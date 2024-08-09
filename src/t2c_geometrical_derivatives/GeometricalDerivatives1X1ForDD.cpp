#include "GeometricalDerivatives1X1ForDD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_dd(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_dd,
                        const size_t idx_op_pp,
                        const size_t idx_op_pf,
                        const size_t idx_op_fp,
                        const size_t idx_op_ff,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PP

        auto to_x_x = pbuffer.data(idx_op_pp + i * 9 + 0);

        auto to_x_y = pbuffer.data(idx_op_pp + i * 9 + 1);

        auto to_x_z = pbuffer.data(idx_op_pp + i * 9 + 2);

        auto to_y_x = pbuffer.data(idx_op_pp + i * 9 + 3);

        auto to_y_y = pbuffer.data(idx_op_pp + i * 9 + 4);

        auto to_y_z = pbuffer.data(idx_op_pp + i * 9 + 5);

        auto to_z_x = pbuffer.data(idx_op_pp + i * 9 + 6);

        auto to_z_y = pbuffer.data(idx_op_pp + i * 9 + 7);

        auto to_z_z = pbuffer.data(idx_op_pp + i * 9 + 8);

        // Set up components of auxiliary buffer : PF

        auto to_x_xxx = pbuffer.data(idx_op_pf + i * 30 + 0);

        auto to_x_xxy = pbuffer.data(idx_op_pf + i * 30 + 1);

        auto to_x_xxz = pbuffer.data(idx_op_pf + i * 30 + 2);

        auto to_x_xyy = pbuffer.data(idx_op_pf + i * 30 + 3);

        auto to_x_xyz = pbuffer.data(idx_op_pf + i * 30 + 4);

        auto to_x_xzz = pbuffer.data(idx_op_pf + i * 30 + 5);

        auto to_x_yyy = pbuffer.data(idx_op_pf + i * 30 + 6);

        auto to_x_yyz = pbuffer.data(idx_op_pf + i * 30 + 7);

        auto to_x_yzz = pbuffer.data(idx_op_pf + i * 30 + 8);

        auto to_x_zzz = pbuffer.data(idx_op_pf + i * 30 + 9);

        auto to_y_xxx = pbuffer.data(idx_op_pf + i * 30 + 10);

        auto to_y_xxy = pbuffer.data(idx_op_pf + i * 30 + 11);

        auto to_y_xxz = pbuffer.data(idx_op_pf + i * 30 + 12);

        auto to_y_xyy = pbuffer.data(idx_op_pf + i * 30 + 13);

        auto to_y_xyz = pbuffer.data(idx_op_pf + i * 30 + 14);

        auto to_y_xzz = pbuffer.data(idx_op_pf + i * 30 + 15);

        auto to_y_yyy = pbuffer.data(idx_op_pf + i * 30 + 16);

        auto to_y_yyz = pbuffer.data(idx_op_pf + i * 30 + 17);

        auto to_y_yzz = pbuffer.data(idx_op_pf + i * 30 + 18);

        auto to_y_zzz = pbuffer.data(idx_op_pf + i * 30 + 19);

        auto to_z_xxx = pbuffer.data(idx_op_pf + i * 30 + 20);

        auto to_z_xxy = pbuffer.data(idx_op_pf + i * 30 + 21);

        auto to_z_xxz = pbuffer.data(idx_op_pf + i * 30 + 22);

        auto to_z_xyy = pbuffer.data(idx_op_pf + i * 30 + 23);

        auto to_z_xyz = pbuffer.data(idx_op_pf + i * 30 + 24);

        auto to_z_xzz = pbuffer.data(idx_op_pf + i * 30 + 25);

        auto to_z_yyy = pbuffer.data(idx_op_pf + i * 30 + 26);

        auto to_z_yyz = pbuffer.data(idx_op_pf + i * 30 + 27);

        auto to_z_yzz = pbuffer.data(idx_op_pf + i * 30 + 28);

        auto to_z_zzz = pbuffer.data(idx_op_pf + i * 30 + 29);

        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up components of auxiliary buffer : FF

        auto to_xxx_xxx = pbuffer.data(idx_op_ff + i * 100 + 0);

        auto to_xxx_xxy = pbuffer.data(idx_op_ff + i * 100 + 1);

        auto to_xxx_xxz = pbuffer.data(idx_op_ff + i * 100 + 2);

        auto to_xxx_xyy = pbuffer.data(idx_op_ff + i * 100 + 3);

        auto to_xxx_xyz = pbuffer.data(idx_op_ff + i * 100 + 4);

        auto to_xxx_xzz = pbuffer.data(idx_op_ff + i * 100 + 5);

        auto to_xxx_yyy = pbuffer.data(idx_op_ff + i * 100 + 6);

        auto to_xxx_yyz = pbuffer.data(idx_op_ff + i * 100 + 7);

        auto to_xxx_yzz = pbuffer.data(idx_op_ff + i * 100 + 8);

        auto to_xxx_zzz = pbuffer.data(idx_op_ff + i * 100 + 9);

        auto to_xxy_xxx = pbuffer.data(idx_op_ff + i * 100 + 10);

        auto to_xxy_xxy = pbuffer.data(idx_op_ff + i * 100 + 11);

        auto to_xxy_xxz = pbuffer.data(idx_op_ff + i * 100 + 12);

        auto to_xxy_xyy = pbuffer.data(idx_op_ff + i * 100 + 13);

        auto to_xxy_xyz = pbuffer.data(idx_op_ff + i * 100 + 14);

        auto to_xxy_xzz = pbuffer.data(idx_op_ff + i * 100 + 15);

        auto to_xxy_yyy = pbuffer.data(idx_op_ff + i * 100 + 16);

        auto to_xxy_yyz = pbuffer.data(idx_op_ff + i * 100 + 17);

        auto to_xxy_yzz = pbuffer.data(idx_op_ff + i * 100 + 18);

        auto to_xxy_zzz = pbuffer.data(idx_op_ff + i * 100 + 19);

        auto to_xxz_xxx = pbuffer.data(idx_op_ff + i * 100 + 20);

        auto to_xxz_xxy = pbuffer.data(idx_op_ff + i * 100 + 21);

        auto to_xxz_xxz = pbuffer.data(idx_op_ff + i * 100 + 22);

        auto to_xxz_xyy = pbuffer.data(idx_op_ff + i * 100 + 23);

        auto to_xxz_xyz = pbuffer.data(idx_op_ff + i * 100 + 24);

        auto to_xxz_xzz = pbuffer.data(idx_op_ff + i * 100 + 25);

        auto to_xxz_yyy = pbuffer.data(idx_op_ff + i * 100 + 26);

        auto to_xxz_yyz = pbuffer.data(idx_op_ff + i * 100 + 27);

        auto to_xxz_yzz = pbuffer.data(idx_op_ff + i * 100 + 28);

        auto to_xxz_zzz = pbuffer.data(idx_op_ff + i * 100 + 29);

        auto to_xyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 30);

        auto to_xyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 31);

        auto to_xyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 32);

        auto to_xyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 33);

        auto to_xyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 34);

        auto to_xyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 35);

        auto to_xyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 36);

        auto to_xyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 37);

        auto to_xyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 38);

        auto to_xyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 39);

        auto to_xyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 40);

        auto to_xyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 41);

        auto to_xyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 42);

        auto to_xyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 43);

        auto to_xyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 44);

        auto to_xyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 45);

        auto to_xyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 46);

        auto to_xyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 47);

        auto to_xyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 48);

        auto to_xyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 49);

        auto to_xzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 50);

        auto to_xzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 51);

        auto to_xzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 52);

        auto to_xzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 53);

        auto to_xzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 54);

        auto to_xzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 55);

        auto to_xzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 56);

        auto to_xzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 57);

        auto to_xzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 58);

        auto to_xzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 59);

        auto to_yyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 60);

        auto to_yyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 61);

        auto to_yyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 62);

        auto to_yyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 63);

        auto to_yyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 64);

        auto to_yyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 65);

        auto to_yyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 66);

        auto to_yyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 67);

        auto to_yyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 68);

        auto to_yyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 69);

        auto to_yyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 70);

        auto to_yyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 71);

        auto to_yyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 72);

        auto to_yyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 73);

        auto to_yyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 74);

        auto to_yyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 75);

        auto to_yyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 76);

        auto to_yyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 77);

        auto to_yyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 78);

        auto to_yyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 79);

        auto to_yzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 80);

        auto to_yzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 81);

        auto to_yzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 82);

        auto to_yzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 83);

        auto to_yzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 84);

        auto to_yzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 85);

        auto to_yzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 86);

        auto to_yzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 87);

        auto to_yzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 88);

        auto to_yzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 89);

        auto to_zzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 90);

        auto to_zzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 91);

        auto to_zzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 92);

        auto to_zzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 93);

        auto to_zzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 94);

        auto to_zzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 95);

        auto to_zzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 96);

        auto to_zzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 97);

        auto to_zzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 98);

        auto to_zzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 99);

        // Set up 0-6 components of targeted buffer : DD

        auto to_x_x_xx_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 0);

        auto to_x_x_xx_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 1);

        auto to_x_x_xx_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 2);

        auto to_x_x_xx_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 3);

        auto to_x_x_xx_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 4);

        auto to_x_x_xx_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_x_x, to_x_x_xx_xx, to_x_x_xx_xy, to_x_x_xx_xz, to_x_x_xx_yy, to_x_x_xx_yz, to_x_x_xx_zz, to_x_xxx, to_x_xxy, to_x_xxz, to_x_xyy, to_x_xyz, to_x_xzz, to_x_y, to_x_z, to_xxx_x, to_xxx_xxx, to_xxx_xxy, to_xxx_xxz, to_xxx_xyy, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xx_xx[k] = 4.0 * to_x_x[k] - 4.0 * to_x_xxx[k] * tke_0 - 4.0 * to_xxx_x[k] * tbe_0 + 4.0 * to_xxx_xxx[k] * tbe_0 * tke_0;

            to_x_x_xx_xy[k] = 2.0 * to_x_y[k] - 4.0 * to_x_xxy[k] * tke_0 - 2.0 * to_xxx_y[k] * tbe_0 + 4.0 * to_xxx_xxy[k] * tbe_0 * tke_0;

            to_x_x_xx_xz[k] = 2.0 * to_x_z[k] - 4.0 * to_x_xxz[k] * tke_0 - 2.0 * to_xxx_z[k] * tbe_0 + 4.0 * to_xxx_xxz[k] * tbe_0 * tke_0;

            to_x_x_xx_yy[k] = -4.0 * to_x_xyy[k] * tke_0 + 4.0 * to_xxx_xyy[k] * tbe_0 * tke_0;

            to_x_x_xx_yz[k] = -4.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xxx_xyz[k] * tbe_0 * tke_0;

            to_x_x_xx_zz[k] = -4.0 * to_x_xzz[k] * tke_0 + 4.0 * to_xxx_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : DD

        auto to_x_x_xy_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 6);

        auto to_x_x_xy_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 7);

        auto to_x_x_xy_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 8);

        auto to_x_x_xy_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 9);

        auto to_x_x_xy_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 10);

        auto to_x_x_xy_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_x_xy_xx, to_x_x_xy_xy, to_x_x_xy_xz, to_x_x_xy_yy, to_x_x_xy_yz, to_x_x_xy_zz, to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, to_y_x, to_y_xxx, to_y_xxy, to_y_xxz, to_y_xyy, to_y_xyz, to_y_xzz, to_y_y, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xy_xx[k] = 2.0 * to_y_x[k] - 2.0 * to_y_xxx[k] * tke_0 - 4.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xy_xy[k] = to_y_y[k] - 2.0 * to_y_xxy[k] * tke_0 - 2.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xy_xz[k] = to_y_z[k] - 2.0 * to_y_xxz[k] * tke_0 - 2.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xy_yy[k] = -2.0 * to_y_xyy[k] * tke_0 + 4.0 * to_xxy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xy_yz[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xy_zz[k] = -2.0 * to_y_xzz[k] * tke_0 + 4.0 * to_xxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : DD

        auto to_x_x_xz_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 12);

        auto to_x_x_xz_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 13);

        auto to_x_x_xz_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 14);

        auto to_x_x_xz_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 15);

        auto to_x_x_xz_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 16);

        auto to_x_x_xz_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_x_xz_xx, to_x_x_xz_xy, to_x_x_xz_xz, to_x_x_xz_yy, to_x_x_xz_yz, to_x_x_xz_zz, to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, to_z_x, to_z_xxx, to_z_xxy, to_z_xxz, to_z_xyy, to_z_xyz, to_z_xzz, to_z_y, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xz_xx[k] = 2.0 * to_z_x[k] - 2.0 * to_z_xxx[k] * tke_0 - 4.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xz_xy[k] = to_z_y[k] - 2.0 * to_z_xxy[k] * tke_0 - 2.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xz_xz[k] = to_z_z[k] - 2.0 * to_z_xxz[k] * tke_0 - 2.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xz_yy[k] = -2.0 * to_z_xyy[k] * tke_0 + 4.0 * to_xxz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xz_yz[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xz_zz[k] = -2.0 * to_z_xzz[k] * tke_0 + 4.0 * to_xxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : DD

        auto to_x_x_yy_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 18);

        auto to_x_x_yy_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 19);

        auto to_x_x_yy_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 20);

        auto to_x_x_yy_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 21);

        auto to_x_x_yy_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 22);

        auto to_x_x_yy_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_x_x_yy_xx, to_x_x_yy_xy, to_x_x_yy_xz, to_x_x_yy_yy, to_x_x_yy_yz, to_x_x_yy_zz, to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yy_xx[k] = -4.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_yy_xy[k] = -2.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_yy_xz[k] = -2.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_yy_yy[k] = 4.0 * to_xyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_yy_yz[k] = 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_yy_zz[k] = 4.0 * to_xyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : DD

        auto to_x_x_yz_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 24);

        auto to_x_x_yz_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 25);

        auto to_x_x_yz_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 26);

        auto to_x_x_yz_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 27);

        auto to_x_x_yz_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 28);

        auto to_x_x_yz_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_x_x_yz_xx, to_x_x_yz_xy, to_x_x_yz_xz, to_x_x_yz_yy, to_x_x_yz_yz, to_x_x_yz_zz, to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yz_xx[k] = -4.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yz_xy[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yz_xz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yz_yy[k] = 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yz_yz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yz_zz[k] = 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : DD

        auto to_x_x_zz_xx = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 30);

        auto to_x_x_zz_xy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 31);

        auto to_x_x_zz_xz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 32);

        auto to_x_x_zz_yy = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 33);

        auto to_x_x_zz_yz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 34);

        auto to_x_x_zz_zz = pbuffer.data(idx_op_geom_101_dd + 0 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_x_x_zz_xx, to_x_x_zz_xy, to_x_x_zz_xz, to_x_x_zz_yy, to_x_x_zz_yz, to_x_x_zz_zz, to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zz_xx[k] = -4.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_zz_xy[k] = -2.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_zz_xz[k] = -2.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_zz_yy[k] = 4.0 * to_xzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_zz_yz[k] = 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_zz_zz[k] = 4.0 * to_xzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : DD

        auto to_x_y_xx_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 0);

        auto to_x_y_xx_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 1);

        auto to_x_y_xx_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 2);

        auto to_x_y_xx_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 3);

        auto to_x_y_xx_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 4);

        auto to_x_y_xx_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_x_x, to_x_xxy, to_x_xyy, to_x_xyz, to_x_y, to_x_y_xx_xx, to_x_y_xx_xy, to_x_y_xx_xz, to_x_y_xx_yy, to_x_y_xx_yz, to_x_y_xx_zz, to_x_yyy, to_x_yyz, to_x_yzz, to_x_z, to_xxx_x, to_xxx_xxy, to_xxx_xyy, to_xxx_xyz, to_xxx_y, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xx_xx[k] = -4.0 * to_x_xxy[k] * tke_0 + 4.0 * to_xxx_xxy[k] * tbe_0 * tke_0;

            to_x_y_xx_xy[k] = 2.0 * to_x_x[k] - 4.0 * to_x_xyy[k] * tke_0 - 2.0 * to_xxx_x[k] * tbe_0 + 4.0 * to_xxx_xyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xz[k] = -4.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xxx_xyz[k] * tbe_0 * tke_0;

            to_x_y_xx_yy[k] = 4.0 * to_x_y[k] - 4.0 * to_x_yyy[k] * tke_0 - 4.0 * to_xxx_y[k] * tbe_0 + 4.0 * to_xxx_yyy[k] * tbe_0 * tke_0;

            to_x_y_xx_yz[k] = 2.0 * to_x_z[k] - 4.0 * to_x_yyz[k] * tke_0 - 2.0 * to_xxx_z[k] * tbe_0 + 4.0 * to_xxx_yyz[k] * tbe_0 * tke_0;

            to_x_y_xx_zz[k] = -4.0 * to_x_yzz[k] * tke_0 + 4.0 * to_xxx_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : DD

        auto to_x_y_xy_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 6);

        auto to_x_y_xy_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 7);

        auto to_x_y_xy_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 8);

        auto to_x_y_xy_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 9);

        auto to_x_y_xy_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 10);

        auto to_x_y_xy_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_y_xy_xx, to_x_y_xy_xy, to_x_y_xy_xz, to_x_y_xy_yy, to_x_y_xy_yz, to_x_y_xy_zz, to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_y_x, to_y_xxy, to_y_xyy, to_y_xyz, to_y_y, to_y_yyy, to_y_yyz, to_y_yzz, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xy_xx[k] = -2.0 * to_y_xxy[k] * tke_0 + 4.0 * to_xxy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xy_xy[k] = to_y_x[k] - 2.0 * to_y_xyy[k] * tke_0 - 2.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xz[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xy_yy[k] = 2.0 * to_y_y[k] - 2.0 * to_y_yyy[k] * tke_0 - 4.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xy_yz[k] = to_y_z[k] - 2.0 * to_y_yyz[k] * tke_0 - 2.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xy_zz[k] = -2.0 * to_y_yzz[k] * tke_0 + 4.0 * to_xxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : DD

        auto to_x_y_xz_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 12);

        auto to_x_y_xz_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 13);

        auto to_x_y_xz_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 14);

        auto to_x_y_xz_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 15);

        auto to_x_y_xz_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 16);

        auto to_x_y_xz_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_y_xz_xx, to_x_y_xz_xy, to_x_y_xz_xz, to_x_y_xz_yy, to_x_y_xz_yz, to_x_y_xz_zz, to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_z_x, to_z_xxy, to_z_xyy, to_z_xyz, to_z_y, to_z_yyy, to_z_yyz, to_z_yzz, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xz_xx[k] = -2.0 * to_z_xxy[k] * tke_0 + 4.0 * to_xxz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xz_xy[k] = to_z_x[k] - 2.0 * to_z_xyy[k] * tke_0 - 2.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xz[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xz_yy[k] = 2.0 * to_z_y[k] - 2.0 * to_z_yyy[k] * tke_0 - 4.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xz_yz[k] = to_z_z[k] - 2.0 * to_z_yyz[k] * tke_0 - 2.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xz_zz[k] = -2.0 * to_z_yzz[k] * tke_0 + 4.0 * to_xxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 54-60 components of targeted buffer : DD

        auto to_x_y_yy_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 18);

        auto to_x_y_yy_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 19);

        auto to_x_y_yy_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 20);

        auto to_x_y_yy_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 21);

        auto to_x_y_yy_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 22);

        auto to_x_y_yy_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_x_y_yy_xx, to_x_y_yy_xy, to_x_y_yy_xz, to_x_y_yy_yy, to_x_y_yy_yz, to_x_y_yy_zz, to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yy_xx[k] = 4.0 * to_xyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_yy_xy[k] = -2.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xz[k] = 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_yy_yy[k] = -4.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_yy_yz[k] = -2.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_yy_zz[k] = 4.0 * to_xyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-66 components of targeted buffer : DD

        auto to_x_y_yz_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 24);

        auto to_x_y_yz_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 25);

        auto to_x_y_yz_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 26);

        auto to_x_y_yz_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 27);

        auto to_x_y_yz_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 28);

        auto to_x_y_yz_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_x_y_yz_xx, to_x_y_yz_xy, to_x_y_yz_xz, to_x_y_yz_yy, to_x_y_yz_yz, to_x_y_yz_zz, to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yz_xx[k] = 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yz_xy[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yz_yy[k] = -4.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yz_yz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yz_zz[k] = 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 66-72 components of targeted buffer : DD

        auto to_x_y_zz_xx = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 30);

        auto to_x_y_zz_xy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 31);

        auto to_x_y_zz_xz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 32);

        auto to_x_y_zz_yy = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 33);

        auto to_x_y_zz_yz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 34);

        auto to_x_y_zz_zz = pbuffer.data(idx_op_geom_101_dd + 1 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_x_y_zz_xx, to_x_y_zz_xy, to_x_y_zz_xz, to_x_y_zz_yy, to_x_y_zz_yz, to_x_y_zz_zz, to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zz_xx[k] = 4.0 * to_xzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_zz_xy[k] = -2.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xz[k] = 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_zz_yy[k] = -4.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_zz_yz[k] = -2.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_zz_zz[k] = 4.0 * to_xzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 72-78 components of targeted buffer : DD

        auto to_x_z_xx_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 0);

        auto to_x_z_xx_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 1);

        auto to_x_z_xx_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 2);

        auto to_x_z_xx_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 3);

        auto to_x_z_xx_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 4);

        auto to_x_z_xx_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_x_x, to_x_xxz, to_x_xyz, to_x_xzz, to_x_y, to_x_yyz, to_x_yzz, to_x_z, to_x_z_xx_xx, to_x_z_xx_xy, to_x_z_xx_xz, to_x_z_xx_yy, to_x_z_xx_yz, to_x_z_xx_zz, to_x_zzz, to_xxx_x, to_xxx_xxz, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xx_xx[k] = -4.0 * to_x_xxz[k] * tke_0 + 4.0 * to_xxx_xxz[k] * tbe_0 * tke_0;

            to_x_z_xx_xy[k] = -4.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xxx_xyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xz[k] = 2.0 * to_x_x[k] - 4.0 * to_x_xzz[k] * tke_0 - 2.0 * to_xxx_x[k] * tbe_0 + 4.0 * to_xxx_xzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yy[k] = -4.0 * to_x_yyz[k] * tke_0 + 4.0 * to_xxx_yyz[k] * tbe_0 * tke_0;

            to_x_z_xx_yz[k] = 2.0 * to_x_y[k] - 4.0 * to_x_yzz[k] * tke_0 - 2.0 * to_xxx_y[k] * tbe_0 + 4.0 * to_xxx_yzz[k] * tbe_0 * tke_0;

            to_x_z_xx_zz[k] = 4.0 * to_x_z[k] - 4.0 * to_x_zzz[k] * tke_0 - 4.0 * to_xxx_z[k] * tbe_0 + 4.0 * to_xxx_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 78-84 components of targeted buffer : DD

        auto to_x_z_xy_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 6);

        auto to_x_z_xy_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 7);

        auto to_x_z_xy_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 8);

        auto to_x_z_xy_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 9);

        auto to_x_z_xy_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 10);

        auto to_x_z_xy_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_z_xy_xx, to_x_z_xy_xy, to_x_z_xy_xz, to_x_z_xy_yy, to_x_z_xy_yz, to_x_z_xy_zz, to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, to_y_x, to_y_xxz, to_y_xyz, to_y_xzz, to_y_y, to_y_yyz, to_y_yzz, to_y_z, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xy_xx[k] = -2.0 * to_y_xxz[k] * tke_0 + 4.0 * to_xxy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xy_xy[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xz[k] = to_y_x[k] - 2.0 * to_y_xzz[k] * tke_0 - 2.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yy[k] = -2.0 * to_y_yyz[k] * tke_0 + 4.0 * to_xxy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xy_yz[k] = to_y_y[k] - 2.0 * to_y_yzz[k] * tke_0 - 2.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xy_zz[k] = 2.0 * to_y_z[k] - 2.0 * to_y_zzz[k] * tke_0 - 4.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 84-90 components of targeted buffer : DD

        auto to_x_z_xz_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 12);

        auto to_x_z_xz_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 13);

        auto to_x_z_xz_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 14);

        auto to_x_z_xz_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 15);

        auto to_x_z_xz_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 16);

        auto to_x_z_xz_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_z_xz_xx, to_x_z_xz_xy, to_x_z_xz_xz, to_x_z_xz_yy, to_x_z_xz_yz, to_x_z_xz_zz, to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, to_z_x, to_z_xxz, to_z_xyz, to_z_xzz, to_z_y, to_z_yyz, to_z_yzz, to_z_z, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xz_xx[k] = -2.0 * to_z_xxz[k] * tke_0 + 4.0 * to_xxz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xz_xy[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xz[k] = to_z_x[k] - 2.0 * to_z_xzz[k] * tke_0 - 2.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yy[k] = -2.0 * to_z_yyz[k] * tke_0 + 4.0 * to_xxz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xz_yz[k] = to_z_y[k] - 2.0 * to_z_yzz[k] * tke_0 - 2.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xz_zz[k] = 2.0 * to_z_z[k] - 2.0 * to_z_zzz[k] * tke_0 - 4.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-96 components of targeted buffer : DD

        auto to_x_z_yy_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 18);

        auto to_x_z_yy_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 19);

        auto to_x_z_yy_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 20);

        auto to_x_z_yy_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 21);

        auto to_x_z_yy_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 22);

        auto to_x_z_yy_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_x_z_yy_xx, to_x_z_yy_xy, to_x_z_yy_xz, to_x_z_yy_yy, to_x_z_yy_yz, to_x_z_yy_zz, to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yy_xx[k] = 4.0 * to_xyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_yy_xy[k] = 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xz[k] = -2.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yy[k] = 4.0 * to_xyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_yy_yz[k] = -2.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_yy_zz[k] = -4.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 96-102 components of targeted buffer : DD

        auto to_x_z_yz_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 24);

        auto to_x_z_yz_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 25);

        auto to_x_z_yz_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 26);

        auto to_x_z_yz_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 27);

        auto to_x_z_yz_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 28);

        auto to_x_z_yz_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_x_z_yz_xx, to_x_z_yz_xy, to_x_z_yz_xz, to_x_z_yz_yy, to_x_z_yz_yz, to_x_z_yz_zz, to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yz_xx[k] = 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yz_xy[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xz[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yy[k] = 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yz_yz[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yz_zz[k] = -4.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 102-108 components of targeted buffer : DD

        auto to_x_z_zz_xx = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 30);

        auto to_x_z_zz_xy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 31);

        auto to_x_z_zz_xz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 32);

        auto to_x_z_zz_yy = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 33);

        auto to_x_z_zz_yz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 34);

        auto to_x_z_zz_zz = pbuffer.data(idx_op_geom_101_dd + 2 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_x_z_zz_xx, to_x_z_zz_xy, to_x_z_zz_xz, to_x_z_zz_yy, to_x_z_zz_yz, to_x_z_zz_zz, to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zz_xx[k] = 4.0 * to_xzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_zz_xy[k] = 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xz[k] = -2.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yy[k] = 4.0 * to_xzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_zz_yz[k] = -2.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_zz_zz[k] = -4.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 108-114 components of targeted buffer : DD

        auto to_y_x_xx_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 0);

        auto to_y_x_xx_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 1);

        auto to_y_x_xx_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 2);

        auto to_y_x_xx_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 3);

        auto to_y_x_xx_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 4);

        auto to_y_x_xx_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, to_y_x_xx_xx, to_y_x_xx_xy, to_y_x_xx_xz, to_y_x_xx_yy, to_y_x_xx_yz, to_y_x_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xx_xx[k] = -4.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xx_xy[k] = -2.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xx_xz[k] = -2.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xx_yy[k] = 4.0 * to_xxy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xx_yz[k] = 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xx_zz[k] = 4.0 * to_xxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 114-120 components of targeted buffer : DD

        auto to_y_x_xy_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 6);

        auto to_y_x_xy_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 7);

        auto to_y_x_xy_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 8);

        auto to_y_x_xy_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 9);

        auto to_y_x_xy_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 10);

        auto to_y_x_xy_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_x, to_x_xxx, to_x_xxy, to_x_xxz, to_x_xyy, to_x_xyz, to_x_xzz, to_x_y, to_x_z, to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, to_y_x_xy_xx, to_y_x_xy_xy, to_y_x_xy_xz, to_y_x_xy_yy, to_y_x_xy_yz, to_y_x_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xy_xx[k] = 2.0 * to_x_x[k] - 2.0 * to_x_xxx[k] * tke_0 - 4.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xy_xy[k] = to_x_y[k] - 2.0 * to_x_xxy[k] * tke_0 - 2.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xy_xz[k] = to_x_z[k] - 2.0 * to_x_xxz[k] * tke_0 - 2.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xy_yy[k] = -2.0 * to_x_xyy[k] * tke_0 + 4.0 * to_xyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xy_yz[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xy_zz[k] = -2.0 * to_x_xzz[k] * tke_0 + 4.0 * to_xyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-126 components of targeted buffer : DD

        auto to_y_x_xz_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 12);

        auto to_y_x_xz_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 13);

        auto to_y_x_xz_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 14);

        auto to_y_x_xz_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 15);

        auto to_y_x_xz_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 16);

        auto to_y_x_xz_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, to_y_x_xz_xx, to_y_x_xz_xy, to_y_x_xz_xz, to_y_x_xz_yy, to_y_x_xz_yz, to_y_x_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xz_xx[k] = -4.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xz_xy[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xz_xz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xz_yy[k] = 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xz_yz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xz_zz[k] = 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 126-132 components of targeted buffer : DD

        auto to_y_x_yy_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 18);

        auto to_y_x_yy_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 19);

        auto to_y_x_yy_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 20);

        auto to_y_x_yy_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 21);

        auto to_y_x_yy_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 22);

        auto to_y_x_yy_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_y_x, to_y_x_yy_xx, to_y_x_yy_xy, to_y_x_yy_xz, to_y_x_yy_yy, to_y_x_yy_yz, to_y_x_yy_zz, to_y_xxx, to_y_xxy, to_y_xxz, to_y_xyy, to_y_xyz, to_y_xzz, to_y_y, to_y_z, to_yyy_x, to_yyy_xxx, to_yyy_xxy, to_yyy_xxz, to_yyy_xyy, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yy_xx[k] = 4.0 * to_y_x[k] - 4.0 * to_y_xxx[k] * tke_0 - 4.0 * to_yyy_x[k] * tbe_0 + 4.0 * to_yyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_yy_xy[k] = 2.0 * to_y_y[k] - 4.0 * to_y_xxy[k] * tke_0 - 2.0 * to_yyy_y[k] * tbe_0 + 4.0 * to_yyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_yy_xz[k] = 2.0 * to_y_z[k] - 4.0 * to_y_xxz[k] * tke_0 - 2.0 * to_yyy_z[k] * tbe_0 + 4.0 * to_yyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_yy_yy[k] = -4.0 * to_y_xyy[k] * tke_0 + 4.0 * to_yyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_yy_yz[k] = -4.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_yy_zz[k] = -4.0 * to_y_xzz[k] * tke_0 + 4.0 * to_yyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 132-138 components of targeted buffer : DD

        auto to_y_x_yz_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 24);

        auto to_y_x_yz_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 25);

        auto to_y_x_yz_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 26);

        auto to_y_x_yz_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 27);

        auto to_y_x_yz_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 28);

        auto to_y_x_yz_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_x_yz_xx, to_y_x_yz_xy, to_y_x_yz_xz, to_y_x_yz_yy, to_y_x_yz_yz, to_y_x_yz_zz, to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, to_z_x, to_z_xxx, to_z_xxy, to_z_xxz, to_z_xyy, to_z_xyz, to_z_xzz, to_z_y, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yz_xx[k] = 2.0 * to_z_x[k] - 2.0 * to_z_xxx[k] * tke_0 - 4.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yz_xy[k] = to_z_y[k] - 2.0 * to_z_xxy[k] * tke_0 - 2.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yz_xz[k] = to_z_z[k] - 2.0 * to_z_xxz[k] * tke_0 - 2.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yz_yy[k] = -2.0 * to_z_xyy[k] * tke_0 + 4.0 * to_yyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yz_yz[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yz_zz[k] = -2.0 * to_z_xzz[k] * tke_0 + 4.0 * to_yyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 138-144 components of targeted buffer : DD

        auto to_y_x_zz_xx = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 30);

        auto to_y_x_zz_xy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 31);

        auto to_y_x_zz_xz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 32);

        auto to_y_x_zz_yy = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 33);

        auto to_y_x_zz_yz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 34);

        auto to_y_x_zz_zz = pbuffer.data(idx_op_geom_101_dd + 3 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_y_x_zz_xx, to_y_x_zz_xy, to_y_x_zz_xz, to_y_x_zz_yy, to_y_x_zz_yz, to_y_x_zz_zz, to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zz_xx[k] = -4.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_zz_xy[k] = -2.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_zz_xz[k] = -2.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_zz_yy[k] = 4.0 * to_yzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_zz_yz[k] = 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_zz_zz[k] = 4.0 * to_yzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 144-150 components of targeted buffer : DD

        auto to_y_y_xx_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 0);

        auto to_y_y_xx_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 1);

        auto to_y_y_xx_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 2);

        auto to_y_y_xx_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 3);

        auto to_y_y_xx_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 4);

        auto to_y_y_xx_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_y_y_xx_xx, to_y_y_xx_xy, to_y_y_xx_xz, to_y_y_xx_yy, to_y_y_xx_yz, to_y_y_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xx_xx[k] = 4.0 * to_xxy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xx_xy[k] = -2.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xz[k] = 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xx_yy[k] = -4.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xx_yz[k] = -2.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xx_zz[k] = 4.0 * to_xxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-156 components of targeted buffer : DD

        auto to_y_y_xy_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 6);

        auto to_y_y_xy_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 7);

        auto to_y_y_xy_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 8);

        auto to_y_y_xy_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 9);

        auto to_y_y_xy_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 10);

        auto to_y_y_xy_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_x, to_x_xxy, to_x_xyy, to_x_xyz, to_x_y, to_x_yyy, to_x_yyz, to_x_yzz, to_x_z, to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_y_y_xy_xx, to_y_y_xy_xy, to_y_y_xy_xz, to_y_y_xy_yy, to_y_y_xy_yz, to_y_y_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xy_xx[k] = -2.0 * to_x_xxy[k] * tke_0 + 4.0 * to_xyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xy_xy[k] = to_x_x[k] - 2.0 * to_x_xyy[k] * tke_0 - 2.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xz[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xy_yy[k] = 2.0 * to_x_y[k] - 2.0 * to_x_yyy[k] * tke_0 - 4.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xy_yz[k] = to_x_z[k] - 2.0 * to_x_yyz[k] * tke_0 - 2.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xy_zz[k] = -2.0 * to_x_yzz[k] * tke_0 + 4.0 * to_xyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 156-162 components of targeted buffer : DD

        auto to_y_y_xz_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 12);

        auto to_y_y_xz_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 13);

        auto to_y_y_xz_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 14);

        auto to_y_y_xz_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 15);

        auto to_y_y_xz_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 16);

        auto to_y_y_xz_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_y_y_xz_xx, to_y_y_xz_xy, to_y_y_xz_xz, to_y_y_xz_yy, to_y_y_xz_yz, to_y_y_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xz_xx[k] = 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xz_xy[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xz_yy[k] = -4.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xz_yz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xz_zz[k] = 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 162-168 components of targeted buffer : DD

        auto to_y_y_yy_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 18);

        auto to_y_y_yy_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 19);

        auto to_y_y_yy_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 20);

        auto to_y_y_yy_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 21);

        auto to_y_y_yy_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 22);

        auto to_y_y_yy_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_y_x, to_y_xxy, to_y_xyy, to_y_xyz, to_y_y, to_y_y_yy_xx, to_y_y_yy_xy, to_y_y_yy_xz, to_y_y_yy_yy, to_y_y_yy_yz, to_y_y_yy_zz, to_y_yyy, to_y_yyz, to_y_yzz, to_y_z, to_yyy_x, to_yyy_xxy, to_yyy_xyy, to_yyy_xyz, to_yyy_y, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yy_xx[k] = -4.0 * to_y_xxy[k] * tke_0 + 4.0 * to_yyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_yy_xy[k] = 2.0 * to_y_x[k] - 4.0 * to_y_xyy[k] * tke_0 - 2.0 * to_yyy_x[k] * tbe_0 + 4.0 * to_yyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xz[k] = -4.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_yy_yy[k] = 4.0 * to_y_y[k] - 4.0 * to_y_yyy[k] * tke_0 - 4.0 * to_yyy_y[k] * tbe_0 + 4.0 * to_yyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_yy_yz[k] = 2.0 * to_y_z[k] - 4.0 * to_y_yyz[k] * tke_0 - 2.0 * to_yyy_z[k] * tbe_0 + 4.0 * to_yyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_yy_zz[k] = -4.0 * to_y_yzz[k] * tke_0 + 4.0 * to_yyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 168-174 components of targeted buffer : DD

        auto to_y_y_yz_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 24);

        auto to_y_y_yz_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 25);

        auto to_y_y_yz_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 26);

        auto to_y_y_yz_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 27);

        auto to_y_y_yz_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 28);

        auto to_y_y_yz_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_y_yz_xx, to_y_y_yz_xy, to_y_y_yz_xz, to_y_y_yz_yy, to_y_y_yz_yz, to_y_y_yz_zz, to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_z_x, to_z_xxy, to_z_xyy, to_z_xyz, to_z_y, to_z_yyy, to_z_yyz, to_z_yzz, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yz_xx[k] = -2.0 * to_z_xxy[k] * tke_0 + 4.0 * to_yyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yz_xy[k] = to_z_x[k] - 2.0 * to_z_xyy[k] * tke_0 - 2.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xz[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yz_yy[k] = 2.0 * to_z_y[k] - 2.0 * to_z_yyy[k] * tke_0 - 4.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yz_yz[k] = to_z_z[k] - 2.0 * to_z_yyz[k] * tke_0 - 2.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yz_zz[k] = -2.0 * to_z_yzz[k] * tke_0 + 4.0 * to_yyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 174-180 components of targeted buffer : DD

        auto to_y_y_zz_xx = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 30);

        auto to_y_y_zz_xy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 31);

        auto to_y_y_zz_xz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 32);

        auto to_y_y_zz_yy = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 33);

        auto to_y_y_zz_yz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 34);

        auto to_y_y_zz_zz = pbuffer.data(idx_op_geom_101_dd + 4 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_y_y_zz_xx, to_y_y_zz_xy, to_y_y_zz_xz, to_y_y_zz_yy, to_y_y_zz_yz, to_y_y_zz_zz, to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zz_xx[k] = 4.0 * to_yzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_zz_xy[k] = -2.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xz[k] = 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_zz_yy[k] = -4.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_zz_yz[k] = -2.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_zz_zz[k] = 4.0 * to_yzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-186 components of targeted buffer : DD

        auto to_y_z_xx_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 0);

        auto to_y_z_xx_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 1);

        auto to_y_z_xx_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 2);

        auto to_y_z_xx_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 3);

        auto to_y_z_xx_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 4);

        auto to_y_z_xx_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, to_y_z_xx_xx, to_y_z_xx_xy, to_y_z_xx_xz, to_y_z_xx_yy, to_y_z_xx_yz, to_y_z_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xx_xx[k] = 4.0 * to_xxy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xx_xy[k] = 4.0 * to_xxy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xz[k] = -2.0 * to_xxy_x[k] * tbe_0 + 4.0 * to_xxy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yy[k] = 4.0 * to_xxy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xx_yz[k] = -2.0 * to_xxy_y[k] * tbe_0 + 4.0 * to_xxy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xx_zz[k] = -4.0 * to_xxy_z[k] * tbe_0 + 4.0 * to_xxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 186-192 components of targeted buffer : DD

        auto to_y_z_xy_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 6);

        auto to_y_z_xy_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 7);

        auto to_y_z_xy_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 8);

        auto to_y_z_xy_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 9);

        auto to_y_z_xy_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 10);

        auto to_y_z_xy_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_x_x, to_x_xxz, to_x_xyz, to_x_xzz, to_x_y, to_x_yyz, to_x_yzz, to_x_z, to_x_zzz, to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, to_y_z_xy_xx, to_y_z_xy_xy, to_y_z_xy_xz, to_y_z_xy_yy, to_y_z_xy_yz, to_y_z_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xy_xx[k] = -2.0 * to_x_xxz[k] * tke_0 + 4.0 * to_xyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xy_xy[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xz[k] = to_x_x[k] - 2.0 * to_x_xzz[k] * tke_0 - 2.0 * to_xyy_x[k] * tbe_0 + 4.0 * to_xyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yy[k] = -2.0 * to_x_yyz[k] * tke_0 + 4.0 * to_xyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xy_yz[k] = to_x_y[k] - 2.0 * to_x_yzz[k] * tke_0 - 2.0 * to_xyy_y[k] * tbe_0 + 4.0 * to_xyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xy_zz[k] = 2.0 * to_x_z[k] - 2.0 * to_x_zzz[k] * tke_0 - 4.0 * to_xyy_z[k] * tbe_0 + 4.0 * to_xyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 192-198 components of targeted buffer : DD

        auto to_y_z_xz_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 12);

        auto to_y_z_xz_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 13);

        auto to_y_z_xz_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 14);

        auto to_y_z_xz_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 15);

        auto to_y_z_xz_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 16);

        auto to_y_z_xz_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, to_y_z_xz_xx, to_y_z_xz_xy, to_y_z_xz_xz, to_y_z_xz_yy, to_y_z_xz_yz, to_y_z_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xz_xx[k] = 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xz_xy[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xz[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yy[k] = 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xz_yz[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xz_zz[k] = -4.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 198-204 components of targeted buffer : DD

        auto to_y_z_yy_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 18);

        auto to_y_z_yy_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 19);

        auto to_y_z_yy_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 20);

        auto to_y_z_yy_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 21);

        auto to_y_z_yy_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 22);

        auto to_y_z_yy_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_y_x, to_y_xxz, to_y_xyz, to_y_xzz, to_y_y, to_y_yyz, to_y_yzz, to_y_z, to_y_z_yy_xx, to_y_z_yy_xy, to_y_z_yy_xz, to_y_z_yy_yy, to_y_z_yy_yz, to_y_z_yy_zz, to_y_zzz, to_yyy_x, to_yyy_xxz, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yy_xx[k] = -4.0 * to_y_xxz[k] * tke_0 + 4.0 * to_yyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_yy_xy[k] = -4.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xz[k] = 2.0 * to_y_x[k] - 4.0 * to_y_xzz[k] * tke_0 - 2.0 * to_yyy_x[k] * tbe_0 + 4.0 * to_yyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yy[k] = -4.0 * to_y_yyz[k] * tke_0 + 4.0 * to_yyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_yy_yz[k] = 2.0 * to_y_y[k] - 4.0 * to_y_yzz[k] * tke_0 - 2.0 * to_yyy_y[k] * tbe_0 + 4.0 * to_yyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_yy_zz[k] = 4.0 * to_y_z[k] - 4.0 * to_y_zzz[k] * tke_0 - 4.0 * to_yyy_z[k] * tbe_0 + 4.0 * to_yyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 204-210 components of targeted buffer : DD

        auto to_y_z_yz_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 24);

        auto to_y_z_yz_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 25);

        auto to_y_z_yz_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 26);

        auto to_y_z_yz_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 27);

        auto to_y_z_yz_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 28);

        auto to_y_z_yz_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_z_yz_xx, to_y_z_yz_xy, to_y_z_yz_xz, to_y_z_yz_yy, to_y_z_yz_yz, to_y_z_yz_zz, to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, to_z_x, to_z_xxz, to_z_xyz, to_z_xzz, to_z_y, to_z_yyz, to_z_yzz, to_z_z, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yz_xx[k] = -2.0 * to_z_xxz[k] * tke_0 + 4.0 * to_yyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yz_xy[k] = -2.0 * to_z_xyz[k] * tke_0 + 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xz[k] = to_z_x[k] - 2.0 * to_z_xzz[k] * tke_0 - 2.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yy[k] = -2.0 * to_z_yyz[k] * tke_0 + 4.0 * to_yyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yz_yz[k] = to_z_y[k] - 2.0 * to_z_yzz[k] * tke_0 - 2.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yz_zz[k] = 2.0 * to_z_z[k] - 2.0 * to_z_zzz[k] * tke_0 - 4.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-216 components of targeted buffer : DD

        auto to_y_z_zz_xx = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 30);

        auto to_y_z_zz_xy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 31);

        auto to_y_z_zz_xz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 32);

        auto to_y_z_zz_yy = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 33);

        auto to_y_z_zz_yz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 34);

        auto to_y_z_zz_zz = pbuffer.data(idx_op_geom_101_dd + 5 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_y_z_zz_xx, to_y_z_zz_xy, to_y_z_zz_xz, to_y_z_zz_yy, to_y_z_zz_yz, to_y_z_zz_zz, to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zz_xx[k] = 4.0 * to_yzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_zz_xy[k] = 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xz[k] = -2.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yy[k] = 4.0 * to_yzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_zz_yz[k] = -2.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_zz_zz[k] = -4.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 216-222 components of targeted buffer : DD

        auto to_z_x_xx_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 0);

        auto to_z_x_xx_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 1);

        auto to_z_x_xx_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 2);

        auto to_z_x_xx_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 3);

        auto to_z_x_xx_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 4);

        auto to_z_x_xx_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, to_z_x_xx_xx, to_z_x_xx_xy, to_z_x_xx_xz, to_z_x_xx_yy, to_z_x_xx_yz, to_z_x_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xx_xx[k] = -4.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xx_xy[k] = -2.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xx_xz[k] = -2.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xx_yy[k] = 4.0 * to_xxz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xx_yz[k] = 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xx_zz[k] = 4.0 * to_xxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 222-228 components of targeted buffer : DD

        auto to_z_x_xy_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 6);

        auto to_z_x_xy_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 7);

        auto to_z_x_xy_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 8);

        auto to_z_x_xy_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 9);

        auto to_z_x_xy_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 10);

        auto to_z_x_xy_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, to_z_x_xy_xx, to_z_x_xy_xy, to_z_x_xy_xz, to_z_x_xy_yy, to_z_x_xy_yz, to_z_x_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xy_xx[k] = -4.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xy_xy[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xy_xz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xy_yy[k] = 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xy_yz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xy_zz[k] = 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 228-234 components of targeted buffer : DD

        auto to_z_x_xz_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 12);

        auto to_z_x_xz_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 13);

        auto to_z_x_xz_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 14);

        auto to_z_x_xz_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 15);

        auto to_z_x_xz_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 16);

        auto to_z_x_xz_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_x, to_x_xxx, to_x_xxy, to_x_xxz, to_x_xyy, to_x_xyz, to_x_xzz, to_x_y, to_x_z, to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, to_z_x_xz_xx, to_z_x_xz_xy, to_z_x_xz_xz, to_z_x_xz_yy, to_z_x_xz_yz, to_z_x_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xz_xx[k] = 2.0 * to_x_x[k] - 2.0 * to_x_xxx[k] * tke_0 - 4.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xz_xy[k] = to_x_y[k] - 2.0 * to_x_xxy[k] * tke_0 - 2.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xz_xz[k] = to_x_z[k] - 2.0 * to_x_xxz[k] * tke_0 - 2.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xz_yy[k] = -2.0 * to_x_xyy[k] * tke_0 + 4.0 * to_xzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xz_yz[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xz_zz[k] = -2.0 * to_x_xzz[k] * tke_0 + 4.0 * to_xzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 234-240 components of targeted buffer : DD

        auto to_z_x_yy_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 18);

        auto to_z_x_yy_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 19);

        auto to_z_x_yy_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 20);

        auto to_z_x_yy_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 21);

        auto to_z_x_yy_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 22);

        auto to_z_x_yy_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, to_z_x_yy_xx, to_z_x_yy_xy, to_z_x_yy_xz, to_z_x_yy_yy, to_z_x_yy_yz, to_z_x_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yy_xx[k] = -4.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yy_xy[k] = -2.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yy_xz[k] = -2.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yy_yy[k] = 4.0 * to_yyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yy_yz[k] = 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yy_zz[k] = 4.0 * to_yyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-246 components of targeted buffer : DD

        auto to_z_x_yz_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 24);

        auto to_z_x_yz_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 25);

        auto to_z_x_yz_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 26);

        auto to_z_x_yz_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 27);

        auto to_z_x_yz_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 28);

        auto to_z_x_yz_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_x, to_y_xxx, to_y_xxy, to_y_xxz, to_y_xyy, to_y_xyz, to_y_xzz, to_y_y, to_y_z, to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, to_z_x_yz_xx, to_z_x_yz_xy, to_z_x_yz_xz, to_z_x_yz_yy, to_z_x_yz_yz, to_z_x_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yz_xx[k] = 2.0 * to_y_x[k] - 2.0 * to_y_xxx[k] * tke_0 - 4.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yz_xy[k] = to_y_y[k] - 2.0 * to_y_xxy[k] * tke_0 - 2.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yz_xz[k] = to_y_z[k] - 2.0 * to_y_xxz[k] * tke_0 - 2.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yz_yy[k] = -2.0 * to_y_xyy[k] * tke_0 + 4.0 * to_yzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yz_yz[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yz_zz[k] = -2.0 * to_y_xzz[k] * tke_0 + 4.0 * to_yzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 246-252 components of targeted buffer : DD

        auto to_z_x_zz_xx = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 30);

        auto to_z_x_zz_xy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 31);

        auto to_z_x_zz_xz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 32);

        auto to_z_x_zz_yy = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 33);

        auto to_z_x_zz_yz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 34);

        auto to_z_x_zz_zz = pbuffer.data(idx_op_geom_101_dd + 6 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_z_x, to_z_x_zz_xx, to_z_x_zz_xy, to_z_x_zz_xz, to_z_x_zz_yy, to_z_x_zz_yz, to_z_x_zz_zz, to_z_xxx, to_z_xxy, to_z_xxz, to_z_xyy, to_z_xyz, to_z_xzz, to_z_y, to_z_z, to_zzz_x, to_zzz_xxx, to_zzz_xxy, to_zzz_xxz, to_zzz_xyy, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zz_xx[k] = 4.0 * to_z_x[k] - 4.0 * to_z_xxx[k] * tke_0 - 4.0 * to_zzz_x[k] * tbe_0 + 4.0 * to_zzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_zz_xy[k] = 2.0 * to_z_y[k] - 4.0 * to_z_xxy[k] * tke_0 - 2.0 * to_zzz_y[k] * tbe_0 + 4.0 * to_zzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_zz_xz[k] = 2.0 * to_z_z[k] - 4.0 * to_z_xxz[k] * tke_0 - 2.0 * to_zzz_z[k] * tbe_0 + 4.0 * to_zzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_zz_yy[k] = -4.0 * to_z_xyy[k] * tke_0 + 4.0 * to_zzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_zz_yz[k] = -4.0 * to_z_xyz[k] * tke_0 + 4.0 * to_zzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_zz_zz[k] = -4.0 * to_z_xzz[k] * tke_0 + 4.0 * to_zzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 252-258 components of targeted buffer : DD

        auto to_z_y_xx_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 0);

        auto to_z_y_xx_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 1);

        auto to_z_y_xx_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 2);

        auto to_z_y_xx_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 3);

        auto to_z_y_xx_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 4);

        auto to_z_y_xx_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_z_y_xx_xx, to_z_y_xx_xy, to_z_y_xx_xz, to_z_y_xx_yy, to_z_y_xx_yz, to_z_y_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xx_xx[k] = 4.0 * to_xxz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xx_xy[k] = -2.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xz[k] = 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xx_yy[k] = -4.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xx_yz[k] = -2.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xx_zz[k] = 4.0 * to_xxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 258-264 components of targeted buffer : DD

        auto to_z_y_xy_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 6);

        auto to_z_y_xy_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 7);

        auto to_z_y_xy_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 8);

        auto to_z_y_xy_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 9);

        auto to_z_y_xy_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 10);

        auto to_z_y_xy_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_z_y_xy_xx, to_z_y_xy_xy, to_z_y_xy_xz, to_z_y_xy_yy, to_z_y_xy_yz, to_z_y_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xy_xx[k] = 4.0 * to_xyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xy_xy[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xz[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xy_yy[k] = -4.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xy_yz[k] = -2.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xy_zz[k] = 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 264-270 components of targeted buffer : DD

        auto to_z_y_xz_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 12);

        auto to_z_y_xz_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 13);

        auto to_z_y_xz_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 14);

        auto to_z_y_xz_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 15);

        auto to_z_y_xz_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 16);

        auto to_z_y_xz_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_x, to_x_xxy, to_x_xyy, to_x_xyz, to_x_y, to_x_yyy, to_x_yyz, to_x_yzz, to_x_z, to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_z_y_xz_xx, to_z_y_xz_xy, to_z_y_xz_xz, to_z_y_xz_yy, to_z_y_xz_yz, to_z_y_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xz_xx[k] = -2.0 * to_x_xxy[k] * tke_0 + 4.0 * to_xzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xz_xy[k] = to_x_x[k] - 2.0 * to_x_xyy[k] * tke_0 - 2.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xz[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xz_yy[k] = 2.0 * to_x_y[k] - 2.0 * to_x_yyy[k] * tke_0 - 4.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xz_yz[k] = to_x_z[k] - 2.0 * to_x_yyz[k] * tke_0 - 2.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xz_zz[k] = -2.0 * to_x_yzz[k] * tke_0 + 4.0 * to_xzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-276 components of targeted buffer : DD

        auto to_z_y_yy_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 18);

        auto to_z_y_yy_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 19);

        auto to_z_y_yy_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 20);

        auto to_z_y_yy_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 21);

        auto to_z_y_yy_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 22);

        auto to_z_y_yy_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_z_y_yy_xx, to_z_y_yy_xy, to_z_y_yy_xz, to_z_y_yy_yy, to_z_y_yy_yz, to_z_y_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yy_xx[k] = 4.0 * to_yyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yy_xy[k] = -2.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xz[k] = 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yy_yy[k] = -4.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yy_yz[k] = -2.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yy_zz[k] = 4.0 * to_yyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 276-282 components of targeted buffer : DD

        auto to_z_y_yz_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 24);

        auto to_z_y_yz_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 25);

        auto to_z_y_yz_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 26);

        auto to_z_y_yz_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 27);

        auto to_z_y_yz_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 28);

        auto to_z_y_yz_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_x, to_y_xxy, to_y_xyy, to_y_xyz, to_y_y, to_y_yyy, to_y_yyz, to_y_yzz, to_y_z, to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_z_y_yz_xx, to_z_y_yz_xy, to_z_y_yz_xz, to_z_y_yz_yy, to_z_y_yz_yz, to_z_y_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yz_xx[k] = -2.0 * to_y_xxy[k] * tke_0 + 4.0 * to_yzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yz_xy[k] = to_y_x[k] - 2.0 * to_y_xyy[k] * tke_0 - 2.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xz[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yz_yy[k] = 2.0 * to_y_y[k] - 2.0 * to_y_yyy[k] * tke_0 - 4.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yz_yz[k] = to_y_z[k] - 2.0 * to_y_yyz[k] * tke_0 - 2.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yz_zz[k] = -2.0 * to_y_yzz[k] * tke_0 + 4.0 * to_yzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 282-288 components of targeted buffer : DD

        auto to_z_y_zz_xx = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 30);

        auto to_z_y_zz_xy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 31);

        auto to_z_y_zz_xz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 32);

        auto to_z_y_zz_yy = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 33);

        auto to_z_y_zz_yz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 34);

        auto to_z_y_zz_zz = pbuffer.data(idx_op_geom_101_dd + 7 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_z_x, to_z_xxy, to_z_xyy, to_z_xyz, to_z_y, to_z_y_zz_xx, to_z_y_zz_xy, to_z_y_zz_xz, to_z_y_zz_yy, to_z_y_zz_yz, to_z_y_zz_zz, to_z_yyy, to_z_yyz, to_z_yzz, to_z_z, to_zzz_x, to_zzz_xxy, to_zzz_xyy, to_zzz_xyz, to_zzz_y, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zz_xx[k] = -4.0 * to_z_xxy[k] * tke_0 + 4.0 * to_zzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_zz_xy[k] = 2.0 * to_z_x[k] - 4.0 * to_z_xyy[k] * tke_0 - 2.0 * to_zzz_x[k] * tbe_0 + 4.0 * to_zzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xz[k] = -4.0 * to_z_xyz[k] * tke_0 + 4.0 * to_zzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_zz_yy[k] = 4.0 * to_z_y[k] - 4.0 * to_z_yyy[k] * tke_0 - 4.0 * to_zzz_y[k] * tbe_0 + 4.0 * to_zzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_zz_yz[k] = 2.0 * to_z_z[k] - 4.0 * to_z_yyz[k] * tke_0 - 2.0 * to_zzz_z[k] * tbe_0 + 4.0 * to_zzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_zz_zz[k] = -4.0 * to_z_yzz[k] * tke_0 + 4.0 * to_zzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 288-294 components of targeted buffer : DD

        auto to_z_z_xx_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 0);

        auto to_z_z_xx_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 1);

        auto to_z_z_xx_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 2);

        auto to_z_z_xx_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 3);

        auto to_z_z_xx_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 4);

        auto to_z_z_xx_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, to_z_z_xx_xx, to_z_z_xx_xy, to_z_z_xx_xz, to_z_z_xx_yy, to_z_z_xx_yz, to_z_z_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xx_xx[k] = 4.0 * to_xxz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xx_xy[k] = 4.0 * to_xxz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xz[k] = -2.0 * to_xxz_x[k] * tbe_0 + 4.0 * to_xxz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yy[k] = 4.0 * to_xxz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xx_yz[k] = -2.0 * to_xxz_y[k] * tbe_0 + 4.0 * to_xxz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xx_zz[k] = -4.0 * to_xxz_z[k] * tbe_0 + 4.0 * to_xxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 294-300 components of targeted buffer : DD

        auto to_z_z_xy_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 6);

        auto to_z_z_xy_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 7);

        auto to_z_z_xy_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 8);

        auto to_z_z_xy_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 9);

        auto to_z_z_xy_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 10);

        auto to_z_z_xy_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, to_z_z_xy_xx, to_z_z_xy_xy, to_z_z_xy_xz, to_z_z_xy_yy, to_z_z_xy_yz, to_z_z_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xy_xx[k] = 4.0 * to_xyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xy_xy[k] = 4.0 * to_xyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xz[k] = -2.0 * to_xyz_x[k] * tbe_0 + 4.0 * to_xyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yy[k] = 4.0 * to_xyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xy_yz[k] = -2.0 * to_xyz_y[k] * tbe_0 + 4.0 * to_xyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xy_zz[k] = -4.0 * to_xyz_z[k] * tbe_0 + 4.0 * to_xyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-306 components of targeted buffer : DD

        auto to_z_z_xz_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 12);

        auto to_z_z_xz_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 13);

        auto to_z_z_xz_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 14);

        auto to_z_z_xz_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 15);

        auto to_z_z_xz_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 16);

        auto to_z_z_xz_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_x_x, to_x_xxz, to_x_xyz, to_x_xzz, to_x_y, to_x_yyz, to_x_yzz, to_x_z, to_x_zzz, to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, to_z_z_xz_xx, to_z_z_xz_xy, to_z_z_xz_xz, to_z_z_xz_yy, to_z_z_xz_yz, to_z_z_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xz_xx[k] = -2.0 * to_x_xxz[k] * tke_0 + 4.0 * to_xzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xz_xy[k] = -2.0 * to_x_xyz[k] * tke_0 + 4.0 * to_xzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xz[k] = to_x_x[k] - 2.0 * to_x_xzz[k] * tke_0 - 2.0 * to_xzz_x[k] * tbe_0 + 4.0 * to_xzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yy[k] = -2.0 * to_x_yyz[k] * tke_0 + 4.0 * to_xzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xz_yz[k] = to_x_y[k] - 2.0 * to_x_yzz[k] * tke_0 - 2.0 * to_xzz_y[k] * tbe_0 + 4.0 * to_xzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xz_zz[k] = 2.0 * to_x_z[k] - 2.0 * to_x_zzz[k] * tke_0 - 4.0 * to_xzz_z[k] * tbe_0 + 4.0 * to_xzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 306-312 components of targeted buffer : DD

        auto to_z_z_yy_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 18);

        auto to_z_z_yy_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 19);

        auto to_z_z_yy_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 20);

        auto to_z_z_yy_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 21);

        auto to_z_z_yy_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 22);

        auto to_z_z_yy_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, to_z_z_yy_xx, to_z_z_yy_xy, to_z_z_yy_xz, to_z_z_yy_yy, to_z_z_yy_yz, to_z_z_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yy_xx[k] = 4.0 * to_yyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yy_xy[k] = 4.0 * to_yyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xz[k] = -2.0 * to_yyz_x[k] * tbe_0 + 4.0 * to_yyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yy[k] = 4.0 * to_yyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yy_yz[k] = -2.0 * to_yyz_y[k] * tbe_0 + 4.0 * to_yyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yy_zz[k] = -4.0 * to_yyz_z[k] * tbe_0 + 4.0 * to_yyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 312-318 components of targeted buffer : DD

        auto to_z_z_yz_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 24);

        auto to_z_z_yz_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 25);

        auto to_z_z_yz_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 26);

        auto to_z_z_yz_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 27);

        auto to_z_z_yz_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 28);

        auto to_z_z_yz_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_y_x, to_y_xxz, to_y_xyz, to_y_xzz, to_y_y, to_y_yyz, to_y_yzz, to_y_z, to_y_zzz, to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, to_z_z_yz_xx, to_z_z_yz_xy, to_z_z_yz_xz, to_z_z_yz_yy, to_z_z_yz_yz, to_z_z_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yz_xx[k] = -2.0 * to_y_xxz[k] * tke_0 + 4.0 * to_yzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yz_xy[k] = -2.0 * to_y_xyz[k] * tke_0 + 4.0 * to_yzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xz[k] = to_y_x[k] - 2.0 * to_y_xzz[k] * tke_0 - 2.0 * to_yzz_x[k] * tbe_0 + 4.0 * to_yzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yy[k] = -2.0 * to_y_yyz[k] * tke_0 + 4.0 * to_yzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yz_yz[k] = to_y_y[k] - 2.0 * to_y_yzz[k] * tke_0 - 2.0 * to_yzz_y[k] * tbe_0 + 4.0 * to_yzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yz_zz[k] = 2.0 * to_y_z[k] - 2.0 * to_y_zzz[k] * tke_0 - 4.0 * to_yzz_z[k] * tbe_0 + 4.0 * to_yzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 318-324 components of targeted buffer : DD

        auto to_z_z_zz_xx = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 30);

        auto to_z_z_zz_xy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 31);

        auto to_z_z_zz_xz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 32);

        auto to_z_z_zz_yy = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 33);

        auto to_z_z_zz_yz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 34);

        auto to_z_z_zz_zz = pbuffer.data(idx_op_geom_101_dd + 8 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_z_x, to_z_xxz, to_z_xyz, to_z_xzz, to_z_y, to_z_yyz, to_z_yzz, to_z_z, to_z_z_zz_xx, to_z_z_zz_xy, to_z_z_zz_xz, to_z_z_zz_yy, to_z_z_zz_yz, to_z_z_zz_zz, to_z_zzz, to_zzz_x, to_zzz_xxz, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zz_xx[k] = -4.0 * to_z_xxz[k] * tke_0 + 4.0 * to_zzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_zz_xy[k] = -4.0 * to_z_xyz[k] * tke_0 + 4.0 * to_zzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xz[k] = 2.0 * to_z_x[k] - 4.0 * to_z_xzz[k] * tke_0 - 2.0 * to_zzz_x[k] * tbe_0 + 4.0 * to_zzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yy[k] = -4.0 * to_z_yyz[k] * tke_0 + 4.0 * to_zzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_zz_yz[k] = 2.0 * to_z_y[k] - 4.0 * to_z_yzz[k] * tke_0 - 2.0 * to_zzz_y[k] * tbe_0 + 4.0 * to_zzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_zz_zz[k] = 4.0 * to_z_z[k] - 4.0 * to_z_zzz[k] * tke_0 - 4.0 * to_zzz_z[k] * tbe_0 + 4.0 * to_zzz_zzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

