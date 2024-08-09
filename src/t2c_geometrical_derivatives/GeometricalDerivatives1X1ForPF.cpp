#include "GeometricalDerivatives1X1ForPF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_pf(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_pf,
                        const size_t idx_op_sd,
                        const size_t idx_op_sg,
                        const size_t idx_op_dd,
                        const size_t idx_op_dg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SD

        auto to_0_xx = pbuffer.data(idx_op_sd + i * 6 + 0);

        auto to_0_xy = pbuffer.data(idx_op_sd + i * 6 + 1);

        auto to_0_xz = pbuffer.data(idx_op_sd + i * 6 + 2);

        auto to_0_yy = pbuffer.data(idx_op_sd + i * 6 + 3);

        auto to_0_yz = pbuffer.data(idx_op_sd + i * 6 + 4);

        auto to_0_zz = pbuffer.data(idx_op_sd + i * 6 + 5);

        // Set up components of auxiliary buffer : SG

        auto to_0_xxxx = pbuffer.data(idx_op_sg + i * 15 + 0);

        auto to_0_xxxy = pbuffer.data(idx_op_sg + i * 15 + 1);

        auto to_0_xxxz = pbuffer.data(idx_op_sg + i * 15 + 2);

        auto to_0_xxyy = pbuffer.data(idx_op_sg + i * 15 + 3);

        auto to_0_xxyz = pbuffer.data(idx_op_sg + i * 15 + 4);

        auto to_0_xxzz = pbuffer.data(idx_op_sg + i * 15 + 5);

        auto to_0_xyyy = pbuffer.data(idx_op_sg + i * 15 + 6);

        auto to_0_xyyz = pbuffer.data(idx_op_sg + i * 15 + 7);

        auto to_0_xyzz = pbuffer.data(idx_op_sg + i * 15 + 8);

        auto to_0_xzzz = pbuffer.data(idx_op_sg + i * 15 + 9);

        auto to_0_yyyy = pbuffer.data(idx_op_sg + i * 15 + 10);

        auto to_0_yyyz = pbuffer.data(idx_op_sg + i * 15 + 11);

        auto to_0_yyzz = pbuffer.data(idx_op_sg + i * 15 + 12);

        auto to_0_yzzz = pbuffer.data(idx_op_sg + i * 15 + 13);

        auto to_0_zzzz = pbuffer.data(idx_op_sg + i * 15 + 14);

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

        // Set up components of auxiliary buffer : DG

        auto to_xx_xxxx = pbuffer.data(idx_op_dg + i * 90 + 0);

        auto to_xx_xxxy = pbuffer.data(idx_op_dg + i * 90 + 1);

        auto to_xx_xxxz = pbuffer.data(idx_op_dg + i * 90 + 2);

        auto to_xx_xxyy = pbuffer.data(idx_op_dg + i * 90 + 3);

        auto to_xx_xxyz = pbuffer.data(idx_op_dg + i * 90 + 4);

        auto to_xx_xxzz = pbuffer.data(idx_op_dg + i * 90 + 5);

        auto to_xx_xyyy = pbuffer.data(idx_op_dg + i * 90 + 6);

        auto to_xx_xyyz = pbuffer.data(idx_op_dg + i * 90 + 7);

        auto to_xx_xyzz = pbuffer.data(idx_op_dg + i * 90 + 8);

        auto to_xx_xzzz = pbuffer.data(idx_op_dg + i * 90 + 9);

        auto to_xx_yyyy = pbuffer.data(idx_op_dg + i * 90 + 10);

        auto to_xx_yyyz = pbuffer.data(idx_op_dg + i * 90 + 11);

        auto to_xx_yyzz = pbuffer.data(idx_op_dg + i * 90 + 12);

        auto to_xx_yzzz = pbuffer.data(idx_op_dg + i * 90 + 13);

        auto to_xx_zzzz = pbuffer.data(idx_op_dg + i * 90 + 14);

        auto to_xy_xxxx = pbuffer.data(idx_op_dg + i * 90 + 15);

        auto to_xy_xxxy = pbuffer.data(idx_op_dg + i * 90 + 16);

        auto to_xy_xxxz = pbuffer.data(idx_op_dg + i * 90 + 17);

        auto to_xy_xxyy = pbuffer.data(idx_op_dg + i * 90 + 18);

        auto to_xy_xxyz = pbuffer.data(idx_op_dg + i * 90 + 19);

        auto to_xy_xxzz = pbuffer.data(idx_op_dg + i * 90 + 20);

        auto to_xy_xyyy = pbuffer.data(idx_op_dg + i * 90 + 21);

        auto to_xy_xyyz = pbuffer.data(idx_op_dg + i * 90 + 22);

        auto to_xy_xyzz = pbuffer.data(idx_op_dg + i * 90 + 23);

        auto to_xy_xzzz = pbuffer.data(idx_op_dg + i * 90 + 24);

        auto to_xy_yyyy = pbuffer.data(idx_op_dg + i * 90 + 25);

        auto to_xy_yyyz = pbuffer.data(idx_op_dg + i * 90 + 26);

        auto to_xy_yyzz = pbuffer.data(idx_op_dg + i * 90 + 27);

        auto to_xy_yzzz = pbuffer.data(idx_op_dg + i * 90 + 28);

        auto to_xy_zzzz = pbuffer.data(idx_op_dg + i * 90 + 29);

        auto to_xz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 30);

        auto to_xz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 31);

        auto to_xz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 32);

        auto to_xz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 33);

        auto to_xz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 34);

        auto to_xz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 35);

        auto to_xz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 36);

        auto to_xz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 37);

        auto to_xz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 38);

        auto to_xz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 39);

        auto to_xz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 40);

        auto to_xz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 41);

        auto to_xz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 42);

        auto to_xz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 43);

        auto to_xz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 44);

        auto to_yy_xxxx = pbuffer.data(idx_op_dg + i * 90 + 45);

        auto to_yy_xxxy = pbuffer.data(idx_op_dg + i * 90 + 46);

        auto to_yy_xxxz = pbuffer.data(idx_op_dg + i * 90 + 47);

        auto to_yy_xxyy = pbuffer.data(idx_op_dg + i * 90 + 48);

        auto to_yy_xxyz = pbuffer.data(idx_op_dg + i * 90 + 49);

        auto to_yy_xxzz = pbuffer.data(idx_op_dg + i * 90 + 50);

        auto to_yy_xyyy = pbuffer.data(idx_op_dg + i * 90 + 51);

        auto to_yy_xyyz = pbuffer.data(idx_op_dg + i * 90 + 52);

        auto to_yy_xyzz = pbuffer.data(idx_op_dg + i * 90 + 53);

        auto to_yy_xzzz = pbuffer.data(idx_op_dg + i * 90 + 54);

        auto to_yy_yyyy = pbuffer.data(idx_op_dg + i * 90 + 55);

        auto to_yy_yyyz = pbuffer.data(idx_op_dg + i * 90 + 56);

        auto to_yy_yyzz = pbuffer.data(idx_op_dg + i * 90 + 57);

        auto to_yy_yzzz = pbuffer.data(idx_op_dg + i * 90 + 58);

        auto to_yy_zzzz = pbuffer.data(idx_op_dg + i * 90 + 59);

        auto to_yz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 60);

        auto to_yz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 61);

        auto to_yz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 62);

        auto to_yz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 63);

        auto to_yz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 64);

        auto to_yz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 65);

        auto to_yz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 66);

        auto to_yz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 67);

        auto to_yz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 68);

        auto to_yz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 69);

        auto to_yz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 70);

        auto to_yz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 71);

        auto to_yz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 72);

        auto to_yz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 73);

        auto to_yz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 74);

        auto to_zz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 75);

        auto to_zz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 76);

        auto to_zz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 77);

        auto to_zz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 78);

        auto to_zz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 79);

        auto to_zz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 80);

        auto to_zz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 81);

        auto to_zz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 82);

        auto to_zz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 83);

        auto to_zz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 84);

        auto to_zz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 85);

        auto to_zz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 86);

        auto to_zz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 87);

        auto to_zz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 88);

        auto to_zz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 89);

        // Set up 0-10 components of targeted buffer : PF

        auto to_x_x_x_xxx = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 0);

        auto to_x_x_x_xxy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 1);

        auto to_x_x_x_xxz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 2);

        auto to_x_x_x_xyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 3);

        auto to_x_x_x_xyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 4);

        auto to_x_x_x_xzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 5);

        auto to_x_x_x_yyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 6);

        auto to_x_x_x_yyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 7);

        auto to_x_x_x_yzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 8);

        auto to_x_x_x_zzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_xx, to_0_xxxx, to_0_xxxy, to_0_xxxz, to_0_xxyy, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yz, to_0_zz, to_x_x_x_xxx, to_x_x_x_xxy, to_x_x_x_xxz, to_x_x_x_xyy, to_x_x_x_xyz, to_x_x_x_xzz, to_x_x_x_yyy, to_x_x_x_yyz, to_x_x_x_yzz, to_x_x_x_zzz, to_xx_xx, to_xx_xxxx, to_xx_xxxy, to_xx_xxxz, to_xx_xxyy, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_x_xxx[k] = 3.0 * to_0_xx[k] - 2.0 * to_0_xxxx[k] * tke_0 - 6.0 * to_xx_xx[k] * tbe_0 + 4.0 * to_xx_xxxx[k] * tbe_0 * tke_0;

            to_x_x_x_xxy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xxxy[k] * tke_0 - 4.0 * to_xx_xy[k] * tbe_0 + 4.0 * to_xx_xxxy[k] * tbe_0 * tke_0;

            to_x_x_x_xxz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xxxz[k] * tke_0 - 4.0 * to_xx_xz[k] * tbe_0 + 4.0 * to_xx_xxxz[k] * tbe_0 * tke_0;

            to_x_x_x_xyy[k] = to_0_yy[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_xx_yy[k] * tbe_0 + 4.0 * to_xx_xxyy[k] * tbe_0 * tke_0;

            to_x_x_x_xyz[k] = to_0_yz[k] - 2.0 * to_0_xxyz[k] * tke_0 - 2.0 * to_xx_yz[k] * tbe_0 + 4.0 * to_xx_xxyz[k] * tbe_0 * tke_0;

            to_x_x_x_xzz[k] = to_0_zz[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_xx_zz[k] * tbe_0 + 4.0 * to_xx_xxzz[k] * tbe_0 * tke_0;

            to_x_x_x_yyy[k] = -2.0 * to_0_xyyy[k] * tke_0 + 4.0 * to_xx_xyyy[k] * tbe_0 * tke_0;

            to_x_x_x_yyz[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_xx_xyyz[k] * tbe_0 * tke_0;

            to_x_x_x_yzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_xx_xyzz[k] * tbe_0 * tke_0;

            to_x_x_x_zzz[k] = -2.0 * to_0_xzzz[k] * tke_0 + 4.0 * to_xx_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : PF

        auto to_x_x_y_xxx = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 10);

        auto to_x_x_y_xxy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 11);

        auto to_x_x_y_xxz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 12);

        auto to_x_x_y_xyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 13);

        auto to_x_x_y_xyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 14);

        auto to_x_x_y_xzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 15);

        auto to_x_x_y_yyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 16);

        auto to_x_x_y_yyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 17);

        auto to_x_x_y_yzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 18);

        auto to_x_x_y_zzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_x_x_y_xxx, to_x_x_y_xxy, to_x_x_y_xxz, to_x_x_y_xyy, to_x_x_y_xyz, to_x_x_y_xzz, to_x_x_y_yyy, to_x_x_y_yyz, to_x_x_y_yzz, to_x_x_y_zzz, to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_y_xxx[k] = -6.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_y_xxy[k] = -4.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_y_xxz[k] = -4.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_y_xyy[k] = -2.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_y_xyz[k] = -2.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_y_xzz[k] = -2.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_y_yyy[k] = 4.0 * to_xy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_y_yyz[k] = 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_y_yzz[k] = 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_y_zzz[k] = 4.0 * to_xy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : PF

        auto to_x_x_z_xxx = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 20);

        auto to_x_x_z_xxy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 21);

        auto to_x_x_z_xxz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 22);

        auto to_x_x_z_xyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 23);

        auto to_x_x_z_xyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 24);

        auto to_x_x_z_xzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 25);

        auto to_x_x_z_yyy = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 26);

        auto to_x_x_z_yyz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 27);

        auto to_x_x_z_yzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 28);

        auto to_x_x_z_zzz = pbuffer.data(idx_op_geom_101_pf + 0 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_x_z_xxx, to_x_x_z_xxy, to_x_x_z_xxz, to_x_x_z_xyy, to_x_x_z_xyz, to_x_x_z_xzz, to_x_x_z_yyy, to_x_x_z_yyz, to_x_x_z_yzz, to_x_x_z_zzz, to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_z_xxx[k] = -6.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_z_xxy[k] = -4.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_z_xxz[k] = -4.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_z_xyy[k] = -2.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_z_xyz[k] = -2.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_z_xzz[k] = -2.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_z_yyy[k] = 4.0 * to_xz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_z_yyz[k] = 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_z_yzz[k] = 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_z_zzz[k] = 4.0 * to_xz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : PF

        auto to_x_y_x_xxx = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 0);

        auto to_x_y_x_xxy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 1);

        auto to_x_y_x_xxz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 2);

        auto to_x_y_x_xyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 3);

        auto to_x_y_x_xyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 4);

        auto to_x_y_x_xzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 5);

        auto to_x_y_x_yyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 6);

        auto to_x_y_x_yyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 7);

        auto to_x_y_x_yzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 8);

        auto to_x_y_x_zzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_xx, to_0_xxxy, to_0_xxyy, to_0_xxyz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_yy, to_0_yyyy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_x_y_x_xxx, to_x_y_x_xxy, to_x_y_x_xxz, to_x_y_x_xyy, to_x_y_x_xyz, to_x_y_x_xzz, to_x_y_x_yyy, to_x_y_x_yyz, to_x_y_x_yzz, to_x_y_x_zzz, to_xx_xx, to_xx_xxxy, to_xx_xxyy, to_xx_xxyz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_yy, to_xx_yyyy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_x_xxx[k] = -2.0 * to_0_xxxy[k] * tke_0 + 4.0 * to_xx_xxxy[k] * tbe_0 * tke_0;

            to_x_y_x_xxy[k] = to_0_xx[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_xx_xx[k] * tbe_0 + 4.0 * to_xx_xxyy[k] * tbe_0 * tke_0;

            to_x_y_x_xxz[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_xx_xxyz[k] * tbe_0 * tke_0;

            to_x_y_x_xyy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xyyy[k] * tke_0 - 4.0 * to_xx_xy[k] * tbe_0 + 4.0 * to_xx_xyyy[k] * tbe_0 * tke_0;

            to_x_y_x_xyz[k] = to_0_xz[k] - 2.0 * to_0_xyyz[k] * tke_0 - 2.0 * to_xx_xz[k] * tbe_0 + 4.0 * to_xx_xyyz[k] * tbe_0 * tke_0;

            to_x_y_x_xzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_xx_xyzz[k] * tbe_0 * tke_0;

            to_x_y_x_yyy[k] = 3.0 * to_0_yy[k] - 2.0 * to_0_yyyy[k] * tke_0 - 6.0 * to_xx_yy[k] * tbe_0 + 4.0 * to_xx_yyyy[k] * tbe_0 * tke_0;

            to_x_y_x_yyz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yyyz[k] * tke_0 - 4.0 * to_xx_yz[k] * tbe_0 + 4.0 * to_xx_yyyz[k] * tbe_0 * tke_0;

            to_x_y_x_yzz[k] = to_0_zz[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_xx_zz[k] * tbe_0 + 4.0 * to_xx_yyzz[k] * tbe_0 * tke_0;

            to_x_y_x_zzz[k] = -2.0 * to_0_yzzz[k] * tke_0 + 4.0 * to_xx_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : PF

        auto to_x_y_y_xxx = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 10);

        auto to_x_y_y_xxy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 11);

        auto to_x_y_y_xxz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 12);

        auto to_x_y_y_xyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 13);

        auto to_x_y_y_xyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 14);

        auto to_x_y_y_xzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 15);

        auto to_x_y_y_yyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 16);

        auto to_x_y_y_yyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 17);

        auto to_x_y_y_yzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 18);

        auto to_x_y_y_zzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_x_y_y_xxx, to_x_y_y_xxy, to_x_y_y_xxz, to_x_y_y_xyy, to_x_y_y_xyz, to_x_y_y_xzz, to_x_y_y_yyy, to_x_y_y_yyz, to_x_y_y_yzz, to_x_y_y_zzz, to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_y_xxx[k] = 4.0 * to_xy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_y_xxy[k] = -2.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_y_xxz[k] = 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_y_xyy[k] = -4.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_y_xyz[k] = -2.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_y_xzz[k] = 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_y_yyy[k] = -6.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_y_yyz[k] = -4.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_y_yzz[k] = -2.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_y_zzz[k] = 4.0 * to_xy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : PF

        auto to_x_y_z_xxx = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 20);

        auto to_x_y_z_xxy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 21);

        auto to_x_y_z_xxz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 22);

        auto to_x_y_z_xyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 23);

        auto to_x_y_z_xyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 24);

        auto to_x_y_z_xzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 25);

        auto to_x_y_z_yyy = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 26);

        auto to_x_y_z_yyz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 27);

        auto to_x_y_z_yzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 28);

        auto to_x_y_z_zzz = pbuffer.data(idx_op_geom_101_pf + 1 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_y_z_xxx, to_x_y_z_xxy, to_x_y_z_xxz, to_x_y_z_xyy, to_x_y_z_xyz, to_x_y_z_xzz, to_x_y_z_yyy, to_x_y_z_yyz, to_x_y_z_yzz, to_x_y_z_zzz, to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_z_xxx[k] = 4.0 * to_xz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_z_xxy[k] = -2.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_z_xxz[k] = 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_z_xyy[k] = -4.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_z_xyz[k] = -2.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_z_xzz[k] = 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_z_yyy[k] = -6.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_z_yyz[k] = -4.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_z_yzz[k] = -2.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_z_zzz[k] = 4.0 * to_xz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : PF

        auto to_x_z_x_xxx = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 0);

        auto to_x_z_x_xxy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 1);

        auto to_x_z_x_xxz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 2);

        auto to_x_z_x_xyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 3);

        auto to_x_z_x_xyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 4);

        auto to_x_z_x_xzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 5);

        auto to_x_z_x_yyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 6);

        auto to_x_z_x_yyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 7);

        auto to_x_z_x_yzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 8);

        auto to_x_z_x_zzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_xx, to_0_xxxz, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_0_zzzz, to_x_z_x_xxx, to_x_z_x_xxy, to_x_z_x_xxz, to_x_z_x_xyy, to_x_z_x_xyz, to_x_z_x_xzz, to_x_z_x_yyy, to_x_z_x_yyz, to_x_z_x_yzz, to_x_z_x_zzz, to_xx_xx, to_xx_xxxz, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_x_xxx[k] = -2.0 * to_0_xxxz[k] * tke_0 + 4.0 * to_xx_xxxz[k] * tbe_0 * tke_0;

            to_x_z_x_xxy[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_xx_xxyz[k] * tbe_0 * tke_0;

            to_x_z_x_xxz[k] = to_0_xx[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_xx_xx[k] * tbe_0 + 4.0 * to_xx_xxzz[k] * tbe_0 * tke_0;

            to_x_z_x_xyy[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_xx_xyyz[k] * tbe_0 * tke_0;

            to_x_z_x_xyz[k] = to_0_xy[k] - 2.0 * to_0_xyzz[k] * tke_0 - 2.0 * to_xx_xy[k] * tbe_0 + 4.0 * to_xx_xyzz[k] * tbe_0 * tke_0;

            to_x_z_x_xzz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xzzz[k] * tke_0 - 4.0 * to_xx_xz[k] * tbe_0 + 4.0 * to_xx_xzzz[k] * tbe_0 * tke_0;

            to_x_z_x_yyy[k] = -2.0 * to_0_yyyz[k] * tke_0 + 4.0 * to_xx_yyyz[k] * tbe_0 * tke_0;

            to_x_z_x_yyz[k] = to_0_yy[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_xx_yy[k] * tbe_0 + 4.0 * to_xx_yyzz[k] * tbe_0 * tke_0;

            to_x_z_x_yzz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yzzz[k] * tke_0 - 4.0 * to_xx_yz[k] * tbe_0 + 4.0 * to_xx_yzzz[k] * tbe_0 * tke_0;

            to_x_z_x_zzz[k] = 3.0 * to_0_zz[k] - 2.0 * to_0_zzzz[k] * tke_0 - 6.0 * to_xx_zz[k] * tbe_0 + 4.0 * to_xx_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : PF

        auto to_x_z_y_xxx = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 10);

        auto to_x_z_y_xxy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 11);

        auto to_x_z_y_xxz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 12);

        auto to_x_z_y_xyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 13);

        auto to_x_z_y_xyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 14);

        auto to_x_z_y_xzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 15);

        auto to_x_z_y_yyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 16);

        auto to_x_z_y_yyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 17);

        auto to_x_z_y_yzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 18);

        auto to_x_z_y_zzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_x_z_y_xxx, to_x_z_y_xxy, to_x_z_y_xxz, to_x_z_y_xyy, to_x_z_y_xyz, to_x_z_y_xzz, to_x_z_y_yyy, to_x_z_y_yyz, to_x_z_y_yzz, to_x_z_y_zzz, to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_y_xxx[k] = 4.0 * to_xy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_y_xxy[k] = 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_y_xxz[k] = -2.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_y_xyy[k] = 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_y_xyz[k] = -2.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_y_xzz[k] = -4.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_y_yyy[k] = 4.0 * to_xy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_y_yyz[k] = -2.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_y_yzz[k] = -4.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_y_zzz[k] = -6.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : PF

        auto to_x_z_z_xxx = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 20);

        auto to_x_z_z_xxy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 21);

        auto to_x_z_z_xxz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 22);

        auto to_x_z_z_xyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 23);

        auto to_x_z_z_xyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 24);

        auto to_x_z_z_xzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 25);

        auto to_x_z_z_yyy = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 26);

        auto to_x_z_z_yyz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 27);

        auto to_x_z_z_yzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 28);

        auto to_x_z_z_zzz = pbuffer.data(idx_op_geom_101_pf + 2 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_z_z_xxx, to_x_z_z_xxy, to_x_z_z_xxz, to_x_z_z_xyy, to_x_z_z_xyz, to_x_z_z_xzz, to_x_z_z_yyy, to_x_z_z_yyz, to_x_z_z_yzz, to_x_z_z_zzz, to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_z_xxx[k] = 4.0 * to_xz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_z_xxy[k] = 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_z_xxz[k] = -2.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_z_xyy[k] = 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_z_xyz[k] = -2.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_z_xzz[k] = -4.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_z_yyy[k] = 4.0 * to_xz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_z_yyz[k] = -2.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_z_yzz[k] = -4.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_z_zzz[k] = -6.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-100 components of targeted buffer : PF

        auto to_y_x_x_xxx = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 0);

        auto to_y_x_x_xxy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 1);

        auto to_y_x_x_xxz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 2);

        auto to_y_x_x_xyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 3);

        auto to_y_x_x_xyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 4);

        auto to_y_x_x_xzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 5);

        auto to_y_x_x_yyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 6);

        auto to_y_x_x_yyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 7);

        auto to_y_x_x_yzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 8);

        auto to_y_x_x_zzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, to_y_x_x_xxx, to_y_x_x_xxy, to_y_x_x_xxz, to_y_x_x_xyy, to_y_x_x_xyz, to_y_x_x_xzz, to_y_x_x_yyy, to_y_x_x_yyz, to_y_x_x_yzz, to_y_x_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_x_xxx[k] = -6.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_x_xxy[k] = -4.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_x_xxz[k] = -4.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_x_xyy[k] = -2.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_x_xyz[k] = -2.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_x_xzz[k] = -2.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_x_yyy[k] = 4.0 * to_xy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_x_yyz[k] = 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_x_yzz[k] = 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_x_zzz[k] = 4.0 * to_xy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 100-110 components of targeted buffer : PF

        auto to_y_x_y_xxx = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 10);

        auto to_y_x_y_xxy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 11);

        auto to_y_x_y_xxz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 12);

        auto to_y_x_y_xyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 13);

        auto to_y_x_y_xyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 14);

        auto to_y_x_y_xzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 15);

        auto to_y_x_y_yyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 16);

        auto to_y_x_y_yyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 17);

        auto to_y_x_y_yzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 18);

        auto to_y_x_y_zzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_xx, to_0_xxxx, to_0_xxxy, to_0_xxxz, to_0_xxyy, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yz, to_0_zz, to_y_x_y_xxx, to_y_x_y_xxy, to_y_x_y_xxz, to_y_x_y_xyy, to_y_x_y_xyz, to_y_x_y_xzz, to_y_x_y_yyy, to_y_x_y_yyz, to_y_x_y_yzz, to_y_x_y_zzz, to_yy_xx, to_yy_xxxx, to_yy_xxxy, to_yy_xxxz, to_yy_xxyy, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_y_xxx[k] = 3.0 * to_0_xx[k] - 2.0 * to_0_xxxx[k] * tke_0 - 6.0 * to_yy_xx[k] * tbe_0 + 4.0 * to_yy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_y_xxy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xxxy[k] * tke_0 - 4.0 * to_yy_xy[k] * tbe_0 + 4.0 * to_yy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_y_xxz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xxxz[k] * tke_0 - 4.0 * to_yy_xz[k] * tbe_0 + 4.0 * to_yy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_y_xyy[k] = to_0_yy[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_yy_yy[k] * tbe_0 + 4.0 * to_yy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_y_xyz[k] = to_0_yz[k] - 2.0 * to_0_xxyz[k] * tke_0 - 2.0 * to_yy_yz[k] * tbe_0 + 4.0 * to_yy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_y_xzz[k] = to_0_zz[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_yy_zz[k] * tbe_0 + 4.0 * to_yy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_y_yyy[k] = -2.0 * to_0_xyyy[k] * tke_0 + 4.0 * to_yy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_y_yyz[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_yy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_y_yzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_yy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_y_zzz[k] = -2.0 * to_0_xzzz[k] * tke_0 + 4.0 * to_yy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 110-120 components of targeted buffer : PF

        auto to_y_x_z_xxx = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 20);

        auto to_y_x_z_xxy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 21);

        auto to_y_x_z_xxz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 22);

        auto to_y_x_z_xyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 23);

        auto to_y_x_z_xyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 24);

        auto to_y_x_z_xzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 25);

        auto to_y_x_z_yyy = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 26);

        auto to_y_x_z_yyz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 27);

        auto to_y_x_z_yzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 28);

        auto to_y_x_z_zzz = pbuffer.data(idx_op_geom_101_pf + 3 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_x_z_xxx, to_y_x_z_xxy, to_y_x_z_xxz, to_y_x_z_xyy, to_y_x_z_xyz, to_y_x_z_xzz, to_y_x_z_yyy, to_y_x_z_yyz, to_y_x_z_yzz, to_y_x_z_zzz, to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_z_xxx[k] = -6.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_z_xxy[k] = -4.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_z_xxz[k] = -4.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_z_xyy[k] = -2.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_z_xyz[k] = -2.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_z_xzz[k] = -2.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_z_yyy[k] = 4.0 * to_yz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_z_yyz[k] = 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_z_yzz[k] = 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_z_zzz[k] = 4.0 * to_yz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-130 components of targeted buffer : PF

        auto to_y_y_x_xxx = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 0);

        auto to_y_y_x_xxy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 1);

        auto to_y_y_x_xxz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 2);

        auto to_y_y_x_xyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 3);

        auto to_y_y_x_xyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 4);

        auto to_y_y_x_xzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 5);

        auto to_y_y_x_yyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 6);

        auto to_y_y_x_yyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 7);

        auto to_y_y_x_yzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 8);

        auto to_y_y_x_zzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_y_y_x_xxx, to_y_y_x_xxy, to_y_y_x_xxz, to_y_y_x_xyy, to_y_y_x_xyz, to_y_y_x_xzz, to_y_y_x_yyy, to_y_y_x_yyz, to_y_y_x_yzz, to_y_y_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_x_xxx[k] = 4.0 * to_xy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_x_xxy[k] = -2.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_x_xxz[k] = 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_x_xyy[k] = -4.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_x_xyz[k] = -2.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_x_xzz[k] = 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_x_yyy[k] = -6.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_x_yyz[k] = -4.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_x_yzz[k] = -2.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_x_zzz[k] = 4.0 * to_xy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 130-140 components of targeted buffer : PF

        auto to_y_y_y_xxx = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 10);

        auto to_y_y_y_xxy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 11);

        auto to_y_y_y_xxz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 12);

        auto to_y_y_y_xyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 13);

        auto to_y_y_y_xyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 14);

        auto to_y_y_y_xzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 15);

        auto to_y_y_y_yyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 16);

        auto to_y_y_y_yyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 17);

        auto to_y_y_y_yzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 18);

        auto to_y_y_y_zzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_xx, to_0_xxxy, to_0_xxyy, to_0_xxyz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_yy, to_0_yyyy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_y_y_y_xxx, to_y_y_y_xxy, to_y_y_y_xxz, to_y_y_y_xyy, to_y_y_y_xyz, to_y_y_y_xzz, to_y_y_y_yyy, to_y_y_y_yyz, to_y_y_y_yzz, to_y_y_y_zzz, to_yy_xx, to_yy_xxxy, to_yy_xxyy, to_yy_xxyz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_yy, to_yy_yyyy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_y_xxx[k] = -2.0 * to_0_xxxy[k] * tke_0 + 4.0 * to_yy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_y_xxy[k] = to_0_xx[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_yy_xx[k] * tbe_0 + 4.0 * to_yy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_y_xxz[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_yy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_y_xyy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xyyy[k] * tke_0 - 4.0 * to_yy_xy[k] * tbe_0 + 4.0 * to_yy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_y_xyz[k] = to_0_xz[k] - 2.0 * to_0_xyyz[k] * tke_0 - 2.0 * to_yy_xz[k] * tbe_0 + 4.0 * to_yy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_y_xzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_yy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_y_yyy[k] = 3.0 * to_0_yy[k] - 2.0 * to_0_yyyy[k] * tke_0 - 6.0 * to_yy_yy[k] * tbe_0 + 4.0 * to_yy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_y_yyz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yyyz[k] * tke_0 - 4.0 * to_yy_yz[k] * tbe_0 + 4.0 * to_yy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_y_yzz[k] = to_0_zz[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_yy_zz[k] * tbe_0 + 4.0 * to_yy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_y_zzz[k] = -2.0 * to_0_yzzz[k] * tke_0 + 4.0 * to_yy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 140-150 components of targeted buffer : PF

        auto to_y_y_z_xxx = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 20);

        auto to_y_y_z_xxy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 21);

        auto to_y_y_z_xxz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 22);

        auto to_y_y_z_xyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 23);

        auto to_y_y_z_xyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 24);

        auto to_y_y_z_xzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 25);

        auto to_y_y_z_yyy = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 26);

        auto to_y_y_z_yyz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 27);

        auto to_y_y_z_yzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 28);

        auto to_y_y_z_zzz = pbuffer.data(idx_op_geom_101_pf + 4 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_y_z_xxx, to_y_y_z_xxy, to_y_y_z_xxz, to_y_y_z_xyy, to_y_y_z_xyz, to_y_y_z_xzz, to_y_y_z_yyy, to_y_y_z_yyz, to_y_y_z_yzz, to_y_y_z_zzz, to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_z_xxx[k] = 4.0 * to_yz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_z_xxy[k] = -2.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_z_xxz[k] = 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_z_xyy[k] = -4.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_z_xyz[k] = -2.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_z_xzz[k] = 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_z_yyy[k] = -6.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_z_yyz[k] = -4.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_z_yzz[k] = -2.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_z_zzz[k] = 4.0 * to_yz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-160 components of targeted buffer : PF

        auto to_y_z_x_xxx = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 0);

        auto to_y_z_x_xxy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 1);

        auto to_y_z_x_xxz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 2);

        auto to_y_z_x_xyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 3);

        auto to_y_z_x_xyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 4);

        auto to_y_z_x_xzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 5);

        auto to_y_z_x_yyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 6);

        auto to_y_z_x_yyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 7);

        auto to_y_z_x_yzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 8);

        auto to_y_z_x_zzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, to_y_z_x_xxx, to_y_z_x_xxy, to_y_z_x_xxz, to_y_z_x_xyy, to_y_z_x_xyz, to_y_z_x_xzz, to_y_z_x_yyy, to_y_z_x_yyz, to_y_z_x_yzz, to_y_z_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_x_xxx[k] = 4.0 * to_xy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_x_xxy[k] = 4.0 * to_xy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_x_xxz[k] = -2.0 * to_xy_xx[k] * tbe_0 + 4.0 * to_xy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_x_xyy[k] = 4.0 * to_xy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_x_xyz[k] = -2.0 * to_xy_xy[k] * tbe_0 + 4.0 * to_xy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_x_xzz[k] = -4.0 * to_xy_xz[k] * tbe_0 + 4.0 * to_xy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_x_yyy[k] = 4.0 * to_xy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_x_yyz[k] = -2.0 * to_xy_yy[k] * tbe_0 + 4.0 * to_xy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_x_yzz[k] = -4.0 * to_xy_yz[k] * tbe_0 + 4.0 * to_xy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_x_zzz[k] = -6.0 * to_xy_zz[k] * tbe_0 + 4.0 * to_xy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 160-170 components of targeted buffer : PF

        auto to_y_z_y_xxx = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 10);

        auto to_y_z_y_xxy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 11);

        auto to_y_z_y_xxz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 12);

        auto to_y_z_y_xyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 13);

        auto to_y_z_y_xyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 14);

        auto to_y_z_y_xzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 15);

        auto to_y_z_y_yyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 16);

        auto to_y_z_y_yyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 17);

        auto to_y_z_y_yzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 18);

        auto to_y_z_y_zzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_xx, to_0_xxxz, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_0_zzzz, to_y_z_y_xxx, to_y_z_y_xxy, to_y_z_y_xxz, to_y_z_y_xyy, to_y_z_y_xyz, to_y_z_y_xzz, to_y_z_y_yyy, to_y_z_y_yyz, to_y_z_y_yzz, to_y_z_y_zzz, to_yy_xx, to_yy_xxxz, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_y_xxx[k] = -2.0 * to_0_xxxz[k] * tke_0 + 4.0 * to_yy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_y_xxy[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_yy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_y_xxz[k] = to_0_xx[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_yy_xx[k] * tbe_0 + 4.0 * to_yy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_y_xyy[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_yy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_y_xyz[k] = to_0_xy[k] - 2.0 * to_0_xyzz[k] * tke_0 - 2.0 * to_yy_xy[k] * tbe_0 + 4.0 * to_yy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_y_xzz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xzzz[k] * tke_0 - 4.0 * to_yy_xz[k] * tbe_0 + 4.0 * to_yy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_y_yyy[k] = -2.0 * to_0_yyyz[k] * tke_0 + 4.0 * to_yy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_y_yyz[k] = to_0_yy[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_yy_yy[k] * tbe_0 + 4.0 * to_yy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_y_yzz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yzzz[k] * tke_0 - 4.0 * to_yy_yz[k] * tbe_0 + 4.0 * to_yy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_y_zzz[k] = 3.0 * to_0_zz[k] - 2.0 * to_0_zzzz[k] * tke_0 - 6.0 * to_yy_zz[k] * tbe_0 + 4.0 * to_yy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 170-180 components of targeted buffer : PF

        auto to_y_z_z_xxx = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 20);

        auto to_y_z_z_xxy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 21);

        auto to_y_z_z_xxz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 22);

        auto to_y_z_z_xyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 23);

        auto to_y_z_z_xyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 24);

        auto to_y_z_z_xzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 25);

        auto to_y_z_z_yyy = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 26);

        auto to_y_z_z_yyz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 27);

        auto to_y_z_z_yzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 28);

        auto to_y_z_z_zzz = pbuffer.data(idx_op_geom_101_pf + 5 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_z_z_xxx, to_y_z_z_xxy, to_y_z_z_xxz, to_y_z_z_xyy, to_y_z_z_xyz, to_y_z_z_xzz, to_y_z_z_yyy, to_y_z_z_yyz, to_y_z_z_yzz, to_y_z_z_zzz, to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_z_xxx[k] = 4.0 * to_yz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_z_xxy[k] = 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_z_xxz[k] = -2.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_z_xyy[k] = 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_z_xyz[k] = -2.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_z_xzz[k] = -4.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_z_yyy[k] = 4.0 * to_yz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_z_yyz[k] = -2.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_z_yzz[k] = -4.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_z_zzz[k] = -6.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-190 components of targeted buffer : PF

        auto to_z_x_x_xxx = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 0);

        auto to_z_x_x_xxy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 1);

        auto to_z_x_x_xxz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 2);

        auto to_z_x_x_xyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 3);

        auto to_z_x_x_xyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 4);

        auto to_z_x_x_xzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 5);

        auto to_z_x_x_yyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 6);

        auto to_z_x_x_yyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 7);

        auto to_z_x_x_yzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 8);

        auto to_z_x_x_zzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, to_z_x_x_xxx, to_z_x_x_xxy, to_z_x_x_xxz, to_z_x_x_xyy, to_z_x_x_xyz, to_z_x_x_xzz, to_z_x_x_yyy, to_z_x_x_yyz, to_z_x_x_yzz, to_z_x_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_x_xxx[k] = -6.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_x_xxy[k] = -4.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_x_xxz[k] = -4.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_x_xyy[k] = -2.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_x_xyz[k] = -2.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_x_xzz[k] = -2.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_x_yyy[k] = 4.0 * to_xz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_x_yyz[k] = 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_x_yzz[k] = 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_x_zzz[k] = 4.0 * to_xz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 190-200 components of targeted buffer : PF

        auto to_z_x_y_xxx = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 10);

        auto to_z_x_y_xxy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 11);

        auto to_z_x_y_xxz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 12);

        auto to_z_x_y_xyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 13);

        auto to_z_x_y_xyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 14);

        auto to_z_x_y_xzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 15);

        auto to_z_x_y_yyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 16);

        auto to_z_x_y_yyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 17);

        auto to_z_x_y_yzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 18);

        auto to_z_x_y_zzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, to_z_x_y_xxx, to_z_x_y_xxy, to_z_x_y_xxz, to_z_x_y_xyy, to_z_x_y_xyz, to_z_x_y_xzz, to_z_x_y_yyy, to_z_x_y_yyz, to_z_x_y_yzz, to_z_x_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_y_xxx[k] = -6.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_y_xxy[k] = -4.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_y_xxz[k] = -4.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_y_xyy[k] = -2.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_y_xyz[k] = -2.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_y_xzz[k] = -2.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_y_yyy[k] = 4.0 * to_yz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_y_yyz[k] = 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_y_yzz[k] = 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_y_zzz[k] = 4.0 * to_yz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 200-210 components of targeted buffer : PF

        auto to_z_x_z_xxx = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 20);

        auto to_z_x_z_xxy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 21);

        auto to_z_x_z_xxz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 22);

        auto to_z_x_z_xyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 23);

        auto to_z_x_z_xyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 24);

        auto to_z_x_z_xzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 25);

        auto to_z_x_z_yyy = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 26);

        auto to_z_x_z_yyz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 27);

        auto to_z_x_z_yzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 28);

        auto to_z_x_z_zzz = pbuffer.data(idx_op_geom_101_pf + 6 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_xx, to_0_xxxx, to_0_xxxy, to_0_xxxz, to_0_xxyy, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yz, to_0_zz, to_z_x_z_xxx, to_z_x_z_xxy, to_z_x_z_xxz, to_z_x_z_xyy, to_z_x_z_xyz, to_z_x_z_xzz, to_z_x_z_yyy, to_z_x_z_yyz, to_z_x_z_yzz, to_z_x_z_zzz, to_zz_xx, to_zz_xxxx, to_zz_xxxy, to_zz_xxxz, to_zz_xxyy, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_z_xxx[k] = 3.0 * to_0_xx[k] - 2.0 * to_0_xxxx[k] * tke_0 - 6.0 * to_zz_xx[k] * tbe_0 + 4.0 * to_zz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_z_xxy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xxxy[k] * tke_0 - 4.0 * to_zz_xy[k] * tbe_0 + 4.0 * to_zz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_z_xxz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xxxz[k] * tke_0 - 4.0 * to_zz_xz[k] * tbe_0 + 4.0 * to_zz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_z_xyy[k] = to_0_yy[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_zz_yy[k] * tbe_0 + 4.0 * to_zz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_z_xyz[k] = to_0_yz[k] - 2.0 * to_0_xxyz[k] * tke_0 - 2.0 * to_zz_yz[k] * tbe_0 + 4.0 * to_zz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_z_xzz[k] = to_0_zz[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_zz_zz[k] * tbe_0 + 4.0 * to_zz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_z_yyy[k] = -2.0 * to_0_xyyy[k] * tke_0 + 4.0 * to_zz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_z_yyz[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_zz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_z_yzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_zz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_z_zzz[k] = -2.0 * to_0_xzzz[k] * tke_0 + 4.0 * to_zz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-220 components of targeted buffer : PF

        auto to_z_y_x_xxx = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 0);

        auto to_z_y_x_xxy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 1);

        auto to_z_y_x_xxz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 2);

        auto to_z_y_x_xyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 3);

        auto to_z_y_x_xyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 4);

        auto to_z_y_x_xzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 5);

        auto to_z_y_x_yyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 6);

        auto to_z_y_x_yyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 7);

        auto to_z_y_x_yzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 8);

        auto to_z_y_x_zzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_z_y_x_xxx, to_z_y_x_xxy, to_z_y_x_xxz, to_z_y_x_xyy, to_z_y_x_xyz, to_z_y_x_xzz, to_z_y_x_yyy, to_z_y_x_yyz, to_z_y_x_yzz, to_z_y_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_x_xxx[k] = 4.0 * to_xz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_x_xxy[k] = -2.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_x_xxz[k] = 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_x_xyy[k] = -4.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_x_xyz[k] = -2.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_x_xzz[k] = 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_x_yyy[k] = -6.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_x_yyz[k] = -4.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_x_yzz[k] = -2.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_x_zzz[k] = 4.0 * to_xz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 220-230 components of targeted buffer : PF

        auto to_z_y_y_xxx = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 10);

        auto to_z_y_y_xxy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 11);

        auto to_z_y_y_xxz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 12);

        auto to_z_y_y_xyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 13);

        auto to_z_y_y_xyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 14);

        auto to_z_y_y_xzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 15);

        auto to_z_y_y_yyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 16);

        auto to_z_y_y_yyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 17);

        auto to_z_y_y_yzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 18);

        auto to_z_y_y_zzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_z_y_y_xxx, to_z_y_y_xxy, to_z_y_y_xxz, to_z_y_y_xyy, to_z_y_y_xyz, to_z_y_y_xzz, to_z_y_y_yyy, to_z_y_y_yyz, to_z_y_y_yzz, to_z_y_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_y_xxx[k] = 4.0 * to_yz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_y_xxy[k] = -2.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_y_xxz[k] = 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_y_xyy[k] = -4.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_y_xyz[k] = -2.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_y_xzz[k] = 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_y_yyy[k] = -6.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_y_yyz[k] = -4.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_y_yzz[k] = -2.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_y_zzz[k] = 4.0 * to_yz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 230-240 components of targeted buffer : PF

        auto to_z_y_z_xxx = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 20);

        auto to_z_y_z_xxy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 21);

        auto to_z_y_z_xxz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 22);

        auto to_z_y_z_xyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 23);

        auto to_z_y_z_xyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 24);

        auto to_z_y_z_xzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 25);

        auto to_z_y_z_yyy = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 26);

        auto to_z_y_z_yyz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 27);

        auto to_z_y_z_yzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 28);

        auto to_z_y_z_zzz = pbuffer.data(idx_op_geom_101_pf + 7 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_xx, to_0_xxxy, to_0_xxyy, to_0_xxyz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_yy, to_0_yyyy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_z_y_z_xxx, to_z_y_z_xxy, to_z_y_z_xxz, to_z_y_z_xyy, to_z_y_z_xyz, to_z_y_z_xzz, to_z_y_z_yyy, to_z_y_z_yyz, to_z_y_z_yzz, to_z_y_z_zzz, to_zz_xx, to_zz_xxxy, to_zz_xxyy, to_zz_xxyz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_yy, to_zz_yyyy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_z_xxx[k] = -2.0 * to_0_xxxy[k] * tke_0 + 4.0 * to_zz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_z_xxy[k] = to_0_xx[k] - 2.0 * to_0_xxyy[k] * tke_0 - 2.0 * to_zz_xx[k] * tbe_0 + 4.0 * to_zz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_z_xxz[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_zz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_z_xyy[k] = 2.0 * to_0_xy[k] - 2.0 * to_0_xyyy[k] * tke_0 - 4.0 * to_zz_xy[k] * tbe_0 + 4.0 * to_zz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_z_xyz[k] = to_0_xz[k] - 2.0 * to_0_xyyz[k] * tke_0 - 2.0 * to_zz_xz[k] * tbe_0 + 4.0 * to_zz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_z_xzz[k] = -2.0 * to_0_xyzz[k] * tke_0 + 4.0 * to_zz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_z_yyy[k] = 3.0 * to_0_yy[k] - 2.0 * to_0_yyyy[k] * tke_0 - 6.0 * to_zz_yy[k] * tbe_0 + 4.0 * to_zz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_z_yyz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yyyz[k] * tke_0 - 4.0 * to_zz_yz[k] * tbe_0 + 4.0 * to_zz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_z_yzz[k] = to_0_zz[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_zz_zz[k] * tbe_0 + 4.0 * to_zz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_z_zzz[k] = -2.0 * to_0_yzzz[k] * tke_0 + 4.0 * to_zz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-250 components of targeted buffer : PF

        auto to_z_z_x_xxx = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 0);

        auto to_z_z_x_xxy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 1);

        auto to_z_z_x_xxz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 2);

        auto to_z_z_x_xyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 3);

        auto to_z_z_x_xyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 4);

        auto to_z_z_x_xzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 5);

        auto to_z_z_x_yyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 6);

        auto to_z_z_x_yyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 7);

        auto to_z_z_x_yzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 8);

        auto to_z_z_x_zzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, to_z_z_x_xxx, to_z_z_x_xxy, to_z_z_x_xxz, to_z_z_x_xyy, to_z_z_x_xyz, to_z_z_x_xzz, to_z_z_x_yyy, to_z_z_x_yyz, to_z_z_x_yzz, to_z_z_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_x_xxx[k] = 4.0 * to_xz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_x_xxy[k] = 4.0 * to_xz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_x_xxz[k] = -2.0 * to_xz_xx[k] * tbe_0 + 4.0 * to_xz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_x_xyy[k] = 4.0 * to_xz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_x_xyz[k] = -2.0 * to_xz_xy[k] * tbe_0 + 4.0 * to_xz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_x_xzz[k] = -4.0 * to_xz_xz[k] * tbe_0 + 4.0 * to_xz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_x_yyy[k] = 4.0 * to_xz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_x_yyz[k] = -2.0 * to_xz_yy[k] * tbe_0 + 4.0 * to_xz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_x_yzz[k] = -4.0 * to_xz_yz[k] * tbe_0 + 4.0 * to_xz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_x_zzz[k] = -6.0 * to_xz_zz[k] * tbe_0 + 4.0 * to_xz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 250-260 components of targeted buffer : PF

        auto to_z_z_y_xxx = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 10);

        auto to_z_z_y_xxy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 11);

        auto to_z_z_y_xxz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 12);

        auto to_z_z_y_xyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 13);

        auto to_z_z_y_xyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 14);

        auto to_z_z_y_xzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 15);

        auto to_z_z_y_yyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 16);

        auto to_z_z_y_yyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 17);

        auto to_z_z_y_yzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 18);

        auto to_z_z_y_zzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, to_z_z_y_xxx, to_z_z_y_xxy, to_z_z_y_xxz, to_z_z_y_xyy, to_z_z_y_xyz, to_z_z_y_xzz, to_z_z_y_yyy, to_z_z_y_yyz, to_z_z_y_yzz, to_z_z_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_y_xxx[k] = 4.0 * to_yz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_y_xxy[k] = 4.0 * to_yz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_y_xxz[k] = -2.0 * to_yz_xx[k] * tbe_0 + 4.0 * to_yz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_y_xyy[k] = 4.0 * to_yz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_y_xyz[k] = -2.0 * to_yz_xy[k] * tbe_0 + 4.0 * to_yz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_y_xzz[k] = -4.0 * to_yz_xz[k] * tbe_0 + 4.0 * to_yz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_y_yyy[k] = 4.0 * to_yz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_y_yyz[k] = -2.0 * to_yz_yy[k] * tbe_0 + 4.0 * to_yz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_y_yzz[k] = -4.0 * to_yz_yz[k] * tbe_0 + 4.0 * to_yz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_y_zzz[k] = -6.0 * to_yz_zz[k] * tbe_0 + 4.0 * to_yz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 260-270 components of targeted buffer : PF

        auto to_z_z_z_xxx = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 20);

        auto to_z_z_z_xxy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 21);

        auto to_z_z_z_xxz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 22);

        auto to_z_z_z_xyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 23);

        auto to_z_z_z_xyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 24);

        auto to_z_z_z_xzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 25);

        auto to_z_z_z_yyy = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 26);

        auto to_z_z_z_yyz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 27);

        auto to_z_z_z_yzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 28);

        auto to_z_z_z_zzz = pbuffer.data(idx_op_geom_101_pf + 8 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_xx, to_0_xxxz, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, to_0_zzzz, to_z_z_z_xxx, to_z_z_z_xxy, to_z_z_z_xxz, to_z_z_z_xyy, to_z_z_z_xyz, to_z_z_z_xzz, to_z_z_z_yyy, to_z_z_z_yyz, to_z_z_z_yzz, to_z_z_z_zzz, to_zz_xx, to_zz_xxxz, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_z_xxx[k] = -2.0 * to_0_xxxz[k] * tke_0 + 4.0 * to_zz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_z_xxy[k] = -2.0 * to_0_xxyz[k] * tke_0 + 4.0 * to_zz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_z_xxz[k] = to_0_xx[k] - 2.0 * to_0_xxzz[k] * tke_0 - 2.0 * to_zz_xx[k] * tbe_0 + 4.0 * to_zz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_z_xyy[k] = -2.0 * to_0_xyyz[k] * tke_0 + 4.0 * to_zz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_z_xyz[k] = to_0_xy[k] - 2.0 * to_0_xyzz[k] * tke_0 - 2.0 * to_zz_xy[k] * tbe_0 + 4.0 * to_zz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_z_xzz[k] = 2.0 * to_0_xz[k] - 2.0 * to_0_xzzz[k] * tke_0 - 4.0 * to_zz_xz[k] * tbe_0 + 4.0 * to_zz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_z_yyy[k] = -2.0 * to_0_yyyz[k] * tke_0 + 4.0 * to_zz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_z_yyz[k] = to_0_yy[k] - 2.0 * to_0_yyzz[k] * tke_0 - 2.0 * to_zz_yy[k] * tbe_0 + 4.0 * to_zz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_z_yzz[k] = 2.0 * to_0_yz[k] - 2.0 * to_0_yzzz[k] * tke_0 - 4.0 * to_zz_yz[k] * tbe_0 + 4.0 * to_zz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_z_zzz[k] = 3.0 * to_0_zz[k] - 2.0 * to_0_zzzz[k] * tke_0 - 6.0 * to_zz_zz[k] * tbe_0 + 4.0 * to_zz_zzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

