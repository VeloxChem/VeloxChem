#include "GeometricalDerivatives0X1ForDF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_df(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_df,
                       const int idx_op_dd,
                       const int idx_op_dg,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up 0-10 components of targeted buffer : DF

        auto to_0_x_xx_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 0);

        auto to_0_x_xx_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 1);

        auto to_0_x_xx_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 2);

        auto to_0_x_xx_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 3);

        auto to_0_x_xx_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 4);

        auto to_0_x_xx_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 5);

        auto to_0_x_xx_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 6);

        auto to_0_x_xx_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 7);

        auto to_0_x_xx_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 8);

        auto to_0_x_xx_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_0_x_xx_xxx, to_0_x_xx_xxy, to_0_x_xx_xxz, to_0_x_xx_xyy, to_0_x_xx_xyz, to_0_x_xx_xzz, to_0_x_xx_yyy, to_0_x_xx_yyz, to_0_x_xx_yzz, to_0_x_xx_zzz, to_xx_xx, to_xx_xxxx, to_xx_xxxy, to_xx_xxxz, to_xx_xxyy, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xx_xxx[k] = -3.0 * to_xx_xx[k] + 2.0 * to_xx_xxxx[k] * tke_0;

            to_0_x_xx_xxy[k] = -2.0 * to_xx_xy[k] + 2.0 * to_xx_xxxy[k] * tke_0;

            to_0_x_xx_xxz[k] = -2.0 * to_xx_xz[k] + 2.0 * to_xx_xxxz[k] * tke_0;

            to_0_x_xx_xyy[k] = -to_xx_yy[k] + 2.0 * to_xx_xxyy[k] * tke_0;

            to_0_x_xx_xyz[k] = -to_xx_yz[k] + 2.0 * to_xx_xxyz[k] * tke_0;

            to_0_x_xx_xzz[k] = -to_xx_zz[k] + 2.0 * to_xx_xxzz[k] * tke_0;

            to_0_x_xx_yyy[k] = 2.0 * to_xx_xyyy[k] * tke_0;

            to_0_x_xx_yyz[k] = 2.0 * to_xx_xyyz[k] * tke_0;

            to_0_x_xx_yzz[k] = 2.0 * to_xx_xyzz[k] * tke_0;

            to_0_x_xx_zzz[k] = 2.0 * to_xx_xzzz[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : DF

        auto to_0_x_xy_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 10);

        auto to_0_x_xy_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 11);

        auto to_0_x_xy_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 12);

        auto to_0_x_xy_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 13);

        auto to_0_x_xy_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 14);

        auto to_0_x_xy_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 15);

        auto to_0_x_xy_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 16);

        auto to_0_x_xy_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 17);

        auto to_0_x_xy_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 18);

        auto to_0_x_xy_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_0_x_xy_xxx, to_0_x_xy_xxy, to_0_x_xy_xxz, to_0_x_xy_xyy, to_0_x_xy_xyz, to_0_x_xy_xzz, to_0_x_xy_yyy, to_0_x_xy_yyz, to_0_x_xy_yzz, to_0_x_xy_zzz, to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xy_xxx[k] = -3.0 * to_xy_xx[k] + 2.0 * to_xy_xxxx[k] * tke_0;

            to_0_x_xy_xxy[k] = -2.0 * to_xy_xy[k] + 2.0 * to_xy_xxxy[k] * tke_0;

            to_0_x_xy_xxz[k] = -2.0 * to_xy_xz[k] + 2.0 * to_xy_xxxz[k] * tke_0;

            to_0_x_xy_xyy[k] = -to_xy_yy[k] + 2.0 * to_xy_xxyy[k] * tke_0;

            to_0_x_xy_xyz[k] = -to_xy_yz[k] + 2.0 * to_xy_xxyz[k] * tke_0;

            to_0_x_xy_xzz[k] = -to_xy_zz[k] + 2.0 * to_xy_xxzz[k] * tke_0;

            to_0_x_xy_yyy[k] = 2.0 * to_xy_xyyy[k] * tke_0;

            to_0_x_xy_yyz[k] = 2.0 * to_xy_xyyz[k] * tke_0;

            to_0_x_xy_yzz[k] = 2.0 * to_xy_xyzz[k] * tke_0;

            to_0_x_xy_zzz[k] = 2.0 * to_xy_xzzz[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : DF

        auto to_0_x_xz_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 20);

        auto to_0_x_xz_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 21);

        auto to_0_x_xz_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 22);

        auto to_0_x_xz_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 23);

        auto to_0_x_xz_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 24);

        auto to_0_x_xz_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 25);

        auto to_0_x_xz_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 26);

        auto to_0_x_xz_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 27);

        auto to_0_x_xz_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 28);

        auto to_0_x_xz_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_x_xz_xxx, to_0_x_xz_xxy, to_0_x_xz_xxz, to_0_x_xz_xyy, to_0_x_xz_xyz, to_0_x_xz_xzz, to_0_x_xz_yyy, to_0_x_xz_yyz, to_0_x_xz_yzz, to_0_x_xz_zzz, to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xz_xxx[k] = -3.0 * to_xz_xx[k] + 2.0 * to_xz_xxxx[k] * tke_0;

            to_0_x_xz_xxy[k] = -2.0 * to_xz_xy[k] + 2.0 * to_xz_xxxy[k] * tke_0;

            to_0_x_xz_xxz[k] = -2.0 * to_xz_xz[k] + 2.0 * to_xz_xxxz[k] * tke_0;

            to_0_x_xz_xyy[k] = -to_xz_yy[k] + 2.0 * to_xz_xxyy[k] * tke_0;

            to_0_x_xz_xyz[k] = -to_xz_yz[k] + 2.0 * to_xz_xxyz[k] * tke_0;

            to_0_x_xz_xzz[k] = -to_xz_zz[k] + 2.0 * to_xz_xxzz[k] * tke_0;

            to_0_x_xz_yyy[k] = 2.0 * to_xz_xyyy[k] * tke_0;

            to_0_x_xz_yyz[k] = 2.0 * to_xz_xyyz[k] * tke_0;

            to_0_x_xz_yzz[k] = 2.0 * to_xz_xyzz[k] * tke_0;

            to_0_x_xz_zzz[k] = 2.0 * to_xz_xzzz[k] * tke_0;
        }

        // Set up 30-40 components of targeted buffer : DF

        auto to_0_x_yy_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 30);

        auto to_0_x_yy_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 31);

        auto to_0_x_yy_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 32);

        auto to_0_x_yy_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 33);

        auto to_0_x_yy_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 34);

        auto to_0_x_yy_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 35);

        auto to_0_x_yy_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 36);

        auto to_0_x_yy_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 37);

        auto to_0_x_yy_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 38);

        auto to_0_x_yy_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_0_x_yy_xxx, to_0_x_yy_xxy, to_0_x_yy_xxz, to_0_x_yy_xyy, to_0_x_yy_xyz, to_0_x_yy_xzz, to_0_x_yy_yyy, to_0_x_yy_yyz, to_0_x_yy_yzz, to_0_x_yy_zzz, to_yy_xx, to_yy_xxxx, to_yy_xxxy, to_yy_xxxz, to_yy_xxyy, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yy_xxx[k] = -3.0 * to_yy_xx[k] + 2.0 * to_yy_xxxx[k] * tke_0;

            to_0_x_yy_xxy[k] = -2.0 * to_yy_xy[k] + 2.0 * to_yy_xxxy[k] * tke_0;

            to_0_x_yy_xxz[k] = -2.0 * to_yy_xz[k] + 2.0 * to_yy_xxxz[k] * tke_0;

            to_0_x_yy_xyy[k] = -to_yy_yy[k] + 2.0 * to_yy_xxyy[k] * tke_0;

            to_0_x_yy_xyz[k] = -to_yy_yz[k] + 2.0 * to_yy_xxyz[k] * tke_0;

            to_0_x_yy_xzz[k] = -to_yy_zz[k] + 2.0 * to_yy_xxzz[k] * tke_0;

            to_0_x_yy_yyy[k] = 2.0 * to_yy_xyyy[k] * tke_0;

            to_0_x_yy_yyz[k] = 2.0 * to_yy_xyyz[k] * tke_0;

            to_0_x_yy_yzz[k] = 2.0 * to_yy_xyzz[k] * tke_0;

            to_0_x_yy_zzz[k] = 2.0 * to_yy_xzzz[k] * tke_0;
        }

        // Set up 40-50 components of targeted buffer : DF

        auto to_0_x_yz_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 40);

        auto to_0_x_yz_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 41);

        auto to_0_x_yz_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 42);

        auto to_0_x_yz_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 43);

        auto to_0_x_yz_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 44);

        auto to_0_x_yz_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 45);

        auto to_0_x_yz_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 46);

        auto to_0_x_yz_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 47);

        auto to_0_x_yz_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 48);

        auto to_0_x_yz_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_0_x_yz_xxx, to_0_x_yz_xxy, to_0_x_yz_xxz, to_0_x_yz_xyy, to_0_x_yz_xyz, to_0_x_yz_xzz, to_0_x_yz_yyy, to_0_x_yz_yyz, to_0_x_yz_yzz, to_0_x_yz_zzz, to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yz_xxx[k] = -3.0 * to_yz_xx[k] + 2.0 * to_yz_xxxx[k] * tke_0;

            to_0_x_yz_xxy[k] = -2.0 * to_yz_xy[k] + 2.0 * to_yz_xxxy[k] * tke_0;

            to_0_x_yz_xxz[k] = -2.0 * to_yz_xz[k] + 2.0 * to_yz_xxxz[k] * tke_0;

            to_0_x_yz_xyy[k] = -to_yz_yy[k] + 2.0 * to_yz_xxyy[k] * tke_0;

            to_0_x_yz_xyz[k] = -to_yz_yz[k] + 2.0 * to_yz_xxyz[k] * tke_0;

            to_0_x_yz_xzz[k] = -to_yz_zz[k] + 2.0 * to_yz_xxzz[k] * tke_0;

            to_0_x_yz_yyy[k] = 2.0 * to_yz_xyyy[k] * tke_0;

            to_0_x_yz_yyz[k] = 2.0 * to_yz_xyyz[k] * tke_0;

            to_0_x_yz_yzz[k] = 2.0 * to_yz_xyzz[k] * tke_0;

            to_0_x_yz_zzz[k] = 2.0 * to_yz_xzzz[k] * tke_0;
        }

        // Set up 50-60 components of targeted buffer : DF

        auto to_0_x_zz_xxx = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 50);

        auto to_0_x_zz_xxy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 51);

        auto to_0_x_zz_xxz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 52);

        auto to_0_x_zz_xyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 53);

        auto to_0_x_zz_xyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 54);

        auto to_0_x_zz_xzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 55);

        auto to_0_x_zz_yyy = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 56);

        auto to_0_x_zz_yyz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 57);

        auto to_0_x_zz_yzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 58);

        auto to_0_x_zz_zzz = pbuffer.data(idx_op_geom_001_df + 0 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_x_zz_xxx, to_0_x_zz_xxy, to_0_x_zz_xxz, to_0_x_zz_xyy, to_0_x_zz_xyz, to_0_x_zz_xzz, to_0_x_zz_yyy, to_0_x_zz_yyz, to_0_x_zz_yzz, to_0_x_zz_zzz, to_zz_xx, to_zz_xxxx, to_zz_xxxy, to_zz_xxxz, to_zz_xxyy, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zz_xxx[k] = -3.0 * to_zz_xx[k] + 2.0 * to_zz_xxxx[k] * tke_0;

            to_0_x_zz_xxy[k] = -2.0 * to_zz_xy[k] + 2.0 * to_zz_xxxy[k] * tke_0;

            to_0_x_zz_xxz[k] = -2.0 * to_zz_xz[k] + 2.0 * to_zz_xxxz[k] * tke_0;

            to_0_x_zz_xyy[k] = -to_zz_yy[k] + 2.0 * to_zz_xxyy[k] * tke_0;

            to_0_x_zz_xyz[k] = -to_zz_yz[k] + 2.0 * to_zz_xxyz[k] * tke_0;

            to_0_x_zz_xzz[k] = -to_zz_zz[k] + 2.0 * to_zz_xxzz[k] * tke_0;

            to_0_x_zz_yyy[k] = 2.0 * to_zz_xyyy[k] * tke_0;

            to_0_x_zz_yyz[k] = 2.0 * to_zz_xyyz[k] * tke_0;

            to_0_x_zz_yzz[k] = 2.0 * to_zz_xyzz[k] * tke_0;

            to_0_x_zz_zzz[k] = 2.0 * to_zz_xzzz[k] * tke_0;
        }

        // Set up 60-70 components of targeted buffer : DF

        auto to_0_y_xx_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 0);

        auto to_0_y_xx_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 1);

        auto to_0_y_xx_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 2);

        auto to_0_y_xx_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 3);

        auto to_0_y_xx_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 4);

        auto to_0_y_xx_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 5);

        auto to_0_y_xx_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 6);

        auto to_0_y_xx_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 7);

        auto to_0_y_xx_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 8);

        auto to_0_y_xx_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_0_y_xx_xxx, to_0_y_xx_xxy, to_0_y_xx_xxz, to_0_y_xx_xyy, to_0_y_xx_xyz, to_0_y_xx_xzz, to_0_y_xx_yyy, to_0_y_xx_yyz, to_0_y_xx_yzz, to_0_y_xx_zzz, to_xx_xx, to_xx_xxxy, to_xx_xxyy, to_xx_xxyz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_yy, to_xx_yyyy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xx_xxx[k] = 2.0 * to_xx_xxxy[k] * tke_0;

            to_0_y_xx_xxy[k] = -to_xx_xx[k] + 2.0 * to_xx_xxyy[k] * tke_0;

            to_0_y_xx_xxz[k] = 2.0 * to_xx_xxyz[k] * tke_0;

            to_0_y_xx_xyy[k] = -2.0 * to_xx_xy[k] + 2.0 * to_xx_xyyy[k] * tke_0;

            to_0_y_xx_xyz[k] = -to_xx_xz[k] + 2.0 * to_xx_xyyz[k] * tke_0;

            to_0_y_xx_xzz[k] = 2.0 * to_xx_xyzz[k] * tke_0;

            to_0_y_xx_yyy[k] = -3.0 * to_xx_yy[k] + 2.0 * to_xx_yyyy[k] * tke_0;

            to_0_y_xx_yyz[k] = -2.0 * to_xx_yz[k] + 2.0 * to_xx_yyyz[k] * tke_0;

            to_0_y_xx_yzz[k] = -to_xx_zz[k] + 2.0 * to_xx_yyzz[k] * tke_0;

            to_0_y_xx_zzz[k] = 2.0 * to_xx_yzzz[k] * tke_0;
        }

        // Set up 70-80 components of targeted buffer : DF

        auto to_0_y_xy_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 10);

        auto to_0_y_xy_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 11);

        auto to_0_y_xy_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 12);

        auto to_0_y_xy_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 13);

        auto to_0_y_xy_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 14);

        auto to_0_y_xy_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 15);

        auto to_0_y_xy_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 16);

        auto to_0_y_xy_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 17);

        auto to_0_y_xy_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 18);

        auto to_0_y_xy_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_0_y_xy_xxx, to_0_y_xy_xxy, to_0_y_xy_xxz, to_0_y_xy_xyy, to_0_y_xy_xyz, to_0_y_xy_xzz, to_0_y_xy_yyy, to_0_y_xy_yyz, to_0_y_xy_yzz, to_0_y_xy_zzz, to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xy_xxx[k] = 2.0 * to_xy_xxxy[k] * tke_0;

            to_0_y_xy_xxy[k] = -to_xy_xx[k] + 2.0 * to_xy_xxyy[k] * tke_0;

            to_0_y_xy_xxz[k] = 2.0 * to_xy_xxyz[k] * tke_0;

            to_0_y_xy_xyy[k] = -2.0 * to_xy_xy[k] + 2.0 * to_xy_xyyy[k] * tke_0;

            to_0_y_xy_xyz[k] = -to_xy_xz[k] + 2.0 * to_xy_xyyz[k] * tke_0;

            to_0_y_xy_xzz[k] = 2.0 * to_xy_xyzz[k] * tke_0;

            to_0_y_xy_yyy[k] = -3.0 * to_xy_yy[k] + 2.0 * to_xy_yyyy[k] * tke_0;

            to_0_y_xy_yyz[k] = -2.0 * to_xy_yz[k] + 2.0 * to_xy_yyyz[k] * tke_0;

            to_0_y_xy_yzz[k] = -to_xy_zz[k] + 2.0 * to_xy_yyzz[k] * tke_0;

            to_0_y_xy_zzz[k] = 2.0 * to_xy_yzzz[k] * tke_0;
        }

        // Set up 80-90 components of targeted buffer : DF

        auto to_0_y_xz_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 20);

        auto to_0_y_xz_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 21);

        auto to_0_y_xz_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 22);

        auto to_0_y_xz_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 23);

        auto to_0_y_xz_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 24);

        auto to_0_y_xz_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 25);

        auto to_0_y_xz_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 26);

        auto to_0_y_xz_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 27);

        auto to_0_y_xz_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 28);

        auto to_0_y_xz_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_y_xz_xxx, to_0_y_xz_xxy, to_0_y_xz_xxz, to_0_y_xz_xyy, to_0_y_xz_xyz, to_0_y_xz_xzz, to_0_y_xz_yyy, to_0_y_xz_yyz, to_0_y_xz_yzz, to_0_y_xz_zzz, to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xz_xxx[k] = 2.0 * to_xz_xxxy[k] * tke_0;

            to_0_y_xz_xxy[k] = -to_xz_xx[k] + 2.0 * to_xz_xxyy[k] * tke_0;

            to_0_y_xz_xxz[k] = 2.0 * to_xz_xxyz[k] * tke_0;

            to_0_y_xz_xyy[k] = -2.0 * to_xz_xy[k] + 2.0 * to_xz_xyyy[k] * tke_0;

            to_0_y_xz_xyz[k] = -to_xz_xz[k] + 2.0 * to_xz_xyyz[k] * tke_0;

            to_0_y_xz_xzz[k] = 2.0 * to_xz_xyzz[k] * tke_0;

            to_0_y_xz_yyy[k] = -3.0 * to_xz_yy[k] + 2.0 * to_xz_yyyy[k] * tke_0;

            to_0_y_xz_yyz[k] = -2.0 * to_xz_yz[k] + 2.0 * to_xz_yyyz[k] * tke_0;

            to_0_y_xz_yzz[k] = -to_xz_zz[k] + 2.0 * to_xz_yyzz[k] * tke_0;

            to_0_y_xz_zzz[k] = 2.0 * to_xz_yzzz[k] * tke_0;
        }

        // Set up 90-100 components of targeted buffer : DF

        auto to_0_y_yy_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 30);

        auto to_0_y_yy_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 31);

        auto to_0_y_yy_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 32);

        auto to_0_y_yy_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 33);

        auto to_0_y_yy_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 34);

        auto to_0_y_yy_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 35);

        auto to_0_y_yy_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 36);

        auto to_0_y_yy_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 37);

        auto to_0_y_yy_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 38);

        auto to_0_y_yy_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_0_y_yy_xxx, to_0_y_yy_xxy, to_0_y_yy_xxz, to_0_y_yy_xyy, to_0_y_yy_xyz, to_0_y_yy_xzz, to_0_y_yy_yyy, to_0_y_yy_yyz, to_0_y_yy_yzz, to_0_y_yy_zzz, to_yy_xx, to_yy_xxxy, to_yy_xxyy, to_yy_xxyz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_yy, to_yy_yyyy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yy_xxx[k] = 2.0 * to_yy_xxxy[k] * tke_0;

            to_0_y_yy_xxy[k] = -to_yy_xx[k] + 2.0 * to_yy_xxyy[k] * tke_0;

            to_0_y_yy_xxz[k] = 2.0 * to_yy_xxyz[k] * tke_0;

            to_0_y_yy_xyy[k] = -2.0 * to_yy_xy[k] + 2.0 * to_yy_xyyy[k] * tke_0;

            to_0_y_yy_xyz[k] = -to_yy_xz[k] + 2.0 * to_yy_xyyz[k] * tke_0;

            to_0_y_yy_xzz[k] = 2.0 * to_yy_xyzz[k] * tke_0;

            to_0_y_yy_yyy[k] = -3.0 * to_yy_yy[k] + 2.0 * to_yy_yyyy[k] * tke_0;

            to_0_y_yy_yyz[k] = -2.0 * to_yy_yz[k] + 2.0 * to_yy_yyyz[k] * tke_0;

            to_0_y_yy_yzz[k] = -to_yy_zz[k] + 2.0 * to_yy_yyzz[k] * tke_0;

            to_0_y_yy_zzz[k] = 2.0 * to_yy_yzzz[k] * tke_0;
        }

        // Set up 100-110 components of targeted buffer : DF

        auto to_0_y_yz_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 40);

        auto to_0_y_yz_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 41);

        auto to_0_y_yz_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 42);

        auto to_0_y_yz_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 43);

        auto to_0_y_yz_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 44);

        auto to_0_y_yz_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 45);

        auto to_0_y_yz_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 46);

        auto to_0_y_yz_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 47);

        auto to_0_y_yz_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 48);

        auto to_0_y_yz_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_0_y_yz_xxx, to_0_y_yz_xxy, to_0_y_yz_xxz, to_0_y_yz_xyy, to_0_y_yz_xyz, to_0_y_yz_xzz, to_0_y_yz_yyy, to_0_y_yz_yyz, to_0_y_yz_yzz, to_0_y_yz_zzz, to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yz_xxx[k] = 2.0 * to_yz_xxxy[k] * tke_0;

            to_0_y_yz_xxy[k] = -to_yz_xx[k] + 2.0 * to_yz_xxyy[k] * tke_0;

            to_0_y_yz_xxz[k] = 2.0 * to_yz_xxyz[k] * tke_0;

            to_0_y_yz_xyy[k] = -2.0 * to_yz_xy[k] + 2.0 * to_yz_xyyy[k] * tke_0;

            to_0_y_yz_xyz[k] = -to_yz_xz[k] + 2.0 * to_yz_xyyz[k] * tke_0;

            to_0_y_yz_xzz[k] = 2.0 * to_yz_xyzz[k] * tke_0;

            to_0_y_yz_yyy[k] = -3.0 * to_yz_yy[k] + 2.0 * to_yz_yyyy[k] * tke_0;

            to_0_y_yz_yyz[k] = -2.0 * to_yz_yz[k] + 2.0 * to_yz_yyyz[k] * tke_0;

            to_0_y_yz_yzz[k] = -to_yz_zz[k] + 2.0 * to_yz_yyzz[k] * tke_0;

            to_0_y_yz_zzz[k] = 2.0 * to_yz_yzzz[k] * tke_0;
        }

        // Set up 110-120 components of targeted buffer : DF

        auto to_0_y_zz_xxx = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 50);

        auto to_0_y_zz_xxy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 51);

        auto to_0_y_zz_xxz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 52);

        auto to_0_y_zz_xyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 53);

        auto to_0_y_zz_xyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 54);

        auto to_0_y_zz_xzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 55);

        auto to_0_y_zz_yyy = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 56);

        auto to_0_y_zz_yyz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 57);

        auto to_0_y_zz_yzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 58);

        auto to_0_y_zz_zzz = pbuffer.data(idx_op_geom_001_df + 1 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_y_zz_xxx, to_0_y_zz_xxy, to_0_y_zz_xxz, to_0_y_zz_xyy, to_0_y_zz_xyz, to_0_y_zz_xzz, to_0_y_zz_yyy, to_0_y_zz_yyz, to_0_y_zz_yzz, to_0_y_zz_zzz, to_zz_xx, to_zz_xxxy, to_zz_xxyy, to_zz_xxyz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_yy, to_zz_yyyy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zz_xxx[k] = 2.0 * to_zz_xxxy[k] * tke_0;

            to_0_y_zz_xxy[k] = -to_zz_xx[k] + 2.0 * to_zz_xxyy[k] * tke_0;

            to_0_y_zz_xxz[k] = 2.0 * to_zz_xxyz[k] * tke_0;

            to_0_y_zz_xyy[k] = -2.0 * to_zz_xy[k] + 2.0 * to_zz_xyyy[k] * tke_0;

            to_0_y_zz_xyz[k] = -to_zz_xz[k] + 2.0 * to_zz_xyyz[k] * tke_0;

            to_0_y_zz_xzz[k] = 2.0 * to_zz_xyzz[k] * tke_0;

            to_0_y_zz_yyy[k] = -3.0 * to_zz_yy[k] + 2.0 * to_zz_yyyy[k] * tke_0;

            to_0_y_zz_yyz[k] = -2.0 * to_zz_yz[k] + 2.0 * to_zz_yyyz[k] * tke_0;

            to_0_y_zz_yzz[k] = -to_zz_zz[k] + 2.0 * to_zz_yyzz[k] * tke_0;

            to_0_y_zz_zzz[k] = 2.0 * to_zz_yzzz[k] * tke_0;
        }

        // Set up 120-130 components of targeted buffer : DF

        auto to_0_z_xx_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 0);

        auto to_0_z_xx_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 1);

        auto to_0_z_xx_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 2);

        auto to_0_z_xx_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 3);

        auto to_0_z_xx_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 4);

        auto to_0_z_xx_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 5);

        auto to_0_z_xx_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 6);

        auto to_0_z_xx_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 7);

        auto to_0_z_xx_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 8);

        auto to_0_z_xx_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_0_z_xx_xxx, to_0_z_xx_xxy, to_0_z_xx_xxz, to_0_z_xx_xyy, to_0_z_xx_xyz, to_0_z_xx_xzz, to_0_z_xx_yyy, to_0_z_xx_yyz, to_0_z_xx_yzz, to_0_z_xx_zzz, to_xx_xx, to_xx_xxxz, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xx_xxx[k] = 2.0 * to_xx_xxxz[k] * tke_0;

            to_0_z_xx_xxy[k] = 2.0 * to_xx_xxyz[k] * tke_0;

            to_0_z_xx_xxz[k] = -to_xx_xx[k] + 2.0 * to_xx_xxzz[k] * tke_0;

            to_0_z_xx_xyy[k] = 2.0 * to_xx_xyyz[k] * tke_0;

            to_0_z_xx_xyz[k] = -to_xx_xy[k] + 2.0 * to_xx_xyzz[k] * tke_0;

            to_0_z_xx_xzz[k] = -2.0 * to_xx_xz[k] + 2.0 * to_xx_xzzz[k] * tke_0;

            to_0_z_xx_yyy[k] = 2.0 * to_xx_yyyz[k] * tke_0;

            to_0_z_xx_yyz[k] = -to_xx_yy[k] + 2.0 * to_xx_yyzz[k] * tke_0;

            to_0_z_xx_yzz[k] = -2.0 * to_xx_yz[k] + 2.0 * to_xx_yzzz[k] * tke_0;

            to_0_z_xx_zzz[k] = -3.0 * to_xx_zz[k] + 2.0 * to_xx_zzzz[k] * tke_0;
        }

        // Set up 130-140 components of targeted buffer : DF

        auto to_0_z_xy_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 10);

        auto to_0_z_xy_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 11);

        auto to_0_z_xy_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 12);

        auto to_0_z_xy_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 13);

        auto to_0_z_xy_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 14);

        auto to_0_z_xy_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 15);

        auto to_0_z_xy_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 16);

        auto to_0_z_xy_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 17);

        auto to_0_z_xy_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 18);

        auto to_0_z_xy_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_0_z_xy_xxx, to_0_z_xy_xxy, to_0_z_xy_xxz, to_0_z_xy_xyy, to_0_z_xy_xyz, to_0_z_xy_xzz, to_0_z_xy_yyy, to_0_z_xy_yyz, to_0_z_xy_yzz, to_0_z_xy_zzz, to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xy_xxx[k] = 2.0 * to_xy_xxxz[k] * tke_0;

            to_0_z_xy_xxy[k] = 2.0 * to_xy_xxyz[k] * tke_0;

            to_0_z_xy_xxz[k] = -to_xy_xx[k] + 2.0 * to_xy_xxzz[k] * tke_0;

            to_0_z_xy_xyy[k] = 2.0 * to_xy_xyyz[k] * tke_0;

            to_0_z_xy_xyz[k] = -to_xy_xy[k] + 2.0 * to_xy_xyzz[k] * tke_0;

            to_0_z_xy_xzz[k] = -2.0 * to_xy_xz[k] + 2.0 * to_xy_xzzz[k] * tke_0;

            to_0_z_xy_yyy[k] = 2.0 * to_xy_yyyz[k] * tke_0;

            to_0_z_xy_yyz[k] = -to_xy_yy[k] + 2.0 * to_xy_yyzz[k] * tke_0;

            to_0_z_xy_yzz[k] = -2.0 * to_xy_yz[k] + 2.0 * to_xy_yzzz[k] * tke_0;

            to_0_z_xy_zzz[k] = -3.0 * to_xy_zz[k] + 2.0 * to_xy_zzzz[k] * tke_0;
        }

        // Set up 140-150 components of targeted buffer : DF

        auto to_0_z_xz_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 20);

        auto to_0_z_xz_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 21);

        auto to_0_z_xz_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 22);

        auto to_0_z_xz_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 23);

        auto to_0_z_xz_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 24);

        auto to_0_z_xz_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 25);

        auto to_0_z_xz_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 26);

        auto to_0_z_xz_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 27);

        auto to_0_z_xz_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 28);

        auto to_0_z_xz_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_z_xz_xxx, to_0_z_xz_xxy, to_0_z_xz_xxz, to_0_z_xz_xyy, to_0_z_xz_xyz, to_0_z_xz_xzz, to_0_z_xz_yyy, to_0_z_xz_yyz, to_0_z_xz_yzz, to_0_z_xz_zzz, to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xz_xxx[k] = 2.0 * to_xz_xxxz[k] * tke_0;

            to_0_z_xz_xxy[k] = 2.0 * to_xz_xxyz[k] * tke_0;

            to_0_z_xz_xxz[k] = -to_xz_xx[k] + 2.0 * to_xz_xxzz[k] * tke_0;

            to_0_z_xz_xyy[k] = 2.0 * to_xz_xyyz[k] * tke_0;

            to_0_z_xz_xyz[k] = -to_xz_xy[k] + 2.0 * to_xz_xyzz[k] * tke_0;

            to_0_z_xz_xzz[k] = -2.0 * to_xz_xz[k] + 2.0 * to_xz_xzzz[k] * tke_0;

            to_0_z_xz_yyy[k] = 2.0 * to_xz_yyyz[k] * tke_0;

            to_0_z_xz_yyz[k] = -to_xz_yy[k] + 2.0 * to_xz_yyzz[k] * tke_0;

            to_0_z_xz_yzz[k] = -2.0 * to_xz_yz[k] + 2.0 * to_xz_yzzz[k] * tke_0;

            to_0_z_xz_zzz[k] = -3.0 * to_xz_zz[k] + 2.0 * to_xz_zzzz[k] * tke_0;
        }

        // Set up 150-160 components of targeted buffer : DF

        auto to_0_z_yy_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 30);

        auto to_0_z_yy_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 31);

        auto to_0_z_yy_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 32);

        auto to_0_z_yy_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 33);

        auto to_0_z_yy_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 34);

        auto to_0_z_yy_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 35);

        auto to_0_z_yy_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 36);

        auto to_0_z_yy_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 37);

        auto to_0_z_yy_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 38);

        auto to_0_z_yy_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_0_z_yy_xxx, to_0_z_yy_xxy, to_0_z_yy_xxz, to_0_z_yy_xyy, to_0_z_yy_xyz, to_0_z_yy_xzz, to_0_z_yy_yyy, to_0_z_yy_yyz, to_0_z_yy_yzz, to_0_z_yy_zzz, to_yy_xx, to_yy_xxxz, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yy_xxx[k] = 2.0 * to_yy_xxxz[k] * tke_0;

            to_0_z_yy_xxy[k] = 2.0 * to_yy_xxyz[k] * tke_0;

            to_0_z_yy_xxz[k] = -to_yy_xx[k] + 2.0 * to_yy_xxzz[k] * tke_0;

            to_0_z_yy_xyy[k] = 2.0 * to_yy_xyyz[k] * tke_0;

            to_0_z_yy_xyz[k] = -to_yy_xy[k] + 2.0 * to_yy_xyzz[k] * tke_0;

            to_0_z_yy_xzz[k] = -2.0 * to_yy_xz[k] + 2.0 * to_yy_xzzz[k] * tke_0;

            to_0_z_yy_yyy[k] = 2.0 * to_yy_yyyz[k] * tke_0;

            to_0_z_yy_yyz[k] = -to_yy_yy[k] + 2.0 * to_yy_yyzz[k] * tke_0;

            to_0_z_yy_yzz[k] = -2.0 * to_yy_yz[k] + 2.0 * to_yy_yzzz[k] * tke_0;

            to_0_z_yy_zzz[k] = -3.0 * to_yy_zz[k] + 2.0 * to_yy_zzzz[k] * tke_0;
        }

        // Set up 160-170 components of targeted buffer : DF

        auto to_0_z_yz_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 40);

        auto to_0_z_yz_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 41);

        auto to_0_z_yz_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 42);

        auto to_0_z_yz_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 43);

        auto to_0_z_yz_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 44);

        auto to_0_z_yz_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 45);

        auto to_0_z_yz_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 46);

        auto to_0_z_yz_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 47);

        auto to_0_z_yz_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 48);

        auto to_0_z_yz_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_0_z_yz_xxx, to_0_z_yz_xxy, to_0_z_yz_xxz, to_0_z_yz_xyy, to_0_z_yz_xyz, to_0_z_yz_xzz, to_0_z_yz_yyy, to_0_z_yz_yyz, to_0_z_yz_yzz, to_0_z_yz_zzz, to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yz_xxx[k] = 2.0 * to_yz_xxxz[k] * tke_0;

            to_0_z_yz_xxy[k] = 2.0 * to_yz_xxyz[k] * tke_0;

            to_0_z_yz_xxz[k] = -to_yz_xx[k] + 2.0 * to_yz_xxzz[k] * tke_0;

            to_0_z_yz_xyy[k] = 2.0 * to_yz_xyyz[k] * tke_0;

            to_0_z_yz_xyz[k] = -to_yz_xy[k] + 2.0 * to_yz_xyzz[k] * tke_0;

            to_0_z_yz_xzz[k] = -2.0 * to_yz_xz[k] + 2.0 * to_yz_xzzz[k] * tke_0;

            to_0_z_yz_yyy[k] = 2.0 * to_yz_yyyz[k] * tke_0;

            to_0_z_yz_yyz[k] = -to_yz_yy[k] + 2.0 * to_yz_yyzz[k] * tke_0;

            to_0_z_yz_yzz[k] = -2.0 * to_yz_yz[k] + 2.0 * to_yz_yzzz[k] * tke_0;

            to_0_z_yz_zzz[k] = -3.0 * to_yz_zz[k] + 2.0 * to_yz_zzzz[k] * tke_0;
        }

        // Set up 170-180 components of targeted buffer : DF

        auto to_0_z_zz_xxx = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 50);

        auto to_0_z_zz_xxy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 51);

        auto to_0_z_zz_xxz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 52);

        auto to_0_z_zz_xyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 53);

        auto to_0_z_zz_xyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 54);

        auto to_0_z_zz_xzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 55);

        auto to_0_z_zz_yyy = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 56);

        auto to_0_z_zz_yyz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 57);

        auto to_0_z_zz_yzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 58);

        auto to_0_z_zz_zzz = pbuffer.data(idx_op_geom_001_df + 2 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_z_zz_xxx, to_0_z_zz_xxy, to_0_z_zz_xxz, to_0_z_zz_xyy, to_0_z_zz_xyz, to_0_z_zz_xzz, to_0_z_zz_yyy, to_0_z_zz_yyz, to_0_z_zz_yzz, to_0_z_zz_zzz, to_zz_xx, to_zz_xxxz, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zz_xxx[k] = 2.0 * to_zz_xxxz[k] * tke_0;

            to_0_z_zz_xxy[k] = 2.0 * to_zz_xxyz[k] * tke_0;

            to_0_z_zz_xxz[k] = -to_zz_xx[k] + 2.0 * to_zz_xxzz[k] * tke_0;

            to_0_z_zz_xyy[k] = 2.0 * to_zz_xyyz[k] * tke_0;

            to_0_z_zz_xyz[k] = -to_zz_xy[k] + 2.0 * to_zz_xyzz[k] * tke_0;

            to_0_z_zz_xzz[k] = -2.0 * to_zz_xz[k] + 2.0 * to_zz_xzzz[k] * tke_0;

            to_0_z_zz_yyy[k] = 2.0 * to_zz_yyyz[k] * tke_0;

            to_0_z_zz_yyz[k] = -to_zz_yy[k] + 2.0 * to_zz_yyzz[k] * tke_0;

            to_0_z_zz_yzz[k] = -2.0 * to_zz_yz[k] + 2.0 * to_zz_yzzz[k] * tke_0;

            to_0_z_zz_zzz[k] = -3.0 * to_zz_zz[k] + 2.0 * to_zz_zzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

