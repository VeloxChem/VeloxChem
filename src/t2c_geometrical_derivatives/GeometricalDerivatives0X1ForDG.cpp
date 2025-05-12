#include "GeometricalDerivatives0X1ForDG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_dg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_dg,
                       const int idx_op_df,
                       const int idx_op_dh,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DF

        auto to_xx_xxx = pbuffer.data(idx_op_df + i * 60 + 0);

        auto to_xx_xxy = pbuffer.data(idx_op_df + i * 60 + 1);

        auto to_xx_xxz = pbuffer.data(idx_op_df + i * 60 + 2);

        auto to_xx_xyy = pbuffer.data(idx_op_df + i * 60 + 3);

        auto to_xx_xyz = pbuffer.data(idx_op_df + i * 60 + 4);

        auto to_xx_xzz = pbuffer.data(idx_op_df + i * 60 + 5);

        auto to_xx_yyy = pbuffer.data(idx_op_df + i * 60 + 6);

        auto to_xx_yyz = pbuffer.data(idx_op_df + i * 60 + 7);

        auto to_xx_yzz = pbuffer.data(idx_op_df + i * 60 + 8);

        auto to_xx_zzz = pbuffer.data(idx_op_df + i * 60 + 9);

        auto to_xy_xxx = pbuffer.data(idx_op_df + i * 60 + 10);

        auto to_xy_xxy = pbuffer.data(idx_op_df + i * 60 + 11);

        auto to_xy_xxz = pbuffer.data(idx_op_df + i * 60 + 12);

        auto to_xy_xyy = pbuffer.data(idx_op_df + i * 60 + 13);

        auto to_xy_xyz = pbuffer.data(idx_op_df + i * 60 + 14);

        auto to_xy_xzz = pbuffer.data(idx_op_df + i * 60 + 15);

        auto to_xy_yyy = pbuffer.data(idx_op_df + i * 60 + 16);

        auto to_xy_yyz = pbuffer.data(idx_op_df + i * 60 + 17);

        auto to_xy_yzz = pbuffer.data(idx_op_df + i * 60 + 18);

        auto to_xy_zzz = pbuffer.data(idx_op_df + i * 60 + 19);

        auto to_xz_xxx = pbuffer.data(idx_op_df + i * 60 + 20);

        auto to_xz_xxy = pbuffer.data(idx_op_df + i * 60 + 21);

        auto to_xz_xxz = pbuffer.data(idx_op_df + i * 60 + 22);

        auto to_xz_xyy = pbuffer.data(idx_op_df + i * 60 + 23);

        auto to_xz_xyz = pbuffer.data(idx_op_df + i * 60 + 24);

        auto to_xz_xzz = pbuffer.data(idx_op_df + i * 60 + 25);

        auto to_xz_yyy = pbuffer.data(idx_op_df + i * 60 + 26);

        auto to_xz_yyz = pbuffer.data(idx_op_df + i * 60 + 27);

        auto to_xz_yzz = pbuffer.data(idx_op_df + i * 60 + 28);

        auto to_xz_zzz = pbuffer.data(idx_op_df + i * 60 + 29);

        auto to_yy_xxx = pbuffer.data(idx_op_df + i * 60 + 30);

        auto to_yy_xxy = pbuffer.data(idx_op_df + i * 60 + 31);

        auto to_yy_xxz = pbuffer.data(idx_op_df + i * 60 + 32);

        auto to_yy_xyy = pbuffer.data(idx_op_df + i * 60 + 33);

        auto to_yy_xyz = pbuffer.data(idx_op_df + i * 60 + 34);

        auto to_yy_xzz = pbuffer.data(idx_op_df + i * 60 + 35);

        auto to_yy_yyy = pbuffer.data(idx_op_df + i * 60 + 36);

        auto to_yy_yyz = pbuffer.data(idx_op_df + i * 60 + 37);

        auto to_yy_yzz = pbuffer.data(idx_op_df + i * 60 + 38);

        auto to_yy_zzz = pbuffer.data(idx_op_df + i * 60 + 39);

        auto to_yz_xxx = pbuffer.data(idx_op_df + i * 60 + 40);

        auto to_yz_xxy = pbuffer.data(idx_op_df + i * 60 + 41);

        auto to_yz_xxz = pbuffer.data(idx_op_df + i * 60 + 42);

        auto to_yz_xyy = pbuffer.data(idx_op_df + i * 60 + 43);

        auto to_yz_xyz = pbuffer.data(idx_op_df + i * 60 + 44);

        auto to_yz_xzz = pbuffer.data(idx_op_df + i * 60 + 45);

        auto to_yz_yyy = pbuffer.data(idx_op_df + i * 60 + 46);

        auto to_yz_yyz = pbuffer.data(idx_op_df + i * 60 + 47);

        auto to_yz_yzz = pbuffer.data(idx_op_df + i * 60 + 48);

        auto to_yz_zzz = pbuffer.data(idx_op_df + i * 60 + 49);

        auto to_zz_xxx = pbuffer.data(idx_op_df + i * 60 + 50);

        auto to_zz_xxy = pbuffer.data(idx_op_df + i * 60 + 51);

        auto to_zz_xxz = pbuffer.data(idx_op_df + i * 60 + 52);

        auto to_zz_xyy = pbuffer.data(idx_op_df + i * 60 + 53);

        auto to_zz_xyz = pbuffer.data(idx_op_df + i * 60 + 54);

        auto to_zz_xzz = pbuffer.data(idx_op_df + i * 60 + 55);

        auto to_zz_yyy = pbuffer.data(idx_op_df + i * 60 + 56);

        auto to_zz_yyz = pbuffer.data(idx_op_df + i * 60 + 57);

        auto to_zz_yzz = pbuffer.data(idx_op_df + i * 60 + 58);

        auto to_zz_zzz = pbuffer.data(idx_op_df + i * 60 + 59);

        // Set up components of auxiliary buffer : DH

        auto to_xx_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 0);

        auto to_xx_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 1);

        auto to_xx_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 2);

        auto to_xx_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 3);

        auto to_xx_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 4);

        auto to_xx_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 5);

        auto to_xx_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 6);

        auto to_xx_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 7);

        auto to_xx_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 8);

        auto to_xx_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 9);

        auto to_xx_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 10);

        auto to_xx_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 11);

        auto to_xx_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 12);

        auto to_xx_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 13);

        auto to_xx_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 14);

        auto to_xx_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 15);

        auto to_xx_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 16);

        auto to_xx_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 17);

        auto to_xx_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 18);

        auto to_xx_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 19);

        auto to_xx_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 20);

        auto to_xy_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 21);

        auto to_xy_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 22);

        auto to_xy_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 23);

        auto to_xy_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 24);

        auto to_xy_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 25);

        auto to_xy_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 26);

        auto to_xy_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 27);

        auto to_xy_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 28);

        auto to_xy_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 29);

        auto to_xy_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 30);

        auto to_xy_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 31);

        auto to_xy_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 32);

        auto to_xy_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 33);

        auto to_xy_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 34);

        auto to_xy_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 35);

        auto to_xy_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 36);

        auto to_xy_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 37);

        auto to_xy_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 38);

        auto to_xy_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 39);

        auto to_xy_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 40);

        auto to_xy_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 41);

        auto to_xz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 42);

        auto to_xz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 43);

        auto to_xz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 44);

        auto to_xz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 45);

        auto to_xz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 46);

        auto to_xz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 47);

        auto to_xz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 48);

        auto to_xz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 49);

        auto to_xz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 50);

        auto to_xz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 51);

        auto to_xz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 52);

        auto to_xz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 53);

        auto to_xz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 54);

        auto to_xz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 55);

        auto to_xz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 56);

        auto to_xz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 57);

        auto to_xz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 58);

        auto to_xz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 59);

        auto to_xz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 60);

        auto to_xz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 61);

        auto to_xz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 62);

        auto to_yy_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 63);

        auto to_yy_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 64);

        auto to_yy_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 65);

        auto to_yy_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 66);

        auto to_yy_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 67);

        auto to_yy_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 68);

        auto to_yy_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 69);

        auto to_yy_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 70);

        auto to_yy_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 71);

        auto to_yy_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 72);

        auto to_yy_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 73);

        auto to_yy_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 74);

        auto to_yy_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 75);

        auto to_yy_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 76);

        auto to_yy_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 77);

        auto to_yy_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 78);

        auto to_yy_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 79);

        auto to_yy_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 80);

        auto to_yy_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 81);

        auto to_yy_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 82);

        auto to_yy_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 83);

        auto to_yz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 84);

        auto to_yz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 85);

        auto to_yz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 86);

        auto to_yz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 87);

        auto to_yz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 88);

        auto to_yz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 89);

        auto to_yz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 90);

        auto to_yz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 91);

        auto to_yz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 92);

        auto to_yz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 93);

        auto to_yz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 94);

        auto to_yz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 95);

        auto to_yz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 96);

        auto to_yz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 97);

        auto to_yz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 98);

        auto to_yz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 99);

        auto to_yz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 100);

        auto to_yz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 101);

        auto to_yz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 102);

        auto to_yz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 103);

        auto to_yz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 104);

        auto to_zz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 105);

        auto to_zz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 106);

        auto to_zz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 107);

        auto to_zz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 108);

        auto to_zz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 109);

        auto to_zz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 110);

        auto to_zz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 111);

        auto to_zz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 112);

        auto to_zz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 113);

        auto to_zz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 114);

        auto to_zz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 115);

        auto to_zz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 116);

        auto to_zz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 117);

        auto to_zz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 118);

        auto to_zz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 119);

        auto to_zz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 120);

        auto to_zz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 121);

        auto to_zz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 122);

        auto to_zz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 123);

        auto to_zz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 124);

        auto to_zz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 125);

        // Set up 0-15 components of targeted buffer : DG

        auto to_0_x_xx_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 0);

        auto to_0_x_xx_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 1);

        auto to_0_x_xx_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 2);

        auto to_0_x_xx_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 3);

        auto to_0_x_xx_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 4);

        auto to_0_x_xx_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 5);

        auto to_0_x_xx_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 6);

        auto to_0_x_xx_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 7);

        auto to_0_x_xx_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 8);

        auto to_0_x_xx_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 9);

        auto to_0_x_xx_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 10);

        auto to_0_x_xx_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 11);

        auto to_0_x_xx_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 12);

        auto to_0_x_xx_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 13);

        auto to_0_x_xx_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 14);

        #pragma omp simd aligned(to_0_x_xx_xxxx, to_0_x_xx_xxxy, to_0_x_xx_xxxz, to_0_x_xx_xxyy, to_0_x_xx_xxyz, to_0_x_xx_xxzz, to_0_x_xx_xyyy, to_0_x_xx_xyyz, to_0_x_xx_xyzz, to_0_x_xx_xzzz, to_0_x_xx_yyyy, to_0_x_xx_yyyz, to_0_x_xx_yyzz, to_0_x_xx_yzzz, to_0_x_xx_zzzz, to_xx_xxx, to_xx_xxxxx, to_xx_xxxxy, to_xx_xxxxz, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xx_xxxx[k] = -4.0 * to_xx_xxx[k] + 2.0 * to_xx_xxxxx[k] * tke_0;

            to_0_x_xx_xxxy[k] = -3.0 * to_xx_xxy[k] + 2.0 * to_xx_xxxxy[k] * tke_0;

            to_0_x_xx_xxxz[k] = -3.0 * to_xx_xxz[k] + 2.0 * to_xx_xxxxz[k] * tke_0;

            to_0_x_xx_xxyy[k] = -2.0 * to_xx_xyy[k] + 2.0 * to_xx_xxxyy[k] * tke_0;

            to_0_x_xx_xxyz[k] = -2.0 * to_xx_xyz[k] + 2.0 * to_xx_xxxyz[k] * tke_0;

            to_0_x_xx_xxzz[k] = -2.0 * to_xx_xzz[k] + 2.0 * to_xx_xxxzz[k] * tke_0;

            to_0_x_xx_xyyy[k] = -to_xx_yyy[k] + 2.0 * to_xx_xxyyy[k] * tke_0;

            to_0_x_xx_xyyz[k] = -to_xx_yyz[k] + 2.0 * to_xx_xxyyz[k] * tke_0;

            to_0_x_xx_xyzz[k] = -to_xx_yzz[k] + 2.0 * to_xx_xxyzz[k] * tke_0;

            to_0_x_xx_xzzz[k] = -to_xx_zzz[k] + 2.0 * to_xx_xxzzz[k] * tke_0;

            to_0_x_xx_yyyy[k] = 2.0 * to_xx_xyyyy[k] * tke_0;

            to_0_x_xx_yyyz[k] = 2.0 * to_xx_xyyyz[k] * tke_0;

            to_0_x_xx_yyzz[k] = 2.0 * to_xx_xyyzz[k] * tke_0;

            to_0_x_xx_yzzz[k] = 2.0 * to_xx_xyzzz[k] * tke_0;

            to_0_x_xx_zzzz[k] = 2.0 * to_xx_xzzzz[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : DG

        auto to_0_x_xy_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 15);

        auto to_0_x_xy_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 16);

        auto to_0_x_xy_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 17);

        auto to_0_x_xy_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 18);

        auto to_0_x_xy_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 19);

        auto to_0_x_xy_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 20);

        auto to_0_x_xy_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 21);

        auto to_0_x_xy_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 22);

        auto to_0_x_xy_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 23);

        auto to_0_x_xy_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 24);

        auto to_0_x_xy_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 25);

        auto to_0_x_xy_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 26);

        auto to_0_x_xy_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 27);

        auto to_0_x_xy_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 28);

        auto to_0_x_xy_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_x_xy_xxxx, to_0_x_xy_xxxy, to_0_x_xy_xxxz, to_0_x_xy_xxyy, to_0_x_xy_xxyz, to_0_x_xy_xxzz, to_0_x_xy_xyyy, to_0_x_xy_xyyz, to_0_x_xy_xyzz, to_0_x_xy_xzzz, to_0_x_xy_yyyy, to_0_x_xy_yyyz, to_0_x_xy_yyzz, to_0_x_xy_yzzz, to_0_x_xy_zzzz, to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xy_xxxx[k] = -4.0 * to_xy_xxx[k] + 2.0 * to_xy_xxxxx[k] * tke_0;

            to_0_x_xy_xxxy[k] = -3.0 * to_xy_xxy[k] + 2.0 * to_xy_xxxxy[k] * tke_0;

            to_0_x_xy_xxxz[k] = -3.0 * to_xy_xxz[k] + 2.0 * to_xy_xxxxz[k] * tke_0;

            to_0_x_xy_xxyy[k] = -2.0 * to_xy_xyy[k] + 2.0 * to_xy_xxxyy[k] * tke_0;

            to_0_x_xy_xxyz[k] = -2.0 * to_xy_xyz[k] + 2.0 * to_xy_xxxyz[k] * tke_0;

            to_0_x_xy_xxzz[k] = -2.0 * to_xy_xzz[k] + 2.0 * to_xy_xxxzz[k] * tke_0;

            to_0_x_xy_xyyy[k] = -to_xy_yyy[k] + 2.0 * to_xy_xxyyy[k] * tke_0;

            to_0_x_xy_xyyz[k] = -to_xy_yyz[k] + 2.0 * to_xy_xxyyz[k] * tke_0;

            to_0_x_xy_xyzz[k] = -to_xy_yzz[k] + 2.0 * to_xy_xxyzz[k] * tke_0;

            to_0_x_xy_xzzz[k] = -to_xy_zzz[k] + 2.0 * to_xy_xxzzz[k] * tke_0;

            to_0_x_xy_yyyy[k] = 2.0 * to_xy_xyyyy[k] * tke_0;

            to_0_x_xy_yyyz[k] = 2.0 * to_xy_xyyyz[k] * tke_0;

            to_0_x_xy_yyzz[k] = 2.0 * to_xy_xyyzz[k] * tke_0;

            to_0_x_xy_yzzz[k] = 2.0 * to_xy_xyzzz[k] * tke_0;

            to_0_x_xy_zzzz[k] = 2.0 * to_xy_xzzzz[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : DG

        auto to_0_x_xz_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 30);

        auto to_0_x_xz_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 31);

        auto to_0_x_xz_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 32);

        auto to_0_x_xz_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 33);

        auto to_0_x_xz_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 34);

        auto to_0_x_xz_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 35);

        auto to_0_x_xz_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 36);

        auto to_0_x_xz_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 37);

        auto to_0_x_xz_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 38);

        auto to_0_x_xz_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 39);

        auto to_0_x_xz_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 40);

        auto to_0_x_xz_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 41);

        auto to_0_x_xz_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 42);

        auto to_0_x_xz_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 43);

        auto to_0_x_xz_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 44);

        #pragma omp simd aligned(to_0_x_xz_xxxx, to_0_x_xz_xxxy, to_0_x_xz_xxxz, to_0_x_xz_xxyy, to_0_x_xz_xxyz, to_0_x_xz_xxzz, to_0_x_xz_xyyy, to_0_x_xz_xyyz, to_0_x_xz_xyzz, to_0_x_xz_xzzz, to_0_x_xz_yyyy, to_0_x_xz_yyyz, to_0_x_xz_yyzz, to_0_x_xz_yzzz, to_0_x_xz_zzzz, to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xz_xxxx[k] = -4.0 * to_xz_xxx[k] + 2.0 * to_xz_xxxxx[k] * tke_0;

            to_0_x_xz_xxxy[k] = -3.0 * to_xz_xxy[k] + 2.0 * to_xz_xxxxy[k] * tke_0;

            to_0_x_xz_xxxz[k] = -3.0 * to_xz_xxz[k] + 2.0 * to_xz_xxxxz[k] * tke_0;

            to_0_x_xz_xxyy[k] = -2.0 * to_xz_xyy[k] + 2.0 * to_xz_xxxyy[k] * tke_0;

            to_0_x_xz_xxyz[k] = -2.0 * to_xz_xyz[k] + 2.0 * to_xz_xxxyz[k] * tke_0;

            to_0_x_xz_xxzz[k] = -2.0 * to_xz_xzz[k] + 2.0 * to_xz_xxxzz[k] * tke_0;

            to_0_x_xz_xyyy[k] = -to_xz_yyy[k] + 2.0 * to_xz_xxyyy[k] * tke_0;

            to_0_x_xz_xyyz[k] = -to_xz_yyz[k] + 2.0 * to_xz_xxyyz[k] * tke_0;

            to_0_x_xz_xyzz[k] = -to_xz_yzz[k] + 2.0 * to_xz_xxyzz[k] * tke_0;

            to_0_x_xz_xzzz[k] = -to_xz_zzz[k] + 2.0 * to_xz_xxzzz[k] * tke_0;

            to_0_x_xz_yyyy[k] = 2.0 * to_xz_xyyyy[k] * tke_0;

            to_0_x_xz_yyyz[k] = 2.0 * to_xz_xyyyz[k] * tke_0;

            to_0_x_xz_yyzz[k] = 2.0 * to_xz_xyyzz[k] * tke_0;

            to_0_x_xz_yzzz[k] = 2.0 * to_xz_xyzzz[k] * tke_0;

            to_0_x_xz_zzzz[k] = 2.0 * to_xz_xzzzz[k] * tke_0;
        }

        // Set up 45-60 components of targeted buffer : DG

        auto to_0_x_yy_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 45);

        auto to_0_x_yy_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 46);

        auto to_0_x_yy_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 47);

        auto to_0_x_yy_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 48);

        auto to_0_x_yy_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 49);

        auto to_0_x_yy_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 50);

        auto to_0_x_yy_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 51);

        auto to_0_x_yy_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 52);

        auto to_0_x_yy_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 53);

        auto to_0_x_yy_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 54);

        auto to_0_x_yy_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 55);

        auto to_0_x_yy_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 56);

        auto to_0_x_yy_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 57);

        auto to_0_x_yy_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 58);

        auto to_0_x_yy_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_x_yy_xxxx, to_0_x_yy_xxxy, to_0_x_yy_xxxz, to_0_x_yy_xxyy, to_0_x_yy_xxyz, to_0_x_yy_xxzz, to_0_x_yy_xyyy, to_0_x_yy_xyyz, to_0_x_yy_xyzz, to_0_x_yy_xzzz, to_0_x_yy_yyyy, to_0_x_yy_yyyz, to_0_x_yy_yyzz, to_0_x_yy_yzzz, to_0_x_yy_zzzz, to_yy_xxx, to_yy_xxxxx, to_yy_xxxxy, to_yy_xxxxz, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yy_xxxx[k] = -4.0 * to_yy_xxx[k] + 2.0 * to_yy_xxxxx[k] * tke_0;

            to_0_x_yy_xxxy[k] = -3.0 * to_yy_xxy[k] + 2.0 * to_yy_xxxxy[k] * tke_0;

            to_0_x_yy_xxxz[k] = -3.0 * to_yy_xxz[k] + 2.0 * to_yy_xxxxz[k] * tke_0;

            to_0_x_yy_xxyy[k] = -2.0 * to_yy_xyy[k] + 2.0 * to_yy_xxxyy[k] * tke_0;

            to_0_x_yy_xxyz[k] = -2.0 * to_yy_xyz[k] + 2.0 * to_yy_xxxyz[k] * tke_0;

            to_0_x_yy_xxzz[k] = -2.0 * to_yy_xzz[k] + 2.0 * to_yy_xxxzz[k] * tke_0;

            to_0_x_yy_xyyy[k] = -to_yy_yyy[k] + 2.0 * to_yy_xxyyy[k] * tke_0;

            to_0_x_yy_xyyz[k] = -to_yy_yyz[k] + 2.0 * to_yy_xxyyz[k] * tke_0;

            to_0_x_yy_xyzz[k] = -to_yy_yzz[k] + 2.0 * to_yy_xxyzz[k] * tke_0;

            to_0_x_yy_xzzz[k] = -to_yy_zzz[k] + 2.0 * to_yy_xxzzz[k] * tke_0;

            to_0_x_yy_yyyy[k] = 2.0 * to_yy_xyyyy[k] * tke_0;

            to_0_x_yy_yyyz[k] = 2.0 * to_yy_xyyyz[k] * tke_0;

            to_0_x_yy_yyzz[k] = 2.0 * to_yy_xyyzz[k] * tke_0;

            to_0_x_yy_yzzz[k] = 2.0 * to_yy_xyzzz[k] * tke_0;

            to_0_x_yy_zzzz[k] = 2.0 * to_yy_xzzzz[k] * tke_0;
        }

        // Set up 60-75 components of targeted buffer : DG

        auto to_0_x_yz_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 60);

        auto to_0_x_yz_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 61);

        auto to_0_x_yz_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 62);

        auto to_0_x_yz_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 63);

        auto to_0_x_yz_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 64);

        auto to_0_x_yz_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 65);

        auto to_0_x_yz_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 66);

        auto to_0_x_yz_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 67);

        auto to_0_x_yz_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 68);

        auto to_0_x_yz_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 69);

        auto to_0_x_yz_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 70);

        auto to_0_x_yz_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 71);

        auto to_0_x_yz_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 72);

        auto to_0_x_yz_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 73);

        auto to_0_x_yz_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 74);

        #pragma omp simd aligned(to_0_x_yz_xxxx, to_0_x_yz_xxxy, to_0_x_yz_xxxz, to_0_x_yz_xxyy, to_0_x_yz_xxyz, to_0_x_yz_xxzz, to_0_x_yz_xyyy, to_0_x_yz_xyyz, to_0_x_yz_xyzz, to_0_x_yz_xzzz, to_0_x_yz_yyyy, to_0_x_yz_yyyz, to_0_x_yz_yyzz, to_0_x_yz_yzzz, to_0_x_yz_zzzz, to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yz_xxxx[k] = -4.0 * to_yz_xxx[k] + 2.0 * to_yz_xxxxx[k] * tke_0;

            to_0_x_yz_xxxy[k] = -3.0 * to_yz_xxy[k] + 2.0 * to_yz_xxxxy[k] * tke_0;

            to_0_x_yz_xxxz[k] = -3.0 * to_yz_xxz[k] + 2.0 * to_yz_xxxxz[k] * tke_0;

            to_0_x_yz_xxyy[k] = -2.0 * to_yz_xyy[k] + 2.0 * to_yz_xxxyy[k] * tke_0;

            to_0_x_yz_xxyz[k] = -2.0 * to_yz_xyz[k] + 2.0 * to_yz_xxxyz[k] * tke_0;

            to_0_x_yz_xxzz[k] = -2.0 * to_yz_xzz[k] + 2.0 * to_yz_xxxzz[k] * tke_0;

            to_0_x_yz_xyyy[k] = -to_yz_yyy[k] + 2.0 * to_yz_xxyyy[k] * tke_0;

            to_0_x_yz_xyyz[k] = -to_yz_yyz[k] + 2.0 * to_yz_xxyyz[k] * tke_0;

            to_0_x_yz_xyzz[k] = -to_yz_yzz[k] + 2.0 * to_yz_xxyzz[k] * tke_0;

            to_0_x_yz_xzzz[k] = -to_yz_zzz[k] + 2.0 * to_yz_xxzzz[k] * tke_0;

            to_0_x_yz_yyyy[k] = 2.0 * to_yz_xyyyy[k] * tke_0;

            to_0_x_yz_yyyz[k] = 2.0 * to_yz_xyyyz[k] * tke_0;

            to_0_x_yz_yyzz[k] = 2.0 * to_yz_xyyzz[k] * tke_0;

            to_0_x_yz_yzzz[k] = 2.0 * to_yz_xyzzz[k] * tke_0;

            to_0_x_yz_zzzz[k] = 2.0 * to_yz_xzzzz[k] * tke_0;
        }

        // Set up 75-90 components of targeted buffer : DG

        auto to_0_x_zz_xxxx = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 75);

        auto to_0_x_zz_xxxy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 76);

        auto to_0_x_zz_xxxz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 77);

        auto to_0_x_zz_xxyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 78);

        auto to_0_x_zz_xxyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 79);

        auto to_0_x_zz_xxzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 80);

        auto to_0_x_zz_xyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 81);

        auto to_0_x_zz_xyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 82);

        auto to_0_x_zz_xyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 83);

        auto to_0_x_zz_xzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 84);

        auto to_0_x_zz_yyyy = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 85);

        auto to_0_x_zz_yyyz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 86);

        auto to_0_x_zz_yyzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 87);

        auto to_0_x_zz_yzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 88);

        auto to_0_x_zz_zzzz = pbuffer.data(idx_op_geom_001_dg + 0 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_x_zz_xxxx, to_0_x_zz_xxxy, to_0_x_zz_xxxz, to_0_x_zz_xxyy, to_0_x_zz_xxyz, to_0_x_zz_xxzz, to_0_x_zz_xyyy, to_0_x_zz_xyyz, to_0_x_zz_xyzz, to_0_x_zz_xzzz, to_0_x_zz_yyyy, to_0_x_zz_yyyz, to_0_x_zz_yyzz, to_0_x_zz_yzzz, to_0_x_zz_zzzz, to_zz_xxx, to_zz_xxxxx, to_zz_xxxxy, to_zz_xxxxz, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zz_xxxx[k] = -4.0 * to_zz_xxx[k] + 2.0 * to_zz_xxxxx[k] * tke_0;

            to_0_x_zz_xxxy[k] = -3.0 * to_zz_xxy[k] + 2.0 * to_zz_xxxxy[k] * tke_0;

            to_0_x_zz_xxxz[k] = -3.0 * to_zz_xxz[k] + 2.0 * to_zz_xxxxz[k] * tke_0;

            to_0_x_zz_xxyy[k] = -2.0 * to_zz_xyy[k] + 2.0 * to_zz_xxxyy[k] * tke_0;

            to_0_x_zz_xxyz[k] = -2.0 * to_zz_xyz[k] + 2.0 * to_zz_xxxyz[k] * tke_0;

            to_0_x_zz_xxzz[k] = -2.0 * to_zz_xzz[k] + 2.0 * to_zz_xxxzz[k] * tke_0;

            to_0_x_zz_xyyy[k] = -to_zz_yyy[k] + 2.0 * to_zz_xxyyy[k] * tke_0;

            to_0_x_zz_xyyz[k] = -to_zz_yyz[k] + 2.0 * to_zz_xxyyz[k] * tke_0;

            to_0_x_zz_xyzz[k] = -to_zz_yzz[k] + 2.0 * to_zz_xxyzz[k] * tke_0;

            to_0_x_zz_xzzz[k] = -to_zz_zzz[k] + 2.0 * to_zz_xxzzz[k] * tke_0;

            to_0_x_zz_yyyy[k] = 2.0 * to_zz_xyyyy[k] * tke_0;

            to_0_x_zz_yyyz[k] = 2.0 * to_zz_xyyyz[k] * tke_0;

            to_0_x_zz_yyzz[k] = 2.0 * to_zz_xyyzz[k] * tke_0;

            to_0_x_zz_yzzz[k] = 2.0 * to_zz_xyzzz[k] * tke_0;

            to_0_x_zz_zzzz[k] = 2.0 * to_zz_xzzzz[k] * tke_0;
        }

        // Set up 90-105 components of targeted buffer : DG

        auto to_0_y_xx_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 0);

        auto to_0_y_xx_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 1);

        auto to_0_y_xx_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 2);

        auto to_0_y_xx_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 3);

        auto to_0_y_xx_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 4);

        auto to_0_y_xx_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 5);

        auto to_0_y_xx_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 6);

        auto to_0_y_xx_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 7);

        auto to_0_y_xx_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 8);

        auto to_0_y_xx_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 9);

        auto to_0_y_xx_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 10);

        auto to_0_y_xx_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 11);

        auto to_0_y_xx_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 12);

        auto to_0_y_xx_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 13);

        auto to_0_y_xx_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 14);

        #pragma omp simd aligned(to_0_y_xx_xxxx, to_0_y_xx_xxxy, to_0_y_xx_xxxz, to_0_y_xx_xxyy, to_0_y_xx_xxyz, to_0_y_xx_xxzz, to_0_y_xx_xyyy, to_0_y_xx_xyyz, to_0_y_xx_xyzz, to_0_y_xx_xzzz, to_0_y_xx_yyyy, to_0_y_xx_yyyz, to_0_y_xx_yyzz, to_0_y_xx_yzzz, to_0_y_xx_zzzz, to_xx_xxx, to_xx_xxxxy, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_yyy, to_xx_yyyyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xx_xxxx[k] = 2.0 * to_xx_xxxxy[k] * tke_0;

            to_0_y_xx_xxxy[k] = -to_xx_xxx[k] + 2.0 * to_xx_xxxyy[k] * tke_0;

            to_0_y_xx_xxxz[k] = 2.0 * to_xx_xxxyz[k] * tke_0;

            to_0_y_xx_xxyy[k] = -2.0 * to_xx_xxy[k] + 2.0 * to_xx_xxyyy[k] * tke_0;

            to_0_y_xx_xxyz[k] = -to_xx_xxz[k] + 2.0 * to_xx_xxyyz[k] * tke_0;

            to_0_y_xx_xxzz[k] = 2.0 * to_xx_xxyzz[k] * tke_0;

            to_0_y_xx_xyyy[k] = -3.0 * to_xx_xyy[k] + 2.0 * to_xx_xyyyy[k] * tke_0;

            to_0_y_xx_xyyz[k] = -2.0 * to_xx_xyz[k] + 2.0 * to_xx_xyyyz[k] * tke_0;

            to_0_y_xx_xyzz[k] = -to_xx_xzz[k] + 2.0 * to_xx_xyyzz[k] * tke_0;

            to_0_y_xx_xzzz[k] = 2.0 * to_xx_xyzzz[k] * tke_0;

            to_0_y_xx_yyyy[k] = -4.0 * to_xx_yyy[k] + 2.0 * to_xx_yyyyy[k] * tke_0;

            to_0_y_xx_yyyz[k] = -3.0 * to_xx_yyz[k] + 2.0 * to_xx_yyyyz[k] * tke_0;

            to_0_y_xx_yyzz[k] = -2.0 * to_xx_yzz[k] + 2.0 * to_xx_yyyzz[k] * tke_0;

            to_0_y_xx_yzzz[k] = -to_xx_zzz[k] + 2.0 * to_xx_yyzzz[k] * tke_0;

            to_0_y_xx_zzzz[k] = 2.0 * to_xx_yzzzz[k] * tke_0;
        }

        // Set up 105-120 components of targeted buffer : DG

        auto to_0_y_xy_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 15);

        auto to_0_y_xy_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 16);

        auto to_0_y_xy_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 17);

        auto to_0_y_xy_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 18);

        auto to_0_y_xy_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 19);

        auto to_0_y_xy_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 20);

        auto to_0_y_xy_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 21);

        auto to_0_y_xy_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 22);

        auto to_0_y_xy_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 23);

        auto to_0_y_xy_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 24);

        auto to_0_y_xy_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 25);

        auto to_0_y_xy_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 26);

        auto to_0_y_xy_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 27);

        auto to_0_y_xy_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 28);

        auto to_0_y_xy_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_y_xy_xxxx, to_0_y_xy_xxxy, to_0_y_xy_xxxz, to_0_y_xy_xxyy, to_0_y_xy_xxyz, to_0_y_xy_xxzz, to_0_y_xy_xyyy, to_0_y_xy_xyyz, to_0_y_xy_xyzz, to_0_y_xy_xzzz, to_0_y_xy_yyyy, to_0_y_xy_yyyz, to_0_y_xy_yyzz, to_0_y_xy_yzzz, to_0_y_xy_zzzz, to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xy_xxxx[k] = 2.0 * to_xy_xxxxy[k] * tke_0;

            to_0_y_xy_xxxy[k] = -to_xy_xxx[k] + 2.0 * to_xy_xxxyy[k] * tke_0;

            to_0_y_xy_xxxz[k] = 2.0 * to_xy_xxxyz[k] * tke_0;

            to_0_y_xy_xxyy[k] = -2.0 * to_xy_xxy[k] + 2.0 * to_xy_xxyyy[k] * tke_0;

            to_0_y_xy_xxyz[k] = -to_xy_xxz[k] + 2.0 * to_xy_xxyyz[k] * tke_0;

            to_0_y_xy_xxzz[k] = 2.0 * to_xy_xxyzz[k] * tke_0;

            to_0_y_xy_xyyy[k] = -3.0 * to_xy_xyy[k] + 2.0 * to_xy_xyyyy[k] * tke_0;

            to_0_y_xy_xyyz[k] = -2.0 * to_xy_xyz[k] + 2.0 * to_xy_xyyyz[k] * tke_0;

            to_0_y_xy_xyzz[k] = -to_xy_xzz[k] + 2.0 * to_xy_xyyzz[k] * tke_0;

            to_0_y_xy_xzzz[k] = 2.0 * to_xy_xyzzz[k] * tke_0;

            to_0_y_xy_yyyy[k] = -4.0 * to_xy_yyy[k] + 2.0 * to_xy_yyyyy[k] * tke_0;

            to_0_y_xy_yyyz[k] = -3.0 * to_xy_yyz[k] + 2.0 * to_xy_yyyyz[k] * tke_0;

            to_0_y_xy_yyzz[k] = -2.0 * to_xy_yzz[k] + 2.0 * to_xy_yyyzz[k] * tke_0;

            to_0_y_xy_yzzz[k] = -to_xy_zzz[k] + 2.0 * to_xy_yyzzz[k] * tke_0;

            to_0_y_xy_zzzz[k] = 2.0 * to_xy_yzzzz[k] * tke_0;
        }

        // Set up 120-135 components of targeted buffer : DG

        auto to_0_y_xz_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 30);

        auto to_0_y_xz_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 31);

        auto to_0_y_xz_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 32);

        auto to_0_y_xz_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 33);

        auto to_0_y_xz_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 34);

        auto to_0_y_xz_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 35);

        auto to_0_y_xz_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 36);

        auto to_0_y_xz_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 37);

        auto to_0_y_xz_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 38);

        auto to_0_y_xz_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 39);

        auto to_0_y_xz_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 40);

        auto to_0_y_xz_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 41);

        auto to_0_y_xz_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 42);

        auto to_0_y_xz_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 43);

        auto to_0_y_xz_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 44);

        #pragma omp simd aligned(to_0_y_xz_xxxx, to_0_y_xz_xxxy, to_0_y_xz_xxxz, to_0_y_xz_xxyy, to_0_y_xz_xxyz, to_0_y_xz_xxzz, to_0_y_xz_xyyy, to_0_y_xz_xyyz, to_0_y_xz_xyzz, to_0_y_xz_xzzz, to_0_y_xz_yyyy, to_0_y_xz_yyyz, to_0_y_xz_yyzz, to_0_y_xz_yzzz, to_0_y_xz_zzzz, to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xz_xxxx[k] = 2.0 * to_xz_xxxxy[k] * tke_0;

            to_0_y_xz_xxxy[k] = -to_xz_xxx[k] + 2.0 * to_xz_xxxyy[k] * tke_0;

            to_0_y_xz_xxxz[k] = 2.0 * to_xz_xxxyz[k] * tke_0;

            to_0_y_xz_xxyy[k] = -2.0 * to_xz_xxy[k] + 2.0 * to_xz_xxyyy[k] * tke_0;

            to_0_y_xz_xxyz[k] = -to_xz_xxz[k] + 2.0 * to_xz_xxyyz[k] * tke_0;

            to_0_y_xz_xxzz[k] = 2.0 * to_xz_xxyzz[k] * tke_0;

            to_0_y_xz_xyyy[k] = -3.0 * to_xz_xyy[k] + 2.0 * to_xz_xyyyy[k] * tke_0;

            to_0_y_xz_xyyz[k] = -2.0 * to_xz_xyz[k] + 2.0 * to_xz_xyyyz[k] * tke_0;

            to_0_y_xz_xyzz[k] = -to_xz_xzz[k] + 2.0 * to_xz_xyyzz[k] * tke_0;

            to_0_y_xz_xzzz[k] = 2.0 * to_xz_xyzzz[k] * tke_0;

            to_0_y_xz_yyyy[k] = -4.0 * to_xz_yyy[k] + 2.0 * to_xz_yyyyy[k] * tke_0;

            to_0_y_xz_yyyz[k] = -3.0 * to_xz_yyz[k] + 2.0 * to_xz_yyyyz[k] * tke_0;

            to_0_y_xz_yyzz[k] = -2.0 * to_xz_yzz[k] + 2.0 * to_xz_yyyzz[k] * tke_0;

            to_0_y_xz_yzzz[k] = -to_xz_zzz[k] + 2.0 * to_xz_yyzzz[k] * tke_0;

            to_0_y_xz_zzzz[k] = 2.0 * to_xz_yzzzz[k] * tke_0;
        }

        // Set up 135-150 components of targeted buffer : DG

        auto to_0_y_yy_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 45);

        auto to_0_y_yy_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 46);

        auto to_0_y_yy_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 47);

        auto to_0_y_yy_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 48);

        auto to_0_y_yy_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 49);

        auto to_0_y_yy_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 50);

        auto to_0_y_yy_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 51);

        auto to_0_y_yy_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 52);

        auto to_0_y_yy_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 53);

        auto to_0_y_yy_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 54);

        auto to_0_y_yy_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 55);

        auto to_0_y_yy_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 56);

        auto to_0_y_yy_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 57);

        auto to_0_y_yy_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 58);

        auto to_0_y_yy_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_y_yy_xxxx, to_0_y_yy_xxxy, to_0_y_yy_xxxz, to_0_y_yy_xxyy, to_0_y_yy_xxyz, to_0_y_yy_xxzz, to_0_y_yy_xyyy, to_0_y_yy_xyyz, to_0_y_yy_xyzz, to_0_y_yy_xzzz, to_0_y_yy_yyyy, to_0_y_yy_yyyz, to_0_y_yy_yyzz, to_0_y_yy_yzzz, to_0_y_yy_zzzz, to_yy_xxx, to_yy_xxxxy, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_yyy, to_yy_yyyyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yy_xxxx[k] = 2.0 * to_yy_xxxxy[k] * tke_0;

            to_0_y_yy_xxxy[k] = -to_yy_xxx[k] + 2.0 * to_yy_xxxyy[k] * tke_0;

            to_0_y_yy_xxxz[k] = 2.0 * to_yy_xxxyz[k] * tke_0;

            to_0_y_yy_xxyy[k] = -2.0 * to_yy_xxy[k] + 2.0 * to_yy_xxyyy[k] * tke_0;

            to_0_y_yy_xxyz[k] = -to_yy_xxz[k] + 2.0 * to_yy_xxyyz[k] * tke_0;

            to_0_y_yy_xxzz[k] = 2.0 * to_yy_xxyzz[k] * tke_0;

            to_0_y_yy_xyyy[k] = -3.0 * to_yy_xyy[k] + 2.0 * to_yy_xyyyy[k] * tke_0;

            to_0_y_yy_xyyz[k] = -2.0 * to_yy_xyz[k] + 2.0 * to_yy_xyyyz[k] * tke_0;

            to_0_y_yy_xyzz[k] = -to_yy_xzz[k] + 2.0 * to_yy_xyyzz[k] * tke_0;

            to_0_y_yy_xzzz[k] = 2.0 * to_yy_xyzzz[k] * tke_0;

            to_0_y_yy_yyyy[k] = -4.0 * to_yy_yyy[k] + 2.0 * to_yy_yyyyy[k] * tke_0;

            to_0_y_yy_yyyz[k] = -3.0 * to_yy_yyz[k] + 2.0 * to_yy_yyyyz[k] * tke_0;

            to_0_y_yy_yyzz[k] = -2.0 * to_yy_yzz[k] + 2.0 * to_yy_yyyzz[k] * tke_0;

            to_0_y_yy_yzzz[k] = -to_yy_zzz[k] + 2.0 * to_yy_yyzzz[k] * tke_0;

            to_0_y_yy_zzzz[k] = 2.0 * to_yy_yzzzz[k] * tke_0;
        }

        // Set up 150-165 components of targeted buffer : DG

        auto to_0_y_yz_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 60);

        auto to_0_y_yz_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 61);

        auto to_0_y_yz_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 62);

        auto to_0_y_yz_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 63);

        auto to_0_y_yz_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 64);

        auto to_0_y_yz_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 65);

        auto to_0_y_yz_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 66);

        auto to_0_y_yz_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 67);

        auto to_0_y_yz_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 68);

        auto to_0_y_yz_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 69);

        auto to_0_y_yz_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 70);

        auto to_0_y_yz_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 71);

        auto to_0_y_yz_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 72);

        auto to_0_y_yz_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 73);

        auto to_0_y_yz_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 74);

        #pragma omp simd aligned(to_0_y_yz_xxxx, to_0_y_yz_xxxy, to_0_y_yz_xxxz, to_0_y_yz_xxyy, to_0_y_yz_xxyz, to_0_y_yz_xxzz, to_0_y_yz_xyyy, to_0_y_yz_xyyz, to_0_y_yz_xyzz, to_0_y_yz_xzzz, to_0_y_yz_yyyy, to_0_y_yz_yyyz, to_0_y_yz_yyzz, to_0_y_yz_yzzz, to_0_y_yz_zzzz, to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yz_xxxx[k] = 2.0 * to_yz_xxxxy[k] * tke_0;

            to_0_y_yz_xxxy[k] = -to_yz_xxx[k] + 2.0 * to_yz_xxxyy[k] * tke_0;

            to_0_y_yz_xxxz[k] = 2.0 * to_yz_xxxyz[k] * tke_0;

            to_0_y_yz_xxyy[k] = -2.0 * to_yz_xxy[k] + 2.0 * to_yz_xxyyy[k] * tke_0;

            to_0_y_yz_xxyz[k] = -to_yz_xxz[k] + 2.0 * to_yz_xxyyz[k] * tke_0;

            to_0_y_yz_xxzz[k] = 2.0 * to_yz_xxyzz[k] * tke_0;

            to_0_y_yz_xyyy[k] = -3.0 * to_yz_xyy[k] + 2.0 * to_yz_xyyyy[k] * tke_0;

            to_0_y_yz_xyyz[k] = -2.0 * to_yz_xyz[k] + 2.0 * to_yz_xyyyz[k] * tke_0;

            to_0_y_yz_xyzz[k] = -to_yz_xzz[k] + 2.0 * to_yz_xyyzz[k] * tke_0;

            to_0_y_yz_xzzz[k] = 2.0 * to_yz_xyzzz[k] * tke_0;

            to_0_y_yz_yyyy[k] = -4.0 * to_yz_yyy[k] + 2.0 * to_yz_yyyyy[k] * tke_0;

            to_0_y_yz_yyyz[k] = -3.0 * to_yz_yyz[k] + 2.0 * to_yz_yyyyz[k] * tke_0;

            to_0_y_yz_yyzz[k] = -2.0 * to_yz_yzz[k] + 2.0 * to_yz_yyyzz[k] * tke_0;

            to_0_y_yz_yzzz[k] = -to_yz_zzz[k] + 2.0 * to_yz_yyzzz[k] * tke_0;

            to_0_y_yz_zzzz[k] = 2.0 * to_yz_yzzzz[k] * tke_0;
        }

        // Set up 165-180 components of targeted buffer : DG

        auto to_0_y_zz_xxxx = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 75);

        auto to_0_y_zz_xxxy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 76);

        auto to_0_y_zz_xxxz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 77);

        auto to_0_y_zz_xxyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 78);

        auto to_0_y_zz_xxyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 79);

        auto to_0_y_zz_xxzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 80);

        auto to_0_y_zz_xyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 81);

        auto to_0_y_zz_xyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 82);

        auto to_0_y_zz_xyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 83);

        auto to_0_y_zz_xzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 84);

        auto to_0_y_zz_yyyy = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 85);

        auto to_0_y_zz_yyyz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 86);

        auto to_0_y_zz_yyzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 87);

        auto to_0_y_zz_yzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 88);

        auto to_0_y_zz_zzzz = pbuffer.data(idx_op_geom_001_dg + 1 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_y_zz_xxxx, to_0_y_zz_xxxy, to_0_y_zz_xxxz, to_0_y_zz_xxyy, to_0_y_zz_xxyz, to_0_y_zz_xxzz, to_0_y_zz_xyyy, to_0_y_zz_xyyz, to_0_y_zz_xyzz, to_0_y_zz_xzzz, to_0_y_zz_yyyy, to_0_y_zz_yyyz, to_0_y_zz_yyzz, to_0_y_zz_yzzz, to_0_y_zz_zzzz, to_zz_xxx, to_zz_xxxxy, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_yyy, to_zz_yyyyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zz_xxxx[k] = 2.0 * to_zz_xxxxy[k] * tke_0;

            to_0_y_zz_xxxy[k] = -to_zz_xxx[k] + 2.0 * to_zz_xxxyy[k] * tke_0;

            to_0_y_zz_xxxz[k] = 2.0 * to_zz_xxxyz[k] * tke_0;

            to_0_y_zz_xxyy[k] = -2.0 * to_zz_xxy[k] + 2.0 * to_zz_xxyyy[k] * tke_0;

            to_0_y_zz_xxyz[k] = -to_zz_xxz[k] + 2.0 * to_zz_xxyyz[k] * tke_0;

            to_0_y_zz_xxzz[k] = 2.0 * to_zz_xxyzz[k] * tke_0;

            to_0_y_zz_xyyy[k] = -3.0 * to_zz_xyy[k] + 2.0 * to_zz_xyyyy[k] * tke_0;

            to_0_y_zz_xyyz[k] = -2.0 * to_zz_xyz[k] + 2.0 * to_zz_xyyyz[k] * tke_0;

            to_0_y_zz_xyzz[k] = -to_zz_xzz[k] + 2.0 * to_zz_xyyzz[k] * tke_0;

            to_0_y_zz_xzzz[k] = 2.0 * to_zz_xyzzz[k] * tke_0;

            to_0_y_zz_yyyy[k] = -4.0 * to_zz_yyy[k] + 2.0 * to_zz_yyyyy[k] * tke_0;

            to_0_y_zz_yyyz[k] = -3.0 * to_zz_yyz[k] + 2.0 * to_zz_yyyyz[k] * tke_0;

            to_0_y_zz_yyzz[k] = -2.0 * to_zz_yzz[k] + 2.0 * to_zz_yyyzz[k] * tke_0;

            to_0_y_zz_yzzz[k] = -to_zz_zzz[k] + 2.0 * to_zz_yyzzz[k] * tke_0;

            to_0_y_zz_zzzz[k] = 2.0 * to_zz_yzzzz[k] * tke_0;
        }

        // Set up 180-195 components of targeted buffer : DG

        auto to_0_z_xx_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 0);

        auto to_0_z_xx_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 1);

        auto to_0_z_xx_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 2);

        auto to_0_z_xx_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 3);

        auto to_0_z_xx_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 4);

        auto to_0_z_xx_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 5);

        auto to_0_z_xx_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 6);

        auto to_0_z_xx_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 7);

        auto to_0_z_xx_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 8);

        auto to_0_z_xx_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 9);

        auto to_0_z_xx_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 10);

        auto to_0_z_xx_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 11);

        auto to_0_z_xx_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 12);

        auto to_0_z_xx_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 13);

        auto to_0_z_xx_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 14);

        #pragma omp simd aligned(to_0_z_xx_xxxx, to_0_z_xx_xxxy, to_0_z_xx_xxxz, to_0_z_xx_xxyy, to_0_z_xx_xxyz, to_0_z_xx_xxzz, to_0_z_xx_xyyy, to_0_z_xx_xyyz, to_0_z_xx_xyzz, to_0_z_xx_xzzz, to_0_z_xx_yyyy, to_0_z_xx_yyyz, to_0_z_xx_yyzz, to_0_z_xx_yzzz, to_0_z_xx_zzzz, to_xx_xxx, to_xx_xxxxz, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xx_xxxx[k] = 2.0 * to_xx_xxxxz[k] * tke_0;

            to_0_z_xx_xxxy[k] = 2.0 * to_xx_xxxyz[k] * tke_0;

            to_0_z_xx_xxxz[k] = -to_xx_xxx[k] + 2.0 * to_xx_xxxzz[k] * tke_0;

            to_0_z_xx_xxyy[k] = 2.0 * to_xx_xxyyz[k] * tke_0;

            to_0_z_xx_xxyz[k] = -to_xx_xxy[k] + 2.0 * to_xx_xxyzz[k] * tke_0;

            to_0_z_xx_xxzz[k] = -2.0 * to_xx_xxz[k] + 2.0 * to_xx_xxzzz[k] * tke_0;

            to_0_z_xx_xyyy[k] = 2.0 * to_xx_xyyyz[k] * tke_0;

            to_0_z_xx_xyyz[k] = -to_xx_xyy[k] + 2.0 * to_xx_xyyzz[k] * tke_0;

            to_0_z_xx_xyzz[k] = -2.0 * to_xx_xyz[k] + 2.0 * to_xx_xyzzz[k] * tke_0;

            to_0_z_xx_xzzz[k] = -3.0 * to_xx_xzz[k] + 2.0 * to_xx_xzzzz[k] * tke_0;

            to_0_z_xx_yyyy[k] = 2.0 * to_xx_yyyyz[k] * tke_0;

            to_0_z_xx_yyyz[k] = -to_xx_yyy[k] + 2.0 * to_xx_yyyzz[k] * tke_0;

            to_0_z_xx_yyzz[k] = -2.0 * to_xx_yyz[k] + 2.0 * to_xx_yyzzz[k] * tke_0;

            to_0_z_xx_yzzz[k] = -3.0 * to_xx_yzz[k] + 2.0 * to_xx_yzzzz[k] * tke_0;

            to_0_z_xx_zzzz[k] = -4.0 * to_xx_zzz[k] + 2.0 * to_xx_zzzzz[k] * tke_0;
        }

        // Set up 195-210 components of targeted buffer : DG

        auto to_0_z_xy_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 15);

        auto to_0_z_xy_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 16);

        auto to_0_z_xy_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 17);

        auto to_0_z_xy_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 18);

        auto to_0_z_xy_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 19);

        auto to_0_z_xy_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 20);

        auto to_0_z_xy_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 21);

        auto to_0_z_xy_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 22);

        auto to_0_z_xy_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 23);

        auto to_0_z_xy_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 24);

        auto to_0_z_xy_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 25);

        auto to_0_z_xy_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 26);

        auto to_0_z_xy_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 27);

        auto to_0_z_xy_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 28);

        auto to_0_z_xy_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_z_xy_xxxx, to_0_z_xy_xxxy, to_0_z_xy_xxxz, to_0_z_xy_xxyy, to_0_z_xy_xxyz, to_0_z_xy_xxzz, to_0_z_xy_xyyy, to_0_z_xy_xyyz, to_0_z_xy_xyzz, to_0_z_xy_xzzz, to_0_z_xy_yyyy, to_0_z_xy_yyyz, to_0_z_xy_yyzz, to_0_z_xy_yzzz, to_0_z_xy_zzzz, to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xy_xxxx[k] = 2.0 * to_xy_xxxxz[k] * tke_0;

            to_0_z_xy_xxxy[k] = 2.0 * to_xy_xxxyz[k] * tke_0;

            to_0_z_xy_xxxz[k] = -to_xy_xxx[k] + 2.0 * to_xy_xxxzz[k] * tke_0;

            to_0_z_xy_xxyy[k] = 2.0 * to_xy_xxyyz[k] * tke_0;

            to_0_z_xy_xxyz[k] = -to_xy_xxy[k] + 2.0 * to_xy_xxyzz[k] * tke_0;

            to_0_z_xy_xxzz[k] = -2.0 * to_xy_xxz[k] + 2.0 * to_xy_xxzzz[k] * tke_0;

            to_0_z_xy_xyyy[k] = 2.0 * to_xy_xyyyz[k] * tke_0;

            to_0_z_xy_xyyz[k] = -to_xy_xyy[k] + 2.0 * to_xy_xyyzz[k] * tke_0;

            to_0_z_xy_xyzz[k] = -2.0 * to_xy_xyz[k] + 2.0 * to_xy_xyzzz[k] * tke_0;

            to_0_z_xy_xzzz[k] = -3.0 * to_xy_xzz[k] + 2.0 * to_xy_xzzzz[k] * tke_0;

            to_0_z_xy_yyyy[k] = 2.0 * to_xy_yyyyz[k] * tke_0;

            to_0_z_xy_yyyz[k] = -to_xy_yyy[k] + 2.0 * to_xy_yyyzz[k] * tke_0;

            to_0_z_xy_yyzz[k] = -2.0 * to_xy_yyz[k] + 2.0 * to_xy_yyzzz[k] * tke_0;

            to_0_z_xy_yzzz[k] = -3.0 * to_xy_yzz[k] + 2.0 * to_xy_yzzzz[k] * tke_0;

            to_0_z_xy_zzzz[k] = -4.0 * to_xy_zzz[k] + 2.0 * to_xy_zzzzz[k] * tke_0;
        }

        // Set up 210-225 components of targeted buffer : DG

        auto to_0_z_xz_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 30);

        auto to_0_z_xz_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 31);

        auto to_0_z_xz_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 32);

        auto to_0_z_xz_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 33);

        auto to_0_z_xz_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 34);

        auto to_0_z_xz_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 35);

        auto to_0_z_xz_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 36);

        auto to_0_z_xz_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 37);

        auto to_0_z_xz_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 38);

        auto to_0_z_xz_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 39);

        auto to_0_z_xz_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 40);

        auto to_0_z_xz_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 41);

        auto to_0_z_xz_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 42);

        auto to_0_z_xz_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 43);

        auto to_0_z_xz_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 44);

        #pragma omp simd aligned(to_0_z_xz_xxxx, to_0_z_xz_xxxy, to_0_z_xz_xxxz, to_0_z_xz_xxyy, to_0_z_xz_xxyz, to_0_z_xz_xxzz, to_0_z_xz_xyyy, to_0_z_xz_xyyz, to_0_z_xz_xyzz, to_0_z_xz_xzzz, to_0_z_xz_yyyy, to_0_z_xz_yyyz, to_0_z_xz_yyzz, to_0_z_xz_yzzz, to_0_z_xz_zzzz, to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xz_xxxx[k] = 2.0 * to_xz_xxxxz[k] * tke_0;

            to_0_z_xz_xxxy[k] = 2.0 * to_xz_xxxyz[k] * tke_0;

            to_0_z_xz_xxxz[k] = -to_xz_xxx[k] + 2.0 * to_xz_xxxzz[k] * tke_0;

            to_0_z_xz_xxyy[k] = 2.0 * to_xz_xxyyz[k] * tke_0;

            to_0_z_xz_xxyz[k] = -to_xz_xxy[k] + 2.0 * to_xz_xxyzz[k] * tke_0;

            to_0_z_xz_xxzz[k] = -2.0 * to_xz_xxz[k] + 2.0 * to_xz_xxzzz[k] * tke_0;

            to_0_z_xz_xyyy[k] = 2.0 * to_xz_xyyyz[k] * tke_0;

            to_0_z_xz_xyyz[k] = -to_xz_xyy[k] + 2.0 * to_xz_xyyzz[k] * tke_0;

            to_0_z_xz_xyzz[k] = -2.0 * to_xz_xyz[k] + 2.0 * to_xz_xyzzz[k] * tke_0;

            to_0_z_xz_xzzz[k] = -3.0 * to_xz_xzz[k] + 2.0 * to_xz_xzzzz[k] * tke_0;

            to_0_z_xz_yyyy[k] = 2.0 * to_xz_yyyyz[k] * tke_0;

            to_0_z_xz_yyyz[k] = -to_xz_yyy[k] + 2.0 * to_xz_yyyzz[k] * tke_0;

            to_0_z_xz_yyzz[k] = -2.0 * to_xz_yyz[k] + 2.0 * to_xz_yyzzz[k] * tke_0;

            to_0_z_xz_yzzz[k] = -3.0 * to_xz_yzz[k] + 2.0 * to_xz_yzzzz[k] * tke_0;

            to_0_z_xz_zzzz[k] = -4.0 * to_xz_zzz[k] + 2.0 * to_xz_zzzzz[k] * tke_0;
        }

        // Set up 225-240 components of targeted buffer : DG

        auto to_0_z_yy_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 45);

        auto to_0_z_yy_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 46);

        auto to_0_z_yy_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 47);

        auto to_0_z_yy_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 48);

        auto to_0_z_yy_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 49);

        auto to_0_z_yy_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 50);

        auto to_0_z_yy_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 51);

        auto to_0_z_yy_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 52);

        auto to_0_z_yy_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 53);

        auto to_0_z_yy_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 54);

        auto to_0_z_yy_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 55);

        auto to_0_z_yy_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 56);

        auto to_0_z_yy_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 57);

        auto to_0_z_yy_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 58);

        auto to_0_z_yy_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_z_yy_xxxx, to_0_z_yy_xxxy, to_0_z_yy_xxxz, to_0_z_yy_xxyy, to_0_z_yy_xxyz, to_0_z_yy_xxzz, to_0_z_yy_xyyy, to_0_z_yy_xyyz, to_0_z_yy_xyzz, to_0_z_yy_xzzz, to_0_z_yy_yyyy, to_0_z_yy_yyyz, to_0_z_yy_yyzz, to_0_z_yy_yzzz, to_0_z_yy_zzzz, to_yy_xxx, to_yy_xxxxz, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yy_xxxx[k] = 2.0 * to_yy_xxxxz[k] * tke_0;

            to_0_z_yy_xxxy[k] = 2.0 * to_yy_xxxyz[k] * tke_0;

            to_0_z_yy_xxxz[k] = -to_yy_xxx[k] + 2.0 * to_yy_xxxzz[k] * tke_0;

            to_0_z_yy_xxyy[k] = 2.0 * to_yy_xxyyz[k] * tke_0;

            to_0_z_yy_xxyz[k] = -to_yy_xxy[k] + 2.0 * to_yy_xxyzz[k] * tke_0;

            to_0_z_yy_xxzz[k] = -2.0 * to_yy_xxz[k] + 2.0 * to_yy_xxzzz[k] * tke_0;

            to_0_z_yy_xyyy[k] = 2.0 * to_yy_xyyyz[k] * tke_0;

            to_0_z_yy_xyyz[k] = -to_yy_xyy[k] + 2.0 * to_yy_xyyzz[k] * tke_0;

            to_0_z_yy_xyzz[k] = -2.0 * to_yy_xyz[k] + 2.0 * to_yy_xyzzz[k] * tke_0;

            to_0_z_yy_xzzz[k] = -3.0 * to_yy_xzz[k] + 2.0 * to_yy_xzzzz[k] * tke_0;

            to_0_z_yy_yyyy[k] = 2.0 * to_yy_yyyyz[k] * tke_0;

            to_0_z_yy_yyyz[k] = -to_yy_yyy[k] + 2.0 * to_yy_yyyzz[k] * tke_0;

            to_0_z_yy_yyzz[k] = -2.0 * to_yy_yyz[k] + 2.0 * to_yy_yyzzz[k] * tke_0;

            to_0_z_yy_yzzz[k] = -3.0 * to_yy_yzz[k] + 2.0 * to_yy_yzzzz[k] * tke_0;

            to_0_z_yy_zzzz[k] = -4.0 * to_yy_zzz[k] + 2.0 * to_yy_zzzzz[k] * tke_0;
        }

        // Set up 240-255 components of targeted buffer : DG

        auto to_0_z_yz_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 60);

        auto to_0_z_yz_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 61);

        auto to_0_z_yz_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 62);

        auto to_0_z_yz_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 63);

        auto to_0_z_yz_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 64);

        auto to_0_z_yz_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 65);

        auto to_0_z_yz_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 66);

        auto to_0_z_yz_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 67);

        auto to_0_z_yz_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 68);

        auto to_0_z_yz_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 69);

        auto to_0_z_yz_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 70);

        auto to_0_z_yz_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 71);

        auto to_0_z_yz_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 72);

        auto to_0_z_yz_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 73);

        auto to_0_z_yz_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 74);

        #pragma omp simd aligned(to_0_z_yz_xxxx, to_0_z_yz_xxxy, to_0_z_yz_xxxz, to_0_z_yz_xxyy, to_0_z_yz_xxyz, to_0_z_yz_xxzz, to_0_z_yz_xyyy, to_0_z_yz_xyyz, to_0_z_yz_xyzz, to_0_z_yz_xzzz, to_0_z_yz_yyyy, to_0_z_yz_yyyz, to_0_z_yz_yyzz, to_0_z_yz_yzzz, to_0_z_yz_zzzz, to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yz_xxxx[k] = 2.0 * to_yz_xxxxz[k] * tke_0;

            to_0_z_yz_xxxy[k] = 2.0 * to_yz_xxxyz[k] * tke_0;

            to_0_z_yz_xxxz[k] = -to_yz_xxx[k] + 2.0 * to_yz_xxxzz[k] * tke_0;

            to_0_z_yz_xxyy[k] = 2.0 * to_yz_xxyyz[k] * tke_0;

            to_0_z_yz_xxyz[k] = -to_yz_xxy[k] + 2.0 * to_yz_xxyzz[k] * tke_0;

            to_0_z_yz_xxzz[k] = -2.0 * to_yz_xxz[k] + 2.0 * to_yz_xxzzz[k] * tke_0;

            to_0_z_yz_xyyy[k] = 2.0 * to_yz_xyyyz[k] * tke_0;

            to_0_z_yz_xyyz[k] = -to_yz_xyy[k] + 2.0 * to_yz_xyyzz[k] * tke_0;

            to_0_z_yz_xyzz[k] = -2.0 * to_yz_xyz[k] + 2.0 * to_yz_xyzzz[k] * tke_0;

            to_0_z_yz_xzzz[k] = -3.0 * to_yz_xzz[k] + 2.0 * to_yz_xzzzz[k] * tke_0;

            to_0_z_yz_yyyy[k] = 2.0 * to_yz_yyyyz[k] * tke_0;

            to_0_z_yz_yyyz[k] = -to_yz_yyy[k] + 2.0 * to_yz_yyyzz[k] * tke_0;

            to_0_z_yz_yyzz[k] = -2.0 * to_yz_yyz[k] + 2.0 * to_yz_yyzzz[k] * tke_0;

            to_0_z_yz_yzzz[k] = -3.0 * to_yz_yzz[k] + 2.0 * to_yz_yzzzz[k] * tke_0;

            to_0_z_yz_zzzz[k] = -4.0 * to_yz_zzz[k] + 2.0 * to_yz_zzzzz[k] * tke_0;
        }

        // Set up 255-270 components of targeted buffer : DG

        auto to_0_z_zz_xxxx = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 75);

        auto to_0_z_zz_xxxy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 76);

        auto to_0_z_zz_xxxz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 77);

        auto to_0_z_zz_xxyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 78);

        auto to_0_z_zz_xxyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 79);

        auto to_0_z_zz_xxzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 80);

        auto to_0_z_zz_xyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 81);

        auto to_0_z_zz_xyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 82);

        auto to_0_z_zz_xyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 83);

        auto to_0_z_zz_xzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 84);

        auto to_0_z_zz_yyyy = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 85);

        auto to_0_z_zz_yyyz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 86);

        auto to_0_z_zz_yyzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 87);

        auto to_0_z_zz_yzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 88);

        auto to_0_z_zz_zzzz = pbuffer.data(idx_op_geom_001_dg + 2 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_z_zz_xxxx, to_0_z_zz_xxxy, to_0_z_zz_xxxz, to_0_z_zz_xxyy, to_0_z_zz_xxyz, to_0_z_zz_xxzz, to_0_z_zz_xyyy, to_0_z_zz_xyyz, to_0_z_zz_xyzz, to_0_z_zz_xzzz, to_0_z_zz_yyyy, to_0_z_zz_yyyz, to_0_z_zz_yyzz, to_0_z_zz_yzzz, to_0_z_zz_zzzz, to_zz_xxx, to_zz_xxxxz, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zz_xxxx[k] = 2.0 * to_zz_xxxxz[k] * tke_0;

            to_0_z_zz_xxxy[k] = 2.0 * to_zz_xxxyz[k] * tke_0;

            to_0_z_zz_xxxz[k] = -to_zz_xxx[k] + 2.0 * to_zz_xxxzz[k] * tke_0;

            to_0_z_zz_xxyy[k] = 2.0 * to_zz_xxyyz[k] * tke_0;

            to_0_z_zz_xxyz[k] = -to_zz_xxy[k] + 2.0 * to_zz_xxyzz[k] * tke_0;

            to_0_z_zz_xxzz[k] = -2.0 * to_zz_xxz[k] + 2.0 * to_zz_xxzzz[k] * tke_0;

            to_0_z_zz_xyyy[k] = 2.0 * to_zz_xyyyz[k] * tke_0;

            to_0_z_zz_xyyz[k] = -to_zz_xyy[k] + 2.0 * to_zz_xyyzz[k] * tke_0;

            to_0_z_zz_xyzz[k] = -2.0 * to_zz_xyz[k] + 2.0 * to_zz_xyzzz[k] * tke_0;

            to_0_z_zz_xzzz[k] = -3.0 * to_zz_xzz[k] + 2.0 * to_zz_xzzzz[k] * tke_0;

            to_0_z_zz_yyyy[k] = 2.0 * to_zz_yyyyz[k] * tke_0;

            to_0_z_zz_yyyz[k] = -to_zz_yyy[k] + 2.0 * to_zz_yyyzz[k] * tke_0;

            to_0_z_zz_yyzz[k] = -2.0 * to_zz_yyz[k] + 2.0 * to_zz_yyzzz[k] * tke_0;

            to_0_z_zz_yzzz[k] = -3.0 * to_zz_yzz[k] + 2.0 * to_zz_yzzzz[k] * tke_0;

            to_0_z_zz_zzzz[k] = -4.0 * to_zz_zzz[k] + 2.0 * to_zz_zzzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

