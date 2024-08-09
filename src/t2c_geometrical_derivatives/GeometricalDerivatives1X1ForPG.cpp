#include "GeometricalDerivatives1X1ForPG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_pg(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_pg,
                        const size_t idx_op_sf,
                        const size_t idx_op_sh,
                        const size_t idx_op_df,
                        const size_t idx_op_dh,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SF

        auto to_0_xxx = pbuffer.data(idx_op_sf + i * 10 + 0);

        auto to_0_xxy = pbuffer.data(idx_op_sf + i * 10 + 1);

        auto to_0_xxz = pbuffer.data(idx_op_sf + i * 10 + 2);

        auto to_0_xyy = pbuffer.data(idx_op_sf + i * 10 + 3);

        auto to_0_xyz = pbuffer.data(idx_op_sf + i * 10 + 4);

        auto to_0_xzz = pbuffer.data(idx_op_sf + i * 10 + 5);

        auto to_0_yyy = pbuffer.data(idx_op_sf + i * 10 + 6);

        auto to_0_yyz = pbuffer.data(idx_op_sf + i * 10 + 7);

        auto to_0_yzz = pbuffer.data(idx_op_sf + i * 10 + 8);

        auto to_0_zzz = pbuffer.data(idx_op_sf + i * 10 + 9);

        // Set up components of auxiliary buffer : SH

        auto to_0_xxxxx = pbuffer.data(idx_op_sh + i * 21 + 0);

        auto to_0_xxxxy = pbuffer.data(idx_op_sh + i * 21 + 1);

        auto to_0_xxxxz = pbuffer.data(idx_op_sh + i * 21 + 2);

        auto to_0_xxxyy = pbuffer.data(idx_op_sh + i * 21 + 3);

        auto to_0_xxxyz = pbuffer.data(idx_op_sh + i * 21 + 4);

        auto to_0_xxxzz = pbuffer.data(idx_op_sh + i * 21 + 5);

        auto to_0_xxyyy = pbuffer.data(idx_op_sh + i * 21 + 6);

        auto to_0_xxyyz = pbuffer.data(idx_op_sh + i * 21 + 7);

        auto to_0_xxyzz = pbuffer.data(idx_op_sh + i * 21 + 8);

        auto to_0_xxzzz = pbuffer.data(idx_op_sh + i * 21 + 9);

        auto to_0_xyyyy = pbuffer.data(idx_op_sh + i * 21 + 10);

        auto to_0_xyyyz = pbuffer.data(idx_op_sh + i * 21 + 11);

        auto to_0_xyyzz = pbuffer.data(idx_op_sh + i * 21 + 12);

        auto to_0_xyzzz = pbuffer.data(idx_op_sh + i * 21 + 13);

        auto to_0_xzzzz = pbuffer.data(idx_op_sh + i * 21 + 14);

        auto to_0_yyyyy = pbuffer.data(idx_op_sh + i * 21 + 15);

        auto to_0_yyyyz = pbuffer.data(idx_op_sh + i * 21 + 16);

        auto to_0_yyyzz = pbuffer.data(idx_op_sh + i * 21 + 17);

        auto to_0_yyzzz = pbuffer.data(idx_op_sh + i * 21 + 18);

        auto to_0_yzzzz = pbuffer.data(idx_op_sh + i * 21 + 19);

        auto to_0_zzzzz = pbuffer.data(idx_op_sh + i * 21 + 20);

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

        // Set up 0-15 components of targeted buffer : PG

        auto to_x_x_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 0);

        auto to_x_x_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 1);

        auto to_x_x_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 2);

        auto to_x_x_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 3);

        auto to_x_x_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 4);

        auto to_x_x_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 5);

        auto to_x_x_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 6);

        auto to_x_x_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 7);

        auto to_x_x_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 8);

        auto to_x_x_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 9);

        auto to_x_x_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 10);

        auto to_x_x_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 11);

        auto to_x_x_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 12);

        auto to_x_x_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 13);

        auto to_x_x_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxx, to_0_xxxxy, to_0_xxxxz, to_0_xxxyy, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyz, to_0_yzz, to_0_zzz, to_x_x_x_xxxx, to_x_x_x_xxxy, to_x_x_x_xxxz, to_x_x_x_xxyy, to_x_x_x_xxyz, to_x_x_x_xxzz, to_x_x_x_xyyy, to_x_x_x_xyyz, to_x_x_x_xyzz, to_x_x_x_xzzz, to_x_x_x_yyyy, to_x_x_x_yyyz, to_x_x_x_yyzz, to_x_x_x_yzzz, to_x_x_x_zzzz, to_xx_xxx, to_xx_xxxxx, to_xx_xxxxy, to_xx_xxxxz, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_x_xxxx[k] = 4.0 * to_0_xxx[k] - 2.0 * to_0_xxxxx[k] * tke_0 - 8.0 * to_xx_xxx[k] * tbe_0 + 4.0 * to_xx_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_x_xxxy[k] = 3.0 * to_0_xxy[k] - 2.0 * to_0_xxxxy[k] * tke_0 - 6.0 * to_xx_xxy[k] * tbe_0 + 4.0 * to_xx_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_x_xxxz[k] = 3.0 * to_0_xxz[k] - 2.0 * to_0_xxxxz[k] * tke_0 - 6.0 * to_xx_xxz[k] * tbe_0 + 4.0 * to_xx_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_x_xxyy[k] = 2.0 * to_0_xyy[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 4.0 * to_xx_xyy[k] * tbe_0 + 4.0 * to_xx_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_x_xxyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xxxyz[k] * tke_0 - 4.0 * to_xx_xyz[k] * tbe_0 + 4.0 * to_xx_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_x_xxzz[k] = 2.0 * to_0_xzz[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 4.0 * to_xx_xzz[k] * tbe_0 + 4.0 * to_xx_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_x_xyyy[k] = to_0_yyy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 2.0 * to_xx_yyy[k] * tbe_0 + 4.0 * to_xx_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_x_xyyz[k] = to_0_yyz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_xx_yyz[k] * tbe_0 + 4.0 * to_xx_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_x_xyzz[k] = to_0_yzz[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_xx_yzz[k] * tbe_0 + 4.0 * to_xx_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_x_xzzz[k] = to_0_zzz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 2.0 * to_xx_zzz[k] * tbe_0 + 4.0 * to_xx_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_x_yyyy[k] = -2.0 * to_0_xyyyy[k] * tke_0 + 4.0 * to_xx_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_x_yyyz[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_xx_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_x_yyzz[k] = -2.0 * to_0_xyyzz[k] * tke_0 + 4.0 * to_xx_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_x_yzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_xx_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_x_zzzz[k] = -2.0 * to_0_xzzzz[k] * tke_0 + 4.0 * to_xx_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : PG

        auto to_x_x_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 15);

        auto to_x_x_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 16);

        auto to_x_x_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 17);

        auto to_x_x_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 18);

        auto to_x_x_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 19);

        auto to_x_x_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 20);

        auto to_x_x_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 21);

        auto to_x_x_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 22);

        auto to_x_x_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 23);

        auto to_x_x_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 24);

        auto to_x_x_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 25);

        auto to_x_x_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 26);

        auto to_x_x_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 27);

        auto to_x_x_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 28);

        auto to_x_x_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_x_x_y_xxxx, to_x_x_y_xxxy, to_x_x_y_xxxz, to_x_x_y_xxyy, to_x_x_y_xxyz, to_x_x_y_xxzz, to_x_x_y_xyyy, to_x_x_y_xyyz, to_x_x_y_xyzz, to_x_x_y_xzzz, to_x_x_y_yyyy, to_x_x_y_yyyz, to_x_x_y_yyzz, to_x_x_y_yzzz, to_x_x_y_zzzz, to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_y_xxxx[k] = -8.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_y_xxxy[k] = -6.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_y_xxxz[k] = -6.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_y_xxyy[k] = -4.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_y_xxyz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_y_xxzz[k] = -4.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_y_xyyy[k] = -2.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_y_xyyz[k] = -2.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_y_xyzz[k] = -2.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_y_xzzz[k] = -2.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_y_yyyy[k] = 4.0 * to_xy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_y_yyyz[k] = 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_y_yyzz[k] = 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_y_yzzz[k] = 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_y_zzzz[k] = 4.0 * to_xy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : PG

        auto to_x_x_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 30);

        auto to_x_x_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 31);

        auto to_x_x_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 32);

        auto to_x_x_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 33);

        auto to_x_x_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 34);

        auto to_x_x_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 35);

        auto to_x_x_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 36);

        auto to_x_x_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 37);

        auto to_x_x_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 38);

        auto to_x_x_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 39);

        auto to_x_x_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 40);

        auto to_x_x_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 41);

        auto to_x_x_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 42);

        auto to_x_x_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 43);

        auto to_x_x_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 0 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_x_x_z_xxxx, to_x_x_z_xxxy, to_x_x_z_xxxz, to_x_x_z_xxyy, to_x_x_z_xxyz, to_x_x_z_xxzz, to_x_x_z_xyyy, to_x_x_z_xyyz, to_x_x_z_xyzz, to_x_x_z_xzzz, to_x_x_z_yyyy, to_x_x_z_yyyz, to_x_x_z_yyzz, to_x_x_z_yzzz, to_x_x_z_zzzz, to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_z_xxxx[k] = -8.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_z_xxxy[k] = -6.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_z_xxxz[k] = -6.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_z_xxyy[k] = -4.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_z_xxyz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_z_xxzz[k] = -4.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_z_xyyy[k] = -2.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_z_xyyz[k] = -2.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_z_xyzz[k] = -2.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_z_xzzz[k] = -2.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_z_yyyy[k] = 4.0 * to_xz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_z_yyyz[k] = 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_z_yyzz[k] = 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_z_yzzz[k] = 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_z_zzzz[k] = 4.0 * to_xz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : PG

        auto to_x_y_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 0);

        auto to_x_y_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 1);

        auto to_x_y_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 2);

        auto to_x_y_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 3);

        auto to_x_y_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 4);

        auto to_x_y_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 5);

        auto to_x_y_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 6);

        auto to_x_y_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 7);

        auto to_x_y_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 8);

        auto to_x_y_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 9);

        auto to_x_y_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 10);

        auto to_x_y_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 11);

        auto to_x_y_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 12);

        auto to_x_y_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 13);

        auto to_x_y_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxy, to_0_xxxyy, to_0_xxxyz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_yyy, to_0_yyyyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_x_y_x_xxxx, to_x_y_x_xxxy, to_x_y_x_xxxz, to_x_y_x_xxyy, to_x_y_x_xxyz, to_x_y_x_xxzz, to_x_y_x_xyyy, to_x_y_x_xyyz, to_x_y_x_xyzz, to_x_y_x_xzzz, to_x_y_x_yyyy, to_x_y_x_yyyz, to_x_y_x_yyzz, to_x_y_x_yzzz, to_x_y_x_zzzz, to_xx_xxx, to_xx_xxxxy, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_yyy, to_xx_yyyyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_x_xxxx[k] = -2.0 * to_0_xxxxy[k] * tke_0 + 4.0 * to_xx_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_x_xxxy[k] = to_0_xxx[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 2.0 * to_xx_xxx[k] * tbe_0 + 4.0 * to_xx_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_x_xxxz[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_xx_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_x_xxyy[k] = 2.0 * to_0_xxy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 4.0 * to_xx_xxy[k] * tbe_0 + 4.0 * to_xx_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_x_xxyz[k] = to_0_xxz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_xx_xxz[k] * tbe_0 + 4.0 * to_xx_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_x_xxzz[k] = -2.0 * to_0_xxyzz[k] * tke_0 + 4.0 * to_xx_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_x_xyyy[k] = 3.0 * to_0_xyy[k] - 2.0 * to_0_xyyyy[k] * tke_0 - 6.0 * to_xx_xyy[k] * tbe_0 + 4.0 * to_xx_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_x_xyyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyyyz[k] * tke_0 - 4.0 * to_xx_xyz[k] * tbe_0 + 4.0 * to_xx_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_x_xyzz[k] = to_0_xzz[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_xx_xzz[k] * tbe_0 + 4.0 * to_xx_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_x_xzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_xx_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_x_yyyy[k] = 4.0 * to_0_yyy[k] - 2.0 * to_0_yyyyy[k] * tke_0 - 8.0 * to_xx_yyy[k] * tbe_0 + 4.0 * to_xx_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_x_yyyz[k] = 3.0 * to_0_yyz[k] - 2.0 * to_0_yyyyz[k] * tke_0 - 6.0 * to_xx_yyz[k] * tbe_0 + 4.0 * to_xx_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_x_yyzz[k] = 2.0 * to_0_yzz[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 4.0 * to_xx_yzz[k] * tbe_0 + 4.0 * to_xx_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_x_yzzz[k] = to_0_zzz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 2.0 * to_xx_zzz[k] * tbe_0 + 4.0 * to_xx_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_x_zzzz[k] = -2.0 * to_0_yzzzz[k] * tke_0 + 4.0 * to_xx_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : PG

        auto to_x_y_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 15);

        auto to_x_y_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 16);

        auto to_x_y_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 17);

        auto to_x_y_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 18);

        auto to_x_y_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 19);

        auto to_x_y_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 20);

        auto to_x_y_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 21);

        auto to_x_y_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 22);

        auto to_x_y_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 23);

        auto to_x_y_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 24);

        auto to_x_y_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 25);

        auto to_x_y_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 26);

        auto to_x_y_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 27);

        auto to_x_y_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 28);

        auto to_x_y_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_x_y_y_xxxx, to_x_y_y_xxxy, to_x_y_y_xxxz, to_x_y_y_xxyy, to_x_y_y_xxyz, to_x_y_y_xxzz, to_x_y_y_xyyy, to_x_y_y_xyyz, to_x_y_y_xyzz, to_x_y_y_xzzz, to_x_y_y_yyyy, to_x_y_y_yyyz, to_x_y_y_yyzz, to_x_y_y_yzzz, to_x_y_y_zzzz, to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_y_xxxx[k] = 4.0 * to_xy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_y_xxxy[k] = -2.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_y_xxxz[k] = 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_y_xxyy[k] = -4.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_y_xxyz[k] = -2.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_y_xxzz[k] = 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_y_xyyy[k] = -6.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_y_xyyz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_y_xyzz[k] = -2.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_y_xzzz[k] = 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_y_yyyy[k] = -8.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_y_yyyz[k] = -6.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_y_yyzz[k] = -4.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_y_yzzz[k] = -2.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_y_zzzz[k] = 4.0 * to_xy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : PG

        auto to_x_y_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 30);

        auto to_x_y_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 31);

        auto to_x_y_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 32);

        auto to_x_y_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 33);

        auto to_x_y_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 34);

        auto to_x_y_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 35);

        auto to_x_y_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 36);

        auto to_x_y_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 37);

        auto to_x_y_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 38);

        auto to_x_y_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 39);

        auto to_x_y_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 40);

        auto to_x_y_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 41);

        auto to_x_y_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 42);

        auto to_x_y_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 43);

        auto to_x_y_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 1 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_x_y_z_xxxx, to_x_y_z_xxxy, to_x_y_z_xxxz, to_x_y_z_xxyy, to_x_y_z_xxyz, to_x_y_z_xxzz, to_x_y_z_xyyy, to_x_y_z_xyyz, to_x_y_z_xyzz, to_x_y_z_xzzz, to_x_y_z_yyyy, to_x_y_z_yyyz, to_x_y_z_yyzz, to_x_y_z_yzzz, to_x_y_z_zzzz, to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_z_xxxx[k] = 4.0 * to_xz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_z_xxxy[k] = -2.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_z_xxxz[k] = 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_z_xxyy[k] = -4.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_z_xxyz[k] = -2.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_z_xxzz[k] = 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_z_xyyy[k] = -6.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_z_xyyz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_z_xyzz[k] = -2.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_z_xzzz[k] = 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_z_yyyy[k] = -8.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_z_yyyz[k] = -6.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_z_yyzz[k] = -4.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_z_yzzz[k] = -2.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_z_zzzz[k] = 4.0 * to_xz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : PG

        auto to_x_z_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 0);

        auto to_x_z_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 1);

        auto to_x_z_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 2);

        auto to_x_z_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 3);

        auto to_x_z_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 4);

        auto to_x_z_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 5);

        auto to_x_z_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 6);

        auto to_x_z_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 7);

        auto to_x_z_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 8);

        auto to_x_z_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 9);

        auto to_x_z_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 10);

        auto to_x_z_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 11);

        auto to_x_z_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 12);

        auto to_x_z_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 13);

        auto to_x_z_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxz, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_0_zzzzz, to_x_z_x_xxxx, to_x_z_x_xxxy, to_x_z_x_xxxz, to_x_z_x_xxyy, to_x_z_x_xxyz, to_x_z_x_xxzz, to_x_z_x_xyyy, to_x_z_x_xyyz, to_x_z_x_xyzz, to_x_z_x_xzzz, to_x_z_x_yyyy, to_x_z_x_yyyz, to_x_z_x_yyzz, to_x_z_x_yzzz, to_x_z_x_zzzz, to_xx_xxx, to_xx_xxxxz, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_x_xxxx[k] = -2.0 * to_0_xxxxz[k] * tke_0 + 4.0 * to_xx_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_x_xxxy[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_xx_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_x_xxxz[k] = to_0_xxx[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 2.0 * to_xx_xxx[k] * tbe_0 + 4.0 * to_xx_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_x_xxyy[k] = -2.0 * to_0_xxyyz[k] * tke_0 + 4.0 * to_xx_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_x_xxyz[k] = to_0_xxy[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_xx_xxy[k] * tbe_0 + 4.0 * to_xx_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_x_xxzz[k] = 2.0 * to_0_xxz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 4.0 * to_xx_xxz[k] * tbe_0 + 4.0 * to_xx_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_x_xyyy[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_xx_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_x_xyyz[k] = to_0_xyy[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_xx_xyy[k] * tbe_0 + 4.0 * to_xx_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_x_xyzz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyzzz[k] * tke_0 - 4.0 * to_xx_xyz[k] * tbe_0 + 4.0 * to_xx_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_x_xzzz[k] = 3.0 * to_0_xzz[k] - 2.0 * to_0_xzzzz[k] * tke_0 - 6.0 * to_xx_xzz[k] * tbe_0 + 4.0 * to_xx_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_x_yyyy[k] = -2.0 * to_0_yyyyz[k] * tke_0 + 4.0 * to_xx_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_x_yyyz[k] = to_0_yyy[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 2.0 * to_xx_yyy[k] * tbe_0 + 4.0 * to_xx_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_x_yyzz[k] = 2.0 * to_0_yyz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 4.0 * to_xx_yyz[k] * tbe_0 + 4.0 * to_xx_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_x_yzzz[k] = 3.0 * to_0_yzz[k] - 2.0 * to_0_yzzzz[k] * tke_0 - 6.0 * to_xx_yzz[k] * tbe_0 + 4.0 * to_xx_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_x_zzzz[k] = 4.0 * to_0_zzz[k] - 2.0 * to_0_zzzzz[k] * tke_0 - 8.0 * to_xx_zzz[k] * tbe_0 + 4.0 * to_xx_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : PG

        auto to_x_z_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 15);

        auto to_x_z_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 16);

        auto to_x_z_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 17);

        auto to_x_z_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 18);

        auto to_x_z_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 19);

        auto to_x_z_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 20);

        auto to_x_z_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 21);

        auto to_x_z_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 22);

        auto to_x_z_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 23);

        auto to_x_z_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 24);

        auto to_x_z_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 25);

        auto to_x_z_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 26);

        auto to_x_z_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 27);

        auto to_x_z_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 28);

        auto to_x_z_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_x_z_y_xxxx, to_x_z_y_xxxy, to_x_z_y_xxxz, to_x_z_y_xxyy, to_x_z_y_xxyz, to_x_z_y_xxzz, to_x_z_y_xyyy, to_x_z_y_xyyz, to_x_z_y_xyzz, to_x_z_y_xzzz, to_x_z_y_yyyy, to_x_z_y_yyyz, to_x_z_y_yyzz, to_x_z_y_yzzz, to_x_z_y_zzzz, to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_y_xxxx[k] = 4.0 * to_xy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_y_xxxy[k] = 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_y_xxxz[k] = -2.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_y_xxyy[k] = 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_y_xxyz[k] = -2.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_y_xxzz[k] = -4.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_y_xyyy[k] = 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_y_xyyz[k] = -2.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_y_xyzz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_y_xzzz[k] = -6.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_y_yyyy[k] = 4.0 * to_xy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_y_yyyz[k] = -2.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_y_yyzz[k] = -4.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_y_yzzz[k] = -6.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_y_zzzz[k] = -8.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : PG

        auto to_x_z_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 30);

        auto to_x_z_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 31);

        auto to_x_z_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 32);

        auto to_x_z_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 33);

        auto to_x_z_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 34);

        auto to_x_z_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 35);

        auto to_x_z_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 36);

        auto to_x_z_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 37);

        auto to_x_z_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 38);

        auto to_x_z_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 39);

        auto to_x_z_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 40);

        auto to_x_z_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 41);

        auto to_x_z_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 42);

        auto to_x_z_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 43);

        auto to_x_z_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 2 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_x_z_z_xxxx, to_x_z_z_xxxy, to_x_z_z_xxxz, to_x_z_z_xxyy, to_x_z_z_xxyz, to_x_z_z_xxzz, to_x_z_z_xyyy, to_x_z_z_xyyz, to_x_z_z_xyzz, to_x_z_z_xzzz, to_x_z_z_yyyy, to_x_z_z_yyyz, to_x_z_z_yyzz, to_x_z_z_yzzz, to_x_z_z_zzzz, to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_z_xxxx[k] = 4.0 * to_xz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_z_xxxy[k] = 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_z_xxxz[k] = -2.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_z_xxyy[k] = 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_z_xxyz[k] = -2.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_z_xxzz[k] = -4.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_z_xyyy[k] = 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_z_xyyz[k] = -2.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_z_xyzz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_z_xzzz[k] = -6.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_z_yyyy[k] = 4.0 * to_xz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_z_yyyz[k] = -2.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_z_yyzz[k] = -4.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_z_yzzz[k] = -6.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_z_zzzz[k] = -8.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 135-150 components of targeted buffer : PG

        auto to_y_x_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 0);

        auto to_y_x_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 1);

        auto to_y_x_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 2);

        auto to_y_x_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 3);

        auto to_y_x_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 4);

        auto to_y_x_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 5);

        auto to_y_x_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 6);

        auto to_y_x_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 7);

        auto to_y_x_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 8);

        auto to_y_x_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 9);

        auto to_y_x_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 10);

        auto to_y_x_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 11);

        auto to_y_x_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 12);

        auto to_y_x_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 13);

        auto to_y_x_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, to_y_x_x_xxxx, to_y_x_x_xxxy, to_y_x_x_xxxz, to_y_x_x_xxyy, to_y_x_x_xxyz, to_y_x_x_xxzz, to_y_x_x_xyyy, to_y_x_x_xyyz, to_y_x_x_xyzz, to_y_x_x_xzzz, to_y_x_x_yyyy, to_y_x_x_yyyz, to_y_x_x_yyzz, to_y_x_x_yzzz, to_y_x_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_x_xxxx[k] = -8.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_x_xxxy[k] = -6.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_x_xxxz[k] = -6.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_x_xxyy[k] = -4.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_x_xxyz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_x_xxzz[k] = -4.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_x_xyyy[k] = -2.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_x_xyyz[k] = -2.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_x_xyzz[k] = -2.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_x_xzzz[k] = -2.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_x_yyyy[k] = 4.0 * to_xy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_x_yyyz[k] = 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_x_yyzz[k] = 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_x_yzzz[k] = 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_x_zzzz[k] = 4.0 * to_xy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-165 components of targeted buffer : PG

        auto to_y_x_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 15);

        auto to_y_x_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 16);

        auto to_y_x_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 17);

        auto to_y_x_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 18);

        auto to_y_x_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 19);

        auto to_y_x_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 20);

        auto to_y_x_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 21);

        auto to_y_x_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 22);

        auto to_y_x_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 23);

        auto to_y_x_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 24);

        auto to_y_x_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 25);

        auto to_y_x_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 26);

        auto to_y_x_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 27);

        auto to_y_x_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 28);

        auto to_y_x_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxx, to_0_xxxxy, to_0_xxxxz, to_0_xxxyy, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyz, to_0_yzz, to_0_zzz, to_y_x_y_xxxx, to_y_x_y_xxxy, to_y_x_y_xxxz, to_y_x_y_xxyy, to_y_x_y_xxyz, to_y_x_y_xxzz, to_y_x_y_xyyy, to_y_x_y_xyyz, to_y_x_y_xyzz, to_y_x_y_xzzz, to_y_x_y_yyyy, to_y_x_y_yyyz, to_y_x_y_yyzz, to_y_x_y_yzzz, to_y_x_y_zzzz, to_yy_xxx, to_yy_xxxxx, to_yy_xxxxy, to_yy_xxxxz, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_y_xxxx[k] = 4.0 * to_0_xxx[k] - 2.0 * to_0_xxxxx[k] * tke_0 - 8.0 * to_yy_xxx[k] * tbe_0 + 4.0 * to_yy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_y_xxxy[k] = 3.0 * to_0_xxy[k] - 2.0 * to_0_xxxxy[k] * tke_0 - 6.0 * to_yy_xxy[k] * tbe_0 + 4.0 * to_yy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_y_xxxz[k] = 3.0 * to_0_xxz[k] - 2.0 * to_0_xxxxz[k] * tke_0 - 6.0 * to_yy_xxz[k] * tbe_0 + 4.0 * to_yy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_y_xxyy[k] = 2.0 * to_0_xyy[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 4.0 * to_yy_xyy[k] * tbe_0 + 4.0 * to_yy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_y_xxyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xxxyz[k] * tke_0 - 4.0 * to_yy_xyz[k] * tbe_0 + 4.0 * to_yy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_y_xxzz[k] = 2.0 * to_0_xzz[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 4.0 * to_yy_xzz[k] * tbe_0 + 4.0 * to_yy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_y_xyyy[k] = to_0_yyy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 2.0 * to_yy_yyy[k] * tbe_0 + 4.0 * to_yy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_y_xyyz[k] = to_0_yyz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_yy_yyz[k] * tbe_0 + 4.0 * to_yy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_y_xyzz[k] = to_0_yzz[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_yy_yzz[k] * tbe_0 + 4.0 * to_yy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_y_xzzz[k] = to_0_zzz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 2.0 * to_yy_zzz[k] * tbe_0 + 4.0 * to_yy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_y_yyyy[k] = -2.0 * to_0_xyyyy[k] * tke_0 + 4.0 * to_yy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_y_yyyz[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_yy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_y_yyzz[k] = -2.0 * to_0_xyyzz[k] * tke_0 + 4.0 * to_yy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_y_yzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_yy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_y_zzzz[k] = -2.0 * to_0_xzzzz[k] * tke_0 + 4.0 * to_yy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 165-180 components of targeted buffer : PG

        auto to_y_x_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 30);

        auto to_y_x_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 31);

        auto to_y_x_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 32);

        auto to_y_x_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 33);

        auto to_y_x_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 34);

        auto to_y_x_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 35);

        auto to_y_x_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 36);

        auto to_y_x_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 37);

        auto to_y_x_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 38);

        auto to_y_x_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 39);

        auto to_y_x_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 40);

        auto to_y_x_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 41);

        auto to_y_x_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 42);

        auto to_y_x_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 43);

        auto to_y_x_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 3 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_y_x_z_xxxx, to_y_x_z_xxxy, to_y_x_z_xxxz, to_y_x_z_xxyy, to_y_x_z_xxyz, to_y_x_z_xxzz, to_y_x_z_xyyy, to_y_x_z_xyyz, to_y_x_z_xyzz, to_y_x_z_xzzz, to_y_x_z_yyyy, to_y_x_z_yyyz, to_y_x_z_yyzz, to_y_x_z_yzzz, to_y_x_z_zzzz, to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_z_xxxx[k] = -8.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_z_xxxy[k] = -6.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_z_xxxz[k] = -6.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_z_xxyy[k] = -4.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_z_xxyz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_z_xxzz[k] = -4.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_z_xyyy[k] = -2.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_z_xyyz[k] = -2.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_z_xyzz[k] = -2.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_z_xzzz[k] = -2.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_z_yyyy[k] = 4.0 * to_yz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_z_yyyz[k] = 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_z_yyzz[k] = 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_z_yzzz[k] = 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_z_zzzz[k] = 4.0 * to_yz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-195 components of targeted buffer : PG

        auto to_y_y_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 0);

        auto to_y_y_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 1);

        auto to_y_y_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 2);

        auto to_y_y_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 3);

        auto to_y_y_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 4);

        auto to_y_y_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 5);

        auto to_y_y_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 6);

        auto to_y_y_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 7);

        auto to_y_y_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 8);

        auto to_y_y_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 9);

        auto to_y_y_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 10);

        auto to_y_y_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 11);

        auto to_y_y_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 12);

        auto to_y_y_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 13);

        auto to_y_y_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_y_y_x_xxxx, to_y_y_x_xxxy, to_y_y_x_xxxz, to_y_y_x_xxyy, to_y_y_x_xxyz, to_y_y_x_xxzz, to_y_y_x_xyyy, to_y_y_x_xyyz, to_y_y_x_xyzz, to_y_y_x_xzzz, to_y_y_x_yyyy, to_y_y_x_yyyz, to_y_y_x_yyzz, to_y_y_x_yzzz, to_y_y_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_x_xxxx[k] = 4.0 * to_xy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_x_xxxy[k] = -2.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_x_xxxz[k] = 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_x_xxyy[k] = -4.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_x_xxyz[k] = -2.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_x_xxzz[k] = 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_x_xyyy[k] = -6.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_x_xyyz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_x_xyzz[k] = -2.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_x_xzzz[k] = 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_x_yyyy[k] = -8.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_x_yyyz[k] = -6.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_x_yyzz[k] = -4.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_x_yzzz[k] = -2.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_x_zzzz[k] = 4.0 * to_xy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 195-210 components of targeted buffer : PG

        auto to_y_y_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 15);

        auto to_y_y_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 16);

        auto to_y_y_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 17);

        auto to_y_y_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 18);

        auto to_y_y_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 19);

        auto to_y_y_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 20);

        auto to_y_y_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 21);

        auto to_y_y_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 22);

        auto to_y_y_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 23);

        auto to_y_y_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 24);

        auto to_y_y_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 25);

        auto to_y_y_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 26);

        auto to_y_y_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 27);

        auto to_y_y_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 28);

        auto to_y_y_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxy, to_0_xxxyy, to_0_xxxyz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_yyy, to_0_yyyyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_y_y_y_xxxx, to_y_y_y_xxxy, to_y_y_y_xxxz, to_y_y_y_xxyy, to_y_y_y_xxyz, to_y_y_y_xxzz, to_y_y_y_xyyy, to_y_y_y_xyyz, to_y_y_y_xyzz, to_y_y_y_xzzz, to_y_y_y_yyyy, to_y_y_y_yyyz, to_y_y_y_yyzz, to_y_y_y_yzzz, to_y_y_y_zzzz, to_yy_xxx, to_yy_xxxxy, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_yyy, to_yy_yyyyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_y_xxxx[k] = -2.0 * to_0_xxxxy[k] * tke_0 + 4.0 * to_yy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_y_xxxy[k] = to_0_xxx[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 2.0 * to_yy_xxx[k] * tbe_0 + 4.0 * to_yy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_y_xxxz[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_yy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_y_xxyy[k] = 2.0 * to_0_xxy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 4.0 * to_yy_xxy[k] * tbe_0 + 4.0 * to_yy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_y_xxyz[k] = to_0_xxz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_yy_xxz[k] * tbe_0 + 4.0 * to_yy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_y_xxzz[k] = -2.0 * to_0_xxyzz[k] * tke_0 + 4.0 * to_yy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_y_xyyy[k] = 3.0 * to_0_xyy[k] - 2.0 * to_0_xyyyy[k] * tke_0 - 6.0 * to_yy_xyy[k] * tbe_0 + 4.0 * to_yy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_y_xyyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyyyz[k] * tke_0 - 4.0 * to_yy_xyz[k] * tbe_0 + 4.0 * to_yy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_y_xyzz[k] = to_0_xzz[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_yy_xzz[k] * tbe_0 + 4.0 * to_yy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_y_xzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_yy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_y_yyyy[k] = 4.0 * to_0_yyy[k] - 2.0 * to_0_yyyyy[k] * tke_0 - 8.0 * to_yy_yyy[k] * tbe_0 + 4.0 * to_yy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_y_yyyz[k] = 3.0 * to_0_yyz[k] - 2.0 * to_0_yyyyz[k] * tke_0 - 6.0 * to_yy_yyz[k] * tbe_0 + 4.0 * to_yy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_y_yyzz[k] = 2.0 * to_0_yzz[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 4.0 * to_yy_yzz[k] * tbe_0 + 4.0 * to_yy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_y_yzzz[k] = to_0_zzz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 2.0 * to_yy_zzz[k] * tbe_0 + 4.0 * to_yy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_y_zzzz[k] = -2.0 * to_0_yzzzz[k] * tke_0 + 4.0 * to_yy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-225 components of targeted buffer : PG

        auto to_y_y_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 30);

        auto to_y_y_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 31);

        auto to_y_y_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 32);

        auto to_y_y_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 33);

        auto to_y_y_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 34);

        auto to_y_y_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 35);

        auto to_y_y_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 36);

        auto to_y_y_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 37);

        auto to_y_y_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 38);

        auto to_y_y_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 39);

        auto to_y_y_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 40);

        auto to_y_y_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 41);

        auto to_y_y_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 42);

        auto to_y_y_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 43);

        auto to_y_y_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 4 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_y_y_z_xxxx, to_y_y_z_xxxy, to_y_y_z_xxxz, to_y_y_z_xxyy, to_y_y_z_xxyz, to_y_y_z_xxzz, to_y_y_z_xyyy, to_y_y_z_xyyz, to_y_y_z_xyzz, to_y_y_z_xzzz, to_y_y_z_yyyy, to_y_y_z_yyyz, to_y_y_z_yyzz, to_y_y_z_yzzz, to_y_y_z_zzzz, to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_z_xxxx[k] = 4.0 * to_yz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_z_xxxy[k] = -2.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_z_xxxz[k] = 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_z_xxyy[k] = -4.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_z_xxyz[k] = -2.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_z_xxzz[k] = 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_z_xyyy[k] = -6.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_z_xyyz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_z_xyzz[k] = -2.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_z_xzzz[k] = 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_z_yyyy[k] = -8.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_z_yyyz[k] = -6.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_z_yyzz[k] = -4.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_z_yzzz[k] = -2.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_z_zzzz[k] = 4.0 * to_yz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 225-240 components of targeted buffer : PG

        auto to_y_z_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 0);

        auto to_y_z_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 1);

        auto to_y_z_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 2);

        auto to_y_z_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 3);

        auto to_y_z_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 4);

        auto to_y_z_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 5);

        auto to_y_z_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 6);

        auto to_y_z_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 7);

        auto to_y_z_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 8);

        auto to_y_z_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 9);

        auto to_y_z_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 10);

        auto to_y_z_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 11);

        auto to_y_z_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 12);

        auto to_y_z_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 13);

        auto to_y_z_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, to_y_z_x_xxxx, to_y_z_x_xxxy, to_y_z_x_xxxz, to_y_z_x_xxyy, to_y_z_x_xxyz, to_y_z_x_xxzz, to_y_z_x_xyyy, to_y_z_x_xyyz, to_y_z_x_xyzz, to_y_z_x_xzzz, to_y_z_x_yyyy, to_y_z_x_yyyz, to_y_z_x_yyzz, to_y_z_x_yzzz, to_y_z_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_x_xxxx[k] = 4.0 * to_xy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_x_xxxy[k] = 4.0 * to_xy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_x_xxxz[k] = -2.0 * to_xy_xxx[k] * tbe_0 + 4.0 * to_xy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_x_xxyy[k] = 4.0 * to_xy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_x_xxyz[k] = -2.0 * to_xy_xxy[k] * tbe_0 + 4.0 * to_xy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_x_xxzz[k] = -4.0 * to_xy_xxz[k] * tbe_0 + 4.0 * to_xy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_x_xyyy[k] = 4.0 * to_xy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_x_xyyz[k] = -2.0 * to_xy_xyy[k] * tbe_0 + 4.0 * to_xy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_x_xyzz[k] = -4.0 * to_xy_xyz[k] * tbe_0 + 4.0 * to_xy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_x_xzzz[k] = -6.0 * to_xy_xzz[k] * tbe_0 + 4.0 * to_xy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_x_yyyy[k] = 4.0 * to_xy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_x_yyyz[k] = -2.0 * to_xy_yyy[k] * tbe_0 + 4.0 * to_xy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_x_yyzz[k] = -4.0 * to_xy_yyz[k] * tbe_0 + 4.0 * to_xy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_x_yzzz[k] = -6.0 * to_xy_yzz[k] * tbe_0 + 4.0 * to_xy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_x_zzzz[k] = -8.0 * to_xy_zzz[k] * tbe_0 + 4.0 * to_xy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-255 components of targeted buffer : PG

        auto to_y_z_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 15);

        auto to_y_z_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 16);

        auto to_y_z_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 17);

        auto to_y_z_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 18);

        auto to_y_z_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 19);

        auto to_y_z_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 20);

        auto to_y_z_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 21);

        auto to_y_z_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 22);

        auto to_y_z_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 23);

        auto to_y_z_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 24);

        auto to_y_z_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 25);

        auto to_y_z_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 26);

        auto to_y_z_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 27);

        auto to_y_z_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 28);

        auto to_y_z_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxz, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_0_zzzzz, to_y_z_y_xxxx, to_y_z_y_xxxy, to_y_z_y_xxxz, to_y_z_y_xxyy, to_y_z_y_xxyz, to_y_z_y_xxzz, to_y_z_y_xyyy, to_y_z_y_xyyz, to_y_z_y_xyzz, to_y_z_y_xzzz, to_y_z_y_yyyy, to_y_z_y_yyyz, to_y_z_y_yyzz, to_y_z_y_yzzz, to_y_z_y_zzzz, to_yy_xxx, to_yy_xxxxz, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_y_xxxx[k] = -2.0 * to_0_xxxxz[k] * tke_0 + 4.0 * to_yy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_y_xxxy[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_yy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_y_xxxz[k] = to_0_xxx[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 2.0 * to_yy_xxx[k] * tbe_0 + 4.0 * to_yy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_y_xxyy[k] = -2.0 * to_0_xxyyz[k] * tke_0 + 4.0 * to_yy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_y_xxyz[k] = to_0_xxy[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_yy_xxy[k] * tbe_0 + 4.0 * to_yy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_y_xxzz[k] = 2.0 * to_0_xxz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 4.0 * to_yy_xxz[k] * tbe_0 + 4.0 * to_yy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_y_xyyy[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_yy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_y_xyyz[k] = to_0_xyy[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_yy_xyy[k] * tbe_0 + 4.0 * to_yy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_y_xyzz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyzzz[k] * tke_0 - 4.0 * to_yy_xyz[k] * tbe_0 + 4.0 * to_yy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_y_xzzz[k] = 3.0 * to_0_xzz[k] - 2.0 * to_0_xzzzz[k] * tke_0 - 6.0 * to_yy_xzz[k] * tbe_0 + 4.0 * to_yy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_y_yyyy[k] = -2.0 * to_0_yyyyz[k] * tke_0 + 4.0 * to_yy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_y_yyyz[k] = to_0_yyy[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 2.0 * to_yy_yyy[k] * tbe_0 + 4.0 * to_yy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_y_yyzz[k] = 2.0 * to_0_yyz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 4.0 * to_yy_yyz[k] * tbe_0 + 4.0 * to_yy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_y_yzzz[k] = 3.0 * to_0_yzz[k] - 2.0 * to_0_yzzzz[k] * tke_0 - 6.0 * to_yy_yzz[k] * tbe_0 + 4.0 * to_yy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_y_zzzz[k] = 4.0 * to_0_zzz[k] - 2.0 * to_0_zzzzz[k] * tke_0 - 8.0 * to_yy_zzz[k] * tbe_0 + 4.0 * to_yy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 255-270 components of targeted buffer : PG

        auto to_y_z_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 30);

        auto to_y_z_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 31);

        auto to_y_z_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 32);

        auto to_y_z_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 33);

        auto to_y_z_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 34);

        auto to_y_z_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 35);

        auto to_y_z_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 36);

        auto to_y_z_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 37);

        auto to_y_z_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 38);

        auto to_y_z_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 39);

        auto to_y_z_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 40);

        auto to_y_z_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 41);

        auto to_y_z_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 42);

        auto to_y_z_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 43);

        auto to_y_z_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 5 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_y_z_z_xxxx, to_y_z_z_xxxy, to_y_z_z_xxxz, to_y_z_z_xxyy, to_y_z_z_xxyz, to_y_z_z_xxzz, to_y_z_z_xyyy, to_y_z_z_xyyz, to_y_z_z_xyzz, to_y_z_z_xzzz, to_y_z_z_yyyy, to_y_z_z_yyyz, to_y_z_z_yyzz, to_y_z_z_yzzz, to_y_z_z_zzzz, to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_z_xxxx[k] = 4.0 * to_yz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_z_xxxy[k] = 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_z_xxxz[k] = -2.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_z_xxyy[k] = 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_z_xxyz[k] = -2.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_z_xxzz[k] = -4.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_z_xyyy[k] = 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_z_xyyz[k] = -2.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_z_xyzz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_z_xzzz[k] = -6.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_z_yyyy[k] = 4.0 * to_yz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_z_yyyz[k] = -2.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_z_yyzz[k] = -4.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_z_yzzz[k] = -6.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_z_zzzz[k] = -8.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-285 components of targeted buffer : PG

        auto to_z_x_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 0);

        auto to_z_x_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 1);

        auto to_z_x_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 2);

        auto to_z_x_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 3);

        auto to_z_x_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 4);

        auto to_z_x_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 5);

        auto to_z_x_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 6);

        auto to_z_x_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 7);

        auto to_z_x_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 8);

        auto to_z_x_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 9);

        auto to_z_x_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 10);

        auto to_z_x_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 11);

        auto to_z_x_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 12);

        auto to_z_x_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 13);

        auto to_z_x_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, to_z_x_x_xxxx, to_z_x_x_xxxy, to_z_x_x_xxxz, to_z_x_x_xxyy, to_z_x_x_xxyz, to_z_x_x_xxzz, to_z_x_x_xyyy, to_z_x_x_xyyz, to_z_x_x_xyzz, to_z_x_x_xzzz, to_z_x_x_yyyy, to_z_x_x_yyyz, to_z_x_x_yyzz, to_z_x_x_yzzz, to_z_x_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_x_xxxx[k] = -8.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_x_xxxy[k] = -6.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_x_xxxz[k] = -6.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_x_xxyy[k] = -4.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_x_xxyz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_x_xxzz[k] = -4.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_x_xyyy[k] = -2.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_x_xyyz[k] = -2.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_x_xyzz[k] = -2.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_x_xzzz[k] = -2.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_x_yyyy[k] = 4.0 * to_xz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_x_yyyz[k] = 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_x_yyzz[k] = 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_x_yzzz[k] = 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_x_zzzz[k] = 4.0 * to_xz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 285-300 components of targeted buffer : PG

        auto to_z_x_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 15);

        auto to_z_x_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 16);

        auto to_z_x_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 17);

        auto to_z_x_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 18);

        auto to_z_x_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 19);

        auto to_z_x_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 20);

        auto to_z_x_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 21);

        auto to_z_x_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 22);

        auto to_z_x_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 23);

        auto to_z_x_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 24);

        auto to_z_x_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 25);

        auto to_z_x_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 26);

        auto to_z_x_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 27);

        auto to_z_x_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 28);

        auto to_z_x_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, to_z_x_y_xxxx, to_z_x_y_xxxy, to_z_x_y_xxxz, to_z_x_y_xxyy, to_z_x_y_xxyz, to_z_x_y_xxzz, to_z_x_y_xyyy, to_z_x_y_xyyz, to_z_x_y_xyzz, to_z_x_y_xzzz, to_z_x_y_yyyy, to_z_x_y_yyyz, to_z_x_y_yyzz, to_z_x_y_yzzz, to_z_x_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_y_xxxx[k] = -8.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_y_xxxy[k] = -6.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_y_xxxz[k] = -6.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_y_xxyy[k] = -4.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_y_xxyz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_y_xxzz[k] = -4.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_y_xyyy[k] = -2.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_y_xyyz[k] = -2.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_y_xyzz[k] = -2.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_y_xzzz[k] = -2.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_y_yyyy[k] = 4.0 * to_yz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_y_yyyz[k] = 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_y_yyzz[k] = 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_y_yzzz[k] = 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_y_zzzz[k] = 4.0 * to_yz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-315 components of targeted buffer : PG

        auto to_z_x_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 30);

        auto to_z_x_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 31);

        auto to_z_x_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 32);

        auto to_z_x_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 33);

        auto to_z_x_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 34);

        auto to_z_x_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 35);

        auto to_z_x_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 36);

        auto to_z_x_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 37);

        auto to_z_x_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 38);

        auto to_z_x_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 39);

        auto to_z_x_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 40);

        auto to_z_x_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 41);

        auto to_z_x_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 42);

        auto to_z_x_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 43);

        auto to_z_x_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 6 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxx, to_0_xxxxy, to_0_xxxxz, to_0_xxxyy, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyz, to_0_yzz, to_0_zzz, to_z_x_z_xxxx, to_z_x_z_xxxy, to_z_x_z_xxxz, to_z_x_z_xxyy, to_z_x_z_xxyz, to_z_x_z_xxzz, to_z_x_z_xyyy, to_z_x_z_xyyz, to_z_x_z_xyzz, to_z_x_z_xzzz, to_z_x_z_yyyy, to_z_x_z_yyyz, to_z_x_z_yyzz, to_z_x_z_yzzz, to_z_x_z_zzzz, to_zz_xxx, to_zz_xxxxx, to_zz_xxxxy, to_zz_xxxxz, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_z_xxxx[k] = 4.0 * to_0_xxx[k] - 2.0 * to_0_xxxxx[k] * tke_0 - 8.0 * to_zz_xxx[k] * tbe_0 + 4.0 * to_zz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_z_xxxy[k] = 3.0 * to_0_xxy[k] - 2.0 * to_0_xxxxy[k] * tke_0 - 6.0 * to_zz_xxy[k] * tbe_0 + 4.0 * to_zz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_z_xxxz[k] = 3.0 * to_0_xxz[k] - 2.0 * to_0_xxxxz[k] * tke_0 - 6.0 * to_zz_xxz[k] * tbe_0 + 4.0 * to_zz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_z_xxyy[k] = 2.0 * to_0_xyy[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 4.0 * to_zz_xyy[k] * tbe_0 + 4.0 * to_zz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_z_xxyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xxxyz[k] * tke_0 - 4.0 * to_zz_xyz[k] * tbe_0 + 4.0 * to_zz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_z_xxzz[k] = 2.0 * to_0_xzz[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 4.0 * to_zz_xzz[k] * tbe_0 + 4.0 * to_zz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_z_xyyy[k] = to_0_yyy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 2.0 * to_zz_yyy[k] * tbe_0 + 4.0 * to_zz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_z_xyyz[k] = to_0_yyz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_zz_yyz[k] * tbe_0 + 4.0 * to_zz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_z_xyzz[k] = to_0_yzz[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_zz_yzz[k] * tbe_0 + 4.0 * to_zz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_z_xzzz[k] = to_0_zzz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 2.0 * to_zz_zzz[k] * tbe_0 + 4.0 * to_zz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_z_yyyy[k] = -2.0 * to_0_xyyyy[k] * tke_0 + 4.0 * to_zz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_z_yyyz[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_zz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_z_yyzz[k] = -2.0 * to_0_xyyzz[k] * tke_0 + 4.0 * to_zz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_z_yzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_zz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_z_zzzz[k] = -2.0 * to_0_xzzzz[k] * tke_0 + 4.0 * to_zz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 315-330 components of targeted buffer : PG

        auto to_z_y_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 0);

        auto to_z_y_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 1);

        auto to_z_y_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 2);

        auto to_z_y_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 3);

        auto to_z_y_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 4);

        auto to_z_y_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 5);

        auto to_z_y_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 6);

        auto to_z_y_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 7);

        auto to_z_y_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 8);

        auto to_z_y_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 9);

        auto to_z_y_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 10);

        auto to_z_y_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 11);

        auto to_z_y_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 12);

        auto to_z_y_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 13);

        auto to_z_y_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_z_y_x_xxxx, to_z_y_x_xxxy, to_z_y_x_xxxz, to_z_y_x_xxyy, to_z_y_x_xxyz, to_z_y_x_xxzz, to_z_y_x_xyyy, to_z_y_x_xyyz, to_z_y_x_xyzz, to_z_y_x_xzzz, to_z_y_x_yyyy, to_z_y_x_yyyz, to_z_y_x_yyzz, to_z_y_x_yzzz, to_z_y_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_x_xxxx[k] = 4.0 * to_xz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_x_xxxy[k] = -2.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_x_xxxz[k] = 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_x_xxyy[k] = -4.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_x_xxyz[k] = -2.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_x_xxzz[k] = 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_x_xyyy[k] = -6.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_x_xyyz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_x_xyzz[k] = -2.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_x_xzzz[k] = 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_x_yyyy[k] = -8.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_x_yyyz[k] = -6.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_x_yyzz[k] = -4.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_x_yzzz[k] = -2.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_x_zzzz[k] = 4.0 * to_xz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-345 components of targeted buffer : PG

        auto to_z_y_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 15);

        auto to_z_y_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 16);

        auto to_z_y_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 17);

        auto to_z_y_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 18);

        auto to_z_y_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 19);

        auto to_z_y_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 20);

        auto to_z_y_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 21);

        auto to_z_y_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 22);

        auto to_z_y_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 23);

        auto to_z_y_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 24);

        auto to_z_y_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 25);

        auto to_z_y_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 26);

        auto to_z_y_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 27);

        auto to_z_y_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 28);

        auto to_z_y_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_z_y_y_xxxx, to_z_y_y_xxxy, to_z_y_y_xxxz, to_z_y_y_xxyy, to_z_y_y_xxyz, to_z_y_y_xxzz, to_z_y_y_xyyy, to_z_y_y_xyyz, to_z_y_y_xyzz, to_z_y_y_xzzz, to_z_y_y_yyyy, to_z_y_y_yyyz, to_z_y_y_yyzz, to_z_y_y_yzzz, to_z_y_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_y_xxxx[k] = 4.0 * to_yz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_y_xxxy[k] = -2.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_y_xxxz[k] = 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_y_xxyy[k] = -4.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_y_xxyz[k] = -2.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_y_xxzz[k] = 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_y_xyyy[k] = -6.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_y_xyyz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_y_xyzz[k] = -2.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_y_xzzz[k] = 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_y_yyyy[k] = -8.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_y_yyyz[k] = -6.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_y_yyzz[k] = -4.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_y_yzzz[k] = -2.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_y_zzzz[k] = 4.0 * to_yz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 345-360 components of targeted buffer : PG

        auto to_z_y_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 30);

        auto to_z_y_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 31);

        auto to_z_y_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 32);

        auto to_z_y_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 33);

        auto to_z_y_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 34);

        auto to_z_y_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 35);

        auto to_z_y_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 36);

        auto to_z_y_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 37);

        auto to_z_y_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 38);

        auto to_z_y_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 39);

        auto to_z_y_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 40);

        auto to_z_y_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 41);

        auto to_z_y_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 42);

        auto to_z_y_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 43);

        auto to_z_y_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 7 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxy, to_0_xxxyy, to_0_xxxyz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_yyy, to_0_yyyyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_z_y_z_xxxx, to_z_y_z_xxxy, to_z_y_z_xxxz, to_z_y_z_xxyy, to_z_y_z_xxyz, to_z_y_z_xxzz, to_z_y_z_xyyy, to_z_y_z_xyyz, to_z_y_z_xyzz, to_z_y_z_xzzz, to_z_y_z_yyyy, to_z_y_z_yyyz, to_z_y_z_yyzz, to_z_y_z_yzzz, to_z_y_z_zzzz, to_zz_xxx, to_zz_xxxxy, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_yyy, to_zz_yyyyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_z_xxxx[k] = -2.0 * to_0_xxxxy[k] * tke_0 + 4.0 * to_zz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_z_xxxy[k] = to_0_xxx[k] - 2.0 * to_0_xxxyy[k] * tke_0 - 2.0 * to_zz_xxx[k] * tbe_0 + 4.0 * to_zz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_z_xxxz[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_zz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_z_xxyy[k] = 2.0 * to_0_xxy[k] - 2.0 * to_0_xxyyy[k] * tke_0 - 4.0 * to_zz_xxy[k] * tbe_0 + 4.0 * to_zz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_z_xxyz[k] = to_0_xxz[k] - 2.0 * to_0_xxyyz[k] * tke_0 - 2.0 * to_zz_xxz[k] * tbe_0 + 4.0 * to_zz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_z_xxzz[k] = -2.0 * to_0_xxyzz[k] * tke_0 + 4.0 * to_zz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_z_xyyy[k] = 3.0 * to_0_xyy[k] - 2.0 * to_0_xyyyy[k] * tke_0 - 6.0 * to_zz_xyy[k] * tbe_0 + 4.0 * to_zz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_z_xyyz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyyyz[k] * tke_0 - 4.0 * to_zz_xyz[k] * tbe_0 + 4.0 * to_zz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_z_xyzz[k] = to_0_xzz[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_zz_xzz[k] * tbe_0 + 4.0 * to_zz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_z_xzzz[k] = -2.0 * to_0_xyzzz[k] * tke_0 + 4.0 * to_zz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_z_yyyy[k] = 4.0 * to_0_yyy[k] - 2.0 * to_0_yyyyy[k] * tke_0 - 8.0 * to_zz_yyy[k] * tbe_0 + 4.0 * to_zz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_z_yyyz[k] = 3.0 * to_0_yyz[k] - 2.0 * to_0_yyyyz[k] * tke_0 - 6.0 * to_zz_yyz[k] * tbe_0 + 4.0 * to_zz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_z_yyzz[k] = 2.0 * to_0_yzz[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 4.0 * to_zz_yzz[k] * tbe_0 + 4.0 * to_zz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_z_yzzz[k] = to_0_zzz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 2.0 * to_zz_zzz[k] * tbe_0 + 4.0 * to_zz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_z_zzzz[k] = -2.0 * to_0_yzzzz[k] * tke_0 + 4.0 * to_zz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-375 components of targeted buffer : PG

        auto to_z_z_x_xxxx = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 0);

        auto to_z_z_x_xxxy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 1);

        auto to_z_z_x_xxxz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 2);

        auto to_z_z_x_xxyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 3);

        auto to_z_z_x_xxyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 4);

        auto to_z_z_x_xxzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 5);

        auto to_z_z_x_xyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 6);

        auto to_z_z_x_xyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 7);

        auto to_z_z_x_xyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 8);

        auto to_z_z_x_xzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 9);

        auto to_z_z_x_yyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 10);

        auto to_z_z_x_yyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 11);

        auto to_z_z_x_yyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 12);

        auto to_z_z_x_yzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 13);

        auto to_z_z_x_zzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, to_z_z_x_xxxx, to_z_z_x_xxxy, to_z_z_x_xxxz, to_z_z_x_xxyy, to_z_z_x_xxyz, to_z_z_x_xxzz, to_z_z_x_xyyy, to_z_z_x_xyyz, to_z_z_x_xyzz, to_z_z_x_xzzz, to_z_z_x_yyyy, to_z_z_x_yyyz, to_z_z_x_yyzz, to_z_z_x_yzzz, to_z_z_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_x_xxxx[k] = 4.0 * to_xz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_x_xxxy[k] = 4.0 * to_xz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_x_xxxz[k] = -2.0 * to_xz_xxx[k] * tbe_0 + 4.0 * to_xz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_x_xxyy[k] = 4.0 * to_xz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_x_xxyz[k] = -2.0 * to_xz_xxy[k] * tbe_0 + 4.0 * to_xz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_x_xxzz[k] = -4.0 * to_xz_xxz[k] * tbe_0 + 4.0 * to_xz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_x_xyyy[k] = 4.0 * to_xz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_x_xyyz[k] = -2.0 * to_xz_xyy[k] * tbe_0 + 4.0 * to_xz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_x_xyzz[k] = -4.0 * to_xz_xyz[k] * tbe_0 + 4.0 * to_xz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_x_xzzz[k] = -6.0 * to_xz_xzz[k] * tbe_0 + 4.0 * to_xz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_x_yyyy[k] = 4.0 * to_xz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_x_yyyz[k] = -2.0 * to_xz_yyy[k] * tbe_0 + 4.0 * to_xz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_x_yyzz[k] = -4.0 * to_xz_yyz[k] * tbe_0 + 4.0 * to_xz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_x_yzzz[k] = -6.0 * to_xz_yzz[k] * tbe_0 + 4.0 * to_xz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_x_zzzz[k] = -8.0 * to_xz_zzz[k] * tbe_0 + 4.0 * to_xz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 375-390 components of targeted buffer : PG

        auto to_z_z_y_xxxx = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 15);

        auto to_z_z_y_xxxy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 16);

        auto to_z_z_y_xxxz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 17);

        auto to_z_z_y_xxyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 18);

        auto to_z_z_y_xxyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 19);

        auto to_z_z_y_xxzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 20);

        auto to_z_z_y_xyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 21);

        auto to_z_z_y_xyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 22);

        auto to_z_z_y_xyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 23);

        auto to_z_z_y_xzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 24);

        auto to_z_z_y_yyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 25);

        auto to_z_z_y_yyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 26);

        auto to_z_z_y_yyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 27);

        auto to_z_z_y_yzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 28);

        auto to_z_z_y_zzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, to_z_z_y_xxxx, to_z_z_y_xxxy, to_z_z_y_xxxz, to_z_z_y_xxyy, to_z_z_y_xxyz, to_z_z_y_xxzz, to_z_z_y_xyyy, to_z_z_y_xyyz, to_z_z_y_xyzz, to_z_z_y_xzzz, to_z_z_y_yyyy, to_z_z_y_yyyz, to_z_z_y_yyzz, to_z_z_y_yzzz, to_z_z_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_y_xxxx[k] = 4.0 * to_yz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_y_xxxy[k] = 4.0 * to_yz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_y_xxxz[k] = -2.0 * to_yz_xxx[k] * tbe_0 + 4.0 * to_yz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_y_xxyy[k] = 4.0 * to_yz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_y_xxyz[k] = -2.0 * to_yz_xxy[k] * tbe_0 + 4.0 * to_yz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_y_xxzz[k] = -4.0 * to_yz_xxz[k] * tbe_0 + 4.0 * to_yz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_y_xyyy[k] = 4.0 * to_yz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_y_xyyz[k] = -2.0 * to_yz_xyy[k] * tbe_0 + 4.0 * to_yz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_y_xyzz[k] = -4.0 * to_yz_xyz[k] * tbe_0 + 4.0 * to_yz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_y_xzzz[k] = -6.0 * to_yz_xzz[k] * tbe_0 + 4.0 * to_yz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_y_yyyy[k] = 4.0 * to_yz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_y_yyyz[k] = -2.0 * to_yz_yyy[k] * tbe_0 + 4.0 * to_yz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_y_yyzz[k] = -4.0 * to_yz_yyz[k] * tbe_0 + 4.0 * to_yz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_y_yzzz[k] = -6.0 * to_yz_yzz[k] * tbe_0 + 4.0 * to_yz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_y_zzzz[k] = -8.0 * to_yz_zzz[k] * tbe_0 + 4.0 * to_yz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-405 components of targeted buffer : PG

        auto to_z_z_z_xxxx = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 30);

        auto to_z_z_z_xxxy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 31);

        auto to_z_z_z_xxxz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 32);

        auto to_z_z_z_xxyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 33);

        auto to_z_z_z_xxyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 34);

        auto to_z_z_z_xxzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 35);

        auto to_z_z_z_xyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 36);

        auto to_z_z_z_xyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 37);

        auto to_z_z_z_xyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 38);

        auto to_z_z_z_xzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 39);

        auto to_z_z_z_yyyy = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 40);

        auto to_z_z_z_yyyz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 41);

        auto to_z_z_z_yyzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 42);

        auto to_z_z_z_yzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 43);

        auto to_z_z_z_zzzz = pbuffer.data(idx_op_geom_101_pg + 8 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxz, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, to_0_zzzzz, to_z_z_z_xxxx, to_z_z_z_xxxy, to_z_z_z_xxxz, to_z_z_z_xxyy, to_z_z_z_xxyz, to_z_z_z_xxzz, to_z_z_z_xyyy, to_z_z_z_xyyz, to_z_z_z_xyzz, to_z_z_z_xzzz, to_z_z_z_yyyy, to_z_z_z_yyyz, to_z_z_z_yyzz, to_z_z_z_yzzz, to_z_z_z_zzzz, to_zz_xxx, to_zz_xxxxz, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_z_xxxx[k] = -2.0 * to_0_xxxxz[k] * tke_0 + 4.0 * to_zz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_z_xxxy[k] = -2.0 * to_0_xxxyz[k] * tke_0 + 4.0 * to_zz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_z_xxxz[k] = to_0_xxx[k] - 2.0 * to_0_xxxzz[k] * tke_0 - 2.0 * to_zz_xxx[k] * tbe_0 + 4.0 * to_zz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_z_xxyy[k] = -2.0 * to_0_xxyyz[k] * tke_0 + 4.0 * to_zz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_z_xxyz[k] = to_0_xxy[k] - 2.0 * to_0_xxyzz[k] * tke_0 - 2.0 * to_zz_xxy[k] * tbe_0 + 4.0 * to_zz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_z_xxzz[k] = 2.0 * to_0_xxz[k] - 2.0 * to_0_xxzzz[k] * tke_0 - 4.0 * to_zz_xxz[k] * tbe_0 + 4.0 * to_zz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_z_xyyy[k] = -2.0 * to_0_xyyyz[k] * tke_0 + 4.0 * to_zz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_z_xyyz[k] = to_0_xyy[k] - 2.0 * to_0_xyyzz[k] * tke_0 - 2.0 * to_zz_xyy[k] * tbe_0 + 4.0 * to_zz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_z_xyzz[k] = 2.0 * to_0_xyz[k] - 2.0 * to_0_xyzzz[k] * tke_0 - 4.0 * to_zz_xyz[k] * tbe_0 + 4.0 * to_zz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_z_xzzz[k] = 3.0 * to_0_xzz[k] - 2.0 * to_0_xzzzz[k] * tke_0 - 6.0 * to_zz_xzz[k] * tbe_0 + 4.0 * to_zz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_z_yyyy[k] = -2.0 * to_0_yyyyz[k] * tke_0 + 4.0 * to_zz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_z_yyyz[k] = to_0_yyy[k] - 2.0 * to_0_yyyzz[k] * tke_0 - 2.0 * to_zz_yyy[k] * tbe_0 + 4.0 * to_zz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_z_yyzz[k] = 2.0 * to_0_yyz[k] - 2.0 * to_0_yyzzz[k] * tke_0 - 4.0 * to_zz_yyz[k] * tbe_0 + 4.0 * to_zz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_z_yzzz[k] = 3.0 * to_0_yzz[k] - 2.0 * to_0_yzzzz[k] * tke_0 - 6.0 * to_zz_yzz[k] * tbe_0 + 4.0 * to_zz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_z_zzzz[k] = 4.0 * to_0_zzz[k] - 2.0 * to_0_zzzzz[k] * tke_0 - 8.0 * to_zz_zzz[k] * tbe_0 + 4.0 * to_zz_zzzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

