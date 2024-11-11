#include "GeometricalDerivatives1X1ForDG.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_dg(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_dg,
                        const size_t              idx_op_pf,
                        const size_t              idx_op_ph,
                        const size_t              idx_op_ff,
                        const size_t              idx_op_fh,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up components of auxiliary buffer : PH

        auto to_x_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 0);

        auto to_x_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 1);

        auto to_x_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 2);

        auto to_x_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 3);

        auto to_x_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 4);

        auto to_x_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 5);

        auto to_x_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 6);

        auto to_x_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 7);

        auto to_x_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 8);

        auto to_x_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 9);

        auto to_x_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 10);

        auto to_x_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 11);

        auto to_x_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 12);

        auto to_x_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 13);

        auto to_x_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 14);

        auto to_x_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 15);

        auto to_x_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 16);

        auto to_x_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 17);

        auto to_x_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 18);

        auto to_x_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 19);

        auto to_x_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 20);

        auto to_y_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 21);

        auto to_y_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 22);

        auto to_y_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 23);

        auto to_y_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 24);

        auto to_y_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 25);

        auto to_y_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 26);

        auto to_y_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 27);

        auto to_y_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 28);

        auto to_y_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 29);

        auto to_y_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 30);

        auto to_y_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 31);

        auto to_y_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 32);

        auto to_y_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 33);

        auto to_y_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 34);

        auto to_y_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 35);

        auto to_y_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 36);

        auto to_y_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 37);

        auto to_y_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 38);

        auto to_y_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 39);

        auto to_y_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 40);

        auto to_y_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 41);

        auto to_z_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 42);

        auto to_z_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 43);

        auto to_z_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 44);

        auto to_z_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 45);

        auto to_z_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 46);

        auto to_z_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 47);

        auto to_z_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 48);

        auto to_z_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 49);

        auto to_z_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 50);

        auto to_z_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 51);

        auto to_z_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 52);

        auto to_z_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 53);

        auto to_z_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 54);

        auto to_z_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 55);

        auto to_z_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 56);

        auto to_z_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 57);

        auto to_z_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 58);

        auto to_z_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 59);

        auto to_z_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 60);

        auto to_z_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 61);

        auto to_z_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 62);

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

        // Set up components of auxiliary buffer : FH

        auto to_xxx_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 0);

        auto to_xxx_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 1);

        auto to_xxx_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 2);

        auto to_xxx_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 3);

        auto to_xxx_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 4);

        auto to_xxx_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 5);

        auto to_xxx_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 6);

        auto to_xxx_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 7);

        auto to_xxx_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 8);

        auto to_xxx_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 9);

        auto to_xxx_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 10);

        auto to_xxx_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 11);

        auto to_xxx_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 12);

        auto to_xxx_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 13);

        auto to_xxx_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 14);

        auto to_xxx_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 15);

        auto to_xxx_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 16);

        auto to_xxx_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 17);

        auto to_xxx_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 18);

        auto to_xxx_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 19);

        auto to_xxx_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 20);

        auto to_xxy_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 21);

        auto to_xxy_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 22);

        auto to_xxy_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 23);

        auto to_xxy_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 24);

        auto to_xxy_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 25);

        auto to_xxy_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 26);

        auto to_xxy_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 27);

        auto to_xxy_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 28);

        auto to_xxy_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 29);

        auto to_xxy_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 30);

        auto to_xxy_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 31);

        auto to_xxy_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 32);

        auto to_xxy_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 33);

        auto to_xxy_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 34);

        auto to_xxy_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 35);

        auto to_xxy_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 36);

        auto to_xxy_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 37);

        auto to_xxy_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 38);

        auto to_xxy_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 39);

        auto to_xxy_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 40);

        auto to_xxy_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 41);

        auto to_xxz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 42);

        auto to_xxz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 43);

        auto to_xxz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 44);

        auto to_xxz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 45);

        auto to_xxz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 46);

        auto to_xxz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 47);

        auto to_xxz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 48);

        auto to_xxz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 49);

        auto to_xxz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 50);

        auto to_xxz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 51);

        auto to_xxz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 52);

        auto to_xxz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 53);

        auto to_xxz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 54);

        auto to_xxz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 55);

        auto to_xxz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 56);

        auto to_xxz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 57);

        auto to_xxz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 58);

        auto to_xxz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 59);

        auto to_xxz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 60);

        auto to_xxz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 61);

        auto to_xxz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 62);

        auto to_xyy_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 63);

        auto to_xyy_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 64);

        auto to_xyy_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 65);

        auto to_xyy_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 66);

        auto to_xyy_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 67);

        auto to_xyy_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 68);

        auto to_xyy_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 69);

        auto to_xyy_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 70);

        auto to_xyy_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 71);

        auto to_xyy_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 72);

        auto to_xyy_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 73);

        auto to_xyy_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 74);

        auto to_xyy_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 75);

        auto to_xyy_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 76);

        auto to_xyy_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 77);

        auto to_xyy_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 78);

        auto to_xyy_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 79);

        auto to_xyy_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 80);

        auto to_xyy_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 81);

        auto to_xyy_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 82);

        auto to_xyy_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 83);

        auto to_xyz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 84);

        auto to_xyz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 85);

        auto to_xyz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 86);

        auto to_xyz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 87);

        auto to_xyz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 88);

        auto to_xyz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 89);

        auto to_xyz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 90);

        auto to_xyz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 91);

        auto to_xyz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 92);

        auto to_xyz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 93);

        auto to_xyz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 94);

        auto to_xyz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 95);

        auto to_xyz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 96);

        auto to_xyz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 97);

        auto to_xyz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 98);

        auto to_xyz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 99);

        auto to_xyz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 100);

        auto to_xyz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 101);

        auto to_xyz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 102);

        auto to_xyz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 103);

        auto to_xyz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 104);

        auto to_xzz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 105);

        auto to_xzz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 106);

        auto to_xzz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 107);

        auto to_xzz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 108);

        auto to_xzz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 109);

        auto to_xzz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 110);

        auto to_xzz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 111);

        auto to_xzz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 112);

        auto to_xzz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 113);

        auto to_xzz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 114);

        auto to_xzz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 115);

        auto to_xzz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 116);

        auto to_xzz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 117);

        auto to_xzz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 118);

        auto to_xzz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 119);

        auto to_xzz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 120);

        auto to_xzz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 121);

        auto to_xzz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 122);

        auto to_xzz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 123);

        auto to_xzz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 124);

        auto to_xzz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 125);

        auto to_yyy_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 126);

        auto to_yyy_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 127);

        auto to_yyy_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 128);

        auto to_yyy_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 129);

        auto to_yyy_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 130);

        auto to_yyy_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 131);

        auto to_yyy_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 132);

        auto to_yyy_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 133);

        auto to_yyy_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 134);

        auto to_yyy_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 135);

        auto to_yyy_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 136);

        auto to_yyy_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 137);

        auto to_yyy_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 138);

        auto to_yyy_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 139);

        auto to_yyy_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 140);

        auto to_yyy_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 141);

        auto to_yyy_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 142);

        auto to_yyy_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 143);

        auto to_yyy_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 144);

        auto to_yyy_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 145);

        auto to_yyy_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 146);

        auto to_yyz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 147);

        auto to_yyz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 148);

        auto to_yyz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 149);

        auto to_yyz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 150);

        auto to_yyz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 151);

        auto to_yyz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 152);

        auto to_yyz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 153);

        auto to_yyz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 154);

        auto to_yyz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 155);

        auto to_yyz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 156);

        auto to_yyz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 157);

        auto to_yyz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 158);

        auto to_yyz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 159);

        auto to_yyz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 160);

        auto to_yyz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 161);

        auto to_yyz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 162);

        auto to_yyz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 163);

        auto to_yyz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 164);

        auto to_yyz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 165);

        auto to_yyz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 166);

        auto to_yyz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 167);

        auto to_yzz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 168);

        auto to_yzz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 169);

        auto to_yzz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 170);

        auto to_yzz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 171);

        auto to_yzz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 172);

        auto to_yzz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 173);

        auto to_yzz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 174);

        auto to_yzz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 175);

        auto to_yzz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 176);

        auto to_yzz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 177);

        auto to_yzz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 178);

        auto to_yzz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 179);

        auto to_yzz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 180);

        auto to_yzz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 181);

        auto to_yzz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 182);

        auto to_yzz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 183);

        auto to_yzz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 184);

        auto to_yzz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 185);

        auto to_yzz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 186);

        auto to_yzz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 187);

        auto to_yzz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 188);

        auto to_zzz_xxxxx = pbuffer.data(idx_op_fh + i * 210 + 189);

        auto to_zzz_xxxxy = pbuffer.data(idx_op_fh + i * 210 + 190);

        auto to_zzz_xxxxz = pbuffer.data(idx_op_fh + i * 210 + 191);

        auto to_zzz_xxxyy = pbuffer.data(idx_op_fh + i * 210 + 192);

        auto to_zzz_xxxyz = pbuffer.data(idx_op_fh + i * 210 + 193);

        auto to_zzz_xxxzz = pbuffer.data(idx_op_fh + i * 210 + 194);

        auto to_zzz_xxyyy = pbuffer.data(idx_op_fh + i * 210 + 195);

        auto to_zzz_xxyyz = pbuffer.data(idx_op_fh + i * 210 + 196);

        auto to_zzz_xxyzz = pbuffer.data(idx_op_fh + i * 210 + 197);

        auto to_zzz_xxzzz = pbuffer.data(idx_op_fh + i * 210 + 198);

        auto to_zzz_xyyyy = pbuffer.data(idx_op_fh + i * 210 + 199);

        auto to_zzz_xyyyz = pbuffer.data(idx_op_fh + i * 210 + 200);

        auto to_zzz_xyyzz = pbuffer.data(idx_op_fh + i * 210 + 201);

        auto to_zzz_xyzzz = pbuffer.data(idx_op_fh + i * 210 + 202);

        auto to_zzz_xzzzz = pbuffer.data(idx_op_fh + i * 210 + 203);

        auto to_zzz_yyyyy = pbuffer.data(idx_op_fh + i * 210 + 204);

        auto to_zzz_yyyyz = pbuffer.data(idx_op_fh + i * 210 + 205);

        auto to_zzz_yyyzz = pbuffer.data(idx_op_fh + i * 210 + 206);

        auto to_zzz_yyzzz = pbuffer.data(idx_op_fh + i * 210 + 207);

        auto to_zzz_yzzzz = pbuffer.data(idx_op_fh + i * 210 + 208);

        auto to_zzz_zzzzz = pbuffer.data(idx_op_fh + i * 210 + 209);

        // Set up 0-15 components of targeted buffer : DG

        auto to_x_x_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 0);

        auto to_x_x_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 1);

        auto to_x_x_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 2);

        auto to_x_x_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 3);

        auto to_x_x_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 4);

        auto to_x_x_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 5);

        auto to_x_x_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 6);

        auto to_x_x_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 7);

        auto to_x_x_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 8);

        auto to_x_x_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 9);

        auto to_x_x_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 10);

        auto to_x_x_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 11);

        auto to_x_x_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 12);

        auto to_x_x_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 13);

        auto to_x_x_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_x_x_xx_xxxx,     \
                             to_x_x_xx_xxxy, \
                             to_x_x_xx_xxxz, \
                             to_x_x_xx_xxyy, \
                             to_x_x_xx_xxyz, \
                             to_x_x_xx_xxzz, \
                             to_x_x_xx_xyyy, \
                             to_x_x_xx_xyyz, \
                             to_x_x_xx_xyzz, \
                             to_x_x_xx_xzzz, \
                             to_x_x_xx_yyyy, \
                             to_x_x_xx_yyyz, \
                             to_x_x_xx_yyzz, \
                             to_x_x_xx_yzzz, \
                             to_x_x_xx_zzzz, \
                             to_x_xxx,       \
                             to_x_xxxxx,     \
                             to_x_xxxxy,     \
                             to_x_xxxxz,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyz,       \
                             to_x_yzz,       \
                             to_x_zzz,       \
                             to_xxx_xxx,     \
                             to_xxx_xxxxx,   \
                             to_xxx_xxxxy,   \
                             to_xxx_xxxxz,   \
                             to_xxx_xxxyy,   \
                             to_xxx_xxxyz,   \
                             to_xxx_xxxzz,   \
                             to_xxx_xxy,     \
                             to_xxx_xxyyy,   \
                             to_xxx_xxyyz,   \
                             to_xxx_xxyzz,   \
                             to_xxx_xxz,     \
                             to_xxx_xxzzz,   \
                             to_xxx_xyy,     \
                             to_xxx_xyyyy,   \
                             to_xxx_xyyyz,   \
                             to_xxx_xyyzz,   \
                             to_xxx_xyz,     \
                             to_xxx_xyzzz,   \
                             to_xxx_xzz,     \
                             to_xxx_xzzzz,   \
                             to_xxx_yyy,     \
                             to_xxx_yyz,     \
                             to_xxx_yzz,     \
                             to_xxx_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xx_xxxx[k] = 8.0 * to_x_xxx[k] - 4.0 * to_x_xxxxx[k] * tke_0 - 8.0 * to_xxx_xxx[k] * tbe_0 + 4.0 * to_xxx_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xx_xxxy[k] = 6.0 * to_x_xxy[k] - 4.0 * to_x_xxxxy[k] * tke_0 - 6.0 * to_xxx_xxy[k] * tbe_0 + 4.0 * to_xxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xx_xxxz[k] = 6.0 * to_x_xxz[k] - 4.0 * to_x_xxxxz[k] * tke_0 - 6.0 * to_xxx_xxz[k] * tbe_0 + 4.0 * to_xxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xx_xxyy[k] = 4.0 * to_x_xyy[k] - 4.0 * to_x_xxxyy[k] * tke_0 - 4.0 * to_xxx_xyy[k] * tbe_0 + 4.0 * to_xxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xx_xxyz[k] = 4.0 * to_x_xyz[k] - 4.0 * to_x_xxxyz[k] * tke_0 - 4.0 * to_xxx_xyz[k] * tbe_0 + 4.0 * to_xxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xx_xxzz[k] = 4.0 * to_x_xzz[k] - 4.0 * to_x_xxxzz[k] * tke_0 - 4.0 * to_xxx_xzz[k] * tbe_0 + 4.0 * to_xxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xx_xyyy[k] = 2.0 * to_x_yyy[k] - 4.0 * to_x_xxyyy[k] * tke_0 - 2.0 * to_xxx_yyy[k] * tbe_0 + 4.0 * to_xxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xx_xyyz[k] = 2.0 * to_x_yyz[k] - 4.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xxx_yyz[k] * tbe_0 + 4.0 * to_xxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xx_xyzz[k] = 2.0 * to_x_yzz[k] - 4.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xxx_yzz[k] * tbe_0 + 4.0 * to_xxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xx_xzzz[k] = 2.0 * to_x_zzz[k] - 4.0 * to_x_xxzzz[k] * tke_0 - 2.0 * to_xxx_zzz[k] * tbe_0 + 4.0 * to_xxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xx_yyyy[k] = -4.0 * to_x_xyyyy[k] * tke_0 + 4.0 * to_xxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xx_yyyz[k] = -4.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xx_yyzz[k] = -4.0 * to_x_xyyzz[k] * tke_0 + 4.0 * to_xxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xx_yzzz[k] = -4.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xx_zzzz[k] = -4.0 * to_x_xzzzz[k] * tke_0 + 4.0 * to_xxx_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : DG

        auto to_x_x_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 15);

        auto to_x_x_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 16);

        auto to_x_x_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 17);

        auto to_x_x_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 18);

        auto to_x_x_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 19);

        auto to_x_x_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 20);

        auto to_x_x_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 21);

        auto to_x_x_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 22);

        auto to_x_x_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 23);

        auto to_x_x_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 24);

        auto to_x_x_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 25);

        auto to_x_x_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 26);

        auto to_x_x_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 27);

        auto to_x_x_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 28);

        auto to_x_x_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_x_xy_xxxx,     \
                             to_x_x_xy_xxxy, \
                             to_x_x_xy_xxxz, \
                             to_x_x_xy_xxyy, \
                             to_x_x_xy_xxyz, \
                             to_x_x_xy_xxzz, \
                             to_x_x_xy_xyyy, \
                             to_x_x_xy_xyyz, \
                             to_x_x_xy_xyzz, \
                             to_x_x_xy_xzzz, \
                             to_x_x_xy_yyyy, \
                             to_x_x_xy_yyyz, \
                             to_x_x_xy_yyzz, \
                             to_x_x_xy_yzzz, \
                             to_x_x_xy_zzzz, \
                             to_xxy_xxx,     \
                             to_xxy_xxxxx,   \
                             to_xxy_xxxxy,   \
                             to_xxy_xxxxz,   \
                             to_xxy_xxxyy,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxxzz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyy,   \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xxzzz,   \
                             to_xxy_xyy,     \
                             to_xxy_xyyyy,   \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_xzzzz,   \
                             to_xxy_yyy,     \
                             to_xxy_yyz,     \
                             to_xxy_yzz,     \
                             to_xxy_zzz,     \
                             to_y_xxx,       \
                             to_y_xxxxx,     \
                             to_y_xxxxy,     \
                             to_y_xxxxz,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyz,       \
                             to_y_yzz,       \
                             to_y_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xy_xxxx[k] = 4.0 * to_y_xxx[k] - 2.0 * to_y_xxxxx[k] * tke_0 - 8.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xy_xxxy[k] = 3.0 * to_y_xxy[k] - 2.0 * to_y_xxxxy[k] * tke_0 - 6.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xy_xxxz[k] = 3.0 * to_y_xxz[k] - 2.0 * to_y_xxxxz[k] * tke_0 - 6.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xy_xxyy[k] = 2.0 * to_y_xyy[k] - 2.0 * to_y_xxxyy[k] * tke_0 - 4.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xy_xxyz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xxxyz[k] * tke_0 - 4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xy_xxzz[k] = 2.0 * to_y_xzz[k] - 2.0 * to_y_xxxzz[k] * tke_0 - 4.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xy_xyyy[k] = to_y_yyy[k] - 2.0 * to_y_xxyyy[k] * tke_0 - 2.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xy_xyyz[k] = to_y_yyz[k] - 2.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xy_xyzz[k] = to_y_yzz[k] - 2.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xy_xzzz[k] = to_y_zzz[k] - 2.0 * to_y_xxzzz[k] * tke_0 - 2.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xy_yyyy[k] = -2.0 * to_y_xyyyy[k] * tke_0 + 4.0 * to_xxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xy_yyyz[k] = -2.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xy_yyzz[k] = -2.0 * to_y_xyyzz[k] * tke_0 + 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xy_yzzz[k] = -2.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xy_zzzz[k] = -2.0 * to_y_xzzzz[k] * tke_0 + 4.0 * to_xxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : DG

        auto to_x_x_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 30);

        auto to_x_x_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 31);

        auto to_x_x_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 32);

        auto to_x_x_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 33);

        auto to_x_x_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 34);

        auto to_x_x_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 35);

        auto to_x_x_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 36);

        auto to_x_x_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 37);

        auto to_x_x_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 38);

        auto to_x_x_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 39);

        auto to_x_x_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 40);

        auto to_x_x_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 41);

        auto to_x_x_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 42);

        auto to_x_x_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 43);

        auto to_x_x_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_x_xz_xxxx,     \
                             to_x_x_xz_xxxy, \
                             to_x_x_xz_xxxz, \
                             to_x_x_xz_xxyy, \
                             to_x_x_xz_xxyz, \
                             to_x_x_xz_xxzz, \
                             to_x_x_xz_xyyy, \
                             to_x_x_xz_xyyz, \
                             to_x_x_xz_xyzz, \
                             to_x_x_xz_xzzz, \
                             to_x_x_xz_yyyy, \
                             to_x_x_xz_yyyz, \
                             to_x_x_xz_yyzz, \
                             to_x_x_xz_yzzz, \
                             to_x_x_xz_zzzz, \
                             to_xxz_xxx,     \
                             to_xxz_xxxxx,   \
                             to_xxz_xxxxy,   \
                             to_xxz_xxxxz,   \
                             to_xxz_xxxyy,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxxzz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyy,   \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xxzzz,   \
                             to_xxz_xyy,     \
                             to_xxz_xyyyy,   \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_xzzzz,   \
                             to_xxz_yyy,     \
                             to_xxz_yyz,     \
                             to_xxz_yzz,     \
                             to_xxz_zzz,     \
                             to_z_xxx,       \
                             to_z_xxxxx,     \
                             to_z_xxxxy,     \
                             to_z_xxxxz,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyz,       \
                             to_z_yzz,       \
                             to_z_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xz_xxxx[k] = 4.0 * to_z_xxx[k] - 2.0 * to_z_xxxxx[k] * tke_0 - 8.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xz_xxxy[k] = 3.0 * to_z_xxy[k] - 2.0 * to_z_xxxxy[k] * tke_0 - 6.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xz_xxxz[k] = 3.0 * to_z_xxz[k] - 2.0 * to_z_xxxxz[k] * tke_0 - 6.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xz_xxyy[k] = 2.0 * to_z_xyy[k] - 2.0 * to_z_xxxyy[k] * tke_0 - 4.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xz_xxyz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xxxyz[k] * tke_0 - 4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xz_xxzz[k] = 2.0 * to_z_xzz[k] - 2.0 * to_z_xxxzz[k] * tke_0 - 4.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xz_xyyy[k] = to_z_yyy[k] - 2.0 * to_z_xxyyy[k] * tke_0 - 2.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xz_xyyz[k] = to_z_yyz[k] - 2.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xz_xyzz[k] = to_z_yzz[k] - 2.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xz_xzzz[k] = to_z_zzz[k] - 2.0 * to_z_xxzzz[k] * tke_0 - 2.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xz_yyyy[k] = -2.0 * to_z_xyyyy[k] * tke_0 + 4.0 * to_xxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xz_yyyz[k] = -2.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xz_yyzz[k] = -2.0 * to_z_xyyzz[k] * tke_0 + 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xz_yzzz[k] = -2.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xz_zzzz[k] = -2.0 * to_z_xzzzz[k] * tke_0 + 4.0 * to_xxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : DG

        auto to_x_x_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 45);

        auto to_x_x_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 46);

        auto to_x_x_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 47);

        auto to_x_x_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 48);

        auto to_x_x_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 49);

        auto to_x_x_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 50);

        auto to_x_x_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 51);

        auto to_x_x_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 52);

        auto to_x_x_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 53);

        auto to_x_x_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 54);

        auto to_x_x_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 55);

        auto to_x_x_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 56);

        auto to_x_x_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 57);

        auto to_x_x_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 58);

        auto to_x_x_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_x_x_yy_xxxx,     \
                             to_x_x_yy_xxxy, \
                             to_x_x_yy_xxxz, \
                             to_x_x_yy_xxyy, \
                             to_x_x_yy_xxyz, \
                             to_x_x_yy_xxzz, \
                             to_x_x_yy_xyyy, \
                             to_x_x_yy_xyyz, \
                             to_x_x_yy_xyzz, \
                             to_x_x_yy_xzzz, \
                             to_x_x_yy_yyyy, \
                             to_x_x_yy_yyyz, \
                             to_x_x_yy_yyzz, \
                             to_x_x_yy_yzzz, \
                             to_x_x_yy_zzzz, \
                             to_xyy_xxx,     \
                             to_xyy_xxxxx,   \
                             to_xyy_xxxxy,   \
                             to_xyy_xxxxz,   \
                             to_xyy_xxxyy,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxxzz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyy,   \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xxzzz,   \
                             to_xyy_xyy,     \
                             to_xyy_xyyyy,   \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_xzzzz,   \
                             to_xyy_yyy,     \
                             to_xyy_yyz,     \
                             to_xyy_yzz,     \
                             to_xyy_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yy_xxxx[k] = -8.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yy_xxxy[k] = -6.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yy_xxxz[k] = -6.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yy_xxyy[k] = -4.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yy_xxyz[k] = -4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yy_xxzz[k] = -4.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yy_xyyy[k] = -2.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yy_xyyz[k] = -2.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yy_xyzz[k] = -2.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yy_xzzz[k] = -2.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yy_yyyy[k] = 4.0 * to_xyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yy_yyyz[k] = 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yy_yyzz[k] = 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yy_yzzz[k] = 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yy_zzzz[k] = 4.0 * to_xyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : DG

        auto to_x_x_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 60);

        auto to_x_x_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 61);

        auto to_x_x_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 62);

        auto to_x_x_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 63);

        auto to_x_x_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 64);

        auto to_x_x_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 65);

        auto to_x_x_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 66);

        auto to_x_x_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 67);

        auto to_x_x_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 68);

        auto to_x_x_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 69);

        auto to_x_x_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 70);

        auto to_x_x_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 71);

        auto to_x_x_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 72);

        auto to_x_x_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 73);

        auto to_x_x_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_x_x_yz_xxxx,     \
                             to_x_x_yz_xxxy, \
                             to_x_x_yz_xxxz, \
                             to_x_x_yz_xxyy, \
                             to_x_x_yz_xxyz, \
                             to_x_x_yz_xxzz, \
                             to_x_x_yz_xyyy, \
                             to_x_x_yz_xyyz, \
                             to_x_x_yz_xyzz, \
                             to_x_x_yz_xzzz, \
                             to_x_x_yz_yyyy, \
                             to_x_x_yz_yyyz, \
                             to_x_x_yz_yyzz, \
                             to_x_x_yz_yzzz, \
                             to_x_x_yz_zzzz, \
                             to_xyz_xxx,     \
                             to_xyz_xxxxx,   \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyz,     \
                             to_xyz_yzz,     \
                             to_xyz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yz_xxxx[k] = -8.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yz_xxxy[k] = -6.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yz_xxxz[k] = -6.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yz_xxyy[k] = -4.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yz_xxyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yz_xxzz[k] = -4.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yz_xyyy[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yz_xyyz[k] = -2.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yz_xyzz[k] = -2.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yz_xzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yz_yyyy[k] = 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yz_yyyz[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yz_yyzz[k] = 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yz_yzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yz_zzzz[k] = 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : DG

        auto to_x_x_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 75);

        auto to_x_x_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 76);

        auto to_x_x_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 77);

        auto to_x_x_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 78);

        auto to_x_x_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 79);

        auto to_x_x_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 80);

        auto to_x_x_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 81);

        auto to_x_x_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 82);

        auto to_x_x_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 83);

        auto to_x_x_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 84);

        auto to_x_x_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 85);

        auto to_x_x_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 86);

        auto to_x_x_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 87);

        auto to_x_x_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 88);

        auto to_x_x_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 0 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_x_x_zz_xxxx,     \
                             to_x_x_zz_xxxy, \
                             to_x_x_zz_xxxz, \
                             to_x_x_zz_xxyy, \
                             to_x_x_zz_xxyz, \
                             to_x_x_zz_xxzz, \
                             to_x_x_zz_xyyy, \
                             to_x_x_zz_xyyz, \
                             to_x_x_zz_xyzz, \
                             to_x_x_zz_xzzz, \
                             to_x_x_zz_yyyy, \
                             to_x_x_zz_yyyz, \
                             to_x_x_zz_yyzz, \
                             to_x_x_zz_yzzz, \
                             to_x_x_zz_zzzz, \
                             to_xzz_xxx,     \
                             to_xzz_xxxxx,   \
                             to_xzz_xxxxy,   \
                             to_xzz_xxxxz,   \
                             to_xzz_xxxyy,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxxzz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyy,   \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xxzzz,   \
                             to_xzz_xyy,     \
                             to_xzz_xyyyy,   \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_xzzzz,   \
                             to_xzz_yyy,     \
                             to_xzz_yyz,     \
                             to_xzz_yzz,     \
                             to_xzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zz_xxxx[k] = -8.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_zz_xxxy[k] = -6.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_zz_xxxz[k] = -6.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_zz_xxyy[k] = -4.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_zz_xxyz[k] = -4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_zz_xxzz[k] = -4.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_zz_xyyy[k] = -2.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_zz_xyyz[k] = -2.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_zz_xyzz[k] = -2.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_zz_xzzz[k] = -2.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_zz_yyyy[k] = 4.0 * to_xzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_zz_yyyz[k] = 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_zz_yyzz[k] = 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_zz_yzzz[k] = 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_zz_zzzz[k] = 4.0 * to_xzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : DG

        auto to_x_y_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 0);

        auto to_x_y_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 1);

        auto to_x_y_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 2);

        auto to_x_y_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 3);

        auto to_x_y_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 4);

        auto to_x_y_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 5);

        auto to_x_y_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 6);

        auto to_x_y_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 7);

        auto to_x_y_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 8);

        auto to_x_y_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 9);

        auto to_x_y_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 10);

        auto to_x_y_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 11);

        auto to_x_y_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 12);

        auto to_x_y_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 13);

        auto to_x_y_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxy,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_y_xx_xxxx, \
                             to_x_y_xx_xxxy, \
                             to_x_y_xx_xxxz, \
                             to_x_y_xx_xxyy, \
                             to_x_y_xx_xxyz, \
                             to_x_y_xx_xxzz, \
                             to_x_y_xx_xyyy, \
                             to_x_y_xx_xyyz, \
                             to_x_y_xx_xyzz, \
                             to_x_y_xx_xzzz, \
                             to_x_y_xx_yyyy, \
                             to_x_y_xx_yyyz, \
                             to_x_y_xx_yyzz, \
                             to_x_y_xx_yzzz, \
                             to_x_y_xx_zzzz, \
                             to_x_yyy,       \
                             to_x_yyyyy,     \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_zzz,       \
                             to_xxx_xxx,     \
                             to_xxx_xxxxy,   \
                             to_xxx_xxxyy,   \
                             to_xxx_xxxyz,   \
                             to_xxx_xxy,     \
                             to_xxx_xxyyy,   \
                             to_xxx_xxyyz,   \
                             to_xxx_xxyzz,   \
                             to_xxx_xxz,     \
                             to_xxx_xyy,     \
                             to_xxx_xyyyy,   \
                             to_xxx_xyyyz,   \
                             to_xxx_xyyzz,   \
                             to_xxx_xyz,     \
                             to_xxx_xyzzz,   \
                             to_xxx_xzz,     \
                             to_xxx_yyy,     \
                             to_xxx_yyyyy,   \
                             to_xxx_yyyyz,   \
                             to_xxx_yyyzz,   \
                             to_xxx_yyz,     \
                             to_xxx_yyzzz,   \
                             to_xxx_yzz,     \
                             to_xxx_yzzzz,   \
                             to_xxx_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xx_xxxx[k] = -4.0 * to_x_xxxxy[k] * tke_0 + 4.0 * to_xxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xx_xxxy[k] = 2.0 * to_x_xxx[k] - 4.0 * to_x_xxxyy[k] * tke_0 - 2.0 * to_xxx_xxx[k] * tbe_0 + 4.0 * to_xxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xxxz[k] = -4.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xx_xxyy[k] = 4.0 * to_x_xxy[k] - 4.0 * to_x_xxyyy[k] * tke_0 - 4.0 * to_xxx_xxy[k] * tbe_0 + 4.0 * to_xxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xxyz[k] = 2.0 * to_x_xxz[k] - 4.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xxx_xxz[k] * tbe_0 + 4.0 * to_xxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xx_xxzz[k] = -4.0 * to_x_xxyzz[k] * tke_0 + 4.0 * to_xxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xx_xyyy[k] = 6.0 * to_x_xyy[k] - 4.0 * to_x_xyyyy[k] * tke_0 - 6.0 * to_xxx_xyy[k] * tbe_0 + 4.0 * to_xxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xyyz[k] = 4.0 * to_x_xyz[k] - 4.0 * to_x_xyyyz[k] * tke_0 - 4.0 * to_xxx_xyz[k] * tbe_0 + 4.0 * to_xxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xx_xyzz[k] = 2.0 * to_x_xzz[k] - 4.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xxx_xzz[k] * tbe_0 + 4.0 * to_xxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xx_xzzz[k] = -4.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xx_yyyy[k] = 8.0 * to_x_yyy[k] - 4.0 * to_x_yyyyy[k] * tke_0 - 8.0 * to_xxx_yyy[k] * tbe_0 + 4.0 * to_xxx_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xx_yyyz[k] = 6.0 * to_x_yyz[k] - 4.0 * to_x_yyyyz[k] * tke_0 - 6.0 * to_xxx_yyz[k] * tbe_0 + 4.0 * to_xxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xx_yyzz[k] = 4.0 * to_x_yzz[k] - 4.0 * to_x_yyyzz[k] * tke_0 - 4.0 * to_xxx_yzz[k] * tbe_0 + 4.0 * to_xxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xx_yzzz[k] = 2.0 * to_x_zzz[k] - 4.0 * to_x_yyzzz[k] * tke_0 - 2.0 * to_xxx_zzz[k] * tbe_0 + 4.0 * to_xxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xx_zzzz[k] = -4.0 * to_x_yzzzz[k] * tke_0 + 4.0 * to_xxx_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : DG

        auto to_x_y_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 15);

        auto to_x_y_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 16);

        auto to_x_y_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 17);

        auto to_x_y_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 18);

        auto to_x_y_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 19);

        auto to_x_y_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 20);

        auto to_x_y_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 21);

        auto to_x_y_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 22);

        auto to_x_y_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 23);

        auto to_x_y_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 24);

        auto to_x_y_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 25);

        auto to_x_y_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 26);

        auto to_x_y_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 27);

        auto to_x_y_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 28);

        auto to_x_y_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_y_xy_xxxx,     \
                             to_x_y_xy_xxxy, \
                             to_x_y_xy_xxxz, \
                             to_x_y_xy_xxyy, \
                             to_x_y_xy_xxyz, \
                             to_x_y_xy_xxzz, \
                             to_x_y_xy_xyyy, \
                             to_x_y_xy_xyyz, \
                             to_x_y_xy_xyzz, \
                             to_x_y_xy_xzzz, \
                             to_x_y_xy_yyyy, \
                             to_x_y_xy_yyyz, \
                             to_x_y_xy_yyzz, \
                             to_x_y_xy_yzzz, \
                             to_x_y_xy_zzzz, \
                             to_xxy_xxx,     \
                             to_xxy_xxxxy,   \
                             to_xxy_xxxyy,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyy,   \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xyy,     \
                             to_xxy_xyyyy,   \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_yyy,     \
                             to_xxy_yyyyy,   \
                             to_xxy_yyyyz,   \
                             to_xxy_yyyzz,   \
                             to_xxy_yyz,     \
                             to_xxy_yyzzz,   \
                             to_xxy_yzz,     \
                             to_xxy_yzzzz,   \
                             to_xxy_zzz,     \
                             to_y_xxx,       \
                             to_y_xxxxy,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_yyy,       \
                             to_y_yyyyy,     \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xy_xxxx[k] = -2.0 * to_y_xxxxy[k] * tke_0 + 4.0 * to_xxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xy_xxxy[k] = to_y_xxx[k] - 2.0 * to_y_xxxyy[k] * tke_0 - 2.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xxxz[k] = -2.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xy_xxyy[k] = 2.0 * to_y_xxy[k] - 2.0 * to_y_xxyyy[k] * tke_0 - 4.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xxyz[k] = to_y_xxz[k] - 2.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xy_xxzz[k] = -2.0 * to_y_xxyzz[k] * tke_0 + 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xy_xyyy[k] = 3.0 * to_y_xyy[k] - 2.0 * to_y_xyyyy[k] * tke_0 - 6.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xyyz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xyyyz[k] * tke_0 - 4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xy_xyzz[k] = to_y_xzz[k] - 2.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xy_xzzz[k] = -2.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xy_yyyy[k] = 4.0 * to_y_yyy[k] - 2.0 * to_y_yyyyy[k] * tke_0 - 8.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xy_yyyz[k] = 3.0 * to_y_yyz[k] - 2.0 * to_y_yyyyz[k] * tke_0 - 6.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xy_yyzz[k] = 2.0 * to_y_yzz[k] - 2.0 * to_y_yyyzz[k] * tke_0 - 4.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xy_yzzz[k] = to_y_zzz[k] - 2.0 * to_y_yyzzz[k] * tke_0 - 2.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xy_zzzz[k] = -2.0 * to_y_yzzzz[k] * tke_0 + 4.0 * to_xxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : DG

        auto to_x_y_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 30);

        auto to_x_y_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 31);

        auto to_x_y_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 32);

        auto to_x_y_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 33);

        auto to_x_y_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 34);

        auto to_x_y_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 35);

        auto to_x_y_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 36);

        auto to_x_y_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 37);

        auto to_x_y_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 38);

        auto to_x_y_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 39);

        auto to_x_y_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 40);

        auto to_x_y_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 41);

        auto to_x_y_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 42);

        auto to_x_y_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 43);

        auto to_x_y_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_y_xz_xxxx,     \
                             to_x_y_xz_xxxy, \
                             to_x_y_xz_xxxz, \
                             to_x_y_xz_xxyy, \
                             to_x_y_xz_xxyz, \
                             to_x_y_xz_xxzz, \
                             to_x_y_xz_xyyy, \
                             to_x_y_xz_xyyz, \
                             to_x_y_xz_xyzz, \
                             to_x_y_xz_xzzz, \
                             to_x_y_xz_yyyy, \
                             to_x_y_xz_yyyz, \
                             to_x_y_xz_yyzz, \
                             to_x_y_xz_yzzz, \
                             to_x_y_xz_zzzz, \
                             to_xxz_xxx,     \
                             to_xxz_xxxxy,   \
                             to_xxz_xxxyy,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyy,   \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xyy,     \
                             to_xxz_xyyyy,   \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_yyy,     \
                             to_xxz_yyyyy,   \
                             to_xxz_yyyyz,   \
                             to_xxz_yyyzz,   \
                             to_xxz_yyz,     \
                             to_xxz_yyzzz,   \
                             to_xxz_yzz,     \
                             to_xxz_yzzzz,   \
                             to_xxz_zzz,     \
                             to_z_xxx,       \
                             to_z_xxxxy,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_yyy,       \
                             to_z_yyyyy,     \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xz_xxxx[k] = -2.0 * to_z_xxxxy[k] * tke_0 + 4.0 * to_xxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xz_xxxy[k] = to_z_xxx[k] - 2.0 * to_z_xxxyy[k] * tke_0 - 2.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xxxz[k] = -2.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xz_xxyy[k] = 2.0 * to_z_xxy[k] - 2.0 * to_z_xxyyy[k] * tke_0 - 4.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xxyz[k] = to_z_xxz[k] - 2.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xz_xxzz[k] = -2.0 * to_z_xxyzz[k] * tke_0 + 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xz_xyyy[k] = 3.0 * to_z_xyy[k] - 2.0 * to_z_xyyyy[k] * tke_0 - 6.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xyyz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xyyyz[k] * tke_0 - 4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xz_xyzz[k] = to_z_xzz[k] - 2.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xz_xzzz[k] = -2.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xz_yyyy[k] = 4.0 * to_z_yyy[k] - 2.0 * to_z_yyyyy[k] * tke_0 - 8.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xz_yyyz[k] = 3.0 * to_z_yyz[k] - 2.0 * to_z_yyyyz[k] * tke_0 - 6.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xz_yyzz[k] = 2.0 * to_z_yzz[k] - 2.0 * to_z_yyyzz[k] * tke_0 - 4.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xz_yzzz[k] = to_z_zzz[k] - 2.0 * to_z_yyzzz[k] * tke_0 - 2.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xz_zzzz[k] = -2.0 * to_z_yzzzz[k] * tke_0 + 4.0 * to_xxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 135-150 components of targeted buffer : DG

        auto to_x_y_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 45);

        auto to_x_y_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 46);

        auto to_x_y_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 47);

        auto to_x_y_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 48);

        auto to_x_y_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 49);

        auto to_x_y_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 50);

        auto to_x_y_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 51);

        auto to_x_y_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 52);

        auto to_x_y_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 53);

        auto to_x_y_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 54);

        auto to_x_y_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 55);

        auto to_x_y_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 56);

        auto to_x_y_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 57);

        auto to_x_y_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 58);

        auto to_x_y_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_x_y_yy_xxxx,     \
                             to_x_y_yy_xxxy, \
                             to_x_y_yy_xxxz, \
                             to_x_y_yy_xxyy, \
                             to_x_y_yy_xxyz, \
                             to_x_y_yy_xxzz, \
                             to_x_y_yy_xyyy, \
                             to_x_y_yy_xyyz, \
                             to_x_y_yy_xyzz, \
                             to_x_y_yy_xzzz, \
                             to_x_y_yy_yyyy, \
                             to_x_y_yy_yyyz, \
                             to_x_y_yy_yyzz, \
                             to_x_y_yy_yzzz, \
                             to_x_y_yy_zzzz, \
                             to_xyy_xxx,     \
                             to_xyy_xxxxy,   \
                             to_xyy_xxxyy,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyy,   \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xyy,     \
                             to_xyy_xyyyy,   \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_yyy,     \
                             to_xyy_yyyyy,   \
                             to_xyy_yyyyz,   \
                             to_xyy_yyyzz,   \
                             to_xyy_yyz,     \
                             to_xyy_yyzzz,   \
                             to_xyy_yzz,     \
                             to_xyy_yzzzz,   \
                             to_xyy_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yy_xxxx[k] = 4.0 * to_xyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yy_xxxy[k] = -2.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xxxz[k] = 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yy_xxyy[k] = -4.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xxyz[k] = -2.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yy_xxzz[k] = 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yy_xyyy[k] = -6.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xyyz[k] = -4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yy_xyzz[k] = -2.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yy_xzzz[k] = 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yy_yyyy[k] = -8.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yy_yyyz[k] = -6.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yy_yyzz[k] = -4.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yy_yzzz[k] = -2.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yy_zzzz[k] = 4.0 * to_xyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-165 components of targeted buffer : DG

        auto to_x_y_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 60);

        auto to_x_y_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 61);

        auto to_x_y_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 62);

        auto to_x_y_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 63);

        auto to_x_y_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 64);

        auto to_x_y_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 65);

        auto to_x_y_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 66);

        auto to_x_y_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 67);

        auto to_x_y_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 68);

        auto to_x_y_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 69);

        auto to_x_y_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 70);

        auto to_x_y_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 71);

        auto to_x_y_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 72);

        auto to_x_y_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 73);

        auto to_x_y_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_x_y_yz_xxxx,     \
                             to_x_y_yz_xxxy, \
                             to_x_y_yz_xxxz, \
                             to_x_y_yz_xxyy, \
                             to_x_y_yz_xxyz, \
                             to_x_y_yz_xxzz, \
                             to_x_y_yz_xyyy, \
                             to_x_y_yz_xyyz, \
                             to_x_y_yz_xyzz, \
                             to_x_y_yz_xzzz, \
                             to_x_y_yz_yyyy, \
                             to_x_y_yz_yyyz, \
                             to_x_y_yz_yyzz, \
                             to_x_y_yz_yzzz, \
                             to_x_y_yz_zzzz, \
                             to_xyz_xxx,     \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_yyy,     \
                             to_xyz_yyyyy,   \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yz_xxxx[k] = 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yz_xxxy[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xxxz[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yz_xxyy[k] = -4.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xxyz[k] = -2.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yz_xxzz[k] = 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yz_xyyy[k] = -6.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xyyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yz_xyzz[k] = -2.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yz_xzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yz_yyyy[k] = -8.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yz_yyyz[k] = -6.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yz_yyzz[k] = -4.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yz_yzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yz_zzzz[k] = 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 165-180 components of targeted buffer : DG

        auto to_x_y_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 75);

        auto to_x_y_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 76);

        auto to_x_y_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 77);

        auto to_x_y_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 78);

        auto to_x_y_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 79);

        auto to_x_y_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 80);

        auto to_x_y_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 81);

        auto to_x_y_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 82);

        auto to_x_y_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 83);

        auto to_x_y_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 84);

        auto to_x_y_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 85);

        auto to_x_y_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 86);

        auto to_x_y_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 87);

        auto to_x_y_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 88);

        auto to_x_y_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 1 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_x_y_zz_xxxx,     \
                             to_x_y_zz_xxxy, \
                             to_x_y_zz_xxxz, \
                             to_x_y_zz_xxyy, \
                             to_x_y_zz_xxyz, \
                             to_x_y_zz_xxzz, \
                             to_x_y_zz_xyyy, \
                             to_x_y_zz_xyyz, \
                             to_x_y_zz_xyzz, \
                             to_x_y_zz_xzzz, \
                             to_x_y_zz_yyyy, \
                             to_x_y_zz_yyyz, \
                             to_x_y_zz_yyzz, \
                             to_x_y_zz_yzzz, \
                             to_x_y_zz_zzzz, \
                             to_xzz_xxx,     \
                             to_xzz_xxxxy,   \
                             to_xzz_xxxyy,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyy,   \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xyy,     \
                             to_xzz_xyyyy,   \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_yyy,     \
                             to_xzz_yyyyy,   \
                             to_xzz_yyyyz,   \
                             to_xzz_yyyzz,   \
                             to_xzz_yyz,     \
                             to_xzz_yyzzz,   \
                             to_xzz_yzz,     \
                             to_xzz_yzzzz,   \
                             to_xzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zz_xxxx[k] = 4.0 * to_xzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_zz_xxxy[k] = -2.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xxxz[k] = 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_zz_xxyy[k] = -4.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xxyz[k] = -2.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_zz_xxzz[k] = 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_zz_xyyy[k] = -6.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xyyz[k] = -4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_zz_xyzz[k] = -2.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_zz_xzzz[k] = 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_zz_yyyy[k] = -8.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_zz_yyyz[k] = -6.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_zz_yyzz[k] = -4.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_zz_yzzz[k] = -2.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_zz_zzzz[k] = 4.0 * to_xzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-195 components of targeted buffer : DG

        auto to_x_z_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 0);

        auto to_x_z_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 1);

        auto to_x_z_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 2);

        auto to_x_z_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 3);

        auto to_x_z_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 4);

        auto to_x_z_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 5);

        auto to_x_z_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 6);

        auto to_x_z_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 7);

        auto to_x_z_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 8);

        auto to_x_z_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 9);

        auto to_x_z_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 10);

        auto to_x_z_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 11);

        auto to_x_z_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 12);

        auto to_x_z_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 13);

        auto to_x_z_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxz,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_z_xx_xxxx, \
                             to_x_z_xx_xxxy, \
                             to_x_z_xx_xxxz, \
                             to_x_z_xx_xxyy, \
                             to_x_z_xx_xxyz, \
                             to_x_z_xx_xxzz, \
                             to_x_z_xx_xyyy, \
                             to_x_z_xx_xyyz, \
                             to_x_z_xx_xyzz, \
                             to_x_z_xx_xzzz, \
                             to_x_z_xx_yyyy, \
                             to_x_z_xx_yyyz, \
                             to_x_z_xx_yyzz, \
                             to_x_z_xx_yzzz, \
                             to_x_z_xx_zzzz, \
                             to_x_zzz,       \
                             to_x_zzzzz,     \
                             to_xxx_xxx,     \
                             to_xxx_xxxxz,   \
                             to_xxx_xxxyz,   \
                             to_xxx_xxxzz,   \
                             to_xxx_xxy,     \
                             to_xxx_xxyyz,   \
                             to_xxx_xxyzz,   \
                             to_xxx_xxz,     \
                             to_xxx_xxzzz,   \
                             to_xxx_xyy,     \
                             to_xxx_xyyyz,   \
                             to_xxx_xyyzz,   \
                             to_xxx_xyz,     \
                             to_xxx_xyzzz,   \
                             to_xxx_xzz,     \
                             to_xxx_xzzzz,   \
                             to_xxx_yyy,     \
                             to_xxx_yyyyz,   \
                             to_xxx_yyyzz,   \
                             to_xxx_yyz,     \
                             to_xxx_yyzzz,   \
                             to_xxx_yzz,     \
                             to_xxx_yzzzz,   \
                             to_xxx_zzz,     \
                             to_xxx_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xx_xxxx[k] = -4.0 * to_x_xxxxz[k] * tke_0 + 4.0 * to_xxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxxy[k] = -4.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxxz[k] = 2.0 * to_x_xxx[k] - 4.0 * to_x_xxxzz[k] * tke_0 - 2.0 * to_xxx_xxx[k] * tbe_0 + 4.0 * to_xxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxyy[k] = -4.0 * to_x_xxyyz[k] * tke_0 + 4.0 * to_xxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxyz[k] = 2.0 * to_x_xxy[k] - 4.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xxx_xxy[k] * tbe_0 + 4.0 * to_xxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxzz[k] = 4.0 * to_x_xxz[k] - 4.0 * to_x_xxzzz[k] * tke_0 - 4.0 * to_xxx_xxz[k] * tbe_0 + 4.0 * to_xxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xyyy[k] = -4.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xyyz[k] = 2.0 * to_x_xyy[k] - 4.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xxx_xyy[k] * tbe_0 + 4.0 * to_xxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xyzz[k] = 4.0 * to_x_xyz[k] - 4.0 * to_x_xyzzz[k] * tke_0 - 4.0 * to_xxx_xyz[k] * tbe_0 + 4.0 * to_xxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xzzz[k] = 6.0 * to_x_xzz[k] - 4.0 * to_x_xzzzz[k] * tke_0 - 6.0 * to_xxx_xzz[k] * tbe_0 + 4.0 * to_xxx_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yyyy[k] = -4.0 * to_x_yyyyz[k] * tke_0 + 4.0 * to_xxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xx_yyyz[k] = 2.0 * to_x_yyy[k] - 4.0 * to_x_yyyzz[k] * tke_0 - 2.0 * to_xxx_yyy[k] * tbe_0 + 4.0 * to_xxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yyzz[k] = 4.0 * to_x_yyz[k] - 4.0 * to_x_yyzzz[k] * tke_0 - 4.0 * to_xxx_yyz[k] * tbe_0 + 4.0 * to_xxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yzzz[k] = 6.0 * to_x_yzz[k] - 4.0 * to_x_yzzzz[k] * tke_0 - 6.0 * to_xxx_yzz[k] * tbe_0 + 4.0 * to_xxx_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_zzzz[k] = 8.0 * to_x_zzz[k] - 4.0 * to_x_zzzzz[k] * tke_0 - 8.0 * to_xxx_zzz[k] * tbe_0 + 4.0 * to_xxx_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 195-210 components of targeted buffer : DG

        auto to_x_z_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 15);

        auto to_x_z_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 16);

        auto to_x_z_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 17);

        auto to_x_z_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 18);

        auto to_x_z_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 19);

        auto to_x_z_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 20);

        auto to_x_z_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 21);

        auto to_x_z_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 22);

        auto to_x_z_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 23);

        auto to_x_z_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 24);

        auto to_x_z_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 25);

        auto to_x_z_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 26);

        auto to_x_z_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 27);

        auto to_x_z_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 28);

        auto to_x_z_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_z_xy_xxxx,     \
                             to_x_z_xy_xxxy, \
                             to_x_z_xy_xxxz, \
                             to_x_z_xy_xxyy, \
                             to_x_z_xy_xxyz, \
                             to_x_z_xy_xxzz, \
                             to_x_z_xy_xyyy, \
                             to_x_z_xy_xyyz, \
                             to_x_z_xy_xyzz, \
                             to_x_z_xy_xzzz, \
                             to_x_z_xy_yyyy, \
                             to_x_z_xy_yyyz, \
                             to_x_z_xy_yyzz, \
                             to_x_z_xy_yzzz, \
                             to_x_z_xy_zzzz, \
                             to_xxy_xxx,     \
                             to_xxy_xxxxz,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxxzz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xxzzz,   \
                             to_xxy_xyy,     \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_xzzzz,   \
                             to_xxy_yyy,     \
                             to_xxy_yyyyz,   \
                             to_xxy_yyyzz,   \
                             to_xxy_yyz,     \
                             to_xxy_yyzzz,   \
                             to_xxy_yzz,     \
                             to_xxy_yzzzz,   \
                             to_xxy_zzz,     \
                             to_xxy_zzzzz,   \
                             to_y_xxx,       \
                             to_y_xxxxz,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_zzz,       \
                             to_y_zzzzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xy_xxxx[k] = -2.0 * to_y_xxxxz[k] * tke_0 + 4.0 * to_xxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxxy[k] = -2.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxxz[k] = to_y_xxx[k] - 2.0 * to_y_xxxzz[k] * tke_0 - 2.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxyy[k] = -2.0 * to_y_xxyyz[k] * tke_0 + 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxyz[k] = to_y_xxy[k] - 2.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxzz[k] = 2.0 * to_y_xxz[k] - 2.0 * to_y_xxzzz[k] * tke_0 - 4.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xyyy[k] = -2.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xyyz[k] = to_y_xyy[k] - 2.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xyzz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xyzzz[k] * tke_0 - 4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xzzz[k] = 3.0 * to_y_xzz[k] - 2.0 * to_y_xzzzz[k] * tke_0 - 6.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yyyy[k] = -2.0 * to_y_yyyyz[k] * tke_0 + 4.0 * to_xxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xy_yyyz[k] = to_y_yyy[k] - 2.0 * to_y_yyyzz[k] * tke_0 - 2.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yyzz[k] = 2.0 * to_y_yyz[k] - 2.0 * to_y_yyzzz[k] * tke_0 - 4.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yzzz[k] = 3.0 * to_y_yzz[k] - 2.0 * to_y_yzzzz[k] * tke_0 - 6.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_zzzz[k] = 4.0 * to_y_zzz[k] - 2.0 * to_y_zzzzz[k] * tke_0 - 8.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-225 components of targeted buffer : DG

        auto to_x_z_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 30);

        auto to_x_z_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 31);

        auto to_x_z_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 32);

        auto to_x_z_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 33);

        auto to_x_z_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 34);

        auto to_x_z_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 35);

        auto to_x_z_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 36);

        auto to_x_z_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 37);

        auto to_x_z_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 38);

        auto to_x_z_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 39);

        auto to_x_z_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 40);

        auto to_x_z_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 41);

        auto to_x_z_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 42);

        auto to_x_z_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 43);

        auto to_x_z_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_z_xz_xxxx,     \
                             to_x_z_xz_xxxy, \
                             to_x_z_xz_xxxz, \
                             to_x_z_xz_xxyy, \
                             to_x_z_xz_xxyz, \
                             to_x_z_xz_xxzz, \
                             to_x_z_xz_xyyy, \
                             to_x_z_xz_xyyz, \
                             to_x_z_xz_xyzz, \
                             to_x_z_xz_xzzz, \
                             to_x_z_xz_yyyy, \
                             to_x_z_xz_yyyz, \
                             to_x_z_xz_yyzz, \
                             to_x_z_xz_yzzz, \
                             to_x_z_xz_zzzz, \
                             to_xxz_xxx,     \
                             to_xxz_xxxxz,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxxzz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xxzzz,   \
                             to_xxz_xyy,     \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_xzzzz,   \
                             to_xxz_yyy,     \
                             to_xxz_yyyyz,   \
                             to_xxz_yyyzz,   \
                             to_xxz_yyz,     \
                             to_xxz_yyzzz,   \
                             to_xxz_yzz,     \
                             to_xxz_yzzzz,   \
                             to_xxz_zzz,     \
                             to_xxz_zzzzz,   \
                             to_z_xxx,       \
                             to_z_xxxxz,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_zzz,       \
                             to_z_zzzzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xz_xxxx[k] = -2.0 * to_z_xxxxz[k] * tke_0 + 4.0 * to_xxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxxy[k] = -2.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxxz[k] = to_z_xxx[k] - 2.0 * to_z_xxxzz[k] * tke_0 - 2.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxyy[k] = -2.0 * to_z_xxyyz[k] * tke_0 + 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxyz[k] = to_z_xxy[k] - 2.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxzz[k] = 2.0 * to_z_xxz[k] - 2.0 * to_z_xxzzz[k] * tke_0 - 4.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xyyy[k] = -2.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xyyz[k] = to_z_xyy[k] - 2.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xyzz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xyzzz[k] * tke_0 - 4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xzzz[k] = 3.0 * to_z_xzz[k] - 2.0 * to_z_xzzzz[k] * tke_0 - 6.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yyyy[k] = -2.0 * to_z_yyyyz[k] * tke_0 + 4.0 * to_xxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xz_yyyz[k] = to_z_yyy[k] - 2.0 * to_z_yyyzz[k] * tke_0 - 2.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yyzz[k] = 2.0 * to_z_yyz[k] - 2.0 * to_z_yyzzz[k] * tke_0 - 4.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yzzz[k] = 3.0 * to_z_yzz[k] - 2.0 * to_z_yzzzz[k] * tke_0 - 6.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_zzzz[k] = 4.0 * to_z_zzz[k] - 2.0 * to_z_zzzzz[k] * tke_0 - 8.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 225-240 components of targeted buffer : DG

        auto to_x_z_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 45);

        auto to_x_z_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 46);

        auto to_x_z_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 47);

        auto to_x_z_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 48);

        auto to_x_z_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 49);

        auto to_x_z_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 50);

        auto to_x_z_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 51);

        auto to_x_z_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 52);

        auto to_x_z_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 53);

        auto to_x_z_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 54);

        auto to_x_z_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 55);

        auto to_x_z_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 56);

        auto to_x_z_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 57);

        auto to_x_z_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 58);

        auto to_x_z_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_x_z_yy_xxxx,     \
                             to_x_z_yy_xxxy, \
                             to_x_z_yy_xxxz, \
                             to_x_z_yy_xxyy, \
                             to_x_z_yy_xxyz, \
                             to_x_z_yy_xxzz, \
                             to_x_z_yy_xyyy, \
                             to_x_z_yy_xyyz, \
                             to_x_z_yy_xyzz, \
                             to_x_z_yy_xzzz, \
                             to_x_z_yy_yyyy, \
                             to_x_z_yy_yyyz, \
                             to_x_z_yy_yyzz, \
                             to_x_z_yy_yzzz, \
                             to_x_z_yy_zzzz, \
                             to_xyy_xxx,     \
                             to_xyy_xxxxz,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxxzz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xxzzz,   \
                             to_xyy_xyy,     \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_xzzzz,   \
                             to_xyy_yyy,     \
                             to_xyy_yyyyz,   \
                             to_xyy_yyyzz,   \
                             to_xyy_yyz,     \
                             to_xyy_yyzzz,   \
                             to_xyy_yzz,     \
                             to_xyy_yzzzz,   \
                             to_xyy_zzz,     \
                             to_xyy_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yy_xxxx[k] = 4.0 * to_xyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxxy[k] = 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxxz[k] = -2.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxyy[k] = 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxyz[k] = -2.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxzz[k] = -4.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xyyy[k] = 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xyyz[k] = -2.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xyzz[k] = -4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xzzz[k] = -6.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yyyy[k] = 4.0 * to_xyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yy_yyyz[k] = -2.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yyzz[k] = -4.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yzzz[k] = -6.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_zzzz[k] = -8.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-255 components of targeted buffer : DG

        auto to_x_z_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 60);

        auto to_x_z_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 61);

        auto to_x_z_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 62);

        auto to_x_z_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 63);

        auto to_x_z_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 64);

        auto to_x_z_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 65);

        auto to_x_z_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 66);

        auto to_x_z_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 67);

        auto to_x_z_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 68);

        auto to_x_z_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 69);

        auto to_x_z_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 70);

        auto to_x_z_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 71);

        auto to_x_z_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 72);

        auto to_x_z_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 73);

        auto to_x_z_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_x_z_yz_xxxx,     \
                             to_x_z_yz_xxxy, \
                             to_x_z_yz_xxxz, \
                             to_x_z_yz_xxyy, \
                             to_x_z_yz_xxyz, \
                             to_x_z_yz_xxzz, \
                             to_x_z_yz_xyyy, \
                             to_x_z_yz_xyyz, \
                             to_x_z_yz_xyzz, \
                             to_x_z_yz_xzzz, \
                             to_x_z_yz_yyyy, \
                             to_x_z_yz_yyyz, \
                             to_x_z_yz_yyzz, \
                             to_x_z_yz_yzzz, \
                             to_x_z_yz_zzzz, \
                             to_xyz_xxx,     \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             to_xyz_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yz_xxxx[k] = 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxxy[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxxz[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxyy[k] = 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxyz[k] = -2.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxzz[k] = -4.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xyyy[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xyyz[k] = -2.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xyzz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xzzz[k] = -6.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yyyy[k] = 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yz_yyyz[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yyzz[k] = -4.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yzzz[k] = -6.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_zzzz[k] = -8.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 255-270 components of targeted buffer : DG

        auto to_x_z_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 75);

        auto to_x_z_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 76);

        auto to_x_z_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 77);

        auto to_x_z_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 78);

        auto to_x_z_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 79);

        auto to_x_z_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 80);

        auto to_x_z_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 81);

        auto to_x_z_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 82);

        auto to_x_z_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 83);

        auto to_x_z_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 84);

        auto to_x_z_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 85);

        auto to_x_z_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 86);

        auto to_x_z_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 87);

        auto to_x_z_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 88);

        auto to_x_z_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 2 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_x_z_zz_xxxx,     \
                             to_x_z_zz_xxxy, \
                             to_x_z_zz_xxxz, \
                             to_x_z_zz_xxyy, \
                             to_x_z_zz_xxyz, \
                             to_x_z_zz_xxzz, \
                             to_x_z_zz_xyyy, \
                             to_x_z_zz_xyyz, \
                             to_x_z_zz_xyzz, \
                             to_x_z_zz_xzzz, \
                             to_x_z_zz_yyyy, \
                             to_x_z_zz_yyyz, \
                             to_x_z_zz_yyzz, \
                             to_x_z_zz_yzzz, \
                             to_x_z_zz_zzzz, \
                             to_xzz_xxx,     \
                             to_xzz_xxxxz,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxxzz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xxzzz,   \
                             to_xzz_xyy,     \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_xzzzz,   \
                             to_xzz_yyy,     \
                             to_xzz_yyyyz,   \
                             to_xzz_yyyzz,   \
                             to_xzz_yyz,     \
                             to_xzz_yyzzz,   \
                             to_xzz_yzz,     \
                             to_xzz_yzzzz,   \
                             to_xzz_zzz,     \
                             to_xzz_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zz_xxxx[k] = 4.0 * to_xzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxxy[k] = 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxxz[k] = -2.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxyy[k] = 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxyz[k] = -2.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxzz[k] = -4.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xyyy[k] = 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xyyz[k] = -2.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xyzz[k] = -4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xzzz[k] = -6.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yyyy[k] = 4.0 * to_xzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_zz_yyyz[k] = -2.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yyzz[k] = -4.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yzzz[k] = -6.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_zzzz[k] = -8.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-285 components of targeted buffer : DG

        auto to_y_x_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 0);

        auto to_y_x_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 1);

        auto to_y_x_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 2);

        auto to_y_x_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 3);

        auto to_y_x_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 4);

        auto to_y_x_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 5);

        auto to_y_x_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 6);

        auto to_y_x_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 7);

        auto to_y_x_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 8);

        auto to_y_x_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 9);

        auto to_y_x_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 10);

        auto to_y_x_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 11);

        auto to_y_x_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 12);

        auto to_y_x_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 13);

        auto to_y_x_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxy_xxx,         \
                             to_xxy_xxxxx,   \
                             to_xxy_xxxxy,   \
                             to_xxy_xxxxz,   \
                             to_xxy_xxxyy,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxxzz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyy,   \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xxzzz,   \
                             to_xxy_xyy,     \
                             to_xxy_xyyyy,   \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_xzzzz,   \
                             to_xxy_yyy,     \
                             to_xxy_yyz,     \
                             to_xxy_yzz,     \
                             to_xxy_zzz,     \
                             to_y_x_xx_xxxx, \
                             to_y_x_xx_xxxy, \
                             to_y_x_xx_xxxz, \
                             to_y_x_xx_xxyy, \
                             to_y_x_xx_xxyz, \
                             to_y_x_xx_xxzz, \
                             to_y_x_xx_xyyy, \
                             to_y_x_xx_xyyz, \
                             to_y_x_xx_xyzz, \
                             to_y_x_xx_xzzz, \
                             to_y_x_xx_yyyy, \
                             to_y_x_xx_yyyz, \
                             to_y_x_xx_yyzz, \
                             to_y_x_xx_yzzz, \
                             to_y_x_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xx_xxxx[k] = -8.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xx_xxxy[k] = -6.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xx_xxxz[k] = -6.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xx_xxyy[k] = -4.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xx_xxyz[k] = -4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xx_xxzz[k] = -4.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xx_xyyy[k] = -2.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xx_xyyz[k] = -2.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xx_xyzz[k] = -2.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xx_xzzz[k] = -2.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xx_yyyy[k] = 4.0 * to_xxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xx_yyyz[k] = 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xx_yyzz[k] = 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xx_yzzz[k] = 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xx_zzzz[k] = 4.0 * to_xxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 285-300 components of targeted buffer : DG

        auto to_y_x_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 15);

        auto to_y_x_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 16);

        auto to_y_x_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 17);

        auto to_y_x_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 18);

        auto to_y_x_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 19);

        auto to_y_x_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 20);

        auto to_y_x_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 21);

        auto to_y_x_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 22);

        auto to_y_x_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 23);

        auto to_y_x_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 24);

        auto to_y_x_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 25);

        auto to_y_x_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 26);

        auto to_y_x_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 27);

        auto to_y_x_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 28);

        auto to_y_x_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxx,     \
                             to_x_xxxxy,     \
                             to_x_xxxxz,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyz,       \
                             to_x_yzz,       \
                             to_x_zzz,       \
                             to_xyy_xxx,     \
                             to_xyy_xxxxx,   \
                             to_xyy_xxxxy,   \
                             to_xyy_xxxxz,   \
                             to_xyy_xxxyy,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxxzz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyy,   \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xxzzz,   \
                             to_xyy_xyy,     \
                             to_xyy_xyyyy,   \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_xzzzz,   \
                             to_xyy_yyy,     \
                             to_xyy_yyz,     \
                             to_xyy_yzz,     \
                             to_xyy_zzz,     \
                             to_y_x_xy_xxxx, \
                             to_y_x_xy_xxxy, \
                             to_y_x_xy_xxxz, \
                             to_y_x_xy_xxyy, \
                             to_y_x_xy_xxyz, \
                             to_y_x_xy_xxzz, \
                             to_y_x_xy_xyyy, \
                             to_y_x_xy_xyyz, \
                             to_y_x_xy_xyzz, \
                             to_y_x_xy_xzzz, \
                             to_y_x_xy_yyyy, \
                             to_y_x_xy_yyyz, \
                             to_y_x_xy_yyzz, \
                             to_y_x_xy_yzzz, \
                             to_y_x_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xy_xxxx[k] = 4.0 * to_x_xxx[k] - 2.0 * to_x_xxxxx[k] * tke_0 - 8.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xy_xxxy[k] = 3.0 * to_x_xxy[k] - 2.0 * to_x_xxxxy[k] * tke_0 - 6.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xy_xxxz[k] = 3.0 * to_x_xxz[k] - 2.0 * to_x_xxxxz[k] * tke_0 - 6.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xy_xxyy[k] = 2.0 * to_x_xyy[k] - 2.0 * to_x_xxxyy[k] * tke_0 - 4.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xy_xxyz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xxxyz[k] * tke_0 - 4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xy_xxzz[k] = 2.0 * to_x_xzz[k] - 2.0 * to_x_xxxzz[k] * tke_0 - 4.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xy_xyyy[k] = to_x_yyy[k] - 2.0 * to_x_xxyyy[k] * tke_0 - 2.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xy_xyyz[k] = to_x_yyz[k] - 2.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xy_xyzz[k] = to_x_yzz[k] - 2.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xy_xzzz[k] = to_x_zzz[k] - 2.0 * to_x_xxzzz[k] * tke_0 - 2.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xy_yyyy[k] = -2.0 * to_x_xyyyy[k] * tke_0 + 4.0 * to_xyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xy_yyyz[k] = -2.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xy_yyzz[k] = -2.0 * to_x_xyyzz[k] * tke_0 + 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xy_yzzz[k] = -2.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xy_zzzz[k] = -2.0 * to_x_xzzzz[k] * tke_0 + 4.0 * to_xyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-315 components of targeted buffer : DG

        auto to_y_x_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 30);

        auto to_y_x_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 31);

        auto to_y_x_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 32);

        auto to_y_x_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 33);

        auto to_y_x_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 34);

        auto to_y_x_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 35);

        auto to_y_x_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 36);

        auto to_y_x_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 37);

        auto to_y_x_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 38);

        auto to_y_x_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 39);

        auto to_y_x_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 40);

        auto to_y_x_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 41);

        auto to_y_x_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 42);

        auto to_y_x_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 43);

        auto to_y_x_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxx,   \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyz,     \
                             to_xyz_yzz,     \
                             to_xyz_zzz,     \
                             to_y_x_xz_xxxx, \
                             to_y_x_xz_xxxy, \
                             to_y_x_xz_xxxz, \
                             to_y_x_xz_xxyy, \
                             to_y_x_xz_xxyz, \
                             to_y_x_xz_xxzz, \
                             to_y_x_xz_xyyy, \
                             to_y_x_xz_xyyz, \
                             to_y_x_xz_xyzz, \
                             to_y_x_xz_xzzz, \
                             to_y_x_xz_yyyy, \
                             to_y_x_xz_yyyz, \
                             to_y_x_xz_yyzz, \
                             to_y_x_xz_yzzz, \
                             to_y_x_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xz_xxxx[k] = -8.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xz_xxxy[k] = -6.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xz_xxxz[k] = -6.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xz_xxyy[k] = -4.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xz_xxyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xz_xxzz[k] = -4.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xz_xyyy[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xz_xyyz[k] = -2.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xz_xyzz[k] = -2.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xz_xzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xz_yyyy[k] = 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xz_yyyz[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xz_yyzz[k] = 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xz_yzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xz_zzzz[k] = 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 315-330 components of targeted buffer : DG

        auto to_y_x_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 45);

        auto to_y_x_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 46);

        auto to_y_x_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 47);

        auto to_y_x_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 48);

        auto to_y_x_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 49);

        auto to_y_x_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 50);

        auto to_y_x_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 51);

        auto to_y_x_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 52);

        auto to_y_x_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 53);

        auto to_y_x_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 54);

        auto to_y_x_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 55);

        auto to_y_x_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 56);

        auto to_y_x_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 57);

        auto to_y_x_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 58);

        auto to_y_x_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_y_x_yy_xxxx,     \
                             to_y_x_yy_xxxy, \
                             to_y_x_yy_xxxz, \
                             to_y_x_yy_xxyy, \
                             to_y_x_yy_xxyz, \
                             to_y_x_yy_xxzz, \
                             to_y_x_yy_xyyy, \
                             to_y_x_yy_xyyz, \
                             to_y_x_yy_xyzz, \
                             to_y_x_yy_xzzz, \
                             to_y_x_yy_yyyy, \
                             to_y_x_yy_yyyz, \
                             to_y_x_yy_yyzz, \
                             to_y_x_yy_yzzz, \
                             to_y_x_yy_zzzz, \
                             to_y_xxx,       \
                             to_y_xxxxx,     \
                             to_y_xxxxy,     \
                             to_y_xxxxz,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyz,       \
                             to_y_yzz,       \
                             to_y_zzz,       \
                             to_yyy_xxx,     \
                             to_yyy_xxxxx,   \
                             to_yyy_xxxxy,   \
                             to_yyy_xxxxz,   \
                             to_yyy_xxxyy,   \
                             to_yyy_xxxyz,   \
                             to_yyy_xxxzz,   \
                             to_yyy_xxy,     \
                             to_yyy_xxyyy,   \
                             to_yyy_xxyyz,   \
                             to_yyy_xxyzz,   \
                             to_yyy_xxz,     \
                             to_yyy_xxzzz,   \
                             to_yyy_xyy,     \
                             to_yyy_xyyyy,   \
                             to_yyy_xyyyz,   \
                             to_yyy_xyyzz,   \
                             to_yyy_xyz,     \
                             to_yyy_xyzzz,   \
                             to_yyy_xzz,     \
                             to_yyy_xzzzz,   \
                             to_yyy_yyy,     \
                             to_yyy_yyz,     \
                             to_yyy_yzz,     \
                             to_yyy_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yy_xxxx[k] = 8.0 * to_y_xxx[k] - 4.0 * to_y_xxxxx[k] * tke_0 - 8.0 * to_yyy_xxx[k] * tbe_0 + 4.0 * to_yyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yy_xxxy[k] = 6.0 * to_y_xxy[k] - 4.0 * to_y_xxxxy[k] * tke_0 - 6.0 * to_yyy_xxy[k] * tbe_0 + 4.0 * to_yyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yy_xxxz[k] = 6.0 * to_y_xxz[k] - 4.0 * to_y_xxxxz[k] * tke_0 - 6.0 * to_yyy_xxz[k] * tbe_0 + 4.0 * to_yyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yy_xxyy[k] = 4.0 * to_y_xyy[k] - 4.0 * to_y_xxxyy[k] * tke_0 - 4.0 * to_yyy_xyy[k] * tbe_0 + 4.0 * to_yyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yy_xxyz[k] = 4.0 * to_y_xyz[k] - 4.0 * to_y_xxxyz[k] * tke_0 - 4.0 * to_yyy_xyz[k] * tbe_0 + 4.0 * to_yyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yy_xxzz[k] = 4.0 * to_y_xzz[k] - 4.0 * to_y_xxxzz[k] * tke_0 - 4.0 * to_yyy_xzz[k] * tbe_0 + 4.0 * to_yyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yy_xyyy[k] = 2.0 * to_y_yyy[k] - 4.0 * to_y_xxyyy[k] * tke_0 - 2.0 * to_yyy_yyy[k] * tbe_0 + 4.0 * to_yyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yy_xyyz[k] = 2.0 * to_y_yyz[k] - 4.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_yyy_yyz[k] * tbe_0 + 4.0 * to_yyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yy_xyzz[k] = 2.0 * to_y_yzz[k] - 4.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_yyy_yzz[k] * tbe_0 + 4.0 * to_yyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yy_xzzz[k] = 2.0 * to_y_zzz[k] - 4.0 * to_y_xxzzz[k] * tke_0 - 2.0 * to_yyy_zzz[k] * tbe_0 + 4.0 * to_yyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yy_yyyy[k] = -4.0 * to_y_xyyyy[k] * tke_0 + 4.0 * to_yyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yy_yyyz[k] = -4.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_yyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yy_yyzz[k] = -4.0 * to_y_xyyzz[k] * tke_0 + 4.0 * to_yyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yy_yzzz[k] = -4.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_yyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yy_zzzz[k] = -4.0 * to_y_xzzzz[k] * tke_0 + 4.0 * to_yyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-345 components of targeted buffer : DG

        auto to_y_x_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 60);

        auto to_y_x_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 61);

        auto to_y_x_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 62);

        auto to_y_x_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 63);

        auto to_y_x_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 64);

        auto to_y_x_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 65);

        auto to_y_x_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 66);

        auto to_y_x_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 67);

        auto to_y_x_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 68);

        auto to_y_x_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 69);

        auto to_y_x_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 70);

        auto to_y_x_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 71);

        auto to_y_x_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 72);

        auto to_y_x_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 73);

        auto to_y_x_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_x_yz_xxxx,     \
                             to_y_x_yz_xxxy, \
                             to_y_x_yz_xxxz, \
                             to_y_x_yz_xxyy, \
                             to_y_x_yz_xxyz, \
                             to_y_x_yz_xxzz, \
                             to_y_x_yz_xyyy, \
                             to_y_x_yz_xyyz, \
                             to_y_x_yz_xyzz, \
                             to_y_x_yz_xzzz, \
                             to_y_x_yz_yyyy, \
                             to_y_x_yz_yyyz, \
                             to_y_x_yz_yyzz, \
                             to_y_x_yz_yzzz, \
                             to_y_x_yz_zzzz, \
                             to_yyz_xxx,     \
                             to_yyz_xxxxx,   \
                             to_yyz_xxxxy,   \
                             to_yyz_xxxxz,   \
                             to_yyz_xxxyy,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxxzz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyy,   \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xxzzz,   \
                             to_yyz_xyy,     \
                             to_yyz_xyyyy,   \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_xzzzz,   \
                             to_yyz_yyy,     \
                             to_yyz_yyz,     \
                             to_yyz_yzz,     \
                             to_yyz_zzz,     \
                             to_z_xxx,       \
                             to_z_xxxxx,     \
                             to_z_xxxxy,     \
                             to_z_xxxxz,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyz,       \
                             to_z_yzz,       \
                             to_z_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yz_xxxx[k] = 4.0 * to_z_xxx[k] - 2.0 * to_z_xxxxx[k] * tke_0 - 8.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yz_xxxy[k] = 3.0 * to_z_xxy[k] - 2.0 * to_z_xxxxy[k] * tke_0 - 6.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yz_xxxz[k] = 3.0 * to_z_xxz[k] - 2.0 * to_z_xxxxz[k] * tke_0 - 6.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yz_xxyy[k] = 2.0 * to_z_xyy[k] - 2.0 * to_z_xxxyy[k] * tke_0 - 4.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yz_xxyz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xxxyz[k] * tke_0 - 4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yz_xxzz[k] = 2.0 * to_z_xzz[k] - 2.0 * to_z_xxxzz[k] * tke_0 - 4.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yz_xyyy[k] = to_z_yyy[k] - 2.0 * to_z_xxyyy[k] * tke_0 - 2.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yz_xyyz[k] = to_z_yyz[k] - 2.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yz_xyzz[k] = to_z_yzz[k] - 2.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yz_xzzz[k] = to_z_zzz[k] - 2.0 * to_z_xxzzz[k] * tke_0 - 2.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yz_yyyy[k] = -2.0 * to_z_xyyyy[k] * tke_0 + 4.0 * to_yyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yz_yyyz[k] = -2.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yz_yyzz[k] = -2.0 * to_z_xyyzz[k] * tke_0 + 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yz_yzzz[k] = -2.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yz_zzzz[k] = -2.0 * to_z_xzzzz[k] * tke_0 + 4.0 * to_yyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 345-360 components of targeted buffer : DG

        auto to_y_x_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 75);

        auto to_y_x_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 76);

        auto to_y_x_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 77);

        auto to_y_x_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 78);

        auto to_y_x_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 79);

        auto to_y_x_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 80);

        auto to_y_x_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 81);

        auto to_y_x_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 82);

        auto to_y_x_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 83);

        auto to_y_x_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 84);

        auto to_y_x_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 85);

        auto to_y_x_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 86);

        auto to_y_x_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 87);

        auto to_y_x_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 88);

        auto to_y_x_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 3 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_y_x_zz_xxxx,     \
                             to_y_x_zz_xxxy, \
                             to_y_x_zz_xxxz, \
                             to_y_x_zz_xxyy, \
                             to_y_x_zz_xxyz, \
                             to_y_x_zz_xxzz, \
                             to_y_x_zz_xyyy, \
                             to_y_x_zz_xyyz, \
                             to_y_x_zz_xyzz, \
                             to_y_x_zz_xzzz, \
                             to_y_x_zz_yyyy, \
                             to_y_x_zz_yyyz, \
                             to_y_x_zz_yyzz, \
                             to_y_x_zz_yzzz, \
                             to_y_x_zz_zzzz, \
                             to_yzz_xxx,     \
                             to_yzz_xxxxx,   \
                             to_yzz_xxxxy,   \
                             to_yzz_xxxxz,   \
                             to_yzz_xxxyy,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxxzz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyy,   \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xxzzz,   \
                             to_yzz_xyy,     \
                             to_yzz_xyyyy,   \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_xzzzz,   \
                             to_yzz_yyy,     \
                             to_yzz_yyz,     \
                             to_yzz_yzz,     \
                             to_yzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zz_xxxx[k] = -8.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_zz_xxxy[k] = -6.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_zz_xxxz[k] = -6.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_zz_xxyy[k] = -4.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_zz_xxyz[k] = -4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_zz_xxzz[k] = -4.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_zz_xyyy[k] = -2.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_zz_xyyz[k] = -2.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_zz_xyzz[k] = -2.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_zz_xzzz[k] = -2.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_zz_yyyy[k] = 4.0 * to_yzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_zz_yyyz[k] = 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_zz_yyzz[k] = 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_zz_yzzz[k] = 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_zz_zzzz[k] = 4.0 * to_yzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-375 components of targeted buffer : DG

        auto to_y_y_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 0);

        auto to_y_y_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 1);

        auto to_y_y_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 2);

        auto to_y_y_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 3);

        auto to_y_y_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 4);

        auto to_y_y_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 5);

        auto to_y_y_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 6);

        auto to_y_y_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 7);

        auto to_y_y_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 8);

        auto to_y_y_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 9);

        auto to_y_y_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 10);

        auto to_y_y_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 11);

        auto to_y_y_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 12);

        auto to_y_y_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 13);

        auto to_y_y_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxy_xxx,         \
                             to_xxy_xxxxy,   \
                             to_xxy_xxxyy,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyy,   \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xyy,     \
                             to_xxy_xyyyy,   \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_yyy,     \
                             to_xxy_yyyyy,   \
                             to_xxy_yyyyz,   \
                             to_xxy_yyyzz,   \
                             to_xxy_yyz,     \
                             to_xxy_yyzzz,   \
                             to_xxy_yzz,     \
                             to_xxy_yzzzz,   \
                             to_xxy_zzz,     \
                             to_y_y_xx_xxxx, \
                             to_y_y_xx_xxxy, \
                             to_y_y_xx_xxxz, \
                             to_y_y_xx_xxyy, \
                             to_y_y_xx_xxyz, \
                             to_y_y_xx_xxzz, \
                             to_y_y_xx_xyyy, \
                             to_y_y_xx_xyyz, \
                             to_y_y_xx_xyzz, \
                             to_y_y_xx_xzzz, \
                             to_y_y_xx_yyyy, \
                             to_y_y_xx_yyyz, \
                             to_y_y_xx_yyzz, \
                             to_y_y_xx_yzzz, \
                             to_y_y_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xx_xxxx[k] = 4.0 * to_xxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xx_xxxy[k] = -2.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xxxz[k] = 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xx_xxyy[k] = -4.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xxyz[k] = -2.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xx_xxzz[k] = 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xx_xyyy[k] = -6.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xyyz[k] = -4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xx_xyzz[k] = -2.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xx_xzzz[k] = 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xx_yyyy[k] = -8.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xx_yyyz[k] = -6.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xx_yyzz[k] = -4.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xx_yzzz[k] = -2.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xx_zzzz[k] = 4.0 * to_xxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 375-390 components of targeted buffer : DG

        auto to_y_y_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 15);

        auto to_y_y_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 16);

        auto to_y_y_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 17);

        auto to_y_y_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 18);

        auto to_y_y_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 19);

        auto to_y_y_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 20);

        auto to_y_y_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 21);

        auto to_y_y_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 22);

        auto to_y_y_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 23);

        auto to_y_y_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 24);

        auto to_y_y_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 25);

        auto to_y_y_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 26);

        auto to_y_y_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 27);

        auto to_y_y_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 28);

        auto to_y_y_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxy,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_yyy,       \
                             to_x_yyyyy,     \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_zzz,       \
                             to_xyy_xxx,     \
                             to_xyy_xxxxy,   \
                             to_xyy_xxxyy,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyy,   \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xyy,     \
                             to_xyy_xyyyy,   \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_yyy,     \
                             to_xyy_yyyyy,   \
                             to_xyy_yyyyz,   \
                             to_xyy_yyyzz,   \
                             to_xyy_yyz,     \
                             to_xyy_yyzzz,   \
                             to_xyy_yzz,     \
                             to_xyy_yzzzz,   \
                             to_xyy_zzz,     \
                             to_y_y_xy_xxxx, \
                             to_y_y_xy_xxxy, \
                             to_y_y_xy_xxxz, \
                             to_y_y_xy_xxyy, \
                             to_y_y_xy_xxyz, \
                             to_y_y_xy_xxzz, \
                             to_y_y_xy_xyyy, \
                             to_y_y_xy_xyyz, \
                             to_y_y_xy_xyzz, \
                             to_y_y_xy_xzzz, \
                             to_y_y_xy_yyyy, \
                             to_y_y_xy_yyyz, \
                             to_y_y_xy_yyzz, \
                             to_y_y_xy_yzzz, \
                             to_y_y_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xy_xxxx[k] = -2.0 * to_x_xxxxy[k] * tke_0 + 4.0 * to_xyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xy_xxxy[k] = to_x_xxx[k] - 2.0 * to_x_xxxyy[k] * tke_0 - 2.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xxxz[k] = -2.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xy_xxyy[k] = 2.0 * to_x_xxy[k] - 2.0 * to_x_xxyyy[k] * tke_0 - 4.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xxyz[k] = to_x_xxz[k] - 2.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xy_xxzz[k] = -2.0 * to_x_xxyzz[k] * tke_0 + 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xy_xyyy[k] = 3.0 * to_x_xyy[k] - 2.0 * to_x_xyyyy[k] * tke_0 - 6.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xyyz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xyyyz[k] * tke_0 - 4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xy_xyzz[k] = to_x_xzz[k] - 2.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xy_xzzz[k] = -2.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xy_yyyy[k] = 4.0 * to_x_yyy[k] - 2.0 * to_x_yyyyy[k] * tke_0 - 8.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xy_yyyz[k] = 3.0 * to_x_yyz[k] - 2.0 * to_x_yyyyz[k] * tke_0 - 6.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xy_yyzz[k] = 2.0 * to_x_yzz[k] - 2.0 * to_x_yyyzz[k] * tke_0 - 4.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xy_yzzz[k] = to_x_zzz[k] - 2.0 * to_x_yyzzz[k] * tke_0 - 2.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xy_zzzz[k] = -2.0 * to_x_yzzzz[k] * tke_0 + 4.0 * to_xyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-405 components of targeted buffer : DG

        auto to_y_y_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 30);

        auto to_y_y_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 31);

        auto to_y_y_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 32);

        auto to_y_y_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 33);

        auto to_y_y_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 34);

        auto to_y_y_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 35);

        auto to_y_y_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 36);

        auto to_y_y_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 37);

        auto to_y_y_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 38);

        auto to_y_y_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 39);

        auto to_y_y_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 40);

        auto to_y_y_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 41);

        auto to_y_y_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 42);

        auto to_y_y_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 43);

        auto to_y_y_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_yyy,     \
                             to_xyz_yyyyy,   \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             to_y_y_xz_xxxx, \
                             to_y_y_xz_xxxy, \
                             to_y_y_xz_xxxz, \
                             to_y_y_xz_xxyy, \
                             to_y_y_xz_xxyz, \
                             to_y_y_xz_xxzz, \
                             to_y_y_xz_xyyy, \
                             to_y_y_xz_xyyz, \
                             to_y_y_xz_xyzz, \
                             to_y_y_xz_xzzz, \
                             to_y_y_xz_yyyy, \
                             to_y_y_xz_yyyz, \
                             to_y_y_xz_yyzz, \
                             to_y_y_xz_yzzz, \
                             to_y_y_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xz_xxxx[k] = 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xz_xxxy[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xxxz[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xz_xxyy[k] = -4.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xxyz[k] = -2.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xz_xxzz[k] = 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xz_xyyy[k] = -6.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xyyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xz_xyzz[k] = -2.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xz_xzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xz_yyyy[k] = -8.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xz_yyyz[k] = -6.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xz_yyzz[k] = -4.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xz_yzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xz_zzzz[k] = 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 405-420 components of targeted buffer : DG

        auto to_y_y_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 45);

        auto to_y_y_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 46);

        auto to_y_y_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 47);

        auto to_y_y_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 48);

        auto to_y_y_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 49);

        auto to_y_y_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 50);

        auto to_y_y_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 51);

        auto to_y_y_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 52);

        auto to_y_y_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 53);

        auto to_y_y_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 54);

        auto to_y_y_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 55);

        auto to_y_y_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 56);

        auto to_y_y_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 57);

        auto to_y_y_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 58);

        auto to_y_y_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_y_xxx,           \
                             to_y_xxxxy,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_y_yy_xxxx, \
                             to_y_y_yy_xxxy, \
                             to_y_y_yy_xxxz, \
                             to_y_y_yy_xxyy, \
                             to_y_y_yy_xxyz, \
                             to_y_y_yy_xxzz, \
                             to_y_y_yy_xyyy, \
                             to_y_y_yy_xyyz, \
                             to_y_y_yy_xyzz, \
                             to_y_y_yy_xzzz, \
                             to_y_y_yy_yyyy, \
                             to_y_y_yy_yyyz, \
                             to_y_y_yy_yyzz, \
                             to_y_y_yy_yzzz, \
                             to_y_y_yy_zzzz, \
                             to_y_yyy,       \
                             to_y_yyyyy,     \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_zzz,       \
                             to_yyy_xxx,     \
                             to_yyy_xxxxy,   \
                             to_yyy_xxxyy,   \
                             to_yyy_xxxyz,   \
                             to_yyy_xxy,     \
                             to_yyy_xxyyy,   \
                             to_yyy_xxyyz,   \
                             to_yyy_xxyzz,   \
                             to_yyy_xxz,     \
                             to_yyy_xyy,     \
                             to_yyy_xyyyy,   \
                             to_yyy_xyyyz,   \
                             to_yyy_xyyzz,   \
                             to_yyy_xyz,     \
                             to_yyy_xyzzz,   \
                             to_yyy_xzz,     \
                             to_yyy_yyy,     \
                             to_yyy_yyyyy,   \
                             to_yyy_yyyyz,   \
                             to_yyy_yyyzz,   \
                             to_yyy_yyz,     \
                             to_yyy_yyzzz,   \
                             to_yyy_yzz,     \
                             to_yyy_yzzzz,   \
                             to_yyy_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yy_xxxx[k] = -4.0 * to_y_xxxxy[k] * tke_0 + 4.0 * to_yyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yy_xxxy[k] = 2.0 * to_y_xxx[k] - 4.0 * to_y_xxxyy[k] * tke_0 - 2.0 * to_yyy_xxx[k] * tbe_0 + 4.0 * to_yyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xxxz[k] = -4.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_yyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yy_xxyy[k] = 4.0 * to_y_xxy[k] - 4.0 * to_y_xxyyy[k] * tke_0 - 4.0 * to_yyy_xxy[k] * tbe_0 + 4.0 * to_yyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xxyz[k] = 2.0 * to_y_xxz[k] - 4.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_yyy_xxz[k] * tbe_0 + 4.0 * to_yyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yy_xxzz[k] = -4.0 * to_y_xxyzz[k] * tke_0 + 4.0 * to_yyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yy_xyyy[k] = 6.0 * to_y_xyy[k] - 4.0 * to_y_xyyyy[k] * tke_0 - 6.0 * to_yyy_xyy[k] * tbe_0 + 4.0 * to_yyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xyyz[k] = 4.0 * to_y_xyz[k] - 4.0 * to_y_xyyyz[k] * tke_0 - 4.0 * to_yyy_xyz[k] * tbe_0 + 4.0 * to_yyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yy_xyzz[k] = 2.0 * to_y_xzz[k] - 4.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_yyy_xzz[k] * tbe_0 + 4.0 * to_yyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yy_xzzz[k] = -4.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_yyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yy_yyyy[k] = 8.0 * to_y_yyy[k] - 4.0 * to_y_yyyyy[k] * tke_0 - 8.0 * to_yyy_yyy[k] * tbe_0 + 4.0 * to_yyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yy_yyyz[k] = 6.0 * to_y_yyz[k] - 4.0 * to_y_yyyyz[k] * tke_0 - 6.0 * to_yyy_yyz[k] * tbe_0 + 4.0 * to_yyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yy_yyzz[k] = 4.0 * to_y_yzz[k] - 4.0 * to_y_yyyzz[k] * tke_0 - 4.0 * to_yyy_yzz[k] * tbe_0 + 4.0 * to_yyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yy_yzzz[k] = 2.0 * to_y_zzz[k] - 4.0 * to_y_yyzzz[k] * tke_0 - 2.0 * to_yyy_zzz[k] * tbe_0 + 4.0 * to_yyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yy_zzzz[k] = -4.0 * to_y_yzzzz[k] * tke_0 + 4.0 * to_yyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-435 components of targeted buffer : DG

        auto to_y_y_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 60);

        auto to_y_y_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 61);

        auto to_y_y_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 62);

        auto to_y_y_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 63);

        auto to_y_y_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 64);

        auto to_y_y_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 65);

        auto to_y_y_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 66);

        auto to_y_y_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 67);

        auto to_y_y_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 68);

        auto to_y_y_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 69);

        auto to_y_y_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 70);

        auto to_y_y_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 71);

        auto to_y_y_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 72);

        auto to_y_y_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 73);

        auto to_y_y_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_y_yz_xxxx,     \
                             to_y_y_yz_xxxy, \
                             to_y_y_yz_xxxz, \
                             to_y_y_yz_xxyy, \
                             to_y_y_yz_xxyz, \
                             to_y_y_yz_xxzz, \
                             to_y_y_yz_xyyy, \
                             to_y_y_yz_xyyz, \
                             to_y_y_yz_xyzz, \
                             to_y_y_yz_xzzz, \
                             to_y_y_yz_yyyy, \
                             to_y_y_yz_yyyz, \
                             to_y_y_yz_yyzz, \
                             to_y_y_yz_yzzz, \
                             to_y_y_yz_zzzz, \
                             to_yyz_xxx,     \
                             to_yyz_xxxxy,   \
                             to_yyz_xxxyy,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyy,   \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xyy,     \
                             to_yyz_xyyyy,   \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_yyy,     \
                             to_yyz_yyyyy,   \
                             to_yyz_yyyyz,   \
                             to_yyz_yyyzz,   \
                             to_yyz_yyz,     \
                             to_yyz_yyzzz,   \
                             to_yyz_yzz,     \
                             to_yyz_yzzzz,   \
                             to_yyz_zzz,     \
                             to_z_xxx,       \
                             to_z_xxxxy,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_yyy,       \
                             to_z_yyyyy,     \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_zzz,       \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yz_xxxx[k] = -2.0 * to_z_xxxxy[k] * tke_0 + 4.0 * to_yyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yz_xxxy[k] = to_z_xxx[k] - 2.0 * to_z_xxxyy[k] * tke_0 - 2.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xxxz[k] = -2.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yz_xxyy[k] = 2.0 * to_z_xxy[k] - 2.0 * to_z_xxyyy[k] * tke_0 - 4.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xxyz[k] = to_z_xxz[k] - 2.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yz_xxzz[k] = -2.0 * to_z_xxyzz[k] * tke_0 + 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yz_xyyy[k] = 3.0 * to_z_xyy[k] - 2.0 * to_z_xyyyy[k] * tke_0 - 6.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xyyz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xyyyz[k] * tke_0 - 4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yz_xyzz[k] = to_z_xzz[k] - 2.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yz_xzzz[k] = -2.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yz_yyyy[k] = 4.0 * to_z_yyy[k] - 2.0 * to_z_yyyyy[k] * tke_0 - 8.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yz_yyyz[k] = 3.0 * to_z_yyz[k] - 2.0 * to_z_yyyyz[k] * tke_0 - 6.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yz_yyzz[k] = 2.0 * to_z_yzz[k] - 2.0 * to_z_yyyzz[k] * tke_0 - 4.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yz_yzzz[k] = to_z_zzz[k] - 2.0 * to_z_yyzzz[k] * tke_0 - 2.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yz_zzzz[k] = -2.0 * to_z_yzzzz[k] * tke_0 + 4.0 * to_yyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 435-450 components of targeted buffer : DG

        auto to_y_y_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 75);

        auto to_y_y_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 76);

        auto to_y_y_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 77);

        auto to_y_y_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 78);

        auto to_y_y_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 79);

        auto to_y_y_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 80);

        auto to_y_y_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 81);

        auto to_y_y_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 82);

        auto to_y_y_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 83);

        auto to_y_y_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 84);

        auto to_y_y_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 85);

        auto to_y_y_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 86);

        auto to_y_y_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 87);

        auto to_y_y_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 88);

        auto to_y_y_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 4 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_y_y_zz_xxxx,     \
                             to_y_y_zz_xxxy, \
                             to_y_y_zz_xxxz, \
                             to_y_y_zz_xxyy, \
                             to_y_y_zz_xxyz, \
                             to_y_y_zz_xxzz, \
                             to_y_y_zz_xyyy, \
                             to_y_y_zz_xyyz, \
                             to_y_y_zz_xyzz, \
                             to_y_y_zz_xzzz, \
                             to_y_y_zz_yyyy, \
                             to_y_y_zz_yyyz, \
                             to_y_y_zz_yyzz, \
                             to_y_y_zz_yzzz, \
                             to_y_y_zz_zzzz, \
                             to_yzz_xxx,     \
                             to_yzz_xxxxy,   \
                             to_yzz_xxxyy,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyy,   \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xyy,     \
                             to_yzz_xyyyy,   \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_yyy,     \
                             to_yzz_yyyyy,   \
                             to_yzz_yyyyz,   \
                             to_yzz_yyyzz,   \
                             to_yzz_yyz,     \
                             to_yzz_yyzzz,   \
                             to_yzz_yzz,     \
                             to_yzz_yzzzz,   \
                             to_yzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zz_xxxx[k] = 4.0 * to_yzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_zz_xxxy[k] = -2.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xxxz[k] = 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_zz_xxyy[k] = -4.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xxyz[k] = -2.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_zz_xxzz[k] = 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_zz_xyyy[k] = -6.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xyyz[k] = -4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_zz_xyzz[k] = -2.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_zz_xzzz[k] = 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_zz_yyyy[k] = -8.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_zz_yyyz[k] = -6.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_zz_yyzz[k] = -4.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_zz_yzzz[k] = -2.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_zz_zzzz[k] = 4.0 * to_yzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-465 components of targeted buffer : DG

        auto to_y_z_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 0);

        auto to_y_z_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 1);

        auto to_y_z_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 2);

        auto to_y_z_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 3);

        auto to_y_z_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 4);

        auto to_y_z_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 5);

        auto to_y_z_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 6);

        auto to_y_z_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 7);

        auto to_y_z_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 8);

        auto to_y_z_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 9);

        auto to_y_z_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 10);

        auto to_y_z_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 11);

        auto to_y_z_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 12);

        auto to_y_z_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 13);

        auto to_y_z_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxy_xxx,         \
                             to_xxy_xxxxz,   \
                             to_xxy_xxxyz,   \
                             to_xxy_xxxzz,   \
                             to_xxy_xxy,     \
                             to_xxy_xxyyz,   \
                             to_xxy_xxyzz,   \
                             to_xxy_xxz,     \
                             to_xxy_xxzzz,   \
                             to_xxy_xyy,     \
                             to_xxy_xyyyz,   \
                             to_xxy_xyyzz,   \
                             to_xxy_xyz,     \
                             to_xxy_xyzzz,   \
                             to_xxy_xzz,     \
                             to_xxy_xzzzz,   \
                             to_xxy_yyy,     \
                             to_xxy_yyyyz,   \
                             to_xxy_yyyzz,   \
                             to_xxy_yyz,     \
                             to_xxy_yyzzz,   \
                             to_xxy_yzz,     \
                             to_xxy_yzzzz,   \
                             to_xxy_zzz,     \
                             to_xxy_zzzzz,   \
                             to_y_z_xx_xxxx, \
                             to_y_z_xx_xxxy, \
                             to_y_z_xx_xxxz, \
                             to_y_z_xx_xxyy, \
                             to_y_z_xx_xxyz, \
                             to_y_z_xx_xxzz, \
                             to_y_z_xx_xyyy, \
                             to_y_z_xx_xyyz, \
                             to_y_z_xx_xyzz, \
                             to_y_z_xx_xzzz, \
                             to_y_z_xx_yyyy, \
                             to_y_z_xx_yyyz, \
                             to_y_z_xx_yyzz, \
                             to_y_z_xx_yzzz, \
                             to_y_z_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xx_xxxx[k] = 4.0 * to_xxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxxy[k] = 4.0 * to_xxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxxz[k] = -2.0 * to_xxy_xxx[k] * tbe_0 + 4.0 * to_xxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxyy[k] = 4.0 * to_xxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxyz[k] = -2.0 * to_xxy_xxy[k] * tbe_0 + 4.0 * to_xxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxzz[k] = -4.0 * to_xxy_xxz[k] * tbe_0 + 4.0 * to_xxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xyyy[k] = 4.0 * to_xxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xyyz[k] = -2.0 * to_xxy_xyy[k] * tbe_0 + 4.0 * to_xxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xyzz[k] = -4.0 * to_xxy_xyz[k] * tbe_0 + 4.0 * to_xxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xzzz[k] = -6.0 * to_xxy_xzz[k] * tbe_0 + 4.0 * to_xxy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yyyy[k] = 4.0 * to_xxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xx_yyyz[k] = -2.0 * to_xxy_yyy[k] * tbe_0 + 4.0 * to_xxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yyzz[k] = -4.0 * to_xxy_yyz[k] * tbe_0 + 4.0 * to_xxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yzzz[k] = -6.0 * to_xxy_yzz[k] * tbe_0 + 4.0 * to_xxy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_zzzz[k] = -8.0 * to_xxy_zzz[k] * tbe_0 + 4.0 * to_xxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 465-480 components of targeted buffer : DG

        auto to_y_z_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 15);

        auto to_y_z_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 16);

        auto to_y_z_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 17);

        auto to_y_z_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 18);

        auto to_y_z_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 19);

        auto to_y_z_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 20);

        auto to_y_z_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 21);

        auto to_y_z_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 22);

        auto to_y_z_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 23);

        auto to_y_z_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 24);

        auto to_y_z_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 25);

        auto to_y_z_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 26);

        auto to_y_z_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 27);

        auto to_y_z_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 28);

        auto to_y_z_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxz,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_zzz,       \
                             to_x_zzzzz,     \
                             to_xyy_xxx,     \
                             to_xyy_xxxxz,   \
                             to_xyy_xxxyz,   \
                             to_xyy_xxxzz,   \
                             to_xyy_xxy,     \
                             to_xyy_xxyyz,   \
                             to_xyy_xxyzz,   \
                             to_xyy_xxz,     \
                             to_xyy_xxzzz,   \
                             to_xyy_xyy,     \
                             to_xyy_xyyyz,   \
                             to_xyy_xyyzz,   \
                             to_xyy_xyz,     \
                             to_xyy_xyzzz,   \
                             to_xyy_xzz,     \
                             to_xyy_xzzzz,   \
                             to_xyy_yyy,     \
                             to_xyy_yyyyz,   \
                             to_xyy_yyyzz,   \
                             to_xyy_yyz,     \
                             to_xyy_yyzzz,   \
                             to_xyy_yzz,     \
                             to_xyy_yzzzz,   \
                             to_xyy_zzz,     \
                             to_xyy_zzzzz,   \
                             to_y_z_xy_xxxx, \
                             to_y_z_xy_xxxy, \
                             to_y_z_xy_xxxz, \
                             to_y_z_xy_xxyy, \
                             to_y_z_xy_xxyz, \
                             to_y_z_xy_xxzz, \
                             to_y_z_xy_xyyy, \
                             to_y_z_xy_xyyz, \
                             to_y_z_xy_xyzz, \
                             to_y_z_xy_xzzz, \
                             to_y_z_xy_yyyy, \
                             to_y_z_xy_yyyz, \
                             to_y_z_xy_yyzz, \
                             to_y_z_xy_yzzz, \
                             to_y_z_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xy_xxxx[k] = -2.0 * to_x_xxxxz[k] * tke_0 + 4.0 * to_xyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxxy[k] = -2.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxxz[k] = to_x_xxx[k] - 2.0 * to_x_xxxzz[k] * tke_0 - 2.0 * to_xyy_xxx[k] * tbe_0 + 4.0 * to_xyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxyy[k] = -2.0 * to_x_xxyyz[k] * tke_0 + 4.0 * to_xyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxyz[k] = to_x_xxy[k] - 2.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xyy_xxy[k] * tbe_0 + 4.0 * to_xyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxzz[k] = 2.0 * to_x_xxz[k] - 2.0 * to_x_xxzzz[k] * tke_0 - 4.0 * to_xyy_xxz[k] * tbe_0 + 4.0 * to_xyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xyyy[k] = -2.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xyyz[k] = to_x_xyy[k] - 2.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xyy_xyy[k] * tbe_0 + 4.0 * to_xyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xyzz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xyzzz[k] * tke_0 - 4.0 * to_xyy_xyz[k] * tbe_0 + 4.0 * to_xyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xzzz[k] = 3.0 * to_x_xzz[k] - 2.0 * to_x_xzzzz[k] * tke_0 - 6.0 * to_xyy_xzz[k] * tbe_0 + 4.0 * to_xyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yyyy[k] = -2.0 * to_x_yyyyz[k] * tke_0 + 4.0 * to_xyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xy_yyyz[k] = to_x_yyy[k] - 2.0 * to_x_yyyzz[k] * tke_0 - 2.0 * to_xyy_yyy[k] * tbe_0 + 4.0 * to_xyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yyzz[k] = 2.0 * to_x_yyz[k] - 2.0 * to_x_yyzzz[k] * tke_0 - 4.0 * to_xyy_yyz[k] * tbe_0 + 4.0 * to_xyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yzzz[k] = 3.0 * to_x_yzz[k] - 2.0 * to_x_yzzzz[k] * tke_0 - 6.0 * to_xyy_yzz[k] * tbe_0 + 4.0 * to_xyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_zzzz[k] = 4.0 * to_x_zzz[k] - 2.0 * to_x_zzzzz[k] * tke_0 - 8.0 * to_xyy_zzz[k] * tbe_0 + 4.0 * to_xyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-495 components of targeted buffer : DG

        auto to_y_z_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 30);

        auto to_y_z_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 31);

        auto to_y_z_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 32);

        auto to_y_z_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 33);

        auto to_y_z_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 34);

        auto to_y_z_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 35);

        auto to_y_z_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 36);

        auto to_y_z_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 37);

        auto to_y_z_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 38);

        auto to_y_z_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 39);

        auto to_y_z_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 40);

        auto to_y_z_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 41);

        auto to_y_z_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 42);

        auto to_y_z_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 43);

        auto to_y_z_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             to_xyz_zzzzz,   \
                             to_y_z_xz_xxxx, \
                             to_y_z_xz_xxxy, \
                             to_y_z_xz_xxxz, \
                             to_y_z_xz_xxyy, \
                             to_y_z_xz_xxyz, \
                             to_y_z_xz_xxzz, \
                             to_y_z_xz_xyyy, \
                             to_y_z_xz_xyyz, \
                             to_y_z_xz_xyzz, \
                             to_y_z_xz_xzzz, \
                             to_y_z_xz_yyyy, \
                             to_y_z_xz_yyyz, \
                             to_y_z_xz_yyzz, \
                             to_y_z_xz_yzzz, \
                             to_y_z_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xz_xxxx[k] = 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxxy[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxxz[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxyy[k] = 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxyz[k] = -2.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxzz[k] = -4.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xyyy[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xyyz[k] = -2.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xyzz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xzzz[k] = -6.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yyyy[k] = 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xz_yyyz[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yyzz[k] = -4.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yzzz[k] = -6.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_zzzz[k] = -8.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 495-510 components of targeted buffer : DG

        auto to_y_z_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 45);

        auto to_y_z_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 46);

        auto to_y_z_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 47);

        auto to_y_z_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 48);

        auto to_y_z_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 49);

        auto to_y_z_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 50);

        auto to_y_z_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 51);

        auto to_y_z_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 52);

        auto to_y_z_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 53);

        auto to_y_z_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 54);

        auto to_y_z_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 55);

        auto to_y_z_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 56);

        auto to_y_z_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 57);

        auto to_y_z_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 58);

        auto to_y_z_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_y_xxx,           \
                             to_y_xxxxz,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_z_yy_xxxx, \
                             to_y_z_yy_xxxy, \
                             to_y_z_yy_xxxz, \
                             to_y_z_yy_xxyy, \
                             to_y_z_yy_xxyz, \
                             to_y_z_yy_xxzz, \
                             to_y_z_yy_xyyy, \
                             to_y_z_yy_xyyz, \
                             to_y_z_yy_xyzz, \
                             to_y_z_yy_xzzz, \
                             to_y_z_yy_yyyy, \
                             to_y_z_yy_yyyz, \
                             to_y_z_yy_yyzz, \
                             to_y_z_yy_yzzz, \
                             to_y_z_yy_zzzz, \
                             to_y_zzz,       \
                             to_y_zzzzz,     \
                             to_yyy_xxx,     \
                             to_yyy_xxxxz,   \
                             to_yyy_xxxyz,   \
                             to_yyy_xxxzz,   \
                             to_yyy_xxy,     \
                             to_yyy_xxyyz,   \
                             to_yyy_xxyzz,   \
                             to_yyy_xxz,     \
                             to_yyy_xxzzz,   \
                             to_yyy_xyy,     \
                             to_yyy_xyyyz,   \
                             to_yyy_xyyzz,   \
                             to_yyy_xyz,     \
                             to_yyy_xyzzz,   \
                             to_yyy_xzz,     \
                             to_yyy_xzzzz,   \
                             to_yyy_yyy,     \
                             to_yyy_yyyyz,   \
                             to_yyy_yyyzz,   \
                             to_yyy_yyz,     \
                             to_yyy_yyzzz,   \
                             to_yyy_yzz,     \
                             to_yyy_yzzzz,   \
                             to_yyy_zzz,     \
                             to_yyy_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yy_xxxx[k] = -4.0 * to_y_xxxxz[k] * tke_0 + 4.0 * to_yyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxxy[k] = -4.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_yyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxxz[k] = 2.0 * to_y_xxx[k] - 4.0 * to_y_xxxzz[k] * tke_0 - 2.0 * to_yyy_xxx[k] * tbe_0 + 4.0 * to_yyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxyy[k] = -4.0 * to_y_xxyyz[k] * tke_0 + 4.0 * to_yyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxyz[k] = 2.0 * to_y_xxy[k] - 4.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_yyy_xxy[k] * tbe_0 + 4.0 * to_yyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxzz[k] = 4.0 * to_y_xxz[k] - 4.0 * to_y_xxzzz[k] * tke_0 - 4.0 * to_yyy_xxz[k] * tbe_0 + 4.0 * to_yyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xyyy[k] = -4.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_yyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xyyz[k] = 2.0 * to_y_xyy[k] - 4.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_yyy_xyy[k] * tbe_0 + 4.0 * to_yyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xyzz[k] = 4.0 * to_y_xyz[k] - 4.0 * to_y_xyzzz[k] * tke_0 - 4.0 * to_yyy_xyz[k] * tbe_0 + 4.0 * to_yyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xzzz[k] = 6.0 * to_y_xzz[k] - 4.0 * to_y_xzzzz[k] * tke_0 - 6.0 * to_yyy_xzz[k] * tbe_0 + 4.0 * to_yyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yyyy[k] = -4.0 * to_y_yyyyz[k] * tke_0 + 4.0 * to_yyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yy_yyyz[k] = 2.0 * to_y_yyy[k] - 4.0 * to_y_yyyzz[k] * tke_0 - 2.0 * to_yyy_yyy[k] * tbe_0 + 4.0 * to_yyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yyzz[k] = 4.0 * to_y_yyz[k] - 4.0 * to_y_yyzzz[k] * tke_0 - 4.0 * to_yyy_yyz[k] * tbe_0 + 4.0 * to_yyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yzzz[k] = 6.0 * to_y_yzz[k] - 4.0 * to_y_yzzzz[k] * tke_0 - 6.0 * to_yyy_yzz[k] * tbe_0 + 4.0 * to_yyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_zzzz[k] = 8.0 * to_y_zzz[k] - 4.0 * to_y_zzzzz[k] * tke_0 - 8.0 * to_yyy_zzz[k] * tbe_0 + 4.0 * to_yyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-525 components of targeted buffer : DG

        auto to_y_z_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 60);

        auto to_y_z_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 61);

        auto to_y_z_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 62);

        auto to_y_z_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 63);

        auto to_y_z_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 64);

        auto to_y_z_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 65);

        auto to_y_z_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 66);

        auto to_y_z_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 67);

        auto to_y_z_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 68);

        auto to_y_z_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 69);

        auto to_y_z_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 70);

        auto to_y_z_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 71);

        auto to_y_z_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 72);

        auto to_y_z_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 73);

        auto to_y_z_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_z_yz_xxxx,     \
                             to_y_z_yz_xxxy, \
                             to_y_z_yz_xxxz, \
                             to_y_z_yz_xxyy, \
                             to_y_z_yz_xxyz, \
                             to_y_z_yz_xxzz, \
                             to_y_z_yz_xyyy, \
                             to_y_z_yz_xyyz, \
                             to_y_z_yz_xyzz, \
                             to_y_z_yz_xzzz, \
                             to_y_z_yz_yyyy, \
                             to_y_z_yz_yyyz, \
                             to_y_z_yz_yyzz, \
                             to_y_z_yz_yzzz, \
                             to_y_z_yz_zzzz, \
                             to_yyz_xxx,     \
                             to_yyz_xxxxz,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxxzz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xxzzz,   \
                             to_yyz_xyy,     \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_xzzzz,   \
                             to_yyz_yyy,     \
                             to_yyz_yyyyz,   \
                             to_yyz_yyyzz,   \
                             to_yyz_yyz,     \
                             to_yyz_yyzzz,   \
                             to_yyz_yzz,     \
                             to_yyz_yzzzz,   \
                             to_yyz_zzz,     \
                             to_yyz_zzzzz,   \
                             to_z_xxx,       \
                             to_z_xxxxz,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_zzz,       \
                             to_z_zzzzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yz_xxxx[k] = -2.0 * to_z_xxxxz[k] * tke_0 + 4.0 * to_yyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxxy[k] = -2.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxxz[k] = to_z_xxx[k] - 2.0 * to_z_xxxzz[k] * tke_0 - 2.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxyy[k] = -2.0 * to_z_xxyyz[k] * tke_0 + 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxyz[k] = to_z_xxy[k] - 2.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxzz[k] = 2.0 * to_z_xxz[k] - 2.0 * to_z_xxzzz[k] * tke_0 - 4.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xyyy[k] = -2.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xyyz[k] = to_z_xyy[k] - 2.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xyzz[k] = 2.0 * to_z_xyz[k] - 2.0 * to_z_xyzzz[k] * tke_0 - 4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xzzz[k] = 3.0 * to_z_xzz[k] - 2.0 * to_z_xzzzz[k] * tke_0 - 6.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yyyy[k] = -2.0 * to_z_yyyyz[k] * tke_0 + 4.0 * to_yyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yz_yyyz[k] = to_z_yyy[k] - 2.0 * to_z_yyyzz[k] * tke_0 - 2.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yyzz[k] = 2.0 * to_z_yyz[k] - 2.0 * to_z_yyzzz[k] * tke_0 - 4.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yzzz[k] = 3.0 * to_z_yzz[k] - 2.0 * to_z_yzzzz[k] * tke_0 - 6.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_zzzz[k] = 4.0 * to_z_zzz[k] - 2.0 * to_z_zzzzz[k] * tke_0 - 8.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 525-540 components of targeted buffer : DG

        auto to_y_z_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 75);

        auto to_y_z_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 76);

        auto to_y_z_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 77);

        auto to_y_z_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 78);

        auto to_y_z_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 79);

        auto to_y_z_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 80);

        auto to_y_z_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 81);

        auto to_y_z_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 82);

        auto to_y_z_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 83);

        auto to_y_z_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 84);

        auto to_y_z_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 85);

        auto to_y_z_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 86);

        auto to_y_z_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 87);

        auto to_y_z_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 88);

        auto to_y_z_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 5 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_y_z_zz_xxxx,     \
                             to_y_z_zz_xxxy, \
                             to_y_z_zz_xxxz, \
                             to_y_z_zz_xxyy, \
                             to_y_z_zz_xxyz, \
                             to_y_z_zz_xxzz, \
                             to_y_z_zz_xyyy, \
                             to_y_z_zz_xyyz, \
                             to_y_z_zz_xyzz, \
                             to_y_z_zz_xzzz, \
                             to_y_z_zz_yyyy, \
                             to_y_z_zz_yyyz, \
                             to_y_z_zz_yyzz, \
                             to_y_z_zz_yzzz, \
                             to_y_z_zz_zzzz, \
                             to_yzz_xxx,     \
                             to_yzz_xxxxz,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxxzz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xxzzz,   \
                             to_yzz_xyy,     \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_xzzzz,   \
                             to_yzz_yyy,     \
                             to_yzz_yyyyz,   \
                             to_yzz_yyyzz,   \
                             to_yzz_yyz,     \
                             to_yzz_yyzzz,   \
                             to_yzz_yzz,     \
                             to_yzz_yzzzz,   \
                             to_yzz_zzz,     \
                             to_yzz_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zz_xxxx[k] = 4.0 * to_yzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxxy[k] = 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxxz[k] = -2.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxyy[k] = 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxyz[k] = -2.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxzz[k] = -4.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xyyy[k] = 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xyyz[k] = -2.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xyzz[k] = -4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xzzz[k] = -6.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yyyy[k] = 4.0 * to_yzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_zz_yyyz[k] = -2.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yyzz[k] = -4.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yzzz[k] = -6.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_zzzz[k] = -8.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-555 components of targeted buffer : DG

        auto to_z_x_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 0);

        auto to_z_x_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 1);

        auto to_z_x_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 2);

        auto to_z_x_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 3);

        auto to_z_x_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 4);

        auto to_z_x_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 5);

        auto to_z_x_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 6);

        auto to_z_x_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 7);

        auto to_z_x_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 8);

        auto to_z_x_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 9);

        auto to_z_x_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 10);

        auto to_z_x_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 11);

        auto to_z_x_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 12);

        auto to_z_x_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 13);

        auto to_z_x_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxz_xxx,         \
                             to_xxz_xxxxx,   \
                             to_xxz_xxxxy,   \
                             to_xxz_xxxxz,   \
                             to_xxz_xxxyy,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxxzz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyy,   \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xxzzz,   \
                             to_xxz_xyy,     \
                             to_xxz_xyyyy,   \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_xzzzz,   \
                             to_xxz_yyy,     \
                             to_xxz_yyz,     \
                             to_xxz_yzz,     \
                             to_xxz_zzz,     \
                             to_z_x_xx_xxxx, \
                             to_z_x_xx_xxxy, \
                             to_z_x_xx_xxxz, \
                             to_z_x_xx_xxyy, \
                             to_z_x_xx_xxyz, \
                             to_z_x_xx_xxzz, \
                             to_z_x_xx_xyyy, \
                             to_z_x_xx_xyyz, \
                             to_z_x_xx_xyzz, \
                             to_z_x_xx_xzzz, \
                             to_z_x_xx_yyyy, \
                             to_z_x_xx_yyyz, \
                             to_z_x_xx_yyzz, \
                             to_z_x_xx_yzzz, \
                             to_z_x_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xx_xxxx[k] = -8.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xx_xxxy[k] = -6.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xx_xxxz[k] = -6.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xx_xxyy[k] = -4.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xx_xxyz[k] = -4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xx_xxzz[k] = -4.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xx_xyyy[k] = -2.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xx_xyyz[k] = -2.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xx_xyzz[k] = -2.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xx_xzzz[k] = -2.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xx_yyyy[k] = 4.0 * to_xxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xx_yyyz[k] = 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xx_yyzz[k] = 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xx_yzzz[k] = 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xx_zzzz[k] = 4.0 * to_xxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 555-570 components of targeted buffer : DG

        auto to_z_x_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 15);

        auto to_z_x_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 16);

        auto to_z_x_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 17);

        auto to_z_x_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 18);

        auto to_z_x_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 19);

        auto to_z_x_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 20);

        auto to_z_x_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 21);

        auto to_z_x_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 22);

        auto to_z_x_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 23);

        auto to_z_x_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 24);

        auto to_z_x_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 25);

        auto to_z_x_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 26);

        auto to_z_x_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 27);

        auto to_z_x_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 28);

        auto to_z_x_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxx,   \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyz,     \
                             to_xyz_yzz,     \
                             to_xyz_zzz,     \
                             to_z_x_xy_xxxx, \
                             to_z_x_xy_xxxy, \
                             to_z_x_xy_xxxz, \
                             to_z_x_xy_xxyy, \
                             to_z_x_xy_xxyz, \
                             to_z_x_xy_xxzz, \
                             to_z_x_xy_xyyy, \
                             to_z_x_xy_xyyz, \
                             to_z_x_xy_xyzz, \
                             to_z_x_xy_xzzz, \
                             to_z_x_xy_yyyy, \
                             to_z_x_xy_yyyz, \
                             to_z_x_xy_yyzz, \
                             to_z_x_xy_yzzz, \
                             to_z_x_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xy_xxxx[k] = -8.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xy_xxxy[k] = -6.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xy_xxxz[k] = -6.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xy_xxyy[k] = -4.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xy_xxyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xy_xxzz[k] = -4.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xy_xyyy[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xy_xyyz[k] = -2.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xy_xyzz[k] = -2.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xy_xzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xy_yyyy[k] = 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xy_yyyz[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xy_yyzz[k] = 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xy_yzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xy_zzzz[k] = 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-585 components of targeted buffer : DG

        auto to_z_x_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 30);

        auto to_z_x_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 31);

        auto to_z_x_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 32);

        auto to_z_x_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 33);

        auto to_z_x_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 34);

        auto to_z_x_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 35);

        auto to_z_x_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 36);

        auto to_z_x_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 37);

        auto to_z_x_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 38);

        auto to_z_x_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 39);

        auto to_z_x_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 40);

        auto to_z_x_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 41);

        auto to_z_x_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 42);

        auto to_z_x_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 43);

        auto to_z_x_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxx,     \
                             to_x_xxxxy,     \
                             to_x_xxxxz,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyz,       \
                             to_x_yzz,       \
                             to_x_zzz,       \
                             to_xzz_xxx,     \
                             to_xzz_xxxxx,   \
                             to_xzz_xxxxy,   \
                             to_xzz_xxxxz,   \
                             to_xzz_xxxyy,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxxzz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyy,   \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xxzzz,   \
                             to_xzz_xyy,     \
                             to_xzz_xyyyy,   \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_xzzzz,   \
                             to_xzz_yyy,     \
                             to_xzz_yyz,     \
                             to_xzz_yzz,     \
                             to_xzz_zzz,     \
                             to_z_x_xz_xxxx, \
                             to_z_x_xz_xxxy, \
                             to_z_x_xz_xxxz, \
                             to_z_x_xz_xxyy, \
                             to_z_x_xz_xxyz, \
                             to_z_x_xz_xxzz, \
                             to_z_x_xz_xyyy, \
                             to_z_x_xz_xyyz, \
                             to_z_x_xz_xyzz, \
                             to_z_x_xz_xzzz, \
                             to_z_x_xz_yyyy, \
                             to_z_x_xz_yyyz, \
                             to_z_x_xz_yyzz, \
                             to_z_x_xz_yzzz, \
                             to_z_x_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xz_xxxx[k] = 4.0 * to_x_xxx[k] - 2.0 * to_x_xxxxx[k] * tke_0 - 8.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xz_xxxy[k] = 3.0 * to_x_xxy[k] - 2.0 * to_x_xxxxy[k] * tke_0 - 6.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xz_xxxz[k] = 3.0 * to_x_xxz[k] - 2.0 * to_x_xxxxz[k] * tke_0 - 6.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xz_xxyy[k] = 2.0 * to_x_xyy[k] - 2.0 * to_x_xxxyy[k] * tke_0 - 4.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xz_xxyz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xxxyz[k] * tke_0 - 4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xz_xxzz[k] = 2.0 * to_x_xzz[k] - 2.0 * to_x_xxxzz[k] * tke_0 - 4.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xz_xyyy[k] = to_x_yyy[k] - 2.0 * to_x_xxyyy[k] * tke_0 - 2.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xz_xyyz[k] = to_x_yyz[k] - 2.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xz_xyzz[k] = to_x_yzz[k] - 2.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xz_xzzz[k] = to_x_zzz[k] - 2.0 * to_x_xxzzz[k] * tke_0 - 2.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xz_yyyy[k] = -2.0 * to_x_xyyyy[k] * tke_0 + 4.0 * to_xzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xz_yyyz[k] = -2.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xz_yyzz[k] = -2.0 * to_x_xyyzz[k] * tke_0 + 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xz_yzzz[k] = -2.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xz_zzzz[k] = -2.0 * to_x_xzzzz[k] * tke_0 + 4.0 * to_xzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 585-600 components of targeted buffer : DG

        auto to_z_x_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 45);

        auto to_z_x_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 46);

        auto to_z_x_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 47);

        auto to_z_x_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 48);

        auto to_z_x_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 49);

        auto to_z_x_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 50);

        auto to_z_x_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 51);

        auto to_z_x_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 52);

        auto to_z_x_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 53);

        auto to_z_x_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 54);

        auto to_z_x_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 55);

        auto to_z_x_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 56);

        auto to_z_x_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 57);

        auto to_z_x_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 58);

        auto to_z_x_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_yyz_xxx,         \
                             to_yyz_xxxxx,   \
                             to_yyz_xxxxy,   \
                             to_yyz_xxxxz,   \
                             to_yyz_xxxyy,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxxzz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyy,   \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xxzzz,   \
                             to_yyz_xyy,     \
                             to_yyz_xyyyy,   \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_xzzzz,   \
                             to_yyz_yyy,     \
                             to_yyz_yyz,     \
                             to_yyz_yzz,     \
                             to_yyz_zzz,     \
                             to_z_x_yy_xxxx, \
                             to_z_x_yy_xxxy, \
                             to_z_x_yy_xxxz, \
                             to_z_x_yy_xxyy, \
                             to_z_x_yy_xxyz, \
                             to_z_x_yy_xxzz, \
                             to_z_x_yy_xyyy, \
                             to_z_x_yy_xyyz, \
                             to_z_x_yy_xyzz, \
                             to_z_x_yy_xzzz, \
                             to_z_x_yy_yyyy, \
                             to_z_x_yy_yyyz, \
                             to_z_x_yy_yyzz, \
                             to_z_x_yy_yzzz, \
                             to_z_x_yy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yy_xxxx[k] = -8.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yy_xxxy[k] = -6.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yy_xxxz[k] = -6.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yy_xxyy[k] = -4.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yy_xxyz[k] = -4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yy_xxzz[k] = -4.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yy_xyyy[k] = -2.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yy_xyyz[k] = -2.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yy_xyzz[k] = -2.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yy_xzzz[k] = -2.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yy_yyyy[k] = 4.0 * to_yyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yy_yyyz[k] = 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yy_yyzz[k] = 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yy_yzzz[k] = 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yy_zzzz[k] = 4.0 * to_yyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-615 components of targeted buffer : DG

        auto to_z_x_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 60);

        auto to_z_x_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 61);

        auto to_z_x_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 62);

        auto to_z_x_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 63);

        auto to_z_x_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 64);

        auto to_z_x_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 65);

        auto to_z_x_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 66);

        auto to_z_x_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 67);

        auto to_z_x_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 68);

        auto to_z_x_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 69);

        auto to_z_x_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 70);

        auto to_z_x_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 71);

        auto to_z_x_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 72);

        auto to_z_x_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 73);

        auto to_z_x_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_xxx,           \
                             to_y_xxxxx,     \
                             to_y_xxxxy,     \
                             to_y_xxxxz,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyz,       \
                             to_y_yzz,       \
                             to_y_zzz,       \
                             to_yzz_xxx,     \
                             to_yzz_xxxxx,   \
                             to_yzz_xxxxy,   \
                             to_yzz_xxxxz,   \
                             to_yzz_xxxyy,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxxzz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyy,   \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xxzzz,   \
                             to_yzz_xyy,     \
                             to_yzz_xyyyy,   \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_xzzzz,   \
                             to_yzz_yyy,     \
                             to_yzz_yyz,     \
                             to_yzz_yzz,     \
                             to_yzz_zzz,     \
                             to_z_x_yz_xxxx, \
                             to_z_x_yz_xxxy, \
                             to_z_x_yz_xxxz, \
                             to_z_x_yz_xxyy, \
                             to_z_x_yz_xxyz, \
                             to_z_x_yz_xxzz, \
                             to_z_x_yz_xyyy, \
                             to_z_x_yz_xyyz, \
                             to_z_x_yz_xyzz, \
                             to_z_x_yz_xzzz, \
                             to_z_x_yz_yyyy, \
                             to_z_x_yz_yyyz, \
                             to_z_x_yz_yyzz, \
                             to_z_x_yz_yzzz, \
                             to_z_x_yz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yz_xxxx[k] = 4.0 * to_y_xxx[k] - 2.0 * to_y_xxxxx[k] * tke_0 - 8.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yz_xxxy[k] = 3.0 * to_y_xxy[k] - 2.0 * to_y_xxxxy[k] * tke_0 - 6.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yz_xxxz[k] = 3.0 * to_y_xxz[k] - 2.0 * to_y_xxxxz[k] * tke_0 - 6.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yz_xxyy[k] = 2.0 * to_y_xyy[k] - 2.0 * to_y_xxxyy[k] * tke_0 - 4.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yz_xxyz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xxxyz[k] * tke_0 - 4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yz_xxzz[k] = 2.0 * to_y_xzz[k] - 2.0 * to_y_xxxzz[k] * tke_0 - 4.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yz_xyyy[k] = to_y_yyy[k] - 2.0 * to_y_xxyyy[k] * tke_0 - 2.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yz_xyyz[k] = to_y_yyz[k] - 2.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yz_xyzz[k] = to_y_yzz[k] - 2.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yz_xzzz[k] = to_y_zzz[k] - 2.0 * to_y_xxzzz[k] * tke_0 - 2.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yz_yyyy[k] = -2.0 * to_y_xyyyy[k] * tke_0 + 4.0 * to_yzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yz_yyyz[k] = -2.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yz_yyzz[k] = -2.0 * to_y_xyyzz[k] * tke_0 + 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yz_yzzz[k] = -2.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yz_zzzz[k] = -2.0 * to_y_xzzzz[k] * tke_0 + 4.0 * to_yzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 615-630 components of targeted buffer : DG

        auto to_z_x_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 75);

        auto to_z_x_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 76);

        auto to_z_x_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 77);

        auto to_z_x_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 78);

        auto to_z_x_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 79);

        auto to_z_x_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 80);

        auto to_z_x_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 81);

        auto to_z_x_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 82);

        auto to_z_x_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 83);

        auto to_z_x_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 84);

        auto to_z_x_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 85);

        auto to_z_x_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 86);

        auto to_z_x_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 87);

        auto to_z_x_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 88);

        auto to_z_x_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 6 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_z_x_zz_xxxx,     \
                             to_z_x_zz_xxxy, \
                             to_z_x_zz_xxxz, \
                             to_z_x_zz_xxyy, \
                             to_z_x_zz_xxyz, \
                             to_z_x_zz_xxzz, \
                             to_z_x_zz_xyyy, \
                             to_z_x_zz_xyyz, \
                             to_z_x_zz_xyzz, \
                             to_z_x_zz_xzzz, \
                             to_z_x_zz_yyyy, \
                             to_z_x_zz_yyyz, \
                             to_z_x_zz_yyzz, \
                             to_z_x_zz_yzzz, \
                             to_z_x_zz_zzzz, \
                             to_z_xxx,       \
                             to_z_xxxxx,     \
                             to_z_xxxxy,     \
                             to_z_xxxxz,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyz,       \
                             to_z_yzz,       \
                             to_z_zzz,       \
                             to_zzz_xxx,     \
                             to_zzz_xxxxx,   \
                             to_zzz_xxxxy,   \
                             to_zzz_xxxxz,   \
                             to_zzz_xxxyy,   \
                             to_zzz_xxxyz,   \
                             to_zzz_xxxzz,   \
                             to_zzz_xxy,     \
                             to_zzz_xxyyy,   \
                             to_zzz_xxyyz,   \
                             to_zzz_xxyzz,   \
                             to_zzz_xxz,     \
                             to_zzz_xxzzz,   \
                             to_zzz_xyy,     \
                             to_zzz_xyyyy,   \
                             to_zzz_xyyyz,   \
                             to_zzz_xyyzz,   \
                             to_zzz_xyz,     \
                             to_zzz_xyzzz,   \
                             to_zzz_xzz,     \
                             to_zzz_xzzzz,   \
                             to_zzz_yyy,     \
                             to_zzz_yyz,     \
                             to_zzz_yzz,     \
                             to_zzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zz_xxxx[k] = 8.0 * to_z_xxx[k] - 4.0 * to_z_xxxxx[k] * tke_0 - 8.0 * to_zzz_xxx[k] * tbe_0 + 4.0 * to_zzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_zz_xxxy[k] = 6.0 * to_z_xxy[k] - 4.0 * to_z_xxxxy[k] * tke_0 - 6.0 * to_zzz_xxy[k] * tbe_0 + 4.0 * to_zzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_zz_xxxz[k] = 6.0 * to_z_xxz[k] - 4.0 * to_z_xxxxz[k] * tke_0 - 6.0 * to_zzz_xxz[k] * tbe_0 + 4.0 * to_zzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_zz_xxyy[k] = 4.0 * to_z_xyy[k] - 4.0 * to_z_xxxyy[k] * tke_0 - 4.0 * to_zzz_xyy[k] * tbe_0 + 4.0 * to_zzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_zz_xxyz[k] = 4.0 * to_z_xyz[k] - 4.0 * to_z_xxxyz[k] * tke_0 - 4.0 * to_zzz_xyz[k] * tbe_0 + 4.0 * to_zzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_zz_xxzz[k] = 4.0 * to_z_xzz[k] - 4.0 * to_z_xxxzz[k] * tke_0 - 4.0 * to_zzz_xzz[k] * tbe_0 + 4.0 * to_zzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_zz_xyyy[k] = 2.0 * to_z_yyy[k] - 4.0 * to_z_xxyyy[k] * tke_0 - 2.0 * to_zzz_yyy[k] * tbe_0 + 4.0 * to_zzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_zz_xyyz[k] = 2.0 * to_z_yyz[k] - 4.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_zzz_yyz[k] * tbe_0 + 4.0 * to_zzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_zz_xyzz[k] = 2.0 * to_z_yzz[k] - 4.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_zzz_yzz[k] * tbe_0 + 4.0 * to_zzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_zz_xzzz[k] = 2.0 * to_z_zzz[k] - 4.0 * to_z_xxzzz[k] * tke_0 - 2.0 * to_zzz_zzz[k] * tbe_0 + 4.0 * to_zzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_zz_yyyy[k] = -4.0 * to_z_xyyyy[k] * tke_0 + 4.0 * to_zzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_zz_yyyz[k] = -4.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_zzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_zz_yyzz[k] = -4.0 * to_z_xyyzz[k] * tke_0 + 4.0 * to_zzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_zz_yzzz[k] = -4.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_zzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_zz_zzzz[k] = -4.0 * to_z_xzzzz[k] * tke_0 + 4.0 * to_zzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-645 components of targeted buffer : DG

        auto to_z_y_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 0);

        auto to_z_y_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 1);

        auto to_z_y_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 2);

        auto to_z_y_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 3);

        auto to_z_y_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 4);

        auto to_z_y_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 5);

        auto to_z_y_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 6);

        auto to_z_y_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 7);

        auto to_z_y_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 8);

        auto to_z_y_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 9);

        auto to_z_y_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 10);

        auto to_z_y_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 11);

        auto to_z_y_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 12);

        auto to_z_y_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 13);

        auto to_z_y_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxz_xxx,         \
                             to_xxz_xxxxy,   \
                             to_xxz_xxxyy,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyy,   \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xyy,     \
                             to_xxz_xyyyy,   \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_yyy,     \
                             to_xxz_yyyyy,   \
                             to_xxz_yyyyz,   \
                             to_xxz_yyyzz,   \
                             to_xxz_yyz,     \
                             to_xxz_yyzzz,   \
                             to_xxz_yzz,     \
                             to_xxz_yzzzz,   \
                             to_xxz_zzz,     \
                             to_z_y_xx_xxxx, \
                             to_z_y_xx_xxxy, \
                             to_z_y_xx_xxxz, \
                             to_z_y_xx_xxyy, \
                             to_z_y_xx_xxyz, \
                             to_z_y_xx_xxzz, \
                             to_z_y_xx_xyyy, \
                             to_z_y_xx_xyyz, \
                             to_z_y_xx_xyzz, \
                             to_z_y_xx_xzzz, \
                             to_z_y_xx_yyyy, \
                             to_z_y_xx_yyyz, \
                             to_z_y_xx_yyzz, \
                             to_z_y_xx_yzzz, \
                             to_z_y_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xx_xxxx[k] = 4.0 * to_xxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xx_xxxy[k] = -2.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xxxz[k] = 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xx_xxyy[k] = -4.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xxyz[k] = -2.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xx_xxzz[k] = 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xx_xyyy[k] = -6.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xyyz[k] = -4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xx_xyzz[k] = -2.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xx_xzzz[k] = 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xx_yyyy[k] = -8.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xx_yyyz[k] = -6.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xx_yyzz[k] = -4.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xx_yzzz[k] = -2.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xx_zzzz[k] = 4.0 * to_xxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 645-660 components of targeted buffer : DG

        auto to_z_y_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 15);

        auto to_z_y_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 16);

        auto to_z_y_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 17);

        auto to_z_y_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 18);

        auto to_z_y_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 19);

        auto to_z_y_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 20);

        auto to_z_y_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 21);

        auto to_z_y_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 22);

        auto to_z_y_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 23);

        auto to_z_y_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 24);

        auto to_z_y_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 25);

        auto to_z_y_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 26);

        auto to_z_y_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 27);

        auto to_z_y_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 28);

        auto to_z_y_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxy,   \
                             to_xyz_xxxyy,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyy,   \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xyy,     \
                             to_xyz_xyyyy,   \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_yyy,     \
                             to_xyz_yyyyy,   \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             to_z_y_xy_xxxx, \
                             to_z_y_xy_xxxy, \
                             to_z_y_xy_xxxz, \
                             to_z_y_xy_xxyy, \
                             to_z_y_xy_xxyz, \
                             to_z_y_xy_xxzz, \
                             to_z_y_xy_xyyy, \
                             to_z_y_xy_xyyz, \
                             to_z_y_xy_xyzz, \
                             to_z_y_xy_xzzz, \
                             to_z_y_xy_yyyy, \
                             to_z_y_xy_yyyz, \
                             to_z_y_xy_yyzz, \
                             to_z_y_xy_yzzz, \
                             to_z_y_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xy_xxxx[k] = 4.0 * to_xyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xy_xxxy[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xxxz[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xy_xxyy[k] = -4.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xxyz[k] = -2.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xy_xxzz[k] = 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xy_xyyy[k] = -6.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xyyz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xy_xyzz[k] = -2.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xy_xzzz[k] = 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xy_yyyy[k] = -8.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xy_yyyz[k] = -6.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xy_yyzz[k] = -4.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xy_yzzz[k] = -2.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xy_zzzz[k] = 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-675 components of targeted buffer : DG

        auto to_z_y_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 30);

        auto to_z_y_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 31);

        auto to_z_y_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 32);

        auto to_z_y_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 33);

        auto to_z_y_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 34);

        auto to_z_y_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 35);

        auto to_z_y_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 36);

        auto to_z_y_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 37);

        auto to_z_y_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 38);

        auto to_z_y_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 39);

        auto to_z_y_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 40);

        auto to_z_y_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 41);

        auto to_z_y_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 42);

        auto to_z_y_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 43);

        auto to_z_y_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxy,     \
                             to_x_xxxyy,     \
                             to_x_xxxyz,     \
                             to_x_xxy,       \
                             to_x_xxyyy,     \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xyy,       \
                             to_x_xyyyy,     \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_yyy,       \
                             to_x_yyyyy,     \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_zzz,       \
                             to_xzz_xxx,     \
                             to_xzz_xxxxy,   \
                             to_xzz_xxxyy,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyy,   \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xyy,     \
                             to_xzz_xyyyy,   \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_yyy,     \
                             to_xzz_yyyyy,   \
                             to_xzz_yyyyz,   \
                             to_xzz_yyyzz,   \
                             to_xzz_yyz,     \
                             to_xzz_yyzzz,   \
                             to_xzz_yzz,     \
                             to_xzz_yzzzz,   \
                             to_xzz_zzz,     \
                             to_z_y_xz_xxxx, \
                             to_z_y_xz_xxxy, \
                             to_z_y_xz_xxxz, \
                             to_z_y_xz_xxyy, \
                             to_z_y_xz_xxyz, \
                             to_z_y_xz_xxzz, \
                             to_z_y_xz_xyyy, \
                             to_z_y_xz_xyyz, \
                             to_z_y_xz_xyzz, \
                             to_z_y_xz_xzzz, \
                             to_z_y_xz_yyyy, \
                             to_z_y_xz_yyyz, \
                             to_z_y_xz_yyzz, \
                             to_z_y_xz_yzzz, \
                             to_z_y_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xz_xxxx[k] = -2.0 * to_x_xxxxy[k] * tke_0 + 4.0 * to_xzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xz_xxxy[k] = to_x_xxx[k] - 2.0 * to_x_xxxyy[k] * tke_0 - 2.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xxxz[k] = -2.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xz_xxyy[k] = 2.0 * to_x_xxy[k] - 2.0 * to_x_xxyyy[k] * tke_0 - 4.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xxyz[k] = to_x_xxz[k] - 2.0 * to_x_xxyyz[k] * tke_0 - 2.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xz_xxzz[k] = -2.0 * to_x_xxyzz[k] * tke_0 + 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xz_xyyy[k] = 3.0 * to_x_xyy[k] - 2.0 * to_x_xyyyy[k] * tke_0 - 6.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xyyz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xyyyz[k] * tke_0 - 4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xz_xyzz[k] = to_x_xzz[k] - 2.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xz_xzzz[k] = -2.0 * to_x_xyzzz[k] * tke_0 + 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xz_yyyy[k] = 4.0 * to_x_yyy[k] - 2.0 * to_x_yyyyy[k] * tke_0 - 8.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xz_yyyz[k] = 3.0 * to_x_yyz[k] - 2.0 * to_x_yyyyz[k] * tke_0 - 6.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xz_yyzz[k] = 2.0 * to_x_yzz[k] - 2.0 * to_x_yyyzz[k] * tke_0 - 4.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xz_yzzz[k] = to_x_zzz[k] - 2.0 * to_x_yyzzz[k] * tke_0 - 2.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xz_zzzz[k] = -2.0 * to_x_yzzzz[k] * tke_0 + 4.0 * to_xzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 675-690 components of targeted buffer : DG

        auto to_z_y_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 45);

        auto to_z_y_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 46);

        auto to_z_y_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 47);

        auto to_z_y_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 48);

        auto to_z_y_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 49);

        auto to_z_y_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 50);

        auto to_z_y_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 51);

        auto to_z_y_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 52);

        auto to_z_y_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 53);

        auto to_z_y_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 54);

        auto to_z_y_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 55);

        auto to_z_y_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 56);

        auto to_z_y_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 57);

        auto to_z_y_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 58);

        auto to_z_y_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_yyz_xxx,         \
                             to_yyz_xxxxy,   \
                             to_yyz_xxxyy,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyy,   \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xyy,     \
                             to_yyz_xyyyy,   \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_yyy,     \
                             to_yyz_yyyyy,   \
                             to_yyz_yyyyz,   \
                             to_yyz_yyyzz,   \
                             to_yyz_yyz,     \
                             to_yyz_yyzzz,   \
                             to_yyz_yzz,     \
                             to_yyz_yzzzz,   \
                             to_yyz_zzz,     \
                             to_z_y_yy_xxxx, \
                             to_z_y_yy_xxxy, \
                             to_z_y_yy_xxxz, \
                             to_z_y_yy_xxyy, \
                             to_z_y_yy_xxyz, \
                             to_z_y_yy_xxzz, \
                             to_z_y_yy_xyyy, \
                             to_z_y_yy_xyyz, \
                             to_z_y_yy_xyzz, \
                             to_z_y_yy_xzzz, \
                             to_z_y_yy_yyyy, \
                             to_z_y_yy_yyyz, \
                             to_z_y_yy_yyzz, \
                             to_z_y_yy_yzzz, \
                             to_z_y_yy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yy_xxxx[k] = 4.0 * to_yyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yy_xxxy[k] = -2.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xxxz[k] = 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yy_xxyy[k] = -4.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xxyz[k] = -2.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yy_xxzz[k] = 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yy_xyyy[k] = -6.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xyyz[k] = -4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yy_xyzz[k] = -2.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yy_xzzz[k] = 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yy_yyyy[k] = -8.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yy_yyyz[k] = -6.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yy_yyzz[k] = -4.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yy_yzzz[k] = -2.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yy_zzzz[k] = 4.0 * to_yyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-705 components of targeted buffer : DG

        auto to_z_y_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 60);

        auto to_z_y_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 61);

        auto to_z_y_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 62);

        auto to_z_y_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 63);

        auto to_z_y_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 64);

        auto to_z_y_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 65);

        auto to_z_y_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 66);

        auto to_z_y_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 67);

        auto to_z_y_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 68);

        auto to_z_y_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 69);

        auto to_z_y_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 70);

        auto to_z_y_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 71);

        auto to_z_y_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 72);

        auto to_z_y_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 73);

        auto to_z_y_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_xxx,           \
                             to_y_xxxxy,     \
                             to_y_xxxyy,     \
                             to_y_xxxyz,     \
                             to_y_xxy,       \
                             to_y_xxyyy,     \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xyy,       \
                             to_y_xyyyy,     \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_yyy,       \
                             to_y_yyyyy,     \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_zzz,       \
                             to_yzz_xxx,     \
                             to_yzz_xxxxy,   \
                             to_yzz_xxxyy,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyy,   \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xyy,     \
                             to_yzz_xyyyy,   \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_yyy,     \
                             to_yzz_yyyyy,   \
                             to_yzz_yyyyz,   \
                             to_yzz_yyyzz,   \
                             to_yzz_yyz,     \
                             to_yzz_yyzzz,   \
                             to_yzz_yzz,     \
                             to_yzz_yzzzz,   \
                             to_yzz_zzz,     \
                             to_z_y_yz_xxxx, \
                             to_z_y_yz_xxxy, \
                             to_z_y_yz_xxxz, \
                             to_z_y_yz_xxyy, \
                             to_z_y_yz_xxyz, \
                             to_z_y_yz_xxzz, \
                             to_z_y_yz_xyyy, \
                             to_z_y_yz_xyyz, \
                             to_z_y_yz_xyzz, \
                             to_z_y_yz_xzzz, \
                             to_z_y_yz_yyyy, \
                             to_z_y_yz_yyyz, \
                             to_z_y_yz_yyzz, \
                             to_z_y_yz_yzzz, \
                             to_z_y_yz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yz_xxxx[k] = -2.0 * to_y_xxxxy[k] * tke_0 + 4.0 * to_yzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yz_xxxy[k] = to_y_xxx[k] - 2.0 * to_y_xxxyy[k] * tke_0 - 2.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xxxz[k] = -2.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yz_xxyy[k] = 2.0 * to_y_xxy[k] - 2.0 * to_y_xxyyy[k] * tke_0 - 4.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xxyz[k] = to_y_xxz[k] - 2.0 * to_y_xxyyz[k] * tke_0 - 2.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yz_xxzz[k] = -2.0 * to_y_xxyzz[k] * tke_0 + 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yz_xyyy[k] = 3.0 * to_y_xyy[k] - 2.0 * to_y_xyyyy[k] * tke_0 - 6.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xyyz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xyyyz[k] * tke_0 - 4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yz_xyzz[k] = to_y_xzz[k] - 2.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yz_xzzz[k] = -2.0 * to_y_xyzzz[k] * tke_0 + 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yz_yyyy[k] = 4.0 * to_y_yyy[k] - 2.0 * to_y_yyyyy[k] * tke_0 - 8.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yz_yyyz[k] = 3.0 * to_y_yyz[k] - 2.0 * to_y_yyyyz[k] * tke_0 - 6.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yz_yyzz[k] = 2.0 * to_y_yzz[k] - 2.0 * to_y_yyyzz[k] * tke_0 - 4.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yz_yzzz[k] = to_y_zzz[k] - 2.0 * to_y_yyzzz[k] * tke_0 - 2.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yz_zzzz[k] = -2.0 * to_y_yzzzz[k] * tke_0 + 4.0 * to_yzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 705-720 components of targeted buffer : DG

        auto to_z_y_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 75);

        auto to_z_y_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 76);

        auto to_z_y_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 77);

        auto to_z_y_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 78);

        auto to_z_y_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 79);

        auto to_z_y_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 80);

        auto to_z_y_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 81);

        auto to_z_y_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 82);

        auto to_z_y_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 83);

        auto to_z_y_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 84);

        auto to_z_y_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 85);

        auto to_z_y_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 86);

        auto to_z_y_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 87);

        auto to_z_y_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 88);

        auto to_z_y_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 7 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_z_xxx,           \
                             to_z_xxxxy,     \
                             to_z_xxxyy,     \
                             to_z_xxxyz,     \
                             to_z_xxy,       \
                             to_z_xxyyy,     \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xyy,       \
                             to_z_xyyyy,     \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_y_zz_xxxx, \
                             to_z_y_zz_xxxy, \
                             to_z_y_zz_xxxz, \
                             to_z_y_zz_xxyy, \
                             to_z_y_zz_xxyz, \
                             to_z_y_zz_xxzz, \
                             to_z_y_zz_xyyy, \
                             to_z_y_zz_xyyz, \
                             to_z_y_zz_xyzz, \
                             to_z_y_zz_xzzz, \
                             to_z_y_zz_yyyy, \
                             to_z_y_zz_yyyz, \
                             to_z_y_zz_yyzz, \
                             to_z_y_zz_yzzz, \
                             to_z_y_zz_zzzz, \
                             to_z_yyy,       \
                             to_z_yyyyy,     \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_zzz,       \
                             to_zzz_xxx,     \
                             to_zzz_xxxxy,   \
                             to_zzz_xxxyy,   \
                             to_zzz_xxxyz,   \
                             to_zzz_xxy,     \
                             to_zzz_xxyyy,   \
                             to_zzz_xxyyz,   \
                             to_zzz_xxyzz,   \
                             to_zzz_xxz,     \
                             to_zzz_xyy,     \
                             to_zzz_xyyyy,   \
                             to_zzz_xyyyz,   \
                             to_zzz_xyyzz,   \
                             to_zzz_xyz,     \
                             to_zzz_xyzzz,   \
                             to_zzz_xzz,     \
                             to_zzz_yyy,     \
                             to_zzz_yyyyy,   \
                             to_zzz_yyyyz,   \
                             to_zzz_yyyzz,   \
                             to_zzz_yyz,     \
                             to_zzz_yyzzz,   \
                             to_zzz_yzz,     \
                             to_zzz_yzzzz,   \
                             to_zzz_zzz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zz_xxxx[k] = -4.0 * to_z_xxxxy[k] * tke_0 + 4.0 * to_zzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_zz_xxxy[k] = 2.0 * to_z_xxx[k] - 4.0 * to_z_xxxyy[k] * tke_0 - 2.0 * to_zzz_xxx[k] * tbe_0 + 4.0 * to_zzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xxxz[k] = -4.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_zzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_zz_xxyy[k] = 4.0 * to_z_xxy[k] - 4.0 * to_z_xxyyy[k] * tke_0 - 4.0 * to_zzz_xxy[k] * tbe_0 + 4.0 * to_zzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xxyz[k] = 2.0 * to_z_xxz[k] - 4.0 * to_z_xxyyz[k] * tke_0 - 2.0 * to_zzz_xxz[k] * tbe_0 + 4.0 * to_zzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_zz_xxzz[k] = -4.0 * to_z_xxyzz[k] * tke_0 + 4.0 * to_zzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_zz_xyyy[k] = 6.0 * to_z_xyy[k] - 4.0 * to_z_xyyyy[k] * tke_0 - 6.0 * to_zzz_xyy[k] * tbe_0 + 4.0 * to_zzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xyyz[k] = 4.0 * to_z_xyz[k] - 4.0 * to_z_xyyyz[k] * tke_0 - 4.0 * to_zzz_xyz[k] * tbe_0 + 4.0 * to_zzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_zz_xyzz[k] = 2.0 * to_z_xzz[k] - 4.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_zzz_xzz[k] * tbe_0 + 4.0 * to_zzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_zz_xzzz[k] = -4.0 * to_z_xyzzz[k] * tke_0 + 4.0 * to_zzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_zz_yyyy[k] = 8.0 * to_z_yyy[k] - 4.0 * to_z_yyyyy[k] * tke_0 - 8.0 * to_zzz_yyy[k] * tbe_0 + 4.0 * to_zzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_zz_yyyz[k] = 6.0 * to_z_yyz[k] - 4.0 * to_z_yyyyz[k] * tke_0 - 6.0 * to_zzz_yyz[k] * tbe_0 + 4.0 * to_zzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_zz_yyzz[k] = 4.0 * to_z_yzz[k] - 4.0 * to_z_yyyzz[k] * tke_0 - 4.0 * to_zzz_yzz[k] * tbe_0 + 4.0 * to_zzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_zz_yzzz[k] = 2.0 * to_z_zzz[k] - 4.0 * to_z_yyzzz[k] * tke_0 - 2.0 * to_zzz_zzz[k] * tbe_0 + 4.0 * to_zzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_zz_zzzz[k] = -4.0 * to_z_yzzzz[k] * tke_0 + 4.0 * to_zzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-735 components of targeted buffer : DG

        auto to_z_z_xx_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 0);

        auto to_z_z_xx_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 1);

        auto to_z_z_xx_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 2);

        auto to_z_z_xx_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 3);

        auto to_z_z_xx_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 4);

        auto to_z_z_xx_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 5);

        auto to_z_z_xx_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 6);

        auto to_z_z_xx_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 7);

        auto to_z_z_xx_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 8);

        auto to_z_z_xx_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 9);

        auto to_z_z_xx_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 10);

        auto to_z_z_xx_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 11);

        auto to_z_z_xx_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 12);

        auto to_z_z_xx_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 13);

        auto to_z_z_xx_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 14);

#pragma omp simd aligned(to_xxz_xxx,         \
                             to_xxz_xxxxz,   \
                             to_xxz_xxxyz,   \
                             to_xxz_xxxzz,   \
                             to_xxz_xxy,     \
                             to_xxz_xxyyz,   \
                             to_xxz_xxyzz,   \
                             to_xxz_xxz,     \
                             to_xxz_xxzzz,   \
                             to_xxz_xyy,     \
                             to_xxz_xyyyz,   \
                             to_xxz_xyyzz,   \
                             to_xxz_xyz,     \
                             to_xxz_xyzzz,   \
                             to_xxz_xzz,     \
                             to_xxz_xzzzz,   \
                             to_xxz_yyy,     \
                             to_xxz_yyyyz,   \
                             to_xxz_yyyzz,   \
                             to_xxz_yyz,     \
                             to_xxz_yyzzz,   \
                             to_xxz_yzz,     \
                             to_xxz_yzzzz,   \
                             to_xxz_zzz,     \
                             to_xxz_zzzzz,   \
                             to_z_z_xx_xxxx, \
                             to_z_z_xx_xxxy, \
                             to_z_z_xx_xxxz, \
                             to_z_z_xx_xxyy, \
                             to_z_z_xx_xxyz, \
                             to_z_z_xx_xxzz, \
                             to_z_z_xx_xyyy, \
                             to_z_z_xx_xyyz, \
                             to_z_z_xx_xyzz, \
                             to_z_z_xx_xzzz, \
                             to_z_z_xx_yyyy, \
                             to_z_z_xx_yyyz, \
                             to_z_z_xx_yyzz, \
                             to_z_z_xx_yzzz, \
                             to_z_z_xx_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xx_xxxx[k] = 4.0 * to_xxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxxy[k] = 4.0 * to_xxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxxz[k] = -2.0 * to_xxz_xxx[k] * tbe_0 + 4.0 * to_xxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxyy[k] = 4.0 * to_xxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxyz[k] = -2.0 * to_xxz_xxy[k] * tbe_0 + 4.0 * to_xxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxzz[k] = -4.0 * to_xxz_xxz[k] * tbe_0 + 4.0 * to_xxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xyyy[k] = 4.0 * to_xxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xyyz[k] = -2.0 * to_xxz_xyy[k] * tbe_0 + 4.0 * to_xxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xyzz[k] = -4.0 * to_xxz_xyz[k] * tbe_0 + 4.0 * to_xxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xzzz[k] = -6.0 * to_xxz_xzz[k] * tbe_0 + 4.0 * to_xxz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yyyy[k] = 4.0 * to_xxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xx_yyyz[k] = -2.0 * to_xxz_yyy[k] * tbe_0 + 4.0 * to_xxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yyzz[k] = -4.0 * to_xxz_yyz[k] * tbe_0 + 4.0 * to_xxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yzzz[k] = -6.0 * to_xxz_yzz[k] * tbe_0 + 4.0 * to_xxz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_zzzz[k] = -8.0 * to_xxz_zzz[k] * tbe_0 + 4.0 * to_xxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 735-750 components of targeted buffer : DG

        auto to_z_z_xy_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 15);

        auto to_z_z_xy_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 16);

        auto to_z_z_xy_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 17);

        auto to_z_z_xy_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 18);

        auto to_z_z_xy_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 19);

        auto to_z_z_xy_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 20);

        auto to_z_z_xy_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 21);

        auto to_z_z_xy_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 22);

        auto to_z_z_xy_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 23);

        auto to_z_z_xy_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 24);

        auto to_z_z_xy_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 25);

        auto to_z_z_xy_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 26);

        auto to_z_z_xy_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 27);

        auto to_z_z_xy_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 28);

        auto to_z_z_xy_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 29);

#pragma omp simd aligned(to_xyz_xxx,         \
                             to_xyz_xxxxz,   \
                             to_xyz_xxxyz,   \
                             to_xyz_xxxzz,   \
                             to_xyz_xxy,     \
                             to_xyz_xxyyz,   \
                             to_xyz_xxyzz,   \
                             to_xyz_xxz,     \
                             to_xyz_xxzzz,   \
                             to_xyz_xyy,     \
                             to_xyz_xyyyz,   \
                             to_xyz_xyyzz,   \
                             to_xyz_xyz,     \
                             to_xyz_xyzzz,   \
                             to_xyz_xzz,     \
                             to_xyz_xzzzz,   \
                             to_xyz_yyy,     \
                             to_xyz_yyyyz,   \
                             to_xyz_yyyzz,   \
                             to_xyz_yyz,     \
                             to_xyz_yyzzz,   \
                             to_xyz_yzz,     \
                             to_xyz_yzzzz,   \
                             to_xyz_zzz,     \
                             to_xyz_zzzzz,   \
                             to_z_z_xy_xxxx, \
                             to_z_z_xy_xxxy, \
                             to_z_z_xy_xxxz, \
                             to_z_z_xy_xxyy, \
                             to_z_z_xy_xxyz, \
                             to_z_z_xy_xxzz, \
                             to_z_z_xy_xyyy, \
                             to_z_z_xy_xyyz, \
                             to_z_z_xy_xyzz, \
                             to_z_z_xy_xzzz, \
                             to_z_z_xy_yyyy, \
                             to_z_z_xy_yyyz, \
                             to_z_z_xy_yyzz, \
                             to_z_z_xy_yzzz, \
                             to_z_z_xy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xy_xxxx[k] = 4.0 * to_xyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxxy[k] = 4.0 * to_xyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxxz[k] = -2.0 * to_xyz_xxx[k] * tbe_0 + 4.0 * to_xyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxyy[k] = 4.0 * to_xyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxyz[k] = -2.0 * to_xyz_xxy[k] * tbe_0 + 4.0 * to_xyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxzz[k] = -4.0 * to_xyz_xxz[k] * tbe_0 + 4.0 * to_xyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xyyy[k] = 4.0 * to_xyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xyyz[k] = -2.0 * to_xyz_xyy[k] * tbe_0 + 4.0 * to_xyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xyzz[k] = -4.0 * to_xyz_xyz[k] * tbe_0 + 4.0 * to_xyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xzzz[k] = -6.0 * to_xyz_xzz[k] * tbe_0 + 4.0 * to_xyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yyyy[k] = 4.0 * to_xyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xy_yyyz[k] = -2.0 * to_xyz_yyy[k] * tbe_0 + 4.0 * to_xyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yyzz[k] = -4.0 * to_xyz_yyz[k] * tbe_0 + 4.0 * to_xyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yzzz[k] = -6.0 * to_xyz_yzz[k] * tbe_0 + 4.0 * to_xyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_zzzz[k] = -8.0 * to_xyz_zzz[k] * tbe_0 + 4.0 * to_xyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-765 components of targeted buffer : DG

        auto to_z_z_xz_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 30);

        auto to_z_z_xz_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 31);

        auto to_z_z_xz_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 32);

        auto to_z_z_xz_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 33);

        auto to_z_z_xz_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 34);

        auto to_z_z_xz_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 35);

        auto to_z_z_xz_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 36);

        auto to_z_z_xz_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 37);

        auto to_z_z_xz_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 38);

        auto to_z_z_xz_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 39);

        auto to_z_z_xz_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 40);

        auto to_z_z_xz_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 41);

        auto to_z_z_xz_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 42);

        auto to_z_z_xz_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 43);

        auto to_z_z_xz_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 44);

#pragma omp simd aligned(to_x_xxx,           \
                             to_x_xxxxz,     \
                             to_x_xxxyz,     \
                             to_x_xxxzz,     \
                             to_x_xxy,       \
                             to_x_xxyyz,     \
                             to_x_xxyzz,     \
                             to_x_xxz,       \
                             to_x_xxzzz,     \
                             to_x_xyy,       \
                             to_x_xyyyz,     \
                             to_x_xyyzz,     \
                             to_x_xyz,       \
                             to_x_xyzzz,     \
                             to_x_xzz,       \
                             to_x_xzzzz,     \
                             to_x_yyy,       \
                             to_x_yyyyz,     \
                             to_x_yyyzz,     \
                             to_x_yyz,       \
                             to_x_yyzzz,     \
                             to_x_yzz,       \
                             to_x_yzzzz,     \
                             to_x_zzz,       \
                             to_x_zzzzz,     \
                             to_xzz_xxx,     \
                             to_xzz_xxxxz,   \
                             to_xzz_xxxyz,   \
                             to_xzz_xxxzz,   \
                             to_xzz_xxy,     \
                             to_xzz_xxyyz,   \
                             to_xzz_xxyzz,   \
                             to_xzz_xxz,     \
                             to_xzz_xxzzz,   \
                             to_xzz_xyy,     \
                             to_xzz_xyyyz,   \
                             to_xzz_xyyzz,   \
                             to_xzz_xyz,     \
                             to_xzz_xyzzz,   \
                             to_xzz_xzz,     \
                             to_xzz_xzzzz,   \
                             to_xzz_yyy,     \
                             to_xzz_yyyyz,   \
                             to_xzz_yyyzz,   \
                             to_xzz_yyz,     \
                             to_xzz_yyzzz,   \
                             to_xzz_yzz,     \
                             to_xzz_yzzzz,   \
                             to_xzz_zzz,     \
                             to_xzz_zzzzz,   \
                             to_z_z_xz_xxxx, \
                             to_z_z_xz_xxxy, \
                             to_z_z_xz_xxxz, \
                             to_z_z_xz_xxyy, \
                             to_z_z_xz_xxyz, \
                             to_z_z_xz_xxzz, \
                             to_z_z_xz_xyyy, \
                             to_z_z_xz_xyyz, \
                             to_z_z_xz_xyzz, \
                             to_z_z_xz_xzzz, \
                             to_z_z_xz_yyyy, \
                             to_z_z_xz_yyyz, \
                             to_z_z_xz_yyzz, \
                             to_z_z_xz_yzzz, \
                             to_z_z_xz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xz_xxxx[k] = -2.0 * to_x_xxxxz[k] * tke_0 + 4.0 * to_xzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxxy[k] = -2.0 * to_x_xxxyz[k] * tke_0 + 4.0 * to_xzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxxz[k] = to_x_xxx[k] - 2.0 * to_x_xxxzz[k] * tke_0 - 2.0 * to_xzz_xxx[k] * tbe_0 + 4.0 * to_xzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxyy[k] = -2.0 * to_x_xxyyz[k] * tke_0 + 4.0 * to_xzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxyz[k] = to_x_xxy[k] - 2.0 * to_x_xxyzz[k] * tke_0 - 2.0 * to_xzz_xxy[k] * tbe_0 + 4.0 * to_xzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxzz[k] = 2.0 * to_x_xxz[k] - 2.0 * to_x_xxzzz[k] * tke_0 - 4.0 * to_xzz_xxz[k] * tbe_0 + 4.0 * to_xzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xyyy[k] = -2.0 * to_x_xyyyz[k] * tke_0 + 4.0 * to_xzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xyyz[k] = to_x_xyy[k] - 2.0 * to_x_xyyzz[k] * tke_0 - 2.0 * to_xzz_xyy[k] * tbe_0 + 4.0 * to_xzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xyzz[k] = 2.0 * to_x_xyz[k] - 2.0 * to_x_xyzzz[k] * tke_0 - 4.0 * to_xzz_xyz[k] * tbe_0 + 4.0 * to_xzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xzzz[k] = 3.0 * to_x_xzz[k] - 2.0 * to_x_xzzzz[k] * tke_0 - 6.0 * to_xzz_xzz[k] * tbe_0 + 4.0 * to_xzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yyyy[k] = -2.0 * to_x_yyyyz[k] * tke_0 + 4.0 * to_xzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xz_yyyz[k] = to_x_yyy[k] - 2.0 * to_x_yyyzz[k] * tke_0 - 2.0 * to_xzz_yyy[k] * tbe_0 + 4.0 * to_xzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yyzz[k] = 2.0 * to_x_yyz[k] - 2.0 * to_x_yyzzz[k] * tke_0 - 4.0 * to_xzz_yyz[k] * tbe_0 + 4.0 * to_xzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yzzz[k] = 3.0 * to_x_yzz[k] - 2.0 * to_x_yzzzz[k] * tke_0 - 6.0 * to_xzz_yzz[k] * tbe_0 + 4.0 * to_xzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_zzzz[k] = 4.0 * to_x_zzz[k] - 2.0 * to_x_zzzzz[k] * tke_0 - 8.0 * to_xzz_zzz[k] * tbe_0 + 4.0 * to_xzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 765-780 components of targeted buffer : DG

        auto to_z_z_yy_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 45);

        auto to_z_z_yy_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 46);

        auto to_z_z_yy_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 47);

        auto to_z_z_yy_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 48);

        auto to_z_z_yy_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 49);

        auto to_z_z_yy_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 50);

        auto to_z_z_yy_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 51);

        auto to_z_z_yy_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 52);

        auto to_z_z_yy_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 53);

        auto to_z_z_yy_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 54);

        auto to_z_z_yy_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 55);

        auto to_z_z_yy_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 56);

        auto to_z_z_yy_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 57);

        auto to_z_z_yy_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 58);

        auto to_z_z_yy_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 59);

#pragma omp simd aligned(to_yyz_xxx,         \
                             to_yyz_xxxxz,   \
                             to_yyz_xxxyz,   \
                             to_yyz_xxxzz,   \
                             to_yyz_xxy,     \
                             to_yyz_xxyyz,   \
                             to_yyz_xxyzz,   \
                             to_yyz_xxz,     \
                             to_yyz_xxzzz,   \
                             to_yyz_xyy,     \
                             to_yyz_xyyyz,   \
                             to_yyz_xyyzz,   \
                             to_yyz_xyz,     \
                             to_yyz_xyzzz,   \
                             to_yyz_xzz,     \
                             to_yyz_xzzzz,   \
                             to_yyz_yyy,     \
                             to_yyz_yyyyz,   \
                             to_yyz_yyyzz,   \
                             to_yyz_yyz,     \
                             to_yyz_yyzzz,   \
                             to_yyz_yzz,     \
                             to_yyz_yzzzz,   \
                             to_yyz_zzz,     \
                             to_yyz_zzzzz,   \
                             to_z_z_yy_xxxx, \
                             to_z_z_yy_xxxy, \
                             to_z_z_yy_xxxz, \
                             to_z_z_yy_xxyy, \
                             to_z_z_yy_xxyz, \
                             to_z_z_yy_xxzz, \
                             to_z_z_yy_xyyy, \
                             to_z_z_yy_xyyz, \
                             to_z_z_yy_xyzz, \
                             to_z_z_yy_xzzz, \
                             to_z_z_yy_yyyy, \
                             to_z_z_yy_yyyz, \
                             to_z_z_yy_yyzz, \
                             to_z_z_yy_yzzz, \
                             to_z_z_yy_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yy_xxxx[k] = 4.0 * to_yyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxxy[k] = 4.0 * to_yyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxxz[k] = -2.0 * to_yyz_xxx[k] * tbe_0 + 4.0 * to_yyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxyy[k] = 4.0 * to_yyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxyz[k] = -2.0 * to_yyz_xxy[k] * tbe_0 + 4.0 * to_yyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxzz[k] = -4.0 * to_yyz_xxz[k] * tbe_0 + 4.0 * to_yyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xyyy[k] = 4.0 * to_yyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xyyz[k] = -2.0 * to_yyz_xyy[k] * tbe_0 + 4.0 * to_yyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xyzz[k] = -4.0 * to_yyz_xyz[k] * tbe_0 + 4.0 * to_yyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xzzz[k] = -6.0 * to_yyz_xzz[k] * tbe_0 + 4.0 * to_yyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yyyy[k] = 4.0 * to_yyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yy_yyyz[k] = -2.0 * to_yyz_yyy[k] * tbe_0 + 4.0 * to_yyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yyzz[k] = -4.0 * to_yyz_yyz[k] * tbe_0 + 4.0 * to_yyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yzzz[k] = -6.0 * to_yyz_yzz[k] * tbe_0 + 4.0 * to_yyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_zzzz[k] = -8.0 * to_yyz_zzz[k] * tbe_0 + 4.0 * to_yyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-795 components of targeted buffer : DG

        auto to_z_z_yz_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 60);

        auto to_z_z_yz_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 61);

        auto to_z_z_yz_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 62);

        auto to_z_z_yz_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 63);

        auto to_z_z_yz_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 64);

        auto to_z_z_yz_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 65);

        auto to_z_z_yz_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 66);

        auto to_z_z_yz_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 67);

        auto to_z_z_yz_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 68);

        auto to_z_z_yz_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 69);

        auto to_z_z_yz_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 70);

        auto to_z_z_yz_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 71);

        auto to_z_z_yz_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 72);

        auto to_z_z_yz_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 73);

        auto to_z_z_yz_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 74);

#pragma omp simd aligned(to_y_xxx,           \
                             to_y_xxxxz,     \
                             to_y_xxxyz,     \
                             to_y_xxxzz,     \
                             to_y_xxy,       \
                             to_y_xxyyz,     \
                             to_y_xxyzz,     \
                             to_y_xxz,       \
                             to_y_xxzzz,     \
                             to_y_xyy,       \
                             to_y_xyyyz,     \
                             to_y_xyyzz,     \
                             to_y_xyz,       \
                             to_y_xyzzz,     \
                             to_y_xzz,       \
                             to_y_xzzzz,     \
                             to_y_yyy,       \
                             to_y_yyyyz,     \
                             to_y_yyyzz,     \
                             to_y_yyz,       \
                             to_y_yyzzz,     \
                             to_y_yzz,       \
                             to_y_yzzzz,     \
                             to_y_zzz,       \
                             to_y_zzzzz,     \
                             to_yzz_xxx,     \
                             to_yzz_xxxxz,   \
                             to_yzz_xxxyz,   \
                             to_yzz_xxxzz,   \
                             to_yzz_xxy,     \
                             to_yzz_xxyyz,   \
                             to_yzz_xxyzz,   \
                             to_yzz_xxz,     \
                             to_yzz_xxzzz,   \
                             to_yzz_xyy,     \
                             to_yzz_xyyyz,   \
                             to_yzz_xyyzz,   \
                             to_yzz_xyz,     \
                             to_yzz_xyzzz,   \
                             to_yzz_xzz,     \
                             to_yzz_xzzzz,   \
                             to_yzz_yyy,     \
                             to_yzz_yyyyz,   \
                             to_yzz_yyyzz,   \
                             to_yzz_yyz,     \
                             to_yzz_yyzzz,   \
                             to_yzz_yzz,     \
                             to_yzz_yzzzz,   \
                             to_yzz_zzz,     \
                             to_yzz_zzzzz,   \
                             to_z_z_yz_xxxx, \
                             to_z_z_yz_xxxy, \
                             to_z_z_yz_xxxz, \
                             to_z_z_yz_xxyy, \
                             to_z_z_yz_xxyz, \
                             to_z_z_yz_xxzz, \
                             to_z_z_yz_xyyy, \
                             to_z_z_yz_xyyz, \
                             to_z_z_yz_xyzz, \
                             to_z_z_yz_xzzz, \
                             to_z_z_yz_yyyy, \
                             to_z_z_yz_yyyz, \
                             to_z_z_yz_yyzz, \
                             to_z_z_yz_yzzz, \
                             to_z_z_yz_zzzz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yz_xxxx[k] = -2.0 * to_y_xxxxz[k] * tke_0 + 4.0 * to_yzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxxy[k] = -2.0 * to_y_xxxyz[k] * tke_0 + 4.0 * to_yzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxxz[k] = to_y_xxx[k] - 2.0 * to_y_xxxzz[k] * tke_0 - 2.0 * to_yzz_xxx[k] * tbe_0 + 4.0 * to_yzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxyy[k] = -2.0 * to_y_xxyyz[k] * tke_0 + 4.0 * to_yzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxyz[k] = to_y_xxy[k] - 2.0 * to_y_xxyzz[k] * tke_0 - 2.0 * to_yzz_xxy[k] * tbe_0 + 4.0 * to_yzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxzz[k] = 2.0 * to_y_xxz[k] - 2.0 * to_y_xxzzz[k] * tke_0 - 4.0 * to_yzz_xxz[k] * tbe_0 + 4.0 * to_yzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xyyy[k] = -2.0 * to_y_xyyyz[k] * tke_0 + 4.0 * to_yzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xyyz[k] = to_y_xyy[k] - 2.0 * to_y_xyyzz[k] * tke_0 - 2.0 * to_yzz_xyy[k] * tbe_0 + 4.0 * to_yzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xyzz[k] = 2.0 * to_y_xyz[k] - 2.0 * to_y_xyzzz[k] * tke_0 - 4.0 * to_yzz_xyz[k] * tbe_0 + 4.0 * to_yzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xzzz[k] = 3.0 * to_y_xzz[k] - 2.0 * to_y_xzzzz[k] * tke_0 - 6.0 * to_yzz_xzz[k] * tbe_0 + 4.0 * to_yzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yyyy[k] = -2.0 * to_y_yyyyz[k] * tke_0 + 4.0 * to_yzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yz_yyyz[k] = to_y_yyy[k] - 2.0 * to_y_yyyzz[k] * tke_0 - 2.0 * to_yzz_yyy[k] * tbe_0 + 4.0 * to_yzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yyzz[k] = 2.0 * to_y_yyz[k] - 2.0 * to_y_yyzzz[k] * tke_0 - 4.0 * to_yzz_yyz[k] * tbe_0 + 4.0 * to_yzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yzzz[k] = 3.0 * to_y_yzz[k] - 2.0 * to_y_yzzzz[k] * tke_0 - 6.0 * to_yzz_yzz[k] * tbe_0 + 4.0 * to_yzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_zzzz[k] = 4.0 * to_y_zzz[k] - 2.0 * to_y_zzzzz[k] * tke_0 - 8.0 * to_yzz_zzz[k] * tbe_0 + 4.0 * to_yzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 795-810 components of targeted buffer : DG

        auto to_z_z_zz_xxxx = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 75);

        auto to_z_z_zz_xxxy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 76);

        auto to_z_z_zz_xxxz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 77);

        auto to_z_z_zz_xxyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 78);

        auto to_z_z_zz_xxyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 79);

        auto to_z_z_zz_xxzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 80);

        auto to_z_z_zz_xyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 81);

        auto to_z_z_zz_xyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 82);

        auto to_z_z_zz_xyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 83);

        auto to_z_z_zz_xzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 84);

        auto to_z_z_zz_yyyy = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 85);

        auto to_z_z_zz_yyyz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 86);

        auto to_z_z_zz_yyzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 87);

        auto to_z_z_zz_yzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 88);

        auto to_z_z_zz_zzzz = pbuffer.data(idx_op_geom_101_dg + 8 * op_comps * 90 + i * 90 + 89);

#pragma omp simd aligned(to_z_xxx,           \
                             to_z_xxxxz,     \
                             to_z_xxxyz,     \
                             to_z_xxxzz,     \
                             to_z_xxy,       \
                             to_z_xxyyz,     \
                             to_z_xxyzz,     \
                             to_z_xxz,       \
                             to_z_xxzzz,     \
                             to_z_xyy,       \
                             to_z_xyyyz,     \
                             to_z_xyyzz,     \
                             to_z_xyz,       \
                             to_z_xyzzz,     \
                             to_z_xzz,       \
                             to_z_xzzzz,     \
                             to_z_yyy,       \
                             to_z_yyyyz,     \
                             to_z_yyyzz,     \
                             to_z_yyz,       \
                             to_z_yyzzz,     \
                             to_z_yzz,       \
                             to_z_yzzzz,     \
                             to_z_z_zz_xxxx, \
                             to_z_z_zz_xxxy, \
                             to_z_z_zz_xxxz, \
                             to_z_z_zz_xxyy, \
                             to_z_z_zz_xxyz, \
                             to_z_z_zz_xxzz, \
                             to_z_z_zz_xyyy, \
                             to_z_z_zz_xyyz, \
                             to_z_z_zz_xyzz, \
                             to_z_z_zz_xzzz, \
                             to_z_z_zz_yyyy, \
                             to_z_z_zz_yyyz, \
                             to_z_z_zz_yyzz, \
                             to_z_z_zz_yzzz, \
                             to_z_z_zz_zzzz, \
                             to_z_zzz,       \
                             to_z_zzzzz,     \
                             to_zzz_xxx,     \
                             to_zzz_xxxxz,   \
                             to_zzz_xxxyz,   \
                             to_zzz_xxxzz,   \
                             to_zzz_xxy,     \
                             to_zzz_xxyyz,   \
                             to_zzz_xxyzz,   \
                             to_zzz_xxz,     \
                             to_zzz_xxzzz,   \
                             to_zzz_xyy,     \
                             to_zzz_xyyyz,   \
                             to_zzz_xyyzz,   \
                             to_zzz_xyz,     \
                             to_zzz_xyzzz,   \
                             to_zzz_xzz,     \
                             to_zzz_xzzzz,   \
                             to_zzz_yyy,     \
                             to_zzz_yyyyz,   \
                             to_zzz_yyyzz,   \
                             to_zzz_yyz,     \
                             to_zzz_yyzzz,   \
                             to_zzz_yzz,     \
                             to_zzz_yzzzz,   \
                             to_zzz_zzz,     \
                             to_zzz_zzzzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zz_xxxx[k] = -4.0 * to_z_xxxxz[k] * tke_0 + 4.0 * to_zzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxxy[k] = -4.0 * to_z_xxxyz[k] * tke_0 + 4.0 * to_zzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxxz[k] = 2.0 * to_z_xxx[k] - 4.0 * to_z_xxxzz[k] * tke_0 - 2.0 * to_zzz_xxx[k] * tbe_0 + 4.0 * to_zzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxyy[k] = -4.0 * to_z_xxyyz[k] * tke_0 + 4.0 * to_zzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxyz[k] = 2.0 * to_z_xxy[k] - 4.0 * to_z_xxyzz[k] * tke_0 - 2.0 * to_zzz_xxy[k] * tbe_0 + 4.0 * to_zzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxzz[k] = 4.0 * to_z_xxz[k] - 4.0 * to_z_xxzzz[k] * tke_0 - 4.0 * to_zzz_xxz[k] * tbe_0 + 4.0 * to_zzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xyyy[k] = -4.0 * to_z_xyyyz[k] * tke_0 + 4.0 * to_zzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xyyz[k] = 2.0 * to_z_xyy[k] - 4.0 * to_z_xyyzz[k] * tke_0 - 2.0 * to_zzz_xyy[k] * tbe_0 + 4.0 * to_zzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xyzz[k] = 4.0 * to_z_xyz[k] - 4.0 * to_z_xyzzz[k] * tke_0 - 4.0 * to_zzz_xyz[k] * tbe_0 + 4.0 * to_zzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xzzz[k] = 6.0 * to_z_xzz[k] - 4.0 * to_z_xzzzz[k] * tke_0 - 6.0 * to_zzz_xzz[k] * tbe_0 + 4.0 * to_zzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yyyy[k] = -4.0 * to_z_yyyyz[k] * tke_0 + 4.0 * to_zzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_zz_yyyz[k] = 2.0 * to_z_yyy[k] - 4.0 * to_z_yyyzz[k] * tke_0 - 2.0 * to_zzz_yyy[k] * tbe_0 + 4.0 * to_zzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yyzz[k] = 4.0 * to_z_yyz[k] - 4.0 * to_z_yyzzz[k] * tke_0 - 4.0 * to_zzz_yyz[k] * tbe_0 + 4.0 * to_zzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yzzz[k] = 6.0 * to_z_yzz[k] - 4.0 * to_z_yzzzz[k] * tke_0 - 6.0 * to_zzz_yzz[k] * tbe_0 + 4.0 * to_zzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_zzzz[k] = 8.0 * to_z_zzz[k] - 4.0 * to_z_zzzzz[k] * tke_0 - 8.0 * to_zzz_zzz[k] * tbe_0 + 4.0 * to_zzz_zzzzz[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
