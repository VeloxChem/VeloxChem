#include "GeometricalDerivatives0X1ForFG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_fg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_fg,
                       const int idx_op_ff,
                       const int idx_op_fh,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up 0-15 components of targeted buffer : FG

        auto to_0_x_xxx_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 0);

        auto to_0_x_xxx_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 1);

        auto to_0_x_xxx_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 2);

        auto to_0_x_xxx_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 3);

        auto to_0_x_xxx_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 4);

        auto to_0_x_xxx_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 5);

        auto to_0_x_xxx_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 6);

        auto to_0_x_xxx_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 7);

        auto to_0_x_xxx_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 8);

        auto to_0_x_xxx_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 9);

        auto to_0_x_xxx_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 10);

        auto to_0_x_xxx_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 11);

        auto to_0_x_xxx_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 12);

        auto to_0_x_xxx_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 13);

        auto to_0_x_xxx_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_0_x_xxx_xxxx, to_0_x_xxx_xxxy, to_0_x_xxx_xxxz, to_0_x_xxx_xxyy, to_0_x_xxx_xxyz, to_0_x_xxx_xxzz, to_0_x_xxx_xyyy, to_0_x_xxx_xyyz, to_0_x_xxx_xyzz, to_0_x_xxx_xzzz, to_0_x_xxx_yyyy, to_0_x_xxx_yyyz, to_0_x_xxx_yyzz, to_0_x_xxx_yzzz, to_0_x_xxx_zzzz, to_xxx_xxx, to_xxx_xxxxx, to_xxx_xxxxy, to_xxx_xxxxz, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxx_xxxx[k] = -4.0 * to_xxx_xxx[k] + 2.0 * to_xxx_xxxxx[k] * tke_0;

            to_0_x_xxx_xxxy[k] = -3.0 * to_xxx_xxy[k] + 2.0 * to_xxx_xxxxy[k] * tke_0;

            to_0_x_xxx_xxxz[k] = -3.0 * to_xxx_xxz[k] + 2.0 * to_xxx_xxxxz[k] * tke_0;

            to_0_x_xxx_xxyy[k] = -2.0 * to_xxx_xyy[k] + 2.0 * to_xxx_xxxyy[k] * tke_0;

            to_0_x_xxx_xxyz[k] = -2.0 * to_xxx_xyz[k] + 2.0 * to_xxx_xxxyz[k] * tke_0;

            to_0_x_xxx_xxzz[k] = -2.0 * to_xxx_xzz[k] + 2.0 * to_xxx_xxxzz[k] * tke_0;

            to_0_x_xxx_xyyy[k] = -to_xxx_yyy[k] + 2.0 * to_xxx_xxyyy[k] * tke_0;

            to_0_x_xxx_xyyz[k] = -to_xxx_yyz[k] + 2.0 * to_xxx_xxyyz[k] * tke_0;

            to_0_x_xxx_xyzz[k] = -to_xxx_yzz[k] + 2.0 * to_xxx_xxyzz[k] * tke_0;

            to_0_x_xxx_xzzz[k] = -to_xxx_zzz[k] + 2.0 * to_xxx_xxzzz[k] * tke_0;

            to_0_x_xxx_yyyy[k] = 2.0 * to_xxx_xyyyy[k] * tke_0;

            to_0_x_xxx_yyyz[k] = 2.0 * to_xxx_xyyyz[k] * tke_0;

            to_0_x_xxx_yyzz[k] = 2.0 * to_xxx_xyyzz[k] * tke_0;

            to_0_x_xxx_yzzz[k] = 2.0 * to_xxx_xyzzz[k] * tke_0;

            to_0_x_xxx_zzzz[k] = 2.0 * to_xxx_xzzzz[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : FG

        auto to_0_x_xxy_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 15);

        auto to_0_x_xxy_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 16);

        auto to_0_x_xxy_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 17);

        auto to_0_x_xxy_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 18);

        auto to_0_x_xxy_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 19);

        auto to_0_x_xxy_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 20);

        auto to_0_x_xxy_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 21);

        auto to_0_x_xxy_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 22);

        auto to_0_x_xxy_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 23);

        auto to_0_x_xxy_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 24);

        auto to_0_x_xxy_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 25);

        auto to_0_x_xxy_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 26);

        auto to_0_x_xxy_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 27);

        auto to_0_x_xxy_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 28);

        auto to_0_x_xxy_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_x_xxy_xxxx, to_0_x_xxy_xxxy, to_0_x_xxy_xxxz, to_0_x_xxy_xxyy, to_0_x_xxy_xxyz, to_0_x_xxy_xxzz, to_0_x_xxy_xyyy, to_0_x_xxy_xyyz, to_0_x_xxy_xyzz, to_0_x_xxy_xzzz, to_0_x_xxy_yyyy, to_0_x_xxy_yyyz, to_0_x_xxy_yyzz, to_0_x_xxy_yzzz, to_0_x_xxy_zzzz, to_xxy_xxx, to_xxy_xxxxx, to_xxy_xxxxy, to_xxy_xxxxz, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxy_xxxx[k] = -4.0 * to_xxy_xxx[k] + 2.0 * to_xxy_xxxxx[k] * tke_0;

            to_0_x_xxy_xxxy[k] = -3.0 * to_xxy_xxy[k] + 2.0 * to_xxy_xxxxy[k] * tke_0;

            to_0_x_xxy_xxxz[k] = -3.0 * to_xxy_xxz[k] + 2.0 * to_xxy_xxxxz[k] * tke_0;

            to_0_x_xxy_xxyy[k] = -2.0 * to_xxy_xyy[k] + 2.0 * to_xxy_xxxyy[k] * tke_0;

            to_0_x_xxy_xxyz[k] = -2.0 * to_xxy_xyz[k] + 2.0 * to_xxy_xxxyz[k] * tke_0;

            to_0_x_xxy_xxzz[k] = -2.0 * to_xxy_xzz[k] + 2.0 * to_xxy_xxxzz[k] * tke_0;

            to_0_x_xxy_xyyy[k] = -to_xxy_yyy[k] + 2.0 * to_xxy_xxyyy[k] * tke_0;

            to_0_x_xxy_xyyz[k] = -to_xxy_yyz[k] + 2.0 * to_xxy_xxyyz[k] * tke_0;

            to_0_x_xxy_xyzz[k] = -to_xxy_yzz[k] + 2.0 * to_xxy_xxyzz[k] * tke_0;

            to_0_x_xxy_xzzz[k] = -to_xxy_zzz[k] + 2.0 * to_xxy_xxzzz[k] * tke_0;

            to_0_x_xxy_yyyy[k] = 2.0 * to_xxy_xyyyy[k] * tke_0;

            to_0_x_xxy_yyyz[k] = 2.0 * to_xxy_xyyyz[k] * tke_0;

            to_0_x_xxy_yyzz[k] = 2.0 * to_xxy_xyyzz[k] * tke_0;

            to_0_x_xxy_yzzz[k] = 2.0 * to_xxy_xyzzz[k] * tke_0;

            to_0_x_xxy_zzzz[k] = 2.0 * to_xxy_xzzzz[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : FG

        auto to_0_x_xxz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 30);

        auto to_0_x_xxz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 31);

        auto to_0_x_xxz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 32);

        auto to_0_x_xxz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 33);

        auto to_0_x_xxz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 34);

        auto to_0_x_xxz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 35);

        auto to_0_x_xxz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 36);

        auto to_0_x_xxz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 37);

        auto to_0_x_xxz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 38);

        auto to_0_x_xxz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 39);

        auto to_0_x_xxz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 40);

        auto to_0_x_xxz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 41);

        auto to_0_x_xxz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 42);

        auto to_0_x_xxz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 43);

        auto to_0_x_xxz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_0_x_xxz_xxxx, to_0_x_xxz_xxxy, to_0_x_xxz_xxxz, to_0_x_xxz_xxyy, to_0_x_xxz_xxyz, to_0_x_xxz_xxzz, to_0_x_xxz_xyyy, to_0_x_xxz_xyyz, to_0_x_xxz_xyzz, to_0_x_xxz_xzzz, to_0_x_xxz_yyyy, to_0_x_xxz_yyyz, to_0_x_xxz_yyzz, to_0_x_xxz_yzzz, to_0_x_xxz_zzzz, to_xxz_xxx, to_xxz_xxxxx, to_xxz_xxxxy, to_xxz_xxxxz, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxz_xxxx[k] = -4.0 * to_xxz_xxx[k] + 2.0 * to_xxz_xxxxx[k] * tke_0;

            to_0_x_xxz_xxxy[k] = -3.0 * to_xxz_xxy[k] + 2.0 * to_xxz_xxxxy[k] * tke_0;

            to_0_x_xxz_xxxz[k] = -3.0 * to_xxz_xxz[k] + 2.0 * to_xxz_xxxxz[k] * tke_0;

            to_0_x_xxz_xxyy[k] = -2.0 * to_xxz_xyy[k] + 2.0 * to_xxz_xxxyy[k] * tke_0;

            to_0_x_xxz_xxyz[k] = -2.0 * to_xxz_xyz[k] + 2.0 * to_xxz_xxxyz[k] * tke_0;

            to_0_x_xxz_xxzz[k] = -2.0 * to_xxz_xzz[k] + 2.0 * to_xxz_xxxzz[k] * tke_0;

            to_0_x_xxz_xyyy[k] = -to_xxz_yyy[k] + 2.0 * to_xxz_xxyyy[k] * tke_0;

            to_0_x_xxz_xyyz[k] = -to_xxz_yyz[k] + 2.0 * to_xxz_xxyyz[k] * tke_0;

            to_0_x_xxz_xyzz[k] = -to_xxz_yzz[k] + 2.0 * to_xxz_xxyzz[k] * tke_0;

            to_0_x_xxz_xzzz[k] = -to_xxz_zzz[k] + 2.0 * to_xxz_xxzzz[k] * tke_0;

            to_0_x_xxz_yyyy[k] = 2.0 * to_xxz_xyyyy[k] * tke_0;

            to_0_x_xxz_yyyz[k] = 2.0 * to_xxz_xyyyz[k] * tke_0;

            to_0_x_xxz_yyzz[k] = 2.0 * to_xxz_xyyzz[k] * tke_0;

            to_0_x_xxz_yzzz[k] = 2.0 * to_xxz_xyzzz[k] * tke_0;

            to_0_x_xxz_zzzz[k] = 2.0 * to_xxz_xzzzz[k] * tke_0;
        }

        // Set up 45-60 components of targeted buffer : FG

        auto to_0_x_xyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 45);

        auto to_0_x_xyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 46);

        auto to_0_x_xyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 47);

        auto to_0_x_xyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 48);

        auto to_0_x_xyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 49);

        auto to_0_x_xyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 50);

        auto to_0_x_xyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 51);

        auto to_0_x_xyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 52);

        auto to_0_x_xyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 53);

        auto to_0_x_xyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 54);

        auto to_0_x_xyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 55);

        auto to_0_x_xyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 56);

        auto to_0_x_xyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 57);

        auto to_0_x_xyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 58);

        auto to_0_x_xyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_x_xyy_xxxx, to_0_x_xyy_xxxy, to_0_x_xyy_xxxz, to_0_x_xyy_xxyy, to_0_x_xyy_xxyz, to_0_x_xyy_xxzz, to_0_x_xyy_xyyy, to_0_x_xyy_xyyz, to_0_x_xyy_xyzz, to_0_x_xyy_xzzz, to_0_x_xyy_yyyy, to_0_x_xyy_yyyz, to_0_x_xyy_yyzz, to_0_x_xyy_yzzz, to_0_x_xyy_zzzz, to_xyy_xxx, to_xyy_xxxxx, to_xyy_xxxxy, to_xyy_xxxxz, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyy_xxxx[k] = -4.0 * to_xyy_xxx[k] + 2.0 * to_xyy_xxxxx[k] * tke_0;

            to_0_x_xyy_xxxy[k] = -3.0 * to_xyy_xxy[k] + 2.0 * to_xyy_xxxxy[k] * tke_0;

            to_0_x_xyy_xxxz[k] = -3.0 * to_xyy_xxz[k] + 2.0 * to_xyy_xxxxz[k] * tke_0;

            to_0_x_xyy_xxyy[k] = -2.0 * to_xyy_xyy[k] + 2.0 * to_xyy_xxxyy[k] * tke_0;

            to_0_x_xyy_xxyz[k] = -2.0 * to_xyy_xyz[k] + 2.0 * to_xyy_xxxyz[k] * tke_0;

            to_0_x_xyy_xxzz[k] = -2.0 * to_xyy_xzz[k] + 2.0 * to_xyy_xxxzz[k] * tke_0;

            to_0_x_xyy_xyyy[k] = -to_xyy_yyy[k] + 2.0 * to_xyy_xxyyy[k] * tke_0;

            to_0_x_xyy_xyyz[k] = -to_xyy_yyz[k] + 2.0 * to_xyy_xxyyz[k] * tke_0;

            to_0_x_xyy_xyzz[k] = -to_xyy_yzz[k] + 2.0 * to_xyy_xxyzz[k] * tke_0;

            to_0_x_xyy_xzzz[k] = -to_xyy_zzz[k] + 2.0 * to_xyy_xxzzz[k] * tke_0;

            to_0_x_xyy_yyyy[k] = 2.0 * to_xyy_xyyyy[k] * tke_0;

            to_0_x_xyy_yyyz[k] = 2.0 * to_xyy_xyyyz[k] * tke_0;

            to_0_x_xyy_yyzz[k] = 2.0 * to_xyy_xyyzz[k] * tke_0;

            to_0_x_xyy_yzzz[k] = 2.0 * to_xyy_xyzzz[k] * tke_0;

            to_0_x_xyy_zzzz[k] = 2.0 * to_xyy_xzzzz[k] * tke_0;
        }

        // Set up 60-75 components of targeted buffer : FG

        auto to_0_x_xyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 60);

        auto to_0_x_xyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 61);

        auto to_0_x_xyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 62);

        auto to_0_x_xyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 63);

        auto to_0_x_xyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 64);

        auto to_0_x_xyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 65);

        auto to_0_x_xyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 66);

        auto to_0_x_xyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 67);

        auto to_0_x_xyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 68);

        auto to_0_x_xyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 69);

        auto to_0_x_xyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 70);

        auto to_0_x_xyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 71);

        auto to_0_x_xyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 72);

        auto to_0_x_xyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 73);

        auto to_0_x_xyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_0_x_xyz_xxxx, to_0_x_xyz_xxxy, to_0_x_xyz_xxxz, to_0_x_xyz_xxyy, to_0_x_xyz_xxyz, to_0_x_xyz_xxzz, to_0_x_xyz_xyyy, to_0_x_xyz_xyyz, to_0_x_xyz_xyzz, to_0_x_xyz_xzzz, to_0_x_xyz_yyyy, to_0_x_xyz_yyyz, to_0_x_xyz_yyzz, to_0_x_xyz_yzzz, to_0_x_xyz_zzzz, to_xyz_xxx, to_xyz_xxxxx, to_xyz_xxxxy, to_xyz_xxxxz, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyz_xxxx[k] = -4.0 * to_xyz_xxx[k] + 2.0 * to_xyz_xxxxx[k] * tke_0;

            to_0_x_xyz_xxxy[k] = -3.0 * to_xyz_xxy[k] + 2.0 * to_xyz_xxxxy[k] * tke_0;

            to_0_x_xyz_xxxz[k] = -3.0 * to_xyz_xxz[k] + 2.0 * to_xyz_xxxxz[k] * tke_0;

            to_0_x_xyz_xxyy[k] = -2.0 * to_xyz_xyy[k] + 2.0 * to_xyz_xxxyy[k] * tke_0;

            to_0_x_xyz_xxyz[k] = -2.0 * to_xyz_xyz[k] + 2.0 * to_xyz_xxxyz[k] * tke_0;

            to_0_x_xyz_xxzz[k] = -2.0 * to_xyz_xzz[k] + 2.0 * to_xyz_xxxzz[k] * tke_0;

            to_0_x_xyz_xyyy[k] = -to_xyz_yyy[k] + 2.0 * to_xyz_xxyyy[k] * tke_0;

            to_0_x_xyz_xyyz[k] = -to_xyz_yyz[k] + 2.0 * to_xyz_xxyyz[k] * tke_0;

            to_0_x_xyz_xyzz[k] = -to_xyz_yzz[k] + 2.0 * to_xyz_xxyzz[k] * tke_0;

            to_0_x_xyz_xzzz[k] = -to_xyz_zzz[k] + 2.0 * to_xyz_xxzzz[k] * tke_0;

            to_0_x_xyz_yyyy[k] = 2.0 * to_xyz_xyyyy[k] * tke_0;

            to_0_x_xyz_yyyz[k] = 2.0 * to_xyz_xyyyz[k] * tke_0;

            to_0_x_xyz_yyzz[k] = 2.0 * to_xyz_xyyzz[k] * tke_0;

            to_0_x_xyz_yzzz[k] = 2.0 * to_xyz_xyzzz[k] * tke_0;

            to_0_x_xyz_zzzz[k] = 2.0 * to_xyz_xzzzz[k] * tke_0;
        }

        // Set up 75-90 components of targeted buffer : FG

        auto to_0_x_xzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 75);

        auto to_0_x_xzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 76);

        auto to_0_x_xzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 77);

        auto to_0_x_xzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 78);

        auto to_0_x_xzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 79);

        auto to_0_x_xzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 80);

        auto to_0_x_xzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 81);

        auto to_0_x_xzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 82);

        auto to_0_x_xzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 83);

        auto to_0_x_xzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 84);

        auto to_0_x_xzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 85);

        auto to_0_x_xzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 86);

        auto to_0_x_xzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 87);

        auto to_0_x_xzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 88);

        auto to_0_x_xzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_x_xzz_xxxx, to_0_x_xzz_xxxy, to_0_x_xzz_xxxz, to_0_x_xzz_xxyy, to_0_x_xzz_xxyz, to_0_x_xzz_xxzz, to_0_x_xzz_xyyy, to_0_x_xzz_xyyz, to_0_x_xzz_xyzz, to_0_x_xzz_xzzz, to_0_x_xzz_yyyy, to_0_x_xzz_yyyz, to_0_x_xzz_yyzz, to_0_x_xzz_yzzz, to_0_x_xzz_zzzz, to_xzz_xxx, to_xzz_xxxxx, to_xzz_xxxxy, to_xzz_xxxxz, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzz_xxxx[k] = -4.0 * to_xzz_xxx[k] + 2.0 * to_xzz_xxxxx[k] * tke_0;

            to_0_x_xzz_xxxy[k] = -3.0 * to_xzz_xxy[k] + 2.0 * to_xzz_xxxxy[k] * tke_0;

            to_0_x_xzz_xxxz[k] = -3.0 * to_xzz_xxz[k] + 2.0 * to_xzz_xxxxz[k] * tke_0;

            to_0_x_xzz_xxyy[k] = -2.0 * to_xzz_xyy[k] + 2.0 * to_xzz_xxxyy[k] * tke_0;

            to_0_x_xzz_xxyz[k] = -2.0 * to_xzz_xyz[k] + 2.0 * to_xzz_xxxyz[k] * tke_0;

            to_0_x_xzz_xxzz[k] = -2.0 * to_xzz_xzz[k] + 2.0 * to_xzz_xxxzz[k] * tke_0;

            to_0_x_xzz_xyyy[k] = -to_xzz_yyy[k] + 2.0 * to_xzz_xxyyy[k] * tke_0;

            to_0_x_xzz_xyyz[k] = -to_xzz_yyz[k] + 2.0 * to_xzz_xxyyz[k] * tke_0;

            to_0_x_xzz_xyzz[k] = -to_xzz_yzz[k] + 2.0 * to_xzz_xxyzz[k] * tke_0;

            to_0_x_xzz_xzzz[k] = -to_xzz_zzz[k] + 2.0 * to_xzz_xxzzz[k] * tke_0;

            to_0_x_xzz_yyyy[k] = 2.0 * to_xzz_xyyyy[k] * tke_0;

            to_0_x_xzz_yyyz[k] = 2.0 * to_xzz_xyyyz[k] * tke_0;

            to_0_x_xzz_yyzz[k] = 2.0 * to_xzz_xyyzz[k] * tke_0;

            to_0_x_xzz_yzzz[k] = 2.0 * to_xzz_xyzzz[k] * tke_0;

            to_0_x_xzz_zzzz[k] = 2.0 * to_xzz_xzzzz[k] * tke_0;
        }

        // Set up 90-105 components of targeted buffer : FG

        auto to_0_x_yyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 90);

        auto to_0_x_yyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 91);

        auto to_0_x_yyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 92);

        auto to_0_x_yyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 93);

        auto to_0_x_yyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 94);

        auto to_0_x_yyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 95);

        auto to_0_x_yyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 96);

        auto to_0_x_yyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 97);

        auto to_0_x_yyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 98);

        auto to_0_x_yyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 99);

        auto to_0_x_yyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 100);

        auto to_0_x_yyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 101);

        auto to_0_x_yyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 102);

        auto to_0_x_yyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 103);

        auto to_0_x_yyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_0_x_yyy_xxxx, to_0_x_yyy_xxxy, to_0_x_yyy_xxxz, to_0_x_yyy_xxyy, to_0_x_yyy_xxyz, to_0_x_yyy_xxzz, to_0_x_yyy_xyyy, to_0_x_yyy_xyyz, to_0_x_yyy_xyzz, to_0_x_yyy_xzzz, to_0_x_yyy_yyyy, to_0_x_yyy_yyyz, to_0_x_yyy_yyzz, to_0_x_yyy_yzzz, to_0_x_yyy_zzzz, to_yyy_xxx, to_yyy_xxxxx, to_yyy_xxxxy, to_yyy_xxxxz, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyy_xxxx[k] = -4.0 * to_yyy_xxx[k] + 2.0 * to_yyy_xxxxx[k] * tke_0;

            to_0_x_yyy_xxxy[k] = -3.0 * to_yyy_xxy[k] + 2.0 * to_yyy_xxxxy[k] * tke_0;

            to_0_x_yyy_xxxz[k] = -3.0 * to_yyy_xxz[k] + 2.0 * to_yyy_xxxxz[k] * tke_0;

            to_0_x_yyy_xxyy[k] = -2.0 * to_yyy_xyy[k] + 2.0 * to_yyy_xxxyy[k] * tke_0;

            to_0_x_yyy_xxyz[k] = -2.0 * to_yyy_xyz[k] + 2.0 * to_yyy_xxxyz[k] * tke_0;

            to_0_x_yyy_xxzz[k] = -2.0 * to_yyy_xzz[k] + 2.0 * to_yyy_xxxzz[k] * tke_0;

            to_0_x_yyy_xyyy[k] = -to_yyy_yyy[k] + 2.0 * to_yyy_xxyyy[k] * tke_0;

            to_0_x_yyy_xyyz[k] = -to_yyy_yyz[k] + 2.0 * to_yyy_xxyyz[k] * tke_0;

            to_0_x_yyy_xyzz[k] = -to_yyy_yzz[k] + 2.0 * to_yyy_xxyzz[k] * tke_0;

            to_0_x_yyy_xzzz[k] = -to_yyy_zzz[k] + 2.0 * to_yyy_xxzzz[k] * tke_0;

            to_0_x_yyy_yyyy[k] = 2.0 * to_yyy_xyyyy[k] * tke_0;

            to_0_x_yyy_yyyz[k] = 2.0 * to_yyy_xyyyz[k] * tke_0;

            to_0_x_yyy_yyzz[k] = 2.0 * to_yyy_xyyzz[k] * tke_0;

            to_0_x_yyy_yzzz[k] = 2.0 * to_yyy_xyzzz[k] * tke_0;

            to_0_x_yyy_zzzz[k] = 2.0 * to_yyy_xzzzz[k] * tke_0;
        }

        // Set up 105-120 components of targeted buffer : FG

        auto to_0_x_yyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 105);

        auto to_0_x_yyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 106);

        auto to_0_x_yyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 107);

        auto to_0_x_yyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 108);

        auto to_0_x_yyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 109);

        auto to_0_x_yyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 110);

        auto to_0_x_yyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 111);

        auto to_0_x_yyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 112);

        auto to_0_x_yyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 113);

        auto to_0_x_yyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 114);

        auto to_0_x_yyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 115);

        auto to_0_x_yyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 116);

        auto to_0_x_yyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 117);

        auto to_0_x_yyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 118);

        auto to_0_x_yyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_x_yyz_xxxx, to_0_x_yyz_xxxy, to_0_x_yyz_xxxz, to_0_x_yyz_xxyy, to_0_x_yyz_xxyz, to_0_x_yyz_xxzz, to_0_x_yyz_xyyy, to_0_x_yyz_xyyz, to_0_x_yyz_xyzz, to_0_x_yyz_xzzz, to_0_x_yyz_yyyy, to_0_x_yyz_yyyz, to_0_x_yyz_yyzz, to_0_x_yyz_yzzz, to_0_x_yyz_zzzz, to_yyz_xxx, to_yyz_xxxxx, to_yyz_xxxxy, to_yyz_xxxxz, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyz_xxxx[k] = -4.0 * to_yyz_xxx[k] + 2.0 * to_yyz_xxxxx[k] * tke_0;

            to_0_x_yyz_xxxy[k] = -3.0 * to_yyz_xxy[k] + 2.0 * to_yyz_xxxxy[k] * tke_0;

            to_0_x_yyz_xxxz[k] = -3.0 * to_yyz_xxz[k] + 2.0 * to_yyz_xxxxz[k] * tke_0;

            to_0_x_yyz_xxyy[k] = -2.0 * to_yyz_xyy[k] + 2.0 * to_yyz_xxxyy[k] * tke_0;

            to_0_x_yyz_xxyz[k] = -2.0 * to_yyz_xyz[k] + 2.0 * to_yyz_xxxyz[k] * tke_0;

            to_0_x_yyz_xxzz[k] = -2.0 * to_yyz_xzz[k] + 2.0 * to_yyz_xxxzz[k] * tke_0;

            to_0_x_yyz_xyyy[k] = -to_yyz_yyy[k] + 2.0 * to_yyz_xxyyy[k] * tke_0;

            to_0_x_yyz_xyyz[k] = -to_yyz_yyz[k] + 2.0 * to_yyz_xxyyz[k] * tke_0;

            to_0_x_yyz_xyzz[k] = -to_yyz_yzz[k] + 2.0 * to_yyz_xxyzz[k] * tke_0;

            to_0_x_yyz_xzzz[k] = -to_yyz_zzz[k] + 2.0 * to_yyz_xxzzz[k] * tke_0;

            to_0_x_yyz_yyyy[k] = 2.0 * to_yyz_xyyyy[k] * tke_0;

            to_0_x_yyz_yyyz[k] = 2.0 * to_yyz_xyyyz[k] * tke_0;

            to_0_x_yyz_yyzz[k] = 2.0 * to_yyz_xyyzz[k] * tke_0;

            to_0_x_yyz_yzzz[k] = 2.0 * to_yyz_xyzzz[k] * tke_0;

            to_0_x_yyz_zzzz[k] = 2.0 * to_yyz_xzzzz[k] * tke_0;
        }

        // Set up 120-135 components of targeted buffer : FG

        auto to_0_x_yzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 120);

        auto to_0_x_yzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 121);

        auto to_0_x_yzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 122);

        auto to_0_x_yzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 123);

        auto to_0_x_yzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 124);

        auto to_0_x_yzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 125);

        auto to_0_x_yzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 126);

        auto to_0_x_yzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 127);

        auto to_0_x_yzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 128);

        auto to_0_x_yzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 129);

        auto to_0_x_yzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 130);

        auto to_0_x_yzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 131);

        auto to_0_x_yzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 132);

        auto to_0_x_yzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 133);

        auto to_0_x_yzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_0_x_yzz_xxxx, to_0_x_yzz_xxxy, to_0_x_yzz_xxxz, to_0_x_yzz_xxyy, to_0_x_yzz_xxyz, to_0_x_yzz_xxzz, to_0_x_yzz_xyyy, to_0_x_yzz_xyyz, to_0_x_yzz_xyzz, to_0_x_yzz_xzzz, to_0_x_yzz_yyyy, to_0_x_yzz_yyyz, to_0_x_yzz_yyzz, to_0_x_yzz_yzzz, to_0_x_yzz_zzzz, to_yzz_xxx, to_yzz_xxxxx, to_yzz_xxxxy, to_yzz_xxxxz, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzz_xxxx[k] = -4.0 * to_yzz_xxx[k] + 2.0 * to_yzz_xxxxx[k] * tke_0;

            to_0_x_yzz_xxxy[k] = -3.0 * to_yzz_xxy[k] + 2.0 * to_yzz_xxxxy[k] * tke_0;

            to_0_x_yzz_xxxz[k] = -3.0 * to_yzz_xxz[k] + 2.0 * to_yzz_xxxxz[k] * tke_0;

            to_0_x_yzz_xxyy[k] = -2.0 * to_yzz_xyy[k] + 2.0 * to_yzz_xxxyy[k] * tke_0;

            to_0_x_yzz_xxyz[k] = -2.0 * to_yzz_xyz[k] + 2.0 * to_yzz_xxxyz[k] * tke_0;

            to_0_x_yzz_xxzz[k] = -2.0 * to_yzz_xzz[k] + 2.0 * to_yzz_xxxzz[k] * tke_0;

            to_0_x_yzz_xyyy[k] = -to_yzz_yyy[k] + 2.0 * to_yzz_xxyyy[k] * tke_0;

            to_0_x_yzz_xyyz[k] = -to_yzz_yyz[k] + 2.0 * to_yzz_xxyyz[k] * tke_0;

            to_0_x_yzz_xyzz[k] = -to_yzz_yzz[k] + 2.0 * to_yzz_xxyzz[k] * tke_0;

            to_0_x_yzz_xzzz[k] = -to_yzz_zzz[k] + 2.0 * to_yzz_xxzzz[k] * tke_0;

            to_0_x_yzz_yyyy[k] = 2.0 * to_yzz_xyyyy[k] * tke_0;

            to_0_x_yzz_yyyz[k] = 2.0 * to_yzz_xyyyz[k] * tke_0;

            to_0_x_yzz_yyzz[k] = 2.0 * to_yzz_xyyzz[k] * tke_0;

            to_0_x_yzz_yzzz[k] = 2.0 * to_yzz_xyzzz[k] * tke_0;

            to_0_x_yzz_zzzz[k] = 2.0 * to_yzz_xzzzz[k] * tke_0;
        }

        // Set up 135-150 components of targeted buffer : FG

        auto to_0_x_zzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 135);

        auto to_0_x_zzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 136);

        auto to_0_x_zzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 137);

        auto to_0_x_zzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 138);

        auto to_0_x_zzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 139);

        auto to_0_x_zzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 140);

        auto to_0_x_zzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 141);

        auto to_0_x_zzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 142);

        auto to_0_x_zzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 143);

        auto to_0_x_zzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 144);

        auto to_0_x_zzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 145);

        auto to_0_x_zzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 146);

        auto to_0_x_zzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 147);

        auto to_0_x_zzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 148);

        auto to_0_x_zzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 0 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_x_zzz_xxxx, to_0_x_zzz_xxxy, to_0_x_zzz_xxxz, to_0_x_zzz_xxyy, to_0_x_zzz_xxyz, to_0_x_zzz_xxzz, to_0_x_zzz_xyyy, to_0_x_zzz_xyyz, to_0_x_zzz_xyzz, to_0_x_zzz_xzzz, to_0_x_zzz_yyyy, to_0_x_zzz_yyyz, to_0_x_zzz_yyzz, to_0_x_zzz_yzzz, to_0_x_zzz_zzzz, to_zzz_xxx, to_zzz_xxxxx, to_zzz_xxxxy, to_zzz_xxxxz, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzz_xxxx[k] = -4.0 * to_zzz_xxx[k] + 2.0 * to_zzz_xxxxx[k] * tke_0;

            to_0_x_zzz_xxxy[k] = -3.0 * to_zzz_xxy[k] + 2.0 * to_zzz_xxxxy[k] * tke_0;

            to_0_x_zzz_xxxz[k] = -3.0 * to_zzz_xxz[k] + 2.0 * to_zzz_xxxxz[k] * tke_0;

            to_0_x_zzz_xxyy[k] = -2.0 * to_zzz_xyy[k] + 2.0 * to_zzz_xxxyy[k] * tke_0;

            to_0_x_zzz_xxyz[k] = -2.0 * to_zzz_xyz[k] + 2.0 * to_zzz_xxxyz[k] * tke_0;

            to_0_x_zzz_xxzz[k] = -2.0 * to_zzz_xzz[k] + 2.0 * to_zzz_xxxzz[k] * tke_0;

            to_0_x_zzz_xyyy[k] = -to_zzz_yyy[k] + 2.0 * to_zzz_xxyyy[k] * tke_0;

            to_0_x_zzz_xyyz[k] = -to_zzz_yyz[k] + 2.0 * to_zzz_xxyyz[k] * tke_0;

            to_0_x_zzz_xyzz[k] = -to_zzz_yzz[k] + 2.0 * to_zzz_xxyzz[k] * tke_0;

            to_0_x_zzz_xzzz[k] = -to_zzz_zzz[k] + 2.0 * to_zzz_xxzzz[k] * tke_0;

            to_0_x_zzz_yyyy[k] = 2.0 * to_zzz_xyyyy[k] * tke_0;

            to_0_x_zzz_yyyz[k] = 2.0 * to_zzz_xyyyz[k] * tke_0;

            to_0_x_zzz_yyzz[k] = 2.0 * to_zzz_xyyzz[k] * tke_0;

            to_0_x_zzz_yzzz[k] = 2.0 * to_zzz_xyzzz[k] * tke_0;

            to_0_x_zzz_zzzz[k] = 2.0 * to_zzz_xzzzz[k] * tke_0;
        }

        // Set up 150-165 components of targeted buffer : FG

        auto to_0_y_xxx_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 0);

        auto to_0_y_xxx_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 1);

        auto to_0_y_xxx_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 2);

        auto to_0_y_xxx_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 3);

        auto to_0_y_xxx_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 4);

        auto to_0_y_xxx_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 5);

        auto to_0_y_xxx_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 6);

        auto to_0_y_xxx_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 7);

        auto to_0_y_xxx_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 8);

        auto to_0_y_xxx_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 9);

        auto to_0_y_xxx_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 10);

        auto to_0_y_xxx_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 11);

        auto to_0_y_xxx_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 12);

        auto to_0_y_xxx_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 13);

        auto to_0_y_xxx_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_0_y_xxx_xxxx, to_0_y_xxx_xxxy, to_0_y_xxx_xxxz, to_0_y_xxx_xxyy, to_0_y_xxx_xxyz, to_0_y_xxx_xxzz, to_0_y_xxx_xyyy, to_0_y_xxx_xyyz, to_0_y_xxx_xyzz, to_0_y_xxx_xzzz, to_0_y_xxx_yyyy, to_0_y_xxx_yyyz, to_0_y_xxx_yyzz, to_0_y_xxx_yzzz, to_0_y_xxx_zzzz, to_xxx_xxx, to_xxx_xxxxy, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_yyy, to_xxx_yyyyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxx_xxxx[k] = 2.0 * to_xxx_xxxxy[k] * tke_0;

            to_0_y_xxx_xxxy[k] = -to_xxx_xxx[k] + 2.0 * to_xxx_xxxyy[k] * tke_0;

            to_0_y_xxx_xxxz[k] = 2.0 * to_xxx_xxxyz[k] * tke_0;

            to_0_y_xxx_xxyy[k] = -2.0 * to_xxx_xxy[k] + 2.0 * to_xxx_xxyyy[k] * tke_0;

            to_0_y_xxx_xxyz[k] = -to_xxx_xxz[k] + 2.0 * to_xxx_xxyyz[k] * tke_0;

            to_0_y_xxx_xxzz[k] = 2.0 * to_xxx_xxyzz[k] * tke_0;

            to_0_y_xxx_xyyy[k] = -3.0 * to_xxx_xyy[k] + 2.0 * to_xxx_xyyyy[k] * tke_0;

            to_0_y_xxx_xyyz[k] = -2.0 * to_xxx_xyz[k] + 2.0 * to_xxx_xyyyz[k] * tke_0;

            to_0_y_xxx_xyzz[k] = -to_xxx_xzz[k] + 2.0 * to_xxx_xyyzz[k] * tke_0;

            to_0_y_xxx_xzzz[k] = 2.0 * to_xxx_xyzzz[k] * tke_0;

            to_0_y_xxx_yyyy[k] = -4.0 * to_xxx_yyy[k] + 2.0 * to_xxx_yyyyy[k] * tke_0;

            to_0_y_xxx_yyyz[k] = -3.0 * to_xxx_yyz[k] + 2.0 * to_xxx_yyyyz[k] * tke_0;

            to_0_y_xxx_yyzz[k] = -2.0 * to_xxx_yzz[k] + 2.0 * to_xxx_yyyzz[k] * tke_0;

            to_0_y_xxx_yzzz[k] = -to_xxx_zzz[k] + 2.0 * to_xxx_yyzzz[k] * tke_0;

            to_0_y_xxx_zzzz[k] = 2.0 * to_xxx_yzzzz[k] * tke_0;
        }

        // Set up 165-180 components of targeted buffer : FG

        auto to_0_y_xxy_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 15);

        auto to_0_y_xxy_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 16);

        auto to_0_y_xxy_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 17);

        auto to_0_y_xxy_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 18);

        auto to_0_y_xxy_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 19);

        auto to_0_y_xxy_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 20);

        auto to_0_y_xxy_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 21);

        auto to_0_y_xxy_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 22);

        auto to_0_y_xxy_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 23);

        auto to_0_y_xxy_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 24);

        auto to_0_y_xxy_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 25);

        auto to_0_y_xxy_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 26);

        auto to_0_y_xxy_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 27);

        auto to_0_y_xxy_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 28);

        auto to_0_y_xxy_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_y_xxy_xxxx, to_0_y_xxy_xxxy, to_0_y_xxy_xxxz, to_0_y_xxy_xxyy, to_0_y_xxy_xxyz, to_0_y_xxy_xxzz, to_0_y_xxy_xyyy, to_0_y_xxy_xyyz, to_0_y_xxy_xyzz, to_0_y_xxy_xzzz, to_0_y_xxy_yyyy, to_0_y_xxy_yyyz, to_0_y_xxy_yyzz, to_0_y_xxy_yzzz, to_0_y_xxy_zzzz, to_xxy_xxx, to_xxy_xxxxy, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_yyy, to_xxy_yyyyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxy_xxxx[k] = 2.0 * to_xxy_xxxxy[k] * tke_0;

            to_0_y_xxy_xxxy[k] = -to_xxy_xxx[k] + 2.0 * to_xxy_xxxyy[k] * tke_0;

            to_0_y_xxy_xxxz[k] = 2.0 * to_xxy_xxxyz[k] * tke_0;

            to_0_y_xxy_xxyy[k] = -2.0 * to_xxy_xxy[k] + 2.0 * to_xxy_xxyyy[k] * tke_0;

            to_0_y_xxy_xxyz[k] = -to_xxy_xxz[k] + 2.0 * to_xxy_xxyyz[k] * tke_0;

            to_0_y_xxy_xxzz[k] = 2.0 * to_xxy_xxyzz[k] * tke_0;

            to_0_y_xxy_xyyy[k] = -3.0 * to_xxy_xyy[k] + 2.0 * to_xxy_xyyyy[k] * tke_0;

            to_0_y_xxy_xyyz[k] = -2.0 * to_xxy_xyz[k] + 2.0 * to_xxy_xyyyz[k] * tke_0;

            to_0_y_xxy_xyzz[k] = -to_xxy_xzz[k] + 2.0 * to_xxy_xyyzz[k] * tke_0;

            to_0_y_xxy_xzzz[k] = 2.0 * to_xxy_xyzzz[k] * tke_0;

            to_0_y_xxy_yyyy[k] = -4.0 * to_xxy_yyy[k] + 2.0 * to_xxy_yyyyy[k] * tke_0;

            to_0_y_xxy_yyyz[k] = -3.0 * to_xxy_yyz[k] + 2.0 * to_xxy_yyyyz[k] * tke_0;

            to_0_y_xxy_yyzz[k] = -2.0 * to_xxy_yzz[k] + 2.0 * to_xxy_yyyzz[k] * tke_0;

            to_0_y_xxy_yzzz[k] = -to_xxy_zzz[k] + 2.0 * to_xxy_yyzzz[k] * tke_0;

            to_0_y_xxy_zzzz[k] = 2.0 * to_xxy_yzzzz[k] * tke_0;
        }

        // Set up 180-195 components of targeted buffer : FG

        auto to_0_y_xxz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 30);

        auto to_0_y_xxz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 31);

        auto to_0_y_xxz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 32);

        auto to_0_y_xxz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 33);

        auto to_0_y_xxz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 34);

        auto to_0_y_xxz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 35);

        auto to_0_y_xxz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 36);

        auto to_0_y_xxz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 37);

        auto to_0_y_xxz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 38);

        auto to_0_y_xxz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 39);

        auto to_0_y_xxz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 40);

        auto to_0_y_xxz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 41);

        auto to_0_y_xxz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 42);

        auto to_0_y_xxz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 43);

        auto to_0_y_xxz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_0_y_xxz_xxxx, to_0_y_xxz_xxxy, to_0_y_xxz_xxxz, to_0_y_xxz_xxyy, to_0_y_xxz_xxyz, to_0_y_xxz_xxzz, to_0_y_xxz_xyyy, to_0_y_xxz_xyyz, to_0_y_xxz_xyzz, to_0_y_xxz_xzzz, to_0_y_xxz_yyyy, to_0_y_xxz_yyyz, to_0_y_xxz_yyzz, to_0_y_xxz_yzzz, to_0_y_xxz_zzzz, to_xxz_xxx, to_xxz_xxxxy, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_yyy, to_xxz_yyyyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxz_xxxx[k] = 2.0 * to_xxz_xxxxy[k] * tke_0;

            to_0_y_xxz_xxxy[k] = -to_xxz_xxx[k] + 2.0 * to_xxz_xxxyy[k] * tke_0;

            to_0_y_xxz_xxxz[k] = 2.0 * to_xxz_xxxyz[k] * tke_0;

            to_0_y_xxz_xxyy[k] = -2.0 * to_xxz_xxy[k] + 2.0 * to_xxz_xxyyy[k] * tke_0;

            to_0_y_xxz_xxyz[k] = -to_xxz_xxz[k] + 2.0 * to_xxz_xxyyz[k] * tke_0;

            to_0_y_xxz_xxzz[k] = 2.0 * to_xxz_xxyzz[k] * tke_0;

            to_0_y_xxz_xyyy[k] = -3.0 * to_xxz_xyy[k] + 2.0 * to_xxz_xyyyy[k] * tke_0;

            to_0_y_xxz_xyyz[k] = -2.0 * to_xxz_xyz[k] + 2.0 * to_xxz_xyyyz[k] * tke_0;

            to_0_y_xxz_xyzz[k] = -to_xxz_xzz[k] + 2.0 * to_xxz_xyyzz[k] * tke_0;

            to_0_y_xxz_xzzz[k] = 2.0 * to_xxz_xyzzz[k] * tke_0;

            to_0_y_xxz_yyyy[k] = -4.0 * to_xxz_yyy[k] + 2.0 * to_xxz_yyyyy[k] * tke_0;

            to_0_y_xxz_yyyz[k] = -3.0 * to_xxz_yyz[k] + 2.0 * to_xxz_yyyyz[k] * tke_0;

            to_0_y_xxz_yyzz[k] = -2.0 * to_xxz_yzz[k] + 2.0 * to_xxz_yyyzz[k] * tke_0;

            to_0_y_xxz_yzzz[k] = -to_xxz_zzz[k] + 2.0 * to_xxz_yyzzz[k] * tke_0;

            to_0_y_xxz_zzzz[k] = 2.0 * to_xxz_yzzzz[k] * tke_0;
        }

        // Set up 195-210 components of targeted buffer : FG

        auto to_0_y_xyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 45);

        auto to_0_y_xyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 46);

        auto to_0_y_xyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 47);

        auto to_0_y_xyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 48);

        auto to_0_y_xyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 49);

        auto to_0_y_xyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 50);

        auto to_0_y_xyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 51);

        auto to_0_y_xyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 52);

        auto to_0_y_xyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 53);

        auto to_0_y_xyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 54);

        auto to_0_y_xyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 55);

        auto to_0_y_xyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 56);

        auto to_0_y_xyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 57);

        auto to_0_y_xyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 58);

        auto to_0_y_xyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_y_xyy_xxxx, to_0_y_xyy_xxxy, to_0_y_xyy_xxxz, to_0_y_xyy_xxyy, to_0_y_xyy_xxyz, to_0_y_xyy_xxzz, to_0_y_xyy_xyyy, to_0_y_xyy_xyyz, to_0_y_xyy_xyzz, to_0_y_xyy_xzzz, to_0_y_xyy_yyyy, to_0_y_xyy_yyyz, to_0_y_xyy_yyzz, to_0_y_xyy_yzzz, to_0_y_xyy_zzzz, to_xyy_xxx, to_xyy_xxxxy, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_yyy, to_xyy_yyyyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyy_xxxx[k] = 2.0 * to_xyy_xxxxy[k] * tke_0;

            to_0_y_xyy_xxxy[k] = -to_xyy_xxx[k] + 2.0 * to_xyy_xxxyy[k] * tke_0;

            to_0_y_xyy_xxxz[k] = 2.0 * to_xyy_xxxyz[k] * tke_0;

            to_0_y_xyy_xxyy[k] = -2.0 * to_xyy_xxy[k] + 2.0 * to_xyy_xxyyy[k] * tke_0;

            to_0_y_xyy_xxyz[k] = -to_xyy_xxz[k] + 2.0 * to_xyy_xxyyz[k] * tke_0;

            to_0_y_xyy_xxzz[k] = 2.0 * to_xyy_xxyzz[k] * tke_0;

            to_0_y_xyy_xyyy[k] = -3.0 * to_xyy_xyy[k] + 2.0 * to_xyy_xyyyy[k] * tke_0;

            to_0_y_xyy_xyyz[k] = -2.0 * to_xyy_xyz[k] + 2.0 * to_xyy_xyyyz[k] * tke_0;

            to_0_y_xyy_xyzz[k] = -to_xyy_xzz[k] + 2.0 * to_xyy_xyyzz[k] * tke_0;

            to_0_y_xyy_xzzz[k] = 2.0 * to_xyy_xyzzz[k] * tke_0;

            to_0_y_xyy_yyyy[k] = -4.0 * to_xyy_yyy[k] + 2.0 * to_xyy_yyyyy[k] * tke_0;

            to_0_y_xyy_yyyz[k] = -3.0 * to_xyy_yyz[k] + 2.0 * to_xyy_yyyyz[k] * tke_0;

            to_0_y_xyy_yyzz[k] = -2.0 * to_xyy_yzz[k] + 2.0 * to_xyy_yyyzz[k] * tke_0;

            to_0_y_xyy_yzzz[k] = -to_xyy_zzz[k] + 2.0 * to_xyy_yyzzz[k] * tke_0;

            to_0_y_xyy_zzzz[k] = 2.0 * to_xyy_yzzzz[k] * tke_0;
        }

        // Set up 210-225 components of targeted buffer : FG

        auto to_0_y_xyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 60);

        auto to_0_y_xyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 61);

        auto to_0_y_xyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 62);

        auto to_0_y_xyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 63);

        auto to_0_y_xyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 64);

        auto to_0_y_xyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 65);

        auto to_0_y_xyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 66);

        auto to_0_y_xyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 67);

        auto to_0_y_xyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 68);

        auto to_0_y_xyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 69);

        auto to_0_y_xyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 70);

        auto to_0_y_xyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 71);

        auto to_0_y_xyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 72);

        auto to_0_y_xyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 73);

        auto to_0_y_xyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_0_y_xyz_xxxx, to_0_y_xyz_xxxy, to_0_y_xyz_xxxz, to_0_y_xyz_xxyy, to_0_y_xyz_xxyz, to_0_y_xyz_xxzz, to_0_y_xyz_xyyy, to_0_y_xyz_xyyz, to_0_y_xyz_xyzz, to_0_y_xyz_xzzz, to_0_y_xyz_yyyy, to_0_y_xyz_yyyz, to_0_y_xyz_yyzz, to_0_y_xyz_yzzz, to_0_y_xyz_zzzz, to_xyz_xxx, to_xyz_xxxxy, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_yyy, to_xyz_yyyyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyz_xxxx[k] = 2.0 * to_xyz_xxxxy[k] * tke_0;

            to_0_y_xyz_xxxy[k] = -to_xyz_xxx[k] + 2.0 * to_xyz_xxxyy[k] * tke_0;

            to_0_y_xyz_xxxz[k] = 2.0 * to_xyz_xxxyz[k] * tke_0;

            to_0_y_xyz_xxyy[k] = -2.0 * to_xyz_xxy[k] + 2.0 * to_xyz_xxyyy[k] * tke_0;

            to_0_y_xyz_xxyz[k] = -to_xyz_xxz[k] + 2.0 * to_xyz_xxyyz[k] * tke_0;

            to_0_y_xyz_xxzz[k] = 2.0 * to_xyz_xxyzz[k] * tke_0;

            to_0_y_xyz_xyyy[k] = -3.0 * to_xyz_xyy[k] + 2.0 * to_xyz_xyyyy[k] * tke_0;

            to_0_y_xyz_xyyz[k] = -2.0 * to_xyz_xyz[k] + 2.0 * to_xyz_xyyyz[k] * tke_0;

            to_0_y_xyz_xyzz[k] = -to_xyz_xzz[k] + 2.0 * to_xyz_xyyzz[k] * tke_0;

            to_0_y_xyz_xzzz[k] = 2.0 * to_xyz_xyzzz[k] * tke_0;

            to_0_y_xyz_yyyy[k] = -4.0 * to_xyz_yyy[k] + 2.0 * to_xyz_yyyyy[k] * tke_0;

            to_0_y_xyz_yyyz[k] = -3.0 * to_xyz_yyz[k] + 2.0 * to_xyz_yyyyz[k] * tke_0;

            to_0_y_xyz_yyzz[k] = -2.0 * to_xyz_yzz[k] + 2.0 * to_xyz_yyyzz[k] * tke_0;

            to_0_y_xyz_yzzz[k] = -to_xyz_zzz[k] + 2.0 * to_xyz_yyzzz[k] * tke_0;

            to_0_y_xyz_zzzz[k] = 2.0 * to_xyz_yzzzz[k] * tke_0;
        }

        // Set up 225-240 components of targeted buffer : FG

        auto to_0_y_xzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 75);

        auto to_0_y_xzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 76);

        auto to_0_y_xzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 77);

        auto to_0_y_xzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 78);

        auto to_0_y_xzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 79);

        auto to_0_y_xzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 80);

        auto to_0_y_xzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 81);

        auto to_0_y_xzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 82);

        auto to_0_y_xzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 83);

        auto to_0_y_xzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 84);

        auto to_0_y_xzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 85);

        auto to_0_y_xzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 86);

        auto to_0_y_xzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 87);

        auto to_0_y_xzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 88);

        auto to_0_y_xzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_y_xzz_xxxx, to_0_y_xzz_xxxy, to_0_y_xzz_xxxz, to_0_y_xzz_xxyy, to_0_y_xzz_xxyz, to_0_y_xzz_xxzz, to_0_y_xzz_xyyy, to_0_y_xzz_xyyz, to_0_y_xzz_xyzz, to_0_y_xzz_xzzz, to_0_y_xzz_yyyy, to_0_y_xzz_yyyz, to_0_y_xzz_yyzz, to_0_y_xzz_yzzz, to_0_y_xzz_zzzz, to_xzz_xxx, to_xzz_xxxxy, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_yyy, to_xzz_yyyyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzz_xxxx[k] = 2.0 * to_xzz_xxxxy[k] * tke_0;

            to_0_y_xzz_xxxy[k] = -to_xzz_xxx[k] + 2.0 * to_xzz_xxxyy[k] * tke_0;

            to_0_y_xzz_xxxz[k] = 2.0 * to_xzz_xxxyz[k] * tke_0;

            to_0_y_xzz_xxyy[k] = -2.0 * to_xzz_xxy[k] + 2.0 * to_xzz_xxyyy[k] * tke_0;

            to_0_y_xzz_xxyz[k] = -to_xzz_xxz[k] + 2.0 * to_xzz_xxyyz[k] * tke_0;

            to_0_y_xzz_xxzz[k] = 2.0 * to_xzz_xxyzz[k] * tke_0;

            to_0_y_xzz_xyyy[k] = -3.0 * to_xzz_xyy[k] + 2.0 * to_xzz_xyyyy[k] * tke_0;

            to_0_y_xzz_xyyz[k] = -2.0 * to_xzz_xyz[k] + 2.0 * to_xzz_xyyyz[k] * tke_0;

            to_0_y_xzz_xyzz[k] = -to_xzz_xzz[k] + 2.0 * to_xzz_xyyzz[k] * tke_0;

            to_0_y_xzz_xzzz[k] = 2.0 * to_xzz_xyzzz[k] * tke_0;

            to_0_y_xzz_yyyy[k] = -4.0 * to_xzz_yyy[k] + 2.0 * to_xzz_yyyyy[k] * tke_0;

            to_0_y_xzz_yyyz[k] = -3.0 * to_xzz_yyz[k] + 2.0 * to_xzz_yyyyz[k] * tke_0;

            to_0_y_xzz_yyzz[k] = -2.0 * to_xzz_yzz[k] + 2.0 * to_xzz_yyyzz[k] * tke_0;

            to_0_y_xzz_yzzz[k] = -to_xzz_zzz[k] + 2.0 * to_xzz_yyzzz[k] * tke_0;

            to_0_y_xzz_zzzz[k] = 2.0 * to_xzz_yzzzz[k] * tke_0;
        }

        // Set up 240-255 components of targeted buffer : FG

        auto to_0_y_yyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 90);

        auto to_0_y_yyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 91);

        auto to_0_y_yyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 92);

        auto to_0_y_yyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 93);

        auto to_0_y_yyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 94);

        auto to_0_y_yyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 95);

        auto to_0_y_yyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 96);

        auto to_0_y_yyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 97);

        auto to_0_y_yyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 98);

        auto to_0_y_yyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 99);

        auto to_0_y_yyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 100);

        auto to_0_y_yyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 101);

        auto to_0_y_yyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 102);

        auto to_0_y_yyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 103);

        auto to_0_y_yyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_0_y_yyy_xxxx, to_0_y_yyy_xxxy, to_0_y_yyy_xxxz, to_0_y_yyy_xxyy, to_0_y_yyy_xxyz, to_0_y_yyy_xxzz, to_0_y_yyy_xyyy, to_0_y_yyy_xyyz, to_0_y_yyy_xyzz, to_0_y_yyy_xzzz, to_0_y_yyy_yyyy, to_0_y_yyy_yyyz, to_0_y_yyy_yyzz, to_0_y_yyy_yzzz, to_0_y_yyy_zzzz, to_yyy_xxx, to_yyy_xxxxy, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_yyy, to_yyy_yyyyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyy_xxxx[k] = 2.0 * to_yyy_xxxxy[k] * tke_0;

            to_0_y_yyy_xxxy[k] = -to_yyy_xxx[k] + 2.0 * to_yyy_xxxyy[k] * tke_0;

            to_0_y_yyy_xxxz[k] = 2.0 * to_yyy_xxxyz[k] * tke_0;

            to_0_y_yyy_xxyy[k] = -2.0 * to_yyy_xxy[k] + 2.0 * to_yyy_xxyyy[k] * tke_0;

            to_0_y_yyy_xxyz[k] = -to_yyy_xxz[k] + 2.0 * to_yyy_xxyyz[k] * tke_0;

            to_0_y_yyy_xxzz[k] = 2.0 * to_yyy_xxyzz[k] * tke_0;

            to_0_y_yyy_xyyy[k] = -3.0 * to_yyy_xyy[k] + 2.0 * to_yyy_xyyyy[k] * tke_0;

            to_0_y_yyy_xyyz[k] = -2.0 * to_yyy_xyz[k] + 2.0 * to_yyy_xyyyz[k] * tke_0;

            to_0_y_yyy_xyzz[k] = -to_yyy_xzz[k] + 2.0 * to_yyy_xyyzz[k] * tke_0;

            to_0_y_yyy_xzzz[k] = 2.0 * to_yyy_xyzzz[k] * tke_0;

            to_0_y_yyy_yyyy[k] = -4.0 * to_yyy_yyy[k] + 2.0 * to_yyy_yyyyy[k] * tke_0;

            to_0_y_yyy_yyyz[k] = -3.0 * to_yyy_yyz[k] + 2.0 * to_yyy_yyyyz[k] * tke_0;

            to_0_y_yyy_yyzz[k] = -2.0 * to_yyy_yzz[k] + 2.0 * to_yyy_yyyzz[k] * tke_0;

            to_0_y_yyy_yzzz[k] = -to_yyy_zzz[k] + 2.0 * to_yyy_yyzzz[k] * tke_0;

            to_0_y_yyy_zzzz[k] = 2.0 * to_yyy_yzzzz[k] * tke_0;
        }

        // Set up 255-270 components of targeted buffer : FG

        auto to_0_y_yyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 105);

        auto to_0_y_yyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 106);

        auto to_0_y_yyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 107);

        auto to_0_y_yyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 108);

        auto to_0_y_yyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 109);

        auto to_0_y_yyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 110);

        auto to_0_y_yyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 111);

        auto to_0_y_yyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 112);

        auto to_0_y_yyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 113);

        auto to_0_y_yyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 114);

        auto to_0_y_yyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 115);

        auto to_0_y_yyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 116);

        auto to_0_y_yyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 117);

        auto to_0_y_yyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 118);

        auto to_0_y_yyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_y_yyz_xxxx, to_0_y_yyz_xxxy, to_0_y_yyz_xxxz, to_0_y_yyz_xxyy, to_0_y_yyz_xxyz, to_0_y_yyz_xxzz, to_0_y_yyz_xyyy, to_0_y_yyz_xyyz, to_0_y_yyz_xyzz, to_0_y_yyz_xzzz, to_0_y_yyz_yyyy, to_0_y_yyz_yyyz, to_0_y_yyz_yyzz, to_0_y_yyz_yzzz, to_0_y_yyz_zzzz, to_yyz_xxx, to_yyz_xxxxy, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_yyy, to_yyz_yyyyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyz_xxxx[k] = 2.0 * to_yyz_xxxxy[k] * tke_0;

            to_0_y_yyz_xxxy[k] = -to_yyz_xxx[k] + 2.0 * to_yyz_xxxyy[k] * tke_0;

            to_0_y_yyz_xxxz[k] = 2.0 * to_yyz_xxxyz[k] * tke_0;

            to_0_y_yyz_xxyy[k] = -2.0 * to_yyz_xxy[k] + 2.0 * to_yyz_xxyyy[k] * tke_0;

            to_0_y_yyz_xxyz[k] = -to_yyz_xxz[k] + 2.0 * to_yyz_xxyyz[k] * tke_0;

            to_0_y_yyz_xxzz[k] = 2.0 * to_yyz_xxyzz[k] * tke_0;

            to_0_y_yyz_xyyy[k] = -3.0 * to_yyz_xyy[k] + 2.0 * to_yyz_xyyyy[k] * tke_0;

            to_0_y_yyz_xyyz[k] = -2.0 * to_yyz_xyz[k] + 2.0 * to_yyz_xyyyz[k] * tke_0;

            to_0_y_yyz_xyzz[k] = -to_yyz_xzz[k] + 2.0 * to_yyz_xyyzz[k] * tke_0;

            to_0_y_yyz_xzzz[k] = 2.0 * to_yyz_xyzzz[k] * tke_0;

            to_0_y_yyz_yyyy[k] = -4.0 * to_yyz_yyy[k] + 2.0 * to_yyz_yyyyy[k] * tke_0;

            to_0_y_yyz_yyyz[k] = -3.0 * to_yyz_yyz[k] + 2.0 * to_yyz_yyyyz[k] * tke_0;

            to_0_y_yyz_yyzz[k] = -2.0 * to_yyz_yzz[k] + 2.0 * to_yyz_yyyzz[k] * tke_0;

            to_0_y_yyz_yzzz[k] = -to_yyz_zzz[k] + 2.0 * to_yyz_yyzzz[k] * tke_0;

            to_0_y_yyz_zzzz[k] = 2.0 * to_yyz_yzzzz[k] * tke_0;
        }

        // Set up 270-285 components of targeted buffer : FG

        auto to_0_y_yzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 120);

        auto to_0_y_yzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 121);

        auto to_0_y_yzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 122);

        auto to_0_y_yzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 123);

        auto to_0_y_yzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 124);

        auto to_0_y_yzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 125);

        auto to_0_y_yzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 126);

        auto to_0_y_yzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 127);

        auto to_0_y_yzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 128);

        auto to_0_y_yzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 129);

        auto to_0_y_yzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 130);

        auto to_0_y_yzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 131);

        auto to_0_y_yzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 132);

        auto to_0_y_yzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 133);

        auto to_0_y_yzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_0_y_yzz_xxxx, to_0_y_yzz_xxxy, to_0_y_yzz_xxxz, to_0_y_yzz_xxyy, to_0_y_yzz_xxyz, to_0_y_yzz_xxzz, to_0_y_yzz_xyyy, to_0_y_yzz_xyyz, to_0_y_yzz_xyzz, to_0_y_yzz_xzzz, to_0_y_yzz_yyyy, to_0_y_yzz_yyyz, to_0_y_yzz_yyzz, to_0_y_yzz_yzzz, to_0_y_yzz_zzzz, to_yzz_xxx, to_yzz_xxxxy, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_yyy, to_yzz_yyyyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzz_xxxx[k] = 2.0 * to_yzz_xxxxy[k] * tke_0;

            to_0_y_yzz_xxxy[k] = -to_yzz_xxx[k] + 2.0 * to_yzz_xxxyy[k] * tke_0;

            to_0_y_yzz_xxxz[k] = 2.0 * to_yzz_xxxyz[k] * tke_0;

            to_0_y_yzz_xxyy[k] = -2.0 * to_yzz_xxy[k] + 2.0 * to_yzz_xxyyy[k] * tke_0;

            to_0_y_yzz_xxyz[k] = -to_yzz_xxz[k] + 2.0 * to_yzz_xxyyz[k] * tke_0;

            to_0_y_yzz_xxzz[k] = 2.0 * to_yzz_xxyzz[k] * tke_0;

            to_0_y_yzz_xyyy[k] = -3.0 * to_yzz_xyy[k] + 2.0 * to_yzz_xyyyy[k] * tke_0;

            to_0_y_yzz_xyyz[k] = -2.0 * to_yzz_xyz[k] + 2.0 * to_yzz_xyyyz[k] * tke_0;

            to_0_y_yzz_xyzz[k] = -to_yzz_xzz[k] + 2.0 * to_yzz_xyyzz[k] * tke_0;

            to_0_y_yzz_xzzz[k] = 2.0 * to_yzz_xyzzz[k] * tke_0;

            to_0_y_yzz_yyyy[k] = -4.0 * to_yzz_yyy[k] + 2.0 * to_yzz_yyyyy[k] * tke_0;

            to_0_y_yzz_yyyz[k] = -3.0 * to_yzz_yyz[k] + 2.0 * to_yzz_yyyyz[k] * tke_0;

            to_0_y_yzz_yyzz[k] = -2.0 * to_yzz_yzz[k] + 2.0 * to_yzz_yyyzz[k] * tke_0;

            to_0_y_yzz_yzzz[k] = -to_yzz_zzz[k] + 2.0 * to_yzz_yyzzz[k] * tke_0;

            to_0_y_yzz_zzzz[k] = 2.0 * to_yzz_yzzzz[k] * tke_0;
        }

        // Set up 285-300 components of targeted buffer : FG

        auto to_0_y_zzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 135);

        auto to_0_y_zzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 136);

        auto to_0_y_zzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 137);

        auto to_0_y_zzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 138);

        auto to_0_y_zzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 139);

        auto to_0_y_zzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 140);

        auto to_0_y_zzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 141);

        auto to_0_y_zzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 142);

        auto to_0_y_zzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 143);

        auto to_0_y_zzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 144);

        auto to_0_y_zzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 145);

        auto to_0_y_zzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 146);

        auto to_0_y_zzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 147);

        auto to_0_y_zzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 148);

        auto to_0_y_zzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 1 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_y_zzz_xxxx, to_0_y_zzz_xxxy, to_0_y_zzz_xxxz, to_0_y_zzz_xxyy, to_0_y_zzz_xxyz, to_0_y_zzz_xxzz, to_0_y_zzz_xyyy, to_0_y_zzz_xyyz, to_0_y_zzz_xyzz, to_0_y_zzz_xzzz, to_0_y_zzz_yyyy, to_0_y_zzz_yyyz, to_0_y_zzz_yyzz, to_0_y_zzz_yzzz, to_0_y_zzz_zzzz, to_zzz_xxx, to_zzz_xxxxy, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_yyy, to_zzz_yyyyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzz_xxxx[k] = 2.0 * to_zzz_xxxxy[k] * tke_0;

            to_0_y_zzz_xxxy[k] = -to_zzz_xxx[k] + 2.0 * to_zzz_xxxyy[k] * tke_0;

            to_0_y_zzz_xxxz[k] = 2.0 * to_zzz_xxxyz[k] * tke_0;

            to_0_y_zzz_xxyy[k] = -2.0 * to_zzz_xxy[k] + 2.0 * to_zzz_xxyyy[k] * tke_0;

            to_0_y_zzz_xxyz[k] = -to_zzz_xxz[k] + 2.0 * to_zzz_xxyyz[k] * tke_0;

            to_0_y_zzz_xxzz[k] = 2.0 * to_zzz_xxyzz[k] * tke_0;

            to_0_y_zzz_xyyy[k] = -3.0 * to_zzz_xyy[k] + 2.0 * to_zzz_xyyyy[k] * tke_0;

            to_0_y_zzz_xyyz[k] = -2.0 * to_zzz_xyz[k] + 2.0 * to_zzz_xyyyz[k] * tke_0;

            to_0_y_zzz_xyzz[k] = -to_zzz_xzz[k] + 2.0 * to_zzz_xyyzz[k] * tke_0;

            to_0_y_zzz_xzzz[k] = 2.0 * to_zzz_xyzzz[k] * tke_0;

            to_0_y_zzz_yyyy[k] = -4.0 * to_zzz_yyy[k] + 2.0 * to_zzz_yyyyy[k] * tke_0;

            to_0_y_zzz_yyyz[k] = -3.0 * to_zzz_yyz[k] + 2.0 * to_zzz_yyyyz[k] * tke_0;

            to_0_y_zzz_yyzz[k] = -2.0 * to_zzz_yzz[k] + 2.0 * to_zzz_yyyzz[k] * tke_0;

            to_0_y_zzz_yzzz[k] = -to_zzz_zzz[k] + 2.0 * to_zzz_yyzzz[k] * tke_0;

            to_0_y_zzz_zzzz[k] = 2.0 * to_zzz_yzzzz[k] * tke_0;
        }

        // Set up 300-315 components of targeted buffer : FG

        auto to_0_z_xxx_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 0);

        auto to_0_z_xxx_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 1);

        auto to_0_z_xxx_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 2);

        auto to_0_z_xxx_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 3);

        auto to_0_z_xxx_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 4);

        auto to_0_z_xxx_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 5);

        auto to_0_z_xxx_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 6);

        auto to_0_z_xxx_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 7);

        auto to_0_z_xxx_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 8);

        auto to_0_z_xxx_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 9);

        auto to_0_z_xxx_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 10);

        auto to_0_z_xxx_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 11);

        auto to_0_z_xxx_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 12);

        auto to_0_z_xxx_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 13);

        auto to_0_z_xxx_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_0_z_xxx_xxxx, to_0_z_xxx_xxxy, to_0_z_xxx_xxxz, to_0_z_xxx_xxyy, to_0_z_xxx_xxyz, to_0_z_xxx_xxzz, to_0_z_xxx_xyyy, to_0_z_xxx_xyyz, to_0_z_xxx_xyzz, to_0_z_xxx_xzzz, to_0_z_xxx_yyyy, to_0_z_xxx_yyyz, to_0_z_xxx_yyzz, to_0_z_xxx_yzzz, to_0_z_xxx_zzzz, to_xxx_xxx, to_xxx_xxxxz, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxx_xxxx[k] = 2.0 * to_xxx_xxxxz[k] * tke_0;

            to_0_z_xxx_xxxy[k] = 2.0 * to_xxx_xxxyz[k] * tke_0;

            to_0_z_xxx_xxxz[k] = -to_xxx_xxx[k] + 2.0 * to_xxx_xxxzz[k] * tke_0;

            to_0_z_xxx_xxyy[k] = 2.0 * to_xxx_xxyyz[k] * tke_0;

            to_0_z_xxx_xxyz[k] = -to_xxx_xxy[k] + 2.0 * to_xxx_xxyzz[k] * tke_0;

            to_0_z_xxx_xxzz[k] = -2.0 * to_xxx_xxz[k] + 2.0 * to_xxx_xxzzz[k] * tke_0;

            to_0_z_xxx_xyyy[k] = 2.0 * to_xxx_xyyyz[k] * tke_0;

            to_0_z_xxx_xyyz[k] = -to_xxx_xyy[k] + 2.0 * to_xxx_xyyzz[k] * tke_0;

            to_0_z_xxx_xyzz[k] = -2.0 * to_xxx_xyz[k] + 2.0 * to_xxx_xyzzz[k] * tke_0;

            to_0_z_xxx_xzzz[k] = -3.0 * to_xxx_xzz[k] + 2.0 * to_xxx_xzzzz[k] * tke_0;

            to_0_z_xxx_yyyy[k] = 2.0 * to_xxx_yyyyz[k] * tke_0;

            to_0_z_xxx_yyyz[k] = -to_xxx_yyy[k] + 2.0 * to_xxx_yyyzz[k] * tke_0;

            to_0_z_xxx_yyzz[k] = -2.0 * to_xxx_yyz[k] + 2.0 * to_xxx_yyzzz[k] * tke_0;

            to_0_z_xxx_yzzz[k] = -3.0 * to_xxx_yzz[k] + 2.0 * to_xxx_yzzzz[k] * tke_0;

            to_0_z_xxx_zzzz[k] = -4.0 * to_xxx_zzz[k] + 2.0 * to_xxx_zzzzz[k] * tke_0;
        }

        // Set up 315-330 components of targeted buffer : FG

        auto to_0_z_xxy_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 15);

        auto to_0_z_xxy_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 16);

        auto to_0_z_xxy_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 17);

        auto to_0_z_xxy_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 18);

        auto to_0_z_xxy_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 19);

        auto to_0_z_xxy_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 20);

        auto to_0_z_xxy_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 21);

        auto to_0_z_xxy_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 22);

        auto to_0_z_xxy_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 23);

        auto to_0_z_xxy_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 24);

        auto to_0_z_xxy_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 25);

        auto to_0_z_xxy_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 26);

        auto to_0_z_xxy_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 27);

        auto to_0_z_xxy_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 28);

        auto to_0_z_xxy_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_z_xxy_xxxx, to_0_z_xxy_xxxy, to_0_z_xxy_xxxz, to_0_z_xxy_xxyy, to_0_z_xxy_xxyz, to_0_z_xxy_xxzz, to_0_z_xxy_xyyy, to_0_z_xxy_xyyz, to_0_z_xxy_xyzz, to_0_z_xxy_xzzz, to_0_z_xxy_yyyy, to_0_z_xxy_yyyz, to_0_z_xxy_yyzz, to_0_z_xxy_yzzz, to_0_z_xxy_zzzz, to_xxy_xxx, to_xxy_xxxxz, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxy_xxxx[k] = 2.0 * to_xxy_xxxxz[k] * tke_0;

            to_0_z_xxy_xxxy[k] = 2.0 * to_xxy_xxxyz[k] * tke_0;

            to_0_z_xxy_xxxz[k] = -to_xxy_xxx[k] + 2.0 * to_xxy_xxxzz[k] * tke_0;

            to_0_z_xxy_xxyy[k] = 2.0 * to_xxy_xxyyz[k] * tke_0;

            to_0_z_xxy_xxyz[k] = -to_xxy_xxy[k] + 2.0 * to_xxy_xxyzz[k] * tke_0;

            to_0_z_xxy_xxzz[k] = -2.0 * to_xxy_xxz[k] + 2.0 * to_xxy_xxzzz[k] * tke_0;

            to_0_z_xxy_xyyy[k] = 2.0 * to_xxy_xyyyz[k] * tke_0;

            to_0_z_xxy_xyyz[k] = -to_xxy_xyy[k] + 2.0 * to_xxy_xyyzz[k] * tke_0;

            to_0_z_xxy_xyzz[k] = -2.0 * to_xxy_xyz[k] + 2.0 * to_xxy_xyzzz[k] * tke_0;

            to_0_z_xxy_xzzz[k] = -3.0 * to_xxy_xzz[k] + 2.0 * to_xxy_xzzzz[k] * tke_0;

            to_0_z_xxy_yyyy[k] = 2.0 * to_xxy_yyyyz[k] * tke_0;

            to_0_z_xxy_yyyz[k] = -to_xxy_yyy[k] + 2.0 * to_xxy_yyyzz[k] * tke_0;

            to_0_z_xxy_yyzz[k] = -2.0 * to_xxy_yyz[k] + 2.0 * to_xxy_yyzzz[k] * tke_0;

            to_0_z_xxy_yzzz[k] = -3.0 * to_xxy_yzz[k] + 2.0 * to_xxy_yzzzz[k] * tke_0;

            to_0_z_xxy_zzzz[k] = -4.0 * to_xxy_zzz[k] + 2.0 * to_xxy_zzzzz[k] * tke_0;
        }

        // Set up 330-345 components of targeted buffer : FG

        auto to_0_z_xxz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 30);

        auto to_0_z_xxz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 31);

        auto to_0_z_xxz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 32);

        auto to_0_z_xxz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 33);

        auto to_0_z_xxz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 34);

        auto to_0_z_xxz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 35);

        auto to_0_z_xxz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 36);

        auto to_0_z_xxz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 37);

        auto to_0_z_xxz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 38);

        auto to_0_z_xxz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 39);

        auto to_0_z_xxz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 40);

        auto to_0_z_xxz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 41);

        auto to_0_z_xxz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 42);

        auto to_0_z_xxz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 43);

        auto to_0_z_xxz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_0_z_xxz_xxxx, to_0_z_xxz_xxxy, to_0_z_xxz_xxxz, to_0_z_xxz_xxyy, to_0_z_xxz_xxyz, to_0_z_xxz_xxzz, to_0_z_xxz_xyyy, to_0_z_xxz_xyyz, to_0_z_xxz_xyzz, to_0_z_xxz_xzzz, to_0_z_xxz_yyyy, to_0_z_xxz_yyyz, to_0_z_xxz_yyzz, to_0_z_xxz_yzzz, to_0_z_xxz_zzzz, to_xxz_xxx, to_xxz_xxxxz, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_xxz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxz_xxxx[k] = 2.0 * to_xxz_xxxxz[k] * tke_0;

            to_0_z_xxz_xxxy[k] = 2.0 * to_xxz_xxxyz[k] * tke_0;

            to_0_z_xxz_xxxz[k] = -to_xxz_xxx[k] + 2.0 * to_xxz_xxxzz[k] * tke_0;

            to_0_z_xxz_xxyy[k] = 2.0 * to_xxz_xxyyz[k] * tke_0;

            to_0_z_xxz_xxyz[k] = -to_xxz_xxy[k] + 2.0 * to_xxz_xxyzz[k] * tke_0;

            to_0_z_xxz_xxzz[k] = -2.0 * to_xxz_xxz[k] + 2.0 * to_xxz_xxzzz[k] * tke_0;

            to_0_z_xxz_xyyy[k] = 2.0 * to_xxz_xyyyz[k] * tke_0;

            to_0_z_xxz_xyyz[k] = -to_xxz_xyy[k] + 2.0 * to_xxz_xyyzz[k] * tke_0;

            to_0_z_xxz_xyzz[k] = -2.0 * to_xxz_xyz[k] + 2.0 * to_xxz_xyzzz[k] * tke_0;

            to_0_z_xxz_xzzz[k] = -3.0 * to_xxz_xzz[k] + 2.0 * to_xxz_xzzzz[k] * tke_0;

            to_0_z_xxz_yyyy[k] = 2.0 * to_xxz_yyyyz[k] * tke_0;

            to_0_z_xxz_yyyz[k] = -to_xxz_yyy[k] + 2.0 * to_xxz_yyyzz[k] * tke_0;

            to_0_z_xxz_yyzz[k] = -2.0 * to_xxz_yyz[k] + 2.0 * to_xxz_yyzzz[k] * tke_0;

            to_0_z_xxz_yzzz[k] = -3.0 * to_xxz_yzz[k] + 2.0 * to_xxz_yzzzz[k] * tke_0;

            to_0_z_xxz_zzzz[k] = -4.0 * to_xxz_zzz[k] + 2.0 * to_xxz_zzzzz[k] * tke_0;
        }

        // Set up 345-360 components of targeted buffer : FG

        auto to_0_z_xyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 45);

        auto to_0_z_xyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 46);

        auto to_0_z_xyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 47);

        auto to_0_z_xyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 48);

        auto to_0_z_xyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 49);

        auto to_0_z_xyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 50);

        auto to_0_z_xyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 51);

        auto to_0_z_xyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 52);

        auto to_0_z_xyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 53);

        auto to_0_z_xyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 54);

        auto to_0_z_xyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 55);

        auto to_0_z_xyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 56);

        auto to_0_z_xyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 57);

        auto to_0_z_xyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 58);

        auto to_0_z_xyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_z_xyy_xxxx, to_0_z_xyy_xxxy, to_0_z_xyy_xxxz, to_0_z_xyy_xxyy, to_0_z_xyy_xxyz, to_0_z_xyy_xxzz, to_0_z_xyy_xyyy, to_0_z_xyy_xyyz, to_0_z_xyy_xyzz, to_0_z_xyy_xzzz, to_0_z_xyy_yyyy, to_0_z_xyy_yyyz, to_0_z_xyy_yyzz, to_0_z_xyy_yzzz, to_0_z_xyy_zzzz, to_xyy_xxx, to_xyy_xxxxz, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyy_xxxx[k] = 2.0 * to_xyy_xxxxz[k] * tke_0;

            to_0_z_xyy_xxxy[k] = 2.0 * to_xyy_xxxyz[k] * tke_0;

            to_0_z_xyy_xxxz[k] = -to_xyy_xxx[k] + 2.0 * to_xyy_xxxzz[k] * tke_0;

            to_0_z_xyy_xxyy[k] = 2.0 * to_xyy_xxyyz[k] * tke_0;

            to_0_z_xyy_xxyz[k] = -to_xyy_xxy[k] + 2.0 * to_xyy_xxyzz[k] * tke_0;

            to_0_z_xyy_xxzz[k] = -2.0 * to_xyy_xxz[k] + 2.0 * to_xyy_xxzzz[k] * tke_0;

            to_0_z_xyy_xyyy[k] = 2.0 * to_xyy_xyyyz[k] * tke_0;

            to_0_z_xyy_xyyz[k] = -to_xyy_xyy[k] + 2.0 * to_xyy_xyyzz[k] * tke_0;

            to_0_z_xyy_xyzz[k] = -2.0 * to_xyy_xyz[k] + 2.0 * to_xyy_xyzzz[k] * tke_0;

            to_0_z_xyy_xzzz[k] = -3.0 * to_xyy_xzz[k] + 2.0 * to_xyy_xzzzz[k] * tke_0;

            to_0_z_xyy_yyyy[k] = 2.0 * to_xyy_yyyyz[k] * tke_0;

            to_0_z_xyy_yyyz[k] = -to_xyy_yyy[k] + 2.0 * to_xyy_yyyzz[k] * tke_0;

            to_0_z_xyy_yyzz[k] = -2.0 * to_xyy_yyz[k] + 2.0 * to_xyy_yyzzz[k] * tke_0;

            to_0_z_xyy_yzzz[k] = -3.0 * to_xyy_yzz[k] + 2.0 * to_xyy_yzzzz[k] * tke_0;

            to_0_z_xyy_zzzz[k] = -4.0 * to_xyy_zzz[k] + 2.0 * to_xyy_zzzzz[k] * tke_0;
        }

        // Set up 360-375 components of targeted buffer : FG

        auto to_0_z_xyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 60);

        auto to_0_z_xyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 61);

        auto to_0_z_xyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 62);

        auto to_0_z_xyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 63);

        auto to_0_z_xyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 64);

        auto to_0_z_xyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 65);

        auto to_0_z_xyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 66);

        auto to_0_z_xyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 67);

        auto to_0_z_xyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 68);

        auto to_0_z_xyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 69);

        auto to_0_z_xyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 70);

        auto to_0_z_xyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 71);

        auto to_0_z_xyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 72);

        auto to_0_z_xyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 73);

        auto to_0_z_xyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_0_z_xyz_xxxx, to_0_z_xyz_xxxy, to_0_z_xyz_xxxz, to_0_z_xyz_xxyy, to_0_z_xyz_xxyz, to_0_z_xyz_xxzz, to_0_z_xyz_xyyy, to_0_z_xyz_xyyz, to_0_z_xyz_xyzz, to_0_z_xyz_xzzz, to_0_z_xyz_yyyy, to_0_z_xyz_yyyz, to_0_z_xyz_yyzz, to_0_z_xyz_yzzz, to_0_z_xyz_zzzz, to_xyz_xxx, to_xyz_xxxxz, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_xyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyz_xxxx[k] = 2.0 * to_xyz_xxxxz[k] * tke_0;

            to_0_z_xyz_xxxy[k] = 2.0 * to_xyz_xxxyz[k] * tke_0;

            to_0_z_xyz_xxxz[k] = -to_xyz_xxx[k] + 2.0 * to_xyz_xxxzz[k] * tke_0;

            to_0_z_xyz_xxyy[k] = 2.0 * to_xyz_xxyyz[k] * tke_0;

            to_0_z_xyz_xxyz[k] = -to_xyz_xxy[k] + 2.0 * to_xyz_xxyzz[k] * tke_0;

            to_0_z_xyz_xxzz[k] = -2.0 * to_xyz_xxz[k] + 2.0 * to_xyz_xxzzz[k] * tke_0;

            to_0_z_xyz_xyyy[k] = 2.0 * to_xyz_xyyyz[k] * tke_0;

            to_0_z_xyz_xyyz[k] = -to_xyz_xyy[k] + 2.0 * to_xyz_xyyzz[k] * tke_0;

            to_0_z_xyz_xyzz[k] = -2.0 * to_xyz_xyz[k] + 2.0 * to_xyz_xyzzz[k] * tke_0;

            to_0_z_xyz_xzzz[k] = -3.0 * to_xyz_xzz[k] + 2.0 * to_xyz_xzzzz[k] * tke_0;

            to_0_z_xyz_yyyy[k] = 2.0 * to_xyz_yyyyz[k] * tke_0;

            to_0_z_xyz_yyyz[k] = -to_xyz_yyy[k] + 2.0 * to_xyz_yyyzz[k] * tke_0;

            to_0_z_xyz_yyzz[k] = -2.0 * to_xyz_yyz[k] + 2.0 * to_xyz_yyzzz[k] * tke_0;

            to_0_z_xyz_yzzz[k] = -3.0 * to_xyz_yzz[k] + 2.0 * to_xyz_yzzzz[k] * tke_0;

            to_0_z_xyz_zzzz[k] = -4.0 * to_xyz_zzz[k] + 2.0 * to_xyz_zzzzz[k] * tke_0;
        }

        // Set up 375-390 components of targeted buffer : FG

        auto to_0_z_xzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 75);

        auto to_0_z_xzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 76);

        auto to_0_z_xzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 77);

        auto to_0_z_xzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 78);

        auto to_0_z_xzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 79);

        auto to_0_z_xzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 80);

        auto to_0_z_xzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 81);

        auto to_0_z_xzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 82);

        auto to_0_z_xzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 83);

        auto to_0_z_xzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 84);

        auto to_0_z_xzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 85);

        auto to_0_z_xzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 86);

        auto to_0_z_xzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 87);

        auto to_0_z_xzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 88);

        auto to_0_z_xzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_z_xzz_xxxx, to_0_z_xzz_xxxy, to_0_z_xzz_xxxz, to_0_z_xzz_xxyy, to_0_z_xzz_xxyz, to_0_z_xzz_xxzz, to_0_z_xzz_xyyy, to_0_z_xzz_xyyz, to_0_z_xzz_xyzz, to_0_z_xzz_xzzz, to_0_z_xzz_yyyy, to_0_z_xzz_yyyz, to_0_z_xzz_yyzz, to_0_z_xzz_yzzz, to_0_z_xzz_zzzz, to_xzz_xxx, to_xzz_xxxxz, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_xzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzz_xxxx[k] = 2.0 * to_xzz_xxxxz[k] * tke_0;

            to_0_z_xzz_xxxy[k] = 2.0 * to_xzz_xxxyz[k] * tke_0;

            to_0_z_xzz_xxxz[k] = -to_xzz_xxx[k] + 2.0 * to_xzz_xxxzz[k] * tke_0;

            to_0_z_xzz_xxyy[k] = 2.0 * to_xzz_xxyyz[k] * tke_0;

            to_0_z_xzz_xxyz[k] = -to_xzz_xxy[k] + 2.0 * to_xzz_xxyzz[k] * tke_0;

            to_0_z_xzz_xxzz[k] = -2.0 * to_xzz_xxz[k] + 2.0 * to_xzz_xxzzz[k] * tke_0;

            to_0_z_xzz_xyyy[k] = 2.0 * to_xzz_xyyyz[k] * tke_0;

            to_0_z_xzz_xyyz[k] = -to_xzz_xyy[k] + 2.0 * to_xzz_xyyzz[k] * tke_0;

            to_0_z_xzz_xyzz[k] = -2.0 * to_xzz_xyz[k] + 2.0 * to_xzz_xyzzz[k] * tke_0;

            to_0_z_xzz_xzzz[k] = -3.0 * to_xzz_xzz[k] + 2.0 * to_xzz_xzzzz[k] * tke_0;

            to_0_z_xzz_yyyy[k] = 2.0 * to_xzz_yyyyz[k] * tke_0;

            to_0_z_xzz_yyyz[k] = -to_xzz_yyy[k] + 2.0 * to_xzz_yyyzz[k] * tke_0;

            to_0_z_xzz_yyzz[k] = -2.0 * to_xzz_yyz[k] + 2.0 * to_xzz_yyzzz[k] * tke_0;

            to_0_z_xzz_yzzz[k] = -3.0 * to_xzz_yzz[k] + 2.0 * to_xzz_yzzzz[k] * tke_0;

            to_0_z_xzz_zzzz[k] = -4.0 * to_xzz_zzz[k] + 2.0 * to_xzz_zzzzz[k] * tke_0;
        }

        // Set up 390-405 components of targeted buffer : FG

        auto to_0_z_yyy_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 90);

        auto to_0_z_yyy_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 91);

        auto to_0_z_yyy_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 92);

        auto to_0_z_yyy_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 93);

        auto to_0_z_yyy_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 94);

        auto to_0_z_yyy_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 95);

        auto to_0_z_yyy_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 96);

        auto to_0_z_yyy_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 97);

        auto to_0_z_yyy_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 98);

        auto to_0_z_yyy_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 99);

        auto to_0_z_yyy_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 100);

        auto to_0_z_yyy_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 101);

        auto to_0_z_yyy_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 102);

        auto to_0_z_yyy_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 103);

        auto to_0_z_yyy_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_0_z_yyy_xxxx, to_0_z_yyy_xxxy, to_0_z_yyy_xxxz, to_0_z_yyy_xxyy, to_0_z_yyy_xxyz, to_0_z_yyy_xxzz, to_0_z_yyy_xyyy, to_0_z_yyy_xyyz, to_0_z_yyy_xyzz, to_0_z_yyy_xzzz, to_0_z_yyy_yyyy, to_0_z_yyy_yyyz, to_0_z_yyy_yyzz, to_0_z_yyy_yzzz, to_0_z_yyy_zzzz, to_yyy_xxx, to_yyy_xxxxz, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyy_xxxx[k] = 2.0 * to_yyy_xxxxz[k] * tke_0;

            to_0_z_yyy_xxxy[k] = 2.0 * to_yyy_xxxyz[k] * tke_0;

            to_0_z_yyy_xxxz[k] = -to_yyy_xxx[k] + 2.0 * to_yyy_xxxzz[k] * tke_0;

            to_0_z_yyy_xxyy[k] = 2.0 * to_yyy_xxyyz[k] * tke_0;

            to_0_z_yyy_xxyz[k] = -to_yyy_xxy[k] + 2.0 * to_yyy_xxyzz[k] * tke_0;

            to_0_z_yyy_xxzz[k] = -2.0 * to_yyy_xxz[k] + 2.0 * to_yyy_xxzzz[k] * tke_0;

            to_0_z_yyy_xyyy[k] = 2.0 * to_yyy_xyyyz[k] * tke_0;

            to_0_z_yyy_xyyz[k] = -to_yyy_xyy[k] + 2.0 * to_yyy_xyyzz[k] * tke_0;

            to_0_z_yyy_xyzz[k] = -2.0 * to_yyy_xyz[k] + 2.0 * to_yyy_xyzzz[k] * tke_0;

            to_0_z_yyy_xzzz[k] = -3.0 * to_yyy_xzz[k] + 2.0 * to_yyy_xzzzz[k] * tke_0;

            to_0_z_yyy_yyyy[k] = 2.0 * to_yyy_yyyyz[k] * tke_0;

            to_0_z_yyy_yyyz[k] = -to_yyy_yyy[k] + 2.0 * to_yyy_yyyzz[k] * tke_0;

            to_0_z_yyy_yyzz[k] = -2.0 * to_yyy_yyz[k] + 2.0 * to_yyy_yyzzz[k] * tke_0;

            to_0_z_yyy_yzzz[k] = -3.0 * to_yyy_yzz[k] + 2.0 * to_yyy_yzzzz[k] * tke_0;

            to_0_z_yyy_zzzz[k] = -4.0 * to_yyy_zzz[k] + 2.0 * to_yyy_zzzzz[k] * tke_0;
        }

        // Set up 405-420 components of targeted buffer : FG

        auto to_0_z_yyz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 105);

        auto to_0_z_yyz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 106);

        auto to_0_z_yyz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 107);

        auto to_0_z_yyz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 108);

        auto to_0_z_yyz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 109);

        auto to_0_z_yyz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 110);

        auto to_0_z_yyz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 111);

        auto to_0_z_yyz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 112);

        auto to_0_z_yyz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 113);

        auto to_0_z_yyz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 114);

        auto to_0_z_yyz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 115);

        auto to_0_z_yyz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 116);

        auto to_0_z_yyz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 117);

        auto to_0_z_yyz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 118);

        auto to_0_z_yyz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_z_yyz_xxxx, to_0_z_yyz_xxxy, to_0_z_yyz_xxxz, to_0_z_yyz_xxyy, to_0_z_yyz_xxyz, to_0_z_yyz_xxzz, to_0_z_yyz_xyyy, to_0_z_yyz_xyyz, to_0_z_yyz_xyzz, to_0_z_yyz_xzzz, to_0_z_yyz_yyyy, to_0_z_yyz_yyyz, to_0_z_yyz_yyzz, to_0_z_yyz_yzzz, to_0_z_yyz_zzzz, to_yyz_xxx, to_yyz_xxxxz, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, to_yyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyz_xxxx[k] = 2.0 * to_yyz_xxxxz[k] * tke_0;

            to_0_z_yyz_xxxy[k] = 2.0 * to_yyz_xxxyz[k] * tke_0;

            to_0_z_yyz_xxxz[k] = -to_yyz_xxx[k] + 2.0 * to_yyz_xxxzz[k] * tke_0;

            to_0_z_yyz_xxyy[k] = 2.0 * to_yyz_xxyyz[k] * tke_0;

            to_0_z_yyz_xxyz[k] = -to_yyz_xxy[k] + 2.0 * to_yyz_xxyzz[k] * tke_0;

            to_0_z_yyz_xxzz[k] = -2.0 * to_yyz_xxz[k] + 2.0 * to_yyz_xxzzz[k] * tke_0;

            to_0_z_yyz_xyyy[k] = 2.0 * to_yyz_xyyyz[k] * tke_0;

            to_0_z_yyz_xyyz[k] = -to_yyz_xyy[k] + 2.0 * to_yyz_xyyzz[k] * tke_0;

            to_0_z_yyz_xyzz[k] = -2.0 * to_yyz_xyz[k] + 2.0 * to_yyz_xyzzz[k] * tke_0;

            to_0_z_yyz_xzzz[k] = -3.0 * to_yyz_xzz[k] + 2.0 * to_yyz_xzzzz[k] * tke_0;

            to_0_z_yyz_yyyy[k] = 2.0 * to_yyz_yyyyz[k] * tke_0;

            to_0_z_yyz_yyyz[k] = -to_yyz_yyy[k] + 2.0 * to_yyz_yyyzz[k] * tke_0;

            to_0_z_yyz_yyzz[k] = -2.0 * to_yyz_yyz[k] + 2.0 * to_yyz_yyzzz[k] * tke_0;

            to_0_z_yyz_yzzz[k] = -3.0 * to_yyz_yzz[k] + 2.0 * to_yyz_yzzzz[k] * tke_0;

            to_0_z_yyz_zzzz[k] = -4.0 * to_yyz_zzz[k] + 2.0 * to_yyz_zzzzz[k] * tke_0;
        }

        // Set up 420-435 components of targeted buffer : FG

        auto to_0_z_yzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 120);

        auto to_0_z_yzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 121);

        auto to_0_z_yzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 122);

        auto to_0_z_yzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 123);

        auto to_0_z_yzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 124);

        auto to_0_z_yzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 125);

        auto to_0_z_yzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 126);

        auto to_0_z_yzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 127);

        auto to_0_z_yzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 128);

        auto to_0_z_yzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 129);

        auto to_0_z_yzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 130);

        auto to_0_z_yzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 131);

        auto to_0_z_yzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 132);

        auto to_0_z_yzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 133);

        auto to_0_z_yzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_0_z_yzz_xxxx, to_0_z_yzz_xxxy, to_0_z_yzz_xxxz, to_0_z_yzz_xxyy, to_0_z_yzz_xxyz, to_0_z_yzz_xxzz, to_0_z_yzz_xyyy, to_0_z_yzz_xyyz, to_0_z_yzz_xyzz, to_0_z_yzz_xzzz, to_0_z_yzz_yyyy, to_0_z_yzz_yyyz, to_0_z_yzz_yyzz, to_0_z_yzz_yzzz, to_0_z_yzz_zzzz, to_yzz_xxx, to_yzz_xxxxz, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, to_yzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzz_xxxx[k] = 2.0 * to_yzz_xxxxz[k] * tke_0;

            to_0_z_yzz_xxxy[k] = 2.0 * to_yzz_xxxyz[k] * tke_0;

            to_0_z_yzz_xxxz[k] = -to_yzz_xxx[k] + 2.0 * to_yzz_xxxzz[k] * tke_0;

            to_0_z_yzz_xxyy[k] = 2.0 * to_yzz_xxyyz[k] * tke_0;

            to_0_z_yzz_xxyz[k] = -to_yzz_xxy[k] + 2.0 * to_yzz_xxyzz[k] * tke_0;

            to_0_z_yzz_xxzz[k] = -2.0 * to_yzz_xxz[k] + 2.0 * to_yzz_xxzzz[k] * tke_0;

            to_0_z_yzz_xyyy[k] = 2.0 * to_yzz_xyyyz[k] * tke_0;

            to_0_z_yzz_xyyz[k] = -to_yzz_xyy[k] + 2.0 * to_yzz_xyyzz[k] * tke_0;

            to_0_z_yzz_xyzz[k] = -2.0 * to_yzz_xyz[k] + 2.0 * to_yzz_xyzzz[k] * tke_0;

            to_0_z_yzz_xzzz[k] = -3.0 * to_yzz_xzz[k] + 2.0 * to_yzz_xzzzz[k] * tke_0;

            to_0_z_yzz_yyyy[k] = 2.0 * to_yzz_yyyyz[k] * tke_0;

            to_0_z_yzz_yyyz[k] = -to_yzz_yyy[k] + 2.0 * to_yzz_yyyzz[k] * tke_0;

            to_0_z_yzz_yyzz[k] = -2.0 * to_yzz_yyz[k] + 2.0 * to_yzz_yyzzz[k] * tke_0;

            to_0_z_yzz_yzzz[k] = -3.0 * to_yzz_yzz[k] + 2.0 * to_yzz_yzzzz[k] * tke_0;

            to_0_z_yzz_zzzz[k] = -4.0 * to_yzz_zzz[k] + 2.0 * to_yzz_zzzzz[k] * tke_0;
        }

        // Set up 435-450 components of targeted buffer : FG

        auto to_0_z_zzz_xxxx = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 135);

        auto to_0_z_zzz_xxxy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 136);

        auto to_0_z_zzz_xxxz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 137);

        auto to_0_z_zzz_xxyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 138);

        auto to_0_z_zzz_xxyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 139);

        auto to_0_z_zzz_xxzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 140);

        auto to_0_z_zzz_xyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 141);

        auto to_0_z_zzz_xyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 142);

        auto to_0_z_zzz_xyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 143);

        auto to_0_z_zzz_xzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 144);

        auto to_0_z_zzz_yyyy = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 145);

        auto to_0_z_zzz_yyyz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 146);

        auto to_0_z_zzz_yyzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 147);

        auto to_0_z_zzz_yzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 148);

        auto to_0_z_zzz_zzzz = pbuffer.data(idx_op_geom_001_fg + 2 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_z_zzz_xxxx, to_0_z_zzz_xxxy, to_0_z_zzz_xxxz, to_0_z_zzz_xxyy, to_0_z_zzz_xxyz, to_0_z_zzz_xxzz, to_0_z_zzz_xyyy, to_0_z_zzz_xyyz, to_0_z_zzz_xyzz, to_0_z_zzz_xzzz, to_0_z_zzz_yyyy, to_0_z_zzz_yyyz, to_0_z_zzz_yyzz, to_0_z_zzz_yzzz, to_0_z_zzz_zzzz, to_zzz_xxx, to_zzz_xxxxz, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, to_zzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzz_xxxx[k] = 2.0 * to_zzz_xxxxz[k] * tke_0;

            to_0_z_zzz_xxxy[k] = 2.0 * to_zzz_xxxyz[k] * tke_0;

            to_0_z_zzz_xxxz[k] = -to_zzz_xxx[k] + 2.0 * to_zzz_xxxzz[k] * tke_0;

            to_0_z_zzz_xxyy[k] = 2.0 * to_zzz_xxyyz[k] * tke_0;

            to_0_z_zzz_xxyz[k] = -to_zzz_xxy[k] + 2.0 * to_zzz_xxyzz[k] * tke_0;

            to_0_z_zzz_xxzz[k] = -2.0 * to_zzz_xxz[k] + 2.0 * to_zzz_xxzzz[k] * tke_0;

            to_0_z_zzz_xyyy[k] = 2.0 * to_zzz_xyyyz[k] * tke_0;

            to_0_z_zzz_xyyz[k] = -to_zzz_xyy[k] + 2.0 * to_zzz_xyyzz[k] * tke_0;

            to_0_z_zzz_xyzz[k] = -2.0 * to_zzz_xyz[k] + 2.0 * to_zzz_xyzzz[k] * tke_0;

            to_0_z_zzz_xzzz[k] = -3.0 * to_zzz_xzz[k] + 2.0 * to_zzz_xzzzz[k] * tke_0;

            to_0_z_zzz_yyyy[k] = 2.0 * to_zzz_yyyyz[k] * tke_0;

            to_0_z_zzz_yyyz[k] = -to_zzz_yyy[k] + 2.0 * to_zzz_yyyzz[k] * tke_0;

            to_0_z_zzz_yyzz[k] = -2.0 * to_zzz_yyz[k] + 2.0 * to_zzz_yyzzz[k] * tke_0;

            to_0_z_zzz_yzzz[k] = -3.0 * to_zzz_yzz[k] + 2.0 * to_zzz_yzzzz[k] * tke_0;

            to_0_z_zzz_zzzz[k] = -4.0 * to_zzz_zzz[k] + 2.0 * to_zzz_zzzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

