#include "GeometricalDerivatives1X1ForGG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_gg(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_gg,
                        const size_t idx_op_ff,
                        const size_t idx_op_fh,
                        const size_t idx_op_hf,
                        const size_t idx_op_hh,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
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

        // Set up components of auxiliary buffer : HF

        auto to_xxxxx_xxx = pbuffer.data(idx_op_hf + i * 210 + 0);

        auto to_xxxxx_xxy = pbuffer.data(idx_op_hf + i * 210 + 1);

        auto to_xxxxx_xxz = pbuffer.data(idx_op_hf + i * 210 + 2);

        auto to_xxxxx_xyy = pbuffer.data(idx_op_hf + i * 210 + 3);

        auto to_xxxxx_xyz = pbuffer.data(idx_op_hf + i * 210 + 4);

        auto to_xxxxx_xzz = pbuffer.data(idx_op_hf + i * 210 + 5);

        auto to_xxxxx_yyy = pbuffer.data(idx_op_hf + i * 210 + 6);

        auto to_xxxxx_yyz = pbuffer.data(idx_op_hf + i * 210 + 7);

        auto to_xxxxx_yzz = pbuffer.data(idx_op_hf + i * 210 + 8);

        auto to_xxxxx_zzz = pbuffer.data(idx_op_hf + i * 210 + 9);

        auto to_xxxxy_xxx = pbuffer.data(idx_op_hf + i * 210 + 10);

        auto to_xxxxy_xxy = pbuffer.data(idx_op_hf + i * 210 + 11);

        auto to_xxxxy_xxz = pbuffer.data(idx_op_hf + i * 210 + 12);

        auto to_xxxxy_xyy = pbuffer.data(idx_op_hf + i * 210 + 13);

        auto to_xxxxy_xyz = pbuffer.data(idx_op_hf + i * 210 + 14);

        auto to_xxxxy_xzz = pbuffer.data(idx_op_hf + i * 210 + 15);

        auto to_xxxxy_yyy = pbuffer.data(idx_op_hf + i * 210 + 16);

        auto to_xxxxy_yyz = pbuffer.data(idx_op_hf + i * 210 + 17);

        auto to_xxxxy_yzz = pbuffer.data(idx_op_hf + i * 210 + 18);

        auto to_xxxxy_zzz = pbuffer.data(idx_op_hf + i * 210 + 19);

        auto to_xxxxz_xxx = pbuffer.data(idx_op_hf + i * 210 + 20);

        auto to_xxxxz_xxy = pbuffer.data(idx_op_hf + i * 210 + 21);

        auto to_xxxxz_xxz = pbuffer.data(idx_op_hf + i * 210 + 22);

        auto to_xxxxz_xyy = pbuffer.data(idx_op_hf + i * 210 + 23);

        auto to_xxxxz_xyz = pbuffer.data(idx_op_hf + i * 210 + 24);

        auto to_xxxxz_xzz = pbuffer.data(idx_op_hf + i * 210 + 25);

        auto to_xxxxz_yyy = pbuffer.data(idx_op_hf + i * 210 + 26);

        auto to_xxxxz_yyz = pbuffer.data(idx_op_hf + i * 210 + 27);

        auto to_xxxxz_yzz = pbuffer.data(idx_op_hf + i * 210 + 28);

        auto to_xxxxz_zzz = pbuffer.data(idx_op_hf + i * 210 + 29);

        auto to_xxxyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 30);

        auto to_xxxyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 31);

        auto to_xxxyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 32);

        auto to_xxxyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 33);

        auto to_xxxyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 34);

        auto to_xxxyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 35);

        auto to_xxxyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 36);

        auto to_xxxyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 37);

        auto to_xxxyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 38);

        auto to_xxxyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 39);

        auto to_xxxyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 40);

        auto to_xxxyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 41);

        auto to_xxxyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 42);

        auto to_xxxyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 43);

        auto to_xxxyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 44);

        auto to_xxxyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 45);

        auto to_xxxyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 46);

        auto to_xxxyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 47);

        auto to_xxxyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 48);

        auto to_xxxyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 49);

        auto to_xxxzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 50);

        auto to_xxxzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 51);

        auto to_xxxzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 52);

        auto to_xxxzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 53);

        auto to_xxxzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 54);

        auto to_xxxzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 55);

        auto to_xxxzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 56);

        auto to_xxxzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 57);

        auto to_xxxzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 58);

        auto to_xxxzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 59);

        auto to_xxyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 60);

        auto to_xxyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 61);

        auto to_xxyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 62);

        auto to_xxyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 63);

        auto to_xxyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 64);

        auto to_xxyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 65);

        auto to_xxyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 66);

        auto to_xxyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 67);

        auto to_xxyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 68);

        auto to_xxyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 69);

        auto to_xxyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 70);

        auto to_xxyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 71);

        auto to_xxyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 72);

        auto to_xxyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 73);

        auto to_xxyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 74);

        auto to_xxyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 75);

        auto to_xxyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 76);

        auto to_xxyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 77);

        auto to_xxyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 78);

        auto to_xxyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 79);

        auto to_xxyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 80);

        auto to_xxyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 81);

        auto to_xxyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 82);

        auto to_xxyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 83);

        auto to_xxyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 84);

        auto to_xxyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 85);

        auto to_xxyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 86);

        auto to_xxyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 87);

        auto to_xxyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 88);

        auto to_xxyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 89);

        auto to_xxzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 90);

        auto to_xxzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 91);

        auto to_xxzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 92);

        auto to_xxzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 93);

        auto to_xxzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 94);

        auto to_xxzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 95);

        auto to_xxzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 96);

        auto to_xxzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 97);

        auto to_xxzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 98);

        auto to_xxzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 99);

        auto to_xyyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 100);

        auto to_xyyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 101);

        auto to_xyyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 102);

        auto to_xyyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 103);

        auto to_xyyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 104);

        auto to_xyyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 105);

        auto to_xyyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 106);

        auto to_xyyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 107);

        auto to_xyyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 108);

        auto to_xyyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 109);

        auto to_xyyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 110);

        auto to_xyyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 111);

        auto to_xyyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 112);

        auto to_xyyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 113);

        auto to_xyyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 114);

        auto to_xyyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 115);

        auto to_xyyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 116);

        auto to_xyyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 117);

        auto to_xyyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 118);

        auto to_xyyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 119);

        auto to_xyyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 120);

        auto to_xyyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 121);

        auto to_xyyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 122);

        auto to_xyyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 123);

        auto to_xyyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 124);

        auto to_xyyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 125);

        auto to_xyyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 126);

        auto to_xyyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 127);

        auto to_xyyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 128);

        auto to_xyyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 129);

        auto to_xyzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 130);

        auto to_xyzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 131);

        auto to_xyzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 132);

        auto to_xyzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 133);

        auto to_xyzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 134);

        auto to_xyzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 135);

        auto to_xyzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 136);

        auto to_xyzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 137);

        auto to_xyzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 138);

        auto to_xyzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 139);

        auto to_xzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 140);

        auto to_xzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 141);

        auto to_xzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 142);

        auto to_xzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 143);

        auto to_xzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 144);

        auto to_xzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 145);

        auto to_xzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 146);

        auto to_xzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 147);

        auto to_xzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 148);

        auto to_xzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 149);

        auto to_yyyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 150);

        auto to_yyyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 151);

        auto to_yyyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 152);

        auto to_yyyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 153);

        auto to_yyyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 154);

        auto to_yyyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 155);

        auto to_yyyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 156);

        auto to_yyyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 157);

        auto to_yyyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 158);

        auto to_yyyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 159);

        auto to_yyyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 160);

        auto to_yyyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 161);

        auto to_yyyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 162);

        auto to_yyyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 163);

        auto to_yyyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 164);

        auto to_yyyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 165);

        auto to_yyyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 166);

        auto to_yyyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 167);

        auto to_yyyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 168);

        auto to_yyyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 169);

        auto to_yyyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 170);

        auto to_yyyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 171);

        auto to_yyyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 172);

        auto to_yyyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 173);

        auto to_yyyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 174);

        auto to_yyyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 175);

        auto to_yyyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 176);

        auto to_yyyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 177);

        auto to_yyyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 178);

        auto to_yyyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 179);

        auto to_yyzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 180);

        auto to_yyzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 181);

        auto to_yyzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 182);

        auto to_yyzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 183);

        auto to_yyzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 184);

        auto to_yyzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 185);

        auto to_yyzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 186);

        auto to_yyzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 187);

        auto to_yyzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 188);

        auto to_yyzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 189);

        auto to_yzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 190);

        auto to_yzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 191);

        auto to_yzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 192);

        auto to_yzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 193);

        auto to_yzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 194);

        auto to_yzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 195);

        auto to_yzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 196);

        auto to_yzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 197);

        auto to_yzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 198);

        auto to_yzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 199);

        auto to_zzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 200);

        auto to_zzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 201);

        auto to_zzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 202);

        auto to_zzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 203);

        auto to_zzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 204);

        auto to_zzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 205);

        auto to_zzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 206);

        auto to_zzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 207);

        auto to_zzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 208);

        auto to_zzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 209);

        // Set up components of auxiliary buffer : HH

        auto to_xxxxx_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 0);

        auto to_xxxxx_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 1);

        auto to_xxxxx_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 2);

        auto to_xxxxx_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 3);

        auto to_xxxxx_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 4);

        auto to_xxxxx_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 5);

        auto to_xxxxx_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 6);

        auto to_xxxxx_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 7);

        auto to_xxxxx_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 8);

        auto to_xxxxx_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 9);

        auto to_xxxxx_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 10);

        auto to_xxxxx_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 11);

        auto to_xxxxx_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 12);

        auto to_xxxxx_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 13);

        auto to_xxxxx_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 14);

        auto to_xxxxx_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 15);

        auto to_xxxxx_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 16);

        auto to_xxxxx_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 17);

        auto to_xxxxx_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 18);

        auto to_xxxxx_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 19);

        auto to_xxxxx_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 20);

        auto to_xxxxy_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 21);

        auto to_xxxxy_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 22);

        auto to_xxxxy_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 23);

        auto to_xxxxy_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 24);

        auto to_xxxxy_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 25);

        auto to_xxxxy_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 26);

        auto to_xxxxy_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 27);

        auto to_xxxxy_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 28);

        auto to_xxxxy_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 29);

        auto to_xxxxy_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 30);

        auto to_xxxxy_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 31);

        auto to_xxxxy_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 32);

        auto to_xxxxy_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 33);

        auto to_xxxxy_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 34);

        auto to_xxxxy_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 35);

        auto to_xxxxy_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 36);

        auto to_xxxxy_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 37);

        auto to_xxxxy_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 38);

        auto to_xxxxy_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 39);

        auto to_xxxxy_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 40);

        auto to_xxxxy_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 41);

        auto to_xxxxz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 42);

        auto to_xxxxz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 43);

        auto to_xxxxz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 44);

        auto to_xxxxz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 45);

        auto to_xxxxz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 46);

        auto to_xxxxz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 47);

        auto to_xxxxz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 48);

        auto to_xxxxz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 49);

        auto to_xxxxz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 50);

        auto to_xxxxz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 51);

        auto to_xxxxz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 52);

        auto to_xxxxz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 53);

        auto to_xxxxz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 54);

        auto to_xxxxz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 55);

        auto to_xxxxz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 56);

        auto to_xxxxz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 57);

        auto to_xxxxz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 58);

        auto to_xxxxz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 59);

        auto to_xxxxz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 60);

        auto to_xxxxz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 61);

        auto to_xxxxz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 62);

        auto to_xxxyy_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 63);

        auto to_xxxyy_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 64);

        auto to_xxxyy_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 65);

        auto to_xxxyy_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 66);

        auto to_xxxyy_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 67);

        auto to_xxxyy_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 68);

        auto to_xxxyy_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 69);

        auto to_xxxyy_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 70);

        auto to_xxxyy_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 71);

        auto to_xxxyy_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 72);

        auto to_xxxyy_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 73);

        auto to_xxxyy_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 74);

        auto to_xxxyy_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 75);

        auto to_xxxyy_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 76);

        auto to_xxxyy_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 77);

        auto to_xxxyy_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 78);

        auto to_xxxyy_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 79);

        auto to_xxxyy_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 80);

        auto to_xxxyy_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 81);

        auto to_xxxyy_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 82);

        auto to_xxxyy_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 83);

        auto to_xxxyz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 84);

        auto to_xxxyz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 85);

        auto to_xxxyz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 86);

        auto to_xxxyz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 87);

        auto to_xxxyz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 88);

        auto to_xxxyz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 89);

        auto to_xxxyz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 90);

        auto to_xxxyz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 91);

        auto to_xxxyz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 92);

        auto to_xxxyz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 93);

        auto to_xxxyz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 94);

        auto to_xxxyz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 95);

        auto to_xxxyz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 96);

        auto to_xxxyz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 97);

        auto to_xxxyz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 98);

        auto to_xxxyz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 99);

        auto to_xxxyz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 100);

        auto to_xxxyz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 101);

        auto to_xxxyz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 102);

        auto to_xxxyz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 103);

        auto to_xxxyz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 104);

        auto to_xxxzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 105);

        auto to_xxxzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 106);

        auto to_xxxzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 107);

        auto to_xxxzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 108);

        auto to_xxxzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 109);

        auto to_xxxzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 110);

        auto to_xxxzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 111);

        auto to_xxxzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 112);

        auto to_xxxzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 113);

        auto to_xxxzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 114);

        auto to_xxxzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 115);

        auto to_xxxzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 116);

        auto to_xxxzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 117);

        auto to_xxxzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 118);

        auto to_xxxzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 119);

        auto to_xxxzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 120);

        auto to_xxxzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 121);

        auto to_xxxzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 122);

        auto to_xxxzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 123);

        auto to_xxxzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 124);

        auto to_xxxzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 125);

        auto to_xxyyy_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 126);

        auto to_xxyyy_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 127);

        auto to_xxyyy_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 128);

        auto to_xxyyy_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 129);

        auto to_xxyyy_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 130);

        auto to_xxyyy_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 131);

        auto to_xxyyy_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 132);

        auto to_xxyyy_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 133);

        auto to_xxyyy_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 134);

        auto to_xxyyy_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 135);

        auto to_xxyyy_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 136);

        auto to_xxyyy_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 137);

        auto to_xxyyy_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 138);

        auto to_xxyyy_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 139);

        auto to_xxyyy_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 140);

        auto to_xxyyy_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 141);

        auto to_xxyyy_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 142);

        auto to_xxyyy_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 143);

        auto to_xxyyy_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 144);

        auto to_xxyyy_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 145);

        auto to_xxyyy_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 146);

        auto to_xxyyz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 147);

        auto to_xxyyz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 148);

        auto to_xxyyz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 149);

        auto to_xxyyz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 150);

        auto to_xxyyz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 151);

        auto to_xxyyz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 152);

        auto to_xxyyz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 153);

        auto to_xxyyz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 154);

        auto to_xxyyz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 155);

        auto to_xxyyz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 156);

        auto to_xxyyz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 157);

        auto to_xxyyz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 158);

        auto to_xxyyz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 159);

        auto to_xxyyz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 160);

        auto to_xxyyz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 161);

        auto to_xxyyz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 162);

        auto to_xxyyz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 163);

        auto to_xxyyz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 164);

        auto to_xxyyz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 165);

        auto to_xxyyz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 166);

        auto to_xxyyz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 167);

        auto to_xxyzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 168);

        auto to_xxyzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 169);

        auto to_xxyzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 170);

        auto to_xxyzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 171);

        auto to_xxyzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 172);

        auto to_xxyzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 173);

        auto to_xxyzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 174);

        auto to_xxyzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 175);

        auto to_xxyzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 176);

        auto to_xxyzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 177);

        auto to_xxyzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 178);

        auto to_xxyzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 179);

        auto to_xxyzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 180);

        auto to_xxyzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 181);

        auto to_xxyzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 182);

        auto to_xxyzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 183);

        auto to_xxyzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 184);

        auto to_xxyzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 185);

        auto to_xxyzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 186);

        auto to_xxyzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 187);

        auto to_xxyzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 188);

        auto to_xxzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 189);

        auto to_xxzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 190);

        auto to_xxzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 191);

        auto to_xxzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 192);

        auto to_xxzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 193);

        auto to_xxzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 194);

        auto to_xxzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 195);

        auto to_xxzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 196);

        auto to_xxzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 197);

        auto to_xxzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 198);

        auto to_xxzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 199);

        auto to_xxzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 200);

        auto to_xxzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 201);

        auto to_xxzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 202);

        auto to_xxzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 203);

        auto to_xxzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 204);

        auto to_xxzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 205);

        auto to_xxzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 206);

        auto to_xxzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 207);

        auto to_xxzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 208);

        auto to_xxzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 209);

        auto to_xyyyy_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 210);

        auto to_xyyyy_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 211);

        auto to_xyyyy_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 212);

        auto to_xyyyy_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 213);

        auto to_xyyyy_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 214);

        auto to_xyyyy_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 215);

        auto to_xyyyy_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 216);

        auto to_xyyyy_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 217);

        auto to_xyyyy_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 218);

        auto to_xyyyy_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 219);

        auto to_xyyyy_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 220);

        auto to_xyyyy_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 221);

        auto to_xyyyy_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 222);

        auto to_xyyyy_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 223);

        auto to_xyyyy_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 224);

        auto to_xyyyy_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 225);

        auto to_xyyyy_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 226);

        auto to_xyyyy_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 227);

        auto to_xyyyy_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 228);

        auto to_xyyyy_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 229);

        auto to_xyyyy_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 230);

        auto to_xyyyz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 231);

        auto to_xyyyz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 232);

        auto to_xyyyz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 233);

        auto to_xyyyz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 234);

        auto to_xyyyz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 235);

        auto to_xyyyz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 236);

        auto to_xyyyz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 237);

        auto to_xyyyz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 238);

        auto to_xyyyz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 239);

        auto to_xyyyz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 240);

        auto to_xyyyz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 241);

        auto to_xyyyz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 242);

        auto to_xyyyz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 243);

        auto to_xyyyz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 244);

        auto to_xyyyz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 245);

        auto to_xyyyz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 246);

        auto to_xyyyz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 247);

        auto to_xyyyz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 248);

        auto to_xyyyz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 249);

        auto to_xyyyz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 250);

        auto to_xyyyz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 251);

        auto to_xyyzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 252);

        auto to_xyyzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 253);

        auto to_xyyzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 254);

        auto to_xyyzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 255);

        auto to_xyyzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 256);

        auto to_xyyzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 257);

        auto to_xyyzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 258);

        auto to_xyyzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 259);

        auto to_xyyzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 260);

        auto to_xyyzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 261);

        auto to_xyyzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 262);

        auto to_xyyzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 263);

        auto to_xyyzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 264);

        auto to_xyyzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 265);

        auto to_xyyzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 266);

        auto to_xyyzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 267);

        auto to_xyyzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 268);

        auto to_xyyzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 269);

        auto to_xyyzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 270);

        auto to_xyyzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 271);

        auto to_xyyzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 272);

        auto to_xyzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 273);

        auto to_xyzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 274);

        auto to_xyzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 275);

        auto to_xyzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 276);

        auto to_xyzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 277);

        auto to_xyzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 278);

        auto to_xyzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 279);

        auto to_xyzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 280);

        auto to_xyzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 281);

        auto to_xyzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 282);

        auto to_xyzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 283);

        auto to_xyzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 284);

        auto to_xyzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 285);

        auto to_xyzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 286);

        auto to_xyzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 287);

        auto to_xyzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 288);

        auto to_xyzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 289);

        auto to_xyzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 290);

        auto to_xyzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 291);

        auto to_xyzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 292);

        auto to_xyzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 293);

        auto to_xzzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 294);

        auto to_xzzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 295);

        auto to_xzzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 296);

        auto to_xzzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 297);

        auto to_xzzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 298);

        auto to_xzzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 299);

        auto to_xzzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 300);

        auto to_xzzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 301);

        auto to_xzzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 302);

        auto to_xzzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 303);

        auto to_xzzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 304);

        auto to_xzzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 305);

        auto to_xzzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 306);

        auto to_xzzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 307);

        auto to_xzzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 308);

        auto to_xzzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 309);

        auto to_xzzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 310);

        auto to_xzzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 311);

        auto to_xzzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 312);

        auto to_xzzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 313);

        auto to_xzzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 314);

        auto to_yyyyy_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 315);

        auto to_yyyyy_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 316);

        auto to_yyyyy_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 317);

        auto to_yyyyy_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 318);

        auto to_yyyyy_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 319);

        auto to_yyyyy_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 320);

        auto to_yyyyy_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 321);

        auto to_yyyyy_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 322);

        auto to_yyyyy_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 323);

        auto to_yyyyy_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 324);

        auto to_yyyyy_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 325);

        auto to_yyyyy_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 326);

        auto to_yyyyy_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 327);

        auto to_yyyyy_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 328);

        auto to_yyyyy_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 329);

        auto to_yyyyy_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 330);

        auto to_yyyyy_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 331);

        auto to_yyyyy_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 332);

        auto to_yyyyy_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 333);

        auto to_yyyyy_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 334);

        auto to_yyyyy_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 335);

        auto to_yyyyz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 336);

        auto to_yyyyz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 337);

        auto to_yyyyz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 338);

        auto to_yyyyz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 339);

        auto to_yyyyz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 340);

        auto to_yyyyz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 341);

        auto to_yyyyz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 342);

        auto to_yyyyz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 343);

        auto to_yyyyz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 344);

        auto to_yyyyz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 345);

        auto to_yyyyz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 346);

        auto to_yyyyz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 347);

        auto to_yyyyz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 348);

        auto to_yyyyz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 349);

        auto to_yyyyz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 350);

        auto to_yyyyz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 351);

        auto to_yyyyz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 352);

        auto to_yyyyz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 353);

        auto to_yyyyz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 354);

        auto to_yyyyz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 355);

        auto to_yyyyz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 356);

        auto to_yyyzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 357);

        auto to_yyyzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 358);

        auto to_yyyzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 359);

        auto to_yyyzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 360);

        auto to_yyyzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 361);

        auto to_yyyzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 362);

        auto to_yyyzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 363);

        auto to_yyyzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 364);

        auto to_yyyzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 365);

        auto to_yyyzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 366);

        auto to_yyyzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 367);

        auto to_yyyzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 368);

        auto to_yyyzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 369);

        auto to_yyyzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 370);

        auto to_yyyzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 371);

        auto to_yyyzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 372);

        auto to_yyyzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 373);

        auto to_yyyzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 374);

        auto to_yyyzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 375);

        auto to_yyyzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 376);

        auto to_yyyzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 377);

        auto to_yyzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 378);

        auto to_yyzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 379);

        auto to_yyzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 380);

        auto to_yyzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 381);

        auto to_yyzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 382);

        auto to_yyzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 383);

        auto to_yyzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 384);

        auto to_yyzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 385);

        auto to_yyzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 386);

        auto to_yyzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 387);

        auto to_yyzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 388);

        auto to_yyzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 389);

        auto to_yyzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 390);

        auto to_yyzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 391);

        auto to_yyzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 392);

        auto to_yyzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 393);

        auto to_yyzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 394);

        auto to_yyzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 395);

        auto to_yyzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 396);

        auto to_yyzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 397);

        auto to_yyzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 398);

        auto to_yzzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 399);

        auto to_yzzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 400);

        auto to_yzzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 401);

        auto to_yzzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 402);

        auto to_yzzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 403);

        auto to_yzzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 404);

        auto to_yzzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 405);

        auto to_yzzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 406);

        auto to_yzzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 407);

        auto to_yzzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 408);

        auto to_yzzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 409);

        auto to_yzzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 410);

        auto to_yzzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 411);

        auto to_yzzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 412);

        auto to_yzzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 413);

        auto to_yzzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 414);

        auto to_yzzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 415);

        auto to_yzzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 416);

        auto to_yzzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 417);

        auto to_yzzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 418);

        auto to_yzzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 419);

        auto to_zzzzz_xxxxx = pbuffer.data(idx_op_hh + i * 441 + 420);

        auto to_zzzzz_xxxxy = pbuffer.data(idx_op_hh + i * 441 + 421);

        auto to_zzzzz_xxxxz = pbuffer.data(idx_op_hh + i * 441 + 422);

        auto to_zzzzz_xxxyy = pbuffer.data(idx_op_hh + i * 441 + 423);

        auto to_zzzzz_xxxyz = pbuffer.data(idx_op_hh + i * 441 + 424);

        auto to_zzzzz_xxxzz = pbuffer.data(idx_op_hh + i * 441 + 425);

        auto to_zzzzz_xxyyy = pbuffer.data(idx_op_hh + i * 441 + 426);

        auto to_zzzzz_xxyyz = pbuffer.data(idx_op_hh + i * 441 + 427);

        auto to_zzzzz_xxyzz = pbuffer.data(idx_op_hh + i * 441 + 428);

        auto to_zzzzz_xxzzz = pbuffer.data(idx_op_hh + i * 441 + 429);

        auto to_zzzzz_xyyyy = pbuffer.data(idx_op_hh + i * 441 + 430);

        auto to_zzzzz_xyyyz = pbuffer.data(idx_op_hh + i * 441 + 431);

        auto to_zzzzz_xyyzz = pbuffer.data(idx_op_hh + i * 441 + 432);

        auto to_zzzzz_xyzzz = pbuffer.data(idx_op_hh + i * 441 + 433);

        auto to_zzzzz_xzzzz = pbuffer.data(idx_op_hh + i * 441 + 434);

        auto to_zzzzz_yyyyy = pbuffer.data(idx_op_hh + i * 441 + 435);

        auto to_zzzzz_yyyyz = pbuffer.data(idx_op_hh + i * 441 + 436);

        auto to_zzzzz_yyyzz = pbuffer.data(idx_op_hh + i * 441 + 437);

        auto to_zzzzz_yyzzz = pbuffer.data(idx_op_hh + i * 441 + 438);

        auto to_zzzzz_yzzzz = pbuffer.data(idx_op_hh + i * 441 + 439);

        auto to_zzzzz_zzzzz = pbuffer.data(idx_op_hh + i * 441 + 440);

        // Set up 0-15 components of targeted buffer : GG

        auto to_x_x_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 0);

        auto to_x_x_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 1);

        auto to_x_x_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 2);

        auto to_x_x_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 3);

        auto to_x_x_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 4);

        auto to_x_x_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 5);

        auto to_x_x_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 6);

        auto to_x_x_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 7);

        auto to_x_x_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 8);

        auto to_x_x_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 9);

        auto to_x_x_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 10);

        auto to_x_x_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 11);

        auto to_x_x_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 12);

        auto to_x_x_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 13);

        auto to_x_x_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_x_x_xxxx_xxxx, to_x_x_xxxx_xxxy, to_x_x_xxxx_xxxz, to_x_x_xxxx_xxyy, to_x_x_xxxx_xxyz, to_x_x_xxxx_xxzz, to_x_x_xxxx_xyyy, to_x_x_xxxx_xyyz, to_x_x_xxxx_xyzz, to_x_x_xxxx_xzzz, to_x_x_xxxx_yyyy, to_x_x_xxxx_yyyz, to_x_x_xxxx_yyzz, to_x_x_xxxx_yzzz, to_x_x_xxxx_zzzz, to_xxx_xxx, to_xxx_xxxxx, to_xxx_xxxxy, to_xxx_xxxxz, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_zzz, to_xxxxx_xxx, to_xxxxx_xxxxx, to_xxxxx_xxxxy, to_xxxxx_xxxxz, to_xxxxx_xxxyy, to_xxxxx_xxxyz, to_xxxxx_xxxzz, to_xxxxx_xxy, to_xxxxx_xxyyy, to_xxxxx_xxyyz, to_xxxxx_xxyzz, to_xxxxx_xxz, to_xxxxx_xxzzz, to_xxxxx_xyy, to_xxxxx_xyyyy, to_xxxxx_xyyyz, to_xxxxx_xyyzz, to_xxxxx_xyz, to_xxxxx_xyzzz, to_xxxxx_xzz, to_xxxxx_xzzzz, to_xxxxx_yyy, to_xxxxx_yyz, to_xxxxx_yzz, to_xxxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxx_xxxx[k] = 16.0 * to_xxx_xxx[k] - 8.0 * to_xxx_xxxxx[k] * tke_0 - 8.0 * to_xxxxx_xxx[k] * tbe_0 + 4.0 * to_xxxxx_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxxy[k] = 12.0 * to_xxx_xxy[k] - 8.0 * to_xxx_xxxxy[k] * tke_0 - 6.0 * to_xxxxx_xxy[k] * tbe_0 + 4.0 * to_xxxxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxxz[k] = 12.0 * to_xxx_xxz[k] - 8.0 * to_xxx_xxxxz[k] * tke_0 - 6.0 * to_xxxxx_xxz[k] * tbe_0 + 4.0 * to_xxxxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxyy[k] = 8.0 * to_xxx_xyy[k] - 8.0 * to_xxx_xxxyy[k] * tke_0 - 4.0 * to_xxxxx_xyy[k] * tbe_0 + 4.0 * to_xxxxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxyz[k] = 8.0 * to_xxx_xyz[k] - 8.0 * to_xxx_xxxyz[k] * tke_0 - 4.0 * to_xxxxx_xyz[k] * tbe_0 + 4.0 * to_xxxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxzz[k] = 8.0 * to_xxx_xzz[k] - 8.0 * to_xxx_xxxzz[k] * tke_0 - 4.0 * to_xxxxx_xzz[k] * tbe_0 + 4.0 * to_xxxxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xyyy[k] = 4.0 * to_xxx_yyy[k] - 8.0 * to_xxx_xxyyy[k] * tke_0 - 2.0 * to_xxxxx_yyy[k] * tbe_0 + 4.0 * to_xxxxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xyyz[k] = 4.0 * to_xxx_yyz[k] - 8.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxxx_yyz[k] * tbe_0 + 4.0 * to_xxxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xyzz[k] = 4.0 * to_xxx_yzz[k] - 8.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxxx_yzz[k] * tbe_0 + 4.0 * to_xxxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xzzz[k] = 4.0 * to_xxx_zzz[k] - 8.0 * to_xxx_xxzzz[k] * tke_0 - 2.0 * to_xxxxx_zzz[k] * tbe_0 + 4.0 * to_xxxxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yyyy[k] = -8.0 * to_xxx_xyyyy[k] * tke_0 + 4.0 * to_xxxxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yyyz[k] = -8.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yyzz[k] = -8.0 * to_xxx_xyyzz[k] * tke_0 + 4.0 * to_xxxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yzzz[k] = -8.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_zzzz[k] = -8.0 * to_xxx_xzzzz[k] * tke_0 + 4.0 * to_xxxxx_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : GG

        auto to_x_x_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 15);

        auto to_x_x_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 16);

        auto to_x_x_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 17);

        auto to_x_x_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 18);

        auto to_x_x_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 19);

        auto to_x_x_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 20);

        auto to_x_x_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 21);

        auto to_x_x_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 22);

        auto to_x_x_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 23);

        auto to_x_x_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 24);

        auto to_x_x_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 25);

        auto to_x_x_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 26);

        auto to_x_x_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 27);

        auto to_x_x_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 28);

        auto to_x_x_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_x_x_xxxy_xxxx, to_x_x_xxxy_xxxy, to_x_x_xxxy_xxxz, to_x_x_xxxy_xxyy, to_x_x_xxxy_xxyz, to_x_x_xxxy_xxzz, to_x_x_xxxy_xyyy, to_x_x_xxxy_xyyz, to_x_x_xxxy_xyzz, to_x_x_xxxy_xzzz, to_x_x_xxxy_yyyy, to_x_x_xxxy_yyyz, to_x_x_xxxy_yyzz, to_x_x_xxxy_yzzz, to_x_x_xxxy_zzzz, to_xxxxy_xxx, to_xxxxy_xxxxx, to_xxxxy_xxxxy, to_xxxxy_xxxxz, to_xxxxy_xxxyy, to_xxxxy_xxxyz, to_xxxxy_xxxzz, to_xxxxy_xxy, to_xxxxy_xxyyy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xxzzz, to_xxxxy_xyy, to_xxxxy_xyyyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_xzzzz, to_xxxxy_yyy, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_zzz, to_xxy_xxx, to_xxy_xxxxx, to_xxy_xxxxy, to_xxy_xxxxz, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxy_xxxx[k] = 12.0 * to_xxy_xxx[k] - 6.0 * to_xxy_xxxxx[k] * tke_0 - 8.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxxy[k] = 9.0 * to_xxy_xxy[k] - 6.0 * to_xxy_xxxxy[k] * tke_0 - 6.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxxz[k] = 9.0 * to_xxy_xxz[k] - 6.0 * to_xxy_xxxxz[k] * tke_0 - 6.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxyy[k] = 6.0 * to_xxy_xyy[k] - 6.0 * to_xxy_xxxyy[k] * tke_0 - 4.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxyz[k] = 6.0 * to_xxy_xyz[k] - 6.0 * to_xxy_xxxyz[k] * tke_0 - 4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxzz[k] = 6.0 * to_xxy_xzz[k] - 6.0 * to_xxy_xxxzz[k] * tke_0 - 4.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xyyy[k] = 3.0 * to_xxy_yyy[k] - 6.0 * to_xxy_xxyyy[k] * tke_0 - 2.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xyyz[k] = 3.0 * to_xxy_yyz[k] - 6.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xyzz[k] = 3.0 * to_xxy_yzz[k] - 6.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xzzz[k] = 3.0 * to_xxy_zzz[k] - 6.0 * to_xxy_xxzzz[k] * tke_0 - 2.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yyyy[k] = -6.0 * to_xxy_xyyyy[k] * tke_0 + 4.0 * to_xxxxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yyyz[k] = -6.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yyzz[k] = -6.0 * to_xxy_xyyzz[k] * tke_0 + 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yzzz[k] = -6.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_zzzz[k] = -6.0 * to_xxy_xzzzz[k] * tke_0 + 4.0 * to_xxxxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : GG

        auto to_x_x_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 30);

        auto to_x_x_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 31);

        auto to_x_x_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 32);

        auto to_x_x_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 33);

        auto to_x_x_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 34);

        auto to_x_x_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 35);

        auto to_x_x_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 36);

        auto to_x_x_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 37);

        auto to_x_x_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 38);

        auto to_x_x_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 39);

        auto to_x_x_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 40);

        auto to_x_x_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 41);

        auto to_x_x_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 42);

        auto to_x_x_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 43);

        auto to_x_x_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_x_x_xxxz_xxxx, to_x_x_xxxz_xxxy, to_x_x_xxxz_xxxz, to_x_x_xxxz_xxyy, to_x_x_xxxz_xxyz, to_x_x_xxxz_xxzz, to_x_x_xxxz_xyyy, to_x_x_xxxz_xyyz, to_x_x_xxxz_xyzz, to_x_x_xxxz_xzzz, to_x_x_xxxz_yyyy, to_x_x_xxxz_yyyz, to_x_x_xxxz_yyzz, to_x_x_xxxz_yzzz, to_x_x_xxxz_zzzz, to_xxxxz_xxx, to_xxxxz_xxxxx, to_xxxxz_xxxxy, to_xxxxz_xxxxz, to_xxxxz_xxxyy, to_xxxxz_xxxyz, to_xxxxz_xxxzz, to_xxxxz_xxy, to_xxxxz_xxyyy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xxzzz, to_xxxxz_xyy, to_xxxxz_xyyyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_xzzzz, to_xxxxz_yyy, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_zzz, to_xxz_xxx, to_xxz_xxxxx, to_xxz_xxxxy, to_xxz_xxxxz, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxz_xxxx[k] = 12.0 * to_xxz_xxx[k] - 6.0 * to_xxz_xxxxx[k] * tke_0 - 8.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxxy[k] = 9.0 * to_xxz_xxy[k] - 6.0 * to_xxz_xxxxy[k] * tke_0 - 6.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxxz[k] = 9.0 * to_xxz_xxz[k] - 6.0 * to_xxz_xxxxz[k] * tke_0 - 6.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxyy[k] = 6.0 * to_xxz_xyy[k] - 6.0 * to_xxz_xxxyy[k] * tke_0 - 4.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxyz[k] = 6.0 * to_xxz_xyz[k] - 6.0 * to_xxz_xxxyz[k] * tke_0 - 4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxzz[k] = 6.0 * to_xxz_xzz[k] - 6.0 * to_xxz_xxxzz[k] * tke_0 - 4.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xyyy[k] = 3.0 * to_xxz_yyy[k] - 6.0 * to_xxz_xxyyy[k] * tke_0 - 2.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xyyz[k] = 3.0 * to_xxz_yyz[k] - 6.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xyzz[k] = 3.0 * to_xxz_yzz[k] - 6.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xzzz[k] = 3.0 * to_xxz_zzz[k] - 6.0 * to_xxz_xxzzz[k] * tke_0 - 2.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yyyy[k] = -6.0 * to_xxz_xyyyy[k] * tke_0 + 4.0 * to_xxxxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yyyz[k] = -6.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yyzz[k] = -6.0 * to_xxz_xyyzz[k] * tke_0 + 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yzzz[k] = -6.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_zzzz[k] = -6.0 * to_xxz_xzzzz[k] * tke_0 + 4.0 * to_xxxxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : GG

        auto to_x_x_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 45);

        auto to_x_x_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 46);

        auto to_x_x_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 47);

        auto to_x_x_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 48);

        auto to_x_x_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 49);

        auto to_x_x_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 50);

        auto to_x_x_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 51);

        auto to_x_x_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 52);

        auto to_x_x_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 53);

        auto to_x_x_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 54);

        auto to_x_x_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 55);

        auto to_x_x_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 56);

        auto to_x_x_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 57);

        auto to_x_x_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 58);

        auto to_x_x_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_x_x_xxyy_xxxx, to_x_x_xxyy_xxxy, to_x_x_xxyy_xxxz, to_x_x_xxyy_xxyy, to_x_x_xxyy_xxyz, to_x_x_xxyy_xxzz, to_x_x_xxyy_xyyy, to_x_x_xxyy_xyyz, to_x_x_xxyy_xyzz, to_x_x_xxyy_xzzz, to_x_x_xxyy_yyyy, to_x_x_xxyy_yyyz, to_x_x_xxyy_yyzz, to_x_x_xxyy_yzzz, to_x_x_xxyy_zzzz, to_xxxyy_xxx, to_xxxyy_xxxxx, to_xxxyy_xxxxy, to_xxxyy_xxxxz, to_xxxyy_xxxyy, to_xxxyy_xxxyz, to_xxxyy_xxxzz, to_xxxyy_xxy, to_xxxyy_xxyyy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xxzzz, to_xxxyy_xyy, to_xxxyy_xyyyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_xzzzz, to_xxxyy_yyy, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_zzz, to_xyy_xxx, to_xyy_xxxxx, to_xyy_xxxxy, to_xyy_xxxxz, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyy_xxxx[k] = 8.0 * to_xyy_xxx[k] - 4.0 * to_xyy_xxxxx[k] * tke_0 - 8.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxxy[k] = 6.0 * to_xyy_xxy[k] - 4.0 * to_xyy_xxxxy[k] * tke_0 - 6.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxxz[k] = 6.0 * to_xyy_xxz[k] - 4.0 * to_xyy_xxxxz[k] * tke_0 - 6.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxyy[k] = 4.0 * to_xyy_xyy[k] - 4.0 * to_xyy_xxxyy[k] * tke_0 - 4.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxyz[k] = 4.0 * to_xyy_xyz[k] - 4.0 * to_xyy_xxxyz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxzz[k] = 4.0 * to_xyy_xzz[k] - 4.0 * to_xyy_xxxzz[k] * tke_0 - 4.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xyyy[k] = 2.0 * to_xyy_yyy[k] - 4.0 * to_xyy_xxyyy[k] * tke_0 - 2.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xyyz[k] = 2.0 * to_xyy_yyz[k] - 4.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xyzz[k] = 2.0 * to_xyy_yzz[k] - 4.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xzzz[k] = 2.0 * to_xyy_zzz[k] - 4.0 * to_xyy_xxzzz[k] * tke_0 - 2.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yyyy[k] = -4.0 * to_xyy_xyyyy[k] * tke_0 + 4.0 * to_xxxyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yyyz[k] = -4.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yyzz[k] = -4.0 * to_xyy_xyyzz[k] * tke_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yzzz[k] = -4.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_zzzz[k] = -4.0 * to_xyy_xzzzz[k] * tke_0 + 4.0 * to_xxxyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : GG

        auto to_x_x_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 60);

        auto to_x_x_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 61);

        auto to_x_x_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 62);

        auto to_x_x_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 63);

        auto to_x_x_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 64);

        auto to_x_x_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 65);

        auto to_x_x_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 66);

        auto to_x_x_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 67);

        auto to_x_x_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 68);

        auto to_x_x_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 69);

        auto to_x_x_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 70);

        auto to_x_x_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 71);

        auto to_x_x_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 72);

        auto to_x_x_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 73);

        auto to_x_x_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_x_x_xxyz_xxxx, to_x_x_xxyz_xxxy, to_x_x_xxyz_xxxz, to_x_x_xxyz_xxyy, to_x_x_xxyz_xxyz, to_x_x_xxyz_xxzz, to_x_x_xxyz_xyyy, to_x_x_xxyz_xyyz, to_x_x_xxyz_xyzz, to_x_x_xxyz_xzzz, to_x_x_xxyz_yyyy, to_x_x_xxyz_yyyz, to_x_x_xxyz_yyzz, to_x_x_xxyz_yzzz, to_x_x_xxyz_zzzz, to_xxxyz_xxx, to_xxxyz_xxxxx, to_xxxyz_xxxxy, to_xxxyz_xxxxz, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_zzz, to_xyz_xxx, to_xyz_xxxxx, to_xyz_xxxxy, to_xyz_xxxxz, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyz_xxxx[k] = 8.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxxx[k] * tke_0 - 8.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxxy[k] = 6.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxxxy[k] * tke_0 - 6.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxxz[k] = 6.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxxxz[k] * tke_0 - 6.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxyy[k] = 4.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 4.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xxxyz[k] * tke_0 - 4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxzz[k] = 4.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 4.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xyyy[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xyyz[k] = 2.0 * to_xyz_yyz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xyzz[k] = 2.0 * to_xyz_yzz[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yyyy[k] = -4.0 * to_xyz_xyyyy[k] * tke_0 + 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yyyz[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yyzz[k] = -4.0 * to_xyz_xyyzz[k] * tke_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_zzzz[k] = -4.0 * to_xyz_xzzzz[k] * tke_0 + 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : GG

        auto to_x_x_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 75);

        auto to_x_x_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 76);

        auto to_x_x_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 77);

        auto to_x_x_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 78);

        auto to_x_x_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 79);

        auto to_x_x_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 80);

        auto to_x_x_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 81);

        auto to_x_x_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 82);

        auto to_x_x_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 83);

        auto to_x_x_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 84);

        auto to_x_x_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 85);

        auto to_x_x_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 86);

        auto to_x_x_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 87);

        auto to_x_x_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 88);

        auto to_x_x_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_x_x_xxzz_xxxx, to_x_x_xxzz_xxxy, to_x_x_xxzz_xxxz, to_x_x_xxzz_xxyy, to_x_x_xxzz_xxyz, to_x_x_xxzz_xxzz, to_x_x_xxzz_xyyy, to_x_x_xxzz_xyyz, to_x_x_xxzz_xyzz, to_x_x_xxzz_xzzz, to_x_x_xxzz_yyyy, to_x_x_xxzz_yyyz, to_x_x_xxzz_yyzz, to_x_x_xxzz_yzzz, to_x_x_xxzz_zzzz, to_xxxzz_xxx, to_xxxzz_xxxxx, to_xxxzz_xxxxy, to_xxxzz_xxxxz, to_xxxzz_xxxyy, to_xxxzz_xxxyz, to_xxxzz_xxxzz, to_xxxzz_xxy, to_xxxzz_xxyyy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xxzzz, to_xxxzz_xyy, to_xxxzz_xyyyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_xzzzz, to_xxxzz_yyy, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_zzz, to_xzz_xxx, to_xzz_xxxxx, to_xzz_xxxxy, to_xzz_xxxxz, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxzz_xxxx[k] = 8.0 * to_xzz_xxx[k] - 4.0 * to_xzz_xxxxx[k] * tke_0 - 8.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxxy[k] = 6.0 * to_xzz_xxy[k] - 4.0 * to_xzz_xxxxy[k] * tke_0 - 6.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxxz[k] = 6.0 * to_xzz_xxz[k] - 4.0 * to_xzz_xxxxz[k] * tke_0 - 6.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxyy[k] = 4.0 * to_xzz_xyy[k] - 4.0 * to_xzz_xxxyy[k] * tke_0 - 4.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxyz[k] = 4.0 * to_xzz_xyz[k] - 4.0 * to_xzz_xxxyz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxzz[k] = 4.0 * to_xzz_xzz[k] - 4.0 * to_xzz_xxxzz[k] * tke_0 - 4.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xyyy[k] = 2.0 * to_xzz_yyy[k] - 4.0 * to_xzz_xxyyy[k] * tke_0 - 2.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xyyz[k] = 2.0 * to_xzz_yyz[k] - 4.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xyzz[k] = 2.0 * to_xzz_yzz[k] - 4.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xzzz[k] = 2.0 * to_xzz_zzz[k] - 4.0 * to_xzz_xxzzz[k] * tke_0 - 2.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yyyy[k] = -4.0 * to_xzz_xyyyy[k] * tke_0 + 4.0 * to_xxxzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yyyz[k] = -4.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yyzz[k] = -4.0 * to_xzz_xyyzz[k] * tke_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yzzz[k] = -4.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_zzzz[k] = -4.0 * to_xzz_xzzzz[k] * tke_0 + 4.0 * to_xxxzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : GG

        auto to_x_x_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 90);

        auto to_x_x_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 91);

        auto to_x_x_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 92);

        auto to_x_x_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 93);

        auto to_x_x_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 94);

        auto to_x_x_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 95);

        auto to_x_x_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 96);

        auto to_x_x_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 97);

        auto to_x_x_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 98);

        auto to_x_x_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 99);

        auto to_x_x_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 100);

        auto to_x_x_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 101);

        auto to_x_x_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 102);

        auto to_x_x_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 103);

        auto to_x_x_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_x_x_xyyy_xxxx, to_x_x_xyyy_xxxy, to_x_x_xyyy_xxxz, to_x_x_xyyy_xxyy, to_x_x_xyyy_xxyz, to_x_x_xyyy_xxzz, to_x_x_xyyy_xyyy, to_x_x_xyyy_xyyz, to_x_x_xyyy_xyzz, to_x_x_xyyy_xzzz, to_x_x_xyyy_yyyy, to_x_x_xyyy_yyyz, to_x_x_xyyy_yyzz, to_x_x_xyyy_yzzz, to_x_x_xyyy_zzzz, to_xxyyy_xxx, to_xxyyy_xxxxx, to_xxyyy_xxxxy, to_xxyyy_xxxxz, to_xxyyy_xxxyy, to_xxyyy_xxxyz, to_xxyyy_xxxzz, to_xxyyy_xxy, to_xxyyy_xxyyy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xxzzz, to_xxyyy_xyy, to_xxyyy_xyyyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_xzzzz, to_xxyyy_yyy, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_zzz, to_yyy_xxx, to_yyy_xxxxx, to_yyy_xxxxy, to_yyy_xxxxz, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyy_xxxx[k] = 4.0 * to_yyy_xxx[k] - 2.0 * to_yyy_xxxxx[k] * tke_0 - 8.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxxy[k] = 3.0 * to_yyy_xxy[k] - 2.0 * to_yyy_xxxxy[k] * tke_0 - 6.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxxz[k] = 3.0 * to_yyy_xxz[k] - 2.0 * to_yyy_xxxxz[k] * tke_0 - 6.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxyy[k] = 2.0 * to_yyy_xyy[k] - 2.0 * to_yyy_xxxyy[k] * tke_0 - 4.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxyz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xxxyz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxzz[k] = 2.0 * to_yyy_xzz[k] - 2.0 * to_yyy_xxxzz[k] * tke_0 - 4.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xyyy[k] = to_yyy_yyy[k] - 2.0 * to_yyy_xxyyy[k] * tke_0 - 2.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xyyz[k] = to_yyy_yyz[k] - 2.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xyzz[k] = to_yyy_yzz[k] - 2.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xzzz[k] = to_yyy_zzz[k] - 2.0 * to_yyy_xxzzz[k] * tke_0 - 2.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yyyy[k] = -2.0 * to_yyy_xyyyy[k] * tke_0 + 4.0 * to_xxyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yyyz[k] = -2.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yyzz[k] = -2.0 * to_yyy_xyyzz[k] * tke_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yzzz[k] = -2.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_zzzz[k] = -2.0 * to_yyy_xzzzz[k] * tke_0 + 4.0 * to_xxyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : GG

        auto to_x_x_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 105);

        auto to_x_x_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 106);

        auto to_x_x_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 107);

        auto to_x_x_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 108);

        auto to_x_x_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 109);

        auto to_x_x_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 110);

        auto to_x_x_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 111);

        auto to_x_x_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 112);

        auto to_x_x_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 113);

        auto to_x_x_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 114);

        auto to_x_x_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 115);

        auto to_x_x_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 116);

        auto to_x_x_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 117);

        auto to_x_x_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 118);

        auto to_x_x_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_x_x_xyyz_xxxx, to_x_x_xyyz_xxxy, to_x_x_xyyz_xxxz, to_x_x_xyyz_xxyy, to_x_x_xyyz_xxyz, to_x_x_xyyz_xxzz, to_x_x_xyyz_xyyy, to_x_x_xyyz_xyyz, to_x_x_xyyz_xyzz, to_x_x_xyyz_xzzz, to_x_x_xyyz_yyyy, to_x_x_xyyz_yyyz, to_x_x_xyyz_yyzz, to_x_x_xyyz_yzzz, to_x_x_xyyz_zzzz, to_xxyyz_xxx, to_xxyyz_xxxxx, to_xxyyz_xxxxy, to_xxyyz_xxxxz, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_zzz, to_yyz_xxx, to_yyz_xxxxx, to_yyz_xxxxy, to_yyz_xxxxz, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyz_xxxx[k] = 4.0 * to_yyz_xxx[k] - 2.0 * to_yyz_xxxxx[k] * tke_0 - 8.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxxy[k] = 3.0 * to_yyz_xxy[k] - 2.0 * to_yyz_xxxxy[k] * tke_0 - 6.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxxz[k] = 3.0 * to_yyz_xxz[k] - 2.0 * to_yyz_xxxxz[k] * tke_0 - 6.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxyy[k] = 2.0 * to_yyz_xyy[k] - 2.0 * to_yyz_xxxyy[k] * tke_0 - 4.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxyz[k] = 2.0 * to_yyz_xyz[k] - 2.0 * to_yyz_xxxyz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxzz[k] = 2.0 * to_yyz_xzz[k] - 2.0 * to_yyz_xxxzz[k] * tke_0 - 4.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xyyy[k] = to_yyz_yyy[k] - 2.0 * to_yyz_xxyyy[k] * tke_0 - 2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xyyz[k] = to_yyz_yyz[k] - 2.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xyzz[k] = to_yyz_yzz[k] - 2.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xzzz[k] = to_yyz_zzz[k] - 2.0 * to_yyz_xxzzz[k] * tke_0 - 2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yyyy[k] = -2.0 * to_yyz_xyyyy[k] * tke_0 + 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yyyz[k] = -2.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yyzz[k] = -2.0 * to_yyz_xyyzz[k] * tke_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yzzz[k] = -2.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_zzzz[k] = -2.0 * to_yyz_xzzzz[k] * tke_0 + 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : GG

        auto to_x_x_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 120);

        auto to_x_x_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 121);

        auto to_x_x_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 122);

        auto to_x_x_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 123);

        auto to_x_x_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 124);

        auto to_x_x_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 125);

        auto to_x_x_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 126);

        auto to_x_x_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 127);

        auto to_x_x_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 128);

        auto to_x_x_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 129);

        auto to_x_x_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 130);

        auto to_x_x_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 131);

        auto to_x_x_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 132);

        auto to_x_x_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 133);

        auto to_x_x_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_x_x_xyzz_xxxx, to_x_x_xyzz_xxxy, to_x_x_xyzz_xxxz, to_x_x_xyzz_xxyy, to_x_x_xyzz_xxyz, to_x_x_xyzz_xxzz, to_x_x_xyzz_xyyy, to_x_x_xyzz_xyyz, to_x_x_xyzz_xyzz, to_x_x_xyzz_xzzz, to_x_x_xyzz_yyyy, to_x_x_xyzz_yyyz, to_x_x_xyzz_yyzz, to_x_x_xyzz_yzzz, to_x_x_xyzz_zzzz, to_xxyzz_xxx, to_xxyzz_xxxxx, to_xxyzz_xxxxy, to_xxyzz_xxxxz, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_zzz, to_yzz_xxx, to_yzz_xxxxx, to_yzz_xxxxy, to_yzz_xxxxz, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyzz_xxxx[k] = 4.0 * to_yzz_xxx[k] - 2.0 * to_yzz_xxxxx[k] * tke_0 - 8.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxxy[k] = 3.0 * to_yzz_xxy[k] - 2.0 * to_yzz_xxxxy[k] * tke_0 - 6.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxxz[k] = 3.0 * to_yzz_xxz[k] - 2.0 * to_yzz_xxxxz[k] * tke_0 - 6.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxyy[k] = 2.0 * to_yzz_xyy[k] - 2.0 * to_yzz_xxxyy[k] * tke_0 - 4.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxyz[k] = 2.0 * to_yzz_xyz[k] - 2.0 * to_yzz_xxxyz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxzz[k] = 2.0 * to_yzz_xzz[k] - 2.0 * to_yzz_xxxzz[k] * tke_0 - 4.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xyyy[k] = to_yzz_yyy[k] - 2.0 * to_yzz_xxyyy[k] * tke_0 - 2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xyyz[k] = to_yzz_yyz[k] - 2.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xyzz[k] = to_yzz_yzz[k] - 2.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xzzz[k] = to_yzz_zzz[k] - 2.0 * to_yzz_xxzzz[k] * tke_0 - 2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yyyy[k] = -2.0 * to_yzz_xyyyy[k] * tke_0 + 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yyyz[k] = -2.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yyzz[k] = -2.0 * to_yzz_xyyzz[k] * tke_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yzzz[k] = -2.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_zzzz[k] = -2.0 * to_yzz_xzzzz[k] * tke_0 + 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 135-150 components of targeted buffer : GG

        auto to_x_x_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 135);

        auto to_x_x_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 136);

        auto to_x_x_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 137);

        auto to_x_x_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 138);

        auto to_x_x_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 139);

        auto to_x_x_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 140);

        auto to_x_x_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 141);

        auto to_x_x_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 142);

        auto to_x_x_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 143);

        auto to_x_x_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 144);

        auto to_x_x_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 145);

        auto to_x_x_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 146);

        auto to_x_x_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 147);

        auto to_x_x_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 148);

        auto to_x_x_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_x_x_xzzz_xxxx, to_x_x_xzzz_xxxy, to_x_x_xzzz_xxxz, to_x_x_xzzz_xxyy, to_x_x_xzzz_xxyz, to_x_x_xzzz_xxzz, to_x_x_xzzz_xyyy, to_x_x_xzzz_xyyz, to_x_x_xzzz_xyzz, to_x_x_xzzz_xzzz, to_x_x_xzzz_yyyy, to_x_x_xzzz_yyyz, to_x_x_xzzz_yyzz, to_x_x_xzzz_yzzz, to_x_x_xzzz_zzzz, to_xxzzz_xxx, to_xxzzz_xxxxx, to_xxzzz_xxxxy, to_xxzzz_xxxxz, to_xxzzz_xxxyy, to_xxzzz_xxxyz, to_xxzzz_xxxzz, to_xxzzz_xxy, to_xxzzz_xxyyy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xxzzz, to_xxzzz_xyy, to_xxzzz_xyyyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_xzzzz, to_xxzzz_yyy, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_zzz, to_zzz_xxx, to_zzz_xxxxx, to_zzz_xxxxy, to_zzz_xxxxz, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzzz_xxxx[k] = 4.0 * to_zzz_xxx[k] - 2.0 * to_zzz_xxxxx[k] * tke_0 - 8.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxxy[k] = 3.0 * to_zzz_xxy[k] - 2.0 * to_zzz_xxxxy[k] * tke_0 - 6.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxxz[k] = 3.0 * to_zzz_xxz[k] - 2.0 * to_zzz_xxxxz[k] * tke_0 - 6.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxyy[k] = 2.0 * to_zzz_xyy[k] - 2.0 * to_zzz_xxxyy[k] * tke_0 - 4.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxyz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xxxyz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxzz[k] = 2.0 * to_zzz_xzz[k] - 2.0 * to_zzz_xxxzz[k] * tke_0 - 4.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xyyy[k] = to_zzz_yyy[k] - 2.0 * to_zzz_xxyyy[k] * tke_0 - 2.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xyyz[k] = to_zzz_yyz[k] - 2.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xyzz[k] = to_zzz_yzz[k] - 2.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xzzz[k] = to_zzz_zzz[k] - 2.0 * to_zzz_xxzzz[k] * tke_0 - 2.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yyyy[k] = -2.0 * to_zzz_xyyyy[k] * tke_0 + 4.0 * to_xxzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yyyz[k] = -2.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yyzz[k] = -2.0 * to_zzz_xyyzz[k] * tke_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yzzz[k] = -2.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_zzzz[k] = -2.0 * to_zzz_xzzzz[k] * tke_0 + 4.0 * to_xxzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-165 components of targeted buffer : GG

        auto to_x_x_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 150);

        auto to_x_x_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 151);

        auto to_x_x_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 152);

        auto to_x_x_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 153);

        auto to_x_x_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 154);

        auto to_x_x_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 155);

        auto to_x_x_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 156);

        auto to_x_x_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 157);

        auto to_x_x_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 158);

        auto to_x_x_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 159);

        auto to_x_x_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 160);

        auto to_x_x_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 161);

        auto to_x_x_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 162);

        auto to_x_x_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 163);

        auto to_x_x_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_x_x_yyyy_xxxx, to_x_x_yyyy_xxxy, to_x_x_yyyy_xxxz, to_x_x_yyyy_xxyy, to_x_x_yyyy_xxyz, to_x_x_yyyy_xxzz, to_x_x_yyyy_xyyy, to_x_x_yyyy_xyyz, to_x_x_yyyy_xyzz, to_x_x_yyyy_xzzz, to_x_x_yyyy_yyyy, to_x_x_yyyy_yyyz, to_x_x_yyyy_yyzz, to_x_x_yyyy_yzzz, to_x_x_yyyy_zzzz, to_xyyyy_xxx, to_xyyyy_xxxxx, to_xyyyy_xxxxy, to_xyyyy_xxxxz, to_xyyyy_xxxyy, to_xyyyy_xxxyz, to_xyyyy_xxxzz, to_xyyyy_xxy, to_xyyyy_xxyyy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xxzzz, to_xyyyy_xyy, to_xyyyy_xyyyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_xzzzz, to_xyyyy_yyy, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyy_xxxx[k] = -8.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxxy[k] = -6.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxxz[k] = -6.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxyy[k] = -4.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxyz[k] = -4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxzz[k] = -4.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xyyy[k] = -2.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xyyz[k] = -2.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xyzz[k] = -2.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xzzz[k] = -2.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yyyy[k] = 4.0 * to_xyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yyyz[k] = 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yyzz[k] = 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yzzz[k] = 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_zzzz[k] = 4.0 * to_xyyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 165-180 components of targeted buffer : GG

        auto to_x_x_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 165);

        auto to_x_x_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 166);

        auto to_x_x_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 167);

        auto to_x_x_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 168);

        auto to_x_x_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 169);

        auto to_x_x_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 170);

        auto to_x_x_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 171);

        auto to_x_x_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 172);

        auto to_x_x_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 173);

        auto to_x_x_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 174);

        auto to_x_x_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 175);

        auto to_x_x_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 176);

        auto to_x_x_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 177);

        auto to_x_x_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 178);

        auto to_x_x_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_x_x_yyyz_xxxx, to_x_x_yyyz_xxxy, to_x_x_yyyz_xxxz, to_x_x_yyyz_xxyy, to_x_x_yyyz_xxyz, to_x_x_yyyz_xxzz, to_x_x_yyyz_xyyy, to_x_x_yyyz_xyyz, to_x_x_yyyz_xyzz, to_x_x_yyyz_xzzz, to_x_x_yyyz_yyyy, to_x_x_yyyz_yyyz, to_x_x_yyyz_yyzz, to_x_x_yyyz_yzzz, to_x_x_yyyz_zzzz, to_xyyyz_xxx, to_xyyyz_xxxxx, to_xyyyz_xxxxy, to_xyyyz_xxxxz, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyz_xxxx[k] = -8.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxxy[k] = -6.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxxz[k] = -6.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxyy[k] = -4.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxyz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxzz[k] = -4.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xyyy[k] = -2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xyyz[k] = -2.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xyzz[k] = -2.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xzzz[k] = -2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yyyy[k] = 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yyyz[k] = 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yyzz[k] = 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yzzz[k] = 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_zzzz[k] = 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-195 components of targeted buffer : GG

        auto to_x_x_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 180);

        auto to_x_x_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 181);

        auto to_x_x_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 182);

        auto to_x_x_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 183);

        auto to_x_x_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 184);

        auto to_x_x_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 185);

        auto to_x_x_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 186);

        auto to_x_x_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 187);

        auto to_x_x_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 188);

        auto to_x_x_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 189);

        auto to_x_x_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 190);

        auto to_x_x_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 191);

        auto to_x_x_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 192);

        auto to_x_x_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 193);

        auto to_x_x_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_x_x_yyzz_xxxx, to_x_x_yyzz_xxxy, to_x_x_yyzz_xxxz, to_x_x_yyzz_xxyy, to_x_x_yyzz_xxyz, to_x_x_yyzz_xxzz, to_x_x_yyzz_xyyy, to_x_x_yyzz_xyyz, to_x_x_yyzz_xyzz, to_x_x_yyzz_xzzz, to_x_x_yyzz_yyyy, to_x_x_yyzz_yyyz, to_x_x_yyzz_yyzz, to_x_x_yyzz_yzzz, to_x_x_yyzz_zzzz, to_xyyzz_xxx, to_xyyzz_xxxxx, to_xyyzz_xxxxy, to_xyyzz_xxxxz, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyzz_xxxx[k] = -8.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxxy[k] = -6.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxxz[k] = -6.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxyy[k] = -4.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxyz[k] = -4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxzz[k] = -4.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xyyy[k] = -2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xyyz[k] = -2.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xyzz[k] = -2.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xzzz[k] = -2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yyyy[k] = 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yyyz[k] = 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yyzz[k] = 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yzzz[k] = 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_zzzz[k] = 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 195-210 components of targeted buffer : GG

        auto to_x_x_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 195);

        auto to_x_x_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 196);

        auto to_x_x_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 197);

        auto to_x_x_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 198);

        auto to_x_x_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 199);

        auto to_x_x_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 200);

        auto to_x_x_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 201);

        auto to_x_x_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 202);

        auto to_x_x_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 203);

        auto to_x_x_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 204);

        auto to_x_x_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 205);

        auto to_x_x_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 206);

        auto to_x_x_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 207);

        auto to_x_x_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 208);

        auto to_x_x_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_x_x_yzzz_xxxx, to_x_x_yzzz_xxxy, to_x_x_yzzz_xxxz, to_x_x_yzzz_xxyy, to_x_x_yzzz_xxyz, to_x_x_yzzz_xxzz, to_x_x_yzzz_xyyy, to_x_x_yzzz_xyyz, to_x_x_yzzz_xyzz, to_x_x_yzzz_xzzz, to_x_x_yzzz_yyyy, to_x_x_yzzz_yyyz, to_x_x_yzzz_yyzz, to_x_x_yzzz_yzzz, to_x_x_yzzz_zzzz, to_xyzzz_xxx, to_xyzzz_xxxxx, to_xyzzz_xxxxy, to_xyzzz_xxxxz, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzzz_xxxx[k] = -8.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxxy[k] = -6.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxxz[k] = -6.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxyy[k] = -4.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxyz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxzz[k] = -4.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xyyy[k] = -2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xyyz[k] = -2.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xyzz[k] = -2.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xzzz[k] = -2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yyyy[k] = 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yyyz[k] = 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yyzz[k] = 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yzzz[k] = 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_zzzz[k] = 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-225 components of targeted buffer : GG

        auto to_x_x_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 210);

        auto to_x_x_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 211);

        auto to_x_x_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 212);

        auto to_x_x_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 213);

        auto to_x_x_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 214);

        auto to_x_x_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 215);

        auto to_x_x_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 216);

        auto to_x_x_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 217);

        auto to_x_x_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 218);

        auto to_x_x_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 219);

        auto to_x_x_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 220);

        auto to_x_x_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 221);

        auto to_x_x_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 222);

        auto to_x_x_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 223);

        auto to_x_x_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 0 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_x_x_zzzz_xxxx, to_x_x_zzzz_xxxy, to_x_x_zzzz_xxxz, to_x_x_zzzz_xxyy, to_x_x_zzzz_xxyz, to_x_x_zzzz_xxzz, to_x_x_zzzz_xyyy, to_x_x_zzzz_xyyz, to_x_x_zzzz_xyzz, to_x_x_zzzz_xzzz, to_x_x_zzzz_yyyy, to_x_x_zzzz_yyyz, to_x_x_zzzz_yyzz, to_x_x_zzzz_yzzz, to_x_x_zzzz_zzzz, to_xzzzz_xxx, to_xzzzz_xxxxx, to_xzzzz_xxxxy, to_xzzzz_xxxxz, to_xzzzz_xxxyy, to_xzzzz_xxxyz, to_xzzzz_xxxzz, to_xzzzz_xxy, to_xzzzz_xxyyy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xxzzz, to_xzzzz_xyy, to_xzzzz_xyyyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_xzzzz, to_xzzzz_yyy, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzzz_xxxx[k] = -8.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxxy[k] = -6.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxxz[k] = -6.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxyy[k] = -4.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxyz[k] = -4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxzz[k] = -4.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xyyy[k] = -2.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xyyz[k] = -2.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xyzz[k] = -2.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xzzz[k] = -2.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yyyy[k] = 4.0 * to_xzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yyyz[k] = 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yyzz[k] = 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yzzz[k] = 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_zzzz[k] = 4.0 * to_xzzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 225-240 components of targeted buffer : GG

        auto to_x_y_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 0);

        auto to_x_y_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 1);

        auto to_x_y_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 2);

        auto to_x_y_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 3);

        auto to_x_y_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 4);

        auto to_x_y_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 5);

        auto to_x_y_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 6);

        auto to_x_y_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 7);

        auto to_x_y_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 8);

        auto to_x_y_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 9);

        auto to_x_y_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 10);

        auto to_x_y_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 11);

        auto to_x_y_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 12);

        auto to_x_y_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 13);

        auto to_x_y_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_x_y_xxxx_xxxx, to_x_y_xxxx_xxxy, to_x_y_xxxx_xxxz, to_x_y_xxxx_xxyy, to_x_y_xxxx_xxyz, to_x_y_xxxx_xxzz, to_x_y_xxxx_xyyy, to_x_y_xxxx_xyyz, to_x_y_xxxx_xyzz, to_x_y_xxxx_xzzz, to_x_y_xxxx_yyyy, to_x_y_xxxx_yyyz, to_x_y_xxxx_yyzz, to_x_y_xxxx_yzzz, to_x_y_xxxx_zzzz, to_xxx_xxx, to_xxx_xxxxy, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_yyy, to_xxx_yyyyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxxxx_xxx, to_xxxxx_xxxxy, to_xxxxx_xxxyy, to_xxxxx_xxxyz, to_xxxxx_xxy, to_xxxxx_xxyyy, to_xxxxx_xxyyz, to_xxxxx_xxyzz, to_xxxxx_xxz, to_xxxxx_xyy, to_xxxxx_xyyyy, to_xxxxx_xyyyz, to_xxxxx_xyyzz, to_xxxxx_xyz, to_xxxxx_xyzzz, to_xxxxx_xzz, to_xxxxx_yyy, to_xxxxx_yyyyy, to_xxxxx_yyyyz, to_xxxxx_yyyzz, to_xxxxx_yyz, to_xxxxx_yyzzz, to_xxxxx_yzz, to_xxxxx_yzzzz, to_xxxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxx_xxxx[k] = -8.0 * to_xxx_xxxxy[k] * tke_0 + 4.0 * to_xxxxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxxy[k] = 4.0 * to_xxx_xxx[k] - 8.0 * to_xxx_xxxyy[k] * tke_0 - 2.0 * to_xxxxx_xxx[k] * tbe_0 + 4.0 * to_xxxxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxxz[k] = -8.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxyy[k] = 8.0 * to_xxx_xxy[k] - 8.0 * to_xxx_xxyyy[k] * tke_0 - 4.0 * to_xxxxx_xxy[k] * tbe_0 + 4.0 * to_xxxxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxyz[k] = 4.0 * to_xxx_xxz[k] - 8.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxxx_xxz[k] * tbe_0 + 4.0 * to_xxxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxzz[k] = -8.0 * to_xxx_xxyzz[k] * tke_0 + 4.0 * to_xxxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xyyy[k] = 12.0 * to_xxx_xyy[k] - 8.0 * to_xxx_xyyyy[k] * tke_0 - 6.0 * to_xxxxx_xyy[k] * tbe_0 + 4.0 * to_xxxxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xyyz[k] = 8.0 * to_xxx_xyz[k] - 8.0 * to_xxx_xyyyz[k] * tke_0 - 4.0 * to_xxxxx_xyz[k] * tbe_0 + 4.0 * to_xxxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xyzz[k] = 4.0 * to_xxx_xzz[k] - 8.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxxx_xzz[k] * tbe_0 + 4.0 * to_xxxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xzzz[k] = -8.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yyyy[k] = 16.0 * to_xxx_yyy[k] - 8.0 * to_xxx_yyyyy[k] * tke_0 - 8.0 * to_xxxxx_yyy[k] * tbe_0 + 4.0 * to_xxxxx_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yyyz[k] = 12.0 * to_xxx_yyz[k] - 8.0 * to_xxx_yyyyz[k] * tke_0 - 6.0 * to_xxxxx_yyz[k] * tbe_0 + 4.0 * to_xxxxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yyzz[k] = 8.0 * to_xxx_yzz[k] - 8.0 * to_xxx_yyyzz[k] * tke_0 - 4.0 * to_xxxxx_yzz[k] * tbe_0 + 4.0 * to_xxxxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yzzz[k] = 4.0 * to_xxx_zzz[k] - 8.0 * to_xxx_yyzzz[k] * tke_0 - 2.0 * to_xxxxx_zzz[k] * tbe_0 + 4.0 * to_xxxxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_zzzz[k] = -8.0 * to_xxx_yzzzz[k] * tke_0 + 4.0 * to_xxxxx_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-255 components of targeted buffer : GG

        auto to_x_y_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 15);

        auto to_x_y_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 16);

        auto to_x_y_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 17);

        auto to_x_y_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 18);

        auto to_x_y_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 19);

        auto to_x_y_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 20);

        auto to_x_y_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 21);

        auto to_x_y_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 22);

        auto to_x_y_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 23);

        auto to_x_y_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 24);

        auto to_x_y_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 25);

        auto to_x_y_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 26);

        auto to_x_y_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 27);

        auto to_x_y_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 28);

        auto to_x_y_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_x_y_xxxy_xxxx, to_x_y_xxxy_xxxy, to_x_y_xxxy_xxxz, to_x_y_xxxy_xxyy, to_x_y_xxxy_xxyz, to_x_y_xxxy_xxzz, to_x_y_xxxy_xyyy, to_x_y_xxxy_xyyz, to_x_y_xxxy_xyzz, to_x_y_xxxy_xzzz, to_x_y_xxxy_yyyy, to_x_y_xxxy_yyyz, to_x_y_xxxy_yyzz, to_x_y_xxxy_yzzz, to_x_y_xxxy_zzzz, to_xxxxy_xxx, to_xxxxy_xxxxy, to_xxxxy_xxxyy, to_xxxxy_xxxyz, to_xxxxy_xxy, to_xxxxy_xxyyy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xyy, to_xxxxy_xyyyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_yyy, to_xxxxy_yyyyy, to_xxxxy_yyyyz, to_xxxxy_yyyzz, to_xxxxy_yyz, to_xxxxy_yyzzz, to_xxxxy_yzz, to_xxxxy_yzzzz, to_xxxxy_zzz, to_xxy_xxx, to_xxy_xxxxy, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_yyy, to_xxy_yyyyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxy_xxxx[k] = -6.0 * to_xxy_xxxxy[k] * tke_0 + 4.0 * to_xxxxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxxy[k] = 3.0 * to_xxy_xxx[k] - 6.0 * to_xxy_xxxyy[k] * tke_0 - 2.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxxz[k] = -6.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxyy[k] = 6.0 * to_xxy_xxy[k] - 6.0 * to_xxy_xxyyy[k] * tke_0 - 4.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxyz[k] = 3.0 * to_xxy_xxz[k] - 6.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxzz[k] = -6.0 * to_xxy_xxyzz[k] * tke_0 + 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xyyy[k] = 9.0 * to_xxy_xyy[k] - 6.0 * to_xxy_xyyyy[k] * tke_0 - 6.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xyyz[k] = 6.0 * to_xxy_xyz[k] - 6.0 * to_xxy_xyyyz[k] * tke_0 - 4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xyzz[k] = 3.0 * to_xxy_xzz[k] - 6.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xzzz[k] = -6.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yyyy[k] = 12.0 * to_xxy_yyy[k] - 6.0 * to_xxy_yyyyy[k] * tke_0 - 8.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yyyz[k] = 9.0 * to_xxy_yyz[k] - 6.0 * to_xxy_yyyyz[k] * tke_0 - 6.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yyzz[k] = 6.0 * to_xxy_yzz[k] - 6.0 * to_xxy_yyyzz[k] * tke_0 - 4.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yzzz[k] = 3.0 * to_xxy_zzz[k] - 6.0 * to_xxy_yyzzz[k] * tke_0 - 2.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_zzzz[k] = -6.0 * to_xxy_yzzzz[k] * tke_0 + 4.0 * to_xxxxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 255-270 components of targeted buffer : GG

        auto to_x_y_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 30);

        auto to_x_y_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 31);

        auto to_x_y_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 32);

        auto to_x_y_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 33);

        auto to_x_y_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 34);

        auto to_x_y_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 35);

        auto to_x_y_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 36);

        auto to_x_y_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 37);

        auto to_x_y_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 38);

        auto to_x_y_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 39);

        auto to_x_y_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 40);

        auto to_x_y_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 41);

        auto to_x_y_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 42);

        auto to_x_y_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 43);

        auto to_x_y_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_x_y_xxxz_xxxx, to_x_y_xxxz_xxxy, to_x_y_xxxz_xxxz, to_x_y_xxxz_xxyy, to_x_y_xxxz_xxyz, to_x_y_xxxz_xxzz, to_x_y_xxxz_xyyy, to_x_y_xxxz_xyyz, to_x_y_xxxz_xyzz, to_x_y_xxxz_xzzz, to_x_y_xxxz_yyyy, to_x_y_xxxz_yyyz, to_x_y_xxxz_yyzz, to_x_y_xxxz_yzzz, to_x_y_xxxz_zzzz, to_xxxxz_xxx, to_xxxxz_xxxxy, to_xxxxz_xxxyy, to_xxxxz_xxxyz, to_xxxxz_xxy, to_xxxxz_xxyyy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xyy, to_xxxxz_xyyyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_yyy, to_xxxxz_yyyyy, to_xxxxz_yyyyz, to_xxxxz_yyyzz, to_xxxxz_yyz, to_xxxxz_yyzzz, to_xxxxz_yzz, to_xxxxz_yzzzz, to_xxxxz_zzz, to_xxz_xxx, to_xxz_xxxxy, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_yyy, to_xxz_yyyyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxz_xxxx[k] = -6.0 * to_xxz_xxxxy[k] * tke_0 + 4.0 * to_xxxxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxxy[k] = 3.0 * to_xxz_xxx[k] - 6.0 * to_xxz_xxxyy[k] * tke_0 - 2.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxxz[k] = -6.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxyy[k] = 6.0 * to_xxz_xxy[k] - 6.0 * to_xxz_xxyyy[k] * tke_0 - 4.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxyz[k] = 3.0 * to_xxz_xxz[k] - 6.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxzz[k] = -6.0 * to_xxz_xxyzz[k] * tke_0 + 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xyyy[k] = 9.0 * to_xxz_xyy[k] - 6.0 * to_xxz_xyyyy[k] * tke_0 - 6.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xyyz[k] = 6.0 * to_xxz_xyz[k] - 6.0 * to_xxz_xyyyz[k] * tke_0 - 4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xyzz[k] = 3.0 * to_xxz_xzz[k] - 6.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xzzz[k] = -6.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yyyy[k] = 12.0 * to_xxz_yyy[k] - 6.0 * to_xxz_yyyyy[k] * tke_0 - 8.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yyyz[k] = 9.0 * to_xxz_yyz[k] - 6.0 * to_xxz_yyyyz[k] * tke_0 - 6.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yyzz[k] = 6.0 * to_xxz_yzz[k] - 6.0 * to_xxz_yyyzz[k] * tke_0 - 4.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yzzz[k] = 3.0 * to_xxz_zzz[k] - 6.0 * to_xxz_yyzzz[k] * tke_0 - 2.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_zzzz[k] = -6.0 * to_xxz_yzzzz[k] * tke_0 + 4.0 * to_xxxxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-285 components of targeted buffer : GG

        auto to_x_y_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 45);

        auto to_x_y_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 46);

        auto to_x_y_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 47);

        auto to_x_y_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 48);

        auto to_x_y_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 49);

        auto to_x_y_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 50);

        auto to_x_y_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 51);

        auto to_x_y_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 52);

        auto to_x_y_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 53);

        auto to_x_y_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 54);

        auto to_x_y_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 55);

        auto to_x_y_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 56);

        auto to_x_y_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 57);

        auto to_x_y_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 58);

        auto to_x_y_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_x_y_xxyy_xxxx, to_x_y_xxyy_xxxy, to_x_y_xxyy_xxxz, to_x_y_xxyy_xxyy, to_x_y_xxyy_xxyz, to_x_y_xxyy_xxzz, to_x_y_xxyy_xyyy, to_x_y_xxyy_xyyz, to_x_y_xxyy_xyzz, to_x_y_xxyy_xzzz, to_x_y_xxyy_yyyy, to_x_y_xxyy_yyyz, to_x_y_xxyy_yyzz, to_x_y_xxyy_yzzz, to_x_y_xxyy_zzzz, to_xxxyy_xxx, to_xxxyy_xxxxy, to_xxxyy_xxxyy, to_xxxyy_xxxyz, to_xxxyy_xxy, to_xxxyy_xxyyy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xyy, to_xxxyy_xyyyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_yyy, to_xxxyy_yyyyy, to_xxxyy_yyyyz, to_xxxyy_yyyzz, to_xxxyy_yyz, to_xxxyy_yyzzz, to_xxxyy_yzz, to_xxxyy_yzzzz, to_xxxyy_zzz, to_xyy_xxx, to_xyy_xxxxy, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_yyy, to_xyy_yyyyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyy_xxxx[k] = -4.0 * to_xyy_xxxxy[k] * tke_0 + 4.0 * to_xxxyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxxy[k] = 2.0 * to_xyy_xxx[k] - 4.0 * to_xyy_xxxyy[k] * tke_0 - 2.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxxz[k] = -4.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxyy[k] = 4.0 * to_xyy_xxy[k] - 4.0 * to_xyy_xxyyy[k] * tke_0 - 4.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxyz[k] = 2.0 * to_xyy_xxz[k] - 4.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxzz[k] = -4.0 * to_xyy_xxyzz[k] * tke_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xyyy[k] = 6.0 * to_xyy_xyy[k] - 4.0 * to_xyy_xyyyy[k] * tke_0 - 6.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xyyz[k] = 4.0 * to_xyy_xyz[k] - 4.0 * to_xyy_xyyyz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xyzz[k] = 2.0 * to_xyy_xzz[k] - 4.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xzzz[k] = -4.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yyyy[k] = 8.0 * to_xyy_yyy[k] - 4.0 * to_xyy_yyyyy[k] * tke_0 - 8.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yyyz[k] = 6.0 * to_xyy_yyz[k] - 4.0 * to_xyy_yyyyz[k] * tke_0 - 6.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yyzz[k] = 4.0 * to_xyy_yzz[k] - 4.0 * to_xyy_yyyzz[k] * tke_0 - 4.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yzzz[k] = 2.0 * to_xyy_zzz[k] - 4.0 * to_xyy_yyzzz[k] * tke_0 - 2.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_zzzz[k] = -4.0 * to_xyy_yzzzz[k] * tke_0 + 4.0 * to_xxxyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 285-300 components of targeted buffer : GG

        auto to_x_y_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 60);

        auto to_x_y_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 61);

        auto to_x_y_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 62);

        auto to_x_y_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 63);

        auto to_x_y_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 64);

        auto to_x_y_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 65);

        auto to_x_y_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 66);

        auto to_x_y_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 67);

        auto to_x_y_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 68);

        auto to_x_y_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 69);

        auto to_x_y_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 70);

        auto to_x_y_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 71);

        auto to_x_y_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 72);

        auto to_x_y_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 73);

        auto to_x_y_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_x_y_xxyz_xxxx, to_x_y_xxyz_xxxy, to_x_y_xxyz_xxxz, to_x_y_xxyz_xxyy, to_x_y_xxyz_xxyz, to_x_y_xxyz_xxzz, to_x_y_xxyz_xyyy, to_x_y_xxyz_xyyz, to_x_y_xxyz_xyzz, to_x_y_xxyz_xzzz, to_x_y_xxyz_yyyy, to_x_y_xxyz_yyyz, to_x_y_xxyz_yyzz, to_x_y_xxyz_yzzz, to_x_y_xxyz_zzzz, to_xxxyz_xxx, to_xxxyz_xxxxy, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_yyy, to_xxxyz_yyyyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_xyz_xxx, to_xyz_xxxxy, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_yyy, to_xyz_yyyyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyz_xxxx[k] = -4.0 * to_xyz_xxxxy[k] * tke_0 + 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxxy[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxxz[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxyy[k] = 4.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 4.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxyz[k] = 2.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxzz[k] = -4.0 * to_xyz_xxyzz[k] * tke_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xyyy[k] = 6.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyyy[k] * tke_0 - 6.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xyyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyyyz[k] * tke_0 - 4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xyzz[k] = 2.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yyyy[k] = 8.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyyy[k] * tke_0 - 8.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yyyz[k] = 6.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyyyz[k] * tke_0 - 6.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yyzz[k] = 4.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 4.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_zzzz[k] = -4.0 * to_xyz_yzzzz[k] * tke_0 + 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-315 components of targeted buffer : GG

        auto to_x_y_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 75);

        auto to_x_y_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 76);

        auto to_x_y_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 77);

        auto to_x_y_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 78);

        auto to_x_y_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 79);

        auto to_x_y_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 80);

        auto to_x_y_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 81);

        auto to_x_y_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 82);

        auto to_x_y_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 83);

        auto to_x_y_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 84);

        auto to_x_y_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 85);

        auto to_x_y_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 86);

        auto to_x_y_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 87);

        auto to_x_y_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 88);

        auto to_x_y_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_x_y_xxzz_xxxx, to_x_y_xxzz_xxxy, to_x_y_xxzz_xxxz, to_x_y_xxzz_xxyy, to_x_y_xxzz_xxyz, to_x_y_xxzz_xxzz, to_x_y_xxzz_xyyy, to_x_y_xxzz_xyyz, to_x_y_xxzz_xyzz, to_x_y_xxzz_xzzz, to_x_y_xxzz_yyyy, to_x_y_xxzz_yyyz, to_x_y_xxzz_yyzz, to_x_y_xxzz_yzzz, to_x_y_xxzz_zzzz, to_xxxzz_xxx, to_xxxzz_xxxxy, to_xxxzz_xxxyy, to_xxxzz_xxxyz, to_xxxzz_xxy, to_xxxzz_xxyyy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xyy, to_xxxzz_xyyyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_yyy, to_xxxzz_yyyyy, to_xxxzz_yyyyz, to_xxxzz_yyyzz, to_xxxzz_yyz, to_xxxzz_yyzzz, to_xxxzz_yzz, to_xxxzz_yzzzz, to_xxxzz_zzz, to_xzz_xxx, to_xzz_xxxxy, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_yyy, to_xzz_yyyyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxzz_xxxx[k] = -4.0 * to_xzz_xxxxy[k] * tke_0 + 4.0 * to_xxxzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxxy[k] = 2.0 * to_xzz_xxx[k] - 4.0 * to_xzz_xxxyy[k] * tke_0 - 2.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxxz[k] = -4.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxyy[k] = 4.0 * to_xzz_xxy[k] - 4.0 * to_xzz_xxyyy[k] * tke_0 - 4.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxyz[k] = 2.0 * to_xzz_xxz[k] - 4.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxzz[k] = -4.0 * to_xzz_xxyzz[k] * tke_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xyyy[k] = 6.0 * to_xzz_xyy[k] - 4.0 * to_xzz_xyyyy[k] * tke_0 - 6.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xyyz[k] = 4.0 * to_xzz_xyz[k] - 4.0 * to_xzz_xyyyz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xyzz[k] = 2.0 * to_xzz_xzz[k] - 4.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xzzz[k] = -4.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yyyy[k] = 8.0 * to_xzz_yyy[k] - 4.0 * to_xzz_yyyyy[k] * tke_0 - 8.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yyyz[k] = 6.0 * to_xzz_yyz[k] - 4.0 * to_xzz_yyyyz[k] * tke_0 - 6.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yyzz[k] = 4.0 * to_xzz_yzz[k] - 4.0 * to_xzz_yyyzz[k] * tke_0 - 4.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yzzz[k] = 2.0 * to_xzz_zzz[k] - 4.0 * to_xzz_yyzzz[k] * tke_0 - 2.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_zzzz[k] = -4.0 * to_xzz_yzzzz[k] * tke_0 + 4.0 * to_xxxzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 315-330 components of targeted buffer : GG

        auto to_x_y_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 90);

        auto to_x_y_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 91);

        auto to_x_y_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 92);

        auto to_x_y_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 93);

        auto to_x_y_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 94);

        auto to_x_y_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 95);

        auto to_x_y_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 96);

        auto to_x_y_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 97);

        auto to_x_y_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 98);

        auto to_x_y_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 99);

        auto to_x_y_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 100);

        auto to_x_y_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 101);

        auto to_x_y_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 102);

        auto to_x_y_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 103);

        auto to_x_y_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_x_y_xyyy_xxxx, to_x_y_xyyy_xxxy, to_x_y_xyyy_xxxz, to_x_y_xyyy_xxyy, to_x_y_xyyy_xxyz, to_x_y_xyyy_xxzz, to_x_y_xyyy_xyyy, to_x_y_xyyy_xyyz, to_x_y_xyyy_xyzz, to_x_y_xyyy_xzzz, to_x_y_xyyy_yyyy, to_x_y_xyyy_yyyz, to_x_y_xyyy_yyzz, to_x_y_xyyy_yzzz, to_x_y_xyyy_zzzz, to_xxyyy_xxx, to_xxyyy_xxxxy, to_xxyyy_xxxyy, to_xxyyy_xxxyz, to_xxyyy_xxy, to_xxyyy_xxyyy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xyy, to_xxyyy_xyyyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_yyy, to_xxyyy_yyyyy, to_xxyyy_yyyyz, to_xxyyy_yyyzz, to_xxyyy_yyz, to_xxyyy_yyzzz, to_xxyyy_yzz, to_xxyyy_yzzzz, to_xxyyy_zzz, to_yyy_xxx, to_yyy_xxxxy, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_yyy, to_yyy_yyyyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyy_xxxx[k] = -2.0 * to_yyy_xxxxy[k] * tke_0 + 4.0 * to_xxyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxxy[k] = to_yyy_xxx[k] - 2.0 * to_yyy_xxxyy[k] * tke_0 - 2.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxxz[k] = -2.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxyy[k] = 2.0 * to_yyy_xxy[k] - 2.0 * to_yyy_xxyyy[k] * tke_0 - 4.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxyz[k] = to_yyy_xxz[k] - 2.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxzz[k] = -2.0 * to_yyy_xxyzz[k] * tke_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xyyy[k] = 3.0 * to_yyy_xyy[k] - 2.0 * to_yyy_xyyyy[k] * tke_0 - 6.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xyyz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xyyyz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xyzz[k] = to_yyy_xzz[k] - 2.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xzzz[k] = -2.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yyyy[k] = 4.0 * to_yyy_yyy[k] - 2.0 * to_yyy_yyyyy[k] * tke_0 - 8.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yyyz[k] = 3.0 * to_yyy_yyz[k] - 2.0 * to_yyy_yyyyz[k] * tke_0 - 6.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yyzz[k] = 2.0 * to_yyy_yzz[k] - 2.0 * to_yyy_yyyzz[k] * tke_0 - 4.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yzzz[k] = to_yyy_zzz[k] - 2.0 * to_yyy_yyzzz[k] * tke_0 - 2.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_zzzz[k] = -2.0 * to_yyy_yzzzz[k] * tke_0 + 4.0 * to_xxyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-345 components of targeted buffer : GG

        auto to_x_y_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 105);

        auto to_x_y_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 106);

        auto to_x_y_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 107);

        auto to_x_y_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 108);

        auto to_x_y_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 109);

        auto to_x_y_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 110);

        auto to_x_y_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 111);

        auto to_x_y_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 112);

        auto to_x_y_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 113);

        auto to_x_y_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 114);

        auto to_x_y_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 115);

        auto to_x_y_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 116);

        auto to_x_y_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 117);

        auto to_x_y_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 118);

        auto to_x_y_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_x_y_xyyz_xxxx, to_x_y_xyyz_xxxy, to_x_y_xyyz_xxxz, to_x_y_xyyz_xxyy, to_x_y_xyyz_xxyz, to_x_y_xyyz_xxzz, to_x_y_xyyz_xyyy, to_x_y_xyyz_xyyz, to_x_y_xyyz_xyzz, to_x_y_xyyz_xzzz, to_x_y_xyyz_yyyy, to_x_y_xyyz_yyyz, to_x_y_xyyz_yyzz, to_x_y_xyyz_yzzz, to_x_y_xyyz_zzzz, to_xxyyz_xxx, to_xxyyz_xxxxy, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_yyy, to_xxyyz_yyyyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_yyz_xxx, to_yyz_xxxxy, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_yyy, to_yyz_yyyyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyz_xxxx[k] = -2.0 * to_yyz_xxxxy[k] * tke_0 + 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxxy[k] = to_yyz_xxx[k] - 2.0 * to_yyz_xxxyy[k] * tke_0 - 2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxxz[k] = -2.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxyy[k] = 2.0 * to_yyz_xxy[k] - 2.0 * to_yyz_xxyyy[k] * tke_0 - 4.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxyz[k] = to_yyz_xxz[k] - 2.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxzz[k] = -2.0 * to_yyz_xxyzz[k] * tke_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xyyy[k] = 3.0 * to_yyz_xyy[k] - 2.0 * to_yyz_xyyyy[k] * tke_0 - 6.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xyyz[k] = 2.0 * to_yyz_xyz[k] - 2.0 * to_yyz_xyyyz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xyzz[k] = to_yyz_xzz[k] - 2.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xzzz[k] = -2.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yyyy[k] = 4.0 * to_yyz_yyy[k] - 2.0 * to_yyz_yyyyy[k] * tke_0 - 8.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yyyz[k] = 3.0 * to_yyz_yyz[k] - 2.0 * to_yyz_yyyyz[k] * tke_0 - 6.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yyzz[k] = 2.0 * to_yyz_yzz[k] - 2.0 * to_yyz_yyyzz[k] * tke_0 - 4.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yzzz[k] = to_yyz_zzz[k] - 2.0 * to_yyz_yyzzz[k] * tke_0 - 2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_zzzz[k] = -2.0 * to_yyz_yzzzz[k] * tke_0 + 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 345-360 components of targeted buffer : GG

        auto to_x_y_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 120);

        auto to_x_y_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 121);

        auto to_x_y_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 122);

        auto to_x_y_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 123);

        auto to_x_y_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 124);

        auto to_x_y_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 125);

        auto to_x_y_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 126);

        auto to_x_y_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 127);

        auto to_x_y_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 128);

        auto to_x_y_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 129);

        auto to_x_y_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 130);

        auto to_x_y_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 131);

        auto to_x_y_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 132);

        auto to_x_y_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 133);

        auto to_x_y_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_x_y_xyzz_xxxx, to_x_y_xyzz_xxxy, to_x_y_xyzz_xxxz, to_x_y_xyzz_xxyy, to_x_y_xyzz_xxyz, to_x_y_xyzz_xxzz, to_x_y_xyzz_xyyy, to_x_y_xyzz_xyyz, to_x_y_xyzz_xyzz, to_x_y_xyzz_xzzz, to_x_y_xyzz_yyyy, to_x_y_xyzz_yyyz, to_x_y_xyzz_yyzz, to_x_y_xyzz_yzzz, to_x_y_xyzz_zzzz, to_xxyzz_xxx, to_xxyzz_xxxxy, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_yyy, to_xxyzz_yyyyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_yzz_xxx, to_yzz_xxxxy, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_yyy, to_yzz_yyyyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyzz_xxxx[k] = -2.0 * to_yzz_xxxxy[k] * tke_0 + 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxxy[k] = to_yzz_xxx[k] - 2.0 * to_yzz_xxxyy[k] * tke_0 - 2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxxz[k] = -2.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxyy[k] = 2.0 * to_yzz_xxy[k] - 2.0 * to_yzz_xxyyy[k] * tke_0 - 4.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxyz[k] = to_yzz_xxz[k] - 2.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxzz[k] = -2.0 * to_yzz_xxyzz[k] * tke_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xyyy[k] = 3.0 * to_yzz_xyy[k] - 2.0 * to_yzz_xyyyy[k] * tke_0 - 6.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xyyz[k] = 2.0 * to_yzz_xyz[k] - 2.0 * to_yzz_xyyyz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xyzz[k] = to_yzz_xzz[k] - 2.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xzzz[k] = -2.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yyyy[k] = 4.0 * to_yzz_yyy[k] - 2.0 * to_yzz_yyyyy[k] * tke_0 - 8.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yyyz[k] = 3.0 * to_yzz_yyz[k] - 2.0 * to_yzz_yyyyz[k] * tke_0 - 6.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yyzz[k] = 2.0 * to_yzz_yzz[k] - 2.0 * to_yzz_yyyzz[k] * tke_0 - 4.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yzzz[k] = to_yzz_zzz[k] - 2.0 * to_yzz_yyzzz[k] * tke_0 - 2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_zzzz[k] = -2.0 * to_yzz_yzzzz[k] * tke_0 + 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-375 components of targeted buffer : GG

        auto to_x_y_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 135);

        auto to_x_y_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 136);

        auto to_x_y_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 137);

        auto to_x_y_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 138);

        auto to_x_y_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 139);

        auto to_x_y_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 140);

        auto to_x_y_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 141);

        auto to_x_y_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 142);

        auto to_x_y_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 143);

        auto to_x_y_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 144);

        auto to_x_y_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 145);

        auto to_x_y_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 146);

        auto to_x_y_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 147);

        auto to_x_y_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 148);

        auto to_x_y_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_x_y_xzzz_xxxx, to_x_y_xzzz_xxxy, to_x_y_xzzz_xxxz, to_x_y_xzzz_xxyy, to_x_y_xzzz_xxyz, to_x_y_xzzz_xxzz, to_x_y_xzzz_xyyy, to_x_y_xzzz_xyyz, to_x_y_xzzz_xyzz, to_x_y_xzzz_xzzz, to_x_y_xzzz_yyyy, to_x_y_xzzz_yyyz, to_x_y_xzzz_yyzz, to_x_y_xzzz_yzzz, to_x_y_xzzz_zzzz, to_xxzzz_xxx, to_xxzzz_xxxxy, to_xxzzz_xxxyy, to_xxzzz_xxxyz, to_xxzzz_xxy, to_xxzzz_xxyyy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xyy, to_xxzzz_xyyyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_yyy, to_xxzzz_yyyyy, to_xxzzz_yyyyz, to_xxzzz_yyyzz, to_xxzzz_yyz, to_xxzzz_yyzzz, to_xxzzz_yzz, to_xxzzz_yzzzz, to_xxzzz_zzz, to_zzz_xxx, to_zzz_xxxxy, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_yyy, to_zzz_yyyyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzzz_xxxx[k] = -2.0 * to_zzz_xxxxy[k] * tke_0 + 4.0 * to_xxzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxxy[k] = to_zzz_xxx[k] - 2.0 * to_zzz_xxxyy[k] * tke_0 - 2.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxxz[k] = -2.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxyy[k] = 2.0 * to_zzz_xxy[k] - 2.0 * to_zzz_xxyyy[k] * tke_0 - 4.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxyz[k] = to_zzz_xxz[k] - 2.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxzz[k] = -2.0 * to_zzz_xxyzz[k] * tke_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xyyy[k] = 3.0 * to_zzz_xyy[k] - 2.0 * to_zzz_xyyyy[k] * tke_0 - 6.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xyyz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xyyyz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xyzz[k] = to_zzz_xzz[k] - 2.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xzzz[k] = -2.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yyyy[k] = 4.0 * to_zzz_yyy[k] - 2.0 * to_zzz_yyyyy[k] * tke_0 - 8.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yyyz[k] = 3.0 * to_zzz_yyz[k] - 2.0 * to_zzz_yyyyz[k] * tke_0 - 6.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yyzz[k] = 2.0 * to_zzz_yzz[k] - 2.0 * to_zzz_yyyzz[k] * tke_0 - 4.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yzzz[k] = to_zzz_zzz[k] - 2.0 * to_zzz_yyzzz[k] * tke_0 - 2.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_zzzz[k] = -2.0 * to_zzz_yzzzz[k] * tke_0 + 4.0 * to_xxzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 375-390 components of targeted buffer : GG

        auto to_x_y_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 150);

        auto to_x_y_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 151);

        auto to_x_y_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 152);

        auto to_x_y_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 153);

        auto to_x_y_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 154);

        auto to_x_y_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 155);

        auto to_x_y_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 156);

        auto to_x_y_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 157);

        auto to_x_y_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 158);

        auto to_x_y_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 159);

        auto to_x_y_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 160);

        auto to_x_y_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 161);

        auto to_x_y_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 162);

        auto to_x_y_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 163);

        auto to_x_y_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_x_y_yyyy_xxxx, to_x_y_yyyy_xxxy, to_x_y_yyyy_xxxz, to_x_y_yyyy_xxyy, to_x_y_yyyy_xxyz, to_x_y_yyyy_xxzz, to_x_y_yyyy_xyyy, to_x_y_yyyy_xyyz, to_x_y_yyyy_xyzz, to_x_y_yyyy_xzzz, to_x_y_yyyy_yyyy, to_x_y_yyyy_yyyz, to_x_y_yyyy_yyzz, to_x_y_yyyy_yzzz, to_x_y_yyyy_zzzz, to_xyyyy_xxx, to_xyyyy_xxxxy, to_xyyyy_xxxyy, to_xyyyy_xxxyz, to_xyyyy_xxy, to_xyyyy_xxyyy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xyy, to_xyyyy_xyyyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_yyy, to_xyyyy_yyyyy, to_xyyyy_yyyyz, to_xyyyy_yyyzz, to_xyyyy_yyz, to_xyyyy_yyzzz, to_xyyyy_yzz, to_xyyyy_yzzzz, to_xyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyy_xxxx[k] = 4.0 * to_xyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxxy[k] = -2.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxxz[k] = 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxyy[k] = -4.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxyz[k] = -2.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxzz[k] = 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xyyy[k] = -6.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xyyz[k] = -4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xyzz[k] = -2.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xzzz[k] = 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yyyy[k] = -8.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yyyz[k] = -6.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yyzz[k] = -4.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yzzz[k] = -2.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_zzzz[k] = 4.0 * to_xyyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-405 components of targeted buffer : GG

        auto to_x_y_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 165);

        auto to_x_y_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 166);

        auto to_x_y_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 167);

        auto to_x_y_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 168);

        auto to_x_y_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 169);

        auto to_x_y_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 170);

        auto to_x_y_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 171);

        auto to_x_y_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 172);

        auto to_x_y_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 173);

        auto to_x_y_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 174);

        auto to_x_y_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 175);

        auto to_x_y_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 176);

        auto to_x_y_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 177);

        auto to_x_y_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 178);

        auto to_x_y_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_x_y_yyyz_xxxx, to_x_y_yyyz_xxxy, to_x_y_yyyz_xxxz, to_x_y_yyyz_xxyy, to_x_y_yyyz_xxyz, to_x_y_yyyz_xxzz, to_x_y_yyyz_xyyy, to_x_y_yyyz_xyyz, to_x_y_yyyz_xyzz, to_x_y_yyyz_xzzz, to_x_y_yyyz_yyyy, to_x_y_yyyz_yyyz, to_x_y_yyyz_yyzz, to_x_y_yyyz_yzzz, to_x_y_yyyz_zzzz, to_xyyyz_xxx, to_xyyyz_xxxxy, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_yyy, to_xyyyz_yyyyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyz_xxxx[k] = 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxxy[k] = -2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxxz[k] = 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxyy[k] = -4.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxyz[k] = -2.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxzz[k] = 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xyyy[k] = -6.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xyyz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xyzz[k] = -2.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xzzz[k] = 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yyyy[k] = -8.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yyyz[k] = -6.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yyzz[k] = -4.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yzzz[k] = -2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_zzzz[k] = 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 405-420 components of targeted buffer : GG

        auto to_x_y_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 180);

        auto to_x_y_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 181);

        auto to_x_y_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 182);

        auto to_x_y_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 183);

        auto to_x_y_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 184);

        auto to_x_y_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 185);

        auto to_x_y_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 186);

        auto to_x_y_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 187);

        auto to_x_y_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 188);

        auto to_x_y_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 189);

        auto to_x_y_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 190);

        auto to_x_y_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 191);

        auto to_x_y_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 192);

        auto to_x_y_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 193);

        auto to_x_y_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_x_y_yyzz_xxxx, to_x_y_yyzz_xxxy, to_x_y_yyzz_xxxz, to_x_y_yyzz_xxyy, to_x_y_yyzz_xxyz, to_x_y_yyzz_xxzz, to_x_y_yyzz_xyyy, to_x_y_yyzz_xyyz, to_x_y_yyzz_xyzz, to_x_y_yyzz_xzzz, to_x_y_yyzz_yyyy, to_x_y_yyzz_yyyz, to_x_y_yyzz_yyzz, to_x_y_yyzz_yzzz, to_x_y_yyzz_zzzz, to_xyyzz_xxx, to_xyyzz_xxxxy, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_yyy, to_xyyzz_yyyyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyzz_xxxx[k] = 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxxy[k] = -2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxxz[k] = 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxyy[k] = -4.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxyz[k] = -2.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxzz[k] = 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xyyy[k] = -6.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xyyz[k] = -4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xyzz[k] = -2.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xzzz[k] = 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yyyy[k] = -8.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yyyz[k] = -6.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yyzz[k] = -4.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yzzz[k] = -2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_zzzz[k] = 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-435 components of targeted buffer : GG

        auto to_x_y_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 195);

        auto to_x_y_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 196);

        auto to_x_y_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 197);

        auto to_x_y_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 198);

        auto to_x_y_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 199);

        auto to_x_y_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 200);

        auto to_x_y_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 201);

        auto to_x_y_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 202);

        auto to_x_y_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 203);

        auto to_x_y_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 204);

        auto to_x_y_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 205);

        auto to_x_y_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 206);

        auto to_x_y_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 207);

        auto to_x_y_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 208);

        auto to_x_y_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_x_y_yzzz_xxxx, to_x_y_yzzz_xxxy, to_x_y_yzzz_xxxz, to_x_y_yzzz_xxyy, to_x_y_yzzz_xxyz, to_x_y_yzzz_xxzz, to_x_y_yzzz_xyyy, to_x_y_yzzz_xyyz, to_x_y_yzzz_xyzz, to_x_y_yzzz_xzzz, to_x_y_yzzz_yyyy, to_x_y_yzzz_yyyz, to_x_y_yzzz_yyzz, to_x_y_yzzz_yzzz, to_x_y_yzzz_zzzz, to_xyzzz_xxx, to_xyzzz_xxxxy, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_yyy, to_xyzzz_yyyyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzzz_xxxx[k] = 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxxy[k] = -2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxxz[k] = 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxyy[k] = -4.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxyz[k] = -2.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxzz[k] = 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xyyy[k] = -6.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xyyz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xyzz[k] = -2.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xzzz[k] = 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yyyy[k] = -8.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yyyz[k] = -6.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yyzz[k] = -4.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yzzz[k] = -2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_zzzz[k] = 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 435-450 components of targeted buffer : GG

        auto to_x_y_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 210);

        auto to_x_y_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 211);

        auto to_x_y_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 212);

        auto to_x_y_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 213);

        auto to_x_y_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 214);

        auto to_x_y_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 215);

        auto to_x_y_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 216);

        auto to_x_y_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 217);

        auto to_x_y_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 218);

        auto to_x_y_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 219);

        auto to_x_y_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 220);

        auto to_x_y_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 221);

        auto to_x_y_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 222);

        auto to_x_y_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 223);

        auto to_x_y_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 1 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_x_y_zzzz_xxxx, to_x_y_zzzz_xxxy, to_x_y_zzzz_xxxz, to_x_y_zzzz_xxyy, to_x_y_zzzz_xxyz, to_x_y_zzzz_xxzz, to_x_y_zzzz_xyyy, to_x_y_zzzz_xyyz, to_x_y_zzzz_xyzz, to_x_y_zzzz_xzzz, to_x_y_zzzz_yyyy, to_x_y_zzzz_yyyz, to_x_y_zzzz_yyzz, to_x_y_zzzz_yzzz, to_x_y_zzzz_zzzz, to_xzzzz_xxx, to_xzzzz_xxxxy, to_xzzzz_xxxyy, to_xzzzz_xxxyz, to_xzzzz_xxy, to_xzzzz_xxyyy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xyy, to_xzzzz_xyyyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_yyy, to_xzzzz_yyyyy, to_xzzzz_yyyyz, to_xzzzz_yyyzz, to_xzzzz_yyz, to_xzzzz_yyzzz, to_xzzzz_yzz, to_xzzzz_yzzzz, to_xzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzzz_xxxx[k] = 4.0 * to_xzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxxy[k] = -2.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxxz[k] = 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxyy[k] = -4.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxyz[k] = -2.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxzz[k] = 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xyyy[k] = -6.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xyyz[k] = -4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xyzz[k] = -2.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xzzz[k] = 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yyyy[k] = -8.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yyyz[k] = -6.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yyzz[k] = -4.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yzzz[k] = -2.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_zzzz[k] = 4.0 * to_xzzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-465 components of targeted buffer : GG

        auto to_x_z_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 0);

        auto to_x_z_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 1);

        auto to_x_z_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 2);

        auto to_x_z_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 3);

        auto to_x_z_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 4);

        auto to_x_z_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 5);

        auto to_x_z_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 6);

        auto to_x_z_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 7);

        auto to_x_z_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 8);

        auto to_x_z_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 9);

        auto to_x_z_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 10);

        auto to_x_z_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 11);

        auto to_x_z_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 12);

        auto to_x_z_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 13);

        auto to_x_z_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_x_z_xxxx_xxxx, to_x_z_xxxx_xxxy, to_x_z_xxxx_xxxz, to_x_z_xxxx_xxyy, to_x_z_xxxx_xxyz, to_x_z_xxxx_xxzz, to_x_z_xxxx_xyyy, to_x_z_xxxx_xyyz, to_x_z_xxxx_xyzz, to_x_z_xxxx_xzzz, to_x_z_xxxx_yyyy, to_x_z_xxxx_yyyz, to_x_z_xxxx_yyzz, to_x_z_xxxx_yzzz, to_x_z_xxxx_zzzz, to_xxx_xxx, to_xxx_xxxxz, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxx_zzzzz, to_xxxxx_xxx, to_xxxxx_xxxxz, to_xxxxx_xxxyz, to_xxxxx_xxxzz, to_xxxxx_xxy, to_xxxxx_xxyyz, to_xxxxx_xxyzz, to_xxxxx_xxz, to_xxxxx_xxzzz, to_xxxxx_xyy, to_xxxxx_xyyyz, to_xxxxx_xyyzz, to_xxxxx_xyz, to_xxxxx_xyzzz, to_xxxxx_xzz, to_xxxxx_xzzzz, to_xxxxx_yyy, to_xxxxx_yyyyz, to_xxxxx_yyyzz, to_xxxxx_yyz, to_xxxxx_yyzzz, to_xxxxx_yzz, to_xxxxx_yzzzz, to_xxxxx_zzz, to_xxxxx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxx_xxxx[k] = -8.0 * to_xxx_xxxxz[k] * tke_0 + 4.0 * to_xxxxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxxy[k] = -8.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxxz[k] = 4.0 * to_xxx_xxx[k] - 8.0 * to_xxx_xxxzz[k] * tke_0 - 2.0 * to_xxxxx_xxx[k] * tbe_0 + 4.0 * to_xxxxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxyy[k] = -8.0 * to_xxx_xxyyz[k] * tke_0 + 4.0 * to_xxxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxyz[k] = 4.0 * to_xxx_xxy[k] - 8.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxxx_xxy[k] * tbe_0 + 4.0 * to_xxxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxzz[k] = 8.0 * to_xxx_xxz[k] - 8.0 * to_xxx_xxzzz[k] * tke_0 - 4.0 * to_xxxxx_xxz[k] * tbe_0 + 4.0 * to_xxxxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xyyy[k] = -8.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xyyz[k] = 4.0 * to_xxx_xyy[k] - 8.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxxx_xyy[k] * tbe_0 + 4.0 * to_xxxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xyzz[k] = 8.0 * to_xxx_xyz[k] - 8.0 * to_xxx_xyzzz[k] * tke_0 - 4.0 * to_xxxxx_xyz[k] * tbe_0 + 4.0 * to_xxxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xzzz[k] = 12.0 * to_xxx_xzz[k] - 8.0 * to_xxx_xzzzz[k] * tke_0 - 6.0 * to_xxxxx_xzz[k] * tbe_0 + 4.0 * to_xxxxx_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yyyy[k] = -8.0 * to_xxx_yyyyz[k] * tke_0 + 4.0 * to_xxxxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yyyz[k] = 4.0 * to_xxx_yyy[k] - 8.0 * to_xxx_yyyzz[k] * tke_0 - 2.0 * to_xxxxx_yyy[k] * tbe_0 + 4.0 * to_xxxxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yyzz[k] = 8.0 * to_xxx_yyz[k] - 8.0 * to_xxx_yyzzz[k] * tke_0 - 4.0 * to_xxxxx_yyz[k] * tbe_0 + 4.0 * to_xxxxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yzzz[k] = 12.0 * to_xxx_yzz[k] - 8.0 * to_xxx_yzzzz[k] * tke_0 - 6.0 * to_xxxxx_yzz[k] * tbe_0 + 4.0 * to_xxxxx_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_zzzz[k] = 16.0 * to_xxx_zzz[k] - 8.0 * to_xxx_zzzzz[k] * tke_0 - 8.0 * to_xxxxx_zzz[k] * tbe_0 + 4.0 * to_xxxxx_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 465-480 components of targeted buffer : GG

        auto to_x_z_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 15);

        auto to_x_z_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 16);

        auto to_x_z_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 17);

        auto to_x_z_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 18);

        auto to_x_z_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 19);

        auto to_x_z_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 20);

        auto to_x_z_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 21);

        auto to_x_z_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 22);

        auto to_x_z_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 23);

        auto to_x_z_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 24);

        auto to_x_z_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 25);

        auto to_x_z_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 26);

        auto to_x_z_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 27);

        auto to_x_z_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 28);

        auto to_x_z_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_x_z_xxxy_xxxx, to_x_z_xxxy_xxxy, to_x_z_xxxy_xxxz, to_x_z_xxxy_xxyy, to_x_z_xxxy_xxyz, to_x_z_xxxy_xxzz, to_x_z_xxxy_xyyy, to_x_z_xxxy_xyyz, to_x_z_xxxy_xyzz, to_x_z_xxxy_xzzz, to_x_z_xxxy_yyyy, to_x_z_xxxy_yyyz, to_x_z_xxxy_yyzz, to_x_z_xxxy_yzzz, to_x_z_xxxy_zzzz, to_xxxxy_xxx, to_xxxxy_xxxxz, to_xxxxy_xxxyz, to_xxxxy_xxxzz, to_xxxxy_xxy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xxzzz, to_xxxxy_xyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_xzzzz, to_xxxxy_yyy, to_xxxxy_yyyyz, to_xxxxy_yyyzz, to_xxxxy_yyz, to_xxxxy_yyzzz, to_xxxxy_yzz, to_xxxxy_yzzzz, to_xxxxy_zzz, to_xxxxy_zzzzz, to_xxy_xxx, to_xxy_xxxxz, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxy_xxxx[k] = -6.0 * to_xxy_xxxxz[k] * tke_0 + 4.0 * to_xxxxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxxy[k] = -6.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxxz[k] = 3.0 * to_xxy_xxx[k] - 6.0 * to_xxy_xxxzz[k] * tke_0 - 2.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxyy[k] = -6.0 * to_xxy_xxyyz[k] * tke_0 + 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxyz[k] = 3.0 * to_xxy_xxy[k] - 6.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxzz[k] = 6.0 * to_xxy_xxz[k] - 6.0 * to_xxy_xxzzz[k] * tke_0 - 4.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xyyy[k] = -6.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xyyz[k] = 3.0 * to_xxy_xyy[k] - 6.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xyzz[k] = 6.0 * to_xxy_xyz[k] - 6.0 * to_xxy_xyzzz[k] * tke_0 - 4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xzzz[k] = 9.0 * to_xxy_xzz[k] - 6.0 * to_xxy_xzzzz[k] * tke_0 - 6.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yyyy[k] = -6.0 * to_xxy_yyyyz[k] * tke_0 + 4.0 * to_xxxxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yyyz[k] = 3.0 * to_xxy_yyy[k] - 6.0 * to_xxy_yyyzz[k] * tke_0 - 2.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yyzz[k] = 6.0 * to_xxy_yyz[k] - 6.0 * to_xxy_yyzzz[k] * tke_0 - 4.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yzzz[k] = 9.0 * to_xxy_yzz[k] - 6.0 * to_xxy_yzzzz[k] * tke_0 - 6.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_zzzz[k] = 12.0 * to_xxy_zzz[k] - 6.0 * to_xxy_zzzzz[k] * tke_0 - 8.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-495 components of targeted buffer : GG

        auto to_x_z_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 30);

        auto to_x_z_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 31);

        auto to_x_z_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 32);

        auto to_x_z_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 33);

        auto to_x_z_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 34);

        auto to_x_z_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 35);

        auto to_x_z_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 36);

        auto to_x_z_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 37);

        auto to_x_z_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 38);

        auto to_x_z_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 39);

        auto to_x_z_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 40);

        auto to_x_z_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 41);

        auto to_x_z_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 42);

        auto to_x_z_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 43);

        auto to_x_z_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_x_z_xxxz_xxxx, to_x_z_xxxz_xxxy, to_x_z_xxxz_xxxz, to_x_z_xxxz_xxyy, to_x_z_xxxz_xxyz, to_x_z_xxxz_xxzz, to_x_z_xxxz_xyyy, to_x_z_xxxz_xyyz, to_x_z_xxxz_xyzz, to_x_z_xxxz_xzzz, to_x_z_xxxz_yyyy, to_x_z_xxxz_yyyz, to_x_z_xxxz_yyzz, to_x_z_xxxz_yzzz, to_x_z_xxxz_zzzz, to_xxxxz_xxx, to_xxxxz_xxxxz, to_xxxxz_xxxyz, to_xxxxz_xxxzz, to_xxxxz_xxy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xxzzz, to_xxxxz_xyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_xzzzz, to_xxxxz_yyy, to_xxxxz_yyyyz, to_xxxxz_yyyzz, to_xxxxz_yyz, to_xxxxz_yyzzz, to_xxxxz_yzz, to_xxxxz_yzzzz, to_xxxxz_zzz, to_xxxxz_zzzzz, to_xxz_xxx, to_xxz_xxxxz, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_xxz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxz_xxxx[k] = -6.0 * to_xxz_xxxxz[k] * tke_0 + 4.0 * to_xxxxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxxy[k] = -6.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxxz[k] = 3.0 * to_xxz_xxx[k] - 6.0 * to_xxz_xxxzz[k] * tke_0 - 2.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxyy[k] = -6.0 * to_xxz_xxyyz[k] * tke_0 + 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxyz[k] = 3.0 * to_xxz_xxy[k] - 6.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxzz[k] = 6.0 * to_xxz_xxz[k] - 6.0 * to_xxz_xxzzz[k] * tke_0 - 4.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xyyy[k] = -6.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xyyz[k] = 3.0 * to_xxz_xyy[k] - 6.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xyzz[k] = 6.0 * to_xxz_xyz[k] - 6.0 * to_xxz_xyzzz[k] * tke_0 - 4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xzzz[k] = 9.0 * to_xxz_xzz[k] - 6.0 * to_xxz_xzzzz[k] * tke_0 - 6.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yyyy[k] = -6.0 * to_xxz_yyyyz[k] * tke_0 + 4.0 * to_xxxxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yyyz[k] = 3.0 * to_xxz_yyy[k] - 6.0 * to_xxz_yyyzz[k] * tke_0 - 2.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yyzz[k] = 6.0 * to_xxz_yyz[k] - 6.0 * to_xxz_yyzzz[k] * tke_0 - 4.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yzzz[k] = 9.0 * to_xxz_yzz[k] - 6.0 * to_xxz_yzzzz[k] * tke_0 - 6.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_zzzz[k] = 12.0 * to_xxz_zzz[k] - 6.0 * to_xxz_zzzzz[k] * tke_0 - 8.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 495-510 components of targeted buffer : GG

        auto to_x_z_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 45);

        auto to_x_z_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 46);

        auto to_x_z_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 47);

        auto to_x_z_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 48);

        auto to_x_z_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 49);

        auto to_x_z_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 50);

        auto to_x_z_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 51);

        auto to_x_z_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 52);

        auto to_x_z_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 53);

        auto to_x_z_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 54);

        auto to_x_z_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 55);

        auto to_x_z_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 56);

        auto to_x_z_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 57);

        auto to_x_z_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 58);

        auto to_x_z_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_x_z_xxyy_xxxx, to_x_z_xxyy_xxxy, to_x_z_xxyy_xxxz, to_x_z_xxyy_xxyy, to_x_z_xxyy_xxyz, to_x_z_xxyy_xxzz, to_x_z_xxyy_xyyy, to_x_z_xxyy_xyyz, to_x_z_xxyy_xyzz, to_x_z_xxyy_xzzz, to_x_z_xxyy_yyyy, to_x_z_xxyy_yyyz, to_x_z_xxyy_yyzz, to_x_z_xxyy_yzzz, to_x_z_xxyy_zzzz, to_xxxyy_xxx, to_xxxyy_xxxxz, to_xxxyy_xxxyz, to_xxxyy_xxxzz, to_xxxyy_xxy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xxzzz, to_xxxyy_xyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_xzzzz, to_xxxyy_yyy, to_xxxyy_yyyyz, to_xxxyy_yyyzz, to_xxxyy_yyz, to_xxxyy_yyzzz, to_xxxyy_yzz, to_xxxyy_yzzzz, to_xxxyy_zzz, to_xxxyy_zzzzz, to_xyy_xxx, to_xyy_xxxxz, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyy_xxxx[k] = -4.0 * to_xyy_xxxxz[k] * tke_0 + 4.0 * to_xxxyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxxy[k] = -4.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxxz[k] = 2.0 * to_xyy_xxx[k] - 4.0 * to_xyy_xxxzz[k] * tke_0 - 2.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxyy[k] = -4.0 * to_xyy_xxyyz[k] * tke_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxyz[k] = 2.0 * to_xyy_xxy[k] - 4.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxzz[k] = 4.0 * to_xyy_xxz[k] - 4.0 * to_xyy_xxzzz[k] * tke_0 - 4.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xyyy[k] = -4.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xyyz[k] = 2.0 * to_xyy_xyy[k] - 4.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xyzz[k] = 4.0 * to_xyy_xyz[k] - 4.0 * to_xyy_xyzzz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xzzz[k] = 6.0 * to_xyy_xzz[k] - 4.0 * to_xyy_xzzzz[k] * tke_0 - 6.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yyyy[k] = -4.0 * to_xyy_yyyyz[k] * tke_0 + 4.0 * to_xxxyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yyyz[k] = 2.0 * to_xyy_yyy[k] - 4.0 * to_xyy_yyyzz[k] * tke_0 - 2.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yyzz[k] = 4.0 * to_xyy_yyz[k] - 4.0 * to_xyy_yyzzz[k] * tke_0 - 4.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yzzz[k] = 6.0 * to_xyy_yzz[k] - 4.0 * to_xyy_yzzzz[k] * tke_0 - 6.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_zzzz[k] = 8.0 * to_xyy_zzz[k] - 4.0 * to_xyy_zzzzz[k] * tke_0 - 8.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-525 components of targeted buffer : GG

        auto to_x_z_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 60);

        auto to_x_z_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 61);

        auto to_x_z_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 62);

        auto to_x_z_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 63);

        auto to_x_z_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 64);

        auto to_x_z_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 65);

        auto to_x_z_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 66);

        auto to_x_z_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 67);

        auto to_x_z_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 68);

        auto to_x_z_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 69);

        auto to_x_z_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 70);

        auto to_x_z_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 71);

        auto to_x_z_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 72);

        auto to_x_z_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 73);

        auto to_x_z_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_x_z_xxyz_xxxx, to_x_z_xxyz_xxxy, to_x_z_xxyz_xxxz, to_x_z_xxyz_xxyy, to_x_z_xxyz_xxyz, to_x_z_xxyz_xxzz, to_x_z_xxyz_xyyy, to_x_z_xxyz_xyyz, to_x_z_xxyz_xyzz, to_x_z_xxyz_xzzz, to_x_z_xxyz_yyyy, to_x_z_xxyz_yyyz, to_x_z_xxyz_yyzz, to_x_z_xxyz_yzzz, to_x_z_xxyz_zzzz, to_xxxyz_xxx, to_xxxyz_xxxxz, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_xxxyz_zzzzz, to_xyz_xxx, to_xyz_xxxxz, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_xyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyz_xxxx[k] = -4.0 * to_xyz_xxxxz[k] * tke_0 + 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxxy[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxxz[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxyy[k] = -4.0 * to_xyz_xxyyz[k] * tke_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxyz[k] = 2.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxzz[k] = 4.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 4.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xyyy[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xyyz[k] = 2.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xyzz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyzzz[k] * tke_0 - 4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xzzz[k] = 6.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xzzzz[k] * tke_0 - 6.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yyyy[k] = -4.0 * to_xyz_yyyyz[k] * tke_0 + 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yyyz[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yyzz[k] = 4.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 4.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yzzz[k] = 6.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yzzzz[k] * tke_0 - 6.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_zzzz[k] = 8.0 * to_xyz_zzz[k] - 4.0 * to_xyz_zzzzz[k] * tke_0 - 8.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 525-540 components of targeted buffer : GG

        auto to_x_z_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 75);

        auto to_x_z_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 76);

        auto to_x_z_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 77);

        auto to_x_z_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 78);

        auto to_x_z_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 79);

        auto to_x_z_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 80);

        auto to_x_z_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 81);

        auto to_x_z_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 82);

        auto to_x_z_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 83);

        auto to_x_z_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 84);

        auto to_x_z_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 85);

        auto to_x_z_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 86);

        auto to_x_z_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 87);

        auto to_x_z_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 88);

        auto to_x_z_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_x_z_xxzz_xxxx, to_x_z_xxzz_xxxy, to_x_z_xxzz_xxxz, to_x_z_xxzz_xxyy, to_x_z_xxzz_xxyz, to_x_z_xxzz_xxzz, to_x_z_xxzz_xyyy, to_x_z_xxzz_xyyz, to_x_z_xxzz_xyzz, to_x_z_xxzz_xzzz, to_x_z_xxzz_yyyy, to_x_z_xxzz_yyyz, to_x_z_xxzz_yyzz, to_x_z_xxzz_yzzz, to_x_z_xxzz_zzzz, to_xxxzz_xxx, to_xxxzz_xxxxz, to_xxxzz_xxxyz, to_xxxzz_xxxzz, to_xxxzz_xxy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xxzzz, to_xxxzz_xyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_xzzzz, to_xxxzz_yyy, to_xxxzz_yyyyz, to_xxxzz_yyyzz, to_xxxzz_yyz, to_xxxzz_yyzzz, to_xxxzz_yzz, to_xxxzz_yzzzz, to_xxxzz_zzz, to_xxxzz_zzzzz, to_xzz_xxx, to_xzz_xxxxz, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_xzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxzz_xxxx[k] = -4.0 * to_xzz_xxxxz[k] * tke_0 + 4.0 * to_xxxzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxxy[k] = -4.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxxz[k] = 2.0 * to_xzz_xxx[k] - 4.0 * to_xzz_xxxzz[k] * tke_0 - 2.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxyy[k] = -4.0 * to_xzz_xxyyz[k] * tke_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxyz[k] = 2.0 * to_xzz_xxy[k] - 4.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxzz[k] = 4.0 * to_xzz_xxz[k] - 4.0 * to_xzz_xxzzz[k] * tke_0 - 4.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xyyy[k] = -4.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xyyz[k] = 2.0 * to_xzz_xyy[k] - 4.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xyzz[k] = 4.0 * to_xzz_xyz[k] - 4.0 * to_xzz_xyzzz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xzzz[k] = 6.0 * to_xzz_xzz[k] - 4.0 * to_xzz_xzzzz[k] * tke_0 - 6.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yyyy[k] = -4.0 * to_xzz_yyyyz[k] * tke_0 + 4.0 * to_xxxzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yyyz[k] = 2.0 * to_xzz_yyy[k] - 4.0 * to_xzz_yyyzz[k] * tke_0 - 2.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yyzz[k] = 4.0 * to_xzz_yyz[k] - 4.0 * to_xzz_yyzzz[k] * tke_0 - 4.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yzzz[k] = 6.0 * to_xzz_yzz[k] - 4.0 * to_xzz_yzzzz[k] * tke_0 - 6.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_zzzz[k] = 8.0 * to_xzz_zzz[k] - 4.0 * to_xzz_zzzzz[k] * tke_0 - 8.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-555 components of targeted buffer : GG

        auto to_x_z_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 90);

        auto to_x_z_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 91);

        auto to_x_z_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 92);

        auto to_x_z_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 93);

        auto to_x_z_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 94);

        auto to_x_z_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 95);

        auto to_x_z_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 96);

        auto to_x_z_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 97);

        auto to_x_z_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 98);

        auto to_x_z_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 99);

        auto to_x_z_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 100);

        auto to_x_z_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 101);

        auto to_x_z_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 102);

        auto to_x_z_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 103);

        auto to_x_z_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_x_z_xyyy_xxxx, to_x_z_xyyy_xxxy, to_x_z_xyyy_xxxz, to_x_z_xyyy_xxyy, to_x_z_xyyy_xxyz, to_x_z_xyyy_xxzz, to_x_z_xyyy_xyyy, to_x_z_xyyy_xyyz, to_x_z_xyyy_xyzz, to_x_z_xyyy_xzzz, to_x_z_xyyy_yyyy, to_x_z_xyyy_yyyz, to_x_z_xyyy_yyzz, to_x_z_xyyy_yzzz, to_x_z_xyyy_zzzz, to_xxyyy_xxx, to_xxyyy_xxxxz, to_xxyyy_xxxyz, to_xxyyy_xxxzz, to_xxyyy_xxy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xxzzz, to_xxyyy_xyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_xzzzz, to_xxyyy_yyy, to_xxyyy_yyyyz, to_xxyyy_yyyzz, to_xxyyy_yyz, to_xxyyy_yyzzz, to_xxyyy_yzz, to_xxyyy_yzzzz, to_xxyyy_zzz, to_xxyyy_zzzzz, to_yyy_xxx, to_yyy_xxxxz, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyy_xxxx[k] = -2.0 * to_yyy_xxxxz[k] * tke_0 + 4.0 * to_xxyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxxy[k] = -2.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxxz[k] = to_yyy_xxx[k] - 2.0 * to_yyy_xxxzz[k] * tke_0 - 2.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxyy[k] = -2.0 * to_yyy_xxyyz[k] * tke_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxyz[k] = to_yyy_xxy[k] - 2.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxzz[k] = 2.0 * to_yyy_xxz[k] - 2.0 * to_yyy_xxzzz[k] * tke_0 - 4.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xyyy[k] = -2.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xyyz[k] = to_yyy_xyy[k] - 2.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xyzz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xyzzz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xzzz[k] = 3.0 * to_yyy_xzz[k] - 2.0 * to_yyy_xzzzz[k] * tke_0 - 6.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yyyy[k] = -2.0 * to_yyy_yyyyz[k] * tke_0 + 4.0 * to_xxyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yyyz[k] = to_yyy_yyy[k] - 2.0 * to_yyy_yyyzz[k] * tke_0 - 2.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yyzz[k] = 2.0 * to_yyy_yyz[k] - 2.0 * to_yyy_yyzzz[k] * tke_0 - 4.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yzzz[k] = 3.0 * to_yyy_yzz[k] - 2.0 * to_yyy_yzzzz[k] * tke_0 - 6.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_zzzz[k] = 4.0 * to_yyy_zzz[k] - 2.0 * to_yyy_zzzzz[k] * tke_0 - 8.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 555-570 components of targeted buffer : GG

        auto to_x_z_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 105);

        auto to_x_z_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 106);

        auto to_x_z_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 107);

        auto to_x_z_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 108);

        auto to_x_z_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 109);

        auto to_x_z_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 110);

        auto to_x_z_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 111);

        auto to_x_z_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 112);

        auto to_x_z_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 113);

        auto to_x_z_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 114);

        auto to_x_z_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 115);

        auto to_x_z_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 116);

        auto to_x_z_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 117);

        auto to_x_z_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 118);

        auto to_x_z_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_x_z_xyyz_xxxx, to_x_z_xyyz_xxxy, to_x_z_xyyz_xxxz, to_x_z_xyyz_xxyy, to_x_z_xyyz_xxyz, to_x_z_xyyz_xxzz, to_x_z_xyyz_xyyy, to_x_z_xyyz_xyyz, to_x_z_xyyz_xyzz, to_x_z_xyyz_xzzz, to_x_z_xyyz_yyyy, to_x_z_xyyz_yyyz, to_x_z_xyyz_yyzz, to_x_z_xyyz_yzzz, to_x_z_xyyz_zzzz, to_xxyyz_xxx, to_xxyyz_xxxxz, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_xxyyz_zzzzz, to_yyz_xxx, to_yyz_xxxxz, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, to_yyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyz_xxxx[k] = -2.0 * to_yyz_xxxxz[k] * tke_0 + 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxxy[k] = -2.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxxz[k] = to_yyz_xxx[k] - 2.0 * to_yyz_xxxzz[k] * tke_0 - 2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxyy[k] = -2.0 * to_yyz_xxyyz[k] * tke_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxyz[k] = to_yyz_xxy[k] - 2.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxzz[k] = 2.0 * to_yyz_xxz[k] - 2.0 * to_yyz_xxzzz[k] * tke_0 - 4.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xyyy[k] = -2.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xyyz[k] = to_yyz_xyy[k] - 2.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xyzz[k] = 2.0 * to_yyz_xyz[k] - 2.0 * to_yyz_xyzzz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xzzz[k] = 3.0 * to_yyz_xzz[k] - 2.0 * to_yyz_xzzzz[k] * tke_0 - 6.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yyyy[k] = -2.0 * to_yyz_yyyyz[k] * tke_0 + 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yyyz[k] = to_yyz_yyy[k] - 2.0 * to_yyz_yyyzz[k] * tke_0 - 2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yyzz[k] = 2.0 * to_yyz_yyz[k] - 2.0 * to_yyz_yyzzz[k] * tke_0 - 4.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yzzz[k] = 3.0 * to_yyz_yzz[k] - 2.0 * to_yyz_yzzzz[k] * tke_0 - 6.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_zzzz[k] = 4.0 * to_yyz_zzz[k] - 2.0 * to_yyz_zzzzz[k] * tke_0 - 8.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-585 components of targeted buffer : GG

        auto to_x_z_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 120);

        auto to_x_z_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 121);

        auto to_x_z_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 122);

        auto to_x_z_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 123);

        auto to_x_z_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 124);

        auto to_x_z_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 125);

        auto to_x_z_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 126);

        auto to_x_z_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 127);

        auto to_x_z_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 128);

        auto to_x_z_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 129);

        auto to_x_z_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 130);

        auto to_x_z_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 131);

        auto to_x_z_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 132);

        auto to_x_z_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 133);

        auto to_x_z_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_x_z_xyzz_xxxx, to_x_z_xyzz_xxxy, to_x_z_xyzz_xxxz, to_x_z_xyzz_xxyy, to_x_z_xyzz_xxyz, to_x_z_xyzz_xxzz, to_x_z_xyzz_xyyy, to_x_z_xyzz_xyyz, to_x_z_xyzz_xyzz, to_x_z_xyzz_xzzz, to_x_z_xyzz_yyyy, to_x_z_xyzz_yyyz, to_x_z_xyzz_yyzz, to_x_z_xyzz_yzzz, to_x_z_xyzz_zzzz, to_xxyzz_xxx, to_xxyzz_xxxxz, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_xxyzz_zzzzz, to_yzz_xxx, to_yzz_xxxxz, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, to_yzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyzz_xxxx[k] = -2.0 * to_yzz_xxxxz[k] * tke_0 + 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxxy[k] = -2.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxxz[k] = to_yzz_xxx[k] - 2.0 * to_yzz_xxxzz[k] * tke_0 - 2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxyy[k] = -2.0 * to_yzz_xxyyz[k] * tke_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxyz[k] = to_yzz_xxy[k] - 2.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxzz[k] = 2.0 * to_yzz_xxz[k] - 2.0 * to_yzz_xxzzz[k] * tke_0 - 4.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xyyy[k] = -2.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xyyz[k] = to_yzz_xyy[k] - 2.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xyzz[k] = 2.0 * to_yzz_xyz[k] - 2.0 * to_yzz_xyzzz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xzzz[k] = 3.0 * to_yzz_xzz[k] - 2.0 * to_yzz_xzzzz[k] * tke_0 - 6.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yyyy[k] = -2.0 * to_yzz_yyyyz[k] * tke_0 + 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yyyz[k] = to_yzz_yyy[k] - 2.0 * to_yzz_yyyzz[k] * tke_0 - 2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yyzz[k] = 2.0 * to_yzz_yyz[k] - 2.0 * to_yzz_yyzzz[k] * tke_0 - 4.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yzzz[k] = 3.0 * to_yzz_yzz[k] - 2.0 * to_yzz_yzzzz[k] * tke_0 - 6.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_zzzz[k] = 4.0 * to_yzz_zzz[k] - 2.0 * to_yzz_zzzzz[k] * tke_0 - 8.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 585-600 components of targeted buffer : GG

        auto to_x_z_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 135);

        auto to_x_z_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 136);

        auto to_x_z_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 137);

        auto to_x_z_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 138);

        auto to_x_z_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 139);

        auto to_x_z_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 140);

        auto to_x_z_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 141);

        auto to_x_z_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 142);

        auto to_x_z_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 143);

        auto to_x_z_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 144);

        auto to_x_z_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 145);

        auto to_x_z_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 146);

        auto to_x_z_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 147);

        auto to_x_z_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 148);

        auto to_x_z_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_x_z_xzzz_xxxx, to_x_z_xzzz_xxxy, to_x_z_xzzz_xxxz, to_x_z_xzzz_xxyy, to_x_z_xzzz_xxyz, to_x_z_xzzz_xxzz, to_x_z_xzzz_xyyy, to_x_z_xzzz_xyyz, to_x_z_xzzz_xyzz, to_x_z_xzzz_xzzz, to_x_z_xzzz_yyyy, to_x_z_xzzz_yyyz, to_x_z_xzzz_yyzz, to_x_z_xzzz_yzzz, to_x_z_xzzz_zzzz, to_xxzzz_xxx, to_xxzzz_xxxxz, to_xxzzz_xxxyz, to_xxzzz_xxxzz, to_xxzzz_xxy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xxzzz, to_xxzzz_xyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_xzzzz, to_xxzzz_yyy, to_xxzzz_yyyyz, to_xxzzz_yyyzz, to_xxzzz_yyz, to_xxzzz_yyzzz, to_xxzzz_yzz, to_xxzzz_yzzzz, to_xxzzz_zzz, to_xxzzz_zzzzz, to_zzz_xxx, to_zzz_xxxxz, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, to_zzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzzz_xxxx[k] = -2.0 * to_zzz_xxxxz[k] * tke_0 + 4.0 * to_xxzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxxy[k] = -2.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxxz[k] = to_zzz_xxx[k] - 2.0 * to_zzz_xxxzz[k] * tke_0 - 2.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxyy[k] = -2.0 * to_zzz_xxyyz[k] * tke_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxyz[k] = to_zzz_xxy[k] - 2.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxzz[k] = 2.0 * to_zzz_xxz[k] - 2.0 * to_zzz_xxzzz[k] * tke_0 - 4.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xyyy[k] = -2.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xyyz[k] = to_zzz_xyy[k] - 2.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xyzz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xyzzz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xzzz[k] = 3.0 * to_zzz_xzz[k] - 2.0 * to_zzz_xzzzz[k] * tke_0 - 6.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yyyy[k] = -2.0 * to_zzz_yyyyz[k] * tke_0 + 4.0 * to_xxzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yyyz[k] = to_zzz_yyy[k] - 2.0 * to_zzz_yyyzz[k] * tke_0 - 2.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yyzz[k] = 2.0 * to_zzz_yyz[k] - 2.0 * to_zzz_yyzzz[k] * tke_0 - 4.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yzzz[k] = 3.0 * to_zzz_yzz[k] - 2.0 * to_zzz_yzzzz[k] * tke_0 - 6.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_zzzz[k] = 4.0 * to_zzz_zzz[k] - 2.0 * to_zzz_zzzzz[k] * tke_0 - 8.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-615 components of targeted buffer : GG

        auto to_x_z_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 150);

        auto to_x_z_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 151);

        auto to_x_z_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 152);

        auto to_x_z_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 153);

        auto to_x_z_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 154);

        auto to_x_z_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 155);

        auto to_x_z_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 156);

        auto to_x_z_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 157);

        auto to_x_z_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 158);

        auto to_x_z_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 159);

        auto to_x_z_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 160);

        auto to_x_z_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 161);

        auto to_x_z_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 162);

        auto to_x_z_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 163);

        auto to_x_z_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_x_z_yyyy_xxxx, to_x_z_yyyy_xxxy, to_x_z_yyyy_xxxz, to_x_z_yyyy_xxyy, to_x_z_yyyy_xxyz, to_x_z_yyyy_xxzz, to_x_z_yyyy_xyyy, to_x_z_yyyy_xyyz, to_x_z_yyyy_xyzz, to_x_z_yyyy_xzzz, to_x_z_yyyy_yyyy, to_x_z_yyyy_yyyz, to_x_z_yyyy_yyzz, to_x_z_yyyy_yzzz, to_x_z_yyyy_zzzz, to_xyyyy_xxx, to_xyyyy_xxxxz, to_xyyyy_xxxyz, to_xyyyy_xxxzz, to_xyyyy_xxy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xxzzz, to_xyyyy_xyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_xzzzz, to_xyyyy_yyy, to_xyyyy_yyyyz, to_xyyyy_yyyzz, to_xyyyy_yyz, to_xyyyy_yyzzz, to_xyyyy_yzz, to_xyyyy_yzzzz, to_xyyyy_zzz, to_xyyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyy_xxxx[k] = 4.0 * to_xyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxxy[k] = 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxxz[k] = -2.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxyy[k] = 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxyz[k] = -2.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxzz[k] = -4.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xyyy[k] = 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xyyz[k] = -2.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xyzz[k] = -4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xzzz[k] = -6.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yyyy[k] = 4.0 * to_xyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yyyz[k] = -2.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yyzz[k] = -4.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yzzz[k] = -6.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_zzzz[k] = -8.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 615-630 components of targeted buffer : GG

        auto to_x_z_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 165);

        auto to_x_z_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 166);

        auto to_x_z_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 167);

        auto to_x_z_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 168);

        auto to_x_z_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 169);

        auto to_x_z_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 170);

        auto to_x_z_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 171);

        auto to_x_z_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 172);

        auto to_x_z_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 173);

        auto to_x_z_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 174);

        auto to_x_z_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 175);

        auto to_x_z_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 176);

        auto to_x_z_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 177);

        auto to_x_z_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 178);

        auto to_x_z_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_x_z_yyyz_xxxx, to_x_z_yyyz_xxxy, to_x_z_yyyz_xxxz, to_x_z_yyyz_xxyy, to_x_z_yyyz_xxyz, to_x_z_yyyz_xxzz, to_x_z_yyyz_xyyy, to_x_z_yyyz_xyyz, to_x_z_yyyz_xyzz, to_x_z_yyyz_xzzz, to_x_z_yyyz_yyyy, to_x_z_yyyz_yyyz, to_x_z_yyyz_yyzz, to_x_z_yyyz_yzzz, to_x_z_yyyz_zzzz, to_xyyyz_xxx, to_xyyyz_xxxxz, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, to_xyyyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyz_xxxx[k] = 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxxy[k] = 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxxz[k] = -2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxyy[k] = 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxyz[k] = -2.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxzz[k] = -4.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xyyy[k] = 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xyyz[k] = -2.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xyzz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xzzz[k] = -6.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yyyy[k] = 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yyyz[k] = -2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yyzz[k] = -4.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yzzz[k] = -6.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_zzzz[k] = -8.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-645 components of targeted buffer : GG

        auto to_x_z_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 180);

        auto to_x_z_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 181);

        auto to_x_z_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 182);

        auto to_x_z_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 183);

        auto to_x_z_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 184);

        auto to_x_z_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 185);

        auto to_x_z_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 186);

        auto to_x_z_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 187);

        auto to_x_z_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 188);

        auto to_x_z_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 189);

        auto to_x_z_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 190);

        auto to_x_z_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 191);

        auto to_x_z_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 192);

        auto to_x_z_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 193);

        auto to_x_z_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_x_z_yyzz_xxxx, to_x_z_yyzz_xxxy, to_x_z_yyzz_xxxz, to_x_z_yyzz_xxyy, to_x_z_yyzz_xxyz, to_x_z_yyzz_xxzz, to_x_z_yyzz_xyyy, to_x_z_yyzz_xyyz, to_x_z_yyzz_xyzz, to_x_z_yyzz_xzzz, to_x_z_yyzz_yyyy, to_x_z_yyzz_yyyz, to_x_z_yyzz_yyzz, to_x_z_yyzz_yzzz, to_x_z_yyzz_zzzz, to_xyyzz_xxx, to_xyyzz_xxxxz, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, to_xyyzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyzz_xxxx[k] = 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxxy[k] = 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxxz[k] = -2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxyy[k] = 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxyz[k] = -2.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxzz[k] = -4.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xyyy[k] = 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xyyz[k] = -2.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xyzz[k] = -4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xzzz[k] = -6.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yyyy[k] = 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yyyz[k] = -2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yyzz[k] = -4.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yzzz[k] = -6.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_zzzz[k] = -8.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 645-660 components of targeted buffer : GG

        auto to_x_z_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 195);

        auto to_x_z_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 196);

        auto to_x_z_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 197);

        auto to_x_z_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 198);

        auto to_x_z_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 199);

        auto to_x_z_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 200);

        auto to_x_z_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 201);

        auto to_x_z_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 202);

        auto to_x_z_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 203);

        auto to_x_z_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 204);

        auto to_x_z_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 205);

        auto to_x_z_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 206);

        auto to_x_z_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 207);

        auto to_x_z_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 208);

        auto to_x_z_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_x_z_yzzz_xxxx, to_x_z_yzzz_xxxy, to_x_z_yzzz_xxxz, to_x_z_yzzz_xxyy, to_x_z_yzzz_xxyz, to_x_z_yzzz_xxzz, to_x_z_yzzz_xyyy, to_x_z_yzzz_xyyz, to_x_z_yzzz_xyzz, to_x_z_yzzz_xzzz, to_x_z_yzzz_yyyy, to_x_z_yzzz_yyyz, to_x_z_yzzz_yyzz, to_x_z_yzzz_yzzz, to_x_z_yzzz_zzzz, to_xyzzz_xxx, to_xyzzz_xxxxz, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, to_xyzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzzz_xxxx[k] = 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxxy[k] = 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxxz[k] = -2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxyy[k] = 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxyz[k] = -2.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxzz[k] = -4.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xyyy[k] = 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xyyz[k] = -2.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xyzz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xzzz[k] = -6.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yyyy[k] = 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yyyz[k] = -2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yyzz[k] = -4.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yzzz[k] = -6.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_zzzz[k] = -8.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-675 components of targeted buffer : GG

        auto to_x_z_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 210);

        auto to_x_z_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 211);

        auto to_x_z_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 212);

        auto to_x_z_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 213);

        auto to_x_z_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 214);

        auto to_x_z_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 215);

        auto to_x_z_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 216);

        auto to_x_z_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 217);

        auto to_x_z_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 218);

        auto to_x_z_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 219);

        auto to_x_z_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 220);

        auto to_x_z_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 221);

        auto to_x_z_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 222);

        auto to_x_z_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 223);

        auto to_x_z_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 2 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_x_z_zzzz_xxxx, to_x_z_zzzz_xxxy, to_x_z_zzzz_xxxz, to_x_z_zzzz_xxyy, to_x_z_zzzz_xxyz, to_x_z_zzzz_xxzz, to_x_z_zzzz_xyyy, to_x_z_zzzz_xyyz, to_x_z_zzzz_xyzz, to_x_z_zzzz_xzzz, to_x_z_zzzz_yyyy, to_x_z_zzzz_yyyz, to_x_z_zzzz_yyzz, to_x_z_zzzz_yzzz, to_x_z_zzzz_zzzz, to_xzzzz_xxx, to_xzzzz_xxxxz, to_xzzzz_xxxyz, to_xzzzz_xxxzz, to_xzzzz_xxy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xxzzz, to_xzzzz_xyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_xzzzz, to_xzzzz_yyy, to_xzzzz_yyyyz, to_xzzzz_yyyzz, to_xzzzz_yyz, to_xzzzz_yyzzz, to_xzzzz_yzz, to_xzzzz_yzzzz, to_xzzzz_zzz, to_xzzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzzz_xxxx[k] = 4.0 * to_xzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxxy[k] = 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxxz[k] = -2.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxyy[k] = 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxyz[k] = -2.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxzz[k] = -4.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xyyy[k] = 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xyyz[k] = -2.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xyzz[k] = -4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xzzz[k] = -6.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yyyy[k] = 4.0 * to_xzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yyyz[k] = -2.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yyzz[k] = -4.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yzzz[k] = -6.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_zzzz[k] = -8.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 675-690 components of targeted buffer : GG

        auto to_y_x_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 0);

        auto to_y_x_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 1);

        auto to_y_x_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 2);

        auto to_y_x_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 3);

        auto to_y_x_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 4);

        auto to_y_x_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 5);

        auto to_y_x_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 6);

        auto to_y_x_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 7);

        auto to_y_x_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 8);

        auto to_y_x_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 9);

        auto to_y_x_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 10);

        auto to_y_x_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 11);

        auto to_y_x_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 12);

        auto to_y_x_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 13);

        auto to_y_x_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxy_xxx, to_xxxxy_xxxxx, to_xxxxy_xxxxy, to_xxxxy_xxxxz, to_xxxxy_xxxyy, to_xxxxy_xxxyz, to_xxxxy_xxxzz, to_xxxxy_xxy, to_xxxxy_xxyyy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xxzzz, to_xxxxy_xyy, to_xxxxy_xyyyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_xzzzz, to_xxxxy_yyy, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_zzz, to_y_x_xxxx_xxxx, to_y_x_xxxx_xxxy, to_y_x_xxxx_xxxz, to_y_x_xxxx_xxyy, to_y_x_xxxx_xxyz, to_y_x_xxxx_xxzz, to_y_x_xxxx_xyyy, to_y_x_xxxx_xyyz, to_y_x_xxxx_xyzz, to_y_x_xxxx_xzzz, to_y_x_xxxx_yyyy, to_y_x_xxxx_yyyz, to_y_x_xxxx_yyzz, to_y_x_xxxx_yzzz, to_y_x_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxx_xxxx[k] = -8.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxxy[k] = -6.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxxz[k] = -6.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxyy[k] = -4.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxyz[k] = -4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxzz[k] = -4.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xyyy[k] = -2.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xyyz[k] = -2.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xyzz[k] = -2.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xzzz[k] = -2.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yyyy[k] = 4.0 * to_xxxxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yyyz[k] = 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yyzz[k] = 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yzzz[k] = 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_zzzz[k] = 4.0 * to_xxxxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-705 components of targeted buffer : GG

        auto to_y_x_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 15);

        auto to_y_x_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 16);

        auto to_y_x_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 17);

        auto to_y_x_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 18);

        auto to_y_x_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 19);

        auto to_y_x_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 20);

        auto to_y_x_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 21);

        auto to_y_x_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 22);

        auto to_y_x_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 23);

        auto to_y_x_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 24);

        auto to_y_x_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 25);

        auto to_y_x_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 26);

        auto to_y_x_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 27);

        auto to_y_x_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 28);

        auto to_y_x_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxx, to_xxx_xxxxy, to_xxx_xxxxz, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_zzz, to_xxxyy_xxx, to_xxxyy_xxxxx, to_xxxyy_xxxxy, to_xxxyy_xxxxz, to_xxxyy_xxxyy, to_xxxyy_xxxyz, to_xxxyy_xxxzz, to_xxxyy_xxy, to_xxxyy_xxyyy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xxzzz, to_xxxyy_xyy, to_xxxyy_xyyyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_xzzzz, to_xxxyy_yyy, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_zzz, to_y_x_xxxy_xxxx, to_y_x_xxxy_xxxy, to_y_x_xxxy_xxxz, to_y_x_xxxy_xxyy, to_y_x_xxxy_xxyz, to_y_x_xxxy_xxzz, to_y_x_xxxy_xyyy, to_y_x_xxxy_xyyz, to_y_x_xxxy_xyzz, to_y_x_xxxy_xzzz, to_y_x_xxxy_yyyy, to_y_x_xxxy_yyyz, to_y_x_xxxy_yyzz, to_y_x_xxxy_yzzz, to_y_x_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxy_xxxx[k] = 4.0 * to_xxx_xxx[k] - 2.0 * to_xxx_xxxxx[k] * tke_0 - 8.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxxy[k] = 3.0 * to_xxx_xxy[k] - 2.0 * to_xxx_xxxxy[k] * tke_0 - 6.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxxz[k] = 3.0 * to_xxx_xxz[k] - 2.0 * to_xxx_xxxxz[k] * tke_0 - 6.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxyy[k] = 2.0 * to_xxx_xyy[k] - 2.0 * to_xxx_xxxyy[k] * tke_0 - 4.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxyz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xxxyz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxzz[k] = 2.0 * to_xxx_xzz[k] - 2.0 * to_xxx_xxxzz[k] * tke_0 - 4.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xyyy[k] = to_xxx_yyy[k] - 2.0 * to_xxx_xxyyy[k] * tke_0 - 2.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xyyz[k] = to_xxx_yyz[k] - 2.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xyzz[k] = to_xxx_yzz[k] - 2.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xzzz[k] = to_xxx_zzz[k] - 2.0 * to_xxx_xxzzz[k] * tke_0 - 2.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yyyy[k] = -2.0 * to_xxx_xyyyy[k] * tke_0 + 4.0 * to_xxxyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yyyz[k] = -2.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yyzz[k] = -2.0 * to_xxx_xyyzz[k] * tke_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yzzz[k] = -2.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_zzzz[k] = -2.0 * to_xxx_xzzzz[k] * tke_0 + 4.0 * to_xxxyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 705-720 components of targeted buffer : GG

        auto to_y_x_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 30);

        auto to_y_x_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 31);

        auto to_y_x_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 32);

        auto to_y_x_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 33);

        auto to_y_x_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 34);

        auto to_y_x_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 35);

        auto to_y_x_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 36);

        auto to_y_x_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 37);

        auto to_y_x_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 38);

        auto to_y_x_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 39);

        auto to_y_x_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 40);

        auto to_y_x_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 41);

        auto to_y_x_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 42);

        auto to_y_x_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 43);

        auto to_y_x_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxx, to_xxxyz_xxxxy, to_xxxyz_xxxxz, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_zzz, to_y_x_xxxz_xxxx, to_y_x_xxxz_xxxy, to_y_x_xxxz_xxxz, to_y_x_xxxz_xxyy, to_y_x_xxxz_xxyz, to_y_x_xxxz_xxzz, to_y_x_xxxz_xyyy, to_y_x_xxxz_xyyz, to_y_x_xxxz_xyzz, to_y_x_xxxz_xzzz, to_y_x_xxxz_yyyy, to_y_x_xxxz_yyyz, to_y_x_xxxz_yyzz, to_y_x_xxxz_yzzz, to_y_x_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxz_xxxx[k] = -8.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxxy[k] = -6.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxxz[k] = -6.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxyy[k] = -4.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxyz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxzz[k] = -4.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xyyy[k] = -2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xyyz[k] = -2.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xyzz[k] = -2.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xzzz[k] = -2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yyyy[k] = 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yyyz[k] = 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yyzz[k] = 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yzzz[k] = 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_zzzz[k] = 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-735 components of targeted buffer : GG

        auto to_y_x_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 45);

        auto to_y_x_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 46);

        auto to_y_x_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 47);

        auto to_y_x_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 48);

        auto to_y_x_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 49);

        auto to_y_x_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 50);

        auto to_y_x_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 51);

        auto to_y_x_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 52);

        auto to_y_x_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 53);

        auto to_y_x_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 54);

        auto to_y_x_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 55);

        auto to_y_x_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 56);

        auto to_y_x_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 57);

        auto to_y_x_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 58);

        auto to_y_x_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxx, to_xxy_xxxxy, to_xxy_xxxxz, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_zzz, to_xxyyy_xxx, to_xxyyy_xxxxx, to_xxyyy_xxxxy, to_xxyyy_xxxxz, to_xxyyy_xxxyy, to_xxyyy_xxxyz, to_xxyyy_xxxzz, to_xxyyy_xxy, to_xxyyy_xxyyy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xxzzz, to_xxyyy_xyy, to_xxyyy_xyyyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_xzzzz, to_xxyyy_yyy, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_zzz, to_y_x_xxyy_xxxx, to_y_x_xxyy_xxxy, to_y_x_xxyy_xxxz, to_y_x_xxyy_xxyy, to_y_x_xxyy_xxyz, to_y_x_xxyy_xxzz, to_y_x_xxyy_xyyy, to_y_x_xxyy_xyyz, to_y_x_xxyy_xyzz, to_y_x_xxyy_xzzz, to_y_x_xxyy_yyyy, to_y_x_xxyy_yyyz, to_y_x_xxyy_yyzz, to_y_x_xxyy_yzzz, to_y_x_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyy_xxxx[k] = 8.0 * to_xxy_xxx[k] - 4.0 * to_xxy_xxxxx[k] * tke_0 - 8.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxxy[k] = 6.0 * to_xxy_xxy[k] - 4.0 * to_xxy_xxxxy[k] * tke_0 - 6.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxxz[k] = 6.0 * to_xxy_xxz[k] - 4.0 * to_xxy_xxxxz[k] * tke_0 - 6.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxyy[k] = 4.0 * to_xxy_xyy[k] - 4.0 * to_xxy_xxxyy[k] * tke_0 - 4.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxyz[k] = 4.0 * to_xxy_xyz[k] - 4.0 * to_xxy_xxxyz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxzz[k] = 4.0 * to_xxy_xzz[k] - 4.0 * to_xxy_xxxzz[k] * tke_0 - 4.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xyyy[k] = 2.0 * to_xxy_yyy[k] - 4.0 * to_xxy_xxyyy[k] * tke_0 - 2.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xyyz[k] = 2.0 * to_xxy_yyz[k] - 4.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xyzz[k] = 2.0 * to_xxy_yzz[k] - 4.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xzzz[k] = 2.0 * to_xxy_zzz[k] - 4.0 * to_xxy_xxzzz[k] * tke_0 - 2.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yyyy[k] = -4.0 * to_xxy_xyyyy[k] * tke_0 + 4.0 * to_xxyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yyyz[k] = -4.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yyzz[k] = -4.0 * to_xxy_xyyzz[k] * tke_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yzzz[k] = -4.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_zzzz[k] = -4.0 * to_xxy_xzzzz[k] * tke_0 + 4.0 * to_xxyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 735-750 components of targeted buffer : GG

        auto to_y_x_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 60);

        auto to_y_x_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 61);

        auto to_y_x_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 62);

        auto to_y_x_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 63);

        auto to_y_x_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 64);

        auto to_y_x_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 65);

        auto to_y_x_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 66);

        auto to_y_x_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 67);

        auto to_y_x_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 68);

        auto to_y_x_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 69);

        auto to_y_x_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 70);

        auto to_y_x_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 71);

        auto to_y_x_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 72);

        auto to_y_x_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 73);

        auto to_y_x_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxx, to_xxyyz_xxxxy, to_xxyyz_xxxxz, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_zzz, to_xxz_xxx, to_xxz_xxxxx, to_xxz_xxxxy, to_xxz_xxxxz, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_zzz, to_y_x_xxyz_xxxx, to_y_x_xxyz_xxxy, to_y_x_xxyz_xxxz, to_y_x_xxyz_xxyy, to_y_x_xxyz_xxyz, to_y_x_xxyz_xxzz, to_y_x_xxyz_xyyy, to_y_x_xxyz_xyyz, to_y_x_xxyz_xyzz, to_y_x_xxyz_xzzz, to_y_x_xxyz_yyyy, to_y_x_xxyz_yyyz, to_y_x_xxyz_yyzz, to_y_x_xxyz_yzzz, to_y_x_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyz_xxxx[k] = 4.0 * to_xxz_xxx[k] - 2.0 * to_xxz_xxxxx[k] * tke_0 - 8.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxxy[k] = 3.0 * to_xxz_xxy[k] - 2.0 * to_xxz_xxxxy[k] * tke_0 - 6.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxxz[k] = 3.0 * to_xxz_xxz[k] - 2.0 * to_xxz_xxxxz[k] * tke_0 - 6.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxyy[k] = 2.0 * to_xxz_xyy[k] - 2.0 * to_xxz_xxxyy[k] * tke_0 - 4.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxyz[k] = 2.0 * to_xxz_xyz[k] - 2.0 * to_xxz_xxxyz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxzz[k] = 2.0 * to_xxz_xzz[k] - 2.0 * to_xxz_xxxzz[k] * tke_0 - 4.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xyyy[k] = to_xxz_yyy[k] - 2.0 * to_xxz_xxyyy[k] * tke_0 - 2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xyyz[k] = to_xxz_yyz[k] - 2.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xyzz[k] = to_xxz_yzz[k] - 2.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xzzz[k] = to_xxz_zzz[k] - 2.0 * to_xxz_xxzzz[k] * tke_0 - 2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yyyy[k] = -2.0 * to_xxz_xyyyy[k] * tke_0 + 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yyyz[k] = -2.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yyzz[k] = -2.0 * to_xxz_xyyzz[k] * tke_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yzzz[k] = -2.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_zzzz[k] = -2.0 * to_xxz_xzzzz[k] * tke_0 + 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-765 components of targeted buffer : GG

        auto to_y_x_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 75);

        auto to_y_x_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 76);

        auto to_y_x_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 77);

        auto to_y_x_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 78);

        auto to_y_x_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 79);

        auto to_y_x_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 80);

        auto to_y_x_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 81);

        auto to_y_x_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 82);

        auto to_y_x_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 83);

        auto to_y_x_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 84);

        auto to_y_x_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 85);

        auto to_y_x_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 86);

        auto to_y_x_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 87);

        auto to_y_x_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 88);

        auto to_y_x_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxyzz_xxx, to_xxyzz_xxxxx, to_xxyzz_xxxxy, to_xxyzz_xxxxz, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_zzz, to_y_x_xxzz_xxxx, to_y_x_xxzz_xxxy, to_y_x_xxzz_xxxz, to_y_x_xxzz_xxyy, to_y_x_xxzz_xxyz, to_y_x_xxzz_xxzz, to_y_x_xxzz_xyyy, to_y_x_xxzz_xyyz, to_y_x_xxzz_xyzz, to_y_x_xxzz_xzzz, to_y_x_xxzz_yyyy, to_y_x_xxzz_yyyz, to_y_x_xxzz_yyzz, to_y_x_xxzz_yzzz, to_y_x_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxzz_xxxx[k] = -8.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxxy[k] = -6.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxxz[k] = -6.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxyy[k] = -4.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxyz[k] = -4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxzz[k] = -4.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xyyy[k] = -2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xyyz[k] = -2.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xyzz[k] = -2.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xzzz[k] = -2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yyyy[k] = 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yyyz[k] = 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yyzz[k] = 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yzzz[k] = 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_zzzz[k] = 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 765-780 components of targeted buffer : GG

        auto to_y_x_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 90);

        auto to_y_x_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 91);

        auto to_y_x_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 92);

        auto to_y_x_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 93);

        auto to_y_x_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 94);

        auto to_y_x_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 95);

        auto to_y_x_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 96);

        auto to_y_x_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 97);

        auto to_y_x_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 98);

        auto to_y_x_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 99);

        auto to_y_x_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 100);

        auto to_y_x_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 101);

        auto to_y_x_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 102);

        auto to_y_x_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 103);

        auto to_y_x_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxx, to_xyy_xxxxy, to_xyy_xxxxz, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_zzz, to_xyyyy_xxx, to_xyyyy_xxxxx, to_xyyyy_xxxxy, to_xyyyy_xxxxz, to_xyyyy_xxxyy, to_xyyyy_xxxyz, to_xyyyy_xxxzz, to_xyyyy_xxy, to_xyyyy_xxyyy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xxzzz, to_xyyyy_xyy, to_xyyyy_xyyyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_xzzzz, to_xyyyy_yyy, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_zzz, to_y_x_xyyy_xxxx, to_y_x_xyyy_xxxy, to_y_x_xyyy_xxxz, to_y_x_xyyy_xxyy, to_y_x_xyyy_xxyz, to_y_x_xyyy_xxzz, to_y_x_xyyy_xyyy, to_y_x_xyyy_xyyz, to_y_x_xyyy_xyzz, to_y_x_xyyy_xzzz, to_y_x_xyyy_yyyy, to_y_x_xyyy_yyyz, to_y_x_xyyy_yyzz, to_y_x_xyyy_yzzz, to_y_x_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyy_xxxx[k] = 12.0 * to_xyy_xxx[k] - 6.0 * to_xyy_xxxxx[k] * tke_0 - 8.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxxy[k] = 9.0 * to_xyy_xxy[k] - 6.0 * to_xyy_xxxxy[k] * tke_0 - 6.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxxz[k] = 9.0 * to_xyy_xxz[k] - 6.0 * to_xyy_xxxxz[k] * tke_0 - 6.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxyy[k] = 6.0 * to_xyy_xyy[k] - 6.0 * to_xyy_xxxyy[k] * tke_0 - 4.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxyz[k] = 6.0 * to_xyy_xyz[k] - 6.0 * to_xyy_xxxyz[k] * tke_0 - 4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxzz[k] = 6.0 * to_xyy_xzz[k] - 6.0 * to_xyy_xxxzz[k] * tke_0 - 4.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xyyy[k] = 3.0 * to_xyy_yyy[k] - 6.0 * to_xyy_xxyyy[k] * tke_0 - 2.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xyyz[k] = 3.0 * to_xyy_yyz[k] - 6.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xyzz[k] = 3.0 * to_xyy_yzz[k] - 6.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xzzz[k] = 3.0 * to_xyy_zzz[k] - 6.0 * to_xyy_xxzzz[k] * tke_0 - 2.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yyyy[k] = -6.0 * to_xyy_xyyyy[k] * tke_0 + 4.0 * to_xyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yyyz[k] = -6.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yyzz[k] = -6.0 * to_xyy_xyyzz[k] * tke_0 + 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yzzz[k] = -6.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_zzzz[k] = -6.0 * to_xyy_xzzzz[k] * tke_0 + 4.0 * to_xyyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-795 components of targeted buffer : GG

        auto to_y_x_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 105);

        auto to_y_x_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 106);

        auto to_y_x_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 107);

        auto to_y_x_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 108);

        auto to_y_x_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 109);

        auto to_y_x_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 110);

        auto to_y_x_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 111);

        auto to_y_x_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 112);

        auto to_y_x_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 113);

        auto to_y_x_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 114);

        auto to_y_x_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 115);

        auto to_y_x_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 116);

        auto to_y_x_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 117);

        auto to_y_x_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 118);

        auto to_y_x_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxx, to_xyyyz_xxxxy, to_xyyyz_xxxxz, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_zzz, to_xyz_xxx, to_xyz_xxxxx, to_xyz_xxxxy, to_xyz_xxxxz, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_zzz, to_y_x_xyyz_xxxx, to_y_x_xyyz_xxxy, to_y_x_xyyz_xxxz, to_y_x_xyyz_xxyy, to_y_x_xyyz_xxyz, to_y_x_xyyz_xxzz, to_y_x_xyyz_xyyy, to_y_x_xyyz_xyyz, to_y_x_xyyz_xyzz, to_y_x_xyyz_xzzz, to_y_x_xyyz_yyyy, to_y_x_xyyz_yyyz, to_y_x_xyyz_yyzz, to_y_x_xyyz_yzzz, to_y_x_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyz_xxxx[k] = 8.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxxx[k] * tke_0 - 8.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxxy[k] = 6.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxxxy[k] * tke_0 - 6.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxxz[k] = 6.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxxxz[k] * tke_0 - 6.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxyy[k] = 4.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 4.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xxxyz[k] * tke_0 - 4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxzz[k] = 4.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 4.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xyyy[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xyyz[k] = 2.0 * to_xyz_yyz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xyzz[k] = 2.0 * to_xyz_yzz[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yyyy[k] = -4.0 * to_xyz_xyyyy[k] * tke_0 + 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yyyz[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yyzz[k] = -4.0 * to_xyz_xyyzz[k] * tke_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_zzzz[k] = -4.0 * to_xyz_xzzzz[k] * tke_0 + 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 795-810 components of targeted buffer : GG

        auto to_y_x_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 120);

        auto to_y_x_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 121);

        auto to_y_x_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 122);

        auto to_y_x_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 123);

        auto to_y_x_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 124);

        auto to_y_x_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 125);

        auto to_y_x_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 126);

        auto to_y_x_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 127);

        auto to_y_x_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 128);

        auto to_y_x_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 129);

        auto to_y_x_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 130);

        auto to_y_x_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 131);

        auto to_y_x_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 132);

        auto to_y_x_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 133);

        auto to_y_x_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyyzz_xxx, to_xyyzz_xxxxx, to_xyyzz_xxxxy, to_xyyzz_xxxxz, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_zzz, to_xzz_xxx, to_xzz_xxxxx, to_xzz_xxxxy, to_xzz_xxxxz, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_zzz, to_y_x_xyzz_xxxx, to_y_x_xyzz_xxxy, to_y_x_xyzz_xxxz, to_y_x_xyzz_xxyy, to_y_x_xyzz_xxyz, to_y_x_xyzz_xxzz, to_y_x_xyzz_xyyy, to_y_x_xyzz_xyyz, to_y_x_xyzz_xyzz, to_y_x_xyzz_xzzz, to_y_x_xyzz_yyyy, to_y_x_xyzz_yyyz, to_y_x_xyzz_yyzz, to_y_x_xyzz_yzzz, to_y_x_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyzz_xxxx[k] = 4.0 * to_xzz_xxx[k] - 2.0 * to_xzz_xxxxx[k] * tke_0 - 8.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxxy[k] = 3.0 * to_xzz_xxy[k] - 2.0 * to_xzz_xxxxy[k] * tke_0 - 6.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxxz[k] = 3.0 * to_xzz_xxz[k] - 2.0 * to_xzz_xxxxz[k] * tke_0 - 6.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxyy[k] = 2.0 * to_xzz_xyy[k] - 2.0 * to_xzz_xxxyy[k] * tke_0 - 4.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxyz[k] = 2.0 * to_xzz_xyz[k] - 2.0 * to_xzz_xxxyz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxzz[k] = 2.0 * to_xzz_xzz[k] - 2.0 * to_xzz_xxxzz[k] * tke_0 - 4.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xyyy[k] = to_xzz_yyy[k] - 2.0 * to_xzz_xxyyy[k] * tke_0 - 2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xyyz[k] = to_xzz_yyz[k] - 2.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xyzz[k] = to_xzz_yzz[k] - 2.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xzzz[k] = to_xzz_zzz[k] - 2.0 * to_xzz_xxzzz[k] * tke_0 - 2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yyyy[k] = -2.0 * to_xzz_xyyyy[k] * tke_0 + 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yyyz[k] = -2.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yyzz[k] = -2.0 * to_xzz_xyyzz[k] * tke_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yzzz[k] = -2.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_zzzz[k] = -2.0 * to_xzz_xzzzz[k] * tke_0 + 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 810-825 components of targeted buffer : GG

        auto to_y_x_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 135);

        auto to_y_x_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 136);

        auto to_y_x_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 137);

        auto to_y_x_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 138);

        auto to_y_x_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 139);

        auto to_y_x_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 140);

        auto to_y_x_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 141);

        auto to_y_x_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 142);

        auto to_y_x_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 143);

        auto to_y_x_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 144);

        auto to_y_x_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 145);

        auto to_y_x_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 146);

        auto to_y_x_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 147);

        auto to_y_x_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 148);

        auto to_y_x_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xyzzz_xxx, to_xyzzz_xxxxx, to_xyzzz_xxxxy, to_xyzzz_xxxxz, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_zzz, to_y_x_xzzz_xxxx, to_y_x_xzzz_xxxy, to_y_x_xzzz_xxxz, to_y_x_xzzz_xxyy, to_y_x_xzzz_xxyz, to_y_x_xzzz_xxzz, to_y_x_xzzz_xyyy, to_y_x_xzzz_xyyz, to_y_x_xzzz_xyzz, to_y_x_xzzz_xzzz, to_y_x_xzzz_yyyy, to_y_x_xzzz_yyyz, to_y_x_xzzz_yyzz, to_y_x_xzzz_yzzz, to_y_x_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzzz_xxxx[k] = -8.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxxy[k] = -6.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxxz[k] = -6.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxyy[k] = -4.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxyz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxzz[k] = -4.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xyyy[k] = -2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xyyz[k] = -2.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xyzz[k] = -2.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xzzz[k] = -2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yyyy[k] = 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yyyz[k] = 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yyzz[k] = 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yzzz[k] = 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_zzzz[k] = 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 825-840 components of targeted buffer : GG

        auto to_y_x_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 150);

        auto to_y_x_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 151);

        auto to_y_x_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 152);

        auto to_y_x_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 153);

        auto to_y_x_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 154);

        auto to_y_x_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 155);

        auto to_y_x_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 156);

        auto to_y_x_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 157);

        auto to_y_x_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 158);

        auto to_y_x_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 159);

        auto to_y_x_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 160);

        auto to_y_x_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 161);

        auto to_y_x_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 162);

        auto to_y_x_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 163);

        auto to_y_x_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_y_x_yyyy_xxxx, to_y_x_yyyy_xxxy, to_y_x_yyyy_xxxz, to_y_x_yyyy_xxyy, to_y_x_yyyy_xxyz, to_y_x_yyyy_xxzz, to_y_x_yyyy_xyyy, to_y_x_yyyy_xyyz, to_y_x_yyyy_xyzz, to_y_x_yyyy_xzzz, to_y_x_yyyy_yyyy, to_y_x_yyyy_yyyz, to_y_x_yyyy_yyzz, to_y_x_yyyy_yzzz, to_y_x_yyyy_zzzz, to_yyy_xxx, to_yyy_xxxxx, to_yyy_xxxxy, to_yyy_xxxxz, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_zzz, to_yyyyy_xxx, to_yyyyy_xxxxx, to_yyyyy_xxxxy, to_yyyyy_xxxxz, to_yyyyy_xxxyy, to_yyyyy_xxxyz, to_yyyyy_xxxzz, to_yyyyy_xxy, to_yyyyy_xxyyy, to_yyyyy_xxyyz, to_yyyyy_xxyzz, to_yyyyy_xxz, to_yyyyy_xxzzz, to_yyyyy_xyy, to_yyyyy_xyyyy, to_yyyyy_xyyyz, to_yyyyy_xyyzz, to_yyyyy_xyz, to_yyyyy_xyzzz, to_yyyyy_xzz, to_yyyyy_xzzzz, to_yyyyy_yyy, to_yyyyy_yyz, to_yyyyy_yzz, to_yyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyy_xxxx[k] = 16.0 * to_yyy_xxx[k] - 8.0 * to_yyy_xxxxx[k] * tke_0 - 8.0 * to_yyyyy_xxx[k] * tbe_0 + 4.0 * to_yyyyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxxy[k] = 12.0 * to_yyy_xxy[k] - 8.0 * to_yyy_xxxxy[k] * tke_0 - 6.0 * to_yyyyy_xxy[k] * tbe_0 + 4.0 * to_yyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxxz[k] = 12.0 * to_yyy_xxz[k] - 8.0 * to_yyy_xxxxz[k] * tke_0 - 6.0 * to_yyyyy_xxz[k] * tbe_0 + 4.0 * to_yyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxyy[k] = 8.0 * to_yyy_xyy[k] - 8.0 * to_yyy_xxxyy[k] * tke_0 - 4.0 * to_yyyyy_xyy[k] * tbe_0 + 4.0 * to_yyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxyz[k] = 8.0 * to_yyy_xyz[k] - 8.0 * to_yyy_xxxyz[k] * tke_0 - 4.0 * to_yyyyy_xyz[k] * tbe_0 + 4.0 * to_yyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxzz[k] = 8.0 * to_yyy_xzz[k] - 8.0 * to_yyy_xxxzz[k] * tke_0 - 4.0 * to_yyyyy_xzz[k] * tbe_0 + 4.0 * to_yyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xyyy[k] = 4.0 * to_yyy_yyy[k] - 8.0 * to_yyy_xxyyy[k] * tke_0 - 2.0 * to_yyyyy_yyy[k] * tbe_0 + 4.0 * to_yyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xyyz[k] = 4.0 * to_yyy_yyz[k] - 8.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_yyyyy_yyz[k] * tbe_0 + 4.0 * to_yyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xyzz[k] = 4.0 * to_yyy_yzz[k] - 8.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_yyyyy_yzz[k] * tbe_0 + 4.0 * to_yyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xzzz[k] = 4.0 * to_yyy_zzz[k] - 8.0 * to_yyy_xxzzz[k] * tke_0 - 2.0 * to_yyyyy_zzz[k] * tbe_0 + 4.0 * to_yyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yyyy[k] = -8.0 * to_yyy_xyyyy[k] * tke_0 + 4.0 * to_yyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yyyz[k] = -8.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_yyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yyzz[k] = -8.0 * to_yyy_xyyzz[k] * tke_0 + 4.0 * to_yyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yzzz[k] = -8.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_yyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_zzzz[k] = -8.0 * to_yyy_xzzzz[k] * tke_0 + 4.0 * to_yyyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 840-855 components of targeted buffer : GG

        auto to_y_x_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 165);

        auto to_y_x_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 166);

        auto to_y_x_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 167);

        auto to_y_x_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 168);

        auto to_y_x_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 169);

        auto to_y_x_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 170);

        auto to_y_x_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 171);

        auto to_y_x_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 172);

        auto to_y_x_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 173);

        auto to_y_x_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 174);

        auto to_y_x_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 175);

        auto to_y_x_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 176);

        auto to_y_x_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 177);

        auto to_y_x_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 178);

        auto to_y_x_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_y_x_yyyz_xxxx, to_y_x_yyyz_xxxy, to_y_x_yyyz_xxxz, to_y_x_yyyz_xxyy, to_y_x_yyyz_xxyz, to_y_x_yyyz_xxzz, to_y_x_yyyz_xyyy, to_y_x_yyyz_xyyz, to_y_x_yyyz_xyzz, to_y_x_yyyz_xzzz, to_y_x_yyyz_yyyy, to_y_x_yyyz_yyyz, to_y_x_yyyz_yyzz, to_y_x_yyyz_yzzz, to_y_x_yyyz_zzzz, to_yyyyz_xxx, to_yyyyz_xxxxx, to_yyyyz_xxxxy, to_yyyyz_xxxxz, to_yyyyz_xxxyy, to_yyyyz_xxxyz, to_yyyyz_xxxzz, to_yyyyz_xxy, to_yyyyz_xxyyy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xxzzz, to_yyyyz_xyy, to_yyyyz_xyyyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_xzzzz, to_yyyyz_yyy, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_zzz, to_yyz_xxx, to_yyz_xxxxx, to_yyz_xxxxy, to_yyz_xxxxz, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyz_xxxx[k] = 12.0 * to_yyz_xxx[k] - 6.0 * to_yyz_xxxxx[k] * tke_0 - 8.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxxy[k] = 9.0 * to_yyz_xxy[k] - 6.0 * to_yyz_xxxxy[k] * tke_0 - 6.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxxz[k] = 9.0 * to_yyz_xxz[k] - 6.0 * to_yyz_xxxxz[k] * tke_0 - 6.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxyy[k] = 6.0 * to_yyz_xyy[k] - 6.0 * to_yyz_xxxyy[k] * tke_0 - 4.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxyz[k] = 6.0 * to_yyz_xyz[k] - 6.0 * to_yyz_xxxyz[k] * tke_0 - 4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxzz[k] = 6.0 * to_yyz_xzz[k] - 6.0 * to_yyz_xxxzz[k] * tke_0 - 4.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xyyy[k] = 3.0 * to_yyz_yyy[k] - 6.0 * to_yyz_xxyyy[k] * tke_0 - 2.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xyyz[k] = 3.0 * to_yyz_yyz[k] - 6.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xyzz[k] = 3.0 * to_yyz_yzz[k] - 6.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xzzz[k] = 3.0 * to_yyz_zzz[k] - 6.0 * to_yyz_xxzzz[k] * tke_0 - 2.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yyyy[k] = -6.0 * to_yyz_xyyyy[k] * tke_0 + 4.0 * to_yyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yyyz[k] = -6.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yyzz[k] = -6.0 * to_yyz_xyyzz[k] * tke_0 + 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yzzz[k] = -6.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_zzzz[k] = -6.0 * to_yyz_xzzzz[k] * tke_0 + 4.0 * to_yyyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 855-870 components of targeted buffer : GG

        auto to_y_x_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 180);

        auto to_y_x_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 181);

        auto to_y_x_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 182);

        auto to_y_x_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 183);

        auto to_y_x_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 184);

        auto to_y_x_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 185);

        auto to_y_x_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 186);

        auto to_y_x_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 187);

        auto to_y_x_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 188);

        auto to_y_x_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 189);

        auto to_y_x_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 190);

        auto to_y_x_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 191);

        auto to_y_x_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 192);

        auto to_y_x_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 193);

        auto to_y_x_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_y_x_yyzz_xxxx, to_y_x_yyzz_xxxy, to_y_x_yyzz_xxxz, to_y_x_yyzz_xxyy, to_y_x_yyzz_xxyz, to_y_x_yyzz_xxzz, to_y_x_yyzz_xyyy, to_y_x_yyzz_xyyz, to_y_x_yyzz_xyzz, to_y_x_yyzz_xzzz, to_y_x_yyzz_yyyy, to_y_x_yyzz_yyyz, to_y_x_yyzz_yyzz, to_y_x_yyzz_yzzz, to_y_x_yyzz_zzzz, to_yyyzz_xxx, to_yyyzz_xxxxx, to_yyyzz_xxxxy, to_yyyzz_xxxxz, to_yyyzz_xxxyy, to_yyyzz_xxxyz, to_yyyzz_xxxzz, to_yyyzz_xxy, to_yyyzz_xxyyy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xxzzz, to_yyyzz_xyy, to_yyyzz_xyyyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_xzzzz, to_yyyzz_yyy, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_zzz, to_yzz_xxx, to_yzz_xxxxx, to_yzz_xxxxy, to_yzz_xxxxz, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyzz_xxxx[k] = 8.0 * to_yzz_xxx[k] - 4.0 * to_yzz_xxxxx[k] * tke_0 - 8.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxxy[k] = 6.0 * to_yzz_xxy[k] - 4.0 * to_yzz_xxxxy[k] * tke_0 - 6.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxxz[k] = 6.0 * to_yzz_xxz[k] - 4.0 * to_yzz_xxxxz[k] * tke_0 - 6.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxyy[k] = 4.0 * to_yzz_xyy[k] - 4.0 * to_yzz_xxxyy[k] * tke_0 - 4.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxyz[k] = 4.0 * to_yzz_xyz[k] - 4.0 * to_yzz_xxxyz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxzz[k] = 4.0 * to_yzz_xzz[k] - 4.0 * to_yzz_xxxzz[k] * tke_0 - 4.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xyyy[k] = 2.0 * to_yzz_yyy[k] - 4.0 * to_yzz_xxyyy[k] * tke_0 - 2.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xyyz[k] = 2.0 * to_yzz_yyz[k] - 4.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xyzz[k] = 2.0 * to_yzz_yzz[k] - 4.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xzzz[k] = 2.0 * to_yzz_zzz[k] - 4.0 * to_yzz_xxzzz[k] * tke_0 - 2.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yyyy[k] = -4.0 * to_yzz_xyyyy[k] * tke_0 + 4.0 * to_yyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yyyz[k] = -4.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yyzz[k] = -4.0 * to_yzz_xyyzz[k] * tke_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yzzz[k] = -4.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_zzzz[k] = -4.0 * to_yzz_xzzzz[k] * tke_0 + 4.0 * to_yyyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 870-885 components of targeted buffer : GG

        auto to_y_x_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 195);

        auto to_y_x_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 196);

        auto to_y_x_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 197);

        auto to_y_x_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 198);

        auto to_y_x_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 199);

        auto to_y_x_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 200);

        auto to_y_x_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 201);

        auto to_y_x_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 202);

        auto to_y_x_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 203);

        auto to_y_x_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 204);

        auto to_y_x_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 205);

        auto to_y_x_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 206);

        auto to_y_x_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 207);

        auto to_y_x_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 208);

        auto to_y_x_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_y_x_yzzz_xxxx, to_y_x_yzzz_xxxy, to_y_x_yzzz_xxxz, to_y_x_yzzz_xxyy, to_y_x_yzzz_xxyz, to_y_x_yzzz_xxzz, to_y_x_yzzz_xyyy, to_y_x_yzzz_xyyz, to_y_x_yzzz_xyzz, to_y_x_yzzz_xzzz, to_y_x_yzzz_yyyy, to_y_x_yzzz_yyyz, to_y_x_yzzz_yyzz, to_y_x_yzzz_yzzz, to_y_x_yzzz_zzzz, to_yyzzz_xxx, to_yyzzz_xxxxx, to_yyzzz_xxxxy, to_yyzzz_xxxxz, to_yyzzz_xxxyy, to_yyzzz_xxxyz, to_yyzzz_xxxzz, to_yyzzz_xxy, to_yyzzz_xxyyy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xxzzz, to_yyzzz_xyy, to_yyzzz_xyyyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_xzzzz, to_yyzzz_yyy, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_zzz, to_zzz_xxx, to_zzz_xxxxx, to_zzz_xxxxy, to_zzz_xxxxz, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzzz_xxxx[k] = 4.0 * to_zzz_xxx[k] - 2.0 * to_zzz_xxxxx[k] * tke_0 - 8.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxxy[k] = 3.0 * to_zzz_xxy[k] - 2.0 * to_zzz_xxxxy[k] * tke_0 - 6.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxxz[k] = 3.0 * to_zzz_xxz[k] - 2.0 * to_zzz_xxxxz[k] * tke_0 - 6.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxyy[k] = 2.0 * to_zzz_xyy[k] - 2.0 * to_zzz_xxxyy[k] * tke_0 - 4.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxyz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xxxyz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxzz[k] = 2.0 * to_zzz_xzz[k] - 2.0 * to_zzz_xxxzz[k] * tke_0 - 4.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xyyy[k] = to_zzz_yyy[k] - 2.0 * to_zzz_xxyyy[k] * tke_0 - 2.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xyyz[k] = to_zzz_yyz[k] - 2.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xyzz[k] = to_zzz_yzz[k] - 2.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xzzz[k] = to_zzz_zzz[k] - 2.0 * to_zzz_xxzzz[k] * tke_0 - 2.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yyyy[k] = -2.0 * to_zzz_xyyyy[k] * tke_0 + 4.0 * to_yyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yyyz[k] = -2.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yyzz[k] = -2.0 * to_zzz_xyyzz[k] * tke_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yzzz[k] = -2.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_zzzz[k] = -2.0 * to_zzz_xzzzz[k] * tke_0 + 4.0 * to_yyzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 885-900 components of targeted buffer : GG

        auto to_y_x_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 210);

        auto to_y_x_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 211);

        auto to_y_x_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 212);

        auto to_y_x_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 213);

        auto to_y_x_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 214);

        auto to_y_x_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 215);

        auto to_y_x_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 216);

        auto to_y_x_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 217);

        auto to_y_x_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 218);

        auto to_y_x_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 219);

        auto to_y_x_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 220);

        auto to_y_x_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 221);

        auto to_y_x_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 222);

        auto to_y_x_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 223);

        auto to_y_x_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 3 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_y_x_zzzz_xxxx, to_y_x_zzzz_xxxy, to_y_x_zzzz_xxxz, to_y_x_zzzz_xxyy, to_y_x_zzzz_xxyz, to_y_x_zzzz_xxzz, to_y_x_zzzz_xyyy, to_y_x_zzzz_xyyz, to_y_x_zzzz_xyzz, to_y_x_zzzz_xzzz, to_y_x_zzzz_yyyy, to_y_x_zzzz_yyyz, to_y_x_zzzz_yyzz, to_y_x_zzzz_yzzz, to_y_x_zzzz_zzzz, to_yzzzz_xxx, to_yzzzz_xxxxx, to_yzzzz_xxxxy, to_yzzzz_xxxxz, to_yzzzz_xxxyy, to_yzzzz_xxxyz, to_yzzzz_xxxzz, to_yzzzz_xxy, to_yzzzz_xxyyy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xxzzz, to_yzzzz_xyy, to_yzzzz_xyyyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_xzzzz, to_yzzzz_yyy, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzzz_xxxx[k] = -8.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxxy[k] = -6.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxxz[k] = -6.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxyy[k] = -4.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxyz[k] = -4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxzz[k] = -4.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xyyy[k] = -2.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xyyz[k] = -2.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xyzz[k] = -2.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xzzz[k] = -2.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yyyy[k] = 4.0 * to_yzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yyyz[k] = 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yyzz[k] = 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yzzz[k] = 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_zzzz[k] = 4.0 * to_yzzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 900-915 components of targeted buffer : GG

        auto to_y_y_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 0);

        auto to_y_y_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 1);

        auto to_y_y_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 2);

        auto to_y_y_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 3);

        auto to_y_y_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 4);

        auto to_y_y_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 5);

        auto to_y_y_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 6);

        auto to_y_y_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 7);

        auto to_y_y_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 8);

        auto to_y_y_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 9);

        auto to_y_y_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 10);

        auto to_y_y_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 11);

        auto to_y_y_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 12);

        auto to_y_y_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 13);

        auto to_y_y_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxy_xxx, to_xxxxy_xxxxy, to_xxxxy_xxxyy, to_xxxxy_xxxyz, to_xxxxy_xxy, to_xxxxy_xxyyy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xyy, to_xxxxy_xyyyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_yyy, to_xxxxy_yyyyy, to_xxxxy_yyyyz, to_xxxxy_yyyzz, to_xxxxy_yyz, to_xxxxy_yyzzz, to_xxxxy_yzz, to_xxxxy_yzzzz, to_xxxxy_zzz, to_y_y_xxxx_xxxx, to_y_y_xxxx_xxxy, to_y_y_xxxx_xxxz, to_y_y_xxxx_xxyy, to_y_y_xxxx_xxyz, to_y_y_xxxx_xxzz, to_y_y_xxxx_xyyy, to_y_y_xxxx_xyyz, to_y_y_xxxx_xyzz, to_y_y_xxxx_xzzz, to_y_y_xxxx_yyyy, to_y_y_xxxx_yyyz, to_y_y_xxxx_yyzz, to_y_y_xxxx_yzzz, to_y_y_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxx_xxxx[k] = 4.0 * to_xxxxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxxy[k] = -2.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxxz[k] = 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxyy[k] = -4.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxyz[k] = -2.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxzz[k] = 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xyyy[k] = -6.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xyyz[k] = -4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xyzz[k] = -2.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xzzz[k] = 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yyyy[k] = -8.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yyyz[k] = -6.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yyzz[k] = -4.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yzzz[k] = -2.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_zzzz[k] = 4.0 * to_xxxxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 915-930 components of targeted buffer : GG

        auto to_y_y_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 15);

        auto to_y_y_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 16);

        auto to_y_y_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 17);

        auto to_y_y_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 18);

        auto to_y_y_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 19);

        auto to_y_y_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 20);

        auto to_y_y_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 21);

        auto to_y_y_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 22);

        auto to_y_y_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 23);

        auto to_y_y_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 24);

        auto to_y_y_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 25);

        auto to_y_y_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 26);

        auto to_y_y_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 27);

        auto to_y_y_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 28);

        auto to_y_y_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxy, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_yyy, to_xxx_yyyyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxxyy_xxx, to_xxxyy_xxxxy, to_xxxyy_xxxyy, to_xxxyy_xxxyz, to_xxxyy_xxy, to_xxxyy_xxyyy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xyy, to_xxxyy_xyyyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_yyy, to_xxxyy_yyyyy, to_xxxyy_yyyyz, to_xxxyy_yyyzz, to_xxxyy_yyz, to_xxxyy_yyzzz, to_xxxyy_yzz, to_xxxyy_yzzzz, to_xxxyy_zzz, to_y_y_xxxy_xxxx, to_y_y_xxxy_xxxy, to_y_y_xxxy_xxxz, to_y_y_xxxy_xxyy, to_y_y_xxxy_xxyz, to_y_y_xxxy_xxzz, to_y_y_xxxy_xyyy, to_y_y_xxxy_xyyz, to_y_y_xxxy_xyzz, to_y_y_xxxy_xzzz, to_y_y_xxxy_yyyy, to_y_y_xxxy_yyyz, to_y_y_xxxy_yyzz, to_y_y_xxxy_yzzz, to_y_y_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxy_xxxx[k] = -2.0 * to_xxx_xxxxy[k] * tke_0 + 4.0 * to_xxxyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxxy[k] = to_xxx_xxx[k] - 2.0 * to_xxx_xxxyy[k] * tke_0 - 2.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxxz[k] = -2.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxyy[k] = 2.0 * to_xxx_xxy[k] - 2.0 * to_xxx_xxyyy[k] * tke_0 - 4.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxyz[k] = to_xxx_xxz[k] - 2.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxzz[k] = -2.0 * to_xxx_xxyzz[k] * tke_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xyyy[k] = 3.0 * to_xxx_xyy[k] - 2.0 * to_xxx_xyyyy[k] * tke_0 - 6.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xyyz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xyyyz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xyzz[k] = to_xxx_xzz[k] - 2.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xzzz[k] = -2.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yyyy[k] = 4.0 * to_xxx_yyy[k] - 2.0 * to_xxx_yyyyy[k] * tke_0 - 8.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yyyz[k] = 3.0 * to_xxx_yyz[k] - 2.0 * to_xxx_yyyyz[k] * tke_0 - 6.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yyzz[k] = 2.0 * to_xxx_yzz[k] - 2.0 * to_xxx_yyyzz[k] * tke_0 - 4.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yzzz[k] = to_xxx_zzz[k] - 2.0 * to_xxx_yyzzz[k] * tke_0 - 2.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_zzzz[k] = -2.0 * to_xxx_yzzzz[k] * tke_0 + 4.0 * to_xxxyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 930-945 components of targeted buffer : GG

        auto to_y_y_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 30);

        auto to_y_y_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 31);

        auto to_y_y_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 32);

        auto to_y_y_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 33);

        auto to_y_y_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 34);

        auto to_y_y_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 35);

        auto to_y_y_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 36);

        auto to_y_y_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 37);

        auto to_y_y_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 38);

        auto to_y_y_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 39);

        auto to_y_y_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 40);

        auto to_y_y_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 41);

        auto to_y_y_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 42);

        auto to_y_y_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 43);

        auto to_y_y_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxy, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_yyy, to_xxxyz_yyyyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_y_y_xxxz_xxxx, to_y_y_xxxz_xxxy, to_y_y_xxxz_xxxz, to_y_y_xxxz_xxyy, to_y_y_xxxz_xxyz, to_y_y_xxxz_xxzz, to_y_y_xxxz_xyyy, to_y_y_xxxz_xyyz, to_y_y_xxxz_xyzz, to_y_y_xxxz_xzzz, to_y_y_xxxz_yyyy, to_y_y_xxxz_yyyz, to_y_y_xxxz_yyzz, to_y_y_xxxz_yzzz, to_y_y_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxz_xxxx[k] = 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxxy[k] = -2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxxz[k] = 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxyy[k] = -4.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxyz[k] = -2.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxzz[k] = 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xyyy[k] = -6.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xyyz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xyzz[k] = -2.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xzzz[k] = 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yyyy[k] = -8.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yyyz[k] = -6.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yyzz[k] = -4.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yzzz[k] = -2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_zzzz[k] = 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 945-960 components of targeted buffer : GG

        auto to_y_y_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 45);

        auto to_y_y_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 46);

        auto to_y_y_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 47);

        auto to_y_y_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 48);

        auto to_y_y_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 49);

        auto to_y_y_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 50);

        auto to_y_y_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 51);

        auto to_y_y_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 52);

        auto to_y_y_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 53);

        auto to_y_y_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 54);

        auto to_y_y_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 55);

        auto to_y_y_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 56);

        auto to_y_y_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 57);

        auto to_y_y_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 58);

        auto to_y_y_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxy, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_yyy, to_xxy_yyyyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxyyy_xxx, to_xxyyy_xxxxy, to_xxyyy_xxxyy, to_xxyyy_xxxyz, to_xxyyy_xxy, to_xxyyy_xxyyy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xyy, to_xxyyy_xyyyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_yyy, to_xxyyy_yyyyy, to_xxyyy_yyyyz, to_xxyyy_yyyzz, to_xxyyy_yyz, to_xxyyy_yyzzz, to_xxyyy_yzz, to_xxyyy_yzzzz, to_xxyyy_zzz, to_y_y_xxyy_xxxx, to_y_y_xxyy_xxxy, to_y_y_xxyy_xxxz, to_y_y_xxyy_xxyy, to_y_y_xxyy_xxyz, to_y_y_xxyy_xxzz, to_y_y_xxyy_xyyy, to_y_y_xxyy_xyyz, to_y_y_xxyy_xyzz, to_y_y_xxyy_xzzz, to_y_y_xxyy_yyyy, to_y_y_xxyy_yyyz, to_y_y_xxyy_yyzz, to_y_y_xxyy_yzzz, to_y_y_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyy_xxxx[k] = -4.0 * to_xxy_xxxxy[k] * tke_0 + 4.0 * to_xxyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxxy[k] = 2.0 * to_xxy_xxx[k] - 4.0 * to_xxy_xxxyy[k] * tke_0 - 2.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxxz[k] = -4.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxyy[k] = 4.0 * to_xxy_xxy[k] - 4.0 * to_xxy_xxyyy[k] * tke_0 - 4.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxyz[k] = 2.0 * to_xxy_xxz[k] - 4.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxzz[k] = -4.0 * to_xxy_xxyzz[k] * tke_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xyyy[k] = 6.0 * to_xxy_xyy[k] - 4.0 * to_xxy_xyyyy[k] * tke_0 - 6.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xyyz[k] = 4.0 * to_xxy_xyz[k] - 4.0 * to_xxy_xyyyz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xyzz[k] = 2.0 * to_xxy_xzz[k] - 4.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xzzz[k] = -4.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yyyy[k] = 8.0 * to_xxy_yyy[k] - 4.0 * to_xxy_yyyyy[k] * tke_0 - 8.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yyyz[k] = 6.0 * to_xxy_yyz[k] - 4.0 * to_xxy_yyyyz[k] * tke_0 - 6.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yyzz[k] = 4.0 * to_xxy_yzz[k] - 4.0 * to_xxy_yyyzz[k] * tke_0 - 4.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yzzz[k] = 2.0 * to_xxy_zzz[k] - 4.0 * to_xxy_yyzzz[k] * tke_0 - 2.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_zzzz[k] = -4.0 * to_xxy_yzzzz[k] * tke_0 + 4.0 * to_xxyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 960-975 components of targeted buffer : GG

        auto to_y_y_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 60);

        auto to_y_y_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 61);

        auto to_y_y_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 62);

        auto to_y_y_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 63);

        auto to_y_y_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 64);

        auto to_y_y_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 65);

        auto to_y_y_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 66);

        auto to_y_y_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 67);

        auto to_y_y_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 68);

        auto to_y_y_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 69);

        auto to_y_y_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 70);

        auto to_y_y_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 71);

        auto to_y_y_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 72);

        auto to_y_y_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 73);

        auto to_y_y_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxy, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_yyy, to_xxyyz_yyyyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_xxz_xxx, to_xxz_xxxxy, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_yyy, to_xxz_yyyyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_y_y_xxyz_xxxx, to_y_y_xxyz_xxxy, to_y_y_xxyz_xxxz, to_y_y_xxyz_xxyy, to_y_y_xxyz_xxyz, to_y_y_xxyz_xxzz, to_y_y_xxyz_xyyy, to_y_y_xxyz_xyyz, to_y_y_xxyz_xyzz, to_y_y_xxyz_xzzz, to_y_y_xxyz_yyyy, to_y_y_xxyz_yyyz, to_y_y_xxyz_yyzz, to_y_y_xxyz_yzzz, to_y_y_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyz_xxxx[k] = -2.0 * to_xxz_xxxxy[k] * tke_0 + 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxxy[k] = to_xxz_xxx[k] - 2.0 * to_xxz_xxxyy[k] * tke_0 - 2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxxz[k] = -2.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxyy[k] = 2.0 * to_xxz_xxy[k] - 2.0 * to_xxz_xxyyy[k] * tke_0 - 4.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxyz[k] = to_xxz_xxz[k] - 2.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxzz[k] = -2.0 * to_xxz_xxyzz[k] * tke_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xyyy[k] = 3.0 * to_xxz_xyy[k] - 2.0 * to_xxz_xyyyy[k] * tke_0 - 6.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xyyz[k] = 2.0 * to_xxz_xyz[k] - 2.0 * to_xxz_xyyyz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xyzz[k] = to_xxz_xzz[k] - 2.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xzzz[k] = -2.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yyyy[k] = 4.0 * to_xxz_yyy[k] - 2.0 * to_xxz_yyyyy[k] * tke_0 - 8.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yyyz[k] = 3.0 * to_xxz_yyz[k] - 2.0 * to_xxz_yyyyz[k] * tke_0 - 6.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yyzz[k] = 2.0 * to_xxz_yzz[k] - 2.0 * to_xxz_yyyzz[k] * tke_0 - 4.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yzzz[k] = to_xxz_zzz[k] - 2.0 * to_xxz_yyzzz[k] * tke_0 - 2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_zzzz[k] = -2.0 * to_xxz_yzzzz[k] * tke_0 + 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 975-990 components of targeted buffer : GG

        auto to_y_y_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 75);

        auto to_y_y_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 76);

        auto to_y_y_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 77);

        auto to_y_y_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 78);

        auto to_y_y_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 79);

        auto to_y_y_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 80);

        auto to_y_y_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 81);

        auto to_y_y_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 82);

        auto to_y_y_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 83);

        auto to_y_y_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 84);

        auto to_y_y_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 85);

        auto to_y_y_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 86);

        auto to_y_y_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 87);

        auto to_y_y_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 88);

        auto to_y_y_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxyzz_xxx, to_xxyzz_xxxxy, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_yyy, to_xxyzz_yyyyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_y_y_xxzz_xxxx, to_y_y_xxzz_xxxy, to_y_y_xxzz_xxxz, to_y_y_xxzz_xxyy, to_y_y_xxzz_xxyz, to_y_y_xxzz_xxzz, to_y_y_xxzz_xyyy, to_y_y_xxzz_xyyz, to_y_y_xxzz_xyzz, to_y_y_xxzz_xzzz, to_y_y_xxzz_yyyy, to_y_y_xxzz_yyyz, to_y_y_xxzz_yyzz, to_y_y_xxzz_yzzz, to_y_y_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxzz_xxxx[k] = 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxxy[k] = -2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxxz[k] = 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxyy[k] = -4.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxyz[k] = -2.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxzz[k] = 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xyyy[k] = -6.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xyyz[k] = -4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xyzz[k] = -2.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xzzz[k] = 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yyyy[k] = -8.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yyyz[k] = -6.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yyzz[k] = -4.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yzzz[k] = -2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_zzzz[k] = 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 990-1005 components of targeted buffer : GG

        auto to_y_y_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 90);

        auto to_y_y_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 91);

        auto to_y_y_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 92);

        auto to_y_y_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 93);

        auto to_y_y_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 94);

        auto to_y_y_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 95);

        auto to_y_y_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 96);

        auto to_y_y_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 97);

        auto to_y_y_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 98);

        auto to_y_y_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 99);

        auto to_y_y_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 100);

        auto to_y_y_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 101);

        auto to_y_y_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 102);

        auto to_y_y_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 103);

        auto to_y_y_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxy, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_yyy, to_xyy_yyyyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyyyy_xxx, to_xyyyy_xxxxy, to_xyyyy_xxxyy, to_xyyyy_xxxyz, to_xyyyy_xxy, to_xyyyy_xxyyy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xyy, to_xyyyy_xyyyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_yyy, to_xyyyy_yyyyy, to_xyyyy_yyyyz, to_xyyyy_yyyzz, to_xyyyy_yyz, to_xyyyy_yyzzz, to_xyyyy_yzz, to_xyyyy_yzzzz, to_xyyyy_zzz, to_y_y_xyyy_xxxx, to_y_y_xyyy_xxxy, to_y_y_xyyy_xxxz, to_y_y_xyyy_xxyy, to_y_y_xyyy_xxyz, to_y_y_xyyy_xxzz, to_y_y_xyyy_xyyy, to_y_y_xyyy_xyyz, to_y_y_xyyy_xyzz, to_y_y_xyyy_xzzz, to_y_y_xyyy_yyyy, to_y_y_xyyy_yyyz, to_y_y_xyyy_yyzz, to_y_y_xyyy_yzzz, to_y_y_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyy_xxxx[k] = -6.0 * to_xyy_xxxxy[k] * tke_0 + 4.0 * to_xyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxxy[k] = 3.0 * to_xyy_xxx[k] - 6.0 * to_xyy_xxxyy[k] * tke_0 - 2.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxxz[k] = -6.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxyy[k] = 6.0 * to_xyy_xxy[k] - 6.0 * to_xyy_xxyyy[k] * tke_0 - 4.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxyz[k] = 3.0 * to_xyy_xxz[k] - 6.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxzz[k] = -6.0 * to_xyy_xxyzz[k] * tke_0 + 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xyyy[k] = 9.0 * to_xyy_xyy[k] - 6.0 * to_xyy_xyyyy[k] * tke_0 - 6.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xyyz[k] = 6.0 * to_xyy_xyz[k] - 6.0 * to_xyy_xyyyz[k] * tke_0 - 4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xyzz[k] = 3.0 * to_xyy_xzz[k] - 6.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xzzz[k] = -6.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yyyy[k] = 12.0 * to_xyy_yyy[k] - 6.0 * to_xyy_yyyyy[k] * tke_0 - 8.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yyyz[k] = 9.0 * to_xyy_yyz[k] - 6.0 * to_xyy_yyyyz[k] * tke_0 - 6.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yyzz[k] = 6.0 * to_xyy_yzz[k] - 6.0 * to_xyy_yyyzz[k] * tke_0 - 4.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yzzz[k] = 3.0 * to_xyy_zzz[k] - 6.0 * to_xyy_yyzzz[k] * tke_0 - 2.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_zzzz[k] = -6.0 * to_xyy_yzzzz[k] * tke_0 + 4.0 * to_xyyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1005-1020 components of targeted buffer : GG

        auto to_y_y_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 105);

        auto to_y_y_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 106);

        auto to_y_y_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 107);

        auto to_y_y_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 108);

        auto to_y_y_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 109);

        auto to_y_y_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 110);

        auto to_y_y_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 111);

        auto to_y_y_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 112);

        auto to_y_y_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 113);

        auto to_y_y_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 114);

        auto to_y_y_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 115);

        auto to_y_y_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 116);

        auto to_y_y_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 117);

        auto to_y_y_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 118);

        auto to_y_y_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxy, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_yyy, to_xyyyz_yyyyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, to_xyz_xxx, to_xyz_xxxxy, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_yyy, to_xyz_yyyyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_y_y_xyyz_xxxx, to_y_y_xyyz_xxxy, to_y_y_xyyz_xxxz, to_y_y_xyyz_xxyy, to_y_y_xyyz_xxyz, to_y_y_xyyz_xxzz, to_y_y_xyyz_xyyy, to_y_y_xyyz_xyyz, to_y_y_xyyz_xyzz, to_y_y_xyyz_xzzz, to_y_y_xyyz_yyyy, to_y_y_xyyz_yyyz, to_y_y_xyyz_yyzz, to_y_y_xyyz_yzzz, to_y_y_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyz_xxxx[k] = -4.0 * to_xyz_xxxxy[k] * tke_0 + 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxxy[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxxz[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxyy[k] = 4.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 4.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxyz[k] = 2.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxzz[k] = -4.0 * to_xyz_xxyzz[k] * tke_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xyyy[k] = 6.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyyy[k] * tke_0 - 6.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xyyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyyyz[k] * tke_0 - 4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xyzz[k] = 2.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yyyy[k] = 8.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyyy[k] * tke_0 - 8.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yyyz[k] = 6.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyyyz[k] * tke_0 - 6.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yyzz[k] = 4.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 4.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_zzzz[k] = -4.0 * to_xyz_yzzzz[k] * tke_0 + 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1020-1035 components of targeted buffer : GG

        auto to_y_y_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 120);

        auto to_y_y_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 121);

        auto to_y_y_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 122);

        auto to_y_y_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 123);

        auto to_y_y_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 124);

        auto to_y_y_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 125);

        auto to_y_y_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 126);

        auto to_y_y_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 127);

        auto to_y_y_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 128);

        auto to_y_y_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 129);

        auto to_y_y_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 130);

        auto to_y_y_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 131);

        auto to_y_y_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 132);

        auto to_y_y_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 133);

        auto to_y_y_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyyzz_xxx, to_xyyzz_xxxxy, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_yyy, to_xyyzz_yyyyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, to_xzz_xxx, to_xzz_xxxxy, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_yyy, to_xzz_yyyyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_y_y_xyzz_xxxx, to_y_y_xyzz_xxxy, to_y_y_xyzz_xxxz, to_y_y_xyzz_xxyy, to_y_y_xyzz_xxyz, to_y_y_xyzz_xxzz, to_y_y_xyzz_xyyy, to_y_y_xyzz_xyyz, to_y_y_xyzz_xyzz, to_y_y_xyzz_xzzz, to_y_y_xyzz_yyyy, to_y_y_xyzz_yyyz, to_y_y_xyzz_yyzz, to_y_y_xyzz_yzzz, to_y_y_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyzz_xxxx[k] = -2.0 * to_xzz_xxxxy[k] * tke_0 + 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxxy[k] = to_xzz_xxx[k] - 2.0 * to_xzz_xxxyy[k] * tke_0 - 2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxxz[k] = -2.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxyy[k] = 2.0 * to_xzz_xxy[k] - 2.0 * to_xzz_xxyyy[k] * tke_0 - 4.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxyz[k] = to_xzz_xxz[k] - 2.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxzz[k] = -2.0 * to_xzz_xxyzz[k] * tke_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xyyy[k] = 3.0 * to_xzz_xyy[k] - 2.0 * to_xzz_xyyyy[k] * tke_0 - 6.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xyyz[k] = 2.0 * to_xzz_xyz[k] - 2.0 * to_xzz_xyyyz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xyzz[k] = to_xzz_xzz[k] - 2.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xzzz[k] = -2.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yyyy[k] = 4.0 * to_xzz_yyy[k] - 2.0 * to_xzz_yyyyy[k] * tke_0 - 8.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yyyz[k] = 3.0 * to_xzz_yyz[k] - 2.0 * to_xzz_yyyyz[k] * tke_0 - 6.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yyzz[k] = 2.0 * to_xzz_yzz[k] - 2.0 * to_xzz_yyyzz[k] * tke_0 - 4.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yzzz[k] = to_xzz_zzz[k] - 2.0 * to_xzz_yyzzz[k] * tke_0 - 2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_zzzz[k] = -2.0 * to_xzz_yzzzz[k] * tke_0 + 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1035-1050 components of targeted buffer : GG

        auto to_y_y_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 135);

        auto to_y_y_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 136);

        auto to_y_y_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 137);

        auto to_y_y_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 138);

        auto to_y_y_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 139);

        auto to_y_y_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 140);

        auto to_y_y_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 141);

        auto to_y_y_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 142);

        auto to_y_y_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 143);

        auto to_y_y_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 144);

        auto to_y_y_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 145);

        auto to_y_y_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 146);

        auto to_y_y_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 147);

        auto to_y_y_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 148);

        auto to_y_y_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xyzzz_xxx, to_xyzzz_xxxxy, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_yyy, to_xyzzz_yyyyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, to_y_y_xzzz_xxxx, to_y_y_xzzz_xxxy, to_y_y_xzzz_xxxz, to_y_y_xzzz_xxyy, to_y_y_xzzz_xxyz, to_y_y_xzzz_xxzz, to_y_y_xzzz_xyyy, to_y_y_xzzz_xyyz, to_y_y_xzzz_xyzz, to_y_y_xzzz_xzzz, to_y_y_xzzz_yyyy, to_y_y_xzzz_yyyz, to_y_y_xzzz_yyzz, to_y_y_xzzz_yzzz, to_y_y_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzzz_xxxx[k] = 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxxy[k] = -2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxxz[k] = 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxyy[k] = -4.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxyz[k] = -2.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxzz[k] = 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xyyy[k] = -6.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xyyz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xyzz[k] = -2.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xzzz[k] = 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yyyy[k] = -8.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yyyz[k] = -6.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yyzz[k] = -4.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yzzz[k] = -2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_zzzz[k] = 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1050-1065 components of targeted buffer : GG

        auto to_y_y_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 150);

        auto to_y_y_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 151);

        auto to_y_y_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 152);

        auto to_y_y_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 153);

        auto to_y_y_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 154);

        auto to_y_y_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 155);

        auto to_y_y_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 156);

        auto to_y_y_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 157);

        auto to_y_y_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 158);

        auto to_y_y_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 159);

        auto to_y_y_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 160);

        auto to_y_y_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 161);

        auto to_y_y_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 162);

        auto to_y_y_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 163);

        auto to_y_y_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_y_y_yyyy_xxxx, to_y_y_yyyy_xxxy, to_y_y_yyyy_xxxz, to_y_y_yyyy_xxyy, to_y_y_yyyy_xxyz, to_y_y_yyyy_xxzz, to_y_y_yyyy_xyyy, to_y_y_yyyy_xyyz, to_y_y_yyyy_xyzz, to_y_y_yyyy_xzzz, to_y_y_yyyy_yyyy, to_y_y_yyyy_yyyz, to_y_y_yyyy_yyzz, to_y_y_yyyy_yzzz, to_y_y_yyyy_zzzz, to_yyy_xxx, to_yyy_xxxxy, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_yyy, to_yyy_yyyyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyyyy_xxx, to_yyyyy_xxxxy, to_yyyyy_xxxyy, to_yyyyy_xxxyz, to_yyyyy_xxy, to_yyyyy_xxyyy, to_yyyyy_xxyyz, to_yyyyy_xxyzz, to_yyyyy_xxz, to_yyyyy_xyy, to_yyyyy_xyyyy, to_yyyyy_xyyyz, to_yyyyy_xyyzz, to_yyyyy_xyz, to_yyyyy_xyzzz, to_yyyyy_xzz, to_yyyyy_yyy, to_yyyyy_yyyyy, to_yyyyy_yyyyz, to_yyyyy_yyyzz, to_yyyyy_yyz, to_yyyyy_yyzzz, to_yyyyy_yzz, to_yyyyy_yzzzz, to_yyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyy_xxxx[k] = -8.0 * to_yyy_xxxxy[k] * tke_0 + 4.0 * to_yyyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxxy[k] = 4.0 * to_yyy_xxx[k] - 8.0 * to_yyy_xxxyy[k] * tke_0 - 2.0 * to_yyyyy_xxx[k] * tbe_0 + 4.0 * to_yyyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxxz[k] = -8.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_yyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxyy[k] = 8.0 * to_yyy_xxy[k] - 8.0 * to_yyy_xxyyy[k] * tke_0 - 4.0 * to_yyyyy_xxy[k] * tbe_0 + 4.0 * to_yyyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxyz[k] = 4.0 * to_yyy_xxz[k] - 8.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_yyyyy_xxz[k] * tbe_0 + 4.0 * to_yyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxzz[k] = -8.0 * to_yyy_xxyzz[k] * tke_0 + 4.0 * to_yyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xyyy[k] = 12.0 * to_yyy_xyy[k] - 8.0 * to_yyy_xyyyy[k] * tke_0 - 6.0 * to_yyyyy_xyy[k] * tbe_0 + 4.0 * to_yyyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xyyz[k] = 8.0 * to_yyy_xyz[k] - 8.0 * to_yyy_xyyyz[k] * tke_0 - 4.0 * to_yyyyy_xyz[k] * tbe_0 + 4.0 * to_yyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xyzz[k] = 4.0 * to_yyy_xzz[k] - 8.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_yyyyy_xzz[k] * tbe_0 + 4.0 * to_yyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xzzz[k] = -8.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_yyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yyyy[k] = 16.0 * to_yyy_yyy[k] - 8.0 * to_yyy_yyyyy[k] * tke_0 - 8.0 * to_yyyyy_yyy[k] * tbe_0 + 4.0 * to_yyyyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yyyz[k] = 12.0 * to_yyy_yyz[k] - 8.0 * to_yyy_yyyyz[k] * tke_0 - 6.0 * to_yyyyy_yyz[k] * tbe_0 + 4.0 * to_yyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yyzz[k] = 8.0 * to_yyy_yzz[k] - 8.0 * to_yyy_yyyzz[k] * tke_0 - 4.0 * to_yyyyy_yzz[k] * tbe_0 + 4.0 * to_yyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yzzz[k] = 4.0 * to_yyy_zzz[k] - 8.0 * to_yyy_yyzzz[k] * tke_0 - 2.0 * to_yyyyy_zzz[k] * tbe_0 + 4.0 * to_yyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_zzzz[k] = -8.0 * to_yyy_yzzzz[k] * tke_0 + 4.0 * to_yyyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1065-1080 components of targeted buffer : GG

        auto to_y_y_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 165);

        auto to_y_y_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 166);

        auto to_y_y_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 167);

        auto to_y_y_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 168);

        auto to_y_y_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 169);

        auto to_y_y_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 170);

        auto to_y_y_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 171);

        auto to_y_y_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 172);

        auto to_y_y_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 173);

        auto to_y_y_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 174);

        auto to_y_y_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 175);

        auto to_y_y_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 176);

        auto to_y_y_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 177);

        auto to_y_y_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 178);

        auto to_y_y_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_y_y_yyyz_xxxx, to_y_y_yyyz_xxxy, to_y_y_yyyz_xxxz, to_y_y_yyyz_xxyy, to_y_y_yyyz_xxyz, to_y_y_yyyz_xxzz, to_y_y_yyyz_xyyy, to_y_y_yyyz_xyyz, to_y_y_yyyz_xyzz, to_y_y_yyyz_xzzz, to_y_y_yyyz_yyyy, to_y_y_yyyz_yyyz, to_y_y_yyyz_yyzz, to_y_y_yyyz_yzzz, to_y_y_yyyz_zzzz, to_yyyyz_xxx, to_yyyyz_xxxxy, to_yyyyz_xxxyy, to_yyyyz_xxxyz, to_yyyyz_xxy, to_yyyyz_xxyyy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xyy, to_yyyyz_xyyyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_yyy, to_yyyyz_yyyyy, to_yyyyz_yyyyz, to_yyyyz_yyyzz, to_yyyyz_yyz, to_yyyyz_yyzzz, to_yyyyz_yzz, to_yyyyz_yzzzz, to_yyyyz_zzz, to_yyz_xxx, to_yyz_xxxxy, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_yyy, to_yyz_yyyyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyz_xxxx[k] = -6.0 * to_yyz_xxxxy[k] * tke_0 + 4.0 * to_yyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxxy[k] = 3.0 * to_yyz_xxx[k] - 6.0 * to_yyz_xxxyy[k] * tke_0 - 2.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxxz[k] = -6.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxyy[k] = 6.0 * to_yyz_xxy[k] - 6.0 * to_yyz_xxyyy[k] * tke_0 - 4.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxyz[k] = 3.0 * to_yyz_xxz[k] - 6.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxzz[k] = -6.0 * to_yyz_xxyzz[k] * tke_0 + 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xyyy[k] = 9.0 * to_yyz_xyy[k] - 6.0 * to_yyz_xyyyy[k] * tke_0 - 6.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xyyz[k] = 6.0 * to_yyz_xyz[k] - 6.0 * to_yyz_xyyyz[k] * tke_0 - 4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xyzz[k] = 3.0 * to_yyz_xzz[k] - 6.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xzzz[k] = -6.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yyyy[k] = 12.0 * to_yyz_yyy[k] - 6.0 * to_yyz_yyyyy[k] * tke_0 - 8.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yyyz[k] = 9.0 * to_yyz_yyz[k] - 6.0 * to_yyz_yyyyz[k] * tke_0 - 6.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yyzz[k] = 6.0 * to_yyz_yzz[k] - 6.0 * to_yyz_yyyzz[k] * tke_0 - 4.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yzzz[k] = 3.0 * to_yyz_zzz[k] - 6.0 * to_yyz_yyzzz[k] * tke_0 - 2.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_zzzz[k] = -6.0 * to_yyz_yzzzz[k] * tke_0 + 4.0 * to_yyyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1080-1095 components of targeted buffer : GG

        auto to_y_y_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 180);

        auto to_y_y_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 181);

        auto to_y_y_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 182);

        auto to_y_y_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 183);

        auto to_y_y_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 184);

        auto to_y_y_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 185);

        auto to_y_y_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 186);

        auto to_y_y_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 187);

        auto to_y_y_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 188);

        auto to_y_y_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 189);

        auto to_y_y_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 190);

        auto to_y_y_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 191);

        auto to_y_y_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 192);

        auto to_y_y_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 193);

        auto to_y_y_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_y_y_yyzz_xxxx, to_y_y_yyzz_xxxy, to_y_y_yyzz_xxxz, to_y_y_yyzz_xxyy, to_y_y_yyzz_xxyz, to_y_y_yyzz_xxzz, to_y_y_yyzz_xyyy, to_y_y_yyzz_xyyz, to_y_y_yyzz_xyzz, to_y_y_yyzz_xzzz, to_y_y_yyzz_yyyy, to_y_y_yyzz_yyyz, to_y_y_yyzz_yyzz, to_y_y_yyzz_yzzz, to_y_y_yyzz_zzzz, to_yyyzz_xxx, to_yyyzz_xxxxy, to_yyyzz_xxxyy, to_yyyzz_xxxyz, to_yyyzz_xxy, to_yyyzz_xxyyy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xyy, to_yyyzz_xyyyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_yyy, to_yyyzz_yyyyy, to_yyyzz_yyyyz, to_yyyzz_yyyzz, to_yyyzz_yyz, to_yyyzz_yyzzz, to_yyyzz_yzz, to_yyyzz_yzzzz, to_yyyzz_zzz, to_yzz_xxx, to_yzz_xxxxy, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_yyy, to_yzz_yyyyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyzz_xxxx[k] = -4.0 * to_yzz_xxxxy[k] * tke_0 + 4.0 * to_yyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxxy[k] = 2.0 * to_yzz_xxx[k] - 4.0 * to_yzz_xxxyy[k] * tke_0 - 2.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxxz[k] = -4.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxyy[k] = 4.0 * to_yzz_xxy[k] - 4.0 * to_yzz_xxyyy[k] * tke_0 - 4.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxyz[k] = 2.0 * to_yzz_xxz[k] - 4.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxzz[k] = -4.0 * to_yzz_xxyzz[k] * tke_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xyyy[k] = 6.0 * to_yzz_xyy[k] - 4.0 * to_yzz_xyyyy[k] * tke_0 - 6.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xyyz[k] = 4.0 * to_yzz_xyz[k] - 4.0 * to_yzz_xyyyz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xyzz[k] = 2.0 * to_yzz_xzz[k] - 4.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xzzz[k] = -4.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yyyy[k] = 8.0 * to_yzz_yyy[k] - 4.0 * to_yzz_yyyyy[k] * tke_0 - 8.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yyyz[k] = 6.0 * to_yzz_yyz[k] - 4.0 * to_yzz_yyyyz[k] * tke_0 - 6.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yyzz[k] = 4.0 * to_yzz_yzz[k] - 4.0 * to_yzz_yyyzz[k] * tke_0 - 4.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yzzz[k] = 2.0 * to_yzz_zzz[k] - 4.0 * to_yzz_yyzzz[k] * tke_0 - 2.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_zzzz[k] = -4.0 * to_yzz_yzzzz[k] * tke_0 + 4.0 * to_yyyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1095-1110 components of targeted buffer : GG

        auto to_y_y_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 195);

        auto to_y_y_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 196);

        auto to_y_y_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 197);

        auto to_y_y_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 198);

        auto to_y_y_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 199);

        auto to_y_y_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 200);

        auto to_y_y_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 201);

        auto to_y_y_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 202);

        auto to_y_y_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 203);

        auto to_y_y_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 204);

        auto to_y_y_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 205);

        auto to_y_y_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 206);

        auto to_y_y_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 207);

        auto to_y_y_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 208);

        auto to_y_y_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_y_y_yzzz_xxxx, to_y_y_yzzz_xxxy, to_y_y_yzzz_xxxz, to_y_y_yzzz_xxyy, to_y_y_yzzz_xxyz, to_y_y_yzzz_xxzz, to_y_y_yzzz_xyyy, to_y_y_yzzz_xyyz, to_y_y_yzzz_xyzz, to_y_y_yzzz_xzzz, to_y_y_yzzz_yyyy, to_y_y_yzzz_yyyz, to_y_y_yzzz_yyzz, to_y_y_yzzz_yzzz, to_y_y_yzzz_zzzz, to_yyzzz_xxx, to_yyzzz_xxxxy, to_yyzzz_xxxyy, to_yyzzz_xxxyz, to_yyzzz_xxy, to_yyzzz_xxyyy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xyy, to_yyzzz_xyyyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_yyy, to_yyzzz_yyyyy, to_yyzzz_yyyyz, to_yyzzz_yyyzz, to_yyzzz_yyz, to_yyzzz_yyzzz, to_yyzzz_yzz, to_yyzzz_yzzzz, to_yyzzz_zzz, to_zzz_xxx, to_zzz_xxxxy, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_yyy, to_zzz_yyyyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzzz_xxxx[k] = -2.0 * to_zzz_xxxxy[k] * tke_0 + 4.0 * to_yyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxxy[k] = to_zzz_xxx[k] - 2.0 * to_zzz_xxxyy[k] * tke_0 - 2.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxxz[k] = -2.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxyy[k] = 2.0 * to_zzz_xxy[k] - 2.0 * to_zzz_xxyyy[k] * tke_0 - 4.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxyz[k] = to_zzz_xxz[k] - 2.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxzz[k] = -2.0 * to_zzz_xxyzz[k] * tke_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xyyy[k] = 3.0 * to_zzz_xyy[k] - 2.0 * to_zzz_xyyyy[k] * tke_0 - 6.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xyyz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xyyyz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xyzz[k] = to_zzz_xzz[k] - 2.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xzzz[k] = -2.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yyyy[k] = 4.0 * to_zzz_yyy[k] - 2.0 * to_zzz_yyyyy[k] * tke_0 - 8.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yyyz[k] = 3.0 * to_zzz_yyz[k] - 2.0 * to_zzz_yyyyz[k] * tke_0 - 6.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yyzz[k] = 2.0 * to_zzz_yzz[k] - 2.0 * to_zzz_yyyzz[k] * tke_0 - 4.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yzzz[k] = to_zzz_zzz[k] - 2.0 * to_zzz_yyzzz[k] * tke_0 - 2.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_zzzz[k] = -2.0 * to_zzz_yzzzz[k] * tke_0 + 4.0 * to_yyzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1110-1125 components of targeted buffer : GG

        auto to_y_y_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 210);

        auto to_y_y_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 211);

        auto to_y_y_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 212);

        auto to_y_y_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 213);

        auto to_y_y_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 214);

        auto to_y_y_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 215);

        auto to_y_y_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 216);

        auto to_y_y_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 217);

        auto to_y_y_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 218);

        auto to_y_y_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 219);

        auto to_y_y_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 220);

        auto to_y_y_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 221);

        auto to_y_y_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 222);

        auto to_y_y_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 223);

        auto to_y_y_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 4 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_y_y_zzzz_xxxx, to_y_y_zzzz_xxxy, to_y_y_zzzz_xxxz, to_y_y_zzzz_xxyy, to_y_y_zzzz_xxyz, to_y_y_zzzz_xxzz, to_y_y_zzzz_xyyy, to_y_y_zzzz_xyyz, to_y_y_zzzz_xyzz, to_y_y_zzzz_xzzz, to_y_y_zzzz_yyyy, to_y_y_zzzz_yyyz, to_y_y_zzzz_yyzz, to_y_y_zzzz_yzzz, to_y_y_zzzz_zzzz, to_yzzzz_xxx, to_yzzzz_xxxxy, to_yzzzz_xxxyy, to_yzzzz_xxxyz, to_yzzzz_xxy, to_yzzzz_xxyyy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xyy, to_yzzzz_xyyyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_yyy, to_yzzzz_yyyyy, to_yzzzz_yyyyz, to_yzzzz_yyyzz, to_yzzzz_yyz, to_yzzzz_yyzzz, to_yzzzz_yzz, to_yzzzz_yzzzz, to_yzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzzz_xxxx[k] = 4.0 * to_yzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxxy[k] = -2.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxxz[k] = 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxyy[k] = -4.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxyz[k] = -2.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxzz[k] = 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xyyy[k] = -6.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xyyz[k] = -4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xyzz[k] = -2.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xzzz[k] = 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yyyy[k] = -8.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yyyz[k] = -6.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yyzz[k] = -4.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yzzz[k] = -2.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_zzzz[k] = 4.0 * to_yzzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1125-1140 components of targeted buffer : GG

        auto to_y_z_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 0);

        auto to_y_z_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 1);

        auto to_y_z_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 2);

        auto to_y_z_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 3);

        auto to_y_z_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 4);

        auto to_y_z_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 5);

        auto to_y_z_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 6);

        auto to_y_z_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 7);

        auto to_y_z_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 8);

        auto to_y_z_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 9);

        auto to_y_z_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 10);

        auto to_y_z_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 11);

        auto to_y_z_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 12);

        auto to_y_z_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 13);

        auto to_y_z_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxy_xxx, to_xxxxy_xxxxz, to_xxxxy_xxxyz, to_xxxxy_xxxzz, to_xxxxy_xxy, to_xxxxy_xxyyz, to_xxxxy_xxyzz, to_xxxxy_xxz, to_xxxxy_xxzzz, to_xxxxy_xyy, to_xxxxy_xyyyz, to_xxxxy_xyyzz, to_xxxxy_xyz, to_xxxxy_xyzzz, to_xxxxy_xzz, to_xxxxy_xzzzz, to_xxxxy_yyy, to_xxxxy_yyyyz, to_xxxxy_yyyzz, to_xxxxy_yyz, to_xxxxy_yyzzz, to_xxxxy_yzz, to_xxxxy_yzzzz, to_xxxxy_zzz, to_xxxxy_zzzzz, to_y_z_xxxx_xxxx, to_y_z_xxxx_xxxy, to_y_z_xxxx_xxxz, to_y_z_xxxx_xxyy, to_y_z_xxxx_xxyz, to_y_z_xxxx_xxzz, to_y_z_xxxx_xyyy, to_y_z_xxxx_xyyz, to_y_z_xxxx_xyzz, to_y_z_xxxx_xzzz, to_y_z_xxxx_yyyy, to_y_z_xxxx_yyyz, to_y_z_xxxx_yyzz, to_y_z_xxxx_yzzz, to_y_z_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxx_xxxx[k] = 4.0 * to_xxxxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxxy[k] = 4.0 * to_xxxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxxz[k] = -2.0 * to_xxxxy_xxx[k] * tbe_0 + 4.0 * to_xxxxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxyy[k] = 4.0 * to_xxxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxyz[k] = -2.0 * to_xxxxy_xxy[k] * tbe_0 + 4.0 * to_xxxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxzz[k] = -4.0 * to_xxxxy_xxz[k] * tbe_0 + 4.0 * to_xxxxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xyyy[k] = 4.0 * to_xxxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xyyz[k] = -2.0 * to_xxxxy_xyy[k] * tbe_0 + 4.0 * to_xxxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xyzz[k] = -4.0 * to_xxxxy_xyz[k] * tbe_0 + 4.0 * to_xxxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xzzz[k] = -6.0 * to_xxxxy_xzz[k] * tbe_0 + 4.0 * to_xxxxy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yyyy[k] = 4.0 * to_xxxxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yyyz[k] = -2.0 * to_xxxxy_yyy[k] * tbe_0 + 4.0 * to_xxxxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yyzz[k] = -4.0 * to_xxxxy_yyz[k] * tbe_0 + 4.0 * to_xxxxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yzzz[k] = -6.0 * to_xxxxy_yzz[k] * tbe_0 + 4.0 * to_xxxxy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_zzzz[k] = -8.0 * to_xxxxy_zzz[k] * tbe_0 + 4.0 * to_xxxxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1140-1155 components of targeted buffer : GG

        auto to_y_z_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 15);

        auto to_y_z_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 16);

        auto to_y_z_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 17);

        auto to_y_z_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 18);

        auto to_y_z_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 19);

        auto to_y_z_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 20);

        auto to_y_z_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 21);

        auto to_y_z_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 22);

        auto to_y_z_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 23);

        auto to_y_z_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 24);

        auto to_y_z_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 25);

        auto to_y_z_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 26);

        auto to_y_z_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 27);

        auto to_y_z_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 28);

        auto to_y_z_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxz, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxx_zzzzz, to_xxxyy_xxx, to_xxxyy_xxxxz, to_xxxyy_xxxyz, to_xxxyy_xxxzz, to_xxxyy_xxy, to_xxxyy_xxyyz, to_xxxyy_xxyzz, to_xxxyy_xxz, to_xxxyy_xxzzz, to_xxxyy_xyy, to_xxxyy_xyyyz, to_xxxyy_xyyzz, to_xxxyy_xyz, to_xxxyy_xyzzz, to_xxxyy_xzz, to_xxxyy_xzzzz, to_xxxyy_yyy, to_xxxyy_yyyyz, to_xxxyy_yyyzz, to_xxxyy_yyz, to_xxxyy_yyzzz, to_xxxyy_yzz, to_xxxyy_yzzzz, to_xxxyy_zzz, to_xxxyy_zzzzz, to_y_z_xxxy_xxxx, to_y_z_xxxy_xxxy, to_y_z_xxxy_xxxz, to_y_z_xxxy_xxyy, to_y_z_xxxy_xxyz, to_y_z_xxxy_xxzz, to_y_z_xxxy_xyyy, to_y_z_xxxy_xyyz, to_y_z_xxxy_xyzz, to_y_z_xxxy_xzzz, to_y_z_xxxy_yyyy, to_y_z_xxxy_yyyz, to_y_z_xxxy_yyzz, to_y_z_xxxy_yzzz, to_y_z_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxy_xxxx[k] = -2.0 * to_xxx_xxxxz[k] * tke_0 + 4.0 * to_xxxyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxxy[k] = -2.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxxz[k] = to_xxx_xxx[k] - 2.0 * to_xxx_xxxzz[k] * tke_0 - 2.0 * to_xxxyy_xxx[k] * tbe_0 + 4.0 * to_xxxyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxyy[k] = -2.0 * to_xxx_xxyyz[k] * tke_0 + 4.0 * to_xxxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxyz[k] = to_xxx_xxy[k] - 2.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxyy_xxy[k] * tbe_0 + 4.0 * to_xxxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxzz[k] = 2.0 * to_xxx_xxz[k] - 2.0 * to_xxx_xxzzz[k] * tke_0 - 4.0 * to_xxxyy_xxz[k] * tbe_0 + 4.0 * to_xxxyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xyyy[k] = -2.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xyyz[k] = to_xxx_xyy[k] - 2.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxyy_xyy[k] * tbe_0 + 4.0 * to_xxxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xyzz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xyzzz[k] * tke_0 - 4.0 * to_xxxyy_xyz[k] * tbe_0 + 4.0 * to_xxxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xzzz[k] = 3.0 * to_xxx_xzz[k] - 2.0 * to_xxx_xzzzz[k] * tke_0 - 6.0 * to_xxxyy_xzz[k] * tbe_0 + 4.0 * to_xxxyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yyyy[k] = -2.0 * to_xxx_yyyyz[k] * tke_0 + 4.0 * to_xxxyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yyyz[k] = to_xxx_yyy[k] - 2.0 * to_xxx_yyyzz[k] * tke_0 - 2.0 * to_xxxyy_yyy[k] * tbe_0 + 4.0 * to_xxxyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yyzz[k] = 2.0 * to_xxx_yyz[k] - 2.0 * to_xxx_yyzzz[k] * tke_0 - 4.0 * to_xxxyy_yyz[k] * tbe_0 + 4.0 * to_xxxyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yzzz[k] = 3.0 * to_xxx_yzz[k] - 2.0 * to_xxx_yzzzz[k] * tke_0 - 6.0 * to_xxxyy_yzz[k] * tbe_0 + 4.0 * to_xxxyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_zzzz[k] = 4.0 * to_xxx_zzz[k] - 2.0 * to_xxx_zzzzz[k] * tke_0 - 8.0 * to_xxxyy_zzz[k] * tbe_0 + 4.0 * to_xxxyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1155-1170 components of targeted buffer : GG

        auto to_y_z_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 30);

        auto to_y_z_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 31);

        auto to_y_z_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 32);

        auto to_y_z_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 33);

        auto to_y_z_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 34);

        auto to_y_z_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 35);

        auto to_y_z_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 36);

        auto to_y_z_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 37);

        auto to_y_z_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 38);

        auto to_y_z_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 39);

        auto to_y_z_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 40);

        auto to_y_z_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 41);

        auto to_y_z_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 42);

        auto to_y_z_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 43);

        auto to_y_z_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxz, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_xxxyz_zzzzz, to_y_z_xxxz_xxxx, to_y_z_xxxz_xxxy, to_y_z_xxxz_xxxz, to_y_z_xxxz_xxyy, to_y_z_xxxz_xxyz, to_y_z_xxxz_xxzz, to_y_z_xxxz_xyyy, to_y_z_xxxz_xyyz, to_y_z_xxxz_xyzz, to_y_z_xxxz_xzzz, to_y_z_xxxz_yyyy, to_y_z_xxxz_yyyz, to_y_z_xxxz_yyzz, to_y_z_xxxz_yzzz, to_y_z_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxz_xxxx[k] = 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxxy[k] = 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxxz[k] = -2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxyy[k] = 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxyz[k] = -2.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxzz[k] = -4.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xyyy[k] = 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xyyz[k] = -2.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xyzz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xzzz[k] = -6.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yyyy[k] = 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yyyz[k] = -2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yyzz[k] = -4.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yzzz[k] = -6.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_zzzz[k] = -8.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1170-1185 components of targeted buffer : GG

        auto to_y_z_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 45);

        auto to_y_z_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 46);

        auto to_y_z_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 47);

        auto to_y_z_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 48);

        auto to_y_z_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 49);

        auto to_y_z_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 50);

        auto to_y_z_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 51);

        auto to_y_z_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 52);

        auto to_y_z_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 53);

        auto to_y_z_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 54);

        auto to_y_z_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 55);

        auto to_y_z_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 56);

        auto to_y_z_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 57);

        auto to_y_z_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 58);

        auto to_y_z_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxz, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxy_zzzzz, to_xxyyy_xxx, to_xxyyy_xxxxz, to_xxyyy_xxxyz, to_xxyyy_xxxzz, to_xxyyy_xxy, to_xxyyy_xxyyz, to_xxyyy_xxyzz, to_xxyyy_xxz, to_xxyyy_xxzzz, to_xxyyy_xyy, to_xxyyy_xyyyz, to_xxyyy_xyyzz, to_xxyyy_xyz, to_xxyyy_xyzzz, to_xxyyy_xzz, to_xxyyy_xzzzz, to_xxyyy_yyy, to_xxyyy_yyyyz, to_xxyyy_yyyzz, to_xxyyy_yyz, to_xxyyy_yyzzz, to_xxyyy_yzz, to_xxyyy_yzzzz, to_xxyyy_zzz, to_xxyyy_zzzzz, to_y_z_xxyy_xxxx, to_y_z_xxyy_xxxy, to_y_z_xxyy_xxxz, to_y_z_xxyy_xxyy, to_y_z_xxyy_xxyz, to_y_z_xxyy_xxzz, to_y_z_xxyy_xyyy, to_y_z_xxyy_xyyz, to_y_z_xxyy_xyzz, to_y_z_xxyy_xzzz, to_y_z_xxyy_yyyy, to_y_z_xxyy_yyyz, to_y_z_xxyy_yyzz, to_y_z_xxyy_yzzz, to_y_z_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyy_xxxx[k] = -4.0 * to_xxy_xxxxz[k] * tke_0 + 4.0 * to_xxyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxxy[k] = -4.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxxz[k] = 2.0 * to_xxy_xxx[k] - 4.0 * to_xxy_xxxzz[k] * tke_0 - 2.0 * to_xxyyy_xxx[k] * tbe_0 + 4.0 * to_xxyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxyy[k] = -4.0 * to_xxy_xxyyz[k] * tke_0 + 4.0 * to_xxyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxyz[k] = 2.0 * to_xxy_xxy[k] - 4.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxyyy_xxy[k] * tbe_0 + 4.0 * to_xxyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxzz[k] = 4.0 * to_xxy_xxz[k] - 4.0 * to_xxy_xxzzz[k] * tke_0 - 4.0 * to_xxyyy_xxz[k] * tbe_0 + 4.0 * to_xxyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xyyy[k] = -4.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xyyz[k] = 2.0 * to_xxy_xyy[k] - 4.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxyyy_xyy[k] * tbe_0 + 4.0 * to_xxyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xyzz[k] = 4.0 * to_xxy_xyz[k] - 4.0 * to_xxy_xyzzz[k] * tke_0 - 4.0 * to_xxyyy_xyz[k] * tbe_0 + 4.0 * to_xxyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xzzz[k] = 6.0 * to_xxy_xzz[k] - 4.0 * to_xxy_xzzzz[k] * tke_0 - 6.0 * to_xxyyy_xzz[k] * tbe_0 + 4.0 * to_xxyyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yyyy[k] = -4.0 * to_xxy_yyyyz[k] * tke_0 + 4.0 * to_xxyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yyyz[k] = 2.0 * to_xxy_yyy[k] - 4.0 * to_xxy_yyyzz[k] * tke_0 - 2.0 * to_xxyyy_yyy[k] * tbe_0 + 4.0 * to_xxyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yyzz[k] = 4.0 * to_xxy_yyz[k] - 4.0 * to_xxy_yyzzz[k] * tke_0 - 4.0 * to_xxyyy_yyz[k] * tbe_0 + 4.0 * to_xxyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yzzz[k] = 6.0 * to_xxy_yzz[k] - 4.0 * to_xxy_yzzzz[k] * tke_0 - 6.0 * to_xxyyy_yzz[k] * tbe_0 + 4.0 * to_xxyyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_zzzz[k] = 8.0 * to_xxy_zzz[k] - 4.0 * to_xxy_zzzzz[k] * tke_0 - 8.0 * to_xxyyy_zzz[k] * tbe_0 + 4.0 * to_xxyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1185-1200 components of targeted buffer : GG

        auto to_y_z_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 60);

        auto to_y_z_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 61);

        auto to_y_z_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 62);

        auto to_y_z_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 63);

        auto to_y_z_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 64);

        auto to_y_z_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 65);

        auto to_y_z_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 66);

        auto to_y_z_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 67);

        auto to_y_z_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 68);

        auto to_y_z_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 69);

        auto to_y_z_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 70);

        auto to_y_z_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 71);

        auto to_y_z_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 72);

        auto to_y_z_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 73);

        auto to_y_z_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxz, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_xxyyz_zzzzz, to_xxz_xxx, to_xxz_xxxxz, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_xxz_zzzzz, to_y_z_xxyz_xxxx, to_y_z_xxyz_xxxy, to_y_z_xxyz_xxxz, to_y_z_xxyz_xxyy, to_y_z_xxyz_xxyz, to_y_z_xxyz_xxzz, to_y_z_xxyz_xyyy, to_y_z_xxyz_xyyz, to_y_z_xxyz_xyzz, to_y_z_xxyz_xzzz, to_y_z_xxyz_yyyy, to_y_z_xxyz_yyyz, to_y_z_xxyz_yyzz, to_y_z_xxyz_yzzz, to_y_z_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyz_xxxx[k] = -2.0 * to_xxz_xxxxz[k] * tke_0 + 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxxy[k] = -2.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxxz[k] = to_xxz_xxx[k] - 2.0 * to_xxz_xxxzz[k] * tke_0 - 2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxyy[k] = -2.0 * to_xxz_xxyyz[k] * tke_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxyz[k] = to_xxz_xxy[k] - 2.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxzz[k] = 2.0 * to_xxz_xxz[k] - 2.0 * to_xxz_xxzzz[k] * tke_0 - 4.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xyyy[k] = -2.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xyyz[k] = to_xxz_xyy[k] - 2.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xyzz[k] = 2.0 * to_xxz_xyz[k] - 2.0 * to_xxz_xyzzz[k] * tke_0 - 4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xzzz[k] = 3.0 * to_xxz_xzz[k] - 2.0 * to_xxz_xzzzz[k] * tke_0 - 6.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yyyy[k] = -2.0 * to_xxz_yyyyz[k] * tke_0 + 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yyyz[k] = to_xxz_yyy[k] - 2.0 * to_xxz_yyyzz[k] * tke_0 - 2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yyzz[k] = 2.0 * to_xxz_yyz[k] - 2.0 * to_xxz_yyzzz[k] * tke_0 - 4.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yzzz[k] = 3.0 * to_xxz_yzz[k] - 2.0 * to_xxz_yzzzz[k] * tke_0 - 6.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_zzzz[k] = 4.0 * to_xxz_zzz[k] - 2.0 * to_xxz_zzzzz[k] * tke_0 - 8.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1200-1215 components of targeted buffer : GG

        auto to_y_z_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 75);

        auto to_y_z_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 76);

        auto to_y_z_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 77);

        auto to_y_z_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 78);

        auto to_y_z_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 79);

        auto to_y_z_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 80);

        auto to_y_z_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 81);

        auto to_y_z_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 82);

        auto to_y_z_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 83);

        auto to_y_z_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 84);

        auto to_y_z_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 85);

        auto to_y_z_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 86);

        auto to_y_z_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 87);

        auto to_y_z_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 88);

        auto to_y_z_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxyzz_xxx, to_xxyzz_xxxxz, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_xxyzz_zzzzz, to_y_z_xxzz_xxxx, to_y_z_xxzz_xxxy, to_y_z_xxzz_xxxz, to_y_z_xxzz_xxyy, to_y_z_xxzz_xxyz, to_y_z_xxzz_xxzz, to_y_z_xxzz_xyyy, to_y_z_xxzz_xyyz, to_y_z_xxzz_xyzz, to_y_z_xxzz_xzzz, to_y_z_xxzz_yyyy, to_y_z_xxzz_yyyz, to_y_z_xxzz_yyzz, to_y_z_xxzz_yzzz, to_y_z_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxzz_xxxx[k] = 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxxy[k] = 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxxz[k] = -2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxyy[k] = 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxyz[k] = -2.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxzz[k] = -4.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xyyy[k] = 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xyyz[k] = -2.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xyzz[k] = -4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xzzz[k] = -6.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yyyy[k] = 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yyyz[k] = -2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yyzz[k] = -4.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yzzz[k] = -6.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_zzzz[k] = -8.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1215-1230 components of targeted buffer : GG

        auto to_y_z_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 90);

        auto to_y_z_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 91);

        auto to_y_z_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 92);

        auto to_y_z_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 93);

        auto to_y_z_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 94);

        auto to_y_z_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 95);

        auto to_y_z_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 96);

        auto to_y_z_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 97);

        auto to_y_z_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 98);

        auto to_y_z_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 99);

        auto to_y_z_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 100);

        auto to_y_z_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 101);

        auto to_y_z_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 102);

        auto to_y_z_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 103);

        auto to_y_z_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxz, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyy_zzzzz, to_xyyyy_xxx, to_xyyyy_xxxxz, to_xyyyy_xxxyz, to_xyyyy_xxxzz, to_xyyyy_xxy, to_xyyyy_xxyyz, to_xyyyy_xxyzz, to_xyyyy_xxz, to_xyyyy_xxzzz, to_xyyyy_xyy, to_xyyyy_xyyyz, to_xyyyy_xyyzz, to_xyyyy_xyz, to_xyyyy_xyzzz, to_xyyyy_xzz, to_xyyyy_xzzzz, to_xyyyy_yyy, to_xyyyy_yyyyz, to_xyyyy_yyyzz, to_xyyyy_yyz, to_xyyyy_yyzzz, to_xyyyy_yzz, to_xyyyy_yzzzz, to_xyyyy_zzz, to_xyyyy_zzzzz, to_y_z_xyyy_xxxx, to_y_z_xyyy_xxxy, to_y_z_xyyy_xxxz, to_y_z_xyyy_xxyy, to_y_z_xyyy_xxyz, to_y_z_xyyy_xxzz, to_y_z_xyyy_xyyy, to_y_z_xyyy_xyyz, to_y_z_xyyy_xyzz, to_y_z_xyyy_xzzz, to_y_z_xyyy_yyyy, to_y_z_xyyy_yyyz, to_y_z_xyyy_yyzz, to_y_z_xyyy_yzzz, to_y_z_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyy_xxxx[k] = -6.0 * to_xyy_xxxxz[k] * tke_0 + 4.0 * to_xyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxxy[k] = -6.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxxz[k] = 3.0 * to_xyy_xxx[k] - 6.0 * to_xyy_xxxzz[k] * tke_0 - 2.0 * to_xyyyy_xxx[k] * tbe_0 + 4.0 * to_xyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxyy[k] = -6.0 * to_xyy_xxyyz[k] * tke_0 + 4.0 * to_xyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxyz[k] = 3.0 * to_xyy_xxy[k] - 6.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xyyyy_xxy[k] * tbe_0 + 4.0 * to_xyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxzz[k] = 6.0 * to_xyy_xxz[k] - 6.0 * to_xyy_xxzzz[k] * tke_0 - 4.0 * to_xyyyy_xxz[k] * tbe_0 + 4.0 * to_xyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xyyy[k] = -6.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xyyz[k] = 3.0 * to_xyy_xyy[k] - 6.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xyyyy_xyy[k] * tbe_0 + 4.0 * to_xyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xyzz[k] = 6.0 * to_xyy_xyz[k] - 6.0 * to_xyy_xyzzz[k] * tke_0 - 4.0 * to_xyyyy_xyz[k] * tbe_0 + 4.0 * to_xyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xzzz[k] = 9.0 * to_xyy_xzz[k] - 6.0 * to_xyy_xzzzz[k] * tke_0 - 6.0 * to_xyyyy_xzz[k] * tbe_0 + 4.0 * to_xyyyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yyyy[k] = -6.0 * to_xyy_yyyyz[k] * tke_0 + 4.0 * to_xyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yyyz[k] = 3.0 * to_xyy_yyy[k] - 6.0 * to_xyy_yyyzz[k] * tke_0 - 2.0 * to_xyyyy_yyy[k] * tbe_0 + 4.0 * to_xyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yyzz[k] = 6.0 * to_xyy_yyz[k] - 6.0 * to_xyy_yyzzz[k] * tke_0 - 4.0 * to_xyyyy_yyz[k] * tbe_0 + 4.0 * to_xyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yzzz[k] = 9.0 * to_xyy_yzz[k] - 6.0 * to_xyy_yzzzz[k] * tke_0 - 6.0 * to_xyyyy_yzz[k] * tbe_0 + 4.0 * to_xyyyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_zzzz[k] = 12.0 * to_xyy_zzz[k] - 6.0 * to_xyy_zzzzz[k] * tke_0 - 8.0 * to_xyyyy_zzz[k] * tbe_0 + 4.0 * to_xyyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1230-1245 components of targeted buffer : GG

        auto to_y_z_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 105);

        auto to_y_z_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 106);

        auto to_y_z_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 107);

        auto to_y_z_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 108);

        auto to_y_z_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 109);

        auto to_y_z_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 110);

        auto to_y_z_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 111);

        auto to_y_z_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 112);

        auto to_y_z_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 113);

        auto to_y_z_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 114);

        auto to_y_z_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 115);

        auto to_y_z_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 116);

        auto to_y_z_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 117);

        auto to_y_z_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 118);

        auto to_y_z_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxz, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, to_xyyyz_zzzzz, to_xyz_xxx, to_xyz_xxxxz, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_xyz_zzzzz, to_y_z_xyyz_xxxx, to_y_z_xyyz_xxxy, to_y_z_xyyz_xxxz, to_y_z_xyyz_xxyy, to_y_z_xyyz_xxyz, to_y_z_xyyz_xxzz, to_y_z_xyyz_xyyy, to_y_z_xyyz_xyyz, to_y_z_xyyz_xyzz, to_y_z_xyyz_xzzz, to_y_z_xyyz_yyyy, to_y_z_xyyz_yyyz, to_y_z_xyyz_yyzz, to_y_z_xyyz_yzzz, to_y_z_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyz_xxxx[k] = -4.0 * to_xyz_xxxxz[k] * tke_0 + 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxxy[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxxz[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxyy[k] = -4.0 * to_xyz_xxyyz[k] * tke_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxyz[k] = 2.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxzz[k] = 4.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 4.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xyyy[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xyyz[k] = 2.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xyzz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyzzz[k] * tke_0 - 4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xzzz[k] = 6.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xzzzz[k] * tke_0 - 6.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yyyy[k] = -4.0 * to_xyz_yyyyz[k] * tke_0 + 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yyyz[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yyzz[k] = 4.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 4.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yzzz[k] = 6.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yzzzz[k] * tke_0 - 6.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_zzzz[k] = 8.0 * to_xyz_zzz[k] - 4.0 * to_xyz_zzzzz[k] * tke_0 - 8.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1245-1260 components of targeted buffer : GG

        auto to_y_z_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 120);

        auto to_y_z_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 121);

        auto to_y_z_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 122);

        auto to_y_z_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 123);

        auto to_y_z_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 124);

        auto to_y_z_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 125);

        auto to_y_z_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 126);

        auto to_y_z_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 127);

        auto to_y_z_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 128);

        auto to_y_z_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 129);

        auto to_y_z_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 130);

        auto to_y_z_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 131);

        auto to_y_z_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 132);

        auto to_y_z_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 133);

        auto to_y_z_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyyzz_xxx, to_xyyzz_xxxxz, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, to_xyyzz_zzzzz, to_xzz_xxx, to_xzz_xxxxz, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_xzz_zzzzz, to_y_z_xyzz_xxxx, to_y_z_xyzz_xxxy, to_y_z_xyzz_xxxz, to_y_z_xyzz_xxyy, to_y_z_xyzz_xxyz, to_y_z_xyzz_xxzz, to_y_z_xyzz_xyyy, to_y_z_xyzz_xyyz, to_y_z_xyzz_xyzz, to_y_z_xyzz_xzzz, to_y_z_xyzz_yyyy, to_y_z_xyzz_yyyz, to_y_z_xyzz_yyzz, to_y_z_xyzz_yzzz, to_y_z_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyzz_xxxx[k] = -2.0 * to_xzz_xxxxz[k] * tke_0 + 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxxy[k] = -2.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxxz[k] = to_xzz_xxx[k] - 2.0 * to_xzz_xxxzz[k] * tke_0 - 2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxyy[k] = -2.0 * to_xzz_xxyyz[k] * tke_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxyz[k] = to_xzz_xxy[k] - 2.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxzz[k] = 2.0 * to_xzz_xxz[k] - 2.0 * to_xzz_xxzzz[k] * tke_0 - 4.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xyyy[k] = -2.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xyyz[k] = to_xzz_xyy[k] - 2.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xyzz[k] = 2.0 * to_xzz_xyz[k] - 2.0 * to_xzz_xyzzz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xzzz[k] = 3.0 * to_xzz_xzz[k] - 2.0 * to_xzz_xzzzz[k] * tke_0 - 6.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yyyy[k] = -2.0 * to_xzz_yyyyz[k] * tke_0 + 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yyyz[k] = to_xzz_yyy[k] - 2.0 * to_xzz_yyyzz[k] * tke_0 - 2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yyzz[k] = 2.0 * to_xzz_yyz[k] - 2.0 * to_xzz_yyzzz[k] * tke_0 - 4.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yzzz[k] = 3.0 * to_xzz_yzz[k] - 2.0 * to_xzz_yzzzz[k] * tke_0 - 6.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_zzzz[k] = 4.0 * to_xzz_zzz[k] - 2.0 * to_xzz_zzzzz[k] * tke_0 - 8.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1260-1275 components of targeted buffer : GG

        auto to_y_z_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 135);

        auto to_y_z_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 136);

        auto to_y_z_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 137);

        auto to_y_z_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 138);

        auto to_y_z_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 139);

        auto to_y_z_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 140);

        auto to_y_z_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 141);

        auto to_y_z_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 142);

        auto to_y_z_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 143);

        auto to_y_z_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 144);

        auto to_y_z_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 145);

        auto to_y_z_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 146);

        auto to_y_z_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 147);

        auto to_y_z_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 148);

        auto to_y_z_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xyzzz_xxx, to_xyzzz_xxxxz, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, to_xyzzz_zzzzz, to_y_z_xzzz_xxxx, to_y_z_xzzz_xxxy, to_y_z_xzzz_xxxz, to_y_z_xzzz_xxyy, to_y_z_xzzz_xxyz, to_y_z_xzzz_xxzz, to_y_z_xzzz_xyyy, to_y_z_xzzz_xyyz, to_y_z_xzzz_xyzz, to_y_z_xzzz_xzzz, to_y_z_xzzz_yyyy, to_y_z_xzzz_yyyz, to_y_z_xzzz_yyzz, to_y_z_xzzz_yzzz, to_y_z_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzzz_xxxx[k] = 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxxy[k] = 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxxz[k] = -2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxyy[k] = 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxyz[k] = -2.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxzz[k] = -4.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xyyy[k] = 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xyyz[k] = -2.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xyzz[k] = -4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xzzz[k] = -6.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yyyy[k] = 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yyyz[k] = -2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yyzz[k] = -4.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yzzz[k] = -6.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_zzzz[k] = -8.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1275-1290 components of targeted buffer : GG

        auto to_y_z_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 150);

        auto to_y_z_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 151);

        auto to_y_z_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 152);

        auto to_y_z_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 153);

        auto to_y_z_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 154);

        auto to_y_z_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 155);

        auto to_y_z_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 156);

        auto to_y_z_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 157);

        auto to_y_z_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 158);

        auto to_y_z_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 159);

        auto to_y_z_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 160);

        auto to_y_z_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 161);

        auto to_y_z_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 162);

        auto to_y_z_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 163);

        auto to_y_z_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_y_z_yyyy_xxxx, to_y_z_yyyy_xxxy, to_y_z_yyyy_xxxz, to_y_z_yyyy_xxyy, to_y_z_yyyy_xxyz, to_y_z_yyyy_xxzz, to_y_z_yyyy_xyyy, to_y_z_yyyy_xyyz, to_y_z_yyyy_xyzz, to_y_z_yyyy_xzzz, to_y_z_yyyy_yyyy, to_y_z_yyyy_yyyz, to_y_z_yyyy_yyzz, to_y_z_yyyy_yzzz, to_y_z_yyyy_zzzz, to_yyy_xxx, to_yyy_xxxxz, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyy_zzzzz, to_yyyyy_xxx, to_yyyyy_xxxxz, to_yyyyy_xxxyz, to_yyyyy_xxxzz, to_yyyyy_xxy, to_yyyyy_xxyyz, to_yyyyy_xxyzz, to_yyyyy_xxz, to_yyyyy_xxzzz, to_yyyyy_xyy, to_yyyyy_xyyyz, to_yyyyy_xyyzz, to_yyyyy_xyz, to_yyyyy_xyzzz, to_yyyyy_xzz, to_yyyyy_xzzzz, to_yyyyy_yyy, to_yyyyy_yyyyz, to_yyyyy_yyyzz, to_yyyyy_yyz, to_yyyyy_yyzzz, to_yyyyy_yzz, to_yyyyy_yzzzz, to_yyyyy_zzz, to_yyyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyy_xxxx[k] = -8.0 * to_yyy_xxxxz[k] * tke_0 + 4.0 * to_yyyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxxy[k] = -8.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_yyyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxxz[k] = 4.0 * to_yyy_xxx[k] - 8.0 * to_yyy_xxxzz[k] * tke_0 - 2.0 * to_yyyyy_xxx[k] * tbe_0 + 4.0 * to_yyyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxyy[k] = -8.0 * to_yyy_xxyyz[k] * tke_0 + 4.0 * to_yyyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxyz[k] = 4.0 * to_yyy_xxy[k] - 8.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_yyyyy_xxy[k] * tbe_0 + 4.0 * to_yyyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxzz[k] = 8.0 * to_yyy_xxz[k] - 8.0 * to_yyy_xxzzz[k] * tke_0 - 4.0 * to_yyyyy_xxz[k] * tbe_0 + 4.0 * to_yyyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xyyy[k] = -8.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_yyyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xyyz[k] = 4.0 * to_yyy_xyy[k] - 8.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_yyyyy_xyy[k] * tbe_0 + 4.0 * to_yyyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xyzz[k] = 8.0 * to_yyy_xyz[k] - 8.0 * to_yyy_xyzzz[k] * tke_0 - 4.0 * to_yyyyy_xyz[k] * tbe_0 + 4.0 * to_yyyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xzzz[k] = 12.0 * to_yyy_xzz[k] - 8.0 * to_yyy_xzzzz[k] * tke_0 - 6.0 * to_yyyyy_xzz[k] * tbe_0 + 4.0 * to_yyyyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yyyy[k] = -8.0 * to_yyy_yyyyz[k] * tke_0 + 4.0 * to_yyyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yyyz[k] = 4.0 * to_yyy_yyy[k] - 8.0 * to_yyy_yyyzz[k] * tke_0 - 2.0 * to_yyyyy_yyy[k] * tbe_0 + 4.0 * to_yyyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yyzz[k] = 8.0 * to_yyy_yyz[k] - 8.0 * to_yyy_yyzzz[k] * tke_0 - 4.0 * to_yyyyy_yyz[k] * tbe_0 + 4.0 * to_yyyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yzzz[k] = 12.0 * to_yyy_yzz[k] - 8.0 * to_yyy_yzzzz[k] * tke_0 - 6.0 * to_yyyyy_yzz[k] * tbe_0 + 4.0 * to_yyyyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_zzzz[k] = 16.0 * to_yyy_zzz[k] - 8.0 * to_yyy_zzzzz[k] * tke_0 - 8.0 * to_yyyyy_zzz[k] * tbe_0 + 4.0 * to_yyyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1290-1305 components of targeted buffer : GG

        auto to_y_z_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 165);

        auto to_y_z_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 166);

        auto to_y_z_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 167);

        auto to_y_z_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 168);

        auto to_y_z_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 169);

        auto to_y_z_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 170);

        auto to_y_z_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 171);

        auto to_y_z_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 172);

        auto to_y_z_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 173);

        auto to_y_z_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 174);

        auto to_y_z_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 175);

        auto to_y_z_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 176);

        auto to_y_z_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 177);

        auto to_y_z_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 178);

        auto to_y_z_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_y_z_yyyz_xxxx, to_y_z_yyyz_xxxy, to_y_z_yyyz_xxxz, to_y_z_yyyz_xxyy, to_y_z_yyyz_xxyz, to_y_z_yyyz_xxzz, to_y_z_yyyz_xyyy, to_y_z_yyyz_xyyz, to_y_z_yyyz_xyzz, to_y_z_yyyz_xzzz, to_y_z_yyyz_yyyy, to_y_z_yyyz_yyyz, to_y_z_yyyz_yyzz, to_y_z_yyyz_yzzz, to_y_z_yyyz_zzzz, to_yyyyz_xxx, to_yyyyz_xxxxz, to_yyyyz_xxxyz, to_yyyyz_xxxzz, to_yyyyz_xxy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xxzzz, to_yyyyz_xyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_xzzzz, to_yyyyz_yyy, to_yyyyz_yyyyz, to_yyyyz_yyyzz, to_yyyyz_yyz, to_yyyyz_yyzzz, to_yyyyz_yzz, to_yyyyz_yzzzz, to_yyyyz_zzz, to_yyyyz_zzzzz, to_yyz_xxx, to_yyz_xxxxz, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, to_yyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyz_xxxx[k] = -6.0 * to_yyz_xxxxz[k] * tke_0 + 4.0 * to_yyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxxy[k] = -6.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxxz[k] = 3.0 * to_yyz_xxx[k] - 6.0 * to_yyz_xxxzz[k] * tke_0 - 2.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxyy[k] = -6.0 * to_yyz_xxyyz[k] * tke_0 + 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxyz[k] = 3.0 * to_yyz_xxy[k] - 6.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxzz[k] = 6.0 * to_yyz_xxz[k] - 6.0 * to_yyz_xxzzz[k] * tke_0 - 4.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xyyy[k] = -6.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xyyz[k] = 3.0 * to_yyz_xyy[k] - 6.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xyzz[k] = 6.0 * to_yyz_xyz[k] - 6.0 * to_yyz_xyzzz[k] * tke_0 - 4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xzzz[k] = 9.0 * to_yyz_xzz[k] - 6.0 * to_yyz_xzzzz[k] * tke_0 - 6.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yyyy[k] = -6.0 * to_yyz_yyyyz[k] * tke_0 + 4.0 * to_yyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yyyz[k] = 3.0 * to_yyz_yyy[k] - 6.0 * to_yyz_yyyzz[k] * tke_0 - 2.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yyzz[k] = 6.0 * to_yyz_yyz[k] - 6.0 * to_yyz_yyzzz[k] * tke_0 - 4.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yzzz[k] = 9.0 * to_yyz_yzz[k] - 6.0 * to_yyz_yzzzz[k] * tke_0 - 6.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_zzzz[k] = 12.0 * to_yyz_zzz[k] - 6.0 * to_yyz_zzzzz[k] * tke_0 - 8.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1305-1320 components of targeted buffer : GG

        auto to_y_z_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 180);

        auto to_y_z_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 181);

        auto to_y_z_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 182);

        auto to_y_z_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 183);

        auto to_y_z_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 184);

        auto to_y_z_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 185);

        auto to_y_z_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 186);

        auto to_y_z_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 187);

        auto to_y_z_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 188);

        auto to_y_z_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 189);

        auto to_y_z_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 190);

        auto to_y_z_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 191);

        auto to_y_z_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 192);

        auto to_y_z_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 193);

        auto to_y_z_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_y_z_yyzz_xxxx, to_y_z_yyzz_xxxy, to_y_z_yyzz_xxxz, to_y_z_yyzz_xxyy, to_y_z_yyzz_xxyz, to_y_z_yyzz_xxzz, to_y_z_yyzz_xyyy, to_y_z_yyzz_xyyz, to_y_z_yyzz_xyzz, to_y_z_yyzz_xzzz, to_y_z_yyzz_yyyy, to_y_z_yyzz_yyyz, to_y_z_yyzz_yyzz, to_y_z_yyzz_yzzz, to_y_z_yyzz_zzzz, to_yyyzz_xxx, to_yyyzz_xxxxz, to_yyyzz_xxxyz, to_yyyzz_xxxzz, to_yyyzz_xxy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xxzzz, to_yyyzz_xyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_xzzzz, to_yyyzz_yyy, to_yyyzz_yyyyz, to_yyyzz_yyyzz, to_yyyzz_yyz, to_yyyzz_yyzzz, to_yyyzz_yzz, to_yyyzz_yzzzz, to_yyyzz_zzz, to_yyyzz_zzzzz, to_yzz_xxx, to_yzz_xxxxz, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, to_yzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyzz_xxxx[k] = -4.0 * to_yzz_xxxxz[k] * tke_0 + 4.0 * to_yyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxxy[k] = -4.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxxz[k] = 2.0 * to_yzz_xxx[k] - 4.0 * to_yzz_xxxzz[k] * tke_0 - 2.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxyy[k] = -4.0 * to_yzz_xxyyz[k] * tke_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxyz[k] = 2.0 * to_yzz_xxy[k] - 4.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxzz[k] = 4.0 * to_yzz_xxz[k] - 4.0 * to_yzz_xxzzz[k] * tke_0 - 4.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xyyy[k] = -4.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xyyz[k] = 2.0 * to_yzz_xyy[k] - 4.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xyzz[k] = 4.0 * to_yzz_xyz[k] - 4.0 * to_yzz_xyzzz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xzzz[k] = 6.0 * to_yzz_xzz[k] - 4.0 * to_yzz_xzzzz[k] * tke_0 - 6.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yyyy[k] = -4.0 * to_yzz_yyyyz[k] * tke_0 + 4.0 * to_yyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yyyz[k] = 2.0 * to_yzz_yyy[k] - 4.0 * to_yzz_yyyzz[k] * tke_0 - 2.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yyzz[k] = 4.0 * to_yzz_yyz[k] - 4.0 * to_yzz_yyzzz[k] * tke_0 - 4.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yzzz[k] = 6.0 * to_yzz_yzz[k] - 4.0 * to_yzz_yzzzz[k] * tke_0 - 6.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_zzzz[k] = 8.0 * to_yzz_zzz[k] - 4.0 * to_yzz_zzzzz[k] * tke_0 - 8.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1320-1335 components of targeted buffer : GG

        auto to_y_z_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 195);

        auto to_y_z_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 196);

        auto to_y_z_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 197);

        auto to_y_z_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 198);

        auto to_y_z_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 199);

        auto to_y_z_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 200);

        auto to_y_z_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 201);

        auto to_y_z_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 202);

        auto to_y_z_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 203);

        auto to_y_z_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 204);

        auto to_y_z_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 205);

        auto to_y_z_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 206);

        auto to_y_z_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 207);

        auto to_y_z_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 208);

        auto to_y_z_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_y_z_yzzz_xxxx, to_y_z_yzzz_xxxy, to_y_z_yzzz_xxxz, to_y_z_yzzz_xxyy, to_y_z_yzzz_xxyz, to_y_z_yzzz_xxzz, to_y_z_yzzz_xyyy, to_y_z_yzzz_xyyz, to_y_z_yzzz_xyzz, to_y_z_yzzz_xzzz, to_y_z_yzzz_yyyy, to_y_z_yzzz_yyyz, to_y_z_yzzz_yyzz, to_y_z_yzzz_yzzz, to_y_z_yzzz_zzzz, to_yyzzz_xxx, to_yyzzz_xxxxz, to_yyzzz_xxxyz, to_yyzzz_xxxzz, to_yyzzz_xxy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xxzzz, to_yyzzz_xyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_xzzzz, to_yyzzz_yyy, to_yyzzz_yyyyz, to_yyzzz_yyyzz, to_yyzzz_yyz, to_yyzzz_yyzzz, to_yyzzz_yzz, to_yyzzz_yzzzz, to_yyzzz_zzz, to_yyzzz_zzzzz, to_zzz_xxx, to_zzz_xxxxz, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, to_zzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzzz_xxxx[k] = -2.0 * to_zzz_xxxxz[k] * tke_0 + 4.0 * to_yyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxxy[k] = -2.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxxz[k] = to_zzz_xxx[k] - 2.0 * to_zzz_xxxzz[k] * tke_0 - 2.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxyy[k] = -2.0 * to_zzz_xxyyz[k] * tke_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxyz[k] = to_zzz_xxy[k] - 2.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxzz[k] = 2.0 * to_zzz_xxz[k] - 2.0 * to_zzz_xxzzz[k] * tke_0 - 4.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xyyy[k] = -2.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xyyz[k] = to_zzz_xyy[k] - 2.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xyzz[k] = 2.0 * to_zzz_xyz[k] - 2.0 * to_zzz_xyzzz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xzzz[k] = 3.0 * to_zzz_xzz[k] - 2.0 * to_zzz_xzzzz[k] * tke_0 - 6.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yyyy[k] = -2.0 * to_zzz_yyyyz[k] * tke_0 + 4.0 * to_yyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yyyz[k] = to_zzz_yyy[k] - 2.0 * to_zzz_yyyzz[k] * tke_0 - 2.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yyzz[k] = 2.0 * to_zzz_yyz[k] - 2.0 * to_zzz_yyzzz[k] * tke_0 - 4.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yzzz[k] = 3.0 * to_zzz_yzz[k] - 2.0 * to_zzz_yzzzz[k] * tke_0 - 6.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_zzzz[k] = 4.0 * to_zzz_zzz[k] - 2.0 * to_zzz_zzzzz[k] * tke_0 - 8.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1335-1350 components of targeted buffer : GG

        auto to_y_z_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 210);

        auto to_y_z_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 211);

        auto to_y_z_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 212);

        auto to_y_z_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 213);

        auto to_y_z_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 214);

        auto to_y_z_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 215);

        auto to_y_z_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 216);

        auto to_y_z_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 217);

        auto to_y_z_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 218);

        auto to_y_z_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 219);

        auto to_y_z_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 220);

        auto to_y_z_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 221);

        auto to_y_z_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 222);

        auto to_y_z_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 223);

        auto to_y_z_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 5 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_y_z_zzzz_xxxx, to_y_z_zzzz_xxxy, to_y_z_zzzz_xxxz, to_y_z_zzzz_xxyy, to_y_z_zzzz_xxyz, to_y_z_zzzz_xxzz, to_y_z_zzzz_xyyy, to_y_z_zzzz_xyyz, to_y_z_zzzz_xyzz, to_y_z_zzzz_xzzz, to_y_z_zzzz_yyyy, to_y_z_zzzz_yyyz, to_y_z_zzzz_yyzz, to_y_z_zzzz_yzzz, to_y_z_zzzz_zzzz, to_yzzzz_xxx, to_yzzzz_xxxxz, to_yzzzz_xxxyz, to_yzzzz_xxxzz, to_yzzzz_xxy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xxzzz, to_yzzzz_xyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_xzzzz, to_yzzzz_yyy, to_yzzzz_yyyyz, to_yzzzz_yyyzz, to_yzzzz_yyz, to_yzzzz_yyzzz, to_yzzzz_yzz, to_yzzzz_yzzzz, to_yzzzz_zzz, to_yzzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzzz_xxxx[k] = 4.0 * to_yzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxxy[k] = 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxxz[k] = -2.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxyy[k] = 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxyz[k] = -2.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxzz[k] = -4.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xyyy[k] = 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xyyz[k] = -2.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xyzz[k] = -4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xzzz[k] = -6.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yyyy[k] = 4.0 * to_yzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yyyz[k] = -2.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yyzz[k] = -4.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yzzz[k] = -6.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_zzzz[k] = -8.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1350-1365 components of targeted buffer : GG

        auto to_z_x_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 0);

        auto to_z_x_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 1);

        auto to_z_x_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 2);

        auto to_z_x_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 3);

        auto to_z_x_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 4);

        auto to_z_x_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 5);

        auto to_z_x_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 6);

        auto to_z_x_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 7);

        auto to_z_x_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 8);

        auto to_z_x_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 9);

        auto to_z_x_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 10);

        auto to_z_x_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 11);

        auto to_z_x_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 12);

        auto to_z_x_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 13);

        auto to_z_x_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxz_xxx, to_xxxxz_xxxxx, to_xxxxz_xxxxy, to_xxxxz_xxxxz, to_xxxxz_xxxyy, to_xxxxz_xxxyz, to_xxxxz_xxxzz, to_xxxxz_xxy, to_xxxxz_xxyyy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xxzzz, to_xxxxz_xyy, to_xxxxz_xyyyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_xzzzz, to_xxxxz_yyy, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_zzz, to_z_x_xxxx_xxxx, to_z_x_xxxx_xxxy, to_z_x_xxxx_xxxz, to_z_x_xxxx_xxyy, to_z_x_xxxx_xxyz, to_z_x_xxxx_xxzz, to_z_x_xxxx_xyyy, to_z_x_xxxx_xyyz, to_z_x_xxxx_xyzz, to_z_x_xxxx_xzzz, to_z_x_xxxx_yyyy, to_z_x_xxxx_yyyz, to_z_x_xxxx_yyzz, to_z_x_xxxx_yzzz, to_z_x_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxx_xxxx[k] = -8.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxxy[k] = -6.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxxz[k] = -6.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxyy[k] = -4.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxyz[k] = -4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxzz[k] = -4.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xyyy[k] = -2.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xyyz[k] = -2.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xyzz[k] = -2.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xzzz[k] = -2.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yyyy[k] = 4.0 * to_xxxxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yyyz[k] = 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yyzz[k] = 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yzzz[k] = 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_zzzz[k] = 4.0 * to_xxxxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1365-1380 components of targeted buffer : GG

        auto to_z_x_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 15);

        auto to_z_x_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 16);

        auto to_z_x_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 17);

        auto to_z_x_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 18);

        auto to_z_x_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 19);

        auto to_z_x_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 20);

        auto to_z_x_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 21);

        auto to_z_x_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 22);

        auto to_z_x_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 23);

        auto to_z_x_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 24);

        auto to_z_x_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 25);

        auto to_z_x_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 26);

        auto to_z_x_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 27);

        auto to_z_x_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 28);

        auto to_z_x_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxx, to_xxxyz_xxxxy, to_xxxyz_xxxxz, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_zzz, to_z_x_xxxy_xxxx, to_z_x_xxxy_xxxy, to_z_x_xxxy_xxxz, to_z_x_xxxy_xxyy, to_z_x_xxxy_xxyz, to_z_x_xxxy_xxzz, to_z_x_xxxy_xyyy, to_z_x_xxxy_xyyz, to_z_x_xxxy_xyzz, to_z_x_xxxy_xzzz, to_z_x_xxxy_yyyy, to_z_x_xxxy_yyyz, to_z_x_xxxy_yyzz, to_z_x_xxxy_yzzz, to_z_x_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxy_xxxx[k] = -8.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxxy[k] = -6.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxxz[k] = -6.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxyy[k] = -4.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxyz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxzz[k] = -4.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xyyy[k] = -2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xyyz[k] = -2.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xyzz[k] = -2.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xzzz[k] = -2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yyyy[k] = 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yyyz[k] = 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yyzz[k] = 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yzzz[k] = 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_zzzz[k] = 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1380-1395 components of targeted buffer : GG

        auto to_z_x_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 30);

        auto to_z_x_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 31);

        auto to_z_x_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 32);

        auto to_z_x_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 33);

        auto to_z_x_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 34);

        auto to_z_x_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 35);

        auto to_z_x_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 36);

        auto to_z_x_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 37);

        auto to_z_x_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 38);

        auto to_z_x_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 39);

        auto to_z_x_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 40);

        auto to_z_x_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 41);

        auto to_z_x_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 42);

        auto to_z_x_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 43);

        auto to_z_x_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxx, to_xxx_xxxxy, to_xxx_xxxxz, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_zzz, to_xxxzz_xxx, to_xxxzz_xxxxx, to_xxxzz_xxxxy, to_xxxzz_xxxxz, to_xxxzz_xxxyy, to_xxxzz_xxxyz, to_xxxzz_xxxzz, to_xxxzz_xxy, to_xxxzz_xxyyy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xxzzz, to_xxxzz_xyy, to_xxxzz_xyyyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_xzzzz, to_xxxzz_yyy, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_zzz, to_z_x_xxxz_xxxx, to_z_x_xxxz_xxxy, to_z_x_xxxz_xxxz, to_z_x_xxxz_xxyy, to_z_x_xxxz_xxyz, to_z_x_xxxz_xxzz, to_z_x_xxxz_xyyy, to_z_x_xxxz_xyyz, to_z_x_xxxz_xyzz, to_z_x_xxxz_xzzz, to_z_x_xxxz_yyyy, to_z_x_xxxz_yyyz, to_z_x_xxxz_yyzz, to_z_x_xxxz_yzzz, to_z_x_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxz_xxxx[k] = 4.0 * to_xxx_xxx[k] - 2.0 * to_xxx_xxxxx[k] * tke_0 - 8.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxxy[k] = 3.0 * to_xxx_xxy[k] - 2.0 * to_xxx_xxxxy[k] * tke_0 - 6.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxxz[k] = 3.0 * to_xxx_xxz[k] - 2.0 * to_xxx_xxxxz[k] * tke_0 - 6.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxyy[k] = 2.0 * to_xxx_xyy[k] - 2.0 * to_xxx_xxxyy[k] * tke_0 - 4.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxyz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xxxyz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxzz[k] = 2.0 * to_xxx_xzz[k] - 2.0 * to_xxx_xxxzz[k] * tke_0 - 4.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xyyy[k] = to_xxx_yyy[k] - 2.0 * to_xxx_xxyyy[k] * tke_0 - 2.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xyyz[k] = to_xxx_yyz[k] - 2.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xyzz[k] = to_xxx_yzz[k] - 2.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xzzz[k] = to_xxx_zzz[k] - 2.0 * to_xxx_xxzzz[k] * tke_0 - 2.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yyyy[k] = -2.0 * to_xxx_xyyyy[k] * tke_0 + 4.0 * to_xxxzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yyyz[k] = -2.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yyzz[k] = -2.0 * to_xxx_xyyzz[k] * tke_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yzzz[k] = -2.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_zzzz[k] = -2.0 * to_xxx_xzzzz[k] * tke_0 + 4.0 * to_xxxzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1395-1410 components of targeted buffer : GG

        auto to_z_x_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 45);

        auto to_z_x_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 46);

        auto to_z_x_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 47);

        auto to_z_x_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 48);

        auto to_z_x_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 49);

        auto to_z_x_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 50);

        auto to_z_x_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 51);

        auto to_z_x_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 52);

        auto to_z_x_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 53);

        auto to_z_x_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 54);

        auto to_z_x_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 55);

        auto to_z_x_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 56);

        auto to_z_x_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 57);

        auto to_z_x_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 58);

        auto to_z_x_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxx, to_xxyyz_xxxxy, to_xxyyz_xxxxz, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_zzz, to_z_x_xxyy_xxxx, to_z_x_xxyy_xxxy, to_z_x_xxyy_xxxz, to_z_x_xxyy_xxyy, to_z_x_xxyy_xxyz, to_z_x_xxyy_xxzz, to_z_x_xxyy_xyyy, to_z_x_xxyy_xyyz, to_z_x_xxyy_xyzz, to_z_x_xxyy_xzzz, to_z_x_xxyy_yyyy, to_z_x_xxyy_yyyz, to_z_x_xxyy_yyzz, to_z_x_xxyy_yzzz, to_z_x_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyy_xxxx[k] = -8.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxxy[k] = -6.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxxz[k] = -6.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxyy[k] = -4.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxyz[k] = -4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxzz[k] = -4.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xyyy[k] = -2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xyyz[k] = -2.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xyzz[k] = -2.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xzzz[k] = -2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yyyy[k] = 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yyyz[k] = 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yyzz[k] = 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yzzz[k] = 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_zzzz[k] = 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1410-1425 components of targeted buffer : GG

        auto to_z_x_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 60);

        auto to_z_x_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 61);

        auto to_z_x_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 62);

        auto to_z_x_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 63);

        auto to_z_x_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 64);

        auto to_z_x_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 65);

        auto to_z_x_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 66);

        auto to_z_x_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 67);

        auto to_z_x_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 68);

        auto to_z_x_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 69);

        auto to_z_x_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 70);

        auto to_z_x_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 71);

        auto to_z_x_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 72);

        auto to_z_x_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 73);

        auto to_z_x_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxx, to_xxy_xxxxy, to_xxy_xxxxz, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_zzz, to_xxyzz_xxx, to_xxyzz_xxxxx, to_xxyzz_xxxxy, to_xxyzz_xxxxz, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_zzz, to_z_x_xxyz_xxxx, to_z_x_xxyz_xxxy, to_z_x_xxyz_xxxz, to_z_x_xxyz_xxyy, to_z_x_xxyz_xxyz, to_z_x_xxyz_xxzz, to_z_x_xxyz_xyyy, to_z_x_xxyz_xyyz, to_z_x_xxyz_xyzz, to_z_x_xxyz_xzzz, to_z_x_xxyz_yyyy, to_z_x_xxyz_yyyz, to_z_x_xxyz_yyzz, to_z_x_xxyz_yzzz, to_z_x_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyz_xxxx[k] = 4.0 * to_xxy_xxx[k] - 2.0 * to_xxy_xxxxx[k] * tke_0 - 8.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxxy[k] = 3.0 * to_xxy_xxy[k] - 2.0 * to_xxy_xxxxy[k] * tke_0 - 6.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxxz[k] = 3.0 * to_xxy_xxz[k] - 2.0 * to_xxy_xxxxz[k] * tke_0 - 6.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxyy[k] = 2.0 * to_xxy_xyy[k] - 2.0 * to_xxy_xxxyy[k] * tke_0 - 4.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxyz[k] = 2.0 * to_xxy_xyz[k] - 2.0 * to_xxy_xxxyz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxzz[k] = 2.0 * to_xxy_xzz[k] - 2.0 * to_xxy_xxxzz[k] * tke_0 - 4.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xyyy[k] = to_xxy_yyy[k] - 2.0 * to_xxy_xxyyy[k] * tke_0 - 2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xyyz[k] = to_xxy_yyz[k] - 2.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xyzz[k] = to_xxy_yzz[k] - 2.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xzzz[k] = to_xxy_zzz[k] - 2.0 * to_xxy_xxzzz[k] * tke_0 - 2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yyyy[k] = -2.0 * to_xxy_xyyyy[k] * tke_0 + 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yyyz[k] = -2.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yyzz[k] = -2.0 * to_xxy_xyyzz[k] * tke_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yzzz[k] = -2.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_zzzz[k] = -2.0 * to_xxy_xzzzz[k] * tke_0 + 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1425-1440 components of targeted buffer : GG

        auto to_z_x_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 75);

        auto to_z_x_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 76);

        auto to_z_x_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 77);

        auto to_z_x_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 78);

        auto to_z_x_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 79);

        auto to_z_x_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 80);

        auto to_z_x_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 81);

        auto to_z_x_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 82);

        auto to_z_x_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 83);

        auto to_z_x_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 84);

        auto to_z_x_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 85);

        auto to_z_x_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 86);

        auto to_z_x_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 87);

        auto to_z_x_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 88);

        auto to_z_x_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxz_xxx, to_xxz_xxxxx, to_xxz_xxxxy, to_xxz_xxxxz, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_zzz, to_xxzzz_xxx, to_xxzzz_xxxxx, to_xxzzz_xxxxy, to_xxzzz_xxxxz, to_xxzzz_xxxyy, to_xxzzz_xxxyz, to_xxzzz_xxxzz, to_xxzzz_xxy, to_xxzzz_xxyyy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xxzzz, to_xxzzz_xyy, to_xxzzz_xyyyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_xzzzz, to_xxzzz_yyy, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_zzz, to_z_x_xxzz_xxxx, to_z_x_xxzz_xxxy, to_z_x_xxzz_xxxz, to_z_x_xxzz_xxyy, to_z_x_xxzz_xxyz, to_z_x_xxzz_xxzz, to_z_x_xxzz_xyyy, to_z_x_xxzz_xyyz, to_z_x_xxzz_xyzz, to_z_x_xxzz_xzzz, to_z_x_xxzz_yyyy, to_z_x_xxzz_yyyz, to_z_x_xxzz_yyzz, to_z_x_xxzz_yzzz, to_z_x_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxzz_xxxx[k] = 8.0 * to_xxz_xxx[k] - 4.0 * to_xxz_xxxxx[k] * tke_0 - 8.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxxy[k] = 6.0 * to_xxz_xxy[k] - 4.0 * to_xxz_xxxxy[k] * tke_0 - 6.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxxz[k] = 6.0 * to_xxz_xxz[k] - 4.0 * to_xxz_xxxxz[k] * tke_0 - 6.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxyy[k] = 4.0 * to_xxz_xyy[k] - 4.0 * to_xxz_xxxyy[k] * tke_0 - 4.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxyz[k] = 4.0 * to_xxz_xyz[k] - 4.0 * to_xxz_xxxyz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxzz[k] = 4.0 * to_xxz_xzz[k] - 4.0 * to_xxz_xxxzz[k] * tke_0 - 4.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xyyy[k] = 2.0 * to_xxz_yyy[k] - 4.0 * to_xxz_xxyyy[k] * tke_0 - 2.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xyyz[k] = 2.0 * to_xxz_yyz[k] - 4.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xyzz[k] = 2.0 * to_xxz_yzz[k] - 4.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xzzz[k] = 2.0 * to_xxz_zzz[k] - 4.0 * to_xxz_xxzzz[k] * tke_0 - 2.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yyyy[k] = -4.0 * to_xxz_xyyyy[k] * tke_0 + 4.0 * to_xxzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yyyz[k] = -4.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yyzz[k] = -4.0 * to_xxz_xyyzz[k] * tke_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yzzz[k] = -4.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_zzzz[k] = -4.0 * to_xxz_xzzzz[k] * tke_0 + 4.0 * to_xxzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1440-1455 components of targeted buffer : GG

        auto to_z_x_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 90);

        auto to_z_x_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 91);

        auto to_z_x_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 92);

        auto to_z_x_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 93);

        auto to_z_x_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 94);

        auto to_z_x_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 95);

        auto to_z_x_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 96);

        auto to_z_x_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 97);

        auto to_z_x_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 98);

        auto to_z_x_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 99);

        auto to_z_x_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 100);

        auto to_z_x_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 101);

        auto to_z_x_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 102);

        auto to_z_x_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 103);

        auto to_z_x_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxx, to_xyyyz_xxxxy, to_xyyyz_xxxxz, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_zzz, to_z_x_xyyy_xxxx, to_z_x_xyyy_xxxy, to_z_x_xyyy_xxxz, to_z_x_xyyy_xxyy, to_z_x_xyyy_xxyz, to_z_x_xyyy_xxzz, to_z_x_xyyy_xyyy, to_z_x_xyyy_xyyz, to_z_x_xyyy_xyzz, to_z_x_xyyy_xzzz, to_z_x_xyyy_yyyy, to_z_x_xyyy_yyyz, to_z_x_xyyy_yyzz, to_z_x_xyyy_yzzz, to_z_x_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyy_xxxx[k] = -8.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxxy[k] = -6.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxxz[k] = -6.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxyy[k] = -4.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxyz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxzz[k] = -4.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xyyy[k] = -2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xyyz[k] = -2.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xyzz[k] = -2.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xzzz[k] = -2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yyyy[k] = 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yyyz[k] = 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yyzz[k] = 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yzzz[k] = 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_zzzz[k] = 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1455-1470 components of targeted buffer : GG

        auto to_z_x_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 105);

        auto to_z_x_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 106);

        auto to_z_x_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 107);

        auto to_z_x_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 108);

        auto to_z_x_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 109);

        auto to_z_x_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 110);

        auto to_z_x_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 111);

        auto to_z_x_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 112);

        auto to_z_x_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 113);

        auto to_z_x_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 114);

        auto to_z_x_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 115);

        auto to_z_x_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 116);

        auto to_z_x_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 117);

        auto to_z_x_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 118);

        auto to_z_x_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxx, to_xyy_xxxxy, to_xyy_xxxxz, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_zzz, to_xyyzz_xxx, to_xyyzz_xxxxx, to_xyyzz_xxxxy, to_xyyzz_xxxxz, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_zzz, to_z_x_xyyz_xxxx, to_z_x_xyyz_xxxy, to_z_x_xyyz_xxxz, to_z_x_xyyz_xxyy, to_z_x_xyyz_xxyz, to_z_x_xyyz_xxzz, to_z_x_xyyz_xyyy, to_z_x_xyyz_xyyz, to_z_x_xyyz_xyzz, to_z_x_xyyz_xzzz, to_z_x_xyyz_yyyy, to_z_x_xyyz_yyyz, to_z_x_xyyz_yyzz, to_z_x_xyyz_yzzz, to_z_x_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyz_xxxx[k] = 4.0 * to_xyy_xxx[k] - 2.0 * to_xyy_xxxxx[k] * tke_0 - 8.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxxy[k] = 3.0 * to_xyy_xxy[k] - 2.0 * to_xyy_xxxxy[k] * tke_0 - 6.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxxz[k] = 3.0 * to_xyy_xxz[k] - 2.0 * to_xyy_xxxxz[k] * tke_0 - 6.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxyy[k] = 2.0 * to_xyy_xyy[k] - 2.0 * to_xyy_xxxyy[k] * tke_0 - 4.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxyz[k] = 2.0 * to_xyy_xyz[k] - 2.0 * to_xyy_xxxyz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxzz[k] = 2.0 * to_xyy_xzz[k] - 2.0 * to_xyy_xxxzz[k] * tke_0 - 4.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xyyy[k] = to_xyy_yyy[k] - 2.0 * to_xyy_xxyyy[k] * tke_0 - 2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xyyz[k] = to_xyy_yyz[k] - 2.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xyzz[k] = to_xyy_yzz[k] - 2.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xzzz[k] = to_xyy_zzz[k] - 2.0 * to_xyy_xxzzz[k] * tke_0 - 2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yyyy[k] = -2.0 * to_xyy_xyyyy[k] * tke_0 + 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yyyz[k] = -2.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yyzz[k] = -2.0 * to_xyy_xyyzz[k] * tke_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yzzz[k] = -2.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_zzzz[k] = -2.0 * to_xyy_xzzzz[k] * tke_0 + 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1470-1485 components of targeted buffer : GG

        auto to_z_x_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 120);

        auto to_z_x_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 121);

        auto to_z_x_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 122);

        auto to_z_x_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 123);

        auto to_z_x_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 124);

        auto to_z_x_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 125);

        auto to_z_x_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 126);

        auto to_z_x_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 127);

        auto to_z_x_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 128);

        auto to_z_x_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 129);

        auto to_z_x_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 130);

        auto to_z_x_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 131);

        auto to_z_x_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 132);

        auto to_z_x_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 133);

        auto to_z_x_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyz_xxx, to_xyz_xxxxx, to_xyz_xxxxy, to_xyz_xxxxz, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_zzz, to_xyzzz_xxx, to_xyzzz_xxxxx, to_xyzzz_xxxxy, to_xyzzz_xxxxz, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_zzz, to_z_x_xyzz_xxxx, to_z_x_xyzz_xxxy, to_z_x_xyzz_xxxz, to_z_x_xyzz_xxyy, to_z_x_xyzz_xxyz, to_z_x_xyzz_xxzz, to_z_x_xyzz_xyyy, to_z_x_xyzz_xyyz, to_z_x_xyzz_xyzz, to_z_x_xyzz_xzzz, to_z_x_xyzz_yyyy, to_z_x_xyzz_yyyz, to_z_x_xyzz_yyzz, to_z_x_xyzz_yzzz, to_z_x_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyzz_xxxx[k] = 8.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxxx[k] * tke_0 - 8.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxxy[k] = 6.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxxxy[k] * tke_0 - 6.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxxz[k] = 6.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxxxz[k] * tke_0 - 6.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxyy[k] = 4.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 4.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xxxyz[k] * tke_0 - 4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxzz[k] = 4.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 4.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xyyy[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xyyz[k] = 2.0 * to_xyz_yyz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xyzz[k] = 2.0 * to_xyz_yzz[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yyyy[k] = -4.0 * to_xyz_xyyyy[k] * tke_0 + 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yyyz[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yyzz[k] = -4.0 * to_xyz_xyyzz[k] * tke_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_zzzz[k] = -4.0 * to_xyz_xzzzz[k] * tke_0 + 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1485-1500 components of targeted buffer : GG

        auto to_z_x_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 135);

        auto to_z_x_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 136);

        auto to_z_x_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 137);

        auto to_z_x_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 138);

        auto to_z_x_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 139);

        auto to_z_x_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 140);

        auto to_z_x_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 141);

        auto to_z_x_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 142);

        auto to_z_x_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 143);

        auto to_z_x_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 144);

        auto to_z_x_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 145);

        auto to_z_x_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 146);

        auto to_z_x_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 147);

        auto to_z_x_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 148);

        auto to_z_x_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xzz_xxx, to_xzz_xxxxx, to_xzz_xxxxy, to_xzz_xxxxz, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_zzz, to_xzzzz_xxx, to_xzzzz_xxxxx, to_xzzzz_xxxxy, to_xzzzz_xxxxz, to_xzzzz_xxxyy, to_xzzzz_xxxyz, to_xzzzz_xxxzz, to_xzzzz_xxy, to_xzzzz_xxyyy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xxzzz, to_xzzzz_xyy, to_xzzzz_xyyyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_xzzzz, to_xzzzz_yyy, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_zzz, to_z_x_xzzz_xxxx, to_z_x_xzzz_xxxy, to_z_x_xzzz_xxxz, to_z_x_xzzz_xxyy, to_z_x_xzzz_xxyz, to_z_x_xzzz_xxzz, to_z_x_xzzz_xyyy, to_z_x_xzzz_xyyz, to_z_x_xzzz_xyzz, to_z_x_xzzz_xzzz, to_z_x_xzzz_yyyy, to_z_x_xzzz_yyyz, to_z_x_xzzz_yyzz, to_z_x_xzzz_yzzz, to_z_x_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzzz_xxxx[k] = 12.0 * to_xzz_xxx[k] - 6.0 * to_xzz_xxxxx[k] * tke_0 - 8.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxxy[k] = 9.0 * to_xzz_xxy[k] - 6.0 * to_xzz_xxxxy[k] * tke_0 - 6.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxxz[k] = 9.0 * to_xzz_xxz[k] - 6.0 * to_xzz_xxxxz[k] * tke_0 - 6.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxyy[k] = 6.0 * to_xzz_xyy[k] - 6.0 * to_xzz_xxxyy[k] * tke_0 - 4.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxyz[k] = 6.0 * to_xzz_xyz[k] - 6.0 * to_xzz_xxxyz[k] * tke_0 - 4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxzz[k] = 6.0 * to_xzz_xzz[k] - 6.0 * to_xzz_xxxzz[k] * tke_0 - 4.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xyyy[k] = 3.0 * to_xzz_yyy[k] - 6.0 * to_xzz_xxyyy[k] * tke_0 - 2.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xyyz[k] = 3.0 * to_xzz_yyz[k] - 6.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xyzz[k] = 3.0 * to_xzz_yzz[k] - 6.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xzzz[k] = 3.0 * to_xzz_zzz[k] - 6.0 * to_xzz_xxzzz[k] * tke_0 - 2.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yyyy[k] = -6.0 * to_xzz_xyyyy[k] * tke_0 + 4.0 * to_xzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yyyz[k] = -6.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yyzz[k] = -6.0 * to_xzz_xyyzz[k] * tke_0 + 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yzzz[k] = -6.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_zzzz[k] = -6.0 * to_xzz_xzzzz[k] * tke_0 + 4.0 * to_xzzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1500-1515 components of targeted buffer : GG

        auto to_z_x_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 150);

        auto to_z_x_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 151);

        auto to_z_x_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 152);

        auto to_z_x_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 153);

        auto to_z_x_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 154);

        auto to_z_x_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 155);

        auto to_z_x_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 156);

        auto to_z_x_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 157);

        auto to_z_x_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 158);

        auto to_z_x_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 159);

        auto to_z_x_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 160);

        auto to_z_x_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 161);

        auto to_z_x_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 162);

        auto to_z_x_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 163);

        auto to_z_x_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_yyyyz_xxx, to_yyyyz_xxxxx, to_yyyyz_xxxxy, to_yyyyz_xxxxz, to_yyyyz_xxxyy, to_yyyyz_xxxyz, to_yyyyz_xxxzz, to_yyyyz_xxy, to_yyyyz_xxyyy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xxzzz, to_yyyyz_xyy, to_yyyyz_xyyyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_xzzzz, to_yyyyz_yyy, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_zzz, to_z_x_yyyy_xxxx, to_z_x_yyyy_xxxy, to_z_x_yyyy_xxxz, to_z_x_yyyy_xxyy, to_z_x_yyyy_xxyz, to_z_x_yyyy_xxzz, to_z_x_yyyy_xyyy, to_z_x_yyyy_xyyz, to_z_x_yyyy_xyzz, to_z_x_yyyy_xzzz, to_z_x_yyyy_yyyy, to_z_x_yyyy_yyyz, to_z_x_yyyy_yyzz, to_z_x_yyyy_yzzz, to_z_x_yyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyy_xxxx[k] = -8.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxxy[k] = -6.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxxz[k] = -6.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxyy[k] = -4.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxyz[k] = -4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxzz[k] = -4.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xyyy[k] = -2.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xyyz[k] = -2.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xyzz[k] = -2.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xzzz[k] = -2.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yyyy[k] = 4.0 * to_yyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yyyz[k] = 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yyzz[k] = 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yzzz[k] = 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_zzzz[k] = 4.0 * to_yyyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1515-1530 components of targeted buffer : GG

        auto to_z_x_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 165);

        auto to_z_x_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 166);

        auto to_z_x_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 167);

        auto to_z_x_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 168);

        auto to_z_x_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 169);

        auto to_z_x_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 170);

        auto to_z_x_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 171);

        auto to_z_x_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 172);

        auto to_z_x_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 173);

        auto to_z_x_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 174);

        auto to_z_x_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 175);

        auto to_z_x_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 176);

        auto to_z_x_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 177);

        auto to_z_x_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 178);

        auto to_z_x_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_yyy_xxx, to_yyy_xxxxx, to_yyy_xxxxy, to_yyy_xxxxz, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_zzz, to_yyyzz_xxx, to_yyyzz_xxxxx, to_yyyzz_xxxxy, to_yyyzz_xxxxz, to_yyyzz_xxxyy, to_yyyzz_xxxyz, to_yyyzz_xxxzz, to_yyyzz_xxy, to_yyyzz_xxyyy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xxzzz, to_yyyzz_xyy, to_yyyzz_xyyyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_xzzzz, to_yyyzz_yyy, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_zzz, to_z_x_yyyz_xxxx, to_z_x_yyyz_xxxy, to_z_x_yyyz_xxxz, to_z_x_yyyz_xxyy, to_z_x_yyyz_xxyz, to_z_x_yyyz_xxzz, to_z_x_yyyz_xyyy, to_z_x_yyyz_xyyz, to_z_x_yyyz_xyzz, to_z_x_yyyz_xzzz, to_z_x_yyyz_yyyy, to_z_x_yyyz_yyyz, to_z_x_yyyz_yyzz, to_z_x_yyyz_yzzz, to_z_x_yyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyz_xxxx[k] = 4.0 * to_yyy_xxx[k] - 2.0 * to_yyy_xxxxx[k] * tke_0 - 8.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxxy[k] = 3.0 * to_yyy_xxy[k] - 2.0 * to_yyy_xxxxy[k] * tke_0 - 6.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxxz[k] = 3.0 * to_yyy_xxz[k] - 2.0 * to_yyy_xxxxz[k] * tke_0 - 6.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxyy[k] = 2.0 * to_yyy_xyy[k] - 2.0 * to_yyy_xxxyy[k] * tke_0 - 4.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxyz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xxxyz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxzz[k] = 2.0 * to_yyy_xzz[k] - 2.0 * to_yyy_xxxzz[k] * tke_0 - 4.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xyyy[k] = to_yyy_yyy[k] - 2.0 * to_yyy_xxyyy[k] * tke_0 - 2.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xyyz[k] = to_yyy_yyz[k] - 2.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xyzz[k] = to_yyy_yzz[k] - 2.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xzzz[k] = to_yyy_zzz[k] - 2.0 * to_yyy_xxzzz[k] * tke_0 - 2.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yyyy[k] = -2.0 * to_yyy_xyyyy[k] * tke_0 + 4.0 * to_yyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yyyz[k] = -2.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yyzz[k] = -2.0 * to_yyy_xyyzz[k] * tke_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yzzz[k] = -2.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_zzzz[k] = -2.0 * to_yyy_xzzzz[k] * tke_0 + 4.0 * to_yyyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1530-1545 components of targeted buffer : GG

        auto to_z_x_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 180);

        auto to_z_x_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 181);

        auto to_z_x_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 182);

        auto to_z_x_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 183);

        auto to_z_x_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 184);

        auto to_z_x_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 185);

        auto to_z_x_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 186);

        auto to_z_x_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 187);

        auto to_z_x_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 188);

        auto to_z_x_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 189);

        auto to_z_x_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 190);

        auto to_z_x_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 191);

        auto to_z_x_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 192);

        auto to_z_x_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 193);

        auto to_z_x_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_yyz_xxx, to_yyz_xxxxx, to_yyz_xxxxy, to_yyz_xxxxz, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_zzz, to_yyzzz_xxx, to_yyzzz_xxxxx, to_yyzzz_xxxxy, to_yyzzz_xxxxz, to_yyzzz_xxxyy, to_yyzzz_xxxyz, to_yyzzz_xxxzz, to_yyzzz_xxy, to_yyzzz_xxyyy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xxzzz, to_yyzzz_xyy, to_yyzzz_xyyyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_xzzzz, to_yyzzz_yyy, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_zzz, to_z_x_yyzz_xxxx, to_z_x_yyzz_xxxy, to_z_x_yyzz_xxxz, to_z_x_yyzz_xxyy, to_z_x_yyzz_xxyz, to_z_x_yyzz_xxzz, to_z_x_yyzz_xyyy, to_z_x_yyzz_xyyz, to_z_x_yyzz_xyzz, to_z_x_yyzz_xzzz, to_z_x_yyzz_yyyy, to_z_x_yyzz_yyyz, to_z_x_yyzz_yyzz, to_z_x_yyzz_yzzz, to_z_x_yyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyzz_xxxx[k] = 8.0 * to_yyz_xxx[k] - 4.0 * to_yyz_xxxxx[k] * tke_0 - 8.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxxy[k] = 6.0 * to_yyz_xxy[k] - 4.0 * to_yyz_xxxxy[k] * tke_0 - 6.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxxz[k] = 6.0 * to_yyz_xxz[k] - 4.0 * to_yyz_xxxxz[k] * tke_0 - 6.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxyy[k] = 4.0 * to_yyz_xyy[k] - 4.0 * to_yyz_xxxyy[k] * tke_0 - 4.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxyz[k] = 4.0 * to_yyz_xyz[k] - 4.0 * to_yyz_xxxyz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxzz[k] = 4.0 * to_yyz_xzz[k] - 4.0 * to_yyz_xxxzz[k] * tke_0 - 4.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xyyy[k] = 2.0 * to_yyz_yyy[k] - 4.0 * to_yyz_xxyyy[k] * tke_0 - 2.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xyyz[k] = 2.0 * to_yyz_yyz[k] - 4.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xyzz[k] = 2.0 * to_yyz_yzz[k] - 4.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xzzz[k] = 2.0 * to_yyz_zzz[k] - 4.0 * to_yyz_xxzzz[k] * tke_0 - 2.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yyyy[k] = -4.0 * to_yyz_xyyyy[k] * tke_0 + 4.0 * to_yyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yyyz[k] = -4.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yyzz[k] = -4.0 * to_yyz_xyyzz[k] * tke_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yzzz[k] = -4.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_zzzz[k] = -4.0 * to_yyz_xzzzz[k] * tke_0 + 4.0 * to_yyzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1545-1560 components of targeted buffer : GG

        auto to_z_x_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 195);

        auto to_z_x_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 196);

        auto to_z_x_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 197);

        auto to_z_x_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 198);

        auto to_z_x_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 199);

        auto to_z_x_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 200);

        auto to_z_x_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 201);

        auto to_z_x_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 202);

        auto to_z_x_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 203);

        auto to_z_x_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 204);

        auto to_z_x_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 205);

        auto to_z_x_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 206);

        auto to_z_x_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 207);

        auto to_z_x_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 208);

        auto to_z_x_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_yzz_xxx, to_yzz_xxxxx, to_yzz_xxxxy, to_yzz_xxxxz, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_zzz, to_yzzzz_xxx, to_yzzzz_xxxxx, to_yzzzz_xxxxy, to_yzzzz_xxxxz, to_yzzzz_xxxyy, to_yzzzz_xxxyz, to_yzzzz_xxxzz, to_yzzzz_xxy, to_yzzzz_xxyyy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xxzzz, to_yzzzz_xyy, to_yzzzz_xyyyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_xzzzz, to_yzzzz_yyy, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_zzz, to_z_x_yzzz_xxxx, to_z_x_yzzz_xxxy, to_z_x_yzzz_xxxz, to_z_x_yzzz_xxyy, to_z_x_yzzz_xxyz, to_z_x_yzzz_xxzz, to_z_x_yzzz_xyyy, to_z_x_yzzz_xyyz, to_z_x_yzzz_xyzz, to_z_x_yzzz_xzzz, to_z_x_yzzz_yyyy, to_z_x_yzzz_yyyz, to_z_x_yzzz_yyzz, to_z_x_yzzz_yzzz, to_z_x_yzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzzz_xxxx[k] = 12.0 * to_yzz_xxx[k] - 6.0 * to_yzz_xxxxx[k] * tke_0 - 8.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxxy[k] = 9.0 * to_yzz_xxy[k] - 6.0 * to_yzz_xxxxy[k] * tke_0 - 6.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxxz[k] = 9.0 * to_yzz_xxz[k] - 6.0 * to_yzz_xxxxz[k] * tke_0 - 6.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxyy[k] = 6.0 * to_yzz_xyy[k] - 6.0 * to_yzz_xxxyy[k] * tke_0 - 4.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxyz[k] = 6.0 * to_yzz_xyz[k] - 6.0 * to_yzz_xxxyz[k] * tke_0 - 4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxzz[k] = 6.0 * to_yzz_xzz[k] - 6.0 * to_yzz_xxxzz[k] * tke_0 - 4.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xyyy[k] = 3.0 * to_yzz_yyy[k] - 6.0 * to_yzz_xxyyy[k] * tke_0 - 2.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xyyz[k] = 3.0 * to_yzz_yyz[k] - 6.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xyzz[k] = 3.0 * to_yzz_yzz[k] - 6.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xzzz[k] = 3.0 * to_yzz_zzz[k] - 6.0 * to_yzz_xxzzz[k] * tke_0 - 2.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yyyy[k] = -6.0 * to_yzz_xyyyy[k] * tke_0 + 4.0 * to_yzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yyyz[k] = -6.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yyzz[k] = -6.0 * to_yzz_xyyzz[k] * tke_0 + 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yzzz[k] = -6.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_zzzz[k] = -6.0 * to_yzz_xzzzz[k] * tke_0 + 4.0 * to_yzzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1560-1575 components of targeted buffer : GG

        auto to_z_x_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 210);

        auto to_z_x_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 211);

        auto to_z_x_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 212);

        auto to_z_x_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 213);

        auto to_z_x_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 214);

        auto to_z_x_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 215);

        auto to_z_x_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 216);

        auto to_z_x_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 217);

        auto to_z_x_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 218);

        auto to_z_x_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 219);

        auto to_z_x_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 220);

        auto to_z_x_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 221);

        auto to_z_x_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 222);

        auto to_z_x_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 223);

        auto to_z_x_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 6 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_z_x_zzzz_xxxx, to_z_x_zzzz_xxxy, to_z_x_zzzz_xxxz, to_z_x_zzzz_xxyy, to_z_x_zzzz_xxyz, to_z_x_zzzz_xxzz, to_z_x_zzzz_xyyy, to_z_x_zzzz_xyyz, to_z_x_zzzz_xyzz, to_z_x_zzzz_xzzz, to_z_x_zzzz_yyyy, to_z_x_zzzz_yyyz, to_z_x_zzzz_yyzz, to_z_x_zzzz_yzzz, to_z_x_zzzz_zzzz, to_zzz_xxx, to_zzz_xxxxx, to_zzz_xxxxy, to_zzz_xxxxz, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_zzz, to_zzzzz_xxx, to_zzzzz_xxxxx, to_zzzzz_xxxxy, to_zzzzz_xxxxz, to_zzzzz_xxxyy, to_zzzzz_xxxyz, to_zzzzz_xxxzz, to_zzzzz_xxy, to_zzzzz_xxyyy, to_zzzzz_xxyyz, to_zzzzz_xxyzz, to_zzzzz_xxz, to_zzzzz_xxzzz, to_zzzzz_xyy, to_zzzzz_xyyyy, to_zzzzz_xyyyz, to_zzzzz_xyyzz, to_zzzzz_xyz, to_zzzzz_xyzzz, to_zzzzz_xzz, to_zzzzz_xzzzz, to_zzzzz_yyy, to_zzzzz_yyz, to_zzzzz_yzz, to_zzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzzz_xxxx[k] = 16.0 * to_zzz_xxx[k] - 8.0 * to_zzz_xxxxx[k] * tke_0 - 8.0 * to_zzzzz_xxx[k] * tbe_0 + 4.0 * to_zzzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxxy[k] = 12.0 * to_zzz_xxy[k] - 8.0 * to_zzz_xxxxy[k] * tke_0 - 6.0 * to_zzzzz_xxy[k] * tbe_0 + 4.0 * to_zzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxxz[k] = 12.0 * to_zzz_xxz[k] - 8.0 * to_zzz_xxxxz[k] * tke_0 - 6.0 * to_zzzzz_xxz[k] * tbe_0 + 4.0 * to_zzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxyy[k] = 8.0 * to_zzz_xyy[k] - 8.0 * to_zzz_xxxyy[k] * tke_0 - 4.0 * to_zzzzz_xyy[k] * tbe_0 + 4.0 * to_zzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxyz[k] = 8.0 * to_zzz_xyz[k] - 8.0 * to_zzz_xxxyz[k] * tke_0 - 4.0 * to_zzzzz_xyz[k] * tbe_0 + 4.0 * to_zzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxzz[k] = 8.0 * to_zzz_xzz[k] - 8.0 * to_zzz_xxxzz[k] * tke_0 - 4.0 * to_zzzzz_xzz[k] * tbe_0 + 4.0 * to_zzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xyyy[k] = 4.0 * to_zzz_yyy[k] - 8.0 * to_zzz_xxyyy[k] * tke_0 - 2.0 * to_zzzzz_yyy[k] * tbe_0 + 4.0 * to_zzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xyyz[k] = 4.0 * to_zzz_yyz[k] - 8.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_zzzzz_yyz[k] * tbe_0 + 4.0 * to_zzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xyzz[k] = 4.0 * to_zzz_yzz[k] - 8.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_zzzzz_yzz[k] * tbe_0 + 4.0 * to_zzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xzzz[k] = 4.0 * to_zzz_zzz[k] - 8.0 * to_zzz_xxzzz[k] * tke_0 - 2.0 * to_zzzzz_zzz[k] * tbe_0 + 4.0 * to_zzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yyyy[k] = -8.0 * to_zzz_xyyyy[k] * tke_0 + 4.0 * to_zzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yyyz[k] = -8.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_zzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yyzz[k] = -8.0 * to_zzz_xyyzz[k] * tke_0 + 4.0 * to_zzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yzzz[k] = -8.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_zzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_zzzz[k] = -8.0 * to_zzz_xzzzz[k] * tke_0 + 4.0 * to_zzzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1575-1590 components of targeted buffer : GG

        auto to_z_y_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 0);

        auto to_z_y_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 1);

        auto to_z_y_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 2);

        auto to_z_y_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 3);

        auto to_z_y_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 4);

        auto to_z_y_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 5);

        auto to_z_y_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 6);

        auto to_z_y_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 7);

        auto to_z_y_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 8);

        auto to_z_y_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 9);

        auto to_z_y_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 10);

        auto to_z_y_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 11);

        auto to_z_y_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 12);

        auto to_z_y_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 13);

        auto to_z_y_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxz_xxx, to_xxxxz_xxxxy, to_xxxxz_xxxyy, to_xxxxz_xxxyz, to_xxxxz_xxy, to_xxxxz_xxyyy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xyy, to_xxxxz_xyyyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_yyy, to_xxxxz_yyyyy, to_xxxxz_yyyyz, to_xxxxz_yyyzz, to_xxxxz_yyz, to_xxxxz_yyzzz, to_xxxxz_yzz, to_xxxxz_yzzzz, to_xxxxz_zzz, to_z_y_xxxx_xxxx, to_z_y_xxxx_xxxy, to_z_y_xxxx_xxxz, to_z_y_xxxx_xxyy, to_z_y_xxxx_xxyz, to_z_y_xxxx_xxzz, to_z_y_xxxx_xyyy, to_z_y_xxxx_xyyz, to_z_y_xxxx_xyzz, to_z_y_xxxx_xzzz, to_z_y_xxxx_yyyy, to_z_y_xxxx_yyyz, to_z_y_xxxx_yyzz, to_z_y_xxxx_yzzz, to_z_y_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxx_xxxx[k] = 4.0 * to_xxxxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxxy[k] = -2.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxxz[k] = 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxyy[k] = -4.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxyz[k] = -2.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxzz[k] = 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xyyy[k] = -6.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xyyz[k] = -4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xyzz[k] = -2.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xzzz[k] = 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yyyy[k] = -8.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yyyz[k] = -6.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yyzz[k] = -4.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yzzz[k] = -2.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_zzzz[k] = 4.0 * to_xxxxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1590-1605 components of targeted buffer : GG

        auto to_z_y_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 15);

        auto to_z_y_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 16);

        auto to_z_y_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 17);

        auto to_z_y_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 18);

        auto to_z_y_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 19);

        auto to_z_y_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 20);

        auto to_z_y_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 21);

        auto to_z_y_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 22);

        auto to_z_y_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 23);

        auto to_z_y_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 24);

        auto to_z_y_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 25);

        auto to_z_y_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 26);

        auto to_z_y_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 27);

        auto to_z_y_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 28);

        auto to_z_y_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxy, to_xxxyz_xxxyy, to_xxxyz_xxxyz, to_xxxyz_xxy, to_xxxyz_xxyyy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyyyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_yyy, to_xxxyz_yyyyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_z_y_xxxy_xxxx, to_z_y_xxxy_xxxy, to_z_y_xxxy_xxxz, to_z_y_xxxy_xxyy, to_z_y_xxxy_xxyz, to_z_y_xxxy_xxzz, to_z_y_xxxy_xyyy, to_z_y_xxxy_xyyz, to_z_y_xxxy_xyzz, to_z_y_xxxy_xzzz, to_z_y_xxxy_yyyy, to_z_y_xxxy_yyyz, to_z_y_xxxy_yyzz, to_z_y_xxxy_yzzz, to_z_y_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxy_xxxx[k] = 4.0 * to_xxxyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxxy[k] = -2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxxz[k] = 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxyy[k] = -4.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxyz[k] = -2.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxzz[k] = 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xyyy[k] = -6.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xyyz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xyzz[k] = -2.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xzzz[k] = 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yyyy[k] = -8.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yyyz[k] = -6.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yyzz[k] = -4.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yzzz[k] = -2.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_zzzz[k] = 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1605-1620 components of targeted buffer : GG

        auto to_z_y_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 30);

        auto to_z_y_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 31);

        auto to_z_y_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 32);

        auto to_z_y_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 33);

        auto to_z_y_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 34);

        auto to_z_y_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 35);

        auto to_z_y_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 36);

        auto to_z_y_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 37);

        auto to_z_y_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 38);

        auto to_z_y_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 39);

        auto to_z_y_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 40);

        auto to_z_y_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 41);

        auto to_z_y_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 42);

        auto to_z_y_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 43);

        auto to_z_y_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxy, to_xxx_xxxyy, to_xxx_xxxyz, to_xxx_xxy, to_xxx_xxyyy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xyy, to_xxx_xyyyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_yyy, to_xxx_yyyyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxxzz_xxx, to_xxxzz_xxxxy, to_xxxzz_xxxyy, to_xxxzz_xxxyz, to_xxxzz_xxy, to_xxxzz_xxyyy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xyy, to_xxxzz_xyyyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_yyy, to_xxxzz_yyyyy, to_xxxzz_yyyyz, to_xxxzz_yyyzz, to_xxxzz_yyz, to_xxxzz_yyzzz, to_xxxzz_yzz, to_xxxzz_yzzzz, to_xxxzz_zzz, to_z_y_xxxz_xxxx, to_z_y_xxxz_xxxy, to_z_y_xxxz_xxxz, to_z_y_xxxz_xxyy, to_z_y_xxxz_xxyz, to_z_y_xxxz_xxzz, to_z_y_xxxz_xyyy, to_z_y_xxxz_xyyz, to_z_y_xxxz_xyzz, to_z_y_xxxz_xzzz, to_z_y_xxxz_yyyy, to_z_y_xxxz_yyyz, to_z_y_xxxz_yyzz, to_z_y_xxxz_yzzz, to_z_y_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxz_xxxx[k] = -2.0 * to_xxx_xxxxy[k] * tke_0 + 4.0 * to_xxxzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxxy[k] = to_xxx_xxx[k] - 2.0 * to_xxx_xxxyy[k] * tke_0 - 2.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxxz[k] = -2.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxyy[k] = 2.0 * to_xxx_xxy[k] - 2.0 * to_xxx_xxyyy[k] * tke_0 - 4.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxyz[k] = to_xxx_xxz[k] - 2.0 * to_xxx_xxyyz[k] * tke_0 - 2.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxzz[k] = -2.0 * to_xxx_xxyzz[k] * tke_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xyyy[k] = 3.0 * to_xxx_xyy[k] - 2.0 * to_xxx_xyyyy[k] * tke_0 - 6.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xyyz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xyyyz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xyzz[k] = to_xxx_xzz[k] - 2.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xzzz[k] = -2.0 * to_xxx_xyzzz[k] * tke_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yyyy[k] = 4.0 * to_xxx_yyy[k] - 2.0 * to_xxx_yyyyy[k] * tke_0 - 8.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yyyz[k] = 3.0 * to_xxx_yyz[k] - 2.0 * to_xxx_yyyyz[k] * tke_0 - 6.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yyzz[k] = 2.0 * to_xxx_yzz[k] - 2.0 * to_xxx_yyyzz[k] * tke_0 - 4.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yzzz[k] = to_xxx_zzz[k] - 2.0 * to_xxx_yyzzz[k] * tke_0 - 2.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_zzzz[k] = -2.0 * to_xxx_yzzzz[k] * tke_0 + 4.0 * to_xxxzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1620-1635 components of targeted buffer : GG

        auto to_z_y_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 45);

        auto to_z_y_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 46);

        auto to_z_y_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 47);

        auto to_z_y_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 48);

        auto to_z_y_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 49);

        auto to_z_y_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 50);

        auto to_z_y_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 51);

        auto to_z_y_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 52);

        auto to_z_y_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 53);

        auto to_z_y_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 54);

        auto to_z_y_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 55);

        auto to_z_y_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 56);

        auto to_z_y_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 57);

        auto to_z_y_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 58);

        auto to_z_y_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxy, to_xxyyz_xxxyy, to_xxyyz_xxxyz, to_xxyyz_xxy, to_xxyyz_xxyyy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyyyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_yyy, to_xxyyz_yyyyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_z_y_xxyy_xxxx, to_z_y_xxyy_xxxy, to_z_y_xxyy_xxxz, to_z_y_xxyy_xxyy, to_z_y_xxyy_xxyz, to_z_y_xxyy_xxzz, to_z_y_xxyy_xyyy, to_z_y_xxyy_xyyz, to_z_y_xxyy_xyzz, to_z_y_xxyy_xzzz, to_z_y_xxyy_yyyy, to_z_y_xxyy_yyyz, to_z_y_xxyy_yyzz, to_z_y_xxyy_yzzz, to_z_y_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyy_xxxx[k] = 4.0 * to_xxyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxxy[k] = -2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxxz[k] = 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxyy[k] = -4.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxyz[k] = -2.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxzz[k] = 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xyyy[k] = -6.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xyyz[k] = -4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xyzz[k] = -2.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xzzz[k] = 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yyyy[k] = -8.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yyyz[k] = -6.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yyzz[k] = -4.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yzzz[k] = -2.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_zzzz[k] = 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1635-1650 components of targeted buffer : GG

        auto to_z_y_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 60);

        auto to_z_y_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 61);

        auto to_z_y_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 62);

        auto to_z_y_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 63);

        auto to_z_y_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 64);

        auto to_z_y_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 65);

        auto to_z_y_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 66);

        auto to_z_y_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 67);

        auto to_z_y_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 68);

        auto to_z_y_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 69);

        auto to_z_y_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 70);

        auto to_z_y_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 71);

        auto to_z_y_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 72);

        auto to_z_y_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 73);

        auto to_z_y_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxy, to_xxy_xxxyy, to_xxy_xxxyz, to_xxy_xxy, to_xxy_xxyyy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xyy, to_xxy_xyyyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_yyy, to_xxy_yyyyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxyzz_xxx, to_xxyzz_xxxxy, to_xxyzz_xxxyy, to_xxyzz_xxxyz, to_xxyzz_xxy, to_xxyzz_xxyyy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyyyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_yyy, to_xxyzz_yyyyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_z_y_xxyz_xxxx, to_z_y_xxyz_xxxy, to_z_y_xxyz_xxxz, to_z_y_xxyz_xxyy, to_z_y_xxyz_xxyz, to_z_y_xxyz_xxzz, to_z_y_xxyz_xyyy, to_z_y_xxyz_xyyz, to_z_y_xxyz_xyzz, to_z_y_xxyz_xzzz, to_z_y_xxyz_yyyy, to_z_y_xxyz_yyyz, to_z_y_xxyz_yyzz, to_z_y_xxyz_yzzz, to_z_y_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyz_xxxx[k] = -2.0 * to_xxy_xxxxy[k] * tke_0 + 4.0 * to_xxyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxxy[k] = to_xxy_xxx[k] - 2.0 * to_xxy_xxxyy[k] * tke_0 - 2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxxz[k] = -2.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxyy[k] = 2.0 * to_xxy_xxy[k] - 2.0 * to_xxy_xxyyy[k] * tke_0 - 4.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxyz[k] = to_xxy_xxz[k] - 2.0 * to_xxy_xxyyz[k] * tke_0 - 2.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxzz[k] = -2.0 * to_xxy_xxyzz[k] * tke_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xyyy[k] = 3.0 * to_xxy_xyy[k] - 2.0 * to_xxy_xyyyy[k] * tke_0 - 6.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xyyz[k] = 2.0 * to_xxy_xyz[k] - 2.0 * to_xxy_xyyyz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xyzz[k] = to_xxy_xzz[k] - 2.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xzzz[k] = -2.0 * to_xxy_xyzzz[k] * tke_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yyyy[k] = 4.0 * to_xxy_yyy[k] - 2.0 * to_xxy_yyyyy[k] * tke_0 - 8.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yyyz[k] = 3.0 * to_xxy_yyz[k] - 2.0 * to_xxy_yyyyz[k] * tke_0 - 6.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yyzz[k] = 2.0 * to_xxy_yzz[k] - 2.0 * to_xxy_yyyzz[k] * tke_0 - 4.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yzzz[k] = to_xxy_zzz[k] - 2.0 * to_xxy_yyzzz[k] * tke_0 - 2.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_zzzz[k] = -2.0 * to_xxy_yzzzz[k] * tke_0 + 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1650-1665 components of targeted buffer : GG

        auto to_z_y_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 75);

        auto to_z_y_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 76);

        auto to_z_y_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 77);

        auto to_z_y_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 78);

        auto to_z_y_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 79);

        auto to_z_y_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 80);

        auto to_z_y_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 81);

        auto to_z_y_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 82);

        auto to_z_y_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 83);

        auto to_z_y_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 84);

        auto to_z_y_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 85);

        auto to_z_y_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 86);

        auto to_z_y_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 87);

        auto to_z_y_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 88);

        auto to_z_y_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxz_xxx, to_xxz_xxxxy, to_xxz_xxxyy, to_xxz_xxxyz, to_xxz_xxy, to_xxz_xxyyy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xyy, to_xxz_xyyyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_yyy, to_xxz_yyyyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_xxzzz_xxx, to_xxzzz_xxxxy, to_xxzzz_xxxyy, to_xxzzz_xxxyz, to_xxzzz_xxy, to_xxzzz_xxyyy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xyy, to_xxzzz_xyyyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_yyy, to_xxzzz_yyyyy, to_xxzzz_yyyyz, to_xxzzz_yyyzz, to_xxzzz_yyz, to_xxzzz_yyzzz, to_xxzzz_yzz, to_xxzzz_yzzzz, to_xxzzz_zzz, to_z_y_xxzz_xxxx, to_z_y_xxzz_xxxy, to_z_y_xxzz_xxxz, to_z_y_xxzz_xxyy, to_z_y_xxzz_xxyz, to_z_y_xxzz_xxzz, to_z_y_xxzz_xyyy, to_z_y_xxzz_xyyz, to_z_y_xxzz_xyzz, to_z_y_xxzz_xzzz, to_z_y_xxzz_yyyy, to_z_y_xxzz_yyyz, to_z_y_xxzz_yyzz, to_z_y_xxzz_yzzz, to_z_y_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxzz_xxxx[k] = -4.0 * to_xxz_xxxxy[k] * tke_0 + 4.0 * to_xxzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxxy[k] = 2.0 * to_xxz_xxx[k] - 4.0 * to_xxz_xxxyy[k] * tke_0 - 2.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxxz[k] = -4.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxyy[k] = 4.0 * to_xxz_xxy[k] - 4.0 * to_xxz_xxyyy[k] * tke_0 - 4.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxyz[k] = 2.0 * to_xxz_xxz[k] - 4.0 * to_xxz_xxyyz[k] * tke_0 - 2.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxzz[k] = -4.0 * to_xxz_xxyzz[k] * tke_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xyyy[k] = 6.0 * to_xxz_xyy[k] - 4.0 * to_xxz_xyyyy[k] * tke_0 - 6.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xyyz[k] = 4.0 * to_xxz_xyz[k] - 4.0 * to_xxz_xyyyz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xyzz[k] = 2.0 * to_xxz_xzz[k] - 4.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xzzz[k] = -4.0 * to_xxz_xyzzz[k] * tke_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yyyy[k] = 8.0 * to_xxz_yyy[k] - 4.0 * to_xxz_yyyyy[k] * tke_0 - 8.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yyyz[k] = 6.0 * to_xxz_yyz[k] - 4.0 * to_xxz_yyyyz[k] * tke_0 - 6.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yyzz[k] = 4.0 * to_xxz_yzz[k] - 4.0 * to_xxz_yyyzz[k] * tke_0 - 4.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yzzz[k] = 2.0 * to_xxz_zzz[k] - 4.0 * to_xxz_yyzzz[k] * tke_0 - 2.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_zzzz[k] = -4.0 * to_xxz_yzzzz[k] * tke_0 + 4.0 * to_xxzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1665-1680 components of targeted buffer : GG

        auto to_z_y_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 90);

        auto to_z_y_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 91);

        auto to_z_y_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 92);

        auto to_z_y_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 93);

        auto to_z_y_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 94);

        auto to_z_y_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 95);

        auto to_z_y_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 96);

        auto to_z_y_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 97);

        auto to_z_y_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 98);

        auto to_z_y_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 99);

        auto to_z_y_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 100);

        auto to_z_y_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 101);

        auto to_z_y_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 102);

        auto to_z_y_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 103);

        auto to_z_y_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxy, to_xyyyz_xxxyy, to_xyyyz_xxxyz, to_xyyyz_xxy, to_xyyyz_xxyyy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyyyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_yyy, to_xyyyz_yyyyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, to_z_y_xyyy_xxxx, to_z_y_xyyy_xxxy, to_z_y_xyyy_xxxz, to_z_y_xyyy_xxyy, to_z_y_xyyy_xxyz, to_z_y_xyyy_xxzz, to_z_y_xyyy_xyyy, to_z_y_xyyy_xyyz, to_z_y_xyyy_xyzz, to_z_y_xyyy_xzzz, to_z_y_xyyy_yyyy, to_z_y_xyyy_yyyz, to_z_y_xyyy_yyzz, to_z_y_xyyy_yzzz, to_z_y_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyy_xxxx[k] = 4.0 * to_xyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxxy[k] = -2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxxz[k] = 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxyy[k] = -4.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxyz[k] = -2.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxzz[k] = 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xyyy[k] = -6.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xyyz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xyzz[k] = -2.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xzzz[k] = 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yyyy[k] = -8.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yyyz[k] = -6.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yyzz[k] = -4.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yzzz[k] = -2.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_zzzz[k] = 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1680-1695 components of targeted buffer : GG

        auto to_z_y_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 105);

        auto to_z_y_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 106);

        auto to_z_y_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 107);

        auto to_z_y_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 108);

        auto to_z_y_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 109);

        auto to_z_y_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 110);

        auto to_z_y_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 111);

        auto to_z_y_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 112);

        auto to_z_y_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 113);

        auto to_z_y_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 114);

        auto to_z_y_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 115);

        auto to_z_y_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 116);

        auto to_z_y_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 117);

        auto to_z_y_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 118);

        auto to_z_y_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxy, to_xyy_xxxyy, to_xyy_xxxyz, to_xyy_xxy, to_xyy_xxyyy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xyy, to_xyy_xyyyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_yyy, to_xyy_yyyyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyyzz_xxx, to_xyyzz_xxxxy, to_xyyzz_xxxyy, to_xyyzz_xxxyz, to_xyyzz_xxy, to_xyyzz_xxyyy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyyyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_yyy, to_xyyzz_yyyyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, to_z_y_xyyz_xxxx, to_z_y_xyyz_xxxy, to_z_y_xyyz_xxxz, to_z_y_xyyz_xxyy, to_z_y_xyyz_xxyz, to_z_y_xyyz_xxzz, to_z_y_xyyz_xyyy, to_z_y_xyyz_xyyz, to_z_y_xyyz_xyzz, to_z_y_xyyz_xzzz, to_z_y_xyyz_yyyy, to_z_y_xyyz_yyyz, to_z_y_xyyz_yyzz, to_z_y_xyyz_yzzz, to_z_y_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyz_xxxx[k] = -2.0 * to_xyy_xxxxy[k] * tke_0 + 4.0 * to_xyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxxy[k] = to_xyy_xxx[k] - 2.0 * to_xyy_xxxyy[k] * tke_0 - 2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxxz[k] = -2.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxyy[k] = 2.0 * to_xyy_xxy[k] - 2.0 * to_xyy_xxyyy[k] * tke_0 - 4.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxyz[k] = to_xyy_xxz[k] - 2.0 * to_xyy_xxyyz[k] * tke_0 - 2.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxzz[k] = -2.0 * to_xyy_xxyzz[k] * tke_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xyyy[k] = 3.0 * to_xyy_xyy[k] - 2.0 * to_xyy_xyyyy[k] * tke_0 - 6.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xyyz[k] = 2.0 * to_xyy_xyz[k] - 2.0 * to_xyy_xyyyz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xyzz[k] = to_xyy_xzz[k] - 2.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xzzz[k] = -2.0 * to_xyy_xyzzz[k] * tke_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yyyy[k] = 4.0 * to_xyy_yyy[k] - 2.0 * to_xyy_yyyyy[k] * tke_0 - 8.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yyyz[k] = 3.0 * to_xyy_yyz[k] - 2.0 * to_xyy_yyyyz[k] * tke_0 - 6.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yyzz[k] = 2.0 * to_xyy_yzz[k] - 2.0 * to_xyy_yyyzz[k] * tke_0 - 4.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yzzz[k] = to_xyy_zzz[k] - 2.0 * to_xyy_yyzzz[k] * tke_0 - 2.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_zzzz[k] = -2.0 * to_xyy_yzzzz[k] * tke_0 + 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1695-1710 components of targeted buffer : GG

        auto to_z_y_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 120);

        auto to_z_y_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 121);

        auto to_z_y_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 122);

        auto to_z_y_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 123);

        auto to_z_y_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 124);

        auto to_z_y_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 125);

        auto to_z_y_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 126);

        auto to_z_y_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 127);

        auto to_z_y_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 128);

        auto to_z_y_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 129);

        auto to_z_y_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 130);

        auto to_z_y_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 131);

        auto to_z_y_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 132);

        auto to_z_y_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 133);

        auto to_z_y_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyz_xxx, to_xyz_xxxxy, to_xyz_xxxyy, to_xyz_xxxyz, to_xyz_xxy, to_xyz_xxyyy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xyy, to_xyz_xyyyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_yyy, to_xyz_yyyyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_xyzzz_xxx, to_xyzzz_xxxxy, to_xyzzz_xxxyy, to_xyzzz_xxxyz, to_xyzzz_xxy, to_xyzzz_xxyyy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyyyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_yyy, to_xyzzz_yyyyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, to_z_y_xyzz_xxxx, to_z_y_xyzz_xxxy, to_z_y_xyzz_xxxz, to_z_y_xyzz_xxyy, to_z_y_xyzz_xxyz, to_z_y_xyzz_xxzz, to_z_y_xyzz_xyyy, to_z_y_xyzz_xyyz, to_z_y_xyzz_xyzz, to_z_y_xyzz_xzzz, to_z_y_xyzz_yyyy, to_z_y_xyzz_yyyz, to_z_y_xyzz_yyzz, to_z_y_xyzz_yzzz, to_z_y_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyzz_xxxx[k] = -4.0 * to_xyz_xxxxy[k] * tke_0 + 4.0 * to_xyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxxy[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxyy[k] * tke_0 - 2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxxz[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxyy[k] = 4.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyyy[k] * tke_0 - 4.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxyz[k] = 2.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxyyz[k] * tke_0 - 2.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxzz[k] = -4.0 * to_xyz_xxyzz[k] * tke_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xyyy[k] = 6.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyyy[k] * tke_0 - 6.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xyyz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyyyz[k] * tke_0 - 4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xyzz[k] = 2.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xzzz[k] = -4.0 * to_xyz_xyzzz[k] * tke_0 + 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yyyy[k] = 8.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyyy[k] * tke_0 - 8.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yyyz[k] = 6.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyyyz[k] * tke_0 - 6.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yyzz[k] = 4.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 4.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yzzz[k] = 2.0 * to_xyz_zzz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 2.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_zzzz[k] = -4.0 * to_xyz_yzzzz[k] * tke_0 + 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1710-1725 components of targeted buffer : GG

        auto to_z_y_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 135);

        auto to_z_y_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 136);

        auto to_z_y_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 137);

        auto to_z_y_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 138);

        auto to_z_y_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 139);

        auto to_z_y_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 140);

        auto to_z_y_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 141);

        auto to_z_y_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 142);

        auto to_z_y_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 143);

        auto to_z_y_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 144);

        auto to_z_y_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 145);

        auto to_z_y_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 146);

        auto to_z_y_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 147);

        auto to_z_y_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 148);

        auto to_z_y_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xzz_xxx, to_xzz_xxxxy, to_xzz_xxxyy, to_xzz_xxxyz, to_xzz_xxy, to_xzz_xxyyy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xyy, to_xzz_xyyyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_yyy, to_xzz_yyyyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_xzzzz_xxx, to_xzzzz_xxxxy, to_xzzzz_xxxyy, to_xzzzz_xxxyz, to_xzzzz_xxy, to_xzzzz_xxyyy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xyy, to_xzzzz_xyyyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_yyy, to_xzzzz_yyyyy, to_xzzzz_yyyyz, to_xzzzz_yyyzz, to_xzzzz_yyz, to_xzzzz_yyzzz, to_xzzzz_yzz, to_xzzzz_yzzzz, to_xzzzz_zzz, to_z_y_xzzz_xxxx, to_z_y_xzzz_xxxy, to_z_y_xzzz_xxxz, to_z_y_xzzz_xxyy, to_z_y_xzzz_xxyz, to_z_y_xzzz_xxzz, to_z_y_xzzz_xyyy, to_z_y_xzzz_xyyz, to_z_y_xzzz_xyzz, to_z_y_xzzz_xzzz, to_z_y_xzzz_yyyy, to_z_y_xzzz_yyyz, to_z_y_xzzz_yyzz, to_z_y_xzzz_yzzz, to_z_y_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzzz_xxxx[k] = -6.0 * to_xzz_xxxxy[k] * tke_0 + 4.0 * to_xzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxxy[k] = 3.0 * to_xzz_xxx[k] - 6.0 * to_xzz_xxxyy[k] * tke_0 - 2.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxxz[k] = -6.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxyy[k] = 6.0 * to_xzz_xxy[k] - 6.0 * to_xzz_xxyyy[k] * tke_0 - 4.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxyz[k] = 3.0 * to_xzz_xxz[k] - 6.0 * to_xzz_xxyyz[k] * tke_0 - 2.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxzz[k] = -6.0 * to_xzz_xxyzz[k] * tke_0 + 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xyyy[k] = 9.0 * to_xzz_xyy[k] - 6.0 * to_xzz_xyyyy[k] * tke_0 - 6.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xyyz[k] = 6.0 * to_xzz_xyz[k] - 6.0 * to_xzz_xyyyz[k] * tke_0 - 4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xyzz[k] = 3.0 * to_xzz_xzz[k] - 6.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xzzz[k] = -6.0 * to_xzz_xyzzz[k] * tke_0 + 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yyyy[k] = 12.0 * to_xzz_yyy[k] - 6.0 * to_xzz_yyyyy[k] * tke_0 - 8.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yyyz[k] = 9.0 * to_xzz_yyz[k] - 6.0 * to_xzz_yyyyz[k] * tke_0 - 6.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yyzz[k] = 6.0 * to_xzz_yzz[k] - 6.0 * to_xzz_yyyzz[k] * tke_0 - 4.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yzzz[k] = 3.0 * to_xzz_zzz[k] - 6.0 * to_xzz_yyzzz[k] * tke_0 - 2.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_zzzz[k] = -6.0 * to_xzz_yzzzz[k] * tke_0 + 4.0 * to_xzzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1725-1740 components of targeted buffer : GG

        auto to_z_y_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 150);

        auto to_z_y_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 151);

        auto to_z_y_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 152);

        auto to_z_y_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 153);

        auto to_z_y_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 154);

        auto to_z_y_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 155);

        auto to_z_y_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 156);

        auto to_z_y_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 157);

        auto to_z_y_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 158);

        auto to_z_y_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 159);

        auto to_z_y_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 160);

        auto to_z_y_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 161);

        auto to_z_y_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 162);

        auto to_z_y_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 163);

        auto to_z_y_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_yyyyz_xxx, to_yyyyz_xxxxy, to_yyyyz_xxxyy, to_yyyyz_xxxyz, to_yyyyz_xxy, to_yyyyz_xxyyy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xyy, to_yyyyz_xyyyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_yyy, to_yyyyz_yyyyy, to_yyyyz_yyyyz, to_yyyyz_yyyzz, to_yyyyz_yyz, to_yyyyz_yyzzz, to_yyyyz_yzz, to_yyyyz_yzzzz, to_yyyyz_zzz, to_z_y_yyyy_xxxx, to_z_y_yyyy_xxxy, to_z_y_yyyy_xxxz, to_z_y_yyyy_xxyy, to_z_y_yyyy_xxyz, to_z_y_yyyy_xxzz, to_z_y_yyyy_xyyy, to_z_y_yyyy_xyyz, to_z_y_yyyy_xyzz, to_z_y_yyyy_xzzz, to_z_y_yyyy_yyyy, to_z_y_yyyy_yyyz, to_z_y_yyyy_yyzz, to_z_y_yyyy_yzzz, to_z_y_yyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyy_xxxx[k] = 4.0 * to_yyyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxxy[k] = -2.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxxz[k] = 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxyy[k] = -4.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxyz[k] = -2.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxzz[k] = 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xyyy[k] = -6.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xyyz[k] = -4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xyzz[k] = -2.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xzzz[k] = 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yyyy[k] = -8.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yyyz[k] = -6.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yyzz[k] = -4.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yzzz[k] = -2.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_zzzz[k] = 4.0 * to_yyyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1740-1755 components of targeted buffer : GG

        auto to_z_y_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 165);

        auto to_z_y_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 166);

        auto to_z_y_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 167);

        auto to_z_y_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 168);

        auto to_z_y_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 169);

        auto to_z_y_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 170);

        auto to_z_y_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 171);

        auto to_z_y_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 172);

        auto to_z_y_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 173);

        auto to_z_y_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 174);

        auto to_z_y_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 175);

        auto to_z_y_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 176);

        auto to_z_y_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 177);

        auto to_z_y_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 178);

        auto to_z_y_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_yyy_xxx, to_yyy_xxxxy, to_yyy_xxxyy, to_yyy_xxxyz, to_yyy_xxy, to_yyy_xxyyy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xyy, to_yyy_xyyyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_yyy, to_yyy_yyyyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyyzz_xxx, to_yyyzz_xxxxy, to_yyyzz_xxxyy, to_yyyzz_xxxyz, to_yyyzz_xxy, to_yyyzz_xxyyy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xyy, to_yyyzz_xyyyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_yyy, to_yyyzz_yyyyy, to_yyyzz_yyyyz, to_yyyzz_yyyzz, to_yyyzz_yyz, to_yyyzz_yyzzz, to_yyyzz_yzz, to_yyyzz_yzzzz, to_yyyzz_zzz, to_z_y_yyyz_xxxx, to_z_y_yyyz_xxxy, to_z_y_yyyz_xxxz, to_z_y_yyyz_xxyy, to_z_y_yyyz_xxyz, to_z_y_yyyz_xxzz, to_z_y_yyyz_xyyy, to_z_y_yyyz_xyyz, to_z_y_yyyz_xyzz, to_z_y_yyyz_xzzz, to_z_y_yyyz_yyyy, to_z_y_yyyz_yyyz, to_z_y_yyyz_yyzz, to_z_y_yyyz_yzzz, to_z_y_yyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyz_xxxx[k] = -2.0 * to_yyy_xxxxy[k] * tke_0 + 4.0 * to_yyyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxxy[k] = to_yyy_xxx[k] - 2.0 * to_yyy_xxxyy[k] * tke_0 - 2.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxxz[k] = -2.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxyy[k] = 2.0 * to_yyy_xxy[k] - 2.0 * to_yyy_xxyyy[k] * tke_0 - 4.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxyz[k] = to_yyy_xxz[k] - 2.0 * to_yyy_xxyyz[k] * tke_0 - 2.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxzz[k] = -2.0 * to_yyy_xxyzz[k] * tke_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xyyy[k] = 3.0 * to_yyy_xyy[k] - 2.0 * to_yyy_xyyyy[k] * tke_0 - 6.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xyyz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xyyyz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xyzz[k] = to_yyy_xzz[k] - 2.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xzzz[k] = -2.0 * to_yyy_xyzzz[k] * tke_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yyyy[k] = 4.0 * to_yyy_yyy[k] - 2.0 * to_yyy_yyyyy[k] * tke_0 - 8.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yyyz[k] = 3.0 * to_yyy_yyz[k] - 2.0 * to_yyy_yyyyz[k] * tke_0 - 6.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yyzz[k] = 2.0 * to_yyy_yzz[k] - 2.0 * to_yyy_yyyzz[k] * tke_0 - 4.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yzzz[k] = to_yyy_zzz[k] - 2.0 * to_yyy_yyzzz[k] * tke_0 - 2.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_zzzz[k] = -2.0 * to_yyy_yzzzz[k] * tke_0 + 4.0 * to_yyyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1755-1770 components of targeted buffer : GG

        auto to_z_y_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 180);

        auto to_z_y_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 181);

        auto to_z_y_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 182);

        auto to_z_y_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 183);

        auto to_z_y_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 184);

        auto to_z_y_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 185);

        auto to_z_y_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 186);

        auto to_z_y_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 187);

        auto to_z_y_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 188);

        auto to_z_y_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 189);

        auto to_z_y_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 190);

        auto to_z_y_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 191);

        auto to_z_y_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 192);

        auto to_z_y_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 193);

        auto to_z_y_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_yyz_xxx, to_yyz_xxxxy, to_yyz_xxxyy, to_yyz_xxxyz, to_yyz_xxy, to_yyz_xxyyy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xyy, to_yyz_xyyyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_yyy, to_yyz_yyyyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, to_yyzzz_xxx, to_yyzzz_xxxxy, to_yyzzz_xxxyy, to_yyzzz_xxxyz, to_yyzzz_xxy, to_yyzzz_xxyyy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xyy, to_yyzzz_xyyyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_yyy, to_yyzzz_yyyyy, to_yyzzz_yyyyz, to_yyzzz_yyyzz, to_yyzzz_yyz, to_yyzzz_yyzzz, to_yyzzz_yzz, to_yyzzz_yzzzz, to_yyzzz_zzz, to_z_y_yyzz_xxxx, to_z_y_yyzz_xxxy, to_z_y_yyzz_xxxz, to_z_y_yyzz_xxyy, to_z_y_yyzz_xxyz, to_z_y_yyzz_xxzz, to_z_y_yyzz_xyyy, to_z_y_yyzz_xyyz, to_z_y_yyzz_xyzz, to_z_y_yyzz_xzzz, to_z_y_yyzz_yyyy, to_z_y_yyzz_yyyz, to_z_y_yyzz_yyzz, to_z_y_yyzz_yzzz, to_z_y_yyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyzz_xxxx[k] = -4.0 * to_yyz_xxxxy[k] * tke_0 + 4.0 * to_yyzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxxy[k] = 2.0 * to_yyz_xxx[k] - 4.0 * to_yyz_xxxyy[k] * tke_0 - 2.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxxz[k] = -4.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxyy[k] = 4.0 * to_yyz_xxy[k] - 4.0 * to_yyz_xxyyy[k] * tke_0 - 4.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxyz[k] = 2.0 * to_yyz_xxz[k] - 4.0 * to_yyz_xxyyz[k] * tke_0 - 2.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxzz[k] = -4.0 * to_yyz_xxyzz[k] * tke_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xyyy[k] = 6.0 * to_yyz_xyy[k] - 4.0 * to_yyz_xyyyy[k] * tke_0 - 6.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xyyz[k] = 4.0 * to_yyz_xyz[k] - 4.0 * to_yyz_xyyyz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xyzz[k] = 2.0 * to_yyz_xzz[k] - 4.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xzzz[k] = -4.0 * to_yyz_xyzzz[k] * tke_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yyyy[k] = 8.0 * to_yyz_yyy[k] - 4.0 * to_yyz_yyyyy[k] * tke_0 - 8.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yyyz[k] = 6.0 * to_yyz_yyz[k] - 4.0 * to_yyz_yyyyz[k] * tke_0 - 6.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yyzz[k] = 4.0 * to_yyz_yzz[k] - 4.0 * to_yyz_yyyzz[k] * tke_0 - 4.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yzzz[k] = 2.0 * to_yyz_zzz[k] - 4.0 * to_yyz_yyzzz[k] * tke_0 - 2.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_zzzz[k] = -4.0 * to_yyz_yzzzz[k] * tke_0 + 4.0 * to_yyzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1770-1785 components of targeted buffer : GG

        auto to_z_y_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 195);

        auto to_z_y_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 196);

        auto to_z_y_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 197);

        auto to_z_y_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 198);

        auto to_z_y_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 199);

        auto to_z_y_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 200);

        auto to_z_y_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 201);

        auto to_z_y_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 202);

        auto to_z_y_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 203);

        auto to_z_y_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 204);

        auto to_z_y_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 205);

        auto to_z_y_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 206);

        auto to_z_y_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 207);

        auto to_z_y_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 208);

        auto to_z_y_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_yzz_xxx, to_yzz_xxxxy, to_yzz_xxxyy, to_yzz_xxxyz, to_yzz_xxy, to_yzz_xxyyy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xyy, to_yzz_xyyyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_yyy, to_yzz_yyyyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, to_yzzzz_xxx, to_yzzzz_xxxxy, to_yzzzz_xxxyy, to_yzzzz_xxxyz, to_yzzzz_xxy, to_yzzzz_xxyyy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xyy, to_yzzzz_xyyyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_yyy, to_yzzzz_yyyyy, to_yzzzz_yyyyz, to_yzzzz_yyyzz, to_yzzzz_yyz, to_yzzzz_yyzzz, to_yzzzz_yzz, to_yzzzz_yzzzz, to_yzzzz_zzz, to_z_y_yzzz_xxxx, to_z_y_yzzz_xxxy, to_z_y_yzzz_xxxz, to_z_y_yzzz_xxyy, to_z_y_yzzz_xxyz, to_z_y_yzzz_xxzz, to_z_y_yzzz_xyyy, to_z_y_yzzz_xyyz, to_z_y_yzzz_xyzz, to_z_y_yzzz_xzzz, to_z_y_yzzz_yyyy, to_z_y_yzzz_yyyz, to_z_y_yzzz_yyzz, to_z_y_yzzz_yzzz, to_z_y_yzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzzz_xxxx[k] = -6.0 * to_yzz_xxxxy[k] * tke_0 + 4.0 * to_yzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxxy[k] = 3.0 * to_yzz_xxx[k] - 6.0 * to_yzz_xxxyy[k] * tke_0 - 2.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxxz[k] = -6.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxyy[k] = 6.0 * to_yzz_xxy[k] - 6.0 * to_yzz_xxyyy[k] * tke_0 - 4.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxyz[k] = 3.0 * to_yzz_xxz[k] - 6.0 * to_yzz_xxyyz[k] * tke_0 - 2.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxzz[k] = -6.0 * to_yzz_xxyzz[k] * tke_0 + 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xyyy[k] = 9.0 * to_yzz_xyy[k] - 6.0 * to_yzz_xyyyy[k] * tke_0 - 6.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xyyz[k] = 6.0 * to_yzz_xyz[k] - 6.0 * to_yzz_xyyyz[k] * tke_0 - 4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xyzz[k] = 3.0 * to_yzz_xzz[k] - 6.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xzzz[k] = -6.0 * to_yzz_xyzzz[k] * tke_0 + 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yyyy[k] = 12.0 * to_yzz_yyy[k] - 6.0 * to_yzz_yyyyy[k] * tke_0 - 8.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yyyz[k] = 9.0 * to_yzz_yyz[k] - 6.0 * to_yzz_yyyyz[k] * tke_0 - 6.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yyzz[k] = 6.0 * to_yzz_yzz[k] - 6.0 * to_yzz_yyyzz[k] * tke_0 - 4.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yzzz[k] = 3.0 * to_yzz_zzz[k] - 6.0 * to_yzz_yyzzz[k] * tke_0 - 2.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_zzzz[k] = -6.0 * to_yzz_yzzzz[k] * tke_0 + 4.0 * to_yzzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1785-1800 components of targeted buffer : GG

        auto to_z_y_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 210);

        auto to_z_y_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 211);

        auto to_z_y_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 212);

        auto to_z_y_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 213);

        auto to_z_y_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 214);

        auto to_z_y_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 215);

        auto to_z_y_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 216);

        auto to_z_y_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 217);

        auto to_z_y_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 218);

        auto to_z_y_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 219);

        auto to_z_y_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 220);

        auto to_z_y_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 221);

        auto to_z_y_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 222);

        auto to_z_y_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 223);

        auto to_z_y_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 7 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_z_y_zzzz_xxxx, to_z_y_zzzz_xxxy, to_z_y_zzzz_xxxz, to_z_y_zzzz_xxyy, to_z_y_zzzz_xxyz, to_z_y_zzzz_xxzz, to_z_y_zzzz_xyyy, to_z_y_zzzz_xyyz, to_z_y_zzzz_xyzz, to_z_y_zzzz_xzzz, to_z_y_zzzz_yyyy, to_z_y_zzzz_yyyz, to_z_y_zzzz_yyzz, to_z_y_zzzz_yzzz, to_z_y_zzzz_zzzz, to_zzz_xxx, to_zzz_xxxxy, to_zzz_xxxyy, to_zzz_xxxyz, to_zzz_xxy, to_zzz_xxyyy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xyy, to_zzz_xyyyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_yyy, to_zzz_yyyyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, to_zzzzz_xxx, to_zzzzz_xxxxy, to_zzzzz_xxxyy, to_zzzzz_xxxyz, to_zzzzz_xxy, to_zzzzz_xxyyy, to_zzzzz_xxyyz, to_zzzzz_xxyzz, to_zzzzz_xxz, to_zzzzz_xyy, to_zzzzz_xyyyy, to_zzzzz_xyyyz, to_zzzzz_xyyzz, to_zzzzz_xyz, to_zzzzz_xyzzz, to_zzzzz_xzz, to_zzzzz_yyy, to_zzzzz_yyyyy, to_zzzzz_yyyyz, to_zzzzz_yyyzz, to_zzzzz_yyz, to_zzzzz_yyzzz, to_zzzzz_yzz, to_zzzzz_yzzzz, to_zzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzzz_xxxx[k] = -8.0 * to_zzz_xxxxy[k] * tke_0 + 4.0 * to_zzzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxxy[k] = 4.0 * to_zzz_xxx[k] - 8.0 * to_zzz_xxxyy[k] * tke_0 - 2.0 * to_zzzzz_xxx[k] * tbe_0 + 4.0 * to_zzzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxxz[k] = -8.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_zzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxyy[k] = 8.0 * to_zzz_xxy[k] - 8.0 * to_zzz_xxyyy[k] * tke_0 - 4.0 * to_zzzzz_xxy[k] * tbe_0 + 4.0 * to_zzzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxyz[k] = 4.0 * to_zzz_xxz[k] - 8.0 * to_zzz_xxyyz[k] * tke_0 - 2.0 * to_zzzzz_xxz[k] * tbe_0 + 4.0 * to_zzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxzz[k] = -8.0 * to_zzz_xxyzz[k] * tke_0 + 4.0 * to_zzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xyyy[k] = 12.0 * to_zzz_xyy[k] - 8.0 * to_zzz_xyyyy[k] * tke_0 - 6.0 * to_zzzzz_xyy[k] * tbe_0 + 4.0 * to_zzzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xyyz[k] = 8.0 * to_zzz_xyz[k] - 8.0 * to_zzz_xyyyz[k] * tke_0 - 4.0 * to_zzzzz_xyz[k] * tbe_0 + 4.0 * to_zzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xyzz[k] = 4.0 * to_zzz_xzz[k] - 8.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_zzzzz_xzz[k] * tbe_0 + 4.0 * to_zzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xzzz[k] = -8.0 * to_zzz_xyzzz[k] * tke_0 + 4.0 * to_zzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yyyy[k] = 16.0 * to_zzz_yyy[k] - 8.0 * to_zzz_yyyyy[k] * tke_0 - 8.0 * to_zzzzz_yyy[k] * tbe_0 + 4.0 * to_zzzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yyyz[k] = 12.0 * to_zzz_yyz[k] - 8.0 * to_zzz_yyyyz[k] * tke_0 - 6.0 * to_zzzzz_yyz[k] * tbe_0 + 4.0 * to_zzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yyzz[k] = 8.0 * to_zzz_yzz[k] - 8.0 * to_zzz_yyyzz[k] * tke_0 - 4.0 * to_zzzzz_yzz[k] * tbe_0 + 4.0 * to_zzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yzzz[k] = 4.0 * to_zzz_zzz[k] - 8.0 * to_zzz_yyzzz[k] * tke_0 - 2.0 * to_zzzzz_zzz[k] * tbe_0 + 4.0 * to_zzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_zzzz[k] = -8.0 * to_zzz_yzzzz[k] * tke_0 + 4.0 * to_zzzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1800-1815 components of targeted buffer : GG

        auto to_z_z_xxxx_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 0);

        auto to_z_z_xxxx_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 1);

        auto to_z_z_xxxx_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 2);

        auto to_z_z_xxxx_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 3);

        auto to_z_z_xxxx_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 4);

        auto to_z_z_xxxx_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 5);

        auto to_z_z_xxxx_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 6);

        auto to_z_z_xxxx_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 7);

        auto to_z_z_xxxx_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 8);

        auto to_z_z_xxxx_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 9);

        auto to_z_z_xxxx_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 10);

        auto to_z_z_xxxx_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 11);

        auto to_z_z_xxxx_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 12);

        auto to_z_z_xxxx_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 13);

        auto to_z_z_xxxx_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_xxxxz_xxx, to_xxxxz_xxxxz, to_xxxxz_xxxyz, to_xxxxz_xxxzz, to_xxxxz_xxy, to_xxxxz_xxyyz, to_xxxxz_xxyzz, to_xxxxz_xxz, to_xxxxz_xxzzz, to_xxxxz_xyy, to_xxxxz_xyyyz, to_xxxxz_xyyzz, to_xxxxz_xyz, to_xxxxz_xyzzz, to_xxxxz_xzz, to_xxxxz_xzzzz, to_xxxxz_yyy, to_xxxxz_yyyyz, to_xxxxz_yyyzz, to_xxxxz_yyz, to_xxxxz_yyzzz, to_xxxxz_yzz, to_xxxxz_yzzzz, to_xxxxz_zzz, to_xxxxz_zzzzz, to_z_z_xxxx_xxxx, to_z_z_xxxx_xxxy, to_z_z_xxxx_xxxz, to_z_z_xxxx_xxyy, to_z_z_xxxx_xxyz, to_z_z_xxxx_xxzz, to_z_z_xxxx_xyyy, to_z_z_xxxx_xyyz, to_z_z_xxxx_xyzz, to_z_z_xxxx_xzzz, to_z_z_xxxx_yyyy, to_z_z_xxxx_yyyz, to_z_z_xxxx_yyzz, to_z_z_xxxx_yzzz, to_z_z_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxx_xxxx[k] = 4.0 * to_xxxxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxxy[k] = 4.0 * to_xxxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxxz[k] = -2.0 * to_xxxxz_xxx[k] * tbe_0 + 4.0 * to_xxxxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxyy[k] = 4.0 * to_xxxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxyz[k] = -2.0 * to_xxxxz_xxy[k] * tbe_0 + 4.0 * to_xxxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxzz[k] = -4.0 * to_xxxxz_xxz[k] * tbe_0 + 4.0 * to_xxxxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xyyy[k] = 4.0 * to_xxxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xyyz[k] = -2.0 * to_xxxxz_xyy[k] * tbe_0 + 4.0 * to_xxxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xyzz[k] = -4.0 * to_xxxxz_xyz[k] * tbe_0 + 4.0 * to_xxxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xzzz[k] = -6.0 * to_xxxxz_xzz[k] * tbe_0 + 4.0 * to_xxxxz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yyyy[k] = 4.0 * to_xxxxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yyyz[k] = -2.0 * to_xxxxz_yyy[k] * tbe_0 + 4.0 * to_xxxxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yyzz[k] = -4.0 * to_xxxxz_yyz[k] * tbe_0 + 4.0 * to_xxxxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yzzz[k] = -6.0 * to_xxxxz_yzz[k] * tbe_0 + 4.0 * to_xxxxz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_zzzz[k] = -8.0 * to_xxxxz_zzz[k] * tbe_0 + 4.0 * to_xxxxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1815-1830 components of targeted buffer : GG

        auto to_z_z_xxxy_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 15);

        auto to_z_z_xxxy_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 16);

        auto to_z_z_xxxy_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 17);

        auto to_z_z_xxxy_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 18);

        auto to_z_z_xxxy_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 19);

        auto to_z_z_xxxy_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 20);

        auto to_z_z_xxxy_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 21);

        auto to_z_z_xxxy_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 22);

        auto to_z_z_xxxy_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 23);

        auto to_z_z_xxxy_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 24);

        auto to_z_z_xxxy_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 25);

        auto to_z_z_xxxy_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 26);

        auto to_z_z_xxxy_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 27);

        auto to_z_z_xxxy_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 28);

        auto to_z_z_xxxy_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_xxxyz_xxx, to_xxxyz_xxxxz, to_xxxyz_xxxyz, to_xxxyz_xxxzz, to_xxxyz_xxy, to_xxxyz_xxyyz, to_xxxyz_xxyzz, to_xxxyz_xxz, to_xxxyz_xxzzz, to_xxxyz_xyy, to_xxxyz_xyyyz, to_xxxyz_xyyzz, to_xxxyz_xyz, to_xxxyz_xyzzz, to_xxxyz_xzz, to_xxxyz_xzzzz, to_xxxyz_yyy, to_xxxyz_yyyyz, to_xxxyz_yyyzz, to_xxxyz_yyz, to_xxxyz_yyzzz, to_xxxyz_yzz, to_xxxyz_yzzzz, to_xxxyz_zzz, to_xxxyz_zzzzz, to_z_z_xxxy_xxxx, to_z_z_xxxy_xxxy, to_z_z_xxxy_xxxz, to_z_z_xxxy_xxyy, to_z_z_xxxy_xxyz, to_z_z_xxxy_xxzz, to_z_z_xxxy_xyyy, to_z_z_xxxy_xyyz, to_z_z_xxxy_xyzz, to_z_z_xxxy_xzzz, to_z_z_xxxy_yyyy, to_z_z_xxxy_yyyz, to_z_z_xxxy_yyzz, to_z_z_xxxy_yzzz, to_z_z_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxy_xxxx[k] = 4.0 * to_xxxyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxxy[k] = 4.0 * to_xxxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxxz[k] = -2.0 * to_xxxyz_xxx[k] * tbe_0 + 4.0 * to_xxxyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxyy[k] = 4.0 * to_xxxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxyz[k] = -2.0 * to_xxxyz_xxy[k] * tbe_0 + 4.0 * to_xxxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxzz[k] = -4.0 * to_xxxyz_xxz[k] * tbe_0 + 4.0 * to_xxxyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xyyy[k] = 4.0 * to_xxxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xyyz[k] = -2.0 * to_xxxyz_xyy[k] * tbe_0 + 4.0 * to_xxxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xyzz[k] = -4.0 * to_xxxyz_xyz[k] * tbe_0 + 4.0 * to_xxxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xzzz[k] = -6.0 * to_xxxyz_xzz[k] * tbe_0 + 4.0 * to_xxxyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yyyy[k] = 4.0 * to_xxxyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yyyz[k] = -2.0 * to_xxxyz_yyy[k] * tbe_0 + 4.0 * to_xxxyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yyzz[k] = -4.0 * to_xxxyz_yyz[k] * tbe_0 + 4.0 * to_xxxyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yzzz[k] = -6.0 * to_xxxyz_yzz[k] * tbe_0 + 4.0 * to_xxxyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_zzzz[k] = -8.0 * to_xxxyz_zzz[k] * tbe_0 + 4.0 * to_xxxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1830-1845 components of targeted buffer : GG

        auto to_z_z_xxxz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 30);

        auto to_z_z_xxxz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 31);

        auto to_z_z_xxxz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 32);

        auto to_z_z_xxxz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 33);

        auto to_z_z_xxxz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 34);

        auto to_z_z_xxxz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 35);

        auto to_z_z_xxxz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 36);

        auto to_z_z_xxxz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 37);

        auto to_z_z_xxxz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 38);

        auto to_z_z_xxxz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 39);

        auto to_z_z_xxxz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 40);

        auto to_z_z_xxxz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 41);

        auto to_z_z_xxxz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 42);

        auto to_z_z_xxxz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 43);

        auto to_z_z_xxxz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_xxx_xxx, to_xxx_xxxxz, to_xxx_xxxyz, to_xxx_xxxzz, to_xxx_xxy, to_xxx_xxyyz, to_xxx_xxyzz, to_xxx_xxz, to_xxx_xxzzz, to_xxx_xyy, to_xxx_xyyyz, to_xxx_xyyzz, to_xxx_xyz, to_xxx_xyzzz, to_xxx_xzz, to_xxx_xzzzz, to_xxx_yyy, to_xxx_yyyyz, to_xxx_yyyzz, to_xxx_yyz, to_xxx_yyzzz, to_xxx_yzz, to_xxx_yzzzz, to_xxx_zzz, to_xxx_zzzzz, to_xxxzz_xxx, to_xxxzz_xxxxz, to_xxxzz_xxxyz, to_xxxzz_xxxzz, to_xxxzz_xxy, to_xxxzz_xxyyz, to_xxxzz_xxyzz, to_xxxzz_xxz, to_xxxzz_xxzzz, to_xxxzz_xyy, to_xxxzz_xyyyz, to_xxxzz_xyyzz, to_xxxzz_xyz, to_xxxzz_xyzzz, to_xxxzz_xzz, to_xxxzz_xzzzz, to_xxxzz_yyy, to_xxxzz_yyyyz, to_xxxzz_yyyzz, to_xxxzz_yyz, to_xxxzz_yyzzz, to_xxxzz_yzz, to_xxxzz_yzzzz, to_xxxzz_zzz, to_xxxzz_zzzzz, to_z_z_xxxz_xxxx, to_z_z_xxxz_xxxy, to_z_z_xxxz_xxxz, to_z_z_xxxz_xxyy, to_z_z_xxxz_xxyz, to_z_z_xxxz_xxzz, to_z_z_xxxz_xyyy, to_z_z_xxxz_xyyz, to_z_z_xxxz_xyzz, to_z_z_xxxz_xzzz, to_z_z_xxxz_yyyy, to_z_z_xxxz_yyyz, to_z_z_xxxz_yyzz, to_z_z_xxxz_yzzz, to_z_z_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxz_xxxx[k] = -2.0 * to_xxx_xxxxz[k] * tke_0 + 4.0 * to_xxxzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxxy[k] = -2.0 * to_xxx_xxxyz[k] * tke_0 + 4.0 * to_xxxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxxz[k] = to_xxx_xxx[k] - 2.0 * to_xxx_xxxzz[k] * tke_0 - 2.0 * to_xxxzz_xxx[k] * tbe_0 + 4.0 * to_xxxzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxyy[k] = -2.0 * to_xxx_xxyyz[k] * tke_0 + 4.0 * to_xxxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxyz[k] = to_xxx_xxy[k] - 2.0 * to_xxx_xxyzz[k] * tke_0 - 2.0 * to_xxxzz_xxy[k] * tbe_0 + 4.0 * to_xxxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxzz[k] = 2.0 * to_xxx_xxz[k] - 2.0 * to_xxx_xxzzz[k] * tke_0 - 4.0 * to_xxxzz_xxz[k] * tbe_0 + 4.0 * to_xxxzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xyyy[k] = -2.0 * to_xxx_xyyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xyyz[k] = to_xxx_xyy[k] - 2.0 * to_xxx_xyyzz[k] * tke_0 - 2.0 * to_xxxzz_xyy[k] * tbe_0 + 4.0 * to_xxxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xyzz[k] = 2.0 * to_xxx_xyz[k] - 2.0 * to_xxx_xyzzz[k] * tke_0 - 4.0 * to_xxxzz_xyz[k] * tbe_0 + 4.0 * to_xxxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xzzz[k] = 3.0 * to_xxx_xzz[k] - 2.0 * to_xxx_xzzzz[k] * tke_0 - 6.0 * to_xxxzz_xzz[k] * tbe_0 + 4.0 * to_xxxzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yyyy[k] = -2.0 * to_xxx_yyyyz[k] * tke_0 + 4.0 * to_xxxzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yyyz[k] = to_xxx_yyy[k] - 2.0 * to_xxx_yyyzz[k] * tke_0 - 2.0 * to_xxxzz_yyy[k] * tbe_0 + 4.0 * to_xxxzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yyzz[k] = 2.0 * to_xxx_yyz[k] - 2.0 * to_xxx_yyzzz[k] * tke_0 - 4.0 * to_xxxzz_yyz[k] * tbe_0 + 4.0 * to_xxxzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yzzz[k] = 3.0 * to_xxx_yzz[k] - 2.0 * to_xxx_yzzzz[k] * tke_0 - 6.0 * to_xxxzz_yzz[k] * tbe_0 + 4.0 * to_xxxzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_zzzz[k] = 4.0 * to_xxx_zzz[k] - 2.0 * to_xxx_zzzzz[k] * tke_0 - 8.0 * to_xxxzz_zzz[k] * tbe_0 + 4.0 * to_xxxzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1845-1860 components of targeted buffer : GG

        auto to_z_z_xxyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 45);

        auto to_z_z_xxyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 46);

        auto to_z_z_xxyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 47);

        auto to_z_z_xxyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 48);

        auto to_z_z_xxyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 49);

        auto to_z_z_xxyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 50);

        auto to_z_z_xxyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 51);

        auto to_z_z_xxyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 52);

        auto to_z_z_xxyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 53);

        auto to_z_z_xxyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 54);

        auto to_z_z_xxyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 55);

        auto to_z_z_xxyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 56);

        auto to_z_z_xxyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 57);

        auto to_z_z_xxyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 58);

        auto to_z_z_xxyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_xxyyz_xxx, to_xxyyz_xxxxz, to_xxyyz_xxxyz, to_xxyyz_xxxzz, to_xxyyz_xxy, to_xxyyz_xxyyz, to_xxyyz_xxyzz, to_xxyyz_xxz, to_xxyyz_xxzzz, to_xxyyz_xyy, to_xxyyz_xyyyz, to_xxyyz_xyyzz, to_xxyyz_xyz, to_xxyyz_xyzzz, to_xxyyz_xzz, to_xxyyz_xzzzz, to_xxyyz_yyy, to_xxyyz_yyyyz, to_xxyyz_yyyzz, to_xxyyz_yyz, to_xxyyz_yyzzz, to_xxyyz_yzz, to_xxyyz_yzzzz, to_xxyyz_zzz, to_xxyyz_zzzzz, to_z_z_xxyy_xxxx, to_z_z_xxyy_xxxy, to_z_z_xxyy_xxxz, to_z_z_xxyy_xxyy, to_z_z_xxyy_xxyz, to_z_z_xxyy_xxzz, to_z_z_xxyy_xyyy, to_z_z_xxyy_xyyz, to_z_z_xxyy_xyzz, to_z_z_xxyy_xzzz, to_z_z_xxyy_yyyy, to_z_z_xxyy_yyyz, to_z_z_xxyy_yyzz, to_z_z_xxyy_yzzz, to_z_z_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyy_xxxx[k] = 4.0 * to_xxyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxxy[k] = 4.0 * to_xxyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxxz[k] = -2.0 * to_xxyyz_xxx[k] * tbe_0 + 4.0 * to_xxyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxyy[k] = 4.0 * to_xxyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxyz[k] = -2.0 * to_xxyyz_xxy[k] * tbe_0 + 4.0 * to_xxyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxzz[k] = -4.0 * to_xxyyz_xxz[k] * tbe_0 + 4.0 * to_xxyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xyyy[k] = 4.0 * to_xxyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xyyz[k] = -2.0 * to_xxyyz_xyy[k] * tbe_0 + 4.0 * to_xxyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xyzz[k] = -4.0 * to_xxyyz_xyz[k] * tbe_0 + 4.0 * to_xxyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xzzz[k] = -6.0 * to_xxyyz_xzz[k] * tbe_0 + 4.0 * to_xxyyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yyyy[k] = 4.0 * to_xxyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yyyz[k] = -2.0 * to_xxyyz_yyy[k] * tbe_0 + 4.0 * to_xxyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yyzz[k] = -4.0 * to_xxyyz_yyz[k] * tbe_0 + 4.0 * to_xxyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yzzz[k] = -6.0 * to_xxyyz_yzz[k] * tbe_0 + 4.0 * to_xxyyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_zzzz[k] = -8.0 * to_xxyyz_zzz[k] * tbe_0 + 4.0 * to_xxyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1860-1875 components of targeted buffer : GG

        auto to_z_z_xxyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 60);

        auto to_z_z_xxyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 61);

        auto to_z_z_xxyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 62);

        auto to_z_z_xxyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 63);

        auto to_z_z_xxyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 64);

        auto to_z_z_xxyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 65);

        auto to_z_z_xxyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 66);

        auto to_z_z_xxyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 67);

        auto to_z_z_xxyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 68);

        auto to_z_z_xxyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 69);

        auto to_z_z_xxyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 70);

        auto to_z_z_xxyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 71);

        auto to_z_z_xxyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 72);

        auto to_z_z_xxyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 73);

        auto to_z_z_xxyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_xxy_xxx, to_xxy_xxxxz, to_xxy_xxxyz, to_xxy_xxxzz, to_xxy_xxy, to_xxy_xxyyz, to_xxy_xxyzz, to_xxy_xxz, to_xxy_xxzzz, to_xxy_xyy, to_xxy_xyyyz, to_xxy_xyyzz, to_xxy_xyz, to_xxy_xyzzz, to_xxy_xzz, to_xxy_xzzzz, to_xxy_yyy, to_xxy_yyyyz, to_xxy_yyyzz, to_xxy_yyz, to_xxy_yyzzz, to_xxy_yzz, to_xxy_yzzzz, to_xxy_zzz, to_xxy_zzzzz, to_xxyzz_xxx, to_xxyzz_xxxxz, to_xxyzz_xxxyz, to_xxyzz_xxxzz, to_xxyzz_xxy, to_xxyzz_xxyyz, to_xxyzz_xxyzz, to_xxyzz_xxz, to_xxyzz_xxzzz, to_xxyzz_xyy, to_xxyzz_xyyyz, to_xxyzz_xyyzz, to_xxyzz_xyz, to_xxyzz_xyzzz, to_xxyzz_xzz, to_xxyzz_xzzzz, to_xxyzz_yyy, to_xxyzz_yyyyz, to_xxyzz_yyyzz, to_xxyzz_yyz, to_xxyzz_yyzzz, to_xxyzz_yzz, to_xxyzz_yzzzz, to_xxyzz_zzz, to_xxyzz_zzzzz, to_z_z_xxyz_xxxx, to_z_z_xxyz_xxxy, to_z_z_xxyz_xxxz, to_z_z_xxyz_xxyy, to_z_z_xxyz_xxyz, to_z_z_xxyz_xxzz, to_z_z_xxyz_xyyy, to_z_z_xxyz_xyyz, to_z_z_xxyz_xyzz, to_z_z_xxyz_xzzz, to_z_z_xxyz_yyyy, to_z_z_xxyz_yyyz, to_z_z_xxyz_yyzz, to_z_z_xxyz_yzzz, to_z_z_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyz_xxxx[k] = -2.0 * to_xxy_xxxxz[k] * tke_0 + 4.0 * to_xxyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxxy[k] = -2.0 * to_xxy_xxxyz[k] * tke_0 + 4.0 * to_xxyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxxz[k] = to_xxy_xxx[k] - 2.0 * to_xxy_xxxzz[k] * tke_0 - 2.0 * to_xxyzz_xxx[k] * tbe_0 + 4.0 * to_xxyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxyy[k] = -2.0 * to_xxy_xxyyz[k] * tke_0 + 4.0 * to_xxyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxyz[k] = to_xxy_xxy[k] - 2.0 * to_xxy_xxyzz[k] * tke_0 - 2.0 * to_xxyzz_xxy[k] * tbe_0 + 4.0 * to_xxyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxzz[k] = 2.0 * to_xxy_xxz[k] - 2.0 * to_xxy_xxzzz[k] * tke_0 - 4.0 * to_xxyzz_xxz[k] * tbe_0 + 4.0 * to_xxyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xyyy[k] = -2.0 * to_xxy_xyyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xyyz[k] = to_xxy_xyy[k] - 2.0 * to_xxy_xyyzz[k] * tke_0 - 2.0 * to_xxyzz_xyy[k] * tbe_0 + 4.0 * to_xxyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xyzz[k] = 2.0 * to_xxy_xyz[k] - 2.0 * to_xxy_xyzzz[k] * tke_0 - 4.0 * to_xxyzz_xyz[k] * tbe_0 + 4.0 * to_xxyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xzzz[k] = 3.0 * to_xxy_xzz[k] - 2.0 * to_xxy_xzzzz[k] * tke_0 - 6.0 * to_xxyzz_xzz[k] * tbe_0 + 4.0 * to_xxyzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yyyy[k] = -2.0 * to_xxy_yyyyz[k] * tke_0 + 4.0 * to_xxyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yyyz[k] = to_xxy_yyy[k] - 2.0 * to_xxy_yyyzz[k] * tke_0 - 2.0 * to_xxyzz_yyy[k] * tbe_0 + 4.0 * to_xxyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yyzz[k] = 2.0 * to_xxy_yyz[k] - 2.0 * to_xxy_yyzzz[k] * tke_0 - 4.0 * to_xxyzz_yyz[k] * tbe_0 + 4.0 * to_xxyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yzzz[k] = 3.0 * to_xxy_yzz[k] - 2.0 * to_xxy_yzzzz[k] * tke_0 - 6.0 * to_xxyzz_yzz[k] * tbe_0 + 4.0 * to_xxyzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_zzzz[k] = 4.0 * to_xxy_zzz[k] - 2.0 * to_xxy_zzzzz[k] * tke_0 - 8.0 * to_xxyzz_zzz[k] * tbe_0 + 4.0 * to_xxyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1875-1890 components of targeted buffer : GG

        auto to_z_z_xxzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 75);

        auto to_z_z_xxzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 76);

        auto to_z_z_xxzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 77);

        auto to_z_z_xxzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 78);

        auto to_z_z_xxzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 79);

        auto to_z_z_xxzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 80);

        auto to_z_z_xxzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 81);

        auto to_z_z_xxzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 82);

        auto to_z_z_xxzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 83);

        auto to_z_z_xxzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 84);

        auto to_z_z_xxzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 85);

        auto to_z_z_xxzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 86);

        auto to_z_z_xxzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 87);

        auto to_z_z_xxzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 88);

        auto to_z_z_xxzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_xxz_xxx, to_xxz_xxxxz, to_xxz_xxxyz, to_xxz_xxxzz, to_xxz_xxy, to_xxz_xxyyz, to_xxz_xxyzz, to_xxz_xxz, to_xxz_xxzzz, to_xxz_xyy, to_xxz_xyyyz, to_xxz_xyyzz, to_xxz_xyz, to_xxz_xyzzz, to_xxz_xzz, to_xxz_xzzzz, to_xxz_yyy, to_xxz_yyyyz, to_xxz_yyyzz, to_xxz_yyz, to_xxz_yyzzz, to_xxz_yzz, to_xxz_yzzzz, to_xxz_zzz, to_xxz_zzzzz, to_xxzzz_xxx, to_xxzzz_xxxxz, to_xxzzz_xxxyz, to_xxzzz_xxxzz, to_xxzzz_xxy, to_xxzzz_xxyyz, to_xxzzz_xxyzz, to_xxzzz_xxz, to_xxzzz_xxzzz, to_xxzzz_xyy, to_xxzzz_xyyyz, to_xxzzz_xyyzz, to_xxzzz_xyz, to_xxzzz_xyzzz, to_xxzzz_xzz, to_xxzzz_xzzzz, to_xxzzz_yyy, to_xxzzz_yyyyz, to_xxzzz_yyyzz, to_xxzzz_yyz, to_xxzzz_yyzzz, to_xxzzz_yzz, to_xxzzz_yzzzz, to_xxzzz_zzz, to_xxzzz_zzzzz, to_z_z_xxzz_xxxx, to_z_z_xxzz_xxxy, to_z_z_xxzz_xxxz, to_z_z_xxzz_xxyy, to_z_z_xxzz_xxyz, to_z_z_xxzz_xxzz, to_z_z_xxzz_xyyy, to_z_z_xxzz_xyyz, to_z_z_xxzz_xyzz, to_z_z_xxzz_xzzz, to_z_z_xxzz_yyyy, to_z_z_xxzz_yyyz, to_z_z_xxzz_yyzz, to_z_z_xxzz_yzzz, to_z_z_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxzz_xxxx[k] = -4.0 * to_xxz_xxxxz[k] * tke_0 + 4.0 * to_xxzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxxy[k] = -4.0 * to_xxz_xxxyz[k] * tke_0 + 4.0 * to_xxzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxxz[k] = 2.0 * to_xxz_xxx[k] - 4.0 * to_xxz_xxxzz[k] * tke_0 - 2.0 * to_xxzzz_xxx[k] * tbe_0 + 4.0 * to_xxzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxyy[k] = -4.0 * to_xxz_xxyyz[k] * tke_0 + 4.0 * to_xxzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxyz[k] = 2.0 * to_xxz_xxy[k] - 4.0 * to_xxz_xxyzz[k] * tke_0 - 2.0 * to_xxzzz_xxy[k] * tbe_0 + 4.0 * to_xxzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxzz[k] = 4.0 * to_xxz_xxz[k] - 4.0 * to_xxz_xxzzz[k] * tke_0 - 4.0 * to_xxzzz_xxz[k] * tbe_0 + 4.0 * to_xxzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xyyy[k] = -4.0 * to_xxz_xyyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xyyz[k] = 2.0 * to_xxz_xyy[k] - 4.0 * to_xxz_xyyzz[k] * tke_0 - 2.0 * to_xxzzz_xyy[k] * tbe_0 + 4.0 * to_xxzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xyzz[k] = 4.0 * to_xxz_xyz[k] - 4.0 * to_xxz_xyzzz[k] * tke_0 - 4.0 * to_xxzzz_xyz[k] * tbe_0 + 4.0 * to_xxzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xzzz[k] = 6.0 * to_xxz_xzz[k] - 4.0 * to_xxz_xzzzz[k] * tke_0 - 6.0 * to_xxzzz_xzz[k] * tbe_0 + 4.0 * to_xxzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yyyy[k] = -4.0 * to_xxz_yyyyz[k] * tke_0 + 4.0 * to_xxzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yyyz[k] = 2.0 * to_xxz_yyy[k] - 4.0 * to_xxz_yyyzz[k] * tke_0 - 2.0 * to_xxzzz_yyy[k] * tbe_0 + 4.0 * to_xxzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yyzz[k] = 4.0 * to_xxz_yyz[k] - 4.0 * to_xxz_yyzzz[k] * tke_0 - 4.0 * to_xxzzz_yyz[k] * tbe_0 + 4.0 * to_xxzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yzzz[k] = 6.0 * to_xxz_yzz[k] - 4.0 * to_xxz_yzzzz[k] * tke_0 - 6.0 * to_xxzzz_yzz[k] * tbe_0 + 4.0 * to_xxzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_zzzz[k] = 8.0 * to_xxz_zzz[k] - 4.0 * to_xxz_zzzzz[k] * tke_0 - 8.0 * to_xxzzz_zzz[k] * tbe_0 + 4.0 * to_xxzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1890-1905 components of targeted buffer : GG

        auto to_z_z_xyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 90);

        auto to_z_z_xyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 91);

        auto to_z_z_xyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 92);

        auto to_z_z_xyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 93);

        auto to_z_z_xyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 94);

        auto to_z_z_xyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 95);

        auto to_z_z_xyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 96);

        auto to_z_z_xyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 97);

        auto to_z_z_xyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 98);

        auto to_z_z_xyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 99);

        auto to_z_z_xyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 100);

        auto to_z_z_xyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 101);

        auto to_z_z_xyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 102);

        auto to_z_z_xyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 103);

        auto to_z_z_xyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_xyyyz_xxx, to_xyyyz_xxxxz, to_xyyyz_xxxyz, to_xyyyz_xxxzz, to_xyyyz_xxy, to_xyyyz_xxyyz, to_xyyyz_xxyzz, to_xyyyz_xxz, to_xyyyz_xxzzz, to_xyyyz_xyy, to_xyyyz_xyyyz, to_xyyyz_xyyzz, to_xyyyz_xyz, to_xyyyz_xyzzz, to_xyyyz_xzz, to_xyyyz_xzzzz, to_xyyyz_yyy, to_xyyyz_yyyyz, to_xyyyz_yyyzz, to_xyyyz_yyz, to_xyyyz_yyzzz, to_xyyyz_yzz, to_xyyyz_yzzzz, to_xyyyz_zzz, to_xyyyz_zzzzz, to_z_z_xyyy_xxxx, to_z_z_xyyy_xxxy, to_z_z_xyyy_xxxz, to_z_z_xyyy_xxyy, to_z_z_xyyy_xxyz, to_z_z_xyyy_xxzz, to_z_z_xyyy_xyyy, to_z_z_xyyy_xyyz, to_z_z_xyyy_xyzz, to_z_z_xyyy_xzzz, to_z_z_xyyy_yyyy, to_z_z_xyyy_yyyz, to_z_z_xyyy_yyzz, to_z_z_xyyy_yzzz, to_z_z_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyy_xxxx[k] = 4.0 * to_xyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxxy[k] = 4.0 * to_xyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxxz[k] = -2.0 * to_xyyyz_xxx[k] * tbe_0 + 4.0 * to_xyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxyy[k] = 4.0 * to_xyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxyz[k] = -2.0 * to_xyyyz_xxy[k] * tbe_0 + 4.0 * to_xyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxzz[k] = -4.0 * to_xyyyz_xxz[k] * tbe_0 + 4.0 * to_xyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xyyy[k] = 4.0 * to_xyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xyyz[k] = -2.0 * to_xyyyz_xyy[k] * tbe_0 + 4.0 * to_xyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xyzz[k] = -4.0 * to_xyyyz_xyz[k] * tbe_0 + 4.0 * to_xyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xzzz[k] = -6.0 * to_xyyyz_xzz[k] * tbe_0 + 4.0 * to_xyyyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yyyy[k] = 4.0 * to_xyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yyyz[k] = -2.0 * to_xyyyz_yyy[k] * tbe_0 + 4.0 * to_xyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yyzz[k] = -4.0 * to_xyyyz_yyz[k] * tbe_0 + 4.0 * to_xyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yzzz[k] = -6.0 * to_xyyyz_yzz[k] * tbe_0 + 4.0 * to_xyyyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_zzzz[k] = -8.0 * to_xyyyz_zzz[k] * tbe_0 + 4.0 * to_xyyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1905-1920 components of targeted buffer : GG

        auto to_z_z_xyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 105);

        auto to_z_z_xyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 106);

        auto to_z_z_xyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 107);

        auto to_z_z_xyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 108);

        auto to_z_z_xyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 109);

        auto to_z_z_xyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 110);

        auto to_z_z_xyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 111);

        auto to_z_z_xyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 112);

        auto to_z_z_xyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 113);

        auto to_z_z_xyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 114);

        auto to_z_z_xyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 115);

        auto to_z_z_xyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 116);

        auto to_z_z_xyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 117);

        auto to_z_z_xyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 118);

        auto to_z_z_xyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_xyy_xxx, to_xyy_xxxxz, to_xyy_xxxyz, to_xyy_xxxzz, to_xyy_xxy, to_xyy_xxyyz, to_xyy_xxyzz, to_xyy_xxz, to_xyy_xxzzz, to_xyy_xyy, to_xyy_xyyyz, to_xyy_xyyzz, to_xyy_xyz, to_xyy_xyzzz, to_xyy_xzz, to_xyy_xzzzz, to_xyy_yyy, to_xyy_yyyyz, to_xyy_yyyzz, to_xyy_yyz, to_xyy_yyzzz, to_xyy_yzz, to_xyy_yzzzz, to_xyy_zzz, to_xyy_zzzzz, to_xyyzz_xxx, to_xyyzz_xxxxz, to_xyyzz_xxxyz, to_xyyzz_xxxzz, to_xyyzz_xxy, to_xyyzz_xxyyz, to_xyyzz_xxyzz, to_xyyzz_xxz, to_xyyzz_xxzzz, to_xyyzz_xyy, to_xyyzz_xyyyz, to_xyyzz_xyyzz, to_xyyzz_xyz, to_xyyzz_xyzzz, to_xyyzz_xzz, to_xyyzz_xzzzz, to_xyyzz_yyy, to_xyyzz_yyyyz, to_xyyzz_yyyzz, to_xyyzz_yyz, to_xyyzz_yyzzz, to_xyyzz_yzz, to_xyyzz_yzzzz, to_xyyzz_zzz, to_xyyzz_zzzzz, to_z_z_xyyz_xxxx, to_z_z_xyyz_xxxy, to_z_z_xyyz_xxxz, to_z_z_xyyz_xxyy, to_z_z_xyyz_xxyz, to_z_z_xyyz_xxzz, to_z_z_xyyz_xyyy, to_z_z_xyyz_xyyz, to_z_z_xyyz_xyzz, to_z_z_xyyz_xzzz, to_z_z_xyyz_yyyy, to_z_z_xyyz_yyyz, to_z_z_xyyz_yyzz, to_z_z_xyyz_yzzz, to_z_z_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyz_xxxx[k] = -2.0 * to_xyy_xxxxz[k] * tke_0 + 4.0 * to_xyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxxy[k] = -2.0 * to_xyy_xxxyz[k] * tke_0 + 4.0 * to_xyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxxz[k] = to_xyy_xxx[k] - 2.0 * to_xyy_xxxzz[k] * tke_0 - 2.0 * to_xyyzz_xxx[k] * tbe_0 + 4.0 * to_xyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxyy[k] = -2.0 * to_xyy_xxyyz[k] * tke_0 + 4.0 * to_xyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxyz[k] = to_xyy_xxy[k] - 2.0 * to_xyy_xxyzz[k] * tke_0 - 2.0 * to_xyyzz_xxy[k] * tbe_0 + 4.0 * to_xyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxzz[k] = 2.0 * to_xyy_xxz[k] - 2.0 * to_xyy_xxzzz[k] * tke_0 - 4.0 * to_xyyzz_xxz[k] * tbe_0 + 4.0 * to_xyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xyyy[k] = -2.0 * to_xyy_xyyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xyyz[k] = to_xyy_xyy[k] - 2.0 * to_xyy_xyyzz[k] * tke_0 - 2.0 * to_xyyzz_xyy[k] * tbe_0 + 4.0 * to_xyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xyzz[k] = 2.0 * to_xyy_xyz[k] - 2.0 * to_xyy_xyzzz[k] * tke_0 - 4.0 * to_xyyzz_xyz[k] * tbe_0 + 4.0 * to_xyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xzzz[k] = 3.0 * to_xyy_xzz[k] - 2.0 * to_xyy_xzzzz[k] * tke_0 - 6.0 * to_xyyzz_xzz[k] * tbe_0 + 4.0 * to_xyyzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yyyy[k] = -2.0 * to_xyy_yyyyz[k] * tke_0 + 4.0 * to_xyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yyyz[k] = to_xyy_yyy[k] - 2.0 * to_xyy_yyyzz[k] * tke_0 - 2.0 * to_xyyzz_yyy[k] * tbe_0 + 4.0 * to_xyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yyzz[k] = 2.0 * to_xyy_yyz[k] - 2.0 * to_xyy_yyzzz[k] * tke_0 - 4.0 * to_xyyzz_yyz[k] * tbe_0 + 4.0 * to_xyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yzzz[k] = 3.0 * to_xyy_yzz[k] - 2.0 * to_xyy_yzzzz[k] * tke_0 - 6.0 * to_xyyzz_yzz[k] * tbe_0 + 4.0 * to_xyyzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_zzzz[k] = 4.0 * to_xyy_zzz[k] - 2.0 * to_xyy_zzzzz[k] * tke_0 - 8.0 * to_xyyzz_zzz[k] * tbe_0 + 4.0 * to_xyyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1920-1935 components of targeted buffer : GG

        auto to_z_z_xyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 120);

        auto to_z_z_xyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 121);

        auto to_z_z_xyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 122);

        auto to_z_z_xyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 123);

        auto to_z_z_xyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 124);

        auto to_z_z_xyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 125);

        auto to_z_z_xyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 126);

        auto to_z_z_xyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 127);

        auto to_z_z_xyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 128);

        auto to_z_z_xyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 129);

        auto to_z_z_xyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 130);

        auto to_z_z_xyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 131);

        auto to_z_z_xyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 132);

        auto to_z_z_xyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 133);

        auto to_z_z_xyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_xyz_xxx, to_xyz_xxxxz, to_xyz_xxxyz, to_xyz_xxxzz, to_xyz_xxy, to_xyz_xxyyz, to_xyz_xxyzz, to_xyz_xxz, to_xyz_xxzzz, to_xyz_xyy, to_xyz_xyyyz, to_xyz_xyyzz, to_xyz_xyz, to_xyz_xyzzz, to_xyz_xzz, to_xyz_xzzzz, to_xyz_yyy, to_xyz_yyyyz, to_xyz_yyyzz, to_xyz_yyz, to_xyz_yyzzz, to_xyz_yzz, to_xyz_yzzzz, to_xyz_zzz, to_xyz_zzzzz, to_xyzzz_xxx, to_xyzzz_xxxxz, to_xyzzz_xxxyz, to_xyzzz_xxxzz, to_xyzzz_xxy, to_xyzzz_xxyyz, to_xyzzz_xxyzz, to_xyzzz_xxz, to_xyzzz_xxzzz, to_xyzzz_xyy, to_xyzzz_xyyyz, to_xyzzz_xyyzz, to_xyzzz_xyz, to_xyzzz_xyzzz, to_xyzzz_xzz, to_xyzzz_xzzzz, to_xyzzz_yyy, to_xyzzz_yyyyz, to_xyzzz_yyyzz, to_xyzzz_yyz, to_xyzzz_yyzzz, to_xyzzz_yzz, to_xyzzz_yzzzz, to_xyzzz_zzz, to_xyzzz_zzzzz, to_z_z_xyzz_xxxx, to_z_z_xyzz_xxxy, to_z_z_xyzz_xxxz, to_z_z_xyzz_xxyy, to_z_z_xyzz_xxyz, to_z_z_xyzz_xxzz, to_z_z_xyzz_xyyy, to_z_z_xyzz_xyyz, to_z_z_xyzz_xyzz, to_z_z_xyzz_xzzz, to_z_z_xyzz_yyyy, to_z_z_xyzz_yyyz, to_z_z_xyzz_yyzz, to_z_z_xyzz_yzzz, to_z_z_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyzz_xxxx[k] = -4.0 * to_xyz_xxxxz[k] * tke_0 + 4.0 * to_xyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxxy[k] = -4.0 * to_xyz_xxxyz[k] * tke_0 + 4.0 * to_xyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxxz[k] = 2.0 * to_xyz_xxx[k] - 4.0 * to_xyz_xxxzz[k] * tke_0 - 2.0 * to_xyzzz_xxx[k] * tbe_0 + 4.0 * to_xyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxyy[k] = -4.0 * to_xyz_xxyyz[k] * tke_0 + 4.0 * to_xyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxyz[k] = 2.0 * to_xyz_xxy[k] - 4.0 * to_xyz_xxyzz[k] * tke_0 - 2.0 * to_xyzzz_xxy[k] * tbe_0 + 4.0 * to_xyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxzz[k] = 4.0 * to_xyz_xxz[k] - 4.0 * to_xyz_xxzzz[k] * tke_0 - 4.0 * to_xyzzz_xxz[k] * tbe_0 + 4.0 * to_xyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xyyy[k] = -4.0 * to_xyz_xyyyz[k] * tke_0 + 4.0 * to_xyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xyyz[k] = 2.0 * to_xyz_xyy[k] - 4.0 * to_xyz_xyyzz[k] * tke_0 - 2.0 * to_xyzzz_xyy[k] * tbe_0 + 4.0 * to_xyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xyzz[k] = 4.0 * to_xyz_xyz[k] - 4.0 * to_xyz_xyzzz[k] * tke_0 - 4.0 * to_xyzzz_xyz[k] * tbe_0 + 4.0 * to_xyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xzzz[k] = 6.0 * to_xyz_xzz[k] - 4.0 * to_xyz_xzzzz[k] * tke_0 - 6.0 * to_xyzzz_xzz[k] * tbe_0 + 4.0 * to_xyzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yyyy[k] = -4.0 * to_xyz_yyyyz[k] * tke_0 + 4.0 * to_xyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yyyz[k] = 2.0 * to_xyz_yyy[k] - 4.0 * to_xyz_yyyzz[k] * tke_0 - 2.0 * to_xyzzz_yyy[k] * tbe_0 + 4.0 * to_xyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yyzz[k] = 4.0 * to_xyz_yyz[k] - 4.0 * to_xyz_yyzzz[k] * tke_0 - 4.0 * to_xyzzz_yyz[k] * tbe_0 + 4.0 * to_xyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yzzz[k] = 6.0 * to_xyz_yzz[k] - 4.0 * to_xyz_yzzzz[k] * tke_0 - 6.0 * to_xyzzz_yzz[k] * tbe_0 + 4.0 * to_xyzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_zzzz[k] = 8.0 * to_xyz_zzz[k] - 4.0 * to_xyz_zzzzz[k] * tke_0 - 8.0 * to_xyzzz_zzz[k] * tbe_0 + 4.0 * to_xyzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1935-1950 components of targeted buffer : GG

        auto to_z_z_xzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 135);

        auto to_z_z_xzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 136);

        auto to_z_z_xzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 137);

        auto to_z_z_xzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 138);

        auto to_z_z_xzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 139);

        auto to_z_z_xzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 140);

        auto to_z_z_xzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 141);

        auto to_z_z_xzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 142);

        auto to_z_z_xzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 143);

        auto to_z_z_xzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 144);

        auto to_z_z_xzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 145);

        auto to_z_z_xzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 146);

        auto to_z_z_xzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 147);

        auto to_z_z_xzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 148);

        auto to_z_z_xzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_xzz_xxx, to_xzz_xxxxz, to_xzz_xxxyz, to_xzz_xxxzz, to_xzz_xxy, to_xzz_xxyyz, to_xzz_xxyzz, to_xzz_xxz, to_xzz_xxzzz, to_xzz_xyy, to_xzz_xyyyz, to_xzz_xyyzz, to_xzz_xyz, to_xzz_xyzzz, to_xzz_xzz, to_xzz_xzzzz, to_xzz_yyy, to_xzz_yyyyz, to_xzz_yyyzz, to_xzz_yyz, to_xzz_yyzzz, to_xzz_yzz, to_xzz_yzzzz, to_xzz_zzz, to_xzz_zzzzz, to_xzzzz_xxx, to_xzzzz_xxxxz, to_xzzzz_xxxyz, to_xzzzz_xxxzz, to_xzzzz_xxy, to_xzzzz_xxyyz, to_xzzzz_xxyzz, to_xzzzz_xxz, to_xzzzz_xxzzz, to_xzzzz_xyy, to_xzzzz_xyyyz, to_xzzzz_xyyzz, to_xzzzz_xyz, to_xzzzz_xyzzz, to_xzzzz_xzz, to_xzzzz_xzzzz, to_xzzzz_yyy, to_xzzzz_yyyyz, to_xzzzz_yyyzz, to_xzzzz_yyz, to_xzzzz_yyzzz, to_xzzzz_yzz, to_xzzzz_yzzzz, to_xzzzz_zzz, to_xzzzz_zzzzz, to_z_z_xzzz_xxxx, to_z_z_xzzz_xxxy, to_z_z_xzzz_xxxz, to_z_z_xzzz_xxyy, to_z_z_xzzz_xxyz, to_z_z_xzzz_xxzz, to_z_z_xzzz_xyyy, to_z_z_xzzz_xyyz, to_z_z_xzzz_xyzz, to_z_z_xzzz_xzzz, to_z_z_xzzz_yyyy, to_z_z_xzzz_yyyz, to_z_z_xzzz_yyzz, to_z_z_xzzz_yzzz, to_z_z_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzzz_xxxx[k] = -6.0 * to_xzz_xxxxz[k] * tke_0 + 4.0 * to_xzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxxy[k] = -6.0 * to_xzz_xxxyz[k] * tke_0 + 4.0 * to_xzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxxz[k] = 3.0 * to_xzz_xxx[k] - 6.0 * to_xzz_xxxzz[k] * tke_0 - 2.0 * to_xzzzz_xxx[k] * tbe_0 + 4.0 * to_xzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxyy[k] = -6.0 * to_xzz_xxyyz[k] * tke_0 + 4.0 * to_xzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxyz[k] = 3.0 * to_xzz_xxy[k] - 6.0 * to_xzz_xxyzz[k] * tke_0 - 2.0 * to_xzzzz_xxy[k] * tbe_0 + 4.0 * to_xzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxzz[k] = 6.0 * to_xzz_xxz[k] - 6.0 * to_xzz_xxzzz[k] * tke_0 - 4.0 * to_xzzzz_xxz[k] * tbe_0 + 4.0 * to_xzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xyyy[k] = -6.0 * to_xzz_xyyyz[k] * tke_0 + 4.0 * to_xzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xyyz[k] = 3.0 * to_xzz_xyy[k] - 6.0 * to_xzz_xyyzz[k] * tke_0 - 2.0 * to_xzzzz_xyy[k] * tbe_0 + 4.0 * to_xzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xyzz[k] = 6.0 * to_xzz_xyz[k] - 6.0 * to_xzz_xyzzz[k] * tke_0 - 4.0 * to_xzzzz_xyz[k] * tbe_0 + 4.0 * to_xzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xzzz[k] = 9.0 * to_xzz_xzz[k] - 6.0 * to_xzz_xzzzz[k] * tke_0 - 6.0 * to_xzzzz_xzz[k] * tbe_0 + 4.0 * to_xzzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yyyy[k] = -6.0 * to_xzz_yyyyz[k] * tke_0 + 4.0 * to_xzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yyyz[k] = 3.0 * to_xzz_yyy[k] - 6.0 * to_xzz_yyyzz[k] * tke_0 - 2.0 * to_xzzzz_yyy[k] * tbe_0 + 4.0 * to_xzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yyzz[k] = 6.0 * to_xzz_yyz[k] - 6.0 * to_xzz_yyzzz[k] * tke_0 - 4.0 * to_xzzzz_yyz[k] * tbe_0 + 4.0 * to_xzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yzzz[k] = 9.0 * to_xzz_yzz[k] - 6.0 * to_xzz_yzzzz[k] * tke_0 - 6.0 * to_xzzzz_yzz[k] * tbe_0 + 4.0 * to_xzzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_zzzz[k] = 12.0 * to_xzz_zzz[k] - 6.0 * to_xzz_zzzzz[k] * tke_0 - 8.0 * to_xzzzz_zzz[k] * tbe_0 + 4.0 * to_xzzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1950-1965 components of targeted buffer : GG

        auto to_z_z_yyyy_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 150);

        auto to_z_z_yyyy_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 151);

        auto to_z_z_yyyy_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 152);

        auto to_z_z_yyyy_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 153);

        auto to_z_z_yyyy_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 154);

        auto to_z_z_yyyy_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 155);

        auto to_z_z_yyyy_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 156);

        auto to_z_z_yyyy_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 157);

        auto to_z_z_yyyy_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 158);

        auto to_z_z_yyyy_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 159);

        auto to_z_z_yyyy_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 160);

        auto to_z_z_yyyy_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 161);

        auto to_z_z_yyyy_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 162);

        auto to_z_z_yyyy_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 163);

        auto to_z_z_yyyy_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_yyyyz_xxx, to_yyyyz_xxxxz, to_yyyyz_xxxyz, to_yyyyz_xxxzz, to_yyyyz_xxy, to_yyyyz_xxyyz, to_yyyyz_xxyzz, to_yyyyz_xxz, to_yyyyz_xxzzz, to_yyyyz_xyy, to_yyyyz_xyyyz, to_yyyyz_xyyzz, to_yyyyz_xyz, to_yyyyz_xyzzz, to_yyyyz_xzz, to_yyyyz_xzzzz, to_yyyyz_yyy, to_yyyyz_yyyyz, to_yyyyz_yyyzz, to_yyyyz_yyz, to_yyyyz_yyzzz, to_yyyyz_yzz, to_yyyyz_yzzzz, to_yyyyz_zzz, to_yyyyz_zzzzz, to_z_z_yyyy_xxxx, to_z_z_yyyy_xxxy, to_z_z_yyyy_xxxz, to_z_z_yyyy_xxyy, to_z_z_yyyy_xxyz, to_z_z_yyyy_xxzz, to_z_z_yyyy_xyyy, to_z_z_yyyy_xyyz, to_z_z_yyyy_xyzz, to_z_z_yyyy_xzzz, to_z_z_yyyy_yyyy, to_z_z_yyyy_yyyz, to_z_z_yyyy_yyzz, to_z_z_yyyy_yzzz, to_z_z_yyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyy_xxxx[k] = 4.0 * to_yyyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxxy[k] = 4.0 * to_yyyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxxz[k] = -2.0 * to_yyyyz_xxx[k] * tbe_0 + 4.0 * to_yyyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxyy[k] = 4.0 * to_yyyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxyz[k] = -2.0 * to_yyyyz_xxy[k] * tbe_0 + 4.0 * to_yyyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxzz[k] = -4.0 * to_yyyyz_xxz[k] * tbe_0 + 4.0 * to_yyyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xyyy[k] = 4.0 * to_yyyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xyyz[k] = -2.0 * to_yyyyz_xyy[k] * tbe_0 + 4.0 * to_yyyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xyzz[k] = -4.0 * to_yyyyz_xyz[k] * tbe_0 + 4.0 * to_yyyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xzzz[k] = -6.0 * to_yyyyz_xzz[k] * tbe_0 + 4.0 * to_yyyyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yyyy[k] = 4.0 * to_yyyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yyyz[k] = -2.0 * to_yyyyz_yyy[k] * tbe_0 + 4.0 * to_yyyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yyzz[k] = -4.0 * to_yyyyz_yyz[k] * tbe_0 + 4.0 * to_yyyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yzzz[k] = -6.0 * to_yyyyz_yzz[k] * tbe_0 + 4.0 * to_yyyyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_zzzz[k] = -8.0 * to_yyyyz_zzz[k] * tbe_0 + 4.0 * to_yyyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1965-1980 components of targeted buffer : GG

        auto to_z_z_yyyz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 165);

        auto to_z_z_yyyz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 166);

        auto to_z_z_yyyz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 167);

        auto to_z_z_yyyz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 168);

        auto to_z_z_yyyz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 169);

        auto to_z_z_yyyz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 170);

        auto to_z_z_yyyz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 171);

        auto to_z_z_yyyz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 172);

        auto to_z_z_yyyz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 173);

        auto to_z_z_yyyz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 174);

        auto to_z_z_yyyz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 175);

        auto to_z_z_yyyz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 176);

        auto to_z_z_yyyz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 177);

        auto to_z_z_yyyz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 178);

        auto to_z_z_yyyz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_yyy_xxx, to_yyy_xxxxz, to_yyy_xxxyz, to_yyy_xxxzz, to_yyy_xxy, to_yyy_xxyyz, to_yyy_xxyzz, to_yyy_xxz, to_yyy_xxzzz, to_yyy_xyy, to_yyy_xyyyz, to_yyy_xyyzz, to_yyy_xyz, to_yyy_xyzzz, to_yyy_xzz, to_yyy_xzzzz, to_yyy_yyy, to_yyy_yyyyz, to_yyy_yyyzz, to_yyy_yyz, to_yyy_yyzzz, to_yyy_yzz, to_yyy_yzzzz, to_yyy_zzz, to_yyy_zzzzz, to_yyyzz_xxx, to_yyyzz_xxxxz, to_yyyzz_xxxyz, to_yyyzz_xxxzz, to_yyyzz_xxy, to_yyyzz_xxyyz, to_yyyzz_xxyzz, to_yyyzz_xxz, to_yyyzz_xxzzz, to_yyyzz_xyy, to_yyyzz_xyyyz, to_yyyzz_xyyzz, to_yyyzz_xyz, to_yyyzz_xyzzz, to_yyyzz_xzz, to_yyyzz_xzzzz, to_yyyzz_yyy, to_yyyzz_yyyyz, to_yyyzz_yyyzz, to_yyyzz_yyz, to_yyyzz_yyzzz, to_yyyzz_yzz, to_yyyzz_yzzzz, to_yyyzz_zzz, to_yyyzz_zzzzz, to_z_z_yyyz_xxxx, to_z_z_yyyz_xxxy, to_z_z_yyyz_xxxz, to_z_z_yyyz_xxyy, to_z_z_yyyz_xxyz, to_z_z_yyyz_xxzz, to_z_z_yyyz_xyyy, to_z_z_yyyz_xyyz, to_z_z_yyyz_xyzz, to_z_z_yyyz_xzzz, to_z_z_yyyz_yyyy, to_z_z_yyyz_yyyz, to_z_z_yyyz_yyzz, to_z_z_yyyz_yzzz, to_z_z_yyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyz_xxxx[k] = -2.0 * to_yyy_xxxxz[k] * tke_0 + 4.0 * to_yyyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxxy[k] = -2.0 * to_yyy_xxxyz[k] * tke_0 + 4.0 * to_yyyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxxz[k] = to_yyy_xxx[k] - 2.0 * to_yyy_xxxzz[k] * tke_0 - 2.0 * to_yyyzz_xxx[k] * tbe_0 + 4.0 * to_yyyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxyy[k] = -2.0 * to_yyy_xxyyz[k] * tke_0 + 4.0 * to_yyyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxyz[k] = to_yyy_xxy[k] - 2.0 * to_yyy_xxyzz[k] * tke_0 - 2.0 * to_yyyzz_xxy[k] * tbe_0 + 4.0 * to_yyyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxzz[k] = 2.0 * to_yyy_xxz[k] - 2.0 * to_yyy_xxzzz[k] * tke_0 - 4.0 * to_yyyzz_xxz[k] * tbe_0 + 4.0 * to_yyyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xyyy[k] = -2.0 * to_yyy_xyyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xyyz[k] = to_yyy_xyy[k] - 2.0 * to_yyy_xyyzz[k] * tke_0 - 2.0 * to_yyyzz_xyy[k] * tbe_0 + 4.0 * to_yyyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xyzz[k] = 2.0 * to_yyy_xyz[k] - 2.0 * to_yyy_xyzzz[k] * tke_0 - 4.0 * to_yyyzz_xyz[k] * tbe_0 + 4.0 * to_yyyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xzzz[k] = 3.0 * to_yyy_xzz[k] - 2.0 * to_yyy_xzzzz[k] * tke_0 - 6.0 * to_yyyzz_xzz[k] * tbe_0 + 4.0 * to_yyyzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yyyy[k] = -2.0 * to_yyy_yyyyz[k] * tke_0 + 4.0 * to_yyyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yyyz[k] = to_yyy_yyy[k] - 2.0 * to_yyy_yyyzz[k] * tke_0 - 2.0 * to_yyyzz_yyy[k] * tbe_0 + 4.0 * to_yyyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yyzz[k] = 2.0 * to_yyy_yyz[k] - 2.0 * to_yyy_yyzzz[k] * tke_0 - 4.0 * to_yyyzz_yyz[k] * tbe_0 + 4.0 * to_yyyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yzzz[k] = 3.0 * to_yyy_yzz[k] - 2.0 * to_yyy_yzzzz[k] * tke_0 - 6.0 * to_yyyzz_yzz[k] * tbe_0 + 4.0 * to_yyyzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_zzzz[k] = 4.0 * to_yyy_zzz[k] - 2.0 * to_yyy_zzzzz[k] * tke_0 - 8.0 * to_yyyzz_zzz[k] * tbe_0 + 4.0 * to_yyyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1980-1995 components of targeted buffer : GG

        auto to_z_z_yyzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 180);

        auto to_z_z_yyzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 181);

        auto to_z_z_yyzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 182);

        auto to_z_z_yyzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 183);

        auto to_z_z_yyzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 184);

        auto to_z_z_yyzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 185);

        auto to_z_z_yyzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 186);

        auto to_z_z_yyzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 187);

        auto to_z_z_yyzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 188);

        auto to_z_z_yyzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 189);

        auto to_z_z_yyzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 190);

        auto to_z_z_yyzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 191);

        auto to_z_z_yyzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 192);

        auto to_z_z_yyzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 193);

        auto to_z_z_yyzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_yyz_xxx, to_yyz_xxxxz, to_yyz_xxxyz, to_yyz_xxxzz, to_yyz_xxy, to_yyz_xxyyz, to_yyz_xxyzz, to_yyz_xxz, to_yyz_xxzzz, to_yyz_xyy, to_yyz_xyyyz, to_yyz_xyyzz, to_yyz_xyz, to_yyz_xyzzz, to_yyz_xzz, to_yyz_xzzzz, to_yyz_yyy, to_yyz_yyyyz, to_yyz_yyyzz, to_yyz_yyz, to_yyz_yyzzz, to_yyz_yzz, to_yyz_yzzzz, to_yyz_zzz, to_yyz_zzzzz, to_yyzzz_xxx, to_yyzzz_xxxxz, to_yyzzz_xxxyz, to_yyzzz_xxxzz, to_yyzzz_xxy, to_yyzzz_xxyyz, to_yyzzz_xxyzz, to_yyzzz_xxz, to_yyzzz_xxzzz, to_yyzzz_xyy, to_yyzzz_xyyyz, to_yyzzz_xyyzz, to_yyzzz_xyz, to_yyzzz_xyzzz, to_yyzzz_xzz, to_yyzzz_xzzzz, to_yyzzz_yyy, to_yyzzz_yyyyz, to_yyzzz_yyyzz, to_yyzzz_yyz, to_yyzzz_yyzzz, to_yyzzz_yzz, to_yyzzz_yzzzz, to_yyzzz_zzz, to_yyzzz_zzzzz, to_z_z_yyzz_xxxx, to_z_z_yyzz_xxxy, to_z_z_yyzz_xxxz, to_z_z_yyzz_xxyy, to_z_z_yyzz_xxyz, to_z_z_yyzz_xxzz, to_z_z_yyzz_xyyy, to_z_z_yyzz_xyyz, to_z_z_yyzz_xyzz, to_z_z_yyzz_xzzz, to_z_z_yyzz_yyyy, to_z_z_yyzz_yyyz, to_z_z_yyzz_yyzz, to_z_z_yyzz_yzzz, to_z_z_yyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyzz_xxxx[k] = -4.0 * to_yyz_xxxxz[k] * tke_0 + 4.0 * to_yyzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxxy[k] = -4.0 * to_yyz_xxxyz[k] * tke_0 + 4.0 * to_yyzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxxz[k] = 2.0 * to_yyz_xxx[k] - 4.0 * to_yyz_xxxzz[k] * tke_0 - 2.0 * to_yyzzz_xxx[k] * tbe_0 + 4.0 * to_yyzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxyy[k] = -4.0 * to_yyz_xxyyz[k] * tke_0 + 4.0 * to_yyzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxyz[k] = 2.0 * to_yyz_xxy[k] - 4.0 * to_yyz_xxyzz[k] * tke_0 - 2.0 * to_yyzzz_xxy[k] * tbe_0 + 4.0 * to_yyzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxzz[k] = 4.0 * to_yyz_xxz[k] - 4.0 * to_yyz_xxzzz[k] * tke_0 - 4.0 * to_yyzzz_xxz[k] * tbe_0 + 4.0 * to_yyzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xyyy[k] = -4.0 * to_yyz_xyyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xyyz[k] = 2.0 * to_yyz_xyy[k] - 4.0 * to_yyz_xyyzz[k] * tke_0 - 2.0 * to_yyzzz_xyy[k] * tbe_0 + 4.0 * to_yyzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xyzz[k] = 4.0 * to_yyz_xyz[k] - 4.0 * to_yyz_xyzzz[k] * tke_0 - 4.0 * to_yyzzz_xyz[k] * tbe_0 + 4.0 * to_yyzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xzzz[k] = 6.0 * to_yyz_xzz[k] - 4.0 * to_yyz_xzzzz[k] * tke_0 - 6.0 * to_yyzzz_xzz[k] * tbe_0 + 4.0 * to_yyzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yyyy[k] = -4.0 * to_yyz_yyyyz[k] * tke_0 + 4.0 * to_yyzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yyyz[k] = 2.0 * to_yyz_yyy[k] - 4.0 * to_yyz_yyyzz[k] * tke_0 - 2.0 * to_yyzzz_yyy[k] * tbe_0 + 4.0 * to_yyzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yyzz[k] = 4.0 * to_yyz_yyz[k] - 4.0 * to_yyz_yyzzz[k] * tke_0 - 4.0 * to_yyzzz_yyz[k] * tbe_0 + 4.0 * to_yyzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yzzz[k] = 6.0 * to_yyz_yzz[k] - 4.0 * to_yyz_yzzzz[k] * tke_0 - 6.0 * to_yyzzz_yzz[k] * tbe_0 + 4.0 * to_yyzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_zzzz[k] = 8.0 * to_yyz_zzz[k] - 4.0 * to_yyz_zzzzz[k] * tke_0 - 8.0 * to_yyzzz_zzz[k] * tbe_0 + 4.0 * to_yyzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1995-2010 components of targeted buffer : GG

        auto to_z_z_yzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 195);

        auto to_z_z_yzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 196);

        auto to_z_z_yzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 197);

        auto to_z_z_yzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 198);

        auto to_z_z_yzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 199);

        auto to_z_z_yzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 200);

        auto to_z_z_yzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 201);

        auto to_z_z_yzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 202);

        auto to_z_z_yzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 203);

        auto to_z_z_yzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 204);

        auto to_z_z_yzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 205);

        auto to_z_z_yzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 206);

        auto to_z_z_yzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 207);

        auto to_z_z_yzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 208);

        auto to_z_z_yzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_yzz_xxx, to_yzz_xxxxz, to_yzz_xxxyz, to_yzz_xxxzz, to_yzz_xxy, to_yzz_xxyyz, to_yzz_xxyzz, to_yzz_xxz, to_yzz_xxzzz, to_yzz_xyy, to_yzz_xyyyz, to_yzz_xyyzz, to_yzz_xyz, to_yzz_xyzzz, to_yzz_xzz, to_yzz_xzzzz, to_yzz_yyy, to_yzz_yyyyz, to_yzz_yyyzz, to_yzz_yyz, to_yzz_yyzzz, to_yzz_yzz, to_yzz_yzzzz, to_yzz_zzz, to_yzz_zzzzz, to_yzzzz_xxx, to_yzzzz_xxxxz, to_yzzzz_xxxyz, to_yzzzz_xxxzz, to_yzzzz_xxy, to_yzzzz_xxyyz, to_yzzzz_xxyzz, to_yzzzz_xxz, to_yzzzz_xxzzz, to_yzzzz_xyy, to_yzzzz_xyyyz, to_yzzzz_xyyzz, to_yzzzz_xyz, to_yzzzz_xyzzz, to_yzzzz_xzz, to_yzzzz_xzzzz, to_yzzzz_yyy, to_yzzzz_yyyyz, to_yzzzz_yyyzz, to_yzzzz_yyz, to_yzzzz_yyzzz, to_yzzzz_yzz, to_yzzzz_yzzzz, to_yzzzz_zzz, to_yzzzz_zzzzz, to_z_z_yzzz_xxxx, to_z_z_yzzz_xxxy, to_z_z_yzzz_xxxz, to_z_z_yzzz_xxyy, to_z_z_yzzz_xxyz, to_z_z_yzzz_xxzz, to_z_z_yzzz_xyyy, to_z_z_yzzz_xyyz, to_z_z_yzzz_xyzz, to_z_z_yzzz_xzzz, to_z_z_yzzz_yyyy, to_z_z_yzzz_yyyz, to_z_z_yzzz_yyzz, to_z_z_yzzz_yzzz, to_z_z_yzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzzz_xxxx[k] = -6.0 * to_yzz_xxxxz[k] * tke_0 + 4.0 * to_yzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxxy[k] = -6.0 * to_yzz_xxxyz[k] * tke_0 + 4.0 * to_yzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxxz[k] = 3.0 * to_yzz_xxx[k] - 6.0 * to_yzz_xxxzz[k] * tke_0 - 2.0 * to_yzzzz_xxx[k] * tbe_0 + 4.0 * to_yzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxyy[k] = -6.0 * to_yzz_xxyyz[k] * tke_0 + 4.0 * to_yzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxyz[k] = 3.0 * to_yzz_xxy[k] - 6.0 * to_yzz_xxyzz[k] * tke_0 - 2.0 * to_yzzzz_xxy[k] * tbe_0 + 4.0 * to_yzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxzz[k] = 6.0 * to_yzz_xxz[k] - 6.0 * to_yzz_xxzzz[k] * tke_0 - 4.0 * to_yzzzz_xxz[k] * tbe_0 + 4.0 * to_yzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xyyy[k] = -6.0 * to_yzz_xyyyz[k] * tke_0 + 4.0 * to_yzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xyyz[k] = 3.0 * to_yzz_xyy[k] - 6.0 * to_yzz_xyyzz[k] * tke_0 - 2.0 * to_yzzzz_xyy[k] * tbe_0 + 4.0 * to_yzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xyzz[k] = 6.0 * to_yzz_xyz[k] - 6.0 * to_yzz_xyzzz[k] * tke_0 - 4.0 * to_yzzzz_xyz[k] * tbe_0 + 4.0 * to_yzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xzzz[k] = 9.0 * to_yzz_xzz[k] - 6.0 * to_yzz_xzzzz[k] * tke_0 - 6.0 * to_yzzzz_xzz[k] * tbe_0 + 4.0 * to_yzzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yyyy[k] = -6.0 * to_yzz_yyyyz[k] * tke_0 + 4.0 * to_yzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yyyz[k] = 3.0 * to_yzz_yyy[k] - 6.0 * to_yzz_yyyzz[k] * tke_0 - 2.0 * to_yzzzz_yyy[k] * tbe_0 + 4.0 * to_yzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yyzz[k] = 6.0 * to_yzz_yyz[k] - 6.0 * to_yzz_yyzzz[k] * tke_0 - 4.0 * to_yzzzz_yyz[k] * tbe_0 + 4.0 * to_yzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yzzz[k] = 9.0 * to_yzz_yzz[k] - 6.0 * to_yzz_yzzzz[k] * tke_0 - 6.0 * to_yzzzz_yzz[k] * tbe_0 + 4.0 * to_yzzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_zzzz[k] = 12.0 * to_yzz_zzz[k] - 6.0 * to_yzz_zzzzz[k] * tke_0 - 8.0 * to_yzzzz_zzz[k] * tbe_0 + 4.0 * to_yzzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 2010-2025 components of targeted buffer : GG

        auto to_z_z_zzzz_xxxx = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 210);

        auto to_z_z_zzzz_xxxy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 211);

        auto to_z_z_zzzz_xxxz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 212);

        auto to_z_z_zzzz_xxyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 213);

        auto to_z_z_zzzz_xxyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 214);

        auto to_z_z_zzzz_xxzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 215);

        auto to_z_z_zzzz_xyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 216);

        auto to_z_z_zzzz_xyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 217);

        auto to_z_z_zzzz_xyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 218);

        auto to_z_z_zzzz_xzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 219);

        auto to_z_z_zzzz_yyyy = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 220);

        auto to_z_z_zzzz_yyyz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 221);

        auto to_z_z_zzzz_yyzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 222);

        auto to_z_z_zzzz_yzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 223);

        auto to_z_z_zzzz_zzzz = pbuffer.data(idx_op_geom_101_gg + 8 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_z_z_zzzz_xxxx, to_z_z_zzzz_xxxy, to_z_z_zzzz_xxxz, to_z_z_zzzz_xxyy, to_z_z_zzzz_xxyz, to_z_z_zzzz_xxzz, to_z_z_zzzz_xyyy, to_z_z_zzzz_xyyz, to_z_z_zzzz_xyzz, to_z_z_zzzz_xzzz, to_z_z_zzzz_yyyy, to_z_z_zzzz_yyyz, to_z_z_zzzz_yyzz, to_z_z_zzzz_yzzz, to_z_z_zzzz_zzzz, to_zzz_xxx, to_zzz_xxxxz, to_zzz_xxxyz, to_zzz_xxxzz, to_zzz_xxy, to_zzz_xxyyz, to_zzz_xxyzz, to_zzz_xxz, to_zzz_xxzzz, to_zzz_xyy, to_zzz_xyyyz, to_zzz_xyyzz, to_zzz_xyz, to_zzz_xyzzz, to_zzz_xzz, to_zzz_xzzzz, to_zzz_yyy, to_zzz_yyyyz, to_zzz_yyyzz, to_zzz_yyz, to_zzz_yyzzz, to_zzz_yzz, to_zzz_yzzzz, to_zzz_zzz, to_zzz_zzzzz, to_zzzzz_xxx, to_zzzzz_xxxxz, to_zzzzz_xxxyz, to_zzzzz_xxxzz, to_zzzzz_xxy, to_zzzzz_xxyyz, to_zzzzz_xxyzz, to_zzzzz_xxz, to_zzzzz_xxzzz, to_zzzzz_xyy, to_zzzzz_xyyyz, to_zzzzz_xyyzz, to_zzzzz_xyz, to_zzzzz_xyzzz, to_zzzzz_xzz, to_zzzzz_xzzzz, to_zzzzz_yyy, to_zzzzz_yyyyz, to_zzzzz_yyyzz, to_zzzzz_yyz, to_zzzzz_yyzzz, to_zzzzz_yzz, to_zzzzz_yzzzz, to_zzzzz_zzz, to_zzzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzzz_xxxx[k] = -8.0 * to_zzz_xxxxz[k] * tke_0 + 4.0 * to_zzzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxxy[k] = -8.0 * to_zzz_xxxyz[k] * tke_0 + 4.0 * to_zzzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxxz[k] = 4.0 * to_zzz_xxx[k] - 8.0 * to_zzz_xxxzz[k] * tke_0 - 2.0 * to_zzzzz_xxx[k] * tbe_0 + 4.0 * to_zzzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxyy[k] = -8.0 * to_zzz_xxyyz[k] * tke_0 + 4.0 * to_zzzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxyz[k] = 4.0 * to_zzz_xxy[k] - 8.0 * to_zzz_xxyzz[k] * tke_0 - 2.0 * to_zzzzz_xxy[k] * tbe_0 + 4.0 * to_zzzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxzz[k] = 8.0 * to_zzz_xxz[k] - 8.0 * to_zzz_xxzzz[k] * tke_0 - 4.0 * to_zzzzz_xxz[k] * tbe_0 + 4.0 * to_zzzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xyyy[k] = -8.0 * to_zzz_xyyyz[k] * tke_0 + 4.0 * to_zzzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xyyz[k] = 4.0 * to_zzz_xyy[k] - 8.0 * to_zzz_xyyzz[k] * tke_0 - 2.0 * to_zzzzz_xyy[k] * tbe_0 + 4.0 * to_zzzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xyzz[k] = 8.0 * to_zzz_xyz[k] - 8.0 * to_zzz_xyzzz[k] * tke_0 - 4.0 * to_zzzzz_xyz[k] * tbe_0 + 4.0 * to_zzzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xzzz[k] = 12.0 * to_zzz_xzz[k] - 8.0 * to_zzz_xzzzz[k] * tke_0 - 6.0 * to_zzzzz_xzz[k] * tbe_0 + 4.0 * to_zzzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yyyy[k] = -8.0 * to_zzz_yyyyz[k] * tke_0 + 4.0 * to_zzzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yyyz[k] = 4.0 * to_zzz_yyy[k] - 8.0 * to_zzz_yyyzz[k] * tke_0 - 2.0 * to_zzzzz_yyy[k] * tbe_0 + 4.0 * to_zzzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yyzz[k] = 8.0 * to_zzz_yyz[k] - 8.0 * to_zzz_yyzzz[k] * tke_0 - 4.0 * to_zzzzz_yyz[k] * tbe_0 + 4.0 * to_zzzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yzzz[k] = 12.0 * to_zzz_yzz[k] - 8.0 * to_zzz_yzzzz[k] * tke_0 - 6.0 * to_zzzzz_yzz[k] * tbe_0 + 4.0 * to_zzzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_zzzz[k] = 16.0 * to_zzz_zzz[k] - 8.0 * to_zzz_zzzzz[k] * tke_0 - 8.0 * to_zzzzz_zzz[k] * tbe_0 + 4.0 * to_zzzzz_zzzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

