#include "T2CHrrABRecGG.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_gg(CSimdArray<double>& cbuffer, 
            const size_t idx_gg,
            const size_t idx_fg,
            const size_t idx_fh,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : FG

    auto t_xxx_xxxx = cbuffer.data(idx_fg);

    auto t_xxx_xxxy = cbuffer.data(idx_fg + 1);

    auto t_xxx_xxxz = cbuffer.data(idx_fg + 2);

    auto t_xxx_xxyy = cbuffer.data(idx_fg + 3);

    auto t_xxx_xxyz = cbuffer.data(idx_fg + 4);

    auto t_xxx_xxzz = cbuffer.data(idx_fg + 5);

    auto t_xxx_xyyy = cbuffer.data(idx_fg + 6);

    auto t_xxx_xyyz = cbuffer.data(idx_fg + 7);

    auto t_xxx_xyzz = cbuffer.data(idx_fg + 8);

    auto t_xxx_xzzz = cbuffer.data(idx_fg + 9);

    auto t_xxx_yyyy = cbuffer.data(idx_fg + 10);

    auto t_xxx_yyyz = cbuffer.data(idx_fg + 11);

    auto t_xxx_yyzz = cbuffer.data(idx_fg + 12);

    auto t_xxx_yzzz = cbuffer.data(idx_fg + 13);

    auto t_xxx_zzzz = cbuffer.data(idx_fg + 14);

    auto t_xxy_xxxx = cbuffer.data(idx_fg + 15);

    auto t_xxy_xxxy = cbuffer.data(idx_fg + 16);

    auto t_xxy_xxxz = cbuffer.data(idx_fg + 17);

    auto t_xxy_xxyy = cbuffer.data(idx_fg + 18);

    auto t_xxy_xxyz = cbuffer.data(idx_fg + 19);

    auto t_xxy_xxzz = cbuffer.data(idx_fg + 20);

    auto t_xxy_xyyy = cbuffer.data(idx_fg + 21);

    auto t_xxy_xyyz = cbuffer.data(idx_fg + 22);

    auto t_xxy_xyzz = cbuffer.data(idx_fg + 23);

    auto t_xxy_xzzz = cbuffer.data(idx_fg + 24);

    auto t_xxy_yyyy = cbuffer.data(idx_fg + 25);

    auto t_xxy_yyyz = cbuffer.data(idx_fg + 26);

    auto t_xxy_yyzz = cbuffer.data(idx_fg + 27);

    auto t_xxy_yzzz = cbuffer.data(idx_fg + 28);

    auto t_xxy_zzzz = cbuffer.data(idx_fg + 29);

    auto t_xxz_xxxx = cbuffer.data(idx_fg + 30);

    auto t_xxz_xxxy = cbuffer.data(idx_fg + 31);

    auto t_xxz_xxxz = cbuffer.data(idx_fg + 32);

    auto t_xxz_xxyy = cbuffer.data(idx_fg + 33);

    auto t_xxz_xxyz = cbuffer.data(idx_fg + 34);

    auto t_xxz_xxzz = cbuffer.data(idx_fg + 35);

    auto t_xxz_xyyy = cbuffer.data(idx_fg + 36);

    auto t_xxz_xyyz = cbuffer.data(idx_fg + 37);

    auto t_xxz_xyzz = cbuffer.data(idx_fg + 38);

    auto t_xxz_xzzz = cbuffer.data(idx_fg + 39);

    auto t_xxz_yyyy = cbuffer.data(idx_fg + 40);

    auto t_xxz_yyyz = cbuffer.data(idx_fg + 41);

    auto t_xxz_yyzz = cbuffer.data(idx_fg + 42);

    auto t_xxz_yzzz = cbuffer.data(idx_fg + 43);

    auto t_xxz_zzzz = cbuffer.data(idx_fg + 44);

    auto t_xyy_xxxx = cbuffer.data(idx_fg + 45);

    auto t_xyy_xxxy = cbuffer.data(idx_fg + 46);

    auto t_xyy_xxxz = cbuffer.data(idx_fg + 47);

    auto t_xyy_xxyy = cbuffer.data(idx_fg + 48);

    auto t_xyy_xxyz = cbuffer.data(idx_fg + 49);

    auto t_xyy_xxzz = cbuffer.data(idx_fg + 50);

    auto t_xyy_xyyy = cbuffer.data(idx_fg + 51);

    auto t_xyy_xyyz = cbuffer.data(idx_fg + 52);

    auto t_xyy_xyzz = cbuffer.data(idx_fg + 53);

    auto t_xyy_xzzz = cbuffer.data(idx_fg + 54);

    auto t_xyy_yyyy = cbuffer.data(idx_fg + 55);

    auto t_xyy_yyyz = cbuffer.data(idx_fg + 56);

    auto t_xyy_yyzz = cbuffer.data(idx_fg + 57);

    auto t_xyy_yzzz = cbuffer.data(idx_fg + 58);

    auto t_xyy_zzzz = cbuffer.data(idx_fg + 59);

    auto t_xyz_xxxx = cbuffer.data(idx_fg + 60);

    auto t_xyz_xxxy = cbuffer.data(idx_fg + 61);

    auto t_xyz_xxxz = cbuffer.data(idx_fg + 62);

    auto t_xyz_xxyy = cbuffer.data(idx_fg + 63);

    auto t_xyz_xxyz = cbuffer.data(idx_fg + 64);

    auto t_xyz_xxzz = cbuffer.data(idx_fg + 65);

    auto t_xyz_xyyy = cbuffer.data(idx_fg + 66);

    auto t_xyz_xyyz = cbuffer.data(idx_fg + 67);

    auto t_xyz_xyzz = cbuffer.data(idx_fg + 68);

    auto t_xyz_xzzz = cbuffer.data(idx_fg + 69);

    auto t_xyz_yyyy = cbuffer.data(idx_fg + 70);

    auto t_xyz_yyyz = cbuffer.data(idx_fg + 71);

    auto t_xyz_yyzz = cbuffer.data(idx_fg + 72);

    auto t_xyz_yzzz = cbuffer.data(idx_fg + 73);

    auto t_xyz_zzzz = cbuffer.data(idx_fg + 74);

    auto t_xzz_xxxx = cbuffer.data(idx_fg + 75);

    auto t_xzz_xxxy = cbuffer.data(idx_fg + 76);

    auto t_xzz_xxxz = cbuffer.data(idx_fg + 77);

    auto t_xzz_xxyy = cbuffer.data(idx_fg + 78);

    auto t_xzz_xxyz = cbuffer.data(idx_fg + 79);

    auto t_xzz_xxzz = cbuffer.data(idx_fg + 80);

    auto t_xzz_xyyy = cbuffer.data(idx_fg + 81);

    auto t_xzz_xyyz = cbuffer.data(idx_fg + 82);

    auto t_xzz_xyzz = cbuffer.data(idx_fg + 83);

    auto t_xzz_xzzz = cbuffer.data(idx_fg + 84);

    auto t_xzz_yyyy = cbuffer.data(idx_fg + 85);

    auto t_xzz_yyyz = cbuffer.data(idx_fg + 86);

    auto t_xzz_yyzz = cbuffer.data(idx_fg + 87);

    auto t_xzz_yzzz = cbuffer.data(idx_fg + 88);

    auto t_xzz_zzzz = cbuffer.data(idx_fg + 89);

    auto t_yyy_xxxx = cbuffer.data(idx_fg + 90);

    auto t_yyy_xxxy = cbuffer.data(idx_fg + 91);

    auto t_yyy_xxxz = cbuffer.data(idx_fg + 92);

    auto t_yyy_xxyy = cbuffer.data(idx_fg + 93);

    auto t_yyy_xxyz = cbuffer.data(idx_fg + 94);

    auto t_yyy_xxzz = cbuffer.data(idx_fg + 95);

    auto t_yyy_xyyy = cbuffer.data(idx_fg + 96);

    auto t_yyy_xyyz = cbuffer.data(idx_fg + 97);

    auto t_yyy_xyzz = cbuffer.data(idx_fg + 98);

    auto t_yyy_xzzz = cbuffer.data(idx_fg + 99);

    auto t_yyy_yyyy = cbuffer.data(idx_fg + 100);

    auto t_yyy_yyyz = cbuffer.data(idx_fg + 101);

    auto t_yyy_yyzz = cbuffer.data(idx_fg + 102);

    auto t_yyy_yzzz = cbuffer.data(idx_fg + 103);

    auto t_yyy_zzzz = cbuffer.data(idx_fg + 104);

    auto t_yyz_xxxx = cbuffer.data(idx_fg + 105);

    auto t_yyz_xxxy = cbuffer.data(idx_fg + 106);

    auto t_yyz_xxxz = cbuffer.data(idx_fg + 107);

    auto t_yyz_xxyy = cbuffer.data(idx_fg + 108);

    auto t_yyz_xxyz = cbuffer.data(idx_fg + 109);

    auto t_yyz_xxzz = cbuffer.data(idx_fg + 110);

    auto t_yyz_xyyy = cbuffer.data(idx_fg + 111);

    auto t_yyz_xyyz = cbuffer.data(idx_fg + 112);

    auto t_yyz_xyzz = cbuffer.data(idx_fg + 113);

    auto t_yyz_xzzz = cbuffer.data(idx_fg + 114);

    auto t_yyz_yyyy = cbuffer.data(idx_fg + 115);

    auto t_yyz_yyyz = cbuffer.data(idx_fg + 116);

    auto t_yyz_yyzz = cbuffer.data(idx_fg + 117);

    auto t_yyz_yzzz = cbuffer.data(idx_fg + 118);

    auto t_yyz_zzzz = cbuffer.data(idx_fg + 119);

    auto t_yzz_xxxx = cbuffer.data(idx_fg + 120);

    auto t_yzz_xxxy = cbuffer.data(idx_fg + 121);

    auto t_yzz_xxxz = cbuffer.data(idx_fg + 122);

    auto t_yzz_xxyy = cbuffer.data(idx_fg + 123);

    auto t_yzz_xxyz = cbuffer.data(idx_fg + 124);

    auto t_yzz_xxzz = cbuffer.data(idx_fg + 125);

    auto t_yzz_xyyy = cbuffer.data(idx_fg + 126);

    auto t_yzz_xyyz = cbuffer.data(idx_fg + 127);

    auto t_yzz_xyzz = cbuffer.data(idx_fg + 128);

    auto t_yzz_xzzz = cbuffer.data(idx_fg + 129);

    auto t_yzz_yyyy = cbuffer.data(idx_fg + 130);

    auto t_yzz_yyyz = cbuffer.data(idx_fg + 131);

    auto t_yzz_yyzz = cbuffer.data(idx_fg + 132);

    auto t_yzz_yzzz = cbuffer.data(idx_fg + 133);

    auto t_yzz_zzzz = cbuffer.data(idx_fg + 134);

    auto t_zzz_xxxx = cbuffer.data(idx_fg + 135);

    auto t_zzz_xxxy = cbuffer.data(idx_fg + 136);

    auto t_zzz_xxxz = cbuffer.data(idx_fg + 137);

    auto t_zzz_xxyy = cbuffer.data(idx_fg + 138);

    auto t_zzz_xxyz = cbuffer.data(idx_fg + 139);

    auto t_zzz_xxzz = cbuffer.data(idx_fg + 140);

    auto t_zzz_xyyy = cbuffer.data(idx_fg + 141);

    auto t_zzz_xyyz = cbuffer.data(idx_fg + 142);

    auto t_zzz_xyzz = cbuffer.data(idx_fg + 143);

    auto t_zzz_xzzz = cbuffer.data(idx_fg + 144);

    auto t_zzz_yyyy = cbuffer.data(idx_fg + 145);

    auto t_zzz_yyyz = cbuffer.data(idx_fg + 146);

    auto t_zzz_yyzz = cbuffer.data(idx_fg + 147);

    auto t_zzz_yzzz = cbuffer.data(idx_fg + 148);

    auto t_zzz_zzzz = cbuffer.data(idx_fg + 149);

    // Set up components of auxiliary buffer : FH

    auto t_xxx_xxxxx = cbuffer.data(idx_fh);

    auto t_xxx_xxxxy = cbuffer.data(idx_fh + 1);

    auto t_xxx_xxxxz = cbuffer.data(idx_fh + 2);

    auto t_xxx_xxxyy = cbuffer.data(idx_fh + 3);

    auto t_xxx_xxxyz = cbuffer.data(idx_fh + 4);

    auto t_xxx_xxxzz = cbuffer.data(idx_fh + 5);

    auto t_xxx_xxyyy = cbuffer.data(idx_fh + 6);

    auto t_xxx_xxyyz = cbuffer.data(idx_fh + 7);

    auto t_xxx_xxyzz = cbuffer.data(idx_fh + 8);

    auto t_xxx_xxzzz = cbuffer.data(idx_fh + 9);

    auto t_xxx_xyyyy = cbuffer.data(idx_fh + 10);

    auto t_xxx_xyyyz = cbuffer.data(idx_fh + 11);

    auto t_xxx_xyyzz = cbuffer.data(idx_fh + 12);

    auto t_xxx_xyzzz = cbuffer.data(idx_fh + 13);

    auto t_xxx_xzzzz = cbuffer.data(idx_fh + 14);

    auto t_xxy_xxxxx = cbuffer.data(idx_fh + 21);

    auto t_xxy_xxxxy = cbuffer.data(idx_fh + 22);

    auto t_xxy_xxxxz = cbuffer.data(idx_fh + 23);

    auto t_xxy_xxxyy = cbuffer.data(idx_fh + 24);

    auto t_xxy_xxxyz = cbuffer.data(idx_fh + 25);

    auto t_xxy_xxxzz = cbuffer.data(idx_fh + 26);

    auto t_xxy_xxyyy = cbuffer.data(idx_fh + 27);

    auto t_xxy_xxyyz = cbuffer.data(idx_fh + 28);

    auto t_xxy_xxyzz = cbuffer.data(idx_fh + 29);

    auto t_xxy_xxzzz = cbuffer.data(idx_fh + 30);

    auto t_xxy_xyyyy = cbuffer.data(idx_fh + 31);

    auto t_xxy_xyyyz = cbuffer.data(idx_fh + 32);

    auto t_xxy_xyyzz = cbuffer.data(idx_fh + 33);

    auto t_xxy_xyzzz = cbuffer.data(idx_fh + 34);

    auto t_xxy_xzzzz = cbuffer.data(idx_fh + 35);

    auto t_xxz_xxxxx = cbuffer.data(idx_fh + 42);

    auto t_xxz_xxxxy = cbuffer.data(idx_fh + 43);

    auto t_xxz_xxxxz = cbuffer.data(idx_fh + 44);

    auto t_xxz_xxxyy = cbuffer.data(idx_fh + 45);

    auto t_xxz_xxxyz = cbuffer.data(idx_fh + 46);

    auto t_xxz_xxxzz = cbuffer.data(idx_fh + 47);

    auto t_xxz_xxyyy = cbuffer.data(idx_fh + 48);

    auto t_xxz_xxyyz = cbuffer.data(idx_fh + 49);

    auto t_xxz_xxyzz = cbuffer.data(idx_fh + 50);

    auto t_xxz_xxzzz = cbuffer.data(idx_fh + 51);

    auto t_xxz_xyyyy = cbuffer.data(idx_fh + 52);

    auto t_xxz_xyyyz = cbuffer.data(idx_fh + 53);

    auto t_xxz_xyyzz = cbuffer.data(idx_fh + 54);

    auto t_xxz_xyzzz = cbuffer.data(idx_fh + 55);

    auto t_xxz_xzzzz = cbuffer.data(idx_fh + 56);

    auto t_xyy_xxxxx = cbuffer.data(idx_fh + 63);

    auto t_xyy_xxxxy = cbuffer.data(idx_fh + 64);

    auto t_xyy_xxxxz = cbuffer.data(idx_fh + 65);

    auto t_xyy_xxxyy = cbuffer.data(idx_fh + 66);

    auto t_xyy_xxxyz = cbuffer.data(idx_fh + 67);

    auto t_xyy_xxxzz = cbuffer.data(idx_fh + 68);

    auto t_xyy_xxyyy = cbuffer.data(idx_fh + 69);

    auto t_xyy_xxyyz = cbuffer.data(idx_fh + 70);

    auto t_xyy_xxyzz = cbuffer.data(idx_fh + 71);

    auto t_xyy_xxzzz = cbuffer.data(idx_fh + 72);

    auto t_xyy_xyyyy = cbuffer.data(idx_fh + 73);

    auto t_xyy_xyyyz = cbuffer.data(idx_fh + 74);

    auto t_xyy_xyyzz = cbuffer.data(idx_fh + 75);

    auto t_xyy_xyzzz = cbuffer.data(idx_fh + 76);

    auto t_xyy_xzzzz = cbuffer.data(idx_fh + 77);

    auto t_xyz_xxxxx = cbuffer.data(idx_fh + 84);

    auto t_xyz_xxxxy = cbuffer.data(idx_fh + 85);

    auto t_xyz_xxxxz = cbuffer.data(idx_fh + 86);

    auto t_xyz_xxxyy = cbuffer.data(idx_fh + 87);

    auto t_xyz_xxxyz = cbuffer.data(idx_fh + 88);

    auto t_xyz_xxxzz = cbuffer.data(idx_fh + 89);

    auto t_xyz_xxyyy = cbuffer.data(idx_fh + 90);

    auto t_xyz_xxyyz = cbuffer.data(idx_fh + 91);

    auto t_xyz_xxyzz = cbuffer.data(idx_fh + 92);

    auto t_xyz_xxzzz = cbuffer.data(idx_fh + 93);

    auto t_xyz_xyyyy = cbuffer.data(idx_fh + 94);

    auto t_xyz_xyyyz = cbuffer.data(idx_fh + 95);

    auto t_xyz_xyyzz = cbuffer.data(idx_fh + 96);

    auto t_xyz_xyzzz = cbuffer.data(idx_fh + 97);

    auto t_xyz_xzzzz = cbuffer.data(idx_fh + 98);

    auto t_xzz_xxxxx = cbuffer.data(idx_fh + 105);

    auto t_xzz_xxxxy = cbuffer.data(idx_fh + 106);

    auto t_xzz_xxxxz = cbuffer.data(idx_fh + 107);

    auto t_xzz_xxxyy = cbuffer.data(idx_fh + 108);

    auto t_xzz_xxxyz = cbuffer.data(idx_fh + 109);

    auto t_xzz_xxxzz = cbuffer.data(idx_fh + 110);

    auto t_xzz_xxyyy = cbuffer.data(idx_fh + 111);

    auto t_xzz_xxyyz = cbuffer.data(idx_fh + 112);

    auto t_xzz_xxyzz = cbuffer.data(idx_fh + 113);

    auto t_xzz_xxzzz = cbuffer.data(idx_fh + 114);

    auto t_xzz_xyyyy = cbuffer.data(idx_fh + 115);

    auto t_xzz_xyyyz = cbuffer.data(idx_fh + 116);

    auto t_xzz_xyyzz = cbuffer.data(idx_fh + 117);

    auto t_xzz_xyzzz = cbuffer.data(idx_fh + 118);

    auto t_xzz_xzzzz = cbuffer.data(idx_fh + 119);

    auto t_yyy_xxxxx = cbuffer.data(idx_fh + 126);

    auto t_yyy_xxxxy = cbuffer.data(idx_fh + 127);

    auto t_yyy_xxxxz = cbuffer.data(idx_fh + 128);

    auto t_yyy_xxxyy = cbuffer.data(idx_fh + 129);

    auto t_yyy_xxxyz = cbuffer.data(idx_fh + 130);

    auto t_yyy_xxxzz = cbuffer.data(idx_fh + 131);

    auto t_yyy_xxyyy = cbuffer.data(idx_fh + 132);

    auto t_yyy_xxyyz = cbuffer.data(idx_fh + 133);

    auto t_yyy_xxyzz = cbuffer.data(idx_fh + 134);

    auto t_yyy_xxzzz = cbuffer.data(idx_fh + 135);

    auto t_yyy_xyyyy = cbuffer.data(idx_fh + 136);

    auto t_yyy_xyyyz = cbuffer.data(idx_fh + 137);

    auto t_yyy_xyyzz = cbuffer.data(idx_fh + 138);

    auto t_yyy_xyzzz = cbuffer.data(idx_fh + 139);

    auto t_yyy_xzzzz = cbuffer.data(idx_fh + 140);

    auto t_yyy_yyyyy = cbuffer.data(idx_fh + 141);

    auto t_yyy_yyyyz = cbuffer.data(idx_fh + 142);

    auto t_yyy_yyyzz = cbuffer.data(idx_fh + 143);

    auto t_yyy_yyzzz = cbuffer.data(idx_fh + 144);

    auto t_yyy_yzzzz = cbuffer.data(idx_fh + 145);

    auto t_yyz_xxxxx = cbuffer.data(idx_fh + 147);

    auto t_yyz_xxxxy = cbuffer.data(idx_fh + 148);

    auto t_yyz_xxxxz = cbuffer.data(idx_fh + 149);

    auto t_yyz_xxxyy = cbuffer.data(idx_fh + 150);

    auto t_yyz_xxxyz = cbuffer.data(idx_fh + 151);

    auto t_yyz_xxxzz = cbuffer.data(idx_fh + 152);

    auto t_yyz_xxyyy = cbuffer.data(idx_fh + 153);

    auto t_yyz_xxyyz = cbuffer.data(idx_fh + 154);

    auto t_yyz_xxyzz = cbuffer.data(idx_fh + 155);

    auto t_yyz_xxzzz = cbuffer.data(idx_fh + 156);

    auto t_yyz_xyyyy = cbuffer.data(idx_fh + 157);

    auto t_yyz_xyyyz = cbuffer.data(idx_fh + 158);

    auto t_yyz_xyyzz = cbuffer.data(idx_fh + 159);

    auto t_yyz_xyzzz = cbuffer.data(idx_fh + 160);

    auto t_yyz_xzzzz = cbuffer.data(idx_fh + 161);

    auto t_yyz_yyyyy = cbuffer.data(idx_fh + 162);

    auto t_yyz_yyyyz = cbuffer.data(idx_fh + 163);

    auto t_yyz_yyyzz = cbuffer.data(idx_fh + 164);

    auto t_yyz_yyzzz = cbuffer.data(idx_fh + 165);

    auto t_yyz_yzzzz = cbuffer.data(idx_fh + 166);

    auto t_yzz_xxxxx = cbuffer.data(idx_fh + 168);

    auto t_yzz_xxxxy = cbuffer.data(idx_fh + 169);

    auto t_yzz_xxxxz = cbuffer.data(idx_fh + 170);

    auto t_yzz_xxxyy = cbuffer.data(idx_fh + 171);

    auto t_yzz_xxxyz = cbuffer.data(idx_fh + 172);

    auto t_yzz_xxxzz = cbuffer.data(idx_fh + 173);

    auto t_yzz_xxyyy = cbuffer.data(idx_fh + 174);

    auto t_yzz_xxyyz = cbuffer.data(idx_fh + 175);

    auto t_yzz_xxyzz = cbuffer.data(idx_fh + 176);

    auto t_yzz_xxzzz = cbuffer.data(idx_fh + 177);

    auto t_yzz_xyyyy = cbuffer.data(idx_fh + 178);

    auto t_yzz_xyyyz = cbuffer.data(idx_fh + 179);

    auto t_yzz_xyyzz = cbuffer.data(idx_fh + 180);

    auto t_yzz_xyzzz = cbuffer.data(idx_fh + 181);

    auto t_yzz_xzzzz = cbuffer.data(idx_fh + 182);

    auto t_yzz_yyyyy = cbuffer.data(idx_fh + 183);

    auto t_yzz_yyyyz = cbuffer.data(idx_fh + 184);

    auto t_yzz_yyyzz = cbuffer.data(idx_fh + 185);

    auto t_yzz_yyzzz = cbuffer.data(idx_fh + 186);

    auto t_yzz_yzzzz = cbuffer.data(idx_fh + 187);

    auto t_zzz_xxxxx = cbuffer.data(idx_fh + 189);

    auto t_zzz_xxxxy = cbuffer.data(idx_fh + 190);

    auto t_zzz_xxxxz = cbuffer.data(idx_fh + 191);

    auto t_zzz_xxxyy = cbuffer.data(idx_fh + 192);

    auto t_zzz_xxxyz = cbuffer.data(idx_fh + 193);

    auto t_zzz_xxxzz = cbuffer.data(idx_fh + 194);

    auto t_zzz_xxyyy = cbuffer.data(idx_fh + 195);

    auto t_zzz_xxyyz = cbuffer.data(idx_fh + 196);

    auto t_zzz_xxyzz = cbuffer.data(idx_fh + 197);

    auto t_zzz_xxzzz = cbuffer.data(idx_fh + 198);

    auto t_zzz_xyyyy = cbuffer.data(idx_fh + 199);

    auto t_zzz_xyyyz = cbuffer.data(idx_fh + 200);

    auto t_zzz_xyyzz = cbuffer.data(idx_fh + 201);

    auto t_zzz_xyzzz = cbuffer.data(idx_fh + 202);

    auto t_zzz_xzzzz = cbuffer.data(idx_fh + 203);

    auto t_zzz_yyyyy = cbuffer.data(idx_fh + 204);

    auto t_zzz_yyyyz = cbuffer.data(idx_fh + 205);

    auto t_zzz_yyyzz = cbuffer.data(idx_fh + 206);

    auto t_zzz_yyzzz = cbuffer.data(idx_fh + 207);

    auto t_zzz_yzzzz = cbuffer.data(idx_fh + 208);

    auto t_zzz_zzzzz = cbuffer.data(idx_fh + 209);

    // Set up components of targeted buffer : GG

    auto t_xxxx_xxxx = cbuffer.data(idx_gg);

    auto t_xxxx_xxxy = cbuffer.data(idx_gg + 1);

    auto t_xxxx_xxxz = cbuffer.data(idx_gg + 2);

    auto t_xxxx_xxyy = cbuffer.data(idx_gg + 3);

    auto t_xxxx_xxyz = cbuffer.data(idx_gg + 4);

    auto t_xxxx_xxzz = cbuffer.data(idx_gg + 5);

    auto t_xxxx_xyyy = cbuffer.data(idx_gg + 6);

    auto t_xxxx_xyyz = cbuffer.data(idx_gg + 7);

    auto t_xxxx_xyzz = cbuffer.data(idx_gg + 8);

    auto t_xxxx_xzzz = cbuffer.data(idx_gg + 9);

    auto t_xxxx_yyyy = cbuffer.data(idx_gg + 10);

    auto t_xxxx_yyyz = cbuffer.data(idx_gg + 11);

    auto t_xxxx_yyzz = cbuffer.data(idx_gg + 12);

    auto t_xxxx_yzzz = cbuffer.data(idx_gg + 13);

    auto t_xxxx_zzzz = cbuffer.data(idx_gg + 14);

    auto t_xxxy_xxxx = cbuffer.data(idx_gg + 15);

    auto t_xxxy_xxxy = cbuffer.data(idx_gg + 16);

    auto t_xxxy_xxxz = cbuffer.data(idx_gg + 17);

    auto t_xxxy_xxyy = cbuffer.data(idx_gg + 18);

    auto t_xxxy_xxyz = cbuffer.data(idx_gg + 19);

    auto t_xxxy_xxzz = cbuffer.data(idx_gg + 20);

    auto t_xxxy_xyyy = cbuffer.data(idx_gg + 21);

    auto t_xxxy_xyyz = cbuffer.data(idx_gg + 22);

    auto t_xxxy_xyzz = cbuffer.data(idx_gg + 23);

    auto t_xxxy_xzzz = cbuffer.data(idx_gg + 24);

    auto t_xxxy_yyyy = cbuffer.data(idx_gg + 25);

    auto t_xxxy_yyyz = cbuffer.data(idx_gg + 26);

    auto t_xxxy_yyzz = cbuffer.data(idx_gg + 27);

    auto t_xxxy_yzzz = cbuffer.data(idx_gg + 28);

    auto t_xxxy_zzzz = cbuffer.data(idx_gg + 29);

    auto t_xxxz_xxxx = cbuffer.data(idx_gg + 30);

    auto t_xxxz_xxxy = cbuffer.data(idx_gg + 31);

    auto t_xxxz_xxxz = cbuffer.data(idx_gg + 32);

    auto t_xxxz_xxyy = cbuffer.data(idx_gg + 33);

    auto t_xxxz_xxyz = cbuffer.data(idx_gg + 34);

    auto t_xxxz_xxzz = cbuffer.data(idx_gg + 35);

    auto t_xxxz_xyyy = cbuffer.data(idx_gg + 36);

    auto t_xxxz_xyyz = cbuffer.data(idx_gg + 37);

    auto t_xxxz_xyzz = cbuffer.data(idx_gg + 38);

    auto t_xxxz_xzzz = cbuffer.data(idx_gg + 39);

    auto t_xxxz_yyyy = cbuffer.data(idx_gg + 40);

    auto t_xxxz_yyyz = cbuffer.data(idx_gg + 41);

    auto t_xxxz_yyzz = cbuffer.data(idx_gg + 42);

    auto t_xxxz_yzzz = cbuffer.data(idx_gg + 43);

    auto t_xxxz_zzzz = cbuffer.data(idx_gg + 44);

    auto t_xxyy_xxxx = cbuffer.data(idx_gg + 45);

    auto t_xxyy_xxxy = cbuffer.data(idx_gg + 46);

    auto t_xxyy_xxxz = cbuffer.data(idx_gg + 47);

    auto t_xxyy_xxyy = cbuffer.data(idx_gg + 48);

    auto t_xxyy_xxyz = cbuffer.data(idx_gg + 49);

    auto t_xxyy_xxzz = cbuffer.data(idx_gg + 50);

    auto t_xxyy_xyyy = cbuffer.data(idx_gg + 51);

    auto t_xxyy_xyyz = cbuffer.data(idx_gg + 52);

    auto t_xxyy_xyzz = cbuffer.data(idx_gg + 53);

    auto t_xxyy_xzzz = cbuffer.data(idx_gg + 54);

    auto t_xxyy_yyyy = cbuffer.data(idx_gg + 55);

    auto t_xxyy_yyyz = cbuffer.data(idx_gg + 56);

    auto t_xxyy_yyzz = cbuffer.data(idx_gg + 57);

    auto t_xxyy_yzzz = cbuffer.data(idx_gg + 58);

    auto t_xxyy_zzzz = cbuffer.data(idx_gg + 59);

    auto t_xxyz_xxxx = cbuffer.data(idx_gg + 60);

    auto t_xxyz_xxxy = cbuffer.data(idx_gg + 61);

    auto t_xxyz_xxxz = cbuffer.data(idx_gg + 62);

    auto t_xxyz_xxyy = cbuffer.data(idx_gg + 63);

    auto t_xxyz_xxyz = cbuffer.data(idx_gg + 64);

    auto t_xxyz_xxzz = cbuffer.data(idx_gg + 65);

    auto t_xxyz_xyyy = cbuffer.data(idx_gg + 66);

    auto t_xxyz_xyyz = cbuffer.data(idx_gg + 67);

    auto t_xxyz_xyzz = cbuffer.data(idx_gg + 68);

    auto t_xxyz_xzzz = cbuffer.data(idx_gg + 69);

    auto t_xxyz_yyyy = cbuffer.data(idx_gg + 70);

    auto t_xxyz_yyyz = cbuffer.data(idx_gg + 71);

    auto t_xxyz_yyzz = cbuffer.data(idx_gg + 72);

    auto t_xxyz_yzzz = cbuffer.data(idx_gg + 73);

    auto t_xxyz_zzzz = cbuffer.data(idx_gg + 74);

    auto t_xxzz_xxxx = cbuffer.data(idx_gg + 75);

    auto t_xxzz_xxxy = cbuffer.data(idx_gg + 76);

    auto t_xxzz_xxxz = cbuffer.data(idx_gg + 77);

    auto t_xxzz_xxyy = cbuffer.data(idx_gg + 78);

    auto t_xxzz_xxyz = cbuffer.data(idx_gg + 79);

    auto t_xxzz_xxzz = cbuffer.data(idx_gg + 80);

    auto t_xxzz_xyyy = cbuffer.data(idx_gg + 81);

    auto t_xxzz_xyyz = cbuffer.data(idx_gg + 82);

    auto t_xxzz_xyzz = cbuffer.data(idx_gg + 83);

    auto t_xxzz_xzzz = cbuffer.data(idx_gg + 84);

    auto t_xxzz_yyyy = cbuffer.data(idx_gg + 85);

    auto t_xxzz_yyyz = cbuffer.data(idx_gg + 86);

    auto t_xxzz_yyzz = cbuffer.data(idx_gg + 87);

    auto t_xxzz_yzzz = cbuffer.data(idx_gg + 88);

    auto t_xxzz_zzzz = cbuffer.data(idx_gg + 89);

    auto t_xyyy_xxxx = cbuffer.data(idx_gg + 90);

    auto t_xyyy_xxxy = cbuffer.data(idx_gg + 91);

    auto t_xyyy_xxxz = cbuffer.data(idx_gg + 92);

    auto t_xyyy_xxyy = cbuffer.data(idx_gg + 93);

    auto t_xyyy_xxyz = cbuffer.data(idx_gg + 94);

    auto t_xyyy_xxzz = cbuffer.data(idx_gg + 95);

    auto t_xyyy_xyyy = cbuffer.data(idx_gg + 96);

    auto t_xyyy_xyyz = cbuffer.data(idx_gg + 97);

    auto t_xyyy_xyzz = cbuffer.data(idx_gg + 98);

    auto t_xyyy_xzzz = cbuffer.data(idx_gg + 99);

    auto t_xyyy_yyyy = cbuffer.data(idx_gg + 100);

    auto t_xyyy_yyyz = cbuffer.data(idx_gg + 101);

    auto t_xyyy_yyzz = cbuffer.data(idx_gg + 102);

    auto t_xyyy_yzzz = cbuffer.data(idx_gg + 103);

    auto t_xyyy_zzzz = cbuffer.data(idx_gg + 104);

    auto t_xyyz_xxxx = cbuffer.data(idx_gg + 105);

    auto t_xyyz_xxxy = cbuffer.data(idx_gg + 106);

    auto t_xyyz_xxxz = cbuffer.data(idx_gg + 107);

    auto t_xyyz_xxyy = cbuffer.data(idx_gg + 108);

    auto t_xyyz_xxyz = cbuffer.data(idx_gg + 109);

    auto t_xyyz_xxzz = cbuffer.data(idx_gg + 110);

    auto t_xyyz_xyyy = cbuffer.data(idx_gg + 111);

    auto t_xyyz_xyyz = cbuffer.data(idx_gg + 112);

    auto t_xyyz_xyzz = cbuffer.data(idx_gg + 113);

    auto t_xyyz_xzzz = cbuffer.data(idx_gg + 114);

    auto t_xyyz_yyyy = cbuffer.data(idx_gg + 115);

    auto t_xyyz_yyyz = cbuffer.data(idx_gg + 116);

    auto t_xyyz_yyzz = cbuffer.data(idx_gg + 117);

    auto t_xyyz_yzzz = cbuffer.data(idx_gg + 118);

    auto t_xyyz_zzzz = cbuffer.data(idx_gg + 119);

    auto t_xyzz_xxxx = cbuffer.data(idx_gg + 120);

    auto t_xyzz_xxxy = cbuffer.data(idx_gg + 121);

    auto t_xyzz_xxxz = cbuffer.data(idx_gg + 122);

    auto t_xyzz_xxyy = cbuffer.data(idx_gg + 123);

    auto t_xyzz_xxyz = cbuffer.data(idx_gg + 124);

    auto t_xyzz_xxzz = cbuffer.data(idx_gg + 125);

    auto t_xyzz_xyyy = cbuffer.data(idx_gg + 126);

    auto t_xyzz_xyyz = cbuffer.data(idx_gg + 127);

    auto t_xyzz_xyzz = cbuffer.data(idx_gg + 128);

    auto t_xyzz_xzzz = cbuffer.data(idx_gg + 129);

    auto t_xyzz_yyyy = cbuffer.data(idx_gg + 130);

    auto t_xyzz_yyyz = cbuffer.data(idx_gg + 131);

    auto t_xyzz_yyzz = cbuffer.data(idx_gg + 132);

    auto t_xyzz_yzzz = cbuffer.data(idx_gg + 133);

    auto t_xyzz_zzzz = cbuffer.data(idx_gg + 134);

    auto t_xzzz_xxxx = cbuffer.data(idx_gg + 135);

    auto t_xzzz_xxxy = cbuffer.data(idx_gg + 136);

    auto t_xzzz_xxxz = cbuffer.data(idx_gg + 137);

    auto t_xzzz_xxyy = cbuffer.data(idx_gg + 138);

    auto t_xzzz_xxyz = cbuffer.data(idx_gg + 139);

    auto t_xzzz_xxzz = cbuffer.data(idx_gg + 140);

    auto t_xzzz_xyyy = cbuffer.data(idx_gg + 141);

    auto t_xzzz_xyyz = cbuffer.data(idx_gg + 142);

    auto t_xzzz_xyzz = cbuffer.data(idx_gg + 143);

    auto t_xzzz_xzzz = cbuffer.data(idx_gg + 144);

    auto t_xzzz_yyyy = cbuffer.data(idx_gg + 145);

    auto t_xzzz_yyyz = cbuffer.data(idx_gg + 146);

    auto t_xzzz_yyzz = cbuffer.data(idx_gg + 147);

    auto t_xzzz_yzzz = cbuffer.data(idx_gg + 148);

    auto t_xzzz_zzzz = cbuffer.data(idx_gg + 149);

    auto t_yyyy_xxxx = cbuffer.data(idx_gg + 150);

    auto t_yyyy_xxxy = cbuffer.data(idx_gg + 151);

    auto t_yyyy_xxxz = cbuffer.data(idx_gg + 152);

    auto t_yyyy_xxyy = cbuffer.data(idx_gg + 153);

    auto t_yyyy_xxyz = cbuffer.data(idx_gg + 154);

    auto t_yyyy_xxzz = cbuffer.data(idx_gg + 155);

    auto t_yyyy_xyyy = cbuffer.data(idx_gg + 156);

    auto t_yyyy_xyyz = cbuffer.data(idx_gg + 157);

    auto t_yyyy_xyzz = cbuffer.data(idx_gg + 158);

    auto t_yyyy_xzzz = cbuffer.data(idx_gg + 159);

    auto t_yyyy_yyyy = cbuffer.data(idx_gg + 160);

    auto t_yyyy_yyyz = cbuffer.data(idx_gg + 161);

    auto t_yyyy_yyzz = cbuffer.data(idx_gg + 162);

    auto t_yyyy_yzzz = cbuffer.data(idx_gg + 163);

    auto t_yyyy_zzzz = cbuffer.data(idx_gg + 164);

    auto t_yyyz_xxxx = cbuffer.data(idx_gg + 165);

    auto t_yyyz_xxxy = cbuffer.data(idx_gg + 166);

    auto t_yyyz_xxxz = cbuffer.data(idx_gg + 167);

    auto t_yyyz_xxyy = cbuffer.data(idx_gg + 168);

    auto t_yyyz_xxyz = cbuffer.data(idx_gg + 169);

    auto t_yyyz_xxzz = cbuffer.data(idx_gg + 170);

    auto t_yyyz_xyyy = cbuffer.data(idx_gg + 171);

    auto t_yyyz_xyyz = cbuffer.data(idx_gg + 172);

    auto t_yyyz_xyzz = cbuffer.data(idx_gg + 173);

    auto t_yyyz_xzzz = cbuffer.data(idx_gg + 174);

    auto t_yyyz_yyyy = cbuffer.data(idx_gg + 175);

    auto t_yyyz_yyyz = cbuffer.data(idx_gg + 176);

    auto t_yyyz_yyzz = cbuffer.data(idx_gg + 177);

    auto t_yyyz_yzzz = cbuffer.data(idx_gg + 178);

    auto t_yyyz_zzzz = cbuffer.data(idx_gg + 179);

    auto t_yyzz_xxxx = cbuffer.data(idx_gg + 180);

    auto t_yyzz_xxxy = cbuffer.data(idx_gg + 181);

    auto t_yyzz_xxxz = cbuffer.data(idx_gg + 182);

    auto t_yyzz_xxyy = cbuffer.data(idx_gg + 183);

    auto t_yyzz_xxyz = cbuffer.data(idx_gg + 184);

    auto t_yyzz_xxzz = cbuffer.data(idx_gg + 185);

    auto t_yyzz_xyyy = cbuffer.data(idx_gg + 186);

    auto t_yyzz_xyyz = cbuffer.data(idx_gg + 187);

    auto t_yyzz_xyzz = cbuffer.data(idx_gg + 188);

    auto t_yyzz_xzzz = cbuffer.data(idx_gg + 189);

    auto t_yyzz_yyyy = cbuffer.data(idx_gg + 190);

    auto t_yyzz_yyyz = cbuffer.data(idx_gg + 191);

    auto t_yyzz_yyzz = cbuffer.data(idx_gg + 192);

    auto t_yyzz_yzzz = cbuffer.data(idx_gg + 193);

    auto t_yyzz_zzzz = cbuffer.data(idx_gg + 194);

    auto t_yzzz_xxxx = cbuffer.data(idx_gg + 195);

    auto t_yzzz_xxxy = cbuffer.data(idx_gg + 196);

    auto t_yzzz_xxxz = cbuffer.data(idx_gg + 197);

    auto t_yzzz_xxyy = cbuffer.data(idx_gg + 198);

    auto t_yzzz_xxyz = cbuffer.data(idx_gg + 199);

    auto t_yzzz_xxzz = cbuffer.data(idx_gg + 200);

    auto t_yzzz_xyyy = cbuffer.data(idx_gg + 201);

    auto t_yzzz_xyyz = cbuffer.data(idx_gg + 202);

    auto t_yzzz_xyzz = cbuffer.data(idx_gg + 203);

    auto t_yzzz_xzzz = cbuffer.data(idx_gg + 204);

    auto t_yzzz_yyyy = cbuffer.data(idx_gg + 205);

    auto t_yzzz_yyyz = cbuffer.data(idx_gg + 206);

    auto t_yzzz_yyzz = cbuffer.data(idx_gg + 207);

    auto t_yzzz_yzzz = cbuffer.data(idx_gg + 208);

    auto t_yzzz_zzzz = cbuffer.data(idx_gg + 209);

    auto t_zzzz_xxxx = cbuffer.data(idx_gg + 210);

    auto t_zzzz_xxxy = cbuffer.data(idx_gg + 211);

    auto t_zzzz_xxxz = cbuffer.data(idx_gg + 212);

    auto t_zzzz_xxyy = cbuffer.data(idx_gg + 213);

    auto t_zzzz_xxyz = cbuffer.data(idx_gg + 214);

    auto t_zzzz_xxzz = cbuffer.data(idx_gg + 215);

    auto t_zzzz_xyyy = cbuffer.data(idx_gg + 216);

    auto t_zzzz_xyyz = cbuffer.data(idx_gg + 217);

    auto t_zzzz_xyzz = cbuffer.data(idx_gg + 218);

    auto t_zzzz_xzzz = cbuffer.data(idx_gg + 219);

    auto t_zzzz_yyyy = cbuffer.data(idx_gg + 220);

    auto t_zzzz_yyyz = cbuffer.data(idx_gg + 221);

    auto t_zzzz_yyzz = cbuffer.data(idx_gg + 222);

    auto t_zzzz_yzzz = cbuffer.data(idx_gg + 223);

    auto t_zzzz_zzzz = cbuffer.data(idx_gg + 224);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxx_xxxx, t_xxx_xxxxx, t_xxx_xxxxy, t_xxx_xxxxz, t_xxx_xxxy, t_xxx_xxxyy, t_xxx_xxxyz, t_xxx_xxxz, t_xxx_xxxzz, t_xxx_xxyy, t_xxx_xxyyy, t_xxx_xxyyz, t_xxx_xxyz, t_xxx_xxyzz, t_xxx_xxzz, t_xxx_xxzzz, t_xxx_xyyy, t_xxx_xyyyy, t_xxx_xyyyz, t_xxx_xyyz, t_xxx_xyyzz, t_xxx_xyzz, t_xxx_xyzzz, t_xxx_xzzz, t_xxx_xzzzz, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, t_xxx_zzzz, t_xxxx_xxxx, t_xxxx_xxxy, t_xxxx_xxxz, t_xxxx_xxyy, t_xxxx_xxyz, t_xxxx_xxzz, t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz, t_xxxx_xzzz, t_xxxx_yyyy, t_xxxx_yyyz, t_xxxx_yyzz, t_xxxx_yzzz, t_xxxx_zzzz, t_xxxy_xxxx, t_xxxy_xxxy, t_xxxy_xxxz, t_xxxy_xxyy, t_xxxy_xxyz, t_xxxy_xxzz, t_xxxy_xyyy, t_xxxy_xyyz, t_xxxy_xyzz, t_xxxy_xzzz, t_xxxy_yyyy, t_xxxy_yyyz, t_xxxy_yyzz, t_xxxy_yzzz, t_xxxy_zzzz, t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz, t_xxxz_xxyy, t_xxxz_xxyz, t_xxxz_xxzz, t_xxxz_xyyy, t_xxxz_xyyz, t_xxxz_xyzz, t_xxxz_xzzz, t_xxxz_yyyy, t_xxxz_yyyz, t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz, t_xxy_xxxx, t_xxy_xxxxx, t_xxy_xxxxy, t_xxy_xxxxz, t_xxy_xxxy, t_xxy_xxxyy, t_xxy_xxxyz, t_xxy_xxxz, t_xxy_xxxzz, t_xxy_xxyy, t_xxy_xxyyy, t_xxy_xxyyz, t_xxy_xxyz, t_xxy_xxyzz, t_xxy_xxzz, t_xxy_xxzzz, t_xxy_xyyy, t_xxy_xyyyy, t_xxy_xyyyz, t_xxy_xyyz, t_xxy_xyyzz, t_xxy_xyzz, t_xxy_xyzzz, t_xxy_xzzz, t_xxy_xzzzz, t_xxy_yyyy, t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz, t_xxyy_xxxx, t_xxyy_xxxy, t_xxyy_xxxz, t_xxyy_xxyy, t_xxyy_xxyz, t_xxyy_xxzz, t_xxyy_xyyy, t_xxyy_xyyz, t_xxyy_xyzz, t_xxyy_xzzz, t_xxyy_yyyy, t_xxyy_yyyz, t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz, t_xxyz_xxxx, t_xxyz_xxxy, t_xxyz_xxxz, t_xxyz_xxyy, t_xxyz_xxyz, t_xxyz_xxzz, t_xxyz_xyyy, t_xxyz_xyyz, t_xxyz_xyzz, t_xxyz_xzzz, t_xxyz_yyyy, t_xxyz_yyyz, t_xxyz_yyzz, t_xxyz_yzzz, t_xxyz_zzzz, t_xxz_xxxx, t_xxz_xxxxx, t_xxz_xxxxy, t_xxz_xxxxz, t_xxz_xxxy, t_xxz_xxxyy, t_xxz_xxxyz, t_xxz_xxxz, t_xxz_xxxzz, t_xxz_xxyy, t_xxz_xxyyy, t_xxz_xxyyz, t_xxz_xxyz, t_xxz_xxyzz, t_xxz_xxzz, t_xxz_xxzzz, t_xxz_xyyy, t_xxz_xyyyy, t_xxz_xyyyz, t_xxz_xyyz, t_xxz_xyyzz, t_xxz_xyzz, t_xxz_xyzzz, t_xxz_xzzz, t_xxz_xzzzz, t_xxz_yyyy, t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz, t_xxzz_xxxx, t_xxzz_xxxy, t_xxzz_xxxz, t_xxzz_xxyy, t_xxzz_xxyz, t_xxzz_xxzz, t_xxzz_xyyy, t_xxzz_xyyz, t_xxzz_xyzz, t_xxzz_xzzz, t_xxzz_yyyy, t_xxzz_yyyz, t_xxzz_yyzz, t_xxzz_yzzz, t_xxzz_zzzz, t_xyy_xxxx, t_xyy_xxxxx, t_xyy_xxxxy, t_xyy_xxxxz, t_xyy_xxxy, t_xyy_xxxyy, t_xyy_xxxyz, t_xyy_xxxz, t_xyy_xxxzz, t_xyy_xxyy, t_xyy_xxyyy, t_xyy_xxyyz, t_xyy_xxyz, t_xyy_xxyzz, t_xyy_xxzz, t_xyy_xxzzz, t_xyy_xyyy, t_xyy_xyyyy, t_xyy_xyyyz, t_xyy_xyyz, t_xyy_xyyzz, t_xyy_xyzz, t_xyy_xyzzz, t_xyy_xzzz, t_xyy_xzzzz, t_xyy_yyyy, t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz, t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz, t_xyyy_xxyy, t_xyyy_xxyz, t_xyyy_xxzz, t_xyyy_xyyy, t_xyyy_xyyz, t_xyyy_xyzz, t_xyyy_xzzz, t_xyyy_yyyy, t_xyyy_yyyz, t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz, t_xyyz_xxxx, t_xyyz_xxxy, t_xyyz_xxxz, t_xyyz_xxyy, t_xyyz_xxyz, t_xyyz_xxzz, t_xyyz_xyyy, t_xyyz_xyyz, t_xyyz_xyzz, t_xyyz_xzzz, t_xyyz_yyyy, t_xyyz_yyyz, t_xyyz_yyzz, t_xyyz_yzzz, t_xyyz_zzzz, t_xyz_xxxx, t_xyz_xxxxx, t_xyz_xxxxy, t_xyz_xxxxz, t_xyz_xxxy, t_xyz_xxxyy, t_xyz_xxxyz, t_xyz_xxxz, t_xyz_xxxzz, t_xyz_xxyy, t_xyz_xxyyy, t_xyz_xxyyz, t_xyz_xxyz, t_xyz_xxyzz, t_xyz_xxzz, t_xyz_xxzzz, t_xyz_xyyy, t_xyz_xyyyy, t_xyz_xyyyz, t_xyz_xyyz, t_xyz_xyyzz, t_xyz_xyzz, t_xyz_xyzzz, t_xyz_xzzz, t_xyz_xzzzz, t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, t_xyz_yzzz, t_xyz_zzzz, t_xyzz_xxxx, t_xyzz_xxxy, t_xyzz_xxxz, t_xyzz_xxyy, t_xyzz_xxyz, t_xyzz_xxzz, t_xyzz_xyyy, t_xyzz_xyyz, t_xyzz_xyzz, t_xyzz_xzzz, t_xyzz_yyyy, t_xyzz_yyyz, t_xyzz_yyzz, t_xyzz_yzzz, t_xyzz_zzzz, t_xzz_xxxx, t_xzz_xxxxx, t_xzz_xxxxy, t_xzz_xxxxz, t_xzz_xxxy, t_xzz_xxxyy, t_xzz_xxxyz, t_xzz_xxxz, t_xzz_xxxzz, t_xzz_xxyy, t_xzz_xxyyy, t_xzz_xxyyz, t_xzz_xxyz, t_xzz_xxyzz, t_xzz_xxzz, t_xzz_xxzzz, t_xzz_xyyy, t_xzz_xyyyy, t_xzz_xyyyz, t_xzz_xyyz, t_xzz_xyyzz, t_xzz_xyzz, t_xzz_xyzzz, t_xzz_xzzz, t_xzz_xzzzz, t_xzz_yyyy, t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, t_xzz_zzzz, t_xzzz_xxxx, t_xzzz_xxxy, t_xzzz_xxxz, t_xzzz_xxyy, t_xzzz_xxyz, t_xzzz_xxzz, t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz, t_xzzz_xzzz, t_xzzz_yyyy, t_xzzz_yyyz, t_xzzz_yyzz, t_xzzz_yzzz, t_xzzz_zzzz, t_yyy_xxxx, t_yyy_xxxxx, t_yyy_xxxxy, t_yyy_xxxxz, t_yyy_xxxy, t_yyy_xxxyy, t_yyy_xxxyz, t_yyy_xxxz, t_yyy_xxxzz, t_yyy_xxyy, t_yyy_xxyyy, t_yyy_xxyyz, t_yyy_xxyz, t_yyy_xxyzz, t_yyy_xxzz, t_yyy_xxzzz, t_yyy_xyyy, t_yyy_xyyyy, t_yyy_xyyyz, t_yyy_xyyz, t_yyy_xyyzz, t_yyy_xyzz, t_yyy_xyzzz, t_yyy_xzzz, t_yyy_xzzzz, t_yyy_yyyy, t_yyy_yyyyy, t_yyy_yyyyz, t_yyy_yyyz, t_yyy_yyyzz, t_yyy_yyzz, t_yyy_yyzzz, t_yyy_yzzz, t_yyy_yzzzz, t_yyy_zzzz, t_yyyy_xxxx, t_yyyy_xxxy, t_yyyy_xxxz, t_yyyy_xxyy, t_yyyy_xxyz, t_yyyy_xxzz, t_yyyy_xyyy, t_yyyy_xyyz, t_yyyy_xyzz, t_yyyy_xzzz, t_yyyy_yyyy, t_yyyy_yyyz, t_yyyy_yyzz, t_yyyy_yzzz, t_yyyy_zzzz, t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz, t_yyyz_xxyy, t_yyyz_xxyz, t_yyyz_xxzz, t_yyyz_xyyy, t_yyyz_xyyz, t_yyyz_xyzz, t_yyyz_xzzz, t_yyyz_yyyy, t_yyyz_yyyz, t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz, t_yyz_xxxx, t_yyz_xxxxx, t_yyz_xxxxy, t_yyz_xxxxz, t_yyz_xxxy, t_yyz_xxxyy, t_yyz_xxxyz, t_yyz_xxxz, t_yyz_xxxzz, t_yyz_xxyy, t_yyz_xxyyy, t_yyz_xxyyz, t_yyz_xxyz, t_yyz_xxyzz, t_yyz_xxzz, t_yyz_xxzzz, t_yyz_xyyy, t_yyz_xyyyy, t_yyz_xyyyz, t_yyz_xyyz, t_yyz_xyyzz, t_yyz_xyzz, t_yyz_xyzzz, t_yyz_xzzz, t_yyz_xzzzz, t_yyz_yyyy, t_yyz_yyyyy, t_yyz_yyyyz, t_yyz_yyyz, t_yyz_yyyzz, t_yyz_yyzz, t_yyz_yyzzz, t_yyz_yzzz, t_yyz_yzzzz, t_yyz_zzzz, t_yyzz_xxxx, t_yyzz_xxxy, t_yyzz_xxxz, t_yyzz_xxyy, t_yyzz_xxyz, t_yyzz_xxzz, t_yyzz_xyyy, t_yyzz_xyyz, t_yyzz_xyzz, t_yyzz_xzzz, t_yyzz_yyyy, t_yyzz_yyyz, t_yyzz_yyzz, t_yyzz_yzzz, t_yyzz_zzzz, t_yzz_xxxx, t_yzz_xxxxx, t_yzz_xxxxy, t_yzz_xxxxz, t_yzz_xxxy, t_yzz_xxxyy, t_yzz_xxxyz, t_yzz_xxxz, t_yzz_xxxzz, t_yzz_xxyy, t_yzz_xxyyy, t_yzz_xxyyz, t_yzz_xxyz, t_yzz_xxyzz, t_yzz_xxzz, t_yzz_xxzzz, t_yzz_xyyy, t_yzz_xyyyy, t_yzz_xyyyz, t_yzz_xyyz, t_yzz_xyyzz, t_yzz_xyzz, t_yzz_xyzzz, t_yzz_xzzz, t_yzz_xzzzz, t_yzz_yyyy, t_yzz_yyyyy, t_yzz_yyyyz, t_yzz_yyyz, t_yzz_yyyzz, t_yzz_yyzz, t_yzz_yyzzz, t_yzz_yzzz, t_yzz_yzzzz, t_yzz_zzzz, t_yzzz_xxxx, t_yzzz_xxxy, t_yzzz_xxxz, t_yzzz_xxyy, t_yzzz_xxyz, t_yzzz_xxzz, t_yzzz_xyyy, t_yzzz_xyyz, t_yzzz_xyzz, t_yzzz_xzzz, t_yzzz_yyyy, t_yzzz_yyyz, t_yzzz_yyzz, t_yzzz_yzzz, t_yzzz_zzzz, t_zzz_xxxx, t_zzz_xxxxx, t_zzz_xxxxy, t_zzz_xxxxz, t_zzz_xxxy, t_zzz_xxxyy, t_zzz_xxxyz, t_zzz_xxxz, t_zzz_xxxzz, t_zzz_xxyy, t_zzz_xxyyy, t_zzz_xxyyz, t_zzz_xxyz, t_zzz_xxyzz, t_zzz_xxzz, t_zzz_xxzzz, t_zzz_xyyy, t_zzz_xyyyy, t_zzz_xyyyz, t_zzz_xyyz, t_zzz_xyyzz, t_zzz_xyzz, t_zzz_xyzzz, t_zzz_xzzz, t_zzz_xzzzz, t_zzz_yyyy, t_zzz_yyyyy, t_zzz_yyyyz, t_zzz_yyyz, t_zzz_yyyzz, t_zzz_yyzz, t_zzz_yyzzz, t_zzz_yzzz, t_zzz_yzzzz, t_zzz_zzzz, t_zzz_zzzzz, t_zzzz_xxxx, t_zzzz_xxxy, t_zzzz_xxxz, t_zzzz_xxyy, t_zzzz_xxyz, t_zzzz_xxzz, t_zzzz_xyyy, t_zzzz_xyyz, t_zzzz_xyzz, t_zzzz_xzzz, t_zzzz_yyyy, t_zzzz_yyyz, t_zzzz_yyzz, t_zzzz_yzzz, t_zzzz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxx_xxxx[i] = -t_xxx_xxxx[i] * ab_x[i] + t_xxx_xxxxx[i];

        t_xxxx_xxxy[i] = -t_xxx_xxxy[i] * ab_x[i] + t_xxx_xxxxy[i];

        t_xxxx_xxxz[i] = -t_xxx_xxxz[i] * ab_x[i] + t_xxx_xxxxz[i];

        t_xxxx_xxyy[i] = -t_xxx_xxyy[i] * ab_x[i] + t_xxx_xxxyy[i];

        t_xxxx_xxyz[i] = -t_xxx_xxyz[i] * ab_x[i] + t_xxx_xxxyz[i];

        t_xxxx_xxzz[i] = -t_xxx_xxzz[i] * ab_x[i] + t_xxx_xxxzz[i];

        t_xxxx_xyyy[i] = -t_xxx_xyyy[i] * ab_x[i] + t_xxx_xxyyy[i];

        t_xxxx_xyyz[i] = -t_xxx_xyyz[i] * ab_x[i] + t_xxx_xxyyz[i];

        t_xxxx_xyzz[i] = -t_xxx_xyzz[i] * ab_x[i] + t_xxx_xxyzz[i];

        t_xxxx_xzzz[i] = -t_xxx_xzzz[i] * ab_x[i] + t_xxx_xxzzz[i];

        t_xxxx_yyyy[i] = -t_xxx_yyyy[i] * ab_x[i] + t_xxx_xyyyy[i];

        t_xxxx_yyyz[i] = -t_xxx_yyyz[i] * ab_x[i] + t_xxx_xyyyz[i];

        t_xxxx_yyzz[i] = -t_xxx_yyzz[i] * ab_x[i] + t_xxx_xyyzz[i];

        t_xxxx_yzzz[i] = -t_xxx_yzzz[i] * ab_x[i] + t_xxx_xyzzz[i];

        t_xxxx_zzzz[i] = -t_xxx_zzzz[i] * ab_x[i] + t_xxx_xzzzz[i];

        t_xxxy_xxxx[i] = -t_xxy_xxxx[i] * ab_x[i] + t_xxy_xxxxx[i];

        t_xxxy_xxxy[i] = -t_xxy_xxxy[i] * ab_x[i] + t_xxy_xxxxy[i];

        t_xxxy_xxxz[i] = -t_xxy_xxxz[i] * ab_x[i] + t_xxy_xxxxz[i];

        t_xxxy_xxyy[i] = -t_xxy_xxyy[i] * ab_x[i] + t_xxy_xxxyy[i];

        t_xxxy_xxyz[i] = -t_xxy_xxyz[i] * ab_x[i] + t_xxy_xxxyz[i];

        t_xxxy_xxzz[i] = -t_xxy_xxzz[i] * ab_x[i] + t_xxy_xxxzz[i];

        t_xxxy_xyyy[i] = -t_xxy_xyyy[i] * ab_x[i] + t_xxy_xxyyy[i];

        t_xxxy_xyyz[i] = -t_xxy_xyyz[i] * ab_x[i] + t_xxy_xxyyz[i];

        t_xxxy_xyzz[i] = -t_xxy_xyzz[i] * ab_x[i] + t_xxy_xxyzz[i];

        t_xxxy_xzzz[i] = -t_xxy_xzzz[i] * ab_x[i] + t_xxy_xxzzz[i];

        t_xxxy_yyyy[i] = -t_xxy_yyyy[i] * ab_x[i] + t_xxy_xyyyy[i];

        t_xxxy_yyyz[i] = -t_xxy_yyyz[i] * ab_x[i] + t_xxy_xyyyz[i];

        t_xxxy_yyzz[i] = -t_xxy_yyzz[i] * ab_x[i] + t_xxy_xyyzz[i];

        t_xxxy_yzzz[i] = -t_xxy_yzzz[i] * ab_x[i] + t_xxy_xyzzz[i];

        t_xxxy_zzzz[i] = -t_xxy_zzzz[i] * ab_x[i] + t_xxy_xzzzz[i];

        t_xxxz_xxxx[i] = -t_xxz_xxxx[i] * ab_x[i] + t_xxz_xxxxx[i];

        t_xxxz_xxxy[i] = -t_xxz_xxxy[i] * ab_x[i] + t_xxz_xxxxy[i];

        t_xxxz_xxxz[i] = -t_xxz_xxxz[i] * ab_x[i] + t_xxz_xxxxz[i];

        t_xxxz_xxyy[i] = -t_xxz_xxyy[i] * ab_x[i] + t_xxz_xxxyy[i];

        t_xxxz_xxyz[i] = -t_xxz_xxyz[i] * ab_x[i] + t_xxz_xxxyz[i];

        t_xxxz_xxzz[i] = -t_xxz_xxzz[i] * ab_x[i] + t_xxz_xxxzz[i];

        t_xxxz_xyyy[i] = -t_xxz_xyyy[i] * ab_x[i] + t_xxz_xxyyy[i];

        t_xxxz_xyyz[i] = -t_xxz_xyyz[i] * ab_x[i] + t_xxz_xxyyz[i];

        t_xxxz_xyzz[i] = -t_xxz_xyzz[i] * ab_x[i] + t_xxz_xxyzz[i];

        t_xxxz_xzzz[i] = -t_xxz_xzzz[i] * ab_x[i] + t_xxz_xxzzz[i];

        t_xxxz_yyyy[i] = -t_xxz_yyyy[i] * ab_x[i] + t_xxz_xyyyy[i];

        t_xxxz_yyyz[i] = -t_xxz_yyyz[i] * ab_x[i] + t_xxz_xyyyz[i];

        t_xxxz_yyzz[i] = -t_xxz_yyzz[i] * ab_x[i] + t_xxz_xyyzz[i];

        t_xxxz_yzzz[i] = -t_xxz_yzzz[i] * ab_x[i] + t_xxz_xyzzz[i];

        t_xxxz_zzzz[i] = -t_xxz_zzzz[i] * ab_x[i] + t_xxz_xzzzz[i];

        t_xxyy_xxxx[i] = -t_xyy_xxxx[i] * ab_x[i] + t_xyy_xxxxx[i];

        t_xxyy_xxxy[i] = -t_xyy_xxxy[i] * ab_x[i] + t_xyy_xxxxy[i];

        t_xxyy_xxxz[i] = -t_xyy_xxxz[i] * ab_x[i] + t_xyy_xxxxz[i];

        t_xxyy_xxyy[i] = -t_xyy_xxyy[i] * ab_x[i] + t_xyy_xxxyy[i];

        t_xxyy_xxyz[i] = -t_xyy_xxyz[i] * ab_x[i] + t_xyy_xxxyz[i];

        t_xxyy_xxzz[i] = -t_xyy_xxzz[i] * ab_x[i] + t_xyy_xxxzz[i];

        t_xxyy_xyyy[i] = -t_xyy_xyyy[i] * ab_x[i] + t_xyy_xxyyy[i];

        t_xxyy_xyyz[i] = -t_xyy_xyyz[i] * ab_x[i] + t_xyy_xxyyz[i];

        t_xxyy_xyzz[i] = -t_xyy_xyzz[i] * ab_x[i] + t_xyy_xxyzz[i];

        t_xxyy_xzzz[i] = -t_xyy_xzzz[i] * ab_x[i] + t_xyy_xxzzz[i];

        t_xxyy_yyyy[i] = -t_xyy_yyyy[i] * ab_x[i] + t_xyy_xyyyy[i];

        t_xxyy_yyyz[i] = -t_xyy_yyyz[i] * ab_x[i] + t_xyy_xyyyz[i];

        t_xxyy_yyzz[i] = -t_xyy_yyzz[i] * ab_x[i] + t_xyy_xyyzz[i];

        t_xxyy_yzzz[i] = -t_xyy_yzzz[i] * ab_x[i] + t_xyy_xyzzz[i];

        t_xxyy_zzzz[i] = -t_xyy_zzzz[i] * ab_x[i] + t_xyy_xzzzz[i];

        t_xxyz_xxxx[i] = -t_xyz_xxxx[i] * ab_x[i] + t_xyz_xxxxx[i];

        t_xxyz_xxxy[i] = -t_xyz_xxxy[i] * ab_x[i] + t_xyz_xxxxy[i];

        t_xxyz_xxxz[i] = -t_xyz_xxxz[i] * ab_x[i] + t_xyz_xxxxz[i];

        t_xxyz_xxyy[i] = -t_xyz_xxyy[i] * ab_x[i] + t_xyz_xxxyy[i];

        t_xxyz_xxyz[i] = -t_xyz_xxyz[i] * ab_x[i] + t_xyz_xxxyz[i];

        t_xxyz_xxzz[i] = -t_xyz_xxzz[i] * ab_x[i] + t_xyz_xxxzz[i];

        t_xxyz_xyyy[i] = -t_xyz_xyyy[i] * ab_x[i] + t_xyz_xxyyy[i];

        t_xxyz_xyyz[i] = -t_xyz_xyyz[i] * ab_x[i] + t_xyz_xxyyz[i];

        t_xxyz_xyzz[i] = -t_xyz_xyzz[i] * ab_x[i] + t_xyz_xxyzz[i];

        t_xxyz_xzzz[i] = -t_xyz_xzzz[i] * ab_x[i] + t_xyz_xxzzz[i];

        t_xxyz_yyyy[i] = -t_xyz_yyyy[i] * ab_x[i] + t_xyz_xyyyy[i];

        t_xxyz_yyyz[i] = -t_xyz_yyyz[i] * ab_x[i] + t_xyz_xyyyz[i];

        t_xxyz_yyzz[i] = -t_xyz_yyzz[i] * ab_x[i] + t_xyz_xyyzz[i];

        t_xxyz_yzzz[i] = -t_xyz_yzzz[i] * ab_x[i] + t_xyz_xyzzz[i];

        t_xxyz_zzzz[i] = -t_xyz_zzzz[i] * ab_x[i] + t_xyz_xzzzz[i];

        t_xxzz_xxxx[i] = -t_xzz_xxxx[i] * ab_x[i] + t_xzz_xxxxx[i];

        t_xxzz_xxxy[i] = -t_xzz_xxxy[i] * ab_x[i] + t_xzz_xxxxy[i];

        t_xxzz_xxxz[i] = -t_xzz_xxxz[i] * ab_x[i] + t_xzz_xxxxz[i];

        t_xxzz_xxyy[i] = -t_xzz_xxyy[i] * ab_x[i] + t_xzz_xxxyy[i];

        t_xxzz_xxyz[i] = -t_xzz_xxyz[i] * ab_x[i] + t_xzz_xxxyz[i];

        t_xxzz_xxzz[i] = -t_xzz_xxzz[i] * ab_x[i] + t_xzz_xxxzz[i];

        t_xxzz_xyyy[i] = -t_xzz_xyyy[i] * ab_x[i] + t_xzz_xxyyy[i];

        t_xxzz_xyyz[i] = -t_xzz_xyyz[i] * ab_x[i] + t_xzz_xxyyz[i];

        t_xxzz_xyzz[i] = -t_xzz_xyzz[i] * ab_x[i] + t_xzz_xxyzz[i];

        t_xxzz_xzzz[i] = -t_xzz_xzzz[i] * ab_x[i] + t_xzz_xxzzz[i];

        t_xxzz_yyyy[i] = -t_xzz_yyyy[i] * ab_x[i] + t_xzz_xyyyy[i];

        t_xxzz_yyyz[i] = -t_xzz_yyyz[i] * ab_x[i] + t_xzz_xyyyz[i];

        t_xxzz_yyzz[i] = -t_xzz_yyzz[i] * ab_x[i] + t_xzz_xyyzz[i];

        t_xxzz_yzzz[i] = -t_xzz_yzzz[i] * ab_x[i] + t_xzz_xyzzz[i];

        t_xxzz_zzzz[i] = -t_xzz_zzzz[i] * ab_x[i] + t_xzz_xzzzz[i];

        t_xyyy_xxxx[i] = -t_yyy_xxxx[i] * ab_x[i] + t_yyy_xxxxx[i];

        t_xyyy_xxxy[i] = -t_yyy_xxxy[i] * ab_x[i] + t_yyy_xxxxy[i];

        t_xyyy_xxxz[i] = -t_yyy_xxxz[i] * ab_x[i] + t_yyy_xxxxz[i];

        t_xyyy_xxyy[i] = -t_yyy_xxyy[i] * ab_x[i] + t_yyy_xxxyy[i];

        t_xyyy_xxyz[i] = -t_yyy_xxyz[i] * ab_x[i] + t_yyy_xxxyz[i];

        t_xyyy_xxzz[i] = -t_yyy_xxzz[i] * ab_x[i] + t_yyy_xxxzz[i];

        t_xyyy_xyyy[i] = -t_yyy_xyyy[i] * ab_x[i] + t_yyy_xxyyy[i];

        t_xyyy_xyyz[i] = -t_yyy_xyyz[i] * ab_x[i] + t_yyy_xxyyz[i];

        t_xyyy_xyzz[i] = -t_yyy_xyzz[i] * ab_x[i] + t_yyy_xxyzz[i];

        t_xyyy_xzzz[i] = -t_yyy_xzzz[i] * ab_x[i] + t_yyy_xxzzz[i];

        t_xyyy_yyyy[i] = -t_yyy_yyyy[i] * ab_x[i] + t_yyy_xyyyy[i];

        t_xyyy_yyyz[i] = -t_yyy_yyyz[i] * ab_x[i] + t_yyy_xyyyz[i];

        t_xyyy_yyzz[i] = -t_yyy_yyzz[i] * ab_x[i] + t_yyy_xyyzz[i];

        t_xyyy_yzzz[i] = -t_yyy_yzzz[i] * ab_x[i] + t_yyy_xyzzz[i];

        t_xyyy_zzzz[i] = -t_yyy_zzzz[i] * ab_x[i] + t_yyy_xzzzz[i];

        t_xyyz_xxxx[i] = -t_yyz_xxxx[i] * ab_x[i] + t_yyz_xxxxx[i];

        t_xyyz_xxxy[i] = -t_yyz_xxxy[i] * ab_x[i] + t_yyz_xxxxy[i];

        t_xyyz_xxxz[i] = -t_yyz_xxxz[i] * ab_x[i] + t_yyz_xxxxz[i];

        t_xyyz_xxyy[i] = -t_yyz_xxyy[i] * ab_x[i] + t_yyz_xxxyy[i];

        t_xyyz_xxyz[i] = -t_yyz_xxyz[i] * ab_x[i] + t_yyz_xxxyz[i];

        t_xyyz_xxzz[i] = -t_yyz_xxzz[i] * ab_x[i] + t_yyz_xxxzz[i];

        t_xyyz_xyyy[i] = -t_yyz_xyyy[i] * ab_x[i] + t_yyz_xxyyy[i];

        t_xyyz_xyyz[i] = -t_yyz_xyyz[i] * ab_x[i] + t_yyz_xxyyz[i];

        t_xyyz_xyzz[i] = -t_yyz_xyzz[i] * ab_x[i] + t_yyz_xxyzz[i];

        t_xyyz_xzzz[i] = -t_yyz_xzzz[i] * ab_x[i] + t_yyz_xxzzz[i];

        t_xyyz_yyyy[i] = -t_yyz_yyyy[i] * ab_x[i] + t_yyz_xyyyy[i];

        t_xyyz_yyyz[i] = -t_yyz_yyyz[i] * ab_x[i] + t_yyz_xyyyz[i];

        t_xyyz_yyzz[i] = -t_yyz_yyzz[i] * ab_x[i] + t_yyz_xyyzz[i];

        t_xyyz_yzzz[i] = -t_yyz_yzzz[i] * ab_x[i] + t_yyz_xyzzz[i];

        t_xyyz_zzzz[i] = -t_yyz_zzzz[i] * ab_x[i] + t_yyz_xzzzz[i];

        t_xyzz_xxxx[i] = -t_yzz_xxxx[i] * ab_x[i] + t_yzz_xxxxx[i];

        t_xyzz_xxxy[i] = -t_yzz_xxxy[i] * ab_x[i] + t_yzz_xxxxy[i];

        t_xyzz_xxxz[i] = -t_yzz_xxxz[i] * ab_x[i] + t_yzz_xxxxz[i];

        t_xyzz_xxyy[i] = -t_yzz_xxyy[i] * ab_x[i] + t_yzz_xxxyy[i];

        t_xyzz_xxyz[i] = -t_yzz_xxyz[i] * ab_x[i] + t_yzz_xxxyz[i];

        t_xyzz_xxzz[i] = -t_yzz_xxzz[i] * ab_x[i] + t_yzz_xxxzz[i];

        t_xyzz_xyyy[i] = -t_yzz_xyyy[i] * ab_x[i] + t_yzz_xxyyy[i];

        t_xyzz_xyyz[i] = -t_yzz_xyyz[i] * ab_x[i] + t_yzz_xxyyz[i];

        t_xyzz_xyzz[i] = -t_yzz_xyzz[i] * ab_x[i] + t_yzz_xxyzz[i];

        t_xyzz_xzzz[i] = -t_yzz_xzzz[i] * ab_x[i] + t_yzz_xxzzz[i];

        t_xyzz_yyyy[i] = -t_yzz_yyyy[i] * ab_x[i] + t_yzz_xyyyy[i];

        t_xyzz_yyyz[i] = -t_yzz_yyyz[i] * ab_x[i] + t_yzz_xyyyz[i];

        t_xyzz_yyzz[i] = -t_yzz_yyzz[i] * ab_x[i] + t_yzz_xyyzz[i];

        t_xyzz_yzzz[i] = -t_yzz_yzzz[i] * ab_x[i] + t_yzz_xyzzz[i];

        t_xyzz_zzzz[i] = -t_yzz_zzzz[i] * ab_x[i] + t_yzz_xzzzz[i];

        t_xzzz_xxxx[i] = -t_zzz_xxxx[i] * ab_x[i] + t_zzz_xxxxx[i];

        t_xzzz_xxxy[i] = -t_zzz_xxxy[i] * ab_x[i] + t_zzz_xxxxy[i];

        t_xzzz_xxxz[i] = -t_zzz_xxxz[i] * ab_x[i] + t_zzz_xxxxz[i];

        t_xzzz_xxyy[i] = -t_zzz_xxyy[i] * ab_x[i] + t_zzz_xxxyy[i];

        t_xzzz_xxyz[i] = -t_zzz_xxyz[i] * ab_x[i] + t_zzz_xxxyz[i];

        t_xzzz_xxzz[i] = -t_zzz_xxzz[i] * ab_x[i] + t_zzz_xxxzz[i];

        t_xzzz_xyyy[i] = -t_zzz_xyyy[i] * ab_x[i] + t_zzz_xxyyy[i];

        t_xzzz_xyyz[i] = -t_zzz_xyyz[i] * ab_x[i] + t_zzz_xxyyz[i];

        t_xzzz_xyzz[i] = -t_zzz_xyzz[i] * ab_x[i] + t_zzz_xxyzz[i];

        t_xzzz_xzzz[i] = -t_zzz_xzzz[i] * ab_x[i] + t_zzz_xxzzz[i];

        t_xzzz_yyyy[i] = -t_zzz_yyyy[i] * ab_x[i] + t_zzz_xyyyy[i];

        t_xzzz_yyyz[i] = -t_zzz_yyyz[i] * ab_x[i] + t_zzz_xyyyz[i];

        t_xzzz_yyzz[i] = -t_zzz_yyzz[i] * ab_x[i] + t_zzz_xyyzz[i];

        t_xzzz_yzzz[i] = -t_zzz_yzzz[i] * ab_x[i] + t_zzz_xyzzz[i];

        t_xzzz_zzzz[i] = -t_zzz_zzzz[i] * ab_x[i] + t_zzz_xzzzz[i];

        t_yyyy_xxxx[i] = -t_yyy_xxxx[i] * ab_y[i] + t_yyy_xxxxy[i];

        t_yyyy_xxxy[i] = -t_yyy_xxxy[i] * ab_y[i] + t_yyy_xxxyy[i];

        t_yyyy_xxxz[i] = -t_yyy_xxxz[i] * ab_y[i] + t_yyy_xxxyz[i];

        t_yyyy_xxyy[i] = -t_yyy_xxyy[i] * ab_y[i] + t_yyy_xxyyy[i];

        t_yyyy_xxyz[i] = -t_yyy_xxyz[i] * ab_y[i] + t_yyy_xxyyz[i];

        t_yyyy_xxzz[i] = -t_yyy_xxzz[i] * ab_y[i] + t_yyy_xxyzz[i];

        t_yyyy_xyyy[i] = -t_yyy_xyyy[i] * ab_y[i] + t_yyy_xyyyy[i];

        t_yyyy_xyyz[i] = -t_yyy_xyyz[i] * ab_y[i] + t_yyy_xyyyz[i];

        t_yyyy_xyzz[i] = -t_yyy_xyzz[i] * ab_y[i] + t_yyy_xyyzz[i];

        t_yyyy_xzzz[i] = -t_yyy_xzzz[i] * ab_y[i] + t_yyy_xyzzz[i];

        t_yyyy_yyyy[i] = -t_yyy_yyyy[i] * ab_y[i] + t_yyy_yyyyy[i];

        t_yyyy_yyyz[i] = -t_yyy_yyyz[i] * ab_y[i] + t_yyy_yyyyz[i];

        t_yyyy_yyzz[i] = -t_yyy_yyzz[i] * ab_y[i] + t_yyy_yyyzz[i];

        t_yyyy_yzzz[i] = -t_yyy_yzzz[i] * ab_y[i] + t_yyy_yyzzz[i];

        t_yyyy_zzzz[i] = -t_yyy_zzzz[i] * ab_y[i] + t_yyy_yzzzz[i];

        t_yyyz_xxxx[i] = -t_yyz_xxxx[i] * ab_y[i] + t_yyz_xxxxy[i];

        t_yyyz_xxxy[i] = -t_yyz_xxxy[i] * ab_y[i] + t_yyz_xxxyy[i];

        t_yyyz_xxxz[i] = -t_yyz_xxxz[i] * ab_y[i] + t_yyz_xxxyz[i];

        t_yyyz_xxyy[i] = -t_yyz_xxyy[i] * ab_y[i] + t_yyz_xxyyy[i];

        t_yyyz_xxyz[i] = -t_yyz_xxyz[i] * ab_y[i] + t_yyz_xxyyz[i];

        t_yyyz_xxzz[i] = -t_yyz_xxzz[i] * ab_y[i] + t_yyz_xxyzz[i];

        t_yyyz_xyyy[i] = -t_yyz_xyyy[i] * ab_y[i] + t_yyz_xyyyy[i];

        t_yyyz_xyyz[i] = -t_yyz_xyyz[i] * ab_y[i] + t_yyz_xyyyz[i];

        t_yyyz_xyzz[i] = -t_yyz_xyzz[i] * ab_y[i] + t_yyz_xyyzz[i];

        t_yyyz_xzzz[i] = -t_yyz_xzzz[i] * ab_y[i] + t_yyz_xyzzz[i];

        t_yyyz_yyyy[i] = -t_yyz_yyyy[i] * ab_y[i] + t_yyz_yyyyy[i];

        t_yyyz_yyyz[i] = -t_yyz_yyyz[i] * ab_y[i] + t_yyz_yyyyz[i];

        t_yyyz_yyzz[i] = -t_yyz_yyzz[i] * ab_y[i] + t_yyz_yyyzz[i];

        t_yyyz_yzzz[i] = -t_yyz_yzzz[i] * ab_y[i] + t_yyz_yyzzz[i];

        t_yyyz_zzzz[i] = -t_yyz_zzzz[i] * ab_y[i] + t_yyz_yzzzz[i];

        t_yyzz_xxxx[i] = -t_yzz_xxxx[i] * ab_y[i] + t_yzz_xxxxy[i];

        t_yyzz_xxxy[i] = -t_yzz_xxxy[i] * ab_y[i] + t_yzz_xxxyy[i];

        t_yyzz_xxxz[i] = -t_yzz_xxxz[i] * ab_y[i] + t_yzz_xxxyz[i];

        t_yyzz_xxyy[i] = -t_yzz_xxyy[i] * ab_y[i] + t_yzz_xxyyy[i];

        t_yyzz_xxyz[i] = -t_yzz_xxyz[i] * ab_y[i] + t_yzz_xxyyz[i];

        t_yyzz_xxzz[i] = -t_yzz_xxzz[i] * ab_y[i] + t_yzz_xxyzz[i];

        t_yyzz_xyyy[i] = -t_yzz_xyyy[i] * ab_y[i] + t_yzz_xyyyy[i];

        t_yyzz_xyyz[i] = -t_yzz_xyyz[i] * ab_y[i] + t_yzz_xyyyz[i];

        t_yyzz_xyzz[i] = -t_yzz_xyzz[i] * ab_y[i] + t_yzz_xyyzz[i];

        t_yyzz_xzzz[i] = -t_yzz_xzzz[i] * ab_y[i] + t_yzz_xyzzz[i];

        t_yyzz_yyyy[i] = -t_yzz_yyyy[i] * ab_y[i] + t_yzz_yyyyy[i];

        t_yyzz_yyyz[i] = -t_yzz_yyyz[i] * ab_y[i] + t_yzz_yyyyz[i];

        t_yyzz_yyzz[i] = -t_yzz_yyzz[i] * ab_y[i] + t_yzz_yyyzz[i];

        t_yyzz_yzzz[i] = -t_yzz_yzzz[i] * ab_y[i] + t_yzz_yyzzz[i];

        t_yyzz_zzzz[i] = -t_yzz_zzzz[i] * ab_y[i] + t_yzz_yzzzz[i];

        t_yzzz_xxxx[i] = -t_zzz_xxxx[i] * ab_y[i] + t_zzz_xxxxy[i];

        t_yzzz_xxxy[i] = -t_zzz_xxxy[i] * ab_y[i] + t_zzz_xxxyy[i];

        t_yzzz_xxxz[i] = -t_zzz_xxxz[i] * ab_y[i] + t_zzz_xxxyz[i];

        t_yzzz_xxyy[i] = -t_zzz_xxyy[i] * ab_y[i] + t_zzz_xxyyy[i];

        t_yzzz_xxyz[i] = -t_zzz_xxyz[i] * ab_y[i] + t_zzz_xxyyz[i];

        t_yzzz_xxzz[i] = -t_zzz_xxzz[i] * ab_y[i] + t_zzz_xxyzz[i];

        t_yzzz_xyyy[i] = -t_zzz_xyyy[i] * ab_y[i] + t_zzz_xyyyy[i];

        t_yzzz_xyyz[i] = -t_zzz_xyyz[i] * ab_y[i] + t_zzz_xyyyz[i];

        t_yzzz_xyzz[i] = -t_zzz_xyzz[i] * ab_y[i] + t_zzz_xyyzz[i];

        t_yzzz_xzzz[i] = -t_zzz_xzzz[i] * ab_y[i] + t_zzz_xyzzz[i];

        t_yzzz_yyyy[i] = -t_zzz_yyyy[i] * ab_y[i] + t_zzz_yyyyy[i];

        t_yzzz_yyyz[i] = -t_zzz_yyyz[i] * ab_y[i] + t_zzz_yyyyz[i];

        t_yzzz_yyzz[i] = -t_zzz_yyzz[i] * ab_y[i] + t_zzz_yyyzz[i];

        t_yzzz_yzzz[i] = -t_zzz_yzzz[i] * ab_y[i] + t_zzz_yyzzz[i];

        t_yzzz_zzzz[i] = -t_zzz_zzzz[i] * ab_y[i] + t_zzz_yzzzz[i];

        t_zzzz_xxxx[i] = -t_zzz_xxxx[i] * ab_z[i] + t_zzz_xxxxz[i];

        t_zzzz_xxxy[i] = -t_zzz_xxxy[i] * ab_z[i] + t_zzz_xxxyz[i];

        t_zzzz_xxxz[i] = -t_zzz_xxxz[i] * ab_z[i] + t_zzz_xxxzz[i];

        t_zzzz_xxyy[i] = -t_zzz_xxyy[i] * ab_z[i] + t_zzz_xxyyz[i];

        t_zzzz_xxyz[i] = -t_zzz_xxyz[i] * ab_z[i] + t_zzz_xxyzz[i];

        t_zzzz_xxzz[i] = -t_zzz_xxzz[i] * ab_z[i] + t_zzz_xxzzz[i];

        t_zzzz_xyyy[i] = -t_zzz_xyyy[i] * ab_z[i] + t_zzz_xyyyz[i];

        t_zzzz_xyyz[i] = -t_zzz_xyyz[i] * ab_z[i] + t_zzz_xyyzz[i];

        t_zzzz_xyzz[i] = -t_zzz_xyzz[i] * ab_z[i] + t_zzz_xyzzz[i];

        t_zzzz_xzzz[i] = -t_zzz_xzzz[i] * ab_z[i] + t_zzz_xzzzz[i];

        t_zzzz_yyyy[i] = -t_zzz_yyyy[i] * ab_z[i] + t_zzz_yyyyz[i];

        t_zzzz_yyyz[i] = -t_zzz_yyyz[i] * ab_z[i] + t_zzz_yyyzz[i];

        t_zzzz_yyzz[i] = -t_zzz_yyzz[i] * ab_z[i] + t_zzz_yyzzz[i];

        t_zzzz_yzzz[i] = -t_zzz_yzzz[i] * ab_z[i] + t_zzz_yzzzz[i];

        t_zzzz_zzzz[i] = -t_zzz_zzzz[i] * ab_z[i] + t_zzz_zzzzz[i];
    }
}

} // t2chrr namespace

