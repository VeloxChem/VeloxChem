#include "KineticEnergyPrimRecHF.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_hf(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_hf,
                            const size_t idx_ovl_ff,
                            const size_t idx_kin_ff,
                            const size_t idx_kin_gd,
                            const size_t idx_kin_gf,
                            const size_t idx_ovl_hf,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ovl_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ovl_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ovl_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ovl_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ovl_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ovl_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ovl_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ovl_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ovl_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ovl_ff + 9);

    auto ts_xxy_xxx = pbuffer.data(idx_ovl_ff + 10);

    auto ts_xxy_xxz = pbuffer.data(idx_ovl_ff + 12);

    auto ts_xxy_xzz = pbuffer.data(idx_ovl_ff + 15);

    auto ts_xxz_xxx = pbuffer.data(idx_ovl_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ovl_ff + 21);

    auto ts_xxz_xyy = pbuffer.data(idx_ovl_ff + 23);

    auto ts_xyy_xxy = pbuffer.data(idx_ovl_ff + 31);

    auto ts_xyy_xyy = pbuffer.data(idx_ovl_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ovl_ff + 34);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ovl_ff + 39);

    auto ts_xzz_xxz = pbuffer.data(idx_ovl_ff + 52);

    auto ts_xzz_xyz = pbuffer.data(idx_ovl_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ovl_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ovl_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ovl_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ovl_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ovl_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ovl_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ovl_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ovl_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ovl_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ovl_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ovl_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ovl_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ovl_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ovl_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ovl_ff + 69);

    auto ts_yyz_xxy = pbuffer.data(idx_ovl_ff + 71);

    auto ts_yyz_xyy = pbuffer.data(idx_ovl_ff + 73);

    auto ts_yyz_yyy = pbuffer.data(idx_ovl_ff + 76);

    auto ts_yzz_xxx = pbuffer.data(idx_ovl_ff + 80);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xyz = pbuffer.data(idx_ovl_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyz = pbuffer.data(idx_ovl_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ovl_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ovl_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ovl_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ovl_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ovl_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ovl_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ovl_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ovl_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ovl_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ovl_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ovl_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ovl_ff + 99);

    // Set up components of auxiliary buffer : FF

    auto tk_xxx_xxx = pbuffer.data(idx_kin_ff);

    auto tk_xxx_xxy = pbuffer.data(idx_kin_ff + 1);

    auto tk_xxx_xxz = pbuffer.data(idx_kin_ff + 2);

    auto tk_xxx_xyy = pbuffer.data(idx_kin_ff + 3);

    auto tk_xxx_xyz = pbuffer.data(idx_kin_ff + 4);

    auto tk_xxx_xzz = pbuffer.data(idx_kin_ff + 5);

    auto tk_xxx_yyy = pbuffer.data(idx_kin_ff + 6);

    auto tk_xxx_yyz = pbuffer.data(idx_kin_ff + 7);

    auto tk_xxx_yzz = pbuffer.data(idx_kin_ff + 8);

    auto tk_xxx_zzz = pbuffer.data(idx_kin_ff + 9);

    auto tk_xxy_xxx = pbuffer.data(idx_kin_ff + 10);

    auto tk_xxy_xxz = pbuffer.data(idx_kin_ff + 12);

    auto tk_xxy_xzz = pbuffer.data(idx_kin_ff + 15);

    auto tk_xxz_xxx = pbuffer.data(idx_kin_ff + 20);

    auto tk_xxz_xxy = pbuffer.data(idx_kin_ff + 21);

    auto tk_xxz_xyy = pbuffer.data(idx_kin_ff + 23);

    auto tk_xyy_xxy = pbuffer.data(idx_kin_ff + 31);

    auto tk_xyy_xyy = pbuffer.data(idx_kin_ff + 33);

    auto tk_xyy_xyz = pbuffer.data(idx_kin_ff + 34);

    auto tk_xyy_yyy = pbuffer.data(idx_kin_ff + 36);

    auto tk_xyy_yyz = pbuffer.data(idx_kin_ff + 37);

    auto tk_xyy_yzz = pbuffer.data(idx_kin_ff + 38);

    auto tk_xyy_zzz = pbuffer.data(idx_kin_ff + 39);

    auto tk_xzz_xxz = pbuffer.data(idx_kin_ff + 52);

    auto tk_xzz_xyz = pbuffer.data(idx_kin_ff + 54);

    auto tk_xzz_xzz = pbuffer.data(idx_kin_ff + 55);

    auto tk_xzz_yyy = pbuffer.data(idx_kin_ff + 56);

    auto tk_xzz_yyz = pbuffer.data(idx_kin_ff + 57);

    auto tk_xzz_yzz = pbuffer.data(idx_kin_ff + 58);

    auto tk_xzz_zzz = pbuffer.data(idx_kin_ff + 59);

    auto tk_yyy_xxx = pbuffer.data(idx_kin_ff + 60);

    auto tk_yyy_xxy = pbuffer.data(idx_kin_ff + 61);

    auto tk_yyy_xxz = pbuffer.data(idx_kin_ff + 62);

    auto tk_yyy_xyy = pbuffer.data(idx_kin_ff + 63);

    auto tk_yyy_xyz = pbuffer.data(idx_kin_ff + 64);

    auto tk_yyy_xzz = pbuffer.data(idx_kin_ff + 65);

    auto tk_yyy_yyy = pbuffer.data(idx_kin_ff + 66);

    auto tk_yyy_yyz = pbuffer.data(idx_kin_ff + 67);

    auto tk_yyy_yzz = pbuffer.data(idx_kin_ff + 68);

    auto tk_yyy_zzz = pbuffer.data(idx_kin_ff + 69);

    auto tk_yyz_xxy = pbuffer.data(idx_kin_ff + 71);

    auto tk_yyz_xyy = pbuffer.data(idx_kin_ff + 73);

    auto tk_yyz_yyy = pbuffer.data(idx_kin_ff + 76);

    auto tk_yzz_xxx = pbuffer.data(idx_kin_ff + 80);

    auto tk_yzz_xxz = pbuffer.data(idx_kin_ff + 82);

    auto tk_yzz_xyz = pbuffer.data(idx_kin_ff + 84);

    auto tk_yzz_xzz = pbuffer.data(idx_kin_ff + 85);

    auto tk_yzz_yyz = pbuffer.data(idx_kin_ff + 87);

    auto tk_yzz_yzz = pbuffer.data(idx_kin_ff + 88);

    auto tk_yzz_zzz = pbuffer.data(idx_kin_ff + 89);

    auto tk_zzz_xxx = pbuffer.data(idx_kin_ff + 90);

    auto tk_zzz_xxy = pbuffer.data(idx_kin_ff + 91);

    auto tk_zzz_xxz = pbuffer.data(idx_kin_ff + 92);

    auto tk_zzz_xyy = pbuffer.data(idx_kin_ff + 93);

    auto tk_zzz_xyz = pbuffer.data(idx_kin_ff + 94);

    auto tk_zzz_xzz = pbuffer.data(idx_kin_ff + 95);

    auto tk_zzz_yyy = pbuffer.data(idx_kin_ff + 96);

    auto tk_zzz_yyz = pbuffer.data(idx_kin_ff + 97);

    auto tk_zzz_yzz = pbuffer.data(idx_kin_ff + 98);

    auto tk_zzz_zzz = pbuffer.data(idx_kin_ff + 99);

    // Set up components of auxiliary buffer : GD

    auto tk_xxxx_xx = pbuffer.data(idx_kin_gd);

    auto tk_xxxx_xy = pbuffer.data(idx_kin_gd + 1);

    auto tk_xxxx_xz = pbuffer.data(idx_kin_gd + 2);

    auto tk_xxxx_yy = pbuffer.data(idx_kin_gd + 3);

    auto tk_xxxx_yz = pbuffer.data(idx_kin_gd + 4);

    auto tk_xxxx_zz = pbuffer.data(idx_kin_gd + 5);

    auto tk_xxxz_xz = pbuffer.data(idx_kin_gd + 14);

    auto tk_xxxz_yz = pbuffer.data(idx_kin_gd + 16);

    auto tk_xxxz_zz = pbuffer.data(idx_kin_gd + 17);

    auto tk_xxyy_xx = pbuffer.data(idx_kin_gd + 18);

    auto tk_xxyy_xy = pbuffer.data(idx_kin_gd + 19);

    auto tk_xxyy_xz = pbuffer.data(idx_kin_gd + 20);

    auto tk_xxyy_yy = pbuffer.data(idx_kin_gd + 21);

    auto tk_xxyy_yz = pbuffer.data(idx_kin_gd + 22);

    auto tk_xxyy_zz = pbuffer.data(idx_kin_gd + 23);

    auto tk_xxzz_xx = pbuffer.data(idx_kin_gd + 30);

    auto tk_xxzz_xy = pbuffer.data(idx_kin_gd + 31);

    auto tk_xxzz_xz = pbuffer.data(idx_kin_gd + 32);

    auto tk_xxzz_yy = pbuffer.data(idx_kin_gd + 33);

    auto tk_xxzz_yz = pbuffer.data(idx_kin_gd + 34);

    auto tk_xxzz_zz = pbuffer.data(idx_kin_gd + 35);

    auto tk_xyyy_xy = pbuffer.data(idx_kin_gd + 37);

    auto tk_xyyy_yy = pbuffer.data(idx_kin_gd + 39);

    auto tk_xyyy_yz = pbuffer.data(idx_kin_gd + 40);

    auto tk_xzzz_xz = pbuffer.data(idx_kin_gd + 56);

    auto tk_xzzz_yz = pbuffer.data(idx_kin_gd + 58);

    auto tk_xzzz_zz = pbuffer.data(idx_kin_gd + 59);

    auto tk_yyyy_xx = pbuffer.data(idx_kin_gd + 60);

    auto tk_yyyy_xy = pbuffer.data(idx_kin_gd + 61);

    auto tk_yyyy_xz = pbuffer.data(idx_kin_gd + 62);

    auto tk_yyyy_yy = pbuffer.data(idx_kin_gd + 63);

    auto tk_yyyy_yz = pbuffer.data(idx_kin_gd + 64);

    auto tk_yyyy_zz = pbuffer.data(idx_kin_gd + 65);

    auto tk_yyyz_xz = pbuffer.data(idx_kin_gd + 68);

    auto tk_yyyz_yz = pbuffer.data(idx_kin_gd + 70);

    auto tk_yyyz_zz = pbuffer.data(idx_kin_gd + 71);

    auto tk_yyzz_xx = pbuffer.data(idx_kin_gd + 72);

    auto tk_yyzz_xy = pbuffer.data(idx_kin_gd + 73);

    auto tk_yyzz_xz = pbuffer.data(idx_kin_gd + 74);

    auto tk_yyzz_yy = pbuffer.data(idx_kin_gd + 75);

    auto tk_yyzz_yz = pbuffer.data(idx_kin_gd + 76);

    auto tk_yyzz_zz = pbuffer.data(idx_kin_gd + 77);

    auto tk_yzzz_xy = pbuffer.data(idx_kin_gd + 79);

    auto tk_yzzz_xz = pbuffer.data(idx_kin_gd + 80);

    auto tk_yzzz_yy = pbuffer.data(idx_kin_gd + 81);

    auto tk_yzzz_yz = pbuffer.data(idx_kin_gd + 82);

    auto tk_yzzz_zz = pbuffer.data(idx_kin_gd + 83);

    auto tk_zzzz_xx = pbuffer.data(idx_kin_gd + 84);

    auto tk_zzzz_xy = pbuffer.data(idx_kin_gd + 85);

    auto tk_zzzz_xz = pbuffer.data(idx_kin_gd + 86);

    auto tk_zzzz_yy = pbuffer.data(idx_kin_gd + 87);

    auto tk_zzzz_yz = pbuffer.data(idx_kin_gd + 88);

    auto tk_zzzz_zz = pbuffer.data(idx_kin_gd + 89);

    // Set up components of auxiliary buffer : GF

    auto tk_xxxx_xxx = pbuffer.data(idx_kin_gf);

    auto tk_xxxx_xxy = pbuffer.data(idx_kin_gf + 1);

    auto tk_xxxx_xxz = pbuffer.data(idx_kin_gf + 2);

    auto tk_xxxx_xyy = pbuffer.data(idx_kin_gf + 3);

    auto tk_xxxx_xyz = pbuffer.data(idx_kin_gf + 4);

    auto tk_xxxx_xzz = pbuffer.data(idx_kin_gf + 5);

    auto tk_xxxx_yyy = pbuffer.data(idx_kin_gf + 6);

    auto tk_xxxx_yyz = pbuffer.data(idx_kin_gf + 7);

    auto tk_xxxx_yzz = pbuffer.data(idx_kin_gf + 8);

    auto tk_xxxx_zzz = pbuffer.data(idx_kin_gf + 9);

    auto tk_xxxy_xxx = pbuffer.data(idx_kin_gf + 10);

    auto tk_xxxy_xxy = pbuffer.data(idx_kin_gf + 11);

    auto tk_xxxy_xxz = pbuffer.data(idx_kin_gf + 12);

    auto tk_xxxy_xyy = pbuffer.data(idx_kin_gf + 13);

    auto tk_xxxy_xzz = pbuffer.data(idx_kin_gf + 15);

    auto tk_xxxy_yyy = pbuffer.data(idx_kin_gf + 16);

    auto tk_xxxz_xxx = pbuffer.data(idx_kin_gf + 20);

    auto tk_xxxz_xxy = pbuffer.data(idx_kin_gf + 21);

    auto tk_xxxz_xxz = pbuffer.data(idx_kin_gf + 22);

    auto tk_xxxz_xyy = pbuffer.data(idx_kin_gf + 23);

    auto tk_xxxz_xyz = pbuffer.data(idx_kin_gf + 24);

    auto tk_xxxz_xzz = pbuffer.data(idx_kin_gf + 25);

    auto tk_xxxz_yyz = pbuffer.data(idx_kin_gf + 27);

    auto tk_xxxz_yzz = pbuffer.data(idx_kin_gf + 28);

    auto tk_xxxz_zzz = pbuffer.data(idx_kin_gf + 29);

    auto tk_xxyy_xxx = pbuffer.data(idx_kin_gf + 30);

    auto tk_xxyy_xxy = pbuffer.data(idx_kin_gf + 31);

    auto tk_xxyy_xxz = pbuffer.data(idx_kin_gf + 32);

    auto tk_xxyy_xyy = pbuffer.data(idx_kin_gf + 33);

    auto tk_xxyy_xyz = pbuffer.data(idx_kin_gf + 34);

    auto tk_xxyy_xzz = pbuffer.data(idx_kin_gf + 35);

    auto tk_xxyy_yyy = pbuffer.data(idx_kin_gf + 36);

    auto tk_xxyy_yyz = pbuffer.data(idx_kin_gf + 37);

    auto tk_xxyy_yzz = pbuffer.data(idx_kin_gf + 38);

    auto tk_xxyy_zzz = pbuffer.data(idx_kin_gf + 39);

    auto tk_xxzz_xxx = pbuffer.data(idx_kin_gf + 50);

    auto tk_xxzz_xxy = pbuffer.data(idx_kin_gf + 51);

    auto tk_xxzz_xxz = pbuffer.data(idx_kin_gf + 52);

    auto tk_xxzz_xyy = pbuffer.data(idx_kin_gf + 53);

    auto tk_xxzz_xyz = pbuffer.data(idx_kin_gf + 54);

    auto tk_xxzz_xzz = pbuffer.data(idx_kin_gf + 55);

    auto tk_xxzz_yyy = pbuffer.data(idx_kin_gf + 56);

    auto tk_xxzz_yyz = pbuffer.data(idx_kin_gf + 57);

    auto tk_xxzz_yzz = pbuffer.data(idx_kin_gf + 58);

    auto tk_xxzz_zzz = pbuffer.data(idx_kin_gf + 59);

    auto tk_xyyy_xxx = pbuffer.data(idx_kin_gf + 60);

    auto tk_xyyy_xxy = pbuffer.data(idx_kin_gf + 61);

    auto tk_xyyy_xyy = pbuffer.data(idx_kin_gf + 63);

    auto tk_xyyy_xyz = pbuffer.data(idx_kin_gf + 64);

    auto tk_xyyy_yyy = pbuffer.data(idx_kin_gf + 66);

    auto tk_xyyy_yyz = pbuffer.data(idx_kin_gf + 67);

    auto tk_xyyy_yzz = pbuffer.data(idx_kin_gf + 68);

    auto tk_xyyy_zzz = pbuffer.data(idx_kin_gf + 69);

    auto tk_xzzz_xxx = pbuffer.data(idx_kin_gf + 90);

    auto tk_xzzz_xxz = pbuffer.data(idx_kin_gf + 92);

    auto tk_xzzz_xyz = pbuffer.data(idx_kin_gf + 94);

    auto tk_xzzz_xzz = pbuffer.data(idx_kin_gf + 95);

    auto tk_xzzz_yyy = pbuffer.data(idx_kin_gf + 96);

    auto tk_xzzz_yyz = pbuffer.data(idx_kin_gf + 97);

    auto tk_xzzz_yzz = pbuffer.data(idx_kin_gf + 98);

    auto tk_xzzz_zzz = pbuffer.data(idx_kin_gf + 99);

    auto tk_yyyy_xxx = pbuffer.data(idx_kin_gf + 100);

    auto tk_yyyy_xxy = pbuffer.data(idx_kin_gf + 101);

    auto tk_yyyy_xxz = pbuffer.data(idx_kin_gf + 102);

    auto tk_yyyy_xyy = pbuffer.data(idx_kin_gf + 103);

    auto tk_yyyy_xyz = pbuffer.data(idx_kin_gf + 104);

    auto tk_yyyy_xzz = pbuffer.data(idx_kin_gf + 105);

    auto tk_yyyy_yyy = pbuffer.data(idx_kin_gf + 106);

    auto tk_yyyy_yyz = pbuffer.data(idx_kin_gf + 107);

    auto tk_yyyy_yzz = pbuffer.data(idx_kin_gf + 108);

    auto tk_yyyy_zzz = pbuffer.data(idx_kin_gf + 109);

    auto tk_yyyz_xxy = pbuffer.data(idx_kin_gf + 111);

    auto tk_yyyz_xxz = pbuffer.data(idx_kin_gf + 112);

    auto tk_yyyz_xyy = pbuffer.data(idx_kin_gf + 113);

    auto tk_yyyz_xyz = pbuffer.data(idx_kin_gf + 114);

    auto tk_yyyz_xzz = pbuffer.data(idx_kin_gf + 115);

    auto tk_yyyz_yyy = pbuffer.data(idx_kin_gf + 116);

    auto tk_yyyz_yyz = pbuffer.data(idx_kin_gf + 117);

    auto tk_yyyz_yzz = pbuffer.data(idx_kin_gf + 118);

    auto tk_yyyz_zzz = pbuffer.data(idx_kin_gf + 119);

    auto tk_yyzz_xxx = pbuffer.data(idx_kin_gf + 120);

    auto tk_yyzz_xxy = pbuffer.data(idx_kin_gf + 121);

    auto tk_yyzz_xxz = pbuffer.data(idx_kin_gf + 122);

    auto tk_yyzz_xyy = pbuffer.data(idx_kin_gf + 123);

    auto tk_yyzz_xyz = pbuffer.data(idx_kin_gf + 124);

    auto tk_yyzz_xzz = pbuffer.data(idx_kin_gf + 125);

    auto tk_yyzz_yyy = pbuffer.data(idx_kin_gf + 126);

    auto tk_yyzz_yyz = pbuffer.data(idx_kin_gf + 127);

    auto tk_yyzz_yzz = pbuffer.data(idx_kin_gf + 128);

    auto tk_yyzz_zzz = pbuffer.data(idx_kin_gf + 129);

    auto tk_yzzz_xxx = pbuffer.data(idx_kin_gf + 130);

    auto tk_yzzz_xxy = pbuffer.data(idx_kin_gf + 131);

    auto tk_yzzz_xxz = pbuffer.data(idx_kin_gf + 132);

    auto tk_yzzz_xyy = pbuffer.data(idx_kin_gf + 133);

    auto tk_yzzz_xyz = pbuffer.data(idx_kin_gf + 134);

    auto tk_yzzz_xzz = pbuffer.data(idx_kin_gf + 135);

    auto tk_yzzz_yyy = pbuffer.data(idx_kin_gf + 136);

    auto tk_yzzz_yyz = pbuffer.data(idx_kin_gf + 137);

    auto tk_yzzz_yzz = pbuffer.data(idx_kin_gf + 138);

    auto tk_yzzz_zzz = pbuffer.data(idx_kin_gf + 139);

    auto tk_zzzz_xxx = pbuffer.data(idx_kin_gf + 140);

    auto tk_zzzz_xxy = pbuffer.data(idx_kin_gf + 141);

    auto tk_zzzz_xxz = pbuffer.data(idx_kin_gf + 142);

    auto tk_zzzz_xyy = pbuffer.data(idx_kin_gf + 143);

    auto tk_zzzz_xyz = pbuffer.data(idx_kin_gf + 144);

    auto tk_zzzz_xzz = pbuffer.data(idx_kin_gf + 145);

    auto tk_zzzz_yyy = pbuffer.data(idx_kin_gf + 146);

    auto tk_zzzz_yyz = pbuffer.data(idx_kin_gf + 147);

    auto tk_zzzz_yzz = pbuffer.data(idx_kin_gf + 148);

    auto tk_zzzz_zzz = pbuffer.data(idx_kin_gf + 149);

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_ovl_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_ovl_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_ovl_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_ovl_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_ovl_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_ovl_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_ovl_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_ovl_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_ovl_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_ovl_hf + 9);

    auto ts_xxxxy_xxx = pbuffer.data(idx_ovl_hf + 10);

    auto ts_xxxxy_xxy = pbuffer.data(idx_ovl_hf + 11);

    auto ts_xxxxy_xxz = pbuffer.data(idx_ovl_hf + 12);

    auto ts_xxxxy_xyy = pbuffer.data(idx_ovl_hf + 13);

    auto ts_xxxxy_xyz = pbuffer.data(idx_ovl_hf + 14);

    auto ts_xxxxy_xzz = pbuffer.data(idx_ovl_hf + 15);

    auto ts_xxxxy_yyy = pbuffer.data(idx_ovl_hf + 16);

    auto ts_xxxxy_yyz = pbuffer.data(idx_ovl_hf + 17);

    auto ts_xxxxy_yzz = pbuffer.data(idx_ovl_hf + 18);

    auto ts_xxxxy_zzz = pbuffer.data(idx_ovl_hf + 19);

    auto ts_xxxxz_xxx = pbuffer.data(idx_ovl_hf + 20);

    auto ts_xxxxz_xxy = pbuffer.data(idx_ovl_hf + 21);

    auto ts_xxxxz_xxz = pbuffer.data(idx_ovl_hf + 22);

    auto ts_xxxxz_xyy = pbuffer.data(idx_ovl_hf + 23);

    auto ts_xxxxz_xyz = pbuffer.data(idx_ovl_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_ovl_hf + 25);

    auto ts_xxxxz_yyy = pbuffer.data(idx_ovl_hf + 26);

    auto ts_xxxxz_yyz = pbuffer.data(idx_ovl_hf + 27);

    auto ts_xxxxz_yzz = pbuffer.data(idx_ovl_hf + 28);

    auto ts_xxxxz_zzz = pbuffer.data(idx_ovl_hf + 29);

    auto ts_xxxyy_xxx = pbuffer.data(idx_ovl_hf + 30);

    auto ts_xxxyy_xxy = pbuffer.data(idx_ovl_hf + 31);

    auto ts_xxxyy_xxz = pbuffer.data(idx_ovl_hf + 32);

    auto ts_xxxyy_xyy = pbuffer.data(idx_ovl_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_ovl_hf + 34);

    auto ts_xxxyy_xzz = pbuffer.data(idx_ovl_hf + 35);

    auto ts_xxxyy_yyy = pbuffer.data(idx_ovl_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_ovl_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_ovl_hf + 38);

    auto ts_xxxyy_zzz = pbuffer.data(idx_ovl_hf + 39);

    auto ts_xxxyz_xxx = pbuffer.data(idx_ovl_hf + 40);

    auto ts_xxxyz_xxy = pbuffer.data(idx_ovl_hf + 41);

    auto ts_xxxyz_xxz = pbuffer.data(idx_ovl_hf + 42);

    auto ts_xxxyz_xyy = pbuffer.data(idx_ovl_hf + 43);

    auto ts_xxxyz_xyz = pbuffer.data(idx_ovl_hf + 44);

    auto ts_xxxyz_xzz = pbuffer.data(idx_ovl_hf + 45);

    auto ts_xxxyz_yyy = pbuffer.data(idx_ovl_hf + 46);

    auto ts_xxxyz_yyz = pbuffer.data(idx_ovl_hf + 47);

    auto ts_xxxyz_yzz = pbuffer.data(idx_ovl_hf + 48);

    auto ts_xxxyz_zzz = pbuffer.data(idx_ovl_hf + 49);

    auto ts_xxxzz_xxx = pbuffer.data(idx_ovl_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_ovl_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_ovl_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_ovl_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_ovl_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_ovl_hf + 55);

    auto ts_xxxzz_yyy = pbuffer.data(idx_ovl_hf + 56);

    auto ts_xxxzz_yyz = pbuffer.data(idx_ovl_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_ovl_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_ovl_hf + 59);

    auto ts_xxyyy_xxx = pbuffer.data(idx_ovl_hf + 60);

    auto ts_xxyyy_xxy = pbuffer.data(idx_ovl_hf + 61);

    auto ts_xxyyy_xxz = pbuffer.data(idx_ovl_hf + 62);

    auto ts_xxyyy_xyy = pbuffer.data(idx_ovl_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_ovl_hf + 64);

    auto ts_xxyyy_xzz = pbuffer.data(idx_ovl_hf + 65);

    auto ts_xxyyy_yyy = pbuffer.data(idx_ovl_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_ovl_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_ovl_hf + 68);

    auto ts_xxyyy_zzz = pbuffer.data(idx_ovl_hf + 69);

    auto ts_xxyyz_xxx = pbuffer.data(idx_ovl_hf + 70);

    auto ts_xxyyz_xxy = pbuffer.data(idx_ovl_hf + 71);

    auto ts_xxyyz_xxz = pbuffer.data(idx_ovl_hf + 72);

    auto ts_xxyyz_xyy = pbuffer.data(idx_ovl_hf + 73);

    auto ts_xxyyz_xyz = pbuffer.data(idx_ovl_hf + 74);

    auto ts_xxyyz_xzz = pbuffer.data(idx_ovl_hf + 75);

    auto ts_xxyyz_yyy = pbuffer.data(idx_ovl_hf + 76);

    auto ts_xxyyz_yyz = pbuffer.data(idx_ovl_hf + 77);

    auto ts_xxyyz_yzz = pbuffer.data(idx_ovl_hf + 78);

    auto ts_xxyyz_zzz = pbuffer.data(idx_ovl_hf + 79);

    auto ts_xxyzz_xxx = pbuffer.data(idx_ovl_hf + 80);

    auto ts_xxyzz_xxy = pbuffer.data(idx_ovl_hf + 81);

    auto ts_xxyzz_xxz = pbuffer.data(idx_ovl_hf + 82);

    auto ts_xxyzz_xyy = pbuffer.data(idx_ovl_hf + 83);

    auto ts_xxyzz_xyz = pbuffer.data(idx_ovl_hf + 84);

    auto ts_xxyzz_xzz = pbuffer.data(idx_ovl_hf + 85);

    auto ts_xxyzz_yyy = pbuffer.data(idx_ovl_hf + 86);

    auto ts_xxyzz_yyz = pbuffer.data(idx_ovl_hf + 87);

    auto ts_xxyzz_yzz = pbuffer.data(idx_ovl_hf + 88);

    auto ts_xxyzz_zzz = pbuffer.data(idx_ovl_hf + 89);

    auto ts_xxzzz_xxx = pbuffer.data(idx_ovl_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_ovl_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_ovl_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_ovl_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_ovl_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_ovl_hf + 95);

    auto ts_xxzzz_yyy = pbuffer.data(idx_ovl_hf + 96);

    auto ts_xxzzz_yyz = pbuffer.data(idx_ovl_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_ovl_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_ovl_hf + 99);

    auto ts_xyyyy_xxx = pbuffer.data(idx_ovl_hf + 100);

    auto ts_xyyyy_xxy = pbuffer.data(idx_ovl_hf + 101);

    auto ts_xyyyy_xxz = pbuffer.data(idx_ovl_hf + 102);

    auto ts_xyyyy_xyy = pbuffer.data(idx_ovl_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_ovl_hf + 104);

    auto ts_xyyyy_xzz = pbuffer.data(idx_ovl_hf + 105);

    auto ts_xyyyy_yyy = pbuffer.data(idx_ovl_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_ovl_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_ovl_hf + 108);

    auto ts_xyyyy_zzz = pbuffer.data(idx_ovl_hf + 109);

    auto ts_xyyyz_xxx = pbuffer.data(idx_ovl_hf + 110);

    auto ts_xyyyz_xxy = pbuffer.data(idx_ovl_hf + 111);

    auto ts_xyyyz_xxz = pbuffer.data(idx_ovl_hf + 112);

    auto ts_xyyyz_xyy = pbuffer.data(idx_ovl_hf + 113);

    auto ts_xyyyz_xyz = pbuffer.data(idx_ovl_hf + 114);

    auto ts_xyyyz_xzz = pbuffer.data(idx_ovl_hf + 115);

    auto ts_xyyyz_yyy = pbuffer.data(idx_ovl_hf + 116);

    auto ts_xyyyz_yyz = pbuffer.data(idx_ovl_hf + 117);

    auto ts_xyyyz_yzz = pbuffer.data(idx_ovl_hf + 118);

    auto ts_xyyyz_zzz = pbuffer.data(idx_ovl_hf + 119);

    auto ts_xyyzz_xxx = pbuffer.data(idx_ovl_hf + 120);

    auto ts_xyyzz_xxy = pbuffer.data(idx_ovl_hf + 121);

    auto ts_xyyzz_xxz = pbuffer.data(idx_ovl_hf + 122);

    auto ts_xyyzz_xyy = pbuffer.data(idx_ovl_hf + 123);

    auto ts_xyyzz_xyz = pbuffer.data(idx_ovl_hf + 124);

    auto ts_xyyzz_xzz = pbuffer.data(idx_ovl_hf + 125);

    auto ts_xyyzz_yyy = pbuffer.data(idx_ovl_hf + 126);

    auto ts_xyyzz_yyz = pbuffer.data(idx_ovl_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_ovl_hf + 128);

    auto ts_xyyzz_zzz = pbuffer.data(idx_ovl_hf + 129);

    auto ts_xyzzz_xxx = pbuffer.data(idx_ovl_hf + 130);

    auto ts_xyzzz_xxy = pbuffer.data(idx_ovl_hf + 131);

    auto ts_xyzzz_xxz = pbuffer.data(idx_ovl_hf + 132);

    auto ts_xyzzz_xyy = pbuffer.data(idx_ovl_hf + 133);

    auto ts_xyzzz_xyz = pbuffer.data(idx_ovl_hf + 134);

    auto ts_xyzzz_xzz = pbuffer.data(idx_ovl_hf + 135);

    auto ts_xyzzz_yyy = pbuffer.data(idx_ovl_hf + 136);

    auto ts_xyzzz_yyz = pbuffer.data(idx_ovl_hf + 137);

    auto ts_xyzzz_yzz = pbuffer.data(idx_ovl_hf + 138);

    auto ts_xyzzz_zzz = pbuffer.data(idx_ovl_hf + 139);

    auto ts_xzzzz_xxx = pbuffer.data(idx_ovl_hf + 140);

    auto ts_xzzzz_xxy = pbuffer.data(idx_ovl_hf + 141);

    auto ts_xzzzz_xxz = pbuffer.data(idx_ovl_hf + 142);

    auto ts_xzzzz_xyy = pbuffer.data(idx_ovl_hf + 143);

    auto ts_xzzzz_xyz = pbuffer.data(idx_ovl_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_ovl_hf + 145);

    auto ts_xzzzz_yyy = pbuffer.data(idx_ovl_hf + 146);

    auto ts_xzzzz_yyz = pbuffer.data(idx_ovl_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_ovl_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_ovl_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_ovl_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_ovl_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_ovl_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_ovl_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_ovl_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_ovl_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_ovl_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_ovl_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_ovl_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_ovl_hf + 159);

    auto ts_yyyyz_xxx = pbuffer.data(idx_ovl_hf + 160);

    auto ts_yyyyz_xxy = pbuffer.data(idx_ovl_hf + 161);

    auto ts_yyyyz_xxz = pbuffer.data(idx_ovl_hf + 162);

    auto ts_yyyyz_xyy = pbuffer.data(idx_ovl_hf + 163);

    auto ts_yyyyz_xyz = pbuffer.data(idx_ovl_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_ovl_hf + 165);

    auto ts_yyyyz_yyy = pbuffer.data(idx_ovl_hf + 166);

    auto ts_yyyyz_yyz = pbuffer.data(idx_ovl_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_ovl_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_ovl_hf + 169);

    auto ts_yyyzz_xxx = pbuffer.data(idx_ovl_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_ovl_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_ovl_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_ovl_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_ovl_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_ovl_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_ovl_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_ovl_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_ovl_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_ovl_hf + 179);

    auto ts_yyzzz_xxx = pbuffer.data(idx_ovl_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_ovl_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_ovl_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_ovl_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_ovl_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_ovl_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_ovl_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_ovl_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_ovl_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_ovl_hf + 189);

    auto ts_yzzzz_xxx = pbuffer.data(idx_ovl_hf + 190);

    auto ts_yzzzz_xxy = pbuffer.data(idx_ovl_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_ovl_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_ovl_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_ovl_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_ovl_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_ovl_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_ovl_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_ovl_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_ovl_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_ovl_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_ovl_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_ovl_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_ovl_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_ovl_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_ovl_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_ovl_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_ovl_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_ovl_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_ovl_hf + 209);

    // Set up 0-10 components of targeted buffer : HF

    auto tk_xxxxx_xxx = pbuffer.data(idx_kin_hf);

    auto tk_xxxxx_xxy = pbuffer.data(idx_kin_hf + 1);

    auto tk_xxxxx_xxz = pbuffer.data(idx_kin_hf + 2);

    auto tk_xxxxx_xyy = pbuffer.data(idx_kin_hf + 3);

    auto tk_xxxxx_xyz = pbuffer.data(idx_kin_hf + 4);

    auto tk_xxxxx_xzz = pbuffer.data(idx_kin_hf + 5);

    auto tk_xxxxx_yyy = pbuffer.data(idx_kin_hf + 6);

    auto tk_xxxxx_yyz = pbuffer.data(idx_kin_hf + 7);

    auto tk_xxxxx_yzz = pbuffer.data(idx_kin_hf + 8);

    auto tk_xxxxx_zzz = pbuffer.data(idx_kin_hf + 9);

    #pragma omp simd aligned(pa_x, tk_xxx_xxx, tk_xxx_xxy, tk_xxx_xxz, tk_xxx_xyy, tk_xxx_xyz, tk_xxx_xzz, tk_xxx_yyy, tk_xxx_yyz, tk_xxx_yzz, tk_xxx_zzz, tk_xxxx_xx, tk_xxxx_xxx, tk_xxxx_xxy, tk_xxxx_xxz, tk_xxxx_xy, tk_xxxx_xyy, tk_xxxx_xyz, tk_xxxx_xz, tk_xxxx_xzz, tk_xxxx_yy, tk_xxxx_yyy, tk_xxxx_yyz, tk_xxxx_yz, tk_xxxx_yzz, tk_xxxx_zz, tk_xxxx_zzz, tk_xxxxx_xxx, tk_xxxxx_xxy, tk_xxxxx_xxz, tk_xxxxx_xyy, tk_xxxxx_xyz, tk_xxxxx_xzz, tk_xxxxx_yyy, tk_xxxxx_yyz, tk_xxxxx_yzz, tk_xxxxx_zzz, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xzz, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yzz, ts_xxx_zzz, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xzz, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yzz, ts_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_xxx[i] = -8.0 * ts_xxx_xxx[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxx[i] * fe_0 + 3.0 * tk_xxxx_xx[i] * fe_0 + tk_xxxx_xxx[i] * pa_x[i] + 2.0 * ts_xxxxx_xxx[i] * fz_0;

        tk_xxxxx_xxy[i] = -8.0 * ts_xxx_xxy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxy[i] * fe_0 + 2.0 * tk_xxxx_xy[i] * fe_0 + tk_xxxx_xxy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxy[i] * fz_0;

        tk_xxxxx_xxz[i] = -8.0 * ts_xxx_xxz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxz[i] * fe_0 + 2.0 * tk_xxxx_xz[i] * fe_0 + tk_xxxx_xxz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxz[i] * fz_0;

        tk_xxxxx_xyy[i] = -8.0 * ts_xxx_xyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyy[i] * fe_0 + tk_xxxx_yy[i] * fe_0 + tk_xxxx_xyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xyy[i] * fz_0;

        tk_xxxxx_xyz[i] = -8.0 * ts_xxx_xyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyz[i] * fe_0 + tk_xxxx_yz[i] * fe_0 + tk_xxxx_xyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyz[i] * fz_0;

        tk_xxxxx_xzz[i] = -8.0 * ts_xxx_xzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xzz[i] * fe_0 + tk_xxxx_zz[i] * fe_0 + tk_xxxx_xzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xzz[i] * fz_0;

        tk_xxxxx_yyy[i] = -8.0 * ts_xxx_yyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyy[i] * fe_0 + tk_xxxx_yyy[i] * pa_x[i] + 2.0 * ts_xxxxx_yyy[i] * fz_0;

        tk_xxxxx_yyz[i] = -8.0 * ts_xxx_yyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyz[i] * fe_0 + tk_xxxx_yyz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyz[i] * fz_0;

        tk_xxxxx_yzz[i] = -8.0 * ts_xxx_yzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yzz[i] * fe_0 + tk_xxxx_yzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yzz[i] * fz_0;

        tk_xxxxx_zzz[i] = -8.0 * ts_xxx_zzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_zzz[i] * fe_0 + tk_xxxx_zzz[i] * pa_x[i] + 2.0 * ts_xxxxx_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : HF

    auto tk_xxxxy_xxx = pbuffer.data(idx_kin_hf + 10);

    auto tk_xxxxy_xxy = pbuffer.data(idx_kin_hf + 11);

    auto tk_xxxxy_xxz = pbuffer.data(idx_kin_hf + 12);

    auto tk_xxxxy_xyy = pbuffer.data(idx_kin_hf + 13);

    auto tk_xxxxy_xyz = pbuffer.data(idx_kin_hf + 14);

    auto tk_xxxxy_xzz = pbuffer.data(idx_kin_hf + 15);

    auto tk_xxxxy_yyy = pbuffer.data(idx_kin_hf + 16);

    auto tk_xxxxy_yyz = pbuffer.data(idx_kin_hf + 17);

    auto tk_xxxxy_yzz = pbuffer.data(idx_kin_hf + 18);

    auto tk_xxxxy_zzz = pbuffer.data(idx_kin_hf + 19);

    #pragma omp simd aligned(pa_y, tk_xxxx_xx, tk_xxxx_xxx, tk_xxxx_xxy, tk_xxxx_xxz, tk_xxxx_xy, tk_xxxx_xyy, tk_xxxx_xyz, tk_xxxx_xz, tk_xxxx_xzz, tk_xxxx_yy, tk_xxxx_yyy, tk_xxxx_yyz, tk_xxxx_yz, tk_xxxx_yzz, tk_xxxx_zz, tk_xxxx_zzz, tk_xxxxy_xxx, tk_xxxxy_xxy, tk_xxxxy_xxz, tk_xxxxy_xyy, tk_xxxxy_xyz, tk_xxxxy_xzz, tk_xxxxy_yyy, tk_xxxxy_yyz, tk_xxxxy_yzz, tk_xxxxy_zzz, ts_xxxxy_xxx, ts_xxxxy_xxy, ts_xxxxy_xxz, ts_xxxxy_xyy, ts_xxxxy_xyz, ts_xxxxy_xzz, ts_xxxxy_yyy, ts_xxxxy_yyz, ts_xxxxy_yzz, ts_xxxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_xxx[i] = tk_xxxx_xxx[i] * pa_y[i] + 2.0 * ts_xxxxy_xxx[i] * fz_0;

        tk_xxxxy_xxy[i] = tk_xxxx_xx[i] * fe_0 + tk_xxxx_xxy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxy[i] * fz_0;

        tk_xxxxy_xxz[i] = tk_xxxx_xxz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxz[i] * fz_0;

        tk_xxxxy_xyy[i] = 2.0 * tk_xxxx_xy[i] * fe_0 + tk_xxxx_xyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xyy[i] * fz_0;

        tk_xxxxy_xyz[i] = tk_xxxx_xz[i] * fe_0 + tk_xxxx_xyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyz[i] * fz_0;

        tk_xxxxy_xzz[i] = tk_xxxx_xzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xzz[i] * fz_0;

        tk_xxxxy_yyy[i] = 3.0 * tk_xxxx_yy[i] * fe_0 + tk_xxxx_yyy[i] * pa_y[i] + 2.0 * ts_xxxxy_yyy[i] * fz_0;

        tk_xxxxy_yyz[i] = 2.0 * tk_xxxx_yz[i] * fe_0 + tk_xxxx_yyz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyz[i] * fz_0;

        tk_xxxxy_yzz[i] = tk_xxxx_zz[i] * fe_0 + tk_xxxx_yzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yzz[i] * fz_0;

        tk_xxxxy_zzz[i] = tk_xxxx_zzz[i] * pa_y[i] + 2.0 * ts_xxxxy_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : HF

    auto tk_xxxxz_xxx = pbuffer.data(idx_kin_hf + 20);

    auto tk_xxxxz_xxy = pbuffer.data(idx_kin_hf + 21);

    auto tk_xxxxz_xxz = pbuffer.data(idx_kin_hf + 22);

    auto tk_xxxxz_xyy = pbuffer.data(idx_kin_hf + 23);

    auto tk_xxxxz_xyz = pbuffer.data(idx_kin_hf + 24);

    auto tk_xxxxz_xzz = pbuffer.data(idx_kin_hf + 25);

    auto tk_xxxxz_yyy = pbuffer.data(idx_kin_hf + 26);

    auto tk_xxxxz_yyz = pbuffer.data(idx_kin_hf + 27);

    auto tk_xxxxz_yzz = pbuffer.data(idx_kin_hf + 28);

    auto tk_xxxxz_zzz = pbuffer.data(idx_kin_hf + 29);

    #pragma omp simd aligned(pa_z, tk_xxxx_xx, tk_xxxx_xxx, tk_xxxx_xxy, tk_xxxx_xxz, tk_xxxx_xy, tk_xxxx_xyy, tk_xxxx_xyz, tk_xxxx_xz, tk_xxxx_xzz, tk_xxxx_yy, tk_xxxx_yyy, tk_xxxx_yyz, tk_xxxx_yz, tk_xxxx_yzz, tk_xxxx_zz, tk_xxxx_zzz, tk_xxxxz_xxx, tk_xxxxz_xxy, tk_xxxxz_xxz, tk_xxxxz_xyy, tk_xxxxz_xyz, tk_xxxxz_xzz, tk_xxxxz_yyy, tk_xxxxz_yyz, tk_xxxxz_yzz, tk_xxxxz_zzz, ts_xxxxz_xxx, ts_xxxxz_xxy, ts_xxxxz_xxz, ts_xxxxz_xyy, ts_xxxxz_xyz, ts_xxxxz_xzz, ts_xxxxz_yyy, ts_xxxxz_yyz, ts_xxxxz_yzz, ts_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_xxx[i] = tk_xxxx_xxx[i] * pa_z[i] + 2.0 * ts_xxxxz_xxx[i] * fz_0;

        tk_xxxxz_xxy[i] = tk_xxxx_xxy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxy[i] * fz_0;

        tk_xxxxz_xxz[i] = tk_xxxx_xx[i] * fe_0 + tk_xxxx_xxz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxz[i] * fz_0;

        tk_xxxxz_xyy[i] = tk_xxxx_xyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xyy[i] * fz_0;

        tk_xxxxz_xyz[i] = tk_xxxx_xy[i] * fe_0 + tk_xxxx_xyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyz[i] * fz_0;

        tk_xxxxz_xzz[i] = 2.0 * tk_xxxx_xz[i] * fe_0 + tk_xxxx_xzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xzz[i] * fz_0;

        tk_xxxxz_yyy[i] = tk_xxxx_yyy[i] * pa_z[i] + 2.0 * ts_xxxxz_yyy[i] * fz_0;

        tk_xxxxz_yyz[i] = tk_xxxx_yy[i] * fe_0 + tk_xxxx_yyz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyz[i] * fz_0;

        tk_xxxxz_yzz[i] = 2.0 * tk_xxxx_yz[i] * fe_0 + tk_xxxx_yzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yzz[i] * fz_0;

        tk_xxxxz_zzz[i] = 3.0 * tk_xxxx_zz[i] * fe_0 + tk_xxxx_zzz[i] * pa_z[i] + 2.0 * ts_xxxxz_zzz[i] * fz_0;
    }

    // Set up 30-40 components of targeted buffer : HF

    auto tk_xxxyy_xxx = pbuffer.data(idx_kin_hf + 30);

    auto tk_xxxyy_xxy = pbuffer.data(idx_kin_hf + 31);

    auto tk_xxxyy_xxz = pbuffer.data(idx_kin_hf + 32);

    auto tk_xxxyy_xyy = pbuffer.data(idx_kin_hf + 33);

    auto tk_xxxyy_xyz = pbuffer.data(idx_kin_hf + 34);

    auto tk_xxxyy_xzz = pbuffer.data(idx_kin_hf + 35);

    auto tk_xxxyy_yyy = pbuffer.data(idx_kin_hf + 36);

    auto tk_xxxyy_yyz = pbuffer.data(idx_kin_hf + 37);

    auto tk_xxxyy_yzz = pbuffer.data(idx_kin_hf + 38);

    auto tk_xxxyy_zzz = pbuffer.data(idx_kin_hf + 39);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxx_xxx, tk_xxx_xxz, tk_xxx_xzz, tk_xxxy_xxx, tk_xxxy_xxz, tk_xxxy_xzz, tk_xxxyy_xxx, tk_xxxyy_xxy, tk_xxxyy_xxz, tk_xxxyy_xyy, tk_xxxyy_xyz, tk_xxxyy_xzz, tk_xxxyy_yyy, tk_xxxyy_yyz, tk_xxxyy_yzz, tk_xxxyy_zzz, tk_xxyy_xxy, tk_xxyy_xy, tk_xxyy_xyy, tk_xxyy_xyz, tk_xxyy_yy, tk_xxyy_yyy, tk_xxyy_yyz, tk_xxyy_yz, tk_xxyy_yzz, tk_xxyy_zzz, tk_xyy_xxy, tk_xyy_xyy, tk_xyy_xyz, tk_xyy_yyy, tk_xyy_yyz, tk_xyy_yzz, tk_xyy_zzz, ts_xxx_xxx, ts_xxx_xxz, ts_xxx_xzz, ts_xxxyy_xxx, ts_xxxyy_xxy, ts_xxxyy_xxz, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_xzz, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yzz, ts_xxxyy_zzz, ts_xyy_xxy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yzz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_xxx[i] = -2.0 * ts_xxx_xxx[i] * fbe_0 * fz_0 + tk_xxx_xxx[i] * fe_0 + tk_xxxy_xxx[i] * pa_y[i] + 2.0 * ts_xxxyy_xxx[i] * fz_0;

        tk_xxxyy_xxy[i] = -4.0 * ts_xyy_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxy[i] * fe_0 + 2.0 * tk_xxyy_xy[i] * fe_0 + tk_xxyy_xxy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxy[i] * fz_0;

        tk_xxxyy_xxz[i] = -2.0 * ts_xxx_xxz[i] * fbe_0 * fz_0 + tk_xxx_xxz[i] * fe_0 + tk_xxxy_xxz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxz[i] * fz_0;

        tk_xxxyy_xyy[i] = -4.0 * ts_xyy_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyy[i] * fe_0 + tk_xxyy_yy[i] * fe_0 + tk_xxyy_xyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xyy[i] * fz_0;

        tk_xxxyy_xyz[i] = -4.0 * ts_xyy_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyz[i] * fe_0 + tk_xxyy_yz[i] * fe_0 + tk_xxyy_xyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyz[i] * fz_0;

        tk_xxxyy_xzz[i] = -2.0 * ts_xxx_xzz[i] * fbe_0 * fz_0 + tk_xxx_xzz[i] * fe_0 + tk_xxxy_xzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xzz[i] * fz_0;

        tk_xxxyy_yyy[i] = -4.0 * ts_xyy_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyy[i] * fe_0 + tk_xxyy_yyy[i] * pa_x[i] + 2.0 * ts_xxxyy_yyy[i] * fz_0;

        tk_xxxyy_yyz[i] = -4.0 * ts_xyy_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyz[i] * fe_0 + tk_xxyy_yyz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyz[i] * fz_0;

        tk_xxxyy_yzz[i] = -4.0 * ts_xyy_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yzz[i] * fe_0 + tk_xxyy_yzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yzz[i] * fz_0;

        tk_xxxyy_zzz[i] = -4.0 * ts_xyy_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_zzz[i] * fe_0 + tk_xxyy_zzz[i] * pa_x[i] + 2.0 * ts_xxxyy_zzz[i] * fz_0;
    }

    // Set up 40-50 components of targeted buffer : HF

    auto tk_xxxyz_xxx = pbuffer.data(idx_kin_hf + 40);

    auto tk_xxxyz_xxy = pbuffer.data(idx_kin_hf + 41);

    auto tk_xxxyz_xxz = pbuffer.data(idx_kin_hf + 42);

    auto tk_xxxyz_xyy = pbuffer.data(idx_kin_hf + 43);

    auto tk_xxxyz_xyz = pbuffer.data(idx_kin_hf + 44);

    auto tk_xxxyz_xzz = pbuffer.data(idx_kin_hf + 45);

    auto tk_xxxyz_yyy = pbuffer.data(idx_kin_hf + 46);

    auto tk_xxxyz_yyz = pbuffer.data(idx_kin_hf + 47);

    auto tk_xxxyz_yzz = pbuffer.data(idx_kin_hf + 48);

    auto tk_xxxyz_zzz = pbuffer.data(idx_kin_hf + 49);

    #pragma omp simd aligned(pa_y, pa_z, tk_xxxy_xxy, tk_xxxy_xyy, tk_xxxy_yyy, tk_xxxyz_xxx, tk_xxxyz_xxy, tk_xxxyz_xxz, tk_xxxyz_xyy, tk_xxxyz_xyz, tk_xxxyz_xzz, tk_xxxyz_yyy, tk_xxxyz_yyz, tk_xxxyz_yzz, tk_xxxyz_zzz, tk_xxxz_xxx, tk_xxxz_xxz, tk_xxxz_xyz, tk_xxxz_xz, tk_xxxz_xzz, tk_xxxz_yyz, tk_xxxz_yz, tk_xxxz_yzz, tk_xxxz_zz, tk_xxxz_zzz, ts_xxxyz_xxx, ts_xxxyz_xxy, ts_xxxyz_xxz, ts_xxxyz_xyy, ts_xxxyz_xyz, ts_xxxyz_xzz, ts_xxxyz_yyy, ts_xxxyz_yyz, ts_xxxyz_yzz, ts_xxxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyz_xxx[i] = tk_xxxz_xxx[i] * pa_y[i] + 2.0 * ts_xxxyz_xxx[i] * fz_0;

        tk_xxxyz_xxy[i] = tk_xxxy_xxy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxy[i] * fz_0;

        tk_xxxyz_xxz[i] = tk_xxxz_xxz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxz[i] * fz_0;

        tk_xxxyz_xyy[i] = tk_xxxy_xyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xyy[i] * fz_0;

        tk_xxxyz_xyz[i] = tk_xxxz_xz[i] * fe_0 + tk_xxxz_xyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyz[i] * fz_0;

        tk_xxxyz_xzz[i] = tk_xxxz_xzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xzz[i] * fz_0;

        tk_xxxyz_yyy[i] = tk_xxxy_yyy[i] * pa_z[i] + 2.0 * ts_xxxyz_yyy[i] * fz_0;

        tk_xxxyz_yyz[i] = 2.0 * tk_xxxz_yz[i] * fe_0 + tk_xxxz_yyz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyz[i] * fz_0;

        tk_xxxyz_yzz[i] = tk_xxxz_zz[i] * fe_0 + tk_xxxz_yzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yzz[i] * fz_0;

        tk_xxxyz_zzz[i] = tk_xxxz_zzz[i] * pa_y[i] + 2.0 * ts_xxxyz_zzz[i] * fz_0;
    }

    // Set up 50-60 components of targeted buffer : HF

    auto tk_xxxzz_xxx = pbuffer.data(idx_kin_hf + 50);

    auto tk_xxxzz_xxy = pbuffer.data(idx_kin_hf + 51);

    auto tk_xxxzz_xxz = pbuffer.data(idx_kin_hf + 52);

    auto tk_xxxzz_xyy = pbuffer.data(idx_kin_hf + 53);

    auto tk_xxxzz_xyz = pbuffer.data(idx_kin_hf + 54);

    auto tk_xxxzz_xzz = pbuffer.data(idx_kin_hf + 55);

    auto tk_xxxzz_yyy = pbuffer.data(idx_kin_hf + 56);

    auto tk_xxxzz_yyz = pbuffer.data(idx_kin_hf + 57);

    auto tk_xxxzz_yzz = pbuffer.data(idx_kin_hf + 58);

    auto tk_xxxzz_zzz = pbuffer.data(idx_kin_hf + 59);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxx_xxx, tk_xxx_xxy, tk_xxx_xyy, tk_xxxz_xxx, tk_xxxz_xxy, tk_xxxz_xyy, tk_xxxzz_xxx, tk_xxxzz_xxy, tk_xxxzz_xxz, tk_xxxzz_xyy, tk_xxxzz_xyz, tk_xxxzz_xzz, tk_xxxzz_yyy, tk_xxxzz_yyz, tk_xxxzz_yzz, tk_xxxzz_zzz, tk_xxzz_xxz, tk_xxzz_xyz, tk_xxzz_xz, tk_xxzz_xzz, tk_xxzz_yyy, tk_xxzz_yyz, tk_xxzz_yz, tk_xxzz_yzz, tk_xxzz_zz, tk_xxzz_zzz, tk_xzz_xxz, tk_xzz_xyz, tk_xzz_xzz, tk_xzz_yyy, tk_xzz_yyz, tk_xzz_yzz, tk_xzz_zzz, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xyy, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xxz, ts_xxxzz_xyy, ts_xxxzz_xyz, ts_xxxzz_xzz, ts_xxxzz_yyy, ts_xxxzz_yyz, ts_xxxzz_yzz, ts_xxxzz_zzz, ts_xzz_xxz, ts_xzz_xyz, ts_xzz_xzz, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yzz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_xxx[i] = -2.0 * ts_xxx_xxx[i] * fbe_0 * fz_0 + tk_xxx_xxx[i] * fe_0 + tk_xxxz_xxx[i] * pa_z[i] + 2.0 * ts_xxxzz_xxx[i] * fz_0;

        tk_xxxzz_xxy[i] = -2.0 * ts_xxx_xxy[i] * fbe_0 * fz_0 + tk_xxx_xxy[i] * fe_0 + tk_xxxz_xxy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxy[i] * fz_0;

        tk_xxxzz_xxz[i] = -4.0 * ts_xzz_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxz[i] * fe_0 + 2.0 * tk_xxzz_xz[i] * fe_0 + tk_xxzz_xxz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxz[i] * fz_0;

        tk_xxxzz_xyy[i] = -2.0 * ts_xxx_xyy[i] * fbe_0 * fz_0 + tk_xxx_xyy[i] * fe_0 + tk_xxxz_xyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xyy[i] * fz_0;

        tk_xxxzz_xyz[i] = -4.0 * ts_xzz_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyz[i] * fe_0 + tk_xxzz_yz[i] * fe_0 + tk_xxzz_xyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyz[i] * fz_0;

        tk_xxxzz_xzz[i] = -4.0 * ts_xzz_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xzz[i] * fe_0 + tk_xxzz_zz[i] * fe_0 + tk_xxzz_xzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xzz[i] * fz_0;

        tk_xxxzz_yyy[i] = -4.0 * ts_xzz_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyy[i] * fe_0 + tk_xxzz_yyy[i] * pa_x[i] + 2.0 * ts_xxxzz_yyy[i] * fz_0;

        tk_xxxzz_yyz[i] = -4.0 * ts_xzz_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyz[i] * fe_0 + tk_xxzz_yyz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyz[i] * fz_0;

        tk_xxxzz_yzz[i] = -4.0 * ts_xzz_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yzz[i] * fe_0 + tk_xxzz_yzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yzz[i] * fz_0;

        tk_xxxzz_zzz[i] = -4.0 * ts_xzz_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_zzz[i] * fe_0 + tk_xxzz_zzz[i] * pa_x[i] + 2.0 * ts_xxxzz_zzz[i] * fz_0;
    }

    // Set up 60-70 components of targeted buffer : HF

    auto tk_xxyyy_xxx = pbuffer.data(idx_kin_hf + 60);

    auto tk_xxyyy_xxy = pbuffer.data(idx_kin_hf + 61);

    auto tk_xxyyy_xxz = pbuffer.data(idx_kin_hf + 62);

    auto tk_xxyyy_xyy = pbuffer.data(idx_kin_hf + 63);

    auto tk_xxyyy_xyz = pbuffer.data(idx_kin_hf + 64);

    auto tk_xxyyy_xzz = pbuffer.data(idx_kin_hf + 65);

    auto tk_xxyyy_yyy = pbuffer.data(idx_kin_hf + 66);

    auto tk_xxyyy_yyz = pbuffer.data(idx_kin_hf + 67);

    auto tk_xxyyy_yzz = pbuffer.data(idx_kin_hf + 68);

    auto tk_xxyyy_zzz = pbuffer.data(idx_kin_hf + 69);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxy_xxx, tk_xxy_xxz, tk_xxy_xzz, tk_xxyy_xxx, tk_xxyy_xxz, tk_xxyy_xzz, tk_xxyyy_xxx, tk_xxyyy_xxy, tk_xxyyy_xxz, tk_xxyyy_xyy, tk_xxyyy_xyz, tk_xxyyy_xzz, tk_xxyyy_yyy, tk_xxyyy_yyz, tk_xxyyy_yzz, tk_xxyyy_zzz, tk_xyyy_xxy, tk_xyyy_xy, tk_xyyy_xyy, tk_xyyy_xyz, tk_xyyy_yy, tk_xyyy_yyy, tk_xyyy_yyz, tk_xyyy_yz, tk_xyyy_yzz, tk_xyyy_zzz, tk_yyy_xxy, tk_yyy_xyy, tk_yyy_xyz, tk_yyy_yyy, tk_yyy_yyz, tk_yyy_yzz, tk_yyy_zzz, ts_xxy_xxx, ts_xxy_xxz, ts_xxy_xzz, ts_xxyyy_xxx, ts_xxyyy_xxy, ts_xxyyy_xxz, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_xzz, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yzz, ts_xxyyy_zzz, ts_yyy_xxy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yzz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_xxx[i] = -4.0 * ts_xxy_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxx[i] * fe_0 + tk_xxyy_xxx[i] * pa_y[i] + 2.0 * ts_xxyyy_xxx[i] * fz_0;

        tk_xxyyy_xxy[i] = -2.0 * ts_yyy_xxy[i] * fbe_0 * fz_0 + tk_yyy_xxy[i] * fe_0 + 2.0 * tk_xyyy_xy[i] * fe_0 + tk_xyyy_xxy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxy[i] * fz_0;

        tk_xxyyy_xxz[i] = -4.0 * ts_xxy_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxz[i] * fe_0 + tk_xxyy_xxz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxz[i] * fz_0;

        tk_xxyyy_xyy[i] = -2.0 * ts_yyy_xyy[i] * fbe_0 * fz_0 + tk_yyy_xyy[i] * fe_0 + tk_xyyy_yy[i] * fe_0 + tk_xyyy_xyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xyy[i] * fz_0;

        tk_xxyyy_xyz[i] = -2.0 * ts_yyy_xyz[i] * fbe_0 * fz_0 + tk_yyy_xyz[i] * fe_0 + tk_xyyy_yz[i] * fe_0 + tk_xyyy_xyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xyz[i] * fz_0;

        tk_xxyyy_xzz[i] = -4.0 * ts_xxy_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xzz[i] * fe_0 + tk_xxyy_xzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xzz[i] * fz_0;

        tk_xxyyy_yyy[i] = -2.0 * ts_yyy_yyy[i] * fbe_0 * fz_0 + tk_yyy_yyy[i] * fe_0 + tk_xyyy_yyy[i] * pa_x[i] + 2.0 * ts_xxyyy_yyy[i] * fz_0;

        tk_xxyyy_yyz[i] = -2.0 * ts_yyy_yyz[i] * fbe_0 * fz_0 + tk_yyy_yyz[i] * fe_0 + tk_xyyy_yyz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyz[i] * fz_0;

        tk_xxyyy_yzz[i] = -2.0 * ts_yyy_yzz[i] * fbe_0 * fz_0 + tk_yyy_yzz[i] * fe_0 + tk_xyyy_yzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yzz[i] * fz_0;

        tk_xxyyy_zzz[i] = -2.0 * ts_yyy_zzz[i] * fbe_0 * fz_0 + tk_yyy_zzz[i] * fe_0 + tk_xyyy_zzz[i] * pa_x[i] + 2.0 * ts_xxyyy_zzz[i] * fz_0;
    }

    // Set up 70-80 components of targeted buffer : HF

    auto tk_xxyyz_xxx = pbuffer.data(idx_kin_hf + 70);

    auto tk_xxyyz_xxy = pbuffer.data(idx_kin_hf + 71);

    auto tk_xxyyz_xxz = pbuffer.data(idx_kin_hf + 72);

    auto tk_xxyyz_xyy = pbuffer.data(idx_kin_hf + 73);

    auto tk_xxyyz_xyz = pbuffer.data(idx_kin_hf + 74);

    auto tk_xxyyz_xzz = pbuffer.data(idx_kin_hf + 75);

    auto tk_xxyyz_yyy = pbuffer.data(idx_kin_hf + 76);

    auto tk_xxyyz_yyz = pbuffer.data(idx_kin_hf + 77);

    auto tk_xxyyz_yzz = pbuffer.data(idx_kin_hf + 78);

    auto tk_xxyyz_zzz = pbuffer.data(idx_kin_hf + 79);

    #pragma omp simd aligned(pa_z, tk_xxyy_xx, tk_xxyy_xxx, tk_xxyy_xxy, tk_xxyy_xxz, tk_xxyy_xy, tk_xxyy_xyy, tk_xxyy_xyz, tk_xxyy_xz, tk_xxyy_xzz, tk_xxyy_yy, tk_xxyy_yyy, tk_xxyy_yyz, tk_xxyy_yz, tk_xxyy_yzz, tk_xxyy_zz, tk_xxyy_zzz, tk_xxyyz_xxx, tk_xxyyz_xxy, tk_xxyyz_xxz, tk_xxyyz_xyy, tk_xxyyz_xyz, tk_xxyyz_xzz, tk_xxyyz_yyy, tk_xxyyz_yyz, tk_xxyyz_yzz, tk_xxyyz_zzz, ts_xxyyz_xxx, ts_xxyyz_xxy, ts_xxyyz_xxz, ts_xxyyz_xyy, ts_xxyyz_xyz, ts_xxyyz_xzz, ts_xxyyz_yyy, ts_xxyyz_yyz, ts_xxyyz_yzz, ts_xxyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_xxx[i] = tk_xxyy_xxx[i] * pa_z[i] + 2.0 * ts_xxyyz_xxx[i] * fz_0;

        tk_xxyyz_xxy[i] = tk_xxyy_xxy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxy[i] * fz_0;

        tk_xxyyz_xxz[i] = tk_xxyy_xx[i] * fe_0 + tk_xxyy_xxz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxz[i] * fz_0;

        tk_xxyyz_xyy[i] = tk_xxyy_xyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xyy[i] * fz_0;

        tk_xxyyz_xyz[i] = tk_xxyy_xy[i] * fe_0 + tk_xxyy_xyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyz[i] * fz_0;

        tk_xxyyz_xzz[i] = 2.0 * tk_xxyy_xz[i] * fe_0 + tk_xxyy_xzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xzz[i] * fz_0;

        tk_xxyyz_yyy[i] = tk_xxyy_yyy[i] * pa_z[i] + 2.0 * ts_xxyyz_yyy[i] * fz_0;

        tk_xxyyz_yyz[i] = tk_xxyy_yy[i] * fe_0 + tk_xxyy_yyz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyz[i] * fz_0;

        tk_xxyyz_yzz[i] = 2.0 * tk_xxyy_yz[i] * fe_0 + tk_xxyy_yzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yzz[i] * fz_0;

        tk_xxyyz_zzz[i] = 3.0 * tk_xxyy_zz[i] * fe_0 + tk_xxyy_zzz[i] * pa_z[i] + 2.0 * ts_xxyyz_zzz[i] * fz_0;
    }

    // Set up 80-90 components of targeted buffer : HF

    auto tk_xxyzz_xxx = pbuffer.data(idx_kin_hf + 80);

    auto tk_xxyzz_xxy = pbuffer.data(idx_kin_hf + 81);

    auto tk_xxyzz_xxz = pbuffer.data(idx_kin_hf + 82);

    auto tk_xxyzz_xyy = pbuffer.data(idx_kin_hf + 83);

    auto tk_xxyzz_xyz = pbuffer.data(idx_kin_hf + 84);

    auto tk_xxyzz_xzz = pbuffer.data(idx_kin_hf + 85);

    auto tk_xxyzz_yyy = pbuffer.data(idx_kin_hf + 86);

    auto tk_xxyzz_yyz = pbuffer.data(idx_kin_hf + 87);

    auto tk_xxyzz_yzz = pbuffer.data(idx_kin_hf + 88);

    auto tk_xxyzz_zzz = pbuffer.data(idx_kin_hf + 89);

    #pragma omp simd aligned(pa_y, tk_xxyzz_xxx, tk_xxyzz_xxy, tk_xxyzz_xxz, tk_xxyzz_xyy, tk_xxyzz_xyz, tk_xxyzz_xzz, tk_xxyzz_yyy, tk_xxyzz_yyz, tk_xxyzz_yzz, tk_xxyzz_zzz, tk_xxzz_xx, tk_xxzz_xxx, tk_xxzz_xxy, tk_xxzz_xxz, tk_xxzz_xy, tk_xxzz_xyy, tk_xxzz_xyz, tk_xxzz_xz, tk_xxzz_xzz, tk_xxzz_yy, tk_xxzz_yyy, tk_xxzz_yyz, tk_xxzz_yz, tk_xxzz_yzz, tk_xxzz_zz, tk_xxzz_zzz, ts_xxyzz_xxx, ts_xxyzz_xxy, ts_xxyzz_xxz, ts_xxyzz_xyy, ts_xxyzz_xyz, ts_xxyzz_xzz, ts_xxyzz_yyy, ts_xxyzz_yyz, ts_xxyzz_yzz, ts_xxyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_xxx[i] = tk_xxzz_xxx[i] * pa_y[i] + 2.0 * ts_xxyzz_xxx[i] * fz_0;

        tk_xxyzz_xxy[i] = tk_xxzz_xx[i] * fe_0 + tk_xxzz_xxy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxy[i] * fz_0;

        tk_xxyzz_xxz[i] = tk_xxzz_xxz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxz[i] * fz_0;

        tk_xxyzz_xyy[i] = 2.0 * tk_xxzz_xy[i] * fe_0 + tk_xxzz_xyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xyy[i] * fz_0;

        tk_xxyzz_xyz[i] = tk_xxzz_xz[i] * fe_0 + tk_xxzz_xyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyz[i] * fz_0;

        tk_xxyzz_xzz[i] = tk_xxzz_xzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xzz[i] * fz_0;

        tk_xxyzz_yyy[i] = 3.0 * tk_xxzz_yy[i] * fe_0 + tk_xxzz_yyy[i] * pa_y[i] + 2.0 * ts_xxyzz_yyy[i] * fz_0;

        tk_xxyzz_yyz[i] = 2.0 * tk_xxzz_yz[i] * fe_0 + tk_xxzz_yyz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyz[i] * fz_0;

        tk_xxyzz_yzz[i] = tk_xxzz_zz[i] * fe_0 + tk_xxzz_yzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yzz[i] * fz_0;

        tk_xxyzz_zzz[i] = tk_xxzz_zzz[i] * pa_y[i] + 2.0 * ts_xxyzz_zzz[i] * fz_0;
    }

    // Set up 90-100 components of targeted buffer : HF

    auto tk_xxzzz_xxx = pbuffer.data(idx_kin_hf + 90);

    auto tk_xxzzz_xxy = pbuffer.data(idx_kin_hf + 91);

    auto tk_xxzzz_xxz = pbuffer.data(idx_kin_hf + 92);

    auto tk_xxzzz_xyy = pbuffer.data(idx_kin_hf + 93);

    auto tk_xxzzz_xyz = pbuffer.data(idx_kin_hf + 94);

    auto tk_xxzzz_xzz = pbuffer.data(idx_kin_hf + 95);

    auto tk_xxzzz_yyy = pbuffer.data(idx_kin_hf + 96);

    auto tk_xxzzz_yyz = pbuffer.data(idx_kin_hf + 97);

    auto tk_xxzzz_yzz = pbuffer.data(idx_kin_hf + 98);

    auto tk_xxzzz_zzz = pbuffer.data(idx_kin_hf + 99);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxz_xxx, tk_xxz_xxy, tk_xxz_xyy, tk_xxzz_xxx, tk_xxzz_xxy, tk_xxzz_xyy, tk_xxzzz_xxx, tk_xxzzz_xxy, tk_xxzzz_xxz, tk_xxzzz_xyy, tk_xxzzz_xyz, tk_xxzzz_xzz, tk_xxzzz_yyy, tk_xxzzz_yyz, tk_xxzzz_yzz, tk_xxzzz_zzz, tk_xzzz_xxz, tk_xzzz_xyz, tk_xzzz_xz, tk_xzzz_xzz, tk_xzzz_yyy, tk_xzzz_yyz, tk_xzzz_yz, tk_xzzz_yzz, tk_xzzz_zz, tk_xzzz_zzz, tk_zzz_xxz, tk_zzz_xyz, tk_zzz_xzz, tk_zzz_yyy, tk_zzz_yyz, tk_zzz_yzz, tk_zzz_zzz, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xyy, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xxz, ts_xxzzz_xyy, ts_xxzzz_xyz, ts_xxzzz_xzz, ts_xxzzz_yyy, ts_xxzzz_yyz, ts_xxzzz_yzz, ts_xxzzz_zzz, ts_zzz_xxz, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_xxx[i] = -4.0 * ts_xxz_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxx[i] * fe_0 + tk_xxzz_xxx[i] * pa_z[i] + 2.0 * ts_xxzzz_xxx[i] * fz_0;

        tk_xxzzz_xxy[i] = -4.0 * ts_xxz_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxy[i] * fe_0 + tk_xxzz_xxy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxy[i] * fz_0;

        tk_xxzzz_xxz[i] = -2.0 * ts_zzz_xxz[i] * fbe_0 * fz_0 + tk_zzz_xxz[i] * fe_0 + 2.0 * tk_xzzz_xz[i] * fe_0 + tk_xzzz_xxz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxz[i] * fz_0;

        tk_xxzzz_xyy[i] = -4.0 * ts_xxz_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xyy[i] * fe_0 + tk_xxzz_xyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xyy[i] * fz_0;

        tk_xxzzz_xyz[i] = -2.0 * ts_zzz_xyz[i] * fbe_0 * fz_0 + tk_zzz_xyz[i] * fe_0 + tk_xzzz_yz[i] * fe_0 + tk_xzzz_xyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xyz[i] * fz_0;

        tk_xxzzz_xzz[i] = -2.0 * ts_zzz_xzz[i] * fbe_0 * fz_0 + tk_zzz_xzz[i] * fe_0 + tk_xzzz_zz[i] * fe_0 + tk_xzzz_xzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xzz[i] * fz_0;

        tk_xxzzz_yyy[i] = -2.0 * ts_zzz_yyy[i] * fbe_0 * fz_0 + tk_zzz_yyy[i] * fe_0 + tk_xzzz_yyy[i] * pa_x[i] + 2.0 * ts_xxzzz_yyy[i] * fz_0;

        tk_xxzzz_yyz[i] = -2.0 * ts_zzz_yyz[i] * fbe_0 * fz_0 + tk_zzz_yyz[i] * fe_0 + tk_xzzz_yyz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyz[i] * fz_0;

        tk_xxzzz_yzz[i] = -2.0 * ts_zzz_yzz[i] * fbe_0 * fz_0 + tk_zzz_yzz[i] * fe_0 + tk_xzzz_yzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yzz[i] * fz_0;

        tk_xxzzz_zzz[i] = -2.0 * ts_zzz_zzz[i] * fbe_0 * fz_0 + tk_zzz_zzz[i] * fe_0 + tk_xzzz_zzz[i] * pa_x[i] + 2.0 * ts_xxzzz_zzz[i] * fz_0;
    }

    // Set up 100-110 components of targeted buffer : HF

    auto tk_xyyyy_xxx = pbuffer.data(idx_kin_hf + 100);

    auto tk_xyyyy_xxy = pbuffer.data(idx_kin_hf + 101);

    auto tk_xyyyy_xxz = pbuffer.data(idx_kin_hf + 102);

    auto tk_xyyyy_xyy = pbuffer.data(idx_kin_hf + 103);

    auto tk_xyyyy_xyz = pbuffer.data(idx_kin_hf + 104);

    auto tk_xyyyy_xzz = pbuffer.data(idx_kin_hf + 105);

    auto tk_xyyyy_yyy = pbuffer.data(idx_kin_hf + 106);

    auto tk_xyyyy_yyz = pbuffer.data(idx_kin_hf + 107);

    auto tk_xyyyy_yzz = pbuffer.data(idx_kin_hf + 108);

    auto tk_xyyyy_zzz = pbuffer.data(idx_kin_hf + 109);

    #pragma omp simd aligned(pa_x, tk_xyyyy_xxx, tk_xyyyy_xxy, tk_xyyyy_xxz, tk_xyyyy_xyy, tk_xyyyy_xyz, tk_xyyyy_xzz, tk_xyyyy_yyy, tk_xyyyy_yyz, tk_xyyyy_yzz, tk_xyyyy_zzz, tk_yyyy_xx, tk_yyyy_xxx, tk_yyyy_xxy, tk_yyyy_xxz, tk_yyyy_xy, tk_yyyy_xyy, tk_yyyy_xyz, tk_yyyy_xz, tk_yyyy_xzz, tk_yyyy_yy, tk_yyyy_yyy, tk_yyyy_yyz, tk_yyyy_yz, tk_yyyy_yzz, tk_yyyy_zz, tk_yyyy_zzz, ts_xyyyy_xxx, ts_xyyyy_xxy, ts_xyyyy_xxz, ts_xyyyy_xyy, ts_xyyyy_xyz, ts_xyyyy_xzz, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yzz, ts_xyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_xxx[i] = 3.0 * tk_yyyy_xx[i] * fe_0 + tk_yyyy_xxx[i] * pa_x[i] + 2.0 * ts_xyyyy_xxx[i] * fz_0;

        tk_xyyyy_xxy[i] = 2.0 * tk_yyyy_xy[i] * fe_0 + tk_yyyy_xxy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxy[i] * fz_0;

        tk_xyyyy_xxz[i] = 2.0 * tk_yyyy_xz[i] * fe_0 + tk_yyyy_xxz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxz[i] * fz_0;

        tk_xyyyy_xyy[i] = tk_yyyy_yy[i] * fe_0 + tk_yyyy_xyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xyy[i] * fz_0;

        tk_xyyyy_xyz[i] = tk_yyyy_yz[i] * fe_0 + tk_yyyy_xyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyz[i] * fz_0;

        tk_xyyyy_xzz[i] = tk_yyyy_zz[i] * fe_0 + tk_yyyy_xzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xzz[i] * fz_0;

        tk_xyyyy_yyy[i] = tk_yyyy_yyy[i] * pa_x[i] + 2.0 * ts_xyyyy_yyy[i] * fz_0;

        tk_xyyyy_yyz[i] = tk_yyyy_yyz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyz[i] * fz_0;

        tk_xyyyy_yzz[i] = tk_yyyy_yzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yzz[i] * fz_0;

        tk_xyyyy_zzz[i] = tk_yyyy_zzz[i] * pa_x[i] + 2.0 * ts_xyyyy_zzz[i] * fz_0;
    }

    // Set up 110-120 components of targeted buffer : HF

    auto tk_xyyyz_xxx = pbuffer.data(idx_kin_hf + 110);

    auto tk_xyyyz_xxy = pbuffer.data(idx_kin_hf + 111);

    auto tk_xyyyz_xxz = pbuffer.data(idx_kin_hf + 112);

    auto tk_xyyyz_xyy = pbuffer.data(idx_kin_hf + 113);

    auto tk_xyyyz_xyz = pbuffer.data(idx_kin_hf + 114);

    auto tk_xyyyz_xzz = pbuffer.data(idx_kin_hf + 115);

    auto tk_xyyyz_yyy = pbuffer.data(idx_kin_hf + 116);

    auto tk_xyyyz_yyz = pbuffer.data(idx_kin_hf + 117);

    auto tk_xyyyz_yzz = pbuffer.data(idx_kin_hf + 118);

    auto tk_xyyyz_zzz = pbuffer.data(idx_kin_hf + 119);

    #pragma omp simd aligned(pa_x, pa_z, tk_xyyy_xxx, tk_xyyy_xxy, tk_xyyy_xyy, tk_xyyyz_xxx, tk_xyyyz_xxy, tk_xyyyz_xxz, tk_xyyyz_xyy, tk_xyyyz_xyz, tk_xyyyz_xzz, tk_xyyyz_yyy, tk_xyyyz_yyz, tk_xyyyz_yzz, tk_xyyyz_zzz, tk_yyyz_xxz, tk_yyyz_xyz, tk_yyyz_xz, tk_yyyz_xzz, tk_yyyz_yyy, tk_yyyz_yyz, tk_yyyz_yz, tk_yyyz_yzz, tk_yyyz_zz, tk_yyyz_zzz, ts_xyyyz_xxx, ts_xyyyz_xxy, ts_xyyyz_xxz, ts_xyyyz_xyy, ts_xyyyz_xyz, ts_xyyyz_xzz, ts_xyyyz_yyy, ts_xyyyz_yyz, ts_xyyyz_yzz, ts_xyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyz_xxx[i] = tk_xyyy_xxx[i] * pa_z[i] + 2.0 * ts_xyyyz_xxx[i] * fz_0;

        tk_xyyyz_xxy[i] = tk_xyyy_xxy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxy[i] * fz_0;

        tk_xyyyz_xxz[i] = 2.0 * tk_yyyz_xz[i] * fe_0 + tk_yyyz_xxz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxz[i] * fz_0;

        tk_xyyyz_xyy[i] = tk_xyyy_xyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xyy[i] * fz_0;

        tk_xyyyz_xyz[i] = tk_yyyz_yz[i] * fe_0 + tk_yyyz_xyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyz[i] * fz_0;

        tk_xyyyz_xzz[i] = tk_yyyz_zz[i] * fe_0 + tk_yyyz_xzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xzz[i] * fz_0;

        tk_xyyyz_yyy[i] = tk_yyyz_yyy[i] * pa_x[i] + 2.0 * ts_xyyyz_yyy[i] * fz_0;

        tk_xyyyz_yyz[i] = tk_yyyz_yyz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyz[i] * fz_0;

        tk_xyyyz_yzz[i] = tk_yyyz_yzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yzz[i] * fz_0;

        tk_xyyyz_zzz[i] = tk_yyyz_zzz[i] * pa_x[i] + 2.0 * ts_xyyyz_zzz[i] * fz_0;
    }

    // Set up 120-130 components of targeted buffer : HF

    auto tk_xyyzz_xxx = pbuffer.data(idx_kin_hf + 120);

    auto tk_xyyzz_xxy = pbuffer.data(idx_kin_hf + 121);

    auto tk_xyyzz_xxz = pbuffer.data(idx_kin_hf + 122);

    auto tk_xyyzz_xyy = pbuffer.data(idx_kin_hf + 123);

    auto tk_xyyzz_xyz = pbuffer.data(idx_kin_hf + 124);

    auto tk_xyyzz_xzz = pbuffer.data(idx_kin_hf + 125);

    auto tk_xyyzz_yyy = pbuffer.data(idx_kin_hf + 126);

    auto tk_xyyzz_yyz = pbuffer.data(idx_kin_hf + 127);

    auto tk_xyyzz_yzz = pbuffer.data(idx_kin_hf + 128);

    auto tk_xyyzz_zzz = pbuffer.data(idx_kin_hf + 129);

    #pragma omp simd aligned(pa_x, tk_xyyzz_xxx, tk_xyyzz_xxy, tk_xyyzz_xxz, tk_xyyzz_xyy, tk_xyyzz_xyz, tk_xyyzz_xzz, tk_xyyzz_yyy, tk_xyyzz_yyz, tk_xyyzz_yzz, tk_xyyzz_zzz, tk_yyzz_xx, tk_yyzz_xxx, tk_yyzz_xxy, tk_yyzz_xxz, tk_yyzz_xy, tk_yyzz_xyy, tk_yyzz_xyz, tk_yyzz_xz, tk_yyzz_xzz, tk_yyzz_yy, tk_yyzz_yyy, tk_yyzz_yyz, tk_yyzz_yz, tk_yyzz_yzz, tk_yyzz_zz, tk_yyzz_zzz, ts_xyyzz_xxx, ts_xyyzz_xxy, ts_xyyzz_xxz, ts_xyyzz_xyy, ts_xyyzz_xyz, ts_xyyzz_xzz, ts_xyyzz_yyy, ts_xyyzz_yyz, ts_xyyzz_yzz, ts_xyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_xxx[i] = 3.0 * tk_yyzz_xx[i] * fe_0 + tk_yyzz_xxx[i] * pa_x[i] + 2.0 * ts_xyyzz_xxx[i] * fz_0;

        tk_xyyzz_xxy[i] = 2.0 * tk_yyzz_xy[i] * fe_0 + tk_yyzz_xxy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxy[i] * fz_0;

        tk_xyyzz_xxz[i] = 2.0 * tk_yyzz_xz[i] * fe_0 + tk_yyzz_xxz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxz[i] * fz_0;

        tk_xyyzz_xyy[i] = tk_yyzz_yy[i] * fe_0 + tk_yyzz_xyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xyy[i] * fz_0;

        tk_xyyzz_xyz[i] = tk_yyzz_yz[i] * fe_0 + tk_yyzz_xyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyz[i] * fz_0;

        tk_xyyzz_xzz[i] = tk_yyzz_zz[i] * fe_0 + tk_yyzz_xzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xzz[i] * fz_0;

        tk_xyyzz_yyy[i] = tk_yyzz_yyy[i] * pa_x[i] + 2.0 * ts_xyyzz_yyy[i] * fz_0;

        tk_xyyzz_yyz[i] = tk_yyzz_yyz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyz[i] * fz_0;

        tk_xyyzz_yzz[i] = tk_yyzz_yzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yzz[i] * fz_0;

        tk_xyyzz_zzz[i] = tk_yyzz_zzz[i] * pa_x[i] + 2.0 * ts_xyyzz_zzz[i] * fz_0;
    }

    // Set up 130-140 components of targeted buffer : HF

    auto tk_xyzzz_xxx = pbuffer.data(idx_kin_hf + 130);

    auto tk_xyzzz_xxy = pbuffer.data(idx_kin_hf + 131);

    auto tk_xyzzz_xxz = pbuffer.data(idx_kin_hf + 132);

    auto tk_xyzzz_xyy = pbuffer.data(idx_kin_hf + 133);

    auto tk_xyzzz_xyz = pbuffer.data(idx_kin_hf + 134);

    auto tk_xyzzz_xzz = pbuffer.data(idx_kin_hf + 135);

    auto tk_xyzzz_yyy = pbuffer.data(idx_kin_hf + 136);

    auto tk_xyzzz_yyz = pbuffer.data(idx_kin_hf + 137);

    auto tk_xyzzz_yzz = pbuffer.data(idx_kin_hf + 138);

    auto tk_xyzzz_zzz = pbuffer.data(idx_kin_hf + 139);

    #pragma omp simd aligned(pa_x, pa_y, tk_xyzzz_xxx, tk_xyzzz_xxy, tk_xyzzz_xxz, tk_xyzzz_xyy, tk_xyzzz_xyz, tk_xyzzz_xzz, tk_xyzzz_yyy, tk_xyzzz_yyz, tk_xyzzz_yzz, tk_xyzzz_zzz, tk_xzzz_xxx, tk_xzzz_xxz, tk_xzzz_xzz, tk_yzzz_xxy, tk_yzzz_xy, tk_yzzz_xyy, tk_yzzz_xyz, tk_yzzz_yy, tk_yzzz_yyy, tk_yzzz_yyz, tk_yzzz_yz, tk_yzzz_yzz, tk_yzzz_zzz, ts_xyzzz_xxx, ts_xyzzz_xxy, ts_xyzzz_xxz, ts_xyzzz_xyy, ts_xyzzz_xyz, ts_xyzzz_xzz, ts_xyzzz_yyy, ts_xyzzz_yyz, ts_xyzzz_yzz, ts_xyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzz_xxx[i] = tk_xzzz_xxx[i] * pa_y[i] + 2.0 * ts_xyzzz_xxx[i] * fz_0;

        tk_xyzzz_xxy[i] = 2.0 * tk_yzzz_xy[i] * fe_0 + tk_yzzz_xxy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxy[i] * fz_0;

        tk_xyzzz_xxz[i] = tk_xzzz_xxz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxz[i] * fz_0;

        tk_xyzzz_xyy[i] = tk_yzzz_yy[i] * fe_0 + tk_yzzz_xyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xyy[i] * fz_0;

        tk_xyzzz_xyz[i] = tk_yzzz_yz[i] * fe_0 + tk_yzzz_xyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyz[i] * fz_0;

        tk_xyzzz_xzz[i] = tk_xzzz_xzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xzz[i] * fz_0;

        tk_xyzzz_yyy[i] = tk_yzzz_yyy[i] * pa_x[i] + 2.0 * ts_xyzzz_yyy[i] * fz_0;

        tk_xyzzz_yyz[i] = tk_yzzz_yyz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyz[i] * fz_0;

        tk_xyzzz_yzz[i] = tk_yzzz_yzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yzz[i] * fz_0;

        tk_xyzzz_zzz[i] = tk_yzzz_zzz[i] * pa_x[i] + 2.0 * ts_xyzzz_zzz[i] * fz_0;
    }

    // Set up 140-150 components of targeted buffer : HF

    auto tk_xzzzz_xxx = pbuffer.data(idx_kin_hf + 140);

    auto tk_xzzzz_xxy = pbuffer.data(idx_kin_hf + 141);

    auto tk_xzzzz_xxz = pbuffer.data(idx_kin_hf + 142);

    auto tk_xzzzz_xyy = pbuffer.data(idx_kin_hf + 143);

    auto tk_xzzzz_xyz = pbuffer.data(idx_kin_hf + 144);

    auto tk_xzzzz_xzz = pbuffer.data(idx_kin_hf + 145);

    auto tk_xzzzz_yyy = pbuffer.data(idx_kin_hf + 146);

    auto tk_xzzzz_yyz = pbuffer.data(idx_kin_hf + 147);

    auto tk_xzzzz_yzz = pbuffer.data(idx_kin_hf + 148);

    auto tk_xzzzz_zzz = pbuffer.data(idx_kin_hf + 149);

    #pragma omp simd aligned(pa_x, tk_xzzzz_xxx, tk_xzzzz_xxy, tk_xzzzz_xxz, tk_xzzzz_xyy, tk_xzzzz_xyz, tk_xzzzz_xzz, tk_xzzzz_yyy, tk_xzzzz_yyz, tk_xzzzz_yzz, tk_xzzzz_zzz, tk_zzzz_xx, tk_zzzz_xxx, tk_zzzz_xxy, tk_zzzz_xxz, tk_zzzz_xy, tk_zzzz_xyy, tk_zzzz_xyz, tk_zzzz_xz, tk_zzzz_xzz, tk_zzzz_yy, tk_zzzz_yyy, tk_zzzz_yyz, tk_zzzz_yz, tk_zzzz_yzz, tk_zzzz_zz, tk_zzzz_zzz, ts_xzzzz_xxx, ts_xzzzz_xxy, ts_xzzzz_xxz, ts_xzzzz_xyy, ts_xzzzz_xyz, ts_xzzzz_xzz, ts_xzzzz_yyy, ts_xzzzz_yyz, ts_xzzzz_yzz, ts_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_xxx[i] = 3.0 * tk_zzzz_xx[i] * fe_0 + tk_zzzz_xxx[i] * pa_x[i] + 2.0 * ts_xzzzz_xxx[i] * fz_0;

        tk_xzzzz_xxy[i] = 2.0 * tk_zzzz_xy[i] * fe_0 + tk_zzzz_xxy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxy[i] * fz_0;

        tk_xzzzz_xxz[i] = 2.0 * tk_zzzz_xz[i] * fe_0 + tk_zzzz_xxz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxz[i] * fz_0;

        tk_xzzzz_xyy[i] = tk_zzzz_yy[i] * fe_0 + tk_zzzz_xyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xyy[i] * fz_0;

        tk_xzzzz_xyz[i] = tk_zzzz_yz[i] * fe_0 + tk_zzzz_xyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyz[i] * fz_0;

        tk_xzzzz_xzz[i] = tk_zzzz_zz[i] * fe_0 + tk_zzzz_xzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xzz[i] * fz_0;

        tk_xzzzz_yyy[i] = tk_zzzz_yyy[i] * pa_x[i] + 2.0 * ts_xzzzz_yyy[i] * fz_0;

        tk_xzzzz_yyz[i] = tk_zzzz_yyz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyz[i] * fz_0;

        tk_xzzzz_yzz[i] = tk_zzzz_yzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yzz[i] * fz_0;

        tk_xzzzz_zzz[i] = tk_zzzz_zzz[i] * pa_x[i] + 2.0 * ts_xzzzz_zzz[i] * fz_0;
    }

    // Set up 150-160 components of targeted buffer : HF

    auto tk_yyyyy_xxx = pbuffer.data(idx_kin_hf + 150);

    auto tk_yyyyy_xxy = pbuffer.data(idx_kin_hf + 151);

    auto tk_yyyyy_xxz = pbuffer.data(idx_kin_hf + 152);

    auto tk_yyyyy_xyy = pbuffer.data(idx_kin_hf + 153);

    auto tk_yyyyy_xyz = pbuffer.data(idx_kin_hf + 154);

    auto tk_yyyyy_xzz = pbuffer.data(idx_kin_hf + 155);

    auto tk_yyyyy_yyy = pbuffer.data(idx_kin_hf + 156);

    auto tk_yyyyy_yyz = pbuffer.data(idx_kin_hf + 157);

    auto tk_yyyyy_yzz = pbuffer.data(idx_kin_hf + 158);

    auto tk_yyyyy_zzz = pbuffer.data(idx_kin_hf + 159);

    #pragma omp simd aligned(pa_y, tk_yyy_xxx, tk_yyy_xxy, tk_yyy_xxz, tk_yyy_xyy, tk_yyy_xyz, tk_yyy_xzz, tk_yyy_yyy, tk_yyy_yyz, tk_yyy_yzz, tk_yyy_zzz, tk_yyyy_xx, tk_yyyy_xxx, tk_yyyy_xxy, tk_yyyy_xxz, tk_yyyy_xy, tk_yyyy_xyy, tk_yyyy_xyz, tk_yyyy_xz, tk_yyyy_xzz, tk_yyyy_yy, tk_yyyy_yyy, tk_yyyy_yyz, tk_yyyy_yz, tk_yyyy_yzz, tk_yyyy_zz, tk_yyyy_zzz, tk_yyyyy_xxx, tk_yyyyy_xxy, tk_yyyyy_xxz, tk_yyyyy_xyy, tk_yyyyy_xyz, tk_yyyyy_xzz, tk_yyyyy_yyy, tk_yyyyy_yyz, tk_yyyyy_yzz, tk_yyyyy_zzz, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xzz, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yzz, ts_yyy_zzz, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xzz, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yzz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_xxx[i] = -8.0 * ts_yyy_xxx[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxx[i] * fe_0 + tk_yyyy_xxx[i] * pa_y[i] + 2.0 * ts_yyyyy_xxx[i] * fz_0;

        tk_yyyyy_xxy[i] = -8.0 * ts_yyy_xxy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxy[i] * fe_0 + tk_yyyy_xx[i] * fe_0 + tk_yyyy_xxy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxy[i] * fz_0;

        tk_yyyyy_xxz[i] = -8.0 * ts_yyy_xxz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxz[i] * fe_0 + tk_yyyy_xxz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxz[i] * fz_0;

        tk_yyyyy_xyy[i] = -8.0 * ts_yyy_xyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyy[i] * fe_0 + 2.0 * tk_yyyy_xy[i] * fe_0 + tk_yyyy_xyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xyy[i] * fz_0;

        tk_yyyyy_xyz[i] = -8.0 * ts_yyy_xyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyz[i] * fe_0 + tk_yyyy_xz[i] * fe_0 + tk_yyyy_xyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyz[i] * fz_0;

        tk_yyyyy_xzz[i] = -8.0 * ts_yyy_xzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xzz[i] * fe_0 + tk_yyyy_xzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xzz[i] * fz_0;

        tk_yyyyy_yyy[i] = -8.0 * ts_yyy_yyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyy[i] * fe_0 + 3.0 * tk_yyyy_yy[i] * fe_0 + tk_yyyy_yyy[i] * pa_y[i] + 2.0 * ts_yyyyy_yyy[i] * fz_0;

        tk_yyyyy_yyz[i] = -8.0 * ts_yyy_yyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyz[i] * fe_0 + 2.0 * tk_yyyy_yz[i] * fe_0 + tk_yyyy_yyz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyz[i] * fz_0;

        tk_yyyyy_yzz[i] = -8.0 * ts_yyy_yzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yzz[i] * fe_0 + tk_yyyy_zz[i] * fe_0 + tk_yyyy_yzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yzz[i] * fz_0;

        tk_yyyyy_zzz[i] = -8.0 * ts_yyy_zzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_zzz[i] * fe_0 + tk_yyyy_zzz[i] * pa_y[i] + 2.0 * ts_yyyyy_zzz[i] * fz_0;
    }

    // Set up 160-170 components of targeted buffer : HF

    auto tk_yyyyz_xxx = pbuffer.data(idx_kin_hf + 160);

    auto tk_yyyyz_xxy = pbuffer.data(idx_kin_hf + 161);

    auto tk_yyyyz_xxz = pbuffer.data(idx_kin_hf + 162);

    auto tk_yyyyz_xyy = pbuffer.data(idx_kin_hf + 163);

    auto tk_yyyyz_xyz = pbuffer.data(idx_kin_hf + 164);

    auto tk_yyyyz_xzz = pbuffer.data(idx_kin_hf + 165);

    auto tk_yyyyz_yyy = pbuffer.data(idx_kin_hf + 166);

    auto tk_yyyyz_yyz = pbuffer.data(idx_kin_hf + 167);

    auto tk_yyyyz_yzz = pbuffer.data(idx_kin_hf + 168);

    auto tk_yyyyz_zzz = pbuffer.data(idx_kin_hf + 169);

    #pragma omp simd aligned(pa_z, tk_yyyy_xx, tk_yyyy_xxx, tk_yyyy_xxy, tk_yyyy_xxz, tk_yyyy_xy, tk_yyyy_xyy, tk_yyyy_xyz, tk_yyyy_xz, tk_yyyy_xzz, tk_yyyy_yy, tk_yyyy_yyy, tk_yyyy_yyz, tk_yyyy_yz, tk_yyyy_yzz, tk_yyyy_zz, tk_yyyy_zzz, tk_yyyyz_xxx, tk_yyyyz_xxy, tk_yyyyz_xxz, tk_yyyyz_xyy, tk_yyyyz_xyz, tk_yyyyz_xzz, tk_yyyyz_yyy, tk_yyyyz_yyz, tk_yyyyz_yzz, tk_yyyyz_zzz, ts_yyyyz_xxx, ts_yyyyz_xxy, ts_yyyyz_xxz, ts_yyyyz_xyy, ts_yyyyz_xyz, ts_yyyyz_xzz, ts_yyyyz_yyy, ts_yyyyz_yyz, ts_yyyyz_yzz, ts_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_xxx[i] = tk_yyyy_xxx[i] * pa_z[i] + 2.0 * ts_yyyyz_xxx[i] * fz_0;

        tk_yyyyz_xxy[i] = tk_yyyy_xxy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxy[i] * fz_0;

        tk_yyyyz_xxz[i] = tk_yyyy_xx[i] * fe_0 + tk_yyyy_xxz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxz[i] * fz_0;

        tk_yyyyz_xyy[i] = tk_yyyy_xyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xyy[i] * fz_0;

        tk_yyyyz_xyz[i] = tk_yyyy_xy[i] * fe_0 + tk_yyyy_xyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyz[i] * fz_0;

        tk_yyyyz_xzz[i] = 2.0 * tk_yyyy_xz[i] * fe_0 + tk_yyyy_xzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xzz[i] * fz_0;

        tk_yyyyz_yyy[i] = tk_yyyy_yyy[i] * pa_z[i] + 2.0 * ts_yyyyz_yyy[i] * fz_0;

        tk_yyyyz_yyz[i] = tk_yyyy_yy[i] * fe_0 + tk_yyyy_yyz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyz[i] * fz_0;

        tk_yyyyz_yzz[i] = 2.0 * tk_yyyy_yz[i] * fe_0 + tk_yyyy_yzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yzz[i] * fz_0;

        tk_yyyyz_zzz[i] = 3.0 * tk_yyyy_zz[i] * fe_0 + tk_yyyy_zzz[i] * pa_z[i] + 2.0 * ts_yyyyz_zzz[i] * fz_0;
    }

    // Set up 170-180 components of targeted buffer : HF

    auto tk_yyyzz_xxx = pbuffer.data(idx_kin_hf + 170);

    auto tk_yyyzz_xxy = pbuffer.data(idx_kin_hf + 171);

    auto tk_yyyzz_xxz = pbuffer.data(idx_kin_hf + 172);

    auto tk_yyyzz_xyy = pbuffer.data(idx_kin_hf + 173);

    auto tk_yyyzz_xyz = pbuffer.data(idx_kin_hf + 174);

    auto tk_yyyzz_xzz = pbuffer.data(idx_kin_hf + 175);

    auto tk_yyyzz_yyy = pbuffer.data(idx_kin_hf + 176);

    auto tk_yyyzz_yyz = pbuffer.data(idx_kin_hf + 177);

    auto tk_yyyzz_yzz = pbuffer.data(idx_kin_hf + 178);

    auto tk_yyyzz_zzz = pbuffer.data(idx_kin_hf + 179);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyy_xxy, tk_yyy_xyy, tk_yyy_yyy, tk_yyyz_xxy, tk_yyyz_xyy, tk_yyyz_yyy, tk_yyyzz_xxx, tk_yyyzz_xxy, tk_yyyzz_xxz, tk_yyyzz_xyy, tk_yyyzz_xyz, tk_yyyzz_xzz, tk_yyyzz_yyy, tk_yyyzz_yyz, tk_yyyzz_yzz, tk_yyyzz_zzz, tk_yyzz_xxx, tk_yyzz_xxz, tk_yyzz_xyz, tk_yyzz_xz, tk_yyzz_xzz, tk_yyzz_yyz, tk_yyzz_yz, tk_yyzz_yzz, tk_yyzz_zz, tk_yyzz_zzz, tk_yzz_xxx, tk_yzz_xxz, tk_yzz_xyz, tk_yzz_xzz, tk_yzz_yyz, tk_yzz_yzz, tk_yzz_zzz, ts_yyy_xxy, ts_yyy_xyy, ts_yyy_yyy, ts_yyyzz_xxx, ts_yyyzz_xxy, ts_yyyzz_xxz, ts_yyyzz_xyy, ts_yyyzz_xyz, ts_yyyzz_xzz, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yzz, ts_yyyzz_zzz, ts_yzz_xxx, ts_yzz_xxz, ts_yzz_xyz, ts_yzz_xzz, ts_yzz_yyz, ts_yzz_yzz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_xxx[i] = -4.0 * ts_yzz_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxx[i] * fe_0 + tk_yyzz_xxx[i] * pa_y[i] + 2.0 * ts_yyyzz_xxx[i] * fz_0;

        tk_yyyzz_xxy[i] = -2.0 * ts_yyy_xxy[i] * fbe_0 * fz_0 + tk_yyy_xxy[i] * fe_0 + tk_yyyz_xxy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxy[i] * fz_0;

        tk_yyyzz_xxz[i] = -4.0 * ts_yzz_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxz[i] * fe_0 + tk_yyzz_xxz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxz[i] * fz_0;

        tk_yyyzz_xyy[i] = -2.0 * ts_yyy_xyy[i] * fbe_0 * fz_0 + tk_yyy_xyy[i] * fe_0 + tk_yyyz_xyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xyy[i] * fz_0;

        tk_yyyzz_xyz[i] = -4.0 * ts_yzz_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyz[i] * fe_0 + tk_yyzz_xz[i] * fe_0 + tk_yyzz_xyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyz[i] * fz_0;

        tk_yyyzz_xzz[i] = -4.0 * ts_yzz_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xzz[i] * fe_0 + tk_yyzz_xzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xzz[i] * fz_0;

        tk_yyyzz_yyy[i] = -2.0 * ts_yyy_yyy[i] * fbe_0 * fz_0 + tk_yyy_yyy[i] * fe_0 + tk_yyyz_yyy[i] * pa_z[i] + 2.0 * ts_yyyzz_yyy[i] * fz_0;

        tk_yyyzz_yyz[i] = -4.0 * ts_yzz_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyz[i] * fe_0 + 2.0 * tk_yyzz_yz[i] * fe_0 + tk_yyzz_yyz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyz[i] * fz_0;

        tk_yyyzz_yzz[i] = -4.0 * ts_yzz_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yzz[i] * fe_0 + tk_yyzz_zz[i] * fe_0 + tk_yyzz_yzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yzz[i] * fz_0;

        tk_yyyzz_zzz[i] = -4.0 * ts_yzz_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_zzz[i] * fe_0 + tk_yyzz_zzz[i] * pa_y[i] + 2.0 * ts_yyyzz_zzz[i] * fz_0;
    }

    // Set up 180-190 components of targeted buffer : HF

    auto tk_yyzzz_xxx = pbuffer.data(idx_kin_hf + 180);

    auto tk_yyzzz_xxy = pbuffer.data(idx_kin_hf + 181);

    auto tk_yyzzz_xxz = pbuffer.data(idx_kin_hf + 182);

    auto tk_yyzzz_xyy = pbuffer.data(idx_kin_hf + 183);

    auto tk_yyzzz_xyz = pbuffer.data(idx_kin_hf + 184);

    auto tk_yyzzz_xzz = pbuffer.data(idx_kin_hf + 185);

    auto tk_yyzzz_yyy = pbuffer.data(idx_kin_hf + 186);

    auto tk_yyzzz_yyz = pbuffer.data(idx_kin_hf + 187);

    auto tk_yyzzz_yzz = pbuffer.data(idx_kin_hf + 188);

    auto tk_yyzzz_zzz = pbuffer.data(idx_kin_hf + 189);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyz_xxy, tk_yyz_xyy, tk_yyz_yyy, tk_yyzz_xxy, tk_yyzz_xyy, tk_yyzz_yyy, tk_yyzzz_xxx, tk_yyzzz_xxy, tk_yyzzz_xxz, tk_yyzzz_xyy, tk_yyzzz_xyz, tk_yyzzz_xzz, tk_yyzzz_yyy, tk_yyzzz_yyz, tk_yyzzz_yzz, tk_yyzzz_zzz, tk_yzzz_xxx, tk_yzzz_xxz, tk_yzzz_xyz, tk_yzzz_xz, tk_yzzz_xzz, tk_yzzz_yyz, tk_yzzz_yz, tk_yzzz_yzz, tk_yzzz_zz, tk_yzzz_zzz, tk_zzz_xxx, tk_zzz_xxz, tk_zzz_xyz, tk_zzz_xzz, tk_zzz_yyz, tk_zzz_yzz, tk_zzz_zzz, ts_yyz_xxy, ts_yyz_xyy, ts_yyz_yyy, ts_yyzzz_xxx, ts_yyzzz_xxy, ts_yyzzz_xxz, ts_yyzzz_xyy, ts_yyzzz_xyz, ts_yyzzz_xzz, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yzz, ts_yyzzz_zzz, ts_zzz_xxx, ts_zzz_xxz, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_xxx[i] = -2.0 * ts_zzz_xxx[i] * fbe_0 * fz_0 + tk_zzz_xxx[i] * fe_0 + tk_yzzz_xxx[i] * pa_y[i] + 2.0 * ts_yyzzz_xxx[i] * fz_0;

        tk_yyzzz_xxy[i] = -4.0 * ts_yyz_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxy[i] * fe_0 + tk_yyzz_xxy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxy[i] * fz_0;

        tk_yyzzz_xxz[i] = -2.0 * ts_zzz_xxz[i] * fbe_0 * fz_0 + tk_zzz_xxz[i] * fe_0 + tk_yzzz_xxz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxz[i] * fz_0;

        tk_yyzzz_xyy[i] = -4.0 * ts_yyz_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xyy[i] * fe_0 + tk_yyzz_xyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xyy[i] * fz_0;

        tk_yyzzz_xyz[i] = -2.0 * ts_zzz_xyz[i] * fbe_0 * fz_0 + tk_zzz_xyz[i] * fe_0 + tk_yzzz_xz[i] * fe_0 + tk_yzzz_xyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyz[i] * fz_0;

        tk_yyzzz_xzz[i] = -2.0 * ts_zzz_xzz[i] * fbe_0 * fz_0 + tk_zzz_xzz[i] * fe_0 + tk_yzzz_xzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xzz[i] * fz_0;

        tk_yyzzz_yyy[i] = -4.0 * ts_yyz_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_yyy[i] * fe_0 + tk_yyzz_yyy[i] * pa_z[i] + 2.0 * ts_yyzzz_yyy[i] * fz_0;

        tk_yyzzz_yyz[i] = -2.0 * ts_zzz_yyz[i] * fbe_0 * fz_0 + tk_zzz_yyz[i] * fe_0 + 2.0 * tk_yzzz_yz[i] * fe_0 + tk_yzzz_yyz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyz[i] * fz_0;

        tk_yyzzz_yzz[i] = -2.0 * ts_zzz_yzz[i] * fbe_0 * fz_0 + tk_zzz_yzz[i] * fe_0 + tk_yzzz_zz[i] * fe_0 + tk_yzzz_yzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yzz[i] * fz_0;

        tk_yyzzz_zzz[i] = -2.0 * ts_zzz_zzz[i] * fbe_0 * fz_0 + tk_zzz_zzz[i] * fe_0 + tk_yzzz_zzz[i] * pa_y[i] + 2.0 * ts_yyzzz_zzz[i] * fz_0;
    }

    // Set up 190-200 components of targeted buffer : HF

    auto tk_yzzzz_xxx = pbuffer.data(idx_kin_hf + 190);

    auto tk_yzzzz_xxy = pbuffer.data(idx_kin_hf + 191);

    auto tk_yzzzz_xxz = pbuffer.data(idx_kin_hf + 192);

    auto tk_yzzzz_xyy = pbuffer.data(idx_kin_hf + 193);

    auto tk_yzzzz_xyz = pbuffer.data(idx_kin_hf + 194);

    auto tk_yzzzz_xzz = pbuffer.data(idx_kin_hf + 195);

    auto tk_yzzzz_yyy = pbuffer.data(idx_kin_hf + 196);

    auto tk_yzzzz_yyz = pbuffer.data(idx_kin_hf + 197);

    auto tk_yzzzz_yzz = pbuffer.data(idx_kin_hf + 198);

    auto tk_yzzzz_zzz = pbuffer.data(idx_kin_hf + 199);

    #pragma omp simd aligned(pa_y, tk_yzzzz_xxx, tk_yzzzz_xxy, tk_yzzzz_xxz, tk_yzzzz_xyy, tk_yzzzz_xyz, tk_yzzzz_xzz, tk_yzzzz_yyy, tk_yzzzz_yyz, tk_yzzzz_yzz, tk_yzzzz_zzz, tk_zzzz_xx, tk_zzzz_xxx, tk_zzzz_xxy, tk_zzzz_xxz, tk_zzzz_xy, tk_zzzz_xyy, tk_zzzz_xyz, tk_zzzz_xz, tk_zzzz_xzz, tk_zzzz_yy, tk_zzzz_yyy, tk_zzzz_yyz, tk_zzzz_yz, tk_zzzz_yzz, tk_zzzz_zz, tk_zzzz_zzz, ts_yzzzz_xxx, ts_yzzzz_xxy, ts_yzzzz_xxz, ts_yzzzz_xyy, ts_yzzzz_xyz, ts_yzzzz_xzz, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yzz, ts_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_xxx[i] = tk_zzzz_xxx[i] * pa_y[i] + 2.0 * ts_yzzzz_xxx[i] * fz_0;

        tk_yzzzz_xxy[i] = tk_zzzz_xx[i] * fe_0 + tk_zzzz_xxy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxy[i] * fz_0;

        tk_yzzzz_xxz[i] = tk_zzzz_xxz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxz[i] * fz_0;

        tk_yzzzz_xyy[i] = 2.0 * tk_zzzz_xy[i] * fe_0 + tk_zzzz_xyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xyy[i] * fz_0;

        tk_yzzzz_xyz[i] = tk_zzzz_xz[i] * fe_0 + tk_zzzz_xyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyz[i] * fz_0;

        tk_yzzzz_xzz[i] = tk_zzzz_xzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xzz[i] * fz_0;

        tk_yzzzz_yyy[i] = 3.0 * tk_zzzz_yy[i] * fe_0 + tk_zzzz_yyy[i] * pa_y[i] + 2.0 * ts_yzzzz_yyy[i] * fz_0;

        tk_yzzzz_yyz[i] = 2.0 * tk_zzzz_yz[i] * fe_0 + tk_zzzz_yyz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyz[i] * fz_0;

        tk_yzzzz_yzz[i] = tk_zzzz_zz[i] * fe_0 + tk_zzzz_yzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yzz[i] * fz_0;

        tk_yzzzz_zzz[i] = tk_zzzz_zzz[i] * pa_y[i] + 2.0 * ts_yzzzz_zzz[i] * fz_0;
    }

    // Set up 200-210 components of targeted buffer : HF

    auto tk_zzzzz_xxx = pbuffer.data(idx_kin_hf + 200);

    auto tk_zzzzz_xxy = pbuffer.data(idx_kin_hf + 201);

    auto tk_zzzzz_xxz = pbuffer.data(idx_kin_hf + 202);

    auto tk_zzzzz_xyy = pbuffer.data(idx_kin_hf + 203);

    auto tk_zzzzz_xyz = pbuffer.data(idx_kin_hf + 204);

    auto tk_zzzzz_xzz = pbuffer.data(idx_kin_hf + 205);

    auto tk_zzzzz_yyy = pbuffer.data(idx_kin_hf + 206);

    auto tk_zzzzz_yyz = pbuffer.data(idx_kin_hf + 207);

    auto tk_zzzzz_yzz = pbuffer.data(idx_kin_hf + 208);

    auto tk_zzzzz_zzz = pbuffer.data(idx_kin_hf + 209);

    #pragma omp simd aligned(pa_z, tk_zzz_xxx, tk_zzz_xxy, tk_zzz_xxz, tk_zzz_xyy, tk_zzz_xyz, tk_zzz_xzz, tk_zzz_yyy, tk_zzz_yyz, tk_zzz_yzz, tk_zzz_zzz, tk_zzzz_xx, tk_zzzz_xxx, tk_zzzz_xxy, tk_zzzz_xxz, tk_zzzz_xy, tk_zzzz_xyy, tk_zzzz_xyz, tk_zzzz_xz, tk_zzzz_xzz, tk_zzzz_yy, tk_zzzz_yyy, tk_zzzz_yyz, tk_zzzz_yz, tk_zzzz_yzz, tk_zzzz_zz, tk_zzzz_zzz, tk_zzzzz_xxx, tk_zzzzz_xxy, tk_zzzzz_xxz, tk_zzzzz_xyy, tk_zzzzz_xyz, tk_zzzzz_xzz, tk_zzzzz_yyy, tk_zzzzz_yyz, tk_zzzzz_yzz, tk_zzzzz_zzz, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xzz, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yzz, ts_zzz_zzz, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xzz, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yzz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_xxx[i] = -8.0 * ts_zzz_xxx[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxx[i] * fe_0 + tk_zzzz_xxx[i] * pa_z[i] + 2.0 * ts_zzzzz_xxx[i] * fz_0;

        tk_zzzzz_xxy[i] = -8.0 * ts_zzz_xxy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxy[i] * fe_0 + tk_zzzz_xxy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxy[i] * fz_0;

        tk_zzzzz_xxz[i] = -8.0 * ts_zzz_xxz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxz[i] * fe_0 + tk_zzzz_xx[i] * fe_0 + tk_zzzz_xxz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxz[i] * fz_0;

        tk_zzzzz_xyy[i] = -8.0 * ts_zzz_xyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyy[i] * fe_0 + tk_zzzz_xyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xyy[i] * fz_0;

        tk_zzzzz_xyz[i] = -8.0 * ts_zzz_xyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyz[i] * fe_0 + tk_zzzz_xy[i] * fe_0 + tk_zzzz_xyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyz[i] * fz_0;

        tk_zzzzz_xzz[i] = -8.0 * ts_zzz_xzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xzz[i] * fe_0 + 2.0 * tk_zzzz_xz[i] * fe_0 + tk_zzzz_xzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xzz[i] * fz_0;

        tk_zzzzz_yyy[i] = -8.0 * ts_zzz_yyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyy[i] * fe_0 + tk_zzzz_yyy[i] * pa_z[i] + 2.0 * ts_zzzzz_yyy[i] * fz_0;

        tk_zzzzz_yyz[i] = -8.0 * ts_zzz_yyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyz[i] * fe_0 + tk_zzzz_yy[i] * fe_0 + tk_zzzz_yyz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyz[i] * fz_0;

        tk_zzzzz_yzz[i] = -8.0 * ts_zzz_yzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yzz[i] * fe_0 + 2.0 * tk_zzzz_yz[i] * fe_0 + tk_zzzz_yzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yzz[i] * fz_0;

        tk_zzzzz_zzz[i] = -8.0 * ts_zzz_zzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_zzz[i] * fe_0 + 3.0 * tk_zzzz_zz[i] * fe_0 + tk_zzzz_zzz[i] * pa_z[i] + 2.0 * ts_zzzzz_zzz[i] * fz_0;
    }

}

} // kinrec namespace

