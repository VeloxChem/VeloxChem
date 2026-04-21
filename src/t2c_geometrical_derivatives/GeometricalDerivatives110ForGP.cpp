#include "GeometricalDerivatives110ForGP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_gp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_gp,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const int idx_op_hs,
                         const int idx_op_hd,
                         const int idx_op_ip,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tr_xxx_0 = pbuffer.data(idx_op_fs);

    auto tr_xxy_0 = pbuffer.data(idx_op_fs + 1);

    auto tr_xxz_0 = pbuffer.data(idx_op_fs + 2);

    auto tr_xyy_0 = pbuffer.data(idx_op_fs + 3);

    auto tr_xyz_0 = pbuffer.data(idx_op_fs + 4);

    auto tr_xzz_0 = pbuffer.data(idx_op_fs + 5);

    auto tr_yyy_0 = pbuffer.data(idx_op_fs + 6);

    auto tr_yyz_0 = pbuffer.data(idx_op_fs + 7);

    auto tr_yzz_0 = pbuffer.data(idx_op_fs + 8);

    auto tr_zzz_0 = pbuffer.data(idx_op_fs + 9);

    // Set up components of auxiliary buffer : FD

    auto tr_xxx_xx = pbuffer.data(idx_op_fd);

    auto tr_xxx_xy = pbuffer.data(idx_op_fd + 1);

    auto tr_xxx_xz = pbuffer.data(idx_op_fd + 2);

    auto tr_xxx_yy = pbuffer.data(idx_op_fd + 3);

    auto tr_xxx_yz = pbuffer.data(idx_op_fd + 4);

    auto tr_xxx_zz = pbuffer.data(idx_op_fd + 5);

    auto tr_xxy_xx = pbuffer.data(idx_op_fd + 6);

    auto tr_xxy_xy = pbuffer.data(idx_op_fd + 7);

    auto tr_xxy_xz = pbuffer.data(idx_op_fd + 8);

    auto tr_xxy_yy = pbuffer.data(idx_op_fd + 9);

    auto tr_xxy_yz = pbuffer.data(idx_op_fd + 10);

    auto tr_xxy_zz = pbuffer.data(idx_op_fd + 11);

    auto tr_xxz_xx = pbuffer.data(idx_op_fd + 12);

    auto tr_xxz_xy = pbuffer.data(idx_op_fd + 13);

    auto tr_xxz_xz = pbuffer.data(idx_op_fd + 14);

    auto tr_xxz_yy = pbuffer.data(idx_op_fd + 15);

    auto tr_xxz_yz = pbuffer.data(idx_op_fd + 16);

    auto tr_xxz_zz = pbuffer.data(idx_op_fd + 17);

    auto tr_xyy_xx = pbuffer.data(idx_op_fd + 18);

    auto tr_xyy_xy = pbuffer.data(idx_op_fd + 19);

    auto tr_xyy_xz = pbuffer.data(idx_op_fd + 20);

    auto tr_xyy_yy = pbuffer.data(idx_op_fd + 21);

    auto tr_xyy_yz = pbuffer.data(idx_op_fd + 22);

    auto tr_xyy_zz = pbuffer.data(idx_op_fd + 23);

    auto tr_xyz_xx = pbuffer.data(idx_op_fd + 24);

    auto tr_xyz_xy = pbuffer.data(idx_op_fd + 25);

    auto tr_xyz_xz = pbuffer.data(idx_op_fd + 26);

    auto tr_xyz_yy = pbuffer.data(idx_op_fd + 27);

    auto tr_xyz_yz = pbuffer.data(idx_op_fd + 28);

    auto tr_xyz_zz = pbuffer.data(idx_op_fd + 29);

    auto tr_xzz_xx = pbuffer.data(idx_op_fd + 30);

    auto tr_xzz_xy = pbuffer.data(idx_op_fd + 31);

    auto tr_xzz_xz = pbuffer.data(idx_op_fd + 32);

    auto tr_xzz_yy = pbuffer.data(idx_op_fd + 33);

    auto tr_xzz_yz = pbuffer.data(idx_op_fd + 34);

    auto tr_xzz_zz = pbuffer.data(idx_op_fd + 35);

    auto tr_yyy_xx = pbuffer.data(idx_op_fd + 36);

    auto tr_yyy_xy = pbuffer.data(idx_op_fd + 37);

    auto tr_yyy_xz = pbuffer.data(idx_op_fd + 38);

    auto tr_yyy_yy = pbuffer.data(idx_op_fd + 39);

    auto tr_yyy_yz = pbuffer.data(idx_op_fd + 40);

    auto tr_yyy_zz = pbuffer.data(idx_op_fd + 41);

    auto tr_yyz_xx = pbuffer.data(idx_op_fd + 42);

    auto tr_yyz_xy = pbuffer.data(idx_op_fd + 43);

    auto tr_yyz_xz = pbuffer.data(idx_op_fd + 44);

    auto tr_yyz_yy = pbuffer.data(idx_op_fd + 45);

    auto tr_yyz_yz = pbuffer.data(idx_op_fd + 46);

    auto tr_yyz_zz = pbuffer.data(idx_op_fd + 47);

    auto tr_yzz_xx = pbuffer.data(idx_op_fd + 48);

    auto tr_yzz_xy = pbuffer.data(idx_op_fd + 49);

    auto tr_yzz_xz = pbuffer.data(idx_op_fd + 50);

    auto tr_yzz_yy = pbuffer.data(idx_op_fd + 51);

    auto tr_yzz_yz = pbuffer.data(idx_op_fd + 52);

    auto tr_yzz_zz = pbuffer.data(idx_op_fd + 53);

    auto tr_zzz_xx = pbuffer.data(idx_op_fd + 54);

    auto tr_zzz_xy = pbuffer.data(idx_op_fd + 55);

    auto tr_zzz_xz = pbuffer.data(idx_op_fd + 56);

    auto tr_zzz_yy = pbuffer.data(idx_op_fd + 57);

    auto tr_zzz_yz = pbuffer.data(idx_op_fd + 58);

    auto tr_zzz_zz = pbuffer.data(idx_op_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto tr_xxxx_x = pbuffer.data(idx_op_gp);

    auto tr_xxxx_y = pbuffer.data(idx_op_gp + 1);

    auto tr_xxxx_z = pbuffer.data(idx_op_gp + 2);

    auto tr_xxxy_x = pbuffer.data(idx_op_gp + 3);

    auto tr_xxxy_y = pbuffer.data(idx_op_gp + 4);

    auto tr_xxxy_z = pbuffer.data(idx_op_gp + 5);

    auto tr_xxxz_x = pbuffer.data(idx_op_gp + 6);

    auto tr_xxxz_y = pbuffer.data(idx_op_gp + 7);

    auto tr_xxxz_z = pbuffer.data(idx_op_gp + 8);

    auto tr_xxyy_x = pbuffer.data(idx_op_gp + 9);

    auto tr_xxyy_y = pbuffer.data(idx_op_gp + 10);

    auto tr_xxyy_z = pbuffer.data(idx_op_gp + 11);

    auto tr_xxyz_x = pbuffer.data(idx_op_gp + 12);

    auto tr_xxyz_y = pbuffer.data(idx_op_gp + 13);

    auto tr_xxyz_z = pbuffer.data(idx_op_gp + 14);

    auto tr_xxzz_x = pbuffer.data(idx_op_gp + 15);

    auto tr_xxzz_y = pbuffer.data(idx_op_gp + 16);

    auto tr_xxzz_z = pbuffer.data(idx_op_gp + 17);

    auto tr_xyyy_x = pbuffer.data(idx_op_gp + 18);

    auto tr_xyyy_y = pbuffer.data(idx_op_gp + 19);

    auto tr_xyyy_z = pbuffer.data(idx_op_gp + 20);

    auto tr_xyyz_x = pbuffer.data(idx_op_gp + 21);

    auto tr_xyyz_y = pbuffer.data(idx_op_gp + 22);

    auto tr_xyyz_z = pbuffer.data(idx_op_gp + 23);

    auto tr_xyzz_x = pbuffer.data(idx_op_gp + 24);

    auto tr_xyzz_y = pbuffer.data(idx_op_gp + 25);

    auto tr_xyzz_z = pbuffer.data(idx_op_gp + 26);

    auto tr_xzzz_x = pbuffer.data(idx_op_gp + 27);

    auto tr_xzzz_y = pbuffer.data(idx_op_gp + 28);

    auto tr_xzzz_z = pbuffer.data(idx_op_gp + 29);

    auto tr_yyyy_x = pbuffer.data(idx_op_gp + 30);

    auto tr_yyyy_y = pbuffer.data(idx_op_gp + 31);

    auto tr_yyyy_z = pbuffer.data(idx_op_gp + 32);

    auto tr_yyyz_x = pbuffer.data(idx_op_gp + 33);

    auto tr_yyyz_y = pbuffer.data(idx_op_gp + 34);

    auto tr_yyyz_z = pbuffer.data(idx_op_gp + 35);

    auto tr_yyzz_x = pbuffer.data(idx_op_gp + 36);

    auto tr_yyzz_y = pbuffer.data(idx_op_gp + 37);

    auto tr_yyzz_z = pbuffer.data(idx_op_gp + 38);

    auto tr_yzzz_x = pbuffer.data(idx_op_gp + 39);

    auto tr_yzzz_y = pbuffer.data(idx_op_gp + 40);

    auto tr_yzzz_z = pbuffer.data(idx_op_gp + 41);

    auto tr_zzzz_x = pbuffer.data(idx_op_gp + 42);

    auto tr_zzzz_y = pbuffer.data(idx_op_gp + 43);

    auto tr_zzzz_z = pbuffer.data(idx_op_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto tr_xxxxx_0 = pbuffer.data(idx_op_hs);

    auto tr_xxxxy_0 = pbuffer.data(idx_op_hs + 1);

    auto tr_xxxxz_0 = pbuffer.data(idx_op_hs + 2);

    auto tr_xxxyy_0 = pbuffer.data(idx_op_hs + 3);

    auto tr_xxxyz_0 = pbuffer.data(idx_op_hs + 4);

    auto tr_xxxzz_0 = pbuffer.data(idx_op_hs + 5);

    auto tr_xxyyy_0 = pbuffer.data(idx_op_hs + 6);

    auto tr_xxyyz_0 = pbuffer.data(idx_op_hs + 7);

    auto tr_xxyzz_0 = pbuffer.data(idx_op_hs + 8);

    auto tr_xxzzz_0 = pbuffer.data(idx_op_hs + 9);

    auto tr_xyyyy_0 = pbuffer.data(idx_op_hs + 10);

    auto tr_xyyyz_0 = pbuffer.data(idx_op_hs + 11);

    auto tr_xyyzz_0 = pbuffer.data(idx_op_hs + 12);

    auto tr_xyzzz_0 = pbuffer.data(idx_op_hs + 13);

    auto tr_xzzzz_0 = pbuffer.data(idx_op_hs + 14);

    auto tr_yyyyy_0 = pbuffer.data(idx_op_hs + 15);

    auto tr_yyyyz_0 = pbuffer.data(idx_op_hs + 16);

    auto tr_yyyzz_0 = pbuffer.data(idx_op_hs + 17);

    auto tr_yyzzz_0 = pbuffer.data(idx_op_hs + 18);

    auto tr_yzzzz_0 = pbuffer.data(idx_op_hs + 19);

    auto tr_zzzzz_0 = pbuffer.data(idx_op_hs + 20);

    // Set up components of auxiliary buffer : HD

    auto tr_xxxxx_xx = pbuffer.data(idx_op_hd);

    auto tr_xxxxx_xy = pbuffer.data(idx_op_hd + 1);

    auto tr_xxxxx_xz = pbuffer.data(idx_op_hd + 2);

    auto tr_xxxxx_yy = pbuffer.data(idx_op_hd + 3);

    auto tr_xxxxx_yz = pbuffer.data(idx_op_hd + 4);

    auto tr_xxxxx_zz = pbuffer.data(idx_op_hd + 5);

    auto tr_xxxxy_xx = pbuffer.data(idx_op_hd + 6);

    auto tr_xxxxy_xy = pbuffer.data(idx_op_hd + 7);

    auto tr_xxxxy_xz = pbuffer.data(idx_op_hd + 8);

    auto tr_xxxxy_yy = pbuffer.data(idx_op_hd + 9);

    auto tr_xxxxy_yz = pbuffer.data(idx_op_hd + 10);

    auto tr_xxxxy_zz = pbuffer.data(idx_op_hd + 11);

    auto tr_xxxxz_xx = pbuffer.data(idx_op_hd + 12);

    auto tr_xxxxz_xy = pbuffer.data(idx_op_hd + 13);

    auto tr_xxxxz_xz = pbuffer.data(idx_op_hd + 14);

    auto tr_xxxxz_yy = pbuffer.data(idx_op_hd + 15);

    auto tr_xxxxz_yz = pbuffer.data(idx_op_hd + 16);

    auto tr_xxxxz_zz = pbuffer.data(idx_op_hd + 17);

    auto tr_xxxyy_xx = pbuffer.data(idx_op_hd + 18);

    auto tr_xxxyy_xy = pbuffer.data(idx_op_hd + 19);

    auto tr_xxxyy_xz = pbuffer.data(idx_op_hd + 20);

    auto tr_xxxyy_yy = pbuffer.data(idx_op_hd + 21);

    auto tr_xxxyy_yz = pbuffer.data(idx_op_hd + 22);

    auto tr_xxxyy_zz = pbuffer.data(idx_op_hd + 23);

    auto tr_xxxyz_xx = pbuffer.data(idx_op_hd + 24);

    auto tr_xxxyz_xy = pbuffer.data(idx_op_hd + 25);

    auto tr_xxxyz_xz = pbuffer.data(idx_op_hd + 26);

    auto tr_xxxyz_yy = pbuffer.data(idx_op_hd + 27);

    auto tr_xxxyz_yz = pbuffer.data(idx_op_hd + 28);

    auto tr_xxxyz_zz = pbuffer.data(idx_op_hd + 29);

    auto tr_xxxzz_xx = pbuffer.data(idx_op_hd + 30);

    auto tr_xxxzz_xy = pbuffer.data(idx_op_hd + 31);

    auto tr_xxxzz_xz = pbuffer.data(idx_op_hd + 32);

    auto tr_xxxzz_yy = pbuffer.data(idx_op_hd + 33);

    auto tr_xxxzz_yz = pbuffer.data(idx_op_hd + 34);

    auto tr_xxxzz_zz = pbuffer.data(idx_op_hd + 35);

    auto tr_xxyyy_xx = pbuffer.data(idx_op_hd + 36);

    auto tr_xxyyy_xy = pbuffer.data(idx_op_hd + 37);

    auto tr_xxyyy_xz = pbuffer.data(idx_op_hd + 38);

    auto tr_xxyyy_yy = pbuffer.data(idx_op_hd + 39);

    auto tr_xxyyy_yz = pbuffer.data(idx_op_hd + 40);

    auto tr_xxyyy_zz = pbuffer.data(idx_op_hd + 41);

    auto tr_xxyyz_xx = pbuffer.data(idx_op_hd + 42);

    auto tr_xxyyz_xy = pbuffer.data(idx_op_hd + 43);

    auto tr_xxyyz_xz = pbuffer.data(idx_op_hd + 44);

    auto tr_xxyyz_yy = pbuffer.data(idx_op_hd + 45);

    auto tr_xxyyz_yz = pbuffer.data(idx_op_hd + 46);

    auto tr_xxyyz_zz = pbuffer.data(idx_op_hd + 47);

    auto tr_xxyzz_xx = pbuffer.data(idx_op_hd + 48);

    auto tr_xxyzz_xy = pbuffer.data(idx_op_hd + 49);

    auto tr_xxyzz_xz = pbuffer.data(idx_op_hd + 50);

    auto tr_xxyzz_yy = pbuffer.data(idx_op_hd + 51);

    auto tr_xxyzz_yz = pbuffer.data(idx_op_hd + 52);

    auto tr_xxyzz_zz = pbuffer.data(idx_op_hd + 53);

    auto tr_xxzzz_xx = pbuffer.data(idx_op_hd + 54);

    auto tr_xxzzz_xy = pbuffer.data(idx_op_hd + 55);

    auto tr_xxzzz_xz = pbuffer.data(idx_op_hd + 56);

    auto tr_xxzzz_yy = pbuffer.data(idx_op_hd + 57);

    auto tr_xxzzz_yz = pbuffer.data(idx_op_hd + 58);

    auto tr_xxzzz_zz = pbuffer.data(idx_op_hd + 59);

    auto tr_xyyyy_xx = pbuffer.data(idx_op_hd + 60);

    auto tr_xyyyy_xy = pbuffer.data(idx_op_hd + 61);

    auto tr_xyyyy_xz = pbuffer.data(idx_op_hd + 62);

    auto tr_xyyyy_yy = pbuffer.data(idx_op_hd + 63);

    auto tr_xyyyy_yz = pbuffer.data(idx_op_hd + 64);

    auto tr_xyyyy_zz = pbuffer.data(idx_op_hd + 65);

    auto tr_xyyyz_xx = pbuffer.data(idx_op_hd + 66);

    auto tr_xyyyz_xy = pbuffer.data(idx_op_hd + 67);

    auto tr_xyyyz_xz = pbuffer.data(idx_op_hd + 68);

    auto tr_xyyyz_yy = pbuffer.data(idx_op_hd + 69);

    auto tr_xyyyz_yz = pbuffer.data(idx_op_hd + 70);

    auto tr_xyyyz_zz = pbuffer.data(idx_op_hd + 71);

    auto tr_xyyzz_xx = pbuffer.data(idx_op_hd + 72);

    auto tr_xyyzz_xy = pbuffer.data(idx_op_hd + 73);

    auto tr_xyyzz_xz = pbuffer.data(idx_op_hd + 74);

    auto tr_xyyzz_yy = pbuffer.data(idx_op_hd + 75);

    auto tr_xyyzz_yz = pbuffer.data(idx_op_hd + 76);

    auto tr_xyyzz_zz = pbuffer.data(idx_op_hd + 77);

    auto tr_xyzzz_xx = pbuffer.data(idx_op_hd + 78);

    auto tr_xyzzz_xy = pbuffer.data(idx_op_hd + 79);

    auto tr_xyzzz_xz = pbuffer.data(idx_op_hd + 80);

    auto tr_xyzzz_yy = pbuffer.data(idx_op_hd + 81);

    auto tr_xyzzz_yz = pbuffer.data(idx_op_hd + 82);

    auto tr_xyzzz_zz = pbuffer.data(idx_op_hd + 83);

    auto tr_xzzzz_xx = pbuffer.data(idx_op_hd + 84);

    auto tr_xzzzz_xy = pbuffer.data(idx_op_hd + 85);

    auto tr_xzzzz_xz = pbuffer.data(idx_op_hd + 86);

    auto tr_xzzzz_yy = pbuffer.data(idx_op_hd + 87);

    auto tr_xzzzz_yz = pbuffer.data(idx_op_hd + 88);

    auto tr_xzzzz_zz = pbuffer.data(idx_op_hd + 89);

    auto tr_yyyyy_xx = pbuffer.data(idx_op_hd + 90);

    auto tr_yyyyy_xy = pbuffer.data(idx_op_hd + 91);

    auto tr_yyyyy_xz = pbuffer.data(idx_op_hd + 92);

    auto tr_yyyyy_yy = pbuffer.data(idx_op_hd + 93);

    auto tr_yyyyy_yz = pbuffer.data(idx_op_hd + 94);

    auto tr_yyyyy_zz = pbuffer.data(idx_op_hd + 95);

    auto tr_yyyyz_xx = pbuffer.data(idx_op_hd + 96);

    auto tr_yyyyz_xy = pbuffer.data(idx_op_hd + 97);

    auto tr_yyyyz_xz = pbuffer.data(idx_op_hd + 98);

    auto tr_yyyyz_yy = pbuffer.data(idx_op_hd + 99);

    auto tr_yyyyz_yz = pbuffer.data(idx_op_hd + 100);

    auto tr_yyyyz_zz = pbuffer.data(idx_op_hd + 101);

    auto tr_yyyzz_xx = pbuffer.data(idx_op_hd + 102);

    auto tr_yyyzz_xy = pbuffer.data(idx_op_hd + 103);

    auto tr_yyyzz_xz = pbuffer.data(idx_op_hd + 104);

    auto tr_yyyzz_yy = pbuffer.data(idx_op_hd + 105);

    auto tr_yyyzz_yz = pbuffer.data(idx_op_hd + 106);

    auto tr_yyyzz_zz = pbuffer.data(idx_op_hd + 107);

    auto tr_yyzzz_xx = pbuffer.data(idx_op_hd + 108);

    auto tr_yyzzz_xy = pbuffer.data(idx_op_hd + 109);

    auto tr_yyzzz_xz = pbuffer.data(idx_op_hd + 110);

    auto tr_yyzzz_yy = pbuffer.data(idx_op_hd + 111);

    auto tr_yyzzz_yz = pbuffer.data(idx_op_hd + 112);

    auto tr_yyzzz_zz = pbuffer.data(idx_op_hd + 113);

    auto tr_yzzzz_xx = pbuffer.data(idx_op_hd + 114);

    auto tr_yzzzz_xy = pbuffer.data(idx_op_hd + 115);

    auto tr_yzzzz_xz = pbuffer.data(idx_op_hd + 116);

    auto tr_yzzzz_yy = pbuffer.data(idx_op_hd + 117);

    auto tr_yzzzz_yz = pbuffer.data(idx_op_hd + 118);

    auto tr_yzzzz_zz = pbuffer.data(idx_op_hd + 119);

    auto tr_zzzzz_xx = pbuffer.data(idx_op_hd + 120);

    auto tr_zzzzz_xy = pbuffer.data(idx_op_hd + 121);

    auto tr_zzzzz_xz = pbuffer.data(idx_op_hd + 122);

    auto tr_zzzzz_yy = pbuffer.data(idx_op_hd + 123);

    auto tr_zzzzz_yz = pbuffer.data(idx_op_hd + 124);

    auto tr_zzzzz_zz = pbuffer.data(idx_op_hd + 125);

    // Set up components of auxiliary buffer : IP

    auto tr_xxxxxx_x = pbuffer.data(idx_op_ip);

    auto tr_xxxxxx_y = pbuffer.data(idx_op_ip + 1);

    auto tr_xxxxxx_z = pbuffer.data(idx_op_ip + 2);

    auto tr_xxxxxy_x = pbuffer.data(idx_op_ip + 3);

    auto tr_xxxxxy_y = pbuffer.data(idx_op_ip + 4);

    auto tr_xxxxxy_z = pbuffer.data(idx_op_ip + 5);

    auto tr_xxxxxz_x = pbuffer.data(idx_op_ip + 6);

    auto tr_xxxxxz_y = pbuffer.data(idx_op_ip + 7);

    auto tr_xxxxxz_z = pbuffer.data(idx_op_ip + 8);

    auto tr_xxxxyy_x = pbuffer.data(idx_op_ip + 9);

    auto tr_xxxxyy_y = pbuffer.data(idx_op_ip + 10);

    auto tr_xxxxyy_z = pbuffer.data(idx_op_ip + 11);

    auto tr_xxxxyz_x = pbuffer.data(idx_op_ip + 12);

    auto tr_xxxxyz_y = pbuffer.data(idx_op_ip + 13);

    auto tr_xxxxyz_z = pbuffer.data(idx_op_ip + 14);

    auto tr_xxxxzz_x = pbuffer.data(idx_op_ip + 15);

    auto tr_xxxxzz_y = pbuffer.data(idx_op_ip + 16);

    auto tr_xxxxzz_z = pbuffer.data(idx_op_ip + 17);

    auto tr_xxxyyy_x = pbuffer.data(idx_op_ip + 18);

    auto tr_xxxyyy_y = pbuffer.data(idx_op_ip + 19);

    auto tr_xxxyyy_z = pbuffer.data(idx_op_ip + 20);

    auto tr_xxxyyz_x = pbuffer.data(idx_op_ip + 21);

    auto tr_xxxyyz_y = pbuffer.data(idx_op_ip + 22);

    auto tr_xxxyyz_z = pbuffer.data(idx_op_ip + 23);

    auto tr_xxxyzz_x = pbuffer.data(idx_op_ip + 24);

    auto tr_xxxyzz_y = pbuffer.data(idx_op_ip + 25);

    auto tr_xxxyzz_z = pbuffer.data(idx_op_ip + 26);

    auto tr_xxxzzz_x = pbuffer.data(idx_op_ip + 27);

    auto tr_xxxzzz_y = pbuffer.data(idx_op_ip + 28);

    auto tr_xxxzzz_z = pbuffer.data(idx_op_ip + 29);

    auto tr_xxyyyy_x = pbuffer.data(idx_op_ip + 30);

    auto tr_xxyyyy_y = pbuffer.data(idx_op_ip + 31);

    auto tr_xxyyyy_z = pbuffer.data(idx_op_ip + 32);

    auto tr_xxyyyz_x = pbuffer.data(idx_op_ip + 33);

    auto tr_xxyyyz_y = pbuffer.data(idx_op_ip + 34);

    auto tr_xxyyyz_z = pbuffer.data(idx_op_ip + 35);

    auto tr_xxyyzz_x = pbuffer.data(idx_op_ip + 36);

    auto tr_xxyyzz_y = pbuffer.data(idx_op_ip + 37);

    auto tr_xxyyzz_z = pbuffer.data(idx_op_ip + 38);

    auto tr_xxyzzz_x = pbuffer.data(idx_op_ip + 39);

    auto tr_xxyzzz_y = pbuffer.data(idx_op_ip + 40);

    auto tr_xxyzzz_z = pbuffer.data(idx_op_ip + 41);

    auto tr_xxzzzz_x = pbuffer.data(idx_op_ip + 42);

    auto tr_xxzzzz_y = pbuffer.data(idx_op_ip + 43);

    auto tr_xxzzzz_z = pbuffer.data(idx_op_ip + 44);

    auto tr_xyyyyy_x = pbuffer.data(idx_op_ip + 45);

    auto tr_xyyyyy_y = pbuffer.data(idx_op_ip + 46);

    auto tr_xyyyyy_z = pbuffer.data(idx_op_ip + 47);

    auto tr_xyyyyz_x = pbuffer.data(idx_op_ip + 48);

    auto tr_xyyyyz_y = pbuffer.data(idx_op_ip + 49);

    auto tr_xyyyyz_z = pbuffer.data(idx_op_ip + 50);

    auto tr_xyyyzz_x = pbuffer.data(idx_op_ip + 51);

    auto tr_xyyyzz_y = pbuffer.data(idx_op_ip + 52);

    auto tr_xyyyzz_z = pbuffer.data(idx_op_ip + 53);

    auto tr_xyyzzz_x = pbuffer.data(idx_op_ip + 54);

    auto tr_xyyzzz_y = pbuffer.data(idx_op_ip + 55);

    auto tr_xyyzzz_z = pbuffer.data(idx_op_ip + 56);

    auto tr_xyzzzz_x = pbuffer.data(idx_op_ip + 57);

    auto tr_xyzzzz_y = pbuffer.data(idx_op_ip + 58);

    auto tr_xyzzzz_z = pbuffer.data(idx_op_ip + 59);

    auto tr_xzzzzz_x = pbuffer.data(idx_op_ip + 60);

    auto tr_xzzzzz_y = pbuffer.data(idx_op_ip + 61);

    auto tr_xzzzzz_z = pbuffer.data(idx_op_ip + 62);

    auto tr_yyyyyy_x = pbuffer.data(idx_op_ip + 63);

    auto tr_yyyyyy_y = pbuffer.data(idx_op_ip + 64);

    auto tr_yyyyyy_z = pbuffer.data(idx_op_ip + 65);

    auto tr_yyyyyz_x = pbuffer.data(idx_op_ip + 66);

    auto tr_yyyyyz_y = pbuffer.data(idx_op_ip + 67);

    auto tr_yyyyyz_z = pbuffer.data(idx_op_ip + 68);

    auto tr_yyyyzz_x = pbuffer.data(idx_op_ip + 69);

    auto tr_yyyyzz_y = pbuffer.data(idx_op_ip + 70);

    auto tr_yyyyzz_z = pbuffer.data(idx_op_ip + 71);

    auto tr_yyyzzz_x = pbuffer.data(idx_op_ip + 72);

    auto tr_yyyzzz_y = pbuffer.data(idx_op_ip + 73);

    auto tr_yyyzzz_z = pbuffer.data(idx_op_ip + 74);

    auto tr_yyzzzz_x = pbuffer.data(idx_op_ip + 75);

    auto tr_yyzzzz_y = pbuffer.data(idx_op_ip + 76);

    auto tr_yyzzzz_z = pbuffer.data(idx_op_ip + 77);

    auto tr_yzzzzz_x = pbuffer.data(idx_op_ip + 78);

    auto tr_yzzzzz_y = pbuffer.data(idx_op_ip + 79);

    auto tr_yzzzzz_z = pbuffer.data(idx_op_ip + 80);

    auto tr_zzzzzz_x = pbuffer.data(idx_op_ip + 81);

    auto tr_zzzzzz_y = pbuffer.data(idx_op_ip + 82);

    auto tr_zzzzzz_z = pbuffer.data(idx_op_ip + 83);

    // Set up 0-3 components of targeted buffer : GP

    auto tr_x_0_x_xxxx_x = pbuffer.data(idx_op_geom_110_gp);

    auto tr_x_0_x_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 1);

    auto tr_x_0_x_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 2);

    #pragma omp simd aligned(tr_x_0_x_xxxx_x, tr_x_0_x_xxxx_y, tr_x_0_x_xxxx_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxx_0, tr_xxxxx_xx, tr_xxxxx_xy, tr_xxxxx_xz, tr_xxxxxx_x, tr_xxxxxx_y, tr_xxxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxx_x[i] = 12.0 * tr_xx_x[i] + 4.0 * tr_xxx_0[i] - 8.0 * tr_xxx_xx[i] * tke_0 - 18.0 * tr_xxxx_x[i] * tbe_0 - 2.0 * tr_xxxxx_0[i] * tbe_0 + 4.0 * tr_xxxxx_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_y[i] = 12.0 * tr_xx_y[i] - 8.0 * tr_xxx_xy[i] * tke_0 - 18.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_z[i] = 12.0 * tr_xx_z[i] - 8.0 * tr_xxx_xz[i] * tke_0 - 18.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_z[i] * tbe_0 * tbe_0;
    }

    // Set up 3-6 components of targeted buffer : GP

    auto tr_x_0_x_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 3);

    auto tr_x_0_x_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 4);

    auto tr_x_0_x_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 5);

    #pragma omp simd aligned(tr_x_0_x_xxxy_x, tr_x_0_x_xxxy_y, tr_x_0_x_xxxy_z, tr_xxxxxy_x, tr_xxxxxy_y, tr_xxxxxy_z, tr_xxxxy_0, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xy_x, tr_xy_y, tr_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxy_x[i] = 6.0 * tr_xy_x[i] + 3.0 * tr_xxy_0[i] - 6.0 * tr_xxy_xx[i] * tke_0 - 14.0 * tr_xxxy_x[i] * tbe_0 - 2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_y[i] = 6.0 * tr_xy_y[i] - 6.0 * tr_xxy_xy[i] * tke_0 - 14.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_z[i] = 6.0 * tr_xy_z[i] - 6.0 * tr_xxy_xz[i] * tke_0 - 14.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 6-9 components of targeted buffer : GP

    auto tr_x_0_x_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 6);

    auto tr_x_0_x_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 7);

    auto tr_x_0_x_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 8);

    #pragma omp simd aligned(tr_x_0_x_xxxz_x, tr_x_0_x_xxxz_y, tr_x_0_x_xxxz_z, tr_xxxxxz_x, tr_xxxxxz_y, tr_xxxxxz_z, tr_xxxxz_0, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxz_x[i] = 6.0 * tr_xz_x[i] + 3.0 * tr_xxz_0[i] - 6.0 * tr_xxz_xx[i] * tke_0 - 14.0 * tr_xxxz_x[i] * tbe_0 - 2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_y[i] = 6.0 * tr_xz_y[i] - 6.0 * tr_xxz_xy[i] * tke_0 - 14.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_z[i] = 6.0 * tr_xz_z[i] - 6.0 * tr_xxz_xz[i] * tke_0 - 14.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 9-12 components of targeted buffer : GP

    auto tr_x_0_x_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 9);

    auto tr_x_0_x_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 10);

    auto tr_x_0_x_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 11);

    #pragma omp simd aligned(tr_x_0_x_xxyy_x, tr_x_0_x_xxyy_y, tr_x_0_x_xxyy_z, tr_xxxxyy_x, tr_xxxxyy_y, tr_xxxxyy_z, tr_xxxyy_0, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxyy_x[i] = 2.0 * tr_yy_x[i] + 2.0 * tr_xyy_0[i] - 4.0 * tr_xyy_xx[i] * tke_0 - 10.0 * tr_xxyy_x[i] * tbe_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_y[i] = 2.0 * tr_yy_y[i] - 4.0 * tr_xyy_xy[i] * tke_0 - 10.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_z[i] = 2.0 * tr_yy_z[i] - 4.0 * tr_xyy_xz[i] * tke_0 - 10.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 12-15 components of targeted buffer : GP

    auto tr_x_0_x_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 12);

    auto tr_x_0_x_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 13);

    auto tr_x_0_x_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 14);

    #pragma omp simd aligned(tr_x_0_x_xxyz_x, tr_x_0_x_xxyz_y, tr_x_0_x_xxyz_z, tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxxyz_0, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxyz_x[i] = 2.0 * tr_yz_x[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_xx[i] * tke_0 - 10.0 * tr_xxyz_x[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_y[i] = 2.0 * tr_yz_y[i] - 4.0 * tr_xyz_xy[i] * tke_0 - 10.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_z[i] = 2.0 * tr_yz_z[i] - 4.0 * tr_xyz_xz[i] * tke_0 - 10.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 15-18 components of targeted buffer : GP

    auto tr_x_0_x_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 15);

    auto tr_x_0_x_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 16);

    auto tr_x_0_x_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 17);

    #pragma omp simd aligned(tr_x_0_x_xxzz_x, tr_x_0_x_xxzz_y, tr_x_0_x_xxzz_z, tr_xxxxzz_x, tr_xxxxzz_y, tr_xxxxzz_z, tr_xxxzz_0, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxzz_x[i] = 2.0 * tr_zz_x[i] + 2.0 * tr_xzz_0[i] - 4.0 * tr_xzz_xx[i] * tke_0 - 10.0 * tr_xxzz_x[i] * tbe_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_y[i] = 2.0 * tr_zz_y[i] - 4.0 * tr_xzz_xy[i] * tke_0 - 10.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_z[i] = 2.0 * tr_zz_z[i] - 4.0 * tr_xzz_xz[i] * tke_0 - 10.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 18-21 components of targeted buffer : GP

    auto tr_x_0_x_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 18);

    auto tr_x_0_x_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 19);

    auto tr_x_0_x_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 20);

    #pragma omp simd aligned(tr_x_0_x_xyyy_x, tr_x_0_x_xyyy_y, tr_x_0_x_xyyy_z, tr_xxxyyy_x, tr_xxxyyy_y, tr_xxxyyy_z, tr_xxyyy_0, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyyy_x[i] = tr_yyy_0[i] - 2.0 * tr_yyy_xx[i] * tke_0 - 6.0 * tr_xyyy_x[i] * tbe_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_y[i] = -2.0 * tr_yyy_xy[i] * tke_0 - 6.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_z[i] = -2.0 * tr_yyy_xz[i] * tke_0 - 6.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 21-24 components of targeted buffer : GP

    auto tr_x_0_x_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 21);

    auto tr_x_0_x_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 22);

    auto tr_x_0_x_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 23);

    #pragma omp simd aligned(tr_x_0_x_xyyz_x, tr_x_0_x_xyyz_y, tr_x_0_x_xyyz_z, tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxyyz_0, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyyz_x[i] = tr_yyz_0[i] - 2.0 * tr_yyz_xx[i] * tke_0 - 6.0 * tr_xyyz_x[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_y[i] = -2.0 * tr_yyz_xy[i] * tke_0 - 6.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_z[i] = -2.0 * tr_yyz_xz[i] * tke_0 - 6.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 24-27 components of targeted buffer : GP

    auto tr_x_0_x_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 24);

    auto tr_x_0_x_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 25);

    auto tr_x_0_x_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 26);

    #pragma omp simd aligned(tr_x_0_x_xyzz_x, tr_x_0_x_xyzz_y, tr_x_0_x_xyzz_z, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xxyzz_0, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyzz_x[i] = tr_yzz_0[i] - 2.0 * tr_yzz_xx[i] * tke_0 - 6.0 * tr_xyzz_x[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_y[i] = -2.0 * tr_yzz_xy[i] * tke_0 - 6.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_z[i] = -2.0 * tr_yzz_xz[i] * tke_0 - 6.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 27-30 components of targeted buffer : GP

    auto tr_x_0_x_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 27);

    auto tr_x_0_x_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 28);

    auto tr_x_0_x_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 29);

    #pragma omp simd aligned(tr_x_0_x_xzzz_x, tr_x_0_x_xzzz_y, tr_x_0_x_xzzz_z, tr_xxxzzz_x, tr_xxxzzz_y, tr_xxxzzz_z, tr_xxzzz_0, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xzzz_x[i] = tr_zzz_0[i] - 2.0 * tr_zzz_xx[i] * tke_0 - 6.0 * tr_xzzz_x[i] * tbe_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_y[i] = -2.0 * tr_zzz_xy[i] * tke_0 - 6.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_z[i] = -2.0 * tr_zzz_xz[i] * tke_0 - 6.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 30-33 components of targeted buffer : GP

    auto tr_x_0_x_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 30);

    auto tr_x_0_x_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 31);

    auto tr_x_0_x_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 32);

    #pragma omp simd aligned(tr_x_0_x_yyyy_x, tr_x_0_x_yyyy_y, tr_x_0_x_yyyy_z, tr_xxyyyy_x, tr_xxyyyy_y, tr_xxyyyy_z, tr_xyyyy_0, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyyy_x[i] = -2.0 * tr_yyyy_x[i] * tbe_0 - 2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_y[i] = -2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_z[i] = -2.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 33-36 components of targeted buffer : GP

    auto tr_x_0_x_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 33);

    auto tr_x_0_x_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 34);

    auto tr_x_0_x_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 35);

    #pragma omp simd aligned(tr_x_0_x_yyyz_x, tr_x_0_x_yyyz_y, tr_x_0_x_yyyz_z, tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xyyyz_0, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyyz_x[i] = -2.0 * tr_yyyz_x[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_y[i] = -2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_z[i] = -2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 36-39 components of targeted buffer : GP

    auto tr_x_0_x_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 36);

    auto tr_x_0_x_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 37);

    auto tr_x_0_x_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 38);

    #pragma omp simd aligned(tr_x_0_x_yyzz_x, tr_x_0_x_yyzz_y, tr_x_0_x_yyzz_z, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xyyzz_0, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyzz_x[i] = -2.0 * tr_yyzz_x[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_y[i] = -2.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_z[i] = -2.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 39-42 components of targeted buffer : GP

    auto tr_x_0_x_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 39);

    auto tr_x_0_x_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 40);

    auto tr_x_0_x_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 41);

    #pragma omp simd aligned(tr_x_0_x_yzzz_x, tr_x_0_x_yzzz_y, tr_x_0_x_yzzz_z, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_xyzzz_0, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yzzz_x[i] = -2.0 * tr_yzzz_x[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_y[i] = -2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_z[i] = -2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 42-45 components of targeted buffer : GP

    auto tr_x_0_x_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 42);

    auto tr_x_0_x_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 43);

    auto tr_x_0_x_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 44);

    #pragma omp simd aligned(tr_x_0_x_zzzz_x, tr_x_0_x_zzzz_y, tr_x_0_x_zzzz_z, tr_xxzzzz_x, tr_xxzzzz_y, tr_xxzzzz_z, tr_xzzzz_0, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zzzz_x[i] = -2.0 * tr_zzzz_x[i] * tbe_0 - 2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_y[i] = -2.0 * tr_zzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_z[i] = -2.0 * tr_zzzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 45-48 components of targeted buffer : GP

    auto tr_x_0_y_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 45);

    auto tr_x_0_y_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 46);

    auto tr_x_0_y_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 47);

    #pragma omp simd aligned(tr_x_0_y_xxxx_x, tr_x_0_y_xxxx_y, tr_x_0_y_xxxx_z, tr_xxx_0, tr_xxx_xy, tr_xxx_yy, tr_xxx_yz, tr_xxxxx_0, tr_xxxxx_xy, tr_xxxxx_yy, tr_xxxxx_yz, tr_xxxxxy_x, tr_xxxxxy_y, tr_xxxxxy_z, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxx_x[i] = -8.0 * tr_xxx_xy[i] * tke_0 - 8.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_y[i] = 4.0 * tr_xxx_0[i] - 8.0 * tr_xxx_yy[i] * tke_0 - 8.0 * tr_xxxy_y[i] * tbe_0 - 2.0 * tr_xxxxx_0[i] * tbe_0 + 4.0 * tr_xxxxx_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_z[i] = -8.0 * tr_xxx_yz[i] * tke_0 - 8.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 48-51 components of targeted buffer : GP

    auto tr_x_0_y_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 48);

    auto tr_x_0_y_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 49);

    auto tr_x_0_y_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 50);

    #pragma omp simd aligned(tr_x_0_y_xxxy_x, tr_x_0_y_xxxy_y, tr_x_0_y_xxxy_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxy_0, tr_xxxxy_xy, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxyy_x, tr_xxxxyy_y, tr_xxxxyy_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxy_x[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xxy_xy[i] * tke_0 - 6.0 * tr_xxyy_x[i] * tbe_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_y[i] = 3.0 * tr_xx_y[i] + 3.0 * tr_xxy_0[i] - 6.0 * tr_xxy_yy[i] * tke_0 - 6.0 * tr_xxyy_y[i] * tbe_0 - 2.0 * tr_xxxx_y[i] * tbe_0 - 2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_z[i] = 3.0 * tr_xx_z[i] - 6.0 * tr_xxy_yz[i] * tke_0 - 6.0 * tr_xxyy_z[i] * tbe_0 - 2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 51-54 components of targeted buffer : GP

    auto tr_x_0_y_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 51);

    auto tr_x_0_y_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 52);

    auto tr_x_0_y_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 53);

    #pragma omp simd aligned(tr_x_0_y_xxxz_x, tr_x_0_y_xxxz_y, tr_x_0_y_xxxz_z, tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxxxz_0, tr_xxxxz_xy, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxz_x[i] = -6.0 * tr_xxz_xy[i] * tke_0 - 6.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_y[i] = 3.0 * tr_xxz_0[i] - 6.0 * tr_xxz_yy[i] * tke_0 - 6.0 * tr_xxyz_y[i] * tbe_0 - 2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_z[i] = -6.0 * tr_xxz_yz[i] * tke_0 - 6.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 54-57 components of targeted buffer : GP

    auto tr_x_0_y_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 54);

    auto tr_x_0_y_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 55);

    auto tr_x_0_y_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 56);

    #pragma omp simd aligned(tr_x_0_y_xxyy_x, tr_x_0_y_xxyy_y, tr_x_0_y_xxyy_z, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyy_0, tr_xxxyy_xy, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyyy_x, tr_xxxyyy_y, tr_xxxyyy_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxyy_x[i] = 4.0 * tr_xy_x[i] - 4.0 * tr_xyy_xy[i] * tke_0 - 4.0 * tr_xyyy_x[i] * tbe_0 - 4.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_y[i] = 4.0 * tr_xy_y[i] + 2.0 * tr_xyy_0[i] - 4.0 * tr_xyy_yy[i] * tke_0 - 4.0 * tr_xyyy_y[i] * tbe_0 - 4.0 * tr_xxxy_y[i] * tbe_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_z[i] = 4.0 * tr_xy_z[i] - 4.0 * tr_xyy_yz[i] * tke_0 - 4.0 * tr_xyyy_z[i] * tbe_0 - 4.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 57-60 components of targeted buffer : GP

    auto tr_x_0_y_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 57);

    auto tr_x_0_y_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 58);

    auto tr_x_0_y_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 59);

    #pragma omp simd aligned(tr_x_0_y_xxyz_x, tr_x_0_y_xxyz_y, tr_x_0_y_xxyz_z, tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxxyz_0, tr_xxxyz_xy, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxyz_x[i] = 2.0 * tr_xz_x[i] - 4.0 * tr_xyz_xy[i] * tke_0 - 4.0 * tr_xyyz_x[i] * tbe_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_y[i] = 2.0 * tr_xz_y[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_yy[i] * tke_0 - 4.0 * tr_xyyz_y[i] * tbe_0 - 2.0 * tr_xxxz_y[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_z[i] = 2.0 * tr_xz_z[i] - 4.0 * tr_xyz_yz[i] * tke_0 - 4.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 60-63 components of targeted buffer : GP

    auto tr_x_0_y_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 60);

    auto tr_x_0_y_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 61);

    auto tr_x_0_y_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 62);

    #pragma omp simd aligned(tr_x_0_y_xxzz_x, tr_x_0_y_xxzz_y, tr_x_0_y_xxzz_z, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xxxzz_0, tr_xxxzz_xy, tr_xxxzz_yy, tr_xxxzz_yz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxzz_x[i] = -4.0 * tr_xzz_xy[i] * tke_0 - 4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_y[i] = 2.0 * tr_xzz_0[i] - 4.0 * tr_xzz_yy[i] * tke_0 - 4.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_z[i] = -4.0 * tr_xzz_yz[i] * tke_0 - 4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 63-66 components of targeted buffer : GP

    auto tr_x_0_y_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 63);

    auto tr_x_0_y_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 64);

    auto tr_x_0_y_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 65);

    #pragma omp simd aligned(tr_x_0_y_xyyy_x, tr_x_0_y_xyyy_y, tr_x_0_y_xyyy_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyy_0, tr_xxyyy_xy, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyyy_x, tr_xxyyyy_y, tr_xxyyyy_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xy, tr_yyy_yy, tr_yyy_yz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyyy_x[i] = 3.0 * tr_yy_x[i] - 2.0 * tr_yyy_xy[i] * tke_0 - 2.0 * tr_yyyy_x[i] * tbe_0 - 6.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_y[i] = 3.0 * tr_yy_y[i] + tr_yyy_0[i] - 2.0 * tr_yyy_yy[i] * tke_0 - 2.0 * tr_yyyy_y[i] * tbe_0 - 6.0 * tr_xxyy_y[i] * tbe_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_z[i] = 3.0 * tr_yy_z[i] - 2.0 * tr_yyy_yz[i] * tke_0 - 2.0 * tr_yyyy_z[i] * tbe_0 - 6.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 66-69 components of targeted buffer : GP

    auto tr_x_0_y_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 66);

    auto tr_x_0_y_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 67);

    auto tr_x_0_y_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 68);

    #pragma omp simd aligned(tr_x_0_y_xyyz_x, tr_x_0_y_xyyz_y, tr_x_0_y_xyyz_z, tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xxyyz_0, tr_xxyyz_xy, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyyz_x[i] = 2.0 * tr_yz_x[i] - 2.0 * tr_yyz_xy[i] * tke_0 - 2.0 * tr_yyyz_x[i] * tbe_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_y[i] = 2.0 * tr_yz_y[i] + tr_yyz_0[i] - 2.0 * tr_yyz_yy[i] * tke_0 - 2.0 * tr_yyyz_y[i] * tbe_0 - 4.0 * tr_xxyz_y[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_z[i] = 2.0 * tr_yz_z[i] - 2.0 * tr_yyz_yz[i] * tke_0 - 2.0 * tr_yyyz_z[i] * tbe_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 69-72 components of targeted buffer : GP

    auto tr_x_0_y_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 69);

    auto tr_x_0_y_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 70);

    auto tr_x_0_y_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 71);

    #pragma omp simd aligned(tr_x_0_y_xyzz_x, tr_x_0_y_xyzz_y, tr_x_0_y_xyzz_z, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xxyzz_0, tr_xxyzz_xy, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyzz_x[i] = tr_zz_x[i] - 2.0 * tr_yzz_xy[i] * tke_0 - 2.0 * tr_yyzz_x[i] * tbe_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_y[i] = tr_zz_y[i] + tr_yzz_0[i] - 2.0 * tr_yzz_yy[i] * tke_0 - 2.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_xxzz_y[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_z[i] = tr_zz_z[i] - 2.0 * tr_yzz_yz[i] * tke_0 - 2.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 72-75 components of targeted buffer : GP

    auto tr_x_0_y_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 72);

    auto tr_x_0_y_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 73);

    auto tr_x_0_y_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 74);

    #pragma omp simd aligned(tr_x_0_y_xzzz_x, tr_x_0_y_xzzz_y, tr_x_0_y_xzzz_z, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_xxzzz_0, tr_xxzzz_xy, tr_xxzzz_yy, tr_xxzzz_yz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_zzz_0, tr_zzz_xy, tr_zzz_yy, tr_zzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xzzz_x[i] = -2.0 * tr_zzz_xy[i] * tke_0 - 2.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_y[i] = tr_zzz_0[i] - 2.0 * tr_zzz_yy[i] * tke_0 - 2.0 * tr_yzzz_y[i] * tbe_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_z[i] = -2.0 * tr_zzz_yz[i] * tke_0 - 2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 75-78 components of targeted buffer : GP

    auto tr_x_0_y_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 75);

    auto tr_x_0_y_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 76);

    auto tr_x_0_y_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 77);

    #pragma omp simd aligned(tr_x_0_y_yyyy_x, tr_x_0_y_yyyy_y, tr_x_0_y_yyyy_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyy_0, tr_xyyyy_xy, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyyy_x, tr_xyyyyy_y, tr_xyyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyyy_x[i] = -8.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_y[i] = -8.0 * tr_xyyy_y[i] * tbe_0 - 2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_z[i] = -8.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 78-81 components of targeted buffer : GP

    auto tr_x_0_y_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 78);

    auto tr_x_0_y_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 79);

    auto tr_x_0_y_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 80);

    #pragma omp simd aligned(tr_x_0_y_yyyz_x, tr_x_0_y_yyyz_y, tr_x_0_y_yyyz_z, tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_xyyyz_0, tr_xyyyz_xy, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyyz_x[i] = -6.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_y[i] = -6.0 * tr_xyyz_y[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_z[i] = -6.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 81-84 components of targeted buffer : GP

    auto tr_x_0_y_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 81);

    auto tr_x_0_y_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 82);

    auto tr_x_0_y_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 83);

    #pragma omp simd aligned(tr_x_0_y_yyzz_x, tr_x_0_y_yyzz_y, tr_x_0_y_yyzz_z, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_xyyzz_0, tr_xyyzz_xy, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyzz_x[i] = -4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_y[i] = -4.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_z[i] = -4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 84-87 components of targeted buffer : GP

    auto tr_x_0_y_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 84);

    auto tr_x_0_y_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 85);

    auto tr_x_0_y_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 86);

    #pragma omp simd aligned(tr_x_0_y_yzzz_x, tr_x_0_y_yzzz_y, tr_x_0_y_yzzz_z, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_xyzzz_0, tr_xyzzz_xy, tr_xyzzz_yy, tr_xyzzz_yz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yzzz_x[i] = -2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_y[i] = -2.0 * tr_xzzz_y[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_z[i] = -2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 87-90 components of targeted buffer : GP

    auto tr_x_0_y_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 87);

    auto tr_x_0_y_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 88);

    auto tr_x_0_y_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 89);

    #pragma omp simd aligned(tr_x_0_y_zzzz_x, tr_x_0_y_zzzz_y, tr_x_0_y_zzzz_z, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_xzzzz_0, tr_xzzzz_xy, tr_xzzzz_yy, tr_xzzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zzzz_x[i] = 4.0 * tr_xzzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_y[i] = -2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_z[i] = 4.0 * tr_xzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 90-93 components of targeted buffer : GP

    auto tr_x_0_z_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 90);

    auto tr_x_0_z_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 91);

    auto tr_x_0_z_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 92);

    #pragma omp simd aligned(tr_x_0_z_xxxx_x, tr_x_0_z_xxxx_y, tr_x_0_z_xxxx_z, tr_xxx_0, tr_xxx_xz, tr_xxx_yz, tr_xxx_zz, tr_xxxxx_0, tr_xxxxx_xz, tr_xxxxx_yz, tr_xxxxx_zz, tr_xxxxxz_x, tr_xxxxxz_y, tr_xxxxxz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxx_x[i] = -8.0 * tr_xxx_xz[i] * tke_0 - 8.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_y[i] = -8.0 * tr_xxx_yz[i] * tke_0 - 8.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_z[i] = 4.0 * tr_xxx_0[i] - 8.0 * tr_xxx_zz[i] * tke_0 - 8.0 * tr_xxxz_z[i] * tbe_0 - 2.0 * tr_xxxxx_0[i] * tbe_0 + 4.0 * tr_xxxxx_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 93-96 components of targeted buffer : GP

    auto tr_x_0_z_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 93);

    auto tr_x_0_z_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 94);

    auto tr_x_0_z_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 95);

    #pragma omp simd aligned(tr_x_0_z_xxxy_x, tr_x_0_z_xxxy_y, tr_x_0_z_xxxy_z, tr_xxxxy_0, tr_xxxxy_xz, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxy_x[i] = -6.0 * tr_xxy_xz[i] * tke_0 - 6.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_y[i] = -6.0 * tr_xxy_yz[i] * tke_0 - 6.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_z[i] = 3.0 * tr_xxy_0[i] - 6.0 * tr_xxy_zz[i] * tke_0 - 6.0 * tr_xxyz_z[i] * tbe_0 - 2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 96-99 components of targeted buffer : GP

    auto tr_x_0_z_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 96);

    auto tr_x_0_z_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 97);

    auto tr_x_0_z_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 98);

    #pragma omp simd aligned(tr_x_0_z_xxxz_x, tr_x_0_z_xxxz_y, tr_x_0_z_xxxz_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxz_0, tr_xxxxz_xz, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxxzz_x, tr_xxxxzz_y, tr_xxxxzz_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxz_x[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xxz_xz[i] * tke_0 - 6.0 * tr_xxzz_x[i] * tbe_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_y[i] = 3.0 * tr_xx_y[i] - 6.0 * tr_xxz_yz[i] * tke_0 - 6.0 * tr_xxzz_y[i] * tbe_0 - 2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_z[i] = 3.0 * tr_xx_z[i] + 3.0 * tr_xxz_0[i] - 6.0 * tr_xxz_zz[i] * tke_0 - 6.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxxx_z[i] * tbe_0 - 2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 99-102 components of targeted buffer : GP

    auto tr_x_0_z_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 99);

    auto tr_x_0_z_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 100);

    auto tr_x_0_z_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 101);

    #pragma omp simd aligned(tr_x_0_z_xxyy_x, tr_x_0_z_xxyy_y, tr_x_0_z_xxyy_z, tr_xxxyy_0, tr_xxxyy_xz, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxyy_x[i] = -4.0 * tr_xyy_xz[i] * tke_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_y[i] = -4.0 * tr_xyy_yz[i] * tke_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_z[i] = 2.0 * tr_xyy_0[i] - 4.0 * tr_xyy_zz[i] * tke_0 - 4.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 102-105 components of targeted buffer : GP

    auto tr_x_0_z_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 102);

    auto tr_x_0_z_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 103);

    auto tr_x_0_z_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 104);

    #pragma omp simd aligned(tr_x_0_z_xxyz_x, tr_x_0_z_xxyz_y, tr_x_0_z_xxyz_z, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyz_0, tr_xxxyz_xz, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxyz_x[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xyz_xz[i] * tke_0 - 4.0 * tr_xyzz_x[i] * tbe_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_y[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xyz_yz[i] * tke_0 - 4.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_z[i] = 2.0 * tr_xy_z[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_zz[i] * tke_0 - 4.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xxxy_z[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 105-108 components of targeted buffer : GP

    auto tr_x_0_z_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 105);

    auto tr_x_0_z_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 106);

    auto tr_x_0_z_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 107);

    #pragma omp simd aligned(tr_x_0_z_xxzz_x, tr_x_0_z_xxzz_y, tr_x_0_z_xxzz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxxzz_0, tr_xxxzz_xz, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxxzzz_x, tr_xxxzzz_y, tr_xxxzzz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxzz_x[i] = 4.0 * tr_xz_x[i] - 4.0 * tr_xzz_xz[i] * tke_0 - 4.0 * tr_xzzz_x[i] * tbe_0 - 4.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_y[i] = 4.0 * tr_xz_y[i] - 4.0 * tr_xzz_yz[i] * tke_0 - 4.0 * tr_xzzz_y[i] * tbe_0 - 4.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_z[i] = 4.0 * tr_xz_z[i] + 2.0 * tr_xzz_0[i] - 4.0 * tr_xzz_zz[i] * tke_0 - 4.0 * tr_xzzz_z[i] * tbe_0 - 4.0 * tr_xxxz_z[i] * tbe_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 108-111 components of targeted buffer : GP

    auto tr_x_0_z_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 108);

    auto tr_x_0_z_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 109);

    auto tr_x_0_z_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 110);

    #pragma omp simd aligned(tr_x_0_z_xyyy_x, tr_x_0_z_xyyy_y, tr_x_0_z_xyyy_z, tr_xxyyy_0, tr_xxyyy_xz, tr_xxyyy_yz, tr_xxyyy_zz, tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_yyy_0, tr_yyy_xz, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyyy_x[i] = -2.0 * tr_yyy_xz[i] * tke_0 - 2.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_y[i] = -2.0 * tr_yyy_yz[i] * tke_0 - 2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_z[i] = tr_yyy_0[i] - 2.0 * tr_yyy_zz[i] * tke_0 - 2.0 * tr_yyyz_z[i] * tbe_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 111-114 components of targeted buffer : GP

    auto tr_x_0_z_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 111);

    auto tr_x_0_z_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 112);

    auto tr_x_0_z_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 113);

    #pragma omp simd aligned(tr_x_0_z_xyyz_x, tr_x_0_z_xyyz_y, tr_x_0_z_xyyz_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyz_0, tr_xxyyz_xz, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyyz_x[i] = tr_yy_x[i] - 2.0 * tr_yyz_xz[i] * tke_0 - 2.0 * tr_yyzz_x[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_y[i] = tr_yy_y[i] - 2.0 * tr_yyz_yz[i] * tke_0 - 2.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_z[i] = tr_yy_z[i] + tr_yyz_0[i] - 2.0 * tr_yyz_zz[i] * tke_0 - 2.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_xxyy_z[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 114-117 components of targeted buffer : GP

    auto tr_x_0_z_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 114);

    auto tr_x_0_z_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 115);

    auto tr_x_0_z_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 116);

    #pragma omp simd aligned(tr_x_0_z_xyzz_x, tr_x_0_z_xyzz_y, tr_x_0_z_xyzz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzz_0, tr_xxyzz_xz, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyzz_x[i] = 2.0 * tr_yz_x[i] - 2.0 * tr_yzz_xz[i] * tke_0 - 2.0 * tr_yzzz_x[i] * tbe_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_y[i] = 2.0 * tr_yz_y[i] - 2.0 * tr_yzz_yz[i] * tke_0 - 2.0 * tr_yzzz_y[i] * tbe_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_z[i] = 2.0 * tr_yz_z[i] + tr_yzz_0[i] - 2.0 * tr_yzz_zz[i] * tke_0 - 2.0 * tr_yzzz_z[i] * tbe_0 - 4.0 * tr_xxyz_z[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 117-120 components of targeted buffer : GP

    auto tr_x_0_z_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 117);

    auto tr_x_0_z_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 118);

    auto tr_x_0_z_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 119);

    #pragma omp simd aligned(tr_x_0_z_xzzz_x, tr_x_0_z_xzzz_y, tr_x_0_z_xzzz_z, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xxzzz_0, tr_xxzzz_xz, tr_xxzzz_yz, tr_xxzzz_zz, tr_xxzzzz_x, tr_xxzzzz_y, tr_xxzzzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xz, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xzzz_x[i] = 3.0 * tr_zz_x[i] - 2.0 * tr_zzz_xz[i] * tke_0 - 2.0 * tr_zzzz_x[i] * tbe_0 - 6.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_y[i] = 3.0 * tr_zz_y[i] - 2.0 * tr_zzz_yz[i] * tke_0 - 2.0 * tr_zzzz_y[i] * tbe_0 - 6.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_z[i] = 3.0 * tr_zz_z[i] + tr_zzz_0[i] - 2.0 * tr_zzz_zz[i] * tke_0 - 2.0 * tr_zzzz_z[i] * tbe_0 - 6.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 120-123 components of targeted buffer : GP

    auto tr_x_0_z_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 120);

    auto tr_x_0_z_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 121);

    auto tr_x_0_z_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 122);

    #pragma omp simd aligned(tr_x_0_z_yyyy_x, tr_x_0_z_yyyy_y, tr_x_0_z_yyyy_z, tr_xyyyy_0, tr_xyyyy_xz, tr_xyyyy_yz, tr_xyyyy_zz, tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyyy_x[i] = 4.0 * tr_xyyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_y[i] = 4.0 * tr_xyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_z[i] = -2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 123-126 components of targeted buffer : GP

    auto tr_x_0_z_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 123);

    auto tr_x_0_z_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 124);

    auto tr_x_0_z_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 125);

    #pragma omp simd aligned(tr_x_0_z_yyyz_x, tr_x_0_z_yyyz_y, tr_x_0_z_yyyz_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyz_0, tr_xyyyz_xz, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyyz_x[i] = -2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_y[i] = -2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_z[i] = -2.0 * tr_xyyy_z[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 126-129 components of targeted buffer : GP

    auto tr_x_0_z_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 126);

    auto tr_x_0_z_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 127);

    auto tr_x_0_z_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 128);

    #pragma omp simd aligned(tr_x_0_z_yyzz_x, tr_x_0_z_yyzz_y, tr_x_0_z_yyzz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzz_0, tr_xyyzz_xz, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyzz_x[i] = -4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_y[i] = -4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_z[i] = -4.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 129-132 components of targeted buffer : GP

    auto tr_x_0_z_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 129);

    auto tr_x_0_z_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 130);

    auto tr_x_0_z_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 131);

    #pragma omp simd aligned(tr_x_0_z_yzzz_x, tr_x_0_z_yzzz_y, tr_x_0_z_yzzz_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzz_0, tr_xyzzz_xz, tr_xyzzz_yz, tr_xyzzz_zz, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yzzz_x[i] = -6.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_y[i] = -6.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_z[i] = -6.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 132-135 components of targeted buffer : GP

    auto tr_x_0_z_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 132);

    auto tr_x_0_z_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 133);

    auto tr_x_0_z_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 134);

    #pragma omp simd aligned(tr_x_0_z_zzzz_x, tr_x_0_z_zzzz_y, tr_x_0_z_zzzz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_xzzzz_0, tr_xzzzz_xz, tr_xzzzz_yz, tr_xzzzz_zz, tr_xzzzzz_x, tr_xzzzzz_y, tr_xzzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zzzz_x[i] = -8.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_y[i] = -8.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_z[i] = -8.0 * tr_xzzz_z[i] * tbe_0 - 2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 135-138 components of targeted buffer : GP

    auto tr_y_0_x_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 135);

    auto tr_y_0_x_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 136);

    auto tr_y_0_x_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 137);

    #pragma omp simd aligned(tr_xxxxxy_x, tr_xxxxxy_y, tr_xxxxxy_z, tr_xxxxy_0, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_y_0_x_xxxx_x, tr_y_0_x_xxxx_y, tr_y_0_x_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxx_x[i] = -8.0 * tr_xxxy_x[i] * tbe_0 - 2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_y[i] = -8.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_z[i] = -8.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 138-141 components of targeted buffer : GP

    auto tr_y_0_x_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 138);

    auto tr_y_0_x_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 139);

    auto tr_y_0_x_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 140);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxyy_x, tr_xxxxyy_y, tr_xxxxyy_z, tr_xxxyy_0, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_y_0_x_xxxy_x, tr_y_0_x_xxxy_y, tr_y_0_x_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxy_x[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xxyy_x[i] * tbe_0 + tr_xxx_0[i] - 2.0 * tr_xxx_xx[i] * tke_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_y[i] = 3.0 * tr_xx_y[i] - 6.0 * tr_xxyy_y[i] * tbe_0 - 2.0 * tr_xxx_xy[i] * tke_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_z[i] = 3.0 * tr_xx_z[i] - 6.0 * tr_xxyy_z[i] * tbe_0 - 2.0 * tr_xxx_xz[i] * tke_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 141-144 components of targeted buffer : GP

    auto tr_y_0_x_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 141);

    auto tr_y_0_x_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 142);

    auto tr_y_0_x_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 143);

    #pragma omp simd aligned(tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxxyz_0, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_y_0_x_xxxz_x, tr_y_0_x_xxxz_y, tr_y_0_x_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxz_x[i] = -6.0 * tr_xxyz_x[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_y[i] = -6.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_z[i] = -6.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 144-147 components of targeted buffer : GP

    auto tr_y_0_x_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 144);

    auto tr_y_0_x_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 145);

    auto tr_y_0_x_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 146);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyyy_x, tr_xxxyyy_y, tr_xxxyyy_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxyyy_0, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_y_0_x_xxyy_x, tr_y_0_x_xxyy_y, tr_y_0_x_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxyy_x[i] = 4.0 * tr_xy_x[i] - 4.0 * tr_xyyy_x[i] * tbe_0 + 2.0 * tr_xxy_0[i] - 4.0 * tr_xxy_xx[i] * tke_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_y[i] = 4.0 * tr_xy_y[i] - 4.0 * tr_xyyy_y[i] * tbe_0 - 4.0 * tr_xxy_xy[i] * tke_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_z[i] = 4.0 * tr_xy_z[i] - 4.0 * tr_xyyy_z[i] * tbe_0 - 4.0 * tr_xxy_xz[i] * tke_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 147-150 components of targeted buffer : GP

    auto tr_y_0_x_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 147);

    auto tr_y_0_x_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 148);

    auto tr_y_0_x_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 149);

    #pragma omp simd aligned(tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxyyz_0, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0_x_xxyz_x, tr_y_0_x_xxyz_y, tr_y_0_x_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxyz_x[i] = 2.0 * tr_xz_x[i] - 4.0 * tr_xyyz_x[i] * tbe_0 + tr_xxz_0[i] - 2.0 * tr_xxz_xx[i] * tke_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_y[i] = 2.0 * tr_xz_y[i] - 4.0 * tr_xyyz_y[i] * tbe_0 - 2.0 * tr_xxz_xy[i] * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_z[i] = 2.0 * tr_xz_z[i] - 4.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xxz_xz[i] * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 150-153 components of targeted buffer : GP

    auto tr_y_0_x_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 150);

    auto tr_y_0_x_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 151);

    auto tr_y_0_x_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 152);

    #pragma omp simd aligned(tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xxyzz_0, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0_x_xxzz_x, tr_y_0_x_xxzz_y, tr_y_0_x_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxzz_x[i] = -4.0 * tr_xyzz_x[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_y[i] = -4.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_z[i] = -4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 153-156 components of targeted buffer : GP

    auto tr_y_0_x_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 153);

    auto tr_y_0_x_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 154);

    auto tr_y_0_x_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 155);

    #pragma omp simd aligned(tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyyy_x, tr_xxyyyy_y, tr_xxyyyy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyyyy_0, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_y_0_x_xyyy_x, tr_y_0_x_xyyy_y, tr_y_0_x_xyyy_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyyy_x[i] = 3.0 * tr_yy_x[i] - 2.0 * tr_yyyy_x[i] * tbe_0 + 3.0 * tr_xyy_0[i] - 6.0 * tr_xyy_xx[i] * tke_0 - 2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_y[i] = 3.0 * tr_yy_y[i] - 2.0 * tr_yyyy_y[i] * tbe_0 - 6.0 * tr_xyy_xy[i] * tke_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_z[i] = 3.0 * tr_yy_z[i] - 2.0 * tr_yyyy_z[i] * tbe_0 - 6.0 * tr_xyy_xz[i] * tke_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 156-159 components of targeted buffer : GP

    auto tr_y_0_x_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 156);

    auto tr_y_0_x_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 157);

    auto tr_y_0_x_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 158);

    #pragma omp simd aligned(tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyyyz_0, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_y_0_x_xyyz_x, tr_y_0_x_xyyz_y, tr_y_0_x_xyyz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyyz_x[i] = 2.0 * tr_yz_x[i] - 2.0 * tr_yyyz_x[i] * tbe_0 + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_xx[i] * tke_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_y[i] = 2.0 * tr_yz_y[i] - 2.0 * tr_yyyz_y[i] * tbe_0 - 4.0 * tr_xyz_xy[i] * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_z[i] = 2.0 * tr_yz_z[i] - 2.0 * tr_yyyz_z[i] * tbe_0 - 4.0 * tr_xyz_xz[i] * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 159-162 components of targeted buffer : GP

    auto tr_y_0_x_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 159);

    auto tr_y_0_x_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 160);

    auto tr_y_0_x_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 161);

    #pragma omp simd aligned(tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xyyzz_0, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_y_0_x_xyzz_x, tr_y_0_x_xyzz_y, tr_y_0_x_xyzz_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyzz_x[i] = tr_zz_x[i] - 2.0 * tr_yyzz_x[i] * tbe_0 + tr_xzz_0[i] - 2.0 * tr_xzz_xx[i] * tke_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_y[i] = tr_zz_y[i] - 2.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_xzz_xy[i] * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_z[i] = tr_zz_z[i] - 2.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_xzz_xz[i] * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 162-165 components of targeted buffer : GP

    auto tr_y_0_x_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 162);

    auto tr_y_0_x_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 163);

    auto tr_y_0_x_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 164);

    #pragma omp simd aligned(tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_xyzzz_0, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_y_0_x_xzzz_x, tr_y_0_x_xzzz_y, tr_y_0_x_xzzz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xzzz_x[i] = -2.0 * tr_yzzz_x[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_y[i] = -2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_z[i] = -2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 165-168 components of targeted buffer : GP

    auto tr_y_0_x_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 165);

    auto tr_y_0_x_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 166);

    auto tr_y_0_x_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 167);

    #pragma omp simd aligned(tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyyy_x, tr_xyyyyy_y, tr_xyyyyy_z, tr_y_0_x_yyyy_x, tr_y_0_x_yyyy_y, tr_y_0_x_yyyy_z, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyyyy_0, tr_yyyyy_xx, tr_yyyyy_xy, tr_yyyyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyyy_x[i] = 4.0 * tr_yyy_0[i] - 8.0 * tr_yyy_xx[i] * tke_0 - 2.0 * tr_yyyyy_0[i] * tbe_0 + 4.0 * tr_yyyyy_xx[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_y[i] = -8.0 * tr_yyy_xy[i] * tke_0 + 4.0 * tr_yyyyy_xy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_z[i] = -8.0 * tr_yyy_xz[i] * tke_0 + 4.0 * tr_yyyyy_xz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 168-171 components of targeted buffer : GP

    auto tr_y_0_x_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 168);

    auto tr_y_0_x_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 169);

    auto tr_y_0_x_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 170);

    #pragma omp simd aligned(tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_y_0_x_yyyz_x, tr_y_0_x_yyyz_y, tr_y_0_x_yyyz_z, tr_yyyyz_0, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyyz_x[i] = 3.0 * tr_yyz_0[i] - 6.0 * tr_yyz_xx[i] * tke_0 - 2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_y[i] = -6.0 * tr_yyz_xy[i] * tke_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_z[i] = -6.0 * tr_yyz_xz[i] * tke_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 171-174 components of targeted buffer : GP

    auto tr_y_0_x_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 171);

    auto tr_y_0_x_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 172);

    auto tr_y_0_x_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 173);

    #pragma omp simd aligned(tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0_x_yyzz_x, tr_y_0_x_yyzz_y, tr_y_0_x_yyzz_z, tr_yyyzz_0, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyzz_x[i] = 2.0 * tr_yzz_0[i] - 4.0 * tr_yzz_xx[i] * tke_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_y[i] = -4.0 * tr_yzz_xy[i] * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_z[i] = -4.0 * tr_yzz_xz[i] * tke_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 174-177 components of targeted buffer : GP

    auto tr_y_0_x_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 174);

    auto tr_y_0_x_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 175);

    auto tr_y_0_x_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 176);

    #pragma omp simd aligned(tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_y_0_x_yzzz_x, tr_y_0_x_yzzz_y, tr_y_0_x_yzzz_z, tr_yyzzz_0, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yzzz_x[i] = tr_zzz_0[i] - 2.0 * tr_zzz_xx[i] * tke_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_y[i] = -2.0 * tr_zzz_xy[i] * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_z[i] = -2.0 * tr_zzz_xz[i] * tke_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 177-180 components of targeted buffer : GP

    auto tr_y_0_x_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 177);

    auto tr_y_0_x_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 178);

    auto tr_y_0_x_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 179);

    #pragma omp simd aligned(tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_y_0_x_zzzz_x, tr_y_0_x_zzzz_y, tr_y_0_x_zzzz_z, tr_yzzzz_0, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zzzz_x[i] = -2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_y[i] = 4.0 * tr_yzzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_z[i] = 4.0 * tr_yzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 180-183 components of targeted buffer : GP

    auto tr_y_0_y_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 180);

    auto tr_y_0_y_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 181);

    auto tr_y_0_y_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 182);

    #pragma omp simd aligned(tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxy_0, tr_xxxxy_xy, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxyy_x, tr_xxxxyy_y, tr_xxxxyy_z, tr_y_0_y_xxxx_x, tr_y_0_y_xxxx_y, tr_y_0_y_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxx_x[i] = -2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_y[i] = -2.0 * tr_xxxx_y[i] * tbe_0 - 2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_z[i] = -2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 183-186 components of targeted buffer : GP

    auto tr_y_0_y_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 183);

    auto tr_y_0_y_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 184);

    auto tr_y_0_y_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 185);

    #pragma omp simd aligned(tr_xxx_0, tr_xxx_xy, tr_xxx_yy, tr_xxx_yz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyy_0, tr_xxxyy_xy, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyyy_x, tr_xxxyyy_y, tr_xxxyyy_z, tr_y_0_y_xxxy_x, tr_y_0_y_xxxy_y, tr_y_0_y_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxy_x[i] = -2.0 * tr_xxx_xy[i] * tke_0 - 6.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_y[i] = tr_xxx_0[i] - 2.0 * tr_xxx_yy[i] * tke_0 - 6.0 * tr_xxxy_y[i] * tbe_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_z[i] = -2.0 * tr_xxx_yz[i] * tke_0 - 6.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 186-189 components of targeted buffer : GP

    auto tr_y_0_y_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 186);

    auto tr_y_0_y_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 187);

    auto tr_y_0_y_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 188);

    #pragma omp simd aligned(tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxxyz_0, tr_xxxyz_xy, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_y_0_y_xxxz_x, tr_y_0_y_xxxz_y, tr_y_0_y_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxz_x[i] = -2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_y[i] = -2.0 * tr_xxxz_y[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_z[i] = -2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 189-192 components of targeted buffer : GP

    auto tr_y_0_y_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 189);

    auto tr_y_0_y_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 190);

    auto tr_y_0_y_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 191);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyy_0, tr_xxyyy_xy, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyyy_x, tr_xxyyyy_y, tr_xxyyyy_z, tr_y_0_y_xxyy_x, tr_y_0_y_xxyy_y, tr_y_0_y_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxyy_x[i] = 2.0 * tr_xx_x[i] - 4.0 * tr_xxy_xy[i] * tke_0 - 10.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_y[i] = 2.0 * tr_xx_y[i] + 2.0 * tr_xxy_0[i] - 4.0 * tr_xxy_yy[i] * tke_0 - 10.0 * tr_xxyy_y[i] * tbe_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_z[i] = 2.0 * tr_xx_z[i] - 4.0 * tr_xxy_yz[i] * tke_0 - 10.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 192-195 components of targeted buffer : GP

    auto tr_y_0_y_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 192);

    auto tr_y_0_y_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 193);

    auto tr_y_0_y_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 194);

    #pragma omp simd aligned(tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xxyyz_0, tr_xxyyz_xy, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, tr_y_0_y_xxyz_x, tr_y_0_y_xxyz_y, tr_y_0_y_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxyz_x[i] = -2.0 * tr_xxz_xy[i] * tke_0 - 6.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_y[i] = tr_xxz_0[i] - 2.0 * tr_xxz_yy[i] * tke_0 - 6.0 * tr_xxyz_y[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_z[i] = -2.0 * tr_xxz_yz[i] * tke_0 - 6.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 195-198 components of targeted buffer : GP

    auto tr_y_0_y_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 195);

    auto tr_y_0_y_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 196);

    auto tr_y_0_y_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 197);

    #pragma omp simd aligned(tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xxyzz_0, tr_xxyzz_xy, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_y_0_y_xxzz_x, tr_y_0_y_xxzz_y, tr_y_0_y_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxzz_x[i] = -2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_y[i] = -2.0 * tr_xxzz_y[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_z[i] = -2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 198-201 components of targeted buffer : GP

    auto tr_y_0_y_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 198);

    auto tr_y_0_y_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 199);

    auto tr_y_0_y_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 200);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyy_0, tr_xyyyy_xy, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyyy_x, tr_xyyyyy_y, tr_xyyyyy_z, tr_y_0_y_xyyy_x, tr_y_0_y_xyyy_y, tr_y_0_y_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyyy_x[i] = 6.0 * tr_xy_x[i] - 6.0 * tr_xyy_xy[i] * tke_0 - 14.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_y[i] = 6.0 * tr_xy_y[i] + 3.0 * tr_xyy_0[i] - 6.0 * tr_xyy_yy[i] * tke_0 - 14.0 * tr_xyyy_y[i] * tbe_0 - 2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_z[i] = 6.0 * tr_xy_z[i] - 6.0 * tr_xyy_yz[i] * tke_0 - 14.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 201-204 components of targeted buffer : GP

    auto tr_y_0_y_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 201);

    auto tr_y_0_y_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 202);

    auto tr_y_0_y_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 203);

    #pragma omp simd aligned(tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_xyyyz_0, tr_xyyyz_xy, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0_y_xyyz_x, tr_y_0_y_xyyz_y, tr_y_0_y_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyyz_x[i] = 2.0 * tr_xz_x[i] - 4.0 * tr_xyz_xy[i] * tke_0 - 10.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_y[i] = 2.0 * tr_xz_y[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_yy[i] * tke_0 - 10.0 * tr_xyyz_y[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_z[i] = 2.0 * tr_xz_z[i] - 4.0 * tr_xyz_yz[i] * tke_0 - 10.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 204-207 components of targeted buffer : GP

    auto tr_y_0_y_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 204);

    auto tr_y_0_y_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 205);

    auto tr_y_0_y_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 206);

    #pragma omp simd aligned(tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_xyyzz_0, tr_xyyzz_xy, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, tr_y_0_y_xyzz_x, tr_y_0_y_xyzz_y, tr_y_0_y_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyzz_x[i] = -2.0 * tr_xzz_xy[i] * tke_0 - 6.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_y[i] = tr_xzz_0[i] - 2.0 * tr_xzz_yy[i] * tke_0 - 6.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_z[i] = -2.0 * tr_xzz_yz[i] * tke_0 - 6.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 207-210 components of targeted buffer : GP

    auto tr_y_0_y_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 207);

    auto tr_y_0_y_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 208);

    auto tr_y_0_y_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 209);

    #pragma omp simd aligned(tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_xyzzz_0, tr_xyzzz_xy, tr_xyzzz_yy, tr_xyzzz_yz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_y_0_y_xzzz_x, tr_y_0_y_xzzz_y, tr_y_0_y_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xzzz_x[i] = -2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_y[i] = -2.0 * tr_xzzz_y[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_z[i] = -2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 210-213 components of targeted buffer : GP

    auto tr_y_0_y_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 210);

    auto tr_y_0_y_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 211);

    auto tr_y_0_y_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 212);

    #pragma omp simd aligned(tr_y_0_y_yyyy_x, tr_y_0_y_yyyy_y, tr_y_0_y_yyyy_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xy, tr_yyy_yy, tr_yyy_yz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyy_0, tr_yyyyy_xy, tr_yyyyy_yy, tr_yyyyy_yz, tr_yyyyyy_x, tr_yyyyyy_y, tr_yyyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyyy_x[i] = 12.0 * tr_yy_x[i] - 8.0 * tr_yyy_xy[i] * tke_0 - 18.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_y[i] = 12.0 * tr_yy_y[i] + 4.0 * tr_yyy_0[i] - 8.0 * tr_yyy_yy[i] * tke_0 - 18.0 * tr_yyyy_y[i] * tbe_0 - 2.0 * tr_yyyyy_0[i] * tbe_0 + 4.0 * tr_yyyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_z[i] = 12.0 * tr_yy_z[i] - 8.0 * tr_yyy_yz[i] * tke_0 - 18.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_yyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 213-216 components of targeted buffer : GP

    auto tr_y_0_y_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 213);

    auto tr_y_0_y_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 214);

    auto tr_y_0_y_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 215);

    #pragma omp simd aligned(tr_y_0_y_yyyz_x, tr_y_0_y_yyyz_y, tr_y_0_y_yyyz_z, tr_yyyyyz_x, tr_yyyyyz_y, tr_yyyyyz_z, tr_yyyyz_0, tr_yyyyz_xy, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyyz_x[i] = 6.0 * tr_yz_x[i] - 6.0 * tr_yyz_xy[i] * tke_0 - 14.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_y[i] = 6.0 * tr_yz_y[i] + 3.0 * tr_yyz_0[i] - 6.0 * tr_yyz_yy[i] * tke_0 - 14.0 * tr_yyyz_y[i] * tbe_0 - 2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_z[i] = 6.0 * tr_yz_z[i] - 6.0 * tr_yyz_yz[i] * tke_0 - 14.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 216-219 components of targeted buffer : GP

    auto tr_y_0_y_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 216);

    auto tr_y_0_y_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 217);

    auto tr_y_0_y_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 218);

    #pragma omp simd aligned(tr_y_0_y_yyzz_x, tr_y_0_y_yyzz_y, tr_y_0_y_yyzz_z, tr_yyyyzz_x, tr_yyyyzz_y, tr_yyyyzz_z, tr_yyyzz_0, tr_yyyzz_xy, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyzz_x[i] = 2.0 * tr_zz_x[i] - 4.0 * tr_yzz_xy[i] * tke_0 - 10.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_y[i] = 2.0 * tr_zz_y[i] + 2.0 * tr_yzz_0[i] - 4.0 * tr_yzz_yy[i] * tke_0 - 10.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_z[i] = 2.0 * tr_zz_z[i] - 4.0 * tr_yzz_yz[i] * tke_0 - 10.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 219-222 components of targeted buffer : GP

    auto tr_y_0_y_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 219);

    auto tr_y_0_y_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 220);

    auto tr_y_0_y_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 221);

    #pragma omp simd aligned(tr_y_0_y_yzzz_x, tr_y_0_y_yzzz_y, tr_y_0_y_yzzz_z, tr_yyyzzz_x, tr_yyyzzz_y, tr_yyyzzz_z, tr_yyzzz_0, tr_yyzzz_xy, tr_yyzzz_yy, tr_yyzzz_yz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_zzz_0, tr_zzz_xy, tr_zzz_yy, tr_zzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yzzz_x[i] = -2.0 * tr_zzz_xy[i] * tke_0 - 6.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_y[i] = tr_zzz_0[i] - 2.0 * tr_zzz_yy[i] * tke_0 - 6.0 * tr_yzzz_y[i] * tbe_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_z[i] = -2.0 * tr_zzz_yz[i] * tke_0 - 6.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 222-225 components of targeted buffer : GP

    auto tr_y_0_y_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 222);

    auto tr_y_0_y_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 223);

    auto tr_y_0_y_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 224);

    #pragma omp simd aligned(tr_y_0_y_zzzz_x, tr_y_0_y_zzzz_y, tr_y_0_y_zzzz_z, tr_yyzzzz_x, tr_yyzzzz_y, tr_yyzzzz_z, tr_yzzzz_0, tr_yzzzz_xy, tr_yzzzz_yy, tr_yzzzz_yz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zzzz_x[i] = -2.0 * tr_zzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_y[i] = -2.0 * tr_zzzz_y[i] * tbe_0 - 2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_z[i] = -2.0 * tr_zzzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 225-228 components of targeted buffer : GP

    auto tr_y_0_z_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 225);

    auto tr_y_0_z_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 226);

    auto tr_y_0_z_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 227);

    #pragma omp simd aligned(tr_xxxxy_0, tr_xxxxy_xz, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_y_0_z_xxxx_x, tr_y_0_z_xxxx_y, tr_y_0_z_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxx_x[i] = 4.0 * tr_xxxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_y[i] = 4.0 * tr_xxxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_z[i] = -2.0 * tr_xxxxy_0[i] * tbe_0 + 4.0 * tr_xxxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 228-231 components of targeted buffer : GP

    auto tr_y_0_z_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 228);

    auto tr_y_0_z_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 229);

    auto tr_y_0_z_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 230);

    #pragma omp simd aligned(tr_xxx_0, tr_xxx_xz, tr_xxx_yz, tr_xxx_zz, tr_xxxyy_0, tr_xxxyy_xz, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_y_0_z_xxxy_x, tr_y_0_z_xxxy_y, tr_y_0_z_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxy_x[i] = -2.0 * tr_xxx_xz[i] * tke_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_y[i] = -2.0 * tr_xxx_yz[i] * tke_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_z[i] = tr_xxx_0[i] - 2.0 * tr_xxx_zz[i] * tke_0 - 2.0 * tr_xxxz_z[i] * tbe_0 - 2.0 * tr_xxxyy_0[i] * tbe_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 231-234 components of targeted buffer : GP

    auto tr_y_0_z_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 231);

    auto tr_y_0_z_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 232);

    auto tr_y_0_z_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 233);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyz_0, tr_xxxyz_xz, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_y_0_z_xxxz_x, tr_y_0_z_xxxz_y, tr_y_0_z_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxz_x[i] = -2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_y[i] = -2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_z[i] = -2.0 * tr_xxxy_z[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 234-237 components of targeted buffer : GP

    auto tr_y_0_z_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 234);

    auto tr_y_0_z_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 235);

    auto tr_y_0_z_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 236);

    #pragma omp simd aligned(tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyyy_0, tr_xxyyy_xz, tr_xxyyy_yz, tr_xxyyy_zz, tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_y_0_z_xxyy_x, tr_y_0_z_xxyy_y, tr_y_0_z_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxyy_x[i] = -4.0 * tr_xxy_xz[i] * tke_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_y[i] = -4.0 * tr_xxy_yz[i] * tke_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_z[i] = 2.0 * tr_xxy_0[i] - 4.0 * tr_xxy_zz[i] * tke_0 - 4.0 * tr_xxyz_z[i] * tbe_0 - 2.0 * tr_xxyyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 237-240 components of targeted buffer : GP

    auto tr_y_0_z_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 237);

    auto tr_y_0_z_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 238);

    auto tr_y_0_z_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 239);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyz_0, tr_xxyyz_xz, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_y_0_z_xxyz_x, tr_y_0_z_xxyz_y, tr_y_0_z_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxyz_x[i] = tr_xx_x[i] - 2.0 * tr_xxz_xz[i] * tke_0 - 2.0 * tr_xxzz_x[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_y[i] = tr_xx_y[i] - 2.0 * tr_xxz_yz[i] * tke_0 - 2.0 * tr_xxzz_y[i] * tbe_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_z[i] = tr_xx_z[i] + tr_xxz_0[i] - 2.0 * tr_xxz_zz[i] * tke_0 - 2.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxyy_z[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 240-243 components of targeted buffer : GP

    auto tr_y_0_z_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 240);

    auto tr_y_0_z_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 241);

    auto tr_y_0_z_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 242);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzz_0, tr_xxyzz_xz, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_y_0_z_xxzz_x, tr_y_0_z_xxzz_y, tr_y_0_z_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxzz_x[i] = -4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_y[i] = -4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_z[i] = -4.0 * tr_xxyz_z[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 243-246 components of targeted buffer : GP

    auto tr_y_0_z_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 243);

    auto tr_y_0_z_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 244);

    auto tr_y_0_z_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 245);

    #pragma omp simd aligned(tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyyy_0, tr_xyyyy_xz, tr_xyyyy_yz, tr_xyyyy_zz, tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_y_0_z_xyyy_x, tr_y_0_z_xyyy_y, tr_y_0_z_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyyy_x[i] = -6.0 * tr_xyy_xz[i] * tke_0 - 6.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_y[i] = -6.0 * tr_xyy_yz[i] * tke_0 - 6.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_z[i] = 3.0 * tr_xyy_0[i] - 6.0 * tr_xyy_zz[i] * tke_0 - 6.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xyyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 246-249 components of targeted buffer : GP

    auto tr_y_0_z_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 246);

    auto tr_y_0_z_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 247);

    auto tr_y_0_z_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 248);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyz_0, tr_xyyyz_xz, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0_z_xyyz_x, tr_y_0_z_xyyz_y, tr_y_0_z_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyyz_x[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xyz_xz[i] * tke_0 - 4.0 * tr_xyzz_x[i] * tbe_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_y[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xyz_yz[i] * tke_0 - 4.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_z[i] = 2.0 * tr_xy_z[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_zz[i] * tke_0 - 4.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xyyy_z[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 249-252 components of targeted buffer : GP

    auto tr_y_0_z_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 249);

    auto tr_y_0_z_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 250);

    auto tr_y_0_z_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 251);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzz_0, tr_xyyzz_xz, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_y_0_z_xyzz_x, tr_y_0_z_xyzz_y, tr_y_0_z_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyzz_x[i] = 2.0 * tr_xz_x[i] - 2.0 * tr_xzz_xz[i] * tke_0 - 2.0 * tr_xzzz_x[i] * tbe_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_y[i] = 2.0 * tr_xz_y[i] - 2.0 * tr_xzz_yz[i] * tke_0 - 2.0 * tr_xzzz_y[i] * tbe_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_z[i] = 2.0 * tr_xz_z[i] + tr_xzz_0[i] - 2.0 * tr_xzz_zz[i] * tke_0 - 2.0 * tr_xzzz_z[i] * tbe_0 - 4.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 252-255 components of targeted buffer : GP

    auto tr_y_0_z_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 252);

    auto tr_y_0_z_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 253);

    auto tr_y_0_z_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 254);

    #pragma omp simd aligned(tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzz_0, tr_xyzzz_xz, tr_xyzzz_yz, tr_xyzzz_zz, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_y_0_z_xzzz_x, tr_y_0_z_xzzz_y, tr_y_0_z_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xzzz_x[i] = -6.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_y[i] = -6.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_z[i] = -6.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 255-258 components of targeted buffer : GP

    auto tr_y_0_z_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 255);

    auto tr_y_0_z_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 256);

    auto tr_y_0_z_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 257);

    #pragma omp simd aligned(tr_y_0_z_yyyy_x, tr_y_0_z_yyyy_y, tr_y_0_z_yyyy_z, tr_yyy_0, tr_yyy_xz, tr_yyy_yz, tr_yyy_zz, tr_yyyyy_0, tr_yyyyy_xz, tr_yyyyy_yz, tr_yyyyy_zz, tr_yyyyyz_x, tr_yyyyyz_y, tr_yyyyyz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyyy_x[i] = -8.0 * tr_yyy_xz[i] * tke_0 - 8.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_y[i] = -8.0 * tr_yyy_yz[i] * tke_0 - 8.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_z[i] = 4.0 * tr_yyy_0[i] - 8.0 * tr_yyy_zz[i] * tke_0 - 8.0 * tr_yyyz_z[i] * tbe_0 - 2.0 * tr_yyyyy_0[i] * tbe_0 + 4.0 * tr_yyyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 258-261 components of targeted buffer : GP

    auto tr_y_0_z_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 258);

    auto tr_y_0_z_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 259);

    auto tr_y_0_z_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 260);

    #pragma omp simd aligned(tr_y_0_z_yyyz_x, tr_y_0_z_yyyz_y, tr_y_0_z_yyyz_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyz_0, tr_yyyyz_xz, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyyyzz_x, tr_yyyyzz_y, tr_yyyyzz_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyyz_x[i] = 3.0 * tr_yy_x[i] - 6.0 * tr_yyz_xz[i] * tke_0 - 6.0 * tr_yyzz_x[i] * tbe_0 - 2.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_y[i] = 3.0 * tr_yy_y[i] - 6.0 * tr_yyz_yz[i] * tke_0 - 6.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_z[i] = 3.0 * tr_yy_z[i] + 3.0 * tr_yyz_0[i] - 6.0 * tr_yyz_zz[i] * tke_0 - 6.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_yyyy_z[i] * tbe_0 - 2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 261-264 components of targeted buffer : GP

    auto tr_y_0_z_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 261);

    auto tr_y_0_z_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 262);

    auto tr_y_0_z_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 263);

    #pragma omp simd aligned(tr_y_0_z_yyzz_x, tr_y_0_z_yyzz_y, tr_y_0_z_yyzz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyyzz_0, tr_yyyzz_xz, tr_yyyzz_yz, tr_yyyzz_zz, tr_yyyzzz_x, tr_yyyzzz_y, tr_yyyzzz_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyzz_x[i] = 4.0 * tr_yz_x[i] - 4.0 * tr_yzz_xz[i] * tke_0 - 4.0 * tr_yzzz_x[i] * tbe_0 - 4.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_y[i] = 4.0 * tr_yz_y[i] - 4.0 * tr_yzz_yz[i] * tke_0 - 4.0 * tr_yzzz_y[i] * tbe_0 - 4.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_z[i] = 4.0 * tr_yz_z[i] + 2.0 * tr_yzz_0[i] - 4.0 * tr_yzz_zz[i] * tke_0 - 4.0 * tr_yzzz_z[i] * tbe_0 - 4.0 * tr_yyyz_z[i] * tbe_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 264-267 components of targeted buffer : GP

    auto tr_y_0_z_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 264);

    auto tr_y_0_z_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 265);

    auto tr_y_0_z_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 266);

    #pragma omp simd aligned(tr_y_0_z_yzzz_x, tr_y_0_z_yzzz_y, tr_y_0_z_yzzz_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yyzzz_0, tr_yyzzz_xz, tr_yyzzz_yz, tr_yyzzz_zz, tr_yyzzzz_x, tr_yyzzzz_y, tr_yyzzzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xz, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yzzz_x[i] = 3.0 * tr_zz_x[i] - 2.0 * tr_zzz_xz[i] * tke_0 - 2.0 * tr_zzzz_x[i] * tbe_0 - 6.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_y[i] = 3.0 * tr_zz_y[i] - 2.0 * tr_zzz_yz[i] * tke_0 - 2.0 * tr_zzzz_y[i] * tbe_0 - 6.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_z[i] = 3.0 * tr_zz_z[i] + tr_zzz_0[i] - 2.0 * tr_zzz_zz[i] * tke_0 - 2.0 * tr_zzzz_z[i] * tbe_0 - 6.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 267-270 components of targeted buffer : GP

    auto tr_y_0_z_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 267);

    auto tr_y_0_z_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 268);

    auto tr_y_0_z_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 269);

    #pragma omp simd aligned(tr_y_0_z_zzzz_x, tr_y_0_z_zzzz_y, tr_y_0_z_zzzz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_yzzzz_0, tr_yzzzz_xz, tr_yzzzz_yz, tr_yzzzz_zz, tr_yzzzzz_x, tr_yzzzzz_y, tr_yzzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zzzz_x[i] = -8.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_y[i] = -8.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_z[i] = -8.0 * tr_yzzz_z[i] * tbe_0 - 2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 270-273 components of targeted buffer : GP

    auto tr_z_0_x_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 270);

    auto tr_z_0_x_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 271);

    auto tr_z_0_x_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 272);

    #pragma omp simd aligned(tr_xxxxxz_x, tr_xxxxxz_y, tr_xxxxxz_z, tr_xxxxz_0, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_z_0_x_xxxx_x, tr_z_0_x_xxxx_y, tr_z_0_x_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxx_x[i] = -8.0 * tr_xxxz_x[i] * tbe_0 - 2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_y[i] = -8.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_z[i] = -8.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 273-276 components of targeted buffer : GP

    auto tr_z_0_x_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 273);

    auto tr_z_0_x_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 274);

    auto tr_z_0_x_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 275);

    #pragma omp simd aligned(tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxxyz_0, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_z_0_x_xxxy_x, tr_z_0_x_xxxy_y, tr_z_0_x_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxy_x[i] = -6.0 * tr_xxyz_x[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_y[i] = -6.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_z[i] = -6.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 276-279 components of targeted buffer : GP

    auto tr_z_0_x_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 276);

    auto tr_z_0_x_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 277);

    auto tr_z_0_x_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 278);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxzz_x, tr_xxxxzz_y, tr_xxxxzz_z, tr_xxxzz_0, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_z_0_x_xxxz_x, tr_z_0_x_xxxz_y, tr_z_0_x_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxz_x[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xxzz_x[i] * tbe_0 + tr_xxx_0[i] - 2.0 * tr_xxx_xx[i] * tke_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_y[i] = 3.0 * tr_xx_y[i] - 6.0 * tr_xxzz_y[i] * tbe_0 - 2.0 * tr_xxx_xy[i] * tke_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_z[i] = 3.0 * tr_xx_z[i] - 6.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxx_xz[i] * tke_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 279-282 components of targeted buffer : GP

    auto tr_z_0_x_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 279);

    auto tr_z_0_x_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 280);

    auto tr_z_0_x_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 281);

    #pragma omp simd aligned(tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxyyz_0, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_z_0_x_xxyy_x, tr_z_0_x_xxyy_y, tr_z_0_x_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxyy_x[i] = -4.0 * tr_xyyz_x[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_y[i] = -4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_z[i] = -4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 282-285 components of targeted buffer : GP

    auto tr_z_0_x_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 282);

    auto tr_z_0_x_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 283);

    auto tr_z_0_x_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 284);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxyzz_0, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_z_0_x_xxyz_x, tr_z_0_x_xxyz_y, tr_z_0_x_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxyz_x[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xyzz_x[i] * tbe_0 + tr_xxy_0[i] - 2.0 * tr_xxy_xx[i] * tke_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_y[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xyzz_y[i] * tbe_0 - 2.0 * tr_xxy_xy[i] * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_z[i] = 2.0 * tr_xy_z[i] - 4.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xxy_xz[i] * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 285-288 components of targeted buffer : GP

    auto tr_z_0_x_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 285);

    auto tr_z_0_x_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 286);

    auto tr_z_0_x_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 287);

    #pragma omp simd aligned(tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxxzzz_x, tr_xxxzzz_y, tr_xxxzzz_z, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxzzz_0, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_z_0_x_xxzz_x, tr_z_0_x_xxzz_y, tr_z_0_x_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxzz_x[i] = 4.0 * tr_xz_x[i] - 4.0 * tr_xzzz_x[i] * tbe_0 + 2.0 * tr_xxz_0[i] - 4.0 * tr_xxz_xx[i] * tke_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_y[i] = 4.0 * tr_xz_y[i] - 4.0 * tr_xzzz_y[i] * tbe_0 - 4.0 * tr_xxz_xy[i] * tke_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_z[i] = 4.0 * tr_xz_z[i] - 4.0 * tr_xzzz_z[i] * tbe_0 - 4.0 * tr_xxz_xz[i] * tke_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 288-291 components of targeted buffer : GP

    auto tr_z_0_x_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 288);

    auto tr_z_0_x_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 289);

    auto tr_z_0_x_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 290);

    #pragma omp simd aligned(tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xyyyz_0, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_z_0_x_xyyy_x, tr_z_0_x_xyyy_y, tr_z_0_x_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyyy_x[i] = -2.0 * tr_yyyz_x[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_y[i] = -2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_z[i] = -2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 291-294 components of targeted buffer : GP

    auto tr_z_0_x_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 291);

    auto tr_z_0_x_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 292);

    auto tr_z_0_x_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 293);

    #pragma omp simd aligned(tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyyzz_0, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_z_0_x_xyyz_x, tr_z_0_x_xyyz_y, tr_z_0_x_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyyz_x[i] = tr_yy_x[i] - 2.0 * tr_yyzz_x[i] * tbe_0 + tr_xyy_0[i] - 2.0 * tr_xyy_xx[i] * tke_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_y[i] = tr_yy_y[i] - 2.0 * tr_yyzz_y[i] * tbe_0 - 2.0 * tr_xyy_xy[i] * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_z[i] = tr_yy_z[i] - 2.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_xyy_xz[i] * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 294-297 components of targeted buffer : GP

    auto tr_z_0_x_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 294);

    auto tr_z_0_x_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 295);

    auto tr_z_0_x_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 296);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyzzz_0, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_z_0_x_xyzz_x, tr_z_0_x_xyzz_y, tr_z_0_x_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyzz_x[i] = 2.0 * tr_yz_x[i] - 2.0 * tr_yzzz_x[i] * tbe_0 + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_xx[i] * tke_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_y[i] = 2.0 * tr_yz_y[i] - 2.0 * tr_yzzz_y[i] * tbe_0 - 4.0 * tr_xyz_xy[i] * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_z[i] = 2.0 * tr_yz_z[i] - 2.0 * tr_yzzz_z[i] * tbe_0 - 4.0 * tr_xyz_xz[i] * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 297-300 components of targeted buffer : GP

    auto tr_z_0_x_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 297);

    auto tr_z_0_x_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 298);

    auto tr_z_0_x_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 299);

    #pragma omp simd aligned(tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xxzzzz_x, tr_xxzzzz_y, tr_xxzzzz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzzzz_0, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_z_0_x_xzzz_x, tr_z_0_x_xzzz_y, tr_z_0_x_xzzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xzzz_x[i] = 3.0 * tr_zz_x[i] - 2.0 * tr_zzzz_x[i] * tbe_0 + 3.0 * tr_xzz_0[i] - 6.0 * tr_xzz_xx[i] * tke_0 - 2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_y[i] = 3.0 * tr_zz_y[i] - 2.0 * tr_zzzz_y[i] * tbe_0 - 6.0 * tr_xzz_xy[i] * tke_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_z[i] = 3.0 * tr_zz_z[i] - 2.0 * tr_zzzz_z[i] * tbe_0 - 6.0 * tr_xzz_xz[i] * tke_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 300-303 components of targeted buffer : GP

    auto tr_z_0_x_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 300);

    auto tr_z_0_x_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 301);

    auto tr_z_0_x_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 302);

    #pragma omp simd aligned(tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_yyyyz_0, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_z_0_x_yyyy_x, tr_z_0_x_yyyy_y, tr_z_0_x_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyyy_x[i] = -2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_y[i] = 4.0 * tr_yyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_z[i] = 4.0 * tr_yyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 303-306 components of targeted buffer : GP

    auto tr_z_0_x_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 303);

    auto tr_z_0_x_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 304);

    auto tr_z_0_x_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 305);

    #pragma omp simd aligned(tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyyzz_0, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_z_0_x_yyyz_x, tr_z_0_x_yyyz_y, tr_z_0_x_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyyz_x[i] = tr_yyy_0[i] - 2.0 * tr_yyy_xx[i] * tke_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_y[i] = -2.0 * tr_yyy_xy[i] * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_z[i] = -2.0 * tr_yyy_xz[i] * tke_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 306-309 components of targeted buffer : GP

    auto tr_z_0_x_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 306);

    auto tr_z_0_x_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 307);

    auto tr_z_0_x_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 308);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyzzz_0, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_z_0_x_yyzz_x, tr_z_0_x_yyzz_y, tr_z_0_x_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyzz_x[i] = 2.0 * tr_yyz_0[i] - 4.0 * tr_yyz_xx[i] * tke_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_y[i] = -4.0 * tr_yyz_xy[i] * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_z[i] = -4.0 * tr_yyz_xz[i] * tke_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 309-312 components of targeted buffer : GP

    auto tr_z_0_x_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 309);

    auto tr_z_0_x_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 310);

    auto tr_z_0_x_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 311);

    #pragma omp simd aligned(tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzzzz_0, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_z_0_x_yzzz_x, tr_z_0_x_yzzz_y, tr_z_0_x_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yzzz_x[i] = 3.0 * tr_yzz_0[i] - 6.0 * tr_yzz_xx[i] * tke_0 - 2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_y[i] = -6.0 * tr_yzz_xy[i] * tke_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_z[i] = -6.0 * tr_yzz_xz[i] * tke_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 312-315 components of targeted buffer : GP

    auto tr_z_0_x_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 312);

    auto tr_z_0_x_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 313);

    auto tr_z_0_x_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 314);

    #pragma omp simd aligned(tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_xzzzzz_x, tr_xzzzzz_y, tr_xzzzzz_z, tr_z_0_x_zzzz_x, tr_z_0_x_zzzz_y, tr_z_0_x_zzzz_z, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzzzz_0, tr_zzzzz_xx, tr_zzzzz_xy, tr_zzzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zzzz_x[i] = 4.0 * tr_zzz_0[i] - 8.0 * tr_zzz_xx[i] * tke_0 - 2.0 * tr_zzzzz_0[i] * tbe_0 + 4.0 * tr_zzzzz_xx[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_y[i] = -8.0 * tr_zzz_xy[i] * tke_0 + 4.0 * tr_zzzzz_xy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_z[i] = -8.0 * tr_zzz_xz[i] * tke_0 + 4.0 * tr_zzzzz_xz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 315-318 components of targeted buffer : GP

    auto tr_z_0_y_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 315);

    auto tr_z_0_y_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 316);

    auto tr_z_0_y_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 317);

    #pragma omp simd aligned(tr_xxxxyz_x, tr_xxxxyz_y, tr_xxxxyz_z, tr_xxxxz_0, tr_xxxxz_xy, tr_xxxxz_yy, tr_xxxxz_yz, tr_z_0_y_xxxx_x, tr_z_0_y_xxxx_y, tr_z_0_y_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxx_x[i] = 4.0 * tr_xxxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_y[i] = -2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_z[i] = 4.0 * tr_xxxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 318-321 components of targeted buffer : GP

    auto tr_z_0_y_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 318);

    auto tr_z_0_y_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 319);

    auto tr_z_0_y_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 320);

    #pragma omp simd aligned(tr_xxxyyz_x, tr_xxxyyz_y, tr_xxxyyz_z, tr_xxxyz_0, tr_xxxyz_xy, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_z_0_y_xxxy_x, tr_z_0_y_xxxy_y, tr_z_0_y_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxy_x[i] = -2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_y[i] = -2.0 * tr_xxxz_y[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_z[i] = -2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 321-324 components of targeted buffer : GP

    auto tr_z_0_y_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 321);

    auto tr_z_0_y_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 322);

    auto tr_z_0_y_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 323);

    #pragma omp simd aligned(tr_xxx_0, tr_xxx_xy, tr_xxx_yy, tr_xxx_yz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_xxxzz_0, tr_xxxzz_xy, tr_xxxzz_yy, tr_xxxzz_yz, tr_z_0_y_xxxz_x, tr_z_0_y_xxxz_y, tr_z_0_y_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxz_x[i] = -2.0 * tr_xxx_xy[i] * tke_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_y[i] = tr_xxx_0[i] - 2.0 * tr_xxx_yy[i] * tke_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_z[i] = -2.0 * tr_xxx_yz[i] * tke_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 324-327 components of targeted buffer : GP

    auto tr_z_0_y_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 324);

    auto tr_z_0_y_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 325);

    auto tr_z_0_y_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 326);

    #pragma omp simd aligned(tr_xxyyyz_x, tr_xxyyyz_y, tr_xxyyyz_z, tr_xxyyz_0, tr_xxyyz_xy, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_z_0_y_xxyy_x, tr_z_0_y_xxyy_y, tr_z_0_y_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxyy_x[i] = -4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_y[i] = -4.0 * tr_xxyz_y[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_z[i] = -4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 327-330 components of targeted buffer : GP

    auto tr_z_0_y_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 327);

    auto tr_z_0_y_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 328);

    auto tr_z_0_y_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 329);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_xxyzz_0, tr_xxyzz_xy, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_z_0_y_xxyz_x, tr_z_0_y_xxyz_y, tr_z_0_y_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxyz_x[i] = tr_xx_x[i] - 2.0 * tr_xxzz_x[i] * tbe_0 - 2.0 * tr_xxy_xy[i] * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_y[i] = tr_xx_y[i] - 2.0 * tr_xxzz_y[i] * tbe_0 + tr_xxy_0[i] - 2.0 * tr_xxy_yy[i] * tke_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_z[i] = tr_xx_z[i] - 2.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxy_yz[i] * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 330-333 components of targeted buffer : GP

    auto tr_z_0_y_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 330);

    auto tr_z_0_y_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 331);

    auto tr_z_0_y_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 332);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, tr_xxzzz_0, tr_xxzzz_xy, tr_xxzzz_yy, tr_xxzzz_yz, tr_z_0_y_xxzz_x, tr_z_0_y_xxzz_y, tr_z_0_y_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxzz_x[i] = -4.0 * tr_xxz_xy[i] * tke_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_y[i] = 2.0 * tr_xxz_0[i] - 4.0 * tr_xxz_yy[i] * tke_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_z[i] = -4.0 * tr_xxz_yz[i] * tke_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 333-336 components of targeted buffer : GP

    auto tr_z_0_y_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 333);

    auto tr_z_0_y_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 334);

    auto tr_z_0_y_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 335);

    #pragma omp simd aligned(tr_xyyyyz_x, tr_xyyyyz_y, tr_xyyyyz_z, tr_xyyyz_0, tr_xyyyz_xy, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_z_0_y_xyyy_x, tr_z_0_y_xyyy_y, tr_z_0_y_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyyy_x[i] = -6.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_y[i] = -6.0 * tr_xyyz_y[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_z[i] = -6.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 336-339 components of targeted buffer : GP

    auto tr_z_0_y_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 336);

    auto tr_z_0_y_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 337);

    auto tr_z_0_y_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 338);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_xyyzz_0, tr_xyyzz_xy, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_z_0_y_xyyz_x, tr_z_0_y_xyyz_y, tr_z_0_y_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyyz_x[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xyzz_x[i] * tbe_0 - 2.0 * tr_xyy_xy[i] * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_y[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xyzz_y[i] * tbe_0 + tr_xyy_0[i] - 2.0 * tr_xyy_yy[i] * tke_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_z[i] = 2.0 * tr_xy_z[i] - 4.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xyy_yz[i] * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 339-342 components of targeted buffer : GP

    auto tr_z_0_y_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 339);

    auto tr_z_0_y_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 340);

    auto tr_z_0_y_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 341);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xyzzz_0, tr_xyzzz_xy, tr_xyzzz_yy, tr_xyzzz_yz, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_z_0_y_xyzz_x, tr_z_0_y_xyzz_y, tr_z_0_y_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyzz_x[i] = 2.0 * tr_xz_x[i] - 2.0 * tr_xzzz_x[i] * tbe_0 - 4.0 * tr_xyz_xy[i] * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_y[i] = 2.0 * tr_xz_y[i] - 2.0 * tr_xzzz_y[i] * tbe_0 + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_yy[i] * tke_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_z[i] = 2.0 * tr_xz_z[i] - 2.0 * tr_xzzz_z[i] * tbe_0 - 4.0 * tr_xyz_yz[i] * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 342-345 components of targeted buffer : GP

    auto tr_z_0_y_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 342);

    auto tr_z_0_y_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 343);

    auto tr_z_0_y_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 344);

    #pragma omp simd aligned(tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, tr_xzzzz_0, tr_xzzzz_xy, tr_xzzzz_yy, tr_xzzzz_yz, tr_z_0_y_xzzz_x, tr_z_0_y_xzzz_y, tr_z_0_y_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xzzz_x[i] = -6.0 * tr_xzz_xy[i] * tke_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_y[i] = 3.0 * tr_xzz_0[i] - 6.0 * tr_xzz_yy[i] * tke_0 - 2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_yy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_z[i] = -6.0 * tr_xzz_yz[i] * tke_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 345-348 components of targeted buffer : GP

    auto tr_z_0_y_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 345);

    auto tr_z_0_y_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 346);

    auto tr_z_0_y_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 347);

    #pragma omp simd aligned(tr_yyyyyz_x, tr_yyyyyz_y, tr_yyyyyz_z, tr_yyyyz_0, tr_yyyyz_xy, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_z_0_y_yyyy_x, tr_z_0_y_yyyy_y, tr_z_0_y_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyyy_x[i] = -8.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_y[i] = -8.0 * tr_yyyz_y[i] * tbe_0 - 2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_z[i] = -8.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 348-351 components of targeted buffer : GP

    auto tr_z_0_y_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 348);

    auto tr_z_0_y_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 349);

    auto tr_z_0_y_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 350);

    #pragma omp simd aligned(tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xy, tr_yyy_yy, tr_yyy_yz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyzz_x, tr_yyyyzz_y, tr_yyyyzz_z, tr_yyyzz_0, tr_yyyzz_xy, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_z_0_y_yyyz_x, tr_z_0_y_yyyz_y, tr_z_0_y_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyyz_x[i] = 3.0 * tr_yy_x[i] - 6.0 * tr_yyzz_x[i] * tbe_0 - 2.0 * tr_yyy_xy[i] * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_y[i] = 3.0 * tr_yy_y[i] - 6.0 * tr_yyzz_y[i] * tbe_0 + tr_yyy_0[i] - 2.0 * tr_yyy_yy[i] * tke_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_z[i] = 3.0 * tr_yy_z[i] - 6.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_yyy_yz[i] * tke_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_yyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 351-354 components of targeted buffer : GP

    auto tr_z_0_y_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 351);

    auto tr_z_0_y_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 352);

    auto tr_z_0_y_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 353);

    #pragma omp simd aligned(tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyyzzz_x, tr_yyyzzz_y, tr_yyyzzz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yyzzz_0, tr_yyzzz_xy, tr_yyzzz_yy, tr_yyzzz_yz, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_z_0_y_yyzz_x, tr_z_0_y_yyzz_y, tr_z_0_y_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyzz_x[i] = 4.0 * tr_yz_x[i] - 4.0 * tr_yzzz_x[i] * tbe_0 - 4.0 * tr_yyz_xy[i] * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_y[i] = 4.0 * tr_yz_y[i] - 4.0 * tr_yzzz_y[i] * tbe_0 + 2.0 * tr_yyz_0[i] - 4.0 * tr_yyz_yy[i] * tke_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_z[i] = 4.0 * tr_yz_z[i] - 4.0 * tr_yzzz_z[i] * tbe_0 - 4.0 * tr_yyz_yz[i] * tke_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 354-357 components of targeted buffer : GP

    auto tr_z_0_y_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 354);

    auto tr_z_0_y_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 355);

    auto tr_z_0_y_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 356);

    #pragma omp simd aligned(tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yyzzzz_x, tr_yyzzzz_y, tr_yyzzzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_yzzzz_0, tr_yzzzz_xy, tr_yzzzz_yy, tr_yzzzz_yz, tr_z_0_y_yzzz_x, tr_z_0_y_yzzz_y, tr_z_0_y_yzzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yzzz_x[i] = 3.0 * tr_zz_x[i] - 2.0 * tr_zzzz_x[i] * tbe_0 - 6.0 * tr_yzz_xy[i] * tke_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_y[i] = 3.0 * tr_zz_y[i] - 2.0 * tr_zzzz_y[i] * tbe_0 + 3.0 * tr_yzz_0[i] - 6.0 * tr_yzz_yy[i] * tke_0 - 2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_yy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_z[i] = 3.0 * tr_zz_z[i] - 2.0 * tr_zzzz_z[i] * tbe_0 - 6.0 * tr_yzz_yz[i] * tke_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 357-360 components of targeted buffer : GP

    auto tr_z_0_y_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 357);

    auto tr_z_0_y_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 358);

    auto tr_z_0_y_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 359);

    #pragma omp simd aligned(tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_yzzzzz_x, tr_yzzzzz_y, tr_yzzzzz_z, tr_z_0_y_zzzz_x, tr_z_0_y_zzzz_y, tr_z_0_y_zzzz_z, tr_zzz_0, tr_zzz_xy, tr_zzz_yy, tr_zzz_yz, tr_zzzzz_0, tr_zzzzz_xy, tr_zzzzz_yy, tr_zzzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zzzz_x[i] = -8.0 * tr_zzz_xy[i] * tke_0 + 4.0 * tr_zzzzz_xy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_y[i] = 4.0 * tr_zzz_0[i] - 8.0 * tr_zzz_yy[i] * tke_0 - 2.0 * tr_zzzzz_0[i] * tbe_0 + 4.0 * tr_zzzzz_yy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_z[i] = -8.0 * tr_zzz_yz[i] * tke_0 + 4.0 * tr_zzzzz_yz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 360-363 components of targeted buffer : GP

    auto tr_z_0_z_xxxx_x = pbuffer.data(idx_op_geom_110_gp + 360);

    auto tr_z_0_z_xxxx_y = pbuffer.data(idx_op_geom_110_gp + 361);

    auto tr_z_0_z_xxxx_z = pbuffer.data(idx_op_geom_110_gp + 362);

    #pragma omp simd aligned(tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxz_0, tr_xxxxz_xz, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxxzz_x, tr_xxxxzz_y, tr_xxxxzz_z, tr_z_0_z_xxxx_x, tr_z_0_z_xxxx_y, tr_z_0_z_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxx_x[i] = -2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_y[i] = -2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_z[i] = -2.0 * tr_xxxx_z[i] * tbe_0 - 2.0 * tr_xxxxz_0[i] * tbe_0 + 4.0 * tr_xxxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 363-366 components of targeted buffer : GP

    auto tr_z_0_z_xxxy_x = pbuffer.data(idx_op_geom_110_gp + 363);

    auto tr_z_0_z_xxxy_y = pbuffer.data(idx_op_geom_110_gp + 364);

    auto tr_z_0_z_xxxy_z = pbuffer.data(idx_op_geom_110_gp + 365);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyz_0, tr_xxxyz_xz, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxyzz_x, tr_xxxyzz_y, tr_xxxyzz_z, tr_z_0_z_xxxy_x, tr_z_0_z_xxxy_y, tr_z_0_z_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxy_x[i] = -2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_y[i] = -2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_z[i] = -2.0 * tr_xxxy_z[i] * tbe_0 - 2.0 * tr_xxxyz_0[i] * tbe_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 366-369 components of targeted buffer : GP

    auto tr_z_0_z_xxxz_x = pbuffer.data(idx_op_geom_110_gp + 366);

    auto tr_z_0_z_xxxz_y = pbuffer.data(idx_op_geom_110_gp + 367);

    auto tr_z_0_z_xxxz_z = pbuffer.data(idx_op_geom_110_gp + 368);

    #pragma omp simd aligned(tr_xxx_0, tr_xxx_xz, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxxzz_0, tr_xxxzz_xz, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxxzzz_x, tr_xxxzzz_y, tr_xxxzzz_z, tr_z_0_z_xxxz_x, tr_z_0_z_xxxz_y, tr_z_0_z_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxz_x[i] = -2.0 * tr_xxx_xz[i] * tke_0 - 6.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_y[i] = -2.0 * tr_xxx_yz[i] * tke_0 - 6.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_z[i] = tr_xxx_0[i] - 2.0 * tr_xxx_zz[i] * tke_0 - 6.0 * tr_xxxz_z[i] * tbe_0 - 2.0 * tr_xxxzz_0[i] * tbe_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 369-372 components of targeted buffer : GP

    auto tr_z_0_z_xxyy_x = pbuffer.data(idx_op_geom_110_gp + 369);

    auto tr_z_0_z_xxyy_y = pbuffer.data(idx_op_geom_110_gp + 370);

    auto tr_z_0_z_xxyy_z = pbuffer.data(idx_op_geom_110_gp + 371);

    #pragma omp simd aligned(tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyz_0, tr_xxyyz_xz, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyyzz_x, tr_xxyyzz_y, tr_xxyyzz_z, tr_z_0_z_xxyy_x, tr_z_0_z_xxyy_y, tr_z_0_z_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxyy_x[i] = -2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_y[i] = -2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_z[i] = -2.0 * tr_xxyy_z[i] * tbe_0 - 2.0 * tr_xxyyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 372-375 components of targeted buffer : GP

    auto tr_z_0_z_xxyz_x = pbuffer.data(idx_op_geom_110_gp + 372);

    auto tr_z_0_z_xxyz_y = pbuffer.data(idx_op_geom_110_gp + 373);

    auto tr_z_0_z_xxyz_z = pbuffer.data(idx_op_geom_110_gp + 374);

    #pragma omp simd aligned(tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzz_0, tr_xxyzz_xz, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxyzzz_x, tr_xxyzzz_y, tr_xxyzzz_z, tr_z_0_z_xxyz_x, tr_z_0_z_xxyz_y, tr_z_0_z_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxyz_x[i] = -2.0 * tr_xxy_xz[i] * tke_0 - 6.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_y[i] = -2.0 * tr_xxy_yz[i] * tke_0 - 6.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_z[i] = tr_xxy_0[i] - 2.0 * tr_xxy_zz[i] * tke_0 - 6.0 * tr_xxyz_z[i] * tbe_0 - 2.0 * tr_xxyzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 375-378 components of targeted buffer : GP

    auto tr_z_0_z_xxzz_x = pbuffer.data(idx_op_geom_110_gp + 375);

    auto tr_z_0_z_xxzz_y = pbuffer.data(idx_op_geom_110_gp + 376);

    auto tr_z_0_z_xxzz_z = pbuffer.data(idx_op_geom_110_gp + 377);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xxzzz_0, tr_xxzzz_xz, tr_xxzzz_yz, tr_xxzzz_zz, tr_xxzzzz_x, tr_xxzzzz_y, tr_xxzzzz_z, tr_z_0_z_xxzz_x, tr_z_0_z_xxzz_y, tr_z_0_z_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxzz_x[i] = 2.0 * tr_xx_x[i] - 4.0 * tr_xxz_xz[i] * tke_0 - 10.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_y[i] = 2.0 * tr_xx_y[i] - 4.0 * tr_xxz_yz[i] * tke_0 - 10.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_z[i] = 2.0 * tr_xx_z[i] + 2.0 * tr_xxz_0[i] - 4.0 * tr_xxz_zz[i] * tke_0 - 10.0 * tr_xxzz_z[i] * tbe_0 - 2.0 * tr_xxzzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 378-381 components of targeted buffer : GP

    auto tr_z_0_z_xyyy_x = pbuffer.data(idx_op_geom_110_gp + 378);

    auto tr_z_0_z_xyyy_y = pbuffer.data(idx_op_geom_110_gp + 379);

    auto tr_z_0_z_xyyy_z = pbuffer.data(idx_op_geom_110_gp + 380);

    #pragma omp simd aligned(tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyz_0, tr_xyyyz_xz, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyyzz_x, tr_xyyyzz_y, tr_xyyyzz_z, tr_z_0_z_xyyy_x, tr_z_0_z_xyyy_y, tr_z_0_z_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyyy_x[i] = -2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_y[i] = -2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_z[i] = -2.0 * tr_xyyy_z[i] * tbe_0 - 2.0 * tr_xyyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 381-384 components of targeted buffer : GP

    auto tr_z_0_z_xyyz_x = pbuffer.data(idx_op_geom_110_gp + 381);

    auto tr_z_0_z_xyyz_y = pbuffer.data(idx_op_geom_110_gp + 382);

    auto tr_z_0_z_xyyz_z = pbuffer.data(idx_op_geom_110_gp + 383);

    #pragma omp simd aligned(tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzz_0, tr_xyyzz_xz, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyyzzz_x, tr_xyyzzz_y, tr_xyyzzz_z, tr_z_0_z_xyyz_x, tr_z_0_z_xyyz_y, tr_z_0_z_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyyz_x[i] = -2.0 * tr_xyy_xz[i] * tke_0 - 6.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_y[i] = -2.0 * tr_xyy_yz[i] * tke_0 - 6.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_z[i] = tr_xyy_0[i] - 2.0 * tr_xyy_zz[i] * tke_0 - 6.0 * tr_xyyz_z[i] * tbe_0 - 2.0 * tr_xyyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 384-387 components of targeted buffer : GP

    auto tr_z_0_z_xyzz_x = pbuffer.data(idx_op_geom_110_gp + 384);

    auto tr_z_0_z_xyzz_y = pbuffer.data(idx_op_geom_110_gp + 385);

    auto tr_z_0_z_xyzz_z = pbuffer.data(idx_op_geom_110_gp + 386);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzz_0, tr_xyzzz_xz, tr_xyzzz_yz, tr_xyzzz_zz, tr_xyzzzz_x, tr_xyzzzz_y, tr_xyzzzz_z, tr_z_0_z_xyzz_x, tr_z_0_z_xyzz_y, tr_z_0_z_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyzz_x[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xyz_xz[i] * tke_0 - 10.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_y[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xyz_yz[i] * tke_0 - 10.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_z[i] = 2.0 * tr_xy_z[i] + 2.0 * tr_xyz_0[i] - 4.0 * tr_xyz_zz[i] * tke_0 - 10.0 * tr_xyzz_z[i] * tbe_0 - 2.0 * tr_xyzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 387-390 components of targeted buffer : GP

    auto tr_z_0_z_xzzz_x = pbuffer.data(idx_op_geom_110_gp + 387);

    auto tr_z_0_z_xzzz_y = pbuffer.data(idx_op_geom_110_gp + 388);

    auto tr_z_0_z_xzzz_z = pbuffer.data(idx_op_geom_110_gp + 389);

    #pragma omp simd aligned(tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_xzzzz_0, tr_xzzzz_xz, tr_xzzzz_yz, tr_xzzzz_zz, tr_xzzzzz_x, tr_xzzzzz_y, tr_xzzzzz_z, tr_z_0_z_xzzz_x, tr_z_0_z_xzzz_y, tr_z_0_z_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xzzz_x[i] = 6.0 * tr_xz_x[i] - 6.0 * tr_xzz_xz[i] * tke_0 - 14.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_y[i] = 6.0 * tr_xz_y[i] - 6.0 * tr_xzz_yz[i] * tke_0 - 14.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_z[i] = 6.0 * tr_xz_z[i] + 3.0 * tr_xzz_0[i] - 6.0 * tr_xzz_zz[i] * tke_0 - 14.0 * tr_xzzz_z[i] * tbe_0 - 2.0 * tr_xzzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 390-393 components of targeted buffer : GP

    auto tr_z_0_z_yyyy_x = pbuffer.data(idx_op_geom_110_gp + 390);

    auto tr_z_0_z_yyyy_y = pbuffer.data(idx_op_geom_110_gp + 391);

    auto tr_z_0_z_yyyy_z = pbuffer.data(idx_op_geom_110_gp + 392);

    #pragma omp simd aligned(tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyz_0, tr_yyyyz_xz, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyyyzz_x, tr_yyyyzz_y, tr_yyyyzz_z, tr_z_0_z_yyyy_x, tr_z_0_z_yyyy_y, tr_z_0_z_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyyy_x[i] = -2.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_y[i] = -2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_z[i] = -2.0 * tr_yyyy_z[i] * tbe_0 - 2.0 * tr_yyyyz_0[i] * tbe_0 + 4.0 * tr_yyyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 393-396 components of targeted buffer : GP

    auto tr_z_0_z_yyyz_x = pbuffer.data(idx_op_geom_110_gp + 393);

    auto tr_z_0_z_yyyz_y = pbuffer.data(idx_op_geom_110_gp + 394);

    auto tr_z_0_z_yyyz_z = pbuffer.data(idx_op_geom_110_gp + 395);

    #pragma omp simd aligned(tr_yyy_0, tr_yyy_xz, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyyzz_0, tr_yyyzz_xz, tr_yyyzz_yz, tr_yyyzz_zz, tr_yyyzzz_x, tr_yyyzzz_y, tr_yyyzzz_z, tr_z_0_z_yyyz_x, tr_z_0_z_yyyz_y, tr_z_0_z_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyyz_x[i] = -2.0 * tr_yyy_xz[i] * tke_0 - 6.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_y[i] = -2.0 * tr_yyy_yz[i] * tke_0 - 6.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_z[i] = tr_yyy_0[i] - 2.0 * tr_yyy_zz[i] * tke_0 - 6.0 * tr_yyyz_z[i] * tbe_0 - 2.0 * tr_yyyzz_0[i] * tbe_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 396-399 components of targeted buffer : GP

    auto tr_z_0_z_yyzz_x = pbuffer.data(idx_op_geom_110_gp + 396);

    auto tr_z_0_z_yyzz_y = pbuffer.data(idx_op_geom_110_gp + 397);

    auto tr_z_0_z_yyzz_z = pbuffer.data(idx_op_geom_110_gp + 398);

    #pragma omp simd aligned(tr_yy_x, tr_yy_y, tr_yy_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yyzzz_0, tr_yyzzz_xz, tr_yyzzz_yz, tr_yyzzz_zz, tr_yyzzzz_x, tr_yyzzzz_y, tr_yyzzzz_z, tr_z_0_z_yyzz_x, tr_z_0_z_yyzz_y, tr_z_0_z_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyzz_x[i] = 2.0 * tr_yy_x[i] - 4.0 * tr_yyz_xz[i] * tke_0 - 10.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_y[i] = 2.0 * tr_yy_y[i] - 4.0 * tr_yyz_yz[i] * tke_0 - 10.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_z[i] = 2.0 * tr_yy_z[i] + 2.0 * tr_yyz_0[i] - 4.0 * tr_yyz_zz[i] * tke_0 - 10.0 * tr_yyzz_z[i] * tbe_0 - 2.0 * tr_yyzzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 399-402 components of targeted buffer : GP

    auto tr_z_0_z_yzzz_x = pbuffer.data(idx_op_geom_110_gp + 399);

    auto tr_z_0_z_yzzz_y = pbuffer.data(idx_op_geom_110_gp + 400);

    auto tr_z_0_z_yzzz_z = pbuffer.data(idx_op_geom_110_gp + 401);

    #pragma omp simd aligned(tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_yzzzz_0, tr_yzzzz_xz, tr_yzzzz_yz, tr_yzzzz_zz, tr_yzzzzz_x, tr_yzzzzz_y, tr_yzzzzz_z, tr_z_0_z_yzzz_x, tr_z_0_z_yzzz_y, tr_z_0_z_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yzzz_x[i] = 6.0 * tr_yz_x[i] - 6.0 * tr_yzz_xz[i] * tke_0 - 14.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_y[i] = 6.0 * tr_yz_y[i] - 6.0 * tr_yzz_yz[i] * tke_0 - 14.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_z[i] = 6.0 * tr_yz_z[i] + 3.0 * tr_yzz_0[i] - 6.0 * tr_yzz_zz[i] * tke_0 - 14.0 * tr_yzzz_z[i] * tbe_0 - 2.0 * tr_yzzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 402-405 components of targeted buffer : GP

    auto tr_z_0_z_zzzz_x = pbuffer.data(idx_op_geom_110_gp + 402);

    auto tr_z_0_z_zzzz_y = pbuffer.data(idx_op_geom_110_gp + 403);

    auto tr_z_0_z_zzzz_z = pbuffer.data(idx_op_geom_110_gp + 404);

    #pragma omp simd aligned(tr_z_0_z_zzzz_x, tr_z_0_z_zzzz_y, tr_z_0_z_zzzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xz, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, tr_zzzzz_0, tr_zzzzz_xz, tr_zzzzz_yz, tr_zzzzz_zz, tr_zzzzzz_x, tr_zzzzzz_y, tr_zzzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zzzz_x[i] = 12.0 * tr_zz_x[i] - 8.0 * tr_zzz_xz[i] * tke_0 - 18.0 * tr_zzzz_x[i] * tbe_0 + 4.0 * tr_zzzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_y[i] = 12.0 * tr_zz_y[i] - 8.0 * tr_zzz_yz[i] * tke_0 - 18.0 * tr_zzzz_y[i] * tbe_0 + 4.0 * tr_zzzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_z[i] = 12.0 * tr_zz_z[i] + 4.0 * tr_zzz_0[i] - 8.0 * tr_zzz_zz[i] * tke_0 - 18.0 * tr_zzzz_z[i] * tbe_0 - 2.0 * tr_zzzzz_0[i] * tbe_0 + 4.0 * tr_zzzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_z[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

