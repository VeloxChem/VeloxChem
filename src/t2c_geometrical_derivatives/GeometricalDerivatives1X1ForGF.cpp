#include "GeometricalDerivatives1X1ForGF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_gf(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_gf,
                        const size_t idx_op_fd,
                        const size_t idx_op_fg,
                        const size_t idx_op_hd,
                        const size_t idx_op_hg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FD

        auto to_xxx_xx = pbuffer.data(idx_op_fd + i * 60 + 0);

        auto to_xxx_xy = pbuffer.data(idx_op_fd + i * 60 + 1);

        auto to_xxx_xz = pbuffer.data(idx_op_fd + i * 60 + 2);

        auto to_xxx_yy = pbuffer.data(idx_op_fd + i * 60 + 3);

        auto to_xxx_yz = pbuffer.data(idx_op_fd + i * 60 + 4);

        auto to_xxx_zz = pbuffer.data(idx_op_fd + i * 60 + 5);

        auto to_xxy_xx = pbuffer.data(idx_op_fd + i * 60 + 6);

        auto to_xxy_xy = pbuffer.data(idx_op_fd + i * 60 + 7);

        auto to_xxy_xz = pbuffer.data(idx_op_fd + i * 60 + 8);

        auto to_xxy_yy = pbuffer.data(idx_op_fd + i * 60 + 9);

        auto to_xxy_yz = pbuffer.data(idx_op_fd + i * 60 + 10);

        auto to_xxy_zz = pbuffer.data(idx_op_fd + i * 60 + 11);

        auto to_xxz_xx = pbuffer.data(idx_op_fd + i * 60 + 12);

        auto to_xxz_xy = pbuffer.data(idx_op_fd + i * 60 + 13);

        auto to_xxz_xz = pbuffer.data(idx_op_fd + i * 60 + 14);

        auto to_xxz_yy = pbuffer.data(idx_op_fd + i * 60 + 15);

        auto to_xxz_yz = pbuffer.data(idx_op_fd + i * 60 + 16);

        auto to_xxz_zz = pbuffer.data(idx_op_fd + i * 60 + 17);

        auto to_xyy_xx = pbuffer.data(idx_op_fd + i * 60 + 18);

        auto to_xyy_xy = pbuffer.data(idx_op_fd + i * 60 + 19);

        auto to_xyy_xz = pbuffer.data(idx_op_fd + i * 60 + 20);

        auto to_xyy_yy = pbuffer.data(idx_op_fd + i * 60 + 21);

        auto to_xyy_yz = pbuffer.data(idx_op_fd + i * 60 + 22);

        auto to_xyy_zz = pbuffer.data(idx_op_fd + i * 60 + 23);

        auto to_xyz_xx = pbuffer.data(idx_op_fd + i * 60 + 24);

        auto to_xyz_xy = pbuffer.data(idx_op_fd + i * 60 + 25);

        auto to_xyz_xz = pbuffer.data(idx_op_fd + i * 60 + 26);

        auto to_xyz_yy = pbuffer.data(idx_op_fd + i * 60 + 27);

        auto to_xyz_yz = pbuffer.data(idx_op_fd + i * 60 + 28);

        auto to_xyz_zz = pbuffer.data(idx_op_fd + i * 60 + 29);

        auto to_xzz_xx = pbuffer.data(idx_op_fd + i * 60 + 30);

        auto to_xzz_xy = pbuffer.data(idx_op_fd + i * 60 + 31);

        auto to_xzz_xz = pbuffer.data(idx_op_fd + i * 60 + 32);

        auto to_xzz_yy = pbuffer.data(idx_op_fd + i * 60 + 33);

        auto to_xzz_yz = pbuffer.data(idx_op_fd + i * 60 + 34);

        auto to_xzz_zz = pbuffer.data(idx_op_fd + i * 60 + 35);

        auto to_yyy_xx = pbuffer.data(idx_op_fd + i * 60 + 36);

        auto to_yyy_xy = pbuffer.data(idx_op_fd + i * 60 + 37);

        auto to_yyy_xz = pbuffer.data(idx_op_fd + i * 60 + 38);

        auto to_yyy_yy = pbuffer.data(idx_op_fd + i * 60 + 39);

        auto to_yyy_yz = pbuffer.data(idx_op_fd + i * 60 + 40);

        auto to_yyy_zz = pbuffer.data(idx_op_fd + i * 60 + 41);

        auto to_yyz_xx = pbuffer.data(idx_op_fd + i * 60 + 42);

        auto to_yyz_xy = pbuffer.data(idx_op_fd + i * 60 + 43);

        auto to_yyz_xz = pbuffer.data(idx_op_fd + i * 60 + 44);

        auto to_yyz_yy = pbuffer.data(idx_op_fd + i * 60 + 45);

        auto to_yyz_yz = pbuffer.data(idx_op_fd + i * 60 + 46);

        auto to_yyz_zz = pbuffer.data(idx_op_fd + i * 60 + 47);

        auto to_yzz_xx = pbuffer.data(idx_op_fd + i * 60 + 48);

        auto to_yzz_xy = pbuffer.data(idx_op_fd + i * 60 + 49);

        auto to_yzz_xz = pbuffer.data(idx_op_fd + i * 60 + 50);

        auto to_yzz_yy = pbuffer.data(idx_op_fd + i * 60 + 51);

        auto to_yzz_yz = pbuffer.data(idx_op_fd + i * 60 + 52);

        auto to_yzz_zz = pbuffer.data(idx_op_fd + i * 60 + 53);

        auto to_zzz_xx = pbuffer.data(idx_op_fd + i * 60 + 54);

        auto to_zzz_xy = pbuffer.data(idx_op_fd + i * 60 + 55);

        auto to_zzz_xz = pbuffer.data(idx_op_fd + i * 60 + 56);

        auto to_zzz_yy = pbuffer.data(idx_op_fd + i * 60 + 57);

        auto to_zzz_yz = pbuffer.data(idx_op_fd + i * 60 + 58);

        auto to_zzz_zz = pbuffer.data(idx_op_fd + i * 60 + 59);

        // Set up components of auxiliary buffer : FG

        auto to_xxx_xxxx = pbuffer.data(idx_op_fg + i * 150 + 0);

        auto to_xxx_xxxy = pbuffer.data(idx_op_fg + i * 150 + 1);

        auto to_xxx_xxxz = pbuffer.data(idx_op_fg + i * 150 + 2);

        auto to_xxx_xxyy = pbuffer.data(idx_op_fg + i * 150 + 3);

        auto to_xxx_xxyz = pbuffer.data(idx_op_fg + i * 150 + 4);

        auto to_xxx_xxzz = pbuffer.data(idx_op_fg + i * 150 + 5);

        auto to_xxx_xyyy = pbuffer.data(idx_op_fg + i * 150 + 6);

        auto to_xxx_xyyz = pbuffer.data(idx_op_fg + i * 150 + 7);

        auto to_xxx_xyzz = pbuffer.data(idx_op_fg + i * 150 + 8);

        auto to_xxx_xzzz = pbuffer.data(idx_op_fg + i * 150 + 9);

        auto to_xxx_yyyy = pbuffer.data(idx_op_fg + i * 150 + 10);

        auto to_xxx_yyyz = pbuffer.data(idx_op_fg + i * 150 + 11);

        auto to_xxx_yyzz = pbuffer.data(idx_op_fg + i * 150 + 12);

        auto to_xxx_yzzz = pbuffer.data(idx_op_fg + i * 150 + 13);

        auto to_xxx_zzzz = pbuffer.data(idx_op_fg + i * 150 + 14);

        auto to_xxy_xxxx = pbuffer.data(idx_op_fg + i * 150 + 15);

        auto to_xxy_xxxy = pbuffer.data(idx_op_fg + i * 150 + 16);

        auto to_xxy_xxxz = pbuffer.data(idx_op_fg + i * 150 + 17);

        auto to_xxy_xxyy = pbuffer.data(idx_op_fg + i * 150 + 18);

        auto to_xxy_xxyz = pbuffer.data(idx_op_fg + i * 150 + 19);

        auto to_xxy_xxzz = pbuffer.data(idx_op_fg + i * 150 + 20);

        auto to_xxy_xyyy = pbuffer.data(idx_op_fg + i * 150 + 21);

        auto to_xxy_xyyz = pbuffer.data(idx_op_fg + i * 150 + 22);

        auto to_xxy_xyzz = pbuffer.data(idx_op_fg + i * 150 + 23);

        auto to_xxy_xzzz = pbuffer.data(idx_op_fg + i * 150 + 24);

        auto to_xxy_yyyy = pbuffer.data(idx_op_fg + i * 150 + 25);

        auto to_xxy_yyyz = pbuffer.data(idx_op_fg + i * 150 + 26);

        auto to_xxy_yyzz = pbuffer.data(idx_op_fg + i * 150 + 27);

        auto to_xxy_yzzz = pbuffer.data(idx_op_fg + i * 150 + 28);

        auto to_xxy_zzzz = pbuffer.data(idx_op_fg + i * 150 + 29);

        auto to_xxz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 30);

        auto to_xxz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 31);

        auto to_xxz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 32);

        auto to_xxz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 33);

        auto to_xxz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 34);

        auto to_xxz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 35);

        auto to_xxz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 36);

        auto to_xxz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 37);

        auto to_xxz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 38);

        auto to_xxz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 39);

        auto to_xxz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 40);

        auto to_xxz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 41);

        auto to_xxz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 42);

        auto to_xxz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 43);

        auto to_xxz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 44);

        auto to_xyy_xxxx = pbuffer.data(idx_op_fg + i * 150 + 45);

        auto to_xyy_xxxy = pbuffer.data(idx_op_fg + i * 150 + 46);

        auto to_xyy_xxxz = pbuffer.data(idx_op_fg + i * 150 + 47);

        auto to_xyy_xxyy = pbuffer.data(idx_op_fg + i * 150 + 48);

        auto to_xyy_xxyz = pbuffer.data(idx_op_fg + i * 150 + 49);

        auto to_xyy_xxzz = pbuffer.data(idx_op_fg + i * 150 + 50);

        auto to_xyy_xyyy = pbuffer.data(idx_op_fg + i * 150 + 51);

        auto to_xyy_xyyz = pbuffer.data(idx_op_fg + i * 150 + 52);

        auto to_xyy_xyzz = pbuffer.data(idx_op_fg + i * 150 + 53);

        auto to_xyy_xzzz = pbuffer.data(idx_op_fg + i * 150 + 54);

        auto to_xyy_yyyy = pbuffer.data(idx_op_fg + i * 150 + 55);

        auto to_xyy_yyyz = pbuffer.data(idx_op_fg + i * 150 + 56);

        auto to_xyy_yyzz = pbuffer.data(idx_op_fg + i * 150 + 57);

        auto to_xyy_yzzz = pbuffer.data(idx_op_fg + i * 150 + 58);

        auto to_xyy_zzzz = pbuffer.data(idx_op_fg + i * 150 + 59);

        auto to_xyz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 60);

        auto to_xyz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 61);

        auto to_xyz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 62);

        auto to_xyz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 63);

        auto to_xyz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 64);

        auto to_xyz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 65);

        auto to_xyz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 66);

        auto to_xyz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 67);

        auto to_xyz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 68);

        auto to_xyz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 69);

        auto to_xyz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 70);

        auto to_xyz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 71);

        auto to_xyz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 72);

        auto to_xyz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 73);

        auto to_xyz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 74);

        auto to_xzz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 75);

        auto to_xzz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 76);

        auto to_xzz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 77);

        auto to_xzz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 78);

        auto to_xzz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 79);

        auto to_xzz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 80);

        auto to_xzz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 81);

        auto to_xzz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 82);

        auto to_xzz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 83);

        auto to_xzz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 84);

        auto to_xzz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 85);

        auto to_xzz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 86);

        auto to_xzz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 87);

        auto to_xzz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 88);

        auto to_xzz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 89);

        auto to_yyy_xxxx = pbuffer.data(idx_op_fg + i * 150 + 90);

        auto to_yyy_xxxy = pbuffer.data(idx_op_fg + i * 150 + 91);

        auto to_yyy_xxxz = pbuffer.data(idx_op_fg + i * 150 + 92);

        auto to_yyy_xxyy = pbuffer.data(idx_op_fg + i * 150 + 93);

        auto to_yyy_xxyz = pbuffer.data(idx_op_fg + i * 150 + 94);

        auto to_yyy_xxzz = pbuffer.data(idx_op_fg + i * 150 + 95);

        auto to_yyy_xyyy = pbuffer.data(idx_op_fg + i * 150 + 96);

        auto to_yyy_xyyz = pbuffer.data(idx_op_fg + i * 150 + 97);

        auto to_yyy_xyzz = pbuffer.data(idx_op_fg + i * 150 + 98);

        auto to_yyy_xzzz = pbuffer.data(idx_op_fg + i * 150 + 99);

        auto to_yyy_yyyy = pbuffer.data(idx_op_fg + i * 150 + 100);

        auto to_yyy_yyyz = pbuffer.data(idx_op_fg + i * 150 + 101);

        auto to_yyy_yyzz = pbuffer.data(idx_op_fg + i * 150 + 102);

        auto to_yyy_yzzz = pbuffer.data(idx_op_fg + i * 150 + 103);

        auto to_yyy_zzzz = pbuffer.data(idx_op_fg + i * 150 + 104);

        auto to_yyz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 105);

        auto to_yyz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 106);

        auto to_yyz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 107);

        auto to_yyz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 108);

        auto to_yyz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 109);

        auto to_yyz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 110);

        auto to_yyz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 111);

        auto to_yyz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 112);

        auto to_yyz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 113);

        auto to_yyz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 114);

        auto to_yyz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 115);

        auto to_yyz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 116);

        auto to_yyz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 117);

        auto to_yyz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 118);

        auto to_yyz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 119);

        auto to_yzz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 120);

        auto to_yzz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 121);

        auto to_yzz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 122);

        auto to_yzz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 123);

        auto to_yzz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 124);

        auto to_yzz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 125);

        auto to_yzz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 126);

        auto to_yzz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 127);

        auto to_yzz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 128);

        auto to_yzz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 129);

        auto to_yzz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 130);

        auto to_yzz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 131);

        auto to_yzz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 132);

        auto to_yzz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 133);

        auto to_yzz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 134);

        auto to_zzz_xxxx = pbuffer.data(idx_op_fg + i * 150 + 135);

        auto to_zzz_xxxy = pbuffer.data(idx_op_fg + i * 150 + 136);

        auto to_zzz_xxxz = pbuffer.data(idx_op_fg + i * 150 + 137);

        auto to_zzz_xxyy = pbuffer.data(idx_op_fg + i * 150 + 138);

        auto to_zzz_xxyz = pbuffer.data(idx_op_fg + i * 150 + 139);

        auto to_zzz_xxzz = pbuffer.data(idx_op_fg + i * 150 + 140);

        auto to_zzz_xyyy = pbuffer.data(idx_op_fg + i * 150 + 141);

        auto to_zzz_xyyz = pbuffer.data(idx_op_fg + i * 150 + 142);

        auto to_zzz_xyzz = pbuffer.data(idx_op_fg + i * 150 + 143);

        auto to_zzz_xzzz = pbuffer.data(idx_op_fg + i * 150 + 144);

        auto to_zzz_yyyy = pbuffer.data(idx_op_fg + i * 150 + 145);

        auto to_zzz_yyyz = pbuffer.data(idx_op_fg + i * 150 + 146);

        auto to_zzz_yyzz = pbuffer.data(idx_op_fg + i * 150 + 147);

        auto to_zzz_yzzz = pbuffer.data(idx_op_fg + i * 150 + 148);

        auto to_zzz_zzzz = pbuffer.data(idx_op_fg + i * 150 + 149);

        // Set up components of auxiliary buffer : HD

        auto to_xxxxx_xx = pbuffer.data(idx_op_hd + i * 126 + 0);

        auto to_xxxxx_xy = pbuffer.data(idx_op_hd + i * 126 + 1);

        auto to_xxxxx_xz = pbuffer.data(idx_op_hd + i * 126 + 2);

        auto to_xxxxx_yy = pbuffer.data(idx_op_hd + i * 126 + 3);

        auto to_xxxxx_yz = pbuffer.data(idx_op_hd + i * 126 + 4);

        auto to_xxxxx_zz = pbuffer.data(idx_op_hd + i * 126 + 5);

        auto to_xxxxy_xx = pbuffer.data(idx_op_hd + i * 126 + 6);

        auto to_xxxxy_xy = pbuffer.data(idx_op_hd + i * 126 + 7);

        auto to_xxxxy_xz = pbuffer.data(idx_op_hd + i * 126 + 8);

        auto to_xxxxy_yy = pbuffer.data(idx_op_hd + i * 126 + 9);

        auto to_xxxxy_yz = pbuffer.data(idx_op_hd + i * 126 + 10);

        auto to_xxxxy_zz = pbuffer.data(idx_op_hd + i * 126 + 11);

        auto to_xxxxz_xx = pbuffer.data(idx_op_hd + i * 126 + 12);

        auto to_xxxxz_xy = pbuffer.data(idx_op_hd + i * 126 + 13);

        auto to_xxxxz_xz = pbuffer.data(idx_op_hd + i * 126 + 14);

        auto to_xxxxz_yy = pbuffer.data(idx_op_hd + i * 126 + 15);

        auto to_xxxxz_yz = pbuffer.data(idx_op_hd + i * 126 + 16);

        auto to_xxxxz_zz = pbuffer.data(idx_op_hd + i * 126 + 17);

        auto to_xxxyy_xx = pbuffer.data(idx_op_hd + i * 126 + 18);

        auto to_xxxyy_xy = pbuffer.data(idx_op_hd + i * 126 + 19);

        auto to_xxxyy_xz = pbuffer.data(idx_op_hd + i * 126 + 20);

        auto to_xxxyy_yy = pbuffer.data(idx_op_hd + i * 126 + 21);

        auto to_xxxyy_yz = pbuffer.data(idx_op_hd + i * 126 + 22);

        auto to_xxxyy_zz = pbuffer.data(idx_op_hd + i * 126 + 23);

        auto to_xxxyz_xx = pbuffer.data(idx_op_hd + i * 126 + 24);

        auto to_xxxyz_xy = pbuffer.data(idx_op_hd + i * 126 + 25);

        auto to_xxxyz_xz = pbuffer.data(idx_op_hd + i * 126 + 26);

        auto to_xxxyz_yy = pbuffer.data(idx_op_hd + i * 126 + 27);

        auto to_xxxyz_yz = pbuffer.data(idx_op_hd + i * 126 + 28);

        auto to_xxxyz_zz = pbuffer.data(idx_op_hd + i * 126 + 29);

        auto to_xxxzz_xx = pbuffer.data(idx_op_hd + i * 126 + 30);

        auto to_xxxzz_xy = pbuffer.data(idx_op_hd + i * 126 + 31);

        auto to_xxxzz_xz = pbuffer.data(idx_op_hd + i * 126 + 32);

        auto to_xxxzz_yy = pbuffer.data(idx_op_hd + i * 126 + 33);

        auto to_xxxzz_yz = pbuffer.data(idx_op_hd + i * 126 + 34);

        auto to_xxxzz_zz = pbuffer.data(idx_op_hd + i * 126 + 35);

        auto to_xxyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 36);

        auto to_xxyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 37);

        auto to_xxyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 38);

        auto to_xxyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 39);

        auto to_xxyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 40);

        auto to_xxyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 41);

        auto to_xxyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 42);

        auto to_xxyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 43);

        auto to_xxyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 44);

        auto to_xxyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 45);

        auto to_xxyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 46);

        auto to_xxyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 47);

        auto to_xxyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 48);

        auto to_xxyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 49);

        auto to_xxyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 50);

        auto to_xxyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 51);

        auto to_xxyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 52);

        auto to_xxyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 53);

        auto to_xxzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 54);

        auto to_xxzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 55);

        auto to_xxzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 56);

        auto to_xxzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 57);

        auto to_xxzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 58);

        auto to_xxzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 59);

        auto to_xyyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 60);

        auto to_xyyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 61);

        auto to_xyyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 62);

        auto to_xyyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 63);

        auto to_xyyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 64);

        auto to_xyyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 65);

        auto to_xyyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 66);

        auto to_xyyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 67);

        auto to_xyyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 68);

        auto to_xyyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 69);

        auto to_xyyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 70);

        auto to_xyyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 71);

        auto to_xyyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 72);

        auto to_xyyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 73);

        auto to_xyyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 74);

        auto to_xyyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 75);

        auto to_xyyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 76);

        auto to_xyyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 77);

        auto to_xyzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 78);

        auto to_xyzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 79);

        auto to_xyzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 80);

        auto to_xyzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 81);

        auto to_xyzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 82);

        auto to_xyzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 83);

        auto to_xzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 84);

        auto to_xzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 85);

        auto to_xzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 86);

        auto to_xzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 87);

        auto to_xzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 88);

        auto to_xzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 89);

        auto to_yyyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 90);

        auto to_yyyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 91);

        auto to_yyyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 92);

        auto to_yyyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 93);

        auto to_yyyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 94);

        auto to_yyyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 95);

        auto to_yyyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 96);

        auto to_yyyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 97);

        auto to_yyyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 98);

        auto to_yyyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 99);

        auto to_yyyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 100);

        auto to_yyyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 101);

        auto to_yyyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 102);

        auto to_yyyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 103);

        auto to_yyyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 104);

        auto to_yyyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 105);

        auto to_yyyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 106);

        auto to_yyyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 107);

        auto to_yyzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 108);

        auto to_yyzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 109);

        auto to_yyzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 110);

        auto to_yyzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 111);

        auto to_yyzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 112);

        auto to_yyzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 113);

        auto to_yzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 114);

        auto to_yzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 115);

        auto to_yzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 116);

        auto to_yzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 117);

        auto to_yzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 118);

        auto to_yzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 119);

        auto to_zzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 120);

        auto to_zzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 121);

        auto to_zzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 122);

        auto to_zzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 123);

        auto to_zzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 124);

        auto to_zzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 125);

        // Set up components of auxiliary buffer : HG

        auto to_xxxxx_xxxx = pbuffer.data(idx_op_hg + i * 315 + 0);

        auto to_xxxxx_xxxy = pbuffer.data(idx_op_hg + i * 315 + 1);

        auto to_xxxxx_xxxz = pbuffer.data(idx_op_hg + i * 315 + 2);

        auto to_xxxxx_xxyy = pbuffer.data(idx_op_hg + i * 315 + 3);

        auto to_xxxxx_xxyz = pbuffer.data(idx_op_hg + i * 315 + 4);

        auto to_xxxxx_xxzz = pbuffer.data(idx_op_hg + i * 315 + 5);

        auto to_xxxxx_xyyy = pbuffer.data(idx_op_hg + i * 315 + 6);

        auto to_xxxxx_xyyz = pbuffer.data(idx_op_hg + i * 315 + 7);

        auto to_xxxxx_xyzz = pbuffer.data(idx_op_hg + i * 315 + 8);

        auto to_xxxxx_xzzz = pbuffer.data(idx_op_hg + i * 315 + 9);

        auto to_xxxxx_yyyy = pbuffer.data(idx_op_hg + i * 315 + 10);

        auto to_xxxxx_yyyz = pbuffer.data(idx_op_hg + i * 315 + 11);

        auto to_xxxxx_yyzz = pbuffer.data(idx_op_hg + i * 315 + 12);

        auto to_xxxxx_yzzz = pbuffer.data(idx_op_hg + i * 315 + 13);

        auto to_xxxxx_zzzz = pbuffer.data(idx_op_hg + i * 315 + 14);

        auto to_xxxxy_xxxx = pbuffer.data(idx_op_hg + i * 315 + 15);

        auto to_xxxxy_xxxy = pbuffer.data(idx_op_hg + i * 315 + 16);

        auto to_xxxxy_xxxz = pbuffer.data(idx_op_hg + i * 315 + 17);

        auto to_xxxxy_xxyy = pbuffer.data(idx_op_hg + i * 315 + 18);

        auto to_xxxxy_xxyz = pbuffer.data(idx_op_hg + i * 315 + 19);

        auto to_xxxxy_xxzz = pbuffer.data(idx_op_hg + i * 315 + 20);

        auto to_xxxxy_xyyy = pbuffer.data(idx_op_hg + i * 315 + 21);

        auto to_xxxxy_xyyz = pbuffer.data(idx_op_hg + i * 315 + 22);

        auto to_xxxxy_xyzz = pbuffer.data(idx_op_hg + i * 315 + 23);

        auto to_xxxxy_xzzz = pbuffer.data(idx_op_hg + i * 315 + 24);

        auto to_xxxxy_yyyy = pbuffer.data(idx_op_hg + i * 315 + 25);

        auto to_xxxxy_yyyz = pbuffer.data(idx_op_hg + i * 315 + 26);

        auto to_xxxxy_yyzz = pbuffer.data(idx_op_hg + i * 315 + 27);

        auto to_xxxxy_yzzz = pbuffer.data(idx_op_hg + i * 315 + 28);

        auto to_xxxxy_zzzz = pbuffer.data(idx_op_hg + i * 315 + 29);

        auto to_xxxxz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 30);

        auto to_xxxxz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 31);

        auto to_xxxxz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 32);

        auto to_xxxxz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 33);

        auto to_xxxxz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 34);

        auto to_xxxxz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 35);

        auto to_xxxxz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 36);

        auto to_xxxxz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 37);

        auto to_xxxxz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 38);

        auto to_xxxxz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 39);

        auto to_xxxxz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 40);

        auto to_xxxxz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 41);

        auto to_xxxxz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 42);

        auto to_xxxxz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 43);

        auto to_xxxxz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 44);

        auto to_xxxyy_xxxx = pbuffer.data(idx_op_hg + i * 315 + 45);

        auto to_xxxyy_xxxy = pbuffer.data(idx_op_hg + i * 315 + 46);

        auto to_xxxyy_xxxz = pbuffer.data(idx_op_hg + i * 315 + 47);

        auto to_xxxyy_xxyy = pbuffer.data(idx_op_hg + i * 315 + 48);

        auto to_xxxyy_xxyz = pbuffer.data(idx_op_hg + i * 315 + 49);

        auto to_xxxyy_xxzz = pbuffer.data(idx_op_hg + i * 315 + 50);

        auto to_xxxyy_xyyy = pbuffer.data(idx_op_hg + i * 315 + 51);

        auto to_xxxyy_xyyz = pbuffer.data(idx_op_hg + i * 315 + 52);

        auto to_xxxyy_xyzz = pbuffer.data(idx_op_hg + i * 315 + 53);

        auto to_xxxyy_xzzz = pbuffer.data(idx_op_hg + i * 315 + 54);

        auto to_xxxyy_yyyy = pbuffer.data(idx_op_hg + i * 315 + 55);

        auto to_xxxyy_yyyz = pbuffer.data(idx_op_hg + i * 315 + 56);

        auto to_xxxyy_yyzz = pbuffer.data(idx_op_hg + i * 315 + 57);

        auto to_xxxyy_yzzz = pbuffer.data(idx_op_hg + i * 315 + 58);

        auto to_xxxyy_zzzz = pbuffer.data(idx_op_hg + i * 315 + 59);

        auto to_xxxyz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 60);

        auto to_xxxyz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 61);

        auto to_xxxyz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 62);

        auto to_xxxyz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 63);

        auto to_xxxyz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 64);

        auto to_xxxyz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 65);

        auto to_xxxyz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 66);

        auto to_xxxyz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 67);

        auto to_xxxyz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 68);

        auto to_xxxyz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 69);

        auto to_xxxyz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 70);

        auto to_xxxyz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 71);

        auto to_xxxyz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 72);

        auto to_xxxyz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 73);

        auto to_xxxyz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 74);

        auto to_xxxzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 75);

        auto to_xxxzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 76);

        auto to_xxxzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 77);

        auto to_xxxzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 78);

        auto to_xxxzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 79);

        auto to_xxxzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 80);

        auto to_xxxzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 81);

        auto to_xxxzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 82);

        auto to_xxxzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 83);

        auto to_xxxzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 84);

        auto to_xxxzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 85);

        auto to_xxxzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 86);

        auto to_xxxzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 87);

        auto to_xxxzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 88);

        auto to_xxxzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 89);

        auto to_xxyyy_xxxx = pbuffer.data(idx_op_hg + i * 315 + 90);

        auto to_xxyyy_xxxy = pbuffer.data(idx_op_hg + i * 315 + 91);

        auto to_xxyyy_xxxz = pbuffer.data(idx_op_hg + i * 315 + 92);

        auto to_xxyyy_xxyy = pbuffer.data(idx_op_hg + i * 315 + 93);

        auto to_xxyyy_xxyz = pbuffer.data(idx_op_hg + i * 315 + 94);

        auto to_xxyyy_xxzz = pbuffer.data(idx_op_hg + i * 315 + 95);

        auto to_xxyyy_xyyy = pbuffer.data(idx_op_hg + i * 315 + 96);

        auto to_xxyyy_xyyz = pbuffer.data(idx_op_hg + i * 315 + 97);

        auto to_xxyyy_xyzz = pbuffer.data(idx_op_hg + i * 315 + 98);

        auto to_xxyyy_xzzz = pbuffer.data(idx_op_hg + i * 315 + 99);

        auto to_xxyyy_yyyy = pbuffer.data(idx_op_hg + i * 315 + 100);

        auto to_xxyyy_yyyz = pbuffer.data(idx_op_hg + i * 315 + 101);

        auto to_xxyyy_yyzz = pbuffer.data(idx_op_hg + i * 315 + 102);

        auto to_xxyyy_yzzz = pbuffer.data(idx_op_hg + i * 315 + 103);

        auto to_xxyyy_zzzz = pbuffer.data(idx_op_hg + i * 315 + 104);

        auto to_xxyyz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 105);

        auto to_xxyyz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 106);

        auto to_xxyyz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 107);

        auto to_xxyyz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 108);

        auto to_xxyyz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 109);

        auto to_xxyyz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 110);

        auto to_xxyyz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 111);

        auto to_xxyyz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 112);

        auto to_xxyyz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 113);

        auto to_xxyyz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 114);

        auto to_xxyyz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 115);

        auto to_xxyyz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 116);

        auto to_xxyyz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 117);

        auto to_xxyyz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 118);

        auto to_xxyyz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 119);

        auto to_xxyzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 120);

        auto to_xxyzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 121);

        auto to_xxyzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 122);

        auto to_xxyzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 123);

        auto to_xxyzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 124);

        auto to_xxyzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 125);

        auto to_xxyzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 126);

        auto to_xxyzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 127);

        auto to_xxyzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 128);

        auto to_xxyzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 129);

        auto to_xxyzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 130);

        auto to_xxyzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 131);

        auto to_xxyzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 132);

        auto to_xxyzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 133);

        auto to_xxyzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 134);

        auto to_xxzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 135);

        auto to_xxzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 136);

        auto to_xxzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 137);

        auto to_xxzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 138);

        auto to_xxzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 139);

        auto to_xxzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 140);

        auto to_xxzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 141);

        auto to_xxzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 142);

        auto to_xxzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 143);

        auto to_xxzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 144);

        auto to_xxzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 145);

        auto to_xxzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 146);

        auto to_xxzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 147);

        auto to_xxzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 148);

        auto to_xxzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 149);

        auto to_xyyyy_xxxx = pbuffer.data(idx_op_hg + i * 315 + 150);

        auto to_xyyyy_xxxy = pbuffer.data(idx_op_hg + i * 315 + 151);

        auto to_xyyyy_xxxz = pbuffer.data(idx_op_hg + i * 315 + 152);

        auto to_xyyyy_xxyy = pbuffer.data(idx_op_hg + i * 315 + 153);

        auto to_xyyyy_xxyz = pbuffer.data(idx_op_hg + i * 315 + 154);

        auto to_xyyyy_xxzz = pbuffer.data(idx_op_hg + i * 315 + 155);

        auto to_xyyyy_xyyy = pbuffer.data(idx_op_hg + i * 315 + 156);

        auto to_xyyyy_xyyz = pbuffer.data(idx_op_hg + i * 315 + 157);

        auto to_xyyyy_xyzz = pbuffer.data(idx_op_hg + i * 315 + 158);

        auto to_xyyyy_xzzz = pbuffer.data(idx_op_hg + i * 315 + 159);

        auto to_xyyyy_yyyy = pbuffer.data(idx_op_hg + i * 315 + 160);

        auto to_xyyyy_yyyz = pbuffer.data(idx_op_hg + i * 315 + 161);

        auto to_xyyyy_yyzz = pbuffer.data(idx_op_hg + i * 315 + 162);

        auto to_xyyyy_yzzz = pbuffer.data(idx_op_hg + i * 315 + 163);

        auto to_xyyyy_zzzz = pbuffer.data(idx_op_hg + i * 315 + 164);

        auto to_xyyyz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 165);

        auto to_xyyyz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 166);

        auto to_xyyyz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 167);

        auto to_xyyyz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 168);

        auto to_xyyyz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 169);

        auto to_xyyyz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 170);

        auto to_xyyyz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 171);

        auto to_xyyyz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 172);

        auto to_xyyyz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 173);

        auto to_xyyyz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 174);

        auto to_xyyyz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 175);

        auto to_xyyyz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 176);

        auto to_xyyyz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 177);

        auto to_xyyyz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 178);

        auto to_xyyyz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 179);

        auto to_xyyzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 180);

        auto to_xyyzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 181);

        auto to_xyyzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 182);

        auto to_xyyzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 183);

        auto to_xyyzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 184);

        auto to_xyyzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 185);

        auto to_xyyzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 186);

        auto to_xyyzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 187);

        auto to_xyyzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 188);

        auto to_xyyzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 189);

        auto to_xyyzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 190);

        auto to_xyyzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 191);

        auto to_xyyzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 192);

        auto to_xyyzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 193);

        auto to_xyyzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 194);

        auto to_xyzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 195);

        auto to_xyzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 196);

        auto to_xyzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 197);

        auto to_xyzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 198);

        auto to_xyzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 199);

        auto to_xyzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 200);

        auto to_xyzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 201);

        auto to_xyzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 202);

        auto to_xyzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 203);

        auto to_xyzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 204);

        auto to_xyzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 205);

        auto to_xyzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 206);

        auto to_xyzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 207);

        auto to_xyzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 208);

        auto to_xyzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 209);

        auto to_xzzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 210);

        auto to_xzzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 211);

        auto to_xzzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 212);

        auto to_xzzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 213);

        auto to_xzzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 214);

        auto to_xzzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 215);

        auto to_xzzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 216);

        auto to_xzzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 217);

        auto to_xzzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 218);

        auto to_xzzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 219);

        auto to_xzzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 220);

        auto to_xzzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 221);

        auto to_xzzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 222);

        auto to_xzzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 223);

        auto to_xzzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 224);

        auto to_yyyyy_xxxx = pbuffer.data(idx_op_hg + i * 315 + 225);

        auto to_yyyyy_xxxy = pbuffer.data(idx_op_hg + i * 315 + 226);

        auto to_yyyyy_xxxz = pbuffer.data(idx_op_hg + i * 315 + 227);

        auto to_yyyyy_xxyy = pbuffer.data(idx_op_hg + i * 315 + 228);

        auto to_yyyyy_xxyz = pbuffer.data(idx_op_hg + i * 315 + 229);

        auto to_yyyyy_xxzz = pbuffer.data(idx_op_hg + i * 315 + 230);

        auto to_yyyyy_xyyy = pbuffer.data(idx_op_hg + i * 315 + 231);

        auto to_yyyyy_xyyz = pbuffer.data(idx_op_hg + i * 315 + 232);

        auto to_yyyyy_xyzz = pbuffer.data(idx_op_hg + i * 315 + 233);

        auto to_yyyyy_xzzz = pbuffer.data(idx_op_hg + i * 315 + 234);

        auto to_yyyyy_yyyy = pbuffer.data(idx_op_hg + i * 315 + 235);

        auto to_yyyyy_yyyz = pbuffer.data(idx_op_hg + i * 315 + 236);

        auto to_yyyyy_yyzz = pbuffer.data(idx_op_hg + i * 315 + 237);

        auto to_yyyyy_yzzz = pbuffer.data(idx_op_hg + i * 315 + 238);

        auto to_yyyyy_zzzz = pbuffer.data(idx_op_hg + i * 315 + 239);

        auto to_yyyyz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 240);

        auto to_yyyyz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 241);

        auto to_yyyyz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 242);

        auto to_yyyyz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 243);

        auto to_yyyyz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 244);

        auto to_yyyyz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 245);

        auto to_yyyyz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 246);

        auto to_yyyyz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 247);

        auto to_yyyyz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 248);

        auto to_yyyyz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 249);

        auto to_yyyyz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 250);

        auto to_yyyyz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 251);

        auto to_yyyyz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 252);

        auto to_yyyyz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 253);

        auto to_yyyyz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 254);

        auto to_yyyzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 255);

        auto to_yyyzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 256);

        auto to_yyyzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 257);

        auto to_yyyzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 258);

        auto to_yyyzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 259);

        auto to_yyyzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 260);

        auto to_yyyzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 261);

        auto to_yyyzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 262);

        auto to_yyyzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 263);

        auto to_yyyzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 264);

        auto to_yyyzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 265);

        auto to_yyyzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 266);

        auto to_yyyzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 267);

        auto to_yyyzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 268);

        auto to_yyyzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 269);

        auto to_yyzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 270);

        auto to_yyzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 271);

        auto to_yyzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 272);

        auto to_yyzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 273);

        auto to_yyzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 274);

        auto to_yyzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 275);

        auto to_yyzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 276);

        auto to_yyzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 277);

        auto to_yyzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 278);

        auto to_yyzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 279);

        auto to_yyzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 280);

        auto to_yyzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 281);

        auto to_yyzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 282);

        auto to_yyzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 283);

        auto to_yyzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 284);

        auto to_yzzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 285);

        auto to_yzzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 286);

        auto to_yzzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 287);

        auto to_yzzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 288);

        auto to_yzzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 289);

        auto to_yzzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 290);

        auto to_yzzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 291);

        auto to_yzzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 292);

        auto to_yzzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 293);

        auto to_yzzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 294);

        auto to_yzzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 295);

        auto to_yzzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 296);

        auto to_yzzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 297);

        auto to_yzzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 298);

        auto to_yzzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 299);

        auto to_zzzzz_xxxx = pbuffer.data(idx_op_hg + i * 315 + 300);

        auto to_zzzzz_xxxy = pbuffer.data(idx_op_hg + i * 315 + 301);

        auto to_zzzzz_xxxz = pbuffer.data(idx_op_hg + i * 315 + 302);

        auto to_zzzzz_xxyy = pbuffer.data(idx_op_hg + i * 315 + 303);

        auto to_zzzzz_xxyz = pbuffer.data(idx_op_hg + i * 315 + 304);

        auto to_zzzzz_xxzz = pbuffer.data(idx_op_hg + i * 315 + 305);

        auto to_zzzzz_xyyy = pbuffer.data(idx_op_hg + i * 315 + 306);

        auto to_zzzzz_xyyz = pbuffer.data(idx_op_hg + i * 315 + 307);

        auto to_zzzzz_xyzz = pbuffer.data(idx_op_hg + i * 315 + 308);

        auto to_zzzzz_xzzz = pbuffer.data(idx_op_hg + i * 315 + 309);

        auto to_zzzzz_yyyy = pbuffer.data(idx_op_hg + i * 315 + 310);

        auto to_zzzzz_yyyz = pbuffer.data(idx_op_hg + i * 315 + 311);

        auto to_zzzzz_yyzz = pbuffer.data(idx_op_hg + i * 315 + 312);

        auto to_zzzzz_yzzz = pbuffer.data(idx_op_hg + i * 315 + 313);

        auto to_zzzzz_zzzz = pbuffer.data(idx_op_hg + i * 315 + 314);

        // Set up 0-10 components of targeted buffer : GF

        auto to_x_x_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 0);

        auto to_x_x_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 1);

        auto to_x_x_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 2);

        auto to_x_x_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 3);

        auto to_x_x_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 4);

        auto to_x_x_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 5);

        auto to_x_x_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 6);

        auto to_x_x_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 7);

        auto to_x_x_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 8);

        auto to_x_x_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_x_x_xxxx_xxx, to_x_x_xxxx_xxy, to_x_x_xxxx_xxz, to_x_x_xxxx_xyy, to_x_x_xxxx_xyz, to_x_x_xxxx_xzz, to_x_x_xxxx_yyy, to_x_x_xxxx_yyz, to_x_x_xxxx_yzz, to_x_x_xxxx_zzz, to_xxx_xx, to_xxx_xxxx, to_xxx_xxxy, to_xxx_xxxz, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yz, to_xxx_zz, to_xxxxx_xx, to_xxxxx_xxxx, to_xxxxx_xxxy, to_xxxxx_xxxz, to_xxxxx_xxyy, to_xxxxx_xxyz, to_xxxxx_xxzz, to_xxxxx_xy, to_xxxxx_xyyy, to_xxxxx_xyyz, to_xxxxx_xyzz, to_xxxxx_xz, to_xxxxx_xzzz, to_xxxxx_yy, to_xxxxx_yz, to_xxxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxx_xxx[k] = 12.0 * to_xxx_xx[k] - 8.0 * to_xxx_xxxx[k] * tke_0 - 6.0 * to_xxxxx_xx[k] * tbe_0 + 4.0 * to_xxxxx_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxy[k] = 8.0 * to_xxx_xy[k] - 8.0 * to_xxx_xxxy[k] * tke_0 - 4.0 * to_xxxxx_xy[k] * tbe_0 + 4.0 * to_xxxxx_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xxz[k] = 8.0 * to_xxx_xz[k] - 8.0 * to_xxx_xxxz[k] * tke_0 - 4.0 * to_xxxxx_xz[k] * tbe_0 + 4.0 * to_xxxxx_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xyy[k] = 4.0 * to_xxx_yy[k] - 8.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxxx_yy[k] * tbe_0 + 4.0 * to_xxxxx_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xyz[k] = 4.0 * to_xxx_yz[k] - 8.0 * to_xxx_xxyz[k] * tke_0 - 2.0 * to_xxxxx_yz[k] * tbe_0 + 4.0 * to_xxxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xzz[k] = 4.0 * to_xxx_zz[k] - 8.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxxx_zz[k] * tbe_0 + 4.0 * to_xxxxx_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yyy[k] = -8.0 * to_xxx_xyyy[k] * tke_0 + 4.0 * to_xxxxx_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yyz[k] = -8.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yzz[k] = -8.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_zzz[k] = -8.0 * to_xxx_xzzz[k] * tke_0 + 4.0 * to_xxxxx_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : GF

        auto to_x_x_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 10);

        auto to_x_x_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 11);

        auto to_x_x_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 12);

        auto to_x_x_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 13);

        auto to_x_x_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 14);

        auto to_x_x_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 15);

        auto to_x_x_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 16);

        auto to_x_x_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 17);

        auto to_x_x_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 18);

        auto to_x_x_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_x_x_xxxy_xxx, to_x_x_xxxy_xxy, to_x_x_xxxy_xxz, to_x_x_xxxy_xyy, to_x_x_xxxy_xyz, to_x_x_xxxy_xzz, to_x_x_xxxy_yyy, to_x_x_xxxy_yyz, to_x_x_xxxy_yzz, to_x_x_xxxy_zzz, to_xxxxy_xx, to_xxxxy_xxxx, to_xxxxy_xxxy, to_xxxxy_xxxz, to_xxxxy_xxyy, to_xxxxy_xxyz, to_xxxxy_xxzz, to_xxxxy_xy, to_xxxxy_xyyy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_xzzz, to_xxxxy_yy, to_xxxxy_yz, to_xxxxy_zz, to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxy_xxx[k] = 9.0 * to_xxy_xx[k] - 6.0 * to_xxy_xxxx[k] * tke_0 - 6.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxy[k] = 6.0 * to_xxy_xy[k] - 6.0 * to_xxy_xxxy[k] * tke_0 - 4.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xxz[k] = 6.0 * to_xxy_xz[k] - 6.0 * to_xxy_xxxz[k] * tke_0 - 4.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xyy[k] = 3.0 * to_xxy_yy[k] - 6.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xyz[k] = 3.0 * to_xxy_yz[k] - 6.0 * to_xxy_xxyz[k] * tke_0 - 2.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xzz[k] = 3.0 * to_xxy_zz[k] - 6.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yyy[k] = -6.0 * to_xxy_xyyy[k] * tke_0 + 4.0 * to_xxxxy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yyz[k] = -6.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yzz[k] = -6.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_zzz[k] = -6.0 * to_xxy_xzzz[k] * tke_0 + 4.0 * to_xxxxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : GF

        auto to_x_x_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 20);

        auto to_x_x_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 21);

        auto to_x_x_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 22);

        auto to_x_x_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 23);

        auto to_x_x_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 24);

        auto to_x_x_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 25);

        auto to_x_x_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 26);

        auto to_x_x_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 27);

        auto to_x_x_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 28);

        auto to_x_x_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_x_xxxz_xxx, to_x_x_xxxz_xxy, to_x_x_xxxz_xxz, to_x_x_xxxz_xyy, to_x_x_xxxz_xyz, to_x_x_xxxz_xzz, to_x_x_xxxz_yyy, to_x_x_xxxz_yyz, to_x_x_xxxz_yzz, to_x_x_xxxz_zzz, to_xxxxz_xx, to_xxxxz_xxxx, to_xxxxz_xxxy, to_xxxxz_xxxz, to_xxxxz_xxyy, to_xxxxz_xxyz, to_xxxxz_xxzz, to_xxxxz_xy, to_xxxxz_xyyy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_xzzz, to_xxxxz_yy, to_xxxxz_yz, to_xxxxz_zz, to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxz_xxx[k] = 9.0 * to_xxz_xx[k] - 6.0 * to_xxz_xxxx[k] * tke_0 - 6.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxy[k] = 6.0 * to_xxz_xy[k] - 6.0 * to_xxz_xxxy[k] * tke_0 - 4.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xxz[k] = 6.0 * to_xxz_xz[k] - 6.0 * to_xxz_xxxz[k] * tke_0 - 4.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xyy[k] = 3.0 * to_xxz_yy[k] - 6.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xyz[k] = 3.0 * to_xxz_yz[k] - 6.0 * to_xxz_xxyz[k] * tke_0 - 2.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xzz[k] = 3.0 * to_xxz_zz[k] - 6.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yyy[k] = -6.0 * to_xxz_xyyy[k] * tke_0 + 4.0 * to_xxxxz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yyz[k] = -6.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yzz[k] = -6.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_zzz[k] = -6.0 * to_xxz_xzzz[k] * tke_0 + 4.0 * to_xxxxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : GF

        auto to_x_x_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 30);

        auto to_x_x_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 31);

        auto to_x_x_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 32);

        auto to_x_x_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 33);

        auto to_x_x_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 34);

        auto to_x_x_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 35);

        auto to_x_x_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 36);

        auto to_x_x_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 37);

        auto to_x_x_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 38);

        auto to_x_x_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_x_x_xxyy_xxx, to_x_x_xxyy_xxy, to_x_x_xxyy_xxz, to_x_x_xxyy_xyy, to_x_x_xxyy_xyz, to_x_x_xxyy_xzz, to_x_x_xxyy_yyy, to_x_x_xxyy_yyz, to_x_x_xxyy_yzz, to_x_x_xxyy_zzz, to_xxxyy_xx, to_xxxyy_xxxx, to_xxxyy_xxxy, to_xxxyy_xxxz, to_xxxyy_xxyy, to_xxxyy_xxyz, to_xxxyy_xxzz, to_xxxyy_xy, to_xxxyy_xyyy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_xzzz, to_xxxyy_yy, to_xxxyy_yz, to_xxxyy_zz, to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyy_xxx[k] = 6.0 * to_xyy_xx[k] - 4.0 * to_xyy_xxxx[k] * tke_0 - 6.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxy[k] = 4.0 * to_xyy_xy[k] - 4.0 * to_xyy_xxxy[k] * tke_0 - 4.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xxz[k] = 4.0 * to_xyy_xz[k] - 4.0 * to_xyy_xxxz[k] * tke_0 - 4.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xyy[k] = 2.0 * to_xyy_yy[k] - 4.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xyz[k] = 2.0 * to_xyy_yz[k] - 4.0 * to_xyy_xxyz[k] * tke_0 - 2.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xzz[k] = 2.0 * to_xyy_zz[k] - 4.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yyy[k] = -4.0 * to_xyy_xyyy[k] * tke_0 + 4.0 * to_xxxyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yyz[k] = -4.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yzz[k] = -4.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_zzz[k] = -4.0 * to_xyy_xzzz[k] * tke_0 + 4.0 * to_xxxyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : GF

        auto to_x_x_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 40);

        auto to_x_x_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 41);

        auto to_x_x_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 42);

        auto to_x_x_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 43);

        auto to_x_x_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 44);

        auto to_x_x_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 45);

        auto to_x_x_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 46);

        auto to_x_x_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 47);

        auto to_x_x_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 48);

        auto to_x_x_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_x_x_xxyz_xxx, to_x_x_xxyz_xxy, to_x_x_xxyz_xxz, to_x_x_xxyz_xyy, to_x_x_xxyz_xyz, to_x_x_xxyz_xzz, to_x_x_xxyz_yyy, to_x_x_xxyz_yyz, to_x_x_xxyz_yzz, to_x_x_xxyz_zzz, to_xxxyz_xx, to_xxxyz_xxxx, to_xxxyz_xxxy, to_xxxyz_xxxz, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yz, to_xxxyz_zz, to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyz_xxx[k] = 6.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxxx[k] * tke_0 - 6.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xxxy[k] * tke_0 - 4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xxz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xxxz[k] * tke_0 - 4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xyy[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xyz[k] = 2.0 * to_xyz_yz[k] - 4.0 * to_xyz_xxyz[k] * tke_0 - 2.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yyy[k] = -4.0 * to_xyz_xyyy[k] * tke_0 + 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yyz[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_zzz[k] = -4.0 * to_xyz_xzzz[k] * tke_0 + 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : GF

        auto to_x_x_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 50);

        auto to_x_x_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 51);

        auto to_x_x_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 52);

        auto to_x_x_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 53);

        auto to_x_x_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 54);

        auto to_x_x_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 55);

        auto to_x_x_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 56);

        auto to_x_x_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 57);

        auto to_x_x_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 58);

        auto to_x_x_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_x_xxzz_xxx, to_x_x_xxzz_xxy, to_x_x_xxzz_xxz, to_x_x_xxzz_xyy, to_x_x_xxzz_xyz, to_x_x_xxzz_xzz, to_x_x_xxzz_yyy, to_x_x_xxzz_yyz, to_x_x_xxzz_yzz, to_x_x_xxzz_zzz, to_xxxzz_xx, to_xxxzz_xxxx, to_xxxzz_xxxy, to_xxxzz_xxxz, to_xxxzz_xxyy, to_xxxzz_xxyz, to_xxxzz_xxzz, to_xxxzz_xy, to_xxxzz_xyyy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_xzzz, to_xxxzz_yy, to_xxxzz_yz, to_xxxzz_zz, to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxzz_xxx[k] = 6.0 * to_xzz_xx[k] - 4.0 * to_xzz_xxxx[k] * tke_0 - 6.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxy[k] = 4.0 * to_xzz_xy[k] - 4.0 * to_xzz_xxxy[k] * tke_0 - 4.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xxz[k] = 4.0 * to_xzz_xz[k] - 4.0 * to_xzz_xxxz[k] * tke_0 - 4.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xyy[k] = 2.0 * to_xzz_yy[k] - 4.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xyz[k] = 2.0 * to_xzz_yz[k] - 4.0 * to_xzz_xxyz[k] * tke_0 - 2.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xzz[k] = 2.0 * to_xzz_zz[k] - 4.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yyy[k] = -4.0 * to_xzz_xyyy[k] * tke_0 + 4.0 * to_xxxzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yyz[k] = -4.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yzz[k] = -4.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_zzz[k] = -4.0 * to_xzz_xzzz[k] * tke_0 + 4.0 * to_xxxzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : GF

        auto to_x_x_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 60);

        auto to_x_x_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 61);

        auto to_x_x_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 62);

        auto to_x_x_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 63);

        auto to_x_x_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 64);

        auto to_x_x_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 65);

        auto to_x_x_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 66);

        auto to_x_x_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 67);

        auto to_x_x_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 68);

        auto to_x_x_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_x_x_xyyy_xxx, to_x_x_xyyy_xxy, to_x_x_xyyy_xxz, to_x_x_xyyy_xyy, to_x_x_xyyy_xyz, to_x_x_xyyy_xzz, to_x_x_xyyy_yyy, to_x_x_xyyy_yyz, to_x_x_xyyy_yzz, to_x_x_xyyy_zzz, to_xxyyy_xx, to_xxyyy_xxxx, to_xxyyy_xxxy, to_xxyyy_xxxz, to_xxyyy_xxyy, to_xxyyy_xxyz, to_xxyyy_xxzz, to_xxyyy_xy, to_xxyyy_xyyy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_xzzz, to_xxyyy_yy, to_xxyyy_yz, to_xxyyy_zz, to_yyy_xx, to_yyy_xxxx, to_yyy_xxxy, to_yyy_xxxz, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyy_xxx[k] = 3.0 * to_yyy_xx[k] - 2.0 * to_yyy_xxxx[k] * tke_0 - 6.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxy[k] = 2.0 * to_yyy_xy[k] - 2.0 * to_yyy_xxxy[k] * tke_0 - 4.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xxz[k] = 2.0 * to_yyy_xz[k] - 2.0 * to_yyy_xxxz[k] * tke_0 - 4.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xyy[k] = to_yyy_yy[k] - 2.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xyz[k] = to_yyy_yz[k] - 2.0 * to_yyy_xxyz[k] * tke_0 - 2.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xzz[k] = to_yyy_zz[k] - 2.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yyy[k] = -2.0 * to_yyy_xyyy[k] * tke_0 + 4.0 * to_xxyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yyz[k] = -2.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yzz[k] = -2.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_zzz[k] = -2.0 * to_yyy_xzzz[k] * tke_0 + 4.0 * to_xxyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : GF

        auto to_x_x_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 70);

        auto to_x_x_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 71);

        auto to_x_x_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 72);

        auto to_x_x_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 73);

        auto to_x_x_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 74);

        auto to_x_x_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 75);

        auto to_x_x_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 76);

        auto to_x_x_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 77);

        auto to_x_x_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 78);

        auto to_x_x_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_x_x_xyyz_xxx, to_x_x_xyyz_xxy, to_x_x_xyyz_xxz, to_x_x_xyyz_xyy, to_x_x_xyyz_xyz, to_x_x_xyyz_xzz, to_x_x_xyyz_yyy, to_x_x_xyyz_yyz, to_x_x_xyyz_yzz, to_x_x_xyyz_zzz, to_xxyyz_xx, to_xxyyz_xxxx, to_xxyyz_xxxy, to_xxyyz_xxxz, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yz, to_xxyyz_zz, to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyz_xxx[k] = 3.0 * to_yyz_xx[k] - 2.0 * to_yyz_xxxx[k] * tke_0 - 6.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxy[k] = 2.0 * to_yyz_xy[k] - 2.0 * to_yyz_xxxy[k] * tke_0 - 4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xxz[k] = 2.0 * to_yyz_xz[k] - 2.0 * to_yyz_xxxz[k] * tke_0 - 4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xyy[k] = to_yyz_yy[k] - 2.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xyz[k] = to_yyz_yz[k] - 2.0 * to_yyz_xxyz[k] * tke_0 - 2.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xzz[k] = to_yyz_zz[k] - 2.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yyy[k] = -2.0 * to_yyz_xyyy[k] * tke_0 + 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yyz[k] = -2.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yzz[k] = -2.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_zzz[k] = -2.0 * to_yyz_xzzz[k] * tke_0 + 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : GF

        auto to_x_x_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 80);

        auto to_x_x_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 81);

        auto to_x_x_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 82);

        auto to_x_x_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 83);

        auto to_x_x_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 84);

        auto to_x_x_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 85);

        auto to_x_x_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 86);

        auto to_x_x_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 87);

        auto to_x_x_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 88);

        auto to_x_x_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_x_xyzz_xxx, to_x_x_xyzz_xxy, to_x_x_xyzz_xxz, to_x_x_xyzz_xyy, to_x_x_xyzz_xyz, to_x_x_xyzz_xzz, to_x_x_xyzz_yyy, to_x_x_xyzz_yyz, to_x_x_xyzz_yzz, to_x_x_xyzz_zzz, to_xxyzz_xx, to_xxyzz_xxxx, to_xxyzz_xxxy, to_xxyzz_xxxz, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yz, to_xxyzz_zz, to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyzz_xxx[k] = 3.0 * to_yzz_xx[k] - 2.0 * to_yzz_xxxx[k] * tke_0 - 6.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxy[k] = 2.0 * to_yzz_xy[k] - 2.0 * to_yzz_xxxy[k] * tke_0 - 4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xxz[k] = 2.0 * to_yzz_xz[k] - 2.0 * to_yzz_xxxz[k] * tke_0 - 4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xyy[k] = to_yzz_yy[k] - 2.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xyz[k] = to_yzz_yz[k] - 2.0 * to_yzz_xxyz[k] * tke_0 - 2.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xzz[k] = to_yzz_zz[k] - 2.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yyy[k] = -2.0 * to_yzz_xyyy[k] * tke_0 + 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yyz[k] = -2.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yzz[k] = -2.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_zzz[k] = -2.0 * to_yzz_xzzz[k] * tke_0 + 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-100 components of targeted buffer : GF

        auto to_x_x_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 90);

        auto to_x_x_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 91);

        auto to_x_x_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 92);

        auto to_x_x_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 93);

        auto to_x_x_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 94);

        auto to_x_x_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 95);

        auto to_x_x_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 96);

        auto to_x_x_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 97);

        auto to_x_x_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 98);

        auto to_x_x_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_x_x_xzzz_xxx, to_x_x_xzzz_xxy, to_x_x_xzzz_xxz, to_x_x_xzzz_xyy, to_x_x_xzzz_xyz, to_x_x_xzzz_xzz, to_x_x_xzzz_yyy, to_x_x_xzzz_yyz, to_x_x_xzzz_yzz, to_x_x_xzzz_zzz, to_xxzzz_xx, to_xxzzz_xxxx, to_xxzzz_xxxy, to_xxzzz_xxxz, to_xxzzz_xxyy, to_xxzzz_xxyz, to_xxzzz_xxzz, to_xxzzz_xy, to_xxzzz_xyyy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_xzzz, to_xxzzz_yy, to_xxzzz_yz, to_xxzzz_zz, to_zzz_xx, to_zzz_xxxx, to_zzz_xxxy, to_zzz_xxxz, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzzz_xxx[k] = 3.0 * to_zzz_xx[k] - 2.0 * to_zzz_xxxx[k] * tke_0 - 6.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxy[k] = 2.0 * to_zzz_xy[k] - 2.0 * to_zzz_xxxy[k] * tke_0 - 4.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xxz[k] = 2.0 * to_zzz_xz[k] - 2.0 * to_zzz_xxxz[k] * tke_0 - 4.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xyy[k] = to_zzz_yy[k] - 2.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xyz[k] = to_zzz_yz[k] - 2.0 * to_zzz_xxyz[k] * tke_0 - 2.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xzz[k] = to_zzz_zz[k] - 2.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yyy[k] = -2.0 * to_zzz_xyyy[k] * tke_0 + 4.0 * to_xxzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yyz[k] = -2.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yzz[k] = -2.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_zzz[k] = -2.0 * to_zzz_xzzz[k] * tke_0 + 4.0 * to_xxzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 100-110 components of targeted buffer : GF

        auto to_x_x_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 100);

        auto to_x_x_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 101);

        auto to_x_x_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 102);

        auto to_x_x_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 103);

        auto to_x_x_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 104);

        auto to_x_x_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 105);

        auto to_x_x_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 106);

        auto to_x_x_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 107);

        auto to_x_x_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 108);

        auto to_x_x_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_x_x_yyyy_xxx, to_x_x_yyyy_xxy, to_x_x_yyyy_xxz, to_x_x_yyyy_xyy, to_x_x_yyyy_xyz, to_x_x_yyyy_xzz, to_x_x_yyyy_yyy, to_x_x_yyyy_yyz, to_x_x_yyyy_yzz, to_x_x_yyyy_zzz, to_xyyyy_xx, to_xyyyy_xxxx, to_xyyyy_xxxy, to_xyyyy_xxxz, to_xyyyy_xxyy, to_xyyyy_xxyz, to_xyyyy_xxzz, to_xyyyy_xy, to_xyyyy_xyyy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_xzzz, to_xyyyy_yy, to_xyyyy_yz, to_xyyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyy_xxx[k] = -6.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxy[k] = -4.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xxz[k] = -4.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xyy[k] = -2.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xyz[k] = -2.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xzz[k] = -2.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yyy[k] = 4.0 * to_xyyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yyz[k] = 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yzz[k] = 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_zzz[k] = 4.0 * to_xyyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 110-120 components of targeted buffer : GF

        auto to_x_x_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 110);

        auto to_x_x_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 111);

        auto to_x_x_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 112);

        auto to_x_x_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 113);

        auto to_x_x_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 114);

        auto to_x_x_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 115);

        auto to_x_x_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 116);

        auto to_x_x_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 117);

        auto to_x_x_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 118);

        auto to_x_x_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_x_yyyz_xxx, to_x_x_yyyz_xxy, to_x_x_yyyz_xxz, to_x_x_yyyz_xyy, to_x_x_yyyz_xyz, to_x_x_yyyz_xzz, to_x_x_yyyz_yyy, to_x_x_yyyz_yyz, to_x_x_yyyz_yzz, to_x_x_yyyz_zzz, to_xyyyz_xx, to_xyyyz_xxxx, to_xyyyz_xxxy, to_xyyyz_xxxz, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yz, to_xyyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyz_xxx[k] = -6.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxy[k] = -4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xxz[k] = -4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xyy[k] = -2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xyz[k] = -2.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xzz[k] = -2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yyy[k] = 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yyz[k] = 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yzz[k] = 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_zzz[k] = 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-130 components of targeted buffer : GF

        auto to_x_x_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 120);

        auto to_x_x_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 121);

        auto to_x_x_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 122);

        auto to_x_x_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 123);

        auto to_x_x_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 124);

        auto to_x_x_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 125);

        auto to_x_x_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 126);

        auto to_x_x_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 127);

        auto to_x_x_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 128);

        auto to_x_x_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_x_x_yyzz_xxx, to_x_x_yyzz_xxy, to_x_x_yyzz_xxz, to_x_x_yyzz_xyy, to_x_x_yyzz_xyz, to_x_x_yyzz_xzz, to_x_x_yyzz_yyy, to_x_x_yyzz_yyz, to_x_x_yyzz_yzz, to_x_x_yyzz_zzz, to_xyyzz_xx, to_xyyzz_xxxx, to_xyyzz_xxxy, to_xyyzz_xxxz, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yz, to_xyyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyzz_xxx[k] = -6.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxy[k] = -4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xxz[k] = -4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xyy[k] = -2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xyz[k] = -2.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xzz[k] = -2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yyy[k] = 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yyz[k] = 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yzz[k] = 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_zzz[k] = 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 130-140 components of targeted buffer : GF

        auto to_x_x_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 130);

        auto to_x_x_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 131);

        auto to_x_x_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 132);

        auto to_x_x_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 133);

        auto to_x_x_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 134);

        auto to_x_x_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 135);

        auto to_x_x_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 136);

        auto to_x_x_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 137);

        auto to_x_x_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 138);

        auto to_x_x_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_x_x_yzzz_xxx, to_x_x_yzzz_xxy, to_x_x_yzzz_xxz, to_x_x_yzzz_xyy, to_x_x_yzzz_xyz, to_x_x_yzzz_xzz, to_x_x_yzzz_yyy, to_x_x_yzzz_yyz, to_x_x_yzzz_yzz, to_x_x_yzzz_zzz, to_xyzzz_xx, to_xyzzz_xxxx, to_xyzzz_xxxy, to_xyzzz_xxxz, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yz, to_xyzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzzz_xxx[k] = -6.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxy[k] = -4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xxz[k] = -4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xyy[k] = -2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xyz[k] = -2.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xzz[k] = -2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yyy[k] = 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yyz[k] = 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yzz[k] = 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_zzz[k] = 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 140-150 components of targeted buffer : GF

        auto to_x_x_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 140);

        auto to_x_x_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 141);

        auto to_x_x_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 142);

        auto to_x_x_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 143);

        auto to_x_x_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 144);

        auto to_x_x_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 145);

        auto to_x_x_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 146);

        auto to_x_x_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 147);

        auto to_x_x_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 148);

        auto to_x_x_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 0 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_x_zzzz_xxx, to_x_x_zzzz_xxy, to_x_x_zzzz_xxz, to_x_x_zzzz_xyy, to_x_x_zzzz_xyz, to_x_x_zzzz_xzz, to_x_x_zzzz_yyy, to_x_x_zzzz_yyz, to_x_x_zzzz_yzz, to_x_x_zzzz_zzz, to_xzzzz_xx, to_xzzzz_xxxx, to_xzzzz_xxxy, to_xzzzz_xxxz, to_xzzzz_xxyy, to_xzzzz_xxyz, to_xzzzz_xxzz, to_xzzzz_xy, to_xzzzz_xyyy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_xzzz, to_xzzzz_yy, to_xzzzz_yz, to_xzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzzz_xxx[k] = -6.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxy[k] = -4.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xxz[k] = -4.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xyy[k] = -2.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xyz[k] = -2.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xzz[k] = -2.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yyy[k] = 4.0 * to_xzzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yyz[k] = 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yzz[k] = 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_zzz[k] = 4.0 * to_xzzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-160 components of targeted buffer : GF

        auto to_x_y_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 0);

        auto to_x_y_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 1);

        auto to_x_y_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 2);

        auto to_x_y_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 3);

        auto to_x_y_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 4);

        auto to_x_y_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 5);

        auto to_x_y_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 6);

        auto to_x_y_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 7);

        auto to_x_y_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 8);

        auto to_x_y_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_x_y_xxxx_xxx, to_x_y_xxxx_xxy, to_x_y_xxxx_xxz, to_x_y_xxxx_xyy, to_x_y_xxxx_xyz, to_x_y_xxxx_xzz, to_x_y_xxxx_yyy, to_x_y_xxxx_yyz, to_x_y_xxxx_yzz, to_x_y_xxxx_zzz, to_xxx_xx, to_xxx_xxxy, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_yy, to_xxx_yyyy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxxxx_xx, to_xxxxx_xxxy, to_xxxxx_xxyy, to_xxxxx_xxyz, to_xxxxx_xy, to_xxxxx_xyyy, to_xxxxx_xyyz, to_xxxxx_xyzz, to_xxxxx_xz, to_xxxxx_yy, to_xxxxx_yyyy, to_xxxxx_yyyz, to_xxxxx_yyzz, to_xxxxx_yz, to_xxxxx_yzzz, to_xxxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxx_xxx[k] = -8.0 * to_xxx_xxxy[k] * tke_0 + 4.0 * to_xxxxx_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxy[k] = 4.0 * to_xxx_xx[k] - 8.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxxx_xx[k] * tbe_0 + 4.0 * to_xxxxx_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xxz[k] = -8.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xyy[k] = 8.0 * to_xxx_xy[k] - 8.0 * to_xxx_xyyy[k] * tke_0 - 4.0 * to_xxxxx_xy[k] * tbe_0 + 4.0 * to_xxxxx_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xyz[k] = 4.0 * to_xxx_xz[k] - 8.0 * to_xxx_xyyz[k] * tke_0 - 2.0 * to_xxxxx_xz[k] * tbe_0 + 4.0 * to_xxxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xzz[k] = -8.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yyy[k] = 12.0 * to_xxx_yy[k] - 8.0 * to_xxx_yyyy[k] * tke_0 - 6.0 * to_xxxxx_yy[k] * tbe_0 + 4.0 * to_xxxxx_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yyz[k] = 8.0 * to_xxx_yz[k] - 8.0 * to_xxx_yyyz[k] * tke_0 - 4.0 * to_xxxxx_yz[k] * tbe_0 + 4.0 * to_xxxxx_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yzz[k] = 4.0 * to_xxx_zz[k] - 8.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxxx_zz[k] * tbe_0 + 4.0 * to_xxxxx_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_zzz[k] = -8.0 * to_xxx_yzzz[k] * tke_0 + 4.0 * to_xxxxx_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 160-170 components of targeted buffer : GF

        auto to_x_y_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 10);

        auto to_x_y_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 11);

        auto to_x_y_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 12);

        auto to_x_y_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 13);

        auto to_x_y_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 14);

        auto to_x_y_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 15);

        auto to_x_y_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 16);

        auto to_x_y_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 17);

        auto to_x_y_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 18);

        auto to_x_y_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_x_y_xxxy_xxx, to_x_y_xxxy_xxy, to_x_y_xxxy_xxz, to_x_y_xxxy_xyy, to_x_y_xxxy_xyz, to_x_y_xxxy_xzz, to_x_y_xxxy_yyy, to_x_y_xxxy_yyz, to_x_y_xxxy_yzz, to_x_y_xxxy_zzz, to_xxxxy_xx, to_xxxxy_xxxy, to_xxxxy_xxyy, to_xxxxy_xxyz, to_xxxxy_xy, to_xxxxy_xyyy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_yy, to_xxxxy_yyyy, to_xxxxy_yyyz, to_xxxxy_yyzz, to_xxxxy_yz, to_xxxxy_yzzz, to_xxxxy_zz, to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxy_xxx[k] = -6.0 * to_xxy_xxxy[k] * tke_0 + 4.0 * to_xxxxy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxy[k] = 3.0 * to_xxy_xx[k] - 6.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xxz[k] = -6.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xyy[k] = 6.0 * to_xxy_xy[k] - 6.0 * to_xxy_xyyy[k] * tke_0 - 4.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xyz[k] = 3.0 * to_xxy_xz[k] - 6.0 * to_xxy_xyyz[k] * tke_0 - 2.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xzz[k] = -6.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yyy[k] = 9.0 * to_xxy_yy[k] - 6.0 * to_xxy_yyyy[k] * tke_0 - 6.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yyz[k] = 6.0 * to_xxy_yz[k] - 6.0 * to_xxy_yyyz[k] * tke_0 - 4.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yzz[k] = 3.0 * to_xxy_zz[k] - 6.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_zzz[k] = -6.0 * to_xxy_yzzz[k] * tke_0 + 4.0 * to_xxxxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 170-180 components of targeted buffer : GF

        auto to_x_y_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 20);

        auto to_x_y_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 21);

        auto to_x_y_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 22);

        auto to_x_y_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 23);

        auto to_x_y_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 24);

        auto to_x_y_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 25);

        auto to_x_y_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 26);

        auto to_x_y_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 27);

        auto to_x_y_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 28);

        auto to_x_y_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_y_xxxz_xxx, to_x_y_xxxz_xxy, to_x_y_xxxz_xxz, to_x_y_xxxz_xyy, to_x_y_xxxz_xyz, to_x_y_xxxz_xzz, to_x_y_xxxz_yyy, to_x_y_xxxz_yyz, to_x_y_xxxz_yzz, to_x_y_xxxz_zzz, to_xxxxz_xx, to_xxxxz_xxxy, to_xxxxz_xxyy, to_xxxxz_xxyz, to_xxxxz_xy, to_xxxxz_xyyy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_yy, to_xxxxz_yyyy, to_xxxxz_yyyz, to_xxxxz_yyzz, to_xxxxz_yz, to_xxxxz_yzzz, to_xxxxz_zz, to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxz_xxx[k] = -6.0 * to_xxz_xxxy[k] * tke_0 + 4.0 * to_xxxxz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxy[k] = 3.0 * to_xxz_xx[k] - 6.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xxz[k] = -6.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xyy[k] = 6.0 * to_xxz_xy[k] - 6.0 * to_xxz_xyyy[k] * tke_0 - 4.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xyz[k] = 3.0 * to_xxz_xz[k] - 6.0 * to_xxz_xyyz[k] * tke_0 - 2.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xzz[k] = -6.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yyy[k] = 9.0 * to_xxz_yy[k] - 6.0 * to_xxz_yyyy[k] * tke_0 - 6.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yyz[k] = 6.0 * to_xxz_yz[k] - 6.0 * to_xxz_yyyz[k] * tke_0 - 4.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yzz[k] = 3.0 * to_xxz_zz[k] - 6.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_zzz[k] = -6.0 * to_xxz_yzzz[k] * tke_0 + 4.0 * to_xxxxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-190 components of targeted buffer : GF

        auto to_x_y_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 30);

        auto to_x_y_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 31);

        auto to_x_y_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 32);

        auto to_x_y_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 33);

        auto to_x_y_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 34);

        auto to_x_y_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 35);

        auto to_x_y_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 36);

        auto to_x_y_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 37);

        auto to_x_y_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 38);

        auto to_x_y_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_x_y_xxyy_xxx, to_x_y_xxyy_xxy, to_x_y_xxyy_xxz, to_x_y_xxyy_xyy, to_x_y_xxyy_xyz, to_x_y_xxyy_xzz, to_x_y_xxyy_yyy, to_x_y_xxyy_yyz, to_x_y_xxyy_yzz, to_x_y_xxyy_zzz, to_xxxyy_xx, to_xxxyy_xxxy, to_xxxyy_xxyy, to_xxxyy_xxyz, to_xxxyy_xy, to_xxxyy_xyyy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_yy, to_xxxyy_yyyy, to_xxxyy_yyyz, to_xxxyy_yyzz, to_xxxyy_yz, to_xxxyy_yzzz, to_xxxyy_zz, to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyy_xxx[k] = -4.0 * to_xyy_xxxy[k] * tke_0 + 4.0 * to_xxxyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxy[k] = 2.0 * to_xyy_xx[k] - 4.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xxz[k] = -4.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xyy[k] = 4.0 * to_xyy_xy[k] - 4.0 * to_xyy_xyyy[k] * tke_0 - 4.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xyz[k] = 2.0 * to_xyy_xz[k] - 4.0 * to_xyy_xyyz[k] * tke_0 - 2.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xzz[k] = -4.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yyy[k] = 6.0 * to_xyy_yy[k] - 4.0 * to_xyy_yyyy[k] * tke_0 - 6.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yyz[k] = 4.0 * to_xyy_yz[k] - 4.0 * to_xyy_yyyz[k] * tke_0 - 4.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yzz[k] = 2.0 * to_xyy_zz[k] - 4.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_zzz[k] = -4.0 * to_xyy_yzzz[k] * tke_0 + 4.0 * to_xxxyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 190-200 components of targeted buffer : GF

        auto to_x_y_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 40);

        auto to_x_y_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 41);

        auto to_x_y_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 42);

        auto to_x_y_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 43);

        auto to_x_y_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 44);

        auto to_x_y_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 45);

        auto to_x_y_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 46);

        auto to_x_y_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 47);

        auto to_x_y_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 48);

        auto to_x_y_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_x_y_xxyz_xxx, to_x_y_xxyz_xxy, to_x_y_xxyz_xxz, to_x_y_xxyz_xyy, to_x_y_xxyz_xyz, to_x_y_xxyz_xzz, to_x_y_xxyz_yyy, to_x_y_xxyz_yyz, to_x_y_xxyz_yzz, to_x_y_xxyz_zzz, to_xxxyz_xx, to_xxxyz_xxxy, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_yy, to_xxxyz_yyyy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyz_xxx[k] = -4.0 * to_xyz_xxxy[k] * tke_0 + 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxy[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xxz[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xyy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyyy[k] * tke_0 - 4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xyz[k] = 2.0 * to_xyz_xz[k] - 4.0 * to_xyz_xyyz[k] * tke_0 - 2.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yyy[k] = 6.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyyy[k] * tke_0 - 6.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yyz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yyyz[k] * tke_0 - 4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_zzz[k] = -4.0 * to_xyz_yzzz[k] * tke_0 + 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 200-210 components of targeted buffer : GF

        auto to_x_y_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 50);

        auto to_x_y_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 51);

        auto to_x_y_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 52);

        auto to_x_y_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 53);

        auto to_x_y_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 54);

        auto to_x_y_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 55);

        auto to_x_y_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 56);

        auto to_x_y_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 57);

        auto to_x_y_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 58);

        auto to_x_y_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_y_xxzz_xxx, to_x_y_xxzz_xxy, to_x_y_xxzz_xxz, to_x_y_xxzz_xyy, to_x_y_xxzz_xyz, to_x_y_xxzz_xzz, to_x_y_xxzz_yyy, to_x_y_xxzz_yyz, to_x_y_xxzz_yzz, to_x_y_xxzz_zzz, to_xxxzz_xx, to_xxxzz_xxxy, to_xxxzz_xxyy, to_xxxzz_xxyz, to_xxxzz_xy, to_xxxzz_xyyy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_yy, to_xxxzz_yyyy, to_xxxzz_yyyz, to_xxxzz_yyzz, to_xxxzz_yz, to_xxxzz_yzzz, to_xxxzz_zz, to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxzz_xxx[k] = -4.0 * to_xzz_xxxy[k] * tke_0 + 4.0 * to_xxxzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxy[k] = 2.0 * to_xzz_xx[k] - 4.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xxz[k] = -4.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xyy[k] = 4.0 * to_xzz_xy[k] - 4.0 * to_xzz_xyyy[k] * tke_0 - 4.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xyz[k] = 2.0 * to_xzz_xz[k] - 4.0 * to_xzz_xyyz[k] * tke_0 - 2.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xzz[k] = -4.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yyy[k] = 6.0 * to_xzz_yy[k] - 4.0 * to_xzz_yyyy[k] * tke_0 - 6.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yyz[k] = 4.0 * to_xzz_yz[k] - 4.0 * to_xzz_yyyz[k] * tke_0 - 4.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yzz[k] = 2.0 * to_xzz_zz[k] - 4.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_zzz[k] = -4.0 * to_xzz_yzzz[k] * tke_0 + 4.0 * to_xxxzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-220 components of targeted buffer : GF

        auto to_x_y_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 60);

        auto to_x_y_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 61);

        auto to_x_y_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 62);

        auto to_x_y_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 63);

        auto to_x_y_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 64);

        auto to_x_y_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 65);

        auto to_x_y_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 66);

        auto to_x_y_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 67);

        auto to_x_y_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 68);

        auto to_x_y_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_x_y_xyyy_xxx, to_x_y_xyyy_xxy, to_x_y_xyyy_xxz, to_x_y_xyyy_xyy, to_x_y_xyyy_xyz, to_x_y_xyyy_xzz, to_x_y_xyyy_yyy, to_x_y_xyyy_yyz, to_x_y_xyyy_yzz, to_x_y_xyyy_zzz, to_xxyyy_xx, to_xxyyy_xxxy, to_xxyyy_xxyy, to_xxyyy_xxyz, to_xxyyy_xy, to_xxyyy_xyyy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_yy, to_xxyyy_yyyy, to_xxyyy_yyyz, to_xxyyy_yyzz, to_xxyyy_yz, to_xxyyy_yzzz, to_xxyyy_zz, to_yyy_xx, to_yyy_xxxy, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_yy, to_yyy_yyyy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyy_xxx[k] = -2.0 * to_yyy_xxxy[k] * tke_0 + 4.0 * to_xxyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxy[k] = to_yyy_xx[k] - 2.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xxz[k] = -2.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xyy[k] = 2.0 * to_yyy_xy[k] - 2.0 * to_yyy_xyyy[k] * tke_0 - 4.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xyz[k] = to_yyy_xz[k] - 2.0 * to_yyy_xyyz[k] * tke_0 - 2.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xzz[k] = -2.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yyy[k] = 3.0 * to_yyy_yy[k] - 2.0 * to_yyy_yyyy[k] * tke_0 - 6.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yyz[k] = 2.0 * to_yyy_yz[k] - 2.0 * to_yyy_yyyz[k] * tke_0 - 4.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yzz[k] = to_yyy_zz[k] - 2.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_zzz[k] = -2.0 * to_yyy_yzzz[k] * tke_0 + 4.0 * to_xxyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 220-230 components of targeted buffer : GF

        auto to_x_y_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 70);

        auto to_x_y_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 71);

        auto to_x_y_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 72);

        auto to_x_y_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 73);

        auto to_x_y_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 74);

        auto to_x_y_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 75);

        auto to_x_y_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 76);

        auto to_x_y_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 77);

        auto to_x_y_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 78);

        auto to_x_y_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_x_y_xyyz_xxx, to_x_y_xyyz_xxy, to_x_y_xyyz_xxz, to_x_y_xyyz_xyy, to_x_y_xyyz_xyz, to_x_y_xyyz_xzz, to_x_y_xyyz_yyy, to_x_y_xyyz_yyz, to_x_y_xyyz_yzz, to_x_y_xyyz_zzz, to_xxyyz_xx, to_xxyyz_xxxy, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_yy, to_xxyyz_yyyy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyz_xxx[k] = -2.0 * to_yyz_xxxy[k] * tke_0 + 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxy[k] = to_yyz_xx[k] - 2.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xxz[k] = -2.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xyy[k] = 2.0 * to_yyz_xy[k] - 2.0 * to_yyz_xyyy[k] * tke_0 - 4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xyz[k] = to_yyz_xz[k] - 2.0 * to_yyz_xyyz[k] * tke_0 - 2.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xzz[k] = -2.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yyy[k] = 3.0 * to_yyz_yy[k] - 2.0 * to_yyz_yyyy[k] * tke_0 - 6.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yyz[k] = 2.0 * to_yyz_yz[k] - 2.0 * to_yyz_yyyz[k] * tke_0 - 4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yzz[k] = to_yyz_zz[k] - 2.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_zzz[k] = -2.0 * to_yyz_yzzz[k] * tke_0 + 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 230-240 components of targeted buffer : GF

        auto to_x_y_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 80);

        auto to_x_y_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 81);

        auto to_x_y_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 82);

        auto to_x_y_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 83);

        auto to_x_y_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 84);

        auto to_x_y_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 85);

        auto to_x_y_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 86);

        auto to_x_y_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 87);

        auto to_x_y_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 88);

        auto to_x_y_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_y_xyzz_xxx, to_x_y_xyzz_xxy, to_x_y_xyzz_xxz, to_x_y_xyzz_xyy, to_x_y_xyzz_xyz, to_x_y_xyzz_xzz, to_x_y_xyzz_yyy, to_x_y_xyzz_yyz, to_x_y_xyzz_yzz, to_x_y_xyzz_zzz, to_xxyzz_xx, to_xxyzz_xxxy, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_yy, to_xxyzz_yyyy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyzz_xxx[k] = -2.0 * to_yzz_xxxy[k] * tke_0 + 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxy[k] = to_yzz_xx[k] - 2.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xxz[k] = -2.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xyy[k] = 2.0 * to_yzz_xy[k] - 2.0 * to_yzz_xyyy[k] * tke_0 - 4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xyz[k] = to_yzz_xz[k] - 2.0 * to_yzz_xyyz[k] * tke_0 - 2.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xzz[k] = -2.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yyy[k] = 3.0 * to_yzz_yy[k] - 2.0 * to_yzz_yyyy[k] * tke_0 - 6.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yyz[k] = 2.0 * to_yzz_yz[k] - 2.0 * to_yzz_yyyz[k] * tke_0 - 4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yzz[k] = to_yzz_zz[k] - 2.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_zzz[k] = -2.0 * to_yzz_yzzz[k] * tke_0 + 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-250 components of targeted buffer : GF

        auto to_x_y_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 90);

        auto to_x_y_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 91);

        auto to_x_y_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 92);

        auto to_x_y_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 93);

        auto to_x_y_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 94);

        auto to_x_y_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 95);

        auto to_x_y_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 96);

        auto to_x_y_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 97);

        auto to_x_y_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 98);

        auto to_x_y_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_x_y_xzzz_xxx, to_x_y_xzzz_xxy, to_x_y_xzzz_xxz, to_x_y_xzzz_xyy, to_x_y_xzzz_xyz, to_x_y_xzzz_xzz, to_x_y_xzzz_yyy, to_x_y_xzzz_yyz, to_x_y_xzzz_yzz, to_x_y_xzzz_zzz, to_xxzzz_xx, to_xxzzz_xxxy, to_xxzzz_xxyy, to_xxzzz_xxyz, to_xxzzz_xy, to_xxzzz_xyyy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_yy, to_xxzzz_yyyy, to_xxzzz_yyyz, to_xxzzz_yyzz, to_xxzzz_yz, to_xxzzz_yzzz, to_xxzzz_zz, to_zzz_xx, to_zzz_xxxy, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_yy, to_zzz_yyyy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzzz_xxx[k] = -2.0 * to_zzz_xxxy[k] * tke_0 + 4.0 * to_xxzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxy[k] = to_zzz_xx[k] - 2.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xxz[k] = -2.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xyy[k] = 2.0 * to_zzz_xy[k] - 2.0 * to_zzz_xyyy[k] * tke_0 - 4.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xyz[k] = to_zzz_xz[k] - 2.0 * to_zzz_xyyz[k] * tke_0 - 2.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xzz[k] = -2.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yyy[k] = 3.0 * to_zzz_yy[k] - 2.0 * to_zzz_yyyy[k] * tke_0 - 6.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yyz[k] = 2.0 * to_zzz_yz[k] - 2.0 * to_zzz_yyyz[k] * tke_0 - 4.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yzz[k] = to_zzz_zz[k] - 2.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_zzz[k] = -2.0 * to_zzz_yzzz[k] * tke_0 + 4.0 * to_xxzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 250-260 components of targeted buffer : GF

        auto to_x_y_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 100);

        auto to_x_y_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 101);

        auto to_x_y_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 102);

        auto to_x_y_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 103);

        auto to_x_y_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 104);

        auto to_x_y_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 105);

        auto to_x_y_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 106);

        auto to_x_y_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 107);

        auto to_x_y_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 108);

        auto to_x_y_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_x_y_yyyy_xxx, to_x_y_yyyy_xxy, to_x_y_yyyy_xxz, to_x_y_yyyy_xyy, to_x_y_yyyy_xyz, to_x_y_yyyy_xzz, to_x_y_yyyy_yyy, to_x_y_yyyy_yyz, to_x_y_yyyy_yzz, to_x_y_yyyy_zzz, to_xyyyy_xx, to_xyyyy_xxxy, to_xyyyy_xxyy, to_xyyyy_xxyz, to_xyyyy_xy, to_xyyyy_xyyy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_yy, to_xyyyy_yyyy, to_xyyyy_yyyz, to_xyyyy_yyzz, to_xyyyy_yz, to_xyyyy_yzzz, to_xyyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyy_xxx[k] = 4.0 * to_xyyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxy[k] = -2.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xxz[k] = 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xyy[k] = -4.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xyz[k] = -2.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xzz[k] = 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yyy[k] = -6.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yyz[k] = -4.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yzz[k] = -2.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_zzz[k] = 4.0 * to_xyyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 260-270 components of targeted buffer : GF

        auto to_x_y_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 110);

        auto to_x_y_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 111);

        auto to_x_y_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 112);

        auto to_x_y_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 113);

        auto to_x_y_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 114);

        auto to_x_y_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 115);

        auto to_x_y_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 116);

        auto to_x_y_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 117);

        auto to_x_y_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 118);

        auto to_x_y_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_y_yyyz_xxx, to_x_y_yyyz_xxy, to_x_y_yyyz_xxz, to_x_y_yyyz_xyy, to_x_y_yyyz_xyz, to_x_y_yyyz_xzz, to_x_y_yyyz_yyy, to_x_y_yyyz_yyz, to_x_y_yyyz_yzz, to_x_y_yyyz_zzz, to_xyyyz_xx, to_xyyyz_xxxy, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_yy, to_xyyyz_yyyy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyz_xxx[k] = 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxy[k] = -2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xxz[k] = 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xyy[k] = -4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xyz[k] = -2.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xzz[k] = 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yyy[k] = -6.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yyz[k] = -4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yzz[k] = -2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_zzz[k] = 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-280 components of targeted buffer : GF

        auto to_x_y_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 120);

        auto to_x_y_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 121);

        auto to_x_y_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 122);

        auto to_x_y_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 123);

        auto to_x_y_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 124);

        auto to_x_y_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 125);

        auto to_x_y_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 126);

        auto to_x_y_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 127);

        auto to_x_y_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 128);

        auto to_x_y_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_x_y_yyzz_xxx, to_x_y_yyzz_xxy, to_x_y_yyzz_xxz, to_x_y_yyzz_xyy, to_x_y_yyzz_xyz, to_x_y_yyzz_xzz, to_x_y_yyzz_yyy, to_x_y_yyzz_yyz, to_x_y_yyzz_yzz, to_x_y_yyzz_zzz, to_xyyzz_xx, to_xyyzz_xxxy, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_yy, to_xyyzz_yyyy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyzz_xxx[k] = 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxy[k] = -2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xxz[k] = 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xyy[k] = -4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xyz[k] = -2.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xzz[k] = 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yyy[k] = -6.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yyz[k] = -4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yzz[k] = -2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_zzz[k] = 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 280-290 components of targeted buffer : GF

        auto to_x_y_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 130);

        auto to_x_y_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 131);

        auto to_x_y_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 132);

        auto to_x_y_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 133);

        auto to_x_y_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 134);

        auto to_x_y_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 135);

        auto to_x_y_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 136);

        auto to_x_y_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 137);

        auto to_x_y_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 138);

        auto to_x_y_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_x_y_yzzz_xxx, to_x_y_yzzz_xxy, to_x_y_yzzz_xxz, to_x_y_yzzz_xyy, to_x_y_yzzz_xyz, to_x_y_yzzz_xzz, to_x_y_yzzz_yyy, to_x_y_yzzz_yyz, to_x_y_yzzz_yzz, to_x_y_yzzz_zzz, to_xyzzz_xx, to_xyzzz_xxxy, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_yy, to_xyzzz_yyyy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzzz_xxx[k] = 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxy[k] = -2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xxz[k] = 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xyy[k] = -4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xyz[k] = -2.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xzz[k] = 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yyy[k] = -6.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yyz[k] = -4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yzz[k] = -2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_zzz[k] = 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 290-300 components of targeted buffer : GF

        auto to_x_y_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 140);

        auto to_x_y_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 141);

        auto to_x_y_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 142);

        auto to_x_y_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 143);

        auto to_x_y_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 144);

        auto to_x_y_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 145);

        auto to_x_y_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 146);

        auto to_x_y_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 147);

        auto to_x_y_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 148);

        auto to_x_y_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 1 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_y_zzzz_xxx, to_x_y_zzzz_xxy, to_x_y_zzzz_xxz, to_x_y_zzzz_xyy, to_x_y_zzzz_xyz, to_x_y_zzzz_xzz, to_x_y_zzzz_yyy, to_x_y_zzzz_yyz, to_x_y_zzzz_yzz, to_x_y_zzzz_zzz, to_xzzzz_xx, to_xzzzz_xxxy, to_xzzzz_xxyy, to_xzzzz_xxyz, to_xzzzz_xy, to_xzzzz_xyyy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_yy, to_xzzzz_yyyy, to_xzzzz_yyyz, to_xzzzz_yyzz, to_xzzzz_yz, to_xzzzz_yzzz, to_xzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzzz_xxx[k] = 4.0 * to_xzzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxy[k] = -2.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xxz[k] = 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xyy[k] = -4.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xyz[k] = -2.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xzz[k] = 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yyy[k] = -6.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yyz[k] = -4.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yzz[k] = -2.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_zzz[k] = 4.0 * to_xzzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-310 components of targeted buffer : GF

        auto to_x_z_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 0);

        auto to_x_z_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 1);

        auto to_x_z_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 2);

        auto to_x_z_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 3);

        auto to_x_z_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 4);

        auto to_x_z_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 5);

        auto to_x_z_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 6);

        auto to_x_z_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 7);

        auto to_x_z_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 8);

        auto to_x_z_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_x_z_xxxx_xxx, to_x_z_xxxx_xxy, to_x_z_xxxx_xxz, to_x_z_xxxx_xyy, to_x_z_xxxx_xyz, to_x_z_xxxx_xzz, to_x_z_xxxx_yyy, to_x_z_xxxx_yyz, to_x_z_xxxx_yzz, to_x_z_xxxx_zzz, to_xxx_xx, to_xxx_xxxz, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxx_zzzz, to_xxxxx_xx, to_xxxxx_xxxz, to_xxxxx_xxyz, to_xxxxx_xxzz, to_xxxxx_xy, to_xxxxx_xyyz, to_xxxxx_xyzz, to_xxxxx_xz, to_xxxxx_xzzz, to_xxxxx_yy, to_xxxxx_yyyz, to_xxxxx_yyzz, to_xxxxx_yz, to_xxxxx_yzzz, to_xxxxx_zz, to_xxxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxx_xxx[k] = -8.0 * to_xxx_xxxz[k] * tke_0 + 4.0 * to_xxxxx_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxy[k] = -8.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xxz[k] = 4.0 * to_xxx_xx[k] - 8.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxxx_xx[k] * tbe_0 + 4.0 * to_xxxxx_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xyy[k] = -8.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xyz[k] = 4.0 * to_xxx_xy[k] - 8.0 * to_xxx_xyzz[k] * tke_0 - 2.0 * to_xxxxx_xy[k] * tbe_0 + 4.0 * to_xxxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xzz[k] = 8.0 * to_xxx_xz[k] - 8.0 * to_xxx_xzzz[k] * tke_0 - 4.0 * to_xxxxx_xz[k] * tbe_0 + 4.0 * to_xxxxx_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yyy[k] = -8.0 * to_xxx_yyyz[k] * tke_0 + 4.0 * to_xxxxx_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yyz[k] = 4.0 * to_xxx_yy[k] - 8.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxxx_yy[k] * tbe_0 + 4.0 * to_xxxxx_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yzz[k] = 8.0 * to_xxx_yz[k] - 8.0 * to_xxx_yzzz[k] * tke_0 - 4.0 * to_xxxxx_yz[k] * tbe_0 + 4.0 * to_xxxxx_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_zzz[k] = 12.0 * to_xxx_zz[k] - 8.0 * to_xxx_zzzz[k] * tke_0 - 6.0 * to_xxxxx_zz[k] * tbe_0 + 4.0 * to_xxxxx_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 310-320 components of targeted buffer : GF

        auto to_x_z_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 10);

        auto to_x_z_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 11);

        auto to_x_z_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 12);

        auto to_x_z_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 13);

        auto to_x_z_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 14);

        auto to_x_z_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 15);

        auto to_x_z_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 16);

        auto to_x_z_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 17);

        auto to_x_z_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 18);

        auto to_x_z_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_x_z_xxxy_xxx, to_x_z_xxxy_xxy, to_x_z_xxxy_xxz, to_x_z_xxxy_xyy, to_x_z_xxxy_xyz, to_x_z_xxxy_xzz, to_x_z_xxxy_yyy, to_x_z_xxxy_yyz, to_x_z_xxxy_yzz, to_x_z_xxxy_zzz, to_xxxxy_xx, to_xxxxy_xxxz, to_xxxxy_xxyz, to_xxxxy_xxzz, to_xxxxy_xy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_xzzz, to_xxxxy_yy, to_xxxxy_yyyz, to_xxxxy_yyzz, to_xxxxy_yz, to_xxxxy_yzzz, to_xxxxy_zz, to_xxxxy_zzzz, to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxy_xxx[k] = -6.0 * to_xxy_xxxz[k] * tke_0 + 4.0 * to_xxxxy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxy[k] = -6.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xxz[k] = 3.0 * to_xxy_xx[k] - 6.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xyy[k] = -6.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xyz[k] = 3.0 * to_xxy_xy[k] - 6.0 * to_xxy_xyzz[k] * tke_0 - 2.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xzz[k] = 6.0 * to_xxy_xz[k] - 6.0 * to_xxy_xzzz[k] * tke_0 - 4.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yyy[k] = -6.0 * to_xxy_yyyz[k] * tke_0 + 4.0 * to_xxxxy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yyz[k] = 3.0 * to_xxy_yy[k] - 6.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yzz[k] = 6.0 * to_xxy_yz[k] - 6.0 * to_xxy_yzzz[k] * tke_0 - 4.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_zzz[k] = 9.0 * to_xxy_zz[k] - 6.0 * to_xxy_zzzz[k] * tke_0 - 6.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 320-330 components of targeted buffer : GF

        auto to_x_z_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 20);

        auto to_x_z_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 21);

        auto to_x_z_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 22);

        auto to_x_z_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 23);

        auto to_x_z_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 24);

        auto to_x_z_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 25);

        auto to_x_z_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 26);

        auto to_x_z_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 27);

        auto to_x_z_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 28);

        auto to_x_z_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_z_xxxz_xxx, to_x_z_xxxz_xxy, to_x_z_xxxz_xxz, to_x_z_xxxz_xyy, to_x_z_xxxz_xyz, to_x_z_xxxz_xzz, to_x_z_xxxz_yyy, to_x_z_xxxz_yyz, to_x_z_xxxz_yzz, to_x_z_xxxz_zzz, to_xxxxz_xx, to_xxxxz_xxxz, to_xxxxz_xxyz, to_xxxxz_xxzz, to_xxxxz_xy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_xzzz, to_xxxxz_yy, to_xxxxz_yyyz, to_xxxxz_yyzz, to_xxxxz_yz, to_xxxxz_yzzz, to_xxxxz_zz, to_xxxxz_zzzz, to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxz_xxx[k] = -6.0 * to_xxz_xxxz[k] * tke_0 + 4.0 * to_xxxxz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxy[k] = -6.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xxz[k] = 3.0 * to_xxz_xx[k] - 6.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xyy[k] = -6.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xyz[k] = 3.0 * to_xxz_xy[k] - 6.0 * to_xxz_xyzz[k] * tke_0 - 2.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xzz[k] = 6.0 * to_xxz_xz[k] - 6.0 * to_xxz_xzzz[k] * tke_0 - 4.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yyy[k] = -6.0 * to_xxz_yyyz[k] * tke_0 + 4.0 * to_xxxxz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yyz[k] = 3.0 * to_xxz_yy[k] - 6.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yzz[k] = 6.0 * to_xxz_yz[k] - 6.0 * to_xxz_yzzz[k] * tke_0 - 4.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_zzz[k] = 9.0 * to_xxz_zz[k] - 6.0 * to_xxz_zzzz[k] * tke_0 - 6.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-340 components of targeted buffer : GF

        auto to_x_z_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 30);

        auto to_x_z_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 31);

        auto to_x_z_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 32);

        auto to_x_z_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 33);

        auto to_x_z_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 34);

        auto to_x_z_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 35);

        auto to_x_z_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 36);

        auto to_x_z_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 37);

        auto to_x_z_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 38);

        auto to_x_z_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_x_z_xxyy_xxx, to_x_z_xxyy_xxy, to_x_z_xxyy_xxz, to_x_z_xxyy_xyy, to_x_z_xxyy_xyz, to_x_z_xxyy_xzz, to_x_z_xxyy_yyy, to_x_z_xxyy_yyz, to_x_z_xxyy_yzz, to_x_z_xxyy_zzz, to_xxxyy_xx, to_xxxyy_xxxz, to_xxxyy_xxyz, to_xxxyy_xxzz, to_xxxyy_xy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_xzzz, to_xxxyy_yy, to_xxxyy_yyyz, to_xxxyy_yyzz, to_xxxyy_yz, to_xxxyy_yzzz, to_xxxyy_zz, to_xxxyy_zzzz, to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyy_xxx[k] = -4.0 * to_xyy_xxxz[k] * tke_0 + 4.0 * to_xxxyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxy[k] = -4.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xxz[k] = 2.0 * to_xyy_xx[k] - 4.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xyy[k] = -4.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xyz[k] = 2.0 * to_xyy_xy[k] - 4.0 * to_xyy_xyzz[k] * tke_0 - 2.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xzz[k] = 4.0 * to_xyy_xz[k] - 4.0 * to_xyy_xzzz[k] * tke_0 - 4.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yyy[k] = -4.0 * to_xyy_yyyz[k] * tke_0 + 4.0 * to_xxxyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yyz[k] = 2.0 * to_xyy_yy[k] - 4.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yzz[k] = 4.0 * to_xyy_yz[k] - 4.0 * to_xyy_yzzz[k] * tke_0 - 4.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_zzz[k] = 6.0 * to_xyy_zz[k] - 4.0 * to_xyy_zzzz[k] * tke_0 - 6.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 340-350 components of targeted buffer : GF

        auto to_x_z_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 40);

        auto to_x_z_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 41);

        auto to_x_z_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 42);

        auto to_x_z_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 43);

        auto to_x_z_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 44);

        auto to_x_z_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 45);

        auto to_x_z_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 46);

        auto to_x_z_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 47);

        auto to_x_z_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 48);

        auto to_x_z_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_x_z_xxyz_xxx, to_x_z_xxyz_xxy, to_x_z_xxyz_xxz, to_x_z_xxyz_xyy, to_x_z_xxyz_xyz, to_x_z_xxyz_xzz, to_x_z_xxyz_yyy, to_x_z_xxyz_yyz, to_x_z_xxyz_yzz, to_x_z_xxyz_zzz, to_xxxyz_xx, to_xxxyz_xxxz, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_xxxyz_zzzz, to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyz_xxx[k] = -4.0 * to_xyz_xxxz[k] * tke_0 + 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxy[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xxz[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xyy[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xyz[k] = 2.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyzz[k] * tke_0 - 2.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xzz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xzzz[k] * tke_0 - 4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yyy[k] = -4.0 * to_xyz_yyyz[k] * tke_0 + 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yyz[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yzz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yzzz[k] * tke_0 - 4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_zzz[k] = 6.0 * to_xyz_zz[k] - 4.0 * to_xyz_zzzz[k] * tke_0 - 6.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 350-360 components of targeted buffer : GF

        auto to_x_z_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 50);

        auto to_x_z_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 51);

        auto to_x_z_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 52);

        auto to_x_z_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 53);

        auto to_x_z_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 54);

        auto to_x_z_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 55);

        auto to_x_z_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 56);

        auto to_x_z_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 57);

        auto to_x_z_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 58);

        auto to_x_z_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_z_xxzz_xxx, to_x_z_xxzz_xxy, to_x_z_xxzz_xxz, to_x_z_xxzz_xyy, to_x_z_xxzz_xyz, to_x_z_xxzz_xzz, to_x_z_xxzz_yyy, to_x_z_xxzz_yyz, to_x_z_xxzz_yzz, to_x_z_xxzz_zzz, to_xxxzz_xx, to_xxxzz_xxxz, to_xxxzz_xxyz, to_xxxzz_xxzz, to_xxxzz_xy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_xzzz, to_xxxzz_yy, to_xxxzz_yyyz, to_xxxzz_yyzz, to_xxxzz_yz, to_xxxzz_yzzz, to_xxxzz_zz, to_xxxzz_zzzz, to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxzz_xxx[k] = -4.0 * to_xzz_xxxz[k] * tke_0 + 4.0 * to_xxxzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxy[k] = -4.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xxz[k] = 2.0 * to_xzz_xx[k] - 4.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xyy[k] = -4.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xyz[k] = 2.0 * to_xzz_xy[k] - 4.0 * to_xzz_xyzz[k] * tke_0 - 2.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xzz[k] = 4.0 * to_xzz_xz[k] - 4.0 * to_xzz_xzzz[k] * tke_0 - 4.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yyy[k] = -4.0 * to_xzz_yyyz[k] * tke_0 + 4.0 * to_xxxzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yyz[k] = 2.0 * to_xzz_yy[k] - 4.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yzz[k] = 4.0 * to_xzz_yz[k] - 4.0 * to_xzz_yzzz[k] * tke_0 - 4.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_zzz[k] = 6.0 * to_xzz_zz[k] - 4.0 * to_xzz_zzzz[k] * tke_0 - 6.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-370 components of targeted buffer : GF

        auto to_x_z_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 60);

        auto to_x_z_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 61);

        auto to_x_z_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 62);

        auto to_x_z_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 63);

        auto to_x_z_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 64);

        auto to_x_z_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 65);

        auto to_x_z_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 66);

        auto to_x_z_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 67);

        auto to_x_z_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 68);

        auto to_x_z_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_x_z_xyyy_xxx, to_x_z_xyyy_xxy, to_x_z_xyyy_xxz, to_x_z_xyyy_xyy, to_x_z_xyyy_xyz, to_x_z_xyyy_xzz, to_x_z_xyyy_yyy, to_x_z_xyyy_yyz, to_x_z_xyyy_yzz, to_x_z_xyyy_zzz, to_xxyyy_xx, to_xxyyy_xxxz, to_xxyyy_xxyz, to_xxyyy_xxzz, to_xxyyy_xy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_xzzz, to_xxyyy_yy, to_xxyyy_yyyz, to_xxyyy_yyzz, to_xxyyy_yz, to_xxyyy_yzzz, to_xxyyy_zz, to_xxyyy_zzzz, to_yyy_xx, to_yyy_xxxz, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyy_xxx[k] = -2.0 * to_yyy_xxxz[k] * tke_0 + 4.0 * to_xxyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxy[k] = -2.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xxz[k] = to_yyy_xx[k] - 2.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xyy[k] = -2.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xyz[k] = to_yyy_xy[k] - 2.0 * to_yyy_xyzz[k] * tke_0 - 2.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xzz[k] = 2.0 * to_yyy_xz[k] - 2.0 * to_yyy_xzzz[k] * tke_0 - 4.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yyy[k] = -2.0 * to_yyy_yyyz[k] * tke_0 + 4.0 * to_xxyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yyz[k] = to_yyy_yy[k] - 2.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yzz[k] = 2.0 * to_yyy_yz[k] - 2.0 * to_yyy_yzzz[k] * tke_0 - 4.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_zzz[k] = 3.0 * to_yyy_zz[k] - 2.0 * to_yyy_zzzz[k] * tke_0 - 6.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 370-380 components of targeted buffer : GF

        auto to_x_z_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 70);

        auto to_x_z_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 71);

        auto to_x_z_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 72);

        auto to_x_z_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 73);

        auto to_x_z_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 74);

        auto to_x_z_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 75);

        auto to_x_z_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 76);

        auto to_x_z_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 77);

        auto to_x_z_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 78);

        auto to_x_z_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_x_z_xyyz_xxx, to_x_z_xyyz_xxy, to_x_z_xyyz_xxz, to_x_z_xyyz_xyy, to_x_z_xyyz_xyz, to_x_z_xyyz_xzz, to_x_z_xyyz_yyy, to_x_z_xyyz_yyz, to_x_z_xyyz_yzz, to_x_z_xyyz_zzz, to_xxyyz_xx, to_xxyyz_xxxz, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_xxyyz_zzzz, to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyz_xxx[k] = -2.0 * to_yyz_xxxz[k] * tke_0 + 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxy[k] = -2.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xxz[k] = to_yyz_xx[k] - 2.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xyy[k] = -2.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xyz[k] = to_yyz_xy[k] - 2.0 * to_yyz_xyzz[k] * tke_0 - 2.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xzz[k] = 2.0 * to_yyz_xz[k] - 2.0 * to_yyz_xzzz[k] * tke_0 - 4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yyy[k] = -2.0 * to_yyz_yyyz[k] * tke_0 + 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yyz[k] = to_yyz_yy[k] - 2.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yzz[k] = 2.0 * to_yyz_yz[k] - 2.0 * to_yyz_yzzz[k] * tke_0 - 4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_zzz[k] = 3.0 * to_yyz_zz[k] - 2.0 * to_yyz_zzzz[k] * tke_0 - 6.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 380-390 components of targeted buffer : GF

        auto to_x_z_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 80);

        auto to_x_z_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 81);

        auto to_x_z_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 82);

        auto to_x_z_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 83);

        auto to_x_z_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 84);

        auto to_x_z_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 85);

        auto to_x_z_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 86);

        auto to_x_z_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 87);

        auto to_x_z_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 88);

        auto to_x_z_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_z_xyzz_xxx, to_x_z_xyzz_xxy, to_x_z_xyzz_xxz, to_x_z_xyzz_xyy, to_x_z_xyzz_xyz, to_x_z_xyzz_xzz, to_x_z_xyzz_yyy, to_x_z_xyzz_yyz, to_x_z_xyzz_yzz, to_x_z_xyzz_zzz, to_xxyzz_xx, to_xxyzz_xxxz, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_xxyzz_zzzz, to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyzz_xxx[k] = -2.0 * to_yzz_xxxz[k] * tke_0 + 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxy[k] = -2.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xxz[k] = to_yzz_xx[k] - 2.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xyy[k] = -2.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xyz[k] = to_yzz_xy[k] - 2.0 * to_yzz_xyzz[k] * tke_0 - 2.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xzz[k] = 2.0 * to_yzz_xz[k] - 2.0 * to_yzz_xzzz[k] * tke_0 - 4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yyy[k] = -2.0 * to_yzz_yyyz[k] * tke_0 + 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yyz[k] = to_yzz_yy[k] - 2.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yzz[k] = 2.0 * to_yzz_yz[k] - 2.0 * to_yzz_yzzz[k] * tke_0 - 4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_zzz[k] = 3.0 * to_yzz_zz[k] - 2.0 * to_yzz_zzzz[k] * tke_0 - 6.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-400 components of targeted buffer : GF

        auto to_x_z_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 90);

        auto to_x_z_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 91);

        auto to_x_z_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 92);

        auto to_x_z_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 93);

        auto to_x_z_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 94);

        auto to_x_z_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 95);

        auto to_x_z_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 96);

        auto to_x_z_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 97);

        auto to_x_z_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 98);

        auto to_x_z_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_x_z_xzzz_xxx, to_x_z_xzzz_xxy, to_x_z_xzzz_xxz, to_x_z_xzzz_xyy, to_x_z_xzzz_xyz, to_x_z_xzzz_xzz, to_x_z_xzzz_yyy, to_x_z_xzzz_yyz, to_x_z_xzzz_yzz, to_x_z_xzzz_zzz, to_xxzzz_xx, to_xxzzz_xxxz, to_xxzzz_xxyz, to_xxzzz_xxzz, to_xxzzz_xy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_xzzz, to_xxzzz_yy, to_xxzzz_yyyz, to_xxzzz_yyzz, to_xxzzz_yz, to_xxzzz_yzzz, to_xxzzz_zz, to_xxzzz_zzzz, to_zzz_xx, to_zzz_xxxz, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzzz_xxx[k] = -2.0 * to_zzz_xxxz[k] * tke_0 + 4.0 * to_xxzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxy[k] = -2.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xxz[k] = to_zzz_xx[k] - 2.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xyy[k] = -2.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xyz[k] = to_zzz_xy[k] - 2.0 * to_zzz_xyzz[k] * tke_0 - 2.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xzz[k] = 2.0 * to_zzz_xz[k] - 2.0 * to_zzz_xzzz[k] * tke_0 - 4.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yyy[k] = -2.0 * to_zzz_yyyz[k] * tke_0 + 4.0 * to_xxzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yyz[k] = to_zzz_yy[k] - 2.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yzz[k] = 2.0 * to_zzz_yz[k] - 2.0 * to_zzz_yzzz[k] * tke_0 - 4.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_zzz[k] = 3.0 * to_zzz_zz[k] - 2.0 * to_zzz_zzzz[k] * tke_0 - 6.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 400-410 components of targeted buffer : GF

        auto to_x_z_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 100);

        auto to_x_z_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 101);

        auto to_x_z_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 102);

        auto to_x_z_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 103);

        auto to_x_z_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 104);

        auto to_x_z_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 105);

        auto to_x_z_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 106);

        auto to_x_z_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 107);

        auto to_x_z_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 108);

        auto to_x_z_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_x_z_yyyy_xxx, to_x_z_yyyy_xxy, to_x_z_yyyy_xxz, to_x_z_yyyy_xyy, to_x_z_yyyy_xyz, to_x_z_yyyy_xzz, to_x_z_yyyy_yyy, to_x_z_yyyy_yyz, to_x_z_yyyy_yzz, to_x_z_yyyy_zzz, to_xyyyy_xx, to_xyyyy_xxxz, to_xyyyy_xxyz, to_xyyyy_xxzz, to_xyyyy_xy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_xzzz, to_xyyyy_yy, to_xyyyy_yyyz, to_xyyyy_yyzz, to_xyyyy_yz, to_xyyyy_yzzz, to_xyyyy_zz, to_xyyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyy_xxx[k] = 4.0 * to_xyyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxy[k] = 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xxz[k] = -2.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xyy[k] = 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xyz[k] = -2.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xzz[k] = -4.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yyy[k] = 4.0 * to_xyyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yyz[k] = -2.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yzz[k] = -4.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_zzz[k] = -6.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 410-420 components of targeted buffer : GF

        auto to_x_z_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 110);

        auto to_x_z_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 111);

        auto to_x_z_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 112);

        auto to_x_z_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 113);

        auto to_x_z_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 114);

        auto to_x_z_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 115);

        auto to_x_z_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 116);

        auto to_x_z_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 117);

        auto to_x_z_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 118);

        auto to_x_z_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_z_yyyz_xxx, to_x_z_yyyz_xxy, to_x_z_yyyz_xxz, to_x_z_yyyz_xyy, to_x_z_yyyz_xyz, to_x_z_yyyz_xzz, to_x_z_yyyz_yyy, to_x_z_yyyz_yyz, to_x_z_yyyz_yzz, to_x_z_yyyz_zzz, to_xyyyz_xx, to_xyyyz_xxxz, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, to_xyyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyz_xxx[k] = 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxy[k] = 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xxz[k] = -2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xyy[k] = 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xyz[k] = -2.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xzz[k] = -4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yyy[k] = 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yyz[k] = -2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yzz[k] = -4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_zzz[k] = -6.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-430 components of targeted buffer : GF

        auto to_x_z_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 120);

        auto to_x_z_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 121);

        auto to_x_z_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 122);

        auto to_x_z_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 123);

        auto to_x_z_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 124);

        auto to_x_z_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 125);

        auto to_x_z_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 126);

        auto to_x_z_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 127);

        auto to_x_z_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 128);

        auto to_x_z_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_x_z_yyzz_xxx, to_x_z_yyzz_xxy, to_x_z_yyzz_xxz, to_x_z_yyzz_xyy, to_x_z_yyzz_xyz, to_x_z_yyzz_xzz, to_x_z_yyzz_yyy, to_x_z_yyzz_yyz, to_x_z_yyzz_yzz, to_x_z_yyzz_zzz, to_xyyzz_xx, to_xyyzz_xxxz, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, to_xyyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyzz_xxx[k] = 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxy[k] = 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xxz[k] = -2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xyy[k] = 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xyz[k] = -2.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xzz[k] = -4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yyy[k] = 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yyz[k] = -2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yzz[k] = -4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_zzz[k] = -6.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 430-440 components of targeted buffer : GF

        auto to_x_z_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 130);

        auto to_x_z_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 131);

        auto to_x_z_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 132);

        auto to_x_z_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 133);

        auto to_x_z_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 134);

        auto to_x_z_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 135);

        auto to_x_z_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 136);

        auto to_x_z_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 137);

        auto to_x_z_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 138);

        auto to_x_z_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_x_z_yzzz_xxx, to_x_z_yzzz_xxy, to_x_z_yzzz_xxz, to_x_z_yzzz_xyy, to_x_z_yzzz_xyz, to_x_z_yzzz_xzz, to_x_z_yzzz_yyy, to_x_z_yzzz_yyz, to_x_z_yzzz_yzz, to_x_z_yzzz_zzz, to_xyzzz_xx, to_xyzzz_xxxz, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, to_xyzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzzz_xxx[k] = 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxy[k] = 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xxz[k] = -2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xyy[k] = 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xyz[k] = -2.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xzz[k] = -4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yyy[k] = 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yyz[k] = -2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yzz[k] = -4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_zzz[k] = -6.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 440-450 components of targeted buffer : GF

        auto to_x_z_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 140);

        auto to_x_z_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 141);

        auto to_x_z_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 142);

        auto to_x_z_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 143);

        auto to_x_z_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 144);

        auto to_x_z_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 145);

        auto to_x_z_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 146);

        auto to_x_z_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 147);

        auto to_x_z_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 148);

        auto to_x_z_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 2 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_z_zzzz_xxx, to_x_z_zzzz_xxy, to_x_z_zzzz_xxz, to_x_z_zzzz_xyy, to_x_z_zzzz_xyz, to_x_z_zzzz_xzz, to_x_z_zzzz_yyy, to_x_z_zzzz_yyz, to_x_z_zzzz_yzz, to_x_z_zzzz_zzz, to_xzzzz_xx, to_xzzzz_xxxz, to_xzzzz_xxyz, to_xzzzz_xxzz, to_xzzzz_xy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_xzzz, to_xzzzz_yy, to_xzzzz_yyyz, to_xzzzz_yyzz, to_xzzzz_yz, to_xzzzz_yzzz, to_xzzzz_zz, to_xzzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzzz_xxx[k] = 4.0 * to_xzzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxy[k] = 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xxz[k] = -2.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xyy[k] = 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xyz[k] = -2.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xzz[k] = -4.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yyy[k] = 4.0 * to_xzzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yyz[k] = -2.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yzz[k] = -4.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_zzz[k] = -6.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-460 components of targeted buffer : GF

        auto to_y_x_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 0);

        auto to_y_x_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 1);

        auto to_y_x_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 2);

        auto to_y_x_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 3);

        auto to_y_x_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 4);

        auto to_y_x_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 5);

        auto to_y_x_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 6);

        auto to_y_x_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 7);

        auto to_y_x_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 8);

        auto to_y_x_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxy_xx, to_xxxxy_xxxx, to_xxxxy_xxxy, to_xxxxy_xxxz, to_xxxxy_xxyy, to_xxxxy_xxyz, to_xxxxy_xxzz, to_xxxxy_xy, to_xxxxy_xyyy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_xzzz, to_xxxxy_yy, to_xxxxy_yz, to_xxxxy_zz, to_y_x_xxxx_xxx, to_y_x_xxxx_xxy, to_y_x_xxxx_xxz, to_y_x_xxxx_xyy, to_y_x_xxxx_xyz, to_y_x_xxxx_xzz, to_y_x_xxxx_yyy, to_y_x_xxxx_yyz, to_y_x_xxxx_yzz, to_y_x_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxx_xxx[k] = -6.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxy[k] = -4.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xxz[k] = -4.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xyy[k] = -2.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xyz[k] = -2.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xzz[k] = -2.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yyy[k] = 4.0 * to_xxxxy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yyz[k] = 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yzz[k] = 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_zzz[k] = 4.0 * to_xxxxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 460-470 components of targeted buffer : GF

        auto to_y_x_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 10);

        auto to_y_x_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 11);

        auto to_y_x_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 12);

        auto to_y_x_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 13);

        auto to_y_x_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 14);

        auto to_y_x_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 15);

        auto to_y_x_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 16);

        auto to_y_x_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 17);

        auto to_y_x_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 18);

        auto to_y_x_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxx, to_xxx_xxxy, to_xxx_xxxz, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yz, to_xxx_zz, to_xxxyy_xx, to_xxxyy_xxxx, to_xxxyy_xxxy, to_xxxyy_xxxz, to_xxxyy_xxyy, to_xxxyy_xxyz, to_xxxyy_xxzz, to_xxxyy_xy, to_xxxyy_xyyy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_xzzz, to_xxxyy_yy, to_xxxyy_yz, to_xxxyy_zz, to_y_x_xxxy_xxx, to_y_x_xxxy_xxy, to_y_x_xxxy_xxz, to_y_x_xxxy_xyy, to_y_x_xxxy_xyz, to_y_x_xxxy_xzz, to_y_x_xxxy_yyy, to_y_x_xxxy_yyz, to_y_x_xxxy_yzz, to_y_x_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxy_xxx[k] = 3.0 * to_xxx_xx[k] - 2.0 * to_xxx_xxxx[k] * tke_0 - 6.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxy[k] = 2.0 * to_xxx_xy[k] - 2.0 * to_xxx_xxxy[k] * tke_0 - 4.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xxz[k] = 2.0 * to_xxx_xz[k] - 2.0 * to_xxx_xxxz[k] * tke_0 - 4.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xyy[k] = to_xxx_yy[k] - 2.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xyz[k] = to_xxx_yz[k] - 2.0 * to_xxx_xxyz[k] * tke_0 - 2.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xzz[k] = to_xxx_zz[k] - 2.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yyy[k] = -2.0 * to_xxx_xyyy[k] * tke_0 + 4.0 * to_xxxyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yyz[k] = -2.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yzz[k] = -2.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_zzz[k] = -2.0 * to_xxx_xzzz[k] * tke_0 + 4.0 * to_xxxyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 470-480 components of targeted buffer : GF

        auto to_y_x_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 20);

        auto to_y_x_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 21);

        auto to_y_x_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 22);

        auto to_y_x_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 23);

        auto to_y_x_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 24);

        auto to_y_x_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 25);

        auto to_y_x_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 26);

        auto to_y_x_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 27);

        auto to_y_x_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 28);

        auto to_y_x_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxx, to_xxxyz_xxxy, to_xxxyz_xxxz, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yz, to_xxxyz_zz, to_y_x_xxxz_xxx, to_y_x_xxxz_xxy, to_y_x_xxxz_xxz, to_y_x_xxxz_xyy, to_y_x_xxxz_xyz, to_y_x_xxxz_xzz, to_y_x_xxxz_yyy, to_y_x_xxxz_yyz, to_y_x_xxxz_yzz, to_y_x_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxz_xxx[k] = -6.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxy[k] = -4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xxz[k] = -4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xyy[k] = -2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xyz[k] = -2.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xzz[k] = -2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yyy[k] = 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yyz[k] = 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yzz[k] = 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_zzz[k] = 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-490 components of targeted buffer : GF

        auto to_y_x_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 30);

        auto to_y_x_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 31);

        auto to_y_x_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 32);

        auto to_y_x_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 33);

        auto to_y_x_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 34);

        auto to_y_x_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 35);

        auto to_y_x_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 36);

        auto to_y_x_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 37);

        auto to_y_x_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 38);

        auto to_y_x_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, to_xxyyy_xx, to_xxyyy_xxxx, to_xxyyy_xxxy, to_xxyyy_xxxz, to_xxyyy_xxyy, to_xxyyy_xxyz, to_xxyyy_xxzz, to_xxyyy_xy, to_xxyyy_xyyy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_xzzz, to_xxyyy_yy, to_xxyyy_yz, to_xxyyy_zz, to_y_x_xxyy_xxx, to_y_x_xxyy_xxy, to_y_x_xxyy_xxz, to_y_x_xxyy_xyy, to_y_x_xxyy_xyz, to_y_x_xxyy_xzz, to_y_x_xxyy_yyy, to_y_x_xxyy_yyz, to_y_x_xxyy_yzz, to_y_x_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyy_xxx[k] = 6.0 * to_xxy_xx[k] - 4.0 * to_xxy_xxxx[k] * tke_0 - 6.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxy[k] = 4.0 * to_xxy_xy[k] - 4.0 * to_xxy_xxxy[k] * tke_0 - 4.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xxz[k] = 4.0 * to_xxy_xz[k] - 4.0 * to_xxy_xxxz[k] * tke_0 - 4.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xyy[k] = 2.0 * to_xxy_yy[k] - 4.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xyz[k] = 2.0 * to_xxy_yz[k] - 4.0 * to_xxy_xxyz[k] * tke_0 - 2.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xzz[k] = 2.0 * to_xxy_zz[k] - 4.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yyy[k] = -4.0 * to_xxy_xyyy[k] * tke_0 + 4.0 * to_xxyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yyz[k] = -4.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yzz[k] = -4.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_zzz[k] = -4.0 * to_xxy_xzzz[k] * tke_0 + 4.0 * to_xxyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 490-500 components of targeted buffer : GF

        auto to_y_x_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 40);

        auto to_y_x_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 41);

        auto to_y_x_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 42);

        auto to_y_x_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 43);

        auto to_y_x_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 44);

        auto to_y_x_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 45);

        auto to_y_x_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 46);

        auto to_y_x_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 47);

        auto to_y_x_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 48);

        auto to_y_x_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxx, to_xxyyz_xxxy, to_xxyyz_xxxz, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yz, to_xxyyz_zz, to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, to_y_x_xxyz_xxx, to_y_x_xxyz_xxy, to_y_x_xxyz_xxz, to_y_x_xxyz_xyy, to_y_x_xxyz_xyz, to_y_x_xxyz_xzz, to_y_x_xxyz_yyy, to_y_x_xxyz_yyz, to_y_x_xxyz_yzz, to_y_x_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyz_xxx[k] = 3.0 * to_xxz_xx[k] - 2.0 * to_xxz_xxxx[k] * tke_0 - 6.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxy[k] = 2.0 * to_xxz_xy[k] - 2.0 * to_xxz_xxxy[k] * tke_0 - 4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xxz[k] = 2.0 * to_xxz_xz[k] - 2.0 * to_xxz_xxxz[k] * tke_0 - 4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xyy[k] = to_xxz_yy[k] - 2.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xyz[k] = to_xxz_yz[k] - 2.0 * to_xxz_xxyz[k] * tke_0 - 2.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xzz[k] = to_xxz_zz[k] - 2.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yyy[k] = -2.0 * to_xxz_xyyy[k] * tke_0 + 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yyz[k] = -2.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yzz[k] = -2.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_zzz[k] = -2.0 * to_xxz_xzzz[k] * tke_0 + 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 500-510 components of targeted buffer : GF

        auto to_y_x_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 50);

        auto to_y_x_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 51);

        auto to_y_x_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 52);

        auto to_y_x_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 53);

        auto to_y_x_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 54);

        auto to_y_x_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 55);

        auto to_y_x_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 56);

        auto to_y_x_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 57);

        auto to_y_x_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 58);

        auto to_y_x_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxyzz_xx, to_xxyzz_xxxx, to_xxyzz_xxxy, to_xxyzz_xxxz, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yz, to_xxyzz_zz, to_y_x_xxzz_xxx, to_y_x_xxzz_xxy, to_y_x_xxzz_xxz, to_y_x_xxzz_xyy, to_y_x_xxzz_xyz, to_y_x_xxzz_xzz, to_y_x_xxzz_yyy, to_y_x_xxzz_yyz, to_y_x_xxzz_yzz, to_y_x_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxzz_xxx[k] = -6.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxy[k] = -4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xxz[k] = -4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xyy[k] = -2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xyz[k] = -2.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xzz[k] = -2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yyy[k] = 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yyz[k] = 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yzz[k] = 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_zzz[k] = 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-520 components of targeted buffer : GF

        auto to_y_x_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 60);

        auto to_y_x_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 61);

        auto to_y_x_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 62);

        auto to_y_x_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 63);

        auto to_y_x_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 64);

        auto to_y_x_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 65);

        auto to_y_x_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 66);

        auto to_y_x_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 67);

        auto to_y_x_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 68);

        auto to_y_x_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, to_xyyyy_xx, to_xyyyy_xxxx, to_xyyyy_xxxy, to_xyyyy_xxxz, to_xyyyy_xxyy, to_xyyyy_xxyz, to_xyyyy_xxzz, to_xyyyy_xy, to_xyyyy_xyyy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_xzzz, to_xyyyy_yy, to_xyyyy_yz, to_xyyyy_zz, to_y_x_xyyy_xxx, to_y_x_xyyy_xxy, to_y_x_xyyy_xxz, to_y_x_xyyy_xyy, to_y_x_xyyy_xyz, to_y_x_xyyy_xzz, to_y_x_xyyy_yyy, to_y_x_xyyy_yyz, to_y_x_xyyy_yzz, to_y_x_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyy_xxx[k] = 9.0 * to_xyy_xx[k] - 6.0 * to_xyy_xxxx[k] * tke_0 - 6.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxy[k] = 6.0 * to_xyy_xy[k] - 6.0 * to_xyy_xxxy[k] * tke_0 - 4.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xxz[k] = 6.0 * to_xyy_xz[k] - 6.0 * to_xyy_xxxz[k] * tke_0 - 4.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xyy[k] = 3.0 * to_xyy_yy[k] - 6.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xyz[k] = 3.0 * to_xyy_yz[k] - 6.0 * to_xyy_xxyz[k] * tke_0 - 2.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xzz[k] = 3.0 * to_xyy_zz[k] - 6.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yyy[k] = -6.0 * to_xyy_xyyy[k] * tke_0 + 4.0 * to_xyyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yyz[k] = -6.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yzz[k] = -6.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_zzz[k] = -6.0 * to_xyy_xzzz[k] * tke_0 + 4.0 * to_xyyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 520-530 components of targeted buffer : GF

        auto to_y_x_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 70);

        auto to_y_x_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 71);

        auto to_y_x_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 72);

        auto to_y_x_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 73);

        auto to_y_x_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 74);

        auto to_y_x_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 75);

        auto to_y_x_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 76);

        auto to_y_x_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 77);

        auto to_y_x_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 78);

        auto to_y_x_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxx, to_xyyyz_xxxy, to_xyyyz_xxxz, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yz, to_xyyyz_zz, to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, to_y_x_xyyz_xxx, to_y_x_xyyz_xxy, to_y_x_xyyz_xxz, to_y_x_xyyz_xyy, to_y_x_xyyz_xyz, to_y_x_xyyz_xzz, to_y_x_xyyz_yyy, to_y_x_xyyz_yyz, to_y_x_xyyz_yzz, to_y_x_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyz_xxx[k] = 6.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxxx[k] * tke_0 - 6.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xxxy[k] * tke_0 - 4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xxz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xxxz[k] * tke_0 - 4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xyy[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xyz[k] = 2.0 * to_xyz_yz[k] - 4.0 * to_xyz_xxyz[k] * tke_0 - 2.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yyy[k] = -4.0 * to_xyz_xyyy[k] * tke_0 + 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yyz[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_zzz[k] = -4.0 * to_xyz_xzzz[k] * tke_0 + 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 530-540 components of targeted buffer : GF

        auto to_y_x_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 80);

        auto to_y_x_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 81);

        auto to_y_x_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 82);

        auto to_y_x_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 83);

        auto to_y_x_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 84);

        auto to_y_x_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 85);

        auto to_y_x_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 86);

        auto to_y_x_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 87);

        auto to_y_x_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 88);

        auto to_y_x_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyyzz_xx, to_xyyzz_xxxx, to_xyyzz_xxxy, to_xyyzz_xxxz, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yz, to_xyyzz_zz, to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, to_y_x_xyzz_xxx, to_y_x_xyzz_xxy, to_y_x_xyzz_xxz, to_y_x_xyzz_xyy, to_y_x_xyzz_xyz, to_y_x_xyzz_xzz, to_y_x_xyzz_yyy, to_y_x_xyzz_yyz, to_y_x_xyzz_yzz, to_y_x_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyzz_xxx[k] = 3.0 * to_xzz_xx[k] - 2.0 * to_xzz_xxxx[k] * tke_0 - 6.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxy[k] = 2.0 * to_xzz_xy[k] - 2.0 * to_xzz_xxxy[k] * tke_0 - 4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xxz[k] = 2.0 * to_xzz_xz[k] - 2.0 * to_xzz_xxxz[k] * tke_0 - 4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xyy[k] = to_xzz_yy[k] - 2.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xyz[k] = to_xzz_yz[k] - 2.0 * to_xzz_xxyz[k] * tke_0 - 2.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xzz[k] = to_xzz_zz[k] - 2.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yyy[k] = -2.0 * to_xzz_xyyy[k] * tke_0 + 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yyz[k] = -2.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yzz[k] = -2.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_zzz[k] = -2.0 * to_xzz_xzzz[k] * tke_0 + 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-550 components of targeted buffer : GF

        auto to_y_x_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 90);

        auto to_y_x_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 91);

        auto to_y_x_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 92);

        auto to_y_x_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 93);

        auto to_y_x_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 94);

        auto to_y_x_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 95);

        auto to_y_x_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 96);

        auto to_y_x_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 97);

        auto to_y_x_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 98);

        auto to_y_x_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xyzzz_xx, to_xyzzz_xxxx, to_xyzzz_xxxy, to_xyzzz_xxxz, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yz, to_xyzzz_zz, to_y_x_xzzz_xxx, to_y_x_xzzz_xxy, to_y_x_xzzz_xxz, to_y_x_xzzz_xyy, to_y_x_xzzz_xyz, to_y_x_xzzz_xzz, to_y_x_xzzz_yyy, to_y_x_xzzz_yyz, to_y_x_xzzz_yzz, to_y_x_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzzz_xxx[k] = -6.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxy[k] = -4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xxz[k] = -4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xyy[k] = -2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xyz[k] = -2.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xzz[k] = -2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yyy[k] = 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yyz[k] = 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yzz[k] = 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_zzz[k] = 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 550-560 components of targeted buffer : GF

        auto to_y_x_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 100);

        auto to_y_x_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 101);

        auto to_y_x_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 102);

        auto to_y_x_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 103);

        auto to_y_x_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 104);

        auto to_y_x_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 105);

        auto to_y_x_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 106);

        auto to_y_x_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 107);

        auto to_y_x_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 108);

        auto to_y_x_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_y_x_yyyy_xxx, to_y_x_yyyy_xxy, to_y_x_yyyy_xxz, to_y_x_yyyy_xyy, to_y_x_yyyy_xyz, to_y_x_yyyy_xzz, to_y_x_yyyy_yyy, to_y_x_yyyy_yyz, to_y_x_yyyy_yzz, to_y_x_yyyy_zzz, to_yyy_xx, to_yyy_xxxx, to_yyy_xxxy, to_yyy_xxxz, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yz, to_yyy_zz, to_yyyyy_xx, to_yyyyy_xxxx, to_yyyyy_xxxy, to_yyyyy_xxxz, to_yyyyy_xxyy, to_yyyyy_xxyz, to_yyyyy_xxzz, to_yyyyy_xy, to_yyyyy_xyyy, to_yyyyy_xyyz, to_yyyyy_xyzz, to_yyyyy_xz, to_yyyyy_xzzz, to_yyyyy_yy, to_yyyyy_yz, to_yyyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyy_xxx[k] = 12.0 * to_yyy_xx[k] - 8.0 * to_yyy_xxxx[k] * tke_0 - 6.0 * to_yyyyy_xx[k] * tbe_0 + 4.0 * to_yyyyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxy[k] = 8.0 * to_yyy_xy[k] - 8.0 * to_yyy_xxxy[k] * tke_0 - 4.0 * to_yyyyy_xy[k] * tbe_0 + 4.0 * to_yyyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xxz[k] = 8.0 * to_yyy_xz[k] - 8.0 * to_yyy_xxxz[k] * tke_0 - 4.0 * to_yyyyy_xz[k] * tbe_0 + 4.0 * to_yyyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xyy[k] = 4.0 * to_yyy_yy[k] - 8.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_yyyyy_yy[k] * tbe_0 + 4.0 * to_yyyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xyz[k] = 4.0 * to_yyy_yz[k] - 8.0 * to_yyy_xxyz[k] * tke_0 - 2.0 * to_yyyyy_yz[k] * tbe_0 + 4.0 * to_yyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xzz[k] = 4.0 * to_yyy_zz[k] - 8.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_yyyyy_zz[k] * tbe_0 + 4.0 * to_yyyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yyy[k] = -8.0 * to_yyy_xyyy[k] * tke_0 + 4.0 * to_yyyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yyz[k] = -8.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_yyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yzz[k] = -8.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_yyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_zzz[k] = -8.0 * to_yyy_xzzz[k] * tke_0 + 4.0 * to_yyyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 560-570 components of targeted buffer : GF

        auto to_y_x_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 110);

        auto to_y_x_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 111);

        auto to_y_x_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 112);

        auto to_y_x_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 113);

        auto to_y_x_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 114);

        auto to_y_x_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 115);

        auto to_y_x_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 116);

        auto to_y_x_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 117);

        auto to_y_x_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 118);

        auto to_y_x_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_x_yyyz_xxx, to_y_x_yyyz_xxy, to_y_x_yyyz_xxz, to_y_x_yyyz_xyy, to_y_x_yyyz_xyz, to_y_x_yyyz_xzz, to_y_x_yyyz_yyy, to_y_x_yyyz_yyz, to_y_x_yyyz_yzz, to_y_x_yyyz_zzz, to_yyyyz_xx, to_yyyyz_xxxx, to_yyyyz_xxxy, to_yyyyz_xxxz, to_yyyyz_xxyy, to_yyyyz_xxyz, to_yyyyz_xxzz, to_yyyyz_xy, to_yyyyz_xyyy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_xzzz, to_yyyyz_yy, to_yyyyz_yz, to_yyyyz_zz, to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyz_xxx[k] = 9.0 * to_yyz_xx[k] - 6.0 * to_yyz_xxxx[k] * tke_0 - 6.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxy[k] = 6.0 * to_yyz_xy[k] - 6.0 * to_yyz_xxxy[k] * tke_0 - 4.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xxz[k] = 6.0 * to_yyz_xz[k] - 6.0 * to_yyz_xxxz[k] * tke_0 - 4.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xyy[k] = 3.0 * to_yyz_yy[k] - 6.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xyz[k] = 3.0 * to_yyz_yz[k] - 6.0 * to_yyz_xxyz[k] * tke_0 - 2.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xzz[k] = 3.0 * to_yyz_zz[k] - 6.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yyy[k] = -6.0 * to_yyz_xyyy[k] * tke_0 + 4.0 * to_yyyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yyz[k] = -6.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yzz[k] = -6.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_zzz[k] = -6.0 * to_yyz_xzzz[k] * tke_0 + 4.0 * to_yyyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-580 components of targeted buffer : GF

        auto to_y_x_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 120);

        auto to_y_x_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 121);

        auto to_y_x_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 122);

        auto to_y_x_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 123);

        auto to_y_x_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 124);

        auto to_y_x_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 125);

        auto to_y_x_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 126);

        auto to_y_x_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 127);

        auto to_y_x_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 128);

        auto to_y_x_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_y_x_yyzz_xxx, to_y_x_yyzz_xxy, to_y_x_yyzz_xxz, to_y_x_yyzz_xyy, to_y_x_yyzz_xyz, to_y_x_yyzz_xzz, to_y_x_yyzz_yyy, to_y_x_yyzz_yyz, to_y_x_yyzz_yzz, to_y_x_yyzz_zzz, to_yyyzz_xx, to_yyyzz_xxxx, to_yyyzz_xxxy, to_yyyzz_xxxz, to_yyyzz_xxyy, to_yyyzz_xxyz, to_yyyzz_xxzz, to_yyyzz_xy, to_yyyzz_xyyy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_xzzz, to_yyyzz_yy, to_yyyzz_yz, to_yyyzz_zz, to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyzz_xxx[k] = 6.0 * to_yzz_xx[k] - 4.0 * to_yzz_xxxx[k] * tke_0 - 6.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxy[k] = 4.0 * to_yzz_xy[k] - 4.0 * to_yzz_xxxy[k] * tke_0 - 4.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xxz[k] = 4.0 * to_yzz_xz[k] - 4.0 * to_yzz_xxxz[k] * tke_0 - 4.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xyy[k] = 2.0 * to_yzz_yy[k] - 4.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xyz[k] = 2.0 * to_yzz_yz[k] - 4.0 * to_yzz_xxyz[k] * tke_0 - 2.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xzz[k] = 2.0 * to_yzz_zz[k] - 4.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yyy[k] = -4.0 * to_yzz_xyyy[k] * tke_0 + 4.0 * to_yyyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yyz[k] = -4.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yzz[k] = -4.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_zzz[k] = -4.0 * to_yzz_xzzz[k] * tke_0 + 4.0 * to_yyyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 580-590 components of targeted buffer : GF

        auto to_y_x_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 130);

        auto to_y_x_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 131);

        auto to_y_x_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 132);

        auto to_y_x_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 133);

        auto to_y_x_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 134);

        auto to_y_x_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 135);

        auto to_y_x_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 136);

        auto to_y_x_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 137);

        auto to_y_x_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 138);

        auto to_y_x_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_y_x_yzzz_xxx, to_y_x_yzzz_xxy, to_y_x_yzzz_xxz, to_y_x_yzzz_xyy, to_y_x_yzzz_xyz, to_y_x_yzzz_xzz, to_y_x_yzzz_yyy, to_y_x_yzzz_yyz, to_y_x_yzzz_yzz, to_y_x_yzzz_zzz, to_yyzzz_xx, to_yyzzz_xxxx, to_yyzzz_xxxy, to_yyzzz_xxxz, to_yyzzz_xxyy, to_yyzzz_xxyz, to_yyzzz_xxzz, to_yyzzz_xy, to_yyzzz_xyyy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_xzzz, to_yyzzz_yy, to_yyzzz_yz, to_yyzzz_zz, to_zzz_xx, to_zzz_xxxx, to_zzz_xxxy, to_zzz_xxxz, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzzz_xxx[k] = 3.0 * to_zzz_xx[k] - 2.0 * to_zzz_xxxx[k] * tke_0 - 6.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxy[k] = 2.0 * to_zzz_xy[k] - 2.0 * to_zzz_xxxy[k] * tke_0 - 4.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xxz[k] = 2.0 * to_zzz_xz[k] - 2.0 * to_zzz_xxxz[k] * tke_0 - 4.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xyy[k] = to_zzz_yy[k] - 2.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xyz[k] = to_zzz_yz[k] - 2.0 * to_zzz_xxyz[k] * tke_0 - 2.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xzz[k] = to_zzz_zz[k] - 2.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yyy[k] = -2.0 * to_zzz_xyyy[k] * tke_0 + 4.0 * to_yyzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yyz[k] = -2.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yzz[k] = -2.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_zzz[k] = -2.0 * to_zzz_xzzz[k] * tke_0 + 4.0 * to_yyzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 590-600 components of targeted buffer : GF

        auto to_y_x_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 140);

        auto to_y_x_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 141);

        auto to_y_x_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 142);

        auto to_y_x_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 143);

        auto to_y_x_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 144);

        auto to_y_x_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 145);

        auto to_y_x_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 146);

        auto to_y_x_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 147);

        auto to_y_x_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 148);

        auto to_y_x_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 3 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_x_zzzz_xxx, to_y_x_zzzz_xxy, to_y_x_zzzz_xxz, to_y_x_zzzz_xyy, to_y_x_zzzz_xyz, to_y_x_zzzz_xzz, to_y_x_zzzz_yyy, to_y_x_zzzz_yyz, to_y_x_zzzz_yzz, to_y_x_zzzz_zzz, to_yzzzz_xx, to_yzzzz_xxxx, to_yzzzz_xxxy, to_yzzzz_xxxz, to_yzzzz_xxyy, to_yzzzz_xxyz, to_yzzzz_xxzz, to_yzzzz_xy, to_yzzzz_xyyy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_xzzz, to_yzzzz_yy, to_yzzzz_yz, to_yzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzzz_xxx[k] = -6.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxy[k] = -4.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xxz[k] = -4.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xyy[k] = -2.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xyz[k] = -2.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xzz[k] = -2.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yyy[k] = 4.0 * to_yzzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yyz[k] = 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yzz[k] = 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_zzz[k] = 4.0 * to_yzzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-610 components of targeted buffer : GF

        auto to_y_y_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 0);

        auto to_y_y_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 1);

        auto to_y_y_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 2);

        auto to_y_y_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 3);

        auto to_y_y_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 4);

        auto to_y_y_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 5);

        auto to_y_y_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 6);

        auto to_y_y_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 7);

        auto to_y_y_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 8);

        auto to_y_y_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxy_xx, to_xxxxy_xxxy, to_xxxxy_xxyy, to_xxxxy_xxyz, to_xxxxy_xy, to_xxxxy_xyyy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_yy, to_xxxxy_yyyy, to_xxxxy_yyyz, to_xxxxy_yyzz, to_xxxxy_yz, to_xxxxy_yzzz, to_xxxxy_zz, to_y_y_xxxx_xxx, to_y_y_xxxx_xxy, to_y_y_xxxx_xxz, to_y_y_xxxx_xyy, to_y_y_xxxx_xyz, to_y_y_xxxx_xzz, to_y_y_xxxx_yyy, to_y_y_xxxx_yyz, to_y_y_xxxx_yzz, to_y_y_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxx_xxx[k] = 4.0 * to_xxxxy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxy[k] = -2.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xxz[k] = 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xyy[k] = -4.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xyz[k] = -2.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xzz[k] = 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yyy[k] = -6.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yyz[k] = -4.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yzz[k] = -2.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_zzz[k] = 4.0 * to_xxxxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 610-620 components of targeted buffer : GF

        auto to_y_y_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 10);

        auto to_y_y_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 11);

        auto to_y_y_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 12);

        auto to_y_y_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 13);

        auto to_y_y_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 14);

        auto to_y_y_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 15);

        auto to_y_y_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 16);

        auto to_y_y_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 17);

        auto to_y_y_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 18);

        auto to_y_y_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxy, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_yy, to_xxx_yyyy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxxyy_xx, to_xxxyy_xxxy, to_xxxyy_xxyy, to_xxxyy_xxyz, to_xxxyy_xy, to_xxxyy_xyyy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_yy, to_xxxyy_yyyy, to_xxxyy_yyyz, to_xxxyy_yyzz, to_xxxyy_yz, to_xxxyy_yzzz, to_xxxyy_zz, to_y_y_xxxy_xxx, to_y_y_xxxy_xxy, to_y_y_xxxy_xxz, to_y_y_xxxy_xyy, to_y_y_xxxy_xyz, to_y_y_xxxy_xzz, to_y_y_xxxy_yyy, to_y_y_xxxy_yyz, to_y_y_xxxy_yzz, to_y_y_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxy_xxx[k] = -2.0 * to_xxx_xxxy[k] * tke_0 + 4.0 * to_xxxyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxy[k] = to_xxx_xx[k] - 2.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xxz[k] = -2.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xyy[k] = 2.0 * to_xxx_xy[k] - 2.0 * to_xxx_xyyy[k] * tke_0 - 4.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xyz[k] = to_xxx_xz[k] - 2.0 * to_xxx_xyyz[k] * tke_0 - 2.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xzz[k] = -2.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yyy[k] = 3.0 * to_xxx_yy[k] - 2.0 * to_xxx_yyyy[k] * tke_0 - 6.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yyz[k] = 2.0 * to_xxx_yz[k] - 2.0 * to_xxx_yyyz[k] * tke_0 - 4.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yzz[k] = to_xxx_zz[k] - 2.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_zzz[k] = -2.0 * to_xxx_yzzz[k] * tke_0 + 4.0 * to_xxxyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 620-630 components of targeted buffer : GF

        auto to_y_y_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 20);

        auto to_y_y_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 21);

        auto to_y_y_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 22);

        auto to_y_y_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 23);

        auto to_y_y_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 24);

        auto to_y_y_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 25);

        auto to_y_y_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 26);

        auto to_y_y_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 27);

        auto to_y_y_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 28);

        auto to_y_y_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxy, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_yy, to_xxxyz_yyyy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_y_y_xxxz_xxx, to_y_y_xxxz_xxy, to_y_y_xxxz_xxz, to_y_y_xxxz_xyy, to_y_y_xxxz_xyz, to_y_y_xxxz_xzz, to_y_y_xxxz_yyy, to_y_y_xxxz_yyz, to_y_y_xxxz_yzz, to_y_y_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxz_xxx[k] = 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxy[k] = -2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xxz[k] = 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xyy[k] = -4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xyz[k] = -2.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xzz[k] = 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yyy[k] = -6.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yyz[k] = -4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yzz[k] = -2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_zzz[k] = 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-640 components of targeted buffer : GF

        auto to_y_y_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 30);

        auto to_y_y_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 31);

        auto to_y_y_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 32);

        auto to_y_y_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 33);

        auto to_y_y_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 34);

        auto to_y_y_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 35);

        auto to_y_y_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 36);

        auto to_y_y_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 37);

        auto to_y_y_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 38);

        auto to_y_y_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxyyy_xx, to_xxyyy_xxxy, to_xxyyy_xxyy, to_xxyyy_xxyz, to_xxyyy_xy, to_xxyyy_xyyy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_yy, to_xxyyy_yyyy, to_xxyyy_yyyz, to_xxyyy_yyzz, to_xxyyy_yz, to_xxyyy_yzzz, to_xxyyy_zz, to_y_y_xxyy_xxx, to_y_y_xxyy_xxy, to_y_y_xxyy_xxz, to_y_y_xxyy_xyy, to_y_y_xxyy_xyz, to_y_y_xxyy_xzz, to_y_y_xxyy_yyy, to_y_y_xxyy_yyz, to_y_y_xxyy_yzz, to_y_y_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyy_xxx[k] = -4.0 * to_xxy_xxxy[k] * tke_0 + 4.0 * to_xxyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxy[k] = 2.0 * to_xxy_xx[k] - 4.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xxz[k] = -4.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xyy[k] = 4.0 * to_xxy_xy[k] - 4.0 * to_xxy_xyyy[k] * tke_0 - 4.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xyz[k] = 2.0 * to_xxy_xz[k] - 4.0 * to_xxy_xyyz[k] * tke_0 - 2.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xzz[k] = -4.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yyy[k] = 6.0 * to_xxy_yy[k] - 4.0 * to_xxy_yyyy[k] * tke_0 - 6.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yyz[k] = 4.0 * to_xxy_yz[k] - 4.0 * to_xxy_yyyz[k] * tke_0 - 4.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yzz[k] = 2.0 * to_xxy_zz[k] - 4.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_zzz[k] = -4.0 * to_xxy_yzzz[k] * tke_0 + 4.0 * to_xxyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 640-650 components of targeted buffer : GF

        auto to_y_y_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 40);

        auto to_y_y_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 41);

        auto to_y_y_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 42);

        auto to_y_y_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 43);

        auto to_y_y_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 44);

        auto to_y_y_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 45);

        auto to_y_y_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 46);

        auto to_y_y_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 47);

        auto to_y_y_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 48);

        auto to_y_y_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxy, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_yy, to_xxyyz_yyyy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_y_y_xxyz_xxx, to_y_y_xxyz_xxy, to_y_y_xxyz_xxz, to_y_y_xxyz_xyy, to_y_y_xxyz_xyz, to_y_y_xxyz_xzz, to_y_y_xxyz_yyy, to_y_y_xxyz_yyz, to_y_y_xxyz_yzz, to_y_y_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyz_xxx[k] = -2.0 * to_xxz_xxxy[k] * tke_0 + 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxy[k] = to_xxz_xx[k] - 2.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xxz[k] = -2.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xyy[k] = 2.0 * to_xxz_xy[k] - 2.0 * to_xxz_xyyy[k] * tke_0 - 4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xyz[k] = to_xxz_xz[k] - 2.0 * to_xxz_xyyz[k] * tke_0 - 2.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xzz[k] = -2.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yyy[k] = 3.0 * to_xxz_yy[k] - 2.0 * to_xxz_yyyy[k] * tke_0 - 6.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yyz[k] = 2.0 * to_xxz_yz[k] - 2.0 * to_xxz_yyyz[k] * tke_0 - 4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yzz[k] = to_xxz_zz[k] - 2.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_zzz[k] = -2.0 * to_xxz_yzzz[k] * tke_0 + 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 650-660 components of targeted buffer : GF

        auto to_y_y_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 50);

        auto to_y_y_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 51);

        auto to_y_y_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 52);

        auto to_y_y_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 53);

        auto to_y_y_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 54);

        auto to_y_y_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 55);

        auto to_y_y_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 56);

        auto to_y_y_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 57);

        auto to_y_y_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 58);

        auto to_y_y_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxyzz_xx, to_xxyzz_xxxy, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_yy, to_xxyzz_yyyy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_y_y_xxzz_xxx, to_y_y_xxzz_xxy, to_y_y_xxzz_xxz, to_y_y_xxzz_xyy, to_y_y_xxzz_xyz, to_y_y_xxzz_xzz, to_y_y_xxzz_yyy, to_y_y_xxzz_yyz, to_y_y_xxzz_yzz, to_y_y_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxzz_xxx[k] = 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxy[k] = -2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xxz[k] = 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xyy[k] = -4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xyz[k] = -2.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xzz[k] = 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yyy[k] = -6.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yyz[k] = -4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yzz[k] = -2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_zzz[k] = 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-670 components of targeted buffer : GF

        auto to_y_y_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 60);

        auto to_y_y_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 61);

        auto to_y_y_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 62);

        auto to_y_y_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 63);

        auto to_y_y_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 64);

        auto to_y_y_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 65);

        auto to_y_y_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 66);

        auto to_y_y_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 67);

        auto to_y_y_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 68);

        auto to_y_y_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyyyy_xx, to_xyyyy_xxxy, to_xyyyy_xxyy, to_xyyyy_xxyz, to_xyyyy_xy, to_xyyyy_xyyy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_yy, to_xyyyy_yyyy, to_xyyyy_yyyz, to_xyyyy_yyzz, to_xyyyy_yz, to_xyyyy_yzzz, to_xyyyy_zz, to_y_y_xyyy_xxx, to_y_y_xyyy_xxy, to_y_y_xyyy_xxz, to_y_y_xyyy_xyy, to_y_y_xyyy_xyz, to_y_y_xyyy_xzz, to_y_y_xyyy_yyy, to_y_y_xyyy_yyz, to_y_y_xyyy_yzz, to_y_y_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyy_xxx[k] = -6.0 * to_xyy_xxxy[k] * tke_0 + 4.0 * to_xyyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxy[k] = 3.0 * to_xyy_xx[k] - 6.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xxz[k] = -6.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xyy[k] = 6.0 * to_xyy_xy[k] - 6.0 * to_xyy_xyyy[k] * tke_0 - 4.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xyz[k] = 3.0 * to_xyy_xz[k] - 6.0 * to_xyy_xyyz[k] * tke_0 - 2.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xzz[k] = -6.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yyy[k] = 9.0 * to_xyy_yy[k] - 6.0 * to_xyy_yyyy[k] * tke_0 - 6.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yyz[k] = 6.0 * to_xyy_yz[k] - 6.0 * to_xyy_yyyz[k] * tke_0 - 4.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yzz[k] = 3.0 * to_xyy_zz[k] - 6.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_zzz[k] = -6.0 * to_xyy_yzzz[k] * tke_0 + 4.0 * to_xyyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 670-680 components of targeted buffer : GF

        auto to_y_y_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 70);

        auto to_y_y_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 71);

        auto to_y_y_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 72);

        auto to_y_y_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 73);

        auto to_y_y_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 74);

        auto to_y_y_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 75);

        auto to_y_y_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 76);

        auto to_y_y_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 77);

        auto to_y_y_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 78);

        auto to_y_y_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxy, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_yy, to_xyyyz_yyyy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_y_y_xyyz_xxx, to_y_y_xyyz_xxy, to_y_y_xyyz_xxz, to_y_y_xyyz_xyy, to_y_y_xyyz_xyz, to_y_y_xyyz_xzz, to_y_y_xyyz_yyy, to_y_y_xyyz_yyz, to_y_y_xyyz_yzz, to_y_y_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyz_xxx[k] = -4.0 * to_xyz_xxxy[k] * tke_0 + 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxy[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xxz[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xyy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyyy[k] * tke_0 - 4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xyz[k] = 2.0 * to_xyz_xz[k] - 4.0 * to_xyz_xyyz[k] * tke_0 - 2.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yyy[k] = 6.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyyy[k] * tke_0 - 6.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yyz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yyyz[k] * tke_0 - 4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_zzz[k] = -4.0 * to_xyz_yzzz[k] * tke_0 + 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 680-690 components of targeted buffer : GF

        auto to_y_y_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 80);

        auto to_y_y_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 81);

        auto to_y_y_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 82);

        auto to_y_y_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 83);

        auto to_y_y_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 84);

        auto to_y_y_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 85);

        auto to_y_y_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 86);

        auto to_y_y_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 87);

        auto to_y_y_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 88);

        auto to_y_y_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyyzz_xx, to_xyyzz_xxxy, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_yy, to_xyyzz_yyyy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_y_y_xyzz_xxx, to_y_y_xyzz_xxy, to_y_y_xyzz_xxz, to_y_y_xyzz_xyy, to_y_y_xyzz_xyz, to_y_y_xyzz_xzz, to_y_y_xyzz_yyy, to_y_y_xyzz_yyz, to_y_y_xyzz_yzz, to_y_y_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyzz_xxx[k] = -2.0 * to_xzz_xxxy[k] * tke_0 + 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxy[k] = to_xzz_xx[k] - 2.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xxz[k] = -2.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xyy[k] = 2.0 * to_xzz_xy[k] - 2.0 * to_xzz_xyyy[k] * tke_0 - 4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xyz[k] = to_xzz_xz[k] - 2.0 * to_xzz_xyyz[k] * tke_0 - 2.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xzz[k] = -2.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yyy[k] = 3.0 * to_xzz_yy[k] - 2.0 * to_xzz_yyyy[k] * tke_0 - 6.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yyz[k] = 2.0 * to_xzz_yz[k] - 2.0 * to_xzz_yyyz[k] * tke_0 - 4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yzz[k] = to_xzz_zz[k] - 2.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_zzz[k] = -2.0 * to_xzz_yzzz[k] * tke_0 + 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-700 components of targeted buffer : GF

        auto to_y_y_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 90);

        auto to_y_y_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 91);

        auto to_y_y_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 92);

        auto to_y_y_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 93);

        auto to_y_y_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 94);

        auto to_y_y_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 95);

        auto to_y_y_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 96);

        auto to_y_y_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 97);

        auto to_y_y_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 98);

        auto to_y_y_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xyzzz_xx, to_xyzzz_xxxy, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_yy, to_xyzzz_yyyy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, to_y_y_xzzz_xxx, to_y_y_xzzz_xxy, to_y_y_xzzz_xxz, to_y_y_xzzz_xyy, to_y_y_xzzz_xyz, to_y_y_xzzz_xzz, to_y_y_xzzz_yyy, to_y_y_xzzz_yyz, to_y_y_xzzz_yzz, to_y_y_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzzz_xxx[k] = 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxy[k] = -2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xxz[k] = 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xyy[k] = -4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xyz[k] = -2.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xzz[k] = 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yyy[k] = -6.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yyz[k] = -4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yzz[k] = -2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_zzz[k] = 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 700-710 components of targeted buffer : GF

        auto to_y_y_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 100);

        auto to_y_y_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 101);

        auto to_y_y_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 102);

        auto to_y_y_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 103);

        auto to_y_y_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 104);

        auto to_y_y_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 105);

        auto to_y_y_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 106);

        auto to_y_y_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 107);

        auto to_y_y_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 108);

        auto to_y_y_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_y_y_yyyy_xxx, to_y_y_yyyy_xxy, to_y_y_yyyy_xxz, to_y_y_yyyy_xyy, to_y_y_yyyy_xyz, to_y_y_yyyy_xzz, to_y_y_yyyy_yyy, to_y_y_yyyy_yyz, to_y_y_yyyy_yzz, to_y_y_yyyy_zzz, to_yyy_xx, to_yyy_xxxy, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_yy, to_yyy_yyyy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyyyy_xx, to_yyyyy_xxxy, to_yyyyy_xxyy, to_yyyyy_xxyz, to_yyyyy_xy, to_yyyyy_xyyy, to_yyyyy_xyyz, to_yyyyy_xyzz, to_yyyyy_xz, to_yyyyy_yy, to_yyyyy_yyyy, to_yyyyy_yyyz, to_yyyyy_yyzz, to_yyyyy_yz, to_yyyyy_yzzz, to_yyyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyy_xxx[k] = -8.0 * to_yyy_xxxy[k] * tke_0 + 4.0 * to_yyyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxy[k] = 4.0 * to_yyy_xx[k] - 8.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_yyyyy_xx[k] * tbe_0 + 4.0 * to_yyyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xxz[k] = -8.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_yyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xyy[k] = 8.0 * to_yyy_xy[k] - 8.0 * to_yyy_xyyy[k] * tke_0 - 4.0 * to_yyyyy_xy[k] * tbe_0 + 4.0 * to_yyyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xyz[k] = 4.0 * to_yyy_xz[k] - 8.0 * to_yyy_xyyz[k] * tke_0 - 2.0 * to_yyyyy_xz[k] * tbe_0 + 4.0 * to_yyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xzz[k] = -8.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_yyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yyy[k] = 12.0 * to_yyy_yy[k] - 8.0 * to_yyy_yyyy[k] * tke_0 - 6.0 * to_yyyyy_yy[k] * tbe_0 + 4.0 * to_yyyyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yyz[k] = 8.0 * to_yyy_yz[k] - 8.0 * to_yyy_yyyz[k] * tke_0 - 4.0 * to_yyyyy_yz[k] * tbe_0 + 4.0 * to_yyyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yzz[k] = 4.0 * to_yyy_zz[k] - 8.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_yyyyy_zz[k] * tbe_0 + 4.0 * to_yyyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_zzz[k] = -8.0 * to_yyy_yzzz[k] * tke_0 + 4.0 * to_yyyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 710-720 components of targeted buffer : GF

        auto to_y_y_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 110);

        auto to_y_y_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 111);

        auto to_y_y_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 112);

        auto to_y_y_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 113);

        auto to_y_y_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 114);

        auto to_y_y_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 115);

        auto to_y_y_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 116);

        auto to_y_y_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 117);

        auto to_y_y_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 118);

        auto to_y_y_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_y_yyyz_xxx, to_y_y_yyyz_xxy, to_y_y_yyyz_xxz, to_y_y_yyyz_xyy, to_y_y_yyyz_xyz, to_y_y_yyyz_xzz, to_y_y_yyyz_yyy, to_y_y_yyyz_yyz, to_y_y_yyyz_yzz, to_y_y_yyyz_zzz, to_yyyyz_xx, to_yyyyz_xxxy, to_yyyyz_xxyy, to_yyyyz_xxyz, to_yyyyz_xy, to_yyyyz_xyyy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_yy, to_yyyyz_yyyy, to_yyyyz_yyyz, to_yyyyz_yyzz, to_yyyyz_yz, to_yyyyz_yzzz, to_yyyyz_zz, to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyz_xxx[k] = -6.0 * to_yyz_xxxy[k] * tke_0 + 4.0 * to_yyyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxy[k] = 3.0 * to_yyz_xx[k] - 6.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xxz[k] = -6.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xyy[k] = 6.0 * to_yyz_xy[k] - 6.0 * to_yyz_xyyy[k] * tke_0 - 4.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xyz[k] = 3.0 * to_yyz_xz[k] - 6.0 * to_yyz_xyyz[k] * tke_0 - 2.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xzz[k] = -6.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yyy[k] = 9.0 * to_yyz_yy[k] - 6.0 * to_yyz_yyyy[k] * tke_0 - 6.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yyz[k] = 6.0 * to_yyz_yz[k] - 6.0 * to_yyz_yyyz[k] * tke_0 - 4.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yzz[k] = 3.0 * to_yyz_zz[k] - 6.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_zzz[k] = -6.0 * to_yyz_yzzz[k] * tke_0 + 4.0 * to_yyyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-730 components of targeted buffer : GF

        auto to_y_y_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 120);

        auto to_y_y_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 121);

        auto to_y_y_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 122);

        auto to_y_y_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 123);

        auto to_y_y_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 124);

        auto to_y_y_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 125);

        auto to_y_y_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 126);

        auto to_y_y_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 127);

        auto to_y_y_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 128);

        auto to_y_y_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_y_y_yyzz_xxx, to_y_y_yyzz_xxy, to_y_y_yyzz_xxz, to_y_y_yyzz_xyy, to_y_y_yyzz_xyz, to_y_y_yyzz_xzz, to_y_y_yyzz_yyy, to_y_y_yyzz_yyz, to_y_y_yyzz_yzz, to_y_y_yyzz_zzz, to_yyyzz_xx, to_yyyzz_xxxy, to_yyyzz_xxyy, to_yyyzz_xxyz, to_yyyzz_xy, to_yyyzz_xyyy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_yy, to_yyyzz_yyyy, to_yyyzz_yyyz, to_yyyzz_yyzz, to_yyyzz_yz, to_yyyzz_yzzz, to_yyyzz_zz, to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyzz_xxx[k] = -4.0 * to_yzz_xxxy[k] * tke_0 + 4.0 * to_yyyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxy[k] = 2.0 * to_yzz_xx[k] - 4.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xxz[k] = -4.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xyy[k] = 4.0 * to_yzz_xy[k] - 4.0 * to_yzz_xyyy[k] * tke_0 - 4.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xyz[k] = 2.0 * to_yzz_xz[k] - 4.0 * to_yzz_xyyz[k] * tke_0 - 2.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xzz[k] = -4.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yyy[k] = 6.0 * to_yzz_yy[k] - 4.0 * to_yzz_yyyy[k] * tke_0 - 6.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yyz[k] = 4.0 * to_yzz_yz[k] - 4.0 * to_yzz_yyyz[k] * tke_0 - 4.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yzz[k] = 2.0 * to_yzz_zz[k] - 4.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_zzz[k] = -4.0 * to_yzz_yzzz[k] * tke_0 + 4.0 * to_yyyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 730-740 components of targeted buffer : GF

        auto to_y_y_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 130);

        auto to_y_y_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 131);

        auto to_y_y_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 132);

        auto to_y_y_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 133);

        auto to_y_y_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 134);

        auto to_y_y_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 135);

        auto to_y_y_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 136);

        auto to_y_y_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 137);

        auto to_y_y_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 138);

        auto to_y_y_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_y_y_yzzz_xxx, to_y_y_yzzz_xxy, to_y_y_yzzz_xxz, to_y_y_yzzz_xyy, to_y_y_yzzz_xyz, to_y_y_yzzz_xzz, to_y_y_yzzz_yyy, to_y_y_yzzz_yyz, to_y_y_yzzz_yzz, to_y_y_yzzz_zzz, to_yyzzz_xx, to_yyzzz_xxxy, to_yyzzz_xxyy, to_yyzzz_xxyz, to_yyzzz_xy, to_yyzzz_xyyy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_yy, to_yyzzz_yyyy, to_yyzzz_yyyz, to_yyzzz_yyzz, to_yyzzz_yz, to_yyzzz_yzzz, to_yyzzz_zz, to_zzz_xx, to_zzz_xxxy, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_yy, to_zzz_yyyy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzzz_xxx[k] = -2.0 * to_zzz_xxxy[k] * tke_0 + 4.0 * to_yyzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxy[k] = to_zzz_xx[k] - 2.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xxz[k] = -2.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xyy[k] = 2.0 * to_zzz_xy[k] - 2.0 * to_zzz_xyyy[k] * tke_0 - 4.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xyz[k] = to_zzz_xz[k] - 2.0 * to_zzz_xyyz[k] * tke_0 - 2.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xzz[k] = -2.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yyy[k] = 3.0 * to_zzz_yy[k] - 2.0 * to_zzz_yyyy[k] * tke_0 - 6.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yyz[k] = 2.0 * to_zzz_yz[k] - 2.0 * to_zzz_yyyz[k] * tke_0 - 4.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yzz[k] = to_zzz_zz[k] - 2.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_zzz[k] = -2.0 * to_zzz_yzzz[k] * tke_0 + 4.0 * to_yyzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 740-750 components of targeted buffer : GF

        auto to_y_y_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 140);

        auto to_y_y_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 141);

        auto to_y_y_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 142);

        auto to_y_y_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 143);

        auto to_y_y_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 144);

        auto to_y_y_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 145);

        auto to_y_y_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 146);

        auto to_y_y_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 147);

        auto to_y_y_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 148);

        auto to_y_y_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 4 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_y_zzzz_xxx, to_y_y_zzzz_xxy, to_y_y_zzzz_xxz, to_y_y_zzzz_xyy, to_y_y_zzzz_xyz, to_y_y_zzzz_xzz, to_y_y_zzzz_yyy, to_y_y_zzzz_yyz, to_y_y_zzzz_yzz, to_y_y_zzzz_zzz, to_yzzzz_xx, to_yzzzz_xxxy, to_yzzzz_xxyy, to_yzzzz_xxyz, to_yzzzz_xy, to_yzzzz_xyyy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_yy, to_yzzzz_yyyy, to_yzzzz_yyyz, to_yzzzz_yyzz, to_yzzzz_yz, to_yzzzz_yzzz, to_yzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzzz_xxx[k] = 4.0 * to_yzzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxy[k] = -2.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xxz[k] = 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xyy[k] = -4.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xyz[k] = -2.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xzz[k] = 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yyy[k] = -6.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yyz[k] = -4.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yzz[k] = -2.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_zzz[k] = 4.0 * to_yzzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-760 components of targeted buffer : GF

        auto to_y_z_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 0);

        auto to_y_z_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 1);

        auto to_y_z_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 2);

        auto to_y_z_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 3);

        auto to_y_z_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 4);

        auto to_y_z_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 5);

        auto to_y_z_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 6);

        auto to_y_z_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 7);

        auto to_y_z_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 8);

        auto to_y_z_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxy_xx, to_xxxxy_xxxz, to_xxxxy_xxyz, to_xxxxy_xxzz, to_xxxxy_xy, to_xxxxy_xyyz, to_xxxxy_xyzz, to_xxxxy_xz, to_xxxxy_xzzz, to_xxxxy_yy, to_xxxxy_yyyz, to_xxxxy_yyzz, to_xxxxy_yz, to_xxxxy_yzzz, to_xxxxy_zz, to_xxxxy_zzzz, to_y_z_xxxx_xxx, to_y_z_xxxx_xxy, to_y_z_xxxx_xxz, to_y_z_xxxx_xyy, to_y_z_xxxx_xyz, to_y_z_xxxx_xzz, to_y_z_xxxx_yyy, to_y_z_xxxx_yyz, to_y_z_xxxx_yzz, to_y_z_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxx_xxx[k] = 4.0 * to_xxxxy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxy[k] = 4.0 * to_xxxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xxz[k] = -2.0 * to_xxxxy_xx[k] * tbe_0 + 4.0 * to_xxxxy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xyy[k] = 4.0 * to_xxxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xyz[k] = -2.0 * to_xxxxy_xy[k] * tbe_0 + 4.0 * to_xxxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xzz[k] = -4.0 * to_xxxxy_xz[k] * tbe_0 + 4.0 * to_xxxxy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yyy[k] = 4.0 * to_xxxxy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yyz[k] = -2.0 * to_xxxxy_yy[k] * tbe_0 + 4.0 * to_xxxxy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yzz[k] = -4.0 * to_xxxxy_yz[k] * tbe_0 + 4.0 * to_xxxxy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_zzz[k] = -6.0 * to_xxxxy_zz[k] * tbe_0 + 4.0 * to_xxxxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 760-770 components of targeted buffer : GF

        auto to_y_z_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 10);

        auto to_y_z_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 11);

        auto to_y_z_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 12);

        auto to_y_z_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 13);

        auto to_y_z_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 14);

        auto to_y_z_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 15);

        auto to_y_z_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 16);

        auto to_y_z_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 17);

        auto to_y_z_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 18);

        auto to_y_z_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxz, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxx_zzzz, to_xxxyy_xx, to_xxxyy_xxxz, to_xxxyy_xxyz, to_xxxyy_xxzz, to_xxxyy_xy, to_xxxyy_xyyz, to_xxxyy_xyzz, to_xxxyy_xz, to_xxxyy_xzzz, to_xxxyy_yy, to_xxxyy_yyyz, to_xxxyy_yyzz, to_xxxyy_yz, to_xxxyy_yzzz, to_xxxyy_zz, to_xxxyy_zzzz, to_y_z_xxxy_xxx, to_y_z_xxxy_xxy, to_y_z_xxxy_xxz, to_y_z_xxxy_xyy, to_y_z_xxxy_xyz, to_y_z_xxxy_xzz, to_y_z_xxxy_yyy, to_y_z_xxxy_yyz, to_y_z_xxxy_yzz, to_y_z_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxy_xxx[k] = -2.0 * to_xxx_xxxz[k] * tke_0 + 4.0 * to_xxxyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxy[k] = -2.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xxz[k] = to_xxx_xx[k] - 2.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxyy_xx[k] * tbe_0 + 4.0 * to_xxxyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xyy[k] = -2.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xyz[k] = to_xxx_xy[k] - 2.0 * to_xxx_xyzz[k] * tke_0 - 2.0 * to_xxxyy_xy[k] * tbe_0 + 4.0 * to_xxxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xzz[k] = 2.0 * to_xxx_xz[k] - 2.0 * to_xxx_xzzz[k] * tke_0 - 4.0 * to_xxxyy_xz[k] * tbe_0 + 4.0 * to_xxxyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yyy[k] = -2.0 * to_xxx_yyyz[k] * tke_0 + 4.0 * to_xxxyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yyz[k] = to_xxx_yy[k] - 2.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxyy_yy[k] * tbe_0 + 4.0 * to_xxxyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yzz[k] = 2.0 * to_xxx_yz[k] - 2.0 * to_xxx_yzzz[k] * tke_0 - 4.0 * to_xxxyy_yz[k] * tbe_0 + 4.0 * to_xxxyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_zzz[k] = 3.0 * to_xxx_zz[k] - 2.0 * to_xxx_zzzz[k] * tke_0 - 6.0 * to_xxxyy_zz[k] * tbe_0 + 4.0 * to_xxxyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 770-780 components of targeted buffer : GF

        auto to_y_z_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 20);

        auto to_y_z_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 21);

        auto to_y_z_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 22);

        auto to_y_z_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 23);

        auto to_y_z_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 24);

        auto to_y_z_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 25);

        auto to_y_z_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 26);

        auto to_y_z_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 27);

        auto to_y_z_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 28);

        auto to_y_z_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxz, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_xxxyz_zzzz, to_y_z_xxxz_xxx, to_y_z_xxxz_xxy, to_y_z_xxxz_xxz, to_y_z_xxxz_xyy, to_y_z_xxxz_xyz, to_y_z_xxxz_xzz, to_y_z_xxxz_yyy, to_y_z_xxxz_yyz, to_y_z_xxxz_yzz, to_y_z_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxz_xxx[k] = 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxy[k] = 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xxz[k] = -2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xyy[k] = 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xyz[k] = -2.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xzz[k] = -4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yyy[k] = 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yyz[k] = -2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yzz[k] = -4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_zzz[k] = -6.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-790 components of targeted buffer : GF

        auto to_y_z_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 30);

        auto to_y_z_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 31);

        auto to_y_z_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 32);

        auto to_y_z_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 33);

        auto to_y_z_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 34);

        auto to_y_z_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 35);

        auto to_y_z_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 36);

        auto to_y_z_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 37);

        auto to_y_z_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 38);

        auto to_y_z_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, to_xxyyy_xx, to_xxyyy_xxxz, to_xxyyy_xxyz, to_xxyyy_xxzz, to_xxyyy_xy, to_xxyyy_xyyz, to_xxyyy_xyzz, to_xxyyy_xz, to_xxyyy_xzzz, to_xxyyy_yy, to_xxyyy_yyyz, to_xxyyy_yyzz, to_xxyyy_yz, to_xxyyy_yzzz, to_xxyyy_zz, to_xxyyy_zzzz, to_y_z_xxyy_xxx, to_y_z_xxyy_xxy, to_y_z_xxyy_xxz, to_y_z_xxyy_xyy, to_y_z_xxyy_xyz, to_y_z_xxyy_xzz, to_y_z_xxyy_yyy, to_y_z_xxyy_yyz, to_y_z_xxyy_yzz, to_y_z_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyy_xxx[k] = -4.0 * to_xxy_xxxz[k] * tke_0 + 4.0 * to_xxyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxy[k] = -4.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xxz[k] = 2.0 * to_xxy_xx[k] - 4.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxyyy_xx[k] * tbe_0 + 4.0 * to_xxyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xyy[k] = -4.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xyz[k] = 2.0 * to_xxy_xy[k] - 4.0 * to_xxy_xyzz[k] * tke_0 - 2.0 * to_xxyyy_xy[k] * tbe_0 + 4.0 * to_xxyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xzz[k] = 4.0 * to_xxy_xz[k] - 4.0 * to_xxy_xzzz[k] * tke_0 - 4.0 * to_xxyyy_xz[k] * tbe_0 + 4.0 * to_xxyyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yyy[k] = -4.0 * to_xxy_yyyz[k] * tke_0 + 4.0 * to_xxyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yyz[k] = 2.0 * to_xxy_yy[k] - 4.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxyyy_yy[k] * tbe_0 + 4.0 * to_xxyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yzz[k] = 4.0 * to_xxy_yz[k] - 4.0 * to_xxy_yzzz[k] * tke_0 - 4.0 * to_xxyyy_yz[k] * tbe_0 + 4.0 * to_xxyyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_zzz[k] = 6.0 * to_xxy_zz[k] - 4.0 * to_xxy_zzzz[k] * tke_0 - 6.0 * to_xxyyy_zz[k] * tbe_0 + 4.0 * to_xxyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 790-800 components of targeted buffer : GF

        auto to_y_z_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 40);

        auto to_y_z_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 41);

        auto to_y_z_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 42);

        auto to_y_z_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 43);

        auto to_y_z_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 44);

        auto to_y_z_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 45);

        auto to_y_z_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 46);

        auto to_y_z_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 47);

        auto to_y_z_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 48);

        auto to_y_z_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxz, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_xxyyz_zzzz, to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, to_y_z_xxyz_xxx, to_y_z_xxyz_xxy, to_y_z_xxyz_xxz, to_y_z_xxyz_xyy, to_y_z_xxyz_xyz, to_y_z_xxyz_xzz, to_y_z_xxyz_yyy, to_y_z_xxyz_yyz, to_y_z_xxyz_yzz, to_y_z_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyz_xxx[k] = -2.0 * to_xxz_xxxz[k] * tke_0 + 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxy[k] = -2.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xxz[k] = to_xxz_xx[k] - 2.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xyy[k] = -2.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xyz[k] = to_xxz_xy[k] - 2.0 * to_xxz_xyzz[k] * tke_0 - 2.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xzz[k] = 2.0 * to_xxz_xz[k] - 2.0 * to_xxz_xzzz[k] * tke_0 - 4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yyy[k] = -2.0 * to_xxz_yyyz[k] * tke_0 + 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yyz[k] = to_xxz_yy[k] - 2.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yzz[k] = 2.0 * to_xxz_yz[k] - 2.0 * to_xxz_yzzz[k] * tke_0 - 4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_zzz[k] = 3.0 * to_xxz_zz[k] - 2.0 * to_xxz_zzzz[k] * tke_0 - 6.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 800-810 components of targeted buffer : GF

        auto to_y_z_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 50);

        auto to_y_z_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 51);

        auto to_y_z_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 52);

        auto to_y_z_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 53);

        auto to_y_z_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 54);

        auto to_y_z_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 55);

        auto to_y_z_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 56);

        auto to_y_z_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 57);

        auto to_y_z_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 58);

        auto to_y_z_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxyzz_xx, to_xxyzz_xxxz, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_xxyzz_zzzz, to_y_z_xxzz_xxx, to_y_z_xxzz_xxy, to_y_z_xxzz_xxz, to_y_z_xxzz_xyy, to_y_z_xxzz_xyz, to_y_z_xxzz_xzz, to_y_z_xxzz_yyy, to_y_z_xxzz_yyz, to_y_z_xxzz_yzz, to_y_z_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxzz_xxx[k] = 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxy[k] = 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xxz[k] = -2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xyy[k] = 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xyz[k] = -2.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xzz[k] = -4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yyy[k] = 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yyz[k] = -2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yzz[k] = -4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_zzz[k] = -6.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 810-820 components of targeted buffer : GF

        auto to_y_z_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 60);

        auto to_y_z_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 61);

        auto to_y_z_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 62);

        auto to_y_z_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 63);

        auto to_y_z_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 64);

        auto to_y_z_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 65);

        auto to_y_z_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 66);

        auto to_y_z_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 67);

        auto to_y_z_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 68);

        auto to_y_z_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, to_xyyyy_xx, to_xyyyy_xxxz, to_xyyyy_xxyz, to_xyyyy_xxzz, to_xyyyy_xy, to_xyyyy_xyyz, to_xyyyy_xyzz, to_xyyyy_xz, to_xyyyy_xzzz, to_xyyyy_yy, to_xyyyy_yyyz, to_xyyyy_yyzz, to_xyyyy_yz, to_xyyyy_yzzz, to_xyyyy_zz, to_xyyyy_zzzz, to_y_z_xyyy_xxx, to_y_z_xyyy_xxy, to_y_z_xyyy_xxz, to_y_z_xyyy_xyy, to_y_z_xyyy_xyz, to_y_z_xyyy_xzz, to_y_z_xyyy_yyy, to_y_z_xyyy_yyz, to_y_z_xyyy_yzz, to_y_z_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyy_xxx[k] = -6.0 * to_xyy_xxxz[k] * tke_0 + 4.0 * to_xyyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxy[k] = -6.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xxz[k] = 3.0 * to_xyy_xx[k] - 6.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xyyyy_xx[k] * tbe_0 + 4.0 * to_xyyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xyy[k] = -6.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xyz[k] = 3.0 * to_xyy_xy[k] - 6.0 * to_xyy_xyzz[k] * tke_0 - 2.0 * to_xyyyy_xy[k] * tbe_0 + 4.0 * to_xyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xzz[k] = 6.0 * to_xyy_xz[k] - 6.0 * to_xyy_xzzz[k] * tke_0 - 4.0 * to_xyyyy_xz[k] * tbe_0 + 4.0 * to_xyyyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yyy[k] = -6.0 * to_xyy_yyyz[k] * tke_0 + 4.0 * to_xyyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yyz[k] = 3.0 * to_xyy_yy[k] - 6.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xyyyy_yy[k] * tbe_0 + 4.0 * to_xyyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yzz[k] = 6.0 * to_xyy_yz[k] - 6.0 * to_xyy_yzzz[k] * tke_0 - 4.0 * to_xyyyy_yz[k] * tbe_0 + 4.0 * to_xyyyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_zzz[k] = 9.0 * to_xyy_zz[k] - 6.0 * to_xyy_zzzz[k] * tke_0 - 6.0 * to_xyyyy_zz[k] * tbe_0 + 4.0 * to_xyyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 820-830 components of targeted buffer : GF

        auto to_y_z_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 70);

        auto to_y_z_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 71);

        auto to_y_z_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 72);

        auto to_y_z_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 73);

        auto to_y_z_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 74);

        auto to_y_z_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 75);

        auto to_y_z_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 76);

        auto to_y_z_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 77);

        auto to_y_z_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 78);

        auto to_y_z_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxz, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, to_xyyyz_zzzz, to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, to_y_z_xyyz_xxx, to_y_z_xyyz_xxy, to_y_z_xyyz_xxz, to_y_z_xyyz_xyy, to_y_z_xyyz_xyz, to_y_z_xyyz_xzz, to_y_z_xyyz_yyy, to_y_z_xyyz_yyz, to_y_z_xyyz_yzz, to_y_z_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyz_xxx[k] = -4.0 * to_xyz_xxxz[k] * tke_0 + 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxy[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xxz[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xyy[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xyz[k] = 2.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyzz[k] * tke_0 - 2.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xzz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xzzz[k] * tke_0 - 4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yyy[k] = -4.0 * to_xyz_yyyz[k] * tke_0 + 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yyz[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yzz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yzzz[k] * tke_0 - 4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_zzz[k] = 6.0 * to_xyz_zz[k] - 4.0 * to_xyz_zzzz[k] * tke_0 - 6.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 830-840 components of targeted buffer : GF

        auto to_y_z_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 80);

        auto to_y_z_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 81);

        auto to_y_z_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 82);

        auto to_y_z_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 83);

        auto to_y_z_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 84);

        auto to_y_z_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 85);

        auto to_y_z_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 86);

        auto to_y_z_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 87);

        auto to_y_z_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 88);

        auto to_y_z_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyyzz_xx, to_xyyzz_xxxz, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, to_xyyzz_zzzz, to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, to_y_z_xyzz_xxx, to_y_z_xyzz_xxy, to_y_z_xyzz_xxz, to_y_z_xyzz_xyy, to_y_z_xyzz_xyz, to_y_z_xyzz_xzz, to_y_z_xyzz_yyy, to_y_z_xyzz_yyz, to_y_z_xyzz_yzz, to_y_z_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyzz_xxx[k] = -2.0 * to_xzz_xxxz[k] * tke_0 + 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxy[k] = -2.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xxz[k] = to_xzz_xx[k] - 2.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xyy[k] = -2.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xyz[k] = to_xzz_xy[k] - 2.0 * to_xzz_xyzz[k] * tke_0 - 2.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xzz[k] = 2.0 * to_xzz_xz[k] - 2.0 * to_xzz_xzzz[k] * tke_0 - 4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yyy[k] = -2.0 * to_xzz_yyyz[k] * tke_0 + 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yyz[k] = to_xzz_yy[k] - 2.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yzz[k] = 2.0 * to_xzz_yz[k] - 2.0 * to_xzz_yzzz[k] * tke_0 - 4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_zzz[k] = 3.0 * to_xzz_zz[k] - 2.0 * to_xzz_zzzz[k] * tke_0 - 6.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 840-850 components of targeted buffer : GF

        auto to_y_z_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 90);

        auto to_y_z_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 91);

        auto to_y_z_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 92);

        auto to_y_z_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 93);

        auto to_y_z_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 94);

        auto to_y_z_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 95);

        auto to_y_z_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 96);

        auto to_y_z_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 97);

        auto to_y_z_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 98);

        auto to_y_z_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xyzzz_xx, to_xyzzz_xxxz, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, to_xyzzz_zzzz, to_y_z_xzzz_xxx, to_y_z_xzzz_xxy, to_y_z_xzzz_xxz, to_y_z_xzzz_xyy, to_y_z_xzzz_xyz, to_y_z_xzzz_xzz, to_y_z_xzzz_yyy, to_y_z_xzzz_yyz, to_y_z_xzzz_yzz, to_y_z_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzzz_xxx[k] = 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxy[k] = 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xxz[k] = -2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xyy[k] = 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xyz[k] = -2.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xzz[k] = -4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yyy[k] = 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yyz[k] = -2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yzz[k] = -4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_zzz[k] = -6.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 850-860 components of targeted buffer : GF

        auto to_y_z_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 100);

        auto to_y_z_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 101);

        auto to_y_z_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 102);

        auto to_y_z_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 103);

        auto to_y_z_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 104);

        auto to_y_z_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 105);

        auto to_y_z_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 106);

        auto to_y_z_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 107);

        auto to_y_z_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 108);

        auto to_y_z_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_y_z_yyyy_xxx, to_y_z_yyyy_xxy, to_y_z_yyyy_xxz, to_y_z_yyyy_xyy, to_y_z_yyyy_xyz, to_y_z_yyyy_xzz, to_y_z_yyyy_yyy, to_y_z_yyyy_yyz, to_y_z_yyyy_yzz, to_y_z_yyyy_zzz, to_yyy_xx, to_yyy_xxxz, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyy_zzzz, to_yyyyy_xx, to_yyyyy_xxxz, to_yyyyy_xxyz, to_yyyyy_xxzz, to_yyyyy_xy, to_yyyyy_xyyz, to_yyyyy_xyzz, to_yyyyy_xz, to_yyyyy_xzzz, to_yyyyy_yy, to_yyyyy_yyyz, to_yyyyy_yyzz, to_yyyyy_yz, to_yyyyy_yzzz, to_yyyyy_zz, to_yyyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyy_xxx[k] = -8.0 * to_yyy_xxxz[k] * tke_0 + 4.0 * to_yyyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxy[k] = -8.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_yyyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xxz[k] = 4.0 * to_yyy_xx[k] - 8.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_yyyyy_xx[k] * tbe_0 + 4.0 * to_yyyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xyy[k] = -8.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_yyyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xyz[k] = 4.0 * to_yyy_xy[k] - 8.0 * to_yyy_xyzz[k] * tke_0 - 2.0 * to_yyyyy_xy[k] * tbe_0 + 4.0 * to_yyyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xzz[k] = 8.0 * to_yyy_xz[k] - 8.0 * to_yyy_xzzz[k] * tke_0 - 4.0 * to_yyyyy_xz[k] * tbe_0 + 4.0 * to_yyyyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yyy[k] = -8.0 * to_yyy_yyyz[k] * tke_0 + 4.0 * to_yyyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yyz[k] = 4.0 * to_yyy_yy[k] - 8.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_yyyyy_yy[k] * tbe_0 + 4.0 * to_yyyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yzz[k] = 8.0 * to_yyy_yz[k] - 8.0 * to_yyy_yzzz[k] * tke_0 - 4.0 * to_yyyyy_yz[k] * tbe_0 + 4.0 * to_yyyyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_zzz[k] = 12.0 * to_yyy_zz[k] - 8.0 * to_yyy_zzzz[k] * tke_0 - 6.0 * to_yyyyy_zz[k] * tbe_0 + 4.0 * to_yyyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 860-870 components of targeted buffer : GF

        auto to_y_z_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 110);

        auto to_y_z_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 111);

        auto to_y_z_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 112);

        auto to_y_z_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 113);

        auto to_y_z_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 114);

        auto to_y_z_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 115);

        auto to_y_z_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 116);

        auto to_y_z_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 117);

        auto to_y_z_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 118);

        auto to_y_z_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_z_yyyz_xxx, to_y_z_yyyz_xxy, to_y_z_yyyz_xxz, to_y_z_yyyz_xyy, to_y_z_yyyz_xyz, to_y_z_yyyz_xzz, to_y_z_yyyz_yyy, to_y_z_yyyz_yyz, to_y_z_yyyz_yzz, to_y_z_yyyz_zzz, to_yyyyz_xx, to_yyyyz_xxxz, to_yyyyz_xxyz, to_yyyyz_xxzz, to_yyyyz_xy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_xzzz, to_yyyyz_yy, to_yyyyz_yyyz, to_yyyyz_yyzz, to_yyyyz_yz, to_yyyyz_yzzz, to_yyyyz_zz, to_yyyyz_zzzz, to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyz_xxx[k] = -6.0 * to_yyz_xxxz[k] * tke_0 + 4.0 * to_yyyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxy[k] = -6.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xxz[k] = 3.0 * to_yyz_xx[k] - 6.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xyy[k] = -6.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xyz[k] = 3.0 * to_yyz_xy[k] - 6.0 * to_yyz_xyzz[k] * tke_0 - 2.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xzz[k] = 6.0 * to_yyz_xz[k] - 6.0 * to_yyz_xzzz[k] * tke_0 - 4.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yyy[k] = -6.0 * to_yyz_yyyz[k] * tke_0 + 4.0 * to_yyyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yyz[k] = 3.0 * to_yyz_yy[k] - 6.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yzz[k] = 6.0 * to_yyz_yz[k] - 6.0 * to_yyz_yzzz[k] * tke_0 - 4.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_zzz[k] = 9.0 * to_yyz_zz[k] - 6.0 * to_yyz_zzzz[k] * tke_0 - 6.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 870-880 components of targeted buffer : GF

        auto to_y_z_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 120);

        auto to_y_z_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 121);

        auto to_y_z_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 122);

        auto to_y_z_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 123);

        auto to_y_z_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 124);

        auto to_y_z_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 125);

        auto to_y_z_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 126);

        auto to_y_z_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 127);

        auto to_y_z_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 128);

        auto to_y_z_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_y_z_yyzz_xxx, to_y_z_yyzz_xxy, to_y_z_yyzz_xxz, to_y_z_yyzz_xyy, to_y_z_yyzz_xyz, to_y_z_yyzz_xzz, to_y_z_yyzz_yyy, to_y_z_yyzz_yyz, to_y_z_yyzz_yzz, to_y_z_yyzz_zzz, to_yyyzz_xx, to_yyyzz_xxxz, to_yyyzz_xxyz, to_yyyzz_xxzz, to_yyyzz_xy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_xzzz, to_yyyzz_yy, to_yyyzz_yyyz, to_yyyzz_yyzz, to_yyyzz_yz, to_yyyzz_yzzz, to_yyyzz_zz, to_yyyzz_zzzz, to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyzz_xxx[k] = -4.0 * to_yzz_xxxz[k] * tke_0 + 4.0 * to_yyyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxy[k] = -4.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xxz[k] = 2.0 * to_yzz_xx[k] - 4.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xyy[k] = -4.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xyz[k] = 2.0 * to_yzz_xy[k] - 4.0 * to_yzz_xyzz[k] * tke_0 - 2.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xzz[k] = 4.0 * to_yzz_xz[k] - 4.0 * to_yzz_xzzz[k] * tke_0 - 4.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yyy[k] = -4.0 * to_yzz_yyyz[k] * tke_0 + 4.0 * to_yyyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yyz[k] = 2.0 * to_yzz_yy[k] - 4.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yzz[k] = 4.0 * to_yzz_yz[k] - 4.0 * to_yzz_yzzz[k] * tke_0 - 4.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_zzz[k] = 6.0 * to_yzz_zz[k] - 4.0 * to_yzz_zzzz[k] * tke_0 - 6.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 880-890 components of targeted buffer : GF

        auto to_y_z_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 130);

        auto to_y_z_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 131);

        auto to_y_z_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 132);

        auto to_y_z_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 133);

        auto to_y_z_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 134);

        auto to_y_z_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 135);

        auto to_y_z_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 136);

        auto to_y_z_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 137);

        auto to_y_z_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 138);

        auto to_y_z_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_y_z_yzzz_xxx, to_y_z_yzzz_xxy, to_y_z_yzzz_xxz, to_y_z_yzzz_xyy, to_y_z_yzzz_xyz, to_y_z_yzzz_xzz, to_y_z_yzzz_yyy, to_y_z_yzzz_yyz, to_y_z_yzzz_yzz, to_y_z_yzzz_zzz, to_yyzzz_xx, to_yyzzz_xxxz, to_yyzzz_xxyz, to_yyzzz_xxzz, to_yyzzz_xy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_xzzz, to_yyzzz_yy, to_yyzzz_yyyz, to_yyzzz_yyzz, to_yyzzz_yz, to_yyzzz_yzzz, to_yyzzz_zz, to_yyzzz_zzzz, to_zzz_xx, to_zzz_xxxz, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzzz_xxx[k] = -2.0 * to_zzz_xxxz[k] * tke_0 + 4.0 * to_yyzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxy[k] = -2.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xxz[k] = to_zzz_xx[k] - 2.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xyy[k] = -2.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xyz[k] = to_zzz_xy[k] - 2.0 * to_zzz_xyzz[k] * tke_0 - 2.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xzz[k] = 2.0 * to_zzz_xz[k] - 2.0 * to_zzz_xzzz[k] * tke_0 - 4.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yyy[k] = -2.0 * to_zzz_yyyz[k] * tke_0 + 4.0 * to_yyzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yyz[k] = to_zzz_yy[k] - 2.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yzz[k] = 2.0 * to_zzz_yz[k] - 2.0 * to_zzz_yzzz[k] * tke_0 - 4.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_zzz[k] = 3.0 * to_zzz_zz[k] - 2.0 * to_zzz_zzzz[k] * tke_0 - 6.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 890-900 components of targeted buffer : GF

        auto to_y_z_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 140);

        auto to_y_z_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 141);

        auto to_y_z_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 142);

        auto to_y_z_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 143);

        auto to_y_z_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 144);

        auto to_y_z_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 145);

        auto to_y_z_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 146);

        auto to_y_z_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 147);

        auto to_y_z_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 148);

        auto to_y_z_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 5 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_z_zzzz_xxx, to_y_z_zzzz_xxy, to_y_z_zzzz_xxz, to_y_z_zzzz_xyy, to_y_z_zzzz_xyz, to_y_z_zzzz_xzz, to_y_z_zzzz_yyy, to_y_z_zzzz_yyz, to_y_z_zzzz_yzz, to_y_z_zzzz_zzz, to_yzzzz_xx, to_yzzzz_xxxz, to_yzzzz_xxyz, to_yzzzz_xxzz, to_yzzzz_xy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_xzzz, to_yzzzz_yy, to_yzzzz_yyyz, to_yzzzz_yyzz, to_yzzzz_yz, to_yzzzz_yzzz, to_yzzzz_zz, to_yzzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzzz_xxx[k] = 4.0 * to_yzzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxy[k] = 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xxz[k] = -2.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xyy[k] = 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xyz[k] = -2.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xzz[k] = -4.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yyy[k] = 4.0 * to_yzzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yyz[k] = -2.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yzz[k] = -4.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_zzz[k] = -6.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 900-910 components of targeted buffer : GF

        auto to_z_x_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 0);

        auto to_z_x_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 1);

        auto to_z_x_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 2);

        auto to_z_x_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 3);

        auto to_z_x_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 4);

        auto to_z_x_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 5);

        auto to_z_x_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 6);

        auto to_z_x_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 7);

        auto to_z_x_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 8);

        auto to_z_x_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxz_xx, to_xxxxz_xxxx, to_xxxxz_xxxy, to_xxxxz_xxxz, to_xxxxz_xxyy, to_xxxxz_xxyz, to_xxxxz_xxzz, to_xxxxz_xy, to_xxxxz_xyyy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_xzzz, to_xxxxz_yy, to_xxxxz_yz, to_xxxxz_zz, to_z_x_xxxx_xxx, to_z_x_xxxx_xxy, to_z_x_xxxx_xxz, to_z_x_xxxx_xyy, to_z_x_xxxx_xyz, to_z_x_xxxx_xzz, to_z_x_xxxx_yyy, to_z_x_xxxx_yyz, to_z_x_xxxx_yzz, to_z_x_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxx_xxx[k] = -6.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxy[k] = -4.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xxz[k] = -4.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xyy[k] = -2.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xyz[k] = -2.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xzz[k] = -2.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yyy[k] = 4.0 * to_xxxxz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yyz[k] = 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yzz[k] = 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_zzz[k] = 4.0 * to_xxxxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 910-920 components of targeted buffer : GF

        auto to_z_x_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 10);

        auto to_z_x_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 11);

        auto to_z_x_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 12);

        auto to_z_x_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 13);

        auto to_z_x_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 14);

        auto to_z_x_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 15);

        auto to_z_x_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 16);

        auto to_z_x_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 17);

        auto to_z_x_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 18);

        auto to_z_x_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxx, to_xxxyz_xxxy, to_xxxyz_xxxz, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yz, to_xxxyz_zz, to_z_x_xxxy_xxx, to_z_x_xxxy_xxy, to_z_x_xxxy_xxz, to_z_x_xxxy_xyy, to_z_x_xxxy_xyz, to_z_x_xxxy_xzz, to_z_x_xxxy_yyy, to_z_x_xxxy_yyz, to_z_x_xxxy_yzz, to_z_x_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxy_xxx[k] = -6.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxy[k] = -4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xxz[k] = -4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xyy[k] = -2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xyz[k] = -2.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xzz[k] = -2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yyy[k] = 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yyz[k] = 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yzz[k] = 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_zzz[k] = 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 920-930 components of targeted buffer : GF

        auto to_z_x_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 20);

        auto to_z_x_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 21);

        auto to_z_x_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 22);

        auto to_z_x_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 23);

        auto to_z_x_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 24);

        auto to_z_x_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 25);

        auto to_z_x_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 26);

        auto to_z_x_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 27);

        auto to_z_x_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 28);

        auto to_z_x_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxx, to_xxx_xxxy, to_xxx_xxxz, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yz, to_xxx_zz, to_xxxzz_xx, to_xxxzz_xxxx, to_xxxzz_xxxy, to_xxxzz_xxxz, to_xxxzz_xxyy, to_xxxzz_xxyz, to_xxxzz_xxzz, to_xxxzz_xy, to_xxxzz_xyyy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_xzzz, to_xxxzz_yy, to_xxxzz_yz, to_xxxzz_zz, to_z_x_xxxz_xxx, to_z_x_xxxz_xxy, to_z_x_xxxz_xxz, to_z_x_xxxz_xyy, to_z_x_xxxz_xyz, to_z_x_xxxz_xzz, to_z_x_xxxz_yyy, to_z_x_xxxz_yyz, to_z_x_xxxz_yzz, to_z_x_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxz_xxx[k] = 3.0 * to_xxx_xx[k] - 2.0 * to_xxx_xxxx[k] * tke_0 - 6.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxy[k] = 2.0 * to_xxx_xy[k] - 2.0 * to_xxx_xxxy[k] * tke_0 - 4.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xxz[k] = 2.0 * to_xxx_xz[k] - 2.0 * to_xxx_xxxz[k] * tke_0 - 4.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xyy[k] = to_xxx_yy[k] - 2.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xyz[k] = to_xxx_yz[k] - 2.0 * to_xxx_xxyz[k] * tke_0 - 2.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xzz[k] = to_xxx_zz[k] - 2.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yyy[k] = -2.0 * to_xxx_xyyy[k] * tke_0 + 4.0 * to_xxxzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yyz[k] = -2.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yzz[k] = -2.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_zzz[k] = -2.0 * to_xxx_xzzz[k] * tke_0 + 4.0 * to_xxxzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 930-940 components of targeted buffer : GF

        auto to_z_x_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 30);

        auto to_z_x_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 31);

        auto to_z_x_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 32);

        auto to_z_x_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 33);

        auto to_z_x_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 34);

        auto to_z_x_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 35);

        auto to_z_x_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 36);

        auto to_z_x_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 37);

        auto to_z_x_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 38);

        auto to_z_x_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxx, to_xxyyz_xxxy, to_xxyyz_xxxz, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yz, to_xxyyz_zz, to_z_x_xxyy_xxx, to_z_x_xxyy_xxy, to_z_x_xxyy_xxz, to_z_x_xxyy_xyy, to_z_x_xxyy_xyz, to_z_x_xxyy_xzz, to_z_x_xxyy_yyy, to_z_x_xxyy_yyz, to_z_x_xxyy_yzz, to_z_x_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyy_xxx[k] = -6.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxy[k] = -4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xxz[k] = -4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xyy[k] = -2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xyz[k] = -2.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xzz[k] = -2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yyy[k] = 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yyz[k] = 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yzz[k] = 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_zzz[k] = 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 940-950 components of targeted buffer : GF

        auto to_z_x_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 40);

        auto to_z_x_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 41);

        auto to_z_x_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 42);

        auto to_z_x_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 43);

        auto to_z_x_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 44);

        auto to_z_x_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 45);

        auto to_z_x_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 46);

        auto to_z_x_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 47);

        auto to_z_x_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 48);

        auto to_z_x_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, to_xxyzz_xx, to_xxyzz_xxxx, to_xxyzz_xxxy, to_xxyzz_xxxz, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yz, to_xxyzz_zz, to_z_x_xxyz_xxx, to_z_x_xxyz_xxy, to_z_x_xxyz_xxz, to_z_x_xxyz_xyy, to_z_x_xxyz_xyz, to_z_x_xxyz_xzz, to_z_x_xxyz_yyy, to_z_x_xxyz_yyz, to_z_x_xxyz_yzz, to_z_x_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyz_xxx[k] = 3.0 * to_xxy_xx[k] - 2.0 * to_xxy_xxxx[k] * tke_0 - 6.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxy[k] = 2.0 * to_xxy_xy[k] - 2.0 * to_xxy_xxxy[k] * tke_0 - 4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xxz[k] = 2.0 * to_xxy_xz[k] - 2.0 * to_xxy_xxxz[k] * tke_0 - 4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xyy[k] = to_xxy_yy[k] - 2.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xyz[k] = to_xxy_yz[k] - 2.0 * to_xxy_xxyz[k] * tke_0 - 2.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xzz[k] = to_xxy_zz[k] - 2.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yyy[k] = -2.0 * to_xxy_xyyy[k] * tke_0 + 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yyz[k] = -2.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yzz[k] = -2.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_zzz[k] = -2.0 * to_xxy_xzzz[k] * tke_0 + 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 950-960 components of targeted buffer : GF

        auto to_z_x_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 50);

        auto to_z_x_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 51);

        auto to_z_x_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 52);

        auto to_z_x_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 53);

        auto to_z_x_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 54);

        auto to_z_x_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 55);

        auto to_z_x_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 56);

        auto to_z_x_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 57);

        auto to_z_x_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 58);

        auto to_z_x_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, to_xxzzz_xx, to_xxzzz_xxxx, to_xxzzz_xxxy, to_xxzzz_xxxz, to_xxzzz_xxyy, to_xxzzz_xxyz, to_xxzzz_xxzz, to_xxzzz_xy, to_xxzzz_xyyy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_xzzz, to_xxzzz_yy, to_xxzzz_yz, to_xxzzz_zz, to_z_x_xxzz_xxx, to_z_x_xxzz_xxy, to_z_x_xxzz_xxz, to_z_x_xxzz_xyy, to_z_x_xxzz_xyz, to_z_x_xxzz_xzz, to_z_x_xxzz_yyy, to_z_x_xxzz_yyz, to_z_x_xxzz_yzz, to_z_x_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxzz_xxx[k] = 6.0 * to_xxz_xx[k] - 4.0 * to_xxz_xxxx[k] * tke_0 - 6.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxy[k] = 4.0 * to_xxz_xy[k] - 4.0 * to_xxz_xxxy[k] * tke_0 - 4.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xxz[k] = 4.0 * to_xxz_xz[k] - 4.0 * to_xxz_xxxz[k] * tke_0 - 4.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xyy[k] = 2.0 * to_xxz_yy[k] - 4.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xyz[k] = 2.0 * to_xxz_yz[k] - 4.0 * to_xxz_xxyz[k] * tke_0 - 2.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xzz[k] = 2.0 * to_xxz_zz[k] - 4.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yyy[k] = -4.0 * to_xxz_xyyy[k] * tke_0 + 4.0 * to_xxzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yyz[k] = -4.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yzz[k] = -4.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_zzz[k] = -4.0 * to_xxz_xzzz[k] * tke_0 + 4.0 * to_xxzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 960-970 components of targeted buffer : GF

        auto to_z_x_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 60);

        auto to_z_x_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 61);

        auto to_z_x_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 62);

        auto to_z_x_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 63);

        auto to_z_x_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 64);

        auto to_z_x_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 65);

        auto to_z_x_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 66);

        auto to_z_x_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 67);

        auto to_z_x_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 68);

        auto to_z_x_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxx, to_xyyyz_xxxy, to_xyyyz_xxxz, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yz, to_xyyyz_zz, to_z_x_xyyy_xxx, to_z_x_xyyy_xxy, to_z_x_xyyy_xxz, to_z_x_xyyy_xyy, to_z_x_xyyy_xyz, to_z_x_xyyy_xzz, to_z_x_xyyy_yyy, to_z_x_xyyy_yyz, to_z_x_xyyy_yzz, to_z_x_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyy_xxx[k] = -6.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxy[k] = -4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xxz[k] = -4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xyy[k] = -2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xyz[k] = -2.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xzz[k] = -2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yyy[k] = 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yyz[k] = 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yzz[k] = 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_zzz[k] = 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 970-980 components of targeted buffer : GF

        auto to_z_x_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 70);

        auto to_z_x_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 71);

        auto to_z_x_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 72);

        auto to_z_x_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 73);

        auto to_z_x_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 74);

        auto to_z_x_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 75);

        auto to_z_x_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 76);

        auto to_z_x_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 77);

        auto to_z_x_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 78);

        auto to_z_x_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, to_xyyzz_xx, to_xyyzz_xxxx, to_xyyzz_xxxy, to_xyyzz_xxxz, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yz, to_xyyzz_zz, to_z_x_xyyz_xxx, to_z_x_xyyz_xxy, to_z_x_xyyz_xxz, to_z_x_xyyz_xyy, to_z_x_xyyz_xyz, to_z_x_xyyz_xzz, to_z_x_xyyz_yyy, to_z_x_xyyz_yyz, to_z_x_xyyz_yzz, to_z_x_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyz_xxx[k] = 3.0 * to_xyy_xx[k] - 2.0 * to_xyy_xxxx[k] * tke_0 - 6.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxy[k] = 2.0 * to_xyy_xy[k] - 2.0 * to_xyy_xxxy[k] * tke_0 - 4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xxz[k] = 2.0 * to_xyy_xz[k] - 2.0 * to_xyy_xxxz[k] * tke_0 - 4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xyy[k] = to_xyy_yy[k] - 2.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xyz[k] = to_xyy_yz[k] - 2.0 * to_xyy_xxyz[k] * tke_0 - 2.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xzz[k] = to_xyy_zz[k] - 2.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yyy[k] = -2.0 * to_xyy_xyyy[k] * tke_0 + 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yyz[k] = -2.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yzz[k] = -2.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_zzz[k] = -2.0 * to_xyy_xzzz[k] * tke_0 + 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 980-990 components of targeted buffer : GF

        auto to_z_x_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 80);

        auto to_z_x_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 81);

        auto to_z_x_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 82);

        auto to_z_x_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 83);

        auto to_z_x_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 84);

        auto to_z_x_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 85);

        auto to_z_x_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 86);

        auto to_z_x_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 87);

        auto to_z_x_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 88);

        auto to_z_x_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, to_xyzzz_xx, to_xyzzz_xxxx, to_xyzzz_xxxy, to_xyzzz_xxxz, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yz, to_xyzzz_zz, to_z_x_xyzz_xxx, to_z_x_xyzz_xxy, to_z_x_xyzz_xxz, to_z_x_xyzz_xyy, to_z_x_xyzz_xyz, to_z_x_xyzz_xzz, to_z_x_xyzz_yyy, to_z_x_xyzz_yyz, to_z_x_xyzz_yzz, to_z_x_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyzz_xxx[k] = 6.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxxx[k] * tke_0 - 6.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xxxy[k] * tke_0 - 4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xxz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xxxz[k] * tke_0 - 4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xyy[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xyz[k] = 2.0 * to_xyz_yz[k] - 4.0 * to_xyz_xxyz[k] * tke_0 - 2.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yyy[k] = -4.0 * to_xyz_xyyy[k] * tke_0 + 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yyz[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_zzz[k] = -4.0 * to_xyz_xzzz[k] * tke_0 + 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 990-1000 components of targeted buffer : GF

        auto to_z_x_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 90);

        auto to_z_x_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 91);

        auto to_z_x_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 92);

        auto to_z_x_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 93);

        auto to_z_x_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 94);

        auto to_z_x_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 95);

        auto to_z_x_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 96);

        auto to_z_x_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 97);

        auto to_z_x_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 98);

        auto to_z_x_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, to_xzzzz_xx, to_xzzzz_xxxx, to_xzzzz_xxxy, to_xzzzz_xxxz, to_xzzzz_xxyy, to_xzzzz_xxyz, to_xzzzz_xxzz, to_xzzzz_xy, to_xzzzz_xyyy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_xzzz, to_xzzzz_yy, to_xzzzz_yz, to_xzzzz_zz, to_z_x_xzzz_xxx, to_z_x_xzzz_xxy, to_z_x_xzzz_xxz, to_z_x_xzzz_xyy, to_z_x_xzzz_xyz, to_z_x_xzzz_xzz, to_z_x_xzzz_yyy, to_z_x_xzzz_yyz, to_z_x_xzzz_yzz, to_z_x_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzzz_xxx[k] = 9.0 * to_xzz_xx[k] - 6.0 * to_xzz_xxxx[k] * tke_0 - 6.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxy[k] = 6.0 * to_xzz_xy[k] - 6.0 * to_xzz_xxxy[k] * tke_0 - 4.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xxz[k] = 6.0 * to_xzz_xz[k] - 6.0 * to_xzz_xxxz[k] * tke_0 - 4.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xyy[k] = 3.0 * to_xzz_yy[k] - 6.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xyz[k] = 3.0 * to_xzz_yz[k] - 6.0 * to_xzz_xxyz[k] * tke_0 - 2.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xzz[k] = 3.0 * to_xzz_zz[k] - 6.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yyy[k] = -6.0 * to_xzz_xyyy[k] * tke_0 + 4.0 * to_xzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yyz[k] = -6.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yzz[k] = -6.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_zzz[k] = -6.0 * to_xzz_xzzz[k] * tke_0 + 4.0 * to_xzzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1000-1010 components of targeted buffer : GF

        auto to_z_x_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 100);

        auto to_z_x_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 101);

        auto to_z_x_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 102);

        auto to_z_x_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 103);

        auto to_z_x_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 104);

        auto to_z_x_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 105);

        auto to_z_x_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 106);

        auto to_z_x_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 107);

        auto to_z_x_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 108);

        auto to_z_x_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_yyyyz_xx, to_yyyyz_xxxx, to_yyyyz_xxxy, to_yyyyz_xxxz, to_yyyyz_xxyy, to_yyyyz_xxyz, to_yyyyz_xxzz, to_yyyyz_xy, to_yyyyz_xyyy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_xzzz, to_yyyyz_yy, to_yyyyz_yz, to_yyyyz_zz, to_z_x_yyyy_xxx, to_z_x_yyyy_xxy, to_z_x_yyyy_xxz, to_z_x_yyyy_xyy, to_z_x_yyyy_xyz, to_z_x_yyyy_xzz, to_z_x_yyyy_yyy, to_z_x_yyyy_yyz, to_z_x_yyyy_yzz, to_z_x_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyy_xxx[k] = -6.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxy[k] = -4.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xxz[k] = -4.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xyy[k] = -2.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xyz[k] = -2.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xzz[k] = -2.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yyy[k] = 4.0 * to_yyyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yyz[k] = 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yzz[k] = 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_zzz[k] = 4.0 * to_yyyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1010-1020 components of targeted buffer : GF

        auto to_z_x_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 110);

        auto to_z_x_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 111);

        auto to_z_x_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 112);

        auto to_z_x_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 113);

        auto to_z_x_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 114);

        auto to_z_x_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 115);

        auto to_z_x_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 116);

        auto to_z_x_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 117);

        auto to_z_x_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 118);

        auto to_z_x_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yyy_xx, to_yyy_xxxx, to_yyy_xxxy, to_yyy_xxxz, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yz, to_yyy_zz, to_yyyzz_xx, to_yyyzz_xxxx, to_yyyzz_xxxy, to_yyyzz_xxxz, to_yyyzz_xxyy, to_yyyzz_xxyz, to_yyyzz_xxzz, to_yyyzz_xy, to_yyyzz_xyyy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_xzzz, to_yyyzz_yy, to_yyyzz_yz, to_yyyzz_zz, to_z_x_yyyz_xxx, to_z_x_yyyz_xxy, to_z_x_yyyz_xxz, to_z_x_yyyz_xyy, to_z_x_yyyz_xyz, to_z_x_yyyz_xzz, to_z_x_yyyz_yyy, to_z_x_yyyz_yyz, to_z_x_yyyz_yzz, to_z_x_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyz_xxx[k] = 3.0 * to_yyy_xx[k] - 2.0 * to_yyy_xxxx[k] * tke_0 - 6.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxy[k] = 2.0 * to_yyy_xy[k] - 2.0 * to_yyy_xxxy[k] * tke_0 - 4.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xxz[k] = 2.0 * to_yyy_xz[k] - 2.0 * to_yyy_xxxz[k] * tke_0 - 4.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xyy[k] = to_yyy_yy[k] - 2.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xyz[k] = to_yyy_yz[k] - 2.0 * to_yyy_xxyz[k] * tke_0 - 2.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xzz[k] = to_yyy_zz[k] - 2.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yyy[k] = -2.0 * to_yyy_xyyy[k] * tke_0 + 4.0 * to_yyyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yyz[k] = -2.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yzz[k] = -2.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_zzz[k] = -2.0 * to_yyy_xzzz[k] * tke_0 + 4.0 * to_yyyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1020-1030 components of targeted buffer : GF

        auto to_z_x_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 120);

        auto to_z_x_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 121);

        auto to_z_x_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 122);

        auto to_z_x_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 123);

        auto to_z_x_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 124);

        auto to_z_x_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 125);

        auto to_z_x_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 126);

        auto to_z_x_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 127);

        auto to_z_x_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 128);

        auto to_z_x_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, to_yyzzz_xx, to_yyzzz_xxxx, to_yyzzz_xxxy, to_yyzzz_xxxz, to_yyzzz_xxyy, to_yyzzz_xxyz, to_yyzzz_xxzz, to_yyzzz_xy, to_yyzzz_xyyy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_xzzz, to_yyzzz_yy, to_yyzzz_yz, to_yyzzz_zz, to_z_x_yyzz_xxx, to_z_x_yyzz_xxy, to_z_x_yyzz_xxz, to_z_x_yyzz_xyy, to_z_x_yyzz_xyz, to_z_x_yyzz_xzz, to_z_x_yyzz_yyy, to_z_x_yyzz_yyz, to_z_x_yyzz_yzz, to_z_x_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyzz_xxx[k] = 6.0 * to_yyz_xx[k] - 4.0 * to_yyz_xxxx[k] * tke_0 - 6.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxy[k] = 4.0 * to_yyz_xy[k] - 4.0 * to_yyz_xxxy[k] * tke_0 - 4.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xxz[k] = 4.0 * to_yyz_xz[k] - 4.0 * to_yyz_xxxz[k] * tke_0 - 4.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xyy[k] = 2.0 * to_yyz_yy[k] - 4.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xyz[k] = 2.0 * to_yyz_yz[k] - 4.0 * to_yyz_xxyz[k] * tke_0 - 2.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xzz[k] = 2.0 * to_yyz_zz[k] - 4.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yyy[k] = -4.0 * to_yyz_xyyy[k] * tke_0 + 4.0 * to_yyzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yyz[k] = -4.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yzz[k] = -4.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_zzz[k] = -4.0 * to_yyz_xzzz[k] * tke_0 + 4.0 * to_yyzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1030-1040 components of targeted buffer : GF

        auto to_z_x_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 130);

        auto to_z_x_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 131);

        auto to_z_x_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 132);

        auto to_z_x_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 133);

        auto to_z_x_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 134);

        auto to_z_x_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 135);

        auto to_z_x_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 136);

        auto to_z_x_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 137);

        auto to_z_x_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 138);

        auto to_z_x_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, to_yzzzz_xx, to_yzzzz_xxxx, to_yzzzz_xxxy, to_yzzzz_xxxz, to_yzzzz_xxyy, to_yzzzz_xxyz, to_yzzzz_xxzz, to_yzzzz_xy, to_yzzzz_xyyy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_xzzz, to_yzzzz_yy, to_yzzzz_yz, to_yzzzz_zz, to_z_x_yzzz_xxx, to_z_x_yzzz_xxy, to_z_x_yzzz_xxz, to_z_x_yzzz_xyy, to_z_x_yzzz_xyz, to_z_x_yzzz_xzz, to_z_x_yzzz_yyy, to_z_x_yzzz_yyz, to_z_x_yzzz_yzz, to_z_x_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzzz_xxx[k] = 9.0 * to_yzz_xx[k] - 6.0 * to_yzz_xxxx[k] * tke_0 - 6.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxy[k] = 6.0 * to_yzz_xy[k] - 6.0 * to_yzz_xxxy[k] * tke_0 - 4.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xxz[k] = 6.0 * to_yzz_xz[k] - 6.0 * to_yzz_xxxz[k] * tke_0 - 4.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xyy[k] = 3.0 * to_yzz_yy[k] - 6.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xyz[k] = 3.0 * to_yzz_yz[k] - 6.0 * to_yzz_xxyz[k] * tke_0 - 2.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xzz[k] = 3.0 * to_yzz_zz[k] - 6.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yyy[k] = -6.0 * to_yzz_xyyy[k] * tke_0 + 4.0 * to_yzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yyz[k] = -6.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yzz[k] = -6.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_zzz[k] = -6.0 * to_yzz_xzzz[k] * tke_0 + 4.0 * to_yzzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1040-1050 components of targeted buffer : GF

        auto to_z_x_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 140);

        auto to_z_x_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 141);

        auto to_z_x_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 142);

        auto to_z_x_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 143);

        auto to_z_x_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 144);

        auto to_z_x_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 145);

        auto to_z_x_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 146);

        auto to_z_x_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 147);

        auto to_z_x_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 148);

        auto to_z_x_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 6 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_x_zzzz_xxx, to_z_x_zzzz_xxy, to_z_x_zzzz_xxz, to_z_x_zzzz_xyy, to_z_x_zzzz_xyz, to_z_x_zzzz_xzz, to_z_x_zzzz_yyy, to_z_x_zzzz_yyz, to_z_x_zzzz_yzz, to_z_x_zzzz_zzz, to_zzz_xx, to_zzz_xxxx, to_zzz_xxxy, to_zzz_xxxz, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yz, to_zzz_zz, to_zzzzz_xx, to_zzzzz_xxxx, to_zzzzz_xxxy, to_zzzzz_xxxz, to_zzzzz_xxyy, to_zzzzz_xxyz, to_zzzzz_xxzz, to_zzzzz_xy, to_zzzzz_xyyy, to_zzzzz_xyyz, to_zzzzz_xyzz, to_zzzzz_xz, to_zzzzz_xzzz, to_zzzzz_yy, to_zzzzz_yz, to_zzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzzz_xxx[k] = 12.0 * to_zzz_xx[k] - 8.0 * to_zzz_xxxx[k] * tke_0 - 6.0 * to_zzzzz_xx[k] * tbe_0 + 4.0 * to_zzzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxy[k] = 8.0 * to_zzz_xy[k] - 8.0 * to_zzz_xxxy[k] * tke_0 - 4.0 * to_zzzzz_xy[k] * tbe_0 + 4.0 * to_zzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xxz[k] = 8.0 * to_zzz_xz[k] - 8.0 * to_zzz_xxxz[k] * tke_0 - 4.0 * to_zzzzz_xz[k] * tbe_0 + 4.0 * to_zzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xyy[k] = 4.0 * to_zzz_yy[k] - 8.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_zzzzz_yy[k] * tbe_0 + 4.0 * to_zzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xyz[k] = 4.0 * to_zzz_yz[k] - 8.0 * to_zzz_xxyz[k] * tke_0 - 2.0 * to_zzzzz_yz[k] * tbe_0 + 4.0 * to_zzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xzz[k] = 4.0 * to_zzz_zz[k] - 8.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_zzzzz_zz[k] * tbe_0 + 4.0 * to_zzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yyy[k] = -8.0 * to_zzz_xyyy[k] * tke_0 + 4.0 * to_zzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yyz[k] = -8.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_zzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yzz[k] = -8.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_zzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_zzz[k] = -8.0 * to_zzz_xzzz[k] * tke_0 + 4.0 * to_zzzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1050-1060 components of targeted buffer : GF

        auto to_z_y_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 0);

        auto to_z_y_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 1);

        auto to_z_y_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 2);

        auto to_z_y_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 3);

        auto to_z_y_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 4);

        auto to_z_y_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 5);

        auto to_z_y_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 6);

        auto to_z_y_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 7);

        auto to_z_y_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 8);

        auto to_z_y_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxz_xx, to_xxxxz_xxxy, to_xxxxz_xxyy, to_xxxxz_xxyz, to_xxxxz_xy, to_xxxxz_xyyy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_yy, to_xxxxz_yyyy, to_xxxxz_yyyz, to_xxxxz_yyzz, to_xxxxz_yz, to_xxxxz_yzzz, to_xxxxz_zz, to_z_y_xxxx_xxx, to_z_y_xxxx_xxy, to_z_y_xxxx_xxz, to_z_y_xxxx_xyy, to_z_y_xxxx_xyz, to_z_y_xxxx_xzz, to_z_y_xxxx_yyy, to_z_y_xxxx_yyz, to_z_y_xxxx_yzz, to_z_y_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxx_xxx[k] = 4.0 * to_xxxxz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxy[k] = -2.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xxz[k] = 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xyy[k] = -4.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xyz[k] = -2.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xzz[k] = 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yyy[k] = -6.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yyz[k] = -4.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yzz[k] = -2.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_zzz[k] = 4.0 * to_xxxxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1060-1070 components of targeted buffer : GF

        auto to_z_y_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 10);

        auto to_z_y_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 11);

        auto to_z_y_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 12);

        auto to_z_y_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 13);

        auto to_z_y_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 14);

        auto to_z_y_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 15);

        auto to_z_y_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 16);

        auto to_z_y_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 17);

        auto to_z_y_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 18);

        auto to_z_y_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxy, to_xxxyz_xxyy, to_xxxyz_xxyz, to_xxxyz_xy, to_xxxyz_xyyy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_yy, to_xxxyz_yyyy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_z_y_xxxy_xxx, to_z_y_xxxy_xxy, to_z_y_xxxy_xxz, to_z_y_xxxy_xyy, to_z_y_xxxy_xyz, to_z_y_xxxy_xzz, to_z_y_xxxy_yyy, to_z_y_xxxy_yyz, to_z_y_xxxy_yzz, to_z_y_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxy_xxx[k] = 4.0 * to_xxxyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxy[k] = -2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xxz[k] = 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xyy[k] = -4.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xyz[k] = -2.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xzz[k] = 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yyy[k] = -6.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yyz[k] = -4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yzz[k] = -2.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_zzz[k] = 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1070-1080 components of targeted buffer : GF

        auto to_z_y_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 20);

        auto to_z_y_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 21);

        auto to_z_y_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 22);

        auto to_z_y_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 23);

        auto to_z_y_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 24);

        auto to_z_y_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 25);

        auto to_z_y_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 26);

        auto to_z_y_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 27);

        auto to_z_y_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 28);

        auto to_z_y_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxy, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_yy, to_xxx_yyyy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxxzz_xx, to_xxxzz_xxxy, to_xxxzz_xxyy, to_xxxzz_xxyz, to_xxxzz_xy, to_xxxzz_xyyy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_yy, to_xxxzz_yyyy, to_xxxzz_yyyz, to_xxxzz_yyzz, to_xxxzz_yz, to_xxxzz_yzzz, to_xxxzz_zz, to_z_y_xxxz_xxx, to_z_y_xxxz_xxy, to_z_y_xxxz_xxz, to_z_y_xxxz_xyy, to_z_y_xxxz_xyz, to_z_y_xxxz_xzz, to_z_y_xxxz_yyy, to_z_y_xxxz_yyz, to_z_y_xxxz_yzz, to_z_y_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxz_xxx[k] = -2.0 * to_xxx_xxxy[k] * tke_0 + 4.0 * to_xxxzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxy[k] = to_xxx_xx[k] - 2.0 * to_xxx_xxyy[k] * tke_0 - 2.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xxz[k] = -2.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xyy[k] = 2.0 * to_xxx_xy[k] - 2.0 * to_xxx_xyyy[k] * tke_0 - 4.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xyz[k] = to_xxx_xz[k] - 2.0 * to_xxx_xyyz[k] * tke_0 - 2.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xzz[k] = -2.0 * to_xxx_xyzz[k] * tke_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yyy[k] = 3.0 * to_xxx_yy[k] - 2.0 * to_xxx_yyyy[k] * tke_0 - 6.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yyz[k] = 2.0 * to_xxx_yz[k] - 2.0 * to_xxx_yyyz[k] * tke_0 - 4.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yzz[k] = to_xxx_zz[k] - 2.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_zzz[k] = -2.0 * to_xxx_yzzz[k] * tke_0 + 4.0 * to_xxxzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1080-1090 components of targeted buffer : GF

        auto to_z_y_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 30);

        auto to_z_y_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 31);

        auto to_z_y_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 32);

        auto to_z_y_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 33);

        auto to_z_y_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 34);

        auto to_z_y_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 35);

        auto to_z_y_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 36);

        auto to_z_y_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 37);

        auto to_z_y_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 38);

        auto to_z_y_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxy, to_xxyyz_xxyy, to_xxyyz_xxyz, to_xxyyz_xy, to_xxyyz_xyyy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_yy, to_xxyyz_yyyy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_z_y_xxyy_xxx, to_z_y_xxyy_xxy, to_z_y_xxyy_xxz, to_z_y_xxyy_xyy, to_z_y_xxyy_xyz, to_z_y_xxyy_xzz, to_z_y_xxyy_yyy, to_z_y_xxyy_yyz, to_z_y_xxyy_yzz, to_z_y_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyy_xxx[k] = 4.0 * to_xxyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxy[k] = -2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xxz[k] = 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xyy[k] = -4.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xyz[k] = -2.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xzz[k] = 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yyy[k] = -6.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yyz[k] = -4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yzz[k] = -2.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_zzz[k] = 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1090-1100 components of targeted buffer : GF

        auto to_z_y_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 40);

        auto to_z_y_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 41);

        auto to_z_y_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 42);

        auto to_z_y_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 43);

        auto to_z_y_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 44);

        auto to_z_y_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 45);

        auto to_z_y_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 46);

        auto to_z_y_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 47);

        auto to_z_y_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 48);

        auto to_z_y_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxyzz_xx, to_xxyzz_xxxy, to_xxyzz_xxyy, to_xxyzz_xxyz, to_xxyzz_xy, to_xxyzz_xyyy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_yy, to_xxyzz_yyyy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_z_y_xxyz_xxx, to_z_y_xxyz_xxy, to_z_y_xxyz_xxz, to_z_y_xxyz_xyy, to_z_y_xxyz_xyz, to_z_y_xxyz_xzz, to_z_y_xxyz_yyy, to_z_y_xxyz_yyz, to_z_y_xxyz_yzz, to_z_y_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyz_xxx[k] = -2.0 * to_xxy_xxxy[k] * tke_0 + 4.0 * to_xxyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxy[k] = to_xxy_xx[k] - 2.0 * to_xxy_xxyy[k] * tke_0 - 2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xxz[k] = -2.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xyy[k] = 2.0 * to_xxy_xy[k] - 2.0 * to_xxy_xyyy[k] * tke_0 - 4.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xyz[k] = to_xxy_xz[k] - 2.0 * to_xxy_xyyz[k] * tke_0 - 2.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xzz[k] = -2.0 * to_xxy_xyzz[k] * tke_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yyy[k] = 3.0 * to_xxy_yy[k] - 2.0 * to_xxy_yyyy[k] * tke_0 - 6.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yyz[k] = 2.0 * to_xxy_yz[k] - 2.0 * to_xxy_yyyz[k] * tke_0 - 4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yzz[k] = to_xxy_zz[k] - 2.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_zzz[k] = -2.0 * to_xxy_yzzz[k] * tke_0 + 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1100-1110 components of targeted buffer : GF

        auto to_z_y_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 50);

        auto to_z_y_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 51);

        auto to_z_y_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 52);

        auto to_z_y_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 53);

        auto to_z_y_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 54);

        auto to_z_y_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 55);

        auto to_z_y_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 56);

        auto to_z_y_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 57);

        auto to_z_y_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 58);

        auto to_z_y_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxzzz_xx, to_xxzzz_xxxy, to_xxzzz_xxyy, to_xxzzz_xxyz, to_xxzzz_xy, to_xxzzz_xyyy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_yy, to_xxzzz_yyyy, to_xxzzz_yyyz, to_xxzzz_yyzz, to_xxzzz_yz, to_xxzzz_yzzz, to_xxzzz_zz, to_z_y_xxzz_xxx, to_z_y_xxzz_xxy, to_z_y_xxzz_xxz, to_z_y_xxzz_xyy, to_z_y_xxzz_xyz, to_z_y_xxzz_xzz, to_z_y_xxzz_yyy, to_z_y_xxzz_yyz, to_z_y_xxzz_yzz, to_z_y_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxzz_xxx[k] = -4.0 * to_xxz_xxxy[k] * tke_0 + 4.0 * to_xxzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxy[k] = 2.0 * to_xxz_xx[k] - 4.0 * to_xxz_xxyy[k] * tke_0 - 2.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xxz[k] = -4.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xyy[k] = 4.0 * to_xxz_xy[k] - 4.0 * to_xxz_xyyy[k] * tke_0 - 4.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xyz[k] = 2.0 * to_xxz_xz[k] - 4.0 * to_xxz_xyyz[k] * tke_0 - 2.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xzz[k] = -4.0 * to_xxz_xyzz[k] * tke_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yyy[k] = 6.0 * to_xxz_yy[k] - 4.0 * to_xxz_yyyy[k] * tke_0 - 6.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yyz[k] = 4.0 * to_xxz_yz[k] - 4.0 * to_xxz_yyyz[k] * tke_0 - 4.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yzz[k] = 2.0 * to_xxz_zz[k] - 4.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_zzz[k] = -4.0 * to_xxz_yzzz[k] * tke_0 + 4.0 * to_xxzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1110-1120 components of targeted buffer : GF

        auto to_z_y_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 60);

        auto to_z_y_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 61);

        auto to_z_y_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 62);

        auto to_z_y_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 63);

        auto to_z_y_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 64);

        auto to_z_y_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 65);

        auto to_z_y_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 66);

        auto to_z_y_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 67);

        auto to_z_y_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 68);

        auto to_z_y_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxy, to_xyyyz_xxyy, to_xyyyz_xxyz, to_xyyyz_xy, to_xyyyz_xyyy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_yy, to_xyyyz_yyyy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, to_z_y_xyyy_xxx, to_z_y_xyyy_xxy, to_z_y_xyyy_xxz, to_z_y_xyyy_xyy, to_z_y_xyyy_xyz, to_z_y_xyyy_xzz, to_z_y_xyyy_yyy, to_z_y_xyyy_yyz, to_z_y_xyyy_yzz, to_z_y_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyy_xxx[k] = 4.0 * to_xyyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxy[k] = -2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xxz[k] = 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xyy[k] = -4.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xyz[k] = -2.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xzz[k] = 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yyy[k] = -6.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yyz[k] = -4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yzz[k] = -2.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_zzz[k] = 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1120-1130 components of targeted buffer : GF

        auto to_z_y_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 70);

        auto to_z_y_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 71);

        auto to_z_y_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 72);

        auto to_z_y_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 73);

        auto to_z_y_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 74);

        auto to_z_y_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 75);

        auto to_z_y_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 76);

        auto to_z_y_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 77);

        auto to_z_y_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 78);

        auto to_z_y_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyyzz_xx, to_xyyzz_xxxy, to_xyyzz_xxyy, to_xyyzz_xxyz, to_xyyzz_xy, to_xyyzz_xyyy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_yy, to_xyyzz_yyyy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, to_z_y_xyyz_xxx, to_z_y_xyyz_xxy, to_z_y_xyyz_xxz, to_z_y_xyyz_xyy, to_z_y_xyyz_xyz, to_z_y_xyyz_xzz, to_z_y_xyyz_yyy, to_z_y_xyyz_yyz, to_z_y_xyyz_yzz, to_z_y_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyz_xxx[k] = -2.0 * to_xyy_xxxy[k] * tke_0 + 4.0 * to_xyyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxy[k] = to_xyy_xx[k] - 2.0 * to_xyy_xxyy[k] * tke_0 - 2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xxz[k] = -2.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xyy[k] = 2.0 * to_xyy_xy[k] - 2.0 * to_xyy_xyyy[k] * tke_0 - 4.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xyz[k] = to_xyy_xz[k] - 2.0 * to_xyy_xyyz[k] * tke_0 - 2.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xzz[k] = -2.0 * to_xyy_xyzz[k] * tke_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yyy[k] = 3.0 * to_xyy_yy[k] - 2.0 * to_xyy_yyyy[k] * tke_0 - 6.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yyz[k] = 2.0 * to_xyy_yz[k] - 2.0 * to_xyy_yyyz[k] * tke_0 - 4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yzz[k] = to_xyy_zz[k] - 2.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_zzz[k] = -2.0 * to_xyy_yzzz[k] * tke_0 + 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1130-1140 components of targeted buffer : GF

        auto to_z_y_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 80);

        auto to_z_y_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 81);

        auto to_z_y_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 82);

        auto to_z_y_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 83);

        auto to_z_y_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 84);

        auto to_z_y_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 85);

        auto to_z_y_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 86);

        auto to_z_y_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 87);

        auto to_z_y_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 88);

        auto to_z_y_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyzzz_xx, to_xyzzz_xxxy, to_xyzzz_xxyy, to_xyzzz_xxyz, to_xyzzz_xy, to_xyzzz_xyyy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_yy, to_xyzzz_yyyy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, to_z_y_xyzz_xxx, to_z_y_xyzz_xxy, to_z_y_xyzz_xxz, to_z_y_xyzz_xyy, to_z_y_xyzz_xyz, to_z_y_xyzz_xzz, to_z_y_xyzz_yyy, to_z_y_xyzz_yyz, to_z_y_xyzz_yzz, to_z_y_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyzz_xxx[k] = -4.0 * to_xyz_xxxy[k] * tke_0 + 4.0 * to_xyzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxy[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxyy[k] * tke_0 - 2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xxz[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xyy[k] = 4.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyyy[k] * tke_0 - 4.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xyz[k] = 2.0 * to_xyz_xz[k] - 4.0 * to_xyz_xyyz[k] * tke_0 - 2.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xzz[k] = -4.0 * to_xyz_xyzz[k] * tke_0 + 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yyy[k] = 6.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyyy[k] * tke_0 - 6.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yyz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yyyz[k] * tke_0 - 4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yzz[k] = 2.0 * to_xyz_zz[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_zzz[k] = -4.0 * to_xyz_yzzz[k] * tke_0 + 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1140-1150 components of targeted buffer : GF

        auto to_z_y_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 90);

        auto to_z_y_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 91);

        auto to_z_y_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 92);

        auto to_z_y_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 93);

        auto to_z_y_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 94);

        auto to_z_y_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 95);

        auto to_z_y_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 96);

        auto to_z_y_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 97);

        auto to_z_y_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 98);

        auto to_z_y_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzzzz_xx, to_xzzzz_xxxy, to_xzzzz_xxyy, to_xzzzz_xxyz, to_xzzzz_xy, to_xzzzz_xyyy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_yy, to_xzzzz_yyyy, to_xzzzz_yyyz, to_xzzzz_yyzz, to_xzzzz_yz, to_xzzzz_yzzz, to_xzzzz_zz, to_z_y_xzzz_xxx, to_z_y_xzzz_xxy, to_z_y_xzzz_xxz, to_z_y_xzzz_xyy, to_z_y_xzzz_xyz, to_z_y_xzzz_xzz, to_z_y_xzzz_yyy, to_z_y_xzzz_yyz, to_z_y_xzzz_yzz, to_z_y_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzzz_xxx[k] = -6.0 * to_xzz_xxxy[k] * tke_0 + 4.0 * to_xzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxy[k] = 3.0 * to_xzz_xx[k] - 6.0 * to_xzz_xxyy[k] * tke_0 - 2.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xxz[k] = -6.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xyy[k] = 6.0 * to_xzz_xy[k] - 6.0 * to_xzz_xyyy[k] * tke_0 - 4.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xyz[k] = 3.0 * to_xzz_xz[k] - 6.0 * to_xzz_xyyz[k] * tke_0 - 2.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xzz[k] = -6.0 * to_xzz_xyzz[k] * tke_0 + 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yyy[k] = 9.0 * to_xzz_yy[k] - 6.0 * to_xzz_yyyy[k] * tke_0 - 6.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yyz[k] = 6.0 * to_xzz_yz[k] - 6.0 * to_xzz_yyyz[k] * tke_0 - 4.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yzz[k] = 3.0 * to_xzz_zz[k] - 6.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_zzz[k] = -6.0 * to_xzz_yzzz[k] * tke_0 + 4.0 * to_xzzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1150-1160 components of targeted buffer : GF

        auto to_z_y_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 100);

        auto to_z_y_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 101);

        auto to_z_y_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 102);

        auto to_z_y_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 103);

        auto to_z_y_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 104);

        auto to_z_y_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 105);

        auto to_z_y_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 106);

        auto to_z_y_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 107);

        auto to_z_y_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 108);

        auto to_z_y_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_yyyyz_xx, to_yyyyz_xxxy, to_yyyyz_xxyy, to_yyyyz_xxyz, to_yyyyz_xy, to_yyyyz_xyyy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_yy, to_yyyyz_yyyy, to_yyyyz_yyyz, to_yyyyz_yyzz, to_yyyyz_yz, to_yyyyz_yzzz, to_yyyyz_zz, to_z_y_yyyy_xxx, to_z_y_yyyy_xxy, to_z_y_yyyy_xxz, to_z_y_yyyy_xyy, to_z_y_yyyy_xyz, to_z_y_yyyy_xzz, to_z_y_yyyy_yyy, to_z_y_yyyy_yyz, to_z_y_yyyy_yzz, to_z_y_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyy_xxx[k] = 4.0 * to_yyyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxy[k] = -2.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xxz[k] = 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xyy[k] = -4.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xyz[k] = -2.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xzz[k] = 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yyy[k] = -6.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yyz[k] = -4.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yzz[k] = -2.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_zzz[k] = 4.0 * to_yyyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1160-1170 components of targeted buffer : GF

        auto to_z_y_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 110);

        auto to_z_y_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 111);

        auto to_z_y_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 112);

        auto to_z_y_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 113);

        auto to_z_y_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 114);

        auto to_z_y_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 115);

        auto to_z_y_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 116);

        auto to_z_y_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 117);

        auto to_z_y_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 118);

        auto to_z_y_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yyy_xx, to_yyy_xxxy, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_yy, to_yyy_yyyy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyyzz_xx, to_yyyzz_xxxy, to_yyyzz_xxyy, to_yyyzz_xxyz, to_yyyzz_xy, to_yyyzz_xyyy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_yy, to_yyyzz_yyyy, to_yyyzz_yyyz, to_yyyzz_yyzz, to_yyyzz_yz, to_yyyzz_yzzz, to_yyyzz_zz, to_z_y_yyyz_xxx, to_z_y_yyyz_xxy, to_z_y_yyyz_xxz, to_z_y_yyyz_xyy, to_z_y_yyyz_xyz, to_z_y_yyyz_xzz, to_z_y_yyyz_yyy, to_z_y_yyyz_yyz, to_z_y_yyyz_yzz, to_z_y_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyz_xxx[k] = -2.0 * to_yyy_xxxy[k] * tke_0 + 4.0 * to_yyyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxy[k] = to_yyy_xx[k] - 2.0 * to_yyy_xxyy[k] * tke_0 - 2.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xxz[k] = -2.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xyy[k] = 2.0 * to_yyy_xy[k] - 2.0 * to_yyy_xyyy[k] * tke_0 - 4.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xyz[k] = to_yyy_xz[k] - 2.0 * to_yyy_xyyz[k] * tke_0 - 2.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xzz[k] = -2.0 * to_yyy_xyzz[k] * tke_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yyy[k] = 3.0 * to_yyy_yy[k] - 2.0 * to_yyy_yyyy[k] * tke_0 - 6.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yyz[k] = 2.0 * to_yyy_yz[k] - 2.0 * to_yyy_yyyz[k] * tke_0 - 4.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yzz[k] = to_yyy_zz[k] - 2.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_zzz[k] = -2.0 * to_yyy_yzzz[k] * tke_0 + 4.0 * to_yyyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1170-1180 components of targeted buffer : GF

        auto to_z_y_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 120);

        auto to_z_y_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 121);

        auto to_z_y_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 122);

        auto to_z_y_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 123);

        auto to_z_y_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 124);

        auto to_z_y_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 125);

        auto to_z_y_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 126);

        auto to_z_y_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 127);

        auto to_z_y_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 128);

        auto to_z_y_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyzzz_xx, to_yyzzz_xxxy, to_yyzzz_xxyy, to_yyzzz_xxyz, to_yyzzz_xy, to_yyzzz_xyyy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_yy, to_yyzzz_yyyy, to_yyzzz_yyyz, to_yyzzz_yyzz, to_yyzzz_yz, to_yyzzz_yzzz, to_yyzzz_zz, to_z_y_yyzz_xxx, to_z_y_yyzz_xxy, to_z_y_yyzz_xxz, to_z_y_yyzz_xyy, to_z_y_yyzz_xyz, to_z_y_yyzz_xzz, to_z_y_yyzz_yyy, to_z_y_yyzz_yyz, to_z_y_yyzz_yzz, to_z_y_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyzz_xxx[k] = -4.0 * to_yyz_xxxy[k] * tke_0 + 4.0 * to_yyzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxy[k] = 2.0 * to_yyz_xx[k] - 4.0 * to_yyz_xxyy[k] * tke_0 - 2.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xxz[k] = -4.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xyy[k] = 4.0 * to_yyz_xy[k] - 4.0 * to_yyz_xyyy[k] * tke_0 - 4.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xyz[k] = 2.0 * to_yyz_xz[k] - 4.0 * to_yyz_xyyz[k] * tke_0 - 2.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xzz[k] = -4.0 * to_yyz_xyzz[k] * tke_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yyy[k] = 6.0 * to_yyz_yy[k] - 4.0 * to_yyz_yyyy[k] * tke_0 - 6.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yyz[k] = 4.0 * to_yyz_yz[k] - 4.0 * to_yyz_yyyz[k] * tke_0 - 4.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yzz[k] = 2.0 * to_yyz_zz[k] - 4.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_zzz[k] = -4.0 * to_yyz_yzzz[k] * tke_0 + 4.0 * to_yyzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1180-1190 components of targeted buffer : GF

        auto to_z_y_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 130);

        auto to_z_y_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 131);

        auto to_z_y_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 132);

        auto to_z_y_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 133);

        auto to_z_y_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 134);

        auto to_z_y_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 135);

        auto to_z_y_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 136);

        auto to_z_y_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 137);

        auto to_z_y_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 138);

        auto to_z_y_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzzzz_xx, to_yzzzz_xxxy, to_yzzzz_xxyy, to_yzzzz_xxyz, to_yzzzz_xy, to_yzzzz_xyyy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_yy, to_yzzzz_yyyy, to_yzzzz_yyyz, to_yzzzz_yyzz, to_yzzzz_yz, to_yzzzz_yzzz, to_yzzzz_zz, to_z_y_yzzz_xxx, to_z_y_yzzz_xxy, to_z_y_yzzz_xxz, to_z_y_yzzz_xyy, to_z_y_yzzz_xyz, to_z_y_yzzz_xzz, to_z_y_yzzz_yyy, to_z_y_yzzz_yyz, to_z_y_yzzz_yzz, to_z_y_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzzz_xxx[k] = -6.0 * to_yzz_xxxy[k] * tke_0 + 4.0 * to_yzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxy[k] = 3.0 * to_yzz_xx[k] - 6.0 * to_yzz_xxyy[k] * tke_0 - 2.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xxz[k] = -6.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xyy[k] = 6.0 * to_yzz_xy[k] - 6.0 * to_yzz_xyyy[k] * tke_0 - 4.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xyz[k] = 3.0 * to_yzz_xz[k] - 6.0 * to_yzz_xyyz[k] * tke_0 - 2.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xzz[k] = -6.0 * to_yzz_xyzz[k] * tke_0 + 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yyy[k] = 9.0 * to_yzz_yy[k] - 6.0 * to_yzz_yyyy[k] * tke_0 - 6.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yyz[k] = 6.0 * to_yzz_yz[k] - 6.0 * to_yzz_yyyz[k] * tke_0 - 4.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yzz[k] = 3.0 * to_yzz_zz[k] - 6.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_zzz[k] = -6.0 * to_yzz_yzzz[k] * tke_0 + 4.0 * to_yzzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1190-1200 components of targeted buffer : GF

        auto to_z_y_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 140);

        auto to_z_y_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 141);

        auto to_z_y_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 142);

        auto to_z_y_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 143);

        auto to_z_y_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 144);

        auto to_z_y_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 145);

        auto to_z_y_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 146);

        auto to_z_y_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 147);

        auto to_z_y_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 148);

        auto to_z_y_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 7 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_y_zzzz_xxx, to_z_y_zzzz_xxy, to_z_y_zzzz_xxz, to_z_y_zzzz_xyy, to_z_y_zzzz_xyz, to_z_y_zzzz_xzz, to_z_y_zzzz_yyy, to_z_y_zzzz_yyz, to_z_y_zzzz_yzz, to_z_y_zzzz_zzz, to_zzz_xx, to_zzz_xxxy, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_yy, to_zzz_yyyy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzzzz_xx, to_zzzzz_xxxy, to_zzzzz_xxyy, to_zzzzz_xxyz, to_zzzzz_xy, to_zzzzz_xyyy, to_zzzzz_xyyz, to_zzzzz_xyzz, to_zzzzz_xz, to_zzzzz_yy, to_zzzzz_yyyy, to_zzzzz_yyyz, to_zzzzz_yyzz, to_zzzzz_yz, to_zzzzz_yzzz, to_zzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzzz_xxx[k] = -8.0 * to_zzz_xxxy[k] * tke_0 + 4.0 * to_zzzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxy[k] = 4.0 * to_zzz_xx[k] - 8.0 * to_zzz_xxyy[k] * tke_0 - 2.0 * to_zzzzz_xx[k] * tbe_0 + 4.0 * to_zzzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xxz[k] = -8.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_zzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xyy[k] = 8.0 * to_zzz_xy[k] - 8.0 * to_zzz_xyyy[k] * tke_0 - 4.0 * to_zzzzz_xy[k] * tbe_0 + 4.0 * to_zzzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xyz[k] = 4.0 * to_zzz_xz[k] - 8.0 * to_zzz_xyyz[k] * tke_0 - 2.0 * to_zzzzz_xz[k] * tbe_0 + 4.0 * to_zzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xzz[k] = -8.0 * to_zzz_xyzz[k] * tke_0 + 4.0 * to_zzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yyy[k] = 12.0 * to_zzz_yy[k] - 8.0 * to_zzz_yyyy[k] * tke_0 - 6.0 * to_zzzzz_yy[k] * tbe_0 + 4.0 * to_zzzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yyz[k] = 8.0 * to_zzz_yz[k] - 8.0 * to_zzz_yyyz[k] * tke_0 - 4.0 * to_zzzzz_yz[k] * tbe_0 + 4.0 * to_zzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yzz[k] = 4.0 * to_zzz_zz[k] - 8.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_zzzzz_zz[k] * tbe_0 + 4.0 * to_zzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_zzz[k] = -8.0 * to_zzz_yzzz[k] * tke_0 + 4.0 * to_zzzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1200-1210 components of targeted buffer : GF

        auto to_z_z_xxxx_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 0);

        auto to_z_z_xxxx_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 1);

        auto to_z_z_xxxx_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 2);

        auto to_z_z_xxxx_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 3);

        auto to_z_z_xxxx_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 4);

        auto to_z_z_xxxx_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 5);

        auto to_z_z_xxxx_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 6);

        auto to_z_z_xxxx_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 7);

        auto to_z_z_xxxx_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 8);

        auto to_z_z_xxxx_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_xxxxz_xx, to_xxxxz_xxxz, to_xxxxz_xxyz, to_xxxxz_xxzz, to_xxxxz_xy, to_xxxxz_xyyz, to_xxxxz_xyzz, to_xxxxz_xz, to_xxxxz_xzzz, to_xxxxz_yy, to_xxxxz_yyyz, to_xxxxz_yyzz, to_xxxxz_yz, to_xxxxz_yzzz, to_xxxxz_zz, to_xxxxz_zzzz, to_z_z_xxxx_xxx, to_z_z_xxxx_xxy, to_z_z_xxxx_xxz, to_z_z_xxxx_xyy, to_z_z_xxxx_xyz, to_z_z_xxxx_xzz, to_z_z_xxxx_yyy, to_z_z_xxxx_yyz, to_z_z_xxxx_yzz, to_z_z_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxx_xxx[k] = 4.0 * to_xxxxz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxy[k] = 4.0 * to_xxxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xxz[k] = -2.0 * to_xxxxz_xx[k] * tbe_0 + 4.0 * to_xxxxz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xyy[k] = 4.0 * to_xxxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xyz[k] = -2.0 * to_xxxxz_xy[k] * tbe_0 + 4.0 * to_xxxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xzz[k] = -4.0 * to_xxxxz_xz[k] * tbe_0 + 4.0 * to_xxxxz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yyy[k] = 4.0 * to_xxxxz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yyz[k] = -2.0 * to_xxxxz_yy[k] * tbe_0 + 4.0 * to_xxxxz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yzz[k] = -4.0 * to_xxxxz_yz[k] * tbe_0 + 4.0 * to_xxxxz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_zzz[k] = -6.0 * to_xxxxz_zz[k] * tbe_0 + 4.0 * to_xxxxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1210-1220 components of targeted buffer : GF

        auto to_z_z_xxxy_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 10);

        auto to_z_z_xxxy_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 11);

        auto to_z_z_xxxy_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 12);

        auto to_z_z_xxxy_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 13);

        auto to_z_z_xxxy_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 14);

        auto to_z_z_xxxy_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 15);

        auto to_z_z_xxxy_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 16);

        auto to_z_z_xxxy_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 17);

        auto to_z_z_xxxy_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 18);

        auto to_z_z_xxxy_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_xxxyz_xx, to_xxxyz_xxxz, to_xxxyz_xxyz, to_xxxyz_xxzz, to_xxxyz_xy, to_xxxyz_xyyz, to_xxxyz_xyzz, to_xxxyz_xz, to_xxxyz_xzzz, to_xxxyz_yy, to_xxxyz_yyyz, to_xxxyz_yyzz, to_xxxyz_yz, to_xxxyz_yzzz, to_xxxyz_zz, to_xxxyz_zzzz, to_z_z_xxxy_xxx, to_z_z_xxxy_xxy, to_z_z_xxxy_xxz, to_z_z_xxxy_xyy, to_z_z_xxxy_xyz, to_z_z_xxxy_xzz, to_z_z_xxxy_yyy, to_z_z_xxxy_yyz, to_z_z_xxxy_yzz, to_z_z_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxy_xxx[k] = 4.0 * to_xxxyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxy[k] = 4.0 * to_xxxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xxz[k] = -2.0 * to_xxxyz_xx[k] * tbe_0 + 4.0 * to_xxxyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xyy[k] = 4.0 * to_xxxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xyz[k] = -2.0 * to_xxxyz_xy[k] * tbe_0 + 4.0 * to_xxxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xzz[k] = -4.0 * to_xxxyz_xz[k] * tbe_0 + 4.0 * to_xxxyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yyy[k] = 4.0 * to_xxxyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yyz[k] = -2.0 * to_xxxyz_yy[k] * tbe_0 + 4.0 * to_xxxyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yzz[k] = -4.0 * to_xxxyz_yz[k] * tbe_0 + 4.0 * to_xxxyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_zzz[k] = -6.0 * to_xxxyz_zz[k] * tbe_0 + 4.0 * to_xxxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1220-1230 components of targeted buffer : GF

        auto to_z_z_xxxz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 20);

        auto to_z_z_xxxz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 21);

        auto to_z_z_xxxz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 22);

        auto to_z_z_xxxz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 23);

        auto to_z_z_xxxz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 24);

        auto to_z_z_xxxz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 25);

        auto to_z_z_xxxz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 26);

        auto to_z_z_xxxz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 27);

        auto to_z_z_xxxz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 28);

        auto to_z_z_xxxz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxx_xx, to_xxx_xxxz, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxx_zzzz, to_xxxzz_xx, to_xxxzz_xxxz, to_xxxzz_xxyz, to_xxxzz_xxzz, to_xxxzz_xy, to_xxxzz_xyyz, to_xxxzz_xyzz, to_xxxzz_xz, to_xxxzz_xzzz, to_xxxzz_yy, to_xxxzz_yyyz, to_xxxzz_yyzz, to_xxxzz_yz, to_xxxzz_yzzz, to_xxxzz_zz, to_xxxzz_zzzz, to_z_z_xxxz_xxx, to_z_z_xxxz_xxy, to_z_z_xxxz_xxz, to_z_z_xxxz_xyy, to_z_z_xxxz_xyz, to_z_z_xxxz_xzz, to_z_z_xxxz_yyy, to_z_z_xxxz_yyz, to_z_z_xxxz_yzz, to_z_z_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxz_xxx[k] = -2.0 * to_xxx_xxxz[k] * tke_0 + 4.0 * to_xxxzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxy[k] = -2.0 * to_xxx_xxyz[k] * tke_0 + 4.0 * to_xxxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xxz[k] = to_xxx_xx[k] - 2.0 * to_xxx_xxzz[k] * tke_0 - 2.0 * to_xxxzz_xx[k] * tbe_0 + 4.0 * to_xxxzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xyy[k] = -2.0 * to_xxx_xyyz[k] * tke_0 + 4.0 * to_xxxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xyz[k] = to_xxx_xy[k] - 2.0 * to_xxx_xyzz[k] * tke_0 - 2.0 * to_xxxzz_xy[k] * tbe_0 + 4.0 * to_xxxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xzz[k] = 2.0 * to_xxx_xz[k] - 2.0 * to_xxx_xzzz[k] * tke_0 - 4.0 * to_xxxzz_xz[k] * tbe_0 + 4.0 * to_xxxzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yyy[k] = -2.0 * to_xxx_yyyz[k] * tke_0 + 4.0 * to_xxxzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yyz[k] = to_xxx_yy[k] - 2.0 * to_xxx_yyzz[k] * tke_0 - 2.0 * to_xxxzz_yy[k] * tbe_0 + 4.0 * to_xxxzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yzz[k] = 2.0 * to_xxx_yz[k] - 2.0 * to_xxx_yzzz[k] * tke_0 - 4.0 * to_xxxzz_yz[k] * tbe_0 + 4.0 * to_xxxzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_zzz[k] = 3.0 * to_xxx_zz[k] - 2.0 * to_xxx_zzzz[k] * tke_0 - 6.0 * to_xxxzz_zz[k] * tbe_0 + 4.0 * to_xxxzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1230-1240 components of targeted buffer : GF

        auto to_z_z_xxyy_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 30);

        auto to_z_z_xxyy_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 31);

        auto to_z_z_xxyy_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 32);

        auto to_z_z_xxyy_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 33);

        auto to_z_z_xxyy_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 34);

        auto to_z_z_xxyy_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 35);

        auto to_z_z_xxyy_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 36);

        auto to_z_z_xxyy_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 37);

        auto to_z_z_xxyy_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 38);

        auto to_z_z_xxyy_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_xxyyz_xx, to_xxyyz_xxxz, to_xxyyz_xxyz, to_xxyyz_xxzz, to_xxyyz_xy, to_xxyyz_xyyz, to_xxyyz_xyzz, to_xxyyz_xz, to_xxyyz_xzzz, to_xxyyz_yy, to_xxyyz_yyyz, to_xxyyz_yyzz, to_xxyyz_yz, to_xxyyz_yzzz, to_xxyyz_zz, to_xxyyz_zzzz, to_z_z_xxyy_xxx, to_z_z_xxyy_xxy, to_z_z_xxyy_xxz, to_z_z_xxyy_xyy, to_z_z_xxyy_xyz, to_z_z_xxyy_xzz, to_z_z_xxyy_yyy, to_z_z_xxyy_yyz, to_z_z_xxyy_yzz, to_z_z_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyy_xxx[k] = 4.0 * to_xxyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxy[k] = 4.0 * to_xxyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xxz[k] = -2.0 * to_xxyyz_xx[k] * tbe_0 + 4.0 * to_xxyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xyy[k] = 4.0 * to_xxyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xyz[k] = -2.0 * to_xxyyz_xy[k] * tbe_0 + 4.0 * to_xxyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xzz[k] = -4.0 * to_xxyyz_xz[k] * tbe_0 + 4.0 * to_xxyyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yyy[k] = 4.0 * to_xxyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yyz[k] = -2.0 * to_xxyyz_yy[k] * tbe_0 + 4.0 * to_xxyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yzz[k] = -4.0 * to_xxyyz_yz[k] * tbe_0 + 4.0 * to_xxyyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_zzz[k] = -6.0 * to_xxyyz_zz[k] * tbe_0 + 4.0 * to_xxyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1240-1250 components of targeted buffer : GF

        auto to_z_z_xxyz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 40);

        auto to_z_z_xxyz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 41);

        auto to_z_z_xxyz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 42);

        auto to_z_z_xxyz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 43);

        auto to_z_z_xxyz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 44);

        auto to_z_z_xxyz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 45);

        auto to_z_z_xxyz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 46);

        auto to_z_z_xxyz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 47);

        auto to_z_z_xxyz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 48);

        auto to_z_z_xxyz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, to_xxyzz_xx, to_xxyzz_xxxz, to_xxyzz_xxyz, to_xxyzz_xxzz, to_xxyzz_xy, to_xxyzz_xyyz, to_xxyzz_xyzz, to_xxyzz_xz, to_xxyzz_xzzz, to_xxyzz_yy, to_xxyzz_yyyz, to_xxyzz_yyzz, to_xxyzz_yz, to_xxyzz_yzzz, to_xxyzz_zz, to_xxyzz_zzzz, to_z_z_xxyz_xxx, to_z_z_xxyz_xxy, to_z_z_xxyz_xxz, to_z_z_xxyz_xyy, to_z_z_xxyz_xyz, to_z_z_xxyz_xzz, to_z_z_xxyz_yyy, to_z_z_xxyz_yyz, to_z_z_xxyz_yzz, to_z_z_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyz_xxx[k] = -2.0 * to_xxy_xxxz[k] * tke_0 + 4.0 * to_xxyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxy[k] = -2.0 * to_xxy_xxyz[k] * tke_0 + 4.0 * to_xxyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xxz[k] = to_xxy_xx[k] - 2.0 * to_xxy_xxzz[k] * tke_0 - 2.0 * to_xxyzz_xx[k] * tbe_0 + 4.0 * to_xxyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xyy[k] = -2.0 * to_xxy_xyyz[k] * tke_0 + 4.0 * to_xxyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xyz[k] = to_xxy_xy[k] - 2.0 * to_xxy_xyzz[k] * tke_0 - 2.0 * to_xxyzz_xy[k] * tbe_0 + 4.0 * to_xxyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xzz[k] = 2.0 * to_xxy_xz[k] - 2.0 * to_xxy_xzzz[k] * tke_0 - 4.0 * to_xxyzz_xz[k] * tbe_0 + 4.0 * to_xxyzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yyy[k] = -2.0 * to_xxy_yyyz[k] * tke_0 + 4.0 * to_xxyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yyz[k] = to_xxy_yy[k] - 2.0 * to_xxy_yyzz[k] * tke_0 - 2.0 * to_xxyzz_yy[k] * tbe_0 + 4.0 * to_xxyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yzz[k] = 2.0 * to_xxy_yz[k] - 2.0 * to_xxy_yzzz[k] * tke_0 - 4.0 * to_xxyzz_yz[k] * tbe_0 + 4.0 * to_xxyzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_zzz[k] = 3.0 * to_xxy_zz[k] - 2.0 * to_xxy_zzzz[k] * tke_0 - 6.0 * to_xxyzz_zz[k] * tbe_0 + 4.0 * to_xxyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1250-1260 components of targeted buffer : GF

        auto to_z_z_xxzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 50);

        auto to_z_z_xxzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 51);

        auto to_z_z_xxzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 52);

        auto to_z_z_xxzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 53);

        auto to_z_z_xxzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 54);

        auto to_z_z_xxzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 55);

        auto to_z_z_xxzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 56);

        auto to_z_z_xxzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 57);

        auto to_z_z_xxzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 58);

        auto to_z_z_xxzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, to_xxzzz_xx, to_xxzzz_xxxz, to_xxzzz_xxyz, to_xxzzz_xxzz, to_xxzzz_xy, to_xxzzz_xyyz, to_xxzzz_xyzz, to_xxzzz_xz, to_xxzzz_xzzz, to_xxzzz_yy, to_xxzzz_yyyz, to_xxzzz_yyzz, to_xxzzz_yz, to_xxzzz_yzzz, to_xxzzz_zz, to_xxzzz_zzzz, to_z_z_xxzz_xxx, to_z_z_xxzz_xxy, to_z_z_xxzz_xxz, to_z_z_xxzz_xyy, to_z_z_xxzz_xyz, to_z_z_xxzz_xzz, to_z_z_xxzz_yyy, to_z_z_xxzz_yyz, to_z_z_xxzz_yzz, to_z_z_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxzz_xxx[k] = -4.0 * to_xxz_xxxz[k] * tke_0 + 4.0 * to_xxzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxy[k] = -4.0 * to_xxz_xxyz[k] * tke_0 + 4.0 * to_xxzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xxz[k] = 2.0 * to_xxz_xx[k] - 4.0 * to_xxz_xxzz[k] * tke_0 - 2.0 * to_xxzzz_xx[k] * tbe_0 + 4.0 * to_xxzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xyy[k] = -4.0 * to_xxz_xyyz[k] * tke_0 + 4.0 * to_xxzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xyz[k] = 2.0 * to_xxz_xy[k] - 4.0 * to_xxz_xyzz[k] * tke_0 - 2.0 * to_xxzzz_xy[k] * tbe_0 + 4.0 * to_xxzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xzz[k] = 4.0 * to_xxz_xz[k] - 4.0 * to_xxz_xzzz[k] * tke_0 - 4.0 * to_xxzzz_xz[k] * tbe_0 + 4.0 * to_xxzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yyy[k] = -4.0 * to_xxz_yyyz[k] * tke_0 + 4.0 * to_xxzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yyz[k] = 2.0 * to_xxz_yy[k] - 4.0 * to_xxz_yyzz[k] * tke_0 - 2.0 * to_xxzzz_yy[k] * tbe_0 + 4.0 * to_xxzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yzz[k] = 4.0 * to_xxz_yz[k] - 4.0 * to_xxz_yzzz[k] * tke_0 - 4.0 * to_xxzzz_yz[k] * tbe_0 + 4.0 * to_xxzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_zzz[k] = 6.0 * to_xxz_zz[k] - 4.0 * to_xxz_zzzz[k] * tke_0 - 6.0 * to_xxzzz_zz[k] * tbe_0 + 4.0 * to_xxzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1260-1270 components of targeted buffer : GF

        auto to_z_z_xyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 60);

        auto to_z_z_xyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 61);

        auto to_z_z_xyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 62);

        auto to_z_z_xyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 63);

        auto to_z_z_xyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 64);

        auto to_z_z_xyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 65);

        auto to_z_z_xyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 66);

        auto to_z_z_xyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 67);

        auto to_z_z_xyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 68);

        auto to_z_z_xyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_xyyyz_xx, to_xyyyz_xxxz, to_xyyyz_xxyz, to_xyyyz_xxzz, to_xyyyz_xy, to_xyyyz_xyyz, to_xyyyz_xyzz, to_xyyyz_xz, to_xyyyz_xzzz, to_xyyyz_yy, to_xyyyz_yyyz, to_xyyyz_yyzz, to_xyyyz_yz, to_xyyyz_yzzz, to_xyyyz_zz, to_xyyyz_zzzz, to_z_z_xyyy_xxx, to_z_z_xyyy_xxy, to_z_z_xyyy_xxz, to_z_z_xyyy_xyy, to_z_z_xyyy_xyz, to_z_z_xyyy_xzz, to_z_z_xyyy_yyy, to_z_z_xyyy_yyz, to_z_z_xyyy_yzz, to_z_z_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyy_xxx[k] = 4.0 * to_xyyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxy[k] = 4.0 * to_xyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xxz[k] = -2.0 * to_xyyyz_xx[k] * tbe_0 + 4.0 * to_xyyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xyy[k] = 4.0 * to_xyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xyz[k] = -2.0 * to_xyyyz_xy[k] * tbe_0 + 4.0 * to_xyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xzz[k] = -4.0 * to_xyyyz_xz[k] * tbe_0 + 4.0 * to_xyyyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yyy[k] = 4.0 * to_xyyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yyz[k] = -2.0 * to_xyyyz_yy[k] * tbe_0 + 4.0 * to_xyyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yzz[k] = -4.0 * to_xyyyz_yz[k] * tbe_0 + 4.0 * to_xyyyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_zzz[k] = -6.0 * to_xyyyz_zz[k] * tbe_0 + 4.0 * to_xyyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1270-1280 components of targeted buffer : GF

        auto to_z_z_xyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 70);

        auto to_z_z_xyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 71);

        auto to_z_z_xyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 72);

        auto to_z_z_xyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 73);

        auto to_z_z_xyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 74);

        auto to_z_z_xyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 75);

        auto to_z_z_xyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 76);

        auto to_z_z_xyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 77);

        auto to_z_z_xyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 78);

        auto to_z_z_xyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, to_xyyzz_xx, to_xyyzz_xxxz, to_xyyzz_xxyz, to_xyyzz_xxzz, to_xyyzz_xy, to_xyyzz_xyyz, to_xyyzz_xyzz, to_xyyzz_xz, to_xyyzz_xzzz, to_xyyzz_yy, to_xyyzz_yyyz, to_xyyzz_yyzz, to_xyyzz_yz, to_xyyzz_yzzz, to_xyyzz_zz, to_xyyzz_zzzz, to_z_z_xyyz_xxx, to_z_z_xyyz_xxy, to_z_z_xyyz_xxz, to_z_z_xyyz_xyy, to_z_z_xyyz_xyz, to_z_z_xyyz_xzz, to_z_z_xyyz_yyy, to_z_z_xyyz_yyz, to_z_z_xyyz_yzz, to_z_z_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyz_xxx[k] = -2.0 * to_xyy_xxxz[k] * tke_0 + 4.0 * to_xyyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxy[k] = -2.0 * to_xyy_xxyz[k] * tke_0 + 4.0 * to_xyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xxz[k] = to_xyy_xx[k] - 2.0 * to_xyy_xxzz[k] * tke_0 - 2.0 * to_xyyzz_xx[k] * tbe_0 + 4.0 * to_xyyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xyy[k] = -2.0 * to_xyy_xyyz[k] * tke_0 + 4.0 * to_xyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xyz[k] = to_xyy_xy[k] - 2.0 * to_xyy_xyzz[k] * tke_0 - 2.0 * to_xyyzz_xy[k] * tbe_0 + 4.0 * to_xyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xzz[k] = 2.0 * to_xyy_xz[k] - 2.0 * to_xyy_xzzz[k] * tke_0 - 4.0 * to_xyyzz_xz[k] * tbe_0 + 4.0 * to_xyyzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yyy[k] = -2.0 * to_xyy_yyyz[k] * tke_0 + 4.0 * to_xyyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yyz[k] = to_xyy_yy[k] - 2.0 * to_xyy_yyzz[k] * tke_0 - 2.0 * to_xyyzz_yy[k] * tbe_0 + 4.0 * to_xyyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yzz[k] = 2.0 * to_xyy_yz[k] - 2.0 * to_xyy_yzzz[k] * tke_0 - 4.0 * to_xyyzz_yz[k] * tbe_0 + 4.0 * to_xyyzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_zzz[k] = 3.0 * to_xyy_zz[k] - 2.0 * to_xyy_zzzz[k] * tke_0 - 6.0 * to_xyyzz_zz[k] * tbe_0 + 4.0 * to_xyyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1280-1290 components of targeted buffer : GF

        auto to_z_z_xyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 80);

        auto to_z_z_xyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 81);

        auto to_z_z_xyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 82);

        auto to_z_z_xyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 83);

        auto to_z_z_xyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 84);

        auto to_z_z_xyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 85);

        auto to_z_z_xyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 86);

        auto to_z_z_xyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 87);

        auto to_z_z_xyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 88);

        auto to_z_z_xyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, to_xyzzz_xx, to_xyzzz_xxxz, to_xyzzz_xxyz, to_xyzzz_xxzz, to_xyzzz_xy, to_xyzzz_xyyz, to_xyzzz_xyzz, to_xyzzz_xz, to_xyzzz_xzzz, to_xyzzz_yy, to_xyzzz_yyyz, to_xyzzz_yyzz, to_xyzzz_yz, to_xyzzz_yzzz, to_xyzzz_zz, to_xyzzz_zzzz, to_z_z_xyzz_xxx, to_z_z_xyzz_xxy, to_z_z_xyzz_xxz, to_z_z_xyzz_xyy, to_z_z_xyzz_xyz, to_z_z_xyzz_xzz, to_z_z_xyzz_yyy, to_z_z_xyzz_yyz, to_z_z_xyzz_yzz, to_z_z_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyzz_xxx[k] = -4.0 * to_xyz_xxxz[k] * tke_0 + 4.0 * to_xyzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxy[k] = -4.0 * to_xyz_xxyz[k] * tke_0 + 4.0 * to_xyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xxz[k] = 2.0 * to_xyz_xx[k] - 4.0 * to_xyz_xxzz[k] * tke_0 - 2.0 * to_xyzzz_xx[k] * tbe_0 + 4.0 * to_xyzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xyy[k] = -4.0 * to_xyz_xyyz[k] * tke_0 + 4.0 * to_xyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xyz[k] = 2.0 * to_xyz_xy[k] - 4.0 * to_xyz_xyzz[k] * tke_0 - 2.0 * to_xyzzz_xy[k] * tbe_0 + 4.0 * to_xyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xzz[k] = 4.0 * to_xyz_xz[k] - 4.0 * to_xyz_xzzz[k] * tke_0 - 4.0 * to_xyzzz_xz[k] * tbe_0 + 4.0 * to_xyzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yyy[k] = -4.0 * to_xyz_yyyz[k] * tke_0 + 4.0 * to_xyzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yyz[k] = 2.0 * to_xyz_yy[k] - 4.0 * to_xyz_yyzz[k] * tke_0 - 2.0 * to_xyzzz_yy[k] * tbe_0 + 4.0 * to_xyzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yzz[k] = 4.0 * to_xyz_yz[k] - 4.0 * to_xyz_yzzz[k] * tke_0 - 4.0 * to_xyzzz_yz[k] * tbe_0 + 4.0 * to_xyzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_zzz[k] = 6.0 * to_xyz_zz[k] - 4.0 * to_xyz_zzzz[k] * tke_0 - 6.0 * to_xyzzz_zz[k] * tbe_0 + 4.0 * to_xyzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1290-1300 components of targeted buffer : GF

        auto to_z_z_xzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 90);

        auto to_z_z_xzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 91);

        auto to_z_z_xzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 92);

        auto to_z_z_xzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 93);

        auto to_z_z_xzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 94);

        auto to_z_z_xzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 95);

        auto to_z_z_xzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 96);

        auto to_z_z_xzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 97);

        auto to_z_z_xzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 98);

        auto to_z_z_xzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, to_xzzzz_xx, to_xzzzz_xxxz, to_xzzzz_xxyz, to_xzzzz_xxzz, to_xzzzz_xy, to_xzzzz_xyyz, to_xzzzz_xyzz, to_xzzzz_xz, to_xzzzz_xzzz, to_xzzzz_yy, to_xzzzz_yyyz, to_xzzzz_yyzz, to_xzzzz_yz, to_xzzzz_yzzz, to_xzzzz_zz, to_xzzzz_zzzz, to_z_z_xzzz_xxx, to_z_z_xzzz_xxy, to_z_z_xzzz_xxz, to_z_z_xzzz_xyy, to_z_z_xzzz_xyz, to_z_z_xzzz_xzz, to_z_z_xzzz_yyy, to_z_z_xzzz_yyz, to_z_z_xzzz_yzz, to_z_z_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzzz_xxx[k] = -6.0 * to_xzz_xxxz[k] * tke_0 + 4.0 * to_xzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxy[k] = -6.0 * to_xzz_xxyz[k] * tke_0 + 4.0 * to_xzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xxz[k] = 3.0 * to_xzz_xx[k] - 6.0 * to_xzz_xxzz[k] * tke_0 - 2.0 * to_xzzzz_xx[k] * tbe_0 + 4.0 * to_xzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xyy[k] = -6.0 * to_xzz_xyyz[k] * tke_0 + 4.0 * to_xzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xyz[k] = 3.0 * to_xzz_xy[k] - 6.0 * to_xzz_xyzz[k] * tke_0 - 2.0 * to_xzzzz_xy[k] * tbe_0 + 4.0 * to_xzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xzz[k] = 6.0 * to_xzz_xz[k] - 6.0 * to_xzz_xzzz[k] * tke_0 - 4.0 * to_xzzzz_xz[k] * tbe_0 + 4.0 * to_xzzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yyy[k] = -6.0 * to_xzz_yyyz[k] * tke_0 + 4.0 * to_xzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yyz[k] = 3.0 * to_xzz_yy[k] - 6.0 * to_xzz_yyzz[k] * tke_0 - 2.0 * to_xzzzz_yy[k] * tbe_0 + 4.0 * to_xzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yzz[k] = 6.0 * to_xzz_yz[k] - 6.0 * to_xzz_yzzz[k] * tke_0 - 4.0 * to_xzzzz_yz[k] * tbe_0 + 4.0 * to_xzzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_zzz[k] = 9.0 * to_xzz_zz[k] - 6.0 * to_xzz_zzzz[k] * tke_0 - 6.0 * to_xzzzz_zz[k] * tbe_0 + 4.0 * to_xzzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1300-1310 components of targeted buffer : GF

        auto to_z_z_yyyy_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 100);

        auto to_z_z_yyyy_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 101);

        auto to_z_z_yyyy_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 102);

        auto to_z_z_yyyy_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 103);

        auto to_z_z_yyyy_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 104);

        auto to_z_z_yyyy_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 105);

        auto to_z_z_yyyy_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 106);

        auto to_z_z_yyyy_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 107);

        auto to_z_z_yyyy_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 108);

        auto to_z_z_yyyy_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_yyyyz_xx, to_yyyyz_xxxz, to_yyyyz_xxyz, to_yyyyz_xxzz, to_yyyyz_xy, to_yyyyz_xyyz, to_yyyyz_xyzz, to_yyyyz_xz, to_yyyyz_xzzz, to_yyyyz_yy, to_yyyyz_yyyz, to_yyyyz_yyzz, to_yyyyz_yz, to_yyyyz_yzzz, to_yyyyz_zz, to_yyyyz_zzzz, to_z_z_yyyy_xxx, to_z_z_yyyy_xxy, to_z_z_yyyy_xxz, to_z_z_yyyy_xyy, to_z_z_yyyy_xyz, to_z_z_yyyy_xzz, to_z_z_yyyy_yyy, to_z_z_yyyy_yyz, to_z_z_yyyy_yzz, to_z_z_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyy_xxx[k] = 4.0 * to_yyyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxy[k] = 4.0 * to_yyyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xxz[k] = -2.0 * to_yyyyz_xx[k] * tbe_0 + 4.0 * to_yyyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xyy[k] = 4.0 * to_yyyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xyz[k] = -2.0 * to_yyyyz_xy[k] * tbe_0 + 4.0 * to_yyyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xzz[k] = -4.0 * to_yyyyz_xz[k] * tbe_0 + 4.0 * to_yyyyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yyy[k] = 4.0 * to_yyyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yyz[k] = -2.0 * to_yyyyz_yy[k] * tbe_0 + 4.0 * to_yyyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yzz[k] = -4.0 * to_yyyyz_yz[k] * tbe_0 + 4.0 * to_yyyyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_zzz[k] = -6.0 * to_yyyyz_zz[k] * tbe_0 + 4.0 * to_yyyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1310-1320 components of targeted buffer : GF

        auto to_z_z_yyyz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 110);

        auto to_z_z_yyyz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 111);

        auto to_z_z_yyyz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 112);

        auto to_z_z_yyyz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 113);

        auto to_z_z_yyyz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 114);

        auto to_z_z_yyyz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 115);

        auto to_z_z_yyyz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 116);

        auto to_z_z_yyyz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 117);

        auto to_z_z_yyyz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 118);

        auto to_z_z_yyyz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yyy_xx, to_yyy_xxxz, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyy_zzzz, to_yyyzz_xx, to_yyyzz_xxxz, to_yyyzz_xxyz, to_yyyzz_xxzz, to_yyyzz_xy, to_yyyzz_xyyz, to_yyyzz_xyzz, to_yyyzz_xz, to_yyyzz_xzzz, to_yyyzz_yy, to_yyyzz_yyyz, to_yyyzz_yyzz, to_yyyzz_yz, to_yyyzz_yzzz, to_yyyzz_zz, to_yyyzz_zzzz, to_z_z_yyyz_xxx, to_z_z_yyyz_xxy, to_z_z_yyyz_xxz, to_z_z_yyyz_xyy, to_z_z_yyyz_xyz, to_z_z_yyyz_xzz, to_z_z_yyyz_yyy, to_z_z_yyyz_yyz, to_z_z_yyyz_yzz, to_z_z_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyz_xxx[k] = -2.0 * to_yyy_xxxz[k] * tke_0 + 4.0 * to_yyyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxy[k] = -2.0 * to_yyy_xxyz[k] * tke_0 + 4.0 * to_yyyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xxz[k] = to_yyy_xx[k] - 2.0 * to_yyy_xxzz[k] * tke_0 - 2.0 * to_yyyzz_xx[k] * tbe_0 + 4.0 * to_yyyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xyy[k] = -2.0 * to_yyy_xyyz[k] * tke_0 + 4.0 * to_yyyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xyz[k] = to_yyy_xy[k] - 2.0 * to_yyy_xyzz[k] * tke_0 - 2.0 * to_yyyzz_xy[k] * tbe_0 + 4.0 * to_yyyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xzz[k] = 2.0 * to_yyy_xz[k] - 2.0 * to_yyy_xzzz[k] * tke_0 - 4.0 * to_yyyzz_xz[k] * tbe_0 + 4.0 * to_yyyzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yyy[k] = -2.0 * to_yyy_yyyz[k] * tke_0 + 4.0 * to_yyyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yyz[k] = to_yyy_yy[k] - 2.0 * to_yyy_yyzz[k] * tke_0 - 2.0 * to_yyyzz_yy[k] * tbe_0 + 4.0 * to_yyyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yzz[k] = 2.0 * to_yyy_yz[k] - 2.0 * to_yyy_yzzz[k] * tke_0 - 4.0 * to_yyyzz_yz[k] * tbe_0 + 4.0 * to_yyyzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_zzz[k] = 3.0 * to_yyy_zz[k] - 2.0 * to_yyy_zzzz[k] * tke_0 - 6.0 * to_yyyzz_zz[k] * tbe_0 + 4.0 * to_yyyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1320-1330 components of targeted buffer : GF

        auto to_z_z_yyzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 120);

        auto to_z_z_yyzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 121);

        auto to_z_z_yyzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 122);

        auto to_z_z_yyzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 123);

        auto to_z_z_yyzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 124);

        auto to_z_z_yyzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 125);

        auto to_z_z_yyzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 126);

        auto to_z_z_yyzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 127);

        auto to_z_z_yyzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 128);

        auto to_z_z_yyzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, to_yyzzz_xx, to_yyzzz_xxxz, to_yyzzz_xxyz, to_yyzzz_xxzz, to_yyzzz_xy, to_yyzzz_xyyz, to_yyzzz_xyzz, to_yyzzz_xz, to_yyzzz_xzzz, to_yyzzz_yy, to_yyzzz_yyyz, to_yyzzz_yyzz, to_yyzzz_yz, to_yyzzz_yzzz, to_yyzzz_zz, to_yyzzz_zzzz, to_z_z_yyzz_xxx, to_z_z_yyzz_xxy, to_z_z_yyzz_xxz, to_z_z_yyzz_xyy, to_z_z_yyzz_xyz, to_z_z_yyzz_xzz, to_z_z_yyzz_yyy, to_z_z_yyzz_yyz, to_z_z_yyzz_yzz, to_z_z_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyzz_xxx[k] = -4.0 * to_yyz_xxxz[k] * tke_0 + 4.0 * to_yyzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxy[k] = -4.0 * to_yyz_xxyz[k] * tke_0 + 4.0 * to_yyzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xxz[k] = 2.0 * to_yyz_xx[k] - 4.0 * to_yyz_xxzz[k] * tke_0 - 2.0 * to_yyzzz_xx[k] * tbe_0 + 4.0 * to_yyzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xyy[k] = -4.0 * to_yyz_xyyz[k] * tke_0 + 4.0 * to_yyzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xyz[k] = 2.0 * to_yyz_xy[k] - 4.0 * to_yyz_xyzz[k] * tke_0 - 2.0 * to_yyzzz_xy[k] * tbe_0 + 4.0 * to_yyzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xzz[k] = 4.0 * to_yyz_xz[k] - 4.0 * to_yyz_xzzz[k] * tke_0 - 4.0 * to_yyzzz_xz[k] * tbe_0 + 4.0 * to_yyzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yyy[k] = -4.0 * to_yyz_yyyz[k] * tke_0 + 4.0 * to_yyzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yyz[k] = 2.0 * to_yyz_yy[k] - 4.0 * to_yyz_yyzz[k] * tke_0 - 2.0 * to_yyzzz_yy[k] * tbe_0 + 4.0 * to_yyzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yzz[k] = 4.0 * to_yyz_yz[k] - 4.0 * to_yyz_yzzz[k] * tke_0 - 4.0 * to_yyzzz_yz[k] * tbe_0 + 4.0 * to_yyzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_zzz[k] = 6.0 * to_yyz_zz[k] - 4.0 * to_yyz_zzzz[k] * tke_0 - 6.0 * to_yyzzz_zz[k] * tbe_0 + 4.0 * to_yyzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1330-1340 components of targeted buffer : GF

        auto to_z_z_yzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 130);

        auto to_z_z_yzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 131);

        auto to_z_z_yzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 132);

        auto to_z_z_yzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 133);

        auto to_z_z_yzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 134);

        auto to_z_z_yzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 135);

        auto to_z_z_yzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 136);

        auto to_z_z_yzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 137);

        auto to_z_z_yzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 138);

        auto to_z_z_yzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, to_yzzzz_xx, to_yzzzz_xxxz, to_yzzzz_xxyz, to_yzzzz_xxzz, to_yzzzz_xy, to_yzzzz_xyyz, to_yzzzz_xyzz, to_yzzzz_xz, to_yzzzz_xzzz, to_yzzzz_yy, to_yzzzz_yyyz, to_yzzzz_yyzz, to_yzzzz_yz, to_yzzzz_yzzz, to_yzzzz_zz, to_yzzzz_zzzz, to_z_z_yzzz_xxx, to_z_z_yzzz_xxy, to_z_z_yzzz_xxz, to_z_z_yzzz_xyy, to_z_z_yzzz_xyz, to_z_z_yzzz_xzz, to_z_z_yzzz_yyy, to_z_z_yzzz_yyz, to_z_z_yzzz_yzz, to_z_z_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzzz_xxx[k] = -6.0 * to_yzz_xxxz[k] * tke_0 + 4.0 * to_yzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxy[k] = -6.0 * to_yzz_xxyz[k] * tke_0 + 4.0 * to_yzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xxz[k] = 3.0 * to_yzz_xx[k] - 6.0 * to_yzz_xxzz[k] * tke_0 - 2.0 * to_yzzzz_xx[k] * tbe_0 + 4.0 * to_yzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xyy[k] = -6.0 * to_yzz_xyyz[k] * tke_0 + 4.0 * to_yzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xyz[k] = 3.0 * to_yzz_xy[k] - 6.0 * to_yzz_xyzz[k] * tke_0 - 2.0 * to_yzzzz_xy[k] * tbe_0 + 4.0 * to_yzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xzz[k] = 6.0 * to_yzz_xz[k] - 6.0 * to_yzz_xzzz[k] * tke_0 - 4.0 * to_yzzzz_xz[k] * tbe_0 + 4.0 * to_yzzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yyy[k] = -6.0 * to_yzz_yyyz[k] * tke_0 + 4.0 * to_yzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yyz[k] = 3.0 * to_yzz_yy[k] - 6.0 * to_yzz_yyzz[k] * tke_0 - 2.0 * to_yzzzz_yy[k] * tbe_0 + 4.0 * to_yzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yzz[k] = 6.0 * to_yzz_yz[k] - 6.0 * to_yzz_yzzz[k] * tke_0 - 4.0 * to_yzzzz_yz[k] * tbe_0 + 4.0 * to_yzzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_zzz[k] = 9.0 * to_yzz_zz[k] - 6.0 * to_yzz_zzzz[k] * tke_0 - 6.0 * to_yzzzz_zz[k] * tbe_0 + 4.0 * to_yzzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1340-1350 components of targeted buffer : GF

        auto to_z_z_zzzz_xxx = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 140);

        auto to_z_z_zzzz_xxy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 141);

        auto to_z_z_zzzz_xxz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 142);

        auto to_z_z_zzzz_xyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 143);

        auto to_z_z_zzzz_xyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 144);

        auto to_z_z_zzzz_xzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 145);

        auto to_z_z_zzzz_yyy = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 146);

        auto to_z_z_zzzz_yyz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 147);

        auto to_z_z_zzzz_yzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 148);

        auto to_z_z_zzzz_zzz = pbuffer.data(idx_op_geom_101_gf + 8 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_z_zzzz_xxx, to_z_z_zzzz_xxy, to_z_z_zzzz_xxz, to_z_z_zzzz_xyy, to_z_z_zzzz_xyz, to_z_z_zzzz_xzz, to_z_z_zzzz_yyy, to_z_z_zzzz_yyz, to_z_z_zzzz_yzz, to_z_z_zzzz_zzz, to_zzz_xx, to_zzz_xxxz, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzz_zzzz, to_zzzzz_xx, to_zzzzz_xxxz, to_zzzzz_xxyz, to_zzzzz_xxzz, to_zzzzz_xy, to_zzzzz_xyyz, to_zzzzz_xyzz, to_zzzzz_xz, to_zzzzz_xzzz, to_zzzzz_yy, to_zzzzz_yyyz, to_zzzzz_yyzz, to_zzzzz_yz, to_zzzzz_yzzz, to_zzzzz_zz, to_zzzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzzz_xxx[k] = -8.0 * to_zzz_xxxz[k] * tke_0 + 4.0 * to_zzzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxy[k] = -8.0 * to_zzz_xxyz[k] * tke_0 + 4.0 * to_zzzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xxz[k] = 4.0 * to_zzz_xx[k] - 8.0 * to_zzz_xxzz[k] * tke_0 - 2.0 * to_zzzzz_xx[k] * tbe_0 + 4.0 * to_zzzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xyy[k] = -8.0 * to_zzz_xyyz[k] * tke_0 + 4.0 * to_zzzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xyz[k] = 4.0 * to_zzz_xy[k] - 8.0 * to_zzz_xyzz[k] * tke_0 - 2.0 * to_zzzzz_xy[k] * tbe_0 + 4.0 * to_zzzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xzz[k] = 8.0 * to_zzz_xz[k] - 8.0 * to_zzz_xzzz[k] * tke_0 - 4.0 * to_zzzzz_xz[k] * tbe_0 + 4.0 * to_zzzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yyy[k] = -8.0 * to_zzz_yyyz[k] * tke_0 + 4.0 * to_zzzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yyz[k] = 4.0 * to_zzz_yy[k] - 8.0 * to_zzz_yyzz[k] * tke_0 - 2.0 * to_zzzzz_yy[k] * tbe_0 + 4.0 * to_zzzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yzz[k] = 8.0 * to_zzz_yz[k] - 8.0 * to_zzz_yzzz[k] * tke_0 - 4.0 * to_zzzzz_yz[k] * tbe_0 + 4.0 * to_zzzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_zzz[k] = 12.0 * to_zzz_zz[k] - 8.0 * to_zzz_zzzz[k] * tke_0 - 6.0 * to_zzzzz_zz[k] * tbe_0 + 4.0 * to_zzzzz_zzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

