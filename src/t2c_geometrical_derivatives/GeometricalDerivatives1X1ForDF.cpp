#include "GeometricalDerivatives1X1ForDF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_df(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_df,
                        const size_t idx_op_pd,
                        const size_t idx_op_pg,
                        const size_t idx_op_fd,
                        const size_t idx_op_fg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PD

        auto to_x_xx = pbuffer.data(idx_op_pd + i * 18 + 0);

        auto to_x_xy = pbuffer.data(idx_op_pd + i * 18 + 1);

        auto to_x_xz = pbuffer.data(idx_op_pd + i * 18 + 2);

        auto to_x_yy = pbuffer.data(idx_op_pd + i * 18 + 3);

        auto to_x_yz = pbuffer.data(idx_op_pd + i * 18 + 4);

        auto to_x_zz = pbuffer.data(idx_op_pd + i * 18 + 5);

        auto to_y_xx = pbuffer.data(idx_op_pd + i * 18 + 6);

        auto to_y_xy = pbuffer.data(idx_op_pd + i * 18 + 7);

        auto to_y_xz = pbuffer.data(idx_op_pd + i * 18 + 8);

        auto to_y_yy = pbuffer.data(idx_op_pd + i * 18 + 9);

        auto to_y_yz = pbuffer.data(idx_op_pd + i * 18 + 10);

        auto to_y_zz = pbuffer.data(idx_op_pd + i * 18 + 11);

        auto to_z_xx = pbuffer.data(idx_op_pd + i * 18 + 12);

        auto to_z_xy = pbuffer.data(idx_op_pd + i * 18 + 13);

        auto to_z_xz = pbuffer.data(idx_op_pd + i * 18 + 14);

        auto to_z_yy = pbuffer.data(idx_op_pd + i * 18 + 15);

        auto to_z_yz = pbuffer.data(idx_op_pd + i * 18 + 16);

        auto to_z_zz = pbuffer.data(idx_op_pd + i * 18 + 17);

        // Set up components of auxiliary buffer : PG

        auto to_x_xxxx = pbuffer.data(idx_op_pg + i * 45 + 0);

        auto to_x_xxxy = pbuffer.data(idx_op_pg + i * 45 + 1);

        auto to_x_xxxz = pbuffer.data(idx_op_pg + i * 45 + 2);

        auto to_x_xxyy = pbuffer.data(idx_op_pg + i * 45 + 3);

        auto to_x_xxyz = pbuffer.data(idx_op_pg + i * 45 + 4);

        auto to_x_xxzz = pbuffer.data(idx_op_pg + i * 45 + 5);

        auto to_x_xyyy = pbuffer.data(idx_op_pg + i * 45 + 6);

        auto to_x_xyyz = pbuffer.data(idx_op_pg + i * 45 + 7);

        auto to_x_xyzz = pbuffer.data(idx_op_pg + i * 45 + 8);

        auto to_x_xzzz = pbuffer.data(idx_op_pg + i * 45 + 9);

        auto to_x_yyyy = pbuffer.data(idx_op_pg + i * 45 + 10);

        auto to_x_yyyz = pbuffer.data(idx_op_pg + i * 45 + 11);

        auto to_x_yyzz = pbuffer.data(idx_op_pg + i * 45 + 12);

        auto to_x_yzzz = pbuffer.data(idx_op_pg + i * 45 + 13);

        auto to_x_zzzz = pbuffer.data(idx_op_pg + i * 45 + 14);

        auto to_y_xxxx = pbuffer.data(idx_op_pg + i * 45 + 15);

        auto to_y_xxxy = pbuffer.data(idx_op_pg + i * 45 + 16);

        auto to_y_xxxz = pbuffer.data(idx_op_pg + i * 45 + 17);

        auto to_y_xxyy = pbuffer.data(idx_op_pg + i * 45 + 18);

        auto to_y_xxyz = pbuffer.data(idx_op_pg + i * 45 + 19);

        auto to_y_xxzz = pbuffer.data(idx_op_pg + i * 45 + 20);

        auto to_y_xyyy = pbuffer.data(idx_op_pg + i * 45 + 21);

        auto to_y_xyyz = pbuffer.data(idx_op_pg + i * 45 + 22);

        auto to_y_xyzz = pbuffer.data(idx_op_pg + i * 45 + 23);

        auto to_y_xzzz = pbuffer.data(idx_op_pg + i * 45 + 24);

        auto to_y_yyyy = pbuffer.data(idx_op_pg + i * 45 + 25);

        auto to_y_yyyz = pbuffer.data(idx_op_pg + i * 45 + 26);

        auto to_y_yyzz = pbuffer.data(idx_op_pg + i * 45 + 27);

        auto to_y_yzzz = pbuffer.data(idx_op_pg + i * 45 + 28);

        auto to_y_zzzz = pbuffer.data(idx_op_pg + i * 45 + 29);

        auto to_z_xxxx = pbuffer.data(idx_op_pg + i * 45 + 30);

        auto to_z_xxxy = pbuffer.data(idx_op_pg + i * 45 + 31);

        auto to_z_xxxz = pbuffer.data(idx_op_pg + i * 45 + 32);

        auto to_z_xxyy = pbuffer.data(idx_op_pg + i * 45 + 33);

        auto to_z_xxyz = pbuffer.data(idx_op_pg + i * 45 + 34);

        auto to_z_xxzz = pbuffer.data(idx_op_pg + i * 45 + 35);

        auto to_z_xyyy = pbuffer.data(idx_op_pg + i * 45 + 36);

        auto to_z_xyyz = pbuffer.data(idx_op_pg + i * 45 + 37);

        auto to_z_xyzz = pbuffer.data(idx_op_pg + i * 45 + 38);

        auto to_z_xzzz = pbuffer.data(idx_op_pg + i * 45 + 39);

        auto to_z_yyyy = pbuffer.data(idx_op_pg + i * 45 + 40);

        auto to_z_yyyz = pbuffer.data(idx_op_pg + i * 45 + 41);

        auto to_z_yyzz = pbuffer.data(idx_op_pg + i * 45 + 42);

        auto to_z_yzzz = pbuffer.data(idx_op_pg + i * 45 + 43);

        auto to_z_zzzz = pbuffer.data(idx_op_pg + i * 45 + 44);

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

        // Set up 0-10 components of targeted buffer : DF

        auto to_x_x_xx_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 0);

        auto to_x_x_xx_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 1);

        auto to_x_x_xx_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 2);

        auto to_x_x_xx_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 3);

        auto to_x_x_xx_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 4);

        auto to_x_x_xx_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 5);

        auto to_x_x_xx_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 6);

        auto to_x_x_xx_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 7);

        auto to_x_x_xx_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 8);

        auto to_x_x_xx_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_x_x_xx_xxx, to_x_x_xx_xxy, to_x_x_xx_xxz, to_x_x_xx_xyy, to_x_x_xx_xyz, to_x_x_xx_xzz, to_x_x_xx_yyy, to_x_x_xx_yyz, to_x_x_xx_yzz, to_x_x_xx_zzz, to_x_xx, to_x_xxxx, to_x_xxxy, to_x_xxxz, to_x_xxyy, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yz, to_x_zz, to_xxx_xx, to_xxx_xxxx, to_xxx_xxxy, to_xxx_xxxz, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xx_xxx[k] = 6.0 * to_x_xx[k] - 4.0 * to_x_xxxx[k] * tke_0 - 6.0 * to_xxx_xx[k] * tbe_0 + 4.0 * to_xxx_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xx_xxy[k] = 4.0 * to_x_xy[k] - 4.0 * to_x_xxxy[k] * tke_0 - 4.0 * to_xxx_xy[k] * tbe_0 + 4.0 * to_xxx_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xx_xxz[k] = 4.0 * to_x_xz[k] - 4.0 * to_x_xxxz[k] * tke_0 - 4.0 * to_xxx_xz[k] * tbe_0 + 4.0 * to_xxx_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xx_xyy[k] = 2.0 * to_x_yy[k] - 4.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xxx_yy[k] * tbe_0 + 4.0 * to_xxx_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xx_xyz[k] = 2.0 * to_x_yz[k] - 4.0 * to_x_xxyz[k] * tke_0 - 2.0 * to_xxx_yz[k] * tbe_0 + 4.0 * to_xxx_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xx_xzz[k] = 2.0 * to_x_zz[k] - 4.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xxx_zz[k] * tbe_0 + 4.0 * to_xxx_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xx_yyy[k] = -4.0 * to_x_xyyy[k] * tke_0 + 4.0 * to_xxx_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xx_yyz[k] = -4.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xxx_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xx_yzz[k] = -4.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xxx_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xx_zzz[k] = -4.0 * to_x_xzzz[k] * tke_0 + 4.0 * to_xxx_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : DF

        auto to_x_x_xy_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 10);

        auto to_x_x_xy_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 11);

        auto to_x_x_xy_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 12);

        auto to_x_x_xy_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 13);

        auto to_x_x_xy_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 14);

        auto to_x_x_xy_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 15);

        auto to_x_x_xy_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 16);

        auto to_x_x_xy_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 17);

        auto to_x_x_xy_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 18);

        auto to_x_x_xy_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_x_xy_xxx, to_x_x_xy_xxy, to_x_x_xy_xxz, to_x_x_xy_xyy, to_x_x_xy_xyz, to_x_x_xy_xzz, to_x_x_xy_yyy, to_x_x_xy_yyz, to_x_x_xy_yzz, to_x_x_xy_zzz, to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, to_y_xx, to_y_xxxx, to_y_xxxy, to_y_xxxz, to_y_xxyy, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xy_xxx[k] = 3.0 * to_y_xx[k] - 2.0 * to_y_xxxx[k] * tke_0 - 6.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xy_xxy[k] = 2.0 * to_y_xy[k] - 2.0 * to_y_xxxy[k] * tke_0 - 4.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xy_xxz[k] = 2.0 * to_y_xz[k] - 2.0 * to_y_xxxz[k] * tke_0 - 4.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xy_xyy[k] = to_y_yy[k] - 2.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xy_xyz[k] = to_y_yz[k] - 2.0 * to_y_xxyz[k] * tke_0 - 2.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xy_xzz[k] = to_y_zz[k] - 2.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xy_yyy[k] = -2.0 * to_y_xyyy[k] * tke_0 + 4.0 * to_xxy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xy_yyz[k] = -2.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xy_yzz[k] = -2.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xy_zzz[k] = -2.0 * to_y_xzzz[k] * tke_0 + 4.0 * to_xxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : DF

        auto to_x_x_xz_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 20);

        auto to_x_x_xz_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 21);

        auto to_x_x_xz_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 22);

        auto to_x_x_xz_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 23);

        auto to_x_x_xz_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 24);

        auto to_x_x_xz_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 25);

        auto to_x_x_xz_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 26);

        auto to_x_x_xz_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 27);

        auto to_x_x_xz_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 28);

        auto to_x_x_xz_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_x_xz_xxx, to_x_x_xz_xxy, to_x_x_xz_xxz, to_x_x_xz_xyy, to_x_x_xz_xyz, to_x_x_xz_xzz, to_x_x_xz_yyy, to_x_x_xz_yyz, to_x_x_xz_yzz, to_x_x_xz_zzz, to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, to_z_xx, to_z_xxxx, to_z_xxxy, to_z_xxxz, to_z_xxyy, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xz_xxx[k] = 3.0 * to_z_xx[k] - 2.0 * to_z_xxxx[k] * tke_0 - 6.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xz_xxy[k] = 2.0 * to_z_xy[k] - 2.0 * to_z_xxxy[k] * tke_0 - 4.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xz_xxz[k] = 2.0 * to_z_xz[k] - 2.0 * to_z_xxxz[k] * tke_0 - 4.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xz_xyy[k] = to_z_yy[k] - 2.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xz_xyz[k] = to_z_yz[k] - 2.0 * to_z_xxyz[k] * tke_0 - 2.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xz_xzz[k] = to_z_zz[k] - 2.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xz_yyy[k] = -2.0 * to_z_xyyy[k] * tke_0 + 4.0 * to_xxz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xz_yyz[k] = -2.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xz_yzz[k] = -2.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xz_zzz[k] = -2.0 * to_z_xzzz[k] * tke_0 + 4.0 * to_xxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : DF

        auto to_x_x_yy_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 30);

        auto to_x_x_yy_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 31);

        auto to_x_x_yy_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 32);

        auto to_x_x_yy_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 33);

        auto to_x_x_yy_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 34);

        auto to_x_x_yy_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 35);

        auto to_x_x_yy_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 36);

        auto to_x_x_yy_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 37);

        auto to_x_x_yy_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 38);

        auto to_x_x_yy_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_x_x_yy_xxx, to_x_x_yy_xxy, to_x_x_yy_xxz, to_x_x_yy_xyy, to_x_x_yy_xyz, to_x_x_yy_xzz, to_x_x_yy_yyy, to_x_x_yy_yyz, to_x_x_yy_yzz, to_x_x_yy_zzz, to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yy_xxx[k] = -6.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yy_xxy[k] = -4.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yy_xxz[k] = -4.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yy_xyy[k] = -2.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yy_xyz[k] = -2.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yy_xzz[k] = -2.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yy_yyy[k] = 4.0 * to_xyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yy_yyz[k] = 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yy_yzz[k] = 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yy_zzz[k] = 4.0 * to_xyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : DF

        auto to_x_x_yz_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 40);

        auto to_x_x_yz_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 41);

        auto to_x_x_yz_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 42);

        auto to_x_x_yz_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 43);

        auto to_x_x_yz_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 44);

        auto to_x_x_yz_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 45);

        auto to_x_x_yz_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 46);

        auto to_x_x_yz_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 47);

        auto to_x_x_yz_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 48);

        auto to_x_x_yz_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_x_x_yz_xxx, to_x_x_yz_xxy, to_x_x_yz_xxz, to_x_x_yz_xyy, to_x_x_yz_xyz, to_x_x_yz_xzz, to_x_x_yz_yyy, to_x_x_yz_yyz, to_x_x_yz_yzz, to_x_x_yz_zzz, to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yz_xxx[k] = -6.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yz_xxy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yz_xxz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yz_xyy[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yz_xyz[k] = -2.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yz_xzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yz_yyy[k] = 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yz_yyz[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yz_yzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yz_zzz[k] = 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : DF

        auto to_x_x_zz_xxx = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 50);

        auto to_x_x_zz_xxy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 51);

        auto to_x_x_zz_xxz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 52);

        auto to_x_x_zz_xyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 53);

        auto to_x_x_zz_xyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 54);

        auto to_x_x_zz_xzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 55);

        auto to_x_x_zz_yyy = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 56);

        auto to_x_x_zz_yyz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 57);

        auto to_x_x_zz_yzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 58);

        auto to_x_x_zz_zzz = pbuffer.data(idx_op_geom_101_df + 0 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_x_zz_xxx, to_x_x_zz_xxy, to_x_x_zz_xxz, to_x_x_zz_xyy, to_x_x_zz_xyz, to_x_x_zz_xzz, to_x_x_zz_yyy, to_x_x_zz_yyz, to_x_x_zz_yzz, to_x_x_zz_zzz, to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zz_xxx[k] = -6.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_zz_xxy[k] = -4.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_zz_xxz[k] = -4.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_zz_xyy[k] = -2.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_zz_xyz[k] = -2.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_zz_xzz[k] = -2.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_zz_yyy[k] = 4.0 * to_xzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_zz_yyz[k] = 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_zz_yzz[k] = 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_zz_zzz[k] = 4.0 * to_xzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : DF

        auto to_x_y_xx_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 0);

        auto to_x_y_xx_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 1);

        auto to_x_y_xx_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 2);

        auto to_x_y_xx_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 3);

        auto to_x_y_xx_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 4);

        auto to_x_y_xx_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 5);

        auto to_x_y_xx_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 6);

        auto to_x_y_xx_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 7);

        auto to_x_y_xx_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 8);

        auto to_x_y_xx_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_x_xx, to_x_xxxy, to_x_xxyy, to_x_xxyz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_y_xx_xxx, to_x_y_xx_xxy, to_x_y_xx_xxz, to_x_y_xx_xyy, to_x_y_xx_xyz, to_x_y_xx_xzz, to_x_y_xx_yyy, to_x_y_xx_yyz, to_x_y_xx_yzz, to_x_y_xx_zzz, to_x_yy, to_x_yyyy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_xxx_xx, to_xxx_xxxy, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_yy, to_xxx_yyyy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xx_xxx[k] = -4.0 * to_x_xxxy[k] * tke_0 + 4.0 * to_xxx_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xx_xxy[k] = 2.0 * to_x_xx[k] - 4.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xxx_xx[k] * tbe_0 + 4.0 * to_xxx_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xxz[k] = -4.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xxx_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xx_xyy[k] = 4.0 * to_x_xy[k] - 4.0 * to_x_xyyy[k] * tke_0 - 4.0 * to_xxx_xy[k] * tbe_0 + 4.0 * to_xxx_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xx_xyz[k] = 2.0 * to_x_xz[k] - 4.0 * to_x_xyyz[k] * tke_0 - 2.0 * to_xxx_xz[k] * tbe_0 + 4.0 * to_xxx_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xx_xzz[k] = -4.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xxx_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xx_yyy[k] = 6.0 * to_x_yy[k] - 4.0 * to_x_yyyy[k] * tke_0 - 6.0 * to_xxx_yy[k] * tbe_0 + 4.0 * to_xxx_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xx_yyz[k] = 4.0 * to_x_yz[k] - 4.0 * to_x_yyyz[k] * tke_0 - 4.0 * to_xxx_yz[k] * tbe_0 + 4.0 * to_xxx_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xx_yzz[k] = 2.0 * to_x_zz[k] - 4.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xxx_zz[k] * tbe_0 + 4.0 * to_xxx_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xx_zzz[k] = -4.0 * to_x_yzzz[k] * tke_0 + 4.0 * to_xxx_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : DF

        auto to_x_y_xy_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 10);

        auto to_x_y_xy_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 11);

        auto to_x_y_xy_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 12);

        auto to_x_y_xy_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 13);

        auto to_x_y_xy_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 14);

        auto to_x_y_xy_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 15);

        auto to_x_y_xy_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 16);

        auto to_x_y_xy_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 17);

        auto to_x_y_xy_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 18);

        auto to_x_y_xy_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_y_xy_xxx, to_x_y_xy_xxy, to_x_y_xy_xxz, to_x_y_xy_xyy, to_x_y_xy_xyz, to_x_y_xy_xzz, to_x_y_xy_yyy, to_x_y_xy_yyz, to_x_y_xy_yzz, to_x_y_xy_zzz, to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_y_xx, to_y_xxxy, to_y_xxyy, to_y_xxyz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_yy, to_y_yyyy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xy_xxx[k] = -2.0 * to_y_xxxy[k] * tke_0 + 4.0 * to_xxy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xy_xxy[k] = to_y_xx[k] - 2.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xxz[k] = -2.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xy_xyy[k] = 2.0 * to_y_xy[k] - 2.0 * to_y_xyyy[k] * tke_0 - 4.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xy_xyz[k] = to_y_xz[k] - 2.0 * to_y_xyyz[k] * tke_0 - 2.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xy_xzz[k] = -2.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xy_yyy[k] = 3.0 * to_y_yy[k] - 2.0 * to_y_yyyy[k] * tke_0 - 6.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xy_yyz[k] = 2.0 * to_y_yz[k] - 2.0 * to_y_yyyz[k] * tke_0 - 4.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xy_yzz[k] = to_y_zz[k] - 2.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xy_zzz[k] = -2.0 * to_y_yzzz[k] * tke_0 + 4.0 * to_xxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : DF

        auto to_x_y_xz_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 20);

        auto to_x_y_xz_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 21);

        auto to_x_y_xz_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 22);

        auto to_x_y_xz_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 23);

        auto to_x_y_xz_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 24);

        auto to_x_y_xz_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 25);

        auto to_x_y_xz_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 26);

        auto to_x_y_xz_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 27);

        auto to_x_y_xz_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 28);

        auto to_x_y_xz_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_y_xz_xxx, to_x_y_xz_xxy, to_x_y_xz_xxz, to_x_y_xz_xyy, to_x_y_xz_xyz, to_x_y_xz_xzz, to_x_y_xz_yyy, to_x_y_xz_yyz, to_x_y_xz_yzz, to_x_y_xz_zzz, to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_z_xx, to_z_xxxy, to_z_xxyy, to_z_xxyz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_yy, to_z_yyyy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xz_xxx[k] = -2.0 * to_z_xxxy[k] * tke_0 + 4.0 * to_xxz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xz_xxy[k] = to_z_xx[k] - 2.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xxz[k] = -2.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xz_xyy[k] = 2.0 * to_z_xy[k] - 2.0 * to_z_xyyy[k] * tke_0 - 4.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xz_xyz[k] = to_z_xz[k] - 2.0 * to_z_xyyz[k] * tke_0 - 2.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xz_xzz[k] = -2.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xz_yyy[k] = 3.0 * to_z_yy[k] - 2.0 * to_z_yyyy[k] * tke_0 - 6.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xz_yyz[k] = 2.0 * to_z_yz[k] - 2.0 * to_z_yyyz[k] * tke_0 - 4.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xz_yzz[k] = to_z_zz[k] - 2.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xz_zzz[k] = -2.0 * to_z_yzzz[k] * tke_0 + 4.0 * to_xxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-100 components of targeted buffer : DF

        auto to_x_y_yy_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 30);

        auto to_x_y_yy_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 31);

        auto to_x_y_yy_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 32);

        auto to_x_y_yy_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 33);

        auto to_x_y_yy_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 34);

        auto to_x_y_yy_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 35);

        auto to_x_y_yy_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 36);

        auto to_x_y_yy_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 37);

        auto to_x_y_yy_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 38);

        auto to_x_y_yy_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_x_y_yy_xxx, to_x_y_yy_xxy, to_x_y_yy_xxz, to_x_y_yy_xyy, to_x_y_yy_xyz, to_x_y_yy_xzz, to_x_y_yy_yyy, to_x_y_yy_yyz, to_x_y_yy_yzz, to_x_y_yy_zzz, to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yy_xxx[k] = 4.0 * to_xyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yy_xxy[k] = -2.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xxz[k] = 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yy_xyy[k] = -4.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yy_xyz[k] = -2.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yy_xzz[k] = 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yy_yyy[k] = -6.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yy_yyz[k] = -4.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yy_yzz[k] = -2.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yy_zzz[k] = 4.0 * to_xyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 100-110 components of targeted buffer : DF

        auto to_x_y_yz_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 40);

        auto to_x_y_yz_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 41);

        auto to_x_y_yz_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 42);

        auto to_x_y_yz_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 43);

        auto to_x_y_yz_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 44);

        auto to_x_y_yz_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 45);

        auto to_x_y_yz_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 46);

        auto to_x_y_yz_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 47);

        auto to_x_y_yz_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 48);

        auto to_x_y_yz_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_x_y_yz_xxx, to_x_y_yz_xxy, to_x_y_yz_xxz, to_x_y_yz_xyy, to_x_y_yz_xyz, to_x_y_yz_xzz, to_x_y_yz_yyy, to_x_y_yz_yyz, to_x_y_yz_yzz, to_x_y_yz_zzz, to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yz_xxx[k] = 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yz_xxy[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xxz[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yz_xyy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yz_xyz[k] = -2.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yz_xzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yz_yyy[k] = -6.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yz_yyz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yz_yzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yz_zzz[k] = 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 110-120 components of targeted buffer : DF

        auto to_x_y_zz_xxx = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 50);

        auto to_x_y_zz_xxy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 51);

        auto to_x_y_zz_xxz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 52);

        auto to_x_y_zz_xyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 53);

        auto to_x_y_zz_xyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 54);

        auto to_x_y_zz_xzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 55);

        auto to_x_y_zz_yyy = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 56);

        auto to_x_y_zz_yyz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 57);

        auto to_x_y_zz_yzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 58);

        auto to_x_y_zz_zzz = pbuffer.data(idx_op_geom_101_df + 1 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_y_zz_xxx, to_x_y_zz_xxy, to_x_y_zz_xxz, to_x_y_zz_xyy, to_x_y_zz_xyz, to_x_y_zz_xzz, to_x_y_zz_yyy, to_x_y_zz_yyz, to_x_y_zz_yzz, to_x_y_zz_zzz, to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zz_xxx[k] = 4.0 * to_xzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_zz_xxy[k] = -2.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xxz[k] = 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_zz_xyy[k] = -4.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_zz_xyz[k] = -2.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_zz_xzz[k] = 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_zz_yyy[k] = -6.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_zz_yyz[k] = -4.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_zz_yzz[k] = -2.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_zz_zzz[k] = 4.0 * to_xzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-130 components of targeted buffer : DF

        auto to_x_z_xx_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 0);

        auto to_x_z_xx_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 1);

        auto to_x_z_xx_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 2);

        auto to_x_z_xx_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 3);

        auto to_x_z_xx_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 4);

        auto to_x_z_xx_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 5);

        auto to_x_z_xx_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 6);

        auto to_x_z_xx_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 7);

        auto to_x_z_xx_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 8);

        auto to_x_z_xx_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_x_xx, to_x_xxxz, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_z_xx_xxx, to_x_z_xx_xxy, to_x_z_xx_xxz, to_x_z_xx_xyy, to_x_z_xx_xyz, to_x_z_xx_xzz, to_x_z_xx_yyy, to_x_z_xx_yyz, to_x_z_xx_yzz, to_x_z_xx_zzz, to_x_zz, to_x_zzzz, to_xxx_xx, to_xxx_xxxz, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xx_xxx[k] = -4.0 * to_x_xxxz[k] * tke_0 + 4.0 * to_xxx_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxy[k] = -4.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xxx_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xxz[k] = 2.0 * to_x_xx[k] - 4.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xxx_xx[k] * tbe_0 + 4.0 * to_xxx_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xyy[k] = -4.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xxx_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xx_xyz[k] = 2.0 * to_x_xy[k] - 4.0 * to_x_xyzz[k] * tke_0 - 2.0 * to_xxx_xy[k] * tbe_0 + 4.0 * to_xxx_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xx_xzz[k] = 4.0 * to_x_xz[k] - 4.0 * to_x_xzzz[k] * tke_0 - 4.0 * to_xxx_xz[k] * tbe_0 + 4.0 * to_xxx_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yyy[k] = -4.0 * to_x_yyyz[k] * tke_0 + 4.0 * to_xxx_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xx_yyz[k] = 2.0 * to_x_yy[k] - 4.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xxx_yy[k] * tbe_0 + 4.0 * to_xxx_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xx_yzz[k] = 4.0 * to_x_yz[k] - 4.0 * to_x_yzzz[k] * tke_0 - 4.0 * to_xxx_yz[k] * tbe_0 + 4.0 * to_xxx_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xx_zzz[k] = 6.0 * to_x_zz[k] - 4.0 * to_x_zzzz[k] * tke_0 - 6.0 * to_xxx_zz[k] * tbe_0 + 4.0 * to_xxx_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 130-140 components of targeted buffer : DF

        auto to_x_z_xy_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 10);

        auto to_x_z_xy_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 11);

        auto to_x_z_xy_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 12);

        auto to_x_z_xy_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 13);

        auto to_x_z_xy_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 14);

        auto to_x_z_xy_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 15);

        auto to_x_z_xy_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 16);

        auto to_x_z_xy_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 17);

        auto to_x_z_xy_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 18);

        auto to_x_z_xy_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_z_xy_xxx, to_x_z_xy_xxy, to_x_z_xy_xxz, to_x_z_xy_xyy, to_x_z_xy_xyz, to_x_z_xy_xzz, to_x_z_xy_yyy, to_x_z_xy_yyz, to_x_z_xy_yzz, to_x_z_xy_zzz, to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, to_y_xx, to_y_xxxz, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, to_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xy_xxx[k] = -2.0 * to_y_xxxz[k] * tke_0 + 4.0 * to_xxy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxy[k] = -2.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xxz[k] = to_y_xx[k] - 2.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xyy[k] = -2.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xy_xyz[k] = to_y_xy[k] - 2.0 * to_y_xyzz[k] * tke_0 - 2.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xy_xzz[k] = 2.0 * to_y_xz[k] - 2.0 * to_y_xzzz[k] * tke_0 - 4.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yyy[k] = -2.0 * to_y_yyyz[k] * tke_0 + 4.0 * to_xxy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xy_yyz[k] = to_y_yy[k] - 2.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xy_yzz[k] = 2.0 * to_y_yz[k] - 2.0 * to_y_yzzz[k] * tke_0 - 4.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xy_zzz[k] = 3.0 * to_y_zz[k] - 2.0 * to_y_zzzz[k] * tke_0 - 6.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 140-150 components of targeted buffer : DF

        auto to_x_z_xz_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 20);

        auto to_x_z_xz_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 21);

        auto to_x_z_xz_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 22);

        auto to_x_z_xz_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 23);

        auto to_x_z_xz_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 24);

        auto to_x_z_xz_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 25);

        auto to_x_z_xz_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 26);

        auto to_x_z_xz_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 27);

        auto to_x_z_xz_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 28);

        auto to_x_z_xz_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_z_xz_xxx, to_x_z_xz_xxy, to_x_z_xz_xxz, to_x_z_xz_xyy, to_x_z_xz_xyz, to_x_z_xz_xzz, to_x_z_xz_yyy, to_x_z_xz_yyz, to_x_z_xz_yzz, to_x_z_xz_zzz, to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, to_z_xx, to_z_xxxz, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, to_z_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xz_xxx[k] = -2.0 * to_z_xxxz[k] * tke_0 + 4.0 * to_xxz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxy[k] = -2.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xxz[k] = to_z_xx[k] - 2.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xyy[k] = -2.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xz_xyz[k] = to_z_xy[k] - 2.0 * to_z_xyzz[k] * tke_0 - 2.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xz_xzz[k] = 2.0 * to_z_xz[k] - 2.0 * to_z_xzzz[k] * tke_0 - 4.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yyy[k] = -2.0 * to_z_yyyz[k] * tke_0 + 4.0 * to_xxz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xz_yyz[k] = to_z_yy[k] - 2.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xz_yzz[k] = 2.0 * to_z_yz[k] - 2.0 * to_z_yzzz[k] * tke_0 - 4.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xz_zzz[k] = 3.0 * to_z_zz[k] - 2.0 * to_z_zzzz[k] * tke_0 - 6.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-160 components of targeted buffer : DF

        auto to_x_z_yy_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 30);

        auto to_x_z_yy_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 31);

        auto to_x_z_yy_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 32);

        auto to_x_z_yy_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 33);

        auto to_x_z_yy_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 34);

        auto to_x_z_yy_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 35);

        auto to_x_z_yy_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 36);

        auto to_x_z_yy_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 37);

        auto to_x_z_yy_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 38);

        auto to_x_z_yy_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_x_z_yy_xxx, to_x_z_yy_xxy, to_x_z_yy_xxz, to_x_z_yy_xyy, to_x_z_yy_xyz, to_x_z_yy_xzz, to_x_z_yy_yyy, to_x_z_yy_yyz, to_x_z_yy_yzz, to_x_z_yy_zzz, to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yy_xxx[k] = 4.0 * to_xyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxy[k] = 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xxz[k] = -2.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xyy[k] = 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yy_xyz[k] = -2.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yy_xzz[k] = -4.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yyy[k] = 4.0 * to_xyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yy_yyz[k] = -2.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yy_yzz[k] = -4.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yy_zzz[k] = -6.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 160-170 components of targeted buffer : DF

        auto to_x_z_yz_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 40);

        auto to_x_z_yz_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 41);

        auto to_x_z_yz_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 42);

        auto to_x_z_yz_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 43);

        auto to_x_z_yz_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 44);

        auto to_x_z_yz_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 45);

        auto to_x_z_yz_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 46);

        auto to_x_z_yz_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 47);

        auto to_x_z_yz_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 48);

        auto to_x_z_yz_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_x_z_yz_xxx, to_x_z_yz_xxy, to_x_z_yz_xxz, to_x_z_yz_xyy, to_x_z_yz_xyz, to_x_z_yz_xzz, to_x_z_yz_yyy, to_x_z_yz_yyz, to_x_z_yz_yzz, to_x_z_yz_zzz, to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yz_xxx[k] = 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxy[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xxz[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xyy[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yz_xyz[k] = -2.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yz_xzz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yyy[k] = 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yz_yyz[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yz_yzz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yz_zzz[k] = -6.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 170-180 components of targeted buffer : DF

        auto to_x_z_zz_xxx = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 50);

        auto to_x_z_zz_xxy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 51);

        auto to_x_z_zz_xxz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 52);

        auto to_x_z_zz_xyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 53);

        auto to_x_z_zz_xyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 54);

        auto to_x_z_zz_xzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 55);

        auto to_x_z_zz_yyy = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 56);

        auto to_x_z_zz_yyz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 57);

        auto to_x_z_zz_yzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 58);

        auto to_x_z_zz_zzz = pbuffer.data(idx_op_geom_101_df + 2 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_z_zz_xxx, to_x_z_zz_xxy, to_x_z_zz_xxz, to_x_z_zz_xyy, to_x_z_zz_xyz, to_x_z_zz_xzz, to_x_z_zz_yyy, to_x_z_zz_yyz, to_x_z_zz_yzz, to_x_z_zz_zzz, to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zz_xxx[k] = 4.0 * to_xzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxy[k] = 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xxz[k] = -2.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xyy[k] = 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_zz_xyz[k] = -2.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_zz_xzz[k] = -4.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yyy[k] = 4.0 * to_xzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_zz_yyz[k] = -2.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_zz_yzz[k] = -4.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_zz_zzz[k] = -6.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-190 components of targeted buffer : DF

        auto to_y_x_xx_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 0);

        auto to_y_x_xx_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 1);

        auto to_y_x_xx_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 2);

        auto to_y_x_xx_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 3);

        auto to_y_x_xx_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 4);

        auto to_y_x_xx_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 5);

        auto to_y_x_xx_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 6);

        auto to_y_x_xx_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 7);

        auto to_y_x_xx_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 8);

        auto to_y_x_xx_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, to_y_x_xx_xxx, to_y_x_xx_xxy, to_y_x_xx_xxz, to_y_x_xx_xyy, to_y_x_xx_xyz, to_y_x_xx_xzz, to_y_x_xx_yyy, to_y_x_xx_yyz, to_y_x_xx_yzz, to_y_x_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xx_xxx[k] = -6.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xx_xxy[k] = -4.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xx_xxz[k] = -4.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xx_xyy[k] = -2.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xx_xyz[k] = -2.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xx_xzz[k] = -2.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xx_yyy[k] = 4.0 * to_xxy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xx_yyz[k] = 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xx_yzz[k] = 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xx_zzz[k] = 4.0 * to_xxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 190-200 components of targeted buffer : DF

        auto to_y_x_xy_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 10);

        auto to_y_x_xy_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 11);

        auto to_y_x_xy_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 12);

        auto to_y_x_xy_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 13);

        auto to_y_x_xy_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 14);

        auto to_y_x_xy_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 15);

        auto to_y_x_xy_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 16);

        auto to_y_x_xy_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 17);

        auto to_y_x_xy_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 18);

        auto to_y_x_xy_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_xx, to_x_xxxx, to_x_xxxy, to_x_xxxz, to_x_xxyy, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yz, to_x_zz, to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, to_y_x_xy_xxx, to_y_x_xy_xxy, to_y_x_xy_xxz, to_y_x_xy_xyy, to_y_x_xy_xyz, to_y_x_xy_xzz, to_y_x_xy_yyy, to_y_x_xy_yyz, to_y_x_xy_yzz, to_y_x_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xy_xxx[k] = 3.0 * to_x_xx[k] - 2.0 * to_x_xxxx[k] * tke_0 - 6.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xy_xxy[k] = 2.0 * to_x_xy[k] - 2.0 * to_x_xxxy[k] * tke_0 - 4.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xy_xxz[k] = 2.0 * to_x_xz[k] - 2.0 * to_x_xxxz[k] * tke_0 - 4.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xy_xyy[k] = to_x_yy[k] - 2.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xy_xyz[k] = to_x_yz[k] - 2.0 * to_x_xxyz[k] * tke_0 - 2.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xy_xzz[k] = to_x_zz[k] - 2.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xy_yyy[k] = -2.0 * to_x_xyyy[k] * tke_0 + 4.0 * to_xyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xy_yyz[k] = -2.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xy_yzz[k] = -2.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xy_zzz[k] = -2.0 * to_x_xzzz[k] * tke_0 + 4.0 * to_xyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 200-210 components of targeted buffer : DF

        auto to_y_x_xz_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 20);

        auto to_y_x_xz_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 21);

        auto to_y_x_xz_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 22);

        auto to_y_x_xz_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 23);

        auto to_y_x_xz_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 24);

        auto to_y_x_xz_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 25);

        auto to_y_x_xz_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 26);

        auto to_y_x_xz_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 27);

        auto to_y_x_xz_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 28);

        auto to_y_x_xz_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, to_y_x_xz_xxx, to_y_x_xz_xxy, to_y_x_xz_xxz, to_y_x_xz_xyy, to_y_x_xz_xyz, to_y_x_xz_xzz, to_y_x_xz_yyy, to_y_x_xz_yyz, to_y_x_xz_yzz, to_y_x_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xz_xxx[k] = -6.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xz_xxy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xz_xxz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xz_xyy[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xz_xyz[k] = -2.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xz_xzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xz_yyy[k] = 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xz_yyz[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xz_yzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xz_zzz[k] = 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-220 components of targeted buffer : DF

        auto to_y_x_yy_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 30);

        auto to_y_x_yy_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 31);

        auto to_y_x_yy_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 32);

        auto to_y_x_yy_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 33);

        auto to_y_x_yy_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 34);

        auto to_y_x_yy_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 35);

        auto to_y_x_yy_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 36);

        auto to_y_x_yy_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 37);

        auto to_y_x_yy_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 38);

        auto to_y_x_yy_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_y_x_yy_xxx, to_y_x_yy_xxy, to_y_x_yy_xxz, to_y_x_yy_xyy, to_y_x_yy_xyz, to_y_x_yy_xzz, to_y_x_yy_yyy, to_y_x_yy_yyz, to_y_x_yy_yzz, to_y_x_yy_zzz, to_y_xx, to_y_xxxx, to_y_xxxy, to_y_xxxz, to_y_xxyy, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yz, to_y_zz, to_yyy_xx, to_yyy_xxxx, to_yyy_xxxy, to_yyy_xxxz, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yy_xxx[k] = 6.0 * to_y_xx[k] - 4.0 * to_y_xxxx[k] * tke_0 - 6.0 * to_yyy_xx[k] * tbe_0 + 4.0 * to_yyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yy_xxy[k] = 4.0 * to_y_xy[k] - 4.0 * to_y_xxxy[k] * tke_0 - 4.0 * to_yyy_xy[k] * tbe_0 + 4.0 * to_yyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yy_xxz[k] = 4.0 * to_y_xz[k] - 4.0 * to_y_xxxz[k] * tke_0 - 4.0 * to_yyy_xz[k] * tbe_0 + 4.0 * to_yyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yy_xyy[k] = 2.0 * to_y_yy[k] - 4.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_yyy_yy[k] * tbe_0 + 4.0 * to_yyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yy_xyz[k] = 2.0 * to_y_yz[k] - 4.0 * to_y_xxyz[k] * tke_0 - 2.0 * to_yyy_yz[k] * tbe_0 + 4.0 * to_yyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yy_xzz[k] = 2.0 * to_y_zz[k] - 4.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_yyy_zz[k] * tbe_0 + 4.0 * to_yyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yy_yyy[k] = -4.0 * to_y_xyyy[k] * tke_0 + 4.0 * to_yyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yy_yyz[k] = -4.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_yyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yy_yzz[k] = -4.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_yyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yy_zzz[k] = -4.0 * to_y_xzzz[k] * tke_0 + 4.0 * to_yyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 220-230 components of targeted buffer : DF

        auto to_y_x_yz_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 40);

        auto to_y_x_yz_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 41);

        auto to_y_x_yz_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 42);

        auto to_y_x_yz_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 43);

        auto to_y_x_yz_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 44);

        auto to_y_x_yz_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 45);

        auto to_y_x_yz_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 46);

        auto to_y_x_yz_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 47);

        auto to_y_x_yz_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 48);

        auto to_y_x_yz_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_x_yz_xxx, to_y_x_yz_xxy, to_y_x_yz_xxz, to_y_x_yz_xyy, to_y_x_yz_xyz, to_y_x_yz_xzz, to_y_x_yz_yyy, to_y_x_yz_yyz, to_y_x_yz_yzz, to_y_x_yz_zzz, to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, to_z_xx, to_z_xxxx, to_z_xxxy, to_z_xxxz, to_z_xxyy, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yz_xxx[k] = 3.0 * to_z_xx[k] - 2.0 * to_z_xxxx[k] * tke_0 - 6.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yz_xxy[k] = 2.0 * to_z_xy[k] - 2.0 * to_z_xxxy[k] * tke_0 - 4.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yz_xxz[k] = 2.0 * to_z_xz[k] - 2.0 * to_z_xxxz[k] * tke_0 - 4.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yz_xyy[k] = to_z_yy[k] - 2.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yz_xyz[k] = to_z_yz[k] - 2.0 * to_z_xxyz[k] * tke_0 - 2.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yz_xzz[k] = to_z_zz[k] - 2.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yz_yyy[k] = -2.0 * to_z_xyyy[k] * tke_0 + 4.0 * to_yyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yz_yyz[k] = -2.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yz_yzz[k] = -2.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yz_zzz[k] = -2.0 * to_z_xzzz[k] * tke_0 + 4.0 * to_yyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 230-240 components of targeted buffer : DF

        auto to_y_x_zz_xxx = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 50);

        auto to_y_x_zz_xxy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 51);

        auto to_y_x_zz_xxz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 52);

        auto to_y_x_zz_xyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 53);

        auto to_y_x_zz_xyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 54);

        auto to_y_x_zz_xzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 55);

        auto to_y_x_zz_yyy = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 56);

        auto to_y_x_zz_yyz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 57);

        auto to_y_x_zz_yzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 58);

        auto to_y_x_zz_zzz = pbuffer.data(idx_op_geom_101_df + 3 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_x_zz_xxx, to_y_x_zz_xxy, to_y_x_zz_xxz, to_y_x_zz_xyy, to_y_x_zz_xyz, to_y_x_zz_xzz, to_y_x_zz_yyy, to_y_x_zz_yyz, to_y_x_zz_yzz, to_y_x_zz_zzz, to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zz_xxx[k] = -6.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_zz_xxy[k] = -4.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_zz_xxz[k] = -4.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_zz_xyy[k] = -2.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_zz_xyz[k] = -2.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_zz_xzz[k] = -2.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_zz_yyy[k] = 4.0 * to_yzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_zz_yyz[k] = 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_zz_yzz[k] = 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_zz_zzz[k] = 4.0 * to_yzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-250 components of targeted buffer : DF

        auto to_y_y_xx_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 0);

        auto to_y_y_xx_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 1);

        auto to_y_y_xx_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 2);

        auto to_y_y_xx_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 3);

        auto to_y_y_xx_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 4);

        auto to_y_y_xx_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 5);

        auto to_y_y_xx_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 6);

        auto to_y_y_xx_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 7);

        auto to_y_y_xx_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 8);

        auto to_y_y_xx_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_y_y_xx_xxx, to_y_y_xx_xxy, to_y_y_xx_xxz, to_y_y_xx_xyy, to_y_y_xx_xyz, to_y_y_xx_xzz, to_y_y_xx_yyy, to_y_y_xx_yyz, to_y_y_xx_yzz, to_y_y_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xx_xxx[k] = 4.0 * to_xxy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xx_xxy[k] = -2.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xxz[k] = 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xx_xyy[k] = -4.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xx_xyz[k] = -2.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xx_xzz[k] = 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xx_yyy[k] = -6.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xx_yyz[k] = -4.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xx_yzz[k] = -2.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xx_zzz[k] = 4.0 * to_xxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 250-260 components of targeted buffer : DF

        auto to_y_y_xy_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 10);

        auto to_y_y_xy_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 11);

        auto to_y_y_xy_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 12);

        auto to_y_y_xy_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 13);

        auto to_y_y_xy_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 14);

        auto to_y_y_xy_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 15);

        auto to_y_y_xy_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 16);

        auto to_y_y_xy_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 17);

        auto to_y_y_xy_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 18);

        auto to_y_y_xy_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_xx, to_x_xxxy, to_x_xxyy, to_x_xxyz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_yy, to_x_yyyy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_y_y_xy_xxx, to_y_y_xy_xxy, to_y_y_xy_xxz, to_y_y_xy_xyy, to_y_y_xy_xyz, to_y_y_xy_xzz, to_y_y_xy_yyy, to_y_y_xy_yyz, to_y_y_xy_yzz, to_y_y_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xy_xxx[k] = -2.0 * to_x_xxxy[k] * tke_0 + 4.0 * to_xyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xy_xxy[k] = to_x_xx[k] - 2.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xxz[k] = -2.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xy_xyy[k] = 2.0 * to_x_xy[k] - 2.0 * to_x_xyyy[k] * tke_0 - 4.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xy_xyz[k] = to_x_xz[k] - 2.0 * to_x_xyyz[k] * tke_0 - 2.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xy_xzz[k] = -2.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xy_yyy[k] = 3.0 * to_x_yy[k] - 2.0 * to_x_yyyy[k] * tke_0 - 6.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xy_yyz[k] = 2.0 * to_x_yz[k] - 2.0 * to_x_yyyz[k] * tke_0 - 4.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xy_yzz[k] = to_x_zz[k] - 2.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xy_zzz[k] = -2.0 * to_x_yzzz[k] * tke_0 + 4.0 * to_xyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 260-270 components of targeted buffer : DF

        auto to_y_y_xz_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 20);

        auto to_y_y_xz_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 21);

        auto to_y_y_xz_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 22);

        auto to_y_y_xz_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 23);

        auto to_y_y_xz_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 24);

        auto to_y_y_xz_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 25);

        auto to_y_y_xz_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 26);

        auto to_y_y_xz_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 27);

        auto to_y_y_xz_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 28);

        auto to_y_y_xz_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_y_y_xz_xxx, to_y_y_xz_xxy, to_y_y_xz_xxz, to_y_y_xz_xyy, to_y_y_xz_xyz, to_y_y_xz_xzz, to_y_y_xz_yyy, to_y_y_xz_yyz, to_y_y_xz_yzz, to_y_y_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xz_xxx[k] = 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xz_xxy[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xxz[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xz_xyy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xz_xyz[k] = -2.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xz_xzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xz_yyy[k] = -6.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xz_yyz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xz_yzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xz_zzz[k] = 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-280 components of targeted buffer : DF

        auto to_y_y_yy_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 30);

        auto to_y_y_yy_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 31);

        auto to_y_y_yy_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 32);

        auto to_y_y_yy_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 33);

        auto to_y_y_yy_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 34);

        auto to_y_y_yy_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 35);

        auto to_y_y_yy_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 36);

        auto to_y_y_yy_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 37);

        auto to_y_y_yy_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 38);

        auto to_y_y_yy_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_y_xx, to_y_xxxy, to_y_xxyy, to_y_xxyz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_y_yy_xxx, to_y_y_yy_xxy, to_y_y_yy_xxz, to_y_y_yy_xyy, to_y_y_yy_xyz, to_y_y_yy_xzz, to_y_y_yy_yyy, to_y_y_yy_yyz, to_y_y_yy_yzz, to_y_y_yy_zzz, to_y_yy, to_y_yyyy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, to_yyy_xx, to_yyy_xxxy, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_yy, to_yyy_yyyy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yy_xxx[k] = -4.0 * to_y_xxxy[k] * tke_0 + 4.0 * to_yyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yy_xxy[k] = 2.0 * to_y_xx[k] - 4.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_yyy_xx[k] * tbe_0 + 4.0 * to_yyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xxz[k] = -4.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_yyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yy_xyy[k] = 4.0 * to_y_xy[k] - 4.0 * to_y_xyyy[k] * tke_0 - 4.0 * to_yyy_xy[k] * tbe_0 + 4.0 * to_yyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yy_xyz[k] = 2.0 * to_y_xz[k] - 4.0 * to_y_xyyz[k] * tke_0 - 2.0 * to_yyy_xz[k] * tbe_0 + 4.0 * to_yyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yy_xzz[k] = -4.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_yyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yy_yyy[k] = 6.0 * to_y_yy[k] - 4.0 * to_y_yyyy[k] * tke_0 - 6.0 * to_yyy_yy[k] * tbe_0 + 4.0 * to_yyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yy_yyz[k] = 4.0 * to_y_yz[k] - 4.0 * to_y_yyyz[k] * tke_0 - 4.0 * to_yyy_yz[k] * tbe_0 + 4.0 * to_yyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yy_yzz[k] = 2.0 * to_y_zz[k] - 4.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_yyy_zz[k] * tbe_0 + 4.0 * to_yyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yy_zzz[k] = -4.0 * to_y_yzzz[k] * tke_0 + 4.0 * to_yyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 280-290 components of targeted buffer : DF

        auto to_y_y_yz_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 40);

        auto to_y_y_yz_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 41);

        auto to_y_y_yz_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 42);

        auto to_y_y_yz_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 43);

        auto to_y_y_yz_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 44);

        auto to_y_y_yz_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 45);

        auto to_y_y_yz_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 46);

        auto to_y_y_yz_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 47);

        auto to_y_y_yz_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 48);

        auto to_y_y_yz_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_y_yz_xxx, to_y_y_yz_xxy, to_y_y_yz_xxz, to_y_y_yz_xyy, to_y_y_yz_xyz, to_y_y_yz_xzz, to_y_y_yz_yyy, to_y_y_yz_yyz, to_y_y_yz_yzz, to_y_y_yz_zzz, to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_z_xx, to_z_xxxy, to_z_xxyy, to_z_xxyz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_yy, to_z_yyyy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yz_xxx[k] = -2.0 * to_z_xxxy[k] * tke_0 + 4.0 * to_yyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yz_xxy[k] = to_z_xx[k] - 2.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xxz[k] = -2.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yz_xyy[k] = 2.0 * to_z_xy[k] - 2.0 * to_z_xyyy[k] * tke_0 - 4.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yz_xyz[k] = to_z_xz[k] - 2.0 * to_z_xyyz[k] * tke_0 - 2.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yz_xzz[k] = -2.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yz_yyy[k] = 3.0 * to_z_yy[k] - 2.0 * to_z_yyyy[k] * tke_0 - 6.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yz_yyz[k] = 2.0 * to_z_yz[k] - 2.0 * to_z_yyyz[k] * tke_0 - 4.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yz_yzz[k] = to_z_zz[k] - 2.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yz_zzz[k] = -2.0 * to_z_yzzz[k] * tke_0 + 4.0 * to_yyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 290-300 components of targeted buffer : DF

        auto to_y_y_zz_xxx = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 50);

        auto to_y_y_zz_xxy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 51);

        auto to_y_y_zz_xxz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 52);

        auto to_y_y_zz_xyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 53);

        auto to_y_y_zz_xyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 54);

        auto to_y_y_zz_xzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 55);

        auto to_y_y_zz_yyy = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 56);

        auto to_y_y_zz_yyz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 57);

        auto to_y_y_zz_yzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 58);

        auto to_y_y_zz_zzz = pbuffer.data(idx_op_geom_101_df + 4 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_y_zz_xxx, to_y_y_zz_xxy, to_y_y_zz_xxz, to_y_y_zz_xyy, to_y_y_zz_xyz, to_y_y_zz_xzz, to_y_y_zz_yyy, to_y_y_zz_yyz, to_y_y_zz_yzz, to_y_y_zz_zzz, to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zz_xxx[k] = 4.0 * to_yzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_zz_xxy[k] = -2.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xxz[k] = 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_zz_xyy[k] = -4.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_zz_xyz[k] = -2.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_zz_xzz[k] = 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_zz_yyy[k] = -6.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_zz_yyz[k] = -4.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_zz_yzz[k] = -2.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_zz_zzz[k] = 4.0 * to_yzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-310 components of targeted buffer : DF

        auto to_y_z_xx_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 0);

        auto to_y_z_xx_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 1);

        auto to_y_z_xx_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 2);

        auto to_y_z_xx_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 3);

        auto to_y_z_xx_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 4);

        auto to_y_z_xx_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 5);

        auto to_y_z_xx_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 6);

        auto to_y_z_xx_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 7);

        auto to_y_z_xx_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 8);

        auto to_y_z_xx_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, to_y_z_xx_xxx, to_y_z_xx_xxy, to_y_z_xx_xxz, to_y_z_xx_xyy, to_y_z_xx_xyz, to_y_z_xx_xzz, to_y_z_xx_yyy, to_y_z_xx_yyz, to_y_z_xx_yzz, to_y_z_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xx_xxx[k] = 4.0 * to_xxy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxy[k] = 4.0 * to_xxy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xxz[k] = -2.0 * to_xxy_xx[k] * tbe_0 + 4.0 * to_xxy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xyy[k] = 4.0 * to_xxy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xx_xyz[k] = -2.0 * to_xxy_xy[k] * tbe_0 + 4.0 * to_xxy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xx_xzz[k] = -4.0 * to_xxy_xz[k] * tbe_0 + 4.0 * to_xxy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yyy[k] = 4.0 * to_xxy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xx_yyz[k] = -2.0 * to_xxy_yy[k] * tbe_0 + 4.0 * to_xxy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xx_yzz[k] = -4.0 * to_xxy_yz[k] * tbe_0 + 4.0 * to_xxy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xx_zzz[k] = -6.0 * to_xxy_zz[k] * tbe_0 + 4.0 * to_xxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 310-320 components of targeted buffer : DF

        auto to_y_z_xy_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 10);

        auto to_y_z_xy_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 11);

        auto to_y_z_xy_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 12);

        auto to_y_z_xy_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 13);

        auto to_y_z_xy_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 14);

        auto to_y_z_xy_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 15);

        auto to_y_z_xy_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 16);

        auto to_y_z_xy_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 17);

        auto to_y_z_xy_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 18);

        auto to_y_z_xy_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_x_xx, to_x_xxxz, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_x_zzzz, to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, to_y_z_xy_xxx, to_y_z_xy_xxy, to_y_z_xy_xxz, to_y_z_xy_xyy, to_y_z_xy_xyz, to_y_z_xy_xzz, to_y_z_xy_yyy, to_y_z_xy_yyz, to_y_z_xy_yzz, to_y_z_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xy_xxx[k] = -2.0 * to_x_xxxz[k] * tke_0 + 4.0 * to_xyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxy[k] = -2.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xxz[k] = to_x_xx[k] - 2.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xyy_xx[k] * tbe_0 + 4.0 * to_xyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xyy[k] = -2.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xy_xyz[k] = to_x_xy[k] - 2.0 * to_x_xyzz[k] * tke_0 - 2.0 * to_xyy_xy[k] * tbe_0 + 4.0 * to_xyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xy_xzz[k] = 2.0 * to_x_xz[k] - 2.0 * to_x_xzzz[k] * tke_0 - 4.0 * to_xyy_xz[k] * tbe_0 + 4.0 * to_xyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yyy[k] = -2.0 * to_x_yyyz[k] * tke_0 + 4.0 * to_xyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xy_yyz[k] = to_x_yy[k] - 2.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xyy_yy[k] * tbe_0 + 4.0 * to_xyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xy_yzz[k] = 2.0 * to_x_yz[k] - 2.0 * to_x_yzzz[k] * tke_0 - 4.0 * to_xyy_yz[k] * tbe_0 + 4.0 * to_xyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xy_zzz[k] = 3.0 * to_x_zz[k] - 2.0 * to_x_zzzz[k] * tke_0 - 6.0 * to_xyy_zz[k] * tbe_0 + 4.0 * to_xyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 320-330 components of targeted buffer : DF

        auto to_y_z_xz_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 20);

        auto to_y_z_xz_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 21);

        auto to_y_z_xz_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 22);

        auto to_y_z_xz_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 23);

        auto to_y_z_xz_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 24);

        auto to_y_z_xz_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 25);

        auto to_y_z_xz_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 26);

        auto to_y_z_xz_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 27);

        auto to_y_z_xz_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 28);

        auto to_y_z_xz_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, to_y_z_xz_xxx, to_y_z_xz_xxy, to_y_z_xz_xxz, to_y_z_xz_xyy, to_y_z_xz_xyz, to_y_z_xz_xzz, to_y_z_xz_yyy, to_y_z_xz_yyz, to_y_z_xz_yzz, to_y_z_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xz_xxx[k] = 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxy[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xxz[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xyy[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xz_xyz[k] = -2.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xz_xzz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yyy[k] = 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xz_yyz[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xz_yzz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xz_zzz[k] = -6.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-340 components of targeted buffer : DF

        auto to_y_z_yy_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 30);

        auto to_y_z_yy_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 31);

        auto to_y_z_yy_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 32);

        auto to_y_z_yy_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 33);

        auto to_y_z_yy_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 34);

        auto to_y_z_yy_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 35);

        auto to_y_z_yy_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 36);

        auto to_y_z_yy_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 37);

        auto to_y_z_yy_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 38);

        auto to_y_z_yy_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_y_xx, to_y_xxxz, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_z_yy_xxx, to_y_z_yy_xxy, to_y_z_yy_xxz, to_y_z_yy_xyy, to_y_z_yy_xyz, to_y_z_yy_xzz, to_y_z_yy_yyy, to_y_z_yy_yyz, to_y_z_yy_yzz, to_y_z_yy_zzz, to_y_zz, to_y_zzzz, to_yyy_xx, to_yyy_xxxz, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yy_xxx[k] = -4.0 * to_y_xxxz[k] * tke_0 + 4.0 * to_yyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxy[k] = -4.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_yyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xxz[k] = 2.0 * to_y_xx[k] - 4.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_yyy_xx[k] * tbe_0 + 4.0 * to_yyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xyy[k] = -4.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_yyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yy_xyz[k] = 2.0 * to_y_xy[k] - 4.0 * to_y_xyzz[k] * tke_0 - 2.0 * to_yyy_xy[k] * tbe_0 + 4.0 * to_yyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yy_xzz[k] = 4.0 * to_y_xz[k] - 4.0 * to_y_xzzz[k] * tke_0 - 4.0 * to_yyy_xz[k] * tbe_0 + 4.0 * to_yyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yyy[k] = -4.0 * to_y_yyyz[k] * tke_0 + 4.0 * to_yyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yy_yyz[k] = 2.0 * to_y_yy[k] - 4.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_yyy_yy[k] * tbe_0 + 4.0 * to_yyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yy_yzz[k] = 4.0 * to_y_yz[k] - 4.0 * to_y_yzzz[k] * tke_0 - 4.0 * to_yyy_yz[k] * tbe_0 + 4.0 * to_yyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yy_zzz[k] = 6.0 * to_y_zz[k] - 4.0 * to_y_zzzz[k] * tke_0 - 6.0 * to_yyy_zz[k] * tbe_0 + 4.0 * to_yyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 340-350 components of targeted buffer : DF

        auto to_y_z_yz_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 40);

        auto to_y_z_yz_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 41);

        auto to_y_z_yz_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 42);

        auto to_y_z_yz_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 43);

        auto to_y_z_yz_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 44);

        auto to_y_z_yz_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 45);

        auto to_y_z_yz_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 46);

        auto to_y_z_yz_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 47);

        auto to_y_z_yz_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 48);

        auto to_y_z_yz_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_z_yz_xxx, to_y_z_yz_xxy, to_y_z_yz_xxz, to_y_z_yz_xyy, to_y_z_yz_xyz, to_y_z_yz_xzz, to_y_z_yz_yyy, to_y_z_yz_yyz, to_y_z_yz_yzz, to_y_z_yz_zzz, to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, to_z_xx, to_z_xxxz, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, to_z_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yz_xxx[k] = -2.0 * to_z_xxxz[k] * tke_0 + 4.0 * to_yyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxy[k] = -2.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xxz[k] = to_z_xx[k] - 2.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xyy[k] = -2.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yz_xyz[k] = to_z_xy[k] - 2.0 * to_z_xyzz[k] * tke_0 - 2.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yz_xzz[k] = 2.0 * to_z_xz[k] - 2.0 * to_z_xzzz[k] * tke_0 - 4.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yyy[k] = -2.0 * to_z_yyyz[k] * tke_0 + 4.0 * to_yyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yz_yyz[k] = to_z_yy[k] - 2.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yz_yzz[k] = 2.0 * to_z_yz[k] - 2.0 * to_z_yzzz[k] * tke_0 - 4.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yz_zzz[k] = 3.0 * to_z_zz[k] - 2.0 * to_z_zzzz[k] * tke_0 - 6.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 350-360 components of targeted buffer : DF

        auto to_y_z_zz_xxx = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 50);

        auto to_y_z_zz_xxy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 51);

        auto to_y_z_zz_xxz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 52);

        auto to_y_z_zz_xyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 53);

        auto to_y_z_zz_xyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 54);

        auto to_y_z_zz_xzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 55);

        auto to_y_z_zz_yyy = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 56);

        auto to_y_z_zz_yyz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 57);

        auto to_y_z_zz_yzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 58);

        auto to_y_z_zz_zzz = pbuffer.data(idx_op_geom_101_df + 5 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_z_zz_xxx, to_y_z_zz_xxy, to_y_z_zz_xxz, to_y_z_zz_xyy, to_y_z_zz_xyz, to_y_z_zz_xzz, to_y_z_zz_yyy, to_y_z_zz_yyz, to_y_z_zz_yzz, to_y_z_zz_zzz, to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zz_xxx[k] = 4.0 * to_yzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxy[k] = 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xxz[k] = -2.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xyy[k] = 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_zz_xyz[k] = -2.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_zz_xzz[k] = -4.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yyy[k] = 4.0 * to_yzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_zz_yyz[k] = -2.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_zz_yzz[k] = -4.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_zz_zzz[k] = -6.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-370 components of targeted buffer : DF

        auto to_z_x_xx_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 0);

        auto to_z_x_xx_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 1);

        auto to_z_x_xx_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 2);

        auto to_z_x_xx_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 3);

        auto to_z_x_xx_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 4);

        auto to_z_x_xx_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 5);

        auto to_z_x_xx_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 6);

        auto to_z_x_xx_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 7);

        auto to_z_x_xx_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 8);

        auto to_z_x_xx_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, to_z_x_xx_xxx, to_z_x_xx_xxy, to_z_x_xx_xxz, to_z_x_xx_xyy, to_z_x_xx_xyz, to_z_x_xx_xzz, to_z_x_xx_yyy, to_z_x_xx_yyz, to_z_x_xx_yzz, to_z_x_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xx_xxx[k] = -6.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xx_xxy[k] = -4.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xx_xxz[k] = -4.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xx_xyy[k] = -2.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xx_xyz[k] = -2.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xx_xzz[k] = -2.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xx_yyy[k] = 4.0 * to_xxz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xx_yyz[k] = 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xx_yzz[k] = 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xx_zzz[k] = 4.0 * to_xxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 370-380 components of targeted buffer : DF

        auto to_z_x_xy_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 10);

        auto to_z_x_xy_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 11);

        auto to_z_x_xy_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 12);

        auto to_z_x_xy_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 13);

        auto to_z_x_xy_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 14);

        auto to_z_x_xy_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 15);

        auto to_z_x_xy_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 16);

        auto to_z_x_xy_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 17);

        auto to_z_x_xy_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 18);

        auto to_z_x_xy_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, to_z_x_xy_xxx, to_z_x_xy_xxy, to_z_x_xy_xxz, to_z_x_xy_xyy, to_z_x_xy_xyz, to_z_x_xy_xzz, to_z_x_xy_yyy, to_z_x_xy_yyz, to_z_x_xy_yzz, to_z_x_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xy_xxx[k] = -6.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xy_xxy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xy_xxz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xy_xyy[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xy_xyz[k] = -2.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xy_xzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xy_yyy[k] = 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xy_yyz[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xy_yzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xy_zzz[k] = 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 380-390 components of targeted buffer : DF

        auto to_z_x_xz_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 20);

        auto to_z_x_xz_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 21);

        auto to_z_x_xz_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 22);

        auto to_z_x_xz_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 23);

        auto to_z_x_xz_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 24);

        auto to_z_x_xz_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 25);

        auto to_z_x_xz_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 26);

        auto to_z_x_xz_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 27);

        auto to_z_x_xz_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 28);

        auto to_z_x_xz_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_xx, to_x_xxxx, to_x_xxxy, to_x_xxxz, to_x_xxyy, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yz, to_x_zz, to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, to_z_x_xz_xxx, to_z_x_xz_xxy, to_z_x_xz_xxz, to_z_x_xz_xyy, to_z_x_xz_xyz, to_z_x_xz_xzz, to_z_x_xz_yyy, to_z_x_xz_yyz, to_z_x_xz_yzz, to_z_x_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xz_xxx[k] = 3.0 * to_x_xx[k] - 2.0 * to_x_xxxx[k] * tke_0 - 6.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xz_xxy[k] = 2.0 * to_x_xy[k] - 2.0 * to_x_xxxy[k] * tke_0 - 4.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xz_xxz[k] = 2.0 * to_x_xz[k] - 2.0 * to_x_xxxz[k] * tke_0 - 4.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xz_xyy[k] = to_x_yy[k] - 2.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xz_xyz[k] = to_x_yz[k] - 2.0 * to_x_xxyz[k] * tke_0 - 2.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xz_xzz[k] = to_x_zz[k] - 2.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xz_yyy[k] = -2.0 * to_x_xyyy[k] * tke_0 + 4.0 * to_xzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xz_yyz[k] = -2.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xz_yzz[k] = -2.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xz_zzz[k] = -2.0 * to_x_xzzz[k] * tke_0 + 4.0 * to_xzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-400 components of targeted buffer : DF

        auto to_z_x_yy_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 30);

        auto to_z_x_yy_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 31);

        auto to_z_x_yy_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 32);

        auto to_z_x_yy_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 33);

        auto to_z_x_yy_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 34);

        auto to_z_x_yy_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 35);

        auto to_z_x_yy_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 36);

        auto to_z_x_yy_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 37);

        auto to_z_x_yy_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 38);

        auto to_z_x_yy_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, to_z_x_yy_xxx, to_z_x_yy_xxy, to_z_x_yy_xxz, to_z_x_yy_xyy, to_z_x_yy_xyz, to_z_x_yy_xzz, to_z_x_yy_yyy, to_z_x_yy_yyz, to_z_x_yy_yzz, to_z_x_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yy_xxx[k] = -6.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yy_xxy[k] = -4.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yy_xxz[k] = -4.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yy_xyy[k] = -2.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yy_xyz[k] = -2.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yy_xzz[k] = -2.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yy_yyy[k] = 4.0 * to_yyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yy_yyz[k] = 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yy_yzz[k] = 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yy_zzz[k] = 4.0 * to_yyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 400-410 components of targeted buffer : DF

        auto to_z_x_yz_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 40);

        auto to_z_x_yz_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 41);

        auto to_z_x_yz_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 42);

        auto to_z_x_yz_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 43);

        auto to_z_x_yz_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 44);

        auto to_z_x_yz_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 45);

        auto to_z_x_yz_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 46);

        auto to_z_x_yz_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 47);

        auto to_z_x_yz_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 48);

        auto to_z_x_yz_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_xx, to_y_xxxx, to_y_xxxy, to_y_xxxz, to_y_xxyy, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yz, to_y_zz, to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, to_z_x_yz_xxx, to_z_x_yz_xxy, to_z_x_yz_xxz, to_z_x_yz_xyy, to_z_x_yz_xyz, to_z_x_yz_xzz, to_z_x_yz_yyy, to_z_x_yz_yyz, to_z_x_yz_yzz, to_z_x_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yz_xxx[k] = 3.0 * to_y_xx[k] - 2.0 * to_y_xxxx[k] * tke_0 - 6.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yz_xxy[k] = 2.0 * to_y_xy[k] - 2.0 * to_y_xxxy[k] * tke_0 - 4.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yz_xxz[k] = 2.0 * to_y_xz[k] - 2.0 * to_y_xxxz[k] * tke_0 - 4.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yz_xyy[k] = to_y_yy[k] - 2.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yz_xyz[k] = to_y_yz[k] - 2.0 * to_y_xxyz[k] * tke_0 - 2.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yz_xzz[k] = to_y_zz[k] - 2.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yz_yyy[k] = -2.0 * to_y_xyyy[k] * tke_0 + 4.0 * to_yzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yz_yyz[k] = -2.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yz_yzz[k] = -2.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yz_zzz[k] = -2.0 * to_y_xzzz[k] * tke_0 + 4.0 * to_yzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 410-420 components of targeted buffer : DF

        auto to_z_x_zz_xxx = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 50);

        auto to_z_x_zz_xxy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 51);

        auto to_z_x_zz_xxz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 52);

        auto to_z_x_zz_xyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 53);

        auto to_z_x_zz_xyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 54);

        auto to_z_x_zz_xzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 55);

        auto to_z_x_zz_yyy = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 56);

        auto to_z_x_zz_yyz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 57);

        auto to_z_x_zz_yzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 58);

        auto to_z_x_zz_zzz = pbuffer.data(idx_op_geom_101_df + 6 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_x_zz_xxx, to_z_x_zz_xxy, to_z_x_zz_xxz, to_z_x_zz_xyy, to_z_x_zz_xyz, to_z_x_zz_xzz, to_z_x_zz_yyy, to_z_x_zz_yyz, to_z_x_zz_yzz, to_z_x_zz_zzz, to_z_xx, to_z_xxxx, to_z_xxxy, to_z_xxxz, to_z_xxyy, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yz, to_z_zz, to_zzz_xx, to_zzz_xxxx, to_zzz_xxxy, to_zzz_xxxz, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zz_xxx[k] = 6.0 * to_z_xx[k] - 4.0 * to_z_xxxx[k] * tke_0 - 6.0 * to_zzz_xx[k] * tbe_0 + 4.0 * to_zzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_zz_xxy[k] = 4.0 * to_z_xy[k] - 4.0 * to_z_xxxy[k] * tke_0 - 4.0 * to_zzz_xy[k] * tbe_0 + 4.0 * to_zzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_zz_xxz[k] = 4.0 * to_z_xz[k] - 4.0 * to_z_xxxz[k] * tke_0 - 4.0 * to_zzz_xz[k] * tbe_0 + 4.0 * to_zzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_zz_xyy[k] = 2.0 * to_z_yy[k] - 4.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_zzz_yy[k] * tbe_0 + 4.0 * to_zzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_zz_xyz[k] = 2.0 * to_z_yz[k] - 4.0 * to_z_xxyz[k] * tke_0 - 2.0 * to_zzz_yz[k] * tbe_0 + 4.0 * to_zzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_zz_xzz[k] = 2.0 * to_z_zz[k] - 4.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_zzz_zz[k] * tbe_0 + 4.0 * to_zzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_zz_yyy[k] = -4.0 * to_z_xyyy[k] * tke_0 + 4.0 * to_zzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_zz_yyz[k] = -4.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_zzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_zz_yzz[k] = -4.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_zzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_zz_zzz[k] = -4.0 * to_z_xzzz[k] * tke_0 + 4.0 * to_zzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-430 components of targeted buffer : DF

        auto to_z_y_xx_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 0);

        auto to_z_y_xx_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 1);

        auto to_z_y_xx_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 2);

        auto to_z_y_xx_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 3);

        auto to_z_y_xx_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 4);

        auto to_z_y_xx_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 5);

        auto to_z_y_xx_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 6);

        auto to_z_y_xx_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 7);

        auto to_z_y_xx_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 8);

        auto to_z_y_xx_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_z_y_xx_xxx, to_z_y_xx_xxy, to_z_y_xx_xxz, to_z_y_xx_xyy, to_z_y_xx_xyz, to_z_y_xx_xzz, to_z_y_xx_yyy, to_z_y_xx_yyz, to_z_y_xx_yzz, to_z_y_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xx_xxx[k] = 4.0 * to_xxz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xx_xxy[k] = -2.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xxz[k] = 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xx_xyy[k] = -4.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xx_xyz[k] = -2.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xx_xzz[k] = 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xx_yyy[k] = -6.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xx_yyz[k] = -4.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xx_yzz[k] = -2.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xx_zzz[k] = 4.0 * to_xxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 430-440 components of targeted buffer : DF

        auto to_z_y_xy_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 10);

        auto to_z_y_xy_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 11);

        auto to_z_y_xy_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 12);

        auto to_z_y_xy_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 13);

        auto to_z_y_xy_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 14);

        auto to_z_y_xy_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 15);

        auto to_z_y_xy_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 16);

        auto to_z_y_xy_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 17);

        auto to_z_y_xy_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 18);

        auto to_z_y_xy_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_z_y_xy_xxx, to_z_y_xy_xxy, to_z_y_xy_xxz, to_z_y_xy_xyy, to_z_y_xy_xyz, to_z_y_xy_xzz, to_z_y_xy_yyy, to_z_y_xy_yyz, to_z_y_xy_yzz, to_z_y_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xy_xxx[k] = 4.0 * to_xyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xy_xxy[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xxz[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xy_xyy[k] = -4.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xy_xyz[k] = -2.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xy_xzz[k] = 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xy_yyy[k] = -6.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xy_yyz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xy_yzz[k] = -2.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xy_zzz[k] = 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 440-450 components of targeted buffer : DF

        auto to_z_y_xz_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 20);

        auto to_z_y_xz_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 21);

        auto to_z_y_xz_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 22);

        auto to_z_y_xz_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 23);

        auto to_z_y_xz_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 24);

        auto to_z_y_xz_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 25);

        auto to_z_y_xz_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 26);

        auto to_z_y_xz_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 27);

        auto to_z_y_xz_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 28);

        auto to_z_y_xz_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_xx, to_x_xxxy, to_x_xxyy, to_x_xxyz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_yy, to_x_yyyy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_z_y_xz_xxx, to_z_y_xz_xxy, to_z_y_xz_xxz, to_z_y_xz_xyy, to_z_y_xz_xyz, to_z_y_xz_xzz, to_z_y_xz_yyy, to_z_y_xz_yyz, to_z_y_xz_yzz, to_z_y_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xz_xxx[k] = -2.0 * to_x_xxxy[k] * tke_0 + 4.0 * to_xzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xz_xxy[k] = to_x_xx[k] - 2.0 * to_x_xxyy[k] * tke_0 - 2.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xxz[k] = -2.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xz_xyy[k] = 2.0 * to_x_xy[k] - 2.0 * to_x_xyyy[k] * tke_0 - 4.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xz_xyz[k] = to_x_xz[k] - 2.0 * to_x_xyyz[k] * tke_0 - 2.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xz_xzz[k] = -2.0 * to_x_xyzz[k] * tke_0 + 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xz_yyy[k] = 3.0 * to_x_yy[k] - 2.0 * to_x_yyyy[k] * tke_0 - 6.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xz_yyz[k] = 2.0 * to_x_yz[k] - 2.0 * to_x_yyyz[k] * tke_0 - 4.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xz_yzz[k] = to_x_zz[k] - 2.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xz_zzz[k] = -2.0 * to_x_yzzz[k] * tke_0 + 4.0 * to_xzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-460 components of targeted buffer : DF

        auto to_z_y_yy_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 30);

        auto to_z_y_yy_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 31);

        auto to_z_y_yy_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 32);

        auto to_z_y_yy_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 33);

        auto to_z_y_yy_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 34);

        auto to_z_y_yy_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 35);

        auto to_z_y_yy_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 36);

        auto to_z_y_yy_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 37);

        auto to_z_y_yy_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 38);

        auto to_z_y_yy_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_z_y_yy_xxx, to_z_y_yy_xxy, to_z_y_yy_xxz, to_z_y_yy_xyy, to_z_y_yy_xyz, to_z_y_yy_xzz, to_z_y_yy_yyy, to_z_y_yy_yyz, to_z_y_yy_yzz, to_z_y_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yy_xxx[k] = 4.0 * to_yyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yy_xxy[k] = -2.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xxz[k] = 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yy_xyy[k] = -4.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yy_xyz[k] = -2.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yy_xzz[k] = 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yy_yyy[k] = -6.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yy_yyz[k] = -4.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yy_yzz[k] = -2.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yy_zzz[k] = 4.0 * to_yyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 460-470 components of targeted buffer : DF

        auto to_z_y_yz_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 40);

        auto to_z_y_yz_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 41);

        auto to_z_y_yz_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 42);

        auto to_z_y_yz_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 43);

        auto to_z_y_yz_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 44);

        auto to_z_y_yz_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 45);

        auto to_z_y_yz_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 46);

        auto to_z_y_yz_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 47);

        auto to_z_y_yz_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 48);

        auto to_z_y_yz_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_xx, to_y_xxxy, to_y_xxyy, to_y_xxyz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_yy, to_y_yyyy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_z_y_yz_xxx, to_z_y_yz_xxy, to_z_y_yz_xxz, to_z_y_yz_xyy, to_z_y_yz_xyz, to_z_y_yz_xzz, to_z_y_yz_yyy, to_z_y_yz_yyz, to_z_y_yz_yzz, to_z_y_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yz_xxx[k] = -2.0 * to_y_xxxy[k] * tke_0 + 4.0 * to_yzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yz_xxy[k] = to_y_xx[k] - 2.0 * to_y_xxyy[k] * tke_0 - 2.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xxz[k] = -2.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yz_xyy[k] = 2.0 * to_y_xy[k] - 2.0 * to_y_xyyy[k] * tke_0 - 4.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yz_xyz[k] = to_y_xz[k] - 2.0 * to_y_xyyz[k] * tke_0 - 2.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yz_xzz[k] = -2.0 * to_y_xyzz[k] * tke_0 + 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yz_yyy[k] = 3.0 * to_y_yy[k] - 2.0 * to_y_yyyy[k] * tke_0 - 6.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yz_yyz[k] = 2.0 * to_y_yz[k] - 2.0 * to_y_yyyz[k] * tke_0 - 4.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yz_yzz[k] = to_y_zz[k] - 2.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yz_zzz[k] = -2.0 * to_y_yzzz[k] * tke_0 + 4.0 * to_yzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 470-480 components of targeted buffer : DF

        auto to_z_y_zz_xxx = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 50);

        auto to_z_y_zz_xxy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 51);

        auto to_z_y_zz_xxz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 52);

        auto to_z_y_zz_xyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 53);

        auto to_z_y_zz_xyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 54);

        auto to_z_y_zz_xzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 55);

        auto to_z_y_zz_yyy = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 56);

        auto to_z_y_zz_yyz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 57);

        auto to_z_y_zz_yzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 58);

        auto to_z_y_zz_zzz = pbuffer.data(idx_op_geom_101_df + 7 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_xx, to_z_xxxy, to_z_xxyy, to_z_xxyz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_y_zz_xxx, to_z_y_zz_xxy, to_z_y_zz_xxz, to_z_y_zz_xyy, to_z_y_zz_xyz, to_z_y_zz_xzz, to_z_y_zz_yyy, to_z_y_zz_yyz, to_z_y_zz_yzz, to_z_y_zz_zzz, to_z_yy, to_z_yyyy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, to_zzz_xx, to_zzz_xxxy, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_yy, to_zzz_yyyy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zz_xxx[k] = -4.0 * to_z_xxxy[k] * tke_0 + 4.0 * to_zzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_zz_xxy[k] = 2.0 * to_z_xx[k] - 4.0 * to_z_xxyy[k] * tke_0 - 2.0 * to_zzz_xx[k] * tbe_0 + 4.0 * to_zzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xxz[k] = -4.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_zzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_zz_xyy[k] = 4.0 * to_z_xy[k] - 4.0 * to_z_xyyy[k] * tke_0 - 4.0 * to_zzz_xy[k] * tbe_0 + 4.0 * to_zzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_zz_xyz[k] = 2.0 * to_z_xz[k] - 4.0 * to_z_xyyz[k] * tke_0 - 2.0 * to_zzz_xz[k] * tbe_0 + 4.0 * to_zzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_zz_xzz[k] = -4.0 * to_z_xyzz[k] * tke_0 + 4.0 * to_zzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_zz_yyy[k] = 6.0 * to_z_yy[k] - 4.0 * to_z_yyyy[k] * tke_0 - 6.0 * to_zzz_yy[k] * tbe_0 + 4.0 * to_zzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_zz_yyz[k] = 4.0 * to_z_yz[k] - 4.0 * to_z_yyyz[k] * tke_0 - 4.0 * to_zzz_yz[k] * tbe_0 + 4.0 * to_zzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_zz_yzz[k] = 2.0 * to_z_zz[k] - 4.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_zzz_zz[k] * tbe_0 + 4.0 * to_zzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_zz_zzz[k] = -4.0 * to_z_yzzz[k] * tke_0 + 4.0 * to_zzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-490 components of targeted buffer : DF

        auto to_z_z_xx_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 0);

        auto to_z_z_xx_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 1);

        auto to_z_z_xx_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 2);

        auto to_z_z_xx_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 3);

        auto to_z_z_xx_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 4);

        auto to_z_z_xx_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 5);

        auto to_z_z_xx_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 6);

        auto to_z_z_xx_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 7);

        auto to_z_z_xx_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 8);

        auto to_z_z_xx_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 9);

        #pragma omp simd aligned(to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, to_z_z_xx_xxx, to_z_z_xx_xxy, to_z_z_xx_xxz, to_z_z_xx_xyy, to_z_z_xx_xyz, to_z_z_xx_xzz, to_z_z_xx_yyy, to_z_z_xx_yyz, to_z_z_xx_yzz, to_z_z_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xx_xxx[k] = 4.0 * to_xxz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxy[k] = 4.0 * to_xxz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xxz[k] = -2.0 * to_xxz_xx[k] * tbe_0 + 4.0 * to_xxz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xyy[k] = 4.0 * to_xxz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xx_xyz[k] = -2.0 * to_xxz_xy[k] * tbe_0 + 4.0 * to_xxz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xx_xzz[k] = -4.0 * to_xxz_xz[k] * tbe_0 + 4.0 * to_xxz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yyy[k] = 4.0 * to_xxz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xx_yyz[k] = -2.0 * to_xxz_yy[k] * tbe_0 + 4.0 * to_xxz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xx_yzz[k] = -4.0 * to_xxz_yz[k] * tbe_0 + 4.0 * to_xxz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xx_zzz[k] = -6.0 * to_xxz_zz[k] * tbe_0 + 4.0 * to_xxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 490-500 components of targeted buffer : DF

        auto to_z_z_xy_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 10);

        auto to_z_z_xy_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 11);

        auto to_z_z_xy_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 12);

        auto to_z_z_xy_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 13);

        auto to_z_z_xy_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 14);

        auto to_z_z_xy_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 15);

        auto to_z_z_xy_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 16);

        auto to_z_z_xy_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 17);

        auto to_z_z_xy_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 18);

        auto to_z_z_xy_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 19);

        #pragma omp simd aligned(to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, to_z_z_xy_xxx, to_z_z_xy_xxy, to_z_z_xy_xxz, to_z_z_xy_xyy, to_z_z_xy_xyz, to_z_z_xy_xzz, to_z_z_xy_yyy, to_z_z_xy_yyz, to_z_z_xy_yzz, to_z_z_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xy_xxx[k] = 4.0 * to_xyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxy[k] = 4.0 * to_xyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xxz[k] = -2.0 * to_xyz_xx[k] * tbe_0 + 4.0 * to_xyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xyy[k] = 4.0 * to_xyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xy_xyz[k] = -2.0 * to_xyz_xy[k] * tbe_0 + 4.0 * to_xyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xy_xzz[k] = -4.0 * to_xyz_xz[k] * tbe_0 + 4.0 * to_xyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yyy[k] = 4.0 * to_xyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xy_yyz[k] = -2.0 * to_xyz_yy[k] * tbe_0 + 4.0 * to_xyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xy_yzz[k] = -4.0 * to_xyz_yz[k] * tbe_0 + 4.0 * to_xyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xy_zzz[k] = -6.0 * to_xyz_zz[k] * tbe_0 + 4.0 * to_xyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 500-510 components of targeted buffer : DF

        auto to_z_z_xz_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 20);

        auto to_z_z_xz_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 21);

        auto to_z_z_xz_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 22);

        auto to_z_z_xz_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 23);

        auto to_z_z_xz_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 24);

        auto to_z_z_xz_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 25);

        auto to_z_z_xz_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 26);

        auto to_z_z_xz_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 27);

        auto to_z_z_xz_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 28);

        auto to_z_z_xz_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_xx, to_x_xxxz, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_x_zzzz, to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, to_z_z_xz_xxx, to_z_z_xz_xxy, to_z_z_xz_xxz, to_z_z_xz_xyy, to_z_z_xz_xyz, to_z_z_xz_xzz, to_z_z_xz_yyy, to_z_z_xz_yyz, to_z_z_xz_yzz, to_z_z_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xz_xxx[k] = -2.0 * to_x_xxxz[k] * tke_0 + 4.0 * to_xzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxy[k] = -2.0 * to_x_xxyz[k] * tke_0 + 4.0 * to_xzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xxz[k] = to_x_xx[k] - 2.0 * to_x_xxzz[k] * tke_0 - 2.0 * to_xzz_xx[k] * tbe_0 + 4.0 * to_xzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xyy[k] = -2.0 * to_x_xyyz[k] * tke_0 + 4.0 * to_xzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xz_xyz[k] = to_x_xy[k] - 2.0 * to_x_xyzz[k] * tke_0 - 2.0 * to_xzz_xy[k] * tbe_0 + 4.0 * to_xzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xz_xzz[k] = 2.0 * to_x_xz[k] - 2.0 * to_x_xzzz[k] * tke_0 - 4.0 * to_xzz_xz[k] * tbe_0 + 4.0 * to_xzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yyy[k] = -2.0 * to_x_yyyz[k] * tke_0 + 4.0 * to_xzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xz_yyz[k] = to_x_yy[k] - 2.0 * to_x_yyzz[k] * tke_0 - 2.0 * to_xzz_yy[k] * tbe_0 + 4.0 * to_xzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xz_yzz[k] = 2.0 * to_x_yz[k] - 2.0 * to_x_yzzz[k] * tke_0 - 4.0 * to_xzz_yz[k] * tbe_0 + 4.0 * to_xzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xz_zzz[k] = 3.0 * to_x_zz[k] - 2.0 * to_x_zzzz[k] * tke_0 - 6.0 * to_xzz_zz[k] * tbe_0 + 4.0 * to_xzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-520 components of targeted buffer : DF

        auto to_z_z_yy_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 30);

        auto to_z_z_yy_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 31);

        auto to_z_z_yy_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 32);

        auto to_z_z_yy_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 33);

        auto to_z_z_yy_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 34);

        auto to_z_z_yy_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 35);

        auto to_z_z_yy_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 36);

        auto to_z_z_yy_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 37);

        auto to_z_z_yy_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 38);

        auto to_z_z_yy_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 39);

        #pragma omp simd aligned(to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, to_z_z_yy_xxx, to_z_z_yy_xxy, to_z_z_yy_xxz, to_z_z_yy_xyy, to_z_z_yy_xyz, to_z_z_yy_xzz, to_z_z_yy_yyy, to_z_z_yy_yyz, to_z_z_yy_yzz, to_z_z_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yy_xxx[k] = 4.0 * to_yyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxy[k] = 4.0 * to_yyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xxz[k] = -2.0 * to_yyz_xx[k] * tbe_0 + 4.0 * to_yyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xyy[k] = 4.0 * to_yyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yy_xyz[k] = -2.0 * to_yyz_xy[k] * tbe_0 + 4.0 * to_yyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yy_xzz[k] = -4.0 * to_yyz_xz[k] * tbe_0 + 4.0 * to_yyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yyy[k] = 4.0 * to_yyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yy_yyz[k] = -2.0 * to_yyz_yy[k] * tbe_0 + 4.0 * to_yyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yy_yzz[k] = -4.0 * to_yyz_yz[k] * tbe_0 + 4.0 * to_yyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yy_zzz[k] = -6.0 * to_yyz_zz[k] * tbe_0 + 4.0 * to_yyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 520-530 components of targeted buffer : DF

        auto to_z_z_yz_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 40);

        auto to_z_z_yz_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 41);

        auto to_z_z_yz_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 42);

        auto to_z_z_yz_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 43);

        auto to_z_z_yz_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 44);

        auto to_z_z_yz_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 45);

        auto to_z_z_yz_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 46);

        auto to_z_z_yz_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 47);

        auto to_z_z_yz_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 48);

        auto to_z_z_yz_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 49);

        #pragma omp simd aligned(to_y_xx, to_y_xxxz, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, to_y_zzzz, to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, to_z_z_yz_xxx, to_z_z_yz_xxy, to_z_z_yz_xxz, to_z_z_yz_xyy, to_z_z_yz_xyz, to_z_z_yz_xzz, to_z_z_yz_yyy, to_z_z_yz_yyz, to_z_z_yz_yzz, to_z_z_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yz_xxx[k] = -2.0 * to_y_xxxz[k] * tke_0 + 4.0 * to_yzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxy[k] = -2.0 * to_y_xxyz[k] * tke_0 + 4.0 * to_yzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xxz[k] = to_y_xx[k] - 2.0 * to_y_xxzz[k] * tke_0 - 2.0 * to_yzz_xx[k] * tbe_0 + 4.0 * to_yzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xyy[k] = -2.0 * to_y_xyyz[k] * tke_0 + 4.0 * to_yzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yz_xyz[k] = to_y_xy[k] - 2.0 * to_y_xyzz[k] * tke_0 - 2.0 * to_yzz_xy[k] * tbe_0 + 4.0 * to_yzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yz_xzz[k] = 2.0 * to_y_xz[k] - 2.0 * to_y_xzzz[k] * tke_0 - 4.0 * to_yzz_xz[k] * tbe_0 + 4.0 * to_yzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yyy[k] = -2.0 * to_y_yyyz[k] * tke_0 + 4.0 * to_yzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yz_yyz[k] = to_y_yy[k] - 2.0 * to_y_yyzz[k] * tke_0 - 2.0 * to_yzz_yy[k] * tbe_0 + 4.0 * to_yzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yz_yzz[k] = 2.0 * to_y_yz[k] - 2.0 * to_y_yzzz[k] * tke_0 - 4.0 * to_yzz_yz[k] * tbe_0 + 4.0 * to_yzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yz_zzz[k] = 3.0 * to_y_zz[k] - 2.0 * to_y_zzzz[k] * tke_0 - 6.0 * to_yzz_zz[k] * tbe_0 + 4.0 * to_yzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 530-540 components of targeted buffer : DF

        auto to_z_z_zz_xxx = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 50);

        auto to_z_z_zz_xxy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 51);

        auto to_z_z_zz_xxz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 52);

        auto to_z_z_zz_xyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 53);

        auto to_z_z_zz_xyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 54);

        auto to_z_z_zz_xzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 55);

        auto to_z_z_zz_yyy = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 56);

        auto to_z_z_zz_yyz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 57);

        auto to_z_z_zz_yzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 58);

        auto to_z_z_zz_zzz = pbuffer.data(idx_op_geom_101_df + 8 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_xx, to_z_xxxz, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_z_zz_xxx, to_z_z_zz_xxy, to_z_z_zz_xxz, to_z_z_zz_xyy, to_z_z_zz_xyz, to_z_z_zz_xzz, to_z_z_zz_yyy, to_z_z_zz_yyz, to_z_z_zz_yzz, to_z_z_zz_zzz, to_z_zz, to_z_zzzz, to_zzz_xx, to_zzz_xxxz, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zz_xxx[k] = -4.0 * to_z_xxxz[k] * tke_0 + 4.0 * to_zzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxy[k] = -4.0 * to_z_xxyz[k] * tke_0 + 4.0 * to_zzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xxz[k] = 2.0 * to_z_xx[k] - 4.0 * to_z_xxzz[k] * tke_0 - 2.0 * to_zzz_xx[k] * tbe_0 + 4.0 * to_zzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xyy[k] = -4.0 * to_z_xyyz[k] * tke_0 + 4.0 * to_zzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_zz_xyz[k] = 2.0 * to_z_xy[k] - 4.0 * to_z_xyzz[k] * tke_0 - 2.0 * to_zzz_xy[k] * tbe_0 + 4.0 * to_zzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_zz_xzz[k] = 4.0 * to_z_xz[k] - 4.0 * to_z_xzzz[k] * tke_0 - 4.0 * to_zzz_xz[k] * tbe_0 + 4.0 * to_zzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yyy[k] = -4.0 * to_z_yyyz[k] * tke_0 + 4.0 * to_zzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_zz_yyz[k] = 2.0 * to_z_yy[k] - 4.0 * to_z_yyzz[k] * tke_0 - 2.0 * to_zzz_yy[k] * tbe_0 + 4.0 * to_zzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_zz_yzz[k] = 4.0 * to_z_yz[k] - 4.0 * to_z_yzzz[k] * tke_0 - 4.0 * to_zzz_yz[k] * tbe_0 + 4.0 * to_zzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_zz_zzz[k] = 6.0 * to_z_zz[k] - 4.0 * to_z_zzzz[k] * tke_0 - 6.0 * to_zzz_zz[k] * tbe_0 + 4.0 * to_zzz_zzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

