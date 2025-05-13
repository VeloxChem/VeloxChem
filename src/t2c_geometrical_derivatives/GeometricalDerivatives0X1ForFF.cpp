#include "GeometricalDerivatives0X1ForFF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_ff(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_ff,
                       const int idx_op_fd,
                       const int idx_op_fg,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
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

        // Set up 0-10 components of targeted buffer : FF

        auto to_0_x_xxx_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 0);

        auto to_0_x_xxx_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 1);

        auto to_0_x_xxx_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 2);

        auto to_0_x_xxx_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 3);

        auto to_0_x_xxx_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 4);

        auto to_0_x_xxx_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 5);

        auto to_0_x_xxx_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 6);

        auto to_0_x_xxx_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 7);

        auto to_0_x_xxx_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 8);

        auto to_0_x_xxx_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_0_x_xxx_xxx, to_0_x_xxx_xxy, to_0_x_xxx_xxz, to_0_x_xxx_xyy, to_0_x_xxx_xyz, to_0_x_xxx_xzz, to_0_x_xxx_yyy, to_0_x_xxx_yyz, to_0_x_xxx_yzz, to_0_x_xxx_zzz, to_xxx_xx, to_xxx_xxxx, to_xxx_xxxy, to_xxx_xxxz, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxx_xxx[k] = -3.0 * to_xxx_xx[k] + 2.0 * to_xxx_xxxx[k] * tke_0;

            to_0_x_xxx_xxy[k] = -2.0 * to_xxx_xy[k] + 2.0 * to_xxx_xxxy[k] * tke_0;

            to_0_x_xxx_xxz[k] = -2.0 * to_xxx_xz[k] + 2.0 * to_xxx_xxxz[k] * tke_0;

            to_0_x_xxx_xyy[k] = -to_xxx_yy[k] + 2.0 * to_xxx_xxyy[k] * tke_0;

            to_0_x_xxx_xyz[k] = -to_xxx_yz[k] + 2.0 * to_xxx_xxyz[k] * tke_0;

            to_0_x_xxx_xzz[k] = -to_xxx_zz[k] + 2.0 * to_xxx_xxzz[k] * tke_0;

            to_0_x_xxx_yyy[k] = 2.0 * to_xxx_xyyy[k] * tke_0;

            to_0_x_xxx_yyz[k] = 2.0 * to_xxx_xyyz[k] * tke_0;

            to_0_x_xxx_yzz[k] = 2.0 * to_xxx_xyzz[k] * tke_0;

            to_0_x_xxx_zzz[k] = 2.0 * to_xxx_xzzz[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : FF

        auto to_0_x_xxy_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 10);

        auto to_0_x_xxy_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 11);

        auto to_0_x_xxy_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 12);

        auto to_0_x_xxy_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 13);

        auto to_0_x_xxy_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 14);

        auto to_0_x_xxy_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 15);

        auto to_0_x_xxy_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 16);

        auto to_0_x_xxy_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 17);

        auto to_0_x_xxy_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 18);

        auto to_0_x_xxy_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_0_x_xxy_xxx, to_0_x_xxy_xxy, to_0_x_xxy_xxz, to_0_x_xxy_xyy, to_0_x_xxy_xyz, to_0_x_xxy_xzz, to_0_x_xxy_yyy, to_0_x_xxy_yyz, to_0_x_xxy_yzz, to_0_x_xxy_zzz, to_xxy_xx, to_xxy_xxxx, to_xxy_xxxy, to_xxy_xxxz, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yz, to_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxy_xxx[k] = -3.0 * to_xxy_xx[k] + 2.0 * to_xxy_xxxx[k] * tke_0;

            to_0_x_xxy_xxy[k] = -2.0 * to_xxy_xy[k] + 2.0 * to_xxy_xxxy[k] * tke_0;

            to_0_x_xxy_xxz[k] = -2.0 * to_xxy_xz[k] + 2.0 * to_xxy_xxxz[k] * tke_0;

            to_0_x_xxy_xyy[k] = -to_xxy_yy[k] + 2.0 * to_xxy_xxyy[k] * tke_0;

            to_0_x_xxy_xyz[k] = -to_xxy_yz[k] + 2.0 * to_xxy_xxyz[k] * tke_0;

            to_0_x_xxy_xzz[k] = -to_xxy_zz[k] + 2.0 * to_xxy_xxzz[k] * tke_0;

            to_0_x_xxy_yyy[k] = 2.0 * to_xxy_xyyy[k] * tke_0;

            to_0_x_xxy_yyz[k] = 2.0 * to_xxy_xyyz[k] * tke_0;

            to_0_x_xxy_yzz[k] = 2.0 * to_xxy_xyzz[k] * tke_0;

            to_0_x_xxy_zzz[k] = 2.0 * to_xxy_xzzz[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : FF

        auto to_0_x_xxz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 20);

        auto to_0_x_xxz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 21);

        auto to_0_x_xxz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 22);

        auto to_0_x_xxz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 23);

        auto to_0_x_xxz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 24);

        auto to_0_x_xxz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 25);

        auto to_0_x_xxz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 26);

        auto to_0_x_xxz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 27);

        auto to_0_x_xxz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 28);

        auto to_0_x_xxz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_0_x_xxz_xxx, to_0_x_xxz_xxy, to_0_x_xxz_xxz, to_0_x_xxz_xyy, to_0_x_xxz_xyz, to_0_x_xxz_xzz, to_0_x_xxz_yyy, to_0_x_xxz_yyz, to_0_x_xxz_yzz, to_0_x_xxz_zzz, to_xxz_xx, to_xxz_xxxx, to_xxz_xxxy, to_xxz_xxxz, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yz, to_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxz_xxx[k] = -3.0 * to_xxz_xx[k] + 2.0 * to_xxz_xxxx[k] * tke_0;

            to_0_x_xxz_xxy[k] = -2.0 * to_xxz_xy[k] + 2.0 * to_xxz_xxxy[k] * tke_0;

            to_0_x_xxz_xxz[k] = -2.0 * to_xxz_xz[k] + 2.0 * to_xxz_xxxz[k] * tke_0;

            to_0_x_xxz_xyy[k] = -to_xxz_yy[k] + 2.0 * to_xxz_xxyy[k] * tke_0;

            to_0_x_xxz_xyz[k] = -to_xxz_yz[k] + 2.0 * to_xxz_xxyz[k] * tke_0;

            to_0_x_xxz_xzz[k] = -to_xxz_zz[k] + 2.0 * to_xxz_xxzz[k] * tke_0;

            to_0_x_xxz_yyy[k] = 2.0 * to_xxz_xyyy[k] * tke_0;

            to_0_x_xxz_yyz[k] = 2.0 * to_xxz_xyyz[k] * tke_0;

            to_0_x_xxz_yzz[k] = 2.0 * to_xxz_xyzz[k] * tke_0;

            to_0_x_xxz_zzz[k] = 2.0 * to_xxz_xzzz[k] * tke_0;
        }

        // Set up 30-40 components of targeted buffer : FF

        auto to_0_x_xyy_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 30);

        auto to_0_x_xyy_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 31);

        auto to_0_x_xyy_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 32);

        auto to_0_x_xyy_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 33);

        auto to_0_x_xyy_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 34);

        auto to_0_x_xyy_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 35);

        auto to_0_x_xyy_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 36);

        auto to_0_x_xyy_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 37);

        auto to_0_x_xyy_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 38);

        auto to_0_x_xyy_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_0_x_xyy_xxx, to_0_x_xyy_xxy, to_0_x_xyy_xxz, to_0_x_xyy_xyy, to_0_x_xyy_xyz, to_0_x_xyy_xzz, to_0_x_xyy_yyy, to_0_x_xyy_yyz, to_0_x_xyy_yzz, to_0_x_xyy_zzz, to_xyy_xx, to_xyy_xxxx, to_xyy_xxxy, to_xyy_xxxz, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyy_xxx[k] = -3.0 * to_xyy_xx[k] + 2.0 * to_xyy_xxxx[k] * tke_0;

            to_0_x_xyy_xxy[k] = -2.0 * to_xyy_xy[k] + 2.0 * to_xyy_xxxy[k] * tke_0;

            to_0_x_xyy_xxz[k] = -2.0 * to_xyy_xz[k] + 2.0 * to_xyy_xxxz[k] * tke_0;

            to_0_x_xyy_xyy[k] = -to_xyy_yy[k] + 2.0 * to_xyy_xxyy[k] * tke_0;

            to_0_x_xyy_xyz[k] = -to_xyy_yz[k] + 2.0 * to_xyy_xxyz[k] * tke_0;

            to_0_x_xyy_xzz[k] = -to_xyy_zz[k] + 2.0 * to_xyy_xxzz[k] * tke_0;

            to_0_x_xyy_yyy[k] = 2.0 * to_xyy_xyyy[k] * tke_0;

            to_0_x_xyy_yyz[k] = 2.0 * to_xyy_xyyz[k] * tke_0;

            to_0_x_xyy_yzz[k] = 2.0 * to_xyy_xyzz[k] * tke_0;

            to_0_x_xyy_zzz[k] = 2.0 * to_xyy_xzzz[k] * tke_0;
        }

        // Set up 40-50 components of targeted buffer : FF

        auto to_0_x_xyz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 40);

        auto to_0_x_xyz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 41);

        auto to_0_x_xyz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 42);

        auto to_0_x_xyz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 43);

        auto to_0_x_xyz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 44);

        auto to_0_x_xyz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 45);

        auto to_0_x_xyz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 46);

        auto to_0_x_xyz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 47);

        auto to_0_x_xyz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 48);

        auto to_0_x_xyz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_0_x_xyz_xxx, to_0_x_xyz_xxy, to_0_x_xyz_xxz, to_0_x_xyz_xyy, to_0_x_xyz_xyz, to_0_x_xyz_xzz, to_0_x_xyz_yyy, to_0_x_xyz_yyz, to_0_x_xyz_yzz, to_0_x_xyz_zzz, to_xyz_xx, to_xyz_xxxx, to_xyz_xxxy, to_xyz_xxxz, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyz_xxx[k] = -3.0 * to_xyz_xx[k] + 2.0 * to_xyz_xxxx[k] * tke_0;

            to_0_x_xyz_xxy[k] = -2.0 * to_xyz_xy[k] + 2.0 * to_xyz_xxxy[k] * tke_0;

            to_0_x_xyz_xxz[k] = -2.0 * to_xyz_xz[k] + 2.0 * to_xyz_xxxz[k] * tke_0;

            to_0_x_xyz_xyy[k] = -to_xyz_yy[k] + 2.0 * to_xyz_xxyy[k] * tke_0;

            to_0_x_xyz_xyz[k] = -to_xyz_yz[k] + 2.0 * to_xyz_xxyz[k] * tke_0;

            to_0_x_xyz_xzz[k] = -to_xyz_zz[k] + 2.0 * to_xyz_xxzz[k] * tke_0;

            to_0_x_xyz_yyy[k] = 2.0 * to_xyz_xyyy[k] * tke_0;

            to_0_x_xyz_yyz[k] = 2.0 * to_xyz_xyyz[k] * tke_0;

            to_0_x_xyz_yzz[k] = 2.0 * to_xyz_xyzz[k] * tke_0;

            to_0_x_xyz_zzz[k] = 2.0 * to_xyz_xzzz[k] * tke_0;
        }

        // Set up 50-60 components of targeted buffer : FF

        auto to_0_x_xzz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 50);

        auto to_0_x_xzz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 51);

        auto to_0_x_xzz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 52);

        auto to_0_x_xzz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 53);

        auto to_0_x_xzz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 54);

        auto to_0_x_xzz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 55);

        auto to_0_x_xzz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 56);

        auto to_0_x_xzz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 57);

        auto to_0_x_xzz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 58);

        auto to_0_x_xzz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_0_x_xzz_xxx, to_0_x_xzz_xxy, to_0_x_xzz_xxz, to_0_x_xzz_xyy, to_0_x_xzz_xyz, to_0_x_xzz_xzz, to_0_x_xzz_yyy, to_0_x_xzz_yyz, to_0_x_xzz_yzz, to_0_x_xzz_zzz, to_xzz_xx, to_xzz_xxxx, to_xzz_xxxy, to_xzz_xxxz, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzz_xxx[k] = -3.0 * to_xzz_xx[k] + 2.0 * to_xzz_xxxx[k] * tke_0;

            to_0_x_xzz_xxy[k] = -2.0 * to_xzz_xy[k] + 2.0 * to_xzz_xxxy[k] * tke_0;

            to_0_x_xzz_xxz[k] = -2.0 * to_xzz_xz[k] + 2.0 * to_xzz_xxxz[k] * tke_0;

            to_0_x_xzz_xyy[k] = -to_xzz_yy[k] + 2.0 * to_xzz_xxyy[k] * tke_0;

            to_0_x_xzz_xyz[k] = -to_xzz_yz[k] + 2.0 * to_xzz_xxyz[k] * tke_0;

            to_0_x_xzz_xzz[k] = -to_xzz_zz[k] + 2.0 * to_xzz_xxzz[k] * tke_0;

            to_0_x_xzz_yyy[k] = 2.0 * to_xzz_xyyy[k] * tke_0;

            to_0_x_xzz_yyz[k] = 2.0 * to_xzz_xyyz[k] * tke_0;

            to_0_x_xzz_yzz[k] = 2.0 * to_xzz_xyzz[k] * tke_0;

            to_0_x_xzz_zzz[k] = 2.0 * to_xzz_xzzz[k] * tke_0;
        }

        // Set up 60-70 components of targeted buffer : FF

        auto to_0_x_yyy_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 60);

        auto to_0_x_yyy_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 61);

        auto to_0_x_yyy_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 62);

        auto to_0_x_yyy_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 63);

        auto to_0_x_yyy_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 64);

        auto to_0_x_yyy_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 65);

        auto to_0_x_yyy_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 66);

        auto to_0_x_yyy_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 67);

        auto to_0_x_yyy_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 68);

        auto to_0_x_yyy_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_0_x_yyy_xxx, to_0_x_yyy_xxy, to_0_x_yyy_xxz, to_0_x_yyy_xyy, to_0_x_yyy_xyz, to_0_x_yyy_xzz, to_0_x_yyy_yyy, to_0_x_yyy_yyz, to_0_x_yyy_yzz, to_0_x_yyy_zzz, to_yyy_xx, to_yyy_xxxx, to_yyy_xxxy, to_yyy_xxxz, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyy_xxx[k] = -3.0 * to_yyy_xx[k] + 2.0 * to_yyy_xxxx[k] * tke_0;

            to_0_x_yyy_xxy[k] = -2.0 * to_yyy_xy[k] + 2.0 * to_yyy_xxxy[k] * tke_0;

            to_0_x_yyy_xxz[k] = -2.0 * to_yyy_xz[k] + 2.0 * to_yyy_xxxz[k] * tke_0;

            to_0_x_yyy_xyy[k] = -to_yyy_yy[k] + 2.0 * to_yyy_xxyy[k] * tke_0;

            to_0_x_yyy_xyz[k] = -to_yyy_yz[k] + 2.0 * to_yyy_xxyz[k] * tke_0;

            to_0_x_yyy_xzz[k] = -to_yyy_zz[k] + 2.0 * to_yyy_xxzz[k] * tke_0;

            to_0_x_yyy_yyy[k] = 2.0 * to_yyy_xyyy[k] * tke_0;

            to_0_x_yyy_yyz[k] = 2.0 * to_yyy_xyyz[k] * tke_0;

            to_0_x_yyy_yzz[k] = 2.0 * to_yyy_xyzz[k] * tke_0;

            to_0_x_yyy_zzz[k] = 2.0 * to_yyy_xzzz[k] * tke_0;
        }

        // Set up 70-80 components of targeted buffer : FF

        auto to_0_x_yyz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 70);

        auto to_0_x_yyz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 71);

        auto to_0_x_yyz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 72);

        auto to_0_x_yyz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 73);

        auto to_0_x_yyz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 74);

        auto to_0_x_yyz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 75);

        auto to_0_x_yyz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 76);

        auto to_0_x_yyz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 77);

        auto to_0_x_yyz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 78);

        auto to_0_x_yyz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_0_x_yyz_xxx, to_0_x_yyz_xxy, to_0_x_yyz_xxz, to_0_x_yyz_xyy, to_0_x_yyz_xyz, to_0_x_yyz_xzz, to_0_x_yyz_yyy, to_0_x_yyz_yyz, to_0_x_yyz_yzz, to_0_x_yyz_zzz, to_yyz_xx, to_yyz_xxxx, to_yyz_xxxy, to_yyz_xxxz, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyz_xxx[k] = -3.0 * to_yyz_xx[k] + 2.0 * to_yyz_xxxx[k] * tke_0;

            to_0_x_yyz_xxy[k] = -2.0 * to_yyz_xy[k] + 2.0 * to_yyz_xxxy[k] * tke_0;

            to_0_x_yyz_xxz[k] = -2.0 * to_yyz_xz[k] + 2.0 * to_yyz_xxxz[k] * tke_0;

            to_0_x_yyz_xyy[k] = -to_yyz_yy[k] + 2.0 * to_yyz_xxyy[k] * tke_0;

            to_0_x_yyz_xyz[k] = -to_yyz_yz[k] + 2.0 * to_yyz_xxyz[k] * tke_0;

            to_0_x_yyz_xzz[k] = -to_yyz_zz[k] + 2.0 * to_yyz_xxzz[k] * tke_0;

            to_0_x_yyz_yyy[k] = 2.0 * to_yyz_xyyy[k] * tke_0;

            to_0_x_yyz_yyz[k] = 2.0 * to_yyz_xyyz[k] * tke_0;

            to_0_x_yyz_yzz[k] = 2.0 * to_yyz_xyzz[k] * tke_0;

            to_0_x_yyz_zzz[k] = 2.0 * to_yyz_xzzz[k] * tke_0;
        }

        // Set up 80-90 components of targeted buffer : FF

        auto to_0_x_yzz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 80);

        auto to_0_x_yzz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 81);

        auto to_0_x_yzz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 82);

        auto to_0_x_yzz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 83);

        auto to_0_x_yzz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 84);

        auto to_0_x_yzz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 85);

        auto to_0_x_yzz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 86);

        auto to_0_x_yzz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 87);

        auto to_0_x_yzz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 88);

        auto to_0_x_yzz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_0_x_yzz_xxx, to_0_x_yzz_xxy, to_0_x_yzz_xxz, to_0_x_yzz_xyy, to_0_x_yzz_xyz, to_0_x_yzz_xzz, to_0_x_yzz_yyy, to_0_x_yzz_yyz, to_0_x_yzz_yzz, to_0_x_yzz_zzz, to_yzz_xx, to_yzz_xxxx, to_yzz_xxxy, to_yzz_xxxz, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzz_xxx[k] = -3.0 * to_yzz_xx[k] + 2.0 * to_yzz_xxxx[k] * tke_0;

            to_0_x_yzz_xxy[k] = -2.0 * to_yzz_xy[k] + 2.0 * to_yzz_xxxy[k] * tke_0;

            to_0_x_yzz_xxz[k] = -2.0 * to_yzz_xz[k] + 2.0 * to_yzz_xxxz[k] * tke_0;

            to_0_x_yzz_xyy[k] = -to_yzz_yy[k] + 2.0 * to_yzz_xxyy[k] * tke_0;

            to_0_x_yzz_xyz[k] = -to_yzz_yz[k] + 2.0 * to_yzz_xxyz[k] * tke_0;

            to_0_x_yzz_xzz[k] = -to_yzz_zz[k] + 2.0 * to_yzz_xxzz[k] * tke_0;

            to_0_x_yzz_yyy[k] = 2.0 * to_yzz_xyyy[k] * tke_0;

            to_0_x_yzz_yyz[k] = 2.0 * to_yzz_xyyz[k] * tke_0;

            to_0_x_yzz_yzz[k] = 2.0 * to_yzz_xyzz[k] * tke_0;

            to_0_x_yzz_zzz[k] = 2.0 * to_yzz_xzzz[k] * tke_0;
        }

        // Set up 90-100 components of targeted buffer : FF

        auto to_0_x_zzz_xxx = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 90);

        auto to_0_x_zzz_xxy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 91);

        auto to_0_x_zzz_xxz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 92);

        auto to_0_x_zzz_xyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 93);

        auto to_0_x_zzz_xyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 94);

        auto to_0_x_zzz_xzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 95);

        auto to_0_x_zzz_yyy = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 96);

        auto to_0_x_zzz_yyz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 97);

        auto to_0_x_zzz_yzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 98);

        auto to_0_x_zzz_zzz = pbuffer.data(idx_op_geom_001_ff + 0 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_0_x_zzz_xxx, to_0_x_zzz_xxy, to_0_x_zzz_xxz, to_0_x_zzz_xyy, to_0_x_zzz_xyz, to_0_x_zzz_xzz, to_0_x_zzz_yyy, to_0_x_zzz_yyz, to_0_x_zzz_yzz, to_0_x_zzz_zzz, to_zzz_xx, to_zzz_xxxx, to_zzz_xxxy, to_zzz_xxxz, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzz_xxx[k] = -3.0 * to_zzz_xx[k] + 2.0 * to_zzz_xxxx[k] * tke_0;

            to_0_x_zzz_xxy[k] = -2.0 * to_zzz_xy[k] + 2.0 * to_zzz_xxxy[k] * tke_0;

            to_0_x_zzz_xxz[k] = -2.0 * to_zzz_xz[k] + 2.0 * to_zzz_xxxz[k] * tke_0;

            to_0_x_zzz_xyy[k] = -to_zzz_yy[k] + 2.0 * to_zzz_xxyy[k] * tke_0;

            to_0_x_zzz_xyz[k] = -to_zzz_yz[k] + 2.0 * to_zzz_xxyz[k] * tke_0;

            to_0_x_zzz_xzz[k] = -to_zzz_zz[k] + 2.0 * to_zzz_xxzz[k] * tke_0;

            to_0_x_zzz_yyy[k] = 2.0 * to_zzz_xyyy[k] * tke_0;

            to_0_x_zzz_yyz[k] = 2.0 * to_zzz_xyyz[k] * tke_0;

            to_0_x_zzz_yzz[k] = 2.0 * to_zzz_xyzz[k] * tke_0;

            to_0_x_zzz_zzz[k] = 2.0 * to_zzz_xzzz[k] * tke_0;
        }

        // Set up 100-110 components of targeted buffer : FF

        auto to_0_y_xxx_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 0);

        auto to_0_y_xxx_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 1);

        auto to_0_y_xxx_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 2);

        auto to_0_y_xxx_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 3);

        auto to_0_y_xxx_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 4);

        auto to_0_y_xxx_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 5);

        auto to_0_y_xxx_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 6);

        auto to_0_y_xxx_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 7);

        auto to_0_y_xxx_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 8);

        auto to_0_y_xxx_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_0_y_xxx_xxx, to_0_y_xxx_xxy, to_0_y_xxx_xxz, to_0_y_xxx_xyy, to_0_y_xxx_xyz, to_0_y_xxx_xzz, to_0_y_xxx_yyy, to_0_y_xxx_yyz, to_0_y_xxx_yzz, to_0_y_xxx_zzz, to_xxx_xx, to_xxx_xxxy, to_xxx_xxyy, to_xxx_xxyz, to_xxx_xy, to_xxx_xyyy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_yy, to_xxx_yyyy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxx_xxx[k] = 2.0 * to_xxx_xxxy[k] * tke_0;

            to_0_y_xxx_xxy[k] = -to_xxx_xx[k] + 2.0 * to_xxx_xxyy[k] * tke_0;

            to_0_y_xxx_xxz[k] = 2.0 * to_xxx_xxyz[k] * tke_0;

            to_0_y_xxx_xyy[k] = -2.0 * to_xxx_xy[k] + 2.0 * to_xxx_xyyy[k] * tke_0;

            to_0_y_xxx_xyz[k] = -to_xxx_xz[k] + 2.0 * to_xxx_xyyz[k] * tke_0;

            to_0_y_xxx_xzz[k] = 2.0 * to_xxx_xyzz[k] * tke_0;

            to_0_y_xxx_yyy[k] = -3.0 * to_xxx_yy[k] + 2.0 * to_xxx_yyyy[k] * tke_0;

            to_0_y_xxx_yyz[k] = -2.0 * to_xxx_yz[k] + 2.0 * to_xxx_yyyz[k] * tke_0;

            to_0_y_xxx_yzz[k] = -to_xxx_zz[k] + 2.0 * to_xxx_yyzz[k] * tke_0;

            to_0_y_xxx_zzz[k] = 2.0 * to_xxx_yzzz[k] * tke_0;
        }

        // Set up 110-120 components of targeted buffer : FF

        auto to_0_y_xxy_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 10);

        auto to_0_y_xxy_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 11);

        auto to_0_y_xxy_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 12);

        auto to_0_y_xxy_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 13);

        auto to_0_y_xxy_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 14);

        auto to_0_y_xxy_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 15);

        auto to_0_y_xxy_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 16);

        auto to_0_y_xxy_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 17);

        auto to_0_y_xxy_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 18);

        auto to_0_y_xxy_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_0_y_xxy_xxx, to_0_y_xxy_xxy, to_0_y_xxy_xxz, to_0_y_xxy_xyy, to_0_y_xxy_xyz, to_0_y_xxy_xzz, to_0_y_xxy_yyy, to_0_y_xxy_yyz, to_0_y_xxy_yzz, to_0_y_xxy_zzz, to_xxy_xx, to_xxy_xxxy, to_xxy_xxyy, to_xxy_xxyz, to_xxy_xy, to_xxy_xyyy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_yy, to_xxy_yyyy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxy_xxx[k] = 2.0 * to_xxy_xxxy[k] * tke_0;

            to_0_y_xxy_xxy[k] = -to_xxy_xx[k] + 2.0 * to_xxy_xxyy[k] * tke_0;

            to_0_y_xxy_xxz[k] = 2.0 * to_xxy_xxyz[k] * tke_0;

            to_0_y_xxy_xyy[k] = -2.0 * to_xxy_xy[k] + 2.0 * to_xxy_xyyy[k] * tke_0;

            to_0_y_xxy_xyz[k] = -to_xxy_xz[k] + 2.0 * to_xxy_xyyz[k] * tke_0;

            to_0_y_xxy_xzz[k] = 2.0 * to_xxy_xyzz[k] * tke_0;

            to_0_y_xxy_yyy[k] = -3.0 * to_xxy_yy[k] + 2.0 * to_xxy_yyyy[k] * tke_0;

            to_0_y_xxy_yyz[k] = -2.0 * to_xxy_yz[k] + 2.0 * to_xxy_yyyz[k] * tke_0;

            to_0_y_xxy_yzz[k] = -to_xxy_zz[k] + 2.0 * to_xxy_yyzz[k] * tke_0;

            to_0_y_xxy_zzz[k] = 2.0 * to_xxy_yzzz[k] * tke_0;
        }

        // Set up 120-130 components of targeted buffer : FF

        auto to_0_y_xxz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 20);

        auto to_0_y_xxz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 21);

        auto to_0_y_xxz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 22);

        auto to_0_y_xxz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 23);

        auto to_0_y_xxz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 24);

        auto to_0_y_xxz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 25);

        auto to_0_y_xxz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 26);

        auto to_0_y_xxz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 27);

        auto to_0_y_xxz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 28);

        auto to_0_y_xxz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_0_y_xxz_xxx, to_0_y_xxz_xxy, to_0_y_xxz_xxz, to_0_y_xxz_xyy, to_0_y_xxz_xyz, to_0_y_xxz_xzz, to_0_y_xxz_yyy, to_0_y_xxz_yyz, to_0_y_xxz_yzz, to_0_y_xxz_zzz, to_xxz_xx, to_xxz_xxxy, to_xxz_xxyy, to_xxz_xxyz, to_xxz_xy, to_xxz_xyyy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_yy, to_xxz_yyyy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxz_xxx[k] = 2.0 * to_xxz_xxxy[k] * tke_0;

            to_0_y_xxz_xxy[k] = -to_xxz_xx[k] + 2.0 * to_xxz_xxyy[k] * tke_0;

            to_0_y_xxz_xxz[k] = 2.0 * to_xxz_xxyz[k] * tke_0;

            to_0_y_xxz_xyy[k] = -2.0 * to_xxz_xy[k] + 2.0 * to_xxz_xyyy[k] * tke_0;

            to_0_y_xxz_xyz[k] = -to_xxz_xz[k] + 2.0 * to_xxz_xyyz[k] * tke_0;

            to_0_y_xxz_xzz[k] = 2.0 * to_xxz_xyzz[k] * tke_0;

            to_0_y_xxz_yyy[k] = -3.0 * to_xxz_yy[k] + 2.0 * to_xxz_yyyy[k] * tke_0;

            to_0_y_xxz_yyz[k] = -2.0 * to_xxz_yz[k] + 2.0 * to_xxz_yyyz[k] * tke_0;

            to_0_y_xxz_yzz[k] = -to_xxz_zz[k] + 2.0 * to_xxz_yyzz[k] * tke_0;

            to_0_y_xxz_zzz[k] = 2.0 * to_xxz_yzzz[k] * tke_0;
        }

        // Set up 130-140 components of targeted buffer : FF

        auto to_0_y_xyy_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 30);

        auto to_0_y_xyy_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 31);

        auto to_0_y_xyy_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 32);

        auto to_0_y_xyy_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 33);

        auto to_0_y_xyy_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 34);

        auto to_0_y_xyy_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 35);

        auto to_0_y_xyy_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 36);

        auto to_0_y_xyy_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 37);

        auto to_0_y_xyy_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 38);

        auto to_0_y_xyy_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_0_y_xyy_xxx, to_0_y_xyy_xxy, to_0_y_xyy_xxz, to_0_y_xyy_xyy, to_0_y_xyy_xyz, to_0_y_xyy_xzz, to_0_y_xyy_yyy, to_0_y_xyy_yyz, to_0_y_xyy_yzz, to_0_y_xyy_zzz, to_xyy_xx, to_xyy_xxxy, to_xyy_xxyy, to_xyy_xxyz, to_xyy_xy, to_xyy_xyyy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_yy, to_xyy_yyyy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyy_xxx[k] = 2.0 * to_xyy_xxxy[k] * tke_0;

            to_0_y_xyy_xxy[k] = -to_xyy_xx[k] + 2.0 * to_xyy_xxyy[k] * tke_0;

            to_0_y_xyy_xxz[k] = 2.0 * to_xyy_xxyz[k] * tke_0;

            to_0_y_xyy_xyy[k] = -2.0 * to_xyy_xy[k] + 2.0 * to_xyy_xyyy[k] * tke_0;

            to_0_y_xyy_xyz[k] = -to_xyy_xz[k] + 2.0 * to_xyy_xyyz[k] * tke_0;

            to_0_y_xyy_xzz[k] = 2.0 * to_xyy_xyzz[k] * tke_0;

            to_0_y_xyy_yyy[k] = -3.0 * to_xyy_yy[k] + 2.0 * to_xyy_yyyy[k] * tke_0;

            to_0_y_xyy_yyz[k] = -2.0 * to_xyy_yz[k] + 2.0 * to_xyy_yyyz[k] * tke_0;

            to_0_y_xyy_yzz[k] = -to_xyy_zz[k] + 2.0 * to_xyy_yyzz[k] * tke_0;

            to_0_y_xyy_zzz[k] = 2.0 * to_xyy_yzzz[k] * tke_0;
        }

        // Set up 140-150 components of targeted buffer : FF

        auto to_0_y_xyz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 40);

        auto to_0_y_xyz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 41);

        auto to_0_y_xyz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 42);

        auto to_0_y_xyz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 43);

        auto to_0_y_xyz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 44);

        auto to_0_y_xyz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 45);

        auto to_0_y_xyz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 46);

        auto to_0_y_xyz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 47);

        auto to_0_y_xyz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 48);

        auto to_0_y_xyz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_0_y_xyz_xxx, to_0_y_xyz_xxy, to_0_y_xyz_xxz, to_0_y_xyz_xyy, to_0_y_xyz_xyz, to_0_y_xyz_xzz, to_0_y_xyz_yyy, to_0_y_xyz_yyz, to_0_y_xyz_yzz, to_0_y_xyz_zzz, to_xyz_xx, to_xyz_xxxy, to_xyz_xxyy, to_xyz_xxyz, to_xyz_xy, to_xyz_xyyy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_yy, to_xyz_yyyy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyz_xxx[k] = 2.0 * to_xyz_xxxy[k] * tke_0;

            to_0_y_xyz_xxy[k] = -to_xyz_xx[k] + 2.0 * to_xyz_xxyy[k] * tke_0;

            to_0_y_xyz_xxz[k] = 2.0 * to_xyz_xxyz[k] * tke_0;

            to_0_y_xyz_xyy[k] = -2.0 * to_xyz_xy[k] + 2.0 * to_xyz_xyyy[k] * tke_0;

            to_0_y_xyz_xyz[k] = -to_xyz_xz[k] + 2.0 * to_xyz_xyyz[k] * tke_0;

            to_0_y_xyz_xzz[k] = 2.0 * to_xyz_xyzz[k] * tke_0;

            to_0_y_xyz_yyy[k] = -3.0 * to_xyz_yy[k] + 2.0 * to_xyz_yyyy[k] * tke_0;

            to_0_y_xyz_yyz[k] = -2.0 * to_xyz_yz[k] + 2.0 * to_xyz_yyyz[k] * tke_0;

            to_0_y_xyz_yzz[k] = -to_xyz_zz[k] + 2.0 * to_xyz_yyzz[k] * tke_0;

            to_0_y_xyz_zzz[k] = 2.0 * to_xyz_yzzz[k] * tke_0;
        }

        // Set up 150-160 components of targeted buffer : FF

        auto to_0_y_xzz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 50);

        auto to_0_y_xzz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 51);

        auto to_0_y_xzz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 52);

        auto to_0_y_xzz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 53);

        auto to_0_y_xzz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 54);

        auto to_0_y_xzz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 55);

        auto to_0_y_xzz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 56);

        auto to_0_y_xzz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 57);

        auto to_0_y_xzz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 58);

        auto to_0_y_xzz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_0_y_xzz_xxx, to_0_y_xzz_xxy, to_0_y_xzz_xxz, to_0_y_xzz_xyy, to_0_y_xzz_xyz, to_0_y_xzz_xzz, to_0_y_xzz_yyy, to_0_y_xzz_yyz, to_0_y_xzz_yzz, to_0_y_xzz_zzz, to_xzz_xx, to_xzz_xxxy, to_xzz_xxyy, to_xzz_xxyz, to_xzz_xy, to_xzz_xyyy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_yy, to_xzz_yyyy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzz_xxx[k] = 2.0 * to_xzz_xxxy[k] * tke_0;

            to_0_y_xzz_xxy[k] = -to_xzz_xx[k] + 2.0 * to_xzz_xxyy[k] * tke_0;

            to_0_y_xzz_xxz[k] = 2.0 * to_xzz_xxyz[k] * tke_0;

            to_0_y_xzz_xyy[k] = -2.0 * to_xzz_xy[k] + 2.0 * to_xzz_xyyy[k] * tke_0;

            to_0_y_xzz_xyz[k] = -to_xzz_xz[k] + 2.0 * to_xzz_xyyz[k] * tke_0;

            to_0_y_xzz_xzz[k] = 2.0 * to_xzz_xyzz[k] * tke_0;

            to_0_y_xzz_yyy[k] = -3.0 * to_xzz_yy[k] + 2.0 * to_xzz_yyyy[k] * tke_0;

            to_0_y_xzz_yyz[k] = -2.0 * to_xzz_yz[k] + 2.0 * to_xzz_yyyz[k] * tke_0;

            to_0_y_xzz_yzz[k] = -to_xzz_zz[k] + 2.0 * to_xzz_yyzz[k] * tke_0;

            to_0_y_xzz_zzz[k] = 2.0 * to_xzz_yzzz[k] * tke_0;
        }

        // Set up 160-170 components of targeted buffer : FF

        auto to_0_y_yyy_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 60);

        auto to_0_y_yyy_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 61);

        auto to_0_y_yyy_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 62);

        auto to_0_y_yyy_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 63);

        auto to_0_y_yyy_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 64);

        auto to_0_y_yyy_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 65);

        auto to_0_y_yyy_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 66);

        auto to_0_y_yyy_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 67);

        auto to_0_y_yyy_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 68);

        auto to_0_y_yyy_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_0_y_yyy_xxx, to_0_y_yyy_xxy, to_0_y_yyy_xxz, to_0_y_yyy_xyy, to_0_y_yyy_xyz, to_0_y_yyy_xzz, to_0_y_yyy_yyy, to_0_y_yyy_yyz, to_0_y_yyy_yzz, to_0_y_yyy_zzz, to_yyy_xx, to_yyy_xxxy, to_yyy_xxyy, to_yyy_xxyz, to_yyy_xy, to_yyy_xyyy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_yy, to_yyy_yyyy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyy_xxx[k] = 2.0 * to_yyy_xxxy[k] * tke_0;

            to_0_y_yyy_xxy[k] = -to_yyy_xx[k] + 2.0 * to_yyy_xxyy[k] * tke_0;

            to_0_y_yyy_xxz[k] = 2.0 * to_yyy_xxyz[k] * tke_0;

            to_0_y_yyy_xyy[k] = -2.0 * to_yyy_xy[k] + 2.0 * to_yyy_xyyy[k] * tke_0;

            to_0_y_yyy_xyz[k] = -to_yyy_xz[k] + 2.0 * to_yyy_xyyz[k] * tke_0;

            to_0_y_yyy_xzz[k] = 2.0 * to_yyy_xyzz[k] * tke_0;

            to_0_y_yyy_yyy[k] = -3.0 * to_yyy_yy[k] + 2.0 * to_yyy_yyyy[k] * tke_0;

            to_0_y_yyy_yyz[k] = -2.0 * to_yyy_yz[k] + 2.0 * to_yyy_yyyz[k] * tke_0;

            to_0_y_yyy_yzz[k] = -to_yyy_zz[k] + 2.0 * to_yyy_yyzz[k] * tke_0;

            to_0_y_yyy_zzz[k] = 2.0 * to_yyy_yzzz[k] * tke_0;
        }

        // Set up 170-180 components of targeted buffer : FF

        auto to_0_y_yyz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 70);

        auto to_0_y_yyz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 71);

        auto to_0_y_yyz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 72);

        auto to_0_y_yyz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 73);

        auto to_0_y_yyz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 74);

        auto to_0_y_yyz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 75);

        auto to_0_y_yyz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 76);

        auto to_0_y_yyz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 77);

        auto to_0_y_yyz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 78);

        auto to_0_y_yyz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_0_y_yyz_xxx, to_0_y_yyz_xxy, to_0_y_yyz_xxz, to_0_y_yyz_xyy, to_0_y_yyz_xyz, to_0_y_yyz_xzz, to_0_y_yyz_yyy, to_0_y_yyz_yyz, to_0_y_yyz_yzz, to_0_y_yyz_zzz, to_yyz_xx, to_yyz_xxxy, to_yyz_xxyy, to_yyz_xxyz, to_yyz_xy, to_yyz_xyyy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_yy, to_yyz_yyyy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyz_xxx[k] = 2.0 * to_yyz_xxxy[k] * tke_0;

            to_0_y_yyz_xxy[k] = -to_yyz_xx[k] + 2.0 * to_yyz_xxyy[k] * tke_0;

            to_0_y_yyz_xxz[k] = 2.0 * to_yyz_xxyz[k] * tke_0;

            to_0_y_yyz_xyy[k] = -2.0 * to_yyz_xy[k] + 2.0 * to_yyz_xyyy[k] * tke_0;

            to_0_y_yyz_xyz[k] = -to_yyz_xz[k] + 2.0 * to_yyz_xyyz[k] * tke_0;

            to_0_y_yyz_xzz[k] = 2.0 * to_yyz_xyzz[k] * tke_0;

            to_0_y_yyz_yyy[k] = -3.0 * to_yyz_yy[k] + 2.0 * to_yyz_yyyy[k] * tke_0;

            to_0_y_yyz_yyz[k] = -2.0 * to_yyz_yz[k] + 2.0 * to_yyz_yyyz[k] * tke_0;

            to_0_y_yyz_yzz[k] = -to_yyz_zz[k] + 2.0 * to_yyz_yyzz[k] * tke_0;

            to_0_y_yyz_zzz[k] = 2.0 * to_yyz_yzzz[k] * tke_0;
        }

        // Set up 180-190 components of targeted buffer : FF

        auto to_0_y_yzz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 80);

        auto to_0_y_yzz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 81);

        auto to_0_y_yzz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 82);

        auto to_0_y_yzz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 83);

        auto to_0_y_yzz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 84);

        auto to_0_y_yzz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 85);

        auto to_0_y_yzz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 86);

        auto to_0_y_yzz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 87);

        auto to_0_y_yzz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 88);

        auto to_0_y_yzz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_0_y_yzz_xxx, to_0_y_yzz_xxy, to_0_y_yzz_xxz, to_0_y_yzz_xyy, to_0_y_yzz_xyz, to_0_y_yzz_xzz, to_0_y_yzz_yyy, to_0_y_yzz_yyz, to_0_y_yzz_yzz, to_0_y_yzz_zzz, to_yzz_xx, to_yzz_xxxy, to_yzz_xxyy, to_yzz_xxyz, to_yzz_xy, to_yzz_xyyy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_yy, to_yzz_yyyy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzz_xxx[k] = 2.0 * to_yzz_xxxy[k] * tke_0;

            to_0_y_yzz_xxy[k] = -to_yzz_xx[k] + 2.0 * to_yzz_xxyy[k] * tke_0;

            to_0_y_yzz_xxz[k] = 2.0 * to_yzz_xxyz[k] * tke_0;

            to_0_y_yzz_xyy[k] = -2.0 * to_yzz_xy[k] + 2.0 * to_yzz_xyyy[k] * tke_0;

            to_0_y_yzz_xyz[k] = -to_yzz_xz[k] + 2.0 * to_yzz_xyyz[k] * tke_0;

            to_0_y_yzz_xzz[k] = 2.0 * to_yzz_xyzz[k] * tke_0;

            to_0_y_yzz_yyy[k] = -3.0 * to_yzz_yy[k] + 2.0 * to_yzz_yyyy[k] * tke_0;

            to_0_y_yzz_yyz[k] = -2.0 * to_yzz_yz[k] + 2.0 * to_yzz_yyyz[k] * tke_0;

            to_0_y_yzz_yzz[k] = -to_yzz_zz[k] + 2.0 * to_yzz_yyzz[k] * tke_0;

            to_0_y_yzz_zzz[k] = 2.0 * to_yzz_yzzz[k] * tke_0;
        }

        // Set up 190-200 components of targeted buffer : FF

        auto to_0_y_zzz_xxx = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 90);

        auto to_0_y_zzz_xxy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 91);

        auto to_0_y_zzz_xxz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 92);

        auto to_0_y_zzz_xyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 93);

        auto to_0_y_zzz_xyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 94);

        auto to_0_y_zzz_xzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 95);

        auto to_0_y_zzz_yyy = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 96);

        auto to_0_y_zzz_yyz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 97);

        auto to_0_y_zzz_yzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 98);

        auto to_0_y_zzz_zzz = pbuffer.data(idx_op_geom_001_ff + 1 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_0_y_zzz_xxx, to_0_y_zzz_xxy, to_0_y_zzz_xxz, to_0_y_zzz_xyy, to_0_y_zzz_xyz, to_0_y_zzz_xzz, to_0_y_zzz_yyy, to_0_y_zzz_yyz, to_0_y_zzz_yzz, to_0_y_zzz_zzz, to_zzz_xx, to_zzz_xxxy, to_zzz_xxyy, to_zzz_xxyz, to_zzz_xy, to_zzz_xyyy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_yy, to_zzz_yyyy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzz_xxx[k] = 2.0 * to_zzz_xxxy[k] * tke_0;

            to_0_y_zzz_xxy[k] = -to_zzz_xx[k] + 2.0 * to_zzz_xxyy[k] * tke_0;

            to_0_y_zzz_xxz[k] = 2.0 * to_zzz_xxyz[k] * tke_0;

            to_0_y_zzz_xyy[k] = -2.0 * to_zzz_xy[k] + 2.0 * to_zzz_xyyy[k] * tke_0;

            to_0_y_zzz_xyz[k] = -to_zzz_xz[k] + 2.0 * to_zzz_xyyz[k] * tke_0;

            to_0_y_zzz_xzz[k] = 2.0 * to_zzz_xyzz[k] * tke_0;

            to_0_y_zzz_yyy[k] = -3.0 * to_zzz_yy[k] + 2.0 * to_zzz_yyyy[k] * tke_0;

            to_0_y_zzz_yyz[k] = -2.0 * to_zzz_yz[k] + 2.0 * to_zzz_yyyz[k] * tke_0;

            to_0_y_zzz_yzz[k] = -to_zzz_zz[k] + 2.0 * to_zzz_yyzz[k] * tke_0;

            to_0_y_zzz_zzz[k] = 2.0 * to_zzz_yzzz[k] * tke_0;
        }

        // Set up 200-210 components of targeted buffer : FF

        auto to_0_z_xxx_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 0);

        auto to_0_z_xxx_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 1);

        auto to_0_z_xxx_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 2);

        auto to_0_z_xxx_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 3);

        auto to_0_z_xxx_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 4);

        auto to_0_z_xxx_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 5);

        auto to_0_z_xxx_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 6);

        auto to_0_z_xxx_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 7);

        auto to_0_z_xxx_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 8);

        auto to_0_z_xxx_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_0_z_xxx_xxx, to_0_z_xxx_xxy, to_0_z_xxx_xxz, to_0_z_xxx_xyy, to_0_z_xxx_xyz, to_0_z_xxx_xzz, to_0_z_xxx_yyy, to_0_z_xxx_yyz, to_0_z_xxx_yzz, to_0_z_xxx_zzz, to_xxx_xx, to_xxx_xxxz, to_xxx_xxyz, to_xxx_xxzz, to_xxx_xy, to_xxx_xyyz, to_xxx_xyzz, to_xxx_xz, to_xxx_xzzz, to_xxx_yy, to_xxx_yyyz, to_xxx_yyzz, to_xxx_yz, to_xxx_yzzz, to_xxx_zz, to_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxx_xxx[k] = 2.0 * to_xxx_xxxz[k] * tke_0;

            to_0_z_xxx_xxy[k] = 2.0 * to_xxx_xxyz[k] * tke_0;

            to_0_z_xxx_xxz[k] = -to_xxx_xx[k] + 2.0 * to_xxx_xxzz[k] * tke_0;

            to_0_z_xxx_xyy[k] = 2.0 * to_xxx_xyyz[k] * tke_0;

            to_0_z_xxx_xyz[k] = -to_xxx_xy[k] + 2.0 * to_xxx_xyzz[k] * tke_0;

            to_0_z_xxx_xzz[k] = -2.0 * to_xxx_xz[k] + 2.0 * to_xxx_xzzz[k] * tke_0;

            to_0_z_xxx_yyy[k] = 2.0 * to_xxx_yyyz[k] * tke_0;

            to_0_z_xxx_yyz[k] = -to_xxx_yy[k] + 2.0 * to_xxx_yyzz[k] * tke_0;

            to_0_z_xxx_yzz[k] = -2.0 * to_xxx_yz[k] + 2.0 * to_xxx_yzzz[k] * tke_0;

            to_0_z_xxx_zzz[k] = -3.0 * to_xxx_zz[k] + 2.0 * to_xxx_zzzz[k] * tke_0;
        }

        // Set up 210-220 components of targeted buffer : FF

        auto to_0_z_xxy_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 10);

        auto to_0_z_xxy_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 11);

        auto to_0_z_xxy_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 12);

        auto to_0_z_xxy_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 13);

        auto to_0_z_xxy_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 14);

        auto to_0_z_xxy_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 15);

        auto to_0_z_xxy_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 16);

        auto to_0_z_xxy_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 17);

        auto to_0_z_xxy_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 18);

        auto to_0_z_xxy_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_0_z_xxy_xxx, to_0_z_xxy_xxy, to_0_z_xxy_xxz, to_0_z_xxy_xyy, to_0_z_xxy_xyz, to_0_z_xxy_xzz, to_0_z_xxy_yyy, to_0_z_xxy_yyz, to_0_z_xxy_yzz, to_0_z_xxy_zzz, to_xxy_xx, to_xxy_xxxz, to_xxy_xxyz, to_xxy_xxzz, to_xxy_xy, to_xxy_xyyz, to_xxy_xyzz, to_xxy_xz, to_xxy_xzzz, to_xxy_yy, to_xxy_yyyz, to_xxy_yyzz, to_xxy_yz, to_xxy_yzzz, to_xxy_zz, to_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxy_xxx[k] = 2.0 * to_xxy_xxxz[k] * tke_0;

            to_0_z_xxy_xxy[k] = 2.0 * to_xxy_xxyz[k] * tke_0;

            to_0_z_xxy_xxz[k] = -to_xxy_xx[k] + 2.0 * to_xxy_xxzz[k] * tke_0;

            to_0_z_xxy_xyy[k] = 2.0 * to_xxy_xyyz[k] * tke_0;

            to_0_z_xxy_xyz[k] = -to_xxy_xy[k] + 2.0 * to_xxy_xyzz[k] * tke_0;

            to_0_z_xxy_xzz[k] = -2.0 * to_xxy_xz[k] + 2.0 * to_xxy_xzzz[k] * tke_0;

            to_0_z_xxy_yyy[k] = 2.0 * to_xxy_yyyz[k] * tke_0;

            to_0_z_xxy_yyz[k] = -to_xxy_yy[k] + 2.0 * to_xxy_yyzz[k] * tke_0;

            to_0_z_xxy_yzz[k] = -2.0 * to_xxy_yz[k] + 2.0 * to_xxy_yzzz[k] * tke_0;

            to_0_z_xxy_zzz[k] = -3.0 * to_xxy_zz[k] + 2.0 * to_xxy_zzzz[k] * tke_0;
        }

        // Set up 220-230 components of targeted buffer : FF

        auto to_0_z_xxz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 20);

        auto to_0_z_xxz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 21);

        auto to_0_z_xxz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 22);

        auto to_0_z_xxz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 23);

        auto to_0_z_xxz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 24);

        auto to_0_z_xxz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 25);

        auto to_0_z_xxz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 26);

        auto to_0_z_xxz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 27);

        auto to_0_z_xxz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 28);

        auto to_0_z_xxz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_0_z_xxz_xxx, to_0_z_xxz_xxy, to_0_z_xxz_xxz, to_0_z_xxz_xyy, to_0_z_xxz_xyz, to_0_z_xxz_xzz, to_0_z_xxz_yyy, to_0_z_xxz_yyz, to_0_z_xxz_yzz, to_0_z_xxz_zzz, to_xxz_xx, to_xxz_xxxz, to_xxz_xxyz, to_xxz_xxzz, to_xxz_xy, to_xxz_xyyz, to_xxz_xyzz, to_xxz_xz, to_xxz_xzzz, to_xxz_yy, to_xxz_yyyz, to_xxz_yyzz, to_xxz_yz, to_xxz_yzzz, to_xxz_zz, to_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxz_xxx[k] = 2.0 * to_xxz_xxxz[k] * tke_0;

            to_0_z_xxz_xxy[k] = 2.0 * to_xxz_xxyz[k] * tke_0;

            to_0_z_xxz_xxz[k] = -to_xxz_xx[k] + 2.0 * to_xxz_xxzz[k] * tke_0;

            to_0_z_xxz_xyy[k] = 2.0 * to_xxz_xyyz[k] * tke_0;

            to_0_z_xxz_xyz[k] = -to_xxz_xy[k] + 2.0 * to_xxz_xyzz[k] * tke_0;

            to_0_z_xxz_xzz[k] = -2.0 * to_xxz_xz[k] + 2.0 * to_xxz_xzzz[k] * tke_0;

            to_0_z_xxz_yyy[k] = 2.0 * to_xxz_yyyz[k] * tke_0;

            to_0_z_xxz_yyz[k] = -to_xxz_yy[k] + 2.0 * to_xxz_yyzz[k] * tke_0;

            to_0_z_xxz_yzz[k] = -2.0 * to_xxz_yz[k] + 2.0 * to_xxz_yzzz[k] * tke_0;

            to_0_z_xxz_zzz[k] = -3.0 * to_xxz_zz[k] + 2.0 * to_xxz_zzzz[k] * tke_0;
        }

        // Set up 230-240 components of targeted buffer : FF

        auto to_0_z_xyy_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 30);

        auto to_0_z_xyy_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 31);

        auto to_0_z_xyy_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 32);

        auto to_0_z_xyy_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 33);

        auto to_0_z_xyy_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 34);

        auto to_0_z_xyy_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 35);

        auto to_0_z_xyy_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 36);

        auto to_0_z_xyy_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 37);

        auto to_0_z_xyy_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 38);

        auto to_0_z_xyy_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_0_z_xyy_xxx, to_0_z_xyy_xxy, to_0_z_xyy_xxz, to_0_z_xyy_xyy, to_0_z_xyy_xyz, to_0_z_xyy_xzz, to_0_z_xyy_yyy, to_0_z_xyy_yyz, to_0_z_xyy_yzz, to_0_z_xyy_zzz, to_xyy_xx, to_xyy_xxxz, to_xyy_xxyz, to_xyy_xxzz, to_xyy_xy, to_xyy_xyyz, to_xyy_xyzz, to_xyy_xz, to_xyy_xzzz, to_xyy_yy, to_xyy_yyyz, to_xyy_yyzz, to_xyy_yz, to_xyy_yzzz, to_xyy_zz, to_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyy_xxx[k] = 2.0 * to_xyy_xxxz[k] * tke_0;

            to_0_z_xyy_xxy[k] = 2.0 * to_xyy_xxyz[k] * tke_0;

            to_0_z_xyy_xxz[k] = -to_xyy_xx[k] + 2.0 * to_xyy_xxzz[k] * tke_0;

            to_0_z_xyy_xyy[k] = 2.0 * to_xyy_xyyz[k] * tke_0;

            to_0_z_xyy_xyz[k] = -to_xyy_xy[k] + 2.0 * to_xyy_xyzz[k] * tke_0;

            to_0_z_xyy_xzz[k] = -2.0 * to_xyy_xz[k] + 2.0 * to_xyy_xzzz[k] * tke_0;

            to_0_z_xyy_yyy[k] = 2.0 * to_xyy_yyyz[k] * tke_0;

            to_0_z_xyy_yyz[k] = -to_xyy_yy[k] + 2.0 * to_xyy_yyzz[k] * tke_0;

            to_0_z_xyy_yzz[k] = -2.0 * to_xyy_yz[k] + 2.0 * to_xyy_yzzz[k] * tke_0;

            to_0_z_xyy_zzz[k] = -3.0 * to_xyy_zz[k] + 2.0 * to_xyy_zzzz[k] * tke_0;
        }

        // Set up 240-250 components of targeted buffer : FF

        auto to_0_z_xyz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 40);

        auto to_0_z_xyz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 41);

        auto to_0_z_xyz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 42);

        auto to_0_z_xyz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 43);

        auto to_0_z_xyz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 44);

        auto to_0_z_xyz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 45);

        auto to_0_z_xyz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 46);

        auto to_0_z_xyz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 47);

        auto to_0_z_xyz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 48);

        auto to_0_z_xyz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_0_z_xyz_xxx, to_0_z_xyz_xxy, to_0_z_xyz_xxz, to_0_z_xyz_xyy, to_0_z_xyz_xyz, to_0_z_xyz_xzz, to_0_z_xyz_yyy, to_0_z_xyz_yyz, to_0_z_xyz_yzz, to_0_z_xyz_zzz, to_xyz_xx, to_xyz_xxxz, to_xyz_xxyz, to_xyz_xxzz, to_xyz_xy, to_xyz_xyyz, to_xyz_xyzz, to_xyz_xz, to_xyz_xzzz, to_xyz_yy, to_xyz_yyyz, to_xyz_yyzz, to_xyz_yz, to_xyz_yzzz, to_xyz_zz, to_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyz_xxx[k] = 2.0 * to_xyz_xxxz[k] * tke_0;

            to_0_z_xyz_xxy[k] = 2.0 * to_xyz_xxyz[k] * tke_0;

            to_0_z_xyz_xxz[k] = -to_xyz_xx[k] + 2.0 * to_xyz_xxzz[k] * tke_0;

            to_0_z_xyz_xyy[k] = 2.0 * to_xyz_xyyz[k] * tke_0;

            to_0_z_xyz_xyz[k] = -to_xyz_xy[k] + 2.0 * to_xyz_xyzz[k] * tke_0;

            to_0_z_xyz_xzz[k] = -2.0 * to_xyz_xz[k] + 2.0 * to_xyz_xzzz[k] * tke_0;

            to_0_z_xyz_yyy[k] = 2.0 * to_xyz_yyyz[k] * tke_0;

            to_0_z_xyz_yyz[k] = -to_xyz_yy[k] + 2.0 * to_xyz_yyzz[k] * tke_0;

            to_0_z_xyz_yzz[k] = -2.0 * to_xyz_yz[k] + 2.0 * to_xyz_yzzz[k] * tke_0;

            to_0_z_xyz_zzz[k] = -3.0 * to_xyz_zz[k] + 2.0 * to_xyz_zzzz[k] * tke_0;
        }

        // Set up 250-260 components of targeted buffer : FF

        auto to_0_z_xzz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 50);

        auto to_0_z_xzz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 51);

        auto to_0_z_xzz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 52);

        auto to_0_z_xzz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 53);

        auto to_0_z_xzz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 54);

        auto to_0_z_xzz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 55);

        auto to_0_z_xzz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 56);

        auto to_0_z_xzz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 57);

        auto to_0_z_xzz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 58);

        auto to_0_z_xzz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_0_z_xzz_xxx, to_0_z_xzz_xxy, to_0_z_xzz_xxz, to_0_z_xzz_xyy, to_0_z_xzz_xyz, to_0_z_xzz_xzz, to_0_z_xzz_yyy, to_0_z_xzz_yyz, to_0_z_xzz_yzz, to_0_z_xzz_zzz, to_xzz_xx, to_xzz_xxxz, to_xzz_xxyz, to_xzz_xxzz, to_xzz_xy, to_xzz_xyyz, to_xzz_xyzz, to_xzz_xz, to_xzz_xzzz, to_xzz_yy, to_xzz_yyyz, to_xzz_yyzz, to_xzz_yz, to_xzz_yzzz, to_xzz_zz, to_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzz_xxx[k] = 2.0 * to_xzz_xxxz[k] * tke_0;

            to_0_z_xzz_xxy[k] = 2.0 * to_xzz_xxyz[k] * tke_0;

            to_0_z_xzz_xxz[k] = -to_xzz_xx[k] + 2.0 * to_xzz_xxzz[k] * tke_0;

            to_0_z_xzz_xyy[k] = 2.0 * to_xzz_xyyz[k] * tke_0;

            to_0_z_xzz_xyz[k] = -to_xzz_xy[k] + 2.0 * to_xzz_xyzz[k] * tke_0;

            to_0_z_xzz_xzz[k] = -2.0 * to_xzz_xz[k] + 2.0 * to_xzz_xzzz[k] * tke_0;

            to_0_z_xzz_yyy[k] = 2.0 * to_xzz_yyyz[k] * tke_0;

            to_0_z_xzz_yyz[k] = -to_xzz_yy[k] + 2.0 * to_xzz_yyzz[k] * tke_0;

            to_0_z_xzz_yzz[k] = -2.0 * to_xzz_yz[k] + 2.0 * to_xzz_yzzz[k] * tke_0;

            to_0_z_xzz_zzz[k] = -3.0 * to_xzz_zz[k] + 2.0 * to_xzz_zzzz[k] * tke_0;
        }

        // Set up 260-270 components of targeted buffer : FF

        auto to_0_z_yyy_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 60);

        auto to_0_z_yyy_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 61);

        auto to_0_z_yyy_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 62);

        auto to_0_z_yyy_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 63);

        auto to_0_z_yyy_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 64);

        auto to_0_z_yyy_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 65);

        auto to_0_z_yyy_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 66);

        auto to_0_z_yyy_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 67);

        auto to_0_z_yyy_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 68);

        auto to_0_z_yyy_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_0_z_yyy_xxx, to_0_z_yyy_xxy, to_0_z_yyy_xxz, to_0_z_yyy_xyy, to_0_z_yyy_xyz, to_0_z_yyy_xzz, to_0_z_yyy_yyy, to_0_z_yyy_yyz, to_0_z_yyy_yzz, to_0_z_yyy_zzz, to_yyy_xx, to_yyy_xxxz, to_yyy_xxyz, to_yyy_xxzz, to_yyy_xy, to_yyy_xyyz, to_yyy_xyzz, to_yyy_xz, to_yyy_xzzz, to_yyy_yy, to_yyy_yyyz, to_yyy_yyzz, to_yyy_yz, to_yyy_yzzz, to_yyy_zz, to_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyy_xxx[k] = 2.0 * to_yyy_xxxz[k] * tke_0;

            to_0_z_yyy_xxy[k] = 2.0 * to_yyy_xxyz[k] * tke_0;

            to_0_z_yyy_xxz[k] = -to_yyy_xx[k] + 2.0 * to_yyy_xxzz[k] * tke_0;

            to_0_z_yyy_xyy[k] = 2.0 * to_yyy_xyyz[k] * tke_0;

            to_0_z_yyy_xyz[k] = -to_yyy_xy[k] + 2.0 * to_yyy_xyzz[k] * tke_0;

            to_0_z_yyy_xzz[k] = -2.0 * to_yyy_xz[k] + 2.0 * to_yyy_xzzz[k] * tke_0;

            to_0_z_yyy_yyy[k] = 2.0 * to_yyy_yyyz[k] * tke_0;

            to_0_z_yyy_yyz[k] = -to_yyy_yy[k] + 2.0 * to_yyy_yyzz[k] * tke_0;

            to_0_z_yyy_yzz[k] = -2.0 * to_yyy_yz[k] + 2.0 * to_yyy_yzzz[k] * tke_0;

            to_0_z_yyy_zzz[k] = -3.0 * to_yyy_zz[k] + 2.0 * to_yyy_zzzz[k] * tke_0;
        }

        // Set up 270-280 components of targeted buffer : FF

        auto to_0_z_yyz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 70);

        auto to_0_z_yyz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 71);

        auto to_0_z_yyz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 72);

        auto to_0_z_yyz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 73);

        auto to_0_z_yyz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 74);

        auto to_0_z_yyz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 75);

        auto to_0_z_yyz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 76);

        auto to_0_z_yyz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 77);

        auto to_0_z_yyz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 78);

        auto to_0_z_yyz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_0_z_yyz_xxx, to_0_z_yyz_xxy, to_0_z_yyz_xxz, to_0_z_yyz_xyy, to_0_z_yyz_xyz, to_0_z_yyz_xzz, to_0_z_yyz_yyy, to_0_z_yyz_yyz, to_0_z_yyz_yzz, to_0_z_yyz_zzz, to_yyz_xx, to_yyz_xxxz, to_yyz_xxyz, to_yyz_xxzz, to_yyz_xy, to_yyz_xyyz, to_yyz_xyzz, to_yyz_xz, to_yyz_xzzz, to_yyz_yy, to_yyz_yyyz, to_yyz_yyzz, to_yyz_yz, to_yyz_yzzz, to_yyz_zz, to_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyz_xxx[k] = 2.0 * to_yyz_xxxz[k] * tke_0;

            to_0_z_yyz_xxy[k] = 2.0 * to_yyz_xxyz[k] * tke_0;

            to_0_z_yyz_xxz[k] = -to_yyz_xx[k] + 2.0 * to_yyz_xxzz[k] * tke_0;

            to_0_z_yyz_xyy[k] = 2.0 * to_yyz_xyyz[k] * tke_0;

            to_0_z_yyz_xyz[k] = -to_yyz_xy[k] + 2.0 * to_yyz_xyzz[k] * tke_0;

            to_0_z_yyz_xzz[k] = -2.0 * to_yyz_xz[k] + 2.0 * to_yyz_xzzz[k] * tke_0;

            to_0_z_yyz_yyy[k] = 2.0 * to_yyz_yyyz[k] * tke_0;

            to_0_z_yyz_yyz[k] = -to_yyz_yy[k] + 2.0 * to_yyz_yyzz[k] * tke_0;

            to_0_z_yyz_yzz[k] = -2.0 * to_yyz_yz[k] + 2.0 * to_yyz_yzzz[k] * tke_0;

            to_0_z_yyz_zzz[k] = -3.0 * to_yyz_zz[k] + 2.0 * to_yyz_zzzz[k] * tke_0;
        }

        // Set up 280-290 components of targeted buffer : FF

        auto to_0_z_yzz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 80);

        auto to_0_z_yzz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 81);

        auto to_0_z_yzz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 82);

        auto to_0_z_yzz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 83);

        auto to_0_z_yzz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 84);

        auto to_0_z_yzz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 85);

        auto to_0_z_yzz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 86);

        auto to_0_z_yzz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 87);

        auto to_0_z_yzz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 88);

        auto to_0_z_yzz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_0_z_yzz_xxx, to_0_z_yzz_xxy, to_0_z_yzz_xxz, to_0_z_yzz_xyy, to_0_z_yzz_xyz, to_0_z_yzz_xzz, to_0_z_yzz_yyy, to_0_z_yzz_yyz, to_0_z_yzz_yzz, to_0_z_yzz_zzz, to_yzz_xx, to_yzz_xxxz, to_yzz_xxyz, to_yzz_xxzz, to_yzz_xy, to_yzz_xyyz, to_yzz_xyzz, to_yzz_xz, to_yzz_xzzz, to_yzz_yy, to_yzz_yyyz, to_yzz_yyzz, to_yzz_yz, to_yzz_yzzz, to_yzz_zz, to_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzz_xxx[k] = 2.0 * to_yzz_xxxz[k] * tke_0;

            to_0_z_yzz_xxy[k] = 2.0 * to_yzz_xxyz[k] * tke_0;

            to_0_z_yzz_xxz[k] = -to_yzz_xx[k] + 2.0 * to_yzz_xxzz[k] * tke_0;

            to_0_z_yzz_xyy[k] = 2.0 * to_yzz_xyyz[k] * tke_0;

            to_0_z_yzz_xyz[k] = -to_yzz_xy[k] + 2.0 * to_yzz_xyzz[k] * tke_0;

            to_0_z_yzz_xzz[k] = -2.0 * to_yzz_xz[k] + 2.0 * to_yzz_xzzz[k] * tke_0;

            to_0_z_yzz_yyy[k] = 2.0 * to_yzz_yyyz[k] * tke_0;

            to_0_z_yzz_yyz[k] = -to_yzz_yy[k] + 2.0 * to_yzz_yyzz[k] * tke_0;

            to_0_z_yzz_yzz[k] = -2.0 * to_yzz_yz[k] + 2.0 * to_yzz_yzzz[k] * tke_0;

            to_0_z_yzz_zzz[k] = -3.0 * to_yzz_zz[k] + 2.0 * to_yzz_zzzz[k] * tke_0;
        }

        // Set up 290-300 components of targeted buffer : FF

        auto to_0_z_zzz_xxx = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 90);

        auto to_0_z_zzz_xxy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 91);

        auto to_0_z_zzz_xxz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 92);

        auto to_0_z_zzz_xyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 93);

        auto to_0_z_zzz_xyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 94);

        auto to_0_z_zzz_xzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 95);

        auto to_0_z_zzz_yyy = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 96);

        auto to_0_z_zzz_yyz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 97);

        auto to_0_z_zzz_yzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 98);

        auto to_0_z_zzz_zzz = pbuffer.data(idx_op_geom_001_ff + 2 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_0_z_zzz_xxx, to_0_z_zzz_xxy, to_0_z_zzz_xxz, to_0_z_zzz_xyy, to_0_z_zzz_xyz, to_0_z_zzz_xzz, to_0_z_zzz_yyy, to_0_z_zzz_yyz, to_0_z_zzz_yzz, to_0_z_zzz_zzz, to_zzz_xx, to_zzz_xxxz, to_zzz_xxyz, to_zzz_xxzz, to_zzz_xy, to_zzz_xyyz, to_zzz_xyzz, to_zzz_xz, to_zzz_xzzz, to_zzz_yy, to_zzz_yyyz, to_zzz_yyzz, to_zzz_yz, to_zzz_yzzz, to_zzz_zz, to_zzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzz_xxx[k] = 2.0 * to_zzz_xxxz[k] * tke_0;

            to_0_z_zzz_xxy[k] = 2.0 * to_zzz_xxyz[k] * tke_0;

            to_0_z_zzz_xxz[k] = -to_zzz_xx[k] + 2.0 * to_zzz_xxzz[k] * tke_0;

            to_0_z_zzz_xyy[k] = 2.0 * to_zzz_xyyz[k] * tke_0;

            to_0_z_zzz_xyz[k] = -to_zzz_xy[k] + 2.0 * to_zzz_xyzz[k] * tke_0;

            to_0_z_zzz_xzz[k] = -2.0 * to_zzz_xz[k] + 2.0 * to_zzz_xzzz[k] * tke_0;

            to_0_z_zzz_yyy[k] = 2.0 * to_zzz_yyyz[k] * tke_0;

            to_0_z_zzz_yyz[k] = -to_zzz_yy[k] + 2.0 * to_zzz_yyzz[k] * tke_0;

            to_0_z_zzz_yzz[k] = -2.0 * to_zzz_yz[k] + 2.0 * to_zzz_yzzz[k] * tke_0;

            to_0_z_zzz_zzz[k] = -3.0 * to_zzz_zz[k] + 2.0 * to_zzz_zzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

