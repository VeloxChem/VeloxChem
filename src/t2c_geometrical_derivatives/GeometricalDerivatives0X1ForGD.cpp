#include "GeometricalDerivatives0X1ForGD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_gd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gd,
                       const int idx_op_gp,
                       const int idx_op_gf,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : GP

        auto to_xxxx_x = pbuffer.data(idx_op_gp + i * 45 + 0);

        auto to_xxxx_y = pbuffer.data(idx_op_gp + i * 45 + 1);

        auto to_xxxx_z = pbuffer.data(idx_op_gp + i * 45 + 2);

        auto to_xxxy_x = pbuffer.data(idx_op_gp + i * 45 + 3);

        auto to_xxxy_y = pbuffer.data(idx_op_gp + i * 45 + 4);

        auto to_xxxy_z = pbuffer.data(idx_op_gp + i * 45 + 5);

        auto to_xxxz_x = pbuffer.data(idx_op_gp + i * 45 + 6);

        auto to_xxxz_y = pbuffer.data(idx_op_gp + i * 45 + 7);

        auto to_xxxz_z = pbuffer.data(idx_op_gp + i * 45 + 8);

        auto to_xxyy_x = pbuffer.data(idx_op_gp + i * 45 + 9);

        auto to_xxyy_y = pbuffer.data(idx_op_gp + i * 45 + 10);

        auto to_xxyy_z = pbuffer.data(idx_op_gp + i * 45 + 11);

        auto to_xxyz_x = pbuffer.data(idx_op_gp + i * 45 + 12);

        auto to_xxyz_y = pbuffer.data(idx_op_gp + i * 45 + 13);

        auto to_xxyz_z = pbuffer.data(idx_op_gp + i * 45 + 14);

        auto to_xxzz_x = pbuffer.data(idx_op_gp + i * 45 + 15);

        auto to_xxzz_y = pbuffer.data(idx_op_gp + i * 45 + 16);

        auto to_xxzz_z = pbuffer.data(idx_op_gp + i * 45 + 17);

        auto to_xyyy_x = pbuffer.data(idx_op_gp + i * 45 + 18);

        auto to_xyyy_y = pbuffer.data(idx_op_gp + i * 45 + 19);

        auto to_xyyy_z = pbuffer.data(idx_op_gp + i * 45 + 20);

        auto to_xyyz_x = pbuffer.data(idx_op_gp + i * 45 + 21);

        auto to_xyyz_y = pbuffer.data(idx_op_gp + i * 45 + 22);

        auto to_xyyz_z = pbuffer.data(idx_op_gp + i * 45 + 23);

        auto to_xyzz_x = pbuffer.data(idx_op_gp + i * 45 + 24);

        auto to_xyzz_y = pbuffer.data(idx_op_gp + i * 45 + 25);

        auto to_xyzz_z = pbuffer.data(idx_op_gp + i * 45 + 26);

        auto to_xzzz_x = pbuffer.data(idx_op_gp + i * 45 + 27);

        auto to_xzzz_y = pbuffer.data(idx_op_gp + i * 45 + 28);

        auto to_xzzz_z = pbuffer.data(idx_op_gp + i * 45 + 29);

        auto to_yyyy_x = pbuffer.data(idx_op_gp + i * 45 + 30);

        auto to_yyyy_y = pbuffer.data(idx_op_gp + i * 45 + 31);

        auto to_yyyy_z = pbuffer.data(idx_op_gp + i * 45 + 32);

        auto to_yyyz_x = pbuffer.data(idx_op_gp + i * 45 + 33);

        auto to_yyyz_y = pbuffer.data(idx_op_gp + i * 45 + 34);

        auto to_yyyz_z = pbuffer.data(idx_op_gp + i * 45 + 35);

        auto to_yyzz_x = pbuffer.data(idx_op_gp + i * 45 + 36);

        auto to_yyzz_y = pbuffer.data(idx_op_gp + i * 45 + 37);

        auto to_yyzz_z = pbuffer.data(idx_op_gp + i * 45 + 38);

        auto to_yzzz_x = pbuffer.data(idx_op_gp + i * 45 + 39);

        auto to_yzzz_y = pbuffer.data(idx_op_gp + i * 45 + 40);

        auto to_yzzz_z = pbuffer.data(idx_op_gp + i * 45 + 41);

        auto to_zzzz_x = pbuffer.data(idx_op_gp + i * 45 + 42);

        auto to_zzzz_y = pbuffer.data(idx_op_gp + i * 45 + 43);

        auto to_zzzz_z = pbuffer.data(idx_op_gp + i * 45 + 44);

        // Set up components of auxiliary buffer : GF

        auto to_xxxx_xxx = pbuffer.data(idx_op_gf + i * 150 + 0);

        auto to_xxxx_xxy = pbuffer.data(idx_op_gf + i * 150 + 1);

        auto to_xxxx_xxz = pbuffer.data(idx_op_gf + i * 150 + 2);

        auto to_xxxx_xyy = pbuffer.data(idx_op_gf + i * 150 + 3);

        auto to_xxxx_xyz = pbuffer.data(idx_op_gf + i * 150 + 4);

        auto to_xxxx_xzz = pbuffer.data(idx_op_gf + i * 150 + 5);

        auto to_xxxx_yyy = pbuffer.data(idx_op_gf + i * 150 + 6);

        auto to_xxxx_yyz = pbuffer.data(idx_op_gf + i * 150 + 7);

        auto to_xxxx_yzz = pbuffer.data(idx_op_gf + i * 150 + 8);

        auto to_xxxx_zzz = pbuffer.data(idx_op_gf + i * 150 + 9);

        auto to_xxxy_xxx = pbuffer.data(idx_op_gf + i * 150 + 10);

        auto to_xxxy_xxy = pbuffer.data(idx_op_gf + i * 150 + 11);

        auto to_xxxy_xxz = pbuffer.data(idx_op_gf + i * 150 + 12);

        auto to_xxxy_xyy = pbuffer.data(idx_op_gf + i * 150 + 13);

        auto to_xxxy_xyz = pbuffer.data(idx_op_gf + i * 150 + 14);

        auto to_xxxy_xzz = pbuffer.data(idx_op_gf + i * 150 + 15);

        auto to_xxxy_yyy = pbuffer.data(idx_op_gf + i * 150 + 16);

        auto to_xxxy_yyz = pbuffer.data(idx_op_gf + i * 150 + 17);

        auto to_xxxy_yzz = pbuffer.data(idx_op_gf + i * 150 + 18);

        auto to_xxxy_zzz = pbuffer.data(idx_op_gf + i * 150 + 19);

        auto to_xxxz_xxx = pbuffer.data(idx_op_gf + i * 150 + 20);

        auto to_xxxz_xxy = pbuffer.data(idx_op_gf + i * 150 + 21);

        auto to_xxxz_xxz = pbuffer.data(idx_op_gf + i * 150 + 22);

        auto to_xxxz_xyy = pbuffer.data(idx_op_gf + i * 150 + 23);

        auto to_xxxz_xyz = pbuffer.data(idx_op_gf + i * 150 + 24);

        auto to_xxxz_xzz = pbuffer.data(idx_op_gf + i * 150 + 25);

        auto to_xxxz_yyy = pbuffer.data(idx_op_gf + i * 150 + 26);

        auto to_xxxz_yyz = pbuffer.data(idx_op_gf + i * 150 + 27);

        auto to_xxxz_yzz = pbuffer.data(idx_op_gf + i * 150 + 28);

        auto to_xxxz_zzz = pbuffer.data(idx_op_gf + i * 150 + 29);

        auto to_xxyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 30);

        auto to_xxyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 31);

        auto to_xxyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 32);

        auto to_xxyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 33);

        auto to_xxyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 34);

        auto to_xxyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 35);

        auto to_xxyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 36);

        auto to_xxyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 37);

        auto to_xxyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 38);

        auto to_xxyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 39);

        auto to_xxyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 40);

        auto to_xxyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 41);

        auto to_xxyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 42);

        auto to_xxyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 43);

        auto to_xxyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 44);

        auto to_xxyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 45);

        auto to_xxyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 46);

        auto to_xxyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 47);

        auto to_xxyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 48);

        auto to_xxyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 49);

        auto to_xxzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 50);

        auto to_xxzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 51);

        auto to_xxzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 52);

        auto to_xxzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 53);

        auto to_xxzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 54);

        auto to_xxzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 55);

        auto to_xxzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 56);

        auto to_xxzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 57);

        auto to_xxzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 58);

        auto to_xxzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 59);

        auto to_xyyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 60);

        auto to_xyyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 61);

        auto to_xyyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 62);

        auto to_xyyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 63);

        auto to_xyyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 64);

        auto to_xyyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 65);

        auto to_xyyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 66);

        auto to_xyyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 67);

        auto to_xyyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 68);

        auto to_xyyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 69);

        auto to_xyyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 70);

        auto to_xyyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 71);

        auto to_xyyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 72);

        auto to_xyyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 73);

        auto to_xyyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 74);

        auto to_xyyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 75);

        auto to_xyyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 76);

        auto to_xyyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 77);

        auto to_xyyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 78);

        auto to_xyyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 79);

        auto to_xyzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 80);

        auto to_xyzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 81);

        auto to_xyzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 82);

        auto to_xyzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 83);

        auto to_xyzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 84);

        auto to_xyzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 85);

        auto to_xyzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 86);

        auto to_xyzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 87);

        auto to_xyzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 88);

        auto to_xyzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 89);

        auto to_xzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 90);

        auto to_xzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 91);

        auto to_xzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 92);

        auto to_xzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 93);

        auto to_xzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 94);

        auto to_xzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 95);

        auto to_xzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 96);

        auto to_xzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 97);

        auto to_xzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 98);

        auto to_xzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 99);

        auto to_yyyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 100);

        auto to_yyyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 101);

        auto to_yyyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 102);

        auto to_yyyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 103);

        auto to_yyyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 104);

        auto to_yyyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 105);

        auto to_yyyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 106);

        auto to_yyyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 107);

        auto to_yyyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 108);

        auto to_yyyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 109);

        auto to_yyyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 110);

        auto to_yyyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 111);

        auto to_yyyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 112);

        auto to_yyyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 113);

        auto to_yyyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 114);

        auto to_yyyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 115);

        auto to_yyyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 116);

        auto to_yyyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 117);

        auto to_yyyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 118);

        auto to_yyyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 119);

        auto to_yyzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 120);

        auto to_yyzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 121);

        auto to_yyzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 122);

        auto to_yyzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 123);

        auto to_yyzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 124);

        auto to_yyzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 125);

        auto to_yyzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 126);

        auto to_yyzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 127);

        auto to_yyzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 128);

        auto to_yyzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 129);

        auto to_yzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 130);

        auto to_yzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 131);

        auto to_yzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 132);

        auto to_yzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 133);

        auto to_yzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 134);

        auto to_yzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 135);

        auto to_yzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 136);

        auto to_yzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 137);

        auto to_yzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 138);

        auto to_yzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 139);

        auto to_zzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 140);

        auto to_zzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 141);

        auto to_zzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 142);

        auto to_zzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 143);

        auto to_zzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 144);

        auto to_zzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 145);

        auto to_zzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 146);

        auto to_zzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 147);

        auto to_zzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 148);

        auto to_zzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 149);

        // Set up 0-6 components of targeted buffer : GD

        auto to_0_x_xxxx_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 0);

        auto to_0_x_xxxx_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 1);

        auto to_0_x_xxxx_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 2);

        auto to_0_x_xxxx_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 3);

        auto to_0_x_xxxx_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 4);

        auto to_0_x_xxxx_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_0_x_xxxx_xx, to_0_x_xxxx_xy, to_0_x_xxxx_xz, to_0_x_xxxx_yy, to_0_x_xxxx_yz, to_0_x_xxxx_zz, to_xxxx_x, to_xxxx_xxx, to_xxxx_xxy, to_xxxx_xxz, to_xxxx_xyy, to_xxxx_xyz, to_xxxx_xzz, to_xxxx_y, to_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxx_xx[k] = -2.0 * to_xxxx_x[k] + 2.0 * to_xxxx_xxx[k] * tke_0;

            to_0_x_xxxx_xy[k] = -to_xxxx_y[k] + 2.0 * to_xxxx_xxy[k] * tke_0;

            to_0_x_xxxx_xz[k] = -to_xxxx_z[k] + 2.0 * to_xxxx_xxz[k] * tke_0;

            to_0_x_xxxx_yy[k] = 2.0 * to_xxxx_xyy[k] * tke_0;

            to_0_x_xxxx_yz[k] = 2.0 * to_xxxx_xyz[k] * tke_0;

            to_0_x_xxxx_zz[k] = 2.0 * to_xxxx_xzz[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : GD

        auto to_0_x_xxxy_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 6);

        auto to_0_x_xxxy_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 7);

        auto to_0_x_xxxy_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 8);

        auto to_0_x_xxxy_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 9);

        auto to_0_x_xxxy_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 10);

        auto to_0_x_xxxy_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_0_x_xxxy_xx, to_0_x_xxxy_xy, to_0_x_xxxy_xz, to_0_x_xxxy_yy, to_0_x_xxxy_yz, to_0_x_xxxy_zz, to_xxxy_x, to_xxxy_xxx, to_xxxy_xxy, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxy_xx[k] = -2.0 * to_xxxy_x[k] + 2.0 * to_xxxy_xxx[k] * tke_0;

            to_0_x_xxxy_xy[k] = -to_xxxy_y[k] + 2.0 * to_xxxy_xxy[k] * tke_0;

            to_0_x_xxxy_xz[k] = -to_xxxy_z[k] + 2.0 * to_xxxy_xxz[k] * tke_0;

            to_0_x_xxxy_yy[k] = 2.0 * to_xxxy_xyy[k] * tke_0;

            to_0_x_xxxy_yz[k] = 2.0 * to_xxxy_xyz[k] * tke_0;

            to_0_x_xxxy_zz[k] = 2.0 * to_xxxy_xzz[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : GD

        auto to_0_x_xxxz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 12);

        auto to_0_x_xxxz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 13);

        auto to_0_x_xxxz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 14);

        auto to_0_x_xxxz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 15);

        auto to_0_x_xxxz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 16);

        auto to_0_x_xxxz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_0_x_xxxz_xx, to_0_x_xxxz_xy, to_0_x_xxxz_xz, to_0_x_xxxz_yy, to_0_x_xxxz_yz, to_0_x_xxxz_zz, to_xxxz_x, to_xxxz_xxx, to_xxxz_xxy, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxz_xx[k] = -2.0 * to_xxxz_x[k] + 2.0 * to_xxxz_xxx[k] * tke_0;

            to_0_x_xxxz_xy[k] = -to_xxxz_y[k] + 2.0 * to_xxxz_xxy[k] * tke_0;

            to_0_x_xxxz_xz[k] = -to_xxxz_z[k] + 2.0 * to_xxxz_xxz[k] * tke_0;

            to_0_x_xxxz_yy[k] = 2.0 * to_xxxz_xyy[k] * tke_0;

            to_0_x_xxxz_yz[k] = 2.0 * to_xxxz_xyz[k] * tke_0;

            to_0_x_xxxz_zz[k] = 2.0 * to_xxxz_xzz[k] * tke_0;
        }

        // Set up 18-24 components of targeted buffer : GD

        auto to_0_x_xxyy_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 18);

        auto to_0_x_xxyy_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 19);

        auto to_0_x_xxyy_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 20);

        auto to_0_x_xxyy_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 21);

        auto to_0_x_xxyy_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 22);

        auto to_0_x_xxyy_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_0_x_xxyy_xx, to_0_x_xxyy_xy, to_0_x_xxyy_xz, to_0_x_xxyy_yy, to_0_x_xxyy_yz, to_0_x_xxyy_zz, to_xxyy_x, to_xxyy_xxx, to_xxyy_xxy, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyy_xx[k] = -2.0 * to_xxyy_x[k] + 2.0 * to_xxyy_xxx[k] * tke_0;

            to_0_x_xxyy_xy[k] = -to_xxyy_y[k] + 2.0 * to_xxyy_xxy[k] * tke_0;

            to_0_x_xxyy_xz[k] = -to_xxyy_z[k] + 2.0 * to_xxyy_xxz[k] * tke_0;

            to_0_x_xxyy_yy[k] = 2.0 * to_xxyy_xyy[k] * tke_0;

            to_0_x_xxyy_yz[k] = 2.0 * to_xxyy_xyz[k] * tke_0;

            to_0_x_xxyy_zz[k] = 2.0 * to_xxyy_xzz[k] * tke_0;
        }

        // Set up 24-30 components of targeted buffer : GD

        auto to_0_x_xxyz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 24);

        auto to_0_x_xxyz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 25);

        auto to_0_x_xxyz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 26);

        auto to_0_x_xxyz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 27);

        auto to_0_x_xxyz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 28);

        auto to_0_x_xxyz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_x_xxyz_xx, to_0_x_xxyz_xy, to_0_x_xxyz_xz, to_0_x_xxyz_yy, to_0_x_xxyz_yz, to_0_x_xxyz_zz, to_xxyz_x, to_xxyz_xxx, to_xxyz_xxy, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyz_xx[k] = -2.0 * to_xxyz_x[k] + 2.0 * to_xxyz_xxx[k] * tke_0;

            to_0_x_xxyz_xy[k] = -to_xxyz_y[k] + 2.0 * to_xxyz_xxy[k] * tke_0;

            to_0_x_xxyz_xz[k] = -to_xxyz_z[k] + 2.0 * to_xxyz_xxz[k] * tke_0;

            to_0_x_xxyz_yy[k] = 2.0 * to_xxyz_xyy[k] * tke_0;

            to_0_x_xxyz_yz[k] = 2.0 * to_xxyz_xyz[k] * tke_0;

            to_0_x_xxyz_zz[k] = 2.0 * to_xxyz_xzz[k] * tke_0;
        }

        // Set up 30-36 components of targeted buffer : GD

        auto to_0_x_xxzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 30);

        auto to_0_x_xxzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 31);

        auto to_0_x_xxzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 32);

        auto to_0_x_xxzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 33);

        auto to_0_x_xxzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 34);

        auto to_0_x_xxzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_0_x_xxzz_xx, to_0_x_xxzz_xy, to_0_x_xxzz_xz, to_0_x_xxzz_yy, to_0_x_xxzz_yz, to_0_x_xxzz_zz, to_xxzz_x, to_xxzz_xxx, to_xxzz_xxy, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxzz_xx[k] = -2.0 * to_xxzz_x[k] + 2.0 * to_xxzz_xxx[k] * tke_0;

            to_0_x_xxzz_xy[k] = -to_xxzz_y[k] + 2.0 * to_xxzz_xxy[k] * tke_0;

            to_0_x_xxzz_xz[k] = -to_xxzz_z[k] + 2.0 * to_xxzz_xxz[k] * tke_0;

            to_0_x_xxzz_yy[k] = 2.0 * to_xxzz_xyy[k] * tke_0;

            to_0_x_xxzz_yz[k] = 2.0 * to_xxzz_xyz[k] * tke_0;

            to_0_x_xxzz_zz[k] = 2.0 * to_xxzz_xzz[k] * tke_0;
        }

        // Set up 36-42 components of targeted buffer : GD

        auto to_0_x_xyyy_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 36);

        auto to_0_x_xyyy_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 37);

        auto to_0_x_xyyy_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 38);

        auto to_0_x_xyyy_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 39);

        auto to_0_x_xyyy_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 40);

        auto to_0_x_xyyy_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_0_x_xyyy_xx, to_0_x_xyyy_xy, to_0_x_xyyy_xz, to_0_x_xyyy_yy, to_0_x_xyyy_yz, to_0_x_xyyy_zz, to_xyyy_x, to_xyyy_xxx, to_xyyy_xxy, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyy_xx[k] = -2.0 * to_xyyy_x[k] + 2.0 * to_xyyy_xxx[k] * tke_0;

            to_0_x_xyyy_xy[k] = -to_xyyy_y[k] + 2.0 * to_xyyy_xxy[k] * tke_0;

            to_0_x_xyyy_xz[k] = -to_xyyy_z[k] + 2.0 * to_xyyy_xxz[k] * tke_0;

            to_0_x_xyyy_yy[k] = 2.0 * to_xyyy_xyy[k] * tke_0;

            to_0_x_xyyy_yz[k] = 2.0 * to_xyyy_xyz[k] * tke_0;

            to_0_x_xyyy_zz[k] = 2.0 * to_xyyy_xzz[k] * tke_0;
        }

        // Set up 42-48 components of targeted buffer : GD

        auto to_0_x_xyyz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 42);

        auto to_0_x_xyyz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 43);

        auto to_0_x_xyyz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 44);

        auto to_0_x_xyyz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 45);

        auto to_0_x_xyyz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 46);

        auto to_0_x_xyyz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_0_x_xyyz_xx, to_0_x_xyyz_xy, to_0_x_xyyz_xz, to_0_x_xyyz_yy, to_0_x_xyyz_yz, to_0_x_xyyz_zz, to_xyyz_x, to_xyyz_xxx, to_xyyz_xxy, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyz_xx[k] = -2.0 * to_xyyz_x[k] + 2.0 * to_xyyz_xxx[k] * tke_0;

            to_0_x_xyyz_xy[k] = -to_xyyz_y[k] + 2.0 * to_xyyz_xxy[k] * tke_0;

            to_0_x_xyyz_xz[k] = -to_xyyz_z[k] + 2.0 * to_xyyz_xxz[k] * tke_0;

            to_0_x_xyyz_yy[k] = 2.0 * to_xyyz_xyy[k] * tke_0;

            to_0_x_xyyz_yz[k] = 2.0 * to_xyyz_xyz[k] * tke_0;

            to_0_x_xyyz_zz[k] = 2.0 * to_xyyz_xzz[k] * tke_0;
        }

        // Set up 48-54 components of targeted buffer : GD

        auto to_0_x_xyzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 48);

        auto to_0_x_xyzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 49);

        auto to_0_x_xyzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 50);

        auto to_0_x_xyzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 51);

        auto to_0_x_xyzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 52);

        auto to_0_x_xyzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_0_x_xyzz_xx, to_0_x_xyzz_xy, to_0_x_xyzz_xz, to_0_x_xyzz_yy, to_0_x_xyzz_yz, to_0_x_xyzz_zz, to_xyzz_x, to_xyzz_xxx, to_xyzz_xxy, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyzz_xx[k] = -2.0 * to_xyzz_x[k] + 2.0 * to_xyzz_xxx[k] * tke_0;

            to_0_x_xyzz_xy[k] = -to_xyzz_y[k] + 2.0 * to_xyzz_xxy[k] * tke_0;

            to_0_x_xyzz_xz[k] = -to_xyzz_z[k] + 2.0 * to_xyzz_xxz[k] * tke_0;

            to_0_x_xyzz_yy[k] = 2.0 * to_xyzz_xyy[k] * tke_0;

            to_0_x_xyzz_yz[k] = 2.0 * to_xyzz_xyz[k] * tke_0;

            to_0_x_xyzz_zz[k] = 2.0 * to_xyzz_xzz[k] * tke_0;
        }

        // Set up 54-60 components of targeted buffer : GD

        auto to_0_x_xzzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 54);

        auto to_0_x_xzzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 55);

        auto to_0_x_xzzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 56);

        auto to_0_x_xzzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 57);

        auto to_0_x_xzzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 58);

        auto to_0_x_xzzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_x_xzzz_xx, to_0_x_xzzz_xy, to_0_x_xzzz_xz, to_0_x_xzzz_yy, to_0_x_xzzz_yz, to_0_x_xzzz_zz, to_xzzz_x, to_xzzz_xxx, to_xzzz_xxy, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzzz_xx[k] = -2.0 * to_xzzz_x[k] + 2.0 * to_xzzz_xxx[k] * tke_0;

            to_0_x_xzzz_xy[k] = -to_xzzz_y[k] + 2.0 * to_xzzz_xxy[k] * tke_0;

            to_0_x_xzzz_xz[k] = -to_xzzz_z[k] + 2.0 * to_xzzz_xxz[k] * tke_0;

            to_0_x_xzzz_yy[k] = 2.0 * to_xzzz_xyy[k] * tke_0;

            to_0_x_xzzz_yz[k] = 2.0 * to_xzzz_xyz[k] * tke_0;

            to_0_x_xzzz_zz[k] = 2.0 * to_xzzz_xzz[k] * tke_0;
        }

        // Set up 60-66 components of targeted buffer : GD

        auto to_0_x_yyyy_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 60);

        auto to_0_x_yyyy_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 61);

        auto to_0_x_yyyy_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 62);

        auto to_0_x_yyyy_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 63);

        auto to_0_x_yyyy_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 64);

        auto to_0_x_yyyy_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_0_x_yyyy_xx, to_0_x_yyyy_xy, to_0_x_yyyy_xz, to_0_x_yyyy_yy, to_0_x_yyyy_yz, to_0_x_yyyy_zz, to_yyyy_x, to_yyyy_xxx, to_yyyy_xxy, to_yyyy_xxz, to_yyyy_xyy, to_yyyy_xyz, to_yyyy_xzz, to_yyyy_y, to_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyy_xx[k] = -2.0 * to_yyyy_x[k] + 2.0 * to_yyyy_xxx[k] * tke_0;

            to_0_x_yyyy_xy[k] = -to_yyyy_y[k] + 2.0 * to_yyyy_xxy[k] * tke_0;

            to_0_x_yyyy_xz[k] = -to_yyyy_z[k] + 2.0 * to_yyyy_xxz[k] * tke_0;

            to_0_x_yyyy_yy[k] = 2.0 * to_yyyy_xyy[k] * tke_0;

            to_0_x_yyyy_yz[k] = 2.0 * to_yyyy_xyz[k] * tke_0;

            to_0_x_yyyy_zz[k] = 2.0 * to_yyyy_xzz[k] * tke_0;
        }

        // Set up 66-72 components of targeted buffer : GD

        auto to_0_x_yyyz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 66);

        auto to_0_x_yyyz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 67);

        auto to_0_x_yyyz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 68);

        auto to_0_x_yyyz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 69);

        auto to_0_x_yyyz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 70);

        auto to_0_x_yyyz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_0_x_yyyz_xx, to_0_x_yyyz_xy, to_0_x_yyyz_xz, to_0_x_yyyz_yy, to_0_x_yyyz_yz, to_0_x_yyyz_zz, to_yyyz_x, to_yyyz_xxx, to_yyyz_xxy, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyz_xx[k] = -2.0 * to_yyyz_x[k] + 2.0 * to_yyyz_xxx[k] * tke_0;

            to_0_x_yyyz_xy[k] = -to_yyyz_y[k] + 2.0 * to_yyyz_xxy[k] * tke_0;

            to_0_x_yyyz_xz[k] = -to_yyyz_z[k] + 2.0 * to_yyyz_xxz[k] * tke_0;

            to_0_x_yyyz_yy[k] = 2.0 * to_yyyz_xyy[k] * tke_0;

            to_0_x_yyyz_yz[k] = 2.0 * to_yyyz_xyz[k] * tke_0;

            to_0_x_yyyz_zz[k] = 2.0 * to_yyyz_xzz[k] * tke_0;
        }

        // Set up 72-78 components of targeted buffer : GD

        auto to_0_x_yyzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 72);

        auto to_0_x_yyzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 73);

        auto to_0_x_yyzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 74);

        auto to_0_x_yyzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 75);

        auto to_0_x_yyzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 76);

        auto to_0_x_yyzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_0_x_yyzz_xx, to_0_x_yyzz_xy, to_0_x_yyzz_xz, to_0_x_yyzz_yy, to_0_x_yyzz_yz, to_0_x_yyzz_zz, to_yyzz_x, to_yyzz_xxx, to_yyzz_xxy, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyzz_xx[k] = -2.0 * to_yyzz_x[k] + 2.0 * to_yyzz_xxx[k] * tke_0;

            to_0_x_yyzz_xy[k] = -to_yyzz_y[k] + 2.0 * to_yyzz_xxy[k] * tke_0;

            to_0_x_yyzz_xz[k] = -to_yyzz_z[k] + 2.0 * to_yyzz_xxz[k] * tke_0;

            to_0_x_yyzz_yy[k] = 2.0 * to_yyzz_xyy[k] * tke_0;

            to_0_x_yyzz_yz[k] = 2.0 * to_yyzz_xyz[k] * tke_0;

            to_0_x_yyzz_zz[k] = 2.0 * to_yyzz_xzz[k] * tke_0;
        }

        // Set up 78-84 components of targeted buffer : GD

        auto to_0_x_yzzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 78);

        auto to_0_x_yzzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 79);

        auto to_0_x_yzzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 80);

        auto to_0_x_yzzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 81);

        auto to_0_x_yzzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 82);

        auto to_0_x_yzzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_0_x_yzzz_xx, to_0_x_yzzz_xy, to_0_x_yzzz_xz, to_0_x_yzzz_yy, to_0_x_yzzz_yz, to_0_x_yzzz_zz, to_yzzz_x, to_yzzz_xxx, to_yzzz_xxy, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzzz_xx[k] = -2.0 * to_yzzz_x[k] + 2.0 * to_yzzz_xxx[k] * tke_0;

            to_0_x_yzzz_xy[k] = -to_yzzz_y[k] + 2.0 * to_yzzz_xxy[k] * tke_0;

            to_0_x_yzzz_xz[k] = -to_yzzz_z[k] + 2.0 * to_yzzz_xxz[k] * tke_0;

            to_0_x_yzzz_yy[k] = 2.0 * to_yzzz_xyy[k] * tke_0;

            to_0_x_yzzz_yz[k] = 2.0 * to_yzzz_xyz[k] * tke_0;

            to_0_x_yzzz_zz[k] = 2.0 * to_yzzz_xzz[k] * tke_0;
        }

        // Set up 84-90 components of targeted buffer : GD

        auto to_0_x_zzzz_xx = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 84);

        auto to_0_x_zzzz_xy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 85);

        auto to_0_x_zzzz_xz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 86);

        auto to_0_x_zzzz_yy = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 87);

        auto to_0_x_zzzz_yz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 88);

        auto to_0_x_zzzz_zz = pbuffer.data(idx_op_geom_001_gd + 0 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_x_zzzz_xx, to_0_x_zzzz_xy, to_0_x_zzzz_xz, to_0_x_zzzz_yy, to_0_x_zzzz_yz, to_0_x_zzzz_zz, to_zzzz_x, to_zzzz_xxx, to_zzzz_xxy, to_zzzz_xxz, to_zzzz_xyy, to_zzzz_xyz, to_zzzz_xzz, to_zzzz_y, to_zzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzzz_xx[k] = -2.0 * to_zzzz_x[k] + 2.0 * to_zzzz_xxx[k] * tke_0;

            to_0_x_zzzz_xy[k] = -to_zzzz_y[k] + 2.0 * to_zzzz_xxy[k] * tke_0;

            to_0_x_zzzz_xz[k] = -to_zzzz_z[k] + 2.0 * to_zzzz_xxz[k] * tke_0;

            to_0_x_zzzz_yy[k] = 2.0 * to_zzzz_xyy[k] * tke_0;

            to_0_x_zzzz_yz[k] = 2.0 * to_zzzz_xyz[k] * tke_0;

            to_0_x_zzzz_zz[k] = 2.0 * to_zzzz_xzz[k] * tke_0;
        }

        // Set up 90-96 components of targeted buffer : GD

        auto to_0_y_xxxx_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 0);

        auto to_0_y_xxxx_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 1);

        auto to_0_y_xxxx_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 2);

        auto to_0_y_xxxx_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 3);

        auto to_0_y_xxxx_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 4);

        auto to_0_y_xxxx_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_0_y_xxxx_xx, to_0_y_xxxx_xy, to_0_y_xxxx_xz, to_0_y_xxxx_yy, to_0_y_xxxx_yz, to_0_y_xxxx_zz, to_xxxx_x, to_xxxx_xxy, to_xxxx_xyy, to_xxxx_xyz, to_xxxx_y, to_xxxx_yyy, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxx_xx[k] = 2.0 * to_xxxx_xxy[k] * tke_0;

            to_0_y_xxxx_xy[k] = -to_xxxx_x[k] + 2.0 * to_xxxx_xyy[k] * tke_0;

            to_0_y_xxxx_xz[k] = 2.0 * to_xxxx_xyz[k] * tke_0;

            to_0_y_xxxx_yy[k] = -2.0 * to_xxxx_y[k] + 2.0 * to_xxxx_yyy[k] * tke_0;

            to_0_y_xxxx_yz[k] = -to_xxxx_z[k] + 2.0 * to_xxxx_yyz[k] * tke_0;

            to_0_y_xxxx_zz[k] = 2.0 * to_xxxx_yzz[k] * tke_0;
        }

        // Set up 96-102 components of targeted buffer : GD

        auto to_0_y_xxxy_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 6);

        auto to_0_y_xxxy_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 7);

        auto to_0_y_xxxy_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 8);

        auto to_0_y_xxxy_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 9);

        auto to_0_y_xxxy_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 10);

        auto to_0_y_xxxy_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_0_y_xxxy_xx, to_0_y_xxxy_xy, to_0_y_xxxy_xz, to_0_y_xxxy_yy, to_0_y_xxxy_yz, to_0_y_xxxy_zz, to_xxxy_x, to_xxxy_xxy, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_y, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxy_xx[k] = 2.0 * to_xxxy_xxy[k] * tke_0;

            to_0_y_xxxy_xy[k] = -to_xxxy_x[k] + 2.0 * to_xxxy_xyy[k] * tke_0;

            to_0_y_xxxy_xz[k] = 2.0 * to_xxxy_xyz[k] * tke_0;

            to_0_y_xxxy_yy[k] = -2.0 * to_xxxy_y[k] + 2.0 * to_xxxy_yyy[k] * tke_0;

            to_0_y_xxxy_yz[k] = -to_xxxy_z[k] + 2.0 * to_xxxy_yyz[k] * tke_0;

            to_0_y_xxxy_zz[k] = 2.0 * to_xxxy_yzz[k] * tke_0;
        }

        // Set up 102-108 components of targeted buffer : GD

        auto to_0_y_xxxz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 12);

        auto to_0_y_xxxz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 13);

        auto to_0_y_xxxz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 14);

        auto to_0_y_xxxz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 15);

        auto to_0_y_xxxz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 16);

        auto to_0_y_xxxz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_0_y_xxxz_xx, to_0_y_xxxz_xy, to_0_y_xxxz_xz, to_0_y_xxxz_yy, to_0_y_xxxz_yz, to_0_y_xxxz_zz, to_xxxz_x, to_xxxz_xxy, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_y, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxz_xx[k] = 2.0 * to_xxxz_xxy[k] * tke_0;

            to_0_y_xxxz_xy[k] = -to_xxxz_x[k] + 2.0 * to_xxxz_xyy[k] * tke_0;

            to_0_y_xxxz_xz[k] = 2.0 * to_xxxz_xyz[k] * tke_0;

            to_0_y_xxxz_yy[k] = -2.0 * to_xxxz_y[k] + 2.0 * to_xxxz_yyy[k] * tke_0;

            to_0_y_xxxz_yz[k] = -to_xxxz_z[k] + 2.0 * to_xxxz_yyz[k] * tke_0;

            to_0_y_xxxz_zz[k] = 2.0 * to_xxxz_yzz[k] * tke_0;
        }

        // Set up 108-114 components of targeted buffer : GD

        auto to_0_y_xxyy_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 18);

        auto to_0_y_xxyy_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 19);

        auto to_0_y_xxyy_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 20);

        auto to_0_y_xxyy_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 21);

        auto to_0_y_xxyy_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 22);

        auto to_0_y_xxyy_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_0_y_xxyy_xx, to_0_y_xxyy_xy, to_0_y_xxyy_xz, to_0_y_xxyy_yy, to_0_y_xxyy_yz, to_0_y_xxyy_zz, to_xxyy_x, to_xxyy_xxy, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_y, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyy_xx[k] = 2.0 * to_xxyy_xxy[k] * tke_0;

            to_0_y_xxyy_xy[k] = -to_xxyy_x[k] + 2.0 * to_xxyy_xyy[k] * tke_0;

            to_0_y_xxyy_xz[k] = 2.0 * to_xxyy_xyz[k] * tke_0;

            to_0_y_xxyy_yy[k] = -2.0 * to_xxyy_y[k] + 2.0 * to_xxyy_yyy[k] * tke_0;

            to_0_y_xxyy_yz[k] = -to_xxyy_z[k] + 2.0 * to_xxyy_yyz[k] * tke_0;

            to_0_y_xxyy_zz[k] = 2.0 * to_xxyy_yzz[k] * tke_0;
        }

        // Set up 114-120 components of targeted buffer : GD

        auto to_0_y_xxyz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 24);

        auto to_0_y_xxyz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 25);

        auto to_0_y_xxyz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 26);

        auto to_0_y_xxyz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 27);

        auto to_0_y_xxyz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 28);

        auto to_0_y_xxyz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_y_xxyz_xx, to_0_y_xxyz_xy, to_0_y_xxyz_xz, to_0_y_xxyz_yy, to_0_y_xxyz_yz, to_0_y_xxyz_zz, to_xxyz_x, to_xxyz_xxy, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_y, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyz_xx[k] = 2.0 * to_xxyz_xxy[k] * tke_0;

            to_0_y_xxyz_xy[k] = -to_xxyz_x[k] + 2.0 * to_xxyz_xyy[k] * tke_0;

            to_0_y_xxyz_xz[k] = 2.0 * to_xxyz_xyz[k] * tke_0;

            to_0_y_xxyz_yy[k] = -2.0 * to_xxyz_y[k] + 2.0 * to_xxyz_yyy[k] * tke_0;

            to_0_y_xxyz_yz[k] = -to_xxyz_z[k] + 2.0 * to_xxyz_yyz[k] * tke_0;

            to_0_y_xxyz_zz[k] = 2.0 * to_xxyz_yzz[k] * tke_0;
        }

        // Set up 120-126 components of targeted buffer : GD

        auto to_0_y_xxzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 30);

        auto to_0_y_xxzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 31);

        auto to_0_y_xxzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 32);

        auto to_0_y_xxzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 33);

        auto to_0_y_xxzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 34);

        auto to_0_y_xxzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_0_y_xxzz_xx, to_0_y_xxzz_xy, to_0_y_xxzz_xz, to_0_y_xxzz_yy, to_0_y_xxzz_yz, to_0_y_xxzz_zz, to_xxzz_x, to_xxzz_xxy, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_y, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxzz_xx[k] = 2.0 * to_xxzz_xxy[k] * tke_0;

            to_0_y_xxzz_xy[k] = -to_xxzz_x[k] + 2.0 * to_xxzz_xyy[k] * tke_0;

            to_0_y_xxzz_xz[k] = 2.0 * to_xxzz_xyz[k] * tke_0;

            to_0_y_xxzz_yy[k] = -2.0 * to_xxzz_y[k] + 2.0 * to_xxzz_yyy[k] * tke_0;

            to_0_y_xxzz_yz[k] = -to_xxzz_z[k] + 2.0 * to_xxzz_yyz[k] * tke_0;

            to_0_y_xxzz_zz[k] = 2.0 * to_xxzz_yzz[k] * tke_0;
        }

        // Set up 126-132 components of targeted buffer : GD

        auto to_0_y_xyyy_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 36);

        auto to_0_y_xyyy_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 37);

        auto to_0_y_xyyy_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 38);

        auto to_0_y_xyyy_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 39);

        auto to_0_y_xyyy_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 40);

        auto to_0_y_xyyy_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_0_y_xyyy_xx, to_0_y_xyyy_xy, to_0_y_xyyy_xz, to_0_y_xyyy_yy, to_0_y_xyyy_yz, to_0_y_xyyy_zz, to_xyyy_x, to_xyyy_xxy, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_y, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyy_xx[k] = 2.0 * to_xyyy_xxy[k] * tke_0;

            to_0_y_xyyy_xy[k] = -to_xyyy_x[k] + 2.0 * to_xyyy_xyy[k] * tke_0;

            to_0_y_xyyy_xz[k] = 2.0 * to_xyyy_xyz[k] * tke_0;

            to_0_y_xyyy_yy[k] = -2.0 * to_xyyy_y[k] + 2.0 * to_xyyy_yyy[k] * tke_0;

            to_0_y_xyyy_yz[k] = -to_xyyy_z[k] + 2.0 * to_xyyy_yyz[k] * tke_0;

            to_0_y_xyyy_zz[k] = 2.0 * to_xyyy_yzz[k] * tke_0;
        }

        // Set up 132-138 components of targeted buffer : GD

        auto to_0_y_xyyz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 42);

        auto to_0_y_xyyz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 43);

        auto to_0_y_xyyz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 44);

        auto to_0_y_xyyz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 45);

        auto to_0_y_xyyz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 46);

        auto to_0_y_xyyz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_0_y_xyyz_xx, to_0_y_xyyz_xy, to_0_y_xyyz_xz, to_0_y_xyyz_yy, to_0_y_xyyz_yz, to_0_y_xyyz_zz, to_xyyz_x, to_xyyz_xxy, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_y, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyz_xx[k] = 2.0 * to_xyyz_xxy[k] * tke_0;

            to_0_y_xyyz_xy[k] = -to_xyyz_x[k] + 2.0 * to_xyyz_xyy[k] * tke_0;

            to_0_y_xyyz_xz[k] = 2.0 * to_xyyz_xyz[k] * tke_0;

            to_0_y_xyyz_yy[k] = -2.0 * to_xyyz_y[k] + 2.0 * to_xyyz_yyy[k] * tke_0;

            to_0_y_xyyz_yz[k] = -to_xyyz_z[k] + 2.0 * to_xyyz_yyz[k] * tke_0;

            to_0_y_xyyz_zz[k] = 2.0 * to_xyyz_yzz[k] * tke_0;
        }

        // Set up 138-144 components of targeted buffer : GD

        auto to_0_y_xyzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 48);

        auto to_0_y_xyzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 49);

        auto to_0_y_xyzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 50);

        auto to_0_y_xyzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 51);

        auto to_0_y_xyzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 52);

        auto to_0_y_xyzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_0_y_xyzz_xx, to_0_y_xyzz_xy, to_0_y_xyzz_xz, to_0_y_xyzz_yy, to_0_y_xyzz_yz, to_0_y_xyzz_zz, to_xyzz_x, to_xyzz_xxy, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_y, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyzz_xx[k] = 2.0 * to_xyzz_xxy[k] * tke_0;

            to_0_y_xyzz_xy[k] = -to_xyzz_x[k] + 2.0 * to_xyzz_xyy[k] * tke_0;

            to_0_y_xyzz_xz[k] = 2.0 * to_xyzz_xyz[k] * tke_0;

            to_0_y_xyzz_yy[k] = -2.0 * to_xyzz_y[k] + 2.0 * to_xyzz_yyy[k] * tke_0;

            to_0_y_xyzz_yz[k] = -to_xyzz_z[k] + 2.0 * to_xyzz_yyz[k] * tke_0;

            to_0_y_xyzz_zz[k] = 2.0 * to_xyzz_yzz[k] * tke_0;
        }

        // Set up 144-150 components of targeted buffer : GD

        auto to_0_y_xzzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 54);

        auto to_0_y_xzzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 55);

        auto to_0_y_xzzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 56);

        auto to_0_y_xzzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 57);

        auto to_0_y_xzzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 58);

        auto to_0_y_xzzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_y_xzzz_xx, to_0_y_xzzz_xy, to_0_y_xzzz_xz, to_0_y_xzzz_yy, to_0_y_xzzz_yz, to_0_y_xzzz_zz, to_xzzz_x, to_xzzz_xxy, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_y, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzzz_xx[k] = 2.0 * to_xzzz_xxy[k] * tke_0;

            to_0_y_xzzz_xy[k] = -to_xzzz_x[k] + 2.0 * to_xzzz_xyy[k] * tke_0;

            to_0_y_xzzz_xz[k] = 2.0 * to_xzzz_xyz[k] * tke_0;

            to_0_y_xzzz_yy[k] = -2.0 * to_xzzz_y[k] + 2.0 * to_xzzz_yyy[k] * tke_0;

            to_0_y_xzzz_yz[k] = -to_xzzz_z[k] + 2.0 * to_xzzz_yyz[k] * tke_0;

            to_0_y_xzzz_zz[k] = 2.0 * to_xzzz_yzz[k] * tke_0;
        }

        // Set up 150-156 components of targeted buffer : GD

        auto to_0_y_yyyy_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 60);

        auto to_0_y_yyyy_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 61);

        auto to_0_y_yyyy_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 62);

        auto to_0_y_yyyy_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 63);

        auto to_0_y_yyyy_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 64);

        auto to_0_y_yyyy_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_0_y_yyyy_xx, to_0_y_yyyy_xy, to_0_y_yyyy_xz, to_0_y_yyyy_yy, to_0_y_yyyy_yz, to_0_y_yyyy_zz, to_yyyy_x, to_yyyy_xxy, to_yyyy_xyy, to_yyyy_xyz, to_yyyy_y, to_yyyy_yyy, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyy_xx[k] = 2.0 * to_yyyy_xxy[k] * tke_0;

            to_0_y_yyyy_xy[k] = -to_yyyy_x[k] + 2.0 * to_yyyy_xyy[k] * tke_0;

            to_0_y_yyyy_xz[k] = 2.0 * to_yyyy_xyz[k] * tke_0;

            to_0_y_yyyy_yy[k] = -2.0 * to_yyyy_y[k] + 2.0 * to_yyyy_yyy[k] * tke_0;

            to_0_y_yyyy_yz[k] = -to_yyyy_z[k] + 2.0 * to_yyyy_yyz[k] * tke_0;

            to_0_y_yyyy_zz[k] = 2.0 * to_yyyy_yzz[k] * tke_0;
        }

        // Set up 156-162 components of targeted buffer : GD

        auto to_0_y_yyyz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 66);

        auto to_0_y_yyyz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 67);

        auto to_0_y_yyyz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 68);

        auto to_0_y_yyyz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 69);

        auto to_0_y_yyyz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 70);

        auto to_0_y_yyyz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_0_y_yyyz_xx, to_0_y_yyyz_xy, to_0_y_yyyz_xz, to_0_y_yyyz_yy, to_0_y_yyyz_yz, to_0_y_yyyz_zz, to_yyyz_x, to_yyyz_xxy, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_y, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyz_xx[k] = 2.0 * to_yyyz_xxy[k] * tke_0;

            to_0_y_yyyz_xy[k] = -to_yyyz_x[k] + 2.0 * to_yyyz_xyy[k] * tke_0;

            to_0_y_yyyz_xz[k] = 2.0 * to_yyyz_xyz[k] * tke_0;

            to_0_y_yyyz_yy[k] = -2.0 * to_yyyz_y[k] + 2.0 * to_yyyz_yyy[k] * tke_0;

            to_0_y_yyyz_yz[k] = -to_yyyz_z[k] + 2.0 * to_yyyz_yyz[k] * tke_0;

            to_0_y_yyyz_zz[k] = 2.0 * to_yyyz_yzz[k] * tke_0;
        }

        // Set up 162-168 components of targeted buffer : GD

        auto to_0_y_yyzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 72);

        auto to_0_y_yyzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 73);

        auto to_0_y_yyzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 74);

        auto to_0_y_yyzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 75);

        auto to_0_y_yyzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 76);

        auto to_0_y_yyzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_0_y_yyzz_xx, to_0_y_yyzz_xy, to_0_y_yyzz_xz, to_0_y_yyzz_yy, to_0_y_yyzz_yz, to_0_y_yyzz_zz, to_yyzz_x, to_yyzz_xxy, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_y, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyzz_xx[k] = 2.0 * to_yyzz_xxy[k] * tke_0;

            to_0_y_yyzz_xy[k] = -to_yyzz_x[k] + 2.0 * to_yyzz_xyy[k] * tke_0;

            to_0_y_yyzz_xz[k] = 2.0 * to_yyzz_xyz[k] * tke_0;

            to_0_y_yyzz_yy[k] = -2.0 * to_yyzz_y[k] + 2.0 * to_yyzz_yyy[k] * tke_0;

            to_0_y_yyzz_yz[k] = -to_yyzz_z[k] + 2.0 * to_yyzz_yyz[k] * tke_0;

            to_0_y_yyzz_zz[k] = 2.0 * to_yyzz_yzz[k] * tke_0;
        }

        // Set up 168-174 components of targeted buffer : GD

        auto to_0_y_yzzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 78);

        auto to_0_y_yzzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 79);

        auto to_0_y_yzzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 80);

        auto to_0_y_yzzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 81);

        auto to_0_y_yzzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 82);

        auto to_0_y_yzzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_0_y_yzzz_xx, to_0_y_yzzz_xy, to_0_y_yzzz_xz, to_0_y_yzzz_yy, to_0_y_yzzz_yz, to_0_y_yzzz_zz, to_yzzz_x, to_yzzz_xxy, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_y, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzzz_xx[k] = 2.0 * to_yzzz_xxy[k] * tke_0;

            to_0_y_yzzz_xy[k] = -to_yzzz_x[k] + 2.0 * to_yzzz_xyy[k] * tke_0;

            to_0_y_yzzz_xz[k] = 2.0 * to_yzzz_xyz[k] * tke_0;

            to_0_y_yzzz_yy[k] = -2.0 * to_yzzz_y[k] + 2.0 * to_yzzz_yyy[k] * tke_0;

            to_0_y_yzzz_yz[k] = -to_yzzz_z[k] + 2.0 * to_yzzz_yyz[k] * tke_0;

            to_0_y_yzzz_zz[k] = 2.0 * to_yzzz_yzz[k] * tke_0;
        }

        // Set up 174-180 components of targeted buffer : GD

        auto to_0_y_zzzz_xx = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 84);

        auto to_0_y_zzzz_xy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 85);

        auto to_0_y_zzzz_xz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 86);

        auto to_0_y_zzzz_yy = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 87);

        auto to_0_y_zzzz_yz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 88);

        auto to_0_y_zzzz_zz = pbuffer.data(idx_op_geom_001_gd + 1 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_y_zzzz_xx, to_0_y_zzzz_xy, to_0_y_zzzz_xz, to_0_y_zzzz_yy, to_0_y_zzzz_yz, to_0_y_zzzz_zz, to_zzzz_x, to_zzzz_xxy, to_zzzz_xyy, to_zzzz_xyz, to_zzzz_y, to_zzzz_yyy, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzzz_xx[k] = 2.0 * to_zzzz_xxy[k] * tke_0;

            to_0_y_zzzz_xy[k] = -to_zzzz_x[k] + 2.0 * to_zzzz_xyy[k] * tke_0;

            to_0_y_zzzz_xz[k] = 2.0 * to_zzzz_xyz[k] * tke_0;

            to_0_y_zzzz_yy[k] = -2.0 * to_zzzz_y[k] + 2.0 * to_zzzz_yyy[k] * tke_0;

            to_0_y_zzzz_yz[k] = -to_zzzz_z[k] + 2.0 * to_zzzz_yyz[k] * tke_0;

            to_0_y_zzzz_zz[k] = 2.0 * to_zzzz_yzz[k] * tke_0;
        }

        // Set up 180-186 components of targeted buffer : GD

        auto to_0_z_xxxx_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 0);

        auto to_0_z_xxxx_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 1);

        auto to_0_z_xxxx_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 2);

        auto to_0_z_xxxx_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 3);

        auto to_0_z_xxxx_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 4);

        auto to_0_z_xxxx_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_0_z_xxxx_xx, to_0_z_xxxx_xy, to_0_z_xxxx_xz, to_0_z_xxxx_yy, to_0_z_xxxx_yz, to_0_z_xxxx_zz, to_xxxx_x, to_xxxx_xxz, to_xxxx_xyz, to_xxxx_xzz, to_xxxx_y, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_z, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxx_xx[k] = 2.0 * to_xxxx_xxz[k] * tke_0;

            to_0_z_xxxx_xy[k] = 2.0 * to_xxxx_xyz[k] * tke_0;

            to_0_z_xxxx_xz[k] = -to_xxxx_x[k] + 2.0 * to_xxxx_xzz[k] * tke_0;

            to_0_z_xxxx_yy[k] = 2.0 * to_xxxx_yyz[k] * tke_0;

            to_0_z_xxxx_yz[k] = -to_xxxx_y[k] + 2.0 * to_xxxx_yzz[k] * tke_0;

            to_0_z_xxxx_zz[k] = -2.0 * to_xxxx_z[k] + 2.0 * to_xxxx_zzz[k] * tke_0;
        }

        // Set up 186-192 components of targeted buffer : GD

        auto to_0_z_xxxy_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 6);

        auto to_0_z_xxxy_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 7);

        auto to_0_z_xxxy_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 8);

        auto to_0_z_xxxy_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 9);

        auto to_0_z_xxxy_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 10);

        auto to_0_z_xxxy_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_0_z_xxxy_xx, to_0_z_xxxy_xy, to_0_z_xxxy_xz, to_0_z_xxxy_yy, to_0_z_xxxy_yz, to_0_z_xxxy_zz, to_xxxy_x, to_xxxy_xxz, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, to_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxy_xx[k] = 2.0 * to_xxxy_xxz[k] * tke_0;

            to_0_z_xxxy_xy[k] = 2.0 * to_xxxy_xyz[k] * tke_0;

            to_0_z_xxxy_xz[k] = -to_xxxy_x[k] + 2.0 * to_xxxy_xzz[k] * tke_0;

            to_0_z_xxxy_yy[k] = 2.0 * to_xxxy_yyz[k] * tke_0;

            to_0_z_xxxy_yz[k] = -to_xxxy_y[k] + 2.0 * to_xxxy_yzz[k] * tke_0;

            to_0_z_xxxy_zz[k] = -2.0 * to_xxxy_z[k] + 2.0 * to_xxxy_zzz[k] * tke_0;
        }

        // Set up 192-198 components of targeted buffer : GD

        auto to_0_z_xxxz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 12);

        auto to_0_z_xxxz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 13);

        auto to_0_z_xxxz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 14);

        auto to_0_z_xxxz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 15);

        auto to_0_z_xxxz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 16);

        auto to_0_z_xxxz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_0_z_xxxz_xx, to_0_z_xxxz_xy, to_0_z_xxxz_xz, to_0_z_xxxz_yy, to_0_z_xxxz_yz, to_0_z_xxxz_zz, to_xxxz_x, to_xxxz_xxz, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, to_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxz_xx[k] = 2.0 * to_xxxz_xxz[k] * tke_0;

            to_0_z_xxxz_xy[k] = 2.0 * to_xxxz_xyz[k] * tke_0;

            to_0_z_xxxz_xz[k] = -to_xxxz_x[k] + 2.0 * to_xxxz_xzz[k] * tke_0;

            to_0_z_xxxz_yy[k] = 2.0 * to_xxxz_yyz[k] * tke_0;

            to_0_z_xxxz_yz[k] = -to_xxxz_y[k] + 2.0 * to_xxxz_yzz[k] * tke_0;

            to_0_z_xxxz_zz[k] = -2.0 * to_xxxz_z[k] + 2.0 * to_xxxz_zzz[k] * tke_0;
        }

        // Set up 198-204 components of targeted buffer : GD

        auto to_0_z_xxyy_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 18);

        auto to_0_z_xxyy_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 19);

        auto to_0_z_xxyy_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 20);

        auto to_0_z_xxyy_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 21);

        auto to_0_z_xxyy_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 22);

        auto to_0_z_xxyy_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_0_z_xxyy_xx, to_0_z_xxyy_xy, to_0_z_xxyy_xz, to_0_z_xxyy_yy, to_0_z_xxyy_yz, to_0_z_xxyy_zz, to_xxyy_x, to_xxyy_xxz, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, to_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyy_xx[k] = 2.0 * to_xxyy_xxz[k] * tke_0;

            to_0_z_xxyy_xy[k] = 2.0 * to_xxyy_xyz[k] * tke_0;

            to_0_z_xxyy_xz[k] = -to_xxyy_x[k] + 2.0 * to_xxyy_xzz[k] * tke_0;

            to_0_z_xxyy_yy[k] = 2.0 * to_xxyy_yyz[k] * tke_0;

            to_0_z_xxyy_yz[k] = -to_xxyy_y[k] + 2.0 * to_xxyy_yzz[k] * tke_0;

            to_0_z_xxyy_zz[k] = -2.0 * to_xxyy_z[k] + 2.0 * to_xxyy_zzz[k] * tke_0;
        }

        // Set up 204-210 components of targeted buffer : GD

        auto to_0_z_xxyz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 24);

        auto to_0_z_xxyz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 25);

        auto to_0_z_xxyz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 26);

        auto to_0_z_xxyz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 27);

        auto to_0_z_xxyz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 28);

        auto to_0_z_xxyz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_0_z_xxyz_xx, to_0_z_xxyz_xy, to_0_z_xxyz_xz, to_0_z_xxyz_yy, to_0_z_xxyz_yz, to_0_z_xxyz_zz, to_xxyz_x, to_xxyz_xxz, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyz_xx[k] = 2.0 * to_xxyz_xxz[k] * tke_0;

            to_0_z_xxyz_xy[k] = 2.0 * to_xxyz_xyz[k] * tke_0;

            to_0_z_xxyz_xz[k] = -to_xxyz_x[k] + 2.0 * to_xxyz_xzz[k] * tke_0;

            to_0_z_xxyz_yy[k] = 2.0 * to_xxyz_yyz[k] * tke_0;

            to_0_z_xxyz_yz[k] = -to_xxyz_y[k] + 2.0 * to_xxyz_yzz[k] * tke_0;

            to_0_z_xxyz_zz[k] = -2.0 * to_xxyz_z[k] + 2.0 * to_xxyz_zzz[k] * tke_0;
        }

        // Set up 210-216 components of targeted buffer : GD

        auto to_0_z_xxzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 30);

        auto to_0_z_xxzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 31);

        auto to_0_z_xxzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 32);

        auto to_0_z_xxzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 33);

        auto to_0_z_xxzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 34);

        auto to_0_z_xxzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_0_z_xxzz_xx, to_0_z_xxzz_xy, to_0_z_xxzz_xz, to_0_z_xxzz_yy, to_0_z_xxzz_yz, to_0_z_xxzz_zz, to_xxzz_x, to_xxzz_xxz, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, to_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxzz_xx[k] = 2.0 * to_xxzz_xxz[k] * tke_0;

            to_0_z_xxzz_xy[k] = 2.0 * to_xxzz_xyz[k] * tke_0;

            to_0_z_xxzz_xz[k] = -to_xxzz_x[k] + 2.0 * to_xxzz_xzz[k] * tke_0;

            to_0_z_xxzz_yy[k] = 2.0 * to_xxzz_yyz[k] * tke_0;

            to_0_z_xxzz_yz[k] = -to_xxzz_y[k] + 2.0 * to_xxzz_yzz[k] * tke_0;

            to_0_z_xxzz_zz[k] = -2.0 * to_xxzz_z[k] + 2.0 * to_xxzz_zzz[k] * tke_0;
        }

        // Set up 216-222 components of targeted buffer : GD

        auto to_0_z_xyyy_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 36);

        auto to_0_z_xyyy_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 37);

        auto to_0_z_xyyy_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 38);

        auto to_0_z_xyyy_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 39);

        auto to_0_z_xyyy_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 40);

        auto to_0_z_xyyy_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_0_z_xyyy_xx, to_0_z_xyyy_xy, to_0_z_xyyy_xz, to_0_z_xyyy_yy, to_0_z_xyyy_yz, to_0_z_xyyy_zz, to_xyyy_x, to_xyyy_xxz, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyy_xx[k] = 2.0 * to_xyyy_xxz[k] * tke_0;

            to_0_z_xyyy_xy[k] = 2.0 * to_xyyy_xyz[k] * tke_0;

            to_0_z_xyyy_xz[k] = -to_xyyy_x[k] + 2.0 * to_xyyy_xzz[k] * tke_0;

            to_0_z_xyyy_yy[k] = 2.0 * to_xyyy_yyz[k] * tke_0;

            to_0_z_xyyy_yz[k] = -to_xyyy_y[k] + 2.0 * to_xyyy_yzz[k] * tke_0;

            to_0_z_xyyy_zz[k] = -2.0 * to_xyyy_z[k] + 2.0 * to_xyyy_zzz[k] * tke_0;
        }

        // Set up 222-228 components of targeted buffer : GD

        auto to_0_z_xyyz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 42);

        auto to_0_z_xyyz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 43);

        auto to_0_z_xyyz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 44);

        auto to_0_z_xyyz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 45);

        auto to_0_z_xyyz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 46);

        auto to_0_z_xyyz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_0_z_xyyz_xx, to_0_z_xyyz_xy, to_0_z_xyyz_xz, to_0_z_xyyz_yy, to_0_z_xyyz_yz, to_0_z_xyyz_zz, to_xyyz_x, to_xyyz_xxz, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyz_xx[k] = 2.0 * to_xyyz_xxz[k] * tke_0;

            to_0_z_xyyz_xy[k] = 2.0 * to_xyyz_xyz[k] * tke_0;

            to_0_z_xyyz_xz[k] = -to_xyyz_x[k] + 2.0 * to_xyyz_xzz[k] * tke_0;

            to_0_z_xyyz_yy[k] = 2.0 * to_xyyz_yyz[k] * tke_0;

            to_0_z_xyyz_yz[k] = -to_xyyz_y[k] + 2.0 * to_xyyz_yzz[k] * tke_0;

            to_0_z_xyyz_zz[k] = -2.0 * to_xyyz_z[k] + 2.0 * to_xyyz_zzz[k] * tke_0;
        }

        // Set up 228-234 components of targeted buffer : GD

        auto to_0_z_xyzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 48);

        auto to_0_z_xyzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 49);

        auto to_0_z_xyzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 50);

        auto to_0_z_xyzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 51);

        auto to_0_z_xyzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 52);

        auto to_0_z_xyzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_0_z_xyzz_xx, to_0_z_xyzz_xy, to_0_z_xyzz_xz, to_0_z_xyzz_yy, to_0_z_xyzz_yz, to_0_z_xyzz_zz, to_xyzz_x, to_xyzz_xxz, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyzz_xx[k] = 2.0 * to_xyzz_xxz[k] * tke_0;

            to_0_z_xyzz_xy[k] = 2.0 * to_xyzz_xyz[k] * tke_0;

            to_0_z_xyzz_xz[k] = -to_xyzz_x[k] + 2.0 * to_xyzz_xzz[k] * tke_0;

            to_0_z_xyzz_yy[k] = 2.0 * to_xyzz_yyz[k] * tke_0;

            to_0_z_xyzz_yz[k] = -to_xyzz_y[k] + 2.0 * to_xyzz_yzz[k] * tke_0;

            to_0_z_xyzz_zz[k] = -2.0 * to_xyzz_z[k] + 2.0 * to_xyzz_zzz[k] * tke_0;
        }

        // Set up 234-240 components of targeted buffer : GD

        auto to_0_z_xzzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 54);

        auto to_0_z_xzzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 55);

        auto to_0_z_xzzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 56);

        auto to_0_z_xzzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 57);

        auto to_0_z_xzzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 58);

        auto to_0_z_xzzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_0_z_xzzz_xx, to_0_z_xzzz_xy, to_0_z_xzzz_xz, to_0_z_xzzz_yy, to_0_z_xzzz_yz, to_0_z_xzzz_zz, to_xzzz_x, to_xzzz_xxz, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzzz_xx[k] = 2.0 * to_xzzz_xxz[k] * tke_0;

            to_0_z_xzzz_xy[k] = 2.0 * to_xzzz_xyz[k] * tke_0;

            to_0_z_xzzz_xz[k] = -to_xzzz_x[k] + 2.0 * to_xzzz_xzz[k] * tke_0;

            to_0_z_xzzz_yy[k] = 2.0 * to_xzzz_yyz[k] * tke_0;

            to_0_z_xzzz_yz[k] = -to_xzzz_y[k] + 2.0 * to_xzzz_yzz[k] * tke_0;

            to_0_z_xzzz_zz[k] = -2.0 * to_xzzz_z[k] + 2.0 * to_xzzz_zzz[k] * tke_0;
        }

        // Set up 240-246 components of targeted buffer : GD

        auto to_0_z_yyyy_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 60);

        auto to_0_z_yyyy_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 61);

        auto to_0_z_yyyy_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 62);

        auto to_0_z_yyyy_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 63);

        auto to_0_z_yyyy_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 64);

        auto to_0_z_yyyy_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_0_z_yyyy_xx, to_0_z_yyyy_xy, to_0_z_yyyy_xz, to_0_z_yyyy_yy, to_0_z_yyyy_yz, to_0_z_yyyy_zz, to_yyyy_x, to_yyyy_xxz, to_yyyy_xyz, to_yyyy_xzz, to_yyyy_y, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_z, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyy_xx[k] = 2.0 * to_yyyy_xxz[k] * tke_0;

            to_0_z_yyyy_xy[k] = 2.0 * to_yyyy_xyz[k] * tke_0;

            to_0_z_yyyy_xz[k] = -to_yyyy_x[k] + 2.0 * to_yyyy_xzz[k] * tke_0;

            to_0_z_yyyy_yy[k] = 2.0 * to_yyyy_yyz[k] * tke_0;

            to_0_z_yyyy_yz[k] = -to_yyyy_y[k] + 2.0 * to_yyyy_yzz[k] * tke_0;

            to_0_z_yyyy_zz[k] = -2.0 * to_yyyy_z[k] + 2.0 * to_yyyy_zzz[k] * tke_0;
        }

        // Set up 246-252 components of targeted buffer : GD

        auto to_0_z_yyyz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 66);

        auto to_0_z_yyyz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 67);

        auto to_0_z_yyyz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 68);

        auto to_0_z_yyyz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 69);

        auto to_0_z_yyyz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 70);

        auto to_0_z_yyyz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_0_z_yyyz_xx, to_0_z_yyyz_xy, to_0_z_yyyz_xz, to_0_z_yyyz_yy, to_0_z_yyyz_yz, to_0_z_yyyz_zz, to_yyyz_x, to_yyyz_xxz, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, to_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyz_xx[k] = 2.0 * to_yyyz_xxz[k] * tke_0;

            to_0_z_yyyz_xy[k] = 2.0 * to_yyyz_xyz[k] * tke_0;

            to_0_z_yyyz_xz[k] = -to_yyyz_x[k] + 2.0 * to_yyyz_xzz[k] * tke_0;

            to_0_z_yyyz_yy[k] = 2.0 * to_yyyz_yyz[k] * tke_0;

            to_0_z_yyyz_yz[k] = -to_yyyz_y[k] + 2.0 * to_yyyz_yzz[k] * tke_0;

            to_0_z_yyyz_zz[k] = -2.0 * to_yyyz_z[k] + 2.0 * to_yyyz_zzz[k] * tke_0;
        }

        // Set up 252-258 components of targeted buffer : GD

        auto to_0_z_yyzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 72);

        auto to_0_z_yyzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 73);

        auto to_0_z_yyzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 74);

        auto to_0_z_yyzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 75);

        auto to_0_z_yyzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 76);

        auto to_0_z_yyzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_0_z_yyzz_xx, to_0_z_yyzz_xy, to_0_z_yyzz_xz, to_0_z_yyzz_yy, to_0_z_yyzz_yz, to_0_z_yyzz_zz, to_yyzz_x, to_yyzz_xxz, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, to_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyzz_xx[k] = 2.0 * to_yyzz_xxz[k] * tke_0;

            to_0_z_yyzz_xy[k] = 2.0 * to_yyzz_xyz[k] * tke_0;

            to_0_z_yyzz_xz[k] = -to_yyzz_x[k] + 2.0 * to_yyzz_xzz[k] * tke_0;

            to_0_z_yyzz_yy[k] = 2.0 * to_yyzz_yyz[k] * tke_0;

            to_0_z_yyzz_yz[k] = -to_yyzz_y[k] + 2.0 * to_yyzz_yzz[k] * tke_0;

            to_0_z_yyzz_zz[k] = -2.0 * to_yyzz_z[k] + 2.0 * to_yyzz_zzz[k] * tke_0;
        }

        // Set up 258-264 components of targeted buffer : GD

        auto to_0_z_yzzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 78);

        auto to_0_z_yzzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 79);

        auto to_0_z_yzzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 80);

        auto to_0_z_yzzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 81);

        auto to_0_z_yzzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 82);

        auto to_0_z_yzzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_0_z_yzzz_xx, to_0_z_yzzz_xy, to_0_z_yzzz_xz, to_0_z_yzzz_yy, to_0_z_yzzz_yz, to_0_z_yzzz_zz, to_yzzz_x, to_yzzz_xxz, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzzz_xx[k] = 2.0 * to_yzzz_xxz[k] * tke_0;

            to_0_z_yzzz_xy[k] = 2.0 * to_yzzz_xyz[k] * tke_0;

            to_0_z_yzzz_xz[k] = -to_yzzz_x[k] + 2.0 * to_yzzz_xzz[k] * tke_0;

            to_0_z_yzzz_yy[k] = 2.0 * to_yzzz_yyz[k] * tke_0;

            to_0_z_yzzz_yz[k] = -to_yzzz_y[k] + 2.0 * to_yzzz_yzz[k] * tke_0;

            to_0_z_yzzz_zz[k] = -2.0 * to_yzzz_z[k] + 2.0 * to_yzzz_zzz[k] * tke_0;
        }

        // Set up 264-270 components of targeted buffer : GD

        auto to_0_z_zzzz_xx = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 84);

        auto to_0_z_zzzz_xy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 85);

        auto to_0_z_zzzz_xz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 86);

        auto to_0_z_zzzz_yy = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 87);

        auto to_0_z_zzzz_yz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 88);

        auto to_0_z_zzzz_zz = pbuffer.data(idx_op_geom_001_gd + 2 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_0_z_zzzz_xx, to_0_z_zzzz_xy, to_0_z_zzzz_xz, to_0_z_zzzz_yy, to_0_z_zzzz_yz, to_0_z_zzzz_zz, to_zzzz_x, to_zzzz_xxz, to_zzzz_xyz, to_zzzz_xzz, to_zzzz_y, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_z, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzzz_xx[k] = 2.0 * to_zzzz_xxz[k] * tke_0;

            to_0_z_zzzz_xy[k] = 2.0 * to_zzzz_xyz[k] * tke_0;

            to_0_z_zzzz_xz[k] = -to_zzzz_x[k] + 2.0 * to_zzzz_xzz[k] * tke_0;

            to_0_z_zzzz_yy[k] = 2.0 * to_zzzz_yyz[k] * tke_0;

            to_0_z_zzzz_yz[k] = -to_zzzz_y[k] + 2.0 * to_zzzz_yzz[k] * tke_0;

            to_0_z_zzzz_zz[k] = -2.0 * to_zzzz_z[k] + 2.0 * to_zzzz_zzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

