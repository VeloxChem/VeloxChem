#include "LocalCorePotentialPrimRecDG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_dg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_dg,
                                  const size_t idx_sg,
                                  const size_t idx_pf,
                                  const size_t idx_pg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxxy = pbuffer.data(idx_sg + 1);

    auto tg_0_xxxz = pbuffer.data(idx_sg + 2);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxyz = pbuffer.data(idx_sg + 4);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xyyz = pbuffer.data(idx_sg + 7);

    auto tg_0_xyzz = pbuffer.data(idx_sg + 8);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyyz = pbuffer.data(idx_sg + 11);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx = pbuffer.data(idx_pf);

    auto tg_x_xxy = pbuffer.data(idx_pf + 1);

    auto tg_x_xxz = pbuffer.data(idx_pf + 2);

    auto tg_x_xyy = pbuffer.data(idx_pf + 3);

    auto tg_x_xyz = pbuffer.data(idx_pf + 4);

    auto tg_x_xzz = pbuffer.data(idx_pf + 5);

    auto tg_x_yyy = pbuffer.data(idx_pf + 6);

    auto tg_x_yyz = pbuffer.data(idx_pf + 7);

    auto tg_x_yzz = pbuffer.data(idx_pf + 8);

    auto tg_x_zzz = pbuffer.data(idx_pf + 9);

    auto tg_y_xxx = pbuffer.data(idx_pf + 10);

    auto tg_y_xxy = pbuffer.data(idx_pf + 11);

    auto tg_y_xxz = pbuffer.data(idx_pf + 12);

    auto tg_y_xyy = pbuffer.data(idx_pf + 13);

    auto tg_y_xyz = pbuffer.data(idx_pf + 14);

    auto tg_y_xzz = pbuffer.data(idx_pf + 15);

    auto tg_y_yyy = pbuffer.data(idx_pf + 16);

    auto tg_y_yyz = pbuffer.data(idx_pf + 17);

    auto tg_y_yzz = pbuffer.data(idx_pf + 18);

    auto tg_y_zzz = pbuffer.data(idx_pf + 19);

    auto tg_z_xxx = pbuffer.data(idx_pf + 20);

    auto tg_z_xxy = pbuffer.data(idx_pf + 21);

    auto tg_z_xxz = pbuffer.data(idx_pf + 22);

    auto tg_z_xyy = pbuffer.data(idx_pf + 23);

    auto tg_z_xyz = pbuffer.data(idx_pf + 24);

    auto tg_z_xzz = pbuffer.data(idx_pf + 25);

    auto tg_z_yyy = pbuffer.data(idx_pf + 26);

    auto tg_z_yyz = pbuffer.data(idx_pf + 27);

    auto tg_z_yzz = pbuffer.data(idx_pf + 28);

    auto tg_z_zzz = pbuffer.data(idx_pf + 29);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx = pbuffer.data(idx_pg);

    auto tg_x_xxxy = pbuffer.data(idx_pg + 1);

    auto tg_x_xxxz = pbuffer.data(idx_pg + 2);

    auto tg_x_xxyy = pbuffer.data(idx_pg + 3);

    auto tg_x_xxyz = pbuffer.data(idx_pg + 4);

    auto tg_x_xxzz = pbuffer.data(idx_pg + 5);

    auto tg_x_xyyy = pbuffer.data(idx_pg + 6);

    auto tg_x_xyyz = pbuffer.data(idx_pg + 7);

    auto tg_x_xyzz = pbuffer.data(idx_pg + 8);

    auto tg_x_xzzz = pbuffer.data(idx_pg + 9);

    auto tg_x_yyyy = pbuffer.data(idx_pg + 10);

    auto tg_x_yyyz = pbuffer.data(idx_pg + 11);

    auto tg_x_yyzz = pbuffer.data(idx_pg + 12);

    auto tg_x_yzzz = pbuffer.data(idx_pg + 13);

    auto tg_x_zzzz = pbuffer.data(idx_pg + 14);

    auto tg_y_xxxx = pbuffer.data(idx_pg + 15);

    auto tg_y_xxxy = pbuffer.data(idx_pg + 16);

    auto tg_y_xxxz = pbuffer.data(idx_pg + 17);

    auto tg_y_xxyy = pbuffer.data(idx_pg + 18);

    auto tg_y_xxyz = pbuffer.data(idx_pg + 19);

    auto tg_y_xxzz = pbuffer.data(idx_pg + 20);

    auto tg_y_xyyy = pbuffer.data(idx_pg + 21);

    auto tg_y_xyyz = pbuffer.data(idx_pg + 22);

    auto tg_y_xyzz = pbuffer.data(idx_pg + 23);

    auto tg_y_xzzz = pbuffer.data(idx_pg + 24);

    auto tg_y_yyyy = pbuffer.data(idx_pg + 25);

    auto tg_y_yyyz = pbuffer.data(idx_pg + 26);

    auto tg_y_yyzz = pbuffer.data(idx_pg + 27);

    auto tg_y_yzzz = pbuffer.data(idx_pg + 28);

    auto tg_y_zzzz = pbuffer.data(idx_pg + 29);

    auto tg_z_xxxx = pbuffer.data(idx_pg + 30);

    auto tg_z_xxxy = pbuffer.data(idx_pg + 31);

    auto tg_z_xxxz = pbuffer.data(idx_pg + 32);

    auto tg_z_xxyy = pbuffer.data(idx_pg + 33);

    auto tg_z_xxyz = pbuffer.data(idx_pg + 34);

    auto tg_z_xxzz = pbuffer.data(idx_pg + 35);

    auto tg_z_xyyy = pbuffer.data(idx_pg + 36);

    auto tg_z_xyyz = pbuffer.data(idx_pg + 37);

    auto tg_z_xyzz = pbuffer.data(idx_pg + 38);

    auto tg_z_xzzz = pbuffer.data(idx_pg + 39);

    auto tg_z_yyyy = pbuffer.data(idx_pg + 40);

    auto tg_z_yyyz = pbuffer.data(idx_pg + 41);

    auto tg_z_yyzz = pbuffer.data(idx_pg + 42);

    auto tg_z_yzzz = pbuffer.data(idx_pg + 43);

    auto tg_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up components of targeted buffer : DG

    auto tg_xx_xxxx = pbuffer.data(idx_dg);

    auto tg_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto tg_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto tg_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto tg_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto tg_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto tg_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto tg_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto tg_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto tg_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto tg_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto tg_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto tg_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto tg_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto tg_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto tg_xy_xxxx = pbuffer.data(idx_dg + 15);

    auto tg_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto tg_xy_xxxz = pbuffer.data(idx_dg + 17);

    auto tg_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto tg_xy_xxyz = pbuffer.data(idx_dg + 19);

    auto tg_xy_xxzz = pbuffer.data(idx_dg + 20);

    auto tg_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto tg_xy_xyyz = pbuffer.data(idx_dg + 22);

    auto tg_xy_xyzz = pbuffer.data(idx_dg + 23);

    auto tg_xy_xzzz = pbuffer.data(idx_dg + 24);

    auto tg_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto tg_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto tg_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto tg_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto tg_xy_zzzz = pbuffer.data(idx_dg + 29);

    auto tg_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto tg_xz_xxxy = pbuffer.data(idx_dg + 31);

    auto tg_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto tg_xz_xxyy = pbuffer.data(idx_dg + 33);

    auto tg_xz_xxyz = pbuffer.data(idx_dg + 34);

    auto tg_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto tg_xz_xyyy = pbuffer.data(idx_dg + 36);

    auto tg_xz_xyyz = pbuffer.data(idx_dg + 37);

    auto tg_xz_xyzz = pbuffer.data(idx_dg + 38);

    auto tg_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto tg_xz_yyyy = pbuffer.data(idx_dg + 40);

    auto tg_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto tg_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto tg_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto tg_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto tg_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto tg_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto tg_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto tg_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto tg_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto tg_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto tg_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto tg_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto tg_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto tg_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto tg_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto tg_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto tg_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto tg_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto tg_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto tg_yz_xxxx = pbuffer.data(idx_dg + 60);

    auto tg_yz_xxxy = pbuffer.data(idx_dg + 61);

    auto tg_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto tg_yz_xxyy = pbuffer.data(idx_dg + 63);

    auto tg_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto tg_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto tg_yz_xyyy = pbuffer.data(idx_dg + 66);

    auto tg_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto tg_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto tg_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto tg_yz_yyyy = pbuffer.data(idx_dg + 70);

    auto tg_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto tg_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto tg_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto tg_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto tg_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto tg_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto tg_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto tg_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto tg_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto tg_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto tg_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto tg_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto tg_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto tg_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto tg_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto tg_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto tg_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto tg_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto tg_zz_zzzz = pbuffer.data(idx_dg + 89);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxxx, tg_0_xxxy, tg_0_xxxz, tg_0_xxyy, tg_0_xxyz, tg_0_xxzz, tg_0_xyyy, tg_0_xyyz, tg_0_xyzz, tg_0_xzzz, tg_0_yyyy, tg_0_yyyz, tg_0_yyzz, tg_0_yzzz, tg_0_zzzz, tg_x_xxx, tg_x_xxxx, tg_x_xxxy, tg_x_xxxz, tg_x_xxy, tg_x_xxyy, tg_x_xxyz, tg_x_xxz, tg_x_xxzz, tg_x_xyy, tg_x_xyyy, tg_x_xyyz, tg_x_xyz, tg_x_xyzz, tg_x_xzz, tg_x_xzzz, tg_x_yyy, tg_x_yyyy, tg_x_yyyz, tg_x_yyz, tg_x_yyzz, tg_x_yzz, tg_x_yzzz, tg_x_zzz, tg_x_zzzz, tg_xx_xxxx, tg_xx_xxxy, tg_xx_xxxz, tg_xx_xxyy, tg_xx_xxyz, tg_xx_xxzz, tg_xx_xyyy, tg_xx_xyyz, tg_xx_xyzz, tg_xx_xzzz, tg_xx_yyyy, tg_xx_yyyz, tg_xx_yyzz, tg_xx_yzzz, tg_xx_zzzz, tg_xy_xxxx, tg_xy_xxxy, tg_xy_xxxz, tg_xy_xxyy, tg_xy_xxyz, tg_xy_xxzz, tg_xy_xyyy, tg_xy_xyyz, tg_xy_xyzz, tg_xy_xzzz, tg_xy_yyyy, tg_xy_yyyz, tg_xy_yyzz, tg_xy_yzzz, tg_xy_zzzz, tg_xz_xxxx, tg_xz_xxxy, tg_xz_xxxz, tg_xz_xxyy, tg_xz_xxyz, tg_xz_xxzz, tg_xz_xyyy, tg_xz_xyyz, tg_xz_xyzz, tg_xz_xzzz, tg_xz_yyyy, tg_xz_yyyz, tg_xz_yyzz, tg_xz_yzzz, tg_xz_zzzz, tg_y_xxx, tg_y_xxxx, tg_y_xxxy, tg_y_xxxz, tg_y_xxy, tg_y_xxyy, tg_y_xxyz, tg_y_xxz, tg_y_xxzz, tg_y_xyy, tg_y_xyyy, tg_y_xyyz, tg_y_xyz, tg_y_xyzz, tg_y_xzz, tg_y_xzzz, tg_y_yyy, tg_y_yyyy, tg_y_yyyz, tg_y_yyz, tg_y_yyzz, tg_y_yzz, tg_y_yzzz, tg_y_zzz, tg_y_zzzz, tg_yy_xxxx, tg_yy_xxxy, tg_yy_xxxz, tg_yy_xxyy, tg_yy_xxyz, tg_yy_xxzz, tg_yy_xyyy, tg_yy_xyyz, tg_yy_xyzz, tg_yy_xzzz, tg_yy_yyyy, tg_yy_yyyz, tg_yy_yyzz, tg_yy_yzzz, tg_yy_zzzz, tg_yz_xxxx, tg_yz_xxxy, tg_yz_xxxz, tg_yz_xxyy, tg_yz_xxyz, tg_yz_xxzz, tg_yz_xyyy, tg_yz_xyyz, tg_yz_xyzz, tg_yz_xzzz, tg_yz_yyyy, tg_yz_yyyz, tg_yz_yyzz, tg_yz_yzzz, tg_yz_zzzz, tg_z_xxx, tg_z_xxxx, tg_z_xxxy, tg_z_xxxz, tg_z_xxy, tg_z_xxyy, tg_z_xxyz, tg_z_xxz, tg_z_xxzz, tg_z_xyy, tg_z_xyyy, tg_z_xyyz, tg_z_xyz, tg_z_xyzz, tg_z_xzz, tg_z_xzzz, tg_z_yyy, tg_z_yyyy, tg_z_yyyz, tg_z_yyz, tg_z_yyzz, tg_z_yzz, tg_z_yzzz, tg_z_zzz, tg_z_zzzz, tg_zz_xxxx, tg_zz_xxxy, tg_zz_xxxz, tg_zz_xxyy, tg_zz_xxyz, tg_zz_xxzz, tg_zz_xyyy, tg_zz_xyyz, tg_zz_xyzz, tg_zz_xzzz, tg_zz_yyyy, tg_zz_yyyz, tg_zz_yyzz, tg_zz_yzzz, tg_zz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_xxxx[i] = tg_0_xxxx[i] * fxi[i] + 4.0 * tg_x_xxx[i] * fxi[i] + tg_x_xxxx[i] * ra_x[i];

        tg_xx_xxxy[i] = tg_0_xxxy[i] * fxi[i] + 3.0 * tg_x_xxy[i] * fxi[i] + tg_x_xxxy[i] * ra_x[i];

        tg_xx_xxxz[i] = tg_0_xxxz[i] * fxi[i] + 3.0 * tg_x_xxz[i] * fxi[i] + tg_x_xxxz[i] * ra_x[i];

        tg_xx_xxyy[i] = tg_0_xxyy[i] * fxi[i] + 2.0 * tg_x_xyy[i] * fxi[i] + tg_x_xxyy[i] * ra_x[i];

        tg_xx_xxyz[i] = tg_0_xxyz[i] * fxi[i] + 2.0 * tg_x_xyz[i] * fxi[i] + tg_x_xxyz[i] * ra_x[i];

        tg_xx_xxzz[i] = tg_0_xxzz[i] * fxi[i] + 2.0 * tg_x_xzz[i] * fxi[i] + tg_x_xxzz[i] * ra_x[i];

        tg_xx_xyyy[i] = tg_0_xyyy[i] * fxi[i] + tg_x_yyy[i] * fxi[i] + tg_x_xyyy[i] * ra_x[i];

        tg_xx_xyyz[i] = tg_0_xyyz[i] * fxi[i] + tg_x_yyz[i] * fxi[i] + tg_x_xyyz[i] * ra_x[i];

        tg_xx_xyzz[i] = tg_0_xyzz[i] * fxi[i] + tg_x_yzz[i] * fxi[i] + tg_x_xyzz[i] * ra_x[i];

        tg_xx_xzzz[i] = tg_0_xzzz[i] * fxi[i] + tg_x_zzz[i] * fxi[i] + tg_x_xzzz[i] * ra_x[i];

        tg_xx_yyyy[i] = tg_0_yyyy[i] * fxi[i] + tg_x_yyyy[i] * ra_x[i];

        tg_xx_yyyz[i] = tg_0_yyyz[i] * fxi[i] + tg_x_yyyz[i] * ra_x[i];

        tg_xx_yyzz[i] = tg_0_yyzz[i] * fxi[i] + tg_x_yyzz[i] * ra_x[i];

        tg_xx_yzzz[i] = tg_0_yzzz[i] * fxi[i] + tg_x_yzzz[i] * ra_x[i];

        tg_xx_zzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_x_zzzz[i] * ra_x[i];

        tg_xy_xxxx[i] = tg_x_xxxx[i] * ra_y[i];

        tg_xy_xxxy[i] = 3.0 * tg_y_xxy[i] * fxi[i] + tg_y_xxxy[i] * ra_x[i];

        tg_xy_xxxz[i] = tg_x_xxxz[i] * ra_y[i];

        tg_xy_xxyy[i] = 2.0 * tg_y_xyy[i] * fxi[i] + tg_y_xxyy[i] * ra_x[i];

        tg_xy_xxyz[i] = 2.0 * tg_y_xyz[i] * fxi[i] + tg_y_xxyz[i] * ra_x[i];

        tg_xy_xxzz[i] = tg_x_xxzz[i] * ra_y[i];

        tg_xy_xyyy[i] = tg_y_yyy[i] * fxi[i] + tg_y_xyyy[i] * ra_x[i];

        tg_xy_xyyz[i] = tg_y_yyz[i] * fxi[i] + tg_y_xyyz[i] * ra_x[i];

        tg_xy_xyzz[i] = tg_y_yzz[i] * fxi[i] + tg_y_xyzz[i] * ra_x[i];

        tg_xy_xzzz[i] = tg_x_xzzz[i] * ra_y[i];

        tg_xy_yyyy[i] = tg_y_yyyy[i] * ra_x[i];

        tg_xy_yyyz[i] = tg_y_yyyz[i] * ra_x[i];

        tg_xy_yyzz[i] = tg_y_yyzz[i] * ra_x[i];

        tg_xy_yzzz[i] = tg_y_yzzz[i] * ra_x[i];

        tg_xy_zzzz[i] = tg_y_zzzz[i] * ra_x[i];

        tg_xz_xxxx[i] = tg_x_xxxx[i] * ra_z[i];

        tg_xz_xxxy[i] = tg_x_xxxy[i] * ra_z[i];

        tg_xz_xxxz[i] = 3.0 * tg_z_xxz[i] * fxi[i] + tg_z_xxxz[i] * ra_x[i];

        tg_xz_xxyy[i] = tg_x_xxyy[i] * ra_z[i];

        tg_xz_xxyz[i] = 2.0 * tg_z_xyz[i] * fxi[i] + tg_z_xxyz[i] * ra_x[i];

        tg_xz_xxzz[i] = 2.0 * tg_z_xzz[i] * fxi[i] + tg_z_xxzz[i] * ra_x[i];

        tg_xz_xyyy[i] = tg_x_xyyy[i] * ra_z[i];

        tg_xz_xyyz[i] = tg_z_yyz[i] * fxi[i] + tg_z_xyyz[i] * ra_x[i];

        tg_xz_xyzz[i] = tg_z_yzz[i] * fxi[i] + tg_z_xyzz[i] * ra_x[i];

        tg_xz_xzzz[i] = tg_z_zzz[i] * fxi[i] + tg_z_xzzz[i] * ra_x[i];

        tg_xz_yyyy[i] = tg_z_yyyy[i] * ra_x[i];

        tg_xz_yyyz[i] = tg_z_yyyz[i] * ra_x[i];

        tg_xz_yyzz[i] = tg_z_yyzz[i] * ra_x[i];

        tg_xz_yzzz[i] = tg_z_yzzz[i] * ra_x[i];

        tg_xz_zzzz[i] = tg_z_zzzz[i] * ra_x[i];

        tg_yy_xxxx[i] = tg_0_xxxx[i] * fxi[i] + tg_y_xxxx[i] * ra_y[i];

        tg_yy_xxxy[i] = tg_0_xxxy[i] * fxi[i] + tg_y_xxx[i] * fxi[i] + tg_y_xxxy[i] * ra_y[i];

        tg_yy_xxxz[i] = tg_0_xxxz[i] * fxi[i] + tg_y_xxxz[i] * ra_y[i];

        tg_yy_xxyy[i] = tg_0_xxyy[i] * fxi[i] + 2.0 * tg_y_xxy[i] * fxi[i] + tg_y_xxyy[i] * ra_y[i];

        tg_yy_xxyz[i] = tg_0_xxyz[i] * fxi[i] + tg_y_xxz[i] * fxi[i] + tg_y_xxyz[i] * ra_y[i];

        tg_yy_xxzz[i] = tg_0_xxzz[i] * fxi[i] + tg_y_xxzz[i] * ra_y[i];

        tg_yy_xyyy[i] = tg_0_xyyy[i] * fxi[i] + 3.0 * tg_y_xyy[i] * fxi[i] + tg_y_xyyy[i] * ra_y[i];

        tg_yy_xyyz[i] = tg_0_xyyz[i] * fxi[i] + 2.0 * tg_y_xyz[i] * fxi[i] + tg_y_xyyz[i] * ra_y[i];

        tg_yy_xyzz[i] = tg_0_xyzz[i] * fxi[i] + tg_y_xzz[i] * fxi[i] + tg_y_xyzz[i] * ra_y[i];

        tg_yy_xzzz[i] = tg_0_xzzz[i] * fxi[i] + tg_y_xzzz[i] * ra_y[i];

        tg_yy_yyyy[i] = tg_0_yyyy[i] * fxi[i] + 4.0 * tg_y_yyy[i] * fxi[i] + tg_y_yyyy[i] * ra_y[i];

        tg_yy_yyyz[i] = tg_0_yyyz[i] * fxi[i] + 3.0 * tg_y_yyz[i] * fxi[i] + tg_y_yyyz[i] * ra_y[i];

        tg_yy_yyzz[i] = tg_0_yyzz[i] * fxi[i] + 2.0 * tg_y_yzz[i] * fxi[i] + tg_y_yyzz[i] * ra_y[i];

        tg_yy_yzzz[i] = tg_0_yzzz[i] * fxi[i] + tg_y_zzz[i] * fxi[i] + tg_y_yzzz[i] * ra_y[i];

        tg_yy_zzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_y_zzzz[i] * ra_y[i];

        tg_yz_xxxx[i] = tg_z_xxxx[i] * ra_y[i];

        tg_yz_xxxy[i] = tg_y_xxxy[i] * ra_z[i];

        tg_yz_xxxz[i] = tg_z_xxxz[i] * ra_y[i];

        tg_yz_xxyy[i] = tg_y_xxyy[i] * ra_z[i];

        tg_yz_xxyz[i] = tg_z_xxz[i] * fxi[i] + tg_z_xxyz[i] * ra_y[i];

        tg_yz_xxzz[i] = tg_z_xxzz[i] * ra_y[i];

        tg_yz_xyyy[i] = tg_y_xyyy[i] * ra_z[i];

        tg_yz_xyyz[i] = 2.0 * tg_z_xyz[i] * fxi[i] + tg_z_xyyz[i] * ra_y[i];

        tg_yz_xyzz[i] = tg_z_xzz[i] * fxi[i] + tg_z_xyzz[i] * ra_y[i];

        tg_yz_xzzz[i] = tg_z_xzzz[i] * ra_y[i];

        tg_yz_yyyy[i] = tg_y_yyyy[i] * ra_z[i];

        tg_yz_yyyz[i] = 3.0 * tg_z_yyz[i] * fxi[i] + tg_z_yyyz[i] * ra_y[i];

        tg_yz_yyzz[i] = 2.0 * tg_z_yzz[i] * fxi[i] + tg_z_yyzz[i] * ra_y[i];

        tg_yz_yzzz[i] = tg_z_zzz[i] * fxi[i] + tg_z_yzzz[i] * ra_y[i];

        tg_yz_zzzz[i] = tg_z_zzzz[i] * ra_y[i];

        tg_zz_xxxx[i] = tg_0_xxxx[i] * fxi[i] + tg_z_xxxx[i] * ra_z[i];

        tg_zz_xxxy[i] = tg_0_xxxy[i] * fxi[i] + tg_z_xxxy[i] * ra_z[i];

        tg_zz_xxxz[i] = tg_0_xxxz[i] * fxi[i] + tg_z_xxx[i] * fxi[i] + tg_z_xxxz[i] * ra_z[i];

        tg_zz_xxyy[i] = tg_0_xxyy[i] * fxi[i] + tg_z_xxyy[i] * ra_z[i];

        tg_zz_xxyz[i] = tg_0_xxyz[i] * fxi[i] + tg_z_xxy[i] * fxi[i] + tg_z_xxyz[i] * ra_z[i];

        tg_zz_xxzz[i] = tg_0_xxzz[i] * fxi[i] + 2.0 * tg_z_xxz[i] * fxi[i] + tg_z_xxzz[i] * ra_z[i];

        tg_zz_xyyy[i] = tg_0_xyyy[i] * fxi[i] + tg_z_xyyy[i] * ra_z[i];

        tg_zz_xyyz[i] = tg_0_xyyz[i] * fxi[i] + tg_z_xyy[i] * fxi[i] + tg_z_xyyz[i] * ra_z[i];

        tg_zz_xyzz[i] = tg_0_xyzz[i] * fxi[i] + 2.0 * tg_z_xyz[i] * fxi[i] + tg_z_xyzz[i] * ra_z[i];

        tg_zz_xzzz[i] = tg_0_xzzz[i] * fxi[i] + 3.0 * tg_z_xzz[i] * fxi[i] + tg_z_xzzz[i] * ra_z[i];

        tg_zz_yyyy[i] = tg_0_yyyy[i] * fxi[i] + tg_z_yyyy[i] * ra_z[i];

        tg_zz_yyyz[i] = tg_0_yyyz[i] * fxi[i] + tg_z_yyy[i] * fxi[i] + tg_z_yyyz[i] * ra_z[i];

        tg_zz_yyzz[i] = tg_0_yyzz[i] * fxi[i] + 2.0 * tg_z_yyz[i] * fxi[i] + tg_z_yyzz[i] * ra_z[i];

        tg_zz_yzzz[i] = tg_0_yzzz[i] * fxi[i] + 3.0 * tg_z_yzz[i] * fxi[i] + tg_z_yzzz[i] * ra_z[i];

        tg_zz_zzzz[i] = tg_0_zzzz[i] * fxi[i] + 4.0 * tg_z_zzz[i] * fxi[i] + tg_z_zzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

