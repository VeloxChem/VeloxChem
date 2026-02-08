#include "T2CHrrABRecFF.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_ff(CSimdArray<double>& cbuffer, 
            const size_t idx_ff,
            const size_t idx_df,
            const size_t idx_dg,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : DF

    auto t_xx_xxx = cbuffer.data(idx_df);

    auto t_xx_xxy = cbuffer.data(idx_df + 1);

    auto t_xx_xxz = cbuffer.data(idx_df + 2);

    auto t_xx_xyy = cbuffer.data(idx_df + 3);

    auto t_xx_xyz = cbuffer.data(idx_df + 4);

    auto t_xx_xzz = cbuffer.data(idx_df + 5);

    auto t_xx_yyy = cbuffer.data(idx_df + 6);

    auto t_xx_yyz = cbuffer.data(idx_df + 7);

    auto t_xx_yzz = cbuffer.data(idx_df + 8);

    auto t_xx_zzz = cbuffer.data(idx_df + 9);

    auto t_xy_xxx = cbuffer.data(idx_df + 10);

    auto t_xy_xxy = cbuffer.data(idx_df + 11);

    auto t_xy_xxz = cbuffer.data(idx_df + 12);

    auto t_xy_xyy = cbuffer.data(idx_df + 13);

    auto t_xy_xyz = cbuffer.data(idx_df + 14);

    auto t_xy_xzz = cbuffer.data(idx_df + 15);

    auto t_xy_yyy = cbuffer.data(idx_df + 16);

    auto t_xy_yyz = cbuffer.data(idx_df + 17);

    auto t_xy_yzz = cbuffer.data(idx_df + 18);

    auto t_xy_zzz = cbuffer.data(idx_df + 19);

    auto t_xz_xxx = cbuffer.data(idx_df + 20);

    auto t_xz_xxy = cbuffer.data(idx_df + 21);

    auto t_xz_xxz = cbuffer.data(idx_df + 22);

    auto t_xz_xyy = cbuffer.data(idx_df + 23);

    auto t_xz_xyz = cbuffer.data(idx_df + 24);

    auto t_xz_xzz = cbuffer.data(idx_df + 25);

    auto t_xz_yyy = cbuffer.data(idx_df + 26);

    auto t_xz_yyz = cbuffer.data(idx_df + 27);

    auto t_xz_yzz = cbuffer.data(idx_df + 28);

    auto t_xz_zzz = cbuffer.data(idx_df + 29);

    auto t_yy_xxx = cbuffer.data(idx_df + 30);

    auto t_yy_xxy = cbuffer.data(idx_df + 31);

    auto t_yy_xxz = cbuffer.data(idx_df + 32);

    auto t_yy_xyy = cbuffer.data(idx_df + 33);

    auto t_yy_xyz = cbuffer.data(idx_df + 34);

    auto t_yy_xzz = cbuffer.data(idx_df + 35);

    auto t_yy_yyy = cbuffer.data(idx_df + 36);

    auto t_yy_yyz = cbuffer.data(idx_df + 37);

    auto t_yy_yzz = cbuffer.data(idx_df + 38);

    auto t_yy_zzz = cbuffer.data(idx_df + 39);

    auto t_yz_xxx = cbuffer.data(idx_df + 40);

    auto t_yz_xxy = cbuffer.data(idx_df + 41);

    auto t_yz_xxz = cbuffer.data(idx_df + 42);

    auto t_yz_xyy = cbuffer.data(idx_df + 43);

    auto t_yz_xyz = cbuffer.data(idx_df + 44);

    auto t_yz_xzz = cbuffer.data(idx_df + 45);

    auto t_yz_yyy = cbuffer.data(idx_df + 46);

    auto t_yz_yyz = cbuffer.data(idx_df + 47);

    auto t_yz_yzz = cbuffer.data(idx_df + 48);

    auto t_yz_zzz = cbuffer.data(idx_df + 49);

    auto t_zz_xxx = cbuffer.data(idx_df + 50);

    auto t_zz_xxy = cbuffer.data(idx_df + 51);

    auto t_zz_xxz = cbuffer.data(idx_df + 52);

    auto t_zz_xyy = cbuffer.data(idx_df + 53);

    auto t_zz_xyz = cbuffer.data(idx_df + 54);

    auto t_zz_xzz = cbuffer.data(idx_df + 55);

    auto t_zz_yyy = cbuffer.data(idx_df + 56);

    auto t_zz_yyz = cbuffer.data(idx_df + 57);

    auto t_zz_yzz = cbuffer.data(idx_df + 58);

    auto t_zz_zzz = cbuffer.data(idx_df + 59);

    // Set up components of auxiliary buffer : DG

    auto t_xx_xxxx = cbuffer.data(idx_dg);

    auto t_xx_xxxy = cbuffer.data(idx_dg + 1);

    auto t_xx_xxxz = cbuffer.data(idx_dg + 2);

    auto t_xx_xxyy = cbuffer.data(idx_dg + 3);

    auto t_xx_xxyz = cbuffer.data(idx_dg + 4);

    auto t_xx_xxzz = cbuffer.data(idx_dg + 5);

    auto t_xx_xyyy = cbuffer.data(idx_dg + 6);

    auto t_xx_xyyz = cbuffer.data(idx_dg + 7);

    auto t_xx_xyzz = cbuffer.data(idx_dg + 8);

    auto t_xx_xzzz = cbuffer.data(idx_dg + 9);

    auto t_xy_xxxx = cbuffer.data(idx_dg + 15);

    auto t_xy_xxxy = cbuffer.data(idx_dg + 16);

    auto t_xy_xxxz = cbuffer.data(idx_dg + 17);

    auto t_xy_xxyy = cbuffer.data(idx_dg + 18);

    auto t_xy_xxyz = cbuffer.data(idx_dg + 19);

    auto t_xy_xxzz = cbuffer.data(idx_dg + 20);

    auto t_xy_xyyy = cbuffer.data(idx_dg + 21);

    auto t_xy_xyyz = cbuffer.data(idx_dg + 22);

    auto t_xy_xyzz = cbuffer.data(idx_dg + 23);

    auto t_xy_xzzz = cbuffer.data(idx_dg + 24);

    auto t_xz_xxxx = cbuffer.data(idx_dg + 30);

    auto t_xz_xxxy = cbuffer.data(idx_dg + 31);

    auto t_xz_xxxz = cbuffer.data(idx_dg + 32);

    auto t_xz_xxyy = cbuffer.data(idx_dg + 33);

    auto t_xz_xxyz = cbuffer.data(idx_dg + 34);

    auto t_xz_xxzz = cbuffer.data(idx_dg + 35);

    auto t_xz_xyyy = cbuffer.data(idx_dg + 36);

    auto t_xz_xyyz = cbuffer.data(idx_dg + 37);

    auto t_xz_xyzz = cbuffer.data(idx_dg + 38);

    auto t_xz_xzzz = cbuffer.data(idx_dg + 39);

    auto t_yy_xxxx = cbuffer.data(idx_dg + 45);

    auto t_yy_xxxy = cbuffer.data(idx_dg + 46);

    auto t_yy_xxxz = cbuffer.data(idx_dg + 47);

    auto t_yy_xxyy = cbuffer.data(idx_dg + 48);

    auto t_yy_xxyz = cbuffer.data(idx_dg + 49);

    auto t_yy_xxzz = cbuffer.data(idx_dg + 50);

    auto t_yy_xyyy = cbuffer.data(idx_dg + 51);

    auto t_yy_xyyz = cbuffer.data(idx_dg + 52);

    auto t_yy_xyzz = cbuffer.data(idx_dg + 53);

    auto t_yy_xzzz = cbuffer.data(idx_dg + 54);

    auto t_yy_yyyy = cbuffer.data(idx_dg + 55);

    auto t_yy_yyyz = cbuffer.data(idx_dg + 56);

    auto t_yy_yyzz = cbuffer.data(idx_dg + 57);

    auto t_yy_yzzz = cbuffer.data(idx_dg + 58);

    auto t_yz_xxxx = cbuffer.data(idx_dg + 60);

    auto t_yz_xxxy = cbuffer.data(idx_dg + 61);

    auto t_yz_xxxz = cbuffer.data(idx_dg + 62);

    auto t_yz_xxyy = cbuffer.data(idx_dg + 63);

    auto t_yz_xxyz = cbuffer.data(idx_dg + 64);

    auto t_yz_xxzz = cbuffer.data(idx_dg + 65);

    auto t_yz_xyyy = cbuffer.data(idx_dg + 66);

    auto t_yz_xyyz = cbuffer.data(idx_dg + 67);

    auto t_yz_xyzz = cbuffer.data(idx_dg + 68);

    auto t_yz_xzzz = cbuffer.data(idx_dg + 69);

    auto t_yz_yyyy = cbuffer.data(idx_dg + 70);

    auto t_yz_yyyz = cbuffer.data(idx_dg + 71);

    auto t_yz_yyzz = cbuffer.data(idx_dg + 72);

    auto t_yz_yzzz = cbuffer.data(idx_dg + 73);

    auto t_zz_xxxx = cbuffer.data(idx_dg + 75);

    auto t_zz_xxxy = cbuffer.data(idx_dg + 76);

    auto t_zz_xxxz = cbuffer.data(idx_dg + 77);

    auto t_zz_xxyy = cbuffer.data(idx_dg + 78);

    auto t_zz_xxyz = cbuffer.data(idx_dg + 79);

    auto t_zz_xxzz = cbuffer.data(idx_dg + 80);

    auto t_zz_xyyy = cbuffer.data(idx_dg + 81);

    auto t_zz_xyyz = cbuffer.data(idx_dg + 82);

    auto t_zz_xyzz = cbuffer.data(idx_dg + 83);

    auto t_zz_xzzz = cbuffer.data(idx_dg + 84);

    auto t_zz_yyyy = cbuffer.data(idx_dg + 85);

    auto t_zz_yyyz = cbuffer.data(idx_dg + 86);

    auto t_zz_yyzz = cbuffer.data(idx_dg + 87);

    auto t_zz_yzzz = cbuffer.data(idx_dg + 88);

    auto t_zz_zzzz = cbuffer.data(idx_dg + 89);

    // Set up components of targeted buffer : FF

    auto t_xxx_xxx = cbuffer.data(idx_ff);

    auto t_xxx_xxy = cbuffer.data(idx_ff + 1);

    auto t_xxx_xxz = cbuffer.data(idx_ff + 2);

    auto t_xxx_xyy = cbuffer.data(idx_ff + 3);

    auto t_xxx_xyz = cbuffer.data(idx_ff + 4);

    auto t_xxx_xzz = cbuffer.data(idx_ff + 5);

    auto t_xxx_yyy = cbuffer.data(idx_ff + 6);

    auto t_xxx_yyz = cbuffer.data(idx_ff + 7);

    auto t_xxx_yzz = cbuffer.data(idx_ff + 8);

    auto t_xxx_zzz = cbuffer.data(idx_ff + 9);

    auto t_xxy_xxx = cbuffer.data(idx_ff + 10);

    auto t_xxy_xxy = cbuffer.data(idx_ff + 11);

    auto t_xxy_xxz = cbuffer.data(idx_ff + 12);

    auto t_xxy_xyy = cbuffer.data(idx_ff + 13);

    auto t_xxy_xyz = cbuffer.data(idx_ff + 14);

    auto t_xxy_xzz = cbuffer.data(idx_ff + 15);

    auto t_xxy_yyy = cbuffer.data(idx_ff + 16);

    auto t_xxy_yyz = cbuffer.data(idx_ff + 17);

    auto t_xxy_yzz = cbuffer.data(idx_ff + 18);

    auto t_xxy_zzz = cbuffer.data(idx_ff + 19);

    auto t_xxz_xxx = cbuffer.data(idx_ff + 20);

    auto t_xxz_xxy = cbuffer.data(idx_ff + 21);

    auto t_xxz_xxz = cbuffer.data(idx_ff + 22);

    auto t_xxz_xyy = cbuffer.data(idx_ff + 23);

    auto t_xxz_xyz = cbuffer.data(idx_ff + 24);

    auto t_xxz_xzz = cbuffer.data(idx_ff + 25);

    auto t_xxz_yyy = cbuffer.data(idx_ff + 26);

    auto t_xxz_yyz = cbuffer.data(idx_ff + 27);

    auto t_xxz_yzz = cbuffer.data(idx_ff + 28);

    auto t_xxz_zzz = cbuffer.data(idx_ff + 29);

    auto t_xyy_xxx = cbuffer.data(idx_ff + 30);

    auto t_xyy_xxy = cbuffer.data(idx_ff + 31);

    auto t_xyy_xxz = cbuffer.data(idx_ff + 32);

    auto t_xyy_xyy = cbuffer.data(idx_ff + 33);

    auto t_xyy_xyz = cbuffer.data(idx_ff + 34);

    auto t_xyy_xzz = cbuffer.data(idx_ff + 35);

    auto t_xyy_yyy = cbuffer.data(idx_ff + 36);

    auto t_xyy_yyz = cbuffer.data(idx_ff + 37);

    auto t_xyy_yzz = cbuffer.data(idx_ff + 38);

    auto t_xyy_zzz = cbuffer.data(idx_ff + 39);

    auto t_xyz_xxx = cbuffer.data(idx_ff + 40);

    auto t_xyz_xxy = cbuffer.data(idx_ff + 41);

    auto t_xyz_xxz = cbuffer.data(idx_ff + 42);

    auto t_xyz_xyy = cbuffer.data(idx_ff + 43);

    auto t_xyz_xyz = cbuffer.data(idx_ff + 44);

    auto t_xyz_xzz = cbuffer.data(idx_ff + 45);

    auto t_xyz_yyy = cbuffer.data(idx_ff + 46);

    auto t_xyz_yyz = cbuffer.data(idx_ff + 47);

    auto t_xyz_yzz = cbuffer.data(idx_ff + 48);

    auto t_xyz_zzz = cbuffer.data(idx_ff + 49);

    auto t_xzz_xxx = cbuffer.data(idx_ff + 50);

    auto t_xzz_xxy = cbuffer.data(idx_ff + 51);

    auto t_xzz_xxz = cbuffer.data(idx_ff + 52);

    auto t_xzz_xyy = cbuffer.data(idx_ff + 53);

    auto t_xzz_xyz = cbuffer.data(idx_ff + 54);

    auto t_xzz_xzz = cbuffer.data(idx_ff + 55);

    auto t_xzz_yyy = cbuffer.data(idx_ff + 56);

    auto t_xzz_yyz = cbuffer.data(idx_ff + 57);

    auto t_xzz_yzz = cbuffer.data(idx_ff + 58);

    auto t_xzz_zzz = cbuffer.data(idx_ff + 59);

    auto t_yyy_xxx = cbuffer.data(idx_ff + 60);

    auto t_yyy_xxy = cbuffer.data(idx_ff + 61);

    auto t_yyy_xxz = cbuffer.data(idx_ff + 62);

    auto t_yyy_xyy = cbuffer.data(idx_ff + 63);

    auto t_yyy_xyz = cbuffer.data(idx_ff + 64);

    auto t_yyy_xzz = cbuffer.data(idx_ff + 65);

    auto t_yyy_yyy = cbuffer.data(idx_ff + 66);

    auto t_yyy_yyz = cbuffer.data(idx_ff + 67);

    auto t_yyy_yzz = cbuffer.data(idx_ff + 68);

    auto t_yyy_zzz = cbuffer.data(idx_ff + 69);

    auto t_yyz_xxx = cbuffer.data(idx_ff + 70);

    auto t_yyz_xxy = cbuffer.data(idx_ff + 71);

    auto t_yyz_xxz = cbuffer.data(idx_ff + 72);

    auto t_yyz_xyy = cbuffer.data(idx_ff + 73);

    auto t_yyz_xyz = cbuffer.data(idx_ff + 74);

    auto t_yyz_xzz = cbuffer.data(idx_ff + 75);

    auto t_yyz_yyy = cbuffer.data(idx_ff + 76);

    auto t_yyz_yyz = cbuffer.data(idx_ff + 77);

    auto t_yyz_yzz = cbuffer.data(idx_ff + 78);

    auto t_yyz_zzz = cbuffer.data(idx_ff + 79);

    auto t_yzz_xxx = cbuffer.data(idx_ff + 80);

    auto t_yzz_xxy = cbuffer.data(idx_ff + 81);

    auto t_yzz_xxz = cbuffer.data(idx_ff + 82);

    auto t_yzz_xyy = cbuffer.data(idx_ff + 83);

    auto t_yzz_xyz = cbuffer.data(idx_ff + 84);

    auto t_yzz_xzz = cbuffer.data(idx_ff + 85);

    auto t_yzz_yyy = cbuffer.data(idx_ff + 86);

    auto t_yzz_yyz = cbuffer.data(idx_ff + 87);

    auto t_yzz_yzz = cbuffer.data(idx_ff + 88);

    auto t_yzz_zzz = cbuffer.data(idx_ff + 89);

    auto t_zzz_xxx = cbuffer.data(idx_ff + 90);

    auto t_zzz_xxy = cbuffer.data(idx_ff + 91);

    auto t_zzz_xxz = cbuffer.data(idx_ff + 92);

    auto t_zzz_xyy = cbuffer.data(idx_ff + 93);

    auto t_zzz_xyz = cbuffer.data(idx_ff + 94);

    auto t_zzz_xzz = cbuffer.data(idx_ff + 95);

    auto t_zzz_yyy = cbuffer.data(idx_ff + 96);

    auto t_zzz_yyz = cbuffer.data(idx_ff + 97);

    auto t_zzz_yzz = cbuffer.data(idx_ff + 98);

    auto t_zzz_zzz = cbuffer.data(idx_ff + 99);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xx_xxx, t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxy, t_xx_xxyy, t_xx_xxyz, t_xx_xxz, t_xx_xxzz, t_xx_xyy, t_xx_xyyy, t_xx_xyyz, t_xx_xyz, t_xx_xyzz, t_xx_xzz, t_xx_xzzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz, t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, t_xxx_xyy, t_xxx_xyz, t_xxx_xzz, t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz, t_xxy_xxx, t_xxy_xxy, t_xxy_xxz, t_xxy_xyy, t_xxy_xyz, t_xxy_xzz, t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, t_xxy_zzz, t_xxz_xxx, t_xxz_xxy, t_xxz_xxz, t_xxz_xyy, t_xxz_xyz, t_xxz_xzz, t_xxz_yyy, t_xxz_yyz, t_xxz_yzz, t_xxz_zzz, t_xy_xxx, t_xy_xxxx, t_xy_xxxy, t_xy_xxxz, t_xy_xxy, t_xy_xxyy, t_xy_xxyz, t_xy_xxz, t_xy_xxzz, t_xy_xyy, t_xy_xyyy, t_xy_xyyz, t_xy_xyz, t_xy_xyzz, t_xy_xzz, t_xy_xzzz, t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz, t_xyy_xxx, t_xyy_xxy, t_xyy_xxz, t_xyy_xyy, t_xyy_xyz, t_xyy_xzz, t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, t_xyy_zzz, t_xyz_xxx, t_xyz_xxy, t_xyz_xxz, t_xyz_xyy, t_xyz_xyz, t_xyz_xzz, t_xyz_yyy, t_xyz_yyz, t_xyz_yzz, t_xyz_zzz, t_xz_xxx, t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxy, t_xz_xxyy, t_xz_xxyz, t_xz_xxz, t_xz_xxzz, t_xz_xyy, t_xz_xyyy, t_xz_xyyz, t_xz_xyz, t_xz_xyzz, t_xz_xzz, t_xz_xzzz, t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz, t_xzz_xxx, t_xzz_xxy, t_xzz_xxz, t_xzz_xyy, t_xzz_xyz, t_xzz_xzz, t_xzz_yyy, t_xzz_yyz, t_xzz_yzz, t_xzz_zzz, t_yy_xxx, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz, t_yy_xxy, t_yy_xxyy, t_yy_xxyz, t_yy_xxz, t_yy_xxzz, t_yy_xyy, t_yy_xyyy, t_yy_xyyz, t_yy_xyz, t_yy_xyzz, t_yy_xzz, t_yy_xzzz, t_yy_yyy, t_yy_yyyy, t_yy_yyyz, t_yy_yyz, t_yy_yyzz, t_yy_yzz, t_yy_yzzz, t_yy_zzz, t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, t_yyy_xyy, t_yyy_xyz, t_yyy_xzz, t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, t_yyy_zzz, t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy, t_yyz_xyz, t_yyz_xzz, t_yyz_yyy, t_yyz_yyz, t_yyz_yzz, t_yyz_zzz, t_yz_xxx, t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, t_yz_xxy, t_yz_xxyy, t_yz_xxyz, t_yz_xxz, t_yz_xxzz, t_yz_xyy, t_yz_xyyy, t_yz_xyyz, t_yz_xyz, t_yz_xyzz, t_yz_xzz, t_yz_xzzz, t_yz_yyy, t_yz_yyyy, t_yz_yyyz, t_yz_yyz, t_yz_yyzz, t_yz_yzz, t_yz_yzzz, t_yz_zzz, t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy, t_yzz_xyz, t_yzz_xzz, t_yzz_yyy, t_yzz_yyz, t_yzz_yzz, t_yzz_zzz, t_zz_xxx, t_zz_xxxx, t_zz_xxxy, t_zz_xxxz, t_zz_xxy, t_zz_xxyy, t_zz_xxyz, t_zz_xxz, t_zz_xxzz, t_zz_xyy, t_zz_xyyy, t_zz_xyyz, t_zz_xyz, t_zz_xyzz, t_zz_xzz, t_zz_xzzz, t_zz_yyy, t_zz_yyyy, t_zz_yyyz, t_zz_yyz, t_zz_yyzz, t_zz_yzz, t_zz_yzzz, t_zz_zzz, t_zz_zzzz, t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, t_zzz_xyy, t_zzz_xyz, t_zzz_xzz, t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, t_zzz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxx_xxx[i] = -t_xx_xxx[i] * ab_x[i] + t_xx_xxxx[i];

        t_xxx_xxy[i] = -t_xx_xxy[i] * ab_x[i] + t_xx_xxxy[i];

        t_xxx_xxz[i] = -t_xx_xxz[i] * ab_x[i] + t_xx_xxxz[i];

        t_xxx_xyy[i] = -t_xx_xyy[i] * ab_x[i] + t_xx_xxyy[i];

        t_xxx_xyz[i] = -t_xx_xyz[i] * ab_x[i] + t_xx_xxyz[i];

        t_xxx_xzz[i] = -t_xx_xzz[i] * ab_x[i] + t_xx_xxzz[i];

        t_xxx_yyy[i] = -t_xx_yyy[i] * ab_x[i] + t_xx_xyyy[i];

        t_xxx_yyz[i] = -t_xx_yyz[i] * ab_x[i] + t_xx_xyyz[i];

        t_xxx_yzz[i] = -t_xx_yzz[i] * ab_x[i] + t_xx_xyzz[i];

        t_xxx_zzz[i] = -t_xx_zzz[i] * ab_x[i] + t_xx_xzzz[i];

        t_xxy_xxx[i] = -t_xy_xxx[i] * ab_x[i] + t_xy_xxxx[i];

        t_xxy_xxy[i] = -t_xy_xxy[i] * ab_x[i] + t_xy_xxxy[i];

        t_xxy_xxz[i] = -t_xy_xxz[i] * ab_x[i] + t_xy_xxxz[i];

        t_xxy_xyy[i] = -t_xy_xyy[i] * ab_x[i] + t_xy_xxyy[i];

        t_xxy_xyz[i] = -t_xy_xyz[i] * ab_x[i] + t_xy_xxyz[i];

        t_xxy_xzz[i] = -t_xy_xzz[i] * ab_x[i] + t_xy_xxzz[i];

        t_xxy_yyy[i] = -t_xy_yyy[i] * ab_x[i] + t_xy_xyyy[i];

        t_xxy_yyz[i] = -t_xy_yyz[i] * ab_x[i] + t_xy_xyyz[i];

        t_xxy_yzz[i] = -t_xy_yzz[i] * ab_x[i] + t_xy_xyzz[i];

        t_xxy_zzz[i] = -t_xy_zzz[i] * ab_x[i] + t_xy_xzzz[i];

        t_xxz_xxx[i] = -t_xz_xxx[i] * ab_x[i] + t_xz_xxxx[i];

        t_xxz_xxy[i] = -t_xz_xxy[i] * ab_x[i] + t_xz_xxxy[i];

        t_xxz_xxz[i] = -t_xz_xxz[i] * ab_x[i] + t_xz_xxxz[i];

        t_xxz_xyy[i] = -t_xz_xyy[i] * ab_x[i] + t_xz_xxyy[i];

        t_xxz_xyz[i] = -t_xz_xyz[i] * ab_x[i] + t_xz_xxyz[i];

        t_xxz_xzz[i] = -t_xz_xzz[i] * ab_x[i] + t_xz_xxzz[i];

        t_xxz_yyy[i] = -t_xz_yyy[i] * ab_x[i] + t_xz_xyyy[i];

        t_xxz_yyz[i] = -t_xz_yyz[i] * ab_x[i] + t_xz_xyyz[i];

        t_xxz_yzz[i] = -t_xz_yzz[i] * ab_x[i] + t_xz_xyzz[i];

        t_xxz_zzz[i] = -t_xz_zzz[i] * ab_x[i] + t_xz_xzzz[i];

        t_xyy_xxx[i] = -t_yy_xxx[i] * ab_x[i] + t_yy_xxxx[i];

        t_xyy_xxy[i] = -t_yy_xxy[i] * ab_x[i] + t_yy_xxxy[i];

        t_xyy_xxz[i] = -t_yy_xxz[i] * ab_x[i] + t_yy_xxxz[i];

        t_xyy_xyy[i] = -t_yy_xyy[i] * ab_x[i] + t_yy_xxyy[i];

        t_xyy_xyz[i] = -t_yy_xyz[i] * ab_x[i] + t_yy_xxyz[i];

        t_xyy_xzz[i] = -t_yy_xzz[i] * ab_x[i] + t_yy_xxzz[i];

        t_xyy_yyy[i] = -t_yy_yyy[i] * ab_x[i] + t_yy_xyyy[i];

        t_xyy_yyz[i] = -t_yy_yyz[i] * ab_x[i] + t_yy_xyyz[i];

        t_xyy_yzz[i] = -t_yy_yzz[i] * ab_x[i] + t_yy_xyzz[i];

        t_xyy_zzz[i] = -t_yy_zzz[i] * ab_x[i] + t_yy_xzzz[i];

        t_xyz_xxx[i] = -t_yz_xxx[i] * ab_x[i] + t_yz_xxxx[i];

        t_xyz_xxy[i] = -t_yz_xxy[i] * ab_x[i] + t_yz_xxxy[i];

        t_xyz_xxz[i] = -t_yz_xxz[i] * ab_x[i] + t_yz_xxxz[i];

        t_xyz_xyy[i] = -t_yz_xyy[i] * ab_x[i] + t_yz_xxyy[i];

        t_xyz_xyz[i] = -t_yz_xyz[i] * ab_x[i] + t_yz_xxyz[i];

        t_xyz_xzz[i] = -t_yz_xzz[i] * ab_x[i] + t_yz_xxzz[i];

        t_xyz_yyy[i] = -t_yz_yyy[i] * ab_x[i] + t_yz_xyyy[i];

        t_xyz_yyz[i] = -t_yz_yyz[i] * ab_x[i] + t_yz_xyyz[i];

        t_xyz_yzz[i] = -t_yz_yzz[i] * ab_x[i] + t_yz_xyzz[i];

        t_xyz_zzz[i] = -t_yz_zzz[i] * ab_x[i] + t_yz_xzzz[i];

        t_xzz_xxx[i] = -t_zz_xxx[i] * ab_x[i] + t_zz_xxxx[i];

        t_xzz_xxy[i] = -t_zz_xxy[i] * ab_x[i] + t_zz_xxxy[i];

        t_xzz_xxz[i] = -t_zz_xxz[i] * ab_x[i] + t_zz_xxxz[i];

        t_xzz_xyy[i] = -t_zz_xyy[i] * ab_x[i] + t_zz_xxyy[i];

        t_xzz_xyz[i] = -t_zz_xyz[i] * ab_x[i] + t_zz_xxyz[i];

        t_xzz_xzz[i] = -t_zz_xzz[i] * ab_x[i] + t_zz_xxzz[i];

        t_xzz_yyy[i] = -t_zz_yyy[i] * ab_x[i] + t_zz_xyyy[i];

        t_xzz_yyz[i] = -t_zz_yyz[i] * ab_x[i] + t_zz_xyyz[i];

        t_xzz_yzz[i] = -t_zz_yzz[i] * ab_x[i] + t_zz_xyzz[i];

        t_xzz_zzz[i] = -t_zz_zzz[i] * ab_x[i] + t_zz_xzzz[i];

        t_yyy_xxx[i] = -t_yy_xxx[i] * ab_y[i] + t_yy_xxxy[i];

        t_yyy_xxy[i] = -t_yy_xxy[i] * ab_y[i] + t_yy_xxyy[i];

        t_yyy_xxz[i] = -t_yy_xxz[i] * ab_y[i] + t_yy_xxyz[i];

        t_yyy_xyy[i] = -t_yy_xyy[i] * ab_y[i] + t_yy_xyyy[i];

        t_yyy_xyz[i] = -t_yy_xyz[i] * ab_y[i] + t_yy_xyyz[i];

        t_yyy_xzz[i] = -t_yy_xzz[i] * ab_y[i] + t_yy_xyzz[i];

        t_yyy_yyy[i] = -t_yy_yyy[i] * ab_y[i] + t_yy_yyyy[i];

        t_yyy_yyz[i] = -t_yy_yyz[i] * ab_y[i] + t_yy_yyyz[i];

        t_yyy_yzz[i] = -t_yy_yzz[i] * ab_y[i] + t_yy_yyzz[i];

        t_yyy_zzz[i] = -t_yy_zzz[i] * ab_y[i] + t_yy_yzzz[i];

        t_yyz_xxx[i] = -t_yz_xxx[i] * ab_y[i] + t_yz_xxxy[i];

        t_yyz_xxy[i] = -t_yz_xxy[i] * ab_y[i] + t_yz_xxyy[i];

        t_yyz_xxz[i] = -t_yz_xxz[i] * ab_y[i] + t_yz_xxyz[i];

        t_yyz_xyy[i] = -t_yz_xyy[i] * ab_y[i] + t_yz_xyyy[i];

        t_yyz_xyz[i] = -t_yz_xyz[i] * ab_y[i] + t_yz_xyyz[i];

        t_yyz_xzz[i] = -t_yz_xzz[i] * ab_y[i] + t_yz_xyzz[i];

        t_yyz_yyy[i] = -t_yz_yyy[i] * ab_y[i] + t_yz_yyyy[i];

        t_yyz_yyz[i] = -t_yz_yyz[i] * ab_y[i] + t_yz_yyyz[i];

        t_yyz_yzz[i] = -t_yz_yzz[i] * ab_y[i] + t_yz_yyzz[i];

        t_yyz_zzz[i] = -t_yz_zzz[i] * ab_y[i] + t_yz_yzzz[i];

        t_yzz_xxx[i] = -t_zz_xxx[i] * ab_y[i] + t_zz_xxxy[i];

        t_yzz_xxy[i] = -t_zz_xxy[i] * ab_y[i] + t_zz_xxyy[i];

        t_yzz_xxz[i] = -t_zz_xxz[i] * ab_y[i] + t_zz_xxyz[i];

        t_yzz_xyy[i] = -t_zz_xyy[i] * ab_y[i] + t_zz_xyyy[i];

        t_yzz_xyz[i] = -t_zz_xyz[i] * ab_y[i] + t_zz_xyyz[i];

        t_yzz_xzz[i] = -t_zz_xzz[i] * ab_y[i] + t_zz_xyzz[i];

        t_yzz_yyy[i] = -t_zz_yyy[i] * ab_y[i] + t_zz_yyyy[i];

        t_yzz_yyz[i] = -t_zz_yyz[i] * ab_y[i] + t_zz_yyyz[i];

        t_yzz_yzz[i] = -t_zz_yzz[i] * ab_y[i] + t_zz_yyzz[i];

        t_yzz_zzz[i] = -t_zz_zzz[i] * ab_y[i] + t_zz_yzzz[i];

        t_zzz_xxx[i] = -t_zz_xxx[i] * ab_z[i] + t_zz_xxxz[i];

        t_zzz_xxy[i] = -t_zz_xxy[i] * ab_z[i] + t_zz_xxyz[i];

        t_zzz_xxz[i] = -t_zz_xxz[i] * ab_z[i] + t_zz_xxzz[i];

        t_zzz_xyy[i] = -t_zz_xyy[i] * ab_z[i] + t_zz_xyyz[i];

        t_zzz_xyz[i] = -t_zz_xyz[i] * ab_z[i] + t_zz_xyzz[i];

        t_zzz_xzz[i] = -t_zz_xzz[i] * ab_z[i] + t_zz_xzzz[i];

        t_zzz_yyy[i] = -t_zz_yyy[i] * ab_z[i] + t_zz_yyyz[i];

        t_zzz_yyz[i] = -t_zz_yyz[i] * ab_z[i] + t_zz_yyzz[i];

        t_zzz_yzz[i] = -t_zz_yzz[i] * ab_z[i] + t_zz_yzzz[i];

        t_zzz_zzz[i] = -t_zz_zzz[i] * ab_z[i] + t_zz_zzzz[i];
    }
}

} // t2chrr namespace

