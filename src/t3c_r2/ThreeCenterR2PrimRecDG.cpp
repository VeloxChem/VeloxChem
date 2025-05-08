#include "ThreeCenterR2PrimRecDG.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_dg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_dg,
                const size_t idx_sg,
                const size_t idx_pf,
                const size_t idx_pg,
                const size_t idx_dd,
                const size_t idx_df,
                const size_t idx_dg,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_sg);

    auto ts_0_xxxy = pbuffer.data(idx_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_pf);

    auto ts_x_xxy = pbuffer.data(idx_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_pf + 9);

    auto ts_y_xxx = pbuffer.data(idx_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_pf + 19);

    auto ts_z_xxx = pbuffer.data(idx_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_pf + 29);

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_pg);

    auto ts_x_xxxy = pbuffer.data(idx_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_dd);

    auto ts_xx_xy = pbuffer.data(idx_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_dd + 5);

    auto ts_xy_xx = pbuffer.data(idx_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_dd + 11);

    auto ts_xz_xx = pbuffer.data(idx_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_dd + 23);

    auto ts_yz_xx = pbuffer.data(idx_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_df);

    auto ts_xx_xxy = pbuffer.data(idx_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_df + 9);

    auto ts_xy_xxx = pbuffer.data(idx_df + 10);

    auto ts_xy_xxy = pbuffer.data(idx_df + 11);

    auto ts_xy_xxz = pbuffer.data(idx_df + 12);

    auto ts_xy_xyy = pbuffer.data(idx_df + 13);

    auto ts_xy_xyz = pbuffer.data(idx_df + 14);

    auto ts_xy_xzz = pbuffer.data(idx_df + 15);

    auto ts_xy_yyy = pbuffer.data(idx_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_df + 18);

    auto ts_xy_zzz = pbuffer.data(idx_df + 19);

    auto ts_xz_xxx = pbuffer.data(idx_df + 20);

    auto ts_xz_xxy = pbuffer.data(idx_df + 21);

    auto ts_xz_xxz = pbuffer.data(idx_df + 22);

    auto ts_xz_xyy = pbuffer.data(idx_df + 23);

    auto ts_xz_xyz = pbuffer.data(idx_df + 24);

    auto ts_xz_xzz = pbuffer.data(idx_df + 25);

    auto ts_xz_yyy = pbuffer.data(idx_df + 26);

    auto ts_xz_yyz = pbuffer.data(idx_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_df + 39);

    auto ts_yz_xxx = pbuffer.data(idx_df + 40);

    auto ts_yz_xxy = pbuffer.data(idx_df + 41);

    auto ts_yz_xxz = pbuffer.data(idx_df + 42);

    auto ts_yz_xyy = pbuffer.data(idx_df + 43);

    auto ts_yz_xyz = pbuffer.data(idx_df + 44);

    auto ts_yz_xzz = pbuffer.data(idx_df + 45);

    auto ts_yz_yyy = pbuffer.data(idx_df + 46);

    auto ts_yz_yyz = pbuffer.data(idx_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto ts_xy_xxxx = pbuffer.data(idx_dg + 15);

    auto ts_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto ts_xy_xxxz = pbuffer.data(idx_dg + 17);

    auto ts_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto ts_xy_xxyz = pbuffer.data(idx_dg + 19);

    auto ts_xy_xxzz = pbuffer.data(idx_dg + 20);

    auto ts_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto ts_xy_xyyz = pbuffer.data(idx_dg + 22);

    auto ts_xy_xyzz = pbuffer.data(idx_dg + 23);

    auto ts_xy_xzzz = pbuffer.data(idx_dg + 24);

    auto ts_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto ts_xy_zzzz = pbuffer.data(idx_dg + 29);

    auto ts_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto ts_xz_xxxy = pbuffer.data(idx_dg + 31);

    auto ts_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto ts_xz_xxyy = pbuffer.data(idx_dg + 33);

    auto ts_xz_xxyz = pbuffer.data(idx_dg + 34);

    auto ts_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto ts_xz_xyyy = pbuffer.data(idx_dg + 36);

    auto ts_xz_xyyz = pbuffer.data(idx_dg + 37);

    auto ts_xz_xyzz = pbuffer.data(idx_dg + 38);

    auto ts_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto ts_xz_yyyy = pbuffer.data(idx_dg + 40);

    auto ts_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto ts_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto ts_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto ts_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto ts_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto ts_yz_xxxx = pbuffer.data(idx_dg + 60);

    auto ts_yz_xxxy = pbuffer.data(idx_dg + 61);

    auto ts_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto ts_yz_xxyy = pbuffer.data(idx_dg + 63);

    auto ts_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto ts_yz_xyyy = pbuffer.data(idx_dg + 66);

    auto ts_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto ts_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto ts_yz_yyyy = pbuffer.data(idx_dg + 70);

    auto ts_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto ts_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto ts_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up 0-15 components of targeted buffer : DG

    auto gr_xx_xxxx = pbuffer.data(idx_g_dg);

    auto gr_xx_xxxy = pbuffer.data(idx_g_dg + 1);

    auto gr_xx_xxxz = pbuffer.data(idx_g_dg + 2);

    auto gr_xx_xxyy = pbuffer.data(idx_g_dg + 3);

    auto gr_xx_xxyz = pbuffer.data(idx_g_dg + 4);

    auto gr_xx_xxzz = pbuffer.data(idx_g_dg + 5);

    auto gr_xx_xyyy = pbuffer.data(idx_g_dg + 6);

    auto gr_xx_xyyz = pbuffer.data(idx_g_dg + 7);

    auto gr_xx_xyzz = pbuffer.data(idx_g_dg + 8);

    auto gr_xx_xzzz = pbuffer.data(idx_g_dg + 9);

    auto gr_xx_yyyy = pbuffer.data(idx_g_dg + 10);

    auto gr_xx_yyyz = pbuffer.data(idx_g_dg + 11);

    auto gr_xx_yyzz = pbuffer.data(idx_g_dg + 12);

    auto gr_xx_yzzz = pbuffer.data(idx_g_dg + 13);

    auto gr_xx_zzzz = pbuffer.data(idx_g_dg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxxx, gr_xx_xxxy, gr_xx_xxxz, gr_xx_xxyy, gr_xx_xxyz, gr_xx_xxzz, gr_xx_xyyy, gr_xx_xyyz, gr_xx_xyzz, gr_xx_xzzz, gr_xx_yyyy, gr_xx_yyyz, gr_xx_yyzz, gr_xx_yzzz, gr_xx_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, ts_xx_xx, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xy, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yy, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xx_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 + 16.0 * ts_x_xxx[i] * gfe_0 + 4.0 * ts_x_xxxx[i] * gfe_0 * gc_x[i] + 12.0 * ts_xx_xx[i] * gfe_0 + 8.0 * ts_xx_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xx_xxxx[i] * gfe_0 + ts_xx_xxxx[i] * rgc2_0;

        gr_xx_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 + 12.0 * ts_x_xxy[i] * gfe_0 + 4.0 * ts_x_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xy[i] * gfe_0 + 6.0 * ts_xx_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xxxy[i] * gfe_0 + ts_xx_xxxy[i] * rgc2_0;

        gr_xx_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 + 12.0 * ts_x_xxz[i] * gfe_0 + 4.0 * ts_x_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xz[i] * gfe_0 + 6.0 * ts_xx_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xxxz[i] * gfe_0 + ts_xx_xxxz[i] * rgc2_0;

        gr_xx_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 + 8.0 * ts_x_xyy[i] * gfe_0 + 4.0 * ts_x_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe_0 + 4.0 * ts_xx_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 + 4.0 * ts_xx_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xxyy[i] * gfe_0 + ts_xx_xxyy[i] * rgc2_0;

        gr_xx_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 + 8.0 * ts_x_xyz[i] * gfe_0 + 4.0 * ts_x_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yz[i] * gfe_0 + 4.0 * ts_xx_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xxyz[i] * gfe_0 + ts_xx_xxyz[i] * rgc2_0;

        gr_xx_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 + 8.0 * ts_x_xzz[i] * gfe_0 + 4.0 * ts_x_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 + 4.0 * ts_xx_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 + 4.0 * ts_xx_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xxzz[i] * gfe_0 + ts_xx_xxzz[i] * rgc2_0;

        gr_xx_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 + 4.0 * ts_x_yyy[i] * gfe_0 + 4.0 * ts_x_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xy[i] * gfe_0 + 6.0 * ts_xx_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xyyy[i] * gfe_0 + ts_xx_xyyy[i] * rgc2_0;

        gr_xx_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 + 4.0 * ts_x_yyz[i] * gfe_0 + 4.0 * ts_x_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe_0 + 4.0 * ts_xx_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xyyz[i] * gfe_0 + ts_xx_xyyz[i] * rgc2_0;

        gr_xx_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 + 4.0 * ts_x_yzz[i] * gfe_0 + 4.0 * ts_x_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_xy[i] * gfe_0 + 4.0 * ts_xx_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xyzz[i] * gfe_0 + ts_xx_xyzz[i] * rgc2_0;

        gr_xx_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 + 4.0 * ts_x_zzz[i] * gfe_0 + 4.0 * ts_x_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_xz[i] * gfe_0 + 6.0 * ts_xx_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xzzz[i] * gfe_0 + ts_xx_xzzz[i] * rgc2_0;

        gr_xx_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 + 4.0 * ts_x_yyyy[i] * gfe_0 * gc_x[i] + 12.0 * ts_xx_yy[i] * gfe_0 + 8.0 * ts_xx_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_yyyy[i] * gfe_0 + ts_xx_yyyy[i] * rgc2_0;

        gr_xx_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 + 4.0 * ts_x_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xx_yz[i] * gfe_0 + 6.0 * ts_xx_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yyyz[i] * gfe_0 + ts_xx_yyyz[i] * rgc2_0;

        gr_xx_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 + 4.0 * ts_x_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 + 4.0 * ts_xx_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_yy[i] * gfe_0 + 4.0 * ts_xx_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yyzz[i] * gfe_0 + ts_xx_yyzz[i] * rgc2_0;

        gr_xx_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 + 4.0 * ts_x_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xx_yz[i] * gfe_0 + 6.0 * ts_xx_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yzzz[i] * gfe_0 + ts_xx_yzzz[i] * rgc2_0;

        gr_xx_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 + 4.0 * ts_x_zzzz[i] * gfe_0 * gc_x[i] + 12.0 * ts_xx_zz[i] * gfe_0 + 8.0 * ts_xx_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_zzzz[i] * gfe_0 + ts_xx_zzzz[i] * rgc2_0;
    }

    // Set up 15-30 components of targeted buffer : DG

    auto gr_xy_xxxx = pbuffer.data(idx_g_dg + 15);

    auto gr_xy_xxxy = pbuffer.data(idx_g_dg + 16);

    auto gr_xy_xxxz = pbuffer.data(idx_g_dg + 17);

    auto gr_xy_xxyy = pbuffer.data(idx_g_dg + 18);

    auto gr_xy_xxyz = pbuffer.data(idx_g_dg + 19);

    auto gr_xy_xxzz = pbuffer.data(idx_g_dg + 20);

    auto gr_xy_xyyy = pbuffer.data(idx_g_dg + 21);

    auto gr_xy_xyyz = pbuffer.data(idx_g_dg + 22);

    auto gr_xy_xyzz = pbuffer.data(idx_g_dg + 23);

    auto gr_xy_xzzz = pbuffer.data(idx_g_dg + 24);

    auto gr_xy_yyyy = pbuffer.data(idx_g_dg + 25);

    auto gr_xy_yyyz = pbuffer.data(idx_g_dg + 26);

    auto gr_xy_yyzz = pbuffer.data(idx_g_dg + 27);

    auto gr_xy_yzzz = pbuffer.data(idx_g_dg + 28);

    auto gr_xy_zzzz = pbuffer.data(idx_g_dg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxxx, gr_xy_xxxy, gr_xy_xxxz, gr_xy_xxyy, gr_xy_xxyz, gr_xy_xxzz, gr_xy_xyyy, gr_xy_xyyz, gr_xy_xyzz, gr_xy_xzzz, gr_xy_yyyy, gr_xy_yyyz, gr_xy_yyzz, gr_xy_yzzz, gr_xy_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, ts_xy_xx, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xy, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yy, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zz, ts_xy_zzz, ts_xy_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xy_xxxx[i] = 8.0 * ts_y_xxx[i] * gfe_0 + 2.0 * ts_y_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_xy_xx[i] * gfe_0 + 8.0 * ts_xy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xy_xxxx[i] * gfe_0 + ts_xy_xxxx[i] * rgc2_0;

        gr_xy_xxxy[i] = 6.0 * ts_y_xxy[i] * gfe_0 + 2.0 * ts_y_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 + 2.0 * ts_x_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_xy[i] * gfe_0 + 6.0 * ts_xy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xxxy[i] * gfe_0 + ts_xy_xxxy[i] * rgc2_0;

        gr_xy_xxxz[i] = 6.0 * ts_y_xxz[i] * gfe_0 + 2.0 * ts_y_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_xz[i] * gfe_0 + 6.0 * ts_xy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xxxz[i] * gfe_0 + ts_xy_xxxz[i] * rgc2_0;

        gr_xy_xxyy[i] = 4.0 * ts_y_xyy[i] * gfe_0 + 2.0 * ts_y_xxyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xxy[i] * gfe_0 + 2.0 * ts_x_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 + 4.0 * ts_xy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xx[i] * gfe_0 + 4.0 * ts_xy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xxyy[i] * gfe_0 + ts_xy_xxyy[i] * rgc2_0;

        gr_xy_xxyz[i] = 4.0 * ts_y_xyz[i] * gfe_0 + 2.0 * ts_y_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxz[i] * gfe_0 + 2.0 * ts_x_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yz[i] * gfe_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xxyz[i] * gfe_0 + ts_xy_xxyz[i] * rgc2_0;

        gr_xy_xxzz[i] = 4.0 * ts_y_xzz[i] * gfe_0 + 2.0 * ts_y_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zz[i] * gfe_0 + 4.0 * ts_xy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xx[i] * gfe_0 + 4.0 * ts_xy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xxzz[i] * gfe_0 + ts_xy_xxzz[i] * rgc2_0;

        gr_xy_xyyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 + 2.0 * ts_y_xyyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xyy[i] * gfe_0 + 2.0 * ts_x_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xy[i] * gfe_0 + 6.0 * ts_xy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xyyy[i] * gfe_0 + ts_xy_xyyy[i] * rgc2_0;

        gr_xy_xyyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 + 2.0 * ts_y_xyyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xyz[i] * gfe_0 + 2.0 * ts_x_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xz[i] * gfe_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xyyz[i] * gfe_0 + ts_xy_xyyz[i] * rgc2_0;

        gr_xy_xyzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 + 2.0 * ts_y_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzz[i] * gfe_0 + 2.0 * ts_x_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xy[i] * gfe_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xyzz[i] * gfe_0 + ts_xy_xyzz[i] * rgc2_0;

        gr_xy_xzzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 + 2.0 * ts_y_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xz[i] * gfe_0 + 6.0 * ts_xy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xzzz[i] * gfe_0 + ts_xy_xzzz[i] * rgc2_0;

        gr_xy_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe_0 * gc_x[i] + 8.0 * ts_x_yyy[i] * gfe_0 + 2.0 * ts_x_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_xy_yy[i] * gfe_0 + 8.0 * ts_xy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_yyyy[i] * gfe_0 + ts_xy_yyyy[i] * rgc2_0;

        gr_xy_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_yyz[i] * gfe_0 + 2.0 * ts_x_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_yz[i] * gfe_0 + 6.0 * ts_xy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yyyz[i] * gfe_0 + ts_xy_yyyz[i] * rgc2_0;

        gr_xy_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_yzz[i] * gfe_0 + 2.0 * ts_x_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zz[i] * gfe_0 + 4.0 * ts_xy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 + 4.0 * ts_xy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yyzz[i] * gfe_0 + ts_xy_yyzz[i] * rgc2_0;

        gr_xy_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe_0 + 2.0 * ts_x_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xy_yz[i] * gfe_0 + 6.0 * ts_xy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yzzz[i] * gfe_0 + ts_xy_yzzz[i] * rgc2_0;

        gr_xy_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_xy_zz[i] * gfe_0 + 8.0 * ts_xy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_zzzz[i] * gfe_0 + ts_xy_zzzz[i] * rgc2_0;
    }

    // Set up 30-45 components of targeted buffer : DG

    auto gr_xz_xxxx = pbuffer.data(idx_g_dg + 30);

    auto gr_xz_xxxy = pbuffer.data(idx_g_dg + 31);

    auto gr_xz_xxxz = pbuffer.data(idx_g_dg + 32);

    auto gr_xz_xxyy = pbuffer.data(idx_g_dg + 33);

    auto gr_xz_xxyz = pbuffer.data(idx_g_dg + 34);

    auto gr_xz_xxzz = pbuffer.data(idx_g_dg + 35);

    auto gr_xz_xyyy = pbuffer.data(idx_g_dg + 36);

    auto gr_xz_xyyz = pbuffer.data(idx_g_dg + 37);

    auto gr_xz_xyzz = pbuffer.data(idx_g_dg + 38);

    auto gr_xz_xzzz = pbuffer.data(idx_g_dg + 39);

    auto gr_xz_yyyy = pbuffer.data(idx_g_dg + 40);

    auto gr_xz_yyyz = pbuffer.data(idx_g_dg + 41);

    auto gr_xz_yyzz = pbuffer.data(idx_g_dg + 42);

    auto gr_xz_yzzz = pbuffer.data(idx_g_dg + 43);

    auto gr_xz_zzzz = pbuffer.data(idx_g_dg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxxx, gr_xz_xxxy, gr_xz_xxxz, gr_xz_xxyy, gr_xz_xxyz, gr_xz_xxzz, gr_xz_xyyy, gr_xz_xyyz, gr_xz_xyzz, gr_xz_xzzz, gr_xz_yyyy, gr_xz_yyyz, gr_xz_yyzz, gr_xz_yzzz, gr_xz_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, ts_xz_xx, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xy, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yy, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zz, ts_xz_zzz, ts_xz_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xz_xxxx[i] = 8.0 * ts_z_xxx[i] * gfe_0 + 2.0 * ts_z_xxxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_xz_xx[i] * gfe_0 + 8.0 * ts_xz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xz_xxxx[i] * gfe_0 + ts_xz_xxxx[i] * rgc2_0;

        gr_xz_xxxy[i] = 6.0 * ts_z_xxy[i] * gfe_0 + 2.0 * ts_z_xxxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_xy[i] * gfe_0 + 6.0 * ts_xz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xxxy[i] * gfe_0 + ts_xz_xxxy[i] * rgc2_0;

        gr_xz_xxxz[i] = 6.0 * ts_z_xxz[i] * gfe_0 + 2.0 * ts_z_xxxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 + 2.0 * ts_x_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_xz[i] * gfe_0 + 6.0 * ts_xz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xxxz[i] * gfe_0 + ts_xz_xxxz[i] * rgc2_0;

        gr_xz_xxyy[i] = 4.0 * ts_z_xyy[i] * gfe_0 + 2.0 * ts_z_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yy[i] * gfe_0 + 4.0 * ts_xz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 + 4.0 * ts_xz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xxyy[i] * gfe_0 + ts_xz_xxyy[i] * rgc2_0;

        gr_xz_xxyz[i] = 4.0 * ts_z_xyz[i] * gfe_0 + 2.0 * ts_z_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxy[i] * gfe_0 + 2.0 * ts_x_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yz[i] * gfe_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xxyz[i] * gfe_0 + ts_xz_xxyz[i] * rgc2_0;

        gr_xz_xxzz[i] = 4.0 * ts_z_xzz[i] * gfe_0 + 2.0 * ts_z_xxzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xxz[i] * gfe_0 + 2.0 * ts_x_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zz[i] * gfe_0 + 4.0 * ts_xz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 + 4.0 * ts_xz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xxzz[i] * gfe_0 + ts_xz_xxzz[i] * rgc2_0;

        gr_xz_xyyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 + 2.0 * ts_z_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xy[i] * gfe_0 + 6.0 * ts_xz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xyyy[i] * gfe_0 + ts_xz_xyyy[i] * rgc2_0;

        gr_xz_xyyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 + 2.0 * ts_z_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xyy[i] * gfe_0 + 2.0 * ts_x_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xz[i] * gfe_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xyyz[i] * gfe_0 + ts_xz_xyyz[i] * rgc2_0;

        gr_xz_xyzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 + 2.0 * ts_z_xyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_xyz[i] * gfe_0 + 2.0 * ts_x_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_xy[i] * gfe_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xyzz[i] * gfe_0 + ts_xz_xyzz[i] * rgc2_0;

        gr_xz_xzzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 + 2.0 * ts_z_xzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xzz[i] * gfe_0 + 2.0 * ts_x_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xz[i] * gfe_0 + 6.0 * ts_xz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xzzz[i] * gfe_0 + ts_xz_xzzz[i] * rgc2_0;

        gr_xz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_xz_yy[i] * gfe_0 + 8.0 * ts_xz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_yyyy[i] * gfe_0 + ts_xz_yyyy[i] * rgc2_0;

        gr_xz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyy[i] * gfe_0 + 2.0 * ts_x_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xz_yz[i] * gfe_0 + 6.0 * ts_xz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yyyz[i] * gfe_0 + ts_xz_yyyz[i] * rgc2_0;

        gr_xz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_yyz[i] * gfe_0 + 2.0 * ts_x_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zz[i] * gfe_0 + 4.0 * ts_xz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_yy[i] * gfe_0 + 4.0 * ts_xz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yyzz[i] * gfe_0 + ts_xz_yyzz[i] * rgc2_0;

        gr_xz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_yzz[i] * gfe_0 + 2.0 * ts_x_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xz_yz[i] * gfe_0 + 6.0 * ts_xz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yzzz[i] * gfe_0 + ts_xz_yzzz[i] * rgc2_0;

        gr_xz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * gc_x[i] + 8.0 * ts_x_zzz[i] * gfe_0 + 2.0 * ts_x_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_xz_zz[i] * gfe_0 + 8.0 * ts_xz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_zzzz[i] * gfe_0 + ts_xz_zzzz[i] * rgc2_0;
    }

    // Set up 45-60 components of targeted buffer : DG

    auto gr_yy_xxxx = pbuffer.data(idx_g_dg + 45);

    auto gr_yy_xxxy = pbuffer.data(idx_g_dg + 46);

    auto gr_yy_xxxz = pbuffer.data(idx_g_dg + 47);

    auto gr_yy_xxyy = pbuffer.data(idx_g_dg + 48);

    auto gr_yy_xxyz = pbuffer.data(idx_g_dg + 49);

    auto gr_yy_xxzz = pbuffer.data(idx_g_dg + 50);

    auto gr_yy_xyyy = pbuffer.data(idx_g_dg + 51);

    auto gr_yy_xyyz = pbuffer.data(idx_g_dg + 52);

    auto gr_yy_xyzz = pbuffer.data(idx_g_dg + 53);

    auto gr_yy_xzzz = pbuffer.data(idx_g_dg + 54);

    auto gr_yy_yyyy = pbuffer.data(idx_g_dg + 55);

    auto gr_yy_yyyz = pbuffer.data(idx_g_dg + 56);

    auto gr_yy_yyzz = pbuffer.data(idx_g_dg + 57);

    auto gr_yy_yzzz = pbuffer.data(idx_g_dg + 58);

    auto gr_yy_zzzz = pbuffer.data(idx_g_dg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxxx, gr_yy_xxxy, gr_yy_xxxz, gr_yy_xxyy, gr_yy_xxyz, gr_yy_xxzz, gr_yy_xyyy, gr_yy_xyyz, gr_yy_xyzz, gr_yy_xzzz, gr_yy_yyyy, gr_yy_yyyz, gr_yy_yyzz, gr_yy_yzzz, gr_yy_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, ts_yy_xx, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xy, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yy, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yy_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 + 4.0 * ts_y_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_yy_xx[i] * gfe_0 + 8.0 * ts_yy_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yy_xxxx[i] * gfe_0 + ts_yy_xxxx[i] * rgc2_0;

        gr_yy_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 + 4.0 * ts_y_xxx[i] * gfe_0 + 4.0 * ts_y_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_xy[i] * gfe_0 + 6.0 * ts_yy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xxxy[i] * gfe_0 + ts_yy_xxxy[i] * rgc2_0;

        gr_yy_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 + 4.0 * ts_y_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_xz[i] * gfe_0 + 6.0 * ts_yy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xxxz[i] * gfe_0 + ts_yy_xxxz[i] * rgc2_0;

        gr_yy_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 + 8.0 * ts_y_xxy[i] * gfe_0 + 4.0 * ts_y_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 + 4.0 * ts_yy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xx[i] * gfe_0 + 4.0 * ts_yy_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xxyy[i] * gfe_0 + ts_yy_xxyy[i] * rgc2_0;

        gr_yy_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 + 4.0 * ts_y_xxz[i] * gfe_0 + 4.0 * ts_y_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yz[i] * gfe_0 + 4.0 * ts_yy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xxyz[i] * gfe_0 + ts_yy_xxyz[i] * rgc2_0;

        gr_yy_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 + 4.0 * ts_y_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zz[i] * gfe_0 + 4.0 * ts_yy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xx[i] * gfe_0 + 4.0 * ts_yy_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xxzz[i] * gfe_0 + ts_yy_xxzz[i] * rgc2_0;

        gr_yy_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 + 12.0 * ts_y_xyy[i] * gfe_0 + 4.0 * ts_y_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yy_xy[i] * gfe_0 + 6.0 * ts_yy_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xyyy[i] * gfe_0 + ts_yy_xyyy[i] * rgc2_0;

        gr_yy_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 + 8.0 * ts_y_xyz[i] * gfe_0 + 4.0 * ts_y_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xz[i] * gfe_0 + 4.0 * ts_yy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xyyz[i] * gfe_0 + ts_yy_xyyz[i] * rgc2_0;

        gr_yy_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 + 4.0 * ts_y_xzz[i] * gfe_0 + 4.0 * ts_y_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xy[i] * gfe_0 + 4.0 * ts_yy_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xyzz[i] * gfe_0 + ts_yy_xyzz[i] * rgc2_0;

        gr_yy_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 + 4.0 * ts_y_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yy_xz[i] * gfe_0 + 6.0 * ts_yy_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xzzz[i] * gfe_0 + ts_yy_xzzz[i] * rgc2_0;

        gr_yy_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 + 16.0 * ts_y_yyy[i] * gfe_0 + 4.0 * ts_y_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_yy_yy[i] * gfe_0 + 8.0 * ts_yy_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_yyyy[i] * gfe_0 + ts_yy_yyyy[i] * rgc2_0;

        gr_yy_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 + 12.0 * ts_y_yyz[i] * gfe_0 + 4.0 * ts_y_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_yz[i] * gfe_0 + 6.0 * ts_yy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yyyz[i] * gfe_0 + ts_yy_yyyz[i] * rgc2_0;

        gr_yy_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 + 8.0 * ts_y_yzz[i] * gfe_0 + 4.0 * ts_y_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zz[i] * gfe_0 + 4.0 * ts_yy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 + 4.0 * ts_yy_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yyzz[i] * gfe_0 + ts_yy_yyzz[i] * rgc2_0;

        gr_yy_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 + 4.0 * ts_y_zzz[i] * gfe_0 + 4.0 * ts_y_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yy_yz[i] * gfe_0 + 6.0 * ts_yy_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yzzz[i] * gfe_0 + ts_yy_yzzz[i] * rgc2_0;

        gr_yy_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 + 4.0 * ts_y_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_yy_zz[i] * gfe_0 + 8.0 * ts_yy_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_zzzz[i] * gfe_0 + ts_yy_zzzz[i] * rgc2_0;
    }

    // Set up 60-75 components of targeted buffer : DG

    auto gr_yz_xxxx = pbuffer.data(idx_g_dg + 60);

    auto gr_yz_xxxy = pbuffer.data(idx_g_dg + 61);

    auto gr_yz_xxxz = pbuffer.data(idx_g_dg + 62);

    auto gr_yz_xxyy = pbuffer.data(idx_g_dg + 63);

    auto gr_yz_xxyz = pbuffer.data(idx_g_dg + 64);

    auto gr_yz_xxzz = pbuffer.data(idx_g_dg + 65);

    auto gr_yz_xyyy = pbuffer.data(idx_g_dg + 66);

    auto gr_yz_xyyz = pbuffer.data(idx_g_dg + 67);

    auto gr_yz_xyzz = pbuffer.data(idx_g_dg + 68);

    auto gr_yz_xzzz = pbuffer.data(idx_g_dg + 69);

    auto gr_yz_yyyy = pbuffer.data(idx_g_dg + 70);

    auto gr_yz_yyyz = pbuffer.data(idx_g_dg + 71);

    auto gr_yz_yyzz = pbuffer.data(idx_g_dg + 72);

    auto gr_yz_yzzz = pbuffer.data(idx_g_dg + 73);

    auto gr_yz_zzzz = pbuffer.data(idx_g_dg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxxx, gr_yz_xxxy, gr_yz_xxxz, gr_yz_xxyy, gr_yz_xxyz, gr_yz_xxzz, gr_yz_xyyy, gr_yz_xyyz, gr_yz_xyzz, gr_yz_xzzz, gr_yz_yyyy, gr_yz_yyyz, gr_yz_yyzz, gr_yz_yzzz, gr_yz_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, ts_yz_xx, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xy, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yy, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zz, ts_yz_zzz, ts_yz_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_yz_xx[i] * gfe_0 + 8.0 * ts_yz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yz_xxxx[i] * gfe_0 + ts_yz_xxxx[i] * rgc2_0;

        gr_yz_xxxy[i] = 2.0 * ts_z_xxx[i] * gfe_0 + 2.0 * ts_z_xxxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_xy[i] * gfe_0 + 6.0 * ts_yz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xxxy[i] * gfe_0 + ts_yz_xxxy[i] * rgc2_0;

        gr_yz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxx[i] * gfe_0 + 2.0 * ts_y_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_xz[i] * gfe_0 + 6.0 * ts_yz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xxxz[i] * gfe_0 + ts_yz_xxxz[i] * rgc2_0;

        gr_yz_xxyy[i] = 4.0 * ts_z_xxy[i] * gfe_0 + 2.0 * ts_z_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yy[i] * gfe_0 + 4.0 * ts_yz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xx[i] * gfe_0 + 4.0 * ts_yz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xxyy[i] * gfe_0 + ts_yz_xxyy[i] * rgc2_0;

        gr_yz_xxyz[i] = 2.0 * ts_z_xxz[i] * gfe_0 + 2.0 * ts_z_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxy[i] * gfe_0 + 2.0 * ts_y_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yz[i] * gfe_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xxyz[i] * gfe_0 + ts_yz_xxyz[i] * rgc2_0;

        gr_yz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_xxz[i] * gfe_0 + 2.0 * ts_y_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zz[i] * gfe_0 + 4.0 * ts_yz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xx[i] * gfe_0 + 4.0 * ts_yz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xxzz[i] * gfe_0 + ts_yz_xxzz[i] * rgc2_0;

        gr_yz_xyyy[i] = 6.0 * ts_z_xyy[i] * gfe_0 + 2.0 * ts_z_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_yz_xy[i] * gfe_0 + 6.0 * ts_yz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xyyy[i] * gfe_0 + ts_yz_xyyy[i] * rgc2_0;

        gr_yz_xyyz[i] = 4.0 * ts_z_xyz[i] * gfe_0 + 2.0 * ts_z_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyy[i] * gfe_0 + 2.0 * ts_y_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xz[i] * gfe_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xyyz[i] * gfe_0 + ts_yz_xyyz[i] * rgc2_0;

        gr_yz_xyzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 + 2.0 * ts_z_xyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_xyz[i] * gfe_0 + 2.0 * ts_y_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_xy[i] * gfe_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xyzz[i] * gfe_0 + ts_yz_xyzz[i] * rgc2_0;

        gr_yz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_xzz[i] * gfe_0 + 2.0 * ts_y_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_yz_xz[i] * gfe_0 + 6.0 * ts_yz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xzzz[i] * gfe_0 + ts_yz_xzzz[i] * rgc2_0;

        gr_yz_yyyy[i] = 8.0 * ts_z_yyy[i] * gfe_0 + 2.0 * ts_z_yyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_yz_yy[i] * gfe_0 + 8.0 * ts_yz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_yyyy[i] * gfe_0 + ts_yz_yyyy[i] * rgc2_0;

        gr_yz_yyyz[i] = 6.0 * ts_z_yyz[i] * gfe_0 + 2.0 * ts_z_yyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyy[i] * gfe_0 + 2.0 * ts_y_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yz_yz[i] * gfe_0 + 6.0 * ts_yz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yyyz[i] * gfe_0 + ts_yz_yyyz[i] * rgc2_0;

        gr_yz_yyzz[i] = 4.0 * ts_z_yzz[i] * gfe_0 + 2.0 * ts_z_yyzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_yyz[i] * gfe_0 + 2.0 * ts_y_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zz[i] * gfe_0 + 4.0 * ts_yz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_yy[i] * gfe_0 + 4.0 * ts_yz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yyzz[i] * gfe_0 + ts_yz_yyzz[i] * rgc2_0;

        gr_yz_yzzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 + 2.0 * ts_z_yzzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_yzz[i] * gfe_0 + 2.0 * ts_y_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yz[i] * gfe_0 + 6.0 * ts_yz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yzzz[i] * gfe_0 + ts_yz_yzzz[i] * rgc2_0;

        gr_yz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * gc_y[i] + 8.0 * ts_y_zzz[i] * gfe_0 + 2.0 * ts_y_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_yz_zz[i] * gfe_0 + 8.0 * ts_yz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_zzzz[i] * gfe_0 + ts_yz_zzzz[i] * rgc2_0;
    }

    // Set up 75-90 components of targeted buffer : DG

    auto gr_zz_xxxx = pbuffer.data(idx_g_dg + 75);

    auto gr_zz_xxxy = pbuffer.data(idx_g_dg + 76);

    auto gr_zz_xxxz = pbuffer.data(idx_g_dg + 77);

    auto gr_zz_xxyy = pbuffer.data(idx_g_dg + 78);

    auto gr_zz_xxyz = pbuffer.data(idx_g_dg + 79);

    auto gr_zz_xxzz = pbuffer.data(idx_g_dg + 80);

    auto gr_zz_xyyy = pbuffer.data(idx_g_dg + 81);

    auto gr_zz_xyyz = pbuffer.data(idx_g_dg + 82);

    auto gr_zz_xyzz = pbuffer.data(idx_g_dg + 83);

    auto gr_zz_xzzz = pbuffer.data(idx_g_dg + 84);

    auto gr_zz_yyyy = pbuffer.data(idx_g_dg + 85);

    auto gr_zz_yyyz = pbuffer.data(idx_g_dg + 86);

    auto gr_zz_yyzz = pbuffer.data(idx_g_dg + 87);

    auto gr_zz_yzzz = pbuffer.data(idx_g_dg + 88);

    auto gr_zz_zzzz = pbuffer.data(idx_g_dg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxxx, gr_zz_xxxy, gr_zz_xxxz, gr_zz_xxyy, gr_zz_xxyz, gr_zz_xxzz, gr_zz_xyyy, gr_zz_xyyz, gr_zz_xyzz, gr_zz_xzzz, gr_zz_yyyy, gr_zz_yyyz, gr_zz_yyzz, gr_zz_yzzz, gr_zz_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, ts_zz_xx, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xy, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yy, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zz_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 + 4.0 * ts_z_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_zz_xx[i] * gfe_0 + 8.0 * ts_zz_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zz_xxxx[i] * gfe_0 + ts_zz_xxxx[i] * rgc2_0;

        gr_zz_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 + 4.0 * ts_z_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_xy[i] * gfe_0 + 6.0 * ts_zz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xxxy[i] * gfe_0 + ts_zz_xxxy[i] * rgc2_0;

        gr_zz_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 + 4.0 * ts_z_xxx[i] * gfe_0 + 4.0 * ts_z_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_xz[i] * gfe_0 + 6.0 * ts_zz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xxxz[i] * gfe_0 + ts_zz_xxxz[i] * rgc2_0;

        gr_zz_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 + 4.0 * ts_z_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yy[i] * gfe_0 + 4.0 * ts_zz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xx[i] * gfe_0 + 4.0 * ts_zz_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xxyy[i] * gfe_0 + ts_zz_xxyy[i] * rgc2_0;

        gr_zz_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 + 4.0 * ts_z_xxy[i] * gfe_0 + 4.0 * ts_z_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yz[i] * gfe_0 + 4.0 * ts_zz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xxyz[i] * gfe_0 + ts_zz_xxyz[i] * rgc2_0;

        gr_zz_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 + 8.0 * ts_z_xxz[i] * gfe_0 + 4.0 * ts_z_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zz[i] * gfe_0 + 4.0 * ts_zz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xx[i] * gfe_0 + 4.0 * ts_zz_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xxzz[i] * gfe_0 + ts_zz_xxzz[i] * rgc2_0;

        gr_zz_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 + 4.0 * ts_z_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_zz_xy[i] * gfe_0 + 6.0 * ts_zz_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xyyy[i] * gfe_0 + ts_zz_xyyy[i] * rgc2_0;

        gr_zz_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 + 4.0 * ts_z_xyy[i] * gfe_0 + 4.0 * ts_z_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xz[i] * gfe_0 + 4.0 * ts_zz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xyyz[i] * gfe_0 + ts_zz_xyyz[i] * rgc2_0;

        gr_zz_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 + 8.0 * ts_z_xyz[i] * gfe_0 + 4.0 * ts_z_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_xy[i] * gfe_0 + 4.0 * ts_zz_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xyzz[i] * gfe_0 + ts_zz_xyzz[i] * rgc2_0;

        gr_zz_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 + 12.0 * ts_z_xzz[i] * gfe_0 + 4.0 * ts_z_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_zz_xz[i] * gfe_0 + 6.0 * ts_zz_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xzzz[i] * gfe_0 + ts_zz_xzzz[i] * rgc2_0;

        gr_zz_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 + 4.0 * ts_z_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_zz_yy[i] * gfe_0 + 8.0 * ts_zz_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_yyyy[i] * gfe_0 + ts_zz_yyyy[i] * rgc2_0;

        gr_zz_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 + 4.0 * ts_z_yyy[i] * gfe_0 + 4.0 * ts_z_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zz_yz[i] * gfe_0 + 6.0 * ts_zz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yyyz[i] * gfe_0 + ts_zz_yyyz[i] * rgc2_0;

        gr_zz_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 + 8.0 * ts_z_yyz[i] * gfe_0 + 4.0 * ts_z_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zz[i] * gfe_0 + 4.0 * ts_zz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_yy[i] * gfe_0 + 4.0 * ts_zz_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yyzz[i] * gfe_0 + ts_zz_yyzz[i] * rgc2_0;

        gr_zz_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 + 12.0 * ts_z_yzz[i] * gfe_0 + 4.0 * ts_z_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_zz_yz[i] * gfe_0 + 6.0 * ts_zz_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yzzz[i] * gfe_0 + ts_zz_yzzz[i] * rgc2_0;

        gr_zz_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 + 16.0 * ts_z_zzz[i] * gfe_0 + 4.0 * ts_z_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_zz_zz[i] * gfe_0 + 8.0 * ts_zz_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_zzzz[i] * gfe_0 + ts_zz_zzzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

