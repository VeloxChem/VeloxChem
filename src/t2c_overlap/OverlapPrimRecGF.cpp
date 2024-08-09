#include "OverlapPrimRecGF.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_gf(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_gf,
                     const size_t              idx_ovl_df,
                     const size_t              idx_ovl_fd,
                     const size_t              idx_ovl_ff,
                     const CSimdArray<double>& factors,
                     const size_t              idx_rpa,
                     const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_ovl_df);

    auto ts_xx_xxy = pbuffer.data(idx_ovl_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_ovl_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_ovl_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_ovl_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_ovl_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_ovl_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_ovl_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_ovl_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_ovl_df + 9);

    auto ts_xy_yyy = pbuffer.data(idx_ovl_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_ovl_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_ovl_df + 18);

    auto ts_xz_yyz = pbuffer.data(idx_ovl_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_ovl_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_ovl_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_ovl_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_ovl_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_ovl_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_ovl_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_ovl_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_ovl_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_ovl_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_ovl_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_ovl_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_ovl_df + 39);

    auto ts_yz_xxz = pbuffer.data(idx_ovl_df + 42);

    auto ts_yz_xzz = pbuffer.data(idx_ovl_df + 45);

    auto ts_yz_yyz = pbuffer.data(idx_ovl_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_ovl_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_ovl_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_ovl_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_ovl_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_ovl_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_ovl_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_ovl_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_ovl_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_ovl_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_ovl_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_ovl_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_ovl_df + 59);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xzz_xz = pbuffer.data(idx_ovl_fd + 32);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

    auto ts_yyz_xz = pbuffer.data(idx_ovl_fd + 44);

    auto ts_yyz_yz = pbuffer.data(idx_ovl_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_ovl_fd + 47);

    auto ts_yzz_xy = pbuffer.data(idx_ovl_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_ovl_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

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

    auto ts_xxy_xxy = pbuffer.data(idx_ovl_ff + 11);

    auto ts_xxy_xxz = pbuffer.data(idx_ovl_ff + 12);

    auto ts_xxy_xyy = pbuffer.data(idx_ovl_ff + 13);

    auto ts_xxy_xzz = pbuffer.data(idx_ovl_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ovl_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ovl_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ovl_ff + 18);

    auto ts_xxz_xxx = pbuffer.data(idx_ovl_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ovl_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ovl_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ovl_ff + 23);

    auto ts_xxz_xyz = pbuffer.data(idx_ovl_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ovl_ff + 25);

    auto ts_xxz_yyz = pbuffer.data(idx_ovl_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ovl_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ovl_ff + 29);

    auto ts_xyy_xxx = pbuffer.data(idx_ovl_ff + 30);

    auto ts_xyy_xxy = pbuffer.data(idx_ovl_ff + 31);

    auto ts_xyy_xyy = pbuffer.data(idx_ovl_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ovl_ff + 34);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ovl_ff + 39);

    auto ts_xyz_yyz = pbuffer.data(idx_ovl_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ovl_ff + 48);

    auto ts_xzz_xxx = pbuffer.data(idx_ovl_ff + 50);

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

    auto ts_yyz_xxz = pbuffer.data(idx_ovl_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ovl_ff + 73);

    auto ts_yyz_xyz = pbuffer.data(idx_ovl_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ovl_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ovl_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ovl_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ovl_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ovl_ff + 79);

    auto ts_yzz_xxx = pbuffer.data(idx_ovl_ff + 80);

    auto ts_yzz_xxy = pbuffer.data(idx_ovl_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ovl_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ovl_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ovl_ff + 86);

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

    // Set up 0-10 components of targeted buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

#pragma omp simd aligned(pa_x,            \
                             ts_xx_xxx,   \
                             ts_xx_xxy,   \
                             ts_xx_xxz,   \
                             ts_xx_xyy,   \
                             ts_xx_xyz,   \
                             ts_xx_xzz,   \
                             ts_xx_yyy,   \
                             ts_xx_yyz,   \
                             ts_xx_yzz,   \
                             ts_xx_zzz,   \
                             ts_xxx_xx,   \
                             ts_xxx_xxx,  \
                             ts_xxx_xxy,  \
                             ts_xxx_xxz,  \
                             ts_xxx_xy,   \
                             ts_xxx_xyy,  \
                             ts_xxx_xyz,  \
                             ts_xxx_xz,   \
                             ts_xxx_xzz,  \
                             ts_xxx_yy,   \
                             ts_xxx_yyy,  \
                             ts_xxx_yyz,  \
                             ts_xxx_yz,   \
                             ts_xxx_yzz,  \
                             ts_xxx_zz,   \
                             ts_xxx_zzz,  \
                             ts_xxxx_xxx, \
                             ts_xxxx_xxy, \
                             ts_xxxx_xxz, \
                             ts_xxxx_xyy, \
                             ts_xxxx_xyz, \
                             ts_xxxx_xzz, \
                             ts_xxxx_yyy, \
                             ts_xxxx_yyz, \
                             ts_xxxx_yzz, \
                             ts_xxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxx_xxx[i] = 3.0 * ts_xx_xxx[i] * fe_0 + 3.0 * ts_xxx_xx[i] * fe_0 + ts_xxx_xxx[i] * pa_x[i];

        ts_xxxx_xxy[i] = 3.0 * ts_xx_xxy[i] * fe_0 + 2.0 * ts_xxx_xy[i] * fe_0 + ts_xxx_xxy[i] * pa_x[i];

        ts_xxxx_xxz[i] = 3.0 * ts_xx_xxz[i] * fe_0 + 2.0 * ts_xxx_xz[i] * fe_0 + ts_xxx_xxz[i] * pa_x[i];

        ts_xxxx_xyy[i] = 3.0 * ts_xx_xyy[i] * fe_0 + ts_xxx_yy[i] * fe_0 + ts_xxx_xyy[i] * pa_x[i];

        ts_xxxx_xyz[i] = 3.0 * ts_xx_xyz[i] * fe_0 + ts_xxx_yz[i] * fe_0 + ts_xxx_xyz[i] * pa_x[i];

        ts_xxxx_xzz[i] = 3.0 * ts_xx_xzz[i] * fe_0 + ts_xxx_zz[i] * fe_0 + ts_xxx_xzz[i] * pa_x[i];

        ts_xxxx_yyy[i] = 3.0 * ts_xx_yyy[i] * fe_0 + ts_xxx_yyy[i] * pa_x[i];

        ts_xxxx_yyz[i] = 3.0 * ts_xx_yyz[i] * fe_0 + ts_xxx_yyz[i] * pa_x[i];

        ts_xxxx_yzz[i] = 3.0 * ts_xx_yzz[i] * fe_0 + ts_xxx_yzz[i] * pa_x[i];

        ts_xxxx_zzz[i] = 3.0 * ts_xx_zzz[i] * fe_0 + ts_xxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto ts_xxxy_xxx = pbuffer.data(idx_ovl_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_ovl_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_ovl_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_ovl_gf + 13);

    auto ts_xxxy_xyz = pbuffer.data(idx_ovl_gf + 14);

    auto ts_xxxy_xzz = pbuffer.data(idx_ovl_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_ovl_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_ovl_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_ovl_gf + 18);

    auto ts_xxxy_zzz = pbuffer.data(idx_ovl_gf + 19);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             ts_xxx_xx,   \
                             ts_xxx_xxx,  \
                             ts_xxx_xxy,  \
                             ts_xxx_xxz,  \
                             ts_xxx_xy,   \
                             ts_xxx_xyy,  \
                             ts_xxx_xyz,  \
                             ts_xxx_xz,   \
                             ts_xxx_xzz,  \
                             ts_xxx_zzz,  \
                             ts_xxxy_xxx, \
                             ts_xxxy_xxy, \
                             ts_xxxy_xxz, \
                             ts_xxxy_xyy, \
                             ts_xxxy_xyz, \
                             ts_xxxy_xzz, \
                             ts_xxxy_yyy, \
                             ts_xxxy_yyz, \
                             ts_xxxy_yzz, \
                             ts_xxxy_zzz, \
                             ts_xxy_yyy,  \
                             ts_xxy_yyz,  \
                             ts_xxy_yzz,  \
                             ts_xy_yyy,   \
                             ts_xy_yyz,   \
                             ts_xy_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxy_xxx[i] = ts_xxx_xxx[i] * pa_y[i];

        ts_xxxy_xxy[i] = ts_xxx_xx[i] * fe_0 + ts_xxx_xxy[i] * pa_y[i];

        ts_xxxy_xxz[i] = ts_xxx_xxz[i] * pa_y[i];

        ts_xxxy_xyy[i] = 2.0 * ts_xxx_xy[i] * fe_0 + ts_xxx_xyy[i] * pa_y[i];

        ts_xxxy_xyz[i] = ts_xxx_xz[i] * fe_0 + ts_xxx_xyz[i] * pa_y[i];

        ts_xxxy_xzz[i] = ts_xxx_xzz[i] * pa_y[i];

        ts_xxxy_yyy[i] = 2.0 * ts_xy_yyy[i] * fe_0 + ts_xxy_yyy[i] * pa_x[i];

        ts_xxxy_yyz[i] = 2.0 * ts_xy_yyz[i] * fe_0 + ts_xxy_yyz[i] * pa_x[i];

        ts_xxxy_yzz[i] = 2.0 * ts_xy_yzz[i] * fe_0 + ts_xxy_yzz[i] * pa_x[i];

        ts_xxxy_zzz[i] = ts_xxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto ts_xxxz_xxx = pbuffer.data(idx_ovl_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_ovl_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_ovl_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_ovl_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_ovl_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_ovl_gf + 25);

    auto ts_xxxz_yyy = pbuffer.data(idx_ovl_gf + 26);

    auto ts_xxxz_yyz = pbuffer.data(idx_ovl_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_ovl_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_ovl_gf + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             ts_xxx_xx,   \
                             ts_xxx_xxx,  \
                             ts_xxx_xxy,  \
                             ts_xxx_xxz,  \
                             ts_xxx_xy,   \
                             ts_xxx_xyy,  \
                             ts_xxx_xyz,  \
                             ts_xxx_xz,   \
                             ts_xxx_xzz,  \
                             ts_xxx_yyy,  \
                             ts_xxxz_xxx, \
                             ts_xxxz_xxy, \
                             ts_xxxz_xxz, \
                             ts_xxxz_xyy, \
                             ts_xxxz_xyz, \
                             ts_xxxz_xzz, \
                             ts_xxxz_yyy, \
                             ts_xxxz_yyz, \
                             ts_xxxz_yzz, \
                             ts_xxxz_zzz, \
                             ts_xxz_yyz,  \
                             ts_xxz_yzz,  \
                             ts_xxz_zzz,  \
                             ts_xz_yyz,   \
                             ts_xz_yzz,   \
                             ts_xz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxz_xxx[i] = ts_xxx_xxx[i] * pa_z[i];

        ts_xxxz_xxy[i] = ts_xxx_xxy[i] * pa_z[i];

        ts_xxxz_xxz[i] = ts_xxx_xx[i] * fe_0 + ts_xxx_xxz[i] * pa_z[i];

        ts_xxxz_xyy[i] = ts_xxx_xyy[i] * pa_z[i];

        ts_xxxz_xyz[i] = ts_xxx_xy[i] * fe_0 + ts_xxx_xyz[i] * pa_z[i];

        ts_xxxz_xzz[i] = 2.0 * ts_xxx_xz[i] * fe_0 + ts_xxx_xzz[i] * pa_z[i];

        ts_xxxz_yyy[i] = ts_xxx_yyy[i] * pa_z[i];

        ts_xxxz_yyz[i] = 2.0 * ts_xz_yyz[i] * fe_0 + ts_xxz_yyz[i] * pa_x[i];

        ts_xxxz_yzz[i] = 2.0 * ts_xz_yzz[i] * fe_0 + ts_xxz_yzz[i] * pa_x[i];

        ts_xxxz_zzz[i] = 2.0 * ts_xz_zzz[i] * fe_0 + ts_xxz_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto ts_xxyy_xxx = pbuffer.data(idx_ovl_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_ovl_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_ovl_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_ovl_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_ovl_gf + 39);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             ts_xx_xxx,   \
                             ts_xx_xxz,   \
                             ts_xx_xzz,   \
                             ts_xxy_xxx,  \
                             ts_xxy_xxz,  \
                             ts_xxy_xzz,  \
                             ts_xxyy_xxx, \
                             ts_xxyy_xxy, \
                             ts_xxyy_xxz, \
                             ts_xxyy_xyy, \
                             ts_xxyy_xyz, \
                             ts_xxyy_xzz, \
                             ts_xxyy_yyy, \
                             ts_xxyy_yyz, \
                             ts_xxyy_yzz, \
                             ts_xxyy_zzz, \
                             ts_xyy_xxy,  \
                             ts_xyy_xy,   \
                             ts_xyy_xyy,  \
                             ts_xyy_xyz,  \
                             ts_xyy_yy,   \
                             ts_xyy_yyy,  \
                             ts_xyy_yyz,  \
                             ts_xyy_yz,   \
                             ts_xyy_yzz,  \
                             ts_xyy_zzz,  \
                             ts_yy_xxy,   \
                             ts_yy_xyy,   \
                             ts_yy_xyz,   \
                             ts_yy_yyy,   \
                             ts_yy_yyz,   \
                             ts_yy_yzz,   \
                             ts_yy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyy_xxx[i] = ts_xx_xxx[i] * fe_0 + ts_xxy_xxx[i] * pa_y[i];

        ts_xxyy_xxy[i] = ts_yy_xxy[i] * fe_0 + 2.0 * ts_xyy_xy[i] * fe_0 + ts_xyy_xxy[i] * pa_x[i];

        ts_xxyy_xxz[i] = ts_xx_xxz[i] * fe_0 + ts_xxy_xxz[i] * pa_y[i];

        ts_xxyy_xyy[i] = ts_yy_xyy[i] * fe_0 + ts_xyy_yy[i] * fe_0 + ts_xyy_xyy[i] * pa_x[i];

        ts_xxyy_xyz[i] = ts_yy_xyz[i] * fe_0 + ts_xyy_yz[i] * fe_0 + ts_xyy_xyz[i] * pa_x[i];

        ts_xxyy_xzz[i] = ts_xx_xzz[i] * fe_0 + ts_xxy_xzz[i] * pa_y[i];

        ts_xxyy_yyy[i] = ts_yy_yyy[i] * fe_0 + ts_xyy_yyy[i] * pa_x[i];

        ts_xxyy_yyz[i] = ts_yy_yyz[i] * fe_0 + ts_xyy_yyz[i] * pa_x[i];

        ts_xxyy_yzz[i] = ts_yy_yzz[i] * fe_0 + ts_xyy_yzz[i] * pa_x[i];

        ts_xxyy_zzz[i] = ts_yy_zzz[i] * fe_0 + ts_xyy_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto ts_xxyz_xxx = pbuffer.data(idx_ovl_gf + 40);

    auto ts_xxyz_xxy = pbuffer.data(idx_ovl_gf + 41);

    auto ts_xxyz_xxz = pbuffer.data(idx_ovl_gf + 42);

    auto ts_xxyz_xyy = pbuffer.data(idx_ovl_gf + 43);

    auto ts_xxyz_xyz = pbuffer.data(idx_ovl_gf + 44);

    auto ts_xxyz_xzz = pbuffer.data(idx_ovl_gf + 45);

    auto ts_xxyz_yyy = pbuffer.data(idx_ovl_gf + 46);

    auto ts_xxyz_yyz = pbuffer.data(idx_ovl_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_ovl_gf + 48);

    auto ts_xxyz_zzz = pbuffer.data(idx_ovl_gf + 49);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             ts_xxy_xxy,  \
                             ts_xxy_xyy,  \
                             ts_xxy_yyy,  \
                             ts_xxyz_xxx, \
                             ts_xxyz_xxy, \
                             ts_xxyz_xxz, \
                             ts_xxyz_xyy, \
                             ts_xxyz_xyz, \
                             ts_xxyz_xzz, \
                             ts_xxyz_yyy, \
                             ts_xxyz_yyz, \
                             ts_xxyz_yzz, \
                             ts_xxyz_zzz, \
                             ts_xxz_xxx,  \
                             ts_xxz_xxz,  \
                             ts_xxz_xyz,  \
                             ts_xxz_xz,   \
                             ts_xxz_xzz,  \
                             ts_xxz_zzz,  \
                             ts_xyz_yyz,  \
                             ts_xyz_yzz,  \
                             ts_yz_yyz,   \
                             ts_yz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyz_xxx[i] = ts_xxz_xxx[i] * pa_y[i];

        ts_xxyz_xxy[i] = ts_xxy_xxy[i] * pa_z[i];

        ts_xxyz_xxz[i] = ts_xxz_xxz[i] * pa_y[i];

        ts_xxyz_xyy[i] = ts_xxy_xyy[i] * pa_z[i];

        ts_xxyz_xyz[i] = ts_xxz_xz[i] * fe_0 + ts_xxz_xyz[i] * pa_y[i];

        ts_xxyz_xzz[i] = ts_xxz_xzz[i] * pa_y[i];

        ts_xxyz_yyy[i] = ts_xxy_yyy[i] * pa_z[i];

        ts_xxyz_yyz[i] = ts_yz_yyz[i] * fe_0 + ts_xyz_yyz[i] * pa_x[i];

        ts_xxyz_yzz[i] = ts_yz_yzz[i] * fe_0 + ts_xyz_yzz[i] * pa_x[i];

        ts_xxyz_zzz[i] = ts_xxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_ovl_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_ovl_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_ovl_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_ovl_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             ts_xx_xxx,   \
                             ts_xx_xxy,   \
                             ts_xx_xyy,   \
                             ts_xxz_xxx,  \
                             ts_xxz_xxy,  \
                             ts_xxz_xyy,  \
                             ts_xxzz_xxx, \
                             ts_xxzz_xxy, \
                             ts_xxzz_xxz, \
                             ts_xxzz_xyy, \
                             ts_xxzz_xyz, \
                             ts_xxzz_xzz, \
                             ts_xxzz_yyy, \
                             ts_xxzz_yyz, \
                             ts_xxzz_yzz, \
                             ts_xxzz_zzz, \
                             ts_xzz_xxz,  \
                             ts_xzz_xyz,  \
                             ts_xzz_xz,   \
                             ts_xzz_xzz,  \
                             ts_xzz_yyy,  \
                             ts_xzz_yyz,  \
                             ts_xzz_yz,   \
                             ts_xzz_yzz,  \
                             ts_xzz_zz,   \
                             ts_xzz_zzz,  \
                             ts_zz_xxz,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyy,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzz_xxx[i] = ts_xx_xxx[i] * fe_0 + ts_xxz_xxx[i] * pa_z[i];

        ts_xxzz_xxy[i] = ts_xx_xxy[i] * fe_0 + ts_xxz_xxy[i] * pa_z[i];

        ts_xxzz_xxz[i] = ts_zz_xxz[i] * fe_0 + 2.0 * ts_xzz_xz[i] * fe_0 + ts_xzz_xxz[i] * pa_x[i];

        ts_xxzz_xyy[i] = ts_xx_xyy[i] * fe_0 + ts_xxz_xyy[i] * pa_z[i];

        ts_xxzz_xyz[i] = ts_zz_xyz[i] * fe_0 + ts_xzz_yz[i] * fe_0 + ts_xzz_xyz[i] * pa_x[i];

        ts_xxzz_xzz[i] = ts_zz_xzz[i] * fe_0 + ts_xzz_zz[i] * fe_0 + ts_xzz_xzz[i] * pa_x[i];

        ts_xxzz_yyy[i] = ts_zz_yyy[i] * fe_0 + ts_xzz_yyy[i] * pa_x[i];

        ts_xxzz_yyz[i] = ts_zz_yyz[i] * fe_0 + ts_xzz_yyz[i] * pa_x[i];

        ts_xxzz_yzz[i] = ts_zz_yzz[i] * fe_0 + ts_xzz_yzz[i] * pa_x[i];

        ts_xxzz_zzz[i] = ts_zz_zzz[i] * fe_0 + ts_xzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto ts_xyyy_xxx = pbuffer.data(idx_ovl_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_ovl_gf + 61);

    auto ts_xyyy_xxz = pbuffer.data(idx_ovl_gf + 62);

    auto ts_xyyy_xyy = pbuffer.data(idx_ovl_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_ovl_gf + 64);

    auto ts_xyyy_xzz = pbuffer.data(idx_ovl_gf + 65);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_ovl_gf + 69);

#pragma omp simd aligned(pa_x,            \
                             ts_xyyy_xxx, \
                             ts_xyyy_xxy, \
                             ts_xyyy_xxz, \
                             ts_xyyy_xyy, \
                             ts_xyyy_xyz, \
                             ts_xyyy_xzz, \
                             ts_xyyy_yyy, \
                             ts_xyyy_yyz, \
                             ts_xyyy_yzz, \
                             ts_xyyy_zzz, \
                             ts_yyy_xx,   \
                             ts_yyy_xxx,  \
                             ts_yyy_xxy,  \
                             ts_yyy_xxz,  \
                             ts_yyy_xy,   \
                             ts_yyy_xyy,  \
                             ts_yyy_xyz,  \
                             ts_yyy_xz,   \
                             ts_yyy_xzz,  \
                             ts_yyy_yy,   \
                             ts_yyy_yyy,  \
                             ts_yyy_yyz,  \
                             ts_yyy_yz,   \
                             ts_yyy_yzz,  \
                             ts_yyy_zz,   \
                             ts_yyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyy_xxx[i] = 3.0 * ts_yyy_xx[i] * fe_0 + ts_yyy_xxx[i] * pa_x[i];

        ts_xyyy_xxy[i] = 2.0 * ts_yyy_xy[i] * fe_0 + ts_yyy_xxy[i] * pa_x[i];

        ts_xyyy_xxz[i] = 2.0 * ts_yyy_xz[i] * fe_0 + ts_yyy_xxz[i] * pa_x[i];

        ts_xyyy_xyy[i] = ts_yyy_yy[i] * fe_0 + ts_yyy_xyy[i] * pa_x[i];

        ts_xyyy_xyz[i] = ts_yyy_yz[i] * fe_0 + ts_yyy_xyz[i] * pa_x[i];

        ts_xyyy_xzz[i] = ts_yyy_zz[i] * fe_0 + ts_yyy_xzz[i] * pa_x[i];

        ts_xyyy_yyy[i] = ts_yyy_yyy[i] * pa_x[i];

        ts_xyyy_yyz[i] = ts_yyy_yyz[i] * pa_x[i];

        ts_xyyy_yzz[i] = ts_yyy_yzz[i] * pa_x[i];

        ts_xyyy_zzz[i] = ts_yyy_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto ts_xyyz_xxx = pbuffer.data(idx_ovl_gf + 70);

    auto ts_xyyz_xxy = pbuffer.data(idx_ovl_gf + 71);

    auto ts_xyyz_xxz = pbuffer.data(idx_ovl_gf + 72);

    auto ts_xyyz_xyy = pbuffer.data(idx_ovl_gf + 73);

    auto ts_xyyz_xyz = pbuffer.data(idx_ovl_gf + 74);

    auto ts_xyyz_xzz = pbuffer.data(idx_ovl_gf + 75);

    auto ts_xyyz_yyy = pbuffer.data(idx_ovl_gf + 76);

    auto ts_xyyz_yyz = pbuffer.data(idx_ovl_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_ovl_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_ovl_gf + 79);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             ts_xyy_xxx,  \
                             ts_xyy_xxy,  \
                             ts_xyy_xyy,  \
                             ts_xyyz_xxx, \
                             ts_xyyz_xxy, \
                             ts_xyyz_xxz, \
                             ts_xyyz_xyy, \
                             ts_xyyz_xyz, \
                             ts_xyyz_xzz, \
                             ts_xyyz_yyy, \
                             ts_xyyz_yyz, \
                             ts_xyyz_yzz, \
                             ts_xyyz_zzz, \
                             ts_yyz_xxz,  \
                             ts_yyz_xyz,  \
                             ts_yyz_xz,   \
                             ts_yyz_xzz,  \
                             ts_yyz_yyy,  \
                             ts_yyz_yyz,  \
                             ts_yyz_yz,   \
                             ts_yyz_yzz,  \
                             ts_yyz_zz,   \
                             ts_yyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyz_xxx[i] = ts_xyy_xxx[i] * pa_z[i];

        ts_xyyz_xxy[i] = ts_xyy_xxy[i] * pa_z[i];

        ts_xyyz_xxz[i] = 2.0 * ts_yyz_xz[i] * fe_0 + ts_yyz_xxz[i] * pa_x[i];

        ts_xyyz_xyy[i] = ts_xyy_xyy[i] * pa_z[i];

        ts_xyyz_xyz[i] = ts_yyz_yz[i] * fe_0 + ts_yyz_xyz[i] * pa_x[i];

        ts_xyyz_xzz[i] = ts_yyz_zz[i] * fe_0 + ts_yyz_xzz[i] * pa_x[i];

        ts_xyyz_yyy[i] = ts_yyz_yyy[i] * pa_x[i];

        ts_xyyz_yyz[i] = ts_yyz_yyz[i] * pa_x[i];

        ts_xyyz_yzz[i] = ts_yyz_yzz[i] * pa_x[i];

        ts_xyyz_zzz[i] = ts_yyz_zzz[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto ts_xyzz_xxx = pbuffer.data(idx_ovl_gf + 80);

    auto ts_xyzz_xxy = pbuffer.data(idx_ovl_gf + 81);

    auto ts_xyzz_xxz = pbuffer.data(idx_ovl_gf + 82);

    auto ts_xyzz_xyy = pbuffer.data(idx_ovl_gf + 83);

    auto ts_xyzz_xyz = pbuffer.data(idx_ovl_gf + 84);

    auto ts_xyzz_xzz = pbuffer.data(idx_ovl_gf + 85);

    auto ts_xyzz_yyy = pbuffer.data(idx_ovl_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_ovl_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_ovl_gf + 88);

    auto ts_xyzz_zzz = pbuffer.data(idx_ovl_gf + 89);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             ts_xyzz_xxx, \
                             ts_xyzz_xxy, \
                             ts_xyzz_xxz, \
                             ts_xyzz_xyy, \
                             ts_xyzz_xyz, \
                             ts_xyzz_xzz, \
                             ts_xyzz_yyy, \
                             ts_xyzz_yyz, \
                             ts_xyzz_yzz, \
                             ts_xyzz_zzz, \
                             ts_xzz_xxx,  \
                             ts_xzz_xxz,  \
                             ts_xzz_xzz,  \
                             ts_yzz_xxy,  \
                             ts_yzz_xy,   \
                             ts_yzz_xyy,  \
                             ts_yzz_xyz,  \
                             ts_yzz_yy,   \
                             ts_yzz_yyy,  \
                             ts_yzz_yyz,  \
                             ts_yzz_yz,   \
                             ts_yzz_yzz,  \
                             ts_yzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzz_xxx[i] = ts_xzz_xxx[i] * pa_y[i];

        ts_xyzz_xxy[i] = 2.0 * ts_yzz_xy[i] * fe_0 + ts_yzz_xxy[i] * pa_x[i];

        ts_xyzz_xxz[i] = ts_xzz_xxz[i] * pa_y[i];

        ts_xyzz_xyy[i] = ts_yzz_yy[i] * fe_0 + ts_yzz_xyy[i] * pa_x[i];

        ts_xyzz_xyz[i] = ts_yzz_yz[i] * fe_0 + ts_yzz_xyz[i] * pa_x[i];

        ts_xyzz_xzz[i] = ts_xzz_xzz[i] * pa_y[i];

        ts_xyzz_yyy[i] = ts_yzz_yyy[i] * pa_x[i];

        ts_xyzz_yyz[i] = ts_yzz_yyz[i] * pa_x[i];

        ts_xyzz_yzz[i] = ts_yzz_yzz[i] * pa_x[i];

        ts_xyzz_zzz[i] = ts_yzz_zzz[i] * pa_x[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto ts_xzzz_xxx = pbuffer.data(idx_ovl_gf + 90);

    auto ts_xzzz_xxy = pbuffer.data(idx_ovl_gf + 91);

    auto ts_xzzz_xxz = pbuffer.data(idx_ovl_gf + 92);

    auto ts_xzzz_xyy = pbuffer.data(idx_ovl_gf + 93);

    auto ts_xzzz_xyz = pbuffer.data(idx_ovl_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_ovl_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_ovl_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

#pragma omp simd aligned(pa_x,            \
                             ts_xzzz_xxx, \
                             ts_xzzz_xxy, \
                             ts_xzzz_xxz, \
                             ts_xzzz_xyy, \
                             ts_xzzz_xyz, \
                             ts_xzzz_xzz, \
                             ts_xzzz_yyy, \
                             ts_xzzz_yyz, \
                             ts_xzzz_yzz, \
                             ts_xzzz_zzz, \
                             ts_zzz_xx,   \
                             ts_zzz_xxx,  \
                             ts_zzz_xxy,  \
                             ts_zzz_xxz,  \
                             ts_zzz_xy,   \
                             ts_zzz_xyy,  \
                             ts_zzz_xyz,  \
                             ts_zzz_xz,   \
                             ts_zzz_xzz,  \
                             ts_zzz_yy,   \
                             ts_zzz_yyy,  \
                             ts_zzz_yyz,  \
                             ts_zzz_yz,   \
                             ts_zzz_yzz,  \
                             ts_zzz_zz,   \
                             ts_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzz_xxx[i] = 3.0 * ts_zzz_xx[i] * fe_0 + ts_zzz_xxx[i] * pa_x[i];

        ts_xzzz_xxy[i] = 2.0 * ts_zzz_xy[i] * fe_0 + ts_zzz_xxy[i] * pa_x[i];

        ts_xzzz_xxz[i] = 2.0 * ts_zzz_xz[i] * fe_0 + ts_zzz_xxz[i] * pa_x[i];

        ts_xzzz_xyy[i] = ts_zzz_yy[i] * fe_0 + ts_zzz_xyy[i] * pa_x[i];

        ts_xzzz_xyz[i] = ts_zzz_yz[i] * fe_0 + ts_zzz_xyz[i] * pa_x[i];

        ts_xzzz_xzz[i] = ts_zzz_zz[i] * fe_0 + ts_zzz_xzz[i] * pa_x[i];

        ts_xzzz_yyy[i] = ts_zzz_yyy[i] * pa_x[i];

        ts_xzzz_yyz[i] = ts_zzz_yyz[i] * pa_x[i];

        ts_xzzz_yzz[i] = ts_zzz_yzz[i] * pa_x[i];

        ts_xzzz_zzz[i] = ts_zzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

#pragma omp simd aligned(pa_y,            \
                             ts_yy_xxx,   \
                             ts_yy_xxy,   \
                             ts_yy_xxz,   \
                             ts_yy_xyy,   \
                             ts_yy_xyz,   \
                             ts_yy_xzz,   \
                             ts_yy_yyy,   \
                             ts_yy_yyz,   \
                             ts_yy_yzz,   \
                             ts_yy_zzz,   \
                             ts_yyy_xx,   \
                             ts_yyy_xxx,  \
                             ts_yyy_xxy,  \
                             ts_yyy_xxz,  \
                             ts_yyy_xy,   \
                             ts_yyy_xyy,  \
                             ts_yyy_xyz,  \
                             ts_yyy_xz,   \
                             ts_yyy_xzz,  \
                             ts_yyy_yy,   \
                             ts_yyy_yyy,  \
                             ts_yyy_yyz,  \
                             ts_yyy_yz,   \
                             ts_yyy_yzz,  \
                             ts_yyy_zz,   \
                             ts_yyy_zzz,  \
                             ts_yyyy_xxx, \
                             ts_yyyy_xxy, \
                             ts_yyyy_xxz, \
                             ts_yyyy_xyy, \
                             ts_yyyy_xyz, \
                             ts_yyyy_xzz, \
                             ts_yyyy_yyy, \
                             ts_yyyy_yyz, \
                             ts_yyyy_yzz, \
                             ts_yyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyy_xxx[i] = 3.0 * ts_yy_xxx[i] * fe_0 + ts_yyy_xxx[i] * pa_y[i];

        ts_yyyy_xxy[i] = 3.0 * ts_yy_xxy[i] * fe_0 + ts_yyy_xx[i] * fe_0 + ts_yyy_xxy[i] * pa_y[i];

        ts_yyyy_xxz[i] = 3.0 * ts_yy_xxz[i] * fe_0 + ts_yyy_xxz[i] * pa_y[i];

        ts_yyyy_xyy[i] = 3.0 * ts_yy_xyy[i] * fe_0 + 2.0 * ts_yyy_xy[i] * fe_0 + ts_yyy_xyy[i] * pa_y[i];

        ts_yyyy_xyz[i] = 3.0 * ts_yy_xyz[i] * fe_0 + ts_yyy_xz[i] * fe_0 + ts_yyy_xyz[i] * pa_y[i];

        ts_yyyy_xzz[i] = 3.0 * ts_yy_xzz[i] * fe_0 + ts_yyy_xzz[i] * pa_y[i];

        ts_yyyy_yyy[i] = 3.0 * ts_yy_yyy[i] * fe_0 + 3.0 * ts_yyy_yy[i] * fe_0 + ts_yyy_yyy[i] * pa_y[i];

        ts_yyyy_yyz[i] = 3.0 * ts_yy_yyz[i] * fe_0 + 2.0 * ts_yyy_yz[i] * fe_0 + ts_yyy_yyz[i] * pa_y[i];

        ts_yyyy_yzz[i] = 3.0 * ts_yy_yzz[i] * fe_0 + ts_yyy_zz[i] * fe_0 + ts_yyy_yzz[i] * pa_y[i];

        ts_yyyy_zzz[i] = 3.0 * ts_yy_zzz[i] * fe_0 + ts_yyy_zzz[i] * pa_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto ts_yyyz_xxx = pbuffer.data(idx_ovl_gf + 110);

    auto ts_yyyz_xxy = pbuffer.data(idx_ovl_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_ovl_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_ovl_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_ovl_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_ovl_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_ovl_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_ovl_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_ovl_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_ovl_gf + 119);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             ts_yyy_xxx,  \
                             ts_yyy_xxy,  \
                             ts_yyy_xy,   \
                             ts_yyy_xyy,  \
                             ts_yyy_xyz,  \
                             ts_yyy_yy,   \
                             ts_yyy_yyy,  \
                             ts_yyy_yyz,  \
                             ts_yyy_yz,   \
                             ts_yyy_yzz,  \
                             ts_yyyz_xxx, \
                             ts_yyyz_xxy, \
                             ts_yyyz_xxz, \
                             ts_yyyz_xyy, \
                             ts_yyyz_xyz, \
                             ts_yyyz_xzz, \
                             ts_yyyz_yyy, \
                             ts_yyyz_yyz, \
                             ts_yyyz_yzz, \
                             ts_yyyz_zzz, \
                             ts_yyz_xxz,  \
                             ts_yyz_xzz,  \
                             ts_yyz_zzz,  \
                             ts_yz_xxz,   \
                             ts_yz_xzz,   \
                             ts_yz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyz_xxx[i] = ts_yyy_xxx[i] * pa_z[i];

        ts_yyyz_xxy[i] = ts_yyy_xxy[i] * pa_z[i];

        ts_yyyz_xxz[i] = 2.0 * ts_yz_xxz[i] * fe_0 + ts_yyz_xxz[i] * pa_y[i];

        ts_yyyz_xyy[i] = ts_yyy_xyy[i] * pa_z[i];

        ts_yyyz_xyz[i] = ts_yyy_xy[i] * fe_0 + ts_yyy_xyz[i] * pa_z[i];

        ts_yyyz_xzz[i] = 2.0 * ts_yz_xzz[i] * fe_0 + ts_yyz_xzz[i] * pa_y[i];

        ts_yyyz_yyy[i] = ts_yyy_yyy[i] * pa_z[i];

        ts_yyyz_yyz[i] = ts_yyy_yy[i] * fe_0 + ts_yyy_yyz[i] * pa_z[i];

        ts_yyyz_yzz[i] = 2.0 * ts_yyy_yz[i] * fe_0 + ts_yyy_yzz[i] * pa_z[i];

        ts_yyyz_zzz[i] = 2.0 * ts_yz_zzz[i] * fe_0 + ts_yyz_zzz[i] * pa_y[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto ts_yyzz_xxx = pbuffer.data(idx_ovl_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_ovl_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_ovl_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             ts_yy_xxy,   \
                             ts_yy_xyy,   \
                             ts_yy_yyy,   \
                             ts_yyz_xxy,  \
                             ts_yyz_xyy,  \
                             ts_yyz_yyy,  \
                             ts_yyzz_xxx, \
                             ts_yyzz_xxy, \
                             ts_yyzz_xxz, \
                             ts_yyzz_xyy, \
                             ts_yyzz_xyz, \
                             ts_yyzz_xzz, \
                             ts_yyzz_yyy, \
                             ts_yyzz_yyz, \
                             ts_yyzz_yzz, \
                             ts_yyzz_zzz, \
                             ts_yzz_xxx,  \
                             ts_yzz_xxz,  \
                             ts_yzz_xyz,  \
                             ts_yzz_xz,   \
                             ts_yzz_xzz,  \
                             ts_yzz_yyz,  \
                             ts_yzz_yz,   \
                             ts_yzz_yzz,  \
                             ts_yzz_zz,   \
                             ts_yzz_zzz,  \
                             ts_zz_xxx,   \
                             ts_zz_xxz,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzz_xxx[i] = ts_zz_xxx[i] * fe_0 + ts_yzz_xxx[i] * pa_y[i];

        ts_yyzz_xxy[i] = ts_yy_xxy[i] * fe_0 + ts_yyz_xxy[i] * pa_z[i];

        ts_yyzz_xxz[i] = ts_zz_xxz[i] * fe_0 + ts_yzz_xxz[i] * pa_y[i];

        ts_yyzz_xyy[i] = ts_yy_xyy[i] * fe_0 + ts_yyz_xyy[i] * pa_z[i];

        ts_yyzz_xyz[i] = ts_zz_xyz[i] * fe_0 + ts_yzz_xz[i] * fe_0 + ts_yzz_xyz[i] * pa_y[i];

        ts_yyzz_xzz[i] = ts_zz_xzz[i] * fe_0 + ts_yzz_xzz[i] * pa_y[i];

        ts_yyzz_yyy[i] = ts_yy_yyy[i] * fe_0 + ts_yyz_yyy[i] * pa_z[i];

        ts_yyzz_yyz[i] = ts_zz_yyz[i] * fe_0 + 2.0 * ts_yzz_yz[i] * fe_0 + ts_yzz_yyz[i] * pa_y[i];

        ts_yyzz_yzz[i] = ts_zz_yzz[i] * fe_0 + ts_yzz_zz[i] * fe_0 + ts_yzz_yzz[i] * pa_y[i];

        ts_yyzz_zzz[i] = ts_zz_zzz[i] * fe_0 + ts_yzz_zzz[i] * pa_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto ts_yzzz_xxx = pbuffer.data(idx_ovl_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_ovl_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_ovl_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_ovl_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_ovl_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

#pragma omp simd aligned(pa_y,            \
                             ts_yzzz_xxx, \
                             ts_yzzz_xxy, \
                             ts_yzzz_xxz, \
                             ts_yzzz_xyy, \
                             ts_yzzz_xyz, \
                             ts_yzzz_xzz, \
                             ts_yzzz_yyy, \
                             ts_yzzz_yyz, \
                             ts_yzzz_yzz, \
                             ts_yzzz_zzz, \
                             ts_zzz_xx,   \
                             ts_zzz_xxx,  \
                             ts_zzz_xxy,  \
                             ts_zzz_xxz,  \
                             ts_zzz_xy,   \
                             ts_zzz_xyy,  \
                             ts_zzz_xyz,  \
                             ts_zzz_xz,   \
                             ts_zzz_xzz,  \
                             ts_zzz_yy,   \
                             ts_zzz_yyy,  \
                             ts_zzz_yyz,  \
                             ts_zzz_yz,   \
                             ts_zzz_yzz,  \
                             ts_zzz_zz,   \
                             ts_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzz_xxx[i] = ts_zzz_xxx[i] * pa_y[i];

        ts_yzzz_xxy[i] = ts_zzz_xx[i] * fe_0 + ts_zzz_xxy[i] * pa_y[i];

        ts_yzzz_xxz[i] = ts_zzz_xxz[i] * pa_y[i];

        ts_yzzz_xyy[i] = 2.0 * ts_zzz_xy[i] * fe_0 + ts_zzz_xyy[i] * pa_y[i];

        ts_yzzz_xyz[i] = ts_zzz_xz[i] * fe_0 + ts_zzz_xyz[i] * pa_y[i];

        ts_yzzz_xzz[i] = ts_zzz_xzz[i] * pa_y[i];

        ts_yzzz_yyy[i] = 3.0 * ts_zzz_yy[i] * fe_0 + ts_zzz_yyy[i] * pa_y[i];

        ts_yzzz_yyz[i] = 2.0 * ts_zzz_yz[i] * fe_0 + ts_zzz_yyz[i] * pa_y[i];

        ts_yzzz_yzz[i] = ts_zzz_zz[i] * fe_0 + ts_zzz_yzz[i] * pa_y[i];

        ts_yzzz_zzz[i] = ts_zzz_zzz[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

#pragma omp simd aligned(pa_z,            \
                             ts_zz_xxx,   \
                             ts_zz_xxy,   \
                             ts_zz_xxz,   \
                             ts_zz_xyy,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyy,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             ts_zzz_xx,   \
                             ts_zzz_xxx,  \
                             ts_zzz_xxy,  \
                             ts_zzz_xxz,  \
                             ts_zzz_xy,   \
                             ts_zzz_xyy,  \
                             ts_zzz_xyz,  \
                             ts_zzz_xz,   \
                             ts_zzz_xzz,  \
                             ts_zzz_yy,   \
                             ts_zzz_yyy,  \
                             ts_zzz_yyz,  \
                             ts_zzz_yz,   \
                             ts_zzz_yzz,  \
                             ts_zzz_zz,   \
                             ts_zzz_zzz,  \
                             ts_zzzz_xxx, \
                             ts_zzzz_xxy, \
                             ts_zzzz_xxz, \
                             ts_zzzz_xyy, \
                             ts_zzzz_xyz, \
                             ts_zzzz_xzz, \
                             ts_zzzz_yyy, \
                             ts_zzzz_yyz, \
                             ts_zzzz_yzz, \
                             ts_zzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzz_xxx[i] = 3.0 * ts_zz_xxx[i] * fe_0 + ts_zzz_xxx[i] * pa_z[i];

        ts_zzzz_xxy[i] = 3.0 * ts_zz_xxy[i] * fe_0 + ts_zzz_xxy[i] * pa_z[i];

        ts_zzzz_xxz[i] = 3.0 * ts_zz_xxz[i] * fe_0 + ts_zzz_xx[i] * fe_0 + ts_zzz_xxz[i] * pa_z[i];

        ts_zzzz_xyy[i] = 3.0 * ts_zz_xyy[i] * fe_0 + ts_zzz_xyy[i] * pa_z[i];

        ts_zzzz_xyz[i] = 3.0 * ts_zz_xyz[i] * fe_0 + ts_zzz_xy[i] * fe_0 + ts_zzz_xyz[i] * pa_z[i];

        ts_zzzz_xzz[i] = 3.0 * ts_zz_xzz[i] * fe_0 + 2.0 * ts_zzz_xz[i] * fe_0 + ts_zzz_xzz[i] * pa_z[i];

        ts_zzzz_yyy[i] = 3.0 * ts_zz_yyy[i] * fe_0 + ts_zzz_yyy[i] * pa_z[i];

        ts_zzzz_yyz[i] = 3.0 * ts_zz_yyz[i] * fe_0 + ts_zzz_yy[i] * fe_0 + ts_zzz_yyz[i] * pa_z[i];

        ts_zzzz_yzz[i] = 3.0 * ts_zz_yzz[i] * fe_0 + 2.0 * ts_zzz_yz[i] * fe_0 + ts_zzz_yzz[i] * pa_z[i];

        ts_zzzz_zzz[i] = 3.0 * ts_zz_zzz[i] * fe_0 + 3.0 * ts_zzz_zz[i] * fe_0 + ts_zzz_zzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
