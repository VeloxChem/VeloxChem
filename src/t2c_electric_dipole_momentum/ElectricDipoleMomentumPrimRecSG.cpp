#include "ElectricDipoleMomentumPrimRecSG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_sg(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_sg,
                                      const size_t              idx_dip_sd,
                                      const size_t              idx_ovl_sf,
                                      const size_t              idx_dip_sf,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpb,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SD

    auto tr_x_0_xx = pbuffer.data(idx_dip_sd);

    auto tr_x_0_yy = pbuffer.data(idx_dip_sd + 3);

    auto tr_x_0_zz = pbuffer.data(idx_dip_sd + 5);

    auto tr_y_0_xx = pbuffer.data(idx_dip_sd + 6);

    auto tr_y_0_xy = pbuffer.data(idx_dip_sd + 7);

    auto tr_y_0_yy = pbuffer.data(idx_dip_sd + 9);

    auto tr_y_0_zz = pbuffer.data(idx_dip_sd + 11);

    auto tr_z_0_xx = pbuffer.data(idx_dip_sd + 12);

    auto tr_z_0_xz = pbuffer.data(idx_dip_sd + 14);

    auto tr_z_0_yy = pbuffer.data(idx_dip_sd + 15);

    auto tr_z_0_yz = pbuffer.data(idx_dip_sd + 16);

    auto tr_z_0_zz = pbuffer.data(idx_dip_sd + 17);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto tr_x_0_xxx = pbuffer.data(idx_dip_sf);

    auto tr_x_0_xxy = pbuffer.data(idx_dip_sf + 1);

    auto tr_x_0_xxz = pbuffer.data(idx_dip_sf + 2);

    auto tr_x_0_xyy = pbuffer.data(idx_dip_sf + 3);

    auto tr_x_0_xzz = pbuffer.data(idx_dip_sf + 5);

    auto tr_x_0_yyy = pbuffer.data(idx_dip_sf + 6);

    auto tr_x_0_yzz = pbuffer.data(idx_dip_sf + 8);

    auto tr_x_0_zzz = pbuffer.data(idx_dip_sf + 9);

    auto tr_y_0_xxx = pbuffer.data(idx_dip_sf + 10);

    auto tr_y_0_xxy = pbuffer.data(idx_dip_sf + 11);

    auto tr_y_0_xyy = pbuffer.data(idx_dip_sf + 13);

    auto tr_y_0_xzz = pbuffer.data(idx_dip_sf + 15);

    auto tr_y_0_yyy = pbuffer.data(idx_dip_sf + 16);

    auto tr_y_0_yyz = pbuffer.data(idx_dip_sf + 17);

    auto tr_y_0_yzz = pbuffer.data(idx_dip_sf + 18);

    auto tr_y_0_zzz = pbuffer.data(idx_dip_sf + 19);

    auto tr_z_0_xxx = pbuffer.data(idx_dip_sf + 20);

    auto tr_z_0_xxz = pbuffer.data(idx_dip_sf + 22);

    auto tr_z_0_xyy = pbuffer.data(idx_dip_sf + 23);

    auto tr_z_0_xzz = pbuffer.data(idx_dip_sf + 25);

    auto tr_z_0_yyy = pbuffer.data(idx_dip_sf + 26);

    auto tr_z_0_yyz = pbuffer.data(idx_dip_sf + 27);

    auto tr_z_0_yzz = pbuffer.data(idx_dip_sf + 28);

    auto tr_z_0_zzz = pbuffer.data(idx_dip_sf + 29);

    // Set up components of targeted buffer : SG

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxyz = pbuffer.data(idx_dip_sg + 4);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

    auto tr_x_0_xyyy = pbuffer.data(idx_dip_sg + 6);

    auto tr_x_0_xyyz = pbuffer.data(idx_dip_sg + 7);

    auto tr_x_0_xyzz = pbuffer.data(idx_dip_sg + 8);

    auto tr_x_0_xzzz = pbuffer.data(idx_dip_sg + 9);

    auto tr_x_0_yyyy = pbuffer.data(idx_dip_sg + 10);

    auto tr_x_0_yyyz = pbuffer.data(idx_dip_sg + 11);

    auto tr_x_0_yyzz = pbuffer.data(idx_dip_sg + 12);

    auto tr_x_0_yzzz = pbuffer.data(idx_dip_sg + 13);

    auto tr_x_0_zzzz = pbuffer.data(idx_dip_sg + 14);

    auto tr_y_0_xxxx = pbuffer.data(idx_dip_sg + 15);

    auto tr_y_0_xxxy = pbuffer.data(idx_dip_sg + 16);

    auto tr_y_0_xxxz = pbuffer.data(idx_dip_sg + 17);

    auto tr_y_0_xxyy = pbuffer.data(idx_dip_sg + 18);

    auto tr_y_0_xxyz = pbuffer.data(idx_dip_sg + 19);

    auto tr_y_0_xxzz = pbuffer.data(idx_dip_sg + 20);

    auto tr_y_0_xyyy = pbuffer.data(idx_dip_sg + 21);

    auto tr_y_0_xyyz = pbuffer.data(idx_dip_sg + 22);

    auto tr_y_0_xyzz = pbuffer.data(idx_dip_sg + 23);

    auto tr_y_0_xzzz = pbuffer.data(idx_dip_sg + 24);

    auto tr_y_0_yyyy = pbuffer.data(idx_dip_sg + 25);

    auto tr_y_0_yyyz = pbuffer.data(idx_dip_sg + 26);

    auto tr_y_0_yyzz = pbuffer.data(idx_dip_sg + 27);

    auto tr_y_0_yzzz = pbuffer.data(idx_dip_sg + 28);

    auto tr_y_0_zzzz = pbuffer.data(idx_dip_sg + 29);

    auto tr_z_0_xxxx = pbuffer.data(idx_dip_sg + 30);

    auto tr_z_0_xxxy = pbuffer.data(idx_dip_sg + 31);

    auto tr_z_0_xxxz = pbuffer.data(idx_dip_sg + 32);

    auto tr_z_0_xxyy = pbuffer.data(idx_dip_sg + 33);

    auto tr_z_0_xxyz = pbuffer.data(idx_dip_sg + 34);

    auto tr_z_0_xxzz = pbuffer.data(idx_dip_sg + 35);

    auto tr_z_0_xyyy = pbuffer.data(idx_dip_sg + 36);

    auto tr_z_0_xyyz = pbuffer.data(idx_dip_sg + 37);

    auto tr_z_0_xyzz = pbuffer.data(idx_dip_sg + 38);

    auto tr_z_0_xzzz = pbuffer.data(idx_dip_sg + 39);

    auto tr_z_0_yyyy = pbuffer.data(idx_dip_sg + 40);

    auto tr_z_0_yyyz = pbuffer.data(idx_dip_sg + 41);

    auto tr_z_0_yyzz = pbuffer.data(idx_dip_sg + 42);

    auto tr_z_0_yzzz = pbuffer.data(idx_dip_sg + 43);

    auto tr_z_0_zzzz = pbuffer.data(idx_dip_sg + 44);

#pragma omp simd aligned(pb_x,            \
                             pb_y,        \
                             pb_z,        \
                             tr_x_0_xx,   \
                             tr_x_0_xxx,  \
                             tr_x_0_xxxx, \
                             tr_x_0_xxxy, \
                             tr_x_0_xxxz, \
                             tr_x_0_xxy,  \
                             tr_x_0_xxyy, \
                             tr_x_0_xxyz, \
                             tr_x_0_xxz,  \
                             tr_x_0_xxzz, \
                             tr_x_0_xyy,  \
                             tr_x_0_xyyy, \
                             tr_x_0_xyyz, \
                             tr_x_0_xyzz, \
                             tr_x_0_xzz,  \
                             tr_x_0_xzzz, \
                             tr_x_0_yy,   \
                             tr_x_0_yyy,  \
                             tr_x_0_yyyy, \
                             tr_x_0_yyyz, \
                             tr_x_0_yyzz, \
                             tr_x_0_yzz,  \
                             tr_x_0_yzzz, \
                             tr_x_0_zz,   \
                             tr_x_0_zzz,  \
                             tr_x_0_zzzz, \
                             tr_y_0_xx,   \
                             tr_y_0_xxx,  \
                             tr_y_0_xxxx, \
                             tr_y_0_xxxy, \
                             tr_y_0_xxxz, \
                             tr_y_0_xxy,  \
                             tr_y_0_xxyy, \
                             tr_y_0_xxyz, \
                             tr_y_0_xxzz, \
                             tr_y_0_xy,   \
                             tr_y_0_xyy,  \
                             tr_y_0_xyyy, \
                             tr_y_0_xyyz, \
                             tr_y_0_xyzz, \
                             tr_y_0_xzz,  \
                             tr_y_0_xzzz, \
                             tr_y_0_yy,   \
                             tr_y_0_yyy,  \
                             tr_y_0_yyyy, \
                             tr_y_0_yyyz, \
                             tr_y_0_yyz,  \
                             tr_y_0_yyzz, \
                             tr_y_0_yzz,  \
                             tr_y_0_yzzz, \
                             tr_y_0_zz,   \
                             tr_y_0_zzz,  \
                             tr_y_0_zzzz, \
                             tr_z_0_xx,   \
                             tr_z_0_xxx,  \
                             tr_z_0_xxxx, \
                             tr_z_0_xxxy, \
                             tr_z_0_xxxz, \
                             tr_z_0_xxyy, \
                             tr_z_0_xxyz, \
                             tr_z_0_xxz,  \
                             tr_z_0_xxzz, \
                             tr_z_0_xyy,  \
                             tr_z_0_xyyy, \
                             tr_z_0_xyyz, \
                             tr_z_0_xyzz, \
                             tr_z_0_xz,   \
                             tr_z_0_xzz,  \
                             tr_z_0_xzzz, \
                             tr_z_0_yy,   \
                             tr_z_0_yyy,  \
                             tr_z_0_yyyy, \
                             tr_z_0_yyyz, \
                             tr_z_0_yyz,  \
                             tr_z_0_yyzz, \
                             tr_z_0_yz,   \
                             tr_z_0_yzz,  \
                             tr_z_0_yzzz, \
                             tr_z_0_zz,   \
                             tr_z_0_zzz,  \
                             tr_z_0_zzzz, \
                             ts_0_xxx,    \
                             ts_0_yyy,    \
                             ts_0_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_xxxx[i] = 3.0 * tr_x_0_xx[i] * fe_0 + ts_0_xxx[i] * fe_0 + tr_x_0_xxx[i] * pb_x[i];

        tr_x_0_xxxy[i] = tr_x_0_xxx[i] * pb_y[i];

        tr_x_0_xxxz[i] = tr_x_0_xxx[i] * pb_z[i];

        tr_x_0_xxyy[i] = tr_x_0_xx[i] * fe_0 + tr_x_0_xxy[i] * pb_y[i];

        tr_x_0_xxyz[i] = tr_x_0_xxz[i] * pb_y[i];

        tr_x_0_xxzz[i] = tr_x_0_xx[i] * fe_0 + tr_x_0_xxz[i] * pb_z[i];

        tr_x_0_xyyy[i] = ts_0_yyy[i] * fe_0 + tr_x_0_yyy[i] * pb_x[i];

        tr_x_0_xyyz[i] = tr_x_0_xyy[i] * pb_z[i];

        tr_x_0_xyzz[i] = tr_x_0_xzz[i] * pb_y[i];

        tr_x_0_xzzz[i] = ts_0_zzz[i] * fe_0 + tr_x_0_zzz[i] * pb_x[i];

        tr_x_0_yyyy[i] = 3.0 * tr_x_0_yy[i] * fe_0 + tr_x_0_yyy[i] * pb_y[i];

        tr_x_0_yyyz[i] = tr_x_0_yyy[i] * pb_z[i];

        tr_x_0_yyzz[i] = tr_x_0_zz[i] * fe_0 + tr_x_0_yzz[i] * pb_y[i];

        tr_x_0_yzzz[i] = tr_x_0_zzz[i] * pb_y[i];

        tr_x_0_zzzz[i] = 3.0 * tr_x_0_zz[i] * fe_0 + tr_x_0_zzz[i] * pb_z[i];

        tr_y_0_xxxx[i] = 3.0 * tr_y_0_xx[i] * fe_0 + tr_y_0_xxx[i] * pb_x[i];

        tr_y_0_xxxy[i] = 2.0 * tr_y_0_xy[i] * fe_0 + tr_y_0_xxy[i] * pb_x[i];

        tr_y_0_xxxz[i] = tr_y_0_xxx[i] * pb_z[i];

        tr_y_0_xxyy[i] = tr_y_0_yy[i] * fe_0 + tr_y_0_xyy[i] * pb_x[i];

        tr_y_0_xxyz[i] = tr_y_0_xxy[i] * pb_z[i];

        tr_y_0_xxzz[i] = tr_y_0_zz[i] * fe_0 + tr_y_0_xzz[i] * pb_x[i];

        tr_y_0_xyyy[i] = tr_y_0_yyy[i] * pb_x[i];

        tr_y_0_xyyz[i] = tr_y_0_yyz[i] * pb_x[i];

        tr_y_0_xyzz[i] = tr_y_0_yzz[i] * pb_x[i];

        tr_y_0_xzzz[i] = tr_y_0_zzz[i] * pb_x[i];

        tr_y_0_yyyy[i] = 3.0 * tr_y_0_yy[i] * fe_0 + ts_0_yyy[i] * fe_0 + tr_y_0_yyy[i] * pb_y[i];

        tr_y_0_yyyz[i] = tr_y_0_yyy[i] * pb_z[i];

        tr_y_0_yyzz[i] = tr_y_0_yy[i] * fe_0 + tr_y_0_yyz[i] * pb_z[i];

        tr_y_0_yzzz[i] = ts_0_zzz[i] * fe_0 + tr_y_0_zzz[i] * pb_y[i];

        tr_y_0_zzzz[i] = 3.0 * tr_y_0_zz[i] * fe_0 + tr_y_0_zzz[i] * pb_z[i];

        tr_z_0_xxxx[i] = 3.0 * tr_z_0_xx[i] * fe_0 + tr_z_0_xxx[i] * pb_x[i];

        tr_z_0_xxxy[i] = tr_z_0_xxx[i] * pb_y[i];

        tr_z_0_xxxz[i] = 2.0 * tr_z_0_xz[i] * fe_0 + tr_z_0_xxz[i] * pb_x[i];

        tr_z_0_xxyy[i] = tr_z_0_yy[i] * fe_0 + tr_z_0_xyy[i] * pb_x[i];

        tr_z_0_xxyz[i] = tr_z_0_xxz[i] * pb_y[i];

        tr_z_0_xxzz[i] = tr_z_0_zz[i] * fe_0 + tr_z_0_xzz[i] * pb_x[i];

        tr_z_0_xyyy[i] = tr_z_0_yyy[i] * pb_x[i];

        tr_z_0_xyyz[i] = tr_z_0_yyz[i] * pb_x[i];

        tr_z_0_xyzz[i] = tr_z_0_yzz[i] * pb_x[i];

        tr_z_0_xzzz[i] = tr_z_0_zzz[i] * pb_x[i];

        tr_z_0_yyyy[i] = 3.0 * tr_z_0_yy[i] * fe_0 + tr_z_0_yyy[i] * pb_y[i];

        tr_z_0_yyyz[i] = 2.0 * tr_z_0_yz[i] * fe_0 + tr_z_0_yyz[i] * pb_y[i];

        tr_z_0_yyzz[i] = tr_z_0_zz[i] * fe_0 + tr_z_0_yzz[i] * pb_y[i];

        tr_z_0_yzzz[i] = tr_z_0_zzz[i] * pb_y[i];

        tr_z_0_zzzz[i] = 3.0 * tr_z_0_zz[i] * fe_0 + ts_0_zzz[i] * fe_0 + tr_z_0_zzz[i] * pb_z[i];
    }
}

}  // namespace diprec
