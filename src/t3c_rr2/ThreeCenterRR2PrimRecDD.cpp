#include "ThreeCenterRR2PrimRecDD.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_dd(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_dd,
                  const size_t idx_pd,
                  const size_t idx_g_pd,
                  const size_t idx_dp,
                  const size_t idx_g_dp,
                  const size_t idx_dd,
                  const size_t idx_g_dd,
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

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

    // Set up components of auxiliary buffer : PD

    auto gr_x_xx = pbuffer.data(idx_g_pd);

    auto gr_x_xy = pbuffer.data(idx_g_pd + 1);

    auto gr_x_xz = pbuffer.data(idx_g_pd + 2);

    auto gr_x_yy = pbuffer.data(idx_g_pd + 3);

    auto gr_x_yz = pbuffer.data(idx_g_pd + 4);

    auto gr_x_zz = pbuffer.data(idx_g_pd + 5);

    auto gr_y_xx = pbuffer.data(idx_g_pd + 6);

    auto gr_y_xy = pbuffer.data(idx_g_pd + 7);

    auto gr_y_xz = pbuffer.data(idx_g_pd + 8);

    auto gr_y_yy = pbuffer.data(idx_g_pd + 9);

    auto gr_y_yz = pbuffer.data(idx_g_pd + 10);

    auto gr_y_zz = pbuffer.data(idx_g_pd + 11);

    auto gr_z_xx = pbuffer.data(idx_g_pd + 12);

    auto gr_z_xy = pbuffer.data(idx_g_pd + 13);

    auto gr_z_xz = pbuffer.data(idx_g_pd + 14);

    auto gr_z_yy = pbuffer.data(idx_g_pd + 15);

    auto gr_z_yz = pbuffer.data(idx_g_pd + 16);

    auto gr_z_zz = pbuffer.data(idx_g_pd + 17);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto gr_xx_x = pbuffer.data(idx_g_dp);

    auto gr_xx_y = pbuffer.data(idx_g_dp + 1);

    auto gr_xx_z = pbuffer.data(idx_g_dp + 2);

    auto gr_xy_x = pbuffer.data(idx_g_dp + 3);

    auto gr_xy_y = pbuffer.data(idx_g_dp + 4);

    auto gr_xy_z = pbuffer.data(idx_g_dp + 5);

    auto gr_xz_x = pbuffer.data(idx_g_dp + 6);

    auto gr_xz_y = pbuffer.data(idx_g_dp + 7);

    auto gr_xz_z = pbuffer.data(idx_g_dp + 8);

    auto gr_yy_x = pbuffer.data(idx_g_dp + 9);

    auto gr_yy_y = pbuffer.data(idx_g_dp + 10);

    auto gr_yy_z = pbuffer.data(idx_g_dp + 11);

    auto gr_yz_x = pbuffer.data(idx_g_dp + 12);

    auto gr_yz_y = pbuffer.data(idx_g_dp + 13);

    auto gr_yz_z = pbuffer.data(idx_g_dp + 14);

    auto gr_zz_x = pbuffer.data(idx_g_dp + 15);

    auto gr_zz_y = pbuffer.data(idx_g_dp + 16);

    auto gr_zz_z = pbuffer.data(idx_g_dp + 17);

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

    // Set up components of auxiliary buffer : DD

    auto gr_xx_xx = pbuffer.data(idx_g_dd);

    auto gr_xx_xy = pbuffer.data(idx_g_dd + 1);

    auto gr_xx_xz = pbuffer.data(idx_g_dd + 2);

    auto gr_xx_yy = pbuffer.data(idx_g_dd + 3);

    auto gr_xx_yz = pbuffer.data(idx_g_dd + 4);

    auto gr_xx_zz = pbuffer.data(idx_g_dd + 5);

    auto gr_xy_xx = pbuffer.data(idx_g_dd + 6);

    auto gr_xy_xy = pbuffer.data(idx_g_dd + 7);

    auto gr_xy_xz = pbuffer.data(idx_g_dd + 8);

    auto gr_xy_yy = pbuffer.data(idx_g_dd + 9);

    auto gr_xy_yz = pbuffer.data(idx_g_dd + 10);

    auto gr_xy_zz = pbuffer.data(idx_g_dd + 11);

    auto gr_xz_xx = pbuffer.data(idx_g_dd + 12);

    auto gr_xz_xy = pbuffer.data(idx_g_dd + 13);

    auto gr_xz_xz = pbuffer.data(idx_g_dd + 14);

    auto gr_xz_yy = pbuffer.data(idx_g_dd + 15);

    auto gr_xz_yz = pbuffer.data(idx_g_dd + 16);

    auto gr_xz_zz = pbuffer.data(idx_g_dd + 17);

    auto gr_yy_xx = pbuffer.data(idx_g_dd + 18);

    auto gr_yy_xy = pbuffer.data(idx_g_dd + 19);

    auto gr_yy_xz = pbuffer.data(idx_g_dd + 20);

    auto gr_yy_yy = pbuffer.data(idx_g_dd + 21);

    auto gr_yy_yz = pbuffer.data(idx_g_dd + 22);

    auto gr_yy_zz = pbuffer.data(idx_g_dd + 23);

    auto gr_yz_xx = pbuffer.data(idx_g_dd + 24);

    auto gr_yz_xy = pbuffer.data(idx_g_dd + 25);

    auto gr_yz_xz = pbuffer.data(idx_g_dd + 26);

    auto gr_yz_yy = pbuffer.data(idx_g_dd + 27);

    auto gr_yz_yz = pbuffer.data(idx_g_dd + 28);

    auto gr_yz_zz = pbuffer.data(idx_g_dd + 29);

    auto gr_zz_xx = pbuffer.data(idx_g_dd + 30);

    auto gr_zz_xy = pbuffer.data(idx_g_dd + 31);

    auto gr_zz_xz = pbuffer.data(idx_g_dd + 32);

    auto gr_zz_yy = pbuffer.data(idx_g_dd + 33);

    auto gr_zz_yz = pbuffer.data(idx_g_dd + 34);

    auto gr_zz_zz = pbuffer.data(idx_g_dd + 35);

    // Set up 0-6 components of targeted buffer : DD

    auto grr_x_xx_xx = pbuffer.data(idx_gr_dd);

    auto grr_x_xx_xy = pbuffer.data(idx_gr_dd + 1);

    auto grr_x_xx_xz = pbuffer.data(idx_gr_dd + 2);

    auto grr_x_xx_yy = pbuffer.data(idx_gr_dd + 3);

    auto grr_x_xx_yz = pbuffer.data(idx_gr_dd + 4);

    auto grr_x_xx_zz = pbuffer.data(idx_gr_dd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_yy, gr_x_yz, gr_x_zz, gr_xx_x, gr_xx_xx, gr_xx_xy, gr_xx_xz, gr_xx_y, gr_xx_yy, gr_xx_yz, gr_xx_z, gr_xx_zz, grr_x_xx_xx, grr_x_xx_xy, grr_x_xx_xz, grr_x_xx_yy, grr_x_xx_yz, grr_x_xx_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xx_xx[i] = 2.0 * ts_x_xx[i] * gfe2_0 + 2.0 * gr_x_xx[i] * gfe_0 + 2.0 * ts_xx_x[i] * gfe2_0 + 2.0 * gr_xx_x[i] * gfe_0 + ts_xx_xx[i] * gfe_0 * gc_x[i] + gr_xx_xx[i] * gc_x[i];

        grr_x_xx_xy[i] = 2.0 * ts_x_xy[i] * gfe2_0 + 2.0 * gr_x_xy[i] * gfe_0 + ts_xx_y[i] * gfe2_0 + gr_xx_y[i] * gfe_0 + ts_xx_xy[i] * gfe_0 * gc_x[i] + gr_xx_xy[i] * gc_x[i];

        grr_x_xx_xz[i] = 2.0 * ts_x_xz[i] * gfe2_0 + 2.0 * gr_x_xz[i] * gfe_0 + ts_xx_z[i] * gfe2_0 + gr_xx_z[i] * gfe_0 + ts_xx_xz[i] * gfe_0 * gc_x[i] + gr_xx_xz[i] * gc_x[i];

        grr_x_xx_yy[i] = 2.0 * ts_x_yy[i] * gfe2_0 + 2.0 * gr_x_yy[i] * gfe_0 + ts_xx_yy[i] * gfe_0 * gc_x[i] + gr_xx_yy[i] * gc_x[i];

        grr_x_xx_yz[i] = 2.0 * ts_x_yz[i] * gfe2_0 + 2.0 * gr_x_yz[i] * gfe_0 + ts_xx_yz[i] * gfe_0 * gc_x[i] + gr_xx_yz[i] * gc_x[i];

        grr_x_xx_zz[i] = 2.0 * ts_x_zz[i] * gfe2_0 + 2.0 * gr_x_zz[i] * gfe_0 + ts_xx_zz[i] * gfe_0 * gc_x[i] + gr_xx_zz[i] * gc_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto grr_x_xy_xx = pbuffer.data(idx_gr_dd + 6);

    auto grr_x_xy_xy = pbuffer.data(idx_gr_dd + 7);

    auto grr_x_xy_xz = pbuffer.data(idx_gr_dd + 8);

    auto grr_x_xy_yy = pbuffer.data(idx_gr_dd + 9);

    auto grr_x_xy_yz = pbuffer.data(idx_gr_dd + 10);

    auto grr_x_xy_zz = pbuffer.data(idx_gr_dd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_x, gr_xy_xx, gr_xy_xy, gr_xy_xz, gr_xy_y, gr_xy_yy, gr_xy_yz, gr_xy_z, gr_xy_zz, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_yy, gr_y_yz, gr_y_zz, grr_x_xy_xx, grr_x_xy_xy, grr_x_xy_xz, grr_x_xy_yy, grr_x_xy_yz, grr_x_xy_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xy_xx[i] = ts_y_xx[i] * gfe2_0 + gr_y_xx[i] * gfe_0 + 2.0 * ts_xy_x[i] * gfe2_0 + 2.0 * gr_xy_x[i] * gfe_0 + ts_xy_xx[i] * gfe_0 * gc_x[i] + gr_xy_xx[i] * gc_x[i];

        grr_x_xy_xy[i] = ts_y_xy[i] * gfe2_0 + gr_y_xy[i] * gfe_0 + ts_xy_y[i] * gfe2_0 + gr_xy_y[i] * gfe_0 + ts_xy_xy[i] * gfe_0 * gc_x[i] + gr_xy_xy[i] * gc_x[i];

        grr_x_xy_xz[i] = ts_y_xz[i] * gfe2_0 + gr_y_xz[i] * gfe_0 + ts_xy_z[i] * gfe2_0 + gr_xy_z[i] * gfe_0 + ts_xy_xz[i] * gfe_0 * gc_x[i] + gr_xy_xz[i] * gc_x[i];

        grr_x_xy_yy[i] = ts_y_yy[i] * gfe2_0 + gr_y_yy[i] * gfe_0 + ts_xy_yy[i] * gfe_0 * gc_x[i] + gr_xy_yy[i] * gc_x[i];

        grr_x_xy_yz[i] = ts_y_yz[i] * gfe2_0 + gr_y_yz[i] * gfe_0 + ts_xy_yz[i] * gfe_0 * gc_x[i] + gr_xy_yz[i] * gc_x[i];

        grr_x_xy_zz[i] = ts_y_zz[i] * gfe2_0 + gr_y_zz[i] * gfe_0 + ts_xy_zz[i] * gfe_0 * gc_x[i] + gr_xy_zz[i] * gc_x[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto grr_x_xz_xx = pbuffer.data(idx_gr_dd + 12);

    auto grr_x_xz_xy = pbuffer.data(idx_gr_dd + 13);

    auto grr_x_xz_xz = pbuffer.data(idx_gr_dd + 14);

    auto grr_x_xz_yy = pbuffer.data(idx_gr_dd + 15);

    auto grr_x_xz_yz = pbuffer.data(idx_gr_dd + 16);

    auto grr_x_xz_zz = pbuffer.data(idx_gr_dd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_x, gr_xz_xx, gr_xz_xy, gr_xz_xz, gr_xz_y, gr_xz_yy, gr_xz_yz, gr_xz_z, gr_xz_zz, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_yy, gr_z_yz, gr_z_zz, grr_x_xz_xx, grr_x_xz_xy, grr_x_xz_xz, grr_x_xz_yy, grr_x_xz_yz, grr_x_xz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xz_xx[i] = ts_z_xx[i] * gfe2_0 + gr_z_xx[i] * gfe_0 + 2.0 * ts_xz_x[i] * gfe2_0 + 2.0 * gr_xz_x[i] * gfe_0 + ts_xz_xx[i] * gfe_0 * gc_x[i] + gr_xz_xx[i] * gc_x[i];

        grr_x_xz_xy[i] = ts_z_xy[i] * gfe2_0 + gr_z_xy[i] * gfe_0 + ts_xz_y[i] * gfe2_0 + gr_xz_y[i] * gfe_0 + ts_xz_xy[i] * gfe_0 * gc_x[i] + gr_xz_xy[i] * gc_x[i];

        grr_x_xz_xz[i] = ts_z_xz[i] * gfe2_0 + gr_z_xz[i] * gfe_0 + ts_xz_z[i] * gfe2_0 + gr_xz_z[i] * gfe_0 + ts_xz_xz[i] * gfe_0 * gc_x[i] + gr_xz_xz[i] * gc_x[i];

        grr_x_xz_yy[i] = ts_z_yy[i] * gfe2_0 + gr_z_yy[i] * gfe_0 + ts_xz_yy[i] * gfe_0 * gc_x[i] + gr_xz_yy[i] * gc_x[i];

        grr_x_xz_yz[i] = ts_z_yz[i] * gfe2_0 + gr_z_yz[i] * gfe_0 + ts_xz_yz[i] * gfe_0 * gc_x[i] + gr_xz_yz[i] * gc_x[i];

        grr_x_xz_zz[i] = ts_z_zz[i] * gfe2_0 + gr_z_zz[i] * gfe_0 + ts_xz_zz[i] * gfe_0 * gc_x[i] + gr_xz_zz[i] * gc_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto grr_x_yy_xx = pbuffer.data(idx_gr_dd + 18);

    auto grr_x_yy_xy = pbuffer.data(idx_gr_dd + 19);

    auto grr_x_yy_xz = pbuffer.data(idx_gr_dd + 20);

    auto grr_x_yy_yy = pbuffer.data(idx_gr_dd + 21);

    auto grr_x_yy_yz = pbuffer.data(idx_gr_dd + 22);

    auto grr_x_yy_zz = pbuffer.data(idx_gr_dd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_x, gr_yy_xx, gr_yy_xy, gr_yy_xz, gr_yy_y, gr_yy_yy, gr_yy_yz, gr_yy_z, gr_yy_zz, grr_x_yy_xx, grr_x_yy_xy, grr_x_yy_xz, grr_x_yy_yy, grr_x_yy_yz, grr_x_yy_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yy_xx[i] = 2.0 * ts_yy_x[i] * gfe2_0 + 2.0 * gr_yy_x[i] * gfe_0 + ts_yy_xx[i] * gfe_0 * gc_x[i] + gr_yy_xx[i] * gc_x[i];

        grr_x_yy_xy[i] = ts_yy_y[i] * gfe2_0 + gr_yy_y[i] * gfe_0 + ts_yy_xy[i] * gfe_0 * gc_x[i] + gr_yy_xy[i] * gc_x[i];

        grr_x_yy_xz[i] = ts_yy_z[i] * gfe2_0 + gr_yy_z[i] * gfe_0 + ts_yy_xz[i] * gfe_0 * gc_x[i] + gr_yy_xz[i] * gc_x[i];

        grr_x_yy_yy[i] = ts_yy_yy[i] * gfe_0 * gc_x[i] + gr_yy_yy[i] * gc_x[i];

        grr_x_yy_yz[i] = ts_yy_yz[i] * gfe_0 * gc_x[i] + gr_yy_yz[i] * gc_x[i];

        grr_x_yy_zz[i] = ts_yy_zz[i] * gfe_0 * gc_x[i] + gr_yy_zz[i] * gc_x[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto grr_x_yz_xx = pbuffer.data(idx_gr_dd + 24);

    auto grr_x_yz_xy = pbuffer.data(idx_gr_dd + 25);

    auto grr_x_yz_xz = pbuffer.data(idx_gr_dd + 26);

    auto grr_x_yz_yy = pbuffer.data(idx_gr_dd + 27);

    auto grr_x_yz_yz = pbuffer.data(idx_gr_dd + 28);

    auto grr_x_yz_zz = pbuffer.data(idx_gr_dd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_x, gr_yz_xx, gr_yz_xy, gr_yz_xz, gr_yz_y, gr_yz_yy, gr_yz_yz, gr_yz_z, gr_yz_zz, grr_x_yz_xx, grr_x_yz_xy, grr_x_yz_xz, grr_x_yz_yy, grr_x_yz_yz, grr_x_yz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yz_xx[i] = 2.0 * ts_yz_x[i] * gfe2_0 + 2.0 * gr_yz_x[i] * gfe_0 + ts_yz_xx[i] * gfe_0 * gc_x[i] + gr_yz_xx[i] * gc_x[i];

        grr_x_yz_xy[i] = ts_yz_y[i] * gfe2_0 + gr_yz_y[i] * gfe_0 + ts_yz_xy[i] * gfe_0 * gc_x[i] + gr_yz_xy[i] * gc_x[i];

        grr_x_yz_xz[i] = ts_yz_z[i] * gfe2_0 + gr_yz_z[i] * gfe_0 + ts_yz_xz[i] * gfe_0 * gc_x[i] + gr_yz_xz[i] * gc_x[i];

        grr_x_yz_yy[i] = ts_yz_yy[i] * gfe_0 * gc_x[i] + gr_yz_yy[i] * gc_x[i];

        grr_x_yz_yz[i] = ts_yz_yz[i] * gfe_0 * gc_x[i] + gr_yz_yz[i] * gc_x[i];

        grr_x_yz_zz[i] = ts_yz_zz[i] * gfe_0 * gc_x[i] + gr_yz_zz[i] * gc_x[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto grr_x_zz_xx = pbuffer.data(idx_gr_dd + 30);

    auto grr_x_zz_xy = pbuffer.data(idx_gr_dd + 31);

    auto grr_x_zz_xz = pbuffer.data(idx_gr_dd + 32);

    auto grr_x_zz_yy = pbuffer.data(idx_gr_dd + 33);

    auto grr_x_zz_yz = pbuffer.data(idx_gr_dd + 34);

    auto grr_x_zz_zz = pbuffer.data(idx_gr_dd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_x, gr_zz_xx, gr_zz_xy, gr_zz_xz, gr_zz_y, gr_zz_yy, gr_zz_yz, gr_zz_z, gr_zz_zz, grr_x_zz_xx, grr_x_zz_xy, grr_x_zz_xz, grr_x_zz_yy, grr_x_zz_yz, grr_x_zz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zz_xx[i] = 2.0 * ts_zz_x[i] * gfe2_0 + 2.0 * gr_zz_x[i] * gfe_0 + ts_zz_xx[i] * gfe_0 * gc_x[i] + gr_zz_xx[i] * gc_x[i];

        grr_x_zz_xy[i] = ts_zz_y[i] * gfe2_0 + gr_zz_y[i] * gfe_0 + ts_zz_xy[i] * gfe_0 * gc_x[i] + gr_zz_xy[i] * gc_x[i];

        grr_x_zz_xz[i] = ts_zz_z[i] * gfe2_0 + gr_zz_z[i] * gfe_0 + ts_zz_xz[i] * gfe_0 * gc_x[i] + gr_zz_xz[i] * gc_x[i];

        grr_x_zz_yy[i] = ts_zz_yy[i] * gfe_0 * gc_x[i] + gr_zz_yy[i] * gc_x[i];

        grr_x_zz_yz[i] = ts_zz_yz[i] * gfe_0 * gc_x[i] + gr_zz_yz[i] * gc_x[i];

        grr_x_zz_zz[i] = ts_zz_zz[i] * gfe_0 * gc_x[i] + gr_zz_zz[i] * gc_x[i];
    }

    // Set up 36-42 components of targeted buffer : DD

    auto grr_y_xx_xx = pbuffer.data(idx_gr_dd + 36);

    auto grr_y_xx_xy = pbuffer.data(idx_gr_dd + 37);

    auto grr_y_xx_xz = pbuffer.data(idx_gr_dd + 38);

    auto grr_y_xx_yy = pbuffer.data(idx_gr_dd + 39);

    auto grr_y_xx_yz = pbuffer.data(idx_gr_dd + 40);

    auto grr_y_xx_zz = pbuffer.data(idx_gr_dd + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_xx, gr_xx_xy, gr_xx_xz, gr_xx_y, gr_xx_yy, gr_xx_yz, gr_xx_z, gr_xx_zz, grr_y_xx_xx, grr_y_xx_xy, grr_y_xx_xz, grr_y_xx_yy, grr_y_xx_yz, grr_y_xx_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xx_xx[i] = ts_xx_xx[i] * gfe_0 * gc_y[i] + gr_xx_xx[i] * gc_y[i];

        grr_y_xx_xy[i] = ts_xx_x[i] * gfe2_0 + gr_xx_x[i] * gfe_0 + ts_xx_xy[i] * gfe_0 * gc_y[i] + gr_xx_xy[i] * gc_y[i];

        grr_y_xx_xz[i] = ts_xx_xz[i] * gfe_0 * gc_y[i] + gr_xx_xz[i] * gc_y[i];

        grr_y_xx_yy[i] = 2.0 * ts_xx_y[i] * gfe2_0 + 2.0 * gr_xx_y[i] * gfe_0 + ts_xx_yy[i] * gfe_0 * gc_y[i] + gr_xx_yy[i] * gc_y[i];

        grr_y_xx_yz[i] = ts_xx_z[i] * gfe2_0 + gr_xx_z[i] * gfe_0 + ts_xx_yz[i] * gfe_0 * gc_y[i] + gr_xx_yz[i] * gc_y[i];

        grr_y_xx_zz[i] = ts_xx_zz[i] * gfe_0 * gc_y[i] + gr_xx_zz[i] * gc_y[i];
    }

    // Set up 42-48 components of targeted buffer : DD

    auto grr_y_xy_xx = pbuffer.data(idx_gr_dd + 42);

    auto grr_y_xy_xy = pbuffer.data(idx_gr_dd + 43);

    auto grr_y_xy_xz = pbuffer.data(idx_gr_dd + 44);

    auto grr_y_xy_yy = pbuffer.data(idx_gr_dd + 45);

    auto grr_y_xy_yz = pbuffer.data(idx_gr_dd + 46);

    auto grr_y_xy_zz = pbuffer.data(idx_gr_dd + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_yy, gr_x_yz, gr_x_zz, gr_xy_x, gr_xy_xx, gr_xy_xy, gr_xy_xz, gr_xy_y, gr_xy_yy, gr_xy_yz, gr_xy_z, gr_xy_zz, grr_y_xy_xx, grr_y_xy_xy, grr_y_xy_xz, grr_y_xy_yy, grr_y_xy_yz, grr_y_xy_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xy_xx[i] = ts_x_xx[i] * gfe2_0 + gr_x_xx[i] * gfe_0 + ts_xy_xx[i] * gfe_0 * gc_y[i] + gr_xy_xx[i] * gc_y[i];

        grr_y_xy_xy[i] = ts_x_xy[i] * gfe2_0 + gr_x_xy[i] * gfe_0 + ts_xy_x[i] * gfe2_0 + gr_xy_x[i] * gfe_0 + ts_xy_xy[i] * gfe_0 * gc_y[i] + gr_xy_xy[i] * gc_y[i];

        grr_y_xy_xz[i] = ts_x_xz[i] * gfe2_0 + gr_x_xz[i] * gfe_0 + ts_xy_xz[i] * gfe_0 * gc_y[i] + gr_xy_xz[i] * gc_y[i];

        grr_y_xy_yy[i] = ts_x_yy[i] * gfe2_0 + gr_x_yy[i] * gfe_0 + 2.0 * ts_xy_y[i] * gfe2_0 + 2.0 * gr_xy_y[i] * gfe_0 + ts_xy_yy[i] * gfe_0 * gc_y[i] + gr_xy_yy[i] * gc_y[i];

        grr_y_xy_yz[i] = ts_x_yz[i] * gfe2_0 + gr_x_yz[i] * gfe_0 + ts_xy_z[i] * gfe2_0 + gr_xy_z[i] * gfe_0 + ts_xy_yz[i] * gfe_0 * gc_y[i] + gr_xy_yz[i] * gc_y[i];

        grr_y_xy_zz[i] = ts_x_zz[i] * gfe2_0 + gr_x_zz[i] * gfe_0 + ts_xy_zz[i] * gfe_0 * gc_y[i] + gr_xy_zz[i] * gc_y[i];
    }

    // Set up 48-54 components of targeted buffer : DD

    auto grr_y_xz_xx = pbuffer.data(idx_gr_dd + 48);

    auto grr_y_xz_xy = pbuffer.data(idx_gr_dd + 49);

    auto grr_y_xz_xz = pbuffer.data(idx_gr_dd + 50);

    auto grr_y_xz_yy = pbuffer.data(idx_gr_dd + 51);

    auto grr_y_xz_yz = pbuffer.data(idx_gr_dd + 52);

    auto grr_y_xz_zz = pbuffer.data(idx_gr_dd + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_x, gr_xz_xx, gr_xz_xy, gr_xz_xz, gr_xz_y, gr_xz_yy, gr_xz_yz, gr_xz_z, gr_xz_zz, grr_y_xz_xx, grr_y_xz_xy, grr_y_xz_xz, grr_y_xz_yy, grr_y_xz_yz, grr_y_xz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xz_xx[i] = ts_xz_xx[i] * gfe_0 * gc_y[i] + gr_xz_xx[i] * gc_y[i];

        grr_y_xz_xy[i] = ts_xz_x[i] * gfe2_0 + gr_xz_x[i] * gfe_0 + ts_xz_xy[i] * gfe_0 * gc_y[i] + gr_xz_xy[i] * gc_y[i];

        grr_y_xz_xz[i] = ts_xz_xz[i] * gfe_0 * gc_y[i] + gr_xz_xz[i] * gc_y[i];

        grr_y_xz_yy[i] = 2.0 * ts_xz_y[i] * gfe2_0 + 2.0 * gr_xz_y[i] * gfe_0 + ts_xz_yy[i] * gfe_0 * gc_y[i] + gr_xz_yy[i] * gc_y[i];

        grr_y_xz_yz[i] = ts_xz_z[i] * gfe2_0 + gr_xz_z[i] * gfe_0 + ts_xz_yz[i] * gfe_0 * gc_y[i] + gr_xz_yz[i] * gc_y[i];

        grr_y_xz_zz[i] = ts_xz_zz[i] * gfe_0 * gc_y[i] + gr_xz_zz[i] * gc_y[i];
    }

    // Set up 54-60 components of targeted buffer : DD

    auto grr_y_yy_xx = pbuffer.data(idx_gr_dd + 54);

    auto grr_y_yy_xy = pbuffer.data(idx_gr_dd + 55);

    auto grr_y_yy_xz = pbuffer.data(idx_gr_dd + 56);

    auto grr_y_yy_yy = pbuffer.data(idx_gr_dd + 57);

    auto grr_y_yy_yz = pbuffer.data(idx_gr_dd + 58);

    auto grr_y_yy_zz = pbuffer.data(idx_gr_dd + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_yy, gr_y_yz, gr_y_zz, gr_yy_x, gr_yy_xx, gr_yy_xy, gr_yy_xz, gr_yy_y, gr_yy_yy, gr_yy_yz, gr_yy_z, gr_yy_zz, grr_y_yy_xx, grr_y_yy_xy, grr_y_yy_xz, grr_y_yy_yy, grr_y_yy_yz, grr_y_yy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yy_xx[i] = 2.0 * ts_y_xx[i] * gfe2_0 + 2.0 * gr_y_xx[i] * gfe_0 + ts_yy_xx[i] * gfe_0 * gc_y[i] + gr_yy_xx[i] * gc_y[i];

        grr_y_yy_xy[i] = 2.0 * ts_y_xy[i] * gfe2_0 + 2.0 * gr_y_xy[i] * gfe_0 + ts_yy_x[i] * gfe2_0 + gr_yy_x[i] * gfe_0 + ts_yy_xy[i] * gfe_0 * gc_y[i] + gr_yy_xy[i] * gc_y[i];

        grr_y_yy_xz[i] = 2.0 * ts_y_xz[i] * gfe2_0 + 2.0 * gr_y_xz[i] * gfe_0 + ts_yy_xz[i] * gfe_0 * gc_y[i] + gr_yy_xz[i] * gc_y[i];

        grr_y_yy_yy[i] = 2.0 * ts_y_yy[i] * gfe2_0 + 2.0 * gr_y_yy[i] * gfe_0 + 2.0 * ts_yy_y[i] * gfe2_0 + 2.0 * gr_yy_y[i] * gfe_0 + ts_yy_yy[i] * gfe_0 * gc_y[i] + gr_yy_yy[i] * gc_y[i];

        grr_y_yy_yz[i] = 2.0 * ts_y_yz[i] * gfe2_0 + 2.0 * gr_y_yz[i] * gfe_0 + ts_yy_z[i] * gfe2_0 + gr_yy_z[i] * gfe_0 + ts_yy_yz[i] * gfe_0 * gc_y[i] + gr_yy_yz[i] * gc_y[i];

        grr_y_yy_zz[i] = 2.0 * ts_y_zz[i] * gfe2_0 + 2.0 * gr_y_zz[i] * gfe_0 + ts_yy_zz[i] * gfe_0 * gc_y[i] + gr_yy_zz[i] * gc_y[i];
    }

    // Set up 60-66 components of targeted buffer : DD

    auto grr_y_yz_xx = pbuffer.data(idx_gr_dd + 60);

    auto grr_y_yz_xy = pbuffer.data(idx_gr_dd + 61);

    auto grr_y_yz_xz = pbuffer.data(idx_gr_dd + 62);

    auto grr_y_yz_yy = pbuffer.data(idx_gr_dd + 63);

    auto grr_y_yz_yz = pbuffer.data(idx_gr_dd + 64);

    auto grr_y_yz_zz = pbuffer.data(idx_gr_dd + 65);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_x, gr_yz_xx, gr_yz_xy, gr_yz_xz, gr_yz_y, gr_yz_yy, gr_yz_yz, gr_yz_z, gr_yz_zz, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_yy, gr_z_yz, gr_z_zz, grr_y_yz_xx, grr_y_yz_xy, grr_y_yz_xz, grr_y_yz_yy, grr_y_yz_yz, grr_y_yz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yz_xx[i] = ts_z_xx[i] * gfe2_0 + gr_z_xx[i] * gfe_0 + ts_yz_xx[i] * gfe_0 * gc_y[i] + gr_yz_xx[i] * gc_y[i];

        grr_y_yz_xy[i] = ts_z_xy[i] * gfe2_0 + gr_z_xy[i] * gfe_0 + ts_yz_x[i] * gfe2_0 + gr_yz_x[i] * gfe_0 + ts_yz_xy[i] * gfe_0 * gc_y[i] + gr_yz_xy[i] * gc_y[i];

        grr_y_yz_xz[i] = ts_z_xz[i] * gfe2_0 + gr_z_xz[i] * gfe_0 + ts_yz_xz[i] * gfe_0 * gc_y[i] + gr_yz_xz[i] * gc_y[i];

        grr_y_yz_yy[i] = ts_z_yy[i] * gfe2_0 + gr_z_yy[i] * gfe_0 + 2.0 * ts_yz_y[i] * gfe2_0 + 2.0 * gr_yz_y[i] * gfe_0 + ts_yz_yy[i] * gfe_0 * gc_y[i] + gr_yz_yy[i] * gc_y[i];

        grr_y_yz_yz[i] = ts_z_yz[i] * gfe2_0 + gr_z_yz[i] * gfe_0 + ts_yz_z[i] * gfe2_0 + gr_yz_z[i] * gfe_0 + ts_yz_yz[i] * gfe_0 * gc_y[i] + gr_yz_yz[i] * gc_y[i];

        grr_y_yz_zz[i] = ts_z_zz[i] * gfe2_0 + gr_z_zz[i] * gfe_0 + ts_yz_zz[i] * gfe_0 * gc_y[i] + gr_yz_zz[i] * gc_y[i];
    }

    // Set up 66-72 components of targeted buffer : DD

    auto grr_y_zz_xx = pbuffer.data(idx_gr_dd + 66);

    auto grr_y_zz_xy = pbuffer.data(idx_gr_dd + 67);

    auto grr_y_zz_xz = pbuffer.data(idx_gr_dd + 68);

    auto grr_y_zz_yy = pbuffer.data(idx_gr_dd + 69);

    auto grr_y_zz_yz = pbuffer.data(idx_gr_dd + 70);

    auto grr_y_zz_zz = pbuffer.data(idx_gr_dd + 71);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_x, gr_zz_xx, gr_zz_xy, gr_zz_xz, gr_zz_y, gr_zz_yy, gr_zz_yz, gr_zz_z, gr_zz_zz, grr_y_zz_xx, grr_y_zz_xy, grr_y_zz_xz, grr_y_zz_yy, grr_y_zz_yz, grr_y_zz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zz_xx[i] = ts_zz_xx[i] * gfe_0 * gc_y[i] + gr_zz_xx[i] * gc_y[i];

        grr_y_zz_xy[i] = ts_zz_x[i] * gfe2_0 + gr_zz_x[i] * gfe_0 + ts_zz_xy[i] * gfe_0 * gc_y[i] + gr_zz_xy[i] * gc_y[i];

        grr_y_zz_xz[i] = ts_zz_xz[i] * gfe_0 * gc_y[i] + gr_zz_xz[i] * gc_y[i];

        grr_y_zz_yy[i] = 2.0 * ts_zz_y[i] * gfe2_0 + 2.0 * gr_zz_y[i] * gfe_0 + ts_zz_yy[i] * gfe_0 * gc_y[i] + gr_zz_yy[i] * gc_y[i];

        grr_y_zz_yz[i] = ts_zz_z[i] * gfe2_0 + gr_zz_z[i] * gfe_0 + ts_zz_yz[i] * gfe_0 * gc_y[i] + gr_zz_yz[i] * gc_y[i];

        grr_y_zz_zz[i] = ts_zz_zz[i] * gfe_0 * gc_y[i] + gr_zz_zz[i] * gc_y[i];
    }

    // Set up 72-78 components of targeted buffer : DD

    auto grr_z_xx_xx = pbuffer.data(idx_gr_dd + 72);

    auto grr_z_xx_xy = pbuffer.data(idx_gr_dd + 73);

    auto grr_z_xx_xz = pbuffer.data(idx_gr_dd + 74);

    auto grr_z_xx_yy = pbuffer.data(idx_gr_dd + 75);

    auto grr_z_xx_yz = pbuffer.data(idx_gr_dd + 76);

    auto grr_z_xx_zz = pbuffer.data(idx_gr_dd + 77);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_xx, gr_xx_xy, gr_xx_xz, gr_xx_y, gr_xx_yy, gr_xx_yz, gr_xx_z, gr_xx_zz, grr_z_xx_xx, grr_z_xx_xy, grr_z_xx_xz, grr_z_xx_yy, grr_z_xx_yz, grr_z_xx_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xx_xx[i] = ts_xx_xx[i] * gfe_0 * gc_z[i] + gr_xx_xx[i] * gc_z[i];

        grr_z_xx_xy[i] = ts_xx_xy[i] * gfe_0 * gc_z[i] + gr_xx_xy[i] * gc_z[i];

        grr_z_xx_xz[i] = ts_xx_x[i] * gfe2_0 + gr_xx_x[i] * gfe_0 + ts_xx_xz[i] * gfe_0 * gc_z[i] + gr_xx_xz[i] * gc_z[i];

        grr_z_xx_yy[i] = ts_xx_yy[i] * gfe_0 * gc_z[i] + gr_xx_yy[i] * gc_z[i];

        grr_z_xx_yz[i] = ts_xx_y[i] * gfe2_0 + gr_xx_y[i] * gfe_0 + ts_xx_yz[i] * gfe_0 * gc_z[i] + gr_xx_yz[i] * gc_z[i];

        grr_z_xx_zz[i] = 2.0 * ts_xx_z[i] * gfe2_0 + 2.0 * gr_xx_z[i] * gfe_0 + ts_xx_zz[i] * gfe_0 * gc_z[i] + gr_xx_zz[i] * gc_z[i];
    }

    // Set up 78-84 components of targeted buffer : DD

    auto grr_z_xy_xx = pbuffer.data(idx_gr_dd + 78);

    auto grr_z_xy_xy = pbuffer.data(idx_gr_dd + 79);

    auto grr_z_xy_xz = pbuffer.data(idx_gr_dd + 80);

    auto grr_z_xy_yy = pbuffer.data(idx_gr_dd + 81);

    auto grr_z_xy_yz = pbuffer.data(idx_gr_dd + 82);

    auto grr_z_xy_zz = pbuffer.data(idx_gr_dd + 83);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_x, gr_xy_xx, gr_xy_xy, gr_xy_xz, gr_xy_y, gr_xy_yy, gr_xy_yz, gr_xy_z, gr_xy_zz, grr_z_xy_xx, grr_z_xy_xy, grr_z_xy_xz, grr_z_xy_yy, grr_z_xy_yz, grr_z_xy_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xy_xx[i] = ts_xy_xx[i] * gfe_0 * gc_z[i] + gr_xy_xx[i] * gc_z[i];

        grr_z_xy_xy[i] = ts_xy_xy[i] * gfe_0 * gc_z[i] + gr_xy_xy[i] * gc_z[i];

        grr_z_xy_xz[i] = ts_xy_x[i] * gfe2_0 + gr_xy_x[i] * gfe_0 + ts_xy_xz[i] * gfe_0 * gc_z[i] + gr_xy_xz[i] * gc_z[i];

        grr_z_xy_yy[i] = ts_xy_yy[i] * gfe_0 * gc_z[i] + gr_xy_yy[i] * gc_z[i];

        grr_z_xy_yz[i] = ts_xy_y[i] * gfe2_0 + gr_xy_y[i] * gfe_0 + ts_xy_yz[i] * gfe_0 * gc_z[i] + gr_xy_yz[i] * gc_z[i];

        grr_z_xy_zz[i] = 2.0 * ts_xy_z[i] * gfe2_0 + 2.0 * gr_xy_z[i] * gfe_0 + ts_xy_zz[i] * gfe_0 * gc_z[i] + gr_xy_zz[i] * gc_z[i];
    }

    // Set up 84-90 components of targeted buffer : DD

    auto grr_z_xz_xx = pbuffer.data(idx_gr_dd + 84);

    auto grr_z_xz_xy = pbuffer.data(idx_gr_dd + 85);

    auto grr_z_xz_xz = pbuffer.data(idx_gr_dd + 86);

    auto grr_z_xz_yy = pbuffer.data(idx_gr_dd + 87);

    auto grr_z_xz_yz = pbuffer.data(idx_gr_dd + 88);

    auto grr_z_xz_zz = pbuffer.data(idx_gr_dd + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_yy, gr_x_yz, gr_x_zz, gr_xz_x, gr_xz_xx, gr_xz_xy, gr_xz_xz, gr_xz_y, gr_xz_yy, gr_xz_yz, gr_xz_z, gr_xz_zz, grr_z_xz_xx, grr_z_xz_xy, grr_z_xz_xz, grr_z_xz_yy, grr_z_xz_yz, grr_z_xz_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xz_xx[i] = ts_x_xx[i] * gfe2_0 + gr_x_xx[i] * gfe_0 + ts_xz_xx[i] * gfe_0 * gc_z[i] + gr_xz_xx[i] * gc_z[i];

        grr_z_xz_xy[i] = ts_x_xy[i] * gfe2_0 + gr_x_xy[i] * gfe_0 + ts_xz_xy[i] * gfe_0 * gc_z[i] + gr_xz_xy[i] * gc_z[i];

        grr_z_xz_xz[i] = ts_x_xz[i] * gfe2_0 + gr_x_xz[i] * gfe_0 + ts_xz_x[i] * gfe2_0 + gr_xz_x[i] * gfe_0 + ts_xz_xz[i] * gfe_0 * gc_z[i] + gr_xz_xz[i] * gc_z[i];

        grr_z_xz_yy[i] = ts_x_yy[i] * gfe2_0 + gr_x_yy[i] * gfe_0 + ts_xz_yy[i] * gfe_0 * gc_z[i] + gr_xz_yy[i] * gc_z[i];

        grr_z_xz_yz[i] = ts_x_yz[i] * gfe2_0 + gr_x_yz[i] * gfe_0 + ts_xz_y[i] * gfe2_0 + gr_xz_y[i] * gfe_0 + ts_xz_yz[i] * gfe_0 * gc_z[i] + gr_xz_yz[i] * gc_z[i];

        grr_z_xz_zz[i] = ts_x_zz[i] * gfe2_0 + gr_x_zz[i] * gfe_0 + 2.0 * ts_xz_z[i] * gfe2_0 + 2.0 * gr_xz_z[i] * gfe_0 + ts_xz_zz[i] * gfe_0 * gc_z[i] + gr_xz_zz[i] * gc_z[i];
    }

    // Set up 90-96 components of targeted buffer : DD

    auto grr_z_yy_xx = pbuffer.data(idx_gr_dd + 90);

    auto grr_z_yy_xy = pbuffer.data(idx_gr_dd + 91);

    auto grr_z_yy_xz = pbuffer.data(idx_gr_dd + 92);

    auto grr_z_yy_yy = pbuffer.data(idx_gr_dd + 93);

    auto grr_z_yy_yz = pbuffer.data(idx_gr_dd + 94);

    auto grr_z_yy_zz = pbuffer.data(idx_gr_dd + 95);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_x, gr_yy_xx, gr_yy_xy, gr_yy_xz, gr_yy_y, gr_yy_yy, gr_yy_yz, gr_yy_z, gr_yy_zz, grr_z_yy_xx, grr_z_yy_xy, grr_z_yy_xz, grr_z_yy_yy, grr_z_yy_yz, grr_z_yy_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yy_xx[i] = ts_yy_xx[i] * gfe_0 * gc_z[i] + gr_yy_xx[i] * gc_z[i];

        grr_z_yy_xy[i] = ts_yy_xy[i] * gfe_0 * gc_z[i] + gr_yy_xy[i] * gc_z[i];

        grr_z_yy_xz[i] = ts_yy_x[i] * gfe2_0 + gr_yy_x[i] * gfe_0 + ts_yy_xz[i] * gfe_0 * gc_z[i] + gr_yy_xz[i] * gc_z[i];

        grr_z_yy_yy[i] = ts_yy_yy[i] * gfe_0 * gc_z[i] + gr_yy_yy[i] * gc_z[i];

        grr_z_yy_yz[i] = ts_yy_y[i] * gfe2_0 + gr_yy_y[i] * gfe_0 + ts_yy_yz[i] * gfe_0 * gc_z[i] + gr_yy_yz[i] * gc_z[i];

        grr_z_yy_zz[i] = 2.0 * ts_yy_z[i] * gfe2_0 + 2.0 * gr_yy_z[i] * gfe_0 + ts_yy_zz[i] * gfe_0 * gc_z[i] + gr_yy_zz[i] * gc_z[i];
    }

    // Set up 96-102 components of targeted buffer : DD

    auto grr_z_yz_xx = pbuffer.data(idx_gr_dd + 96);

    auto grr_z_yz_xy = pbuffer.data(idx_gr_dd + 97);

    auto grr_z_yz_xz = pbuffer.data(idx_gr_dd + 98);

    auto grr_z_yz_yy = pbuffer.data(idx_gr_dd + 99);

    auto grr_z_yz_yz = pbuffer.data(idx_gr_dd + 100);

    auto grr_z_yz_zz = pbuffer.data(idx_gr_dd + 101);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_yy, gr_y_yz, gr_y_zz, gr_yz_x, gr_yz_xx, gr_yz_xy, gr_yz_xz, gr_yz_y, gr_yz_yy, gr_yz_yz, gr_yz_z, gr_yz_zz, grr_z_yz_xx, grr_z_yz_xy, grr_z_yz_xz, grr_z_yz_yy, grr_z_yz_yz, grr_z_yz_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yz_xx[i] = ts_y_xx[i] * gfe2_0 + gr_y_xx[i] * gfe_0 + ts_yz_xx[i] * gfe_0 * gc_z[i] + gr_yz_xx[i] * gc_z[i];

        grr_z_yz_xy[i] = ts_y_xy[i] * gfe2_0 + gr_y_xy[i] * gfe_0 + ts_yz_xy[i] * gfe_0 * gc_z[i] + gr_yz_xy[i] * gc_z[i];

        grr_z_yz_xz[i] = ts_y_xz[i] * gfe2_0 + gr_y_xz[i] * gfe_0 + ts_yz_x[i] * gfe2_0 + gr_yz_x[i] * gfe_0 + ts_yz_xz[i] * gfe_0 * gc_z[i] + gr_yz_xz[i] * gc_z[i];

        grr_z_yz_yy[i] = ts_y_yy[i] * gfe2_0 + gr_y_yy[i] * gfe_0 + ts_yz_yy[i] * gfe_0 * gc_z[i] + gr_yz_yy[i] * gc_z[i];

        grr_z_yz_yz[i] = ts_y_yz[i] * gfe2_0 + gr_y_yz[i] * gfe_0 + ts_yz_y[i] * gfe2_0 + gr_yz_y[i] * gfe_0 + ts_yz_yz[i] * gfe_0 * gc_z[i] + gr_yz_yz[i] * gc_z[i];

        grr_z_yz_zz[i] = ts_y_zz[i] * gfe2_0 + gr_y_zz[i] * gfe_0 + 2.0 * ts_yz_z[i] * gfe2_0 + 2.0 * gr_yz_z[i] * gfe_0 + ts_yz_zz[i] * gfe_0 * gc_z[i] + gr_yz_zz[i] * gc_z[i];
    }

    // Set up 102-108 components of targeted buffer : DD

    auto grr_z_zz_xx = pbuffer.data(idx_gr_dd + 102);

    auto grr_z_zz_xy = pbuffer.data(idx_gr_dd + 103);

    auto grr_z_zz_xz = pbuffer.data(idx_gr_dd + 104);

    auto grr_z_zz_yy = pbuffer.data(idx_gr_dd + 105);

    auto grr_z_zz_yz = pbuffer.data(idx_gr_dd + 106);

    auto grr_z_zz_zz = pbuffer.data(idx_gr_dd + 107);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_yy, gr_z_yz, gr_z_zz, gr_zz_x, gr_zz_xx, gr_zz_xy, gr_zz_xz, gr_zz_y, gr_zz_yy, gr_zz_yz, gr_zz_z, gr_zz_zz, grr_z_zz_xx, grr_z_zz_xy, grr_z_zz_xz, grr_z_zz_yy, grr_z_zz_yz, grr_z_zz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zz_xx[i] = 2.0 * ts_z_xx[i] * gfe2_0 + 2.0 * gr_z_xx[i] * gfe_0 + ts_zz_xx[i] * gfe_0 * gc_z[i] + gr_zz_xx[i] * gc_z[i];

        grr_z_zz_xy[i] = 2.0 * ts_z_xy[i] * gfe2_0 + 2.0 * gr_z_xy[i] * gfe_0 + ts_zz_xy[i] * gfe_0 * gc_z[i] + gr_zz_xy[i] * gc_z[i];

        grr_z_zz_xz[i] = 2.0 * ts_z_xz[i] * gfe2_0 + 2.0 * gr_z_xz[i] * gfe_0 + ts_zz_x[i] * gfe2_0 + gr_zz_x[i] * gfe_0 + ts_zz_xz[i] * gfe_0 * gc_z[i] + gr_zz_xz[i] * gc_z[i];

        grr_z_zz_yy[i] = 2.0 * ts_z_yy[i] * gfe2_0 + 2.0 * gr_z_yy[i] * gfe_0 + ts_zz_yy[i] * gfe_0 * gc_z[i] + gr_zz_yy[i] * gc_z[i];

        grr_z_zz_yz[i] = 2.0 * ts_z_yz[i] * gfe2_0 + 2.0 * gr_z_yz[i] * gfe_0 + ts_zz_y[i] * gfe2_0 + gr_zz_y[i] * gfe_0 + ts_zz_yz[i] * gfe_0 * gc_z[i] + gr_zz_yz[i] * gc_z[i];

        grr_z_zz_zz[i] = 2.0 * ts_z_zz[i] * gfe2_0 + 2.0 * gr_z_zz[i] * gfe_0 + 2.0 * ts_zz_z[i] * gfe2_0 + 2.0 * gr_zz_z[i] * gfe_0 + ts_zz_zz[i] * gfe_0 * gc_z[i] + gr_zz_zz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

