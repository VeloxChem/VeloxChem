#include "ThreeCenterOverlapPrimRecPD.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_pd(CSimdArray<double>& pbuffer, 
                     const size_t idx_pd,
                     const size_t idx_sp,
                     const size_t idx_sd,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up 0-6 components of targeted buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    #pragma omp simd aligned(ga_x, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_x_xx[i] = 2.0 * ts_0_x[i] * gfe_0 + ts_0_xx[i] * ga_x[i];

        ts_x_xy[i] = ts_0_y[i] * gfe_0 + ts_0_xy[i] * ga_x[i];

        ts_x_xz[i] = ts_0_z[i] * gfe_0 + ts_0_xz[i] * ga_x[i];

        ts_x_yy[i] = ts_0_yy[i] * ga_x[i];

        ts_x_yz[i] = ts_0_yz[i] * ga_x[i];

        ts_x_zz[i] = ts_0_zz[i] * ga_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    #pragma omp simd aligned(ga_y, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_y_xx[i] = ts_0_xx[i] * ga_y[i];

        ts_y_xy[i] = ts_0_x[i] * gfe_0 + ts_0_xy[i] * ga_y[i];

        ts_y_xz[i] = ts_0_xz[i] * ga_y[i];

        ts_y_yy[i] = 2.0 * ts_0_y[i] * gfe_0 + ts_0_yy[i] * ga_y[i];

        ts_y_yz[i] = ts_0_z[i] * gfe_0 + ts_0_yz[i] * ga_y[i];

        ts_y_zz[i] = ts_0_zz[i] * ga_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

    #pragma omp simd aligned(ga_z, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_z_xx[i] = ts_0_xx[i] * ga_z[i];

        ts_z_xy[i] = ts_0_xy[i] * ga_z[i];

        ts_z_xz[i] = ts_0_x[i] * gfe_0 + ts_0_xz[i] * ga_z[i];

        ts_z_yy[i] = ts_0_yy[i] * ga_z[i];

        ts_z_yz[i] = ts_0_y[i] * gfe_0 + ts_0_yz[i] * ga_z[i];

        ts_z_zz[i] = 2.0 * ts_0_z[i] * gfe_0 + ts_0_zz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

