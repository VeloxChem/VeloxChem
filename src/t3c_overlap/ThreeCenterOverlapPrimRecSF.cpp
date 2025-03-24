#include "ThreeCenterOverlapPrimRecSF.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_sf(CSimdArray<double>& pbuffer, 
                     const size_t idx_sf,
                     const size_t idx_sp,
                     const size_t idx_sd,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GB) distances

    auto gb_x = factors.data(idx_rgb);

    auto gb_y = factors.data(idx_rgb + 1);

    auto gb_z = factors.data(idx_rgb + 2);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of targeted buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xxy = pbuffer.data(idx_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    #pragma omp simd aligned(gb_x, gb_y, gb_z, ts_0_x, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_y, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_z, ts_0_zz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_0_xxx[i] = 2.0 * ts_0_x[i] * gfe_0 + ts_0_xx[i] * gb_x[i];

        ts_0_xxy[i] = ts_0_xx[i] * gb_y[i];

        ts_0_xxz[i] = ts_0_xx[i] * gb_z[i];

        ts_0_xyy[i] = ts_0_yy[i] * gb_x[i];

        ts_0_xyz[i] = ts_0_yz[i] * gb_x[i];

        ts_0_xzz[i] = ts_0_zz[i] * gb_x[i];

        ts_0_yyy[i] = 2.0 * ts_0_y[i] * gfe_0 + ts_0_yy[i] * gb_y[i];

        ts_0_yyz[i] = ts_0_yy[i] * gb_z[i];

        ts_0_yzz[i] = ts_0_zz[i] * gb_y[i];

        ts_0_zzz[i] = 2.0 * ts_0_z[i] * gfe_0 + ts_0_zz[i] * gb_z[i];
    }
}

} // t3ovlrec namespace

