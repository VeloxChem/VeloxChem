#include "ThreeCenterOverlapPrimRecSG.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_sg(CSimdArray<double>& pbuffer, 
                     const size_t idx_sg,
                     const size_t idx_sd,
                     const size_t idx_sf,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xxz = pbuffer.data(idx_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of targeted buffer : SG

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

    #pragma omp simd aligned(gb_x, gb_y, gb_z, ts_0_xx, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yy, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zz, ts_0_zzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_0_xxxx[i] = 3.0 * ts_0_xx[i] * gfe_0 + ts_0_xxx[i] * gb_x[i];

        ts_0_xxxy[i] = ts_0_xxx[i] * gb_y[i];

        ts_0_xxxz[i] = ts_0_xxx[i] * gb_z[i];

        ts_0_xxyy[i] = ts_0_yy[i] * gfe_0 + ts_0_xyy[i] * gb_x[i];

        ts_0_xxyz[i] = ts_0_xxz[i] * gb_y[i];

        ts_0_xxzz[i] = ts_0_zz[i] * gfe_0 + ts_0_xzz[i] * gb_x[i];

        ts_0_xyyy[i] = ts_0_yyy[i] * gb_x[i];

        ts_0_xyyz[i] = ts_0_yyz[i] * gb_x[i];

        ts_0_xyzz[i] = ts_0_yzz[i] * gb_x[i];

        ts_0_xzzz[i] = ts_0_zzz[i] * gb_x[i];

        ts_0_yyyy[i] = 3.0 * ts_0_yy[i] * gfe_0 + ts_0_yyy[i] * gb_y[i];

        ts_0_yyyz[i] = ts_0_yyy[i] * gb_z[i];

        ts_0_yyzz[i] = ts_0_zz[i] * gfe_0 + ts_0_yzz[i] * gb_y[i];

        ts_0_yzzz[i] = ts_0_zzz[i] * gb_y[i];

        ts_0_zzzz[i] = 3.0 * ts_0_zz[i] * gfe_0 + ts_0_zzz[i] * gb_z[i];
    }
}

} // t3ovlrec namespace

