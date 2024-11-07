#include "OverlapPrimRecSF.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_sf(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_sf,
                     const size_t              idx_ovl_sp,
                     const size_t              idx_ovl_sd,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of targeted buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xxy = pbuffer.data(idx_ovl_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_ovl_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_ovl_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_ovl_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

#pragma omp simd aligned(pb_x,         \
                             pb_y,     \
                             pb_z,     \
                             ts_0_x,   \
                             ts_0_xx,  \
                             ts_0_xxx, \
                             ts_0_xxy, \
                             ts_0_xxz, \
                             ts_0_xyy, \
                             ts_0_xyz, \
                             ts_0_xzz, \
                             ts_0_y,   \
                             ts_0_yy,  \
                             ts_0_yyy, \
                             ts_0_yyz, \
                             ts_0_yz,  \
                             ts_0_yzz, \
                             ts_0_z,   \
                             ts_0_zz,  \
                             ts_0_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xxx[i] = 2.0 * ts_0_x[i] * fe_0 + ts_0_xx[i] * pb_x[i];

        ts_0_xxy[i] = ts_0_xx[i] * pb_y[i];

        ts_0_xxz[i] = ts_0_xx[i] * pb_z[i];

        ts_0_xyy[i] = ts_0_yy[i] * pb_x[i];

        ts_0_xyz[i] = ts_0_yz[i] * pb_x[i];

        ts_0_xzz[i] = ts_0_zz[i] * pb_x[i];

        ts_0_yyy[i] = 2.0 * ts_0_y[i] * fe_0 + ts_0_yy[i] * pb_y[i];

        ts_0_yyz[i] = ts_0_yy[i] * pb_z[i];

        ts_0_yzz[i] = ts_0_zz[i] * pb_y[i];

        ts_0_zzz[i] = 2.0 * ts_0_z[i] * fe_0 + ts_0_zz[i] * pb_z[i];
    }
}

}  // namespace ovlrec
