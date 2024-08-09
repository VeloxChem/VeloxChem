#include "OverlapPrimRecSD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_sd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_sd,
                     const size_t              idx_ovl_ss,
                     const size_t              idx_ovl_sp,
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

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of targeted buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

#pragma omp simd aligned(pb_x,        \
                             pb_y,    \
                             pb_z,    \
                             ts_0_0,  \
                             ts_0_x,  \
                             ts_0_xx, \
                             ts_0_xy, \
                             ts_0_xz, \
                             ts_0_y,  \
                             ts_0_yy, \
                             ts_0_yz, \
                             ts_0_z,  \
                             ts_0_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xx[i] = ts_0_0[i] * fe_0 + ts_0_x[i] * pb_x[i];

        ts_0_xy[i] = ts_0_y[i] * pb_x[i];

        ts_0_xz[i] = ts_0_z[i] * pb_x[i];

        ts_0_yy[i] = ts_0_0[i] * fe_0 + ts_0_y[i] * pb_y[i];

        ts_0_yz[i] = ts_0_z[i] * pb_y[i];

        ts_0_zz[i] = ts_0_0[i] * fe_0 + ts_0_z[i] * pb_z[i];
    }
}

}  // namespace ovlrec
