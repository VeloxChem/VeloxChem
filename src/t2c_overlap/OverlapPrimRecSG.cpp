#include "OverlapPrimRecSG.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_sg(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_sg,
                     const size_t              idx_ovl_sd,
                     const size_t              idx_ovl_sf,
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

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xxz = pbuffer.data(idx_ovl_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_ovl_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of targeted buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxxy = pbuffer.data(idx_ovl_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_ovl_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_ovl_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_ovl_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_ovl_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_ovl_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

#pragma omp simd aligned(pb_x,          \
                             pb_y,      \
                             pb_z,      \
                             ts_0_xx,   \
                             ts_0_xxx,  \
                             ts_0_xxxx, \
                             ts_0_xxxy, \
                             ts_0_xxxz, \
                             ts_0_xxyy, \
                             ts_0_xxyz, \
                             ts_0_xxz,  \
                             ts_0_xxzz, \
                             ts_0_xyy,  \
                             ts_0_xyyy, \
                             ts_0_xyyz, \
                             ts_0_xyzz, \
                             ts_0_xzz,  \
                             ts_0_xzzz, \
                             ts_0_yy,   \
                             ts_0_yyy,  \
                             ts_0_yyyy, \
                             ts_0_yyyz, \
                             ts_0_yyz,  \
                             ts_0_yyzz, \
                             ts_0_yzz,  \
                             ts_0_yzzz, \
                             ts_0_zz,   \
                             ts_0_zzz,  \
                             ts_0_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xxxx[i] = 3.0 * ts_0_xx[i] * fe_0 + ts_0_xxx[i] * pb_x[i];

        ts_0_xxxy[i] = ts_0_xxx[i] * pb_y[i];

        ts_0_xxxz[i] = ts_0_xxx[i] * pb_z[i];

        ts_0_xxyy[i] = ts_0_yy[i] * fe_0 + ts_0_xyy[i] * pb_x[i];

        ts_0_xxyz[i] = ts_0_xxz[i] * pb_y[i];

        ts_0_xxzz[i] = ts_0_zz[i] * fe_0 + ts_0_xzz[i] * pb_x[i];

        ts_0_xyyy[i] = ts_0_yyy[i] * pb_x[i];

        ts_0_xyyz[i] = ts_0_yyz[i] * pb_x[i];

        ts_0_xyzz[i] = ts_0_yzz[i] * pb_x[i];

        ts_0_xzzz[i] = ts_0_zzz[i] * pb_x[i];

        ts_0_yyyy[i] = 3.0 * ts_0_yy[i] * fe_0 + ts_0_yyy[i] * pb_y[i];

        ts_0_yyyz[i] = ts_0_yyy[i] * pb_z[i];

        ts_0_yyzz[i] = ts_0_zz[i] * fe_0 + ts_0_yzz[i] * pb_y[i];

        ts_0_yzzz[i] = ts_0_zzz[i] * pb_y[i];

        ts_0_zzzz[i] = 3.0 * ts_0_zz[i] * fe_0 + ts_0_zzz[i] * pb_z[i];
    }
}

}  // namespace ovlrec
