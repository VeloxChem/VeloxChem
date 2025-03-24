#include "ThreeCenterOverlapPrimRecGS.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_gs(CSimdArray<double>& pbuffer, 
                     const size_t idx_gs,
                     const size_t idx_ds,
                     const size_t idx_fs,
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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_xxz_0 = pbuffer.data(idx_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_fs + 3);

    auto ts_xzz_0 = pbuffer.data(idx_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of targeted buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_gs);

    auto ts_xxxy_0 = pbuffer.data(idx_gs + 1);

    auto ts_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto ts_xxyz_0 = pbuffer.data(idx_gs + 4);

    auto ts_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto ts_xyyz_0 = pbuffer.data(idx_gs + 7);

    auto ts_xyzz_0 = pbuffer.data(idx_gs + 8);

    auto ts_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_gs + 14);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xx_0, ts_xxx_0, ts_xxxx_0, ts_xxxy_0, ts_xxxz_0, ts_xxyy_0, ts_xxyz_0, ts_xxz_0, ts_xxzz_0, ts_xyy_0, ts_xyyy_0, ts_xyyz_0, ts_xyzz_0, ts_xzz_0, ts_xzzz_0, ts_yy_0, ts_yyy_0, ts_yyyy_0, ts_yyyz_0, ts_yyz_0, ts_yyzz_0, ts_yzz_0, ts_yzzz_0, ts_zz_0, ts_zzz_0, ts_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxx_0[i] = 3.0 * ts_xx_0[i] * gfe_0 + ts_xxx_0[i] * ga_x[i];

        ts_xxxy_0[i] = ts_xxx_0[i] * ga_y[i];

        ts_xxxz_0[i] = ts_xxx_0[i] * ga_z[i];

        ts_xxyy_0[i] = ts_yy_0[i] * gfe_0 + ts_xyy_0[i] * ga_x[i];

        ts_xxyz_0[i] = ts_xxz_0[i] * ga_y[i];

        ts_xxzz_0[i] = ts_zz_0[i] * gfe_0 + ts_xzz_0[i] * ga_x[i];

        ts_xyyy_0[i] = ts_yyy_0[i] * ga_x[i];

        ts_xyyz_0[i] = ts_yyz_0[i] * ga_x[i];

        ts_xyzz_0[i] = ts_yzz_0[i] * ga_x[i];

        ts_xzzz_0[i] = ts_zzz_0[i] * ga_x[i];

        ts_yyyy_0[i] = 3.0 * ts_yy_0[i] * gfe_0 + ts_yyy_0[i] * ga_y[i];

        ts_yyyz_0[i] = ts_yyy_0[i] * ga_z[i];

        ts_yyzz_0[i] = ts_zz_0[i] * gfe_0 + ts_yzz_0[i] * ga_y[i];

        ts_yzzz_0[i] = ts_zzz_0[i] * ga_y[i];

        ts_zzzz_0[i] = 3.0 * ts_zz_0[i] * gfe_0 + ts_zzz_0[i] * ga_z[i];
    }
}

} // t3ovlrec namespace

