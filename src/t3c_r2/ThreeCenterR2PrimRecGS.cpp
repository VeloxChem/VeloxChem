#include "ThreeCenterR2PrimRecGS.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_gs(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gs,
                const size_t idx_ds,
                const size_t idx_fs,
                const size_t idx_gs,
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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_xxy_0 = pbuffer.data(idx_fs + 1);

    auto ts_xxz_0 = pbuffer.data(idx_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_fs + 3);

    auto ts_xyz_0 = pbuffer.data(idx_fs + 4);

    auto ts_xzz_0 = pbuffer.data(idx_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : GS

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

    // Set up components of targeted buffer : GS

    auto gr_xxxx_0 = pbuffer.data(idx_g_gs);

    auto gr_xxxy_0 = pbuffer.data(idx_g_gs + 1);

    auto gr_xxxz_0 = pbuffer.data(idx_g_gs + 2);

    auto gr_xxyy_0 = pbuffer.data(idx_g_gs + 3);

    auto gr_xxyz_0 = pbuffer.data(idx_g_gs + 4);

    auto gr_xxzz_0 = pbuffer.data(idx_g_gs + 5);

    auto gr_xyyy_0 = pbuffer.data(idx_g_gs + 6);

    auto gr_xyyz_0 = pbuffer.data(idx_g_gs + 7);

    auto gr_xyzz_0 = pbuffer.data(idx_g_gs + 8);

    auto gr_xzzz_0 = pbuffer.data(idx_g_gs + 9);

    auto gr_yyyy_0 = pbuffer.data(idx_g_gs + 10);

    auto gr_yyyz_0 = pbuffer.data(idx_g_gs + 11);

    auto gr_yyzz_0 = pbuffer.data(idx_g_gs + 12);

    auto gr_yzzz_0 = pbuffer.data(idx_g_gs + 13);

    auto gr_zzzz_0 = pbuffer.data(idx_g_gs + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_0, gr_xxxy_0, gr_xxxz_0, gr_xxyy_0, gr_xxyz_0, gr_xxzz_0, gr_xyyy_0, gr_xyyz_0, gr_xyzz_0, gr_xzzz_0, gr_yyyy_0, gr_yyyz_0, gr_yyzz_0, gr_yzzz_0, gr_zzzz_0, ts_xx_0, ts_xxx_0, ts_xxxx_0, ts_xxxy_0, ts_xxxz_0, ts_xxy_0, ts_xxyy_0, ts_xxyz_0, ts_xxz_0, ts_xxzz_0, ts_xy_0, ts_xyy_0, ts_xyyy_0, ts_xyyz_0, ts_xyz_0, ts_xyzz_0, ts_xz_0, ts_xzz_0, ts_xzzz_0, ts_yy_0, ts_yyy_0, ts_yyyy_0, ts_yyyz_0, ts_yyz_0, ts_yyzz_0, ts_yz_0, ts_yzz_0, ts_yzzz_0, ts_zz_0, ts_zzz_0, ts_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxxx_0[i] = 12.0 * ts_xx_0[i] * gfe2_0 + 8.0 * ts_xxx_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxx_0[i] * gfe_0 + ts_xxxx_0[i] * rgc2_0;

        gr_xxxy_0[i] = 6.0 * ts_xy_0[i] * gfe2_0 + 6.0 * ts_xxy_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_0[i] * gfe_0 + ts_xxxy_0[i] * rgc2_0;

        gr_xxxz_0[i] = 6.0 * ts_xz_0[i] * gfe2_0 + 6.0 * ts_xxz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_0[i] * gfe_0 + ts_xxxz_0[i] * rgc2_0;

        gr_xxyy_0[i] = 2.0 * ts_yy_0[i] * gfe2_0 + 4.0 * ts_xyy_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe2_0 + 4.0 * ts_xxy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_0[i] * gfe_0 + ts_xxyy_0[i] * rgc2_0;

        gr_xxyz_0[i] = 2.0 * ts_yz_0[i] * gfe2_0 + 4.0 * ts_xyz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_0[i] * gfe_0 + ts_xxyz_0[i] * rgc2_0;

        gr_xxzz_0[i] = 2.0 * ts_zz_0[i] * gfe2_0 + 4.0 * ts_xzz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe2_0 + 4.0 * ts_xxz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_0[i] * gfe_0 + ts_xxzz_0[i] * rgc2_0;

        gr_xyyy_0[i] = 2.0 * ts_yyy_0[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_0[i] * gfe2_0 + 6.0 * ts_xyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_0[i] * gfe_0 + ts_xyyy_0[i] * rgc2_0;

        gr_xyyz_0[i] = 2.0 * ts_yyz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_0[i] * gfe2_0 + 4.0 * ts_xyz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_0[i] * gfe_0 + ts_xyyz_0[i] * rgc2_0;

        gr_xyzz_0[i] = 2.0 * ts_yzz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe2_0 + 4.0 * ts_xyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_0[i] * gfe_0 + ts_xyzz_0[i] * rgc2_0;

        gr_xzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_0[i] * gfe2_0 + 6.0 * ts_xzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_0[i] * gfe_0 + ts_xzzz_0[i] * rgc2_0;

        gr_yyyy_0[i] = 12.0 * ts_yy_0[i] * gfe2_0 + 8.0 * ts_yyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_0[i] * gfe_0 + ts_yyyy_0[i] * rgc2_0;

        gr_yyyz_0[i] = 6.0 * ts_yz_0[i] * gfe2_0 + 6.0 * ts_yyz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_0[i] * gfe_0 + ts_yyyz_0[i] * rgc2_0;

        gr_yyzz_0[i] = 2.0 * ts_zz_0[i] * gfe2_0 + 4.0 * ts_yzz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe2_0 + 4.0 * ts_yyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_0[i] * gfe_0 + ts_yyzz_0[i] * rgc2_0;

        gr_yzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_0[i] * gfe2_0 + 6.0 * ts_yzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_0[i] * gfe_0 + ts_yzzz_0[i] * rgc2_0;

        gr_zzzz_0[i] = 12.0 * ts_zz_0[i] * gfe2_0 + 8.0 * ts_zzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_0[i] * gfe_0 + ts_zzzz_0[i] * rgc2_0;
    }
}

} // t3r2rec namespace

