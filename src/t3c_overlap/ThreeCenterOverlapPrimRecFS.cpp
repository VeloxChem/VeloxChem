#include "ThreeCenterOverlapPrimRecFS.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_fs(CSimdArray<double>& pbuffer, 
                     const size_t idx_fs,
                     const size_t idx_ps,
                     const size_t idx_ds,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of targeted buffer : FS

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

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_x_0, ts_xx_0, ts_xxx_0, ts_xxy_0, ts_xxz_0, ts_xyy_0, ts_xyz_0, ts_xzz_0, ts_y_0, ts_yy_0, ts_yyy_0, ts_yyz_0, ts_yz_0, ts_yzz_0, ts_z_0, ts_zz_0, ts_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxx_0[i] = 2.0 * ts_x_0[i] * gfe_0 + ts_xx_0[i] * ga_x[i];

        ts_xxy_0[i] = ts_xx_0[i] * ga_y[i];

        ts_xxz_0[i] = ts_xx_0[i] * ga_z[i];

        ts_xyy_0[i] = ts_yy_0[i] * ga_x[i];

        ts_xyz_0[i] = ts_yz_0[i] * ga_x[i];

        ts_xzz_0[i] = ts_zz_0[i] * ga_x[i];

        ts_yyy_0[i] = 2.0 * ts_y_0[i] * gfe_0 + ts_yy_0[i] * ga_y[i];

        ts_yyz_0[i] = ts_yy_0[i] * ga_z[i];

        ts_yzz_0[i] = ts_zz_0[i] * ga_y[i];

        ts_zzz_0[i] = 2.0 * ts_z_0[i] * gfe_0 + ts_zz_0[i] * ga_z[i];
    }
}

} // t3ovlrec namespace

