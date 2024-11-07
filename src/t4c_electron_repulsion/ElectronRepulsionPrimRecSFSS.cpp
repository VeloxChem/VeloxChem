#include "ElectronRepulsionPrimRecSFSS.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfss(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfss,
                                  size_t                idx_eri_0_spss,
                                  size_t                idx_eri_1_spss,
                                  size_t                idx_eri_0_sdss,
                                  size_t                idx_eri_1_sdss,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SPSS

    auto g_0_x_0_0_0 = pbuffer.data(idx_eri_0_spss);

    auto g_0_y_0_0_0 = pbuffer.data(idx_eri_0_spss + 1);

    auto g_0_z_0_0_0 = pbuffer.data(idx_eri_0_spss + 2);

    /// Set up components of auxilary buffer : SPSS

    auto g_0_x_0_0_1 = pbuffer.data(idx_eri_1_spss);

    auto g_0_y_0_0_1 = pbuffer.data(idx_eri_1_spss + 1);

    auto g_0_z_0_0_1 = pbuffer.data(idx_eri_1_spss + 2);

    /// Set up components of auxilary buffer : SDSS

    auto g_0_xx_0_0_0 = pbuffer.data(idx_eri_0_sdss);

    auto g_0_yy_0_0_0 = pbuffer.data(idx_eri_0_sdss + 3);

    auto g_0_yz_0_0_0 = pbuffer.data(idx_eri_0_sdss + 4);

    auto g_0_zz_0_0_0 = pbuffer.data(idx_eri_0_sdss + 5);

    /// Set up components of auxilary buffer : SDSS

    auto g_0_xx_0_0_1 = pbuffer.data(idx_eri_1_sdss);

    auto g_0_yy_0_0_1 = pbuffer.data(idx_eri_1_sdss + 3);

    auto g_0_yz_0_0_1 = pbuffer.data(idx_eri_1_sdss + 4);

    auto g_0_zz_0_0_1 = pbuffer.data(idx_eri_1_sdss + 5);

    /// Set up components of targeted buffer : SFSS

    auto g_0_xxx_0_0_0 = pbuffer.data(idx_eri_0_sfss);

    auto g_0_xxy_0_0_0 = pbuffer.data(idx_eri_0_sfss + 1);

    auto g_0_xxz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 2);

    auto g_0_xyy_0_0_0 = pbuffer.data(idx_eri_0_sfss + 3);

    auto g_0_xyz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 4);

    auto g_0_xzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 5);

    auto g_0_yyy_0_0_0 = pbuffer.data(idx_eri_0_sfss + 6);

    auto g_0_yyz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 7);

    auto g_0_yzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 8);

    auto g_0_zzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 9);

#pragma omp simd aligned(g_0_x_0_0_0,       \
                             g_0_x_0_0_1,   \
                             g_0_xx_0_0_0,  \
                             g_0_xx_0_0_1,  \
                             g_0_xxx_0_0_0, \
                             g_0_xxy_0_0_0, \
                             g_0_xxz_0_0_0, \
                             g_0_xyy_0_0_0, \
                             g_0_xyz_0_0_0, \
                             g_0_xzz_0_0_0, \
                             g_0_y_0_0_0,   \
                             g_0_y_0_0_1,   \
                             g_0_yy_0_0_0,  \
                             g_0_yy_0_0_1,  \
                             g_0_yyy_0_0_0, \
                             g_0_yyz_0_0_0, \
                             g_0_yz_0_0_0,  \
                             g_0_yz_0_0_1,  \
                             g_0_yzz_0_0_0, \
                             g_0_z_0_0_0,   \
                             g_0_z_0_0_1,   \
                             g_0_zz_0_0_0,  \
                             g_0_zz_0_0_1,  \
                             g_0_zzz_0_0_0, \
                             wp_x,          \
                             wp_y,          \
                             wp_z,          \
                             c_exps,        \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxx_0_0_0[i] = 2.0 * g_0_x_0_0_0[i] * fi_ab_0 - 2.0 * g_0_x_0_0_1[i] * fti_ab_0 + g_0_xx_0_0_0[i] * pb_x + g_0_xx_0_0_1[i] * wp_x[i];

        g_0_xxy_0_0_0[i] = g_0_xx_0_0_0[i] * pb_y + g_0_xx_0_0_1[i] * wp_y[i];

        g_0_xxz_0_0_0[i] = g_0_xx_0_0_0[i] * pb_z + g_0_xx_0_0_1[i] * wp_z[i];

        g_0_xyy_0_0_0[i] = g_0_yy_0_0_0[i] * pb_x + g_0_yy_0_0_1[i] * wp_x[i];

        g_0_xyz_0_0_0[i] = g_0_yz_0_0_0[i] * pb_x + g_0_yz_0_0_1[i] * wp_x[i];

        g_0_xzz_0_0_0[i] = g_0_zz_0_0_0[i] * pb_x + g_0_zz_0_0_1[i] * wp_x[i];

        g_0_yyy_0_0_0[i] = 2.0 * g_0_y_0_0_0[i] * fi_ab_0 - 2.0 * g_0_y_0_0_1[i] * fti_ab_0 + g_0_yy_0_0_0[i] * pb_y + g_0_yy_0_0_1[i] * wp_y[i];

        g_0_yyz_0_0_0[i] = g_0_yy_0_0_0[i] * pb_z + g_0_yy_0_0_1[i] * wp_z[i];

        g_0_yzz_0_0_0[i] = g_0_zz_0_0_0[i] * pb_y + g_0_zz_0_0_1[i] * wp_y[i];

        g_0_zzz_0_0_0[i] = 2.0 * g_0_z_0_0_0[i] * fi_ab_0 - 2.0 * g_0_z_0_0_1[i] * fti_ab_0 + g_0_zz_0_0_0[i] * pb_z + g_0_zz_0_0_1[i] * wp_z[i];
    }
}

}  // namespace erirec
