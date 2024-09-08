#include "ElectronRepulsionPrimRecSPSF.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_spsf(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_spsf,
                                  size_t                idx_eri_1_sssd,
                                  size_t                idx_eri_0_sssf,
                                  size_t                idx_eri_1_sssf,
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

    /// Set up components of auxilary buffer : SSSD

    auto g_0_0_0_xx_1 = pbuffer.data(idx_eri_1_sssd);

    auto g_0_0_0_xy_1 = pbuffer.data(idx_eri_1_sssd + 1);

    auto g_0_0_0_xz_1 = pbuffer.data(idx_eri_1_sssd + 2);

    auto g_0_0_0_yy_1 = pbuffer.data(idx_eri_1_sssd + 3);

    auto g_0_0_0_yz_1 = pbuffer.data(idx_eri_1_sssd + 4);

    auto g_0_0_0_zz_1 = pbuffer.data(idx_eri_1_sssd + 5);

    /// Set up components of auxilary buffer : SSSF

    auto g_0_0_0_xxx_0 = pbuffer.data(idx_eri_0_sssf);

    auto g_0_0_0_xxy_0 = pbuffer.data(idx_eri_0_sssf + 1);

    auto g_0_0_0_xxz_0 = pbuffer.data(idx_eri_0_sssf + 2);

    auto g_0_0_0_xyy_0 = pbuffer.data(idx_eri_0_sssf + 3);

    auto g_0_0_0_xyz_0 = pbuffer.data(idx_eri_0_sssf + 4);

    auto g_0_0_0_xzz_0 = pbuffer.data(idx_eri_0_sssf + 5);

    auto g_0_0_0_yyy_0 = pbuffer.data(idx_eri_0_sssf + 6);

    auto g_0_0_0_yyz_0 = pbuffer.data(idx_eri_0_sssf + 7);

    auto g_0_0_0_yzz_0 = pbuffer.data(idx_eri_0_sssf + 8);

    auto g_0_0_0_zzz_0 = pbuffer.data(idx_eri_0_sssf + 9);

    /// Set up components of auxilary buffer : SSSF

    auto g_0_0_0_xxx_1 = pbuffer.data(idx_eri_1_sssf);

    auto g_0_0_0_xxy_1 = pbuffer.data(idx_eri_1_sssf + 1);

    auto g_0_0_0_xxz_1 = pbuffer.data(idx_eri_1_sssf + 2);

    auto g_0_0_0_xyy_1 = pbuffer.data(idx_eri_1_sssf + 3);

    auto g_0_0_0_xyz_1 = pbuffer.data(idx_eri_1_sssf + 4);

    auto g_0_0_0_xzz_1 = pbuffer.data(idx_eri_1_sssf + 5);

    auto g_0_0_0_yyy_1 = pbuffer.data(idx_eri_1_sssf + 6);

    auto g_0_0_0_yyz_1 = pbuffer.data(idx_eri_1_sssf + 7);

    auto g_0_0_0_yzz_1 = pbuffer.data(idx_eri_1_sssf + 8);

    auto g_0_0_0_zzz_1 = pbuffer.data(idx_eri_1_sssf + 9);

    /// Set up 0-10 components of targeted buffer : SPSF

    auto g_0_x_0_xxx_0 = pbuffer.data(idx_eri_0_spsf);

    auto g_0_x_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 1);

    auto g_0_x_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 2);

    auto g_0_x_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 3);

    auto g_0_x_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 4);

    auto g_0_x_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 5);

    auto g_0_x_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 6);

    auto g_0_x_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 7);

    auto g_0_x_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 8);

    auto g_0_x_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 9);

#pragma omp simd aligned(g_0_0_0_xx_1,      \
                             g_0_0_0_xxx_0, \
                             g_0_0_0_xxx_1, \
                             g_0_0_0_xxy_0, \
                             g_0_0_0_xxy_1, \
                             g_0_0_0_xxz_0, \
                             g_0_0_0_xxz_1, \
                             g_0_0_0_xy_1,  \
                             g_0_0_0_xyy_0, \
                             g_0_0_0_xyy_1, \
                             g_0_0_0_xyz_0, \
                             g_0_0_0_xyz_1, \
                             g_0_0_0_xz_1,  \
                             g_0_0_0_xzz_0, \
                             g_0_0_0_xzz_1, \
                             g_0_0_0_yy_1,  \
                             g_0_0_0_yyy_0, \
                             g_0_0_0_yyy_1, \
                             g_0_0_0_yyz_0, \
                             g_0_0_0_yyz_1, \
                             g_0_0_0_yz_1,  \
                             g_0_0_0_yzz_0, \
                             g_0_0_0_yzz_1, \
                             g_0_0_0_zz_1,  \
                             g_0_0_0_zzz_0, \
                             g_0_0_0_zzz_1, \
                             g_0_x_0_xxx_0, \
                             g_0_x_0_xxy_0, \
                             g_0_x_0_xxz_0, \
                             g_0_x_0_xyy_0, \
                             g_0_x_0_xyz_0, \
                             g_0_x_0_xzz_0, \
                             g_0_x_0_yyy_0, \
                             g_0_x_0_yyz_0, \
                             g_0_x_0_yzz_0, \
                             g_0_x_0_zzz_0, \
                             wp_x : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxx_0[i] = 3.0 * g_0_0_0_xx_1[i] * fi_abcd_0 + g_0_0_0_xxx_0[i] * pb_x + g_0_0_0_xxx_1[i] * wp_x[i];

        g_0_x_0_xxy_0[i] = 2.0 * g_0_0_0_xy_1[i] * fi_abcd_0 + g_0_0_0_xxy_0[i] * pb_x + g_0_0_0_xxy_1[i] * wp_x[i];

        g_0_x_0_xxz_0[i] = 2.0 * g_0_0_0_xz_1[i] * fi_abcd_0 + g_0_0_0_xxz_0[i] * pb_x + g_0_0_0_xxz_1[i] * wp_x[i];

        g_0_x_0_xyy_0[i] = g_0_0_0_yy_1[i] * fi_abcd_0 + g_0_0_0_xyy_0[i] * pb_x + g_0_0_0_xyy_1[i] * wp_x[i];

        g_0_x_0_xyz_0[i] = g_0_0_0_yz_1[i] * fi_abcd_0 + g_0_0_0_xyz_0[i] * pb_x + g_0_0_0_xyz_1[i] * wp_x[i];

        g_0_x_0_xzz_0[i] = g_0_0_0_zz_1[i] * fi_abcd_0 + g_0_0_0_xzz_0[i] * pb_x + g_0_0_0_xzz_1[i] * wp_x[i];

        g_0_x_0_yyy_0[i] = g_0_0_0_yyy_0[i] * pb_x + g_0_0_0_yyy_1[i] * wp_x[i];

        g_0_x_0_yyz_0[i] = g_0_0_0_yyz_0[i] * pb_x + g_0_0_0_yyz_1[i] * wp_x[i];

        g_0_x_0_yzz_0[i] = g_0_0_0_yzz_0[i] * pb_x + g_0_0_0_yzz_1[i] * wp_x[i];

        g_0_x_0_zzz_0[i] = g_0_0_0_zzz_0[i] * pb_x + g_0_0_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : SPSF

    auto g_0_y_0_xxx_0 = pbuffer.data(idx_eri_0_spsf + 10);

    auto g_0_y_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 11);

    auto g_0_y_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 12);

    auto g_0_y_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 13);

    auto g_0_y_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 14);

    auto g_0_y_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 15);

    auto g_0_y_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 16);

    auto g_0_y_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 17);

    auto g_0_y_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 18);

    auto g_0_y_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 19);

#pragma omp simd aligned(g_0_0_0_xx_1,      \
                             g_0_0_0_xxx_0, \
                             g_0_0_0_xxx_1, \
                             g_0_0_0_xxy_0, \
                             g_0_0_0_xxy_1, \
                             g_0_0_0_xxz_0, \
                             g_0_0_0_xxz_1, \
                             g_0_0_0_xy_1,  \
                             g_0_0_0_xyy_0, \
                             g_0_0_0_xyy_1, \
                             g_0_0_0_xyz_0, \
                             g_0_0_0_xyz_1, \
                             g_0_0_0_xz_1,  \
                             g_0_0_0_xzz_0, \
                             g_0_0_0_xzz_1, \
                             g_0_0_0_yy_1,  \
                             g_0_0_0_yyy_0, \
                             g_0_0_0_yyy_1, \
                             g_0_0_0_yyz_0, \
                             g_0_0_0_yyz_1, \
                             g_0_0_0_yz_1,  \
                             g_0_0_0_yzz_0, \
                             g_0_0_0_yzz_1, \
                             g_0_0_0_zz_1,  \
                             g_0_0_0_zzz_0, \
                             g_0_0_0_zzz_1, \
                             g_0_y_0_xxx_0, \
                             g_0_y_0_xxy_0, \
                             g_0_y_0_xxz_0, \
                             g_0_y_0_xyy_0, \
                             g_0_y_0_xyz_0, \
                             g_0_y_0_xzz_0, \
                             g_0_y_0_yyy_0, \
                             g_0_y_0_yyz_0, \
                             g_0_y_0_yzz_0, \
                             g_0_y_0_zzz_0, \
                             wp_y : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxx_0[i] = g_0_0_0_xxx_0[i] * pb_y + g_0_0_0_xxx_1[i] * wp_y[i];

        g_0_y_0_xxy_0[i] = g_0_0_0_xx_1[i] * fi_abcd_0 + g_0_0_0_xxy_0[i] * pb_y + g_0_0_0_xxy_1[i] * wp_y[i];

        g_0_y_0_xxz_0[i] = g_0_0_0_xxz_0[i] * pb_y + g_0_0_0_xxz_1[i] * wp_y[i];

        g_0_y_0_xyy_0[i] = 2.0 * g_0_0_0_xy_1[i] * fi_abcd_0 + g_0_0_0_xyy_0[i] * pb_y + g_0_0_0_xyy_1[i] * wp_y[i];

        g_0_y_0_xyz_0[i] = g_0_0_0_xz_1[i] * fi_abcd_0 + g_0_0_0_xyz_0[i] * pb_y + g_0_0_0_xyz_1[i] * wp_y[i];

        g_0_y_0_xzz_0[i] = g_0_0_0_xzz_0[i] * pb_y + g_0_0_0_xzz_1[i] * wp_y[i];

        g_0_y_0_yyy_0[i] = 3.0 * g_0_0_0_yy_1[i] * fi_abcd_0 + g_0_0_0_yyy_0[i] * pb_y + g_0_0_0_yyy_1[i] * wp_y[i];

        g_0_y_0_yyz_0[i] = 2.0 * g_0_0_0_yz_1[i] * fi_abcd_0 + g_0_0_0_yyz_0[i] * pb_y + g_0_0_0_yyz_1[i] * wp_y[i];

        g_0_y_0_yzz_0[i] = g_0_0_0_zz_1[i] * fi_abcd_0 + g_0_0_0_yzz_0[i] * pb_y + g_0_0_0_yzz_1[i] * wp_y[i];

        g_0_y_0_zzz_0[i] = g_0_0_0_zzz_0[i] * pb_y + g_0_0_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : SPSF

    auto g_0_z_0_xxx_0 = pbuffer.data(idx_eri_0_spsf + 20);

    auto g_0_z_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 21);

    auto g_0_z_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 22);

    auto g_0_z_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 23);

    auto g_0_z_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 24);

    auto g_0_z_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 25);

    auto g_0_z_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 26);

    auto g_0_z_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 27);

    auto g_0_z_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 28);

    auto g_0_z_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 29);

#pragma omp simd aligned(g_0_0_0_xx_1,      \
                             g_0_0_0_xxx_0, \
                             g_0_0_0_xxx_1, \
                             g_0_0_0_xxy_0, \
                             g_0_0_0_xxy_1, \
                             g_0_0_0_xxz_0, \
                             g_0_0_0_xxz_1, \
                             g_0_0_0_xy_1,  \
                             g_0_0_0_xyy_0, \
                             g_0_0_0_xyy_1, \
                             g_0_0_0_xyz_0, \
                             g_0_0_0_xyz_1, \
                             g_0_0_0_xz_1,  \
                             g_0_0_0_xzz_0, \
                             g_0_0_0_xzz_1, \
                             g_0_0_0_yy_1,  \
                             g_0_0_0_yyy_0, \
                             g_0_0_0_yyy_1, \
                             g_0_0_0_yyz_0, \
                             g_0_0_0_yyz_1, \
                             g_0_0_0_yz_1,  \
                             g_0_0_0_yzz_0, \
                             g_0_0_0_yzz_1, \
                             g_0_0_0_zz_1,  \
                             g_0_0_0_zzz_0, \
                             g_0_0_0_zzz_1, \
                             g_0_z_0_xxx_0, \
                             g_0_z_0_xxy_0, \
                             g_0_z_0_xxz_0, \
                             g_0_z_0_xyy_0, \
                             g_0_z_0_xyz_0, \
                             g_0_z_0_xzz_0, \
                             g_0_z_0_yyy_0, \
                             g_0_z_0_yyz_0, \
                             g_0_z_0_yzz_0, \
                             g_0_z_0_zzz_0, \
                             wp_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxx_0[i] = g_0_0_0_xxx_0[i] * pb_z + g_0_0_0_xxx_1[i] * wp_z[i];

        g_0_z_0_xxy_0[i] = g_0_0_0_xxy_0[i] * pb_z + g_0_0_0_xxy_1[i] * wp_z[i];

        g_0_z_0_xxz_0[i] = g_0_0_0_xx_1[i] * fi_abcd_0 + g_0_0_0_xxz_0[i] * pb_z + g_0_0_0_xxz_1[i] * wp_z[i];

        g_0_z_0_xyy_0[i] = g_0_0_0_xyy_0[i] * pb_z + g_0_0_0_xyy_1[i] * wp_z[i];

        g_0_z_0_xyz_0[i] = g_0_0_0_xy_1[i] * fi_abcd_0 + g_0_0_0_xyz_0[i] * pb_z + g_0_0_0_xyz_1[i] * wp_z[i];

        g_0_z_0_xzz_0[i] = 2.0 * g_0_0_0_xz_1[i] * fi_abcd_0 + g_0_0_0_xzz_0[i] * pb_z + g_0_0_0_xzz_1[i] * wp_z[i];

        g_0_z_0_yyy_0[i] = g_0_0_0_yyy_0[i] * pb_z + g_0_0_0_yyy_1[i] * wp_z[i];

        g_0_z_0_yyz_0[i] = g_0_0_0_yy_1[i] * fi_abcd_0 + g_0_0_0_yyz_0[i] * pb_z + g_0_0_0_yyz_1[i] * wp_z[i];

        g_0_z_0_yzz_0[i] = 2.0 * g_0_0_0_yz_1[i] * fi_abcd_0 + g_0_0_0_yzz_0[i] * pb_z + g_0_0_0_yzz_1[i] * wp_z[i];

        g_0_z_0_zzz_0[i] = 3.0 * g_0_0_0_zz_1[i] * fi_abcd_0 + g_0_0_0_zzz_0[i] * pb_z + g_0_0_0_zzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
