#include "TwoCenterElectronRepulsionPrimRecGF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gf,
                                const size_t idx_eri_0_df,
                                const size_t idx_eri_1_df,
                                const size_t idx_eri_1_fd,
                                const size_t idx_eri_1_ff,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DF

    auto g_xx_xxx_0 = pbuffer.data(idx_eri_0_df);

    auto g_xx_xxy_0 = pbuffer.data(idx_eri_0_df + 1);

    auto g_xx_xxz_0 = pbuffer.data(idx_eri_0_df + 2);

    auto g_xx_xyy_0 = pbuffer.data(idx_eri_0_df + 3);

    auto g_xx_xyz_0 = pbuffer.data(idx_eri_0_df + 4);

    auto g_xx_xzz_0 = pbuffer.data(idx_eri_0_df + 5);

    auto g_xx_yyy_0 = pbuffer.data(idx_eri_0_df + 6);

    auto g_xx_yyz_0 = pbuffer.data(idx_eri_0_df + 7);

    auto g_xx_yzz_0 = pbuffer.data(idx_eri_0_df + 8);

    auto g_xx_zzz_0 = pbuffer.data(idx_eri_0_df + 9);

    auto g_yy_xxx_0 = pbuffer.data(idx_eri_0_df + 30);

    auto g_yy_xxy_0 = pbuffer.data(idx_eri_0_df + 31);

    auto g_yy_xxz_0 = pbuffer.data(idx_eri_0_df + 32);

    auto g_yy_xyy_0 = pbuffer.data(idx_eri_0_df + 33);

    auto g_yy_xyz_0 = pbuffer.data(idx_eri_0_df + 34);

    auto g_yy_xzz_0 = pbuffer.data(idx_eri_0_df + 35);

    auto g_yy_yyy_0 = pbuffer.data(idx_eri_0_df + 36);

    auto g_yy_yyz_0 = pbuffer.data(idx_eri_0_df + 37);

    auto g_yy_yzz_0 = pbuffer.data(idx_eri_0_df + 38);

    auto g_yy_zzz_0 = pbuffer.data(idx_eri_0_df + 39);

    auto g_zz_xxx_0 = pbuffer.data(idx_eri_0_df + 50);

    auto g_zz_xxy_0 = pbuffer.data(idx_eri_0_df + 51);

    auto g_zz_xxz_0 = pbuffer.data(idx_eri_0_df + 52);

    auto g_zz_xyy_0 = pbuffer.data(idx_eri_0_df + 53);

    auto g_zz_xyz_0 = pbuffer.data(idx_eri_0_df + 54);

    auto g_zz_xzz_0 = pbuffer.data(idx_eri_0_df + 55);

    auto g_zz_yyy_0 = pbuffer.data(idx_eri_0_df + 56);

    auto g_zz_yyz_0 = pbuffer.data(idx_eri_0_df + 57);

    auto g_zz_yzz_0 = pbuffer.data(idx_eri_0_df + 58);

    auto g_zz_zzz_0 = pbuffer.data(idx_eri_0_df + 59);

    // Set up components of auxiliary buffer : DF

    auto g_xx_xxx_1 = pbuffer.data(idx_eri_1_df);

    auto g_xx_xxy_1 = pbuffer.data(idx_eri_1_df + 1);

    auto g_xx_xxz_1 = pbuffer.data(idx_eri_1_df + 2);

    auto g_xx_xyy_1 = pbuffer.data(idx_eri_1_df + 3);

    auto g_xx_xyz_1 = pbuffer.data(idx_eri_1_df + 4);

    auto g_xx_xzz_1 = pbuffer.data(idx_eri_1_df + 5);

    auto g_xx_yyy_1 = pbuffer.data(idx_eri_1_df + 6);

    auto g_xx_yyz_1 = pbuffer.data(idx_eri_1_df + 7);

    auto g_xx_yzz_1 = pbuffer.data(idx_eri_1_df + 8);

    auto g_xx_zzz_1 = pbuffer.data(idx_eri_1_df + 9);

    auto g_yy_xxx_1 = pbuffer.data(idx_eri_1_df + 30);

    auto g_yy_xxy_1 = pbuffer.data(idx_eri_1_df + 31);

    auto g_yy_xxz_1 = pbuffer.data(idx_eri_1_df + 32);

    auto g_yy_xyy_1 = pbuffer.data(idx_eri_1_df + 33);

    auto g_yy_xyz_1 = pbuffer.data(idx_eri_1_df + 34);

    auto g_yy_xzz_1 = pbuffer.data(idx_eri_1_df + 35);

    auto g_yy_yyy_1 = pbuffer.data(idx_eri_1_df + 36);

    auto g_yy_yyz_1 = pbuffer.data(idx_eri_1_df + 37);

    auto g_yy_yzz_1 = pbuffer.data(idx_eri_1_df + 38);

    auto g_yy_zzz_1 = pbuffer.data(idx_eri_1_df + 39);

    auto g_zz_xxx_1 = pbuffer.data(idx_eri_1_df + 50);

    auto g_zz_xxy_1 = pbuffer.data(idx_eri_1_df + 51);

    auto g_zz_xxz_1 = pbuffer.data(idx_eri_1_df + 52);

    auto g_zz_xyy_1 = pbuffer.data(idx_eri_1_df + 53);

    auto g_zz_xyz_1 = pbuffer.data(idx_eri_1_df + 54);

    auto g_zz_xzz_1 = pbuffer.data(idx_eri_1_df + 55);

    auto g_zz_yyy_1 = pbuffer.data(idx_eri_1_df + 56);

    auto g_zz_yyz_1 = pbuffer.data(idx_eri_1_df + 57);

    auto g_zz_yzz_1 = pbuffer.data(idx_eri_1_df + 58);

    auto g_zz_zzz_1 = pbuffer.data(idx_eri_1_df + 59);

    // Set up components of auxiliary buffer : FD

    auto g_xxx_xx_1 = pbuffer.data(idx_eri_1_fd);

    auto g_xxx_xy_1 = pbuffer.data(idx_eri_1_fd + 1);

    auto g_xxx_xz_1 = pbuffer.data(idx_eri_1_fd + 2);

    auto g_xxx_yy_1 = pbuffer.data(idx_eri_1_fd + 3);

    auto g_xxx_yz_1 = pbuffer.data(idx_eri_1_fd + 4);

    auto g_xxx_zz_1 = pbuffer.data(idx_eri_1_fd + 5);

    auto g_xxz_xz_1 = pbuffer.data(idx_eri_1_fd + 14);

    auto g_xxz_yz_1 = pbuffer.data(idx_eri_1_fd + 16);

    auto g_xxz_zz_1 = pbuffer.data(idx_eri_1_fd + 17);

    auto g_xyy_xy_1 = pbuffer.data(idx_eri_1_fd + 19);

    auto g_xyy_yy_1 = pbuffer.data(idx_eri_1_fd + 21);

    auto g_xyy_yz_1 = pbuffer.data(idx_eri_1_fd + 22);

    auto g_xzz_xz_1 = pbuffer.data(idx_eri_1_fd + 32);

    auto g_xzz_yz_1 = pbuffer.data(idx_eri_1_fd + 34);

    auto g_xzz_zz_1 = pbuffer.data(idx_eri_1_fd + 35);

    auto g_yyy_xx_1 = pbuffer.data(idx_eri_1_fd + 36);

    auto g_yyy_xy_1 = pbuffer.data(idx_eri_1_fd + 37);

    auto g_yyy_xz_1 = pbuffer.data(idx_eri_1_fd + 38);

    auto g_yyy_yy_1 = pbuffer.data(idx_eri_1_fd + 39);

    auto g_yyy_yz_1 = pbuffer.data(idx_eri_1_fd + 40);

    auto g_yyy_zz_1 = pbuffer.data(idx_eri_1_fd + 41);

    auto g_yyz_xz_1 = pbuffer.data(idx_eri_1_fd + 44);

    auto g_yyz_yz_1 = pbuffer.data(idx_eri_1_fd + 46);

    auto g_yyz_zz_1 = pbuffer.data(idx_eri_1_fd + 47);

    auto g_yzz_xy_1 = pbuffer.data(idx_eri_1_fd + 49);

    auto g_yzz_xz_1 = pbuffer.data(idx_eri_1_fd + 50);

    auto g_yzz_yy_1 = pbuffer.data(idx_eri_1_fd + 51);

    auto g_yzz_yz_1 = pbuffer.data(idx_eri_1_fd + 52);

    auto g_yzz_zz_1 = pbuffer.data(idx_eri_1_fd + 53);

    auto g_zzz_xx_1 = pbuffer.data(idx_eri_1_fd + 54);

    auto g_zzz_xy_1 = pbuffer.data(idx_eri_1_fd + 55);

    auto g_zzz_xz_1 = pbuffer.data(idx_eri_1_fd + 56);

    auto g_zzz_yy_1 = pbuffer.data(idx_eri_1_fd + 57);

    auto g_zzz_yz_1 = pbuffer.data(idx_eri_1_fd + 58);

    auto g_zzz_zz_1 = pbuffer.data(idx_eri_1_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto g_xxx_xxx_1 = pbuffer.data(idx_eri_1_ff);

    auto g_xxx_xxy_1 = pbuffer.data(idx_eri_1_ff + 1);

    auto g_xxx_xxz_1 = pbuffer.data(idx_eri_1_ff + 2);

    auto g_xxx_xyy_1 = pbuffer.data(idx_eri_1_ff + 3);

    auto g_xxx_xyz_1 = pbuffer.data(idx_eri_1_ff + 4);

    auto g_xxx_xzz_1 = pbuffer.data(idx_eri_1_ff + 5);

    auto g_xxx_yyy_1 = pbuffer.data(idx_eri_1_ff + 6);

    auto g_xxx_yyz_1 = pbuffer.data(idx_eri_1_ff + 7);

    auto g_xxx_yzz_1 = pbuffer.data(idx_eri_1_ff + 8);

    auto g_xxx_zzz_1 = pbuffer.data(idx_eri_1_ff + 9);

    auto g_xxy_xxx_1 = pbuffer.data(idx_eri_1_ff + 10);

    auto g_xxy_xxy_1 = pbuffer.data(idx_eri_1_ff + 11);

    auto g_xxy_xxz_1 = pbuffer.data(idx_eri_1_ff + 12);

    auto g_xxy_xyy_1 = pbuffer.data(idx_eri_1_ff + 13);

    auto g_xxy_xzz_1 = pbuffer.data(idx_eri_1_ff + 15);

    auto g_xxy_yyy_1 = pbuffer.data(idx_eri_1_ff + 16);

    auto g_xxz_xxx_1 = pbuffer.data(idx_eri_1_ff + 20);

    auto g_xxz_xxy_1 = pbuffer.data(idx_eri_1_ff + 21);

    auto g_xxz_xxz_1 = pbuffer.data(idx_eri_1_ff + 22);

    auto g_xxz_xyy_1 = pbuffer.data(idx_eri_1_ff + 23);

    auto g_xxz_xyz_1 = pbuffer.data(idx_eri_1_ff + 24);

    auto g_xxz_xzz_1 = pbuffer.data(idx_eri_1_ff + 25);

    auto g_xxz_yyz_1 = pbuffer.data(idx_eri_1_ff + 27);

    auto g_xxz_yzz_1 = pbuffer.data(idx_eri_1_ff + 28);

    auto g_xxz_zzz_1 = pbuffer.data(idx_eri_1_ff + 29);

    auto g_xyy_xxx_1 = pbuffer.data(idx_eri_1_ff + 30);

    auto g_xyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 31);

    auto g_xyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 33);

    auto g_xyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 34);

    auto g_xyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 36);

    auto g_xyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 37);

    auto g_xyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 38);

    auto g_xyy_zzz_1 = pbuffer.data(idx_eri_1_ff + 39);

    auto g_xzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 50);

    auto g_xzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 52);

    auto g_xzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 54);

    auto g_xzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 55);

    auto g_xzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 56);

    auto g_xzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 57);

    auto g_xzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 58);

    auto g_xzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 59);

    auto g_yyy_xxx_1 = pbuffer.data(idx_eri_1_ff + 60);

    auto g_yyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 61);

    auto g_yyy_xxz_1 = pbuffer.data(idx_eri_1_ff + 62);

    auto g_yyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 63);

    auto g_yyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 64);

    auto g_yyy_xzz_1 = pbuffer.data(idx_eri_1_ff + 65);

    auto g_yyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 66);

    auto g_yyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 67);

    auto g_yyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 68);

    auto g_yyy_zzz_1 = pbuffer.data(idx_eri_1_ff + 69);

    auto g_yyz_xxy_1 = pbuffer.data(idx_eri_1_ff + 71);

    auto g_yyz_xxz_1 = pbuffer.data(idx_eri_1_ff + 72);

    auto g_yyz_xyy_1 = pbuffer.data(idx_eri_1_ff + 73);

    auto g_yyz_xyz_1 = pbuffer.data(idx_eri_1_ff + 74);

    auto g_yyz_xzz_1 = pbuffer.data(idx_eri_1_ff + 75);

    auto g_yyz_yyy_1 = pbuffer.data(idx_eri_1_ff + 76);

    auto g_yyz_yyz_1 = pbuffer.data(idx_eri_1_ff + 77);

    auto g_yyz_yzz_1 = pbuffer.data(idx_eri_1_ff + 78);

    auto g_yyz_zzz_1 = pbuffer.data(idx_eri_1_ff + 79);

    auto g_yzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 80);

    auto g_yzz_xxy_1 = pbuffer.data(idx_eri_1_ff + 81);

    auto g_yzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 82);

    auto g_yzz_xyy_1 = pbuffer.data(idx_eri_1_ff + 83);

    auto g_yzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 84);

    auto g_yzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 85);

    auto g_yzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 86);

    auto g_yzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 87);

    auto g_yzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 88);

    auto g_yzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 89);

    auto g_zzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 90);

    auto g_zzz_xxy_1 = pbuffer.data(idx_eri_1_ff + 91);

    auto g_zzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 92);

    auto g_zzz_xyy_1 = pbuffer.data(idx_eri_1_ff + 93);

    auto g_zzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 94);

    auto g_zzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 95);

    auto g_zzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 96);

    auto g_zzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 97);

    auto g_zzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 98);

    auto g_zzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 99);

    // Set up 0-10 components of targeted buffer : GF

    auto g_xxxx_xxx_0 = pbuffer.data(idx_eri_0_gf);

    auto g_xxxx_xxy_0 = pbuffer.data(idx_eri_0_gf + 1);

    auto g_xxxx_xxz_0 = pbuffer.data(idx_eri_0_gf + 2);

    auto g_xxxx_xyy_0 = pbuffer.data(idx_eri_0_gf + 3);

    auto g_xxxx_xyz_0 = pbuffer.data(idx_eri_0_gf + 4);

    auto g_xxxx_xzz_0 = pbuffer.data(idx_eri_0_gf + 5);

    auto g_xxxx_yyy_0 = pbuffer.data(idx_eri_0_gf + 6);

    auto g_xxxx_yyz_0 = pbuffer.data(idx_eri_0_gf + 7);

    auto g_xxxx_yzz_0 = pbuffer.data(idx_eri_0_gf + 8);

    auto g_xxxx_zzz_0 = pbuffer.data(idx_eri_0_gf + 9);

    #pragma omp simd aligned(g_xx_xxx_0, g_xx_xxx_1, g_xx_xxy_0, g_xx_xxy_1, g_xx_xxz_0, g_xx_xxz_1, g_xx_xyy_0, g_xx_xyy_1, g_xx_xyz_0, g_xx_xyz_1, g_xx_xzz_0, g_xx_xzz_1, g_xx_yyy_0, g_xx_yyy_1, g_xx_yyz_0, g_xx_yyz_1, g_xx_yzz_0, g_xx_yzz_1, g_xx_zzz_0, g_xx_zzz_1, g_xxx_xx_1, g_xxx_xxx_1, g_xxx_xxy_1, g_xxx_xxz_1, g_xxx_xy_1, g_xxx_xyy_1, g_xxx_xyz_1, g_xxx_xz_1, g_xxx_xzz_1, g_xxx_yy_1, g_xxx_yyy_1, g_xxx_yyz_1, g_xxx_yz_1, g_xxx_yzz_1, g_xxx_zz_1, g_xxx_zzz_1, g_xxxx_xxx_0, g_xxxx_xxy_0, g_xxxx_xxz_0, g_xxxx_xyy_0, g_xxxx_xyz_0, g_xxxx_xzz_0, g_xxxx_yyy_0, g_xxxx_yyz_0, g_xxxx_yzz_0, g_xxxx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_xxx_0[i] = 3.0 * g_xx_xxx_0[i] * fbe_0 - 3.0 * g_xx_xxx_1[i] * fz_be_0 + 3.0 * g_xxx_xx_1[i] * fe_0 + g_xxx_xxx_1[i] * pa_x[i];

        g_xxxx_xxy_0[i] = 3.0 * g_xx_xxy_0[i] * fbe_0 - 3.0 * g_xx_xxy_1[i] * fz_be_0 + 2.0 * g_xxx_xy_1[i] * fe_0 + g_xxx_xxy_1[i] * pa_x[i];

        g_xxxx_xxz_0[i] = 3.0 * g_xx_xxz_0[i] * fbe_0 - 3.0 * g_xx_xxz_1[i] * fz_be_0 + 2.0 * g_xxx_xz_1[i] * fe_0 + g_xxx_xxz_1[i] * pa_x[i];

        g_xxxx_xyy_0[i] = 3.0 * g_xx_xyy_0[i] * fbe_0 - 3.0 * g_xx_xyy_1[i] * fz_be_0 + g_xxx_yy_1[i] * fe_0 + g_xxx_xyy_1[i] * pa_x[i];

        g_xxxx_xyz_0[i] = 3.0 * g_xx_xyz_0[i] * fbe_0 - 3.0 * g_xx_xyz_1[i] * fz_be_0 + g_xxx_yz_1[i] * fe_0 + g_xxx_xyz_1[i] * pa_x[i];

        g_xxxx_xzz_0[i] = 3.0 * g_xx_xzz_0[i] * fbe_0 - 3.0 * g_xx_xzz_1[i] * fz_be_0 + g_xxx_zz_1[i] * fe_0 + g_xxx_xzz_1[i] * pa_x[i];

        g_xxxx_yyy_0[i] = 3.0 * g_xx_yyy_0[i] * fbe_0 - 3.0 * g_xx_yyy_1[i] * fz_be_0 + g_xxx_yyy_1[i] * pa_x[i];

        g_xxxx_yyz_0[i] = 3.0 * g_xx_yyz_0[i] * fbe_0 - 3.0 * g_xx_yyz_1[i] * fz_be_0 + g_xxx_yyz_1[i] * pa_x[i];

        g_xxxx_yzz_0[i] = 3.0 * g_xx_yzz_0[i] * fbe_0 - 3.0 * g_xx_yzz_1[i] * fz_be_0 + g_xxx_yzz_1[i] * pa_x[i];

        g_xxxx_zzz_0[i] = 3.0 * g_xx_zzz_0[i] * fbe_0 - 3.0 * g_xx_zzz_1[i] * fz_be_0 + g_xxx_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto g_xxxy_xxx_0 = pbuffer.data(idx_eri_0_gf + 10);

    auto g_xxxy_xxy_0 = pbuffer.data(idx_eri_0_gf + 11);

    auto g_xxxy_xxz_0 = pbuffer.data(idx_eri_0_gf + 12);

    auto g_xxxy_xyy_0 = pbuffer.data(idx_eri_0_gf + 13);

    auto g_xxxy_xyz_0 = pbuffer.data(idx_eri_0_gf + 14);

    auto g_xxxy_xzz_0 = pbuffer.data(idx_eri_0_gf + 15);

    auto g_xxxy_yyy_0 = pbuffer.data(idx_eri_0_gf + 16);

    auto g_xxxy_yyz_0 = pbuffer.data(idx_eri_0_gf + 17);

    auto g_xxxy_yzz_0 = pbuffer.data(idx_eri_0_gf + 18);

    auto g_xxxy_zzz_0 = pbuffer.data(idx_eri_0_gf + 19);

    #pragma omp simd aligned(g_xxx_xx_1, g_xxx_xxx_1, g_xxx_xxy_1, g_xxx_xxz_1, g_xxx_xy_1, g_xxx_xyy_1, g_xxx_xyz_1, g_xxx_xz_1, g_xxx_xzz_1, g_xxx_yy_1, g_xxx_yyy_1, g_xxx_yyz_1, g_xxx_yz_1, g_xxx_yzz_1, g_xxx_zz_1, g_xxx_zzz_1, g_xxxy_xxx_0, g_xxxy_xxy_0, g_xxxy_xxz_0, g_xxxy_xyy_0, g_xxxy_xyz_0, g_xxxy_xzz_0, g_xxxy_yyy_0, g_xxxy_yyz_0, g_xxxy_yzz_0, g_xxxy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_xxx_0[i] = g_xxx_xxx_1[i] * pa_y[i];

        g_xxxy_xxy_0[i] = g_xxx_xx_1[i] * fe_0 + g_xxx_xxy_1[i] * pa_y[i];

        g_xxxy_xxz_0[i] = g_xxx_xxz_1[i] * pa_y[i];

        g_xxxy_xyy_0[i] = 2.0 * g_xxx_xy_1[i] * fe_0 + g_xxx_xyy_1[i] * pa_y[i];

        g_xxxy_xyz_0[i] = g_xxx_xz_1[i] * fe_0 + g_xxx_xyz_1[i] * pa_y[i];

        g_xxxy_xzz_0[i] = g_xxx_xzz_1[i] * pa_y[i];

        g_xxxy_yyy_0[i] = 3.0 * g_xxx_yy_1[i] * fe_0 + g_xxx_yyy_1[i] * pa_y[i];

        g_xxxy_yyz_0[i] = 2.0 * g_xxx_yz_1[i] * fe_0 + g_xxx_yyz_1[i] * pa_y[i];

        g_xxxy_yzz_0[i] = g_xxx_zz_1[i] * fe_0 + g_xxx_yzz_1[i] * pa_y[i];

        g_xxxy_zzz_0[i] = g_xxx_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto g_xxxz_xxx_0 = pbuffer.data(idx_eri_0_gf + 20);

    auto g_xxxz_xxy_0 = pbuffer.data(idx_eri_0_gf + 21);

    auto g_xxxz_xxz_0 = pbuffer.data(idx_eri_0_gf + 22);

    auto g_xxxz_xyy_0 = pbuffer.data(idx_eri_0_gf + 23);

    auto g_xxxz_xyz_0 = pbuffer.data(idx_eri_0_gf + 24);

    auto g_xxxz_xzz_0 = pbuffer.data(idx_eri_0_gf + 25);

    auto g_xxxz_yyy_0 = pbuffer.data(idx_eri_0_gf + 26);

    auto g_xxxz_yyz_0 = pbuffer.data(idx_eri_0_gf + 27);

    auto g_xxxz_yzz_0 = pbuffer.data(idx_eri_0_gf + 28);

    auto g_xxxz_zzz_0 = pbuffer.data(idx_eri_0_gf + 29);

    #pragma omp simd aligned(g_xxx_xx_1, g_xxx_xxx_1, g_xxx_xxy_1, g_xxx_xxz_1, g_xxx_xy_1, g_xxx_xyy_1, g_xxx_xyz_1, g_xxx_xz_1, g_xxx_xzz_1, g_xxx_yy_1, g_xxx_yyy_1, g_xxx_yyz_1, g_xxx_yz_1, g_xxx_yzz_1, g_xxx_zz_1, g_xxx_zzz_1, g_xxxz_xxx_0, g_xxxz_xxy_0, g_xxxz_xxz_0, g_xxxz_xyy_0, g_xxxz_xyz_0, g_xxxz_xzz_0, g_xxxz_yyy_0, g_xxxz_yyz_0, g_xxxz_yzz_0, g_xxxz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_xxx_0[i] = g_xxx_xxx_1[i] * pa_z[i];

        g_xxxz_xxy_0[i] = g_xxx_xxy_1[i] * pa_z[i];

        g_xxxz_xxz_0[i] = g_xxx_xx_1[i] * fe_0 + g_xxx_xxz_1[i] * pa_z[i];

        g_xxxz_xyy_0[i] = g_xxx_xyy_1[i] * pa_z[i];

        g_xxxz_xyz_0[i] = g_xxx_xy_1[i] * fe_0 + g_xxx_xyz_1[i] * pa_z[i];

        g_xxxz_xzz_0[i] = 2.0 * g_xxx_xz_1[i] * fe_0 + g_xxx_xzz_1[i] * pa_z[i];

        g_xxxz_yyy_0[i] = g_xxx_yyy_1[i] * pa_z[i];

        g_xxxz_yyz_0[i] = g_xxx_yy_1[i] * fe_0 + g_xxx_yyz_1[i] * pa_z[i];

        g_xxxz_yzz_0[i] = 2.0 * g_xxx_yz_1[i] * fe_0 + g_xxx_yzz_1[i] * pa_z[i];

        g_xxxz_zzz_0[i] = 3.0 * g_xxx_zz_1[i] * fe_0 + g_xxx_zzz_1[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto g_xxyy_xxx_0 = pbuffer.data(idx_eri_0_gf + 30);

    auto g_xxyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 31);

    auto g_xxyy_xxz_0 = pbuffer.data(idx_eri_0_gf + 32);

    auto g_xxyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 33);

    auto g_xxyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 34);

    auto g_xxyy_xzz_0 = pbuffer.data(idx_eri_0_gf + 35);

    auto g_xxyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 36);

    auto g_xxyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 37);

    auto g_xxyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 38);

    auto g_xxyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 39);

    #pragma omp simd aligned(g_xx_xxx_0, g_xx_xxx_1, g_xx_xxz_0, g_xx_xxz_1, g_xx_xzz_0, g_xx_xzz_1, g_xxy_xxx_1, g_xxy_xxz_1, g_xxy_xzz_1, g_xxyy_xxx_0, g_xxyy_xxy_0, g_xxyy_xxz_0, g_xxyy_xyy_0, g_xxyy_xyz_0, g_xxyy_xzz_0, g_xxyy_yyy_0, g_xxyy_yyz_0, g_xxyy_yzz_0, g_xxyy_zzz_0, g_xyy_xxy_1, g_xyy_xy_1, g_xyy_xyy_1, g_xyy_xyz_1, g_xyy_yy_1, g_xyy_yyy_1, g_xyy_yyz_1, g_xyy_yz_1, g_xyy_yzz_1, g_xyy_zzz_1, g_yy_xxy_0, g_yy_xxy_1, g_yy_xyy_0, g_yy_xyy_1, g_yy_xyz_0, g_yy_xyz_1, g_yy_yyy_0, g_yy_yyy_1, g_yy_yyz_0, g_yy_yyz_1, g_yy_yzz_0, g_yy_yzz_1, g_yy_zzz_0, g_yy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyy_xxx_0[i] = g_xx_xxx_0[i] * fbe_0 - g_xx_xxx_1[i] * fz_be_0 + g_xxy_xxx_1[i] * pa_y[i];

        g_xxyy_xxy_0[i] = g_yy_xxy_0[i] * fbe_0 - g_yy_xxy_1[i] * fz_be_0 + 2.0 * g_xyy_xy_1[i] * fe_0 + g_xyy_xxy_1[i] * pa_x[i];

        g_xxyy_xxz_0[i] = g_xx_xxz_0[i] * fbe_0 - g_xx_xxz_1[i] * fz_be_0 + g_xxy_xxz_1[i] * pa_y[i];

        g_xxyy_xyy_0[i] = g_yy_xyy_0[i] * fbe_0 - g_yy_xyy_1[i] * fz_be_0 + g_xyy_yy_1[i] * fe_0 + g_xyy_xyy_1[i] * pa_x[i];

        g_xxyy_xyz_0[i] = g_yy_xyz_0[i] * fbe_0 - g_yy_xyz_1[i] * fz_be_0 + g_xyy_yz_1[i] * fe_0 + g_xyy_xyz_1[i] * pa_x[i];

        g_xxyy_xzz_0[i] = g_xx_xzz_0[i] * fbe_0 - g_xx_xzz_1[i] * fz_be_0 + g_xxy_xzz_1[i] * pa_y[i];

        g_xxyy_yyy_0[i] = g_yy_yyy_0[i] * fbe_0 - g_yy_yyy_1[i] * fz_be_0 + g_xyy_yyy_1[i] * pa_x[i];

        g_xxyy_yyz_0[i] = g_yy_yyz_0[i] * fbe_0 - g_yy_yyz_1[i] * fz_be_0 + g_xyy_yyz_1[i] * pa_x[i];

        g_xxyy_yzz_0[i] = g_yy_yzz_0[i] * fbe_0 - g_yy_yzz_1[i] * fz_be_0 + g_xyy_yzz_1[i] * pa_x[i];

        g_xxyy_zzz_0[i] = g_yy_zzz_0[i] * fbe_0 - g_yy_zzz_1[i] * fz_be_0 + g_xyy_zzz_1[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto g_xxyz_xxx_0 = pbuffer.data(idx_eri_0_gf + 40);

    auto g_xxyz_xxy_0 = pbuffer.data(idx_eri_0_gf + 41);

    auto g_xxyz_xxz_0 = pbuffer.data(idx_eri_0_gf + 42);

    auto g_xxyz_xyy_0 = pbuffer.data(idx_eri_0_gf + 43);

    auto g_xxyz_xyz_0 = pbuffer.data(idx_eri_0_gf + 44);

    auto g_xxyz_xzz_0 = pbuffer.data(idx_eri_0_gf + 45);

    auto g_xxyz_yyy_0 = pbuffer.data(idx_eri_0_gf + 46);

    auto g_xxyz_yyz_0 = pbuffer.data(idx_eri_0_gf + 47);

    auto g_xxyz_yzz_0 = pbuffer.data(idx_eri_0_gf + 48);

    auto g_xxyz_zzz_0 = pbuffer.data(idx_eri_0_gf + 49);

    #pragma omp simd aligned(g_xxy_xxy_1, g_xxy_xyy_1, g_xxy_yyy_1, g_xxyz_xxx_0, g_xxyz_xxy_0, g_xxyz_xxz_0, g_xxyz_xyy_0, g_xxyz_xyz_0, g_xxyz_xzz_0, g_xxyz_yyy_0, g_xxyz_yyz_0, g_xxyz_yzz_0, g_xxyz_zzz_0, g_xxz_xxx_1, g_xxz_xxz_1, g_xxz_xyz_1, g_xxz_xz_1, g_xxz_xzz_1, g_xxz_yyz_1, g_xxz_yz_1, g_xxz_yzz_1, g_xxz_zz_1, g_xxz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyz_xxx_0[i] = g_xxz_xxx_1[i] * pa_y[i];

        g_xxyz_xxy_0[i] = g_xxy_xxy_1[i] * pa_z[i];

        g_xxyz_xxz_0[i] = g_xxz_xxz_1[i] * pa_y[i];

        g_xxyz_xyy_0[i] = g_xxy_xyy_1[i] * pa_z[i];

        g_xxyz_xyz_0[i] = g_xxz_xz_1[i] * fe_0 + g_xxz_xyz_1[i] * pa_y[i];

        g_xxyz_xzz_0[i] = g_xxz_xzz_1[i] * pa_y[i];

        g_xxyz_yyy_0[i] = g_xxy_yyy_1[i] * pa_z[i];

        g_xxyz_yyz_0[i] = 2.0 * g_xxz_yz_1[i] * fe_0 + g_xxz_yyz_1[i] * pa_y[i];

        g_xxyz_yzz_0[i] = g_xxz_zz_1[i] * fe_0 + g_xxz_yzz_1[i] * pa_y[i];

        g_xxyz_zzz_0[i] = g_xxz_zzz_1[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto g_xxzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 50);

    auto g_xxzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 51);

    auto g_xxzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 52);

    auto g_xxzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 53);

    auto g_xxzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 54);

    auto g_xxzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 55);

    auto g_xxzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 56);

    auto g_xxzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 57);

    auto g_xxzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 58);

    auto g_xxzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 59);

    #pragma omp simd aligned(g_xx_xxx_0, g_xx_xxx_1, g_xx_xxy_0, g_xx_xxy_1, g_xx_xyy_0, g_xx_xyy_1, g_xxz_xxx_1, g_xxz_xxy_1, g_xxz_xyy_1, g_xxzz_xxx_0, g_xxzz_xxy_0, g_xxzz_xxz_0, g_xxzz_xyy_0, g_xxzz_xyz_0, g_xxzz_xzz_0, g_xxzz_yyy_0, g_xxzz_yyz_0, g_xxzz_yzz_0, g_xxzz_zzz_0, g_xzz_xxz_1, g_xzz_xyz_1, g_xzz_xz_1, g_xzz_xzz_1, g_xzz_yyy_1, g_xzz_yyz_1, g_xzz_yz_1, g_xzz_yzz_1, g_xzz_zz_1, g_xzz_zzz_1, g_zz_xxz_0, g_zz_xxz_1, g_zz_xyz_0, g_zz_xyz_1, g_zz_xzz_0, g_zz_xzz_1, g_zz_yyy_0, g_zz_yyy_1, g_zz_yyz_0, g_zz_yyz_1, g_zz_yzz_0, g_zz_yzz_1, g_zz_zzz_0, g_zz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzz_xxx_0[i] = g_xx_xxx_0[i] * fbe_0 - g_xx_xxx_1[i] * fz_be_0 + g_xxz_xxx_1[i] * pa_z[i];

        g_xxzz_xxy_0[i] = g_xx_xxy_0[i] * fbe_0 - g_xx_xxy_1[i] * fz_be_0 + g_xxz_xxy_1[i] * pa_z[i];

        g_xxzz_xxz_0[i] = g_zz_xxz_0[i] * fbe_0 - g_zz_xxz_1[i] * fz_be_0 + 2.0 * g_xzz_xz_1[i] * fe_0 + g_xzz_xxz_1[i] * pa_x[i];

        g_xxzz_xyy_0[i] = g_xx_xyy_0[i] * fbe_0 - g_xx_xyy_1[i] * fz_be_0 + g_xxz_xyy_1[i] * pa_z[i];

        g_xxzz_xyz_0[i] = g_zz_xyz_0[i] * fbe_0 - g_zz_xyz_1[i] * fz_be_0 + g_xzz_yz_1[i] * fe_0 + g_xzz_xyz_1[i] * pa_x[i];

        g_xxzz_xzz_0[i] = g_zz_xzz_0[i] * fbe_0 - g_zz_xzz_1[i] * fz_be_0 + g_xzz_zz_1[i] * fe_0 + g_xzz_xzz_1[i] * pa_x[i];

        g_xxzz_yyy_0[i] = g_zz_yyy_0[i] * fbe_0 - g_zz_yyy_1[i] * fz_be_0 + g_xzz_yyy_1[i] * pa_x[i];

        g_xxzz_yyz_0[i] = g_zz_yyz_0[i] * fbe_0 - g_zz_yyz_1[i] * fz_be_0 + g_xzz_yyz_1[i] * pa_x[i];

        g_xxzz_yzz_0[i] = g_zz_yzz_0[i] * fbe_0 - g_zz_yzz_1[i] * fz_be_0 + g_xzz_yzz_1[i] * pa_x[i];

        g_xxzz_zzz_0[i] = g_zz_zzz_0[i] * fbe_0 - g_zz_zzz_1[i] * fz_be_0 + g_xzz_zzz_1[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto g_xyyy_xxx_0 = pbuffer.data(idx_eri_0_gf + 60);

    auto g_xyyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 61);

    auto g_xyyy_xxz_0 = pbuffer.data(idx_eri_0_gf + 62);

    auto g_xyyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 63);

    auto g_xyyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 64);

    auto g_xyyy_xzz_0 = pbuffer.data(idx_eri_0_gf + 65);

    auto g_xyyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 66);

    auto g_xyyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 67);

    auto g_xyyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 68);

    auto g_xyyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 69);

    #pragma omp simd aligned(g_xyyy_xxx_0, g_xyyy_xxy_0, g_xyyy_xxz_0, g_xyyy_xyy_0, g_xyyy_xyz_0, g_xyyy_xzz_0, g_xyyy_yyy_0, g_xyyy_yyz_0, g_xyyy_yzz_0, g_xyyy_zzz_0, g_yyy_xx_1, g_yyy_xxx_1, g_yyy_xxy_1, g_yyy_xxz_1, g_yyy_xy_1, g_yyy_xyy_1, g_yyy_xyz_1, g_yyy_xz_1, g_yyy_xzz_1, g_yyy_yy_1, g_yyy_yyy_1, g_yyy_yyz_1, g_yyy_yz_1, g_yyy_yzz_1, g_yyy_zz_1, g_yyy_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_xxx_0[i] = 3.0 * g_yyy_xx_1[i] * fe_0 + g_yyy_xxx_1[i] * pa_x[i];

        g_xyyy_xxy_0[i] = 2.0 * g_yyy_xy_1[i] * fe_0 + g_yyy_xxy_1[i] * pa_x[i];

        g_xyyy_xxz_0[i] = 2.0 * g_yyy_xz_1[i] * fe_0 + g_yyy_xxz_1[i] * pa_x[i];

        g_xyyy_xyy_0[i] = g_yyy_yy_1[i] * fe_0 + g_yyy_xyy_1[i] * pa_x[i];

        g_xyyy_xyz_0[i] = g_yyy_yz_1[i] * fe_0 + g_yyy_xyz_1[i] * pa_x[i];

        g_xyyy_xzz_0[i] = g_yyy_zz_1[i] * fe_0 + g_yyy_xzz_1[i] * pa_x[i];

        g_xyyy_yyy_0[i] = g_yyy_yyy_1[i] * pa_x[i];

        g_xyyy_yyz_0[i] = g_yyy_yyz_1[i] * pa_x[i];

        g_xyyy_yzz_0[i] = g_yyy_yzz_1[i] * pa_x[i];

        g_xyyy_zzz_0[i] = g_yyy_zzz_1[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto g_xyyz_xxx_0 = pbuffer.data(idx_eri_0_gf + 70);

    auto g_xyyz_xxy_0 = pbuffer.data(idx_eri_0_gf + 71);

    auto g_xyyz_xxz_0 = pbuffer.data(idx_eri_0_gf + 72);

    auto g_xyyz_xyy_0 = pbuffer.data(idx_eri_0_gf + 73);

    auto g_xyyz_xyz_0 = pbuffer.data(idx_eri_0_gf + 74);

    auto g_xyyz_xzz_0 = pbuffer.data(idx_eri_0_gf + 75);

    auto g_xyyz_yyy_0 = pbuffer.data(idx_eri_0_gf + 76);

    auto g_xyyz_yyz_0 = pbuffer.data(idx_eri_0_gf + 77);

    auto g_xyyz_yzz_0 = pbuffer.data(idx_eri_0_gf + 78);

    auto g_xyyz_zzz_0 = pbuffer.data(idx_eri_0_gf + 79);

    #pragma omp simd aligned(g_xyy_xxx_1, g_xyy_xxy_1, g_xyy_xyy_1, g_xyyz_xxx_0, g_xyyz_xxy_0, g_xyyz_xxz_0, g_xyyz_xyy_0, g_xyyz_xyz_0, g_xyyz_xzz_0, g_xyyz_yyy_0, g_xyyz_yyz_0, g_xyyz_yzz_0, g_xyyz_zzz_0, g_yyz_xxz_1, g_yyz_xyz_1, g_yyz_xz_1, g_yyz_xzz_1, g_yyz_yyy_1, g_yyz_yyz_1, g_yyz_yz_1, g_yyz_yzz_1, g_yyz_zz_1, g_yyz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyz_xxx_0[i] = g_xyy_xxx_1[i] * pa_z[i];

        g_xyyz_xxy_0[i] = g_xyy_xxy_1[i] * pa_z[i];

        g_xyyz_xxz_0[i] = 2.0 * g_yyz_xz_1[i] * fe_0 + g_yyz_xxz_1[i] * pa_x[i];

        g_xyyz_xyy_0[i] = g_xyy_xyy_1[i] * pa_z[i];

        g_xyyz_xyz_0[i] = g_yyz_yz_1[i] * fe_0 + g_yyz_xyz_1[i] * pa_x[i];

        g_xyyz_xzz_0[i] = g_yyz_zz_1[i] * fe_0 + g_yyz_xzz_1[i] * pa_x[i];

        g_xyyz_yyy_0[i] = g_yyz_yyy_1[i] * pa_x[i];

        g_xyyz_yyz_0[i] = g_yyz_yyz_1[i] * pa_x[i];

        g_xyyz_yzz_0[i] = g_yyz_yzz_1[i] * pa_x[i];

        g_xyyz_zzz_0[i] = g_yyz_zzz_1[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto g_xyzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 80);

    auto g_xyzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 81);

    auto g_xyzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 82);

    auto g_xyzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 83);

    auto g_xyzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 84);

    auto g_xyzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 85);

    auto g_xyzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 86);

    auto g_xyzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 87);

    auto g_xyzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 88);

    auto g_xyzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 89);

    #pragma omp simd aligned(g_xyzz_xxx_0, g_xyzz_xxy_0, g_xyzz_xxz_0, g_xyzz_xyy_0, g_xyzz_xyz_0, g_xyzz_xzz_0, g_xyzz_yyy_0, g_xyzz_yyz_0, g_xyzz_yzz_0, g_xyzz_zzz_0, g_xzz_xxx_1, g_xzz_xxz_1, g_xzz_xzz_1, g_yzz_xxy_1, g_yzz_xy_1, g_yzz_xyy_1, g_yzz_xyz_1, g_yzz_yy_1, g_yzz_yyy_1, g_yzz_yyz_1, g_yzz_yz_1, g_yzz_yzz_1, g_yzz_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzz_xxx_0[i] = g_xzz_xxx_1[i] * pa_y[i];

        g_xyzz_xxy_0[i] = 2.0 * g_yzz_xy_1[i] * fe_0 + g_yzz_xxy_1[i] * pa_x[i];

        g_xyzz_xxz_0[i] = g_xzz_xxz_1[i] * pa_y[i];

        g_xyzz_xyy_0[i] = g_yzz_yy_1[i] * fe_0 + g_yzz_xyy_1[i] * pa_x[i];

        g_xyzz_xyz_0[i] = g_yzz_yz_1[i] * fe_0 + g_yzz_xyz_1[i] * pa_x[i];

        g_xyzz_xzz_0[i] = g_xzz_xzz_1[i] * pa_y[i];

        g_xyzz_yyy_0[i] = g_yzz_yyy_1[i] * pa_x[i];

        g_xyzz_yyz_0[i] = g_yzz_yyz_1[i] * pa_x[i];

        g_xyzz_yzz_0[i] = g_yzz_yzz_1[i] * pa_x[i];

        g_xyzz_zzz_0[i] = g_yzz_zzz_1[i] * pa_x[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto g_xzzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 90);

    auto g_xzzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 91);

    auto g_xzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 92);

    auto g_xzzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 93);

    auto g_xzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 94);

    auto g_xzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 95);

    auto g_xzzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 96);

    auto g_xzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 97);

    auto g_xzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 98);

    auto g_xzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 99);

    #pragma omp simd aligned(g_xzzz_xxx_0, g_xzzz_xxy_0, g_xzzz_xxz_0, g_xzzz_xyy_0, g_xzzz_xyz_0, g_xzzz_xzz_0, g_xzzz_yyy_0, g_xzzz_yyz_0, g_xzzz_yzz_0, g_xzzz_zzz_0, g_zzz_xx_1, g_zzz_xxx_1, g_zzz_xxy_1, g_zzz_xxz_1, g_zzz_xy_1, g_zzz_xyy_1, g_zzz_xyz_1, g_zzz_xz_1, g_zzz_xzz_1, g_zzz_yy_1, g_zzz_yyy_1, g_zzz_yyz_1, g_zzz_yz_1, g_zzz_yzz_1, g_zzz_zz_1, g_zzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_xxx_0[i] = 3.0 * g_zzz_xx_1[i] * fe_0 + g_zzz_xxx_1[i] * pa_x[i];

        g_xzzz_xxy_0[i] = 2.0 * g_zzz_xy_1[i] * fe_0 + g_zzz_xxy_1[i] * pa_x[i];

        g_xzzz_xxz_0[i] = 2.0 * g_zzz_xz_1[i] * fe_0 + g_zzz_xxz_1[i] * pa_x[i];

        g_xzzz_xyy_0[i] = g_zzz_yy_1[i] * fe_0 + g_zzz_xyy_1[i] * pa_x[i];

        g_xzzz_xyz_0[i] = g_zzz_yz_1[i] * fe_0 + g_zzz_xyz_1[i] * pa_x[i];

        g_xzzz_xzz_0[i] = g_zzz_zz_1[i] * fe_0 + g_zzz_xzz_1[i] * pa_x[i];

        g_xzzz_yyy_0[i] = g_zzz_yyy_1[i] * pa_x[i];

        g_xzzz_yyz_0[i] = g_zzz_yyz_1[i] * pa_x[i];

        g_xzzz_yzz_0[i] = g_zzz_yzz_1[i] * pa_x[i];

        g_xzzz_zzz_0[i] = g_zzz_zzz_1[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto g_yyyy_xxx_0 = pbuffer.data(idx_eri_0_gf + 100);

    auto g_yyyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 101);

    auto g_yyyy_xxz_0 = pbuffer.data(idx_eri_0_gf + 102);

    auto g_yyyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 103);

    auto g_yyyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 104);

    auto g_yyyy_xzz_0 = pbuffer.data(idx_eri_0_gf + 105);

    auto g_yyyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 106);

    auto g_yyyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 107);

    auto g_yyyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 108);

    auto g_yyyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 109);

    #pragma omp simd aligned(g_yy_xxx_0, g_yy_xxx_1, g_yy_xxy_0, g_yy_xxy_1, g_yy_xxz_0, g_yy_xxz_1, g_yy_xyy_0, g_yy_xyy_1, g_yy_xyz_0, g_yy_xyz_1, g_yy_xzz_0, g_yy_xzz_1, g_yy_yyy_0, g_yy_yyy_1, g_yy_yyz_0, g_yy_yyz_1, g_yy_yzz_0, g_yy_yzz_1, g_yy_zzz_0, g_yy_zzz_1, g_yyy_xx_1, g_yyy_xxx_1, g_yyy_xxy_1, g_yyy_xxz_1, g_yyy_xy_1, g_yyy_xyy_1, g_yyy_xyz_1, g_yyy_xz_1, g_yyy_xzz_1, g_yyy_yy_1, g_yyy_yyy_1, g_yyy_yyz_1, g_yyy_yz_1, g_yyy_yzz_1, g_yyy_zz_1, g_yyy_zzz_1, g_yyyy_xxx_0, g_yyyy_xxy_0, g_yyyy_xxz_0, g_yyyy_xyy_0, g_yyyy_xyz_0, g_yyyy_xzz_0, g_yyyy_yyy_0, g_yyyy_yyz_0, g_yyyy_yzz_0, g_yyyy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_xxx_0[i] = 3.0 * g_yy_xxx_0[i] * fbe_0 - 3.0 * g_yy_xxx_1[i] * fz_be_0 + g_yyy_xxx_1[i] * pa_y[i];

        g_yyyy_xxy_0[i] = 3.0 * g_yy_xxy_0[i] * fbe_0 - 3.0 * g_yy_xxy_1[i] * fz_be_0 + g_yyy_xx_1[i] * fe_0 + g_yyy_xxy_1[i] * pa_y[i];

        g_yyyy_xxz_0[i] = 3.0 * g_yy_xxz_0[i] * fbe_0 - 3.0 * g_yy_xxz_1[i] * fz_be_0 + g_yyy_xxz_1[i] * pa_y[i];

        g_yyyy_xyy_0[i] = 3.0 * g_yy_xyy_0[i] * fbe_0 - 3.0 * g_yy_xyy_1[i] * fz_be_0 + 2.0 * g_yyy_xy_1[i] * fe_0 + g_yyy_xyy_1[i] * pa_y[i];

        g_yyyy_xyz_0[i] = 3.0 * g_yy_xyz_0[i] * fbe_0 - 3.0 * g_yy_xyz_1[i] * fz_be_0 + g_yyy_xz_1[i] * fe_0 + g_yyy_xyz_1[i] * pa_y[i];

        g_yyyy_xzz_0[i] = 3.0 * g_yy_xzz_0[i] * fbe_0 - 3.0 * g_yy_xzz_1[i] * fz_be_0 + g_yyy_xzz_1[i] * pa_y[i];

        g_yyyy_yyy_0[i] = 3.0 * g_yy_yyy_0[i] * fbe_0 - 3.0 * g_yy_yyy_1[i] * fz_be_0 + 3.0 * g_yyy_yy_1[i] * fe_0 + g_yyy_yyy_1[i] * pa_y[i];

        g_yyyy_yyz_0[i] = 3.0 * g_yy_yyz_0[i] * fbe_0 - 3.0 * g_yy_yyz_1[i] * fz_be_0 + 2.0 * g_yyy_yz_1[i] * fe_0 + g_yyy_yyz_1[i] * pa_y[i];

        g_yyyy_yzz_0[i] = 3.0 * g_yy_yzz_0[i] * fbe_0 - 3.0 * g_yy_yzz_1[i] * fz_be_0 + g_yyy_zz_1[i] * fe_0 + g_yyy_yzz_1[i] * pa_y[i];

        g_yyyy_zzz_0[i] = 3.0 * g_yy_zzz_0[i] * fbe_0 - 3.0 * g_yy_zzz_1[i] * fz_be_0 + g_yyy_zzz_1[i] * pa_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto g_yyyz_xxx_0 = pbuffer.data(idx_eri_0_gf + 110);

    auto g_yyyz_xxy_0 = pbuffer.data(idx_eri_0_gf + 111);

    auto g_yyyz_xxz_0 = pbuffer.data(idx_eri_0_gf + 112);

    auto g_yyyz_xyy_0 = pbuffer.data(idx_eri_0_gf + 113);

    auto g_yyyz_xyz_0 = pbuffer.data(idx_eri_0_gf + 114);

    auto g_yyyz_xzz_0 = pbuffer.data(idx_eri_0_gf + 115);

    auto g_yyyz_yyy_0 = pbuffer.data(idx_eri_0_gf + 116);

    auto g_yyyz_yyz_0 = pbuffer.data(idx_eri_0_gf + 117);

    auto g_yyyz_yzz_0 = pbuffer.data(idx_eri_0_gf + 118);

    auto g_yyyz_zzz_0 = pbuffer.data(idx_eri_0_gf + 119);

    #pragma omp simd aligned(g_yyy_xx_1, g_yyy_xxx_1, g_yyy_xxy_1, g_yyy_xxz_1, g_yyy_xy_1, g_yyy_xyy_1, g_yyy_xyz_1, g_yyy_xz_1, g_yyy_xzz_1, g_yyy_yy_1, g_yyy_yyy_1, g_yyy_yyz_1, g_yyy_yz_1, g_yyy_yzz_1, g_yyy_zz_1, g_yyy_zzz_1, g_yyyz_xxx_0, g_yyyz_xxy_0, g_yyyz_xxz_0, g_yyyz_xyy_0, g_yyyz_xyz_0, g_yyyz_xzz_0, g_yyyz_yyy_0, g_yyyz_yyz_0, g_yyyz_yzz_0, g_yyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_xxx_0[i] = g_yyy_xxx_1[i] * pa_z[i];

        g_yyyz_xxy_0[i] = g_yyy_xxy_1[i] * pa_z[i];

        g_yyyz_xxz_0[i] = g_yyy_xx_1[i] * fe_0 + g_yyy_xxz_1[i] * pa_z[i];

        g_yyyz_xyy_0[i] = g_yyy_xyy_1[i] * pa_z[i];

        g_yyyz_xyz_0[i] = g_yyy_xy_1[i] * fe_0 + g_yyy_xyz_1[i] * pa_z[i];

        g_yyyz_xzz_0[i] = 2.0 * g_yyy_xz_1[i] * fe_0 + g_yyy_xzz_1[i] * pa_z[i];

        g_yyyz_yyy_0[i] = g_yyy_yyy_1[i] * pa_z[i];

        g_yyyz_yyz_0[i] = g_yyy_yy_1[i] * fe_0 + g_yyy_yyz_1[i] * pa_z[i];

        g_yyyz_yzz_0[i] = 2.0 * g_yyy_yz_1[i] * fe_0 + g_yyy_yzz_1[i] * pa_z[i];

        g_yyyz_zzz_0[i] = 3.0 * g_yyy_zz_1[i] * fe_0 + g_yyy_zzz_1[i] * pa_z[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto g_yyzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 120);

    auto g_yyzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 121);

    auto g_yyzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 122);

    auto g_yyzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 123);

    auto g_yyzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 124);

    auto g_yyzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 125);

    auto g_yyzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 126);

    auto g_yyzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 127);

    auto g_yyzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 128);

    auto g_yyzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 129);

    #pragma omp simd aligned(g_yy_xxy_0, g_yy_xxy_1, g_yy_xyy_0, g_yy_xyy_1, g_yy_yyy_0, g_yy_yyy_1, g_yyz_xxy_1, g_yyz_xyy_1, g_yyz_yyy_1, g_yyzz_xxx_0, g_yyzz_xxy_0, g_yyzz_xxz_0, g_yyzz_xyy_0, g_yyzz_xyz_0, g_yyzz_xzz_0, g_yyzz_yyy_0, g_yyzz_yyz_0, g_yyzz_yzz_0, g_yyzz_zzz_0, g_yzz_xxx_1, g_yzz_xxz_1, g_yzz_xyz_1, g_yzz_xz_1, g_yzz_xzz_1, g_yzz_yyz_1, g_yzz_yz_1, g_yzz_yzz_1, g_yzz_zz_1, g_yzz_zzz_1, g_zz_xxx_0, g_zz_xxx_1, g_zz_xxz_0, g_zz_xxz_1, g_zz_xyz_0, g_zz_xyz_1, g_zz_xzz_0, g_zz_xzz_1, g_zz_yyz_0, g_zz_yyz_1, g_zz_yzz_0, g_zz_yzz_1, g_zz_zzz_0, g_zz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzz_xxx_0[i] = g_zz_xxx_0[i] * fbe_0 - g_zz_xxx_1[i] * fz_be_0 + g_yzz_xxx_1[i] * pa_y[i];

        g_yyzz_xxy_0[i] = g_yy_xxy_0[i] * fbe_0 - g_yy_xxy_1[i] * fz_be_0 + g_yyz_xxy_1[i] * pa_z[i];

        g_yyzz_xxz_0[i] = g_zz_xxz_0[i] * fbe_0 - g_zz_xxz_1[i] * fz_be_0 + g_yzz_xxz_1[i] * pa_y[i];

        g_yyzz_xyy_0[i] = g_yy_xyy_0[i] * fbe_0 - g_yy_xyy_1[i] * fz_be_0 + g_yyz_xyy_1[i] * pa_z[i];

        g_yyzz_xyz_0[i] = g_zz_xyz_0[i] * fbe_0 - g_zz_xyz_1[i] * fz_be_0 + g_yzz_xz_1[i] * fe_0 + g_yzz_xyz_1[i] * pa_y[i];

        g_yyzz_xzz_0[i] = g_zz_xzz_0[i] * fbe_0 - g_zz_xzz_1[i] * fz_be_0 + g_yzz_xzz_1[i] * pa_y[i];

        g_yyzz_yyy_0[i] = g_yy_yyy_0[i] * fbe_0 - g_yy_yyy_1[i] * fz_be_0 + g_yyz_yyy_1[i] * pa_z[i];

        g_yyzz_yyz_0[i] = g_zz_yyz_0[i] * fbe_0 - g_zz_yyz_1[i] * fz_be_0 + 2.0 * g_yzz_yz_1[i] * fe_0 + g_yzz_yyz_1[i] * pa_y[i];

        g_yyzz_yzz_0[i] = g_zz_yzz_0[i] * fbe_0 - g_zz_yzz_1[i] * fz_be_0 + g_yzz_zz_1[i] * fe_0 + g_yzz_yzz_1[i] * pa_y[i];

        g_yyzz_zzz_0[i] = g_zz_zzz_0[i] * fbe_0 - g_zz_zzz_1[i] * fz_be_0 + g_yzz_zzz_1[i] * pa_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto g_yzzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 130);

    auto g_yzzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 131);

    auto g_yzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 132);

    auto g_yzzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 133);

    auto g_yzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 134);

    auto g_yzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 135);

    auto g_yzzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 136);

    auto g_yzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 137);

    auto g_yzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 138);

    auto g_yzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 139);

    #pragma omp simd aligned(g_yzzz_xxx_0, g_yzzz_xxy_0, g_yzzz_xxz_0, g_yzzz_xyy_0, g_yzzz_xyz_0, g_yzzz_xzz_0, g_yzzz_yyy_0, g_yzzz_yyz_0, g_yzzz_yzz_0, g_yzzz_zzz_0, g_zzz_xx_1, g_zzz_xxx_1, g_zzz_xxy_1, g_zzz_xxz_1, g_zzz_xy_1, g_zzz_xyy_1, g_zzz_xyz_1, g_zzz_xz_1, g_zzz_xzz_1, g_zzz_yy_1, g_zzz_yyy_1, g_zzz_yyz_1, g_zzz_yz_1, g_zzz_yzz_1, g_zzz_zz_1, g_zzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_xxx_0[i] = g_zzz_xxx_1[i] * pa_y[i];

        g_yzzz_xxy_0[i] = g_zzz_xx_1[i] * fe_0 + g_zzz_xxy_1[i] * pa_y[i];

        g_yzzz_xxz_0[i] = g_zzz_xxz_1[i] * pa_y[i];

        g_yzzz_xyy_0[i] = 2.0 * g_zzz_xy_1[i] * fe_0 + g_zzz_xyy_1[i] * pa_y[i];

        g_yzzz_xyz_0[i] = g_zzz_xz_1[i] * fe_0 + g_zzz_xyz_1[i] * pa_y[i];

        g_yzzz_xzz_0[i] = g_zzz_xzz_1[i] * pa_y[i];

        g_yzzz_yyy_0[i] = 3.0 * g_zzz_yy_1[i] * fe_0 + g_zzz_yyy_1[i] * pa_y[i];

        g_yzzz_yyz_0[i] = 2.0 * g_zzz_yz_1[i] * fe_0 + g_zzz_yyz_1[i] * pa_y[i];

        g_yzzz_yzz_0[i] = g_zzz_zz_1[i] * fe_0 + g_zzz_yzz_1[i] * pa_y[i];

        g_yzzz_zzz_0[i] = g_zzz_zzz_1[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto g_zzzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 140);

    auto g_zzzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 141);

    auto g_zzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 142);

    auto g_zzzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 143);

    auto g_zzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 144);

    auto g_zzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 145);

    auto g_zzzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 146);

    auto g_zzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 147);

    auto g_zzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 148);

    auto g_zzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 149);

    #pragma omp simd aligned(g_zz_xxx_0, g_zz_xxx_1, g_zz_xxy_0, g_zz_xxy_1, g_zz_xxz_0, g_zz_xxz_1, g_zz_xyy_0, g_zz_xyy_1, g_zz_xyz_0, g_zz_xyz_1, g_zz_xzz_0, g_zz_xzz_1, g_zz_yyy_0, g_zz_yyy_1, g_zz_yyz_0, g_zz_yyz_1, g_zz_yzz_0, g_zz_yzz_1, g_zz_zzz_0, g_zz_zzz_1, g_zzz_xx_1, g_zzz_xxx_1, g_zzz_xxy_1, g_zzz_xxz_1, g_zzz_xy_1, g_zzz_xyy_1, g_zzz_xyz_1, g_zzz_xz_1, g_zzz_xzz_1, g_zzz_yy_1, g_zzz_yyy_1, g_zzz_yyz_1, g_zzz_yz_1, g_zzz_yzz_1, g_zzz_zz_1, g_zzz_zzz_1, g_zzzz_xxx_0, g_zzzz_xxy_0, g_zzzz_xxz_0, g_zzzz_xyy_0, g_zzzz_xyz_0, g_zzzz_xzz_0, g_zzzz_yyy_0, g_zzzz_yyz_0, g_zzzz_yzz_0, g_zzzz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_xxx_0[i] = 3.0 * g_zz_xxx_0[i] * fbe_0 - 3.0 * g_zz_xxx_1[i] * fz_be_0 + g_zzz_xxx_1[i] * pa_z[i];

        g_zzzz_xxy_0[i] = 3.0 * g_zz_xxy_0[i] * fbe_0 - 3.0 * g_zz_xxy_1[i] * fz_be_0 + g_zzz_xxy_1[i] * pa_z[i];

        g_zzzz_xxz_0[i] = 3.0 * g_zz_xxz_0[i] * fbe_0 - 3.0 * g_zz_xxz_1[i] * fz_be_0 + g_zzz_xx_1[i] * fe_0 + g_zzz_xxz_1[i] * pa_z[i];

        g_zzzz_xyy_0[i] = 3.0 * g_zz_xyy_0[i] * fbe_0 - 3.0 * g_zz_xyy_1[i] * fz_be_0 + g_zzz_xyy_1[i] * pa_z[i];

        g_zzzz_xyz_0[i] = 3.0 * g_zz_xyz_0[i] * fbe_0 - 3.0 * g_zz_xyz_1[i] * fz_be_0 + g_zzz_xy_1[i] * fe_0 + g_zzz_xyz_1[i] * pa_z[i];

        g_zzzz_xzz_0[i] = 3.0 * g_zz_xzz_0[i] * fbe_0 - 3.0 * g_zz_xzz_1[i] * fz_be_0 + 2.0 * g_zzz_xz_1[i] * fe_0 + g_zzz_xzz_1[i] * pa_z[i];

        g_zzzz_yyy_0[i] = 3.0 * g_zz_yyy_0[i] * fbe_0 - 3.0 * g_zz_yyy_1[i] * fz_be_0 + g_zzz_yyy_1[i] * pa_z[i];

        g_zzzz_yyz_0[i] = 3.0 * g_zz_yyz_0[i] * fbe_0 - 3.0 * g_zz_yyz_1[i] * fz_be_0 + g_zzz_yy_1[i] * fe_0 + g_zzz_yyz_1[i] * pa_z[i];

        g_zzzz_yzz_0[i] = 3.0 * g_zz_yzz_0[i] * fbe_0 - 3.0 * g_zz_yzz_1[i] * fz_be_0 + 2.0 * g_zzz_yz_1[i] * fe_0 + g_zzz_yzz_1[i] * pa_z[i];

        g_zzzz_zzz_0[i] = 3.0 * g_zz_zzz_0[i] * fbe_0 - 3.0 * g_zz_zzz_1[i] * fz_be_0 + 3.0 * g_zzz_zz_1[i] * fe_0 + g_zzz_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

