#include "ProjectedCorePotentialPrimRecSGForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_sg_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_sg_s_0_0_0,
                                        const size_t idx_sd_s_0_0_0,
                                        const size_t idx_sf_s_0_0_0,
                                        const size_t idx_sd_s_0_1_0,
                                        const size_t idx_sf_s_0_1_0,
                                        const int m,
                                        const size_t idx_sd_s_0_0_1,
                                        const size_t idx_sf_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up B center coordinates

    auto rb_x = factors.data(idx_b);

    auto rb_y = factors.data(idx_b + 1);

    auto rb_z = factors.data(idx_b + 2);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0);

    auto tg_0_xy_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 1);

    auto tg_0_xz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 2);

    auto tg_0_yy_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 3);

    auto tg_0_yz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 4);

    auto tg_0_zz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 5);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0);

    auto tg_0_xxy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 1);

    auto tg_0_xxz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 2);

    auto tg_0_xyy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 3);

    auto tg_0_xyz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 4);

    auto tg_0_xzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 5);

    auto tg_0_yyy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 6);

    auto tg_0_yyz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 7);

    auto tg_0_yzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 8);

    auto tg_0_zzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 9);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0);

    auto tg_0_xy_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0 + 1);

    auto tg_0_xz_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0 + 2);

    auto tg_0_yy_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0 + 3);

    auto tg_0_yz_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0 + 4);

    auto tg_0_zz_s_0_1_0 = pbuffer.data(idx_sd_s_0_1_0 + 5);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0);

    auto tg_0_xxy_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 1);

    auto tg_0_xxz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 2);

    auto tg_0_xyy_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 3);

    auto tg_0_xyz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 4);

    auto tg_0_xzz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 5);

    auto tg_0_yyy_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 6);

    auto tg_0_yyz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 7);

    auto tg_0_yzz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 8);

    auto tg_0_zzz_s_0_1_0 = pbuffer.data(idx_sf_s_0_1_0 + 9);

    // Set up components of targeted buffer : SG

    auto tg_0_xxxx_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0);

    auto tg_0_xxxy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 1);

    auto tg_0_xxxz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 2);

    auto tg_0_xxyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 3);

    auto tg_0_xxyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 4);

    auto tg_0_xxzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 5);

    auto tg_0_xyyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 6);

    auto tg_0_xyyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 7);

    auto tg_0_xyzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 8);

    auto tg_0_xzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 9);

    auto tg_0_yyyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 10);

    auto tg_0_yyyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 11);

    auto tg_0_yyzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 12);

    auto tg_0_yzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 13);

    auto tg_0_zzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 14);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xx_s_0_0_0, tg_0_xx_s_0_1_0, tg_0_xxx_s_0_0_0, tg_0_xxx_s_0_1_0, tg_0_xxxx_s_0_0_0, tg_0_xxxy_s_0_0_0, tg_0_xxxz_s_0_0_0, tg_0_xxyy_s_0_0_0, tg_0_xxyz_s_0_0_0, tg_0_xxz_s_0_0_0, tg_0_xxz_s_0_1_0, tg_0_xxzz_s_0_0_0, tg_0_xyy_s_0_0_0, tg_0_xyy_s_0_1_0, tg_0_xyyy_s_0_0_0, tg_0_xyyz_s_0_0_0, tg_0_xyzz_s_0_0_0, tg_0_xzz_s_0_0_0, tg_0_xzz_s_0_1_0, tg_0_xzzz_s_0_0_0, tg_0_yy_s_0_0_0, tg_0_yy_s_0_1_0, tg_0_yyy_s_0_0_0, tg_0_yyy_s_0_1_0, tg_0_yyyy_s_0_0_0, tg_0_yyyz_s_0_0_0, tg_0_yyz_s_0_0_0, tg_0_yyz_s_0_1_0, tg_0_yyzz_s_0_0_0, tg_0_yzz_s_0_0_0, tg_0_yzz_s_0_1_0, tg_0_yzzz_s_0_0_0, tg_0_zz_s_0_0_0, tg_0_zz_s_0_1_0, tg_0_zzz_s_0_0_0, tg_0_zzz_s_0_1_0, tg_0_zzzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double fbz_0 = -(a_exp + c_exp) * fzi_0;

        const double fazi_0 = a_exp * fzi_0;

        const double fb_0 = b_exps[i];

        tg_0_xxxx_s_0_0_0[i] = 3.0 * tg_0_xx_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_xx_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxx_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxx_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxy_s_0_0_0[i] = 2.0 * tg_0_xxx_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxx_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxz_s_0_0_0[i] = 2.0 * tg_0_xxx_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxx_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxyy_s_0_0_0[i] = tg_0_yy_s_0_0_0[i] * fzi_0 + tg_0_yy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyz_s_0_0_0[i] = 2.0 * tg_0_xxz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxzz_s_0_0_0[i] = tg_0_zz_s_0_0_0[i] * fzi_0 + tg_0_zz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyy_s_0_0_0[i] = 2.0 * tg_0_yyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyz_s_0_0_0[i] = 2.0 * tg_0_yyz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyzz_s_0_0_0[i] = 2.0 * tg_0_yzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xzzz_s_0_0_0[i] = 2.0 * tg_0_zzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_yyyy_s_0_0_0[i] = 3.0 * tg_0_yy_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_yy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yyy_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyy_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyyz_s_0_0_0[i] = 2.0 * tg_0_yyy_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyy_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_yyzz_s_0_0_0[i] = tg_0_zz_s_0_0_0[i] * fzi_0 + tg_0_zz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yzzz_s_0_0_0[i] = 2.0 * tg_0_zzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_zzzz_s_0_0_0[i] = 3.0 * tg_0_zz_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_zz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_zzz_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzz_s_0_0_0[i] * rb_z[i] * fbz_0;
    }

    if (m > 0)
    {
        const double fm_0 = (double)m;

        // Set up components of auxiliary buffer : SD

        auto tg_0_xx_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1);

        auto tg_0_xy_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 1);

        auto tg_0_xz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 2);

        auto tg_0_yy_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 3);

        auto tg_0_yz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 4);

        auto tg_0_zz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 5);

        // Set up components of auxiliary buffer : SF

        auto tg_0_xxx_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1);

        auto tg_0_xxy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 1);

        auto tg_0_xxz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 2);

        auto tg_0_xyy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 3);

        auto tg_0_xyz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 4);

        auto tg_0_xzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 5);

        auto tg_0_yyy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 6);

        auto tg_0_yyz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 7);

        auto tg_0_yzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 8);

        auto tg_0_zzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 9);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xx_s_0_0_1, tg_0_xxx_s_0_0_1, tg_0_xxxx_s_0_0_0, tg_0_xxxy_s_0_0_0, tg_0_xxxz_s_0_0_0, tg_0_xxyy_s_0_0_0, tg_0_xxyz_s_0_0_0, tg_0_xxz_s_0_0_1, tg_0_xxzz_s_0_0_0, tg_0_xyy_s_0_0_1, tg_0_xyyy_s_0_0_0, tg_0_xyyz_s_0_0_0, tg_0_xyzz_s_0_0_0, tg_0_xzz_s_0_0_1, tg_0_xzzz_s_0_0_0, tg_0_yy_s_0_0_1, tg_0_yyy_s_0_0_1, tg_0_yyyy_s_0_0_0, tg_0_yyyz_s_0_0_0, tg_0_yyz_s_0_0_1, tg_0_yyzz_s_0_0_0, tg_0_yzz_s_0_0_1, tg_0_yzzz_s_0_0_0, tg_0_zz_s_0_0_1, tg_0_zzz_s_0_0_1, tg_0_zzzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fbi_0 = 1.0 / b_exps[i];

            tg_0_xxxx_s_0_0_0[i] = 3.0 / 2.0 * tg_0_xx_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxx_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxy_s_0_0_0[i] = tg_0_xxx_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxz_s_0_0_0[i] = tg_0_xxx_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyz_s_0_0_0[i] = tg_0_xxz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyy_s_0_0_0[i] = tg_0_yyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyz_s_0_0_0[i] = tg_0_yyz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyzz_s_0_0_0[i] = tg_0_yzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xzzz_s_0_0_0[i] = tg_0_zzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_yyyy_s_0_0_0[i] = 3.0 / 2.0 * tg_0_yy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyy_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyyz_s_0_0_0[i] = tg_0_yyy_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_yyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yzzz_s_0_0_0[i] = tg_0_zzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_zzzz_s_0_0_0[i] = 3.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_zzz_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;
        }
    }
}

} // t2pecp namespace

