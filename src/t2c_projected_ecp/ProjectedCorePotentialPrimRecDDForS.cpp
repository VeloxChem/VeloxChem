#include "ProjectedCorePotentialPrimRecDDForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_dd_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_dd_s_0_0_0,
                                        const size_t idx_sd_s_0_0_0,
                                        const size_t idx_pd_s_0_0_0,
                                        const size_t idx_sd_s_1_0_0,
                                        const size_t idx_pd_s_1_0_0,
                                        const int p,
                                        const size_t idx_sd_s_0_0_1,
                                        const size_t idx_pd_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0);

    auto tg_0_xy_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 1);

    auto tg_0_xz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 2);

    auto tg_0_yy_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 3);

    auto tg_0_yz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 4);

    auto tg_0_zz_s_0_0_0 = pbuffer.data(idx_sd_s_0_0_0 + 5);

    // Set up components of auxiliary buffer : PD

    auto tg_x_xx_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0);

    auto tg_x_xy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 1);

    auto tg_x_xz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 2);

    auto tg_x_yy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 3);

    auto tg_x_yz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 4);

    auto tg_x_zz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 5);

    auto tg_y_xx_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 6);

    auto tg_y_xy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 7);

    auto tg_y_xz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 8);

    auto tg_y_yy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 9);

    auto tg_y_yz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 10);

    auto tg_y_zz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 11);

    auto tg_z_xx_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 12);

    auto tg_z_xy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 13);

    auto tg_z_xz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 14);

    auto tg_z_yy_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 15);

    auto tg_z_yz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 16);

    auto tg_z_zz_s_0_0_0 = pbuffer.data(idx_pd_s_0_0_0 + 17);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0);

    auto tg_0_xy_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0 + 1);

    auto tg_0_xz_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0 + 2);

    auto tg_0_yy_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0 + 3);

    auto tg_0_yz_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0 + 4);

    auto tg_0_zz_s_1_0_0 = pbuffer.data(idx_sd_s_1_0_0 + 5);

    // Set up components of auxiliary buffer : PD

    auto tg_x_xx_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0);

    auto tg_x_xy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 1);

    auto tg_x_xz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 2);

    auto tg_x_yy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 3);

    auto tg_x_yz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 4);

    auto tg_x_zz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 5);

    auto tg_y_xx_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 6);

    auto tg_y_xy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 7);

    auto tg_y_xz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 8);

    auto tg_y_yy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 9);

    auto tg_y_yz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 10);

    auto tg_y_zz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 11);

    auto tg_z_xx_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 12);

    auto tg_z_xy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 13);

    auto tg_z_xz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 14);

    auto tg_z_yy_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 15);

    auto tg_z_yz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 16);

    auto tg_z_zz_s_1_0_0 = pbuffer.data(idx_pd_s_1_0_0 + 17);

    // Set up components of targeted buffer : DD

    auto tg_xx_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0);

    auto tg_xx_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 1);

    auto tg_xx_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 2);

    auto tg_xx_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 3);

    auto tg_xx_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 4);

    auto tg_xx_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 5);

    auto tg_xy_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 6);

    auto tg_xy_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 7);

    auto tg_xy_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 8);

    auto tg_xy_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 9);

    auto tg_xy_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 10);

    auto tg_xy_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 11);

    auto tg_xz_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 12);

    auto tg_xz_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 13);

    auto tg_xz_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 14);

    auto tg_xz_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 15);

    auto tg_xz_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 16);

    auto tg_xz_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 17);

    auto tg_yy_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 18);

    auto tg_yy_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 19);

    auto tg_yy_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 20);

    auto tg_yy_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 21);

    auto tg_yy_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 22);

    auto tg_yy_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 23);

    auto tg_yz_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 24);

    auto tg_yz_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 25);

    auto tg_yz_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 26);

    auto tg_yz_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 27);

    auto tg_yz_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 28);

    auto tg_yz_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 29);

    auto tg_zz_xx_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 30);

    auto tg_zz_xy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 31);

    auto tg_zz_xz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 32);

    auto tg_zz_yy_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 33);

    auto tg_zz_yz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 34);

    auto tg_zz_zz_s_0_0_0 = pbuffer.data(idx_dd_s_0_0_0 + 35);

    #pragma omp simd aligned(b_exps, tg_0_xx_s_0_0_0, tg_0_xx_s_1_0_0, tg_0_xy_s_0_0_0, tg_0_xy_s_1_0_0, tg_0_xz_s_0_0_0, tg_0_xz_s_1_0_0, tg_0_yy_s_0_0_0, tg_0_yy_s_1_0_0, tg_0_yz_s_0_0_0, tg_0_yz_s_1_0_0, tg_0_zz_s_0_0_0, tg_0_zz_s_1_0_0, tg_x_xx_s_0_0_0, tg_x_xx_s_1_0_0, tg_x_xy_s_0_0_0, tg_x_xy_s_1_0_0, tg_x_xz_s_0_0_0, tg_x_xz_s_1_0_0, tg_x_yy_s_0_0_0, tg_x_yy_s_1_0_0, tg_x_yz_s_0_0_0, tg_x_yz_s_1_0_0, tg_x_zz_s_0_0_0, tg_x_zz_s_1_0_0, tg_xx_xx_s_0_0_0, tg_xx_xy_s_0_0_0, tg_xx_xz_s_0_0_0, tg_xx_yy_s_0_0_0, tg_xx_yz_s_0_0_0, tg_xx_zz_s_0_0_0, tg_xy_xx_s_0_0_0, tg_xy_xy_s_0_0_0, tg_xy_xz_s_0_0_0, tg_xy_yy_s_0_0_0, tg_xy_yz_s_0_0_0, tg_xy_zz_s_0_0_0, tg_xz_xx_s_0_0_0, tg_xz_xy_s_0_0_0, tg_xz_xz_s_0_0_0, tg_xz_yy_s_0_0_0, tg_xz_yz_s_0_0_0, tg_xz_zz_s_0_0_0, tg_y_xx_s_0_0_0, tg_y_xx_s_1_0_0, tg_y_xy_s_0_0_0, tg_y_xy_s_1_0_0, tg_y_xz_s_0_0_0, tg_y_xz_s_1_0_0, tg_y_yy_s_0_0_0, tg_y_yy_s_1_0_0, tg_y_yz_s_0_0_0, tg_y_yz_s_1_0_0, tg_y_zz_s_0_0_0, tg_y_zz_s_1_0_0, tg_yy_xx_s_0_0_0, tg_yy_xy_s_0_0_0, tg_yy_xz_s_0_0_0, tg_yy_yy_s_0_0_0, tg_yy_yz_s_0_0_0, tg_yy_zz_s_0_0_0, tg_yz_xx_s_0_0_0, tg_yz_xy_s_0_0_0, tg_yz_xz_s_0_0_0, tg_yz_yy_s_0_0_0, tg_yz_yz_s_0_0_0, tg_yz_zz_s_0_0_0, tg_z_xx_s_0_0_0, tg_z_xx_s_1_0_0, tg_z_xy_s_0_0_0, tg_z_xy_s_1_0_0, tg_z_xz_s_0_0_0, tg_z_xz_s_1_0_0, tg_z_yy_s_0_0_0, tg_z_yy_s_1_0_0, tg_z_yz_s_0_0_0, tg_z_yz_s_1_0_0, tg_z_zz_s_0_0_0, tg_z_zz_s_1_0_0, tg_zz_xx_s_0_0_0, tg_zz_xy_s_0_0_0, tg_zz_xz_s_0_0_0, tg_zz_yy_s_0_0_0, tg_zz_yz_s_0_0_0, tg_zz_zz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xx_xx_s_0_0_0[i] = tg_0_xx_s_0_0_0[i] * fzi_0 + tg_0_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xy_s_0_0_0[i] = tg_0_xy_s_0_0_0[i] * fzi_0 + tg_0_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xz_s_0_0_0[i] = tg_0_xz_s_0_0_0[i] * fzi_0 + tg_0_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yy_s_0_0_0[i] = tg_0_yy_s_0_0_0[i] * fzi_0 + tg_0_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yz_s_0_0_0[i] = tg_0_yz_s_0_0_0[i] * fzi_0 + tg_0_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_zz_s_0_0_0[i] = tg_0_zz_s_0_0_0[i] * fzi_0 + tg_0_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xx_s_0_0_0[i] = 2.0 * tg_y_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xy_s_0_0_0[i] = 2.0 * tg_y_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xz_s_0_0_0[i] = 2.0 * tg_y_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yy_s_0_0_0[i] = 2.0 * tg_y_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yz_s_0_0_0[i] = 2.0 * tg_y_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_zz_s_0_0_0[i] = 2.0 * tg_y_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xx_s_0_0_0[i] = 2.0 * tg_z_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xy_s_0_0_0[i] = 2.0 * tg_z_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xz_s_0_0_0[i] = 2.0 * tg_z_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yy_s_0_0_0[i] = 2.0 * tg_z_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yz_s_0_0_0[i] = 2.0 * tg_z_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_zz_s_0_0_0[i] = 2.0 * tg_z_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_zz_s_0_0_0[i] * a_x * faz_0;

        tg_yy_xx_s_0_0_0[i] = tg_0_xx_s_0_0_0[i] * fzi_0 + tg_0_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xy_s_0_0_0[i] = tg_0_xy_s_0_0_0[i] * fzi_0 + tg_0_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xz_s_0_0_0[i] = tg_0_xz_s_0_0_0[i] * fzi_0 + tg_0_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yy_s_0_0_0[i] = tg_0_yy_s_0_0_0[i] * fzi_0 + tg_0_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yz_s_0_0_0[i] = tg_0_yz_s_0_0_0[i] * fzi_0 + tg_0_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_zz_s_0_0_0[i] = tg_0_zz_s_0_0_0[i] * fzi_0 + tg_0_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_zz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xx_s_0_0_0[i] = 2.0 * tg_z_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xy_s_0_0_0[i] = 2.0 * tg_z_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xz_s_0_0_0[i] = 2.0 * tg_z_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yy_s_0_0_0[i] = 2.0 * tg_z_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yz_s_0_0_0[i] = 2.0 * tg_z_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_zz_s_0_0_0[i] = 2.0 * tg_z_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_zz_s_0_0_0[i] * a_y * faz_0;

        tg_zz_xx_s_0_0_0[i] = tg_0_xx_s_0_0_0[i] * fzi_0 + tg_0_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xx_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xy_s_0_0_0[i] = tg_0_xy_s_0_0_0[i] * fzi_0 + tg_0_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xz_s_0_0_0[i] = tg_0_xz_s_0_0_0[i] * fzi_0 + tg_0_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yy_s_0_0_0[i] = tg_0_yy_s_0_0_0[i] * fzi_0 + tg_0_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yz_s_0_0_0[i] = tg_0_yz_s_0_0_0[i] * fzi_0 + tg_0_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_zz_s_0_0_0[i] = tg_0_zz_s_0_0_0[i] * fzi_0 + tg_0_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_zz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SD

        auto tg_0_xx_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1);

        auto tg_0_xy_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 1);

        auto tg_0_xz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 2);

        auto tg_0_yy_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 3);

        auto tg_0_yz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 4);

        auto tg_0_zz_s_0_0_1 = pbuffer.data(idx_sd_s_0_0_1 + 5);

        // Set up components of auxiliary buffer : PD

        auto tg_x_xx_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1);

        auto tg_x_xy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 1);

        auto tg_x_xz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 2);

        auto tg_x_yy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 3);

        auto tg_x_yz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 4);

        auto tg_x_zz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 5);

        auto tg_y_xx_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 6);

        auto tg_y_xy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 7);

        auto tg_y_xz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 8);

        auto tg_y_yy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 9);

        auto tg_y_yz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 10);

        auto tg_y_zz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 11);

        auto tg_z_xx_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 12);

        auto tg_z_xy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 13);

        auto tg_z_xz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 14);

        auto tg_z_yy_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 15);

        auto tg_z_yz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 16);

        auto tg_z_zz_s_0_0_1 = pbuffer.data(idx_pd_s_0_0_1 + 17);

        #pragma omp simd aligned(b_exps, tg_0_xx_s_0_0_1, tg_0_xy_s_0_0_1, tg_0_xz_s_0_0_1, tg_0_yy_s_0_0_1, tg_0_yz_s_0_0_1, tg_0_zz_s_0_0_1, tg_x_xx_s_0_0_1, tg_x_xy_s_0_0_1, tg_x_xz_s_0_0_1, tg_x_yy_s_0_0_1, tg_x_yz_s_0_0_1, tg_x_zz_s_0_0_1, tg_xx_xx_s_0_0_0, tg_xx_xy_s_0_0_0, tg_xx_xz_s_0_0_0, tg_xx_yy_s_0_0_0, tg_xx_yz_s_0_0_0, tg_xx_zz_s_0_0_0, tg_xy_xx_s_0_0_0, tg_xy_xy_s_0_0_0, tg_xy_xz_s_0_0_0, tg_xy_yy_s_0_0_0, tg_xy_yz_s_0_0_0, tg_xy_zz_s_0_0_0, tg_xz_xx_s_0_0_0, tg_xz_xy_s_0_0_0, tg_xz_xz_s_0_0_0, tg_xz_yy_s_0_0_0, tg_xz_yz_s_0_0_0, tg_xz_zz_s_0_0_0, tg_y_xx_s_0_0_1, tg_y_xy_s_0_0_1, tg_y_xz_s_0_0_1, tg_y_yy_s_0_0_1, tg_y_yz_s_0_0_1, tg_y_zz_s_0_0_1, tg_yy_xx_s_0_0_0, tg_yy_xy_s_0_0_0, tg_yy_xz_s_0_0_0, tg_yy_yy_s_0_0_0, tg_yy_yz_s_0_0_0, tg_yy_zz_s_0_0_0, tg_yz_xx_s_0_0_0, tg_yz_xy_s_0_0_0, tg_yz_xz_s_0_0_0, tg_yz_yy_s_0_0_0, tg_yz_yz_s_0_0_0, tg_yz_zz_s_0_0_0, tg_z_xx_s_0_0_1, tg_z_xy_s_0_0_1, tg_z_xz_s_0_0_1, tg_z_yy_s_0_0_1, tg_z_yz_s_0_0_1, tg_z_zz_s_0_0_1, tg_zz_xx_s_0_0_0, tg_zz_xy_s_0_0_0, tg_zz_xz_s_0_0_0, tg_zz_yy_s_0_0_0, tg_zz_yz_s_0_0_0, tg_zz_zz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xx_s_0_0_0[i] = tg_y_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xy_s_0_0_0[i] = tg_y_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xz_s_0_0_0[i] = tg_y_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yy_s_0_0_0[i] = tg_y_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yz_s_0_0_0[i] = tg_y_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_zz_s_0_0_0[i] = tg_y_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xx_s_0_0_0[i] = tg_z_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xy_s_0_0_0[i] = tg_z_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xz_s_0_0_0[i] = tg_z_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yy_s_0_0_0[i] = tg_z_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yz_s_0_0_0[i] = tg_z_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_zz_s_0_0_0[i] = tg_z_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xx_s_0_0_0[i] = tg_z_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xy_s_0_0_0[i] = tg_z_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xz_s_0_0_0[i] = tg_z_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yy_s_0_0_0[i] = tg_z_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yz_s_0_0_0[i] = tg_z_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_zz_s_0_0_0[i] = tg_z_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

