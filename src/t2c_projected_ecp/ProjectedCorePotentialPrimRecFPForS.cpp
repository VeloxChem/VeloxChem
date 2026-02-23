#include "ProjectedCorePotentialPrimRecFPForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_fp_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_fp_s_0_0_0,
                                        const size_t idx_pp_s_0_0_0,
                                        const size_t idx_dp_s_0_0_0,
                                        const size_t idx_pp_s_1_0_0,
                                        const size_t idx_dp_s_1_0_0,
                                        const int p,
                                        const size_t idx_pp_s_0_0_1,
                                        const size_t idx_dp_s_0_0_1,
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

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0);

    auto tg_x_y_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 1);

    auto tg_x_z_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 2);

    auto tg_y_x_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 3);

    auto tg_y_y_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 4);

    auto tg_y_z_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 5);

    auto tg_z_x_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 6);

    auto tg_z_y_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 7);

    auto tg_z_z_s_0_0_0 = pbuffer.data(idx_pp_s_0_0_0 + 8);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0);

    auto tg_xx_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 1);

    auto tg_xx_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 2);

    auto tg_xy_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 3);

    auto tg_xy_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 4);

    auto tg_xy_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 5);

    auto tg_xz_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 6);

    auto tg_xz_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 7);

    auto tg_xz_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 8);

    auto tg_yy_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 9);

    auto tg_yy_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 10);

    auto tg_yy_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 11);

    auto tg_yz_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 12);

    auto tg_yz_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 13);

    auto tg_yz_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 14);

    auto tg_zz_x_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 15);

    auto tg_zz_y_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 16);

    auto tg_zz_z_s_0_0_0 = pbuffer.data(idx_dp_s_0_0_0 + 17);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0);

    auto tg_x_y_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 1);

    auto tg_x_z_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 2);

    auto tg_y_x_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 3);

    auto tg_y_y_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 4);

    auto tg_y_z_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 5);

    auto tg_z_x_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 6);

    auto tg_z_y_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 7);

    auto tg_z_z_s_1_0_0 = pbuffer.data(idx_pp_s_1_0_0 + 8);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0);

    auto tg_xx_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 1);

    auto tg_xx_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 2);

    auto tg_xy_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 3);

    auto tg_xy_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 4);

    auto tg_xy_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 5);

    auto tg_xz_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 6);

    auto tg_xz_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 7);

    auto tg_xz_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 8);

    auto tg_yy_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 9);

    auto tg_yy_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 10);

    auto tg_yy_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 11);

    auto tg_yz_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 12);

    auto tg_yz_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 13);

    auto tg_yz_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 14);

    auto tg_zz_x_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 15);

    auto tg_zz_y_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 16);

    auto tg_zz_z_s_1_0_0 = pbuffer.data(idx_dp_s_1_0_0 + 17);

    // Set up components of targeted buffer : FP

    auto tg_xxx_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0);

    auto tg_xxx_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 1);

    auto tg_xxx_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 2);

    auto tg_xxy_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 3);

    auto tg_xxy_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 4);

    auto tg_xxy_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 5);

    auto tg_xxz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 6);

    auto tg_xxz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 7);

    auto tg_xxz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 8);

    auto tg_xyy_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 9);

    auto tg_xyy_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 10);

    auto tg_xyy_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 11);

    auto tg_xyz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 12);

    auto tg_xyz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 13);

    auto tg_xyz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 14);

    auto tg_xzz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 15);

    auto tg_xzz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 16);

    auto tg_xzz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 17);

    auto tg_yyy_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 18);

    auto tg_yyy_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 19);

    auto tg_yyy_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 20);

    auto tg_yyz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 21);

    auto tg_yyz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 22);

    auto tg_yyz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 23);

    auto tg_yzz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 24);

    auto tg_yzz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 25);

    auto tg_yzz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 26);

    auto tg_zzz_x_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 27);

    auto tg_zzz_y_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 28);

    auto tg_zzz_z_s_0_0_0 = pbuffer.data(idx_fp_s_0_0_0 + 29);

    #pragma omp simd aligned(b_exps, tg_x_x_s_0_0_0, tg_x_x_s_1_0_0, tg_x_y_s_0_0_0, tg_x_y_s_1_0_0, tg_x_z_s_0_0_0, tg_x_z_s_1_0_0, tg_xx_x_s_0_0_0, tg_xx_x_s_1_0_0, tg_xx_y_s_0_0_0, tg_xx_y_s_1_0_0, tg_xx_z_s_0_0_0, tg_xx_z_s_1_0_0, tg_xxx_x_s_0_0_0, tg_xxx_y_s_0_0_0, tg_xxx_z_s_0_0_0, tg_xxy_x_s_0_0_0, tg_xxy_y_s_0_0_0, tg_xxy_z_s_0_0_0, tg_xxz_x_s_0_0_0, tg_xxz_y_s_0_0_0, tg_xxz_z_s_0_0_0, tg_xyy_x_s_0_0_0, tg_xyy_y_s_0_0_0, tg_xyy_z_s_0_0_0, tg_xyz_x_s_0_0_0, tg_xyz_y_s_0_0_0, tg_xyz_z_s_0_0_0, tg_xzz_x_s_0_0_0, tg_xzz_y_s_0_0_0, tg_xzz_z_s_0_0_0, tg_y_x_s_0_0_0, tg_y_x_s_1_0_0, tg_y_y_s_0_0_0, tg_y_y_s_1_0_0, tg_y_z_s_0_0_0, tg_y_z_s_1_0_0, tg_yy_x_s_0_0_0, tg_yy_x_s_1_0_0, tg_yy_y_s_0_0_0, tg_yy_y_s_1_0_0, tg_yy_z_s_0_0_0, tg_yy_z_s_1_0_0, tg_yyy_x_s_0_0_0, tg_yyy_y_s_0_0_0, tg_yyy_z_s_0_0_0, tg_yyz_x_s_0_0_0, tg_yyz_y_s_0_0_0, tg_yyz_z_s_0_0_0, tg_yz_x_s_0_0_0, tg_yz_x_s_1_0_0, tg_yz_y_s_0_0_0, tg_yz_y_s_1_0_0, tg_yz_z_s_0_0_0, tg_yz_z_s_1_0_0, tg_yzz_x_s_0_0_0, tg_yzz_y_s_0_0_0, tg_yzz_z_s_0_0_0, tg_z_x_s_0_0_0, tg_z_x_s_1_0_0, tg_z_y_s_0_0_0, tg_z_y_s_1_0_0, tg_z_z_s_0_0_0, tg_z_z_s_1_0_0, tg_zz_x_s_0_0_0, tg_zz_x_s_1_0_0, tg_zz_y_s_0_0_0, tg_zz_y_s_1_0_0, tg_zz_z_s_0_0_0, tg_zz_z_s_1_0_0, tg_zzz_x_s_0_0_0, tg_zzz_y_s_0_0_0, tg_zzz_z_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xxx_x_s_0_0_0[i] = 2.0 * tg_x_x_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_x_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_x_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_x_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_y_s_0_0_0[i] = 2.0 * tg_x_y_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_y_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_y_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_y_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_z_s_0_0_0[i] = 2.0 * tg_x_z_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_z_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_z_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_z_s_0_0_0[i] * a_x * faz_0;

        tg_xxy_x_s_0_0_0[i] = 2.0 * tg_xx_x_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_x_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_y_s_0_0_0[i] = 2.0 * tg_xx_y_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_y_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_z_s_0_0_0[i] = 2.0 * tg_xx_z_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_z_s_0_0_0[i] * a_y * faz_0;

        tg_xxz_x_s_0_0_0[i] = 2.0 * tg_xx_x_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_x_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_y_s_0_0_0[i] = 2.0 * tg_xx_y_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_y_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_z_s_0_0_0[i] = 2.0 * tg_xx_z_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_z_s_0_0_0[i] * a_z * faz_0;

        tg_xyy_x_s_0_0_0[i] = 2.0 * tg_yy_x_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_x_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_y_s_0_0_0[i] = 2.0 * tg_yy_y_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_y_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_z_s_0_0_0[i] = 2.0 * tg_yy_z_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_z_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_x_s_0_0_0[i] = 2.0 * tg_yz_x_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_x_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_y_s_0_0_0[i] = 2.0 * tg_yz_y_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_y_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_z_s_0_0_0[i] = 2.0 * tg_yz_z_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_z_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_x_s_0_0_0[i] = 2.0 * tg_zz_x_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_x_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_y_s_0_0_0[i] = 2.0 * tg_zz_y_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_y_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_z_s_0_0_0[i] = 2.0 * tg_zz_z_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_z_s_0_0_0[i] * a_x * faz_0;

        tg_yyy_x_s_0_0_0[i] = 2.0 * tg_y_x_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_x_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_x_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_x_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_y_s_0_0_0[i] = 2.0 * tg_y_y_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_y_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_y_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_y_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_z_s_0_0_0[i] = 2.0 * tg_y_z_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_z_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_z_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_z_s_0_0_0[i] * a_y * faz_0;

        tg_yyz_x_s_0_0_0[i] = 2.0 * tg_yy_x_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_x_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_y_s_0_0_0[i] = 2.0 * tg_yy_y_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_y_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_z_s_0_0_0[i] = 2.0 * tg_yy_z_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_z_s_0_0_0[i] * a_z * faz_0;

        tg_yzz_x_s_0_0_0[i] = 2.0 * tg_zz_x_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_x_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_y_s_0_0_0[i] = 2.0 * tg_zz_y_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_y_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_z_s_0_0_0[i] = 2.0 * tg_zz_z_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_z_s_0_0_0[i] * a_y * faz_0;

        tg_zzz_x_s_0_0_0[i] = 2.0 * tg_z_x_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_x_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_x_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_x_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_y_s_0_0_0[i] = 2.0 * tg_z_y_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_y_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_y_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_y_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_z_s_0_0_0[i] = 2.0 * tg_z_z_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_z_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_z_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_z_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : PP

        auto tg_x_x_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1);

        auto tg_x_y_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 1);

        auto tg_x_z_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 2);

        auto tg_y_x_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 3);

        auto tg_y_y_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 4);

        auto tg_y_z_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 5);

        auto tg_z_x_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 6);

        auto tg_z_y_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 7);

        auto tg_z_z_s_0_0_1 = pbuffer.data(idx_pp_s_0_0_1 + 8);

        // Set up components of auxiliary buffer : DP

        auto tg_xx_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1);

        auto tg_xx_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 1);

        auto tg_xx_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 2);

        auto tg_xy_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 3);

        auto tg_xy_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 4);

        auto tg_xy_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 5);

        auto tg_xz_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 6);

        auto tg_xz_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 7);

        auto tg_xz_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 8);

        auto tg_yy_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 9);

        auto tg_yy_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 10);

        auto tg_yy_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 11);

        auto tg_yz_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 12);

        auto tg_yz_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 13);

        auto tg_yz_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 14);

        auto tg_zz_x_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 15);

        auto tg_zz_y_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 16);

        auto tg_zz_z_s_0_0_1 = pbuffer.data(idx_dp_s_0_0_1 + 17);

        #pragma omp simd aligned(b_exps, tg_x_x_s_0_0_1, tg_x_y_s_0_0_1, tg_x_z_s_0_0_1, tg_xx_x_s_0_0_1, tg_xx_y_s_0_0_1, tg_xx_z_s_0_0_1, tg_xxx_x_s_0_0_0, tg_xxx_y_s_0_0_0, tg_xxx_z_s_0_0_0, tg_xxy_x_s_0_0_0, tg_xxy_y_s_0_0_0, tg_xxy_z_s_0_0_0, tg_xxz_x_s_0_0_0, tg_xxz_y_s_0_0_0, tg_xxz_z_s_0_0_0, tg_xyy_x_s_0_0_0, tg_xyy_y_s_0_0_0, tg_xyy_z_s_0_0_0, tg_xyz_x_s_0_0_0, tg_xyz_y_s_0_0_0, tg_xyz_z_s_0_0_0, tg_xzz_x_s_0_0_0, tg_xzz_y_s_0_0_0, tg_xzz_z_s_0_0_0, tg_y_x_s_0_0_1, tg_y_y_s_0_0_1, tg_y_z_s_0_0_1, tg_yy_x_s_0_0_1, tg_yy_y_s_0_0_1, tg_yy_z_s_0_0_1, tg_yyy_x_s_0_0_0, tg_yyy_y_s_0_0_0, tg_yyy_z_s_0_0_0, tg_yyz_x_s_0_0_0, tg_yyz_y_s_0_0_0, tg_yyz_z_s_0_0_0, tg_yz_x_s_0_0_1, tg_yz_y_s_0_0_1, tg_yz_z_s_0_0_1, tg_yzz_x_s_0_0_0, tg_yzz_y_s_0_0_0, tg_yzz_z_s_0_0_0, tg_z_x_s_0_0_1, tg_z_y_s_0_0_1, tg_z_z_s_0_0_1, tg_zz_x_s_0_0_1, tg_zz_y_s_0_0_1, tg_zz_z_s_0_0_1, tg_zzz_x_s_0_0_0, tg_zzz_y_s_0_0_0, tg_zzz_z_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxx_x_s_0_0_0[i] = tg_x_x_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_x_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_y_s_0_0_0[i] = tg_x_y_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_y_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_z_s_0_0_0[i] = tg_x_z_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_z_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxy_x_s_0_0_0[i] = tg_xx_x_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_y_s_0_0_0[i] = tg_xx_y_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_z_s_0_0_0[i] = tg_xx_z_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxz_x_s_0_0_0[i] = tg_xx_x_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_y_s_0_0_0[i] = tg_xx_y_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_z_s_0_0_0[i] = tg_xx_z_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xyy_x_s_0_0_0[i] = tg_yy_x_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_y_s_0_0_0[i] = tg_yy_y_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_z_s_0_0_0[i] = tg_yy_z_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_x_s_0_0_0[i] = tg_yz_x_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_y_s_0_0_0[i] = tg_yz_y_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_z_s_0_0_0[i] = tg_yz_z_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_x_s_0_0_0[i] = tg_zz_x_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_y_s_0_0_0[i] = tg_zz_y_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_z_s_0_0_0[i] = tg_zz_z_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyy_x_s_0_0_0[i] = tg_y_x_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_x_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_y_s_0_0_0[i] = tg_y_y_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_y_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_z_s_0_0_0[i] = tg_y_z_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_z_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyz_x_s_0_0_0[i] = tg_yy_x_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_y_s_0_0_0[i] = tg_yy_y_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_z_s_0_0_0[i] = tg_yy_z_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yzz_x_s_0_0_0[i] = tg_zz_x_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_y_s_0_0_0[i] = tg_zz_y_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_z_s_0_0_0[i] = tg_zz_z_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzz_x_s_0_0_0[i] = tg_z_x_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_x_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_y_s_0_0_0[i] = tg_z_y_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_y_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_z_s_0_0_0[i] = tg_z_z_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_z_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

