#include "ElectronRepulsionGeom0010ContrRecXXHP.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhp(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhp,
                                            const size_t idx_xxgp,
                                            const size_t idx_geom_10_xxgp,
                                            const size_t idx_geom_10_xxgd,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSGP

            const auto gp_off = idx_xxgp + (i * bcomps + j) * 45;

            auto g_xxxx_x = cbuffer.data(gp_off + 0);

            auto g_xxxx_y = cbuffer.data(gp_off + 1);

            auto g_xxxx_z = cbuffer.data(gp_off + 2);

            auto g_yyyy_x = cbuffer.data(gp_off + 30);

            auto g_yyyy_y = cbuffer.data(gp_off + 31);

            auto g_yyyy_z = cbuffer.data(gp_off + 32);

            auto g_zzzz_x = cbuffer.data(gp_off + 42);

            auto g_zzzz_y = cbuffer.data(gp_off + 43);

            auto g_zzzz_z = cbuffer.data(gp_off + 44);

            /// Set up components of auxilary buffer : SSGP

            const auto gp_geom_10_off = idx_geom_10_xxgp + (i * bcomps + j) * 45;

            auto g_x_0_xxxx_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xyyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xyyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xyyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xyyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xyzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xyzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_yyyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_yyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_yyyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_yyyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_yyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_yyyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_yyzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_yyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_yyzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_yzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_yzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_yzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_zzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_zzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_zzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_y_0_xxxx_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 0);

            auto g_y_0_xxxx_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 1);

            auto g_y_0_xxxx_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 2);

            auto g_y_0_xxxy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 3);

            auto g_y_0_xxxy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 4);

            auto g_y_0_xxxy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 5);

            auto g_y_0_xxxz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 6);

            auto g_y_0_xxxz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 7);

            auto g_y_0_xxxz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 8);

            auto g_y_0_xxyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 9);

            auto g_y_0_xxyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 10);

            auto g_y_0_xxyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 11);

            auto g_y_0_xxyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 12);

            auto g_y_0_xxyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 13);

            auto g_y_0_xxyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 14);

            auto g_y_0_xxzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 15);

            auto g_y_0_xxzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 16);

            auto g_y_0_xxzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 17);

            auto g_y_0_xyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 18);

            auto g_y_0_xyyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 19);

            auto g_y_0_xyyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 20);

            auto g_y_0_xyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 21);

            auto g_y_0_xyyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 22);

            auto g_y_0_xyyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 23);

            auto g_y_0_xyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 24);

            auto g_y_0_xyzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 25);

            auto g_y_0_xyzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 26);

            auto g_y_0_xzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 27);

            auto g_y_0_xzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 28);

            auto g_y_0_xzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 29);

            auto g_y_0_yyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 30);

            auto g_y_0_yyyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 31);

            auto g_y_0_yyyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 32);

            auto g_y_0_yyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 33);

            auto g_y_0_yyyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 34);

            auto g_y_0_yyyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 35);

            auto g_y_0_yyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 36);

            auto g_y_0_yyzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 37);

            auto g_y_0_yyzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 38);

            auto g_y_0_yzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 39);

            auto g_y_0_yzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 40);

            auto g_y_0_yzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 41);

            auto g_y_0_zzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 42);

            auto g_y_0_zzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 43);

            auto g_y_0_zzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 44);

            auto g_z_0_xxxx_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_z_0_xxxx_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 1);

            auto g_z_0_xxxx_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 2);

            auto g_z_0_xxxy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 3);

            auto g_z_0_xxxy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 4);

            auto g_z_0_xxxy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 5);

            auto g_z_0_xxxz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_z_0_xxxz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 7);

            auto g_z_0_xxxz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 8);

            auto g_z_0_xxyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 9);

            auto g_z_0_xxyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 10);

            auto g_z_0_xxyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 11);

            auto g_z_0_xxyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_z_0_xxyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 13);

            auto g_z_0_xxyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 14);

            auto g_z_0_xxzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 15);

            auto g_z_0_xxzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 16);

            auto g_z_0_xxzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 17);

            auto g_z_0_xyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_z_0_xyyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 19);

            auto g_z_0_xyyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 20);

            auto g_z_0_xyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 21);

            auto g_z_0_xyyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 22);

            auto g_z_0_xyyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 23);

            auto g_z_0_xyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_z_0_xyzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 25);

            auto g_z_0_xyzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 26);

            auto g_z_0_xzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 27);

            auto g_z_0_xzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 28);

            auto g_z_0_xzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 29);

            auto g_z_0_yyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_z_0_yyyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_z_0_yyyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 32);

            auto g_z_0_yyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 33);

            auto g_z_0_yyyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 34);

            auto g_z_0_yyyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 35);

            auto g_z_0_yyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_z_0_yyzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_z_0_yyzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 38);

            auto g_z_0_yzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 39);

            auto g_z_0_yzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 40);

            auto g_z_0_yzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 41);

            auto g_z_0_zzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_z_0_zzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_z_0_zzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 44);

            /// Set up components of auxilary buffer : SSGD

            const auto gd_geom_10_off = idx_geom_10_xxgd + (i * bcomps + j) * 90;

            auto g_x_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_y_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 2);

            auto g_y_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_y_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 7);

            auto g_y_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 8);

            auto g_y_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_y_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 13);

            auto g_y_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 14);

            auto g_y_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_y_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 19);

            auto g_y_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 20);

            auto g_y_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_y_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 25);

            auto g_y_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 26);

            auto g_y_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_y_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_y_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 32);

            auto g_y_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_y_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_y_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 38);

            auto g_y_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_y_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_y_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 44);

            auto g_y_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 48);

            auto g_y_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 49);

            auto g_y_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 50);

            auto g_y_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 54);

            auto g_y_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 55);

            auto g_y_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 56);

            auto g_y_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 60);

            auto g_y_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 61);

            auto g_y_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 62);

            auto g_y_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 63);

            auto g_y_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 64);

            auto g_y_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 65);

            auto g_y_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 66);

            auto g_y_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 67);

            auto g_y_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 68);

            auto g_y_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 70);

            auto g_y_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 71);

            auto g_y_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 72);

            auto g_y_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 73);

            auto g_y_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 74);

            auto g_y_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 76);

            auto g_y_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 77);

            auto g_y_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 78);

            auto g_y_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 79);

            auto g_y_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 80);

            auto g_y_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 82);

            auto g_y_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 83);

            auto g_y_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 84);

            auto g_y_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 85);

            auto g_y_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 86);

            auto g_y_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 88);

            auto g_y_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 89);

            auto g_z_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 2);

            auto g_z_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 6);

            auto g_z_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 7);

            auto g_z_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 8);

            auto g_z_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 12);

            auto g_z_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 13);

            auto g_z_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 14);

            auto g_z_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 18);

            auto g_z_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 19);

            auto g_z_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 20);

            auto g_z_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 24);

            auto g_z_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 25);

            auto g_z_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 26);

            auto g_z_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 30);

            auto g_z_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 31);

            auto g_z_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 32);

            auto g_z_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 36);

            auto g_z_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 37);

            auto g_z_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 38);

            auto g_z_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 42);

            auto g_z_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 43);

            auto g_z_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 44);

            auto g_z_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 48);

            auto g_z_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 49);

            auto g_z_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 50);

            auto g_z_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 54);

            auto g_z_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 55);

            auto g_z_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 56);

            auto g_z_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 60);

            auto g_z_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 61);

            auto g_z_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 62);

            auto g_z_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 63);

            auto g_z_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 64);

            auto g_z_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 66);

            auto g_z_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 67);

            auto g_z_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 68);

            auto g_z_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 69);

            auto g_z_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 70);

            auto g_z_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 72);

            auto g_z_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 73);

            auto g_z_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 74);

            auto g_z_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 75);

            auto g_z_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 76);

            auto g_z_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 78);

            auto g_z_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 79);

            auto g_z_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 80);

            auto g_z_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 81);

            auto g_z_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 82);

            auto g_z_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 84);

            auto g_z_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 85);

            auto g_z_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 86);

            auto g_z_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 87);

            auto g_z_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 88);

            auto g_z_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 89);

            /// set up bra offset for contr_buffer_xxhp

            const auto hp_geom_10_off = idx_geom_10_xxhp + (i * bcomps + j) * 63;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_x, g_x_0_xxxx_xx, g_x_0_xxxx_xy, g_x_0_xxxx_xz, g_x_0_xxxx_y, g_x_0_xxxx_z, g_x_0_xxxxx_x, g_x_0_xxxxx_y, g_x_0_xxxxx_z, g_xxxx_x, g_xxxx_y, g_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_x[k] = -g_xxxx_x[k] - g_x_0_xxxx_x[k] * cd_x[k] + g_x_0_xxxx_xx[k];

                g_x_0_xxxxx_y[k] = -g_xxxx_y[k] - g_x_0_xxxx_y[k] * cd_x[k] + g_x_0_xxxx_xy[k];

                g_x_0_xxxxx_z[k] = -g_xxxx_z[k] - g_x_0_xxxx_z[k] * cd_x[k] + g_x_0_xxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_x, g_x_0_xxxx_xy, g_x_0_xxxx_y, g_x_0_xxxx_yy, g_x_0_xxxx_yz, g_x_0_xxxx_z, g_x_0_xxxxy_x, g_x_0_xxxxy_y, g_x_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_x[k] = -g_x_0_xxxx_x[k] * cd_y[k] + g_x_0_xxxx_xy[k];

                g_x_0_xxxxy_y[k] = -g_x_0_xxxx_y[k] * cd_y[k] + g_x_0_xxxx_yy[k];

                g_x_0_xxxxy_z[k] = -g_x_0_xxxx_z[k] * cd_y[k] + g_x_0_xxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_x, g_x_0_xxxx_xz, g_x_0_xxxx_y, g_x_0_xxxx_yz, g_x_0_xxxx_z, g_x_0_xxxx_zz, g_x_0_xxxxz_x, g_x_0_xxxxz_y, g_x_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_x[k] = -g_x_0_xxxx_x[k] * cd_z[k] + g_x_0_xxxx_xz[k];

                g_x_0_xxxxz_y[k] = -g_x_0_xxxx_y[k] * cd_z[k] + g_x_0_xxxx_yz[k];

                g_x_0_xxxxz_z[k] = -g_x_0_xxxx_z[k] * cd_z[k] + g_x_0_xxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_x, g_x_0_xxxy_xy, g_x_0_xxxy_y, g_x_0_xxxy_yy, g_x_0_xxxy_yz, g_x_0_xxxy_z, g_x_0_xxxyy_x, g_x_0_xxxyy_y, g_x_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_x[k] = -g_x_0_xxxy_x[k] * cd_y[k] + g_x_0_xxxy_xy[k];

                g_x_0_xxxyy_y[k] = -g_x_0_xxxy_y[k] * cd_y[k] + g_x_0_xxxy_yy[k];

                g_x_0_xxxyy_z[k] = -g_x_0_xxxy_z[k] * cd_y[k] + g_x_0_xxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_x, g_x_0_xxxyz_y, g_x_0_xxxyz_z, g_x_0_xxxz_x, g_x_0_xxxz_xy, g_x_0_xxxz_y, g_x_0_xxxz_yy, g_x_0_xxxz_yz, g_x_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_x[k] = -g_x_0_xxxz_x[k] * cd_y[k] + g_x_0_xxxz_xy[k];

                g_x_0_xxxyz_y[k] = -g_x_0_xxxz_y[k] * cd_y[k] + g_x_0_xxxz_yy[k];

                g_x_0_xxxyz_z[k] = -g_x_0_xxxz_z[k] * cd_y[k] + g_x_0_xxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_x, g_x_0_xxxz_xz, g_x_0_xxxz_y, g_x_0_xxxz_yz, g_x_0_xxxz_z, g_x_0_xxxz_zz, g_x_0_xxxzz_x, g_x_0_xxxzz_y, g_x_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_x[k] = -g_x_0_xxxz_x[k] * cd_z[k] + g_x_0_xxxz_xz[k];

                g_x_0_xxxzz_y[k] = -g_x_0_xxxz_y[k] * cd_z[k] + g_x_0_xxxz_yz[k];

                g_x_0_xxxzz_z[k] = -g_x_0_xxxz_z[k] * cd_z[k] + g_x_0_xxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_x, g_x_0_xxyy_xy, g_x_0_xxyy_y, g_x_0_xxyy_yy, g_x_0_xxyy_yz, g_x_0_xxyy_z, g_x_0_xxyyy_x, g_x_0_xxyyy_y, g_x_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_x[k] = -g_x_0_xxyy_x[k] * cd_y[k] + g_x_0_xxyy_xy[k];

                g_x_0_xxyyy_y[k] = -g_x_0_xxyy_y[k] * cd_y[k] + g_x_0_xxyy_yy[k];

                g_x_0_xxyyy_z[k] = -g_x_0_xxyy_z[k] * cd_y[k] + g_x_0_xxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_x, g_x_0_xxyyz_y, g_x_0_xxyyz_z, g_x_0_xxyz_x, g_x_0_xxyz_xy, g_x_0_xxyz_y, g_x_0_xxyz_yy, g_x_0_xxyz_yz, g_x_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_x[k] = -g_x_0_xxyz_x[k] * cd_y[k] + g_x_0_xxyz_xy[k];

                g_x_0_xxyyz_y[k] = -g_x_0_xxyz_y[k] * cd_y[k] + g_x_0_xxyz_yy[k];

                g_x_0_xxyyz_z[k] = -g_x_0_xxyz_z[k] * cd_y[k] + g_x_0_xxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_x, g_x_0_xxyzz_y, g_x_0_xxyzz_z, g_x_0_xxzz_x, g_x_0_xxzz_xy, g_x_0_xxzz_y, g_x_0_xxzz_yy, g_x_0_xxzz_yz, g_x_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_x[k] = -g_x_0_xxzz_x[k] * cd_y[k] + g_x_0_xxzz_xy[k];

                g_x_0_xxyzz_y[k] = -g_x_0_xxzz_y[k] * cd_y[k] + g_x_0_xxzz_yy[k];

                g_x_0_xxyzz_z[k] = -g_x_0_xxzz_z[k] * cd_y[k] + g_x_0_xxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_x, g_x_0_xxzz_xz, g_x_0_xxzz_y, g_x_0_xxzz_yz, g_x_0_xxzz_z, g_x_0_xxzz_zz, g_x_0_xxzzz_x, g_x_0_xxzzz_y, g_x_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_x[k] = -g_x_0_xxzz_x[k] * cd_z[k] + g_x_0_xxzz_xz[k];

                g_x_0_xxzzz_y[k] = -g_x_0_xxzz_y[k] * cd_z[k] + g_x_0_xxzz_yz[k];

                g_x_0_xxzzz_z[k] = -g_x_0_xxzz_z[k] * cd_z[k] + g_x_0_xxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_x, g_x_0_xyyy_xy, g_x_0_xyyy_y, g_x_0_xyyy_yy, g_x_0_xyyy_yz, g_x_0_xyyy_z, g_x_0_xyyyy_x, g_x_0_xyyyy_y, g_x_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_x[k] = -g_x_0_xyyy_x[k] * cd_y[k] + g_x_0_xyyy_xy[k];

                g_x_0_xyyyy_y[k] = -g_x_0_xyyy_y[k] * cd_y[k] + g_x_0_xyyy_yy[k];

                g_x_0_xyyyy_z[k] = -g_x_0_xyyy_z[k] * cd_y[k] + g_x_0_xyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_x, g_x_0_xyyyz_y, g_x_0_xyyyz_z, g_x_0_xyyz_x, g_x_0_xyyz_xy, g_x_0_xyyz_y, g_x_0_xyyz_yy, g_x_0_xyyz_yz, g_x_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_x[k] = -g_x_0_xyyz_x[k] * cd_y[k] + g_x_0_xyyz_xy[k];

                g_x_0_xyyyz_y[k] = -g_x_0_xyyz_y[k] * cd_y[k] + g_x_0_xyyz_yy[k];

                g_x_0_xyyyz_z[k] = -g_x_0_xyyz_z[k] * cd_y[k] + g_x_0_xyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_x, g_x_0_xyyzz_y, g_x_0_xyyzz_z, g_x_0_xyzz_x, g_x_0_xyzz_xy, g_x_0_xyzz_y, g_x_0_xyzz_yy, g_x_0_xyzz_yz, g_x_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_x[k] = -g_x_0_xyzz_x[k] * cd_y[k] + g_x_0_xyzz_xy[k];

                g_x_0_xyyzz_y[k] = -g_x_0_xyzz_y[k] * cd_y[k] + g_x_0_xyzz_yy[k];

                g_x_0_xyyzz_z[k] = -g_x_0_xyzz_z[k] * cd_y[k] + g_x_0_xyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_x, g_x_0_xyzzz_y, g_x_0_xyzzz_z, g_x_0_xzzz_x, g_x_0_xzzz_xy, g_x_0_xzzz_y, g_x_0_xzzz_yy, g_x_0_xzzz_yz, g_x_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_x[k] = -g_x_0_xzzz_x[k] * cd_y[k] + g_x_0_xzzz_xy[k];

                g_x_0_xyzzz_y[k] = -g_x_0_xzzz_y[k] * cd_y[k] + g_x_0_xzzz_yy[k];

                g_x_0_xyzzz_z[k] = -g_x_0_xzzz_z[k] * cd_y[k] + g_x_0_xzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_x, g_x_0_xzzz_xz, g_x_0_xzzz_y, g_x_0_xzzz_yz, g_x_0_xzzz_z, g_x_0_xzzz_zz, g_x_0_xzzzz_x, g_x_0_xzzzz_y, g_x_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_x[k] = -g_x_0_xzzz_x[k] * cd_z[k] + g_x_0_xzzz_xz[k];

                g_x_0_xzzzz_y[k] = -g_x_0_xzzz_y[k] * cd_z[k] + g_x_0_xzzz_yz[k];

                g_x_0_xzzzz_z[k] = -g_x_0_xzzz_z[k] * cd_z[k] + g_x_0_xzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_x, g_x_0_yyyy_xy, g_x_0_yyyy_y, g_x_0_yyyy_yy, g_x_0_yyyy_yz, g_x_0_yyyy_z, g_x_0_yyyyy_x, g_x_0_yyyyy_y, g_x_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_x[k] = -g_x_0_yyyy_x[k] * cd_y[k] + g_x_0_yyyy_xy[k];

                g_x_0_yyyyy_y[k] = -g_x_0_yyyy_y[k] * cd_y[k] + g_x_0_yyyy_yy[k];

                g_x_0_yyyyy_z[k] = -g_x_0_yyyy_z[k] * cd_y[k] + g_x_0_yyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_x, g_x_0_yyyyz_y, g_x_0_yyyyz_z, g_x_0_yyyz_x, g_x_0_yyyz_xy, g_x_0_yyyz_y, g_x_0_yyyz_yy, g_x_0_yyyz_yz, g_x_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_x[k] = -g_x_0_yyyz_x[k] * cd_y[k] + g_x_0_yyyz_xy[k];

                g_x_0_yyyyz_y[k] = -g_x_0_yyyz_y[k] * cd_y[k] + g_x_0_yyyz_yy[k];

                g_x_0_yyyyz_z[k] = -g_x_0_yyyz_z[k] * cd_y[k] + g_x_0_yyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_x, g_x_0_yyyzz_y, g_x_0_yyyzz_z, g_x_0_yyzz_x, g_x_0_yyzz_xy, g_x_0_yyzz_y, g_x_0_yyzz_yy, g_x_0_yyzz_yz, g_x_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_x[k] = -g_x_0_yyzz_x[k] * cd_y[k] + g_x_0_yyzz_xy[k];

                g_x_0_yyyzz_y[k] = -g_x_0_yyzz_y[k] * cd_y[k] + g_x_0_yyzz_yy[k];

                g_x_0_yyyzz_z[k] = -g_x_0_yyzz_z[k] * cd_y[k] + g_x_0_yyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_x, g_x_0_yyzzz_y, g_x_0_yyzzz_z, g_x_0_yzzz_x, g_x_0_yzzz_xy, g_x_0_yzzz_y, g_x_0_yzzz_yy, g_x_0_yzzz_yz, g_x_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_x[k] = -g_x_0_yzzz_x[k] * cd_y[k] + g_x_0_yzzz_xy[k];

                g_x_0_yyzzz_y[k] = -g_x_0_yzzz_y[k] * cd_y[k] + g_x_0_yzzz_yy[k];

                g_x_0_yyzzz_z[k] = -g_x_0_yzzz_z[k] * cd_y[k] + g_x_0_yzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_x, g_x_0_yzzzz_y, g_x_0_yzzzz_z, g_x_0_zzzz_x, g_x_0_zzzz_xy, g_x_0_zzzz_y, g_x_0_zzzz_yy, g_x_0_zzzz_yz, g_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_x[k] = -g_x_0_zzzz_x[k] * cd_y[k] + g_x_0_zzzz_xy[k];

                g_x_0_yzzzz_y[k] = -g_x_0_zzzz_y[k] * cd_y[k] + g_x_0_zzzz_yy[k];

                g_x_0_yzzzz_z[k] = -g_x_0_zzzz_z[k] * cd_y[k] + g_x_0_zzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_x, g_x_0_zzzz_xz, g_x_0_zzzz_y, g_x_0_zzzz_yz, g_x_0_zzzz_z, g_x_0_zzzz_zz, g_x_0_zzzzz_x, g_x_0_zzzzz_y, g_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_x[k] = -g_x_0_zzzz_x[k] * cd_z[k] + g_x_0_zzzz_xz[k];

                g_x_0_zzzzz_y[k] = -g_x_0_zzzz_y[k] * cd_z[k] + g_x_0_zzzz_yz[k];

                g_x_0_zzzzz_z[k] = -g_x_0_zzzz_z[k] * cd_z[k] + g_x_0_zzzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_x, g_y_0_xxxx_xx, g_y_0_xxxx_xy, g_y_0_xxxx_xz, g_y_0_xxxx_y, g_y_0_xxxx_z, g_y_0_xxxxx_x, g_y_0_xxxxx_y, g_y_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_x[k] = -g_y_0_xxxx_x[k] * cd_x[k] + g_y_0_xxxx_xx[k];

                g_y_0_xxxxx_y[k] = -g_y_0_xxxx_y[k] * cd_x[k] + g_y_0_xxxx_xy[k];

                g_y_0_xxxxx_z[k] = -g_y_0_xxxx_z[k] * cd_x[k] + g_y_0_xxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 3);

            auto g_y_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 4);

            auto g_y_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_x, g_y_0_xxxxy_y, g_y_0_xxxxy_z, g_y_0_xxxy_x, g_y_0_xxxy_xx, g_y_0_xxxy_xy, g_y_0_xxxy_xz, g_y_0_xxxy_y, g_y_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_x[k] = -g_y_0_xxxy_x[k] * cd_x[k] + g_y_0_xxxy_xx[k];

                g_y_0_xxxxy_y[k] = -g_y_0_xxxy_y[k] * cd_x[k] + g_y_0_xxxy_xy[k];

                g_y_0_xxxxy_z[k] = -g_y_0_xxxy_z[k] * cd_x[k] + g_y_0_xxxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 6);

            auto g_y_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 7);

            auto g_y_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_x, g_y_0_xxxxz_y, g_y_0_xxxxz_z, g_y_0_xxxz_x, g_y_0_xxxz_xx, g_y_0_xxxz_xy, g_y_0_xxxz_xz, g_y_0_xxxz_y, g_y_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_x[k] = -g_y_0_xxxz_x[k] * cd_x[k] + g_y_0_xxxz_xx[k];

                g_y_0_xxxxz_y[k] = -g_y_0_xxxz_y[k] * cd_x[k] + g_y_0_xxxz_xy[k];

                g_y_0_xxxxz_z[k] = -g_y_0_xxxz_z[k] * cd_x[k] + g_y_0_xxxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 9);

            auto g_y_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 10);

            auto g_y_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_x, g_y_0_xxxyy_y, g_y_0_xxxyy_z, g_y_0_xxyy_x, g_y_0_xxyy_xx, g_y_0_xxyy_xy, g_y_0_xxyy_xz, g_y_0_xxyy_y, g_y_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_x[k] = -g_y_0_xxyy_x[k] * cd_x[k] + g_y_0_xxyy_xx[k];

                g_y_0_xxxyy_y[k] = -g_y_0_xxyy_y[k] * cd_x[k] + g_y_0_xxyy_xy[k];

                g_y_0_xxxyy_z[k] = -g_y_0_xxyy_z[k] * cd_x[k] + g_y_0_xxyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 12);

            auto g_y_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 13);

            auto g_y_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_x, g_y_0_xxxyz_y, g_y_0_xxxyz_z, g_y_0_xxyz_x, g_y_0_xxyz_xx, g_y_0_xxyz_xy, g_y_0_xxyz_xz, g_y_0_xxyz_y, g_y_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_x[k] = -g_y_0_xxyz_x[k] * cd_x[k] + g_y_0_xxyz_xx[k];

                g_y_0_xxxyz_y[k] = -g_y_0_xxyz_y[k] * cd_x[k] + g_y_0_xxyz_xy[k];

                g_y_0_xxxyz_z[k] = -g_y_0_xxyz_z[k] * cd_x[k] + g_y_0_xxyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 15);

            auto g_y_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 16);

            auto g_y_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_x, g_y_0_xxxzz_y, g_y_0_xxxzz_z, g_y_0_xxzz_x, g_y_0_xxzz_xx, g_y_0_xxzz_xy, g_y_0_xxzz_xz, g_y_0_xxzz_y, g_y_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_x[k] = -g_y_0_xxzz_x[k] * cd_x[k] + g_y_0_xxzz_xx[k];

                g_y_0_xxxzz_y[k] = -g_y_0_xxzz_y[k] * cd_x[k] + g_y_0_xxzz_xy[k];

                g_y_0_xxxzz_z[k] = -g_y_0_xxzz_z[k] * cd_x[k] + g_y_0_xxzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 18);

            auto g_y_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 19);

            auto g_y_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_x, g_y_0_xxyyy_y, g_y_0_xxyyy_z, g_y_0_xyyy_x, g_y_0_xyyy_xx, g_y_0_xyyy_xy, g_y_0_xyyy_xz, g_y_0_xyyy_y, g_y_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_x[k] = -g_y_0_xyyy_x[k] * cd_x[k] + g_y_0_xyyy_xx[k];

                g_y_0_xxyyy_y[k] = -g_y_0_xyyy_y[k] * cd_x[k] + g_y_0_xyyy_xy[k];

                g_y_0_xxyyy_z[k] = -g_y_0_xyyy_z[k] * cd_x[k] + g_y_0_xyyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 21);

            auto g_y_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 22);

            auto g_y_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_x, g_y_0_xxyyz_y, g_y_0_xxyyz_z, g_y_0_xyyz_x, g_y_0_xyyz_xx, g_y_0_xyyz_xy, g_y_0_xyyz_xz, g_y_0_xyyz_y, g_y_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_x[k] = -g_y_0_xyyz_x[k] * cd_x[k] + g_y_0_xyyz_xx[k];

                g_y_0_xxyyz_y[k] = -g_y_0_xyyz_y[k] * cd_x[k] + g_y_0_xyyz_xy[k];

                g_y_0_xxyyz_z[k] = -g_y_0_xyyz_z[k] * cd_x[k] + g_y_0_xyyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 24);

            auto g_y_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 25);

            auto g_y_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_x, g_y_0_xxyzz_y, g_y_0_xxyzz_z, g_y_0_xyzz_x, g_y_0_xyzz_xx, g_y_0_xyzz_xy, g_y_0_xyzz_xz, g_y_0_xyzz_y, g_y_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_x[k] = -g_y_0_xyzz_x[k] * cd_x[k] + g_y_0_xyzz_xx[k];

                g_y_0_xxyzz_y[k] = -g_y_0_xyzz_y[k] * cd_x[k] + g_y_0_xyzz_xy[k];

                g_y_0_xxyzz_z[k] = -g_y_0_xyzz_z[k] * cd_x[k] + g_y_0_xyzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 27);

            auto g_y_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 28);

            auto g_y_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_x, g_y_0_xxzzz_y, g_y_0_xxzzz_z, g_y_0_xzzz_x, g_y_0_xzzz_xx, g_y_0_xzzz_xy, g_y_0_xzzz_xz, g_y_0_xzzz_y, g_y_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_x[k] = -g_y_0_xzzz_x[k] * cd_x[k] + g_y_0_xzzz_xx[k];

                g_y_0_xxzzz_y[k] = -g_y_0_xzzz_y[k] * cd_x[k] + g_y_0_xzzz_xy[k];

                g_y_0_xxzzz_z[k] = -g_y_0_xzzz_z[k] * cd_x[k] + g_y_0_xzzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 30);

            auto g_y_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 31);

            auto g_y_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_x, g_y_0_xyyyy_y, g_y_0_xyyyy_z, g_y_0_yyyy_x, g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_y, g_y_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_x[k] = -g_y_0_yyyy_x[k] * cd_x[k] + g_y_0_yyyy_xx[k];

                g_y_0_xyyyy_y[k] = -g_y_0_yyyy_y[k] * cd_x[k] + g_y_0_yyyy_xy[k];

                g_y_0_xyyyy_z[k] = -g_y_0_yyyy_z[k] * cd_x[k] + g_y_0_yyyy_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 33);

            auto g_y_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 34);

            auto g_y_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_x, g_y_0_xyyyz_y, g_y_0_xyyyz_z, g_y_0_yyyz_x, g_y_0_yyyz_xx, g_y_0_yyyz_xy, g_y_0_yyyz_xz, g_y_0_yyyz_y, g_y_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_x[k] = -g_y_0_yyyz_x[k] * cd_x[k] + g_y_0_yyyz_xx[k];

                g_y_0_xyyyz_y[k] = -g_y_0_yyyz_y[k] * cd_x[k] + g_y_0_yyyz_xy[k];

                g_y_0_xyyyz_z[k] = -g_y_0_yyyz_z[k] * cd_x[k] + g_y_0_yyyz_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 36);

            auto g_y_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 37);

            auto g_y_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_x, g_y_0_xyyzz_y, g_y_0_xyyzz_z, g_y_0_yyzz_x, g_y_0_yyzz_xx, g_y_0_yyzz_xy, g_y_0_yyzz_xz, g_y_0_yyzz_y, g_y_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_x[k] = -g_y_0_yyzz_x[k] * cd_x[k] + g_y_0_yyzz_xx[k];

                g_y_0_xyyzz_y[k] = -g_y_0_yyzz_y[k] * cd_x[k] + g_y_0_yyzz_xy[k];

                g_y_0_xyyzz_z[k] = -g_y_0_yyzz_z[k] * cd_x[k] + g_y_0_yyzz_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 39);

            auto g_y_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 40);

            auto g_y_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_x, g_y_0_xyzzz_y, g_y_0_xyzzz_z, g_y_0_yzzz_x, g_y_0_yzzz_xx, g_y_0_yzzz_xy, g_y_0_yzzz_xz, g_y_0_yzzz_y, g_y_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_x[k] = -g_y_0_yzzz_x[k] * cd_x[k] + g_y_0_yzzz_xx[k];

                g_y_0_xyzzz_y[k] = -g_y_0_yzzz_y[k] * cd_x[k] + g_y_0_yzzz_xy[k];

                g_y_0_xyzzz_z[k] = -g_y_0_yzzz_z[k] * cd_x[k] + g_y_0_yzzz_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 42);

            auto g_y_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 43);

            auto g_y_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_x, g_y_0_xzzzz_y, g_y_0_xzzzz_z, g_y_0_zzzz_x, g_y_0_zzzz_xx, g_y_0_zzzz_xy, g_y_0_zzzz_xz, g_y_0_zzzz_y, g_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_x[k] = -g_y_0_zzzz_x[k] * cd_x[k] + g_y_0_zzzz_xx[k];

                g_y_0_xzzzz_y[k] = -g_y_0_zzzz_y[k] * cd_x[k] + g_y_0_zzzz_xy[k];

                g_y_0_xzzzz_z[k] = -g_y_0_zzzz_z[k] * cd_x[k] + g_y_0_zzzz_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 45);

            auto g_y_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 46);

            auto g_y_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_x, g_y_0_yyyy_xy, g_y_0_yyyy_y, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_z, g_y_0_yyyyy_x, g_y_0_yyyyy_y, g_y_0_yyyyy_z, g_yyyy_x, g_yyyy_y, g_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_x[k] = -g_yyyy_x[k] - g_y_0_yyyy_x[k] * cd_y[k] + g_y_0_yyyy_xy[k];

                g_y_0_yyyyy_y[k] = -g_yyyy_y[k] - g_y_0_yyyy_y[k] * cd_y[k] + g_y_0_yyyy_yy[k];

                g_y_0_yyyyy_z[k] = -g_yyyy_z[k] - g_y_0_yyyy_z[k] * cd_y[k] + g_y_0_yyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 48);

            auto g_y_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 49);

            auto g_y_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_x, g_y_0_yyyy_xz, g_y_0_yyyy_y, g_y_0_yyyy_yz, g_y_0_yyyy_z, g_y_0_yyyy_zz, g_y_0_yyyyz_x, g_y_0_yyyyz_y, g_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_x[k] = -g_y_0_yyyy_x[k] * cd_z[k] + g_y_0_yyyy_xz[k];

                g_y_0_yyyyz_y[k] = -g_y_0_yyyy_y[k] * cd_z[k] + g_y_0_yyyy_yz[k];

                g_y_0_yyyyz_z[k] = -g_y_0_yyyy_z[k] * cd_z[k] + g_y_0_yyyy_zz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 51);

            auto g_y_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 52);

            auto g_y_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_x, g_y_0_yyyz_xz, g_y_0_yyyz_y, g_y_0_yyyz_yz, g_y_0_yyyz_z, g_y_0_yyyz_zz, g_y_0_yyyzz_x, g_y_0_yyyzz_y, g_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_x[k] = -g_y_0_yyyz_x[k] * cd_z[k] + g_y_0_yyyz_xz[k];

                g_y_0_yyyzz_y[k] = -g_y_0_yyyz_y[k] * cd_z[k] + g_y_0_yyyz_yz[k];

                g_y_0_yyyzz_z[k] = -g_y_0_yyyz_z[k] * cd_z[k] + g_y_0_yyyz_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 54);

            auto g_y_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 55);

            auto g_y_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_x, g_y_0_yyzz_xz, g_y_0_yyzz_y, g_y_0_yyzz_yz, g_y_0_yyzz_z, g_y_0_yyzz_zz, g_y_0_yyzzz_x, g_y_0_yyzzz_y, g_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_x[k] = -g_y_0_yyzz_x[k] * cd_z[k] + g_y_0_yyzz_xz[k];

                g_y_0_yyzzz_y[k] = -g_y_0_yyzz_y[k] * cd_z[k] + g_y_0_yyzz_yz[k];

                g_y_0_yyzzz_z[k] = -g_y_0_yyzz_z[k] * cd_z[k] + g_y_0_yyzz_zz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 57);

            auto g_y_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 58);

            auto g_y_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_x, g_y_0_yzzz_xz, g_y_0_yzzz_y, g_y_0_yzzz_yz, g_y_0_yzzz_z, g_y_0_yzzz_zz, g_y_0_yzzzz_x, g_y_0_yzzzz_y, g_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_x[k] = -g_y_0_yzzz_x[k] * cd_z[k] + g_y_0_yzzz_xz[k];

                g_y_0_yzzzz_y[k] = -g_y_0_yzzz_y[k] * cd_z[k] + g_y_0_yzzz_yz[k];

                g_y_0_yzzzz_z[k] = -g_y_0_yzzz_z[k] * cd_z[k] + g_y_0_yzzz_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 60);

            auto g_y_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 61);

            auto g_y_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_x, g_y_0_zzzz_xz, g_y_0_zzzz_y, g_y_0_zzzz_yz, g_y_0_zzzz_z, g_y_0_zzzz_zz, g_y_0_zzzzz_x, g_y_0_zzzzz_y, g_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_x[k] = -g_y_0_zzzz_x[k] * cd_z[k] + g_y_0_zzzz_xz[k];

                g_y_0_zzzzz_y[k] = -g_y_0_zzzz_y[k] * cd_z[k] + g_y_0_zzzz_yz[k];

                g_y_0_zzzzz_z[k] = -g_y_0_zzzz_z[k] * cd_z[k] + g_y_0_zzzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_x, g_z_0_xxxx_xx, g_z_0_xxxx_xy, g_z_0_xxxx_xz, g_z_0_xxxx_y, g_z_0_xxxx_z, g_z_0_xxxxx_x, g_z_0_xxxxx_y, g_z_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_x[k] = -g_z_0_xxxx_x[k] * cd_x[k] + g_z_0_xxxx_xx[k];

                g_z_0_xxxxx_y[k] = -g_z_0_xxxx_y[k] * cd_x[k] + g_z_0_xxxx_xy[k];

                g_z_0_xxxxx_z[k] = -g_z_0_xxxx_z[k] * cd_x[k] + g_z_0_xxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 3);

            auto g_z_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 4);

            auto g_z_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_x, g_z_0_xxxxy_y, g_z_0_xxxxy_z, g_z_0_xxxy_x, g_z_0_xxxy_xx, g_z_0_xxxy_xy, g_z_0_xxxy_xz, g_z_0_xxxy_y, g_z_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_x[k] = -g_z_0_xxxy_x[k] * cd_x[k] + g_z_0_xxxy_xx[k];

                g_z_0_xxxxy_y[k] = -g_z_0_xxxy_y[k] * cd_x[k] + g_z_0_xxxy_xy[k];

                g_z_0_xxxxy_z[k] = -g_z_0_xxxy_z[k] * cd_x[k] + g_z_0_xxxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 6);

            auto g_z_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 7);

            auto g_z_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_x, g_z_0_xxxxz_y, g_z_0_xxxxz_z, g_z_0_xxxz_x, g_z_0_xxxz_xx, g_z_0_xxxz_xy, g_z_0_xxxz_xz, g_z_0_xxxz_y, g_z_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_x[k] = -g_z_0_xxxz_x[k] * cd_x[k] + g_z_0_xxxz_xx[k];

                g_z_0_xxxxz_y[k] = -g_z_0_xxxz_y[k] * cd_x[k] + g_z_0_xxxz_xy[k];

                g_z_0_xxxxz_z[k] = -g_z_0_xxxz_z[k] * cd_x[k] + g_z_0_xxxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 9);

            auto g_z_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 10);

            auto g_z_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_x, g_z_0_xxxyy_y, g_z_0_xxxyy_z, g_z_0_xxyy_x, g_z_0_xxyy_xx, g_z_0_xxyy_xy, g_z_0_xxyy_xz, g_z_0_xxyy_y, g_z_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_x[k] = -g_z_0_xxyy_x[k] * cd_x[k] + g_z_0_xxyy_xx[k];

                g_z_0_xxxyy_y[k] = -g_z_0_xxyy_y[k] * cd_x[k] + g_z_0_xxyy_xy[k];

                g_z_0_xxxyy_z[k] = -g_z_0_xxyy_z[k] * cd_x[k] + g_z_0_xxyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 12);

            auto g_z_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 13);

            auto g_z_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_x, g_z_0_xxxyz_y, g_z_0_xxxyz_z, g_z_0_xxyz_x, g_z_0_xxyz_xx, g_z_0_xxyz_xy, g_z_0_xxyz_xz, g_z_0_xxyz_y, g_z_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_x[k] = -g_z_0_xxyz_x[k] * cd_x[k] + g_z_0_xxyz_xx[k];

                g_z_0_xxxyz_y[k] = -g_z_0_xxyz_y[k] * cd_x[k] + g_z_0_xxyz_xy[k];

                g_z_0_xxxyz_z[k] = -g_z_0_xxyz_z[k] * cd_x[k] + g_z_0_xxyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 15);

            auto g_z_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 16);

            auto g_z_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_x, g_z_0_xxxzz_y, g_z_0_xxxzz_z, g_z_0_xxzz_x, g_z_0_xxzz_xx, g_z_0_xxzz_xy, g_z_0_xxzz_xz, g_z_0_xxzz_y, g_z_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_x[k] = -g_z_0_xxzz_x[k] * cd_x[k] + g_z_0_xxzz_xx[k];

                g_z_0_xxxzz_y[k] = -g_z_0_xxzz_y[k] * cd_x[k] + g_z_0_xxzz_xy[k];

                g_z_0_xxxzz_z[k] = -g_z_0_xxzz_z[k] * cd_x[k] + g_z_0_xxzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 18);

            auto g_z_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 19);

            auto g_z_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_x, g_z_0_xxyyy_y, g_z_0_xxyyy_z, g_z_0_xyyy_x, g_z_0_xyyy_xx, g_z_0_xyyy_xy, g_z_0_xyyy_xz, g_z_0_xyyy_y, g_z_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_x[k] = -g_z_0_xyyy_x[k] * cd_x[k] + g_z_0_xyyy_xx[k];

                g_z_0_xxyyy_y[k] = -g_z_0_xyyy_y[k] * cd_x[k] + g_z_0_xyyy_xy[k];

                g_z_0_xxyyy_z[k] = -g_z_0_xyyy_z[k] * cd_x[k] + g_z_0_xyyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 21);

            auto g_z_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 22);

            auto g_z_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_x, g_z_0_xxyyz_y, g_z_0_xxyyz_z, g_z_0_xyyz_x, g_z_0_xyyz_xx, g_z_0_xyyz_xy, g_z_0_xyyz_xz, g_z_0_xyyz_y, g_z_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_x[k] = -g_z_0_xyyz_x[k] * cd_x[k] + g_z_0_xyyz_xx[k];

                g_z_0_xxyyz_y[k] = -g_z_0_xyyz_y[k] * cd_x[k] + g_z_0_xyyz_xy[k];

                g_z_0_xxyyz_z[k] = -g_z_0_xyyz_z[k] * cd_x[k] + g_z_0_xyyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 24);

            auto g_z_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 25);

            auto g_z_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_x, g_z_0_xxyzz_y, g_z_0_xxyzz_z, g_z_0_xyzz_x, g_z_0_xyzz_xx, g_z_0_xyzz_xy, g_z_0_xyzz_xz, g_z_0_xyzz_y, g_z_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_x[k] = -g_z_0_xyzz_x[k] * cd_x[k] + g_z_0_xyzz_xx[k];

                g_z_0_xxyzz_y[k] = -g_z_0_xyzz_y[k] * cd_x[k] + g_z_0_xyzz_xy[k];

                g_z_0_xxyzz_z[k] = -g_z_0_xyzz_z[k] * cd_x[k] + g_z_0_xyzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 27);

            auto g_z_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 28);

            auto g_z_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_x, g_z_0_xxzzz_y, g_z_0_xxzzz_z, g_z_0_xzzz_x, g_z_0_xzzz_xx, g_z_0_xzzz_xy, g_z_0_xzzz_xz, g_z_0_xzzz_y, g_z_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_x[k] = -g_z_0_xzzz_x[k] * cd_x[k] + g_z_0_xzzz_xx[k];

                g_z_0_xxzzz_y[k] = -g_z_0_xzzz_y[k] * cd_x[k] + g_z_0_xzzz_xy[k];

                g_z_0_xxzzz_z[k] = -g_z_0_xzzz_z[k] * cd_x[k] + g_z_0_xzzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 30);

            auto g_z_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 31);

            auto g_z_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_x, g_z_0_xyyyy_y, g_z_0_xyyyy_z, g_z_0_yyyy_x, g_z_0_yyyy_xx, g_z_0_yyyy_xy, g_z_0_yyyy_xz, g_z_0_yyyy_y, g_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_x[k] = -g_z_0_yyyy_x[k] * cd_x[k] + g_z_0_yyyy_xx[k];

                g_z_0_xyyyy_y[k] = -g_z_0_yyyy_y[k] * cd_x[k] + g_z_0_yyyy_xy[k];

                g_z_0_xyyyy_z[k] = -g_z_0_yyyy_z[k] * cd_x[k] + g_z_0_yyyy_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 33);

            auto g_z_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 34);

            auto g_z_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_x, g_z_0_xyyyz_y, g_z_0_xyyyz_z, g_z_0_yyyz_x, g_z_0_yyyz_xx, g_z_0_yyyz_xy, g_z_0_yyyz_xz, g_z_0_yyyz_y, g_z_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_x[k] = -g_z_0_yyyz_x[k] * cd_x[k] + g_z_0_yyyz_xx[k];

                g_z_0_xyyyz_y[k] = -g_z_0_yyyz_y[k] * cd_x[k] + g_z_0_yyyz_xy[k];

                g_z_0_xyyyz_z[k] = -g_z_0_yyyz_z[k] * cd_x[k] + g_z_0_yyyz_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 36);

            auto g_z_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 37);

            auto g_z_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_x, g_z_0_xyyzz_y, g_z_0_xyyzz_z, g_z_0_yyzz_x, g_z_0_yyzz_xx, g_z_0_yyzz_xy, g_z_0_yyzz_xz, g_z_0_yyzz_y, g_z_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_x[k] = -g_z_0_yyzz_x[k] * cd_x[k] + g_z_0_yyzz_xx[k];

                g_z_0_xyyzz_y[k] = -g_z_0_yyzz_y[k] * cd_x[k] + g_z_0_yyzz_xy[k];

                g_z_0_xyyzz_z[k] = -g_z_0_yyzz_z[k] * cd_x[k] + g_z_0_yyzz_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 39);

            auto g_z_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 40);

            auto g_z_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_x, g_z_0_xyzzz_y, g_z_0_xyzzz_z, g_z_0_yzzz_x, g_z_0_yzzz_xx, g_z_0_yzzz_xy, g_z_0_yzzz_xz, g_z_0_yzzz_y, g_z_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_x[k] = -g_z_0_yzzz_x[k] * cd_x[k] + g_z_0_yzzz_xx[k];

                g_z_0_xyzzz_y[k] = -g_z_0_yzzz_y[k] * cd_x[k] + g_z_0_yzzz_xy[k];

                g_z_0_xyzzz_z[k] = -g_z_0_yzzz_z[k] * cd_x[k] + g_z_0_yzzz_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 42);

            auto g_z_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 43);

            auto g_z_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_x, g_z_0_xzzzz_y, g_z_0_xzzzz_z, g_z_0_zzzz_x, g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_y, g_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_x[k] = -g_z_0_zzzz_x[k] * cd_x[k] + g_z_0_zzzz_xx[k];

                g_z_0_xzzzz_y[k] = -g_z_0_zzzz_y[k] * cd_x[k] + g_z_0_zzzz_xy[k];

                g_z_0_xzzzz_z[k] = -g_z_0_zzzz_z[k] * cd_x[k] + g_z_0_zzzz_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 45);

            auto g_z_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 46);

            auto g_z_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_x, g_z_0_yyyy_xy, g_z_0_yyyy_y, g_z_0_yyyy_yy, g_z_0_yyyy_yz, g_z_0_yyyy_z, g_z_0_yyyyy_x, g_z_0_yyyyy_y, g_z_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_x[k] = -g_z_0_yyyy_x[k] * cd_y[k] + g_z_0_yyyy_xy[k];

                g_z_0_yyyyy_y[k] = -g_z_0_yyyy_y[k] * cd_y[k] + g_z_0_yyyy_yy[k];

                g_z_0_yyyyy_z[k] = -g_z_0_yyyy_z[k] * cd_y[k] + g_z_0_yyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 48);

            auto g_z_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 49);

            auto g_z_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_x, g_z_0_yyyyz_y, g_z_0_yyyyz_z, g_z_0_yyyz_x, g_z_0_yyyz_xy, g_z_0_yyyz_y, g_z_0_yyyz_yy, g_z_0_yyyz_yz, g_z_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_x[k] = -g_z_0_yyyz_x[k] * cd_y[k] + g_z_0_yyyz_xy[k];

                g_z_0_yyyyz_y[k] = -g_z_0_yyyz_y[k] * cd_y[k] + g_z_0_yyyz_yy[k];

                g_z_0_yyyyz_z[k] = -g_z_0_yyyz_z[k] * cd_y[k] + g_z_0_yyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 51);

            auto g_z_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 52);

            auto g_z_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_x, g_z_0_yyyzz_y, g_z_0_yyyzz_z, g_z_0_yyzz_x, g_z_0_yyzz_xy, g_z_0_yyzz_y, g_z_0_yyzz_yy, g_z_0_yyzz_yz, g_z_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_x[k] = -g_z_0_yyzz_x[k] * cd_y[k] + g_z_0_yyzz_xy[k];

                g_z_0_yyyzz_y[k] = -g_z_0_yyzz_y[k] * cd_y[k] + g_z_0_yyzz_yy[k];

                g_z_0_yyyzz_z[k] = -g_z_0_yyzz_z[k] * cd_y[k] + g_z_0_yyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 54);

            auto g_z_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 55);

            auto g_z_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_x, g_z_0_yyzzz_y, g_z_0_yyzzz_z, g_z_0_yzzz_x, g_z_0_yzzz_xy, g_z_0_yzzz_y, g_z_0_yzzz_yy, g_z_0_yzzz_yz, g_z_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_x[k] = -g_z_0_yzzz_x[k] * cd_y[k] + g_z_0_yzzz_xy[k];

                g_z_0_yyzzz_y[k] = -g_z_0_yzzz_y[k] * cd_y[k] + g_z_0_yzzz_yy[k];

                g_z_0_yyzzz_z[k] = -g_z_0_yzzz_z[k] * cd_y[k] + g_z_0_yzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 57);

            auto g_z_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 58);

            auto g_z_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_x, g_z_0_yzzzz_y, g_z_0_yzzzz_z, g_z_0_zzzz_x, g_z_0_zzzz_xy, g_z_0_zzzz_y, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_x[k] = -g_z_0_zzzz_x[k] * cd_y[k] + g_z_0_zzzz_xy[k];

                g_z_0_yzzzz_y[k] = -g_z_0_zzzz_y[k] * cd_y[k] + g_z_0_zzzz_yy[k];

                g_z_0_yzzzz_z[k] = -g_z_0_zzzz_z[k] * cd_y[k] + g_z_0_zzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 60);

            auto g_z_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 61);

            auto g_z_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_x, g_z_0_zzzz_xz, g_z_0_zzzz_y, g_z_0_zzzz_yz, g_z_0_zzzz_z, g_z_0_zzzz_zz, g_z_0_zzzzz_x, g_z_0_zzzzz_y, g_z_0_zzzzz_z, g_zzzz_x, g_zzzz_y, g_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_x[k] = -g_zzzz_x[k] - g_z_0_zzzz_x[k] * cd_z[k] + g_z_0_zzzz_xz[k];

                g_z_0_zzzzz_y[k] = -g_zzzz_y[k] - g_z_0_zzzz_y[k] * cd_z[k] + g_z_0_zzzz_yz[k];

                g_z_0_zzzzz_z[k] = -g_zzzz_z[k] - g_z_0_zzzz_z[k] * cd_z[k] + g_z_0_zzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

