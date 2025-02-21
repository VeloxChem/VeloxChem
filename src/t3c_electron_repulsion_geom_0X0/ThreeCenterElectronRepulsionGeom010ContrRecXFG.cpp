#include "ThreeCenterElectronRepulsionGeom010ContrRecXFG.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xfg(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfg,
                                        const size_t idx_xdg,
                                        const size_t idx_geom_10_xdg,
                                        const size_t idx_geom_10_xdh,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SDG

        const auto dg_off = idx_xdg + i * 90;

        auto g_xx_xxxx = cbuffer.data(dg_off + 0);

        auto g_xx_xxxy = cbuffer.data(dg_off + 1);

        auto g_xx_xxxz = cbuffer.data(dg_off + 2);

        auto g_xx_xxyy = cbuffer.data(dg_off + 3);

        auto g_xx_xxyz = cbuffer.data(dg_off + 4);

        auto g_xx_xxzz = cbuffer.data(dg_off + 5);

        auto g_xx_xyyy = cbuffer.data(dg_off + 6);

        auto g_xx_xyyz = cbuffer.data(dg_off + 7);

        auto g_xx_xyzz = cbuffer.data(dg_off + 8);

        auto g_xx_xzzz = cbuffer.data(dg_off + 9);

        auto g_xx_yyyy = cbuffer.data(dg_off + 10);

        auto g_xx_yyyz = cbuffer.data(dg_off + 11);

        auto g_xx_yyzz = cbuffer.data(dg_off + 12);

        auto g_xx_yzzz = cbuffer.data(dg_off + 13);

        auto g_xx_zzzz = cbuffer.data(dg_off + 14);

        auto g_yy_xxxx = cbuffer.data(dg_off + 45);

        auto g_yy_xxxy = cbuffer.data(dg_off + 46);

        auto g_yy_xxxz = cbuffer.data(dg_off + 47);

        auto g_yy_xxyy = cbuffer.data(dg_off + 48);

        auto g_yy_xxyz = cbuffer.data(dg_off + 49);

        auto g_yy_xxzz = cbuffer.data(dg_off + 50);

        auto g_yy_xyyy = cbuffer.data(dg_off + 51);

        auto g_yy_xyyz = cbuffer.data(dg_off + 52);

        auto g_yy_xyzz = cbuffer.data(dg_off + 53);

        auto g_yy_xzzz = cbuffer.data(dg_off + 54);

        auto g_yy_yyyy = cbuffer.data(dg_off + 55);

        auto g_yy_yyyz = cbuffer.data(dg_off + 56);

        auto g_yy_yyzz = cbuffer.data(dg_off + 57);

        auto g_yy_yzzz = cbuffer.data(dg_off + 58);

        auto g_yy_zzzz = cbuffer.data(dg_off + 59);

        auto g_zz_xxxx = cbuffer.data(dg_off + 75);

        auto g_zz_xxxy = cbuffer.data(dg_off + 76);

        auto g_zz_xxxz = cbuffer.data(dg_off + 77);

        auto g_zz_xxyy = cbuffer.data(dg_off + 78);

        auto g_zz_xxyz = cbuffer.data(dg_off + 79);

        auto g_zz_xxzz = cbuffer.data(dg_off + 80);

        auto g_zz_xyyy = cbuffer.data(dg_off + 81);

        auto g_zz_xyyz = cbuffer.data(dg_off + 82);

        auto g_zz_xyzz = cbuffer.data(dg_off + 83);

        auto g_zz_xzzz = cbuffer.data(dg_off + 84);

        auto g_zz_yyyy = cbuffer.data(dg_off + 85);

        auto g_zz_yyyz = cbuffer.data(dg_off + 86);

        auto g_zz_yyzz = cbuffer.data(dg_off + 87);

        auto g_zz_yzzz = cbuffer.data(dg_off + 88);

        auto g_zz_zzzz = cbuffer.data(dg_off + 89);

        /// Set up components of auxilary buffer : SDG

        const auto dg_geom_10_off = idx_geom_10_xdg + i * 90;

        auto g_x_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 30);

        auto g_x_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 35);

        auto g_x_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 41);

        auto g_x_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 42);

        auto g_x_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 44);

        auto g_x_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 45);

        auto g_x_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 46);

        auto g_x_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 47);

        auto g_x_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 48);

        auto g_x_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 49);

        auto g_x_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 50);

        auto g_x_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 51);

        auto g_x_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 52);

        auto g_x_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 53);

        auto g_x_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 54);

        auto g_x_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 55);

        auto g_x_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 56);

        auto g_x_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 57);

        auto g_x_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 58);

        auto g_x_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 59);

        auto g_x_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 60);

        auto g_x_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 61);

        auto g_x_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 62);

        auto g_x_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 63);

        auto g_x_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 64);

        auto g_x_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 65);

        auto g_x_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 66);

        auto g_x_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 67);

        auto g_x_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 68);

        auto g_x_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 69);

        auto g_x_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 70);

        auto g_x_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 71);

        auto g_x_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 72);

        auto g_x_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 73);

        auto g_x_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 74);

        auto g_x_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps + 75);

        auto g_x_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps + 76);

        auto g_x_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps + 77);

        auto g_x_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 78);

        auto g_x_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 79);

        auto g_x_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 80);

        auto g_x_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 81);

        auto g_x_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 82);

        auto g_x_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 83);

        auto g_x_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 84);

        auto g_x_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps + 85);

        auto g_x_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps + 86);

        auto g_x_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 87);

        auto g_x_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 88);

        auto g_x_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps + 89);

        auto g_y_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 0);

        auto g_y_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 1);

        auto g_y_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 2);

        auto g_y_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 3);

        auto g_y_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 4);

        auto g_y_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 5);

        auto g_y_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 6);

        auto g_y_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 7);

        auto g_y_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 8);

        auto g_y_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 9);

        auto g_y_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 10);

        auto g_y_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 11);

        auto g_y_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 12);

        auto g_y_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 13);

        auto g_y_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 14);

        auto g_y_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 15);

        auto g_y_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 16);

        auto g_y_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 17);

        auto g_y_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 18);

        auto g_y_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 19);

        auto g_y_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 20);

        auto g_y_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 21);

        auto g_y_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 22);

        auto g_y_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 23);

        auto g_y_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 24);

        auto g_y_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 25);

        auto g_y_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 26);

        auto g_y_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 27);

        auto g_y_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 28);

        auto g_y_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 29);

        auto g_y_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 30);

        auto g_y_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 31);

        auto g_y_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 32);

        auto g_y_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 33);

        auto g_y_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 34);

        auto g_y_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 35);

        auto g_y_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 36);

        auto g_y_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 37);

        auto g_y_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 38);

        auto g_y_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 39);

        auto g_y_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 40);

        auto g_y_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 41);

        auto g_y_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 42);

        auto g_y_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 43);

        auto g_y_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 44);

        auto g_y_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 45);

        auto g_y_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 46);

        auto g_y_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 47);

        auto g_y_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 48);

        auto g_y_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 49);

        auto g_y_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 50);

        auto g_y_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 51);

        auto g_y_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 52);

        auto g_y_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 53);

        auto g_y_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 54);

        auto g_y_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 55);

        auto g_y_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 56);

        auto g_y_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 57);

        auto g_y_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 58);

        auto g_y_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 59);

        auto g_y_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 60);

        auto g_y_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 61);

        auto g_y_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 62);

        auto g_y_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 63);

        auto g_y_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 64);

        auto g_y_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 65);

        auto g_y_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 66);

        auto g_y_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 67);

        auto g_y_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 68);

        auto g_y_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 69);

        auto g_y_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 70);

        auto g_y_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 71);

        auto g_y_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 72);

        auto g_y_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 73);

        auto g_y_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 74);

        auto g_y_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps + 75);

        auto g_y_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps + 76);

        auto g_y_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps + 77);

        auto g_y_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 78);

        auto g_y_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 79);

        auto g_y_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 80);

        auto g_y_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 81);

        auto g_y_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 82);

        auto g_y_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 83);

        auto g_y_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 84);

        auto g_y_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps + 85);

        auto g_y_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps + 86);

        auto g_y_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 87);

        auto g_y_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 88);

        auto g_y_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps + 89);

        auto g_z_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 0);

        auto g_z_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 1);

        auto g_z_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 2);

        auto g_z_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 3);

        auto g_z_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 4);

        auto g_z_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 5);

        auto g_z_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 6);

        auto g_z_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 7);

        auto g_z_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 8);

        auto g_z_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 9);

        auto g_z_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 10);

        auto g_z_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 11);

        auto g_z_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 12);

        auto g_z_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 13);

        auto g_z_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 14);

        auto g_z_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 15);

        auto g_z_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 16);

        auto g_z_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 17);

        auto g_z_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 18);

        auto g_z_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 19);

        auto g_z_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 20);

        auto g_z_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 21);

        auto g_z_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 22);

        auto g_z_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 23);

        auto g_z_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 24);

        auto g_z_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 25);

        auto g_z_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 26);

        auto g_z_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 27);

        auto g_z_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 28);

        auto g_z_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 29);

        auto g_z_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 30);

        auto g_z_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 31);

        auto g_z_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 32);

        auto g_z_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 33);

        auto g_z_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 34);

        auto g_z_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 35);

        auto g_z_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 36);

        auto g_z_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 37);

        auto g_z_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 38);

        auto g_z_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 39);

        auto g_z_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 40);

        auto g_z_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 41);

        auto g_z_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 42);

        auto g_z_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 43);

        auto g_z_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 44);

        auto g_z_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 45);

        auto g_z_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 46);

        auto g_z_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 47);

        auto g_z_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 48);

        auto g_z_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 49);

        auto g_z_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 50);

        auto g_z_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 51);

        auto g_z_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 52);

        auto g_z_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 53);

        auto g_z_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 54);

        auto g_z_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 55);

        auto g_z_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 56);

        auto g_z_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 57);

        auto g_z_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 58);

        auto g_z_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 59);

        auto g_z_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 60);

        auto g_z_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 61);

        auto g_z_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 62);

        auto g_z_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 63);

        auto g_z_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 64);

        auto g_z_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 65);

        auto g_z_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 66);

        auto g_z_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 67);

        auto g_z_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 68);

        auto g_z_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 69);

        auto g_z_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 70);

        auto g_z_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 71);

        auto g_z_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 72);

        auto g_z_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 73);

        auto g_z_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 74);

        auto g_z_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps + 75);

        auto g_z_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps + 76);

        auto g_z_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps + 77);

        auto g_z_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 78);

        auto g_z_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 79);

        auto g_z_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 80);

        auto g_z_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 81);

        auto g_z_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 82);

        auto g_z_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 83);

        auto g_z_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 84);

        auto g_z_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps + 85);

        auto g_z_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps + 86);

        auto g_z_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 87);

        auto g_z_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 88);

        auto g_z_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps + 89);

        /// Set up components of auxilary buffer : SDH

        const auto dh_geom_10_off = idx_geom_10_xdh + i * 126;

        auto g_x_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 48);

        auto g_x_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 55);

        auto g_x_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 56);

        auto g_x_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 59);

        auto g_x_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 60);

        auto g_x_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 61);

        auto g_x_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 62);

        auto g_x_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 64);

        auto g_x_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 66);

        auto g_x_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 67);

        auto g_x_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 69);

        auto g_x_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 70);

        auto g_x_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 71);

        auto g_x_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 73);

        auto g_x_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 74);

        auto g_x_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 75);

        auto g_x_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 76);

        auto g_x_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 78);

        auto g_x_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 79);

        auto g_x_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 80);

        auto g_x_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 81);

        auto g_x_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 82);

        auto g_x_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 85);

        auto g_x_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 87);

        auto g_x_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 88);

        auto g_x_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 90);

        auto g_x_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 91);

        auto g_x_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 92);

        auto g_x_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 94);

        auto g_x_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 95);

        auto g_x_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 96);

        auto g_x_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 97);

        auto g_x_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 99);

        auto g_x_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 100);

        auto g_x_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 101);

        auto g_x_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 102);

        auto g_x_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 103);

        auto g_x_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 106);

        auto g_x_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 107);

        auto g_x_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 108);

        auto g_x_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 109);

        auto g_x_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 110);

        auto g_x_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 111);

        auto g_x_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 112);

        auto g_x_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 113);

        auto g_x_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 114);

        auto g_x_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 115);

        auto g_x_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 116);

        auto g_x_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 117);

        auto g_x_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 118);

        auto g_x_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 119);

        auto g_x_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 120);

        auto g_x_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 121);

        auto g_x_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 122);

        auto g_x_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 123);

        auto g_x_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 124);

        auto g_x_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 125);

        auto g_y_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 0);

        auto g_y_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 1);

        auto g_y_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 2);

        auto g_y_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 3);

        auto g_y_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 4);

        auto g_y_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 5);

        auto g_y_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 6);

        auto g_y_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 7);

        auto g_y_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 8);

        auto g_y_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 9);

        auto g_y_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 10);

        auto g_y_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 11);

        auto g_y_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 12);

        auto g_y_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 13);

        auto g_y_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 14);

        auto g_y_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 21);

        auto g_y_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 22);

        auto g_y_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 23);

        auto g_y_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 24);

        auto g_y_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 25);

        auto g_y_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 26);

        auto g_y_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 27);

        auto g_y_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 28);

        auto g_y_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 29);

        auto g_y_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 30);

        auto g_y_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 31);

        auto g_y_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 32);

        auto g_y_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 33);

        auto g_y_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 34);

        auto g_y_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 35);

        auto g_y_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 42);

        auto g_y_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 43);

        auto g_y_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 44);

        auto g_y_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 45);

        auto g_y_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 46);

        auto g_y_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 47);

        auto g_y_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 48);

        auto g_y_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 49);

        auto g_y_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 50);

        auto g_y_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 51);

        auto g_y_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 52);

        auto g_y_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 53);

        auto g_y_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 54);

        auto g_y_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 55);

        auto g_y_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 56);

        auto g_y_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 63);

        auto g_y_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 64);

        auto g_y_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 65);

        auto g_y_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 66);

        auto g_y_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 67);

        auto g_y_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 68);

        auto g_y_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 69);

        auto g_y_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 70);

        auto g_y_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 71);

        auto g_y_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 72);

        auto g_y_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 73);

        auto g_y_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 74);

        auto g_y_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 75);

        auto g_y_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 76);

        auto g_y_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 77);

        auto g_y_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 78);

        auto g_y_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 79);

        auto g_y_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 80);

        auto g_y_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 81);

        auto g_y_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 82);

        auto g_y_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 83);

        auto g_y_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 84);

        auto g_y_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 85);

        auto g_y_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 86);

        auto g_y_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 87);

        auto g_y_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 88);

        auto g_y_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 89);

        auto g_y_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 90);

        auto g_y_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 91);

        auto g_y_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 92);

        auto g_y_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 93);

        auto g_y_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 94);

        auto g_y_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 95);

        auto g_y_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 96);

        auto g_y_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 97);

        auto g_y_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 98);

        auto g_y_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 100);

        auto g_y_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 101);

        auto g_y_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 102);

        auto g_y_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 103);

        auto g_y_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 104);

        auto g_y_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 105);

        auto g_y_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 106);

        auto g_y_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 107);

        auto g_y_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 108);

        auto g_y_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 109);

        auto g_y_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 110);

        auto g_y_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 111);

        auto g_y_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 112);

        auto g_y_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 113);

        auto g_y_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 114);

        auto g_y_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 115);

        auto g_y_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 116);

        auto g_y_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 117);

        auto g_y_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 118);

        auto g_y_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 119);

        auto g_y_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 121);

        auto g_y_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 122);

        auto g_y_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 123);

        auto g_y_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 124);

        auto g_y_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 125);

        auto g_z_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 0);

        auto g_z_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 1);

        auto g_z_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 2);

        auto g_z_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 3);

        auto g_z_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 4);

        auto g_z_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 5);

        auto g_z_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 6);

        auto g_z_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 7);

        auto g_z_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 8);

        auto g_z_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 9);

        auto g_z_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 10);

        auto g_z_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 11);

        auto g_z_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 12);

        auto g_z_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 13);

        auto g_z_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 14);

        auto g_z_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 21);

        auto g_z_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 22);

        auto g_z_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 23);

        auto g_z_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 24);

        auto g_z_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 25);

        auto g_z_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 26);

        auto g_z_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 27);

        auto g_z_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 28);

        auto g_z_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 29);

        auto g_z_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 30);

        auto g_z_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 31);

        auto g_z_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 32);

        auto g_z_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 33);

        auto g_z_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 34);

        auto g_z_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 35);

        auto g_z_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 42);

        auto g_z_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 43);

        auto g_z_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 44);

        auto g_z_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 45);

        auto g_z_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 46);

        auto g_z_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 47);

        auto g_z_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 48);

        auto g_z_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 49);

        auto g_z_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 50);

        auto g_z_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 51);

        auto g_z_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 52);

        auto g_z_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 53);

        auto g_z_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 54);

        auto g_z_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 55);

        auto g_z_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 56);

        auto g_z_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 63);

        auto g_z_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 64);

        auto g_z_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 65);

        auto g_z_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 66);

        auto g_z_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 67);

        auto g_z_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 68);

        auto g_z_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 69);

        auto g_z_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 70);

        auto g_z_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 71);

        auto g_z_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 72);

        auto g_z_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 73);

        auto g_z_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 74);

        auto g_z_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 75);

        auto g_z_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 76);

        auto g_z_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 77);

        auto g_z_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 78);

        auto g_z_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 79);

        auto g_z_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 80);

        auto g_z_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 81);

        auto g_z_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 82);

        auto g_z_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 84);

        auto g_z_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 85);

        auto g_z_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 86);

        auto g_z_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 87);

        auto g_z_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 88);

        auto g_z_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 89);

        auto g_z_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 90);

        auto g_z_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 91);

        auto g_z_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 92);

        auto g_z_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 93);

        auto g_z_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 94);

        auto g_z_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 95);

        auto g_z_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 96);

        auto g_z_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 97);

        auto g_z_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 98);

        auto g_z_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 99);

        auto g_z_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 100);

        auto g_z_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 101);

        auto g_z_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 102);

        auto g_z_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 103);

        auto g_z_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 105);

        auto g_z_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 106);

        auto g_z_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 107);

        auto g_z_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 108);

        auto g_z_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 109);

        auto g_z_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 110);

        auto g_z_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 111);

        auto g_z_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 112);

        auto g_z_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 113);

        auto g_z_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 114);

        auto g_z_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 115);

        auto g_z_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 116);

        auto g_z_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 117);

        auto g_z_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 118);

        auto g_z_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 119);

        auto g_z_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 120);

        auto g_z_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 121);

        auto g_z_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 122);

        auto g_z_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 123);

        auto g_z_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 124);

        auto g_z_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 125);

        /// set up bra offset for contr_buffer_xxfg

        const auto fg_geom_10_off = idx_geom_10_xfg + i * 150;

        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_x_0_xx_xxxx, g_x_0_xx_xxxxx, g_x_0_xx_xxxxy, g_x_0_xx_xxxxz, g_x_0_xx_xxxy, g_x_0_xx_xxxyy, g_x_0_xx_xxxyz, g_x_0_xx_xxxz, g_x_0_xx_xxxzz, g_x_0_xx_xxyy, g_x_0_xx_xxyyy, g_x_0_xx_xxyyz, g_x_0_xx_xxyz, g_x_0_xx_xxyzz, g_x_0_xx_xxzz, g_x_0_xx_xxzzz, g_x_0_xx_xyyy, g_x_0_xx_xyyyy, g_x_0_xx_xyyyz, g_x_0_xx_xyyz, g_x_0_xx_xyyzz, g_x_0_xx_xyzz, g_x_0_xx_xyzzz, g_x_0_xx_xzzz, g_x_0_xx_xzzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyzz, g_x_0_xx_yzzz, g_x_0_xx_zzzz, g_x_0_xxx_xxxx, g_x_0_xxx_xxxy, g_x_0_xxx_xxxz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyz, g_x_0_xxx_xxzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyz, g_x_0_xxx_xyzz, g_x_0_xxx_xzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyz, g_x_0_xxx_yyzz, g_x_0_xxx_yzzz, g_x_0_xxx_zzzz, g_xx_xxxx, g_xx_xxxy, g_xx_xxxz, g_xx_xxyy, g_xx_xxyz, g_xx_xxzz, g_xx_xyyy, g_xx_xyyz, g_xx_xyzz, g_xx_xzzz, g_xx_yyyy, g_xx_yyyz, g_xx_yyzz, g_xx_yzzz, g_xx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_xxxx[k] = -g_xx_xxxx[k] - g_x_0_xx_xxxx[k] * cd_x[k] + g_x_0_xx_xxxxx[k];

            g_x_0_xxx_xxxy[k] = -g_xx_xxxy[k] - g_x_0_xx_xxxy[k] * cd_x[k] + g_x_0_xx_xxxxy[k];

            g_x_0_xxx_xxxz[k] = -g_xx_xxxz[k] - g_x_0_xx_xxxz[k] * cd_x[k] + g_x_0_xx_xxxxz[k];

            g_x_0_xxx_xxyy[k] = -g_xx_xxyy[k] - g_x_0_xx_xxyy[k] * cd_x[k] + g_x_0_xx_xxxyy[k];

            g_x_0_xxx_xxyz[k] = -g_xx_xxyz[k] - g_x_0_xx_xxyz[k] * cd_x[k] + g_x_0_xx_xxxyz[k];

            g_x_0_xxx_xxzz[k] = -g_xx_xxzz[k] - g_x_0_xx_xxzz[k] * cd_x[k] + g_x_0_xx_xxxzz[k];

            g_x_0_xxx_xyyy[k] = -g_xx_xyyy[k] - g_x_0_xx_xyyy[k] * cd_x[k] + g_x_0_xx_xxyyy[k];

            g_x_0_xxx_xyyz[k] = -g_xx_xyyz[k] - g_x_0_xx_xyyz[k] * cd_x[k] + g_x_0_xx_xxyyz[k];

            g_x_0_xxx_xyzz[k] = -g_xx_xyzz[k] - g_x_0_xx_xyzz[k] * cd_x[k] + g_x_0_xx_xxyzz[k];

            g_x_0_xxx_xzzz[k] = -g_xx_xzzz[k] - g_x_0_xx_xzzz[k] * cd_x[k] + g_x_0_xx_xxzzz[k];

            g_x_0_xxx_yyyy[k] = -g_xx_yyyy[k] - g_x_0_xx_yyyy[k] * cd_x[k] + g_x_0_xx_xyyyy[k];

            g_x_0_xxx_yyyz[k] = -g_xx_yyyz[k] - g_x_0_xx_yyyz[k] * cd_x[k] + g_x_0_xx_xyyyz[k];

            g_x_0_xxx_yyzz[k] = -g_xx_yyzz[k] - g_x_0_xx_yyzz[k] * cd_x[k] + g_x_0_xx_xyyzz[k];

            g_x_0_xxx_yzzz[k] = -g_xx_yzzz[k] - g_x_0_xx_yzzz[k] * cd_x[k] + g_x_0_xx_xyzzz[k];

            g_x_0_xxx_zzzz[k] = -g_xx_zzzz[k] - g_x_0_xx_zzzz[k] * cd_x[k] + g_x_0_xx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_xx_xxxx, g_x_0_xx_xxxxy, g_x_0_xx_xxxy, g_x_0_xx_xxxyy, g_x_0_xx_xxxyz, g_x_0_xx_xxxz, g_x_0_xx_xxyy, g_x_0_xx_xxyyy, g_x_0_xx_xxyyz, g_x_0_xx_xxyz, g_x_0_xx_xxyzz, g_x_0_xx_xxzz, g_x_0_xx_xyyy, g_x_0_xx_xyyyy, g_x_0_xx_xyyyz, g_x_0_xx_xyyz, g_x_0_xx_xyyzz, g_x_0_xx_xyzz, g_x_0_xx_xyzzz, g_x_0_xx_xzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyyy, g_x_0_xx_yyyyz, g_x_0_xx_yyyz, g_x_0_xx_yyyzz, g_x_0_xx_yyzz, g_x_0_xx_yyzzz, g_x_0_xx_yzzz, g_x_0_xx_yzzzz, g_x_0_xx_zzzz, g_x_0_xxy_xxxx, g_x_0_xxy_xxxy, g_x_0_xxy_xxxz, g_x_0_xxy_xxyy, g_x_0_xxy_xxyz, g_x_0_xxy_xxzz, g_x_0_xxy_xyyy, g_x_0_xxy_xyyz, g_x_0_xxy_xyzz, g_x_0_xxy_xzzz, g_x_0_xxy_yyyy, g_x_0_xxy_yyyz, g_x_0_xxy_yyzz, g_x_0_xxy_yzzz, g_x_0_xxy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_xxxx[k] = -g_x_0_xx_xxxx[k] * cd_y[k] + g_x_0_xx_xxxxy[k];

            g_x_0_xxy_xxxy[k] = -g_x_0_xx_xxxy[k] * cd_y[k] + g_x_0_xx_xxxyy[k];

            g_x_0_xxy_xxxz[k] = -g_x_0_xx_xxxz[k] * cd_y[k] + g_x_0_xx_xxxyz[k];

            g_x_0_xxy_xxyy[k] = -g_x_0_xx_xxyy[k] * cd_y[k] + g_x_0_xx_xxyyy[k];

            g_x_0_xxy_xxyz[k] = -g_x_0_xx_xxyz[k] * cd_y[k] + g_x_0_xx_xxyyz[k];

            g_x_0_xxy_xxzz[k] = -g_x_0_xx_xxzz[k] * cd_y[k] + g_x_0_xx_xxyzz[k];

            g_x_0_xxy_xyyy[k] = -g_x_0_xx_xyyy[k] * cd_y[k] + g_x_0_xx_xyyyy[k];

            g_x_0_xxy_xyyz[k] = -g_x_0_xx_xyyz[k] * cd_y[k] + g_x_0_xx_xyyyz[k];

            g_x_0_xxy_xyzz[k] = -g_x_0_xx_xyzz[k] * cd_y[k] + g_x_0_xx_xyyzz[k];

            g_x_0_xxy_xzzz[k] = -g_x_0_xx_xzzz[k] * cd_y[k] + g_x_0_xx_xyzzz[k];

            g_x_0_xxy_yyyy[k] = -g_x_0_xx_yyyy[k] * cd_y[k] + g_x_0_xx_yyyyy[k];

            g_x_0_xxy_yyyz[k] = -g_x_0_xx_yyyz[k] * cd_y[k] + g_x_0_xx_yyyyz[k];

            g_x_0_xxy_yyzz[k] = -g_x_0_xx_yyzz[k] * cd_y[k] + g_x_0_xx_yyyzz[k];

            g_x_0_xxy_yzzz[k] = -g_x_0_xx_yzzz[k] * cd_y[k] + g_x_0_xx_yyzzz[k];

            g_x_0_xxy_zzzz[k] = -g_x_0_xx_zzzz[k] * cd_y[k] + g_x_0_xx_yzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 44);

        #pragma omp simd aligned(cd_z, g_x_0_xx_xxxx, g_x_0_xx_xxxxz, g_x_0_xx_xxxy, g_x_0_xx_xxxyz, g_x_0_xx_xxxz, g_x_0_xx_xxxzz, g_x_0_xx_xxyy, g_x_0_xx_xxyyz, g_x_0_xx_xxyz, g_x_0_xx_xxyzz, g_x_0_xx_xxzz, g_x_0_xx_xxzzz, g_x_0_xx_xyyy, g_x_0_xx_xyyyz, g_x_0_xx_xyyz, g_x_0_xx_xyyzz, g_x_0_xx_xyzz, g_x_0_xx_xyzzz, g_x_0_xx_xzzz, g_x_0_xx_xzzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyyz, g_x_0_xx_yyyz, g_x_0_xx_yyyzz, g_x_0_xx_yyzz, g_x_0_xx_yyzzz, g_x_0_xx_yzzz, g_x_0_xx_yzzzz, g_x_0_xx_zzzz, g_x_0_xx_zzzzz, g_x_0_xxz_xxxx, g_x_0_xxz_xxxy, g_x_0_xxz_xxxz, g_x_0_xxz_xxyy, g_x_0_xxz_xxyz, g_x_0_xxz_xxzz, g_x_0_xxz_xyyy, g_x_0_xxz_xyyz, g_x_0_xxz_xyzz, g_x_0_xxz_xzzz, g_x_0_xxz_yyyy, g_x_0_xxz_yyyz, g_x_0_xxz_yyzz, g_x_0_xxz_yzzz, g_x_0_xxz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_xxxx[k] = -g_x_0_xx_xxxx[k] * cd_z[k] + g_x_0_xx_xxxxz[k];

            g_x_0_xxz_xxxy[k] = -g_x_0_xx_xxxy[k] * cd_z[k] + g_x_0_xx_xxxyz[k];

            g_x_0_xxz_xxxz[k] = -g_x_0_xx_xxxz[k] * cd_z[k] + g_x_0_xx_xxxzz[k];

            g_x_0_xxz_xxyy[k] = -g_x_0_xx_xxyy[k] * cd_z[k] + g_x_0_xx_xxyyz[k];

            g_x_0_xxz_xxyz[k] = -g_x_0_xx_xxyz[k] * cd_z[k] + g_x_0_xx_xxyzz[k];

            g_x_0_xxz_xxzz[k] = -g_x_0_xx_xxzz[k] * cd_z[k] + g_x_0_xx_xxzzz[k];

            g_x_0_xxz_xyyy[k] = -g_x_0_xx_xyyy[k] * cd_z[k] + g_x_0_xx_xyyyz[k];

            g_x_0_xxz_xyyz[k] = -g_x_0_xx_xyyz[k] * cd_z[k] + g_x_0_xx_xyyzz[k];

            g_x_0_xxz_xyzz[k] = -g_x_0_xx_xyzz[k] * cd_z[k] + g_x_0_xx_xyzzz[k];

            g_x_0_xxz_xzzz[k] = -g_x_0_xx_xzzz[k] * cd_z[k] + g_x_0_xx_xzzzz[k];

            g_x_0_xxz_yyyy[k] = -g_x_0_xx_yyyy[k] * cd_z[k] + g_x_0_xx_yyyyz[k];

            g_x_0_xxz_yyyz[k] = -g_x_0_xx_yyyz[k] * cd_z[k] + g_x_0_xx_yyyzz[k];

            g_x_0_xxz_yyzz[k] = -g_x_0_xx_yyzz[k] * cd_z[k] + g_x_0_xx_yyzzz[k];

            g_x_0_xxz_yzzz[k] = -g_x_0_xx_yzzz[k] * cd_z[k] + g_x_0_xx_yzzzz[k];

            g_x_0_xxz_zzzz[k] = -g_x_0_xx_zzzz[k] * cd_z[k] + g_x_0_xx_zzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_y, g_x_0_xy_xxxx, g_x_0_xy_xxxxy, g_x_0_xy_xxxy, g_x_0_xy_xxxyy, g_x_0_xy_xxxyz, g_x_0_xy_xxxz, g_x_0_xy_xxyy, g_x_0_xy_xxyyy, g_x_0_xy_xxyyz, g_x_0_xy_xxyz, g_x_0_xy_xxyzz, g_x_0_xy_xxzz, g_x_0_xy_xyyy, g_x_0_xy_xyyyy, g_x_0_xy_xyyyz, g_x_0_xy_xyyz, g_x_0_xy_xyyzz, g_x_0_xy_xyzz, g_x_0_xy_xyzzz, g_x_0_xy_xzzz, g_x_0_xy_yyyy, g_x_0_xy_yyyyy, g_x_0_xy_yyyyz, g_x_0_xy_yyyz, g_x_0_xy_yyyzz, g_x_0_xy_yyzz, g_x_0_xy_yyzzz, g_x_0_xy_yzzz, g_x_0_xy_yzzzz, g_x_0_xy_zzzz, g_x_0_xyy_xxxx, g_x_0_xyy_xxxy, g_x_0_xyy_xxxz, g_x_0_xyy_xxyy, g_x_0_xyy_xxyz, g_x_0_xyy_xxzz, g_x_0_xyy_xyyy, g_x_0_xyy_xyyz, g_x_0_xyy_xyzz, g_x_0_xyy_xzzz, g_x_0_xyy_yyyy, g_x_0_xyy_yyyz, g_x_0_xyy_yyzz, g_x_0_xyy_yzzz, g_x_0_xyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_xxxx[k] = -g_x_0_xy_xxxx[k] * cd_y[k] + g_x_0_xy_xxxxy[k];

            g_x_0_xyy_xxxy[k] = -g_x_0_xy_xxxy[k] * cd_y[k] + g_x_0_xy_xxxyy[k];

            g_x_0_xyy_xxxz[k] = -g_x_0_xy_xxxz[k] * cd_y[k] + g_x_0_xy_xxxyz[k];

            g_x_0_xyy_xxyy[k] = -g_x_0_xy_xxyy[k] * cd_y[k] + g_x_0_xy_xxyyy[k];

            g_x_0_xyy_xxyz[k] = -g_x_0_xy_xxyz[k] * cd_y[k] + g_x_0_xy_xxyyz[k];

            g_x_0_xyy_xxzz[k] = -g_x_0_xy_xxzz[k] * cd_y[k] + g_x_0_xy_xxyzz[k];

            g_x_0_xyy_xyyy[k] = -g_x_0_xy_xyyy[k] * cd_y[k] + g_x_0_xy_xyyyy[k];

            g_x_0_xyy_xyyz[k] = -g_x_0_xy_xyyz[k] * cd_y[k] + g_x_0_xy_xyyyz[k];

            g_x_0_xyy_xyzz[k] = -g_x_0_xy_xyzz[k] * cd_y[k] + g_x_0_xy_xyyzz[k];

            g_x_0_xyy_xzzz[k] = -g_x_0_xy_xzzz[k] * cd_y[k] + g_x_0_xy_xyzzz[k];

            g_x_0_xyy_yyyy[k] = -g_x_0_xy_yyyy[k] * cd_y[k] + g_x_0_xy_yyyyy[k];

            g_x_0_xyy_yyyz[k] = -g_x_0_xy_yyyz[k] * cd_y[k] + g_x_0_xy_yyyyz[k];

            g_x_0_xyy_yyzz[k] = -g_x_0_xy_yyzz[k] * cd_y[k] + g_x_0_xy_yyyzz[k];

            g_x_0_xyy_yzzz[k] = -g_x_0_xy_yzzz[k] * cd_y[k] + g_x_0_xy_yyzzz[k];

            g_x_0_xyy_zzzz[k] = -g_x_0_xy_zzzz[k] * cd_y[k] + g_x_0_xy_yzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 74);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_xxxx, g_x_0_xyz_xxxy, g_x_0_xyz_xxxz, g_x_0_xyz_xxyy, g_x_0_xyz_xxyz, g_x_0_xyz_xxzz, g_x_0_xyz_xyyy, g_x_0_xyz_xyyz, g_x_0_xyz_xyzz, g_x_0_xyz_xzzz, g_x_0_xyz_yyyy, g_x_0_xyz_yyyz, g_x_0_xyz_yyzz, g_x_0_xyz_yzzz, g_x_0_xyz_zzzz, g_x_0_xz_xxxx, g_x_0_xz_xxxxy, g_x_0_xz_xxxy, g_x_0_xz_xxxyy, g_x_0_xz_xxxyz, g_x_0_xz_xxxz, g_x_0_xz_xxyy, g_x_0_xz_xxyyy, g_x_0_xz_xxyyz, g_x_0_xz_xxyz, g_x_0_xz_xxyzz, g_x_0_xz_xxzz, g_x_0_xz_xyyy, g_x_0_xz_xyyyy, g_x_0_xz_xyyyz, g_x_0_xz_xyyz, g_x_0_xz_xyyzz, g_x_0_xz_xyzz, g_x_0_xz_xyzzz, g_x_0_xz_xzzz, g_x_0_xz_yyyy, g_x_0_xz_yyyyy, g_x_0_xz_yyyyz, g_x_0_xz_yyyz, g_x_0_xz_yyyzz, g_x_0_xz_yyzz, g_x_0_xz_yyzzz, g_x_0_xz_yzzz, g_x_0_xz_yzzzz, g_x_0_xz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_xxxx[k] = -g_x_0_xz_xxxx[k] * cd_y[k] + g_x_0_xz_xxxxy[k];

            g_x_0_xyz_xxxy[k] = -g_x_0_xz_xxxy[k] * cd_y[k] + g_x_0_xz_xxxyy[k];

            g_x_0_xyz_xxxz[k] = -g_x_0_xz_xxxz[k] * cd_y[k] + g_x_0_xz_xxxyz[k];

            g_x_0_xyz_xxyy[k] = -g_x_0_xz_xxyy[k] * cd_y[k] + g_x_0_xz_xxyyy[k];

            g_x_0_xyz_xxyz[k] = -g_x_0_xz_xxyz[k] * cd_y[k] + g_x_0_xz_xxyyz[k];

            g_x_0_xyz_xxzz[k] = -g_x_0_xz_xxzz[k] * cd_y[k] + g_x_0_xz_xxyzz[k];

            g_x_0_xyz_xyyy[k] = -g_x_0_xz_xyyy[k] * cd_y[k] + g_x_0_xz_xyyyy[k];

            g_x_0_xyz_xyyz[k] = -g_x_0_xz_xyyz[k] * cd_y[k] + g_x_0_xz_xyyyz[k];

            g_x_0_xyz_xyzz[k] = -g_x_0_xz_xyzz[k] * cd_y[k] + g_x_0_xz_xyyzz[k];

            g_x_0_xyz_xzzz[k] = -g_x_0_xz_xzzz[k] * cd_y[k] + g_x_0_xz_xyzzz[k];

            g_x_0_xyz_yyyy[k] = -g_x_0_xz_yyyy[k] * cd_y[k] + g_x_0_xz_yyyyy[k];

            g_x_0_xyz_yyyz[k] = -g_x_0_xz_yyyz[k] * cd_y[k] + g_x_0_xz_yyyyz[k];

            g_x_0_xyz_yyzz[k] = -g_x_0_xz_yyzz[k] * cd_y[k] + g_x_0_xz_yyyzz[k];

            g_x_0_xyz_yzzz[k] = -g_x_0_xz_yzzz[k] * cd_y[k] + g_x_0_xz_yyzzz[k];

            g_x_0_xyz_zzzz[k] = -g_x_0_xz_zzzz[k] * cd_y[k] + g_x_0_xz_yzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 83);

        auto g_x_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_x_0_xz_xxxx, g_x_0_xz_xxxxz, g_x_0_xz_xxxy, g_x_0_xz_xxxyz, g_x_0_xz_xxxz, g_x_0_xz_xxxzz, g_x_0_xz_xxyy, g_x_0_xz_xxyyz, g_x_0_xz_xxyz, g_x_0_xz_xxyzz, g_x_0_xz_xxzz, g_x_0_xz_xxzzz, g_x_0_xz_xyyy, g_x_0_xz_xyyyz, g_x_0_xz_xyyz, g_x_0_xz_xyyzz, g_x_0_xz_xyzz, g_x_0_xz_xyzzz, g_x_0_xz_xzzz, g_x_0_xz_xzzzz, g_x_0_xz_yyyy, g_x_0_xz_yyyyz, g_x_0_xz_yyyz, g_x_0_xz_yyyzz, g_x_0_xz_yyzz, g_x_0_xz_yyzzz, g_x_0_xz_yzzz, g_x_0_xz_yzzzz, g_x_0_xz_zzzz, g_x_0_xz_zzzzz, g_x_0_xzz_xxxx, g_x_0_xzz_xxxy, g_x_0_xzz_xxxz, g_x_0_xzz_xxyy, g_x_0_xzz_xxyz, g_x_0_xzz_xxzz, g_x_0_xzz_xyyy, g_x_0_xzz_xyyz, g_x_0_xzz_xyzz, g_x_0_xzz_xzzz, g_x_0_xzz_yyyy, g_x_0_xzz_yyyz, g_x_0_xzz_yyzz, g_x_0_xzz_yzzz, g_x_0_xzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_xxxx[k] = -g_x_0_xz_xxxx[k] * cd_z[k] + g_x_0_xz_xxxxz[k];

            g_x_0_xzz_xxxy[k] = -g_x_0_xz_xxxy[k] * cd_z[k] + g_x_0_xz_xxxyz[k];

            g_x_0_xzz_xxxz[k] = -g_x_0_xz_xxxz[k] * cd_z[k] + g_x_0_xz_xxxzz[k];

            g_x_0_xzz_xxyy[k] = -g_x_0_xz_xxyy[k] * cd_z[k] + g_x_0_xz_xxyyz[k];

            g_x_0_xzz_xxyz[k] = -g_x_0_xz_xxyz[k] * cd_z[k] + g_x_0_xz_xxyzz[k];

            g_x_0_xzz_xxzz[k] = -g_x_0_xz_xxzz[k] * cd_z[k] + g_x_0_xz_xxzzz[k];

            g_x_0_xzz_xyyy[k] = -g_x_0_xz_xyyy[k] * cd_z[k] + g_x_0_xz_xyyyz[k];

            g_x_0_xzz_xyyz[k] = -g_x_0_xz_xyyz[k] * cd_z[k] + g_x_0_xz_xyyzz[k];

            g_x_0_xzz_xyzz[k] = -g_x_0_xz_xyzz[k] * cd_z[k] + g_x_0_xz_xyzzz[k];

            g_x_0_xzz_xzzz[k] = -g_x_0_xz_xzzz[k] * cd_z[k] + g_x_0_xz_xzzzz[k];

            g_x_0_xzz_yyyy[k] = -g_x_0_xz_yyyy[k] * cd_z[k] + g_x_0_xz_yyyyz[k];

            g_x_0_xzz_yyyz[k] = -g_x_0_xz_yyyz[k] * cd_z[k] + g_x_0_xz_yyyzz[k];

            g_x_0_xzz_yyzz[k] = -g_x_0_xz_yyzz[k] * cd_z[k] + g_x_0_xz_yyzzz[k];

            g_x_0_xzz_yzzz[k] = -g_x_0_xz_yzzz[k] * cd_z[k] + g_x_0_xz_yzzzz[k];

            g_x_0_xzz_zzzz[k] = -g_x_0_xz_zzzz[k] * cd_z[k] + g_x_0_xz_zzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 99);

        auto g_x_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 100);

        auto g_x_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 101);

        auto g_x_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 102);

        auto g_x_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 103);

        auto g_x_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_x_0_yy_xxxx, g_x_0_yy_xxxxy, g_x_0_yy_xxxy, g_x_0_yy_xxxyy, g_x_0_yy_xxxyz, g_x_0_yy_xxxz, g_x_0_yy_xxyy, g_x_0_yy_xxyyy, g_x_0_yy_xxyyz, g_x_0_yy_xxyz, g_x_0_yy_xxyzz, g_x_0_yy_xxzz, g_x_0_yy_xyyy, g_x_0_yy_xyyyy, g_x_0_yy_xyyyz, g_x_0_yy_xyyz, g_x_0_yy_xyyzz, g_x_0_yy_xyzz, g_x_0_yy_xyzzz, g_x_0_yy_xzzz, g_x_0_yy_yyyy, g_x_0_yy_yyyyy, g_x_0_yy_yyyyz, g_x_0_yy_yyyz, g_x_0_yy_yyyzz, g_x_0_yy_yyzz, g_x_0_yy_yyzzz, g_x_0_yy_yzzz, g_x_0_yy_yzzzz, g_x_0_yy_zzzz, g_x_0_yyy_xxxx, g_x_0_yyy_xxxy, g_x_0_yyy_xxxz, g_x_0_yyy_xxyy, g_x_0_yyy_xxyz, g_x_0_yyy_xxzz, g_x_0_yyy_xyyy, g_x_0_yyy_xyyz, g_x_0_yyy_xyzz, g_x_0_yyy_xzzz, g_x_0_yyy_yyyy, g_x_0_yyy_yyyz, g_x_0_yyy_yyzz, g_x_0_yyy_yzzz, g_x_0_yyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_xxxx[k] = -g_x_0_yy_xxxx[k] * cd_y[k] + g_x_0_yy_xxxxy[k];

            g_x_0_yyy_xxxy[k] = -g_x_0_yy_xxxy[k] * cd_y[k] + g_x_0_yy_xxxyy[k];

            g_x_0_yyy_xxxz[k] = -g_x_0_yy_xxxz[k] * cd_y[k] + g_x_0_yy_xxxyz[k];

            g_x_0_yyy_xxyy[k] = -g_x_0_yy_xxyy[k] * cd_y[k] + g_x_0_yy_xxyyy[k];

            g_x_0_yyy_xxyz[k] = -g_x_0_yy_xxyz[k] * cd_y[k] + g_x_0_yy_xxyyz[k];

            g_x_0_yyy_xxzz[k] = -g_x_0_yy_xxzz[k] * cd_y[k] + g_x_0_yy_xxyzz[k];

            g_x_0_yyy_xyyy[k] = -g_x_0_yy_xyyy[k] * cd_y[k] + g_x_0_yy_xyyyy[k];

            g_x_0_yyy_xyyz[k] = -g_x_0_yy_xyyz[k] * cd_y[k] + g_x_0_yy_xyyyz[k];

            g_x_0_yyy_xyzz[k] = -g_x_0_yy_xyzz[k] * cd_y[k] + g_x_0_yy_xyyzz[k];

            g_x_0_yyy_xzzz[k] = -g_x_0_yy_xzzz[k] * cd_y[k] + g_x_0_yy_xyzzz[k];

            g_x_0_yyy_yyyy[k] = -g_x_0_yy_yyyy[k] * cd_y[k] + g_x_0_yy_yyyyy[k];

            g_x_0_yyy_yyyz[k] = -g_x_0_yy_yyyz[k] * cd_y[k] + g_x_0_yy_yyyyz[k];

            g_x_0_yyy_yyzz[k] = -g_x_0_yy_yyzz[k] * cd_y[k] + g_x_0_yy_yyyzz[k];

            g_x_0_yyy_yzzz[k] = -g_x_0_yy_yzzz[k] * cd_y[k] + g_x_0_yy_yyzzz[k];

            g_x_0_yyy_zzzz[k] = -g_x_0_yy_zzzz[k] * cd_y[k] + g_x_0_yy_yzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 105);

        auto g_x_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 106);

        auto g_x_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 107);

        auto g_x_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 108);

        auto g_x_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 109);

        auto g_x_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 110);

        auto g_x_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 111);

        auto g_x_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 112);

        auto g_x_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 113);

        auto g_x_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 114);

        auto g_x_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 115);

        auto g_x_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 116);

        auto g_x_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 117);

        auto g_x_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 118);

        auto g_x_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 119);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_xxxx, g_x_0_yyz_xxxy, g_x_0_yyz_xxxz, g_x_0_yyz_xxyy, g_x_0_yyz_xxyz, g_x_0_yyz_xxzz, g_x_0_yyz_xyyy, g_x_0_yyz_xyyz, g_x_0_yyz_xyzz, g_x_0_yyz_xzzz, g_x_0_yyz_yyyy, g_x_0_yyz_yyyz, g_x_0_yyz_yyzz, g_x_0_yyz_yzzz, g_x_0_yyz_zzzz, g_x_0_yz_xxxx, g_x_0_yz_xxxxy, g_x_0_yz_xxxy, g_x_0_yz_xxxyy, g_x_0_yz_xxxyz, g_x_0_yz_xxxz, g_x_0_yz_xxyy, g_x_0_yz_xxyyy, g_x_0_yz_xxyyz, g_x_0_yz_xxyz, g_x_0_yz_xxyzz, g_x_0_yz_xxzz, g_x_0_yz_xyyy, g_x_0_yz_xyyyy, g_x_0_yz_xyyyz, g_x_0_yz_xyyz, g_x_0_yz_xyyzz, g_x_0_yz_xyzz, g_x_0_yz_xyzzz, g_x_0_yz_xzzz, g_x_0_yz_yyyy, g_x_0_yz_yyyyy, g_x_0_yz_yyyyz, g_x_0_yz_yyyz, g_x_0_yz_yyyzz, g_x_0_yz_yyzz, g_x_0_yz_yyzzz, g_x_0_yz_yzzz, g_x_0_yz_yzzzz, g_x_0_yz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_xxxx[k] = -g_x_0_yz_xxxx[k] * cd_y[k] + g_x_0_yz_xxxxy[k];

            g_x_0_yyz_xxxy[k] = -g_x_0_yz_xxxy[k] * cd_y[k] + g_x_0_yz_xxxyy[k];

            g_x_0_yyz_xxxz[k] = -g_x_0_yz_xxxz[k] * cd_y[k] + g_x_0_yz_xxxyz[k];

            g_x_0_yyz_xxyy[k] = -g_x_0_yz_xxyy[k] * cd_y[k] + g_x_0_yz_xxyyy[k];

            g_x_0_yyz_xxyz[k] = -g_x_0_yz_xxyz[k] * cd_y[k] + g_x_0_yz_xxyyz[k];

            g_x_0_yyz_xxzz[k] = -g_x_0_yz_xxzz[k] * cd_y[k] + g_x_0_yz_xxyzz[k];

            g_x_0_yyz_xyyy[k] = -g_x_0_yz_xyyy[k] * cd_y[k] + g_x_0_yz_xyyyy[k];

            g_x_0_yyz_xyyz[k] = -g_x_0_yz_xyyz[k] * cd_y[k] + g_x_0_yz_xyyyz[k];

            g_x_0_yyz_xyzz[k] = -g_x_0_yz_xyzz[k] * cd_y[k] + g_x_0_yz_xyyzz[k];

            g_x_0_yyz_xzzz[k] = -g_x_0_yz_xzzz[k] * cd_y[k] + g_x_0_yz_xyzzz[k];

            g_x_0_yyz_yyyy[k] = -g_x_0_yz_yyyy[k] * cd_y[k] + g_x_0_yz_yyyyy[k];

            g_x_0_yyz_yyyz[k] = -g_x_0_yz_yyyz[k] * cd_y[k] + g_x_0_yz_yyyyz[k];

            g_x_0_yyz_yyzz[k] = -g_x_0_yz_yyzz[k] * cd_y[k] + g_x_0_yz_yyyzz[k];

            g_x_0_yyz_yzzz[k] = -g_x_0_yz_yzzz[k] * cd_y[k] + g_x_0_yz_yyzzz[k];

            g_x_0_yyz_zzzz[k] = -g_x_0_yz_zzzz[k] * cd_y[k] + g_x_0_yz_yzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 120);

        auto g_x_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 121);

        auto g_x_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 122);

        auto g_x_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 123);

        auto g_x_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 124);

        auto g_x_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 125);

        auto g_x_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 126);

        auto g_x_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 127);

        auto g_x_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 128);

        auto g_x_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 129);

        auto g_x_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 130);

        auto g_x_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 131);

        auto g_x_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 132);

        auto g_x_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 133);

        auto g_x_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 134);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_xxxx, g_x_0_yzz_xxxy, g_x_0_yzz_xxxz, g_x_0_yzz_xxyy, g_x_0_yzz_xxyz, g_x_0_yzz_xxzz, g_x_0_yzz_xyyy, g_x_0_yzz_xyyz, g_x_0_yzz_xyzz, g_x_0_yzz_xzzz, g_x_0_yzz_yyyy, g_x_0_yzz_yyyz, g_x_0_yzz_yyzz, g_x_0_yzz_yzzz, g_x_0_yzz_zzzz, g_x_0_zz_xxxx, g_x_0_zz_xxxxy, g_x_0_zz_xxxy, g_x_0_zz_xxxyy, g_x_0_zz_xxxyz, g_x_0_zz_xxxz, g_x_0_zz_xxyy, g_x_0_zz_xxyyy, g_x_0_zz_xxyyz, g_x_0_zz_xxyz, g_x_0_zz_xxyzz, g_x_0_zz_xxzz, g_x_0_zz_xyyy, g_x_0_zz_xyyyy, g_x_0_zz_xyyyz, g_x_0_zz_xyyz, g_x_0_zz_xyyzz, g_x_0_zz_xyzz, g_x_0_zz_xyzzz, g_x_0_zz_xzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyyy, g_x_0_zz_yyyyz, g_x_0_zz_yyyz, g_x_0_zz_yyyzz, g_x_0_zz_yyzz, g_x_0_zz_yyzzz, g_x_0_zz_yzzz, g_x_0_zz_yzzzz, g_x_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_xxxx[k] = -g_x_0_zz_xxxx[k] * cd_y[k] + g_x_0_zz_xxxxy[k];

            g_x_0_yzz_xxxy[k] = -g_x_0_zz_xxxy[k] * cd_y[k] + g_x_0_zz_xxxyy[k];

            g_x_0_yzz_xxxz[k] = -g_x_0_zz_xxxz[k] * cd_y[k] + g_x_0_zz_xxxyz[k];

            g_x_0_yzz_xxyy[k] = -g_x_0_zz_xxyy[k] * cd_y[k] + g_x_0_zz_xxyyy[k];

            g_x_0_yzz_xxyz[k] = -g_x_0_zz_xxyz[k] * cd_y[k] + g_x_0_zz_xxyyz[k];

            g_x_0_yzz_xxzz[k] = -g_x_0_zz_xxzz[k] * cd_y[k] + g_x_0_zz_xxyzz[k];

            g_x_0_yzz_xyyy[k] = -g_x_0_zz_xyyy[k] * cd_y[k] + g_x_0_zz_xyyyy[k];

            g_x_0_yzz_xyyz[k] = -g_x_0_zz_xyyz[k] * cd_y[k] + g_x_0_zz_xyyyz[k];

            g_x_0_yzz_xyzz[k] = -g_x_0_zz_xyzz[k] * cd_y[k] + g_x_0_zz_xyyzz[k];

            g_x_0_yzz_xzzz[k] = -g_x_0_zz_xzzz[k] * cd_y[k] + g_x_0_zz_xyzzz[k];

            g_x_0_yzz_yyyy[k] = -g_x_0_zz_yyyy[k] * cd_y[k] + g_x_0_zz_yyyyy[k];

            g_x_0_yzz_yyyz[k] = -g_x_0_zz_yyyz[k] * cd_y[k] + g_x_0_zz_yyyyz[k];

            g_x_0_yzz_yyzz[k] = -g_x_0_zz_yyzz[k] * cd_y[k] + g_x_0_zz_yyyzz[k];

            g_x_0_yzz_yzzz[k] = -g_x_0_zz_yzzz[k] * cd_y[k] + g_x_0_zz_yyzzz[k];

            g_x_0_yzz_zzzz[k] = -g_x_0_zz_zzzz[k] * cd_y[k] + g_x_0_zz_yzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps  + 135);

        auto g_x_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 136);

        auto g_x_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 137);

        auto g_x_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 138);

        auto g_x_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 139);

        auto g_x_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 140);

        auto g_x_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 141);

        auto g_x_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 142);

        auto g_x_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 143);

        auto g_x_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 144);

        auto g_x_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps  + 145);

        auto g_x_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 146);

        auto g_x_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 147);

        auto g_x_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 148);

        auto g_x_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps  + 149);

        #pragma omp simd aligned(cd_z, g_x_0_zz_xxxx, g_x_0_zz_xxxxz, g_x_0_zz_xxxy, g_x_0_zz_xxxyz, g_x_0_zz_xxxz, g_x_0_zz_xxxzz, g_x_0_zz_xxyy, g_x_0_zz_xxyyz, g_x_0_zz_xxyz, g_x_0_zz_xxyzz, g_x_0_zz_xxzz, g_x_0_zz_xxzzz, g_x_0_zz_xyyy, g_x_0_zz_xyyyz, g_x_0_zz_xyyz, g_x_0_zz_xyyzz, g_x_0_zz_xyzz, g_x_0_zz_xyzzz, g_x_0_zz_xzzz, g_x_0_zz_xzzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyyz, g_x_0_zz_yyyz, g_x_0_zz_yyyzz, g_x_0_zz_yyzz, g_x_0_zz_yyzzz, g_x_0_zz_yzzz, g_x_0_zz_yzzzz, g_x_0_zz_zzzz, g_x_0_zz_zzzzz, g_x_0_zzz_xxxx, g_x_0_zzz_xxxy, g_x_0_zzz_xxxz, g_x_0_zzz_xxyy, g_x_0_zzz_xxyz, g_x_0_zzz_xxzz, g_x_0_zzz_xyyy, g_x_0_zzz_xyyz, g_x_0_zzz_xyzz, g_x_0_zzz_xzzz, g_x_0_zzz_yyyy, g_x_0_zzz_yyyz, g_x_0_zzz_yyzz, g_x_0_zzz_yzzz, g_x_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_xxxx[k] = -g_x_0_zz_xxxx[k] * cd_z[k] + g_x_0_zz_xxxxz[k];

            g_x_0_zzz_xxxy[k] = -g_x_0_zz_xxxy[k] * cd_z[k] + g_x_0_zz_xxxyz[k];

            g_x_0_zzz_xxxz[k] = -g_x_0_zz_xxxz[k] * cd_z[k] + g_x_0_zz_xxxzz[k];

            g_x_0_zzz_xxyy[k] = -g_x_0_zz_xxyy[k] * cd_z[k] + g_x_0_zz_xxyyz[k];

            g_x_0_zzz_xxyz[k] = -g_x_0_zz_xxyz[k] * cd_z[k] + g_x_0_zz_xxyzz[k];

            g_x_0_zzz_xxzz[k] = -g_x_0_zz_xxzz[k] * cd_z[k] + g_x_0_zz_xxzzz[k];

            g_x_0_zzz_xyyy[k] = -g_x_0_zz_xyyy[k] * cd_z[k] + g_x_0_zz_xyyyz[k];

            g_x_0_zzz_xyyz[k] = -g_x_0_zz_xyyz[k] * cd_z[k] + g_x_0_zz_xyyzz[k];

            g_x_0_zzz_xyzz[k] = -g_x_0_zz_xyzz[k] * cd_z[k] + g_x_0_zz_xyzzz[k];

            g_x_0_zzz_xzzz[k] = -g_x_0_zz_xzzz[k] * cd_z[k] + g_x_0_zz_xzzzz[k];

            g_x_0_zzz_yyyy[k] = -g_x_0_zz_yyyy[k] * cd_z[k] + g_x_0_zz_yyyyz[k];

            g_x_0_zzz_yyyz[k] = -g_x_0_zz_yyyz[k] * cd_z[k] + g_x_0_zz_yyyzz[k];

            g_x_0_zzz_yyzz[k] = -g_x_0_zz_yyzz[k] * cd_z[k] + g_x_0_zz_yyzzz[k];

            g_x_0_zzz_yzzz[k] = -g_x_0_zz_yzzz[k] * cd_z[k] + g_x_0_zz_yzzzz[k];

            g_x_0_zzz_zzzz[k] = -g_x_0_zz_zzzz[k] * cd_z[k] + g_x_0_zz_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 0);

        auto g_y_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 1);

        auto g_y_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 2);

        auto g_y_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 3);

        auto g_y_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 4);

        auto g_y_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 5);

        auto g_y_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 6);

        auto g_y_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 7);

        auto g_y_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 8);

        auto g_y_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 9);

        auto g_y_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 10);

        auto g_y_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 11);

        auto g_y_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 12);

        auto g_y_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 13);

        auto g_y_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_y_0_xx_xxxx, g_y_0_xx_xxxxx, g_y_0_xx_xxxxy, g_y_0_xx_xxxxz, g_y_0_xx_xxxy, g_y_0_xx_xxxyy, g_y_0_xx_xxxyz, g_y_0_xx_xxxz, g_y_0_xx_xxxzz, g_y_0_xx_xxyy, g_y_0_xx_xxyyy, g_y_0_xx_xxyyz, g_y_0_xx_xxyz, g_y_0_xx_xxyzz, g_y_0_xx_xxzz, g_y_0_xx_xxzzz, g_y_0_xx_xyyy, g_y_0_xx_xyyyy, g_y_0_xx_xyyyz, g_y_0_xx_xyyz, g_y_0_xx_xyyzz, g_y_0_xx_xyzz, g_y_0_xx_xyzzz, g_y_0_xx_xzzz, g_y_0_xx_xzzzz, g_y_0_xx_yyyy, g_y_0_xx_yyyz, g_y_0_xx_yyzz, g_y_0_xx_yzzz, g_y_0_xx_zzzz, g_y_0_xxx_xxxx, g_y_0_xxx_xxxy, g_y_0_xxx_xxxz, g_y_0_xxx_xxyy, g_y_0_xxx_xxyz, g_y_0_xxx_xxzz, g_y_0_xxx_xyyy, g_y_0_xxx_xyyz, g_y_0_xxx_xyzz, g_y_0_xxx_xzzz, g_y_0_xxx_yyyy, g_y_0_xxx_yyyz, g_y_0_xxx_yyzz, g_y_0_xxx_yzzz, g_y_0_xxx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_xxxx[k] = -g_y_0_xx_xxxx[k] * cd_x[k] + g_y_0_xx_xxxxx[k];

            g_y_0_xxx_xxxy[k] = -g_y_0_xx_xxxy[k] * cd_x[k] + g_y_0_xx_xxxxy[k];

            g_y_0_xxx_xxxz[k] = -g_y_0_xx_xxxz[k] * cd_x[k] + g_y_0_xx_xxxxz[k];

            g_y_0_xxx_xxyy[k] = -g_y_0_xx_xxyy[k] * cd_x[k] + g_y_0_xx_xxxyy[k];

            g_y_0_xxx_xxyz[k] = -g_y_0_xx_xxyz[k] * cd_x[k] + g_y_0_xx_xxxyz[k];

            g_y_0_xxx_xxzz[k] = -g_y_0_xx_xxzz[k] * cd_x[k] + g_y_0_xx_xxxzz[k];

            g_y_0_xxx_xyyy[k] = -g_y_0_xx_xyyy[k] * cd_x[k] + g_y_0_xx_xxyyy[k];

            g_y_0_xxx_xyyz[k] = -g_y_0_xx_xyyz[k] * cd_x[k] + g_y_0_xx_xxyyz[k];

            g_y_0_xxx_xyzz[k] = -g_y_0_xx_xyzz[k] * cd_x[k] + g_y_0_xx_xxyzz[k];

            g_y_0_xxx_xzzz[k] = -g_y_0_xx_xzzz[k] * cd_x[k] + g_y_0_xx_xxzzz[k];

            g_y_0_xxx_yyyy[k] = -g_y_0_xx_yyyy[k] * cd_x[k] + g_y_0_xx_xyyyy[k];

            g_y_0_xxx_yyyz[k] = -g_y_0_xx_yyyz[k] * cd_x[k] + g_y_0_xx_xyyyz[k];

            g_y_0_xxx_yyzz[k] = -g_y_0_xx_yyzz[k] * cd_x[k] + g_y_0_xx_xyyzz[k];

            g_y_0_xxx_yzzz[k] = -g_y_0_xx_yzzz[k] * cd_x[k] + g_y_0_xx_xyzzz[k];

            g_y_0_xxx_zzzz[k] = -g_y_0_xx_zzzz[k] * cd_x[k] + g_y_0_xx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 15);

        auto g_y_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 16);

        auto g_y_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 17);

        auto g_y_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 18);

        auto g_y_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 19);

        auto g_y_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 20);

        auto g_y_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 21);

        auto g_y_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 22);

        auto g_y_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 23);

        auto g_y_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 24);

        auto g_y_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 25);

        auto g_y_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 26);

        auto g_y_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 27);

        auto g_y_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 28);

        auto g_y_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_xxxx, g_y_0_xxy_xxxy, g_y_0_xxy_xxxz, g_y_0_xxy_xxyy, g_y_0_xxy_xxyz, g_y_0_xxy_xxzz, g_y_0_xxy_xyyy, g_y_0_xxy_xyyz, g_y_0_xxy_xyzz, g_y_0_xxy_xzzz, g_y_0_xxy_yyyy, g_y_0_xxy_yyyz, g_y_0_xxy_yyzz, g_y_0_xxy_yzzz, g_y_0_xxy_zzzz, g_y_0_xy_xxxx, g_y_0_xy_xxxxx, g_y_0_xy_xxxxy, g_y_0_xy_xxxxz, g_y_0_xy_xxxy, g_y_0_xy_xxxyy, g_y_0_xy_xxxyz, g_y_0_xy_xxxz, g_y_0_xy_xxxzz, g_y_0_xy_xxyy, g_y_0_xy_xxyyy, g_y_0_xy_xxyyz, g_y_0_xy_xxyz, g_y_0_xy_xxyzz, g_y_0_xy_xxzz, g_y_0_xy_xxzzz, g_y_0_xy_xyyy, g_y_0_xy_xyyyy, g_y_0_xy_xyyyz, g_y_0_xy_xyyz, g_y_0_xy_xyyzz, g_y_0_xy_xyzz, g_y_0_xy_xyzzz, g_y_0_xy_xzzz, g_y_0_xy_xzzzz, g_y_0_xy_yyyy, g_y_0_xy_yyyz, g_y_0_xy_yyzz, g_y_0_xy_yzzz, g_y_0_xy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_xxxx[k] = -g_y_0_xy_xxxx[k] * cd_x[k] + g_y_0_xy_xxxxx[k];

            g_y_0_xxy_xxxy[k] = -g_y_0_xy_xxxy[k] * cd_x[k] + g_y_0_xy_xxxxy[k];

            g_y_0_xxy_xxxz[k] = -g_y_0_xy_xxxz[k] * cd_x[k] + g_y_0_xy_xxxxz[k];

            g_y_0_xxy_xxyy[k] = -g_y_0_xy_xxyy[k] * cd_x[k] + g_y_0_xy_xxxyy[k];

            g_y_0_xxy_xxyz[k] = -g_y_0_xy_xxyz[k] * cd_x[k] + g_y_0_xy_xxxyz[k];

            g_y_0_xxy_xxzz[k] = -g_y_0_xy_xxzz[k] * cd_x[k] + g_y_0_xy_xxxzz[k];

            g_y_0_xxy_xyyy[k] = -g_y_0_xy_xyyy[k] * cd_x[k] + g_y_0_xy_xxyyy[k];

            g_y_0_xxy_xyyz[k] = -g_y_0_xy_xyyz[k] * cd_x[k] + g_y_0_xy_xxyyz[k];

            g_y_0_xxy_xyzz[k] = -g_y_0_xy_xyzz[k] * cd_x[k] + g_y_0_xy_xxyzz[k];

            g_y_0_xxy_xzzz[k] = -g_y_0_xy_xzzz[k] * cd_x[k] + g_y_0_xy_xxzzz[k];

            g_y_0_xxy_yyyy[k] = -g_y_0_xy_yyyy[k] * cd_x[k] + g_y_0_xy_xyyyy[k];

            g_y_0_xxy_yyyz[k] = -g_y_0_xy_yyyz[k] * cd_x[k] + g_y_0_xy_xyyyz[k];

            g_y_0_xxy_yyzz[k] = -g_y_0_xy_yyzz[k] * cd_x[k] + g_y_0_xy_xyyzz[k];

            g_y_0_xxy_yzzz[k] = -g_y_0_xy_yzzz[k] * cd_x[k] + g_y_0_xy_xyzzz[k];

            g_y_0_xxy_zzzz[k] = -g_y_0_xy_zzzz[k] * cd_x[k] + g_y_0_xy_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 30);

        auto g_y_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 31);

        auto g_y_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 32);

        auto g_y_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 33);

        auto g_y_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 34);

        auto g_y_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 35);

        auto g_y_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 36);

        auto g_y_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 37);

        auto g_y_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 38);

        auto g_y_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 39);

        auto g_y_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 40);

        auto g_y_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 41);

        auto g_y_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 42);

        auto g_y_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 43);

        auto g_y_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_xxxx, g_y_0_xxz_xxxy, g_y_0_xxz_xxxz, g_y_0_xxz_xxyy, g_y_0_xxz_xxyz, g_y_0_xxz_xxzz, g_y_0_xxz_xyyy, g_y_0_xxz_xyyz, g_y_0_xxz_xyzz, g_y_0_xxz_xzzz, g_y_0_xxz_yyyy, g_y_0_xxz_yyyz, g_y_0_xxz_yyzz, g_y_0_xxz_yzzz, g_y_0_xxz_zzzz, g_y_0_xz_xxxx, g_y_0_xz_xxxxx, g_y_0_xz_xxxxy, g_y_0_xz_xxxxz, g_y_0_xz_xxxy, g_y_0_xz_xxxyy, g_y_0_xz_xxxyz, g_y_0_xz_xxxz, g_y_0_xz_xxxzz, g_y_0_xz_xxyy, g_y_0_xz_xxyyy, g_y_0_xz_xxyyz, g_y_0_xz_xxyz, g_y_0_xz_xxyzz, g_y_0_xz_xxzz, g_y_0_xz_xxzzz, g_y_0_xz_xyyy, g_y_0_xz_xyyyy, g_y_0_xz_xyyyz, g_y_0_xz_xyyz, g_y_0_xz_xyyzz, g_y_0_xz_xyzz, g_y_0_xz_xyzzz, g_y_0_xz_xzzz, g_y_0_xz_xzzzz, g_y_0_xz_yyyy, g_y_0_xz_yyyz, g_y_0_xz_yyzz, g_y_0_xz_yzzz, g_y_0_xz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_xxxx[k] = -g_y_0_xz_xxxx[k] * cd_x[k] + g_y_0_xz_xxxxx[k];

            g_y_0_xxz_xxxy[k] = -g_y_0_xz_xxxy[k] * cd_x[k] + g_y_0_xz_xxxxy[k];

            g_y_0_xxz_xxxz[k] = -g_y_0_xz_xxxz[k] * cd_x[k] + g_y_0_xz_xxxxz[k];

            g_y_0_xxz_xxyy[k] = -g_y_0_xz_xxyy[k] * cd_x[k] + g_y_0_xz_xxxyy[k];

            g_y_0_xxz_xxyz[k] = -g_y_0_xz_xxyz[k] * cd_x[k] + g_y_0_xz_xxxyz[k];

            g_y_0_xxz_xxzz[k] = -g_y_0_xz_xxzz[k] * cd_x[k] + g_y_0_xz_xxxzz[k];

            g_y_0_xxz_xyyy[k] = -g_y_0_xz_xyyy[k] * cd_x[k] + g_y_0_xz_xxyyy[k];

            g_y_0_xxz_xyyz[k] = -g_y_0_xz_xyyz[k] * cd_x[k] + g_y_0_xz_xxyyz[k];

            g_y_0_xxz_xyzz[k] = -g_y_0_xz_xyzz[k] * cd_x[k] + g_y_0_xz_xxyzz[k];

            g_y_0_xxz_xzzz[k] = -g_y_0_xz_xzzz[k] * cd_x[k] + g_y_0_xz_xxzzz[k];

            g_y_0_xxz_yyyy[k] = -g_y_0_xz_yyyy[k] * cd_x[k] + g_y_0_xz_xyyyy[k];

            g_y_0_xxz_yyyz[k] = -g_y_0_xz_yyyz[k] * cd_x[k] + g_y_0_xz_xyyyz[k];

            g_y_0_xxz_yyzz[k] = -g_y_0_xz_yyzz[k] * cd_x[k] + g_y_0_xz_xyyzz[k];

            g_y_0_xxz_yzzz[k] = -g_y_0_xz_yzzz[k] * cd_x[k] + g_y_0_xz_xyzzz[k];

            g_y_0_xxz_zzzz[k] = -g_y_0_xz_zzzz[k] * cd_x[k] + g_y_0_xz_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 45);

        auto g_y_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 46);

        auto g_y_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 47);

        auto g_y_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 48);

        auto g_y_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 49);

        auto g_y_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 50);

        auto g_y_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 51);

        auto g_y_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 52);

        auto g_y_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 53);

        auto g_y_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 54);

        auto g_y_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 55);

        auto g_y_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 56);

        auto g_y_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 57);

        auto g_y_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 58);

        auto g_y_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_xxxx, g_y_0_xyy_xxxy, g_y_0_xyy_xxxz, g_y_0_xyy_xxyy, g_y_0_xyy_xxyz, g_y_0_xyy_xxzz, g_y_0_xyy_xyyy, g_y_0_xyy_xyyz, g_y_0_xyy_xyzz, g_y_0_xyy_xzzz, g_y_0_xyy_yyyy, g_y_0_xyy_yyyz, g_y_0_xyy_yyzz, g_y_0_xyy_yzzz, g_y_0_xyy_zzzz, g_y_0_yy_xxxx, g_y_0_yy_xxxxx, g_y_0_yy_xxxxy, g_y_0_yy_xxxxz, g_y_0_yy_xxxy, g_y_0_yy_xxxyy, g_y_0_yy_xxxyz, g_y_0_yy_xxxz, g_y_0_yy_xxxzz, g_y_0_yy_xxyy, g_y_0_yy_xxyyy, g_y_0_yy_xxyyz, g_y_0_yy_xxyz, g_y_0_yy_xxyzz, g_y_0_yy_xxzz, g_y_0_yy_xxzzz, g_y_0_yy_xyyy, g_y_0_yy_xyyyy, g_y_0_yy_xyyyz, g_y_0_yy_xyyz, g_y_0_yy_xyyzz, g_y_0_yy_xyzz, g_y_0_yy_xyzzz, g_y_0_yy_xzzz, g_y_0_yy_xzzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_xxxx[k] = -g_y_0_yy_xxxx[k] * cd_x[k] + g_y_0_yy_xxxxx[k];

            g_y_0_xyy_xxxy[k] = -g_y_0_yy_xxxy[k] * cd_x[k] + g_y_0_yy_xxxxy[k];

            g_y_0_xyy_xxxz[k] = -g_y_0_yy_xxxz[k] * cd_x[k] + g_y_0_yy_xxxxz[k];

            g_y_0_xyy_xxyy[k] = -g_y_0_yy_xxyy[k] * cd_x[k] + g_y_0_yy_xxxyy[k];

            g_y_0_xyy_xxyz[k] = -g_y_0_yy_xxyz[k] * cd_x[k] + g_y_0_yy_xxxyz[k];

            g_y_0_xyy_xxzz[k] = -g_y_0_yy_xxzz[k] * cd_x[k] + g_y_0_yy_xxxzz[k];

            g_y_0_xyy_xyyy[k] = -g_y_0_yy_xyyy[k] * cd_x[k] + g_y_0_yy_xxyyy[k];

            g_y_0_xyy_xyyz[k] = -g_y_0_yy_xyyz[k] * cd_x[k] + g_y_0_yy_xxyyz[k];

            g_y_0_xyy_xyzz[k] = -g_y_0_yy_xyzz[k] * cd_x[k] + g_y_0_yy_xxyzz[k];

            g_y_0_xyy_xzzz[k] = -g_y_0_yy_xzzz[k] * cd_x[k] + g_y_0_yy_xxzzz[k];

            g_y_0_xyy_yyyy[k] = -g_y_0_yy_yyyy[k] * cd_x[k] + g_y_0_yy_xyyyy[k];

            g_y_0_xyy_yyyz[k] = -g_y_0_yy_yyyz[k] * cd_x[k] + g_y_0_yy_xyyyz[k];

            g_y_0_xyy_yyzz[k] = -g_y_0_yy_yyzz[k] * cd_x[k] + g_y_0_yy_xyyzz[k];

            g_y_0_xyy_yzzz[k] = -g_y_0_yy_yzzz[k] * cd_x[k] + g_y_0_yy_xyzzz[k];

            g_y_0_xyy_zzzz[k] = -g_y_0_yy_zzzz[k] * cd_x[k] + g_y_0_yy_xzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 60);

        auto g_y_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 61);

        auto g_y_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 62);

        auto g_y_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 63);

        auto g_y_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 64);

        auto g_y_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 65);

        auto g_y_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 66);

        auto g_y_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 67);

        auto g_y_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 68);

        auto g_y_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 69);

        auto g_y_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 70);

        auto g_y_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 71);

        auto g_y_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 72);

        auto g_y_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 73);

        auto g_y_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 74);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_xxxx, g_y_0_xyz_xxxy, g_y_0_xyz_xxxz, g_y_0_xyz_xxyy, g_y_0_xyz_xxyz, g_y_0_xyz_xxzz, g_y_0_xyz_xyyy, g_y_0_xyz_xyyz, g_y_0_xyz_xyzz, g_y_0_xyz_xzzz, g_y_0_xyz_yyyy, g_y_0_xyz_yyyz, g_y_0_xyz_yyzz, g_y_0_xyz_yzzz, g_y_0_xyz_zzzz, g_y_0_yz_xxxx, g_y_0_yz_xxxxx, g_y_0_yz_xxxxy, g_y_0_yz_xxxxz, g_y_0_yz_xxxy, g_y_0_yz_xxxyy, g_y_0_yz_xxxyz, g_y_0_yz_xxxz, g_y_0_yz_xxxzz, g_y_0_yz_xxyy, g_y_0_yz_xxyyy, g_y_0_yz_xxyyz, g_y_0_yz_xxyz, g_y_0_yz_xxyzz, g_y_0_yz_xxzz, g_y_0_yz_xxzzz, g_y_0_yz_xyyy, g_y_0_yz_xyyyy, g_y_0_yz_xyyyz, g_y_0_yz_xyyz, g_y_0_yz_xyyzz, g_y_0_yz_xyzz, g_y_0_yz_xyzzz, g_y_0_yz_xzzz, g_y_0_yz_xzzzz, g_y_0_yz_yyyy, g_y_0_yz_yyyz, g_y_0_yz_yyzz, g_y_0_yz_yzzz, g_y_0_yz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_xxxx[k] = -g_y_0_yz_xxxx[k] * cd_x[k] + g_y_0_yz_xxxxx[k];

            g_y_0_xyz_xxxy[k] = -g_y_0_yz_xxxy[k] * cd_x[k] + g_y_0_yz_xxxxy[k];

            g_y_0_xyz_xxxz[k] = -g_y_0_yz_xxxz[k] * cd_x[k] + g_y_0_yz_xxxxz[k];

            g_y_0_xyz_xxyy[k] = -g_y_0_yz_xxyy[k] * cd_x[k] + g_y_0_yz_xxxyy[k];

            g_y_0_xyz_xxyz[k] = -g_y_0_yz_xxyz[k] * cd_x[k] + g_y_0_yz_xxxyz[k];

            g_y_0_xyz_xxzz[k] = -g_y_0_yz_xxzz[k] * cd_x[k] + g_y_0_yz_xxxzz[k];

            g_y_0_xyz_xyyy[k] = -g_y_0_yz_xyyy[k] * cd_x[k] + g_y_0_yz_xxyyy[k];

            g_y_0_xyz_xyyz[k] = -g_y_0_yz_xyyz[k] * cd_x[k] + g_y_0_yz_xxyyz[k];

            g_y_0_xyz_xyzz[k] = -g_y_0_yz_xyzz[k] * cd_x[k] + g_y_0_yz_xxyzz[k];

            g_y_0_xyz_xzzz[k] = -g_y_0_yz_xzzz[k] * cd_x[k] + g_y_0_yz_xxzzz[k];

            g_y_0_xyz_yyyy[k] = -g_y_0_yz_yyyy[k] * cd_x[k] + g_y_0_yz_xyyyy[k];

            g_y_0_xyz_yyyz[k] = -g_y_0_yz_yyyz[k] * cd_x[k] + g_y_0_yz_xyyyz[k];

            g_y_0_xyz_yyzz[k] = -g_y_0_yz_yyzz[k] * cd_x[k] + g_y_0_yz_xyyzz[k];

            g_y_0_xyz_yzzz[k] = -g_y_0_yz_yzzz[k] * cd_x[k] + g_y_0_yz_xyzzz[k];

            g_y_0_xyz_zzzz[k] = -g_y_0_yz_zzzz[k] * cd_x[k] + g_y_0_yz_xzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 75);

        auto g_y_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 76);

        auto g_y_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 77);

        auto g_y_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 78);

        auto g_y_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 79);

        auto g_y_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 80);

        auto g_y_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 81);

        auto g_y_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 82);

        auto g_y_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 83);

        auto g_y_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 84);

        auto g_y_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 85);

        auto g_y_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 86);

        auto g_y_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 87);

        auto g_y_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 88);

        auto g_y_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 89);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_xxxx, g_y_0_xzz_xxxy, g_y_0_xzz_xxxz, g_y_0_xzz_xxyy, g_y_0_xzz_xxyz, g_y_0_xzz_xxzz, g_y_0_xzz_xyyy, g_y_0_xzz_xyyz, g_y_0_xzz_xyzz, g_y_0_xzz_xzzz, g_y_0_xzz_yyyy, g_y_0_xzz_yyyz, g_y_0_xzz_yyzz, g_y_0_xzz_yzzz, g_y_0_xzz_zzzz, g_y_0_zz_xxxx, g_y_0_zz_xxxxx, g_y_0_zz_xxxxy, g_y_0_zz_xxxxz, g_y_0_zz_xxxy, g_y_0_zz_xxxyy, g_y_0_zz_xxxyz, g_y_0_zz_xxxz, g_y_0_zz_xxxzz, g_y_0_zz_xxyy, g_y_0_zz_xxyyy, g_y_0_zz_xxyyz, g_y_0_zz_xxyz, g_y_0_zz_xxyzz, g_y_0_zz_xxzz, g_y_0_zz_xxzzz, g_y_0_zz_xyyy, g_y_0_zz_xyyyy, g_y_0_zz_xyyyz, g_y_0_zz_xyyz, g_y_0_zz_xyyzz, g_y_0_zz_xyzz, g_y_0_zz_xyzzz, g_y_0_zz_xzzz, g_y_0_zz_xzzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyz, g_y_0_zz_yyzz, g_y_0_zz_yzzz, g_y_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_xxxx[k] = -g_y_0_zz_xxxx[k] * cd_x[k] + g_y_0_zz_xxxxx[k];

            g_y_0_xzz_xxxy[k] = -g_y_0_zz_xxxy[k] * cd_x[k] + g_y_0_zz_xxxxy[k];

            g_y_0_xzz_xxxz[k] = -g_y_0_zz_xxxz[k] * cd_x[k] + g_y_0_zz_xxxxz[k];

            g_y_0_xzz_xxyy[k] = -g_y_0_zz_xxyy[k] * cd_x[k] + g_y_0_zz_xxxyy[k];

            g_y_0_xzz_xxyz[k] = -g_y_0_zz_xxyz[k] * cd_x[k] + g_y_0_zz_xxxyz[k];

            g_y_0_xzz_xxzz[k] = -g_y_0_zz_xxzz[k] * cd_x[k] + g_y_0_zz_xxxzz[k];

            g_y_0_xzz_xyyy[k] = -g_y_0_zz_xyyy[k] * cd_x[k] + g_y_0_zz_xxyyy[k];

            g_y_0_xzz_xyyz[k] = -g_y_0_zz_xyyz[k] * cd_x[k] + g_y_0_zz_xxyyz[k];

            g_y_0_xzz_xyzz[k] = -g_y_0_zz_xyzz[k] * cd_x[k] + g_y_0_zz_xxyzz[k];

            g_y_0_xzz_xzzz[k] = -g_y_0_zz_xzzz[k] * cd_x[k] + g_y_0_zz_xxzzz[k];

            g_y_0_xzz_yyyy[k] = -g_y_0_zz_yyyy[k] * cd_x[k] + g_y_0_zz_xyyyy[k];

            g_y_0_xzz_yyyz[k] = -g_y_0_zz_yyyz[k] * cd_x[k] + g_y_0_zz_xyyyz[k];

            g_y_0_xzz_yyzz[k] = -g_y_0_zz_yyzz[k] * cd_x[k] + g_y_0_zz_xyyzz[k];

            g_y_0_xzz_yzzz[k] = -g_y_0_zz_yzzz[k] * cd_x[k] + g_y_0_zz_xyzzz[k];

            g_y_0_xzz_zzzz[k] = -g_y_0_zz_zzzz[k] * cd_x[k] + g_y_0_zz_xzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 90);

        auto g_y_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 91);

        auto g_y_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 92);

        auto g_y_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 93);

        auto g_y_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 94);

        auto g_y_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 95);

        auto g_y_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 96);

        auto g_y_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 97);

        auto g_y_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 98);

        auto g_y_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 99);

        auto g_y_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 100);

        auto g_y_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 101);

        auto g_y_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 102);

        auto g_y_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 103);

        auto g_y_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_y_0_yy_xxxx, g_y_0_yy_xxxxy, g_y_0_yy_xxxy, g_y_0_yy_xxxyy, g_y_0_yy_xxxyz, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyyy, g_y_0_yy_xxyyz, g_y_0_yy_xxyz, g_y_0_yy_xxyzz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyyy, g_y_0_yy_xyyyz, g_y_0_yy_xyyz, g_y_0_yy_xyyzz, g_y_0_yy_xyzz, g_y_0_yy_xyzzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyyy, g_y_0_yy_yyyyz, g_y_0_yy_yyyz, g_y_0_yy_yyyzz, g_y_0_yy_yyzz, g_y_0_yy_yyzzz, g_y_0_yy_yzzz, g_y_0_yy_yzzzz, g_y_0_yy_zzzz, g_y_0_yyy_xxxx, g_y_0_yyy_xxxy, g_y_0_yyy_xxxz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyz, g_y_0_yyy_xxzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyz, g_y_0_yyy_xyzz, g_y_0_yyy_xzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyz, g_y_0_yyy_yyzz, g_y_0_yyy_yzzz, g_y_0_yyy_zzzz, g_yy_xxxx, g_yy_xxxy, g_yy_xxxz, g_yy_xxyy, g_yy_xxyz, g_yy_xxzz, g_yy_xyyy, g_yy_xyyz, g_yy_xyzz, g_yy_xzzz, g_yy_yyyy, g_yy_yyyz, g_yy_yyzz, g_yy_yzzz, g_yy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_xxxx[k] = -g_yy_xxxx[k] - g_y_0_yy_xxxx[k] * cd_y[k] + g_y_0_yy_xxxxy[k];

            g_y_0_yyy_xxxy[k] = -g_yy_xxxy[k] - g_y_0_yy_xxxy[k] * cd_y[k] + g_y_0_yy_xxxyy[k];

            g_y_0_yyy_xxxz[k] = -g_yy_xxxz[k] - g_y_0_yy_xxxz[k] * cd_y[k] + g_y_0_yy_xxxyz[k];

            g_y_0_yyy_xxyy[k] = -g_yy_xxyy[k] - g_y_0_yy_xxyy[k] * cd_y[k] + g_y_0_yy_xxyyy[k];

            g_y_0_yyy_xxyz[k] = -g_yy_xxyz[k] - g_y_0_yy_xxyz[k] * cd_y[k] + g_y_0_yy_xxyyz[k];

            g_y_0_yyy_xxzz[k] = -g_yy_xxzz[k] - g_y_0_yy_xxzz[k] * cd_y[k] + g_y_0_yy_xxyzz[k];

            g_y_0_yyy_xyyy[k] = -g_yy_xyyy[k] - g_y_0_yy_xyyy[k] * cd_y[k] + g_y_0_yy_xyyyy[k];

            g_y_0_yyy_xyyz[k] = -g_yy_xyyz[k] - g_y_0_yy_xyyz[k] * cd_y[k] + g_y_0_yy_xyyyz[k];

            g_y_0_yyy_xyzz[k] = -g_yy_xyzz[k] - g_y_0_yy_xyzz[k] * cd_y[k] + g_y_0_yy_xyyzz[k];

            g_y_0_yyy_xzzz[k] = -g_yy_xzzz[k] - g_y_0_yy_xzzz[k] * cd_y[k] + g_y_0_yy_xyzzz[k];

            g_y_0_yyy_yyyy[k] = -g_yy_yyyy[k] - g_y_0_yy_yyyy[k] * cd_y[k] + g_y_0_yy_yyyyy[k];

            g_y_0_yyy_yyyz[k] = -g_yy_yyyz[k] - g_y_0_yy_yyyz[k] * cd_y[k] + g_y_0_yy_yyyyz[k];

            g_y_0_yyy_yyzz[k] = -g_yy_yyzz[k] - g_y_0_yy_yyzz[k] * cd_y[k] + g_y_0_yy_yyyzz[k];

            g_y_0_yyy_yzzz[k] = -g_yy_yzzz[k] - g_y_0_yy_yzzz[k] * cd_y[k] + g_y_0_yy_yyzzz[k];

            g_y_0_yyy_zzzz[k] = -g_yy_zzzz[k] - g_y_0_yy_zzzz[k] * cd_y[k] + g_y_0_yy_yzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 105);

        auto g_y_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 106);

        auto g_y_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 107);

        auto g_y_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 108);

        auto g_y_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 109);

        auto g_y_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 110);

        auto g_y_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 111);

        auto g_y_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 112);

        auto g_y_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 113);

        auto g_y_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 114);

        auto g_y_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 115);

        auto g_y_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 116);

        auto g_y_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 117);

        auto g_y_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 118);

        auto g_y_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 119);

        #pragma omp simd aligned(cd_z, g_y_0_yy_xxxx, g_y_0_yy_xxxxz, g_y_0_yy_xxxy, g_y_0_yy_xxxyz, g_y_0_yy_xxxz, g_y_0_yy_xxxzz, g_y_0_yy_xxyy, g_y_0_yy_xxyyz, g_y_0_yy_xxyz, g_y_0_yy_xxyzz, g_y_0_yy_xxzz, g_y_0_yy_xxzzz, g_y_0_yy_xyyy, g_y_0_yy_xyyyz, g_y_0_yy_xyyz, g_y_0_yy_xyyzz, g_y_0_yy_xyzz, g_y_0_yy_xyzzz, g_y_0_yy_xzzz, g_y_0_yy_xzzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyyz, g_y_0_yy_yyyz, g_y_0_yy_yyyzz, g_y_0_yy_yyzz, g_y_0_yy_yyzzz, g_y_0_yy_yzzz, g_y_0_yy_yzzzz, g_y_0_yy_zzzz, g_y_0_yy_zzzzz, g_y_0_yyz_xxxx, g_y_0_yyz_xxxy, g_y_0_yyz_xxxz, g_y_0_yyz_xxyy, g_y_0_yyz_xxyz, g_y_0_yyz_xxzz, g_y_0_yyz_xyyy, g_y_0_yyz_xyyz, g_y_0_yyz_xyzz, g_y_0_yyz_xzzz, g_y_0_yyz_yyyy, g_y_0_yyz_yyyz, g_y_0_yyz_yyzz, g_y_0_yyz_yzzz, g_y_0_yyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_xxxx[k] = -g_y_0_yy_xxxx[k] * cd_z[k] + g_y_0_yy_xxxxz[k];

            g_y_0_yyz_xxxy[k] = -g_y_0_yy_xxxy[k] * cd_z[k] + g_y_0_yy_xxxyz[k];

            g_y_0_yyz_xxxz[k] = -g_y_0_yy_xxxz[k] * cd_z[k] + g_y_0_yy_xxxzz[k];

            g_y_0_yyz_xxyy[k] = -g_y_0_yy_xxyy[k] * cd_z[k] + g_y_0_yy_xxyyz[k];

            g_y_0_yyz_xxyz[k] = -g_y_0_yy_xxyz[k] * cd_z[k] + g_y_0_yy_xxyzz[k];

            g_y_0_yyz_xxzz[k] = -g_y_0_yy_xxzz[k] * cd_z[k] + g_y_0_yy_xxzzz[k];

            g_y_0_yyz_xyyy[k] = -g_y_0_yy_xyyy[k] * cd_z[k] + g_y_0_yy_xyyyz[k];

            g_y_0_yyz_xyyz[k] = -g_y_0_yy_xyyz[k] * cd_z[k] + g_y_0_yy_xyyzz[k];

            g_y_0_yyz_xyzz[k] = -g_y_0_yy_xyzz[k] * cd_z[k] + g_y_0_yy_xyzzz[k];

            g_y_0_yyz_xzzz[k] = -g_y_0_yy_xzzz[k] * cd_z[k] + g_y_0_yy_xzzzz[k];

            g_y_0_yyz_yyyy[k] = -g_y_0_yy_yyyy[k] * cd_z[k] + g_y_0_yy_yyyyz[k];

            g_y_0_yyz_yyyz[k] = -g_y_0_yy_yyyz[k] * cd_z[k] + g_y_0_yy_yyyzz[k];

            g_y_0_yyz_yyzz[k] = -g_y_0_yy_yyzz[k] * cd_z[k] + g_y_0_yy_yyzzz[k];

            g_y_0_yyz_yzzz[k] = -g_y_0_yy_yzzz[k] * cd_z[k] + g_y_0_yy_yzzzz[k];

            g_y_0_yyz_zzzz[k] = -g_y_0_yy_zzzz[k] * cd_z[k] + g_y_0_yy_zzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 120);

        auto g_y_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 121);

        auto g_y_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 122);

        auto g_y_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 123);

        auto g_y_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 124);

        auto g_y_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 125);

        auto g_y_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 126);

        auto g_y_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 127);

        auto g_y_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 128);

        auto g_y_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 129);

        auto g_y_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 130);

        auto g_y_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 131);

        auto g_y_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 132);

        auto g_y_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 133);

        auto g_y_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 134);

        #pragma omp simd aligned(cd_z, g_y_0_yz_xxxx, g_y_0_yz_xxxxz, g_y_0_yz_xxxy, g_y_0_yz_xxxyz, g_y_0_yz_xxxz, g_y_0_yz_xxxzz, g_y_0_yz_xxyy, g_y_0_yz_xxyyz, g_y_0_yz_xxyz, g_y_0_yz_xxyzz, g_y_0_yz_xxzz, g_y_0_yz_xxzzz, g_y_0_yz_xyyy, g_y_0_yz_xyyyz, g_y_0_yz_xyyz, g_y_0_yz_xyyzz, g_y_0_yz_xyzz, g_y_0_yz_xyzzz, g_y_0_yz_xzzz, g_y_0_yz_xzzzz, g_y_0_yz_yyyy, g_y_0_yz_yyyyz, g_y_0_yz_yyyz, g_y_0_yz_yyyzz, g_y_0_yz_yyzz, g_y_0_yz_yyzzz, g_y_0_yz_yzzz, g_y_0_yz_yzzzz, g_y_0_yz_zzzz, g_y_0_yz_zzzzz, g_y_0_yzz_xxxx, g_y_0_yzz_xxxy, g_y_0_yzz_xxxz, g_y_0_yzz_xxyy, g_y_0_yzz_xxyz, g_y_0_yzz_xxzz, g_y_0_yzz_xyyy, g_y_0_yzz_xyyz, g_y_0_yzz_xyzz, g_y_0_yzz_xzzz, g_y_0_yzz_yyyy, g_y_0_yzz_yyyz, g_y_0_yzz_yyzz, g_y_0_yzz_yzzz, g_y_0_yzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_xxxx[k] = -g_y_0_yz_xxxx[k] * cd_z[k] + g_y_0_yz_xxxxz[k];

            g_y_0_yzz_xxxy[k] = -g_y_0_yz_xxxy[k] * cd_z[k] + g_y_0_yz_xxxyz[k];

            g_y_0_yzz_xxxz[k] = -g_y_0_yz_xxxz[k] * cd_z[k] + g_y_0_yz_xxxzz[k];

            g_y_0_yzz_xxyy[k] = -g_y_0_yz_xxyy[k] * cd_z[k] + g_y_0_yz_xxyyz[k];

            g_y_0_yzz_xxyz[k] = -g_y_0_yz_xxyz[k] * cd_z[k] + g_y_0_yz_xxyzz[k];

            g_y_0_yzz_xxzz[k] = -g_y_0_yz_xxzz[k] * cd_z[k] + g_y_0_yz_xxzzz[k];

            g_y_0_yzz_xyyy[k] = -g_y_0_yz_xyyy[k] * cd_z[k] + g_y_0_yz_xyyyz[k];

            g_y_0_yzz_xyyz[k] = -g_y_0_yz_xyyz[k] * cd_z[k] + g_y_0_yz_xyyzz[k];

            g_y_0_yzz_xyzz[k] = -g_y_0_yz_xyzz[k] * cd_z[k] + g_y_0_yz_xyzzz[k];

            g_y_0_yzz_xzzz[k] = -g_y_0_yz_xzzz[k] * cd_z[k] + g_y_0_yz_xzzzz[k];

            g_y_0_yzz_yyyy[k] = -g_y_0_yz_yyyy[k] * cd_z[k] + g_y_0_yz_yyyyz[k];

            g_y_0_yzz_yyyz[k] = -g_y_0_yz_yyyz[k] * cd_z[k] + g_y_0_yz_yyyzz[k];

            g_y_0_yzz_yyzz[k] = -g_y_0_yz_yyzz[k] * cd_z[k] + g_y_0_yz_yyzzz[k];

            g_y_0_yzz_yzzz[k] = -g_y_0_yz_yzzz[k] * cd_z[k] + g_y_0_yz_yzzzz[k];

            g_y_0_yzz_zzzz[k] = -g_y_0_yz_zzzz[k] * cd_z[k] + g_y_0_yz_zzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps  + 135);

        auto g_y_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 136);

        auto g_y_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 137);

        auto g_y_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 138);

        auto g_y_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 139);

        auto g_y_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 140);

        auto g_y_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 141);

        auto g_y_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 142);

        auto g_y_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 143);

        auto g_y_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 144);

        auto g_y_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps  + 145);

        auto g_y_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 146);

        auto g_y_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 147);

        auto g_y_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 148);

        auto g_y_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps  + 149);

        #pragma omp simd aligned(cd_z, g_y_0_zz_xxxx, g_y_0_zz_xxxxz, g_y_0_zz_xxxy, g_y_0_zz_xxxyz, g_y_0_zz_xxxz, g_y_0_zz_xxxzz, g_y_0_zz_xxyy, g_y_0_zz_xxyyz, g_y_0_zz_xxyz, g_y_0_zz_xxyzz, g_y_0_zz_xxzz, g_y_0_zz_xxzzz, g_y_0_zz_xyyy, g_y_0_zz_xyyyz, g_y_0_zz_xyyz, g_y_0_zz_xyyzz, g_y_0_zz_xyzz, g_y_0_zz_xyzzz, g_y_0_zz_xzzz, g_y_0_zz_xzzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyyz, g_y_0_zz_yyyz, g_y_0_zz_yyyzz, g_y_0_zz_yyzz, g_y_0_zz_yyzzz, g_y_0_zz_yzzz, g_y_0_zz_yzzzz, g_y_0_zz_zzzz, g_y_0_zz_zzzzz, g_y_0_zzz_xxxx, g_y_0_zzz_xxxy, g_y_0_zzz_xxxz, g_y_0_zzz_xxyy, g_y_0_zzz_xxyz, g_y_0_zzz_xxzz, g_y_0_zzz_xyyy, g_y_0_zzz_xyyz, g_y_0_zzz_xyzz, g_y_0_zzz_xzzz, g_y_0_zzz_yyyy, g_y_0_zzz_yyyz, g_y_0_zzz_yyzz, g_y_0_zzz_yzzz, g_y_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_xxxx[k] = -g_y_0_zz_xxxx[k] * cd_z[k] + g_y_0_zz_xxxxz[k];

            g_y_0_zzz_xxxy[k] = -g_y_0_zz_xxxy[k] * cd_z[k] + g_y_0_zz_xxxyz[k];

            g_y_0_zzz_xxxz[k] = -g_y_0_zz_xxxz[k] * cd_z[k] + g_y_0_zz_xxxzz[k];

            g_y_0_zzz_xxyy[k] = -g_y_0_zz_xxyy[k] * cd_z[k] + g_y_0_zz_xxyyz[k];

            g_y_0_zzz_xxyz[k] = -g_y_0_zz_xxyz[k] * cd_z[k] + g_y_0_zz_xxyzz[k];

            g_y_0_zzz_xxzz[k] = -g_y_0_zz_xxzz[k] * cd_z[k] + g_y_0_zz_xxzzz[k];

            g_y_0_zzz_xyyy[k] = -g_y_0_zz_xyyy[k] * cd_z[k] + g_y_0_zz_xyyyz[k];

            g_y_0_zzz_xyyz[k] = -g_y_0_zz_xyyz[k] * cd_z[k] + g_y_0_zz_xyyzz[k];

            g_y_0_zzz_xyzz[k] = -g_y_0_zz_xyzz[k] * cd_z[k] + g_y_0_zz_xyzzz[k];

            g_y_0_zzz_xzzz[k] = -g_y_0_zz_xzzz[k] * cd_z[k] + g_y_0_zz_xzzzz[k];

            g_y_0_zzz_yyyy[k] = -g_y_0_zz_yyyy[k] * cd_z[k] + g_y_0_zz_yyyyz[k];

            g_y_0_zzz_yyyz[k] = -g_y_0_zz_yyyz[k] * cd_z[k] + g_y_0_zz_yyyzz[k];

            g_y_0_zzz_yyzz[k] = -g_y_0_zz_yyzz[k] * cd_z[k] + g_y_0_zz_yyzzz[k];

            g_y_0_zzz_yzzz[k] = -g_y_0_zz_yzzz[k] * cd_z[k] + g_y_0_zz_yzzzz[k];

            g_y_0_zzz_zzzz[k] = -g_y_0_zz_zzzz[k] * cd_z[k] + g_y_0_zz_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 0);

        auto g_z_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 1);

        auto g_z_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 2);

        auto g_z_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 3);

        auto g_z_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 4);

        auto g_z_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 5);

        auto g_z_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 6);

        auto g_z_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 7);

        auto g_z_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 8);

        auto g_z_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 9);

        auto g_z_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 10);

        auto g_z_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 11);

        auto g_z_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 12);

        auto g_z_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 13);

        auto g_z_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_z_0_xx_xxxx, g_z_0_xx_xxxxx, g_z_0_xx_xxxxy, g_z_0_xx_xxxxz, g_z_0_xx_xxxy, g_z_0_xx_xxxyy, g_z_0_xx_xxxyz, g_z_0_xx_xxxz, g_z_0_xx_xxxzz, g_z_0_xx_xxyy, g_z_0_xx_xxyyy, g_z_0_xx_xxyyz, g_z_0_xx_xxyz, g_z_0_xx_xxyzz, g_z_0_xx_xxzz, g_z_0_xx_xxzzz, g_z_0_xx_xyyy, g_z_0_xx_xyyyy, g_z_0_xx_xyyyz, g_z_0_xx_xyyz, g_z_0_xx_xyyzz, g_z_0_xx_xyzz, g_z_0_xx_xyzzz, g_z_0_xx_xzzz, g_z_0_xx_xzzzz, g_z_0_xx_yyyy, g_z_0_xx_yyyz, g_z_0_xx_yyzz, g_z_0_xx_yzzz, g_z_0_xx_zzzz, g_z_0_xxx_xxxx, g_z_0_xxx_xxxy, g_z_0_xxx_xxxz, g_z_0_xxx_xxyy, g_z_0_xxx_xxyz, g_z_0_xxx_xxzz, g_z_0_xxx_xyyy, g_z_0_xxx_xyyz, g_z_0_xxx_xyzz, g_z_0_xxx_xzzz, g_z_0_xxx_yyyy, g_z_0_xxx_yyyz, g_z_0_xxx_yyzz, g_z_0_xxx_yzzz, g_z_0_xxx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_xxxx[k] = -g_z_0_xx_xxxx[k] * cd_x[k] + g_z_0_xx_xxxxx[k];

            g_z_0_xxx_xxxy[k] = -g_z_0_xx_xxxy[k] * cd_x[k] + g_z_0_xx_xxxxy[k];

            g_z_0_xxx_xxxz[k] = -g_z_0_xx_xxxz[k] * cd_x[k] + g_z_0_xx_xxxxz[k];

            g_z_0_xxx_xxyy[k] = -g_z_0_xx_xxyy[k] * cd_x[k] + g_z_0_xx_xxxyy[k];

            g_z_0_xxx_xxyz[k] = -g_z_0_xx_xxyz[k] * cd_x[k] + g_z_0_xx_xxxyz[k];

            g_z_0_xxx_xxzz[k] = -g_z_0_xx_xxzz[k] * cd_x[k] + g_z_0_xx_xxxzz[k];

            g_z_0_xxx_xyyy[k] = -g_z_0_xx_xyyy[k] * cd_x[k] + g_z_0_xx_xxyyy[k];

            g_z_0_xxx_xyyz[k] = -g_z_0_xx_xyyz[k] * cd_x[k] + g_z_0_xx_xxyyz[k];

            g_z_0_xxx_xyzz[k] = -g_z_0_xx_xyzz[k] * cd_x[k] + g_z_0_xx_xxyzz[k];

            g_z_0_xxx_xzzz[k] = -g_z_0_xx_xzzz[k] * cd_x[k] + g_z_0_xx_xxzzz[k];

            g_z_0_xxx_yyyy[k] = -g_z_0_xx_yyyy[k] * cd_x[k] + g_z_0_xx_xyyyy[k];

            g_z_0_xxx_yyyz[k] = -g_z_0_xx_yyyz[k] * cd_x[k] + g_z_0_xx_xyyyz[k];

            g_z_0_xxx_yyzz[k] = -g_z_0_xx_yyzz[k] * cd_x[k] + g_z_0_xx_xyyzz[k];

            g_z_0_xxx_yzzz[k] = -g_z_0_xx_yzzz[k] * cd_x[k] + g_z_0_xx_xyzzz[k];

            g_z_0_xxx_zzzz[k] = -g_z_0_xx_zzzz[k] * cd_x[k] + g_z_0_xx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 15);

        auto g_z_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 16);

        auto g_z_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 17);

        auto g_z_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 18);

        auto g_z_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 19);

        auto g_z_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 20);

        auto g_z_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 21);

        auto g_z_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 22);

        auto g_z_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 23);

        auto g_z_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 24);

        auto g_z_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 25);

        auto g_z_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 26);

        auto g_z_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 27);

        auto g_z_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 28);

        auto g_z_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_xxxx, g_z_0_xxy_xxxy, g_z_0_xxy_xxxz, g_z_0_xxy_xxyy, g_z_0_xxy_xxyz, g_z_0_xxy_xxzz, g_z_0_xxy_xyyy, g_z_0_xxy_xyyz, g_z_0_xxy_xyzz, g_z_0_xxy_xzzz, g_z_0_xxy_yyyy, g_z_0_xxy_yyyz, g_z_0_xxy_yyzz, g_z_0_xxy_yzzz, g_z_0_xxy_zzzz, g_z_0_xy_xxxx, g_z_0_xy_xxxxx, g_z_0_xy_xxxxy, g_z_0_xy_xxxxz, g_z_0_xy_xxxy, g_z_0_xy_xxxyy, g_z_0_xy_xxxyz, g_z_0_xy_xxxz, g_z_0_xy_xxxzz, g_z_0_xy_xxyy, g_z_0_xy_xxyyy, g_z_0_xy_xxyyz, g_z_0_xy_xxyz, g_z_0_xy_xxyzz, g_z_0_xy_xxzz, g_z_0_xy_xxzzz, g_z_0_xy_xyyy, g_z_0_xy_xyyyy, g_z_0_xy_xyyyz, g_z_0_xy_xyyz, g_z_0_xy_xyyzz, g_z_0_xy_xyzz, g_z_0_xy_xyzzz, g_z_0_xy_xzzz, g_z_0_xy_xzzzz, g_z_0_xy_yyyy, g_z_0_xy_yyyz, g_z_0_xy_yyzz, g_z_0_xy_yzzz, g_z_0_xy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_xxxx[k] = -g_z_0_xy_xxxx[k] * cd_x[k] + g_z_0_xy_xxxxx[k];

            g_z_0_xxy_xxxy[k] = -g_z_0_xy_xxxy[k] * cd_x[k] + g_z_0_xy_xxxxy[k];

            g_z_0_xxy_xxxz[k] = -g_z_0_xy_xxxz[k] * cd_x[k] + g_z_0_xy_xxxxz[k];

            g_z_0_xxy_xxyy[k] = -g_z_0_xy_xxyy[k] * cd_x[k] + g_z_0_xy_xxxyy[k];

            g_z_0_xxy_xxyz[k] = -g_z_0_xy_xxyz[k] * cd_x[k] + g_z_0_xy_xxxyz[k];

            g_z_0_xxy_xxzz[k] = -g_z_0_xy_xxzz[k] * cd_x[k] + g_z_0_xy_xxxzz[k];

            g_z_0_xxy_xyyy[k] = -g_z_0_xy_xyyy[k] * cd_x[k] + g_z_0_xy_xxyyy[k];

            g_z_0_xxy_xyyz[k] = -g_z_0_xy_xyyz[k] * cd_x[k] + g_z_0_xy_xxyyz[k];

            g_z_0_xxy_xyzz[k] = -g_z_0_xy_xyzz[k] * cd_x[k] + g_z_0_xy_xxyzz[k];

            g_z_0_xxy_xzzz[k] = -g_z_0_xy_xzzz[k] * cd_x[k] + g_z_0_xy_xxzzz[k];

            g_z_0_xxy_yyyy[k] = -g_z_0_xy_yyyy[k] * cd_x[k] + g_z_0_xy_xyyyy[k];

            g_z_0_xxy_yyyz[k] = -g_z_0_xy_yyyz[k] * cd_x[k] + g_z_0_xy_xyyyz[k];

            g_z_0_xxy_yyzz[k] = -g_z_0_xy_yyzz[k] * cd_x[k] + g_z_0_xy_xyyzz[k];

            g_z_0_xxy_yzzz[k] = -g_z_0_xy_yzzz[k] * cd_x[k] + g_z_0_xy_xyzzz[k];

            g_z_0_xxy_zzzz[k] = -g_z_0_xy_zzzz[k] * cd_x[k] + g_z_0_xy_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 30);

        auto g_z_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 31);

        auto g_z_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 32);

        auto g_z_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 33);

        auto g_z_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 34);

        auto g_z_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 35);

        auto g_z_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 36);

        auto g_z_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 37);

        auto g_z_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 38);

        auto g_z_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 39);

        auto g_z_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 40);

        auto g_z_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 41);

        auto g_z_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 42);

        auto g_z_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 43);

        auto g_z_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_xxxx, g_z_0_xxz_xxxy, g_z_0_xxz_xxxz, g_z_0_xxz_xxyy, g_z_0_xxz_xxyz, g_z_0_xxz_xxzz, g_z_0_xxz_xyyy, g_z_0_xxz_xyyz, g_z_0_xxz_xyzz, g_z_0_xxz_xzzz, g_z_0_xxz_yyyy, g_z_0_xxz_yyyz, g_z_0_xxz_yyzz, g_z_0_xxz_yzzz, g_z_0_xxz_zzzz, g_z_0_xz_xxxx, g_z_0_xz_xxxxx, g_z_0_xz_xxxxy, g_z_0_xz_xxxxz, g_z_0_xz_xxxy, g_z_0_xz_xxxyy, g_z_0_xz_xxxyz, g_z_0_xz_xxxz, g_z_0_xz_xxxzz, g_z_0_xz_xxyy, g_z_0_xz_xxyyy, g_z_0_xz_xxyyz, g_z_0_xz_xxyz, g_z_0_xz_xxyzz, g_z_0_xz_xxzz, g_z_0_xz_xxzzz, g_z_0_xz_xyyy, g_z_0_xz_xyyyy, g_z_0_xz_xyyyz, g_z_0_xz_xyyz, g_z_0_xz_xyyzz, g_z_0_xz_xyzz, g_z_0_xz_xyzzz, g_z_0_xz_xzzz, g_z_0_xz_xzzzz, g_z_0_xz_yyyy, g_z_0_xz_yyyz, g_z_0_xz_yyzz, g_z_0_xz_yzzz, g_z_0_xz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_xxxx[k] = -g_z_0_xz_xxxx[k] * cd_x[k] + g_z_0_xz_xxxxx[k];

            g_z_0_xxz_xxxy[k] = -g_z_0_xz_xxxy[k] * cd_x[k] + g_z_0_xz_xxxxy[k];

            g_z_0_xxz_xxxz[k] = -g_z_0_xz_xxxz[k] * cd_x[k] + g_z_0_xz_xxxxz[k];

            g_z_0_xxz_xxyy[k] = -g_z_0_xz_xxyy[k] * cd_x[k] + g_z_0_xz_xxxyy[k];

            g_z_0_xxz_xxyz[k] = -g_z_0_xz_xxyz[k] * cd_x[k] + g_z_0_xz_xxxyz[k];

            g_z_0_xxz_xxzz[k] = -g_z_0_xz_xxzz[k] * cd_x[k] + g_z_0_xz_xxxzz[k];

            g_z_0_xxz_xyyy[k] = -g_z_0_xz_xyyy[k] * cd_x[k] + g_z_0_xz_xxyyy[k];

            g_z_0_xxz_xyyz[k] = -g_z_0_xz_xyyz[k] * cd_x[k] + g_z_0_xz_xxyyz[k];

            g_z_0_xxz_xyzz[k] = -g_z_0_xz_xyzz[k] * cd_x[k] + g_z_0_xz_xxyzz[k];

            g_z_0_xxz_xzzz[k] = -g_z_0_xz_xzzz[k] * cd_x[k] + g_z_0_xz_xxzzz[k];

            g_z_0_xxz_yyyy[k] = -g_z_0_xz_yyyy[k] * cd_x[k] + g_z_0_xz_xyyyy[k];

            g_z_0_xxz_yyyz[k] = -g_z_0_xz_yyyz[k] * cd_x[k] + g_z_0_xz_xyyyz[k];

            g_z_0_xxz_yyzz[k] = -g_z_0_xz_yyzz[k] * cd_x[k] + g_z_0_xz_xyyzz[k];

            g_z_0_xxz_yzzz[k] = -g_z_0_xz_yzzz[k] * cd_x[k] + g_z_0_xz_xyzzz[k];

            g_z_0_xxz_zzzz[k] = -g_z_0_xz_zzzz[k] * cd_x[k] + g_z_0_xz_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 45);

        auto g_z_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 46);

        auto g_z_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 47);

        auto g_z_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 48);

        auto g_z_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 49);

        auto g_z_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 50);

        auto g_z_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 51);

        auto g_z_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 52);

        auto g_z_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 53);

        auto g_z_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 54);

        auto g_z_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 55);

        auto g_z_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 56);

        auto g_z_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 57);

        auto g_z_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 58);

        auto g_z_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_xxxx, g_z_0_xyy_xxxy, g_z_0_xyy_xxxz, g_z_0_xyy_xxyy, g_z_0_xyy_xxyz, g_z_0_xyy_xxzz, g_z_0_xyy_xyyy, g_z_0_xyy_xyyz, g_z_0_xyy_xyzz, g_z_0_xyy_xzzz, g_z_0_xyy_yyyy, g_z_0_xyy_yyyz, g_z_0_xyy_yyzz, g_z_0_xyy_yzzz, g_z_0_xyy_zzzz, g_z_0_yy_xxxx, g_z_0_yy_xxxxx, g_z_0_yy_xxxxy, g_z_0_yy_xxxxz, g_z_0_yy_xxxy, g_z_0_yy_xxxyy, g_z_0_yy_xxxyz, g_z_0_yy_xxxz, g_z_0_yy_xxxzz, g_z_0_yy_xxyy, g_z_0_yy_xxyyy, g_z_0_yy_xxyyz, g_z_0_yy_xxyz, g_z_0_yy_xxyzz, g_z_0_yy_xxzz, g_z_0_yy_xxzzz, g_z_0_yy_xyyy, g_z_0_yy_xyyyy, g_z_0_yy_xyyyz, g_z_0_yy_xyyz, g_z_0_yy_xyyzz, g_z_0_yy_xyzz, g_z_0_yy_xyzzz, g_z_0_yy_xzzz, g_z_0_yy_xzzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyzz, g_z_0_yy_yzzz, g_z_0_yy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_xxxx[k] = -g_z_0_yy_xxxx[k] * cd_x[k] + g_z_0_yy_xxxxx[k];

            g_z_0_xyy_xxxy[k] = -g_z_0_yy_xxxy[k] * cd_x[k] + g_z_0_yy_xxxxy[k];

            g_z_0_xyy_xxxz[k] = -g_z_0_yy_xxxz[k] * cd_x[k] + g_z_0_yy_xxxxz[k];

            g_z_0_xyy_xxyy[k] = -g_z_0_yy_xxyy[k] * cd_x[k] + g_z_0_yy_xxxyy[k];

            g_z_0_xyy_xxyz[k] = -g_z_0_yy_xxyz[k] * cd_x[k] + g_z_0_yy_xxxyz[k];

            g_z_0_xyy_xxzz[k] = -g_z_0_yy_xxzz[k] * cd_x[k] + g_z_0_yy_xxxzz[k];

            g_z_0_xyy_xyyy[k] = -g_z_0_yy_xyyy[k] * cd_x[k] + g_z_0_yy_xxyyy[k];

            g_z_0_xyy_xyyz[k] = -g_z_0_yy_xyyz[k] * cd_x[k] + g_z_0_yy_xxyyz[k];

            g_z_0_xyy_xyzz[k] = -g_z_0_yy_xyzz[k] * cd_x[k] + g_z_0_yy_xxyzz[k];

            g_z_0_xyy_xzzz[k] = -g_z_0_yy_xzzz[k] * cd_x[k] + g_z_0_yy_xxzzz[k];

            g_z_0_xyy_yyyy[k] = -g_z_0_yy_yyyy[k] * cd_x[k] + g_z_0_yy_xyyyy[k];

            g_z_0_xyy_yyyz[k] = -g_z_0_yy_yyyz[k] * cd_x[k] + g_z_0_yy_xyyyz[k];

            g_z_0_xyy_yyzz[k] = -g_z_0_yy_yyzz[k] * cd_x[k] + g_z_0_yy_xyyzz[k];

            g_z_0_xyy_yzzz[k] = -g_z_0_yy_yzzz[k] * cd_x[k] + g_z_0_yy_xyzzz[k];

            g_z_0_xyy_zzzz[k] = -g_z_0_yy_zzzz[k] * cd_x[k] + g_z_0_yy_xzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 60);

        auto g_z_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 61);

        auto g_z_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 62);

        auto g_z_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 63);

        auto g_z_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 64);

        auto g_z_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 65);

        auto g_z_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 66);

        auto g_z_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 67);

        auto g_z_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 68);

        auto g_z_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 69);

        auto g_z_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 70);

        auto g_z_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 71);

        auto g_z_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 72);

        auto g_z_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 73);

        auto g_z_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 74);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_xxxx, g_z_0_xyz_xxxy, g_z_0_xyz_xxxz, g_z_0_xyz_xxyy, g_z_0_xyz_xxyz, g_z_0_xyz_xxzz, g_z_0_xyz_xyyy, g_z_0_xyz_xyyz, g_z_0_xyz_xyzz, g_z_0_xyz_xzzz, g_z_0_xyz_yyyy, g_z_0_xyz_yyyz, g_z_0_xyz_yyzz, g_z_0_xyz_yzzz, g_z_0_xyz_zzzz, g_z_0_yz_xxxx, g_z_0_yz_xxxxx, g_z_0_yz_xxxxy, g_z_0_yz_xxxxz, g_z_0_yz_xxxy, g_z_0_yz_xxxyy, g_z_0_yz_xxxyz, g_z_0_yz_xxxz, g_z_0_yz_xxxzz, g_z_0_yz_xxyy, g_z_0_yz_xxyyy, g_z_0_yz_xxyyz, g_z_0_yz_xxyz, g_z_0_yz_xxyzz, g_z_0_yz_xxzz, g_z_0_yz_xxzzz, g_z_0_yz_xyyy, g_z_0_yz_xyyyy, g_z_0_yz_xyyyz, g_z_0_yz_xyyz, g_z_0_yz_xyyzz, g_z_0_yz_xyzz, g_z_0_yz_xyzzz, g_z_0_yz_xzzz, g_z_0_yz_xzzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyzz, g_z_0_yz_yzzz, g_z_0_yz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_xxxx[k] = -g_z_0_yz_xxxx[k] * cd_x[k] + g_z_0_yz_xxxxx[k];

            g_z_0_xyz_xxxy[k] = -g_z_0_yz_xxxy[k] * cd_x[k] + g_z_0_yz_xxxxy[k];

            g_z_0_xyz_xxxz[k] = -g_z_0_yz_xxxz[k] * cd_x[k] + g_z_0_yz_xxxxz[k];

            g_z_0_xyz_xxyy[k] = -g_z_0_yz_xxyy[k] * cd_x[k] + g_z_0_yz_xxxyy[k];

            g_z_0_xyz_xxyz[k] = -g_z_0_yz_xxyz[k] * cd_x[k] + g_z_0_yz_xxxyz[k];

            g_z_0_xyz_xxzz[k] = -g_z_0_yz_xxzz[k] * cd_x[k] + g_z_0_yz_xxxzz[k];

            g_z_0_xyz_xyyy[k] = -g_z_0_yz_xyyy[k] * cd_x[k] + g_z_0_yz_xxyyy[k];

            g_z_0_xyz_xyyz[k] = -g_z_0_yz_xyyz[k] * cd_x[k] + g_z_0_yz_xxyyz[k];

            g_z_0_xyz_xyzz[k] = -g_z_0_yz_xyzz[k] * cd_x[k] + g_z_0_yz_xxyzz[k];

            g_z_0_xyz_xzzz[k] = -g_z_0_yz_xzzz[k] * cd_x[k] + g_z_0_yz_xxzzz[k];

            g_z_0_xyz_yyyy[k] = -g_z_0_yz_yyyy[k] * cd_x[k] + g_z_0_yz_xyyyy[k];

            g_z_0_xyz_yyyz[k] = -g_z_0_yz_yyyz[k] * cd_x[k] + g_z_0_yz_xyyyz[k];

            g_z_0_xyz_yyzz[k] = -g_z_0_yz_yyzz[k] * cd_x[k] + g_z_0_yz_xyyzz[k];

            g_z_0_xyz_yzzz[k] = -g_z_0_yz_yzzz[k] * cd_x[k] + g_z_0_yz_xyzzz[k];

            g_z_0_xyz_zzzz[k] = -g_z_0_yz_zzzz[k] * cd_x[k] + g_z_0_yz_xzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 75);

        auto g_z_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 76);

        auto g_z_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 77);

        auto g_z_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 78);

        auto g_z_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 79);

        auto g_z_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 80);

        auto g_z_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 81);

        auto g_z_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 82);

        auto g_z_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 83);

        auto g_z_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 84);

        auto g_z_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 85);

        auto g_z_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 86);

        auto g_z_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 87);

        auto g_z_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 88);

        auto g_z_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 89);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_xxxx, g_z_0_xzz_xxxy, g_z_0_xzz_xxxz, g_z_0_xzz_xxyy, g_z_0_xzz_xxyz, g_z_0_xzz_xxzz, g_z_0_xzz_xyyy, g_z_0_xzz_xyyz, g_z_0_xzz_xyzz, g_z_0_xzz_xzzz, g_z_0_xzz_yyyy, g_z_0_xzz_yyyz, g_z_0_xzz_yyzz, g_z_0_xzz_yzzz, g_z_0_xzz_zzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxxx, g_z_0_zz_xxxxy, g_z_0_zz_xxxxz, g_z_0_zz_xxxy, g_z_0_zz_xxxyy, g_z_0_zz_xxxyz, g_z_0_zz_xxxz, g_z_0_zz_xxxzz, g_z_0_zz_xxyy, g_z_0_zz_xxyyy, g_z_0_zz_xxyyz, g_z_0_zz_xxyz, g_z_0_zz_xxyzz, g_z_0_zz_xxzz, g_z_0_zz_xxzzz, g_z_0_zz_xyyy, g_z_0_zz_xyyyy, g_z_0_zz_xyyyz, g_z_0_zz_xyyz, g_z_0_zz_xyyzz, g_z_0_zz_xyzz, g_z_0_zz_xyzzz, g_z_0_zz_xzzz, g_z_0_zz_xzzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_xxxx[k] = -g_z_0_zz_xxxx[k] * cd_x[k] + g_z_0_zz_xxxxx[k];

            g_z_0_xzz_xxxy[k] = -g_z_0_zz_xxxy[k] * cd_x[k] + g_z_0_zz_xxxxy[k];

            g_z_0_xzz_xxxz[k] = -g_z_0_zz_xxxz[k] * cd_x[k] + g_z_0_zz_xxxxz[k];

            g_z_0_xzz_xxyy[k] = -g_z_0_zz_xxyy[k] * cd_x[k] + g_z_0_zz_xxxyy[k];

            g_z_0_xzz_xxyz[k] = -g_z_0_zz_xxyz[k] * cd_x[k] + g_z_0_zz_xxxyz[k];

            g_z_0_xzz_xxzz[k] = -g_z_0_zz_xxzz[k] * cd_x[k] + g_z_0_zz_xxxzz[k];

            g_z_0_xzz_xyyy[k] = -g_z_0_zz_xyyy[k] * cd_x[k] + g_z_0_zz_xxyyy[k];

            g_z_0_xzz_xyyz[k] = -g_z_0_zz_xyyz[k] * cd_x[k] + g_z_0_zz_xxyyz[k];

            g_z_0_xzz_xyzz[k] = -g_z_0_zz_xyzz[k] * cd_x[k] + g_z_0_zz_xxyzz[k];

            g_z_0_xzz_xzzz[k] = -g_z_0_zz_xzzz[k] * cd_x[k] + g_z_0_zz_xxzzz[k];

            g_z_0_xzz_yyyy[k] = -g_z_0_zz_yyyy[k] * cd_x[k] + g_z_0_zz_xyyyy[k];

            g_z_0_xzz_yyyz[k] = -g_z_0_zz_yyyz[k] * cd_x[k] + g_z_0_zz_xyyyz[k];

            g_z_0_xzz_yyzz[k] = -g_z_0_zz_yyzz[k] * cd_x[k] + g_z_0_zz_xyyzz[k];

            g_z_0_xzz_yzzz[k] = -g_z_0_zz_yzzz[k] * cd_x[k] + g_z_0_zz_xyzzz[k];

            g_z_0_xzz_zzzz[k] = -g_z_0_zz_zzzz[k] * cd_x[k] + g_z_0_zz_xzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 90);

        auto g_z_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 91);

        auto g_z_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 92);

        auto g_z_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 93);

        auto g_z_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 94);

        auto g_z_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 95);

        auto g_z_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 96);

        auto g_z_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 97);

        auto g_z_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 98);

        auto g_z_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 99);

        auto g_z_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 100);

        auto g_z_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 101);

        auto g_z_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 102);

        auto g_z_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 103);

        auto g_z_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_z_0_yy_xxxx, g_z_0_yy_xxxxy, g_z_0_yy_xxxy, g_z_0_yy_xxxyy, g_z_0_yy_xxxyz, g_z_0_yy_xxxz, g_z_0_yy_xxyy, g_z_0_yy_xxyyy, g_z_0_yy_xxyyz, g_z_0_yy_xxyz, g_z_0_yy_xxyzz, g_z_0_yy_xxzz, g_z_0_yy_xyyy, g_z_0_yy_xyyyy, g_z_0_yy_xyyyz, g_z_0_yy_xyyz, g_z_0_yy_xyyzz, g_z_0_yy_xyzz, g_z_0_yy_xyzzz, g_z_0_yy_xzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyyy, g_z_0_yy_yyyyz, g_z_0_yy_yyyz, g_z_0_yy_yyyzz, g_z_0_yy_yyzz, g_z_0_yy_yyzzz, g_z_0_yy_yzzz, g_z_0_yy_yzzzz, g_z_0_yy_zzzz, g_z_0_yyy_xxxx, g_z_0_yyy_xxxy, g_z_0_yyy_xxxz, g_z_0_yyy_xxyy, g_z_0_yyy_xxyz, g_z_0_yyy_xxzz, g_z_0_yyy_xyyy, g_z_0_yyy_xyyz, g_z_0_yyy_xyzz, g_z_0_yyy_xzzz, g_z_0_yyy_yyyy, g_z_0_yyy_yyyz, g_z_0_yyy_yyzz, g_z_0_yyy_yzzz, g_z_0_yyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_xxxx[k] = -g_z_0_yy_xxxx[k] * cd_y[k] + g_z_0_yy_xxxxy[k];

            g_z_0_yyy_xxxy[k] = -g_z_0_yy_xxxy[k] * cd_y[k] + g_z_0_yy_xxxyy[k];

            g_z_0_yyy_xxxz[k] = -g_z_0_yy_xxxz[k] * cd_y[k] + g_z_0_yy_xxxyz[k];

            g_z_0_yyy_xxyy[k] = -g_z_0_yy_xxyy[k] * cd_y[k] + g_z_0_yy_xxyyy[k];

            g_z_0_yyy_xxyz[k] = -g_z_0_yy_xxyz[k] * cd_y[k] + g_z_0_yy_xxyyz[k];

            g_z_0_yyy_xxzz[k] = -g_z_0_yy_xxzz[k] * cd_y[k] + g_z_0_yy_xxyzz[k];

            g_z_0_yyy_xyyy[k] = -g_z_0_yy_xyyy[k] * cd_y[k] + g_z_0_yy_xyyyy[k];

            g_z_0_yyy_xyyz[k] = -g_z_0_yy_xyyz[k] * cd_y[k] + g_z_0_yy_xyyyz[k];

            g_z_0_yyy_xyzz[k] = -g_z_0_yy_xyzz[k] * cd_y[k] + g_z_0_yy_xyyzz[k];

            g_z_0_yyy_xzzz[k] = -g_z_0_yy_xzzz[k] * cd_y[k] + g_z_0_yy_xyzzz[k];

            g_z_0_yyy_yyyy[k] = -g_z_0_yy_yyyy[k] * cd_y[k] + g_z_0_yy_yyyyy[k];

            g_z_0_yyy_yyyz[k] = -g_z_0_yy_yyyz[k] * cd_y[k] + g_z_0_yy_yyyyz[k];

            g_z_0_yyy_yyzz[k] = -g_z_0_yy_yyzz[k] * cd_y[k] + g_z_0_yy_yyyzz[k];

            g_z_0_yyy_yzzz[k] = -g_z_0_yy_yzzz[k] * cd_y[k] + g_z_0_yy_yyzzz[k];

            g_z_0_yyy_zzzz[k] = -g_z_0_yy_zzzz[k] * cd_y[k] + g_z_0_yy_yzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 105);

        auto g_z_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 106);

        auto g_z_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 107);

        auto g_z_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 108);

        auto g_z_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 109);

        auto g_z_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 110);

        auto g_z_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 111);

        auto g_z_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 112);

        auto g_z_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 113);

        auto g_z_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 114);

        auto g_z_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 115);

        auto g_z_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 116);

        auto g_z_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 117);

        auto g_z_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 118);

        auto g_z_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 119);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_xxxx, g_z_0_yyz_xxxy, g_z_0_yyz_xxxz, g_z_0_yyz_xxyy, g_z_0_yyz_xxyz, g_z_0_yyz_xxzz, g_z_0_yyz_xyyy, g_z_0_yyz_xyyz, g_z_0_yyz_xyzz, g_z_0_yyz_xzzz, g_z_0_yyz_yyyy, g_z_0_yyz_yyyz, g_z_0_yyz_yyzz, g_z_0_yyz_yzzz, g_z_0_yyz_zzzz, g_z_0_yz_xxxx, g_z_0_yz_xxxxy, g_z_0_yz_xxxy, g_z_0_yz_xxxyy, g_z_0_yz_xxxyz, g_z_0_yz_xxxz, g_z_0_yz_xxyy, g_z_0_yz_xxyyy, g_z_0_yz_xxyyz, g_z_0_yz_xxyz, g_z_0_yz_xxyzz, g_z_0_yz_xxzz, g_z_0_yz_xyyy, g_z_0_yz_xyyyy, g_z_0_yz_xyyyz, g_z_0_yz_xyyz, g_z_0_yz_xyyzz, g_z_0_yz_xyzz, g_z_0_yz_xyzzz, g_z_0_yz_xzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyyy, g_z_0_yz_yyyyz, g_z_0_yz_yyyz, g_z_0_yz_yyyzz, g_z_0_yz_yyzz, g_z_0_yz_yyzzz, g_z_0_yz_yzzz, g_z_0_yz_yzzzz, g_z_0_yz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_xxxx[k] = -g_z_0_yz_xxxx[k] * cd_y[k] + g_z_0_yz_xxxxy[k];

            g_z_0_yyz_xxxy[k] = -g_z_0_yz_xxxy[k] * cd_y[k] + g_z_0_yz_xxxyy[k];

            g_z_0_yyz_xxxz[k] = -g_z_0_yz_xxxz[k] * cd_y[k] + g_z_0_yz_xxxyz[k];

            g_z_0_yyz_xxyy[k] = -g_z_0_yz_xxyy[k] * cd_y[k] + g_z_0_yz_xxyyy[k];

            g_z_0_yyz_xxyz[k] = -g_z_0_yz_xxyz[k] * cd_y[k] + g_z_0_yz_xxyyz[k];

            g_z_0_yyz_xxzz[k] = -g_z_0_yz_xxzz[k] * cd_y[k] + g_z_0_yz_xxyzz[k];

            g_z_0_yyz_xyyy[k] = -g_z_0_yz_xyyy[k] * cd_y[k] + g_z_0_yz_xyyyy[k];

            g_z_0_yyz_xyyz[k] = -g_z_0_yz_xyyz[k] * cd_y[k] + g_z_0_yz_xyyyz[k];

            g_z_0_yyz_xyzz[k] = -g_z_0_yz_xyzz[k] * cd_y[k] + g_z_0_yz_xyyzz[k];

            g_z_0_yyz_xzzz[k] = -g_z_0_yz_xzzz[k] * cd_y[k] + g_z_0_yz_xyzzz[k];

            g_z_0_yyz_yyyy[k] = -g_z_0_yz_yyyy[k] * cd_y[k] + g_z_0_yz_yyyyy[k];

            g_z_0_yyz_yyyz[k] = -g_z_0_yz_yyyz[k] * cd_y[k] + g_z_0_yz_yyyyz[k];

            g_z_0_yyz_yyzz[k] = -g_z_0_yz_yyzz[k] * cd_y[k] + g_z_0_yz_yyyzz[k];

            g_z_0_yyz_yzzz[k] = -g_z_0_yz_yzzz[k] * cd_y[k] + g_z_0_yz_yyzzz[k];

            g_z_0_yyz_zzzz[k] = -g_z_0_yz_zzzz[k] * cd_y[k] + g_z_0_yz_yzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 120);

        auto g_z_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 121);

        auto g_z_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 122);

        auto g_z_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 123);

        auto g_z_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 124);

        auto g_z_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 125);

        auto g_z_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 126);

        auto g_z_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 127);

        auto g_z_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 128);

        auto g_z_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 129);

        auto g_z_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 130);

        auto g_z_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 131);

        auto g_z_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 132);

        auto g_z_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 133);

        auto g_z_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 134);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_xxxx, g_z_0_yzz_xxxy, g_z_0_yzz_xxxz, g_z_0_yzz_xxyy, g_z_0_yzz_xxyz, g_z_0_yzz_xxzz, g_z_0_yzz_xyyy, g_z_0_yzz_xyyz, g_z_0_yzz_xyzz, g_z_0_yzz_xzzz, g_z_0_yzz_yyyy, g_z_0_yzz_yyyz, g_z_0_yzz_yyzz, g_z_0_yzz_yzzz, g_z_0_yzz_zzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxxy, g_z_0_zz_xxxy, g_z_0_zz_xxxyy, g_z_0_zz_xxxyz, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyyy, g_z_0_zz_xxyyz, g_z_0_zz_xxyz, g_z_0_zz_xxyzz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyyy, g_z_0_zz_xyyyz, g_z_0_zz_xyyz, g_z_0_zz_xyyzz, g_z_0_zz_xyzz, g_z_0_zz_xyzzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyyy, g_z_0_zz_yyyyz, g_z_0_zz_yyyz, g_z_0_zz_yyyzz, g_z_0_zz_yyzz, g_z_0_zz_yyzzz, g_z_0_zz_yzzz, g_z_0_zz_yzzzz, g_z_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_xxxx[k] = -g_z_0_zz_xxxx[k] * cd_y[k] + g_z_0_zz_xxxxy[k];

            g_z_0_yzz_xxxy[k] = -g_z_0_zz_xxxy[k] * cd_y[k] + g_z_0_zz_xxxyy[k];

            g_z_0_yzz_xxxz[k] = -g_z_0_zz_xxxz[k] * cd_y[k] + g_z_0_zz_xxxyz[k];

            g_z_0_yzz_xxyy[k] = -g_z_0_zz_xxyy[k] * cd_y[k] + g_z_0_zz_xxyyy[k];

            g_z_0_yzz_xxyz[k] = -g_z_0_zz_xxyz[k] * cd_y[k] + g_z_0_zz_xxyyz[k];

            g_z_0_yzz_xxzz[k] = -g_z_0_zz_xxzz[k] * cd_y[k] + g_z_0_zz_xxyzz[k];

            g_z_0_yzz_xyyy[k] = -g_z_0_zz_xyyy[k] * cd_y[k] + g_z_0_zz_xyyyy[k];

            g_z_0_yzz_xyyz[k] = -g_z_0_zz_xyyz[k] * cd_y[k] + g_z_0_zz_xyyyz[k];

            g_z_0_yzz_xyzz[k] = -g_z_0_zz_xyzz[k] * cd_y[k] + g_z_0_zz_xyyzz[k];

            g_z_0_yzz_xzzz[k] = -g_z_0_zz_xzzz[k] * cd_y[k] + g_z_0_zz_xyzzz[k];

            g_z_0_yzz_yyyy[k] = -g_z_0_zz_yyyy[k] * cd_y[k] + g_z_0_zz_yyyyy[k];

            g_z_0_yzz_yyyz[k] = -g_z_0_zz_yyyz[k] * cd_y[k] + g_z_0_zz_yyyyz[k];

            g_z_0_yzz_yyzz[k] = -g_z_0_zz_yyzz[k] * cd_y[k] + g_z_0_zz_yyyzz[k];

            g_z_0_yzz_yzzz[k] = -g_z_0_zz_yzzz[k] * cd_y[k] + g_z_0_zz_yyzzz[k];

            g_z_0_yzz_zzzz[k] = -g_z_0_zz_zzzz[k] * cd_y[k] + g_z_0_zz_yzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps  + 135);

        auto g_z_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 136);

        auto g_z_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 137);

        auto g_z_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 138);

        auto g_z_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 139);

        auto g_z_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 140);

        auto g_z_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 141);

        auto g_z_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 142);

        auto g_z_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 143);

        auto g_z_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 144);

        auto g_z_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps  + 145);

        auto g_z_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 146);

        auto g_z_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 147);

        auto g_z_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 148);

        auto g_z_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps  + 149);

        #pragma omp simd aligned(cd_z, g_z_0_zz_xxxx, g_z_0_zz_xxxxz, g_z_0_zz_xxxy, g_z_0_zz_xxxyz, g_z_0_zz_xxxz, g_z_0_zz_xxxzz, g_z_0_zz_xxyy, g_z_0_zz_xxyyz, g_z_0_zz_xxyz, g_z_0_zz_xxyzz, g_z_0_zz_xxzz, g_z_0_zz_xxzzz, g_z_0_zz_xyyy, g_z_0_zz_xyyyz, g_z_0_zz_xyyz, g_z_0_zz_xyyzz, g_z_0_zz_xyzz, g_z_0_zz_xyzzz, g_z_0_zz_xzzz, g_z_0_zz_xzzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyyz, g_z_0_zz_yyyz, g_z_0_zz_yyyzz, g_z_0_zz_yyzz, g_z_0_zz_yyzzz, g_z_0_zz_yzzz, g_z_0_zz_yzzzz, g_z_0_zz_zzzz, g_z_0_zz_zzzzz, g_z_0_zzz_xxxx, g_z_0_zzz_xxxy, g_z_0_zzz_xxxz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyz, g_z_0_zzz_xxzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyz, g_z_0_zzz_xyzz, g_z_0_zzz_xzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyz, g_z_0_zzz_yyzz, g_z_0_zzz_yzzz, g_z_0_zzz_zzzz, g_zz_xxxx, g_zz_xxxy, g_zz_xxxz, g_zz_xxyy, g_zz_xxyz, g_zz_xxzz, g_zz_xyyy, g_zz_xyyz, g_zz_xyzz, g_zz_xzzz, g_zz_yyyy, g_zz_yyyz, g_zz_yyzz, g_zz_yzzz, g_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_xxxx[k] = -g_zz_xxxx[k] - g_z_0_zz_xxxx[k] * cd_z[k] + g_z_0_zz_xxxxz[k];

            g_z_0_zzz_xxxy[k] = -g_zz_xxxy[k] - g_z_0_zz_xxxy[k] * cd_z[k] + g_z_0_zz_xxxyz[k];

            g_z_0_zzz_xxxz[k] = -g_zz_xxxz[k] - g_z_0_zz_xxxz[k] * cd_z[k] + g_z_0_zz_xxxzz[k];

            g_z_0_zzz_xxyy[k] = -g_zz_xxyy[k] - g_z_0_zz_xxyy[k] * cd_z[k] + g_z_0_zz_xxyyz[k];

            g_z_0_zzz_xxyz[k] = -g_zz_xxyz[k] - g_z_0_zz_xxyz[k] * cd_z[k] + g_z_0_zz_xxyzz[k];

            g_z_0_zzz_xxzz[k] = -g_zz_xxzz[k] - g_z_0_zz_xxzz[k] * cd_z[k] + g_z_0_zz_xxzzz[k];

            g_z_0_zzz_xyyy[k] = -g_zz_xyyy[k] - g_z_0_zz_xyyy[k] * cd_z[k] + g_z_0_zz_xyyyz[k];

            g_z_0_zzz_xyyz[k] = -g_zz_xyyz[k] - g_z_0_zz_xyyz[k] * cd_z[k] + g_z_0_zz_xyyzz[k];

            g_z_0_zzz_xyzz[k] = -g_zz_xyzz[k] - g_z_0_zz_xyzz[k] * cd_z[k] + g_z_0_zz_xyzzz[k];

            g_z_0_zzz_xzzz[k] = -g_zz_xzzz[k] - g_z_0_zz_xzzz[k] * cd_z[k] + g_z_0_zz_xzzzz[k];

            g_z_0_zzz_yyyy[k] = -g_zz_yyyy[k] - g_z_0_zz_yyyy[k] * cd_z[k] + g_z_0_zz_yyyyz[k];

            g_z_0_zzz_yyyz[k] = -g_zz_yyyz[k] - g_z_0_zz_yyyz[k] * cd_z[k] + g_z_0_zz_yyyzz[k];

            g_z_0_zzz_yyzz[k] = -g_zz_yyzz[k] - g_z_0_zz_yyzz[k] * cd_z[k] + g_z_0_zz_yyzzz[k];

            g_z_0_zzz_yzzz[k] = -g_zz_yzzz[k] - g_z_0_zz_yzzz[k] * cd_z[k] + g_z_0_zz_yzzzz[k];

            g_z_0_zzz_zzzz[k] = -g_zz_zzzz[k] - g_z_0_zz_zzzz[k] * cd_z[k] + g_z_0_zz_zzzzz[k];
        }
    }
}

} // t3ceri namespace

