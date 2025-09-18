#include "ElectronRepulsionGeom0010ContrRecXXHD.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhd(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhd,
                                            const size_t idx_xxgd,
                                            const size_t idx_geom_10_xxgd,
                                            const size_t idx_geom_10_xxgf,
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
            /// Set up components of auxilary buffer : SSGD

            const auto gd_off = idx_xxgd + (i * bcomps + j) * 90;

            auto g_xxxx_xx = cbuffer.data(gd_off + 0);

            auto g_xxxx_xy = cbuffer.data(gd_off + 1);

            auto g_xxxx_xz = cbuffer.data(gd_off + 2);

            auto g_xxxx_yy = cbuffer.data(gd_off + 3);

            auto g_xxxx_yz = cbuffer.data(gd_off + 4);

            auto g_xxxx_zz = cbuffer.data(gd_off + 5);

            auto g_xxxy_xx = cbuffer.data(gd_off + 6);

            auto g_xxxy_xy = cbuffer.data(gd_off + 7);

            auto g_xxxy_xz = cbuffer.data(gd_off + 8);

            auto g_xxxy_yy = cbuffer.data(gd_off + 9);

            auto g_xxxy_yz = cbuffer.data(gd_off + 10);

            auto g_xxxy_zz = cbuffer.data(gd_off + 11);

            auto g_xxxz_xx = cbuffer.data(gd_off + 12);

            auto g_xxxz_xy = cbuffer.data(gd_off + 13);

            auto g_xxxz_xz = cbuffer.data(gd_off + 14);

            auto g_xxxz_yy = cbuffer.data(gd_off + 15);

            auto g_xxxz_yz = cbuffer.data(gd_off + 16);

            auto g_xxxz_zz = cbuffer.data(gd_off + 17);

            auto g_xxyy_xx = cbuffer.data(gd_off + 18);

            auto g_xxyy_xy = cbuffer.data(gd_off + 19);

            auto g_xxyy_xz = cbuffer.data(gd_off + 20);

            auto g_xxyy_yy = cbuffer.data(gd_off + 21);

            auto g_xxyy_yz = cbuffer.data(gd_off + 22);

            auto g_xxyy_zz = cbuffer.data(gd_off + 23);

            auto g_xxyz_xx = cbuffer.data(gd_off + 24);

            auto g_xxyz_xy = cbuffer.data(gd_off + 25);

            auto g_xxyz_xz = cbuffer.data(gd_off + 26);

            auto g_xxyz_yy = cbuffer.data(gd_off + 27);

            auto g_xxyz_yz = cbuffer.data(gd_off + 28);

            auto g_xxyz_zz = cbuffer.data(gd_off + 29);

            auto g_xxzz_xx = cbuffer.data(gd_off + 30);

            auto g_xxzz_xy = cbuffer.data(gd_off + 31);

            auto g_xxzz_xz = cbuffer.data(gd_off + 32);

            auto g_xxzz_yy = cbuffer.data(gd_off + 33);

            auto g_xxzz_yz = cbuffer.data(gd_off + 34);

            auto g_xxzz_zz = cbuffer.data(gd_off + 35);

            auto g_xyyy_xx = cbuffer.data(gd_off + 36);

            auto g_xyyy_xy = cbuffer.data(gd_off + 37);

            auto g_xyyy_xz = cbuffer.data(gd_off + 38);

            auto g_xyyy_yy = cbuffer.data(gd_off + 39);

            auto g_xyyy_yz = cbuffer.data(gd_off + 40);

            auto g_xyyy_zz = cbuffer.data(gd_off + 41);

            auto g_xyyz_xx = cbuffer.data(gd_off + 42);

            auto g_xyyz_xy = cbuffer.data(gd_off + 43);

            auto g_xyyz_xz = cbuffer.data(gd_off + 44);

            auto g_xyyz_yy = cbuffer.data(gd_off + 45);

            auto g_xyyz_yz = cbuffer.data(gd_off + 46);

            auto g_xyyz_zz = cbuffer.data(gd_off + 47);

            auto g_xyzz_xx = cbuffer.data(gd_off + 48);

            auto g_xyzz_xy = cbuffer.data(gd_off + 49);

            auto g_xyzz_xz = cbuffer.data(gd_off + 50);

            auto g_xyzz_yy = cbuffer.data(gd_off + 51);

            auto g_xyzz_yz = cbuffer.data(gd_off + 52);

            auto g_xyzz_zz = cbuffer.data(gd_off + 53);

            auto g_xzzz_xx = cbuffer.data(gd_off + 54);

            auto g_xzzz_xy = cbuffer.data(gd_off + 55);

            auto g_xzzz_xz = cbuffer.data(gd_off + 56);

            auto g_xzzz_yy = cbuffer.data(gd_off + 57);

            auto g_xzzz_yz = cbuffer.data(gd_off + 58);

            auto g_xzzz_zz = cbuffer.data(gd_off + 59);

            auto g_yyyy_xx = cbuffer.data(gd_off + 60);

            auto g_yyyy_xy = cbuffer.data(gd_off + 61);

            auto g_yyyy_xz = cbuffer.data(gd_off + 62);

            auto g_yyyy_yy = cbuffer.data(gd_off + 63);

            auto g_yyyy_yz = cbuffer.data(gd_off + 64);

            auto g_yyyy_zz = cbuffer.data(gd_off + 65);

            auto g_yyyz_xx = cbuffer.data(gd_off + 66);

            auto g_yyyz_xy = cbuffer.data(gd_off + 67);

            auto g_yyyz_xz = cbuffer.data(gd_off + 68);

            auto g_yyyz_yy = cbuffer.data(gd_off + 69);

            auto g_yyyz_yz = cbuffer.data(gd_off + 70);

            auto g_yyyz_zz = cbuffer.data(gd_off + 71);

            auto g_yyzz_xx = cbuffer.data(gd_off + 72);

            auto g_yyzz_xy = cbuffer.data(gd_off + 73);

            auto g_yyzz_xz = cbuffer.data(gd_off + 74);

            auto g_yyzz_yy = cbuffer.data(gd_off + 75);

            auto g_yyzz_yz = cbuffer.data(gd_off + 76);

            auto g_yyzz_zz = cbuffer.data(gd_off + 77);

            auto g_yzzz_xx = cbuffer.data(gd_off + 78);

            auto g_yzzz_xy = cbuffer.data(gd_off + 79);

            auto g_yzzz_xz = cbuffer.data(gd_off + 80);

            auto g_yzzz_yy = cbuffer.data(gd_off + 81);

            auto g_yzzz_yz = cbuffer.data(gd_off + 82);

            auto g_yzzz_zz = cbuffer.data(gd_off + 83);

            auto g_zzzz_xx = cbuffer.data(gd_off + 84);

            auto g_zzzz_xy = cbuffer.data(gd_off + 85);

            auto g_zzzz_xz = cbuffer.data(gd_off + 86);

            auto g_zzzz_yy = cbuffer.data(gd_off + 87);

            auto g_zzzz_yz = cbuffer.data(gd_off + 88);

            auto g_zzzz_zz = cbuffer.data(gd_off + 89);

            /// Set up components of auxilary buffer : SSGD

            const auto gd_geom_10_off = idx_geom_10_xxgd + (i * bcomps + j) * 90;

            auto g_x_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_y_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 2);

            auto g_y_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 3);

            auto g_y_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 4);

            auto g_y_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 5);

            auto g_y_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_y_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 7);

            auto g_y_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 8);

            auto g_y_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 9);

            auto g_y_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 10);

            auto g_y_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 11);

            auto g_y_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_y_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 13);

            auto g_y_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 14);

            auto g_y_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 15);

            auto g_y_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 16);

            auto g_y_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 17);

            auto g_y_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_y_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 19);

            auto g_y_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 20);

            auto g_y_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 21);

            auto g_y_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 22);

            auto g_y_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 23);

            auto g_y_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_y_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 25);

            auto g_y_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 26);

            auto g_y_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 27);

            auto g_y_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 28);

            auto g_y_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 29);

            auto g_y_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_y_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_y_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 32);

            auto g_y_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 33);

            auto g_y_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 34);

            auto g_y_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 35);

            auto g_y_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_y_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_y_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 38);

            auto g_y_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 39);

            auto g_y_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 40);

            auto g_y_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 41);

            auto g_y_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_y_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_y_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 44);

            auto g_y_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 45);

            auto g_y_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 46);

            auto g_y_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 47);

            auto g_y_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 48);

            auto g_y_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 49);

            auto g_y_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 50);

            auto g_y_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 51);

            auto g_y_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 52);

            auto g_y_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 53);

            auto g_y_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 54);

            auto g_y_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 55);

            auto g_y_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 56);

            auto g_y_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 57);

            auto g_y_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 58);

            auto g_y_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 59);

            auto g_y_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 60);

            auto g_y_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 61);

            auto g_y_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 62);

            auto g_y_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 63);

            auto g_y_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 64);

            auto g_y_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 65);

            auto g_y_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 66);

            auto g_y_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 67);

            auto g_y_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 68);

            auto g_y_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 69);

            auto g_y_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 70);

            auto g_y_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 71);

            auto g_y_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 72);

            auto g_y_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 73);

            auto g_y_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 74);

            auto g_y_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 75);

            auto g_y_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 76);

            auto g_y_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 77);

            auto g_y_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 78);

            auto g_y_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 79);

            auto g_y_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 80);

            auto g_y_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 81);

            auto g_y_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 82);

            auto g_y_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 83);

            auto g_y_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 84);

            auto g_y_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 85);

            auto g_y_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 86);

            auto g_y_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 87);

            auto g_y_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 88);

            auto g_y_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps * bcomps + 89);

            auto g_z_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 2);

            auto g_z_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 3);

            auto g_z_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 4);

            auto g_z_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 5);

            auto g_z_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 6);

            auto g_z_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 7);

            auto g_z_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 8);

            auto g_z_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 9);

            auto g_z_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 10);

            auto g_z_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 11);

            auto g_z_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 12);

            auto g_z_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 13);

            auto g_z_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 14);

            auto g_z_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 15);

            auto g_z_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 16);

            auto g_z_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 17);

            auto g_z_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 18);

            auto g_z_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 19);

            auto g_z_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 20);

            auto g_z_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 21);

            auto g_z_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 22);

            auto g_z_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 23);

            auto g_z_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 24);

            auto g_z_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 25);

            auto g_z_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 26);

            auto g_z_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 27);

            auto g_z_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 28);

            auto g_z_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 29);

            auto g_z_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 30);

            auto g_z_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 31);

            auto g_z_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 32);

            auto g_z_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 33);

            auto g_z_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 34);

            auto g_z_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 35);

            auto g_z_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 36);

            auto g_z_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 37);

            auto g_z_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 38);

            auto g_z_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 39);

            auto g_z_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 40);

            auto g_z_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 41);

            auto g_z_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 42);

            auto g_z_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 43);

            auto g_z_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 44);

            auto g_z_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 45);

            auto g_z_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 46);

            auto g_z_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 47);

            auto g_z_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 48);

            auto g_z_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 49);

            auto g_z_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 50);

            auto g_z_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 51);

            auto g_z_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 52);

            auto g_z_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 53);

            auto g_z_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 54);

            auto g_z_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 55);

            auto g_z_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 56);

            auto g_z_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 57);

            auto g_z_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 58);

            auto g_z_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 59);

            auto g_z_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 60);

            auto g_z_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 61);

            auto g_z_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 62);

            auto g_z_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 63);

            auto g_z_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 64);

            auto g_z_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 65);

            auto g_z_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 66);

            auto g_z_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 67);

            auto g_z_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 68);

            auto g_z_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 69);

            auto g_z_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 70);

            auto g_z_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 71);

            auto g_z_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 72);

            auto g_z_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 73);

            auto g_z_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 74);

            auto g_z_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 75);

            auto g_z_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 76);

            auto g_z_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 77);

            auto g_z_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 78);

            auto g_z_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 79);

            auto g_z_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 80);

            auto g_z_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 81);

            auto g_z_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 82);

            auto g_z_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 83);

            auto g_z_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 84);

            auto g_z_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 85);

            auto g_z_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 86);

            auto g_z_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 87);

            auto g_z_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 88);

            auto g_z_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps * bcomps + 89);

            /// Set up components of auxilary buffer : SSGF

            const auto gf_geom_10_off = idx_geom_10_xxgf + (i * bcomps + j) * 150;

            auto g_x_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_y_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 5);

            auto g_y_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 6);

            auto g_y_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 7);

            auto g_y_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 8);

            auto g_y_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 9);

            auto g_y_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 10);

            auto g_y_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 11);

            auto g_y_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 12);

            auto g_y_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 13);

            auto g_y_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 14);

            auto g_y_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 15);

            auto g_y_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 16);

            auto g_y_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 17);

            auto g_y_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 18);

            auto g_y_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 19);

            auto g_y_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 20);

            auto g_y_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 21);

            auto g_y_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 22);

            auto g_y_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 23);

            auto g_y_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 24);

            auto g_y_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 25);

            auto g_y_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 26);

            auto g_y_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 27);

            auto g_y_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 28);

            auto g_y_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 29);

            auto g_y_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 30);

            auto g_y_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 31);

            auto g_y_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 32);

            auto g_y_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 33);

            auto g_y_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 34);

            auto g_y_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 35);

            auto g_y_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 36);

            auto g_y_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 37);

            auto g_y_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 38);

            auto g_y_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 39);

            auto g_y_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 40);

            auto g_y_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 41);

            auto g_y_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 42);

            auto g_y_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 43);

            auto g_y_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 44);

            auto g_y_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 45);

            auto g_y_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 46);

            auto g_y_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 47);

            auto g_y_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 48);

            auto g_y_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 49);

            auto g_y_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 50);

            auto g_y_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 51);

            auto g_y_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 52);

            auto g_y_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 53);

            auto g_y_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 54);

            auto g_y_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 55);

            auto g_y_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 56);

            auto g_y_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 57);

            auto g_y_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 58);

            auto g_y_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 59);

            auto g_y_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 60);

            auto g_y_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 61);

            auto g_y_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 62);

            auto g_y_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 63);

            auto g_y_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 64);

            auto g_y_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 65);

            auto g_y_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 66);

            auto g_y_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 67);

            auto g_y_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 68);

            auto g_y_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 69);

            auto g_y_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 70);

            auto g_y_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 71);

            auto g_y_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 72);

            auto g_y_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 73);

            auto g_y_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 74);

            auto g_y_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 75);

            auto g_y_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 76);

            auto g_y_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 77);

            auto g_y_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 78);

            auto g_y_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 79);

            auto g_y_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 80);

            auto g_y_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 81);

            auto g_y_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 82);

            auto g_y_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 83);

            auto g_y_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 84);

            auto g_y_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 85);

            auto g_y_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 86);

            auto g_y_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 87);

            auto g_y_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 88);

            auto g_y_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 89);

            auto g_y_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 90);

            auto g_y_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 91);

            auto g_y_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 92);

            auto g_y_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 93);

            auto g_y_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 94);

            auto g_y_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 95);

            auto g_y_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 96);

            auto g_y_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 97);

            auto g_y_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 98);

            auto g_y_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 99);

            auto g_y_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 100);

            auto g_y_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 101);

            auto g_y_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 102);

            auto g_y_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 103);

            auto g_y_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 104);

            auto g_y_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 105);

            auto g_y_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 106);

            auto g_y_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 107);

            auto g_y_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 108);

            auto g_y_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 109);

            auto g_y_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 110);

            auto g_y_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 111);

            auto g_y_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 112);

            auto g_y_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 113);

            auto g_y_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 114);

            auto g_y_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 115);

            auto g_y_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 116);

            auto g_y_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 117);

            auto g_y_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 118);

            auto g_y_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 119);

            auto g_y_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 120);

            auto g_y_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 121);

            auto g_y_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 122);

            auto g_y_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 123);

            auto g_y_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 124);

            auto g_y_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 125);

            auto g_y_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 126);

            auto g_y_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 127);

            auto g_y_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 128);

            auto g_y_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 129);

            auto g_y_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 130);

            auto g_y_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 131);

            auto g_y_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 132);

            auto g_y_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 133);

            auto g_y_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 134);

            auto g_y_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 135);

            auto g_y_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 136);

            auto g_y_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 137);

            auto g_y_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 138);

            auto g_y_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 139);

            auto g_y_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 140);

            auto g_y_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 141);

            auto g_y_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 142);

            auto g_y_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 143);

            auto g_y_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 144);

            auto g_y_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 145);

            auto g_y_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 146);

            auto g_y_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 147);

            auto g_y_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 148);

            auto g_y_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 150 * acomps * bcomps + 149);

            auto g_z_0_xxxx_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 5);

            auto g_z_0_xxxx_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 6);

            auto g_z_0_xxxx_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 7);

            auto g_z_0_xxxx_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 8);

            auto g_z_0_xxxx_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 9);

            auto g_z_0_xxxy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 10);

            auto g_z_0_xxxy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 11);

            auto g_z_0_xxxy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 12);

            auto g_z_0_xxxy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 13);

            auto g_z_0_xxxy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 14);

            auto g_z_0_xxxy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 15);

            auto g_z_0_xxxy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 16);

            auto g_z_0_xxxy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 17);

            auto g_z_0_xxxy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 18);

            auto g_z_0_xxxy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 19);

            auto g_z_0_xxxz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 20);

            auto g_z_0_xxxz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 21);

            auto g_z_0_xxxz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 22);

            auto g_z_0_xxxz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 23);

            auto g_z_0_xxxz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 24);

            auto g_z_0_xxxz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 25);

            auto g_z_0_xxxz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 26);

            auto g_z_0_xxxz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 27);

            auto g_z_0_xxxz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 28);

            auto g_z_0_xxxz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 29);

            auto g_z_0_xxyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 30);

            auto g_z_0_xxyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 31);

            auto g_z_0_xxyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 32);

            auto g_z_0_xxyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 33);

            auto g_z_0_xxyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 34);

            auto g_z_0_xxyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 35);

            auto g_z_0_xxyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 36);

            auto g_z_0_xxyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 37);

            auto g_z_0_xxyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 38);

            auto g_z_0_xxyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 39);

            auto g_z_0_xxyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 40);

            auto g_z_0_xxyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 41);

            auto g_z_0_xxyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 42);

            auto g_z_0_xxyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 43);

            auto g_z_0_xxyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 44);

            auto g_z_0_xxyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 45);

            auto g_z_0_xxyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 46);

            auto g_z_0_xxyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 47);

            auto g_z_0_xxyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 48);

            auto g_z_0_xxyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 49);

            auto g_z_0_xxzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 50);

            auto g_z_0_xxzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 51);

            auto g_z_0_xxzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 52);

            auto g_z_0_xxzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 53);

            auto g_z_0_xxzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 54);

            auto g_z_0_xxzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 55);

            auto g_z_0_xxzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 56);

            auto g_z_0_xxzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 57);

            auto g_z_0_xxzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 58);

            auto g_z_0_xxzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 59);

            auto g_z_0_xyyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 60);

            auto g_z_0_xyyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 61);

            auto g_z_0_xyyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 62);

            auto g_z_0_xyyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 63);

            auto g_z_0_xyyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 64);

            auto g_z_0_xyyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 65);

            auto g_z_0_xyyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 66);

            auto g_z_0_xyyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 67);

            auto g_z_0_xyyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 68);

            auto g_z_0_xyyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 69);

            auto g_z_0_xyyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 70);

            auto g_z_0_xyyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 71);

            auto g_z_0_xyyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 72);

            auto g_z_0_xyyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 73);

            auto g_z_0_xyyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 74);

            auto g_z_0_xyyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 75);

            auto g_z_0_xyyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 76);

            auto g_z_0_xyyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 77);

            auto g_z_0_xyyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 78);

            auto g_z_0_xyyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 79);

            auto g_z_0_xyzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 80);

            auto g_z_0_xyzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 81);

            auto g_z_0_xyzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 82);

            auto g_z_0_xyzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 83);

            auto g_z_0_xyzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 84);

            auto g_z_0_xyzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 85);

            auto g_z_0_xyzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 86);

            auto g_z_0_xyzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 87);

            auto g_z_0_xyzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 88);

            auto g_z_0_xyzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 89);

            auto g_z_0_xzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 90);

            auto g_z_0_xzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 91);

            auto g_z_0_xzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 92);

            auto g_z_0_xzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 93);

            auto g_z_0_xzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 94);

            auto g_z_0_xzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 95);

            auto g_z_0_xzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 96);

            auto g_z_0_xzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 97);

            auto g_z_0_xzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 98);

            auto g_z_0_xzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 99);

            auto g_z_0_yyyy_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 100);

            auto g_z_0_yyyy_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 101);

            auto g_z_0_yyyy_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 102);

            auto g_z_0_yyyy_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 103);

            auto g_z_0_yyyy_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 104);

            auto g_z_0_yyyy_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 105);

            auto g_z_0_yyyy_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 106);

            auto g_z_0_yyyy_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 107);

            auto g_z_0_yyyy_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 108);

            auto g_z_0_yyyy_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 109);

            auto g_z_0_yyyz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 110);

            auto g_z_0_yyyz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 111);

            auto g_z_0_yyyz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 112);

            auto g_z_0_yyyz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 113);

            auto g_z_0_yyyz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 114);

            auto g_z_0_yyyz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 115);

            auto g_z_0_yyyz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 116);

            auto g_z_0_yyyz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 117);

            auto g_z_0_yyyz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 118);

            auto g_z_0_yyyz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 119);

            auto g_z_0_yyzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 120);

            auto g_z_0_yyzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 121);

            auto g_z_0_yyzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 122);

            auto g_z_0_yyzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 123);

            auto g_z_0_yyzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 124);

            auto g_z_0_yyzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 125);

            auto g_z_0_yyzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 126);

            auto g_z_0_yyzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 127);

            auto g_z_0_yyzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 128);

            auto g_z_0_yyzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 129);

            auto g_z_0_yzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 130);

            auto g_z_0_yzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 131);

            auto g_z_0_yzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 132);

            auto g_z_0_yzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 133);

            auto g_z_0_yzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 134);

            auto g_z_0_yzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 135);

            auto g_z_0_yzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 136);

            auto g_z_0_yzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 137);

            auto g_z_0_yzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 138);

            auto g_z_0_yzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 139);

            auto g_z_0_zzzz_xxx = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 140);

            auto g_z_0_zzzz_xxy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 141);

            auto g_z_0_zzzz_xxz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 142);

            auto g_z_0_zzzz_xyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 143);

            auto g_z_0_zzzz_xyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 144);

            auto g_z_0_zzzz_xzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 145);

            auto g_z_0_zzzz_yyy = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 146);

            auto g_z_0_zzzz_yyz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 147);

            auto g_z_0_zzzz_yzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 148);

            auto g_z_0_zzzz_zzz = cbuffer.data(gf_geom_10_off + 300 * acomps * bcomps + 149);

            /// set up bra offset for contr_buffer_xxhd

            const auto hd_geom_10_off = idx_geom_10_xxhd + (i * bcomps + j) * 126;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_xx, g_x_0_xxxx_xxx, g_x_0_xxxx_xxy, g_x_0_xxxx_xxz, g_x_0_xxxx_xy, g_x_0_xxxx_xyy, g_x_0_xxxx_xyz, g_x_0_xxxx_xz, g_x_0_xxxx_xzz, g_x_0_xxxx_yy, g_x_0_xxxx_yz, g_x_0_xxxx_zz, g_x_0_xxxxx_xx, g_x_0_xxxxx_xy, g_x_0_xxxxx_xz, g_x_0_xxxxx_yy, g_x_0_xxxxx_yz, g_x_0_xxxxx_zz, g_xxxx_xx, g_xxxx_xy, g_xxxx_xz, g_xxxx_yy, g_xxxx_yz, g_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xx[k] = -g_xxxx_xx[k] - g_x_0_xxxx_xx[k] * cd_x[k] + g_x_0_xxxx_xxx[k];

                g_x_0_xxxxx_xy[k] = -g_xxxx_xy[k] - g_x_0_xxxx_xy[k] * cd_x[k] + g_x_0_xxxx_xxy[k];

                g_x_0_xxxxx_xz[k] = -g_xxxx_xz[k] - g_x_0_xxxx_xz[k] * cd_x[k] + g_x_0_xxxx_xxz[k];

                g_x_0_xxxxx_yy[k] = -g_xxxx_yy[k] - g_x_0_xxxx_yy[k] * cd_x[k] + g_x_0_xxxx_xyy[k];

                g_x_0_xxxxx_yz[k] = -g_xxxx_yz[k] - g_x_0_xxxx_yz[k] * cd_x[k] + g_x_0_xxxx_xyz[k];

                g_x_0_xxxxx_zz[k] = -g_xxxx_zz[k] - g_x_0_xxxx_zz[k] * cd_x[k] + g_x_0_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_xx, g_x_0_xxxx_xxy, g_x_0_xxxx_xy, g_x_0_xxxx_xyy, g_x_0_xxxx_xyz, g_x_0_xxxx_xz, g_x_0_xxxx_yy, g_x_0_xxxx_yyy, g_x_0_xxxx_yyz, g_x_0_xxxx_yz, g_x_0_xxxx_yzz, g_x_0_xxxx_zz, g_x_0_xxxxy_xx, g_x_0_xxxxy_xy, g_x_0_xxxxy_xz, g_x_0_xxxxy_yy, g_x_0_xxxxy_yz, g_x_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xx[k] = -g_x_0_xxxx_xx[k] * cd_y[k] + g_x_0_xxxx_xxy[k];

                g_x_0_xxxxy_xy[k] = -g_x_0_xxxx_xy[k] * cd_y[k] + g_x_0_xxxx_xyy[k];

                g_x_0_xxxxy_xz[k] = -g_x_0_xxxx_xz[k] * cd_y[k] + g_x_0_xxxx_xyz[k];

                g_x_0_xxxxy_yy[k] = -g_x_0_xxxx_yy[k] * cd_y[k] + g_x_0_xxxx_yyy[k];

                g_x_0_xxxxy_yz[k] = -g_x_0_xxxx_yz[k] * cd_y[k] + g_x_0_xxxx_yyz[k];

                g_x_0_xxxxy_zz[k] = -g_x_0_xxxx_zz[k] * cd_y[k] + g_x_0_xxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_xx, g_x_0_xxxx_xxz, g_x_0_xxxx_xy, g_x_0_xxxx_xyz, g_x_0_xxxx_xz, g_x_0_xxxx_xzz, g_x_0_xxxx_yy, g_x_0_xxxx_yyz, g_x_0_xxxx_yz, g_x_0_xxxx_yzz, g_x_0_xxxx_zz, g_x_0_xxxx_zzz, g_x_0_xxxxz_xx, g_x_0_xxxxz_xy, g_x_0_xxxxz_xz, g_x_0_xxxxz_yy, g_x_0_xxxxz_yz, g_x_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xx[k] = -g_x_0_xxxx_xx[k] * cd_z[k] + g_x_0_xxxx_xxz[k];

                g_x_0_xxxxz_xy[k] = -g_x_0_xxxx_xy[k] * cd_z[k] + g_x_0_xxxx_xyz[k];

                g_x_0_xxxxz_xz[k] = -g_x_0_xxxx_xz[k] * cd_z[k] + g_x_0_xxxx_xzz[k];

                g_x_0_xxxxz_yy[k] = -g_x_0_xxxx_yy[k] * cd_z[k] + g_x_0_xxxx_yyz[k];

                g_x_0_xxxxz_yz[k] = -g_x_0_xxxx_yz[k] * cd_z[k] + g_x_0_xxxx_yzz[k];

                g_x_0_xxxxz_zz[k] = -g_x_0_xxxx_zz[k] * cd_z[k] + g_x_0_xxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_xx, g_x_0_xxxy_xxy, g_x_0_xxxy_xy, g_x_0_xxxy_xyy, g_x_0_xxxy_xyz, g_x_0_xxxy_xz, g_x_0_xxxy_yy, g_x_0_xxxy_yyy, g_x_0_xxxy_yyz, g_x_0_xxxy_yz, g_x_0_xxxy_yzz, g_x_0_xxxy_zz, g_x_0_xxxyy_xx, g_x_0_xxxyy_xy, g_x_0_xxxyy_xz, g_x_0_xxxyy_yy, g_x_0_xxxyy_yz, g_x_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xx[k] = -g_x_0_xxxy_xx[k] * cd_y[k] + g_x_0_xxxy_xxy[k];

                g_x_0_xxxyy_xy[k] = -g_x_0_xxxy_xy[k] * cd_y[k] + g_x_0_xxxy_xyy[k];

                g_x_0_xxxyy_xz[k] = -g_x_0_xxxy_xz[k] * cd_y[k] + g_x_0_xxxy_xyz[k];

                g_x_0_xxxyy_yy[k] = -g_x_0_xxxy_yy[k] * cd_y[k] + g_x_0_xxxy_yyy[k];

                g_x_0_xxxyy_yz[k] = -g_x_0_xxxy_yz[k] * cd_y[k] + g_x_0_xxxy_yyz[k];

                g_x_0_xxxyy_zz[k] = -g_x_0_xxxy_zz[k] * cd_y[k] + g_x_0_xxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_xx, g_x_0_xxxyz_xy, g_x_0_xxxyz_xz, g_x_0_xxxyz_yy, g_x_0_xxxyz_yz, g_x_0_xxxyz_zz, g_x_0_xxxz_xx, g_x_0_xxxz_xxy, g_x_0_xxxz_xy, g_x_0_xxxz_xyy, g_x_0_xxxz_xyz, g_x_0_xxxz_xz, g_x_0_xxxz_yy, g_x_0_xxxz_yyy, g_x_0_xxxz_yyz, g_x_0_xxxz_yz, g_x_0_xxxz_yzz, g_x_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xx[k] = -g_x_0_xxxz_xx[k] * cd_y[k] + g_x_0_xxxz_xxy[k];

                g_x_0_xxxyz_xy[k] = -g_x_0_xxxz_xy[k] * cd_y[k] + g_x_0_xxxz_xyy[k];

                g_x_0_xxxyz_xz[k] = -g_x_0_xxxz_xz[k] * cd_y[k] + g_x_0_xxxz_xyz[k];

                g_x_0_xxxyz_yy[k] = -g_x_0_xxxz_yy[k] * cd_y[k] + g_x_0_xxxz_yyy[k];

                g_x_0_xxxyz_yz[k] = -g_x_0_xxxz_yz[k] * cd_y[k] + g_x_0_xxxz_yyz[k];

                g_x_0_xxxyz_zz[k] = -g_x_0_xxxz_zz[k] * cd_y[k] + g_x_0_xxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_xx, g_x_0_xxxz_xxz, g_x_0_xxxz_xy, g_x_0_xxxz_xyz, g_x_0_xxxz_xz, g_x_0_xxxz_xzz, g_x_0_xxxz_yy, g_x_0_xxxz_yyz, g_x_0_xxxz_yz, g_x_0_xxxz_yzz, g_x_0_xxxz_zz, g_x_0_xxxz_zzz, g_x_0_xxxzz_xx, g_x_0_xxxzz_xy, g_x_0_xxxzz_xz, g_x_0_xxxzz_yy, g_x_0_xxxzz_yz, g_x_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xx[k] = -g_x_0_xxxz_xx[k] * cd_z[k] + g_x_0_xxxz_xxz[k];

                g_x_0_xxxzz_xy[k] = -g_x_0_xxxz_xy[k] * cd_z[k] + g_x_0_xxxz_xyz[k];

                g_x_0_xxxzz_xz[k] = -g_x_0_xxxz_xz[k] * cd_z[k] + g_x_0_xxxz_xzz[k];

                g_x_0_xxxzz_yy[k] = -g_x_0_xxxz_yy[k] * cd_z[k] + g_x_0_xxxz_yyz[k];

                g_x_0_xxxzz_yz[k] = -g_x_0_xxxz_yz[k] * cd_z[k] + g_x_0_xxxz_yzz[k];

                g_x_0_xxxzz_zz[k] = -g_x_0_xxxz_zz[k] * cd_z[k] + g_x_0_xxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_xx, g_x_0_xxyy_xxy, g_x_0_xxyy_xy, g_x_0_xxyy_xyy, g_x_0_xxyy_xyz, g_x_0_xxyy_xz, g_x_0_xxyy_yy, g_x_0_xxyy_yyy, g_x_0_xxyy_yyz, g_x_0_xxyy_yz, g_x_0_xxyy_yzz, g_x_0_xxyy_zz, g_x_0_xxyyy_xx, g_x_0_xxyyy_xy, g_x_0_xxyyy_xz, g_x_0_xxyyy_yy, g_x_0_xxyyy_yz, g_x_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xx[k] = -g_x_0_xxyy_xx[k] * cd_y[k] + g_x_0_xxyy_xxy[k];

                g_x_0_xxyyy_xy[k] = -g_x_0_xxyy_xy[k] * cd_y[k] + g_x_0_xxyy_xyy[k];

                g_x_0_xxyyy_xz[k] = -g_x_0_xxyy_xz[k] * cd_y[k] + g_x_0_xxyy_xyz[k];

                g_x_0_xxyyy_yy[k] = -g_x_0_xxyy_yy[k] * cd_y[k] + g_x_0_xxyy_yyy[k];

                g_x_0_xxyyy_yz[k] = -g_x_0_xxyy_yz[k] * cd_y[k] + g_x_0_xxyy_yyz[k];

                g_x_0_xxyyy_zz[k] = -g_x_0_xxyy_zz[k] * cd_y[k] + g_x_0_xxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_xx, g_x_0_xxyyz_xy, g_x_0_xxyyz_xz, g_x_0_xxyyz_yy, g_x_0_xxyyz_yz, g_x_0_xxyyz_zz, g_x_0_xxyz_xx, g_x_0_xxyz_xxy, g_x_0_xxyz_xy, g_x_0_xxyz_xyy, g_x_0_xxyz_xyz, g_x_0_xxyz_xz, g_x_0_xxyz_yy, g_x_0_xxyz_yyy, g_x_0_xxyz_yyz, g_x_0_xxyz_yz, g_x_0_xxyz_yzz, g_x_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xx[k] = -g_x_0_xxyz_xx[k] * cd_y[k] + g_x_0_xxyz_xxy[k];

                g_x_0_xxyyz_xy[k] = -g_x_0_xxyz_xy[k] * cd_y[k] + g_x_0_xxyz_xyy[k];

                g_x_0_xxyyz_xz[k] = -g_x_0_xxyz_xz[k] * cd_y[k] + g_x_0_xxyz_xyz[k];

                g_x_0_xxyyz_yy[k] = -g_x_0_xxyz_yy[k] * cd_y[k] + g_x_0_xxyz_yyy[k];

                g_x_0_xxyyz_yz[k] = -g_x_0_xxyz_yz[k] * cd_y[k] + g_x_0_xxyz_yyz[k];

                g_x_0_xxyyz_zz[k] = -g_x_0_xxyz_zz[k] * cd_y[k] + g_x_0_xxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_xx, g_x_0_xxyzz_xy, g_x_0_xxyzz_xz, g_x_0_xxyzz_yy, g_x_0_xxyzz_yz, g_x_0_xxyzz_zz, g_x_0_xxzz_xx, g_x_0_xxzz_xxy, g_x_0_xxzz_xy, g_x_0_xxzz_xyy, g_x_0_xxzz_xyz, g_x_0_xxzz_xz, g_x_0_xxzz_yy, g_x_0_xxzz_yyy, g_x_0_xxzz_yyz, g_x_0_xxzz_yz, g_x_0_xxzz_yzz, g_x_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xx[k] = -g_x_0_xxzz_xx[k] * cd_y[k] + g_x_0_xxzz_xxy[k];

                g_x_0_xxyzz_xy[k] = -g_x_0_xxzz_xy[k] * cd_y[k] + g_x_0_xxzz_xyy[k];

                g_x_0_xxyzz_xz[k] = -g_x_0_xxzz_xz[k] * cd_y[k] + g_x_0_xxzz_xyz[k];

                g_x_0_xxyzz_yy[k] = -g_x_0_xxzz_yy[k] * cd_y[k] + g_x_0_xxzz_yyy[k];

                g_x_0_xxyzz_yz[k] = -g_x_0_xxzz_yz[k] * cd_y[k] + g_x_0_xxzz_yyz[k];

                g_x_0_xxyzz_zz[k] = -g_x_0_xxzz_zz[k] * cd_y[k] + g_x_0_xxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_xx, g_x_0_xxzz_xxz, g_x_0_xxzz_xy, g_x_0_xxzz_xyz, g_x_0_xxzz_xz, g_x_0_xxzz_xzz, g_x_0_xxzz_yy, g_x_0_xxzz_yyz, g_x_0_xxzz_yz, g_x_0_xxzz_yzz, g_x_0_xxzz_zz, g_x_0_xxzz_zzz, g_x_0_xxzzz_xx, g_x_0_xxzzz_xy, g_x_0_xxzzz_xz, g_x_0_xxzzz_yy, g_x_0_xxzzz_yz, g_x_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xx[k] = -g_x_0_xxzz_xx[k] * cd_z[k] + g_x_0_xxzz_xxz[k];

                g_x_0_xxzzz_xy[k] = -g_x_0_xxzz_xy[k] * cd_z[k] + g_x_0_xxzz_xyz[k];

                g_x_0_xxzzz_xz[k] = -g_x_0_xxzz_xz[k] * cd_z[k] + g_x_0_xxzz_xzz[k];

                g_x_0_xxzzz_yy[k] = -g_x_0_xxzz_yy[k] * cd_z[k] + g_x_0_xxzz_yyz[k];

                g_x_0_xxzzz_yz[k] = -g_x_0_xxzz_yz[k] * cd_z[k] + g_x_0_xxzz_yzz[k];

                g_x_0_xxzzz_zz[k] = -g_x_0_xxzz_zz[k] * cd_z[k] + g_x_0_xxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_xx, g_x_0_xyyy_xxy, g_x_0_xyyy_xy, g_x_0_xyyy_xyy, g_x_0_xyyy_xyz, g_x_0_xyyy_xz, g_x_0_xyyy_yy, g_x_0_xyyy_yyy, g_x_0_xyyy_yyz, g_x_0_xyyy_yz, g_x_0_xyyy_yzz, g_x_0_xyyy_zz, g_x_0_xyyyy_xx, g_x_0_xyyyy_xy, g_x_0_xyyyy_xz, g_x_0_xyyyy_yy, g_x_0_xyyyy_yz, g_x_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xx[k] = -g_x_0_xyyy_xx[k] * cd_y[k] + g_x_0_xyyy_xxy[k];

                g_x_0_xyyyy_xy[k] = -g_x_0_xyyy_xy[k] * cd_y[k] + g_x_0_xyyy_xyy[k];

                g_x_0_xyyyy_xz[k] = -g_x_0_xyyy_xz[k] * cd_y[k] + g_x_0_xyyy_xyz[k];

                g_x_0_xyyyy_yy[k] = -g_x_0_xyyy_yy[k] * cd_y[k] + g_x_0_xyyy_yyy[k];

                g_x_0_xyyyy_yz[k] = -g_x_0_xyyy_yz[k] * cd_y[k] + g_x_0_xyyy_yyz[k];

                g_x_0_xyyyy_zz[k] = -g_x_0_xyyy_zz[k] * cd_y[k] + g_x_0_xyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_xx, g_x_0_xyyyz_xy, g_x_0_xyyyz_xz, g_x_0_xyyyz_yy, g_x_0_xyyyz_yz, g_x_0_xyyyz_zz, g_x_0_xyyz_xx, g_x_0_xyyz_xxy, g_x_0_xyyz_xy, g_x_0_xyyz_xyy, g_x_0_xyyz_xyz, g_x_0_xyyz_xz, g_x_0_xyyz_yy, g_x_0_xyyz_yyy, g_x_0_xyyz_yyz, g_x_0_xyyz_yz, g_x_0_xyyz_yzz, g_x_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xx[k] = -g_x_0_xyyz_xx[k] * cd_y[k] + g_x_0_xyyz_xxy[k];

                g_x_0_xyyyz_xy[k] = -g_x_0_xyyz_xy[k] * cd_y[k] + g_x_0_xyyz_xyy[k];

                g_x_0_xyyyz_xz[k] = -g_x_0_xyyz_xz[k] * cd_y[k] + g_x_0_xyyz_xyz[k];

                g_x_0_xyyyz_yy[k] = -g_x_0_xyyz_yy[k] * cd_y[k] + g_x_0_xyyz_yyy[k];

                g_x_0_xyyyz_yz[k] = -g_x_0_xyyz_yz[k] * cd_y[k] + g_x_0_xyyz_yyz[k];

                g_x_0_xyyyz_zz[k] = -g_x_0_xyyz_zz[k] * cd_y[k] + g_x_0_xyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_xx, g_x_0_xyyzz_xy, g_x_0_xyyzz_xz, g_x_0_xyyzz_yy, g_x_0_xyyzz_yz, g_x_0_xyyzz_zz, g_x_0_xyzz_xx, g_x_0_xyzz_xxy, g_x_0_xyzz_xy, g_x_0_xyzz_xyy, g_x_0_xyzz_xyz, g_x_0_xyzz_xz, g_x_0_xyzz_yy, g_x_0_xyzz_yyy, g_x_0_xyzz_yyz, g_x_0_xyzz_yz, g_x_0_xyzz_yzz, g_x_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xx[k] = -g_x_0_xyzz_xx[k] * cd_y[k] + g_x_0_xyzz_xxy[k];

                g_x_0_xyyzz_xy[k] = -g_x_0_xyzz_xy[k] * cd_y[k] + g_x_0_xyzz_xyy[k];

                g_x_0_xyyzz_xz[k] = -g_x_0_xyzz_xz[k] * cd_y[k] + g_x_0_xyzz_xyz[k];

                g_x_0_xyyzz_yy[k] = -g_x_0_xyzz_yy[k] * cd_y[k] + g_x_0_xyzz_yyy[k];

                g_x_0_xyyzz_yz[k] = -g_x_0_xyzz_yz[k] * cd_y[k] + g_x_0_xyzz_yyz[k];

                g_x_0_xyyzz_zz[k] = -g_x_0_xyzz_zz[k] * cd_y[k] + g_x_0_xyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_xx, g_x_0_xyzzz_xy, g_x_0_xyzzz_xz, g_x_0_xyzzz_yy, g_x_0_xyzzz_yz, g_x_0_xyzzz_zz, g_x_0_xzzz_xx, g_x_0_xzzz_xxy, g_x_0_xzzz_xy, g_x_0_xzzz_xyy, g_x_0_xzzz_xyz, g_x_0_xzzz_xz, g_x_0_xzzz_yy, g_x_0_xzzz_yyy, g_x_0_xzzz_yyz, g_x_0_xzzz_yz, g_x_0_xzzz_yzz, g_x_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xx[k] = -g_x_0_xzzz_xx[k] * cd_y[k] + g_x_0_xzzz_xxy[k];

                g_x_0_xyzzz_xy[k] = -g_x_0_xzzz_xy[k] * cd_y[k] + g_x_0_xzzz_xyy[k];

                g_x_0_xyzzz_xz[k] = -g_x_0_xzzz_xz[k] * cd_y[k] + g_x_0_xzzz_xyz[k];

                g_x_0_xyzzz_yy[k] = -g_x_0_xzzz_yy[k] * cd_y[k] + g_x_0_xzzz_yyy[k];

                g_x_0_xyzzz_yz[k] = -g_x_0_xzzz_yz[k] * cd_y[k] + g_x_0_xzzz_yyz[k];

                g_x_0_xyzzz_zz[k] = -g_x_0_xzzz_zz[k] * cd_y[k] + g_x_0_xzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_xx, g_x_0_xzzz_xxz, g_x_0_xzzz_xy, g_x_0_xzzz_xyz, g_x_0_xzzz_xz, g_x_0_xzzz_xzz, g_x_0_xzzz_yy, g_x_0_xzzz_yyz, g_x_0_xzzz_yz, g_x_0_xzzz_yzz, g_x_0_xzzz_zz, g_x_0_xzzz_zzz, g_x_0_xzzzz_xx, g_x_0_xzzzz_xy, g_x_0_xzzzz_xz, g_x_0_xzzzz_yy, g_x_0_xzzzz_yz, g_x_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xx[k] = -g_x_0_xzzz_xx[k] * cd_z[k] + g_x_0_xzzz_xxz[k];

                g_x_0_xzzzz_xy[k] = -g_x_0_xzzz_xy[k] * cd_z[k] + g_x_0_xzzz_xyz[k];

                g_x_0_xzzzz_xz[k] = -g_x_0_xzzz_xz[k] * cd_z[k] + g_x_0_xzzz_xzz[k];

                g_x_0_xzzzz_yy[k] = -g_x_0_xzzz_yy[k] * cd_z[k] + g_x_0_xzzz_yyz[k];

                g_x_0_xzzzz_yz[k] = -g_x_0_xzzz_yz[k] * cd_z[k] + g_x_0_xzzz_yzz[k];

                g_x_0_xzzzz_zz[k] = -g_x_0_xzzz_zz[k] * cd_z[k] + g_x_0_xzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_xx, g_x_0_yyyy_xxy, g_x_0_yyyy_xy, g_x_0_yyyy_xyy, g_x_0_yyyy_xyz, g_x_0_yyyy_xz, g_x_0_yyyy_yy, g_x_0_yyyy_yyy, g_x_0_yyyy_yyz, g_x_0_yyyy_yz, g_x_0_yyyy_yzz, g_x_0_yyyy_zz, g_x_0_yyyyy_xx, g_x_0_yyyyy_xy, g_x_0_yyyyy_xz, g_x_0_yyyyy_yy, g_x_0_yyyyy_yz, g_x_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xx[k] = -g_x_0_yyyy_xx[k] * cd_y[k] + g_x_0_yyyy_xxy[k];

                g_x_0_yyyyy_xy[k] = -g_x_0_yyyy_xy[k] * cd_y[k] + g_x_0_yyyy_xyy[k];

                g_x_0_yyyyy_xz[k] = -g_x_0_yyyy_xz[k] * cd_y[k] + g_x_0_yyyy_xyz[k];

                g_x_0_yyyyy_yy[k] = -g_x_0_yyyy_yy[k] * cd_y[k] + g_x_0_yyyy_yyy[k];

                g_x_0_yyyyy_yz[k] = -g_x_0_yyyy_yz[k] * cd_y[k] + g_x_0_yyyy_yyz[k];

                g_x_0_yyyyy_zz[k] = -g_x_0_yyyy_zz[k] * cd_y[k] + g_x_0_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_xx, g_x_0_yyyyz_xy, g_x_0_yyyyz_xz, g_x_0_yyyyz_yy, g_x_0_yyyyz_yz, g_x_0_yyyyz_zz, g_x_0_yyyz_xx, g_x_0_yyyz_xxy, g_x_0_yyyz_xy, g_x_0_yyyz_xyy, g_x_0_yyyz_xyz, g_x_0_yyyz_xz, g_x_0_yyyz_yy, g_x_0_yyyz_yyy, g_x_0_yyyz_yyz, g_x_0_yyyz_yz, g_x_0_yyyz_yzz, g_x_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xx[k] = -g_x_0_yyyz_xx[k] * cd_y[k] + g_x_0_yyyz_xxy[k];

                g_x_0_yyyyz_xy[k] = -g_x_0_yyyz_xy[k] * cd_y[k] + g_x_0_yyyz_xyy[k];

                g_x_0_yyyyz_xz[k] = -g_x_0_yyyz_xz[k] * cd_y[k] + g_x_0_yyyz_xyz[k];

                g_x_0_yyyyz_yy[k] = -g_x_0_yyyz_yy[k] * cd_y[k] + g_x_0_yyyz_yyy[k];

                g_x_0_yyyyz_yz[k] = -g_x_0_yyyz_yz[k] * cd_y[k] + g_x_0_yyyz_yyz[k];

                g_x_0_yyyyz_zz[k] = -g_x_0_yyyz_zz[k] * cd_y[k] + g_x_0_yyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_xx, g_x_0_yyyzz_xy, g_x_0_yyyzz_xz, g_x_0_yyyzz_yy, g_x_0_yyyzz_yz, g_x_0_yyyzz_zz, g_x_0_yyzz_xx, g_x_0_yyzz_xxy, g_x_0_yyzz_xy, g_x_0_yyzz_xyy, g_x_0_yyzz_xyz, g_x_0_yyzz_xz, g_x_0_yyzz_yy, g_x_0_yyzz_yyy, g_x_0_yyzz_yyz, g_x_0_yyzz_yz, g_x_0_yyzz_yzz, g_x_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xx[k] = -g_x_0_yyzz_xx[k] * cd_y[k] + g_x_0_yyzz_xxy[k];

                g_x_0_yyyzz_xy[k] = -g_x_0_yyzz_xy[k] * cd_y[k] + g_x_0_yyzz_xyy[k];

                g_x_0_yyyzz_xz[k] = -g_x_0_yyzz_xz[k] * cd_y[k] + g_x_0_yyzz_xyz[k];

                g_x_0_yyyzz_yy[k] = -g_x_0_yyzz_yy[k] * cd_y[k] + g_x_0_yyzz_yyy[k];

                g_x_0_yyyzz_yz[k] = -g_x_0_yyzz_yz[k] * cd_y[k] + g_x_0_yyzz_yyz[k];

                g_x_0_yyyzz_zz[k] = -g_x_0_yyzz_zz[k] * cd_y[k] + g_x_0_yyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_xx, g_x_0_yyzzz_xy, g_x_0_yyzzz_xz, g_x_0_yyzzz_yy, g_x_0_yyzzz_yz, g_x_0_yyzzz_zz, g_x_0_yzzz_xx, g_x_0_yzzz_xxy, g_x_0_yzzz_xy, g_x_0_yzzz_xyy, g_x_0_yzzz_xyz, g_x_0_yzzz_xz, g_x_0_yzzz_yy, g_x_0_yzzz_yyy, g_x_0_yzzz_yyz, g_x_0_yzzz_yz, g_x_0_yzzz_yzz, g_x_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xx[k] = -g_x_0_yzzz_xx[k] * cd_y[k] + g_x_0_yzzz_xxy[k];

                g_x_0_yyzzz_xy[k] = -g_x_0_yzzz_xy[k] * cd_y[k] + g_x_0_yzzz_xyy[k];

                g_x_0_yyzzz_xz[k] = -g_x_0_yzzz_xz[k] * cd_y[k] + g_x_0_yzzz_xyz[k];

                g_x_0_yyzzz_yy[k] = -g_x_0_yzzz_yy[k] * cd_y[k] + g_x_0_yzzz_yyy[k];

                g_x_0_yyzzz_yz[k] = -g_x_0_yzzz_yz[k] * cd_y[k] + g_x_0_yzzz_yyz[k];

                g_x_0_yyzzz_zz[k] = -g_x_0_yzzz_zz[k] * cd_y[k] + g_x_0_yzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_xx, g_x_0_yzzzz_xy, g_x_0_yzzzz_xz, g_x_0_yzzzz_yy, g_x_0_yzzzz_yz, g_x_0_yzzzz_zz, g_x_0_zzzz_xx, g_x_0_zzzz_xxy, g_x_0_zzzz_xy, g_x_0_zzzz_xyy, g_x_0_zzzz_xyz, g_x_0_zzzz_xz, g_x_0_zzzz_yy, g_x_0_zzzz_yyy, g_x_0_zzzz_yyz, g_x_0_zzzz_yz, g_x_0_zzzz_yzz, g_x_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xx[k] = -g_x_0_zzzz_xx[k] * cd_y[k] + g_x_0_zzzz_xxy[k];

                g_x_0_yzzzz_xy[k] = -g_x_0_zzzz_xy[k] * cd_y[k] + g_x_0_zzzz_xyy[k];

                g_x_0_yzzzz_xz[k] = -g_x_0_zzzz_xz[k] * cd_y[k] + g_x_0_zzzz_xyz[k];

                g_x_0_yzzzz_yy[k] = -g_x_0_zzzz_yy[k] * cd_y[k] + g_x_0_zzzz_yyy[k];

                g_x_0_yzzzz_yz[k] = -g_x_0_zzzz_yz[k] * cd_y[k] + g_x_0_zzzz_yyz[k];

                g_x_0_yzzzz_zz[k] = -g_x_0_zzzz_zz[k] * cd_y[k] + g_x_0_zzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_xx, g_x_0_zzzz_xxz, g_x_0_zzzz_xy, g_x_0_zzzz_xyz, g_x_0_zzzz_xz, g_x_0_zzzz_xzz, g_x_0_zzzz_yy, g_x_0_zzzz_yyz, g_x_0_zzzz_yz, g_x_0_zzzz_yzz, g_x_0_zzzz_zz, g_x_0_zzzz_zzz, g_x_0_zzzzz_xx, g_x_0_zzzzz_xy, g_x_0_zzzzz_xz, g_x_0_zzzzz_yy, g_x_0_zzzzz_yz, g_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xx[k] = -g_x_0_zzzz_xx[k] * cd_z[k] + g_x_0_zzzz_xxz[k];

                g_x_0_zzzzz_xy[k] = -g_x_0_zzzz_xy[k] * cd_z[k] + g_x_0_zzzz_xyz[k];

                g_x_0_zzzzz_xz[k] = -g_x_0_zzzz_xz[k] * cd_z[k] + g_x_0_zzzz_xzz[k];

                g_x_0_zzzzz_yy[k] = -g_x_0_zzzz_yy[k] * cd_z[k] + g_x_0_zzzz_yyz[k];

                g_x_0_zzzzz_yz[k] = -g_x_0_zzzz_yz[k] * cd_z[k] + g_x_0_zzzz_yzz[k];

                g_x_0_zzzzz_zz[k] = -g_x_0_zzzz_zz[k] * cd_z[k] + g_x_0_zzzz_zzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_xx, g_y_0_xxxx_xxx, g_y_0_xxxx_xxy, g_y_0_xxxx_xxz, g_y_0_xxxx_xy, g_y_0_xxxx_xyy, g_y_0_xxxx_xyz, g_y_0_xxxx_xz, g_y_0_xxxx_xzz, g_y_0_xxxx_yy, g_y_0_xxxx_yz, g_y_0_xxxx_zz, g_y_0_xxxxx_xx, g_y_0_xxxxx_xy, g_y_0_xxxxx_xz, g_y_0_xxxxx_yy, g_y_0_xxxxx_yz, g_y_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xx[k] = -g_y_0_xxxx_xx[k] * cd_x[k] + g_y_0_xxxx_xxx[k];

                g_y_0_xxxxx_xy[k] = -g_y_0_xxxx_xy[k] * cd_x[k] + g_y_0_xxxx_xxy[k];

                g_y_0_xxxxx_xz[k] = -g_y_0_xxxx_xz[k] * cd_x[k] + g_y_0_xxxx_xxz[k];

                g_y_0_xxxxx_yy[k] = -g_y_0_xxxx_yy[k] * cd_x[k] + g_y_0_xxxx_xyy[k];

                g_y_0_xxxxx_yz[k] = -g_y_0_xxxx_yz[k] * cd_x[k] + g_y_0_xxxx_xyz[k];

                g_y_0_xxxxx_zz[k] = -g_y_0_xxxx_zz[k] * cd_x[k] + g_y_0_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 6);

            auto g_y_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 7);

            auto g_y_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 8);

            auto g_y_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 9);

            auto g_y_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 10);

            auto g_y_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_xx, g_y_0_xxxxy_xy, g_y_0_xxxxy_xz, g_y_0_xxxxy_yy, g_y_0_xxxxy_yz, g_y_0_xxxxy_zz, g_y_0_xxxy_xx, g_y_0_xxxy_xxx, g_y_0_xxxy_xxy, g_y_0_xxxy_xxz, g_y_0_xxxy_xy, g_y_0_xxxy_xyy, g_y_0_xxxy_xyz, g_y_0_xxxy_xz, g_y_0_xxxy_xzz, g_y_0_xxxy_yy, g_y_0_xxxy_yz, g_y_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xx[k] = -g_y_0_xxxy_xx[k] * cd_x[k] + g_y_0_xxxy_xxx[k];

                g_y_0_xxxxy_xy[k] = -g_y_0_xxxy_xy[k] * cd_x[k] + g_y_0_xxxy_xxy[k];

                g_y_0_xxxxy_xz[k] = -g_y_0_xxxy_xz[k] * cd_x[k] + g_y_0_xxxy_xxz[k];

                g_y_0_xxxxy_yy[k] = -g_y_0_xxxy_yy[k] * cd_x[k] + g_y_0_xxxy_xyy[k];

                g_y_0_xxxxy_yz[k] = -g_y_0_xxxy_yz[k] * cd_x[k] + g_y_0_xxxy_xyz[k];

                g_y_0_xxxxy_zz[k] = -g_y_0_xxxy_zz[k] * cd_x[k] + g_y_0_xxxy_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 12);

            auto g_y_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 13);

            auto g_y_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 14);

            auto g_y_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 15);

            auto g_y_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 16);

            auto g_y_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_xx, g_y_0_xxxxz_xy, g_y_0_xxxxz_xz, g_y_0_xxxxz_yy, g_y_0_xxxxz_yz, g_y_0_xxxxz_zz, g_y_0_xxxz_xx, g_y_0_xxxz_xxx, g_y_0_xxxz_xxy, g_y_0_xxxz_xxz, g_y_0_xxxz_xy, g_y_0_xxxz_xyy, g_y_0_xxxz_xyz, g_y_0_xxxz_xz, g_y_0_xxxz_xzz, g_y_0_xxxz_yy, g_y_0_xxxz_yz, g_y_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xx[k] = -g_y_0_xxxz_xx[k] * cd_x[k] + g_y_0_xxxz_xxx[k];

                g_y_0_xxxxz_xy[k] = -g_y_0_xxxz_xy[k] * cd_x[k] + g_y_0_xxxz_xxy[k];

                g_y_0_xxxxz_xz[k] = -g_y_0_xxxz_xz[k] * cd_x[k] + g_y_0_xxxz_xxz[k];

                g_y_0_xxxxz_yy[k] = -g_y_0_xxxz_yy[k] * cd_x[k] + g_y_0_xxxz_xyy[k];

                g_y_0_xxxxz_yz[k] = -g_y_0_xxxz_yz[k] * cd_x[k] + g_y_0_xxxz_xyz[k];

                g_y_0_xxxxz_zz[k] = -g_y_0_xxxz_zz[k] * cd_x[k] + g_y_0_xxxz_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 18);

            auto g_y_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 19);

            auto g_y_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 20);

            auto g_y_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 21);

            auto g_y_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 22);

            auto g_y_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_xx, g_y_0_xxxyy_xy, g_y_0_xxxyy_xz, g_y_0_xxxyy_yy, g_y_0_xxxyy_yz, g_y_0_xxxyy_zz, g_y_0_xxyy_xx, g_y_0_xxyy_xxx, g_y_0_xxyy_xxy, g_y_0_xxyy_xxz, g_y_0_xxyy_xy, g_y_0_xxyy_xyy, g_y_0_xxyy_xyz, g_y_0_xxyy_xz, g_y_0_xxyy_xzz, g_y_0_xxyy_yy, g_y_0_xxyy_yz, g_y_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xx[k] = -g_y_0_xxyy_xx[k] * cd_x[k] + g_y_0_xxyy_xxx[k];

                g_y_0_xxxyy_xy[k] = -g_y_0_xxyy_xy[k] * cd_x[k] + g_y_0_xxyy_xxy[k];

                g_y_0_xxxyy_xz[k] = -g_y_0_xxyy_xz[k] * cd_x[k] + g_y_0_xxyy_xxz[k];

                g_y_0_xxxyy_yy[k] = -g_y_0_xxyy_yy[k] * cd_x[k] + g_y_0_xxyy_xyy[k];

                g_y_0_xxxyy_yz[k] = -g_y_0_xxyy_yz[k] * cd_x[k] + g_y_0_xxyy_xyz[k];

                g_y_0_xxxyy_zz[k] = -g_y_0_xxyy_zz[k] * cd_x[k] + g_y_0_xxyy_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 24);

            auto g_y_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 25);

            auto g_y_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 26);

            auto g_y_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 27);

            auto g_y_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 28);

            auto g_y_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_xx, g_y_0_xxxyz_xy, g_y_0_xxxyz_xz, g_y_0_xxxyz_yy, g_y_0_xxxyz_yz, g_y_0_xxxyz_zz, g_y_0_xxyz_xx, g_y_0_xxyz_xxx, g_y_0_xxyz_xxy, g_y_0_xxyz_xxz, g_y_0_xxyz_xy, g_y_0_xxyz_xyy, g_y_0_xxyz_xyz, g_y_0_xxyz_xz, g_y_0_xxyz_xzz, g_y_0_xxyz_yy, g_y_0_xxyz_yz, g_y_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xx[k] = -g_y_0_xxyz_xx[k] * cd_x[k] + g_y_0_xxyz_xxx[k];

                g_y_0_xxxyz_xy[k] = -g_y_0_xxyz_xy[k] * cd_x[k] + g_y_0_xxyz_xxy[k];

                g_y_0_xxxyz_xz[k] = -g_y_0_xxyz_xz[k] * cd_x[k] + g_y_0_xxyz_xxz[k];

                g_y_0_xxxyz_yy[k] = -g_y_0_xxyz_yy[k] * cd_x[k] + g_y_0_xxyz_xyy[k];

                g_y_0_xxxyz_yz[k] = -g_y_0_xxyz_yz[k] * cd_x[k] + g_y_0_xxyz_xyz[k];

                g_y_0_xxxyz_zz[k] = -g_y_0_xxyz_zz[k] * cd_x[k] + g_y_0_xxyz_xzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 30);

            auto g_y_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 31);

            auto g_y_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 32);

            auto g_y_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 33);

            auto g_y_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 34);

            auto g_y_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_xx, g_y_0_xxxzz_xy, g_y_0_xxxzz_xz, g_y_0_xxxzz_yy, g_y_0_xxxzz_yz, g_y_0_xxxzz_zz, g_y_0_xxzz_xx, g_y_0_xxzz_xxx, g_y_0_xxzz_xxy, g_y_0_xxzz_xxz, g_y_0_xxzz_xy, g_y_0_xxzz_xyy, g_y_0_xxzz_xyz, g_y_0_xxzz_xz, g_y_0_xxzz_xzz, g_y_0_xxzz_yy, g_y_0_xxzz_yz, g_y_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xx[k] = -g_y_0_xxzz_xx[k] * cd_x[k] + g_y_0_xxzz_xxx[k];

                g_y_0_xxxzz_xy[k] = -g_y_0_xxzz_xy[k] * cd_x[k] + g_y_0_xxzz_xxy[k];

                g_y_0_xxxzz_xz[k] = -g_y_0_xxzz_xz[k] * cd_x[k] + g_y_0_xxzz_xxz[k];

                g_y_0_xxxzz_yy[k] = -g_y_0_xxzz_yy[k] * cd_x[k] + g_y_0_xxzz_xyy[k];

                g_y_0_xxxzz_yz[k] = -g_y_0_xxzz_yz[k] * cd_x[k] + g_y_0_xxzz_xyz[k];

                g_y_0_xxxzz_zz[k] = -g_y_0_xxzz_zz[k] * cd_x[k] + g_y_0_xxzz_xzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 36);

            auto g_y_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 37);

            auto g_y_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 38);

            auto g_y_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 39);

            auto g_y_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 40);

            auto g_y_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_xx, g_y_0_xxyyy_xy, g_y_0_xxyyy_xz, g_y_0_xxyyy_yy, g_y_0_xxyyy_yz, g_y_0_xxyyy_zz, g_y_0_xyyy_xx, g_y_0_xyyy_xxx, g_y_0_xyyy_xxy, g_y_0_xyyy_xxz, g_y_0_xyyy_xy, g_y_0_xyyy_xyy, g_y_0_xyyy_xyz, g_y_0_xyyy_xz, g_y_0_xyyy_xzz, g_y_0_xyyy_yy, g_y_0_xyyy_yz, g_y_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xx[k] = -g_y_0_xyyy_xx[k] * cd_x[k] + g_y_0_xyyy_xxx[k];

                g_y_0_xxyyy_xy[k] = -g_y_0_xyyy_xy[k] * cd_x[k] + g_y_0_xyyy_xxy[k];

                g_y_0_xxyyy_xz[k] = -g_y_0_xyyy_xz[k] * cd_x[k] + g_y_0_xyyy_xxz[k];

                g_y_0_xxyyy_yy[k] = -g_y_0_xyyy_yy[k] * cd_x[k] + g_y_0_xyyy_xyy[k];

                g_y_0_xxyyy_yz[k] = -g_y_0_xyyy_yz[k] * cd_x[k] + g_y_0_xyyy_xyz[k];

                g_y_0_xxyyy_zz[k] = -g_y_0_xyyy_zz[k] * cd_x[k] + g_y_0_xyyy_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 42);

            auto g_y_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 43);

            auto g_y_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 44);

            auto g_y_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 45);

            auto g_y_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 46);

            auto g_y_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_xx, g_y_0_xxyyz_xy, g_y_0_xxyyz_xz, g_y_0_xxyyz_yy, g_y_0_xxyyz_yz, g_y_0_xxyyz_zz, g_y_0_xyyz_xx, g_y_0_xyyz_xxx, g_y_0_xyyz_xxy, g_y_0_xyyz_xxz, g_y_0_xyyz_xy, g_y_0_xyyz_xyy, g_y_0_xyyz_xyz, g_y_0_xyyz_xz, g_y_0_xyyz_xzz, g_y_0_xyyz_yy, g_y_0_xyyz_yz, g_y_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xx[k] = -g_y_0_xyyz_xx[k] * cd_x[k] + g_y_0_xyyz_xxx[k];

                g_y_0_xxyyz_xy[k] = -g_y_0_xyyz_xy[k] * cd_x[k] + g_y_0_xyyz_xxy[k];

                g_y_0_xxyyz_xz[k] = -g_y_0_xyyz_xz[k] * cd_x[k] + g_y_0_xyyz_xxz[k];

                g_y_0_xxyyz_yy[k] = -g_y_0_xyyz_yy[k] * cd_x[k] + g_y_0_xyyz_xyy[k];

                g_y_0_xxyyz_yz[k] = -g_y_0_xyyz_yz[k] * cd_x[k] + g_y_0_xyyz_xyz[k];

                g_y_0_xxyyz_zz[k] = -g_y_0_xyyz_zz[k] * cd_x[k] + g_y_0_xyyz_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 48);

            auto g_y_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 49);

            auto g_y_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 50);

            auto g_y_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 51);

            auto g_y_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 52);

            auto g_y_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_xx, g_y_0_xxyzz_xy, g_y_0_xxyzz_xz, g_y_0_xxyzz_yy, g_y_0_xxyzz_yz, g_y_0_xxyzz_zz, g_y_0_xyzz_xx, g_y_0_xyzz_xxx, g_y_0_xyzz_xxy, g_y_0_xyzz_xxz, g_y_0_xyzz_xy, g_y_0_xyzz_xyy, g_y_0_xyzz_xyz, g_y_0_xyzz_xz, g_y_0_xyzz_xzz, g_y_0_xyzz_yy, g_y_0_xyzz_yz, g_y_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xx[k] = -g_y_0_xyzz_xx[k] * cd_x[k] + g_y_0_xyzz_xxx[k];

                g_y_0_xxyzz_xy[k] = -g_y_0_xyzz_xy[k] * cd_x[k] + g_y_0_xyzz_xxy[k];

                g_y_0_xxyzz_xz[k] = -g_y_0_xyzz_xz[k] * cd_x[k] + g_y_0_xyzz_xxz[k];

                g_y_0_xxyzz_yy[k] = -g_y_0_xyzz_yy[k] * cd_x[k] + g_y_0_xyzz_xyy[k];

                g_y_0_xxyzz_yz[k] = -g_y_0_xyzz_yz[k] * cd_x[k] + g_y_0_xyzz_xyz[k];

                g_y_0_xxyzz_zz[k] = -g_y_0_xyzz_zz[k] * cd_x[k] + g_y_0_xyzz_xzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 54);

            auto g_y_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 55);

            auto g_y_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 56);

            auto g_y_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 57);

            auto g_y_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 58);

            auto g_y_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_xx, g_y_0_xxzzz_xy, g_y_0_xxzzz_xz, g_y_0_xxzzz_yy, g_y_0_xxzzz_yz, g_y_0_xxzzz_zz, g_y_0_xzzz_xx, g_y_0_xzzz_xxx, g_y_0_xzzz_xxy, g_y_0_xzzz_xxz, g_y_0_xzzz_xy, g_y_0_xzzz_xyy, g_y_0_xzzz_xyz, g_y_0_xzzz_xz, g_y_0_xzzz_xzz, g_y_0_xzzz_yy, g_y_0_xzzz_yz, g_y_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xx[k] = -g_y_0_xzzz_xx[k] * cd_x[k] + g_y_0_xzzz_xxx[k];

                g_y_0_xxzzz_xy[k] = -g_y_0_xzzz_xy[k] * cd_x[k] + g_y_0_xzzz_xxy[k];

                g_y_0_xxzzz_xz[k] = -g_y_0_xzzz_xz[k] * cd_x[k] + g_y_0_xzzz_xxz[k];

                g_y_0_xxzzz_yy[k] = -g_y_0_xzzz_yy[k] * cd_x[k] + g_y_0_xzzz_xyy[k];

                g_y_0_xxzzz_yz[k] = -g_y_0_xzzz_yz[k] * cd_x[k] + g_y_0_xzzz_xyz[k];

                g_y_0_xxzzz_zz[k] = -g_y_0_xzzz_zz[k] * cd_x[k] + g_y_0_xzzz_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 60);

            auto g_y_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 61);

            auto g_y_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 62);

            auto g_y_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 63);

            auto g_y_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 64);

            auto g_y_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_xx, g_y_0_xyyyy_xy, g_y_0_xyyyy_xz, g_y_0_xyyyy_yy, g_y_0_xyyyy_yz, g_y_0_xyyyy_zz, g_y_0_yyyy_xx, g_y_0_yyyy_xxx, g_y_0_yyyy_xxy, g_y_0_yyyy_xxz, g_y_0_yyyy_xy, g_y_0_yyyy_xyy, g_y_0_yyyy_xyz, g_y_0_yyyy_xz, g_y_0_yyyy_xzz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xx[k] = -g_y_0_yyyy_xx[k] * cd_x[k] + g_y_0_yyyy_xxx[k];

                g_y_0_xyyyy_xy[k] = -g_y_0_yyyy_xy[k] * cd_x[k] + g_y_0_yyyy_xxy[k];

                g_y_0_xyyyy_xz[k] = -g_y_0_yyyy_xz[k] * cd_x[k] + g_y_0_yyyy_xxz[k];

                g_y_0_xyyyy_yy[k] = -g_y_0_yyyy_yy[k] * cd_x[k] + g_y_0_yyyy_xyy[k];

                g_y_0_xyyyy_yz[k] = -g_y_0_yyyy_yz[k] * cd_x[k] + g_y_0_yyyy_xyz[k];

                g_y_0_xyyyy_zz[k] = -g_y_0_yyyy_zz[k] * cd_x[k] + g_y_0_yyyy_xzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 66);

            auto g_y_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 67);

            auto g_y_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 68);

            auto g_y_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 69);

            auto g_y_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 70);

            auto g_y_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_xx, g_y_0_xyyyz_xy, g_y_0_xyyyz_xz, g_y_0_xyyyz_yy, g_y_0_xyyyz_yz, g_y_0_xyyyz_zz, g_y_0_yyyz_xx, g_y_0_yyyz_xxx, g_y_0_yyyz_xxy, g_y_0_yyyz_xxz, g_y_0_yyyz_xy, g_y_0_yyyz_xyy, g_y_0_yyyz_xyz, g_y_0_yyyz_xz, g_y_0_yyyz_xzz, g_y_0_yyyz_yy, g_y_0_yyyz_yz, g_y_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xx[k] = -g_y_0_yyyz_xx[k] * cd_x[k] + g_y_0_yyyz_xxx[k];

                g_y_0_xyyyz_xy[k] = -g_y_0_yyyz_xy[k] * cd_x[k] + g_y_0_yyyz_xxy[k];

                g_y_0_xyyyz_xz[k] = -g_y_0_yyyz_xz[k] * cd_x[k] + g_y_0_yyyz_xxz[k];

                g_y_0_xyyyz_yy[k] = -g_y_0_yyyz_yy[k] * cd_x[k] + g_y_0_yyyz_xyy[k];

                g_y_0_xyyyz_yz[k] = -g_y_0_yyyz_yz[k] * cd_x[k] + g_y_0_yyyz_xyz[k];

                g_y_0_xyyyz_zz[k] = -g_y_0_yyyz_zz[k] * cd_x[k] + g_y_0_yyyz_xzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 72);

            auto g_y_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 73);

            auto g_y_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 74);

            auto g_y_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 75);

            auto g_y_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 76);

            auto g_y_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_xx, g_y_0_xyyzz_xy, g_y_0_xyyzz_xz, g_y_0_xyyzz_yy, g_y_0_xyyzz_yz, g_y_0_xyyzz_zz, g_y_0_yyzz_xx, g_y_0_yyzz_xxx, g_y_0_yyzz_xxy, g_y_0_yyzz_xxz, g_y_0_yyzz_xy, g_y_0_yyzz_xyy, g_y_0_yyzz_xyz, g_y_0_yyzz_xz, g_y_0_yyzz_xzz, g_y_0_yyzz_yy, g_y_0_yyzz_yz, g_y_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xx[k] = -g_y_0_yyzz_xx[k] * cd_x[k] + g_y_0_yyzz_xxx[k];

                g_y_0_xyyzz_xy[k] = -g_y_0_yyzz_xy[k] * cd_x[k] + g_y_0_yyzz_xxy[k];

                g_y_0_xyyzz_xz[k] = -g_y_0_yyzz_xz[k] * cd_x[k] + g_y_0_yyzz_xxz[k];

                g_y_0_xyyzz_yy[k] = -g_y_0_yyzz_yy[k] * cd_x[k] + g_y_0_yyzz_xyy[k];

                g_y_0_xyyzz_yz[k] = -g_y_0_yyzz_yz[k] * cd_x[k] + g_y_0_yyzz_xyz[k];

                g_y_0_xyyzz_zz[k] = -g_y_0_yyzz_zz[k] * cd_x[k] + g_y_0_yyzz_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 78);

            auto g_y_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 79);

            auto g_y_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 80);

            auto g_y_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 81);

            auto g_y_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 82);

            auto g_y_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_xx, g_y_0_xyzzz_xy, g_y_0_xyzzz_xz, g_y_0_xyzzz_yy, g_y_0_xyzzz_yz, g_y_0_xyzzz_zz, g_y_0_yzzz_xx, g_y_0_yzzz_xxx, g_y_0_yzzz_xxy, g_y_0_yzzz_xxz, g_y_0_yzzz_xy, g_y_0_yzzz_xyy, g_y_0_yzzz_xyz, g_y_0_yzzz_xz, g_y_0_yzzz_xzz, g_y_0_yzzz_yy, g_y_0_yzzz_yz, g_y_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xx[k] = -g_y_0_yzzz_xx[k] * cd_x[k] + g_y_0_yzzz_xxx[k];

                g_y_0_xyzzz_xy[k] = -g_y_0_yzzz_xy[k] * cd_x[k] + g_y_0_yzzz_xxy[k];

                g_y_0_xyzzz_xz[k] = -g_y_0_yzzz_xz[k] * cd_x[k] + g_y_0_yzzz_xxz[k];

                g_y_0_xyzzz_yy[k] = -g_y_0_yzzz_yy[k] * cd_x[k] + g_y_0_yzzz_xyy[k];

                g_y_0_xyzzz_yz[k] = -g_y_0_yzzz_yz[k] * cd_x[k] + g_y_0_yzzz_xyz[k];

                g_y_0_xyzzz_zz[k] = -g_y_0_yzzz_zz[k] * cd_x[k] + g_y_0_yzzz_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 84);

            auto g_y_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 85);

            auto g_y_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 86);

            auto g_y_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 87);

            auto g_y_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 88);

            auto g_y_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_xx, g_y_0_xzzzz_xy, g_y_0_xzzzz_xz, g_y_0_xzzzz_yy, g_y_0_xzzzz_yz, g_y_0_xzzzz_zz, g_y_0_zzzz_xx, g_y_0_zzzz_xxx, g_y_0_zzzz_xxy, g_y_0_zzzz_xxz, g_y_0_zzzz_xy, g_y_0_zzzz_xyy, g_y_0_zzzz_xyz, g_y_0_zzzz_xz, g_y_0_zzzz_xzz, g_y_0_zzzz_yy, g_y_0_zzzz_yz, g_y_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xx[k] = -g_y_0_zzzz_xx[k] * cd_x[k] + g_y_0_zzzz_xxx[k];

                g_y_0_xzzzz_xy[k] = -g_y_0_zzzz_xy[k] * cd_x[k] + g_y_0_zzzz_xxy[k];

                g_y_0_xzzzz_xz[k] = -g_y_0_zzzz_xz[k] * cd_x[k] + g_y_0_zzzz_xxz[k];

                g_y_0_xzzzz_yy[k] = -g_y_0_zzzz_yy[k] * cd_x[k] + g_y_0_zzzz_xyy[k];

                g_y_0_xzzzz_yz[k] = -g_y_0_zzzz_yz[k] * cd_x[k] + g_y_0_zzzz_xyz[k];

                g_y_0_xzzzz_zz[k] = -g_y_0_zzzz_zz[k] * cd_x[k] + g_y_0_zzzz_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 90);

            auto g_y_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 91);

            auto g_y_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 92);

            auto g_y_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 93);

            auto g_y_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 94);

            auto g_y_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_xx, g_y_0_yyyy_xxy, g_y_0_yyyy_xy, g_y_0_yyyy_xyy, g_y_0_yyyy_xyz, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yyy, g_y_0_yyyy_yyz, g_y_0_yyyy_yz, g_y_0_yyyy_yzz, g_y_0_yyyy_zz, g_y_0_yyyyy_xx, g_y_0_yyyyy_xy, g_y_0_yyyyy_xz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yz, g_y_0_yyyyy_zz, g_yyyy_xx, g_yyyy_xy, g_yyyy_xz, g_yyyy_yy, g_yyyy_yz, g_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xx[k] = -g_yyyy_xx[k] - g_y_0_yyyy_xx[k] * cd_y[k] + g_y_0_yyyy_xxy[k];

                g_y_0_yyyyy_xy[k] = -g_yyyy_xy[k] - g_y_0_yyyy_xy[k] * cd_y[k] + g_y_0_yyyy_xyy[k];

                g_y_0_yyyyy_xz[k] = -g_yyyy_xz[k] - g_y_0_yyyy_xz[k] * cd_y[k] + g_y_0_yyyy_xyz[k];

                g_y_0_yyyyy_yy[k] = -g_yyyy_yy[k] - g_y_0_yyyy_yy[k] * cd_y[k] + g_y_0_yyyy_yyy[k];

                g_y_0_yyyyy_yz[k] = -g_yyyy_yz[k] - g_y_0_yyyy_yz[k] * cd_y[k] + g_y_0_yyyy_yyz[k];

                g_y_0_yyyyy_zz[k] = -g_yyyy_zz[k] - g_y_0_yyyy_zz[k] * cd_y[k] + g_y_0_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 96);

            auto g_y_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 97);

            auto g_y_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 98);

            auto g_y_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 99);

            auto g_y_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 100);

            auto g_y_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_xx, g_y_0_yyyy_xxz, g_y_0_yyyy_xy, g_y_0_yyyy_xyz, g_y_0_yyyy_xz, g_y_0_yyyy_xzz, g_y_0_yyyy_yy, g_y_0_yyyy_yyz, g_y_0_yyyy_yz, g_y_0_yyyy_yzz, g_y_0_yyyy_zz, g_y_0_yyyy_zzz, g_y_0_yyyyz_xx, g_y_0_yyyyz_xy, g_y_0_yyyyz_xz, g_y_0_yyyyz_yy, g_y_0_yyyyz_yz, g_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xx[k] = -g_y_0_yyyy_xx[k] * cd_z[k] + g_y_0_yyyy_xxz[k];

                g_y_0_yyyyz_xy[k] = -g_y_0_yyyy_xy[k] * cd_z[k] + g_y_0_yyyy_xyz[k];

                g_y_0_yyyyz_xz[k] = -g_y_0_yyyy_xz[k] * cd_z[k] + g_y_0_yyyy_xzz[k];

                g_y_0_yyyyz_yy[k] = -g_y_0_yyyy_yy[k] * cd_z[k] + g_y_0_yyyy_yyz[k];

                g_y_0_yyyyz_yz[k] = -g_y_0_yyyy_yz[k] * cd_z[k] + g_y_0_yyyy_yzz[k];

                g_y_0_yyyyz_zz[k] = -g_y_0_yyyy_zz[k] * cd_z[k] + g_y_0_yyyy_zzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 102);

            auto g_y_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 103);

            auto g_y_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 104);

            auto g_y_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 105);

            auto g_y_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 106);

            auto g_y_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_xx, g_y_0_yyyz_xxz, g_y_0_yyyz_xy, g_y_0_yyyz_xyz, g_y_0_yyyz_xz, g_y_0_yyyz_xzz, g_y_0_yyyz_yy, g_y_0_yyyz_yyz, g_y_0_yyyz_yz, g_y_0_yyyz_yzz, g_y_0_yyyz_zz, g_y_0_yyyz_zzz, g_y_0_yyyzz_xx, g_y_0_yyyzz_xy, g_y_0_yyyzz_xz, g_y_0_yyyzz_yy, g_y_0_yyyzz_yz, g_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xx[k] = -g_y_0_yyyz_xx[k] * cd_z[k] + g_y_0_yyyz_xxz[k];

                g_y_0_yyyzz_xy[k] = -g_y_0_yyyz_xy[k] * cd_z[k] + g_y_0_yyyz_xyz[k];

                g_y_0_yyyzz_xz[k] = -g_y_0_yyyz_xz[k] * cd_z[k] + g_y_0_yyyz_xzz[k];

                g_y_0_yyyzz_yy[k] = -g_y_0_yyyz_yy[k] * cd_z[k] + g_y_0_yyyz_yyz[k];

                g_y_0_yyyzz_yz[k] = -g_y_0_yyyz_yz[k] * cd_z[k] + g_y_0_yyyz_yzz[k];

                g_y_0_yyyzz_zz[k] = -g_y_0_yyyz_zz[k] * cd_z[k] + g_y_0_yyyz_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 108);

            auto g_y_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 109);

            auto g_y_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 110);

            auto g_y_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 111);

            auto g_y_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 112);

            auto g_y_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_xx, g_y_0_yyzz_xxz, g_y_0_yyzz_xy, g_y_0_yyzz_xyz, g_y_0_yyzz_xz, g_y_0_yyzz_xzz, g_y_0_yyzz_yy, g_y_0_yyzz_yyz, g_y_0_yyzz_yz, g_y_0_yyzz_yzz, g_y_0_yyzz_zz, g_y_0_yyzz_zzz, g_y_0_yyzzz_xx, g_y_0_yyzzz_xy, g_y_0_yyzzz_xz, g_y_0_yyzzz_yy, g_y_0_yyzzz_yz, g_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xx[k] = -g_y_0_yyzz_xx[k] * cd_z[k] + g_y_0_yyzz_xxz[k];

                g_y_0_yyzzz_xy[k] = -g_y_0_yyzz_xy[k] * cd_z[k] + g_y_0_yyzz_xyz[k];

                g_y_0_yyzzz_xz[k] = -g_y_0_yyzz_xz[k] * cd_z[k] + g_y_0_yyzz_xzz[k];

                g_y_0_yyzzz_yy[k] = -g_y_0_yyzz_yy[k] * cd_z[k] + g_y_0_yyzz_yyz[k];

                g_y_0_yyzzz_yz[k] = -g_y_0_yyzz_yz[k] * cd_z[k] + g_y_0_yyzz_yzz[k];

                g_y_0_yyzzz_zz[k] = -g_y_0_yyzz_zz[k] * cd_z[k] + g_y_0_yyzz_zzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 114);

            auto g_y_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 115);

            auto g_y_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 116);

            auto g_y_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 117);

            auto g_y_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 118);

            auto g_y_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_xx, g_y_0_yzzz_xxz, g_y_0_yzzz_xy, g_y_0_yzzz_xyz, g_y_0_yzzz_xz, g_y_0_yzzz_xzz, g_y_0_yzzz_yy, g_y_0_yzzz_yyz, g_y_0_yzzz_yz, g_y_0_yzzz_yzz, g_y_0_yzzz_zz, g_y_0_yzzz_zzz, g_y_0_yzzzz_xx, g_y_0_yzzzz_xy, g_y_0_yzzzz_xz, g_y_0_yzzzz_yy, g_y_0_yzzzz_yz, g_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xx[k] = -g_y_0_yzzz_xx[k] * cd_z[k] + g_y_0_yzzz_xxz[k];

                g_y_0_yzzzz_xy[k] = -g_y_0_yzzz_xy[k] * cd_z[k] + g_y_0_yzzz_xyz[k];

                g_y_0_yzzzz_xz[k] = -g_y_0_yzzz_xz[k] * cd_z[k] + g_y_0_yzzz_xzz[k];

                g_y_0_yzzzz_yy[k] = -g_y_0_yzzz_yy[k] * cd_z[k] + g_y_0_yzzz_yyz[k];

                g_y_0_yzzzz_yz[k] = -g_y_0_yzzz_yz[k] * cd_z[k] + g_y_0_yzzz_yzz[k];

                g_y_0_yzzzz_zz[k] = -g_y_0_yzzz_zz[k] * cd_z[k] + g_y_0_yzzz_zzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 120);

            auto g_y_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 121);

            auto g_y_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 122);

            auto g_y_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 123);

            auto g_y_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 124);

            auto g_y_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_xx, g_y_0_zzzz_xxz, g_y_0_zzzz_xy, g_y_0_zzzz_xyz, g_y_0_zzzz_xz, g_y_0_zzzz_xzz, g_y_0_zzzz_yy, g_y_0_zzzz_yyz, g_y_0_zzzz_yz, g_y_0_zzzz_yzz, g_y_0_zzzz_zz, g_y_0_zzzz_zzz, g_y_0_zzzzz_xx, g_y_0_zzzzz_xy, g_y_0_zzzzz_xz, g_y_0_zzzzz_yy, g_y_0_zzzzz_yz, g_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xx[k] = -g_y_0_zzzz_xx[k] * cd_z[k] + g_y_0_zzzz_xxz[k];

                g_y_0_zzzzz_xy[k] = -g_y_0_zzzz_xy[k] * cd_z[k] + g_y_0_zzzz_xyz[k];

                g_y_0_zzzzz_xz[k] = -g_y_0_zzzz_xz[k] * cd_z[k] + g_y_0_zzzz_xzz[k];

                g_y_0_zzzzz_yy[k] = -g_y_0_zzzz_yy[k] * cd_z[k] + g_y_0_zzzz_yyz[k];

                g_y_0_zzzzz_yz[k] = -g_y_0_zzzz_yz[k] * cd_z[k] + g_y_0_zzzz_yzz[k];

                g_y_0_zzzzz_zz[k] = -g_y_0_zzzz_zz[k] * cd_z[k] + g_y_0_zzzz_zzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_xx, g_z_0_xxxx_xxx, g_z_0_xxxx_xxy, g_z_0_xxxx_xxz, g_z_0_xxxx_xy, g_z_0_xxxx_xyy, g_z_0_xxxx_xyz, g_z_0_xxxx_xz, g_z_0_xxxx_xzz, g_z_0_xxxx_yy, g_z_0_xxxx_yz, g_z_0_xxxx_zz, g_z_0_xxxxx_xx, g_z_0_xxxxx_xy, g_z_0_xxxxx_xz, g_z_0_xxxxx_yy, g_z_0_xxxxx_yz, g_z_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xx[k] = -g_z_0_xxxx_xx[k] * cd_x[k] + g_z_0_xxxx_xxx[k];

                g_z_0_xxxxx_xy[k] = -g_z_0_xxxx_xy[k] * cd_x[k] + g_z_0_xxxx_xxy[k];

                g_z_0_xxxxx_xz[k] = -g_z_0_xxxx_xz[k] * cd_x[k] + g_z_0_xxxx_xxz[k];

                g_z_0_xxxxx_yy[k] = -g_z_0_xxxx_yy[k] * cd_x[k] + g_z_0_xxxx_xyy[k];

                g_z_0_xxxxx_yz[k] = -g_z_0_xxxx_yz[k] * cd_x[k] + g_z_0_xxxx_xyz[k];

                g_z_0_xxxxx_zz[k] = -g_z_0_xxxx_zz[k] * cd_x[k] + g_z_0_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 6);

            auto g_z_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 7);

            auto g_z_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 8);

            auto g_z_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 9);

            auto g_z_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 10);

            auto g_z_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_xx, g_z_0_xxxxy_xy, g_z_0_xxxxy_xz, g_z_0_xxxxy_yy, g_z_0_xxxxy_yz, g_z_0_xxxxy_zz, g_z_0_xxxy_xx, g_z_0_xxxy_xxx, g_z_0_xxxy_xxy, g_z_0_xxxy_xxz, g_z_0_xxxy_xy, g_z_0_xxxy_xyy, g_z_0_xxxy_xyz, g_z_0_xxxy_xz, g_z_0_xxxy_xzz, g_z_0_xxxy_yy, g_z_0_xxxy_yz, g_z_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xx[k] = -g_z_0_xxxy_xx[k] * cd_x[k] + g_z_0_xxxy_xxx[k];

                g_z_0_xxxxy_xy[k] = -g_z_0_xxxy_xy[k] * cd_x[k] + g_z_0_xxxy_xxy[k];

                g_z_0_xxxxy_xz[k] = -g_z_0_xxxy_xz[k] * cd_x[k] + g_z_0_xxxy_xxz[k];

                g_z_0_xxxxy_yy[k] = -g_z_0_xxxy_yy[k] * cd_x[k] + g_z_0_xxxy_xyy[k];

                g_z_0_xxxxy_yz[k] = -g_z_0_xxxy_yz[k] * cd_x[k] + g_z_0_xxxy_xyz[k];

                g_z_0_xxxxy_zz[k] = -g_z_0_xxxy_zz[k] * cd_x[k] + g_z_0_xxxy_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 12);

            auto g_z_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 13);

            auto g_z_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 14);

            auto g_z_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 15);

            auto g_z_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 16);

            auto g_z_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_xx, g_z_0_xxxxz_xy, g_z_0_xxxxz_xz, g_z_0_xxxxz_yy, g_z_0_xxxxz_yz, g_z_0_xxxxz_zz, g_z_0_xxxz_xx, g_z_0_xxxz_xxx, g_z_0_xxxz_xxy, g_z_0_xxxz_xxz, g_z_0_xxxz_xy, g_z_0_xxxz_xyy, g_z_0_xxxz_xyz, g_z_0_xxxz_xz, g_z_0_xxxz_xzz, g_z_0_xxxz_yy, g_z_0_xxxz_yz, g_z_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xx[k] = -g_z_0_xxxz_xx[k] * cd_x[k] + g_z_0_xxxz_xxx[k];

                g_z_0_xxxxz_xy[k] = -g_z_0_xxxz_xy[k] * cd_x[k] + g_z_0_xxxz_xxy[k];

                g_z_0_xxxxz_xz[k] = -g_z_0_xxxz_xz[k] * cd_x[k] + g_z_0_xxxz_xxz[k];

                g_z_0_xxxxz_yy[k] = -g_z_0_xxxz_yy[k] * cd_x[k] + g_z_0_xxxz_xyy[k];

                g_z_0_xxxxz_yz[k] = -g_z_0_xxxz_yz[k] * cd_x[k] + g_z_0_xxxz_xyz[k];

                g_z_0_xxxxz_zz[k] = -g_z_0_xxxz_zz[k] * cd_x[k] + g_z_0_xxxz_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 18);

            auto g_z_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 19);

            auto g_z_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 20);

            auto g_z_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 21);

            auto g_z_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 22);

            auto g_z_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_xx, g_z_0_xxxyy_xy, g_z_0_xxxyy_xz, g_z_0_xxxyy_yy, g_z_0_xxxyy_yz, g_z_0_xxxyy_zz, g_z_0_xxyy_xx, g_z_0_xxyy_xxx, g_z_0_xxyy_xxy, g_z_0_xxyy_xxz, g_z_0_xxyy_xy, g_z_0_xxyy_xyy, g_z_0_xxyy_xyz, g_z_0_xxyy_xz, g_z_0_xxyy_xzz, g_z_0_xxyy_yy, g_z_0_xxyy_yz, g_z_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xx[k] = -g_z_0_xxyy_xx[k] * cd_x[k] + g_z_0_xxyy_xxx[k];

                g_z_0_xxxyy_xy[k] = -g_z_0_xxyy_xy[k] * cd_x[k] + g_z_0_xxyy_xxy[k];

                g_z_0_xxxyy_xz[k] = -g_z_0_xxyy_xz[k] * cd_x[k] + g_z_0_xxyy_xxz[k];

                g_z_0_xxxyy_yy[k] = -g_z_0_xxyy_yy[k] * cd_x[k] + g_z_0_xxyy_xyy[k];

                g_z_0_xxxyy_yz[k] = -g_z_0_xxyy_yz[k] * cd_x[k] + g_z_0_xxyy_xyz[k];

                g_z_0_xxxyy_zz[k] = -g_z_0_xxyy_zz[k] * cd_x[k] + g_z_0_xxyy_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 24);

            auto g_z_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 25);

            auto g_z_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 26);

            auto g_z_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 27);

            auto g_z_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 28);

            auto g_z_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_xx, g_z_0_xxxyz_xy, g_z_0_xxxyz_xz, g_z_0_xxxyz_yy, g_z_0_xxxyz_yz, g_z_0_xxxyz_zz, g_z_0_xxyz_xx, g_z_0_xxyz_xxx, g_z_0_xxyz_xxy, g_z_0_xxyz_xxz, g_z_0_xxyz_xy, g_z_0_xxyz_xyy, g_z_0_xxyz_xyz, g_z_0_xxyz_xz, g_z_0_xxyz_xzz, g_z_0_xxyz_yy, g_z_0_xxyz_yz, g_z_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xx[k] = -g_z_0_xxyz_xx[k] * cd_x[k] + g_z_0_xxyz_xxx[k];

                g_z_0_xxxyz_xy[k] = -g_z_0_xxyz_xy[k] * cd_x[k] + g_z_0_xxyz_xxy[k];

                g_z_0_xxxyz_xz[k] = -g_z_0_xxyz_xz[k] * cd_x[k] + g_z_0_xxyz_xxz[k];

                g_z_0_xxxyz_yy[k] = -g_z_0_xxyz_yy[k] * cd_x[k] + g_z_0_xxyz_xyy[k];

                g_z_0_xxxyz_yz[k] = -g_z_0_xxyz_yz[k] * cd_x[k] + g_z_0_xxyz_xyz[k];

                g_z_0_xxxyz_zz[k] = -g_z_0_xxyz_zz[k] * cd_x[k] + g_z_0_xxyz_xzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 30);

            auto g_z_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 31);

            auto g_z_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 32);

            auto g_z_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 33);

            auto g_z_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 34);

            auto g_z_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_xx, g_z_0_xxxzz_xy, g_z_0_xxxzz_xz, g_z_0_xxxzz_yy, g_z_0_xxxzz_yz, g_z_0_xxxzz_zz, g_z_0_xxzz_xx, g_z_0_xxzz_xxx, g_z_0_xxzz_xxy, g_z_0_xxzz_xxz, g_z_0_xxzz_xy, g_z_0_xxzz_xyy, g_z_0_xxzz_xyz, g_z_0_xxzz_xz, g_z_0_xxzz_xzz, g_z_0_xxzz_yy, g_z_0_xxzz_yz, g_z_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xx[k] = -g_z_0_xxzz_xx[k] * cd_x[k] + g_z_0_xxzz_xxx[k];

                g_z_0_xxxzz_xy[k] = -g_z_0_xxzz_xy[k] * cd_x[k] + g_z_0_xxzz_xxy[k];

                g_z_0_xxxzz_xz[k] = -g_z_0_xxzz_xz[k] * cd_x[k] + g_z_0_xxzz_xxz[k];

                g_z_0_xxxzz_yy[k] = -g_z_0_xxzz_yy[k] * cd_x[k] + g_z_0_xxzz_xyy[k];

                g_z_0_xxxzz_yz[k] = -g_z_0_xxzz_yz[k] * cd_x[k] + g_z_0_xxzz_xyz[k];

                g_z_0_xxxzz_zz[k] = -g_z_0_xxzz_zz[k] * cd_x[k] + g_z_0_xxzz_xzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 36);

            auto g_z_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 37);

            auto g_z_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 38);

            auto g_z_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 39);

            auto g_z_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 40);

            auto g_z_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_xx, g_z_0_xxyyy_xy, g_z_0_xxyyy_xz, g_z_0_xxyyy_yy, g_z_0_xxyyy_yz, g_z_0_xxyyy_zz, g_z_0_xyyy_xx, g_z_0_xyyy_xxx, g_z_0_xyyy_xxy, g_z_0_xyyy_xxz, g_z_0_xyyy_xy, g_z_0_xyyy_xyy, g_z_0_xyyy_xyz, g_z_0_xyyy_xz, g_z_0_xyyy_xzz, g_z_0_xyyy_yy, g_z_0_xyyy_yz, g_z_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xx[k] = -g_z_0_xyyy_xx[k] * cd_x[k] + g_z_0_xyyy_xxx[k];

                g_z_0_xxyyy_xy[k] = -g_z_0_xyyy_xy[k] * cd_x[k] + g_z_0_xyyy_xxy[k];

                g_z_0_xxyyy_xz[k] = -g_z_0_xyyy_xz[k] * cd_x[k] + g_z_0_xyyy_xxz[k];

                g_z_0_xxyyy_yy[k] = -g_z_0_xyyy_yy[k] * cd_x[k] + g_z_0_xyyy_xyy[k];

                g_z_0_xxyyy_yz[k] = -g_z_0_xyyy_yz[k] * cd_x[k] + g_z_0_xyyy_xyz[k];

                g_z_0_xxyyy_zz[k] = -g_z_0_xyyy_zz[k] * cd_x[k] + g_z_0_xyyy_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 42);

            auto g_z_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 43);

            auto g_z_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 44);

            auto g_z_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 45);

            auto g_z_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 46);

            auto g_z_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_xx, g_z_0_xxyyz_xy, g_z_0_xxyyz_xz, g_z_0_xxyyz_yy, g_z_0_xxyyz_yz, g_z_0_xxyyz_zz, g_z_0_xyyz_xx, g_z_0_xyyz_xxx, g_z_0_xyyz_xxy, g_z_0_xyyz_xxz, g_z_0_xyyz_xy, g_z_0_xyyz_xyy, g_z_0_xyyz_xyz, g_z_0_xyyz_xz, g_z_0_xyyz_xzz, g_z_0_xyyz_yy, g_z_0_xyyz_yz, g_z_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xx[k] = -g_z_0_xyyz_xx[k] * cd_x[k] + g_z_0_xyyz_xxx[k];

                g_z_0_xxyyz_xy[k] = -g_z_0_xyyz_xy[k] * cd_x[k] + g_z_0_xyyz_xxy[k];

                g_z_0_xxyyz_xz[k] = -g_z_0_xyyz_xz[k] * cd_x[k] + g_z_0_xyyz_xxz[k];

                g_z_0_xxyyz_yy[k] = -g_z_0_xyyz_yy[k] * cd_x[k] + g_z_0_xyyz_xyy[k];

                g_z_0_xxyyz_yz[k] = -g_z_0_xyyz_yz[k] * cd_x[k] + g_z_0_xyyz_xyz[k];

                g_z_0_xxyyz_zz[k] = -g_z_0_xyyz_zz[k] * cd_x[k] + g_z_0_xyyz_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 48);

            auto g_z_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 49);

            auto g_z_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 50);

            auto g_z_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 51);

            auto g_z_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 52);

            auto g_z_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_xx, g_z_0_xxyzz_xy, g_z_0_xxyzz_xz, g_z_0_xxyzz_yy, g_z_0_xxyzz_yz, g_z_0_xxyzz_zz, g_z_0_xyzz_xx, g_z_0_xyzz_xxx, g_z_0_xyzz_xxy, g_z_0_xyzz_xxz, g_z_0_xyzz_xy, g_z_0_xyzz_xyy, g_z_0_xyzz_xyz, g_z_0_xyzz_xz, g_z_0_xyzz_xzz, g_z_0_xyzz_yy, g_z_0_xyzz_yz, g_z_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xx[k] = -g_z_0_xyzz_xx[k] * cd_x[k] + g_z_0_xyzz_xxx[k];

                g_z_0_xxyzz_xy[k] = -g_z_0_xyzz_xy[k] * cd_x[k] + g_z_0_xyzz_xxy[k];

                g_z_0_xxyzz_xz[k] = -g_z_0_xyzz_xz[k] * cd_x[k] + g_z_0_xyzz_xxz[k];

                g_z_0_xxyzz_yy[k] = -g_z_0_xyzz_yy[k] * cd_x[k] + g_z_0_xyzz_xyy[k];

                g_z_0_xxyzz_yz[k] = -g_z_0_xyzz_yz[k] * cd_x[k] + g_z_0_xyzz_xyz[k];

                g_z_0_xxyzz_zz[k] = -g_z_0_xyzz_zz[k] * cd_x[k] + g_z_0_xyzz_xzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 54);

            auto g_z_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 55);

            auto g_z_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 56);

            auto g_z_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 57);

            auto g_z_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 58);

            auto g_z_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_xx, g_z_0_xxzzz_xy, g_z_0_xxzzz_xz, g_z_0_xxzzz_yy, g_z_0_xxzzz_yz, g_z_0_xxzzz_zz, g_z_0_xzzz_xx, g_z_0_xzzz_xxx, g_z_0_xzzz_xxy, g_z_0_xzzz_xxz, g_z_0_xzzz_xy, g_z_0_xzzz_xyy, g_z_0_xzzz_xyz, g_z_0_xzzz_xz, g_z_0_xzzz_xzz, g_z_0_xzzz_yy, g_z_0_xzzz_yz, g_z_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xx[k] = -g_z_0_xzzz_xx[k] * cd_x[k] + g_z_0_xzzz_xxx[k];

                g_z_0_xxzzz_xy[k] = -g_z_0_xzzz_xy[k] * cd_x[k] + g_z_0_xzzz_xxy[k];

                g_z_0_xxzzz_xz[k] = -g_z_0_xzzz_xz[k] * cd_x[k] + g_z_0_xzzz_xxz[k];

                g_z_0_xxzzz_yy[k] = -g_z_0_xzzz_yy[k] * cd_x[k] + g_z_0_xzzz_xyy[k];

                g_z_0_xxzzz_yz[k] = -g_z_0_xzzz_yz[k] * cd_x[k] + g_z_0_xzzz_xyz[k];

                g_z_0_xxzzz_zz[k] = -g_z_0_xzzz_zz[k] * cd_x[k] + g_z_0_xzzz_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 60);

            auto g_z_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 61);

            auto g_z_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 62);

            auto g_z_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 63);

            auto g_z_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 64);

            auto g_z_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_xx, g_z_0_xyyyy_xy, g_z_0_xyyyy_xz, g_z_0_xyyyy_yy, g_z_0_xyyyy_yz, g_z_0_xyyyy_zz, g_z_0_yyyy_xx, g_z_0_yyyy_xxx, g_z_0_yyyy_xxy, g_z_0_yyyy_xxz, g_z_0_yyyy_xy, g_z_0_yyyy_xyy, g_z_0_yyyy_xyz, g_z_0_yyyy_xz, g_z_0_yyyy_xzz, g_z_0_yyyy_yy, g_z_0_yyyy_yz, g_z_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xx[k] = -g_z_0_yyyy_xx[k] * cd_x[k] + g_z_0_yyyy_xxx[k];

                g_z_0_xyyyy_xy[k] = -g_z_0_yyyy_xy[k] * cd_x[k] + g_z_0_yyyy_xxy[k];

                g_z_0_xyyyy_xz[k] = -g_z_0_yyyy_xz[k] * cd_x[k] + g_z_0_yyyy_xxz[k];

                g_z_0_xyyyy_yy[k] = -g_z_0_yyyy_yy[k] * cd_x[k] + g_z_0_yyyy_xyy[k];

                g_z_0_xyyyy_yz[k] = -g_z_0_yyyy_yz[k] * cd_x[k] + g_z_0_yyyy_xyz[k];

                g_z_0_xyyyy_zz[k] = -g_z_0_yyyy_zz[k] * cd_x[k] + g_z_0_yyyy_xzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 66);

            auto g_z_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 67);

            auto g_z_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 68);

            auto g_z_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 69);

            auto g_z_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 70);

            auto g_z_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_xx, g_z_0_xyyyz_xy, g_z_0_xyyyz_xz, g_z_0_xyyyz_yy, g_z_0_xyyyz_yz, g_z_0_xyyyz_zz, g_z_0_yyyz_xx, g_z_0_yyyz_xxx, g_z_0_yyyz_xxy, g_z_0_yyyz_xxz, g_z_0_yyyz_xy, g_z_0_yyyz_xyy, g_z_0_yyyz_xyz, g_z_0_yyyz_xz, g_z_0_yyyz_xzz, g_z_0_yyyz_yy, g_z_0_yyyz_yz, g_z_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xx[k] = -g_z_0_yyyz_xx[k] * cd_x[k] + g_z_0_yyyz_xxx[k];

                g_z_0_xyyyz_xy[k] = -g_z_0_yyyz_xy[k] * cd_x[k] + g_z_0_yyyz_xxy[k];

                g_z_0_xyyyz_xz[k] = -g_z_0_yyyz_xz[k] * cd_x[k] + g_z_0_yyyz_xxz[k];

                g_z_0_xyyyz_yy[k] = -g_z_0_yyyz_yy[k] * cd_x[k] + g_z_0_yyyz_xyy[k];

                g_z_0_xyyyz_yz[k] = -g_z_0_yyyz_yz[k] * cd_x[k] + g_z_0_yyyz_xyz[k];

                g_z_0_xyyyz_zz[k] = -g_z_0_yyyz_zz[k] * cd_x[k] + g_z_0_yyyz_xzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 72);

            auto g_z_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 73);

            auto g_z_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 74);

            auto g_z_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 75);

            auto g_z_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 76);

            auto g_z_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_xx, g_z_0_xyyzz_xy, g_z_0_xyyzz_xz, g_z_0_xyyzz_yy, g_z_0_xyyzz_yz, g_z_0_xyyzz_zz, g_z_0_yyzz_xx, g_z_0_yyzz_xxx, g_z_0_yyzz_xxy, g_z_0_yyzz_xxz, g_z_0_yyzz_xy, g_z_0_yyzz_xyy, g_z_0_yyzz_xyz, g_z_0_yyzz_xz, g_z_0_yyzz_xzz, g_z_0_yyzz_yy, g_z_0_yyzz_yz, g_z_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xx[k] = -g_z_0_yyzz_xx[k] * cd_x[k] + g_z_0_yyzz_xxx[k];

                g_z_0_xyyzz_xy[k] = -g_z_0_yyzz_xy[k] * cd_x[k] + g_z_0_yyzz_xxy[k];

                g_z_0_xyyzz_xz[k] = -g_z_0_yyzz_xz[k] * cd_x[k] + g_z_0_yyzz_xxz[k];

                g_z_0_xyyzz_yy[k] = -g_z_0_yyzz_yy[k] * cd_x[k] + g_z_0_yyzz_xyy[k];

                g_z_0_xyyzz_yz[k] = -g_z_0_yyzz_yz[k] * cd_x[k] + g_z_0_yyzz_xyz[k];

                g_z_0_xyyzz_zz[k] = -g_z_0_yyzz_zz[k] * cd_x[k] + g_z_0_yyzz_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 78);

            auto g_z_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 79);

            auto g_z_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 80);

            auto g_z_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 81);

            auto g_z_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 82);

            auto g_z_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_xx, g_z_0_xyzzz_xy, g_z_0_xyzzz_xz, g_z_0_xyzzz_yy, g_z_0_xyzzz_yz, g_z_0_xyzzz_zz, g_z_0_yzzz_xx, g_z_0_yzzz_xxx, g_z_0_yzzz_xxy, g_z_0_yzzz_xxz, g_z_0_yzzz_xy, g_z_0_yzzz_xyy, g_z_0_yzzz_xyz, g_z_0_yzzz_xz, g_z_0_yzzz_xzz, g_z_0_yzzz_yy, g_z_0_yzzz_yz, g_z_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xx[k] = -g_z_0_yzzz_xx[k] * cd_x[k] + g_z_0_yzzz_xxx[k];

                g_z_0_xyzzz_xy[k] = -g_z_0_yzzz_xy[k] * cd_x[k] + g_z_0_yzzz_xxy[k];

                g_z_0_xyzzz_xz[k] = -g_z_0_yzzz_xz[k] * cd_x[k] + g_z_0_yzzz_xxz[k];

                g_z_0_xyzzz_yy[k] = -g_z_0_yzzz_yy[k] * cd_x[k] + g_z_0_yzzz_xyy[k];

                g_z_0_xyzzz_yz[k] = -g_z_0_yzzz_yz[k] * cd_x[k] + g_z_0_yzzz_xyz[k];

                g_z_0_xyzzz_zz[k] = -g_z_0_yzzz_zz[k] * cd_x[k] + g_z_0_yzzz_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 84);

            auto g_z_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 85);

            auto g_z_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 86);

            auto g_z_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 87);

            auto g_z_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 88);

            auto g_z_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_xx, g_z_0_xzzzz_xy, g_z_0_xzzzz_xz, g_z_0_xzzzz_yy, g_z_0_xzzzz_yz, g_z_0_xzzzz_zz, g_z_0_zzzz_xx, g_z_0_zzzz_xxx, g_z_0_zzzz_xxy, g_z_0_zzzz_xxz, g_z_0_zzzz_xy, g_z_0_zzzz_xyy, g_z_0_zzzz_xyz, g_z_0_zzzz_xz, g_z_0_zzzz_xzz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xx[k] = -g_z_0_zzzz_xx[k] * cd_x[k] + g_z_0_zzzz_xxx[k];

                g_z_0_xzzzz_xy[k] = -g_z_0_zzzz_xy[k] * cd_x[k] + g_z_0_zzzz_xxy[k];

                g_z_0_xzzzz_xz[k] = -g_z_0_zzzz_xz[k] * cd_x[k] + g_z_0_zzzz_xxz[k];

                g_z_0_xzzzz_yy[k] = -g_z_0_zzzz_yy[k] * cd_x[k] + g_z_0_zzzz_xyy[k];

                g_z_0_xzzzz_yz[k] = -g_z_0_zzzz_yz[k] * cd_x[k] + g_z_0_zzzz_xyz[k];

                g_z_0_xzzzz_zz[k] = -g_z_0_zzzz_zz[k] * cd_x[k] + g_z_0_zzzz_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 90);

            auto g_z_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 91);

            auto g_z_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 92);

            auto g_z_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 93);

            auto g_z_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 94);

            auto g_z_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_xx, g_z_0_yyyy_xxy, g_z_0_yyyy_xy, g_z_0_yyyy_xyy, g_z_0_yyyy_xyz, g_z_0_yyyy_xz, g_z_0_yyyy_yy, g_z_0_yyyy_yyy, g_z_0_yyyy_yyz, g_z_0_yyyy_yz, g_z_0_yyyy_yzz, g_z_0_yyyy_zz, g_z_0_yyyyy_xx, g_z_0_yyyyy_xy, g_z_0_yyyyy_xz, g_z_0_yyyyy_yy, g_z_0_yyyyy_yz, g_z_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xx[k] = -g_z_0_yyyy_xx[k] * cd_y[k] + g_z_0_yyyy_xxy[k];

                g_z_0_yyyyy_xy[k] = -g_z_0_yyyy_xy[k] * cd_y[k] + g_z_0_yyyy_xyy[k];

                g_z_0_yyyyy_xz[k] = -g_z_0_yyyy_xz[k] * cd_y[k] + g_z_0_yyyy_xyz[k];

                g_z_0_yyyyy_yy[k] = -g_z_0_yyyy_yy[k] * cd_y[k] + g_z_0_yyyy_yyy[k];

                g_z_0_yyyyy_yz[k] = -g_z_0_yyyy_yz[k] * cd_y[k] + g_z_0_yyyy_yyz[k];

                g_z_0_yyyyy_zz[k] = -g_z_0_yyyy_zz[k] * cd_y[k] + g_z_0_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 96);

            auto g_z_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 97);

            auto g_z_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 98);

            auto g_z_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 99);

            auto g_z_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 100);

            auto g_z_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_xx, g_z_0_yyyyz_xy, g_z_0_yyyyz_xz, g_z_0_yyyyz_yy, g_z_0_yyyyz_yz, g_z_0_yyyyz_zz, g_z_0_yyyz_xx, g_z_0_yyyz_xxy, g_z_0_yyyz_xy, g_z_0_yyyz_xyy, g_z_0_yyyz_xyz, g_z_0_yyyz_xz, g_z_0_yyyz_yy, g_z_0_yyyz_yyy, g_z_0_yyyz_yyz, g_z_0_yyyz_yz, g_z_0_yyyz_yzz, g_z_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xx[k] = -g_z_0_yyyz_xx[k] * cd_y[k] + g_z_0_yyyz_xxy[k];

                g_z_0_yyyyz_xy[k] = -g_z_0_yyyz_xy[k] * cd_y[k] + g_z_0_yyyz_xyy[k];

                g_z_0_yyyyz_xz[k] = -g_z_0_yyyz_xz[k] * cd_y[k] + g_z_0_yyyz_xyz[k];

                g_z_0_yyyyz_yy[k] = -g_z_0_yyyz_yy[k] * cd_y[k] + g_z_0_yyyz_yyy[k];

                g_z_0_yyyyz_yz[k] = -g_z_0_yyyz_yz[k] * cd_y[k] + g_z_0_yyyz_yyz[k];

                g_z_0_yyyyz_zz[k] = -g_z_0_yyyz_zz[k] * cd_y[k] + g_z_0_yyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 102);

            auto g_z_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 103);

            auto g_z_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 104);

            auto g_z_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 105);

            auto g_z_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 106);

            auto g_z_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_xx, g_z_0_yyyzz_xy, g_z_0_yyyzz_xz, g_z_0_yyyzz_yy, g_z_0_yyyzz_yz, g_z_0_yyyzz_zz, g_z_0_yyzz_xx, g_z_0_yyzz_xxy, g_z_0_yyzz_xy, g_z_0_yyzz_xyy, g_z_0_yyzz_xyz, g_z_0_yyzz_xz, g_z_0_yyzz_yy, g_z_0_yyzz_yyy, g_z_0_yyzz_yyz, g_z_0_yyzz_yz, g_z_0_yyzz_yzz, g_z_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xx[k] = -g_z_0_yyzz_xx[k] * cd_y[k] + g_z_0_yyzz_xxy[k];

                g_z_0_yyyzz_xy[k] = -g_z_0_yyzz_xy[k] * cd_y[k] + g_z_0_yyzz_xyy[k];

                g_z_0_yyyzz_xz[k] = -g_z_0_yyzz_xz[k] * cd_y[k] + g_z_0_yyzz_xyz[k];

                g_z_0_yyyzz_yy[k] = -g_z_0_yyzz_yy[k] * cd_y[k] + g_z_0_yyzz_yyy[k];

                g_z_0_yyyzz_yz[k] = -g_z_0_yyzz_yz[k] * cd_y[k] + g_z_0_yyzz_yyz[k];

                g_z_0_yyyzz_zz[k] = -g_z_0_yyzz_zz[k] * cd_y[k] + g_z_0_yyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 108);

            auto g_z_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 109);

            auto g_z_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 110);

            auto g_z_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 111);

            auto g_z_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 112);

            auto g_z_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_xx, g_z_0_yyzzz_xy, g_z_0_yyzzz_xz, g_z_0_yyzzz_yy, g_z_0_yyzzz_yz, g_z_0_yyzzz_zz, g_z_0_yzzz_xx, g_z_0_yzzz_xxy, g_z_0_yzzz_xy, g_z_0_yzzz_xyy, g_z_0_yzzz_xyz, g_z_0_yzzz_xz, g_z_0_yzzz_yy, g_z_0_yzzz_yyy, g_z_0_yzzz_yyz, g_z_0_yzzz_yz, g_z_0_yzzz_yzz, g_z_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xx[k] = -g_z_0_yzzz_xx[k] * cd_y[k] + g_z_0_yzzz_xxy[k];

                g_z_0_yyzzz_xy[k] = -g_z_0_yzzz_xy[k] * cd_y[k] + g_z_0_yzzz_xyy[k];

                g_z_0_yyzzz_xz[k] = -g_z_0_yzzz_xz[k] * cd_y[k] + g_z_0_yzzz_xyz[k];

                g_z_0_yyzzz_yy[k] = -g_z_0_yzzz_yy[k] * cd_y[k] + g_z_0_yzzz_yyy[k];

                g_z_0_yyzzz_yz[k] = -g_z_0_yzzz_yz[k] * cd_y[k] + g_z_0_yzzz_yyz[k];

                g_z_0_yyzzz_zz[k] = -g_z_0_yzzz_zz[k] * cd_y[k] + g_z_0_yzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 114);

            auto g_z_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 115);

            auto g_z_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 116);

            auto g_z_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 117);

            auto g_z_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 118);

            auto g_z_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_xx, g_z_0_yzzzz_xy, g_z_0_yzzzz_xz, g_z_0_yzzzz_yy, g_z_0_yzzzz_yz, g_z_0_yzzzz_zz, g_z_0_zzzz_xx, g_z_0_zzzz_xxy, g_z_0_zzzz_xy, g_z_0_zzzz_xyy, g_z_0_zzzz_xyz, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yyy, g_z_0_zzzz_yyz, g_z_0_zzzz_yz, g_z_0_zzzz_yzz, g_z_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xx[k] = -g_z_0_zzzz_xx[k] * cd_y[k] + g_z_0_zzzz_xxy[k];

                g_z_0_yzzzz_xy[k] = -g_z_0_zzzz_xy[k] * cd_y[k] + g_z_0_zzzz_xyy[k];

                g_z_0_yzzzz_xz[k] = -g_z_0_zzzz_xz[k] * cd_y[k] + g_z_0_zzzz_xyz[k];

                g_z_0_yzzzz_yy[k] = -g_z_0_zzzz_yy[k] * cd_y[k] + g_z_0_zzzz_yyy[k];

                g_z_0_yzzzz_yz[k] = -g_z_0_zzzz_yz[k] * cd_y[k] + g_z_0_zzzz_yyz[k];

                g_z_0_yzzzz_zz[k] = -g_z_0_zzzz_zz[k] * cd_y[k] + g_z_0_zzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 120);

            auto g_z_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 121);

            auto g_z_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 122);

            auto g_z_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 123);

            auto g_z_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 124);

            auto g_z_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_xx, g_z_0_zzzz_xxz, g_z_0_zzzz_xy, g_z_0_zzzz_xyz, g_z_0_zzzz_xz, g_z_0_zzzz_xzz, g_z_0_zzzz_yy, g_z_0_zzzz_yyz, g_z_0_zzzz_yz, g_z_0_zzzz_yzz, g_z_0_zzzz_zz, g_z_0_zzzz_zzz, g_z_0_zzzzz_xx, g_z_0_zzzzz_xy, g_z_0_zzzzz_xz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_zz, g_zzzz_xx, g_zzzz_xy, g_zzzz_xz, g_zzzz_yy, g_zzzz_yz, g_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xx[k] = -g_zzzz_xx[k] - g_z_0_zzzz_xx[k] * cd_z[k] + g_z_0_zzzz_xxz[k];

                g_z_0_zzzzz_xy[k] = -g_zzzz_xy[k] - g_z_0_zzzz_xy[k] * cd_z[k] + g_z_0_zzzz_xyz[k];

                g_z_0_zzzzz_xz[k] = -g_zzzz_xz[k] - g_z_0_zzzz_xz[k] * cd_z[k] + g_z_0_zzzz_xzz[k];

                g_z_0_zzzzz_yy[k] = -g_zzzz_yy[k] - g_z_0_zzzz_yy[k] * cd_z[k] + g_z_0_zzzz_yyz[k];

                g_z_0_zzzzz_yz[k] = -g_zzzz_yz[k] - g_z_0_zzzz_yz[k] * cd_z[k] + g_z_0_zzzz_yzz[k];

                g_z_0_zzzzz_zz[k] = -g_zzzz_zz[k] - g_z_0_zzzz_zz[k] * cd_z[k] + g_z_0_zzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

