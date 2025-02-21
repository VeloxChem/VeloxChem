#include "ThreeCenterElectronRepulsionGeom010ContrRecXFF.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xff(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xff,
                                        const size_t idx_xdf,
                                        const size_t idx_geom_10_xdf,
                                        const size_t idx_geom_10_xdg,
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
        /// Set up components of auxilary buffer : SDF

        const auto df_off = idx_xdf + i * 60;

        auto g_xx_xxx = cbuffer.data(df_off + 0);

        auto g_xx_xxy = cbuffer.data(df_off + 1);

        auto g_xx_xxz = cbuffer.data(df_off + 2);

        auto g_xx_xyy = cbuffer.data(df_off + 3);

        auto g_xx_xyz = cbuffer.data(df_off + 4);

        auto g_xx_xzz = cbuffer.data(df_off + 5);

        auto g_xx_yyy = cbuffer.data(df_off + 6);

        auto g_xx_yyz = cbuffer.data(df_off + 7);

        auto g_xx_yzz = cbuffer.data(df_off + 8);

        auto g_xx_zzz = cbuffer.data(df_off + 9);

        auto g_xy_xxx = cbuffer.data(df_off + 10);

        auto g_xy_xxy = cbuffer.data(df_off + 11);

        auto g_xy_xxz = cbuffer.data(df_off + 12);

        auto g_xy_xyy = cbuffer.data(df_off + 13);

        auto g_xy_xyz = cbuffer.data(df_off + 14);

        auto g_xy_xzz = cbuffer.data(df_off + 15);

        auto g_xy_yyy = cbuffer.data(df_off + 16);

        auto g_xy_yyz = cbuffer.data(df_off + 17);

        auto g_xy_yzz = cbuffer.data(df_off + 18);

        auto g_xy_zzz = cbuffer.data(df_off + 19);

        auto g_xz_xxx = cbuffer.data(df_off + 20);

        auto g_xz_xxy = cbuffer.data(df_off + 21);

        auto g_xz_xxz = cbuffer.data(df_off + 22);

        auto g_xz_xyy = cbuffer.data(df_off + 23);

        auto g_xz_xyz = cbuffer.data(df_off + 24);

        auto g_xz_xzz = cbuffer.data(df_off + 25);

        auto g_xz_yyy = cbuffer.data(df_off + 26);

        auto g_xz_yyz = cbuffer.data(df_off + 27);

        auto g_xz_yzz = cbuffer.data(df_off + 28);

        auto g_xz_zzz = cbuffer.data(df_off + 29);

        auto g_yy_xxx = cbuffer.data(df_off + 30);

        auto g_yy_xxy = cbuffer.data(df_off + 31);

        auto g_yy_xxz = cbuffer.data(df_off + 32);

        auto g_yy_xyy = cbuffer.data(df_off + 33);

        auto g_yy_xyz = cbuffer.data(df_off + 34);

        auto g_yy_xzz = cbuffer.data(df_off + 35);

        auto g_yy_yyy = cbuffer.data(df_off + 36);

        auto g_yy_yyz = cbuffer.data(df_off + 37);

        auto g_yy_yzz = cbuffer.data(df_off + 38);

        auto g_yy_zzz = cbuffer.data(df_off + 39);

        auto g_yz_xxx = cbuffer.data(df_off + 40);

        auto g_yz_xxy = cbuffer.data(df_off + 41);

        auto g_yz_xxz = cbuffer.data(df_off + 42);

        auto g_yz_xyy = cbuffer.data(df_off + 43);

        auto g_yz_xyz = cbuffer.data(df_off + 44);

        auto g_yz_xzz = cbuffer.data(df_off + 45);

        auto g_yz_yyy = cbuffer.data(df_off + 46);

        auto g_yz_yyz = cbuffer.data(df_off + 47);

        auto g_yz_yzz = cbuffer.data(df_off + 48);

        auto g_yz_zzz = cbuffer.data(df_off + 49);

        auto g_zz_xxx = cbuffer.data(df_off + 50);

        auto g_zz_xxy = cbuffer.data(df_off + 51);

        auto g_zz_xxz = cbuffer.data(df_off + 52);

        auto g_zz_xyy = cbuffer.data(df_off + 53);

        auto g_zz_xyz = cbuffer.data(df_off + 54);

        auto g_zz_xzz = cbuffer.data(df_off + 55);

        auto g_zz_yyy = cbuffer.data(df_off + 56);

        auto g_zz_yyz = cbuffer.data(df_off + 57);

        auto g_zz_yzz = cbuffer.data(df_off + 58);

        auto g_zz_zzz = cbuffer.data(df_off + 59);

        /// Set up components of auxilary buffer : SDF

        const auto df_geom_10_off = idx_geom_10_xdf + i * 60;

        auto g_x_0_xx_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xy_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xy_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xy_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xy_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 29);

        auto g_x_0_yy_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 30);

        auto g_x_0_yy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 31);

        auto g_x_0_yy_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 32);

        auto g_x_0_yy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 33);

        auto g_x_0_yy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 34);

        auto g_x_0_yy_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 35);

        auto g_x_0_yy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 36);

        auto g_x_0_yy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 37);

        auto g_x_0_yy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 38);

        auto g_x_0_yy_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 39);

        auto g_x_0_yz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 40);

        auto g_x_0_yz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 41);

        auto g_x_0_yz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 42);

        auto g_x_0_yz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 43);

        auto g_x_0_yz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 44);

        auto g_x_0_yz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 45);

        auto g_x_0_yz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 46);

        auto g_x_0_yz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 47);

        auto g_x_0_yz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 48);

        auto g_x_0_yz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 49);

        auto g_x_0_zz_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 50);

        auto g_x_0_zz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 51);

        auto g_x_0_zz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 52);

        auto g_x_0_zz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 53);

        auto g_x_0_zz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 54);

        auto g_x_0_zz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 55);

        auto g_x_0_zz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 56);

        auto g_x_0_zz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 57);

        auto g_x_0_zz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 58);

        auto g_x_0_zz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 59);

        auto g_y_0_xx_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 0);

        auto g_y_0_xx_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 1);

        auto g_y_0_xx_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 2);

        auto g_y_0_xx_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 3);

        auto g_y_0_xx_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 4);

        auto g_y_0_xx_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 5);

        auto g_y_0_xx_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 6);

        auto g_y_0_xx_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 7);

        auto g_y_0_xx_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 8);

        auto g_y_0_xx_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 9);

        auto g_y_0_xy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 10);

        auto g_y_0_xy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 11);

        auto g_y_0_xy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 12);

        auto g_y_0_xy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 13);

        auto g_y_0_xy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 14);

        auto g_y_0_xy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 15);

        auto g_y_0_xy_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 16);

        auto g_y_0_xy_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 17);

        auto g_y_0_xy_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 18);

        auto g_y_0_xy_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 19);

        auto g_y_0_xz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 20);

        auto g_y_0_xz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 21);

        auto g_y_0_xz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 22);

        auto g_y_0_xz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 23);

        auto g_y_0_xz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 24);

        auto g_y_0_xz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 25);

        auto g_y_0_xz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 26);

        auto g_y_0_xz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 27);

        auto g_y_0_xz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 28);

        auto g_y_0_xz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 29);

        auto g_y_0_yy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 30);

        auto g_y_0_yy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 31);

        auto g_y_0_yy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 32);

        auto g_y_0_yy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 33);

        auto g_y_0_yy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 34);

        auto g_y_0_yy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 35);

        auto g_y_0_yy_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 36);

        auto g_y_0_yy_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 37);

        auto g_y_0_yy_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 38);

        auto g_y_0_yy_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 39);

        auto g_y_0_yz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 40);

        auto g_y_0_yz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 41);

        auto g_y_0_yz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 42);

        auto g_y_0_yz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 43);

        auto g_y_0_yz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 44);

        auto g_y_0_yz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 45);

        auto g_y_0_yz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 46);

        auto g_y_0_yz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 47);

        auto g_y_0_yz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 48);

        auto g_y_0_yz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 49);

        auto g_y_0_zz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 50);

        auto g_y_0_zz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 51);

        auto g_y_0_zz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 52);

        auto g_y_0_zz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 53);

        auto g_y_0_zz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 54);

        auto g_y_0_zz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 55);

        auto g_y_0_zz_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 56);

        auto g_y_0_zz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 57);

        auto g_y_0_zz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 58);

        auto g_y_0_zz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 59);

        auto g_z_0_xx_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 0);

        auto g_z_0_xx_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 1);

        auto g_z_0_xx_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 2);

        auto g_z_0_xx_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 3);

        auto g_z_0_xx_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 4);

        auto g_z_0_xx_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 5);

        auto g_z_0_xx_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 6);

        auto g_z_0_xx_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 7);

        auto g_z_0_xx_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 8);

        auto g_z_0_xx_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 9);

        auto g_z_0_xy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 10);

        auto g_z_0_xy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 11);

        auto g_z_0_xy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 12);

        auto g_z_0_xy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 13);

        auto g_z_0_xy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 14);

        auto g_z_0_xy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 15);

        auto g_z_0_xy_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 16);

        auto g_z_0_xy_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 17);

        auto g_z_0_xy_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 18);

        auto g_z_0_xy_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 19);

        auto g_z_0_xz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 20);

        auto g_z_0_xz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 21);

        auto g_z_0_xz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 22);

        auto g_z_0_xz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 23);

        auto g_z_0_xz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 24);

        auto g_z_0_xz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 25);

        auto g_z_0_xz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 26);

        auto g_z_0_xz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 27);

        auto g_z_0_xz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 28);

        auto g_z_0_xz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 29);

        auto g_z_0_yy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 30);

        auto g_z_0_yy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 31);

        auto g_z_0_yy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 32);

        auto g_z_0_yy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 33);

        auto g_z_0_yy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 34);

        auto g_z_0_yy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 35);

        auto g_z_0_yy_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 36);

        auto g_z_0_yy_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 37);

        auto g_z_0_yy_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 38);

        auto g_z_0_yy_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 39);

        auto g_z_0_yz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 40);

        auto g_z_0_yz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 41);

        auto g_z_0_yz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 42);

        auto g_z_0_yz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 43);

        auto g_z_0_yz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 44);

        auto g_z_0_yz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 45);

        auto g_z_0_yz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 46);

        auto g_z_0_yz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 47);

        auto g_z_0_yz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 48);

        auto g_z_0_yz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 49);

        auto g_z_0_zz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 50);

        auto g_z_0_zz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 51);

        auto g_z_0_zz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 52);

        auto g_z_0_zz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 53);

        auto g_z_0_zz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 54);

        auto g_z_0_zz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 55);

        auto g_z_0_zz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 56);

        auto g_z_0_zz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 57);

        auto g_z_0_zz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 58);

        auto g_z_0_zz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 59);

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

        /// set up bra offset for contr_buffer_xxff

        const auto ff_geom_10_off = idx_geom_10_xff + i * 100;

        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_x_0_xx_xxx, g_x_0_xx_xxxx, g_x_0_xx_xxxy, g_x_0_xx_xxxz, g_x_0_xx_xxy, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxz, g_x_0_xx_xxzz, g_x_0_xx_xyy, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyz, g_x_0_xx_xyzz, g_x_0_xx_xzz, g_x_0_xx_xzzz, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yzz, g_x_0_xx_zzz, g_x_0_xxx_xxx, g_x_0_xxx_xxy, g_x_0_xxx_xxz, g_x_0_xxx_xyy, g_x_0_xxx_xyz, g_x_0_xxx_xzz, g_x_0_xxx_yyy, g_x_0_xxx_yyz, g_x_0_xxx_yzz, g_x_0_xxx_zzz, g_xx_xxx, g_xx_xxy, g_xx_xxz, g_xx_xyy, g_xx_xyz, g_xx_xzz, g_xx_yyy, g_xx_yyz, g_xx_yzz, g_xx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_xxx[k] = -g_xx_xxx[k] - g_x_0_xx_xxx[k] * cd_x[k] + g_x_0_xx_xxxx[k];

            g_x_0_xxx_xxy[k] = -g_xx_xxy[k] - g_x_0_xx_xxy[k] * cd_x[k] + g_x_0_xx_xxxy[k];

            g_x_0_xxx_xxz[k] = -g_xx_xxz[k] - g_x_0_xx_xxz[k] * cd_x[k] + g_x_0_xx_xxxz[k];

            g_x_0_xxx_xyy[k] = -g_xx_xyy[k] - g_x_0_xx_xyy[k] * cd_x[k] + g_x_0_xx_xxyy[k];

            g_x_0_xxx_xyz[k] = -g_xx_xyz[k] - g_x_0_xx_xyz[k] * cd_x[k] + g_x_0_xx_xxyz[k];

            g_x_0_xxx_xzz[k] = -g_xx_xzz[k] - g_x_0_xx_xzz[k] * cd_x[k] + g_x_0_xx_xxzz[k];

            g_x_0_xxx_yyy[k] = -g_xx_yyy[k] - g_x_0_xx_yyy[k] * cd_x[k] + g_x_0_xx_xyyy[k];

            g_x_0_xxx_yyz[k] = -g_xx_yyz[k] - g_x_0_xx_yyz[k] * cd_x[k] + g_x_0_xx_xyyz[k];

            g_x_0_xxx_yzz[k] = -g_xx_yzz[k] - g_x_0_xx_yzz[k] * cd_x[k] + g_x_0_xx_xyzz[k];

            g_x_0_xxx_zzz[k] = -g_xx_zzz[k] - g_x_0_xx_zzz[k] * cd_x[k] + g_x_0_xx_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 19);

        #pragma omp simd aligned(cd_y, g_x_0_xx_xxx, g_x_0_xx_xxxy, g_x_0_xx_xxy, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxz, g_x_0_xx_xyy, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyz, g_x_0_xx_xyzz, g_x_0_xx_xzz, g_x_0_xx_yyy, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyz, g_x_0_xx_yyzz, g_x_0_xx_yzz, g_x_0_xx_yzzz, g_x_0_xx_zzz, g_x_0_xxy_xxx, g_x_0_xxy_xxy, g_x_0_xxy_xxz, g_x_0_xxy_xyy, g_x_0_xxy_xyz, g_x_0_xxy_xzz, g_x_0_xxy_yyy, g_x_0_xxy_yyz, g_x_0_xxy_yzz, g_x_0_xxy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_xxx[k] = -g_x_0_xx_xxx[k] * cd_y[k] + g_x_0_xx_xxxy[k];

            g_x_0_xxy_xxy[k] = -g_x_0_xx_xxy[k] * cd_y[k] + g_x_0_xx_xxyy[k];

            g_x_0_xxy_xxz[k] = -g_x_0_xx_xxz[k] * cd_y[k] + g_x_0_xx_xxyz[k];

            g_x_0_xxy_xyy[k] = -g_x_0_xx_xyy[k] * cd_y[k] + g_x_0_xx_xyyy[k];

            g_x_0_xxy_xyz[k] = -g_x_0_xx_xyz[k] * cd_y[k] + g_x_0_xx_xyyz[k];

            g_x_0_xxy_xzz[k] = -g_x_0_xx_xzz[k] * cd_y[k] + g_x_0_xx_xyzz[k];

            g_x_0_xxy_yyy[k] = -g_x_0_xx_yyy[k] * cd_y[k] + g_x_0_xx_yyyy[k];

            g_x_0_xxy_yyz[k] = -g_x_0_xx_yyz[k] * cd_y[k] + g_x_0_xx_yyyz[k];

            g_x_0_xxy_yzz[k] = -g_x_0_xx_yzz[k] * cd_y[k] + g_x_0_xx_yyzz[k];

            g_x_0_xxy_zzz[k] = -g_x_0_xx_zzz[k] * cd_y[k] + g_x_0_xx_yzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_x_0_xx_xxx, g_x_0_xx_xxxz, g_x_0_xx_xxy, g_x_0_xx_xxyz, g_x_0_xx_xxz, g_x_0_xx_xxzz, g_x_0_xx_xyy, g_x_0_xx_xyyz, g_x_0_xx_xyz, g_x_0_xx_xyzz, g_x_0_xx_xzz, g_x_0_xx_xzzz, g_x_0_xx_yyy, g_x_0_xx_yyyz, g_x_0_xx_yyz, g_x_0_xx_yyzz, g_x_0_xx_yzz, g_x_0_xx_yzzz, g_x_0_xx_zzz, g_x_0_xx_zzzz, g_x_0_xxz_xxx, g_x_0_xxz_xxy, g_x_0_xxz_xxz, g_x_0_xxz_xyy, g_x_0_xxz_xyz, g_x_0_xxz_xzz, g_x_0_xxz_yyy, g_x_0_xxz_yyz, g_x_0_xxz_yzz, g_x_0_xxz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_xxx[k] = -g_x_0_xx_xxx[k] * cd_z[k] + g_x_0_xx_xxxz[k];

            g_x_0_xxz_xxy[k] = -g_x_0_xx_xxy[k] * cd_z[k] + g_x_0_xx_xxyz[k];

            g_x_0_xxz_xxz[k] = -g_x_0_xx_xxz[k] * cd_z[k] + g_x_0_xx_xxzz[k];

            g_x_0_xxz_xyy[k] = -g_x_0_xx_xyy[k] * cd_z[k] + g_x_0_xx_xyyz[k];

            g_x_0_xxz_xyz[k] = -g_x_0_xx_xyz[k] * cd_z[k] + g_x_0_xx_xyzz[k];

            g_x_0_xxz_xzz[k] = -g_x_0_xx_xzz[k] * cd_z[k] + g_x_0_xx_xzzz[k];

            g_x_0_xxz_yyy[k] = -g_x_0_xx_yyy[k] * cd_z[k] + g_x_0_xx_yyyz[k];

            g_x_0_xxz_yyz[k] = -g_x_0_xx_yyz[k] * cd_z[k] + g_x_0_xx_yyzz[k];

            g_x_0_xxz_yzz[k] = -g_x_0_xx_yzz[k] * cd_z[k] + g_x_0_xx_yzzz[k];

            g_x_0_xxz_zzz[k] = -g_x_0_xx_zzz[k] * cd_z[k] + g_x_0_xx_zzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 39);

        #pragma omp simd aligned(cd_y, g_x_0_xy_xxx, g_x_0_xy_xxxy, g_x_0_xy_xxy, g_x_0_xy_xxyy, g_x_0_xy_xxyz, g_x_0_xy_xxz, g_x_0_xy_xyy, g_x_0_xy_xyyy, g_x_0_xy_xyyz, g_x_0_xy_xyz, g_x_0_xy_xyzz, g_x_0_xy_xzz, g_x_0_xy_yyy, g_x_0_xy_yyyy, g_x_0_xy_yyyz, g_x_0_xy_yyz, g_x_0_xy_yyzz, g_x_0_xy_yzz, g_x_0_xy_yzzz, g_x_0_xy_zzz, g_x_0_xyy_xxx, g_x_0_xyy_xxy, g_x_0_xyy_xxz, g_x_0_xyy_xyy, g_x_0_xyy_xyz, g_x_0_xyy_xzz, g_x_0_xyy_yyy, g_x_0_xyy_yyz, g_x_0_xyy_yzz, g_x_0_xyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_xxx[k] = -g_x_0_xy_xxx[k] * cd_y[k] + g_x_0_xy_xxxy[k];

            g_x_0_xyy_xxy[k] = -g_x_0_xy_xxy[k] * cd_y[k] + g_x_0_xy_xxyy[k];

            g_x_0_xyy_xxz[k] = -g_x_0_xy_xxz[k] * cd_y[k] + g_x_0_xy_xxyz[k];

            g_x_0_xyy_xyy[k] = -g_x_0_xy_xyy[k] * cd_y[k] + g_x_0_xy_xyyy[k];

            g_x_0_xyy_xyz[k] = -g_x_0_xy_xyz[k] * cd_y[k] + g_x_0_xy_xyyz[k];

            g_x_0_xyy_xzz[k] = -g_x_0_xy_xzz[k] * cd_y[k] + g_x_0_xy_xyzz[k];

            g_x_0_xyy_yyy[k] = -g_x_0_xy_yyy[k] * cd_y[k] + g_x_0_xy_yyyy[k];

            g_x_0_xyy_yyz[k] = -g_x_0_xy_yyz[k] * cd_y[k] + g_x_0_xy_yyyz[k];

            g_x_0_xyy_yzz[k] = -g_x_0_xy_yzz[k] * cd_y[k] + g_x_0_xy_yyzz[k];

            g_x_0_xyy_zzz[k] = -g_x_0_xy_zzz[k] * cd_y[k] + g_x_0_xy_yzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 49);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_xxx, g_x_0_xyz_xxy, g_x_0_xyz_xxz, g_x_0_xyz_xyy, g_x_0_xyz_xyz, g_x_0_xyz_xzz, g_x_0_xyz_yyy, g_x_0_xyz_yyz, g_x_0_xyz_yzz, g_x_0_xyz_zzz, g_x_0_xz_xxx, g_x_0_xz_xxxy, g_x_0_xz_xxy, g_x_0_xz_xxyy, g_x_0_xz_xxyz, g_x_0_xz_xxz, g_x_0_xz_xyy, g_x_0_xz_xyyy, g_x_0_xz_xyyz, g_x_0_xz_xyz, g_x_0_xz_xyzz, g_x_0_xz_xzz, g_x_0_xz_yyy, g_x_0_xz_yyyy, g_x_0_xz_yyyz, g_x_0_xz_yyz, g_x_0_xz_yyzz, g_x_0_xz_yzz, g_x_0_xz_yzzz, g_x_0_xz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_xxx[k] = -g_x_0_xz_xxx[k] * cd_y[k] + g_x_0_xz_xxxy[k];

            g_x_0_xyz_xxy[k] = -g_x_0_xz_xxy[k] * cd_y[k] + g_x_0_xz_xxyy[k];

            g_x_0_xyz_xxz[k] = -g_x_0_xz_xxz[k] * cd_y[k] + g_x_0_xz_xxyz[k];

            g_x_0_xyz_xyy[k] = -g_x_0_xz_xyy[k] * cd_y[k] + g_x_0_xz_xyyy[k];

            g_x_0_xyz_xyz[k] = -g_x_0_xz_xyz[k] * cd_y[k] + g_x_0_xz_xyyz[k];

            g_x_0_xyz_xzz[k] = -g_x_0_xz_xzz[k] * cd_y[k] + g_x_0_xz_xyzz[k];

            g_x_0_xyz_yyy[k] = -g_x_0_xz_yyy[k] * cd_y[k] + g_x_0_xz_yyyy[k];

            g_x_0_xyz_yyz[k] = -g_x_0_xz_yyz[k] * cd_y[k] + g_x_0_xz_yyyz[k];

            g_x_0_xyz_yzz[k] = -g_x_0_xz_yzz[k] * cd_y[k] + g_x_0_xz_yyzz[k];

            g_x_0_xyz_zzz[k] = -g_x_0_xz_zzz[k] * cd_y[k] + g_x_0_xz_yzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_x_0_xz_xxx, g_x_0_xz_xxxz, g_x_0_xz_xxy, g_x_0_xz_xxyz, g_x_0_xz_xxz, g_x_0_xz_xxzz, g_x_0_xz_xyy, g_x_0_xz_xyyz, g_x_0_xz_xyz, g_x_0_xz_xyzz, g_x_0_xz_xzz, g_x_0_xz_xzzz, g_x_0_xz_yyy, g_x_0_xz_yyyz, g_x_0_xz_yyz, g_x_0_xz_yyzz, g_x_0_xz_yzz, g_x_0_xz_yzzz, g_x_0_xz_zzz, g_x_0_xz_zzzz, g_x_0_xzz_xxx, g_x_0_xzz_xxy, g_x_0_xzz_xxz, g_x_0_xzz_xyy, g_x_0_xzz_xyz, g_x_0_xzz_xzz, g_x_0_xzz_yyy, g_x_0_xzz_yyz, g_x_0_xzz_yzz, g_x_0_xzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_xxx[k] = -g_x_0_xz_xxx[k] * cd_z[k] + g_x_0_xz_xxxz[k];

            g_x_0_xzz_xxy[k] = -g_x_0_xz_xxy[k] * cd_z[k] + g_x_0_xz_xxyz[k];

            g_x_0_xzz_xxz[k] = -g_x_0_xz_xxz[k] * cd_z[k] + g_x_0_xz_xxzz[k];

            g_x_0_xzz_xyy[k] = -g_x_0_xz_xyy[k] * cd_z[k] + g_x_0_xz_xyyz[k];

            g_x_0_xzz_xyz[k] = -g_x_0_xz_xyz[k] * cd_z[k] + g_x_0_xz_xyzz[k];

            g_x_0_xzz_xzz[k] = -g_x_0_xz_xzz[k] * cd_z[k] + g_x_0_xz_xzzz[k];

            g_x_0_xzz_yyy[k] = -g_x_0_xz_yyy[k] * cd_z[k] + g_x_0_xz_yyyz[k];

            g_x_0_xzz_yyz[k] = -g_x_0_xz_yyz[k] * cd_z[k] + g_x_0_xz_yyzz[k];

            g_x_0_xzz_yzz[k] = -g_x_0_xz_yzz[k] * cd_z[k] + g_x_0_xz_yzzz[k];

            g_x_0_xzz_zzz[k] = -g_x_0_xz_zzz[k] * cd_z[k] + g_x_0_xz_zzzz[k];
        }

        /// Set up 60-70 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 69);

        #pragma omp simd aligned(cd_y, g_x_0_yy_xxx, g_x_0_yy_xxxy, g_x_0_yy_xxy, g_x_0_yy_xxyy, g_x_0_yy_xxyz, g_x_0_yy_xxz, g_x_0_yy_xyy, g_x_0_yy_xyyy, g_x_0_yy_xyyz, g_x_0_yy_xyz, g_x_0_yy_xyzz, g_x_0_yy_xzz, g_x_0_yy_yyy, g_x_0_yy_yyyy, g_x_0_yy_yyyz, g_x_0_yy_yyz, g_x_0_yy_yyzz, g_x_0_yy_yzz, g_x_0_yy_yzzz, g_x_0_yy_zzz, g_x_0_yyy_xxx, g_x_0_yyy_xxy, g_x_0_yyy_xxz, g_x_0_yyy_xyy, g_x_0_yyy_xyz, g_x_0_yyy_xzz, g_x_0_yyy_yyy, g_x_0_yyy_yyz, g_x_0_yyy_yzz, g_x_0_yyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_xxx[k] = -g_x_0_yy_xxx[k] * cd_y[k] + g_x_0_yy_xxxy[k];

            g_x_0_yyy_xxy[k] = -g_x_0_yy_xxy[k] * cd_y[k] + g_x_0_yy_xxyy[k];

            g_x_0_yyy_xxz[k] = -g_x_0_yy_xxz[k] * cd_y[k] + g_x_0_yy_xxyz[k];

            g_x_0_yyy_xyy[k] = -g_x_0_yy_xyy[k] * cd_y[k] + g_x_0_yy_xyyy[k];

            g_x_0_yyy_xyz[k] = -g_x_0_yy_xyz[k] * cd_y[k] + g_x_0_yy_xyyz[k];

            g_x_0_yyy_xzz[k] = -g_x_0_yy_xzz[k] * cd_y[k] + g_x_0_yy_xyzz[k];

            g_x_0_yyy_yyy[k] = -g_x_0_yy_yyy[k] * cd_y[k] + g_x_0_yy_yyyy[k];

            g_x_0_yyy_yyz[k] = -g_x_0_yy_yyz[k] * cd_y[k] + g_x_0_yy_yyyz[k];

            g_x_0_yyy_yzz[k] = -g_x_0_yy_yzz[k] * cd_y[k] + g_x_0_yy_yyzz[k];

            g_x_0_yyy_zzz[k] = -g_x_0_yy_zzz[k] * cd_y[k] + g_x_0_yy_yzzz[k];
        }

        /// Set up 70-80 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 74);

        auto g_x_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 79);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_xxx, g_x_0_yyz_xxy, g_x_0_yyz_xxz, g_x_0_yyz_xyy, g_x_0_yyz_xyz, g_x_0_yyz_xzz, g_x_0_yyz_yyy, g_x_0_yyz_yyz, g_x_0_yyz_yzz, g_x_0_yyz_zzz, g_x_0_yz_xxx, g_x_0_yz_xxxy, g_x_0_yz_xxy, g_x_0_yz_xxyy, g_x_0_yz_xxyz, g_x_0_yz_xxz, g_x_0_yz_xyy, g_x_0_yz_xyyy, g_x_0_yz_xyyz, g_x_0_yz_xyz, g_x_0_yz_xyzz, g_x_0_yz_xzz, g_x_0_yz_yyy, g_x_0_yz_yyyy, g_x_0_yz_yyyz, g_x_0_yz_yyz, g_x_0_yz_yyzz, g_x_0_yz_yzz, g_x_0_yz_yzzz, g_x_0_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_xxx[k] = -g_x_0_yz_xxx[k] * cd_y[k] + g_x_0_yz_xxxy[k];

            g_x_0_yyz_xxy[k] = -g_x_0_yz_xxy[k] * cd_y[k] + g_x_0_yz_xxyy[k];

            g_x_0_yyz_xxz[k] = -g_x_0_yz_xxz[k] * cd_y[k] + g_x_0_yz_xxyz[k];

            g_x_0_yyz_xyy[k] = -g_x_0_yz_xyy[k] * cd_y[k] + g_x_0_yz_xyyy[k];

            g_x_0_yyz_xyz[k] = -g_x_0_yz_xyz[k] * cd_y[k] + g_x_0_yz_xyyz[k];

            g_x_0_yyz_xzz[k] = -g_x_0_yz_xzz[k] * cd_y[k] + g_x_0_yz_xyzz[k];

            g_x_0_yyz_yyy[k] = -g_x_0_yz_yyy[k] * cd_y[k] + g_x_0_yz_yyyy[k];

            g_x_0_yyz_yyz[k] = -g_x_0_yz_yyz[k] * cd_y[k] + g_x_0_yz_yyyz[k];

            g_x_0_yyz_yzz[k] = -g_x_0_yz_yzz[k] * cd_y[k] + g_x_0_yz_yyzz[k];

            g_x_0_yyz_zzz[k] = -g_x_0_yz_zzz[k] * cd_y[k] + g_x_0_yz_yzzz[k];
        }

        /// Set up 80-90 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 83);

        auto g_x_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 89);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_xxx, g_x_0_yzz_xxy, g_x_0_yzz_xxz, g_x_0_yzz_xyy, g_x_0_yzz_xyz, g_x_0_yzz_xzz, g_x_0_yzz_yyy, g_x_0_yzz_yyz, g_x_0_yzz_yzz, g_x_0_yzz_zzz, g_x_0_zz_xxx, g_x_0_zz_xxxy, g_x_0_zz_xxy, g_x_0_zz_xxyy, g_x_0_zz_xxyz, g_x_0_zz_xxz, g_x_0_zz_xyy, g_x_0_zz_xyyy, g_x_0_zz_xyyz, g_x_0_zz_xyz, g_x_0_zz_xyzz, g_x_0_zz_xzz, g_x_0_zz_yyy, g_x_0_zz_yyyy, g_x_0_zz_yyyz, g_x_0_zz_yyz, g_x_0_zz_yyzz, g_x_0_zz_yzz, g_x_0_zz_yzzz, g_x_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_xxx[k] = -g_x_0_zz_xxx[k] * cd_y[k] + g_x_0_zz_xxxy[k];

            g_x_0_yzz_xxy[k] = -g_x_0_zz_xxy[k] * cd_y[k] + g_x_0_zz_xxyy[k];

            g_x_0_yzz_xxz[k] = -g_x_0_zz_xxz[k] * cd_y[k] + g_x_0_zz_xxyz[k];

            g_x_0_yzz_xyy[k] = -g_x_0_zz_xyy[k] * cd_y[k] + g_x_0_zz_xyyy[k];

            g_x_0_yzz_xyz[k] = -g_x_0_zz_xyz[k] * cd_y[k] + g_x_0_zz_xyyz[k];

            g_x_0_yzz_xzz[k] = -g_x_0_zz_xzz[k] * cd_y[k] + g_x_0_zz_xyzz[k];

            g_x_0_yzz_yyy[k] = -g_x_0_zz_yyy[k] * cd_y[k] + g_x_0_zz_yyyy[k];

            g_x_0_yzz_yyz[k] = -g_x_0_zz_yyz[k] * cd_y[k] + g_x_0_zz_yyyz[k];

            g_x_0_yzz_yzz[k] = -g_x_0_zz_yzz[k] * cd_y[k] + g_x_0_zz_yyzz[k];

            g_x_0_yzz_zzz[k] = -g_x_0_zz_zzz[k] * cd_y[k] + g_x_0_zz_yzzz[k];
        }

        /// Set up 90-100 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps  + 99);

        #pragma omp simd aligned(cd_z, g_x_0_zz_xxx, g_x_0_zz_xxxz, g_x_0_zz_xxy, g_x_0_zz_xxyz, g_x_0_zz_xxz, g_x_0_zz_xxzz, g_x_0_zz_xyy, g_x_0_zz_xyyz, g_x_0_zz_xyz, g_x_0_zz_xyzz, g_x_0_zz_xzz, g_x_0_zz_xzzz, g_x_0_zz_yyy, g_x_0_zz_yyyz, g_x_0_zz_yyz, g_x_0_zz_yyzz, g_x_0_zz_yzz, g_x_0_zz_yzzz, g_x_0_zz_zzz, g_x_0_zz_zzzz, g_x_0_zzz_xxx, g_x_0_zzz_xxy, g_x_0_zzz_xxz, g_x_0_zzz_xyy, g_x_0_zzz_xyz, g_x_0_zzz_xzz, g_x_0_zzz_yyy, g_x_0_zzz_yyz, g_x_0_zzz_yzz, g_x_0_zzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_xxx[k] = -g_x_0_zz_xxx[k] * cd_z[k] + g_x_0_zz_xxxz[k];

            g_x_0_zzz_xxy[k] = -g_x_0_zz_xxy[k] * cd_z[k] + g_x_0_zz_xxyz[k];

            g_x_0_zzz_xxz[k] = -g_x_0_zz_xxz[k] * cd_z[k] + g_x_0_zz_xxzz[k];

            g_x_0_zzz_xyy[k] = -g_x_0_zz_xyy[k] * cd_z[k] + g_x_0_zz_xyyz[k];

            g_x_0_zzz_xyz[k] = -g_x_0_zz_xyz[k] * cd_z[k] + g_x_0_zz_xyzz[k];

            g_x_0_zzz_xzz[k] = -g_x_0_zz_xzz[k] * cd_z[k] + g_x_0_zz_xzzz[k];

            g_x_0_zzz_yyy[k] = -g_x_0_zz_yyy[k] * cd_z[k] + g_x_0_zz_yyyz[k];

            g_x_0_zzz_yyz[k] = -g_x_0_zz_yyz[k] * cd_z[k] + g_x_0_zz_yyzz[k];

            g_x_0_zzz_yzz[k] = -g_x_0_zz_yzz[k] * cd_z[k] + g_x_0_zz_yzzz[k];

            g_x_0_zzz_zzz[k] = -g_x_0_zz_zzz[k] * cd_z[k] + g_x_0_zz_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 0);

        auto g_y_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 1);

        auto g_y_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 2);

        auto g_y_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 3);

        auto g_y_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 4);

        auto g_y_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 5);

        auto g_y_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 6);

        auto g_y_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 7);

        auto g_y_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 8);

        auto g_y_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_y_0_xx_xxx, g_y_0_xx_xxxx, g_y_0_xx_xxxy, g_y_0_xx_xxxz, g_y_0_xx_xxy, g_y_0_xx_xxyy, g_y_0_xx_xxyz, g_y_0_xx_xxz, g_y_0_xx_xxzz, g_y_0_xx_xyy, g_y_0_xx_xyyy, g_y_0_xx_xyyz, g_y_0_xx_xyz, g_y_0_xx_xyzz, g_y_0_xx_xzz, g_y_0_xx_xzzz, g_y_0_xx_yyy, g_y_0_xx_yyz, g_y_0_xx_yzz, g_y_0_xx_zzz, g_y_0_xxx_xxx, g_y_0_xxx_xxy, g_y_0_xxx_xxz, g_y_0_xxx_xyy, g_y_0_xxx_xyz, g_y_0_xxx_xzz, g_y_0_xxx_yyy, g_y_0_xxx_yyz, g_y_0_xxx_yzz, g_y_0_xxx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_xxx[k] = -g_y_0_xx_xxx[k] * cd_x[k] + g_y_0_xx_xxxx[k];

            g_y_0_xxx_xxy[k] = -g_y_0_xx_xxy[k] * cd_x[k] + g_y_0_xx_xxxy[k];

            g_y_0_xxx_xxz[k] = -g_y_0_xx_xxz[k] * cd_x[k] + g_y_0_xx_xxxz[k];

            g_y_0_xxx_xyy[k] = -g_y_0_xx_xyy[k] * cd_x[k] + g_y_0_xx_xxyy[k];

            g_y_0_xxx_xyz[k] = -g_y_0_xx_xyz[k] * cd_x[k] + g_y_0_xx_xxyz[k];

            g_y_0_xxx_xzz[k] = -g_y_0_xx_xzz[k] * cd_x[k] + g_y_0_xx_xxzz[k];

            g_y_0_xxx_yyy[k] = -g_y_0_xx_yyy[k] * cd_x[k] + g_y_0_xx_xyyy[k];

            g_y_0_xxx_yyz[k] = -g_y_0_xx_yyz[k] * cd_x[k] + g_y_0_xx_xyyz[k];

            g_y_0_xxx_yzz[k] = -g_y_0_xx_yzz[k] * cd_x[k] + g_y_0_xx_xyzz[k];

            g_y_0_xxx_zzz[k] = -g_y_0_xx_zzz[k] * cd_x[k] + g_y_0_xx_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 10);

        auto g_y_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 11);

        auto g_y_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 12);

        auto g_y_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 13);

        auto g_y_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 14);

        auto g_y_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 15);

        auto g_y_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 16);

        auto g_y_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 17);

        auto g_y_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 18);

        auto g_y_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 19);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_xxx, g_y_0_xxy_xxy, g_y_0_xxy_xxz, g_y_0_xxy_xyy, g_y_0_xxy_xyz, g_y_0_xxy_xzz, g_y_0_xxy_yyy, g_y_0_xxy_yyz, g_y_0_xxy_yzz, g_y_0_xxy_zzz, g_y_0_xy_xxx, g_y_0_xy_xxxx, g_y_0_xy_xxxy, g_y_0_xy_xxxz, g_y_0_xy_xxy, g_y_0_xy_xxyy, g_y_0_xy_xxyz, g_y_0_xy_xxz, g_y_0_xy_xxzz, g_y_0_xy_xyy, g_y_0_xy_xyyy, g_y_0_xy_xyyz, g_y_0_xy_xyz, g_y_0_xy_xyzz, g_y_0_xy_xzz, g_y_0_xy_xzzz, g_y_0_xy_yyy, g_y_0_xy_yyz, g_y_0_xy_yzz, g_y_0_xy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_xxx[k] = -g_y_0_xy_xxx[k] * cd_x[k] + g_y_0_xy_xxxx[k];

            g_y_0_xxy_xxy[k] = -g_y_0_xy_xxy[k] * cd_x[k] + g_y_0_xy_xxxy[k];

            g_y_0_xxy_xxz[k] = -g_y_0_xy_xxz[k] * cd_x[k] + g_y_0_xy_xxxz[k];

            g_y_0_xxy_xyy[k] = -g_y_0_xy_xyy[k] * cd_x[k] + g_y_0_xy_xxyy[k];

            g_y_0_xxy_xyz[k] = -g_y_0_xy_xyz[k] * cd_x[k] + g_y_0_xy_xxyz[k];

            g_y_0_xxy_xzz[k] = -g_y_0_xy_xzz[k] * cd_x[k] + g_y_0_xy_xxzz[k];

            g_y_0_xxy_yyy[k] = -g_y_0_xy_yyy[k] * cd_x[k] + g_y_0_xy_xyyy[k];

            g_y_0_xxy_yyz[k] = -g_y_0_xy_yyz[k] * cd_x[k] + g_y_0_xy_xyyz[k];

            g_y_0_xxy_yzz[k] = -g_y_0_xy_yzz[k] * cd_x[k] + g_y_0_xy_xyzz[k];

            g_y_0_xxy_zzz[k] = -g_y_0_xy_zzz[k] * cd_x[k] + g_y_0_xy_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 20);

        auto g_y_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 21);

        auto g_y_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 22);

        auto g_y_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 23);

        auto g_y_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 24);

        auto g_y_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 25);

        auto g_y_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 26);

        auto g_y_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 27);

        auto g_y_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 28);

        auto g_y_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_xxx, g_y_0_xxz_xxy, g_y_0_xxz_xxz, g_y_0_xxz_xyy, g_y_0_xxz_xyz, g_y_0_xxz_xzz, g_y_0_xxz_yyy, g_y_0_xxz_yyz, g_y_0_xxz_yzz, g_y_0_xxz_zzz, g_y_0_xz_xxx, g_y_0_xz_xxxx, g_y_0_xz_xxxy, g_y_0_xz_xxxz, g_y_0_xz_xxy, g_y_0_xz_xxyy, g_y_0_xz_xxyz, g_y_0_xz_xxz, g_y_0_xz_xxzz, g_y_0_xz_xyy, g_y_0_xz_xyyy, g_y_0_xz_xyyz, g_y_0_xz_xyz, g_y_0_xz_xyzz, g_y_0_xz_xzz, g_y_0_xz_xzzz, g_y_0_xz_yyy, g_y_0_xz_yyz, g_y_0_xz_yzz, g_y_0_xz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_xxx[k] = -g_y_0_xz_xxx[k] * cd_x[k] + g_y_0_xz_xxxx[k];

            g_y_0_xxz_xxy[k] = -g_y_0_xz_xxy[k] * cd_x[k] + g_y_0_xz_xxxy[k];

            g_y_0_xxz_xxz[k] = -g_y_0_xz_xxz[k] * cd_x[k] + g_y_0_xz_xxxz[k];

            g_y_0_xxz_xyy[k] = -g_y_0_xz_xyy[k] * cd_x[k] + g_y_0_xz_xxyy[k];

            g_y_0_xxz_xyz[k] = -g_y_0_xz_xyz[k] * cd_x[k] + g_y_0_xz_xxyz[k];

            g_y_0_xxz_xzz[k] = -g_y_0_xz_xzz[k] * cd_x[k] + g_y_0_xz_xxzz[k];

            g_y_0_xxz_yyy[k] = -g_y_0_xz_yyy[k] * cd_x[k] + g_y_0_xz_xyyy[k];

            g_y_0_xxz_yyz[k] = -g_y_0_xz_yyz[k] * cd_x[k] + g_y_0_xz_xyyz[k];

            g_y_0_xxz_yzz[k] = -g_y_0_xz_yzz[k] * cd_x[k] + g_y_0_xz_xyzz[k];

            g_y_0_xxz_zzz[k] = -g_y_0_xz_zzz[k] * cd_x[k] + g_y_0_xz_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 30);

        auto g_y_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 31);

        auto g_y_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 32);

        auto g_y_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 33);

        auto g_y_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 34);

        auto g_y_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 35);

        auto g_y_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 36);

        auto g_y_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 37);

        auto g_y_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 38);

        auto g_y_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 39);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_xxx, g_y_0_xyy_xxy, g_y_0_xyy_xxz, g_y_0_xyy_xyy, g_y_0_xyy_xyz, g_y_0_xyy_xzz, g_y_0_xyy_yyy, g_y_0_xyy_yyz, g_y_0_xyy_yzz, g_y_0_xyy_zzz, g_y_0_yy_xxx, g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxy, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxz, g_y_0_yy_xxzz, g_y_0_yy_xyy, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyz, g_y_0_yy_xyzz, g_y_0_yy_xzz, g_y_0_yy_xzzz, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yzz, g_y_0_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_xxx[k] = -g_y_0_yy_xxx[k] * cd_x[k] + g_y_0_yy_xxxx[k];

            g_y_0_xyy_xxy[k] = -g_y_0_yy_xxy[k] * cd_x[k] + g_y_0_yy_xxxy[k];

            g_y_0_xyy_xxz[k] = -g_y_0_yy_xxz[k] * cd_x[k] + g_y_0_yy_xxxz[k];

            g_y_0_xyy_xyy[k] = -g_y_0_yy_xyy[k] * cd_x[k] + g_y_0_yy_xxyy[k];

            g_y_0_xyy_xyz[k] = -g_y_0_yy_xyz[k] * cd_x[k] + g_y_0_yy_xxyz[k];

            g_y_0_xyy_xzz[k] = -g_y_0_yy_xzz[k] * cd_x[k] + g_y_0_yy_xxzz[k];

            g_y_0_xyy_yyy[k] = -g_y_0_yy_yyy[k] * cd_x[k] + g_y_0_yy_xyyy[k];

            g_y_0_xyy_yyz[k] = -g_y_0_yy_yyz[k] * cd_x[k] + g_y_0_yy_xyyz[k];

            g_y_0_xyy_yzz[k] = -g_y_0_yy_yzz[k] * cd_x[k] + g_y_0_yy_xyzz[k];

            g_y_0_xyy_zzz[k] = -g_y_0_yy_zzz[k] * cd_x[k] + g_y_0_yy_xzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 40);

        auto g_y_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 41);

        auto g_y_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 42);

        auto g_y_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 43);

        auto g_y_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 44);

        auto g_y_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 45);

        auto g_y_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 46);

        auto g_y_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 47);

        auto g_y_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 48);

        auto g_y_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 49);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_xxx, g_y_0_xyz_xxy, g_y_0_xyz_xxz, g_y_0_xyz_xyy, g_y_0_xyz_xyz, g_y_0_xyz_xzz, g_y_0_xyz_yyy, g_y_0_xyz_yyz, g_y_0_xyz_yzz, g_y_0_xyz_zzz, g_y_0_yz_xxx, g_y_0_yz_xxxx, g_y_0_yz_xxxy, g_y_0_yz_xxxz, g_y_0_yz_xxy, g_y_0_yz_xxyy, g_y_0_yz_xxyz, g_y_0_yz_xxz, g_y_0_yz_xxzz, g_y_0_yz_xyy, g_y_0_yz_xyyy, g_y_0_yz_xyyz, g_y_0_yz_xyz, g_y_0_yz_xyzz, g_y_0_yz_xzz, g_y_0_yz_xzzz, g_y_0_yz_yyy, g_y_0_yz_yyz, g_y_0_yz_yzz, g_y_0_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_xxx[k] = -g_y_0_yz_xxx[k] * cd_x[k] + g_y_0_yz_xxxx[k];

            g_y_0_xyz_xxy[k] = -g_y_0_yz_xxy[k] * cd_x[k] + g_y_0_yz_xxxy[k];

            g_y_0_xyz_xxz[k] = -g_y_0_yz_xxz[k] * cd_x[k] + g_y_0_yz_xxxz[k];

            g_y_0_xyz_xyy[k] = -g_y_0_yz_xyy[k] * cd_x[k] + g_y_0_yz_xxyy[k];

            g_y_0_xyz_xyz[k] = -g_y_0_yz_xyz[k] * cd_x[k] + g_y_0_yz_xxyz[k];

            g_y_0_xyz_xzz[k] = -g_y_0_yz_xzz[k] * cd_x[k] + g_y_0_yz_xxzz[k];

            g_y_0_xyz_yyy[k] = -g_y_0_yz_yyy[k] * cd_x[k] + g_y_0_yz_xyyy[k];

            g_y_0_xyz_yyz[k] = -g_y_0_yz_yyz[k] * cd_x[k] + g_y_0_yz_xyyz[k];

            g_y_0_xyz_yzz[k] = -g_y_0_yz_yzz[k] * cd_x[k] + g_y_0_yz_xyzz[k];

            g_y_0_xyz_zzz[k] = -g_y_0_yz_zzz[k] * cd_x[k] + g_y_0_yz_xzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 50);

        auto g_y_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 51);

        auto g_y_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 52);

        auto g_y_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 53);

        auto g_y_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 54);

        auto g_y_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 55);

        auto g_y_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 56);

        auto g_y_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 57);

        auto g_y_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 58);

        auto g_y_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_xxx, g_y_0_xzz_xxy, g_y_0_xzz_xxz, g_y_0_xzz_xyy, g_y_0_xzz_xyz, g_y_0_xzz_xzz, g_y_0_xzz_yyy, g_y_0_xzz_yyz, g_y_0_xzz_yzz, g_y_0_xzz_zzz, g_y_0_zz_xxx, g_y_0_zz_xxxx, g_y_0_zz_xxxy, g_y_0_zz_xxxz, g_y_0_zz_xxy, g_y_0_zz_xxyy, g_y_0_zz_xxyz, g_y_0_zz_xxz, g_y_0_zz_xxzz, g_y_0_zz_xyy, g_y_0_zz_xyyy, g_y_0_zz_xyyz, g_y_0_zz_xyz, g_y_0_zz_xyzz, g_y_0_zz_xzz, g_y_0_zz_xzzz, g_y_0_zz_yyy, g_y_0_zz_yyz, g_y_0_zz_yzz, g_y_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_xxx[k] = -g_y_0_zz_xxx[k] * cd_x[k] + g_y_0_zz_xxxx[k];

            g_y_0_xzz_xxy[k] = -g_y_0_zz_xxy[k] * cd_x[k] + g_y_0_zz_xxxy[k];

            g_y_0_xzz_xxz[k] = -g_y_0_zz_xxz[k] * cd_x[k] + g_y_0_zz_xxxz[k];

            g_y_0_xzz_xyy[k] = -g_y_0_zz_xyy[k] * cd_x[k] + g_y_0_zz_xxyy[k];

            g_y_0_xzz_xyz[k] = -g_y_0_zz_xyz[k] * cd_x[k] + g_y_0_zz_xxyz[k];

            g_y_0_xzz_xzz[k] = -g_y_0_zz_xzz[k] * cd_x[k] + g_y_0_zz_xxzz[k];

            g_y_0_xzz_yyy[k] = -g_y_0_zz_yyy[k] * cd_x[k] + g_y_0_zz_xyyy[k];

            g_y_0_xzz_yyz[k] = -g_y_0_zz_yyz[k] * cd_x[k] + g_y_0_zz_xyyz[k];

            g_y_0_xzz_yzz[k] = -g_y_0_zz_yzz[k] * cd_x[k] + g_y_0_zz_xyzz[k];

            g_y_0_xzz_zzz[k] = -g_y_0_zz_zzz[k] * cd_x[k] + g_y_0_zz_xzzz[k];
        }

        /// Set up 60-70 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 60);

        auto g_y_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 61);

        auto g_y_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 62);

        auto g_y_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 63);

        auto g_y_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 64);

        auto g_y_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 65);

        auto g_y_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 66);

        auto g_y_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 67);

        auto g_y_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 68);

        auto g_y_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 69);

        #pragma omp simd aligned(cd_y, g_y_0_yy_xxx, g_y_0_yy_xxxy, g_y_0_yy_xxy, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxz, g_y_0_yy_xyy, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyz, g_y_0_yy_xyzz, g_y_0_yy_xzz, g_y_0_yy_yyy, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyz, g_y_0_yy_yyzz, g_y_0_yy_yzz, g_y_0_yy_yzzz, g_y_0_yy_zzz, g_y_0_yyy_xxx, g_y_0_yyy_xxy, g_y_0_yyy_xxz, g_y_0_yyy_xyy, g_y_0_yyy_xyz, g_y_0_yyy_xzz, g_y_0_yyy_yyy, g_y_0_yyy_yyz, g_y_0_yyy_yzz, g_y_0_yyy_zzz, g_yy_xxx, g_yy_xxy, g_yy_xxz, g_yy_xyy, g_yy_xyz, g_yy_xzz, g_yy_yyy, g_yy_yyz, g_yy_yzz, g_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_xxx[k] = -g_yy_xxx[k] - g_y_0_yy_xxx[k] * cd_y[k] + g_y_0_yy_xxxy[k];

            g_y_0_yyy_xxy[k] = -g_yy_xxy[k] - g_y_0_yy_xxy[k] * cd_y[k] + g_y_0_yy_xxyy[k];

            g_y_0_yyy_xxz[k] = -g_yy_xxz[k] - g_y_0_yy_xxz[k] * cd_y[k] + g_y_0_yy_xxyz[k];

            g_y_0_yyy_xyy[k] = -g_yy_xyy[k] - g_y_0_yy_xyy[k] * cd_y[k] + g_y_0_yy_xyyy[k];

            g_y_0_yyy_xyz[k] = -g_yy_xyz[k] - g_y_0_yy_xyz[k] * cd_y[k] + g_y_0_yy_xyyz[k];

            g_y_0_yyy_xzz[k] = -g_yy_xzz[k] - g_y_0_yy_xzz[k] * cd_y[k] + g_y_0_yy_xyzz[k];

            g_y_0_yyy_yyy[k] = -g_yy_yyy[k] - g_y_0_yy_yyy[k] * cd_y[k] + g_y_0_yy_yyyy[k];

            g_y_0_yyy_yyz[k] = -g_yy_yyz[k] - g_y_0_yy_yyz[k] * cd_y[k] + g_y_0_yy_yyyz[k];

            g_y_0_yyy_yzz[k] = -g_yy_yzz[k] - g_y_0_yy_yzz[k] * cd_y[k] + g_y_0_yy_yyzz[k];

            g_y_0_yyy_zzz[k] = -g_yy_zzz[k] - g_y_0_yy_zzz[k] * cd_y[k] + g_y_0_yy_yzzz[k];
        }

        /// Set up 70-80 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 70);

        auto g_y_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 71);

        auto g_y_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 72);

        auto g_y_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 73);

        auto g_y_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 74);

        auto g_y_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 75);

        auto g_y_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 76);

        auto g_y_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 77);

        auto g_y_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 78);

        auto g_y_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 79);

        #pragma omp simd aligned(cd_z, g_y_0_yy_xxx, g_y_0_yy_xxxz, g_y_0_yy_xxy, g_y_0_yy_xxyz, g_y_0_yy_xxz, g_y_0_yy_xxzz, g_y_0_yy_xyy, g_y_0_yy_xyyz, g_y_0_yy_xyz, g_y_0_yy_xyzz, g_y_0_yy_xzz, g_y_0_yy_xzzz, g_y_0_yy_yyy, g_y_0_yy_yyyz, g_y_0_yy_yyz, g_y_0_yy_yyzz, g_y_0_yy_yzz, g_y_0_yy_yzzz, g_y_0_yy_zzz, g_y_0_yy_zzzz, g_y_0_yyz_xxx, g_y_0_yyz_xxy, g_y_0_yyz_xxz, g_y_0_yyz_xyy, g_y_0_yyz_xyz, g_y_0_yyz_xzz, g_y_0_yyz_yyy, g_y_0_yyz_yyz, g_y_0_yyz_yzz, g_y_0_yyz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_xxx[k] = -g_y_0_yy_xxx[k] * cd_z[k] + g_y_0_yy_xxxz[k];

            g_y_0_yyz_xxy[k] = -g_y_0_yy_xxy[k] * cd_z[k] + g_y_0_yy_xxyz[k];

            g_y_0_yyz_xxz[k] = -g_y_0_yy_xxz[k] * cd_z[k] + g_y_0_yy_xxzz[k];

            g_y_0_yyz_xyy[k] = -g_y_0_yy_xyy[k] * cd_z[k] + g_y_0_yy_xyyz[k];

            g_y_0_yyz_xyz[k] = -g_y_0_yy_xyz[k] * cd_z[k] + g_y_0_yy_xyzz[k];

            g_y_0_yyz_xzz[k] = -g_y_0_yy_xzz[k] * cd_z[k] + g_y_0_yy_xzzz[k];

            g_y_0_yyz_yyy[k] = -g_y_0_yy_yyy[k] * cd_z[k] + g_y_0_yy_yyyz[k];

            g_y_0_yyz_yyz[k] = -g_y_0_yy_yyz[k] * cd_z[k] + g_y_0_yy_yyzz[k];

            g_y_0_yyz_yzz[k] = -g_y_0_yy_yzz[k] * cd_z[k] + g_y_0_yy_yzzz[k];

            g_y_0_yyz_zzz[k] = -g_y_0_yy_zzz[k] * cd_z[k] + g_y_0_yy_zzzz[k];
        }

        /// Set up 80-90 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 80);

        auto g_y_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 81);

        auto g_y_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 82);

        auto g_y_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 83);

        auto g_y_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 84);

        auto g_y_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 85);

        auto g_y_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 86);

        auto g_y_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 87);

        auto g_y_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 88);

        auto g_y_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_y_0_yz_xxx, g_y_0_yz_xxxz, g_y_0_yz_xxy, g_y_0_yz_xxyz, g_y_0_yz_xxz, g_y_0_yz_xxzz, g_y_0_yz_xyy, g_y_0_yz_xyyz, g_y_0_yz_xyz, g_y_0_yz_xyzz, g_y_0_yz_xzz, g_y_0_yz_xzzz, g_y_0_yz_yyy, g_y_0_yz_yyyz, g_y_0_yz_yyz, g_y_0_yz_yyzz, g_y_0_yz_yzz, g_y_0_yz_yzzz, g_y_0_yz_zzz, g_y_0_yz_zzzz, g_y_0_yzz_xxx, g_y_0_yzz_xxy, g_y_0_yzz_xxz, g_y_0_yzz_xyy, g_y_0_yzz_xyz, g_y_0_yzz_xzz, g_y_0_yzz_yyy, g_y_0_yzz_yyz, g_y_0_yzz_yzz, g_y_0_yzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_xxx[k] = -g_y_0_yz_xxx[k] * cd_z[k] + g_y_0_yz_xxxz[k];

            g_y_0_yzz_xxy[k] = -g_y_0_yz_xxy[k] * cd_z[k] + g_y_0_yz_xxyz[k];

            g_y_0_yzz_xxz[k] = -g_y_0_yz_xxz[k] * cd_z[k] + g_y_0_yz_xxzz[k];

            g_y_0_yzz_xyy[k] = -g_y_0_yz_xyy[k] * cd_z[k] + g_y_0_yz_xyyz[k];

            g_y_0_yzz_xyz[k] = -g_y_0_yz_xyz[k] * cd_z[k] + g_y_0_yz_xyzz[k];

            g_y_0_yzz_xzz[k] = -g_y_0_yz_xzz[k] * cd_z[k] + g_y_0_yz_xzzz[k];

            g_y_0_yzz_yyy[k] = -g_y_0_yz_yyy[k] * cd_z[k] + g_y_0_yz_yyyz[k];

            g_y_0_yzz_yyz[k] = -g_y_0_yz_yyz[k] * cd_z[k] + g_y_0_yz_yyzz[k];

            g_y_0_yzz_yzz[k] = -g_y_0_yz_yzz[k] * cd_z[k] + g_y_0_yz_yzzz[k];

            g_y_0_yzz_zzz[k] = -g_y_0_yz_zzz[k] * cd_z[k] + g_y_0_yz_zzzz[k];
        }

        /// Set up 90-100 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps  + 90);

        auto g_y_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 91);

        auto g_y_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 92);

        auto g_y_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 93);

        auto g_y_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 94);

        auto g_y_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 95);

        auto g_y_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps  + 96);

        auto g_y_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 97);

        auto g_y_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 98);

        auto g_y_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps  + 99);

        #pragma omp simd aligned(cd_z, g_y_0_zz_xxx, g_y_0_zz_xxxz, g_y_0_zz_xxy, g_y_0_zz_xxyz, g_y_0_zz_xxz, g_y_0_zz_xxzz, g_y_0_zz_xyy, g_y_0_zz_xyyz, g_y_0_zz_xyz, g_y_0_zz_xyzz, g_y_0_zz_xzz, g_y_0_zz_xzzz, g_y_0_zz_yyy, g_y_0_zz_yyyz, g_y_0_zz_yyz, g_y_0_zz_yyzz, g_y_0_zz_yzz, g_y_0_zz_yzzz, g_y_0_zz_zzz, g_y_0_zz_zzzz, g_y_0_zzz_xxx, g_y_0_zzz_xxy, g_y_0_zzz_xxz, g_y_0_zzz_xyy, g_y_0_zzz_xyz, g_y_0_zzz_xzz, g_y_0_zzz_yyy, g_y_0_zzz_yyz, g_y_0_zzz_yzz, g_y_0_zzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_xxx[k] = -g_y_0_zz_xxx[k] * cd_z[k] + g_y_0_zz_xxxz[k];

            g_y_0_zzz_xxy[k] = -g_y_0_zz_xxy[k] * cd_z[k] + g_y_0_zz_xxyz[k];

            g_y_0_zzz_xxz[k] = -g_y_0_zz_xxz[k] * cd_z[k] + g_y_0_zz_xxzz[k];

            g_y_0_zzz_xyy[k] = -g_y_0_zz_xyy[k] * cd_z[k] + g_y_0_zz_xyyz[k];

            g_y_0_zzz_xyz[k] = -g_y_0_zz_xyz[k] * cd_z[k] + g_y_0_zz_xyzz[k];

            g_y_0_zzz_xzz[k] = -g_y_0_zz_xzz[k] * cd_z[k] + g_y_0_zz_xzzz[k];

            g_y_0_zzz_yyy[k] = -g_y_0_zz_yyy[k] * cd_z[k] + g_y_0_zz_yyyz[k];

            g_y_0_zzz_yyz[k] = -g_y_0_zz_yyz[k] * cd_z[k] + g_y_0_zz_yyzz[k];

            g_y_0_zzz_yzz[k] = -g_y_0_zz_yzz[k] * cd_z[k] + g_y_0_zz_yzzz[k];

            g_y_0_zzz_zzz[k] = -g_y_0_zz_zzz[k] * cd_z[k] + g_y_0_zz_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 0);

        auto g_z_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 1);

        auto g_z_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 2);

        auto g_z_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 3);

        auto g_z_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 4);

        auto g_z_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 5);

        auto g_z_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 6);

        auto g_z_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 7);

        auto g_z_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 8);

        auto g_z_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_z_0_xx_xxx, g_z_0_xx_xxxx, g_z_0_xx_xxxy, g_z_0_xx_xxxz, g_z_0_xx_xxy, g_z_0_xx_xxyy, g_z_0_xx_xxyz, g_z_0_xx_xxz, g_z_0_xx_xxzz, g_z_0_xx_xyy, g_z_0_xx_xyyy, g_z_0_xx_xyyz, g_z_0_xx_xyz, g_z_0_xx_xyzz, g_z_0_xx_xzz, g_z_0_xx_xzzz, g_z_0_xx_yyy, g_z_0_xx_yyz, g_z_0_xx_yzz, g_z_0_xx_zzz, g_z_0_xxx_xxx, g_z_0_xxx_xxy, g_z_0_xxx_xxz, g_z_0_xxx_xyy, g_z_0_xxx_xyz, g_z_0_xxx_xzz, g_z_0_xxx_yyy, g_z_0_xxx_yyz, g_z_0_xxx_yzz, g_z_0_xxx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_xxx[k] = -g_z_0_xx_xxx[k] * cd_x[k] + g_z_0_xx_xxxx[k];

            g_z_0_xxx_xxy[k] = -g_z_0_xx_xxy[k] * cd_x[k] + g_z_0_xx_xxxy[k];

            g_z_0_xxx_xxz[k] = -g_z_0_xx_xxz[k] * cd_x[k] + g_z_0_xx_xxxz[k];

            g_z_0_xxx_xyy[k] = -g_z_0_xx_xyy[k] * cd_x[k] + g_z_0_xx_xxyy[k];

            g_z_0_xxx_xyz[k] = -g_z_0_xx_xyz[k] * cd_x[k] + g_z_0_xx_xxyz[k];

            g_z_0_xxx_xzz[k] = -g_z_0_xx_xzz[k] * cd_x[k] + g_z_0_xx_xxzz[k];

            g_z_0_xxx_yyy[k] = -g_z_0_xx_yyy[k] * cd_x[k] + g_z_0_xx_xyyy[k];

            g_z_0_xxx_yyz[k] = -g_z_0_xx_yyz[k] * cd_x[k] + g_z_0_xx_xyyz[k];

            g_z_0_xxx_yzz[k] = -g_z_0_xx_yzz[k] * cd_x[k] + g_z_0_xx_xyzz[k];

            g_z_0_xxx_zzz[k] = -g_z_0_xx_zzz[k] * cd_x[k] + g_z_0_xx_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 10);

        auto g_z_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 11);

        auto g_z_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 12);

        auto g_z_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 13);

        auto g_z_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 14);

        auto g_z_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 15);

        auto g_z_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 16);

        auto g_z_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 17);

        auto g_z_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 18);

        auto g_z_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 19);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_xxx, g_z_0_xxy_xxy, g_z_0_xxy_xxz, g_z_0_xxy_xyy, g_z_0_xxy_xyz, g_z_0_xxy_xzz, g_z_0_xxy_yyy, g_z_0_xxy_yyz, g_z_0_xxy_yzz, g_z_0_xxy_zzz, g_z_0_xy_xxx, g_z_0_xy_xxxx, g_z_0_xy_xxxy, g_z_0_xy_xxxz, g_z_0_xy_xxy, g_z_0_xy_xxyy, g_z_0_xy_xxyz, g_z_0_xy_xxz, g_z_0_xy_xxzz, g_z_0_xy_xyy, g_z_0_xy_xyyy, g_z_0_xy_xyyz, g_z_0_xy_xyz, g_z_0_xy_xyzz, g_z_0_xy_xzz, g_z_0_xy_xzzz, g_z_0_xy_yyy, g_z_0_xy_yyz, g_z_0_xy_yzz, g_z_0_xy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_xxx[k] = -g_z_0_xy_xxx[k] * cd_x[k] + g_z_0_xy_xxxx[k];

            g_z_0_xxy_xxy[k] = -g_z_0_xy_xxy[k] * cd_x[k] + g_z_0_xy_xxxy[k];

            g_z_0_xxy_xxz[k] = -g_z_0_xy_xxz[k] * cd_x[k] + g_z_0_xy_xxxz[k];

            g_z_0_xxy_xyy[k] = -g_z_0_xy_xyy[k] * cd_x[k] + g_z_0_xy_xxyy[k];

            g_z_0_xxy_xyz[k] = -g_z_0_xy_xyz[k] * cd_x[k] + g_z_0_xy_xxyz[k];

            g_z_0_xxy_xzz[k] = -g_z_0_xy_xzz[k] * cd_x[k] + g_z_0_xy_xxzz[k];

            g_z_0_xxy_yyy[k] = -g_z_0_xy_yyy[k] * cd_x[k] + g_z_0_xy_xyyy[k];

            g_z_0_xxy_yyz[k] = -g_z_0_xy_yyz[k] * cd_x[k] + g_z_0_xy_xyyz[k];

            g_z_0_xxy_yzz[k] = -g_z_0_xy_yzz[k] * cd_x[k] + g_z_0_xy_xyzz[k];

            g_z_0_xxy_zzz[k] = -g_z_0_xy_zzz[k] * cd_x[k] + g_z_0_xy_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 20);

        auto g_z_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 21);

        auto g_z_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 22);

        auto g_z_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 23);

        auto g_z_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 24);

        auto g_z_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 25);

        auto g_z_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 26);

        auto g_z_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 27);

        auto g_z_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 28);

        auto g_z_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_xxx, g_z_0_xxz_xxy, g_z_0_xxz_xxz, g_z_0_xxz_xyy, g_z_0_xxz_xyz, g_z_0_xxz_xzz, g_z_0_xxz_yyy, g_z_0_xxz_yyz, g_z_0_xxz_yzz, g_z_0_xxz_zzz, g_z_0_xz_xxx, g_z_0_xz_xxxx, g_z_0_xz_xxxy, g_z_0_xz_xxxz, g_z_0_xz_xxy, g_z_0_xz_xxyy, g_z_0_xz_xxyz, g_z_0_xz_xxz, g_z_0_xz_xxzz, g_z_0_xz_xyy, g_z_0_xz_xyyy, g_z_0_xz_xyyz, g_z_0_xz_xyz, g_z_0_xz_xyzz, g_z_0_xz_xzz, g_z_0_xz_xzzz, g_z_0_xz_yyy, g_z_0_xz_yyz, g_z_0_xz_yzz, g_z_0_xz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_xxx[k] = -g_z_0_xz_xxx[k] * cd_x[k] + g_z_0_xz_xxxx[k];

            g_z_0_xxz_xxy[k] = -g_z_0_xz_xxy[k] * cd_x[k] + g_z_0_xz_xxxy[k];

            g_z_0_xxz_xxz[k] = -g_z_0_xz_xxz[k] * cd_x[k] + g_z_0_xz_xxxz[k];

            g_z_0_xxz_xyy[k] = -g_z_0_xz_xyy[k] * cd_x[k] + g_z_0_xz_xxyy[k];

            g_z_0_xxz_xyz[k] = -g_z_0_xz_xyz[k] * cd_x[k] + g_z_0_xz_xxyz[k];

            g_z_0_xxz_xzz[k] = -g_z_0_xz_xzz[k] * cd_x[k] + g_z_0_xz_xxzz[k];

            g_z_0_xxz_yyy[k] = -g_z_0_xz_yyy[k] * cd_x[k] + g_z_0_xz_xyyy[k];

            g_z_0_xxz_yyz[k] = -g_z_0_xz_yyz[k] * cd_x[k] + g_z_0_xz_xyyz[k];

            g_z_0_xxz_yzz[k] = -g_z_0_xz_yzz[k] * cd_x[k] + g_z_0_xz_xyzz[k];

            g_z_0_xxz_zzz[k] = -g_z_0_xz_zzz[k] * cd_x[k] + g_z_0_xz_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 30);

        auto g_z_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 31);

        auto g_z_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 32);

        auto g_z_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 33);

        auto g_z_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 34);

        auto g_z_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 35);

        auto g_z_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 36);

        auto g_z_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 37);

        auto g_z_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 38);

        auto g_z_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 39);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_xxx, g_z_0_xyy_xxy, g_z_0_xyy_xxz, g_z_0_xyy_xyy, g_z_0_xyy_xyz, g_z_0_xyy_xzz, g_z_0_xyy_yyy, g_z_0_xyy_yyz, g_z_0_xyy_yzz, g_z_0_xyy_zzz, g_z_0_yy_xxx, g_z_0_yy_xxxx, g_z_0_yy_xxxy, g_z_0_yy_xxxz, g_z_0_yy_xxy, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxz, g_z_0_yy_xxzz, g_z_0_yy_xyy, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyz, g_z_0_yy_xyzz, g_z_0_yy_xzz, g_z_0_yy_xzzz, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yzz, g_z_0_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_xxx[k] = -g_z_0_yy_xxx[k] * cd_x[k] + g_z_0_yy_xxxx[k];

            g_z_0_xyy_xxy[k] = -g_z_0_yy_xxy[k] * cd_x[k] + g_z_0_yy_xxxy[k];

            g_z_0_xyy_xxz[k] = -g_z_0_yy_xxz[k] * cd_x[k] + g_z_0_yy_xxxz[k];

            g_z_0_xyy_xyy[k] = -g_z_0_yy_xyy[k] * cd_x[k] + g_z_0_yy_xxyy[k];

            g_z_0_xyy_xyz[k] = -g_z_0_yy_xyz[k] * cd_x[k] + g_z_0_yy_xxyz[k];

            g_z_0_xyy_xzz[k] = -g_z_0_yy_xzz[k] * cd_x[k] + g_z_0_yy_xxzz[k];

            g_z_0_xyy_yyy[k] = -g_z_0_yy_yyy[k] * cd_x[k] + g_z_0_yy_xyyy[k];

            g_z_0_xyy_yyz[k] = -g_z_0_yy_yyz[k] * cd_x[k] + g_z_0_yy_xyyz[k];

            g_z_0_xyy_yzz[k] = -g_z_0_yy_yzz[k] * cd_x[k] + g_z_0_yy_xyzz[k];

            g_z_0_xyy_zzz[k] = -g_z_0_yy_zzz[k] * cd_x[k] + g_z_0_yy_xzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 40);

        auto g_z_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 41);

        auto g_z_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 42);

        auto g_z_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 43);

        auto g_z_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 44);

        auto g_z_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 45);

        auto g_z_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 46);

        auto g_z_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 47);

        auto g_z_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 48);

        auto g_z_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 49);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_xxx, g_z_0_xyz_xxy, g_z_0_xyz_xxz, g_z_0_xyz_xyy, g_z_0_xyz_xyz, g_z_0_xyz_xzz, g_z_0_xyz_yyy, g_z_0_xyz_yyz, g_z_0_xyz_yzz, g_z_0_xyz_zzz, g_z_0_yz_xxx, g_z_0_yz_xxxx, g_z_0_yz_xxxy, g_z_0_yz_xxxz, g_z_0_yz_xxy, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxz, g_z_0_yz_xxzz, g_z_0_yz_xyy, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyz, g_z_0_yz_xyzz, g_z_0_yz_xzz, g_z_0_yz_xzzz, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yzz, g_z_0_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_xxx[k] = -g_z_0_yz_xxx[k] * cd_x[k] + g_z_0_yz_xxxx[k];

            g_z_0_xyz_xxy[k] = -g_z_0_yz_xxy[k] * cd_x[k] + g_z_0_yz_xxxy[k];

            g_z_0_xyz_xxz[k] = -g_z_0_yz_xxz[k] * cd_x[k] + g_z_0_yz_xxxz[k];

            g_z_0_xyz_xyy[k] = -g_z_0_yz_xyy[k] * cd_x[k] + g_z_0_yz_xxyy[k];

            g_z_0_xyz_xyz[k] = -g_z_0_yz_xyz[k] * cd_x[k] + g_z_0_yz_xxyz[k];

            g_z_0_xyz_xzz[k] = -g_z_0_yz_xzz[k] * cd_x[k] + g_z_0_yz_xxzz[k];

            g_z_0_xyz_yyy[k] = -g_z_0_yz_yyy[k] * cd_x[k] + g_z_0_yz_xyyy[k];

            g_z_0_xyz_yyz[k] = -g_z_0_yz_yyz[k] * cd_x[k] + g_z_0_yz_xyyz[k];

            g_z_0_xyz_yzz[k] = -g_z_0_yz_yzz[k] * cd_x[k] + g_z_0_yz_xyzz[k];

            g_z_0_xyz_zzz[k] = -g_z_0_yz_zzz[k] * cd_x[k] + g_z_0_yz_xzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 50);

        auto g_z_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 51);

        auto g_z_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 52);

        auto g_z_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 53);

        auto g_z_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 54);

        auto g_z_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 55);

        auto g_z_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 56);

        auto g_z_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 57);

        auto g_z_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 58);

        auto g_z_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_xxx, g_z_0_xzz_xxy, g_z_0_xzz_xxz, g_z_0_xzz_xyy, g_z_0_xzz_xyz, g_z_0_xzz_xzz, g_z_0_xzz_yyy, g_z_0_xzz_yyz, g_z_0_xzz_yzz, g_z_0_xzz_zzz, g_z_0_zz_xxx, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxy, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxz, g_z_0_zz_xxzz, g_z_0_zz_xyy, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyz, g_z_0_zz_xyzz, g_z_0_zz_xzz, g_z_0_zz_xzzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_xxx[k] = -g_z_0_zz_xxx[k] * cd_x[k] + g_z_0_zz_xxxx[k];

            g_z_0_xzz_xxy[k] = -g_z_0_zz_xxy[k] * cd_x[k] + g_z_0_zz_xxxy[k];

            g_z_0_xzz_xxz[k] = -g_z_0_zz_xxz[k] * cd_x[k] + g_z_0_zz_xxxz[k];

            g_z_0_xzz_xyy[k] = -g_z_0_zz_xyy[k] * cd_x[k] + g_z_0_zz_xxyy[k];

            g_z_0_xzz_xyz[k] = -g_z_0_zz_xyz[k] * cd_x[k] + g_z_0_zz_xxyz[k];

            g_z_0_xzz_xzz[k] = -g_z_0_zz_xzz[k] * cd_x[k] + g_z_0_zz_xxzz[k];

            g_z_0_xzz_yyy[k] = -g_z_0_zz_yyy[k] * cd_x[k] + g_z_0_zz_xyyy[k];

            g_z_0_xzz_yyz[k] = -g_z_0_zz_yyz[k] * cd_x[k] + g_z_0_zz_xyyz[k];

            g_z_0_xzz_yzz[k] = -g_z_0_zz_yzz[k] * cd_x[k] + g_z_0_zz_xyzz[k];

            g_z_0_xzz_zzz[k] = -g_z_0_zz_zzz[k] * cd_x[k] + g_z_0_zz_xzzz[k];
        }

        /// Set up 60-70 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 60);

        auto g_z_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 61);

        auto g_z_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 62);

        auto g_z_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 63);

        auto g_z_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 64);

        auto g_z_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 65);

        auto g_z_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 66);

        auto g_z_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 67);

        auto g_z_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 68);

        auto g_z_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 69);

        #pragma omp simd aligned(cd_y, g_z_0_yy_xxx, g_z_0_yy_xxxy, g_z_0_yy_xxy, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxz, g_z_0_yy_xyy, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyz, g_z_0_yy_xyzz, g_z_0_yy_xzz, g_z_0_yy_yyy, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyz, g_z_0_yy_yyzz, g_z_0_yy_yzz, g_z_0_yy_yzzz, g_z_0_yy_zzz, g_z_0_yyy_xxx, g_z_0_yyy_xxy, g_z_0_yyy_xxz, g_z_0_yyy_xyy, g_z_0_yyy_xyz, g_z_0_yyy_xzz, g_z_0_yyy_yyy, g_z_0_yyy_yyz, g_z_0_yyy_yzz, g_z_0_yyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_xxx[k] = -g_z_0_yy_xxx[k] * cd_y[k] + g_z_0_yy_xxxy[k];

            g_z_0_yyy_xxy[k] = -g_z_0_yy_xxy[k] * cd_y[k] + g_z_0_yy_xxyy[k];

            g_z_0_yyy_xxz[k] = -g_z_0_yy_xxz[k] * cd_y[k] + g_z_0_yy_xxyz[k];

            g_z_0_yyy_xyy[k] = -g_z_0_yy_xyy[k] * cd_y[k] + g_z_0_yy_xyyy[k];

            g_z_0_yyy_xyz[k] = -g_z_0_yy_xyz[k] * cd_y[k] + g_z_0_yy_xyyz[k];

            g_z_0_yyy_xzz[k] = -g_z_0_yy_xzz[k] * cd_y[k] + g_z_0_yy_xyzz[k];

            g_z_0_yyy_yyy[k] = -g_z_0_yy_yyy[k] * cd_y[k] + g_z_0_yy_yyyy[k];

            g_z_0_yyy_yyz[k] = -g_z_0_yy_yyz[k] * cd_y[k] + g_z_0_yy_yyyz[k];

            g_z_0_yyy_yzz[k] = -g_z_0_yy_yzz[k] * cd_y[k] + g_z_0_yy_yyzz[k];

            g_z_0_yyy_zzz[k] = -g_z_0_yy_zzz[k] * cd_y[k] + g_z_0_yy_yzzz[k];
        }

        /// Set up 70-80 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 70);

        auto g_z_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 71);

        auto g_z_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 72);

        auto g_z_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 73);

        auto g_z_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 74);

        auto g_z_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 75);

        auto g_z_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 76);

        auto g_z_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 77);

        auto g_z_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 78);

        auto g_z_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 79);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_xxx, g_z_0_yyz_xxy, g_z_0_yyz_xxz, g_z_0_yyz_xyy, g_z_0_yyz_xyz, g_z_0_yyz_xzz, g_z_0_yyz_yyy, g_z_0_yyz_yyz, g_z_0_yyz_yzz, g_z_0_yyz_zzz, g_z_0_yz_xxx, g_z_0_yz_xxxy, g_z_0_yz_xxy, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxz, g_z_0_yz_xyy, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyz, g_z_0_yz_xyzz, g_z_0_yz_xzz, g_z_0_yz_yyy, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyz, g_z_0_yz_yyzz, g_z_0_yz_yzz, g_z_0_yz_yzzz, g_z_0_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_xxx[k] = -g_z_0_yz_xxx[k] * cd_y[k] + g_z_0_yz_xxxy[k];

            g_z_0_yyz_xxy[k] = -g_z_0_yz_xxy[k] * cd_y[k] + g_z_0_yz_xxyy[k];

            g_z_0_yyz_xxz[k] = -g_z_0_yz_xxz[k] * cd_y[k] + g_z_0_yz_xxyz[k];

            g_z_0_yyz_xyy[k] = -g_z_0_yz_xyy[k] * cd_y[k] + g_z_0_yz_xyyy[k];

            g_z_0_yyz_xyz[k] = -g_z_0_yz_xyz[k] * cd_y[k] + g_z_0_yz_xyyz[k];

            g_z_0_yyz_xzz[k] = -g_z_0_yz_xzz[k] * cd_y[k] + g_z_0_yz_xyzz[k];

            g_z_0_yyz_yyy[k] = -g_z_0_yz_yyy[k] * cd_y[k] + g_z_0_yz_yyyy[k];

            g_z_0_yyz_yyz[k] = -g_z_0_yz_yyz[k] * cd_y[k] + g_z_0_yz_yyyz[k];

            g_z_0_yyz_yzz[k] = -g_z_0_yz_yzz[k] * cd_y[k] + g_z_0_yz_yyzz[k];

            g_z_0_yyz_zzz[k] = -g_z_0_yz_zzz[k] * cd_y[k] + g_z_0_yz_yzzz[k];
        }

        /// Set up 80-90 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 80);

        auto g_z_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 81);

        auto g_z_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 82);

        auto g_z_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 83);

        auto g_z_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 84);

        auto g_z_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 85);

        auto g_z_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 86);

        auto g_z_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 87);

        auto g_z_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 88);

        auto g_z_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 89);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_xxx, g_z_0_yzz_xxy, g_z_0_yzz_xxz, g_z_0_yzz_xyy, g_z_0_yzz_xyz, g_z_0_yzz_xzz, g_z_0_yzz_yyy, g_z_0_yzz_yyz, g_z_0_yzz_yzz, g_z_0_yzz_zzz, g_z_0_zz_xxx, g_z_0_zz_xxxy, g_z_0_zz_xxy, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyz, g_z_0_zz_xyzz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyz, g_z_0_zz_yyzz, g_z_0_zz_yzz, g_z_0_zz_yzzz, g_z_0_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_xxx[k] = -g_z_0_zz_xxx[k] * cd_y[k] + g_z_0_zz_xxxy[k];

            g_z_0_yzz_xxy[k] = -g_z_0_zz_xxy[k] * cd_y[k] + g_z_0_zz_xxyy[k];

            g_z_0_yzz_xxz[k] = -g_z_0_zz_xxz[k] * cd_y[k] + g_z_0_zz_xxyz[k];

            g_z_0_yzz_xyy[k] = -g_z_0_zz_xyy[k] * cd_y[k] + g_z_0_zz_xyyy[k];

            g_z_0_yzz_xyz[k] = -g_z_0_zz_xyz[k] * cd_y[k] + g_z_0_zz_xyyz[k];

            g_z_0_yzz_xzz[k] = -g_z_0_zz_xzz[k] * cd_y[k] + g_z_0_zz_xyzz[k];

            g_z_0_yzz_yyy[k] = -g_z_0_zz_yyy[k] * cd_y[k] + g_z_0_zz_yyyy[k];

            g_z_0_yzz_yyz[k] = -g_z_0_zz_yyz[k] * cd_y[k] + g_z_0_zz_yyyz[k];

            g_z_0_yzz_yzz[k] = -g_z_0_zz_yzz[k] * cd_y[k] + g_z_0_zz_yyzz[k];

            g_z_0_yzz_zzz[k] = -g_z_0_zz_zzz[k] * cd_y[k] + g_z_0_zz_yzzz[k];
        }

        /// Set up 90-100 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps  + 90);

        auto g_z_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 91);

        auto g_z_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 92);

        auto g_z_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 93);

        auto g_z_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 94);

        auto g_z_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 95);

        auto g_z_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps  + 96);

        auto g_z_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 97);

        auto g_z_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 98);

        auto g_z_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps  + 99);

        #pragma omp simd aligned(cd_z, g_z_0_zz_xxx, g_z_0_zz_xxxz, g_z_0_zz_xxy, g_z_0_zz_xxyz, g_z_0_zz_xxz, g_z_0_zz_xxzz, g_z_0_zz_xyy, g_z_0_zz_xyyz, g_z_0_zz_xyz, g_z_0_zz_xyzz, g_z_0_zz_xzz, g_z_0_zz_xzzz, g_z_0_zz_yyy, g_z_0_zz_yyyz, g_z_0_zz_yyz, g_z_0_zz_yyzz, g_z_0_zz_yzz, g_z_0_zz_yzzz, g_z_0_zz_zzz, g_z_0_zz_zzzz, g_z_0_zzz_xxx, g_z_0_zzz_xxy, g_z_0_zzz_xxz, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xzz, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yzz, g_z_0_zzz_zzz, g_zz_xxx, g_zz_xxy, g_zz_xxz, g_zz_xyy, g_zz_xyz, g_zz_xzz, g_zz_yyy, g_zz_yyz, g_zz_yzz, g_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_xxx[k] = -g_zz_xxx[k] - g_z_0_zz_xxx[k] * cd_z[k] + g_z_0_zz_xxxz[k];

            g_z_0_zzz_xxy[k] = -g_zz_xxy[k] - g_z_0_zz_xxy[k] * cd_z[k] + g_z_0_zz_xxyz[k];

            g_z_0_zzz_xxz[k] = -g_zz_xxz[k] - g_z_0_zz_xxz[k] * cd_z[k] + g_z_0_zz_xxzz[k];

            g_z_0_zzz_xyy[k] = -g_zz_xyy[k] - g_z_0_zz_xyy[k] * cd_z[k] + g_z_0_zz_xyyz[k];

            g_z_0_zzz_xyz[k] = -g_zz_xyz[k] - g_z_0_zz_xyz[k] * cd_z[k] + g_z_0_zz_xyzz[k];

            g_z_0_zzz_xzz[k] = -g_zz_xzz[k] - g_z_0_zz_xzz[k] * cd_z[k] + g_z_0_zz_xzzz[k];

            g_z_0_zzz_yyy[k] = -g_zz_yyy[k] - g_z_0_zz_yyy[k] * cd_z[k] + g_z_0_zz_yyyz[k];

            g_z_0_zzz_yyz[k] = -g_zz_yyz[k] - g_z_0_zz_yyz[k] * cd_z[k] + g_z_0_zz_yyzz[k];

            g_z_0_zzz_yzz[k] = -g_zz_yzz[k] - g_z_0_zz_yzz[k] * cd_z[k] + g_z_0_zz_yzzz[k];

            g_z_0_zzz_zzz[k] = -g_zz_zzz[k] - g_z_0_zz_zzz[k] * cd_z[k] + g_z_0_zz_zzzz[k];
        }
    }
}

} // t3ceri namespace

