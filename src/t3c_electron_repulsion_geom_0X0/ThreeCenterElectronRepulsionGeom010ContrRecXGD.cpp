#include "ThreeCenterElectronRepulsionGeom010ContrRecXGD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xgd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xgd,
                                        const size_t idx_xfd,
                                        const size_t idx_geom_10_xfd,
                                        const size_t idx_geom_10_xff,
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
        /// Set up components of auxilary buffer : SFD

        const auto fd_off = idx_xfd + i * 60;

        auto g_xxx_xx = cbuffer.data(fd_off + 0);

        auto g_xxx_xy = cbuffer.data(fd_off + 1);

        auto g_xxx_xz = cbuffer.data(fd_off + 2);

        auto g_xxx_yy = cbuffer.data(fd_off + 3);

        auto g_xxx_yz = cbuffer.data(fd_off + 4);

        auto g_xxx_zz = cbuffer.data(fd_off + 5);

        auto g_xxy_xx = cbuffer.data(fd_off + 6);

        auto g_xxy_xy = cbuffer.data(fd_off + 7);

        auto g_xxy_xz = cbuffer.data(fd_off + 8);

        auto g_xxy_yy = cbuffer.data(fd_off + 9);

        auto g_xxy_yz = cbuffer.data(fd_off + 10);

        auto g_xxy_zz = cbuffer.data(fd_off + 11);

        auto g_xxz_xx = cbuffer.data(fd_off + 12);

        auto g_xxz_xy = cbuffer.data(fd_off + 13);

        auto g_xxz_xz = cbuffer.data(fd_off + 14);

        auto g_xxz_yy = cbuffer.data(fd_off + 15);

        auto g_xxz_yz = cbuffer.data(fd_off + 16);

        auto g_xxz_zz = cbuffer.data(fd_off + 17);

        auto g_xyy_xx = cbuffer.data(fd_off + 18);

        auto g_xyy_xy = cbuffer.data(fd_off + 19);

        auto g_xyy_xz = cbuffer.data(fd_off + 20);

        auto g_xyy_yy = cbuffer.data(fd_off + 21);

        auto g_xyy_yz = cbuffer.data(fd_off + 22);

        auto g_xyy_zz = cbuffer.data(fd_off + 23);

        auto g_xyz_xx = cbuffer.data(fd_off + 24);

        auto g_xyz_xy = cbuffer.data(fd_off + 25);

        auto g_xyz_xz = cbuffer.data(fd_off + 26);

        auto g_xyz_yy = cbuffer.data(fd_off + 27);

        auto g_xyz_yz = cbuffer.data(fd_off + 28);

        auto g_xyz_zz = cbuffer.data(fd_off + 29);

        auto g_xzz_xx = cbuffer.data(fd_off + 30);

        auto g_xzz_xy = cbuffer.data(fd_off + 31);

        auto g_xzz_xz = cbuffer.data(fd_off + 32);

        auto g_xzz_yy = cbuffer.data(fd_off + 33);

        auto g_xzz_yz = cbuffer.data(fd_off + 34);

        auto g_xzz_zz = cbuffer.data(fd_off + 35);

        auto g_yyy_xx = cbuffer.data(fd_off + 36);

        auto g_yyy_xy = cbuffer.data(fd_off + 37);

        auto g_yyy_xz = cbuffer.data(fd_off + 38);

        auto g_yyy_yy = cbuffer.data(fd_off + 39);

        auto g_yyy_yz = cbuffer.data(fd_off + 40);

        auto g_yyy_zz = cbuffer.data(fd_off + 41);

        auto g_yyz_xx = cbuffer.data(fd_off + 42);

        auto g_yyz_xy = cbuffer.data(fd_off + 43);

        auto g_yyz_xz = cbuffer.data(fd_off + 44);

        auto g_yyz_yy = cbuffer.data(fd_off + 45);

        auto g_yyz_yz = cbuffer.data(fd_off + 46);

        auto g_yyz_zz = cbuffer.data(fd_off + 47);

        auto g_yzz_xx = cbuffer.data(fd_off + 48);

        auto g_yzz_xy = cbuffer.data(fd_off + 49);

        auto g_yzz_xz = cbuffer.data(fd_off + 50);

        auto g_yzz_yy = cbuffer.data(fd_off + 51);

        auto g_yzz_yz = cbuffer.data(fd_off + 52);

        auto g_yzz_zz = cbuffer.data(fd_off + 53);

        auto g_zzz_xx = cbuffer.data(fd_off + 54);

        auto g_zzz_xy = cbuffer.data(fd_off + 55);

        auto g_zzz_xz = cbuffer.data(fd_off + 56);

        auto g_zzz_yy = cbuffer.data(fd_off + 57);

        auto g_zzz_yz = cbuffer.data(fd_off + 58);

        auto g_zzz_zz = cbuffer.data(fd_off + 59);

        /// Set up components of auxilary buffer : SFD

        const auto fd_geom_10_off = idx_geom_10_xfd + i * 60;

        auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xxy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xxy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xxy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xxz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xyy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xyy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xyy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xyz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xyz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xyz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 30);

        auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 35);

        auto g_x_0_yyy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 36);

        auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 37);

        auto g_x_0_yyy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 38);

        auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 39);

        auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 40);

        auto g_x_0_yyy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 41);

        auto g_x_0_yyz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 42);

        auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 43);

        auto g_x_0_yyz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 44);

        auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 45);

        auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 46);

        auto g_x_0_yyz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 47);

        auto g_x_0_yzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 48);

        auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 49);

        auto g_x_0_yzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 50);

        auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 51);

        auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 52);

        auto g_x_0_yzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 53);

        auto g_x_0_zzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps + 54);

        auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps + 55);

        auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps + 56);

        auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps + 57);

        auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps + 58);

        auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps + 59);

        auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 0);

        auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 1);

        auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 2);

        auto g_y_0_xxx_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 3);

        auto g_y_0_xxx_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 4);

        auto g_y_0_xxx_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 5);

        auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 6);

        auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 7);

        auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 8);

        auto g_y_0_xxy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 9);

        auto g_y_0_xxy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 10);

        auto g_y_0_xxy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 11);

        auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 12);

        auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 13);

        auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 14);

        auto g_y_0_xxz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 15);

        auto g_y_0_xxz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 16);

        auto g_y_0_xxz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 17);

        auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 18);

        auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 19);

        auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 20);

        auto g_y_0_xyy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 21);

        auto g_y_0_xyy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 22);

        auto g_y_0_xyy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 23);

        auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 24);

        auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 25);

        auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 26);

        auto g_y_0_xyz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 27);

        auto g_y_0_xyz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 28);

        auto g_y_0_xyz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 29);

        auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 30);

        auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 31);

        auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 32);

        auto g_y_0_xzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 33);

        auto g_y_0_xzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 34);

        auto g_y_0_xzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 35);

        auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 36);

        auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 37);

        auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 38);

        auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 39);

        auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 40);

        auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 41);

        auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 42);

        auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 43);

        auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 44);

        auto g_y_0_yyz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 45);

        auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 46);

        auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 47);

        auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 48);

        auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 49);

        auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 50);

        auto g_y_0_yzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 51);

        auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 52);

        auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 53);

        auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps + 54);

        auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps + 55);

        auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps + 56);

        auto g_y_0_zzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps + 57);

        auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps + 58);

        auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps + 59);

        auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 0);

        auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 1);

        auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 2);

        auto g_z_0_xxx_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 3);

        auto g_z_0_xxx_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 4);

        auto g_z_0_xxx_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 5);

        auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 6);

        auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 7);

        auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 8);

        auto g_z_0_xxy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 9);

        auto g_z_0_xxy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 10);

        auto g_z_0_xxy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 11);

        auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 12);

        auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 13);

        auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 14);

        auto g_z_0_xxz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 15);

        auto g_z_0_xxz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 16);

        auto g_z_0_xxz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 17);

        auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 18);

        auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 19);

        auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 20);

        auto g_z_0_xyy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 21);

        auto g_z_0_xyy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 22);

        auto g_z_0_xyy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 23);

        auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 24);

        auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 25);

        auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 26);

        auto g_z_0_xyz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 27);

        auto g_z_0_xyz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 28);

        auto g_z_0_xyz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 29);

        auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 30);

        auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 31);

        auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 32);

        auto g_z_0_xzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 33);

        auto g_z_0_xzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 34);

        auto g_z_0_xzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 35);

        auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 36);

        auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 37);

        auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 38);

        auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 39);

        auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 40);

        auto g_z_0_yyy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 41);

        auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 42);

        auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 43);

        auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 44);

        auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 45);

        auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 46);

        auto g_z_0_yyz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 47);

        auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 48);

        auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 49);

        auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 50);

        auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 51);

        auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 52);

        auto g_z_0_yzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 53);

        auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps + 54);

        auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps + 55);

        auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps + 56);

        auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps + 57);

        auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps + 58);

        auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps + 59);

        /// Set up components of auxilary buffer : SFF

        const auto ff_geom_10_off = idx_geom_10_xff + i * 100;

        auto g_x_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 30);

        auto g_x_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 35);

        auto g_x_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 41);

        auto g_x_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 42);

        auto g_x_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 48);

        auto g_x_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 55);

        auto g_x_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 56);

        auto g_x_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 59);

        auto g_x_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 60);

        auto g_x_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 61);

        auto g_x_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 62);

        auto g_x_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 63);

        auto g_x_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 64);

        auto g_x_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 65);

        auto g_x_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 66);

        auto g_x_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 67);

        auto g_x_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 68);

        auto g_x_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 69);

        auto g_x_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 70);

        auto g_x_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 71);

        auto g_x_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 72);

        auto g_x_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 73);

        auto g_x_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 74);

        auto g_x_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 75);

        auto g_x_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 76);

        auto g_x_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 77);

        auto g_x_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 78);

        auto g_x_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 79);

        auto g_x_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 80);

        auto g_x_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 81);

        auto g_x_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 82);

        auto g_x_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 83);

        auto g_x_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 84);

        auto g_x_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 85);

        auto g_x_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 86);

        auto g_x_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 87);

        auto g_x_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 88);

        auto g_x_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 89);

        auto g_x_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps + 90);

        auto g_x_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps + 91);

        auto g_x_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps + 92);

        auto g_x_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 93);

        auto g_x_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 94);

        auto g_x_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 95);

        auto g_x_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps + 96);

        auto g_x_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps + 97);

        auto g_x_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 98);

        auto g_x_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps + 99);

        auto g_y_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 0);

        auto g_y_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 1);

        auto g_y_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 2);

        auto g_y_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 3);

        auto g_y_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 4);

        auto g_y_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 5);

        auto g_y_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 6);

        auto g_y_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 7);

        auto g_y_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 8);

        auto g_y_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 9);

        auto g_y_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 10);

        auto g_y_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 11);

        auto g_y_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 12);

        auto g_y_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 13);

        auto g_y_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 14);

        auto g_y_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 15);

        auto g_y_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 16);

        auto g_y_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 17);

        auto g_y_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 18);

        auto g_y_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 19);

        auto g_y_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 20);

        auto g_y_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 21);

        auto g_y_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 22);

        auto g_y_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 23);

        auto g_y_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 24);

        auto g_y_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 25);

        auto g_y_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 26);

        auto g_y_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 27);

        auto g_y_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 28);

        auto g_y_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 29);

        auto g_y_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 30);

        auto g_y_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 31);

        auto g_y_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 32);

        auto g_y_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 33);

        auto g_y_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 34);

        auto g_y_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 35);

        auto g_y_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 36);

        auto g_y_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 37);

        auto g_y_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 38);

        auto g_y_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 39);

        auto g_y_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 40);

        auto g_y_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 41);

        auto g_y_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 42);

        auto g_y_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 43);

        auto g_y_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 44);

        auto g_y_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 45);

        auto g_y_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 46);

        auto g_y_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 47);

        auto g_y_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 48);

        auto g_y_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 49);

        auto g_y_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 50);

        auto g_y_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 51);

        auto g_y_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 52);

        auto g_y_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 53);

        auto g_y_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 54);

        auto g_y_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 55);

        auto g_y_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 56);

        auto g_y_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 57);

        auto g_y_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 58);

        auto g_y_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 59);

        auto g_y_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 60);

        auto g_y_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 61);

        auto g_y_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 62);

        auto g_y_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 63);

        auto g_y_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 64);

        auto g_y_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 65);

        auto g_y_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 66);

        auto g_y_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 67);

        auto g_y_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 68);

        auto g_y_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 69);

        auto g_y_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 70);

        auto g_y_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 71);

        auto g_y_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 72);

        auto g_y_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 73);

        auto g_y_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 74);

        auto g_y_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 75);

        auto g_y_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 76);

        auto g_y_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 77);

        auto g_y_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 78);

        auto g_y_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 79);

        auto g_y_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 80);

        auto g_y_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 81);

        auto g_y_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 82);

        auto g_y_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 83);

        auto g_y_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 84);

        auto g_y_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 85);

        auto g_y_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 86);

        auto g_y_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 87);

        auto g_y_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 88);

        auto g_y_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 89);

        auto g_y_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps + 90);

        auto g_y_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps + 91);

        auto g_y_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps + 92);

        auto g_y_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 93);

        auto g_y_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 94);

        auto g_y_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 95);

        auto g_y_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps + 96);

        auto g_y_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps + 97);

        auto g_y_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 98);

        auto g_y_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps + 99);

        auto g_z_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 0);

        auto g_z_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 1);

        auto g_z_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 2);

        auto g_z_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 3);

        auto g_z_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 4);

        auto g_z_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 5);

        auto g_z_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 6);

        auto g_z_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 7);

        auto g_z_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 8);

        auto g_z_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 9);

        auto g_z_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 10);

        auto g_z_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 11);

        auto g_z_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 12);

        auto g_z_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 13);

        auto g_z_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 14);

        auto g_z_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 15);

        auto g_z_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 16);

        auto g_z_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 17);

        auto g_z_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 18);

        auto g_z_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 19);

        auto g_z_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 20);

        auto g_z_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 21);

        auto g_z_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 22);

        auto g_z_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 23);

        auto g_z_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 24);

        auto g_z_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 25);

        auto g_z_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 26);

        auto g_z_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 27);

        auto g_z_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 28);

        auto g_z_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 29);

        auto g_z_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 30);

        auto g_z_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 31);

        auto g_z_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 32);

        auto g_z_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 33);

        auto g_z_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 34);

        auto g_z_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 35);

        auto g_z_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 36);

        auto g_z_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 37);

        auto g_z_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 38);

        auto g_z_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 39);

        auto g_z_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 40);

        auto g_z_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 41);

        auto g_z_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 42);

        auto g_z_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 43);

        auto g_z_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 44);

        auto g_z_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 45);

        auto g_z_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 46);

        auto g_z_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 47);

        auto g_z_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 48);

        auto g_z_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 49);

        auto g_z_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 50);

        auto g_z_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 51);

        auto g_z_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 52);

        auto g_z_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 53);

        auto g_z_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 54);

        auto g_z_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 55);

        auto g_z_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 56);

        auto g_z_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 57);

        auto g_z_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 58);

        auto g_z_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 59);

        auto g_z_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 60);

        auto g_z_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 61);

        auto g_z_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 62);

        auto g_z_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 63);

        auto g_z_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 64);

        auto g_z_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 65);

        auto g_z_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 66);

        auto g_z_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 67);

        auto g_z_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 68);

        auto g_z_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 69);

        auto g_z_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 70);

        auto g_z_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 71);

        auto g_z_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 72);

        auto g_z_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 73);

        auto g_z_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 74);

        auto g_z_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 75);

        auto g_z_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 76);

        auto g_z_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 77);

        auto g_z_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 78);

        auto g_z_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 79);

        auto g_z_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 80);

        auto g_z_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 81);

        auto g_z_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 82);

        auto g_z_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 83);

        auto g_z_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 84);

        auto g_z_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 85);

        auto g_z_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 86);

        auto g_z_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 87);

        auto g_z_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 88);

        auto g_z_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 89);

        auto g_z_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps + 90);

        auto g_z_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps + 91);

        auto g_z_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps + 92);

        auto g_z_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 93);

        auto g_z_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 94);

        auto g_z_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 95);

        auto g_z_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps + 96);

        auto g_z_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps + 97);

        auto g_z_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 98);

        auto g_z_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps + 99);

        /// set up bra offset for contr_buffer_xxgd

        const auto gd_geom_10_off = idx_geom_10_xgd + i * 90;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_x_0_xxx_xx, g_x_0_xxx_xxx, g_x_0_xxx_xxy, g_x_0_xxx_xxz, g_x_0_xxx_xy, g_x_0_xxx_xyy, g_x_0_xxx_xyz, g_x_0_xxx_xz, g_x_0_xxx_xzz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, g_x_0_xxxx_xx, g_x_0_xxxx_xy, g_x_0_xxxx_xz, g_x_0_xxxx_yy, g_x_0_xxxx_yz, g_x_0_xxxx_zz, g_xxx_xx, g_xxx_xy, g_xxx_xz, g_xxx_yy, g_xxx_yz, g_xxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxx_xx[k] = -g_xxx_xx[k] - g_x_0_xxx_xx[k] * cd_x[k] + g_x_0_xxx_xxx[k];

            g_x_0_xxxx_xy[k] = -g_xxx_xy[k] - g_x_0_xxx_xy[k] * cd_x[k] + g_x_0_xxx_xxy[k];

            g_x_0_xxxx_xz[k] = -g_xxx_xz[k] - g_x_0_xxx_xz[k] * cd_x[k] + g_x_0_xxx_xxz[k];

            g_x_0_xxxx_yy[k] = -g_xxx_yy[k] - g_x_0_xxx_yy[k] * cd_x[k] + g_x_0_xxx_xyy[k];

            g_x_0_xxxx_yz[k] = -g_xxx_yz[k] - g_x_0_xxx_yz[k] * cd_x[k] + g_x_0_xxx_xyz[k];

            g_x_0_xxxx_zz[k] = -g_xxx_zz[k] - g_x_0_xxx_zz[k] * cd_x[k] + g_x_0_xxx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_xxx_xx, g_x_0_xxx_xxy, g_x_0_xxx_xy, g_x_0_xxx_xyy, g_x_0_xxx_xyz, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yyy, g_x_0_xxx_yyz, g_x_0_xxx_yz, g_x_0_xxx_yzz, g_x_0_xxx_zz, g_x_0_xxxy_xx, g_x_0_xxxy_xy, g_x_0_xxxy_xz, g_x_0_xxxy_yy, g_x_0_xxxy_yz, g_x_0_xxxy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxy_xx[k] = -g_x_0_xxx_xx[k] * cd_y[k] + g_x_0_xxx_xxy[k];

            g_x_0_xxxy_xy[k] = -g_x_0_xxx_xy[k] * cd_y[k] + g_x_0_xxx_xyy[k];

            g_x_0_xxxy_xz[k] = -g_x_0_xxx_xz[k] * cd_y[k] + g_x_0_xxx_xyz[k];

            g_x_0_xxxy_yy[k] = -g_x_0_xxx_yy[k] * cd_y[k] + g_x_0_xxx_yyy[k];

            g_x_0_xxxy_yz[k] = -g_x_0_xxx_yz[k] * cd_y[k] + g_x_0_xxx_yyz[k];

            g_x_0_xxxy_zz[k] = -g_x_0_xxx_zz[k] * cd_y[k] + g_x_0_xxx_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_xxx_xx, g_x_0_xxx_xxz, g_x_0_xxx_xy, g_x_0_xxx_xyz, g_x_0_xxx_xz, g_x_0_xxx_xzz, g_x_0_xxx_yy, g_x_0_xxx_yyz, g_x_0_xxx_yz, g_x_0_xxx_yzz, g_x_0_xxx_zz, g_x_0_xxx_zzz, g_x_0_xxxz_xx, g_x_0_xxxz_xy, g_x_0_xxxz_xz, g_x_0_xxxz_yy, g_x_0_xxxz_yz, g_x_0_xxxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxz_xx[k] = -g_x_0_xxx_xx[k] * cd_z[k] + g_x_0_xxx_xxz[k];

            g_x_0_xxxz_xy[k] = -g_x_0_xxx_xy[k] * cd_z[k] + g_x_0_xxx_xyz[k];

            g_x_0_xxxz_xz[k] = -g_x_0_xxx_xz[k] * cd_z[k] + g_x_0_xxx_xzz[k];

            g_x_0_xxxz_yy[k] = -g_x_0_xxx_yy[k] * cd_z[k] + g_x_0_xxx_yyz[k];

            g_x_0_xxxz_yz[k] = -g_x_0_xxx_yz[k] * cd_z[k] + g_x_0_xxx_yzz[k];

            g_x_0_xxxz_zz[k] = -g_x_0_xxx_zz[k] * cd_z[k] + g_x_0_xxx_zzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_x_0_xxy_xx, g_x_0_xxy_xxy, g_x_0_xxy_xy, g_x_0_xxy_xyy, g_x_0_xxy_xyz, g_x_0_xxy_xz, g_x_0_xxy_yy, g_x_0_xxy_yyy, g_x_0_xxy_yyz, g_x_0_xxy_yz, g_x_0_xxy_yzz, g_x_0_xxy_zz, g_x_0_xxyy_xx, g_x_0_xxyy_xy, g_x_0_xxyy_xz, g_x_0_xxyy_yy, g_x_0_xxyy_yz, g_x_0_xxyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyy_xx[k] = -g_x_0_xxy_xx[k] * cd_y[k] + g_x_0_xxy_xxy[k];

            g_x_0_xxyy_xy[k] = -g_x_0_xxy_xy[k] * cd_y[k] + g_x_0_xxy_xyy[k];

            g_x_0_xxyy_xz[k] = -g_x_0_xxy_xz[k] * cd_y[k] + g_x_0_xxy_xyz[k];

            g_x_0_xxyy_yy[k] = -g_x_0_xxy_yy[k] * cd_y[k] + g_x_0_xxy_yyy[k];

            g_x_0_xxyy_yz[k] = -g_x_0_xxy_yz[k] * cd_y[k] + g_x_0_xxy_yyz[k];

            g_x_0_xxyy_zz[k] = -g_x_0_xxy_zz[k] * cd_y[k] + g_x_0_xxy_yzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_xxyz_xx, g_x_0_xxyz_xy, g_x_0_xxyz_xz, g_x_0_xxyz_yy, g_x_0_xxyz_yz, g_x_0_xxyz_zz, g_x_0_xxz_xx, g_x_0_xxz_xxy, g_x_0_xxz_xy, g_x_0_xxz_xyy, g_x_0_xxz_xyz, g_x_0_xxz_xz, g_x_0_xxz_yy, g_x_0_xxz_yyy, g_x_0_xxz_yyz, g_x_0_xxz_yz, g_x_0_xxz_yzz, g_x_0_xxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyz_xx[k] = -g_x_0_xxz_xx[k] * cd_y[k] + g_x_0_xxz_xxy[k];

            g_x_0_xxyz_xy[k] = -g_x_0_xxz_xy[k] * cd_y[k] + g_x_0_xxz_xyy[k];

            g_x_0_xxyz_xz[k] = -g_x_0_xxz_xz[k] * cd_y[k] + g_x_0_xxz_xyz[k];

            g_x_0_xxyz_yy[k] = -g_x_0_xxz_yy[k] * cd_y[k] + g_x_0_xxz_yyy[k];

            g_x_0_xxyz_yz[k] = -g_x_0_xxz_yz[k] * cd_y[k] + g_x_0_xxz_yyz[k];

            g_x_0_xxyz_zz[k] = -g_x_0_xxz_zz[k] * cd_y[k] + g_x_0_xxz_yzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_x_0_xxz_xx, g_x_0_xxz_xxz, g_x_0_xxz_xy, g_x_0_xxz_xyz, g_x_0_xxz_xz, g_x_0_xxz_xzz, g_x_0_xxz_yy, g_x_0_xxz_yyz, g_x_0_xxz_yz, g_x_0_xxz_yzz, g_x_0_xxz_zz, g_x_0_xxz_zzz, g_x_0_xxzz_xx, g_x_0_xxzz_xy, g_x_0_xxzz_xz, g_x_0_xxzz_yy, g_x_0_xxzz_yz, g_x_0_xxzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxzz_xx[k] = -g_x_0_xxz_xx[k] * cd_z[k] + g_x_0_xxz_xxz[k];

            g_x_0_xxzz_xy[k] = -g_x_0_xxz_xy[k] * cd_z[k] + g_x_0_xxz_xyz[k];

            g_x_0_xxzz_xz[k] = -g_x_0_xxz_xz[k] * cd_z[k] + g_x_0_xxz_xzz[k];

            g_x_0_xxzz_yy[k] = -g_x_0_xxz_yy[k] * cd_z[k] + g_x_0_xxz_yyz[k];

            g_x_0_xxzz_yz[k] = -g_x_0_xxz_yz[k] * cd_z[k] + g_x_0_xxz_yzz[k];

            g_x_0_xxzz_zz[k] = -g_x_0_xxz_zz[k] * cd_z[k] + g_x_0_xxz_zzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_x_0_xyy_xx, g_x_0_xyy_xxy, g_x_0_xyy_xy, g_x_0_xyy_xyy, g_x_0_xyy_xyz, g_x_0_xyy_xz, g_x_0_xyy_yy, g_x_0_xyy_yyy, g_x_0_xyy_yyz, g_x_0_xyy_yz, g_x_0_xyy_yzz, g_x_0_xyy_zz, g_x_0_xyyy_xx, g_x_0_xyyy_xy, g_x_0_xyyy_xz, g_x_0_xyyy_yy, g_x_0_xyyy_yz, g_x_0_xyyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyy_xx[k] = -g_x_0_xyy_xx[k] * cd_y[k] + g_x_0_xyy_xxy[k];

            g_x_0_xyyy_xy[k] = -g_x_0_xyy_xy[k] * cd_y[k] + g_x_0_xyy_xyy[k];

            g_x_0_xyyy_xz[k] = -g_x_0_xyy_xz[k] * cd_y[k] + g_x_0_xyy_xyz[k];

            g_x_0_xyyy_yy[k] = -g_x_0_xyy_yy[k] * cd_y[k] + g_x_0_xyy_yyy[k];

            g_x_0_xyyy_yz[k] = -g_x_0_xyy_yz[k] * cd_y[k] + g_x_0_xyy_yyz[k];

            g_x_0_xyyy_zz[k] = -g_x_0_xyy_zz[k] * cd_y[k] + g_x_0_xyy_yzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 47);

        #pragma omp simd aligned(cd_y, g_x_0_xyyz_xx, g_x_0_xyyz_xy, g_x_0_xyyz_xz, g_x_0_xyyz_yy, g_x_0_xyyz_yz, g_x_0_xyyz_zz, g_x_0_xyz_xx, g_x_0_xyz_xxy, g_x_0_xyz_xy, g_x_0_xyz_xyy, g_x_0_xyz_xyz, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yyy, g_x_0_xyz_yyz, g_x_0_xyz_yz, g_x_0_xyz_yzz, g_x_0_xyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyz_xx[k] = -g_x_0_xyz_xx[k] * cd_y[k] + g_x_0_xyz_xxy[k];

            g_x_0_xyyz_xy[k] = -g_x_0_xyz_xy[k] * cd_y[k] + g_x_0_xyz_xyy[k];

            g_x_0_xyyz_xz[k] = -g_x_0_xyz_xz[k] * cd_y[k] + g_x_0_xyz_xyz[k];

            g_x_0_xyyz_yy[k] = -g_x_0_xyz_yy[k] * cd_y[k] + g_x_0_xyz_yyy[k];

            g_x_0_xyyz_yz[k] = -g_x_0_xyz_yz[k] * cd_y[k] + g_x_0_xyz_yyz[k];

            g_x_0_xyyz_zz[k] = -g_x_0_xyz_zz[k] * cd_y[k] + g_x_0_xyz_yzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 53);

        #pragma omp simd aligned(cd_y, g_x_0_xyzz_xx, g_x_0_xyzz_xy, g_x_0_xyzz_xz, g_x_0_xyzz_yy, g_x_0_xyzz_yz, g_x_0_xyzz_zz, g_x_0_xzz_xx, g_x_0_xzz_xxy, g_x_0_xzz_xy, g_x_0_xzz_xyy, g_x_0_xzz_xyz, g_x_0_xzz_xz, g_x_0_xzz_yy, g_x_0_xzz_yyy, g_x_0_xzz_yyz, g_x_0_xzz_yz, g_x_0_xzz_yzz, g_x_0_xzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyzz_xx[k] = -g_x_0_xzz_xx[k] * cd_y[k] + g_x_0_xzz_xxy[k];

            g_x_0_xyzz_xy[k] = -g_x_0_xzz_xy[k] * cd_y[k] + g_x_0_xzz_xyy[k];

            g_x_0_xyzz_xz[k] = -g_x_0_xzz_xz[k] * cd_y[k] + g_x_0_xzz_xyz[k];

            g_x_0_xyzz_yy[k] = -g_x_0_xzz_yy[k] * cd_y[k] + g_x_0_xzz_yyy[k];

            g_x_0_xyzz_yz[k] = -g_x_0_xzz_yz[k] * cd_y[k] + g_x_0_xzz_yyz[k];

            g_x_0_xyzz_zz[k] = -g_x_0_xzz_zz[k] * cd_y[k] + g_x_0_xzz_yzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_x_0_xzz_xx, g_x_0_xzz_xxz, g_x_0_xzz_xy, g_x_0_xzz_xyz, g_x_0_xzz_xz, g_x_0_xzz_xzz, g_x_0_xzz_yy, g_x_0_xzz_yyz, g_x_0_xzz_yz, g_x_0_xzz_yzz, g_x_0_xzz_zz, g_x_0_xzz_zzz, g_x_0_xzzz_xx, g_x_0_xzzz_xy, g_x_0_xzzz_xz, g_x_0_xzzz_yy, g_x_0_xzzz_yz, g_x_0_xzzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzzz_xx[k] = -g_x_0_xzz_xx[k] * cd_z[k] + g_x_0_xzz_xxz[k];

            g_x_0_xzzz_xy[k] = -g_x_0_xzz_xy[k] * cd_z[k] + g_x_0_xzz_xyz[k];

            g_x_0_xzzz_xz[k] = -g_x_0_xzz_xz[k] * cd_z[k] + g_x_0_xzz_xzz[k];

            g_x_0_xzzz_yy[k] = -g_x_0_xzz_yy[k] * cd_z[k] + g_x_0_xzz_yyz[k];

            g_x_0_xzzz_yz[k] = -g_x_0_xzz_yz[k] * cd_z[k] + g_x_0_xzz_yzz[k];

            g_x_0_xzzz_zz[k] = -g_x_0_xzz_zz[k] * cd_z[k] + g_x_0_xzz_zzz[k];
        }

        /// Set up 60-66 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 65);

        #pragma omp simd aligned(cd_y, g_x_0_yyy_xx, g_x_0_yyy_xxy, g_x_0_yyy_xy, g_x_0_yyy_xyy, g_x_0_yyy_xyz, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yyy, g_x_0_yyy_yyz, g_x_0_yyy_yz, g_x_0_yyy_yzz, g_x_0_yyy_zz, g_x_0_yyyy_xx, g_x_0_yyyy_xy, g_x_0_yyyy_xz, g_x_0_yyyy_yy, g_x_0_yyyy_yz, g_x_0_yyyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyy_xx[k] = -g_x_0_yyy_xx[k] * cd_y[k] + g_x_0_yyy_xxy[k];

            g_x_0_yyyy_xy[k] = -g_x_0_yyy_xy[k] * cd_y[k] + g_x_0_yyy_xyy[k];

            g_x_0_yyyy_xz[k] = -g_x_0_yyy_xz[k] * cd_y[k] + g_x_0_yyy_xyz[k];

            g_x_0_yyyy_yy[k] = -g_x_0_yyy_yy[k] * cd_y[k] + g_x_0_yyy_yyy[k];

            g_x_0_yyyy_yz[k] = -g_x_0_yyy_yz[k] * cd_y[k] + g_x_0_yyy_yyz[k];

            g_x_0_yyyy_zz[k] = -g_x_0_yyy_zz[k] * cd_y[k] + g_x_0_yyy_yzz[k];
        }

        /// Set up 66-72 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 71);

        #pragma omp simd aligned(cd_y, g_x_0_yyyz_xx, g_x_0_yyyz_xy, g_x_0_yyyz_xz, g_x_0_yyyz_yy, g_x_0_yyyz_yz, g_x_0_yyyz_zz, g_x_0_yyz_xx, g_x_0_yyz_xxy, g_x_0_yyz_xy, g_x_0_yyz_xyy, g_x_0_yyz_xyz, g_x_0_yyz_xz, g_x_0_yyz_yy, g_x_0_yyz_yyy, g_x_0_yyz_yyz, g_x_0_yyz_yz, g_x_0_yyz_yzz, g_x_0_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyz_xx[k] = -g_x_0_yyz_xx[k] * cd_y[k] + g_x_0_yyz_xxy[k];

            g_x_0_yyyz_xy[k] = -g_x_0_yyz_xy[k] * cd_y[k] + g_x_0_yyz_xyy[k];

            g_x_0_yyyz_xz[k] = -g_x_0_yyz_xz[k] * cd_y[k] + g_x_0_yyz_xyz[k];

            g_x_0_yyyz_yy[k] = -g_x_0_yyz_yy[k] * cd_y[k] + g_x_0_yyz_yyy[k];

            g_x_0_yyyz_yz[k] = -g_x_0_yyz_yz[k] * cd_y[k] + g_x_0_yyz_yyz[k];

            g_x_0_yyyz_zz[k] = -g_x_0_yyz_zz[k] * cd_y[k] + g_x_0_yyz_yzz[k];
        }

        /// Set up 72-78 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 74);

        auto g_x_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 77);

        #pragma omp simd aligned(cd_y, g_x_0_yyzz_xx, g_x_0_yyzz_xy, g_x_0_yyzz_xz, g_x_0_yyzz_yy, g_x_0_yyzz_yz, g_x_0_yyzz_zz, g_x_0_yzz_xx, g_x_0_yzz_xxy, g_x_0_yzz_xy, g_x_0_yzz_xyy, g_x_0_yzz_xyz, g_x_0_yzz_xz, g_x_0_yzz_yy, g_x_0_yzz_yyy, g_x_0_yzz_yyz, g_x_0_yzz_yz, g_x_0_yzz_yzz, g_x_0_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyzz_xx[k] = -g_x_0_yzz_xx[k] * cd_y[k] + g_x_0_yzz_xxy[k];

            g_x_0_yyzz_xy[k] = -g_x_0_yzz_xy[k] * cd_y[k] + g_x_0_yzz_xyy[k];

            g_x_0_yyzz_xz[k] = -g_x_0_yzz_xz[k] * cd_y[k] + g_x_0_yzz_xyz[k];

            g_x_0_yyzz_yy[k] = -g_x_0_yzz_yy[k] * cd_y[k] + g_x_0_yzz_yyy[k];

            g_x_0_yyzz_yz[k] = -g_x_0_yzz_yz[k] * cd_y[k] + g_x_0_yzz_yyz[k];

            g_x_0_yyzz_zz[k] = -g_x_0_yzz_zz[k] * cd_y[k] + g_x_0_yzz_yzz[k];
        }

        /// Set up 78-84 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_x_0_yzzz_xx, g_x_0_yzzz_xy, g_x_0_yzzz_xz, g_x_0_yzzz_yy, g_x_0_yzzz_yz, g_x_0_yzzz_zz, g_x_0_zzz_xx, g_x_0_zzz_xxy, g_x_0_zzz_xy, g_x_0_zzz_xyy, g_x_0_zzz_xyz, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yyy, g_x_0_zzz_yyz, g_x_0_zzz_yz, g_x_0_zzz_yzz, g_x_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzzz_xx[k] = -g_x_0_zzz_xx[k] * cd_y[k] + g_x_0_zzz_xxy[k];

            g_x_0_yzzz_xy[k] = -g_x_0_zzz_xy[k] * cd_y[k] + g_x_0_zzz_xyy[k];

            g_x_0_yzzz_xz[k] = -g_x_0_zzz_xz[k] * cd_y[k] + g_x_0_zzz_xyz[k];

            g_x_0_yzzz_yy[k] = -g_x_0_zzz_yy[k] * cd_y[k] + g_x_0_zzz_yyy[k];

            g_x_0_yzzz_yz[k] = -g_x_0_zzz_yz[k] * cd_y[k] + g_x_0_zzz_yyz[k];

            g_x_0_yzzz_zz[k] = -g_x_0_zzz_zz[k] * cd_y[k] + g_x_0_zzz_yzz[k];
        }

        /// Set up 84-90 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 0 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_x_0_zzz_xx, g_x_0_zzz_xxz, g_x_0_zzz_xy, g_x_0_zzz_xyz, g_x_0_zzz_xz, g_x_0_zzz_xzz, g_x_0_zzz_yy, g_x_0_zzz_yyz, g_x_0_zzz_yz, g_x_0_zzz_yzz, g_x_0_zzz_zz, g_x_0_zzz_zzz, g_x_0_zzzz_xx, g_x_0_zzzz_xy, g_x_0_zzzz_xz, g_x_0_zzzz_yy, g_x_0_zzzz_yz, g_x_0_zzzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzzz_xx[k] = -g_x_0_zzz_xx[k] * cd_z[k] + g_x_0_zzz_xxz[k];

            g_x_0_zzzz_xy[k] = -g_x_0_zzz_xy[k] * cd_z[k] + g_x_0_zzz_xyz[k];

            g_x_0_zzzz_xz[k] = -g_x_0_zzz_xz[k] * cd_z[k] + g_x_0_zzz_xzz[k];

            g_x_0_zzzz_yy[k] = -g_x_0_zzz_yy[k] * cd_z[k] + g_x_0_zzz_yyz[k];

            g_x_0_zzzz_yz[k] = -g_x_0_zzz_yz[k] * cd_z[k] + g_x_0_zzz_yzz[k];

            g_x_0_zzzz_zz[k] = -g_x_0_zzz_zz[k] * cd_z[k] + g_x_0_zzz_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 0);

        auto g_y_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 1);

        auto g_y_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 2);

        auto g_y_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 3);

        auto g_y_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 4);

        auto g_y_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xxx_xx, g_y_0_xxx_xxx, g_y_0_xxx_xxy, g_y_0_xxx_xxz, g_y_0_xxx_xy, g_y_0_xxx_xyy, g_y_0_xxx_xyz, g_y_0_xxx_xz, g_y_0_xxx_xzz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz, g_y_0_xxxx_xx, g_y_0_xxxx_xy, g_y_0_xxxx_xz, g_y_0_xxxx_yy, g_y_0_xxxx_yz, g_y_0_xxxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxx_xx[k] = -g_y_0_xxx_xx[k] * cd_x[k] + g_y_0_xxx_xxx[k];

            g_y_0_xxxx_xy[k] = -g_y_0_xxx_xy[k] * cd_x[k] + g_y_0_xxx_xxy[k];

            g_y_0_xxxx_xz[k] = -g_y_0_xxx_xz[k] * cd_x[k] + g_y_0_xxx_xxz[k];

            g_y_0_xxxx_yy[k] = -g_y_0_xxx_yy[k] * cd_x[k] + g_y_0_xxx_xyy[k];

            g_y_0_xxxx_yz[k] = -g_y_0_xxx_yz[k] * cd_x[k] + g_y_0_xxx_xyz[k];

            g_y_0_xxxx_zz[k] = -g_y_0_xxx_zz[k] * cd_x[k] + g_y_0_xxx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 6);

        auto g_y_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 7);

        auto g_y_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 8);

        auto g_y_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 9);

        auto g_y_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 10);

        auto g_y_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_y_0_xxxy_xx, g_y_0_xxxy_xy, g_y_0_xxxy_xz, g_y_0_xxxy_yy, g_y_0_xxxy_yz, g_y_0_xxxy_zz, g_y_0_xxy_xx, g_y_0_xxy_xxx, g_y_0_xxy_xxy, g_y_0_xxy_xxz, g_y_0_xxy_xy, g_y_0_xxy_xyy, g_y_0_xxy_xyz, g_y_0_xxy_xz, g_y_0_xxy_xzz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxy_xx[k] = -g_y_0_xxy_xx[k] * cd_x[k] + g_y_0_xxy_xxx[k];

            g_y_0_xxxy_xy[k] = -g_y_0_xxy_xy[k] * cd_x[k] + g_y_0_xxy_xxy[k];

            g_y_0_xxxy_xz[k] = -g_y_0_xxy_xz[k] * cd_x[k] + g_y_0_xxy_xxz[k];

            g_y_0_xxxy_yy[k] = -g_y_0_xxy_yy[k] * cd_x[k] + g_y_0_xxy_xyy[k];

            g_y_0_xxxy_yz[k] = -g_y_0_xxy_yz[k] * cd_x[k] + g_y_0_xxy_xyz[k];

            g_y_0_xxxy_zz[k] = -g_y_0_xxy_zz[k] * cd_x[k] + g_y_0_xxy_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 12);

        auto g_y_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 13);

        auto g_y_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 14);

        auto g_y_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 15);

        auto g_y_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 16);

        auto g_y_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_y_0_xxxz_xx, g_y_0_xxxz_xy, g_y_0_xxxz_xz, g_y_0_xxxz_yy, g_y_0_xxxz_yz, g_y_0_xxxz_zz, g_y_0_xxz_xx, g_y_0_xxz_xxx, g_y_0_xxz_xxy, g_y_0_xxz_xxz, g_y_0_xxz_xy, g_y_0_xxz_xyy, g_y_0_xxz_xyz, g_y_0_xxz_xz, g_y_0_xxz_xzz, g_y_0_xxz_yy, g_y_0_xxz_yz, g_y_0_xxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxz_xx[k] = -g_y_0_xxz_xx[k] * cd_x[k] + g_y_0_xxz_xxx[k];

            g_y_0_xxxz_xy[k] = -g_y_0_xxz_xy[k] * cd_x[k] + g_y_0_xxz_xxy[k];

            g_y_0_xxxz_xz[k] = -g_y_0_xxz_xz[k] * cd_x[k] + g_y_0_xxz_xxz[k];

            g_y_0_xxxz_yy[k] = -g_y_0_xxz_yy[k] * cd_x[k] + g_y_0_xxz_xyy[k];

            g_y_0_xxxz_yz[k] = -g_y_0_xxz_yz[k] * cd_x[k] + g_y_0_xxz_xyz[k];

            g_y_0_xxxz_zz[k] = -g_y_0_xxz_zz[k] * cd_x[k] + g_y_0_xxz_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 18);

        auto g_y_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 19);

        auto g_y_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 20);

        auto g_y_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 21);

        auto g_y_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 22);

        auto g_y_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 23);

        #pragma omp simd aligned(cd_x, g_y_0_xxyy_xx, g_y_0_xxyy_xy, g_y_0_xxyy_xz, g_y_0_xxyy_yy, g_y_0_xxyy_yz, g_y_0_xxyy_zz, g_y_0_xyy_xx, g_y_0_xyy_xxx, g_y_0_xyy_xxy, g_y_0_xyy_xxz, g_y_0_xyy_xy, g_y_0_xyy_xyy, g_y_0_xyy_xyz, g_y_0_xyy_xz, g_y_0_xyy_xzz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyy_xx[k] = -g_y_0_xyy_xx[k] * cd_x[k] + g_y_0_xyy_xxx[k];

            g_y_0_xxyy_xy[k] = -g_y_0_xyy_xy[k] * cd_x[k] + g_y_0_xyy_xxy[k];

            g_y_0_xxyy_xz[k] = -g_y_0_xyy_xz[k] * cd_x[k] + g_y_0_xyy_xxz[k];

            g_y_0_xxyy_yy[k] = -g_y_0_xyy_yy[k] * cd_x[k] + g_y_0_xyy_xyy[k];

            g_y_0_xxyy_yz[k] = -g_y_0_xyy_yz[k] * cd_x[k] + g_y_0_xyy_xyz[k];

            g_y_0_xxyy_zz[k] = -g_y_0_xyy_zz[k] * cd_x[k] + g_y_0_xyy_xzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 24);

        auto g_y_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 25);

        auto g_y_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 26);

        auto g_y_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 27);

        auto g_y_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 28);

        auto g_y_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xxyz_xx, g_y_0_xxyz_xy, g_y_0_xxyz_xz, g_y_0_xxyz_yy, g_y_0_xxyz_yz, g_y_0_xxyz_zz, g_y_0_xyz_xx, g_y_0_xyz_xxx, g_y_0_xyz_xxy, g_y_0_xyz_xxz, g_y_0_xyz_xy, g_y_0_xyz_xyy, g_y_0_xyz_xyz, g_y_0_xyz_xz, g_y_0_xyz_xzz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyz_xx[k] = -g_y_0_xyz_xx[k] * cd_x[k] + g_y_0_xyz_xxx[k];

            g_y_0_xxyz_xy[k] = -g_y_0_xyz_xy[k] * cd_x[k] + g_y_0_xyz_xxy[k];

            g_y_0_xxyz_xz[k] = -g_y_0_xyz_xz[k] * cd_x[k] + g_y_0_xyz_xxz[k];

            g_y_0_xxyz_yy[k] = -g_y_0_xyz_yy[k] * cd_x[k] + g_y_0_xyz_xyy[k];

            g_y_0_xxyz_yz[k] = -g_y_0_xyz_yz[k] * cd_x[k] + g_y_0_xyz_xyz[k];

            g_y_0_xxyz_zz[k] = -g_y_0_xyz_zz[k] * cd_x[k] + g_y_0_xyz_xzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 30);

        auto g_y_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 31);

        auto g_y_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 32);

        auto g_y_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 33);

        auto g_y_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 34);

        auto g_y_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_y_0_xxzz_xx, g_y_0_xxzz_xy, g_y_0_xxzz_xz, g_y_0_xxzz_yy, g_y_0_xxzz_yz, g_y_0_xxzz_zz, g_y_0_xzz_xx, g_y_0_xzz_xxx, g_y_0_xzz_xxy, g_y_0_xzz_xxz, g_y_0_xzz_xy, g_y_0_xzz_xyy, g_y_0_xzz_xyz, g_y_0_xzz_xz, g_y_0_xzz_xzz, g_y_0_xzz_yy, g_y_0_xzz_yz, g_y_0_xzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxzz_xx[k] = -g_y_0_xzz_xx[k] * cd_x[k] + g_y_0_xzz_xxx[k];

            g_y_0_xxzz_xy[k] = -g_y_0_xzz_xy[k] * cd_x[k] + g_y_0_xzz_xxy[k];

            g_y_0_xxzz_xz[k] = -g_y_0_xzz_xz[k] * cd_x[k] + g_y_0_xzz_xxz[k];

            g_y_0_xxzz_yy[k] = -g_y_0_xzz_yy[k] * cd_x[k] + g_y_0_xzz_xyy[k];

            g_y_0_xxzz_yz[k] = -g_y_0_xzz_yz[k] * cd_x[k] + g_y_0_xzz_xyz[k];

            g_y_0_xxzz_zz[k] = -g_y_0_xzz_zz[k] * cd_x[k] + g_y_0_xzz_xzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 36);

        auto g_y_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 37);

        auto g_y_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 38);

        auto g_y_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 39);

        auto g_y_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 40);

        auto g_y_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_y_0_xyyy_xx, g_y_0_xyyy_xy, g_y_0_xyyy_xz, g_y_0_xyyy_yy, g_y_0_xyyy_yz, g_y_0_xyyy_zz, g_y_0_yyy_xx, g_y_0_yyy_xxx, g_y_0_yyy_xxy, g_y_0_yyy_xxz, g_y_0_yyy_xy, g_y_0_yyy_xyy, g_y_0_yyy_xyz, g_y_0_yyy_xz, g_y_0_yyy_xzz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyy_xx[k] = -g_y_0_yyy_xx[k] * cd_x[k] + g_y_0_yyy_xxx[k];

            g_y_0_xyyy_xy[k] = -g_y_0_yyy_xy[k] * cd_x[k] + g_y_0_yyy_xxy[k];

            g_y_0_xyyy_xz[k] = -g_y_0_yyy_xz[k] * cd_x[k] + g_y_0_yyy_xxz[k];

            g_y_0_xyyy_yy[k] = -g_y_0_yyy_yy[k] * cd_x[k] + g_y_0_yyy_xyy[k];

            g_y_0_xyyy_yz[k] = -g_y_0_yyy_yz[k] * cd_x[k] + g_y_0_yyy_xyz[k];

            g_y_0_xyyy_zz[k] = -g_y_0_yyy_zz[k] * cd_x[k] + g_y_0_yyy_xzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 42);

        auto g_y_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 43);

        auto g_y_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 44);

        auto g_y_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 45);

        auto g_y_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 46);

        auto g_y_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 47);

        #pragma omp simd aligned(cd_x, g_y_0_xyyz_xx, g_y_0_xyyz_xy, g_y_0_xyyz_xz, g_y_0_xyyz_yy, g_y_0_xyyz_yz, g_y_0_xyyz_zz, g_y_0_yyz_xx, g_y_0_yyz_xxx, g_y_0_yyz_xxy, g_y_0_yyz_xxz, g_y_0_yyz_xy, g_y_0_yyz_xyy, g_y_0_yyz_xyz, g_y_0_yyz_xz, g_y_0_yyz_xzz, g_y_0_yyz_yy, g_y_0_yyz_yz, g_y_0_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyz_xx[k] = -g_y_0_yyz_xx[k] * cd_x[k] + g_y_0_yyz_xxx[k];

            g_y_0_xyyz_xy[k] = -g_y_0_yyz_xy[k] * cd_x[k] + g_y_0_yyz_xxy[k];

            g_y_0_xyyz_xz[k] = -g_y_0_yyz_xz[k] * cd_x[k] + g_y_0_yyz_xxz[k];

            g_y_0_xyyz_yy[k] = -g_y_0_yyz_yy[k] * cd_x[k] + g_y_0_yyz_xyy[k];

            g_y_0_xyyz_yz[k] = -g_y_0_yyz_yz[k] * cd_x[k] + g_y_0_yyz_xyz[k];

            g_y_0_xyyz_zz[k] = -g_y_0_yyz_zz[k] * cd_x[k] + g_y_0_yyz_xzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 48);

        auto g_y_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 49);

        auto g_y_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 50);

        auto g_y_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 51);

        auto g_y_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 52);

        auto g_y_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 53);

        #pragma omp simd aligned(cd_x, g_y_0_xyzz_xx, g_y_0_xyzz_xy, g_y_0_xyzz_xz, g_y_0_xyzz_yy, g_y_0_xyzz_yz, g_y_0_xyzz_zz, g_y_0_yzz_xx, g_y_0_yzz_xxx, g_y_0_yzz_xxy, g_y_0_yzz_xxz, g_y_0_yzz_xy, g_y_0_yzz_xyy, g_y_0_yzz_xyz, g_y_0_yzz_xz, g_y_0_yzz_xzz, g_y_0_yzz_yy, g_y_0_yzz_yz, g_y_0_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyzz_xx[k] = -g_y_0_yzz_xx[k] * cd_x[k] + g_y_0_yzz_xxx[k];

            g_y_0_xyzz_xy[k] = -g_y_0_yzz_xy[k] * cd_x[k] + g_y_0_yzz_xxy[k];

            g_y_0_xyzz_xz[k] = -g_y_0_yzz_xz[k] * cd_x[k] + g_y_0_yzz_xxz[k];

            g_y_0_xyzz_yy[k] = -g_y_0_yzz_yy[k] * cd_x[k] + g_y_0_yzz_xyy[k];

            g_y_0_xyzz_yz[k] = -g_y_0_yzz_yz[k] * cd_x[k] + g_y_0_yzz_xyz[k];

            g_y_0_xyzz_zz[k] = -g_y_0_yzz_zz[k] * cd_x[k] + g_y_0_yzz_xzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 54);

        auto g_y_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 55);

        auto g_y_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 56);

        auto g_y_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 57);

        auto g_y_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 58);

        auto g_y_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_y_0_xzzz_xx, g_y_0_xzzz_xy, g_y_0_xzzz_xz, g_y_0_xzzz_yy, g_y_0_xzzz_yz, g_y_0_xzzz_zz, g_y_0_zzz_xx, g_y_0_zzz_xxx, g_y_0_zzz_xxy, g_y_0_zzz_xxz, g_y_0_zzz_xy, g_y_0_zzz_xyy, g_y_0_zzz_xyz, g_y_0_zzz_xz, g_y_0_zzz_xzz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzzz_xx[k] = -g_y_0_zzz_xx[k] * cd_x[k] + g_y_0_zzz_xxx[k];

            g_y_0_xzzz_xy[k] = -g_y_0_zzz_xy[k] * cd_x[k] + g_y_0_zzz_xxy[k];

            g_y_0_xzzz_xz[k] = -g_y_0_zzz_xz[k] * cd_x[k] + g_y_0_zzz_xxz[k];

            g_y_0_xzzz_yy[k] = -g_y_0_zzz_yy[k] * cd_x[k] + g_y_0_zzz_xyy[k];

            g_y_0_xzzz_yz[k] = -g_y_0_zzz_yz[k] * cd_x[k] + g_y_0_zzz_xyz[k];

            g_y_0_xzzz_zz[k] = -g_y_0_zzz_zz[k] * cd_x[k] + g_y_0_zzz_xzz[k];
        }

        /// Set up 60-66 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 60);

        auto g_y_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 61);

        auto g_y_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 62);

        auto g_y_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 63);

        auto g_y_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 64);

        auto g_y_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 65);

        #pragma omp simd aligned(cd_y, g_y_0_yyy_xx, g_y_0_yyy_xxy, g_y_0_yyy_xy, g_y_0_yyy_xyy, g_y_0_yyy_xyz, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yyy, g_y_0_yyy_yyz, g_y_0_yyy_yz, g_y_0_yyy_yzz, g_y_0_yyy_zz, g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz, g_yyy_xx, g_yyy_xy, g_yyy_xz, g_yyy_yy, g_yyy_yz, g_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyy_xx[k] = -g_yyy_xx[k] - g_y_0_yyy_xx[k] * cd_y[k] + g_y_0_yyy_xxy[k];

            g_y_0_yyyy_xy[k] = -g_yyy_xy[k] - g_y_0_yyy_xy[k] * cd_y[k] + g_y_0_yyy_xyy[k];

            g_y_0_yyyy_xz[k] = -g_yyy_xz[k] - g_y_0_yyy_xz[k] * cd_y[k] + g_y_0_yyy_xyz[k];

            g_y_0_yyyy_yy[k] = -g_yyy_yy[k] - g_y_0_yyy_yy[k] * cd_y[k] + g_y_0_yyy_yyy[k];

            g_y_0_yyyy_yz[k] = -g_yyy_yz[k] - g_y_0_yyy_yz[k] * cd_y[k] + g_y_0_yyy_yyz[k];

            g_y_0_yyyy_zz[k] = -g_yyy_zz[k] - g_y_0_yyy_zz[k] * cd_y[k] + g_y_0_yyy_yzz[k];
        }

        /// Set up 66-72 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 66);

        auto g_y_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 67);

        auto g_y_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 68);

        auto g_y_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 69);

        auto g_y_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 70);

        auto g_y_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 71);

        #pragma omp simd aligned(cd_z, g_y_0_yyy_xx, g_y_0_yyy_xxz, g_y_0_yyy_xy, g_y_0_yyy_xyz, g_y_0_yyy_xz, g_y_0_yyy_xzz, g_y_0_yyy_yy, g_y_0_yyy_yyz, g_y_0_yyy_yz, g_y_0_yyy_yzz, g_y_0_yyy_zz, g_y_0_yyy_zzz, g_y_0_yyyz_xx, g_y_0_yyyz_xy, g_y_0_yyyz_xz, g_y_0_yyyz_yy, g_y_0_yyyz_yz, g_y_0_yyyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyz_xx[k] = -g_y_0_yyy_xx[k] * cd_z[k] + g_y_0_yyy_xxz[k];

            g_y_0_yyyz_xy[k] = -g_y_0_yyy_xy[k] * cd_z[k] + g_y_0_yyy_xyz[k];

            g_y_0_yyyz_xz[k] = -g_y_0_yyy_xz[k] * cd_z[k] + g_y_0_yyy_xzz[k];

            g_y_0_yyyz_yy[k] = -g_y_0_yyy_yy[k] * cd_z[k] + g_y_0_yyy_yyz[k];

            g_y_0_yyyz_yz[k] = -g_y_0_yyy_yz[k] * cd_z[k] + g_y_0_yyy_yzz[k];

            g_y_0_yyyz_zz[k] = -g_y_0_yyy_zz[k] * cd_z[k] + g_y_0_yyy_zzz[k];
        }

        /// Set up 72-78 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 72);

        auto g_y_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 73);

        auto g_y_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 74);

        auto g_y_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 75);

        auto g_y_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 76);

        auto g_y_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 77);

        #pragma omp simd aligned(cd_z, g_y_0_yyz_xx, g_y_0_yyz_xxz, g_y_0_yyz_xy, g_y_0_yyz_xyz, g_y_0_yyz_xz, g_y_0_yyz_xzz, g_y_0_yyz_yy, g_y_0_yyz_yyz, g_y_0_yyz_yz, g_y_0_yyz_yzz, g_y_0_yyz_zz, g_y_0_yyz_zzz, g_y_0_yyzz_xx, g_y_0_yyzz_xy, g_y_0_yyzz_xz, g_y_0_yyzz_yy, g_y_0_yyzz_yz, g_y_0_yyzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyzz_xx[k] = -g_y_0_yyz_xx[k] * cd_z[k] + g_y_0_yyz_xxz[k];

            g_y_0_yyzz_xy[k] = -g_y_0_yyz_xy[k] * cd_z[k] + g_y_0_yyz_xyz[k];

            g_y_0_yyzz_xz[k] = -g_y_0_yyz_xz[k] * cd_z[k] + g_y_0_yyz_xzz[k];

            g_y_0_yyzz_yy[k] = -g_y_0_yyz_yy[k] * cd_z[k] + g_y_0_yyz_yyz[k];

            g_y_0_yyzz_yz[k] = -g_y_0_yyz_yz[k] * cd_z[k] + g_y_0_yyz_yzz[k];

            g_y_0_yyzz_zz[k] = -g_y_0_yyz_zz[k] * cd_z[k] + g_y_0_yyz_zzz[k];
        }

        /// Set up 78-84 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 78);

        auto g_y_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 79);

        auto g_y_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 80);

        auto g_y_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 81);

        auto g_y_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 82);

        auto g_y_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 83);

        #pragma omp simd aligned(cd_z, g_y_0_yzz_xx, g_y_0_yzz_xxz, g_y_0_yzz_xy, g_y_0_yzz_xyz, g_y_0_yzz_xz, g_y_0_yzz_xzz, g_y_0_yzz_yy, g_y_0_yzz_yyz, g_y_0_yzz_yz, g_y_0_yzz_yzz, g_y_0_yzz_zz, g_y_0_yzz_zzz, g_y_0_yzzz_xx, g_y_0_yzzz_xy, g_y_0_yzzz_xz, g_y_0_yzzz_yy, g_y_0_yzzz_yz, g_y_0_yzzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzzz_xx[k] = -g_y_0_yzz_xx[k] * cd_z[k] + g_y_0_yzz_xxz[k];

            g_y_0_yzzz_xy[k] = -g_y_0_yzz_xy[k] * cd_z[k] + g_y_0_yzz_xyz[k];

            g_y_0_yzzz_xz[k] = -g_y_0_yzz_xz[k] * cd_z[k] + g_y_0_yzz_xzz[k];

            g_y_0_yzzz_yy[k] = -g_y_0_yzz_yy[k] * cd_z[k] + g_y_0_yzz_yyz[k];

            g_y_0_yzzz_yz[k] = -g_y_0_yzz_yz[k] * cd_z[k] + g_y_0_yzz_yzz[k];

            g_y_0_yzzz_zz[k] = -g_y_0_yzz_zz[k] * cd_z[k] + g_y_0_yzz_zzz[k];
        }

        /// Set up 84-90 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 90 * acomps  + 84);

        auto g_y_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 85);

        auto g_y_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 86);

        auto g_y_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 90 * acomps  + 87);

        auto g_y_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 88);

        auto g_y_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 90 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_y_0_zzz_xx, g_y_0_zzz_xxz, g_y_0_zzz_xy, g_y_0_zzz_xyz, g_y_0_zzz_xz, g_y_0_zzz_xzz, g_y_0_zzz_yy, g_y_0_zzz_yyz, g_y_0_zzz_yz, g_y_0_zzz_yzz, g_y_0_zzz_zz, g_y_0_zzz_zzz, g_y_0_zzzz_xx, g_y_0_zzzz_xy, g_y_0_zzzz_xz, g_y_0_zzzz_yy, g_y_0_zzzz_yz, g_y_0_zzzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzzz_xx[k] = -g_y_0_zzz_xx[k] * cd_z[k] + g_y_0_zzz_xxz[k];

            g_y_0_zzzz_xy[k] = -g_y_0_zzz_xy[k] * cd_z[k] + g_y_0_zzz_xyz[k];

            g_y_0_zzzz_xz[k] = -g_y_0_zzz_xz[k] * cd_z[k] + g_y_0_zzz_xzz[k];

            g_y_0_zzzz_yy[k] = -g_y_0_zzz_yy[k] * cd_z[k] + g_y_0_zzz_yyz[k];

            g_y_0_zzzz_yz[k] = -g_y_0_zzz_yz[k] * cd_z[k] + g_y_0_zzz_yzz[k];

            g_y_0_zzzz_zz[k] = -g_y_0_zzz_zz[k] * cd_z[k] + g_y_0_zzz_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 0);

        auto g_z_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 1);

        auto g_z_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 2);

        auto g_z_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 3);

        auto g_z_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 4);

        auto g_z_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xxx_xx, g_z_0_xxx_xxx, g_z_0_xxx_xxy, g_z_0_xxx_xxz, g_z_0_xxx_xy, g_z_0_xxx_xyy, g_z_0_xxx_xyz, g_z_0_xxx_xz, g_z_0_xxx_xzz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz, g_z_0_xxxx_xx, g_z_0_xxxx_xy, g_z_0_xxxx_xz, g_z_0_xxxx_yy, g_z_0_xxxx_yz, g_z_0_xxxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxx_xx[k] = -g_z_0_xxx_xx[k] * cd_x[k] + g_z_0_xxx_xxx[k];

            g_z_0_xxxx_xy[k] = -g_z_0_xxx_xy[k] * cd_x[k] + g_z_0_xxx_xxy[k];

            g_z_0_xxxx_xz[k] = -g_z_0_xxx_xz[k] * cd_x[k] + g_z_0_xxx_xxz[k];

            g_z_0_xxxx_yy[k] = -g_z_0_xxx_yy[k] * cd_x[k] + g_z_0_xxx_xyy[k];

            g_z_0_xxxx_yz[k] = -g_z_0_xxx_yz[k] * cd_x[k] + g_z_0_xxx_xyz[k];

            g_z_0_xxxx_zz[k] = -g_z_0_xxx_zz[k] * cd_x[k] + g_z_0_xxx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 6);

        auto g_z_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 7);

        auto g_z_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 8);

        auto g_z_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 9);

        auto g_z_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 10);

        auto g_z_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_z_0_xxxy_xx, g_z_0_xxxy_xy, g_z_0_xxxy_xz, g_z_0_xxxy_yy, g_z_0_xxxy_yz, g_z_0_xxxy_zz, g_z_0_xxy_xx, g_z_0_xxy_xxx, g_z_0_xxy_xxy, g_z_0_xxy_xxz, g_z_0_xxy_xy, g_z_0_xxy_xyy, g_z_0_xxy_xyz, g_z_0_xxy_xz, g_z_0_xxy_xzz, g_z_0_xxy_yy, g_z_0_xxy_yz, g_z_0_xxy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxy_xx[k] = -g_z_0_xxy_xx[k] * cd_x[k] + g_z_0_xxy_xxx[k];

            g_z_0_xxxy_xy[k] = -g_z_0_xxy_xy[k] * cd_x[k] + g_z_0_xxy_xxy[k];

            g_z_0_xxxy_xz[k] = -g_z_0_xxy_xz[k] * cd_x[k] + g_z_0_xxy_xxz[k];

            g_z_0_xxxy_yy[k] = -g_z_0_xxy_yy[k] * cd_x[k] + g_z_0_xxy_xyy[k];

            g_z_0_xxxy_yz[k] = -g_z_0_xxy_yz[k] * cd_x[k] + g_z_0_xxy_xyz[k];

            g_z_0_xxxy_zz[k] = -g_z_0_xxy_zz[k] * cd_x[k] + g_z_0_xxy_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 12);

        auto g_z_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 13);

        auto g_z_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 14);

        auto g_z_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 15);

        auto g_z_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 16);

        auto g_z_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_z_0_xxxz_xx, g_z_0_xxxz_xy, g_z_0_xxxz_xz, g_z_0_xxxz_yy, g_z_0_xxxz_yz, g_z_0_xxxz_zz, g_z_0_xxz_xx, g_z_0_xxz_xxx, g_z_0_xxz_xxy, g_z_0_xxz_xxz, g_z_0_xxz_xy, g_z_0_xxz_xyy, g_z_0_xxz_xyz, g_z_0_xxz_xz, g_z_0_xxz_xzz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxz_xx[k] = -g_z_0_xxz_xx[k] * cd_x[k] + g_z_0_xxz_xxx[k];

            g_z_0_xxxz_xy[k] = -g_z_0_xxz_xy[k] * cd_x[k] + g_z_0_xxz_xxy[k];

            g_z_0_xxxz_xz[k] = -g_z_0_xxz_xz[k] * cd_x[k] + g_z_0_xxz_xxz[k];

            g_z_0_xxxz_yy[k] = -g_z_0_xxz_yy[k] * cd_x[k] + g_z_0_xxz_xyy[k];

            g_z_0_xxxz_yz[k] = -g_z_0_xxz_yz[k] * cd_x[k] + g_z_0_xxz_xyz[k];

            g_z_0_xxxz_zz[k] = -g_z_0_xxz_zz[k] * cd_x[k] + g_z_0_xxz_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 18);

        auto g_z_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 19);

        auto g_z_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 20);

        auto g_z_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 21);

        auto g_z_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 22);

        auto g_z_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 23);

        #pragma omp simd aligned(cd_x, g_z_0_xxyy_xx, g_z_0_xxyy_xy, g_z_0_xxyy_xz, g_z_0_xxyy_yy, g_z_0_xxyy_yz, g_z_0_xxyy_zz, g_z_0_xyy_xx, g_z_0_xyy_xxx, g_z_0_xyy_xxy, g_z_0_xyy_xxz, g_z_0_xyy_xy, g_z_0_xyy_xyy, g_z_0_xyy_xyz, g_z_0_xyy_xz, g_z_0_xyy_xzz, g_z_0_xyy_yy, g_z_0_xyy_yz, g_z_0_xyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyy_xx[k] = -g_z_0_xyy_xx[k] * cd_x[k] + g_z_0_xyy_xxx[k];

            g_z_0_xxyy_xy[k] = -g_z_0_xyy_xy[k] * cd_x[k] + g_z_0_xyy_xxy[k];

            g_z_0_xxyy_xz[k] = -g_z_0_xyy_xz[k] * cd_x[k] + g_z_0_xyy_xxz[k];

            g_z_0_xxyy_yy[k] = -g_z_0_xyy_yy[k] * cd_x[k] + g_z_0_xyy_xyy[k];

            g_z_0_xxyy_yz[k] = -g_z_0_xyy_yz[k] * cd_x[k] + g_z_0_xyy_xyz[k];

            g_z_0_xxyy_zz[k] = -g_z_0_xyy_zz[k] * cd_x[k] + g_z_0_xyy_xzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 24);

        auto g_z_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 25);

        auto g_z_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 26);

        auto g_z_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 27);

        auto g_z_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 28);

        auto g_z_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xxyz_xx, g_z_0_xxyz_xy, g_z_0_xxyz_xz, g_z_0_xxyz_yy, g_z_0_xxyz_yz, g_z_0_xxyz_zz, g_z_0_xyz_xx, g_z_0_xyz_xxx, g_z_0_xyz_xxy, g_z_0_xyz_xxz, g_z_0_xyz_xy, g_z_0_xyz_xyy, g_z_0_xyz_xyz, g_z_0_xyz_xz, g_z_0_xyz_xzz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyz_xx[k] = -g_z_0_xyz_xx[k] * cd_x[k] + g_z_0_xyz_xxx[k];

            g_z_0_xxyz_xy[k] = -g_z_0_xyz_xy[k] * cd_x[k] + g_z_0_xyz_xxy[k];

            g_z_0_xxyz_xz[k] = -g_z_0_xyz_xz[k] * cd_x[k] + g_z_0_xyz_xxz[k];

            g_z_0_xxyz_yy[k] = -g_z_0_xyz_yy[k] * cd_x[k] + g_z_0_xyz_xyy[k];

            g_z_0_xxyz_yz[k] = -g_z_0_xyz_yz[k] * cd_x[k] + g_z_0_xyz_xyz[k];

            g_z_0_xxyz_zz[k] = -g_z_0_xyz_zz[k] * cd_x[k] + g_z_0_xyz_xzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 30);

        auto g_z_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 31);

        auto g_z_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 32);

        auto g_z_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 33);

        auto g_z_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 34);

        auto g_z_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_z_0_xxzz_xx, g_z_0_xxzz_xy, g_z_0_xxzz_xz, g_z_0_xxzz_yy, g_z_0_xxzz_yz, g_z_0_xxzz_zz, g_z_0_xzz_xx, g_z_0_xzz_xxx, g_z_0_xzz_xxy, g_z_0_xzz_xxz, g_z_0_xzz_xy, g_z_0_xzz_xyy, g_z_0_xzz_xyz, g_z_0_xzz_xz, g_z_0_xzz_xzz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxzz_xx[k] = -g_z_0_xzz_xx[k] * cd_x[k] + g_z_0_xzz_xxx[k];

            g_z_0_xxzz_xy[k] = -g_z_0_xzz_xy[k] * cd_x[k] + g_z_0_xzz_xxy[k];

            g_z_0_xxzz_xz[k] = -g_z_0_xzz_xz[k] * cd_x[k] + g_z_0_xzz_xxz[k];

            g_z_0_xxzz_yy[k] = -g_z_0_xzz_yy[k] * cd_x[k] + g_z_0_xzz_xyy[k];

            g_z_0_xxzz_yz[k] = -g_z_0_xzz_yz[k] * cd_x[k] + g_z_0_xzz_xyz[k];

            g_z_0_xxzz_zz[k] = -g_z_0_xzz_zz[k] * cd_x[k] + g_z_0_xzz_xzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 36);

        auto g_z_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 37);

        auto g_z_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 38);

        auto g_z_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 39);

        auto g_z_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 40);

        auto g_z_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_z_0_xyyy_xx, g_z_0_xyyy_xy, g_z_0_xyyy_xz, g_z_0_xyyy_yy, g_z_0_xyyy_yz, g_z_0_xyyy_zz, g_z_0_yyy_xx, g_z_0_yyy_xxx, g_z_0_yyy_xxy, g_z_0_yyy_xxz, g_z_0_yyy_xy, g_z_0_yyy_xyy, g_z_0_yyy_xyz, g_z_0_yyy_xz, g_z_0_yyy_xzz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyy_xx[k] = -g_z_0_yyy_xx[k] * cd_x[k] + g_z_0_yyy_xxx[k];

            g_z_0_xyyy_xy[k] = -g_z_0_yyy_xy[k] * cd_x[k] + g_z_0_yyy_xxy[k];

            g_z_0_xyyy_xz[k] = -g_z_0_yyy_xz[k] * cd_x[k] + g_z_0_yyy_xxz[k];

            g_z_0_xyyy_yy[k] = -g_z_0_yyy_yy[k] * cd_x[k] + g_z_0_yyy_xyy[k];

            g_z_0_xyyy_yz[k] = -g_z_0_yyy_yz[k] * cd_x[k] + g_z_0_yyy_xyz[k];

            g_z_0_xyyy_zz[k] = -g_z_0_yyy_zz[k] * cd_x[k] + g_z_0_yyy_xzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 42);

        auto g_z_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 43);

        auto g_z_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 44);

        auto g_z_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 45);

        auto g_z_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 46);

        auto g_z_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 47);

        #pragma omp simd aligned(cd_x, g_z_0_xyyz_xx, g_z_0_xyyz_xy, g_z_0_xyyz_xz, g_z_0_xyyz_yy, g_z_0_xyyz_yz, g_z_0_xyyz_zz, g_z_0_yyz_xx, g_z_0_yyz_xxx, g_z_0_yyz_xxy, g_z_0_yyz_xxz, g_z_0_yyz_xy, g_z_0_yyz_xyy, g_z_0_yyz_xyz, g_z_0_yyz_xz, g_z_0_yyz_xzz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyz_xx[k] = -g_z_0_yyz_xx[k] * cd_x[k] + g_z_0_yyz_xxx[k];

            g_z_0_xyyz_xy[k] = -g_z_0_yyz_xy[k] * cd_x[k] + g_z_0_yyz_xxy[k];

            g_z_0_xyyz_xz[k] = -g_z_0_yyz_xz[k] * cd_x[k] + g_z_0_yyz_xxz[k];

            g_z_0_xyyz_yy[k] = -g_z_0_yyz_yy[k] * cd_x[k] + g_z_0_yyz_xyy[k];

            g_z_0_xyyz_yz[k] = -g_z_0_yyz_yz[k] * cd_x[k] + g_z_0_yyz_xyz[k];

            g_z_0_xyyz_zz[k] = -g_z_0_yyz_zz[k] * cd_x[k] + g_z_0_yyz_xzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 48);

        auto g_z_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 49);

        auto g_z_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 50);

        auto g_z_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 51);

        auto g_z_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 52);

        auto g_z_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 53);

        #pragma omp simd aligned(cd_x, g_z_0_xyzz_xx, g_z_0_xyzz_xy, g_z_0_xyzz_xz, g_z_0_xyzz_yy, g_z_0_xyzz_yz, g_z_0_xyzz_zz, g_z_0_yzz_xx, g_z_0_yzz_xxx, g_z_0_yzz_xxy, g_z_0_yzz_xxz, g_z_0_yzz_xy, g_z_0_yzz_xyy, g_z_0_yzz_xyz, g_z_0_yzz_xz, g_z_0_yzz_xzz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyzz_xx[k] = -g_z_0_yzz_xx[k] * cd_x[k] + g_z_0_yzz_xxx[k];

            g_z_0_xyzz_xy[k] = -g_z_0_yzz_xy[k] * cd_x[k] + g_z_0_yzz_xxy[k];

            g_z_0_xyzz_xz[k] = -g_z_0_yzz_xz[k] * cd_x[k] + g_z_0_yzz_xxz[k];

            g_z_0_xyzz_yy[k] = -g_z_0_yzz_yy[k] * cd_x[k] + g_z_0_yzz_xyy[k];

            g_z_0_xyzz_yz[k] = -g_z_0_yzz_yz[k] * cd_x[k] + g_z_0_yzz_xyz[k];

            g_z_0_xyzz_zz[k] = -g_z_0_yzz_zz[k] * cd_x[k] + g_z_0_yzz_xzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 54);

        auto g_z_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 55);

        auto g_z_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 56);

        auto g_z_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 57);

        auto g_z_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 58);

        auto g_z_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_z_0_xzzz_xx, g_z_0_xzzz_xy, g_z_0_xzzz_xz, g_z_0_xzzz_yy, g_z_0_xzzz_yz, g_z_0_xzzz_zz, g_z_0_zzz_xx, g_z_0_zzz_xxx, g_z_0_zzz_xxy, g_z_0_zzz_xxz, g_z_0_zzz_xy, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xz, g_z_0_zzz_xzz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzzz_xx[k] = -g_z_0_zzz_xx[k] * cd_x[k] + g_z_0_zzz_xxx[k];

            g_z_0_xzzz_xy[k] = -g_z_0_zzz_xy[k] * cd_x[k] + g_z_0_zzz_xxy[k];

            g_z_0_xzzz_xz[k] = -g_z_0_zzz_xz[k] * cd_x[k] + g_z_0_zzz_xxz[k];

            g_z_0_xzzz_yy[k] = -g_z_0_zzz_yy[k] * cd_x[k] + g_z_0_zzz_xyy[k];

            g_z_0_xzzz_yz[k] = -g_z_0_zzz_yz[k] * cd_x[k] + g_z_0_zzz_xyz[k];

            g_z_0_xzzz_zz[k] = -g_z_0_zzz_zz[k] * cd_x[k] + g_z_0_zzz_xzz[k];
        }

        /// Set up 60-66 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 60);

        auto g_z_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 61);

        auto g_z_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 62);

        auto g_z_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 63);

        auto g_z_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 64);

        auto g_z_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 65);

        #pragma omp simd aligned(cd_y, g_z_0_yyy_xx, g_z_0_yyy_xxy, g_z_0_yyy_xy, g_z_0_yyy_xyy, g_z_0_yyy_xyz, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yyy, g_z_0_yyy_yyz, g_z_0_yyy_yz, g_z_0_yyy_yzz, g_z_0_yyy_zz, g_z_0_yyyy_xx, g_z_0_yyyy_xy, g_z_0_yyyy_xz, g_z_0_yyyy_yy, g_z_0_yyyy_yz, g_z_0_yyyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyy_xx[k] = -g_z_0_yyy_xx[k] * cd_y[k] + g_z_0_yyy_xxy[k];

            g_z_0_yyyy_xy[k] = -g_z_0_yyy_xy[k] * cd_y[k] + g_z_0_yyy_xyy[k];

            g_z_0_yyyy_xz[k] = -g_z_0_yyy_xz[k] * cd_y[k] + g_z_0_yyy_xyz[k];

            g_z_0_yyyy_yy[k] = -g_z_0_yyy_yy[k] * cd_y[k] + g_z_0_yyy_yyy[k];

            g_z_0_yyyy_yz[k] = -g_z_0_yyy_yz[k] * cd_y[k] + g_z_0_yyy_yyz[k];

            g_z_0_yyyy_zz[k] = -g_z_0_yyy_zz[k] * cd_y[k] + g_z_0_yyy_yzz[k];
        }

        /// Set up 66-72 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 66);

        auto g_z_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 67);

        auto g_z_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 68);

        auto g_z_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 69);

        auto g_z_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 70);

        auto g_z_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 71);

        #pragma omp simd aligned(cd_y, g_z_0_yyyz_xx, g_z_0_yyyz_xy, g_z_0_yyyz_xz, g_z_0_yyyz_yy, g_z_0_yyyz_yz, g_z_0_yyyz_zz, g_z_0_yyz_xx, g_z_0_yyz_xxy, g_z_0_yyz_xy, g_z_0_yyz_xyy, g_z_0_yyz_xyz, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yyy, g_z_0_yyz_yyz, g_z_0_yyz_yz, g_z_0_yyz_yzz, g_z_0_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyz_xx[k] = -g_z_0_yyz_xx[k] * cd_y[k] + g_z_0_yyz_xxy[k];

            g_z_0_yyyz_xy[k] = -g_z_0_yyz_xy[k] * cd_y[k] + g_z_0_yyz_xyy[k];

            g_z_0_yyyz_xz[k] = -g_z_0_yyz_xz[k] * cd_y[k] + g_z_0_yyz_xyz[k];

            g_z_0_yyyz_yy[k] = -g_z_0_yyz_yy[k] * cd_y[k] + g_z_0_yyz_yyy[k];

            g_z_0_yyyz_yz[k] = -g_z_0_yyz_yz[k] * cd_y[k] + g_z_0_yyz_yyz[k];

            g_z_0_yyyz_zz[k] = -g_z_0_yyz_zz[k] * cd_y[k] + g_z_0_yyz_yzz[k];
        }

        /// Set up 72-78 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 72);

        auto g_z_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 73);

        auto g_z_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 74);

        auto g_z_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 75);

        auto g_z_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 76);

        auto g_z_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 77);

        #pragma omp simd aligned(cd_y, g_z_0_yyzz_xx, g_z_0_yyzz_xy, g_z_0_yyzz_xz, g_z_0_yyzz_yy, g_z_0_yyzz_yz, g_z_0_yyzz_zz, g_z_0_yzz_xx, g_z_0_yzz_xxy, g_z_0_yzz_xy, g_z_0_yzz_xyy, g_z_0_yzz_xyz, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yyy, g_z_0_yzz_yyz, g_z_0_yzz_yz, g_z_0_yzz_yzz, g_z_0_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyzz_xx[k] = -g_z_0_yzz_xx[k] * cd_y[k] + g_z_0_yzz_xxy[k];

            g_z_0_yyzz_xy[k] = -g_z_0_yzz_xy[k] * cd_y[k] + g_z_0_yzz_xyy[k];

            g_z_0_yyzz_xz[k] = -g_z_0_yzz_xz[k] * cd_y[k] + g_z_0_yzz_xyz[k];

            g_z_0_yyzz_yy[k] = -g_z_0_yzz_yy[k] * cd_y[k] + g_z_0_yzz_yyy[k];

            g_z_0_yyzz_yz[k] = -g_z_0_yzz_yz[k] * cd_y[k] + g_z_0_yzz_yyz[k];

            g_z_0_yyzz_zz[k] = -g_z_0_yzz_zz[k] * cd_y[k] + g_z_0_yzz_yzz[k];
        }

        /// Set up 78-84 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 78);

        auto g_z_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 79);

        auto g_z_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 80);

        auto g_z_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 81);

        auto g_z_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 82);

        auto g_z_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_z_0_yzzz_xx, g_z_0_yzzz_xy, g_z_0_yzzz_xz, g_z_0_yzzz_yy, g_z_0_yzzz_yz, g_z_0_yzzz_zz, g_z_0_zzz_xx, g_z_0_zzz_xxy, g_z_0_zzz_xy, g_z_0_zzz_xyy, g_z_0_zzz_xyz, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yz, g_z_0_zzz_yzz, g_z_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzzz_xx[k] = -g_z_0_zzz_xx[k] * cd_y[k] + g_z_0_zzz_xxy[k];

            g_z_0_yzzz_xy[k] = -g_z_0_zzz_xy[k] * cd_y[k] + g_z_0_zzz_xyy[k];

            g_z_0_yzzz_xz[k] = -g_z_0_zzz_xz[k] * cd_y[k] + g_z_0_zzz_xyz[k];

            g_z_0_yzzz_yy[k] = -g_z_0_zzz_yy[k] * cd_y[k] + g_z_0_zzz_yyy[k];

            g_z_0_yzzz_yz[k] = -g_z_0_zzz_yz[k] * cd_y[k] + g_z_0_zzz_yyz[k];

            g_z_0_yzzz_zz[k] = -g_z_0_zzz_zz[k] * cd_y[k] + g_z_0_zzz_yzz[k];
        }

        /// Set up 84-90 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 180 * acomps  + 84);

        auto g_z_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 85);

        auto g_z_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 86);

        auto g_z_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 180 * acomps  + 87);

        auto g_z_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 88);

        auto g_z_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 180 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_z_0_zzz_xx, g_z_0_zzz_xxz, g_z_0_zzz_xy, g_z_0_zzz_xyz, g_z_0_zzz_xz, g_z_0_zzz_xzz, g_z_0_zzz_yy, g_z_0_zzz_yyz, g_z_0_zzz_yz, g_z_0_zzz_yzz, g_z_0_zzz_zz, g_z_0_zzz_zzz, g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz, g_zzz_xx, g_zzz_xy, g_zzz_xz, g_zzz_yy, g_zzz_yz, g_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzzz_xx[k] = -g_zzz_xx[k] - g_z_0_zzz_xx[k] * cd_z[k] + g_z_0_zzz_xxz[k];

            g_z_0_zzzz_xy[k] = -g_zzz_xy[k] - g_z_0_zzz_xy[k] * cd_z[k] + g_z_0_zzz_xyz[k];

            g_z_0_zzzz_xz[k] = -g_zzz_xz[k] - g_z_0_zzz_xz[k] * cd_z[k] + g_z_0_zzz_xzz[k];

            g_z_0_zzzz_yy[k] = -g_zzz_yy[k] - g_z_0_zzz_yy[k] * cd_z[k] + g_z_0_zzz_yyz[k];

            g_z_0_zzzz_yz[k] = -g_zzz_yz[k] - g_z_0_zzz_yz[k] * cd_z[k] + g_z_0_zzz_yzz[k];

            g_z_0_zzzz_zz[k] = -g_zzz_zz[k] - g_z_0_zzz_zz[k] * cd_z[k] + g_z_0_zzz_zzz[k];
        }
    }
}

} // t3ceri namespace

