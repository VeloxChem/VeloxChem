#include "ElectronRepulsionGeom0010ContrRecXXGP.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxgp(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgp,
                                            const size_t idx_xxfp,
                                            const size_t idx_geom_10_xxfp,
                                            const size_t idx_geom_10_xxfd,
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
            /// Set up components of auxilary buffer : SSFP

            const auto fp_off = idx_xxfp + (i * bcomps + j) * 30;

            auto g_xxx_x = cbuffer.data(fp_off + 0);

            auto g_xxx_y = cbuffer.data(fp_off + 1);

            auto g_xxx_z = cbuffer.data(fp_off + 2);

            auto g_yyy_x = cbuffer.data(fp_off + 18);

            auto g_yyy_y = cbuffer.data(fp_off + 19);

            auto g_yyy_z = cbuffer.data(fp_off + 20);

            auto g_zzz_x = cbuffer.data(fp_off + 27);

            auto g_zzz_y = cbuffer.data(fp_off + 28);

            auto g_zzz_z = cbuffer.data(fp_off + 29);

            /// Set up components of auxilary buffer : SSFP

            const auto fp_geom_10_off = idx_geom_10_xxfp + (i * bcomps + j) * 30;

            auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxy_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxy_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_yyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_yyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_yyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_yyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_yzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_yzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_zzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 0);

            auto g_y_0_xxx_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 1);

            auto g_y_0_xxx_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 2);

            auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 3);

            auto g_y_0_xxy_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 4);

            auto g_y_0_xxy_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 5);

            auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 6);

            auto g_y_0_xxz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 7);

            auto g_y_0_xxz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 8);

            auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 9);

            auto g_y_0_xyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 10);

            auto g_y_0_xyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 11);

            auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 12);

            auto g_y_0_xyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 13);

            auto g_y_0_xyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 14);

            auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 15);

            auto g_y_0_xzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 16);

            auto g_y_0_xzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 17);

            auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 18);

            auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 19);

            auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 20);

            auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 21);

            auto g_y_0_yyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 22);

            auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 23);

            auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 24);

            auto g_y_0_yzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 25);

            auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 26);

            auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 27);

            auto g_y_0_zzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 28);

            auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 29);

            auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 0);

            auto g_z_0_xxx_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 1);

            auto g_z_0_xxx_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 2);

            auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 3);

            auto g_z_0_xxy_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 4);

            auto g_z_0_xxy_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 5);

            auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 6);

            auto g_z_0_xxz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 7);

            auto g_z_0_xxz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 8);

            auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 9);

            auto g_z_0_xyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 10);

            auto g_z_0_xyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 11);

            auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 12);

            auto g_z_0_xyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 13);

            auto g_z_0_xyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 14);

            auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 15);

            auto g_z_0_xzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 16);

            auto g_z_0_xzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 17);

            auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 18);

            auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 19);

            auto g_z_0_yyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 20);

            auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 21);

            auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 22);

            auto g_z_0_yyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 23);

            auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 24);

            auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 25);

            auto g_z_0_yzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 26);

            auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 27);

            auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 28);

            auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 29);

            /// Set up components of auxilary buffer : SSFD

            const auto fd_geom_10_off = idx_geom_10_xxfd + (i * bcomps + j) * 60;

            auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 0);

            auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 1);

            auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 2);

            auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 6);

            auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 7);

            auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 8);

            auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 12);

            auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 13);

            auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 14);

            auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 18);

            auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 19);

            auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 20);

            auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 24);

            auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 25);

            auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 26);

            auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 30);

            auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 31);

            auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 32);

            auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 36);

            auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 37);

            auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 38);

            auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 39);

            auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 40);

            auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 41);

            auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 42);

            auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 43);

            auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 44);

            auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 46);

            auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 47);

            auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 48);

            auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 49);

            auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 50);

            auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 52);

            auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 53);

            auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 54);

            auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 55);

            auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 56);

            auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 58);

            auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps * bcomps + 59);

            auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 0);

            auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 1);

            auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 2);

            auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 6);

            auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 7);

            auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 8);

            auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 12);

            auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 13);

            auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 14);

            auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 18);

            auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 19);

            auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 20);

            auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 24);

            auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 25);

            auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 26);

            auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 30);

            auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 31);

            auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 32);

            auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 36);

            auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 37);

            auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 38);

            auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 39);

            auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 40);

            auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 42);

            auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 43);

            auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 44);

            auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 45);

            auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 46);

            auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 48);

            auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 49);

            auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 50);

            auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 51);

            auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 52);

            auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 54);

            auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 55);

            auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 56);

            auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 57);

            auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 58);

            auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps * bcomps + 59);

            /// set up bra offset for contr_buffer_xxgp

            const auto gp_geom_10_off = idx_geom_10_xxgp + (i * bcomps + j) * 45;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_x_0_xxx_x, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_y, g_x_0_xxx_z, g_x_0_xxxx_x, g_x_0_xxxx_y, g_x_0_xxxx_z, g_xxx_x, g_xxx_y, g_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_x[k] = -g_xxx_x[k] - g_x_0_xxx_x[k] * cd_x[k] + g_x_0_xxx_xx[k];

                g_x_0_xxxx_y[k] = -g_xxx_y[k] - g_x_0_xxx_y[k] * cd_x[k] + g_x_0_xxx_xy[k];

                g_x_0_xxxx_z[k] = -g_xxx_z[k] - g_x_0_xxx_z[k] * cd_x[k] + g_x_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_y, g_x_0_xxx_x, g_x_0_xxx_xy, g_x_0_xxx_y, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_z, g_x_0_xxxy_x, g_x_0_xxxy_y, g_x_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_x[k] = -g_x_0_xxx_x[k] * cd_y[k] + g_x_0_xxx_xy[k];

                g_x_0_xxxy_y[k] = -g_x_0_xxx_y[k] * cd_y[k] + g_x_0_xxx_yy[k];

                g_x_0_xxxy_z[k] = -g_x_0_xxx_z[k] * cd_y[k] + g_x_0_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_z, g_x_0_xxx_x, g_x_0_xxx_xz, g_x_0_xxx_y, g_x_0_xxx_yz, g_x_0_xxx_z, g_x_0_xxx_zz, g_x_0_xxxz_x, g_x_0_xxxz_y, g_x_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_x[k] = -g_x_0_xxx_x[k] * cd_z[k] + g_x_0_xxx_xz[k];

                g_x_0_xxxz_y[k] = -g_x_0_xxx_y[k] * cd_z[k] + g_x_0_xxx_yz[k];

                g_x_0_xxxz_z[k] = -g_x_0_xxx_z[k] * cd_z[k] + g_x_0_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxy_x, g_x_0_xxy_xy, g_x_0_xxy_y, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_z, g_x_0_xxyy_x, g_x_0_xxyy_y, g_x_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_x[k] = -g_x_0_xxy_x[k] * cd_y[k] + g_x_0_xxy_xy[k];

                g_x_0_xxyy_y[k] = -g_x_0_xxy_y[k] * cd_y[k] + g_x_0_xxy_yy[k];

                g_x_0_xxyy_z[k] = -g_x_0_xxy_z[k] * cd_y[k] + g_x_0_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_y, g_x_0_xxyz_x, g_x_0_xxyz_y, g_x_0_xxyz_z, g_x_0_xxz_x, g_x_0_xxz_xy, g_x_0_xxz_y, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_x[k] = -g_x_0_xxz_x[k] * cd_y[k] + g_x_0_xxz_xy[k];

                g_x_0_xxyz_y[k] = -g_x_0_xxz_y[k] * cd_y[k] + g_x_0_xxz_yy[k];

                g_x_0_xxyz_z[k] = -g_x_0_xxz_z[k] * cd_y[k] + g_x_0_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_x_0_xxz_x, g_x_0_xxz_xz, g_x_0_xxz_y, g_x_0_xxz_yz, g_x_0_xxz_z, g_x_0_xxz_zz, g_x_0_xxzz_x, g_x_0_xxzz_y, g_x_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_x[k] = -g_x_0_xxz_x[k] * cd_z[k] + g_x_0_xxz_xz[k];

                g_x_0_xxzz_y[k] = -g_x_0_xxz_y[k] * cd_z[k] + g_x_0_xxz_yz[k];

                g_x_0_xxzz_z[k] = -g_x_0_xxz_z[k] * cd_z[k] + g_x_0_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xyyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_y, g_x_0_xyy_x, g_x_0_xyy_xy, g_x_0_xyy_y, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_z, g_x_0_xyyy_x, g_x_0_xyyy_y, g_x_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_x[k] = -g_x_0_xyy_x[k] * cd_y[k] + g_x_0_xyy_xy[k];

                g_x_0_xyyy_y[k] = -g_x_0_xyy_y[k] * cd_y[k] + g_x_0_xyy_yy[k];

                g_x_0_xyyy_z[k] = -g_x_0_xyy_z[k] * cd_y[k] + g_x_0_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xyyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_xyyz_x, g_x_0_xyyz_y, g_x_0_xyyz_z, g_x_0_xyz_x, g_x_0_xyz_xy, g_x_0_xyz_y, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_x[k] = -g_x_0_xyz_x[k] * cd_y[k] + g_x_0_xyz_xy[k];

                g_x_0_xyyz_y[k] = -g_x_0_xyz_y[k] * cd_y[k] + g_x_0_xyz_yy[k];

                g_x_0_xyyz_z[k] = -g_x_0_xyz_z[k] * cd_y[k] + g_x_0_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xyzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_y, g_x_0_xyzz_x, g_x_0_xyzz_y, g_x_0_xyzz_z, g_x_0_xzz_x, g_x_0_xzz_xy, g_x_0_xzz_y, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_x[k] = -g_x_0_xzz_x[k] * cd_y[k] + g_x_0_xzz_xy[k];

                g_x_0_xyzz_y[k] = -g_x_0_xzz_y[k] * cd_y[k] + g_x_0_xzz_yy[k];

                g_x_0_xyzz_z[k] = -g_x_0_xzz_z[k] * cd_y[k] + g_x_0_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_xzz_x, g_x_0_xzz_xz, g_x_0_xzz_y, g_x_0_xzz_yz, g_x_0_xzz_z, g_x_0_xzz_zz, g_x_0_xzzz_x, g_x_0_xzzz_y, g_x_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_x[k] = -g_x_0_xzz_x[k] * cd_z[k] + g_x_0_xzz_xz[k];

                g_x_0_xzzz_y[k] = -g_x_0_xzz_y[k] * cd_z[k] + g_x_0_xzz_yz[k];

                g_x_0_xzzz_z[k] = -g_x_0_xzz_z[k] * cd_z[k] + g_x_0_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_yyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_yyyy_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_y, g_x_0_yyy_x, g_x_0_yyy_xy, g_x_0_yyy_y, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_z, g_x_0_yyyy_x, g_x_0_yyyy_y, g_x_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_x[k] = -g_x_0_yyy_x[k] * cd_y[k] + g_x_0_yyy_xy[k];

                g_x_0_yyyy_y[k] = -g_x_0_yyy_y[k] * cd_y[k] + g_x_0_yyy_yy[k];

                g_x_0_yyyy_z[k] = -g_x_0_yyy_z[k] * cd_y[k] + g_x_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_yyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_yyyz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_y, g_x_0_yyyz_x, g_x_0_yyyz_y, g_x_0_yyyz_z, g_x_0_yyz_x, g_x_0_yyz_xy, g_x_0_yyz_y, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_x[k] = -g_x_0_yyz_x[k] * cd_y[k] + g_x_0_yyz_xy[k];

                g_x_0_yyyz_y[k] = -g_x_0_yyz_y[k] * cd_y[k] + g_x_0_yyz_yy[k];

                g_x_0_yyyz_z[k] = -g_x_0_yyz_z[k] * cd_y[k] + g_x_0_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_yyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_yyzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_y, g_x_0_yyzz_x, g_x_0_yyzz_y, g_x_0_yyzz_z, g_x_0_yzz_x, g_x_0_yzz_xy, g_x_0_yzz_y, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_x[k] = -g_x_0_yzz_x[k] * cd_y[k] + g_x_0_yzz_xy[k];

                g_x_0_yyzz_y[k] = -g_x_0_yzz_y[k] * cd_y[k] + g_x_0_yzz_yy[k];

                g_x_0_yyzz_z[k] = -g_x_0_yzz_z[k] * cd_y[k] + g_x_0_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_yzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_yzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_yzzz_x, g_x_0_yzzz_y, g_x_0_yzzz_z, g_x_0_zzz_x, g_x_0_zzz_xy, g_x_0_zzz_y, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_x[k] = -g_x_0_zzz_x[k] * cd_y[k] + g_x_0_zzz_xy[k];

                g_x_0_yzzz_y[k] = -g_x_0_zzz_y[k] * cd_y[k] + g_x_0_zzz_yy[k];

                g_x_0_yzzz_z[k] = -g_x_0_zzz_z[k] * cd_y[k] + g_x_0_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_zzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_zzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_zzz_x, g_x_0_zzz_xz, g_x_0_zzz_y, g_x_0_zzz_yz, g_x_0_zzz_z, g_x_0_zzz_zz, g_x_0_zzzz_x, g_x_0_zzzz_y, g_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_x[k] = -g_x_0_zzz_x[k] * cd_z[k] + g_x_0_zzz_xz[k];

                g_x_0_zzzz_y[k] = -g_x_0_zzz_y[k] * cd_z[k] + g_x_0_zzz_yz[k];

                g_x_0_zzzz_z[k] = -g_x_0_zzz_z[k] * cd_z[k] + g_x_0_zzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 0);

            auto g_y_0_xxxx_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 1);

            auto g_y_0_xxxx_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxx_x, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_y, g_y_0_xxx_z, g_y_0_xxxx_x, g_y_0_xxxx_y, g_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_x[k] = -g_y_0_xxx_x[k] * cd_x[k] + g_y_0_xxx_xx[k];

                g_y_0_xxxx_y[k] = -g_y_0_xxx_y[k] * cd_x[k] + g_y_0_xxx_xy[k];

                g_y_0_xxxx_z[k] = -g_y_0_xxx_z[k] * cd_x[k] + g_y_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 3);

            auto g_y_0_xxxy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 4);

            auto g_y_0_xxxy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxy_x, g_y_0_xxxy_y, g_y_0_xxxy_z, g_y_0_xxy_x, g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_y, g_y_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_x[k] = -g_y_0_xxy_x[k] * cd_x[k] + g_y_0_xxy_xx[k];

                g_y_0_xxxy_y[k] = -g_y_0_xxy_y[k] * cd_x[k] + g_y_0_xxy_xy[k];

                g_y_0_xxxy_z[k] = -g_y_0_xxy_z[k] * cd_x[k] + g_y_0_xxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 6);

            auto g_y_0_xxxz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 7);

            auto g_y_0_xxxz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xxxz_x, g_y_0_xxxz_y, g_y_0_xxxz_z, g_y_0_xxz_x, g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_y, g_y_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_x[k] = -g_y_0_xxz_x[k] * cd_x[k] + g_y_0_xxz_xx[k];

                g_y_0_xxxz_y[k] = -g_y_0_xxz_y[k] * cd_x[k] + g_y_0_xxz_xy[k];

                g_y_0_xxxz_z[k] = -g_y_0_xxz_z[k] * cd_x[k] + g_y_0_xxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 9);

            auto g_y_0_xxyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 10);

            auto g_y_0_xxyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxyy_x, g_y_0_xxyy_y, g_y_0_xxyy_z, g_y_0_xyy_x, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_y, g_y_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_x[k] = -g_y_0_xyy_x[k] * cd_x[k] + g_y_0_xyy_xx[k];

                g_y_0_xxyy_y[k] = -g_y_0_xyy_y[k] * cd_x[k] + g_y_0_xyy_xy[k];

                g_y_0_xxyy_z[k] = -g_y_0_xyy_z[k] * cd_x[k] + g_y_0_xyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 12);

            auto g_y_0_xxyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 13);

            auto g_y_0_xxyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxyz_x, g_y_0_xxyz_y, g_y_0_xxyz_z, g_y_0_xyz_x, g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_y, g_y_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_x[k] = -g_y_0_xyz_x[k] * cd_x[k] + g_y_0_xyz_xx[k];

                g_y_0_xxyz_y[k] = -g_y_0_xyz_y[k] * cd_x[k] + g_y_0_xyz_xy[k];

                g_y_0_xxyz_z[k] = -g_y_0_xyz_z[k] * cd_x[k] + g_y_0_xyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 15);

            auto g_y_0_xxzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 16);

            auto g_y_0_xxzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xxzz_x, g_y_0_xxzz_y, g_y_0_xxzz_z, g_y_0_xzz_x, g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_y, g_y_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_x[k] = -g_y_0_xzz_x[k] * cd_x[k] + g_y_0_xzz_xx[k];

                g_y_0_xxzz_y[k] = -g_y_0_xzz_y[k] * cd_x[k] + g_y_0_xzz_xy[k];

                g_y_0_xxzz_z[k] = -g_y_0_xzz_z[k] * cd_x[k] + g_y_0_xzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 18);

            auto g_y_0_xyyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 19);

            auto g_y_0_xyyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xyyy_x, g_y_0_xyyy_y, g_y_0_xyyy_z, g_y_0_yyy_x, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_y, g_y_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_x[k] = -g_y_0_yyy_x[k] * cd_x[k] + g_y_0_yyy_xx[k];

                g_y_0_xyyy_y[k] = -g_y_0_yyy_y[k] * cd_x[k] + g_y_0_yyy_xy[k];

                g_y_0_xyyy_z[k] = -g_y_0_yyy_z[k] * cd_x[k] + g_y_0_yyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 21);

            auto g_y_0_xyyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 22);

            auto g_y_0_xyyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_y_0_xyyz_x, g_y_0_xyyz_y, g_y_0_xyyz_z, g_y_0_yyz_x, g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_y, g_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_x[k] = -g_y_0_yyz_x[k] * cd_x[k] + g_y_0_yyz_xx[k];

                g_y_0_xyyz_y[k] = -g_y_0_yyz_y[k] * cd_x[k] + g_y_0_yyz_xy[k];

                g_y_0_xyyz_z[k] = -g_y_0_yyz_z[k] * cd_x[k] + g_y_0_yyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 24);

            auto g_y_0_xyzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 25);

            auto g_y_0_xyzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_y_0_xyzz_x, g_y_0_xyzz_y, g_y_0_xyzz_z, g_y_0_yzz_x, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_y, g_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_x[k] = -g_y_0_yzz_x[k] * cd_x[k] + g_y_0_yzz_xx[k];

                g_y_0_xyzz_y[k] = -g_y_0_yzz_y[k] * cd_x[k] + g_y_0_yzz_xy[k];

                g_y_0_xyzz_z[k] = -g_y_0_yzz_z[k] * cd_x[k] + g_y_0_yzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 27);

            auto g_y_0_xzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 28);

            auto g_y_0_xzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xzzz_x, g_y_0_xzzz_y, g_y_0_xzzz_z, g_y_0_zzz_x, g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_y, g_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_x[k] = -g_y_0_zzz_x[k] * cd_x[k] + g_y_0_zzz_xx[k];

                g_y_0_xzzz_y[k] = -g_y_0_zzz_y[k] * cd_x[k] + g_y_0_zzz_xy[k];

                g_y_0_xzzz_z[k] = -g_y_0_zzz_z[k] * cd_x[k] + g_y_0_zzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 30);

            auto g_y_0_yyyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 31);

            auto g_y_0_yyyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_y, g_y_0_yyy_x, g_y_0_yyy_xy, g_y_0_yyy_y, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_z, g_y_0_yyyy_x, g_y_0_yyyy_y, g_y_0_yyyy_z, g_yyy_x, g_yyy_y, g_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_x[k] = -g_yyy_x[k] - g_y_0_yyy_x[k] * cd_y[k] + g_y_0_yyy_xy[k];

                g_y_0_yyyy_y[k] = -g_yyy_y[k] - g_y_0_yyy_y[k] * cd_y[k] + g_y_0_yyy_yy[k];

                g_y_0_yyyy_z[k] = -g_yyy_z[k] - g_y_0_yyy_z[k] * cd_y[k] + g_y_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 33);

            auto g_y_0_yyyz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 34);

            auto g_y_0_yyyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_z, g_y_0_yyy_x, g_y_0_yyy_xz, g_y_0_yyy_y, g_y_0_yyy_yz, g_y_0_yyy_z, g_y_0_yyy_zz, g_y_0_yyyz_x, g_y_0_yyyz_y, g_y_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_x[k] = -g_y_0_yyy_x[k] * cd_z[k] + g_y_0_yyy_xz[k];

                g_y_0_yyyz_y[k] = -g_y_0_yyy_y[k] * cd_z[k] + g_y_0_yyy_yz[k];

                g_y_0_yyyz_z[k] = -g_y_0_yyy_z[k] * cd_z[k] + g_y_0_yyy_zz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 36);

            auto g_y_0_yyzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 37);

            auto g_y_0_yyzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_z, g_y_0_yyz_x, g_y_0_yyz_xz, g_y_0_yyz_y, g_y_0_yyz_yz, g_y_0_yyz_z, g_y_0_yyz_zz, g_y_0_yyzz_x, g_y_0_yyzz_y, g_y_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_x[k] = -g_y_0_yyz_x[k] * cd_z[k] + g_y_0_yyz_xz[k];

                g_y_0_yyzz_y[k] = -g_y_0_yyz_y[k] * cd_z[k] + g_y_0_yyz_yz[k];

                g_y_0_yyzz_z[k] = -g_y_0_yyz_z[k] * cd_z[k] + g_y_0_yyz_zz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 39);

            auto g_y_0_yzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 40);

            auto g_y_0_yzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_z, g_y_0_yzz_x, g_y_0_yzz_xz, g_y_0_yzz_y, g_y_0_yzz_yz, g_y_0_yzz_z, g_y_0_yzz_zz, g_y_0_yzzz_x, g_y_0_yzzz_y, g_y_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_x[k] = -g_y_0_yzz_x[k] * cd_z[k] + g_y_0_yzz_xz[k];

                g_y_0_yzzz_y[k] = -g_y_0_yzz_y[k] * cd_z[k] + g_y_0_yzz_yz[k];

                g_y_0_yzzz_z[k] = -g_y_0_yzz_z[k] * cd_z[k] + g_y_0_yzz_zz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 42);

            auto g_y_0_zzzz_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 43);

            auto g_y_0_zzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_y_0_zzz_x, g_y_0_zzz_xz, g_y_0_zzz_y, g_y_0_zzz_yz, g_y_0_zzz_z, g_y_0_zzz_zz, g_y_0_zzzz_x, g_y_0_zzzz_y, g_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_x[k] = -g_y_0_zzz_x[k] * cd_z[k] + g_y_0_zzz_xz[k];

                g_y_0_zzzz_y[k] = -g_y_0_zzz_y[k] * cd_z[k] + g_y_0_zzz_yz[k];

                g_y_0_zzzz_z[k] = -g_y_0_zzz_z[k] * cd_z[k] + g_y_0_zzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_z_0_xxxx_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 1);

            auto g_z_0_xxxx_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxx_x, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_y, g_z_0_xxx_z, g_z_0_xxxx_x, g_z_0_xxxx_y, g_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_x[k] = -g_z_0_xxx_x[k] * cd_x[k] + g_z_0_xxx_xx[k];

                g_z_0_xxxx_y[k] = -g_z_0_xxx_y[k] * cd_x[k] + g_z_0_xxx_xy[k];

                g_z_0_xxxx_z[k] = -g_z_0_xxx_z[k] * cd_x[k] + g_z_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 3);

            auto g_z_0_xxxy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 4);

            auto g_z_0_xxxy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxy_x, g_z_0_xxxy_y, g_z_0_xxxy_z, g_z_0_xxy_x, g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_y, g_z_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_x[k] = -g_z_0_xxy_x[k] * cd_x[k] + g_z_0_xxy_xx[k];

                g_z_0_xxxy_y[k] = -g_z_0_xxy_y[k] * cd_x[k] + g_z_0_xxy_xy[k];

                g_z_0_xxxy_z[k] = -g_z_0_xxy_z[k] * cd_x[k] + g_z_0_xxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_z_0_xxxz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 7);

            auto g_z_0_xxxz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xxxz_x, g_z_0_xxxz_y, g_z_0_xxxz_z, g_z_0_xxz_x, g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_y, g_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_x[k] = -g_z_0_xxz_x[k] * cd_x[k] + g_z_0_xxz_xx[k];

                g_z_0_xxxz_y[k] = -g_z_0_xxz_y[k] * cd_x[k] + g_z_0_xxz_xy[k];

                g_z_0_xxxz_z[k] = -g_z_0_xxz_z[k] * cd_x[k] + g_z_0_xxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 9);

            auto g_z_0_xxyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 10);

            auto g_z_0_xxyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxyy_x, g_z_0_xxyy_y, g_z_0_xxyy_z, g_z_0_xyy_x, g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_y, g_z_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_x[k] = -g_z_0_xyy_x[k] * cd_x[k] + g_z_0_xyy_xx[k];

                g_z_0_xxyy_y[k] = -g_z_0_xyy_y[k] * cd_x[k] + g_z_0_xyy_xy[k];

                g_z_0_xxyy_z[k] = -g_z_0_xyy_z[k] * cd_x[k] + g_z_0_xyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_z_0_xxyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 13);

            auto g_z_0_xxyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxyz_x, g_z_0_xxyz_y, g_z_0_xxyz_z, g_z_0_xyz_x, g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_y, g_z_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_x[k] = -g_z_0_xyz_x[k] * cd_x[k] + g_z_0_xyz_xx[k];

                g_z_0_xxyz_y[k] = -g_z_0_xyz_y[k] * cd_x[k] + g_z_0_xyz_xy[k];

                g_z_0_xxyz_z[k] = -g_z_0_xyz_z[k] * cd_x[k] + g_z_0_xyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 15);

            auto g_z_0_xxzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 16);

            auto g_z_0_xxzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xxzz_x, g_z_0_xxzz_y, g_z_0_xxzz_z, g_z_0_xzz_x, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_y, g_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_x[k] = -g_z_0_xzz_x[k] * cd_x[k] + g_z_0_xzz_xx[k];

                g_z_0_xxzz_y[k] = -g_z_0_xzz_y[k] * cd_x[k] + g_z_0_xzz_xy[k];

                g_z_0_xxzz_z[k] = -g_z_0_xzz_z[k] * cd_x[k] + g_z_0_xzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_z_0_xyyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 19);

            auto g_z_0_xyyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xyyy_x, g_z_0_xyyy_y, g_z_0_xyyy_z, g_z_0_yyy_x, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_y, g_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_x[k] = -g_z_0_yyy_x[k] * cd_x[k] + g_z_0_yyy_xx[k];

                g_z_0_xyyy_y[k] = -g_z_0_yyy_y[k] * cd_x[k] + g_z_0_yyy_xy[k];

                g_z_0_xyyy_z[k] = -g_z_0_yyy_z[k] * cd_x[k] + g_z_0_yyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 21);

            auto g_z_0_xyyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 22);

            auto g_z_0_xyyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_z_0_xyyz_x, g_z_0_xyyz_y, g_z_0_xyyz_z, g_z_0_yyz_x, g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_y, g_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_x[k] = -g_z_0_yyz_x[k] * cd_x[k] + g_z_0_yyz_xx[k];

                g_z_0_xyyz_y[k] = -g_z_0_yyz_y[k] * cd_x[k] + g_z_0_yyz_xy[k];

                g_z_0_xyyz_z[k] = -g_z_0_yyz_z[k] * cd_x[k] + g_z_0_yyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_z_0_xyzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 25);

            auto g_z_0_xyzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_z_0_xyzz_x, g_z_0_xyzz_y, g_z_0_xyzz_z, g_z_0_yzz_x, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_y, g_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_x[k] = -g_z_0_yzz_x[k] * cd_x[k] + g_z_0_yzz_xx[k];

                g_z_0_xyzz_y[k] = -g_z_0_yzz_y[k] * cd_x[k] + g_z_0_yzz_xy[k];

                g_z_0_xyzz_z[k] = -g_z_0_yzz_z[k] * cd_x[k] + g_z_0_yzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 27);

            auto g_z_0_xzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 28);

            auto g_z_0_xzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xzzz_x, g_z_0_xzzz_y, g_z_0_xzzz_z, g_z_0_zzz_x, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_y, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_x[k] = -g_z_0_zzz_x[k] * cd_x[k] + g_z_0_zzz_xx[k];

                g_z_0_xzzz_y[k] = -g_z_0_zzz_y[k] * cd_x[k] + g_z_0_zzz_xy[k];

                g_z_0_xzzz_z[k] = -g_z_0_zzz_z[k] * cd_x[k] + g_z_0_zzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_z_0_yyyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_z_0_yyyy_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_y, g_z_0_yyy_x, g_z_0_yyy_xy, g_z_0_yyy_y, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_z, g_z_0_yyyy_x, g_z_0_yyyy_y, g_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_x[k] = -g_z_0_yyy_x[k] * cd_y[k] + g_z_0_yyy_xy[k];

                g_z_0_yyyy_y[k] = -g_z_0_yyy_y[k] * cd_y[k] + g_z_0_yyy_yy[k];

                g_z_0_yyyy_z[k] = -g_z_0_yyy_z[k] * cd_y[k] + g_z_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 33);

            auto g_z_0_yyyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 34);

            auto g_z_0_yyyz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_y, g_z_0_yyyz_x, g_z_0_yyyz_y, g_z_0_yyyz_z, g_z_0_yyz_x, g_z_0_yyz_xy, g_z_0_yyz_y, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_x[k] = -g_z_0_yyz_x[k] * cd_y[k] + g_z_0_yyz_xy[k];

                g_z_0_yyyz_y[k] = -g_z_0_yyz_y[k] * cd_y[k] + g_z_0_yyz_yy[k];

                g_z_0_yyyz_z[k] = -g_z_0_yyz_z[k] * cd_y[k] + g_z_0_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_z_0_yyzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_z_0_yyzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_y, g_z_0_yyzz_x, g_z_0_yyzz_y, g_z_0_yyzz_z, g_z_0_yzz_x, g_z_0_yzz_xy, g_z_0_yzz_y, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_x[k] = -g_z_0_yzz_x[k] * cd_y[k] + g_z_0_yzz_xy[k];

                g_z_0_yyzz_y[k] = -g_z_0_yzz_y[k] * cd_y[k] + g_z_0_yzz_yy[k];

                g_z_0_yyzz_z[k] = -g_z_0_yzz_z[k] * cd_y[k] + g_z_0_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 39);

            auto g_z_0_yzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 40);

            auto g_z_0_yzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_z_0_yzzz_x, g_z_0_yzzz_y, g_z_0_yzzz_z, g_z_0_zzz_x, g_z_0_zzz_xy, g_z_0_zzz_y, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_x[k] = -g_z_0_zzz_x[k] * cd_y[k] + g_z_0_zzz_xy[k];

                g_z_0_yzzz_y[k] = -g_z_0_zzz_y[k] * cd_y[k] + g_z_0_zzz_yy[k];

                g_z_0_yzzz_z[k] = -g_z_0_zzz_z[k] * cd_y[k] + g_z_0_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_z_0_zzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_z_0_zzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_z_0_zzz_x, g_z_0_zzz_xz, g_z_0_zzz_y, g_z_0_zzz_yz, g_z_0_zzz_z, g_z_0_zzz_zz, g_z_0_zzzz_x, g_z_0_zzzz_y, g_z_0_zzzz_z, g_zzz_x, g_zzz_y, g_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_x[k] = -g_zzz_x[k] - g_z_0_zzz_x[k] * cd_z[k] + g_z_0_zzz_xz[k];

                g_z_0_zzzz_y[k] = -g_zzz_y[k] - g_z_0_zzz_y[k] * cd_z[k] + g_z_0_zzz_yz[k];

                g_z_0_zzzz_z[k] = -g_zzz_z[k] - g_z_0_zzz_z[k] * cd_z[k] + g_z_0_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

