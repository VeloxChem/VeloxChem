#include "ElectronRepulsionGeom0010ContrRecXXGF.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxgf(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgf,
                                            const size_t idx_xxff,
                                            const size_t idx_geom_10_xxff,
                                            const size_t idx_geom_10_xxfg,
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
            /// Set up components of auxilary buffer : SSFF

            const auto ff_off = idx_xxff + (i * bcomps + j) * 100;

            auto g_xxx_xxx = cbuffer.data(ff_off + 0);

            auto g_xxx_xxy = cbuffer.data(ff_off + 1);

            auto g_xxx_xxz = cbuffer.data(ff_off + 2);

            auto g_xxx_xyy = cbuffer.data(ff_off + 3);

            auto g_xxx_xyz = cbuffer.data(ff_off + 4);

            auto g_xxx_xzz = cbuffer.data(ff_off + 5);

            auto g_xxx_yyy = cbuffer.data(ff_off + 6);

            auto g_xxx_yyz = cbuffer.data(ff_off + 7);

            auto g_xxx_yzz = cbuffer.data(ff_off + 8);

            auto g_xxx_zzz = cbuffer.data(ff_off + 9);

            auto g_xxy_xxx = cbuffer.data(ff_off + 10);

            auto g_xxy_xxy = cbuffer.data(ff_off + 11);

            auto g_xxy_xxz = cbuffer.data(ff_off + 12);

            auto g_xxy_xyy = cbuffer.data(ff_off + 13);

            auto g_xxy_xyz = cbuffer.data(ff_off + 14);

            auto g_xxy_xzz = cbuffer.data(ff_off + 15);

            auto g_xxy_yyy = cbuffer.data(ff_off + 16);

            auto g_xxy_yyz = cbuffer.data(ff_off + 17);

            auto g_xxy_yzz = cbuffer.data(ff_off + 18);

            auto g_xxy_zzz = cbuffer.data(ff_off + 19);

            auto g_xxz_xxx = cbuffer.data(ff_off + 20);

            auto g_xxz_xxy = cbuffer.data(ff_off + 21);

            auto g_xxz_xxz = cbuffer.data(ff_off + 22);

            auto g_xxz_xyy = cbuffer.data(ff_off + 23);

            auto g_xxz_xyz = cbuffer.data(ff_off + 24);

            auto g_xxz_xzz = cbuffer.data(ff_off + 25);

            auto g_xxz_yyy = cbuffer.data(ff_off + 26);

            auto g_xxz_yyz = cbuffer.data(ff_off + 27);

            auto g_xxz_yzz = cbuffer.data(ff_off + 28);

            auto g_xxz_zzz = cbuffer.data(ff_off + 29);

            auto g_xyy_xxx = cbuffer.data(ff_off + 30);

            auto g_xyy_xxy = cbuffer.data(ff_off + 31);

            auto g_xyy_xxz = cbuffer.data(ff_off + 32);

            auto g_xyy_xyy = cbuffer.data(ff_off + 33);

            auto g_xyy_xyz = cbuffer.data(ff_off + 34);

            auto g_xyy_xzz = cbuffer.data(ff_off + 35);

            auto g_xyy_yyy = cbuffer.data(ff_off + 36);

            auto g_xyy_yyz = cbuffer.data(ff_off + 37);

            auto g_xyy_yzz = cbuffer.data(ff_off + 38);

            auto g_xyy_zzz = cbuffer.data(ff_off + 39);

            auto g_xyz_xxx = cbuffer.data(ff_off + 40);

            auto g_xyz_xxy = cbuffer.data(ff_off + 41);

            auto g_xyz_xxz = cbuffer.data(ff_off + 42);

            auto g_xyz_xyy = cbuffer.data(ff_off + 43);

            auto g_xyz_xyz = cbuffer.data(ff_off + 44);

            auto g_xyz_xzz = cbuffer.data(ff_off + 45);

            auto g_xyz_yyy = cbuffer.data(ff_off + 46);

            auto g_xyz_yyz = cbuffer.data(ff_off + 47);

            auto g_xyz_yzz = cbuffer.data(ff_off + 48);

            auto g_xyz_zzz = cbuffer.data(ff_off + 49);

            auto g_xzz_xxx = cbuffer.data(ff_off + 50);

            auto g_xzz_xxy = cbuffer.data(ff_off + 51);

            auto g_xzz_xxz = cbuffer.data(ff_off + 52);

            auto g_xzz_xyy = cbuffer.data(ff_off + 53);

            auto g_xzz_xyz = cbuffer.data(ff_off + 54);

            auto g_xzz_xzz = cbuffer.data(ff_off + 55);

            auto g_xzz_yyy = cbuffer.data(ff_off + 56);

            auto g_xzz_yyz = cbuffer.data(ff_off + 57);

            auto g_xzz_yzz = cbuffer.data(ff_off + 58);

            auto g_xzz_zzz = cbuffer.data(ff_off + 59);

            auto g_yyy_xxx = cbuffer.data(ff_off + 60);

            auto g_yyy_xxy = cbuffer.data(ff_off + 61);

            auto g_yyy_xxz = cbuffer.data(ff_off + 62);

            auto g_yyy_xyy = cbuffer.data(ff_off + 63);

            auto g_yyy_xyz = cbuffer.data(ff_off + 64);

            auto g_yyy_xzz = cbuffer.data(ff_off + 65);

            auto g_yyy_yyy = cbuffer.data(ff_off + 66);

            auto g_yyy_yyz = cbuffer.data(ff_off + 67);

            auto g_yyy_yzz = cbuffer.data(ff_off + 68);

            auto g_yyy_zzz = cbuffer.data(ff_off + 69);

            auto g_yyz_xxx = cbuffer.data(ff_off + 70);

            auto g_yyz_xxy = cbuffer.data(ff_off + 71);

            auto g_yyz_xxz = cbuffer.data(ff_off + 72);

            auto g_yyz_xyy = cbuffer.data(ff_off + 73);

            auto g_yyz_xyz = cbuffer.data(ff_off + 74);

            auto g_yyz_xzz = cbuffer.data(ff_off + 75);

            auto g_yyz_yyy = cbuffer.data(ff_off + 76);

            auto g_yyz_yyz = cbuffer.data(ff_off + 77);

            auto g_yyz_yzz = cbuffer.data(ff_off + 78);

            auto g_yyz_zzz = cbuffer.data(ff_off + 79);

            auto g_yzz_xxx = cbuffer.data(ff_off + 80);

            auto g_yzz_xxy = cbuffer.data(ff_off + 81);

            auto g_yzz_xxz = cbuffer.data(ff_off + 82);

            auto g_yzz_xyy = cbuffer.data(ff_off + 83);

            auto g_yzz_xyz = cbuffer.data(ff_off + 84);

            auto g_yzz_xzz = cbuffer.data(ff_off + 85);

            auto g_yzz_yyy = cbuffer.data(ff_off + 86);

            auto g_yzz_yyz = cbuffer.data(ff_off + 87);

            auto g_yzz_yzz = cbuffer.data(ff_off + 88);

            auto g_yzz_zzz = cbuffer.data(ff_off + 89);

            auto g_zzz_xxx = cbuffer.data(ff_off + 90);

            auto g_zzz_xxy = cbuffer.data(ff_off + 91);

            auto g_zzz_xxz = cbuffer.data(ff_off + 92);

            auto g_zzz_xyy = cbuffer.data(ff_off + 93);

            auto g_zzz_xyz = cbuffer.data(ff_off + 94);

            auto g_zzz_xzz = cbuffer.data(ff_off + 95);

            auto g_zzz_yyy = cbuffer.data(ff_off + 96);

            auto g_zzz_yyz = cbuffer.data(ff_off + 97);

            auto g_zzz_yzz = cbuffer.data(ff_off + 98);

            auto g_zzz_zzz = cbuffer.data(ff_off + 99);

            /// Set up components of auxilary buffer : SSFF

            const auto ff_geom_10_off = idx_geom_10_xxff + (i * bcomps + j) * 100;

            auto g_x_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_y_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 0);

            auto g_y_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 1);

            auto g_y_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 2);

            auto g_y_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 3);

            auto g_y_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 4);

            auto g_y_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 5);

            auto g_y_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 6);

            auto g_y_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 7);

            auto g_y_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 8);

            auto g_y_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 9);

            auto g_y_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 10);

            auto g_y_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 11);

            auto g_y_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 12);

            auto g_y_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 13);

            auto g_y_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 14);

            auto g_y_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 15);

            auto g_y_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 16);

            auto g_y_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 17);

            auto g_y_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 18);

            auto g_y_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 19);

            auto g_y_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 20);

            auto g_y_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 21);

            auto g_y_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 22);

            auto g_y_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 23);

            auto g_y_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 24);

            auto g_y_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 25);

            auto g_y_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 26);

            auto g_y_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 27);

            auto g_y_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 28);

            auto g_y_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 29);

            auto g_y_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 30);

            auto g_y_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 31);

            auto g_y_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 32);

            auto g_y_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 33);

            auto g_y_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 34);

            auto g_y_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 35);

            auto g_y_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 36);

            auto g_y_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 37);

            auto g_y_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 38);

            auto g_y_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 39);

            auto g_y_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 40);

            auto g_y_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 41);

            auto g_y_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 42);

            auto g_y_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 43);

            auto g_y_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 44);

            auto g_y_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 45);

            auto g_y_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 46);

            auto g_y_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 47);

            auto g_y_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 48);

            auto g_y_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 49);

            auto g_y_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 50);

            auto g_y_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 51);

            auto g_y_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 52);

            auto g_y_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 53);

            auto g_y_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 54);

            auto g_y_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 55);

            auto g_y_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 56);

            auto g_y_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 57);

            auto g_y_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 58);

            auto g_y_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 59);

            auto g_y_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 60);

            auto g_y_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 61);

            auto g_y_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 62);

            auto g_y_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 63);

            auto g_y_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 64);

            auto g_y_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 65);

            auto g_y_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 66);

            auto g_y_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 67);

            auto g_y_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 68);

            auto g_y_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 69);

            auto g_y_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 70);

            auto g_y_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 71);

            auto g_y_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 72);

            auto g_y_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 73);

            auto g_y_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 74);

            auto g_y_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 75);

            auto g_y_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 76);

            auto g_y_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 77);

            auto g_y_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 78);

            auto g_y_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 79);

            auto g_y_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 80);

            auto g_y_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 81);

            auto g_y_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 82);

            auto g_y_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 83);

            auto g_y_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 84);

            auto g_y_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 85);

            auto g_y_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 86);

            auto g_y_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 87);

            auto g_y_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 88);

            auto g_y_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 89);

            auto g_y_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 90);

            auto g_y_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 91);

            auto g_y_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 92);

            auto g_y_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 93);

            auto g_y_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 94);

            auto g_y_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 95);

            auto g_y_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 96);

            auto g_y_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 97);

            auto g_y_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 98);

            auto g_y_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 100 * acomps * bcomps + 99);

            auto g_z_0_xxx_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 0);

            auto g_z_0_xxx_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 1);

            auto g_z_0_xxx_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 2);

            auto g_z_0_xxx_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 3);

            auto g_z_0_xxx_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 4);

            auto g_z_0_xxx_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 5);

            auto g_z_0_xxx_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 6);

            auto g_z_0_xxx_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 7);

            auto g_z_0_xxx_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 8);

            auto g_z_0_xxx_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 9);

            auto g_z_0_xxy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 10);

            auto g_z_0_xxy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 11);

            auto g_z_0_xxy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 12);

            auto g_z_0_xxy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 13);

            auto g_z_0_xxy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 14);

            auto g_z_0_xxy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 15);

            auto g_z_0_xxy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 16);

            auto g_z_0_xxy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 17);

            auto g_z_0_xxy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 18);

            auto g_z_0_xxy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 19);

            auto g_z_0_xxz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 20);

            auto g_z_0_xxz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 21);

            auto g_z_0_xxz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 22);

            auto g_z_0_xxz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 23);

            auto g_z_0_xxz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 24);

            auto g_z_0_xxz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 25);

            auto g_z_0_xxz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 26);

            auto g_z_0_xxz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 27);

            auto g_z_0_xxz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 28);

            auto g_z_0_xxz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 29);

            auto g_z_0_xyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 30);

            auto g_z_0_xyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 31);

            auto g_z_0_xyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 32);

            auto g_z_0_xyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 33);

            auto g_z_0_xyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 34);

            auto g_z_0_xyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 35);

            auto g_z_0_xyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 36);

            auto g_z_0_xyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 37);

            auto g_z_0_xyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 38);

            auto g_z_0_xyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 39);

            auto g_z_0_xyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 40);

            auto g_z_0_xyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 41);

            auto g_z_0_xyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 42);

            auto g_z_0_xyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 43);

            auto g_z_0_xyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 44);

            auto g_z_0_xyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 45);

            auto g_z_0_xyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 46);

            auto g_z_0_xyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 47);

            auto g_z_0_xyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 48);

            auto g_z_0_xyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 49);

            auto g_z_0_xzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 50);

            auto g_z_0_xzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 51);

            auto g_z_0_xzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 52);

            auto g_z_0_xzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 53);

            auto g_z_0_xzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 54);

            auto g_z_0_xzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 55);

            auto g_z_0_xzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 56);

            auto g_z_0_xzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 57);

            auto g_z_0_xzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 58);

            auto g_z_0_xzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 59);

            auto g_z_0_yyy_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 60);

            auto g_z_0_yyy_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 61);

            auto g_z_0_yyy_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 62);

            auto g_z_0_yyy_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 63);

            auto g_z_0_yyy_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 64);

            auto g_z_0_yyy_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 65);

            auto g_z_0_yyy_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 66);

            auto g_z_0_yyy_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 67);

            auto g_z_0_yyy_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 68);

            auto g_z_0_yyy_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 69);

            auto g_z_0_yyz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 70);

            auto g_z_0_yyz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 71);

            auto g_z_0_yyz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 72);

            auto g_z_0_yyz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 73);

            auto g_z_0_yyz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 74);

            auto g_z_0_yyz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 75);

            auto g_z_0_yyz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 76);

            auto g_z_0_yyz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 77);

            auto g_z_0_yyz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 78);

            auto g_z_0_yyz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 79);

            auto g_z_0_yzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 80);

            auto g_z_0_yzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 81);

            auto g_z_0_yzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 82);

            auto g_z_0_yzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 83);

            auto g_z_0_yzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 84);

            auto g_z_0_yzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 85);

            auto g_z_0_yzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 86);

            auto g_z_0_yzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 87);

            auto g_z_0_yzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 88);

            auto g_z_0_yzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 89);

            auto g_z_0_zzz_xxx = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 90);

            auto g_z_0_zzz_xxy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 91);

            auto g_z_0_zzz_xxz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 92);

            auto g_z_0_zzz_xyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 93);

            auto g_z_0_zzz_xyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 94);

            auto g_z_0_zzz_xzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 95);

            auto g_z_0_zzz_yyy = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 96);

            auto g_z_0_zzz_yyz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 97);

            auto g_z_0_zzz_yzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 98);

            auto g_z_0_zzz_zzz = cbuffer.data(ff_geom_10_off + 200 * acomps * bcomps + 99);

            /// Set up components of auxilary buffer : SSFG

            const auto fg_geom_10_off = idx_geom_10_xxfg + (i * bcomps + j) * 150;

            auto g_x_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_y_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 0);

            auto g_y_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 1);

            auto g_y_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 2);

            auto g_y_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 3);

            auto g_y_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 4);

            auto g_y_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 5);

            auto g_y_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 6);

            auto g_y_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 7);

            auto g_y_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 8);

            auto g_y_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 9);

            auto g_y_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 10);

            auto g_y_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 11);

            auto g_y_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 12);

            auto g_y_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 13);

            auto g_y_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 14);

            auto g_y_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 15);

            auto g_y_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 16);

            auto g_y_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 17);

            auto g_y_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 18);

            auto g_y_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 19);

            auto g_y_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 20);

            auto g_y_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 21);

            auto g_y_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 22);

            auto g_y_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 23);

            auto g_y_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 24);

            auto g_y_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 25);

            auto g_y_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 26);

            auto g_y_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 27);

            auto g_y_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 28);

            auto g_y_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 29);

            auto g_y_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 30);

            auto g_y_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 31);

            auto g_y_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 32);

            auto g_y_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 33);

            auto g_y_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 34);

            auto g_y_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 35);

            auto g_y_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 36);

            auto g_y_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 37);

            auto g_y_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 38);

            auto g_y_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 39);

            auto g_y_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 40);

            auto g_y_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 41);

            auto g_y_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 42);

            auto g_y_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 43);

            auto g_y_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 44);

            auto g_y_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 45);

            auto g_y_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 46);

            auto g_y_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 47);

            auto g_y_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 48);

            auto g_y_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 49);

            auto g_y_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 50);

            auto g_y_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 51);

            auto g_y_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 52);

            auto g_y_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 53);

            auto g_y_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 54);

            auto g_y_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 55);

            auto g_y_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 56);

            auto g_y_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 57);

            auto g_y_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 58);

            auto g_y_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 59);

            auto g_y_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 60);

            auto g_y_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 61);

            auto g_y_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 62);

            auto g_y_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 63);

            auto g_y_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 64);

            auto g_y_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 65);

            auto g_y_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 66);

            auto g_y_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 67);

            auto g_y_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 68);

            auto g_y_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 69);

            auto g_y_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 70);

            auto g_y_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 71);

            auto g_y_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 72);

            auto g_y_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 73);

            auto g_y_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 74);

            auto g_y_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 75);

            auto g_y_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 76);

            auto g_y_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 77);

            auto g_y_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 78);

            auto g_y_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 79);

            auto g_y_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 80);

            auto g_y_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 81);

            auto g_y_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 82);

            auto g_y_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 83);

            auto g_y_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 84);

            auto g_y_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 85);

            auto g_y_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 86);

            auto g_y_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 87);

            auto g_y_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 88);

            auto g_y_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 89);

            auto g_y_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 90);

            auto g_y_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 91);

            auto g_y_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 92);

            auto g_y_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 93);

            auto g_y_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 94);

            auto g_y_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 95);

            auto g_y_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 96);

            auto g_y_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 97);

            auto g_y_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 98);

            auto g_y_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 99);

            auto g_y_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 100);

            auto g_y_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 101);

            auto g_y_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 102);

            auto g_y_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 103);

            auto g_y_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 104);

            auto g_y_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 105);

            auto g_y_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 106);

            auto g_y_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 107);

            auto g_y_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 108);

            auto g_y_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 109);

            auto g_y_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 110);

            auto g_y_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 111);

            auto g_y_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 112);

            auto g_y_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 113);

            auto g_y_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 114);

            auto g_y_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 115);

            auto g_y_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 116);

            auto g_y_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 117);

            auto g_y_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 118);

            auto g_y_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 119);

            auto g_y_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 120);

            auto g_y_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 121);

            auto g_y_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 122);

            auto g_y_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 123);

            auto g_y_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 124);

            auto g_y_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 125);

            auto g_y_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 126);

            auto g_y_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 127);

            auto g_y_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 128);

            auto g_y_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 129);

            auto g_y_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 130);

            auto g_y_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 131);

            auto g_y_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 132);

            auto g_y_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 133);

            auto g_y_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 134);

            auto g_y_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 135);

            auto g_y_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 136);

            auto g_y_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 137);

            auto g_y_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 138);

            auto g_y_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 139);

            auto g_y_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 140);

            auto g_y_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 141);

            auto g_y_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 142);

            auto g_y_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 143);

            auto g_y_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 144);

            auto g_y_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 145);

            auto g_y_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 146);

            auto g_y_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 147);

            auto g_y_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 148);

            auto g_y_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps * bcomps + 149);

            auto g_z_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 0);

            auto g_z_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 1);

            auto g_z_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 2);

            auto g_z_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 3);

            auto g_z_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 4);

            auto g_z_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 5);

            auto g_z_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 6);

            auto g_z_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 7);

            auto g_z_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 8);

            auto g_z_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 9);

            auto g_z_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 10);

            auto g_z_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 11);

            auto g_z_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 12);

            auto g_z_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 13);

            auto g_z_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 14);

            auto g_z_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 15);

            auto g_z_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 16);

            auto g_z_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 17);

            auto g_z_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 18);

            auto g_z_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 19);

            auto g_z_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 20);

            auto g_z_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 21);

            auto g_z_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 22);

            auto g_z_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 23);

            auto g_z_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 24);

            auto g_z_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 25);

            auto g_z_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 26);

            auto g_z_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 27);

            auto g_z_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 28);

            auto g_z_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 29);

            auto g_z_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 30);

            auto g_z_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 31);

            auto g_z_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 32);

            auto g_z_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 33);

            auto g_z_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 34);

            auto g_z_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 35);

            auto g_z_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 36);

            auto g_z_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 37);

            auto g_z_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 38);

            auto g_z_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 39);

            auto g_z_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 40);

            auto g_z_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 41);

            auto g_z_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 42);

            auto g_z_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 43);

            auto g_z_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 44);

            auto g_z_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 45);

            auto g_z_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 46);

            auto g_z_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 47);

            auto g_z_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 48);

            auto g_z_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 49);

            auto g_z_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 50);

            auto g_z_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 51);

            auto g_z_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 52);

            auto g_z_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 53);

            auto g_z_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 54);

            auto g_z_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 55);

            auto g_z_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 56);

            auto g_z_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 57);

            auto g_z_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 58);

            auto g_z_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 59);

            auto g_z_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 60);

            auto g_z_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 61);

            auto g_z_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 62);

            auto g_z_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 63);

            auto g_z_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 64);

            auto g_z_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 65);

            auto g_z_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 66);

            auto g_z_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 67);

            auto g_z_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 68);

            auto g_z_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 69);

            auto g_z_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 70);

            auto g_z_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 71);

            auto g_z_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 72);

            auto g_z_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 73);

            auto g_z_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 74);

            auto g_z_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 75);

            auto g_z_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 76);

            auto g_z_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 77);

            auto g_z_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 78);

            auto g_z_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 79);

            auto g_z_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 80);

            auto g_z_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 81);

            auto g_z_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 82);

            auto g_z_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 83);

            auto g_z_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 84);

            auto g_z_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 85);

            auto g_z_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 86);

            auto g_z_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 87);

            auto g_z_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 88);

            auto g_z_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 89);

            auto g_z_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 90);

            auto g_z_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 91);

            auto g_z_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 92);

            auto g_z_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 93);

            auto g_z_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 94);

            auto g_z_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 95);

            auto g_z_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 96);

            auto g_z_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 97);

            auto g_z_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 98);

            auto g_z_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 99);

            auto g_z_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 100);

            auto g_z_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 101);

            auto g_z_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 102);

            auto g_z_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 103);

            auto g_z_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 104);

            auto g_z_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 105);

            auto g_z_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 106);

            auto g_z_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 107);

            auto g_z_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 108);

            auto g_z_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 109);

            auto g_z_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 110);

            auto g_z_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 111);

            auto g_z_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 112);

            auto g_z_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 113);

            auto g_z_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 114);

            auto g_z_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 115);

            auto g_z_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 116);

            auto g_z_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 117);

            auto g_z_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 118);

            auto g_z_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 119);

            auto g_z_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 120);

            auto g_z_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 121);

            auto g_z_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 122);

            auto g_z_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 123);

            auto g_z_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 124);

            auto g_z_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 125);

            auto g_z_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 126);

            auto g_z_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 127);

            auto g_z_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 128);

            auto g_z_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 129);

            auto g_z_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 130);

            auto g_z_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 131);

            auto g_z_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 132);

            auto g_z_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 133);

            auto g_z_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 134);

            auto g_z_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 135);

            auto g_z_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 136);

            auto g_z_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 137);

            auto g_z_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 138);

            auto g_z_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 139);

            auto g_z_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 140);

            auto g_z_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 141);

            auto g_z_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 142);

            auto g_z_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 143);

            auto g_z_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 144);

            auto g_z_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 145);

            auto g_z_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 146);

            auto g_z_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 147);

            auto g_z_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 148);

            auto g_z_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps * bcomps + 149);

            /// set up bra offset for contr_buffer_xxgf

            const auto gf_geom_10_off = idx_geom_10_xxgf + (i * bcomps + j) * 150;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_x_0_xxx_xxx, g_x_0_xxx_xxxx, g_x_0_xxx_xxxy, g_x_0_xxx_xxxz, g_x_0_xxx_xxy, g_x_0_xxx_xxyy, g_x_0_xxx_xxyz, g_x_0_xxx_xxz, g_x_0_xxx_xxzz, g_x_0_xxx_xyy, g_x_0_xxx_xyyy, g_x_0_xxx_xyyz, g_x_0_xxx_xyz, g_x_0_xxx_xyzz, g_x_0_xxx_xzz, g_x_0_xxx_xzzz, g_x_0_xxx_yyy, g_x_0_xxx_yyz, g_x_0_xxx_yzz, g_x_0_xxx_zzz, g_x_0_xxxx_xxx, g_x_0_xxxx_xxy, g_x_0_xxxx_xxz, g_x_0_xxxx_xyy, g_x_0_xxxx_xyz, g_x_0_xxxx_xzz, g_x_0_xxxx_yyy, g_x_0_xxxx_yyz, g_x_0_xxxx_yzz, g_x_0_xxxx_zzz, g_xxx_xxx, g_xxx_xxy, g_xxx_xxz, g_xxx_xyy, g_xxx_xyz, g_xxx_xzz, g_xxx_yyy, g_xxx_yyz, g_xxx_yzz, g_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_xxx[k] = -g_xxx_xxx[k] - g_x_0_xxx_xxx[k] * cd_x[k] + g_x_0_xxx_xxxx[k];

                g_x_0_xxxx_xxy[k] = -g_xxx_xxy[k] - g_x_0_xxx_xxy[k] * cd_x[k] + g_x_0_xxx_xxxy[k];

                g_x_0_xxxx_xxz[k] = -g_xxx_xxz[k] - g_x_0_xxx_xxz[k] * cd_x[k] + g_x_0_xxx_xxxz[k];

                g_x_0_xxxx_xyy[k] = -g_xxx_xyy[k] - g_x_0_xxx_xyy[k] * cd_x[k] + g_x_0_xxx_xxyy[k];

                g_x_0_xxxx_xyz[k] = -g_xxx_xyz[k] - g_x_0_xxx_xyz[k] * cd_x[k] + g_x_0_xxx_xxyz[k];

                g_x_0_xxxx_xzz[k] = -g_xxx_xzz[k] - g_x_0_xxx_xzz[k] * cd_x[k] + g_x_0_xxx_xxzz[k];

                g_x_0_xxxx_yyy[k] = -g_xxx_yyy[k] - g_x_0_xxx_yyy[k] * cd_x[k] + g_x_0_xxx_xyyy[k];

                g_x_0_xxxx_yyz[k] = -g_xxx_yyz[k] - g_x_0_xxx_yyz[k] * cd_x[k] + g_x_0_xxx_xyyz[k];

                g_x_0_xxxx_yzz[k] = -g_xxx_yzz[k] - g_x_0_xxx_yzz[k] * cd_x[k] + g_x_0_xxx_xyzz[k];

                g_x_0_xxxx_zzz[k] = -g_xxx_zzz[k] - g_x_0_xxx_zzz[k] * cd_x[k] + g_x_0_xxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xxx_xxx, g_x_0_xxx_xxxy, g_x_0_xxx_xxy, g_x_0_xxx_xxyy, g_x_0_xxx_xxyz, g_x_0_xxx_xxz, g_x_0_xxx_xyy, g_x_0_xxx_xyyy, g_x_0_xxx_xyyz, g_x_0_xxx_xyz, g_x_0_xxx_xyzz, g_x_0_xxx_xzz, g_x_0_xxx_yyy, g_x_0_xxx_yyyy, g_x_0_xxx_yyyz, g_x_0_xxx_yyz, g_x_0_xxx_yyzz, g_x_0_xxx_yzz, g_x_0_xxx_yzzz, g_x_0_xxx_zzz, g_x_0_xxxy_xxx, g_x_0_xxxy_xxy, g_x_0_xxxy_xxz, g_x_0_xxxy_xyy, g_x_0_xxxy_xyz, g_x_0_xxxy_xzz, g_x_0_xxxy_yyy, g_x_0_xxxy_yyz, g_x_0_xxxy_yzz, g_x_0_xxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_xxx[k] = -g_x_0_xxx_xxx[k] * cd_y[k] + g_x_0_xxx_xxxy[k];

                g_x_0_xxxy_xxy[k] = -g_x_0_xxx_xxy[k] * cd_y[k] + g_x_0_xxx_xxyy[k];

                g_x_0_xxxy_xxz[k] = -g_x_0_xxx_xxz[k] * cd_y[k] + g_x_0_xxx_xxyz[k];

                g_x_0_xxxy_xyy[k] = -g_x_0_xxx_xyy[k] * cd_y[k] + g_x_0_xxx_xyyy[k];

                g_x_0_xxxy_xyz[k] = -g_x_0_xxx_xyz[k] * cd_y[k] + g_x_0_xxx_xyyz[k];

                g_x_0_xxxy_xzz[k] = -g_x_0_xxx_xzz[k] * cd_y[k] + g_x_0_xxx_xyzz[k];

                g_x_0_xxxy_yyy[k] = -g_x_0_xxx_yyy[k] * cd_y[k] + g_x_0_xxx_yyyy[k];

                g_x_0_xxxy_yyz[k] = -g_x_0_xxx_yyz[k] * cd_y[k] + g_x_0_xxx_yyyz[k];

                g_x_0_xxxy_yzz[k] = -g_x_0_xxx_yzz[k] * cd_y[k] + g_x_0_xxx_yyzz[k];

                g_x_0_xxxy_zzz[k] = -g_x_0_xxx_zzz[k] * cd_y[k] + g_x_0_xxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_xxx_xxx, g_x_0_xxx_xxxz, g_x_0_xxx_xxy, g_x_0_xxx_xxyz, g_x_0_xxx_xxz, g_x_0_xxx_xxzz, g_x_0_xxx_xyy, g_x_0_xxx_xyyz, g_x_0_xxx_xyz, g_x_0_xxx_xyzz, g_x_0_xxx_xzz, g_x_0_xxx_xzzz, g_x_0_xxx_yyy, g_x_0_xxx_yyyz, g_x_0_xxx_yyz, g_x_0_xxx_yyzz, g_x_0_xxx_yzz, g_x_0_xxx_yzzz, g_x_0_xxx_zzz, g_x_0_xxx_zzzz, g_x_0_xxxz_xxx, g_x_0_xxxz_xxy, g_x_0_xxxz_xxz, g_x_0_xxxz_xyy, g_x_0_xxxz_xyz, g_x_0_xxxz_xzz, g_x_0_xxxz_yyy, g_x_0_xxxz_yyz, g_x_0_xxxz_yzz, g_x_0_xxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_xxx[k] = -g_x_0_xxx_xxx[k] * cd_z[k] + g_x_0_xxx_xxxz[k];

                g_x_0_xxxz_xxy[k] = -g_x_0_xxx_xxy[k] * cd_z[k] + g_x_0_xxx_xxyz[k];

                g_x_0_xxxz_xxz[k] = -g_x_0_xxx_xxz[k] * cd_z[k] + g_x_0_xxx_xxzz[k];

                g_x_0_xxxz_xyy[k] = -g_x_0_xxx_xyy[k] * cd_z[k] + g_x_0_xxx_xyyz[k];

                g_x_0_xxxz_xyz[k] = -g_x_0_xxx_xyz[k] * cd_z[k] + g_x_0_xxx_xyzz[k];

                g_x_0_xxxz_xzz[k] = -g_x_0_xxx_xzz[k] * cd_z[k] + g_x_0_xxx_xzzz[k];

                g_x_0_xxxz_yyy[k] = -g_x_0_xxx_yyy[k] * cd_z[k] + g_x_0_xxx_yyyz[k];

                g_x_0_xxxz_yyz[k] = -g_x_0_xxx_yyz[k] * cd_z[k] + g_x_0_xxx_yyzz[k];

                g_x_0_xxxz_yzz[k] = -g_x_0_xxx_yzz[k] * cd_z[k] + g_x_0_xxx_yzzz[k];

                g_x_0_xxxz_zzz[k] = -g_x_0_xxx_zzz[k] * cd_z[k] + g_x_0_xxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xxy_xxx, g_x_0_xxy_xxxy, g_x_0_xxy_xxy, g_x_0_xxy_xxyy, g_x_0_xxy_xxyz, g_x_0_xxy_xxz, g_x_0_xxy_xyy, g_x_0_xxy_xyyy, g_x_0_xxy_xyyz, g_x_0_xxy_xyz, g_x_0_xxy_xyzz, g_x_0_xxy_xzz, g_x_0_xxy_yyy, g_x_0_xxy_yyyy, g_x_0_xxy_yyyz, g_x_0_xxy_yyz, g_x_0_xxy_yyzz, g_x_0_xxy_yzz, g_x_0_xxy_yzzz, g_x_0_xxy_zzz, g_x_0_xxyy_xxx, g_x_0_xxyy_xxy, g_x_0_xxyy_xxz, g_x_0_xxyy_xyy, g_x_0_xxyy_xyz, g_x_0_xxyy_xzz, g_x_0_xxyy_yyy, g_x_0_xxyy_yyz, g_x_0_xxyy_yzz, g_x_0_xxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_xxx[k] = -g_x_0_xxy_xxx[k] * cd_y[k] + g_x_0_xxy_xxxy[k];

                g_x_0_xxyy_xxy[k] = -g_x_0_xxy_xxy[k] * cd_y[k] + g_x_0_xxy_xxyy[k];

                g_x_0_xxyy_xxz[k] = -g_x_0_xxy_xxz[k] * cd_y[k] + g_x_0_xxy_xxyz[k];

                g_x_0_xxyy_xyy[k] = -g_x_0_xxy_xyy[k] * cd_y[k] + g_x_0_xxy_xyyy[k];

                g_x_0_xxyy_xyz[k] = -g_x_0_xxy_xyz[k] * cd_y[k] + g_x_0_xxy_xyyz[k];

                g_x_0_xxyy_xzz[k] = -g_x_0_xxy_xzz[k] * cd_y[k] + g_x_0_xxy_xyzz[k];

                g_x_0_xxyy_yyy[k] = -g_x_0_xxy_yyy[k] * cd_y[k] + g_x_0_xxy_yyyy[k];

                g_x_0_xxyy_yyz[k] = -g_x_0_xxy_yyz[k] * cd_y[k] + g_x_0_xxy_yyyz[k];

                g_x_0_xxyy_yzz[k] = -g_x_0_xxy_yzz[k] * cd_y[k] + g_x_0_xxy_yyzz[k];

                g_x_0_xxyy_zzz[k] = -g_x_0_xxy_zzz[k] * cd_y[k] + g_x_0_xxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xxyz_xxx, g_x_0_xxyz_xxy, g_x_0_xxyz_xxz, g_x_0_xxyz_xyy, g_x_0_xxyz_xyz, g_x_0_xxyz_xzz, g_x_0_xxyz_yyy, g_x_0_xxyz_yyz, g_x_0_xxyz_yzz, g_x_0_xxyz_zzz, g_x_0_xxz_xxx, g_x_0_xxz_xxxy, g_x_0_xxz_xxy, g_x_0_xxz_xxyy, g_x_0_xxz_xxyz, g_x_0_xxz_xxz, g_x_0_xxz_xyy, g_x_0_xxz_xyyy, g_x_0_xxz_xyyz, g_x_0_xxz_xyz, g_x_0_xxz_xyzz, g_x_0_xxz_xzz, g_x_0_xxz_yyy, g_x_0_xxz_yyyy, g_x_0_xxz_yyyz, g_x_0_xxz_yyz, g_x_0_xxz_yyzz, g_x_0_xxz_yzz, g_x_0_xxz_yzzz, g_x_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_xxx[k] = -g_x_0_xxz_xxx[k] * cd_y[k] + g_x_0_xxz_xxxy[k];

                g_x_0_xxyz_xxy[k] = -g_x_0_xxz_xxy[k] * cd_y[k] + g_x_0_xxz_xxyy[k];

                g_x_0_xxyz_xxz[k] = -g_x_0_xxz_xxz[k] * cd_y[k] + g_x_0_xxz_xxyz[k];

                g_x_0_xxyz_xyy[k] = -g_x_0_xxz_xyy[k] * cd_y[k] + g_x_0_xxz_xyyy[k];

                g_x_0_xxyz_xyz[k] = -g_x_0_xxz_xyz[k] * cd_y[k] + g_x_0_xxz_xyyz[k];

                g_x_0_xxyz_xzz[k] = -g_x_0_xxz_xzz[k] * cd_y[k] + g_x_0_xxz_xyzz[k];

                g_x_0_xxyz_yyy[k] = -g_x_0_xxz_yyy[k] * cd_y[k] + g_x_0_xxz_yyyy[k];

                g_x_0_xxyz_yyz[k] = -g_x_0_xxz_yyz[k] * cd_y[k] + g_x_0_xxz_yyyz[k];

                g_x_0_xxyz_yzz[k] = -g_x_0_xxz_yzz[k] * cd_y[k] + g_x_0_xxz_yyzz[k];

                g_x_0_xxyz_zzz[k] = -g_x_0_xxz_zzz[k] * cd_y[k] + g_x_0_xxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_xxz_xxx, g_x_0_xxz_xxxz, g_x_0_xxz_xxy, g_x_0_xxz_xxyz, g_x_0_xxz_xxz, g_x_0_xxz_xxzz, g_x_0_xxz_xyy, g_x_0_xxz_xyyz, g_x_0_xxz_xyz, g_x_0_xxz_xyzz, g_x_0_xxz_xzz, g_x_0_xxz_xzzz, g_x_0_xxz_yyy, g_x_0_xxz_yyyz, g_x_0_xxz_yyz, g_x_0_xxz_yyzz, g_x_0_xxz_yzz, g_x_0_xxz_yzzz, g_x_0_xxz_zzz, g_x_0_xxz_zzzz, g_x_0_xxzz_xxx, g_x_0_xxzz_xxy, g_x_0_xxzz_xxz, g_x_0_xxzz_xyy, g_x_0_xxzz_xyz, g_x_0_xxzz_xzz, g_x_0_xxzz_yyy, g_x_0_xxzz_yyz, g_x_0_xxzz_yzz, g_x_0_xxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_xxx[k] = -g_x_0_xxz_xxx[k] * cd_z[k] + g_x_0_xxz_xxxz[k];

                g_x_0_xxzz_xxy[k] = -g_x_0_xxz_xxy[k] * cd_z[k] + g_x_0_xxz_xxyz[k];

                g_x_0_xxzz_xxz[k] = -g_x_0_xxz_xxz[k] * cd_z[k] + g_x_0_xxz_xxzz[k];

                g_x_0_xxzz_xyy[k] = -g_x_0_xxz_xyy[k] * cd_z[k] + g_x_0_xxz_xyyz[k];

                g_x_0_xxzz_xyz[k] = -g_x_0_xxz_xyz[k] * cd_z[k] + g_x_0_xxz_xyzz[k];

                g_x_0_xxzz_xzz[k] = -g_x_0_xxz_xzz[k] * cd_z[k] + g_x_0_xxz_xzzz[k];

                g_x_0_xxzz_yyy[k] = -g_x_0_xxz_yyy[k] * cd_z[k] + g_x_0_xxz_yyyz[k];

                g_x_0_xxzz_yyz[k] = -g_x_0_xxz_yyz[k] * cd_z[k] + g_x_0_xxz_yyzz[k];

                g_x_0_xxzz_yzz[k] = -g_x_0_xxz_yzz[k] * cd_z[k] + g_x_0_xxz_yzzz[k];

                g_x_0_xxzz_zzz[k] = -g_x_0_xxz_zzz[k] * cd_z[k] + g_x_0_xxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xyy_xxx, g_x_0_xyy_xxxy, g_x_0_xyy_xxy, g_x_0_xyy_xxyy, g_x_0_xyy_xxyz, g_x_0_xyy_xxz, g_x_0_xyy_xyy, g_x_0_xyy_xyyy, g_x_0_xyy_xyyz, g_x_0_xyy_xyz, g_x_0_xyy_xyzz, g_x_0_xyy_xzz, g_x_0_xyy_yyy, g_x_0_xyy_yyyy, g_x_0_xyy_yyyz, g_x_0_xyy_yyz, g_x_0_xyy_yyzz, g_x_0_xyy_yzz, g_x_0_xyy_yzzz, g_x_0_xyy_zzz, g_x_0_xyyy_xxx, g_x_0_xyyy_xxy, g_x_0_xyyy_xxz, g_x_0_xyyy_xyy, g_x_0_xyyy_xyz, g_x_0_xyyy_xzz, g_x_0_xyyy_yyy, g_x_0_xyyy_yyz, g_x_0_xyyy_yzz, g_x_0_xyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_xxx[k] = -g_x_0_xyy_xxx[k] * cd_y[k] + g_x_0_xyy_xxxy[k];

                g_x_0_xyyy_xxy[k] = -g_x_0_xyy_xxy[k] * cd_y[k] + g_x_0_xyy_xxyy[k];

                g_x_0_xyyy_xxz[k] = -g_x_0_xyy_xxz[k] * cd_y[k] + g_x_0_xyy_xxyz[k];

                g_x_0_xyyy_xyy[k] = -g_x_0_xyy_xyy[k] * cd_y[k] + g_x_0_xyy_xyyy[k];

                g_x_0_xyyy_xyz[k] = -g_x_0_xyy_xyz[k] * cd_y[k] + g_x_0_xyy_xyyz[k];

                g_x_0_xyyy_xzz[k] = -g_x_0_xyy_xzz[k] * cd_y[k] + g_x_0_xyy_xyzz[k];

                g_x_0_xyyy_yyy[k] = -g_x_0_xyy_yyy[k] * cd_y[k] + g_x_0_xyy_yyyy[k];

                g_x_0_xyyy_yyz[k] = -g_x_0_xyy_yyz[k] * cd_y[k] + g_x_0_xyy_yyyz[k];

                g_x_0_xyyy_yzz[k] = -g_x_0_xyy_yzz[k] * cd_y[k] + g_x_0_xyy_yyzz[k];

                g_x_0_xyyy_zzz[k] = -g_x_0_xyy_zzz[k] * cd_y[k] + g_x_0_xyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xyyz_xxx, g_x_0_xyyz_xxy, g_x_0_xyyz_xxz, g_x_0_xyyz_xyy, g_x_0_xyyz_xyz, g_x_0_xyyz_xzz, g_x_0_xyyz_yyy, g_x_0_xyyz_yyz, g_x_0_xyyz_yzz, g_x_0_xyyz_zzz, g_x_0_xyz_xxx, g_x_0_xyz_xxxy, g_x_0_xyz_xxy, g_x_0_xyz_xxyy, g_x_0_xyz_xxyz, g_x_0_xyz_xxz, g_x_0_xyz_xyy, g_x_0_xyz_xyyy, g_x_0_xyz_xyyz, g_x_0_xyz_xyz, g_x_0_xyz_xyzz, g_x_0_xyz_xzz, g_x_0_xyz_yyy, g_x_0_xyz_yyyy, g_x_0_xyz_yyyz, g_x_0_xyz_yyz, g_x_0_xyz_yyzz, g_x_0_xyz_yzz, g_x_0_xyz_yzzz, g_x_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_xxx[k] = -g_x_0_xyz_xxx[k] * cd_y[k] + g_x_0_xyz_xxxy[k];

                g_x_0_xyyz_xxy[k] = -g_x_0_xyz_xxy[k] * cd_y[k] + g_x_0_xyz_xxyy[k];

                g_x_0_xyyz_xxz[k] = -g_x_0_xyz_xxz[k] * cd_y[k] + g_x_0_xyz_xxyz[k];

                g_x_0_xyyz_xyy[k] = -g_x_0_xyz_xyy[k] * cd_y[k] + g_x_0_xyz_xyyy[k];

                g_x_0_xyyz_xyz[k] = -g_x_0_xyz_xyz[k] * cd_y[k] + g_x_0_xyz_xyyz[k];

                g_x_0_xyyz_xzz[k] = -g_x_0_xyz_xzz[k] * cd_y[k] + g_x_0_xyz_xyzz[k];

                g_x_0_xyyz_yyy[k] = -g_x_0_xyz_yyy[k] * cd_y[k] + g_x_0_xyz_yyyy[k];

                g_x_0_xyyz_yyz[k] = -g_x_0_xyz_yyz[k] * cd_y[k] + g_x_0_xyz_yyyz[k];

                g_x_0_xyyz_yzz[k] = -g_x_0_xyz_yzz[k] * cd_y[k] + g_x_0_xyz_yyzz[k];

                g_x_0_xyyz_zzz[k] = -g_x_0_xyz_zzz[k] * cd_y[k] + g_x_0_xyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xyzz_xxx, g_x_0_xyzz_xxy, g_x_0_xyzz_xxz, g_x_0_xyzz_xyy, g_x_0_xyzz_xyz, g_x_0_xyzz_xzz, g_x_0_xyzz_yyy, g_x_0_xyzz_yyz, g_x_0_xyzz_yzz, g_x_0_xyzz_zzz, g_x_0_xzz_xxx, g_x_0_xzz_xxxy, g_x_0_xzz_xxy, g_x_0_xzz_xxyy, g_x_0_xzz_xxyz, g_x_0_xzz_xxz, g_x_0_xzz_xyy, g_x_0_xzz_xyyy, g_x_0_xzz_xyyz, g_x_0_xzz_xyz, g_x_0_xzz_xyzz, g_x_0_xzz_xzz, g_x_0_xzz_yyy, g_x_0_xzz_yyyy, g_x_0_xzz_yyyz, g_x_0_xzz_yyz, g_x_0_xzz_yyzz, g_x_0_xzz_yzz, g_x_0_xzz_yzzz, g_x_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_xxx[k] = -g_x_0_xzz_xxx[k] * cd_y[k] + g_x_0_xzz_xxxy[k];

                g_x_0_xyzz_xxy[k] = -g_x_0_xzz_xxy[k] * cd_y[k] + g_x_0_xzz_xxyy[k];

                g_x_0_xyzz_xxz[k] = -g_x_0_xzz_xxz[k] * cd_y[k] + g_x_0_xzz_xxyz[k];

                g_x_0_xyzz_xyy[k] = -g_x_0_xzz_xyy[k] * cd_y[k] + g_x_0_xzz_xyyy[k];

                g_x_0_xyzz_xyz[k] = -g_x_0_xzz_xyz[k] * cd_y[k] + g_x_0_xzz_xyyz[k];

                g_x_0_xyzz_xzz[k] = -g_x_0_xzz_xzz[k] * cd_y[k] + g_x_0_xzz_xyzz[k];

                g_x_0_xyzz_yyy[k] = -g_x_0_xzz_yyy[k] * cd_y[k] + g_x_0_xzz_yyyy[k];

                g_x_0_xyzz_yyz[k] = -g_x_0_xzz_yyz[k] * cd_y[k] + g_x_0_xzz_yyyz[k];

                g_x_0_xyzz_yzz[k] = -g_x_0_xzz_yzz[k] * cd_y[k] + g_x_0_xzz_yyzz[k];

                g_x_0_xyzz_zzz[k] = -g_x_0_xzz_zzz[k] * cd_y[k] + g_x_0_xzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_xzz_xxx, g_x_0_xzz_xxxz, g_x_0_xzz_xxy, g_x_0_xzz_xxyz, g_x_0_xzz_xxz, g_x_0_xzz_xxzz, g_x_0_xzz_xyy, g_x_0_xzz_xyyz, g_x_0_xzz_xyz, g_x_0_xzz_xyzz, g_x_0_xzz_xzz, g_x_0_xzz_xzzz, g_x_0_xzz_yyy, g_x_0_xzz_yyyz, g_x_0_xzz_yyz, g_x_0_xzz_yyzz, g_x_0_xzz_yzz, g_x_0_xzz_yzzz, g_x_0_xzz_zzz, g_x_0_xzz_zzzz, g_x_0_xzzz_xxx, g_x_0_xzzz_xxy, g_x_0_xzzz_xxz, g_x_0_xzzz_xyy, g_x_0_xzzz_xyz, g_x_0_xzzz_xzz, g_x_0_xzzz_yyy, g_x_0_xzzz_yyz, g_x_0_xzzz_yzz, g_x_0_xzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_xxx[k] = -g_x_0_xzz_xxx[k] * cd_z[k] + g_x_0_xzz_xxxz[k];

                g_x_0_xzzz_xxy[k] = -g_x_0_xzz_xxy[k] * cd_z[k] + g_x_0_xzz_xxyz[k];

                g_x_0_xzzz_xxz[k] = -g_x_0_xzz_xxz[k] * cd_z[k] + g_x_0_xzz_xxzz[k];

                g_x_0_xzzz_xyy[k] = -g_x_0_xzz_xyy[k] * cd_z[k] + g_x_0_xzz_xyyz[k];

                g_x_0_xzzz_xyz[k] = -g_x_0_xzz_xyz[k] * cd_z[k] + g_x_0_xzz_xyzz[k];

                g_x_0_xzzz_xzz[k] = -g_x_0_xzz_xzz[k] * cd_z[k] + g_x_0_xzz_xzzz[k];

                g_x_0_xzzz_yyy[k] = -g_x_0_xzz_yyy[k] * cd_z[k] + g_x_0_xzz_yyyz[k];

                g_x_0_xzzz_yyz[k] = -g_x_0_xzz_yyz[k] * cd_z[k] + g_x_0_xzz_yyzz[k];

                g_x_0_xzzz_yzz[k] = -g_x_0_xzz_yzz[k] * cd_z[k] + g_x_0_xzz_yzzz[k];

                g_x_0_xzzz_zzz[k] = -g_x_0_xzz_zzz[k] * cd_z[k] + g_x_0_xzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yyy_xxx, g_x_0_yyy_xxxy, g_x_0_yyy_xxy, g_x_0_yyy_xxyy, g_x_0_yyy_xxyz, g_x_0_yyy_xxz, g_x_0_yyy_xyy, g_x_0_yyy_xyyy, g_x_0_yyy_xyyz, g_x_0_yyy_xyz, g_x_0_yyy_xyzz, g_x_0_yyy_xzz, g_x_0_yyy_yyy, g_x_0_yyy_yyyy, g_x_0_yyy_yyyz, g_x_0_yyy_yyz, g_x_0_yyy_yyzz, g_x_0_yyy_yzz, g_x_0_yyy_yzzz, g_x_0_yyy_zzz, g_x_0_yyyy_xxx, g_x_0_yyyy_xxy, g_x_0_yyyy_xxz, g_x_0_yyyy_xyy, g_x_0_yyyy_xyz, g_x_0_yyyy_xzz, g_x_0_yyyy_yyy, g_x_0_yyyy_yyz, g_x_0_yyyy_yzz, g_x_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_xxx[k] = -g_x_0_yyy_xxx[k] * cd_y[k] + g_x_0_yyy_xxxy[k];

                g_x_0_yyyy_xxy[k] = -g_x_0_yyy_xxy[k] * cd_y[k] + g_x_0_yyy_xxyy[k];

                g_x_0_yyyy_xxz[k] = -g_x_0_yyy_xxz[k] * cd_y[k] + g_x_0_yyy_xxyz[k];

                g_x_0_yyyy_xyy[k] = -g_x_0_yyy_xyy[k] * cd_y[k] + g_x_0_yyy_xyyy[k];

                g_x_0_yyyy_xyz[k] = -g_x_0_yyy_xyz[k] * cd_y[k] + g_x_0_yyy_xyyz[k];

                g_x_0_yyyy_xzz[k] = -g_x_0_yyy_xzz[k] * cd_y[k] + g_x_0_yyy_xyzz[k];

                g_x_0_yyyy_yyy[k] = -g_x_0_yyy_yyy[k] * cd_y[k] + g_x_0_yyy_yyyy[k];

                g_x_0_yyyy_yyz[k] = -g_x_0_yyy_yyz[k] * cd_y[k] + g_x_0_yyy_yyyz[k];

                g_x_0_yyyy_yzz[k] = -g_x_0_yyy_yzz[k] * cd_y[k] + g_x_0_yyy_yyzz[k];

                g_x_0_yyyy_zzz[k] = -g_x_0_yyy_zzz[k] * cd_y[k] + g_x_0_yyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yyyz_xxx, g_x_0_yyyz_xxy, g_x_0_yyyz_xxz, g_x_0_yyyz_xyy, g_x_0_yyyz_xyz, g_x_0_yyyz_xzz, g_x_0_yyyz_yyy, g_x_0_yyyz_yyz, g_x_0_yyyz_yzz, g_x_0_yyyz_zzz, g_x_0_yyz_xxx, g_x_0_yyz_xxxy, g_x_0_yyz_xxy, g_x_0_yyz_xxyy, g_x_0_yyz_xxyz, g_x_0_yyz_xxz, g_x_0_yyz_xyy, g_x_0_yyz_xyyy, g_x_0_yyz_xyyz, g_x_0_yyz_xyz, g_x_0_yyz_xyzz, g_x_0_yyz_xzz, g_x_0_yyz_yyy, g_x_0_yyz_yyyy, g_x_0_yyz_yyyz, g_x_0_yyz_yyz, g_x_0_yyz_yyzz, g_x_0_yyz_yzz, g_x_0_yyz_yzzz, g_x_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_xxx[k] = -g_x_0_yyz_xxx[k] * cd_y[k] + g_x_0_yyz_xxxy[k];

                g_x_0_yyyz_xxy[k] = -g_x_0_yyz_xxy[k] * cd_y[k] + g_x_0_yyz_xxyy[k];

                g_x_0_yyyz_xxz[k] = -g_x_0_yyz_xxz[k] * cd_y[k] + g_x_0_yyz_xxyz[k];

                g_x_0_yyyz_xyy[k] = -g_x_0_yyz_xyy[k] * cd_y[k] + g_x_0_yyz_xyyy[k];

                g_x_0_yyyz_xyz[k] = -g_x_0_yyz_xyz[k] * cd_y[k] + g_x_0_yyz_xyyz[k];

                g_x_0_yyyz_xzz[k] = -g_x_0_yyz_xzz[k] * cd_y[k] + g_x_0_yyz_xyzz[k];

                g_x_0_yyyz_yyy[k] = -g_x_0_yyz_yyy[k] * cd_y[k] + g_x_0_yyz_yyyy[k];

                g_x_0_yyyz_yyz[k] = -g_x_0_yyz_yyz[k] * cd_y[k] + g_x_0_yyz_yyyz[k];

                g_x_0_yyyz_yzz[k] = -g_x_0_yyz_yzz[k] * cd_y[k] + g_x_0_yyz_yyzz[k];

                g_x_0_yyyz_zzz[k] = -g_x_0_yyz_zzz[k] * cd_y[k] + g_x_0_yyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yyzz_xxx, g_x_0_yyzz_xxy, g_x_0_yyzz_xxz, g_x_0_yyzz_xyy, g_x_0_yyzz_xyz, g_x_0_yyzz_xzz, g_x_0_yyzz_yyy, g_x_0_yyzz_yyz, g_x_0_yyzz_yzz, g_x_0_yyzz_zzz, g_x_0_yzz_xxx, g_x_0_yzz_xxxy, g_x_0_yzz_xxy, g_x_0_yzz_xxyy, g_x_0_yzz_xxyz, g_x_0_yzz_xxz, g_x_0_yzz_xyy, g_x_0_yzz_xyyy, g_x_0_yzz_xyyz, g_x_0_yzz_xyz, g_x_0_yzz_xyzz, g_x_0_yzz_xzz, g_x_0_yzz_yyy, g_x_0_yzz_yyyy, g_x_0_yzz_yyyz, g_x_0_yzz_yyz, g_x_0_yzz_yyzz, g_x_0_yzz_yzz, g_x_0_yzz_yzzz, g_x_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_xxx[k] = -g_x_0_yzz_xxx[k] * cd_y[k] + g_x_0_yzz_xxxy[k];

                g_x_0_yyzz_xxy[k] = -g_x_0_yzz_xxy[k] * cd_y[k] + g_x_0_yzz_xxyy[k];

                g_x_0_yyzz_xxz[k] = -g_x_0_yzz_xxz[k] * cd_y[k] + g_x_0_yzz_xxyz[k];

                g_x_0_yyzz_xyy[k] = -g_x_0_yzz_xyy[k] * cd_y[k] + g_x_0_yzz_xyyy[k];

                g_x_0_yyzz_xyz[k] = -g_x_0_yzz_xyz[k] * cd_y[k] + g_x_0_yzz_xyyz[k];

                g_x_0_yyzz_xzz[k] = -g_x_0_yzz_xzz[k] * cd_y[k] + g_x_0_yzz_xyzz[k];

                g_x_0_yyzz_yyy[k] = -g_x_0_yzz_yyy[k] * cd_y[k] + g_x_0_yzz_yyyy[k];

                g_x_0_yyzz_yyz[k] = -g_x_0_yzz_yyz[k] * cd_y[k] + g_x_0_yzz_yyyz[k];

                g_x_0_yyzz_yzz[k] = -g_x_0_yzz_yzz[k] * cd_y[k] + g_x_0_yzz_yyzz[k];

                g_x_0_yyzz_zzz[k] = -g_x_0_yzz_zzz[k] * cd_y[k] + g_x_0_yzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yzzz_xxx, g_x_0_yzzz_xxy, g_x_0_yzzz_xxz, g_x_0_yzzz_xyy, g_x_0_yzzz_xyz, g_x_0_yzzz_xzz, g_x_0_yzzz_yyy, g_x_0_yzzz_yyz, g_x_0_yzzz_yzz, g_x_0_yzzz_zzz, g_x_0_zzz_xxx, g_x_0_zzz_xxxy, g_x_0_zzz_xxy, g_x_0_zzz_xxyy, g_x_0_zzz_xxyz, g_x_0_zzz_xxz, g_x_0_zzz_xyy, g_x_0_zzz_xyyy, g_x_0_zzz_xyyz, g_x_0_zzz_xyz, g_x_0_zzz_xyzz, g_x_0_zzz_xzz, g_x_0_zzz_yyy, g_x_0_zzz_yyyy, g_x_0_zzz_yyyz, g_x_0_zzz_yyz, g_x_0_zzz_yyzz, g_x_0_zzz_yzz, g_x_0_zzz_yzzz, g_x_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_xxx[k] = -g_x_0_zzz_xxx[k] * cd_y[k] + g_x_0_zzz_xxxy[k];

                g_x_0_yzzz_xxy[k] = -g_x_0_zzz_xxy[k] * cd_y[k] + g_x_0_zzz_xxyy[k];

                g_x_0_yzzz_xxz[k] = -g_x_0_zzz_xxz[k] * cd_y[k] + g_x_0_zzz_xxyz[k];

                g_x_0_yzzz_xyy[k] = -g_x_0_zzz_xyy[k] * cd_y[k] + g_x_0_zzz_xyyy[k];

                g_x_0_yzzz_xyz[k] = -g_x_0_zzz_xyz[k] * cd_y[k] + g_x_0_zzz_xyyz[k];

                g_x_0_yzzz_xzz[k] = -g_x_0_zzz_xzz[k] * cd_y[k] + g_x_0_zzz_xyzz[k];

                g_x_0_yzzz_yyy[k] = -g_x_0_zzz_yyy[k] * cd_y[k] + g_x_0_zzz_yyyy[k];

                g_x_0_yzzz_yyz[k] = -g_x_0_zzz_yyz[k] * cd_y[k] + g_x_0_zzz_yyyz[k];

                g_x_0_yzzz_yzz[k] = -g_x_0_zzz_yzz[k] * cd_y[k] + g_x_0_zzz_yyzz[k];

                g_x_0_yzzz_zzz[k] = -g_x_0_zzz_zzz[k] * cd_y[k] + g_x_0_zzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_zzz_xxx, g_x_0_zzz_xxxz, g_x_0_zzz_xxy, g_x_0_zzz_xxyz, g_x_0_zzz_xxz, g_x_0_zzz_xxzz, g_x_0_zzz_xyy, g_x_0_zzz_xyyz, g_x_0_zzz_xyz, g_x_0_zzz_xyzz, g_x_0_zzz_xzz, g_x_0_zzz_xzzz, g_x_0_zzz_yyy, g_x_0_zzz_yyyz, g_x_0_zzz_yyz, g_x_0_zzz_yyzz, g_x_0_zzz_yzz, g_x_0_zzz_yzzz, g_x_0_zzz_zzz, g_x_0_zzz_zzzz, g_x_0_zzzz_xxx, g_x_0_zzzz_xxy, g_x_0_zzzz_xxz, g_x_0_zzzz_xyy, g_x_0_zzzz_xyz, g_x_0_zzzz_xzz, g_x_0_zzzz_yyy, g_x_0_zzzz_yyz, g_x_0_zzzz_yzz, g_x_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_xxx[k] = -g_x_0_zzz_xxx[k] * cd_z[k] + g_x_0_zzz_xxxz[k];

                g_x_0_zzzz_xxy[k] = -g_x_0_zzz_xxy[k] * cd_z[k] + g_x_0_zzz_xxyz[k];

                g_x_0_zzzz_xxz[k] = -g_x_0_zzz_xxz[k] * cd_z[k] + g_x_0_zzz_xxzz[k];

                g_x_0_zzzz_xyy[k] = -g_x_0_zzz_xyy[k] * cd_z[k] + g_x_0_zzz_xyyz[k];

                g_x_0_zzzz_xyz[k] = -g_x_0_zzz_xyz[k] * cd_z[k] + g_x_0_zzz_xyzz[k];

                g_x_0_zzzz_xzz[k] = -g_x_0_zzz_xzz[k] * cd_z[k] + g_x_0_zzz_xzzz[k];

                g_x_0_zzzz_yyy[k] = -g_x_0_zzz_yyy[k] * cd_z[k] + g_x_0_zzz_yyyz[k];

                g_x_0_zzzz_yyz[k] = -g_x_0_zzz_yyz[k] * cd_z[k] + g_x_0_zzz_yyzz[k];

                g_x_0_zzzz_yzz[k] = -g_x_0_zzz_yzz[k] * cd_z[k] + g_x_0_zzz_yzzz[k];

                g_x_0_zzzz_zzz[k] = -g_x_0_zzz_zzz[k] * cd_z[k] + g_x_0_zzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxx_xxx, g_y_0_xxx_xxxx, g_y_0_xxx_xxxy, g_y_0_xxx_xxxz, g_y_0_xxx_xxy, g_y_0_xxx_xxyy, g_y_0_xxx_xxyz, g_y_0_xxx_xxz, g_y_0_xxx_xxzz, g_y_0_xxx_xyy, g_y_0_xxx_xyyy, g_y_0_xxx_xyyz, g_y_0_xxx_xyz, g_y_0_xxx_xyzz, g_y_0_xxx_xzz, g_y_0_xxx_xzzz, g_y_0_xxx_yyy, g_y_0_xxx_yyz, g_y_0_xxx_yzz, g_y_0_xxx_zzz, g_y_0_xxxx_xxx, g_y_0_xxxx_xxy, g_y_0_xxxx_xxz, g_y_0_xxxx_xyy, g_y_0_xxxx_xyz, g_y_0_xxxx_xzz, g_y_0_xxxx_yyy, g_y_0_xxxx_yyz, g_y_0_xxxx_yzz, g_y_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_xxx[k] = -g_y_0_xxx_xxx[k] * cd_x[k] + g_y_0_xxx_xxxx[k];

                g_y_0_xxxx_xxy[k] = -g_y_0_xxx_xxy[k] * cd_x[k] + g_y_0_xxx_xxxy[k];

                g_y_0_xxxx_xxz[k] = -g_y_0_xxx_xxz[k] * cd_x[k] + g_y_0_xxx_xxxz[k];

                g_y_0_xxxx_xyy[k] = -g_y_0_xxx_xyy[k] * cd_x[k] + g_y_0_xxx_xxyy[k];

                g_y_0_xxxx_xyz[k] = -g_y_0_xxx_xyz[k] * cd_x[k] + g_y_0_xxx_xxyz[k];

                g_y_0_xxxx_xzz[k] = -g_y_0_xxx_xzz[k] * cd_x[k] + g_y_0_xxx_xxzz[k];

                g_y_0_xxxx_yyy[k] = -g_y_0_xxx_yyy[k] * cd_x[k] + g_y_0_xxx_xyyy[k];

                g_y_0_xxxx_yyz[k] = -g_y_0_xxx_yyz[k] * cd_x[k] + g_y_0_xxx_xyyz[k];

                g_y_0_xxxx_yzz[k] = -g_y_0_xxx_yzz[k] * cd_x[k] + g_y_0_xxx_xyzz[k];

                g_y_0_xxxx_zzz[k] = -g_y_0_xxx_zzz[k] * cd_x[k] + g_y_0_xxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxxy_xxx, g_y_0_xxxy_xxy, g_y_0_xxxy_xxz, g_y_0_xxxy_xyy, g_y_0_xxxy_xyz, g_y_0_xxxy_xzz, g_y_0_xxxy_yyy, g_y_0_xxxy_yyz, g_y_0_xxxy_yzz, g_y_0_xxxy_zzz, g_y_0_xxy_xxx, g_y_0_xxy_xxxx, g_y_0_xxy_xxxy, g_y_0_xxy_xxxz, g_y_0_xxy_xxy, g_y_0_xxy_xxyy, g_y_0_xxy_xxyz, g_y_0_xxy_xxz, g_y_0_xxy_xxzz, g_y_0_xxy_xyy, g_y_0_xxy_xyyy, g_y_0_xxy_xyyz, g_y_0_xxy_xyz, g_y_0_xxy_xyzz, g_y_0_xxy_xzz, g_y_0_xxy_xzzz, g_y_0_xxy_yyy, g_y_0_xxy_yyz, g_y_0_xxy_yzz, g_y_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_xxx[k] = -g_y_0_xxy_xxx[k] * cd_x[k] + g_y_0_xxy_xxxx[k];

                g_y_0_xxxy_xxy[k] = -g_y_0_xxy_xxy[k] * cd_x[k] + g_y_0_xxy_xxxy[k];

                g_y_0_xxxy_xxz[k] = -g_y_0_xxy_xxz[k] * cd_x[k] + g_y_0_xxy_xxxz[k];

                g_y_0_xxxy_xyy[k] = -g_y_0_xxy_xyy[k] * cd_x[k] + g_y_0_xxy_xxyy[k];

                g_y_0_xxxy_xyz[k] = -g_y_0_xxy_xyz[k] * cd_x[k] + g_y_0_xxy_xxyz[k];

                g_y_0_xxxy_xzz[k] = -g_y_0_xxy_xzz[k] * cd_x[k] + g_y_0_xxy_xxzz[k];

                g_y_0_xxxy_yyy[k] = -g_y_0_xxy_yyy[k] * cd_x[k] + g_y_0_xxy_xyyy[k];

                g_y_0_xxxy_yyz[k] = -g_y_0_xxy_yyz[k] * cd_x[k] + g_y_0_xxy_xyyz[k];

                g_y_0_xxxy_yzz[k] = -g_y_0_xxy_yzz[k] * cd_x[k] + g_y_0_xxy_xyzz[k];

                g_y_0_xxxy_zzz[k] = -g_y_0_xxy_zzz[k] * cd_x[k] + g_y_0_xxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxxz_xxx, g_y_0_xxxz_xxy, g_y_0_xxxz_xxz, g_y_0_xxxz_xyy, g_y_0_xxxz_xyz, g_y_0_xxxz_xzz, g_y_0_xxxz_yyy, g_y_0_xxxz_yyz, g_y_0_xxxz_yzz, g_y_0_xxxz_zzz, g_y_0_xxz_xxx, g_y_0_xxz_xxxx, g_y_0_xxz_xxxy, g_y_0_xxz_xxxz, g_y_0_xxz_xxy, g_y_0_xxz_xxyy, g_y_0_xxz_xxyz, g_y_0_xxz_xxz, g_y_0_xxz_xxzz, g_y_0_xxz_xyy, g_y_0_xxz_xyyy, g_y_0_xxz_xyyz, g_y_0_xxz_xyz, g_y_0_xxz_xyzz, g_y_0_xxz_xzz, g_y_0_xxz_xzzz, g_y_0_xxz_yyy, g_y_0_xxz_yyz, g_y_0_xxz_yzz, g_y_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_xxx[k] = -g_y_0_xxz_xxx[k] * cd_x[k] + g_y_0_xxz_xxxx[k];

                g_y_0_xxxz_xxy[k] = -g_y_0_xxz_xxy[k] * cd_x[k] + g_y_0_xxz_xxxy[k];

                g_y_0_xxxz_xxz[k] = -g_y_0_xxz_xxz[k] * cd_x[k] + g_y_0_xxz_xxxz[k];

                g_y_0_xxxz_xyy[k] = -g_y_0_xxz_xyy[k] * cd_x[k] + g_y_0_xxz_xxyy[k];

                g_y_0_xxxz_xyz[k] = -g_y_0_xxz_xyz[k] * cd_x[k] + g_y_0_xxz_xxyz[k];

                g_y_0_xxxz_xzz[k] = -g_y_0_xxz_xzz[k] * cd_x[k] + g_y_0_xxz_xxzz[k];

                g_y_0_xxxz_yyy[k] = -g_y_0_xxz_yyy[k] * cd_x[k] + g_y_0_xxz_xyyy[k];

                g_y_0_xxxz_yyz[k] = -g_y_0_xxz_yyz[k] * cd_x[k] + g_y_0_xxz_xyyz[k];

                g_y_0_xxxz_yzz[k] = -g_y_0_xxz_yzz[k] * cd_x[k] + g_y_0_xxz_xyzz[k];

                g_y_0_xxxz_zzz[k] = -g_y_0_xxz_zzz[k] * cd_x[k] + g_y_0_xxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxyy_xxx, g_y_0_xxyy_xxy, g_y_0_xxyy_xxz, g_y_0_xxyy_xyy, g_y_0_xxyy_xyz, g_y_0_xxyy_xzz, g_y_0_xxyy_yyy, g_y_0_xxyy_yyz, g_y_0_xxyy_yzz, g_y_0_xxyy_zzz, g_y_0_xyy_xxx, g_y_0_xyy_xxxx, g_y_0_xyy_xxxy, g_y_0_xyy_xxxz, g_y_0_xyy_xxy, g_y_0_xyy_xxyy, g_y_0_xyy_xxyz, g_y_0_xyy_xxz, g_y_0_xyy_xxzz, g_y_0_xyy_xyy, g_y_0_xyy_xyyy, g_y_0_xyy_xyyz, g_y_0_xyy_xyz, g_y_0_xyy_xyzz, g_y_0_xyy_xzz, g_y_0_xyy_xzzz, g_y_0_xyy_yyy, g_y_0_xyy_yyz, g_y_0_xyy_yzz, g_y_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_xxx[k] = -g_y_0_xyy_xxx[k] * cd_x[k] + g_y_0_xyy_xxxx[k];

                g_y_0_xxyy_xxy[k] = -g_y_0_xyy_xxy[k] * cd_x[k] + g_y_0_xyy_xxxy[k];

                g_y_0_xxyy_xxz[k] = -g_y_0_xyy_xxz[k] * cd_x[k] + g_y_0_xyy_xxxz[k];

                g_y_0_xxyy_xyy[k] = -g_y_0_xyy_xyy[k] * cd_x[k] + g_y_0_xyy_xxyy[k];

                g_y_0_xxyy_xyz[k] = -g_y_0_xyy_xyz[k] * cd_x[k] + g_y_0_xyy_xxyz[k];

                g_y_0_xxyy_xzz[k] = -g_y_0_xyy_xzz[k] * cd_x[k] + g_y_0_xyy_xxzz[k];

                g_y_0_xxyy_yyy[k] = -g_y_0_xyy_yyy[k] * cd_x[k] + g_y_0_xyy_xyyy[k];

                g_y_0_xxyy_yyz[k] = -g_y_0_xyy_yyz[k] * cd_x[k] + g_y_0_xyy_xyyz[k];

                g_y_0_xxyy_yzz[k] = -g_y_0_xyy_yzz[k] * cd_x[k] + g_y_0_xyy_xyzz[k];

                g_y_0_xxyy_zzz[k] = -g_y_0_xyy_zzz[k] * cd_x[k] + g_y_0_xyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxyz_xxx, g_y_0_xxyz_xxy, g_y_0_xxyz_xxz, g_y_0_xxyz_xyy, g_y_0_xxyz_xyz, g_y_0_xxyz_xzz, g_y_0_xxyz_yyy, g_y_0_xxyz_yyz, g_y_0_xxyz_yzz, g_y_0_xxyz_zzz, g_y_0_xyz_xxx, g_y_0_xyz_xxxx, g_y_0_xyz_xxxy, g_y_0_xyz_xxxz, g_y_0_xyz_xxy, g_y_0_xyz_xxyy, g_y_0_xyz_xxyz, g_y_0_xyz_xxz, g_y_0_xyz_xxzz, g_y_0_xyz_xyy, g_y_0_xyz_xyyy, g_y_0_xyz_xyyz, g_y_0_xyz_xyz, g_y_0_xyz_xyzz, g_y_0_xyz_xzz, g_y_0_xyz_xzzz, g_y_0_xyz_yyy, g_y_0_xyz_yyz, g_y_0_xyz_yzz, g_y_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_xxx[k] = -g_y_0_xyz_xxx[k] * cd_x[k] + g_y_0_xyz_xxxx[k];

                g_y_0_xxyz_xxy[k] = -g_y_0_xyz_xxy[k] * cd_x[k] + g_y_0_xyz_xxxy[k];

                g_y_0_xxyz_xxz[k] = -g_y_0_xyz_xxz[k] * cd_x[k] + g_y_0_xyz_xxxz[k];

                g_y_0_xxyz_xyy[k] = -g_y_0_xyz_xyy[k] * cd_x[k] + g_y_0_xyz_xxyy[k];

                g_y_0_xxyz_xyz[k] = -g_y_0_xyz_xyz[k] * cd_x[k] + g_y_0_xyz_xxyz[k];

                g_y_0_xxyz_xzz[k] = -g_y_0_xyz_xzz[k] * cd_x[k] + g_y_0_xyz_xxzz[k];

                g_y_0_xxyz_yyy[k] = -g_y_0_xyz_yyy[k] * cd_x[k] + g_y_0_xyz_xyyy[k];

                g_y_0_xxyz_yyz[k] = -g_y_0_xyz_yyz[k] * cd_x[k] + g_y_0_xyz_xyyz[k];

                g_y_0_xxyz_yzz[k] = -g_y_0_xyz_yzz[k] * cd_x[k] + g_y_0_xyz_xyzz[k];

                g_y_0_xxyz_zzz[k] = -g_y_0_xyz_zzz[k] * cd_x[k] + g_y_0_xyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxzz_xxx, g_y_0_xxzz_xxy, g_y_0_xxzz_xxz, g_y_0_xxzz_xyy, g_y_0_xxzz_xyz, g_y_0_xxzz_xzz, g_y_0_xxzz_yyy, g_y_0_xxzz_yyz, g_y_0_xxzz_yzz, g_y_0_xxzz_zzz, g_y_0_xzz_xxx, g_y_0_xzz_xxxx, g_y_0_xzz_xxxy, g_y_0_xzz_xxxz, g_y_0_xzz_xxy, g_y_0_xzz_xxyy, g_y_0_xzz_xxyz, g_y_0_xzz_xxz, g_y_0_xzz_xxzz, g_y_0_xzz_xyy, g_y_0_xzz_xyyy, g_y_0_xzz_xyyz, g_y_0_xzz_xyz, g_y_0_xzz_xyzz, g_y_0_xzz_xzz, g_y_0_xzz_xzzz, g_y_0_xzz_yyy, g_y_0_xzz_yyz, g_y_0_xzz_yzz, g_y_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_xxx[k] = -g_y_0_xzz_xxx[k] * cd_x[k] + g_y_0_xzz_xxxx[k];

                g_y_0_xxzz_xxy[k] = -g_y_0_xzz_xxy[k] * cd_x[k] + g_y_0_xzz_xxxy[k];

                g_y_0_xxzz_xxz[k] = -g_y_0_xzz_xxz[k] * cd_x[k] + g_y_0_xzz_xxxz[k];

                g_y_0_xxzz_xyy[k] = -g_y_0_xzz_xyy[k] * cd_x[k] + g_y_0_xzz_xxyy[k];

                g_y_0_xxzz_xyz[k] = -g_y_0_xzz_xyz[k] * cd_x[k] + g_y_0_xzz_xxyz[k];

                g_y_0_xxzz_xzz[k] = -g_y_0_xzz_xzz[k] * cd_x[k] + g_y_0_xzz_xxzz[k];

                g_y_0_xxzz_yyy[k] = -g_y_0_xzz_yyy[k] * cd_x[k] + g_y_0_xzz_xyyy[k];

                g_y_0_xxzz_yyz[k] = -g_y_0_xzz_yyz[k] * cd_x[k] + g_y_0_xzz_xyyz[k];

                g_y_0_xxzz_yzz[k] = -g_y_0_xzz_yzz[k] * cd_x[k] + g_y_0_xzz_xyzz[k];

                g_y_0_xxzz_zzz[k] = -g_y_0_xzz_zzz[k] * cd_x[k] + g_y_0_xzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xyyy_xxx, g_y_0_xyyy_xxy, g_y_0_xyyy_xxz, g_y_0_xyyy_xyy, g_y_0_xyyy_xyz, g_y_0_xyyy_xzz, g_y_0_xyyy_yyy, g_y_0_xyyy_yyz, g_y_0_xyyy_yzz, g_y_0_xyyy_zzz, g_y_0_yyy_xxx, g_y_0_yyy_xxxx, g_y_0_yyy_xxxy, g_y_0_yyy_xxxz, g_y_0_yyy_xxy, g_y_0_yyy_xxyy, g_y_0_yyy_xxyz, g_y_0_yyy_xxz, g_y_0_yyy_xxzz, g_y_0_yyy_xyy, g_y_0_yyy_xyyy, g_y_0_yyy_xyyz, g_y_0_yyy_xyz, g_y_0_yyy_xyzz, g_y_0_yyy_xzz, g_y_0_yyy_xzzz, g_y_0_yyy_yyy, g_y_0_yyy_yyz, g_y_0_yyy_yzz, g_y_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_xxx[k] = -g_y_0_yyy_xxx[k] * cd_x[k] + g_y_0_yyy_xxxx[k];

                g_y_0_xyyy_xxy[k] = -g_y_0_yyy_xxy[k] * cd_x[k] + g_y_0_yyy_xxxy[k];

                g_y_0_xyyy_xxz[k] = -g_y_0_yyy_xxz[k] * cd_x[k] + g_y_0_yyy_xxxz[k];

                g_y_0_xyyy_xyy[k] = -g_y_0_yyy_xyy[k] * cd_x[k] + g_y_0_yyy_xxyy[k];

                g_y_0_xyyy_xyz[k] = -g_y_0_yyy_xyz[k] * cd_x[k] + g_y_0_yyy_xxyz[k];

                g_y_0_xyyy_xzz[k] = -g_y_0_yyy_xzz[k] * cd_x[k] + g_y_0_yyy_xxzz[k];

                g_y_0_xyyy_yyy[k] = -g_y_0_yyy_yyy[k] * cd_x[k] + g_y_0_yyy_xyyy[k];

                g_y_0_xyyy_yyz[k] = -g_y_0_yyy_yyz[k] * cd_x[k] + g_y_0_yyy_xyyz[k];

                g_y_0_xyyy_yzz[k] = -g_y_0_yyy_yzz[k] * cd_x[k] + g_y_0_yyy_xyzz[k];

                g_y_0_xyyy_zzz[k] = -g_y_0_yyy_zzz[k] * cd_x[k] + g_y_0_yyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xyyz_xxx, g_y_0_xyyz_xxy, g_y_0_xyyz_xxz, g_y_0_xyyz_xyy, g_y_0_xyyz_xyz, g_y_0_xyyz_xzz, g_y_0_xyyz_yyy, g_y_0_xyyz_yyz, g_y_0_xyyz_yzz, g_y_0_xyyz_zzz, g_y_0_yyz_xxx, g_y_0_yyz_xxxx, g_y_0_yyz_xxxy, g_y_0_yyz_xxxz, g_y_0_yyz_xxy, g_y_0_yyz_xxyy, g_y_0_yyz_xxyz, g_y_0_yyz_xxz, g_y_0_yyz_xxzz, g_y_0_yyz_xyy, g_y_0_yyz_xyyy, g_y_0_yyz_xyyz, g_y_0_yyz_xyz, g_y_0_yyz_xyzz, g_y_0_yyz_xzz, g_y_0_yyz_xzzz, g_y_0_yyz_yyy, g_y_0_yyz_yyz, g_y_0_yyz_yzz, g_y_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_xxx[k] = -g_y_0_yyz_xxx[k] * cd_x[k] + g_y_0_yyz_xxxx[k];

                g_y_0_xyyz_xxy[k] = -g_y_0_yyz_xxy[k] * cd_x[k] + g_y_0_yyz_xxxy[k];

                g_y_0_xyyz_xxz[k] = -g_y_0_yyz_xxz[k] * cd_x[k] + g_y_0_yyz_xxxz[k];

                g_y_0_xyyz_xyy[k] = -g_y_0_yyz_xyy[k] * cd_x[k] + g_y_0_yyz_xxyy[k];

                g_y_0_xyyz_xyz[k] = -g_y_0_yyz_xyz[k] * cd_x[k] + g_y_0_yyz_xxyz[k];

                g_y_0_xyyz_xzz[k] = -g_y_0_yyz_xzz[k] * cd_x[k] + g_y_0_yyz_xxzz[k];

                g_y_0_xyyz_yyy[k] = -g_y_0_yyz_yyy[k] * cd_x[k] + g_y_0_yyz_xyyy[k];

                g_y_0_xyyz_yyz[k] = -g_y_0_yyz_yyz[k] * cd_x[k] + g_y_0_yyz_xyyz[k];

                g_y_0_xyyz_yzz[k] = -g_y_0_yyz_yzz[k] * cd_x[k] + g_y_0_yyz_xyzz[k];

                g_y_0_xyyz_zzz[k] = -g_y_0_yyz_zzz[k] * cd_x[k] + g_y_0_yyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xyzz_xxx, g_y_0_xyzz_xxy, g_y_0_xyzz_xxz, g_y_0_xyzz_xyy, g_y_0_xyzz_xyz, g_y_0_xyzz_xzz, g_y_0_xyzz_yyy, g_y_0_xyzz_yyz, g_y_0_xyzz_yzz, g_y_0_xyzz_zzz, g_y_0_yzz_xxx, g_y_0_yzz_xxxx, g_y_0_yzz_xxxy, g_y_0_yzz_xxxz, g_y_0_yzz_xxy, g_y_0_yzz_xxyy, g_y_0_yzz_xxyz, g_y_0_yzz_xxz, g_y_0_yzz_xxzz, g_y_0_yzz_xyy, g_y_0_yzz_xyyy, g_y_0_yzz_xyyz, g_y_0_yzz_xyz, g_y_0_yzz_xyzz, g_y_0_yzz_xzz, g_y_0_yzz_xzzz, g_y_0_yzz_yyy, g_y_0_yzz_yyz, g_y_0_yzz_yzz, g_y_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_xxx[k] = -g_y_0_yzz_xxx[k] * cd_x[k] + g_y_0_yzz_xxxx[k];

                g_y_0_xyzz_xxy[k] = -g_y_0_yzz_xxy[k] * cd_x[k] + g_y_0_yzz_xxxy[k];

                g_y_0_xyzz_xxz[k] = -g_y_0_yzz_xxz[k] * cd_x[k] + g_y_0_yzz_xxxz[k];

                g_y_0_xyzz_xyy[k] = -g_y_0_yzz_xyy[k] * cd_x[k] + g_y_0_yzz_xxyy[k];

                g_y_0_xyzz_xyz[k] = -g_y_0_yzz_xyz[k] * cd_x[k] + g_y_0_yzz_xxyz[k];

                g_y_0_xyzz_xzz[k] = -g_y_0_yzz_xzz[k] * cd_x[k] + g_y_0_yzz_xxzz[k];

                g_y_0_xyzz_yyy[k] = -g_y_0_yzz_yyy[k] * cd_x[k] + g_y_0_yzz_xyyy[k];

                g_y_0_xyzz_yyz[k] = -g_y_0_yzz_yyz[k] * cd_x[k] + g_y_0_yzz_xyyz[k];

                g_y_0_xyzz_yzz[k] = -g_y_0_yzz_yzz[k] * cd_x[k] + g_y_0_yzz_xyzz[k];

                g_y_0_xyzz_zzz[k] = -g_y_0_yzz_zzz[k] * cd_x[k] + g_y_0_yzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xzzz_xxx, g_y_0_xzzz_xxy, g_y_0_xzzz_xxz, g_y_0_xzzz_xyy, g_y_0_xzzz_xyz, g_y_0_xzzz_xzz, g_y_0_xzzz_yyy, g_y_0_xzzz_yyz, g_y_0_xzzz_yzz, g_y_0_xzzz_zzz, g_y_0_zzz_xxx, g_y_0_zzz_xxxx, g_y_0_zzz_xxxy, g_y_0_zzz_xxxz, g_y_0_zzz_xxy, g_y_0_zzz_xxyy, g_y_0_zzz_xxyz, g_y_0_zzz_xxz, g_y_0_zzz_xxzz, g_y_0_zzz_xyy, g_y_0_zzz_xyyy, g_y_0_zzz_xyyz, g_y_0_zzz_xyz, g_y_0_zzz_xyzz, g_y_0_zzz_xzz, g_y_0_zzz_xzzz, g_y_0_zzz_yyy, g_y_0_zzz_yyz, g_y_0_zzz_yzz, g_y_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_xxx[k] = -g_y_0_zzz_xxx[k] * cd_x[k] + g_y_0_zzz_xxxx[k];

                g_y_0_xzzz_xxy[k] = -g_y_0_zzz_xxy[k] * cd_x[k] + g_y_0_zzz_xxxy[k];

                g_y_0_xzzz_xxz[k] = -g_y_0_zzz_xxz[k] * cd_x[k] + g_y_0_zzz_xxxz[k];

                g_y_0_xzzz_xyy[k] = -g_y_0_zzz_xyy[k] * cd_x[k] + g_y_0_zzz_xxyy[k];

                g_y_0_xzzz_xyz[k] = -g_y_0_zzz_xyz[k] * cd_x[k] + g_y_0_zzz_xxyz[k];

                g_y_0_xzzz_xzz[k] = -g_y_0_zzz_xzz[k] * cd_x[k] + g_y_0_zzz_xxzz[k];

                g_y_0_xzzz_yyy[k] = -g_y_0_zzz_yyy[k] * cd_x[k] + g_y_0_zzz_xyyy[k];

                g_y_0_xzzz_yyz[k] = -g_y_0_zzz_yyz[k] * cd_x[k] + g_y_0_zzz_xyyz[k];

                g_y_0_xzzz_yzz[k] = -g_y_0_zzz_yzz[k] * cd_x[k] + g_y_0_zzz_xyzz[k];

                g_y_0_xzzz_zzz[k] = -g_y_0_zzz_zzz[k] * cd_x[k] + g_y_0_zzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_y_0_yyy_xxx, g_y_0_yyy_xxxy, g_y_0_yyy_xxy, g_y_0_yyy_xxyy, g_y_0_yyy_xxyz, g_y_0_yyy_xxz, g_y_0_yyy_xyy, g_y_0_yyy_xyyy, g_y_0_yyy_xyyz, g_y_0_yyy_xyz, g_y_0_yyy_xyzz, g_y_0_yyy_xzz, g_y_0_yyy_yyy, g_y_0_yyy_yyyy, g_y_0_yyy_yyyz, g_y_0_yyy_yyz, g_y_0_yyy_yyzz, g_y_0_yyy_yzz, g_y_0_yyy_yzzz, g_y_0_yyy_zzz, g_y_0_yyyy_xxx, g_y_0_yyyy_xxy, g_y_0_yyyy_xxz, g_y_0_yyyy_xyy, g_y_0_yyyy_xyz, g_y_0_yyyy_xzz, g_y_0_yyyy_yyy, g_y_0_yyyy_yyz, g_y_0_yyyy_yzz, g_y_0_yyyy_zzz, g_yyy_xxx, g_yyy_xxy, g_yyy_xxz, g_yyy_xyy, g_yyy_xyz, g_yyy_xzz, g_yyy_yyy, g_yyy_yyz, g_yyy_yzz, g_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_xxx[k] = -g_yyy_xxx[k] - g_y_0_yyy_xxx[k] * cd_y[k] + g_y_0_yyy_xxxy[k];

                g_y_0_yyyy_xxy[k] = -g_yyy_xxy[k] - g_y_0_yyy_xxy[k] * cd_y[k] + g_y_0_yyy_xxyy[k];

                g_y_0_yyyy_xxz[k] = -g_yyy_xxz[k] - g_y_0_yyy_xxz[k] * cd_y[k] + g_y_0_yyy_xxyz[k];

                g_y_0_yyyy_xyy[k] = -g_yyy_xyy[k] - g_y_0_yyy_xyy[k] * cd_y[k] + g_y_0_yyy_xyyy[k];

                g_y_0_yyyy_xyz[k] = -g_yyy_xyz[k] - g_y_0_yyy_xyz[k] * cd_y[k] + g_y_0_yyy_xyyz[k];

                g_y_0_yyyy_xzz[k] = -g_yyy_xzz[k] - g_y_0_yyy_xzz[k] * cd_y[k] + g_y_0_yyy_xyzz[k];

                g_y_0_yyyy_yyy[k] = -g_yyy_yyy[k] - g_y_0_yyy_yyy[k] * cd_y[k] + g_y_0_yyy_yyyy[k];

                g_y_0_yyyy_yyz[k] = -g_yyy_yyz[k] - g_y_0_yyy_yyz[k] * cd_y[k] + g_y_0_yyy_yyyz[k];

                g_y_0_yyyy_yzz[k] = -g_yyy_yzz[k] - g_y_0_yyy_yzz[k] * cd_y[k] + g_y_0_yyy_yyzz[k];

                g_y_0_yyyy_zzz[k] = -g_yyy_zzz[k] - g_y_0_yyy_zzz[k] * cd_y[k] + g_y_0_yyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_yyy_xxx, g_y_0_yyy_xxxz, g_y_0_yyy_xxy, g_y_0_yyy_xxyz, g_y_0_yyy_xxz, g_y_0_yyy_xxzz, g_y_0_yyy_xyy, g_y_0_yyy_xyyz, g_y_0_yyy_xyz, g_y_0_yyy_xyzz, g_y_0_yyy_xzz, g_y_0_yyy_xzzz, g_y_0_yyy_yyy, g_y_0_yyy_yyyz, g_y_0_yyy_yyz, g_y_0_yyy_yyzz, g_y_0_yyy_yzz, g_y_0_yyy_yzzz, g_y_0_yyy_zzz, g_y_0_yyy_zzzz, g_y_0_yyyz_xxx, g_y_0_yyyz_xxy, g_y_0_yyyz_xxz, g_y_0_yyyz_xyy, g_y_0_yyyz_xyz, g_y_0_yyyz_xzz, g_y_0_yyyz_yyy, g_y_0_yyyz_yyz, g_y_0_yyyz_yzz, g_y_0_yyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_xxx[k] = -g_y_0_yyy_xxx[k] * cd_z[k] + g_y_0_yyy_xxxz[k];

                g_y_0_yyyz_xxy[k] = -g_y_0_yyy_xxy[k] * cd_z[k] + g_y_0_yyy_xxyz[k];

                g_y_0_yyyz_xxz[k] = -g_y_0_yyy_xxz[k] * cd_z[k] + g_y_0_yyy_xxzz[k];

                g_y_0_yyyz_xyy[k] = -g_y_0_yyy_xyy[k] * cd_z[k] + g_y_0_yyy_xyyz[k];

                g_y_0_yyyz_xyz[k] = -g_y_0_yyy_xyz[k] * cd_z[k] + g_y_0_yyy_xyzz[k];

                g_y_0_yyyz_xzz[k] = -g_y_0_yyy_xzz[k] * cd_z[k] + g_y_0_yyy_xzzz[k];

                g_y_0_yyyz_yyy[k] = -g_y_0_yyy_yyy[k] * cd_z[k] + g_y_0_yyy_yyyz[k];

                g_y_0_yyyz_yyz[k] = -g_y_0_yyy_yyz[k] * cd_z[k] + g_y_0_yyy_yyzz[k];

                g_y_0_yyyz_yzz[k] = -g_y_0_yyy_yzz[k] * cd_z[k] + g_y_0_yyy_yzzz[k];

                g_y_0_yyyz_zzz[k] = -g_y_0_yyy_zzz[k] * cd_z[k] + g_y_0_yyy_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_yyz_xxx, g_y_0_yyz_xxxz, g_y_0_yyz_xxy, g_y_0_yyz_xxyz, g_y_0_yyz_xxz, g_y_0_yyz_xxzz, g_y_0_yyz_xyy, g_y_0_yyz_xyyz, g_y_0_yyz_xyz, g_y_0_yyz_xyzz, g_y_0_yyz_xzz, g_y_0_yyz_xzzz, g_y_0_yyz_yyy, g_y_0_yyz_yyyz, g_y_0_yyz_yyz, g_y_0_yyz_yyzz, g_y_0_yyz_yzz, g_y_0_yyz_yzzz, g_y_0_yyz_zzz, g_y_0_yyz_zzzz, g_y_0_yyzz_xxx, g_y_0_yyzz_xxy, g_y_0_yyzz_xxz, g_y_0_yyzz_xyy, g_y_0_yyzz_xyz, g_y_0_yyzz_xzz, g_y_0_yyzz_yyy, g_y_0_yyzz_yyz, g_y_0_yyzz_yzz, g_y_0_yyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_xxx[k] = -g_y_0_yyz_xxx[k] * cd_z[k] + g_y_0_yyz_xxxz[k];

                g_y_0_yyzz_xxy[k] = -g_y_0_yyz_xxy[k] * cd_z[k] + g_y_0_yyz_xxyz[k];

                g_y_0_yyzz_xxz[k] = -g_y_0_yyz_xxz[k] * cd_z[k] + g_y_0_yyz_xxzz[k];

                g_y_0_yyzz_xyy[k] = -g_y_0_yyz_xyy[k] * cd_z[k] + g_y_0_yyz_xyyz[k];

                g_y_0_yyzz_xyz[k] = -g_y_0_yyz_xyz[k] * cd_z[k] + g_y_0_yyz_xyzz[k];

                g_y_0_yyzz_xzz[k] = -g_y_0_yyz_xzz[k] * cd_z[k] + g_y_0_yyz_xzzz[k];

                g_y_0_yyzz_yyy[k] = -g_y_0_yyz_yyy[k] * cd_z[k] + g_y_0_yyz_yyyz[k];

                g_y_0_yyzz_yyz[k] = -g_y_0_yyz_yyz[k] * cd_z[k] + g_y_0_yyz_yyzz[k];

                g_y_0_yyzz_yzz[k] = -g_y_0_yyz_yzz[k] * cd_z[k] + g_y_0_yyz_yzzz[k];

                g_y_0_yyzz_zzz[k] = -g_y_0_yyz_zzz[k] * cd_z[k] + g_y_0_yyz_zzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_yzz_xxx, g_y_0_yzz_xxxz, g_y_0_yzz_xxy, g_y_0_yzz_xxyz, g_y_0_yzz_xxz, g_y_0_yzz_xxzz, g_y_0_yzz_xyy, g_y_0_yzz_xyyz, g_y_0_yzz_xyz, g_y_0_yzz_xyzz, g_y_0_yzz_xzz, g_y_0_yzz_xzzz, g_y_0_yzz_yyy, g_y_0_yzz_yyyz, g_y_0_yzz_yyz, g_y_0_yzz_yyzz, g_y_0_yzz_yzz, g_y_0_yzz_yzzz, g_y_0_yzz_zzz, g_y_0_yzz_zzzz, g_y_0_yzzz_xxx, g_y_0_yzzz_xxy, g_y_0_yzzz_xxz, g_y_0_yzzz_xyy, g_y_0_yzzz_xyz, g_y_0_yzzz_xzz, g_y_0_yzzz_yyy, g_y_0_yzzz_yyz, g_y_0_yzzz_yzz, g_y_0_yzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_xxx[k] = -g_y_0_yzz_xxx[k] * cd_z[k] + g_y_0_yzz_xxxz[k];

                g_y_0_yzzz_xxy[k] = -g_y_0_yzz_xxy[k] * cd_z[k] + g_y_0_yzz_xxyz[k];

                g_y_0_yzzz_xxz[k] = -g_y_0_yzz_xxz[k] * cd_z[k] + g_y_0_yzz_xxzz[k];

                g_y_0_yzzz_xyy[k] = -g_y_0_yzz_xyy[k] * cd_z[k] + g_y_0_yzz_xyyz[k];

                g_y_0_yzzz_xyz[k] = -g_y_0_yzz_xyz[k] * cd_z[k] + g_y_0_yzz_xyzz[k];

                g_y_0_yzzz_xzz[k] = -g_y_0_yzz_xzz[k] * cd_z[k] + g_y_0_yzz_xzzz[k];

                g_y_0_yzzz_yyy[k] = -g_y_0_yzz_yyy[k] * cd_z[k] + g_y_0_yzz_yyyz[k];

                g_y_0_yzzz_yyz[k] = -g_y_0_yzz_yyz[k] * cd_z[k] + g_y_0_yzz_yyzz[k];

                g_y_0_yzzz_yzz[k] = -g_y_0_yzz_yzz[k] * cd_z[k] + g_y_0_yzz_yzzz[k];

                g_y_0_yzzz_zzz[k] = -g_y_0_yzz_zzz[k] * cd_z[k] + g_y_0_yzz_zzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_zzz_xxx, g_y_0_zzz_xxxz, g_y_0_zzz_xxy, g_y_0_zzz_xxyz, g_y_0_zzz_xxz, g_y_0_zzz_xxzz, g_y_0_zzz_xyy, g_y_0_zzz_xyyz, g_y_0_zzz_xyz, g_y_0_zzz_xyzz, g_y_0_zzz_xzz, g_y_0_zzz_xzzz, g_y_0_zzz_yyy, g_y_0_zzz_yyyz, g_y_0_zzz_yyz, g_y_0_zzz_yyzz, g_y_0_zzz_yzz, g_y_0_zzz_yzzz, g_y_0_zzz_zzz, g_y_0_zzz_zzzz, g_y_0_zzzz_xxx, g_y_0_zzzz_xxy, g_y_0_zzzz_xxz, g_y_0_zzzz_xyy, g_y_0_zzzz_xyz, g_y_0_zzzz_xzz, g_y_0_zzzz_yyy, g_y_0_zzzz_yyz, g_y_0_zzzz_yzz, g_y_0_zzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_xxx[k] = -g_y_0_zzz_xxx[k] * cd_z[k] + g_y_0_zzz_xxxz[k];

                g_y_0_zzzz_xxy[k] = -g_y_0_zzz_xxy[k] * cd_z[k] + g_y_0_zzz_xxyz[k];

                g_y_0_zzzz_xxz[k] = -g_y_0_zzz_xxz[k] * cd_z[k] + g_y_0_zzz_xxzz[k];

                g_y_0_zzzz_xyy[k] = -g_y_0_zzz_xyy[k] * cd_z[k] + g_y_0_zzz_xyyz[k];

                g_y_0_zzzz_xyz[k] = -g_y_0_zzz_xyz[k] * cd_z[k] + g_y_0_zzz_xyzz[k];

                g_y_0_zzzz_xzz[k] = -g_y_0_zzz_xzz[k] * cd_z[k] + g_y_0_zzz_xzzz[k];

                g_y_0_zzzz_yyy[k] = -g_y_0_zzz_yyy[k] * cd_z[k] + g_y_0_zzz_yyyz[k];

                g_y_0_zzzz_yyz[k] = -g_y_0_zzz_yyz[k] * cd_z[k] + g_y_0_zzz_yyzz[k];

                g_y_0_zzzz_yzz[k] = -g_y_0_zzz_yzz[k] * cd_z[k] + g_y_0_zzz_yzzz[k];

                g_y_0_zzzz_zzz[k] = -g_y_0_zzz_zzz[k] * cd_z[k] + g_y_0_zzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxx_xxx, g_z_0_xxx_xxxx, g_z_0_xxx_xxxy, g_z_0_xxx_xxxz, g_z_0_xxx_xxy, g_z_0_xxx_xxyy, g_z_0_xxx_xxyz, g_z_0_xxx_xxz, g_z_0_xxx_xxzz, g_z_0_xxx_xyy, g_z_0_xxx_xyyy, g_z_0_xxx_xyyz, g_z_0_xxx_xyz, g_z_0_xxx_xyzz, g_z_0_xxx_xzz, g_z_0_xxx_xzzz, g_z_0_xxx_yyy, g_z_0_xxx_yyz, g_z_0_xxx_yzz, g_z_0_xxx_zzz, g_z_0_xxxx_xxx, g_z_0_xxxx_xxy, g_z_0_xxxx_xxz, g_z_0_xxxx_xyy, g_z_0_xxxx_xyz, g_z_0_xxxx_xzz, g_z_0_xxxx_yyy, g_z_0_xxxx_yyz, g_z_0_xxxx_yzz, g_z_0_xxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_xxx[k] = -g_z_0_xxx_xxx[k] * cd_x[k] + g_z_0_xxx_xxxx[k];

                g_z_0_xxxx_xxy[k] = -g_z_0_xxx_xxy[k] * cd_x[k] + g_z_0_xxx_xxxy[k];

                g_z_0_xxxx_xxz[k] = -g_z_0_xxx_xxz[k] * cd_x[k] + g_z_0_xxx_xxxz[k];

                g_z_0_xxxx_xyy[k] = -g_z_0_xxx_xyy[k] * cd_x[k] + g_z_0_xxx_xxyy[k];

                g_z_0_xxxx_xyz[k] = -g_z_0_xxx_xyz[k] * cd_x[k] + g_z_0_xxx_xxyz[k];

                g_z_0_xxxx_xzz[k] = -g_z_0_xxx_xzz[k] * cd_x[k] + g_z_0_xxx_xxzz[k];

                g_z_0_xxxx_yyy[k] = -g_z_0_xxx_yyy[k] * cd_x[k] + g_z_0_xxx_xyyy[k];

                g_z_0_xxxx_yyz[k] = -g_z_0_xxx_yyz[k] * cd_x[k] + g_z_0_xxx_xyyz[k];

                g_z_0_xxxx_yzz[k] = -g_z_0_xxx_yzz[k] * cd_x[k] + g_z_0_xxx_xyzz[k];

                g_z_0_xxxx_zzz[k] = -g_z_0_xxx_zzz[k] * cd_x[k] + g_z_0_xxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxxy_xxx, g_z_0_xxxy_xxy, g_z_0_xxxy_xxz, g_z_0_xxxy_xyy, g_z_0_xxxy_xyz, g_z_0_xxxy_xzz, g_z_0_xxxy_yyy, g_z_0_xxxy_yyz, g_z_0_xxxy_yzz, g_z_0_xxxy_zzz, g_z_0_xxy_xxx, g_z_0_xxy_xxxx, g_z_0_xxy_xxxy, g_z_0_xxy_xxxz, g_z_0_xxy_xxy, g_z_0_xxy_xxyy, g_z_0_xxy_xxyz, g_z_0_xxy_xxz, g_z_0_xxy_xxzz, g_z_0_xxy_xyy, g_z_0_xxy_xyyy, g_z_0_xxy_xyyz, g_z_0_xxy_xyz, g_z_0_xxy_xyzz, g_z_0_xxy_xzz, g_z_0_xxy_xzzz, g_z_0_xxy_yyy, g_z_0_xxy_yyz, g_z_0_xxy_yzz, g_z_0_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_xxx[k] = -g_z_0_xxy_xxx[k] * cd_x[k] + g_z_0_xxy_xxxx[k];

                g_z_0_xxxy_xxy[k] = -g_z_0_xxy_xxy[k] * cd_x[k] + g_z_0_xxy_xxxy[k];

                g_z_0_xxxy_xxz[k] = -g_z_0_xxy_xxz[k] * cd_x[k] + g_z_0_xxy_xxxz[k];

                g_z_0_xxxy_xyy[k] = -g_z_0_xxy_xyy[k] * cd_x[k] + g_z_0_xxy_xxyy[k];

                g_z_0_xxxy_xyz[k] = -g_z_0_xxy_xyz[k] * cd_x[k] + g_z_0_xxy_xxyz[k];

                g_z_0_xxxy_xzz[k] = -g_z_0_xxy_xzz[k] * cd_x[k] + g_z_0_xxy_xxzz[k];

                g_z_0_xxxy_yyy[k] = -g_z_0_xxy_yyy[k] * cd_x[k] + g_z_0_xxy_xyyy[k];

                g_z_0_xxxy_yyz[k] = -g_z_0_xxy_yyz[k] * cd_x[k] + g_z_0_xxy_xyyz[k];

                g_z_0_xxxy_yzz[k] = -g_z_0_xxy_yzz[k] * cd_x[k] + g_z_0_xxy_xyzz[k];

                g_z_0_xxxy_zzz[k] = -g_z_0_xxy_zzz[k] * cd_x[k] + g_z_0_xxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxxz_xxx, g_z_0_xxxz_xxy, g_z_0_xxxz_xxz, g_z_0_xxxz_xyy, g_z_0_xxxz_xyz, g_z_0_xxxz_xzz, g_z_0_xxxz_yyy, g_z_0_xxxz_yyz, g_z_0_xxxz_yzz, g_z_0_xxxz_zzz, g_z_0_xxz_xxx, g_z_0_xxz_xxxx, g_z_0_xxz_xxxy, g_z_0_xxz_xxxz, g_z_0_xxz_xxy, g_z_0_xxz_xxyy, g_z_0_xxz_xxyz, g_z_0_xxz_xxz, g_z_0_xxz_xxzz, g_z_0_xxz_xyy, g_z_0_xxz_xyyy, g_z_0_xxz_xyyz, g_z_0_xxz_xyz, g_z_0_xxz_xyzz, g_z_0_xxz_xzz, g_z_0_xxz_xzzz, g_z_0_xxz_yyy, g_z_0_xxz_yyz, g_z_0_xxz_yzz, g_z_0_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_xxx[k] = -g_z_0_xxz_xxx[k] * cd_x[k] + g_z_0_xxz_xxxx[k];

                g_z_0_xxxz_xxy[k] = -g_z_0_xxz_xxy[k] * cd_x[k] + g_z_0_xxz_xxxy[k];

                g_z_0_xxxz_xxz[k] = -g_z_0_xxz_xxz[k] * cd_x[k] + g_z_0_xxz_xxxz[k];

                g_z_0_xxxz_xyy[k] = -g_z_0_xxz_xyy[k] * cd_x[k] + g_z_0_xxz_xxyy[k];

                g_z_0_xxxz_xyz[k] = -g_z_0_xxz_xyz[k] * cd_x[k] + g_z_0_xxz_xxyz[k];

                g_z_0_xxxz_xzz[k] = -g_z_0_xxz_xzz[k] * cd_x[k] + g_z_0_xxz_xxzz[k];

                g_z_0_xxxz_yyy[k] = -g_z_0_xxz_yyy[k] * cd_x[k] + g_z_0_xxz_xyyy[k];

                g_z_0_xxxz_yyz[k] = -g_z_0_xxz_yyz[k] * cd_x[k] + g_z_0_xxz_xyyz[k];

                g_z_0_xxxz_yzz[k] = -g_z_0_xxz_yzz[k] * cd_x[k] + g_z_0_xxz_xyzz[k];

                g_z_0_xxxz_zzz[k] = -g_z_0_xxz_zzz[k] * cd_x[k] + g_z_0_xxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxyy_xxx, g_z_0_xxyy_xxy, g_z_0_xxyy_xxz, g_z_0_xxyy_xyy, g_z_0_xxyy_xyz, g_z_0_xxyy_xzz, g_z_0_xxyy_yyy, g_z_0_xxyy_yyz, g_z_0_xxyy_yzz, g_z_0_xxyy_zzz, g_z_0_xyy_xxx, g_z_0_xyy_xxxx, g_z_0_xyy_xxxy, g_z_0_xyy_xxxz, g_z_0_xyy_xxy, g_z_0_xyy_xxyy, g_z_0_xyy_xxyz, g_z_0_xyy_xxz, g_z_0_xyy_xxzz, g_z_0_xyy_xyy, g_z_0_xyy_xyyy, g_z_0_xyy_xyyz, g_z_0_xyy_xyz, g_z_0_xyy_xyzz, g_z_0_xyy_xzz, g_z_0_xyy_xzzz, g_z_0_xyy_yyy, g_z_0_xyy_yyz, g_z_0_xyy_yzz, g_z_0_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_xxx[k] = -g_z_0_xyy_xxx[k] * cd_x[k] + g_z_0_xyy_xxxx[k];

                g_z_0_xxyy_xxy[k] = -g_z_0_xyy_xxy[k] * cd_x[k] + g_z_0_xyy_xxxy[k];

                g_z_0_xxyy_xxz[k] = -g_z_0_xyy_xxz[k] * cd_x[k] + g_z_0_xyy_xxxz[k];

                g_z_0_xxyy_xyy[k] = -g_z_0_xyy_xyy[k] * cd_x[k] + g_z_0_xyy_xxyy[k];

                g_z_0_xxyy_xyz[k] = -g_z_0_xyy_xyz[k] * cd_x[k] + g_z_0_xyy_xxyz[k];

                g_z_0_xxyy_xzz[k] = -g_z_0_xyy_xzz[k] * cd_x[k] + g_z_0_xyy_xxzz[k];

                g_z_0_xxyy_yyy[k] = -g_z_0_xyy_yyy[k] * cd_x[k] + g_z_0_xyy_xyyy[k];

                g_z_0_xxyy_yyz[k] = -g_z_0_xyy_yyz[k] * cd_x[k] + g_z_0_xyy_xyyz[k];

                g_z_0_xxyy_yzz[k] = -g_z_0_xyy_yzz[k] * cd_x[k] + g_z_0_xyy_xyzz[k];

                g_z_0_xxyy_zzz[k] = -g_z_0_xyy_zzz[k] * cd_x[k] + g_z_0_xyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxyz_xxx, g_z_0_xxyz_xxy, g_z_0_xxyz_xxz, g_z_0_xxyz_xyy, g_z_0_xxyz_xyz, g_z_0_xxyz_xzz, g_z_0_xxyz_yyy, g_z_0_xxyz_yyz, g_z_0_xxyz_yzz, g_z_0_xxyz_zzz, g_z_0_xyz_xxx, g_z_0_xyz_xxxx, g_z_0_xyz_xxxy, g_z_0_xyz_xxxz, g_z_0_xyz_xxy, g_z_0_xyz_xxyy, g_z_0_xyz_xxyz, g_z_0_xyz_xxz, g_z_0_xyz_xxzz, g_z_0_xyz_xyy, g_z_0_xyz_xyyy, g_z_0_xyz_xyyz, g_z_0_xyz_xyz, g_z_0_xyz_xyzz, g_z_0_xyz_xzz, g_z_0_xyz_xzzz, g_z_0_xyz_yyy, g_z_0_xyz_yyz, g_z_0_xyz_yzz, g_z_0_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_xxx[k] = -g_z_0_xyz_xxx[k] * cd_x[k] + g_z_0_xyz_xxxx[k];

                g_z_0_xxyz_xxy[k] = -g_z_0_xyz_xxy[k] * cd_x[k] + g_z_0_xyz_xxxy[k];

                g_z_0_xxyz_xxz[k] = -g_z_0_xyz_xxz[k] * cd_x[k] + g_z_0_xyz_xxxz[k];

                g_z_0_xxyz_xyy[k] = -g_z_0_xyz_xyy[k] * cd_x[k] + g_z_0_xyz_xxyy[k];

                g_z_0_xxyz_xyz[k] = -g_z_0_xyz_xyz[k] * cd_x[k] + g_z_0_xyz_xxyz[k];

                g_z_0_xxyz_xzz[k] = -g_z_0_xyz_xzz[k] * cd_x[k] + g_z_0_xyz_xxzz[k];

                g_z_0_xxyz_yyy[k] = -g_z_0_xyz_yyy[k] * cd_x[k] + g_z_0_xyz_xyyy[k];

                g_z_0_xxyz_yyz[k] = -g_z_0_xyz_yyz[k] * cd_x[k] + g_z_0_xyz_xyyz[k];

                g_z_0_xxyz_yzz[k] = -g_z_0_xyz_yzz[k] * cd_x[k] + g_z_0_xyz_xyzz[k];

                g_z_0_xxyz_zzz[k] = -g_z_0_xyz_zzz[k] * cd_x[k] + g_z_0_xyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxzz_xxx, g_z_0_xxzz_xxy, g_z_0_xxzz_xxz, g_z_0_xxzz_xyy, g_z_0_xxzz_xyz, g_z_0_xxzz_xzz, g_z_0_xxzz_yyy, g_z_0_xxzz_yyz, g_z_0_xxzz_yzz, g_z_0_xxzz_zzz, g_z_0_xzz_xxx, g_z_0_xzz_xxxx, g_z_0_xzz_xxxy, g_z_0_xzz_xxxz, g_z_0_xzz_xxy, g_z_0_xzz_xxyy, g_z_0_xzz_xxyz, g_z_0_xzz_xxz, g_z_0_xzz_xxzz, g_z_0_xzz_xyy, g_z_0_xzz_xyyy, g_z_0_xzz_xyyz, g_z_0_xzz_xyz, g_z_0_xzz_xyzz, g_z_0_xzz_xzz, g_z_0_xzz_xzzz, g_z_0_xzz_yyy, g_z_0_xzz_yyz, g_z_0_xzz_yzz, g_z_0_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_xxx[k] = -g_z_0_xzz_xxx[k] * cd_x[k] + g_z_0_xzz_xxxx[k];

                g_z_0_xxzz_xxy[k] = -g_z_0_xzz_xxy[k] * cd_x[k] + g_z_0_xzz_xxxy[k];

                g_z_0_xxzz_xxz[k] = -g_z_0_xzz_xxz[k] * cd_x[k] + g_z_0_xzz_xxxz[k];

                g_z_0_xxzz_xyy[k] = -g_z_0_xzz_xyy[k] * cd_x[k] + g_z_0_xzz_xxyy[k];

                g_z_0_xxzz_xyz[k] = -g_z_0_xzz_xyz[k] * cd_x[k] + g_z_0_xzz_xxyz[k];

                g_z_0_xxzz_xzz[k] = -g_z_0_xzz_xzz[k] * cd_x[k] + g_z_0_xzz_xxzz[k];

                g_z_0_xxzz_yyy[k] = -g_z_0_xzz_yyy[k] * cd_x[k] + g_z_0_xzz_xyyy[k];

                g_z_0_xxzz_yyz[k] = -g_z_0_xzz_yyz[k] * cd_x[k] + g_z_0_xzz_xyyz[k];

                g_z_0_xxzz_yzz[k] = -g_z_0_xzz_yzz[k] * cd_x[k] + g_z_0_xzz_xyzz[k];

                g_z_0_xxzz_zzz[k] = -g_z_0_xzz_zzz[k] * cd_x[k] + g_z_0_xzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xyyy_xxx, g_z_0_xyyy_xxy, g_z_0_xyyy_xxz, g_z_0_xyyy_xyy, g_z_0_xyyy_xyz, g_z_0_xyyy_xzz, g_z_0_xyyy_yyy, g_z_0_xyyy_yyz, g_z_0_xyyy_yzz, g_z_0_xyyy_zzz, g_z_0_yyy_xxx, g_z_0_yyy_xxxx, g_z_0_yyy_xxxy, g_z_0_yyy_xxxz, g_z_0_yyy_xxy, g_z_0_yyy_xxyy, g_z_0_yyy_xxyz, g_z_0_yyy_xxz, g_z_0_yyy_xxzz, g_z_0_yyy_xyy, g_z_0_yyy_xyyy, g_z_0_yyy_xyyz, g_z_0_yyy_xyz, g_z_0_yyy_xyzz, g_z_0_yyy_xzz, g_z_0_yyy_xzzz, g_z_0_yyy_yyy, g_z_0_yyy_yyz, g_z_0_yyy_yzz, g_z_0_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_xxx[k] = -g_z_0_yyy_xxx[k] * cd_x[k] + g_z_0_yyy_xxxx[k];

                g_z_0_xyyy_xxy[k] = -g_z_0_yyy_xxy[k] * cd_x[k] + g_z_0_yyy_xxxy[k];

                g_z_0_xyyy_xxz[k] = -g_z_0_yyy_xxz[k] * cd_x[k] + g_z_0_yyy_xxxz[k];

                g_z_0_xyyy_xyy[k] = -g_z_0_yyy_xyy[k] * cd_x[k] + g_z_0_yyy_xxyy[k];

                g_z_0_xyyy_xyz[k] = -g_z_0_yyy_xyz[k] * cd_x[k] + g_z_0_yyy_xxyz[k];

                g_z_0_xyyy_xzz[k] = -g_z_0_yyy_xzz[k] * cd_x[k] + g_z_0_yyy_xxzz[k];

                g_z_0_xyyy_yyy[k] = -g_z_0_yyy_yyy[k] * cd_x[k] + g_z_0_yyy_xyyy[k];

                g_z_0_xyyy_yyz[k] = -g_z_0_yyy_yyz[k] * cd_x[k] + g_z_0_yyy_xyyz[k];

                g_z_0_xyyy_yzz[k] = -g_z_0_yyy_yzz[k] * cd_x[k] + g_z_0_yyy_xyzz[k];

                g_z_0_xyyy_zzz[k] = -g_z_0_yyy_zzz[k] * cd_x[k] + g_z_0_yyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xyyz_xxx, g_z_0_xyyz_xxy, g_z_0_xyyz_xxz, g_z_0_xyyz_xyy, g_z_0_xyyz_xyz, g_z_0_xyyz_xzz, g_z_0_xyyz_yyy, g_z_0_xyyz_yyz, g_z_0_xyyz_yzz, g_z_0_xyyz_zzz, g_z_0_yyz_xxx, g_z_0_yyz_xxxx, g_z_0_yyz_xxxy, g_z_0_yyz_xxxz, g_z_0_yyz_xxy, g_z_0_yyz_xxyy, g_z_0_yyz_xxyz, g_z_0_yyz_xxz, g_z_0_yyz_xxzz, g_z_0_yyz_xyy, g_z_0_yyz_xyyy, g_z_0_yyz_xyyz, g_z_0_yyz_xyz, g_z_0_yyz_xyzz, g_z_0_yyz_xzz, g_z_0_yyz_xzzz, g_z_0_yyz_yyy, g_z_0_yyz_yyz, g_z_0_yyz_yzz, g_z_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_xxx[k] = -g_z_0_yyz_xxx[k] * cd_x[k] + g_z_0_yyz_xxxx[k];

                g_z_0_xyyz_xxy[k] = -g_z_0_yyz_xxy[k] * cd_x[k] + g_z_0_yyz_xxxy[k];

                g_z_0_xyyz_xxz[k] = -g_z_0_yyz_xxz[k] * cd_x[k] + g_z_0_yyz_xxxz[k];

                g_z_0_xyyz_xyy[k] = -g_z_0_yyz_xyy[k] * cd_x[k] + g_z_0_yyz_xxyy[k];

                g_z_0_xyyz_xyz[k] = -g_z_0_yyz_xyz[k] * cd_x[k] + g_z_0_yyz_xxyz[k];

                g_z_0_xyyz_xzz[k] = -g_z_0_yyz_xzz[k] * cd_x[k] + g_z_0_yyz_xxzz[k];

                g_z_0_xyyz_yyy[k] = -g_z_0_yyz_yyy[k] * cd_x[k] + g_z_0_yyz_xyyy[k];

                g_z_0_xyyz_yyz[k] = -g_z_0_yyz_yyz[k] * cd_x[k] + g_z_0_yyz_xyyz[k];

                g_z_0_xyyz_yzz[k] = -g_z_0_yyz_yzz[k] * cd_x[k] + g_z_0_yyz_xyzz[k];

                g_z_0_xyyz_zzz[k] = -g_z_0_yyz_zzz[k] * cd_x[k] + g_z_0_yyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xyzz_xxx, g_z_0_xyzz_xxy, g_z_0_xyzz_xxz, g_z_0_xyzz_xyy, g_z_0_xyzz_xyz, g_z_0_xyzz_xzz, g_z_0_xyzz_yyy, g_z_0_xyzz_yyz, g_z_0_xyzz_yzz, g_z_0_xyzz_zzz, g_z_0_yzz_xxx, g_z_0_yzz_xxxx, g_z_0_yzz_xxxy, g_z_0_yzz_xxxz, g_z_0_yzz_xxy, g_z_0_yzz_xxyy, g_z_0_yzz_xxyz, g_z_0_yzz_xxz, g_z_0_yzz_xxzz, g_z_0_yzz_xyy, g_z_0_yzz_xyyy, g_z_0_yzz_xyyz, g_z_0_yzz_xyz, g_z_0_yzz_xyzz, g_z_0_yzz_xzz, g_z_0_yzz_xzzz, g_z_0_yzz_yyy, g_z_0_yzz_yyz, g_z_0_yzz_yzz, g_z_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_xxx[k] = -g_z_0_yzz_xxx[k] * cd_x[k] + g_z_0_yzz_xxxx[k];

                g_z_0_xyzz_xxy[k] = -g_z_0_yzz_xxy[k] * cd_x[k] + g_z_0_yzz_xxxy[k];

                g_z_0_xyzz_xxz[k] = -g_z_0_yzz_xxz[k] * cd_x[k] + g_z_0_yzz_xxxz[k];

                g_z_0_xyzz_xyy[k] = -g_z_0_yzz_xyy[k] * cd_x[k] + g_z_0_yzz_xxyy[k];

                g_z_0_xyzz_xyz[k] = -g_z_0_yzz_xyz[k] * cd_x[k] + g_z_0_yzz_xxyz[k];

                g_z_0_xyzz_xzz[k] = -g_z_0_yzz_xzz[k] * cd_x[k] + g_z_0_yzz_xxzz[k];

                g_z_0_xyzz_yyy[k] = -g_z_0_yzz_yyy[k] * cd_x[k] + g_z_0_yzz_xyyy[k];

                g_z_0_xyzz_yyz[k] = -g_z_0_yzz_yyz[k] * cd_x[k] + g_z_0_yzz_xyyz[k];

                g_z_0_xyzz_yzz[k] = -g_z_0_yzz_yzz[k] * cd_x[k] + g_z_0_yzz_xyzz[k];

                g_z_0_xyzz_zzz[k] = -g_z_0_yzz_zzz[k] * cd_x[k] + g_z_0_yzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xzzz_xxx, g_z_0_xzzz_xxy, g_z_0_xzzz_xxz, g_z_0_xzzz_xyy, g_z_0_xzzz_xyz, g_z_0_xzzz_xzz, g_z_0_xzzz_yyy, g_z_0_xzzz_yyz, g_z_0_xzzz_yzz, g_z_0_xzzz_zzz, g_z_0_zzz_xxx, g_z_0_zzz_xxxx, g_z_0_zzz_xxxy, g_z_0_zzz_xxxz, g_z_0_zzz_xxy, g_z_0_zzz_xxyy, g_z_0_zzz_xxyz, g_z_0_zzz_xxz, g_z_0_zzz_xxzz, g_z_0_zzz_xyy, g_z_0_zzz_xyyy, g_z_0_zzz_xyyz, g_z_0_zzz_xyz, g_z_0_zzz_xyzz, g_z_0_zzz_xzz, g_z_0_zzz_xzzz, g_z_0_zzz_yyy, g_z_0_zzz_yyz, g_z_0_zzz_yzz, g_z_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_xxx[k] = -g_z_0_zzz_xxx[k] * cd_x[k] + g_z_0_zzz_xxxx[k];

                g_z_0_xzzz_xxy[k] = -g_z_0_zzz_xxy[k] * cd_x[k] + g_z_0_zzz_xxxy[k];

                g_z_0_xzzz_xxz[k] = -g_z_0_zzz_xxz[k] * cd_x[k] + g_z_0_zzz_xxxz[k];

                g_z_0_xzzz_xyy[k] = -g_z_0_zzz_xyy[k] * cd_x[k] + g_z_0_zzz_xxyy[k];

                g_z_0_xzzz_xyz[k] = -g_z_0_zzz_xyz[k] * cd_x[k] + g_z_0_zzz_xxyz[k];

                g_z_0_xzzz_xzz[k] = -g_z_0_zzz_xzz[k] * cd_x[k] + g_z_0_zzz_xxzz[k];

                g_z_0_xzzz_yyy[k] = -g_z_0_zzz_yyy[k] * cd_x[k] + g_z_0_zzz_xyyy[k];

                g_z_0_xzzz_yyz[k] = -g_z_0_zzz_yyz[k] * cd_x[k] + g_z_0_zzz_xyyz[k];

                g_z_0_xzzz_yzz[k] = -g_z_0_zzz_yzz[k] * cd_x[k] + g_z_0_zzz_xyzz[k];

                g_z_0_xzzz_zzz[k] = -g_z_0_zzz_zzz[k] * cd_x[k] + g_z_0_zzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yyy_xxx, g_z_0_yyy_xxxy, g_z_0_yyy_xxy, g_z_0_yyy_xxyy, g_z_0_yyy_xxyz, g_z_0_yyy_xxz, g_z_0_yyy_xyy, g_z_0_yyy_xyyy, g_z_0_yyy_xyyz, g_z_0_yyy_xyz, g_z_0_yyy_xyzz, g_z_0_yyy_xzz, g_z_0_yyy_yyy, g_z_0_yyy_yyyy, g_z_0_yyy_yyyz, g_z_0_yyy_yyz, g_z_0_yyy_yyzz, g_z_0_yyy_yzz, g_z_0_yyy_yzzz, g_z_0_yyy_zzz, g_z_0_yyyy_xxx, g_z_0_yyyy_xxy, g_z_0_yyyy_xxz, g_z_0_yyyy_xyy, g_z_0_yyyy_xyz, g_z_0_yyyy_xzz, g_z_0_yyyy_yyy, g_z_0_yyyy_yyz, g_z_0_yyyy_yzz, g_z_0_yyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_xxx[k] = -g_z_0_yyy_xxx[k] * cd_y[k] + g_z_0_yyy_xxxy[k];

                g_z_0_yyyy_xxy[k] = -g_z_0_yyy_xxy[k] * cd_y[k] + g_z_0_yyy_xxyy[k];

                g_z_0_yyyy_xxz[k] = -g_z_0_yyy_xxz[k] * cd_y[k] + g_z_0_yyy_xxyz[k];

                g_z_0_yyyy_xyy[k] = -g_z_0_yyy_xyy[k] * cd_y[k] + g_z_0_yyy_xyyy[k];

                g_z_0_yyyy_xyz[k] = -g_z_0_yyy_xyz[k] * cd_y[k] + g_z_0_yyy_xyyz[k];

                g_z_0_yyyy_xzz[k] = -g_z_0_yyy_xzz[k] * cd_y[k] + g_z_0_yyy_xyzz[k];

                g_z_0_yyyy_yyy[k] = -g_z_0_yyy_yyy[k] * cd_y[k] + g_z_0_yyy_yyyy[k];

                g_z_0_yyyy_yyz[k] = -g_z_0_yyy_yyz[k] * cd_y[k] + g_z_0_yyy_yyyz[k];

                g_z_0_yyyy_yzz[k] = -g_z_0_yyy_yzz[k] * cd_y[k] + g_z_0_yyy_yyzz[k];

                g_z_0_yyyy_zzz[k] = -g_z_0_yyy_zzz[k] * cd_y[k] + g_z_0_yyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yyyz_xxx, g_z_0_yyyz_xxy, g_z_0_yyyz_xxz, g_z_0_yyyz_xyy, g_z_0_yyyz_xyz, g_z_0_yyyz_xzz, g_z_0_yyyz_yyy, g_z_0_yyyz_yyz, g_z_0_yyyz_yzz, g_z_0_yyyz_zzz, g_z_0_yyz_xxx, g_z_0_yyz_xxxy, g_z_0_yyz_xxy, g_z_0_yyz_xxyy, g_z_0_yyz_xxyz, g_z_0_yyz_xxz, g_z_0_yyz_xyy, g_z_0_yyz_xyyy, g_z_0_yyz_xyyz, g_z_0_yyz_xyz, g_z_0_yyz_xyzz, g_z_0_yyz_xzz, g_z_0_yyz_yyy, g_z_0_yyz_yyyy, g_z_0_yyz_yyyz, g_z_0_yyz_yyz, g_z_0_yyz_yyzz, g_z_0_yyz_yzz, g_z_0_yyz_yzzz, g_z_0_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_xxx[k] = -g_z_0_yyz_xxx[k] * cd_y[k] + g_z_0_yyz_xxxy[k];

                g_z_0_yyyz_xxy[k] = -g_z_0_yyz_xxy[k] * cd_y[k] + g_z_0_yyz_xxyy[k];

                g_z_0_yyyz_xxz[k] = -g_z_0_yyz_xxz[k] * cd_y[k] + g_z_0_yyz_xxyz[k];

                g_z_0_yyyz_xyy[k] = -g_z_0_yyz_xyy[k] * cd_y[k] + g_z_0_yyz_xyyy[k];

                g_z_0_yyyz_xyz[k] = -g_z_0_yyz_xyz[k] * cd_y[k] + g_z_0_yyz_xyyz[k];

                g_z_0_yyyz_xzz[k] = -g_z_0_yyz_xzz[k] * cd_y[k] + g_z_0_yyz_xyzz[k];

                g_z_0_yyyz_yyy[k] = -g_z_0_yyz_yyy[k] * cd_y[k] + g_z_0_yyz_yyyy[k];

                g_z_0_yyyz_yyz[k] = -g_z_0_yyz_yyz[k] * cd_y[k] + g_z_0_yyz_yyyz[k];

                g_z_0_yyyz_yzz[k] = -g_z_0_yyz_yzz[k] * cd_y[k] + g_z_0_yyz_yyzz[k];

                g_z_0_yyyz_zzz[k] = -g_z_0_yyz_zzz[k] * cd_y[k] + g_z_0_yyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yyzz_xxx, g_z_0_yyzz_xxy, g_z_0_yyzz_xxz, g_z_0_yyzz_xyy, g_z_0_yyzz_xyz, g_z_0_yyzz_xzz, g_z_0_yyzz_yyy, g_z_0_yyzz_yyz, g_z_0_yyzz_yzz, g_z_0_yyzz_zzz, g_z_0_yzz_xxx, g_z_0_yzz_xxxy, g_z_0_yzz_xxy, g_z_0_yzz_xxyy, g_z_0_yzz_xxyz, g_z_0_yzz_xxz, g_z_0_yzz_xyy, g_z_0_yzz_xyyy, g_z_0_yzz_xyyz, g_z_0_yzz_xyz, g_z_0_yzz_xyzz, g_z_0_yzz_xzz, g_z_0_yzz_yyy, g_z_0_yzz_yyyy, g_z_0_yzz_yyyz, g_z_0_yzz_yyz, g_z_0_yzz_yyzz, g_z_0_yzz_yzz, g_z_0_yzz_yzzz, g_z_0_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_xxx[k] = -g_z_0_yzz_xxx[k] * cd_y[k] + g_z_0_yzz_xxxy[k];

                g_z_0_yyzz_xxy[k] = -g_z_0_yzz_xxy[k] * cd_y[k] + g_z_0_yzz_xxyy[k];

                g_z_0_yyzz_xxz[k] = -g_z_0_yzz_xxz[k] * cd_y[k] + g_z_0_yzz_xxyz[k];

                g_z_0_yyzz_xyy[k] = -g_z_0_yzz_xyy[k] * cd_y[k] + g_z_0_yzz_xyyy[k];

                g_z_0_yyzz_xyz[k] = -g_z_0_yzz_xyz[k] * cd_y[k] + g_z_0_yzz_xyyz[k];

                g_z_0_yyzz_xzz[k] = -g_z_0_yzz_xzz[k] * cd_y[k] + g_z_0_yzz_xyzz[k];

                g_z_0_yyzz_yyy[k] = -g_z_0_yzz_yyy[k] * cd_y[k] + g_z_0_yzz_yyyy[k];

                g_z_0_yyzz_yyz[k] = -g_z_0_yzz_yyz[k] * cd_y[k] + g_z_0_yzz_yyyz[k];

                g_z_0_yyzz_yzz[k] = -g_z_0_yzz_yzz[k] * cd_y[k] + g_z_0_yzz_yyzz[k];

                g_z_0_yyzz_zzz[k] = -g_z_0_yzz_zzz[k] * cd_y[k] + g_z_0_yzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yzzz_xxx, g_z_0_yzzz_xxy, g_z_0_yzzz_xxz, g_z_0_yzzz_xyy, g_z_0_yzzz_xyz, g_z_0_yzzz_xzz, g_z_0_yzzz_yyy, g_z_0_yzzz_yyz, g_z_0_yzzz_yzz, g_z_0_yzzz_zzz, g_z_0_zzz_xxx, g_z_0_zzz_xxxy, g_z_0_zzz_xxy, g_z_0_zzz_xxyy, g_z_0_zzz_xxyz, g_z_0_zzz_xxz, g_z_0_zzz_xyy, g_z_0_zzz_xyyy, g_z_0_zzz_xyyz, g_z_0_zzz_xyz, g_z_0_zzz_xyzz, g_z_0_zzz_xzz, g_z_0_zzz_yyy, g_z_0_zzz_yyyy, g_z_0_zzz_yyyz, g_z_0_zzz_yyz, g_z_0_zzz_yyzz, g_z_0_zzz_yzz, g_z_0_zzz_yzzz, g_z_0_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_xxx[k] = -g_z_0_zzz_xxx[k] * cd_y[k] + g_z_0_zzz_xxxy[k];

                g_z_0_yzzz_xxy[k] = -g_z_0_zzz_xxy[k] * cd_y[k] + g_z_0_zzz_xxyy[k];

                g_z_0_yzzz_xxz[k] = -g_z_0_zzz_xxz[k] * cd_y[k] + g_z_0_zzz_xxyz[k];

                g_z_0_yzzz_xyy[k] = -g_z_0_zzz_xyy[k] * cd_y[k] + g_z_0_zzz_xyyy[k];

                g_z_0_yzzz_xyz[k] = -g_z_0_zzz_xyz[k] * cd_y[k] + g_z_0_zzz_xyyz[k];

                g_z_0_yzzz_xzz[k] = -g_z_0_zzz_xzz[k] * cd_y[k] + g_z_0_zzz_xyzz[k];

                g_z_0_yzzz_yyy[k] = -g_z_0_zzz_yyy[k] * cd_y[k] + g_z_0_zzz_yyyy[k];

                g_z_0_yzzz_yyz[k] = -g_z_0_zzz_yyz[k] * cd_y[k] + g_z_0_zzz_yyyz[k];

                g_z_0_yzzz_yzz[k] = -g_z_0_zzz_yzz[k] * cd_y[k] + g_z_0_zzz_yyzz[k];

                g_z_0_yzzz_zzz[k] = -g_z_0_zzz_zzz[k] * cd_y[k] + g_z_0_zzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_z_0_zzz_xxx, g_z_0_zzz_xxxz, g_z_0_zzz_xxy, g_z_0_zzz_xxyz, g_z_0_zzz_xxz, g_z_0_zzz_xxzz, g_z_0_zzz_xyy, g_z_0_zzz_xyyz, g_z_0_zzz_xyz, g_z_0_zzz_xyzz, g_z_0_zzz_xzz, g_z_0_zzz_xzzz, g_z_0_zzz_yyy, g_z_0_zzz_yyyz, g_z_0_zzz_yyz, g_z_0_zzz_yyzz, g_z_0_zzz_yzz, g_z_0_zzz_yzzz, g_z_0_zzz_zzz, g_z_0_zzz_zzzz, g_z_0_zzzz_xxx, g_z_0_zzzz_xxy, g_z_0_zzzz_xxz, g_z_0_zzzz_xyy, g_z_0_zzzz_xyz, g_z_0_zzzz_xzz, g_z_0_zzzz_yyy, g_z_0_zzzz_yyz, g_z_0_zzzz_yzz, g_z_0_zzzz_zzz, g_zzz_xxx, g_zzz_xxy, g_zzz_xxz, g_zzz_xyy, g_zzz_xyz, g_zzz_xzz, g_zzz_yyy, g_zzz_yyz, g_zzz_yzz, g_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_xxx[k] = -g_zzz_xxx[k] - g_z_0_zzz_xxx[k] * cd_z[k] + g_z_0_zzz_xxxz[k];

                g_z_0_zzzz_xxy[k] = -g_zzz_xxy[k] - g_z_0_zzz_xxy[k] * cd_z[k] + g_z_0_zzz_xxyz[k];

                g_z_0_zzzz_xxz[k] = -g_zzz_xxz[k] - g_z_0_zzz_xxz[k] * cd_z[k] + g_z_0_zzz_xxzz[k];

                g_z_0_zzzz_xyy[k] = -g_zzz_xyy[k] - g_z_0_zzz_xyy[k] * cd_z[k] + g_z_0_zzz_xyyz[k];

                g_z_0_zzzz_xyz[k] = -g_zzz_xyz[k] - g_z_0_zzz_xyz[k] * cd_z[k] + g_z_0_zzz_xyzz[k];

                g_z_0_zzzz_xzz[k] = -g_zzz_xzz[k] - g_z_0_zzz_xzz[k] * cd_z[k] + g_z_0_zzz_xzzz[k];

                g_z_0_zzzz_yyy[k] = -g_zzz_yyy[k] - g_z_0_zzz_yyy[k] * cd_z[k] + g_z_0_zzz_yyyz[k];

                g_z_0_zzzz_yyz[k] = -g_zzz_yyz[k] - g_z_0_zzz_yyz[k] * cd_z[k] + g_z_0_zzz_yyzz[k];

                g_z_0_zzzz_yzz[k] = -g_zzz_yzz[k] - g_z_0_zzz_yzz[k] * cd_z[k] + g_z_0_zzz_yzzz[k];

                g_z_0_zzzz_zzz[k] = -g_zzz_zzz[k] - g_z_0_zzz_zzz[k] * cd_z[k] + g_z_0_zzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

