#include "ElectronRepulsionGeom1100ContrRecFSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_fsxx,
                                            const size_t idx_geom_01_dsxx,
                                            const size_t idx_geom_10_dsxx,
                                            const size_t idx_geom_11_dsxx,
                                            const size_t idx_geom_11_dpxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_01_off = idx_geom_01_dsxx + i * dcomps + j;

            auto g_0_x_xx_0 = cbuffer.data(ds_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xy_0 = cbuffer.data(ds_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xz_0 = cbuffer.data(ds_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_yy_0 = cbuffer.data(ds_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_yz_0 = cbuffer.data(ds_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_zz_0 = cbuffer.data(ds_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_y_xx_0 = cbuffer.data(ds_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_y_xy_0 = cbuffer.data(ds_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_y_xz_0 = cbuffer.data(ds_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_yy_0 = cbuffer.data(ds_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_yz_0 = cbuffer.data(ds_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_zz_0 = cbuffer.data(ds_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_z_xx_0 = cbuffer.data(ds_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_z_xy_0 = cbuffer.data(ds_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_z_xz_0 = cbuffer.data(ds_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_z_yy_0 = cbuffer.data(ds_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_z_yz_0 = cbuffer.data(ds_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_z_zz_0 = cbuffer.data(ds_geom_01_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_10_off = idx_geom_10_dsxx + i * dcomps + j;

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 5 * ccomps * dcomps);

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * ccomps * dcomps);

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 7 * ccomps * dcomps);

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 11 * ccomps * dcomps);

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * ccomps * dcomps);

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 13 * ccomps * dcomps);

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 14 * ccomps * dcomps);

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 15 * ccomps * dcomps);

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 16 * ccomps * dcomps);

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_11_off = idx_geom_11_dsxx + i * dcomps + j;

            auto g_x_x_xx_0 = cbuffer.data(ds_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xy_0 = cbuffer.data(ds_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xz_0 = cbuffer.data(ds_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_yy_0 = cbuffer.data(ds_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_yz_0 = cbuffer.data(ds_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_zz_0 = cbuffer.data(ds_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_y_xx_0 = cbuffer.data(ds_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_y_xy_0 = cbuffer.data(ds_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_y_xz_0 = cbuffer.data(ds_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_y_yy_0 = cbuffer.data(ds_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_y_yz_0 = cbuffer.data(ds_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_y_zz_0 = cbuffer.data(ds_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_z_xx_0 = cbuffer.data(ds_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_z_xy_0 = cbuffer.data(ds_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_z_xz_0 = cbuffer.data(ds_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_z_yy_0 = cbuffer.data(ds_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_z_yz_0 = cbuffer.data(ds_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_z_zz_0 = cbuffer.data(ds_geom_11_off + 17 * ccomps * dcomps);

            auto g_y_x_xx_0 = cbuffer.data(ds_geom_11_off + 18 * ccomps * dcomps);

            auto g_y_x_xy_0 = cbuffer.data(ds_geom_11_off + 19 * ccomps * dcomps);

            auto g_y_x_xz_0 = cbuffer.data(ds_geom_11_off + 20 * ccomps * dcomps);

            auto g_y_x_yy_0 = cbuffer.data(ds_geom_11_off + 21 * ccomps * dcomps);

            auto g_y_x_yz_0 = cbuffer.data(ds_geom_11_off + 22 * ccomps * dcomps);

            auto g_y_x_zz_0 = cbuffer.data(ds_geom_11_off + 23 * ccomps * dcomps);

            auto g_y_y_xx_0 = cbuffer.data(ds_geom_11_off + 24 * ccomps * dcomps);

            auto g_y_y_xy_0 = cbuffer.data(ds_geom_11_off + 25 * ccomps * dcomps);

            auto g_y_y_xz_0 = cbuffer.data(ds_geom_11_off + 26 * ccomps * dcomps);

            auto g_y_y_yy_0 = cbuffer.data(ds_geom_11_off + 27 * ccomps * dcomps);

            auto g_y_y_yz_0 = cbuffer.data(ds_geom_11_off + 28 * ccomps * dcomps);

            auto g_y_y_zz_0 = cbuffer.data(ds_geom_11_off + 29 * ccomps * dcomps);

            auto g_y_z_xx_0 = cbuffer.data(ds_geom_11_off + 30 * ccomps * dcomps);

            auto g_y_z_xy_0 = cbuffer.data(ds_geom_11_off + 31 * ccomps * dcomps);

            auto g_y_z_xz_0 = cbuffer.data(ds_geom_11_off + 32 * ccomps * dcomps);

            auto g_y_z_yy_0 = cbuffer.data(ds_geom_11_off + 33 * ccomps * dcomps);

            auto g_y_z_yz_0 = cbuffer.data(ds_geom_11_off + 34 * ccomps * dcomps);

            auto g_y_z_zz_0 = cbuffer.data(ds_geom_11_off + 35 * ccomps * dcomps);

            auto g_z_x_xx_0 = cbuffer.data(ds_geom_11_off + 36 * ccomps * dcomps);

            auto g_z_x_xy_0 = cbuffer.data(ds_geom_11_off + 37 * ccomps * dcomps);

            auto g_z_x_xz_0 = cbuffer.data(ds_geom_11_off + 38 * ccomps * dcomps);

            auto g_z_x_yy_0 = cbuffer.data(ds_geom_11_off + 39 * ccomps * dcomps);

            auto g_z_x_yz_0 = cbuffer.data(ds_geom_11_off + 40 * ccomps * dcomps);

            auto g_z_x_zz_0 = cbuffer.data(ds_geom_11_off + 41 * ccomps * dcomps);

            auto g_z_y_xx_0 = cbuffer.data(ds_geom_11_off + 42 * ccomps * dcomps);

            auto g_z_y_xy_0 = cbuffer.data(ds_geom_11_off + 43 * ccomps * dcomps);

            auto g_z_y_xz_0 = cbuffer.data(ds_geom_11_off + 44 * ccomps * dcomps);

            auto g_z_y_yy_0 = cbuffer.data(ds_geom_11_off + 45 * ccomps * dcomps);

            auto g_z_y_yz_0 = cbuffer.data(ds_geom_11_off + 46 * ccomps * dcomps);

            auto g_z_y_zz_0 = cbuffer.data(ds_geom_11_off + 47 * ccomps * dcomps);

            auto g_z_z_xx_0 = cbuffer.data(ds_geom_11_off + 48 * ccomps * dcomps);

            auto g_z_z_xy_0 = cbuffer.data(ds_geom_11_off + 49 * ccomps * dcomps);

            auto g_z_z_xz_0 = cbuffer.data(ds_geom_11_off + 50 * ccomps * dcomps);

            auto g_z_z_yy_0 = cbuffer.data(ds_geom_11_off + 51 * ccomps * dcomps);

            auto g_z_z_yz_0 = cbuffer.data(ds_geom_11_off + 52 * ccomps * dcomps);

            auto g_z_z_zz_0 = cbuffer.data(ds_geom_11_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_11_off = idx_geom_11_dpxx + i * dcomps + j;

            auto g_x_x_xx_x = cbuffer.data(dp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_y = cbuffer.data(dp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_z = cbuffer.data(dp_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xy_x = cbuffer.data(dp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xy_y = cbuffer.data(dp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xy_z = cbuffer.data(dp_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xz_x = cbuffer.data(dp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xz_y = cbuffer.data(dp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xz_z = cbuffer.data(dp_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_yy_x = cbuffer.data(dp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_yy_y = cbuffer.data(dp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_yy_z = cbuffer.data(dp_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_yz_x = cbuffer.data(dp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_yz_y = cbuffer.data(dp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_yz_z = cbuffer.data(dp_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_zz_x = cbuffer.data(dp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_zz_y = cbuffer.data(dp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_zz_z = cbuffer.data(dp_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_xx_x = cbuffer.data(dp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_xx_y = cbuffer.data(dp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_xx_z = cbuffer.data(dp_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_xy_x = cbuffer.data(dp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_xy_y = cbuffer.data(dp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_xy_z = cbuffer.data(dp_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_xz_x = cbuffer.data(dp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_xz_y = cbuffer.data(dp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_xz_z = cbuffer.data(dp_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_yy_x = cbuffer.data(dp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_yy_y = cbuffer.data(dp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_yy_z = cbuffer.data(dp_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_yz_x = cbuffer.data(dp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_yz_y = cbuffer.data(dp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_yz_z = cbuffer.data(dp_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_zz_x = cbuffer.data(dp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_zz_y = cbuffer.data(dp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_zz_z = cbuffer.data(dp_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_xx_x = cbuffer.data(dp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_xx_y = cbuffer.data(dp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_xx_z = cbuffer.data(dp_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_xy_x = cbuffer.data(dp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_xy_y = cbuffer.data(dp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_xy_z = cbuffer.data(dp_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_xz_x = cbuffer.data(dp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_xz_y = cbuffer.data(dp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_xz_z = cbuffer.data(dp_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_z_yy_x = cbuffer.data(dp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_z_yy_y = cbuffer.data(dp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_z_yy_z = cbuffer.data(dp_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_z_yz_x = cbuffer.data(dp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_z_yz_y = cbuffer.data(dp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_z_yz_z = cbuffer.data(dp_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_z_zz_x = cbuffer.data(dp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_z_zz_y = cbuffer.data(dp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_z_zz_z = cbuffer.data(dp_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_xx_x = cbuffer.data(dp_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_xx_y = cbuffer.data(dp_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_xx_z = cbuffer.data(dp_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_xy_x = cbuffer.data(dp_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_xy_y = cbuffer.data(dp_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_xy_z = cbuffer.data(dp_geom_11_off + 59 * ccomps * dcomps);

            auto g_y_x_xz_x = cbuffer.data(dp_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_x_xz_y = cbuffer.data(dp_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_x_xz_z = cbuffer.data(dp_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_x_yy_x = cbuffer.data(dp_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_x_yy_y = cbuffer.data(dp_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_x_yy_z = cbuffer.data(dp_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_x_yz_x = cbuffer.data(dp_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_x_yz_y = cbuffer.data(dp_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_x_yz_z = cbuffer.data(dp_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_x_zz_x = cbuffer.data(dp_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_x_zz_y = cbuffer.data(dp_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_x_zz_z = cbuffer.data(dp_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_xx_x = cbuffer.data(dp_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_xx_y = cbuffer.data(dp_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_xx_z = cbuffer.data(dp_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_y_xy_x = cbuffer.data(dp_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_y_xy_y = cbuffer.data(dp_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_y_xy_z = cbuffer.data(dp_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_y_xz_x = cbuffer.data(dp_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_y_xz_y = cbuffer.data(dp_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_y_xz_z = cbuffer.data(dp_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_y_yy_x = cbuffer.data(dp_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_y_yy_y = cbuffer.data(dp_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_y_yy_z = cbuffer.data(dp_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_y_yz_x = cbuffer.data(dp_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_y_yz_y = cbuffer.data(dp_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_y_yz_z = cbuffer.data(dp_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_y_zz_x = cbuffer.data(dp_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_y_zz_y = cbuffer.data(dp_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_y_zz_z = cbuffer.data(dp_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_z_xx_x = cbuffer.data(dp_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_z_xx_y = cbuffer.data(dp_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_z_xx_z = cbuffer.data(dp_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_z_xy_x = cbuffer.data(dp_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_z_xy_y = cbuffer.data(dp_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_z_xy_z = cbuffer.data(dp_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_z_xz_x = cbuffer.data(dp_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_z_xz_y = cbuffer.data(dp_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_z_xz_z = cbuffer.data(dp_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_z_yy_x = cbuffer.data(dp_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_z_yy_y = cbuffer.data(dp_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_z_yy_z = cbuffer.data(dp_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_z_yz_x = cbuffer.data(dp_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_z_yz_y = cbuffer.data(dp_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_z_yz_z = cbuffer.data(dp_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_z_zz_x = cbuffer.data(dp_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_z_zz_y = cbuffer.data(dp_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_z_zz_z = cbuffer.data(dp_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_x_xx_x = cbuffer.data(dp_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_x_xx_y = cbuffer.data(dp_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_x_xx_z = cbuffer.data(dp_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_x_xy_x = cbuffer.data(dp_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_x_xy_y = cbuffer.data(dp_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_x_xy_z = cbuffer.data(dp_geom_11_off + 113 * ccomps * dcomps);

            auto g_z_x_xz_x = cbuffer.data(dp_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_x_xz_y = cbuffer.data(dp_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_x_xz_z = cbuffer.data(dp_geom_11_off + 116 * ccomps * dcomps);

            auto g_z_x_yy_x = cbuffer.data(dp_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_x_yy_y = cbuffer.data(dp_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_x_yy_z = cbuffer.data(dp_geom_11_off + 119 * ccomps * dcomps);

            auto g_z_x_yz_x = cbuffer.data(dp_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_x_yz_y = cbuffer.data(dp_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_x_yz_z = cbuffer.data(dp_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_x_zz_x = cbuffer.data(dp_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_x_zz_y = cbuffer.data(dp_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_x_zz_z = cbuffer.data(dp_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_y_xx_x = cbuffer.data(dp_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_y_xx_y = cbuffer.data(dp_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_y_xx_z = cbuffer.data(dp_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_y_xy_x = cbuffer.data(dp_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_y_xy_y = cbuffer.data(dp_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_y_xy_z = cbuffer.data(dp_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_y_xz_x = cbuffer.data(dp_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_y_xz_y = cbuffer.data(dp_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_y_xz_z = cbuffer.data(dp_geom_11_off + 134 * ccomps * dcomps);

            auto g_z_y_yy_x = cbuffer.data(dp_geom_11_off + 135 * ccomps * dcomps);

            auto g_z_y_yy_y = cbuffer.data(dp_geom_11_off + 136 * ccomps * dcomps);

            auto g_z_y_yy_z = cbuffer.data(dp_geom_11_off + 137 * ccomps * dcomps);

            auto g_z_y_yz_x = cbuffer.data(dp_geom_11_off + 138 * ccomps * dcomps);

            auto g_z_y_yz_y = cbuffer.data(dp_geom_11_off + 139 * ccomps * dcomps);

            auto g_z_y_yz_z = cbuffer.data(dp_geom_11_off + 140 * ccomps * dcomps);

            auto g_z_y_zz_x = cbuffer.data(dp_geom_11_off + 141 * ccomps * dcomps);

            auto g_z_y_zz_y = cbuffer.data(dp_geom_11_off + 142 * ccomps * dcomps);

            auto g_z_y_zz_z = cbuffer.data(dp_geom_11_off + 143 * ccomps * dcomps);

            auto g_z_z_xx_x = cbuffer.data(dp_geom_11_off + 144 * ccomps * dcomps);

            auto g_z_z_xx_y = cbuffer.data(dp_geom_11_off + 145 * ccomps * dcomps);

            auto g_z_z_xx_z = cbuffer.data(dp_geom_11_off + 146 * ccomps * dcomps);

            auto g_z_z_xy_x = cbuffer.data(dp_geom_11_off + 147 * ccomps * dcomps);

            auto g_z_z_xy_y = cbuffer.data(dp_geom_11_off + 148 * ccomps * dcomps);

            auto g_z_z_xy_z = cbuffer.data(dp_geom_11_off + 149 * ccomps * dcomps);

            auto g_z_z_xz_x = cbuffer.data(dp_geom_11_off + 150 * ccomps * dcomps);

            auto g_z_z_xz_y = cbuffer.data(dp_geom_11_off + 151 * ccomps * dcomps);

            auto g_z_z_xz_z = cbuffer.data(dp_geom_11_off + 152 * ccomps * dcomps);

            auto g_z_z_yy_x = cbuffer.data(dp_geom_11_off + 153 * ccomps * dcomps);

            auto g_z_z_yy_y = cbuffer.data(dp_geom_11_off + 154 * ccomps * dcomps);

            auto g_z_z_yy_z = cbuffer.data(dp_geom_11_off + 155 * ccomps * dcomps);

            auto g_z_z_yz_x = cbuffer.data(dp_geom_11_off + 156 * ccomps * dcomps);

            auto g_z_z_yz_y = cbuffer.data(dp_geom_11_off + 157 * ccomps * dcomps);

            auto g_z_z_yz_z = cbuffer.data(dp_geom_11_off + 158 * ccomps * dcomps);

            auto g_z_z_zz_x = cbuffer.data(dp_geom_11_off + 159 * ccomps * dcomps);

            auto g_z_z_zz_y = cbuffer.data(dp_geom_11_off + 160 * ccomps * dcomps);

            auto g_z_z_zz_z = cbuffer.data(dp_geom_11_off + 161 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fsxx

            const auto fs_geom_11_off = idx_geom_11_fsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxx_0 = cbuffer.data(fs_geom_11_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_0, g_x_0_xx_0, g_x_x_xx_0, g_x_x_xx_x, g_x_x_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxx_0[k] = -g_0_x_xx_0[k] + g_x_0_xx_0[k] - g_x_x_xx_0[k] * ab_x + g_x_x_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxy_0 = cbuffer.data(fs_geom_11_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_0, g_x_x_xx_y, g_x_x_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxy_0[k] = -g_x_x_xx_0[k] * ab_y + g_x_x_xx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxz_0 = cbuffer.data(fs_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_0, g_x_x_xx_z, g_x_x_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxz_0[k] = -g_x_x_xx_0[k] * ab_z + g_x_x_xx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyy_0 = cbuffer.data(fs_geom_11_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xy_0, g_x_x_xy_y, g_x_x_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyy_0[k] = -g_x_x_xy_0[k] * ab_y + g_x_x_xy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyz_0 = cbuffer.data(fs_geom_11_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyz_0, g_x_x_xz_0, g_x_x_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyz_0[k] = -g_x_x_xz_0[k] * ab_y + g_x_x_xz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzz_0 = cbuffer.data(fs_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xz_0, g_x_x_xz_z, g_x_x_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzz_0[k] = -g_x_x_xz_0[k] * ab_z + g_x_x_xz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyy_0 = cbuffer.data(fs_geom_11_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yy_0, g_x_x_yy_y, g_x_x_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyy_0[k] = -g_x_x_yy_0[k] * ab_y + g_x_x_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyz_0 = cbuffer.data(fs_geom_11_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyz_0, g_x_x_yz_0, g_x_x_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyz_0[k] = -g_x_x_yz_0[k] * ab_y + g_x_x_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzz_0 = cbuffer.data(fs_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzz_0, g_x_x_zz_0, g_x_x_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzz_0[k] = -g_x_x_zz_0[k] * ab_y + g_x_x_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzz_0 = cbuffer.data(fs_geom_11_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zz_0, g_x_x_zz_z, g_x_x_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzz_0[k] = -g_x_x_zz_0[k] * ab_z + g_x_x_zz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxx_0 = cbuffer.data(fs_geom_11_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_0, g_x_y_xx_0, g_x_y_xx_x, g_x_y_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxx_0[k] = -g_0_y_xx_0[k] - g_x_y_xx_0[k] * ab_x + g_x_y_xx_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxy_0 = cbuffer.data(fs_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_0, g_x_y_xxy_0, g_x_y_xy_0, g_x_y_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxy_0[k] = -g_0_y_xy_0[k] - g_x_y_xy_0[k] * ab_x + g_x_y_xy_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxz_0 = cbuffer.data(fs_geom_11_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xx_0, g_x_y_xx_z, g_x_y_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxz_0[k] = -g_x_y_xx_0[k] * ab_z + g_x_y_xx_z[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyy_0 = cbuffer.data(fs_geom_11_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_0, g_x_y_xyy_0, g_x_y_yy_0, g_x_y_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyy_0[k] = -g_0_y_yy_0[k] - g_x_y_yy_0[k] * ab_x + g_x_y_yy_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyz_0 = cbuffer.data(fs_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xy_0, g_x_y_xy_z, g_x_y_xyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyz_0[k] = -g_x_y_xy_0[k] * ab_z + g_x_y_xy_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzz_0 = cbuffer.data(fs_geom_11_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xz_0, g_x_y_xz_z, g_x_y_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzz_0[k] = -g_x_y_xz_0[k] * ab_z + g_x_y_xz_z[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyy_0 = cbuffer.data(fs_geom_11_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_0, g_x_y_yy_0, g_x_y_yy_y, g_x_y_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyy_0[k] = g_x_0_yy_0[k] - g_x_y_yy_0[k] * ab_y + g_x_y_yy_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyz_0 = cbuffer.data(fs_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yy_0, g_x_y_yy_z, g_x_y_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyz_0[k] = -g_x_y_yy_0[k] * ab_z + g_x_y_yy_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzz_0 = cbuffer.data(fs_geom_11_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yz_0, g_x_y_yz_z, g_x_y_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzz_0[k] = -g_x_y_yz_0[k] * ab_z + g_x_y_yz_z[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzz_0 = cbuffer.data(fs_geom_11_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zz_0, g_x_y_zz_z, g_x_y_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzz_0[k] = -g_x_y_zz_0[k] * ab_z + g_x_y_zz_z[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxx_0 = cbuffer.data(fs_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_0, g_x_z_xx_0, g_x_z_xx_x, g_x_z_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxx_0[k] = -g_0_z_xx_0[k] - g_x_z_xx_0[k] * ab_x + g_x_z_xx_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxy_0 = cbuffer.data(fs_geom_11_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xx_0, g_x_z_xx_y, g_x_z_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxy_0[k] = -g_x_z_xx_0[k] * ab_y + g_x_z_xx_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxz_0 = cbuffer.data(fs_geom_11_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_0, g_x_z_xxz_0, g_x_z_xz_0, g_x_z_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxz_0[k] = -g_0_z_xz_0[k] - g_x_z_xz_0[k] * ab_x + g_x_z_xz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyy_0 = cbuffer.data(fs_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xy_0, g_x_z_xy_y, g_x_z_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyy_0[k] = -g_x_z_xy_0[k] * ab_y + g_x_z_xy_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyz_0 = cbuffer.data(fs_geom_11_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyz_0, g_x_z_xz_0, g_x_z_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyz_0[k] = -g_x_z_xz_0[k] * ab_y + g_x_z_xz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzz_0 = cbuffer.data(fs_geom_11_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_0, g_x_z_xzz_0, g_x_z_zz_0, g_x_z_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzz_0[k] = -g_0_z_zz_0[k] - g_x_z_zz_0[k] * ab_x + g_x_z_zz_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyy_0 = cbuffer.data(fs_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yy_0, g_x_z_yy_y, g_x_z_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyy_0[k] = -g_x_z_yy_0[k] * ab_y + g_x_z_yy_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyz_0 = cbuffer.data(fs_geom_11_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyz_0, g_x_z_yz_0, g_x_z_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyz_0[k] = -g_x_z_yz_0[k] * ab_y + g_x_z_yz_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzz_0 = cbuffer.data(fs_geom_11_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzz_0, g_x_z_zz_0, g_x_z_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzz_0[k] = -g_x_z_zz_0[k] * ab_y + g_x_z_zz_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzz_0 = cbuffer.data(fs_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_0, g_x_z_zz_0, g_x_z_zz_z, g_x_z_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzz_0[k] = g_x_0_zz_0[k] - g_x_z_zz_0[k] * ab_z + g_x_z_zz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxx_0 = cbuffer.data(fs_geom_11_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_0, g_y_x_xx_0, g_y_x_xx_x, g_y_x_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxx_0[k] = g_y_0_xx_0[k] - g_y_x_xx_0[k] * ab_x + g_y_x_xx_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxy_0 = cbuffer.data(fs_geom_11_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_0, g_y_x_xxy_0, g_y_x_xy_0, g_y_x_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxy_0[k] = g_y_0_xy_0[k] - g_y_x_xy_0[k] * ab_x + g_y_x_xy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxz_0 = cbuffer.data(fs_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xx_0, g_y_x_xx_z, g_y_x_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxz_0[k] = -g_y_x_xx_0[k] * ab_z + g_y_x_xx_z[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyy_0 = cbuffer.data(fs_geom_11_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_0, g_y_x_xyy_0, g_y_x_yy_0, g_y_x_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyy_0[k] = g_y_0_yy_0[k] - g_y_x_yy_0[k] * ab_x + g_y_x_yy_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyz_0 = cbuffer.data(fs_geom_11_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xy_0, g_y_x_xy_z, g_y_x_xyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyz_0[k] = -g_y_x_xy_0[k] * ab_z + g_y_x_xy_z[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzz_0 = cbuffer.data(fs_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xz_0, g_y_x_xz_z, g_y_x_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzz_0[k] = -g_y_x_xz_0[k] * ab_z + g_y_x_xz_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyy_0 = cbuffer.data(fs_geom_11_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_0, g_y_x_yy_0, g_y_x_yy_y, g_y_x_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyy_0[k] = -g_0_x_yy_0[k] - g_y_x_yy_0[k] * ab_y + g_y_x_yy_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyz_0 = cbuffer.data(fs_geom_11_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yy_0, g_y_x_yy_z, g_y_x_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyz_0[k] = -g_y_x_yy_0[k] * ab_z + g_y_x_yy_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzz_0 = cbuffer.data(fs_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yz_0, g_y_x_yz_z, g_y_x_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzz_0[k] = -g_y_x_yz_0[k] * ab_z + g_y_x_yz_z[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzz_0 = cbuffer.data(fs_geom_11_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zz_0, g_y_x_zz_z, g_y_x_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzz_0[k] = -g_y_x_zz_0[k] * ab_z + g_y_x_zz_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxx_0 = cbuffer.data(fs_geom_11_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xx_0, g_y_y_xx_x, g_y_y_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxx_0[k] = -g_y_y_xx_0[k] * ab_x + g_y_y_xx_x[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxy_0 = cbuffer.data(fs_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxy_0, g_y_y_xy_0, g_y_y_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxy_0[k] = -g_y_y_xy_0[k] * ab_x + g_y_y_xy_x[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxz_0 = cbuffer.data(fs_geom_11_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxz_0, g_y_y_xz_0, g_y_y_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxz_0[k] = -g_y_y_xz_0[k] * ab_x + g_y_y_xz_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyy_0 = cbuffer.data(fs_geom_11_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyy_0, g_y_y_yy_0, g_y_y_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyy_0[k] = -g_y_y_yy_0[k] * ab_x + g_y_y_yy_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyz_0 = cbuffer.data(fs_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyz_0, g_y_y_yz_0, g_y_y_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyz_0[k] = -g_y_y_yz_0[k] * ab_x + g_y_y_yz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzz_0 = cbuffer.data(fs_geom_11_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzz_0, g_y_y_zz_0, g_y_y_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzz_0[k] = -g_y_y_zz_0[k] * ab_x + g_y_y_zz_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyy_0 = cbuffer.data(fs_geom_11_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_0, g_y_0_yy_0, g_y_y_yy_0, g_y_y_yy_y, g_y_y_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyy_0[k] = -g_0_y_yy_0[k] + g_y_0_yy_0[k] - g_y_y_yy_0[k] * ab_y + g_y_y_yy_y[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyz_0 = cbuffer.data(fs_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yy_0, g_y_y_yy_z, g_y_y_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyz_0[k] = -g_y_y_yy_0[k] * ab_z + g_y_y_yy_z[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzz_0 = cbuffer.data(fs_geom_11_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yz_0, g_y_y_yz_z, g_y_y_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzz_0[k] = -g_y_y_yz_0[k] * ab_z + g_y_y_yz_z[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzz_0 = cbuffer.data(fs_geom_11_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zz_0, g_y_y_zz_z, g_y_y_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzz_0[k] = -g_y_y_zz_0[k] * ab_z + g_y_y_zz_z[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxx_0 = cbuffer.data(fs_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xx_0, g_y_z_xx_x, g_y_z_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxx_0[k] = -g_y_z_xx_0[k] * ab_x + g_y_z_xx_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxy_0 = cbuffer.data(fs_geom_11_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxy_0, g_y_z_xy_0, g_y_z_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxy_0[k] = -g_y_z_xy_0[k] * ab_x + g_y_z_xy_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxz_0 = cbuffer.data(fs_geom_11_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxz_0, g_y_z_xz_0, g_y_z_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxz_0[k] = -g_y_z_xz_0[k] * ab_x + g_y_z_xz_x[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyy_0 = cbuffer.data(fs_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyy_0, g_y_z_yy_0, g_y_z_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyy_0[k] = -g_y_z_yy_0[k] * ab_x + g_y_z_yy_x[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyz_0 = cbuffer.data(fs_geom_11_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyz_0, g_y_z_yz_0, g_y_z_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyz_0[k] = -g_y_z_yz_0[k] * ab_x + g_y_z_yz_x[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzz_0 = cbuffer.data(fs_geom_11_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzz_0, g_y_z_zz_0, g_y_z_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzz_0[k] = -g_y_z_zz_0[k] * ab_x + g_y_z_zz_x[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyy_0 = cbuffer.data(fs_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_0, g_y_z_yy_0, g_y_z_yy_y, g_y_z_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyy_0[k] = -g_0_z_yy_0[k] - g_y_z_yy_0[k] * ab_y + g_y_z_yy_y[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyz_0 = cbuffer.data(fs_geom_11_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_0, g_y_z_yyz_0, g_y_z_yz_0, g_y_z_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyz_0[k] = -g_0_z_yz_0[k] - g_y_z_yz_0[k] * ab_y + g_y_z_yz_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzz_0 = cbuffer.data(fs_geom_11_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_0, g_y_z_yzz_0, g_y_z_zz_0, g_y_z_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzz_0[k] = -g_0_z_zz_0[k] - g_y_z_zz_0[k] * ab_y + g_y_z_zz_y[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzz_0 = cbuffer.data(fs_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_0, g_y_z_zz_0, g_y_z_zz_z, g_y_z_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzz_0[k] = g_y_0_zz_0[k] - g_y_z_zz_0[k] * ab_z + g_y_z_zz_z[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxx_0 = cbuffer.data(fs_geom_11_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_0, g_z_x_xx_0, g_z_x_xx_x, g_z_x_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxx_0[k] = g_z_0_xx_0[k] - g_z_x_xx_0[k] * ab_x + g_z_x_xx_x[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxy_0 = cbuffer.data(fs_geom_11_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xx_0, g_z_x_xx_y, g_z_x_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxy_0[k] = -g_z_x_xx_0[k] * ab_y + g_z_x_xx_y[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxz_0 = cbuffer.data(fs_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_0, g_z_x_xxz_0, g_z_x_xz_0, g_z_x_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxz_0[k] = g_z_0_xz_0[k] - g_z_x_xz_0[k] * ab_x + g_z_x_xz_x[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyy_0 = cbuffer.data(fs_geom_11_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xy_0, g_z_x_xy_y, g_z_x_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyy_0[k] = -g_z_x_xy_0[k] * ab_y + g_z_x_xy_y[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyz_0 = cbuffer.data(fs_geom_11_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyz_0, g_z_x_xz_0, g_z_x_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyz_0[k] = -g_z_x_xz_0[k] * ab_y + g_z_x_xz_y[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzz_0 = cbuffer.data(fs_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_0, g_z_x_xzz_0, g_z_x_zz_0, g_z_x_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzz_0[k] = g_z_0_zz_0[k] - g_z_x_zz_0[k] * ab_x + g_z_x_zz_x[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyy_0 = cbuffer.data(fs_geom_11_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yy_0, g_z_x_yy_y, g_z_x_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyy_0[k] = -g_z_x_yy_0[k] * ab_y + g_z_x_yy_y[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyz_0 = cbuffer.data(fs_geom_11_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyz_0, g_z_x_yz_0, g_z_x_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyz_0[k] = -g_z_x_yz_0[k] * ab_y + g_z_x_yz_y[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzz_0 = cbuffer.data(fs_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzz_0, g_z_x_zz_0, g_z_x_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzz_0[k] = -g_z_x_zz_0[k] * ab_y + g_z_x_zz_y[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzz_0 = cbuffer.data(fs_geom_11_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_0, g_z_x_zz_0, g_z_x_zz_z, g_z_x_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzz_0[k] = -g_0_x_zz_0[k] - g_z_x_zz_0[k] * ab_z + g_z_x_zz_z[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxx_0 = cbuffer.data(fs_geom_11_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xx_0, g_z_y_xx_x, g_z_y_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxx_0[k] = -g_z_y_xx_0[k] * ab_x + g_z_y_xx_x[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxy_0 = cbuffer.data(fs_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxy_0, g_z_y_xy_0, g_z_y_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxy_0[k] = -g_z_y_xy_0[k] * ab_x + g_z_y_xy_x[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxz_0 = cbuffer.data(fs_geom_11_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxz_0, g_z_y_xz_0, g_z_y_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxz_0[k] = -g_z_y_xz_0[k] * ab_x + g_z_y_xz_x[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyy_0 = cbuffer.data(fs_geom_11_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyy_0, g_z_y_yy_0, g_z_y_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyy_0[k] = -g_z_y_yy_0[k] * ab_x + g_z_y_yy_x[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyz_0 = cbuffer.data(fs_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyz_0, g_z_y_yz_0, g_z_y_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyz_0[k] = -g_z_y_yz_0[k] * ab_x + g_z_y_yz_x[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzz_0 = cbuffer.data(fs_geom_11_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzz_0, g_z_y_zz_0, g_z_y_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzz_0[k] = -g_z_y_zz_0[k] * ab_x + g_z_y_zz_x[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyy_0 = cbuffer.data(fs_geom_11_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_0, g_z_y_yy_0, g_z_y_yy_y, g_z_y_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyy_0[k] = g_z_0_yy_0[k] - g_z_y_yy_0[k] * ab_y + g_z_y_yy_y[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyz_0 = cbuffer.data(fs_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_0, g_z_y_yyz_0, g_z_y_yz_0, g_z_y_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyz_0[k] = g_z_0_yz_0[k] - g_z_y_yz_0[k] * ab_y + g_z_y_yz_y[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzz_0 = cbuffer.data(fs_geom_11_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_0, g_z_y_yzz_0, g_z_y_zz_0, g_z_y_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzz_0[k] = g_z_0_zz_0[k] - g_z_y_zz_0[k] * ab_y + g_z_y_zz_y[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzz_0 = cbuffer.data(fs_geom_11_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_0, g_z_y_zz_0, g_z_y_zz_z, g_z_y_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzz_0[k] = -g_0_y_zz_0[k] - g_z_y_zz_0[k] * ab_z + g_z_y_zz_z[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxx_0 = cbuffer.data(fs_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xx_0, g_z_z_xx_x, g_z_z_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxx_0[k] = -g_z_z_xx_0[k] * ab_x + g_z_z_xx_x[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxy_0 = cbuffer.data(fs_geom_11_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxy_0, g_z_z_xy_0, g_z_z_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxy_0[k] = -g_z_z_xy_0[k] * ab_x + g_z_z_xy_x[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxz_0 = cbuffer.data(fs_geom_11_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxz_0, g_z_z_xz_0, g_z_z_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxz_0[k] = -g_z_z_xz_0[k] * ab_x + g_z_z_xz_x[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyy_0 = cbuffer.data(fs_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyy_0, g_z_z_yy_0, g_z_z_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyy_0[k] = -g_z_z_yy_0[k] * ab_x + g_z_z_yy_x[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyz_0 = cbuffer.data(fs_geom_11_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyz_0, g_z_z_yz_0, g_z_z_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyz_0[k] = -g_z_z_yz_0[k] * ab_x + g_z_z_yz_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzz_0 = cbuffer.data(fs_geom_11_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzz_0, g_z_z_zz_0, g_z_z_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzz_0[k] = -g_z_z_zz_0[k] * ab_x + g_z_z_zz_x[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyy_0 = cbuffer.data(fs_geom_11_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yy_0, g_z_z_yy_y, g_z_z_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyy_0[k] = -g_z_z_yy_0[k] * ab_y + g_z_z_yy_y[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyz_0 = cbuffer.data(fs_geom_11_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyz_0, g_z_z_yz_0, g_z_z_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyz_0[k] = -g_z_z_yz_0[k] * ab_y + g_z_z_yz_y[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzz_0 = cbuffer.data(fs_geom_11_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzz_0, g_z_z_zz_0, g_z_z_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzz_0[k] = -g_z_z_zz_0[k] * ab_y + g_z_z_zz_y[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzz_0 = cbuffer.data(fs_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_0, g_z_0_zz_0, g_z_z_zz_0, g_z_z_zz_z, g_z_z_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzz_0[k] = -g_0_z_zz_0[k] + g_z_0_zz_0[k] - g_z_z_zz_0[k] * ab_z + g_z_z_zz_z[k];
            }
        }
    }
}

} // erirec namespace

