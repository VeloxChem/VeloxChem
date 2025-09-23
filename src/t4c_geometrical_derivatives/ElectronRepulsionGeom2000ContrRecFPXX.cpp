#include "ElectronRepulsionGeom2000ContrRecFPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_fpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_fpxx,
                                            const size_t idx_geom_10_dpxx,
                                            const size_t idx_geom_20_dpxx,
                                            const size_t idx_geom_20_ddxx,
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
            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_10_off = idx_geom_10_dpxx + i * dcomps + j;

            auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xy_x = cbuffer.data(dp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xy_z = cbuffer.data(dp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xz_x = cbuffer.data(dp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_yy_x = cbuffer.data(dp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_yy_z = cbuffer.data(dp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_yz_x = cbuffer.data(dp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_yz_z = cbuffer.data(dp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_zz_x = cbuffer.data(dp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_xx_y = cbuffer.data(dp_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_xx_z = cbuffer.data(dp_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_xy_y = cbuffer.data(dp_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_xy_z = cbuffer.data(dp_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_xz_y = cbuffer.data(dp_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_xz_z = cbuffer.data(dp_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_yz_y = cbuffer.data(dp_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_zz_y = cbuffer.data(dp_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_xx_y = cbuffer.data(dp_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_xx_z = cbuffer.data(dp_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_xy_y = cbuffer.data(dp_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_xy_z = cbuffer.data(dp_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_xz_y = cbuffer.data(dp_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_xz_z = cbuffer.data(dp_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_yy_z = cbuffer.data(dp_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_yz_z = cbuffer.data(dp_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_20_off = idx_geom_20_dpxx + i * dcomps + j;

            auto g_xx_0_xx_x = cbuffer.data(dp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_y = cbuffer.data(dp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_z = cbuffer.data(dp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xy_x = cbuffer.data(dp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xy_y = cbuffer.data(dp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xy_z = cbuffer.data(dp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xz_x = cbuffer.data(dp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xz_y = cbuffer.data(dp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xz_z = cbuffer.data(dp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_yy_x = cbuffer.data(dp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_yy_y = cbuffer.data(dp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_yy_z = cbuffer.data(dp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_yz_x = cbuffer.data(dp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_yz_y = cbuffer.data(dp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_yz_z = cbuffer.data(dp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_zz_x = cbuffer.data(dp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_zz_y = cbuffer.data(dp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_zz_z = cbuffer.data(dp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_xx_x = cbuffer.data(dp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_xx_y = cbuffer.data(dp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xy_0_xx_z = cbuffer.data(dp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_xy_x = cbuffer.data(dp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_xy_y = cbuffer.data(dp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_xy_z = cbuffer.data(dp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_xz_x = cbuffer.data(dp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_xz_y = cbuffer.data(dp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_xz_z = cbuffer.data(dp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_yy_x = cbuffer.data(dp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_yy_y = cbuffer.data(dp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_yy_z = cbuffer.data(dp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_yz_x = cbuffer.data(dp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_yz_y = cbuffer.data(dp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_yz_z = cbuffer.data(dp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_zz_x = cbuffer.data(dp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_zz_y = cbuffer.data(dp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_zz_z = cbuffer.data(dp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xz_0_xx_x = cbuffer.data(dp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xz_0_xx_y = cbuffer.data(dp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xz_0_xx_z = cbuffer.data(dp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xz_0_xy_x = cbuffer.data(dp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xz_0_xy_y = cbuffer.data(dp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xz_0_xy_z = cbuffer.data(dp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_xz_x = cbuffer.data(dp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_xz_y = cbuffer.data(dp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_xz_z = cbuffer.data(dp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xz_0_yy_x = cbuffer.data(dp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xz_0_yy_y = cbuffer.data(dp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xz_0_yy_z = cbuffer.data(dp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xz_0_yz_x = cbuffer.data(dp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xz_0_yz_y = cbuffer.data(dp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xz_0_yz_z = cbuffer.data(dp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xz_0_zz_x = cbuffer.data(dp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xz_0_zz_y = cbuffer.data(dp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xz_0_zz_z = cbuffer.data(dp_geom_20_off + 53 * ccomps * dcomps);

            auto g_yy_0_xx_x = cbuffer.data(dp_geom_20_off + 54 * ccomps * dcomps);

            auto g_yy_0_xx_y = cbuffer.data(dp_geom_20_off + 55 * ccomps * dcomps);

            auto g_yy_0_xx_z = cbuffer.data(dp_geom_20_off + 56 * ccomps * dcomps);

            auto g_yy_0_xy_x = cbuffer.data(dp_geom_20_off + 57 * ccomps * dcomps);

            auto g_yy_0_xy_y = cbuffer.data(dp_geom_20_off + 58 * ccomps * dcomps);

            auto g_yy_0_xy_z = cbuffer.data(dp_geom_20_off + 59 * ccomps * dcomps);

            auto g_yy_0_xz_x = cbuffer.data(dp_geom_20_off + 60 * ccomps * dcomps);

            auto g_yy_0_xz_y = cbuffer.data(dp_geom_20_off + 61 * ccomps * dcomps);

            auto g_yy_0_xz_z = cbuffer.data(dp_geom_20_off + 62 * ccomps * dcomps);

            auto g_yy_0_yy_x = cbuffer.data(dp_geom_20_off + 63 * ccomps * dcomps);

            auto g_yy_0_yy_y = cbuffer.data(dp_geom_20_off + 64 * ccomps * dcomps);

            auto g_yy_0_yy_z = cbuffer.data(dp_geom_20_off + 65 * ccomps * dcomps);

            auto g_yy_0_yz_x = cbuffer.data(dp_geom_20_off + 66 * ccomps * dcomps);

            auto g_yy_0_yz_y = cbuffer.data(dp_geom_20_off + 67 * ccomps * dcomps);

            auto g_yy_0_yz_z = cbuffer.data(dp_geom_20_off + 68 * ccomps * dcomps);

            auto g_yy_0_zz_x = cbuffer.data(dp_geom_20_off + 69 * ccomps * dcomps);

            auto g_yy_0_zz_y = cbuffer.data(dp_geom_20_off + 70 * ccomps * dcomps);

            auto g_yy_0_zz_z = cbuffer.data(dp_geom_20_off + 71 * ccomps * dcomps);

            auto g_yz_0_xx_x = cbuffer.data(dp_geom_20_off + 72 * ccomps * dcomps);

            auto g_yz_0_xx_y = cbuffer.data(dp_geom_20_off + 73 * ccomps * dcomps);

            auto g_yz_0_xx_z = cbuffer.data(dp_geom_20_off + 74 * ccomps * dcomps);

            auto g_yz_0_xy_x = cbuffer.data(dp_geom_20_off + 75 * ccomps * dcomps);

            auto g_yz_0_xy_y = cbuffer.data(dp_geom_20_off + 76 * ccomps * dcomps);

            auto g_yz_0_xy_z = cbuffer.data(dp_geom_20_off + 77 * ccomps * dcomps);

            auto g_yz_0_xz_x = cbuffer.data(dp_geom_20_off + 78 * ccomps * dcomps);

            auto g_yz_0_xz_y = cbuffer.data(dp_geom_20_off + 79 * ccomps * dcomps);

            auto g_yz_0_xz_z = cbuffer.data(dp_geom_20_off + 80 * ccomps * dcomps);

            auto g_yz_0_yy_x = cbuffer.data(dp_geom_20_off + 81 * ccomps * dcomps);

            auto g_yz_0_yy_y = cbuffer.data(dp_geom_20_off + 82 * ccomps * dcomps);

            auto g_yz_0_yy_z = cbuffer.data(dp_geom_20_off + 83 * ccomps * dcomps);

            auto g_yz_0_yz_x = cbuffer.data(dp_geom_20_off + 84 * ccomps * dcomps);

            auto g_yz_0_yz_y = cbuffer.data(dp_geom_20_off + 85 * ccomps * dcomps);

            auto g_yz_0_yz_z = cbuffer.data(dp_geom_20_off + 86 * ccomps * dcomps);

            auto g_yz_0_zz_x = cbuffer.data(dp_geom_20_off + 87 * ccomps * dcomps);

            auto g_yz_0_zz_y = cbuffer.data(dp_geom_20_off + 88 * ccomps * dcomps);

            auto g_yz_0_zz_z = cbuffer.data(dp_geom_20_off + 89 * ccomps * dcomps);

            auto g_zz_0_xx_x = cbuffer.data(dp_geom_20_off + 90 * ccomps * dcomps);

            auto g_zz_0_xx_y = cbuffer.data(dp_geom_20_off + 91 * ccomps * dcomps);

            auto g_zz_0_xx_z = cbuffer.data(dp_geom_20_off + 92 * ccomps * dcomps);

            auto g_zz_0_xy_x = cbuffer.data(dp_geom_20_off + 93 * ccomps * dcomps);

            auto g_zz_0_xy_y = cbuffer.data(dp_geom_20_off + 94 * ccomps * dcomps);

            auto g_zz_0_xy_z = cbuffer.data(dp_geom_20_off + 95 * ccomps * dcomps);

            auto g_zz_0_xz_x = cbuffer.data(dp_geom_20_off + 96 * ccomps * dcomps);

            auto g_zz_0_xz_y = cbuffer.data(dp_geom_20_off + 97 * ccomps * dcomps);

            auto g_zz_0_xz_z = cbuffer.data(dp_geom_20_off + 98 * ccomps * dcomps);

            auto g_zz_0_yy_x = cbuffer.data(dp_geom_20_off + 99 * ccomps * dcomps);

            auto g_zz_0_yy_y = cbuffer.data(dp_geom_20_off + 100 * ccomps * dcomps);

            auto g_zz_0_yy_z = cbuffer.data(dp_geom_20_off + 101 * ccomps * dcomps);

            auto g_zz_0_yz_x = cbuffer.data(dp_geom_20_off + 102 * ccomps * dcomps);

            auto g_zz_0_yz_y = cbuffer.data(dp_geom_20_off + 103 * ccomps * dcomps);

            auto g_zz_0_yz_z = cbuffer.data(dp_geom_20_off + 104 * ccomps * dcomps);

            auto g_zz_0_zz_x = cbuffer.data(dp_geom_20_off + 105 * ccomps * dcomps);

            auto g_zz_0_zz_y = cbuffer.data(dp_geom_20_off + 106 * ccomps * dcomps);

            auto g_zz_0_zz_z = cbuffer.data(dp_geom_20_off + 107 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DDSS

            const auto dd_geom_20_off = idx_geom_20_ddxx + i * dcomps + j;

            auto g_xx_0_xx_xx = cbuffer.data(dd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xy = cbuffer.data(dd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xz = cbuffer.data(dd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_yy = cbuffer.data(dd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_yz = cbuffer.data(dd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_zz = cbuffer.data(dd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xy_xx = cbuffer.data(dd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xy_xy = cbuffer.data(dd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xy_xz = cbuffer.data(dd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xy_yy = cbuffer.data(dd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xy_yz = cbuffer.data(dd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xy_zz = cbuffer.data(dd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xz_xx = cbuffer.data(dd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xz_xy = cbuffer.data(dd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xz_xz = cbuffer.data(dd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xz_yy = cbuffer.data(dd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xz_yz = cbuffer.data(dd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xz_zz = cbuffer.data(dd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_yy_xx = cbuffer.data(dd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_yy_xy = cbuffer.data(dd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_yy_xz = cbuffer.data(dd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_yy_yy = cbuffer.data(dd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_yy_yz = cbuffer.data(dd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_yy_zz = cbuffer.data(dd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_yz_xx = cbuffer.data(dd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_yz_xy = cbuffer.data(dd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_yz_xz = cbuffer.data(dd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_yz_yy = cbuffer.data(dd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_yz_yz = cbuffer.data(dd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_yz_zz = cbuffer.data(dd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_zz_xx = cbuffer.data(dd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_zz_xy = cbuffer.data(dd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_zz_xz = cbuffer.data(dd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_zz_yy = cbuffer.data(dd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_zz_yz = cbuffer.data(dd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_zz_zz = cbuffer.data(dd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_xx_xx = cbuffer.data(dd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_xx_xy = cbuffer.data(dd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_xx_xz = cbuffer.data(dd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_xx_yy = cbuffer.data(dd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_xx_yz = cbuffer.data(dd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_xx_zz = cbuffer.data(dd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_xy_xx = cbuffer.data(dd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_xy_xy = cbuffer.data(dd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_xy_xz = cbuffer.data(dd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_xy_yy = cbuffer.data(dd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_xy_yz = cbuffer.data(dd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_xy_zz = cbuffer.data(dd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_xz_xx = cbuffer.data(dd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_xz_xy = cbuffer.data(dd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_xz_xz = cbuffer.data(dd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_xz_yy = cbuffer.data(dd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_xz_yz = cbuffer.data(dd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_xz_zz = cbuffer.data(dd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_yy_xx = cbuffer.data(dd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_yy_xy = cbuffer.data(dd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_yy_xz = cbuffer.data(dd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_yy_yy = cbuffer.data(dd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_yy_yz = cbuffer.data(dd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_yy_zz = cbuffer.data(dd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_yz_xx = cbuffer.data(dd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_yz_xy = cbuffer.data(dd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_yz_xz = cbuffer.data(dd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_yz_yy = cbuffer.data(dd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_yz_yz = cbuffer.data(dd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_yz_zz = cbuffer.data(dd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_zz_xx = cbuffer.data(dd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_zz_xy = cbuffer.data(dd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_zz_xz = cbuffer.data(dd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_zz_yy = cbuffer.data(dd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_zz_yz = cbuffer.data(dd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_zz_zz = cbuffer.data(dd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_xx_xx = cbuffer.data(dd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_xx_xy = cbuffer.data(dd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_xx_xz = cbuffer.data(dd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_xx_yy = cbuffer.data(dd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_xx_yz = cbuffer.data(dd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_xx_zz = cbuffer.data(dd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_xy_xx = cbuffer.data(dd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_xy_xy = cbuffer.data(dd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_xy_xz = cbuffer.data(dd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_xy_yy = cbuffer.data(dd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_xy_yz = cbuffer.data(dd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_xy_zz = cbuffer.data(dd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xz_0_xz_xx = cbuffer.data(dd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_xz_xy = cbuffer.data(dd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_xz_xz = cbuffer.data(dd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xz_0_xz_yy = cbuffer.data(dd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_xz_yz = cbuffer.data(dd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_xz_zz = cbuffer.data(dd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xz_0_yy_xx = cbuffer.data(dd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_yy_xy = cbuffer.data(dd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_yy_xz = cbuffer.data(dd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xz_0_yy_yy = cbuffer.data(dd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_yy_yz = cbuffer.data(dd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_yy_zz = cbuffer.data(dd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xz_0_yz_xx = cbuffer.data(dd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_yz_xy = cbuffer.data(dd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_yz_xz = cbuffer.data(dd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xz_0_yz_yy = cbuffer.data(dd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_yz_yz = cbuffer.data(dd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_yz_zz = cbuffer.data(dd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xz_0_zz_xx = cbuffer.data(dd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_zz_xy = cbuffer.data(dd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_zz_xz = cbuffer.data(dd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xz_0_zz_yy = cbuffer.data(dd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_zz_yz = cbuffer.data(dd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_zz_zz = cbuffer.data(dd_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_xx_xx = cbuffer.data(dd_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_xx_xy = cbuffer.data(dd_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_xx_xz = cbuffer.data(dd_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_xx_yy = cbuffer.data(dd_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_xx_yz = cbuffer.data(dd_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_xx_zz = cbuffer.data(dd_geom_20_off + 113 * ccomps * dcomps);

            auto g_yy_0_xy_xx = cbuffer.data(dd_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_xy_xy = cbuffer.data(dd_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_xy_xz = cbuffer.data(dd_geom_20_off + 116 * ccomps * dcomps);

            auto g_yy_0_xy_yy = cbuffer.data(dd_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_xy_yz = cbuffer.data(dd_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_xy_zz = cbuffer.data(dd_geom_20_off + 119 * ccomps * dcomps);

            auto g_yy_0_xz_xx = cbuffer.data(dd_geom_20_off + 120 * ccomps * dcomps);

            auto g_yy_0_xz_xy = cbuffer.data(dd_geom_20_off + 121 * ccomps * dcomps);

            auto g_yy_0_xz_xz = cbuffer.data(dd_geom_20_off + 122 * ccomps * dcomps);

            auto g_yy_0_xz_yy = cbuffer.data(dd_geom_20_off + 123 * ccomps * dcomps);

            auto g_yy_0_xz_yz = cbuffer.data(dd_geom_20_off + 124 * ccomps * dcomps);

            auto g_yy_0_xz_zz = cbuffer.data(dd_geom_20_off + 125 * ccomps * dcomps);

            auto g_yy_0_yy_xx = cbuffer.data(dd_geom_20_off + 126 * ccomps * dcomps);

            auto g_yy_0_yy_xy = cbuffer.data(dd_geom_20_off + 127 * ccomps * dcomps);

            auto g_yy_0_yy_xz = cbuffer.data(dd_geom_20_off + 128 * ccomps * dcomps);

            auto g_yy_0_yy_yy = cbuffer.data(dd_geom_20_off + 129 * ccomps * dcomps);

            auto g_yy_0_yy_yz = cbuffer.data(dd_geom_20_off + 130 * ccomps * dcomps);

            auto g_yy_0_yy_zz = cbuffer.data(dd_geom_20_off + 131 * ccomps * dcomps);

            auto g_yy_0_yz_xx = cbuffer.data(dd_geom_20_off + 132 * ccomps * dcomps);

            auto g_yy_0_yz_xy = cbuffer.data(dd_geom_20_off + 133 * ccomps * dcomps);

            auto g_yy_0_yz_xz = cbuffer.data(dd_geom_20_off + 134 * ccomps * dcomps);

            auto g_yy_0_yz_yy = cbuffer.data(dd_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_yz_yz = cbuffer.data(dd_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_yz_zz = cbuffer.data(dd_geom_20_off + 137 * ccomps * dcomps);

            auto g_yy_0_zz_xx = cbuffer.data(dd_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_zz_xy = cbuffer.data(dd_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_zz_xz = cbuffer.data(dd_geom_20_off + 140 * ccomps * dcomps);

            auto g_yy_0_zz_yy = cbuffer.data(dd_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_zz_yz = cbuffer.data(dd_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_zz_zz = cbuffer.data(dd_geom_20_off + 143 * ccomps * dcomps);

            auto g_yz_0_xx_xx = cbuffer.data(dd_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_xx_xy = cbuffer.data(dd_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_xx_xz = cbuffer.data(dd_geom_20_off + 146 * ccomps * dcomps);

            auto g_yz_0_xx_yy = cbuffer.data(dd_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_xx_yz = cbuffer.data(dd_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_xx_zz = cbuffer.data(dd_geom_20_off + 149 * ccomps * dcomps);

            auto g_yz_0_xy_xx = cbuffer.data(dd_geom_20_off + 150 * ccomps * dcomps);

            auto g_yz_0_xy_xy = cbuffer.data(dd_geom_20_off + 151 * ccomps * dcomps);

            auto g_yz_0_xy_xz = cbuffer.data(dd_geom_20_off + 152 * ccomps * dcomps);

            auto g_yz_0_xy_yy = cbuffer.data(dd_geom_20_off + 153 * ccomps * dcomps);

            auto g_yz_0_xy_yz = cbuffer.data(dd_geom_20_off + 154 * ccomps * dcomps);

            auto g_yz_0_xy_zz = cbuffer.data(dd_geom_20_off + 155 * ccomps * dcomps);

            auto g_yz_0_xz_xx = cbuffer.data(dd_geom_20_off + 156 * ccomps * dcomps);

            auto g_yz_0_xz_xy = cbuffer.data(dd_geom_20_off + 157 * ccomps * dcomps);

            auto g_yz_0_xz_xz = cbuffer.data(dd_geom_20_off + 158 * ccomps * dcomps);

            auto g_yz_0_xz_yy = cbuffer.data(dd_geom_20_off + 159 * ccomps * dcomps);

            auto g_yz_0_xz_yz = cbuffer.data(dd_geom_20_off + 160 * ccomps * dcomps);

            auto g_yz_0_xz_zz = cbuffer.data(dd_geom_20_off + 161 * ccomps * dcomps);

            auto g_yz_0_yy_xx = cbuffer.data(dd_geom_20_off + 162 * ccomps * dcomps);

            auto g_yz_0_yy_xy = cbuffer.data(dd_geom_20_off + 163 * ccomps * dcomps);

            auto g_yz_0_yy_xz = cbuffer.data(dd_geom_20_off + 164 * ccomps * dcomps);

            auto g_yz_0_yy_yy = cbuffer.data(dd_geom_20_off + 165 * ccomps * dcomps);

            auto g_yz_0_yy_yz = cbuffer.data(dd_geom_20_off + 166 * ccomps * dcomps);

            auto g_yz_0_yy_zz = cbuffer.data(dd_geom_20_off + 167 * ccomps * dcomps);

            auto g_yz_0_yz_xx = cbuffer.data(dd_geom_20_off + 168 * ccomps * dcomps);

            auto g_yz_0_yz_xy = cbuffer.data(dd_geom_20_off + 169 * ccomps * dcomps);

            auto g_yz_0_yz_xz = cbuffer.data(dd_geom_20_off + 170 * ccomps * dcomps);

            auto g_yz_0_yz_yy = cbuffer.data(dd_geom_20_off + 171 * ccomps * dcomps);

            auto g_yz_0_yz_yz = cbuffer.data(dd_geom_20_off + 172 * ccomps * dcomps);

            auto g_yz_0_yz_zz = cbuffer.data(dd_geom_20_off + 173 * ccomps * dcomps);

            auto g_yz_0_zz_xx = cbuffer.data(dd_geom_20_off + 174 * ccomps * dcomps);

            auto g_yz_0_zz_xy = cbuffer.data(dd_geom_20_off + 175 * ccomps * dcomps);

            auto g_yz_0_zz_xz = cbuffer.data(dd_geom_20_off + 176 * ccomps * dcomps);

            auto g_yz_0_zz_yy = cbuffer.data(dd_geom_20_off + 177 * ccomps * dcomps);

            auto g_yz_0_zz_yz = cbuffer.data(dd_geom_20_off + 178 * ccomps * dcomps);

            auto g_yz_0_zz_zz = cbuffer.data(dd_geom_20_off + 179 * ccomps * dcomps);

            auto g_zz_0_xx_xx = cbuffer.data(dd_geom_20_off + 180 * ccomps * dcomps);

            auto g_zz_0_xx_xy = cbuffer.data(dd_geom_20_off + 181 * ccomps * dcomps);

            auto g_zz_0_xx_xz = cbuffer.data(dd_geom_20_off + 182 * ccomps * dcomps);

            auto g_zz_0_xx_yy = cbuffer.data(dd_geom_20_off + 183 * ccomps * dcomps);

            auto g_zz_0_xx_yz = cbuffer.data(dd_geom_20_off + 184 * ccomps * dcomps);

            auto g_zz_0_xx_zz = cbuffer.data(dd_geom_20_off + 185 * ccomps * dcomps);

            auto g_zz_0_xy_xx = cbuffer.data(dd_geom_20_off + 186 * ccomps * dcomps);

            auto g_zz_0_xy_xy = cbuffer.data(dd_geom_20_off + 187 * ccomps * dcomps);

            auto g_zz_0_xy_xz = cbuffer.data(dd_geom_20_off + 188 * ccomps * dcomps);

            auto g_zz_0_xy_yy = cbuffer.data(dd_geom_20_off + 189 * ccomps * dcomps);

            auto g_zz_0_xy_yz = cbuffer.data(dd_geom_20_off + 190 * ccomps * dcomps);

            auto g_zz_0_xy_zz = cbuffer.data(dd_geom_20_off + 191 * ccomps * dcomps);

            auto g_zz_0_xz_xx = cbuffer.data(dd_geom_20_off + 192 * ccomps * dcomps);

            auto g_zz_0_xz_xy = cbuffer.data(dd_geom_20_off + 193 * ccomps * dcomps);

            auto g_zz_0_xz_xz = cbuffer.data(dd_geom_20_off + 194 * ccomps * dcomps);

            auto g_zz_0_xz_yy = cbuffer.data(dd_geom_20_off + 195 * ccomps * dcomps);

            auto g_zz_0_xz_yz = cbuffer.data(dd_geom_20_off + 196 * ccomps * dcomps);

            auto g_zz_0_xz_zz = cbuffer.data(dd_geom_20_off + 197 * ccomps * dcomps);

            auto g_zz_0_yy_xx = cbuffer.data(dd_geom_20_off + 198 * ccomps * dcomps);

            auto g_zz_0_yy_xy = cbuffer.data(dd_geom_20_off + 199 * ccomps * dcomps);

            auto g_zz_0_yy_xz = cbuffer.data(dd_geom_20_off + 200 * ccomps * dcomps);

            auto g_zz_0_yy_yy = cbuffer.data(dd_geom_20_off + 201 * ccomps * dcomps);

            auto g_zz_0_yy_yz = cbuffer.data(dd_geom_20_off + 202 * ccomps * dcomps);

            auto g_zz_0_yy_zz = cbuffer.data(dd_geom_20_off + 203 * ccomps * dcomps);

            auto g_zz_0_yz_xx = cbuffer.data(dd_geom_20_off + 204 * ccomps * dcomps);

            auto g_zz_0_yz_xy = cbuffer.data(dd_geom_20_off + 205 * ccomps * dcomps);

            auto g_zz_0_yz_xz = cbuffer.data(dd_geom_20_off + 206 * ccomps * dcomps);

            auto g_zz_0_yz_yy = cbuffer.data(dd_geom_20_off + 207 * ccomps * dcomps);

            auto g_zz_0_yz_yz = cbuffer.data(dd_geom_20_off + 208 * ccomps * dcomps);

            auto g_zz_0_yz_zz = cbuffer.data(dd_geom_20_off + 209 * ccomps * dcomps);

            auto g_zz_0_zz_xx = cbuffer.data(dd_geom_20_off + 210 * ccomps * dcomps);

            auto g_zz_0_zz_xy = cbuffer.data(dd_geom_20_off + 211 * ccomps * dcomps);

            auto g_zz_0_zz_xz = cbuffer.data(dd_geom_20_off + 212 * ccomps * dcomps);

            auto g_zz_0_zz_yy = cbuffer.data(dd_geom_20_off + 213 * ccomps * dcomps);

            auto g_zz_0_zz_yz = cbuffer.data(dd_geom_20_off + 214 * ccomps * dcomps);

            auto g_zz_0_zz_zz = cbuffer.data(dd_geom_20_off + 215 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fpxx

            const auto fp_geom_20_off = idx_geom_20_fpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxx_x = cbuffer.data(fp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_y = cbuffer.data(fp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_z = cbuffer.data(fp_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_x, g_x_0_xx_y, g_x_0_xx_z, g_xx_0_xx_x, g_xx_0_xx_xx, g_xx_0_xx_xy, g_xx_0_xx_xz, g_xx_0_xx_y, g_xx_0_xx_z, g_xx_0_xxx_x, g_xx_0_xxx_y, g_xx_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxx_x[k] = -2.0 * g_x_0_xx_x[k] - g_xx_0_xx_x[k] * ab_x + g_xx_0_xx_xx[k];

                g_xx_0_xxx_y[k] = -2.0 * g_x_0_xx_y[k] - g_xx_0_xx_y[k] * ab_x + g_xx_0_xx_xy[k];

                g_xx_0_xxx_z[k] = -2.0 * g_x_0_xx_z[k] - g_xx_0_xx_z[k] * ab_x + g_xx_0_xx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxy_x = cbuffer.data(fp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxy_y = cbuffer.data(fp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxy_z = cbuffer.data(fp_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_x, g_xx_0_xx_xy, g_xx_0_xx_y, g_xx_0_xx_yy, g_xx_0_xx_yz, g_xx_0_xx_z, g_xx_0_xxy_x, g_xx_0_xxy_y, g_xx_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxy_x[k] = -g_xx_0_xx_x[k] * ab_y + g_xx_0_xx_xy[k];

                g_xx_0_xxy_y[k] = -g_xx_0_xx_y[k] * ab_y + g_xx_0_xx_yy[k];

                g_xx_0_xxy_z[k] = -g_xx_0_xx_z[k] * ab_y + g_xx_0_xx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxz_x = cbuffer.data(fp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxz_y = cbuffer.data(fp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxz_z = cbuffer.data(fp_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_x, g_xx_0_xx_xz, g_xx_0_xx_y, g_xx_0_xx_yz, g_xx_0_xx_z, g_xx_0_xx_zz, g_xx_0_xxz_x, g_xx_0_xxz_y, g_xx_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxz_x[k] = -g_xx_0_xx_x[k] * ab_z + g_xx_0_xx_xz[k];

                g_xx_0_xxz_y[k] = -g_xx_0_xx_y[k] * ab_z + g_xx_0_xx_yz[k];

                g_xx_0_xxz_z[k] = -g_xx_0_xx_z[k] * ab_z + g_xx_0_xx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyy_x = cbuffer.data(fp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xyy_y = cbuffer.data(fp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xyy_z = cbuffer.data(fp_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xy_x, g_xx_0_xy_xy, g_xx_0_xy_y, g_xx_0_xy_yy, g_xx_0_xy_yz, g_xx_0_xy_z, g_xx_0_xyy_x, g_xx_0_xyy_y, g_xx_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyy_x[k] = -g_xx_0_xy_x[k] * ab_y + g_xx_0_xy_xy[k];

                g_xx_0_xyy_y[k] = -g_xx_0_xy_y[k] * ab_y + g_xx_0_xy_yy[k];

                g_xx_0_xyy_z[k] = -g_xx_0_xy_z[k] * ab_y + g_xx_0_xy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyz_x = cbuffer.data(fp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xyz_y = cbuffer.data(fp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xyz_z = cbuffer.data(fp_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyz_x, g_xx_0_xyz_y, g_xx_0_xyz_z, g_xx_0_xz_x, g_xx_0_xz_xy, g_xx_0_xz_y, g_xx_0_xz_yy, g_xx_0_xz_yz, g_xx_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyz_x[k] = -g_xx_0_xz_x[k] * ab_y + g_xx_0_xz_xy[k];

                g_xx_0_xyz_y[k] = -g_xx_0_xz_y[k] * ab_y + g_xx_0_xz_yy[k];

                g_xx_0_xyz_z[k] = -g_xx_0_xz_z[k] * ab_y + g_xx_0_xz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzz_x = cbuffer.data(fp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xzz_y = cbuffer.data(fp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xzz_z = cbuffer.data(fp_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xz_x, g_xx_0_xz_xz, g_xx_0_xz_y, g_xx_0_xz_yz, g_xx_0_xz_z, g_xx_0_xz_zz, g_xx_0_xzz_x, g_xx_0_xzz_y, g_xx_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzz_x[k] = -g_xx_0_xz_x[k] * ab_z + g_xx_0_xz_xz[k];

                g_xx_0_xzz_y[k] = -g_xx_0_xz_y[k] * ab_z + g_xx_0_xz_yz[k];

                g_xx_0_xzz_z[k] = -g_xx_0_xz_z[k] * ab_z + g_xx_0_xz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyy_x = cbuffer.data(fp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_yyy_y = cbuffer.data(fp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_yyy_z = cbuffer.data(fp_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yy_x, g_xx_0_yy_xy, g_xx_0_yy_y, g_xx_0_yy_yy, g_xx_0_yy_yz, g_xx_0_yy_z, g_xx_0_yyy_x, g_xx_0_yyy_y, g_xx_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyy_x[k] = -g_xx_0_yy_x[k] * ab_y + g_xx_0_yy_xy[k];

                g_xx_0_yyy_y[k] = -g_xx_0_yy_y[k] * ab_y + g_xx_0_yy_yy[k];

                g_xx_0_yyy_z[k] = -g_xx_0_yy_z[k] * ab_y + g_xx_0_yy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyz_x = cbuffer.data(fp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_yyz_y = cbuffer.data(fp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_yyz_z = cbuffer.data(fp_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyz_x, g_xx_0_yyz_y, g_xx_0_yyz_z, g_xx_0_yz_x, g_xx_0_yz_xy, g_xx_0_yz_y, g_xx_0_yz_yy, g_xx_0_yz_yz, g_xx_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyz_x[k] = -g_xx_0_yz_x[k] * ab_y + g_xx_0_yz_xy[k];

                g_xx_0_yyz_y[k] = -g_xx_0_yz_y[k] * ab_y + g_xx_0_yz_yy[k];

                g_xx_0_yyz_z[k] = -g_xx_0_yz_z[k] * ab_y + g_xx_0_yz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzz_x = cbuffer.data(fp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_yzz_y = cbuffer.data(fp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_yzz_z = cbuffer.data(fp_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzz_x, g_xx_0_yzz_y, g_xx_0_yzz_z, g_xx_0_zz_x, g_xx_0_zz_xy, g_xx_0_zz_y, g_xx_0_zz_yy, g_xx_0_zz_yz, g_xx_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzz_x[k] = -g_xx_0_zz_x[k] * ab_y + g_xx_0_zz_xy[k];

                g_xx_0_yzz_y[k] = -g_xx_0_zz_y[k] * ab_y + g_xx_0_zz_yy[k];

                g_xx_0_yzz_z[k] = -g_xx_0_zz_z[k] * ab_y + g_xx_0_zz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzz_x = cbuffer.data(fp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_zzz_y = cbuffer.data(fp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_zzz_z = cbuffer.data(fp_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zz_x, g_xx_0_zz_xz, g_xx_0_zz_y, g_xx_0_zz_yz, g_xx_0_zz_z, g_xx_0_zz_zz, g_xx_0_zzz_x, g_xx_0_zzz_y, g_xx_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzz_x[k] = -g_xx_0_zz_x[k] * ab_z + g_xx_0_zz_xz[k];

                g_xx_0_zzz_y[k] = -g_xx_0_zz_y[k] * ab_z + g_xx_0_zz_yz[k];

                g_xx_0_zzz_z[k] = -g_xx_0_zz_z[k] * ab_z + g_xx_0_zz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxx_x = cbuffer.data(fp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_xxx_y = cbuffer.data(fp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_xxx_z = cbuffer.data(fp_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_x, g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_y, g_xy_0_xx_z, g_xy_0_xxx_x, g_xy_0_xxx_y, g_xy_0_xxx_z, g_y_0_xx_x, g_y_0_xx_y, g_y_0_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxx_x[k] = -g_y_0_xx_x[k] - g_xy_0_xx_x[k] * ab_x + g_xy_0_xx_xx[k];

                g_xy_0_xxx_y[k] = -g_y_0_xx_y[k] - g_xy_0_xx_y[k] * ab_x + g_xy_0_xx_xy[k];

                g_xy_0_xxx_z[k] = -g_y_0_xx_z[k] - g_xy_0_xx_z[k] * ab_x + g_xy_0_xx_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxy_x = cbuffer.data(fp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_xxy_y = cbuffer.data(fp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_xxy_z = cbuffer.data(fp_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_x, g_xy_0_xxy_y, g_xy_0_xxy_z, g_xy_0_xy_x, g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_y, g_xy_0_xy_z, g_y_0_xy_x, g_y_0_xy_y, g_y_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxy_x[k] = -g_y_0_xy_x[k] - g_xy_0_xy_x[k] * ab_x + g_xy_0_xy_xx[k];

                g_xy_0_xxy_y[k] = -g_y_0_xy_y[k] - g_xy_0_xy_y[k] * ab_x + g_xy_0_xy_xy[k];

                g_xy_0_xxy_z[k] = -g_y_0_xy_z[k] - g_xy_0_xy_z[k] * ab_x + g_xy_0_xy_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxz_x = cbuffer.data(fp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_xxz_y = cbuffer.data(fp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_xxz_z = cbuffer.data(fp_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_x, g_xy_0_xx_xz, g_xy_0_xx_y, g_xy_0_xx_yz, g_xy_0_xx_z, g_xy_0_xx_zz, g_xy_0_xxz_x, g_xy_0_xxz_y, g_xy_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxz_x[k] = -g_xy_0_xx_x[k] * ab_z + g_xy_0_xx_xz[k];

                g_xy_0_xxz_y[k] = -g_xy_0_xx_y[k] * ab_z + g_xy_0_xx_yz[k];

                g_xy_0_xxz_z[k] = -g_xy_0_xx_z[k] * ab_z + g_xy_0_xx_zz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyy_x = cbuffer.data(fp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_xyy_y = cbuffer.data(fp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_xyy_z = cbuffer.data(fp_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_x, g_xy_0_xyy_y, g_xy_0_xyy_z, g_xy_0_yy_x, g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_y, g_xy_0_yy_z, g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyy_x[k] = -g_y_0_yy_x[k] - g_xy_0_yy_x[k] * ab_x + g_xy_0_yy_xx[k];

                g_xy_0_xyy_y[k] = -g_y_0_yy_y[k] - g_xy_0_yy_y[k] * ab_x + g_xy_0_yy_xy[k];

                g_xy_0_xyy_z[k] = -g_y_0_yy_z[k] - g_xy_0_yy_z[k] * ab_x + g_xy_0_yy_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyz_x = cbuffer.data(fp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_xyz_y = cbuffer.data(fp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_xyz_z = cbuffer.data(fp_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_x, g_xy_0_xy_xz, g_xy_0_xy_y, g_xy_0_xy_yz, g_xy_0_xy_z, g_xy_0_xy_zz, g_xy_0_xyz_x, g_xy_0_xyz_y, g_xy_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyz_x[k] = -g_xy_0_xy_x[k] * ab_z + g_xy_0_xy_xz[k];

                g_xy_0_xyz_y[k] = -g_xy_0_xy_y[k] * ab_z + g_xy_0_xy_yz[k];

                g_xy_0_xyz_z[k] = -g_xy_0_xy_z[k] * ab_z + g_xy_0_xy_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzz_x = cbuffer.data(fp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_xzz_y = cbuffer.data(fp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_xzz_z = cbuffer.data(fp_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xz_x, g_xy_0_xz_xz, g_xy_0_xz_y, g_xy_0_xz_yz, g_xy_0_xz_z, g_xy_0_xz_zz, g_xy_0_xzz_x, g_xy_0_xzz_y, g_xy_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzz_x[k] = -g_xy_0_xz_x[k] * ab_z + g_xy_0_xz_xz[k];

                g_xy_0_xzz_y[k] = -g_xy_0_xz_y[k] * ab_z + g_xy_0_xz_yz[k];

                g_xy_0_xzz_z[k] = -g_xy_0_xz_z[k] * ab_z + g_xy_0_xz_zz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyy_x = cbuffer.data(fp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_yyy_y = cbuffer.data(fp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_yyy_z = cbuffer.data(fp_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_x, g_x_0_yy_y, g_x_0_yy_z, g_xy_0_yy_x, g_xy_0_yy_xy, g_xy_0_yy_y, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_z, g_xy_0_yyy_x, g_xy_0_yyy_y, g_xy_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyy_x[k] = -g_x_0_yy_x[k] - g_xy_0_yy_x[k] * ab_y + g_xy_0_yy_xy[k];

                g_xy_0_yyy_y[k] = -g_x_0_yy_y[k] - g_xy_0_yy_y[k] * ab_y + g_xy_0_yy_yy[k];

                g_xy_0_yyy_z[k] = -g_x_0_yy_z[k] - g_xy_0_yy_z[k] * ab_y + g_xy_0_yy_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyz_x = cbuffer.data(fp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_yyz_y = cbuffer.data(fp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_yyz_z = cbuffer.data(fp_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yy_x, g_xy_0_yy_xz, g_xy_0_yy_y, g_xy_0_yy_yz, g_xy_0_yy_z, g_xy_0_yy_zz, g_xy_0_yyz_x, g_xy_0_yyz_y, g_xy_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyz_x[k] = -g_xy_0_yy_x[k] * ab_z + g_xy_0_yy_xz[k];

                g_xy_0_yyz_y[k] = -g_xy_0_yy_y[k] * ab_z + g_xy_0_yy_yz[k];

                g_xy_0_yyz_z[k] = -g_xy_0_yy_z[k] * ab_z + g_xy_0_yy_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzz_x = cbuffer.data(fp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_yzz_y = cbuffer.data(fp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_yzz_z = cbuffer.data(fp_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yz_x, g_xy_0_yz_xz, g_xy_0_yz_y, g_xy_0_yz_yz, g_xy_0_yz_z, g_xy_0_yz_zz, g_xy_0_yzz_x, g_xy_0_yzz_y, g_xy_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzz_x[k] = -g_xy_0_yz_x[k] * ab_z + g_xy_0_yz_xz[k];

                g_xy_0_yzz_y[k] = -g_xy_0_yz_y[k] * ab_z + g_xy_0_yz_yz[k];

                g_xy_0_yzz_z[k] = -g_xy_0_yz_z[k] * ab_z + g_xy_0_yz_zz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzz_x = cbuffer.data(fp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_zzz_y = cbuffer.data(fp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_zzz_z = cbuffer.data(fp_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zz_x, g_xy_0_zz_xz, g_xy_0_zz_y, g_xy_0_zz_yz, g_xy_0_zz_z, g_xy_0_zz_zz, g_xy_0_zzz_x, g_xy_0_zzz_y, g_xy_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzz_x[k] = -g_xy_0_zz_x[k] * ab_z + g_xy_0_zz_xz[k];

                g_xy_0_zzz_y[k] = -g_xy_0_zz_y[k] * ab_z + g_xy_0_zz_yz[k];

                g_xy_0_zzz_z[k] = -g_xy_0_zz_z[k] * ab_z + g_xy_0_zz_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxx_x = cbuffer.data(fp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_xxx_y = cbuffer.data(fp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_xxx_z = cbuffer.data(fp_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_x, g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_y, g_xz_0_xx_z, g_xz_0_xxx_x, g_xz_0_xxx_y, g_xz_0_xxx_z, g_z_0_xx_x, g_z_0_xx_y, g_z_0_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxx_x[k] = -g_z_0_xx_x[k] - g_xz_0_xx_x[k] * ab_x + g_xz_0_xx_xx[k];

                g_xz_0_xxx_y[k] = -g_z_0_xx_y[k] - g_xz_0_xx_y[k] * ab_x + g_xz_0_xx_xy[k];

                g_xz_0_xxx_z[k] = -g_z_0_xx_z[k] - g_xz_0_xx_z[k] * ab_x + g_xz_0_xx_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxy_x = cbuffer.data(fp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_xxy_y = cbuffer.data(fp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_xxy_z = cbuffer.data(fp_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_x, g_xz_0_xx_xy, g_xz_0_xx_y, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_z, g_xz_0_xxy_x, g_xz_0_xxy_y, g_xz_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxy_x[k] = -g_xz_0_xx_x[k] * ab_y + g_xz_0_xx_xy[k];

                g_xz_0_xxy_y[k] = -g_xz_0_xx_y[k] * ab_y + g_xz_0_xx_yy[k];

                g_xz_0_xxy_z[k] = -g_xz_0_xx_z[k] * ab_y + g_xz_0_xx_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxz_x = cbuffer.data(fp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_xxz_y = cbuffer.data(fp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_xxz_z = cbuffer.data(fp_geom_20_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxz_x, g_xz_0_xxz_y, g_xz_0_xxz_z, g_xz_0_xz_x, g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_y, g_xz_0_xz_z, g_z_0_xz_x, g_z_0_xz_y, g_z_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxz_x[k] = -g_z_0_xz_x[k] - g_xz_0_xz_x[k] * ab_x + g_xz_0_xz_xx[k];

                g_xz_0_xxz_y[k] = -g_z_0_xz_y[k] - g_xz_0_xz_y[k] * ab_x + g_xz_0_xz_xy[k];

                g_xz_0_xxz_z[k] = -g_z_0_xz_z[k] - g_xz_0_xz_z[k] * ab_x + g_xz_0_xz_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyy_x = cbuffer.data(fp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xz_0_xyy_y = cbuffer.data(fp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_xyy_z = cbuffer.data(fp_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xy_x, g_xz_0_xy_xy, g_xz_0_xy_y, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_z, g_xz_0_xyy_x, g_xz_0_xyy_y, g_xz_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyy_x[k] = -g_xz_0_xy_x[k] * ab_y + g_xz_0_xy_xy[k];

                g_xz_0_xyy_y[k] = -g_xz_0_xy_y[k] * ab_y + g_xz_0_xy_yy[k];

                g_xz_0_xyy_z[k] = -g_xz_0_xy_z[k] * ab_y + g_xz_0_xy_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyz_x = cbuffer.data(fp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_xyz_y = cbuffer.data(fp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_xyz_z = cbuffer.data(fp_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyz_x, g_xz_0_xyz_y, g_xz_0_xyz_z, g_xz_0_xz_x, g_xz_0_xz_xy, g_xz_0_xz_y, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyz_x[k] = -g_xz_0_xz_x[k] * ab_y + g_xz_0_xz_xy[k];

                g_xz_0_xyz_y[k] = -g_xz_0_xz_y[k] * ab_y + g_xz_0_xz_yy[k];

                g_xz_0_xyz_z[k] = -g_xz_0_xz_z[k] * ab_y + g_xz_0_xz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzz_x = cbuffer.data(fp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_xzz_y = cbuffer.data(fp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_xzz_z = cbuffer.data(fp_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzz_x, g_xz_0_xzz_y, g_xz_0_xzz_z, g_xz_0_zz_x, g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_y, g_xz_0_zz_z, g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzz_x[k] = -g_z_0_zz_x[k] - g_xz_0_zz_x[k] * ab_x + g_xz_0_zz_xx[k];

                g_xz_0_xzz_y[k] = -g_z_0_zz_y[k] - g_xz_0_zz_y[k] * ab_x + g_xz_0_zz_xy[k];

                g_xz_0_xzz_z[k] = -g_z_0_zz_z[k] - g_xz_0_zz_z[k] * ab_x + g_xz_0_zz_xz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyy_x = cbuffer.data(fp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_yyy_y = cbuffer.data(fp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_yyy_z = cbuffer.data(fp_geom_20_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yy_x, g_xz_0_yy_xy, g_xz_0_yy_y, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_z, g_xz_0_yyy_x, g_xz_0_yyy_y, g_xz_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyy_x[k] = -g_xz_0_yy_x[k] * ab_y + g_xz_0_yy_xy[k];

                g_xz_0_yyy_y[k] = -g_xz_0_yy_y[k] * ab_y + g_xz_0_yy_yy[k];

                g_xz_0_yyy_z[k] = -g_xz_0_yy_z[k] * ab_y + g_xz_0_yy_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyz_x = cbuffer.data(fp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_yyz_y = cbuffer.data(fp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_yyz_z = cbuffer.data(fp_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyz_x, g_xz_0_yyz_y, g_xz_0_yyz_z, g_xz_0_yz_x, g_xz_0_yz_xy, g_xz_0_yz_y, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyz_x[k] = -g_xz_0_yz_x[k] * ab_y + g_xz_0_yz_xy[k];

                g_xz_0_yyz_y[k] = -g_xz_0_yz_y[k] * ab_y + g_xz_0_yz_yy[k];

                g_xz_0_yyz_z[k] = -g_xz_0_yz_z[k] * ab_y + g_xz_0_yz_yz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzz_x = cbuffer.data(fp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_yzz_y = cbuffer.data(fp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_yzz_z = cbuffer.data(fp_geom_20_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzz_x, g_xz_0_yzz_y, g_xz_0_yzz_z, g_xz_0_zz_x, g_xz_0_zz_xy, g_xz_0_zz_y, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzz_x[k] = -g_xz_0_zz_x[k] * ab_y + g_xz_0_zz_xy[k];

                g_xz_0_yzz_y[k] = -g_xz_0_zz_y[k] * ab_y + g_xz_0_zz_yy[k];

                g_xz_0_yzz_z[k] = -g_xz_0_zz_z[k] * ab_y + g_xz_0_zz_yz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzz_x = cbuffer.data(fp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_zzz_y = cbuffer.data(fp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_zzz_z = cbuffer.data(fp_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_x, g_x_0_zz_y, g_x_0_zz_z, g_xz_0_zz_x, g_xz_0_zz_xz, g_xz_0_zz_y, g_xz_0_zz_yz, g_xz_0_zz_z, g_xz_0_zz_zz, g_xz_0_zzz_x, g_xz_0_zzz_y, g_xz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzz_x[k] = -g_x_0_zz_x[k] - g_xz_0_zz_x[k] * ab_z + g_xz_0_zz_xz[k];

                g_xz_0_zzz_y[k] = -g_x_0_zz_y[k] - g_xz_0_zz_y[k] * ab_z + g_xz_0_zz_yz[k];

                g_xz_0_zzz_z[k] = -g_x_0_zz_z[k] - g_xz_0_zz_z[k] * ab_z + g_xz_0_zz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxx_x = cbuffer.data(fp_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_xxx_y = cbuffer.data(fp_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_xxx_z = cbuffer.data(fp_geom_20_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xx_x, g_yy_0_xx_xx, g_yy_0_xx_xy, g_yy_0_xx_xz, g_yy_0_xx_y, g_yy_0_xx_z, g_yy_0_xxx_x, g_yy_0_xxx_y, g_yy_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxx_x[k] = -g_yy_0_xx_x[k] * ab_x + g_yy_0_xx_xx[k];

                g_yy_0_xxx_y[k] = -g_yy_0_xx_y[k] * ab_x + g_yy_0_xx_xy[k];

                g_yy_0_xxx_z[k] = -g_yy_0_xx_z[k] * ab_x + g_yy_0_xx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxy_x = cbuffer.data(fp_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_xxy_y = cbuffer.data(fp_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_xxy_z = cbuffer.data(fp_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxy_x, g_yy_0_xxy_y, g_yy_0_xxy_z, g_yy_0_xy_x, g_yy_0_xy_xx, g_yy_0_xy_xy, g_yy_0_xy_xz, g_yy_0_xy_y, g_yy_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxy_x[k] = -g_yy_0_xy_x[k] * ab_x + g_yy_0_xy_xx[k];

                g_yy_0_xxy_y[k] = -g_yy_0_xy_y[k] * ab_x + g_yy_0_xy_xy[k];

                g_yy_0_xxy_z[k] = -g_yy_0_xy_z[k] * ab_x + g_yy_0_xy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxz_x = cbuffer.data(fp_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_xxz_y = cbuffer.data(fp_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_xxz_z = cbuffer.data(fp_geom_20_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxz_x, g_yy_0_xxz_y, g_yy_0_xxz_z, g_yy_0_xz_x, g_yy_0_xz_xx, g_yy_0_xz_xy, g_yy_0_xz_xz, g_yy_0_xz_y, g_yy_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxz_x[k] = -g_yy_0_xz_x[k] * ab_x + g_yy_0_xz_xx[k];

                g_yy_0_xxz_y[k] = -g_yy_0_xz_y[k] * ab_x + g_yy_0_xz_xy[k];

                g_yy_0_xxz_z[k] = -g_yy_0_xz_z[k] * ab_x + g_yy_0_xz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyy_x = cbuffer.data(fp_geom_20_off + 99 * ccomps * dcomps);

            auto g_yy_0_xyy_y = cbuffer.data(fp_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_xyy_z = cbuffer.data(fp_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyy_x, g_yy_0_xyy_y, g_yy_0_xyy_z, g_yy_0_yy_x, g_yy_0_yy_xx, g_yy_0_yy_xy, g_yy_0_yy_xz, g_yy_0_yy_y, g_yy_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyy_x[k] = -g_yy_0_yy_x[k] * ab_x + g_yy_0_yy_xx[k];

                g_yy_0_xyy_y[k] = -g_yy_0_yy_y[k] * ab_x + g_yy_0_yy_xy[k];

                g_yy_0_xyy_z[k] = -g_yy_0_yy_z[k] * ab_x + g_yy_0_yy_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyz_x = cbuffer.data(fp_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_xyz_y = cbuffer.data(fp_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_xyz_z = cbuffer.data(fp_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyz_x, g_yy_0_xyz_y, g_yy_0_xyz_z, g_yy_0_yz_x, g_yy_0_yz_xx, g_yy_0_yz_xy, g_yy_0_yz_xz, g_yy_0_yz_y, g_yy_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyz_x[k] = -g_yy_0_yz_x[k] * ab_x + g_yy_0_yz_xx[k];

                g_yy_0_xyz_y[k] = -g_yy_0_yz_y[k] * ab_x + g_yy_0_yz_xy[k];

                g_yy_0_xyz_z[k] = -g_yy_0_yz_z[k] * ab_x + g_yy_0_yz_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzz_x = cbuffer.data(fp_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_xzz_y = cbuffer.data(fp_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_xzz_z = cbuffer.data(fp_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzz_x, g_yy_0_xzz_y, g_yy_0_xzz_z, g_yy_0_zz_x, g_yy_0_zz_xx, g_yy_0_zz_xy, g_yy_0_zz_xz, g_yy_0_zz_y, g_yy_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzz_x[k] = -g_yy_0_zz_x[k] * ab_x + g_yy_0_zz_xx[k];

                g_yy_0_xzz_y[k] = -g_yy_0_zz_y[k] * ab_x + g_yy_0_zz_xy[k];

                g_yy_0_xzz_z[k] = -g_yy_0_zz_z[k] * ab_x + g_yy_0_zz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyy_x = cbuffer.data(fp_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_yyy_y = cbuffer.data(fp_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_yyy_z = cbuffer.data(fp_geom_20_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z, g_yy_0_yy_x, g_yy_0_yy_xy, g_yy_0_yy_y, g_yy_0_yy_yy, g_yy_0_yy_yz, g_yy_0_yy_z, g_yy_0_yyy_x, g_yy_0_yyy_y, g_yy_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyy_x[k] = -2.0 * g_y_0_yy_x[k] - g_yy_0_yy_x[k] * ab_y + g_yy_0_yy_xy[k];

                g_yy_0_yyy_y[k] = -2.0 * g_y_0_yy_y[k] - g_yy_0_yy_y[k] * ab_y + g_yy_0_yy_yy[k];

                g_yy_0_yyy_z[k] = -2.0 * g_y_0_yy_z[k] - g_yy_0_yy_z[k] * ab_y + g_yy_0_yy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyz_x = cbuffer.data(fp_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_yyz_y = cbuffer.data(fp_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_yyz_z = cbuffer.data(fp_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yy_x, g_yy_0_yy_xz, g_yy_0_yy_y, g_yy_0_yy_yz, g_yy_0_yy_z, g_yy_0_yy_zz, g_yy_0_yyz_x, g_yy_0_yyz_y, g_yy_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyz_x[k] = -g_yy_0_yy_x[k] * ab_z + g_yy_0_yy_xz[k];

                g_yy_0_yyz_y[k] = -g_yy_0_yy_y[k] * ab_z + g_yy_0_yy_yz[k];

                g_yy_0_yyz_z[k] = -g_yy_0_yy_z[k] * ab_z + g_yy_0_yy_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzz_x = cbuffer.data(fp_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_yzz_y = cbuffer.data(fp_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_yzz_z = cbuffer.data(fp_geom_20_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yz_x, g_yy_0_yz_xz, g_yy_0_yz_y, g_yy_0_yz_yz, g_yy_0_yz_z, g_yy_0_yz_zz, g_yy_0_yzz_x, g_yy_0_yzz_y, g_yy_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzz_x[k] = -g_yy_0_yz_x[k] * ab_z + g_yy_0_yz_xz[k];

                g_yy_0_yzz_y[k] = -g_yy_0_yz_y[k] * ab_z + g_yy_0_yz_yz[k];

                g_yy_0_yzz_z[k] = -g_yy_0_yz_z[k] * ab_z + g_yy_0_yz_zz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzz_x = cbuffer.data(fp_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_zzz_y = cbuffer.data(fp_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_zzz_z = cbuffer.data(fp_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zz_x, g_yy_0_zz_xz, g_yy_0_zz_y, g_yy_0_zz_yz, g_yy_0_zz_z, g_yy_0_zz_zz, g_yy_0_zzz_x, g_yy_0_zzz_y, g_yy_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzz_x[k] = -g_yy_0_zz_x[k] * ab_z + g_yy_0_zz_xz[k];

                g_yy_0_zzz_y[k] = -g_yy_0_zz_y[k] * ab_z + g_yy_0_zz_yz[k];

                g_yy_0_zzz_z[k] = -g_yy_0_zz_z[k] * ab_z + g_yy_0_zz_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxx_x = cbuffer.data(fp_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_xxx_y = cbuffer.data(fp_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_xxx_z = cbuffer.data(fp_geom_20_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xx_x, g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_y, g_yz_0_xx_z, g_yz_0_xxx_x, g_yz_0_xxx_y, g_yz_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxx_x[k] = -g_yz_0_xx_x[k] * ab_x + g_yz_0_xx_xx[k];

                g_yz_0_xxx_y[k] = -g_yz_0_xx_y[k] * ab_x + g_yz_0_xx_xy[k];

                g_yz_0_xxx_z[k] = -g_yz_0_xx_z[k] * ab_x + g_yz_0_xx_xz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxy_x = cbuffer.data(fp_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_xxy_y = cbuffer.data(fp_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_xxy_z = cbuffer.data(fp_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxy_x, g_yz_0_xxy_y, g_yz_0_xxy_z, g_yz_0_xy_x, g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_y, g_yz_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxy_x[k] = -g_yz_0_xy_x[k] * ab_x + g_yz_0_xy_xx[k];

                g_yz_0_xxy_y[k] = -g_yz_0_xy_y[k] * ab_x + g_yz_0_xy_xy[k];

                g_yz_0_xxy_z[k] = -g_yz_0_xy_z[k] * ab_x + g_yz_0_xy_xz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxz_x = cbuffer.data(fp_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_xxz_y = cbuffer.data(fp_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_xxz_z = cbuffer.data(fp_geom_20_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxz_x, g_yz_0_xxz_y, g_yz_0_xxz_z, g_yz_0_xz_x, g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_y, g_yz_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxz_x[k] = -g_yz_0_xz_x[k] * ab_x + g_yz_0_xz_xx[k];

                g_yz_0_xxz_y[k] = -g_yz_0_xz_y[k] * ab_x + g_yz_0_xz_xy[k];

                g_yz_0_xxz_z[k] = -g_yz_0_xz_z[k] * ab_x + g_yz_0_xz_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyy_x = cbuffer.data(fp_geom_20_off + 129 * ccomps * dcomps);

            auto g_yz_0_xyy_y = cbuffer.data(fp_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_xyy_z = cbuffer.data(fp_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyy_x, g_yz_0_xyy_y, g_yz_0_xyy_z, g_yz_0_yy_x, g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_y, g_yz_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyy_x[k] = -g_yz_0_yy_x[k] * ab_x + g_yz_0_yy_xx[k];

                g_yz_0_xyy_y[k] = -g_yz_0_yy_y[k] * ab_x + g_yz_0_yy_xy[k];

                g_yz_0_xyy_z[k] = -g_yz_0_yy_z[k] * ab_x + g_yz_0_yy_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyz_x = cbuffer.data(fp_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_xyz_y = cbuffer.data(fp_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_xyz_z = cbuffer.data(fp_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyz_x, g_yz_0_xyz_y, g_yz_0_xyz_z, g_yz_0_yz_x, g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_y, g_yz_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyz_x[k] = -g_yz_0_yz_x[k] * ab_x + g_yz_0_yz_xx[k];

                g_yz_0_xyz_y[k] = -g_yz_0_yz_y[k] * ab_x + g_yz_0_yz_xy[k];

                g_yz_0_xyz_z[k] = -g_yz_0_yz_z[k] * ab_x + g_yz_0_yz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzz_x = cbuffer.data(fp_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_xzz_y = cbuffer.data(fp_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_xzz_z = cbuffer.data(fp_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzz_x, g_yz_0_xzz_y, g_yz_0_xzz_z, g_yz_0_zz_x, g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_y, g_yz_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzz_x[k] = -g_yz_0_zz_x[k] * ab_x + g_yz_0_zz_xx[k];

                g_yz_0_xzz_y[k] = -g_yz_0_zz_y[k] * ab_x + g_yz_0_zz_xy[k];

                g_yz_0_xzz_z[k] = -g_yz_0_zz_z[k] * ab_x + g_yz_0_zz_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyy_x = cbuffer.data(fp_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_yyy_y = cbuffer.data(fp_geom_20_off + 139 * ccomps * dcomps);

            auto g_yz_0_yyy_z = cbuffer.data(fp_geom_20_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yy_x, g_yz_0_yy_xy, g_yz_0_yy_y, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_z, g_yz_0_yyy_x, g_yz_0_yyy_y, g_yz_0_yyy_z, g_z_0_yy_x, g_z_0_yy_y, g_z_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyy_x[k] = -g_z_0_yy_x[k] - g_yz_0_yy_x[k] * ab_y + g_yz_0_yy_xy[k];

                g_yz_0_yyy_y[k] = -g_z_0_yy_y[k] - g_yz_0_yy_y[k] * ab_y + g_yz_0_yy_yy[k];

                g_yz_0_yyy_z[k] = -g_z_0_yy_z[k] - g_yz_0_yy_z[k] * ab_y + g_yz_0_yy_yz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyz_x = cbuffer.data(fp_geom_20_off + 141 * ccomps * dcomps);

            auto g_yz_0_yyz_y = cbuffer.data(fp_geom_20_off + 142 * ccomps * dcomps);

            auto g_yz_0_yyz_z = cbuffer.data(fp_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyz_x, g_yz_0_yyz_y, g_yz_0_yyz_z, g_yz_0_yz_x, g_yz_0_yz_xy, g_yz_0_yz_y, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_z, g_z_0_yz_x, g_z_0_yz_y, g_z_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyz_x[k] = -g_z_0_yz_x[k] - g_yz_0_yz_x[k] * ab_y + g_yz_0_yz_xy[k];

                g_yz_0_yyz_y[k] = -g_z_0_yz_y[k] - g_yz_0_yz_y[k] * ab_y + g_yz_0_yz_yy[k];

                g_yz_0_yyz_z[k] = -g_z_0_yz_z[k] - g_yz_0_yz_z[k] * ab_y + g_yz_0_yz_yz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzz_x = cbuffer.data(fp_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_yzz_y = cbuffer.data(fp_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_yzz_z = cbuffer.data(fp_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzz_x, g_yz_0_yzz_y, g_yz_0_yzz_z, g_yz_0_zz_x, g_yz_0_zz_xy, g_yz_0_zz_y, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_z, g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzz_x[k] = -g_z_0_zz_x[k] - g_yz_0_zz_x[k] * ab_y + g_yz_0_zz_xy[k];

                g_yz_0_yzz_y[k] = -g_z_0_zz_y[k] - g_yz_0_zz_y[k] * ab_y + g_yz_0_zz_yy[k];

                g_yz_0_yzz_z[k] = -g_z_0_zz_z[k] - g_yz_0_zz_z[k] * ab_y + g_yz_0_zz_yz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzz_x = cbuffer.data(fp_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_zzz_y = cbuffer.data(fp_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_zzz_z = cbuffer.data(fp_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_x, g_y_0_zz_y, g_y_0_zz_z, g_yz_0_zz_x, g_yz_0_zz_xz, g_yz_0_zz_y, g_yz_0_zz_yz, g_yz_0_zz_z, g_yz_0_zz_zz, g_yz_0_zzz_x, g_yz_0_zzz_y, g_yz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzz_x[k] = -g_y_0_zz_x[k] - g_yz_0_zz_x[k] * ab_z + g_yz_0_zz_xz[k];

                g_yz_0_zzz_y[k] = -g_y_0_zz_y[k] - g_yz_0_zz_y[k] * ab_z + g_yz_0_zz_yz[k];

                g_yz_0_zzz_z[k] = -g_y_0_zz_z[k] - g_yz_0_zz_z[k] * ab_z + g_yz_0_zz_zz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxx_x = cbuffer.data(fp_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_xxx_y = cbuffer.data(fp_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_xxx_z = cbuffer.data(fp_geom_20_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xx_x, g_zz_0_xx_xx, g_zz_0_xx_xy, g_zz_0_xx_xz, g_zz_0_xx_y, g_zz_0_xx_z, g_zz_0_xxx_x, g_zz_0_xxx_y, g_zz_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxx_x[k] = -g_zz_0_xx_x[k] * ab_x + g_zz_0_xx_xx[k];

                g_zz_0_xxx_y[k] = -g_zz_0_xx_y[k] * ab_x + g_zz_0_xx_xy[k];

                g_zz_0_xxx_z[k] = -g_zz_0_xx_z[k] * ab_x + g_zz_0_xx_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxy_x = cbuffer.data(fp_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_xxy_y = cbuffer.data(fp_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_xxy_z = cbuffer.data(fp_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxy_x, g_zz_0_xxy_y, g_zz_0_xxy_z, g_zz_0_xy_x, g_zz_0_xy_xx, g_zz_0_xy_xy, g_zz_0_xy_xz, g_zz_0_xy_y, g_zz_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxy_x[k] = -g_zz_0_xy_x[k] * ab_x + g_zz_0_xy_xx[k];

                g_zz_0_xxy_y[k] = -g_zz_0_xy_y[k] * ab_x + g_zz_0_xy_xy[k];

                g_zz_0_xxy_z[k] = -g_zz_0_xy_z[k] * ab_x + g_zz_0_xy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxz_x = cbuffer.data(fp_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_xxz_y = cbuffer.data(fp_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_xxz_z = cbuffer.data(fp_geom_20_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxz_x, g_zz_0_xxz_y, g_zz_0_xxz_z, g_zz_0_xz_x, g_zz_0_xz_xx, g_zz_0_xz_xy, g_zz_0_xz_xz, g_zz_0_xz_y, g_zz_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxz_x[k] = -g_zz_0_xz_x[k] * ab_x + g_zz_0_xz_xx[k];

                g_zz_0_xxz_y[k] = -g_zz_0_xz_y[k] * ab_x + g_zz_0_xz_xy[k];

                g_zz_0_xxz_z[k] = -g_zz_0_xz_z[k] * ab_x + g_zz_0_xz_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyy_x = cbuffer.data(fp_geom_20_off + 159 * ccomps * dcomps);

            auto g_zz_0_xyy_y = cbuffer.data(fp_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_xyy_z = cbuffer.data(fp_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyy_x, g_zz_0_xyy_y, g_zz_0_xyy_z, g_zz_0_yy_x, g_zz_0_yy_xx, g_zz_0_yy_xy, g_zz_0_yy_xz, g_zz_0_yy_y, g_zz_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyy_x[k] = -g_zz_0_yy_x[k] * ab_x + g_zz_0_yy_xx[k];

                g_zz_0_xyy_y[k] = -g_zz_0_yy_y[k] * ab_x + g_zz_0_yy_xy[k];

                g_zz_0_xyy_z[k] = -g_zz_0_yy_z[k] * ab_x + g_zz_0_yy_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyz_x = cbuffer.data(fp_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_xyz_y = cbuffer.data(fp_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_xyz_z = cbuffer.data(fp_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyz_x, g_zz_0_xyz_y, g_zz_0_xyz_z, g_zz_0_yz_x, g_zz_0_yz_xx, g_zz_0_yz_xy, g_zz_0_yz_xz, g_zz_0_yz_y, g_zz_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyz_x[k] = -g_zz_0_yz_x[k] * ab_x + g_zz_0_yz_xx[k];

                g_zz_0_xyz_y[k] = -g_zz_0_yz_y[k] * ab_x + g_zz_0_yz_xy[k];

                g_zz_0_xyz_z[k] = -g_zz_0_yz_z[k] * ab_x + g_zz_0_yz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzz_x = cbuffer.data(fp_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_xzz_y = cbuffer.data(fp_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_xzz_z = cbuffer.data(fp_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzz_x, g_zz_0_xzz_y, g_zz_0_xzz_z, g_zz_0_zz_x, g_zz_0_zz_xx, g_zz_0_zz_xy, g_zz_0_zz_xz, g_zz_0_zz_y, g_zz_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzz_x[k] = -g_zz_0_zz_x[k] * ab_x + g_zz_0_zz_xx[k];

                g_zz_0_xzz_y[k] = -g_zz_0_zz_y[k] * ab_x + g_zz_0_zz_xy[k];

                g_zz_0_xzz_z[k] = -g_zz_0_zz_z[k] * ab_x + g_zz_0_zz_xz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyy_x = cbuffer.data(fp_geom_20_off + 168 * ccomps * dcomps);

            auto g_zz_0_yyy_y = cbuffer.data(fp_geom_20_off + 169 * ccomps * dcomps);

            auto g_zz_0_yyy_z = cbuffer.data(fp_geom_20_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yy_x, g_zz_0_yy_xy, g_zz_0_yy_y, g_zz_0_yy_yy, g_zz_0_yy_yz, g_zz_0_yy_z, g_zz_0_yyy_x, g_zz_0_yyy_y, g_zz_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyy_x[k] = -g_zz_0_yy_x[k] * ab_y + g_zz_0_yy_xy[k];

                g_zz_0_yyy_y[k] = -g_zz_0_yy_y[k] * ab_y + g_zz_0_yy_yy[k];

                g_zz_0_yyy_z[k] = -g_zz_0_yy_z[k] * ab_y + g_zz_0_yy_yz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyz_x = cbuffer.data(fp_geom_20_off + 171 * ccomps * dcomps);

            auto g_zz_0_yyz_y = cbuffer.data(fp_geom_20_off + 172 * ccomps * dcomps);

            auto g_zz_0_yyz_z = cbuffer.data(fp_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyz_x, g_zz_0_yyz_y, g_zz_0_yyz_z, g_zz_0_yz_x, g_zz_0_yz_xy, g_zz_0_yz_y, g_zz_0_yz_yy, g_zz_0_yz_yz, g_zz_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyz_x[k] = -g_zz_0_yz_x[k] * ab_y + g_zz_0_yz_xy[k];

                g_zz_0_yyz_y[k] = -g_zz_0_yz_y[k] * ab_y + g_zz_0_yz_yy[k];

                g_zz_0_yyz_z[k] = -g_zz_0_yz_z[k] * ab_y + g_zz_0_yz_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzz_x = cbuffer.data(fp_geom_20_off + 174 * ccomps * dcomps);

            auto g_zz_0_yzz_y = cbuffer.data(fp_geom_20_off + 175 * ccomps * dcomps);

            auto g_zz_0_yzz_z = cbuffer.data(fp_geom_20_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzz_x, g_zz_0_yzz_y, g_zz_0_yzz_z, g_zz_0_zz_x, g_zz_0_zz_xy, g_zz_0_zz_y, g_zz_0_zz_yy, g_zz_0_zz_yz, g_zz_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzz_x[k] = -g_zz_0_zz_x[k] * ab_y + g_zz_0_zz_xy[k];

                g_zz_0_yzz_y[k] = -g_zz_0_zz_y[k] * ab_y + g_zz_0_zz_yy[k];

                g_zz_0_yzz_z[k] = -g_zz_0_zz_z[k] * ab_y + g_zz_0_zz_yz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzz_x = cbuffer.data(fp_geom_20_off + 177 * ccomps * dcomps);

            auto g_zz_0_zzz_y = cbuffer.data(fp_geom_20_off + 178 * ccomps * dcomps);

            auto g_zz_0_zzz_z = cbuffer.data(fp_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, g_zz_0_zz_x, g_zz_0_zz_xz, g_zz_0_zz_y, g_zz_0_zz_yz, g_zz_0_zz_z, g_zz_0_zz_zz, g_zz_0_zzz_x, g_zz_0_zzz_y, g_zz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzz_x[k] = -2.0 * g_z_0_zz_x[k] - g_zz_0_zz_x[k] * ab_z + g_zz_0_zz_xz[k];

                g_zz_0_zzz_y[k] = -2.0 * g_z_0_zz_y[k] - g_zz_0_zz_y[k] * ab_z + g_zz_0_zz_yz[k];

                g_zz_0_zzz_z[k] = -2.0 * g_z_0_zz_z[k] - g_zz_0_zz_z[k] * ab_z + g_zz_0_zz_zz[k];
            }
        }
    }
}

} // erirec namespace

