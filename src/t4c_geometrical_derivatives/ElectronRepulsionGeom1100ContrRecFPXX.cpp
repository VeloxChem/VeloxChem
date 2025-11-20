#include "ElectronRepulsionGeom1100ContrRecFPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_fpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_fpxx,
                                            const size_t idx_geom_01_dpxx,
                                            const size_t idx_geom_10_dpxx,
                                            const size_t idx_geom_11_dpxx,
                                            const size_t idx_geom_11_ddxx,
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

            const auto dp_geom_01_off = idx_geom_01_dpxx + i * dcomps + j;

            auto g_0_x_xx_x = cbuffer.data(dp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_y = cbuffer.data(dp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_z = cbuffer.data(dp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xy_x = cbuffer.data(dp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xy_y = cbuffer.data(dp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xy_z = cbuffer.data(dp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xz_x = cbuffer.data(dp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xz_y = cbuffer.data(dp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xz_z = cbuffer.data(dp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_yy_x = cbuffer.data(dp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_yy_y = cbuffer.data(dp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_yy_z = cbuffer.data(dp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_yz_x = cbuffer.data(dp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_yz_y = cbuffer.data(dp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_yz_z = cbuffer.data(dp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_zz_x = cbuffer.data(dp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_zz_y = cbuffer.data(dp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_zz_z = cbuffer.data(dp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_xx_x = cbuffer.data(dp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_xx_y = cbuffer.data(dp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_xx_z = cbuffer.data(dp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_xy_x = cbuffer.data(dp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_xy_y = cbuffer.data(dp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_xy_z = cbuffer.data(dp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_xz_x = cbuffer.data(dp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_xz_y = cbuffer.data(dp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_xz_z = cbuffer.data(dp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_yy_x = cbuffer.data(dp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_yy_y = cbuffer.data(dp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_yy_z = cbuffer.data(dp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_yz_x = cbuffer.data(dp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_yz_y = cbuffer.data(dp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_yz_z = cbuffer.data(dp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_zz_x = cbuffer.data(dp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_zz_y = cbuffer.data(dp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_zz_z = cbuffer.data(dp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_xx_x = cbuffer.data(dp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_xx_y = cbuffer.data(dp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_xx_z = cbuffer.data(dp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_xy_x = cbuffer.data(dp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_xy_y = cbuffer.data(dp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_xy_z = cbuffer.data(dp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_xz_x = cbuffer.data(dp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_xz_y = cbuffer.data(dp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_xz_z = cbuffer.data(dp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_yy_x = cbuffer.data(dp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_yy_y = cbuffer.data(dp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_yy_z = cbuffer.data(dp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_z_yz_x = cbuffer.data(dp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_yz_y = cbuffer.data(dp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_yz_z = cbuffer.data(dp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_zz_x = cbuffer.data(dp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_zz_y = cbuffer.data(dp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_zz_z = cbuffer.data(dp_geom_01_off + 53 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : DDSS

            const auto dd_geom_11_off = idx_geom_11_ddxx + i * dcomps + j;

            auto g_x_x_xx_xx = cbuffer.data(dd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xy = cbuffer.data(dd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xz = cbuffer.data(dd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_yy = cbuffer.data(dd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_yz = cbuffer.data(dd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_zz = cbuffer.data(dd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xy_xx = cbuffer.data(dd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xy_xy = cbuffer.data(dd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xy_xz = cbuffer.data(dd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xy_yy = cbuffer.data(dd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xy_yz = cbuffer.data(dd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xy_zz = cbuffer.data(dd_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xz_xx = cbuffer.data(dd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xz_xy = cbuffer.data(dd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xz_xz = cbuffer.data(dd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xz_yy = cbuffer.data(dd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xz_yz = cbuffer.data(dd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xz_zz = cbuffer.data(dd_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_yy_xx = cbuffer.data(dd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_yy_xy = cbuffer.data(dd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_yy_xz = cbuffer.data(dd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_yy_yy = cbuffer.data(dd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_yy_yz = cbuffer.data(dd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_yy_zz = cbuffer.data(dd_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_yz_xx = cbuffer.data(dd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_yz_xy = cbuffer.data(dd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_yz_xz = cbuffer.data(dd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_yz_yy = cbuffer.data(dd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_yz_yz = cbuffer.data(dd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_yz_zz = cbuffer.data(dd_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_zz_xx = cbuffer.data(dd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_zz_xy = cbuffer.data(dd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_zz_xz = cbuffer.data(dd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_zz_yy = cbuffer.data(dd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_zz_yz = cbuffer.data(dd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_zz_zz = cbuffer.data(dd_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_xx_xx = cbuffer.data(dd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_xx_xy = cbuffer.data(dd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_xx_xz = cbuffer.data(dd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_xx_yy = cbuffer.data(dd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_xx_yz = cbuffer.data(dd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_xx_zz = cbuffer.data(dd_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_xy_xx = cbuffer.data(dd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_xy_xy = cbuffer.data(dd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_xy_xz = cbuffer.data(dd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_xy_yy = cbuffer.data(dd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xy_yz = cbuffer.data(dd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xy_zz = cbuffer.data(dd_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_xz_xx = cbuffer.data(dd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_xz_xy = cbuffer.data(dd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_xz_xz = cbuffer.data(dd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_xz_yy = cbuffer.data(dd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_xz_yz = cbuffer.data(dd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_xz_zz = cbuffer.data(dd_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_yy_xx = cbuffer.data(dd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_yy_xy = cbuffer.data(dd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_yy_xz = cbuffer.data(dd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_yy_yy = cbuffer.data(dd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_yy_yz = cbuffer.data(dd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_yy_zz = cbuffer.data(dd_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_yz_xx = cbuffer.data(dd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_yz_xy = cbuffer.data(dd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_yz_xz = cbuffer.data(dd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_yz_yy = cbuffer.data(dd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_yz_yz = cbuffer.data(dd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_yz_zz = cbuffer.data(dd_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_zz_xx = cbuffer.data(dd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_zz_xy = cbuffer.data(dd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_zz_xz = cbuffer.data(dd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_zz_yy = cbuffer.data(dd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_zz_yz = cbuffer.data(dd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_zz_zz = cbuffer.data(dd_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_xx_xx = cbuffer.data(dd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_xx_xy = cbuffer.data(dd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_xx_xz = cbuffer.data(dd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_xx_yy = cbuffer.data(dd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_xx_yz = cbuffer.data(dd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_xx_zz = cbuffer.data(dd_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_xy_xx = cbuffer.data(dd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_xy_xy = cbuffer.data(dd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_xy_xz = cbuffer.data(dd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_xy_yy = cbuffer.data(dd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_xy_yz = cbuffer.data(dd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_xy_zz = cbuffer.data(dd_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_z_xz_xx = cbuffer.data(dd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_xz_xy = cbuffer.data(dd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_xz_xz = cbuffer.data(dd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_xz_yy = cbuffer.data(dd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_xz_yz = cbuffer.data(dd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_xz_zz = cbuffer.data(dd_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_z_yy_xx = cbuffer.data(dd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_yy_xy = cbuffer.data(dd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_yy_xz = cbuffer.data(dd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_yy_yy = cbuffer.data(dd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_yy_yz = cbuffer.data(dd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_yy_zz = cbuffer.data(dd_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_z_yz_xx = cbuffer.data(dd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_yz_xy = cbuffer.data(dd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_yz_xz = cbuffer.data(dd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_yz_yy = cbuffer.data(dd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_yz_yz = cbuffer.data(dd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_yz_zz = cbuffer.data(dd_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_z_zz_xx = cbuffer.data(dd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_zz_xy = cbuffer.data(dd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_zz_xz = cbuffer.data(dd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_zz_yy = cbuffer.data(dd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_zz_yz = cbuffer.data(dd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_zz_zz = cbuffer.data(dd_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_xx_xx = cbuffer.data(dd_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_xx_xy = cbuffer.data(dd_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_xx_xz = cbuffer.data(dd_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_xx_yy = cbuffer.data(dd_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_xx_yz = cbuffer.data(dd_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_xx_zz = cbuffer.data(dd_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_x_xy_xx = cbuffer.data(dd_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_xy_xy = cbuffer.data(dd_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_xy_xz = cbuffer.data(dd_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_xy_yy = cbuffer.data(dd_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_xy_yz = cbuffer.data(dd_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_xy_zz = cbuffer.data(dd_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_x_xz_xx = cbuffer.data(dd_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_x_xz_xy = cbuffer.data(dd_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_x_xz_xz = cbuffer.data(dd_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_x_xz_yy = cbuffer.data(dd_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_x_xz_yz = cbuffer.data(dd_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_x_xz_zz = cbuffer.data(dd_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_x_yy_xx = cbuffer.data(dd_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_x_yy_xy = cbuffer.data(dd_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_x_yy_xz = cbuffer.data(dd_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_x_yy_yy = cbuffer.data(dd_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_x_yy_yz = cbuffer.data(dd_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_x_yy_zz = cbuffer.data(dd_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_x_yz_xx = cbuffer.data(dd_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_x_yz_xy = cbuffer.data(dd_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_x_yz_xz = cbuffer.data(dd_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_yz_yy = cbuffer.data(dd_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_yz_yz = cbuffer.data(dd_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_yz_zz = cbuffer.data(dd_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_x_zz_xx = cbuffer.data(dd_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_zz_xy = cbuffer.data(dd_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_zz_xz = cbuffer.data(dd_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_zz_yy = cbuffer.data(dd_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_zz_yz = cbuffer.data(dd_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_zz_zz = cbuffer.data(dd_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_y_xx_xx = cbuffer.data(dd_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_xx_xy = cbuffer.data(dd_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_xx_xz = cbuffer.data(dd_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_xx_yy = cbuffer.data(dd_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_xx_yz = cbuffer.data(dd_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_xx_zz = cbuffer.data(dd_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_y_xy_xx = cbuffer.data(dd_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_y_xy_xy = cbuffer.data(dd_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_y_xy_xz = cbuffer.data(dd_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_y_xy_yy = cbuffer.data(dd_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_y_xy_yz = cbuffer.data(dd_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_y_xy_zz = cbuffer.data(dd_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_y_xz_xx = cbuffer.data(dd_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_y_xz_xy = cbuffer.data(dd_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_y_xz_xz = cbuffer.data(dd_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_y_xz_yy = cbuffer.data(dd_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_y_xz_yz = cbuffer.data(dd_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_y_xz_zz = cbuffer.data(dd_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_y_yy_xx = cbuffer.data(dd_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_y_yy_xy = cbuffer.data(dd_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_y_yy_xz = cbuffer.data(dd_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_y_yy_yy = cbuffer.data(dd_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_y_yy_yz = cbuffer.data(dd_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_y_yy_zz = cbuffer.data(dd_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_y_yz_xx = cbuffer.data(dd_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_y_yz_xy = cbuffer.data(dd_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_y_yz_xz = cbuffer.data(dd_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_y_yz_yy = cbuffer.data(dd_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_y_yz_yz = cbuffer.data(dd_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_y_yz_zz = cbuffer.data(dd_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_y_zz_xx = cbuffer.data(dd_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_y_zz_xy = cbuffer.data(dd_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_y_zz_xz = cbuffer.data(dd_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_y_zz_yy = cbuffer.data(dd_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_y_zz_yz = cbuffer.data(dd_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_y_zz_zz = cbuffer.data(dd_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_z_xx_xx = cbuffer.data(dd_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_z_xx_xy = cbuffer.data(dd_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_z_xx_xz = cbuffer.data(dd_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_z_xx_yy = cbuffer.data(dd_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_z_xx_yz = cbuffer.data(dd_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_z_xx_zz = cbuffer.data(dd_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_z_xy_xx = cbuffer.data(dd_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_z_xy_xy = cbuffer.data(dd_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_z_xy_xz = cbuffer.data(dd_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_z_xy_yy = cbuffer.data(dd_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_z_xy_yz = cbuffer.data(dd_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_z_xy_zz = cbuffer.data(dd_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_z_xz_xx = cbuffer.data(dd_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_z_xz_xy = cbuffer.data(dd_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_z_xz_xz = cbuffer.data(dd_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_z_xz_yy = cbuffer.data(dd_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_z_xz_yz = cbuffer.data(dd_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_z_xz_zz = cbuffer.data(dd_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_z_yy_xx = cbuffer.data(dd_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_z_yy_xy = cbuffer.data(dd_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_z_yy_xz = cbuffer.data(dd_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_z_yy_yy = cbuffer.data(dd_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_z_yy_yz = cbuffer.data(dd_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_z_yy_zz = cbuffer.data(dd_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_z_yz_xx = cbuffer.data(dd_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_z_yz_xy = cbuffer.data(dd_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_z_yz_xz = cbuffer.data(dd_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_z_yz_yy = cbuffer.data(dd_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_z_yz_yz = cbuffer.data(dd_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_z_yz_zz = cbuffer.data(dd_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_z_zz_xx = cbuffer.data(dd_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_z_zz_xy = cbuffer.data(dd_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_z_zz_xz = cbuffer.data(dd_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_z_zz_yy = cbuffer.data(dd_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_z_zz_yz = cbuffer.data(dd_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_z_zz_zz = cbuffer.data(dd_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_x_xx_xx = cbuffer.data(dd_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_x_xx_xy = cbuffer.data(dd_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_x_xx_xz = cbuffer.data(dd_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_x_xx_yy = cbuffer.data(dd_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_x_xx_yz = cbuffer.data(dd_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_x_xx_zz = cbuffer.data(dd_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_x_xy_xx = cbuffer.data(dd_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_x_xy_xy = cbuffer.data(dd_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_x_xy_xz = cbuffer.data(dd_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_x_xy_yy = cbuffer.data(dd_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_x_xy_yz = cbuffer.data(dd_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_x_xy_zz = cbuffer.data(dd_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_x_xz_xx = cbuffer.data(dd_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_x_xz_xy = cbuffer.data(dd_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_x_xz_xz = cbuffer.data(dd_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_x_xz_yy = cbuffer.data(dd_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_x_xz_yz = cbuffer.data(dd_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_x_xz_zz = cbuffer.data(dd_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_x_yy_xx = cbuffer.data(dd_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_x_yy_xy = cbuffer.data(dd_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_x_yy_xz = cbuffer.data(dd_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_x_yy_yy = cbuffer.data(dd_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_x_yy_yz = cbuffer.data(dd_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_x_yy_zz = cbuffer.data(dd_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_x_yz_xx = cbuffer.data(dd_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_x_yz_xy = cbuffer.data(dd_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_x_yz_xz = cbuffer.data(dd_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_x_yz_yy = cbuffer.data(dd_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_x_yz_yz = cbuffer.data(dd_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_x_yz_zz = cbuffer.data(dd_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_x_zz_xx = cbuffer.data(dd_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_x_zz_xy = cbuffer.data(dd_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_x_zz_xz = cbuffer.data(dd_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_x_zz_yy = cbuffer.data(dd_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_x_zz_yz = cbuffer.data(dd_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_x_zz_zz = cbuffer.data(dd_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_y_xx_xx = cbuffer.data(dd_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_y_xx_xy = cbuffer.data(dd_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_y_xx_xz = cbuffer.data(dd_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_y_xx_yy = cbuffer.data(dd_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_y_xx_yz = cbuffer.data(dd_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_y_xx_zz = cbuffer.data(dd_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_y_xy_xx = cbuffer.data(dd_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_y_xy_xy = cbuffer.data(dd_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_y_xy_xz = cbuffer.data(dd_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_y_xy_yy = cbuffer.data(dd_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_y_xy_yz = cbuffer.data(dd_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_y_xy_zz = cbuffer.data(dd_geom_11_off + 263 * ccomps * dcomps);

            auto g_z_y_xz_xx = cbuffer.data(dd_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_y_xz_xy = cbuffer.data(dd_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_y_xz_xz = cbuffer.data(dd_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_y_xz_yy = cbuffer.data(dd_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_y_xz_yz = cbuffer.data(dd_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_y_xz_zz = cbuffer.data(dd_geom_11_off + 269 * ccomps * dcomps);

            auto g_z_y_yy_xx = cbuffer.data(dd_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_y_yy_xy = cbuffer.data(dd_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_y_yy_xz = cbuffer.data(dd_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_y_yy_yy = cbuffer.data(dd_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_y_yy_yz = cbuffer.data(dd_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_y_yy_zz = cbuffer.data(dd_geom_11_off + 275 * ccomps * dcomps);

            auto g_z_y_yz_xx = cbuffer.data(dd_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_y_yz_xy = cbuffer.data(dd_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_y_yz_xz = cbuffer.data(dd_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_y_yz_yy = cbuffer.data(dd_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_y_yz_yz = cbuffer.data(dd_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_y_yz_zz = cbuffer.data(dd_geom_11_off + 281 * ccomps * dcomps);

            auto g_z_y_zz_xx = cbuffer.data(dd_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_y_zz_xy = cbuffer.data(dd_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_y_zz_xz = cbuffer.data(dd_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_y_zz_yy = cbuffer.data(dd_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_y_zz_yz = cbuffer.data(dd_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_y_zz_zz = cbuffer.data(dd_geom_11_off + 287 * ccomps * dcomps);

            auto g_z_z_xx_xx = cbuffer.data(dd_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_z_xx_xy = cbuffer.data(dd_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_z_xx_xz = cbuffer.data(dd_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_z_xx_yy = cbuffer.data(dd_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_z_xx_yz = cbuffer.data(dd_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_z_xx_zz = cbuffer.data(dd_geom_11_off + 293 * ccomps * dcomps);

            auto g_z_z_xy_xx = cbuffer.data(dd_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_z_xy_xy = cbuffer.data(dd_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_z_xy_xz = cbuffer.data(dd_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_z_xy_yy = cbuffer.data(dd_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_z_xy_yz = cbuffer.data(dd_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_z_xy_zz = cbuffer.data(dd_geom_11_off + 299 * ccomps * dcomps);

            auto g_z_z_xz_xx = cbuffer.data(dd_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_z_xz_xy = cbuffer.data(dd_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_z_xz_xz = cbuffer.data(dd_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_z_xz_yy = cbuffer.data(dd_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_z_xz_yz = cbuffer.data(dd_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_z_xz_zz = cbuffer.data(dd_geom_11_off + 305 * ccomps * dcomps);

            auto g_z_z_yy_xx = cbuffer.data(dd_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_z_yy_xy = cbuffer.data(dd_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_z_yy_xz = cbuffer.data(dd_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_z_yy_yy = cbuffer.data(dd_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_z_yy_yz = cbuffer.data(dd_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_z_yy_zz = cbuffer.data(dd_geom_11_off + 311 * ccomps * dcomps);

            auto g_z_z_yz_xx = cbuffer.data(dd_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_z_yz_xy = cbuffer.data(dd_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_z_yz_xz = cbuffer.data(dd_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_z_yz_yy = cbuffer.data(dd_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_z_yz_yz = cbuffer.data(dd_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_z_yz_zz = cbuffer.data(dd_geom_11_off + 317 * ccomps * dcomps);

            auto g_z_z_zz_xx = cbuffer.data(dd_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_z_zz_xy = cbuffer.data(dd_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_z_zz_xz = cbuffer.data(dd_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_z_zz_yy = cbuffer.data(dd_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_z_zz_yz = cbuffer.data(dd_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_z_zz_zz = cbuffer.data(dd_geom_11_off + 323 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fpxx

            const auto fp_geom_11_off = idx_geom_11_fpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxx_x = cbuffer.data(fp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxx_y = cbuffer.data(fp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxx_z = cbuffer.data(fp_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_x, g_0_x_xx_y, g_0_x_xx_z, g_x_0_xx_x, g_x_0_xx_y, g_x_0_xx_z, g_x_x_xx_x, g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_y, g_x_x_xx_z, g_x_x_xxx_x, g_x_x_xxx_y, g_x_x_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxx_x[k] = -g_0_x_xx_x[k] + g_x_0_xx_x[k] - g_x_x_xx_x[k] * ab_x + g_x_x_xx_xx[k];

                g_x_x_xxx_y[k] = -g_0_x_xx_y[k] + g_x_0_xx_y[k] - g_x_x_xx_y[k] * ab_x + g_x_x_xx_xy[k];

                g_x_x_xxx_z[k] = -g_0_x_xx_z[k] + g_x_0_xx_z[k] - g_x_x_xx_z[k] * ab_x + g_x_x_xx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxy_x = cbuffer.data(fp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxy_y = cbuffer.data(fp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxy_z = cbuffer.data(fp_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_x, g_x_x_xx_xy, g_x_x_xx_y, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_z, g_x_x_xxy_x, g_x_x_xxy_y, g_x_x_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxy_x[k] = -g_x_x_xx_x[k] * ab_y + g_x_x_xx_xy[k];

                g_x_x_xxy_y[k] = -g_x_x_xx_y[k] * ab_y + g_x_x_xx_yy[k];

                g_x_x_xxy_z[k] = -g_x_x_xx_z[k] * ab_y + g_x_x_xx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxz_x = cbuffer.data(fp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxz_y = cbuffer.data(fp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxz_z = cbuffer.data(fp_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_x, g_x_x_xx_xz, g_x_x_xx_y, g_x_x_xx_yz, g_x_x_xx_z, g_x_x_xx_zz, g_x_x_xxz_x, g_x_x_xxz_y, g_x_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxz_x[k] = -g_x_x_xx_x[k] * ab_z + g_x_x_xx_xz[k];

                g_x_x_xxz_y[k] = -g_x_x_xx_y[k] * ab_z + g_x_x_xx_yz[k];

                g_x_x_xxz_z[k] = -g_x_x_xx_z[k] * ab_z + g_x_x_xx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyy_x = cbuffer.data(fp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xyy_y = cbuffer.data(fp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xyy_z = cbuffer.data(fp_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xy_x, g_x_x_xy_xy, g_x_x_xy_y, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_z, g_x_x_xyy_x, g_x_x_xyy_y, g_x_x_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyy_x[k] = -g_x_x_xy_x[k] * ab_y + g_x_x_xy_xy[k];

                g_x_x_xyy_y[k] = -g_x_x_xy_y[k] * ab_y + g_x_x_xy_yy[k];

                g_x_x_xyy_z[k] = -g_x_x_xy_z[k] * ab_y + g_x_x_xy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyz_x = cbuffer.data(fp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xyz_y = cbuffer.data(fp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xyz_z = cbuffer.data(fp_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyz_x, g_x_x_xyz_y, g_x_x_xyz_z, g_x_x_xz_x, g_x_x_xz_xy, g_x_x_xz_y, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyz_x[k] = -g_x_x_xz_x[k] * ab_y + g_x_x_xz_xy[k];

                g_x_x_xyz_y[k] = -g_x_x_xz_y[k] * ab_y + g_x_x_xz_yy[k];

                g_x_x_xyz_z[k] = -g_x_x_xz_z[k] * ab_y + g_x_x_xz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzz_x = cbuffer.data(fp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xzz_y = cbuffer.data(fp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xzz_z = cbuffer.data(fp_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xz_x, g_x_x_xz_xz, g_x_x_xz_y, g_x_x_xz_yz, g_x_x_xz_z, g_x_x_xz_zz, g_x_x_xzz_x, g_x_x_xzz_y, g_x_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzz_x[k] = -g_x_x_xz_x[k] * ab_z + g_x_x_xz_xz[k];

                g_x_x_xzz_y[k] = -g_x_x_xz_y[k] * ab_z + g_x_x_xz_yz[k];

                g_x_x_xzz_z[k] = -g_x_x_xz_z[k] * ab_z + g_x_x_xz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyy_x = cbuffer.data(fp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_yyy_y = cbuffer.data(fp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_yyy_z = cbuffer.data(fp_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yy_x, g_x_x_yy_xy, g_x_x_yy_y, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_z, g_x_x_yyy_x, g_x_x_yyy_y, g_x_x_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyy_x[k] = -g_x_x_yy_x[k] * ab_y + g_x_x_yy_xy[k];

                g_x_x_yyy_y[k] = -g_x_x_yy_y[k] * ab_y + g_x_x_yy_yy[k];

                g_x_x_yyy_z[k] = -g_x_x_yy_z[k] * ab_y + g_x_x_yy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyz_x = cbuffer.data(fp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_yyz_y = cbuffer.data(fp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_yyz_z = cbuffer.data(fp_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyz_x, g_x_x_yyz_y, g_x_x_yyz_z, g_x_x_yz_x, g_x_x_yz_xy, g_x_x_yz_y, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyz_x[k] = -g_x_x_yz_x[k] * ab_y + g_x_x_yz_xy[k];

                g_x_x_yyz_y[k] = -g_x_x_yz_y[k] * ab_y + g_x_x_yz_yy[k];

                g_x_x_yyz_z[k] = -g_x_x_yz_z[k] * ab_y + g_x_x_yz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzz_x = cbuffer.data(fp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_yzz_y = cbuffer.data(fp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_yzz_z = cbuffer.data(fp_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzz_x, g_x_x_yzz_y, g_x_x_yzz_z, g_x_x_zz_x, g_x_x_zz_xy, g_x_x_zz_y, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzz_x[k] = -g_x_x_zz_x[k] * ab_y + g_x_x_zz_xy[k];

                g_x_x_yzz_y[k] = -g_x_x_zz_y[k] * ab_y + g_x_x_zz_yy[k];

                g_x_x_yzz_z[k] = -g_x_x_zz_z[k] * ab_y + g_x_x_zz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzz_x = cbuffer.data(fp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_zzz_y = cbuffer.data(fp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_zzz_z = cbuffer.data(fp_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zz_x, g_x_x_zz_xz, g_x_x_zz_y, g_x_x_zz_yz, g_x_x_zz_z, g_x_x_zz_zz, g_x_x_zzz_x, g_x_x_zzz_y, g_x_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzz_x[k] = -g_x_x_zz_x[k] * ab_z + g_x_x_zz_xz[k];

                g_x_x_zzz_y[k] = -g_x_x_zz_y[k] * ab_z + g_x_x_zz_yz[k];

                g_x_x_zzz_z[k] = -g_x_x_zz_z[k] * ab_z + g_x_x_zz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxx_x = cbuffer.data(fp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_xxx_y = cbuffer.data(fp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_xxx_z = cbuffer.data(fp_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_x, g_0_y_xx_y, g_0_y_xx_z, g_x_y_xx_x, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_y, g_x_y_xx_z, g_x_y_xxx_x, g_x_y_xxx_y, g_x_y_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxx_x[k] = -g_0_y_xx_x[k] - g_x_y_xx_x[k] * ab_x + g_x_y_xx_xx[k];

                g_x_y_xxx_y[k] = -g_0_y_xx_y[k] - g_x_y_xx_y[k] * ab_x + g_x_y_xx_xy[k];

                g_x_y_xxx_z[k] = -g_0_y_xx_z[k] - g_x_y_xx_z[k] * ab_x + g_x_y_xx_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxy_x = cbuffer.data(fp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_xxy_y = cbuffer.data(fp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_xxy_z = cbuffer.data(fp_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_x_y_xxy_x, g_x_y_xxy_y, g_x_y_xxy_z, g_x_y_xy_x, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_y, g_x_y_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxy_x[k] = -g_0_y_xy_x[k] - g_x_y_xy_x[k] * ab_x + g_x_y_xy_xx[k];

                g_x_y_xxy_y[k] = -g_0_y_xy_y[k] - g_x_y_xy_y[k] * ab_x + g_x_y_xy_xy[k];

                g_x_y_xxy_z[k] = -g_0_y_xy_z[k] - g_x_y_xy_z[k] * ab_x + g_x_y_xy_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxz_x = cbuffer.data(fp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_xxz_y = cbuffer.data(fp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_xxz_z = cbuffer.data(fp_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xx_x, g_x_y_xx_xz, g_x_y_xx_y, g_x_y_xx_yz, g_x_y_xx_z, g_x_y_xx_zz, g_x_y_xxz_x, g_x_y_xxz_y, g_x_y_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxz_x[k] = -g_x_y_xx_x[k] * ab_z + g_x_y_xx_xz[k];

                g_x_y_xxz_y[k] = -g_x_y_xx_y[k] * ab_z + g_x_y_xx_yz[k];

                g_x_y_xxz_z[k] = -g_x_y_xx_z[k] * ab_z + g_x_y_xx_zz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyy_x = cbuffer.data(fp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_xyy_y = cbuffer.data(fp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_xyy_z = cbuffer.data(fp_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_x_y_xyy_x, g_x_y_xyy_y, g_x_y_xyy_z, g_x_y_yy_x, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_y, g_x_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyy_x[k] = -g_0_y_yy_x[k] - g_x_y_yy_x[k] * ab_x + g_x_y_yy_xx[k];

                g_x_y_xyy_y[k] = -g_0_y_yy_y[k] - g_x_y_yy_y[k] * ab_x + g_x_y_yy_xy[k];

                g_x_y_xyy_z[k] = -g_0_y_yy_z[k] - g_x_y_yy_z[k] * ab_x + g_x_y_yy_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyz_x = cbuffer.data(fp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_xyz_y = cbuffer.data(fp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_xyz_z = cbuffer.data(fp_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xy_x, g_x_y_xy_xz, g_x_y_xy_y, g_x_y_xy_yz, g_x_y_xy_z, g_x_y_xy_zz, g_x_y_xyz_x, g_x_y_xyz_y, g_x_y_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyz_x[k] = -g_x_y_xy_x[k] * ab_z + g_x_y_xy_xz[k];

                g_x_y_xyz_y[k] = -g_x_y_xy_y[k] * ab_z + g_x_y_xy_yz[k];

                g_x_y_xyz_z[k] = -g_x_y_xy_z[k] * ab_z + g_x_y_xy_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzz_x = cbuffer.data(fp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xzz_y = cbuffer.data(fp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xzz_z = cbuffer.data(fp_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xz_x, g_x_y_xz_xz, g_x_y_xz_y, g_x_y_xz_yz, g_x_y_xz_z, g_x_y_xz_zz, g_x_y_xzz_x, g_x_y_xzz_y, g_x_y_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzz_x[k] = -g_x_y_xz_x[k] * ab_z + g_x_y_xz_xz[k];

                g_x_y_xzz_y[k] = -g_x_y_xz_y[k] * ab_z + g_x_y_xz_yz[k];

                g_x_y_xzz_z[k] = -g_x_y_xz_z[k] * ab_z + g_x_y_xz_zz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyy_x = cbuffer.data(fp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_yyy_y = cbuffer.data(fp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_yyy_z = cbuffer.data(fp_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_x, g_x_0_yy_y, g_x_0_yy_z, g_x_y_yy_x, g_x_y_yy_xy, g_x_y_yy_y, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_z, g_x_y_yyy_x, g_x_y_yyy_y, g_x_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyy_x[k] = g_x_0_yy_x[k] - g_x_y_yy_x[k] * ab_y + g_x_y_yy_xy[k];

                g_x_y_yyy_y[k] = g_x_0_yy_y[k] - g_x_y_yy_y[k] * ab_y + g_x_y_yy_yy[k];

                g_x_y_yyy_z[k] = g_x_0_yy_z[k] - g_x_y_yy_z[k] * ab_y + g_x_y_yy_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyz_x = cbuffer.data(fp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_yyz_y = cbuffer.data(fp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_yyz_z = cbuffer.data(fp_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yy_x, g_x_y_yy_xz, g_x_y_yy_y, g_x_y_yy_yz, g_x_y_yy_z, g_x_y_yy_zz, g_x_y_yyz_x, g_x_y_yyz_y, g_x_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyz_x[k] = -g_x_y_yy_x[k] * ab_z + g_x_y_yy_xz[k];

                g_x_y_yyz_y[k] = -g_x_y_yy_y[k] * ab_z + g_x_y_yy_yz[k];

                g_x_y_yyz_z[k] = -g_x_y_yy_z[k] * ab_z + g_x_y_yy_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzz_x = cbuffer.data(fp_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_yzz_y = cbuffer.data(fp_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_yzz_z = cbuffer.data(fp_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yz_x, g_x_y_yz_xz, g_x_y_yz_y, g_x_y_yz_yz, g_x_y_yz_z, g_x_y_yz_zz, g_x_y_yzz_x, g_x_y_yzz_y, g_x_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzz_x[k] = -g_x_y_yz_x[k] * ab_z + g_x_y_yz_xz[k];

                g_x_y_yzz_y[k] = -g_x_y_yz_y[k] * ab_z + g_x_y_yz_yz[k];

                g_x_y_yzz_z[k] = -g_x_y_yz_z[k] * ab_z + g_x_y_yz_zz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzz_x = cbuffer.data(fp_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_zzz_y = cbuffer.data(fp_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_zzz_z = cbuffer.data(fp_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zz_x, g_x_y_zz_xz, g_x_y_zz_y, g_x_y_zz_yz, g_x_y_zz_z, g_x_y_zz_zz, g_x_y_zzz_x, g_x_y_zzz_y, g_x_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzz_x[k] = -g_x_y_zz_x[k] * ab_z + g_x_y_zz_xz[k];

                g_x_y_zzz_y[k] = -g_x_y_zz_y[k] * ab_z + g_x_y_zz_yz[k];

                g_x_y_zzz_z[k] = -g_x_y_zz_z[k] * ab_z + g_x_y_zz_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxx_x = cbuffer.data(fp_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_xxx_y = cbuffer.data(fp_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_xxx_z = cbuffer.data(fp_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_x, g_0_z_xx_y, g_0_z_xx_z, g_x_z_xx_x, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_y, g_x_z_xx_z, g_x_z_xxx_x, g_x_z_xxx_y, g_x_z_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxx_x[k] = -g_0_z_xx_x[k] - g_x_z_xx_x[k] * ab_x + g_x_z_xx_xx[k];

                g_x_z_xxx_y[k] = -g_0_z_xx_y[k] - g_x_z_xx_y[k] * ab_x + g_x_z_xx_xy[k];

                g_x_z_xxx_z[k] = -g_0_z_xx_z[k] - g_x_z_xx_z[k] * ab_x + g_x_z_xx_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxy_x = cbuffer.data(fp_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_xxy_y = cbuffer.data(fp_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_xxy_z = cbuffer.data(fp_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xx_x, g_x_z_xx_xy, g_x_z_xx_y, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_z, g_x_z_xxy_x, g_x_z_xxy_y, g_x_z_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxy_x[k] = -g_x_z_xx_x[k] * ab_y + g_x_z_xx_xy[k];

                g_x_z_xxy_y[k] = -g_x_z_xx_y[k] * ab_y + g_x_z_xx_yy[k];

                g_x_z_xxy_z[k] = -g_x_z_xx_z[k] * ab_y + g_x_z_xx_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxz_x = cbuffer.data(fp_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_xxz_y = cbuffer.data(fp_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_xxz_z = cbuffer.data(fp_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_x_z_xxz_x, g_x_z_xxz_y, g_x_z_xxz_z, g_x_z_xz_x, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_y, g_x_z_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxz_x[k] = -g_0_z_xz_x[k] - g_x_z_xz_x[k] * ab_x + g_x_z_xz_xx[k];

                g_x_z_xxz_y[k] = -g_0_z_xz_y[k] - g_x_z_xz_y[k] * ab_x + g_x_z_xz_xy[k];

                g_x_z_xxz_z[k] = -g_0_z_xz_z[k] - g_x_z_xz_z[k] * ab_x + g_x_z_xz_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyy_x = cbuffer.data(fp_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_z_xyy_y = cbuffer.data(fp_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_xyy_z = cbuffer.data(fp_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xy_x, g_x_z_xy_xy, g_x_z_xy_y, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_z, g_x_z_xyy_x, g_x_z_xyy_y, g_x_z_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyy_x[k] = -g_x_z_xy_x[k] * ab_y + g_x_z_xy_xy[k];

                g_x_z_xyy_y[k] = -g_x_z_xy_y[k] * ab_y + g_x_z_xy_yy[k];

                g_x_z_xyy_z[k] = -g_x_z_xy_z[k] * ab_y + g_x_z_xy_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyz_x = cbuffer.data(fp_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_xyz_y = cbuffer.data(fp_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_xyz_z = cbuffer.data(fp_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyz_x, g_x_z_xyz_y, g_x_z_xyz_z, g_x_z_xz_x, g_x_z_xz_xy, g_x_z_xz_y, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyz_x[k] = -g_x_z_xz_x[k] * ab_y + g_x_z_xz_xy[k];

                g_x_z_xyz_y[k] = -g_x_z_xz_y[k] * ab_y + g_x_z_xz_yy[k];

                g_x_z_xyz_z[k] = -g_x_z_xz_z[k] * ab_y + g_x_z_xz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzz_x = cbuffer.data(fp_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_xzz_y = cbuffer.data(fp_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_xzz_z = cbuffer.data(fp_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_x_z_xzz_x, g_x_z_xzz_y, g_x_z_xzz_z, g_x_z_zz_x, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_y, g_x_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzz_x[k] = -g_0_z_zz_x[k] - g_x_z_zz_x[k] * ab_x + g_x_z_zz_xx[k];

                g_x_z_xzz_y[k] = -g_0_z_zz_y[k] - g_x_z_zz_y[k] * ab_x + g_x_z_zz_xy[k];

                g_x_z_xzz_z[k] = -g_0_z_zz_z[k] - g_x_z_zz_z[k] * ab_x + g_x_z_zz_xz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyy_x = cbuffer.data(fp_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_yyy_y = cbuffer.data(fp_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_yyy_z = cbuffer.data(fp_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yy_x, g_x_z_yy_xy, g_x_z_yy_y, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_z, g_x_z_yyy_x, g_x_z_yyy_y, g_x_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyy_x[k] = -g_x_z_yy_x[k] * ab_y + g_x_z_yy_xy[k];

                g_x_z_yyy_y[k] = -g_x_z_yy_y[k] * ab_y + g_x_z_yy_yy[k];

                g_x_z_yyy_z[k] = -g_x_z_yy_z[k] * ab_y + g_x_z_yy_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyz_x = cbuffer.data(fp_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_yyz_y = cbuffer.data(fp_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_yyz_z = cbuffer.data(fp_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyz_x, g_x_z_yyz_y, g_x_z_yyz_z, g_x_z_yz_x, g_x_z_yz_xy, g_x_z_yz_y, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyz_x[k] = -g_x_z_yz_x[k] * ab_y + g_x_z_yz_xy[k];

                g_x_z_yyz_y[k] = -g_x_z_yz_y[k] * ab_y + g_x_z_yz_yy[k];

                g_x_z_yyz_z[k] = -g_x_z_yz_z[k] * ab_y + g_x_z_yz_yz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzz_x = cbuffer.data(fp_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_yzz_y = cbuffer.data(fp_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_yzz_z = cbuffer.data(fp_geom_11_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzz_x, g_x_z_yzz_y, g_x_z_yzz_z, g_x_z_zz_x, g_x_z_zz_xy, g_x_z_zz_y, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzz_x[k] = -g_x_z_zz_x[k] * ab_y + g_x_z_zz_xy[k];

                g_x_z_yzz_y[k] = -g_x_z_zz_y[k] * ab_y + g_x_z_zz_yy[k];

                g_x_z_yzz_z[k] = -g_x_z_zz_z[k] * ab_y + g_x_z_zz_yz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzz_x = cbuffer.data(fp_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_zzz_y = cbuffer.data(fp_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_zzz_z = cbuffer.data(fp_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_x, g_x_0_zz_y, g_x_0_zz_z, g_x_z_zz_x, g_x_z_zz_xz, g_x_z_zz_y, g_x_z_zz_yz, g_x_z_zz_z, g_x_z_zz_zz, g_x_z_zzz_x, g_x_z_zzz_y, g_x_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzz_x[k] = g_x_0_zz_x[k] - g_x_z_zz_x[k] * ab_z + g_x_z_zz_xz[k];

                g_x_z_zzz_y[k] = g_x_0_zz_y[k] - g_x_z_zz_y[k] * ab_z + g_x_z_zz_yz[k];

                g_x_z_zzz_z[k] = g_x_0_zz_z[k] - g_x_z_zz_z[k] * ab_z + g_x_z_zz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxx_x = cbuffer.data(fp_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_xxx_y = cbuffer.data(fp_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_xxx_z = cbuffer.data(fp_geom_11_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_x, g_y_0_xx_y, g_y_0_xx_z, g_y_x_xx_x, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_y, g_y_x_xx_z, g_y_x_xxx_x, g_y_x_xxx_y, g_y_x_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxx_x[k] = g_y_0_xx_x[k] - g_y_x_xx_x[k] * ab_x + g_y_x_xx_xx[k];

                g_y_x_xxx_y[k] = g_y_0_xx_y[k] - g_y_x_xx_y[k] * ab_x + g_y_x_xx_xy[k];

                g_y_x_xxx_z[k] = g_y_0_xx_z[k] - g_y_x_xx_z[k] * ab_x + g_y_x_xx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxy_x = cbuffer.data(fp_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_xxy_y = cbuffer.data(fp_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_xxy_z = cbuffer.data(fp_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_x, g_y_0_xy_y, g_y_0_xy_z, g_y_x_xxy_x, g_y_x_xxy_y, g_y_x_xxy_z, g_y_x_xy_x, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_y, g_y_x_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxy_x[k] = g_y_0_xy_x[k] - g_y_x_xy_x[k] * ab_x + g_y_x_xy_xx[k];

                g_y_x_xxy_y[k] = g_y_0_xy_y[k] - g_y_x_xy_y[k] * ab_x + g_y_x_xy_xy[k];

                g_y_x_xxy_z[k] = g_y_0_xy_z[k] - g_y_x_xy_z[k] * ab_x + g_y_x_xy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxz_x = cbuffer.data(fp_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_xxz_y = cbuffer.data(fp_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_xxz_z = cbuffer.data(fp_geom_11_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xx_x, g_y_x_xx_xz, g_y_x_xx_y, g_y_x_xx_yz, g_y_x_xx_z, g_y_x_xx_zz, g_y_x_xxz_x, g_y_x_xxz_y, g_y_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxz_x[k] = -g_y_x_xx_x[k] * ab_z + g_y_x_xx_xz[k];

                g_y_x_xxz_y[k] = -g_y_x_xx_y[k] * ab_z + g_y_x_xx_yz[k];

                g_y_x_xxz_z[k] = -g_y_x_xx_z[k] * ab_z + g_y_x_xx_zz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyy_x = cbuffer.data(fp_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_x_xyy_y = cbuffer.data(fp_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_xyy_z = cbuffer.data(fp_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z, g_y_x_xyy_x, g_y_x_xyy_y, g_y_x_xyy_z, g_y_x_yy_x, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_y, g_y_x_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyy_x[k] = g_y_0_yy_x[k] - g_y_x_yy_x[k] * ab_x + g_y_x_yy_xx[k];

                g_y_x_xyy_y[k] = g_y_0_yy_y[k] - g_y_x_yy_y[k] * ab_x + g_y_x_yy_xy[k];

                g_y_x_xyy_z[k] = g_y_0_yy_z[k] - g_y_x_yy_z[k] * ab_x + g_y_x_yy_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyz_x = cbuffer.data(fp_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_xyz_y = cbuffer.data(fp_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_xyz_z = cbuffer.data(fp_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xy_x, g_y_x_xy_xz, g_y_x_xy_y, g_y_x_xy_yz, g_y_x_xy_z, g_y_x_xy_zz, g_y_x_xyz_x, g_y_x_xyz_y, g_y_x_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyz_x[k] = -g_y_x_xy_x[k] * ab_z + g_y_x_xy_xz[k];

                g_y_x_xyz_y[k] = -g_y_x_xy_y[k] * ab_z + g_y_x_xy_yz[k];

                g_y_x_xyz_z[k] = -g_y_x_xy_z[k] * ab_z + g_y_x_xy_zz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzz_x = cbuffer.data(fp_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_xzz_y = cbuffer.data(fp_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_xzz_z = cbuffer.data(fp_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xz_x, g_y_x_xz_xz, g_y_x_xz_y, g_y_x_xz_yz, g_y_x_xz_z, g_y_x_xz_zz, g_y_x_xzz_x, g_y_x_xzz_y, g_y_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzz_x[k] = -g_y_x_xz_x[k] * ab_z + g_y_x_xz_xz[k];

                g_y_x_xzz_y[k] = -g_y_x_xz_y[k] * ab_z + g_y_x_xz_yz[k];

                g_y_x_xzz_z[k] = -g_y_x_xz_z[k] * ab_z + g_y_x_xz_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyy_x = cbuffer.data(fp_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_yyy_y = cbuffer.data(fp_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_yyy_z = cbuffer.data(fp_geom_11_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_x, g_0_x_yy_y, g_0_x_yy_z, g_y_x_yy_x, g_y_x_yy_xy, g_y_x_yy_y, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_z, g_y_x_yyy_x, g_y_x_yyy_y, g_y_x_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyy_x[k] = -g_0_x_yy_x[k] - g_y_x_yy_x[k] * ab_y + g_y_x_yy_xy[k];

                g_y_x_yyy_y[k] = -g_0_x_yy_y[k] - g_y_x_yy_y[k] * ab_y + g_y_x_yy_yy[k];

                g_y_x_yyy_z[k] = -g_0_x_yy_z[k] - g_y_x_yy_z[k] * ab_y + g_y_x_yy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyz_x = cbuffer.data(fp_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_yyz_y = cbuffer.data(fp_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_yyz_z = cbuffer.data(fp_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yy_x, g_y_x_yy_xz, g_y_x_yy_y, g_y_x_yy_yz, g_y_x_yy_z, g_y_x_yy_zz, g_y_x_yyz_x, g_y_x_yyz_y, g_y_x_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyz_x[k] = -g_y_x_yy_x[k] * ab_z + g_y_x_yy_xz[k];

                g_y_x_yyz_y[k] = -g_y_x_yy_y[k] * ab_z + g_y_x_yy_yz[k];

                g_y_x_yyz_z[k] = -g_y_x_yy_z[k] * ab_z + g_y_x_yy_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzz_x = cbuffer.data(fp_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_yzz_y = cbuffer.data(fp_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_yzz_z = cbuffer.data(fp_geom_11_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yz_x, g_y_x_yz_xz, g_y_x_yz_y, g_y_x_yz_yz, g_y_x_yz_z, g_y_x_yz_zz, g_y_x_yzz_x, g_y_x_yzz_y, g_y_x_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzz_x[k] = -g_y_x_yz_x[k] * ab_z + g_y_x_yz_xz[k];

                g_y_x_yzz_y[k] = -g_y_x_yz_y[k] * ab_z + g_y_x_yz_yz[k];

                g_y_x_yzz_z[k] = -g_y_x_yz_z[k] * ab_z + g_y_x_yz_zz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzz_x = cbuffer.data(fp_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_zzz_y = cbuffer.data(fp_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_zzz_z = cbuffer.data(fp_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zz_x, g_y_x_zz_xz, g_y_x_zz_y, g_y_x_zz_yz, g_y_x_zz_z, g_y_x_zz_zz, g_y_x_zzz_x, g_y_x_zzz_y, g_y_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzz_x[k] = -g_y_x_zz_x[k] * ab_z + g_y_x_zz_xz[k];

                g_y_x_zzz_y[k] = -g_y_x_zz_y[k] * ab_z + g_y_x_zz_yz[k];

                g_y_x_zzz_z[k] = -g_y_x_zz_z[k] * ab_z + g_y_x_zz_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxx_x = cbuffer.data(fp_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_xxx_y = cbuffer.data(fp_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_xxx_z = cbuffer.data(fp_geom_11_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xx_x, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_y, g_y_y_xx_z, g_y_y_xxx_x, g_y_y_xxx_y, g_y_y_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxx_x[k] = -g_y_y_xx_x[k] * ab_x + g_y_y_xx_xx[k];

                g_y_y_xxx_y[k] = -g_y_y_xx_y[k] * ab_x + g_y_y_xx_xy[k];

                g_y_y_xxx_z[k] = -g_y_y_xx_z[k] * ab_x + g_y_y_xx_xz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxy_x = cbuffer.data(fp_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_xxy_y = cbuffer.data(fp_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_xxy_z = cbuffer.data(fp_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxy_x, g_y_y_xxy_y, g_y_y_xxy_z, g_y_y_xy_x, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_y, g_y_y_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxy_x[k] = -g_y_y_xy_x[k] * ab_x + g_y_y_xy_xx[k];

                g_y_y_xxy_y[k] = -g_y_y_xy_y[k] * ab_x + g_y_y_xy_xy[k];

                g_y_y_xxy_z[k] = -g_y_y_xy_z[k] * ab_x + g_y_y_xy_xz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxz_x = cbuffer.data(fp_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_xxz_y = cbuffer.data(fp_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_xxz_z = cbuffer.data(fp_geom_11_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxz_x, g_y_y_xxz_y, g_y_y_xxz_z, g_y_y_xz_x, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_y, g_y_y_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxz_x[k] = -g_y_y_xz_x[k] * ab_x + g_y_y_xz_xx[k];

                g_y_y_xxz_y[k] = -g_y_y_xz_y[k] * ab_x + g_y_y_xz_xy[k];

                g_y_y_xxz_z[k] = -g_y_y_xz_z[k] * ab_x + g_y_y_xz_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyy_x = cbuffer.data(fp_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_y_xyy_y = cbuffer.data(fp_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_xyy_z = cbuffer.data(fp_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyy_x, g_y_y_xyy_y, g_y_y_xyy_z, g_y_y_yy_x, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_y, g_y_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyy_x[k] = -g_y_y_yy_x[k] * ab_x + g_y_y_yy_xx[k];

                g_y_y_xyy_y[k] = -g_y_y_yy_y[k] * ab_x + g_y_y_yy_xy[k];

                g_y_y_xyy_z[k] = -g_y_y_yy_z[k] * ab_x + g_y_y_yy_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyz_x = cbuffer.data(fp_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_xyz_y = cbuffer.data(fp_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_xyz_z = cbuffer.data(fp_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyz_x, g_y_y_xyz_y, g_y_y_xyz_z, g_y_y_yz_x, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_y, g_y_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyz_x[k] = -g_y_y_yz_x[k] * ab_x + g_y_y_yz_xx[k];

                g_y_y_xyz_y[k] = -g_y_y_yz_y[k] * ab_x + g_y_y_yz_xy[k];

                g_y_y_xyz_z[k] = -g_y_y_yz_z[k] * ab_x + g_y_y_yz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzz_x = cbuffer.data(fp_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_xzz_y = cbuffer.data(fp_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_xzz_z = cbuffer.data(fp_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzz_x, g_y_y_xzz_y, g_y_y_xzz_z, g_y_y_zz_x, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_y, g_y_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzz_x[k] = -g_y_y_zz_x[k] * ab_x + g_y_y_zz_xx[k];

                g_y_y_xzz_y[k] = -g_y_y_zz_y[k] * ab_x + g_y_y_zz_xy[k];

                g_y_y_xzz_z[k] = -g_y_y_zz_z[k] * ab_x + g_y_y_zz_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyy_x = cbuffer.data(fp_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_yyy_y = cbuffer.data(fp_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_y_yyy_z = cbuffer.data(fp_geom_11_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z, g_y_y_yy_x, g_y_y_yy_xy, g_y_y_yy_y, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_z, g_y_y_yyy_x, g_y_y_yyy_y, g_y_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyy_x[k] = -g_0_y_yy_x[k] + g_y_0_yy_x[k] - g_y_y_yy_x[k] * ab_y + g_y_y_yy_xy[k];

                g_y_y_yyy_y[k] = -g_0_y_yy_y[k] + g_y_0_yy_y[k] - g_y_y_yy_y[k] * ab_y + g_y_y_yy_yy[k];

                g_y_y_yyy_z[k] = -g_0_y_yy_z[k] + g_y_0_yy_z[k] - g_y_y_yy_z[k] * ab_y + g_y_y_yy_yz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyz_x = cbuffer.data(fp_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_y_yyz_y = cbuffer.data(fp_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_y_yyz_z = cbuffer.data(fp_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yy_x, g_y_y_yy_xz, g_y_y_yy_y, g_y_y_yy_yz, g_y_y_yy_z, g_y_y_yy_zz, g_y_y_yyz_x, g_y_y_yyz_y, g_y_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyz_x[k] = -g_y_y_yy_x[k] * ab_z + g_y_y_yy_xz[k];

                g_y_y_yyz_y[k] = -g_y_y_yy_y[k] * ab_z + g_y_y_yy_yz[k];

                g_y_y_yyz_z[k] = -g_y_y_yy_z[k] * ab_z + g_y_y_yy_zz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzz_x = cbuffer.data(fp_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_yzz_y = cbuffer.data(fp_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_yzz_z = cbuffer.data(fp_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yz_x, g_y_y_yz_xz, g_y_y_yz_y, g_y_y_yz_yz, g_y_y_yz_z, g_y_y_yz_zz, g_y_y_yzz_x, g_y_y_yzz_y, g_y_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzz_x[k] = -g_y_y_yz_x[k] * ab_z + g_y_y_yz_xz[k];

                g_y_y_yzz_y[k] = -g_y_y_yz_y[k] * ab_z + g_y_y_yz_yz[k];

                g_y_y_yzz_z[k] = -g_y_y_yz_z[k] * ab_z + g_y_y_yz_zz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzz_x = cbuffer.data(fp_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_zzz_y = cbuffer.data(fp_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_zzz_z = cbuffer.data(fp_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zz_x, g_y_y_zz_xz, g_y_y_zz_y, g_y_y_zz_yz, g_y_y_zz_z, g_y_y_zz_zz, g_y_y_zzz_x, g_y_y_zzz_y, g_y_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzz_x[k] = -g_y_y_zz_x[k] * ab_z + g_y_y_zz_xz[k];

                g_y_y_zzz_y[k] = -g_y_y_zz_y[k] * ab_z + g_y_y_zz_yz[k];

                g_y_y_zzz_z[k] = -g_y_y_zz_z[k] * ab_z + g_y_y_zz_zz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxx_x = cbuffer.data(fp_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_xxx_y = cbuffer.data(fp_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_xxx_z = cbuffer.data(fp_geom_11_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xx_x, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_y, g_y_z_xx_z, g_y_z_xxx_x, g_y_z_xxx_y, g_y_z_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxx_x[k] = -g_y_z_xx_x[k] * ab_x + g_y_z_xx_xx[k];

                g_y_z_xxx_y[k] = -g_y_z_xx_y[k] * ab_x + g_y_z_xx_xy[k];

                g_y_z_xxx_z[k] = -g_y_z_xx_z[k] * ab_x + g_y_z_xx_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxy_x = cbuffer.data(fp_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_xxy_y = cbuffer.data(fp_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_xxy_z = cbuffer.data(fp_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxy_x, g_y_z_xxy_y, g_y_z_xxy_z, g_y_z_xy_x, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_y, g_y_z_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxy_x[k] = -g_y_z_xy_x[k] * ab_x + g_y_z_xy_xx[k];

                g_y_z_xxy_y[k] = -g_y_z_xy_y[k] * ab_x + g_y_z_xy_xy[k];

                g_y_z_xxy_z[k] = -g_y_z_xy_z[k] * ab_x + g_y_z_xy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxz_x = cbuffer.data(fp_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_xxz_y = cbuffer.data(fp_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_xxz_z = cbuffer.data(fp_geom_11_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxz_x, g_y_z_xxz_y, g_y_z_xxz_z, g_y_z_xz_x, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_y, g_y_z_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxz_x[k] = -g_y_z_xz_x[k] * ab_x + g_y_z_xz_xx[k];

                g_y_z_xxz_y[k] = -g_y_z_xz_y[k] * ab_x + g_y_z_xz_xy[k];

                g_y_z_xxz_z[k] = -g_y_z_xz_z[k] * ab_x + g_y_z_xz_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyy_x = cbuffer.data(fp_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_z_xyy_y = cbuffer.data(fp_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_xyy_z = cbuffer.data(fp_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyy_x, g_y_z_xyy_y, g_y_z_xyy_z, g_y_z_yy_x, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_y, g_y_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyy_x[k] = -g_y_z_yy_x[k] * ab_x + g_y_z_yy_xx[k];

                g_y_z_xyy_y[k] = -g_y_z_yy_y[k] * ab_x + g_y_z_yy_xy[k];

                g_y_z_xyy_z[k] = -g_y_z_yy_z[k] * ab_x + g_y_z_yy_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyz_x = cbuffer.data(fp_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_xyz_y = cbuffer.data(fp_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_xyz_z = cbuffer.data(fp_geom_11_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyz_x, g_y_z_xyz_y, g_y_z_xyz_z, g_y_z_yz_x, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_y, g_y_z_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyz_x[k] = -g_y_z_yz_x[k] * ab_x + g_y_z_yz_xx[k];

                g_y_z_xyz_y[k] = -g_y_z_yz_y[k] * ab_x + g_y_z_yz_xy[k];

                g_y_z_xyz_z[k] = -g_y_z_yz_z[k] * ab_x + g_y_z_yz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzz_x = cbuffer.data(fp_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_xzz_y = cbuffer.data(fp_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_xzz_z = cbuffer.data(fp_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzz_x, g_y_z_xzz_y, g_y_z_xzz_z, g_y_z_zz_x, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_y, g_y_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzz_x[k] = -g_y_z_zz_x[k] * ab_x + g_y_z_zz_xx[k];

                g_y_z_xzz_y[k] = -g_y_z_zz_y[k] * ab_x + g_y_z_zz_xy[k];

                g_y_z_xzz_z[k] = -g_y_z_zz_z[k] * ab_x + g_y_z_zz_xz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyy_x = cbuffer.data(fp_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_z_yyy_y = cbuffer.data(fp_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_z_yyy_z = cbuffer.data(fp_geom_11_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_x, g_0_z_yy_y, g_0_z_yy_z, g_y_z_yy_x, g_y_z_yy_xy, g_y_z_yy_y, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_z, g_y_z_yyy_x, g_y_z_yyy_y, g_y_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyy_x[k] = -g_0_z_yy_x[k] - g_y_z_yy_x[k] * ab_y + g_y_z_yy_xy[k];

                g_y_z_yyy_y[k] = -g_0_z_yy_y[k] - g_y_z_yy_y[k] * ab_y + g_y_z_yy_yy[k];

                g_y_z_yyy_z[k] = -g_0_z_yy_z[k] - g_y_z_yy_z[k] * ab_y + g_y_z_yy_yz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyz_x = cbuffer.data(fp_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_z_yyz_y = cbuffer.data(fp_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_z_yyz_z = cbuffer.data(fp_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_y_z_yyz_x, g_y_z_yyz_y, g_y_z_yyz_z, g_y_z_yz_x, g_y_z_yz_xy, g_y_z_yz_y, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyz_x[k] = -g_0_z_yz_x[k] - g_y_z_yz_x[k] * ab_y + g_y_z_yz_xy[k];

                g_y_z_yyz_y[k] = -g_0_z_yz_y[k] - g_y_z_yz_y[k] * ab_y + g_y_z_yz_yy[k];

                g_y_z_yyz_z[k] = -g_0_z_yz_z[k] - g_y_z_yz_z[k] * ab_y + g_y_z_yz_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzz_x = cbuffer.data(fp_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_z_yzz_y = cbuffer.data(fp_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_z_yzz_z = cbuffer.data(fp_geom_11_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_y_z_yzz_x, g_y_z_yzz_y, g_y_z_yzz_z, g_y_z_zz_x, g_y_z_zz_xy, g_y_z_zz_y, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzz_x[k] = -g_0_z_zz_x[k] - g_y_z_zz_x[k] * ab_y + g_y_z_zz_xy[k];

                g_y_z_yzz_y[k] = -g_0_z_zz_y[k] - g_y_z_zz_y[k] * ab_y + g_y_z_zz_yy[k];

                g_y_z_yzz_z[k] = -g_0_z_zz_z[k] - g_y_z_zz_z[k] * ab_y + g_y_z_zz_yz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzz_x = cbuffer.data(fp_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_z_zzz_y = cbuffer.data(fp_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_z_zzz_z = cbuffer.data(fp_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_x, g_y_0_zz_y, g_y_0_zz_z, g_y_z_zz_x, g_y_z_zz_xz, g_y_z_zz_y, g_y_z_zz_yz, g_y_z_zz_z, g_y_z_zz_zz, g_y_z_zzz_x, g_y_z_zzz_y, g_y_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzz_x[k] = g_y_0_zz_x[k] - g_y_z_zz_x[k] * ab_z + g_y_z_zz_xz[k];

                g_y_z_zzz_y[k] = g_y_0_zz_y[k] - g_y_z_zz_y[k] * ab_z + g_y_z_zz_yz[k];

                g_y_z_zzz_z[k] = g_y_0_zz_z[k] - g_y_z_zz_z[k] * ab_z + g_y_z_zz_zz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxx_x = cbuffer.data(fp_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_xxx_y = cbuffer.data(fp_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_xxx_z = cbuffer.data(fp_geom_11_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_x, g_z_0_xx_y, g_z_0_xx_z, g_z_x_xx_x, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_y, g_z_x_xx_z, g_z_x_xxx_x, g_z_x_xxx_y, g_z_x_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxx_x[k] = g_z_0_xx_x[k] - g_z_x_xx_x[k] * ab_x + g_z_x_xx_xx[k];

                g_z_x_xxx_y[k] = g_z_0_xx_y[k] - g_z_x_xx_y[k] * ab_x + g_z_x_xx_xy[k];

                g_z_x_xxx_z[k] = g_z_0_xx_z[k] - g_z_x_xx_z[k] * ab_x + g_z_x_xx_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxy_x = cbuffer.data(fp_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_xxy_y = cbuffer.data(fp_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_xxy_z = cbuffer.data(fp_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xx_x, g_z_x_xx_xy, g_z_x_xx_y, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_z, g_z_x_xxy_x, g_z_x_xxy_y, g_z_x_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxy_x[k] = -g_z_x_xx_x[k] * ab_y + g_z_x_xx_xy[k];

                g_z_x_xxy_y[k] = -g_z_x_xx_y[k] * ab_y + g_z_x_xx_yy[k];

                g_z_x_xxy_z[k] = -g_z_x_xx_z[k] * ab_y + g_z_x_xx_yz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxz_x = cbuffer.data(fp_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_xxz_y = cbuffer.data(fp_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_xxz_z = cbuffer.data(fp_geom_11_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_x, g_z_0_xz_y, g_z_0_xz_z, g_z_x_xxz_x, g_z_x_xxz_y, g_z_x_xxz_z, g_z_x_xz_x, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_y, g_z_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxz_x[k] = g_z_0_xz_x[k] - g_z_x_xz_x[k] * ab_x + g_z_x_xz_xx[k];

                g_z_x_xxz_y[k] = g_z_0_xz_y[k] - g_z_x_xz_y[k] * ab_x + g_z_x_xz_xy[k];

                g_z_x_xxz_z[k] = g_z_0_xz_z[k] - g_z_x_xz_z[k] * ab_x + g_z_x_xz_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyy_x = cbuffer.data(fp_geom_11_off + 189 * ccomps * dcomps);

            auto g_z_x_xyy_y = cbuffer.data(fp_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_xyy_z = cbuffer.data(fp_geom_11_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xy_x, g_z_x_xy_xy, g_z_x_xy_y, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_z, g_z_x_xyy_x, g_z_x_xyy_y, g_z_x_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyy_x[k] = -g_z_x_xy_x[k] * ab_y + g_z_x_xy_xy[k];

                g_z_x_xyy_y[k] = -g_z_x_xy_y[k] * ab_y + g_z_x_xy_yy[k];

                g_z_x_xyy_z[k] = -g_z_x_xy_z[k] * ab_y + g_z_x_xy_yz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyz_x = cbuffer.data(fp_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_xyz_y = cbuffer.data(fp_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_xyz_z = cbuffer.data(fp_geom_11_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyz_x, g_z_x_xyz_y, g_z_x_xyz_z, g_z_x_xz_x, g_z_x_xz_xy, g_z_x_xz_y, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyz_x[k] = -g_z_x_xz_x[k] * ab_y + g_z_x_xz_xy[k];

                g_z_x_xyz_y[k] = -g_z_x_xz_y[k] * ab_y + g_z_x_xz_yy[k];

                g_z_x_xyz_z[k] = -g_z_x_xz_z[k] * ab_y + g_z_x_xz_yz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzz_x = cbuffer.data(fp_geom_11_off + 195 * ccomps * dcomps);

            auto g_z_x_xzz_y = cbuffer.data(fp_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_x_xzz_z = cbuffer.data(fp_geom_11_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, g_z_x_xzz_x, g_z_x_xzz_y, g_z_x_xzz_z, g_z_x_zz_x, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_y, g_z_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzz_x[k] = g_z_0_zz_x[k] - g_z_x_zz_x[k] * ab_x + g_z_x_zz_xx[k];

                g_z_x_xzz_y[k] = g_z_0_zz_y[k] - g_z_x_zz_y[k] * ab_x + g_z_x_zz_xy[k];

                g_z_x_xzz_z[k] = g_z_0_zz_z[k] - g_z_x_zz_z[k] * ab_x + g_z_x_zz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyy_x = cbuffer.data(fp_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_x_yyy_y = cbuffer.data(fp_geom_11_off + 199 * ccomps * dcomps);

            auto g_z_x_yyy_z = cbuffer.data(fp_geom_11_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yy_x, g_z_x_yy_xy, g_z_x_yy_y, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_z, g_z_x_yyy_x, g_z_x_yyy_y, g_z_x_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyy_x[k] = -g_z_x_yy_x[k] * ab_y + g_z_x_yy_xy[k];

                g_z_x_yyy_y[k] = -g_z_x_yy_y[k] * ab_y + g_z_x_yy_yy[k];

                g_z_x_yyy_z[k] = -g_z_x_yy_z[k] * ab_y + g_z_x_yy_yz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyz_x = cbuffer.data(fp_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_x_yyz_y = cbuffer.data(fp_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_x_yyz_z = cbuffer.data(fp_geom_11_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyz_x, g_z_x_yyz_y, g_z_x_yyz_z, g_z_x_yz_x, g_z_x_yz_xy, g_z_x_yz_y, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyz_x[k] = -g_z_x_yz_x[k] * ab_y + g_z_x_yz_xy[k];

                g_z_x_yyz_y[k] = -g_z_x_yz_y[k] * ab_y + g_z_x_yz_yy[k];

                g_z_x_yyz_z[k] = -g_z_x_yz_z[k] * ab_y + g_z_x_yz_yz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzz_x = cbuffer.data(fp_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_x_yzz_y = cbuffer.data(fp_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_x_yzz_z = cbuffer.data(fp_geom_11_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzz_x, g_z_x_yzz_y, g_z_x_yzz_z, g_z_x_zz_x, g_z_x_zz_xy, g_z_x_zz_y, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzz_x[k] = -g_z_x_zz_x[k] * ab_y + g_z_x_zz_xy[k];

                g_z_x_yzz_y[k] = -g_z_x_zz_y[k] * ab_y + g_z_x_zz_yy[k];

                g_z_x_yzz_z[k] = -g_z_x_zz_z[k] * ab_y + g_z_x_zz_yz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzz_x = cbuffer.data(fp_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_x_zzz_y = cbuffer.data(fp_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_x_zzz_z = cbuffer.data(fp_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_x, g_0_x_zz_y, g_0_x_zz_z, g_z_x_zz_x, g_z_x_zz_xz, g_z_x_zz_y, g_z_x_zz_yz, g_z_x_zz_z, g_z_x_zz_zz, g_z_x_zzz_x, g_z_x_zzz_y, g_z_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzz_x[k] = -g_0_x_zz_x[k] - g_z_x_zz_x[k] * ab_z + g_z_x_zz_xz[k];

                g_z_x_zzz_y[k] = -g_0_x_zz_y[k] - g_z_x_zz_y[k] * ab_z + g_z_x_zz_yz[k];

                g_z_x_zzz_z[k] = -g_0_x_zz_z[k] - g_z_x_zz_z[k] * ab_z + g_z_x_zz_zz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxx_x = cbuffer.data(fp_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_xxx_y = cbuffer.data(fp_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_xxx_z = cbuffer.data(fp_geom_11_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xx_x, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_y, g_z_y_xx_z, g_z_y_xxx_x, g_z_y_xxx_y, g_z_y_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxx_x[k] = -g_z_y_xx_x[k] * ab_x + g_z_y_xx_xx[k];

                g_z_y_xxx_y[k] = -g_z_y_xx_y[k] * ab_x + g_z_y_xx_xy[k];

                g_z_y_xxx_z[k] = -g_z_y_xx_z[k] * ab_x + g_z_y_xx_xz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxy_x = cbuffer.data(fp_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_xxy_y = cbuffer.data(fp_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_xxy_z = cbuffer.data(fp_geom_11_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxy_x, g_z_y_xxy_y, g_z_y_xxy_z, g_z_y_xy_x, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_y, g_z_y_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxy_x[k] = -g_z_y_xy_x[k] * ab_x + g_z_y_xy_xx[k];

                g_z_y_xxy_y[k] = -g_z_y_xy_y[k] * ab_x + g_z_y_xy_xy[k];

                g_z_y_xxy_z[k] = -g_z_y_xy_z[k] * ab_x + g_z_y_xy_xz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxz_x = cbuffer.data(fp_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_xxz_y = cbuffer.data(fp_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_xxz_z = cbuffer.data(fp_geom_11_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxz_x, g_z_y_xxz_y, g_z_y_xxz_z, g_z_y_xz_x, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_y, g_z_y_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxz_x[k] = -g_z_y_xz_x[k] * ab_x + g_z_y_xz_xx[k];

                g_z_y_xxz_y[k] = -g_z_y_xz_y[k] * ab_x + g_z_y_xz_xy[k];

                g_z_y_xxz_z[k] = -g_z_y_xz_z[k] * ab_x + g_z_y_xz_xz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyy_x = cbuffer.data(fp_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_y_xyy_y = cbuffer.data(fp_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_xyy_z = cbuffer.data(fp_geom_11_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyy_x, g_z_y_xyy_y, g_z_y_xyy_z, g_z_y_yy_x, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_y, g_z_y_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyy_x[k] = -g_z_y_yy_x[k] * ab_x + g_z_y_yy_xx[k];

                g_z_y_xyy_y[k] = -g_z_y_yy_y[k] * ab_x + g_z_y_yy_xy[k];

                g_z_y_xyy_z[k] = -g_z_y_yy_z[k] * ab_x + g_z_y_yy_xz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyz_x = cbuffer.data(fp_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_xyz_y = cbuffer.data(fp_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_y_xyz_z = cbuffer.data(fp_geom_11_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyz_x, g_z_y_xyz_y, g_z_y_xyz_z, g_z_y_yz_x, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_y, g_z_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyz_x[k] = -g_z_y_yz_x[k] * ab_x + g_z_y_yz_xx[k];

                g_z_y_xyz_y[k] = -g_z_y_yz_y[k] * ab_x + g_z_y_yz_xy[k];

                g_z_y_xyz_z[k] = -g_z_y_yz_z[k] * ab_x + g_z_y_yz_xz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzz_x = cbuffer.data(fp_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_y_xzz_y = cbuffer.data(fp_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_y_xzz_z = cbuffer.data(fp_geom_11_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzz_x, g_z_y_xzz_y, g_z_y_xzz_z, g_z_y_zz_x, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_y, g_z_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzz_x[k] = -g_z_y_zz_x[k] * ab_x + g_z_y_zz_xx[k];

                g_z_y_xzz_y[k] = -g_z_y_zz_y[k] * ab_x + g_z_y_zz_xy[k];

                g_z_y_xzz_z[k] = -g_z_y_zz_z[k] * ab_x + g_z_y_zz_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyy_x = cbuffer.data(fp_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_y_yyy_y = cbuffer.data(fp_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_y_yyy_z = cbuffer.data(fp_geom_11_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_x, g_z_0_yy_y, g_z_0_yy_z, g_z_y_yy_x, g_z_y_yy_xy, g_z_y_yy_y, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_z, g_z_y_yyy_x, g_z_y_yyy_y, g_z_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyy_x[k] = g_z_0_yy_x[k] - g_z_y_yy_x[k] * ab_y + g_z_y_yy_xy[k];

                g_z_y_yyy_y[k] = g_z_0_yy_y[k] - g_z_y_yy_y[k] * ab_y + g_z_y_yy_yy[k];

                g_z_y_yyy_z[k] = g_z_0_yy_z[k] - g_z_y_yy_z[k] * ab_y + g_z_y_yy_yz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyz_x = cbuffer.data(fp_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_y_yyz_y = cbuffer.data(fp_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_y_yyz_z = cbuffer.data(fp_geom_11_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_x, g_z_0_yz_y, g_z_0_yz_z, g_z_y_yyz_x, g_z_y_yyz_y, g_z_y_yyz_z, g_z_y_yz_x, g_z_y_yz_xy, g_z_y_yz_y, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyz_x[k] = g_z_0_yz_x[k] - g_z_y_yz_x[k] * ab_y + g_z_y_yz_xy[k];

                g_z_y_yyz_y[k] = g_z_0_yz_y[k] - g_z_y_yz_y[k] * ab_y + g_z_y_yz_yy[k];

                g_z_y_yyz_z[k] = g_z_0_yz_z[k] - g_z_y_yz_z[k] * ab_y + g_z_y_yz_yz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzz_x = cbuffer.data(fp_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_y_yzz_y = cbuffer.data(fp_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_y_yzz_z = cbuffer.data(fp_geom_11_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, g_z_y_yzz_x, g_z_y_yzz_y, g_z_y_yzz_z, g_z_y_zz_x, g_z_y_zz_xy, g_z_y_zz_y, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzz_x[k] = g_z_0_zz_x[k] - g_z_y_zz_x[k] * ab_y + g_z_y_zz_xy[k];

                g_z_y_yzz_y[k] = g_z_0_zz_y[k] - g_z_y_zz_y[k] * ab_y + g_z_y_zz_yy[k];

                g_z_y_yzz_z[k] = g_z_0_zz_z[k] - g_z_y_zz_z[k] * ab_y + g_z_y_zz_yz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzz_x = cbuffer.data(fp_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_y_zzz_y = cbuffer.data(fp_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_y_zzz_z = cbuffer.data(fp_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_x, g_0_y_zz_y, g_0_y_zz_z, g_z_y_zz_x, g_z_y_zz_xz, g_z_y_zz_y, g_z_y_zz_yz, g_z_y_zz_z, g_z_y_zz_zz, g_z_y_zzz_x, g_z_y_zzz_y, g_z_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzz_x[k] = -g_0_y_zz_x[k] - g_z_y_zz_x[k] * ab_z + g_z_y_zz_xz[k];

                g_z_y_zzz_y[k] = -g_0_y_zz_y[k] - g_z_y_zz_y[k] * ab_z + g_z_y_zz_yz[k];

                g_z_y_zzz_z[k] = -g_0_y_zz_z[k] - g_z_y_zz_z[k] * ab_z + g_z_y_zz_zz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxx_x = cbuffer.data(fp_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_xxx_y = cbuffer.data(fp_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_xxx_z = cbuffer.data(fp_geom_11_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xx_x, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_y, g_z_z_xx_z, g_z_z_xxx_x, g_z_z_xxx_y, g_z_z_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxx_x[k] = -g_z_z_xx_x[k] * ab_x + g_z_z_xx_xx[k];

                g_z_z_xxx_y[k] = -g_z_z_xx_y[k] * ab_x + g_z_z_xx_xy[k];

                g_z_z_xxx_z[k] = -g_z_z_xx_z[k] * ab_x + g_z_z_xx_xz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxy_x = cbuffer.data(fp_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_xxy_y = cbuffer.data(fp_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_xxy_z = cbuffer.data(fp_geom_11_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxy_x, g_z_z_xxy_y, g_z_z_xxy_z, g_z_z_xy_x, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_y, g_z_z_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxy_x[k] = -g_z_z_xy_x[k] * ab_x + g_z_z_xy_xx[k];

                g_z_z_xxy_y[k] = -g_z_z_xy_y[k] * ab_x + g_z_z_xy_xy[k];

                g_z_z_xxy_z[k] = -g_z_z_xy_z[k] * ab_x + g_z_z_xy_xz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxz_x = cbuffer.data(fp_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_xxz_y = cbuffer.data(fp_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_xxz_z = cbuffer.data(fp_geom_11_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxz_x, g_z_z_xxz_y, g_z_z_xxz_z, g_z_z_xz_x, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_y, g_z_z_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxz_x[k] = -g_z_z_xz_x[k] * ab_x + g_z_z_xz_xx[k];

                g_z_z_xxz_y[k] = -g_z_z_xz_y[k] * ab_x + g_z_z_xz_xy[k];

                g_z_z_xxz_z[k] = -g_z_z_xz_z[k] * ab_x + g_z_z_xz_xz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyy_x = cbuffer.data(fp_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_z_xyy_y = cbuffer.data(fp_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_xyy_z = cbuffer.data(fp_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyy_x, g_z_z_xyy_y, g_z_z_xyy_z, g_z_z_yy_x, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_y, g_z_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyy_x[k] = -g_z_z_yy_x[k] * ab_x + g_z_z_yy_xx[k];

                g_z_z_xyy_y[k] = -g_z_z_yy_y[k] * ab_x + g_z_z_yy_xy[k];

                g_z_z_xyy_z[k] = -g_z_z_yy_z[k] * ab_x + g_z_z_yy_xz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyz_x = cbuffer.data(fp_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_z_xyz_y = cbuffer.data(fp_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_z_xyz_z = cbuffer.data(fp_geom_11_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyz_x, g_z_z_xyz_y, g_z_z_xyz_z, g_z_z_yz_x, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_y, g_z_z_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyz_x[k] = -g_z_z_yz_x[k] * ab_x + g_z_z_yz_xx[k];

                g_z_z_xyz_y[k] = -g_z_z_yz_y[k] * ab_x + g_z_z_yz_xy[k];

                g_z_z_xyz_z[k] = -g_z_z_yz_z[k] * ab_x + g_z_z_yz_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzz_x = cbuffer.data(fp_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_z_xzz_y = cbuffer.data(fp_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_z_xzz_z = cbuffer.data(fp_geom_11_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzz_x, g_z_z_xzz_y, g_z_z_xzz_z, g_z_z_zz_x, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_y, g_z_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzz_x[k] = -g_z_z_zz_x[k] * ab_x + g_z_z_zz_xx[k];

                g_z_z_xzz_y[k] = -g_z_z_zz_y[k] * ab_x + g_z_z_zz_xy[k];

                g_z_z_xzz_z[k] = -g_z_z_zz_z[k] * ab_x + g_z_z_zz_xz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyy_x = cbuffer.data(fp_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_z_yyy_y = cbuffer.data(fp_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_z_yyy_z = cbuffer.data(fp_geom_11_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yy_x, g_z_z_yy_xy, g_z_z_yy_y, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_z, g_z_z_yyy_x, g_z_z_yyy_y, g_z_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyy_x[k] = -g_z_z_yy_x[k] * ab_y + g_z_z_yy_xy[k];

                g_z_z_yyy_y[k] = -g_z_z_yy_y[k] * ab_y + g_z_z_yy_yy[k];

                g_z_z_yyy_z[k] = -g_z_z_yy_z[k] * ab_y + g_z_z_yy_yz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyz_x = cbuffer.data(fp_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_z_yyz_y = cbuffer.data(fp_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_z_yyz_z = cbuffer.data(fp_geom_11_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyz_x, g_z_z_yyz_y, g_z_z_yyz_z, g_z_z_yz_x, g_z_z_yz_xy, g_z_z_yz_y, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyz_x[k] = -g_z_z_yz_x[k] * ab_y + g_z_z_yz_xy[k];

                g_z_z_yyz_y[k] = -g_z_z_yz_y[k] * ab_y + g_z_z_yz_yy[k];

                g_z_z_yyz_z[k] = -g_z_z_yz_z[k] * ab_y + g_z_z_yz_yz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzz_x = cbuffer.data(fp_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_z_yzz_y = cbuffer.data(fp_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_z_yzz_z = cbuffer.data(fp_geom_11_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzz_x, g_z_z_yzz_y, g_z_z_yzz_z, g_z_z_zz_x, g_z_z_zz_xy, g_z_z_zz_y, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzz_x[k] = -g_z_z_zz_x[k] * ab_y + g_z_z_zz_xy[k];

                g_z_z_yzz_y[k] = -g_z_z_zz_y[k] * ab_y + g_z_z_zz_yy[k];

                g_z_z_yzz_z[k] = -g_z_z_zz_z[k] * ab_y + g_z_z_zz_yz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzz_x = cbuffer.data(fp_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_z_zzz_y = cbuffer.data(fp_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_z_zzz_z = cbuffer.data(fp_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, g_z_z_zz_x, g_z_z_zz_xz, g_z_z_zz_y, g_z_z_zz_yz, g_z_z_zz_z, g_z_z_zz_zz, g_z_z_zzz_x, g_z_z_zzz_y, g_z_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzz_x[k] = -g_0_z_zz_x[k] + g_z_0_zz_x[k] - g_z_z_zz_x[k] * ab_z + g_z_z_zz_xz[k];

                g_z_z_zzz_y[k] = -g_0_z_zz_y[k] + g_z_0_zz_y[k] - g_z_z_zz_y[k] * ab_z + g_z_z_zz_yz[k];

                g_z_z_zzz_z[k] = -g_0_z_zz_z[k] + g_z_0_zz_z[k] - g_z_z_zz_z[k] * ab_z + g_z_z_zz_zz[k];
            }
        }
    }
}

} // erirec namespace

