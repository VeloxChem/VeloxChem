#include "ElectronRepulsionGeom2000ContrRecHPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_hpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_hpxx,
                                            const size_t idx_geom_10_gpxx,
                                            const size_t idx_geom_20_gpxx,
                                            const size_t idx_geom_20_gdxx,
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
            /// Set up components of auxilary buffer : GPSS

            const auto gp_geom_10_off = idx_geom_10_gpxx + i * dcomps + j;

            auto g_x_0_xxxx_x = cbuffer.data(gp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_y = cbuffer.data(gp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_z = cbuffer.data(gp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxy_x = cbuffer.data(gp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxy_y = cbuffer.data(gp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxy_z = cbuffer.data(gp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxz_x = cbuffer.data(gp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxz_y = cbuffer.data(gp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxz_z = cbuffer.data(gp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxyy_x = cbuffer.data(gp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxyy_y = cbuffer.data(gp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxyy_z = cbuffer.data(gp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxyz_x = cbuffer.data(gp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxyz_y = cbuffer.data(gp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxyz_z = cbuffer.data(gp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxzz_x = cbuffer.data(gp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxzz_y = cbuffer.data(gp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxzz_z = cbuffer.data(gp_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xyyy_x = cbuffer.data(gp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xyyy_y = cbuffer.data(gp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xyyy_z = cbuffer.data(gp_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xyyz_x = cbuffer.data(gp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xyyz_y = cbuffer.data(gp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xyyz_z = cbuffer.data(gp_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xyzz_x = cbuffer.data(gp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xyzz_y = cbuffer.data(gp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xyzz_z = cbuffer.data(gp_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xzzz_x = cbuffer.data(gp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xzzz_y = cbuffer.data(gp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xzzz_z = cbuffer.data(gp_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_yyyy_x = cbuffer.data(gp_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_yyyy_y = cbuffer.data(gp_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_yyyy_z = cbuffer.data(gp_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_yyyz_x = cbuffer.data(gp_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_yyyz_y = cbuffer.data(gp_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_yyyz_z = cbuffer.data(gp_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_yyzz_x = cbuffer.data(gp_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yyzz_y = cbuffer.data(gp_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yyzz_z = cbuffer.data(gp_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_yzzz_x = cbuffer.data(gp_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yzzz_y = cbuffer.data(gp_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yzzz_z = cbuffer.data(gp_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_zzzz_x = cbuffer.data(gp_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_zzzz_y = cbuffer.data(gp_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_zzzz_z = cbuffer.data(gp_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_xxxx_x = cbuffer.data(gp_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xxxx_y = cbuffer.data(gp_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xxxx_z = cbuffer.data(gp_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_xxxy_x = cbuffer.data(gp_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_xxxy_y = cbuffer.data(gp_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_xxxy_z = cbuffer.data(gp_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_xxxz_x = cbuffer.data(gp_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_xxxz_y = cbuffer.data(gp_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_xxxz_z = cbuffer.data(gp_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_xxyy_x = cbuffer.data(gp_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_xxyy_y = cbuffer.data(gp_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_xxyy_z = cbuffer.data(gp_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_xxyz_x = cbuffer.data(gp_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_xxyz_y = cbuffer.data(gp_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_xxyz_z = cbuffer.data(gp_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_xxzz_x = cbuffer.data(gp_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xxzz_y = cbuffer.data(gp_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xxzz_z = cbuffer.data(gp_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xyyy_x = cbuffer.data(gp_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xyyy_y = cbuffer.data(gp_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xyyy_z = cbuffer.data(gp_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xyyz_x = cbuffer.data(gp_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xyyz_y = cbuffer.data(gp_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xyyz_z = cbuffer.data(gp_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xyzz_x = cbuffer.data(gp_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xyzz_y = cbuffer.data(gp_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xyzz_z = cbuffer.data(gp_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xzzz_x = cbuffer.data(gp_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xzzz_y = cbuffer.data(gp_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xzzz_z = cbuffer.data(gp_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_yyyy_x = cbuffer.data(gp_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_yyyy_y = cbuffer.data(gp_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_yyyy_z = cbuffer.data(gp_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_yyyz_x = cbuffer.data(gp_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_yyyz_y = cbuffer.data(gp_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_yyyz_z = cbuffer.data(gp_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_yyzz_x = cbuffer.data(gp_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_yyzz_y = cbuffer.data(gp_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_yyzz_z = cbuffer.data(gp_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_yzzz_x = cbuffer.data(gp_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_yzzz_y = cbuffer.data(gp_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_yzzz_z = cbuffer.data(gp_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_zzzz_x = cbuffer.data(gp_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_zzzz_y = cbuffer.data(gp_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_zzzz_z = cbuffer.data(gp_geom_10_off + 89 * ccomps * dcomps);

            auto g_z_0_xxxx_x = cbuffer.data(gp_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_xxxx_y = cbuffer.data(gp_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_xxxx_z = cbuffer.data(gp_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_xxxy_x = cbuffer.data(gp_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_xxxy_y = cbuffer.data(gp_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_xxxy_z = cbuffer.data(gp_geom_10_off + 95 * ccomps * dcomps);

            auto g_z_0_xxxz_x = cbuffer.data(gp_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_xxxz_y = cbuffer.data(gp_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_xxxz_z = cbuffer.data(gp_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_xxyy_x = cbuffer.data(gp_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_xxyy_y = cbuffer.data(gp_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_xxyy_z = cbuffer.data(gp_geom_10_off + 101 * ccomps * dcomps);

            auto g_z_0_xxyz_x = cbuffer.data(gp_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_xxyz_y = cbuffer.data(gp_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_xxyz_z = cbuffer.data(gp_geom_10_off + 104 * ccomps * dcomps);

            auto g_z_0_xxzz_x = cbuffer.data(gp_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_xxzz_y = cbuffer.data(gp_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_xxzz_z = cbuffer.data(gp_geom_10_off + 107 * ccomps * dcomps);

            auto g_z_0_xyyy_x = cbuffer.data(gp_geom_10_off + 108 * ccomps * dcomps);

            auto g_z_0_xyyy_y = cbuffer.data(gp_geom_10_off + 109 * ccomps * dcomps);

            auto g_z_0_xyyy_z = cbuffer.data(gp_geom_10_off + 110 * ccomps * dcomps);

            auto g_z_0_xyyz_x = cbuffer.data(gp_geom_10_off + 111 * ccomps * dcomps);

            auto g_z_0_xyyz_y = cbuffer.data(gp_geom_10_off + 112 * ccomps * dcomps);

            auto g_z_0_xyyz_z = cbuffer.data(gp_geom_10_off + 113 * ccomps * dcomps);

            auto g_z_0_xyzz_x = cbuffer.data(gp_geom_10_off + 114 * ccomps * dcomps);

            auto g_z_0_xyzz_y = cbuffer.data(gp_geom_10_off + 115 * ccomps * dcomps);

            auto g_z_0_xyzz_z = cbuffer.data(gp_geom_10_off + 116 * ccomps * dcomps);

            auto g_z_0_xzzz_x = cbuffer.data(gp_geom_10_off + 117 * ccomps * dcomps);

            auto g_z_0_xzzz_y = cbuffer.data(gp_geom_10_off + 118 * ccomps * dcomps);

            auto g_z_0_xzzz_z = cbuffer.data(gp_geom_10_off + 119 * ccomps * dcomps);

            auto g_z_0_yyyy_x = cbuffer.data(gp_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_yyyy_y = cbuffer.data(gp_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_yyyy_z = cbuffer.data(gp_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_yyyz_x = cbuffer.data(gp_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_yyyz_y = cbuffer.data(gp_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_yyyz_z = cbuffer.data(gp_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_yyzz_x = cbuffer.data(gp_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_yyzz_y = cbuffer.data(gp_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_yyzz_z = cbuffer.data(gp_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_yzzz_x = cbuffer.data(gp_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_yzzz_y = cbuffer.data(gp_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_yzzz_z = cbuffer.data(gp_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_zzzz_x = cbuffer.data(gp_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_zzzz_y = cbuffer.data(gp_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_zzzz_z = cbuffer.data(gp_geom_10_off + 134 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GPSS

            const auto gp_geom_20_off = idx_geom_20_gpxx + i * dcomps + j;

            auto g_xx_0_xxxx_x = cbuffer.data(gp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_y = cbuffer.data(gp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_z = cbuffer.data(gp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxy_x = cbuffer.data(gp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxy_y = cbuffer.data(gp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxy_z = cbuffer.data(gp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxz_x = cbuffer.data(gp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxz_y = cbuffer.data(gp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxz_z = cbuffer.data(gp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxyy_x = cbuffer.data(gp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxyy_y = cbuffer.data(gp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxyy_z = cbuffer.data(gp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxyz_x = cbuffer.data(gp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxyz_y = cbuffer.data(gp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxyz_z = cbuffer.data(gp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxzz_x = cbuffer.data(gp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxzz_y = cbuffer.data(gp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxzz_z = cbuffer.data(gp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xyyy_x = cbuffer.data(gp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xyyy_y = cbuffer.data(gp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xyyy_z = cbuffer.data(gp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xyyz_x = cbuffer.data(gp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xyyz_y = cbuffer.data(gp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xyyz_z = cbuffer.data(gp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xyzz_x = cbuffer.data(gp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xyzz_y = cbuffer.data(gp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xyzz_z = cbuffer.data(gp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xzzz_x = cbuffer.data(gp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xzzz_y = cbuffer.data(gp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xzzz_z = cbuffer.data(gp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_yyyy_x = cbuffer.data(gp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_yyyy_y = cbuffer.data(gp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_yyyy_z = cbuffer.data(gp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_yyyz_x = cbuffer.data(gp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_yyyz_y = cbuffer.data(gp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_yyyz_z = cbuffer.data(gp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_yyzz_x = cbuffer.data(gp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_yyzz_y = cbuffer.data(gp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_yyzz_z = cbuffer.data(gp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_yzzz_x = cbuffer.data(gp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_yzzz_y = cbuffer.data(gp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_yzzz_z = cbuffer.data(gp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_zzzz_x = cbuffer.data(gp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_zzzz_y = cbuffer.data(gp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_zzzz_z = cbuffer.data(gp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_xxxx_x = cbuffer.data(gp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_xxxx_y = cbuffer.data(gp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_xxxx_z = cbuffer.data(gp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_xxxy_x = cbuffer.data(gp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_xxxy_y = cbuffer.data(gp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_xxxy_z = cbuffer.data(gp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_xxxz_x = cbuffer.data(gp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_xxxz_y = cbuffer.data(gp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_xxxz_z = cbuffer.data(gp_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_xxyy_x = cbuffer.data(gp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_xxyy_y = cbuffer.data(gp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_xxyy_z = cbuffer.data(gp_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_xxyz_x = cbuffer.data(gp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_xxyz_y = cbuffer.data(gp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_xxyz_z = cbuffer.data(gp_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_xxzz_x = cbuffer.data(gp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_xxzz_y = cbuffer.data(gp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_xxzz_z = cbuffer.data(gp_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_xyyy_x = cbuffer.data(gp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xyyy_y = cbuffer.data(gp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xyyy_z = cbuffer.data(gp_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_xyyz_x = cbuffer.data(gp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xyyz_y = cbuffer.data(gp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xyyz_z = cbuffer.data(gp_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_xyzz_x = cbuffer.data(gp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xyzz_y = cbuffer.data(gp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xyzz_z = cbuffer.data(gp_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_xzzz_x = cbuffer.data(gp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xzzz_y = cbuffer.data(gp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xzzz_z = cbuffer.data(gp_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_yyyy_x = cbuffer.data(gp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_yyyy_y = cbuffer.data(gp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_yyyy_z = cbuffer.data(gp_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_yyyz_x = cbuffer.data(gp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_yyyz_y = cbuffer.data(gp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_yyyz_z = cbuffer.data(gp_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_yyzz_x = cbuffer.data(gp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_yyzz_y = cbuffer.data(gp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_yyzz_z = cbuffer.data(gp_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_yzzz_x = cbuffer.data(gp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_yzzz_y = cbuffer.data(gp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_yzzz_z = cbuffer.data(gp_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_zzzz_x = cbuffer.data(gp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_zzzz_y = cbuffer.data(gp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_zzzz_z = cbuffer.data(gp_geom_20_off + 89 * ccomps * dcomps);

            auto g_xz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 92 * ccomps * dcomps);

            auto g_xz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 95 * ccomps * dcomps);

            auto g_xz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 98 * ccomps * dcomps);

            auto g_xz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 101 * ccomps * dcomps);

            auto g_xz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 104 * ccomps * dcomps);

            auto g_xz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 107 * ccomps * dcomps);

            auto g_xz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 108 * ccomps * dcomps);

            auto g_xz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 109 * ccomps * dcomps);

            auto g_xz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 110 * ccomps * dcomps);

            auto g_xz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 111 * ccomps * dcomps);

            auto g_xz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 112 * ccomps * dcomps);

            auto g_xz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 113 * ccomps * dcomps);

            auto g_xz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 114 * ccomps * dcomps);

            auto g_xz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 115 * ccomps * dcomps);

            auto g_xz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 116 * ccomps * dcomps);

            auto g_xz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 117 * ccomps * dcomps);

            auto g_xz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 118 * ccomps * dcomps);

            auto g_xz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 119 * ccomps * dcomps);

            auto g_xz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 122 * ccomps * dcomps);

            auto g_xz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 134 * ccomps * dcomps);

            auto g_yy_0_xxxx_x = cbuffer.data(gp_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_xxxx_y = cbuffer.data(gp_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_xxxx_z = cbuffer.data(gp_geom_20_off + 137 * ccomps * dcomps);

            auto g_yy_0_xxxy_x = cbuffer.data(gp_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_xxxy_y = cbuffer.data(gp_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_xxxy_z = cbuffer.data(gp_geom_20_off + 140 * ccomps * dcomps);

            auto g_yy_0_xxxz_x = cbuffer.data(gp_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_xxxz_y = cbuffer.data(gp_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_xxxz_z = cbuffer.data(gp_geom_20_off + 143 * ccomps * dcomps);

            auto g_yy_0_xxyy_x = cbuffer.data(gp_geom_20_off + 144 * ccomps * dcomps);

            auto g_yy_0_xxyy_y = cbuffer.data(gp_geom_20_off + 145 * ccomps * dcomps);

            auto g_yy_0_xxyy_z = cbuffer.data(gp_geom_20_off + 146 * ccomps * dcomps);

            auto g_yy_0_xxyz_x = cbuffer.data(gp_geom_20_off + 147 * ccomps * dcomps);

            auto g_yy_0_xxyz_y = cbuffer.data(gp_geom_20_off + 148 * ccomps * dcomps);

            auto g_yy_0_xxyz_z = cbuffer.data(gp_geom_20_off + 149 * ccomps * dcomps);

            auto g_yy_0_xxzz_x = cbuffer.data(gp_geom_20_off + 150 * ccomps * dcomps);

            auto g_yy_0_xxzz_y = cbuffer.data(gp_geom_20_off + 151 * ccomps * dcomps);

            auto g_yy_0_xxzz_z = cbuffer.data(gp_geom_20_off + 152 * ccomps * dcomps);

            auto g_yy_0_xyyy_x = cbuffer.data(gp_geom_20_off + 153 * ccomps * dcomps);

            auto g_yy_0_xyyy_y = cbuffer.data(gp_geom_20_off + 154 * ccomps * dcomps);

            auto g_yy_0_xyyy_z = cbuffer.data(gp_geom_20_off + 155 * ccomps * dcomps);

            auto g_yy_0_xyyz_x = cbuffer.data(gp_geom_20_off + 156 * ccomps * dcomps);

            auto g_yy_0_xyyz_y = cbuffer.data(gp_geom_20_off + 157 * ccomps * dcomps);

            auto g_yy_0_xyyz_z = cbuffer.data(gp_geom_20_off + 158 * ccomps * dcomps);

            auto g_yy_0_xyzz_x = cbuffer.data(gp_geom_20_off + 159 * ccomps * dcomps);

            auto g_yy_0_xyzz_y = cbuffer.data(gp_geom_20_off + 160 * ccomps * dcomps);

            auto g_yy_0_xyzz_z = cbuffer.data(gp_geom_20_off + 161 * ccomps * dcomps);

            auto g_yy_0_xzzz_x = cbuffer.data(gp_geom_20_off + 162 * ccomps * dcomps);

            auto g_yy_0_xzzz_y = cbuffer.data(gp_geom_20_off + 163 * ccomps * dcomps);

            auto g_yy_0_xzzz_z = cbuffer.data(gp_geom_20_off + 164 * ccomps * dcomps);

            auto g_yy_0_yyyy_x = cbuffer.data(gp_geom_20_off + 165 * ccomps * dcomps);

            auto g_yy_0_yyyy_y = cbuffer.data(gp_geom_20_off + 166 * ccomps * dcomps);

            auto g_yy_0_yyyy_z = cbuffer.data(gp_geom_20_off + 167 * ccomps * dcomps);

            auto g_yy_0_yyyz_x = cbuffer.data(gp_geom_20_off + 168 * ccomps * dcomps);

            auto g_yy_0_yyyz_y = cbuffer.data(gp_geom_20_off + 169 * ccomps * dcomps);

            auto g_yy_0_yyyz_z = cbuffer.data(gp_geom_20_off + 170 * ccomps * dcomps);

            auto g_yy_0_yyzz_x = cbuffer.data(gp_geom_20_off + 171 * ccomps * dcomps);

            auto g_yy_0_yyzz_y = cbuffer.data(gp_geom_20_off + 172 * ccomps * dcomps);

            auto g_yy_0_yyzz_z = cbuffer.data(gp_geom_20_off + 173 * ccomps * dcomps);

            auto g_yy_0_yzzz_x = cbuffer.data(gp_geom_20_off + 174 * ccomps * dcomps);

            auto g_yy_0_yzzz_y = cbuffer.data(gp_geom_20_off + 175 * ccomps * dcomps);

            auto g_yy_0_yzzz_z = cbuffer.data(gp_geom_20_off + 176 * ccomps * dcomps);

            auto g_yy_0_zzzz_x = cbuffer.data(gp_geom_20_off + 177 * ccomps * dcomps);

            auto g_yy_0_zzzz_y = cbuffer.data(gp_geom_20_off + 178 * ccomps * dcomps);

            auto g_yy_0_zzzz_z = cbuffer.data(gp_geom_20_off + 179 * ccomps * dcomps);

            auto g_yz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 180 * ccomps * dcomps);

            auto g_yz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 181 * ccomps * dcomps);

            auto g_yz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 182 * ccomps * dcomps);

            auto g_yz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 183 * ccomps * dcomps);

            auto g_yz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 184 * ccomps * dcomps);

            auto g_yz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 185 * ccomps * dcomps);

            auto g_yz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 186 * ccomps * dcomps);

            auto g_yz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 187 * ccomps * dcomps);

            auto g_yz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 188 * ccomps * dcomps);

            auto g_yz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 189 * ccomps * dcomps);

            auto g_yz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 190 * ccomps * dcomps);

            auto g_yz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 191 * ccomps * dcomps);

            auto g_yz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 192 * ccomps * dcomps);

            auto g_yz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 193 * ccomps * dcomps);

            auto g_yz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 194 * ccomps * dcomps);

            auto g_yz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 195 * ccomps * dcomps);

            auto g_yz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 196 * ccomps * dcomps);

            auto g_yz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 197 * ccomps * dcomps);

            auto g_yz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 198 * ccomps * dcomps);

            auto g_yz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 199 * ccomps * dcomps);

            auto g_yz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 200 * ccomps * dcomps);

            auto g_yz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 201 * ccomps * dcomps);

            auto g_yz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 202 * ccomps * dcomps);

            auto g_yz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 203 * ccomps * dcomps);

            auto g_yz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 204 * ccomps * dcomps);

            auto g_yz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 205 * ccomps * dcomps);

            auto g_yz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 206 * ccomps * dcomps);

            auto g_yz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 207 * ccomps * dcomps);

            auto g_yz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 208 * ccomps * dcomps);

            auto g_yz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 209 * ccomps * dcomps);

            auto g_yz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 210 * ccomps * dcomps);

            auto g_yz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 211 * ccomps * dcomps);

            auto g_yz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 212 * ccomps * dcomps);

            auto g_yz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 213 * ccomps * dcomps);

            auto g_yz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 214 * ccomps * dcomps);

            auto g_yz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 215 * ccomps * dcomps);

            auto g_yz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 216 * ccomps * dcomps);

            auto g_yz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 217 * ccomps * dcomps);

            auto g_yz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 218 * ccomps * dcomps);

            auto g_yz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 219 * ccomps * dcomps);

            auto g_yz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 220 * ccomps * dcomps);

            auto g_yz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 221 * ccomps * dcomps);

            auto g_yz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 222 * ccomps * dcomps);

            auto g_yz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 223 * ccomps * dcomps);

            auto g_yz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 224 * ccomps * dcomps);

            auto g_zz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 225 * ccomps * dcomps);

            auto g_zz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 226 * ccomps * dcomps);

            auto g_zz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 227 * ccomps * dcomps);

            auto g_zz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 228 * ccomps * dcomps);

            auto g_zz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 229 * ccomps * dcomps);

            auto g_zz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 230 * ccomps * dcomps);

            auto g_zz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 231 * ccomps * dcomps);

            auto g_zz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 232 * ccomps * dcomps);

            auto g_zz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 233 * ccomps * dcomps);

            auto g_zz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 234 * ccomps * dcomps);

            auto g_zz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 235 * ccomps * dcomps);

            auto g_zz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 236 * ccomps * dcomps);

            auto g_zz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 237 * ccomps * dcomps);

            auto g_zz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 238 * ccomps * dcomps);

            auto g_zz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 239 * ccomps * dcomps);

            auto g_zz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 240 * ccomps * dcomps);

            auto g_zz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 241 * ccomps * dcomps);

            auto g_zz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 242 * ccomps * dcomps);

            auto g_zz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 243 * ccomps * dcomps);

            auto g_zz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 244 * ccomps * dcomps);

            auto g_zz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 245 * ccomps * dcomps);

            auto g_zz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 246 * ccomps * dcomps);

            auto g_zz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 247 * ccomps * dcomps);

            auto g_zz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 248 * ccomps * dcomps);

            auto g_zz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 249 * ccomps * dcomps);

            auto g_zz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 250 * ccomps * dcomps);

            auto g_zz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 251 * ccomps * dcomps);

            auto g_zz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 252 * ccomps * dcomps);

            auto g_zz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 253 * ccomps * dcomps);

            auto g_zz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 254 * ccomps * dcomps);

            auto g_zz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 255 * ccomps * dcomps);

            auto g_zz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 256 * ccomps * dcomps);

            auto g_zz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 257 * ccomps * dcomps);

            auto g_zz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 258 * ccomps * dcomps);

            auto g_zz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 259 * ccomps * dcomps);

            auto g_zz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 260 * ccomps * dcomps);

            auto g_zz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 261 * ccomps * dcomps);

            auto g_zz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 262 * ccomps * dcomps);

            auto g_zz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 263 * ccomps * dcomps);

            auto g_zz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 264 * ccomps * dcomps);

            auto g_zz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 265 * ccomps * dcomps);

            auto g_zz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 266 * ccomps * dcomps);

            auto g_zz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 267 * ccomps * dcomps);

            auto g_zz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 268 * ccomps * dcomps);

            auto g_zz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GDSS

            const auto gd_geom_20_off = idx_geom_20_gdxx + i * dcomps + j;

            auto g_xx_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 191 * ccomps * dcomps);

            auto g_xz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 197 * ccomps * dcomps);

            auto g_xz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 269 * ccomps * dcomps);

            auto g_yy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 275 * ccomps * dcomps);

            auto g_yy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 281 * ccomps * dcomps);

            auto g_yy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 284 * ccomps * dcomps);

            auto g_yy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 287 * ccomps * dcomps);

            auto g_yy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 293 * ccomps * dcomps);

            auto g_yy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 359 * ccomps * dcomps);

            auto g_yz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 365 * ccomps * dcomps);

            auto g_yz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 371 * ccomps * dcomps);

            auto g_yz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 374 * ccomps * dcomps);

            auto g_yz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 377 * ccomps * dcomps);

            auto g_yz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 383 * ccomps * dcomps);

            auto g_yz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 389 * ccomps * dcomps);

            auto g_yz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 395 * ccomps * dcomps);

            auto g_yz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 419 * ccomps * dcomps);

            auto g_yz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 449 * ccomps * dcomps);

            auto g_zz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 455 * ccomps * dcomps);

            auto g_zz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 461 * ccomps * dcomps);

            auto g_zz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 464 * ccomps * dcomps);

            auto g_zz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 467 * ccomps * dcomps);

            auto g_zz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 473 * ccomps * dcomps);

            auto g_zz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 479 * ccomps * dcomps);

            auto g_zz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 485 * ccomps * dcomps);

            auto g_zz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 491 * ccomps * dcomps);

            auto g_zz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 494 * ccomps * dcomps);

            auto g_zz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 497 * ccomps * dcomps);

            auto g_zz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 509 * ccomps * dcomps);

            auto g_zz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 539 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hpxx

            const auto hp_geom_20_off = idx_geom_20_hpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_x, g_x_0_xxxx_y, g_x_0_xxxx_z, g_xx_0_xxxx_x, g_xx_0_xxxx_xx, g_xx_0_xxxx_xy, g_xx_0_xxxx_xz, g_xx_0_xxxx_y, g_xx_0_xxxx_z, g_xx_0_xxxxx_x, g_xx_0_xxxxx_y, g_xx_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxx_x[k] = -2.0 * g_x_0_xxxx_x[k] - g_xx_0_xxxx_x[k] * ab_x + g_xx_0_xxxx_xx[k];

                g_xx_0_xxxxx_y[k] = -2.0 * g_x_0_xxxx_y[k] - g_xx_0_xxxx_y[k] * ab_x + g_xx_0_xxxx_xy[k];

                g_xx_0_xxxxx_z[k] = -2.0 * g_x_0_xxxx_z[k] - g_xx_0_xxxx_z[k] * ab_x + g_xx_0_xxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_x, g_xx_0_xxxx_xy, g_xx_0_xxxx_y, g_xx_0_xxxx_yy, g_xx_0_xxxx_yz, g_xx_0_xxxx_z, g_xx_0_xxxxy_x, g_xx_0_xxxxy_y, g_xx_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxy_x[k] = -g_xx_0_xxxx_x[k] * ab_y + g_xx_0_xxxx_xy[k];

                g_xx_0_xxxxy_y[k] = -g_xx_0_xxxx_y[k] * ab_y + g_xx_0_xxxx_yy[k];

                g_xx_0_xxxxy_z[k] = -g_xx_0_xxxx_z[k] * ab_y + g_xx_0_xxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_x, g_xx_0_xxxx_xz, g_xx_0_xxxx_y, g_xx_0_xxxx_yz, g_xx_0_xxxx_z, g_xx_0_xxxx_zz, g_xx_0_xxxxz_x, g_xx_0_xxxxz_y, g_xx_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxz_x[k] = -g_xx_0_xxxx_x[k] * ab_z + g_xx_0_xxxx_xz[k];

                g_xx_0_xxxxz_y[k] = -g_xx_0_xxxx_y[k] * ab_z + g_xx_0_xxxx_yz[k];

                g_xx_0_xxxxz_z[k] = -g_xx_0_xxxx_z[k] * ab_z + g_xx_0_xxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxy_x, g_xx_0_xxxy_xy, g_xx_0_xxxy_y, g_xx_0_xxxy_yy, g_xx_0_xxxy_yz, g_xx_0_xxxy_z, g_xx_0_xxxyy_x, g_xx_0_xxxyy_y, g_xx_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyy_x[k] = -g_xx_0_xxxy_x[k] * ab_y + g_xx_0_xxxy_xy[k];

                g_xx_0_xxxyy_y[k] = -g_xx_0_xxxy_y[k] * ab_y + g_xx_0_xxxy_yy[k];

                g_xx_0_xxxyy_z[k] = -g_xx_0_xxxy_z[k] * ab_y + g_xx_0_xxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyz_x, g_xx_0_xxxyz_y, g_xx_0_xxxyz_z, g_xx_0_xxxz_x, g_xx_0_xxxz_xy, g_xx_0_xxxz_y, g_xx_0_xxxz_yy, g_xx_0_xxxz_yz, g_xx_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyz_x[k] = -g_xx_0_xxxz_x[k] * ab_y + g_xx_0_xxxz_xy[k];

                g_xx_0_xxxyz_y[k] = -g_xx_0_xxxz_y[k] * ab_y + g_xx_0_xxxz_yy[k];

                g_xx_0_xxxyz_z[k] = -g_xx_0_xxxz_z[k] * ab_y + g_xx_0_xxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxz_x, g_xx_0_xxxz_xz, g_xx_0_xxxz_y, g_xx_0_xxxz_yz, g_xx_0_xxxz_z, g_xx_0_xxxz_zz, g_xx_0_xxxzz_x, g_xx_0_xxxzz_y, g_xx_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzz_x[k] = -g_xx_0_xxxz_x[k] * ab_z + g_xx_0_xxxz_xz[k];

                g_xx_0_xxxzz_y[k] = -g_xx_0_xxxz_y[k] * ab_z + g_xx_0_xxxz_yz[k];

                g_xx_0_xxxzz_z[k] = -g_xx_0_xxxz_z[k] * ab_z + g_xx_0_xxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyy_x, g_xx_0_xxyy_xy, g_xx_0_xxyy_y, g_xx_0_xxyy_yy, g_xx_0_xxyy_yz, g_xx_0_xxyy_z, g_xx_0_xxyyy_x, g_xx_0_xxyyy_y, g_xx_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyy_x[k] = -g_xx_0_xxyy_x[k] * ab_y + g_xx_0_xxyy_xy[k];

                g_xx_0_xxyyy_y[k] = -g_xx_0_xxyy_y[k] * ab_y + g_xx_0_xxyy_yy[k];

                g_xx_0_xxyyy_z[k] = -g_xx_0_xxyy_z[k] * ab_y + g_xx_0_xxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyz_x, g_xx_0_xxyyz_y, g_xx_0_xxyyz_z, g_xx_0_xxyz_x, g_xx_0_xxyz_xy, g_xx_0_xxyz_y, g_xx_0_xxyz_yy, g_xx_0_xxyz_yz, g_xx_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyz_x[k] = -g_xx_0_xxyz_x[k] * ab_y + g_xx_0_xxyz_xy[k];

                g_xx_0_xxyyz_y[k] = -g_xx_0_xxyz_y[k] * ab_y + g_xx_0_xxyz_yy[k];

                g_xx_0_xxyyz_z[k] = -g_xx_0_xxyz_z[k] * ab_y + g_xx_0_xxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzz_x, g_xx_0_xxyzz_y, g_xx_0_xxyzz_z, g_xx_0_xxzz_x, g_xx_0_xxzz_xy, g_xx_0_xxzz_y, g_xx_0_xxzz_yy, g_xx_0_xxzz_yz, g_xx_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzz_x[k] = -g_xx_0_xxzz_x[k] * ab_y + g_xx_0_xxzz_xy[k];

                g_xx_0_xxyzz_y[k] = -g_xx_0_xxzz_y[k] * ab_y + g_xx_0_xxzz_yy[k];

                g_xx_0_xxyzz_z[k] = -g_xx_0_xxzz_z[k] * ab_y + g_xx_0_xxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzz_x, g_xx_0_xxzz_xz, g_xx_0_xxzz_y, g_xx_0_xxzz_yz, g_xx_0_xxzz_z, g_xx_0_xxzz_zz, g_xx_0_xxzzz_x, g_xx_0_xxzzz_y, g_xx_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzz_x[k] = -g_xx_0_xxzz_x[k] * ab_z + g_xx_0_xxzz_xz[k];

                g_xx_0_xxzzz_y[k] = -g_xx_0_xxzz_y[k] * ab_z + g_xx_0_xxzz_yz[k];

                g_xx_0_xxzzz_z[k] = -g_xx_0_xxzz_z[k] * ab_z + g_xx_0_xxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyy_x, g_xx_0_xyyy_xy, g_xx_0_xyyy_y, g_xx_0_xyyy_yy, g_xx_0_xyyy_yz, g_xx_0_xyyy_z, g_xx_0_xyyyy_x, g_xx_0_xyyyy_y, g_xx_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyy_x[k] = -g_xx_0_xyyy_x[k] * ab_y + g_xx_0_xyyy_xy[k];

                g_xx_0_xyyyy_y[k] = -g_xx_0_xyyy_y[k] * ab_y + g_xx_0_xyyy_yy[k];

                g_xx_0_xyyyy_z[k] = -g_xx_0_xyyy_z[k] * ab_y + g_xx_0_xyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyz_x, g_xx_0_xyyyz_y, g_xx_0_xyyyz_z, g_xx_0_xyyz_x, g_xx_0_xyyz_xy, g_xx_0_xyyz_y, g_xx_0_xyyz_yy, g_xx_0_xyyz_yz, g_xx_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyz_x[k] = -g_xx_0_xyyz_x[k] * ab_y + g_xx_0_xyyz_xy[k];

                g_xx_0_xyyyz_y[k] = -g_xx_0_xyyz_y[k] * ab_y + g_xx_0_xyyz_yy[k];

                g_xx_0_xyyyz_z[k] = -g_xx_0_xyyz_z[k] * ab_y + g_xx_0_xyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzz_x, g_xx_0_xyyzz_y, g_xx_0_xyyzz_z, g_xx_0_xyzz_x, g_xx_0_xyzz_xy, g_xx_0_xyzz_y, g_xx_0_xyzz_yy, g_xx_0_xyzz_yz, g_xx_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzz_x[k] = -g_xx_0_xyzz_x[k] * ab_y + g_xx_0_xyzz_xy[k];

                g_xx_0_xyyzz_y[k] = -g_xx_0_xyzz_y[k] * ab_y + g_xx_0_xyzz_yy[k];

                g_xx_0_xyyzz_z[k] = -g_xx_0_xyzz_z[k] * ab_y + g_xx_0_xyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzz_x, g_xx_0_xyzzz_y, g_xx_0_xyzzz_z, g_xx_0_xzzz_x, g_xx_0_xzzz_xy, g_xx_0_xzzz_y, g_xx_0_xzzz_yy, g_xx_0_xzzz_yz, g_xx_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzz_x[k] = -g_xx_0_xzzz_x[k] * ab_y + g_xx_0_xzzz_xy[k];

                g_xx_0_xyzzz_y[k] = -g_xx_0_xzzz_y[k] * ab_y + g_xx_0_xzzz_yy[k];

                g_xx_0_xyzzz_z[k] = -g_xx_0_xzzz_z[k] * ab_y + g_xx_0_xzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzz_x, g_xx_0_xzzz_xz, g_xx_0_xzzz_y, g_xx_0_xzzz_yz, g_xx_0_xzzz_z, g_xx_0_xzzz_zz, g_xx_0_xzzzz_x, g_xx_0_xzzzz_y, g_xx_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzz_x[k] = -g_xx_0_xzzz_x[k] * ab_z + g_xx_0_xzzz_xz[k];

                g_xx_0_xzzzz_y[k] = -g_xx_0_xzzz_y[k] * ab_z + g_xx_0_xzzz_yz[k];

                g_xx_0_xzzzz_z[k] = -g_xx_0_xzzz_z[k] * ab_z + g_xx_0_xzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyy_x, g_xx_0_yyyy_xy, g_xx_0_yyyy_y, g_xx_0_yyyy_yy, g_xx_0_yyyy_yz, g_xx_0_yyyy_z, g_xx_0_yyyyy_x, g_xx_0_yyyyy_y, g_xx_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyy_x[k] = -g_xx_0_yyyy_x[k] * ab_y + g_xx_0_yyyy_xy[k];

                g_xx_0_yyyyy_y[k] = -g_xx_0_yyyy_y[k] * ab_y + g_xx_0_yyyy_yy[k];

                g_xx_0_yyyyy_z[k] = -g_xx_0_yyyy_z[k] * ab_y + g_xx_0_yyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyz_x, g_xx_0_yyyyz_y, g_xx_0_yyyyz_z, g_xx_0_yyyz_x, g_xx_0_yyyz_xy, g_xx_0_yyyz_y, g_xx_0_yyyz_yy, g_xx_0_yyyz_yz, g_xx_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyz_x[k] = -g_xx_0_yyyz_x[k] * ab_y + g_xx_0_yyyz_xy[k];

                g_xx_0_yyyyz_y[k] = -g_xx_0_yyyz_y[k] * ab_y + g_xx_0_yyyz_yy[k];

                g_xx_0_yyyyz_z[k] = -g_xx_0_yyyz_z[k] * ab_y + g_xx_0_yyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzz_x, g_xx_0_yyyzz_y, g_xx_0_yyyzz_z, g_xx_0_yyzz_x, g_xx_0_yyzz_xy, g_xx_0_yyzz_y, g_xx_0_yyzz_yy, g_xx_0_yyzz_yz, g_xx_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzz_x[k] = -g_xx_0_yyzz_x[k] * ab_y + g_xx_0_yyzz_xy[k];

                g_xx_0_yyyzz_y[k] = -g_xx_0_yyzz_y[k] * ab_y + g_xx_0_yyzz_yy[k];

                g_xx_0_yyyzz_z[k] = -g_xx_0_yyzz_z[k] * ab_y + g_xx_0_yyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzz_x, g_xx_0_yyzzz_y, g_xx_0_yyzzz_z, g_xx_0_yzzz_x, g_xx_0_yzzz_xy, g_xx_0_yzzz_y, g_xx_0_yzzz_yy, g_xx_0_yzzz_yz, g_xx_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzz_x[k] = -g_xx_0_yzzz_x[k] * ab_y + g_xx_0_yzzz_xy[k];

                g_xx_0_yyzzz_y[k] = -g_xx_0_yzzz_y[k] * ab_y + g_xx_0_yzzz_yy[k];

                g_xx_0_yyzzz_z[k] = -g_xx_0_yzzz_z[k] * ab_y + g_xx_0_yzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzz_x, g_xx_0_yzzzz_y, g_xx_0_yzzzz_z, g_xx_0_zzzz_x, g_xx_0_zzzz_xy, g_xx_0_zzzz_y, g_xx_0_zzzz_yy, g_xx_0_zzzz_yz, g_xx_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzz_x[k] = -g_xx_0_zzzz_x[k] * ab_y + g_xx_0_zzzz_xy[k];

                g_xx_0_yzzzz_y[k] = -g_xx_0_zzzz_y[k] * ab_y + g_xx_0_zzzz_yy[k];

                g_xx_0_yzzzz_z[k] = -g_xx_0_zzzz_z[k] * ab_y + g_xx_0_zzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzz_x, g_xx_0_zzzz_xz, g_xx_0_zzzz_y, g_xx_0_zzzz_yz, g_xx_0_zzzz_z, g_xx_0_zzzz_zz, g_xx_0_zzzzz_x, g_xx_0_zzzzz_y, g_xx_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzz_x[k] = -g_xx_0_zzzz_x[k] * ab_z + g_xx_0_zzzz_xz[k];

                g_xx_0_zzzzz_y[k] = -g_xx_0_zzzz_y[k] * ab_z + g_xx_0_zzzz_yz[k];

                g_xx_0_zzzzz_z[k] = -g_xx_0_zzzz_z[k] * ab_z + g_xx_0_zzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_x, g_xy_0_xxxx_xx, g_xy_0_xxxx_xy, g_xy_0_xxxx_xz, g_xy_0_xxxx_y, g_xy_0_xxxx_z, g_xy_0_xxxxx_x, g_xy_0_xxxxx_y, g_xy_0_xxxxx_z, g_y_0_xxxx_x, g_y_0_xxxx_y, g_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxx_x[k] = -g_y_0_xxxx_x[k] - g_xy_0_xxxx_x[k] * ab_x + g_xy_0_xxxx_xx[k];

                g_xy_0_xxxxx_y[k] = -g_y_0_xxxx_y[k] - g_xy_0_xxxx_y[k] * ab_x + g_xy_0_xxxx_xy[k];

                g_xy_0_xxxxx_z[k] = -g_y_0_xxxx_z[k] - g_xy_0_xxxx_z[k] * ab_x + g_xy_0_xxxx_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_x, g_xy_0_xxxxy_y, g_xy_0_xxxxy_z, g_xy_0_xxxy_x, g_xy_0_xxxy_xx, g_xy_0_xxxy_xy, g_xy_0_xxxy_xz, g_xy_0_xxxy_y, g_xy_0_xxxy_z, g_y_0_xxxy_x, g_y_0_xxxy_y, g_y_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxy_x[k] = -g_y_0_xxxy_x[k] - g_xy_0_xxxy_x[k] * ab_x + g_xy_0_xxxy_xx[k];

                g_xy_0_xxxxy_y[k] = -g_y_0_xxxy_y[k] - g_xy_0_xxxy_y[k] * ab_x + g_xy_0_xxxy_xy[k];

                g_xy_0_xxxxy_z[k] = -g_y_0_xxxy_z[k] - g_xy_0_xxxy_z[k] * ab_x + g_xy_0_xxxy_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_x, g_xy_0_xxxx_xz, g_xy_0_xxxx_y, g_xy_0_xxxx_yz, g_xy_0_xxxx_z, g_xy_0_xxxx_zz, g_xy_0_xxxxz_x, g_xy_0_xxxxz_y, g_xy_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxz_x[k] = -g_xy_0_xxxx_x[k] * ab_z + g_xy_0_xxxx_xz[k];

                g_xy_0_xxxxz_y[k] = -g_xy_0_xxxx_y[k] * ab_z + g_xy_0_xxxx_yz[k];

                g_xy_0_xxxxz_z[k] = -g_xy_0_xxxx_z[k] * ab_z + g_xy_0_xxxx_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_x, g_xy_0_xxxyy_y, g_xy_0_xxxyy_z, g_xy_0_xxyy_x, g_xy_0_xxyy_xx, g_xy_0_xxyy_xy, g_xy_0_xxyy_xz, g_xy_0_xxyy_y, g_xy_0_xxyy_z, g_y_0_xxyy_x, g_y_0_xxyy_y, g_y_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyy_x[k] = -g_y_0_xxyy_x[k] - g_xy_0_xxyy_x[k] * ab_x + g_xy_0_xxyy_xx[k];

                g_xy_0_xxxyy_y[k] = -g_y_0_xxyy_y[k] - g_xy_0_xxyy_y[k] * ab_x + g_xy_0_xxyy_xy[k];

                g_xy_0_xxxyy_z[k] = -g_y_0_xxyy_z[k] - g_xy_0_xxyy_z[k] * ab_x + g_xy_0_xxyy_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_x, g_xy_0_xxxy_xz, g_xy_0_xxxy_y, g_xy_0_xxxy_yz, g_xy_0_xxxy_z, g_xy_0_xxxy_zz, g_xy_0_xxxyz_x, g_xy_0_xxxyz_y, g_xy_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyz_x[k] = -g_xy_0_xxxy_x[k] * ab_z + g_xy_0_xxxy_xz[k];

                g_xy_0_xxxyz_y[k] = -g_xy_0_xxxy_y[k] * ab_z + g_xy_0_xxxy_yz[k];

                g_xy_0_xxxyz_z[k] = -g_xy_0_xxxy_z[k] * ab_z + g_xy_0_xxxy_zz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxz_x, g_xy_0_xxxz_xz, g_xy_0_xxxz_y, g_xy_0_xxxz_yz, g_xy_0_xxxz_z, g_xy_0_xxxz_zz, g_xy_0_xxxzz_x, g_xy_0_xxxzz_y, g_xy_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzz_x[k] = -g_xy_0_xxxz_x[k] * ab_z + g_xy_0_xxxz_xz[k];

                g_xy_0_xxxzz_y[k] = -g_xy_0_xxxz_y[k] * ab_z + g_xy_0_xxxz_yz[k];

                g_xy_0_xxxzz_z[k] = -g_xy_0_xxxz_z[k] * ab_z + g_xy_0_xxxz_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_x, g_xy_0_xxyyy_y, g_xy_0_xxyyy_z, g_xy_0_xyyy_x, g_xy_0_xyyy_xx, g_xy_0_xyyy_xy, g_xy_0_xyyy_xz, g_xy_0_xyyy_y, g_xy_0_xyyy_z, g_y_0_xyyy_x, g_y_0_xyyy_y, g_y_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyy_x[k] = -g_y_0_xyyy_x[k] - g_xy_0_xyyy_x[k] * ab_x + g_xy_0_xyyy_xx[k];

                g_xy_0_xxyyy_y[k] = -g_y_0_xyyy_y[k] - g_xy_0_xyyy_y[k] * ab_x + g_xy_0_xyyy_xy[k];

                g_xy_0_xxyyy_z[k] = -g_y_0_xyyy_z[k] - g_xy_0_xyyy_z[k] * ab_x + g_xy_0_xyyy_xz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_x, g_xy_0_xxyy_xz, g_xy_0_xxyy_y, g_xy_0_xxyy_yz, g_xy_0_xxyy_z, g_xy_0_xxyy_zz, g_xy_0_xxyyz_x, g_xy_0_xxyyz_y, g_xy_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyz_x[k] = -g_xy_0_xxyy_x[k] * ab_z + g_xy_0_xxyy_xz[k];

                g_xy_0_xxyyz_y[k] = -g_xy_0_xxyy_y[k] * ab_z + g_xy_0_xxyy_yz[k];

                g_xy_0_xxyyz_z[k] = -g_xy_0_xxyy_z[k] * ab_z + g_xy_0_xxyy_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyz_x, g_xy_0_xxyz_xz, g_xy_0_xxyz_y, g_xy_0_xxyz_yz, g_xy_0_xxyz_z, g_xy_0_xxyz_zz, g_xy_0_xxyzz_x, g_xy_0_xxyzz_y, g_xy_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzz_x[k] = -g_xy_0_xxyz_x[k] * ab_z + g_xy_0_xxyz_xz[k];

                g_xy_0_xxyzz_y[k] = -g_xy_0_xxyz_y[k] * ab_z + g_xy_0_xxyz_yz[k];

                g_xy_0_xxyzz_z[k] = -g_xy_0_xxyz_z[k] * ab_z + g_xy_0_xxyz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzz_x, g_xy_0_xxzz_xz, g_xy_0_xxzz_y, g_xy_0_xxzz_yz, g_xy_0_xxzz_z, g_xy_0_xxzz_zz, g_xy_0_xxzzz_x, g_xy_0_xxzzz_y, g_xy_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzz_x[k] = -g_xy_0_xxzz_x[k] * ab_z + g_xy_0_xxzz_xz[k];

                g_xy_0_xxzzz_y[k] = -g_xy_0_xxzz_y[k] * ab_z + g_xy_0_xxzz_yz[k];

                g_xy_0_xxzzz_z[k] = -g_xy_0_xxzz_z[k] * ab_z + g_xy_0_xxzz_zz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_x, g_xy_0_xyyyy_y, g_xy_0_xyyyy_z, g_xy_0_yyyy_x, g_xy_0_yyyy_xx, g_xy_0_yyyy_xy, g_xy_0_yyyy_xz, g_xy_0_yyyy_y, g_xy_0_yyyy_z, g_y_0_yyyy_x, g_y_0_yyyy_y, g_y_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyy_x[k] = -g_y_0_yyyy_x[k] - g_xy_0_yyyy_x[k] * ab_x + g_xy_0_yyyy_xx[k];

                g_xy_0_xyyyy_y[k] = -g_y_0_yyyy_y[k] - g_xy_0_yyyy_y[k] * ab_x + g_xy_0_yyyy_xy[k];

                g_xy_0_xyyyy_z[k] = -g_y_0_yyyy_z[k] - g_xy_0_yyyy_z[k] * ab_x + g_xy_0_yyyy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_x, g_xy_0_xyyy_xz, g_xy_0_xyyy_y, g_xy_0_xyyy_yz, g_xy_0_xyyy_z, g_xy_0_xyyy_zz, g_xy_0_xyyyz_x, g_xy_0_xyyyz_y, g_xy_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyz_x[k] = -g_xy_0_xyyy_x[k] * ab_z + g_xy_0_xyyy_xz[k];

                g_xy_0_xyyyz_y[k] = -g_xy_0_xyyy_y[k] * ab_z + g_xy_0_xyyy_yz[k];

                g_xy_0_xyyyz_z[k] = -g_xy_0_xyyy_z[k] * ab_z + g_xy_0_xyyy_zz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyz_x, g_xy_0_xyyz_xz, g_xy_0_xyyz_y, g_xy_0_xyyz_yz, g_xy_0_xyyz_z, g_xy_0_xyyz_zz, g_xy_0_xyyzz_x, g_xy_0_xyyzz_y, g_xy_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzz_x[k] = -g_xy_0_xyyz_x[k] * ab_z + g_xy_0_xyyz_xz[k];

                g_xy_0_xyyzz_y[k] = -g_xy_0_xyyz_y[k] * ab_z + g_xy_0_xyyz_yz[k];

                g_xy_0_xyyzz_z[k] = -g_xy_0_xyyz_z[k] * ab_z + g_xy_0_xyyz_zz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzz_x, g_xy_0_xyzz_xz, g_xy_0_xyzz_y, g_xy_0_xyzz_yz, g_xy_0_xyzz_z, g_xy_0_xyzz_zz, g_xy_0_xyzzz_x, g_xy_0_xyzzz_y, g_xy_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzz_x[k] = -g_xy_0_xyzz_x[k] * ab_z + g_xy_0_xyzz_xz[k];

                g_xy_0_xyzzz_y[k] = -g_xy_0_xyzz_y[k] * ab_z + g_xy_0_xyzz_yz[k];

                g_xy_0_xyzzz_z[k] = -g_xy_0_xyzz_z[k] * ab_z + g_xy_0_xyzz_zz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzz_x, g_xy_0_xzzz_xz, g_xy_0_xzzz_y, g_xy_0_xzzz_yz, g_xy_0_xzzz_z, g_xy_0_xzzz_zz, g_xy_0_xzzzz_x, g_xy_0_xzzzz_y, g_xy_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzz_x[k] = -g_xy_0_xzzz_x[k] * ab_z + g_xy_0_xzzz_xz[k];

                g_xy_0_xzzzz_y[k] = -g_xy_0_xzzz_y[k] * ab_z + g_xy_0_xzzz_yz[k];

                g_xy_0_xzzzz_z[k] = -g_xy_0_xzzz_z[k] * ab_z + g_xy_0_xzzz_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_x, g_x_0_yyyy_y, g_x_0_yyyy_z, g_xy_0_yyyy_x, g_xy_0_yyyy_xy, g_xy_0_yyyy_y, g_xy_0_yyyy_yy, g_xy_0_yyyy_yz, g_xy_0_yyyy_z, g_xy_0_yyyyy_x, g_xy_0_yyyyy_y, g_xy_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyy_x[k] = -g_x_0_yyyy_x[k] - g_xy_0_yyyy_x[k] * ab_y + g_xy_0_yyyy_xy[k];

                g_xy_0_yyyyy_y[k] = -g_x_0_yyyy_y[k] - g_xy_0_yyyy_y[k] * ab_y + g_xy_0_yyyy_yy[k];

                g_xy_0_yyyyy_z[k] = -g_x_0_yyyy_z[k] - g_xy_0_yyyy_z[k] * ab_y + g_xy_0_yyyy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyy_x, g_xy_0_yyyy_xz, g_xy_0_yyyy_y, g_xy_0_yyyy_yz, g_xy_0_yyyy_z, g_xy_0_yyyy_zz, g_xy_0_yyyyz_x, g_xy_0_yyyyz_y, g_xy_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyz_x[k] = -g_xy_0_yyyy_x[k] * ab_z + g_xy_0_yyyy_xz[k];

                g_xy_0_yyyyz_y[k] = -g_xy_0_yyyy_y[k] * ab_z + g_xy_0_yyyy_yz[k];

                g_xy_0_yyyyz_z[k] = -g_xy_0_yyyy_z[k] * ab_z + g_xy_0_yyyy_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyz_x, g_xy_0_yyyz_xz, g_xy_0_yyyz_y, g_xy_0_yyyz_yz, g_xy_0_yyyz_z, g_xy_0_yyyz_zz, g_xy_0_yyyzz_x, g_xy_0_yyyzz_y, g_xy_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzz_x[k] = -g_xy_0_yyyz_x[k] * ab_z + g_xy_0_yyyz_xz[k];

                g_xy_0_yyyzz_y[k] = -g_xy_0_yyyz_y[k] * ab_z + g_xy_0_yyyz_yz[k];

                g_xy_0_yyyzz_z[k] = -g_xy_0_yyyz_z[k] * ab_z + g_xy_0_yyyz_zz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzz_x, g_xy_0_yyzz_xz, g_xy_0_yyzz_y, g_xy_0_yyzz_yz, g_xy_0_yyzz_z, g_xy_0_yyzz_zz, g_xy_0_yyzzz_x, g_xy_0_yyzzz_y, g_xy_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzz_x[k] = -g_xy_0_yyzz_x[k] * ab_z + g_xy_0_yyzz_xz[k];

                g_xy_0_yyzzz_y[k] = -g_xy_0_yyzz_y[k] * ab_z + g_xy_0_yyzz_yz[k];

                g_xy_0_yyzzz_z[k] = -g_xy_0_yyzz_z[k] * ab_z + g_xy_0_yyzz_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzz_x, g_xy_0_yzzz_xz, g_xy_0_yzzz_y, g_xy_0_yzzz_yz, g_xy_0_yzzz_z, g_xy_0_yzzz_zz, g_xy_0_yzzzz_x, g_xy_0_yzzzz_y, g_xy_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzz_x[k] = -g_xy_0_yzzz_x[k] * ab_z + g_xy_0_yzzz_xz[k];

                g_xy_0_yzzzz_y[k] = -g_xy_0_yzzz_y[k] * ab_z + g_xy_0_yzzz_yz[k];

                g_xy_0_yzzzz_z[k] = -g_xy_0_yzzz_z[k] * ab_z + g_xy_0_yzzz_zz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzz_x, g_xy_0_zzzz_xz, g_xy_0_zzzz_y, g_xy_0_zzzz_yz, g_xy_0_zzzz_z, g_xy_0_zzzz_zz, g_xy_0_zzzzz_x, g_xy_0_zzzzz_y, g_xy_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzz_x[k] = -g_xy_0_zzzz_x[k] * ab_z + g_xy_0_zzzz_xz[k];

                g_xy_0_zzzzz_y[k] = -g_xy_0_zzzz_y[k] * ab_z + g_xy_0_zzzz_yz[k];

                g_xy_0_zzzzz_z[k] = -g_xy_0_zzzz_z[k] * ab_z + g_xy_0_zzzz_zz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_x, g_xz_0_xxxx_xx, g_xz_0_xxxx_xy, g_xz_0_xxxx_xz, g_xz_0_xxxx_y, g_xz_0_xxxx_z, g_xz_0_xxxxx_x, g_xz_0_xxxxx_y, g_xz_0_xxxxx_z, g_z_0_xxxx_x, g_z_0_xxxx_y, g_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxx_x[k] = -g_z_0_xxxx_x[k] - g_xz_0_xxxx_x[k] * ab_x + g_xz_0_xxxx_xx[k];

                g_xz_0_xxxxx_y[k] = -g_z_0_xxxx_y[k] - g_xz_0_xxxx_y[k] * ab_x + g_xz_0_xxxx_xy[k];

                g_xz_0_xxxxx_z[k] = -g_z_0_xxxx_z[k] - g_xz_0_xxxx_z[k] * ab_x + g_xz_0_xxxx_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_x, g_xz_0_xxxx_xy, g_xz_0_xxxx_y, g_xz_0_xxxx_yy, g_xz_0_xxxx_yz, g_xz_0_xxxx_z, g_xz_0_xxxxy_x, g_xz_0_xxxxy_y, g_xz_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxy_x[k] = -g_xz_0_xxxx_x[k] * ab_y + g_xz_0_xxxx_xy[k];

                g_xz_0_xxxxy_y[k] = -g_xz_0_xxxx_y[k] * ab_y + g_xz_0_xxxx_yy[k];

                g_xz_0_xxxxy_z[k] = -g_xz_0_xxxx_z[k] * ab_y + g_xz_0_xxxx_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxz_x, g_xz_0_xxxxz_y, g_xz_0_xxxxz_z, g_xz_0_xxxz_x, g_xz_0_xxxz_xx, g_xz_0_xxxz_xy, g_xz_0_xxxz_xz, g_xz_0_xxxz_y, g_xz_0_xxxz_z, g_z_0_xxxz_x, g_z_0_xxxz_y, g_z_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxz_x[k] = -g_z_0_xxxz_x[k] - g_xz_0_xxxz_x[k] * ab_x + g_xz_0_xxxz_xx[k];

                g_xz_0_xxxxz_y[k] = -g_z_0_xxxz_y[k] - g_xz_0_xxxz_y[k] * ab_x + g_xz_0_xxxz_xy[k];

                g_xz_0_xxxxz_z[k] = -g_z_0_xxxz_z[k] - g_xz_0_xxxz_z[k] * ab_x + g_xz_0_xxxz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxy_x, g_xz_0_xxxy_xy, g_xz_0_xxxy_y, g_xz_0_xxxy_yy, g_xz_0_xxxy_yz, g_xz_0_xxxy_z, g_xz_0_xxxyy_x, g_xz_0_xxxyy_y, g_xz_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyy_x[k] = -g_xz_0_xxxy_x[k] * ab_y + g_xz_0_xxxy_xy[k];

                g_xz_0_xxxyy_y[k] = -g_xz_0_xxxy_y[k] * ab_y + g_xz_0_xxxy_yy[k];

                g_xz_0_xxxyy_z[k] = -g_xz_0_xxxy_z[k] * ab_y + g_xz_0_xxxy_yz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyz_x, g_xz_0_xxxyz_y, g_xz_0_xxxyz_z, g_xz_0_xxxz_x, g_xz_0_xxxz_xy, g_xz_0_xxxz_y, g_xz_0_xxxz_yy, g_xz_0_xxxz_yz, g_xz_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyz_x[k] = -g_xz_0_xxxz_x[k] * ab_y + g_xz_0_xxxz_xy[k];

                g_xz_0_xxxyz_y[k] = -g_xz_0_xxxz_y[k] * ab_y + g_xz_0_xxxz_yy[k];

                g_xz_0_xxxyz_z[k] = -g_xz_0_xxxz_z[k] * ab_y + g_xz_0_xxxz_yz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzz_x, g_xz_0_xxxzz_y, g_xz_0_xxxzz_z, g_xz_0_xxzz_x, g_xz_0_xxzz_xx, g_xz_0_xxzz_xy, g_xz_0_xxzz_xz, g_xz_0_xxzz_y, g_xz_0_xxzz_z, g_z_0_xxzz_x, g_z_0_xxzz_y, g_z_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzz_x[k] = -g_z_0_xxzz_x[k] - g_xz_0_xxzz_x[k] * ab_x + g_xz_0_xxzz_xx[k];

                g_xz_0_xxxzz_y[k] = -g_z_0_xxzz_y[k] - g_xz_0_xxzz_y[k] * ab_x + g_xz_0_xxzz_xy[k];

                g_xz_0_xxxzz_z[k] = -g_z_0_xxzz_z[k] - g_xz_0_xxzz_z[k] * ab_x + g_xz_0_xxzz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyy_x, g_xz_0_xxyy_xy, g_xz_0_xxyy_y, g_xz_0_xxyy_yy, g_xz_0_xxyy_yz, g_xz_0_xxyy_z, g_xz_0_xxyyy_x, g_xz_0_xxyyy_y, g_xz_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyy_x[k] = -g_xz_0_xxyy_x[k] * ab_y + g_xz_0_xxyy_xy[k];

                g_xz_0_xxyyy_y[k] = -g_xz_0_xxyy_y[k] * ab_y + g_xz_0_xxyy_yy[k];

                g_xz_0_xxyyy_z[k] = -g_xz_0_xxyy_z[k] * ab_y + g_xz_0_xxyy_yz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyz_x, g_xz_0_xxyyz_y, g_xz_0_xxyyz_z, g_xz_0_xxyz_x, g_xz_0_xxyz_xy, g_xz_0_xxyz_y, g_xz_0_xxyz_yy, g_xz_0_xxyz_yz, g_xz_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyz_x[k] = -g_xz_0_xxyz_x[k] * ab_y + g_xz_0_xxyz_xy[k];

                g_xz_0_xxyyz_y[k] = -g_xz_0_xxyz_y[k] * ab_y + g_xz_0_xxyz_yy[k];

                g_xz_0_xxyyz_z[k] = -g_xz_0_xxyz_z[k] * ab_y + g_xz_0_xxyz_yz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzz_x, g_xz_0_xxyzz_y, g_xz_0_xxyzz_z, g_xz_0_xxzz_x, g_xz_0_xxzz_xy, g_xz_0_xxzz_y, g_xz_0_xxzz_yy, g_xz_0_xxzz_yz, g_xz_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzz_x[k] = -g_xz_0_xxzz_x[k] * ab_y + g_xz_0_xxzz_xy[k];

                g_xz_0_xxyzz_y[k] = -g_xz_0_xxzz_y[k] * ab_y + g_xz_0_xxzz_yy[k];

                g_xz_0_xxyzz_z[k] = -g_xz_0_xxzz_z[k] * ab_y + g_xz_0_xxzz_yz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzz_x, g_xz_0_xxzzz_y, g_xz_0_xxzzz_z, g_xz_0_xzzz_x, g_xz_0_xzzz_xx, g_xz_0_xzzz_xy, g_xz_0_xzzz_xz, g_xz_0_xzzz_y, g_xz_0_xzzz_z, g_z_0_xzzz_x, g_z_0_xzzz_y, g_z_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzz_x[k] = -g_z_0_xzzz_x[k] - g_xz_0_xzzz_x[k] * ab_x + g_xz_0_xzzz_xx[k];

                g_xz_0_xxzzz_y[k] = -g_z_0_xzzz_y[k] - g_xz_0_xzzz_y[k] * ab_x + g_xz_0_xzzz_xy[k];

                g_xz_0_xxzzz_z[k] = -g_z_0_xzzz_z[k] - g_xz_0_xzzz_z[k] * ab_x + g_xz_0_xzzz_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyy_x, g_xz_0_xyyy_xy, g_xz_0_xyyy_y, g_xz_0_xyyy_yy, g_xz_0_xyyy_yz, g_xz_0_xyyy_z, g_xz_0_xyyyy_x, g_xz_0_xyyyy_y, g_xz_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyy_x[k] = -g_xz_0_xyyy_x[k] * ab_y + g_xz_0_xyyy_xy[k];

                g_xz_0_xyyyy_y[k] = -g_xz_0_xyyy_y[k] * ab_y + g_xz_0_xyyy_yy[k];

                g_xz_0_xyyyy_z[k] = -g_xz_0_xyyy_z[k] * ab_y + g_xz_0_xyyy_yz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyz_x, g_xz_0_xyyyz_y, g_xz_0_xyyyz_z, g_xz_0_xyyz_x, g_xz_0_xyyz_xy, g_xz_0_xyyz_y, g_xz_0_xyyz_yy, g_xz_0_xyyz_yz, g_xz_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyz_x[k] = -g_xz_0_xyyz_x[k] * ab_y + g_xz_0_xyyz_xy[k];

                g_xz_0_xyyyz_y[k] = -g_xz_0_xyyz_y[k] * ab_y + g_xz_0_xyyz_yy[k];

                g_xz_0_xyyyz_z[k] = -g_xz_0_xyyz_z[k] * ab_y + g_xz_0_xyyz_yz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzz_x, g_xz_0_xyyzz_y, g_xz_0_xyyzz_z, g_xz_0_xyzz_x, g_xz_0_xyzz_xy, g_xz_0_xyzz_y, g_xz_0_xyzz_yy, g_xz_0_xyzz_yz, g_xz_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzz_x[k] = -g_xz_0_xyzz_x[k] * ab_y + g_xz_0_xyzz_xy[k];

                g_xz_0_xyyzz_y[k] = -g_xz_0_xyzz_y[k] * ab_y + g_xz_0_xyzz_yy[k];

                g_xz_0_xyyzz_z[k] = -g_xz_0_xyzz_z[k] * ab_y + g_xz_0_xyzz_yz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzz_x, g_xz_0_xyzzz_y, g_xz_0_xyzzz_z, g_xz_0_xzzz_x, g_xz_0_xzzz_xy, g_xz_0_xzzz_y, g_xz_0_xzzz_yy, g_xz_0_xzzz_yz, g_xz_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzz_x[k] = -g_xz_0_xzzz_x[k] * ab_y + g_xz_0_xzzz_xy[k];

                g_xz_0_xyzzz_y[k] = -g_xz_0_xzzz_y[k] * ab_y + g_xz_0_xzzz_yy[k];

                g_xz_0_xyzzz_z[k] = -g_xz_0_xzzz_z[k] * ab_y + g_xz_0_xzzz_yz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzz_x, g_xz_0_xzzzz_y, g_xz_0_xzzzz_z, g_xz_0_zzzz_x, g_xz_0_zzzz_xx, g_xz_0_zzzz_xy, g_xz_0_zzzz_xz, g_xz_0_zzzz_y, g_xz_0_zzzz_z, g_z_0_zzzz_x, g_z_0_zzzz_y, g_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzz_x[k] = -g_z_0_zzzz_x[k] - g_xz_0_zzzz_x[k] * ab_x + g_xz_0_zzzz_xx[k];

                g_xz_0_xzzzz_y[k] = -g_z_0_zzzz_y[k] - g_xz_0_zzzz_y[k] * ab_x + g_xz_0_zzzz_xy[k];

                g_xz_0_xzzzz_z[k] = -g_z_0_zzzz_z[k] - g_xz_0_zzzz_z[k] * ab_x + g_xz_0_zzzz_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyy_x, g_xz_0_yyyy_xy, g_xz_0_yyyy_y, g_xz_0_yyyy_yy, g_xz_0_yyyy_yz, g_xz_0_yyyy_z, g_xz_0_yyyyy_x, g_xz_0_yyyyy_y, g_xz_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyy_x[k] = -g_xz_0_yyyy_x[k] * ab_y + g_xz_0_yyyy_xy[k];

                g_xz_0_yyyyy_y[k] = -g_xz_0_yyyy_y[k] * ab_y + g_xz_0_yyyy_yy[k];

                g_xz_0_yyyyy_z[k] = -g_xz_0_yyyy_z[k] * ab_y + g_xz_0_yyyy_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyz_x, g_xz_0_yyyyz_y, g_xz_0_yyyyz_z, g_xz_0_yyyz_x, g_xz_0_yyyz_xy, g_xz_0_yyyz_y, g_xz_0_yyyz_yy, g_xz_0_yyyz_yz, g_xz_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyz_x[k] = -g_xz_0_yyyz_x[k] * ab_y + g_xz_0_yyyz_xy[k];

                g_xz_0_yyyyz_y[k] = -g_xz_0_yyyz_y[k] * ab_y + g_xz_0_yyyz_yy[k];

                g_xz_0_yyyyz_z[k] = -g_xz_0_yyyz_z[k] * ab_y + g_xz_0_yyyz_yz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzz_x, g_xz_0_yyyzz_y, g_xz_0_yyyzz_z, g_xz_0_yyzz_x, g_xz_0_yyzz_xy, g_xz_0_yyzz_y, g_xz_0_yyzz_yy, g_xz_0_yyzz_yz, g_xz_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzz_x[k] = -g_xz_0_yyzz_x[k] * ab_y + g_xz_0_yyzz_xy[k];

                g_xz_0_yyyzz_y[k] = -g_xz_0_yyzz_y[k] * ab_y + g_xz_0_yyzz_yy[k];

                g_xz_0_yyyzz_z[k] = -g_xz_0_yyzz_z[k] * ab_y + g_xz_0_yyzz_yz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzz_x, g_xz_0_yyzzz_y, g_xz_0_yyzzz_z, g_xz_0_yzzz_x, g_xz_0_yzzz_xy, g_xz_0_yzzz_y, g_xz_0_yzzz_yy, g_xz_0_yzzz_yz, g_xz_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzz_x[k] = -g_xz_0_yzzz_x[k] * ab_y + g_xz_0_yzzz_xy[k];

                g_xz_0_yyzzz_y[k] = -g_xz_0_yzzz_y[k] * ab_y + g_xz_0_yzzz_yy[k];

                g_xz_0_yyzzz_z[k] = -g_xz_0_yzzz_z[k] * ab_y + g_xz_0_yzzz_yz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzz_x, g_xz_0_yzzzz_y, g_xz_0_yzzzz_z, g_xz_0_zzzz_x, g_xz_0_zzzz_xy, g_xz_0_zzzz_y, g_xz_0_zzzz_yy, g_xz_0_zzzz_yz, g_xz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzz_x[k] = -g_xz_0_zzzz_x[k] * ab_y + g_xz_0_zzzz_xy[k];

                g_xz_0_yzzzz_y[k] = -g_xz_0_zzzz_y[k] * ab_y + g_xz_0_zzzz_yy[k];

                g_xz_0_yzzzz_z[k] = -g_xz_0_zzzz_z[k] * ab_y + g_xz_0_zzzz_yz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_x, g_x_0_zzzz_y, g_x_0_zzzz_z, g_xz_0_zzzz_x, g_xz_0_zzzz_xz, g_xz_0_zzzz_y, g_xz_0_zzzz_yz, g_xz_0_zzzz_z, g_xz_0_zzzz_zz, g_xz_0_zzzzz_x, g_xz_0_zzzzz_y, g_xz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzz_x[k] = -g_x_0_zzzz_x[k] - g_xz_0_zzzz_x[k] * ab_z + g_xz_0_zzzz_xz[k];

                g_xz_0_zzzzz_y[k] = -g_x_0_zzzz_y[k] - g_xz_0_zzzz_y[k] * ab_z + g_xz_0_zzzz_yz[k];

                g_xz_0_zzzzz_z[k] = -g_x_0_zzzz_z[k] - g_xz_0_zzzz_z[k] * ab_z + g_xz_0_zzzz_zz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxx_x, g_yy_0_xxxx_xx, g_yy_0_xxxx_xy, g_yy_0_xxxx_xz, g_yy_0_xxxx_y, g_yy_0_xxxx_z, g_yy_0_xxxxx_x, g_yy_0_xxxxx_y, g_yy_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxx_x[k] = -g_yy_0_xxxx_x[k] * ab_x + g_yy_0_xxxx_xx[k];

                g_yy_0_xxxxx_y[k] = -g_yy_0_xxxx_y[k] * ab_x + g_yy_0_xxxx_xy[k];

                g_yy_0_xxxxx_z[k] = -g_yy_0_xxxx_z[k] * ab_x + g_yy_0_xxxx_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxy_x, g_yy_0_xxxxy_y, g_yy_0_xxxxy_z, g_yy_0_xxxy_x, g_yy_0_xxxy_xx, g_yy_0_xxxy_xy, g_yy_0_xxxy_xz, g_yy_0_xxxy_y, g_yy_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxy_x[k] = -g_yy_0_xxxy_x[k] * ab_x + g_yy_0_xxxy_xx[k];

                g_yy_0_xxxxy_y[k] = -g_yy_0_xxxy_y[k] * ab_x + g_yy_0_xxxy_xy[k];

                g_yy_0_xxxxy_z[k] = -g_yy_0_xxxy_z[k] * ab_x + g_yy_0_xxxy_xz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxz_x, g_yy_0_xxxxz_y, g_yy_0_xxxxz_z, g_yy_0_xxxz_x, g_yy_0_xxxz_xx, g_yy_0_xxxz_xy, g_yy_0_xxxz_xz, g_yy_0_xxxz_y, g_yy_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxz_x[k] = -g_yy_0_xxxz_x[k] * ab_x + g_yy_0_xxxz_xx[k];

                g_yy_0_xxxxz_y[k] = -g_yy_0_xxxz_y[k] * ab_x + g_yy_0_xxxz_xy[k];

                g_yy_0_xxxxz_z[k] = -g_yy_0_xxxz_z[k] * ab_x + g_yy_0_xxxz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyy_x, g_yy_0_xxxyy_y, g_yy_0_xxxyy_z, g_yy_0_xxyy_x, g_yy_0_xxyy_xx, g_yy_0_xxyy_xy, g_yy_0_xxyy_xz, g_yy_0_xxyy_y, g_yy_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyy_x[k] = -g_yy_0_xxyy_x[k] * ab_x + g_yy_0_xxyy_xx[k];

                g_yy_0_xxxyy_y[k] = -g_yy_0_xxyy_y[k] * ab_x + g_yy_0_xxyy_xy[k];

                g_yy_0_xxxyy_z[k] = -g_yy_0_xxyy_z[k] * ab_x + g_yy_0_xxyy_xz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyz_x, g_yy_0_xxxyz_y, g_yy_0_xxxyz_z, g_yy_0_xxyz_x, g_yy_0_xxyz_xx, g_yy_0_xxyz_xy, g_yy_0_xxyz_xz, g_yy_0_xxyz_y, g_yy_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyz_x[k] = -g_yy_0_xxyz_x[k] * ab_x + g_yy_0_xxyz_xx[k];

                g_yy_0_xxxyz_y[k] = -g_yy_0_xxyz_y[k] * ab_x + g_yy_0_xxyz_xy[k];

                g_yy_0_xxxyz_z[k] = -g_yy_0_xxyz_z[k] * ab_x + g_yy_0_xxyz_xz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzz_x, g_yy_0_xxxzz_y, g_yy_0_xxxzz_z, g_yy_0_xxzz_x, g_yy_0_xxzz_xx, g_yy_0_xxzz_xy, g_yy_0_xxzz_xz, g_yy_0_xxzz_y, g_yy_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzz_x[k] = -g_yy_0_xxzz_x[k] * ab_x + g_yy_0_xxzz_xx[k];

                g_yy_0_xxxzz_y[k] = -g_yy_0_xxzz_y[k] * ab_x + g_yy_0_xxzz_xy[k];

                g_yy_0_xxxzz_z[k] = -g_yy_0_xxzz_z[k] * ab_x + g_yy_0_xxzz_xz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyy_x, g_yy_0_xxyyy_y, g_yy_0_xxyyy_z, g_yy_0_xyyy_x, g_yy_0_xyyy_xx, g_yy_0_xyyy_xy, g_yy_0_xyyy_xz, g_yy_0_xyyy_y, g_yy_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyy_x[k] = -g_yy_0_xyyy_x[k] * ab_x + g_yy_0_xyyy_xx[k];

                g_yy_0_xxyyy_y[k] = -g_yy_0_xyyy_y[k] * ab_x + g_yy_0_xyyy_xy[k];

                g_yy_0_xxyyy_z[k] = -g_yy_0_xyyy_z[k] * ab_x + g_yy_0_xyyy_xz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyz_x, g_yy_0_xxyyz_y, g_yy_0_xxyyz_z, g_yy_0_xyyz_x, g_yy_0_xyyz_xx, g_yy_0_xyyz_xy, g_yy_0_xyyz_xz, g_yy_0_xyyz_y, g_yy_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyz_x[k] = -g_yy_0_xyyz_x[k] * ab_x + g_yy_0_xyyz_xx[k];

                g_yy_0_xxyyz_y[k] = -g_yy_0_xyyz_y[k] * ab_x + g_yy_0_xyyz_xy[k];

                g_yy_0_xxyyz_z[k] = -g_yy_0_xyyz_z[k] * ab_x + g_yy_0_xyyz_xz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzz_x, g_yy_0_xxyzz_y, g_yy_0_xxyzz_z, g_yy_0_xyzz_x, g_yy_0_xyzz_xx, g_yy_0_xyzz_xy, g_yy_0_xyzz_xz, g_yy_0_xyzz_y, g_yy_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzz_x[k] = -g_yy_0_xyzz_x[k] * ab_x + g_yy_0_xyzz_xx[k];

                g_yy_0_xxyzz_y[k] = -g_yy_0_xyzz_y[k] * ab_x + g_yy_0_xyzz_xy[k];

                g_yy_0_xxyzz_z[k] = -g_yy_0_xyzz_z[k] * ab_x + g_yy_0_xyzz_xz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzz_x, g_yy_0_xxzzz_y, g_yy_0_xxzzz_z, g_yy_0_xzzz_x, g_yy_0_xzzz_xx, g_yy_0_xzzz_xy, g_yy_0_xzzz_xz, g_yy_0_xzzz_y, g_yy_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzz_x[k] = -g_yy_0_xzzz_x[k] * ab_x + g_yy_0_xzzz_xx[k];

                g_yy_0_xxzzz_y[k] = -g_yy_0_xzzz_y[k] * ab_x + g_yy_0_xzzz_xy[k];

                g_yy_0_xxzzz_z[k] = -g_yy_0_xzzz_z[k] * ab_x + g_yy_0_xzzz_xz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyy_x, g_yy_0_xyyyy_y, g_yy_0_xyyyy_z, g_yy_0_yyyy_x, g_yy_0_yyyy_xx, g_yy_0_yyyy_xy, g_yy_0_yyyy_xz, g_yy_0_yyyy_y, g_yy_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyy_x[k] = -g_yy_0_yyyy_x[k] * ab_x + g_yy_0_yyyy_xx[k];

                g_yy_0_xyyyy_y[k] = -g_yy_0_yyyy_y[k] * ab_x + g_yy_0_yyyy_xy[k];

                g_yy_0_xyyyy_z[k] = -g_yy_0_yyyy_z[k] * ab_x + g_yy_0_yyyy_xz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyz_x, g_yy_0_xyyyz_y, g_yy_0_xyyyz_z, g_yy_0_yyyz_x, g_yy_0_yyyz_xx, g_yy_0_yyyz_xy, g_yy_0_yyyz_xz, g_yy_0_yyyz_y, g_yy_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyz_x[k] = -g_yy_0_yyyz_x[k] * ab_x + g_yy_0_yyyz_xx[k];

                g_yy_0_xyyyz_y[k] = -g_yy_0_yyyz_y[k] * ab_x + g_yy_0_yyyz_xy[k];

                g_yy_0_xyyyz_z[k] = -g_yy_0_yyyz_z[k] * ab_x + g_yy_0_yyyz_xz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzz_x, g_yy_0_xyyzz_y, g_yy_0_xyyzz_z, g_yy_0_yyzz_x, g_yy_0_yyzz_xx, g_yy_0_yyzz_xy, g_yy_0_yyzz_xz, g_yy_0_yyzz_y, g_yy_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzz_x[k] = -g_yy_0_yyzz_x[k] * ab_x + g_yy_0_yyzz_xx[k];

                g_yy_0_xyyzz_y[k] = -g_yy_0_yyzz_y[k] * ab_x + g_yy_0_yyzz_xy[k];

                g_yy_0_xyyzz_z[k] = -g_yy_0_yyzz_z[k] * ab_x + g_yy_0_yyzz_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzz_x, g_yy_0_xyzzz_y, g_yy_0_xyzzz_z, g_yy_0_yzzz_x, g_yy_0_yzzz_xx, g_yy_0_yzzz_xy, g_yy_0_yzzz_xz, g_yy_0_yzzz_y, g_yy_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzz_x[k] = -g_yy_0_yzzz_x[k] * ab_x + g_yy_0_yzzz_xx[k];

                g_yy_0_xyzzz_y[k] = -g_yy_0_yzzz_y[k] * ab_x + g_yy_0_yzzz_xy[k];

                g_yy_0_xyzzz_z[k] = -g_yy_0_yzzz_z[k] * ab_x + g_yy_0_yzzz_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzz_x, g_yy_0_xzzzz_y, g_yy_0_xzzzz_z, g_yy_0_zzzz_x, g_yy_0_zzzz_xx, g_yy_0_zzzz_xy, g_yy_0_zzzz_xz, g_yy_0_zzzz_y, g_yy_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzz_x[k] = -g_yy_0_zzzz_x[k] * ab_x + g_yy_0_zzzz_xx[k];

                g_yy_0_xzzzz_y[k] = -g_yy_0_zzzz_y[k] * ab_x + g_yy_0_zzzz_xy[k];

                g_yy_0_xzzzz_z[k] = -g_yy_0_zzzz_z[k] * ab_x + g_yy_0_zzzz_xz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_x, g_y_0_yyyy_y, g_y_0_yyyy_z, g_yy_0_yyyy_x, g_yy_0_yyyy_xy, g_yy_0_yyyy_y, g_yy_0_yyyy_yy, g_yy_0_yyyy_yz, g_yy_0_yyyy_z, g_yy_0_yyyyy_x, g_yy_0_yyyyy_y, g_yy_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyy_x[k] = -2.0 * g_y_0_yyyy_x[k] - g_yy_0_yyyy_x[k] * ab_y + g_yy_0_yyyy_xy[k];

                g_yy_0_yyyyy_y[k] = -2.0 * g_y_0_yyyy_y[k] - g_yy_0_yyyy_y[k] * ab_y + g_yy_0_yyyy_yy[k];

                g_yy_0_yyyyy_z[k] = -2.0 * g_y_0_yyyy_z[k] - g_yy_0_yyyy_z[k] * ab_y + g_yy_0_yyyy_yz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyy_x, g_yy_0_yyyy_xz, g_yy_0_yyyy_y, g_yy_0_yyyy_yz, g_yy_0_yyyy_z, g_yy_0_yyyy_zz, g_yy_0_yyyyz_x, g_yy_0_yyyyz_y, g_yy_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyz_x[k] = -g_yy_0_yyyy_x[k] * ab_z + g_yy_0_yyyy_xz[k];

                g_yy_0_yyyyz_y[k] = -g_yy_0_yyyy_y[k] * ab_z + g_yy_0_yyyy_yz[k];

                g_yy_0_yyyyz_z[k] = -g_yy_0_yyyy_z[k] * ab_z + g_yy_0_yyyy_zz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 240 * ccomps * dcomps);

            auto g_yy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 241 * ccomps * dcomps);

            auto g_yy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyz_x, g_yy_0_yyyz_xz, g_yy_0_yyyz_y, g_yy_0_yyyz_yz, g_yy_0_yyyz_z, g_yy_0_yyyz_zz, g_yy_0_yyyzz_x, g_yy_0_yyyzz_y, g_yy_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzz_x[k] = -g_yy_0_yyyz_x[k] * ab_z + g_yy_0_yyyz_xz[k];

                g_yy_0_yyyzz_y[k] = -g_yy_0_yyyz_y[k] * ab_z + g_yy_0_yyyz_yz[k];

                g_yy_0_yyyzz_z[k] = -g_yy_0_yyyz_z[k] * ab_z + g_yy_0_yyyz_zz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 243 * ccomps * dcomps);

            auto g_yy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 244 * ccomps * dcomps);

            auto g_yy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzz_x, g_yy_0_yyzz_xz, g_yy_0_yyzz_y, g_yy_0_yyzz_yz, g_yy_0_yyzz_z, g_yy_0_yyzz_zz, g_yy_0_yyzzz_x, g_yy_0_yyzzz_y, g_yy_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzz_x[k] = -g_yy_0_yyzz_x[k] * ab_z + g_yy_0_yyzz_xz[k];

                g_yy_0_yyzzz_y[k] = -g_yy_0_yyzz_y[k] * ab_z + g_yy_0_yyzz_yz[k];

                g_yy_0_yyzzz_z[k] = -g_yy_0_yyzz_z[k] * ab_z + g_yy_0_yyzz_zz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 246 * ccomps * dcomps);

            auto g_yy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 247 * ccomps * dcomps);

            auto g_yy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzz_x, g_yy_0_yzzz_xz, g_yy_0_yzzz_y, g_yy_0_yzzz_yz, g_yy_0_yzzz_z, g_yy_0_yzzz_zz, g_yy_0_yzzzz_x, g_yy_0_yzzzz_y, g_yy_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzz_x[k] = -g_yy_0_yzzz_x[k] * ab_z + g_yy_0_yzzz_xz[k];

                g_yy_0_yzzzz_y[k] = -g_yy_0_yzzz_y[k] * ab_z + g_yy_0_yzzz_yz[k];

                g_yy_0_yzzzz_z[k] = -g_yy_0_yzzz_z[k] * ab_z + g_yy_0_yzzz_zz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 249 * ccomps * dcomps);

            auto g_yy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 250 * ccomps * dcomps);

            auto g_yy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzz_x, g_yy_0_zzzz_xz, g_yy_0_zzzz_y, g_yy_0_zzzz_yz, g_yy_0_zzzz_z, g_yy_0_zzzz_zz, g_yy_0_zzzzz_x, g_yy_0_zzzzz_y, g_yy_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzz_x[k] = -g_yy_0_zzzz_x[k] * ab_z + g_yy_0_zzzz_xz[k];

                g_yy_0_zzzzz_y[k] = -g_yy_0_zzzz_y[k] * ab_z + g_yy_0_zzzz_yz[k];

                g_yy_0_zzzzz_z[k] = -g_yy_0_zzzz_z[k] * ab_z + g_yy_0_zzzz_zz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxx_x, g_yz_0_xxxx_xx, g_yz_0_xxxx_xy, g_yz_0_xxxx_xz, g_yz_0_xxxx_y, g_yz_0_xxxx_z, g_yz_0_xxxxx_x, g_yz_0_xxxxx_y, g_yz_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxx_x[k] = -g_yz_0_xxxx_x[k] * ab_x + g_yz_0_xxxx_xx[k];

                g_yz_0_xxxxx_y[k] = -g_yz_0_xxxx_y[k] * ab_x + g_yz_0_xxxx_xy[k];

                g_yz_0_xxxxx_z[k] = -g_yz_0_xxxx_z[k] * ab_x + g_yz_0_xxxx_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxy_x, g_yz_0_xxxxy_y, g_yz_0_xxxxy_z, g_yz_0_xxxy_x, g_yz_0_xxxy_xx, g_yz_0_xxxy_xy, g_yz_0_xxxy_xz, g_yz_0_xxxy_y, g_yz_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxy_x[k] = -g_yz_0_xxxy_x[k] * ab_x + g_yz_0_xxxy_xx[k];

                g_yz_0_xxxxy_y[k] = -g_yz_0_xxxy_y[k] * ab_x + g_yz_0_xxxy_xy[k];

                g_yz_0_xxxxy_z[k] = -g_yz_0_xxxy_z[k] * ab_x + g_yz_0_xxxy_xz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxz_x, g_yz_0_xxxxz_y, g_yz_0_xxxxz_z, g_yz_0_xxxz_x, g_yz_0_xxxz_xx, g_yz_0_xxxz_xy, g_yz_0_xxxz_xz, g_yz_0_xxxz_y, g_yz_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxz_x[k] = -g_yz_0_xxxz_x[k] * ab_x + g_yz_0_xxxz_xx[k];

                g_yz_0_xxxxz_y[k] = -g_yz_0_xxxz_y[k] * ab_x + g_yz_0_xxxz_xy[k];

                g_yz_0_xxxxz_z[k] = -g_yz_0_xxxz_z[k] * ab_x + g_yz_0_xxxz_xz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyy_x, g_yz_0_xxxyy_y, g_yz_0_xxxyy_z, g_yz_0_xxyy_x, g_yz_0_xxyy_xx, g_yz_0_xxyy_xy, g_yz_0_xxyy_xz, g_yz_0_xxyy_y, g_yz_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyy_x[k] = -g_yz_0_xxyy_x[k] * ab_x + g_yz_0_xxyy_xx[k];

                g_yz_0_xxxyy_y[k] = -g_yz_0_xxyy_y[k] * ab_x + g_yz_0_xxyy_xy[k];

                g_yz_0_xxxyy_z[k] = -g_yz_0_xxyy_z[k] * ab_x + g_yz_0_xxyy_xz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyz_x, g_yz_0_xxxyz_y, g_yz_0_xxxyz_z, g_yz_0_xxyz_x, g_yz_0_xxyz_xx, g_yz_0_xxyz_xy, g_yz_0_xxyz_xz, g_yz_0_xxyz_y, g_yz_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyz_x[k] = -g_yz_0_xxyz_x[k] * ab_x + g_yz_0_xxyz_xx[k];

                g_yz_0_xxxyz_y[k] = -g_yz_0_xxyz_y[k] * ab_x + g_yz_0_xxyz_xy[k];

                g_yz_0_xxxyz_z[k] = -g_yz_0_xxyz_z[k] * ab_x + g_yz_0_xxyz_xz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzz_x, g_yz_0_xxxzz_y, g_yz_0_xxxzz_z, g_yz_0_xxzz_x, g_yz_0_xxzz_xx, g_yz_0_xxzz_xy, g_yz_0_xxzz_xz, g_yz_0_xxzz_y, g_yz_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzz_x[k] = -g_yz_0_xxzz_x[k] * ab_x + g_yz_0_xxzz_xx[k];

                g_yz_0_xxxzz_y[k] = -g_yz_0_xxzz_y[k] * ab_x + g_yz_0_xxzz_xy[k];

                g_yz_0_xxxzz_z[k] = -g_yz_0_xxzz_z[k] * ab_x + g_yz_0_xxzz_xz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyy_x, g_yz_0_xxyyy_y, g_yz_0_xxyyy_z, g_yz_0_xyyy_x, g_yz_0_xyyy_xx, g_yz_0_xyyy_xy, g_yz_0_xyyy_xz, g_yz_0_xyyy_y, g_yz_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyy_x[k] = -g_yz_0_xyyy_x[k] * ab_x + g_yz_0_xyyy_xx[k];

                g_yz_0_xxyyy_y[k] = -g_yz_0_xyyy_y[k] * ab_x + g_yz_0_xyyy_xy[k];

                g_yz_0_xxyyy_z[k] = -g_yz_0_xyyy_z[k] * ab_x + g_yz_0_xyyy_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyz_x, g_yz_0_xxyyz_y, g_yz_0_xxyyz_z, g_yz_0_xyyz_x, g_yz_0_xyyz_xx, g_yz_0_xyyz_xy, g_yz_0_xyyz_xz, g_yz_0_xyyz_y, g_yz_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyz_x[k] = -g_yz_0_xyyz_x[k] * ab_x + g_yz_0_xyyz_xx[k];

                g_yz_0_xxyyz_y[k] = -g_yz_0_xyyz_y[k] * ab_x + g_yz_0_xyyz_xy[k];

                g_yz_0_xxyyz_z[k] = -g_yz_0_xyyz_z[k] * ab_x + g_yz_0_xyyz_xz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzz_x, g_yz_0_xxyzz_y, g_yz_0_xxyzz_z, g_yz_0_xyzz_x, g_yz_0_xyzz_xx, g_yz_0_xyzz_xy, g_yz_0_xyzz_xz, g_yz_0_xyzz_y, g_yz_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzz_x[k] = -g_yz_0_xyzz_x[k] * ab_x + g_yz_0_xyzz_xx[k];

                g_yz_0_xxyzz_y[k] = -g_yz_0_xyzz_y[k] * ab_x + g_yz_0_xyzz_xy[k];

                g_yz_0_xxyzz_z[k] = -g_yz_0_xyzz_z[k] * ab_x + g_yz_0_xyzz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzz_x, g_yz_0_xxzzz_y, g_yz_0_xxzzz_z, g_yz_0_xzzz_x, g_yz_0_xzzz_xx, g_yz_0_xzzz_xy, g_yz_0_xzzz_xz, g_yz_0_xzzz_y, g_yz_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzz_x[k] = -g_yz_0_xzzz_x[k] * ab_x + g_yz_0_xzzz_xx[k];

                g_yz_0_xxzzz_y[k] = -g_yz_0_xzzz_y[k] * ab_x + g_yz_0_xzzz_xy[k];

                g_yz_0_xxzzz_z[k] = -g_yz_0_xzzz_z[k] * ab_x + g_yz_0_xzzz_xz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyy_x, g_yz_0_xyyyy_y, g_yz_0_xyyyy_z, g_yz_0_yyyy_x, g_yz_0_yyyy_xx, g_yz_0_yyyy_xy, g_yz_0_yyyy_xz, g_yz_0_yyyy_y, g_yz_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyy_x[k] = -g_yz_0_yyyy_x[k] * ab_x + g_yz_0_yyyy_xx[k];

                g_yz_0_xyyyy_y[k] = -g_yz_0_yyyy_y[k] * ab_x + g_yz_0_yyyy_xy[k];

                g_yz_0_xyyyy_z[k] = -g_yz_0_yyyy_z[k] * ab_x + g_yz_0_yyyy_xz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyz_x, g_yz_0_xyyyz_y, g_yz_0_xyyyz_z, g_yz_0_yyyz_x, g_yz_0_yyyz_xx, g_yz_0_yyyz_xy, g_yz_0_yyyz_xz, g_yz_0_yyyz_y, g_yz_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyz_x[k] = -g_yz_0_yyyz_x[k] * ab_x + g_yz_0_yyyz_xx[k];

                g_yz_0_xyyyz_y[k] = -g_yz_0_yyyz_y[k] * ab_x + g_yz_0_yyyz_xy[k];

                g_yz_0_xyyyz_z[k] = -g_yz_0_yyyz_z[k] * ab_x + g_yz_0_yyyz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzz_x, g_yz_0_xyyzz_y, g_yz_0_xyyzz_z, g_yz_0_yyzz_x, g_yz_0_yyzz_xx, g_yz_0_yyzz_xy, g_yz_0_yyzz_xz, g_yz_0_yyzz_y, g_yz_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzz_x[k] = -g_yz_0_yyzz_x[k] * ab_x + g_yz_0_yyzz_xx[k];

                g_yz_0_xyyzz_y[k] = -g_yz_0_yyzz_y[k] * ab_x + g_yz_0_yyzz_xy[k];

                g_yz_0_xyyzz_z[k] = -g_yz_0_yyzz_z[k] * ab_x + g_yz_0_yyzz_xz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzz_x, g_yz_0_xyzzz_y, g_yz_0_xyzzz_z, g_yz_0_yzzz_x, g_yz_0_yzzz_xx, g_yz_0_yzzz_xy, g_yz_0_yzzz_xz, g_yz_0_yzzz_y, g_yz_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzz_x[k] = -g_yz_0_yzzz_x[k] * ab_x + g_yz_0_yzzz_xx[k];

                g_yz_0_xyzzz_y[k] = -g_yz_0_yzzz_y[k] * ab_x + g_yz_0_yzzz_xy[k];

                g_yz_0_xyzzz_z[k] = -g_yz_0_yzzz_z[k] * ab_x + g_yz_0_yzzz_xz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzz_x, g_yz_0_xzzzz_y, g_yz_0_xzzzz_z, g_yz_0_zzzz_x, g_yz_0_zzzz_xx, g_yz_0_zzzz_xy, g_yz_0_zzzz_xz, g_yz_0_zzzz_y, g_yz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzz_x[k] = -g_yz_0_zzzz_x[k] * ab_x + g_yz_0_zzzz_xx[k];

                g_yz_0_xzzzz_y[k] = -g_yz_0_zzzz_y[k] * ab_x + g_yz_0_zzzz_xy[k];

                g_yz_0_xzzzz_z[k] = -g_yz_0_zzzz_z[k] * ab_x + g_yz_0_zzzz_xz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyy_x, g_yz_0_yyyy_xy, g_yz_0_yyyy_y, g_yz_0_yyyy_yy, g_yz_0_yyyy_yz, g_yz_0_yyyy_z, g_yz_0_yyyyy_x, g_yz_0_yyyyy_y, g_yz_0_yyyyy_z, g_z_0_yyyy_x, g_z_0_yyyy_y, g_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyy_x[k] = -g_z_0_yyyy_x[k] - g_yz_0_yyyy_x[k] * ab_y + g_yz_0_yyyy_xy[k];

                g_yz_0_yyyyy_y[k] = -g_z_0_yyyy_y[k] - g_yz_0_yyyy_y[k] * ab_y + g_yz_0_yyyy_yy[k];

                g_yz_0_yyyyy_z[k] = -g_z_0_yyyy_z[k] - g_yz_0_yyyy_z[k] * ab_y + g_yz_0_yyyy_yz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 300 * ccomps * dcomps);

            auto g_yz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 301 * ccomps * dcomps);

            auto g_yz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyz_x, g_yz_0_yyyyz_y, g_yz_0_yyyyz_z, g_yz_0_yyyz_x, g_yz_0_yyyz_xy, g_yz_0_yyyz_y, g_yz_0_yyyz_yy, g_yz_0_yyyz_yz, g_yz_0_yyyz_z, g_z_0_yyyz_x, g_z_0_yyyz_y, g_z_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyz_x[k] = -g_z_0_yyyz_x[k] - g_yz_0_yyyz_x[k] * ab_y + g_yz_0_yyyz_xy[k];

                g_yz_0_yyyyz_y[k] = -g_z_0_yyyz_y[k] - g_yz_0_yyyz_y[k] * ab_y + g_yz_0_yyyz_yy[k];

                g_yz_0_yyyyz_z[k] = -g_z_0_yyyz_z[k] - g_yz_0_yyyz_z[k] * ab_y + g_yz_0_yyyz_yz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 303 * ccomps * dcomps);

            auto g_yz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 304 * ccomps * dcomps);

            auto g_yz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzz_x, g_yz_0_yyyzz_y, g_yz_0_yyyzz_z, g_yz_0_yyzz_x, g_yz_0_yyzz_xy, g_yz_0_yyzz_y, g_yz_0_yyzz_yy, g_yz_0_yyzz_yz, g_yz_0_yyzz_z, g_z_0_yyzz_x, g_z_0_yyzz_y, g_z_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzz_x[k] = -g_z_0_yyzz_x[k] - g_yz_0_yyzz_x[k] * ab_y + g_yz_0_yyzz_xy[k];

                g_yz_0_yyyzz_y[k] = -g_z_0_yyzz_y[k] - g_yz_0_yyzz_y[k] * ab_y + g_yz_0_yyzz_yy[k];

                g_yz_0_yyyzz_z[k] = -g_z_0_yyzz_z[k] - g_yz_0_yyzz_z[k] * ab_y + g_yz_0_yyzz_yz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 306 * ccomps * dcomps);

            auto g_yz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 307 * ccomps * dcomps);

            auto g_yz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzz_x, g_yz_0_yyzzz_y, g_yz_0_yyzzz_z, g_yz_0_yzzz_x, g_yz_0_yzzz_xy, g_yz_0_yzzz_y, g_yz_0_yzzz_yy, g_yz_0_yzzz_yz, g_yz_0_yzzz_z, g_z_0_yzzz_x, g_z_0_yzzz_y, g_z_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzz_x[k] = -g_z_0_yzzz_x[k] - g_yz_0_yzzz_x[k] * ab_y + g_yz_0_yzzz_xy[k];

                g_yz_0_yyzzz_y[k] = -g_z_0_yzzz_y[k] - g_yz_0_yzzz_y[k] * ab_y + g_yz_0_yzzz_yy[k];

                g_yz_0_yyzzz_z[k] = -g_z_0_yzzz_z[k] - g_yz_0_yzzz_z[k] * ab_y + g_yz_0_yzzz_yz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 309 * ccomps * dcomps);

            auto g_yz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 310 * ccomps * dcomps);

            auto g_yz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzz_x, g_yz_0_yzzzz_y, g_yz_0_yzzzz_z, g_yz_0_zzzz_x, g_yz_0_zzzz_xy, g_yz_0_zzzz_y, g_yz_0_zzzz_yy, g_yz_0_zzzz_yz, g_yz_0_zzzz_z, g_z_0_zzzz_x, g_z_0_zzzz_y, g_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzz_x[k] = -g_z_0_zzzz_x[k] - g_yz_0_zzzz_x[k] * ab_y + g_yz_0_zzzz_xy[k];

                g_yz_0_yzzzz_y[k] = -g_z_0_zzzz_y[k] - g_yz_0_zzzz_y[k] * ab_y + g_yz_0_zzzz_yy[k];

                g_yz_0_yzzzz_z[k] = -g_z_0_zzzz_z[k] - g_yz_0_zzzz_z[k] * ab_y + g_yz_0_zzzz_yz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 312 * ccomps * dcomps);

            auto g_yz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 313 * ccomps * dcomps);

            auto g_yz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_x, g_y_0_zzzz_y, g_y_0_zzzz_z, g_yz_0_zzzz_x, g_yz_0_zzzz_xz, g_yz_0_zzzz_y, g_yz_0_zzzz_yz, g_yz_0_zzzz_z, g_yz_0_zzzz_zz, g_yz_0_zzzzz_x, g_yz_0_zzzzz_y, g_yz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzz_x[k] = -g_y_0_zzzz_x[k] - g_yz_0_zzzz_x[k] * ab_z + g_yz_0_zzzz_xz[k];

                g_yz_0_zzzzz_y[k] = -g_y_0_zzzz_y[k] - g_yz_0_zzzz_y[k] * ab_z + g_yz_0_zzzz_yz[k];

                g_yz_0_zzzzz_z[k] = -g_y_0_zzzz_z[k] - g_yz_0_zzzz_z[k] * ab_z + g_yz_0_zzzz_zz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxx_x, g_zz_0_xxxx_xx, g_zz_0_xxxx_xy, g_zz_0_xxxx_xz, g_zz_0_xxxx_y, g_zz_0_xxxx_z, g_zz_0_xxxxx_x, g_zz_0_xxxxx_y, g_zz_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxx_x[k] = -g_zz_0_xxxx_x[k] * ab_x + g_zz_0_xxxx_xx[k];

                g_zz_0_xxxxx_y[k] = -g_zz_0_xxxx_y[k] * ab_x + g_zz_0_xxxx_xy[k];

                g_zz_0_xxxxx_z[k] = -g_zz_0_xxxx_z[k] * ab_x + g_zz_0_xxxx_xz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxy_x, g_zz_0_xxxxy_y, g_zz_0_xxxxy_z, g_zz_0_xxxy_x, g_zz_0_xxxy_xx, g_zz_0_xxxy_xy, g_zz_0_xxxy_xz, g_zz_0_xxxy_y, g_zz_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxy_x[k] = -g_zz_0_xxxy_x[k] * ab_x + g_zz_0_xxxy_xx[k];

                g_zz_0_xxxxy_y[k] = -g_zz_0_xxxy_y[k] * ab_x + g_zz_0_xxxy_xy[k];

                g_zz_0_xxxxy_z[k] = -g_zz_0_xxxy_z[k] * ab_x + g_zz_0_xxxy_xz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxz_x, g_zz_0_xxxxz_y, g_zz_0_xxxxz_z, g_zz_0_xxxz_x, g_zz_0_xxxz_xx, g_zz_0_xxxz_xy, g_zz_0_xxxz_xz, g_zz_0_xxxz_y, g_zz_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxz_x[k] = -g_zz_0_xxxz_x[k] * ab_x + g_zz_0_xxxz_xx[k];

                g_zz_0_xxxxz_y[k] = -g_zz_0_xxxz_y[k] * ab_x + g_zz_0_xxxz_xy[k];

                g_zz_0_xxxxz_z[k] = -g_zz_0_xxxz_z[k] * ab_x + g_zz_0_xxxz_xz[k];
            }

            /// Set up 324-327 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 326 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyy_x, g_zz_0_xxxyy_y, g_zz_0_xxxyy_z, g_zz_0_xxyy_x, g_zz_0_xxyy_xx, g_zz_0_xxyy_xy, g_zz_0_xxyy_xz, g_zz_0_xxyy_y, g_zz_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyy_x[k] = -g_zz_0_xxyy_x[k] * ab_x + g_zz_0_xxyy_xx[k];

                g_zz_0_xxxyy_y[k] = -g_zz_0_xxyy_y[k] * ab_x + g_zz_0_xxyy_xy[k];

                g_zz_0_xxxyy_z[k] = -g_zz_0_xxyy_z[k] * ab_x + g_zz_0_xxyy_xz[k];
            }

            /// Set up 327-330 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyz_x, g_zz_0_xxxyz_y, g_zz_0_xxxyz_z, g_zz_0_xxyz_x, g_zz_0_xxyz_xx, g_zz_0_xxyz_xy, g_zz_0_xxyz_xz, g_zz_0_xxyz_y, g_zz_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyz_x[k] = -g_zz_0_xxyz_x[k] * ab_x + g_zz_0_xxyz_xx[k];

                g_zz_0_xxxyz_y[k] = -g_zz_0_xxyz_y[k] * ab_x + g_zz_0_xxyz_xy[k];

                g_zz_0_xxxyz_z[k] = -g_zz_0_xxyz_z[k] * ab_x + g_zz_0_xxyz_xz[k];
            }

            /// Set up 330-333 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 332 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzz_x, g_zz_0_xxxzz_y, g_zz_0_xxxzz_z, g_zz_0_xxzz_x, g_zz_0_xxzz_xx, g_zz_0_xxzz_xy, g_zz_0_xxzz_xz, g_zz_0_xxzz_y, g_zz_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzz_x[k] = -g_zz_0_xxzz_x[k] * ab_x + g_zz_0_xxzz_xx[k];

                g_zz_0_xxxzz_y[k] = -g_zz_0_xxzz_y[k] * ab_x + g_zz_0_xxzz_xy[k];

                g_zz_0_xxxzz_z[k] = -g_zz_0_xxzz_z[k] * ab_x + g_zz_0_xxzz_xz[k];
            }

            /// Set up 333-336 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyy_x, g_zz_0_xxyyy_y, g_zz_0_xxyyy_z, g_zz_0_xyyy_x, g_zz_0_xyyy_xx, g_zz_0_xyyy_xy, g_zz_0_xyyy_xz, g_zz_0_xyyy_y, g_zz_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyy_x[k] = -g_zz_0_xyyy_x[k] * ab_x + g_zz_0_xyyy_xx[k];

                g_zz_0_xxyyy_y[k] = -g_zz_0_xyyy_y[k] * ab_x + g_zz_0_xyyy_xy[k];

                g_zz_0_xxyyy_z[k] = -g_zz_0_xyyy_z[k] * ab_x + g_zz_0_xyyy_xz[k];
            }

            /// Set up 336-339 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 338 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyz_x, g_zz_0_xxyyz_y, g_zz_0_xxyyz_z, g_zz_0_xyyz_x, g_zz_0_xyyz_xx, g_zz_0_xyyz_xy, g_zz_0_xyyz_xz, g_zz_0_xyyz_y, g_zz_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyz_x[k] = -g_zz_0_xyyz_x[k] * ab_x + g_zz_0_xyyz_xx[k];

                g_zz_0_xxyyz_y[k] = -g_zz_0_xyyz_y[k] * ab_x + g_zz_0_xyyz_xy[k];

                g_zz_0_xxyyz_z[k] = -g_zz_0_xyyz_z[k] * ab_x + g_zz_0_xyyz_xz[k];
            }

            /// Set up 339-342 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzz_x, g_zz_0_xxyzz_y, g_zz_0_xxyzz_z, g_zz_0_xyzz_x, g_zz_0_xyzz_xx, g_zz_0_xyzz_xy, g_zz_0_xyzz_xz, g_zz_0_xyzz_y, g_zz_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzz_x[k] = -g_zz_0_xyzz_x[k] * ab_x + g_zz_0_xyzz_xx[k];

                g_zz_0_xxyzz_y[k] = -g_zz_0_xyzz_y[k] * ab_x + g_zz_0_xyzz_xy[k];

                g_zz_0_xxyzz_z[k] = -g_zz_0_xyzz_z[k] * ab_x + g_zz_0_xyzz_xz[k];
            }

            /// Set up 342-345 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzz_x, g_zz_0_xxzzz_y, g_zz_0_xxzzz_z, g_zz_0_xzzz_x, g_zz_0_xzzz_xx, g_zz_0_xzzz_xy, g_zz_0_xzzz_xz, g_zz_0_xzzz_y, g_zz_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzz_x[k] = -g_zz_0_xzzz_x[k] * ab_x + g_zz_0_xzzz_xx[k];

                g_zz_0_xxzzz_y[k] = -g_zz_0_xzzz_y[k] * ab_x + g_zz_0_xzzz_xy[k];

                g_zz_0_xxzzz_z[k] = -g_zz_0_xzzz_z[k] * ab_x + g_zz_0_xzzz_xz[k];
            }

            /// Set up 345-348 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyy_x, g_zz_0_xyyyy_y, g_zz_0_xyyyy_z, g_zz_0_yyyy_x, g_zz_0_yyyy_xx, g_zz_0_yyyy_xy, g_zz_0_yyyy_xz, g_zz_0_yyyy_y, g_zz_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyy_x[k] = -g_zz_0_yyyy_x[k] * ab_x + g_zz_0_yyyy_xx[k];

                g_zz_0_xyyyy_y[k] = -g_zz_0_yyyy_y[k] * ab_x + g_zz_0_yyyy_xy[k];

                g_zz_0_xyyyy_z[k] = -g_zz_0_yyyy_z[k] * ab_x + g_zz_0_yyyy_xz[k];
            }

            /// Set up 348-351 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 350 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyz_x, g_zz_0_xyyyz_y, g_zz_0_xyyyz_z, g_zz_0_yyyz_x, g_zz_0_yyyz_xx, g_zz_0_yyyz_xy, g_zz_0_yyyz_xz, g_zz_0_yyyz_y, g_zz_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyz_x[k] = -g_zz_0_yyyz_x[k] * ab_x + g_zz_0_yyyz_xx[k];

                g_zz_0_xyyyz_y[k] = -g_zz_0_yyyz_y[k] * ab_x + g_zz_0_yyyz_xy[k];

                g_zz_0_xyyyz_z[k] = -g_zz_0_yyyz_z[k] * ab_x + g_zz_0_yyyz_xz[k];
            }

            /// Set up 351-354 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzz_x, g_zz_0_xyyzz_y, g_zz_0_xyyzz_z, g_zz_0_yyzz_x, g_zz_0_yyzz_xx, g_zz_0_yyzz_xy, g_zz_0_yyzz_xz, g_zz_0_yyzz_y, g_zz_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzz_x[k] = -g_zz_0_yyzz_x[k] * ab_x + g_zz_0_yyzz_xx[k];

                g_zz_0_xyyzz_y[k] = -g_zz_0_yyzz_y[k] * ab_x + g_zz_0_yyzz_xy[k];

                g_zz_0_xyyzz_z[k] = -g_zz_0_yyzz_z[k] * ab_x + g_zz_0_yyzz_xz[k];
            }

            /// Set up 354-357 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzz_x, g_zz_0_xyzzz_y, g_zz_0_xyzzz_z, g_zz_0_yzzz_x, g_zz_0_yzzz_xx, g_zz_0_yzzz_xy, g_zz_0_yzzz_xz, g_zz_0_yzzz_y, g_zz_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzz_x[k] = -g_zz_0_yzzz_x[k] * ab_x + g_zz_0_yzzz_xx[k];

                g_zz_0_xyzzz_y[k] = -g_zz_0_yzzz_y[k] * ab_x + g_zz_0_yzzz_xy[k];

                g_zz_0_xyzzz_z[k] = -g_zz_0_yzzz_z[k] * ab_x + g_zz_0_yzzz_xz[k];
            }

            /// Set up 357-360 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzz_x, g_zz_0_xzzzz_y, g_zz_0_xzzzz_z, g_zz_0_zzzz_x, g_zz_0_zzzz_xx, g_zz_0_zzzz_xy, g_zz_0_zzzz_xz, g_zz_0_zzzz_y, g_zz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzz_x[k] = -g_zz_0_zzzz_x[k] * ab_x + g_zz_0_zzzz_xx[k];

                g_zz_0_xzzzz_y[k] = -g_zz_0_zzzz_y[k] * ab_x + g_zz_0_zzzz_xy[k];

                g_zz_0_xzzzz_z[k] = -g_zz_0_zzzz_z[k] * ab_x + g_zz_0_zzzz_xz[k];
            }

            /// Set up 360-363 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 360 * ccomps * dcomps);

            auto g_zz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 361 * ccomps * dcomps);

            auto g_zz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 362 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyy_x, g_zz_0_yyyy_xy, g_zz_0_yyyy_y, g_zz_0_yyyy_yy, g_zz_0_yyyy_yz, g_zz_0_yyyy_z, g_zz_0_yyyyy_x, g_zz_0_yyyyy_y, g_zz_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyy_x[k] = -g_zz_0_yyyy_x[k] * ab_y + g_zz_0_yyyy_xy[k];

                g_zz_0_yyyyy_y[k] = -g_zz_0_yyyy_y[k] * ab_y + g_zz_0_yyyy_yy[k];

                g_zz_0_yyyyy_z[k] = -g_zz_0_yyyy_z[k] * ab_y + g_zz_0_yyyy_yz[k];
            }

            /// Set up 363-366 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 363 * ccomps * dcomps);

            auto g_zz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 364 * ccomps * dcomps);

            auto g_zz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyz_x, g_zz_0_yyyyz_y, g_zz_0_yyyyz_z, g_zz_0_yyyz_x, g_zz_0_yyyz_xy, g_zz_0_yyyz_y, g_zz_0_yyyz_yy, g_zz_0_yyyz_yz, g_zz_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyz_x[k] = -g_zz_0_yyyz_x[k] * ab_y + g_zz_0_yyyz_xy[k];

                g_zz_0_yyyyz_y[k] = -g_zz_0_yyyz_y[k] * ab_y + g_zz_0_yyyz_yy[k];

                g_zz_0_yyyyz_z[k] = -g_zz_0_yyyz_z[k] * ab_y + g_zz_0_yyyz_yz[k];
            }

            /// Set up 366-369 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 366 * ccomps * dcomps);

            auto g_zz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 367 * ccomps * dcomps);

            auto g_zz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 368 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzz_x, g_zz_0_yyyzz_y, g_zz_0_yyyzz_z, g_zz_0_yyzz_x, g_zz_0_yyzz_xy, g_zz_0_yyzz_y, g_zz_0_yyzz_yy, g_zz_0_yyzz_yz, g_zz_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzz_x[k] = -g_zz_0_yyzz_x[k] * ab_y + g_zz_0_yyzz_xy[k];

                g_zz_0_yyyzz_y[k] = -g_zz_0_yyzz_y[k] * ab_y + g_zz_0_yyzz_yy[k];

                g_zz_0_yyyzz_z[k] = -g_zz_0_yyzz_z[k] * ab_y + g_zz_0_yyzz_yz[k];
            }

            /// Set up 369-372 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 369 * ccomps * dcomps);

            auto g_zz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 370 * ccomps * dcomps);

            auto g_zz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzz_x, g_zz_0_yyzzz_y, g_zz_0_yyzzz_z, g_zz_0_yzzz_x, g_zz_0_yzzz_xy, g_zz_0_yzzz_y, g_zz_0_yzzz_yy, g_zz_0_yzzz_yz, g_zz_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzz_x[k] = -g_zz_0_yzzz_x[k] * ab_y + g_zz_0_yzzz_xy[k];

                g_zz_0_yyzzz_y[k] = -g_zz_0_yzzz_y[k] * ab_y + g_zz_0_yzzz_yy[k];

                g_zz_0_yyzzz_z[k] = -g_zz_0_yzzz_z[k] * ab_y + g_zz_0_yzzz_yz[k];
            }

            /// Set up 372-375 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 372 * ccomps * dcomps);

            auto g_zz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 373 * ccomps * dcomps);

            auto g_zz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzz_x, g_zz_0_yzzzz_y, g_zz_0_yzzzz_z, g_zz_0_zzzz_x, g_zz_0_zzzz_xy, g_zz_0_zzzz_y, g_zz_0_zzzz_yy, g_zz_0_zzzz_yz, g_zz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzz_x[k] = -g_zz_0_zzzz_x[k] * ab_y + g_zz_0_zzzz_xy[k];

                g_zz_0_yzzzz_y[k] = -g_zz_0_zzzz_y[k] * ab_y + g_zz_0_zzzz_yy[k];

                g_zz_0_yzzzz_z[k] = -g_zz_0_zzzz_z[k] * ab_y + g_zz_0_zzzz_yz[k];
            }

            /// Set up 375-378 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 375 * ccomps * dcomps);

            auto g_zz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 376 * ccomps * dcomps);

            auto g_zz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_x, g_z_0_zzzz_y, g_z_0_zzzz_z, g_zz_0_zzzz_x, g_zz_0_zzzz_xz, g_zz_0_zzzz_y, g_zz_0_zzzz_yz, g_zz_0_zzzz_z, g_zz_0_zzzz_zz, g_zz_0_zzzzz_x, g_zz_0_zzzzz_y, g_zz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzz_x[k] = -2.0 * g_z_0_zzzz_x[k] - g_zz_0_zzzz_x[k] * ab_z + g_zz_0_zzzz_xz[k];

                g_zz_0_zzzzz_y[k] = -2.0 * g_z_0_zzzz_y[k] - g_zz_0_zzzz_y[k] * ab_z + g_zz_0_zzzz_yz[k];

                g_zz_0_zzzzz_z[k] = -2.0 * g_z_0_zzzz_z[k] - g_zz_0_zzzz_z[k] * ab_z + g_zz_0_zzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

