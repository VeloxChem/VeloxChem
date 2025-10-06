#include "ElectronRepulsionGeom0100ContrRecHPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_hpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hpxx,
                                            const size_t idx_gpxx,
                                            const size_t idx_geom_01_gpxx,
                                            const size_t idx_geom_01_gdxx,
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

            const auto gp_off = idx_gpxx + i * dcomps + j;

            auto g_xxxx_x = cbuffer.data(gp_off + 0 * ccomps * dcomps);

            auto g_xxxx_y = cbuffer.data(gp_off + 1 * ccomps * dcomps);

            auto g_xxxx_z = cbuffer.data(gp_off + 2 * ccomps * dcomps);

            auto g_xxxy_x = cbuffer.data(gp_off + 3 * ccomps * dcomps);

            auto g_xxxy_y = cbuffer.data(gp_off + 4 * ccomps * dcomps);

            auto g_xxxy_z = cbuffer.data(gp_off + 5 * ccomps * dcomps);

            auto g_xxxz_x = cbuffer.data(gp_off + 6 * ccomps * dcomps);

            auto g_xxxz_y = cbuffer.data(gp_off + 7 * ccomps * dcomps);

            auto g_xxxz_z = cbuffer.data(gp_off + 8 * ccomps * dcomps);

            auto g_xxyy_x = cbuffer.data(gp_off + 9 * ccomps * dcomps);

            auto g_xxyy_y = cbuffer.data(gp_off + 10 * ccomps * dcomps);

            auto g_xxyy_z = cbuffer.data(gp_off + 11 * ccomps * dcomps);

            auto g_xxyz_x = cbuffer.data(gp_off + 12 * ccomps * dcomps);

            auto g_xxyz_y = cbuffer.data(gp_off + 13 * ccomps * dcomps);

            auto g_xxyz_z = cbuffer.data(gp_off + 14 * ccomps * dcomps);

            auto g_xxzz_x = cbuffer.data(gp_off + 15 * ccomps * dcomps);

            auto g_xxzz_y = cbuffer.data(gp_off + 16 * ccomps * dcomps);

            auto g_xxzz_z = cbuffer.data(gp_off + 17 * ccomps * dcomps);

            auto g_xyyy_x = cbuffer.data(gp_off + 18 * ccomps * dcomps);

            auto g_xyyy_y = cbuffer.data(gp_off + 19 * ccomps * dcomps);

            auto g_xyyy_z = cbuffer.data(gp_off + 20 * ccomps * dcomps);

            auto g_xyyz_x = cbuffer.data(gp_off + 21 * ccomps * dcomps);

            auto g_xyyz_y = cbuffer.data(gp_off + 22 * ccomps * dcomps);

            auto g_xyyz_z = cbuffer.data(gp_off + 23 * ccomps * dcomps);

            auto g_xyzz_x = cbuffer.data(gp_off + 24 * ccomps * dcomps);

            auto g_xyzz_y = cbuffer.data(gp_off + 25 * ccomps * dcomps);

            auto g_xyzz_z = cbuffer.data(gp_off + 26 * ccomps * dcomps);

            auto g_xzzz_x = cbuffer.data(gp_off + 27 * ccomps * dcomps);

            auto g_xzzz_y = cbuffer.data(gp_off + 28 * ccomps * dcomps);

            auto g_xzzz_z = cbuffer.data(gp_off + 29 * ccomps * dcomps);

            auto g_yyyy_x = cbuffer.data(gp_off + 30 * ccomps * dcomps);

            auto g_yyyy_y = cbuffer.data(gp_off + 31 * ccomps * dcomps);

            auto g_yyyy_z = cbuffer.data(gp_off + 32 * ccomps * dcomps);

            auto g_yyyz_x = cbuffer.data(gp_off + 33 * ccomps * dcomps);

            auto g_yyyz_y = cbuffer.data(gp_off + 34 * ccomps * dcomps);

            auto g_yyyz_z = cbuffer.data(gp_off + 35 * ccomps * dcomps);

            auto g_yyzz_x = cbuffer.data(gp_off + 36 * ccomps * dcomps);

            auto g_yyzz_y = cbuffer.data(gp_off + 37 * ccomps * dcomps);

            auto g_yyzz_z = cbuffer.data(gp_off + 38 * ccomps * dcomps);

            auto g_yzzz_x = cbuffer.data(gp_off + 39 * ccomps * dcomps);

            auto g_yzzz_y = cbuffer.data(gp_off + 40 * ccomps * dcomps);

            auto g_yzzz_z = cbuffer.data(gp_off + 41 * ccomps * dcomps);

            auto g_zzzz_x = cbuffer.data(gp_off + 42 * ccomps * dcomps);

            auto g_zzzz_y = cbuffer.data(gp_off + 43 * ccomps * dcomps);

            auto g_zzzz_z = cbuffer.data(gp_off + 44 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GPSS

            const auto gp_geom_01_off = idx_geom_01_gpxx + i * dcomps + j;

            auto g_0_x_xxxx_x = cbuffer.data(gp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_y = cbuffer.data(gp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_z = cbuffer.data(gp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxy_x = cbuffer.data(gp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxy_y = cbuffer.data(gp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxy_z = cbuffer.data(gp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxz_x = cbuffer.data(gp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxz_y = cbuffer.data(gp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxz_z = cbuffer.data(gp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxyy_x = cbuffer.data(gp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxyy_y = cbuffer.data(gp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxyy_z = cbuffer.data(gp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxyz_x = cbuffer.data(gp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxyz_y = cbuffer.data(gp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxyz_z = cbuffer.data(gp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxzz_x = cbuffer.data(gp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxzz_y = cbuffer.data(gp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxzz_z = cbuffer.data(gp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xyyy_x = cbuffer.data(gp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xyyy_y = cbuffer.data(gp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xyyy_z = cbuffer.data(gp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xyyz_x = cbuffer.data(gp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xyyz_y = cbuffer.data(gp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xyyz_z = cbuffer.data(gp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xyzz_x = cbuffer.data(gp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xyzz_y = cbuffer.data(gp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xyzz_z = cbuffer.data(gp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xzzz_x = cbuffer.data(gp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xzzz_y = cbuffer.data(gp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xzzz_z = cbuffer.data(gp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_yyyy_x = cbuffer.data(gp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_yyyy_y = cbuffer.data(gp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_yyyy_z = cbuffer.data(gp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_yyyz_x = cbuffer.data(gp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_yyyz_y = cbuffer.data(gp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_yyyz_z = cbuffer.data(gp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_yyzz_x = cbuffer.data(gp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_yyzz_y = cbuffer.data(gp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_yyzz_z = cbuffer.data(gp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_yzzz_x = cbuffer.data(gp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_yzzz_y = cbuffer.data(gp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_yzzz_z = cbuffer.data(gp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_zzzz_x = cbuffer.data(gp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_zzzz_y = cbuffer.data(gp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_zzzz_z = cbuffer.data(gp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_xxxx_x = cbuffer.data(gp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_xxxx_y = cbuffer.data(gp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_xxxx_z = cbuffer.data(gp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_xxxy_x = cbuffer.data(gp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_xxxy_y = cbuffer.data(gp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_xxxy_z = cbuffer.data(gp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_xxxz_x = cbuffer.data(gp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_xxxz_y = cbuffer.data(gp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_xxxz_z = cbuffer.data(gp_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_xxyy_x = cbuffer.data(gp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_xxyy_y = cbuffer.data(gp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_xxyy_z = cbuffer.data(gp_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_xxyz_x = cbuffer.data(gp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_xxyz_y = cbuffer.data(gp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_xxyz_z = cbuffer.data(gp_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_xxzz_x = cbuffer.data(gp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_xxzz_y = cbuffer.data(gp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_xxzz_z = cbuffer.data(gp_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_xyyy_x = cbuffer.data(gp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xyyy_y = cbuffer.data(gp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xyyy_z = cbuffer.data(gp_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_xyyz_x = cbuffer.data(gp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xyyz_y = cbuffer.data(gp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xyyz_z = cbuffer.data(gp_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_xyzz_x = cbuffer.data(gp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xyzz_y = cbuffer.data(gp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xyzz_z = cbuffer.data(gp_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_xzzz_x = cbuffer.data(gp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xzzz_y = cbuffer.data(gp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xzzz_z = cbuffer.data(gp_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_yyyy_x = cbuffer.data(gp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_yyyy_y = cbuffer.data(gp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_yyyy_z = cbuffer.data(gp_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_yyyz_x = cbuffer.data(gp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_yyyz_y = cbuffer.data(gp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_yyyz_z = cbuffer.data(gp_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_yyzz_x = cbuffer.data(gp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_yyzz_y = cbuffer.data(gp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_yyzz_z = cbuffer.data(gp_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_yzzz_x = cbuffer.data(gp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_yzzz_y = cbuffer.data(gp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_yzzz_z = cbuffer.data(gp_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_zzzz_x = cbuffer.data(gp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_zzzz_y = cbuffer.data(gp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_zzzz_z = cbuffer.data(gp_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_z_xxxx_x = cbuffer.data(gp_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_z_xxxx_y = cbuffer.data(gp_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_z_xxxx_z = cbuffer.data(gp_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_z_xxxy_x = cbuffer.data(gp_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_z_xxxy_y = cbuffer.data(gp_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_z_xxxy_z = cbuffer.data(gp_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_z_xxxz_x = cbuffer.data(gp_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_z_xxxz_y = cbuffer.data(gp_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_z_xxxz_z = cbuffer.data(gp_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_z_xxyy_x = cbuffer.data(gp_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_z_xxyy_y = cbuffer.data(gp_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_z_xxyy_z = cbuffer.data(gp_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_z_xxyz_x = cbuffer.data(gp_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_z_xxyz_y = cbuffer.data(gp_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_z_xxyz_z = cbuffer.data(gp_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_z_xxzz_x = cbuffer.data(gp_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_z_xxzz_y = cbuffer.data(gp_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_z_xxzz_z = cbuffer.data(gp_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_z_xyyy_x = cbuffer.data(gp_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_z_xyyy_y = cbuffer.data(gp_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_z_xyyy_z = cbuffer.data(gp_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_z_xyyz_x = cbuffer.data(gp_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_z_xyyz_y = cbuffer.data(gp_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_z_xyyz_z = cbuffer.data(gp_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_z_xyzz_x = cbuffer.data(gp_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_z_xyzz_y = cbuffer.data(gp_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_z_xyzz_z = cbuffer.data(gp_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_z_xzzz_x = cbuffer.data(gp_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_z_xzzz_y = cbuffer.data(gp_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_z_xzzz_z = cbuffer.data(gp_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_z_yyyy_x = cbuffer.data(gp_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_yyyy_y = cbuffer.data(gp_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_yyyy_z = cbuffer.data(gp_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_z_yyyz_x = cbuffer.data(gp_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_yyyz_y = cbuffer.data(gp_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_yyyz_z = cbuffer.data(gp_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_yyzz_x = cbuffer.data(gp_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_yyzz_y = cbuffer.data(gp_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_yyzz_z = cbuffer.data(gp_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_yzzz_x = cbuffer.data(gp_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_yzzz_y = cbuffer.data(gp_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_yzzz_z = cbuffer.data(gp_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_zzzz_x = cbuffer.data(gp_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_zzzz_y = cbuffer.data(gp_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_zzzz_z = cbuffer.data(gp_geom_01_off + 134 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GDSS

            const auto gd_geom_01_off = idx_geom_01_gdxx + i * dcomps + j;

            auto g_0_x_xxxx_xx = cbuffer.data(gd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xy = cbuffer.data(gd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xz = cbuffer.data(gd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_yy = cbuffer.data(gd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_yz = cbuffer.data(gd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_zz = cbuffer.data(gd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxy_xx = cbuffer.data(gd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxy_xy = cbuffer.data(gd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxy_xz = cbuffer.data(gd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxy_yy = cbuffer.data(gd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxy_yz = cbuffer.data(gd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxy_zz = cbuffer.data(gd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxz_xx = cbuffer.data(gd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxz_xy = cbuffer.data(gd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxz_xz = cbuffer.data(gd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxz_yy = cbuffer.data(gd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxz_yz = cbuffer.data(gd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxz_zz = cbuffer.data(gd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxyy_xx = cbuffer.data(gd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxyy_xy = cbuffer.data(gd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxyy_xz = cbuffer.data(gd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxyy_yy = cbuffer.data(gd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxyy_yz = cbuffer.data(gd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxyy_zz = cbuffer.data(gd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxyz_xx = cbuffer.data(gd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxyz_xy = cbuffer.data(gd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxyz_xz = cbuffer.data(gd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxyz_yy = cbuffer.data(gd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxyz_yz = cbuffer.data(gd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxyz_zz = cbuffer.data(gd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxzz_xx = cbuffer.data(gd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxzz_xy = cbuffer.data(gd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxzz_xz = cbuffer.data(gd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxzz_yy = cbuffer.data(gd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxzz_yz = cbuffer.data(gd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxzz_zz = cbuffer.data(gd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xyyy_xx = cbuffer.data(gd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyyy_xy = cbuffer.data(gd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyyy_xz = cbuffer.data(gd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyyy_yy = cbuffer.data(gd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyyy_yz = cbuffer.data(gd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyyy_zz = cbuffer.data(gd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xyyz_xx = cbuffer.data(gd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xyyz_xy = cbuffer.data(gd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xyyz_xz = cbuffer.data(gd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyyz_yy = cbuffer.data(gd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyyz_yz = cbuffer.data(gd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyyz_zz = cbuffer.data(gd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyzz_xx = cbuffer.data(gd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyzz_xy = cbuffer.data(gd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyzz_xz = cbuffer.data(gd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyzz_yy = cbuffer.data(gd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyzz_yz = cbuffer.data(gd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyzz_zz = cbuffer.data(gd_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xzzz_xx = cbuffer.data(gd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xzzz_xy = cbuffer.data(gd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xzzz_xz = cbuffer.data(gd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xzzz_yy = cbuffer.data(gd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xzzz_yz = cbuffer.data(gd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xzzz_zz = cbuffer.data(gd_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yyyy_xx = cbuffer.data(gd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yyyy_xy = cbuffer.data(gd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yyyy_xz = cbuffer.data(gd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyyy_yy = cbuffer.data(gd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyyy_yz = cbuffer.data(gd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyyy_zz = cbuffer.data(gd_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yyyz_xx = cbuffer.data(gd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyyz_xy = cbuffer.data(gd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyyz_xz = cbuffer.data(gd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyyz_yy = cbuffer.data(gd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyyz_yz = cbuffer.data(gd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyyz_zz = cbuffer.data(gd_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yyzz_xx = cbuffer.data(gd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyzz_xy = cbuffer.data(gd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyzz_xz = cbuffer.data(gd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyzz_yy = cbuffer.data(gd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyzz_yz = cbuffer.data(gd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyzz_zz = cbuffer.data(gd_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yzzz_xx = cbuffer.data(gd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yzzz_xy = cbuffer.data(gd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzzz_xz = cbuffer.data(gd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yzzz_yy = cbuffer.data(gd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yzzz_yz = cbuffer.data(gd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yzzz_zz = cbuffer.data(gd_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_zzzz_xx = cbuffer.data(gd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zzzz_xy = cbuffer.data(gd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zzzz_xz = cbuffer.data(gd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zzzz_yy = cbuffer.data(gd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zzzz_yz = cbuffer.data(gd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zzzz_zz = cbuffer.data(gd_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xxxx_xx = cbuffer.data(gd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxxx_xy = cbuffer.data(gd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxxx_xz = cbuffer.data(gd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xxxx_yy = cbuffer.data(gd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xxxx_yz = cbuffer.data(gd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xxxx_zz = cbuffer.data(gd_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xxxy_xx = cbuffer.data(gd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xxxy_xy = cbuffer.data(gd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xxxy_xz = cbuffer.data(gd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xxxy_yy = cbuffer.data(gd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxxy_yz = cbuffer.data(gd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxxy_zz = cbuffer.data(gd_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xxxz_xx = cbuffer.data(gd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxxz_xy = cbuffer.data(gd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxxz_xz = cbuffer.data(gd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxxz_yy = cbuffer.data(gd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxxz_yz = cbuffer.data(gd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxxz_zz = cbuffer.data(gd_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xxyy_xx = cbuffer.data(gd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxyy_xy = cbuffer.data(gd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxyy_xz = cbuffer.data(gd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxyy_yy = cbuffer.data(gd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxyy_yz = cbuffer.data(gd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxyy_zz = cbuffer.data(gd_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xxyz_xx = cbuffer.data(gd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxyz_xy = cbuffer.data(gd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxyz_xz = cbuffer.data(gd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxyz_yy = cbuffer.data(gd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxyz_yz = cbuffer.data(gd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxyz_zz = cbuffer.data(gd_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xxzz_xx = cbuffer.data(gd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxzz_xy = cbuffer.data(gd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxzz_xz = cbuffer.data(gd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxzz_yy = cbuffer.data(gd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxzz_yz = cbuffer.data(gd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxzz_zz = cbuffer.data(gd_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xyyy_xx = cbuffer.data(gd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xyyy_xy = cbuffer.data(gd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xyyy_xz = cbuffer.data(gd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xyyy_yy = cbuffer.data(gd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyyy_yz = cbuffer.data(gd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyyy_zz = cbuffer.data(gd_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xyyz_xx = cbuffer.data(gd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyyz_xy = cbuffer.data(gd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyyz_xz = cbuffer.data(gd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyyz_yy = cbuffer.data(gd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyyz_yz = cbuffer.data(gd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyyz_zz = cbuffer.data(gd_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xyzz_xx = cbuffer.data(gd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyzz_xy = cbuffer.data(gd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyzz_xz = cbuffer.data(gd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyzz_yy = cbuffer.data(gd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyzz_yz = cbuffer.data(gd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyzz_zz = cbuffer.data(gd_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xzzz_xx = cbuffer.data(gd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xzzz_xy = cbuffer.data(gd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xzzz_xz = cbuffer.data(gd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xzzz_yy = cbuffer.data(gd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xzzz_yz = cbuffer.data(gd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xzzz_zz = cbuffer.data(gd_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yyyy_xx = cbuffer.data(gd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yyyy_xy = cbuffer.data(gd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yyyy_xz = cbuffer.data(gd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yyyy_yy = cbuffer.data(gd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yyyy_yz = cbuffer.data(gd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yyyy_zz = cbuffer.data(gd_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yyyz_xx = cbuffer.data(gd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yyyz_xy = cbuffer.data(gd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yyyz_xz = cbuffer.data(gd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yyyz_yy = cbuffer.data(gd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyyz_yz = cbuffer.data(gd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyyz_zz = cbuffer.data(gd_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yyzz_xx = cbuffer.data(gd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yyzz_xy = cbuffer.data(gd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yyzz_xz = cbuffer.data(gd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_yyzz_yy = cbuffer.data(gd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_yyzz_yz = cbuffer.data(gd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_yyzz_zz = cbuffer.data(gd_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_yzzz_xx = cbuffer.data(gd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_yzzz_xy = cbuffer.data(gd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_yzzz_xz = cbuffer.data(gd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_yzzz_yy = cbuffer.data(gd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_yzzz_yz = cbuffer.data(gd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_yzzz_zz = cbuffer.data(gd_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_zzzz_xx = cbuffer.data(gd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zzzz_xy = cbuffer.data(gd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zzzz_xz = cbuffer.data(gd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zzzz_yy = cbuffer.data(gd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zzzz_yz = cbuffer.data(gd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zzzz_zz = cbuffer.data(gd_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xxxx_xx = cbuffer.data(gd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xxxx_xy = cbuffer.data(gd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xxxx_xz = cbuffer.data(gd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xxxx_yy = cbuffer.data(gd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xxxx_yz = cbuffer.data(gd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xxxx_zz = cbuffer.data(gd_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xxxy_xx = cbuffer.data(gd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xxxy_xy = cbuffer.data(gd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xxxy_xz = cbuffer.data(gd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xxxy_yy = cbuffer.data(gd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xxxy_yz = cbuffer.data(gd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xxxy_zz = cbuffer.data(gd_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xxxz_xx = cbuffer.data(gd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xxxz_xy = cbuffer.data(gd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xxxz_xz = cbuffer.data(gd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xxxz_yy = cbuffer.data(gd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xxxz_yz = cbuffer.data(gd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xxxz_zz = cbuffer.data(gd_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xxyy_xx = cbuffer.data(gd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xxyy_xy = cbuffer.data(gd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxyy_xz = cbuffer.data(gd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxyy_yy = cbuffer.data(gd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxyy_yz = cbuffer.data(gd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxyy_zz = cbuffer.data(gd_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xxyz_xx = cbuffer.data(gd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxyz_xy = cbuffer.data(gd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxyz_xz = cbuffer.data(gd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxyz_yy = cbuffer.data(gd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxyz_yz = cbuffer.data(gd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxyz_zz = cbuffer.data(gd_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xxzz_xx = cbuffer.data(gd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxzz_xy = cbuffer.data(gd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxzz_xz = cbuffer.data(gd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xxzz_yy = cbuffer.data(gd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xxzz_yz = cbuffer.data(gd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xxzz_zz = cbuffer.data(gd_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xyyy_xx = cbuffer.data(gd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xyyy_xy = cbuffer.data(gd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xyyy_xz = cbuffer.data(gd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xyyy_yy = cbuffer.data(gd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xyyy_yz = cbuffer.data(gd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xyyy_zz = cbuffer.data(gd_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xyyz_xx = cbuffer.data(gd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xyyz_xy = cbuffer.data(gd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xyyz_xz = cbuffer.data(gd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xyyz_yy = cbuffer.data(gd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xyyz_yz = cbuffer.data(gd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xyyz_zz = cbuffer.data(gd_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_xyzz_xx = cbuffer.data(gd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xyzz_xy = cbuffer.data(gd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xyzz_xz = cbuffer.data(gd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_xyzz_yy = cbuffer.data(gd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xyzz_yz = cbuffer.data(gd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xyzz_zz = cbuffer.data(gd_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_xzzz_xx = cbuffer.data(gd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xzzz_xy = cbuffer.data(gd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xzzz_xz = cbuffer.data(gd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_xzzz_yy = cbuffer.data(gd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xzzz_yz = cbuffer.data(gd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xzzz_zz = cbuffer.data(gd_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yyyy_xx = cbuffer.data(gd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yyyy_xy = cbuffer.data(gd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yyyy_xz = cbuffer.data(gd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yyyy_yy = cbuffer.data(gd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yyyy_yz = cbuffer.data(gd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yyyy_zz = cbuffer.data(gd_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yyyz_xx = cbuffer.data(gd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yyyz_xy = cbuffer.data(gd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yyyz_xz = cbuffer.data(gd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yyyz_yy = cbuffer.data(gd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yyyz_yz = cbuffer.data(gd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yyyz_zz = cbuffer.data(gd_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_yyzz_xx = cbuffer.data(gd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yyzz_xy = cbuffer.data(gd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yyzz_xz = cbuffer.data(gd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_yyzz_yy = cbuffer.data(gd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_yyzz_yz = cbuffer.data(gd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_yyzz_zz = cbuffer.data(gd_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_yzzz_xx = cbuffer.data(gd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_yzzz_xy = cbuffer.data(gd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_yzzz_xz = cbuffer.data(gd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_yzzz_yy = cbuffer.data(gd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_yzzz_yz = cbuffer.data(gd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_yzzz_zz = cbuffer.data(gd_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_zzzz_xx = cbuffer.data(gd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zzzz_xy = cbuffer.data(gd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zzzz_xz = cbuffer.data(gd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zzzz_yy = cbuffer.data(gd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zzzz_yz = cbuffer.data(gd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zzzz_zz = cbuffer.data(gd_geom_01_off + 269 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hpxx

            const auto hp_geom_01_off = idx_geom_01_hpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxx_x = cbuffer.data(hp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_y = cbuffer.data(hp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_z = cbuffer.data(hp_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_x, g_0_x_xxxx_xx, g_0_x_xxxx_xy, g_0_x_xxxx_xz, g_0_x_xxxx_y, g_0_x_xxxx_z, g_0_x_xxxxx_x, g_0_x_xxxxx_y, g_0_x_xxxxx_z, g_xxxx_x, g_xxxx_y, g_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxx_x[k] = g_xxxx_x[k] - g_0_x_xxxx_x[k] * ab_x + g_0_x_xxxx_xx[k];

                g_0_x_xxxxx_y[k] = g_xxxx_y[k] - g_0_x_xxxx_y[k] * ab_x + g_0_x_xxxx_xy[k];

                g_0_x_xxxxx_z[k] = g_xxxx_z[k] - g_0_x_xxxx_z[k] * ab_x + g_0_x_xxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxy_x = cbuffer.data(hp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxy_y = cbuffer.data(hp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxy_z = cbuffer.data(hp_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_x, g_0_x_xxxx_xy, g_0_x_xxxx_y, g_0_x_xxxx_yy, g_0_x_xxxx_yz, g_0_x_xxxx_z, g_0_x_xxxxy_x, g_0_x_xxxxy_y, g_0_x_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxy_x[k] = -g_0_x_xxxx_x[k] * ab_y + g_0_x_xxxx_xy[k];

                g_0_x_xxxxy_y[k] = -g_0_x_xxxx_y[k] * ab_y + g_0_x_xxxx_yy[k];

                g_0_x_xxxxy_z[k] = -g_0_x_xxxx_z[k] * ab_y + g_0_x_xxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxz_x = cbuffer.data(hp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxz_y = cbuffer.data(hp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxz_z = cbuffer.data(hp_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_x, g_0_x_xxxx_xz, g_0_x_xxxx_y, g_0_x_xxxx_yz, g_0_x_xxxx_z, g_0_x_xxxx_zz, g_0_x_xxxxz_x, g_0_x_xxxxz_y, g_0_x_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxz_x[k] = -g_0_x_xxxx_x[k] * ab_z + g_0_x_xxxx_xz[k];

                g_0_x_xxxxz_y[k] = -g_0_x_xxxx_y[k] * ab_z + g_0_x_xxxx_yz[k];

                g_0_x_xxxxz_z[k] = -g_0_x_xxxx_z[k] * ab_z + g_0_x_xxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyy_x = cbuffer.data(hp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxyy_y = cbuffer.data(hp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxyy_z = cbuffer.data(hp_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxy_x, g_0_x_xxxy_xy, g_0_x_xxxy_y, g_0_x_xxxy_yy, g_0_x_xxxy_yz, g_0_x_xxxy_z, g_0_x_xxxyy_x, g_0_x_xxxyy_y, g_0_x_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyy_x[k] = -g_0_x_xxxy_x[k] * ab_y + g_0_x_xxxy_xy[k];

                g_0_x_xxxyy_y[k] = -g_0_x_xxxy_y[k] * ab_y + g_0_x_xxxy_yy[k];

                g_0_x_xxxyy_z[k] = -g_0_x_xxxy_z[k] * ab_y + g_0_x_xxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyz_x = cbuffer.data(hp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxyz_y = cbuffer.data(hp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxyz_z = cbuffer.data(hp_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyz_x, g_0_x_xxxyz_y, g_0_x_xxxyz_z, g_0_x_xxxz_x, g_0_x_xxxz_xy, g_0_x_xxxz_y, g_0_x_xxxz_yy, g_0_x_xxxz_yz, g_0_x_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyz_x[k] = -g_0_x_xxxz_x[k] * ab_y + g_0_x_xxxz_xy[k];

                g_0_x_xxxyz_y[k] = -g_0_x_xxxz_y[k] * ab_y + g_0_x_xxxz_yy[k];

                g_0_x_xxxyz_z[k] = -g_0_x_xxxz_z[k] * ab_y + g_0_x_xxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzz_x = cbuffer.data(hp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxzz_y = cbuffer.data(hp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxzz_z = cbuffer.data(hp_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxz_x, g_0_x_xxxz_xz, g_0_x_xxxz_y, g_0_x_xxxz_yz, g_0_x_xxxz_z, g_0_x_xxxz_zz, g_0_x_xxxzz_x, g_0_x_xxxzz_y, g_0_x_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzz_x[k] = -g_0_x_xxxz_x[k] * ab_z + g_0_x_xxxz_xz[k];

                g_0_x_xxxzz_y[k] = -g_0_x_xxxz_y[k] * ab_z + g_0_x_xxxz_yz[k];

                g_0_x_xxxzz_z[k] = -g_0_x_xxxz_z[k] * ab_z + g_0_x_xxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyy_x = cbuffer.data(hp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxyyy_y = cbuffer.data(hp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxyyy_z = cbuffer.data(hp_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyy_x, g_0_x_xxyy_xy, g_0_x_xxyy_y, g_0_x_xxyy_yy, g_0_x_xxyy_yz, g_0_x_xxyy_z, g_0_x_xxyyy_x, g_0_x_xxyyy_y, g_0_x_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyy_x[k] = -g_0_x_xxyy_x[k] * ab_y + g_0_x_xxyy_xy[k];

                g_0_x_xxyyy_y[k] = -g_0_x_xxyy_y[k] * ab_y + g_0_x_xxyy_yy[k];

                g_0_x_xxyyy_z[k] = -g_0_x_xxyy_z[k] * ab_y + g_0_x_xxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyz_x = cbuffer.data(hp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxyyz_y = cbuffer.data(hp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxyyz_z = cbuffer.data(hp_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyz_x, g_0_x_xxyyz_y, g_0_x_xxyyz_z, g_0_x_xxyz_x, g_0_x_xxyz_xy, g_0_x_xxyz_y, g_0_x_xxyz_yy, g_0_x_xxyz_yz, g_0_x_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyz_x[k] = -g_0_x_xxyz_x[k] * ab_y + g_0_x_xxyz_xy[k];

                g_0_x_xxyyz_y[k] = -g_0_x_xxyz_y[k] * ab_y + g_0_x_xxyz_yy[k];

                g_0_x_xxyyz_z[k] = -g_0_x_xxyz_z[k] * ab_y + g_0_x_xxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzz_x = cbuffer.data(hp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxyzz_y = cbuffer.data(hp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxyzz_z = cbuffer.data(hp_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzz_x, g_0_x_xxyzz_y, g_0_x_xxyzz_z, g_0_x_xxzz_x, g_0_x_xxzz_xy, g_0_x_xxzz_y, g_0_x_xxzz_yy, g_0_x_xxzz_yz, g_0_x_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzz_x[k] = -g_0_x_xxzz_x[k] * ab_y + g_0_x_xxzz_xy[k];

                g_0_x_xxyzz_y[k] = -g_0_x_xxzz_y[k] * ab_y + g_0_x_xxzz_yy[k];

                g_0_x_xxyzz_z[k] = -g_0_x_xxzz_z[k] * ab_y + g_0_x_xxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzz_x = cbuffer.data(hp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxzzz_y = cbuffer.data(hp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxzzz_z = cbuffer.data(hp_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzz_x, g_0_x_xxzz_xz, g_0_x_xxzz_y, g_0_x_xxzz_yz, g_0_x_xxzz_z, g_0_x_xxzz_zz, g_0_x_xxzzz_x, g_0_x_xxzzz_y, g_0_x_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzz_x[k] = -g_0_x_xxzz_x[k] * ab_z + g_0_x_xxzz_xz[k];

                g_0_x_xxzzz_y[k] = -g_0_x_xxzz_y[k] * ab_z + g_0_x_xxzz_yz[k];

                g_0_x_xxzzz_z[k] = -g_0_x_xxzz_z[k] * ab_z + g_0_x_xxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyy_x = cbuffer.data(hp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xyyyy_y = cbuffer.data(hp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xyyyy_z = cbuffer.data(hp_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyy_x, g_0_x_xyyy_xy, g_0_x_xyyy_y, g_0_x_xyyy_yy, g_0_x_xyyy_yz, g_0_x_xyyy_z, g_0_x_xyyyy_x, g_0_x_xyyyy_y, g_0_x_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyy_x[k] = -g_0_x_xyyy_x[k] * ab_y + g_0_x_xyyy_xy[k];

                g_0_x_xyyyy_y[k] = -g_0_x_xyyy_y[k] * ab_y + g_0_x_xyyy_yy[k];

                g_0_x_xyyyy_z[k] = -g_0_x_xyyy_z[k] * ab_y + g_0_x_xyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyz_x = cbuffer.data(hp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xyyyz_y = cbuffer.data(hp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xyyyz_z = cbuffer.data(hp_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyz_x, g_0_x_xyyyz_y, g_0_x_xyyyz_z, g_0_x_xyyz_x, g_0_x_xyyz_xy, g_0_x_xyyz_y, g_0_x_xyyz_yy, g_0_x_xyyz_yz, g_0_x_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyz_x[k] = -g_0_x_xyyz_x[k] * ab_y + g_0_x_xyyz_xy[k];

                g_0_x_xyyyz_y[k] = -g_0_x_xyyz_y[k] * ab_y + g_0_x_xyyz_yy[k];

                g_0_x_xyyyz_z[k] = -g_0_x_xyyz_z[k] * ab_y + g_0_x_xyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzz_x = cbuffer.data(hp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyyzz_y = cbuffer.data(hp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyyzz_z = cbuffer.data(hp_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzz_x, g_0_x_xyyzz_y, g_0_x_xyyzz_z, g_0_x_xyzz_x, g_0_x_xyzz_xy, g_0_x_xyzz_y, g_0_x_xyzz_yy, g_0_x_xyzz_yz, g_0_x_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzz_x[k] = -g_0_x_xyzz_x[k] * ab_y + g_0_x_xyzz_xy[k];

                g_0_x_xyyzz_y[k] = -g_0_x_xyzz_y[k] * ab_y + g_0_x_xyzz_yy[k];

                g_0_x_xyyzz_z[k] = -g_0_x_xyzz_z[k] * ab_y + g_0_x_xyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzz_x = cbuffer.data(hp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyzzz_y = cbuffer.data(hp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyzzz_z = cbuffer.data(hp_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzz_x, g_0_x_xyzzz_y, g_0_x_xyzzz_z, g_0_x_xzzz_x, g_0_x_xzzz_xy, g_0_x_xzzz_y, g_0_x_xzzz_yy, g_0_x_xzzz_yz, g_0_x_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzz_x[k] = -g_0_x_xzzz_x[k] * ab_y + g_0_x_xzzz_xy[k];

                g_0_x_xyzzz_y[k] = -g_0_x_xzzz_y[k] * ab_y + g_0_x_xzzz_yy[k];

                g_0_x_xyzzz_z[k] = -g_0_x_xzzz_z[k] * ab_y + g_0_x_xzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzz_x = cbuffer.data(hp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xzzzz_y = cbuffer.data(hp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xzzzz_z = cbuffer.data(hp_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzz_x, g_0_x_xzzz_xz, g_0_x_xzzz_y, g_0_x_xzzz_yz, g_0_x_xzzz_z, g_0_x_xzzz_zz, g_0_x_xzzzz_x, g_0_x_xzzzz_y, g_0_x_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzz_x[k] = -g_0_x_xzzz_x[k] * ab_z + g_0_x_xzzz_xz[k];

                g_0_x_xzzzz_y[k] = -g_0_x_xzzz_y[k] * ab_z + g_0_x_xzzz_yz[k];

                g_0_x_xzzzz_z[k] = -g_0_x_xzzz_z[k] * ab_z + g_0_x_xzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyy_x = cbuffer.data(hp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yyyyy_y = cbuffer.data(hp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yyyyy_z = cbuffer.data(hp_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_x, g_0_x_yyyy_xy, g_0_x_yyyy_y, g_0_x_yyyy_yy, g_0_x_yyyy_yz, g_0_x_yyyy_z, g_0_x_yyyyy_x, g_0_x_yyyyy_y, g_0_x_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyy_x[k] = -g_0_x_yyyy_x[k] * ab_y + g_0_x_yyyy_xy[k];

                g_0_x_yyyyy_y[k] = -g_0_x_yyyy_y[k] * ab_y + g_0_x_yyyy_yy[k];

                g_0_x_yyyyy_z[k] = -g_0_x_yyyy_z[k] * ab_y + g_0_x_yyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyz_x = cbuffer.data(hp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yyyyz_y = cbuffer.data(hp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yyyyz_z = cbuffer.data(hp_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyz_x, g_0_x_yyyyz_y, g_0_x_yyyyz_z, g_0_x_yyyz_x, g_0_x_yyyz_xy, g_0_x_yyyz_y, g_0_x_yyyz_yy, g_0_x_yyyz_yz, g_0_x_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyz_x[k] = -g_0_x_yyyz_x[k] * ab_y + g_0_x_yyyz_xy[k];

                g_0_x_yyyyz_y[k] = -g_0_x_yyyz_y[k] * ab_y + g_0_x_yyyz_yy[k];

                g_0_x_yyyyz_z[k] = -g_0_x_yyyz_z[k] * ab_y + g_0_x_yyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzz_x = cbuffer.data(hp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yyyzz_y = cbuffer.data(hp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yyyzz_z = cbuffer.data(hp_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzz_x, g_0_x_yyyzz_y, g_0_x_yyyzz_z, g_0_x_yyzz_x, g_0_x_yyzz_xy, g_0_x_yyzz_y, g_0_x_yyzz_yy, g_0_x_yyzz_yz, g_0_x_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzz_x[k] = -g_0_x_yyzz_x[k] * ab_y + g_0_x_yyzz_xy[k];

                g_0_x_yyyzz_y[k] = -g_0_x_yyzz_y[k] * ab_y + g_0_x_yyzz_yy[k];

                g_0_x_yyyzz_z[k] = -g_0_x_yyzz_z[k] * ab_y + g_0_x_yyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzz_x = cbuffer.data(hp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_yyzzz_y = cbuffer.data(hp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_yyzzz_z = cbuffer.data(hp_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzz_x, g_0_x_yyzzz_y, g_0_x_yyzzz_z, g_0_x_yzzz_x, g_0_x_yzzz_xy, g_0_x_yzzz_y, g_0_x_yzzz_yy, g_0_x_yzzz_yz, g_0_x_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzz_x[k] = -g_0_x_yzzz_x[k] * ab_y + g_0_x_yzzz_xy[k];

                g_0_x_yyzzz_y[k] = -g_0_x_yzzz_y[k] * ab_y + g_0_x_yzzz_yy[k];

                g_0_x_yyzzz_z[k] = -g_0_x_yzzz_z[k] * ab_y + g_0_x_yzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzz_x = cbuffer.data(hp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_yzzzz_y = cbuffer.data(hp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_yzzzz_z = cbuffer.data(hp_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzz_x, g_0_x_yzzzz_y, g_0_x_yzzzz_z, g_0_x_zzzz_x, g_0_x_zzzz_xy, g_0_x_zzzz_y, g_0_x_zzzz_yy, g_0_x_zzzz_yz, g_0_x_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzz_x[k] = -g_0_x_zzzz_x[k] * ab_y + g_0_x_zzzz_xy[k];

                g_0_x_yzzzz_y[k] = -g_0_x_zzzz_y[k] * ab_y + g_0_x_zzzz_yy[k];

                g_0_x_yzzzz_z[k] = -g_0_x_zzzz_z[k] * ab_y + g_0_x_zzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzz_x = cbuffer.data(hp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_zzzzz_y = cbuffer.data(hp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_zzzzz_z = cbuffer.data(hp_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_x, g_0_x_zzzz_xz, g_0_x_zzzz_y, g_0_x_zzzz_yz, g_0_x_zzzz_z, g_0_x_zzzz_zz, g_0_x_zzzzz_x, g_0_x_zzzzz_y, g_0_x_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzz_x[k] = -g_0_x_zzzz_x[k] * ab_z + g_0_x_zzzz_xz[k];

                g_0_x_zzzzz_y[k] = -g_0_x_zzzz_y[k] * ab_z + g_0_x_zzzz_yz[k];

                g_0_x_zzzzz_z[k] = -g_0_x_zzzz_z[k] * ab_z + g_0_x_zzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxx_x = cbuffer.data(hp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xxxxx_y = cbuffer.data(hp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xxxxx_z = cbuffer.data(hp_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_x, g_0_y_xxxx_xx, g_0_y_xxxx_xy, g_0_y_xxxx_xz, g_0_y_xxxx_y, g_0_y_xxxx_z, g_0_y_xxxxx_x, g_0_y_xxxxx_y, g_0_y_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxx_x[k] = -g_0_y_xxxx_x[k] * ab_x + g_0_y_xxxx_xx[k];

                g_0_y_xxxxx_y[k] = -g_0_y_xxxx_y[k] * ab_x + g_0_y_xxxx_xy[k];

                g_0_y_xxxxx_z[k] = -g_0_y_xxxx_z[k] * ab_x + g_0_y_xxxx_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxy_x = cbuffer.data(hp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xxxxy_y = cbuffer.data(hp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xxxxy_z = cbuffer.data(hp_geom_01_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxy_x, g_0_y_xxxxy_y, g_0_y_xxxxy_z, g_0_y_xxxy_x, g_0_y_xxxy_xx, g_0_y_xxxy_xy, g_0_y_xxxy_xz, g_0_y_xxxy_y, g_0_y_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxy_x[k] = -g_0_y_xxxy_x[k] * ab_x + g_0_y_xxxy_xx[k];

                g_0_y_xxxxy_y[k] = -g_0_y_xxxy_y[k] * ab_x + g_0_y_xxxy_xy[k];

                g_0_y_xxxxy_z[k] = -g_0_y_xxxy_z[k] * ab_x + g_0_y_xxxy_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxz_x = cbuffer.data(hp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xxxxz_y = cbuffer.data(hp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xxxxz_z = cbuffer.data(hp_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxz_x, g_0_y_xxxxz_y, g_0_y_xxxxz_z, g_0_y_xxxz_x, g_0_y_xxxz_xx, g_0_y_xxxz_xy, g_0_y_xxxz_xz, g_0_y_xxxz_y, g_0_y_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxz_x[k] = -g_0_y_xxxz_x[k] * ab_x + g_0_y_xxxz_xx[k];

                g_0_y_xxxxz_y[k] = -g_0_y_xxxz_y[k] * ab_x + g_0_y_xxxz_xy[k];

                g_0_y_xxxxz_z[k] = -g_0_y_xxxz_z[k] * ab_x + g_0_y_xxxz_xz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyy_x = cbuffer.data(hp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xxxyy_y = cbuffer.data(hp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xxxyy_z = cbuffer.data(hp_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyy_x, g_0_y_xxxyy_y, g_0_y_xxxyy_z, g_0_y_xxyy_x, g_0_y_xxyy_xx, g_0_y_xxyy_xy, g_0_y_xxyy_xz, g_0_y_xxyy_y, g_0_y_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyy_x[k] = -g_0_y_xxyy_x[k] * ab_x + g_0_y_xxyy_xx[k];

                g_0_y_xxxyy_y[k] = -g_0_y_xxyy_y[k] * ab_x + g_0_y_xxyy_xy[k];

                g_0_y_xxxyy_z[k] = -g_0_y_xxyy_z[k] * ab_x + g_0_y_xxyy_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyz_x = cbuffer.data(hp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_xxxyz_y = cbuffer.data(hp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_xxxyz_z = cbuffer.data(hp_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyz_x, g_0_y_xxxyz_y, g_0_y_xxxyz_z, g_0_y_xxyz_x, g_0_y_xxyz_xx, g_0_y_xxyz_xy, g_0_y_xxyz_xz, g_0_y_xxyz_y, g_0_y_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyz_x[k] = -g_0_y_xxyz_x[k] * ab_x + g_0_y_xxyz_xx[k];

                g_0_y_xxxyz_y[k] = -g_0_y_xxyz_y[k] * ab_x + g_0_y_xxyz_xy[k];

                g_0_y_xxxyz_z[k] = -g_0_y_xxyz_z[k] * ab_x + g_0_y_xxyz_xz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzz_x = cbuffer.data(hp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_xxxzz_y = cbuffer.data(hp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_xxxzz_z = cbuffer.data(hp_geom_01_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzz_x, g_0_y_xxxzz_y, g_0_y_xxxzz_z, g_0_y_xxzz_x, g_0_y_xxzz_xx, g_0_y_xxzz_xy, g_0_y_xxzz_xz, g_0_y_xxzz_y, g_0_y_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzz_x[k] = -g_0_y_xxzz_x[k] * ab_x + g_0_y_xxzz_xx[k];

                g_0_y_xxxzz_y[k] = -g_0_y_xxzz_y[k] * ab_x + g_0_y_xxzz_xy[k];

                g_0_y_xxxzz_z[k] = -g_0_y_xxzz_z[k] * ab_x + g_0_y_xxzz_xz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyy_x = cbuffer.data(hp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_xxyyy_y = cbuffer.data(hp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_xxyyy_z = cbuffer.data(hp_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyy_x, g_0_y_xxyyy_y, g_0_y_xxyyy_z, g_0_y_xyyy_x, g_0_y_xyyy_xx, g_0_y_xyyy_xy, g_0_y_xyyy_xz, g_0_y_xyyy_y, g_0_y_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyy_x[k] = -g_0_y_xyyy_x[k] * ab_x + g_0_y_xyyy_xx[k];

                g_0_y_xxyyy_y[k] = -g_0_y_xyyy_y[k] * ab_x + g_0_y_xyyy_xy[k];

                g_0_y_xxyyy_z[k] = -g_0_y_xyyy_z[k] * ab_x + g_0_y_xyyy_xz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyz_x = cbuffer.data(hp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xxyyz_y = cbuffer.data(hp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xxyyz_z = cbuffer.data(hp_geom_01_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyz_x, g_0_y_xxyyz_y, g_0_y_xxyyz_z, g_0_y_xyyz_x, g_0_y_xyyz_xx, g_0_y_xyyz_xy, g_0_y_xyyz_xz, g_0_y_xyyz_y, g_0_y_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyz_x[k] = -g_0_y_xyyz_x[k] * ab_x + g_0_y_xyyz_xx[k];

                g_0_y_xxyyz_y[k] = -g_0_y_xyyz_y[k] * ab_x + g_0_y_xyyz_xy[k];

                g_0_y_xxyyz_z[k] = -g_0_y_xyyz_z[k] * ab_x + g_0_y_xyyz_xz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzz_x = cbuffer.data(hp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xxyzz_y = cbuffer.data(hp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xxyzz_z = cbuffer.data(hp_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzz_x, g_0_y_xxyzz_y, g_0_y_xxyzz_z, g_0_y_xyzz_x, g_0_y_xyzz_xx, g_0_y_xyzz_xy, g_0_y_xyzz_xz, g_0_y_xyzz_y, g_0_y_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzz_x[k] = -g_0_y_xyzz_x[k] * ab_x + g_0_y_xyzz_xx[k];

                g_0_y_xxyzz_y[k] = -g_0_y_xyzz_y[k] * ab_x + g_0_y_xyzz_xy[k];

                g_0_y_xxyzz_z[k] = -g_0_y_xyzz_z[k] * ab_x + g_0_y_xyzz_xz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzz_x = cbuffer.data(hp_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxzzz_y = cbuffer.data(hp_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxzzz_z = cbuffer.data(hp_geom_01_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzz_x, g_0_y_xxzzz_y, g_0_y_xxzzz_z, g_0_y_xzzz_x, g_0_y_xzzz_xx, g_0_y_xzzz_xy, g_0_y_xzzz_xz, g_0_y_xzzz_y, g_0_y_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzz_x[k] = -g_0_y_xzzz_x[k] * ab_x + g_0_y_xzzz_xx[k];

                g_0_y_xxzzz_y[k] = -g_0_y_xzzz_y[k] * ab_x + g_0_y_xzzz_xy[k];

                g_0_y_xxzzz_z[k] = -g_0_y_xzzz_z[k] * ab_x + g_0_y_xzzz_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyy_x = cbuffer.data(hp_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xyyyy_y = cbuffer.data(hp_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xyyyy_z = cbuffer.data(hp_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyy_x, g_0_y_xyyyy_y, g_0_y_xyyyy_z, g_0_y_yyyy_x, g_0_y_yyyy_xx, g_0_y_yyyy_xy, g_0_y_yyyy_xz, g_0_y_yyyy_y, g_0_y_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyy_x[k] = -g_0_y_yyyy_x[k] * ab_x + g_0_y_yyyy_xx[k];

                g_0_y_xyyyy_y[k] = -g_0_y_yyyy_y[k] * ab_x + g_0_y_yyyy_xy[k];

                g_0_y_xyyyy_z[k] = -g_0_y_yyyy_z[k] * ab_x + g_0_y_yyyy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyz_x = cbuffer.data(hp_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xyyyz_y = cbuffer.data(hp_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xyyyz_z = cbuffer.data(hp_geom_01_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyz_x, g_0_y_xyyyz_y, g_0_y_xyyyz_z, g_0_y_yyyz_x, g_0_y_yyyz_xx, g_0_y_yyyz_xy, g_0_y_yyyz_xz, g_0_y_yyyz_y, g_0_y_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyz_x[k] = -g_0_y_yyyz_x[k] * ab_x + g_0_y_yyyz_xx[k];

                g_0_y_xyyyz_y[k] = -g_0_y_yyyz_y[k] * ab_x + g_0_y_yyyz_xy[k];

                g_0_y_xyyyz_z[k] = -g_0_y_yyyz_z[k] * ab_x + g_0_y_yyyz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzz_x = cbuffer.data(hp_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xyyzz_y = cbuffer.data(hp_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xyyzz_z = cbuffer.data(hp_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzz_x, g_0_y_xyyzz_y, g_0_y_xyyzz_z, g_0_y_yyzz_x, g_0_y_yyzz_xx, g_0_y_yyzz_xy, g_0_y_yyzz_xz, g_0_y_yyzz_y, g_0_y_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzz_x[k] = -g_0_y_yyzz_x[k] * ab_x + g_0_y_yyzz_xx[k];

                g_0_y_xyyzz_y[k] = -g_0_y_yyzz_y[k] * ab_x + g_0_y_yyzz_xy[k];

                g_0_y_xyyzz_z[k] = -g_0_y_yyzz_z[k] * ab_x + g_0_y_yyzz_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzz_x = cbuffer.data(hp_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xyzzz_y = cbuffer.data(hp_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xyzzz_z = cbuffer.data(hp_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzz_x, g_0_y_xyzzz_y, g_0_y_xyzzz_z, g_0_y_yzzz_x, g_0_y_yzzz_xx, g_0_y_yzzz_xy, g_0_y_yzzz_xz, g_0_y_yzzz_y, g_0_y_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzz_x[k] = -g_0_y_yzzz_x[k] * ab_x + g_0_y_yzzz_xx[k];

                g_0_y_xyzzz_y[k] = -g_0_y_yzzz_y[k] * ab_x + g_0_y_yzzz_xy[k];

                g_0_y_xyzzz_z[k] = -g_0_y_yzzz_z[k] * ab_x + g_0_y_yzzz_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzz_x = cbuffer.data(hp_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xzzzz_y = cbuffer.data(hp_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xzzzz_z = cbuffer.data(hp_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzz_x, g_0_y_xzzzz_y, g_0_y_xzzzz_z, g_0_y_zzzz_x, g_0_y_zzzz_xx, g_0_y_zzzz_xy, g_0_y_zzzz_xz, g_0_y_zzzz_y, g_0_y_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzz_x[k] = -g_0_y_zzzz_x[k] * ab_x + g_0_y_zzzz_xx[k];

                g_0_y_xzzzz_y[k] = -g_0_y_zzzz_y[k] * ab_x + g_0_y_zzzz_xy[k];

                g_0_y_xzzzz_z[k] = -g_0_y_zzzz_z[k] * ab_x + g_0_y_zzzz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyy_x = cbuffer.data(hp_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_yyyyy_y = cbuffer.data(hp_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_yyyyy_z = cbuffer.data(hp_geom_01_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_x, g_0_y_yyyy_xy, g_0_y_yyyy_y, g_0_y_yyyy_yy, g_0_y_yyyy_yz, g_0_y_yyyy_z, g_0_y_yyyyy_x, g_0_y_yyyyy_y, g_0_y_yyyyy_z, g_yyyy_x, g_yyyy_y, g_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyy_x[k] = g_yyyy_x[k] - g_0_y_yyyy_x[k] * ab_y + g_0_y_yyyy_xy[k];

                g_0_y_yyyyy_y[k] = g_yyyy_y[k] - g_0_y_yyyy_y[k] * ab_y + g_0_y_yyyy_yy[k];

                g_0_y_yyyyy_z[k] = g_yyyy_z[k] - g_0_y_yyyy_z[k] * ab_y + g_0_y_yyyy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyz_x = cbuffer.data(hp_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_yyyyz_y = cbuffer.data(hp_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_yyyyz_z = cbuffer.data(hp_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_x, g_0_y_yyyy_xz, g_0_y_yyyy_y, g_0_y_yyyy_yz, g_0_y_yyyy_z, g_0_y_yyyy_zz, g_0_y_yyyyz_x, g_0_y_yyyyz_y, g_0_y_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyz_x[k] = -g_0_y_yyyy_x[k] * ab_z + g_0_y_yyyy_xz[k];

                g_0_y_yyyyz_y[k] = -g_0_y_yyyy_y[k] * ab_z + g_0_y_yyyy_yz[k];

                g_0_y_yyyyz_z[k] = -g_0_y_yyyy_z[k] * ab_z + g_0_y_yyyy_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzz_x = cbuffer.data(hp_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_yyyzz_y = cbuffer.data(hp_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_yyyzz_z = cbuffer.data(hp_geom_01_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyz_x, g_0_y_yyyz_xz, g_0_y_yyyz_y, g_0_y_yyyz_yz, g_0_y_yyyz_z, g_0_y_yyyz_zz, g_0_y_yyyzz_x, g_0_y_yyyzz_y, g_0_y_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzz_x[k] = -g_0_y_yyyz_x[k] * ab_z + g_0_y_yyyz_xz[k];

                g_0_y_yyyzz_y[k] = -g_0_y_yyyz_y[k] * ab_z + g_0_y_yyyz_yz[k];

                g_0_y_yyyzz_z[k] = -g_0_y_yyyz_z[k] * ab_z + g_0_y_yyyz_zz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzz_x = cbuffer.data(hp_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_yyzzz_y = cbuffer.data(hp_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_yyzzz_z = cbuffer.data(hp_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzz_x, g_0_y_yyzz_xz, g_0_y_yyzz_y, g_0_y_yyzz_yz, g_0_y_yyzz_z, g_0_y_yyzz_zz, g_0_y_yyzzz_x, g_0_y_yyzzz_y, g_0_y_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzz_x[k] = -g_0_y_yyzz_x[k] * ab_z + g_0_y_yyzz_xz[k];

                g_0_y_yyzzz_y[k] = -g_0_y_yyzz_y[k] * ab_z + g_0_y_yyzz_yz[k];

                g_0_y_yyzzz_z[k] = -g_0_y_yyzz_z[k] * ab_z + g_0_y_yyzz_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzz_x = cbuffer.data(hp_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_yzzzz_y = cbuffer.data(hp_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_yzzzz_z = cbuffer.data(hp_geom_01_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzz_x, g_0_y_yzzz_xz, g_0_y_yzzz_y, g_0_y_yzzz_yz, g_0_y_yzzz_z, g_0_y_yzzz_zz, g_0_y_yzzzz_x, g_0_y_yzzzz_y, g_0_y_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzz_x[k] = -g_0_y_yzzz_x[k] * ab_z + g_0_y_yzzz_xz[k];

                g_0_y_yzzzz_y[k] = -g_0_y_yzzz_y[k] * ab_z + g_0_y_yzzz_yz[k];

                g_0_y_yzzzz_z[k] = -g_0_y_yzzz_z[k] * ab_z + g_0_y_yzzz_zz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzz_x = cbuffer.data(hp_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_zzzzz_y = cbuffer.data(hp_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_zzzzz_z = cbuffer.data(hp_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_x, g_0_y_zzzz_xz, g_0_y_zzzz_y, g_0_y_zzzz_yz, g_0_y_zzzz_z, g_0_y_zzzz_zz, g_0_y_zzzzz_x, g_0_y_zzzzz_y, g_0_y_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzz_x[k] = -g_0_y_zzzz_x[k] * ab_z + g_0_y_zzzz_xz[k];

                g_0_y_zzzzz_y[k] = -g_0_y_zzzz_y[k] * ab_z + g_0_y_zzzz_yz[k];

                g_0_y_zzzzz_z[k] = -g_0_y_zzzz_z[k] * ab_z + g_0_y_zzzz_zz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxx_x = cbuffer.data(hp_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_xxxxx_y = cbuffer.data(hp_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_xxxxx_z = cbuffer.data(hp_geom_01_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_x, g_0_z_xxxx_xx, g_0_z_xxxx_xy, g_0_z_xxxx_xz, g_0_z_xxxx_y, g_0_z_xxxx_z, g_0_z_xxxxx_x, g_0_z_xxxxx_y, g_0_z_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxx_x[k] = -g_0_z_xxxx_x[k] * ab_x + g_0_z_xxxx_xx[k];

                g_0_z_xxxxx_y[k] = -g_0_z_xxxx_y[k] * ab_x + g_0_z_xxxx_xy[k];

                g_0_z_xxxxx_z[k] = -g_0_z_xxxx_z[k] * ab_x + g_0_z_xxxx_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxy_x = cbuffer.data(hp_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_xxxxy_y = cbuffer.data(hp_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_xxxxy_z = cbuffer.data(hp_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxy_x, g_0_z_xxxxy_y, g_0_z_xxxxy_z, g_0_z_xxxy_x, g_0_z_xxxy_xx, g_0_z_xxxy_xy, g_0_z_xxxy_xz, g_0_z_xxxy_y, g_0_z_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxy_x[k] = -g_0_z_xxxy_x[k] * ab_x + g_0_z_xxxy_xx[k];

                g_0_z_xxxxy_y[k] = -g_0_z_xxxy_y[k] * ab_x + g_0_z_xxxy_xy[k];

                g_0_z_xxxxy_z[k] = -g_0_z_xxxy_z[k] * ab_x + g_0_z_xxxy_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxz_x = cbuffer.data(hp_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_xxxxz_y = cbuffer.data(hp_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_xxxxz_z = cbuffer.data(hp_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxz_x, g_0_z_xxxxz_y, g_0_z_xxxxz_z, g_0_z_xxxz_x, g_0_z_xxxz_xx, g_0_z_xxxz_xy, g_0_z_xxxz_xz, g_0_z_xxxz_y, g_0_z_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxz_x[k] = -g_0_z_xxxz_x[k] * ab_x + g_0_z_xxxz_xx[k];

                g_0_z_xxxxz_y[k] = -g_0_z_xxxz_y[k] * ab_x + g_0_z_xxxz_xy[k];

                g_0_z_xxxxz_z[k] = -g_0_z_xxxz_z[k] * ab_x + g_0_z_xxxz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyy_x = cbuffer.data(hp_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_xxxyy_y = cbuffer.data(hp_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_xxxyy_z = cbuffer.data(hp_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyy_x, g_0_z_xxxyy_y, g_0_z_xxxyy_z, g_0_z_xxyy_x, g_0_z_xxyy_xx, g_0_z_xxyy_xy, g_0_z_xxyy_xz, g_0_z_xxyy_y, g_0_z_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyy_x[k] = -g_0_z_xxyy_x[k] * ab_x + g_0_z_xxyy_xx[k];

                g_0_z_xxxyy_y[k] = -g_0_z_xxyy_y[k] * ab_x + g_0_z_xxyy_xy[k];

                g_0_z_xxxyy_z[k] = -g_0_z_xxyy_z[k] * ab_x + g_0_z_xxyy_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyz_x = cbuffer.data(hp_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_xxxyz_y = cbuffer.data(hp_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_xxxyz_z = cbuffer.data(hp_geom_01_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyz_x, g_0_z_xxxyz_y, g_0_z_xxxyz_z, g_0_z_xxyz_x, g_0_z_xxyz_xx, g_0_z_xxyz_xy, g_0_z_xxyz_xz, g_0_z_xxyz_y, g_0_z_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyz_x[k] = -g_0_z_xxyz_x[k] * ab_x + g_0_z_xxyz_xx[k];

                g_0_z_xxxyz_y[k] = -g_0_z_xxyz_y[k] * ab_x + g_0_z_xxyz_xy[k];

                g_0_z_xxxyz_z[k] = -g_0_z_xxyz_z[k] * ab_x + g_0_z_xxyz_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzz_x = cbuffer.data(hp_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_xxxzz_y = cbuffer.data(hp_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_xxxzz_z = cbuffer.data(hp_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzz_x, g_0_z_xxxzz_y, g_0_z_xxxzz_z, g_0_z_xxzz_x, g_0_z_xxzz_xx, g_0_z_xxzz_xy, g_0_z_xxzz_xz, g_0_z_xxzz_y, g_0_z_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzz_x[k] = -g_0_z_xxzz_x[k] * ab_x + g_0_z_xxzz_xx[k];

                g_0_z_xxxzz_y[k] = -g_0_z_xxzz_y[k] * ab_x + g_0_z_xxzz_xy[k];

                g_0_z_xxxzz_z[k] = -g_0_z_xxzz_z[k] * ab_x + g_0_z_xxzz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyy_x = cbuffer.data(hp_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_xxyyy_y = cbuffer.data(hp_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_xxyyy_z = cbuffer.data(hp_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyy_x, g_0_z_xxyyy_y, g_0_z_xxyyy_z, g_0_z_xyyy_x, g_0_z_xyyy_xx, g_0_z_xyyy_xy, g_0_z_xyyy_xz, g_0_z_xyyy_y, g_0_z_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyy_x[k] = -g_0_z_xyyy_x[k] * ab_x + g_0_z_xyyy_xx[k];

                g_0_z_xxyyy_y[k] = -g_0_z_xyyy_y[k] * ab_x + g_0_z_xyyy_xy[k];

                g_0_z_xxyyy_z[k] = -g_0_z_xyyy_z[k] * ab_x + g_0_z_xyyy_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyz_x = cbuffer.data(hp_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_xxyyz_y = cbuffer.data(hp_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_xxyyz_z = cbuffer.data(hp_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyz_x, g_0_z_xxyyz_y, g_0_z_xxyyz_z, g_0_z_xyyz_x, g_0_z_xyyz_xx, g_0_z_xyyz_xy, g_0_z_xyyz_xz, g_0_z_xyyz_y, g_0_z_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyz_x[k] = -g_0_z_xyyz_x[k] * ab_x + g_0_z_xyyz_xx[k];

                g_0_z_xxyyz_y[k] = -g_0_z_xyyz_y[k] * ab_x + g_0_z_xyyz_xy[k];

                g_0_z_xxyyz_z[k] = -g_0_z_xyyz_z[k] * ab_x + g_0_z_xyyz_xz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzz_x = cbuffer.data(hp_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_xxyzz_y = cbuffer.data(hp_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_xxyzz_z = cbuffer.data(hp_geom_01_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzz_x, g_0_z_xxyzz_y, g_0_z_xxyzz_z, g_0_z_xyzz_x, g_0_z_xyzz_xx, g_0_z_xyzz_xy, g_0_z_xyzz_xz, g_0_z_xyzz_y, g_0_z_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzz_x[k] = -g_0_z_xyzz_x[k] * ab_x + g_0_z_xyzz_xx[k];

                g_0_z_xxyzz_y[k] = -g_0_z_xyzz_y[k] * ab_x + g_0_z_xyzz_xy[k];

                g_0_z_xxyzz_z[k] = -g_0_z_xyzz_z[k] * ab_x + g_0_z_xyzz_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzz_x = cbuffer.data(hp_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_xxzzz_y = cbuffer.data(hp_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_xxzzz_z = cbuffer.data(hp_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzz_x, g_0_z_xxzzz_y, g_0_z_xxzzz_z, g_0_z_xzzz_x, g_0_z_xzzz_xx, g_0_z_xzzz_xy, g_0_z_xzzz_xz, g_0_z_xzzz_y, g_0_z_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzz_x[k] = -g_0_z_xzzz_x[k] * ab_x + g_0_z_xzzz_xx[k];

                g_0_z_xxzzz_y[k] = -g_0_z_xzzz_y[k] * ab_x + g_0_z_xzzz_xy[k];

                g_0_z_xxzzz_z[k] = -g_0_z_xzzz_z[k] * ab_x + g_0_z_xzzz_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyy_x = cbuffer.data(hp_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_xyyyy_y = cbuffer.data(hp_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_xyyyy_z = cbuffer.data(hp_geom_01_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyy_x, g_0_z_xyyyy_y, g_0_z_xyyyy_z, g_0_z_yyyy_x, g_0_z_yyyy_xx, g_0_z_yyyy_xy, g_0_z_yyyy_xz, g_0_z_yyyy_y, g_0_z_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyy_x[k] = -g_0_z_yyyy_x[k] * ab_x + g_0_z_yyyy_xx[k];

                g_0_z_xyyyy_y[k] = -g_0_z_yyyy_y[k] * ab_x + g_0_z_yyyy_xy[k];

                g_0_z_xyyyy_z[k] = -g_0_z_yyyy_z[k] * ab_x + g_0_z_yyyy_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyz_x = cbuffer.data(hp_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_xyyyz_y = cbuffer.data(hp_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_xyyyz_z = cbuffer.data(hp_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyz_x, g_0_z_xyyyz_y, g_0_z_xyyyz_z, g_0_z_yyyz_x, g_0_z_yyyz_xx, g_0_z_yyyz_xy, g_0_z_yyyz_xz, g_0_z_yyyz_y, g_0_z_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyz_x[k] = -g_0_z_yyyz_x[k] * ab_x + g_0_z_yyyz_xx[k];

                g_0_z_xyyyz_y[k] = -g_0_z_yyyz_y[k] * ab_x + g_0_z_yyyz_xy[k];

                g_0_z_xyyyz_z[k] = -g_0_z_yyyz_z[k] * ab_x + g_0_z_yyyz_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzz_x = cbuffer.data(hp_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_xyyzz_y = cbuffer.data(hp_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_xyyzz_z = cbuffer.data(hp_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzz_x, g_0_z_xyyzz_y, g_0_z_xyyzz_z, g_0_z_yyzz_x, g_0_z_yyzz_xx, g_0_z_yyzz_xy, g_0_z_yyzz_xz, g_0_z_yyzz_y, g_0_z_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzz_x[k] = -g_0_z_yyzz_x[k] * ab_x + g_0_z_yyzz_xx[k];

                g_0_z_xyyzz_y[k] = -g_0_z_yyzz_y[k] * ab_x + g_0_z_yyzz_xy[k];

                g_0_z_xyyzz_z[k] = -g_0_z_yyzz_z[k] * ab_x + g_0_z_yyzz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzz_x = cbuffer.data(hp_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_xyzzz_y = cbuffer.data(hp_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_xyzzz_z = cbuffer.data(hp_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzz_x, g_0_z_xyzzz_y, g_0_z_xyzzz_z, g_0_z_yzzz_x, g_0_z_yzzz_xx, g_0_z_yzzz_xy, g_0_z_yzzz_xz, g_0_z_yzzz_y, g_0_z_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzz_x[k] = -g_0_z_yzzz_x[k] * ab_x + g_0_z_yzzz_xx[k];

                g_0_z_xyzzz_y[k] = -g_0_z_yzzz_y[k] * ab_x + g_0_z_yzzz_xy[k];

                g_0_z_xyzzz_z[k] = -g_0_z_yzzz_z[k] * ab_x + g_0_z_yzzz_xz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzz_x = cbuffer.data(hp_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_xzzzz_y = cbuffer.data(hp_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_xzzzz_z = cbuffer.data(hp_geom_01_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzz_x, g_0_z_xzzzz_y, g_0_z_xzzzz_z, g_0_z_zzzz_x, g_0_z_zzzz_xx, g_0_z_zzzz_xy, g_0_z_zzzz_xz, g_0_z_zzzz_y, g_0_z_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzz_x[k] = -g_0_z_zzzz_x[k] * ab_x + g_0_z_zzzz_xx[k];

                g_0_z_xzzzz_y[k] = -g_0_z_zzzz_y[k] * ab_x + g_0_z_zzzz_xy[k];

                g_0_z_xzzzz_z[k] = -g_0_z_zzzz_z[k] * ab_x + g_0_z_zzzz_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyy_x = cbuffer.data(hp_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_yyyyy_y = cbuffer.data(hp_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_yyyyy_z = cbuffer.data(hp_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_x, g_0_z_yyyy_xy, g_0_z_yyyy_y, g_0_z_yyyy_yy, g_0_z_yyyy_yz, g_0_z_yyyy_z, g_0_z_yyyyy_x, g_0_z_yyyyy_y, g_0_z_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyy_x[k] = -g_0_z_yyyy_x[k] * ab_y + g_0_z_yyyy_xy[k];

                g_0_z_yyyyy_y[k] = -g_0_z_yyyy_y[k] * ab_y + g_0_z_yyyy_yy[k];

                g_0_z_yyyyy_z[k] = -g_0_z_yyyy_z[k] * ab_y + g_0_z_yyyy_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyz_x = cbuffer.data(hp_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_yyyyz_y = cbuffer.data(hp_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_yyyyz_z = cbuffer.data(hp_geom_01_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyz_x, g_0_z_yyyyz_y, g_0_z_yyyyz_z, g_0_z_yyyz_x, g_0_z_yyyz_xy, g_0_z_yyyz_y, g_0_z_yyyz_yy, g_0_z_yyyz_yz, g_0_z_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyz_x[k] = -g_0_z_yyyz_x[k] * ab_y + g_0_z_yyyz_xy[k];

                g_0_z_yyyyz_y[k] = -g_0_z_yyyz_y[k] * ab_y + g_0_z_yyyz_yy[k];

                g_0_z_yyyyz_z[k] = -g_0_z_yyyz_z[k] * ab_y + g_0_z_yyyz_yz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzz_x = cbuffer.data(hp_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_yyyzz_y = cbuffer.data(hp_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_yyyzz_z = cbuffer.data(hp_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzz_x, g_0_z_yyyzz_y, g_0_z_yyyzz_z, g_0_z_yyzz_x, g_0_z_yyzz_xy, g_0_z_yyzz_y, g_0_z_yyzz_yy, g_0_z_yyzz_yz, g_0_z_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzz_x[k] = -g_0_z_yyzz_x[k] * ab_y + g_0_z_yyzz_xy[k];

                g_0_z_yyyzz_y[k] = -g_0_z_yyzz_y[k] * ab_y + g_0_z_yyzz_yy[k];

                g_0_z_yyyzz_z[k] = -g_0_z_yyzz_z[k] * ab_y + g_0_z_yyzz_yz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzz_x = cbuffer.data(hp_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_yyzzz_y = cbuffer.data(hp_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_yyzzz_z = cbuffer.data(hp_geom_01_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzz_x, g_0_z_yyzzz_y, g_0_z_yyzzz_z, g_0_z_yzzz_x, g_0_z_yzzz_xy, g_0_z_yzzz_y, g_0_z_yzzz_yy, g_0_z_yzzz_yz, g_0_z_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzz_x[k] = -g_0_z_yzzz_x[k] * ab_y + g_0_z_yzzz_xy[k];

                g_0_z_yyzzz_y[k] = -g_0_z_yzzz_y[k] * ab_y + g_0_z_yzzz_yy[k];

                g_0_z_yyzzz_z[k] = -g_0_z_yzzz_z[k] * ab_y + g_0_z_yzzz_yz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzz_x = cbuffer.data(hp_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_yzzzz_y = cbuffer.data(hp_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_yzzzz_z = cbuffer.data(hp_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzz_x, g_0_z_yzzzz_y, g_0_z_yzzzz_z, g_0_z_zzzz_x, g_0_z_zzzz_xy, g_0_z_zzzz_y, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzz_x[k] = -g_0_z_zzzz_x[k] * ab_y + g_0_z_zzzz_xy[k];

                g_0_z_yzzzz_y[k] = -g_0_z_zzzz_y[k] * ab_y + g_0_z_zzzz_yy[k];

                g_0_z_yzzzz_z[k] = -g_0_z_zzzz_z[k] * ab_y + g_0_z_zzzz_yz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzz_x = cbuffer.data(hp_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_zzzzz_y = cbuffer.data(hp_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_zzzzz_z = cbuffer.data(hp_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_x, g_0_z_zzzz_xz, g_0_z_zzzz_y, g_0_z_zzzz_yz, g_0_z_zzzz_z, g_0_z_zzzz_zz, g_0_z_zzzzz_x, g_0_z_zzzzz_y, g_0_z_zzzzz_z, g_zzzz_x, g_zzzz_y, g_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzz_x[k] = g_zzzz_x[k] - g_0_z_zzzz_x[k] * ab_z + g_0_z_zzzz_xz[k];

                g_0_z_zzzzz_y[k] = g_zzzz_y[k] - g_0_z_zzzz_y[k] * ab_z + g_0_z_zzzz_yz[k];

                g_0_z_zzzzz_z[k] = g_zzzz_z[k] - g_0_z_zzzz_z[k] * ab_z + g_0_z_zzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

