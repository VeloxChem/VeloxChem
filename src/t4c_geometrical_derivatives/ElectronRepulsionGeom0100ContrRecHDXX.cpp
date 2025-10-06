#include "ElectronRepulsionGeom0100ContrRecHDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_hdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hdxx,
                                            const size_t idx_gdxx,
                                            const size_t idx_geom_01_gdxx,
                                            const size_t idx_geom_01_gfxx,
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
            /// Set up components of auxilary buffer : GDSS

            const auto gd_off = idx_gdxx + i * dcomps + j;

            auto g_xxxx_xx = cbuffer.data(gd_off + 0 * ccomps * dcomps);

            auto g_xxxx_xy = cbuffer.data(gd_off + 1 * ccomps * dcomps);

            auto g_xxxx_xz = cbuffer.data(gd_off + 2 * ccomps * dcomps);

            auto g_xxxx_yy = cbuffer.data(gd_off + 3 * ccomps * dcomps);

            auto g_xxxx_yz = cbuffer.data(gd_off + 4 * ccomps * dcomps);

            auto g_xxxx_zz = cbuffer.data(gd_off + 5 * ccomps * dcomps);

            auto g_xxxy_xx = cbuffer.data(gd_off + 6 * ccomps * dcomps);

            auto g_xxxy_xy = cbuffer.data(gd_off + 7 * ccomps * dcomps);

            auto g_xxxy_xz = cbuffer.data(gd_off + 8 * ccomps * dcomps);

            auto g_xxxy_yy = cbuffer.data(gd_off + 9 * ccomps * dcomps);

            auto g_xxxy_yz = cbuffer.data(gd_off + 10 * ccomps * dcomps);

            auto g_xxxy_zz = cbuffer.data(gd_off + 11 * ccomps * dcomps);

            auto g_xxxz_xx = cbuffer.data(gd_off + 12 * ccomps * dcomps);

            auto g_xxxz_xy = cbuffer.data(gd_off + 13 * ccomps * dcomps);

            auto g_xxxz_xz = cbuffer.data(gd_off + 14 * ccomps * dcomps);

            auto g_xxxz_yy = cbuffer.data(gd_off + 15 * ccomps * dcomps);

            auto g_xxxz_yz = cbuffer.data(gd_off + 16 * ccomps * dcomps);

            auto g_xxxz_zz = cbuffer.data(gd_off + 17 * ccomps * dcomps);

            auto g_xxyy_xx = cbuffer.data(gd_off + 18 * ccomps * dcomps);

            auto g_xxyy_xy = cbuffer.data(gd_off + 19 * ccomps * dcomps);

            auto g_xxyy_xz = cbuffer.data(gd_off + 20 * ccomps * dcomps);

            auto g_xxyy_yy = cbuffer.data(gd_off + 21 * ccomps * dcomps);

            auto g_xxyy_yz = cbuffer.data(gd_off + 22 * ccomps * dcomps);

            auto g_xxyy_zz = cbuffer.data(gd_off + 23 * ccomps * dcomps);

            auto g_xxyz_xx = cbuffer.data(gd_off + 24 * ccomps * dcomps);

            auto g_xxyz_xy = cbuffer.data(gd_off + 25 * ccomps * dcomps);

            auto g_xxyz_xz = cbuffer.data(gd_off + 26 * ccomps * dcomps);

            auto g_xxyz_yy = cbuffer.data(gd_off + 27 * ccomps * dcomps);

            auto g_xxyz_yz = cbuffer.data(gd_off + 28 * ccomps * dcomps);

            auto g_xxyz_zz = cbuffer.data(gd_off + 29 * ccomps * dcomps);

            auto g_xxzz_xx = cbuffer.data(gd_off + 30 * ccomps * dcomps);

            auto g_xxzz_xy = cbuffer.data(gd_off + 31 * ccomps * dcomps);

            auto g_xxzz_xz = cbuffer.data(gd_off + 32 * ccomps * dcomps);

            auto g_xxzz_yy = cbuffer.data(gd_off + 33 * ccomps * dcomps);

            auto g_xxzz_yz = cbuffer.data(gd_off + 34 * ccomps * dcomps);

            auto g_xxzz_zz = cbuffer.data(gd_off + 35 * ccomps * dcomps);

            auto g_xyyy_xx = cbuffer.data(gd_off + 36 * ccomps * dcomps);

            auto g_xyyy_xy = cbuffer.data(gd_off + 37 * ccomps * dcomps);

            auto g_xyyy_xz = cbuffer.data(gd_off + 38 * ccomps * dcomps);

            auto g_xyyy_yy = cbuffer.data(gd_off + 39 * ccomps * dcomps);

            auto g_xyyy_yz = cbuffer.data(gd_off + 40 * ccomps * dcomps);

            auto g_xyyy_zz = cbuffer.data(gd_off + 41 * ccomps * dcomps);

            auto g_xyyz_xx = cbuffer.data(gd_off + 42 * ccomps * dcomps);

            auto g_xyyz_xy = cbuffer.data(gd_off + 43 * ccomps * dcomps);

            auto g_xyyz_xz = cbuffer.data(gd_off + 44 * ccomps * dcomps);

            auto g_xyyz_yy = cbuffer.data(gd_off + 45 * ccomps * dcomps);

            auto g_xyyz_yz = cbuffer.data(gd_off + 46 * ccomps * dcomps);

            auto g_xyyz_zz = cbuffer.data(gd_off + 47 * ccomps * dcomps);

            auto g_xyzz_xx = cbuffer.data(gd_off + 48 * ccomps * dcomps);

            auto g_xyzz_xy = cbuffer.data(gd_off + 49 * ccomps * dcomps);

            auto g_xyzz_xz = cbuffer.data(gd_off + 50 * ccomps * dcomps);

            auto g_xyzz_yy = cbuffer.data(gd_off + 51 * ccomps * dcomps);

            auto g_xyzz_yz = cbuffer.data(gd_off + 52 * ccomps * dcomps);

            auto g_xyzz_zz = cbuffer.data(gd_off + 53 * ccomps * dcomps);

            auto g_xzzz_xx = cbuffer.data(gd_off + 54 * ccomps * dcomps);

            auto g_xzzz_xy = cbuffer.data(gd_off + 55 * ccomps * dcomps);

            auto g_xzzz_xz = cbuffer.data(gd_off + 56 * ccomps * dcomps);

            auto g_xzzz_yy = cbuffer.data(gd_off + 57 * ccomps * dcomps);

            auto g_xzzz_yz = cbuffer.data(gd_off + 58 * ccomps * dcomps);

            auto g_xzzz_zz = cbuffer.data(gd_off + 59 * ccomps * dcomps);

            auto g_yyyy_xx = cbuffer.data(gd_off + 60 * ccomps * dcomps);

            auto g_yyyy_xy = cbuffer.data(gd_off + 61 * ccomps * dcomps);

            auto g_yyyy_xz = cbuffer.data(gd_off + 62 * ccomps * dcomps);

            auto g_yyyy_yy = cbuffer.data(gd_off + 63 * ccomps * dcomps);

            auto g_yyyy_yz = cbuffer.data(gd_off + 64 * ccomps * dcomps);

            auto g_yyyy_zz = cbuffer.data(gd_off + 65 * ccomps * dcomps);

            auto g_yyyz_xx = cbuffer.data(gd_off + 66 * ccomps * dcomps);

            auto g_yyyz_xy = cbuffer.data(gd_off + 67 * ccomps * dcomps);

            auto g_yyyz_xz = cbuffer.data(gd_off + 68 * ccomps * dcomps);

            auto g_yyyz_yy = cbuffer.data(gd_off + 69 * ccomps * dcomps);

            auto g_yyyz_yz = cbuffer.data(gd_off + 70 * ccomps * dcomps);

            auto g_yyyz_zz = cbuffer.data(gd_off + 71 * ccomps * dcomps);

            auto g_yyzz_xx = cbuffer.data(gd_off + 72 * ccomps * dcomps);

            auto g_yyzz_xy = cbuffer.data(gd_off + 73 * ccomps * dcomps);

            auto g_yyzz_xz = cbuffer.data(gd_off + 74 * ccomps * dcomps);

            auto g_yyzz_yy = cbuffer.data(gd_off + 75 * ccomps * dcomps);

            auto g_yyzz_yz = cbuffer.data(gd_off + 76 * ccomps * dcomps);

            auto g_yyzz_zz = cbuffer.data(gd_off + 77 * ccomps * dcomps);

            auto g_yzzz_xx = cbuffer.data(gd_off + 78 * ccomps * dcomps);

            auto g_yzzz_xy = cbuffer.data(gd_off + 79 * ccomps * dcomps);

            auto g_yzzz_xz = cbuffer.data(gd_off + 80 * ccomps * dcomps);

            auto g_yzzz_yy = cbuffer.data(gd_off + 81 * ccomps * dcomps);

            auto g_yzzz_yz = cbuffer.data(gd_off + 82 * ccomps * dcomps);

            auto g_yzzz_zz = cbuffer.data(gd_off + 83 * ccomps * dcomps);

            auto g_zzzz_xx = cbuffer.data(gd_off + 84 * ccomps * dcomps);

            auto g_zzzz_xy = cbuffer.data(gd_off + 85 * ccomps * dcomps);

            auto g_zzzz_xz = cbuffer.data(gd_off + 86 * ccomps * dcomps);

            auto g_zzzz_yy = cbuffer.data(gd_off + 87 * ccomps * dcomps);

            auto g_zzzz_yz = cbuffer.data(gd_off + 88 * ccomps * dcomps);

            auto g_zzzz_zz = cbuffer.data(gd_off + 89 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : GFSS

            const auto gf_geom_01_off = idx_geom_01_gfxx + i * dcomps + j;

            auto g_0_x_xxxx_xxx = cbuffer.data(gf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxy = cbuffer.data(gf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxz = cbuffer.data(gf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xyy = cbuffer.data(gf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xyz = cbuffer.data(gf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xzz = cbuffer.data(gf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_yyy = cbuffer.data(gf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_yyz = cbuffer.data(gf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_yzz = cbuffer.data(gf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_zzz = cbuffer.data(gf_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxy_xxx = cbuffer.data(gf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxy_xxy = cbuffer.data(gf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxy_xxz = cbuffer.data(gf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxy_xyy = cbuffer.data(gf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxy_xyz = cbuffer.data(gf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxy_xzz = cbuffer.data(gf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxy_yyy = cbuffer.data(gf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxy_yyz = cbuffer.data(gf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxy_yzz = cbuffer.data(gf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxy_zzz = cbuffer.data(gf_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxz_xxx = cbuffer.data(gf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxz_xxy = cbuffer.data(gf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxz_xxz = cbuffer.data(gf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxz_xyy = cbuffer.data(gf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxz_xyz = cbuffer.data(gf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxz_xzz = cbuffer.data(gf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxz_yyy = cbuffer.data(gf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxz_yyz = cbuffer.data(gf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxz_yzz = cbuffer.data(gf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxz_zzz = cbuffer.data(gf_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxyy_xxx = cbuffer.data(gf_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxyy_xxy = cbuffer.data(gf_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxyy_xxz = cbuffer.data(gf_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxyy_xyy = cbuffer.data(gf_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxyy_xyz = cbuffer.data(gf_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxyy_xzz = cbuffer.data(gf_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxyy_yyy = cbuffer.data(gf_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyy_yyz = cbuffer.data(gf_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyy_yzz = cbuffer.data(gf_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyy_zzz = cbuffer.data(gf_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyz_xxx = cbuffer.data(gf_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyz_xxy = cbuffer.data(gf_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxyz_xxz = cbuffer.data(gf_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxyz_xyy = cbuffer.data(gf_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxyz_xyz = cbuffer.data(gf_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxyz_xzz = cbuffer.data(gf_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyz_yyy = cbuffer.data(gf_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyz_yyz = cbuffer.data(gf_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxyz_yzz = cbuffer.data(gf_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyz_zzz = cbuffer.data(gf_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxzz_xxx = cbuffer.data(gf_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxzz_xxy = cbuffer.data(gf_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxzz_xxz = cbuffer.data(gf_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxzz_xyy = cbuffer.data(gf_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxzz_xyz = cbuffer.data(gf_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxzz_xzz = cbuffer.data(gf_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxzz_yyy = cbuffer.data(gf_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxzz_yyz = cbuffer.data(gf_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxzz_yzz = cbuffer.data(gf_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxzz_zzz = cbuffer.data(gf_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xyyy_xxx = cbuffer.data(gf_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyyy_xxy = cbuffer.data(gf_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyyy_xxz = cbuffer.data(gf_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyyy_xyy = cbuffer.data(gf_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyy_xyz = cbuffer.data(gf_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyy_xzz = cbuffer.data(gf_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyyy_yyy = cbuffer.data(gf_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyy_yyz = cbuffer.data(gf_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyy_yzz = cbuffer.data(gf_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyyy_zzz = cbuffer.data(gf_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyyz_xxx = cbuffer.data(gf_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyz_xxy = cbuffer.data(gf_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyyz_xxz = cbuffer.data(gf_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyz_xyy = cbuffer.data(gf_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyz_xyz = cbuffer.data(gf_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyyz_xzz = cbuffer.data(gf_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyz_yyy = cbuffer.data(gf_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyz_yyz = cbuffer.data(gf_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyyz_yzz = cbuffer.data(gf_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyyz_zzz = cbuffer.data(gf_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyzz_xxx = cbuffer.data(gf_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyzz_xxy = cbuffer.data(gf_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyzz_xxz = cbuffer.data(gf_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyzz_xyy = cbuffer.data(gf_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xyzz_xyz = cbuffer.data(gf_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xyzz_xzz = cbuffer.data(gf_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xyzz_yyy = cbuffer.data(gf_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xyzz_yyz = cbuffer.data(gf_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xyzz_yzz = cbuffer.data(gf_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xyzz_zzz = cbuffer.data(gf_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xzzz_xxx = cbuffer.data(gf_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xzzz_xxy = cbuffer.data(gf_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xzzz_xxz = cbuffer.data(gf_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xzzz_xyy = cbuffer.data(gf_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xzzz_xyz = cbuffer.data(gf_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xzzz_xzz = cbuffer.data(gf_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xzzz_yyy = cbuffer.data(gf_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xzzz_yyz = cbuffer.data(gf_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xzzz_yzz = cbuffer.data(gf_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xzzz_zzz = cbuffer.data(gf_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyyy_xxx = cbuffer.data(gf_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyyy_xxy = cbuffer.data(gf_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyyy_xxz = cbuffer.data(gf_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyyy_xyy = cbuffer.data(gf_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyyy_xyz = cbuffer.data(gf_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyyy_xzz = cbuffer.data(gf_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyyy_yyy = cbuffer.data(gf_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyyy_yyz = cbuffer.data(gf_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyyy_yzz = cbuffer.data(gf_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyyy_zzz = cbuffer.data(gf_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyyz_xxx = cbuffer.data(gf_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyyz_xxy = cbuffer.data(gf_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyyz_xxz = cbuffer.data(gf_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyyz_xyy = cbuffer.data(gf_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yyyz_xyz = cbuffer.data(gf_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yyyz_xzz = cbuffer.data(gf_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yyyz_yyy = cbuffer.data(gf_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yyyz_yyz = cbuffer.data(gf_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yyyz_yzz = cbuffer.data(gf_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yyyz_zzz = cbuffer.data(gf_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_yyzz_xxx = cbuffer.data(gf_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yyzz_xxy = cbuffer.data(gf_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yyzz_xxz = cbuffer.data(gf_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yyzz_xyy = cbuffer.data(gf_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yyzz_xyz = cbuffer.data(gf_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yyzz_xzz = cbuffer.data(gf_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyzz_yyy = cbuffer.data(gf_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyzz_yyz = cbuffer.data(gf_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyzz_yzz = cbuffer.data(gf_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyzz_zzz = cbuffer.data(gf_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yzzz_xxx = cbuffer.data(gf_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yzzz_xxy = cbuffer.data(gf_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yzzz_xxz = cbuffer.data(gf_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yzzz_xyy = cbuffer.data(gf_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yzzz_xyz = cbuffer.data(gf_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yzzz_xzz = cbuffer.data(gf_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yzzz_yyy = cbuffer.data(gf_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yzzz_yyz = cbuffer.data(gf_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yzzz_yzz = cbuffer.data(gf_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yzzz_zzz = cbuffer.data(gf_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_zzzz_xxx = cbuffer.data(gf_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_zzzz_xxy = cbuffer.data(gf_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_zzzz_xxz = cbuffer.data(gf_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_zzzz_xyy = cbuffer.data(gf_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_zzzz_xyz = cbuffer.data(gf_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_zzzz_xzz = cbuffer.data(gf_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_zzzz_yyy = cbuffer.data(gf_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_zzzz_yyz = cbuffer.data(gf_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_zzzz_yzz = cbuffer.data(gf_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_zzzz_zzz = cbuffer.data(gf_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xxxx_xxx = cbuffer.data(gf_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxx_xxy = cbuffer.data(gf_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxx_xxz = cbuffer.data(gf_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxxx_xyy = cbuffer.data(gf_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxxx_xyz = cbuffer.data(gf_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxxx_xzz = cbuffer.data(gf_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxxx_yyy = cbuffer.data(gf_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxxx_yyz = cbuffer.data(gf_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxxx_yzz = cbuffer.data(gf_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxxx_zzz = cbuffer.data(gf_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxxy_xxx = cbuffer.data(gf_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxxy_xxy = cbuffer.data(gf_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxxy_xxz = cbuffer.data(gf_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxxy_xyy = cbuffer.data(gf_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxxy_xyz = cbuffer.data(gf_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxxy_xzz = cbuffer.data(gf_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxxy_yyy = cbuffer.data(gf_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxxy_yyz = cbuffer.data(gf_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxxy_yzz = cbuffer.data(gf_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxxy_zzz = cbuffer.data(gf_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxxz_xxx = cbuffer.data(gf_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxxz_xxy = cbuffer.data(gf_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxxz_xxz = cbuffer.data(gf_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxxz_xyy = cbuffer.data(gf_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxxz_xyz = cbuffer.data(gf_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxxz_xzz = cbuffer.data(gf_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxxz_yyy = cbuffer.data(gf_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxxz_yyz = cbuffer.data(gf_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxxz_yzz = cbuffer.data(gf_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxxz_zzz = cbuffer.data(gf_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxyy_xxx = cbuffer.data(gf_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxyy_xxy = cbuffer.data(gf_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxyy_xxz = cbuffer.data(gf_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxyy_xyy = cbuffer.data(gf_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxyy_xyz = cbuffer.data(gf_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxyy_xzz = cbuffer.data(gf_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxyy_yyy = cbuffer.data(gf_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxyy_yyz = cbuffer.data(gf_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxyy_yzz = cbuffer.data(gf_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxyy_zzz = cbuffer.data(gf_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxyz_xxx = cbuffer.data(gf_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxyz_xxy = cbuffer.data(gf_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxyz_xxz = cbuffer.data(gf_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxyz_xyy = cbuffer.data(gf_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxyz_xyz = cbuffer.data(gf_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xxyz_xzz = cbuffer.data(gf_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xxyz_yyy = cbuffer.data(gf_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xxyz_yyz = cbuffer.data(gf_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xxyz_yzz = cbuffer.data(gf_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xxyz_zzz = cbuffer.data(gf_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xxzz_xxx = cbuffer.data(gf_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xxzz_xxy = cbuffer.data(gf_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xxzz_xxz = cbuffer.data(gf_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xxzz_xyy = cbuffer.data(gf_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xxzz_xyz = cbuffer.data(gf_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xxzz_xzz = cbuffer.data(gf_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xxzz_yyy = cbuffer.data(gf_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xxzz_yyz = cbuffer.data(gf_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xxzz_yzz = cbuffer.data(gf_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xxzz_zzz = cbuffer.data(gf_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xyyy_xxx = cbuffer.data(gf_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xyyy_xxy = cbuffer.data(gf_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xyyy_xxz = cbuffer.data(gf_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xyyy_xyy = cbuffer.data(gf_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xyyy_xyz = cbuffer.data(gf_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xyyy_xzz = cbuffer.data(gf_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xyyy_yyy = cbuffer.data(gf_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xyyy_yyz = cbuffer.data(gf_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xyyy_yzz = cbuffer.data(gf_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xyyy_zzz = cbuffer.data(gf_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xyyz_xxx = cbuffer.data(gf_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xyyz_xxy = cbuffer.data(gf_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xyyz_xxz = cbuffer.data(gf_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xyyz_xyy = cbuffer.data(gf_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xyyz_xyz = cbuffer.data(gf_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xyyz_xzz = cbuffer.data(gf_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xyyz_yyy = cbuffer.data(gf_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xyyz_yyz = cbuffer.data(gf_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xyyz_yzz = cbuffer.data(gf_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xyyz_zzz = cbuffer.data(gf_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xyzz_xxx = cbuffer.data(gf_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xyzz_xxy = cbuffer.data(gf_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xyzz_xxz = cbuffer.data(gf_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xyzz_xyy = cbuffer.data(gf_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xyzz_xyz = cbuffer.data(gf_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xyzz_xzz = cbuffer.data(gf_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xyzz_yyy = cbuffer.data(gf_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xyzz_yyz = cbuffer.data(gf_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xyzz_yzz = cbuffer.data(gf_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xyzz_zzz = cbuffer.data(gf_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xzzz_xxx = cbuffer.data(gf_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xzzz_xxy = cbuffer.data(gf_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xzzz_xxz = cbuffer.data(gf_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xzzz_xyy = cbuffer.data(gf_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xzzz_xyz = cbuffer.data(gf_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xzzz_xzz = cbuffer.data(gf_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xzzz_yyy = cbuffer.data(gf_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xzzz_yyz = cbuffer.data(gf_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xzzz_yzz = cbuffer.data(gf_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xzzz_zzz = cbuffer.data(gf_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_yyyy_xxx = cbuffer.data(gf_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_yyyy_xxy = cbuffer.data(gf_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_yyyy_xxz = cbuffer.data(gf_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_yyyy_xyy = cbuffer.data(gf_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_yyyy_xyz = cbuffer.data(gf_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_yyyy_xzz = cbuffer.data(gf_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_yyyy_yyy = cbuffer.data(gf_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_yyyy_yyz = cbuffer.data(gf_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_yyyy_yzz = cbuffer.data(gf_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_yyyy_zzz = cbuffer.data(gf_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_yyyz_xxx = cbuffer.data(gf_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_yyyz_xxy = cbuffer.data(gf_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_yyyz_xxz = cbuffer.data(gf_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_yyyz_xyy = cbuffer.data(gf_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_yyyz_xyz = cbuffer.data(gf_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_yyyz_xzz = cbuffer.data(gf_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_yyyz_yyy = cbuffer.data(gf_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_yyyz_yyz = cbuffer.data(gf_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_yyyz_yzz = cbuffer.data(gf_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_yyyz_zzz = cbuffer.data(gf_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_yyzz_xxx = cbuffer.data(gf_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_yyzz_xxy = cbuffer.data(gf_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_yyzz_xxz = cbuffer.data(gf_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_yyzz_xyy = cbuffer.data(gf_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_yyzz_xyz = cbuffer.data(gf_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_yyzz_xzz = cbuffer.data(gf_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_yyzz_yyy = cbuffer.data(gf_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_yyzz_yyz = cbuffer.data(gf_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_yyzz_yzz = cbuffer.data(gf_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_yyzz_zzz = cbuffer.data(gf_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_yzzz_xxx = cbuffer.data(gf_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_yzzz_xxy = cbuffer.data(gf_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_yzzz_xxz = cbuffer.data(gf_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_yzzz_xyy = cbuffer.data(gf_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_yzzz_xyz = cbuffer.data(gf_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_yzzz_xzz = cbuffer.data(gf_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_yzzz_yyy = cbuffer.data(gf_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_yzzz_yyz = cbuffer.data(gf_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_yzzz_yzz = cbuffer.data(gf_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_yzzz_zzz = cbuffer.data(gf_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_zzzz_xxx = cbuffer.data(gf_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_zzzz_xxy = cbuffer.data(gf_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_zzzz_xxz = cbuffer.data(gf_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_zzzz_xyy = cbuffer.data(gf_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_zzzz_xyz = cbuffer.data(gf_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_zzzz_xzz = cbuffer.data(gf_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_zzzz_yyy = cbuffer.data(gf_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_zzzz_yyz = cbuffer.data(gf_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_zzzz_yzz = cbuffer.data(gf_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_zzzz_zzz = cbuffer.data(gf_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xxxx_xxx = cbuffer.data(gf_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxxx_xxy = cbuffer.data(gf_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxxx_xxz = cbuffer.data(gf_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxxx_xyy = cbuffer.data(gf_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxxx_xyz = cbuffer.data(gf_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxxx_xzz = cbuffer.data(gf_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxxx_yyy = cbuffer.data(gf_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxxx_yyz = cbuffer.data(gf_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxxx_yzz = cbuffer.data(gf_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxxx_zzz = cbuffer.data(gf_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxxy_xxx = cbuffer.data(gf_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxxy_xxy = cbuffer.data(gf_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xxxy_xxz = cbuffer.data(gf_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xxxy_xyy = cbuffer.data(gf_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xxxy_xyz = cbuffer.data(gf_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xxxy_xzz = cbuffer.data(gf_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xxxy_yyy = cbuffer.data(gf_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xxxy_yyz = cbuffer.data(gf_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xxxy_yzz = cbuffer.data(gf_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xxxy_zzz = cbuffer.data(gf_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xxxz_xxx = cbuffer.data(gf_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xxxz_xxy = cbuffer.data(gf_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xxxz_xxz = cbuffer.data(gf_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xxxz_xyy = cbuffer.data(gf_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xxxz_xyz = cbuffer.data(gf_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xxxz_xzz = cbuffer.data(gf_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xxxz_yyy = cbuffer.data(gf_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xxxz_yyz = cbuffer.data(gf_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xxxz_yzz = cbuffer.data(gf_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xxxz_zzz = cbuffer.data(gf_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_xxyy_xxx = cbuffer.data(gf_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xxyy_xxy = cbuffer.data(gf_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xxyy_xxz = cbuffer.data(gf_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xxyy_xyy = cbuffer.data(gf_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xxyy_xyz = cbuffer.data(gf_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xxyy_xzz = cbuffer.data(gf_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxyy_yyy = cbuffer.data(gf_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxyy_yyz = cbuffer.data(gf_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxyy_yzz = cbuffer.data(gf_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxyy_zzz = cbuffer.data(gf_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxyz_xxx = cbuffer.data(gf_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxyz_xxy = cbuffer.data(gf_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxyz_xxz = cbuffer.data(gf_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxyz_xyy = cbuffer.data(gf_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxyz_xyz = cbuffer.data(gf_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xxyz_xzz = cbuffer.data(gf_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xxyz_yyy = cbuffer.data(gf_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xxyz_yyz = cbuffer.data(gf_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xxyz_yzz = cbuffer.data(gf_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xxyz_zzz = cbuffer.data(gf_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xxzz_xxx = cbuffer.data(gf_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xxzz_xxy = cbuffer.data(gf_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xxzz_xxz = cbuffer.data(gf_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xxzz_xyy = cbuffer.data(gf_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xxzz_xyz = cbuffer.data(gf_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xxzz_xzz = cbuffer.data(gf_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xxzz_yyy = cbuffer.data(gf_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xxzz_yyz = cbuffer.data(gf_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xxzz_yzz = cbuffer.data(gf_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xxzz_zzz = cbuffer.data(gf_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_xyyy_xxx = cbuffer.data(gf_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xyyy_xxy = cbuffer.data(gf_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xyyy_xxz = cbuffer.data(gf_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xyyy_xyy = cbuffer.data(gf_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xyyy_xyz = cbuffer.data(gf_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xyyy_xzz = cbuffer.data(gf_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xyyy_yyy = cbuffer.data(gf_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xyyy_yyz = cbuffer.data(gf_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xyyy_yzz = cbuffer.data(gf_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xyyy_zzz = cbuffer.data(gf_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xyyz_xxx = cbuffer.data(gf_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xyyz_xxy = cbuffer.data(gf_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xyyz_xxz = cbuffer.data(gf_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xyyz_xyy = cbuffer.data(gf_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xyyz_xyz = cbuffer.data(gf_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xyyz_xzz = cbuffer.data(gf_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xyyz_yyy = cbuffer.data(gf_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xyyz_yyz = cbuffer.data(gf_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xyyz_yzz = cbuffer.data(gf_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xyyz_zzz = cbuffer.data(gf_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xyzz_xxx = cbuffer.data(gf_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xyzz_xxy = cbuffer.data(gf_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xyzz_xxz = cbuffer.data(gf_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xyzz_xyy = cbuffer.data(gf_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xyzz_xyz = cbuffer.data(gf_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xyzz_xzz = cbuffer.data(gf_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xyzz_yyy = cbuffer.data(gf_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xyzz_yyz = cbuffer.data(gf_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xyzz_yzz = cbuffer.data(gf_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xyzz_zzz = cbuffer.data(gf_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_xzzz_xxx = cbuffer.data(gf_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_xzzz_xxy = cbuffer.data(gf_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_xzzz_xxz = cbuffer.data(gf_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_xzzz_xyy = cbuffer.data(gf_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_xzzz_xyz = cbuffer.data(gf_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_xzzz_xzz = cbuffer.data(gf_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_xzzz_yyy = cbuffer.data(gf_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_xzzz_yyz = cbuffer.data(gf_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_xzzz_yzz = cbuffer.data(gf_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_xzzz_zzz = cbuffer.data(gf_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_yyyy_xxx = cbuffer.data(gf_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_yyyy_xxy = cbuffer.data(gf_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_yyyy_xxz = cbuffer.data(gf_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_yyyy_xyy = cbuffer.data(gf_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_yyyy_xyz = cbuffer.data(gf_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_yyyy_xzz = cbuffer.data(gf_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_yyyy_yyy = cbuffer.data(gf_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_yyyy_yyz = cbuffer.data(gf_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_yyyy_yzz = cbuffer.data(gf_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_yyyy_zzz = cbuffer.data(gf_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_yyyz_xxx = cbuffer.data(gf_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_yyyz_xxy = cbuffer.data(gf_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_yyyz_xxz = cbuffer.data(gf_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_yyyz_xyy = cbuffer.data(gf_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_yyyz_xyz = cbuffer.data(gf_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_yyyz_xzz = cbuffer.data(gf_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_yyyz_yyy = cbuffer.data(gf_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_yyyz_yyz = cbuffer.data(gf_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_yyyz_yzz = cbuffer.data(gf_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_yyyz_zzz = cbuffer.data(gf_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_yyzz_xxx = cbuffer.data(gf_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_yyzz_xxy = cbuffer.data(gf_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_yyzz_xxz = cbuffer.data(gf_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_yyzz_xyy = cbuffer.data(gf_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_yyzz_xyz = cbuffer.data(gf_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_yyzz_xzz = cbuffer.data(gf_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_yyzz_yyy = cbuffer.data(gf_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_yyzz_yyz = cbuffer.data(gf_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_yyzz_yzz = cbuffer.data(gf_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_yyzz_zzz = cbuffer.data(gf_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_yzzz_xxx = cbuffer.data(gf_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_yzzz_xxy = cbuffer.data(gf_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_yzzz_xxz = cbuffer.data(gf_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_yzzz_xyy = cbuffer.data(gf_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_yzzz_xyz = cbuffer.data(gf_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_yzzz_xzz = cbuffer.data(gf_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_yzzz_yyy = cbuffer.data(gf_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_yzzz_yyz = cbuffer.data(gf_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_yzzz_yzz = cbuffer.data(gf_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_yzzz_zzz = cbuffer.data(gf_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_zzzz_xxx = cbuffer.data(gf_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_zzzz_xxy = cbuffer.data(gf_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_zzzz_xxz = cbuffer.data(gf_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_zzzz_xyy = cbuffer.data(gf_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_zzzz_xyz = cbuffer.data(gf_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_zzzz_xzz = cbuffer.data(gf_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_zzzz_yyy = cbuffer.data(gf_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_zzzz_yyz = cbuffer.data(gf_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_zzzz_yzz = cbuffer.data(gf_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_zzzz_zzz = cbuffer.data(gf_geom_01_off + 449 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hdxx

            const auto hd_geom_01_off = idx_geom_01_hdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxx_xx = cbuffer.data(hd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xy = cbuffer.data(hd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xz = cbuffer.data(hd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_yy = cbuffer.data(hd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_yz = cbuffer.data(hd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_zz = cbuffer.data(hd_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xx, g_0_x_xxxx_xxx, g_0_x_xxxx_xxy, g_0_x_xxxx_xxz, g_0_x_xxxx_xy, g_0_x_xxxx_xyy, g_0_x_xxxx_xyz, g_0_x_xxxx_xz, g_0_x_xxxx_xzz, g_0_x_xxxx_yy, g_0_x_xxxx_yz, g_0_x_xxxx_zz, g_0_x_xxxxx_xx, g_0_x_xxxxx_xy, g_0_x_xxxxx_xz, g_0_x_xxxxx_yy, g_0_x_xxxxx_yz, g_0_x_xxxxx_zz, g_xxxx_xx, g_xxxx_xy, g_xxxx_xz, g_xxxx_yy, g_xxxx_yz, g_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxx_xx[k] = g_xxxx_xx[k] - g_0_x_xxxx_xx[k] * ab_x + g_0_x_xxxx_xxx[k];

                g_0_x_xxxxx_xy[k] = g_xxxx_xy[k] - g_0_x_xxxx_xy[k] * ab_x + g_0_x_xxxx_xxy[k];

                g_0_x_xxxxx_xz[k] = g_xxxx_xz[k] - g_0_x_xxxx_xz[k] * ab_x + g_0_x_xxxx_xxz[k];

                g_0_x_xxxxx_yy[k] = g_xxxx_yy[k] - g_0_x_xxxx_yy[k] * ab_x + g_0_x_xxxx_xyy[k];

                g_0_x_xxxxx_yz[k] = g_xxxx_yz[k] - g_0_x_xxxx_yz[k] * ab_x + g_0_x_xxxx_xyz[k];

                g_0_x_xxxxx_zz[k] = g_xxxx_zz[k] - g_0_x_xxxx_zz[k] * ab_x + g_0_x_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxy_xx = cbuffer.data(hd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxy_xy = cbuffer.data(hd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxy_xz = cbuffer.data(hd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxy_yy = cbuffer.data(hd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxy_yz = cbuffer.data(hd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxy_zz = cbuffer.data(hd_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xx, g_0_x_xxxx_xxy, g_0_x_xxxx_xy, g_0_x_xxxx_xyy, g_0_x_xxxx_xyz, g_0_x_xxxx_xz, g_0_x_xxxx_yy, g_0_x_xxxx_yyy, g_0_x_xxxx_yyz, g_0_x_xxxx_yz, g_0_x_xxxx_yzz, g_0_x_xxxx_zz, g_0_x_xxxxy_xx, g_0_x_xxxxy_xy, g_0_x_xxxxy_xz, g_0_x_xxxxy_yy, g_0_x_xxxxy_yz, g_0_x_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxy_xx[k] = -g_0_x_xxxx_xx[k] * ab_y + g_0_x_xxxx_xxy[k];

                g_0_x_xxxxy_xy[k] = -g_0_x_xxxx_xy[k] * ab_y + g_0_x_xxxx_xyy[k];

                g_0_x_xxxxy_xz[k] = -g_0_x_xxxx_xz[k] * ab_y + g_0_x_xxxx_xyz[k];

                g_0_x_xxxxy_yy[k] = -g_0_x_xxxx_yy[k] * ab_y + g_0_x_xxxx_yyy[k];

                g_0_x_xxxxy_yz[k] = -g_0_x_xxxx_yz[k] * ab_y + g_0_x_xxxx_yyz[k];

                g_0_x_xxxxy_zz[k] = -g_0_x_xxxx_zz[k] * ab_y + g_0_x_xxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxz_xx = cbuffer.data(hd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxz_xy = cbuffer.data(hd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxz_xz = cbuffer.data(hd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxz_yy = cbuffer.data(hd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxz_yz = cbuffer.data(hd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxz_zz = cbuffer.data(hd_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xx, g_0_x_xxxx_xxz, g_0_x_xxxx_xy, g_0_x_xxxx_xyz, g_0_x_xxxx_xz, g_0_x_xxxx_xzz, g_0_x_xxxx_yy, g_0_x_xxxx_yyz, g_0_x_xxxx_yz, g_0_x_xxxx_yzz, g_0_x_xxxx_zz, g_0_x_xxxx_zzz, g_0_x_xxxxz_xx, g_0_x_xxxxz_xy, g_0_x_xxxxz_xz, g_0_x_xxxxz_yy, g_0_x_xxxxz_yz, g_0_x_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxz_xx[k] = -g_0_x_xxxx_xx[k] * ab_z + g_0_x_xxxx_xxz[k];

                g_0_x_xxxxz_xy[k] = -g_0_x_xxxx_xy[k] * ab_z + g_0_x_xxxx_xyz[k];

                g_0_x_xxxxz_xz[k] = -g_0_x_xxxx_xz[k] * ab_z + g_0_x_xxxx_xzz[k];

                g_0_x_xxxxz_yy[k] = -g_0_x_xxxx_yy[k] * ab_z + g_0_x_xxxx_yyz[k];

                g_0_x_xxxxz_yz[k] = -g_0_x_xxxx_yz[k] * ab_z + g_0_x_xxxx_yzz[k];

                g_0_x_xxxxz_zz[k] = -g_0_x_xxxx_zz[k] * ab_z + g_0_x_xxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyy_xx = cbuffer.data(hd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxyy_xy = cbuffer.data(hd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxyy_xz = cbuffer.data(hd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxyy_yy = cbuffer.data(hd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxyy_yz = cbuffer.data(hd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxyy_zz = cbuffer.data(hd_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxy_xx, g_0_x_xxxy_xxy, g_0_x_xxxy_xy, g_0_x_xxxy_xyy, g_0_x_xxxy_xyz, g_0_x_xxxy_xz, g_0_x_xxxy_yy, g_0_x_xxxy_yyy, g_0_x_xxxy_yyz, g_0_x_xxxy_yz, g_0_x_xxxy_yzz, g_0_x_xxxy_zz, g_0_x_xxxyy_xx, g_0_x_xxxyy_xy, g_0_x_xxxyy_xz, g_0_x_xxxyy_yy, g_0_x_xxxyy_yz, g_0_x_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyy_xx[k] = -g_0_x_xxxy_xx[k] * ab_y + g_0_x_xxxy_xxy[k];

                g_0_x_xxxyy_xy[k] = -g_0_x_xxxy_xy[k] * ab_y + g_0_x_xxxy_xyy[k];

                g_0_x_xxxyy_xz[k] = -g_0_x_xxxy_xz[k] * ab_y + g_0_x_xxxy_xyz[k];

                g_0_x_xxxyy_yy[k] = -g_0_x_xxxy_yy[k] * ab_y + g_0_x_xxxy_yyy[k];

                g_0_x_xxxyy_yz[k] = -g_0_x_xxxy_yz[k] * ab_y + g_0_x_xxxy_yyz[k];

                g_0_x_xxxyy_zz[k] = -g_0_x_xxxy_zz[k] * ab_y + g_0_x_xxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyz_xx = cbuffer.data(hd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxyz_xy = cbuffer.data(hd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxyz_xz = cbuffer.data(hd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxyz_yy = cbuffer.data(hd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxyz_yz = cbuffer.data(hd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxyz_zz = cbuffer.data(hd_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyz_xx, g_0_x_xxxyz_xy, g_0_x_xxxyz_xz, g_0_x_xxxyz_yy, g_0_x_xxxyz_yz, g_0_x_xxxyz_zz, g_0_x_xxxz_xx, g_0_x_xxxz_xxy, g_0_x_xxxz_xy, g_0_x_xxxz_xyy, g_0_x_xxxz_xyz, g_0_x_xxxz_xz, g_0_x_xxxz_yy, g_0_x_xxxz_yyy, g_0_x_xxxz_yyz, g_0_x_xxxz_yz, g_0_x_xxxz_yzz, g_0_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyz_xx[k] = -g_0_x_xxxz_xx[k] * ab_y + g_0_x_xxxz_xxy[k];

                g_0_x_xxxyz_xy[k] = -g_0_x_xxxz_xy[k] * ab_y + g_0_x_xxxz_xyy[k];

                g_0_x_xxxyz_xz[k] = -g_0_x_xxxz_xz[k] * ab_y + g_0_x_xxxz_xyz[k];

                g_0_x_xxxyz_yy[k] = -g_0_x_xxxz_yy[k] * ab_y + g_0_x_xxxz_yyy[k];

                g_0_x_xxxyz_yz[k] = -g_0_x_xxxz_yz[k] * ab_y + g_0_x_xxxz_yyz[k];

                g_0_x_xxxyz_zz[k] = -g_0_x_xxxz_zz[k] * ab_y + g_0_x_xxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzz_xx = cbuffer.data(hd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxzz_xy = cbuffer.data(hd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxzz_xz = cbuffer.data(hd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxzz_yy = cbuffer.data(hd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxzz_yz = cbuffer.data(hd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxzz_zz = cbuffer.data(hd_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxz_xx, g_0_x_xxxz_xxz, g_0_x_xxxz_xy, g_0_x_xxxz_xyz, g_0_x_xxxz_xz, g_0_x_xxxz_xzz, g_0_x_xxxz_yy, g_0_x_xxxz_yyz, g_0_x_xxxz_yz, g_0_x_xxxz_yzz, g_0_x_xxxz_zz, g_0_x_xxxz_zzz, g_0_x_xxxzz_xx, g_0_x_xxxzz_xy, g_0_x_xxxzz_xz, g_0_x_xxxzz_yy, g_0_x_xxxzz_yz, g_0_x_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzz_xx[k] = -g_0_x_xxxz_xx[k] * ab_z + g_0_x_xxxz_xxz[k];

                g_0_x_xxxzz_xy[k] = -g_0_x_xxxz_xy[k] * ab_z + g_0_x_xxxz_xyz[k];

                g_0_x_xxxzz_xz[k] = -g_0_x_xxxz_xz[k] * ab_z + g_0_x_xxxz_xzz[k];

                g_0_x_xxxzz_yy[k] = -g_0_x_xxxz_yy[k] * ab_z + g_0_x_xxxz_yyz[k];

                g_0_x_xxxzz_yz[k] = -g_0_x_xxxz_yz[k] * ab_z + g_0_x_xxxz_yzz[k];

                g_0_x_xxxzz_zz[k] = -g_0_x_xxxz_zz[k] * ab_z + g_0_x_xxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyy_xx = cbuffer.data(hd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyyy_xy = cbuffer.data(hd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyyy_xz = cbuffer.data(hd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyyy_yy = cbuffer.data(hd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyyy_yz = cbuffer.data(hd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyyy_zz = cbuffer.data(hd_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyy_xx, g_0_x_xxyy_xxy, g_0_x_xxyy_xy, g_0_x_xxyy_xyy, g_0_x_xxyy_xyz, g_0_x_xxyy_xz, g_0_x_xxyy_yy, g_0_x_xxyy_yyy, g_0_x_xxyy_yyz, g_0_x_xxyy_yz, g_0_x_xxyy_yzz, g_0_x_xxyy_zz, g_0_x_xxyyy_xx, g_0_x_xxyyy_xy, g_0_x_xxyyy_xz, g_0_x_xxyyy_yy, g_0_x_xxyyy_yz, g_0_x_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyy_xx[k] = -g_0_x_xxyy_xx[k] * ab_y + g_0_x_xxyy_xxy[k];

                g_0_x_xxyyy_xy[k] = -g_0_x_xxyy_xy[k] * ab_y + g_0_x_xxyy_xyy[k];

                g_0_x_xxyyy_xz[k] = -g_0_x_xxyy_xz[k] * ab_y + g_0_x_xxyy_xyz[k];

                g_0_x_xxyyy_yy[k] = -g_0_x_xxyy_yy[k] * ab_y + g_0_x_xxyy_yyy[k];

                g_0_x_xxyyy_yz[k] = -g_0_x_xxyy_yz[k] * ab_y + g_0_x_xxyy_yyz[k];

                g_0_x_xxyyy_zz[k] = -g_0_x_xxyy_zz[k] * ab_y + g_0_x_xxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyz_xx = cbuffer.data(hd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxyyz_xy = cbuffer.data(hd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxyyz_xz = cbuffer.data(hd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxyyz_yy = cbuffer.data(hd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyyz_yz = cbuffer.data(hd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyyz_zz = cbuffer.data(hd_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyz_xx, g_0_x_xxyyz_xy, g_0_x_xxyyz_xz, g_0_x_xxyyz_yy, g_0_x_xxyyz_yz, g_0_x_xxyyz_zz, g_0_x_xxyz_xx, g_0_x_xxyz_xxy, g_0_x_xxyz_xy, g_0_x_xxyz_xyy, g_0_x_xxyz_xyz, g_0_x_xxyz_xz, g_0_x_xxyz_yy, g_0_x_xxyz_yyy, g_0_x_xxyz_yyz, g_0_x_xxyz_yz, g_0_x_xxyz_yzz, g_0_x_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyz_xx[k] = -g_0_x_xxyz_xx[k] * ab_y + g_0_x_xxyz_xxy[k];

                g_0_x_xxyyz_xy[k] = -g_0_x_xxyz_xy[k] * ab_y + g_0_x_xxyz_xyy[k];

                g_0_x_xxyyz_xz[k] = -g_0_x_xxyz_xz[k] * ab_y + g_0_x_xxyz_xyz[k];

                g_0_x_xxyyz_yy[k] = -g_0_x_xxyz_yy[k] * ab_y + g_0_x_xxyz_yyy[k];

                g_0_x_xxyyz_yz[k] = -g_0_x_xxyz_yz[k] * ab_y + g_0_x_xxyz_yyz[k];

                g_0_x_xxyyz_zz[k] = -g_0_x_xxyz_zz[k] * ab_y + g_0_x_xxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzz_xx = cbuffer.data(hd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyzz_xy = cbuffer.data(hd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxyzz_xz = cbuffer.data(hd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxyzz_yy = cbuffer.data(hd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxyzz_yz = cbuffer.data(hd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxyzz_zz = cbuffer.data(hd_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzz_xx, g_0_x_xxyzz_xy, g_0_x_xxyzz_xz, g_0_x_xxyzz_yy, g_0_x_xxyzz_yz, g_0_x_xxyzz_zz, g_0_x_xxzz_xx, g_0_x_xxzz_xxy, g_0_x_xxzz_xy, g_0_x_xxzz_xyy, g_0_x_xxzz_xyz, g_0_x_xxzz_xz, g_0_x_xxzz_yy, g_0_x_xxzz_yyy, g_0_x_xxzz_yyz, g_0_x_xxzz_yz, g_0_x_xxzz_yzz, g_0_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzz_xx[k] = -g_0_x_xxzz_xx[k] * ab_y + g_0_x_xxzz_xxy[k];

                g_0_x_xxyzz_xy[k] = -g_0_x_xxzz_xy[k] * ab_y + g_0_x_xxzz_xyy[k];

                g_0_x_xxyzz_xz[k] = -g_0_x_xxzz_xz[k] * ab_y + g_0_x_xxzz_xyz[k];

                g_0_x_xxyzz_yy[k] = -g_0_x_xxzz_yy[k] * ab_y + g_0_x_xxzz_yyy[k];

                g_0_x_xxyzz_yz[k] = -g_0_x_xxzz_yz[k] * ab_y + g_0_x_xxzz_yyz[k];

                g_0_x_xxyzz_zz[k] = -g_0_x_xxzz_zz[k] * ab_y + g_0_x_xxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzz_xx = cbuffer.data(hd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxzzz_xy = cbuffer.data(hd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxzzz_xz = cbuffer.data(hd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxzzz_yy = cbuffer.data(hd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxzzz_yz = cbuffer.data(hd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxzzz_zz = cbuffer.data(hd_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzz_xx, g_0_x_xxzz_xxz, g_0_x_xxzz_xy, g_0_x_xxzz_xyz, g_0_x_xxzz_xz, g_0_x_xxzz_xzz, g_0_x_xxzz_yy, g_0_x_xxzz_yyz, g_0_x_xxzz_yz, g_0_x_xxzz_yzz, g_0_x_xxzz_zz, g_0_x_xxzz_zzz, g_0_x_xxzzz_xx, g_0_x_xxzzz_xy, g_0_x_xxzzz_xz, g_0_x_xxzzz_yy, g_0_x_xxzzz_yz, g_0_x_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzz_xx[k] = -g_0_x_xxzz_xx[k] * ab_z + g_0_x_xxzz_xxz[k];

                g_0_x_xxzzz_xy[k] = -g_0_x_xxzz_xy[k] * ab_z + g_0_x_xxzz_xyz[k];

                g_0_x_xxzzz_xz[k] = -g_0_x_xxzz_xz[k] * ab_z + g_0_x_xxzz_xzz[k];

                g_0_x_xxzzz_yy[k] = -g_0_x_xxzz_yy[k] * ab_z + g_0_x_xxzz_yyz[k];

                g_0_x_xxzzz_yz[k] = -g_0_x_xxzz_yz[k] * ab_z + g_0_x_xxzz_yzz[k];

                g_0_x_xxzzz_zz[k] = -g_0_x_xxzz_zz[k] * ab_z + g_0_x_xxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyy_xx = cbuffer.data(hd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyyyy_xy = cbuffer.data(hd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyyyy_xz = cbuffer.data(hd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyyyy_yy = cbuffer.data(hd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyyy_yz = cbuffer.data(hd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyyy_zz = cbuffer.data(hd_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyy_xx, g_0_x_xyyy_xxy, g_0_x_xyyy_xy, g_0_x_xyyy_xyy, g_0_x_xyyy_xyz, g_0_x_xyyy_xz, g_0_x_xyyy_yy, g_0_x_xyyy_yyy, g_0_x_xyyy_yyz, g_0_x_xyyy_yz, g_0_x_xyyy_yzz, g_0_x_xyyy_zz, g_0_x_xyyyy_xx, g_0_x_xyyyy_xy, g_0_x_xyyyy_xz, g_0_x_xyyyy_yy, g_0_x_xyyyy_yz, g_0_x_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyy_xx[k] = -g_0_x_xyyy_xx[k] * ab_y + g_0_x_xyyy_xxy[k];

                g_0_x_xyyyy_xy[k] = -g_0_x_xyyy_xy[k] * ab_y + g_0_x_xyyy_xyy[k];

                g_0_x_xyyyy_xz[k] = -g_0_x_xyyy_xz[k] * ab_y + g_0_x_xyyy_xyz[k];

                g_0_x_xyyyy_yy[k] = -g_0_x_xyyy_yy[k] * ab_y + g_0_x_xyyy_yyy[k];

                g_0_x_xyyyy_yz[k] = -g_0_x_xyyy_yz[k] * ab_y + g_0_x_xyyy_yyz[k];

                g_0_x_xyyyy_zz[k] = -g_0_x_xyyy_zz[k] * ab_y + g_0_x_xyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyz_xx = cbuffer.data(hd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyyz_xy = cbuffer.data(hd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyyz_xz = cbuffer.data(hd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyyyz_yy = cbuffer.data(hd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyyyz_yz = cbuffer.data(hd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyyz_zz = cbuffer.data(hd_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyz_xx, g_0_x_xyyyz_xy, g_0_x_xyyyz_xz, g_0_x_xyyyz_yy, g_0_x_xyyyz_yz, g_0_x_xyyyz_zz, g_0_x_xyyz_xx, g_0_x_xyyz_xxy, g_0_x_xyyz_xy, g_0_x_xyyz_xyy, g_0_x_xyyz_xyz, g_0_x_xyyz_xz, g_0_x_xyyz_yy, g_0_x_xyyz_yyy, g_0_x_xyyz_yyz, g_0_x_xyyz_yz, g_0_x_xyyz_yzz, g_0_x_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyz_xx[k] = -g_0_x_xyyz_xx[k] * ab_y + g_0_x_xyyz_xxy[k];

                g_0_x_xyyyz_xy[k] = -g_0_x_xyyz_xy[k] * ab_y + g_0_x_xyyz_xyy[k];

                g_0_x_xyyyz_xz[k] = -g_0_x_xyyz_xz[k] * ab_y + g_0_x_xyyz_xyz[k];

                g_0_x_xyyyz_yy[k] = -g_0_x_xyyz_yy[k] * ab_y + g_0_x_xyyz_yyy[k];

                g_0_x_xyyyz_yz[k] = -g_0_x_xyyz_yz[k] * ab_y + g_0_x_xyyz_yyz[k];

                g_0_x_xyyyz_zz[k] = -g_0_x_xyyz_zz[k] * ab_y + g_0_x_xyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzz_xx = cbuffer.data(hd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyzz_xy = cbuffer.data(hd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyzz_xz = cbuffer.data(hd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyyzz_yy = cbuffer.data(hd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyzz_yz = cbuffer.data(hd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyzz_zz = cbuffer.data(hd_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzz_xx, g_0_x_xyyzz_xy, g_0_x_xyyzz_xz, g_0_x_xyyzz_yy, g_0_x_xyyzz_yz, g_0_x_xyyzz_zz, g_0_x_xyzz_xx, g_0_x_xyzz_xxy, g_0_x_xyzz_xy, g_0_x_xyzz_xyy, g_0_x_xyzz_xyz, g_0_x_xyzz_xz, g_0_x_xyzz_yy, g_0_x_xyzz_yyy, g_0_x_xyzz_yyz, g_0_x_xyzz_yz, g_0_x_xyzz_yzz, g_0_x_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzz_xx[k] = -g_0_x_xyzz_xx[k] * ab_y + g_0_x_xyzz_xxy[k];

                g_0_x_xyyzz_xy[k] = -g_0_x_xyzz_xy[k] * ab_y + g_0_x_xyzz_xyy[k];

                g_0_x_xyyzz_xz[k] = -g_0_x_xyzz_xz[k] * ab_y + g_0_x_xyzz_xyz[k];

                g_0_x_xyyzz_yy[k] = -g_0_x_xyzz_yy[k] * ab_y + g_0_x_xyzz_yyy[k];

                g_0_x_xyyzz_yz[k] = -g_0_x_xyzz_yz[k] * ab_y + g_0_x_xyzz_yyz[k];

                g_0_x_xyyzz_zz[k] = -g_0_x_xyzz_zz[k] * ab_y + g_0_x_xyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzz_xx = cbuffer.data(hd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyzzz_xy = cbuffer.data(hd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyzzz_xz = cbuffer.data(hd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyzzz_yy = cbuffer.data(hd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyzzz_yz = cbuffer.data(hd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyzzz_zz = cbuffer.data(hd_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzz_xx, g_0_x_xyzzz_xy, g_0_x_xyzzz_xz, g_0_x_xyzzz_yy, g_0_x_xyzzz_yz, g_0_x_xyzzz_zz, g_0_x_xzzz_xx, g_0_x_xzzz_xxy, g_0_x_xzzz_xy, g_0_x_xzzz_xyy, g_0_x_xzzz_xyz, g_0_x_xzzz_xz, g_0_x_xzzz_yy, g_0_x_xzzz_yyy, g_0_x_xzzz_yyz, g_0_x_xzzz_yz, g_0_x_xzzz_yzz, g_0_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzz_xx[k] = -g_0_x_xzzz_xx[k] * ab_y + g_0_x_xzzz_xxy[k];

                g_0_x_xyzzz_xy[k] = -g_0_x_xzzz_xy[k] * ab_y + g_0_x_xzzz_xyy[k];

                g_0_x_xyzzz_xz[k] = -g_0_x_xzzz_xz[k] * ab_y + g_0_x_xzzz_xyz[k];

                g_0_x_xyzzz_yy[k] = -g_0_x_xzzz_yy[k] * ab_y + g_0_x_xzzz_yyy[k];

                g_0_x_xyzzz_yz[k] = -g_0_x_xzzz_yz[k] * ab_y + g_0_x_xzzz_yyz[k];

                g_0_x_xyzzz_zz[k] = -g_0_x_xzzz_zz[k] * ab_y + g_0_x_xzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzz_xx = cbuffer.data(hd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzzzz_xy = cbuffer.data(hd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzzzz_xz = cbuffer.data(hd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzzzz_yy = cbuffer.data(hd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzzzz_yz = cbuffer.data(hd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzzzz_zz = cbuffer.data(hd_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzz_xx, g_0_x_xzzz_xxz, g_0_x_xzzz_xy, g_0_x_xzzz_xyz, g_0_x_xzzz_xz, g_0_x_xzzz_xzz, g_0_x_xzzz_yy, g_0_x_xzzz_yyz, g_0_x_xzzz_yz, g_0_x_xzzz_yzz, g_0_x_xzzz_zz, g_0_x_xzzz_zzz, g_0_x_xzzzz_xx, g_0_x_xzzzz_xy, g_0_x_xzzzz_xz, g_0_x_xzzzz_yy, g_0_x_xzzzz_yz, g_0_x_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzz_xx[k] = -g_0_x_xzzz_xx[k] * ab_z + g_0_x_xzzz_xxz[k];

                g_0_x_xzzzz_xy[k] = -g_0_x_xzzz_xy[k] * ab_z + g_0_x_xzzz_xyz[k];

                g_0_x_xzzzz_xz[k] = -g_0_x_xzzz_xz[k] * ab_z + g_0_x_xzzz_xzz[k];

                g_0_x_xzzzz_yy[k] = -g_0_x_xzzz_yy[k] * ab_z + g_0_x_xzzz_yyz[k];

                g_0_x_xzzzz_yz[k] = -g_0_x_xzzz_yz[k] * ab_z + g_0_x_xzzz_yzz[k];

                g_0_x_xzzzz_zz[k] = -g_0_x_xzzz_zz[k] * ab_z + g_0_x_xzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyy_xx = cbuffer.data(hd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyyyy_xy = cbuffer.data(hd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyyyy_xz = cbuffer.data(hd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyyyy_yy = cbuffer.data(hd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyyyy_yz = cbuffer.data(hd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyyyy_zz = cbuffer.data(hd_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_xx, g_0_x_yyyy_xxy, g_0_x_yyyy_xy, g_0_x_yyyy_xyy, g_0_x_yyyy_xyz, g_0_x_yyyy_xz, g_0_x_yyyy_yy, g_0_x_yyyy_yyy, g_0_x_yyyy_yyz, g_0_x_yyyy_yz, g_0_x_yyyy_yzz, g_0_x_yyyy_zz, g_0_x_yyyyy_xx, g_0_x_yyyyy_xy, g_0_x_yyyyy_xz, g_0_x_yyyyy_yy, g_0_x_yyyyy_yz, g_0_x_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyy_xx[k] = -g_0_x_yyyy_xx[k] * ab_y + g_0_x_yyyy_xxy[k];

                g_0_x_yyyyy_xy[k] = -g_0_x_yyyy_xy[k] * ab_y + g_0_x_yyyy_xyy[k];

                g_0_x_yyyyy_xz[k] = -g_0_x_yyyy_xz[k] * ab_y + g_0_x_yyyy_xyz[k];

                g_0_x_yyyyy_yy[k] = -g_0_x_yyyy_yy[k] * ab_y + g_0_x_yyyy_yyy[k];

                g_0_x_yyyyy_yz[k] = -g_0_x_yyyy_yz[k] * ab_y + g_0_x_yyyy_yyz[k];

                g_0_x_yyyyy_zz[k] = -g_0_x_yyyy_zz[k] * ab_y + g_0_x_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyz_xx = cbuffer.data(hd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyyyz_xy = cbuffer.data(hd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyyyz_xz = cbuffer.data(hd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyyyz_yy = cbuffer.data(hd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyyyz_yz = cbuffer.data(hd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyyyz_zz = cbuffer.data(hd_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyz_xx, g_0_x_yyyyz_xy, g_0_x_yyyyz_xz, g_0_x_yyyyz_yy, g_0_x_yyyyz_yz, g_0_x_yyyyz_zz, g_0_x_yyyz_xx, g_0_x_yyyz_xxy, g_0_x_yyyz_xy, g_0_x_yyyz_xyy, g_0_x_yyyz_xyz, g_0_x_yyyz_xz, g_0_x_yyyz_yy, g_0_x_yyyz_yyy, g_0_x_yyyz_yyz, g_0_x_yyyz_yz, g_0_x_yyyz_yzz, g_0_x_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyz_xx[k] = -g_0_x_yyyz_xx[k] * ab_y + g_0_x_yyyz_xxy[k];

                g_0_x_yyyyz_xy[k] = -g_0_x_yyyz_xy[k] * ab_y + g_0_x_yyyz_xyy[k];

                g_0_x_yyyyz_xz[k] = -g_0_x_yyyz_xz[k] * ab_y + g_0_x_yyyz_xyz[k];

                g_0_x_yyyyz_yy[k] = -g_0_x_yyyz_yy[k] * ab_y + g_0_x_yyyz_yyy[k];

                g_0_x_yyyyz_yz[k] = -g_0_x_yyyz_yz[k] * ab_y + g_0_x_yyyz_yyz[k];

                g_0_x_yyyyz_zz[k] = -g_0_x_yyyz_zz[k] * ab_y + g_0_x_yyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzz_xx = cbuffer.data(hd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyyzz_xy = cbuffer.data(hd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyyzz_xz = cbuffer.data(hd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyyzz_yy = cbuffer.data(hd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyyzz_yz = cbuffer.data(hd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyyzz_zz = cbuffer.data(hd_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzz_xx, g_0_x_yyyzz_xy, g_0_x_yyyzz_xz, g_0_x_yyyzz_yy, g_0_x_yyyzz_yz, g_0_x_yyyzz_zz, g_0_x_yyzz_xx, g_0_x_yyzz_xxy, g_0_x_yyzz_xy, g_0_x_yyzz_xyy, g_0_x_yyzz_xyz, g_0_x_yyzz_xz, g_0_x_yyzz_yy, g_0_x_yyzz_yyy, g_0_x_yyzz_yyz, g_0_x_yyzz_yz, g_0_x_yyzz_yzz, g_0_x_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzz_xx[k] = -g_0_x_yyzz_xx[k] * ab_y + g_0_x_yyzz_xxy[k];

                g_0_x_yyyzz_xy[k] = -g_0_x_yyzz_xy[k] * ab_y + g_0_x_yyzz_xyy[k];

                g_0_x_yyyzz_xz[k] = -g_0_x_yyzz_xz[k] * ab_y + g_0_x_yyzz_xyz[k];

                g_0_x_yyyzz_yy[k] = -g_0_x_yyzz_yy[k] * ab_y + g_0_x_yyzz_yyy[k];

                g_0_x_yyyzz_yz[k] = -g_0_x_yyzz_yz[k] * ab_y + g_0_x_yyzz_yyz[k];

                g_0_x_yyyzz_zz[k] = -g_0_x_yyzz_zz[k] * ab_y + g_0_x_yyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzz_xx = cbuffer.data(hd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyzzz_xy = cbuffer.data(hd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyzzz_xz = cbuffer.data(hd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyzzz_yy = cbuffer.data(hd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyzzz_yz = cbuffer.data(hd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyzzz_zz = cbuffer.data(hd_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzz_xx, g_0_x_yyzzz_xy, g_0_x_yyzzz_xz, g_0_x_yyzzz_yy, g_0_x_yyzzz_yz, g_0_x_yyzzz_zz, g_0_x_yzzz_xx, g_0_x_yzzz_xxy, g_0_x_yzzz_xy, g_0_x_yzzz_xyy, g_0_x_yzzz_xyz, g_0_x_yzzz_xz, g_0_x_yzzz_yy, g_0_x_yzzz_yyy, g_0_x_yzzz_yyz, g_0_x_yzzz_yz, g_0_x_yzzz_yzz, g_0_x_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzz_xx[k] = -g_0_x_yzzz_xx[k] * ab_y + g_0_x_yzzz_xxy[k];

                g_0_x_yyzzz_xy[k] = -g_0_x_yzzz_xy[k] * ab_y + g_0_x_yzzz_xyy[k];

                g_0_x_yyzzz_xz[k] = -g_0_x_yzzz_xz[k] * ab_y + g_0_x_yzzz_xyz[k];

                g_0_x_yyzzz_yy[k] = -g_0_x_yzzz_yy[k] * ab_y + g_0_x_yzzz_yyy[k];

                g_0_x_yyzzz_yz[k] = -g_0_x_yzzz_yz[k] * ab_y + g_0_x_yzzz_yyz[k];

                g_0_x_yyzzz_zz[k] = -g_0_x_yzzz_zz[k] * ab_y + g_0_x_yzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzz_xx = cbuffer.data(hd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yzzzz_xy = cbuffer.data(hd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yzzzz_xz = cbuffer.data(hd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yzzzz_yy = cbuffer.data(hd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yzzzz_yz = cbuffer.data(hd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yzzzz_zz = cbuffer.data(hd_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzz_xx, g_0_x_yzzzz_xy, g_0_x_yzzzz_xz, g_0_x_yzzzz_yy, g_0_x_yzzzz_yz, g_0_x_yzzzz_zz, g_0_x_zzzz_xx, g_0_x_zzzz_xxy, g_0_x_zzzz_xy, g_0_x_zzzz_xyy, g_0_x_zzzz_xyz, g_0_x_zzzz_xz, g_0_x_zzzz_yy, g_0_x_zzzz_yyy, g_0_x_zzzz_yyz, g_0_x_zzzz_yz, g_0_x_zzzz_yzz, g_0_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzz_xx[k] = -g_0_x_zzzz_xx[k] * ab_y + g_0_x_zzzz_xxy[k];

                g_0_x_yzzzz_xy[k] = -g_0_x_zzzz_xy[k] * ab_y + g_0_x_zzzz_xyy[k];

                g_0_x_yzzzz_xz[k] = -g_0_x_zzzz_xz[k] * ab_y + g_0_x_zzzz_xyz[k];

                g_0_x_yzzzz_yy[k] = -g_0_x_zzzz_yy[k] * ab_y + g_0_x_zzzz_yyy[k];

                g_0_x_yzzzz_yz[k] = -g_0_x_zzzz_yz[k] * ab_y + g_0_x_zzzz_yyz[k];

                g_0_x_yzzzz_zz[k] = -g_0_x_zzzz_zz[k] * ab_y + g_0_x_zzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzz_xx = cbuffer.data(hd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_zzzzz_xy = cbuffer.data(hd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_zzzzz_xz = cbuffer.data(hd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_zzzzz_yy = cbuffer.data(hd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_zzzzz_yz = cbuffer.data(hd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_zzzzz_zz = cbuffer.data(hd_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_xx, g_0_x_zzzz_xxz, g_0_x_zzzz_xy, g_0_x_zzzz_xyz, g_0_x_zzzz_xz, g_0_x_zzzz_xzz, g_0_x_zzzz_yy, g_0_x_zzzz_yyz, g_0_x_zzzz_yz, g_0_x_zzzz_yzz, g_0_x_zzzz_zz, g_0_x_zzzz_zzz, g_0_x_zzzzz_xx, g_0_x_zzzzz_xy, g_0_x_zzzzz_xz, g_0_x_zzzzz_yy, g_0_x_zzzzz_yz, g_0_x_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzz_xx[k] = -g_0_x_zzzz_xx[k] * ab_z + g_0_x_zzzz_xxz[k];

                g_0_x_zzzzz_xy[k] = -g_0_x_zzzz_xy[k] * ab_z + g_0_x_zzzz_xyz[k];

                g_0_x_zzzzz_xz[k] = -g_0_x_zzzz_xz[k] * ab_z + g_0_x_zzzz_xzz[k];

                g_0_x_zzzzz_yy[k] = -g_0_x_zzzz_yy[k] * ab_z + g_0_x_zzzz_yyz[k];

                g_0_x_zzzzz_yz[k] = -g_0_x_zzzz_yz[k] * ab_z + g_0_x_zzzz_yzz[k];

                g_0_x_zzzzz_zz[k] = -g_0_x_zzzz_zz[k] * ab_z + g_0_x_zzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxx_xx = cbuffer.data(hd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxxxx_xy = cbuffer.data(hd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxxxx_xz = cbuffer.data(hd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xxxxx_yy = cbuffer.data(hd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xxxxx_yz = cbuffer.data(hd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xxxxx_zz = cbuffer.data(hd_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_xx, g_0_y_xxxx_xxx, g_0_y_xxxx_xxy, g_0_y_xxxx_xxz, g_0_y_xxxx_xy, g_0_y_xxxx_xyy, g_0_y_xxxx_xyz, g_0_y_xxxx_xz, g_0_y_xxxx_xzz, g_0_y_xxxx_yy, g_0_y_xxxx_yz, g_0_y_xxxx_zz, g_0_y_xxxxx_xx, g_0_y_xxxxx_xy, g_0_y_xxxxx_xz, g_0_y_xxxxx_yy, g_0_y_xxxxx_yz, g_0_y_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxx_xx[k] = -g_0_y_xxxx_xx[k] * ab_x + g_0_y_xxxx_xxx[k];

                g_0_y_xxxxx_xy[k] = -g_0_y_xxxx_xy[k] * ab_x + g_0_y_xxxx_xxy[k];

                g_0_y_xxxxx_xz[k] = -g_0_y_xxxx_xz[k] * ab_x + g_0_y_xxxx_xxz[k];

                g_0_y_xxxxx_yy[k] = -g_0_y_xxxx_yy[k] * ab_x + g_0_y_xxxx_xyy[k];

                g_0_y_xxxxx_yz[k] = -g_0_y_xxxx_yz[k] * ab_x + g_0_y_xxxx_xyz[k];

                g_0_y_xxxxx_zz[k] = -g_0_y_xxxx_zz[k] * ab_x + g_0_y_xxxx_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxy_xx = cbuffer.data(hd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xxxxy_xy = cbuffer.data(hd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xxxxy_xz = cbuffer.data(hd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xxxxy_yy = cbuffer.data(hd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xxxxy_yz = cbuffer.data(hd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xxxxy_zz = cbuffer.data(hd_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxy_xx, g_0_y_xxxxy_xy, g_0_y_xxxxy_xz, g_0_y_xxxxy_yy, g_0_y_xxxxy_yz, g_0_y_xxxxy_zz, g_0_y_xxxy_xx, g_0_y_xxxy_xxx, g_0_y_xxxy_xxy, g_0_y_xxxy_xxz, g_0_y_xxxy_xy, g_0_y_xxxy_xyy, g_0_y_xxxy_xyz, g_0_y_xxxy_xz, g_0_y_xxxy_xzz, g_0_y_xxxy_yy, g_0_y_xxxy_yz, g_0_y_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxy_xx[k] = -g_0_y_xxxy_xx[k] * ab_x + g_0_y_xxxy_xxx[k];

                g_0_y_xxxxy_xy[k] = -g_0_y_xxxy_xy[k] * ab_x + g_0_y_xxxy_xxy[k];

                g_0_y_xxxxy_xz[k] = -g_0_y_xxxy_xz[k] * ab_x + g_0_y_xxxy_xxz[k];

                g_0_y_xxxxy_yy[k] = -g_0_y_xxxy_yy[k] * ab_x + g_0_y_xxxy_xyy[k];

                g_0_y_xxxxy_yz[k] = -g_0_y_xxxy_yz[k] * ab_x + g_0_y_xxxy_xyz[k];

                g_0_y_xxxxy_zz[k] = -g_0_y_xxxy_zz[k] * ab_x + g_0_y_xxxy_xzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxz_xx = cbuffer.data(hd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xxxxz_xy = cbuffer.data(hd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xxxxz_xz = cbuffer.data(hd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xxxxz_yy = cbuffer.data(hd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xxxxz_yz = cbuffer.data(hd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xxxxz_zz = cbuffer.data(hd_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxz_xx, g_0_y_xxxxz_xy, g_0_y_xxxxz_xz, g_0_y_xxxxz_yy, g_0_y_xxxxz_yz, g_0_y_xxxxz_zz, g_0_y_xxxz_xx, g_0_y_xxxz_xxx, g_0_y_xxxz_xxy, g_0_y_xxxz_xxz, g_0_y_xxxz_xy, g_0_y_xxxz_xyy, g_0_y_xxxz_xyz, g_0_y_xxxz_xz, g_0_y_xxxz_xzz, g_0_y_xxxz_yy, g_0_y_xxxz_yz, g_0_y_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxz_xx[k] = -g_0_y_xxxz_xx[k] * ab_x + g_0_y_xxxz_xxx[k];

                g_0_y_xxxxz_xy[k] = -g_0_y_xxxz_xy[k] * ab_x + g_0_y_xxxz_xxy[k];

                g_0_y_xxxxz_xz[k] = -g_0_y_xxxz_xz[k] * ab_x + g_0_y_xxxz_xxz[k];

                g_0_y_xxxxz_yy[k] = -g_0_y_xxxz_yy[k] * ab_x + g_0_y_xxxz_xyy[k];

                g_0_y_xxxxz_yz[k] = -g_0_y_xxxz_yz[k] * ab_x + g_0_y_xxxz_xyz[k];

                g_0_y_xxxxz_zz[k] = -g_0_y_xxxz_zz[k] * ab_x + g_0_y_xxxz_xzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyy_xx = cbuffer.data(hd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xxxyy_xy = cbuffer.data(hd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xxxyy_xz = cbuffer.data(hd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xxxyy_yy = cbuffer.data(hd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xxxyy_yz = cbuffer.data(hd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xxxyy_zz = cbuffer.data(hd_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyy_xx, g_0_y_xxxyy_xy, g_0_y_xxxyy_xz, g_0_y_xxxyy_yy, g_0_y_xxxyy_yz, g_0_y_xxxyy_zz, g_0_y_xxyy_xx, g_0_y_xxyy_xxx, g_0_y_xxyy_xxy, g_0_y_xxyy_xxz, g_0_y_xxyy_xy, g_0_y_xxyy_xyy, g_0_y_xxyy_xyz, g_0_y_xxyy_xz, g_0_y_xxyy_xzz, g_0_y_xxyy_yy, g_0_y_xxyy_yz, g_0_y_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyy_xx[k] = -g_0_y_xxyy_xx[k] * ab_x + g_0_y_xxyy_xxx[k];

                g_0_y_xxxyy_xy[k] = -g_0_y_xxyy_xy[k] * ab_x + g_0_y_xxyy_xxy[k];

                g_0_y_xxxyy_xz[k] = -g_0_y_xxyy_xz[k] * ab_x + g_0_y_xxyy_xxz[k];

                g_0_y_xxxyy_yy[k] = -g_0_y_xxyy_yy[k] * ab_x + g_0_y_xxyy_xyy[k];

                g_0_y_xxxyy_yz[k] = -g_0_y_xxyy_yz[k] * ab_x + g_0_y_xxyy_xyz[k];

                g_0_y_xxxyy_zz[k] = -g_0_y_xxyy_zz[k] * ab_x + g_0_y_xxyy_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyz_xx = cbuffer.data(hd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxyz_xy = cbuffer.data(hd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxyz_xz = cbuffer.data(hd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxxyz_yy = cbuffer.data(hd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxxyz_yz = cbuffer.data(hd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxxyz_zz = cbuffer.data(hd_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyz_xx, g_0_y_xxxyz_xy, g_0_y_xxxyz_xz, g_0_y_xxxyz_yy, g_0_y_xxxyz_yz, g_0_y_xxxyz_zz, g_0_y_xxyz_xx, g_0_y_xxyz_xxx, g_0_y_xxyz_xxy, g_0_y_xxyz_xxz, g_0_y_xxyz_xy, g_0_y_xxyz_xyy, g_0_y_xxyz_xyz, g_0_y_xxyz_xz, g_0_y_xxyz_xzz, g_0_y_xxyz_yy, g_0_y_xxyz_yz, g_0_y_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyz_xx[k] = -g_0_y_xxyz_xx[k] * ab_x + g_0_y_xxyz_xxx[k];

                g_0_y_xxxyz_xy[k] = -g_0_y_xxyz_xy[k] * ab_x + g_0_y_xxyz_xxy[k];

                g_0_y_xxxyz_xz[k] = -g_0_y_xxyz_xz[k] * ab_x + g_0_y_xxyz_xxz[k];

                g_0_y_xxxyz_yy[k] = -g_0_y_xxyz_yy[k] * ab_x + g_0_y_xxyz_xyy[k];

                g_0_y_xxxyz_yz[k] = -g_0_y_xxyz_yz[k] * ab_x + g_0_y_xxyz_xyz[k];

                g_0_y_xxxyz_zz[k] = -g_0_y_xxyz_zz[k] * ab_x + g_0_y_xxyz_xzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzz_xx = cbuffer.data(hd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxxzz_xy = cbuffer.data(hd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxxzz_xz = cbuffer.data(hd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxxzz_yy = cbuffer.data(hd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxxzz_yz = cbuffer.data(hd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxxzz_zz = cbuffer.data(hd_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzz_xx, g_0_y_xxxzz_xy, g_0_y_xxxzz_xz, g_0_y_xxxzz_yy, g_0_y_xxxzz_yz, g_0_y_xxxzz_zz, g_0_y_xxzz_xx, g_0_y_xxzz_xxx, g_0_y_xxzz_xxy, g_0_y_xxzz_xxz, g_0_y_xxzz_xy, g_0_y_xxzz_xyy, g_0_y_xxzz_xyz, g_0_y_xxzz_xz, g_0_y_xxzz_xzz, g_0_y_xxzz_yy, g_0_y_xxzz_yz, g_0_y_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzz_xx[k] = -g_0_y_xxzz_xx[k] * ab_x + g_0_y_xxzz_xxx[k];

                g_0_y_xxxzz_xy[k] = -g_0_y_xxzz_xy[k] * ab_x + g_0_y_xxzz_xxy[k];

                g_0_y_xxxzz_xz[k] = -g_0_y_xxzz_xz[k] * ab_x + g_0_y_xxzz_xxz[k];

                g_0_y_xxxzz_yy[k] = -g_0_y_xxzz_yy[k] * ab_x + g_0_y_xxzz_xyy[k];

                g_0_y_xxxzz_yz[k] = -g_0_y_xxzz_yz[k] * ab_x + g_0_y_xxzz_xyz[k];

                g_0_y_xxxzz_zz[k] = -g_0_y_xxzz_zz[k] * ab_x + g_0_y_xxzz_xzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyy_xx = cbuffer.data(hd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxyyy_xy = cbuffer.data(hd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxyyy_xz = cbuffer.data(hd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxyyy_yy = cbuffer.data(hd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxyyy_yz = cbuffer.data(hd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxyyy_zz = cbuffer.data(hd_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyy_xx, g_0_y_xxyyy_xy, g_0_y_xxyyy_xz, g_0_y_xxyyy_yy, g_0_y_xxyyy_yz, g_0_y_xxyyy_zz, g_0_y_xyyy_xx, g_0_y_xyyy_xxx, g_0_y_xyyy_xxy, g_0_y_xyyy_xxz, g_0_y_xyyy_xy, g_0_y_xyyy_xyy, g_0_y_xyyy_xyz, g_0_y_xyyy_xz, g_0_y_xyyy_xzz, g_0_y_xyyy_yy, g_0_y_xyyy_yz, g_0_y_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyy_xx[k] = -g_0_y_xyyy_xx[k] * ab_x + g_0_y_xyyy_xxx[k];

                g_0_y_xxyyy_xy[k] = -g_0_y_xyyy_xy[k] * ab_x + g_0_y_xyyy_xxy[k];

                g_0_y_xxyyy_xz[k] = -g_0_y_xyyy_xz[k] * ab_x + g_0_y_xyyy_xxz[k];

                g_0_y_xxyyy_yy[k] = -g_0_y_xyyy_yy[k] * ab_x + g_0_y_xyyy_xyy[k];

                g_0_y_xxyyy_yz[k] = -g_0_y_xyyy_yz[k] * ab_x + g_0_y_xyyy_xyz[k];

                g_0_y_xxyyy_zz[k] = -g_0_y_xyyy_zz[k] * ab_x + g_0_y_xyyy_xzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyz_xx = cbuffer.data(hd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxyyz_xy = cbuffer.data(hd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxyyz_xz = cbuffer.data(hd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxyyz_yy = cbuffer.data(hd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxyyz_yz = cbuffer.data(hd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxyyz_zz = cbuffer.data(hd_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyz_xx, g_0_y_xxyyz_xy, g_0_y_xxyyz_xz, g_0_y_xxyyz_yy, g_0_y_xxyyz_yz, g_0_y_xxyyz_zz, g_0_y_xyyz_xx, g_0_y_xyyz_xxx, g_0_y_xyyz_xxy, g_0_y_xyyz_xxz, g_0_y_xyyz_xy, g_0_y_xyyz_xyy, g_0_y_xyyz_xyz, g_0_y_xyyz_xz, g_0_y_xyyz_xzz, g_0_y_xyyz_yy, g_0_y_xyyz_yz, g_0_y_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyz_xx[k] = -g_0_y_xyyz_xx[k] * ab_x + g_0_y_xyyz_xxx[k];

                g_0_y_xxyyz_xy[k] = -g_0_y_xyyz_xy[k] * ab_x + g_0_y_xyyz_xxy[k];

                g_0_y_xxyyz_xz[k] = -g_0_y_xyyz_xz[k] * ab_x + g_0_y_xyyz_xxz[k];

                g_0_y_xxyyz_yy[k] = -g_0_y_xyyz_yy[k] * ab_x + g_0_y_xyyz_xyy[k];

                g_0_y_xxyyz_yz[k] = -g_0_y_xyyz_yz[k] * ab_x + g_0_y_xyyz_xyz[k];

                g_0_y_xxyyz_zz[k] = -g_0_y_xyyz_zz[k] * ab_x + g_0_y_xyyz_xzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzz_xx = cbuffer.data(hd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxyzz_xy = cbuffer.data(hd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxyzz_xz = cbuffer.data(hd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxyzz_yy = cbuffer.data(hd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxyzz_yz = cbuffer.data(hd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxyzz_zz = cbuffer.data(hd_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzz_xx, g_0_y_xxyzz_xy, g_0_y_xxyzz_xz, g_0_y_xxyzz_yy, g_0_y_xxyzz_yz, g_0_y_xxyzz_zz, g_0_y_xyzz_xx, g_0_y_xyzz_xxx, g_0_y_xyzz_xxy, g_0_y_xyzz_xxz, g_0_y_xyzz_xy, g_0_y_xyzz_xyy, g_0_y_xyzz_xyz, g_0_y_xyzz_xz, g_0_y_xyzz_xzz, g_0_y_xyzz_yy, g_0_y_xyzz_yz, g_0_y_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzz_xx[k] = -g_0_y_xyzz_xx[k] * ab_x + g_0_y_xyzz_xxx[k];

                g_0_y_xxyzz_xy[k] = -g_0_y_xyzz_xy[k] * ab_x + g_0_y_xyzz_xxy[k];

                g_0_y_xxyzz_xz[k] = -g_0_y_xyzz_xz[k] * ab_x + g_0_y_xyzz_xxz[k];

                g_0_y_xxyzz_yy[k] = -g_0_y_xyzz_yy[k] * ab_x + g_0_y_xyzz_xyy[k];

                g_0_y_xxyzz_yz[k] = -g_0_y_xyzz_yz[k] * ab_x + g_0_y_xyzz_xyz[k];

                g_0_y_xxyzz_zz[k] = -g_0_y_xyzz_zz[k] * ab_x + g_0_y_xyzz_xzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzz_xx = cbuffer.data(hd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxzzz_xy = cbuffer.data(hd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxzzz_xz = cbuffer.data(hd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxzzz_yy = cbuffer.data(hd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxzzz_yz = cbuffer.data(hd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxzzz_zz = cbuffer.data(hd_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzz_xx, g_0_y_xxzzz_xy, g_0_y_xxzzz_xz, g_0_y_xxzzz_yy, g_0_y_xxzzz_yz, g_0_y_xxzzz_zz, g_0_y_xzzz_xx, g_0_y_xzzz_xxx, g_0_y_xzzz_xxy, g_0_y_xzzz_xxz, g_0_y_xzzz_xy, g_0_y_xzzz_xyy, g_0_y_xzzz_xyz, g_0_y_xzzz_xz, g_0_y_xzzz_xzz, g_0_y_xzzz_yy, g_0_y_xzzz_yz, g_0_y_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzz_xx[k] = -g_0_y_xzzz_xx[k] * ab_x + g_0_y_xzzz_xxx[k];

                g_0_y_xxzzz_xy[k] = -g_0_y_xzzz_xy[k] * ab_x + g_0_y_xzzz_xxy[k];

                g_0_y_xxzzz_xz[k] = -g_0_y_xzzz_xz[k] * ab_x + g_0_y_xzzz_xxz[k];

                g_0_y_xxzzz_yy[k] = -g_0_y_xzzz_yy[k] * ab_x + g_0_y_xzzz_xyy[k];

                g_0_y_xxzzz_yz[k] = -g_0_y_xzzz_yz[k] * ab_x + g_0_y_xzzz_xyz[k];

                g_0_y_xxzzz_zz[k] = -g_0_y_xzzz_zz[k] * ab_x + g_0_y_xzzz_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyy_xx = cbuffer.data(hd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xyyyy_xy = cbuffer.data(hd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xyyyy_xz = cbuffer.data(hd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xyyyy_yy = cbuffer.data(hd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xyyyy_yz = cbuffer.data(hd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xyyyy_zz = cbuffer.data(hd_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyy_xx, g_0_y_xyyyy_xy, g_0_y_xyyyy_xz, g_0_y_xyyyy_yy, g_0_y_xyyyy_yz, g_0_y_xyyyy_zz, g_0_y_yyyy_xx, g_0_y_yyyy_xxx, g_0_y_yyyy_xxy, g_0_y_yyyy_xxz, g_0_y_yyyy_xy, g_0_y_yyyy_xyy, g_0_y_yyyy_xyz, g_0_y_yyyy_xz, g_0_y_yyyy_xzz, g_0_y_yyyy_yy, g_0_y_yyyy_yz, g_0_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyy_xx[k] = -g_0_y_yyyy_xx[k] * ab_x + g_0_y_yyyy_xxx[k];

                g_0_y_xyyyy_xy[k] = -g_0_y_yyyy_xy[k] * ab_x + g_0_y_yyyy_xxy[k];

                g_0_y_xyyyy_xz[k] = -g_0_y_yyyy_xz[k] * ab_x + g_0_y_yyyy_xxz[k];

                g_0_y_xyyyy_yy[k] = -g_0_y_yyyy_yy[k] * ab_x + g_0_y_yyyy_xyy[k];

                g_0_y_xyyyy_yz[k] = -g_0_y_yyyy_yz[k] * ab_x + g_0_y_yyyy_xyz[k];

                g_0_y_xyyyy_zz[k] = -g_0_y_yyyy_zz[k] * ab_x + g_0_y_yyyy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyz_xx = cbuffer.data(hd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xyyyz_xy = cbuffer.data(hd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xyyyz_xz = cbuffer.data(hd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xyyyz_yy = cbuffer.data(hd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyyyz_yz = cbuffer.data(hd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyyyz_zz = cbuffer.data(hd_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyz_xx, g_0_y_xyyyz_xy, g_0_y_xyyyz_xz, g_0_y_xyyyz_yy, g_0_y_xyyyz_yz, g_0_y_xyyyz_zz, g_0_y_yyyz_xx, g_0_y_yyyz_xxx, g_0_y_yyyz_xxy, g_0_y_yyyz_xxz, g_0_y_yyyz_xy, g_0_y_yyyz_xyy, g_0_y_yyyz_xyz, g_0_y_yyyz_xz, g_0_y_yyyz_xzz, g_0_y_yyyz_yy, g_0_y_yyyz_yz, g_0_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyz_xx[k] = -g_0_y_yyyz_xx[k] * ab_x + g_0_y_yyyz_xxx[k];

                g_0_y_xyyyz_xy[k] = -g_0_y_yyyz_xy[k] * ab_x + g_0_y_yyyz_xxy[k];

                g_0_y_xyyyz_xz[k] = -g_0_y_yyyz_xz[k] * ab_x + g_0_y_yyyz_xxz[k];

                g_0_y_xyyyz_yy[k] = -g_0_y_yyyz_yy[k] * ab_x + g_0_y_yyyz_xyy[k];

                g_0_y_xyyyz_yz[k] = -g_0_y_yyyz_yz[k] * ab_x + g_0_y_yyyz_xyz[k];

                g_0_y_xyyyz_zz[k] = -g_0_y_yyyz_zz[k] * ab_x + g_0_y_yyyz_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzz_xx = cbuffer.data(hd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyyzz_xy = cbuffer.data(hd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyyzz_xz = cbuffer.data(hd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyyzz_yy = cbuffer.data(hd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyyzz_yz = cbuffer.data(hd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyyzz_zz = cbuffer.data(hd_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzz_xx, g_0_y_xyyzz_xy, g_0_y_xyyzz_xz, g_0_y_xyyzz_yy, g_0_y_xyyzz_yz, g_0_y_xyyzz_zz, g_0_y_yyzz_xx, g_0_y_yyzz_xxx, g_0_y_yyzz_xxy, g_0_y_yyzz_xxz, g_0_y_yyzz_xy, g_0_y_yyzz_xyy, g_0_y_yyzz_xyz, g_0_y_yyzz_xz, g_0_y_yyzz_xzz, g_0_y_yyzz_yy, g_0_y_yyzz_yz, g_0_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzz_xx[k] = -g_0_y_yyzz_xx[k] * ab_x + g_0_y_yyzz_xxx[k];

                g_0_y_xyyzz_xy[k] = -g_0_y_yyzz_xy[k] * ab_x + g_0_y_yyzz_xxy[k];

                g_0_y_xyyzz_xz[k] = -g_0_y_yyzz_xz[k] * ab_x + g_0_y_yyzz_xxz[k];

                g_0_y_xyyzz_yy[k] = -g_0_y_yyzz_yy[k] * ab_x + g_0_y_yyzz_xyy[k];

                g_0_y_xyyzz_yz[k] = -g_0_y_yyzz_yz[k] * ab_x + g_0_y_yyzz_xyz[k];

                g_0_y_xyyzz_zz[k] = -g_0_y_yyzz_zz[k] * ab_x + g_0_y_yyzz_xzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzz_xx = cbuffer.data(hd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyzzz_xy = cbuffer.data(hd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyzzz_xz = cbuffer.data(hd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyzzz_yy = cbuffer.data(hd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyzzz_yz = cbuffer.data(hd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyzzz_zz = cbuffer.data(hd_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzz_xx, g_0_y_xyzzz_xy, g_0_y_xyzzz_xz, g_0_y_xyzzz_yy, g_0_y_xyzzz_yz, g_0_y_xyzzz_zz, g_0_y_yzzz_xx, g_0_y_yzzz_xxx, g_0_y_yzzz_xxy, g_0_y_yzzz_xxz, g_0_y_yzzz_xy, g_0_y_yzzz_xyy, g_0_y_yzzz_xyz, g_0_y_yzzz_xz, g_0_y_yzzz_xzz, g_0_y_yzzz_yy, g_0_y_yzzz_yz, g_0_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzz_xx[k] = -g_0_y_yzzz_xx[k] * ab_x + g_0_y_yzzz_xxx[k];

                g_0_y_xyzzz_xy[k] = -g_0_y_yzzz_xy[k] * ab_x + g_0_y_yzzz_xxy[k];

                g_0_y_xyzzz_xz[k] = -g_0_y_yzzz_xz[k] * ab_x + g_0_y_yzzz_xxz[k];

                g_0_y_xyzzz_yy[k] = -g_0_y_yzzz_yy[k] * ab_x + g_0_y_yzzz_xyy[k];

                g_0_y_xyzzz_yz[k] = -g_0_y_yzzz_yz[k] * ab_x + g_0_y_yzzz_xyz[k];

                g_0_y_xyzzz_zz[k] = -g_0_y_yzzz_zz[k] * ab_x + g_0_y_yzzz_xzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzz_xx = cbuffer.data(hd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xzzzz_xy = cbuffer.data(hd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xzzzz_xz = cbuffer.data(hd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xzzzz_yy = cbuffer.data(hd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xzzzz_yz = cbuffer.data(hd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xzzzz_zz = cbuffer.data(hd_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzz_xx, g_0_y_xzzzz_xy, g_0_y_xzzzz_xz, g_0_y_xzzzz_yy, g_0_y_xzzzz_yz, g_0_y_xzzzz_zz, g_0_y_zzzz_xx, g_0_y_zzzz_xxx, g_0_y_zzzz_xxy, g_0_y_zzzz_xxz, g_0_y_zzzz_xy, g_0_y_zzzz_xyy, g_0_y_zzzz_xyz, g_0_y_zzzz_xz, g_0_y_zzzz_xzz, g_0_y_zzzz_yy, g_0_y_zzzz_yz, g_0_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzz_xx[k] = -g_0_y_zzzz_xx[k] * ab_x + g_0_y_zzzz_xxx[k];

                g_0_y_xzzzz_xy[k] = -g_0_y_zzzz_xy[k] * ab_x + g_0_y_zzzz_xxy[k];

                g_0_y_xzzzz_xz[k] = -g_0_y_zzzz_xz[k] * ab_x + g_0_y_zzzz_xxz[k];

                g_0_y_xzzzz_yy[k] = -g_0_y_zzzz_yy[k] * ab_x + g_0_y_zzzz_xyy[k];

                g_0_y_xzzzz_yz[k] = -g_0_y_zzzz_yz[k] * ab_x + g_0_y_zzzz_xyz[k];

                g_0_y_xzzzz_zz[k] = -g_0_y_zzzz_zz[k] * ab_x + g_0_y_zzzz_xzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyy_xx = cbuffer.data(hd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_yyyyy_xy = cbuffer.data(hd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_yyyyy_xz = cbuffer.data(hd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_yyyyy_yy = cbuffer.data(hd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_yyyyy_yz = cbuffer.data(hd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_yyyyy_zz = cbuffer.data(hd_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xx, g_0_y_yyyy_xxy, g_0_y_yyyy_xy, g_0_y_yyyy_xyy, g_0_y_yyyy_xyz, g_0_y_yyyy_xz, g_0_y_yyyy_yy, g_0_y_yyyy_yyy, g_0_y_yyyy_yyz, g_0_y_yyyy_yz, g_0_y_yyyy_yzz, g_0_y_yyyy_zz, g_0_y_yyyyy_xx, g_0_y_yyyyy_xy, g_0_y_yyyyy_xz, g_0_y_yyyyy_yy, g_0_y_yyyyy_yz, g_0_y_yyyyy_zz, g_yyyy_xx, g_yyyy_xy, g_yyyy_xz, g_yyyy_yy, g_yyyy_yz, g_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyy_xx[k] = g_yyyy_xx[k] - g_0_y_yyyy_xx[k] * ab_y + g_0_y_yyyy_xxy[k];

                g_0_y_yyyyy_xy[k] = g_yyyy_xy[k] - g_0_y_yyyy_xy[k] * ab_y + g_0_y_yyyy_xyy[k];

                g_0_y_yyyyy_xz[k] = g_yyyy_xz[k] - g_0_y_yyyy_xz[k] * ab_y + g_0_y_yyyy_xyz[k];

                g_0_y_yyyyy_yy[k] = g_yyyy_yy[k] - g_0_y_yyyy_yy[k] * ab_y + g_0_y_yyyy_yyy[k];

                g_0_y_yyyyy_yz[k] = g_yyyy_yz[k] - g_0_y_yyyy_yz[k] * ab_y + g_0_y_yyyy_yyz[k];

                g_0_y_yyyyy_zz[k] = g_yyyy_zz[k] - g_0_y_yyyy_zz[k] * ab_y + g_0_y_yyyy_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyz_xx = cbuffer.data(hd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_yyyyz_xy = cbuffer.data(hd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_yyyyz_xz = cbuffer.data(hd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_yyyyz_yy = cbuffer.data(hd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_yyyyz_yz = cbuffer.data(hd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_yyyyz_zz = cbuffer.data(hd_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xx, g_0_y_yyyy_xxz, g_0_y_yyyy_xy, g_0_y_yyyy_xyz, g_0_y_yyyy_xz, g_0_y_yyyy_xzz, g_0_y_yyyy_yy, g_0_y_yyyy_yyz, g_0_y_yyyy_yz, g_0_y_yyyy_yzz, g_0_y_yyyy_zz, g_0_y_yyyy_zzz, g_0_y_yyyyz_xx, g_0_y_yyyyz_xy, g_0_y_yyyyz_xz, g_0_y_yyyyz_yy, g_0_y_yyyyz_yz, g_0_y_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyz_xx[k] = -g_0_y_yyyy_xx[k] * ab_z + g_0_y_yyyy_xxz[k];

                g_0_y_yyyyz_xy[k] = -g_0_y_yyyy_xy[k] * ab_z + g_0_y_yyyy_xyz[k];

                g_0_y_yyyyz_xz[k] = -g_0_y_yyyy_xz[k] * ab_z + g_0_y_yyyy_xzz[k];

                g_0_y_yyyyz_yy[k] = -g_0_y_yyyy_yy[k] * ab_z + g_0_y_yyyy_yyz[k];

                g_0_y_yyyyz_yz[k] = -g_0_y_yyyy_yz[k] * ab_z + g_0_y_yyyy_yzz[k];

                g_0_y_yyyyz_zz[k] = -g_0_y_yyyy_zz[k] * ab_z + g_0_y_yyyy_zzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzz_xx = cbuffer.data(hd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_yyyzz_xy = cbuffer.data(hd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_yyyzz_xz = cbuffer.data(hd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_yyyzz_yy = cbuffer.data(hd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_yyyzz_yz = cbuffer.data(hd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_yyyzz_zz = cbuffer.data(hd_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyz_xx, g_0_y_yyyz_xxz, g_0_y_yyyz_xy, g_0_y_yyyz_xyz, g_0_y_yyyz_xz, g_0_y_yyyz_xzz, g_0_y_yyyz_yy, g_0_y_yyyz_yyz, g_0_y_yyyz_yz, g_0_y_yyyz_yzz, g_0_y_yyyz_zz, g_0_y_yyyz_zzz, g_0_y_yyyzz_xx, g_0_y_yyyzz_xy, g_0_y_yyyzz_xz, g_0_y_yyyzz_yy, g_0_y_yyyzz_yz, g_0_y_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzz_xx[k] = -g_0_y_yyyz_xx[k] * ab_z + g_0_y_yyyz_xxz[k];

                g_0_y_yyyzz_xy[k] = -g_0_y_yyyz_xy[k] * ab_z + g_0_y_yyyz_xyz[k];

                g_0_y_yyyzz_xz[k] = -g_0_y_yyyz_xz[k] * ab_z + g_0_y_yyyz_xzz[k];

                g_0_y_yyyzz_yy[k] = -g_0_y_yyyz_yy[k] * ab_z + g_0_y_yyyz_yyz[k];

                g_0_y_yyyzz_yz[k] = -g_0_y_yyyz_yz[k] * ab_z + g_0_y_yyyz_yzz[k];

                g_0_y_yyyzz_zz[k] = -g_0_y_yyyz_zz[k] * ab_z + g_0_y_yyyz_zzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzz_xx = cbuffer.data(hd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_yyzzz_xy = cbuffer.data(hd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_yyzzz_xz = cbuffer.data(hd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_yyzzz_yy = cbuffer.data(hd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_yyzzz_yz = cbuffer.data(hd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_yyzzz_zz = cbuffer.data(hd_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzz_xx, g_0_y_yyzz_xxz, g_0_y_yyzz_xy, g_0_y_yyzz_xyz, g_0_y_yyzz_xz, g_0_y_yyzz_xzz, g_0_y_yyzz_yy, g_0_y_yyzz_yyz, g_0_y_yyzz_yz, g_0_y_yyzz_yzz, g_0_y_yyzz_zz, g_0_y_yyzz_zzz, g_0_y_yyzzz_xx, g_0_y_yyzzz_xy, g_0_y_yyzzz_xz, g_0_y_yyzzz_yy, g_0_y_yyzzz_yz, g_0_y_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzz_xx[k] = -g_0_y_yyzz_xx[k] * ab_z + g_0_y_yyzz_xxz[k];

                g_0_y_yyzzz_xy[k] = -g_0_y_yyzz_xy[k] * ab_z + g_0_y_yyzz_xyz[k];

                g_0_y_yyzzz_xz[k] = -g_0_y_yyzz_xz[k] * ab_z + g_0_y_yyzz_xzz[k];

                g_0_y_yyzzz_yy[k] = -g_0_y_yyzz_yy[k] * ab_z + g_0_y_yyzz_yyz[k];

                g_0_y_yyzzz_yz[k] = -g_0_y_yyzz_yz[k] * ab_z + g_0_y_yyzz_yzz[k];

                g_0_y_yyzzz_zz[k] = -g_0_y_yyzz_zz[k] * ab_z + g_0_y_yyzz_zzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzz_xx = cbuffer.data(hd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yzzzz_xy = cbuffer.data(hd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yzzzz_xz = cbuffer.data(hd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yzzzz_yy = cbuffer.data(hd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yzzzz_yz = cbuffer.data(hd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yzzzz_zz = cbuffer.data(hd_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzz_xx, g_0_y_yzzz_xxz, g_0_y_yzzz_xy, g_0_y_yzzz_xyz, g_0_y_yzzz_xz, g_0_y_yzzz_xzz, g_0_y_yzzz_yy, g_0_y_yzzz_yyz, g_0_y_yzzz_yz, g_0_y_yzzz_yzz, g_0_y_yzzz_zz, g_0_y_yzzz_zzz, g_0_y_yzzzz_xx, g_0_y_yzzzz_xy, g_0_y_yzzzz_xz, g_0_y_yzzzz_yy, g_0_y_yzzzz_yz, g_0_y_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzz_xx[k] = -g_0_y_yzzz_xx[k] * ab_z + g_0_y_yzzz_xxz[k];

                g_0_y_yzzzz_xy[k] = -g_0_y_yzzz_xy[k] * ab_z + g_0_y_yzzz_xyz[k];

                g_0_y_yzzzz_xz[k] = -g_0_y_yzzz_xz[k] * ab_z + g_0_y_yzzz_xzz[k];

                g_0_y_yzzzz_yy[k] = -g_0_y_yzzz_yy[k] * ab_z + g_0_y_yzzz_yyz[k];

                g_0_y_yzzzz_yz[k] = -g_0_y_yzzz_yz[k] * ab_z + g_0_y_yzzz_yzz[k];

                g_0_y_yzzzz_zz[k] = -g_0_y_yzzz_zz[k] * ab_z + g_0_y_yzzz_zzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzz_xx = cbuffer.data(hd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_zzzzz_xy = cbuffer.data(hd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_zzzzz_xz = cbuffer.data(hd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_zzzzz_yy = cbuffer.data(hd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_zzzzz_yz = cbuffer.data(hd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_zzzzz_zz = cbuffer.data(hd_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_xx, g_0_y_zzzz_xxz, g_0_y_zzzz_xy, g_0_y_zzzz_xyz, g_0_y_zzzz_xz, g_0_y_zzzz_xzz, g_0_y_zzzz_yy, g_0_y_zzzz_yyz, g_0_y_zzzz_yz, g_0_y_zzzz_yzz, g_0_y_zzzz_zz, g_0_y_zzzz_zzz, g_0_y_zzzzz_xx, g_0_y_zzzzz_xy, g_0_y_zzzzz_xz, g_0_y_zzzzz_yy, g_0_y_zzzzz_yz, g_0_y_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzz_xx[k] = -g_0_y_zzzz_xx[k] * ab_z + g_0_y_zzzz_xxz[k];

                g_0_y_zzzzz_xy[k] = -g_0_y_zzzz_xy[k] * ab_z + g_0_y_zzzz_xyz[k];

                g_0_y_zzzzz_xz[k] = -g_0_y_zzzz_xz[k] * ab_z + g_0_y_zzzz_xzz[k];

                g_0_y_zzzzz_yy[k] = -g_0_y_zzzz_yy[k] * ab_z + g_0_y_zzzz_yyz[k];

                g_0_y_zzzzz_yz[k] = -g_0_y_zzzz_yz[k] * ab_z + g_0_y_zzzz_yzz[k];

                g_0_y_zzzzz_zz[k] = -g_0_y_zzzz_zz[k] * ab_z + g_0_y_zzzz_zzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxx_xx = cbuffer.data(hd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xxxxx_xy = cbuffer.data(hd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xxxxx_xz = cbuffer.data(hd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xxxxx_yy = cbuffer.data(hd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xxxxx_yz = cbuffer.data(hd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xxxxx_zz = cbuffer.data(hd_geom_01_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_xx, g_0_z_xxxx_xxx, g_0_z_xxxx_xxy, g_0_z_xxxx_xxz, g_0_z_xxxx_xy, g_0_z_xxxx_xyy, g_0_z_xxxx_xyz, g_0_z_xxxx_xz, g_0_z_xxxx_xzz, g_0_z_xxxx_yy, g_0_z_xxxx_yz, g_0_z_xxxx_zz, g_0_z_xxxxx_xx, g_0_z_xxxxx_xy, g_0_z_xxxxx_xz, g_0_z_xxxxx_yy, g_0_z_xxxxx_yz, g_0_z_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxx_xx[k] = -g_0_z_xxxx_xx[k] * ab_x + g_0_z_xxxx_xxx[k];

                g_0_z_xxxxx_xy[k] = -g_0_z_xxxx_xy[k] * ab_x + g_0_z_xxxx_xxy[k];

                g_0_z_xxxxx_xz[k] = -g_0_z_xxxx_xz[k] * ab_x + g_0_z_xxxx_xxz[k];

                g_0_z_xxxxx_yy[k] = -g_0_z_xxxx_yy[k] * ab_x + g_0_z_xxxx_xyy[k];

                g_0_z_xxxxx_yz[k] = -g_0_z_xxxx_yz[k] * ab_x + g_0_z_xxxx_xyz[k];

                g_0_z_xxxxx_zz[k] = -g_0_z_xxxx_zz[k] * ab_x + g_0_z_xxxx_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxy_xx = cbuffer.data(hd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xxxxy_xy = cbuffer.data(hd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xxxxy_xz = cbuffer.data(hd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_xxxxy_yy = cbuffer.data(hd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xxxxy_yz = cbuffer.data(hd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xxxxy_zz = cbuffer.data(hd_geom_01_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxy_xx, g_0_z_xxxxy_xy, g_0_z_xxxxy_xz, g_0_z_xxxxy_yy, g_0_z_xxxxy_yz, g_0_z_xxxxy_zz, g_0_z_xxxy_xx, g_0_z_xxxy_xxx, g_0_z_xxxy_xxy, g_0_z_xxxy_xxz, g_0_z_xxxy_xy, g_0_z_xxxy_xyy, g_0_z_xxxy_xyz, g_0_z_xxxy_xz, g_0_z_xxxy_xzz, g_0_z_xxxy_yy, g_0_z_xxxy_yz, g_0_z_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxy_xx[k] = -g_0_z_xxxy_xx[k] * ab_x + g_0_z_xxxy_xxx[k];

                g_0_z_xxxxy_xy[k] = -g_0_z_xxxy_xy[k] * ab_x + g_0_z_xxxy_xxy[k];

                g_0_z_xxxxy_xz[k] = -g_0_z_xxxy_xz[k] * ab_x + g_0_z_xxxy_xxz[k];

                g_0_z_xxxxy_yy[k] = -g_0_z_xxxy_yy[k] * ab_x + g_0_z_xxxy_xyy[k];

                g_0_z_xxxxy_yz[k] = -g_0_z_xxxy_yz[k] * ab_x + g_0_z_xxxy_xyz[k];

                g_0_z_xxxxy_zz[k] = -g_0_z_xxxy_zz[k] * ab_x + g_0_z_xxxy_xzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxz_xx = cbuffer.data(hd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xxxxz_xy = cbuffer.data(hd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xxxxz_xz = cbuffer.data(hd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_xxxxz_yy = cbuffer.data(hd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xxxxz_yz = cbuffer.data(hd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xxxxz_zz = cbuffer.data(hd_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxz_xx, g_0_z_xxxxz_xy, g_0_z_xxxxz_xz, g_0_z_xxxxz_yy, g_0_z_xxxxz_yz, g_0_z_xxxxz_zz, g_0_z_xxxz_xx, g_0_z_xxxz_xxx, g_0_z_xxxz_xxy, g_0_z_xxxz_xxz, g_0_z_xxxz_xy, g_0_z_xxxz_xyy, g_0_z_xxxz_xyz, g_0_z_xxxz_xz, g_0_z_xxxz_xzz, g_0_z_xxxz_yy, g_0_z_xxxz_yz, g_0_z_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxz_xx[k] = -g_0_z_xxxz_xx[k] * ab_x + g_0_z_xxxz_xxx[k];

                g_0_z_xxxxz_xy[k] = -g_0_z_xxxz_xy[k] * ab_x + g_0_z_xxxz_xxy[k];

                g_0_z_xxxxz_xz[k] = -g_0_z_xxxz_xz[k] * ab_x + g_0_z_xxxz_xxz[k];

                g_0_z_xxxxz_yy[k] = -g_0_z_xxxz_yy[k] * ab_x + g_0_z_xxxz_xyy[k];

                g_0_z_xxxxz_yz[k] = -g_0_z_xxxz_yz[k] * ab_x + g_0_z_xxxz_xyz[k];

                g_0_z_xxxxz_zz[k] = -g_0_z_xxxz_zz[k] * ab_x + g_0_z_xxxz_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyy_xx = cbuffer.data(hd_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xxxyy_xy = cbuffer.data(hd_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xxxyy_xz = cbuffer.data(hd_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_xxxyy_yy = cbuffer.data(hd_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xxxyy_yz = cbuffer.data(hd_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xxxyy_zz = cbuffer.data(hd_geom_01_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyy_xx, g_0_z_xxxyy_xy, g_0_z_xxxyy_xz, g_0_z_xxxyy_yy, g_0_z_xxxyy_yz, g_0_z_xxxyy_zz, g_0_z_xxyy_xx, g_0_z_xxyy_xxx, g_0_z_xxyy_xxy, g_0_z_xxyy_xxz, g_0_z_xxyy_xy, g_0_z_xxyy_xyy, g_0_z_xxyy_xyz, g_0_z_xxyy_xz, g_0_z_xxyy_xzz, g_0_z_xxyy_yy, g_0_z_xxyy_yz, g_0_z_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyy_xx[k] = -g_0_z_xxyy_xx[k] * ab_x + g_0_z_xxyy_xxx[k];

                g_0_z_xxxyy_xy[k] = -g_0_z_xxyy_xy[k] * ab_x + g_0_z_xxyy_xxy[k];

                g_0_z_xxxyy_xz[k] = -g_0_z_xxyy_xz[k] * ab_x + g_0_z_xxyy_xxz[k];

                g_0_z_xxxyy_yy[k] = -g_0_z_xxyy_yy[k] * ab_x + g_0_z_xxyy_xyy[k];

                g_0_z_xxxyy_yz[k] = -g_0_z_xxyy_yz[k] * ab_x + g_0_z_xxyy_xyz[k];

                g_0_z_xxxyy_zz[k] = -g_0_z_xxyy_zz[k] * ab_x + g_0_z_xxyy_xzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyz_xx = cbuffer.data(hd_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xxxyz_xy = cbuffer.data(hd_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xxxyz_xz = cbuffer.data(hd_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_xxxyz_yy = cbuffer.data(hd_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xxxyz_yz = cbuffer.data(hd_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xxxyz_zz = cbuffer.data(hd_geom_01_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyz_xx, g_0_z_xxxyz_xy, g_0_z_xxxyz_xz, g_0_z_xxxyz_yy, g_0_z_xxxyz_yz, g_0_z_xxxyz_zz, g_0_z_xxyz_xx, g_0_z_xxyz_xxx, g_0_z_xxyz_xxy, g_0_z_xxyz_xxz, g_0_z_xxyz_xy, g_0_z_xxyz_xyy, g_0_z_xxyz_xyz, g_0_z_xxyz_xz, g_0_z_xxyz_xzz, g_0_z_xxyz_yy, g_0_z_xxyz_yz, g_0_z_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyz_xx[k] = -g_0_z_xxyz_xx[k] * ab_x + g_0_z_xxyz_xxx[k];

                g_0_z_xxxyz_xy[k] = -g_0_z_xxyz_xy[k] * ab_x + g_0_z_xxyz_xxy[k];

                g_0_z_xxxyz_xz[k] = -g_0_z_xxyz_xz[k] * ab_x + g_0_z_xxyz_xxz[k];

                g_0_z_xxxyz_yy[k] = -g_0_z_xxyz_yy[k] * ab_x + g_0_z_xxyz_xyy[k];

                g_0_z_xxxyz_yz[k] = -g_0_z_xxyz_yz[k] * ab_x + g_0_z_xxyz_xyz[k];

                g_0_z_xxxyz_zz[k] = -g_0_z_xxyz_zz[k] * ab_x + g_0_z_xxyz_xzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzz_xx = cbuffer.data(hd_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xxxzz_xy = cbuffer.data(hd_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xxxzz_xz = cbuffer.data(hd_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_xxxzz_yy = cbuffer.data(hd_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xxxzz_yz = cbuffer.data(hd_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xxxzz_zz = cbuffer.data(hd_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzz_xx, g_0_z_xxxzz_xy, g_0_z_xxxzz_xz, g_0_z_xxxzz_yy, g_0_z_xxxzz_yz, g_0_z_xxxzz_zz, g_0_z_xxzz_xx, g_0_z_xxzz_xxx, g_0_z_xxzz_xxy, g_0_z_xxzz_xxz, g_0_z_xxzz_xy, g_0_z_xxzz_xyy, g_0_z_xxzz_xyz, g_0_z_xxzz_xz, g_0_z_xxzz_xzz, g_0_z_xxzz_yy, g_0_z_xxzz_yz, g_0_z_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzz_xx[k] = -g_0_z_xxzz_xx[k] * ab_x + g_0_z_xxzz_xxx[k];

                g_0_z_xxxzz_xy[k] = -g_0_z_xxzz_xy[k] * ab_x + g_0_z_xxzz_xxy[k];

                g_0_z_xxxzz_xz[k] = -g_0_z_xxzz_xz[k] * ab_x + g_0_z_xxzz_xxz[k];

                g_0_z_xxxzz_yy[k] = -g_0_z_xxzz_yy[k] * ab_x + g_0_z_xxzz_xyy[k];

                g_0_z_xxxzz_yz[k] = -g_0_z_xxzz_yz[k] * ab_x + g_0_z_xxzz_xyz[k];

                g_0_z_xxxzz_zz[k] = -g_0_z_xxzz_zz[k] * ab_x + g_0_z_xxzz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyy_xx = cbuffer.data(hd_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xxyyy_xy = cbuffer.data(hd_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xxyyy_xz = cbuffer.data(hd_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_xxyyy_yy = cbuffer.data(hd_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xxyyy_yz = cbuffer.data(hd_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xxyyy_zz = cbuffer.data(hd_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyy_xx, g_0_z_xxyyy_xy, g_0_z_xxyyy_xz, g_0_z_xxyyy_yy, g_0_z_xxyyy_yz, g_0_z_xxyyy_zz, g_0_z_xyyy_xx, g_0_z_xyyy_xxx, g_0_z_xyyy_xxy, g_0_z_xyyy_xxz, g_0_z_xyyy_xy, g_0_z_xyyy_xyy, g_0_z_xyyy_xyz, g_0_z_xyyy_xz, g_0_z_xyyy_xzz, g_0_z_xyyy_yy, g_0_z_xyyy_yz, g_0_z_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyy_xx[k] = -g_0_z_xyyy_xx[k] * ab_x + g_0_z_xyyy_xxx[k];

                g_0_z_xxyyy_xy[k] = -g_0_z_xyyy_xy[k] * ab_x + g_0_z_xyyy_xxy[k];

                g_0_z_xxyyy_xz[k] = -g_0_z_xyyy_xz[k] * ab_x + g_0_z_xyyy_xxz[k];

                g_0_z_xxyyy_yy[k] = -g_0_z_xyyy_yy[k] * ab_x + g_0_z_xyyy_xyy[k];

                g_0_z_xxyyy_yz[k] = -g_0_z_xyyy_yz[k] * ab_x + g_0_z_xyyy_xyz[k];

                g_0_z_xxyyy_zz[k] = -g_0_z_xyyy_zz[k] * ab_x + g_0_z_xyyy_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyz_xx = cbuffer.data(hd_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xxyyz_xy = cbuffer.data(hd_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xxyyz_xz = cbuffer.data(hd_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_xxyyz_yy = cbuffer.data(hd_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xxyyz_yz = cbuffer.data(hd_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xxyyz_zz = cbuffer.data(hd_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyz_xx, g_0_z_xxyyz_xy, g_0_z_xxyyz_xz, g_0_z_xxyyz_yy, g_0_z_xxyyz_yz, g_0_z_xxyyz_zz, g_0_z_xyyz_xx, g_0_z_xyyz_xxx, g_0_z_xyyz_xxy, g_0_z_xyyz_xxz, g_0_z_xyyz_xy, g_0_z_xyyz_xyy, g_0_z_xyyz_xyz, g_0_z_xyyz_xz, g_0_z_xyyz_xzz, g_0_z_xyyz_yy, g_0_z_xyyz_yz, g_0_z_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyz_xx[k] = -g_0_z_xyyz_xx[k] * ab_x + g_0_z_xyyz_xxx[k];

                g_0_z_xxyyz_xy[k] = -g_0_z_xyyz_xy[k] * ab_x + g_0_z_xyyz_xxy[k];

                g_0_z_xxyyz_xz[k] = -g_0_z_xyyz_xz[k] * ab_x + g_0_z_xyyz_xxz[k];

                g_0_z_xxyyz_yy[k] = -g_0_z_xyyz_yy[k] * ab_x + g_0_z_xyyz_xyy[k];

                g_0_z_xxyyz_yz[k] = -g_0_z_xyyz_yz[k] * ab_x + g_0_z_xyyz_xyz[k];

                g_0_z_xxyyz_zz[k] = -g_0_z_xyyz_zz[k] * ab_x + g_0_z_xyyz_xzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzz_xx = cbuffer.data(hd_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxyzz_xy = cbuffer.data(hd_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxyzz_xz = cbuffer.data(hd_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxyzz_yy = cbuffer.data(hd_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxyzz_yz = cbuffer.data(hd_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxyzz_zz = cbuffer.data(hd_geom_01_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzz_xx, g_0_z_xxyzz_xy, g_0_z_xxyzz_xz, g_0_z_xxyzz_yy, g_0_z_xxyzz_yz, g_0_z_xxyzz_zz, g_0_z_xyzz_xx, g_0_z_xyzz_xxx, g_0_z_xyzz_xxy, g_0_z_xyzz_xxz, g_0_z_xyzz_xy, g_0_z_xyzz_xyy, g_0_z_xyzz_xyz, g_0_z_xyzz_xz, g_0_z_xyzz_xzz, g_0_z_xyzz_yy, g_0_z_xyzz_yz, g_0_z_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzz_xx[k] = -g_0_z_xyzz_xx[k] * ab_x + g_0_z_xyzz_xxx[k];

                g_0_z_xxyzz_xy[k] = -g_0_z_xyzz_xy[k] * ab_x + g_0_z_xyzz_xxy[k];

                g_0_z_xxyzz_xz[k] = -g_0_z_xyzz_xz[k] * ab_x + g_0_z_xyzz_xxz[k];

                g_0_z_xxyzz_yy[k] = -g_0_z_xyzz_yy[k] * ab_x + g_0_z_xyzz_xyy[k];

                g_0_z_xxyzz_yz[k] = -g_0_z_xyzz_yz[k] * ab_x + g_0_z_xyzz_xyz[k];

                g_0_z_xxyzz_zz[k] = -g_0_z_xyzz_zz[k] * ab_x + g_0_z_xyzz_xzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzz_xx = cbuffer.data(hd_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxzzz_xy = cbuffer.data(hd_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxzzz_xz = cbuffer.data(hd_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxzzz_yy = cbuffer.data(hd_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxzzz_yz = cbuffer.data(hd_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxzzz_zz = cbuffer.data(hd_geom_01_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzz_xx, g_0_z_xxzzz_xy, g_0_z_xxzzz_xz, g_0_z_xxzzz_yy, g_0_z_xxzzz_yz, g_0_z_xxzzz_zz, g_0_z_xzzz_xx, g_0_z_xzzz_xxx, g_0_z_xzzz_xxy, g_0_z_xzzz_xxz, g_0_z_xzzz_xy, g_0_z_xzzz_xyy, g_0_z_xzzz_xyz, g_0_z_xzzz_xz, g_0_z_xzzz_xzz, g_0_z_xzzz_yy, g_0_z_xzzz_yz, g_0_z_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzz_xx[k] = -g_0_z_xzzz_xx[k] * ab_x + g_0_z_xzzz_xxx[k];

                g_0_z_xxzzz_xy[k] = -g_0_z_xzzz_xy[k] * ab_x + g_0_z_xzzz_xxy[k];

                g_0_z_xxzzz_xz[k] = -g_0_z_xzzz_xz[k] * ab_x + g_0_z_xzzz_xxz[k];

                g_0_z_xxzzz_yy[k] = -g_0_z_xzzz_yy[k] * ab_x + g_0_z_xzzz_xyy[k];

                g_0_z_xxzzz_yz[k] = -g_0_z_xzzz_yz[k] * ab_x + g_0_z_xzzz_xyz[k];

                g_0_z_xxzzz_zz[k] = -g_0_z_xzzz_zz[k] * ab_x + g_0_z_xzzz_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyy_xx = cbuffer.data(hd_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xyyyy_xy = cbuffer.data(hd_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xyyyy_xz = cbuffer.data(hd_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xyyyy_yy = cbuffer.data(hd_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xyyyy_yz = cbuffer.data(hd_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xyyyy_zz = cbuffer.data(hd_geom_01_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyy_xx, g_0_z_xyyyy_xy, g_0_z_xyyyy_xz, g_0_z_xyyyy_yy, g_0_z_xyyyy_yz, g_0_z_xyyyy_zz, g_0_z_yyyy_xx, g_0_z_yyyy_xxx, g_0_z_yyyy_xxy, g_0_z_yyyy_xxz, g_0_z_yyyy_xy, g_0_z_yyyy_xyy, g_0_z_yyyy_xyz, g_0_z_yyyy_xz, g_0_z_yyyy_xzz, g_0_z_yyyy_yy, g_0_z_yyyy_yz, g_0_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyy_xx[k] = -g_0_z_yyyy_xx[k] * ab_x + g_0_z_yyyy_xxx[k];

                g_0_z_xyyyy_xy[k] = -g_0_z_yyyy_xy[k] * ab_x + g_0_z_yyyy_xxy[k];

                g_0_z_xyyyy_xz[k] = -g_0_z_yyyy_xz[k] * ab_x + g_0_z_yyyy_xxz[k];

                g_0_z_xyyyy_yy[k] = -g_0_z_yyyy_yy[k] * ab_x + g_0_z_yyyy_xyy[k];

                g_0_z_xyyyy_yz[k] = -g_0_z_yyyy_yz[k] * ab_x + g_0_z_yyyy_xyz[k];

                g_0_z_xyyyy_zz[k] = -g_0_z_yyyy_zz[k] * ab_x + g_0_z_yyyy_xzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyz_xx = cbuffer.data(hd_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xyyyz_xy = cbuffer.data(hd_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xyyyz_xz = cbuffer.data(hd_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xyyyz_yy = cbuffer.data(hd_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xyyyz_yz = cbuffer.data(hd_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xyyyz_zz = cbuffer.data(hd_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyz_xx, g_0_z_xyyyz_xy, g_0_z_xyyyz_xz, g_0_z_xyyyz_yy, g_0_z_xyyyz_yz, g_0_z_xyyyz_zz, g_0_z_yyyz_xx, g_0_z_yyyz_xxx, g_0_z_yyyz_xxy, g_0_z_yyyz_xxz, g_0_z_yyyz_xy, g_0_z_yyyz_xyy, g_0_z_yyyz_xyz, g_0_z_yyyz_xz, g_0_z_yyyz_xzz, g_0_z_yyyz_yy, g_0_z_yyyz_yz, g_0_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyz_xx[k] = -g_0_z_yyyz_xx[k] * ab_x + g_0_z_yyyz_xxx[k];

                g_0_z_xyyyz_xy[k] = -g_0_z_yyyz_xy[k] * ab_x + g_0_z_yyyz_xxy[k];

                g_0_z_xyyyz_xz[k] = -g_0_z_yyyz_xz[k] * ab_x + g_0_z_yyyz_xxz[k];

                g_0_z_xyyyz_yy[k] = -g_0_z_yyyz_yy[k] * ab_x + g_0_z_yyyz_xyy[k];

                g_0_z_xyyyz_yz[k] = -g_0_z_yyyz_yz[k] * ab_x + g_0_z_yyyz_xyz[k];

                g_0_z_xyyyz_zz[k] = -g_0_z_yyyz_zz[k] * ab_x + g_0_z_yyyz_xzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzz_xx = cbuffer.data(hd_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xyyzz_xy = cbuffer.data(hd_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xyyzz_xz = cbuffer.data(hd_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xyyzz_yy = cbuffer.data(hd_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xyyzz_yz = cbuffer.data(hd_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xyyzz_zz = cbuffer.data(hd_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzz_xx, g_0_z_xyyzz_xy, g_0_z_xyyzz_xz, g_0_z_xyyzz_yy, g_0_z_xyyzz_yz, g_0_z_xyyzz_zz, g_0_z_yyzz_xx, g_0_z_yyzz_xxx, g_0_z_yyzz_xxy, g_0_z_yyzz_xxz, g_0_z_yyzz_xy, g_0_z_yyzz_xyy, g_0_z_yyzz_xyz, g_0_z_yyzz_xz, g_0_z_yyzz_xzz, g_0_z_yyzz_yy, g_0_z_yyzz_yz, g_0_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzz_xx[k] = -g_0_z_yyzz_xx[k] * ab_x + g_0_z_yyzz_xxx[k];

                g_0_z_xyyzz_xy[k] = -g_0_z_yyzz_xy[k] * ab_x + g_0_z_yyzz_xxy[k];

                g_0_z_xyyzz_xz[k] = -g_0_z_yyzz_xz[k] * ab_x + g_0_z_yyzz_xxz[k];

                g_0_z_xyyzz_yy[k] = -g_0_z_yyzz_yy[k] * ab_x + g_0_z_yyzz_xyy[k];

                g_0_z_xyyzz_yz[k] = -g_0_z_yyzz_yz[k] * ab_x + g_0_z_yyzz_xyz[k];

                g_0_z_xyyzz_zz[k] = -g_0_z_yyzz_zz[k] * ab_x + g_0_z_yyzz_xzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzz_xx = cbuffer.data(hd_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xyzzz_xy = cbuffer.data(hd_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xyzzz_xz = cbuffer.data(hd_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xyzzz_yy = cbuffer.data(hd_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xyzzz_yz = cbuffer.data(hd_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xyzzz_zz = cbuffer.data(hd_geom_01_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzz_xx, g_0_z_xyzzz_xy, g_0_z_xyzzz_xz, g_0_z_xyzzz_yy, g_0_z_xyzzz_yz, g_0_z_xyzzz_zz, g_0_z_yzzz_xx, g_0_z_yzzz_xxx, g_0_z_yzzz_xxy, g_0_z_yzzz_xxz, g_0_z_yzzz_xy, g_0_z_yzzz_xyy, g_0_z_yzzz_xyz, g_0_z_yzzz_xz, g_0_z_yzzz_xzz, g_0_z_yzzz_yy, g_0_z_yzzz_yz, g_0_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzz_xx[k] = -g_0_z_yzzz_xx[k] * ab_x + g_0_z_yzzz_xxx[k];

                g_0_z_xyzzz_xy[k] = -g_0_z_yzzz_xy[k] * ab_x + g_0_z_yzzz_xxy[k];

                g_0_z_xyzzz_xz[k] = -g_0_z_yzzz_xz[k] * ab_x + g_0_z_yzzz_xxz[k];

                g_0_z_xyzzz_yy[k] = -g_0_z_yzzz_yy[k] * ab_x + g_0_z_yzzz_xyy[k];

                g_0_z_xyzzz_yz[k] = -g_0_z_yzzz_yz[k] * ab_x + g_0_z_yzzz_xyz[k];

                g_0_z_xyzzz_zz[k] = -g_0_z_yzzz_zz[k] * ab_x + g_0_z_yzzz_xzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzz_xx = cbuffer.data(hd_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xzzzz_xy = cbuffer.data(hd_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xzzzz_xz = cbuffer.data(hd_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xzzzz_yy = cbuffer.data(hd_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xzzzz_yz = cbuffer.data(hd_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xzzzz_zz = cbuffer.data(hd_geom_01_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzz_xx, g_0_z_xzzzz_xy, g_0_z_xzzzz_xz, g_0_z_xzzzz_yy, g_0_z_xzzzz_yz, g_0_z_xzzzz_zz, g_0_z_zzzz_xx, g_0_z_zzzz_xxx, g_0_z_zzzz_xxy, g_0_z_zzzz_xxz, g_0_z_zzzz_xy, g_0_z_zzzz_xyy, g_0_z_zzzz_xyz, g_0_z_zzzz_xz, g_0_z_zzzz_xzz, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzz_xx[k] = -g_0_z_zzzz_xx[k] * ab_x + g_0_z_zzzz_xxx[k];

                g_0_z_xzzzz_xy[k] = -g_0_z_zzzz_xy[k] * ab_x + g_0_z_zzzz_xxy[k];

                g_0_z_xzzzz_xz[k] = -g_0_z_zzzz_xz[k] * ab_x + g_0_z_zzzz_xxz[k];

                g_0_z_xzzzz_yy[k] = -g_0_z_zzzz_yy[k] * ab_x + g_0_z_zzzz_xyy[k];

                g_0_z_xzzzz_yz[k] = -g_0_z_zzzz_yz[k] * ab_x + g_0_z_zzzz_xyz[k];

                g_0_z_xzzzz_zz[k] = -g_0_z_zzzz_zz[k] * ab_x + g_0_z_zzzz_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyy_xx = cbuffer.data(hd_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_yyyyy_xy = cbuffer.data(hd_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_yyyyy_xz = cbuffer.data(hd_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_yyyyy_yy = cbuffer.data(hd_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_yyyyy_yz = cbuffer.data(hd_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_yyyyy_zz = cbuffer.data(hd_geom_01_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_xx, g_0_z_yyyy_xxy, g_0_z_yyyy_xy, g_0_z_yyyy_xyy, g_0_z_yyyy_xyz, g_0_z_yyyy_xz, g_0_z_yyyy_yy, g_0_z_yyyy_yyy, g_0_z_yyyy_yyz, g_0_z_yyyy_yz, g_0_z_yyyy_yzz, g_0_z_yyyy_zz, g_0_z_yyyyy_xx, g_0_z_yyyyy_xy, g_0_z_yyyyy_xz, g_0_z_yyyyy_yy, g_0_z_yyyyy_yz, g_0_z_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyy_xx[k] = -g_0_z_yyyy_xx[k] * ab_y + g_0_z_yyyy_xxy[k];

                g_0_z_yyyyy_xy[k] = -g_0_z_yyyy_xy[k] * ab_y + g_0_z_yyyy_xyy[k];

                g_0_z_yyyyy_xz[k] = -g_0_z_yyyy_xz[k] * ab_y + g_0_z_yyyy_xyz[k];

                g_0_z_yyyyy_yy[k] = -g_0_z_yyyy_yy[k] * ab_y + g_0_z_yyyy_yyy[k];

                g_0_z_yyyyy_yz[k] = -g_0_z_yyyy_yz[k] * ab_y + g_0_z_yyyy_yyz[k];

                g_0_z_yyyyy_zz[k] = -g_0_z_yyyy_zz[k] * ab_y + g_0_z_yyyy_yzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyz_xx = cbuffer.data(hd_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_yyyyz_xy = cbuffer.data(hd_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_yyyyz_xz = cbuffer.data(hd_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_yyyyz_yy = cbuffer.data(hd_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_yyyyz_yz = cbuffer.data(hd_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_yyyyz_zz = cbuffer.data(hd_geom_01_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyz_xx, g_0_z_yyyyz_xy, g_0_z_yyyyz_xz, g_0_z_yyyyz_yy, g_0_z_yyyyz_yz, g_0_z_yyyyz_zz, g_0_z_yyyz_xx, g_0_z_yyyz_xxy, g_0_z_yyyz_xy, g_0_z_yyyz_xyy, g_0_z_yyyz_xyz, g_0_z_yyyz_xz, g_0_z_yyyz_yy, g_0_z_yyyz_yyy, g_0_z_yyyz_yyz, g_0_z_yyyz_yz, g_0_z_yyyz_yzz, g_0_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyz_xx[k] = -g_0_z_yyyz_xx[k] * ab_y + g_0_z_yyyz_xxy[k];

                g_0_z_yyyyz_xy[k] = -g_0_z_yyyz_xy[k] * ab_y + g_0_z_yyyz_xyy[k];

                g_0_z_yyyyz_xz[k] = -g_0_z_yyyz_xz[k] * ab_y + g_0_z_yyyz_xyz[k];

                g_0_z_yyyyz_yy[k] = -g_0_z_yyyz_yy[k] * ab_y + g_0_z_yyyz_yyy[k];

                g_0_z_yyyyz_yz[k] = -g_0_z_yyyz_yz[k] * ab_y + g_0_z_yyyz_yyz[k];

                g_0_z_yyyyz_zz[k] = -g_0_z_yyyz_zz[k] * ab_y + g_0_z_yyyz_yzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzz_xx = cbuffer.data(hd_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_yyyzz_xy = cbuffer.data(hd_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_yyyzz_xz = cbuffer.data(hd_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_yyyzz_yy = cbuffer.data(hd_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_yyyzz_yz = cbuffer.data(hd_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_yyyzz_zz = cbuffer.data(hd_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzz_xx, g_0_z_yyyzz_xy, g_0_z_yyyzz_xz, g_0_z_yyyzz_yy, g_0_z_yyyzz_yz, g_0_z_yyyzz_zz, g_0_z_yyzz_xx, g_0_z_yyzz_xxy, g_0_z_yyzz_xy, g_0_z_yyzz_xyy, g_0_z_yyzz_xyz, g_0_z_yyzz_xz, g_0_z_yyzz_yy, g_0_z_yyzz_yyy, g_0_z_yyzz_yyz, g_0_z_yyzz_yz, g_0_z_yyzz_yzz, g_0_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzz_xx[k] = -g_0_z_yyzz_xx[k] * ab_y + g_0_z_yyzz_xxy[k];

                g_0_z_yyyzz_xy[k] = -g_0_z_yyzz_xy[k] * ab_y + g_0_z_yyzz_xyy[k];

                g_0_z_yyyzz_xz[k] = -g_0_z_yyzz_xz[k] * ab_y + g_0_z_yyzz_xyz[k];

                g_0_z_yyyzz_yy[k] = -g_0_z_yyzz_yy[k] * ab_y + g_0_z_yyzz_yyy[k];

                g_0_z_yyyzz_yz[k] = -g_0_z_yyzz_yz[k] * ab_y + g_0_z_yyzz_yyz[k];

                g_0_z_yyyzz_zz[k] = -g_0_z_yyzz_zz[k] * ab_y + g_0_z_yyzz_yzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzz_xx = cbuffer.data(hd_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_yyzzz_xy = cbuffer.data(hd_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_yyzzz_xz = cbuffer.data(hd_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_yyzzz_yy = cbuffer.data(hd_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_yyzzz_yz = cbuffer.data(hd_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_yyzzz_zz = cbuffer.data(hd_geom_01_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzz_xx, g_0_z_yyzzz_xy, g_0_z_yyzzz_xz, g_0_z_yyzzz_yy, g_0_z_yyzzz_yz, g_0_z_yyzzz_zz, g_0_z_yzzz_xx, g_0_z_yzzz_xxy, g_0_z_yzzz_xy, g_0_z_yzzz_xyy, g_0_z_yzzz_xyz, g_0_z_yzzz_xz, g_0_z_yzzz_yy, g_0_z_yzzz_yyy, g_0_z_yzzz_yyz, g_0_z_yzzz_yz, g_0_z_yzzz_yzz, g_0_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzz_xx[k] = -g_0_z_yzzz_xx[k] * ab_y + g_0_z_yzzz_xxy[k];

                g_0_z_yyzzz_xy[k] = -g_0_z_yzzz_xy[k] * ab_y + g_0_z_yzzz_xyy[k];

                g_0_z_yyzzz_xz[k] = -g_0_z_yzzz_xz[k] * ab_y + g_0_z_yzzz_xyz[k];

                g_0_z_yyzzz_yy[k] = -g_0_z_yzzz_yy[k] * ab_y + g_0_z_yzzz_yyy[k];

                g_0_z_yyzzz_yz[k] = -g_0_z_yzzz_yz[k] * ab_y + g_0_z_yzzz_yyz[k];

                g_0_z_yyzzz_zz[k] = -g_0_z_yzzz_zz[k] * ab_y + g_0_z_yzzz_yzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzz_xx = cbuffer.data(hd_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_yzzzz_xy = cbuffer.data(hd_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_yzzzz_xz = cbuffer.data(hd_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_yzzzz_yy = cbuffer.data(hd_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_yzzzz_yz = cbuffer.data(hd_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_yzzzz_zz = cbuffer.data(hd_geom_01_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzz_xx, g_0_z_yzzzz_xy, g_0_z_yzzzz_xz, g_0_z_yzzzz_yy, g_0_z_yzzzz_yz, g_0_z_yzzzz_zz, g_0_z_zzzz_xx, g_0_z_zzzz_xxy, g_0_z_zzzz_xy, g_0_z_zzzz_xyy, g_0_z_zzzz_xyz, g_0_z_zzzz_xz, g_0_z_zzzz_yy, g_0_z_zzzz_yyy, g_0_z_zzzz_yyz, g_0_z_zzzz_yz, g_0_z_zzzz_yzz, g_0_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzz_xx[k] = -g_0_z_zzzz_xx[k] * ab_y + g_0_z_zzzz_xxy[k];

                g_0_z_yzzzz_xy[k] = -g_0_z_zzzz_xy[k] * ab_y + g_0_z_zzzz_xyy[k];

                g_0_z_yzzzz_xz[k] = -g_0_z_zzzz_xz[k] * ab_y + g_0_z_zzzz_xyz[k];

                g_0_z_yzzzz_yy[k] = -g_0_z_zzzz_yy[k] * ab_y + g_0_z_zzzz_yyy[k];

                g_0_z_yzzzz_yz[k] = -g_0_z_zzzz_yz[k] * ab_y + g_0_z_zzzz_yyz[k];

                g_0_z_yzzzz_zz[k] = -g_0_z_zzzz_zz[k] * ab_y + g_0_z_zzzz_yzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzz_xx = cbuffer.data(hd_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_zzzzz_xy = cbuffer.data(hd_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_zzzzz_xz = cbuffer.data(hd_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_zzzzz_yy = cbuffer.data(hd_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_zzzzz_yz = cbuffer.data(hd_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_zzzzz_zz = cbuffer.data(hd_geom_01_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xx, g_0_z_zzzz_xxz, g_0_z_zzzz_xy, g_0_z_zzzz_xyz, g_0_z_zzzz_xz, g_0_z_zzzz_xzz, g_0_z_zzzz_yy, g_0_z_zzzz_yyz, g_0_z_zzzz_yz, g_0_z_zzzz_yzz, g_0_z_zzzz_zz, g_0_z_zzzz_zzz, g_0_z_zzzzz_xx, g_0_z_zzzzz_xy, g_0_z_zzzzz_xz, g_0_z_zzzzz_yy, g_0_z_zzzzz_yz, g_0_z_zzzzz_zz, g_zzzz_xx, g_zzzz_xy, g_zzzz_xz, g_zzzz_yy, g_zzzz_yz, g_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzz_xx[k] = g_zzzz_xx[k] - g_0_z_zzzz_xx[k] * ab_z + g_0_z_zzzz_xxz[k];

                g_0_z_zzzzz_xy[k] = g_zzzz_xy[k] - g_0_z_zzzz_xy[k] * ab_z + g_0_z_zzzz_xyz[k];

                g_0_z_zzzzz_xz[k] = g_zzzz_xz[k] - g_0_z_zzzz_xz[k] * ab_z + g_0_z_zzzz_xzz[k];

                g_0_z_zzzzz_yy[k] = g_zzzz_yy[k] - g_0_z_zzzz_yy[k] * ab_z + g_0_z_zzzz_yyz[k];

                g_0_z_zzzzz_yz[k] = g_zzzz_yz[k] - g_0_z_zzzz_yz[k] * ab_z + g_0_z_zzzz_yzz[k];

                g_0_z_zzzzz_zz[k] = g_zzzz_zz[k] - g_0_z_zzzz_zz[k] * ab_z + g_0_z_zzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

