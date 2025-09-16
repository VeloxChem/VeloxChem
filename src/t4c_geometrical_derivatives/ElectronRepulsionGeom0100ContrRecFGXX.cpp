#include "ElectronRepulsionGeom0100ContrRecFGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_fgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_fgxx,
                                            const size_t idx_dgxx,
                                            const size_t idx_geom_01_dgxx,
                                            const size_t idx_geom_01_dhxx,
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
            /// Set up components of auxilary buffer : DGSS

            const auto dg_off = idx_dgxx + i * dcomps + j;

            auto g_xx_xxxx = cbuffer.data(dg_off + 0 * ccomps * dcomps);

            auto g_xx_xxxy = cbuffer.data(dg_off + 1 * ccomps * dcomps);

            auto g_xx_xxxz = cbuffer.data(dg_off + 2 * ccomps * dcomps);

            auto g_xx_xxyy = cbuffer.data(dg_off + 3 * ccomps * dcomps);

            auto g_xx_xxyz = cbuffer.data(dg_off + 4 * ccomps * dcomps);

            auto g_xx_xxzz = cbuffer.data(dg_off + 5 * ccomps * dcomps);

            auto g_xx_xyyy = cbuffer.data(dg_off + 6 * ccomps * dcomps);

            auto g_xx_xyyz = cbuffer.data(dg_off + 7 * ccomps * dcomps);

            auto g_xx_xyzz = cbuffer.data(dg_off + 8 * ccomps * dcomps);

            auto g_xx_xzzz = cbuffer.data(dg_off + 9 * ccomps * dcomps);

            auto g_xx_yyyy = cbuffer.data(dg_off + 10 * ccomps * dcomps);

            auto g_xx_yyyz = cbuffer.data(dg_off + 11 * ccomps * dcomps);

            auto g_xx_yyzz = cbuffer.data(dg_off + 12 * ccomps * dcomps);

            auto g_xx_yzzz = cbuffer.data(dg_off + 13 * ccomps * dcomps);

            auto g_xx_zzzz = cbuffer.data(dg_off + 14 * ccomps * dcomps);

            auto g_xy_xxxx = cbuffer.data(dg_off + 15 * ccomps * dcomps);

            auto g_xy_xxxy = cbuffer.data(dg_off + 16 * ccomps * dcomps);

            auto g_xy_xxxz = cbuffer.data(dg_off + 17 * ccomps * dcomps);

            auto g_xy_xxyy = cbuffer.data(dg_off + 18 * ccomps * dcomps);

            auto g_xy_xxyz = cbuffer.data(dg_off + 19 * ccomps * dcomps);

            auto g_xy_xxzz = cbuffer.data(dg_off + 20 * ccomps * dcomps);

            auto g_xy_xyyy = cbuffer.data(dg_off + 21 * ccomps * dcomps);

            auto g_xy_xyyz = cbuffer.data(dg_off + 22 * ccomps * dcomps);

            auto g_xy_xyzz = cbuffer.data(dg_off + 23 * ccomps * dcomps);

            auto g_xy_xzzz = cbuffer.data(dg_off + 24 * ccomps * dcomps);

            auto g_xy_yyyy = cbuffer.data(dg_off + 25 * ccomps * dcomps);

            auto g_xy_yyyz = cbuffer.data(dg_off + 26 * ccomps * dcomps);

            auto g_xy_yyzz = cbuffer.data(dg_off + 27 * ccomps * dcomps);

            auto g_xy_yzzz = cbuffer.data(dg_off + 28 * ccomps * dcomps);

            auto g_xy_zzzz = cbuffer.data(dg_off + 29 * ccomps * dcomps);

            auto g_xz_xxxx = cbuffer.data(dg_off + 30 * ccomps * dcomps);

            auto g_xz_xxxy = cbuffer.data(dg_off + 31 * ccomps * dcomps);

            auto g_xz_xxxz = cbuffer.data(dg_off + 32 * ccomps * dcomps);

            auto g_xz_xxyy = cbuffer.data(dg_off + 33 * ccomps * dcomps);

            auto g_xz_xxyz = cbuffer.data(dg_off + 34 * ccomps * dcomps);

            auto g_xz_xxzz = cbuffer.data(dg_off + 35 * ccomps * dcomps);

            auto g_xz_xyyy = cbuffer.data(dg_off + 36 * ccomps * dcomps);

            auto g_xz_xyyz = cbuffer.data(dg_off + 37 * ccomps * dcomps);

            auto g_xz_xyzz = cbuffer.data(dg_off + 38 * ccomps * dcomps);

            auto g_xz_xzzz = cbuffer.data(dg_off + 39 * ccomps * dcomps);

            auto g_xz_yyyy = cbuffer.data(dg_off + 40 * ccomps * dcomps);

            auto g_xz_yyyz = cbuffer.data(dg_off + 41 * ccomps * dcomps);

            auto g_xz_yyzz = cbuffer.data(dg_off + 42 * ccomps * dcomps);

            auto g_xz_yzzz = cbuffer.data(dg_off + 43 * ccomps * dcomps);

            auto g_xz_zzzz = cbuffer.data(dg_off + 44 * ccomps * dcomps);

            auto g_yy_xxxx = cbuffer.data(dg_off + 45 * ccomps * dcomps);

            auto g_yy_xxxy = cbuffer.data(dg_off + 46 * ccomps * dcomps);

            auto g_yy_xxxz = cbuffer.data(dg_off + 47 * ccomps * dcomps);

            auto g_yy_xxyy = cbuffer.data(dg_off + 48 * ccomps * dcomps);

            auto g_yy_xxyz = cbuffer.data(dg_off + 49 * ccomps * dcomps);

            auto g_yy_xxzz = cbuffer.data(dg_off + 50 * ccomps * dcomps);

            auto g_yy_xyyy = cbuffer.data(dg_off + 51 * ccomps * dcomps);

            auto g_yy_xyyz = cbuffer.data(dg_off + 52 * ccomps * dcomps);

            auto g_yy_xyzz = cbuffer.data(dg_off + 53 * ccomps * dcomps);

            auto g_yy_xzzz = cbuffer.data(dg_off + 54 * ccomps * dcomps);

            auto g_yy_yyyy = cbuffer.data(dg_off + 55 * ccomps * dcomps);

            auto g_yy_yyyz = cbuffer.data(dg_off + 56 * ccomps * dcomps);

            auto g_yy_yyzz = cbuffer.data(dg_off + 57 * ccomps * dcomps);

            auto g_yy_yzzz = cbuffer.data(dg_off + 58 * ccomps * dcomps);

            auto g_yy_zzzz = cbuffer.data(dg_off + 59 * ccomps * dcomps);

            auto g_yz_xxxx = cbuffer.data(dg_off + 60 * ccomps * dcomps);

            auto g_yz_xxxy = cbuffer.data(dg_off + 61 * ccomps * dcomps);

            auto g_yz_xxxz = cbuffer.data(dg_off + 62 * ccomps * dcomps);

            auto g_yz_xxyy = cbuffer.data(dg_off + 63 * ccomps * dcomps);

            auto g_yz_xxyz = cbuffer.data(dg_off + 64 * ccomps * dcomps);

            auto g_yz_xxzz = cbuffer.data(dg_off + 65 * ccomps * dcomps);

            auto g_yz_xyyy = cbuffer.data(dg_off + 66 * ccomps * dcomps);

            auto g_yz_xyyz = cbuffer.data(dg_off + 67 * ccomps * dcomps);

            auto g_yz_xyzz = cbuffer.data(dg_off + 68 * ccomps * dcomps);

            auto g_yz_xzzz = cbuffer.data(dg_off + 69 * ccomps * dcomps);

            auto g_yz_yyyy = cbuffer.data(dg_off + 70 * ccomps * dcomps);

            auto g_yz_yyyz = cbuffer.data(dg_off + 71 * ccomps * dcomps);

            auto g_yz_yyzz = cbuffer.data(dg_off + 72 * ccomps * dcomps);

            auto g_yz_yzzz = cbuffer.data(dg_off + 73 * ccomps * dcomps);

            auto g_yz_zzzz = cbuffer.data(dg_off + 74 * ccomps * dcomps);

            auto g_zz_xxxx = cbuffer.data(dg_off + 75 * ccomps * dcomps);

            auto g_zz_xxxy = cbuffer.data(dg_off + 76 * ccomps * dcomps);

            auto g_zz_xxxz = cbuffer.data(dg_off + 77 * ccomps * dcomps);

            auto g_zz_xxyy = cbuffer.data(dg_off + 78 * ccomps * dcomps);

            auto g_zz_xxyz = cbuffer.data(dg_off + 79 * ccomps * dcomps);

            auto g_zz_xxzz = cbuffer.data(dg_off + 80 * ccomps * dcomps);

            auto g_zz_xyyy = cbuffer.data(dg_off + 81 * ccomps * dcomps);

            auto g_zz_xyyz = cbuffer.data(dg_off + 82 * ccomps * dcomps);

            auto g_zz_xyzz = cbuffer.data(dg_off + 83 * ccomps * dcomps);

            auto g_zz_xzzz = cbuffer.data(dg_off + 84 * ccomps * dcomps);

            auto g_zz_yyyy = cbuffer.data(dg_off + 85 * ccomps * dcomps);

            auto g_zz_yyyz = cbuffer.data(dg_off + 86 * ccomps * dcomps);

            auto g_zz_yyzz = cbuffer.data(dg_off + 87 * ccomps * dcomps);

            auto g_zz_yzzz = cbuffer.data(dg_off + 88 * ccomps * dcomps);

            auto g_zz_zzzz = cbuffer.data(dg_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DGSS

            const auto dg_geom_01_off = idx_geom_01_dgxx + i * dcomps + j;

            auto g_0_x_xx_xxxx = cbuffer.data(dg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxy = cbuffer.data(dg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxz = cbuffer.data(dg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxyy = cbuffer.data(dg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxyz = cbuffer.data(dg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxzz = cbuffer.data(dg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xyyy = cbuffer.data(dg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xyyz = cbuffer.data(dg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xyzz = cbuffer.data(dg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xzzz = cbuffer.data(dg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_yyyy = cbuffer.data(dg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_yyyz = cbuffer.data(dg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_yyzz = cbuffer.data(dg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_yzzz = cbuffer.data(dg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_zzzz = cbuffer.data(dg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xy_xxxx = cbuffer.data(dg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xy_xxxy = cbuffer.data(dg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xy_xxxz = cbuffer.data(dg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xy_xxyy = cbuffer.data(dg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xy_xxyz = cbuffer.data(dg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xy_xxzz = cbuffer.data(dg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xy_xyyy = cbuffer.data(dg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xy_xyyz = cbuffer.data(dg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xy_xyzz = cbuffer.data(dg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xy_xzzz = cbuffer.data(dg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xy_yyyy = cbuffer.data(dg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xy_yyyz = cbuffer.data(dg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xy_yyzz = cbuffer.data(dg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xy_yzzz = cbuffer.data(dg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xy_zzzz = cbuffer.data(dg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xz_xxxx = cbuffer.data(dg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xz_xxxy = cbuffer.data(dg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xz_xxxz = cbuffer.data(dg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xz_xxyy = cbuffer.data(dg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xz_xxyz = cbuffer.data(dg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xz_xxzz = cbuffer.data(dg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xz_xyyy = cbuffer.data(dg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xz_xyyz = cbuffer.data(dg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xz_xyzz = cbuffer.data(dg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xz_xzzz = cbuffer.data(dg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xz_yyyy = cbuffer.data(dg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xz_yyyz = cbuffer.data(dg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xz_yyzz = cbuffer.data(dg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xz_yzzz = cbuffer.data(dg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xz_zzzz = cbuffer.data(dg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yy_xxxx = cbuffer.data(dg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yy_xxxy = cbuffer.data(dg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yy_xxxz = cbuffer.data(dg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yy_xxyy = cbuffer.data(dg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yy_xxyz = cbuffer.data(dg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yy_xxzz = cbuffer.data(dg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_yy_xyyy = cbuffer.data(dg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yy_xyyz = cbuffer.data(dg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yy_xyzz = cbuffer.data(dg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_yy_xzzz = cbuffer.data(dg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_yy_yyyy = cbuffer.data(dg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_yy_yyyz = cbuffer.data(dg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_yy_yyzz = cbuffer.data(dg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_yy_yzzz = cbuffer.data(dg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_yy_zzzz = cbuffer.data(dg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yz_xxxx = cbuffer.data(dg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yz_xxxy = cbuffer.data(dg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yz_xxxz = cbuffer.data(dg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yz_xxyy = cbuffer.data(dg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yz_xxyz = cbuffer.data(dg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yz_xxzz = cbuffer.data(dg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yz_xyyy = cbuffer.data(dg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yz_xyyz = cbuffer.data(dg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yz_xyzz = cbuffer.data(dg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yz_xzzz = cbuffer.data(dg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yz_yyyy = cbuffer.data(dg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yz_yyyz = cbuffer.data(dg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yz_yyzz = cbuffer.data(dg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yz_yzzz = cbuffer.data(dg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yz_zzzz = cbuffer.data(dg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_zz_xxxx = cbuffer.data(dg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_zz_xxxy = cbuffer.data(dg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_zz_xxxz = cbuffer.data(dg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_zz_xxyy = cbuffer.data(dg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_zz_xxyz = cbuffer.data(dg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_zz_xxzz = cbuffer.data(dg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_zz_xyyy = cbuffer.data(dg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_zz_xyyz = cbuffer.data(dg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_zz_xyzz = cbuffer.data(dg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_zz_xzzz = cbuffer.data(dg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zz_yyyy = cbuffer.data(dg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zz_yyyz = cbuffer.data(dg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zz_yyzz = cbuffer.data(dg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zz_yzzz = cbuffer.data(dg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zz_zzzz = cbuffer.data(dg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xx_xxxx = cbuffer.data(dg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xx_xxxy = cbuffer.data(dg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xx_xxxz = cbuffer.data(dg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xx_xxyy = cbuffer.data(dg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xx_xxyz = cbuffer.data(dg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xx_xxzz = cbuffer.data(dg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xx_xyyy = cbuffer.data(dg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xx_xyyz = cbuffer.data(dg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xx_xyzz = cbuffer.data(dg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xx_xzzz = cbuffer.data(dg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xx_yyyy = cbuffer.data(dg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xx_yyyz = cbuffer.data(dg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xx_yyzz = cbuffer.data(dg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xx_yzzz = cbuffer.data(dg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xx_zzzz = cbuffer.data(dg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xy_xxxx = cbuffer.data(dg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xy_xxxy = cbuffer.data(dg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xy_xxxz = cbuffer.data(dg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xy_xxyy = cbuffer.data(dg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xy_xxyz = cbuffer.data(dg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xy_xxzz = cbuffer.data(dg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xy_xyyy = cbuffer.data(dg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xy_xyyz = cbuffer.data(dg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xy_xyzz = cbuffer.data(dg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xy_xzzz = cbuffer.data(dg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xy_yyyy = cbuffer.data(dg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xy_yyyz = cbuffer.data(dg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xy_yyzz = cbuffer.data(dg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xy_yzzz = cbuffer.data(dg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xy_zzzz = cbuffer.data(dg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xz_xxxx = cbuffer.data(dg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xz_xxxy = cbuffer.data(dg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xz_xxxz = cbuffer.data(dg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xz_xxyy = cbuffer.data(dg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xz_xxyz = cbuffer.data(dg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xz_xxzz = cbuffer.data(dg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xz_xyyy = cbuffer.data(dg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xz_xyyz = cbuffer.data(dg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xz_xyzz = cbuffer.data(dg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xz_xzzz = cbuffer.data(dg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xz_yyyy = cbuffer.data(dg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xz_yyyz = cbuffer.data(dg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xz_yyzz = cbuffer.data(dg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xz_yzzz = cbuffer.data(dg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xz_zzzz = cbuffer.data(dg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_yy_xxxx = cbuffer.data(dg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_yy_xxxy = cbuffer.data(dg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_yy_xxxz = cbuffer.data(dg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_yy_xxyy = cbuffer.data(dg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_yy_xxyz = cbuffer.data(dg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_yy_xxzz = cbuffer.data(dg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_yy_xyyy = cbuffer.data(dg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_yy_xyyz = cbuffer.data(dg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_yy_xyzz = cbuffer.data(dg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_yy_xzzz = cbuffer.data(dg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_yy_yyyy = cbuffer.data(dg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_yy_yyyz = cbuffer.data(dg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_yy_yyzz = cbuffer.data(dg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_yy_yzzz = cbuffer.data(dg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_yy_zzzz = cbuffer.data(dg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yz_xxxx = cbuffer.data(dg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yz_xxxy = cbuffer.data(dg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yz_xxxz = cbuffer.data(dg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yz_xxyy = cbuffer.data(dg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yz_xxyz = cbuffer.data(dg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yz_xxzz = cbuffer.data(dg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yz_xyyy = cbuffer.data(dg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yz_xyyz = cbuffer.data(dg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yz_xyzz = cbuffer.data(dg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yz_xzzz = cbuffer.data(dg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yz_yyyy = cbuffer.data(dg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yz_yyyz = cbuffer.data(dg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yz_yyzz = cbuffer.data(dg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yz_yzzz = cbuffer.data(dg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yz_zzzz = cbuffer.data(dg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_zz_xxxx = cbuffer.data(dg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_zz_xxxy = cbuffer.data(dg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_zz_xxxz = cbuffer.data(dg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_zz_xxyy = cbuffer.data(dg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_zz_xxyz = cbuffer.data(dg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_zz_xxzz = cbuffer.data(dg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_zz_xyyy = cbuffer.data(dg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_zz_xyyz = cbuffer.data(dg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_zz_xyzz = cbuffer.data(dg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_zz_xzzz = cbuffer.data(dg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zz_yyyy = cbuffer.data(dg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zz_yyyz = cbuffer.data(dg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zz_yyzz = cbuffer.data(dg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zz_yzzz = cbuffer.data(dg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zz_zzzz = cbuffer.data(dg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xx_xxxx = cbuffer.data(dg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xx_xxxy = cbuffer.data(dg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xx_xxxz = cbuffer.data(dg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xx_xxyy = cbuffer.data(dg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xx_xxyz = cbuffer.data(dg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xx_xxzz = cbuffer.data(dg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xx_xyyy = cbuffer.data(dg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xx_xyyz = cbuffer.data(dg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xx_xyzz = cbuffer.data(dg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xx_xzzz = cbuffer.data(dg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xx_yyyy = cbuffer.data(dg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xx_yyyz = cbuffer.data(dg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xx_yyzz = cbuffer.data(dg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xx_yzzz = cbuffer.data(dg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xx_zzzz = cbuffer.data(dg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xy_xxxx = cbuffer.data(dg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xy_xxxy = cbuffer.data(dg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xy_xxxz = cbuffer.data(dg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xy_xxyy = cbuffer.data(dg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xy_xxyz = cbuffer.data(dg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xy_xxzz = cbuffer.data(dg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xy_xyyy = cbuffer.data(dg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xy_xyyz = cbuffer.data(dg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xy_xyzz = cbuffer.data(dg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xy_xzzz = cbuffer.data(dg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xy_yyyy = cbuffer.data(dg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xy_yyyz = cbuffer.data(dg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xy_yyzz = cbuffer.data(dg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xy_yzzz = cbuffer.data(dg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xy_zzzz = cbuffer.data(dg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xz_xxxx = cbuffer.data(dg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xz_xxxy = cbuffer.data(dg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xz_xxxz = cbuffer.data(dg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xz_xxyy = cbuffer.data(dg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xz_xxyz = cbuffer.data(dg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xz_xxzz = cbuffer.data(dg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xz_xyyy = cbuffer.data(dg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xz_xyyz = cbuffer.data(dg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xz_xyzz = cbuffer.data(dg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xz_xzzz = cbuffer.data(dg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xz_yyyy = cbuffer.data(dg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xz_yyyz = cbuffer.data(dg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xz_yyzz = cbuffer.data(dg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xz_yzzz = cbuffer.data(dg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xz_zzzz = cbuffer.data(dg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_yy_xxxx = cbuffer.data(dg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_yy_xxxy = cbuffer.data(dg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_yy_xxxz = cbuffer.data(dg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_yy_xxyy = cbuffer.data(dg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_yy_xxyz = cbuffer.data(dg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_yy_xxzz = cbuffer.data(dg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_yy_xyyy = cbuffer.data(dg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_yy_xyyz = cbuffer.data(dg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_yy_xyzz = cbuffer.data(dg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_yy_xzzz = cbuffer.data(dg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_yy_yyyy = cbuffer.data(dg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_yy_yyyz = cbuffer.data(dg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_yy_yyzz = cbuffer.data(dg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_yy_yzzz = cbuffer.data(dg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_yy_zzzz = cbuffer.data(dg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yz_xxxx = cbuffer.data(dg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yz_xxxy = cbuffer.data(dg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yz_xxxz = cbuffer.data(dg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yz_xxyy = cbuffer.data(dg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yz_xxyz = cbuffer.data(dg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yz_xxzz = cbuffer.data(dg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yz_xyyy = cbuffer.data(dg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yz_xyyz = cbuffer.data(dg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yz_xyzz = cbuffer.data(dg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yz_xzzz = cbuffer.data(dg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yz_yyyy = cbuffer.data(dg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yz_yyyz = cbuffer.data(dg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_yz_yyzz = cbuffer.data(dg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yz_yzzz = cbuffer.data(dg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yz_zzzz = cbuffer.data(dg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_zz_xxxx = cbuffer.data(dg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_zz_xxxy = cbuffer.data(dg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_zz_xxxz = cbuffer.data(dg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_zz_xxyy = cbuffer.data(dg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_zz_xxyz = cbuffer.data(dg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_zz_xxzz = cbuffer.data(dg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_zz_xyyy = cbuffer.data(dg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_zz_xyyz = cbuffer.data(dg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_zz_xyzz = cbuffer.data(dg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_zz_xzzz = cbuffer.data(dg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zz_yyyy = cbuffer.data(dg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zz_yyyz = cbuffer.data(dg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zz_yyzz = cbuffer.data(dg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zz_yzzz = cbuffer.data(dg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zz_zzzz = cbuffer.data(dg_geom_01_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_01_off = idx_geom_01_dhxx + i * dcomps + j;

            auto g_0_x_xx_xxxxx = cbuffer.data(dh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxxy = cbuffer.data(dh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxxz = cbuffer.data(dh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxxyy = cbuffer.data(dh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxxyz = cbuffer.data(dh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxxzz = cbuffer.data(dh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xxyyy = cbuffer.data(dh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xxyyz = cbuffer.data(dh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xxyzz = cbuffer.data(dh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xxzzz = cbuffer.data(dh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_xyyyy = cbuffer.data(dh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_xyyyz = cbuffer.data(dh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_xyyzz = cbuffer.data(dh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_xyzzz = cbuffer.data(dh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_xzzzz = cbuffer.data(dh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xx_yyyyy = cbuffer.data(dh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xx_yyyyz = cbuffer.data(dh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xx_yyyzz = cbuffer.data(dh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xx_yyzzz = cbuffer.data(dh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xx_yzzzz = cbuffer.data(dh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xx_zzzzz = cbuffer.data(dh_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xy_xxxxx = cbuffer.data(dh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xy_xxxxy = cbuffer.data(dh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xy_xxxxz = cbuffer.data(dh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xy_xxxyy = cbuffer.data(dh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xy_xxxyz = cbuffer.data(dh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xy_xxxzz = cbuffer.data(dh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xy_xxyyy = cbuffer.data(dh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xy_xxyyz = cbuffer.data(dh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xy_xxyzz = cbuffer.data(dh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xy_xxzzz = cbuffer.data(dh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xy_xyyyy = cbuffer.data(dh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xy_xyyyz = cbuffer.data(dh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xy_xyyzz = cbuffer.data(dh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xy_xyzzz = cbuffer.data(dh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xy_xzzzz = cbuffer.data(dh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xy_yyyyy = cbuffer.data(dh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xy_yyyyz = cbuffer.data(dh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xy_yyyzz = cbuffer.data(dh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xy_yyzzz = cbuffer.data(dh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xy_yzzzz = cbuffer.data(dh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xy_zzzzz = cbuffer.data(dh_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xz_xxxxx = cbuffer.data(dh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xz_xxxxy = cbuffer.data(dh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xz_xxxxz = cbuffer.data(dh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xz_xxxyy = cbuffer.data(dh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xz_xxxyz = cbuffer.data(dh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xz_xxxzz = cbuffer.data(dh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xz_xxyyy = cbuffer.data(dh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xz_xxyyz = cbuffer.data(dh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xz_xxyzz = cbuffer.data(dh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xz_xxzzz = cbuffer.data(dh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xz_xyyyy = cbuffer.data(dh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xz_xyyyz = cbuffer.data(dh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xz_xyyzz = cbuffer.data(dh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xz_xyzzz = cbuffer.data(dh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xz_xzzzz = cbuffer.data(dh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xz_yyyyy = cbuffer.data(dh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xz_yyyyz = cbuffer.data(dh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xz_yyyzz = cbuffer.data(dh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xz_yyzzz = cbuffer.data(dh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xz_yzzzz = cbuffer.data(dh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xz_zzzzz = cbuffer.data(dh_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yy_xxxxx = cbuffer.data(dh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yy_xxxxy = cbuffer.data(dh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yy_xxxxz = cbuffer.data(dh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yy_xxxyy = cbuffer.data(dh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yy_xxxyz = cbuffer.data(dh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yy_xxxzz = cbuffer.data(dh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yy_xxyyy = cbuffer.data(dh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yy_xxyyz = cbuffer.data(dh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yy_xxyzz = cbuffer.data(dh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yy_xxzzz = cbuffer.data(dh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yy_xyyyy = cbuffer.data(dh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yy_xyyyz = cbuffer.data(dh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yy_xyyzz = cbuffer.data(dh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yy_xyzzz = cbuffer.data(dh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yy_xzzzz = cbuffer.data(dh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yy_yyyyy = cbuffer.data(dh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yy_yyyyz = cbuffer.data(dh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yy_yyyzz = cbuffer.data(dh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yy_yyzzz = cbuffer.data(dh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yy_yzzzz = cbuffer.data(dh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yy_zzzzz = cbuffer.data(dh_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_yz_xxxxx = cbuffer.data(dh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_yz_xxxxy = cbuffer.data(dh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_yz_xxxxz = cbuffer.data(dh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_yz_xxxyy = cbuffer.data(dh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_yz_xxxyz = cbuffer.data(dh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_yz_xxxzz = cbuffer.data(dh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yz_xxyyy = cbuffer.data(dh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yz_xxyyz = cbuffer.data(dh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yz_xxyzz = cbuffer.data(dh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yz_xxzzz = cbuffer.data(dh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yz_xyyyy = cbuffer.data(dh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yz_xyyyz = cbuffer.data(dh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yz_xyyzz = cbuffer.data(dh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yz_xyzzz = cbuffer.data(dh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yz_xzzzz = cbuffer.data(dh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yz_yyyyy = cbuffer.data(dh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yz_yyyyz = cbuffer.data(dh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yz_yyyzz = cbuffer.data(dh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yz_yyzzz = cbuffer.data(dh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yz_yzzzz = cbuffer.data(dh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yz_zzzzz = cbuffer.data(dh_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_zz_xxxxx = cbuffer.data(dh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_zz_xxxxy = cbuffer.data(dh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_zz_xxxxz = cbuffer.data(dh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_zz_xxxyy = cbuffer.data(dh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_zz_xxxyz = cbuffer.data(dh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_zz_xxxzz = cbuffer.data(dh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_zz_xxyyy = cbuffer.data(dh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_zz_xxyyz = cbuffer.data(dh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_zz_xxyzz = cbuffer.data(dh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_zz_xxzzz = cbuffer.data(dh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_zz_xyyyy = cbuffer.data(dh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_zz_xyyyz = cbuffer.data(dh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_zz_xyyzz = cbuffer.data(dh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_zz_xyzzz = cbuffer.data(dh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_zz_xzzzz = cbuffer.data(dh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_zz_yyyyy = cbuffer.data(dh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_zz_yyyyz = cbuffer.data(dh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_zz_yyyzz = cbuffer.data(dh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_zz_yyzzz = cbuffer.data(dh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_zz_yzzzz = cbuffer.data(dh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_zz_zzzzz = cbuffer.data(dh_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xx_xxxxx = cbuffer.data(dh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xx_xxxxy = cbuffer.data(dh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xx_xxxxz = cbuffer.data(dh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xx_xxxyy = cbuffer.data(dh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xx_xxxyz = cbuffer.data(dh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xx_xxxzz = cbuffer.data(dh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xx_xxyyy = cbuffer.data(dh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xx_xxyyz = cbuffer.data(dh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xx_xxyzz = cbuffer.data(dh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xx_xxzzz = cbuffer.data(dh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xx_xyyyy = cbuffer.data(dh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xx_xyyyz = cbuffer.data(dh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xx_xyyzz = cbuffer.data(dh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xx_xyzzz = cbuffer.data(dh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xx_xzzzz = cbuffer.data(dh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xx_yyyyy = cbuffer.data(dh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xx_yyyyz = cbuffer.data(dh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xx_yyyzz = cbuffer.data(dh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xx_yyzzz = cbuffer.data(dh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xx_yzzzz = cbuffer.data(dh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xx_zzzzz = cbuffer.data(dh_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xy_xxxxx = cbuffer.data(dh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xy_xxxxy = cbuffer.data(dh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xy_xxxxz = cbuffer.data(dh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xy_xxxyy = cbuffer.data(dh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xy_xxxyz = cbuffer.data(dh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xy_xxxzz = cbuffer.data(dh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xy_xxyyy = cbuffer.data(dh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xy_xxyyz = cbuffer.data(dh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xy_xxyzz = cbuffer.data(dh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xy_xxzzz = cbuffer.data(dh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xy_xyyyy = cbuffer.data(dh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xy_xyyyz = cbuffer.data(dh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xy_xyyzz = cbuffer.data(dh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xy_xyzzz = cbuffer.data(dh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xy_xzzzz = cbuffer.data(dh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xy_yyyyy = cbuffer.data(dh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xy_yyyyz = cbuffer.data(dh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xy_yyyzz = cbuffer.data(dh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xy_yyzzz = cbuffer.data(dh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xy_yzzzz = cbuffer.data(dh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xy_zzzzz = cbuffer.data(dh_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xz_xxxxx = cbuffer.data(dh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xz_xxxxy = cbuffer.data(dh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xz_xxxxz = cbuffer.data(dh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xz_xxxyy = cbuffer.data(dh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xz_xxxyz = cbuffer.data(dh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xz_xxxzz = cbuffer.data(dh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xz_xxyyy = cbuffer.data(dh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xz_xxyyz = cbuffer.data(dh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xz_xxyzz = cbuffer.data(dh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xz_xxzzz = cbuffer.data(dh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xz_xyyyy = cbuffer.data(dh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xz_xyyyz = cbuffer.data(dh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xz_xyyzz = cbuffer.data(dh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xz_xyzzz = cbuffer.data(dh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xz_xzzzz = cbuffer.data(dh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xz_yyyyy = cbuffer.data(dh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xz_yyyyz = cbuffer.data(dh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xz_yyyzz = cbuffer.data(dh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xz_yyzzz = cbuffer.data(dh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xz_yzzzz = cbuffer.data(dh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xz_zzzzz = cbuffer.data(dh_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_yy_xxxxx = cbuffer.data(dh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_yy_xxxxy = cbuffer.data(dh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_yy_xxxxz = cbuffer.data(dh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_yy_xxxyy = cbuffer.data(dh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_yy_xxxyz = cbuffer.data(dh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_yy_xxxzz = cbuffer.data(dh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_yy_xxyyy = cbuffer.data(dh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_yy_xxyyz = cbuffer.data(dh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_yy_xxyzz = cbuffer.data(dh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_yy_xxzzz = cbuffer.data(dh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_yy_xyyyy = cbuffer.data(dh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_yy_xyyyz = cbuffer.data(dh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_yy_xyyzz = cbuffer.data(dh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_yy_xyzzz = cbuffer.data(dh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_yy_xzzzz = cbuffer.data(dh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_yy_yyyyy = cbuffer.data(dh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_yy_yyyyz = cbuffer.data(dh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_yy_yyyzz = cbuffer.data(dh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_yy_yyzzz = cbuffer.data(dh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_yy_yzzzz = cbuffer.data(dh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_yy_zzzzz = cbuffer.data(dh_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_yz_xxxxx = cbuffer.data(dh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_yz_xxxxy = cbuffer.data(dh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_yz_xxxxz = cbuffer.data(dh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_yz_xxxyy = cbuffer.data(dh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_yz_xxxyz = cbuffer.data(dh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_yz_xxxzz = cbuffer.data(dh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_yz_xxyyy = cbuffer.data(dh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_yz_xxyyz = cbuffer.data(dh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_yz_xxyzz = cbuffer.data(dh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_yz_xxzzz = cbuffer.data(dh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_yz_xyyyy = cbuffer.data(dh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_yz_xyyyz = cbuffer.data(dh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_yz_xyyzz = cbuffer.data(dh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_yz_xyzzz = cbuffer.data(dh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_yz_xzzzz = cbuffer.data(dh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_yz_yyyyy = cbuffer.data(dh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_yz_yyyyz = cbuffer.data(dh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_yz_yyyzz = cbuffer.data(dh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_yz_yyzzz = cbuffer.data(dh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_yz_yzzzz = cbuffer.data(dh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_yz_zzzzz = cbuffer.data(dh_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_zz_xxxxx = cbuffer.data(dh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_zz_xxxxy = cbuffer.data(dh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_zz_xxxxz = cbuffer.data(dh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_zz_xxxyy = cbuffer.data(dh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_zz_xxxyz = cbuffer.data(dh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_zz_xxxzz = cbuffer.data(dh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_zz_xxyyy = cbuffer.data(dh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_zz_xxyyz = cbuffer.data(dh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_zz_xxyzz = cbuffer.data(dh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_zz_xxzzz = cbuffer.data(dh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_zz_xyyyy = cbuffer.data(dh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_zz_xyyyz = cbuffer.data(dh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_zz_xyyzz = cbuffer.data(dh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_zz_xyzzz = cbuffer.data(dh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_zz_xzzzz = cbuffer.data(dh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_zz_yyyyy = cbuffer.data(dh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_zz_yyyyz = cbuffer.data(dh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_zz_yyyzz = cbuffer.data(dh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_zz_yyzzz = cbuffer.data(dh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_zz_yzzzz = cbuffer.data(dh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_zz_zzzzz = cbuffer.data(dh_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_xx_xxxxx = cbuffer.data(dh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xx_xxxxy = cbuffer.data(dh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xx_xxxxz = cbuffer.data(dh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xx_xxxyy = cbuffer.data(dh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xx_xxxyz = cbuffer.data(dh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xx_xxxzz = cbuffer.data(dh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xx_xxyyy = cbuffer.data(dh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xx_xxyyz = cbuffer.data(dh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xx_xxyzz = cbuffer.data(dh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_xx_xxzzz = cbuffer.data(dh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xx_xyyyy = cbuffer.data(dh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xx_xyyyz = cbuffer.data(dh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_xx_xyyzz = cbuffer.data(dh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xx_xyzzz = cbuffer.data(dh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xx_xzzzz = cbuffer.data(dh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_xx_yyyyy = cbuffer.data(dh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xx_yyyyz = cbuffer.data(dh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xx_yyyzz = cbuffer.data(dh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_xx_yyzzz = cbuffer.data(dh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xx_yzzzz = cbuffer.data(dh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xx_zzzzz = cbuffer.data(dh_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_xy_xxxxx = cbuffer.data(dh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xy_xxxxy = cbuffer.data(dh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xy_xxxxz = cbuffer.data(dh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_xy_xxxyy = cbuffer.data(dh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xy_xxxyz = cbuffer.data(dh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xy_xxxzz = cbuffer.data(dh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_xy_xxyyy = cbuffer.data(dh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xy_xxyyz = cbuffer.data(dh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xy_xxyzz = cbuffer.data(dh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_xy_xxzzz = cbuffer.data(dh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xy_xyyyy = cbuffer.data(dh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xy_xyyyz = cbuffer.data(dh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_xy_xyyzz = cbuffer.data(dh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xy_xyzzz = cbuffer.data(dh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xy_xzzzz = cbuffer.data(dh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_xy_yyyyy = cbuffer.data(dh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xy_yyyyz = cbuffer.data(dh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xy_yyyzz = cbuffer.data(dh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_xy_yyzzz = cbuffer.data(dh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xy_yzzzz = cbuffer.data(dh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xy_zzzzz = cbuffer.data(dh_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_xz_xxxxx = cbuffer.data(dh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xz_xxxxy = cbuffer.data(dh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xz_xxxxz = cbuffer.data(dh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_xz_xxxyy = cbuffer.data(dh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xz_xxxyz = cbuffer.data(dh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xz_xxxzz = cbuffer.data(dh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xz_xxyyy = cbuffer.data(dh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xz_xxyyz = cbuffer.data(dh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xz_xxyzz = cbuffer.data(dh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xz_xxzzz = cbuffer.data(dh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xz_xyyyy = cbuffer.data(dh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xz_xyyyz = cbuffer.data(dh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xz_xyyzz = cbuffer.data(dh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xz_xyzzz = cbuffer.data(dh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xz_xzzzz = cbuffer.data(dh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xz_yyyyy = cbuffer.data(dh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xz_yyyyz = cbuffer.data(dh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xz_yyyzz = cbuffer.data(dh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xz_yyzzz = cbuffer.data(dh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xz_yzzzz = cbuffer.data(dh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xz_zzzzz = cbuffer.data(dh_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_yy_xxxxx = cbuffer.data(dh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_yy_xxxxy = cbuffer.data(dh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_yy_xxxxz = cbuffer.data(dh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_yy_xxxyy = cbuffer.data(dh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_yy_xxxyz = cbuffer.data(dh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_yy_xxxzz = cbuffer.data(dh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_yy_xxyyy = cbuffer.data(dh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_yy_xxyyz = cbuffer.data(dh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_yy_xxyzz = cbuffer.data(dh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_yy_xxzzz = cbuffer.data(dh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_yy_xyyyy = cbuffer.data(dh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_yy_xyyyz = cbuffer.data(dh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_yy_xyyzz = cbuffer.data(dh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_yy_xyzzz = cbuffer.data(dh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_yy_xzzzz = cbuffer.data(dh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_yy_yyyyy = cbuffer.data(dh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_yy_yyyyz = cbuffer.data(dh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_yy_yyyzz = cbuffer.data(dh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_yy_yyzzz = cbuffer.data(dh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_yy_yzzzz = cbuffer.data(dh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_yy_zzzzz = cbuffer.data(dh_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_yz_xxxxx = cbuffer.data(dh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_yz_xxxxy = cbuffer.data(dh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_yz_xxxxz = cbuffer.data(dh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_yz_xxxyy = cbuffer.data(dh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_yz_xxxyz = cbuffer.data(dh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_yz_xxxzz = cbuffer.data(dh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_yz_xxyyy = cbuffer.data(dh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_yz_xxyyz = cbuffer.data(dh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_yz_xxyzz = cbuffer.data(dh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_yz_xxzzz = cbuffer.data(dh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_yz_xyyyy = cbuffer.data(dh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_yz_xyyyz = cbuffer.data(dh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_yz_xyyzz = cbuffer.data(dh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_yz_xyzzz = cbuffer.data(dh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_yz_xzzzz = cbuffer.data(dh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_yz_yyyyy = cbuffer.data(dh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_yz_yyyyz = cbuffer.data(dh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_yz_yyyzz = cbuffer.data(dh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_yz_yyzzz = cbuffer.data(dh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_yz_yzzzz = cbuffer.data(dh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_yz_zzzzz = cbuffer.data(dh_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_zz_xxxxx = cbuffer.data(dh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_zz_xxxxy = cbuffer.data(dh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_zz_xxxxz = cbuffer.data(dh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_zz_xxxyy = cbuffer.data(dh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_zz_xxxyz = cbuffer.data(dh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_zz_xxxzz = cbuffer.data(dh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_zz_xxyyy = cbuffer.data(dh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_zz_xxyyz = cbuffer.data(dh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_zz_xxyzz = cbuffer.data(dh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_zz_xxzzz = cbuffer.data(dh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_zz_xyyyy = cbuffer.data(dh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_zz_xyyyz = cbuffer.data(dh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_zz_xyyzz = cbuffer.data(dh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_zz_xyzzz = cbuffer.data(dh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_zz_xzzzz = cbuffer.data(dh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_zz_yyyyy = cbuffer.data(dh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_zz_yyyyz = cbuffer.data(dh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_zz_yyyzz = cbuffer.data(dh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_zz_yyzzz = cbuffer.data(dh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_zz_yzzzz = cbuffer.data(dh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_zz_zzzzz = cbuffer.data(dh_geom_01_off + 377 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fgxx

            const auto fg_geom_01_off = idx_geom_01_fgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxx_xxxx = cbuffer.data(fg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxy = cbuffer.data(fg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxz = cbuffer.data(fg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxyy = cbuffer.data(fg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxyz = cbuffer.data(fg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxzz = cbuffer.data(fg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xyyy = cbuffer.data(fg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xyyz = cbuffer.data(fg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xyzz = cbuffer.data(fg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xzzz = cbuffer.data(fg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_yyyy = cbuffer.data(fg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_yyyz = cbuffer.data(fg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_yyzz = cbuffer.data(fg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_yzzz = cbuffer.data(fg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_zzzz = cbuffer.data(fg_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_xxxx, g_0_x_xx_xxxxx, g_0_x_xx_xxxxy, g_0_x_xx_xxxxz, g_0_x_xx_xxxy, g_0_x_xx_xxxyy, g_0_x_xx_xxxyz, g_0_x_xx_xxxz, g_0_x_xx_xxxzz, g_0_x_xx_xxyy, g_0_x_xx_xxyyy, g_0_x_xx_xxyyz, g_0_x_xx_xxyz, g_0_x_xx_xxyzz, g_0_x_xx_xxzz, g_0_x_xx_xxzzz, g_0_x_xx_xyyy, g_0_x_xx_xyyyy, g_0_x_xx_xyyyz, g_0_x_xx_xyyz, g_0_x_xx_xyyzz, g_0_x_xx_xyzz, g_0_x_xx_xyzzz, g_0_x_xx_xzzz, g_0_x_xx_xzzzz, g_0_x_xx_yyyy, g_0_x_xx_yyyz, g_0_x_xx_yyzz, g_0_x_xx_yzzz, g_0_x_xx_zzzz, g_0_x_xxx_xxxx, g_0_x_xxx_xxxy, g_0_x_xxx_xxxz, g_0_x_xxx_xxyy, g_0_x_xxx_xxyz, g_0_x_xxx_xxzz, g_0_x_xxx_xyyy, g_0_x_xxx_xyyz, g_0_x_xxx_xyzz, g_0_x_xxx_xzzz, g_0_x_xxx_yyyy, g_0_x_xxx_yyyz, g_0_x_xxx_yyzz, g_0_x_xxx_yzzz, g_0_x_xxx_zzzz, g_xx_xxxx, g_xx_xxxy, g_xx_xxxz, g_xx_xxyy, g_xx_xxyz, g_xx_xxzz, g_xx_xyyy, g_xx_xyyz, g_xx_xyzz, g_xx_xzzz, g_xx_yyyy, g_xx_yyyz, g_xx_yyzz, g_xx_yzzz, g_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxx_xxxx[k] = g_xx_xxxx[k] - g_0_x_xx_xxxx[k] * ab_x + g_0_x_xx_xxxxx[k];

                g_0_x_xxx_xxxy[k] = g_xx_xxxy[k] - g_0_x_xx_xxxy[k] * ab_x + g_0_x_xx_xxxxy[k];

                g_0_x_xxx_xxxz[k] = g_xx_xxxz[k] - g_0_x_xx_xxxz[k] * ab_x + g_0_x_xx_xxxxz[k];

                g_0_x_xxx_xxyy[k] = g_xx_xxyy[k] - g_0_x_xx_xxyy[k] * ab_x + g_0_x_xx_xxxyy[k];

                g_0_x_xxx_xxyz[k] = g_xx_xxyz[k] - g_0_x_xx_xxyz[k] * ab_x + g_0_x_xx_xxxyz[k];

                g_0_x_xxx_xxzz[k] = g_xx_xxzz[k] - g_0_x_xx_xxzz[k] * ab_x + g_0_x_xx_xxxzz[k];

                g_0_x_xxx_xyyy[k] = g_xx_xyyy[k] - g_0_x_xx_xyyy[k] * ab_x + g_0_x_xx_xxyyy[k];

                g_0_x_xxx_xyyz[k] = g_xx_xyyz[k] - g_0_x_xx_xyyz[k] * ab_x + g_0_x_xx_xxyyz[k];

                g_0_x_xxx_xyzz[k] = g_xx_xyzz[k] - g_0_x_xx_xyzz[k] * ab_x + g_0_x_xx_xxyzz[k];

                g_0_x_xxx_xzzz[k] = g_xx_xzzz[k] - g_0_x_xx_xzzz[k] * ab_x + g_0_x_xx_xxzzz[k];

                g_0_x_xxx_yyyy[k] = g_xx_yyyy[k] - g_0_x_xx_yyyy[k] * ab_x + g_0_x_xx_xyyyy[k];

                g_0_x_xxx_yyyz[k] = g_xx_yyyz[k] - g_0_x_xx_yyyz[k] * ab_x + g_0_x_xx_xyyyz[k];

                g_0_x_xxx_yyzz[k] = g_xx_yyzz[k] - g_0_x_xx_yyzz[k] * ab_x + g_0_x_xx_xyyzz[k];

                g_0_x_xxx_yzzz[k] = g_xx_yzzz[k] - g_0_x_xx_yzzz[k] * ab_x + g_0_x_xx_xyzzz[k];

                g_0_x_xxx_zzzz[k] = g_xx_zzzz[k] - g_0_x_xx_zzzz[k] * ab_x + g_0_x_xx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxy_xxxx = cbuffer.data(fg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxy_xxxy = cbuffer.data(fg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxy_xxxz = cbuffer.data(fg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxy_xxyy = cbuffer.data(fg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxy_xxyz = cbuffer.data(fg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxy_xxzz = cbuffer.data(fg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxy_xyyy = cbuffer.data(fg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxy_xyyz = cbuffer.data(fg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxy_xyzz = cbuffer.data(fg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxy_xzzz = cbuffer.data(fg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxy_yyyy = cbuffer.data(fg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxy_yyyz = cbuffer.data(fg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxy_yyzz = cbuffer.data(fg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_yzzz = cbuffer.data(fg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_zzzz = cbuffer.data(fg_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_xxxx, g_0_x_xx_xxxxy, g_0_x_xx_xxxy, g_0_x_xx_xxxyy, g_0_x_xx_xxxyz, g_0_x_xx_xxxz, g_0_x_xx_xxyy, g_0_x_xx_xxyyy, g_0_x_xx_xxyyz, g_0_x_xx_xxyz, g_0_x_xx_xxyzz, g_0_x_xx_xxzz, g_0_x_xx_xyyy, g_0_x_xx_xyyyy, g_0_x_xx_xyyyz, g_0_x_xx_xyyz, g_0_x_xx_xyyzz, g_0_x_xx_xyzz, g_0_x_xx_xyzzz, g_0_x_xx_xzzz, g_0_x_xx_yyyy, g_0_x_xx_yyyyy, g_0_x_xx_yyyyz, g_0_x_xx_yyyz, g_0_x_xx_yyyzz, g_0_x_xx_yyzz, g_0_x_xx_yyzzz, g_0_x_xx_yzzz, g_0_x_xx_yzzzz, g_0_x_xx_zzzz, g_0_x_xxy_xxxx, g_0_x_xxy_xxxy, g_0_x_xxy_xxxz, g_0_x_xxy_xxyy, g_0_x_xxy_xxyz, g_0_x_xxy_xxzz, g_0_x_xxy_xyyy, g_0_x_xxy_xyyz, g_0_x_xxy_xyzz, g_0_x_xxy_xzzz, g_0_x_xxy_yyyy, g_0_x_xxy_yyyz, g_0_x_xxy_yyzz, g_0_x_xxy_yzzz, g_0_x_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxy_xxxx[k] = -g_0_x_xx_xxxx[k] * ab_y + g_0_x_xx_xxxxy[k];

                g_0_x_xxy_xxxy[k] = -g_0_x_xx_xxxy[k] * ab_y + g_0_x_xx_xxxyy[k];

                g_0_x_xxy_xxxz[k] = -g_0_x_xx_xxxz[k] * ab_y + g_0_x_xx_xxxyz[k];

                g_0_x_xxy_xxyy[k] = -g_0_x_xx_xxyy[k] * ab_y + g_0_x_xx_xxyyy[k];

                g_0_x_xxy_xxyz[k] = -g_0_x_xx_xxyz[k] * ab_y + g_0_x_xx_xxyyz[k];

                g_0_x_xxy_xxzz[k] = -g_0_x_xx_xxzz[k] * ab_y + g_0_x_xx_xxyzz[k];

                g_0_x_xxy_xyyy[k] = -g_0_x_xx_xyyy[k] * ab_y + g_0_x_xx_xyyyy[k];

                g_0_x_xxy_xyyz[k] = -g_0_x_xx_xyyz[k] * ab_y + g_0_x_xx_xyyyz[k];

                g_0_x_xxy_xyzz[k] = -g_0_x_xx_xyzz[k] * ab_y + g_0_x_xx_xyyzz[k];

                g_0_x_xxy_xzzz[k] = -g_0_x_xx_xzzz[k] * ab_y + g_0_x_xx_xyzzz[k];

                g_0_x_xxy_yyyy[k] = -g_0_x_xx_yyyy[k] * ab_y + g_0_x_xx_yyyyy[k];

                g_0_x_xxy_yyyz[k] = -g_0_x_xx_yyyz[k] * ab_y + g_0_x_xx_yyyyz[k];

                g_0_x_xxy_yyzz[k] = -g_0_x_xx_yyzz[k] * ab_y + g_0_x_xx_yyyzz[k];

                g_0_x_xxy_yzzz[k] = -g_0_x_xx_yzzz[k] * ab_y + g_0_x_xx_yyzzz[k];

                g_0_x_xxy_zzzz[k] = -g_0_x_xx_zzzz[k] * ab_y + g_0_x_xx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxz_xxxx = cbuffer.data(fg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxz_xxxy = cbuffer.data(fg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxz_xxxz = cbuffer.data(fg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxz_xxyy = cbuffer.data(fg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxz_xxyz = cbuffer.data(fg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxz_xxzz = cbuffer.data(fg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxz_xyyy = cbuffer.data(fg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxz_xyyz = cbuffer.data(fg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxz_xyzz = cbuffer.data(fg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxz_xzzz = cbuffer.data(fg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxz_yyyy = cbuffer.data(fg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxz_yyyz = cbuffer.data(fg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxz_yyzz = cbuffer.data(fg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxz_yzzz = cbuffer.data(fg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxz_zzzz = cbuffer.data(fg_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_xxxx, g_0_x_xx_xxxxz, g_0_x_xx_xxxy, g_0_x_xx_xxxyz, g_0_x_xx_xxxz, g_0_x_xx_xxxzz, g_0_x_xx_xxyy, g_0_x_xx_xxyyz, g_0_x_xx_xxyz, g_0_x_xx_xxyzz, g_0_x_xx_xxzz, g_0_x_xx_xxzzz, g_0_x_xx_xyyy, g_0_x_xx_xyyyz, g_0_x_xx_xyyz, g_0_x_xx_xyyzz, g_0_x_xx_xyzz, g_0_x_xx_xyzzz, g_0_x_xx_xzzz, g_0_x_xx_xzzzz, g_0_x_xx_yyyy, g_0_x_xx_yyyyz, g_0_x_xx_yyyz, g_0_x_xx_yyyzz, g_0_x_xx_yyzz, g_0_x_xx_yyzzz, g_0_x_xx_yzzz, g_0_x_xx_yzzzz, g_0_x_xx_zzzz, g_0_x_xx_zzzzz, g_0_x_xxz_xxxx, g_0_x_xxz_xxxy, g_0_x_xxz_xxxz, g_0_x_xxz_xxyy, g_0_x_xxz_xxyz, g_0_x_xxz_xxzz, g_0_x_xxz_xyyy, g_0_x_xxz_xyyz, g_0_x_xxz_xyzz, g_0_x_xxz_xzzz, g_0_x_xxz_yyyy, g_0_x_xxz_yyyz, g_0_x_xxz_yyzz, g_0_x_xxz_yzzz, g_0_x_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxz_xxxx[k] = -g_0_x_xx_xxxx[k] * ab_z + g_0_x_xx_xxxxz[k];

                g_0_x_xxz_xxxy[k] = -g_0_x_xx_xxxy[k] * ab_z + g_0_x_xx_xxxyz[k];

                g_0_x_xxz_xxxz[k] = -g_0_x_xx_xxxz[k] * ab_z + g_0_x_xx_xxxzz[k];

                g_0_x_xxz_xxyy[k] = -g_0_x_xx_xxyy[k] * ab_z + g_0_x_xx_xxyyz[k];

                g_0_x_xxz_xxyz[k] = -g_0_x_xx_xxyz[k] * ab_z + g_0_x_xx_xxyzz[k];

                g_0_x_xxz_xxzz[k] = -g_0_x_xx_xxzz[k] * ab_z + g_0_x_xx_xxzzz[k];

                g_0_x_xxz_xyyy[k] = -g_0_x_xx_xyyy[k] * ab_z + g_0_x_xx_xyyyz[k];

                g_0_x_xxz_xyyz[k] = -g_0_x_xx_xyyz[k] * ab_z + g_0_x_xx_xyyzz[k];

                g_0_x_xxz_xyzz[k] = -g_0_x_xx_xyzz[k] * ab_z + g_0_x_xx_xyzzz[k];

                g_0_x_xxz_xzzz[k] = -g_0_x_xx_xzzz[k] * ab_z + g_0_x_xx_xzzzz[k];

                g_0_x_xxz_yyyy[k] = -g_0_x_xx_yyyy[k] * ab_z + g_0_x_xx_yyyyz[k];

                g_0_x_xxz_yyyz[k] = -g_0_x_xx_yyyz[k] * ab_z + g_0_x_xx_yyyzz[k];

                g_0_x_xxz_yyzz[k] = -g_0_x_xx_yyzz[k] * ab_z + g_0_x_xx_yyzzz[k];

                g_0_x_xxz_yzzz[k] = -g_0_x_xx_yzzz[k] * ab_z + g_0_x_xx_yzzzz[k];

                g_0_x_xxz_zzzz[k] = -g_0_x_xx_zzzz[k] * ab_z + g_0_x_xx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyy_xxxx = cbuffer.data(fg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyy_xxxy = cbuffer.data(fg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyy_xxxz = cbuffer.data(fg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyy_xxyy = cbuffer.data(fg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyy_xxyz = cbuffer.data(fg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyy_xxzz = cbuffer.data(fg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyy_xyyy = cbuffer.data(fg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyy_xyyz = cbuffer.data(fg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyy_xyzz = cbuffer.data(fg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xyy_xzzz = cbuffer.data(fg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xyy_yyyy = cbuffer.data(fg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xyy_yyyz = cbuffer.data(fg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xyy_yyzz = cbuffer.data(fg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xyy_yzzz = cbuffer.data(fg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xyy_zzzz = cbuffer.data(fg_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xy_xxxx, g_0_x_xy_xxxxy, g_0_x_xy_xxxy, g_0_x_xy_xxxyy, g_0_x_xy_xxxyz, g_0_x_xy_xxxz, g_0_x_xy_xxyy, g_0_x_xy_xxyyy, g_0_x_xy_xxyyz, g_0_x_xy_xxyz, g_0_x_xy_xxyzz, g_0_x_xy_xxzz, g_0_x_xy_xyyy, g_0_x_xy_xyyyy, g_0_x_xy_xyyyz, g_0_x_xy_xyyz, g_0_x_xy_xyyzz, g_0_x_xy_xyzz, g_0_x_xy_xyzzz, g_0_x_xy_xzzz, g_0_x_xy_yyyy, g_0_x_xy_yyyyy, g_0_x_xy_yyyyz, g_0_x_xy_yyyz, g_0_x_xy_yyyzz, g_0_x_xy_yyzz, g_0_x_xy_yyzzz, g_0_x_xy_yzzz, g_0_x_xy_yzzzz, g_0_x_xy_zzzz, g_0_x_xyy_xxxx, g_0_x_xyy_xxxy, g_0_x_xyy_xxxz, g_0_x_xyy_xxyy, g_0_x_xyy_xxyz, g_0_x_xyy_xxzz, g_0_x_xyy_xyyy, g_0_x_xyy_xyyz, g_0_x_xyy_xyzz, g_0_x_xyy_xzzz, g_0_x_xyy_yyyy, g_0_x_xyy_yyyz, g_0_x_xyy_yyzz, g_0_x_xyy_yzzz, g_0_x_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyy_xxxx[k] = -g_0_x_xy_xxxx[k] * ab_y + g_0_x_xy_xxxxy[k];

                g_0_x_xyy_xxxy[k] = -g_0_x_xy_xxxy[k] * ab_y + g_0_x_xy_xxxyy[k];

                g_0_x_xyy_xxxz[k] = -g_0_x_xy_xxxz[k] * ab_y + g_0_x_xy_xxxyz[k];

                g_0_x_xyy_xxyy[k] = -g_0_x_xy_xxyy[k] * ab_y + g_0_x_xy_xxyyy[k];

                g_0_x_xyy_xxyz[k] = -g_0_x_xy_xxyz[k] * ab_y + g_0_x_xy_xxyyz[k];

                g_0_x_xyy_xxzz[k] = -g_0_x_xy_xxzz[k] * ab_y + g_0_x_xy_xxyzz[k];

                g_0_x_xyy_xyyy[k] = -g_0_x_xy_xyyy[k] * ab_y + g_0_x_xy_xyyyy[k];

                g_0_x_xyy_xyyz[k] = -g_0_x_xy_xyyz[k] * ab_y + g_0_x_xy_xyyyz[k];

                g_0_x_xyy_xyzz[k] = -g_0_x_xy_xyzz[k] * ab_y + g_0_x_xy_xyyzz[k];

                g_0_x_xyy_xzzz[k] = -g_0_x_xy_xzzz[k] * ab_y + g_0_x_xy_xyzzz[k];

                g_0_x_xyy_yyyy[k] = -g_0_x_xy_yyyy[k] * ab_y + g_0_x_xy_yyyyy[k];

                g_0_x_xyy_yyyz[k] = -g_0_x_xy_yyyz[k] * ab_y + g_0_x_xy_yyyyz[k];

                g_0_x_xyy_yyzz[k] = -g_0_x_xy_yyzz[k] * ab_y + g_0_x_xy_yyyzz[k];

                g_0_x_xyy_yzzz[k] = -g_0_x_xy_yzzz[k] * ab_y + g_0_x_xy_yyzzz[k];

                g_0_x_xyy_zzzz[k] = -g_0_x_xy_zzzz[k] * ab_y + g_0_x_xy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyz_xxxx = cbuffer.data(fg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyz_xxxy = cbuffer.data(fg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyz_xxxz = cbuffer.data(fg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyz_xxyy = cbuffer.data(fg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyz_xxyz = cbuffer.data(fg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyz_xxzz = cbuffer.data(fg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyz_xyyy = cbuffer.data(fg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyz_xyyz = cbuffer.data(fg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyz_xyzz = cbuffer.data(fg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyz_xzzz = cbuffer.data(fg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyz_yyyy = cbuffer.data(fg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyz_yyyz = cbuffer.data(fg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyz_yyzz = cbuffer.data(fg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyz_yzzz = cbuffer.data(fg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyz_zzzz = cbuffer.data(fg_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyz_xxxx, g_0_x_xyz_xxxy, g_0_x_xyz_xxxz, g_0_x_xyz_xxyy, g_0_x_xyz_xxyz, g_0_x_xyz_xxzz, g_0_x_xyz_xyyy, g_0_x_xyz_xyyz, g_0_x_xyz_xyzz, g_0_x_xyz_xzzz, g_0_x_xyz_yyyy, g_0_x_xyz_yyyz, g_0_x_xyz_yyzz, g_0_x_xyz_yzzz, g_0_x_xyz_zzzz, g_0_x_xz_xxxx, g_0_x_xz_xxxxy, g_0_x_xz_xxxy, g_0_x_xz_xxxyy, g_0_x_xz_xxxyz, g_0_x_xz_xxxz, g_0_x_xz_xxyy, g_0_x_xz_xxyyy, g_0_x_xz_xxyyz, g_0_x_xz_xxyz, g_0_x_xz_xxyzz, g_0_x_xz_xxzz, g_0_x_xz_xyyy, g_0_x_xz_xyyyy, g_0_x_xz_xyyyz, g_0_x_xz_xyyz, g_0_x_xz_xyyzz, g_0_x_xz_xyzz, g_0_x_xz_xyzzz, g_0_x_xz_xzzz, g_0_x_xz_yyyy, g_0_x_xz_yyyyy, g_0_x_xz_yyyyz, g_0_x_xz_yyyz, g_0_x_xz_yyyzz, g_0_x_xz_yyzz, g_0_x_xz_yyzzz, g_0_x_xz_yzzz, g_0_x_xz_yzzzz, g_0_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyz_xxxx[k] = -g_0_x_xz_xxxx[k] * ab_y + g_0_x_xz_xxxxy[k];

                g_0_x_xyz_xxxy[k] = -g_0_x_xz_xxxy[k] * ab_y + g_0_x_xz_xxxyy[k];

                g_0_x_xyz_xxxz[k] = -g_0_x_xz_xxxz[k] * ab_y + g_0_x_xz_xxxyz[k];

                g_0_x_xyz_xxyy[k] = -g_0_x_xz_xxyy[k] * ab_y + g_0_x_xz_xxyyy[k];

                g_0_x_xyz_xxyz[k] = -g_0_x_xz_xxyz[k] * ab_y + g_0_x_xz_xxyyz[k];

                g_0_x_xyz_xxzz[k] = -g_0_x_xz_xxzz[k] * ab_y + g_0_x_xz_xxyzz[k];

                g_0_x_xyz_xyyy[k] = -g_0_x_xz_xyyy[k] * ab_y + g_0_x_xz_xyyyy[k];

                g_0_x_xyz_xyyz[k] = -g_0_x_xz_xyyz[k] * ab_y + g_0_x_xz_xyyyz[k];

                g_0_x_xyz_xyzz[k] = -g_0_x_xz_xyzz[k] * ab_y + g_0_x_xz_xyyzz[k];

                g_0_x_xyz_xzzz[k] = -g_0_x_xz_xzzz[k] * ab_y + g_0_x_xz_xyzzz[k];

                g_0_x_xyz_yyyy[k] = -g_0_x_xz_yyyy[k] * ab_y + g_0_x_xz_yyyyy[k];

                g_0_x_xyz_yyyz[k] = -g_0_x_xz_yyyz[k] * ab_y + g_0_x_xz_yyyyz[k];

                g_0_x_xyz_yyzz[k] = -g_0_x_xz_yyzz[k] * ab_y + g_0_x_xz_yyyzz[k];

                g_0_x_xyz_yzzz[k] = -g_0_x_xz_yzzz[k] * ab_y + g_0_x_xz_yyzzz[k];

                g_0_x_xyz_zzzz[k] = -g_0_x_xz_zzzz[k] * ab_y + g_0_x_xz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzz_xxxx = cbuffer.data(fg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xzz_xxxy = cbuffer.data(fg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xzz_xxxz = cbuffer.data(fg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xzz_xxyy = cbuffer.data(fg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xzz_xxyz = cbuffer.data(fg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xzz_xxzz = cbuffer.data(fg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xzz_xyyy = cbuffer.data(fg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xzz_xyyz = cbuffer.data(fg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xzz_xyzz = cbuffer.data(fg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xzz_xzzz = cbuffer.data(fg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzz_yyyy = cbuffer.data(fg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzz_yyyz = cbuffer.data(fg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzz_yyzz = cbuffer.data(fg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzz_yzzz = cbuffer.data(fg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzz_zzzz = cbuffer.data(fg_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xz_xxxx, g_0_x_xz_xxxxz, g_0_x_xz_xxxy, g_0_x_xz_xxxyz, g_0_x_xz_xxxz, g_0_x_xz_xxxzz, g_0_x_xz_xxyy, g_0_x_xz_xxyyz, g_0_x_xz_xxyz, g_0_x_xz_xxyzz, g_0_x_xz_xxzz, g_0_x_xz_xxzzz, g_0_x_xz_xyyy, g_0_x_xz_xyyyz, g_0_x_xz_xyyz, g_0_x_xz_xyyzz, g_0_x_xz_xyzz, g_0_x_xz_xyzzz, g_0_x_xz_xzzz, g_0_x_xz_xzzzz, g_0_x_xz_yyyy, g_0_x_xz_yyyyz, g_0_x_xz_yyyz, g_0_x_xz_yyyzz, g_0_x_xz_yyzz, g_0_x_xz_yyzzz, g_0_x_xz_yzzz, g_0_x_xz_yzzzz, g_0_x_xz_zzzz, g_0_x_xz_zzzzz, g_0_x_xzz_xxxx, g_0_x_xzz_xxxy, g_0_x_xzz_xxxz, g_0_x_xzz_xxyy, g_0_x_xzz_xxyz, g_0_x_xzz_xxzz, g_0_x_xzz_xyyy, g_0_x_xzz_xyyz, g_0_x_xzz_xyzz, g_0_x_xzz_xzzz, g_0_x_xzz_yyyy, g_0_x_xzz_yyyz, g_0_x_xzz_yyzz, g_0_x_xzz_yzzz, g_0_x_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzz_xxxx[k] = -g_0_x_xz_xxxx[k] * ab_z + g_0_x_xz_xxxxz[k];

                g_0_x_xzz_xxxy[k] = -g_0_x_xz_xxxy[k] * ab_z + g_0_x_xz_xxxyz[k];

                g_0_x_xzz_xxxz[k] = -g_0_x_xz_xxxz[k] * ab_z + g_0_x_xz_xxxzz[k];

                g_0_x_xzz_xxyy[k] = -g_0_x_xz_xxyy[k] * ab_z + g_0_x_xz_xxyyz[k];

                g_0_x_xzz_xxyz[k] = -g_0_x_xz_xxyz[k] * ab_z + g_0_x_xz_xxyzz[k];

                g_0_x_xzz_xxzz[k] = -g_0_x_xz_xxzz[k] * ab_z + g_0_x_xz_xxzzz[k];

                g_0_x_xzz_xyyy[k] = -g_0_x_xz_xyyy[k] * ab_z + g_0_x_xz_xyyyz[k];

                g_0_x_xzz_xyyz[k] = -g_0_x_xz_xyyz[k] * ab_z + g_0_x_xz_xyyzz[k];

                g_0_x_xzz_xyzz[k] = -g_0_x_xz_xyzz[k] * ab_z + g_0_x_xz_xyzzz[k];

                g_0_x_xzz_xzzz[k] = -g_0_x_xz_xzzz[k] * ab_z + g_0_x_xz_xzzzz[k];

                g_0_x_xzz_yyyy[k] = -g_0_x_xz_yyyy[k] * ab_z + g_0_x_xz_yyyyz[k];

                g_0_x_xzz_yyyz[k] = -g_0_x_xz_yyyz[k] * ab_z + g_0_x_xz_yyyzz[k];

                g_0_x_xzz_yyzz[k] = -g_0_x_xz_yyzz[k] * ab_z + g_0_x_xz_yyzzz[k];

                g_0_x_xzz_yzzz[k] = -g_0_x_xz_yzzz[k] * ab_z + g_0_x_xz_yzzzz[k];

                g_0_x_xzz_zzzz[k] = -g_0_x_xz_zzzz[k] * ab_z + g_0_x_xz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyy_xxxx = cbuffer.data(fg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyy_xxxy = cbuffer.data(fg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyy_xxxz = cbuffer.data(fg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyy_xxyy = cbuffer.data(fg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyy_xxyz = cbuffer.data(fg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyy_xxzz = cbuffer.data(fg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yyy_xyyy = cbuffer.data(fg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyy_xyyz = cbuffer.data(fg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyy_xyzz = cbuffer.data(fg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyy_xzzz = cbuffer.data(fg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyy_yyyy = cbuffer.data(fg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyy_yyyz = cbuffer.data(fg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyy_yyzz = cbuffer.data(fg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyy_yzzz = cbuffer.data(fg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyy_zzzz = cbuffer.data(fg_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_xxxx, g_0_x_yy_xxxxy, g_0_x_yy_xxxy, g_0_x_yy_xxxyy, g_0_x_yy_xxxyz, g_0_x_yy_xxxz, g_0_x_yy_xxyy, g_0_x_yy_xxyyy, g_0_x_yy_xxyyz, g_0_x_yy_xxyz, g_0_x_yy_xxyzz, g_0_x_yy_xxzz, g_0_x_yy_xyyy, g_0_x_yy_xyyyy, g_0_x_yy_xyyyz, g_0_x_yy_xyyz, g_0_x_yy_xyyzz, g_0_x_yy_xyzz, g_0_x_yy_xyzzz, g_0_x_yy_xzzz, g_0_x_yy_yyyy, g_0_x_yy_yyyyy, g_0_x_yy_yyyyz, g_0_x_yy_yyyz, g_0_x_yy_yyyzz, g_0_x_yy_yyzz, g_0_x_yy_yyzzz, g_0_x_yy_yzzz, g_0_x_yy_yzzzz, g_0_x_yy_zzzz, g_0_x_yyy_xxxx, g_0_x_yyy_xxxy, g_0_x_yyy_xxxz, g_0_x_yyy_xxyy, g_0_x_yyy_xxyz, g_0_x_yyy_xxzz, g_0_x_yyy_xyyy, g_0_x_yyy_xyyz, g_0_x_yyy_xyzz, g_0_x_yyy_xzzz, g_0_x_yyy_yyyy, g_0_x_yyy_yyyz, g_0_x_yyy_yyzz, g_0_x_yyy_yzzz, g_0_x_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyy_xxxx[k] = -g_0_x_yy_xxxx[k] * ab_y + g_0_x_yy_xxxxy[k];

                g_0_x_yyy_xxxy[k] = -g_0_x_yy_xxxy[k] * ab_y + g_0_x_yy_xxxyy[k];

                g_0_x_yyy_xxxz[k] = -g_0_x_yy_xxxz[k] * ab_y + g_0_x_yy_xxxyz[k];

                g_0_x_yyy_xxyy[k] = -g_0_x_yy_xxyy[k] * ab_y + g_0_x_yy_xxyyy[k];

                g_0_x_yyy_xxyz[k] = -g_0_x_yy_xxyz[k] * ab_y + g_0_x_yy_xxyyz[k];

                g_0_x_yyy_xxzz[k] = -g_0_x_yy_xxzz[k] * ab_y + g_0_x_yy_xxyzz[k];

                g_0_x_yyy_xyyy[k] = -g_0_x_yy_xyyy[k] * ab_y + g_0_x_yy_xyyyy[k];

                g_0_x_yyy_xyyz[k] = -g_0_x_yy_xyyz[k] * ab_y + g_0_x_yy_xyyyz[k];

                g_0_x_yyy_xyzz[k] = -g_0_x_yy_xyzz[k] * ab_y + g_0_x_yy_xyyzz[k];

                g_0_x_yyy_xzzz[k] = -g_0_x_yy_xzzz[k] * ab_y + g_0_x_yy_xyzzz[k];

                g_0_x_yyy_yyyy[k] = -g_0_x_yy_yyyy[k] * ab_y + g_0_x_yy_yyyyy[k];

                g_0_x_yyy_yyyz[k] = -g_0_x_yy_yyyz[k] * ab_y + g_0_x_yy_yyyyz[k];

                g_0_x_yyy_yyzz[k] = -g_0_x_yy_yyzz[k] * ab_y + g_0_x_yy_yyyzz[k];

                g_0_x_yyy_yzzz[k] = -g_0_x_yy_yzzz[k] * ab_y + g_0_x_yy_yyzzz[k];

                g_0_x_yyy_zzzz[k] = -g_0_x_yy_zzzz[k] * ab_y + g_0_x_yy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyz_xxxx = cbuffer.data(fg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyz_xxxy = cbuffer.data(fg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyz_xxxz = cbuffer.data(fg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyz_xxyy = cbuffer.data(fg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyz_xxyz = cbuffer.data(fg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyz_xxzz = cbuffer.data(fg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyz_xyyy = cbuffer.data(fg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyz_xyyz = cbuffer.data(fg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyz_xyzz = cbuffer.data(fg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yyz_xzzz = cbuffer.data(fg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yyz_yyyy = cbuffer.data(fg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yyz_yyyz = cbuffer.data(fg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yyz_yyzz = cbuffer.data(fg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yyz_yzzz = cbuffer.data(fg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yyz_zzzz = cbuffer.data(fg_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyz_xxxx, g_0_x_yyz_xxxy, g_0_x_yyz_xxxz, g_0_x_yyz_xxyy, g_0_x_yyz_xxyz, g_0_x_yyz_xxzz, g_0_x_yyz_xyyy, g_0_x_yyz_xyyz, g_0_x_yyz_xyzz, g_0_x_yyz_xzzz, g_0_x_yyz_yyyy, g_0_x_yyz_yyyz, g_0_x_yyz_yyzz, g_0_x_yyz_yzzz, g_0_x_yyz_zzzz, g_0_x_yz_xxxx, g_0_x_yz_xxxxy, g_0_x_yz_xxxy, g_0_x_yz_xxxyy, g_0_x_yz_xxxyz, g_0_x_yz_xxxz, g_0_x_yz_xxyy, g_0_x_yz_xxyyy, g_0_x_yz_xxyyz, g_0_x_yz_xxyz, g_0_x_yz_xxyzz, g_0_x_yz_xxzz, g_0_x_yz_xyyy, g_0_x_yz_xyyyy, g_0_x_yz_xyyyz, g_0_x_yz_xyyz, g_0_x_yz_xyyzz, g_0_x_yz_xyzz, g_0_x_yz_xyzzz, g_0_x_yz_xzzz, g_0_x_yz_yyyy, g_0_x_yz_yyyyy, g_0_x_yz_yyyyz, g_0_x_yz_yyyz, g_0_x_yz_yyyzz, g_0_x_yz_yyzz, g_0_x_yz_yyzzz, g_0_x_yz_yzzz, g_0_x_yz_yzzzz, g_0_x_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyz_xxxx[k] = -g_0_x_yz_xxxx[k] * ab_y + g_0_x_yz_xxxxy[k];

                g_0_x_yyz_xxxy[k] = -g_0_x_yz_xxxy[k] * ab_y + g_0_x_yz_xxxyy[k];

                g_0_x_yyz_xxxz[k] = -g_0_x_yz_xxxz[k] * ab_y + g_0_x_yz_xxxyz[k];

                g_0_x_yyz_xxyy[k] = -g_0_x_yz_xxyy[k] * ab_y + g_0_x_yz_xxyyy[k];

                g_0_x_yyz_xxyz[k] = -g_0_x_yz_xxyz[k] * ab_y + g_0_x_yz_xxyyz[k];

                g_0_x_yyz_xxzz[k] = -g_0_x_yz_xxzz[k] * ab_y + g_0_x_yz_xxyzz[k];

                g_0_x_yyz_xyyy[k] = -g_0_x_yz_xyyy[k] * ab_y + g_0_x_yz_xyyyy[k];

                g_0_x_yyz_xyyz[k] = -g_0_x_yz_xyyz[k] * ab_y + g_0_x_yz_xyyyz[k];

                g_0_x_yyz_xyzz[k] = -g_0_x_yz_xyzz[k] * ab_y + g_0_x_yz_xyyzz[k];

                g_0_x_yyz_xzzz[k] = -g_0_x_yz_xzzz[k] * ab_y + g_0_x_yz_xyzzz[k];

                g_0_x_yyz_yyyy[k] = -g_0_x_yz_yyyy[k] * ab_y + g_0_x_yz_yyyyy[k];

                g_0_x_yyz_yyyz[k] = -g_0_x_yz_yyyz[k] * ab_y + g_0_x_yz_yyyyz[k];

                g_0_x_yyz_yyzz[k] = -g_0_x_yz_yyzz[k] * ab_y + g_0_x_yz_yyyzz[k];

                g_0_x_yyz_yzzz[k] = -g_0_x_yz_yzzz[k] * ab_y + g_0_x_yz_yyzzz[k];

                g_0_x_yyz_zzzz[k] = -g_0_x_yz_zzzz[k] * ab_y + g_0_x_yz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzz_xxxx = cbuffer.data(fg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yzz_xxxy = cbuffer.data(fg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yzz_xxxz = cbuffer.data(fg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yzz_xxyy = cbuffer.data(fg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yzz_xxyz = cbuffer.data(fg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yzz_xxzz = cbuffer.data(fg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yzz_xyyy = cbuffer.data(fg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yzz_xyyz = cbuffer.data(fg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yzz_xyzz = cbuffer.data(fg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yzz_xzzz = cbuffer.data(fg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yzz_yyyy = cbuffer.data(fg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yzz_yyyz = cbuffer.data(fg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yzz_yyzz = cbuffer.data(fg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yzz_yzzz = cbuffer.data(fg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yzz_zzzz = cbuffer.data(fg_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzz_xxxx, g_0_x_yzz_xxxy, g_0_x_yzz_xxxz, g_0_x_yzz_xxyy, g_0_x_yzz_xxyz, g_0_x_yzz_xxzz, g_0_x_yzz_xyyy, g_0_x_yzz_xyyz, g_0_x_yzz_xyzz, g_0_x_yzz_xzzz, g_0_x_yzz_yyyy, g_0_x_yzz_yyyz, g_0_x_yzz_yyzz, g_0_x_yzz_yzzz, g_0_x_yzz_zzzz, g_0_x_zz_xxxx, g_0_x_zz_xxxxy, g_0_x_zz_xxxy, g_0_x_zz_xxxyy, g_0_x_zz_xxxyz, g_0_x_zz_xxxz, g_0_x_zz_xxyy, g_0_x_zz_xxyyy, g_0_x_zz_xxyyz, g_0_x_zz_xxyz, g_0_x_zz_xxyzz, g_0_x_zz_xxzz, g_0_x_zz_xyyy, g_0_x_zz_xyyyy, g_0_x_zz_xyyyz, g_0_x_zz_xyyz, g_0_x_zz_xyyzz, g_0_x_zz_xyzz, g_0_x_zz_xyzzz, g_0_x_zz_xzzz, g_0_x_zz_yyyy, g_0_x_zz_yyyyy, g_0_x_zz_yyyyz, g_0_x_zz_yyyz, g_0_x_zz_yyyzz, g_0_x_zz_yyzz, g_0_x_zz_yyzzz, g_0_x_zz_yzzz, g_0_x_zz_yzzzz, g_0_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzz_xxxx[k] = -g_0_x_zz_xxxx[k] * ab_y + g_0_x_zz_xxxxy[k];

                g_0_x_yzz_xxxy[k] = -g_0_x_zz_xxxy[k] * ab_y + g_0_x_zz_xxxyy[k];

                g_0_x_yzz_xxxz[k] = -g_0_x_zz_xxxz[k] * ab_y + g_0_x_zz_xxxyz[k];

                g_0_x_yzz_xxyy[k] = -g_0_x_zz_xxyy[k] * ab_y + g_0_x_zz_xxyyy[k];

                g_0_x_yzz_xxyz[k] = -g_0_x_zz_xxyz[k] * ab_y + g_0_x_zz_xxyyz[k];

                g_0_x_yzz_xxzz[k] = -g_0_x_zz_xxzz[k] * ab_y + g_0_x_zz_xxyzz[k];

                g_0_x_yzz_xyyy[k] = -g_0_x_zz_xyyy[k] * ab_y + g_0_x_zz_xyyyy[k];

                g_0_x_yzz_xyyz[k] = -g_0_x_zz_xyyz[k] * ab_y + g_0_x_zz_xyyyz[k];

                g_0_x_yzz_xyzz[k] = -g_0_x_zz_xyzz[k] * ab_y + g_0_x_zz_xyyzz[k];

                g_0_x_yzz_xzzz[k] = -g_0_x_zz_xzzz[k] * ab_y + g_0_x_zz_xyzzz[k];

                g_0_x_yzz_yyyy[k] = -g_0_x_zz_yyyy[k] * ab_y + g_0_x_zz_yyyyy[k];

                g_0_x_yzz_yyyz[k] = -g_0_x_zz_yyyz[k] * ab_y + g_0_x_zz_yyyyz[k];

                g_0_x_yzz_yyzz[k] = -g_0_x_zz_yyzz[k] * ab_y + g_0_x_zz_yyyzz[k];

                g_0_x_yzz_yzzz[k] = -g_0_x_zz_yzzz[k] * ab_y + g_0_x_zz_yyzzz[k];

                g_0_x_yzz_zzzz[k] = -g_0_x_zz_zzzz[k] * ab_y + g_0_x_zz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzz_xxxx = cbuffer.data(fg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_zzz_xxxy = cbuffer.data(fg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_zzz_xxxz = cbuffer.data(fg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_zzz_xxyy = cbuffer.data(fg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_zzz_xxyz = cbuffer.data(fg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_zzz_xxzz = cbuffer.data(fg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_zzz_xyyy = cbuffer.data(fg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_zzz_xyyz = cbuffer.data(fg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_zzz_xyzz = cbuffer.data(fg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_zzz_xzzz = cbuffer.data(fg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_zzz_yyyy = cbuffer.data(fg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_zzz_yyyz = cbuffer.data(fg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_zzz_yyzz = cbuffer.data(fg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_zzz_yzzz = cbuffer.data(fg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_zzz_zzzz = cbuffer.data(fg_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_xxxx, g_0_x_zz_xxxxz, g_0_x_zz_xxxy, g_0_x_zz_xxxyz, g_0_x_zz_xxxz, g_0_x_zz_xxxzz, g_0_x_zz_xxyy, g_0_x_zz_xxyyz, g_0_x_zz_xxyz, g_0_x_zz_xxyzz, g_0_x_zz_xxzz, g_0_x_zz_xxzzz, g_0_x_zz_xyyy, g_0_x_zz_xyyyz, g_0_x_zz_xyyz, g_0_x_zz_xyyzz, g_0_x_zz_xyzz, g_0_x_zz_xyzzz, g_0_x_zz_xzzz, g_0_x_zz_xzzzz, g_0_x_zz_yyyy, g_0_x_zz_yyyyz, g_0_x_zz_yyyz, g_0_x_zz_yyyzz, g_0_x_zz_yyzz, g_0_x_zz_yyzzz, g_0_x_zz_yzzz, g_0_x_zz_yzzzz, g_0_x_zz_zzzz, g_0_x_zz_zzzzz, g_0_x_zzz_xxxx, g_0_x_zzz_xxxy, g_0_x_zzz_xxxz, g_0_x_zzz_xxyy, g_0_x_zzz_xxyz, g_0_x_zzz_xxzz, g_0_x_zzz_xyyy, g_0_x_zzz_xyyz, g_0_x_zzz_xyzz, g_0_x_zzz_xzzz, g_0_x_zzz_yyyy, g_0_x_zzz_yyyz, g_0_x_zzz_yyzz, g_0_x_zzz_yzzz, g_0_x_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzz_xxxx[k] = -g_0_x_zz_xxxx[k] * ab_z + g_0_x_zz_xxxxz[k];

                g_0_x_zzz_xxxy[k] = -g_0_x_zz_xxxy[k] * ab_z + g_0_x_zz_xxxyz[k];

                g_0_x_zzz_xxxz[k] = -g_0_x_zz_xxxz[k] * ab_z + g_0_x_zz_xxxzz[k];

                g_0_x_zzz_xxyy[k] = -g_0_x_zz_xxyy[k] * ab_z + g_0_x_zz_xxyyz[k];

                g_0_x_zzz_xxyz[k] = -g_0_x_zz_xxyz[k] * ab_z + g_0_x_zz_xxyzz[k];

                g_0_x_zzz_xxzz[k] = -g_0_x_zz_xxzz[k] * ab_z + g_0_x_zz_xxzzz[k];

                g_0_x_zzz_xyyy[k] = -g_0_x_zz_xyyy[k] * ab_z + g_0_x_zz_xyyyz[k];

                g_0_x_zzz_xyyz[k] = -g_0_x_zz_xyyz[k] * ab_z + g_0_x_zz_xyyzz[k];

                g_0_x_zzz_xyzz[k] = -g_0_x_zz_xyzz[k] * ab_z + g_0_x_zz_xyzzz[k];

                g_0_x_zzz_xzzz[k] = -g_0_x_zz_xzzz[k] * ab_z + g_0_x_zz_xzzzz[k];

                g_0_x_zzz_yyyy[k] = -g_0_x_zz_yyyy[k] * ab_z + g_0_x_zz_yyyyz[k];

                g_0_x_zzz_yyyz[k] = -g_0_x_zz_yyyz[k] * ab_z + g_0_x_zz_yyyzz[k];

                g_0_x_zzz_yyzz[k] = -g_0_x_zz_yyzz[k] * ab_z + g_0_x_zz_yyzzz[k];

                g_0_x_zzz_yzzz[k] = -g_0_x_zz_yzzz[k] * ab_z + g_0_x_zz_yzzzz[k];

                g_0_x_zzz_zzzz[k] = -g_0_x_zz_zzzz[k] * ab_z + g_0_x_zz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxx_xxxx = cbuffer.data(fg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxx_xxxy = cbuffer.data(fg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxx_xxxz = cbuffer.data(fg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxx_xxyy = cbuffer.data(fg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxx_xxyz = cbuffer.data(fg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxx_xxzz = cbuffer.data(fg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxx_xyyy = cbuffer.data(fg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxx_xyyz = cbuffer.data(fg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxx_xyzz = cbuffer.data(fg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxx_xzzz = cbuffer.data(fg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxx_yyyy = cbuffer.data(fg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxx_yyyz = cbuffer.data(fg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxx_yyzz = cbuffer.data(fg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxx_yzzz = cbuffer.data(fg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxx_zzzz = cbuffer.data(fg_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_xxxx, g_0_y_xx_xxxxx, g_0_y_xx_xxxxy, g_0_y_xx_xxxxz, g_0_y_xx_xxxy, g_0_y_xx_xxxyy, g_0_y_xx_xxxyz, g_0_y_xx_xxxz, g_0_y_xx_xxxzz, g_0_y_xx_xxyy, g_0_y_xx_xxyyy, g_0_y_xx_xxyyz, g_0_y_xx_xxyz, g_0_y_xx_xxyzz, g_0_y_xx_xxzz, g_0_y_xx_xxzzz, g_0_y_xx_xyyy, g_0_y_xx_xyyyy, g_0_y_xx_xyyyz, g_0_y_xx_xyyz, g_0_y_xx_xyyzz, g_0_y_xx_xyzz, g_0_y_xx_xyzzz, g_0_y_xx_xzzz, g_0_y_xx_xzzzz, g_0_y_xx_yyyy, g_0_y_xx_yyyz, g_0_y_xx_yyzz, g_0_y_xx_yzzz, g_0_y_xx_zzzz, g_0_y_xxx_xxxx, g_0_y_xxx_xxxy, g_0_y_xxx_xxxz, g_0_y_xxx_xxyy, g_0_y_xxx_xxyz, g_0_y_xxx_xxzz, g_0_y_xxx_xyyy, g_0_y_xxx_xyyz, g_0_y_xxx_xyzz, g_0_y_xxx_xzzz, g_0_y_xxx_yyyy, g_0_y_xxx_yyyz, g_0_y_xxx_yyzz, g_0_y_xxx_yzzz, g_0_y_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxx_xxxx[k] = -g_0_y_xx_xxxx[k] * ab_x + g_0_y_xx_xxxxx[k];

                g_0_y_xxx_xxxy[k] = -g_0_y_xx_xxxy[k] * ab_x + g_0_y_xx_xxxxy[k];

                g_0_y_xxx_xxxz[k] = -g_0_y_xx_xxxz[k] * ab_x + g_0_y_xx_xxxxz[k];

                g_0_y_xxx_xxyy[k] = -g_0_y_xx_xxyy[k] * ab_x + g_0_y_xx_xxxyy[k];

                g_0_y_xxx_xxyz[k] = -g_0_y_xx_xxyz[k] * ab_x + g_0_y_xx_xxxyz[k];

                g_0_y_xxx_xxzz[k] = -g_0_y_xx_xxzz[k] * ab_x + g_0_y_xx_xxxzz[k];

                g_0_y_xxx_xyyy[k] = -g_0_y_xx_xyyy[k] * ab_x + g_0_y_xx_xxyyy[k];

                g_0_y_xxx_xyyz[k] = -g_0_y_xx_xyyz[k] * ab_x + g_0_y_xx_xxyyz[k];

                g_0_y_xxx_xyzz[k] = -g_0_y_xx_xyzz[k] * ab_x + g_0_y_xx_xxyzz[k];

                g_0_y_xxx_xzzz[k] = -g_0_y_xx_xzzz[k] * ab_x + g_0_y_xx_xxzzz[k];

                g_0_y_xxx_yyyy[k] = -g_0_y_xx_yyyy[k] * ab_x + g_0_y_xx_xyyyy[k];

                g_0_y_xxx_yyyz[k] = -g_0_y_xx_yyyz[k] * ab_x + g_0_y_xx_xyyyz[k];

                g_0_y_xxx_yyzz[k] = -g_0_y_xx_yyzz[k] * ab_x + g_0_y_xx_xyyzz[k];

                g_0_y_xxx_yzzz[k] = -g_0_y_xx_yzzz[k] * ab_x + g_0_y_xx_xyzzz[k];

                g_0_y_xxx_zzzz[k] = -g_0_y_xx_zzzz[k] * ab_x + g_0_y_xx_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxy_xxxx = cbuffer.data(fg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxy_xxxy = cbuffer.data(fg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxy_xxxz = cbuffer.data(fg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxy_xxyy = cbuffer.data(fg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxy_xxyz = cbuffer.data(fg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxy_xxzz = cbuffer.data(fg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxy_xyyy = cbuffer.data(fg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxy_xyyz = cbuffer.data(fg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxy_xyzz = cbuffer.data(fg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxy_xzzz = cbuffer.data(fg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxy_yyyy = cbuffer.data(fg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxy_yyyz = cbuffer.data(fg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxy_yyzz = cbuffer.data(fg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxy_yzzz = cbuffer.data(fg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxy_zzzz = cbuffer.data(fg_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxy_xxxx, g_0_y_xxy_xxxy, g_0_y_xxy_xxxz, g_0_y_xxy_xxyy, g_0_y_xxy_xxyz, g_0_y_xxy_xxzz, g_0_y_xxy_xyyy, g_0_y_xxy_xyyz, g_0_y_xxy_xyzz, g_0_y_xxy_xzzz, g_0_y_xxy_yyyy, g_0_y_xxy_yyyz, g_0_y_xxy_yyzz, g_0_y_xxy_yzzz, g_0_y_xxy_zzzz, g_0_y_xy_xxxx, g_0_y_xy_xxxxx, g_0_y_xy_xxxxy, g_0_y_xy_xxxxz, g_0_y_xy_xxxy, g_0_y_xy_xxxyy, g_0_y_xy_xxxyz, g_0_y_xy_xxxz, g_0_y_xy_xxxzz, g_0_y_xy_xxyy, g_0_y_xy_xxyyy, g_0_y_xy_xxyyz, g_0_y_xy_xxyz, g_0_y_xy_xxyzz, g_0_y_xy_xxzz, g_0_y_xy_xxzzz, g_0_y_xy_xyyy, g_0_y_xy_xyyyy, g_0_y_xy_xyyyz, g_0_y_xy_xyyz, g_0_y_xy_xyyzz, g_0_y_xy_xyzz, g_0_y_xy_xyzzz, g_0_y_xy_xzzz, g_0_y_xy_xzzzz, g_0_y_xy_yyyy, g_0_y_xy_yyyz, g_0_y_xy_yyzz, g_0_y_xy_yzzz, g_0_y_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxy_xxxx[k] = -g_0_y_xy_xxxx[k] * ab_x + g_0_y_xy_xxxxx[k];

                g_0_y_xxy_xxxy[k] = -g_0_y_xy_xxxy[k] * ab_x + g_0_y_xy_xxxxy[k];

                g_0_y_xxy_xxxz[k] = -g_0_y_xy_xxxz[k] * ab_x + g_0_y_xy_xxxxz[k];

                g_0_y_xxy_xxyy[k] = -g_0_y_xy_xxyy[k] * ab_x + g_0_y_xy_xxxyy[k];

                g_0_y_xxy_xxyz[k] = -g_0_y_xy_xxyz[k] * ab_x + g_0_y_xy_xxxyz[k];

                g_0_y_xxy_xxzz[k] = -g_0_y_xy_xxzz[k] * ab_x + g_0_y_xy_xxxzz[k];

                g_0_y_xxy_xyyy[k] = -g_0_y_xy_xyyy[k] * ab_x + g_0_y_xy_xxyyy[k];

                g_0_y_xxy_xyyz[k] = -g_0_y_xy_xyyz[k] * ab_x + g_0_y_xy_xxyyz[k];

                g_0_y_xxy_xyzz[k] = -g_0_y_xy_xyzz[k] * ab_x + g_0_y_xy_xxyzz[k];

                g_0_y_xxy_xzzz[k] = -g_0_y_xy_xzzz[k] * ab_x + g_0_y_xy_xxzzz[k];

                g_0_y_xxy_yyyy[k] = -g_0_y_xy_yyyy[k] * ab_x + g_0_y_xy_xyyyy[k];

                g_0_y_xxy_yyyz[k] = -g_0_y_xy_yyyz[k] * ab_x + g_0_y_xy_xyyyz[k];

                g_0_y_xxy_yyzz[k] = -g_0_y_xy_yyzz[k] * ab_x + g_0_y_xy_xyyzz[k];

                g_0_y_xxy_yzzz[k] = -g_0_y_xy_yzzz[k] * ab_x + g_0_y_xy_xyzzz[k];

                g_0_y_xxy_zzzz[k] = -g_0_y_xy_zzzz[k] * ab_x + g_0_y_xy_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxz_xxxx = cbuffer.data(fg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxz_xxxy = cbuffer.data(fg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxz_xxxz = cbuffer.data(fg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxz_xxyy = cbuffer.data(fg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxz_xxyz = cbuffer.data(fg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxz_xxzz = cbuffer.data(fg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxz_xyyy = cbuffer.data(fg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxz_xyyz = cbuffer.data(fg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxz_xyzz = cbuffer.data(fg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxz_xzzz = cbuffer.data(fg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxz_yyyy = cbuffer.data(fg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxz_yyyz = cbuffer.data(fg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxz_yyzz = cbuffer.data(fg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxz_yzzz = cbuffer.data(fg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxz_zzzz = cbuffer.data(fg_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxz_xxxx, g_0_y_xxz_xxxy, g_0_y_xxz_xxxz, g_0_y_xxz_xxyy, g_0_y_xxz_xxyz, g_0_y_xxz_xxzz, g_0_y_xxz_xyyy, g_0_y_xxz_xyyz, g_0_y_xxz_xyzz, g_0_y_xxz_xzzz, g_0_y_xxz_yyyy, g_0_y_xxz_yyyz, g_0_y_xxz_yyzz, g_0_y_xxz_yzzz, g_0_y_xxz_zzzz, g_0_y_xz_xxxx, g_0_y_xz_xxxxx, g_0_y_xz_xxxxy, g_0_y_xz_xxxxz, g_0_y_xz_xxxy, g_0_y_xz_xxxyy, g_0_y_xz_xxxyz, g_0_y_xz_xxxz, g_0_y_xz_xxxzz, g_0_y_xz_xxyy, g_0_y_xz_xxyyy, g_0_y_xz_xxyyz, g_0_y_xz_xxyz, g_0_y_xz_xxyzz, g_0_y_xz_xxzz, g_0_y_xz_xxzzz, g_0_y_xz_xyyy, g_0_y_xz_xyyyy, g_0_y_xz_xyyyz, g_0_y_xz_xyyz, g_0_y_xz_xyyzz, g_0_y_xz_xyzz, g_0_y_xz_xyzzz, g_0_y_xz_xzzz, g_0_y_xz_xzzzz, g_0_y_xz_yyyy, g_0_y_xz_yyyz, g_0_y_xz_yyzz, g_0_y_xz_yzzz, g_0_y_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxz_xxxx[k] = -g_0_y_xz_xxxx[k] * ab_x + g_0_y_xz_xxxxx[k];

                g_0_y_xxz_xxxy[k] = -g_0_y_xz_xxxy[k] * ab_x + g_0_y_xz_xxxxy[k];

                g_0_y_xxz_xxxz[k] = -g_0_y_xz_xxxz[k] * ab_x + g_0_y_xz_xxxxz[k];

                g_0_y_xxz_xxyy[k] = -g_0_y_xz_xxyy[k] * ab_x + g_0_y_xz_xxxyy[k];

                g_0_y_xxz_xxyz[k] = -g_0_y_xz_xxyz[k] * ab_x + g_0_y_xz_xxxyz[k];

                g_0_y_xxz_xxzz[k] = -g_0_y_xz_xxzz[k] * ab_x + g_0_y_xz_xxxzz[k];

                g_0_y_xxz_xyyy[k] = -g_0_y_xz_xyyy[k] * ab_x + g_0_y_xz_xxyyy[k];

                g_0_y_xxz_xyyz[k] = -g_0_y_xz_xyyz[k] * ab_x + g_0_y_xz_xxyyz[k];

                g_0_y_xxz_xyzz[k] = -g_0_y_xz_xyzz[k] * ab_x + g_0_y_xz_xxyzz[k];

                g_0_y_xxz_xzzz[k] = -g_0_y_xz_xzzz[k] * ab_x + g_0_y_xz_xxzzz[k];

                g_0_y_xxz_yyyy[k] = -g_0_y_xz_yyyy[k] * ab_x + g_0_y_xz_xyyyy[k];

                g_0_y_xxz_yyyz[k] = -g_0_y_xz_yyyz[k] * ab_x + g_0_y_xz_xyyyz[k];

                g_0_y_xxz_yyzz[k] = -g_0_y_xz_yyzz[k] * ab_x + g_0_y_xz_xyyzz[k];

                g_0_y_xxz_yzzz[k] = -g_0_y_xz_yzzz[k] * ab_x + g_0_y_xz_xyzzz[k];

                g_0_y_xxz_zzzz[k] = -g_0_y_xz_zzzz[k] * ab_x + g_0_y_xz_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyy_xxxx = cbuffer.data(fg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyy_xxxy = cbuffer.data(fg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyy_xxxz = cbuffer.data(fg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xyy_xxyy = cbuffer.data(fg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyy_xxyz = cbuffer.data(fg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyy_xxzz = cbuffer.data(fg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyy_xyyy = cbuffer.data(fg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyy_xyyz = cbuffer.data(fg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyy_xyzz = cbuffer.data(fg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xyy_xzzz = cbuffer.data(fg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyy_yyyy = cbuffer.data(fg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyy_yyyz = cbuffer.data(fg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyy_yyzz = cbuffer.data(fg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyy_yzzz = cbuffer.data(fg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyy_zzzz = cbuffer.data(fg_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyy_xxxx, g_0_y_xyy_xxxy, g_0_y_xyy_xxxz, g_0_y_xyy_xxyy, g_0_y_xyy_xxyz, g_0_y_xyy_xxzz, g_0_y_xyy_xyyy, g_0_y_xyy_xyyz, g_0_y_xyy_xyzz, g_0_y_xyy_xzzz, g_0_y_xyy_yyyy, g_0_y_xyy_yyyz, g_0_y_xyy_yyzz, g_0_y_xyy_yzzz, g_0_y_xyy_zzzz, g_0_y_yy_xxxx, g_0_y_yy_xxxxx, g_0_y_yy_xxxxy, g_0_y_yy_xxxxz, g_0_y_yy_xxxy, g_0_y_yy_xxxyy, g_0_y_yy_xxxyz, g_0_y_yy_xxxz, g_0_y_yy_xxxzz, g_0_y_yy_xxyy, g_0_y_yy_xxyyy, g_0_y_yy_xxyyz, g_0_y_yy_xxyz, g_0_y_yy_xxyzz, g_0_y_yy_xxzz, g_0_y_yy_xxzzz, g_0_y_yy_xyyy, g_0_y_yy_xyyyy, g_0_y_yy_xyyyz, g_0_y_yy_xyyz, g_0_y_yy_xyyzz, g_0_y_yy_xyzz, g_0_y_yy_xyzzz, g_0_y_yy_xzzz, g_0_y_yy_xzzzz, g_0_y_yy_yyyy, g_0_y_yy_yyyz, g_0_y_yy_yyzz, g_0_y_yy_yzzz, g_0_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyy_xxxx[k] = -g_0_y_yy_xxxx[k] * ab_x + g_0_y_yy_xxxxx[k];

                g_0_y_xyy_xxxy[k] = -g_0_y_yy_xxxy[k] * ab_x + g_0_y_yy_xxxxy[k];

                g_0_y_xyy_xxxz[k] = -g_0_y_yy_xxxz[k] * ab_x + g_0_y_yy_xxxxz[k];

                g_0_y_xyy_xxyy[k] = -g_0_y_yy_xxyy[k] * ab_x + g_0_y_yy_xxxyy[k];

                g_0_y_xyy_xxyz[k] = -g_0_y_yy_xxyz[k] * ab_x + g_0_y_yy_xxxyz[k];

                g_0_y_xyy_xxzz[k] = -g_0_y_yy_xxzz[k] * ab_x + g_0_y_yy_xxxzz[k];

                g_0_y_xyy_xyyy[k] = -g_0_y_yy_xyyy[k] * ab_x + g_0_y_yy_xxyyy[k];

                g_0_y_xyy_xyyz[k] = -g_0_y_yy_xyyz[k] * ab_x + g_0_y_yy_xxyyz[k];

                g_0_y_xyy_xyzz[k] = -g_0_y_yy_xyzz[k] * ab_x + g_0_y_yy_xxyzz[k];

                g_0_y_xyy_xzzz[k] = -g_0_y_yy_xzzz[k] * ab_x + g_0_y_yy_xxzzz[k];

                g_0_y_xyy_yyyy[k] = -g_0_y_yy_yyyy[k] * ab_x + g_0_y_yy_xyyyy[k];

                g_0_y_xyy_yyyz[k] = -g_0_y_yy_yyyz[k] * ab_x + g_0_y_yy_xyyyz[k];

                g_0_y_xyy_yyzz[k] = -g_0_y_yy_yyzz[k] * ab_x + g_0_y_yy_xyyzz[k];

                g_0_y_xyy_yzzz[k] = -g_0_y_yy_yzzz[k] * ab_x + g_0_y_yy_xyzzz[k];

                g_0_y_xyy_zzzz[k] = -g_0_y_yy_zzzz[k] * ab_x + g_0_y_yy_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyz_xxxx = cbuffer.data(fg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xyz_xxxy = cbuffer.data(fg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xyz_xxxz = cbuffer.data(fg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xyz_xxyy = cbuffer.data(fg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xyz_xxyz = cbuffer.data(fg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xyz_xxzz = cbuffer.data(fg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xyz_xyyy = cbuffer.data(fg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xyz_xyyz = cbuffer.data(fg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xyz_xyzz = cbuffer.data(fg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xyz_xzzz = cbuffer.data(fg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xyz_yyyy = cbuffer.data(fg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xyz_yyyz = cbuffer.data(fg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xyz_yyzz = cbuffer.data(fg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xyz_yzzz = cbuffer.data(fg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xyz_zzzz = cbuffer.data(fg_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyz_xxxx, g_0_y_xyz_xxxy, g_0_y_xyz_xxxz, g_0_y_xyz_xxyy, g_0_y_xyz_xxyz, g_0_y_xyz_xxzz, g_0_y_xyz_xyyy, g_0_y_xyz_xyyz, g_0_y_xyz_xyzz, g_0_y_xyz_xzzz, g_0_y_xyz_yyyy, g_0_y_xyz_yyyz, g_0_y_xyz_yyzz, g_0_y_xyz_yzzz, g_0_y_xyz_zzzz, g_0_y_yz_xxxx, g_0_y_yz_xxxxx, g_0_y_yz_xxxxy, g_0_y_yz_xxxxz, g_0_y_yz_xxxy, g_0_y_yz_xxxyy, g_0_y_yz_xxxyz, g_0_y_yz_xxxz, g_0_y_yz_xxxzz, g_0_y_yz_xxyy, g_0_y_yz_xxyyy, g_0_y_yz_xxyyz, g_0_y_yz_xxyz, g_0_y_yz_xxyzz, g_0_y_yz_xxzz, g_0_y_yz_xxzzz, g_0_y_yz_xyyy, g_0_y_yz_xyyyy, g_0_y_yz_xyyyz, g_0_y_yz_xyyz, g_0_y_yz_xyyzz, g_0_y_yz_xyzz, g_0_y_yz_xyzzz, g_0_y_yz_xzzz, g_0_y_yz_xzzzz, g_0_y_yz_yyyy, g_0_y_yz_yyyz, g_0_y_yz_yyzz, g_0_y_yz_yzzz, g_0_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyz_xxxx[k] = -g_0_y_yz_xxxx[k] * ab_x + g_0_y_yz_xxxxx[k];

                g_0_y_xyz_xxxy[k] = -g_0_y_yz_xxxy[k] * ab_x + g_0_y_yz_xxxxy[k];

                g_0_y_xyz_xxxz[k] = -g_0_y_yz_xxxz[k] * ab_x + g_0_y_yz_xxxxz[k];

                g_0_y_xyz_xxyy[k] = -g_0_y_yz_xxyy[k] * ab_x + g_0_y_yz_xxxyy[k];

                g_0_y_xyz_xxyz[k] = -g_0_y_yz_xxyz[k] * ab_x + g_0_y_yz_xxxyz[k];

                g_0_y_xyz_xxzz[k] = -g_0_y_yz_xxzz[k] * ab_x + g_0_y_yz_xxxzz[k];

                g_0_y_xyz_xyyy[k] = -g_0_y_yz_xyyy[k] * ab_x + g_0_y_yz_xxyyy[k];

                g_0_y_xyz_xyyz[k] = -g_0_y_yz_xyyz[k] * ab_x + g_0_y_yz_xxyyz[k];

                g_0_y_xyz_xyzz[k] = -g_0_y_yz_xyzz[k] * ab_x + g_0_y_yz_xxyzz[k];

                g_0_y_xyz_xzzz[k] = -g_0_y_yz_xzzz[k] * ab_x + g_0_y_yz_xxzzz[k];

                g_0_y_xyz_yyyy[k] = -g_0_y_yz_yyyy[k] * ab_x + g_0_y_yz_xyyyy[k];

                g_0_y_xyz_yyyz[k] = -g_0_y_yz_yyyz[k] * ab_x + g_0_y_yz_xyyyz[k];

                g_0_y_xyz_yyzz[k] = -g_0_y_yz_yyzz[k] * ab_x + g_0_y_yz_xyyzz[k];

                g_0_y_xyz_yzzz[k] = -g_0_y_yz_yzzz[k] * ab_x + g_0_y_yz_xyzzz[k];

                g_0_y_xyz_zzzz[k] = -g_0_y_yz_zzzz[k] * ab_x + g_0_y_yz_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzz_xxxx = cbuffer.data(fg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xzz_xxxy = cbuffer.data(fg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xzz_xxxz = cbuffer.data(fg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xzz_xxyy = cbuffer.data(fg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xzz_xxyz = cbuffer.data(fg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xzz_xxzz = cbuffer.data(fg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xzz_xyyy = cbuffer.data(fg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xzz_xyyz = cbuffer.data(fg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xzz_xyzz = cbuffer.data(fg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xzz_xzzz = cbuffer.data(fg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xzz_yyyy = cbuffer.data(fg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xzz_yyyz = cbuffer.data(fg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xzz_yyzz = cbuffer.data(fg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xzz_yzzz = cbuffer.data(fg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xzz_zzzz = cbuffer.data(fg_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzz_xxxx, g_0_y_xzz_xxxy, g_0_y_xzz_xxxz, g_0_y_xzz_xxyy, g_0_y_xzz_xxyz, g_0_y_xzz_xxzz, g_0_y_xzz_xyyy, g_0_y_xzz_xyyz, g_0_y_xzz_xyzz, g_0_y_xzz_xzzz, g_0_y_xzz_yyyy, g_0_y_xzz_yyyz, g_0_y_xzz_yyzz, g_0_y_xzz_yzzz, g_0_y_xzz_zzzz, g_0_y_zz_xxxx, g_0_y_zz_xxxxx, g_0_y_zz_xxxxy, g_0_y_zz_xxxxz, g_0_y_zz_xxxy, g_0_y_zz_xxxyy, g_0_y_zz_xxxyz, g_0_y_zz_xxxz, g_0_y_zz_xxxzz, g_0_y_zz_xxyy, g_0_y_zz_xxyyy, g_0_y_zz_xxyyz, g_0_y_zz_xxyz, g_0_y_zz_xxyzz, g_0_y_zz_xxzz, g_0_y_zz_xxzzz, g_0_y_zz_xyyy, g_0_y_zz_xyyyy, g_0_y_zz_xyyyz, g_0_y_zz_xyyz, g_0_y_zz_xyyzz, g_0_y_zz_xyzz, g_0_y_zz_xyzzz, g_0_y_zz_xzzz, g_0_y_zz_xzzzz, g_0_y_zz_yyyy, g_0_y_zz_yyyz, g_0_y_zz_yyzz, g_0_y_zz_yzzz, g_0_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzz_xxxx[k] = -g_0_y_zz_xxxx[k] * ab_x + g_0_y_zz_xxxxx[k];

                g_0_y_xzz_xxxy[k] = -g_0_y_zz_xxxy[k] * ab_x + g_0_y_zz_xxxxy[k];

                g_0_y_xzz_xxxz[k] = -g_0_y_zz_xxxz[k] * ab_x + g_0_y_zz_xxxxz[k];

                g_0_y_xzz_xxyy[k] = -g_0_y_zz_xxyy[k] * ab_x + g_0_y_zz_xxxyy[k];

                g_0_y_xzz_xxyz[k] = -g_0_y_zz_xxyz[k] * ab_x + g_0_y_zz_xxxyz[k];

                g_0_y_xzz_xxzz[k] = -g_0_y_zz_xxzz[k] * ab_x + g_0_y_zz_xxxzz[k];

                g_0_y_xzz_xyyy[k] = -g_0_y_zz_xyyy[k] * ab_x + g_0_y_zz_xxyyy[k];

                g_0_y_xzz_xyyz[k] = -g_0_y_zz_xyyz[k] * ab_x + g_0_y_zz_xxyyz[k];

                g_0_y_xzz_xyzz[k] = -g_0_y_zz_xyzz[k] * ab_x + g_0_y_zz_xxyzz[k];

                g_0_y_xzz_xzzz[k] = -g_0_y_zz_xzzz[k] * ab_x + g_0_y_zz_xxzzz[k];

                g_0_y_xzz_yyyy[k] = -g_0_y_zz_yyyy[k] * ab_x + g_0_y_zz_xyyyy[k];

                g_0_y_xzz_yyyz[k] = -g_0_y_zz_yyyz[k] * ab_x + g_0_y_zz_xyyyz[k];

                g_0_y_xzz_yyzz[k] = -g_0_y_zz_yyzz[k] * ab_x + g_0_y_zz_xyyzz[k];

                g_0_y_xzz_yzzz[k] = -g_0_y_zz_yzzz[k] * ab_x + g_0_y_zz_xyzzz[k];

                g_0_y_xzz_zzzz[k] = -g_0_y_zz_zzzz[k] * ab_x + g_0_y_zz_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyy_xxxx = cbuffer.data(fg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yyy_xxxy = cbuffer.data(fg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yyy_xxxz = cbuffer.data(fg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yyy_xxyy = cbuffer.data(fg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yyy_xxyz = cbuffer.data(fg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yyy_xxzz = cbuffer.data(fg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_yyy_xyyy = cbuffer.data(fg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_yyy_xyyz = cbuffer.data(fg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_yyy_xyzz = cbuffer.data(fg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_yyy_xzzz = cbuffer.data(fg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_yyy_yyyy = cbuffer.data(fg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_yyy_yyyz = cbuffer.data(fg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_yyy_yyzz = cbuffer.data(fg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_yyy_yzzz = cbuffer.data(fg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_yyy_zzzz = cbuffer.data(fg_geom_01_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxxx, g_0_y_yy_xxxxy, g_0_y_yy_xxxy, g_0_y_yy_xxxyy, g_0_y_yy_xxxyz, g_0_y_yy_xxxz, g_0_y_yy_xxyy, g_0_y_yy_xxyyy, g_0_y_yy_xxyyz, g_0_y_yy_xxyz, g_0_y_yy_xxyzz, g_0_y_yy_xxzz, g_0_y_yy_xyyy, g_0_y_yy_xyyyy, g_0_y_yy_xyyyz, g_0_y_yy_xyyz, g_0_y_yy_xyyzz, g_0_y_yy_xyzz, g_0_y_yy_xyzzz, g_0_y_yy_xzzz, g_0_y_yy_yyyy, g_0_y_yy_yyyyy, g_0_y_yy_yyyyz, g_0_y_yy_yyyz, g_0_y_yy_yyyzz, g_0_y_yy_yyzz, g_0_y_yy_yyzzz, g_0_y_yy_yzzz, g_0_y_yy_yzzzz, g_0_y_yy_zzzz, g_0_y_yyy_xxxx, g_0_y_yyy_xxxy, g_0_y_yyy_xxxz, g_0_y_yyy_xxyy, g_0_y_yyy_xxyz, g_0_y_yyy_xxzz, g_0_y_yyy_xyyy, g_0_y_yyy_xyyz, g_0_y_yyy_xyzz, g_0_y_yyy_xzzz, g_0_y_yyy_yyyy, g_0_y_yyy_yyyz, g_0_y_yyy_yyzz, g_0_y_yyy_yzzz, g_0_y_yyy_zzzz, g_yy_xxxx, g_yy_xxxy, g_yy_xxxz, g_yy_xxyy, g_yy_xxyz, g_yy_xxzz, g_yy_xyyy, g_yy_xyyz, g_yy_xyzz, g_yy_xzzz, g_yy_yyyy, g_yy_yyyz, g_yy_yyzz, g_yy_yzzz, g_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyy_xxxx[k] = g_yy_xxxx[k] - g_0_y_yy_xxxx[k] * ab_y + g_0_y_yy_xxxxy[k];

                g_0_y_yyy_xxxy[k] = g_yy_xxxy[k] - g_0_y_yy_xxxy[k] * ab_y + g_0_y_yy_xxxyy[k];

                g_0_y_yyy_xxxz[k] = g_yy_xxxz[k] - g_0_y_yy_xxxz[k] * ab_y + g_0_y_yy_xxxyz[k];

                g_0_y_yyy_xxyy[k] = g_yy_xxyy[k] - g_0_y_yy_xxyy[k] * ab_y + g_0_y_yy_xxyyy[k];

                g_0_y_yyy_xxyz[k] = g_yy_xxyz[k] - g_0_y_yy_xxyz[k] * ab_y + g_0_y_yy_xxyyz[k];

                g_0_y_yyy_xxzz[k] = g_yy_xxzz[k] - g_0_y_yy_xxzz[k] * ab_y + g_0_y_yy_xxyzz[k];

                g_0_y_yyy_xyyy[k] = g_yy_xyyy[k] - g_0_y_yy_xyyy[k] * ab_y + g_0_y_yy_xyyyy[k];

                g_0_y_yyy_xyyz[k] = g_yy_xyyz[k] - g_0_y_yy_xyyz[k] * ab_y + g_0_y_yy_xyyyz[k];

                g_0_y_yyy_xyzz[k] = g_yy_xyzz[k] - g_0_y_yy_xyzz[k] * ab_y + g_0_y_yy_xyyzz[k];

                g_0_y_yyy_xzzz[k] = g_yy_xzzz[k] - g_0_y_yy_xzzz[k] * ab_y + g_0_y_yy_xyzzz[k];

                g_0_y_yyy_yyyy[k] = g_yy_yyyy[k] - g_0_y_yy_yyyy[k] * ab_y + g_0_y_yy_yyyyy[k];

                g_0_y_yyy_yyyz[k] = g_yy_yyyz[k] - g_0_y_yy_yyyz[k] * ab_y + g_0_y_yy_yyyyz[k];

                g_0_y_yyy_yyzz[k] = g_yy_yyzz[k] - g_0_y_yy_yyzz[k] * ab_y + g_0_y_yy_yyyzz[k];

                g_0_y_yyy_yzzz[k] = g_yy_yzzz[k] - g_0_y_yy_yzzz[k] * ab_y + g_0_y_yy_yyzzz[k];

                g_0_y_yyy_zzzz[k] = g_yy_zzzz[k] - g_0_y_yy_zzzz[k] * ab_y + g_0_y_yy_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyz_xxxx = cbuffer.data(fg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_yyz_xxxy = cbuffer.data(fg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_yyz_xxxz = cbuffer.data(fg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_yyz_xxyy = cbuffer.data(fg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_yyz_xxyz = cbuffer.data(fg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_yyz_xxzz = cbuffer.data(fg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_yyz_xyyy = cbuffer.data(fg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_yyz_xyyz = cbuffer.data(fg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_yyz_xyzz = cbuffer.data(fg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_yyz_xzzz = cbuffer.data(fg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_yyz_yyyy = cbuffer.data(fg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_yyz_yyyz = cbuffer.data(fg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_yyz_yyzz = cbuffer.data(fg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_yyz_yzzz = cbuffer.data(fg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_yyz_zzzz = cbuffer.data(fg_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxxx, g_0_y_yy_xxxxz, g_0_y_yy_xxxy, g_0_y_yy_xxxyz, g_0_y_yy_xxxz, g_0_y_yy_xxxzz, g_0_y_yy_xxyy, g_0_y_yy_xxyyz, g_0_y_yy_xxyz, g_0_y_yy_xxyzz, g_0_y_yy_xxzz, g_0_y_yy_xxzzz, g_0_y_yy_xyyy, g_0_y_yy_xyyyz, g_0_y_yy_xyyz, g_0_y_yy_xyyzz, g_0_y_yy_xyzz, g_0_y_yy_xyzzz, g_0_y_yy_xzzz, g_0_y_yy_xzzzz, g_0_y_yy_yyyy, g_0_y_yy_yyyyz, g_0_y_yy_yyyz, g_0_y_yy_yyyzz, g_0_y_yy_yyzz, g_0_y_yy_yyzzz, g_0_y_yy_yzzz, g_0_y_yy_yzzzz, g_0_y_yy_zzzz, g_0_y_yy_zzzzz, g_0_y_yyz_xxxx, g_0_y_yyz_xxxy, g_0_y_yyz_xxxz, g_0_y_yyz_xxyy, g_0_y_yyz_xxyz, g_0_y_yyz_xxzz, g_0_y_yyz_xyyy, g_0_y_yyz_xyyz, g_0_y_yyz_xyzz, g_0_y_yyz_xzzz, g_0_y_yyz_yyyy, g_0_y_yyz_yyyz, g_0_y_yyz_yyzz, g_0_y_yyz_yzzz, g_0_y_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyz_xxxx[k] = -g_0_y_yy_xxxx[k] * ab_z + g_0_y_yy_xxxxz[k];

                g_0_y_yyz_xxxy[k] = -g_0_y_yy_xxxy[k] * ab_z + g_0_y_yy_xxxyz[k];

                g_0_y_yyz_xxxz[k] = -g_0_y_yy_xxxz[k] * ab_z + g_0_y_yy_xxxzz[k];

                g_0_y_yyz_xxyy[k] = -g_0_y_yy_xxyy[k] * ab_z + g_0_y_yy_xxyyz[k];

                g_0_y_yyz_xxyz[k] = -g_0_y_yy_xxyz[k] * ab_z + g_0_y_yy_xxyzz[k];

                g_0_y_yyz_xxzz[k] = -g_0_y_yy_xxzz[k] * ab_z + g_0_y_yy_xxzzz[k];

                g_0_y_yyz_xyyy[k] = -g_0_y_yy_xyyy[k] * ab_z + g_0_y_yy_xyyyz[k];

                g_0_y_yyz_xyyz[k] = -g_0_y_yy_xyyz[k] * ab_z + g_0_y_yy_xyyzz[k];

                g_0_y_yyz_xyzz[k] = -g_0_y_yy_xyzz[k] * ab_z + g_0_y_yy_xyzzz[k];

                g_0_y_yyz_xzzz[k] = -g_0_y_yy_xzzz[k] * ab_z + g_0_y_yy_xzzzz[k];

                g_0_y_yyz_yyyy[k] = -g_0_y_yy_yyyy[k] * ab_z + g_0_y_yy_yyyyz[k];

                g_0_y_yyz_yyyz[k] = -g_0_y_yy_yyyz[k] * ab_z + g_0_y_yy_yyyzz[k];

                g_0_y_yyz_yyzz[k] = -g_0_y_yy_yyzz[k] * ab_z + g_0_y_yy_yyzzz[k];

                g_0_y_yyz_yzzz[k] = -g_0_y_yy_yzzz[k] * ab_z + g_0_y_yy_yzzzz[k];

                g_0_y_yyz_zzzz[k] = -g_0_y_yy_zzzz[k] * ab_z + g_0_y_yy_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzz_xxxx = cbuffer.data(fg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_yzz_xxxy = cbuffer.data(fg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_yzz_xxxz = cbuffer.data(fg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_yzz_xxyy = cbuffer.data(fg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_yzz_xxyz = cbuffer.data(fg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_yzz_xxzz = cbuffer.data(fg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_yzz_xyyy = cbuffer.data(fg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_yzz_xyyz = cbuffer.data(fg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_yzz_xyzz = cbuffer.data(fg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_yzz_xzzz = cbuffer.data(fg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_yzz_yyyy = cbuffer.data(fg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_yzz_yyyz = cbuffer.data(fg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_yzz_yyzz = cbuffer.data(fg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_yzz_yzzz = cbuffer.data(fg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_yzz_zzzz = cbuffer.data(fg_geom_01_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yz_xxxx, g_0_y_yz_xxxxz, g_0_y_yz_xxxy, g_0_y_yz_xxxyz, g_0_y_yz_xxxz, g_0_y_yz_xxxzz, g_0_y_yz_xxyy, g_0_y_yz_xxyyz, g_0_y_yz_xxyz, g_0_y_yz_xxyzz, g_0_y_yz_xxzz, g_0_y_yz_xxzzz, g_0_y_yz_xyyy, g_0_y_yz_xyyyz, g_0_y_yz_xyyz, g_0_y_yz_xyyzz, g_0_y_yz_xyzz, g_0_y_yz_xyzzz, g_0_y_yz_xzzz, g_0_y_yz_xzzzz, g_0_y_yz_yyyy, g_0_y_yz_yyyyz, g_0_y_yz_yyyz, g_0_y_yz_yyyzz, g_0_y_yz_yyzz, g_0_y_yz_yyzzz, g_0_y_yz_yzzz, g_0_y_yz_yzzzz, g_0_y_yz_zzzz, g_0_y_yz_zzzzz, g_0_y_yzz_xxxx, g_0_y_yzz_xxxy, g_0_y_yzz_xxxz, g_0_y_yzz_xxyy, g_0_y_yzz_xxyz, g_0_y_yzz_xxzz, g_0_y_yzz_xyyy, g_0_y_yzz_xyyz, g_0_y_yzz_xyzz, g_0_y_yzz_xzzz, g_0_y_yzz_yyyy, g_0_y_yzz_yyyz, g_0_y_yzz_yyzz, g_0_y_yzz_yzzz, g_0_y_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzz_xxxx[k] = -g_0_y_yz_xxxx[k] * ab_z + g_0_y_yz_xxxxz[k];

                g_0_y_yzz_xxxy[k] = -g_0_y_yz_xxxy[k] * ab_z + g_0_y_yz_xxxyz[k];

                g_0_y_yzz_xxxz[k] = -g_0_y_yz_xxxz[k] * ab_z + g_0_y_yz_xxxzz[k];

                g_0_y_yzz_xxyy[k] = -g_0_y_yz_xxyy[k] * ab_z + g_0_y_yz_xxyyz[k];

                g_0_y_yzz_xxyz[k] = -g_0_y_yz_xxyz[k] * ab_z + g_0_y_yz_xxyzz[k];

                g_0_y_yzz_xxzz[k] = -g_0_y_yz_xxzz[k] * ab_z + g_0_y_yz_xxzzz[k];

                g_0_y_yzz_xyyy[k] = -g_0_y_yz_xyyy[k] * ab_z + g_0_y_yz_xyyyz[k];

                g_0_y_yzz_xyyz[k] = -g_0_y_yz_xyyz[k] * ab_z + g_0_y_yz_xyyzz[k];

                g_0_y_yzz_xyzz[k] = -g_0_y_yz_xyzz[k] * ab_z + g_0_y_yz_xyzzz[k];

                g_0_y_yzz_xzzz[k] = -g_0_y_yz_xzzz[k] * ab_z + g_0_y_yz_xzzzz[k];

                g_0_y_yzz_yyyy[k] = -g_0_y_yz_yyyy[k] * ab_z + g_0_y_yz_yyyyz[k];

                g_0_y_yzz_yyyz[k] = -g_0_y_yz_yyyz[k] * ab_z + g_0_y_yz_yyyzz[k];

                g_0_y_yzz_yyzz[k] = -g_0_y_yz_yyzz[k] * ab_z + g_0_y_yz_yyzzz[k];

                g_0_y_yzz_yzzz[k] = -g_0_y_yz_yzzz[k] * ab_z + g_0_y_yz_yzzzz[k];

                g_0_y_yzz_zzzz[k] = -g_0_y_yz_zzzz[k] * ab_z + g_0_y_yz_zzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzz_xxxx = cbuffer.data(fg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_zzz_xxxy = cbuffer.data(fg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_zzz_xxxz = cbuffer.data(fg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_zzz_xxyy = cbuffer.data(fg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_zzz_xxyz = cbuffer.data(fg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_zzz_xxzz = cbuffer.data(fg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_zzz_xyyy = cbuffer.data(fg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_zzz_xyyz = cbuffer.data(fg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_zzz_xyzz = cbuffer.data(fg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_zzz_xzzz = cbuffer.data(fg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_zzz_yyyy = cbuffer.data(fg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_zzz_yyyz = cbuffer.data(fg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_zzz_yyzz = cbuffer.data(fg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_zzz_yzzz = cbuffer.data(fg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_zzz_zzzz = cbuffer.data(fg_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_xxxx, g_0_y_zz_xxxxz, g_0_y_zz_xxxy, g_0_y_zz_xxxyz, g_0_y_zz_xxxz, g_0_y_zz_xxxzz, g_0_y_zz_xxyy, g_0_y_zz_xxyyz, g_0_y_zz_xxyz, g_0_y_zz_xxyzz, g_0_y_zz_xxzz, g_0_y_zz_xxzzz, g_0_y_zz_xyyy, g_0_y_zz_xyyyz, g_0_y_zz_xyyz, g_0_y_zz_xyyzz, g_0_y_zz_xyzz, g_0_y_zz_xyzzz, g_0_y_zz_xzzz, g_0_y_zz_xzzzz, g_0_y_zz_yyyy, g_0_y_zz_yyyyz, g_0_y_zz_yyyz, g_0_y_zz_yyyzz, g_0_y_zz_yyzz, g_0_y_zz_yyzzz, g_0_y_zz_yzzz, g_0_y_zz_yzzzz, g_0_y_zz_zzzz, g_0_y_zz_zzzzz, g_0_y_zzz_xxxx, g_0_y_zzz_xxxy, g_0_y_zzz_xxxz, g_0_y_zzz_xxyy, g_0_y_zzz_xxyz, g_0_y_zzz_xxzz, g_0_y_zzz_xyyy, g_0_y_zzz_xyyz, g_0_y_zzz_xyzz, g_0_y_zzz_xzzz, g_0_y_zzz_yyyy, g_0_y_zzz_yyyz, g_0_y_zzz_yyzz, g_0_y_zzz_yzzz, g_0_y_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzz_xxxx[k] = -g_0_y_zz_xxxx[k] * ab_z + g_0_y_zz_xxxxz[k];

                g_0_y_zzz_xxxy[k] = -g_0_y_zz_xxxy[k] * ab_z + g_0_y_zz_xxxyz[k];

                g_0_y_zzz_xxxz[k] = -g_0_y_zz_xxxz[k] * ab_z + g_0_y_zz_xxxzz[k];

                g_0_y_zzz_xxyy[k] = -g_0_y_zz_xxyy[k] * ab_z + g_0_y_zz_xxyyz[k];

                g_0_y_zzz_xxyz[k] = -g_0_y_zz_xxyz[k] * ab_z + g_0_y_zz_xxyzz[k];

                g_0_y_zzz_xxzz[k] = -g_0_y_zz_xxzz[k] * ab_z + g_0_y_zz_xxzzz[k];

                g_0_y_zzz_xyyy[k] = -g_0_y_zz_xyyy[k] * ab_z + g_0_y_zz_xyyyz[k];

                g_0_y_zzz_xyyz[k] = -g_0_y_zz_xyyz[k] * ab_z + g_0_y_zz_xyyzz[k];

                g_0_y_zzz_xyzz[k] = -g_0_y_zz_xyzz[k] * ab_z + g_0_y_zz_xyzzz[k];

                g_0_y_zzz_xzzz[k] = -g_0_y_zz_xzzz[k] * ab_z + g_0_y_zz_xzzzz[k];

                g_0_y_zzz_yyyy[k] = -g_0_y_zz_yyyy[k] * ab_z + g_0_y_zz_yyyyz[k];

                g_0_y_zzz_yyyz[k] = -g_0_y_zz_yyyz[k] * ab_z + g_0_y_zz_yyyzz[k];

                g_0_y_zzz_yyzz[k] = -g_0_y_zz_yyzz[k] * ab_z + g_0_y_zz_yyzzz[k];

                g_0_y_zzz_yzzz[k] = -g_0_y_zz_yzzz[k] * ab_z + g_0_y_zz_yzzzz[k];

                g_0_y_zzz_zzzz[k] = -g_0_y_zz_zzzz[k] * ab_z + g_0_y_zz_zzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxx_xxxx = cbuffer.data(fg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxx_xxxy = cbuffer.data(fg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxx_xxxz = cbuffer.data(fg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxx_xxyy = cbuffer.data(fg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxx_xxyz = cbuffer.data(fg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxx_xxzz = cbuffer.data(fg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxx_xyyy = cbuffer.data(fg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxx_xyyz = cbuffer.data(fg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxx_xyzz = cbuffer.data(fg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxx_xzzz = cbuffer.data(fg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxx_yyyy = cbuffer.data(fg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxx_yyyz = cbuffer.data(fg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xxx_yyzz = cbuffer.data(fg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xxx_yzzz = cbuffer.data(fg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xxx_zzzz = cbuffer.data(fg_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_xxxx, g_0_z_xx_xxxxx, g_0_z_xx_xxxxy, g_0_z_xx_xxxxz, g_0_z_xx_xxxy, g_0_z_xx_xxxyy, g_0_z_xx_xxxyz, g_0_z_xx_xxxz, g_0_z_xx_xxxzz, g_0_z_xx_xxyy, g_0_z_xx_xxyyy, g_0_z_xx_xxyyz, g_0_z_xx_xxyz, g_0_z_xx_xxyzz, g_0_z_xx_xxzz, g_0_z_xx_xxzzz, g_0_z_xx_xyyy, g_0_z_xx_xyyyy, g_0_z_xx_xyyyz, g_0_z_xx_xyyz, g_0_z_xx_xyyzz, g_0_z_xx_xyzz, g_0_z_xx_xyzzz, g_0_z_xx_xzzz, g_0_z_xx_xzzzz, g_0_z_xx_yyyy, g_0_z_xx_yyyz, g_0_z_xx_yyzz, g_0_z_xx_yzzz, g_0_z_xx_zzzz, g_0_z_xxx_xxxx, g_0_z_xxx_xxxy, g_0_z_xxx_xxxz, g_0_z_xxx_xxyy, g_0_z_xxx_xxyz, g_0_z_xxx_xxzz, g_0_z_xxx_xyyy, g_0_z_xxx_xyyz, g_0_z_xxx_xyzz, g_0_z_xxx_xzzz, g_0_z_xxx_yyyy, g_0_z_xxx_yyyz, g_0_z_xxx_yyzz, g_0_z_xxx_yzzz, g_0_z_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxx_xxxx[k] = -g_0_z_xx_xxxx[k] * ab_x + g_0_z_xx_xxxxx[k];

                g_0_z_xxx_xxxy[k] = -g_0_z_xx_xxxy[k] * ab_x + g_0_z_xx_xxxxy[k];

                g_0_z_xxx_xxxz[k] = -g_0_z_xx_xxxz[k] * ab_x + g_0_z_xx_xxxxz[k];

                g_0_z_xxx_xxyy[k] = -g_0_z_xx_xxyy[k] * ab_x + g_0_z_xx_xxxyy[k];

                g_0_z_xxx_xxyz[k] = -g_0_z_xx_xxyz[k] * ab_x + g_0_z_xx_xxxyz[k];

                g_0_z_xxx_xxzz[k] = -g_0_z_xx_xxzz[k] * ab_x + g_0_z_xx_xxxzz[k];

                g_0_z_xxx_xyyy[k] = -g_0_z_xx_xyyy[k] * ab_x + g_0_z_xx_xxyyy[k];

                g_0_z_xxx_xyyz[k] = -g_0_z_xx_xyyz[k] * ab_x + g_0_z_xx_xxyyz[k];

                g_0_z_xxx_xyzz[k] = -g_0_z_xx_xyzz[k] * ab_x + g_0_z_xx_xxyzz[k];

                g_0_z_xxx_xzzz[k] = -g_0_z_xx_xzzz[k] * ab_x + g_0_z_xx_xxzzz[k];

                g_0_z_xxx_yyyy[k] = -g_0_z_xx_yyyy[k] * ab_x + g_0_z_xx_xyyyy[k];

                g_0_z_xxx_yyyz[k] = -g_0_z_xx_yyyz[k] * ab_x + g_0_z_xx_xyyyz[k];

                g_0_z_xxx_yyzz[k] = -g_0_z_xx_yyzz[k] * ab_x + g_0_z_xx_xyyzz[k];

                g_0_z_xxx_yzzz[k] = -g_0_z_xx_yzzz[k] * ab_x + g_0_z_xx_xyzzz[k];

                g_0_z_xxx_zzzz[k] = -g_0_z_xx_zzzz[k] * ab_x + g_0_z_xx_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxy_xxxx = cbuffer.data(fg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xxy_xxxy = cbuffer.data(fg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xxy_xxxz = cbuffer.data(fg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xxy_xxyy = cbuffer.data(fg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xxy_xxyz = cbuffer.data(fg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xxy_xxzz = cbuffer.data(fg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xxy_xyyy = cbuffer.data(fg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xxy_xyyz = cbuffer.data(fg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xxy_xyzz = cbuffer.data(fg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xxy_xzzz = cbuffer.data(fg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xxy_yyyy = cbuffer.data(fg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xxy_yyyz = cbuffer.data(fg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xxy_yyzz = cbuffer.data(fg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xxy_yzzz = cbuffer.data(fg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xxy_zzzz = cbuffer.data(fg_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxy_xxxx, g_0_z_xxy_xxxy, g_0_z_xxy_xxxz, g_0_z_xxy_xxyy, g_0_z_xxy_xxyz, g_0_z_xxy_xxzz, g_0_z_xxy_xyyy, g_0_z_xxy_xyyz, g_0_z_xxy_xyzz, g_0_z_xxy_xzzz, g_0_z_xxy_yyyy, g_0_z_xxy_yyyz, g_0_z_xxy_yyzz, g_0_z_xxy_yzzz, g_0_z_xxy_zzzz, g_0_z_xy_xxxx, g_0_z_xy_xxxxx, g_0_z_xy_xxxxy, g_0_z_xy_xxxxz, g_0_z_xy_xxxy, g_0_z_xy_xxxyy, g_0_z_xy_xxxyz, g_0_z_xy_xxxz, g_0_z_xy_xxxzz, g_0_z_xy_xxyy, g_0_z_xy_xxyyy, g_0_z_xy_xxyyz, g_0_z_xy_xxyz, g_0_z_xy_xxyzz, g_0_z_xy_xxzz, g_0_z_xy_xxzzz, g_0_z_xy_xyyy, g_0_z_xy_xyyyy, g_0_z_xy_xyyyz, g_0_z_xy_xyyz, g_0_z_xy_xyyzz, g_0_z_xy_xyzz, g_0_z_xy_xyzzz, g_0_z_xy_xzzz, g_0_z_xy_xzzzz, g_0_z_xy_yyyy, g_0_z_xy_yyyz, g_0_z_xy_yyzz, g_0_z_xy_yzzz, g_0_z_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxy_xxxx[k] = -g_0_z_xy_xxxx[k] * ab_x + g_0_z_xy_xxxxx[k];

                g_0_z_xxy_xxxy[k] = -g_0_z_xy_xxxy[k] * ab_x + g_0_z_xy_xxxxy[k];

                g_0_z_xxy_xxxz[k] = -g_0_z_xy_xxxz[k] * ab_x + g_0_z_xy_xxxxz[k];

                g_0_z_xxy_xxyy[k] = -g_0_z_xy_xxyy[k] * ab_x + g_0_z_xy_xxxyy[k];

                g_0_z_xxy_xxyz[k] = -g_0_z_xy_xxyz[k] * ab_x + g_0_z_xy_xxxyz[k];

                g_0_z_xxy_xxzz[k] = -g_0_z_xy_xxzz[k] * ab_x + g_0_z_xy_xxxzz[k];

                g_0_z_xxy_xyyy[k] = -g_0_z_xy_xyyy[k] * ab_x + g_0_z_xy_xxyyy[k];

                g_0_z_xxy_xyyz[k] = -g_0_z_xy_xyyz[k] * ab_x + g_0_z_xy_xxyyz[k];

                g_0_z_xxy_xyzz[k] = -g_0_z_xy_xyzz[k] * ab_x + g_0_z_xy_xxyzz[k];

                g_0_z_xxy_xzzz[k] = -g_0_z_xy_xzzz[k] * ab_x + g_0_z_xy_xxzzz[k];

                g_0_z_xxy_yyyy[k] = -g_0_z_xy_yyyy[k] * ab_x + g_0_z_xy_xyyyy[k];

                g_0_z_xxy_yyyz[k] = -g_0_z_xy_yyyz[k] * ab_x + g_0_z_xy_xyyyz[k];

                g_0_z_xxy_yyzz[k] = -g_0_z_xy_yyzz[k] * ab_x + g_0_z_xy_xyyzz[k];

                g_0_z_xxy_yzzz[k] = -g_0_z_xy_yzzz[k] * ab_x + g_0_z_xy_xyzzz[k];

                g_0_z_xxy_zzzz[k] = -g_0_z_xy_zzzz[k] * ab_x + g_0_z_xy_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxz_xxxx = cbuffer.data(fg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xxz_xxxy = cbuffer.data(fg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xxz_xxxz = cbuffer.data(fg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xxz_xxyy = cbuffer.data(fg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xxz_xxyz = cbuffer.data(fg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xxz_xxzz = cbuffer.data(fg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxz_xyyy = cbuffer.data(fg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxz_xyyz = cbuffer.data(fg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxz_xyzz = cbuffer.data(fg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxz_xzzz = cbuffer.data(fg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxz_yyyy = cbuffer.data(fg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxz_yyyz = cbuffer.data(fg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxz_yyzz = cbuffer.data(fg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxz_yzzz = cbuffer.data(fg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxz_zzzz = cbuffer.data(fg_geom_01_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxz_xxxx, g_0_z_xxz_xxxy, g_0_z_xxz_xxxz, g_0_z_xxz_xxyy, g_0_z_xxz_xxyz, g_0_z_xxz_xxzz, g_0_z_xxz_xyyy, g_0_z_xxz_xyyz, g_0_z_xxz_xyzz, g_0_z_xxz_xzzz, g_0_z_xxz_yyyy, g_0_z_xxz_yyyz, g_0_z_xxz_yyzz, g_0_z_xxz_yzzz, g_0_z_xxz_zzzz, g_0_z_xz_xxxx, g_0_z_xz_xxxxx, g_0_z_xz_xxxxy, g_0_z_xz_xxxxz, g_0_z_xz_xxxy, g_0_z_xz_xxxyy, g_0_z_xz_xxxyz, g_0_z_xz_xxxz, g_0_z_xz_xxxzz, g_0_z_xz_xxyy, g_0_z_xz_xxyyy, g_0_z_xz_xxyyz, g_0_z_xz_xxyz, g_0_z_xz_xxyzz, g_0_z_xz_xxzz, g_0_z_xz_xxzzz, g_0_z_xz_xyyy, g_0_z_xz_xyyyy, g_0_z_xz_xyyyz, g_0_z_xz_xyyz, g_0_z_xz_xyyzz, g_0_z_xz_xyzz, g_0_z_xz_xyzzz, g_0_z_xz_xzzz, g_0_z_xz_xzzzz, g_0_z_xz_yyyy, g_0_z_xz_yyyz, g_0_z_xz_yyzz, g_0_z_xz_yzzz, g_0_z_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxz_xxxx[k] = -g_0_z_xz_xxxx[k] * ab_x + g_0_z_xz_xxxxx[k];

                g_0_z_xxz_xxxy[k] = -g_0_z_xz_xxxy[k] * ab_x + g_0_z_xz_xxxxy[k];

                g_0_z_xxz_xxxz[k] = -g_0_z_xz_xxxz[k] * ab_x + g_0_z_xz_xxxxz[k];

                g_0_z_xxz_xxyy[k] = -g_0_z_xz_xxyy[k] * ab_x + g_0_z_xz_xxxyy[k];

                g_0_z_xxz_xxyz[k] = -g_0_z_xz_xxyz[k] * ab_x + g_0_z_xz_xxxyz[k];

                g_0_z_xxz_xxzz[k] = -g_0_z_xz_xxzz[k] * ab_x + g_0_z_xz_xxxzz[k];

                g_0_z_xxz_xyyy[k] = -g_0_z_xz_xyyy[k] * ab_x + g_0_z_xz_xxyyy[k];

                g_0_z_xxz_xyyz[k] = -g_0_z_xz_xyyz[k] * ab_x + g_0_z_xz_xxyyz[k];

                g_0_z_xxz_xyzz[k] = -g_0_z_xz_xyzz[k] * ab_x + g_0_z_xz_xxyzz[k];

                g_0_z_xxz_xzzz[k] = -g_0_z_xz_xzzz[k] * ab_x + g_0_z_xz_xxzzz[k];

                g_0_z_xxz_yyyy[k] = -g_0_z_xz_yyyy[k] * ab_x + g_0_z_xz_xyyyy[k];

                g_0_z_xxz_yyyz[k] = -g_0_z_xz_yyyz[k] * ab_x + g_0_z_xz_xyyyz[k];

                g_0_z_xxz_yyzz[k] = -g_0_z_xz_yyzz[k] * ab_x + g_0_z_xz_xyyzz[k];

                g_0_z_xxz_yzzz[k] = -g_0_z_xz_yzzz[k] * ab_x + g_0_z_xz_xyzzz[k];

                g_0_z_xxz_zzzz[k] = -g_0_z_xz_zzzz[k] * ab_x + g_0_z_xz_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyy_xxxx = cbuffer.data(fg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xyy_xxxy = cbuffer.data(fg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xyy_xxxz = cbuffer.data(fg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xyy_xxyy = cbuffer.data(fg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xyy_xxyz = cbuffer.data(fg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xyy_xxzz = cbuffer.data(fg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xyy_xyyy = cbuffer.data(fg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xyy_xyyz = cbuffer.data(fg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xyy_xyzz = cbuffer.data(fg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xyy_xzzz = cbuffer.data(fg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xyy_yyyy = cbuffer.data(fg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xyy_yyyz = cbuffer.data(fg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xyy_yyzz = cbuffer.data(fg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xyy_yzzz = cbuffer.data(fg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xyy_zzzz = cbuffer.data(fg_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyy_xxxx, g_0_z_xyy_xxxy, g_0_z_xyy_xxxz, g_0_z_xyy_xxyy, g_0_z_xyy_xxyz, g_0_z_xyy_xxzz, g_0_z_xyy_xyyy, g_0_z_xyy_xyyz, g_0_z_xyy_xyzz, g_0_z_xyy_xzzz, g_0_z_xyy_yyyy, g_0_z_xyy_yyyz, g_0_z_xyy_yyzz, g_0_z_xyy_yzzz, g_0_z_xyy_zzzz, g_0_z_yy_xxxx, g_0_z_yy_xxxxx, g_0_z_yy_xxxxy, g_0_z_yy_xxxxz, g_0_z_yy_xxxy, g_0_z_yy_xxxyy, g_0_z_yy_xxxyz, g_0_z_yy_xxxz, g_0_z_yy_xxxzz, g_0_z_yy_xxyy, g_0_z_yy_xxyyy, g_0_z_yy_xxyyz, g_0_z_yy_xxyz, g_0_z_yy_xxyzz, g_0_z_yy_xxzz, g_0_z_yy_xxzzz, g_0_z_yy_xyyy, g_0_z_yy_xyyyy, g_0_z_yy_xyyyz, g_0_z_yy_xyyz, g_0_z_yy_xyyzz, g_0_z_yy_xyzz, g_0_z_yy_xyzzz, g_0_z_yy_xzzz, g_0_z_yy_xzzzz, g_0_z_yy_yyyy, g_0_z_yy_yyyz, g_0_z_yy_yyzz, g_0_z_yy_yzzz, g_0_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyy_xxxx[k] = -g_0_z_yy_xxxx[k] * ab_x + g_0_z_yy_xxxxx[k];

                g_0_z_xyy_xxxy[k] = -g_0_z_yy_xxxy[k] * ab_x + g_0_z_yy_xxxxy[k];

                g_0_z_xyy_xxxz[k] = -g_0_z_yy_xxxz[k] * ab_x + g_0_z_yy_xxxxz[k];

                g_0_z_xyy_xxyy[k] = -g_0_z_yy_xxyy[k] * ab_x + g_0_z_yy_xxxyy[k];

                g_0_z_xyy_xxyz[k] = -g_0_z_yy_xxyz[k] * ab_x + g_0_z_yy_xxxyz[k];

                g_0_z_xyy_xxzz[k] = -g_0_z_yy_xxzz[k] * ab_x + g_0_z_yy_xxxzz[k];

                g_0_z_xyy_xyyy[k] = -g_0_z_yy_xyyy[k] * ab_x + g_0_z_yy_xxyyy[k];

                g_0_z_xyy_xyyz[k] = -g_0_z_yy_xyyz[k] * ab_x + g_0_z_yy_xxyyz[k];

                g_0_z_xyy_xyzz[k] = -g_0_z_yy_xyzz[k] * ab_x + g_0_z_yy_xxyzz[k];

                g_0_z_xyy_xzzz[k] = -g_0_z_yy_xzzz[k] * ab_x + g_0_z_yy_xxzzz[k];

                g_0_z_xyy_yyyy[k] = -g_0_z_yy_yyyy[k] * ab_x + g_0_z_yy_xyyyy[k];

                g_0_z_xyy_yyyz[k] = -g_0_z_yy_yyyz[k] * ab_x + g_0_z_yy_xyyyz[k];

                g_0_z_xyy_yyzz[k] = -g_0_z_yy_yyzz[k] * ab_x + g_0_z_yy_xyyzz[k];

                g_0_z_xyy_yzzz[k] = -g_0_z_yy_yzzz[k] * ab_x + g_0_z_yy_xyzzz[k];

                g_0_z_xyy_zzzz[k] = -g_0_z_yy_zzzz[k] * ab_x + g_0_z_yy_xzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyz_xxxx = cbuffer.data(fg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xyz_xxxy = cbuffer.data(fg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xyz_xxxz = cbuffer.data(fg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xyz_xxyy = cbuffer.data(fg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xyz_xxyz = cbuffer.data(fg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xyz_xxzz = cbuffer.data(fg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xyz_xyyy = cbuffer.data(fg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xyz_xyyz = cbuffer.data(fg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xyz_xyzz = cbuffer.data(fg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xyz_xzzz = cbuffer.data(fg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xyz_yyyy = cbuffer.data(fg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xyz_yyyz = cbuffer.data(fg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xyz_yyzz = cbuffer.data(fg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xyz_yzzz = cbuffer.data(fg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xyz_zzzz = cbuffer.data(fg_geom_01_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyz_xxxx, g_0_z_xyz_xxxy, g_0_z_xyz_xxxz, g_0_z_xyz_xxyy, g_0_z_xyz_xxyz, g_0_z_xyz_xxzz, g_0_z_xyz_xyyy, g_0_z_xyz_xyyz, g_0_z_xyz_xyzz, g_0_z_xyz_xzzz, g_0_z_xyz_yyyy, g_0_z_xyz_yyyz, g_0_z_xyz_yyzz, g_0_z_xyz_yzzz, g_0_z_xyz_zzzz, g_0_z_yz_xxxx, g_0_z_yz_xxxxx, g_0_z_yz_xxxxy, g_0_z_yz_xxxxz, g_0_z_yz_xxxy, g_0_z_yz_xxxyy, g_0_z_yz_xxxyz, g_0_z_yz_xxxz, g_0_z_yz_xxxzz, g_0_z_yz_xxyy, g_0_z_yz_xxyyy, g_0_z_yz_xxyyz, g_0_z_yz_xxyz, g_0_z_yz_xxyzz, g_0_z_yz_xxzz, g_0_z_yz_xxzzz, g_0_z_yz_xyyy, g_0_z_yz_xyyyy, g_0_z_yz_xyyyz, g_0_z_yz_xyyz, g_0_z_yz_xyyzz, g_0_z_yz_xyzz, g_0_z_yz_xyzzz, g_0_z_yz_xzzz, g_0_z_yz_xzzzz, g_0_z_yz_yyyy, g_0_z_yz_yyyz, g_0_z_yz_yyzz, g_0_z_yz_yzzz, g_0_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyz_xxxx[k] = -g_0_z_yz_xxxx[k] * ab_x + g_0_z_yz_xxxxx[k];

                g_0_z_xyz_xxxy[k] = -g_0_z_yz_xxxy[k] * ab_x + g_0_z_yz_xxxxy[k];

                g_0_z_xyz_xxxz[k] = -g_0_z_yz_xxxz[k] * ab_x + g_0_z_yz_xxxxz[k];

                g_0_z_xyz_xxyy[k] = -g_0_z_yz_xxyy[k] * ab_x + g_0_z_yz_xxxyy[k];

                g_0_z_xyz_xxyz[k] = -g_0_z_yz_xxyz[k] * ab_x + g_0_z_yz_xxxyz[k];

                g_0_z_xyz_xxzz[k] = -g_0_z_yz_xxzz[k] * ab_x + g_0_z_yz_xxxzz[k];

                g_0_z_xyz_xyyy[k] = -g_0_z_yz_xyyy[k] * ab_x + g_0_z_yz_xxyyy[k];

                g_0_z_xyz_xyyz[k] = -g_0_z_yz_xyyz[k] * ab_x + g_0_z_yz_xxyyz[k];

                g_0_z_xyz_xyzz[k] = -g_0_z_yz_xyzz[k] * ab_x + g_0_z_yz_xxyzz[k];

                g_0_z_xyz_xzzz[k] = -g_0_z_yz_xzzz[k] * ab_x + g_0_z_yz_xxzzz[k];

                g_0_z_xyz_yyyy[k] = -g_0_z_yz_yyyy[k] * ab_x + g_0_z_yz_xyyyy[k];

                g_0_z_xyz_yyyz[k] = -g_0_z_yz_yyyz[k] * ab_x + g_0_z_yz_xyyyz[k];

                g_0_z_xyz_yyzz[k] = -g_0_z_yz_yyzz[k] * ab_x + g_0_z_yz_xyyzz[k];

                g_0_z_xyz_yzzz[k] = -g_0_z_yz_yzzz[k] * ab_x + g_0_z_yz_xyzzz[k];

                g_0_z_xyz_zzzz[k] = -g_0_z_yz_zzzz[k] * ab_x + g_0_z_yz_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzz_xxxx = cbuffer.data(fg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xzz_xxxy = cbuffer.data(fg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xzz_xxxz = cbuffer.data(fg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xzz_xxyy = cbuffer.data(fg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xzz_xxyz = cbuffer.data(fg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xzz_xxzz = cbuffer.data(fg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xzz_xyyy = cbuffer.data(fg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xzz_xyyz = cbuffer.data(fg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xzz_xyzz = cbuffer.data(fg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xzz_xzzz = cbuffer.data(fg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xzz_yyyy = cbuffer.data(fg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xzz_yyyz = cbuffer.data(fg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xzz_yyzz = cbuffer.data(fg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xzz_yzzz = cbuffer.data(fg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xzz_zzzz = cbuffer.data(fg_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzz_xxxx, g_0_z_xzz_xxxy, g_0_z_xzz_xxxz, g_0_z_xzz_xxyy, g_0_z_xzz_xxyz, g_0_z_xzz_xxzz, g_0_z_xzz_xyyy, g_0_z_xzz_xyyz, g_0_z_xzz_xyzz, g_0_z_xzz_xzzz, g_0_z_xzz_yyyy, g_0_z_xzz_yyyz, g_0_z_xzz_yyzz, g_0_z_xzz_yzzz, g_0_z_xzz_zzzz, g_0_z_zz_xxxx, g_0_z_zz_xxxxx, g_0_z_zz_xxxxy, g_0_z_zz_xxxxz, g_0_z_zz_xxxy, g_0_z_zz_xxxyy, g_0_z_zz_xxxyz, g_0_z_zz_xxxz, g_0_z_zz_xxxzz, g_0_z_zz_xxyy, g_0_z_zz_xxyyy, g_0_z_zz_xxyyz, g_0_z_zz_xxyz, g_0_z_zz_xxyzz, g_0_z_zz_xxzz, g_0_z_zz_xxzzz, g_0_z_zz_xyyy, g_0_z_zz_xyyyy, g_0_z_zz_xyyyz, g_0_z_zz_xyyz, g_0_z_zz_xyyzz, g_0_z_zz_xyzz, g_0_z_zz_xyzzz, g_0_z_zz_xzzz, g_0_z_zz_xzzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyz, g_0_z_zz_yyzz, g_0_z_zz_yzzz, g_0_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzz_xxxx[k] = -g_0_z_zz_xxxx[k] * ab_x + g_0_z_zz_xxxxx[k];

                g_0_z_xzz_xxxy[k] = -g_0_z_zz_xxxy[k] * ab_x + g_0_z_zz_xxxxy[k];

                g_0_z_xzz_xxxz[k] = -g_0_z_zz_xxxz[k] * ab_x + g_0_z_zz_xxxxz[k];

                g_0_z_xzz_xxyy[k] = -g_0_z_zz_xxyy[k] * ab_x + g_0_z_zz_xxxyy[k];

                g_0_z_xzz_xxyz[k] = -g_0_z_zz_xxyz[k] * ab_x + g_0_z_zz_xxxyz[k];

                g_0_z_xzz_xxzz[k] = -g_0_z_zz_xxzz[k] * ab_x + g_0_z_zz_xxxzz[k];

                g_0_z_xzz_xyyy[k] = -g_0_z_zz_xyyy[k] * ab_x + g_0_z_zz_xxyyy[k];

                g_0_z_xzz_xyyz[k] = -g_0_z_zz_xyyz[k] * ab_x + g_0_z_zz_xxyyz[k];

                g_0_z_xzz_xyzz[k] = -g_0_z_zz_xyzz[k] * ab_x + g_0_z_zz_xxyzz[k];

                g_0_z_xzz_xzzz[k] = -g_0_z_zz_xzzz[k] * ab_x + g_0_z_zz_xxzzz[k];

                g_0_z_xzz_yyyy[k] = -g_0_z_zz_yyyy[k] * ab_x + g_0_z_zz_xyyyy[k];

                g_0_z_xzz_yyyz[k] = -g_0_z_zz_yyyz[k] * ab_x + g_0_z_zz_xyyyz[k];

                g_0_z_xzz_yyzz[k] = -g_0_z_zz_yyzz[k] * ab_x + g_0_z_zz_xyyzz[k];

                g_0_z_xzz_yzzz[k] = -g_0_z_zz_yzzz[k] * ab_x + g_0_z_zz_xyzzz[k];

                g_0_z_xzz_zzzz[k] = -g_0_z_zz_zzzz[k] * ab_x + g_0_z_zz_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyy_xxxx = cbuffer.data(fg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_yyy_xxxy = cbuffer.data(fg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_yyy_xxxz = cbuffer.data(fg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_yyy_xxyy = cbuffer.data(fg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_yyy_xxyz = cbuffer.data(fg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_yyy_xxzz = cbuffer.data(fg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_yyy_xyyy = cbuffer.data(fg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_yyy_xyyz = cbuffer.data(fg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_yyy_xyzz = cbuffer.data(fg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_yyy_xzzz = cbuffer.data(fg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_yyy_yyyy = cbuffer.data(fg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_yyy_yyyz = cbuffer.data(fg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_yyy_yyzz = cbuffer.data(fg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_yyy_yzzz = cbuffer.data(fg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_yyy_zzzz = cbuffer.data(fg_geom_01_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_xxxx, g_0_z_yy_xxxxy, g_0_z_yy_xxxy, g_0_z_yy_xxxyy, g_0_z_yy_xxxyz, g_0_z_yy_xxxz, g_0_z_yy_xxyy, g_0_z_yy_xxyyy, g_0_z_yy_xxyyz, g_0_z_yy_xxyz, g_0_z_yy_xxyzz, g_0_z_yy_xxzz, g_0_z_yy_xyyy, g_0_z_yy_xyyyy, g_0_z_yy_xyyyz, g_0_z_yy_xyyz, g_0_z_yy_xyyzz, g_0_z_yy_xyzz, g_0_z_yy_xyzzz, g_0_z_yy_xzzz, g_0_z_yy_yyyy, g_0_z_yy_yyyyy, g_0_z_yy_yyyyz, g_0_z_yy_yyyz, g_0_z_yy_yyyzz, g_0_z_yy_yyzz, g_0_z_yy_yyzzz, g_0_z_yy_yzzz, g_0_z_yy_yzzzz, g_0_z_yy_zzzz, g_0_z_yyy_xxxx, g_0_z_yyy_xxxy, g_0_z_yyy_xxxz, g_0_z_yyy_xxyy, g_0_z_yyy_xxyz, g_0_z_yyy_xxzz, g_0_z_yyy_xyyy, g_0_z_yyy_xyyz, g_0_z_yyy_xyzz, g_0_z_yyy_xzzz, g_0_z_yyy_yyyy, g_0_z_yyy_yyyz, g_0_z_yyy_yyzz, g_0_z_yyy_yzzz, g_0_z_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyy_xxxx[k] = -g_0_z_yy_xxxx[k] * ab_y + g_0_z_yy_xxxxy[k];

                g_0_z_yyy_xxxy[k] = -g_0_z_yy_xxxy[k] * ab_y + g_0_z_yy_xxxyy[k];

                g_0_z_yyy_xxxz[k] = -g_0_z_yy_xxxz[k] * ab_y + g_0_z_yy_xxxyz[k];

                g_0_z_yyy_xxyy[k] = -g_0_z_yy_xxyy[k] * ab_y + g_0_z_yy_xxyyy[k];

                g_0_z_yyy_xxyz[k] = -g_0_z_yy_xxyz[k] * ab_y + g_0_z_yy_xxyyz[k];

                g_0_z_yyy_xxzz[k] = -g_0_z_yy_xxzz[k] * ab_y + g_0_z_yy_xxyzz[k];

                g_0_z_yyy_xyyy[k] = -g_0_z_yy_xyyy[k] * ab_y + g_0_z_yy_xyyyy[k];

                g_0_z_yyy_xyyz[k] = -g_0_z_yy_xyyz[k] * ab_y + g_0_z_yy_xyyyz[k];

                g_0_z_yyy_xyzz[k] = -g_0_z_yy_xyzz[k] * ab_y + g_0_z_yy_xyyzz[k];

                g_0_z_yyy_xzzz[k] = -g_0_z_yy_xzzz[k] * ab_y + g_0_z_yy_xyzzz[k];

                g_0_z_yyy_yyyy[k] = -g_0_z_yy_yyyy[k] * ab_y + g_0_z_yy_yyyyy[k];

                g_0_z_yyy_yyyz[k] = -g_0_z_yy_yyyz[k] * ab_y + g_0_z_yy_yyyyz[k];

                g_0_z_yyy_yyzz[k] = -g_0_z_yy_yyzz[k] * ab_y + g_0_z_yy_yyyzz[k];

                g_0_z_yyy_yzzz[k] = -g_0_z_yy_yzzz[k] * ab_y + g_0_z_yy_yyzzz[k];

                g_0_z_yyy_zzzz[k] = -g_0_z_yy_zzzz[k] * ab_y + g_0_z_yy_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyz_xxxx = cbuffer.data(fg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_yyz_xxxy = cbuffer.data(fg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_yyz_xxxz = cbuffer.data(fg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_yyz_xxyy = cbuffer.data(fg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_yyz_xxyz = cbuffer.data(fg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_yyz_xxzz = cbuffer.data(fg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_yyz_xyyy = cbuffer.data(fg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_yyz_xyyz = cbuffer.data(fg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_yyz_xyzz = cbuffer.data(fg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_yyz_xzzz = cbuffer.data(fg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_yyz_yyyy = cbuffer.data(fg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_yyz_yyyz = cbuffer.data(fg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_yyz_yyzz = cbuffer.data(fg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_yyz_yzzz = cbuffer.data(fg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_yyz_zzzz = cbuffer.data(fg_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyz_xxxx, g_0_z_yyz_xxxy, g_0_z_yyz_xxxz, g_0_z_yyz_xxyy, g_0_z_yyz_xxyz, g_0_z_yyz_xxzz, g_0_z_yyz_xyyy, g_0_z_yyz_xyyz, g_0_z_yyz_xyzz, g_0_z_yyz_xzzz, g_0_z_yyz_yyyy, g_0_z_yyz_yyyz, g_0_z_yyz_yyzz, g_0_z_yyz_yzzz, g_0_z_yyz_zzzz, g_0_z_yz_xxxx, g_0_z_yz_xxxxy, g_0_z_yz_xxxy, g_0_z_yz_xxxyy, g_0_z_yz_xxxyz, g_0_z_yz_xxxz, g_0_z_yz_xxyy, g_0_z_yz_xxyyy, g_0_z_yz_xxyyz, g_0_z_yz_xxyz, g_0_z_yz_xxyzz, g_0_z_yz_xxzz, g_0_z_yz_xyyy, g_0_z_yz_xyyyy, g_0_z_yz_xyyyz, g_0_z_yz_xyyz, g_0_z_yz_xyyzz, g_0_z_yz_xyzz, g_0_z_yz_xyzzz, g_0_z_yz_xzzz, g_0_z_yz_yyyy, g_0_z_yz_yyyyy, g_0_z_yz_yyyyz, g_0_z_yz_yyyz, g_0_z_yz_yyyzz, g_0_z_yz_yyzz, g_0_z_yz_yyzzz, g_0_z_yz_yzzz, g_0_z_yz_yzzzz, g_0_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyz_xxxx[k] = -g_0_z_yz_xxxx[k] * ab_y + g_0_z_yz_xxxxy[k];

                g_0_z_yyz_xxxy[k] = -g_0_z_yz_xxxy[k] * ab_y + g_0_z_yz_xxxyy[k];

                g_0_z_yyz_xxxz[k] = -g_0_z_yz_xxxz[k] * ab_y + g_0_z_yz_xxxyz[k];

                g_0_z_yyz_xxyy[k] = -g_0_z_yz_xxyy[k] * ab_y + g_0_z_yz_xxyyy[k];

                g_0_z_yyz_xxyz[k] = -g_0_z_yz_xxyz[k] * ab_y + g_0_z_yz_xxyyz[k];

                g_0_z_yyz_xxzz[k] = -g_0_z_yz_xxzz[k] * ab_y + g_0_z_yz_xxyzz[k];

                g_0_z_yyz_xyyy[k] = -g_0_z_yz_xyyy[k] * ab_y + g_0_z_yz_xyyyy[k];

                g_0_z_yyz_xyyz[k] = -g_0_z_yz_xyyz[k] * ab_y + g_0_z_yz_xyyyz[k];

                g_0_z_yyz_xyzz[k] = -g_0_z_yz_xyzz[k] * ab_y + g_0_z_yz_xyyzz[k];

                g_0_z_yyz_xzzz[k] = -g_0_z_yz_xzzz[k] * ab_y + g_0_z_yz_xyzzz[k];

                g_0_z_yyz_yyyy[k] = -g_0_z_yz_yyyy[k] * ab_y + g_0_z_yz_yyyyy[k];

                g_0_z_yyz_yyyz[k] = -g_0_z_yz_yyyz[k] * ab_y + g_0_z_yz_yyyyz[k];

                g_0_z_yyz_yyzz[k] = -g_0_z_yz_yyzz[k] * ab_y + g_0_z_yz_yyyzz[k];

                g_0_z_yyz_yzzz[k] = -g_0_z_yz_yzzz[k] * ab_y + g_0_z_yz_yyzzz[k];

                g_0_z_yyz_zzzz[k] = -g_0_z_yz_zzzz[k] * ab_y + g_0_z_yz_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzz_xxxx = cbuffer.data(fg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_yzz_xxxy = cbuffer.data(fg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_yzz_xxxz = cbuffer.data(fg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_yzz_xxyy = cbuffer.data(fg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_yzz_xxyz = cbuffer.data(fg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_yzz_xxzz = cbuffer.data(fg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_yzz_xyyy = cbuffer.data(fg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_yzz_xyyz = cbuffer.data(fg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_yzz_xyzz = cbuffer.data(fg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_yzz_xzzz = cbuffer.data(fg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_yzz_yyyy = cbuffer.data(fg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_yzz_yyyz = cbuffer.data(fg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_yzz_yyzz = cbuffer.data(fg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_yzz_yzzz = cbuffer.data(fg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_yzz_zzzz = cbuffer.data(fg_geom_01_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzz_xxxx, g_0_z_yzz_xxxy, g_0_z_yzz_xxxz, g_0_z_yzz_xxyy, g_0_z_yzz_xxyz, g_0_z_yzz_xxzz, g_0_z_yzz_xyyy, g_0_z_yzz_xyyz, g_0_z_yzz_xyzz, g_0_z_yzz_xzzz, g_0_z_yzz_yyyy, g_0_z_yzz_yyyz, g_0_z_yzz_yyzz, g_0_z_yzz_yzzz, g_0_z_yzz_zzzz, g_0_z_zz_xxxx, g_0_z_zz_xxxxy, g_0_z_zz_xxxy, g_0_z_zz_xxxyy, g_0_z_zz_xxxyz, g_0_z_zz_xxxz, g_0_z_zz_xxyy, g_0_z_zz_xxyyy, g_0_z_zz_xxyyz, g_0_z_zz_xxyz, g_0_z_zz_xxyzz, g_0_z_zz_xxzz, g_0_z_zz_xyyy, g_0_z_zz_xyyyy, g_0_z_zz_xyyyz, g_0_z_zz_xyyz, g_0_z_zz_xyyzz, g_0_z_zz_xyzz, g_0_z_zz_xyzzz, g_0_z_zz_xzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyyy, g_0_z_zz_yyyyz, g_0_z_zz_yyyz, g_0_z_zz_yyyzz, g_0_z_zz_yyzz, g_0_z_zz_yyzzz, g_0_z_zz_yzzz, g_0_z_zz_yzzzz, g_0_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzz_xxxx[k] = -g_0_z_zz_xxxx[k] * ab_y + g_0_z_zz_xxxxy[k];

                g_0_z_yzz_xxxy[k] = -g_0_z_zz_xxxy[k] * ab_y + g_0_z_zz_xxxyy[k];

                g_0_z_yzz_xxxz[k] = -g_0_z_zz_xxxz[k] * ab_y + g_0_z_zz_xxxyz[k];

                g_0_z_yzz_xxyy[k] = -g_0_z_zz_xxyy[k] * ab_y + g_0_z_zz_xxyyy[k];

                g_0_z_yzz_xxyz[k] = -g_0_z_zz_xxyz[k] * ab_y + g_0_z_zz_xxyyz[k];

                g_0_z_yzz_xxzz[k] = -g_0_z_zz_xxzz[k] * ab_y + g_0_z_zz_xxyzz[k];

                g_0_z_yzz_xyyy[k] = -g_0_z_zz_xyyy[k] * ab_y + g_0_z_zz_xyyyy[k];

                g_0_z_yzz_xyyz[k] = -g_0_z_zz_xyyz[k] * ab_y + g_0_z_zz_xyyyz[k];

                g_0_z_yzz_xyzz[k] = -g_0_z_zz_xyzz[k] * ab_y + g_0_z_zz_xyyzz[k];

                g_0_z_yzz_xzzz[k] = -g_0_z_zz_xzzz[k] * ab_y + g_0_z_zz_xyzzz[k];

                g_0_z_yzz_yyyy[k] = -g_0_z_zz_yyyy[k] * ab_y + g_0_z_zz_yyyyy[k];

                g_0_z_yzz_yyyz[k] = -g_0_z_zz_yyyz[k] * ab_y + g_0_z_zz_yyyyz[k];

                g_0_z_yzz_yyzz[k] = -g_0_z_zz_yyzz[k] * ab_y + g_0_z_zz_yyyzz[k];

                g_0_z_yzz_yzzz[k] = -g_0_z_zz_yzzz[k] * ab_y + g_0_z_zz_yyzzz[k];

                g_0_z_yzz_zzzz[k] = -g_0_z_zz_zzzz[k] * ab_y + g_0_z_zz_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzz_xxxx = cbuffer.data(fg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_zzz_xxxy = cbuffer.data(fg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_zzz_xxxz = cbuffer.data(fg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_zzz_xxyy = cbuffer.data(fg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_zzz_xxyz = cbuffer.data(fg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_zzz_xxzz = cbuffer.data(fg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_zzz_xyyy = cbuffer.data(fg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_zzz_xyyz = cbuffer.data(fg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_zzz_xyzz = cbuffer.data(fg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_zzz_xzzz = cbuffer.data(fg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_zzz_yyyy = cbuffer.data(fg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_zzz_yyyz = cbuffer.data(fg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_zzz_yyzz = cbuffer.data(fg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_zzz_yzzz = cbuffer.data(fg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_zzz_zzzz = cbuffer.data(fg_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxxx, g_0_z_zz_xxxxz, g_0_z_zz_xxxy, g_0_z_zz_xxxyz, g_0_z_zz_xxxz, g_0_z_zz_xxxzz, g_0_z_zz_xxyy, g_0_z_zz_xxyyz, g_0_z_zz_xxyz, g_0_z_zz_xxyzz, g_0_z_zz_xxzz, g_0_z_zz_xxzzz, g_0_z_zz_xyyy, g_0_z_zz_xyyyz, g_0_z_zz_xyyz, g_0_z_zz_xyyzz, g_0_z_zz_xyzz, g_0_z_zz_xyzzz, g_0_z_zz_xzzz, g_0_z_zz_xzzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyyz, g_0_z_zz_yyyz, g_0_z_zz_yyyzz, g_0_z_zz_yyzz, g_0_z_zz_yyzzz, g_0_z_zz_yzzz, g_0_z_zz_yzzzz, g_0_z_zz_zzzz, g_0_z_zz_zzzzz, g_0_z_zzz_xxxx, g_0_z_zzz_xxxy, g_0_z_zzz_xxxz, g_0_z_zzz_xxyy, g_0_z_zzz_xxyz, g_0_z_zzz_xxzz, g_0_z_zzz_xyyy, g_0_z_zzz_xyyz, g_0_z_zzz_xyzz, g_0_z_zzz_xzzz, g_0_z_zzz_yyyy, g_0_z_zzz_yyyz, g_0_z_zzz_yyzz, g_0_z_zzz_yzzz, g_0_z_zzz_zzzz, g_zz_xxxx, g_zz_xxxy, g_zz_xxxz, g_zz_xxyy, g_zz_xxyz, g_zz_xxzz, g_zz_xyyy, g_zz_xyyz, g_zz_xyzz, g_zz_xzzz, g_zz_yyyy, g_zz_yyyz, g_zz_yyzz, g_zz_yzzz, g_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzz_xxxx[k] = g_zz_xxxx[k] - g_0_z_zz_xxxx[k] * ab_z + g_0_z_zz_xxxxz[k];

                g_0_z_zzz_xxxy[k] = g_zz_xxxy[k] - g_0_z_zz_xxxy[k] * ab_z + g_0_z_zz_xxxyz[k];

                g_0_z_zzz_xxxz[k] = g_zz_xxxz[k] - g_0_z_zz_xxxz[k] * ab_z + g_0_z_zz_xxxzz[k];

                g_0_z_zzz_xxyy[k] = g_zz_xxyy[k] - g_0_z_zz_xxyy[k] * ab_z + g_0_z_zz_xxyyz[k];

                g_0_z_zzz_xxyz[k] = g_zz_xxyz[k] - g_0_z_zz_xxyz[k] * ab_z + g_0_z_zz_xxyzz[k];

                g_0_z_zzz_xxzz[k] = g_zz_xxzz[k] - g_0_z_zz_xxzz[k] * ab_z + g_0_z_zz_xxzzz[k];

                g_0_z_zzz_xyyy[k] = g_zz_xyyy[k] - g_0_z_zz_xyyy[k] * ab_z + g_0_z_zz_xyyyz[k];

                g_0_z_zzz_xyyz[k] = g_zz_xyyz[k] - g_0_z_zz_xyyz[k] * ab_z + g_0_z_zz_xyyzz[k];

                g_0_z_zzz_xyzz[k] = g_zz_xyzz[k] - g_0_z_zz_xyzz[k] * ab_z + g_0_z_zz_xyzzz[k];

                g_0_z_zzz_xzzz[k] = g_zz_xzzz[k] - g_0_z_zz_xzzz[k] * ab_z + g_0_z_zz_xzzzz[k];

                g_0_z_zzz_yyyy[k] = g_zz_yyyy[k] - g_0_z_zz_yyyy[k] * ab_z + g_0_z_zz_yyyyz[k];

                g_0_z_zzz_yyyz[k] = g_zz_yyyz[k] - g_0_z_zz_yyyz[k] * ab_z + g_0_z_zz_yyyzz[k];

                g_0_z_zzz_yyzz[k] = g_zz_yyzz[k] - g_0_z_zz_yyzz[k] * ab_z + g_0_z_zz_yyzzz[k];

                g_0_z_zzz_yzzz[k] = g_zz_yzzz[k] - g_0_z_zz_yzzz[k] * ab_z + g_0_z_zz_yzzzz[k];

                g_0_z_zzz_zzzz[k] = g_zz_zzzz[k] - g_0_z_zz_zzzz[k] * ab_z + g_0_z_zz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

