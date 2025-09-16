#include "ElectronRepulsionGeom2000ContrRecDGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_dgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dgxx,
                                            const size_t idx_geom_10_pgxx,
                                            const size_t idx_geom_20_pgxx,
                                            const size_t idx_geom_20_phxx,
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
            /// Set up components of auxilary buffer : PGSS

            const auto pg_geom_10_off = idx_geom_10_pgxx + i * dcomps + j;

            auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_y_xxxx = cbuffer.data(pg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_xxxz = cbuffer.data(pg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_y_xxzz = cbuffer.data(pg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_xzzz = cbuffer.data(pg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_zzzz = cbuffer.data(pg_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_z_xxxx = cbuffer.data(pg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_x_yyyy = cbuffer.data(pg_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_x_yyyz = cbuffer.data(pg_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_x_yyzz = cbuffer.data(pg_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_x_yzzz = cbuffer.data(pg_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_x_zzzz = cbuffer.data(pg_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_z_yyyy = cbuffer.data(pg_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 89 * ccomps * dcomps);

            auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 95 * ccomps * dcomps);

            auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_x_yyyy = cbuffer.data(pg_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_x_yyyz = cbuffer.data(pg_geom_10_off + 101 * ccomps * dcomps);

            auto g_z_0_x_yyzz = cbuffer.data(pg_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_x_yzzz = cbuffer.data(pg_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_x_zzzz = cbuffer.data(pg_geom_10_off + 104 * ccomps * dcomps);

            auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 107 * ccomps * dcomps);

            auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 108 * ccomps * dcomps);

            auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 109 * ccomps * dcomps);

            auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 110 * ccomps * dcomps);

            auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 111 * ccomps * dcomps);

            auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 112 * ccomps * dcomps);

            auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 113 * ccomps * dcomps);

            auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 114 * ccomps * dcomps);

            auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 115 * ccomps * dcomps);

            auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 116 * ccomps * dcomps);

            auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 117 * ccomps * dcomps);

            auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 118 * ccomps * dcomps);

            auto g_z_0_y_zzzz = cbuffer.data(pg_geom_10_off + 119 * ccomps * dcomps);

            auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 134 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PGSS

            const auto pg_geom_20_off = idx_geom_20_pgxx + i * dcomps + j;

            auto g_xx_0_x_xxxx = cbuffer.data(pg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxy = cbuffer.data(pg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxz = cbuffer.data(pg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxyy = cbuffer.data(pg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxyz = cbuffer.data(pg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxzz = cbuffer.data(pg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xyyy = cbuffer.data(pg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xyyz = cbuffer.data(pg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xyzz = cbuffer.data(pg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xzzz = cbuffer.data(pg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_yyyy = cbuffer.data(pg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_yyyz = cbuffer.data(pg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_yyzz = cbuffer.data(pg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_yzzz = cbuffer.data(pg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_zzzz = cbuffer.data(pg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_y_xxxx = cbuffer.data(pg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_y_xxxy = cbuffer.data(pg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_y_xxxz = cbuffer.data(pg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_y_xxyy = cbuffer.data(pg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_y_xxyz = cbuffer.data(pg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_y_xxzz = cbuffer.data(pg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_y_xyyy = cbuffer.data(pg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_y_xyyz = cbuffer.data(pg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_y_xyzz = cbuffer.data(pg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_y_xzzz = cbuffer.data(pg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_y_yyyy = cbuffer.data(pg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_y_yyyz = cbuffer.data(pg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_y_yyzz = cbuffer.data(pg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_y_yzzz = cbuffer.data(pg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_y_zzzz = cbuffer.data(pg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_z_xxxx = cbuffer.data(pg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_z_xxxy = cbuffer.data(pg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_z_xxxz = cbuffer.data(pg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_z_xxyy = cbuffer.data(pg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_z_xxyz = cbuffer.data(pg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_z_xxzz = cbuffer.data(pg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_z_xyyy = cbuffer.data(pg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_z_xyyz = cbuffer.data(pg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_z_xyzz = cbuffer.data(pg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_z_xzzz = cbuffer.data(pg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_z_yyyy = cbuffer.data(pg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_z_yyyz = cbuffer.data(pg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_z_yyzz = cbuffer.data(pg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_z_yzzz = cbuffer.data(pg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_z_zzzz = cbuffer.data(pg_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_x_xxxx = cbuffer.data(pg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_x_xxxy = cbuffer.data(pg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_x_xxxz = cbuffer.data(pg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_x_xxyy = cbuffer.data(pg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_x_xxyz = cbuffer.data(pg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_x_xxzz = cbuffer.data(pg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_x_xyyy = cbuffer.data(pg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_x_xyyz = cbuffer.data(pg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_x_xyzz = cbuffer.data(pg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_x_xzzz = cbuffer.data(pg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_x_yyyy = cbuffer.data(pg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_x_yyyz = cbuffer.data(pg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_x_yyzz = cbuffer.data(pg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_x_yzzz = cbuffer.data(pg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_x_zzzz = cbuffer.data(pg_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_y_xxxx = cbuffer.data(pg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_y_xxxy = cbuffer.data(pg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_y_xxxz = cbuffer.data(pg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_y_xxyy = cbuffer.data(pg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_y_xxyz = cbuffer.data(pg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_y_xxzz = cbuffer.data(pg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_y_xyyy = cbuffer.data(pg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_y_xyyz = cbuffer.data(pg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_y_xyzz = cbuffer.data(pg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_y_xzzz = cbuffer.data(pg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_y_yyyy = cbuffer.data(pg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_y_yyyz = cbuffer.data(pg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_y_yyzz = cbuffer.data(pg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_y_yzzz = cbuffer.data(pg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_y_zzzz = cbuffer.data(pg_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_z_xxxx = cbuffer.data(pg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_z_xxxy = cbuffer.data(pg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_z_xxxz = cbuffer.data(pg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_z_xxyy = cbuffer.data(pg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_z_xxyz = cbuffer.data(pg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_z_xxzz = cbuffer.data(pg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_z_xyyy = cbuffer.data(pg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_z_xyyz = cbuffer.data(pg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_z_xyzz = cbuffer.data(pg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_z_xzzz = cbuffer.data(pg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_z_yyyy = cbuffer.data(pg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_z_yyyz = cbuffer.data(pg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_z_yyzz = cbuffer.data(pg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_z_yzzz = cbuffer.data(pg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_z_zzzz = cbuffer.data(pg_geom_20_off + 89 * ccomps * dcomps);

            auto g_xz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 104 * ccomps * dcomps);

            auto g_xz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 119 * ccomps * dcomps);

            auto g_xz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 134 * ccomps * dcomps);

            auto g_yy_0_x_xxxx = cbuffer.data(pg_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_x_xxxy = cbuffer.data(pg_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_x_xxxz = cbuffer.data(pg_geom_20_off + 137 * ccomps * dcomps);

            auto g_yy_0_x_xxyy = cbuffer.data(pg_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_x_xxyz = cbuffer.data(pg_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_x_xxzz = cbuffer.data(pg_geom_20_off + 140 * ccomps * dcomps);

            auto g_yy_0_x_xyyy = cbuffer.data(pg_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_x_xyyz = cbuffer.data(pg_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_x_xyzz = cbuffer.data(pg_geom_20_off + 143 * ccomps * dcomps);

            auto g_yy_0_x_xzzz = cbuffer.data(pg_geom_20_off + 144 * ccomps * dcomps);

            auto g_yy_0_x_yyyy = cbuffer.data(pg_geom_20_off + 145 * ccomps * dcomps);

            auto g_yy_0_x_yyyz = cbuffer.data(pg_geom_20_off + 146 * ccomps * dcomps);

            auto g_yy_0_x_yyzz = cbuffer.data(pg_geom_20_off + 147 * ccomps * dcomps);

            auto g_yy_0_x_yzzz = cbuffer.data(pg_geom_20_off + 148 * ccomps * dcomps);

            auto g_yy_0_x_zzzz = cbuffer.data(pg_geom_20_off + 149 * ccomps * dcomps);

            auto g_yy_0_y_xxxx = cbuffer.data(pg_geom_20_off + 150 * ccomps * dcomps);

            auto g_yy_0_y_xxxy = cbuffer.data(pg_geom_20_off + 151 * ccomps * dcomps);

            auto g_yy_0_y_xxxz = cbuffer.data(pg_geom_20_off + 152 * ccomps * dcomps);

            auto g_yy_0_y_xxyy = cbuffer.data(pg_geom_20_off + 153 * ccomps * dcomps);

            auto g_yy_0_y_xxyz = cbuffer.data(pg_geom_20_off + 154 * ccomps * dcomps);

            auto g_yy_0_y_xxzz = cbuffer.data(pg_geom_20_off + 155 * ccomps * dcomps);

            auto g_yy_0_y_xyyy = cbuffer.data(pg_geom_20_off + 156 * ccomps * dcomps);

            auto g_yy_0_y_xyyz = cbuffer.data(pg_geom_20_off + 157 * ccomps * dcomps);

            auto g_yy_0_y_xyzz = cbuffer.data(pg_geom_20_off + 158 * ccomps * dcomps);

            auto g_yy_0_y_xzzz = cbuffer.data(pg_geom_20_off + 159 * ccomps * dcomps);

            auto g_yy_0_y_yyyy = cbuffer.data(pg_geom_20_off + 160 * ccomps * dcomps);

            auto g_yy_0_y_yyyz = cbuffer.data(pg_geom_20_off + 161 * ccomps * dcomps);

            auto g_yy_0_y_yyzz = cbuffer.data(pg_geom_20_off + 162 * ccomps * dcomps);

            auto g_yy_0_y_yzzz = cbuffer.data(pg_geom_20_off + 163 * ccomps * dcomps);

            auto g_yy_0_y_zzzz = cbuffer.data(pg_geom_20_off + 164 * ccomps * dcomps);

            auto g_yy_0_z_xxxx = cbuffer.data(pg_geom_20_off + 165 * ccomps * dcomps);

            auto g_yy_0_z_xxxy = cbuffer.data(pg_geom_20_off + 166 * ccomps * dcomps);

            auto g_yy_0_z_xxxz = cbuffer.data(pg_geom_20_off + 167 * ccomps * dcomps);

            auto g_yy_0_z_xxyy = cbuffer.data(pg_geom_20_off + 168 * ccomps * dcomps);

            auto g_yy_0_z_xxyz = cbuffer.data(pg_geom_20_off + 169 * ccomps * dcomps);

            auto g_yy_0_z_xxzz = cbuffer.data(pg_geom_20_off + 170 * ccomps * dcomps);

            auto g_yy_0_z_xyyy = cbuffer.data(pg_geom_20_off + 171 * ccomps * dcomps);

            auto g_yy_0_z_xyyz = cbuffer.data(pg_geom_20_off + 172 * ccomps * dcomps);

            auto g_yy_0_z_xyzz = cbuffer.data(pg_geom_20_off + 173 * ccomps * dcomps);

            auto g_yy_0_z_xzzz = cbuffer.data(pg_geom_20_off + 174 * ccomps * dcomps);

            auto g_yy_0_z_yyyy = cbuffer.data(pg_geom_20_off + 175 * ccomps * dcomps);

            auto g_yy_0_z_yyyz = cbuffer.data(pg_geom_20_off + 176 * ccomps * dcomps);

            auto g_yy_0_z_yyzz = cbuffer.data(pg_geom_20_off + 177 * ccomps * dcomps);

            auto g_yy_0_z_yzzz = cbuffer.data(pg_geom_20_off + 178 * ccomps * dcomps);

            auto g_yy_0_z_zzzz = cbuffer.data(pg_geom_20_off + 179 * ccomps * dcomps);

            auto g_yz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 180 * ccomps * dcomps);

            auto g_yz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 181 * ccomps * dcomps);

            auto g_yz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 182 * ccomps * dcomps);

            auto g_yz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 183 * ccomps * dcomps);

            auto g_yz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 184 * ccomps * dcomps);

            auto g_yz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 185 * ccomps * dcomps);

            auto g_yz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 186 * ccomps * dcomps);

            auto g_yz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 187 * ccomps * dcomps);

            auto g_yz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 188 * ccomps * dcomps);

            auto g_yz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 189 * ccomps * dcomps);

            auto g_yz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 190 * ccomps * dcomps);

            auto g_yz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 191 * ccomps * dcomps);

            auto g_yz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 192 * ccomps * dcomps);

            auto g_yz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 193 * ccomps * dcomps);

            auto g_yz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 194 * ccomps * dcomps);

            auto g_yz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 195 * ccomps * dcomps);

            auto g_yz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 196 * ccomps * dcomps);

            auto g_yz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 197 * ccomps * dcomps);

            auto g_yz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 198 * ccomps * dcomps);

            auto g_yz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 199 * ccomps * dcomps);

            auto g_yz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 200 * ccomps * dcomps);

            auto g_yz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 201 * ccomps * dcomps);

            auto g_yz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 202 * ccomps * dcomps);

            auto g_yz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 203 * ccomps * dcomps);

            auto g_yz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 204 * ccomps * dcomps);

            auto g_yz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 205 * ccomps * dcomps);

            auto g_yz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 206 * ccomps * dcomps);

            auto g_yz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 207 * ccomps * dcomps);

            auto g_yz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 208 * ccomps * dcomps);

            auto g_yz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 209 * ccomps * dcomps);

            auto g_yz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 210 * ccomps * dcomps);

            auto g_yz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 211 * ccomps * dcomps);

            auto g_yz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 212 * ccomps * dcomps);

            auto g_yz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 213 * ccomps * dcomps);

            auto g_yz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 214 * ccomps * dcomps);

            auto g_yz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 215 * ccomps * dcomps);

            auto g_yz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 216 * ccomps * dcomps);

            auto g_yz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 217 * ccomps * dcomps);

            auto g_yz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 218 * ccomps * dcomps);

            auto g_yz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 219 * ccomps * dcomps);

            auto g_yz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 220 * ccomps * dcomps);

            auto g_yz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 221 * ccomps * dcomps);

            auto g_yz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 222 * ccomps * dcomps);

            auto g_yz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 223 * ccomps * dcomps);

            auto g_yz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 224 * ccomps * dcomps);

            auto g_zz_0_x_xxxx = cbuffer.data(pg_geom_20_off + 225 * ccomps * dcomps);

            auto g_zz_0_x_xxxy = cbuffer.data(pg_geom_20_off + 226 * ccomps * dcomps);

            auto g_zz_0_x_xxxz = cbuffer.data(pg_geom_20_off + 227 * ccomps * dcomps);

            auto g_zz_0_x_xxyy = cbuffer.data(pg_geom_20_off + 228 * ccomps * dcomps);

            auto g_zz_0_x_xxyz = cbuffer.data(pg_geom_20_off + 229 * ccomps * dcomps);

            auto g_zz_0_x_xxzz = cbuffer.data(pg_geom_20_off + 230 * ccomps * dcomps);

            auto g_zz_0_x_xyyy = cbuffer.data(pg_geom_20_off + 231 * ccomps * dcomps);

            auto g_zz_0_x_xyyz = cbuffer.data(pg_geom_20_off + 232 * ccomps * dcomps);

            auto g_zz_0_x_xyzz = cbuffer.data(pg_geom_20_off + 233 * ccomps * dcomps);

            auto g_zz_0_x_xzzz = cbuffer.data(pg_geom_20_off + 234 * ccomps * dcomps);

            auto g_zz_0_x_yyyy = cbuffer.data(pg_geom_20_off + 235 * ccomps * dcomps);

            auto g_zz_0_x_yyyz = cbuffer.data(pg_geom_20_off + 236 * ccomps * dcomps);

            auto g_zz_0_x_yyzz = cbuffer.data(pg_geom_20_off + 237 * ccomps * dcomps);

            auto g_zz_0_x_yzzz = cbuffer.data(pg_geom_20_off + 238 * ccomps * dcomps);

            auto g_zz_0_x_zzzz = cbuffer.data(pg_geom_20_off + 239 * ccomps * dcomps);

            auto g_zz_0_y_xxxx = cbuffer.data(pg_geom_20_off + 240 * ccomps * dcomps);

            auto g_zz_0_y_xxxy = cbuffer.data(pg_geom_20_off + 241 * ccomps * dcomps);

            auto g_zz_0_y_xxxz = cbuffer.data(pg_geom_20_off + 242 * ccomps * dcomps);

            auto g_zz_0_y_xxyy = cbuffer.data(pg_geom_20_off + 243 * ccomps * dcomps);

            auto g_zz_0_y_xxyz = cbuffer.data(pg_geom_20_off + 244 * ccomps * dcomps);

            auto g_zz_0_y_xxzz = cbuffer.data(pg_geom_20_off + 245 * ccomps * dcomps);

            auto g_zz_0_y_xyyy = cbuffer.data(pg_geom_20_off + 246 * ccomps * dcomps);

            auto g_zz_0_y_xyyz = cbuffer.data(pg_geom_20_off + 247 * ccomps * dcomps);

            auto g_zz_0_y_xyzz = cbuffer.data(pg_geom_20_off + 248 * ccomps * dcomps);

            auto g_zz_0_y_xzzz = cbuffer.data(pg_geom_20_off + 249 * ccomps * dcomps);

            auto g_zz_0_y_yyyy = cbuffer.data(pg_geom_20_off + 250 * ccomps * dcomps);

            auto g_zz_0_y_yyyz = cbuffer.data(pg_geom_20_off + 251 * ccomps * dcomps);

            auto g_zz_0_y_yyzz = cbuffer.data(pg_geom_20_off + 252 * ccomps * dcomps);

            auto g_zz_0_y_yzzz = cbuffer.data(pg_geom_20_off + 253 * ccomps * dcomps);

            auto g_zz_0_y_zzzz = cbuffer.data(pg_geom_20_off + 254 * ccomps * dcomps);

            auto g_zz_0_z_xxxx = cbuffer.data(pg_geom_20_off + 255 * ccomps * dcomps);

            auto g_zz_0_z_xxxy = cbuffer.data(pg_geom_20_off + 256 * ccomps * dcomps);

            auto g_zz_0_z_xxxz = cbuffer.data(pg_geom_20_off + 257 * ccomps * dcomps);

            auto g_zz_0_z_xxyy = cbuffer.data(pg_geom_20_off + 258 * ccomps * dcomps);

            auto g_zz_0_z_xxyz = cbuffer.data(pg_geom_20_off + 259 * ccomps * dcomps);

            auto g_zz_0_z_xxzz = cbuffer.data(pg_geom_20_off + 260 * ccomps * dcomps);

            auto g_zz_0_z_xyyy = cbuffer.data(pg_geom_20_off + 261 * ccomps * dcomps);

            auto g_zz_0_z_xyyz = cbuffer.data(pg_geom_20_off + 262 * ccomps * dcomps);

            auto g_zz_0_z_xyzz = cbuffer.data(pg_geom_20_off + 263 * ccomps * dcomps);

            auto g_zz_0_z_xzzz = cbuffer.data(pg_geom_20_off + 264 * ccomps * dcomps);

            auto g_zz_0_z_yyyy = cbuffer.data(pg_geom_20_off + 265 * ccomps * dcomps);

            auto g_zz_0_z_yyyz = cbuffer.data(pg_geom_20_off + 266 * ccomps * dcomps);

            auto g_zz_0_z_yyzz = cbuffer.data(pg_geom_20_off + 267 * ccomps * dcomps);

            auto g_zz_0_z_yzzz = cbuffer.data(pg_geom_20_off + 268 * ccomps * dcomps);

            auto g_zz_0_z_zzzz = cbuffer.data(pg_geom_20_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_20_off = idx_geom_20_phxx + i * dcomps + j;

            auto g_xx_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 149 * ccomps * dcomps);

            auto g_xz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 209 * ccomps * dcomps);

            auto g_yy_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 239 * ccomps * dcomps);

            auto g_yy_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 240 * ccomps * dcomps);

            auto g_yy_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 241 * ccomps * dcomps);

            auto g_yy_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 242 * ccomps * dcomps);

            auto g_yy_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 243 * ccomps * dcomps);

            auto g_yy_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 244 * ccomps * dcomps);

            auto g_yy_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 245 * ccomps * dcomps);

            auto g_yy_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 246 * ccomps * dcomps);

            auto g_yy_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 247 * ccomps * dcomps);

            auto g_yy_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 248 * ccomps * dcomps);

            auto g_yy_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 249 * ccomps * dcomps);

            auto g_yy_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 250 * ccomps * dcomps);

            auto g_yy_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 269 * ccomps * dcomps);

            auto g_yz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 299 * ccomps * dcomps);

            auto g_yz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 300 * ccomps * dcomps);

            auto g_yz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 301 * ccomps * dcomps);

            auto g_yz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 302 * ccomps * dcomps);

            auto g_yz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 303 * ccomps * dcomps);

            auto g_yz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 304 * ccomps * dcomps);

            auto g_yz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 305 * ccomps * dcomps);

            auto g_yz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 306 * ccomps * dcomps);

            auto g_yz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 307 * ccomps * dcomps);

            auto g_yz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 308 * ccomps * dcomps);

            auto g_yz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 309 * ccomps * dcomps);

            auto g_yz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 310 * ccomps * dcomps);

            auto g_yz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 311 * ccomps * dcomps);

            auto g_yz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 312 * ccomps * dcomps);

            auto g_yz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 313 * ccomps * dcomps);

            auto g_yz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 329 * ccomps * dcomps);

            auto g_zz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 359 * ccomps * dcomps);

            auto g_zz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 360 * ccomps * dcomps);

            auto g_zz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 361 * ccomps * dcomps);

            auto g_zz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 362 * ccomps * dcomps);

            auto g_zz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 363 * ccomps * dcomps);

            auto g_zz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 364 * ccomps * dcomps);

            auto g_zz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 365 * ccomps * dcomps);

            auto g_zz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 366 * ccomps * dcomps);

            auto g_zz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 367 * ccomps * dcomps);

            auto g_zz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 368 * ccomps * dcomps);

            auto g_zz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 369 * ccomps * dcomps);

            auto g_zz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 370 * ccomps * dcomps);

            auto g_zz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 371 * ccomps * dcomps);

            auto g_zz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 372 * ccomps * dcomps);

            auto g_zz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 373 * ccomps * dcomps);

            auto g_zz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 374 * ccomps * dcomps);

            auto g_zz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 375 * ccomps * dcomps);

            auto g_zz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 376 * ccomps * dcomps);

            auto g_zz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 377 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dgxx

            const auto dg_geom_20_off = idx_geom_20_dgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz, g_xx_0_x_xxxx, g_xx_0_x_xxxxx, g_xx_0_x_xxxxy, g_xx_0_x_xxxxz, g_xx_0_x_xxxy, g_xx_0_x_xxxyy, g_xx_0_x_xxxyz, g_xx_0_x_xxxz, g_xx_0_x_xxxzz, g_xx_0_x_xxyy, g_xx_0_x_xxyyy, g_xx_0_x_xxyyz, g_xx_0_x_xxyz, g_xx_0_x_xxyzz, g_xx_0_x_xxzz, g_xx_0_x_xxzzz, g_xx_0_x_xyyy, g_xx_0_x_xyyyy, g_xx_0_x_xyyyz, g_xx_0_x_xyyz, g_xx_0_x_xyyzz, g_xx_0_x_xyzz, g_xx_0_x_xyzzz, g_xx_0_x_xzzz, g_xx_0_x_xzzzz, g_xx_0_x_yyyy, g_xx_0_x_yyyz, g_xx_0_x_yyzz, g_xx_0_x_yzzz, g_xx_0_x_zzzz, g_xx_0_xx_xxxx, g_xx_0_xx_xxxy, g_xx_0_xx_xxxz, g_xx_0_xx_xxyy, g_xx_0_xx_xxyz, g_xx_0_xx_xxzz, g_xx_0_xx_xyyy, g_xx_0_xx_xyyz, g_xx_0_xx_xyzz, g_xx_0_xx_xzzz, g_xx_0_xx_yyyy, g_xx_0_xx_yyyz, g_xx_0_xx_yyzz, g_xx_0_xx_yzzz, g_xx_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xx_xxxx[k] = -2.0 * g_x_0_x_xxxx[k] - g_xx_0_x_xxxx[k] * ab_x + g_xx_0_x_xxxxx[k];

                g_xx_0_xx_xxxy[k] = -2.0 * g_x_0_x_xxxy[k] - g_xx_0_x_xxxy[k] * ab_x + g_xx_0_x_xxxxy[k];

                g_xx_0_xx_xxxz[k] = -2.0 * g_x_0_x_xxxz[k] - g_xx_0_x_xxxz[k] * ab_x + g_xx_0_x_xxxxz[k];

                g_xx_0_xx_xxyy[k] = -2.0 * g_x_0_x_xxyy[k] - g_xx_0_x_xxyy[k] * ab_x + g_xx_0_x_xxxyy[k];

                g_xx_0_xx_xxyz[k] = -2.0 * g_x_0_x_xxyz[k] - g_xx_0_x_xxyz[k] * ab_x + g_xx_0_x_xxxyz[k];

                g_xx_0_xx_xxzz[k] = -2.0 * g_x_0_x_xxzz[k] - g_xx_0_x_xxzz[k] * ab_x + g_xx_0_x_xxxzz[k];

                g_xx_0_xx_xyyy[k] = -2.0 * g_x_0_x_xyyy[k] - g_xx_0_x_xyyy[k] * ab_x + g_xx_0_x_xxyyy[k];

                g_xx_0_xx_xyyz[k] = -2.0 * g_x_0_x_xyyz[k] - g_xx_0_x_xyyz[k] * ab_x + g_xx_0_x_xxyyz[k];

                g_xx_0_xx_xyzz[k] = -2.0 * g_x_0_x_xyzz[k] - g_xx_0_x_xyzz[k] * ab_x + g_xx_0_x_xxyzz[k];

                g_xx_0_xx_xzzz[k] = -2.0 * g_x_0_x_xzzz[k] - g_xx_0_x_xzzz[k] * ab_x + g_xx_0_x_xxzzz[k];

                g_xx_0_xx_yyyy[k] = -2.0 * g_x_0_x_yyyy[k] - g_xx_0_x_yyyy[k] * ab_x + g_xx_0_x_xyyyy[k];

                g_xx_0_xx_yyyz[k] = -2.0 * g_x_0_x_yyyz[k] - g_xx_0_x_yyyz[k] * ab_x + g_xx_0_x_xyyyz[k];

                g_xx_0_xx_yyzz[k] = -2.0 * g_x_0_x_yyzz[k] - g_xx_0_x_yyzz[k] * ab_x + g_xx_0_x_xyyzz[k];

                g_xx_0_xx_yzzz[k] = -2.0 * g_x_0_x_yzzz[k] - g_xx_0_x_yzzz[k] * ab_x + g_xx_0_x_xyzzz[k];

                g_xx_0_xx_zzzz[k] = -2.0 * g_x_0_x_zzzz[k] - g_xx_0_x_zzzz[k] * ab_x + g_xx_0_x_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxx, g_xx_0_x_xxxxy, g_xx_0_x_xxxy, g_xx_0_x_xxxyy, g_xx_0_x_xxxyz, g_xx_0_x_xxxz, g_xx_0_x_xxyy, g_xx_0_x_xxyyy, g_xx_0_x_xxyyz, g_xx_0_x_xxyz, g_xx_0_x_xxyzz, g_xx_0_x_xxzz, g_xx_0_x_xyyy, g_xx_0_x_xyyyy, g_xx_0_x_xyyyz, g_xx_0_x_xyyz, g_xx_0_x_xyyzz, g_xx_0_x_xyzz, g_xx_0_x_xyzzz, g_xx_0_x_xzzz, g_xx_0_x_yyyy, g_xx_0_x_yyyyy, g_xx_0_x_yyyyz, g_xx_0_x_yyyz, g_xx_0_x_yyyzz, g_xx_0_x_yyzz, g_xx_0_x_yyzzz, g_xx_0_x_yzzz, g_xx_0_x_yzzzz, g_xx_0_x_zzzz, g_xx_0_xy_xxxx, g_xx_0_xy_xxxy, g_xx_0_xy_xxxz, g_xx_0_xy_xxyy, g_xx_0_xy_xxyz, g_xx_0_xy_xxzz, g_xx_0_xy_xyyy, g_xx_0_xy_xyyz, g_xx_0_xy_xyzz, g_xx_0_xy_xzzz, g_xx_0_xy_yyyy, g_xx_0_xy_yyyz, g_xx_0_xy_yyzz, g_xx_0_xy_yzzz, g_xx_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xy_xxxx[k] = -g_xx_0_x_xxxx[k] * ab_y + g_xx_0_x_xxxxy[k];

                g_xx_0_xy_xxxy[k] = -g_xx_0_x_xxxy[k] * ab_y + g_xx_0_x_xxxyy[k];

                g_xx_0_xy_xxxz[k] = -g_xx_0_x_xxxz[k] * ab_y + g_xx_0_x_xxxyz[k];

                g_xx_0_xy_xxyy[k] = -g_xx_0_x_xxyy[k] * ab_y + g_xx_0_x_xxyyy[k];

                g_xx_0_xy_xxyz[k] = -g_xx_0_x_xxyz[k] * ab_y + g_xx_0_x_xxyyz[k];

                g_xx_0_xy_xxzz[k] = -g_xx_0_x_xxzz[k] * ab_y + g_xx_0_x_xxyzz[k];

                g_xx_0_xy_xyyy[k] = -g_xx_0_x_xyyy[k] * ab_y + g_xx_0_x_xyyyy[k];

                g_xx_0_xy_xyyz[k] = -g_xx_0_x_xyyz[k] * ab_y + g_xx_0_x_xyyyz[k];

                g_xx_0_xy_xyzz[k] = -g_xx_0_x_xyzz[k] * ab_y + g_xx_0_x_xyyzz[k];

                g_xx_0_xy_xzzz[k] = -g_xx_0_x_xzzz[k] * ab_y + g_xx_0_x_xyzzz[k];

                g_xx_0_xy_yyyy[k] = -g_xx_0_x_yyyy[k] * ab_y + g_xx_0_x_yyyyy[k];

                g_xx_0_xy_yyyz[k] = -g_xx_0_x_yyyz[k] * ab_y + g_xx_0_x_yyyyz[k];

                g_xx_0_xy_yyzz[k] = -g_xx_0_x_yyzz[k] * ab_y + g_xx_0_x_yyyzz[k];

                g_xx_0_xy_yzzz[k] = -g_xx_0_x_yzzz[k] * ab_y + g_xx_0_x_yyzzz[k];

                g_xx_0_xy_zzzz[k] = -g_xx_0_x_zzzz[k] * ab_y + g_xx_0_x_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxx, g_xx_0_x_xxxxz, g_xx_0_x_xxxy, g_xx_0_x_xxxyz, g_xx_0_x_xxxz, g_xx_0_x_xxxzz, g_xx_0_x_xxyy, g_xx_0_x_xxyyz, g_xx_0_x_xxyz, g_xx_0_x_xxyzz, g_xx_0_x_xxzz, g_xx_0_x_xxzzz, g_xx_0_x_xyyy, g_xx_0_x_xyyyz, g_xx_0_x_xyyz, g_xx_0_x_xyyzz, g_xx_0_x_xyzz, g_xx_0_x_xyzzz, g_xx_0_x_xzzz, g_xx_0_x_xzzzz, g_xx_0_x_yyyy, g_xx_0_x_yyyyz, g_xx_0_x_yyyz, g_xx_0_x_yyyzz, g_xx_0_x_yyzz, g_xx_0_x_yyzzz, g_xx_0_x_yzzz, g_xx_0_x_yzzzz, g_xx_0_x_zzzz, g_xx_0_x_zzzzz, g_xx_0_xz_xxxx, g_xx_0_xz_xxxy, g_xx_0_xz_xxxz, g_xx_0_xz_xxyy, g_xx_0_xz_xxyz, g_xx_0_xz_xxzz, g_xx_0_xz_xyyy, g_xx_0_xz_xyyz, g_xx_0_xz_xyzz, g_xx_0_xz_xzzz, g_xx_0_xz_yyyy, g_xx_0_xz_yyyz, g_xx_0_xz_yyzz, g_xx_0_xz_yzzz, g_xx_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xz_xxxx[k] = -g_xx_0_x_xxxx[k] * ab_z + g_xx_0_x_xxxxz[k];

                g_xx_0_xz_xxxy[k] = -g_xx_0_x_xxxy[k] * ab_z + g_xx_0_x_xxxyz[k];

                g_xx_0_xz_xxxz[k] = -g_xx_0_x_xxxz[k] * ab_z + g_xx_0_x_xxxzz[k];

                g_xx_0_xz_xxyy[k] = -g_xx_0_x_xxyy[k] * ab_z + g_xx_0_x_xxyyz[k];

                g_xx_0_xz_xxyz[k] = -g_xx_0_x_xxyz[k] * ab_z + g_xx_0_x_xxyzz[k];

                g_xx_0_xz_xxzz[k] = -g_xx_0_x_xxzz[k] * ab_z + g_xx_0_x_xxzzz[k];

                g_xx_0_xz_xyyy[k] = -g_xx_0_x_xyyy[k] * ab_z + g_xx_0_x_xyyyz[k];

                g_xx_0_xz_xyyz[k] = -g_xx_0_x_xyyz[k] * ab_z + g_xx_0_x_xyyzz[k];

                g_xx_0_xz_xyzz[k] = -g_xx_0_x_xyzz[k] * ab_z + g_xx_0_x_xyzzz[k];

                g_xx_0_xz_xzzz[k] = -g_xx_0_x_xzzz[k] * ab_z + g_xx_0_x_xzzzz[k];

                g_xx_0_xz_yyyy[k] = -g_xx_0_x_yyyy[k] * ab_z + g_xx_0_x_yyyyz[k];

                g_xx_0_xz_yyyz[k] = -g_xx_0_x_yyyz[k] * ab_z + g_xx_0_x_yyyzz[k];

                g_xx_0_xz_yyzz[k] = -g_xx_0_x_yyzz[k] * ab_z + g_xx_0_x_yyzzz[k];

                g_xx_0_xz_yzzz[k] = -g_xx_0_x_yzzz[k] * ab_z + g_xx_0_x_yzzzz[k];

                g_xx_0_xz_zzzz[k] = -g_xx_0_x_zzzz[k] * ab_z + g_xx_0_x_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_y_xxxx, g_xx_0_y_xxxxy, g_xx_0_y_xxxy, g_xx_0_y_xxxyy, g_xx_0_y_xxxyz, g_xx_0_y_xxxz, g_xx_0_y_xxyy, g_xx_0_y_xxyyy, g_xx_0_y_xxyyz, g_xx_0_y_xxyz, g_xx_0_y_xxyzz, g_xx_0_y_xxzz, g_xx_0_y_xyyy, g_xx_0_y_xyyyy, g_xx_0_y_xyyyz, g_xx_0_y_xyyz, g_xx_0_y_xyyzz, g_xx_0_y_xyzz, g_xx_0_y_xyzzz, g_xx_0_y_xzzz, g_xx_0_y_yyyy, g_xx_0_y_yyyyy, g_xx_0_y_yyyyz, g_xx_0_y_yyyz, g_xx_0_y_yyyzz, g_xx_0_y_yyzz, g_xx_0_y_yyzzz, g_xx_0_y_yzzz, g_xx_0_y_yzzzz, g_xx_0_y_zzzz, g_xx_0_yy_xxxx, g_xx_0_yy_xxxy, g_xx_0_yy_xxxz, g_xx_0_yy_xxyy, g_xx_0_yy_xxyz, g_xx_0_yy_xxzz, g_xx_0_yy_xyyy, g_xx_0_yy_xyyz, g_xx_0_yy_xyzz, g_xx_0_yy_xzzz, g_xx_0_yy_yyyy, g_xx_0_yy_yyyz, g_xx_0_yy_yyzz, g_xx_0_yy_yzzz, g_xx_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yy_xxxx[k] = -g_xx_0_y_xxxx[k] * ab_y + g_xx_0_y_xxxxy[k];

                g_xx_0_yy_xxxy[k] = -g_xx_0_y_xxxy[k] * ab_y + g_xx_0_y_xxxyy[k];

                g_xx_0_yy_xxxz[k] = -g_xx_0_y_xxxz[k] * ab_y + g_xx_0_y_xxxyz[k];

                g_xx_0_yy_xxyy[k] = -g_xx_0_y_xxyy[k] * ab_y + g_xx_0_y_xxyyy[k];

                g_xx_0_yy_xxyz[k] = -g_xx_0_y_xxyz[k] * ab_y + g_xx_0_y_xxyyz[k];

                g_xx_0_yy_xxzz[k] = -g_xx_0_y_xxzz[k] * ab_y + g_xx_0_y_xxyzz[k];

                g_xx_0_yy_xyyy[k] = -g_xx_0_y_xyyy[k] * ab_y + g_xx_0_y_xyyyy[k];

                g_xx_0_yy_xyyz[k] = -g_xx_0_y_xyyz[k] * ab_y + g_xx_0_y_xyyyz[k];

                g_xx_0_yy_xyzz[k] = -g_xx_0_y_xyzz[k] * ab_y + g_xx_0_y_xyyzz[k];

                g_xx_0_yy_xzzz[k] = -g_xx_0_y_xzzz[k] * ab_y + g_xx_0_y_xyzzz[k];

                g_xx_0_yy_yyyy[k] = -g_xx_0_y_yyyy[k] * ab_y + g_xx_0_y_yyyyy[k];

                g_xx_0_yy_yyyz[k] = -g_xx_0_y_yyyz[k] * ab_y + g_xx_0_y_yyyyz[k];

                g_xx_0_yy_yyzz[k] = -g_xx_0_y_yyzz[k] * ab_y + g_xx_0_y_yyyzz[k];

                g_xx_0_yy_yzzz[k] = -g_xx_0_y_yzzz[k] * ab_y + g_xx_0_y_yyzzz[k];

                g_xx_0_yy_zzzz[k] = -g_xx_0_y_zzzz[k] * ab_y + g_xx_0_y_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yz_xxxx, g_xx_0_yz_xxxy, g_xx_0_yz_xxxz, g_xx_0_yz_xxyy, g_xx_0_yz_xxyz, g_xx_0_yz_xxzz, g_xx_0_yz_xyyy, g_xx_0_yz_xyyz, g_xx_0_yz_xyzz, g_xx_0_yz_xzzz, g_xx_0_yz_yyyy, g_xx_0_yz_yyyz, g_xx_0_yz_yyzz, g_xx_0_yz_yzzz, g_xx_0_yz_zzzz, g_xx_0_z_xxxx, g_xx_0_z_xxxxy, g_xx_0_z_xxxy, g_xx_0_z_xxxyy, g_xx_0_z_xxxyz, g_xx_0_z_xxxz, g_xx_0_z_xxyy, g_xx_0_z_xxyyy, g_xx_0_z_xxyyz, g_xx_0_z_xxyz, g_xx_0_z_xxyzz, g_xx_0_z_xxzz, g_xx_0_z_xyyy, g_xx_0_z_xyyyy, g_xx_0_z_xyyyz, g_xx_0_z_xyyz, g_xx_0_z_xyyzz, g_xx_0_z_xyzz, g_xx_0_z_xyzzz, g_xx_0_z_xzzz, g_xx_0_z_yyyy, g_xx_0_z_yyyyy, g_xx_0_z_yyyyz, g_xx_0_z_yyyz, g_xx_0_z_yyyzz, g_xx_0_z_yyzz, g_xx_0_z_yyzzz, g_xx_0_z_yzzz, g_xx_0_z_yzzzz, g_xx_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yz_xxxx[k] = -g_xx_0_z_xxxx[k] * ab_y + g_xx_0_z_xxxxy[k];

                g_xx_0_yz_xxxy[k] = -g_xx_0_z_xxxy[k] * ab_y + g_xx_0_z_xxxyy[k];

                g_xx_0_yz_xxxz[k] = -g_xx_0_z_xxxz[k] * ab_y + g_xx_0_z_xxxyz[k];

                g_xx_0_yz_xxyy[k] = -g_xx_0_z_xxyy[k] * ab_y + g_xx_0_z_xxyyy[k];

                g_xx_0_yz_xxyz[k] = -g_xx_0_z_xxyz[k] * ab_y + g_xx_0_z_xxyyz[k];

                g_xx_0_yz_xxzz[k] = -g_xx_0_z_xxzz[k] * ab_y + g_xx_0_z_xxyzz[k];

                g_xx_0_yz_xyyy[k] = -g_xx_0_z_xyyy[k] * ab_y + g_xx_0_z_xyyyy[k];

                g_xx_0_yz_xyyz[k] = -g_xx_0_z_xyyz[k] * ab_y + g_xx_0_z_xyyyz[k];

                g_xx_0_yz_xyzz[k] = -g_xx_0_z_xyzz[k] * ab_y + g_xx_0_z_xyyzz[k];

                g_xx_0_yz_xzzz[k] = -g_xx_0_z_xzzz[k] * ab_y + g_xx_0_z_xyzzz[k];

                g_xx_0_yz_yyyy[k] = -g_xx_0_z_yyyy[k] * ab_y + g_xx_0_z_yyyyy[k];

                g_xx_0_yz_yyyz[k] = -g_xx_0_z_yyyz[k] * ab_y + g_xx_0_z_yyyyz[k];

                g_xx_0_yz_yyzz[k] = -g_xx_0_z_yyzz[k] * ab_y + g_xx_0_z_yyyzz[k];

                g_xx_0_yz_yzzz[k] = -g_xx_0_z_yzzz[k] * ab_y + g_xx_0_z_yyzzz[k];

                g_xx_0_yz_zzzz[k] = -g_xx_0_z_zzzz[k] * ab_y + g_xx_0_z_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_z_xxxx, g_xx_0_z_xxxxz, g_xx_0_z_xxxy, g_xx_0_z_xxxyz, g_xx_0_z_xxxz, g_xx_0_z_xxxzz, g_xx_0_z_xxyy, g_xx_0_z_xxyyz, g_xx_0_z_xxyz, g_xx_0_z_xxyzz, g_xx_0_z_xxzz, g_xx_0_z_xxzzz, g_xx_0_z_xyyy, g_xx_0_z_xyyyz, g_xx_0_z_xyyz, g_xx_0_z_xyyzz, g_xx_0_z_xyzz, g_xx_0_z_xyzzz, g_xx_0_z_xzzz, g_xx_0_z_xzzzz, g_xx_0_z_yyyy, g_xx_0_z_yyyyz, g_xx_0_z_yyyz, g_xx_0_z_yyyzz, g_xx_0_z_yyzz, g_xx_0_z_yyzzz, g_xx_0_z_yzzz, g_xx_0_z_yzzzz, g_xx_0_z_zzzz, g_xx_0_z_zzzzz, g_xx_0_zz_xxxx, g_xx_0_zz_xxxy, g_xx_0_zz_xxxz, g_xx_0_zz_xxyy, g_xx_0_zz_xxyz, g_xx_0_zz_xxzz, g_xx_0_zz_xyyy, g_xx_0_zz_xyyz, g_xx_0_zz_xyzz, g_xx_0_zz_xzzz, g_xx_0_zz_yyyy, g_xx_0_zz_yyyz, g_xx_0_zz_yyzz, g_xx_0_zz_yzzz, g_xx_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zz_xxxx[k] = -g_xx_0_z_xxxx[k] * ab_z + g_xx_0_z_xxxxz[k];

                g_xx_0_zz_xxxy[k] = -g_xx_0_z_xxxy[k] * ab_z + g_xx_0_z_xxxyz[k];

                g_xx_0_zz_xxxz[k] = -g_xx_0_z_xxxz[k] * ab_z + g_xx_0_z_xxxzz[k];

                g_xx_0_zz_xxyy[k] = -g_xx_0_z_xxyy[k] * ab_z + g_xx_0_z_xxyyz[k];

                g_xx_0_zz_xxyz[k] = -g_xx_0_z_xxyz[k] * ab_z + g_xx_0_z_xxyzz[k];

                g_xx_0_zz_xxzz[k] = -g_xx_0_z_xxzz[k] * ab_z + g_xx_0_z_xxzzz[k];

                g_xx_0_zz_xyyy[k] = -g_xx_0_z_xyyy[k] * ab_z + g_xx_0_z_xyyyz[k];

                g_xx_0_zz_xyyz[k] = -g_xx_0_z_xyyz[k] * ab_z + g_xx_0_z_xyyzz[k];

                g_xx_0_zz_xyzz[k] = -g_xx_0_z_xyzz[k] * ab_z + g_xx_0_z_xyzzz[k];

                g_xx_0_zz_xzzz[k] = -g_xx_0_z_xzzz[k] * ab_z + g_xx_0_z_xzzzz[k];

                g_xx_0_zz_yyyy[k] = -g_xx_0_z_yyyy[k] * ab_z + g_xx_0_z_yyyyz[k];

                g_xx_0_zz_yyyz[k] = -g_xx_0_z_yyyz[k] * ab_z + g_xx_0_z_yyyzz[k];

                g_xx_0_zz_yyzz[k] = -g_xx_0_z_yyzz[k] * ab_z + g_xx_0_z_yyzzz[k];

                g_xx_0_zz_yzzz[k] = -g_xx_0_z_yzzz[k] * ab_z + g_xx_0_z_yzzzz[k];

                g_xx_0_zz_zzzz[k] = -g_xx_0_z_zzzz[k] * ab_z + g_xx_0_z_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxx, g_xy_0_x_xxxxx, g_xy_0_x_xxxxy, g_xy_0_x_xxxxz, g_xy_0_x_xxxy, g_xy_0_x_xxxyy, g_xy_0_x_xxxyz, g_xy_0_x_xxxz, g_xy_0_x_xxxzz, g_xy_0_x_xxyy, g_xy_0_x_xxyyy, g_xy_0_x_xxyyz, g_xy_0_x_xxyz, g_xy_0_x_xxyzz, g_xy_0_x_xxzz, g_xy_0_x_xxzzz, g_xy_0_x_xyyy, g_xy_0_x_xyyyy, g_xy_0_x_xyyyz, g_xy_0_x_xyyz, g_xy_0_x_xyyzz, g_xy_0_x_xyzz, g_xy_0_x_xyzzz, g_xy_0_x_xzzz, g_xy_0_x_xzzzz, g_xy_0_x_yyyy, g_xy_0_x_yyyz, g_xy_0_x_yyzz, g_xy_0_x_yzzz, g_xy_0_x_zzzz, g_xy_0_xx_xxxx, g_xy_0_xx_xxxy, g_xy_0_xx_xxxz, g_xy_0_xx_xxyy, g_xy_0_xx_xxyz, g_xy_0_xx_xxzz, g_xy_0_xx_xyyy, g_xy_0_xx_xyyz, g_xy_0_xx_xyzz, g_xy_0_xx_xzzz, g_xy_0_xx_yyyy, g_xy_0_xx_yyyz, g_xy_0_xx_yyzz, g_xy_0_xx_yzzz, g_xy_0_xx_zzzz, g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxzz, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyzz, g_y_0_x_xzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xx_xxxx[k] = -g_y_0_x_xxxx[k] - g_xy_0_x_xxxx[k] * ab_x + g_xy_0_x_xxxxx[k];

                g_xy_0_xx_xxxy[k] = -g_y_0_x_xxxy[k] - g_xy_0_x_xxxy[k] * ab_x + g_xy_0_x_xxxxy[k];

                g_xy_0_xx_xxxz[k] = -g_y_0_x_xxxz[k] - g_xy_0_x_xxxz[k] * ab_x + g_xy_0_x_xxxxz[k];

                g_xy_0_xx_xxyy[k] = -g_y_0_x_xxyy[k] - g_xy_0_x_xxyy[k] * ab_x + g_xy_0_x_xxxyy[k];

                g_xy_0_xx_xxyz[k] = -g_y_0_x_xxyz[k] - g_xy_0_x_xxyz[k] * ab_x + g_xy_0_x_xxxyz[k];

                g_xy_0_xx_xxzz[k] = -g_y_0_x_xxzz[k] - g_xy_0_x_xxzz[k] * ab_x + g_xy_0_x_xxxzz[k];

                g_xy_0_xx_xyyy[k] = -g_y_0_x_xyyy[k] - g_xy_0_x_xyyy[k] * ab_x + g_xy_0_x_xxyyy[k];

                g_xy_0_xx_xyyz[k] = -g_y_0_x_xyyz[k] - g_xy_0_x_xyyz[k] * ab_x + g_xy_0_x_xxyyz[k];

                g_xy_0_xx_xyzz[k] = -g_y_0_x_xyzz[k] - g_xy_0_x_xyzz[k] * ab_x + g_xy_0_x_xxyzz[k];

                g_xy_0_xx_xzzz[k] = -g_y_0_x_xzzz[k] - g_xy_0_x_xzzz[k] * ab_x + g_xy_0_x_xxzzz[k];

                g_xy_0_xx_yyyy[k] = -g_y_0_x_yyyy[k] - g_xy_0_x_yyyy[k] * ab_x + g_xy_0_x_xyyyy[k];

                g_xy_0_xx_yyyz[k] = -g_y_0_x_yyyz[k] - g_xy_0_x_yyyz[k] * ab_x + g_xy_0_x_xyyyz[k];

                g_xy_0_xx_yyzz[k] = -g_y_0_x_yyzz[k] - g_xy_0_x_yyzz[k] * ab_x + g_xy_0_x_xyyzz[k];

                g_xy_0_xx_yzzz[k] = -g_y_0_x_yzzz[k] - g_xy_0_x_yzzz[k] * ab_x + g_xy_0_x_xyzzz[k];

                g_xy_0_xx_zzzz[k] = -g_y_0_x_zzzz[k] - g_xy_0_x_zzzz[k] * ab_x + g_xy_0_x_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_xxxx, g_xy_0_xy_xxxy, g_xy_0_xy_xxxz, g_xy_0_xy_xxyy, g_xy_0_xy_xxyz, g_xy_0_xy_xxzz, g_xy_0_xy_xyyy, g_xy_0_xy_xyyz, g_xy_0_xy_xyzz, g_xy_0_xy_xzzz, g_xy_0_xy_yyyy, g_xy_0_xy_yyyz, g_xy_0_xy_yyzz, g_xy_0_xy_yzzz, g_xy_0_xy_zzzz, g_xy_0_y_xxxx, g_xy_0_y_xxxxx, g_xy_0_y_xxxxy, g_xy_0_y_xxxxz, g_xy_0_y_xxxy, g_xy_0_y_xxxyy, g_xy_0_y_xxxyz, g_xy_0_y_xxxz, g_xy_0_y_xxxzz, g_xy_0_y_xxyy, g_xy_0_y_xxyyy, g_xy_0_y_xxyyz, g_xy_0_y_xxyz, g_xy_0_y_xxyzz, g_xy_0_y_xxzz, g_xy_0_y_xxzzz, g_xy_0_y_xyyy, g_xy_0_y_xyyyy, g_xy_0_y_xyyyz, g_xy_0_y_xyyz, g_xy_0_y_xyyzz, g_xy_0_y_xyzz, g_xy_0_y_xyzzz, g_xy_0_y_xzzz, g_xy_0_y_xzzzz, g_xy_0_y_yyyy, g_xy_0_y_yyyz, g_xy_0_y_yyzz, g_xy_0_y_yzzz, g_xy_0_y_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xy_xxxx[k] = -g_y_0_y_xxxx[k] - g_xy_0_y_xxxx[k] * ab_x + g_xy_0_y_xxxxx[k];

                g_xy_0_xy_xxxy[k] = -g_y_0_y_xxxy[k] - g_xy_0_y_xxxy[k] * ab_x + g_xy_0_y_xxxxy[k];

                g_xy_0_xy_xxxz[k] = -g_y_0_y_xxxz[k] - g_xy_0_y_xxxz[k] * ab_x + g_xy_0_y_xxxxz[k];

                g_xy_0_xy_xxyy[k] = -g_y_0_y_xxyy[k] - g_xy_0_y_xxyy[k] * ab_x + g_xy_0_y_xxxyy[k];

                g_xy_0_xy_xxyz[k] = -g_y_0_y_xxyz[k] - g_xy_0_y_xxyz[k] * ab_x + g_xy_0_y_xxxyz[k];

                g_xy_0_xy_xxzz[k] = -g_y_0_y_xxzz[k] - g_xy_0_y_xxzz[k] * ab_x + g_xy_0_y_xxxzz[k];

                g_xy_0_xy_xyyy[k] = -g_y_0_y_xyyy[k] - g_xy_0_y_xyyy[k] * ab_x + g_xy_0_y_xxyyy[k];

                g_xy_0_xy_xyyz[k] = -g_y_0_y_xyyz[k] - g_xy_0_y_xyyz[k] * ab_x + g_xy_0_y_xxyyz[k];

                g_xy_0_xy_xyzz[k] = -g_y_0_y_xyzz[k] - g_xy_0_y_xyzz[k] * ab_x + g_xy_0_y_xxyzz[k];

                g_xy_0_xy_xzzz[k] = -g_y_0_y_xzzz[k] - g_xy_0_y_xzzz[k] * ab_x + g_xy_0_y_xxzzz[k];

                g_xy_0_xy_yyyy[k] = -g_y_0_y_yyyy[k] - g_xy_0_y_yyyy[k] * ab_x + g_xy_0_y_xyyyy[k];

                g_xy_0_xy_yyyz[k] = -g_y_0_y_yyyz[k] - g_xy_0_y_yyyz[k] * ab_x + g_xy_0_y_xyyyz[k];

                g_xy_0_xy_yyzz[k] = -g_y_0_y_yyzz[k] - g_xy_0_y_yyzz[k] * ab_x + g_xy_0_y_xyyzz[k];

                g_xy_0_xy_yzzz[k] = -g_y_0_y_yzzz[k] - g_xy_0_y_yzzz[k] * ab_x + g_xy_0_y_xyzzz[k];

                g_xy_0_xy_zzzz[k] = -g_y_0_y_zzzz[k] - g_xy_0_y_zzzz[k] * ab_x + g_xy_0_y_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxx, g_xy_0_x_xxxxz, g_xy_0_x_xxxy, g_xy_0_x_xxxyz, g_xy_0_x_xxxz, g_xy_0_x_xxxzz, g_xy_0_x_xxyy, g_xy_0_x_xxyyz, g_xy_0_x_xxyz, g_xy_0_x_xxyzz, g_xy_0_x_xxzz, g_xy_0_x_xxzzz, g_xy_0_x_xyyy, g_xy_0_x_xyyyz, g_xy_0_x_xyyz, g_xy_0_x_xyyzz, g_xy_0_x_xyzz, g_xy_0_x_xyzzz, g_xy_0_x_xzzz, g_xy_0_x_xzzzz, g_xy_0_x_yyyy, g_xy_0_x_yyyyz, g_xy_0_x_yyyz, g_xy_0_x_yyyzz, g_xy_0_x_yyzz, g_xy_0_x_yyzzz, g_xy_0_x_yzzz, g_xy_0_x_yzzzz, g_xy_0_x_zzzz, g_xy_0_x_zzzzz, g_xy_0_xz_xxxx, g_xy_0_xz_xxxy, g_xy_0_xz_xxxz, g_xy_0_xz_xxyy, g_xy_0_xz_xxyz, g_xy_0_xz_xxzz, g_xy_0_xz_xyyy, g_xy_0_xz_xyyz, g_xy_0_xz_xyzz, g_xy_0_xz_xzzz, g_xy_0_xz_yyyy, g_xy_0_xz_yyyz, g_xy_0_xz_yyzz, g_xy_0_xz_yzzz, g_xy_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xz_xxxx[k] = -g_xy_0_x_xxxx[k] * ab_z + g_xy_0_x_xxxxz[k];

                g_xy_0_xz_xxxy[k] = -g_xy_0_x_xxxy[k] * ab_z + g_xy_0_x_xxxyz[k];

                g_xy_0_xz_xxxz[k] = -g_xy_0_x_xxxz[k] * ab_z + g_xy_0_x_xxxzz[k];

                g_xy_0_xz_xxyy[k] = -g_xy_0_x_xxyy[k] * ab_z + g_xy_0_x_xxyyz[k];

                g_xy_0_xz_xxyz[k] = -g_xy_0_x_xxyz[k] * ab_z + g_xy_0_x_xxyzz[k];

                g_xy_0_xz_xxzz[k] = -g_xy_0_x_xxzz[k] * ab_z + g_xy_0_x_xxzzz[k];

                g_xy_0_xz_xyyy[k] = -g_xy_0_x_xyyy[k] * ab_z + g_xy_0_x_xyyyz[k];

                g_xy_0_xz_xyyz[k] = -g_xy_0_x_xyyz[k] * ab_z + g_xy_0_x_xyyzz[k];

                g_xy_0_xz_xyzz[k] = -g_xy_0_x_xyzz[k] * ab_z + g_xy_0_x_xyzzz[k];

                g_xy_0_xz_xzzz[k] = -g_xy_0_x_xzzz[k] * ab_z + g_xy_0_x_xzzzz[k];

                g_xy_0_xz_yyyy[k] = -g_xy_0_x_yyyy[k] * ab_z + g_xy_0_x_yyyyz[k];

                g_xy_0_xz_yyyz[k] = -g_xy_0_x_yyyz[k] * ab_z + g_xy_0_x_yyyzz[k];

                g_xy_0_xz_yyzz[k] = -g_xy_0_x_yyzz[k] * ab_z + g_xy_0_x_yyzzz[k];

                g_xy_0_xz_yzzz[k] = -g_xy_0_x_yzzz[k] * ab_z + g_xy_0_x_yzzzz[k];

                g_xy_0_xz_zzzz[k] = -g_xy_0_x_zzzz[k] * ab_z + g_xy_0_x_zzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxx, g_x_0_y_xxxy, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyzz, g_x_0_y_yzzz, g_x_0_y_zzzz, g_xy_0_y_xxxx, g_xy_0_y_xxxxy, g_xy_0_y_xxxy, g_xy_0_y_xxxyy, g_xy_0_y_xxxyz, g_xy_0_y_xxxz, g_xy_0_y_xxyy, g_xy_0_y_xxyyy, g_xy_0_y_xxyyz, g_xy_0_y_xxyz, g_xy_0_y_xxyzz, g_xy_0_y_xxzz, g_xy_0_y_xyyy, g_xy_0_y_xyyyy, g_xy_0_y_xyyyz, g_xy_0_y_xyyz, g_xy_0_y_xyyzz, g_xy_0_y_xyzz, g_xy_0_y_xyzzz, g_xy_0_y_xzzz, g_xy_0_y_yyyy, g_xy_0_y_yyyyy, g_xy_0_y_yyyyz, g_xy_0_y_yyyz, g_xy_0_y_yyyzz, g_xy_0_y_yyzz, g_xy_0_y_yyzzz, g_xy_0_y_yzzz, g_xy_0_y_yzzzz, g_xy_0_y_zzzz, g_xy_0_yy_xxxx, g_xy_0_yy_xxxy, g_xy_0_yy_xxxz, g_xy_0_yy_xxyy, g_xy_0_yy_xxyz, g_xy_0_yy_xxzz, g_xy_0_yy_xyyy, g_xy_0_yy_xyyz, g_xy_0_yy_xyzz, g_xy_0_yy_xzzz, g_xy_0_yy_yyyy, g_xy_0_yy_yyyz, g_xy_0_yy_yyzz, g_xy_0_yy_yzzz, g_xy_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yy_xxxx[k] = -g_x_0_y_xxxx[k] - g_xy_0_y_xxxx[k] * ab_y + g_xy_0_y_xxxxy[k];

                g_xy_0_yy_xxxy[k] = -g_x_0_y_xxxy[k] - g_xy_0_y_xxxy[k] * ab_y + g_xy_0_y_xxxyy[k];

                g_xy_0_yy_xxxz[k] = -g_x_0_y_xxxz[k] - g_xy_0_y_xxxz[k] * ab_y + g_xy_0_y_xxxyz[k];

                g_xy_0_yy_xxyy[k] = -g_x_0_y_xxyy[k] - g_xy_0_y_xxyy[k] * ab_y + g_xy_0_y_xxyyy[k];

                g_xy_0_yy_xxyz[k] = -g_x_0_y_xxyz[k] - g_xy_0_y_xxyz[k] * ab_y + g_xy_0_y_xxyyz[k];

                g_xy_0_yy_xxzz[k] = -g_x_0_y_xxzz[k] - g_xy_0_y_xxzz[k] * ab_y + g_xy_0_y_xxyzz[k];

                g_xy_0_yy_xyyy[k] = -g_x_0_y_xyyy[k] - g_xy_0_y_xyyy[k] * ab_y + g_xy_0_y_xyyyy[k];

                g_xy_0_yy_xyyz[k] = -g_x_0_y_xyyz[k] - g_xy_0_y_xyyz[k] * ab_y + g_xy_0_y_xyyyz[k];

                g_xy_0_yy_xyzz[k] = -g_x_0_y_xyzz[k] - g_xy_0_y_xyzz[k] * ab_y + g_xy_0_y_xyyzz[k];

                g_xy_0_yy_xzzz[k] = -g_x_0_y_xzzz[k] - g_xy_0_y_xzzz[k] * ab_y + g_xy_0_y_xyzzz[k];

                g_xy_0_yy_yyyy[k] = -g_x_0_y_yyyy[k] - g_xy_0_y_yyyy[k] * ab_y + g_xy_0_y_yyyyy[k];

                g_xy_0_yy_yyyz[k] = -g_x_0_y_yyyz[k] - g_xy_0_y_yyyz[k] * ab_y + g_xy_0_y_yyyyz[k];

                g_xy_0_yy_yyzz[k] = -g_x_0_y_yyzz[k] - g_xy_0_y_yyzz[k] * ab_y + g_xy_0_y_yyyzz[k];

                g_xy_0_yy_yzzz[k] = -g_x_0_y_yzzz[k] - g_xy_0_y_yzzz[k] * ab_y + g_xy_0_y_yyzzz[k];

                g_xy_0_yy_zzzz[k] = -g_x_0_y_zzzz[k] - g_xy_0_y_zzzz[k] * ab_y + g_xy_0_y_yzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_y_xxxx, g_xy_0_y_xxxxz, g_xy_0_y_xxxy, g_xy_0_y_xxxyz, g_xy_0_y_xxxz, g_xy_0_y_xxxzz, g_xy_0_y_xxyy, g_xy_0_y_xxyyz, g_xy_0_y_xxyz, g_xy_0_y_xxyzz, g_xy_0_y_xxzz, g_xy_0_y_xxzzz, g_xy_0_y_xyyy, g_xy_0_y_xyyyz, g_xy_0_y_xyyz, g_xy_0_y_xyyzz, g_xy_0_y_xyzz, g_xy_0_y_xyzzz, g_xy_0_y_xzzz, g_xy_0_y_xzzzz, g_xy_0_y_yyyy, g_xy_0_y_yyyyz, g_xy_0_y_yyyz, g_xy_0_y_yyyzz, g_xy_0_y_yyzz, g_xy_0_y_yyzzz, g_xy_0_y_yzzz, g_xy_0_y_yzzzz, g_xy_0_y_zzzz, g_xy_0_y_zzzzz, g_xy_0_yz_xxxx, g_xy_0_yz_xxxy, g_xy_0_yz_xxxz, g_xy_0_yz_xxyy, g_xy_0_yz_xxyz, g_xy_0_yz_xxzz, g_xy_0_yz_xyyy, g_xy_0_yz_xyyz, g_xy_0_yz_xyzz, g_xy_0_yz_xzzz, g_xy_0_yz_yyyy, g_xy_0_yz_yyyz, g_xy_0_yz_yyzz, g_xy_0_yz_yzzz, g_xy_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yz_xxxx[k] = -g_xy_0_y_xxxx[k] * ab_z + g_xy_0_y_xxxxz[k];

                g_xy_0_yz_xxxy[k] = -g_xy_0_y_xxxy[k] * ab_z + g_xy_0_y_xxxyz[k];

                g_xy_0_yz_xxxz[k] = -g_xy_0_y_xxxz[k] * ab_z + g_xy_0_y_xxxzz[k];

                g_xy_0_yz_xxyy[k] = -g_xy_0_y_xxyy[k] * ab_z + g_xy_0_y_xxyyz[k];

                g_xy_0_yz_xxyz[k] = -g_xy_0_y_xxyz[k] * ab_z + g_xy_0_y_xxyzz[k];

                g_xy_0_yz_xxzz[k] = -g_xy_0_y_xxzz[k] * ab_z + g_xy_0_y_xxzzz[k];

                g_xy_0_yz_xyyy[k] = -g_xy_0_y_xyyy[k] * ab_z + g_xy_0_y_xyyyz[k];

                g_xy_0_yz_xyyz[k] = -g_xy_0_y_xyyz[k] * ab_z + g_xy_0_y_xyyzz[k];

                g_xy_0_yz_xyzz[k] = -g_xy_0_y_xyzz[k] * ab_z + g_xy_0_y_xyzzz[k];

                g_xy_0_yz_xzzz[k] = -g_xy_0_y_xzzz[k] * ab_z + g_xy_0_y_xzzzz[k];

                g_xy_0_yz_yyyy[k] = -g_xy_0_y_yyyy[k] * ab_z + g_xy_0_y_yyyyz[k];

                g_xy_0_yz_yyyz[k] = -g_xy_0_y_yyyz[k] * ab_z + g_xy_0_y_yyyzz[k];

                g_xy_0_yz_yyzz[k] = -g_xy_0_y_yyzz[k] * ab_z + g_xy_0_y_yyzzz[k];

                g_xy_0_yz_yzzz[k] = -g_xy_0_y_yzzz[k] * ab_z + g_xy_0_y_yzzzz[k];

                g_xy_0_yz_zzzz[k] = -g_xy_0_y_zzzz[k] * ab_z + g_xy_0_y_zzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_z_xxxx, g_xy_0_z_xxxxz, g_xy_0_z_xxxy, g_xy_0_z_xxxyz, g_xy_0_z_xxxz, g_xy_0_z_xxxzz, g_xy_0_z_xxyy, g_xy_0_z_xxyyz, g_xy_0_z_xxyz, g_xy_0_z_xxyzz, g_xy_0_z_xxzz, g_xy_0_z_xxzzz, g_xy_0_z_xyyy, g_xy_0_z_xyyyz, g_xy_0_z_xyyz, g_xy_0_z_xyyzz, g_xy_0_z_xyzz, g_xy_0_z_xyzzz, g_xy_0_z_xzzz, g_xy_0_z_xzzzz, g_xy_0_z_yyyy, g_xy_0_z_yyyyz, g_xy_0_z_yyyz, g_xy_0_z_yyyzz, g_xy_0_z_yyzz, g_xy_0_z_yyzzz, g_xy_0_z_yzzz, g_xy_0_z_yzzzz, g_xy_0_z_zzzz, g_xy_0_z_zzzzz, g_xy_0_zz_xxxx, g_xy_0_zz_xxxy, g_xy_0_zz_xxxz, g_xy_0_zz_xxyy, g_xy_0_zz_xxyz, g_xy_0_zz_xxzz, g_xy_0_zz_xyyy, g_xy_0_zz_xyyz, g_xy_0_zz_xyzz, g_xy_0_zz_xzzz, g_xy_0_zz_yyyy, g_xy_0_zz_yyyz, g_xy_0_zz_yyzz, g_xy_0_zz_yzzz, g_xy_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zz_xxxx[k] = -g_xy_0_z_xxxx[k] * ab_z + g_xy_0_z_xxxxz[k];

                g_xy_0_zz_xxxy[k] = -g_xy_0_z_xxxy[k] * ab_z + g_xy_0_z_xxxyz[k];

                g_xy_0_zz_xxxz[k] = -g_xy_0_z_xxxz[k] * ab_z + g_xy_0_z_xxxzz[k];

                g_xy_0_zz_xxyy[k] = -g_xy_0_z_xxyy[k] * ab_z + g_xy_0_z_xxyyz[k];

                g_xy_0_zz_xxyz[k] = -g_xy_0_z_xxyz[k] * ab_z + g_xy_0_z_xxyzz[k];

                g_xy_0_zz_xxzz[k] = -g_xy_0_z_xxzz[k] * ab_z + g_xy_0_z_xxzzz[k];

                g_xy_0_zz_xyyy[k] = -g_xy_0_z_xyyy[k] * ab_z + g_xy_0_z_xyyyz[k];

                g_xy_0_zz_xyyz[k] = -g_xy_0_z_xyyz[k] * ab_z + g_xy_0_z_xyyzz[k];

                g_xy_0_zz_xyzz[k] = -g_xy_0_z_xyzz[k] * ab_z + g_xy_0_z_xyzzz[k];

                g_xy_0_zz_xzzz[k] = -g_xy_0_z_xzzz[k] * ab_z + g_xy_0_z_xzzzz[k];

                g_xy_0_zz_yyyy[k] = -g_xy_0_z_yyyy[k] * ab_z + g_xy_0_z_yyyyz[k];

                g_xy_0_zz_yyyz[k] = -g_xy_0_z_yyyz[k] * ab_z + g_xy_0_z_yyyzz[k];

                g_xy_0_zz_yyzz[k] = -g_xy_0_z_yyzz[k] * ab_z + g_xy_0_z_yyzzz[k];

                g_xy_0_zz_yzzz[k] = -g_xy_0_z_yzzz[k] * ab_z + g_xy_0_z_yzzzz[k];

                g_xy_0_zz_zzzz[k] = -g_xy_0_z_zzzz[k] * ab_z + g_xy_0_z_zzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxx, g_xz_0_x_xxxxx, g_xz_0_x_xxxxy, g_xz_0_x_xxxxz, g_xz_0_x_xxxy, g_xz_0_x_xxxyy, g_xz_0_x_xxxyz, g_xz_0_x_xxxz, g_xz_0_x_xxxzz, g_xz_0_x_xxyy, g_xz_0_x_xxyyy, g_xz_0_x_xxyyz, g_xz_0_x_xxyz, g_xz_0_x_xxyzz, g_xz_0_x_xxzz, g_xz_0_x_xxzzz, g_xz_0_x_xyyy, g_xz_0_x_xyyyy, g_xz_0_x_xyyyz, g_xz_0_x_xyyz, g_xz_0_x_xyyzz, g_xz_0_x_xyzz, g_xz_0_x_xyzzz, g_xz_0_x_xzzz, g_xz_0_x_xzzzz, g_xz_0_x_yyyy, g_xz_0_x_yyyz, g_xz_0_x_yyzz, g_xz_0_x_yzzz, g_xz_0_x_zzzz, g_xz_0_xx_xxxx, g_xz_0_xx_xxxy, g_xz_0_xx_xxxz, g_xz_0_xx_xxyy, g_xz_0_xx_xxyz, g_xz_0_xx_xxzz, g_xz_0_xx_xyyy, g_xz_0_xx_xyyz, g_xz_0_xx_xyzz, g_xz_0_xx_xzzz, g_xz_0_xx_yyyy, g_xz_0_xx_yyyz, g_xz_0_xx_yyzz, g_xz_0_xx_yzzz, g_xz_0_xx_zzzz, g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxzz, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyzz, g_z_0_x_xzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xx_xxxx[k] = -g_z_0_x_xxxx[k] - g_xz_0_x_xxxx[k] * ab_x + g_xz_0_x_xxxxx[k];

                g_xz_0_xx_xxxy[k] = -g_z_0_x_xxxy[k] - g_xz_0_x_xxxy[k] * ab_x + g_xz_0_x_xxxxy[k];

                g_xz_0_xx_xxxz[k] = -g_z_0_x_xxxz[k] - g_xz_0_x_xxxz[k] * ab_x + g_xz_0_x_xxxxz[k];

                g_xz_0_xx_xxyy[k] = -g_z_0_x_xxyy[k] - g_xz_0_x_xxyy[k] * ab_x + g_xz_0_x_xxxyy[k];

                g_xz_0_xx_xxyz[k] = -g_z_0_x_xxyz[k] - g_xz_0_x_xxyz[k] * ab_x + g_xz_0_x_xxxyz[k];

                g_xz_0_xx_xxzz[k] = -g_z_0_x_xxzz[k] - g_xz_0_x_xxzz[k] * ab_x + g_xz_0_x_xxxzz[k];

                g_xz_0_xx_xyyy[k] = -g_z_0_x_xyyy[k] - g_xz_0_x_xyyy[k] * ab_x + g_xz_0_x_xxyyy[k];

                g_xz_0_xx_xyyz[k] = -g_z_0_x_xyyz[k] - g_xz_0_x_xyyz[k] * ab_x + g_xz_0_x_xxyyz[k];

                g_xz_0_xx_xyzz[k] = -g_z_0_x_xyzz[k] - g_xz_0_x_xyzz[k] * ab_x + g_xz_0_x_xxyzz[k];

                g_xz_0_xx_xzzz[k] = -g_z_0_x_xzzz[k] - g_xz_0_x_xzzz[k] * ab_x + g_xz_0_x_xxzzz[k];

                g_xz_0_xx_yyyy[k] = -g_z_0_x_yyyy[k] - g_xz_0_x_yyyy[k] * ab_x + g_xz_0_x_xyyyy[k];

                g_xz_0_xx_yyyz[k] = -g_z_0_x_yyyz[k] - g_xz_0_x_yyyz[k] * ab_x + g_xz_0_x_xyyyz[k];

                g_xz_0_xx_yyzz[k] = -g_z_0_x_yyzz[k] - g_xz_0_x_yyzz[k] * ab_x + g_xz_0_x_xyyzz[k];

                g_xz_0_xx_yzzz[k] = -g_z_0_x_yzzz[k] - g_xz_0_x_yzzz[k] * ab_x + g_xz_0_x_xyzzz[k];

                g_xz_0_xx_zzzz[k] = -g_z_0_x_zzzz[k] - g_xz_0_x_zzzz[k] * ab_x + g_xz_0_x_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxx, g_xz_0_x_xxxxy, g_xz_0_x_xxxy, g_xz_0_x_xxxyy, g_xz_0_x_xxxyz, g_xz_0_x_xxxz, g_xz_0_x_xxyy, g_xz_0_x_xxyyy, g_xz_0_x_xxyyz, g_xz_0_x_xxyz, g_xz_0_x_xxyzz, g_xz_0_x_xxzz, g_xz_0_x_xyyy, g_xz_0_x_xyyyy, g_xz_0_x_xyyyz, g_xz_0_x_xyyz, g_xz_0_x_xyyzz, g_xz_0_x_xyzz, g_xz_0_x_xyzzz, g_xz_0_x_xzzz, g_xz_0_x_yyyy, g_xz_0_x_yyyyy, g_xz_0_x_yyyyz, g_xz_0_x_yyyz, g_xz_0_x_yyyzz, g_xz_0_x_yyzz, g_xz_0_x_yyzzz, g_xz_0_x_yzzz, g_xz_0_x_yzzzz, g_xz_0_x_zzzz, g_xz_0_xy_xxxx, g_xz_0_xy_xxxy, g_xz_0_xy_xxxz, g_xz_0_xy_xxyy, g_xz_0_xy_xxyz, g_xz_0_xy_xxzz, g_xz_0_xy_xyyy, g_xz_0_xy_xyyz, g_xz_0_xy_xyzz, g_xz_0_xy_xzzz, g_xz_0_xy_yyyy, g_xz_0_xy_yyyz, g_xz_0_xy_yyzz, g_xz_0_xy_yzzz, g_xz_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xy_xxxx[k] = -g_xz_0_x_xxxx[k] * ab_y + g_xz_0_x_xxxxy[k];

                g_xz_0_xy_xxxy[k] = -g_xz_0_x_xxxy[k] * ab_y + g_xz_0_x_xxxyy[k];

                g_xz_0_xy_xxxz[k] = -g_xz_0_x_xxxz[k] * ab_y + g_xz_0_x_xxxyz[k];

                g_xz_0_xy_xxyy[k] = -g_xz_0_x_xxyy[k] * ab_y + g_xz_0_x_xxyyy[k];

                g_xz_0_xy_xxyz[k] = -g_xz_0_x_xxyz[k] * ab_y + g_xz_0_x_xxyyz[k];

                g_xz_0_xy_xxzz[k] = -g_xz_0_x_xxzz[k] * ab_y + g_xz_0_x_xxyzz[k];

                g_xz_0_xy_xyyy[k] = -g_xz_0_x_xyyy[k] * ab_y + g_xz_0_x_xyyyy[k];

                g_xz_0_xy_xyyz[k] = -g_xz_0_x_xyyz[k] * ab_y + g_xz_0_x_xyyyz[k];

                g_xz_0_xy_xyzz[k] = -g_xz_0_x_xyzz[k] * ab_y + g_xz_0_x_xyyzz[k];

                g_xz_0_xy_xzzz[k] = -g_xz_0_x_xzzz[k] * ab_y + g_xz_0_x_xyzzz[k];

                g_xz_0_xy_yyyy[k] = -g_xz_0_x_yyyy[k] * ab_y + g_xz_0_x_yyyyy[k];

                g_xz_0_xy_yyyz[k] = -g_xz_0_x_yyyz[k] * ab_y + g_xz_0_x_yyyyz[k];

                g_xz_0_xy_yyzz[k] = -g_xz_0_x_yyzz[k] * ab_y + g_xz_0_x_yyyzz[k];

                g_xz_0_xy_yzzz[k] = -g_xz_0_x_yzzz[k] * ab_y + g_xz_0_x_yyzzz[k];

                g_xz_0_xy_zzzz[k] = -g_xz_0_x_zzzz[k] * ab_y + g_xz_0_x_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xz_xxxx, g_xz_0_xz_xxxy, g_xz_0_xz_xxxz, g_xz_0_xz_xxyy, g_xz_0_xz_xxyz, g_xz_0_xz_xxzz, g_xz_0_xz_xyyy, g_xz_0_xz_xyyz, g_xz_0_xz_xyzz, g_xz_0_xz_xzzz, g_xz_0_xz_yyyy, g_xz_0_xz_yyyz, g_xz_0_xz_yyzz, g_xz_0_xz_yzzz, g_xz_0_xz_zzzz, g_xz_0_z_xxxx, g_xz_0_z_xxxxx, g_xz_0_z_xxxxy, g_xz_0_z_xxxxz, g_xz_0_z_xxxy, g_xz_0_z_xxxyy, g_xz_0_z_xxxyz, g_xz_0_z_xxxz, g_xz_0_z_xxxzz, g_xz_0_z_xxyy, g_xz_0_z_xxyyy, g_xz_0_z_xxyyz, g_xz_0_z_xxyz, g_xz_0_z_xxyzz, g_xz_0_z_xxzz, g_xz_0_z_xxzzz, g_xz_0_z_xyyy, g_xz_0_z_xyyyy, g_xz_0_z_xyyyz, g_xz_0_z_xyyz, g_xz_0_z_xyyzz, g_xz_0_z_xyzz, g_xz_0_z_xyzzz, g_xz_0_z_xzzz, g_xz_0_z_xzzzz, g_xz_0_z_yyyy, g_xz_0_z_yyyz, g_xz_0_z_yyzz, g_xz_0_z_yzzz, g_xz_0_z_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xz_xxxx[k] = -g_z_0_z_xxxx[k] - g_xz_0_z_xxxx[k] * ab_x + g_xz_0_z_xxxxx[k];

                g_xz_0_xz_xxxy[k] = -g_z_0_z_xxxy[k] - g_xz_0_z_xxxy[k] * ab_x + g_xz_0_z_xxxxy[k];

                g_xz_0_xz_xxxz[k] = -g_z_0_z_xxxz[k] - g_xz_0_z_xxxz[k] * ab_x + g_xz_0_z_xxxxz[k];

                g_xz_0_xz_xxyy[k] = -g_z_0_z_xxyy[k] - g_xz_0_z_xxyy[k] * ab_x + g_xz_0_z_xxxyy[k];

                g_xz_0_xz_xxyz[k] = -g_z_0_z_xxyz[k] - g_xz_0_z_xxyz[k] * ab_x + g_xz_0_z_xxxyz[k];

                g_xz_0_xz_xxzz[k] = -g_z_0_z_xxzz[k] - g_xz_0_z_xxzz[k] * ab_x + g_xz_0_z_xxxzz[k];

                g_xz_0_xz_xyyy[k] = -g_z_0_z_xyyy[k] - g_xz_0_z_xyyy[k] * ab_x + g_xz_0_z_xxyyy[k];

                g_xz_0_xz_xyyz[k] = -g_z_0_z_xyyz[k] - g_xz_0_z_xyyz[k] * ab_x + g_xz_0_z_xxyyz[k];

                g_xz_0_xz_xyzz[k] = -g_z_0_z_xyzz[k] - g_xz_0_z_xyzz[k] * ab_x + g_xz_0_z_xxyzz[k];

                g_xz_0_xz_xzzz[k] = -g_z_0_z_xzzz[k] - g_xz_0_z_xzzz[k] * ab_x + g_xz_0_z_xxzzz[k];

                g_xz_0_xz_yyyy[k] = -g_z_0_z_yyyy[k] - g_xz_0_z_yyyy[k] * ab_x + g_xz_0_z_xyyyy[k];

                g_xz_0_xz_yyyz[k] = -g_z_0_z_yyyz[k] - g_xz_0_z_yyyz[k] * ab_x + g_xz_0_z_xyyyz[k];

                g_xz_0_xz_yyzz[k] = -g_z_0_z_yyzz[k] - g_xz_0_z_yyzz[k] * ab_x + g_xz_0_z_xyyzz[k];

                g_xz_0_xz_yzzz[k] = -g_z_0_z_yzzz[k] - g_xz_0_z_yzzz[k] * ab_x + g_xz_0_z_xyzzz[k];

                g_xz_0_xz_zzzz[k] = -g_z_0_z_zzzz[k] - g_xz_0_z_zzzz[k] * ab_x + g_xz_0_z_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_y_xxxx, g_xz_0_y_xxxxy, g_xz_0_y_xxxy, g_xz_0_y_xxxyy, g_xz_0_y_xxxyz, g_xz_0_y_xxxz, g_xz_0_y_xxyy, g_xz_0_y_xxyyy, g_xz_0_y_xxyyz, g_xz_0_y_xxyz, g_xz_0_y_xxyzz, g_xz_0_y_xxzz, g_xz_0_y_xyyy, g_xz_0_y_xyyyy, g_xz_0_y_xyyyz, g_xz_0_y_xyyz, g_xz_0_y_xyyzz, g_xz_0_y_xyzz, g_xz_0_y_xyzzz, g_xz_0_y_xzzz, g_xz_0_y_yyyy, g_xz_0_y_yyyyy, g_xz_0_y_yyyyz, g_xz_0_y_yyyz, g_xz_0_y_yyyzz, g_xz_0_y_yyzz, g_xz_0_y_yyzzz, g_xz_0_y_yzzz, g_xz_0_y_yzzzz, g_xz_0_y_zzzz, g_xz_0_yy_xxxx, g_xz_0_yy_xxxy, g_xz_0_yy_xxxz, g_xz_0_yy_xxyy, g_xz_0_yy_xxyz, g_xz_0_yy_xxzz, g_xz_0_yy_xyyy, g_xz_0_yy_xyyz, g_xz_0_yy_xyzz, g_xz_0_yy_xzzz, g_xz_0_yy_yyyy, g_xz_0_yy_yyyz, g_xz_0_yy_yyzz, g_xz_0_yy_yzzz, g_xz_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yy_xxxx[k] = -g_xz_0_y_xxxx[k] * ab_y + g_xz_0_y_xxxxy[k];

                g_xz_0_yy_xxxy[k] = -g_xz_0_y_xxxy[k] * ab_y + g_xz_0_y_xxxyy[k];

                g_xz_0_yy_xxxz[k] = -g_xz_0_y_xxxz[k] * ab_y + g_xz_0_y_xxxyz[k];

                g_xz_0_yy_xxyy[k] = -g_xz_0_y_xxyy[k] * ab_y + g_xz_0_y_xxyyy[k];

                g_xz_0_yy_xxyz[k] = -g_xz_0_y_xxyz[k] * ab_y + g_xz_0_y_xxyyz[k];

                g_xz_0_yy_xxzz[k] = -g_xz_0_y_xxzz[k] * ab_y + g_xz_0_y_xxyzz[k];

                g_xz_0_yy_xyyy[k] = -g_xz_0_y_xyyy[k] * ab_y + g_xz_0_y_xyyyy[k];

                g_xz_0_yy_xyyz[k] = -g_xz_0_y_xyyz[k] * ab_y + g_xz_0_y_xyyyz[k];

                g_xz_0_yy_xyzz[k] = -g_xz_0_y_xyzz[k] * ab_y + g_xz_0_y_xyyzz[k];

                g_xz_0_yy_xzzz[k] = -g_xz_0_y_xzzz[k] * ab_y + g_xz_0_y_xyzzz[k];

                g_xz_0_yy_yyyy[k] = -g_xz_0_y_yyyy[k] * ab_y + g_xz_0_y_yyyyy[k];

                g_xz_0_yy_yyyz[k] = -g_xz_0_y_yyyz[k] * ab_y + g_xz_0_y_yyyyz[k];

                g_xz_0_yy_yyzz[k] = -g_xz_0_y_yyzz[k] * ab_y + g_xz_0_y_yyyzz[k];

                g_xz_0_yy_yzzz[k] = -g_xz_0_y_yzzz[k] * ab_y + g_xz_0_y_yyzzz[k];

                g_xz_0_yy_zzzz[k] = -g_xz_0_y_zzzz[k] * ab_y + g_xz_0_y_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yz_xxxx, g_xz_0_yz_xxxy, g_xz_0_yz_xxxz, g_xz_0_yz_xxyy, g_xz_0_yz_xxyz, g_xz_0_yz_xxzz, g_xz_0_yz_xyyy, g_xz_0_yz_xyyz, g_xz_0_yz_xyzz, g_xz_0_yz_xzzz, g_xz_0_yz_yyyy, g_xz_0_yz_yyyz, g_xz_0_yz_yyzz, g_xz_0_yz_yzzz, g_xz_0_yz_zzzz, g_xz_0_z_xxxx, g_xz_0_z_xxxxy, g_xz_0_z_xxxy, g_xz_0_z_xxxyy, g_xz_0_z_xxxyz, g_xz_0_z_xxxz, g_xz_0_z_xxyy, g_xz_0_z_xxyyy, g_xz_0_z_xxyyz, g_xz_0_z_xxyz, g_xz_0_z_xxyzz, g_xz_0_z_xxzz, g_xz_0_z_xyyy, g_xz_0_z_xyyyy, g_xz_0_z_xyyyz, g_xz_0_z_xyyz, g_xz_0_z_xyyzz, g_xz_0_z_xyzz, g_xz_0_z_xyzzz, g_xz_0_z_xzzz, g_xz_0_z_yyyy, g_xz_0_z_yyyyy, g_xz_0_z_yyyyz, g_xz_0_z_yyyz, g_xz_0_z_yyyzz, g_xz_0_z_yyzz, g_xz_0_z_yyzzz, g_xz_0_z_yzzz, g_xz_0_z_yzzzz, g_xz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yz_xxxx[k] = -g_xz_0_z_xxxx[k] * ab_y + g_xz_0_z_xxxxy[k];

                g_xz_0_yz_xxxy[k] = -g_xz_0_z_xxxy[k] * ab_y + g_xz_0_z_xxxyy[k];

                g_xz_0_yz_xxxz[k] = -g_xz_0_z_xxxz[k] * ab_y + g_xz_0_z_xxxyz[k];

                g_xz_0_yz_xxyy[k] = -g_xz_0_z_xxyy[k] * ab_y + g_xz_0_z_xxyyy[k];

                g_xz_0_yz_xxyz[k] = -g_xz_0_z_xxyz[k] * ab_y + g_xz_0_z_xxyyz[k];

                g_xz_0_yz_xxzz[k] = -g_xz_0_z_xxzz[k] * ab_y + g_xz_0_z_xxyzz[k];

                g_xz_0_yz_xyyy[k] = -g_xz_0_z_xyyy[k] * ab_y + g_xz_0_z_xyyyy[k];

                g_xz_0_yz_xyyz[k] = -g_xz_0_z_xyyz[k] * ab_y + g_xz_0_z_xyyyz[k];

                g_xz_0_yz_xyzz[k] = -g_xz_0_z_xyzz[k] * ab_y + g_xz_0_z_xyyzz[k];

                g_xz_0_yz_xzzz[k] = -g_xz_0_z_xzzz[k] * ab_y + g_xz_0_z_xyzzz[k];

                g_xz_0_yz_yyyy[k] = -g_xz_0_z_yyyy[k] * ab_y + g_xz_0_z_yyyyy[k];

                g_xz_0_yz_yyyz[k] = -g_xz_0_z_yyyz[k] * ab_y + g_xz_0_z_yyyyz[k];

                g_xz_0_yz_yyzz[k] = -g_xz_0_z_yyzz[k] * ab_y + g_xz_0_z_yyyzz[k];

                g_xz_0_yz_yzzz[k] = -g_xz_0_z_yzzz[k] * ab_y + g_xz_0_z_yyzzz[k];

                g_xz_0_yz_zzzz[k] = -g_xz_0_z_zzzz[k] * ab_y + g_xz_0_z_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxx, g_x_0_z_xxxy, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyzz, g_x_0_z_yzzz, g_x_0_z_zzzz, g_xz_0_z_xxxx, g_xz_0_z_xxxxz, g_xz_0_z_xxxy, g_xz_0_z_xxxyz, g_xz_0_z_xxxz, g_xz_0_z_xxxzz, g_xz_0_z_xxyy, g_xz_0_z_xxyyz, g_xz_0_z_xxyz, g_xz_0_z_xxyzz, g_xz_0_z_xxzz, g_xz_0_z_xxzzz, g_xz_0_z_xyyy, g_xz_0_z_xyyyz, g_xz_0_z_xyyz, g_xz_0_z_xyyzz, g_xz_0_z_xyzz, g_xz_0_z_xyzzz, g_xz_0_z_xzzz, g_xz_0_z_xzzzz, g_xz_0_z_yyyy, g_xz_0_z_yyyyz, g_xz_0_z_yyyz, g_xz_0_z_yyyzz, g_xz_0_z_yyzz, g_xz_0_z_yyzzz, g_xz_0_z_yzzz, g_xz_0_z_yzzzz, g_xz_0_z_zzzz, g_xz_0_z_zzzzz, g_xz_0_zz_xxxx, g_xz_0_zz_xxxy, g_xz_0_zz_xxxz, g_xz_0_zz_xxyy, g_xz_0_zz_xxyz, g_xz_0_zz_xxzz, g_xz_0_zz_xyyy, g_xz_0_zz_xyyz, g_xz_0_zz_xyzz, g_xz_0_zz_xzzz, g_xz_0_zz_yyyy, g_xz_0_zz_yyyz, g_xz_0_zz_yyzz, g_xz_0_zz_yzzz, g_xz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zz_xxxx[k] = -g_x_0_z_xxxx[k] - g_xz_0_z_xxxx[k] * ab_z + g_xz_0_z_xxxxz[k];

                g_xz_0_zz_xxxy[k] = -g_x_0_z_xxxy[k] - g_xz_0_z_xxxy[k] * ab_z + g_xz_0_z_xxxyz[k];

                g_xz_0_zz_xxxz[k] = -g_x_0_z_xxxz[k] - g_xz_0_z_xxxz[k] * ab_z + g_xz_0_z_xxxzz[k];

                g_xz_0_zz_xxyy[k] = -g_x_0_z_xxyy[k] - g_xz_0_z_xxyy[k] * ab_z + g_xz_0_z_xxyyz[k];

                g_xz_0_zz_xxyz[k] = -g_x_0_z_xxyz[k] - g_xz_0_z_xxyz[k] * ab_z + g_xz_0_z_xxyzz[k];

                g_xz_0_zz_xxzz[k] = -g_x_0_z_xxzz[k] - g_xz_0_z_xxzz[k] * ab_z + g_xz_0_z_xxzzz[k];

                g_xz_0_zz_xyyy[k] = -g_x_0_z_xyyy[k] - g_xz_0_z_xyyy[k] * ab_z + g_xz_0_z_xyyyz[k];

                g_xz_0_zz_xyyz[k] = -g_x_0_z_xyyz[k] - g_xz_0_z_xyyz[k] * ab_z + g_xz_0_z_xyyzz[k];

                g_xz_0_zz_xyzz[k] = -g_x_0_z_xyzz[k] - g_xz_0_z_xyzz[k] * ab_z + g_xz_0_z_xyzzz[k];

                g_xz_0_zz_xzzz[k] = -g_x_0_z_xzzz[k] - g_xz_0_z_xzzz[k] * ab_z + g_xz_0_z_xzzzz[k];

                g_xz_0_zz_yyyy[k] = -g_x_0_z_yyyy[k] - g_xz_0_z_yyyy[k] * ab_z + g_xz_0_z_yyyyz[k];

                g_xz_0_zz_yyyz[k] = -g_x_0_z_yyyz[k] - g_xz_0_z_yyyz[k] * ab_z + g_xz_0_z_yyyzz[k];

                g_xz_0_zz_yyzz[k] = -g_x_0_z_yyzz[k] - g_xz_0_z_yyzz[k] * ab_z + g_xz_0_z_yyzzz[k];

                g_xz_0_zz_yzzz[k] = -g_x_0_z_yzzz[k] - g_xz_0_z_yzzz[k] * ab_z + g_xz_0_z_yzzzz[k];

                g_xz_0_zz_zzzz[k] = -g_x_0_z_zzzz[k] - g_xz_0_z_zzzz[k] * ab_z + g_xz_0_z_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 275 * ccomps * dcomps);

            auto g_yy_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 281 * ccomps * dcomps);

            auto g_yy_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_x_xxxx, g_yy_0_x_xxxxx, g_yy_0_x_xxxxy, g_yy_0_x_xxxxz, g_yy_0_x_xxxy, g_yy_0_x_xxxyy, g_yy_0_x_xxxyz, g_yy_0_x_xxxz, g_yy_0_x_xxxzz, g_yy_0_x_xxyy, g_yy_0_x_xxyyy, g_yy_0_x_xxyyz, g_yy_0_x_xxyz, g_yy_0_x_xxyzz, g_yy_0_x_xxzz, g_yy_0_x_xxzzz, g_yy_0_x_xyyy, g_yy_0_x_xyyyy, g_yy_0_x_xyyyz, g_yy_0_x_xyyz, g_yy_0_x_xyyzz, g_yy_0_x_xyzz, g_yy_0_x_xyzzz, g_yy_0_x_xzzz, g_yy_0_x_xzzzz, g_yy_0_x_yyyy, g_yy_0_x_yyyz, g_yy_0_x_yyzz, g_yy_0_x_yzzz, g_yy_0_x_zzzz, g_yy_0_xx_xxxx, g_yy_0_xx_xxxy, g_yy_0_xx_xxxz, g_yy_0_xx_xxyy, g_yy_0_xx_xxyz, g_yy_0_xx_xxzz, g_yy_0_xx_xyyy, g_yy_0_xx_xyyz, g_yy_0_xx_xyzz, g_yy_0_xx_xzzz, g_yy_0_xx_yyyy, g_yy_0_xx_yyyz, g_yy_0_xx_yyzz, g_yy_0_xx_yzzz, g_yy_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xx_xxxx[k] = -g_yy_0_x_xxxx[k] * ab_x + g_yy_0_x_xxxxx[k];

                g_yy_0_xx_xxxy[k] = -g_yy_0_x_xxxy[k] * ab_x + g_yy_0_x_xxxxy[k];

                g_yy_0_xx_xxxz[k] = -g_yy_0_x_xxxz[k] * ab_x + g_yy_0_x_xxxxz[k];

                g_yy_0_xx_xxyy[k] = -g_yy_0_x_xxyy[k] * ab_x + g_yy_0_x_xxxyy[k];

                g_yy_0_xx_xxyz[k] = -g_yy_0_x_xxyz[k] * ab_x + g_yy_0_x_xxxyz[k];

                g_yy_0_xx_xxzz[k] = -g_yy_0_x_xxzz[k] * ab_x + g_yy_0_x_xxxzz[k];

                g_yy_0_xx_xyyy[k] = -g_yy_0_x_xyyy[k] * ab_x + g_yy_0_x_xxyyy[k];

                g_yy_0_xx_xyyz[k] = -g_yy_0_x_xyyz[k] * ab_x + g_yy_0_x_xxyyz[k];

                g_yy_0_xx_xyzz[k] = -g_yy_0_x_xyzz[k] * ab_x + g_yy_0_x_xxyzz[k];

                g_yy_0_xx_xzzz[k] = -g_yy_0_x_xzzz[k] * ab_x + g_yy_0_x_xxzzz[k];

                g_yy_0_xx_yyyy[k] = -g_yy_0_x_yyyy[k] * ab_x + g_yy_0_x_xyyyy[k];

                g_yy_0_xx_yyyz[k] = -g_yy_0_x_yyyz[k] * ab_x + g_yy_0_x_xyyyz[k];

                g_yy_0_xx_yyzz[k] = -g_yy_0_x_yyzz[k] * ab_x + g_yy_0_x_xyyzz[k];

                g_yy_0_xx_yzzz[k] = -g_yy_0_x_yzzz[k] * ab_x + g_yy_0_x_xyzzz[k];

                g_yy_0_xx_zzzz[k] = -g_yy_0_x_zzzz[k] * ab_x + g_yy_0_x_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 287 * ccomps * dcomps);

            auto g_yy_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 293 * ccomps * dcomps);

            auto g_yy_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xy_xxxx, g_yy_0_xy_xxxy, g_yy_0_xy_xxxz, g_yy_0_xy_xxyy, g_yy_0_xy_xxyz, g_yy_0_xy_xxzz, g_yy_0_xy_xyyy, g_yy_0_xy_xyyz, g_yy_0_xy_xyzz, g_yy_0_xy_xzzz, g_yy_0_xy_yyyy, g_yy_0_xy_yyyz, g_yy_0_xy_yyzz, g_yy_0_xy_yzzz, g_yy_0_xy_zzzz, g_yy_0_y_xxxx, g_yy_0_y_xxxxx, g_yy_0_y_xxxxy, g_yy_0_y_xxxxz, g_yy_0_y_xxxy, g_yy_0_y_xxxyy, g_yy_0_y_xxxyz, g_yy_0_y_xxxz, g_yy_0_y_xxxzz, g_yy_0_y_xxyy, g_yy_0_y_xxyyy, g_yy_0_y_xxyyz, g_yy_0_y_xxyz, g_yy_0_y_xxyzz, g_yy_0_y_xxzz, g_yy_0_y_xxzzz, g_yy_0_y_xyyy, g_yy_0_y_xyyyy, g_yy_0_y_xyyyz, g_yy_0_y_xyyz, g_yy_0_y_xyyzz, g_yy_0_y_xyzz, g_yy_0_y_xyzzz, g_yy_0_y_xzzz, g_yy_0_y_xzzzz, g_yy_0_y_yyyy, g_yy_0_y_yyyz, g_yy_0_y_yyzz, g_yy_0_y_yzzz, g_yy_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xy_xxxx[k] = -g_yy_0_y_xxxx[k] * ab_x + g_yy_0_y_xxxxx[k];

                g_yy_0_xy_xxxy[k] = -g_yy_0_y_xxxy[k] * ab_x + g_yy_0_y_xxxxy[k];

                g_yy_0_xy_xxxz[k] = -g_yy_0_y_xxxz[k] * ab_x + g_yy_0_y_xxxxz[k];

                g_yy_0_xy_xxyy[k] = -g_yy_0_y_xxyy[k] * ab_x + g_yy_0_y_xxxyy[k];

                g_yy_0_xy_xxyz[k] = -g_yy_0_y_xxyz[k] * ab_x + g_yy_0_y_xxxyz[k];

                g_yy_0_xy_xxzz[k] = -g_yy_0_y_xxzz[k] * ab_x + g_yy_0_y_xxxzz[k];

                g_yy_0_xy_xyyy[k] = -g_yy_0_y_xyyy[k] * ab_x + g_yy_0_y_xxyyy[k];

                g_yy_0_xy_xyyz[k] = -g_yy_0_y_xyyz[k] * ab_x + g_yy_0_y_xxyyz[k];

                g_yy_0_xy_xyzz[k] = -g_yy_0_y_xyzz[k] * ab_x + g_yy_0_y_xxyzz[k];

                g_yy_0_xy_xzzz[k] = -g_yy_0_y_xzzz[k] * ab_x + g_yy_0_y_xxzzz[k];

                g_yy_0_xy_yyyy[k] = -g_yy_0_y_yyyy[k] * ab_x + g_yy_0_y_xyyyy[k];

                g_yy_0_xy_yyyz[k] = -g_yy_0_y_yyyz[k] * ab_x + g_yy_0_y_xyyyz[k];

                g_yy_0_xy_yyzz[k] = -g_yy_0_y_yyzz[k] * ab_x + g_yy_0_y_xyyzz[k];

                g_yy_0_xy_yzzz[k] = -g_yy_0_y_yzzz[k] * ab_x + g_yy_0_y_xyzzz[k];

                g_yy_0_xy_zzzz[k] = -g_yy_0_y_zzzz[k] * ab_x + g_yy_0_y_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xz_xxxx, g_yy_0_xz_xxxy, g_yy_0_xz_xxxz, g_yy_0_xz_xxyy, g_yy_0_xz_xxyz, g_yy_0_xz_xxzz, g_yy_0_xz_xyyy, g_yy_0_xz_xyyz, g_yy_0_xz_xyzz, g_yy_0_xz_xzzz, g_yy_0_xz_yyyy, g_yy_0_xz_yyyz, g_yy_0_xz_yyzz, g_yy_0_xz_yzzz, g_yy_0_xz_zzzz, g_yy_0_z_xxxx, g_yy_0_z_xxxxx, g_yy_0_z_xxxxy, g_yy_0_z_xxxxz, g_yy_0_z_xxxy, g_yy_0_z_xxxyy, g_yy_0_z_xxxyz, g_yy_0_z_xxxz, g_yy_0_z_xxxzz, g_yy_0_z_xxyy, g_yy_0_z_xxyyy, g_yy_0_z_xxyyz, g_yy_0_z_xxyz, g_yy_0_z_xxyzz, g_yy_0_z_xxzz, g_yy_0_z_xxzzz, g_yy_0_z_xyyy, g_yy_0_z_xyyyy, g_yy_0_z_xyyyz, g_yy_0_z_xyyz, g_yy_0_z_xyyzz, g_yy_0_z_xyzz, g_yy_0_z_xyzzz, g_yy_0_z_xzzz, g_yy_0_z_xzzzz, g_yy_0_z_yyyy, g_yy_0_z_yyyz, g_yy_0_z_yyzz, g_yy_0_z_yzzz, g_yy_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xz_xxxx[k] = -g_yy_0_z_xxxx[k] * ab_x + g_yy_0_z_xxxxx[k];

                g_yy_0_xz_xxxy[k] = -g_yy_0_z_xxxy[k] * ab_x + g_yy_0_z_xxxxy[k];

                g_yy_0_xz_xxxz[k] = -g_yy_0_z_xxxz[k] * ab_x + g_yy_0_z_xxxxz[k];

                g_yy_0_xz_xxyy[k] = -g_yy_0_z_xxyy[k] * ab_x + g_yy_0_z_xxxyy[k];

                g_yy_0_xz_xxyz[k] = -g_yy_0_z_xxyz[k] * ab_x + g_yy_0_z_xxxyz[k];

                g_yy_0_xz_xxzz[k] = -g_yy_0_z_xxzz[k] * ab_x + g_yy_0_z_xxxzz[k];

                g_yy_0_xz_xyyy[k] = -g_yy_0_z_xyyy[k] * ab_x + g_yy_0_z_xxyyy[k];

                g_yy_0_xz_xyyz[k] = -g_yy_0_z_xyyz[k] * ab_x + g_yy_0_z_xxyyz[k];

                g_yy_0_xz_xyzz[k] = -g_yy_0_z_xyzz[k] * ab_x + g_yy_0_z_xxyzz[k];

                g_yy_0_xz_xzzz[k] = -g_yy_0_z_xzzz[k] * ab_x + g_yy_0_z_xxzzz[k];

                g_yy_0_xz_yyyy[k] = -g_yy_0_z_yyyy[k] * ab_x + g_yy_0_z_xyyyy[k];

                g_yy_0_xz_yyyz[k] = -g_yy_0_z_yyyz[k] * ab_x + g_yy_0_z_xyyyz[k];

                g_yy_0_xz_yyzz[k] = -g_yy_0_z_yyzz[k] * ab_x + g_yy_0_z_xyyzz[k];

                g_yy_0_xz_yzzz[k] = -g_yy_0_z_yzzz[k] * ab_x + g_yy_0_z_xyzzz[k];

                g_yy_0_xz_zzzz[k] = -g_yy_0_z_zzzz[k] * ab_x + g_yy_0_z_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz, g_yy_0_y_xxxx, g_yy_0_y_xxxxy, g_yy_0_y_xxxy, g_yy_0_y_xxxyy, g_yy_0_y_xxxyz, g_yy_0_y_xxxz, g_yy_0_y_xxyy, g_yy_0_y_xxyyy, g_yy_0_y_xxyyz, g_yy_0_y_xxyz, g_yy_0_y_xxyzz, g_yy_0_y_xxzz, g_yy_0_y_xyyy, g_yy_0_y_xyyyy, g_yy_0_y_xyyyz, g_yy_0_y_xyyz, g_yy_0_y_xyyzz, g_yy_0_y_xyzz, g_yy_0_y_xyzzz, g_yy_0_y_xzzz, g_yy_0_y_yyyy, g_yy_0_y_yyyyy, g_yy_0_y_yyyyz, g_yy_0_y_yyyz, g_yy_0_y_yyyzz, g_yy_0_y_yyzz, g_yy_0_y_yyzzz, g_yy_0_y_yzzz, g_yy_0_y_yzzzz, g_yy_0_y_zzzz, g_yy_0_yy_xxxx, g_yy_0_yy_xxxy, g_yy_0_yy_xxxz, g_yy_0_yy_xxyy, g_yy_0_yy_xxyz, g_yy_0_yy_xxzz, g_yy_0_yy_xyyy, g_yy_0_yy_xyyz, g_yy_0_yy_xyzz, g_yy_0_yy_xzzz, g_yy_0_yy_yyyy, g_yy_0_yy_yyyz, g_yy_0_yy_yyzz, g_yy_0_yy_yzzz, g_yy_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yy_xxxx[k] = -2.0 * g_y_0_y_xxxx[k] - g_yy_0_y_xxxx[k] * ab_y + g_yy_0_y_xxxxy[k];

                g_yy_0_yy_xxxy[k] = -2.0 * g_y_0_y_xxxy[k] - g_yy_0_y_xxxy[k] * ab_y + g_yy_0_y_xxxyy[k];

                g_yy_0_yy_xxxz[k] = -2.0 * g_y_0_y_xxxz[k] - g_yy_0_y_xxxz[k] * ab_y + g_yy_0_y_xxxyz[k];

                g_yy_0_yy_xxyy[k] = -2.0 * g_y_0_y_xxyy[k] - g_yy_0_y_xxyy[k] * ab_y + g_yy_0_y_xxyyy[k];

                g_yy_0_yy_xxyz[k] = -2.0 * g_y_0_y_xxyz[k] - g_yy_0_y_xxyz[k] * ab_y + g_yy_0_y_xxyyz[k];

                g_yy_0_yy_xxzz[k] = -2.0 * g_y_0_y_xxzz[k] - g_yy_0_y_xxzz[k] * ab_y + g_yy_0_y_xxyzz[k];

                g_yy_0_yy_xyyy[k] = -2.0 * g_y_0_y_xyyy[k] - g_yy_0_y_xyyy[k] * ab_y + g_yy_0_y_xyyyy[k];

                g_yy_0_yy_xyyz[k] = -2.0 * g_y_0_y_xyyz[k] - g_yy_0_y_xyyz[k] * ab_y + g_yy_0_y_xyyyz[k];

                g_yy_0_yy_xyzz[k] = -2.0 * g_y_0_y_xyzz[k] - g_yy_0_y_xyzz[k] * ab_y + g_yy_0_y_xyyzz[k];

                g_yy_0_yy_xzzz[k] = -2.0 * g_y_0_y_xzzz[k] - g_yy_0_y_xzzz[k] * ab_y + g_yy_0_y_xyzzz[k];

                g_yy_0_yy_yyyy[k] = -2.0 * g_y_0_y_yyyy[k] - g_yy_0_y_yyyy[k] * ab_y + g_yy_0_y_yyyyy[k];

                g_yy_0_yy_yyyz[k] = -2.0 * g_y_0_y_yyyz[k] - g_yy_0_y_yyyz[k] * ab_y + g_yy_0_y_yyyyz[k];

                g_yy_0_yy_yyzz[k] = -2.0 * g_y_0_y_yyzz[k] - g_yy_0_y_yyzz[k] * ab_y + g_yy_0_y_yyyzz[k];

                g_yy_0_yy_yzzz[k] = -2.0 * g_y_0_y_yzzz[k] - g_yy_0_y_yzzz[k] * ab_y + g_yy_0_y_yyzzz[k];

                g_yy_0_yy_zzzz[k] = -2.0 * g_y_0_y_zzzz[k] - g_yy_0_y_zzzz[k] * ab_y + g_yy_0_y_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_y_xxxx, g_yy_0_y_xxxxz, g_yy_0_y_xxxy, g_yy_0_y_xxxyz, g_yy_0_y_xxxz, g_yy_0_y_xxxzz, g_yy_0_y_xxyy, g_yy_0_y_xxyyz, g_yy_0_y_xxyz, g_yy_0_y_xxyzz, g_yy_0_y_xxzz, g_yy_0_y_xxzzz, g_yy_0_y_xyyy, g_yy_0_y_xyyyz, g_yy_0_y_xyyz, g_yy_0_y_xyyzz, g_yy_0_y_xyzz, g_yy_0_y_xyzzz, g_yy_0_y_xzzz, g_yy_0_y_xzzzz, g_yy_0_y_yyyy, g_yy_0_y_yyyyz, g_yy_0_y_yyyz, g_yy_0_y_yyyzz, g_yy_0_y_yyzz, g_yy_0_y_yyzzz, g_yy_0_y_yzzz, g_yy_0_y_yzzzz, g_yy_0_y_zzzz, g_yy_0_y_zzzzz, g_yy_0_yz_xxxx, g_yy_0_yz_xxxy, g_yy_0_yz_xxxz, g_yy_0_yz_xxyy, g_yy_0_yz_xxyz, g_yy_0_yz_xxzz, g_yy_0_yz_xyyy, g_yy_0_yz_xyyz, g_yy_0_yz_xyzz, g_yy_0_yz_xzzz, g_yy_0_yz_yyyy, g_yy_0_yz_yyyz, g_yy_0_yz_yyzz, g_yy_0_yz_yzzz, g_yy_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yz_xxxx[k] = -g_yy_0_y_xxxx[k] * ab_z + g_yy_0_y_xxxxz[k];

                g_yy_0_yz_xxxy[k] = -g_yy_0_y_xxxy[k] * ab_z + g_yy_0_y_xxxyz[k];

                g_yy_0_yz_xxxz[k] = -g_yy_0_y_xxxz[k] * ab_z + g_yy_0_y_xxxzz[k];

                g_yy_0_yz_xxyy[k] = -g_yy_0_y_xxyy[k] * ab_z + g_yy_0_y_xxyyz[k];

                g_yy_0_yz_xxyz[k] = -g_yy_0_y_xxyz[k] * ab_z + g_yy_0_y_xxyzz[k];

                g_yy_0_yz_xxzz[k] = -g_yy_0_y_xxzz[k] * ab_z + g_yy_0_y_xxzzz[k];

                g_yy_0_yz_xyyy[k] = -g_yy_0_y_xyyy[k] * ab_z + g_yy_0_y_xyyyz[k];

                g_yy_0_yz_xyyz[k] = -g_yy_0_y_xyyz[k] * ab_z + g_yy_0_y_xyyzz[k];

                g_yy_0_yz_xyzz[k] = -g_yy_0_y_xyzz[k] * ab_z + g_yy_0_y_xyzzz[k];

                g_yy_0_yz_xzzz[k] = -g_yy_0_y_xzzz[k] * ab_z + g_yy_0_y_xzzzz[k];

                g_yy_0_yz_yyyy[k] = -g_yy_0_y_yyyy[k] * ab_z + g_yy_0_y_yyyyz[k];

                g_yy_0_yz_yyyz[k] = -g_yy_0_y_yyyz[k] * ab_z + g_yy_0_y_yyyzz[k];

                g_yy_0_yz_yyzz[k] = -g_yy_0_y_yyzz[k] * ab_z + g_yy_0_y_yyzzz[k];

                g_yy_0_yz_yzzz[k] = -g_yy_0_y_yzzz[k] * ab_z + g_yy_0_y_yzzzz[k];

                g_yy_0_yz_zzzz[k] = -g_yy_0_y_zzzz[k] * ab_z + g_yy_0_y_zzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_z_xxxx, g_yy_0_z_xxxxz, g_yy_0_z_xxxy, g_yy_0_z_xxxyz, g_yy_0_z_xxxz, g_yy_0_z_xxxzz, g_yy_0_z_xxyy, g_yy_0_z_xxyyz, g_yy_0_z_xxyz, g_yy_0_z_xxyzz, g_yy_0_z_xxzz, g_yy_0_z_xxzzz, g_yy_0_z_xyyy, g_yy_0_z_xyyyz, g_yy_0_z_xyyz, g_yy_0_z_xyyzz, g_yy_0_z_xyzz, g_yy_0_z_xyzzz, g_yy_0_z_xzzz, g_yy_0_z_xzzzz, g_yy_0_z_yyyy, g_yy_0_z_yyyyz, g_yy_0_z_yyyz, g_yy_0_z_yyyzz, g_yy_0_z_yyzz, g_yy_0_z_yyzzz, g_yy_0_z_yzzz, g_yy_0_z_yzzzz, g_yy_0_z_zzzz, g_yy_0_z_zzzzz, g_yy_0_zz_xxxx, g_yy_0_zz_xxxy, g_yy_0_zz_xxxz, g_yy_0_zz_xxyy, g_yy_0_zz_xxyz, g_yy_0_zz_xxzz, g_yy_0_zz_xyyy, g_yy_0_zz_xyyz, g_yy_0_zz_xyzz, g_yy_0_zz_xzzz, g_yy_0_zz_yyyy, g_yy_0_zz_yyyz, g_yy_0_zz_yyzz, g_yy_0_zz_yzzz, g_yy_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zz_xxxx[k] = -g_yy_0_z_xxxx[k] * ab_z + g_yy_0_z_xxxxz[k];

                g_yy_0_zz_xxxy[k] = -g_yy_0_z_xxxy[k] * ab_z + g_yy_0_z_xxxyz[k];

                g_yy_0_zz_xxxz[k] = -g_yy_0_z_xxxz[k] * ab_z + g_yy_0_z_xxxzz[k];

                g_yy_0_zz_xxyy[k] = -g_yy_0_z_xxyy[k] * ab_z + g_yy_0_z_xxyyz[k];

                g_yy_0_zz_xxyz[k] = -g_yy_0_z_xxyz[k] * ab_z + g_yy_0_z_xxyzz[k];

                g_yy_0_zz_xxzz[k] = -g_yy_0_z_xxzz[k] * ab_z + g_yy_0_z_xxzzz[k];

                g_yy_0_zz_xyyy[k] = -g_yy_0_z_xyyy[k] * ab_z + g_yy_0_z_xyyyz[k];

                g_yy_0_zz_xyyz[k] = -g_yy_0_z_xyyz[k] * ab_z + g_yy_0_z_xyyzz[k];

                g_yy_0_zz_xyzz[k] = -g_yy_0_z_xyzz[k] * ab_z + g_yy_0_z_xyzzz[k];

                g_yy_0_zz_xzzz[k] = -g_yy_0_z_xzzz[k] * ab_z + g_yy_0_z_xzzzz[k];

                g_yy_0_zz_yyyy[k] = -g_yy_0_z_yyyy[k] * ab_z + g_yy_0_z_yyyyz[k];

                g_yy_0_zz_yyyz[k] = -g_yy_0_z_yyyz[k] * ab_z + g_yy_0_z_yyyzz[k];

                g_yy_0_zz_yyzz[k] = -g_yy_0_z_yyzz[k] * ab_z + g_yy_0_z_yyzzz[k];

                g_yy_0_zz_yzzz[k] = -g_yy_0_z_yzzz[k] * ab_z + g_yy_0_z_yzzzz[k];

                g_yy_0_zz_zzzz[k] = -g_yy_0_z_zzzz[k] * ab_z + g_yy_0_z_zzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 365 * ccomps * dcomps);

            auto g_yz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 371 * ccomps * dcomps);

            auto g_yz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_x_xxxx, g_yz_0_x_xxxxx, g_yz_0_x_xxxxy, g_yz_0_x_xxxxz, g_yz_0_x_xxxy, g_yz_0_x_xxxyy, g_yz_0_x_xxxyz, g_yz_0_x_xxxz, g_yz_0_x_xxxzz, g_yz_0_x_xxyy, g_yz_0_x_xxyyy, g_yz_0_x_xxyyz, g_yz_0_x_xxyz, g_yz_0_x_xxyzz, g_yz_0_x_xxzz, g_yz_0_x_xxzzz, g_yz_0_x_xyyy, g_yz_0_x_xyyyy, g_yz_0_x_xyyyz, g_yz_0_x_xyyz, g_yz_0_x_xyyzz, g_yz_0_x_xyzz, g_yz_0_x_xyzzz, g_yz_0_x_xzzz, g_yz_0_x_xzzzz, g_yz_0_x_yyyy, g_yz_0_x_yyyz, g_yz_0_x_yyzz, g_yz_0_x_yzzz, g_yz_0_x_zzzz, g_yz_0_xx_xxxx, g_yz_0_xx_xxxy, g_yz_0_xx_xxxz, g_yz_0_xx_xxyy, g_yz_0_xx_xxyz, g_yz_0_xx_xxzz, g_yz_0_xx_xyyy, g_yz_0_xx_xyyz, g_yz_0_xx_xyzz, g_yz_0_xx_xzzz, g_yz_0_xx_yyyy, g_yz_0_xx_yyyz, g_yz_0_xx_yyzz, g_yz_0_xx_yzzz, g_yz_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xx_xxxx[k] = -g_yz_0_x_xxxx[k] * ab_x + g_yz_0_x_xxxxx[k];

                g_yz_0_xx_xxxy[k] = -g_yz_0_x_xxxy[k] * ab_x + g_yz_0_x_xxxxy[k];

                g_yz_0_xx_xxxz[k] = -g_yz_0_x_xxxz[k] * ab_x + g_yz_0_x_xxxxz[k];

                g_yz_0_xx_xxyy[k] = -g_yz_0_x_xxyy[k] * ab_x + g_yz_0_x_xxxyy[k];

                g_yz_0_xx_xxyz[k] = -g_yz_0_x_xxyz[k] * ab_x + g_yz_0_x_xxxyz[k];

                g_yz_0_xx_xxzz[k] = -g_yz_0_x_xxzz[k] * ab_x + g_yz_0_x_xxxzz[k];

                g_yz_0_xx_xyyy[k] = -g_yz_0_x_xyyy[k] * ab_x + g_yz_0_x_xxyyy[k];

                g_yz_0_xx_xyyz[k] = -g_yz_0_x_xyyz[k] * ab_x + g_yz_0_x_xxyyz[k];

                g_yz_0_xx_xyzz[k] = -g_yz_0_x_xyzz[k] * ab_x + g_yz_0_x_xxyzz[k];

                g_yz_0_xx_xzzz[k] = -g_yz_0_x_xzzz[k] * ab_x + g_yz_0_x_xxzzz[k];

                g_yz_0_xx_yyyy[k] = -g_yz_0_x_yyyy[k] * ab_x + g_yz_0_x_xyyyy[k];

                g_yz_0_xx_yyyz[k] = -g_yz_0_x_yyyz[k] * ab_x + g_yz_0_x_xyyyz[k];

                g_yz_0_xx_yyzz[k] = -g_yz_0_x_yyzz[k] * ab_x + g_yz_0_x_xyyzz[k];

                g_yz_0_xx_yzzz[k] = -g_yz_0_x_yzzz[k] * ab_x + g_yz_0_x_xyzzz[k];

                g_yz_0_xx_zzzz[k] = -g_yz_0_x_zzzz[k] * ab_x + g_yz_0_x_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 377 * ccomps * dcomps);

            auto g_yz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 383 * ccomps * dcomps);

            auto g_yz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xy_xxxx, g_yz_0_xy_xxxy, g_yz_0_xy_xxxz, g_yz_0_xy_xxyy, g_yz_0_xy_xxyz, g_yz_0_xy_xxzz, g_yz_0_xy_xyyy, g_yz_0_xy_xyyz, g_yz_0_xy_xyzz, g_yz_0_xy_xzzz, g_yz_0_xy_yyyy, g_yz_0_xy_yyyz, g_yz_0_xy_yyzz, g_yz_0_xy_yzzz, g_yz_0_xy_zzzz, g_yz_0_y_xxxx, g_yz_0_y_xxxxx, g_yz_0_y_xxxxy, g_yz_0_y_xxxxz, g_yz_0_y_xxxy, g_yz_0_y_xxxyy, g_yz_0_y_xxxyz, g_yz_0_y_xxxz, g_yz_0_y_xxxzz, g_yz_0_y_xxyy, g_yz_0_y_xxyyy, g_yz_0_y_xxyyz, g_yz_0_y_xxyz, g_yz_0_y_xxyzz, g_yz_0_y_xxzz, g_yz_0_y_xxzzz, g_yz_0_y_xyyy, g_yz_0_y_xyyyy, g_yz_0_y_xyyyz, g_yz_0_y_xyyz, g_yz_0_y_xyyzz, g_yz_0_y_xyzz, g_yz_0_y_xyzzz, g_yz_0_y_xzzz, g_yz_0_y_xzzzz, g_yz_0_y_yyyy, g_yz_0_y_yyyz, g_yz_0_y_yyzz, g_yz_0_y_yzzz, g_yz_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xy_xxxx[k] = -g_yz_0_y_xxxx[k] * ab_x + g_yz_0_y_xxxxx[k];

                g_yz_0_xy_xxxy[k] = -g_yz_0_y_xxxy[k] * ab_x + g_yz_0_y_xxxxy[k];

                g_yz_0_xy_xxxz[k] = -g_yz_0_y_xxxz[k] * ab_x + g_yz_0_y_xxxxz[k];

                g_yz_0_xy_xxyy[k] = -g_yz_0_y_xxyy[k] * ab_x + g_yz_0_y_xxxyy[k];

                g_yz_0_xy_xxyz[k] = -g_yz_0_y_xxyz[k] * ab_x + g_yz_0_y_xxxyz[k];

                g_yz_0_xy_xxzz[k] = -g_yz_0_y_xxzz[k] * ab_x + g_yz_0_y_xxxzz[k];

                g_yz_0_xy_xyyy[k] = -g_yz_0_y_xyyy[k] * ab_x + g_yz_0_y_xxyyy[k];

                g_yz_0_xy_xyyz[k] = -g_yz_0_y_xyyz[k] * ab_x + g_yz_0_y_xxyyz[k];

                g_yz_0_xy_xyzz[k] = -g_yz_0_y_xyzz[k] * ab_x + g_yz_0_y_xxyzz[k];

                g_yz_0_xy_xzzz[k] = -g_yz_0_y_xzzz[k] * ab_x + g_yz_0_y_xxzzz[k];

                g_yz_0_xy_yyyy[k] = -g_yz_0_y_yyyy[k] * ab_x + g_yz_0_y_xyyyy[k];

                g_yz_0_xy_yyyz[k] = -g_yz_0_y_yyyz[k] * ab_x + g_yz_0_y_xyyyz[k];

                g_yz_0_xy_yyzz[k] = -g_yz_0_y_yyzz[k] * ab_x + g_yz_0_y_xyyzz[k];

                g_yz_0_xy_yzzz[k] = -g_yz_0_y_yzzz[k] * ab_x + g_yz_0_y_xyzzz[k];

                g_yz_0_xy_zzzz[k] = -g_yz_0_y_zzzz[k] * ab_x + g_yz_0_y_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 395 * ccomps * dcomps);

            auto g_yz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xz_xxxx, g_yz_0_xz_xxxy, g_yz_0_xz_xxxz, g_yz_0_xz_xxyy, g_yz_0_xz_xxyz, g_yz_0_xz_xxzz, g_yz_0_xz_xyyy, g_yz_0_xz_xyyz, g_yz_0_xz_xyzz, g_yz_0_xz_xzzz, g_yz_0_xz_yyyy, g_yz_0_xz_yyyz, g_yz_0_xz_yyzz, g_yz_0_xz_yzzz, g_yz_0_xz_zzzz, g_yz_0_z_xxxx, g_yz_0_z_xxxxx, g_yz_0_z_xxxxy, g_yz_0_z_xxxxz, g_yz_0_z_xxxy, g_yz_0_z_xxxyy, g_yz_0_z_xxxyz, g_yz_0_z_xxxz, g_yz_0_z_xxxzz, g_yz_0_z_xxyy, g_yz_0_z_xxyyy, g_yz_0_z_xxyyz, g_yz_0_z_xxyz, g_yz_0_z_xxyzz, g_yz_0_z_xxzz, g_yz_0_z_xxzzz, g_yz_0_z_xyyy, g_yz_0_z_xyyyy, g_yz_0_z_xyyyz, g_yz_0_z_xyyz, g_yz_0_z_xyyzz, g_yz_0_z_xyzz, g_yz_0_z_xyzzz, g_yz_0_z_xzzz, g_yz_0_z_xzzzz, g_yz_0_z_yyyy, g_yz_0_z_yyyz, g_yz_0_z_yyzz, g_yz_0_z_yzzz, g_yz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xz_xxxx[k] = -g_yz_0_z_xxxx[k] * ab_x + g_yz_0_z_xxxxx[k];

                g_yz_0_xz_xxxy[k] = -g_yz_0_z_xxxy[k] * ab_x + g_yz_0_z_xxxxy[k];

                g_yz_0_xz_xxxz[k] = -g_yz_0_z_xxxz[k] * ab_x + g_yz_0_z_xxxxz[k];

                g_yz_0_xz_xxyy[k] = -g_yz_0_z_xxyy[k] * ab_x + g_yz_0_z_xxxyy[k];

                g_yz_0_xz_xxyz[k] = -g_yz_0_z_xxyz[k] * ab_x + g_yz_0_z_xxxyz[k];

                g_yz_0_xz_xxzz[k] = -g_yz_0_z_xxzz[k] * ab_x + g_yz_0_z_xxxzz[k];

                g_yz_0_xz_xyyy[k] = -g_yz_0_z_xyyy[k] * ab_x + g_yz_0_z_xxyyy[k];

                g_yz_0_xz_xyyz[k] = -g_yz_0_z_xyyz[k] * ab_x + g_yz_0_z_xxyyz[k];

                g_yz_0_xz_xyzz[k] = -g_yz_0_z_xyzz[k] * ab_x + g_yz_0_z_xxyzz[k];

                g_yz_0_xz_xzzz[k] = -g_yz_0_z_xzzz[k] * ab_x + g_yz_0_z_xxzzz[k];

                g_yz_0_xz_yyyy[k] = -g_yz_0_z_yyyy[k] * ab_x + g_yz_0_z_xyyyy[k];

                g_yz_0_xz_yyyz[k] = -g_yz_0_z_yyyz[k] * ab_x + g_yz_0_z_xyyyz[k];

                g_yz_0_xz_yyzz[k] = -g_yz_0_z_yyzz[k] * ab_x + g_yz_0_z_xyyzz[k];

                g_yz_0_xz_yzzz[k] = -g_yz_0_z_yzzz[k] * ab_x + g_yz_0_z_xyzzz[k];

                g_yz_0_xz_zzzz[k] = -g_yz_0_z_zzzz[k] * ab_x + g_yz_0_z_xzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_y_xxxx, g_yz_0_y_xxxxy, g_yz_0_y_xxxy, g_yz_0_y_xxxyy, g_yz_0_y_xxxyz, g_yz_0_y_xxxz, g_yz_0_y_xxyy, g_yz_0_y_xxyyy, g_yz_0_y_xxyyz, g_yz_0_y_xxyz, g_yz_0_y_xxyzz, g_yz_0_y_xxzz, g_yz_0_y_xyyy, g_yz_0_y_xyyyy, g_yz_0_y_xyyyz, g_yz_0_y_xyyz, g_yz_0_y_xyyzz, g_yz_0_y_xyzz, g_yz_0_y_xyzzz, g_yz_0_y_xzzz, g_yz_0_y_yyyy, g_yz_0_y_yyyyy, g_yz_0_y_yyyyz, g_yz_0_y_yyyz, g_yz_0_y_yyyzz, g_yz_0_y_yyzz, g_yz_0_y_yyzzz, g_yz_0_y_yzzz, g_yz_0_y_yzzzz, g_yz_0_y_zzzz, g_yz_0_yy_xxxx, g_yz_0_yy_xxxy, g_yz_0_yy_xxxz, g_yz_0_yy_xxyy, g_yz_0_yy_xxyz, g_yz_0_yy_xxzz, g_yz_0_yy_xyyy, g_yz_0_yy_xyyz, g_yz_0_yy_xyzz, g_yz_0_yy_xzzz, g_yz_0_yy_yyyy, g_yz_0_yy_yyyz, g_yz_0_yy_yyzz, g_yz_0_yy_yzzz, g_yz_0_yy_zzzz, g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yy_xxxx[k] = -g_z_0_y_xxxx[k] - g_yz_0_y_xxxx[k] * ab_y + g_yz_0_y_xxxxy[k];

                g_yz_0_yy_xxxy[k] = -g_z_0_y_xxxy[k] - g_yz_0_y_xxxy[k] * ab_y + g_yz_0_y_xxxyy[k];

                g_yz_0_yy_xxxz[k] = -g_z_0_y_xxxz[k] - g_yz_0_y_xxxz[k] * ab_y + g_yz_0_y_xxxyz[k];

                g_yz_0_yy_xxyy[k] = -g_z_0_y_xxyy[k] - g_yz_0_y_xxyy[k] * ab_y + g_yz_0_y_xxyyy[k];

                g_yz_0_yy_xxyz[k] = -g_z_0_y_xxyz[k] - g_yz_0_y_xxyz[k] * ab_y + g_yz_0_y_xxyyz[k];

                g_yz_0_yy_xxzz[k] = -g_z_0_y_xxzz[k] - g_yz_0_y_xxzz[k] * ab_y + g_yz_0_y_xxyzz[k];

                g_yz_0_yy_xyyy[k] = -g_z_0_y_xyyy[k] - g_yz_0_y_xyyy[k] * ab_y + g_yz_0_y_xyyyy[k];

                g_yz_0_yy_xyyz[k] = -g_z_0_y_xyyz[k] - g_yz_0_y_xyyz[k] * ab_y + g_yz_0_y_xyyyz[k];

                g_yz_0_yy_xyzz[k] = -g_z_0_y_xyzz[k] - g_yz_0_y_xyzz[k] * ab_y + g_yz_0_y_xyyzz[k];

                g_yz_0_yy_xzzz[k] = -g_z_0_y_xzzz[k] - g_yz_0_y_xzzz[k] * ab_y + g_yz_0_y_xyzzz[k];

                g_yz_0_yy_yyyy[k] = -g_z_0_y_yyyy[k] - g_yz_0_y_yyyy[k] * ab_y + g_yz_0_y_yyyyy[k];

                g_yz_0_yy_yyyz[k] = -g_z_0_y_yyyz[k] - g_yz_0_y_yyyz[k] * ab_y + g_yz_0_y_yyyyz[k];

                g_yz_0_yy_yyzz[k] = -g_z_0_y_yyzz[k] - g_yz_0_y_yyzz[k] * ab_y + g_yz_0_y_yyyzz[k];

                g_yz_0_yy_yzzz[k] = -g_z_0_y_yzzz[k] - g_yz_0_y_yzzz[k] * ab_y + g_yz_0_y_yyzzz[k];

                g_yz_0_yy_zzzz[k] = -g_z_0_y_zzzz[k] - g_yz_0_y_zzzz[k] * ab_y + g_yz_0_y_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yz_xxxx, g_yz_0_yz_xxxy, g_yz_0_yz_xxxz, g_yz_0_yz_xxyy, g_yz_0_yz_xxyz, g_yz_0_yz_xxzz, g_yz_0_yz_xyyy, g_yz_0_yz_xyyz, g_yz_0_yz_xyzz, g_yz_0_yz_xzzz, g_yz_0_yz_yyyy, g_yz_0_yz_yyyz, g_yz_0_yz_yyzz, g_yz_0_yz_yzzz, g_yz_0_yz_zzzz, g_yz_0_z_xxxx, g_yz_0_z_xxxxy, g_yz_0_z_xxxy, g_yz_0_z_xxxyy, g_yz_0_z_xxxyz, g_yz_0_z_xxxz, g_yz_0_z_xxyy, g_yz_0_z_xxyyy, g_yz_0_z_xxyyz, g_yz_0_z_xxyz, g_yz_0_z_xxyzz, g_yz_0_z_xxzz, g_yz_0_z_xyyy, g_yz_0_z_xyyyy, g_yz_0_z_xyyyz, g_yz_0_z_xyyz, g_yz_0_z_xyyzz, g_yz_0_z_xyzz, g_yz_0_z_xyzzz, g_yz_0_z_xzzz, g_yz_0_z_yyyy, g_yz_0_z_yyyyy, g_yz_0_z_yyyyz, g_yz_0_z_yyyz, g_yz_0_z_yyyzz, g_yz_0_z_yyzz, g_yz_0_z_yyzzz, g_yz_0_z_yzzz, g_yz_0_z_yzzzz, g_yz_0_z_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yz_xxxx[k] = -g_z_0_z_xxxx[k] - g_yz_0_z_xxxx[k] * ab_y + g_yz_0_z_xxxxy[k];

                g_yz_0_yz_xxxy[k] = -g_z_0_z_xxxy[k] - g_yz_0_z_xxxy[k] * ab_y + g_yz_0_z_xxxyy[k];

                g_yz_0_yz_xxxz[k] = -g_z_0_z_xxxz[k] - g_yz_0_z_xxxz[k] * ab_y + g_yz_0_z_xxxyz[k];

                g_yz_0_yz_xxyy[k] = -g_z_0_z_xxyy[k] - g_yz_0_z_xxyy[k] * ab_y + g_yz_0_z_xxyyy[k];

                g_yz_0_yz_xxyz[k] = -g_z_0_z_xxyz[k] - g_yz_0_z_xxyz[k] * ab_y + g_yz_0_z_xxyyz[k];

                g_yz_0_yz_xxzz[k] = -g_z_0_z_xxzz[k] - g_yz_0_z_xxzz[k] * ab_y + g_yz_0_z_xxyzz[k];

                g_yz_0_yz_xyyy[k] = -g_z_0_z_xyyy[k] - g_yz_0_z_xyyy[k] * ab_y + g_yz_0_z_xyyyy[k];

                g_yz_0_yz_xyyz[k] = -g_z_0_z_xyyz[k] - g_yz_0_z_xyyz[k] * ab_y + g_yz_0_z_xyyyz[k];

                g_yz_0_yz_xyzz[k] = -g_z_0_z_xyzz[k] - g_yz_0_z_xyzz[k] * ab_y + g_yz_0_z_xyyzz[k];

                g_yz_0_yz_xzzz[k] = -g_z_0_z_xzzz[k] - g_yz_0_z_xzzz[k] * ab_y + g_yz_0_z_xyzzz[k];

                g_yz_0_yz_yyyy[k] = -g_z_0_z_yyyy[k] - g_yz_0_z_yyyy[k] * ab_y + g_yz_0_z_yyyyy[k];

                g_yz_0_yz_yyyz[k] = -g_z_0_z_yyyz[k] - g_yz_0_z_yyyz[k] * ab_y + g_yz_0_z_yyyyz[k];

                g_yz_0_yz_yyzz[k] = -g_z_0_z_yyzz[k] - g_yz_0_z_yyzz[k] * ab_y + g_yz_0_z_yyyzz[k];

                g_yz_0_yz_yzzz[k] = -g_z_0_z_yzzz[k] - g_yz_0_z_yzzz[k] * ab_y + g_yz_0_z_yyzzz[k];

                g_yz_0_yz_zzzz[k] = -g_z_0_z_zzzz[k] - g_yz_0_z_zzzz[k] * ab_y + g_yz_0_z_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxzz, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyzz, g_y_0_z_xzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz, g_yz_0_z_xxxx, g_yz_0_z_xxxxz, g_yz_0_z_xxxy, g_yz_0_z_xxxyz, g_yz_0_z_xxxz, g_yz_0_z_xxxzz, g_yz_0_z_xxyy, g_yz_0_z_xxyyz, g_yz_0_z_xxyz, g_yz_0_z_xxyzz, g_yz_0_z_xxzz, g_yz_0_z_xxzzz, g_yz_0_z_xyyy, g_yz_0_z_xyyyz, g_yz_0_z_xyyz, g_yz_0_z_xyyzz, g_yz_0_z_xyzz, g_yz_0_z_xyzzz, g_yz_0_z_xzzz, g_yz_0_z_xzzzz, g_yz_0_z_yyyy, g_yz_0_z_yyyyz, g_yz_0_z_yyyz, g_yz_0_z_yyyzz, g_yz_0_z_yyzz, g_yz_0_z_yyzzz, g_yz_0_z_yzzz, g_yz_0_z_yzzzz, g_yz_0_z_zzzz, g_yz_0_z_zzzzz, g_yz_0_zz_xxxx, g_yz_0_zz_xxxy, g_yz_0_zz_xxxz, g_yz_0_zz_xxyy, g_yz_0_zz_xxyz, g_yz_0_zz_xxzz, g_yz_0_zz_xyyy, g_yz_0_zz_xyyz, g_yz_0_zz_xyzz, g_yz_0_zz_xzzz, g_yz_0_zz_yyyy, g_yz_0_zz_yyyz, g_yz_0_zz_yyzz, g_yz_0_zz_yzzz, g_yz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zz_xxxx[k] = -g_y_0_z_xxxx[k] - g_yz_0_z_xxxx[k] * ab_z + g_yz_0_z_xxxxz[k];

                g_yz_0_zz_xxxy[k] = -g_y_0_z_xxxy[k] - g_yz_0_z_xxxy[k] * ab_z + g_yz_0_z_xxxyz[k];

                g_yz_0_zz_xxxz[k] = -g_y_0_z_xxxz[k] - g_yz_0_z_xxxz[k] * ab_z + g_yz_0_z_xxxzz[k];

                g_yz_0_zz_xxyy[k] = -g_y_0_z_xxyy[k] - g_yz_0_z_xxyy[k] * ab_z + g_yz_0_z_xxyyz[k];

                g_yz_0_zz_xxyz[k] = -g_y_0_z_xxyz[k] - g_yz_0_z_xxyz[k] * ab_z + g_yz_0_z_xxyzz[k];

                g_yz_0_zz_xxzz[k] = -g_y_0_z_xxzz[k] - g_yz_0_z_xxzz[k] * ab_z + g_yz_0_z_xxzzz[k];

                g_yz_0_zz_xyyy[k] = -g_y_0_z_xyyy[k] - g_yz_0_z_xyyy[k] * ab_z + g_yz_0_z_xyyyz[k];

                g_yz_0_zz_xyyz[k] = -g_y_0_z_xyyz[k] - g_yz_0_z_xyyz[k] * ab_z + g_yz_0_z_xyyzz[k];

                g_yz_0_zz_xyzz[k] = -g_y_0_z_xyzz[k] - g_yz_0_z_xyzz[k] * ab_z + g_yz_0_z_xyzzz[k];

                g_yz_0_zz_xzzz[k] = -g_y_0_z_xzzz[k] - g_yz_0_z_xzzz[k] * ab_z + g_yz_0_z_xzzzz[k];

                g_yz_0_zz_yyyy[k] = -g_y_0_z_yyyy[k] - g_yz_0_z_yyyy[k] * ab_z + g_yz_0_z_yyyyz[k];

                g_yz_0_zz_yyyz[k] = -g_y_0_z_yyyz[k] - g_yz_0_z_yyyz[k] * ab_z + g_yz_0_z_yyyzz[k];

                g_yz_0_zz_yyzz[k] = -g_y_0_z_yyzz[k] - g_yz_0_z_yyzz[k] * ab_z + g_yz_0_z_yyzzz[k];

                g_yz_0_zz_yzzz[k] = -g_y_0_z_yzzz[k] - g_yz_0_z_yzzz[k] * ab_z + g_yz_0_z_yzzzz[k];

                g_yz_0_zz_zzzz[k] = -g_y_0_z_zzzz[k] - g_yz_0_z_zzzz[k] * ab_z + g_yz_0_z_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 455 * ccomps * dcomps);

            auto g_zz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 461 * ccomps * dcomps);

            auto g_zz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_x_xxxx, g_zz_0_x_xxxxx, g_zz_0_x_xxxxy, g_zz_0_x_xxxxz, g_zz_0_x_xxxy, g_zz_0_x_xxxyy, g_zz_0_x_xxxyz, g_zz_0_x_xxxz, g_zz_0_x_xxxzz, g_zz_0_x_xxyy, g_zz_0_x_xxyyy, g_zz_0_x_xxyyz, g_zz_0_x_xxyz, g_zz_0_x_xxyzz, g_zz_0_x_xxzz, g_zz_0_x_xxzzz, g_zz_0_x_xyyy, g_zz_0_x_xyyyy, g_zz_0_x_xyyyz, g_zz_0_x_xyyz, g_zz_0_x_xyyzz, g_zz_0_x_xyzz, g_zz_0_x_xyzzz, g_zz_0_x_xzzz, g_zz_0_x_xzzzz, g_zz_0_x_yyyy, g_zz_0_x_yyyz, g_zz_0_x_yyzz, g_zz_0_x_yzzz, g_zz_0_x_zzzz, g_zz_0_xx_xxxx, g_zz_0_xx_xxxy, g_zz_0_xx_xxxz, g_zz_0_xx_xxyy, g_zz_0_xx_xxyz, g_zz_0_xx_xxzz, g_zz_0_xx_xyyy, g_zz_0_xx_xyyz, g_zz_0_xx_xyzz, g_zz_0_xx_xzzz, g_zz_0_xx_yyyy, g_zz_0_xx_yyyz, g_zz_0_xx_yyzz, g_zz_0_xx_yzzz, g_zz_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xx_xxxx[k] = -g_zz_0_x_xxxx[k] * ab_x + g_zz_0_x_xxxxx[k];

                g_zz_0_xx_xxxy[k] = -g_zz_0_x_xxxy[k] * ab_x + g_zz_0_x_xxxxy[k];

                g_zz_0_xx_xxxz[k] = -g_zz_0_x_xxxz[k] * ab_x + g_zz_0_x_xxxxz[k];

                g_zz_0_xx_xxyy[k] = -g_zz_0_x_xxyy[k] * ab_x + g_zz_0_x_xxxyy[k];

                g_zz_0_xx_xxyz[k] = -g_zz_0_x_xxyz[k] * ab_x + g_zz_0_x_xxxyz[k];

                g_zz_0_xx_xxzz[k] = -g_zz_0_x_xxzz[k] * ab_x + g_zz_0_x_xxxzz[k];

                g_zz_0_xx_xyyy[k] = -g_zz_0_x_xyyy[k] * ab_x + g_zz_0_x_xxyyy[k];

                g_zz_0_xx_xyyz[k] = -g_zz_0_x_xyyz[k] * ab_x + g_zz_0_x_xxyyz[k];

                g_zz_0_xx_xyzz[k] = -g_zz_0_x_xyzz[k] * ab_x + g_zz_0_x_xxyzz[k];

                g_zz_0_xx_xzzz[k] = -g_zz_0_x_xzzz[k] * ab_x + g_zz_0_x_xxzzz[k];

                g_zz_0_xx_yyyy[k] = -g_zz_0_x_yyyy[k] * ab_x + g_zz_0_x_xyyyy[k];

                g_zz_0_xx_yyyz[k] = -g_zz_0_x_yyyz[k] * ab_x + g_zz_0_x_xyyyz[k];

                g_zz_0_xx_yyzz[k] = -g_zz_0_x_yyzz[k] * ab_x + g_zz_0_x_xyyzz[k];

                g_zz_0_xx_yzzz[k] = -g_zz_0_x_yzzz[k] * ab_x + g_zz_0_x_xyzzz[k];

                g_zz_0_xx_zzzz[k] = -g_zz_0_x_zzzz[k] * ab_x + g_zz_0_x_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 467 * ccomps * dcomps);

            auto g_zz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 473 * ccomps * dcomps);

            auto g_zz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xy_xxxx, g_zz_0_xy_xxxy, g_zz_0_xy_xxxz, g_zz_0_xy_xxyy, g_zz_0_xy_xxyz, g_zz_0_xy_xxzz, g_zz_0_xy_xyyy, g_zz_0_xy_xyyz, g_zz_0_xy_xyzz, g_zz_0_xy_xzzz, g_zz_0_xy_yyyy, g_zz_0_xy_yyyz, g_zz_0_xy_yyzz, g_zz_0_xy_yzzz, g_zz_0_xy_zzzz, g_zz_0_y_xxxx, g_zz_0_y_xxxxx, g_zz_0_y_xxxxy, g_zz_0_y_xxxxz, g_zz_0_y_xxxy, g_zz_0_y_xxxyy, g_zz_0_y_xxxyz, g_zz_0_y_xxxz, g_zz_0_y_xxxzz, g_zz_0_y_xxyy, g_zz_0_y_xxyyy, g_zz_0_y_xxyyz, g_zz_0_y_xxyz, g_zz_0_y_xxyzz, g_zz_0_y_xxzz, g_zz_0_y_xxzzz, g_zz_0_y_xyyy, g_zz_0_y_xyyyy, g_zz_0_y_xyyyz, g_zz_0_y_xyyz, g_zz_0_y_xyyzz, g_zz_0_y_xyzz, g_zz_0_y_xyzzz, g_zz_0_y_xzzz, g_zz_0_y_xzzzz, g_zz_0_y_yyyy, g_zz_0_y_yyyz, g_zz_0_y_yyzz, g_zz_0_y_yzzz, g_zz_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xy_xxxx[k] = -g_zz_0_y_xxxx[k] * ab_x + g_zz_0_y_xxxxx[k];

                g_zz_0_xy_xxxy[k] = -g_zz_0_y_xxxy[k] * ab_x + g_zz_0_y_xxxxy[k];

                g_zz_0_xy_xxxz[k] = -g_zz_0_y_xxxz[k] * ab_x + g_zz_0_y_xxxxz[k];

                g_zz_0_xy_xxyy[k] = -g_zz_0_y_xxyy[k] * ab_x + g_zz_0_y_xxxyy[k];

                g_zz_0_xy_xxyz[k] = -g_zz_0_y_xxyz[k] * ab_x + g_zz_0_y_xxxyz[k];

                g_zz_0_xy_xxzz[k] = -g_zz_0_y_xxzz[k] * ab_x + g_zz_0_y_xxxzz[k];

                g_zz_0_xy_xyyy[k] = -g_zz_0_y_xyyy[k] * ab_x + g_zz_0_y_xxyyy[k];

                g_zz_0_xy_xyyz[k] = -g_zz_0_y_xyyz[k] * ab_x + g_zz_0_y_xxyyz[k];

                g_zz_0_xy_xyzz[k] = -g_zz_0_y_xyzz[k] * ab_x + g_zz_0_y_xxyzz[k];

                g_zz_0_xy_xzzz[k] = -g_zz_0_y_xzzz[k] * ab_x + g_zz_0_y_xxzzz[k];

                g_zz_0_xy_yyyy[k] = -g_zz_0_y_yyyy[k] * ab_x + g_zz_0_y_xyyyy[k];

                g_zz_0_xy_yyyz[k] = -g_zz_0_y_yyyz[k] * ab_x + g_zz_0_y_xyyyz[k];

                g_zz_0_xy_yyzz[k] = -g_zz_0_y_yyzz[k] * ab_x + g_zz_0_y_xyyzz[k];

                g_zz_0_xy_yzzz[k] = -g_zz_0_y_yzzz[k] * ab_x + g_zz_0_y_xyzzz[k];

                g_zz_0_xy_zzzz[k] = -g_zz_0_y_zzzz[k] * ab_x + g_zz_0_y_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 485 * ccomps * dcomps);

            auto g_zz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 491 * ccomps * dcomps);

            auto g_zz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xz_xxxx, g_zz_0_xz_xxxy, g_zz_0_xz_xxxz, g_zz_0_xz_xxyy, g_zz_0_xz_xxyz, g_zz_0_xz_xxzz, g_zz_0_xz_xyyy, g_zz_0_xz_xyyz, g_zz_0_xz_xyzz, g_zz_0_xz_xzzz, g_zz_0_xz_yyyy, g_zz_0_xz_yyyz, g_zz_0_xz_yyzz, g_zz_0_xz_yzzz, g_zz_0_xz_zzzz, g_zz_0_z_xxxx, g_zz_0_z_xxxxx, g_zz_0_z_xxxxy, g_zz_0_z_xxxxz, g_zz_0_z_xxxy, g_zz_0_z_xxxyy, g_zz_0_z_xxxyz, g_zz_0_z_xxxz, g_zz_0_z_xxxzz, g_zz_0_z_xxyy, g_zz_0_z_xxyyy, g_zz_0_z_xxyyz, g_zz_0_z_xxyz, g_zz_0_z_xxyzz, g_zz_0_z_xxzz, g_zz_0_z_xxzzz, g_zz_0_z_xyyy, g_zz_0_z_xyyyy, g_zz_0_z_xyyyz, g_zz_0_z_xyyz, g_zz_0_z_xyyzz, g_zz_0_z_xyzz, g_zz_0_z_xyzzz, g_zz_0_z_xzzz, g_zz_0_z_xzzzz, g_zz_0_z_yyyy, g_zz_0_z_yyyz, g_zz_0_z_yyzz, g_zz_0_z_yzzz, g_zz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xz_xxxx[k] = -g_zz_0_z_xxxx[k] * ab_x + g_zz_0_z_xxxxx[k];

                g_zz_0_xz_xxxy[k] = -g_zz_0_z_xxxy[k] * ab_x + g_zz_0_z_xxxxy[k];

                g_zz_0_xz_xxxz[k] = -g_zz_0_z_xxxz[k] * ab_x + g_zz_0_z_xxxxz[k];

                g_zz_0_xz_xxyy[k] = -g_zz_0_z_xxyy[k] * ab_x + g_zz_0_z_xxxyy[k];

                g_zz_0_xz_xxyz[k] = -g_zz_0_z_xxyz[k] * ab_x + g_zz_0_z_xxxyz[k];

                g_zz_0_xz_xxzz[k] = -g_zz_0_z_xxzz[k] * ab_x + g_zz_0_z_xxxzz[k];

                g_zz_0_xz_xyyy[k] = -g_zz_0_z_xyyy[k] * ab_x + g_zz_0_z_xxyyy[k];

                g_zz_0_xz_xyyz[k] = -g_zz_0_z_xyyz[k] * ab_x + g_zz_0_z_xxyyz[k];

                g_zz_0_xz_xyzz[k] = -g_zz_0_z_xyzz[k] * ab_x + g_zz_0_z_xxyzz[k];

                g_zz_0_xz_xzzz[k] = -g_zz_0_z_xzzz[k] * ab_x + g_zz_0_z_xxzzz[k];

                g_zz_0_xz_yyyy[k] = -g_zz_0_z_yyyy[k] * ab_x + g_zz_0_z_xyyyy[k];

                g_zz_0_xz_yyyz[k] = -g_zz_0_z_yyyz[k] * ab_x + g_zz_0_z_xyyyz[k];

                g_zz_0_xz_yyzz[k] = -g_zz_0_z_yyzz[k] * ab_x + g_zz_0_z_xyyzz[k];

                g_zz_0_xz_yzzz[k] = -g_zz_0_z_yzzz[k] * ab_x + g_zz_0_z_xyzzz[k];

                g_zz_0_xz_zzzz[k] = -g_zz_0_z_zzzz[k] * ab_x + g_zz_0_z_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 497 * ccomps * dcomps);

            auto g_zz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_y_xxxx, g_zz_0_y_xxxxy, g_zz_0_y_xxxy, g_zz_0_y_xxxyy, g_zz_0_y_xxxyz, g_zz_0_y_xxxz, g_zz_0_y_xxyy, g_zz_0_y_xxyyy, g_zz_0_y_xxyyz, g_zz_0_y_xxyz, g_zz_0_y_xxyzz, g_zz_0_y_xxzz, g_zz_0_y_xyyy, g_zz_0_y_xyyyy, g_zz_0_y_xyyyz, g_zz_0_y_xyyz, g_zz_0_y_xyyzz, g_zz_0_y_xyzz, g_zz_0_y_xyzzz, g_zz_0_y_xzzz, g_zz_0_y_yyyy, g_zz_0_y_yyyyy, g_zz_0_y_yyyyz, g_zz_0_y_yyyz, g_zz_0_y_yyyzz, g_zz_0_y_yyzz, g_zz_0_y_yyzzz, g_zz_0_y_yzzz, g_zz_0_y_yzzzz, g_zz_0_y_zzzz, g_zz_0_yy_xxxx, g_zz_0_yy_xxxy, g_zz_0_yy_xxxz, g_zz_0_yy_xxyy, g_zz_0_yy_xxyz, g_zz_0_yy_xxzz, g_zz_0_yy_xyyy, g_zz_0_yy_xyyz, g_zz_0_yy_xyzz, g_zz_0_yy_xzzz, g_zz_0_yy_yyyy, g_zz_0_yy_yyyz, g_zz_0_yy_yyzz, g_zz_0_yy_yzzz, g_zz_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yy_xxxx[k] = -g_zz_0_y_xxxx[k] * ab_y + g_zz_0_y_xxxxy[k];

                g_zz_0_yy_xxxy[k] = -g_zz_0_y_xxxy[k] * ab_y + g_zz_0_y_xxxyy[k];

                g_zz_0_yy_xxxz[k] = -g_zz_0_y_xxxz[k] * ab_y + g_zz_0_y_xxxyz[k];

                g_zz_0_yy_xxyy[k] = -g_zz_0_y_xxyy[k] * ab_y + g_zz_0_y_xxyyy[k];

                g_zz_0_yy_xxyz[k] = -g_zz_0_y_xxyz[k] * ab_y + g_zz_0_y_xxyyz[k];

                g_zz_0_yy_xxzz[k] = -g_zz_0_y_xxzz[k] * ab_y + g_zz_0_y_xxyzz[k];

                g_zz_0_yy_xyyy[k] = -g_zz_0_y_xyyy[k] * ab_y + g_zz_0_y_xyyyy[k];

                g_zz_0_yy_xyyz[k] = -g_zz_0_y_xyyz[k] * ab_y + g_zz_0_y_xyyyz[k];

                g_zz_0_yy_xyzz[k] = -g_zz_0_y_xyzz[k] * ab_y + g_zz_0_y_xyyzz[k];

                g_zz_0_yy_xzzz[k] = -g_zz_0_y_xzzz[k] * ab_y + g_zz_0_y_xyzzz[k];

                g_zz_0_yy_yyyy[k] = -g_zz_0_y_yyyy[k] * ab_y + g_zz_0_y_yyyyy[k];

                g_zz_0_yy_yyyz[k] = -g_zz_0_y_yyyz[k] * ab_y + g_zz_0_y_yyyyz[k];

                g_zz_0_yy_yyzz[k] = -g_zz_0_y_yyzz[k] * ab_y + g_zz_0_y_yyyzz[k];

                g_zz_0_yy_yzzz[k] = -g_zz_0_y_yzzz[k] * ab_y + g_zz_0_y_yyzzz[k];

                g_zz_0_yy_zzzz[k] = -g_zz_0_y_zzzz[k] * ab_y + g_zz_0_y_yzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yz_xxxx, g_zz_0_yz_xxxy, g_zz_0_yz_xxxz, g_zz_0_yz_xxyy, g_zz_0_yz_xxyz, g_zz_0_yz_xxzz, g_zz_0_yz_xyyy, g_zz_0_yz_xyyz, g_zz_0_yz_xyzz, g_zz_0_yz_xzzz, g_zz_0_yz_yyyy, g_zz_0_yz_yyyz, g_zz_0_yz_yyzz, g_zz_0_yz_yzzz, g_zz_0_yz_zzzz, g_zz_0_z_xxxx, g_zz_0_z_xxxxy, g_zz_0_z_xxxy, g_zz_0_z_xxxyy, g_zz_0_z_xxxyz, g_zz_0_z_xxxz, g_zz_0_z_xxyy, g_zz_0_z_xxyyy, g_zz_0_z_xxyyz, g_zz_0_z_xxyz, g_zz_0_z_xxyzz, g_zz_0_z_xxzz, g_zz_0_z_xyyy, g_zz_0_z_xyyyy, g_zz_0_z_xyyyz, g_zz_0_z_xyyz, g_zz_0_z_xyyzz, g_zz_0_z_xyzz, g_zz_0_z_xyzzz, g_zz_0_z_xzzz, g_zz_0_z_yyyy, g_zz_0_z_yyyyy, g_zz_0_z_yyyyz, g_zz_0_z_yyyz, g_zz_0_z_yyyzz, g_zz_0_z_yyzz, g_zz_0_z_yyzzz, g_zz_0_z_yzzz, g_zz_0_z_yzzzz, g_zz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yz_xxxx[k] = -g_zz_0_z_xxxx[k] * ab_y + g_zz_0_z_xxxxy[k];

                g_zz_0_yz_xxxy[k] = -g_zz_0_z_xxxy[k] * ab_y + g_zz_0_z_xxxyy[k];

                g_zz_0_yz_xxxz[k] = -g_zz_0_z_xxxz[k] * ab_y + g_zz_0_z_xxxyz[k];

                g_zz_0_yz_xxyy[k] = -g_zz_0_z_xxyy[k] * ab_y + g_zz_0_z_xxyyy[k];

                g_zz_0_yz_xxyz[k] = -g_zz_0_z_xxyz[k] * ab_y + g_zz_0_z_xxyyz[k];

                g_zz_0_yz_xxzz[k] = -g_zz_0_z_xxzz[k] * ab_y + g_zz_0_z_xxyzz[k];

                g_zz_0_yz_xyyy[k] = -g_zz_0_z_xyyy[k] * ab_y + g_zz_0_z_xyyyy[k];

                g_zz_0_yz_xyyz[k] = -g_zz_0_z_xyyz[k] * ab_y + g_zz_0_z_xyyyz[k];

                g_zz_0_yz_xyzz[k] = -g_zz_0_z_xyzz[k] * ab_y + g_zz_0_z_xyyzz[k];

                g_zz_0_yz_xzzz[k] = -g_zz_0_z_xzzz[k] * ab_y + g_zz_0_z_xyzzz[k];

                g_zz_0_yz_yyyy[k] = -g_zz_0_z_yyyy[k] * ab_y + g_zz_0_z_yyyyy[k];

                g_zz_0_yz_yyyz[k] = -g_zz_0_z_yyyz[k] * ab_y + g_zz_0_z_yyyyz[k];

                g_zz_0_yz_yyzz[k] = -g_zz_0_z_yyzz[k] * ab_y + g_zz_0_z_yyyzz[k];

                g_zz_0_yz_yzzz[k] = -g_zz_0_z_yzzz[k] * ab_y + g_zz_0_z_yyzzz[k];

                g_zz_0_yz_zzzz[k] = -g_zz_0_z_zzzz[k] * ab_y + g_zz_0_z_yzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz, g_zz_0_z_xxxx, g_zz_0_z_xxxxz, g_zz_0_z_xxxy, g_zz_0_z_xxxyz, g_zz_0_z_xxxz, g_zz_0_z_xxxzz, g_zz_0_z_xxyy, g_zz_0_z_xxyyz, g_zz_0_z_xxyz, g_zz_0_z_xxyzz, g_zz_0_z_xxzz, g_zz_0_z_xxzzz, g_zz_0_z_xyyy, g_zz_0_z_xyyyz, g_zz_0_z_xyyz, g_zz_0_z_xyyzz, g_zz_0_z_xyzz, g_zz_0_z_xyzzz, g_zz_0_z_xzzz, g_zz_0_z_xzzzz, g_zz_0_z_yyyy, g_zz_0_z_yyyyz, g_zz_0_z_yyyz, g_zz_0_z_yyyzz, g_zz_0_z_yyzz, g_zz_0_z_yyzzz, g_zz_0_z_yzzz, g_zz_0_z_yzzzz, g_zz_0_z_zzzz, g_zz_0_z_zzzzz, g_zz_0_zz_xxxx, g_zz_0_zz_xxxy, g_zz_0_zz_xxxz, g_zz_0_zz_xxyy, g_zz_0_zz_xxyz, g_zz_0_zz_xxzz, g_zz_0_zz_xyyy, g_zz_0_zz_xyyz, g_zz_0_zz_xyzz, g_zz_0_zz_xzzz, g_zz_0_zz_yyyy, g_zz_0_zz_yyyz, g_zz_0_zz_yyzz, g_zz_0_zz_yzzz, g_zz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zz_xxxx[k] = -2.0 * g_z_0_z_xxxx[k] - g_zz_0_z_xxxx[k] * ab_z + g_zz_0_z_xxxxz[k];

                g_zz_0_zz_xxxy[k] = -2.0 * g_z_0_z_xxxy[k] - g_zz_0_z_xxxy[k] * ab_z + g_zz_0_z_xxxyz[k];

                g_zz_0_zz_xxxz[k] = -2.0 * g_z_0_z_xxxz[k] - g_zz_0_z_xxxz[k] * ab_z + g_zz_0_z_xxxzz[k];

                g_zz_0_zz_xxyy[k] = -2.0 * g_z_0_z_xxyy[k] - g_zz_0_z_xxyy[k] * ab_z + g_zz_0_z_xxyyz[k];

                g_zz_0_zz_xxyz[k] = -2.0 * g_z_0_z_xxyz[k] - g_zz_0_z_xxyz[k] * ab_z + g_zz_0_z_xxyzz[k];

                g_zz_0_zz_xxzz[k] = -2.0 * g_z_0_z_xxzz[k] - g_zz_0_z_xxzz[k] * ab_z + g_zz_0_z_xxzzz[k];

                g_zz_0_zz_xyyy[k] = -2.0 * g_z_0_z_xyyy[k] - g_zz_0_z_xyyy[k] * ab_z + g_zz_0_z_xyyyz[k];

                g_zz_0_zz_xyyz[k] = -2.0 * g_z_0_z_xyyz[k] - g_zz_0_z_xyyz[k] * ab_z + g_zz_0_z_xyyzz[k];

                g_zz_0_zz_xyzz[k] = -2.0 * g_z_0_z_xyzz[k] - g_zz_0_z_xyzz[k] * ab_z + g_zz_0_z_xyzzz[k];

                g_zz_0_zz_xzzz[k] = -2.0 * g_z_0_z_xzzz[k] - g_zz_0_z_xzzz[k] * ab_z + g_zz_0_z_xzzzz[k];

                g_zz_0_zz_yyyy[k] = -2.0 * g_z_0_z_yyyy[k] - g_zz_0_z_yyyy[k] * ab_z + g_zz_0_z_yyyyz[k];

                g_zz_0_zz_yyyz[k] = -2.0 * g_z_0_z_yyyz[k] - g_zz_0_z_yyyz[k] * ab_z + g_zz_0_z_yyyzz[k];

                g_zz_0_zz_yyzz[k] = -2.0 * g_z_0_z_yyzz[k] - g_zz_0_z_yyzz[k] * ab_z + g_zz_0_z_yyzzz[k];

                g_zz_0_zz_yzzz[k] = -2.0 * g_z_0_z_yzzz[k] - g_zz_0_z_yzzz[k] * ab_z + g_zz_0_z_yzzzz[k];

                g_zz_0_zz_zzzz[k] = -2.0 * g_z_0_z_zzzz[k] - g_zz_0_z_zzzz[k] * ab_z + g_zz_0_z_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

