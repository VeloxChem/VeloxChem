#include "ElectronRepulsionGeom1100ContrRecDGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_dgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_dgxx,
                                            const size_t idx_geom_01_pgxx,
                                            const size_t idx_geom_10_pgxx,
                                            const size_t idx_geom_11_pgxx,
                                            const size_t idx_geom_11_phxx,
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

            const auto pg_geom_01_off = idx_geom_01_pgxx + i * dcomps + j;

            auto g_0_x_x_xxxx = cbuffer.data(pg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxy = cbuffer.data(pg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxz = cbuffer.data(pg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxyy = cbuffer.data(pg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxyz = cbuffer.data(pg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxzz = cbuffer.data(pg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xyyy = cbuffer.data(pg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xyyz = cbuffer.data(pg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xyzz = cbuffer.data(pg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xzzz = cbuffer.data(pg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_yyyy = cbuffer.data(pg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_yyyz = cbuffer.data(pg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_yyzz = cbuffer.data(pg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_yzzz = cbuffer.data(pg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_zzzz = cbuffer.data(pg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_y_xxxx = cbuffer.data(pg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_y_xxxy = cbuffer.data(pg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_y_xxxz = cbuffer.data(pg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_y_xxyy = cbuffer.data(pg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_y_xxyz = cbuffer.data(pg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_y_xxzz = cbuffer.data(pg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_y_xyyy = cbuffer.data(pg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_y_xyyz = cbuffer.data(pg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_y_xyzz = cbuffer.data(pg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_y_xzzz = cbuffer.data(pg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_y_yyyy = cbuffer.data(pg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_y_yyyz = cbuffer.data(pg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_y_yyzz = cbuffer.data(pg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_y_yzzz = cbuffer.data(pg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_y_zzzz = cbuffer.data(pg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_z_xxxx = cbuffer.data(pg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_z_xxxy = cbuffer.data(pg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_z_xxxz = cbuffer.data(pg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_z_xxyy = cbuffer.data(pg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_z_xxyz = cbuffer.data(pg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_z_xxzz = cbuffer.data(pg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_z_xyyy = cbuffer.data(pg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_z_xyyz = cbuffer.data(pg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_z_xyzz = cbuffer.data(pg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_z_xzzz = cbuffer.data(pg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_z_yyyy = cbuffer.data(pg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_z_yyyz = cbuffer.data(pg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_z_yyzz = cbuffer.data(pg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_z_yzzz = cbuffer.data(pg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_z_zzzz = cbuffer.data(pg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_x_xxxx = cbuffer.data(pg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_x_xxxy = cbuffer.data(pg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_x_xxxz = cbuffer.data(pg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_x_xxyy = cbuffer.data(pg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_x_xxyz = cbuffer.data(pg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_x_xxzz = cbuffer.data(pg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_x_xyyy = cbuffer.data(pg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_x_xyyz = cbuffer.data(pg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_x_xyzz = cbuffer.data(pg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_x_xzzz = cbuffer.data(pg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_x_yyyy = cbuffer.data(pg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_x_yyyz = cbuffer.data(pg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_x_yyzz = cbuffer.data(pg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_x_yzzz = cbuffer.data(pg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_x_zzzz = cbuffer.data(pg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_y_xxxx = cbuffer.data(pg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_y_xxxy = cbuffer.data(pg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_y_xxxz = cbuffer.data(pg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_y_xxyy = cbuffer.data(pg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_y_xxyz = cbuffer.data(pg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_y_xxzz = cbuffer.data(pg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_y_xyyy = cbuffer.data(pg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_y_xyyz = cbuffer.data(pg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_y_xyzz = cbuffer.data(pg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_y_xzzz = cbuffer.data(pg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_y_yyyy = cbuffer.data(pg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_y_yyyz = cbuffer.data(pg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_y_yyzz = cbuffer.data(pg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_y_yzzz = cbuffer.data(pg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_y_zzzz = cbuffer.data(pg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_z_xxxx = cbuffer.data(pg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_z_xxxy = cbuffer.data(pg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_z_xxxz = cbuffer.data(pg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_z_xxyy = cbuffer.data(pg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_z_xxyz = cbuffer.data(pg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_z_xxzz = cbuffer.data(pg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_z_xyyy = cbuffer.data(pg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_z_xyyz = cbuffer.data(pg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_z_xyzz = cbuffer.data(pg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_z_xzzz = cbuffer.data(pg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_z_yyyy = cbuffer.data(pg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_z_yyyz = cbuffer.data(pg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_z_yyzz = cbuffer.data(pg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_z_yzzz = cbuffer.data(pg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_z_zzzz = cbuffer.data(pg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_z_x_xxxx = cbuffer.data(pg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_z_x_xxxy = cbuffer.data(pg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_z_x_xxxz = cbuffer.data(pg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_z_x_xxyy = cbuffer.data(pg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_z_x_xxyz = cbuffer.data(pg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_z_x_xxzz = cbuffer.data(pg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_z_x_xyyy = cbuffer.data(pg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_z_x_xyyz = cbuffer.data(pg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_z_x_xyzz = cbuffer.data(pg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_z_x_xzzz = cbuffer.data(pg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_z_x_yyyy = cbuffer.data(pg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_z_x_yyyz = cbuffer.data(pg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_z_x_yyzz = cbuffer.data(pg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_z_x_yzzz = cbuffer.data(pg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_z_x_zzzz = cbuffer.data(pg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_z_y_xxxx = cbuffer.data(pg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_z_y_xxxy = cbuffer.data(pg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_z_y_xxxz = cbuffer.data(pg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_z_y_xxyy = cbuffer.data(pg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_z_y_xxyz = cbuffer.data(pg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_z_y_xxzz = cbuffer.data(pg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_z_y_xyyy = cbuffer.data(pg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_z_y_xyyz = cbuffer.data(pg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_z_y_xyzz = cbuffer.data(pg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_z_y_xzzz = cbuffer.data(pg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_z_y_yyyy = cbuffer.data(pg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_z_y_yyyz = cbuffer.data(pg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_z_y_yyzz = cbuffer.data(pg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_z_y_yzzz = cbuffer.data(pg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_z_y_zzzz = cbuffer.data(pg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_z_z_xxxx = cbuffer.data(pg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_z_xxxy = cbuffer.data(pg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_z_xxxz = cbuffer.data(pg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_z_z_xxyy = cbuffer.data(pg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_z_xxyz = cbuffer.data(pg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_z_xxzz = cbuffer.data(pg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_z_xyyy = cbuffer.data(pg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_z_xyyz = cbuffer.data(pg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_z_xyzz = cbuffer.data(pg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_z_xzzz = cbuffer.data(pg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_z_yyyy = cbuffer.data(pg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_z_yyyz = cbuffer.data(pg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_z_yyzz = cbuffer.data(pg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_z_yzzz = cbuffer.data(pg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_z_zzzz = cbuffer.data(pg_geom_01_off + 134 * ccomps * dcomps);

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

            const auto pg_geom_11_off = idx_geom_11_pgxx + i * dcomps + j;

            auto g_x_x_x_xxxx = cbuffer.data(pg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxy = cbuffer.data(pg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxz = cbuffer.data(pg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxyy = cbuffer.data(pg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxyz = cbuffer.data(pg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxzz = cbuffer.data(pg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xyyy = cbuffer.data(pg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xyyz = cbuffer.data(pg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xyzz = cbuffer.data(pg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xzzz = cbuffer.data(pg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_yyyy = cbuffer.data(pg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_yyyz = cbuffer.data(pg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_yyzz = cbuffer.data(pg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_yzzz = cbuffer.data(pg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_zzzz = cbuffer.data(pg_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_y_xxxx = cbuffer.data(pg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_y_xxxy = cbuffer.data(pg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_y_xxxz = cbuffer.data(pg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_y_xxyy = cbuffer.data(pg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_y_xxyz = cbuffer.data(pg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_y_xxzz = cbuffer.data(pg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_y_xyyy = cbuffer.data(pg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_y_xyyz = cbuffer.data(pg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_y_xyzz = cbuffer.data(pg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_y_xzzz = cbuffer.data(pg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_y_yyyy = cbuffer.data(pg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_y_yyyz = cbuffer.data(pg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_y_yyzz = cbuffer.data(pg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_y_yzzz = cbuffer.data(pg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_zzzz = cbuffer.data(pg_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_z_xxxx = cbuffer.data(pg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_z_xxxy = cbuffer.data(pg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_z_xxxz = cbuffer.data(pg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_z_xxyy = cbuffer.data(pg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_z_xxyz = cbuffer.data(pg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_z_xxzz = cbuffer.data(pg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_z_xyyy = cbuffer.data(pg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_z_xyyz = cbuffer.data(pg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_z_xyzz = cbuffer.data(pg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_z_xzzz = cbuffer.data(pg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_z_yyyy = cbuffer.data(pg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_z_yyyz = cbuffer.data(pg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_z_yyzz = cbuffer.data(pg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_z_yzzz = cbuffer.data(pg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_z_zzzz = cbuffer.data(pg_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_x_xxxx = cbuffer.data(pg_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_x_xxxy = cbuffer.data(pg_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_x_xxxz = cbuffer.data(pg_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_x_xxyy = cbuffer.data(pg_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_x_xxyz = cbuffer.data(pg_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_x_xxzz = cbuffer.data(pg_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_x_xyyy = cbuffer.data(pg_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_x_xyyz = cbuffer.data(pg_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_x_xyzz = cbuffer.data(pg_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_x_xzzz = cbuffer.data(pg_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_x_yyyy = cbuffer.data(pg_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_x_yyyz = cbuffer.data(pg_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_x_yyzz = cbuffer.data(pg_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_x_yzzz = cbuffer.data(pg_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_x_zzzz = cbuffer.data(pg_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_y_xxxx = cbuffer.data(pg_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_y_xxxy = cbuffer.data(pg_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_y_xxxz = cbuffer.data(pg_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_y_xxyy = cbuffer.data(pg_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_y_xxyz = cbuffer.data(pg_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_y_xxzz = cbuffer.data(pg_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_y_xyyy = cbuffer.data(pg_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_y_xyyz = cbuffer.data(pg_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_y_xyzz = cbuffer.data(pg_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_y_xzzz = cbuffer.data(pg_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_y_yyyy = cbuffer.data(pg_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_y_yyyz = cbuffer.data(pg_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_y_yyzz = cbuffer.data(pg_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_y_yzzz = cbuffer.data(pg_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_y_zzzz = cbuffer.data(pg_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_z_xxxx = cbuffer.data(pg_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_z_xxxy = cbuffer.data(pg_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_z_xxxz = cbuffer.data(pg_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_z_xxyy = cbuffer.data(pg_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_z_xxyz = cbuffer.data(pg_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_z_xxzz = cbuffer.data(pg_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_z_xyyy = cbuffer.data(pg_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_z_xyyz = cbuffer.data(pg_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_z_xyzz = cbuffer.data(pg_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_z_xzzz = cbuffer.data(pg_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_z_yyyy = cbuffer.data(pg_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_z_yyyz = cbuffer.data(pg_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_z_yyzz = cbuffer.data(pg_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_z_yzzz = cbuffer.data(pg_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_z_zzzz = cbuffer.data(pg_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_z_x_xxxx = cbuffer.data(pg_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_x_xxxy = cbuffer.data(pg_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_x_xxxz = cbuffer.data(pg_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_x_xxyy = cbuffer.data(pg_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_x_xxyz = cbuffer.data(pg_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_x_xxzz = cbuffer.data(pg_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_z_x_xyyy = cbuffer.data(pg_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_x_xyyz = cbuffer.data(pg_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_x_xyzz = cbuffer.data(pg_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_x_xzzz = cbuffer.data(pg_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_x_yyyy = cbuffer.data(pg_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_x_yyyz = cbuffer.data(pg_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_z_x_yyzz = cbuffer.data(pg_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_x_yzzz = cbuffer.data(pg_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_x_zzzz = cbuffer.data(pg_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_y_xxxx = cbuffer.data(pg_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_y_xxxy = cbuffer.data(pg_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_y_xxxz = cbuffer.data(pg_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_z_y_xxyy = cbuffer.data(pg_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_z_y_xxyz = cbuffer.data(pg_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_z_y_xxzz = cbuffer.data(pg_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_z_y_xyyy = cbuffer.data(pg_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_z_y_xyyz = cbuffer.data(pg_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_z_y_xyzz = cbuffer.data(pg_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_z_y_xzzz = cbuffer.data(pg_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_z_y_yyyy = cbuffer.data(pg_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_z_y_yyyz = cbuffer.data(pg_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_z_y_yyzz = cbuffer.data(pg_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_z_y_yzzz = cbuffer.data(pg_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_z_y_zzzz = cbuffer.data(pg_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_z_z_xxxx = cbuffer.data(pg_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_z_xxxy = cbuffer.data(pg_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_z_xxxz = cbuffer.data(pg_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_z_xxyy = cbuffer.data(pg_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_z_xxyz = cbuffer.data(pg_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_z_xxzz = cbuffer.data(pg_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_z_xyyy = cbuffer.data(pg_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_z_xyyz = cbuffer.data(pg_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_z_xyzz = cbuffer.data(pg_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_z_xzzz = cbuffer.data(pg_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_z_yyyy = cbuffer.data(pg_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_z_yyyz = cbuffer.data(pg_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_z_yyzz = cbuffer.data(pg_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_z_yzzz = cbuffer.data(pg_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_z_zzzz = cbuffer.data(pg_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_x_xxxx = cbuffer.data(pg_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_x_xxxy = cbuffer.data(pg_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_x_xxxz = cbuffer.data(pg_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_x_x_xxyy = cbuffer.data(pg_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_x_xxyz = cbuffer.data(pg_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_x_xxzz = cbuffer.data(pg_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_x_xyyy = cbuffer.data(pg_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_x_xyyz = cbuffer.data(pg_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_x_xyzz = cbuffer.data(pg_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_x_x_xzzz = cbuffer.data(pg_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_x_x_yyyy = cbuffer.data(pg_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_x_x_yyyz = cbuffer.data(pg_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_x_x_yyzz = cbuffer.data(pg_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_x_x_yzzz = cbuffer.data(pg_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_x_x_zzzz = cbuffer.data(pg_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_x_y_xxxx = cbuffer.data(pg_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_x_y_xxxy = cbuffer.data(pg_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_x_y_xxxz = cbuffer.data(pg_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_x_y_xxyy = cbuffer.data(pg_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_x_y_xxyz = cbuffer.data(pg_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_x_y_xxzz = cbuffer.data(pg_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_x_y_xyyy = cbuffer.data(pg_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_x_y_xyyz = cbuffer.data(pg_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_x_y_xyzz = cbuffer.data(pg_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_x_y_xzzz = cbuffer.data(pg_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_x_y_yyyy = cbuffer.data(pg_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_x_y_yyyz = cbuffer.data(pg_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_x_y_yyzz = cbuffer.data(pg_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_x_y_yzzz = cbuffer.data(pg_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_x_y_zzzz = cbuffer.data(pg_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_x_z_xxxx = cbuffer.data(pg_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_x_z_xxxy = cbuffer.data(pg_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_x_z_xxxz = cbuffer.data(pg_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_x_z_xxyy = cbuffer.data(pg_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_x_z_xxyz = cbuffer.data(pg_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_x_z_xxzz = cbuffer.data(pg_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_x_z_xyyy = cbuffer.data(pg_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_x_z_xyyz = cbuffer.data(pg_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_x_z_xyzz = cbuffer.data(pg_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_x_z_xzzz = cbuffer.data(pg_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_x_z_yyyy = cbuffer.data(pg_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_x_z_yyyz = cbuffer.data(pg_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_x_z_yyzz = cbuffer.data(pg_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_x_z_yzzz = cbuffer.data(pg_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_x_z_zzzz = cbuffer.data(pg_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_y_x_xxxx = cbuffer.data(pg_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_y_x_xxxy = cbuffer.data(pg_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_y_x_xxxz = cbuffer.data(pg_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_y_x_xxyy = cbuffer.data(pg_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_y_x_xxyz = cbuffer.data(pg_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_y_x_xxzz = cbuffer.data(pg_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_y_x_xyyy = cbuffer.data(pg_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_y_x_xyyz = cbuffer.data(pg_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_y_x_xyzz = cbuffer.data(pg_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_y_x_xzzz = cbuffer.data(pg_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_y_x_yyyy = cbuffer.data(pg_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_y_x_yyyz = cbuffer.data(pg_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_y_x_yyzz = cbuffer.data(pg_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_y_x_yzzz = cbuffer.data(pg_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_y_x_zzzz = cbuffer.data(pg_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_y_y_xxxx = cbuffer.data(pg_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_y_y_xxxy = cbuffer.data(pg_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_y_y_xxxz = cbuffer.data(pg_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_y_y_xxyy = cbuffer.data(pg_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_y_y_xxyz = cbuffer.data(pg_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_y_y_xxzz = cbuffer.data(pg_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_y_y_xyyy = cbuffer.data(pg_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_y_y_xyyz = cbuffer.data(pg_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_y_y_xyzz = cbuffer.data(pg_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_y_y_xzzz = cbuffer.data(pg_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_y_y_yyyy = cbuffer.data(pg_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_y_y_yyyz = cbuffer.data(pg_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_y_y_yyzz = cbuffer.data(pg_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_y_y_yzzz = cbuffer.data(pg_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_y_y_zzzz = cbuffer.data(pg_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_y_z_xxxx = cbuffer.data(pg_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_y_z_xxxy = cbuffer.data(pg_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_y_z_xxxz = cbuffer.data(pg_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_y_z_xxyy = cbuffer.data(pg_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_y_z_xxyz = cbuffer.data(pg_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_y_z_xxzz = cbuffer.data(pg_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_y_z_xyyy = cbuffer.data(pg_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_y_z_xyyz = cbuffer.data(pg_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_y_z_xyzz = cbuffer.data(pg_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_y_z_xzzz = cbuffer.data(pg_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_y_z_yyyy = cbuffer.data(pg_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_y_z_yyyz = cbuffer.data(pg_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_y_z_yyzz = cbuffer.data(pg_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_y_z_yzzz = cbuffer.data(pg_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_y_z_zzzz = cbuffer.data(pg_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_z_x_xxxx = cbuffer.data(pg_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_z_x_xxxy = cbuffer.data(pg_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_z_x_xxxz = cbuffer.data(pg_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_z_x_xxyy = cbuffer.data(pg_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_z_x_xxyz = cbuffer.data(pg_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_z_x_xxzz = cbuffer.data(pg_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_z_x_xyyy = cbuffer.data(pg_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_z_x_xyyz = cbuffer.data(pg_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_z_x_xyzz = cbuffer.data(pg_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_z_x_xzzz = cbuffer.data(pg_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_z_x_yyyy = cbuffer.data(pg_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_z_x_yyyz = cbuffer.data(pg_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_z_x_yyzz = cbuffer.data(pg_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_z_x_yzzz = cbuffer.data(pg_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_z_x_zzzz = cbuffer.data(pg_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_z_y_xxxx = cbuffer.data(pg_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_z_y_xxxy = cbuffer.data(pg_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_z_y_xxxz = cbuffer.data(pg_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_z_y_xxyy = cbuffer.data(pg_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_z_y_xxyz = cbuffer.data(pg_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_z_y_xxzz = cbuffer.data(pg_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_z_y_xyyy = cbuffer.data(pg_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_z_y_xyyz = cbuffer.data(pg_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_z_y_xyzz = cbuffer.data(pg_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_z_y_xzzz = cbuffer.data(pg_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_z_y_yyyy = cbuffer.data(pg_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_z_y_yyyz = cbuffer.data(pg_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_z_y_yyzz = cbuffer.data(pg_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_z_y_yzzz = cbuffer.data(pg_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_z_y_zzzz = cbuffer.data(pg_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_z_z_xxxx = cbuffer.data(pg_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_z_z_xxxy = cbuffer.data(pg_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_z_z_xxxz = cbuffer.data(pg_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_z_z_xxyy = cbuffer.data(pg_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_z_z_xxyz = cbuffer.data(pg_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_z_z_xxzz = cbuffer.data(pg_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_z_z_xyyy = cbuffer.data(pg_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_z_z_xyyz = cbuffer.data(pg_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_z_z_xyzz = cbuffer.data(pg_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_z_z_xzzz = cbuffer.data(pg_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_z_z_yyyy = cbuffer.data(pg_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_z_z_yyyz = cbuffer.data(pg_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_z_z_yyzz = cbuffer.data(pg_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_z_z_yzzz = cbuffer.data(pg_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_z_z_zzzz = cbuffer.data(pg_geom_11_off + 269 * ccomps * dcomps);

            auto g_z_x_x_xxxx = cbuffer.data(pg_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_x_x_xxxy = cbuffer.data(pg_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_x_x_xxxz = cbuffer.data(pg_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_x_x_xxyy = cbuffer.data(pg_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_x_x_xxyz = cbuffer.data(pg_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_x_x_xxzz = cbuffer.data(pg_geom_11_off + 275 * ccomps * dcomps);

            auto g_z_x_x_xyyy = cbuffer.data(pg_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_x_x_xyyz = cbuffer.data(pg_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_x_x_xyzz = cbuffer.data(pg_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_x_x_xzzz = cbuffer.data(pg_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_x_x_yyyy = cbuffer.data(pg_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_x_x_yyyz = cbuffer.data(pg_geom_11_off + 281 * ccomps * dcomps);

            auto g_z_x_x_yyzz = cbuffer.data(pg_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_x_x_yzzz = cbuffer.data(pg_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_x_x_zzzz = cbuffer.data(pg_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_x_y_xxxx = cbuffer.data(pg_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_x_y_xxxy = cbuffer.data(pg_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_x_y_xxxz = cbuffer.data(pg_geom_11_off + 287 * ccomps * dcomps);

            auto g_z_x_y_xxyy = cbuffer.data(pg_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_x_y_xxyz = cbuffer.data(pg_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_x_y_xxzz = cbuffer.data(pg_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_x_y_xyyy = cbuffer.data(pg_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_x_y_xyyz = cbuffer.data(pg_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_x_y_xyzz = cbuffer.data(pg_geom_11_off + 293 * ccomps * dcomps);

            auto g_z_x_y_xzzz = cbuffer.data(pg_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_x_y_yyyy = cbuffer.data(pg_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_x_y_yyyz = cbuffer.data(pg_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_x_y_yyzz = cbuffer.data(pg_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_x_y_yzzz = cbuffer.data(pg_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_x_y_zzzz = cbuffer.data(pg_geom_11_off + 299 * ccomps * dcomps);

            auto g_z_x_z_xxxx = cbuffer.data(pg_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_x_z_xxxy = cbuffer.data(pg_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_x_z_xxxz = cbuffer.data(pg_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_x_z_xxyy = cbuffer.data(pg_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_x_z_xxyz = cbuffer.data(pg_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_x_z_xxzz = cbuffer.data(pg_geom_11_off + 305 * ccomps * dcomps);

            auto g_z_x_z_xyyy = cbuffer.data(pg_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_x_z_xyyz = cbuffer.data(pg_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_x_z_xyzz = cbuffer.data(pg_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_x_z_xzzz = cbuffer.data(pg_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_x_z_yyyy = cbuffer.data(pg_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_x_z_yyyz = cbuffer.data(pg_geom_11_off + 311 * ccomps * dcomps);

            auto g_z_x_z_yyzz = cbuffer.data(pg_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_x_z_yzzz = cbuffer.data(pg_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_x_z_zzzz = cbuffer.data(pg_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_y_x_xxxx = cbuffer.data(pg_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_y_x_xxxy = cbuffer.data(pg_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_y_x_xxxz = cbuffer.data(pg_geom_11_off + 317 * ccomps * dcomps);

            auto g_z_y_x_xxyy = cbuffer.data(pg_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_y_x_xxyz = cbuffer.data(pg_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_y_x_xxzz = cbuffer.data(pg_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_y_x_xyyy = cbuffer.data(pg_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_y_x_xyyz = cbuffer.data(pg_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_y_x_xyzz = cbuffer.data(pg_geom_11_off + 323 * ccomps * dcomps);

            auto g_z_y_x_xzzz = cbuffer.data(pg_geom_11_off + 324 * ccomps * dcomps);

            auto g_z_y_x_yyyy = cbuffer.data(pg_geom_11_off + 325 * ccomps * dcomps);

            auto g_z_y_x_yyyz = cbuffer.data(pg_geom_11_off + 326 * ccomps * dcomps);

            auto g_z_y_x_yyzz = cbuffer.data(pg_geom_11_off + 327 * ccomps * dcomps);

            auto g_z_y_x_yzzz = cbuffer.data(pg_geom_11_off + 328 * ccomps * dcomps);

            auto g_z_y_x_zzzz = cbuffer.data(pg_geom_11_off + 329 * ccomps * dcomps);

            auto g_z_y_y_xxxx = cbuffer.data(pg_geom_11_off + 330 * ccomps * dcomps);

            auto g_z_y_y_xxxy = cbuffer.data(pg_geom_11_off + 331 * ccomps * dcomps);

            auto g_z_y_y_xxxz = cbuffer.data(pg_geom_11_off + 332 * ccomps * dcomps);

            auto g_z_y_y_xxyy = cbuffer.data(pg_geom_11_off + 333 * ccomps * dcomps);

            auto g_z_y_y_xxyz = cbuffer.data(pg_geom_11_off + 334 * ccomps * dcomps);

            auto g_z_y_y_xxzz = cbuffer.data(pg_geom_11_off + 335 * ccomps * dcomps);

            auto g_z_y_y_xyyy = cbuffer.data(pg_geom_11_off + 336 * ccomps * dcomps);

            auto g_z_y_y_xyyz = cbuffer.data(pg_geom_11_off + 337 * ccomps * dcomps);

            auto g_z_y_y_xyzz = cbuffer.data(pg_geom_11_off + 338 * ccomps * dcomps);

            auto g_z_y_y_xzzz = cbuffer.data(pg_geom_11_off + 339 * ccomps * dcomps);

            auto g_z_y_y_yyyy = cbuffer.data(pg_geom_11_off + 340 * ccomps * dcomps);

            auto g_z_y_y_yyyz = cbuffer.data(pg_geom_11_off + 341 * ccomps * dcomps);

            auto g_z_y_y_yyzz = cbuffer.data(pg_geom_11_off + 342 * ccomps * dcomps);

            auto g_z_y_y_yzzz = cbuffer.data(pg_geom_11_off + 343 * ccomps * dcomps);

            auto g_z_y_y_zzzz = cbuffer.data(pg_geom_11_off + 344 * ccomps * dcomps);

            auto g_z_y_z_xxxx = cbuffer.data(pg_geom_11_off + 345 * ccomps * dcomps);

            auto g_z_y_z_xxxy = cbuffer.data(pg_geom_11_off + 346 * ccomps * dcomps);

            auto g_z_y_z_xxxz = cbuffer.data(pg_geom_11_off + 347 * ccomps * dcomps);

            auto g_z_y_z_xxyy = cbuffer.data(pg_geom_11_off + 348 * ccomps * dcomps);

            auto g_z_y_z_xxyz = cbuffer.data(pg_geom_11_off + 349 * ccomps * dcomps);

            auto g_z_y_z_xxzz = cbuffer.data(pg_geom_11_off + 350 * ccomps * dcomps);

            auto g_z_y_z_xyyy = cbuffer.data(pg_geom_11_off + 351 * ccomps * dcomps);

            auto g_z_y_z_xyyz = cbuffer.data(pg_geom_11_off + 352 * ccomps * dcomps);

            auto g_z_y_z_xyzz = cbuffer.data(pg_geom_11_off + 353 * ccomps * dcomps);

            auto g_z_y_z_xzzz = cbuffer.data(pg_geom_11_off + 354 * ccomps * dcomps);

            auto g_z_y_z_yyyy = cbuffer.data(pg_geom_11_off + 355 * ccomps * dcomps);

            auto g_z_y_z_yyyz = cbuffer.data(pg_geom_11_off + 356 * ccomps * dcomps);

            auto g_z_y_z_yyzz = cbuffer.data(pg_geom_11_off + 357 * ccomps * dcomps);

            auto g_z_y_z_yzzz = cbuffer.data(pg_geom_11_off + 358 * ccomps * dcomps);

            auto g_z_y_z_zzzz = cbuffer.data(pg_geom_11_off + 359 * ccomps * dcomps);

            auto g_z_z_x_xxxx = cbuffer.data(pg_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_z_x_xxxy = cbuffer.data(pg_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_z_x_xxxz = cbuffer.data(pg_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_z_x_xxyy = cbuffer.data(pg_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_z_x_xxyz = cbuffer.data(pg_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_z_x_xxzz = cbuffer.data(pg_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_z_x_xyyy = cbuffer.data(pg_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_z_x_xyyz = cbuffer.data(pg_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_z_x_xyzz = cbuffer.data(pg_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_z_x_xzzz = cbuffer.data(pg_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_z_x_yyyy = cbuffer.data(pg_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_z_x_yyyz = cbuffer.data(pg_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_z_x_yyzz = cbuffer.data(pg_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_z_x_yzzz = cbuffer.data(pg_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_z_x_zzzz = cbuffer.data(pg_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_z_y_xxxx = cbuffer.data(pg_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_z_y_xxxy = cbuffer.data(pg_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_z_y_xxxz = cbuffer.data(pg_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_z_y_xxyy = cbuffer.data(pg_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_z_y_xxyz = cbuffer.data(pg_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_z_y_xxzz = cbuffer.data(pg_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_z_y_xyyy = cbuffer.data(pg_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_z_y_xyyz = cbuffer.data(pg_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_z_y_xyzz = cbuffer.data(pg_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_z_y_xzzz = cbuffer.data(pg_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_z_y_yyyy = cbuffer.data(pg_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_z_y_yyyz = cbuffer.data(pg_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_z_y_yyzz = cbuffer.data(pg_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_z_y_yzzz = cbuffer.data(pg_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_z_y_zzzz = cbuffer.data(pg_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_z_z_xxxx = cbuffer.data(pg_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_z_z_xxxy = cbuffer.data(pg_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_z_z_xxxz = cbuffer.data(pg_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_z_z_xxyy = cbuffer.data(pg_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_z_z_xxyz = cbuffer.data(pg_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_z_z_xxzz = cbuffer.data(pg_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_z_z_xyyy = cbuffer.data(pg_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_z_z_xyyz = cbuffer.data(pg_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_z_z_xyzz = cbuffer.data(pg_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_z_z_xzzz = cbuffer.data(pg_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_z_z_yyyy = cbuffer.data(pg_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_z_z_yyyz = cbuffer.data(pg_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_z_z_yyzz = cbuffer.data(pg_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_z_z_yzzz = cbuffer.data(pg_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_z_z_zzzz = cbuffer.data(pg_geom_11_off + 404 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_11_off = idx_geom_11_phxx + i * dcomps + j;

            auto g_x_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 409 * ccomps * dcomps);

            auto g_z_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 419 * ccomps * dcomps);

            auto g_z_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 429 * ccomps * dcomps);

            auto g_z_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 439 * ccomps * dcomps);

            auto g_z_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 440 * ccomps * dcomps);

            auto g_z_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 449 * ccomps * dcomps);

            auto g_z_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 459 * ccomps * dcomps);

            auto g_z_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 461 * ccomps * dcomps);

            auto g_z_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 469 * ccomps * dcomps);

            auto g_z_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 479 * ccomps * dcomps);

            auto g_z_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 482 * ccomps * dcomps);

            auto g_z_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 489 * ccomps * dcomps);

            auto g_z_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 499 * ccomps * dcomps);

            auto g_z_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 503 * ccomps * dcomps);

            auto g_z_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 566 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dgxx

            const auto dg_geom_11_off = idx_geom_11_dgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxx, g_0_x_x_xxxy, g_0_x_x_xxxz, g_0_x_x_xxyy, g_0_x_x_xxyz, g_0_x_x_xxzz, g_0_x_x_xyyy, g_0_x_x_xyyz, g_0_x_x_xyzz, g_0_x_x_xzzz, g_0_x_x_yyyy, g_0_x_x_yyyz, g_0_x_x_yyzz, g_0_x_x_yzzz, g_0_x_x_zzzz, g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz, g_x_x_x_xxxx, g_x_x_x_xxxxx, g_x_x_x_xxxxy, g_x_x_x_xxxxz, g_x_x_x_xxxy, g_x_x_x_xxxyy, g_x_x_x_xxxyz, g_x_x_x_xxxz, g_x_x_x_xxxzz, g_x_x_x_xxyy, g_x_x_x_xxyyy, g_x_x_x_xxyyz, g_x_x_x_xxyz, g_x_x_x_xxyzz, g_x_x_x_xxzz, g_x_x_x_xxzzz, g_x_x_x_xyyy, g_x_x_x_xyyyy, g_x_x_x_xyyyz, g_x_x_x_xyyz, g_x_x_x_xyyzz, g_x_x_x_xyzz, g_x_x_x_xyzzz, g_x_x_x_xzzz, g_x_x_x_xzzzz, g_x_x_x_yyyy, g_x_x_x_yyyz, g_x_x_x_yyzz, g_x_x_x_yzzz, g_x_x_x_zzzz, g_x_x_xx_xxxx, g_x_x_xx_xxxy, g_x_x_xx_xxxz, g_x_x_xx_xxyy, g_x_x_xx_xxyz, g_x_x_xx_xxzz, g_x_x_xx_xyyy, g_x_x_xx_xyyz, g_x_x_xx_xyzz, g_x_x_xx_xzzz, g_x_x_xx_yyyy, g_x_x_xx_yyyz, g_x_x_xx_yyzz, g_x_x_xx_yzzz, g_x_x_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_xxxx[k] = -g_0_x_x_xxxx[k] + g_x_0_x_xxxx[k] - g_x_x_x_xxxx[k] * ab_x + g_x_x_x_xxxxx[k];

                g_x_x_xx_xxxy[k] = -g_0_x_x_xxxy[k] + g_x_0_x_xxxy[k] - g_x_x_x_xxxy[k] * ab_x + g_x_x_x_xxxxy[k];

                g_x_x_xx_xxxz[k] = -g_0_x_x_xxxz[k] + g_x_0_x_xxxz[k] - g_x_x_x_xxxz[k] * ab_x + g_x_x_x_xxxxz[k];

                g_x_x_xx_xxyy[k] = -g_0_x_x_xxyy[k] + g_x_0_x_xxyy[k] - g_x_x_x_xxyy[k] * ab_x + g_x_x_x_xxxyy[k];

                g_x_x_xx_xxyz[k] = -g_0_x_x_xxyz[k] + g_x_0_x_xxyz[k] - g_x_x_x_xxyz[k] * ab_x + g_x_x_x_xxxyz[k];

                g_x_x_xx_xxzz[k] = -g_0_x_x_xxzz[k] + g_x_0_x_xxzz[k] - g_x_x_x_xxzz[k] * ab_x + g_x_x_x_xxxzz[k];

                g_x_x_xx_xyyy[k] = -g_0_x_x_xyyy[k] + g_x_0_x_xyyy[k] - g_x_x_x_xyyy[k] * ab_x + g_x_x_x_xxyyy[k];

                g_x_x_xx_xyyz[k] = -g_0_x_x_xyyz[k] + g_x_0_x_xyyz[k] - g_x_x_x_xyyz[k] * ab_x + g_x_x_x_xxyyz[k];

                g_x_x_xx_xyzz[k] = -g_0_x_x_xyzz[k] + g_x_0_x_xyzz[k] - g_x_x_x_xyzz[k] * ab_x + g_x_x_x_xxyzz[k];

                g_x_x_xx_xzzz[k] = -g_0_x_x_xzzz[k] + g_x_0_x_xzzz[k] - g_x_x_x_xzzz[k] * ab_x + g_x_x_x_xxzzz[k];

                g_x_x_xx_yyyy[k] = -g_0_x_x_yyyy[k] + g_x_0_x_yyyy[k] - g_x_x_x_yyyy[k] * ab_x + g_x_x_x_xyyyy[k];

                g_x_x_xx_yyyz[k] = -g_0_x_x_yyyz[k] + g_x_0_x_yyyz[k] - g_x_x_x_yyyz[k] * ab_x + g_x_x_x_xyyyz[k];

                g_x_x_xx_yyzz[k] = -g_0_x_x_yyzz[k] + g_x_0_x_yyzz[k] - g_x_x_x_yyzz[k] * ab_x + g_x_x_x_xyyzz[k];

                g_x_x_xx_yzzz[k] = -g_0_x_x_yzzz[k] + g_x_0_x_yzzz[k] - g_x_x_x_yzzz[k] * ab_x + g_x_x_x_xyzzz[k];

                g_x_x_xx_zzzz[k] = -g_0_x_x_zzzz[k] + g_x_0_x_zzzz[k] - g_x_x_x_zzzz[k] * ab_x + g_x_x_x_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxx, g_x_x_x_xxxxy, g_x_x_x_xxxy, g_x_x_x_xxxyy, g_x_x_x_xxxyz, g_x_x_x_xxxz, g_x_x_x_xxyy, g_x_x_x_xxyyy, g_x_x_x_xxyyz, g_x_x_x_xxyz, g_x_x_x_xxyzz, g_x_x_x_xxzz, g_x_x_x_xyyy, g_x_x_x_xyyyy, g_x_x_x_xyyyz, g_x_x_x_xyyz, g_x_x_x_xyyzz, g_x_x_x_xyzz, g_x_x_x_xyzzz, g_x_x_x_xzzz, g_x_x_x_yyyy, g_x_x_x_yyyyy, g_x_x_x_yyyyz, g_x_x_x_yyyz, g_x_x_x_yyyzz, g_x_x_x_yyzz, g_x_x_x_yyzzz, g_x_x_x_yzzz, g_x_x_x_yzzzz, g_x_x_x_zzzz, g_x_x_xy_xxxx, g_x_x_xy_xxxy, g_x_x_xy_xxxz, g_x_x_xy_xxyy, g_x_x_xy_xxyz, g_x_x_xy_xxzz, g_x_x_xy_xyyy, g_x_x_xy_xyyz, g_x_x_xy_xyzz, g_x_x_xy_xzzz, g_x_x_xy_yyyy, g_x_x_xy_yyyz, g_x_x_xy_yyzz, g_x_x_xy_yzzz, g_x_x_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_xxxx[k] = -g_x_x_x_xxxx[k] * ab_y + g_x_x_x_xxxxy[k];

                g_x_x_xy_xxxy[k] = -g_x_x_x_xxxy[k] * ab_y + g_x_x_x_xxxyy[k];

                g_x_x_xy_xxxz[k] = -g_x_x_x_xxxz[k] * ab_y + g_x_x_x_xxxyz[k];

                g_x_x_xy_xxyy[k] = -g_x_x_x_xxyy[k] * ab_y + g_x_x_x_xxyyy[k];

                g_x_x_xy_xxyz[k] = -g_x_x_x_xxyz[k] * ab_y + g_x_x_x_xxyyz[k];

                g_x_x_xy_xxzz[k] = -g_x_x_x_xxzz[k] * ab_y + g_x_x_x_xxyzz[k];

                g_x_x_xy_xyyy[k] = -g_x_x_x_xyyy[k] * ab_y + g_x_x_x_xyyyy[k];

                g_x_x_xy_xyyz[k] = -g_x_x_x_xyyz[k] * ab_y + g_x_x_x_xyyyz[k];

                g_x_x_xy_xyzz[k] = -g_x_x_x_xyzz[k] * ab_y + g_x_x_x_xyyzz[k];

                g_x_x_xy_xzzz[k] = -g_x_x_x_xzzz[k] * ab_y + g_x_x_x_xyzzz[k];

                g_x_x_xy_yyyy[k] = -g_x_x_x_yyyy[k] * ab_y + g_x_x_x_yyyyy[k];

                g_x_x_xy_yyyz[k] = -g_x_x_x_yyyz[k] * ab_y + g_x_x_x_yyyyz[k];

                g_x_x_xy_yyzz[k] = -g_x_x_x_yyzz[k] * ab_y + g_x_x_x_yyyzz[k];

                g_x_x_xy_yzzz[k] = -g_x_x_x_yzzz[k] * ab_y + g_x_x_x_yyzzz[k];

                g_x_x_xy_zzzz[k] = -g_x_x_x_zzzz[k] * ab_y + g_x_x_x_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxx, g_x_x_x_xxxxz, g_x_x_x_xxxy, g_x_x_x_xxxyz, g_x_x_x_xxxz, g_x_x_x_xxxzz, g_x_x_x_xxyy, g_x_x_x_xxyyz, g_x_x_x_xxyz, g_x_x_x_xxyzz, g_x_x_x_xxzz, g_x_x_x_xxzzz, g_x_x_x_xyyy, g_x_x_x_xyyyz, g_x_x_x_xyyz, g_x_x_x_xyyzz, g_x_x_x_xyzz, g_x_x_x_xyzzz, g_x_x_x_xzzz, g_x_x_x_xzzzz, g_x_x_x_yyyy, g_x_x_x_yyyyz, g_x_x_x_yyyz, g_x_x_x_yyyzz, g_x_x_x_yyzz, g_x_x_x_yyzzz, g_x_x_x_yzzz, g_x_x_x_yzzzz, g_x_x_x_zzzz, g_x_x_x_zzzzz, g_x_x_xz_xxxx, g_x_x_xz_xxxy, g_x_x_xz_xxxz, g_x_x_xz_xxyy, g_x_x_xz_xxyz, g_x_x_xz_xxzz, g_x_x_xz_xyyy, g_x_x_xz_xyyz, g_x_x_xz_xyzz, g_x_x_xz_xzzz, g_x_x_xz_yyyy, g_x_x_xz_yyyz, g_x_x_xz_yyzz, g_x_x_xz_yzzz, g_x_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_xxxx[k] = -g_x_x_x_xxxx[k] * ab_z + g_x_x_x_xxxxz[k];

                g_x_x_xz_xxxy[k] = -g_x_x_x_xxxy[k] * ab_z + g_x_x_x_xxxyz[k];

                g_x_x_xz_xxxz[k] = -g_x_x_x_xxxz[k] * ab_z + g_x_x_x_xxxzz[k];

                g_x_x_xz_xxyy[k] = -g_x_x_x_xxyy[k] * ab_z + g_x_x_x_xxyyz[k];

                g_x_x_xz_xxyz[k] = -g_x_x_x_xxyz[k] * ab_z + g_x_x_x_xxyzz[k];

                g_x_x_xz_xxzz[k] = -g_x_x_x_xxzz[k] * ab_z + g_x_x_x_xxzzz[k];

                g_x_x_xz_xyyy[k] = -g_x_x_x_xyyy[k] * ab_z + g_x_x_x_xyyyz[k];

                g_x_x_xz_xyyz[k] = -g_x_x_x_xyyz[k] * ab_z + g_x_x_x_xyyzz[k];

                g_x_x_xz_xyzz[k] = -g_x_x_x_xyzz[k] * ab_z + g_x_x_x_xyzzz[k];

                g_x_x_xz_xzzz[k] = -g_x_x_x_xzzz[k] * ab_z + g_x_x_x_xzzzz[k];

                g_x_x_xz_yyyy[k] = -g_x_x_x_yyyy[k] * ab_z + g_x_x_x_yyyyz[k];

                g_x_x_xz_yyyz[k] = -g_x_x_x_yyyz[k] * ab_z + g_x_x_x_yyyzz[k];

                g_x_x_xz_yyzz[k] = -g_x_x_x_yyzz[k] * ab_z + g_x_x_x_yyzzz[k];

                g_x_x_xz_yzzz[k] = -g_x_x_x_yzzz[k] * ab_z + g_x_x_x_yzzzz[k];

                g_x_x_xz_zzzz[k] = -g_x_x_x_zzzz[k] * ab_z + g_x_x_x_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_xxxx, g_x_x_y_xxxxy, g_x_x_y_xxxy, g_x_x_y_xxxyy, g_x_x_y_xxxyz, g_x_x_y_xxxz, g_x_x_y_xxyy, g_x_x_y_xxyyy, g_x_x_y_xxyyz, g_x_x_y_xxyz, g_x_x_y_xxyzz, g_x_x_y_xxzz, g_x_x_y_xyyy, g_x_x_y_xyyyy, g_x_x_y_xyyyz, g_x_x_y_xyyz, g_x_x_y_xyyzz, g_x_x_y_xyzz, g_x_x_y_xyzzz, g_x_x_y_xzzz, g_x_x_y_yyyy, g_x_x_y_yyyyy, g_x_x_y_yyyyz, g_x_x_y_yyyz, g_x_x_y_yyyzz, g_x_x_y_yyzz, g_x_x_y_yyzzz, g_x_x_y_yzzz, g_x_x_y_yzzzz, g_x_x_y_zzzz, g_x_x_yy_xxxx, g_x_x_yy_xxxy, g_x_x_yy_xxxz, g_x_x_yy_xxyy, g_x_x_yy_xxyz, g_x_x_yy_xxzz, g_x_x_yy_xyyy, g_x_x_yy_xyyz, g_x_x_yy_xyzz, g_x_x_yy_xzzz, g_x_x_yy_yyyy, g_x_x_yy_yyyz, g_x_x_yy_yyzz, g_x_x_yy_yzzz, g_x_x_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_xxxx[k] = -g_x_x_y_xxxx[k] * ab_y + g_x_x_y_xxxxy[k];

                g_x_x_yy_xxxy[k] = -g_x_x_y_xxxy[k] * ab_y + g_x_x_y_xxxyy[k];

                g_x_x_yy_xxxz[k] = -g_x_x_y_xxxz[k] * ab_y + g_x_x_y_xxxyz[k];

                g_x_x_yy_xxyy[k] = -g_x_x_y_xxyy[k] * ab_y + g_x_x_y_xxyyy[k];

                g_x_x_yy_xxyz[k] = -g_x_x_y_xxyz[k] * ab_y + g_x_x_y_xxyyz[k];

                g_x_x_yy_xxzz[k] = -g_x_x_y_xxzz[k] * ab_y + g_x_x_y_xxyzz[k];

                g_x_x_yy_xyyy[k] = -g_x_x_y_xyyy[k] * ab_y + g_x_x_y_xyyyy[k];

                g_x_x_yy_xyyz[k] = -g_x_x_y_xyyz[k] * ab_y + g_x_x_y_xyyyz[k];

                g_x_x_yy_xyzz[k] = -g_x_x_y_xyzz[k] * ab_y + g_x_x_y_xyyzz[k];

                g_x_x_yy_xzzz[k] = -g_x_x_y_xzzz[k] * ab_y + g_x_x_y_xyzzz[k];

                g_x_x_yy_yyyy[k] = -g_x_x_y_yyyy[k] * ab_y + g_x_x_y_yyyyy[k];

                g_x_x_yy_yyyz[k] = -g_x_x_y_yyyz[k] * ab_y + g_x_x_y_yyyyz[k];

                g_x_x_yy_yyzz[k] = -g_x_x_y_yyzz[k] * ab_y + g_x_x_y_yyyzz[k];

                g_x_x_yy_yzzz[k] = -g_x_x_y_yzzz[k] * ab_y + g_x_x_y_yyzzz[k];

                g_x_x_yy_zzzz[k] = -g_x_x_y_zzzz[k] * ab_y + g_x_x_y_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_xxxx, g_x_x_yz_xxxy, g_x_x_yz_xxxz, g_x_x_yz_xxyy, g_x_x_yz_xxyz, g_x_x_yz_xxzz, g_x_x_yz_xyyy, g_x_x_yz_xyyz, g_x_x_yz_xyzz, g_x_x_yz_xzzz, g_x_x_yz_yyyy, g_x_x_yz_yyyz, g_x_x_yz_yyzz, g_x_x_yz_yzzz, g_x_x_yz_zzzz, g_x_x_z_xxxx, g_x_x_z_xxxxy, g_x_x_z_xxxy, g_x_x_z_xxxyy, g_x_x_z_xxxyz, g_x_x_z_xxxz, g_x_x_z_xxyy, g_x_x_z_xxyyy, g_x_x_z_xxyyz, g_x_x_z_xxyz, g_x_x_z_xxyzz, g_x_x_z_xxzz, g_x_x_z_xyyy, g_x_x_z_xyyyy, g_x_x_z_xyyyz, g_x_x_z_xyyz, g_x_x_z_xyyzz, g_x_x_z_xyzz, g_x_x_z_xyzzz, g_x_x_z_xzzz, g_x_x_z_yyyy, g_x_x_z_yyyyy, g_x_x_z_yyyyz, g_x_x_z_yyyz, g_x_x_z_yyyzz, g_x_x_z_yyzz, g_x_x_z_yyzzz, g_x_x_z_yzzz, g_x_x_z_yzzzz, g_x_x_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_xxxx[k] = -g_x_x_z_xxxx[k] * ab_y + g_x_x_z_xxxxy[k];

                g_x_x_yz_xxxy[k] = -g_x_x_z_xxxy[k] * ab_y + g_x_x_z_xxxyy[k];

                g_x_x_yz_xxxz[k] = -g_x_x_z_xxxz[k] * ab_y + g_x_x_z_xxxyz[k];

                g_x_x_yz_xxyy[k] = -g_x_x_z_xxyy[k] * ab_y + g_x_x_z_xxyyy[k];

                g_x_x_yz_xxyz[k] = -g_x_x_z_xxyz[k] * ab_y + g_x_x_z_xxyyz[k];

                g_x_x_yz_xxzz[k] = -g_x_x_z_xxzz[k] * ab_y + g_x_x_z_xxyzz[k];

                g_x_x_yz_xyyy[k] = -g_x_x_z_xyyy[k] * ab_y + g_x_x_z_xyyyy[k];

                g_x_x_yz_xyyz[k] = -g_x_x_z_xyyz[k] * ab_y + g_x_x_z_xyyyz[k];

                g_x_x_yz_xyzz[k] = -g_x_x_z_xyzz[k] * ab_y + g_x_x_z_xyyzz[k];

                g_x_x_yz_xzzz[k] = -g_x_x_z_xzzz[k] * ab_y + g_x_x_z_xyzzz[k];

                g_x_x_yz_yyyy[k] = -g_x_x_z_yyyy[k] * ab_y + g_x_x_z_yyyyy[k];

                g_x_x_yz_yyyz[k] = -g_x_x_z_yyyz[k] * ab_y + g_x_x_z_yyyyz[k];

                g_x_x_yz_yyzz[k] = -g_x_x_z_yyzz[k] * ab_y + g_x_x_z_yyyzz[k];

                g_x_x_yz_yzzz[k] = -g_x_x_z_yzzz[k] * ab_y + g_x_x_z_yyzzz[k];

                g_x_x_yz_zzzz[k] = -g_x_x_z_zzzz[k] * ab_y + g_x_x_z_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_xxxx, g_x_x_z_xxxxz, g_x_x_z_xxxy, g_x_x_z_xxxyz, g_x_x_z_xxxz, g_x_x_z_xxxzz, g_x_x_z_xxyy, g_x_x_z_xxyyz, g_x_x_z_xxyz, g_x_x_z_xxyzz, g_x_x_z_xxzz, g_x_x_z_xxzzz, g_x_x_z_xyyy, g_x_x_z_xyyyz, g_x_x_z_xyyz, g_x_x_z_xyyzz, g_x_x_z_xyzz, g_x_x_z_xyzzz, g_x_x_z_xzzz, g_x_x_z_xzzzz, g_x_x_z_yyyy, g_x_x_z_yyyyz, g_x_x_z_yyyz, g_x_x_z_yyyzz, g_x_x_z_yyzz, g_x_x_z_yyzzz, g_x_x_z_yzzz, g_x_x_z_yzzzz, g_x_x_z_zzzz, g_x_x_z_zzzzz, g_x_x_zz_xxxx, g_x_x_zz_xxxy, g_x_x_zz_xxxz, g_x_x_zz_xxyy, g_x_x_zz_xxyz, g_x_x_zz_xxzz, g_x_x_zz_xyyy, g_x_x_zz_xyyz, g_x_x_zz_xyzz, g_x_x_zz_xzzz, g_x_x_zz_yyyy, g_x_x_zz_yyyz, g_x_x_zz_yyzz, g_x_x_zz_yzzz, g_x_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_xxxx[k] = -g_x_x_z_xxxx[k] * ab_z + g_x_x_z_xxxxz[k];

                g_x_x_zz_xxxy[k] = -g_x_x_z_xxxy[k] * ab_z + g_x_x_z_xxxyz[k];

                g_x_x_zz_xxxz[k] = -g_x_x_z_xxxz[k] * ab_z + g_x_x_z_xxxzz[k];

                g_x_x_zz_xxyy[k] = -g_x_x_z_xxyy[k] * ab_z + g_x_x_z_xxyyz[k];

                g_x_x_zz_xxyz[k] = -g_x_x_z_xxyz[k] * ab_z + g_x_x_z_xxyzz[k];

                g_x_x_zz_xxzz[k] = -g_x_x_z_xxzz[k] * ab_z + g_x_x_z_xxzzz[k];

                g_x_x_zz_xyyy[k] = -g_x_x_z_xyyy[k] * ab_z + g_x_x_z_xyyyz[k];

                g_x_x_zz_xyyz[k] = -g_x_x_z_xyyz[k] * ab_z + g_x_x_z_xyyzz[k];

                g_x_x_zz_xyzz[k] = -g_x_x_z_xyzz[k] * ab_z + g_x_x_z_xyzzz[k];

                g_x_x_zz_xzzz[k] = -g_x_x_z_xzzz[k] * ab_z + g_x_x_z_xzzzz[k];

                g_x_x_zz_yyyy[k] = -g_x_x_z_yyyy[k] * ab_z + g_x_x_z_yyyyz[k];

                g_x_x_zz_yyyz[k] = -g_x_x_z_yyyz[k] * ab_z + g_x_x_z_yyyzz[k];

                g_x_x_zz_yyzz[k] = -g_x_x_z_yyzz[k] * ab_z + g_x_x_z_yyzzz[k];

                g_x_x_zz_yzzz[k] = -g_x_x_z_yzzz[k] * ab_z + g_x_x_z_yzzzz[k];

                g_x_x_zz_zzzz[k] = -g_x_x_z_zzzz[k] * ab_z + g_x_x_z_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxxx, g_0_y_x_xxxy, g_0_y_x_xxxz, g_0_y_x_xxyy, g_0_y_x_xxyz, g_0_y_x_xxzz, g_0_y_x_xyyy, g_0_y_x_xyyz, g_0_y_x_xyzz, g_0_y_x_xzzz, g_0_y_x_yyyy, g_0_y_x_yyyz, g_0_y_x_yyzz, g_0_y_x_yzzz, g_0_y_x_zzzz, g_x_y_x_xxxx, g_x_y_x_xxxxx, g_x_y_x_xxxxy, g_x_y_x_xxxxz, g_x_y_x_xxxy, g_x_y_x_xxxyy, g_x_y_x_xxxyz, g_x_y_x_xxxz, g_x_y_x_xxxzz, g_x_y_x_xxyy, g_x_y_x_xxyyy, g_x_y_x_xxyyz, g_x_y_x_xxyz, g_x_y_x_xxyzz, g_x_y_x_xxzz, g_x_y_x_xxzzz, g_x_y_x_xyyy, g_x_y_x_xyyyy, g_x_y_x_xyyyz, g_x_y_x_xyyz, g_x_y_x_xyyzz, g_x_y_x_xyzz, g_x_y_x_xyzzz, g_x_y_x_xzzz, g_x_y_x_xzzzz, g_x_y_x_yyyy, g_x_y_x_yyyz, g_x_y_x_yyzz, g_x_y_x_yzzz, g_x_y_x_zzzz, g_x_y_xx_xxxx, g_x_y_xx_xxxy, g_x_y_xx_xxxz, g_x_y_xx_xxyy, g_x_y_xx_xxyz, g_x_y_xx_xxzz, g_x_y_xx_xyyy, g_x_y_xx_xyyz, g_x_y_xx_xyzz, g_x_y_xx_xzzz, g_x_y_xx_yyyy, g_x_y_xx_yyyz, g_x_y_xx_yyzz, g_x_y_xx_yzzz, g_x_y_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_xxxx[k] = -g_0_y_x_xxxx[k] - g_x_y_x_xxxx[k] * ab_x + g_x_y_x_xxxxx[k];

                g_x_y_xx_xxxy[k] = -g_0_y_x_xxxy[k] - g_x_y_x_xxxy[k] * ab_x + g_x_y_x_xxxxy[k];

                g_x_y_xx_xxxz[k] = -g_0_y_x_xxxz[k] - g_x_y_x_xxxz[k] * ab_x + g_x_y_x_xxxxz[k];

                g_x_y_xx_xxyy[k] = -g_0_y_x_xxyy[k] - g_x_y_x_xxyy[k] * ab_x + g_x_y_x_xxxyy[k];

                g_x_y_xx_xxyz[k] = -g_0_y_x_xxyz[k] - g_x_y_x_xxyz[k] * ab_x + g_x_y_x_xxxyz[k];

                g_x_y_xx_xxzz[k] = -g_0_y_x_xxzz[k] - g_x_y_x_xxzz[k] * ab_x + g_x_y_x_xxxzz[k];

                g_x_y_xx_xyyy[k] = -g_0_y_x_xyyy[k] - g_x_y_x_xyyy[k] * ab_x + g_x_y_x_xxyyy[k];

                g_x_y_xx_xyyz[k] = -g_0_y_x_xyyz[k] - g_x_y_x_xyyz[k] * ab_x + g_x_y_x_xxyyz[k];

                g_x_y_xx_xyzz[k] = -g_0_y_x_xyzz[k] - g_x_y_x_xyzz[k] * ab_x + g_x_y_x_xxyzz[k];

                g_x_y_xx_xzzz[k] = -g_0_y_x_xzzz[k] - g_x_y_x_xzzz[k] * ab_x + g_x_y_x_xxzzz[k];

                g_x_y_xx_yyyy[k] = -g_0_y_x_yyyy[k] - g_x_y_x_yyyy[k] * ab_x + g_x_y_x_xyyyy[k];

                g_x_y_xx_yyyz[k] = -g_0_y_x_yyyz[k] - g_x_y_x_yyyz[k] * ab_x + g_x_y_x_xyyyz[k];

                g_x_y_xx_yyzz[k] = -g_0_y_x_yyzz[k] - g_x_y_x_yyzz[k] * ab_x + g_x_y_x_xyyzz[k];

                g_x_y_xx_yzzz[k] = -g_0_y_x_yzzz[k] - g_x_y_x_yzzz[k] * ab_x + g_x_y_x_xyzzz[k];

                g_x_y_xx_zzzz[k] = -g_0_y_x_zzzz[k] - g_x_y_x_zzzz[k] * ab_x + g_x_y_x_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxx, g_0_y_y_xxxy, g_0_y_y_xxxz, g_0_y_y_xxyy, g_0_y_y_xxyz, g_0_y_y_xxzz, g_0_y_y_xyyy, g_0_y_y_xyyz, g_0_y_y_xyzz, g_0_y_y_xzzz, g_0_y_y_yyyy, g_0_y_y_yyyz, g_0_y_y_yyzz, g_0_y_y_yzzz, g_0_y_y_zzzz, g_x_y_xy_xxxx, g_x_y_xy_xxxy, g_x_y_xy_xxxz, g_x_y_xy_xxyy, g_x_y_xy_xxyz, g_x_y_xy_xxzz, g_x_y_xy_xyyy, g_x_y_xy_xyyz, g_x_y_xy_xyzz, g_x_y_xy_xzzz, g_x_y_xy_yyyy, g_x_y_xy_yyyz, g_x_y_xy_yyzz, g_x_y_xy_yzzz, g_x_y_xy_zzzz, g_x_y_y_xxxx, g_x_y_y_xxxxx, g_x_y_y_xxxxy, g_x_y_y_xxxxz, g_x_y_y_xxxy, g_x_y_y_xxxyy, g_x_y_y_xxxyz, g_x_y_y_xxxz, g_x_y_y_xxxzz, g_x_y_y_xxyy, g_x_y_y_xxyyy, g_x_y_y_xxyyz, g_x_y_y_xxyz, g_x_y_y_xxyzz, g_x_y_y_xxzz, g_x_y_y_xxzzz, g_x_y_y_xyyy, g_x_y_y_xyyyy, g_x_y_y_xyyyz, g_x_y_y_xyyz, g_x_y_y_xyyzz, g_x_y_y_xyzz, g_x_y_y_xyzzz, g_x_y_y_xzzz, g_x_y_y_xzzzz, g_x_y_y_yyyy, g_x_y_y_yyyz, g_x_y_y_yyzz, g_x_y_y_yzzz, g_x_y_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_xxxx[k] = -g_0_y_y_xxxx[k] - g_x_y_y_xxxx[k] * ab_x + g_x_y_y_xxxxx[k];

                g_x_y_xy_xxxy[k] = -g_0_y_y_xxxy[k] - g_x_y_y_xxxy[k] * ab_x + g_x_y_y_xxxxy[k];

                g_x_y_xy_xxxz[k] = -g_0_y_y_xxxz[k] - g_x_y_y_xxxz[k] * ab_x + g_x_y_y_xxxxz[k];

                g_x_y_xy_xxyy[k] = -g_0_y_y_xxyy[k] - g_x_y_y_xxyy[k] * ab_x + g_x_y_y_xxxyy[k];

                g_x_y_xy_xxyz[k] = -g_0_y_y_xxyz[k] - g_x_y_y_xxyz[k] * ab_x + g_x_y_y_xxxyz[k];

                g_x_y_xy_xxzz[k] = -g_0_y_y_xxzz[k] - g_x_y_y_xxzz[k] * ab_x + g_x_y_y_xxxzz[k];

                g_x_y_xy_xyyy[k] = -g_0_y_y_xyyy[k] - g_x_y_y_xyyy[k] * ab_x + g_x_y_y_xxyyy[k];

                g_x_y_xy_xyyz[k] = -g_0_y_y_xyyz[k] - g_x_y_y_xyyz[k] * ab_x + g_x_y_y_xxyyz[k];

                g_x_y_xy_xyzz[k] = -g_0_y_y_xyzz[k] - g_x_y_y_xyzz[k] * ab_x + g_x_y_y_xxyzz[k];

                g_x_y_xy_xzzz[k] = -g_0_y_y_xzzz[k] - g_x_y_y_xzzz[k] * ab_x + g_x_y_y_xxzzz[k];

                g_x_y_xy_yyyy[k] = -g_0_y_y_yyyy[k] - g_x_y_y_yyyy[k] * ab_x + g_x_y_y_xyyyy[k];

                g_x_y_xy_yyyz[k] = -g_0_y_y_yyyz[k] - g_x_y_y_yyyz[k] * ab_x + g_x_y_y_xyyyz[k];

                g_x_y_xy_yyzz[k] = -g_0_y_y_yyzz[k] - g_x_y_y_yyzz[k] * ab_x + g_x_y_y_xyyzz[k];

                g_x_y_xy_yzzz[k] = -g_0_y_y_yzzz[k] - g_x_y_y_yzzz[k] * ab_x + g_x_y_y_xyzzz[k];

                g_x_y_xy_zzzz[k] = -g_0_y_y_zzzz[k] - g_x_y_y_zzzz[k] * ab_x + g_x_y_y_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_xxxx, g_x_y_x_xxxxz, g_x_y_x_xxxy, g_x_y_x_xxxyz, g_x_y_x_xxxz, g_x_y_x_xxxzz, g_x_y_x_xxyy, g_x_y_x_xxyyz, g_x_y_x_xxyz, g_x_y_x_xxyzz, g_x_y_x_xxzz, g_x_y_x_xxzzz, g_x_y_x_xyyy, g_x_y_x_xyyyz, g_x_y_x_xyyz, g_x_y_x_xyyzz, g_x_y_x_xyzz, g_x_y_x_xyzzz, g_x_y_x_xzzz, g_x_y_x_xzzzz, g_x_y_x_yyyy, g_x_y_x_yyyyz, g_x_y_x_yyyz, g_x_y_x_yyyzz, g_x_y_x_yyzz, g_x_y_x_yyzzz, g_x_y_x_yzzz, g_x_y_x_yzzzz, g_x_y_x_zzzz, g_x_y_x_zzzzz, g_x_y_xz_xxxx, g_x_y_xz_xxxy, g_x_y_xz_xxxz, g_x_y_xz_xxyy, g_x_y_xz_xxyz, g_x_y_xz_xxzz, g_x_y_xz_xyyy, g_x_y_xz_xyyz, g_x_y_xz_xyzz, g_x_y_xz_xzzz, g_x_y_xz_yyyy, g_x_y_xz_yyyz, g_x_y_xz_yyzz, g_x_y_xz_yzzz, g_x_y_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_xxxx[k] = -g_x_y_x_xxxx[k] * ab_z + g_x_y_x_xxxxz[k];

                g_x_y_xz_xxxy[k] = -g_x_y_x_xxxy[k] * ab_z + g_x_y_x_xxxyz[k];

                g_x_y_xz_xxxz[k] = -g_x_y_x_xxxz[k] * ab_z + g_x_y_x_xxxzz[k];

                g_x_y_xz_xxyy[k] = -g_x_y_x_xxyy[k] * ab_z + g_x_y_x_xxyyz[k];

                g_x_y_xz_xxyz[k] = -g_x_y_x_xxyz[k] * ab_z + g_x_y_x_xxyzz[k];

                g_x_y_xz_xxzz[k] = -g_x_y_x_xxzz[k] * ab_z + g_x_y_x_xxzzz[k];

                g_x_y_xz_xyyy[k] = -g_x_y_x_xyyy[k] * ab_z + g_x_y_x_xyyyz[k];

                g_x_y_xz_xyyz[k] = -g_x_y_x_xyyz[k] * ab_z + g_x_y_x_xyyzz[k];

                g_x_y_xz_xyzz[k] = -g_x_y_x_xyzz[k] * ab_z + g_x_y_x_xyzzz[k];

                g_x_y_xz_xzzz[k] = -g_x_y_x_xzzz[k] * ab_z + g_x_y_x_xzzzz[k];

                g_x_y_xz_yyyy[k] = -g_x_y_x_yyyy[k] * ab_z + g_x_y_x_yyyyz[k];

                g_x_y_xz_yyyz[k] = -g_x_y_x_yyyz[k] * ab_z + g_x_y_x_yyyzz[k];

                g_x_y_xz_yyzz[k] = -g_x_y_x_yyzz[k] * ab_z + g_x_y_x_yyzzz[k];

                g_x_y_xz_yzzz[k] = -g_x_y_x_yzzz[k] * ab_z + g_x_y_x_yzzzz[k];

                g_x_y_xz_zzzz[k] = -g_x_y_x_zzzz[k] * ab_z + g_x_y_x_zzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxx, g_x_0_y_xxxy, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyzz, g_x_0_y_yzzz, g_x_0_y_zzzz, g_x_y_y_xxxx, g_x_y_y_xxxxy, g_x_y_y_xxxy, g_x_y_y_xxxyy, g_x_y_y_xxxyz, g_x_y_y_xxxz, g_x_y_y_xxyy, g_x_y_y_xxyyy, g_x_y_y_xxyyz, g_x_y_y_xxyz, g_x_y_y_xxyzz, g_x_y_y_xxzz, g_x_y_y_xyyy, g_x_y_y_xyyyy, g_x_y_y_xyyyz, g_x_y_y_xyyz, g_x_y_y_xyyzz, g_x_y_y_xyzz, g_x_y_y_xyzzz, g_x_y_y_xzzz, g_x_y_y_yyyy, g_x_y_y_yyyyy, g_x_y_y_yyyyz, g_x_y_y_yyyz, g_x_y_y_yyyzz, g_x_y_y_yyzz, g_x_y_y_yyzzz, g_x_y_y_yzzz, g_x_y_y_yzzzz, g_x_y_y_zzzz, g_x_y_yy_xxxx, g_x_y_yy_xxxy, g_x_y_yy_xxxz, g_x_y_yy_xxyy, g_x_y_yy_xxyz, g_x_y_yy_xxzz, g_x_y_yy_xyyy, g_x_y_yy_xyyz, g_x_y_yy_xyzz, g_x_y_yy_xzzz, g_x_y_yy_yyyy, g_x_y_yy_yyyz, g_x_y_yy_yyzz, g_x_y_yy_yzzz, g_x_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_xxxx[k] = g_x_0_y_xxxx[k] - g_x_y_y_xxxx[k] * ab_y + g_x_y_y_xxxxy[k];

                g_x_y_yy_xxxy[k] = g_x_0_y_xxxy[k] - g_x_y_y_xxxy[k] * ab_y + g_x_y_y_xxxyy[k];

                g_x_y_yy_xxxz[k] = g_x_0_y_xxxz[k] - g_x_y_y_xxxz[k] * ab_y + g_x_y_y_xxxyz[k];

                g_x_y_yy_xxyy[k] = g_x_0_y_xxyy[k] - g_x_y_y_xxyy[k] * ab_y + g_x_y_y_xxyyy[k];

                g_x_y_yy_xxyz[k] = g_x_0_y_xxyz[k] - g_x_y_y_xxyz[k] * ab_y + g_x_y_y_xxyyz[k];

                g_x_y_yy_xxzz[k] = g_x_0_y_xxzz[k] - g_x_y_y_xxzz[k] * ab_y + g_x_y_y_xxyzz[k];

                g_x_y_yy_xyyy[k] = g_x_0_y_xyyy[k] - g_x_y_y_xyyy[k] * ab_y + g_x_y_y_xyyyy[k];

                g_x_y_yy_xyyz[k] = g_x_0_y_xyyz[k] - g_x_y_y_xyyz[k] * ab_y + g_x_y_y_xyyyz[k];

                g_x_y_yy_xyzz[k] = g_x_0_y_xyzz[k] - g_x_y_y_xyzz[k] * ab_y + g_x_y_y_xyyzz[k];

                g_x_y_yy_xzzz[k] = g_x_0_y_xzzz[k] - g_x_y_y_xzzz[k] * ab_y + g_x_y_y_xyzzz[k];

                g_x_y_yy_yyyy[k] = g_x_0_y_yyyy[k] - g_x_y_y_yyyy[k] * ab_y + g_x_y_y_yyyyy[k];

                g_x_y_yy_yyyz[k] = g_x_0_y_yyyz[k] - g_x_y_y_yyyz[k] * ab_y + g_x_y_y_yyyyz[k];

                g_x_y_yy_yyzz[k] = g_x_0_y_yyzz[k] - g_x_y_y_yyzz[k] * ab_y + g_x_y_y_yyyzz[k];

                g_x_y_yy_yzzz[k] = g_x_0_y_yzzz[k] - g_x_y_y_yzzz[k] * ab_y + g_x_y_y_yyzzz[k];

                g_x_y_yy_zzzz[k] = g_x_0_y_zzzz[k] - g_x_y_y_zzzz[k] * ab_y + g_x_y_y_yzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_xxxx, g_x_y_y_xxxxz, g_x_y_y_xxxy, g_x_y_y_xxxyz, g_x_y_y_xxxz, g_x_y_y_xxxzz, g_x_y_y_xxyy, g_x_y_y_xxyyz, g_x_y_y_xxyz, g_x_y_y_xxyzz, g_x_y_y_xxzz, g_x_y_y_xxzzz, g_x_y_y_xyyy, g_x_y_y_xyyyz, g_x_y_y_xyyz, g_x_y_y_xyyzz, g_x_y_y_xyzz, g_x_y_y_xyzzz, g_x_y_y_xzzz, g_x_y_y_xzzzz, g_x_y_y_yyyy, g_x_y_y_yyyyz, g_x_y_y_yyyz, g_x_y_y_yyyzz, g_x_y_y_yyzz, g_x_y_y_yyzzz, g_x_y_y_yzzz, g_x_y_y_yzzzz, g_x_y_y_zzzz, g_x_y_y_zzzzz, g_x_y_yz_xxxx, g_x_y_yz_xxxy, g_x_y_yz_xxxz, g_x_y_yz_xxyy, g_x_y_yz_xxyz, g_x_y_yz_xxzz, g_x_y_yz_xyyy, g_x_y_yz_xyyz, g_x_y_yz_xyzz, g_x_y_yz_xzzz, g_x_y_yz_yyyy, g_x_y_yz_yyyz, g_x_y_yz_yyzz, g_x_y_yz_yzzz, g_x_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_xxxx[k] = -g_x_y_y_xxxx[k] * ab_z + g_x_y_y_xxxxz[k];

                g_x_y_yz_xxxy[k] = -g_x_y_y_xxxy[k] * ab_z + g_x_y_y_xxxyz[k];

                g_x_y_yz_xxxz[k] = -g_x_y_y_xxxz[k] * ab_z + g_x_y_y_xxxzz[k];

                g_x_y_yz_xxyy[k] = -g_x_y_y_xxyy[k] * ab_z + g_x_y_y_xxyyz[k];

                g_x_y_yz_xxyz[k] = -g_x_y_y_xxyz[k] * ab_z + g_x_y_y_xxyzz[k];

                g_x_y_yz_xxzz[k] = -g_x_y_y_xxzz[k] * ab_z + g_x_y_y_xxzzz[k];

                g_x_y_yz_xyyy[k] = -g_x_y_y_xyyy[k] * ab_z + g_x_y_y_xyyyz[k];

                g_x_y_yz_xyyz[k] = -g_x_y_y_xyyz[k] * ab_z + g_x_y_y_xyyzz[k];

                g_x_y_yz_xyzz[k] = -g_x_y_y_xyzz[k] * ab_z + g_x_y_y_xyzzz[k];

                g_x_y_yz_xzzz[k] = -g_x_y_y_xzzz[k] * ab_z + g_x_y_y_xzzzz[k];

                g_x_y_yz_yyyy[k] = -g_x_y_y_yyyy[k] * ab_z + g_x_y_y_yyyyz[k];

                g_x_y_yz_yyyz[k] = -g_x_y_y_yyyz[k] * ab_z + g_x_y_y_yyyzz[k];

                g_x_y_yz_yyzz[k] = -g_x_y_y_yyzz[k] * ab_z + g_x_y_y_yyzzz[k];

                g_x_y_yz_yzzz[k] = -g_x_y_y_yzzz[k] * ab_z + g_x_y_y_yzzzz[k];

                g_x_y_yz_zzzz[k] = -g_x_y_y_zzzz[k] * ab_z + g_x_y_y_zzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_xxxx, g_x_y_z_xxxxz, g_x_y_z_xxxy, g_x_y_z_xxxyz, g_x_y_z_xxxz, g_x_y_z_xxxzz, g_x_y_z_xxyy, g_x_y_z_xxyyz, g_x_y_z_xxyz, g_x_y_z_xxyzz, g_x_y_z_xxzz, g_x_y_z_xxzzz, g_x_y_z_xyyy, g_x_y_z_xyyyz, g_x_y_z_xyyz, g_x_y_z_xyyzz, g_x_y_z_xyzz, g_x_y_z_xyzzz, g_x_y_z_xzzz, g_x_y_z_xzzzz, g_x_y_z_yyyy, g_x_y_z_yyyyz, g_x_y_z_yyyz, g_x_y_z_yyyzz, g_x_y_z_yyzz, g_x_y_z_yyzzz, g_x_y_z_yzzz, g_x_y_z_yzzzz, g_x_y_z_zzzz, g_x_y_z_zzzzz, g_x_y_zz_xxxx, g_x_y_zz_xxxy, g_x_y_zz_xxxz, g_x_y_zz_xxyy, g_x_y_zz_xxyz, g_x_y_zz_xxzz, g_x_y_zz_xyyy, g_x_y_zz_xyyz, g_x_y_zz_xyzz, g_x_y_zz_xzzz, g_x_y_zz_yyyy, g_x_y_zz_yyyz, g_x_y_zz_yyzz, g_x_y_zz_yzzz, g_x_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_xxxx[k] = -g_x_y_z_xxxx[k] * ab_z + g_x_y_z_xxxxz[k];

                g_x_y_zz_xxxy[k] = -g_x_y_z_xxxy[k] * ab_z + g_x_y_z_xxxyz[k];

                g_x_y_zz_xxxz[k] = -g_x_y_z_xxxz[k] * ab_z + g_x_y_z_xxxzz[k];

                g_x_y_zz_xxyy[k] = -g_x_y_z_xxyy[k] * ab_z + g_x_y_z_xxyyz[k];

                g_x_y_zz_xxyz[k] = -g_x_y_z_xxyz[k] * ab_z + g_x_y_z_xxyzz[k];

                g_x_y_zz_xxzz[k] = -g_x_y_z_xxzz[k] * ab_z + g_x_y_z_xxzzz[k];

                g_x_y_zz_xyyy[k] = -g_x_y_z_xyyy[k] * ab_z + g_x_y_z_xyyyz[k];

                g_x_y_zz_xyyz[k] = -g_x_y_z_xyyz[k] * ab_z + g_x_y_z_xyyzz[k];

                g_x_y_zz_xyzz[k] = -g_x_y_z_xyzz[k] * ab_z + g_x_y_z_xyzzz[k];

                g_x_y_zz_xzzz[k] = -g_x_y_z_xzzz[k] * ab_z + g_x_y_z_xzzzz[k];

                g_x_y_zz_yyyy[k] = -g_x_y_z_yyyy[k] * ab_z + g_x_y_z_yyyyz[k];

                g_x_y_zz_yyyz[k] = -g_x_y_z_yyyz[k] * ab_z + g_x_y_z_yyyzz[k];

                g_x_y_zz_yyzz[k] = -g_x_y_z_yyzz[k] * ab_z + g_x_y_z_yyzzz[k];

                g_x_y_zz_yzzz[k] = -g_x_y_z_yzzz[k] * ab_z + g_x_y_z_yzzzz[k];

                g_x_y_zz_zzzz[k] = -g_x_y_z_zzzz[k] * ab_z + g_x_y_z_zzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxxx, g_0_z_x_xxxy, g_0_z_x_xxxz, g_0_z_x_xxyy, g_0_z_x_xxyz, g_0_z_x_xxzz, g_0_z_x_xyyy, g_0_z_x_xyyz, g_0_z_x_xyzz, g_0_z_x_xzzz, g_0_z_x_yyyy, g_0_z_x_yyyz, g_0_z_x_yyzz, g_0_z_x_yzzz, g_0_z_x_zzzz, g_x_z_x_xxxx, g_x_z_x_xxxxx, g_x_z_x_xxxxy, g_x_z_x_xxxxz, g_x_z_x_xxxy, g_x_z_x_xxxyy, g_x_z_x_xxxyz, g_x_z_x_xxxz, g_x_z_x_xxxzz, g_x_z_x_xxyy, g_x_z_x_xxyyy, g_x_z_x_xxyyz, g_x_z_x_xxyz, g_x_z_x_xxyzz, g_x_z_x_xxzz, g_x_z_x_xxzzz, g_x_z_x_xyyy, g_x_z_x_xyyyy, g_x_z_x_xyyyz, g_x_z_x_xyyz, g_x_z_x_xyyzz, g_x_z_x_xyzz, g_x_z_x_xyzzz, g_x_z_x_xzzz, g_x_z_x_xzzzz, g_x_z_x_yyyy, g_x_z_x_yyyz, g_x_z_x_yyzz, g_x_z_x_yzzz, g_x_z_x_zzzz, g_x_z_xx_xxxx, g_x_z_xx_xxxy, g_x_z_xx_xxxz, g_x_z_xx_xxyy, g_x_z_xx_xxyz, g_x_z_xx_xxzz, g_x_z_xx_xyyy, g_x_z_xx_xyyz, g_x_z_xx_xyzz, g_x_z_xx_xzzz, g_x_z_xx_yyyy, g_x_z_xx_yyyz, g_x_z_xx_yyzz, g_x_z_xx_yzzz, g_x_z_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_xxxx[k] = -g_0_z_x_xxxx[k] - g_x_z_x_xxxx[k] * ab_x + g_x_z_x_xxxxx[k];

                g_x_z_xx_xxxy[k] = -g_0_z_x_xxxy[k] - g_x_z_x_xxxy[k] * ab_x + g_x_z_x_xxxxy[k];

                g_x_z_xx_xxxz[k] = -g_0_z_x_xxxz[k] - g_x_z_x_xxxz[k] * ab_x + g_x_z_x_xxxxz[k];

                g_x_z_xx_xxyy[k] = -g_0_z_x_xxyy[k] - g_x_z_x_xxyy[k] * ab_x + g_x_z_x_xxxyy[k];

                g_x_z_xx_xxyz[k] = -g_0_z_x_xxyz[k] - g_x_z_x_xxyz[k] * ab_x + g_x_z_x_xxxyz[k];

                g_x_z_xx_xxzz[k] = -g_0_z_x_xxzz[k] - g_x_z_x_xxzz[k] * ab_x + g_x_z_x_xxxzz[k];

                g_x_z_xx_xyyy[k] = -g_0_z_x_xyyy[k] - g_x_z_x_xyyy[k] * ab_x + g_x_z_x_xxyyy[k];

                g_x_z_xx_xyyz[k] = -g_0_z_x_xyyz[k] - g_x_z_x_xyyz[k] * ab_x + g_x_z_x_xxyyz[k];

                g_x_z_xx_xyzz[k] = -g_0_z_x_xyzz[k] - g_x_z_x_xyzz[k] * ab_x + g_x_z_x_xxyzz[k];

                g_x_z_xx_xzzz[k] = -g_0_z_x_xzzz[k] - g_x_z_x_xzzz[k] * ab_x + g_x_z_x_xxzzz[k];

                g_x_z_xx_yyyy[k] = -g_0_z_x_yyyy[k] - g_x_z_x_yyyy[k] * ab_x + g_x_z_x_xyyyy[k];

                g_x_z_xx_yyyz[k] = -g_0_z_x_yyyz[k] - g_x_z_x_yyyz[k] * ab_x + g_x_z_x_xyyyz[k];

                g_x_z_xx_yyzz[k] = -g_0_z_x_yyzz[k] - g_x_z_x_yyzz[k] * ab_x + g_x_z_x_xyyzz[k];

                g_x_z_xx_yzzz[k] = -g_0_z_x_yzzz[k] - g_x_z_x_yzzz[k] * ab_x + g_x_z_x_xyzzz[k];

                g_x_z_xx_zzzz[k] = -g_0_z_x_zzzz[k] - g_x_z_x_zzzz[k] * ab_x + g_x_z_x_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_xxxx, g_x_z_x_xxxxy, g_x_z_x_xxxy, g_x_z_x_xxxyy, g_x_z_x_xxxyz, g_x_z_x_xxxz, g_x_z_x_xxyy, g_x_z_x_xxyyy, g_x_z_x_xxyyz, g_x_z_x_xxyz, g_x_z_x_xxyzz, g_x_z_x_xxzz, g_x_z_x_xyyy, g_x_z_x_xyyyy, g_x_z_x_xyyyz, g_x_z_x_xyyz, g_x_z_x_xyyzz, g_x_z_x_xyzz, g_x_z_x_xyzzz, g_x_z_x_xzzz, g_x_z_x_yyyy, g_x_z_x_yyyyy, g_x_z_x_yyyyz, g_x_z_x_yyyz, g_x_z_x_yyyzz, g_x_z_x_yyzz, g_x_z_x_yyzzz, g_x_z_x_yzzz, g_x_z_x_yzzzz, g_x_z_x_zzzz, g_x_z_xy_xxxx, g_x_z_xy_xxxy, g_x_z_xy_xxxz, g_x_z_xy_xxyy, g_x_z_xy_xxyz, g_x_z_xy_xxzz, g_x_z_xy_xyyy, g_x_z_xy_xyyz, g_x_z_xy_xyzz, g_x_z_xy_xzzz, g_x_z_xy_yyyy, g_x_z_xy_yyyz, g_x_z_xy_yyzz, g_x_z_xy_yzzz, g_x_z_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_xxxx[k] = -g_x_z_x_xxxx[k] * ab_y + g_x_z_x_xxxxy[k];

                g_x_z_xy_xxxy[k] = -g_x_z_x_xxxy[k] * ab_y + g_x_z_x_xxxyy[k];

                g_x_z_xy_xxxz[k] = -g_x_z_x_xxxz[k] * ab_y + g_x_z_x_xxxyz[k];

                g_x_z_xy_xxyy[k] = -g_x_z_x_xxyy[k] * ab_y + g_x_z_x_xxyyy[k];

                g_x_z_xy_xxyz[k] = -g_x_z_x_xxyz[k] * ab_y + g_x_z_x_xxyyz[k];

                g_x_z_xy_xxzz[k] = -g_x_z_x_xxzz[k] * ab_y + g_x_z_x_xxyzz[k];

                g_x_z_xy_xyyy[k] = -g_x_z_x_xyyy[k] * ab_y + g_x_z_x_xyyyy[k];

                g_x_z_xy_xyyz[k] = -g_x_z_x_xyyz[k] * ab_y + g_x_z_x_xyyyz[k];

                g_x_z_xy_xyzz[k] = -g_x_z_x_xyzz[k] * ab_y + g_x_z_x_xyyzz[k];

                g_x_z_xy_xzzz[k] = -g_x_z_x_xzzz[k] * ab_y + g_x_z_x_xyzzz[k];

                g_x_z_xy_yyyy[k] = -g_x_z_x_yyyy[k] * ab_y + g_x_z_x_yyyyy[k];

                g_x_z_xy_yyyz[k] = -g_x_z_x_yyyz[k] * ab_y + g_x_z_x_yyyyz[k];

                g_x_z_xy_yyzz[k] = -g_x_z_x_yyzz[k] * ab_y + g_x_z_x_yyyzz[k];

                g_x_z_xy_yzzz[k] = -g_x_z_x_yzzz[k] * ab_y + g_x_z_x_yyzzz[k];

                g_x_z_xy_zzzz[k] = -g_x_z_x_zzzz[k] * ab_y + g_x_z_x_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxx, g_0_z_z_xxxy, g_0_z_z_xxxz, g_0_z_z_xxyy, g_0_z_z_xxyz, g_0_z_z_xxzz, g_0_z_z_xyyy, g_0_z_z_xyyz, g_0_z_z_xyzz, g_0_z_z_xzzz, g_0_z_z_yyyy, g_0_z_z_yyyz, g_0_z_z_yyzz, g_0_z_z_yzzz, g_0_z_z_zzzz, g_x_z_xz_xxxx, g_x_z_xz_xxxy, g_x_z_xz_xxxz, g_x_z_xz_xxyy, g_x_z_xz_xxyz, g_x_z_xz_xxzz, g_x_z_xz_xyyy, g_x_z_xz_xyyz, g_x_z_xz_xyzz, g_x_z_xz_xzzz, g_x_z_xz_yyyy, g_x_z_xz_yyyz, g_x_z_xz_yyzz, g_x_z_xz_yzzz, g_x_z_xz_zzzz, g_x_z_z_xxxx, g_x_z_z_xxxxx, g_x_z_z_xxxxy, g_x_z_z_xxxxz, g_x_z_z_xxxy, g_x_z_z_xxxyy, g_x_z_z_xxxyz, g_x_z_z_xxxz, g_x_z_z_xxxzz, g_x_z_z_xxyy, g_x_z_z_xxyyy, g_x_z_z_xxyyz, g_x_z_z_xxyz, g_x_z_z_xxyzz, g_x_z_z_xxzz, g_x_z_z_xxzzz, g_x_z_z_xyyy, g_x_z_z_xyyyy, g_x_z_z_xyyyz, g_x_z_z_xyyz, g_x_z_z_xyyzz, g_x_z_z_xyzz, g_x_z_z_xyzzz, g_x_z_z_xzzz, g_x_z_z_xzzzz, g_x_z_z_yyyy, g_x_z_z_yyyz, g_x_z_z_yyzz, g_x_z_z_yzzz, g_x_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_xxxx[k] = -g_0_z_z_xxxx[k] - g_x_z_z_xxxx[k] * ab_x + g_x_z_z_xxxxx[k];

                g_x_z_xz_xxxy[k] = -g_0_z_z_xxxy[k] - g_x_z_z_xxxy[k] * ab_x + g_x_z_z_xxxxy[k];

                g_x_z_xz_xxxz[k] = -g_0_z_z_xxxz[k] - g_x_z_z_xxxz[k] * ab_x + g_x_z_z_xxxxz[k];

                g_x_z_xz_xxyy[k] = -g_0_z_z_xxyy[k] - g_x_z_z_xxyy[k] * ab_x + g_x_z_z_xxxyy[k];

                g_x_z_xz_xxyz[k] = -g_0_z_z_xxyz[k] - g_x_z_z_xxyz[k] * ab_x + g_x_z_z_xxxyz[k];

                g_x_z_xz_xxzz[k] = -g_0_z_z_xxzz[k] - g_x_z_z_xxzz[k] * ab_x + g_x_z_z_xxxzz[k];

                g_x_z_xz_xyyy[k] = -g_0_z_z_xyyy[k] - g_x_z_z_xyyy[k] * ab_x + g_x_z_z_xxyyy[k];

                g_x_z_xz_xyyz[k] = -g_0_z_z_xyyz[k] - g_x_z_z_xyyz[k] * ab_x + g_x_z_z_xxyyz[k];

                g_x_z_xz_xyzz[k] = -g_0_z_z_xyzz[k] - g_x_z_z_xyzz[k] * ab_x + g_x_z_z_xxyzz[k];

                g_x_z_xz_xzzz[k] = -g_0_z_z_xzzz[k] - g_x_z_z_xzzz[k] * ab_x + g_x_z_z_xxzzz[k];

                g_x_z_xz_yyyy[k] = -g_0_z_z_yyyy[k] - g_x_z_z_yyyy[k] * ab_x + g_x_z_z_xyyyy[k];

                g_x_z_xz_yyyz[k] = -g_0_z_z_yyyz[k] - g_x_z_z_yyyz[k] * ab_x + g_x_z_z_xyyyz[k];

                g_x_z_xz_yyzz[k] = -g_0_z_z_yyzz[k] - g_x_z_z_yyzz[k] * ab_x + g_x_z_z_xyyzz[k];

                g_x_z_xz_yzzz[k] = -g_0_z_z_yzzz[k] - g_x_z_z_yzzz[k] * ab_x + g_x_z_z_xyzzz[k];

                g_x_z_xz_zzzz[k] = -g_0_z_z_zzzz[k] - g_x_z_z_zzzz[k] * ab_x + g_x_z_z_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_xxxx, g_x_z_y_xxxxy, g_x_z_y_xxxy, g_x_z_y_xxxyy, g_x_z_y_xxxyz, g_x_z_y_xxxz, g_x_z_y_xxyy, g_x_z_y_xxyyy, g_x_z_y_xxyyz, g_x_z_y_xxyz, g_x_z_y_xxyzz, g_x_z_y_xxzz, g_x_z_y_xyyy, g_x_z_y_xyyyy, g_x_z_y_xyyyz, g_x_z_y_xyyz, g_x_z_y_xyyzz, g_x_z_y_xyzz, g_x_z_y_xyzzz, g_x_z_y_xzzz, g_x_z_y_yyyy, g_x_z_y_yyyyy, g_x_z_y_yyyyz, g_x_z_y_yyyz, g_x_z_y_yyyzz, g_x_z_y_yyzz, g_x_z_y_yyzzz, g_x_z_y_yzzz, g_x_z_y_yzzzz, g_x_z_y_zzzz, g_x_z_yy_xxxx, g_x_z_yy_xxxy, g_x_z_yy_xxxz, g_x_z_yy_xxyy, g_x_z_yy_xxyz, g_x_z_yy_xxzz, g_x_z_yy_xyyy, g_x_z_yy_xyyz, g_x_z_yy_xyzz, g_x_z_yy_xzzz, g_x_z_yy_yyyy, g_x_z_yy_yyyz, g_x_z_yy_yyzz, g_x_z_yy_yzzz, g_x_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_xxxx[k] = -g_x_z_y_xxxx[k] * ab_y + g_x_z_y_xxxxy[k];

                g_x_z_yy_xxxy[k] = -g_x_z_y_xxxy[k] * ab_y + g_x_z_y_xxxyy[k];

                g_x_z_yy_xxxz[k] = -g_x_z_y_xxxz[k] * ab_y + g_x_z_y_xxxyz[k];

                g_x_z_yy_xxyy[k] = -g_x_z_y_xxyy[k] * ab_y + g_x_z_y_xxyyy[k];

                g_x_z_yy_xxyz[k] = -g_x_z_y_xxyz[k] * ab_y + g_x_z_y_xxyyz[k];

                g_x_z_yy_xxzz[k] = -g_x_z_y_xxzz[k] * ab_y + g_x_z_y_xxyzz[k];

                g_x_z_yy_xyyy[k] = -g_x_z_y_xyyy[k] * ab_y + g_x_z_y_xyyyy[k];

                g_x_z_yy_xyyz[k] = -g_x_z_y_xyyz[k] * ab_y + g_x_z_y_xyyyz[k];

                g_x_z_yy_xyzz[k] = -g_x_z_y_xyzz[k] * ab_y + g_x_z_y_xyyzz[k];

                g_x_z_yy_xzzz[k] = -g_x_z_y_xzzz[k] * ab_y + g_x_z_y_xyzzz[k];

                g_x_z_yy_yyyy[k] = -g_x_z_y_yyyy[k] * ab_y + g_x_z_y_yyyyy[k];

                g_x_z_yy_yyyz[k] = -g_x_z_y_yyyz[k] * ab_y + g_x_z_y_yyyyz[k];

                g_x_z_yy_yyzz[k] = -g_x_z_y_yyzz[k] * ab_y + g_x_z_y_yyyzz[k];

                g_x_z_yy_yzzz[k] = -g_x_z_y_yzzz[k] * ab_y + g_x_z_y_yyzzz[k];

                g_x_z_yy_zzzz[k] = -g_x_z_y_zzzz[k] * ab_y + g_x_z_y_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_xxxx, g_x_z_yz_xxxy, g_x_z_yz_xxxz, g_x_z_yz_xxyy, g_x_z_yz_xxyz, g_x_z_yz_xxzz, g_x_z_yz_xyyy, g_x_z_yz_xyyz, g_x_z_yz_xyzz, g_x_z_yz_xzzz, g_x_z_yz_yyyy, g_x_z_yz_yyyz, g_x_z_yz_yyzz, g_x_z_yz_yzzz, g_x_z_yz_zzzz, g_x_z_z_xxxx, g_x_z_z_xxxxy, g_x_z_z_xxxy, g_x_z_z_xxxyy, g_x_z_z_xxxyz, g_x_z_z_xxxz, g_x_z_z_xxyy, g_x_z_z_xxyyy, g_x_z_z_xxyyz, g_x_z_z_xxyz, g_x_z_z_xxyzz, g_x_z_z_xxzz, g_x_z_z_xyyy, g_x_z_z_xyyyy, g_x_z_z_xyyyz, g_x_z_z_xyyz, g_x_z_z_xyyzz, g_x_z_z_xyzz, g_x_z_z_xyzzz, g_x_z_z_xzzz, g_x_z_z_yyyy, g_x_z_z_yyyyy, g_x_z_z_yyyyz, g_x_z_z_yyyz, g_x_z_z_yyyzz, g_x_z_z_yyzz, g_x_z_z_yyzzz, g_x_z_z_yzzz, g_x_z_z_yzzzz, g_x_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_xxxx[k] = -g_x_z_z_xxxx[k] * ab_y + g_x_z_z_xxxxy[k];

                g_x_z_yz_xxxy[k] = -g_x_z_z_xxxy[k] * ab_y + g_x_z_z_xxxyy[k];

                g_x_z_yz_xxxz[k] = -g_x_z_z_xxxz[k] * ab_y + g_x_z_z_xxxyz[k];

                g_x_z_yz_xxyy[k] = -g_x_z_z_xxyy[k] * ab_y + g_x_z_z_xxyyy[k];

                g_x_z_yz_xxyz[k] = -g_x_z_z_xxyz[k] * ab_y + g_x_z_z_xxyyz[k];

                g_x_z_yz_xxzz[k] = -g_x_z_z_xxzz[k] * ab_y + g_x_z_z_xxyzz[k];

                g_x_z_yz_xyyy[k] = -g_x_z_z_xyyy[k] * ab_y + g_x_z_z_xyyyy[k];

                g_x_z_yz_xyyz[k] = -g_x_z_z_xyyz[k] * ab_y + g_x_z_z_xyyyz[k];

                g_x_z_yz_xyzz[k] = -g_x_z_z_xyzz[k] * ab_y + g_x_z_z_xyyzz[k];

                g_x_z_yz_xzzz[k] = -g_x_z_z_xzzz[k] * ab_y + g_x_z_z_xyzzz[k];

                g_x_z_yz_yyyy[k] = -g_x_z_z_yyyy[k] * ab_y + g_x_z_z_yyyyy[k];

                g_x_z_yz_yyyz[k] = -g_x_z_z_yyyz[k] * ab_y + g_x_z_z_yyyyz[k];

                g_x_z_yz_yyzz[k] = -g_x_z_z_yyzz[k] * ab_y + g_x_z_z_yyyzz[k];

                g_x_z_yz_yzzz[k] = -g_x_z_z_yzzz[k] * ab_y + g_x_z_z_yyzzz[k];

                g_x_z_yz_zzzz[k] = -g_x_z_z_zzzz[k] * ab_y + g_x_z_z_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxx, g_x_0_z_xxxy, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyzz, g_x_0_z_yzzz, g_x_0_z_zzzz, g_x_z_z_xxxx, g_x_z_z_xxxxz, g_x_z_z_xxxy, g_x_z_z_xxxyz, g_x_z_z_xxxz, g_x_z_z_xxxzz, g_x_z_z_xxyy, g_x_z_z_xxyyz, g_x_z_z_xxyz, g_x_z_z_xxyzz, g_x_z_z_xxzz, g_x_z_z_xxzzz, g_x_z_z_xyyy, g_x_z_z_xyyyz, g_x_z_z_xyyz, g_x_z_z_xyyzz, g_x_z_z_xyzz, g_x_z_z_xyzzz, g_x_z_z_xzzz, g_x_z_z_xzzzz, g_x_z_z_yyyy, g_x_z_z_yyyyz, g_x_z_z_yyyz, g_x_z_z_yyyzz, g_x_z_z_yyzz, g_x_z_z_yyzzz, g_x_z_z_yzzz, g_x_z_z_yzzzz, g_x_z_z_zzzz, g_x_z_z_zzzzz, g_x_z_zz_xxxx, g_x_z_zz_xxxy, g_x_z_zz_xxxz, g_x_z_zz_xxyy, g_x_z_zz_xxyz, g_x_z_zz_xxzz, g_x_z_zz_xyyy, g_x_z_zz_xyyz, g_x_z_zz_xyzz, g_x_z_zz_xzzz, g_x_z_zz_yyyy, g_x_z_zz_yyyz, g_x_z_zz_yyzz, g_x_z_zz_yzzz, g_x_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_xxxx[k] = g_x_0_z_xxxx[k] - g_x_z_z_xxxx[k] * ab_z + g_x_z_z_xxxxz[k];

                g_x_z_zz_xxxy[k] = g_x_0_z_xxxy[k] - g_x_z_z_xxxy[k] * ab_z + g_x_z_z_xxxyz[k];

                g_x_z_zz_xxxz[k] = g_x_0_z_xxxz[k] - g_x_z_z_xxxz[k] * ab_z + g_x_z_z_xxxzz[k];

                g_x_z_zz_xxyy[k] = g_x_0_z_xxyy[k] - g_x_z_z_xxyy[k] * ab_z + g_x_z_z_xxyyz[k];

                g_x_z_zz_xxyz[k] = g_x_0_z_xxyz[k] - g_x_z_z_xxyz[k] * ab_z + g_x_z_z_xxyzz[k];

                g_x_z_zz_xxzz[k] = g_x_0_z_xxzz[k] - g_x_z_z_xxzz[k] * ab_z + g_x_z_z_xxzzz[k];

                g_x_z_zz_xyyy[k] = g_x_0_z_xyyy[k] - g_x_z_z_xyyy[k] * ab_z + g_x_z_z_xyyyz[k];

                g_x_z_zz_xyyz[k] = g_x_0_z_xyyz[k] - g_x_z_z_xyyz[k] * ab_z + g_x_z_z_xyyzz[k];

                g_x_z_zz_xyzz[k] = g_x_0_z_xyzz[k] - g_x_z_z_xyzz[k] * ab_z + g_x_z_z_xyzzz[k];

                g_x_z_zz_xzzz[k] = g_x_0_z_xzzz[k] - g_x_z_z_xzzz[k] * ab_z + g_x_z_z_xzzzz[k];

                g_x_z_zz_yyyy[k] = g_x_0_z_yyyy[k] - g_x_z_z_yyyy[k] * ab_z + g_x_z_z_yyyyz[k];

                g_x_z_zz_yyyz[k] = g_x_0_z_yyyz[k] - g_x_z_z_yyyz[k] * ab_z + g_x_z_z_yyyzz[k];

                g_x_z_zz_yyzz[k] = g_x_0_z_yyzz[k] - g_x_z_z_yyzz[k] * ab_z + g_x_z_z_yyzzz[k];

                g_x_z_zz_yzzz[k] = g_x_0_z_yzzz[k] - g_x_z_z_yzzz[k] * ab_z + g_x_z_z_yzzzz[k];

                g_x_z_zz_zzzz[k] = g_x_0_z_zzzz[k] - g_x_z_z_zzzz[k] * ab_z + g_x_z_z_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxzz, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyzz, g_y_0_x_xzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz, g_y_x_x_xxxx, g_y_x_x_xxxxx, g_y_x_x_xxxxy, g_y_x_x_xxxxz, g_y_x_x_xxxy, g_y_x_x_xxxyy, g_y_x_x_xxxyz, g_y_x_x_xxxz, g_y_x_x_xxxzz, g_y_x_x_xxyy, g_y_x_x_xxyyy, g_y_x_x_xxyyz, g_y_x_x_xxyz, g_y_x_x_xxyzz, g_y_x_x_xxzz, g_y_x_x_xxzzz, g_y_x_x_xyyy, g_y_x_x_xyyyy, g_y_x_x_xyyyz, g_y_x_x_xyyz, g_y_x_x_xyyzz, g_y_x_x_xyzz, g_y_x_x_xyzzz, g_y_x_x_xzzz, g_y_x_x_xzzzz, g_y_x_x_yyyy, g_y_x_x_yyyz, g_y_x_x_yyzz, g_y_x_x_yzzz, g_y_x_x_zzzz, g_y_x_xx_xxxx, g_y_x_xx_xxxy, g_y_x_xx_xxxz, g_y_x_xx_xxyy, g_y_x_xx_xxyz, g_y_x_xx_xxzz, g_y_x_xx_xyyy, g_y_x_xx_xyyz, g_y_x_xx_xyzz, g_y_x_xx_xzzz, g_y_x_xx_yyyy, g_y_x_xx_yyyz, g_y_x_xx_yyzz, g_y_x_xx_yzzz, g_y_x_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_xxxx[k] = g_y_0_x_xxxx[k] - g_y_x_x_xxxx[k] * ab_x + g_y_x_x_xxxxx[k];

                g_y_x_xx_xxxy[k] = g_y_0_x_xxxy[k] - g_y_x_x_xxxy[k] * ab_x + g_y_x_x_xxxxy[k];

                g_y_x_xx_xxxz[k] = g_y_0_x_xxxz[k] - g_y_x_x_xxxz[k] * ab_x + g_y_x_x_xxxxz[k];

                g_y_x_xx_xxyy[k] = g_y_0_x_xxyy[k] - g_y_x_x_xxyy[k] * ab_x + g_y_x_x_xxxyy[k];

                g_y_x_xx_xxyz[k] = g_y_0_x_xxyz[k] - g_y_x_x_xxyz[k] * ab_x + g_y_x_x_xxxyz[k];

                g_y_x_xx_xxzz[k] = g_y_0_x_xxzz[k] - g_y_x_x_xxzz[k] * ab_x + g_y_x_x_xxxzz[k];

                g_y_x_xx_xyyy[k] = g_y_0_x_xyyy[k] - g_y_x_x_xyyy[k] * ab_x + g_y_x_x_xxyyy[k];

                g_y_x_xx_xyyz[k] = g_y_0_x_xyyz[k] - g_y_x_x_xyyz[k] * ab_x + g_y_x_x_xxyyz[k];

                g_y_x_xx_xyzz[k] = g_y_0_x_xyzz[k] - g_y_x_x_xyzz[k] * ab_x + g_y_x_x_xxyzz[k];

                g_y_x_xx_xzzz[k] = g_y_0_x_xzzz[k] - g_y_x_x_xzzz[k] * ab_x + g_y_x_x_xxzzz[k];

                g_y_x_xx_yyyy[k] = g_y_0_x_yyyy[k] - g_y_x_x_yyyy[k] * ab_x + g_y_x_x_xyyyy[k];

                g_y_x_xx_yyyz[k] = g_y_0_x_yyyz[k] - g_y_x_x_yyyz[k] * ab_x + g_y_x_x_xyyyz[k];

                g_y_x_xx_yyzz[k] = g_y_0_x_yyzz[k] - g_y_x_x_yyzz[k] * ab_x + g_y_x_x_xyyzz[k];

                g_y_x_xx_yzzz[k] = g_y_0_x_yzzz[k] - g_y_x_x_yzzz[k] * ab_x + g_y_x_x_xyzzz[k];

                g_y_x_xx_zzzz[k] = g_y_0_x_zzzz[k] - g_y_x_x_zzzz[k] * ab_x + g_y_x_x_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz, g_y_x_xy_xxxx, g_y_x_xy_xxxy, g_y_x_xy_xxxz, g_y_x_xy_xxyy, g_y_x_xy_xxyz, g_y_x_xy_xxzz, g_y_x_xy_xyyy, g_y_x_xy_xyyz, g_y_x_xy_xyzz, g_y_x_xy_xzzz, g_y_x_xy_yyyy, g_y_x_xy_yyyz, g_y_x_xy_yyzz, g_y_x_xy_yzzz, g_y_x_xy_zzzz, g_y_x_y_xxxx, g_y_x_y_xxxxx, g_y_x_y_xxxxy, g_y_x_y_xxxxz, g_y_x_y_xxxy, g_y_x_y_xxxyy, g_y_x_y_xxxyz, g_y_x_y_xxxz, g_y_x_y_xxxzz, g_y_x_y_xxyy, g_y_x_y_xxyyy, g_y_x_y_xxyyz, g_y_x_y_xxyz, g_y_x_y_xxyzz, g_y_x_y_xxzz, g_y_x_y_xxzzz, g_y_x_y_xyyy, g_y_x_y_xyyyy, g_y_x_y_xyyyz, g_y_x_y_xyyz, g_y_x_y_xyyzz, g_y_x_y_xyzz, g_y_x_y_xyzzz, g_y_x_y_xzzz, g_y_x_y_xzzzz, g_y_x_y_yyyy, g_y_x_y_yyyz, g_y_x_y_yyzz, g_y_x_y_yzzz, g_y_x_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_xxxx[k] = g_y_0_y_xxxx[k] - g_y_x_y_xxxx[k] * ab_x + g_y_x_y_xxxxx[k];

                g_y_x_xy_xxxy[k] = g_y_0_y_xxxy[k] - g_y_x_y_xxxy[k] * ab_x + g_y_x_y_xxxxy[k];

                g_y_x_xy_xxxz[k] = g_y_0_y_xxxz[k] - g_y_x_y_xxxz[k] * ab_x + g_y_x_y_xxxxz[k];

                g_y_x_xy_xxyy[k] = g_y_0_y_xxyy[k] - g_y_x_y_xxyy[k] * ab_x + g_y_x_y_xxxyy[k];

                g_y_x_xy_xxyz[k] = g_y_0_y_xxyz[k] - g_y_x_y_xxyz[k] * ab_x + g_y_x_y_xxxyz[k];

                g_y_x_xy_xxzz[k] = g_y_0_y_xxzz[k] - g_y_x_y_xxzz[k] * ab_x + g_y_x_y_xxxzz[k];

                g_y_x_xy_xyyy[k] = g_y_0_y_xyyy[k] - g_y_x_y_xyyy[k] * ab_x + g_y_x_y_xxyyy[k];

                g_y_x_xy_xyyz[k] = g_y_0_y_xyyz[k] - g_y_x_y_xyyz[k] * ab_x + g_y_x_y_xxyyz[k];

                g_y_x_xy_xyzz[k] = g_y_0_y_xyzz[k] - g_y_x_y_xyzz[k] * ab_x + g_y_x_y_xxyzz[k];

                g_y_x_xy_xzzz[k] = g_y_0_y_xzzz[k] - g_y_x_y_xzzz[k] * ab_x + g_y_x_y_xxzzz[k];

                g_y_x_xy_yyyy[k] = g_y_0_y_yyyy[k] - g_y_x_y_yyyy[k] * ab_x + g_y_x_y_xyyyy[k];

                g_y_x_xy_yyyz[k] = g_y_0_y_yyyz[k] - g_y_x_y_yyyz[k] * ab_x + g_y_x_y_xyyyz[k];

                g_y_x_xy_yyzz[k] = g_y_0_y_yyzz[k] - g_y_x_y_yyzz[k] * ab_x + g_y_x_y_xyyzz[k];

                g_y_x_xy_yzzz[k] = g_y_0_y_yzzz[k] - g_y_x_y_yzzz[k] * ab_x + g_y_x_y_xyzzz[k];

                g_y_x_xy_zzzz[k] = g_y_0_y_zzzz[k] - g_y_x_y_zzzz[k] * ab_x + g_y_x_y_xzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_xxxx, g_y_x_x_xxxxz, g_y_x_x_xxxy, g_y_x_x_xxxyz, g_y_x_x_xxxz, g_y_x_x_xxxzz, g_y_x_x_xxyy, g_y_x_x_xxyyz, g_y_x_x_xxyz, g_y_x_x_xxyzz, g_y_x_x_xxzz, g_y_x_x_xxzzz, g_y_x_x_xyyy, g_y_x_x_xyyyz, g_y_x_x_xyyz, g_y_x_x_xyyzz, g_y_x_x_xyzz, g_y_x_x_xyzzz, g_y_x_x_xzzz, g_y_x_x_xzzzz, g_y_x_x_yyyy, g_y_x_x_yyyyz, g_y_x_x_yyyz, g_y_x_x_yyyzz, g_y_x_x_yyzz, g_y_x_x_yyzzz, g_y_x_x_yzzz, g_y_x_x_yzzzz, g_y_x_x_zzzz, g_y_x_x_zzzzz, g_y_x_xz_xxxx, g_y_x_xz_xxxy, g_y_x_xz_xxxz, g_y_x_xz_xxyy, g_y_x_xz_xxyz, g_y_x_xz_xxzz, g_y_x_xz_xyyy, g_y_x_xz_xyyz, g_y_x_xz_xyzz, g_y_x_xz_xzzz, g_y_x_xz_yyyy, g_y_x_xz_yyyz, g_y_x_xz_yyzz, g_y_x_xz_yzzz, g_y_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_xxxx[k] = -g_y_x_x_xxxx[k] * ab_z + g_y_x_x_xxxxz[k];

                g_y_x_xz_xxxy[k] = -g_y_x_x_xxxy[k] * ab_z + g_y_x_x_xxxyz[k];

                g_y_x_xz_xxxz[k] = -g_y_x_x_xxxz[k] * ab_z + g_y_x_x_xxxzz[k];

                g_y_x_xz_xxyy[k] = -g_y_x_x_xxyy[k] * ab_z + g_y_x_x_xxyyz[k];

                g_y_x_xz_xxyz[k] = -g_y_x_x_xxyz[k] * ab_z + g_y_x_x_xxyzz[k];

                g_y_x_xz_xxzz[k] = -g_y_x_x_xxzz[k] * ab_z + g_y_x_x_xxzzz[k];

                g_y_x_xz_xyyy[k] = -g_y_x_x_xyyy[k] * ab_z + g_y_x_x_xyyyz[k];

                g_y_x_xz_xyyz[k] = -g_y_x_x_xyyz[k] * ab_z + g_y_x_x_xyyzz[k];

                g_y_x_xz_xyzz[k] = -g_y_x_x_xyzz[k] * ab_z + g_y_x_x_xyzzz[k];

                g_y_x_xz_xzzz[k] = -g_y_x_x_xzzz[k] * ab_z + g_y_x_x_xzzzz[k];

                g_y_x_xz_yyyy[k] = -g_y_x_x_yyyy[k] * ab_z + g_y_x_x_yyyyz[k];

                g_y_x_xz_yyyz[k] = -g_y_x_x_yyyz[k] * ab_z + g_y_x_x_yyyzz[k];

                g_y_x_xz_yyzz[k] = -g_y_x_x_yyzz[k] * ab_z + g_y_x_x_yyzzz[k];

                g_y_x_xz_yzzz[k] = -g_y_x_x_yzzz[k] * ab_z + g_y_x_x_yzzzz[k];

                g_y_x_xz_zzzz[k] = -g_y_x_x_zzzz[k] * ab_z + g_y_x_x_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxxx, g_0_x_y_xxxy, g_0_x_y_xxxz, g_0_x_y_xxyy, g_0_x_y_xxyz, g_0_x_y_xxzz, g_0_x_y_xyyy, g_0_x_y_xyyz, g_0_x_y_xyzz, g_0_x_y_xzzz, g_0_x_y_yyyy, g_0_x_y_yyyz, g_0_x_y_yyzz, g_0_x_y_yzzz, g_0_x_y_zzzz, g_y_x_y_xxxx, g_y_x_y_xxxxy, g_y_x_y_xxxy, g_y_x_y_xxxyy, g_y_x_y_xxxyz, g_y_x_y_xxxz, g_y_x_y_xxyy, g_y_x_y_xxyyy, g_y_x_y_xxyyz, g_y_x_y_xxyz, g_y_x_y_xxyzz, g_y_x_y_xxzz, g_y_x_y_xyyy, g_y_x_y_xyyyy, g_y_x_y_xyyyz, g_y_x_y_xyyz, g_y_x_y_xyyzz, g_y_x_y_xyzz, g_y_x_y_xyzzz, g_y_x_y_xzzz, g_y_x_y_yyyy, g_y_x_y_yyyyy, g_y_x_y_yyyyz, g_y_x_y_yyyz, g_y_x_y_yyyzz, g_y_x_y_yyzz, g_y_x_y_yyzzz, g_y_x_y_yzzz, g_y_x_y_yzzzz, g_y_x_y_zzzz, g_y_x_yy_xxxx, g_y_x_yy_xxxy, g_y_x_yy_xxxz, g_y_x_yy_xxyy, g_y_x_yy_xxyz, g_y_x_yy_xxzz, g_y_x_yy_xyyy, g_y_x_yy_xyyz, g_y_x_yy_xyzz, g_y_x_yy_xzzz, g_y_x_yy_yyyy, g_y_x_yy_yyyz, g_y_x_yy_yyzz, g_y_x_yy_yzzz, g_y_x_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_xxxx[k] = -g_0_x_y_xxxx[k] - g_y_x_y_xxxx[k] * ab_y + g_y_x_y_xxxxy[k];

                g_y_x_yy_xxxy[k] = -g_0_x_y_xxxy[k] - g_y_x_y_xxxy[k] * ab_y + g_y_x_y_xxxyy[k];

                g_y_x_yy_xxxz[k] = -g_0_x_y_xxxz[k] - g_y_x_y_xxxz[k] * ab_y + g_y_x_y_xxxyz[k];

                g_y_x_yy_xxyy[k] = -g_0_x_y_xxyy[k] - g_y_x_y_xxyy[k] * ab_y + g_y_x_y_xxyyy[k];

                g_y_x_yy_xxyz[k] = -g_0_x_y_xxyz[k] - g_y_x_y_xxyz[k] * ab_y + g_y_x_y_xxyyz[k];

                g_y_x_yy_xxzz[k] = -g_0_x_y_xxzz[k] - g_y_x_y_xxzz[k] * ab_y + g_y_x_y_xxyzz[k];

                g_y_x_yy_xyyy[k] = -g_0_x_y_xyyy[k] - g_y_x_y_xyyy[k] * ab_y + g_y_x_y_xyyyy[k];

                g_y_x_yy_xyyz[k] = -g_0_x_y_xyyz[k] - g_y_x_y_xyyz[k] * ab_y + g_y_x_y_xyyyz[k];

                g_y_x_yy_xyzz[k] = -g_0_x_y_xyzz[k] - g_y_x_y_xyzz[k] * ab_y + g_y_x_y_xyyzz[k];

                g_y_x_yy_xzzz[k] = -g_0_x_y_xzzz[k] - g_y_x_y_xzzz[k] * ab_y + g_y_x_y_xyzzz[k];

                g_y_x_yy_yyyy[k] = -g_0_x_y_yyyy[k] - g_y_x_y_yyyy[k] * ab_y + g_y_x_y_yyyyy[k];

                g_y_x_yy_yyyz[k] = -g_0_x_y_yyyz[k] - g_y_x_y_yyyz[k] * ab_y + g_y_x_y_yyyyz[k];

                g_y_x_yy_yyzz[k] = -g_0_x_y_yyzz[k] - g_y_x_y_yyzz[k] * ab_y + g_y_x_y_yyyzz[k];

                g_y_x_yy_yzzz[k] = -g_0_x_y_yzzz[k] - g_y_x_y_yzzz[k] * ab_y + g_y_x_y_yyzzz[k];

                g_y_x_yy_zzzz[k] = -g_0_x_y_zzzz[k] - g_y_x_y_zzzz[k] * ab_y + g_y_x_y_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_xxxx, g_y_x_y_xxxxz, g_y_x_y_xxxy, g_y_x_y_xxxyz, g_y_x_y_xxxz, g_y_x_y_xxxzz, g_y_x_y_xxyy, g_y_x_y_xxyyz, g_y_x_y_xxyz, g_y_x_y_xxyzz, g_y_x_y_xxzz, g_y_x_y_xxzzz, g_y_x_y_xyyy, g_y_x_y_xyyyz, g_y_x_y_xyyz, g_y_x_y_xyyzz, g_y_x_y_xyzz, g_y_x_y_xyzzz, g_y_x_y_xzzz, g_y_x_y_xzzzz, g_y_x_y_yyyy, g_y_x_y_yyyyz, g_y_x_y_yyyz, g_y_x_y_yyyzz, g_y_x_y_yyzz, g_y_x_y_yyzzz, g_y_x_y_yzzz, g_y_x_y_yzzzz, g_y_x_y_zzzz, g_y_x_y_zzzzz, g_y_x_yz_xxxx, g_y_x_yz_xxxy, g_y_x_yz_xxxz, g_y_x_yz_xxyy, g_y_x_yz_xxyz, g_y_x_yz_xxzz, g_y_x_yz_xyyy, g_y_x_yz_xyyz, g_y_x_yz_xyzz, g_y_x_yz_xzzz, g_y_x_yz_yyyy, g_y_x_yz_yyyz, g_y_x_yz_yyzz, g_y_x_yz_yzzz, g_y_x_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_xxxx[k] = -g_y_x_y_xxxx[k] * ab_z + g_y_x_y_xxxxz[k];

                g_y_x_yz_xxxy[k] = -g_y_x_y_xxxy[k] * ab_z + g_y_x_y_xxxyz[k];

                g_y_x_yz_xxxz[k] = -g_y_x_y_xxxz[k] * ab_z + g_y_x_y_xxxzz[k];

                g_y_x_yz_xxyy[k] = -g_y_x_y_xxyy[k] * ab_z + g_y_x_y_xxyyz[k];

                g_y_x_yz_xxyz[k] = -g_y_x_y_xxyz[k] * ab_z + g_y_x_y_xxyzz[k];

                g_y_x_yz_xxzz[k] = -g_y_x_y_xxzz[k] * ab_z + g_y_x_y_xxzzz[k];

                g_y_x_yz_xyyy[k] = -g_y_x_y_xyyy[k] * ab_z + g_y_x_y_xyyyz[k];

                g_y_x_yz_xyyz[k] = -g_y_x_y_xyyz[k] * ab_z + g_y_x_y_xyyzz[k];

                g_y_x_yz_xyzz[k] = -g_y_x_y_xyzz[k] * ab_z + g_y_x_y_xyzzz[k];

                g_y_x_yz_xzzz[k] = -g_y_x_y_xzzz[k] * ab_z + g_y_x_y_xzzzz[k];

                g_y_x_yz_yyyy[k] = -g_y_x_y_yyyy[k] * ab_z + g_y_x_y_yyyyz[k];

                g_y_x_yz_yyyz[k] = -g_y_x_y_yyyz[k] * ab_z + g_y_x_y_yyyzz[k];

                g_y_x_yz_yyzz[k] = -g_y_x_y_yyzz[k] * ab_z + g_y_x_y_yyzzz[k];

                g_y_x_yz_yzzz[k] = -g_y_x_y_yzzz[k] * ab_z + g_y_x_y_yzzzz[k];

                g_y_x_yz_zzzz[k] = -g_y_x_y_zzzz[k] * ab_z + g_y_x_y_zzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_xxxx, g_y_x_z_xxxxz, g_y_x_z_xxxy, g_y_x_z_xxxyz, g_y_x_z_xxxz, g_y_x_z_xxxzz, g_y_x_z_xxyy, g_y_x_z_xxyyz, g_y_x_z_xxyz, g_y_x_z_xxyzz, g_y_x_z_xxzz, g_y_x_z_xxzzz, g_y_x_z_xyyy, g_y_x_z_xyyyz, g_y_x_z_xyyz, g_y_x_z_xyyzz, g_y_x_z_xyzz, g_y_x_z_xyzzz, g_y_x_z_xzzz, g_y_x_z_xzzzz, g_y_x_z_yyyy, g_y_x_z_yyyyz, g_y_x_z_yyyz, g_y_x_z_yyyzz, g_y_x_z_yyzz, g_y_x_z_yyzzz, g_y_x_z_yzzz, g_y_x_z_yzzzz, g_y_x_z_zzzz, g_y_x_z_zzzzz, g_y_x_zz_xxxx, g_y_x_zz_xxxy, g_y_x_zz_xxxz, g_y_x_zz_xxyy, g_y_x_zz_xxyz, g_y_x_zz_xxzz, g_y_x_zz_xyyy, g_y_x_zz_xyyz, g_y_x_zz_xyzz, g_y_x_zz_xzzz, g_y_x_zz_yyyy, g_y_x_zz_yyyz, g_y_x_zz_yyzz, g_y_x_zz_yzzz, g_y_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_xxxx[k] = -g_y_x_z_xxxx[k] * ab_z + g_y_x_z_xxxxz[k];

                g_y_x_zz_xxxy[k] = -g_y_x_z_xxxy[k] * ab_z + g_y_x_z_xxxyz[k];

                g_y_x_zz_xxxz[k] = -g_y_x_z_xxxz[k] * ab_z + g_y_x_z_xxxzz[k];

                g_y_x_zz_xxyy[k] = -g_y_x_z_xxyy[k] * ab_z + g_y_x_z_xxyyz[k];

                g_y_x_zz_xxyz[k] = -g_y_x_z_xxyz[k] * ab_z + g_y_x_z_xxyzz[k];

                g_y_x_zz_xxzz[k] = -g_y_x_z_xxzz[k] * ab_z + g_y_x_z_xxzzz[k];

                g_y_x_zz_xyyy[k] = -g_y_x_z_xyyy[k] * ab_z + g_y_x_z_xyyyz[k];

                g_y_x_zz_xyyz[k] = -g_y_x_z_xyyz[k] * ab_z + g_y_x_z_xyyzz[k];

                g_y_x_zz_xyzz[k] = -g_y_x_z_xyzz[k] * ab_z + g_y_x_z_xyzzz[k];

                g_y_x_zz_xzzz[k] = -g_y_x_z_xzzz[k] * ab_z + g_y_x_z_xzzzz[k];

                g_y_x_zz_yyyy[k] = -g_y_x_z_yyyy[k] * ab_z + g_y_x_z_yyyyz[k];

                g_y_x_zz_yyyz[k] = -g_y_x_z_yyyz[k] * ab_z + g_y_x_z_yyyzz[k];

                g_y_x_zz_yyzz[k] = -g_y_x_z_yyzz[k] * ab_z + g_y_x_z_yyzzz[k];

                g_y_x_zz_yzzz[k] = -g_y_x_z_yzzz[k] * ab_z + g_y_x_z_yzzzz[k];

                g_y_x_zz_zzzz[k] = -g_y_x_z_zzzz[k] * ab_z + g_y_x_z_zzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_xxxx, g_y_y_x_xxxxx, g_y_y_x_xxxxy, g_y_y_x_xxxxz, g_y_y_x_xxxy, g_y_y_x_xxxyy, g_y_y_x_xxxyz, g_y_y_x_xxxz, g_y_y_x_xxxzz, g_y_y_x_xxyy, g_y_y_x_xxyyy, g_y_y_x_xxyyz, g_y_y_x_xxyz, g_y_y_x_xxyzz, g_y_y_x_xxzz, g_y_y_x_xxzzz, g_y_y_x_xyyy, g_y_y_x_xyyyy, g_y_y_x_xyyyz, g_y_y_x_xyyz, g_y_y_x_xyyzz, g_y_y_x_xyzz, g_y_y_x_xyzzz, g_y_y_x_xzzz, g_y_y_x_xzzzz, g_y_y_x_yyyy, g_y_y_x_yyyz, g_y_y_x_yyzz, g_y_y_x_yzzz, g_y_y_x_zzzz, g_y_y_xx_xxxx, g_y_y_xx_xxxy, g_y_y_xx_xxxz, g_y_y_xx_xxyy, g_y_y_xx_xxyz, g_y_y_xx_xxzz, g_y_y_xx_xyyy, g_y_y_xx_xyyz, g_y_y_xx_xyzz, g_y_y_xx_xzzz, g_y_y_xx_yyyy, g_y_y_xx_yyyz, g_y_y_xx_yyzz, g_y_y_xx_yzzz, g_y_y_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_xxxx[k] = -g_y_y_x_xxxx[k] * ab_x + g_y_y_x_xxxxx[k];

                g_y_y_xx_xxxy[k] = -g_y_y_x_xxxy[k] * ab_x + g_y_y_x_xxxxy[k];

                g_y_y_xx_xxxz[k] = -g_y_y_x_xxxz[k] * ab_x + g_y_y_x_xxxxz[k];

                g_y_y_xx_xxyy[k] = -g_y_y_x_xxyy[k] * ab_x + g_y_y_x_xxxyy[k];

                g_y_y_xx_xxyz[k] = -g_y_y_x_xxyz[k] * ab_x + g_y_y_x_xxxyz[k];

                g_y_y_xx_xxzz[k] = -g_y_y_x_xxzz[k] * ab_x + g_y_y_x_xxxzz[k];

                g_y_y_xx_xyyy[k] = -g_y_y_x_xyyy[k] * ab_x + g_y_y_x_xxyyy[k];

                g_y_y_xx_xyyz[k] = -g_y_y_x_xyyz[k] * ab_x + g_y_y_x_xxyyz[k];

                g_y_y_xx_xyzz[k] = -g_y_y_x_xyzz[k] * ab_x + g_y_y_x_xxyzz[k];

                g_y_y_xx_xzzz[k] = -g_y_y_x_xzzz[k] * ab_x + g_y_y_x_xxzzz[k];

                g_y_y_xx_yyyy[k] = -g_y_y_x_yyyy[k] * ab_x + g_y_y_x_xyyyy[k];

                g_y_y_xx_yyyz[k] = -g_y_y_x_yyyz[k] * ab_x + g_y_y_x_xyyyz[k];

                g_y_y_xx_yyzz[k] = -g_y_y_x_yyzz[k] * ab_x + g_y_y_x_xyyzz[k];

                g_y_y_xx_yzzz[k] = -g_y_y_x_yzzz[k] * ab_x + g_y_y_x_xyzzz[k];

                g_y_y_xx_zzzz[k] = -g_y_y_x_zzzz[k] * ab_x + g_y_y_x_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_xxxx, g_y_y_xy_xxxy, g_y_y_xy_xxxz, g_y_y_xy_xxyy, g_y_y_xy_xxyz, g_y_y_xy_xxzz, g_y_y_xy_xyyy, g_y_y_xy_xyyz, g_y_y_xy_xyzz, g_y_y_xy_xzzz, g_y_y_xy_yyyy, g_y_y_xy_yyyz, g_y_y_xy_yyzz, g_y_y_xy_yzzz, g_y_y_xy_zzzz, g_y_y_y_xxxx, g_y_y_y_xxxxx, g_y_y_y_xxxxy, g_y_y_y_xxxxz, g_y_y_y_xxxy, g_y_y_y_xxxyy, g_y_y_y_xxxyz, g_y_y_y_xxxz, g_y_y_y_xxxzz, g_y_y_y_xxyy, g_y_y_y_xxyyy, g_y_y_y_xxyyz, g_y_y_y_xxyz, g_y_y_y_xxyzz, g_y_y_y_xxzz, g_y_y_y_xxzzz, g_y_y_y_xyyy, g_y_y_y_xyyyy, g_y_y_y_xyyyz, g_y_y_y_xyyz, g_y_y_y_xyyzz, g_y_y_y_xyzz, g_y_y_y_xyzzz, g_y_y_y_xzzz, g_y_y_y_xzzzz, g_y_y_y_yyyy, g_y_y_y_yyyz, g_y_y_y_yyzz, g_y_y_y_yzzz, g_y_y_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_xxxx[k] = -g_y_y_y_xxxx[k] * ab_x + g_y_y_y_xxxxx[k];

                g_y_y_xy_xxxy[k] = -g_y_y_y_xxxy[k] * ab_x + g_y_y_y_xxxxy[k];

                g_y_y_xy_xxxz[k] = -g_y_y_y_xxxz[k] * ab_x + g_y_y_y_xxxxz[k];

                g_y_y_xy_xxyy[k] = -g_y_y_y_xxyy[k] * ab_x + g_y_y_y_xxxyy[k];

                g_y_y_xy_xxyz[k] = -g_y_y_y_xxyz[k] * ab_x + g_y_y_y_xxxyz[k];

                g_y_y_xy_xxzz[k] = -g_y_y_y_xxzz[k] * ab_x + g_y_y_y_xxxzz[k];

                g_y_y_xy_xyyy[k] = -g_y_y_y_xyyy[k] * ab_x + g_y_y_y_xxyyy[k];

                g_y_y_xy_xyyz[k] = -g_y_y_y_xyyz[k] * ab_x + g_y_y_y_xxyyz[k];

                g_y_y_xy_xyzz[k] = -g_y_y_y_xyzz[k] * ab_x + g_y_y_y_xxyzz[k];

                g_y_y_xy_xzzz[k] = -g_y_y_y_xzzz[k] * ab_x + g_y_y_y_xxzzz[k];

                g_y_y_xy_yyyy[k] = -g_y_y_y_yyyy[k] * ab_x + g_y_y_y_xyyyy[k];

                g_y_y_xy_yyyz[k] = -g_y_y_y_yyyz[k] * ab_x + g_y_y_y_xyyyz[k];

                g_y_y_xy_yyzz[k] = -g_y_y_y_yyzz[k] * ab_x + g_y_y_y_xyyzz[k];

                g_y_y_xy_yzzz[k] = -g_y_y_y_yzzz[k] * ab_x + g_y_y_y_xyzzz[k];

                g_y_y_xy_zzzz[k] = -g_y_y_y_zzzz[k] * ab_x + g_y_y_y_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_xxxx, g_y_y_xz_xxxy, g_y_y_xz_xxxz, g_y_y_xz_xxyy, g_y_y_xz_xxyz, g_y_y_xz_xxzz, g_y_y_xz_xyyy, g_y_y_xz_xyyz, g_y_y_xz_xyzz, g_y_y_xz_xzzz, g_y_y_xz_yyyy, g_y_y_xz_yyyz, g_y_y_xz_yyzz, g_y_y_xz_yzzz, g_y_y_xz_zzzz, g_y_y_z_xxxx, g_y_y_z_xxxxx, g_y_y_z_xxxxy, g_y_y_z_xxxxz, g_y_y_z_xxxy, g_y_y_z_xxxyy, g_y_y_z_xxxyz, g_y_y_z_xxxz, g_y_y_z_xxxzz, g_y_y_z_xxyy, g_y_y_z_xxyyy, g_y_y_z_xxyyz, g_y_y_z_xxyz, g_y_y_z_xxyzz, g_y_y_z_xxzz, g_y_y_z_xxzzz, g_y_y_z_xyyy, g_y_y_z_xyyyy, g_y_y_z_xyyyz, g_y_y_z_xyyz, g_y_y_z_xyyzz, g_y_y_z_xyzz, g_y_y_z_xyzzz, g_y_y_z_xzzz, g_y_y_z_xzzzz, g_y_y_z_yyyy, g_y_y_z_yyyz, g_y_y_z_yyzz, g_y_y_z_yzzz, g_y_y_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_xxxx[k] = -g_y_y_z_xxxx[k] * ab_x + g_y_y_z_xxxxx[k];

                g_y_y_xz_xxxy[k] = -g_y_y_z_xxxy[k] * ab_x + g_y_y_z_xxxxy[k];

                g_y_y_xz_xxxz[k] = -g_y_y_z_xxxz[k] * ab_x + g_y_y_z_xxxxz[k];

                g_y_y_xz_xxyy[k] = -g_y_y_z_xxyy[k] * ab_x + g_y_y_z_xxxyy[k];

                g_y_y_xz_xxyz[k] = -g_y_y_z_xxyz[k] * ab_x + g_y_y_z_xxxyz[k];

                g_y_y_xz_xxzz[k] = -g_y_y_z_xxzz[k] * ab_x + g_y_y_z_xxxzz[k];

                g_y_y_xz_xyyy[k] = -g_y_y_z_xyyy[k] * ab_x + g_y_y_z_xxyyy[k];

                g_y_y_xz_xyyz[k] = -g_y_y_z_xyyz[k] * ab_x + g_y_y_z_xxyyz[k];

                g_y_y_xz_xyzz[k] = -g_y_y_z_xyzz[k] * ab_x + g_y_y_z_xxyzz[k];

                g_y_y_xz_xzzz[k] = -g_y_y_z_xzzz[k] * ab_x + g_y_y_z_xxzzz[k];

                g_y_y_xz_yyyy[k] = -g_y_y_z_yyyy[k] * ab_x + g_y_y_z_xyyyy[k];

                g_y_y_xz_yyyz[k] = -g_y_y_z_yyyz[k] * ab_x + g_y_y_z_xyyyz[k];

                g_y_y_xz_yyzz[k] = -g_y_y_z_yyzz[k] * ab_x + g_y_y_z_xyyzz[k];

                g_y_y_xz_yzzz[k] = -g_y_y_z_yzzz[k] * ab_x + g_y_y_z_xyzzz[k];

                g_y_y_xz_zzzz[k] = -g_y_y_z_zzzz[k] * ab_x + g_y_y_z_xzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxx, g_0_y_y_xxxy, g_0_y_y_xxxz, g_0_y_y_xxyy, g_0_y_y_xxyz, g_0_y_y_xxzz, g_0_y_y_xyyy, g_0_y_y_xyyz, g_0_y_y_xyzz, g_0_y_y_xzzz, g_0_y_y_yyyy, g_0_y_y_yyyz, g_0_y_y_yyzz, g_0_y_y_yzzz, g_0_y_y_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz, g_y_y_y_xxxx, g_y_y_y_xxxxy, g_y_y_y_xxxy, g_y_y_y_xxxyy, g_y_y_y_xxxyz, g_y_y_y_xxxz, g_y_y_y_xxyy, g_y_y_y_xxyyy, g_y_y_y_xxyyz, g_y_y_y_xxyz, g_y_y_y_xxyzz, g_y_y_y_xxzz, g_y_y_y_xyyy, g_y_y_y_xyyyy, g_y_y_y_xyyyz, g_y_y_y_xyyz, g_y_y_y_xyyzz, g_y_y_y_xyzz, g_y_y_y_xyzzz, g_y_y_y_xzzz, g_y_y_y_yyyy, g_y_y_y_yyyyy, g_y_y_y_yyyyz, g_y_y_y_yyyz, g_y_y_y_yyyzz, g_y_y_y_yyzz, g_y_y_y_yyzzz, g_y_y_y_yzzz, g_y_y_y_yzzzz, g_y_y_y_zzzz, g_y_y_yy_xxxx, g_y_y_yy_xxxy, g_y_y_yy_xxxz, g_y_y_yy_xxyy, g_y_y_yy_xxyz, g_y_y_yy_xxzz, g_y_y_yy_xyyy, g_y_y_yy_xyyz, g_y_y_yy_xyzz, g_y_y_yy_xzzz, g_y_y_yy_yyyy, g_y_y_yy_yyyz, g_y_y_yy_yyzz, g_y_y_yy_yzzz, g_y_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_xxxx[k] = -g_0_y_y_xxxx[k] + g_y_0_y_xxxx[k] - g_y_y_y_xxxx[k] * ab_y + g_y_y_y_xxxxy[k];

                g_y_y_yy_xxxy[k] = -g_0_y_y_xxxy[k] + g_y_0_y_xxxy[k] - g_y_y_y_xxxy[k] * ab_y + g_y_y_y_xxxyy[k];

                g_y_y_yy_xxxz[k] = -g_0_y_y_xxxz[k] + g_y_0_y_xxxz[k] - g_y_y_y_xxxz[k] * ab_y + g_y_y_y_xxxyz[k];

                g_y_y_yy_xxyy[k] = -g_0_y_y_xxyy[k] + g_y_0_y_xxyy[k] - g_y_y_y_xxyy[k] * ab_y + g_y_y_y_xxyyy[k];

                g_y_y_yy_xxyz[k] = -g_0_y_y_xxyz[k] + g_y_0_y_xxyz[k] - g_y_y_y_xxyz[k] * ab_y + g_y_y_y_xxyyz[k];

                g_y_y_yy_xxzz[k] = -g_0_y_y_xxzz[k] + g_y_0_y_xxzz[k] - g_y_y_y_xxzz[k] * ab_y + g_y_y_y_xxyzz[k];

                g_y_y_yy_xyyy[k] = -g_0_y_y_xyyy[k] + g_y_0_y_xyyy[k] - g_y_y_y_xyyy[k] * ab_y + g_y_y_y_xyyyy[k];

                g_y_y_yy_xyyz[k] = -g_0_y_y_xyyz[k] + g_y_0_y_xyyz[k] - g_y_y_y_xyyz[k] * ab_y + g_y_y_y_xyyyz[k];

                g_y_y_yy_xyzz[k] = -g_0_y_y_xyzz[k] + g_y_0_y_xyzz[k] - g_y_y_y_xyzz[k] * ab_y + g_y_y_y_xyyzz[k];

                g_y_y_yy_xzzz[k] = -g_0_y_y_xzzz[k] + g_y_0_y_xzzz[k] - g_y_y_y_xzzz[k] * ab_y + g_y_y_y_xyzzz[k];

                g_y_y_yy_yyyy[k] = -g_0_y_y_yyyy[k] + g_y_0_y_yyyy[k] - g_y_y_y_yyyy[k] * ab_y + g_y_y_y_yyyyy[k];

                g_y_y_yy_yyyz[k] = -g_0_y_y_yyyz[k] + g_y_0_y_yyyz[k] - g_y_y_y_yyyz[k] * ab_y + g_y_y_y_yyyyz[k];

                g_y_y_yy_yyzz[k] = -g_0_y_y_yyzz[k] + g_y_0_y_yyzz[k] - g_y_y_y_yyzz[k] * ab_y + g_y_y_y_yyyzz[k];

                g_y_y_yy_yzzz[k] = -g_0_y_y_yzzz[k] + g_y_0_y_yzzz[k] - g_y_y_y_yzzz[k] * ab_y + g_y_y_y_yyzzz[k];

                g_y_y_yy_zzzz[k] = -g_0_y_y_zzzz[k] + g_y_0_y_zzzz[k] - g_y_y_y_zzzz[k] * ab_y + g_y_y_y_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_xxxx, g_y_y_y_xxxxz, g_y_y_y_xxxy, g_y_y_y_xxxyz, g_y_y_y_xxxz, g_y_y_y_xxxzz, g_y_y_y_xxyy, g_y_y_y_xxyyz, g_y_y_y_xxyz, g_y_y_y_xxyzz, g_y_y_y_xxzz, g_y_y_y_xxzzz, g_y_y_y_xyyy, g_y_y_y_xyyyz, g_y_y_y_xyyz, g_y_y_y_xyyzz, g_y_y_y_xyzz, g_y_y_y_xyzzz, g_y_y_y_xzzz, g_y_y_y_xzzzz, g_y_y_y_yyyy, g_y_y_y_yyyyz, g_y_y_y_yyyz, g_y_y_y_yyyzz, g_y_y_y_yyzz, g_y_y_y_yyzzz, g_y_y_y_yzzz, g_y_y_y_yzzzz, g_y_y_y_zzzz, g_y_y_y_zzzzz, g_y_y_yz_xxxx, g_y_y_yz_xxxy, g_y_y_yz_xxxz, g_y_y_yz_xxyy, g_y_y_yz_xxyz, g_y_y_yz_xxzz, g_y_y_yz_xyyy, g_y_y_yz_xyyz, g_y_y_yz_xyzz, g_y_y_yz_xzzz, g_y_y_yz_yyyy, g_y_y_yz_yyyz, g_y_y_yz_yyzz, g_y_y_yz_yzzz, g_y_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_xxxx[k] = -g_y_y_y_xxxx[k] * ab_z + g_y_y_y_xxxxz[k];

                g_y_y_yz_xxxy[k] = -g_y_y_y_xxxy[k] * ab_z + g_y_y_y_xxxyz[k];

                g_y_y_yz_xxxz[k] = -g_y_y_y_xxxz[k] * ab_z + g_y_y_y_xxxzz[k];

                g_y_y_yz_xxyy[k] = -g_y_y_y_xxyy[k] * ab_z + g_y_y_y_xxyyz[k];

                g_y_y_yz_xxyz[k] = -g_y_y_y_xxyz[k] * ab_z + g_y_y_y_xxyzz[k];

                g_y_y_yz_xxzz[k] = -g_y_y_y_xxzz[k] * ab_z + g_y_y_y_xxzzz[k];

                g_y_y_yz_xyyy[k] = -g_y_y_y_xyyy[k] * ab_z + g_y_y_y_xyyyz[k];

                g_y_y_yz_xyyz[k] = -g_y_y_y_xyyz[k] * ab_z + g_y_y_y_xyyzz[k];

                g_y_y_yz_xyzz[k] = -g_y_y_y_xyzz[k] * ab_z + g_y_y_y_xyzzz[k];

                g_y_y_yz_xzzz[k] = -g_y_y_y_xzzz[k] * ab_z + g_y_y_y_xzzzz[k];

                g_y_y_yz_yyyy[k] = -g_y_y_y_yyyy[k] * ab_z + g_y_y_y_yyyyz[k];

                g_y_y_yz_yyyz[k] = -g_y_y_y_yyyz[k] * ab_z + g_y_y_y_yyyzz[k];

                g_y_y_yz_yyzz[k] = -g_y_y_y_yyzz[k] * ab_z + g_y_y_y_yyzzz[k];

                g_y_y_yz_yzzz[k] = -g_y_y_y_yzzz[k] * ab_z + g_y_y_y_yzzzz[k];

                g_y_y_yz_zzzz[k] = -g_y_y_y_zzzz[k] * ab_z + g_y_y_y_zzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_xxxx, g_y_y_z_xxxxz, g_y_y_z_xxxy, g_y_y_z_xxxyz, g_y_y_z_xxxz, g_y_y_z_xxxzz, g_y_y_z_xxyy, g_y_y_z_xxyyz, g_y_y_z_xxyz, g_y_y_z_xxyzz, g_y_y_z_xxzz, g_y_y_z_xxzzz, g_y_y_z_xyyy, g_y_y_z_xyyyz, g_y_y_z_xyyz, g_y_y_z_xyyzz, g_y_y_z_xyzz, g_y_y_z_xyzzz, g_y_y_z_xzzz, g_y_y_z_xzzzz, g_y_y_z_yyyy, g_y_y_z_yyyyz, g_y_y_z_yyyz, g_y_y_z_yyyzz, g_y_y_z_yyzz, g_y_y_z_yyzzz, g_y_y_z_yzzz, g_y_y_z_yzzzz, g_y_y_z_zzzz, g_y_y_z_zzzzz, g_y_y_zz_xxxx, g_y_y_zz_xxxy, g_y_y_zz_xxxz, g_y_y_zz_xxyy, g_y_y_zz_xxyz, g_y_y_zz_xxzz, g_y_y_zz_xyyy, g_y_y_zz_xyyz, g_y_y_zz_xyzz, g_y_y_zz_xzzz, g_y_y_zz_yyyy, g_y_y_zz_yyyz, g_y_y_zz_yyzz, g_y_y_zz_yzzz, g_y_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_xxxx[k] = -g_y_y_z_xxxx[k] * ab_z + g_y_y_z_xxxxz[k];

                g_y_y_zz_xxxy[k] = -g_y_y_z_xxxy[k] * ab_z + g_y_y_z_xxxyz[k];

                g_y_y_zz_xxxz[k] = -g_y_y_z_xxxz[k] * ab_z + g_y_y_z_xxxzz[k];

                g_y_y_zz_xxyy[k] = -g_y_y_z_xxyy[k] * ab_z + g_y_y_z_xxyyz[k];

                g_y_y_zz_xxyz[k] = -g_y_y_z_xxyz[k] * ab_z + g_y_y_z_xxyzz[k];

                g_y_y_zz_xxzz[k] = -g_y_y_z_xxzz[k] * ab_z + g_y_y_z_xxzzz[k];

                g_y_y_zz_xyyy[k] = -g_y_y_z_xyyy[k] * ab_z + g_y_y_z_xyyyz[k];

                g_y_y_zz_xyyz[k] = -g_y_y_z_xyyz[k] * ab_z + g_y_y_z_xyyzz[k];

                g_y_y_zz_xyzz[k] = -g_y_y_z_xyzz[k] * ab_z + g_y_y_z_xyzzz[k];

                g_y_y_zz_xzzz[k] = -g_y_y_z_xzzz[k] * ab_z + g_y_y_z_xzzzz[k];

                g_y_y_zz_yyyy[k] = -g_y_y_z_yyyy[k] * ab_z + g_y_y_z_yyyyz[k];

                g_y_y_zz_yyyz[k] = -g_y_y_z_yyyz[k] * ab_z + g_y_y_z_yyyzz[k];

                g_y_y_zz_yyzz[k] = -g_y_y_z_yyzz[k] * ab_z + g_y_y_z_yyzzz[k];

                g_y_y_zz_yzzz[k] = -g_y_y_z_yzzz[k] * ab_z + g_y_y_z_yzzzz[k];

                g_y_y_zz_zzzz[k] = -g_y_y_z_zzzz[k] * ab_z + g_y_y_z_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_xxxx, g_y_z_x_xxxxx, g_y_z_x_xxxxy, g_y_z_x_xxxxz, g_y_z_x_xxxy, g_y_z_x_xxxyy, g_y_z_x_xxxyz, g_y_z_x_xxxz, g_y_z_x_xxxzz, g_y_z_x_xxyy, g_y_z_x_xxyyy, g_y_z_x_xxyyz, g_y_z_x_xxyz, g_y_z_x_xxyzz, g_y_z_x_xxzz, g_y_z_x_xxzzz, g_y_z_x_xyyy, g_y_z_x_xyyyy, g_y_z_x_xyyyz, g_y_z_x_xyyz, g_y_z_x_xyyzz, g_y_z_x_xyzz, g_y_z_x_xyzzz, g_y_z_x_xzzz, g_y_z_x_xzzzz, g_y_z_x_yyyy, g_y_z_x_yyyz, g_y_z_x_yyzz, g_y_z_x_yzzz, g_y_z_x_zzzz, g_y_z_xx_xxxx, g_y_z_xx_xxxy, g_y_z_xx_xxxz, g_y_z_xx_xxyy, g_y_z_xx_xxyz, g_y_z_xx_xxzz, g_y_z_xx_xyyy, g_y_z_xx_xyyz, g_y_z_xx_xyzz, g_y_z_xx_xzzz, g_y_z_xx_yyyy, g_y_z_xx_yyyz, g_y_z_xx_yyzz, g_y_z_xx_yzzz, g_y_z_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_xxxx[k] = -g_y_z_x_xxxx[k] * ab_x + g_y_z_x_xxxxx[k];

                g_y_z_xx_xxxy[k] = -g_y_z_x_xxxy[k] * ab_x + g_y_z_x_xxxxy[k];

                g_y_z_xx_xxxz[k] = -g_y_z_x_xxxz[k] * ab_x + g_y_z_x_xxxxz[k];

                g_y_z_xx_xxyy[k] = -g_y_z_x_xxyy[k] * ab_x + g_y_z_x_xxxyy[k];

                g_y_z_xx_xxyz[k] = -g_y_z_x_xxyz[k] * ab_x + g_y_z_x_xxxyz[k];

                g_y_z_xx_xxzz[k] = -g_y_z_x_xxzz[k] * ab_x + g_y_z_x_xxxzz[k];

                g_y_z_xx_xyyy[k] = -g_y_z_x_xyyy[k] * ab_x + g_y_z_x_xxyyy[k];

                g_y_z_xx_xyyz[k] = -g_y_z_x_xyyz[k] * ab_x + g_y_z_x_xxyyz[k];

                g_y_z_xx_xyzz[k] = -g_y_z_x_xyzz[k] * ab_x + g_y_z_x_xxyzz[k];

                g_y_z_xx_xzzz[k] = -g_y_z_x_xzzz[k] * ab_x + g_y_z_x_xxzzz[k];

                g_y_z_xx_yyyy[k] = -g_y_z_x_yyyy[k] * ab_x + g_y_z_x_xyyyy[k];

                g_y_z_xx_yyyz[k] = -g_y_z_x_yyyz[k] * ab_x + g_y_z_x_xyyyz[k];

                g_y_z_xx_yyzz[k] = -g_y_z_x_yyzz[k] * ab_x + g_y_z_x_xyyzz[k];

                g_y_z_xx_yzzz[k] = -g_y_z_x_yzzz[k] * ab_x + g_y_z_x_xyzzz[k];

                g_y_z_xx_zzzz[k] = -g_y_z_x_zzzz[k] * ab_x + g_y_z_x_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_xxxx, g_y_z_xy_xxxy, g_y_z_xy_xxxz, g_y_z_xy_xxyy, g_y_z_xy_xxyz, g_y_z_xy_xxzz, g_y_z_xy_xyyy, g_y_z_xy_xyyz, g_y_z_xy_xyzz, g_y_z_xy_xzzz, g_y_z_xy_yyyy, g_y_z_xy_yyyz, g_y_z_xy_yyzz, g_y_z_xy_yzzz, g_y_z_xy_zzzz, g_y_z_y_xxxx, g_y_z_y_xxxxx, g_y_z_y_xxxxy, g_y_z_y_xxxxz, g_y_z_y_xxxy, g_y_z_y_xxxyy, g_y_z_y_xxxyz, g_y_z_y_xxxz, g_y_z_y_xxxzz, g_y_z_y_xxyy, g_y_z_y_xxyyy, g_y_z_y_xxyyz, g_y_z_y_xxyz, g_y_z_y_xxyzz, g_y_z_y_xxzz, g_y_z_y_xxzzz, g_y_z_y_xyyy, g_y_z_y_xyyyy, g_y_z_y_xyyyz, g_y_z_y_xyyz, g_y_z_y_xyyzz, g_y_z_y_xyzz, g_y_z_y_xyzzz, g_y_z_y_xzzz, g_y_z_y_xzzzz, g_y_z_y_yyyy, g_y_z_y_yyyz, g_y_z_y_yyzz, g_y_z_y_yzzz, g_y_z_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_xxxx[k] = -g_y_z_y_xxxx[k] * ab_x + g_y_z_y_xxxxx[k];

                g_y_z_xy_xxxy[k] = -g_y_z_y_xxxy[k] * ab_x + g_y_z_y_xxxxy[k];

                g_y_z_xy_xxxz[k] = -g_y_z_y_xxxz[k] * ab_x + g_y_z_y_xxxxz[k];

                g_y_z_xy_xxyy[k] = -g_y_z_y_xxyy[k] * ab_x + g_y_z_y_xxxyy[k];

                g_y_z_xy_xxyz[k] = -g_y_z_y_xxyz[k] * ab_x + g_y_z_y_xxxyz[k];

                g_y_z_xy_xxzz[k] = -g_y_z_y_xxzz[k] * ab_x + g_y_z_y_xxxzz[k];

                g_y_z_xy_xyyy[k] = -g_y_z_y_xyyy[k] * ab_x + g_y_z_y_xxyyy[k];

                g_y_z_xy_xyyz[k] = -g_y_z_y_xyyz[k] * ab_x + g_y_z_y_xxyyz[k];

                g_y_z_xy_xyzz[k] = -g_y_z_y_xyzz[k] * ab_x + g_y_z_y_xxyzz[k];

                g_y_z_xy_xzzz[k] = -g_y_z_y_xzzz[k] * ab_x + g_y_z_y_xxzzz[k];

                g_y_z_xy_yyyy[k] = -g_y_z_y_yyyy[k] * ab_x + g_y_z_y_xyyyy[k];

                g_y_z_xy_yyyz[k] = -g_y_z_y_yyyz[k] * ab_x + g_y_z_y_xyyyz[k];

                g_y_z_xy_yyzz[k] = -g_y_z_y_yyzz[k] * ab_x + g_y_z_y_xyyzz[k];

                g_y_z_xy_yzzz[k] = -g_y_z_y_yzzz[k] * ab_x + g_y_z_y_xyzzz[k];

                g_y_z_xy_zzzz[k] = -g_y_z_y_zzzz[k] * ab_x + g_y_z_y_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_xxxx, g_y_z_xz_xxxy, g_y_z_xz_xxxz, g_y_z_xz_xxyy, g_y_z_xz_xxyz, g_y_z_xz_xxzz, g_y_z_xz_xyyy, g_y_z_xz_xyyz, g_y_z_xz_xyzz, g_y_z_xz_xzzz, g_y_z_xz_yyyy, g_y_z_xz_yyyz, g_y_z_xz_yyzz, g_y_z_xz_yzzz, g_y_z_xz_zzzz, g_y_z_z_xxxx, g_y_z_z_xxxxx, g_y_z_z_xxxxy, g_y_z_z_xxxxz, g_y_z_z_xxxy, g_y_z_z_xxxyy, g_y_z_z_xxxyz, g_y_z_z_xxxz, g_y_z_z_xxxzz, g_y_z_z_xxyy, g_y_z_z_xxyyy, g_y_z_z_xxyyz, g_y_z_z_xxyz, g_y_z_z_xxyzz, g_y_z_z_xxzz, g_y_z_z_xxzzz, g_y_z_z_xyyy, g_y_z_z_xyyyy, g_y_z_z_xyyyz, g_y_z_z_xyyz, g_y_z_z_xyyzz, g_y_z_z_xyzz, g_y_z_z_xyzzz, g_y_z_z_xzzz, g_y_z_z_xzzzz, g_y_z_z_yyyy, g_y_z_z_yyyz, g_y_z_z_yyzz, g_y_z_z_yzzz, g_y_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_xxxx[k] = -g_y_z_z_xxxx[k] * ab_x + g_y_z_z_xxxxx[k];

                g_y_z_xz_xxxy[k] = -g_y_z_z_xxxy[k] * ab_x + g_y_z_z_xxxxy[k];

                g_y_z_xz_xxxz[k] = -g_y_z_z_xxxz[k] * ab_x + g_y_z_z_xxxxz[k];

                g_y_z_xz_xxyy[k] = -g_y_z_z_xxyy[k] * ab_x + g_y_z_z_xxxyy[k];

                g_y_z_xz_xxyz[k] = -g_y_z_z_xxyz[k] * ab_x + g_y_z_z_xxxyz[k];

                g_y_z_xz_xxzz[k] = -g_y_z_z_xxzz[k] * ab_x + g_y_z_z_xxxzz[k];

                g_y_z_xz_xyyy[k] = -g_y_z_z_xyyy[k] * ab_x + g_y_z_z_xxyyy[k];

                g_y_z_xz_xyyz[k] = -g_y_z_z_xyyz[k] * ab_x + g_y_z_z_xxyyz[k];

                g_y_z_xz_xyzz[k] = -g_y_z_z_xyzz[k] * ab_x + g_y_z_z_xxyzz[k];

                g_y_z_xz_xzzz[k] = -g_y_z_z_xzzz[k] * ab_x + g_y_z_z_xxzzz[k];

                g_y_z_xz_yyyy[k] = -g_y_z_z_yyyy[k] * ab_x + g_y_z_z_xyyyy[k];

                g_y_z_xz_yyyz[k] = -g_y_z_z_yyyz[k] * ab_x + g_y_z_z_xyyyz[k];

                g_y_z_xz_yyzz[k] = -g_y_z_z_yyzz[k] * ab_x + g_y_z_z_xyyzz[k];

                g_y_z_xz_yzzz[k] = -g_y_z_z_yzzz[k] * ab_x + g_y_z_z_xyzzz[k];

                g_y_z_xz_zzzz[k] = -g_y_z_z_zzzz[k] * ab_x + g_y_z_z_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxxx, g_0_z_y_xxxy, g_0_z_y_xxxz, g_0_z_y_xxyy, g_0_z_y_xxyz, g_0_z_y_xxzz, g_0_z_y_xyyy, g_0_z_y_xyyz, g_0_z_y_xyzz, g_0_z_y_xzzz, g_0_z_y_yyyy, g_0_z_y_yyyz, g_0_z_y_yyzz, g_0_z_y_yzzz, g_0_z_y_zzzz, g_y_z_y_xxxx, g_y_z_y_xxxxy, g_y_z_y_xxxy, g_y_z_y_xxxyy, g_y_z_y_xxxyz, g_y_z_y_xxxz, g_y_z_y_xxyy, g_y_z_y_xxyyy, g_y_z_y_xxyyz, g_y_z_y_xxyz, g_y_z_y_xxyzz, g_y_z_y_xxzz, g_y_z_y_xyyy, g_y_z_y_xyyyy, g_y_z_y_xyyyz, g_y_z_y_xyyz, g_y_z_y_xyyzz, g_y_z_y_xyzz, g_y_z_y_xyzzz, g_y_z_y_xzzz, g_y_z_y_yyyy, g_y_z_y_yyyyy, g_y_z_y_yyyyz, g_y_z_y_yyyz, g_y_z_y_yyyzz, g_y_z_y_yyzz, g_y_z_y_yyzzz, g_y_z_y_yzzz, g_y_z_y_yzzzz, g_y_z_y_zzzz, g_y_z_yy_xxxx, g_y_z_yy_xxxy, g_y_z_yy_xxxz, g_y_z_yy_xxyy, g_y_z_yy_xxyz, g_y_z_yy_xxzz, g_y_z_yy_xyyy, g_y_z_yy_xyyz, g_y_z_yy_xyzz, g_y_z_yy_xzzz, g_y_z_yy_yyyy, g_y_z_yy_yyyz, g_y_z_yy_yyzz, g_y_z_yy_yzzz, g_y_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_xxxx[k] = -g_0_z_y_xxxx[k] - g_y_z_y_xxxx[k] * ab_y + g_y_z_y_xxxxy[k];

                g_y_z_yy_xxxy[k] = -g_0_z_y_xxxy[k] - g_y_z_y_xxxy[k] * ab_y + g_y_z_y_xxxyy[k];

                g_y_z_yy_xxxz[k] = -g_0_z_y_xxxz[k] - g_y_z_y_xxxz[k] * ab_y + g_y_z_y_xxxyz[k];

                g_y_z_yy_xxyy[k] = -g_0_z_y_xxyy[k] - g_y_z_y_xxyy[k] * ab_y + g_y_z_y_xxyyy[k];

                g_y_z_yy_xxyz[k] = -g_0_z_y_xxyz[k] - g_y_z_y_xxyz[k] * ab_y + g_y_z_y_xxyyz[k];

                g_y_z_yy_xxzz[k] = -g_0_z_y_xxzz[k] - g_y_z_y_xxzz[k] * ab_y + g_y_z_y_xxyzz[k];

                g_y_z_yy_xyyy[k] = -g_0_z_y_xyyy[k] - g_y_z_y_xyyy[k] * ab_y + g_y_z_y_xyyyy[k];

                g_y_z_yy_xyyz[k] = -g_0_z_y_xyyz[k] - g_y_z_y_xyyz[k] * ab_y + g_y_z_y_xyyyz[k];

                g_y_z_yy_xyzz[k] = -g_0_z_y_xyzz[k] - g_y_z_y_xyzz[k] * ab_y + g_y_z_y_xyyzz[k];

                g_y_z_yy_xzzz[k] = -g_0_z_y_xzzz[k] - g_y_z_y_xzzz[k] * ab_y + g_y_z_y_xyzzz[k];

                g_y_z_yy_yyyy[k] = -g_0_z_y_yyyy[k] - g_y_z_y_yyyy[k] * ab_y + g_y_z_y_yyyyy[k];

                g_y_z_yy_yyyz[k] = -g_0_z_y_yyyz[k] - g_y_z_y_yyyz[k] * ab_y + g_y_z_y_yyyyz[k];

                g_y_z_yy_yyzz[k] = -g_0_z_y_yyzz[k] - g_y_z_y_yyzz[k] * ab_y + g_y_z_y_yyyzz[k];

                g_y_z_yy_yzzz[k] = -g_0_z_y_yzzz[k] - g_y_z_y_yzzz[k] * ab_y + g_y_z_y_yyzzz[k];

                g_y_z_yy_zzzz[k] = -g_0_z_y_zzzz[k] - g_y_z_y_zzzz[k] * ab_y + g_y_z_y_yzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxx, g_0_z_z_xxxy, g_0_z_z_xxxz, g_0_z_z_xxyy, g_0_z_z_xxyz, g_0_z_z_xxzz, g_0_z_z_xyyy, g_0_z_z_xyyz, g_0_z_z_xyzz, g_0_z_z_xzzz, g_0_z_z_yyyy, g_0_z_z_yyyz, g_0_z_z_yyzz, g_0_z_z_yzzz, g_0_z_z_zzzz, g_y_z_yz_xxxx, g_y_z_yz_xxxy, g_y_z_yz_xxxz, g_y_z_yz_xxyy, g_y_z_yz_xxyz, g_y_z_yz_xxzz, g_y_z_yz_xyyy, g_y_z_yz_xyyz, g_y_z_yz_xyzz, g_y_z_yz_xzzz, g_y_z_yz_yyyy, g_y_z_yz_yyyz, g_y_z_yz_yyzz, g_y_z_yz_yzzz, g_y_z_yz_zzzz, g_y_z_z_xxxx, g_y_z_z_xxxxy, g_y_z_z_xxxy, g_y_z_z_xxxyy, g_y_z_z_xxxyz, g_y_z_z_xxxz, g_y_z_z_xxyy, g_y_z_z_xxyyy, g_y_z_z_xxyyz, g_y_z_z_xxyz, g_y_z_z_xxyzz, g_y_z_z_xxzz, g_y_z_z_xyyy, g_y_z_z_xyyyy, g_y_z_z_xyyyz, g_y_z_z_xyyz, g_y_z_z_xyyzz, g_y_z_z_xyzz, g_y_z_z_xyzzz, g_y_z_z_xzzz, g_y_z_z_yyyy, g_y_z_z_yyyyy, g_y_z_z_yyyyz, g_y_z_z_yyyz, g_y_z_z_yyyzz, g_y_z_z_yyzz, g_y_z_z_yyzzz, g_y_z_z_yzzz, g_y_z_z_yzzzz, g_y_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_xxxx[k] = -g_0_z_z_xxxx[k] - g_y_z_z_xxxx[k] * ab_y + g_y_z_z_xxxxy[k];

                g_y_z_yz_xxxy[k] = -g_0_z_z_xxxy[k] - g_y_z_z_xxxy[k] * ab_y + g_y_z_z_xxxyy[k];

                g_y_z_yz_xxxz[k] = -g_0_z_z_xxxz[k] - g_y_z_z_xxxz[k] * ab_y + g_y_z_z_xxxyz[k];

                g_y_z_yz_xxyy[k] = -g_0_z_z_xxyy[k] - g_y_z_z_xxyy[k] * ab_y + g_y_z_z_xxyyy[k];

                g_y_z_yz_xxyz[k] = -g_0_z_z_xxyz[k] - g_y_z_z_xxyz[k] * ab_y + g_y_z_z_xxyyz[k];

                g_y_z_yz_xxzz[k] = -g_0_z_z_xxzz[k] - g_y_z_z_xxzz[k] * ab_y + g_y_z_z_xxyzz[k];

                g_y_z_yz_xyyy[k] = -g_0_z_z_xyyy[k] - g_y_z_z_xyyy[k] * ab_y + g_y_z_z_xyyyy[k];

                g_y_z_yz_xyyz[k] = -g_0_z_z_xyyz[k] - g_y_z_z_xyyz[k] * ab_y + g_y_z_z_xyyyz[k];

                g_y_z_yz_xyzz[k] = -g_0_z_z_xyzz[k] - g_y_z_z_xyzz[k] * ab_y + g_y_z_z_xyyzz[k];

                g_y_z_yz_xzzz[k] = -g_0_z_z_xzzz[k] - g_y_z_z_xzzz[k] * ab_y + g_y_z_z_xyzzz[k];

                g_y_z_yz_yyyy[k] = -g_0_z_z_yyyy[k] - g_y_z_z_yyyy[k] * ab_y + g_y_z_z_yyyyy[k];

                g_y_z_yz_yyyz[k] = -g_0_z_z_yyyz[k] - g_y_z_z_yyyz[k] * ab_y + g_y_z_z_yyyyz[k];

                g_y_z_yz_yyzz[k] = -g_0_z_z_yyzz[k] - g_y_z_z_yyzz[k] * ab_y + g_y_z_z_yyyzz[k];

                g_y_z_yz_yzzz[k] = -g_0_z_z_yzzz[k] - g_y_z_z_yzzz[k] * ab_y + g_y_z_z_yyzzz[k];

                g_y_z_yz_zzzz[k] = -g_0_z_z_zzzz[k] - g_y_z_z_zzzz[k] * ab_y + g_y_z_z_yzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxzz, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyzz, g_y_0_z_xzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz, g_y_z_z_xxxx, g_y_z_z_xxxxz, g_y_z_z_xxxy, g_y_z_z_xxxyz, g_y_z_z_xxxz, g_y_z_z_xxxzz, g_y_z_z_xxyy, g_y_z_z_xxyyz, g_y_z_z_xxyz, g_y_z_z_xxyzz, g_y_z_z_xxzz, g_y_z_z_xxzzz, g_y_z_z_xyyy, g_y_z_z_xyyyz, g_y_z_z_xyyz, g_y_z_z_xyyzz, g_y_z_z_xyzz, g_y_z_z_xyzzz, g_y_z_z_xzzz, g_y_z_z_xzzzz, g_y_z_z_yyyy, g_y_z_z_yyyyz, g_y_z_z_yyyz, g_y_z_z_yyyzz, g_y_z_z_yyzz, g_y_z_z_yyzzz, g_y_z_z_yzzz, g_y_z_z_yzzzz, g_y_z_z_zzzz, g_y_z_z_zzzzz, g_y_z_zz_xxxx, g_y_z_zz_xxxy, g_y_z_zz_xxxz, g_y_z_zz_xxyy, g_y_z_zz_xxyz, g_y_z_zz_xxzz, g_y_z_zz_xyyy, g_y_z_zz_xyyz, g_y_z_zz_xyzz, g_y_z_zz_xzzz, g_y_z_zz_yyyy, g_y_z_zz_yyyz, g_y_z_zz_yyzz, g_y_z_zz_yzzz, g_y_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_xxxx[k] = g_y_0_z_xxxx[k] - g_y_z_z_xxxx[k] * ab_z + g_y_z_z_xxxxz[k];

                g_y_z_zz_xxxy[k] = g_y_0_z_xxxy[k] - g_y_z_z_xxxy[k] * ab_z + g_y_z_z_xxxyz[k];

                g_y_z_zz_xxxz[k] = g_y_0_z_xxxz[k] - g_y_z_z_xxxz[k] * ab_z + g_y_z_z_xxxzz[k];

                g_y_z_zz_xxyy[k] = g_y_0_z_xxyy[k] - g_y_z_z_xxyy[k] * ab_z + g_y_z_z_xxyyz[k];

                g_y_z_zz_xxyz[k] = g_y_0_z_xxyz[k] - g_y_z_z_xxyz[k] * ab_z + g_y_z_z_xxyzz[k];

                g_y_z_zz_xxzz[k] = g_y_0_z_xxzz[k] - g_y_z_z_xxzz[k] * ab_z + g_y_z_z_xxzzz[k];

                g_y_z_zz_xyyy[k] = g_y_0_z_xyyy[k] - g_y_z_z_xyyy[k] * ab_z + g_y_z_z_xyyyz[k];

                g_y_z_zz_xyyz[k] = g_y_0_z_xyyz[k] - g_y_z_z_xyyz[k] * ab_z + g_y_z_z_xyyzz[k];

                g_y_z_zz_xyzz[k] = g_y_0_z_xyzz[k] - g_y_z_z_xyzz[k] * ab_z + g_y_z_z_xyzzz[k];

                g_y_z_zz_xzzz[k] = g_y_0_z_xzzz[k] - g_y_z_z_xzzz[k] * ab_z + g_y_z_z_xzzzz[k];

                g_y_z_zz_yyyy[k] = g_y_0_z_yyyy[k] - g_y_z_z_yyyy[k] * ab_z + g_y_z_z_yyyyz[k];

                g_y_z_zz_yyyz[k] = g_y_0_z_yyyz[k] - g_y_z_z_yyyz[k] * ab_z + g_y_z_z_yyyzz[k];

                g_y_z_zz_yyzz[k] = g_y_0_z_yyzz[k] - g_y_z_z_yyzz[k] * ab_z + g_y_z_z_yyzzz[k];

                g_y_z_zz_yzzz[k] = g_y_0_z_yzzz[k] - g_y_z_z_yzzz[k] * ab_z + g_y_z_z_yzzzz[k];

                g_y_z_zz_zzzz[k] = g_y_0_z_zzzz[k] - g_y_z_z_zzzz[k] * ab_z + g_y_z_z_zzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxzz, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyzz, g_z_0_x_xzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz, g_z_x_x_xxxx, g_z_x_x_xxxxx, g_z_x_x_xxxxy, g_z_x_x_xxxxz, g_z_x_x_xxxy, g_z_x_x_xxxyy, g_z_x_x_xxxyz, g_z_x_x_xxxz, g_z_x_x_xxxzz, g_z_x_x_xxyy, g_z_x_x_xxyyy, g_z_x_x_xxyyz, g_z_x_x_xxyz, g_z_x_x_xxyzz, g_z_x_x_xxzz, g_z_x_x_xxzzz, g_z_x_x_xyyy, g_z_x_x_xyyyy, g_z_x_x_xyyyz, g_z_x_x_xyyz, g_z_x_x_xyyzz, g_z_x_x_xyzz, g_z_x_x_xyzzz, g_z_x_x_xzzz, g_z_x_x_xzzzz, g_z_x_x_yyyy, g_z_x_x_yyyz, g_z_x_x_yyzz, g_z_x_x_yzzz, g_z_x_x_zzzz, g_z_x_xx_xxxx, g_z_x_xx_xxxy, g_z_x_xx_xxxz, g_z_x_xx_xxyy, g_z_x_xx_xxyz, g_z_x_xx_xxzz, g_z_x_xx_xyyy, g_z_x_xx_xyyz, g_z_x_xx_xyzz, g_z_x_xx_xzzz, g_z_x_xx_yyyy, g_z_x_xx_yyyz, g_z_x_xx_yyzz, g_z_x_xx_yzzz, g_z_x_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_xxxx[k] = g_z_0_x_xxxx[k] - g_z_x_x_xxxx[k] * ab_x + g_z_x_x_xxxxx[k];

                g_z_x_xx_xxxy[k] = g_z_0_x_xxxy[k] - g_z_x_x_xxxy[k] * ab_x + g_z_x_x_xxxxy[k];

                g_z_x_xx_xxxz[k] = g_z_0_x_xxxz[k] - g_z_x_x_xxxz[k] * ab_x + g_z_x_x_xxxxz[k];

                g_z_x_xx_xxyy[k] = g_z_0_x_xxyy[k] - g_z_x_x_xxyy[k] * ab_x + g_z_x_x_xxxyy[k];

                g_z_x_xx_xxyz[k] = g_z_0_x_xxyz[k] - g_z_x_x_xxyz[k] * ab_x + g_z_x_x_xxxyz[k];

                g_z_x_xx_xxzz[k] = g_z_0_x_xxzz[k] - g_z_x_x_xxzz[k] * ab_x + g_z_x_x_xxxzz[k];

                g_z_x_xx_xyyy[k] = g_z_0_x_xyyy[k] - g_z_x_x_xyyy[k] * ab_x + g_z_x_x_xxyyy[k];

                g_z_x_xx_xyyz[k] = g_z_0_x_xyyz[k] - g_z_x_x_xyyz[k] * ab_x + g_z_x_x_xxyyz[k];

                g_z_x_xx_xyzz[k] = g_z_0_x_xyzz[k] - g_z_x_x_xyzz[k] * ab_x + g_z_x_x_xxyzz[k];

                g_z_x_xx_xzzz[k] = g_z_0_x_xzzz[k] - g_z_x_x_xzzz[k] * ab_x + g_z_x_x_xxzzz[k];

                g_z_x_xx_yyyy[k] = g_z_0_x_yyyy[k] - g_z_x_x_yyyy[k] * ab_x + g_z_x_x_xyyyy[k];

                g_z_x_xx_yyyz[k] = g_z_0_x_yyyz[k] - g_z_x_x_yyyz[k] * ab_x + g_z_x_x_xyyyz[k];

                g_z_x_xx_yyzz[k] = g_z_0_x_yyzz[k] - g_z_x_x_yyzz[k] * ab_x + g_z_x_x_xyyzz[k];

                g_z_x_xx_yzzz[k] = g_z_0_x_yzzz[k] - g_z_x_x_yzzz[k] * ab_x + g_z_x_x_xyzzz[k];

                g_z_x_xx_zzzz[k] = g_z_0_x_zzzz[k] - g_z_x_x_zzzz[k] * ab_x + g_z_x_x_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 566 * ccomps * dcomps);

            auto g_z_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 567 * ccomps * dcomps);

            auto g_z_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 568 * ccomps * dcomps);

            auto g_z_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_xxxx, g_z_x_x_xxxxy, g_z_x_x_xxxy, g_z_x_x_xxxyy, g_z_x_x_xxxyz, g_z_x_x_xxxz, g_z_x_x_xxyy, g_z_x_x_xxyyy, g_z_x_x_xxyyz, g_z_x_x_xxyz, g_z_x_x_xxyzz, g_z_x_x_xxzz, g_z_x_x_xyyy, g_z_x_x_xyyyy, g_z_x_x_xyyyz, g_z_x_x_xyyz, g_z_x_x_xyyzz, g_z_x_x_xyzz, g_z_x_x_xyzzz, g_z_x_x_xzzz, g_z_x_x_yyyy, g_z_x_x_yyyyy, g_z_x_x_yyyyz, g_z_x_x_yyyz, g_z_x_x_yyyzz, g_z_x_x_yyzz, g_z_x_x_yyzzz, g_z_x_x_yzzz, g_z_x_x_yzzzz, g_z_x_x_zzzz, g_z_x_xy_xxxx, g_z_x_xy_xxxy, g_z_x_xy_xxxz, g_z_x_xy_xxyy, g_z_x_xy_xxyz, g_z_x_xy_xxzz, g_z_x_xy_xyyy, g_z_x_xy_xyyz, g_z_x_xy_xyzz, g_z_x_xy_xzzz, g_z_x_xy_yyyy, g_z_x_xy_yyyz, g_z_x_xy_yyzz, g_z_x_xy_yzzz, g_z_x_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_xxxx[k] = -g_z_x_x_xxxx[k] * ab_y + g_z_x_x_xxxxy[k];

                g_z_x_xy_xxxy[k] = -g_z_x_x_xxxy[k] * ab_y + g_z_x_x_xxxyy[k];

                g_z_x_xy_xxxz[k] = -g_z_x_x_xxxz[k] * ab_y + g_z_x_x_xxxyz[k];

                g_z_x_xy_xxyy[k] = -g_z_x_x_xxyy[k] * ab_y + g_z_x_x_xxyyy[k];

                g_z_x_xy_xxyz[k] = -g_z_x_x_xxyz[k] * ab_y + g_z_x_x_xxyyz[k];

                g_z_x_xy_xxzz[k] = -g_z_x_x_xxzz[k] * ab_y + g_z_x_x_xxyzz[k];

                g_z_x_xy_xyyy[k] = -g_z_x_x_xyyy[k] * ab_y + g_z_x_x_xyyyy[k];

                g_z_x_xy_xyyz[k] = -g_z_x_x_xyyz[k] * ab_y + g_z_x_x_xyyyz[k];

                g_z_x_xy_xyzz[k] = -g_z_x_x_xyzz[k] * ab_y + g_z_x_x_xyyzz[k];

                g_z_x_xy_xzzz[k] = -g_z_x_x_xzzz[k] * ab_y + g_z_x_x_xyzzz[k];

                g_z_x_xy_yyyy[k] = -g_z_x_x_yyyy[k] * ab_y + g_z_x_x_yyyyy[k];

                g_z_x_xy_yyyz[k] = -g_z_x_x_yyyz[k] * ab_y + g_z_x_x_yyyyz[k];

                g_z_x_xy_yyzz[k] = -g_z_x_x_yyzz[k] * ab_y + g_z_x_x_yyyzz[k];

                g_z_x_xy_yzzz[k] = -g_z_x_x_yzzz[k] * ab_y + g_z_x_x_yyzzz[k];

                g_z_x_xy_zzzz[k] = -g_z_x_x_zzzz[k] * ab_y + g_z_x_x_yzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 570 * ccomps * dcomps);

            auto g_z_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 571 * ccomps * dcomps);

            auto g_z_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 572 * ccomps * dcomps);

            auto g_z_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 573 * ccomps * dcomps);

            auto g_z_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 574 * ccomps * dcomps);

            auto g_z_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 575 * ccomps * dcomps);

            auto g_z_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 576 * ccomps * dcomps);

            auto g_z_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 577 * ccomps * dcomps);

            auto g_z_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 578 * ccomps * dcomps);

            auto g_z_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 579 * ccomps * dcomps);

            auto g_z_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 580 * ccomps * dcomps);

            auto g_z_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 581 * ccomps * dcomps);

            auto g_z_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 582 * ccomps * dcomps);

            auto g_z_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 583 * ccomps * dcomps);

            auto g_z_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz, g_z_x_xz_xxxx, g_z_x_xz_xxxy, g_z_x_xz_xxxz, g_z_x_xz_xxyy, g_z_x_xz_xxyz, g_z_x_xz_xxzz, g_z_x_xz_xyyy, g_z_x_xz_xyyz, g_z_x_xz_xyzz, g_z_x_xz_xzzz, g_z_x_xz_yyyy, g_z_x_xz_yyyz, g_z_x_xz_yyzz, g_z_x_xz_yzzz, g_z_x_xz_zzzz, g_z_x_z_xxxx, g_z_x_z_xxxxx, g_z_x_z_xxxxy, g_z_x_z_xxxxz, g_z_x_z_xxxy, g_z_x_z_xxxyy, g_z_x_z_xxxyz, g_z_x_z_xxxz, g_z_x_z_xxxzz, g_z_x_z_xxyy, g_z_x_z_xxyyy, g_z_x_z_xxyyz, g_z_x_z_xxyz, g_z_x_z_xxyzz, g_z_x_z_xxzz, g_z_x_z_xxzzz, g_z_x_z_xyyy, g_z_x_z_xyyyy, g_z_x_z_xyyyz, g_z_x_z_xyyz, g_z_x_z_xyyzz, g_z_x_z_xyzz, g_z_x_z_xyzzz, g_z_x_z_xzzz, g_z_x_z_xzzzz, g_z_x_z_yyyy, g_z_x_z_yyyz, g_z_x_z_yyzz, g_z_x_z_yzzz, g_z_x_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_xxxx[k] = g_z_0_z_xxxx[k] - g_z_x_z_xxxx[k] * ab_x + g_z_x_z_xxxxx[k];

                g_z_x_xz_xxxy[k] = g_z_0_z_xxxy[k] - g_z_x_z_xxxy[k] * ab_x + g_z_x_z_xxxxy[k];

                g_z_x_xz_xxxz[k] = g_z_0_z_xxxz[k] - g_z_x_z_xxxz[k] * ab_x + g_z_x_z_xxxxz[k];

                g_z_x_xz_xxyy[k] = g_z_0_z_xxyy[k] - g_z_x_z_xxyy[k] * ab_x + g_z_x_z_xxxyy[k];

                g_z_x_xz_xxyz[k] = g_z_0_z_xxyz[k] - g_z_x_z_xxyz[k] * ab_x + g_z_x_z_xxxyz[k];

                g_z_x_xz_xxzz[k] = g_z_0_z_xxzz[k] - g_z_x_z_xxzz[k] * ab_x + g_z_x_z_xxxzz[k];

                g_z_x_xz_xyyy[k] = g_z_0_z_xyyy[k] - g_z_x_z_xyyy[k] * ab_x + g_z_x_z_xxyyy[k];

                g_z_x_xz_xyyz[k] = g_z_0_z_xyyz[k] - g_z_x_z_xyyz[k] * ab_x + g_z_x_z_xxyyz[k];

                g_z_x_xz_xyzz[k] = g_z_0_z_xyzz[k] - g_z_x_z_xyzz[k] * ab_x + g_z_x_z_xxyzz[k];

                g_z_x_xz_xzzz[k] = g_z_0_z_xzzz[k] - g_z_x_z_xzzz[k] * ab_x + g_z_x_z_xxzzz[k];

                g_z_x_xz_yyyy[k] = g_z_0_z_yyyy[k] - g_z_x_z_yyyy[k] * ab_x + g_z_x_z_xyyyy[k];

                g_z_x_xz_yyyz[k] = g_z_0_z_yyyz[k] - g_z_x_z_yyyz[k] * ab_x + g_z_x_z_xyyyz[k];

                g_z_x_xz_yyzz[k] = g_z_0_z_yyzz[k] - g_z_x_z_yyzz[k] * ab_x + g_z_x_z_xyyzz[k];

                g_z_x_xz_yzzz[k] = g_z_0_z_yzzz[k] - g_z_x_z_yzzz[k] * ab_x + g_z_x_z_xyzzz[k];

                g_z_x_xz_zzzz[k] = g_z_0_z_zzzz[k] - g_z_x_z_zzzz[k] * ab_x + g_z_x_z_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 585 * ccomps * dcomps);

            auto g_z_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 586 * ccomps * dcomps);

            auto g_z_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 587 * ccomps * dcomps);

            auto g_z_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 588 * ccomps * dcomps);

            auto g_z_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 589 * ccomps * dcomps);

            auto g_z_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 590 * ccomps * dcomps);

            auto g_z_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 591 * ccomps * dcomps);

            auto g_z_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 592 * ccomps * dcomps);

            auto g_z_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 593 * ccomps * dcomps);

            auto g_z_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 594 * ccomps * dcomps);

            auto g_z_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 595 * ccomps * dcomps);

            auto g_z_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 596 * ccomps * dcomps);

            auto g_z_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 597 * ccomps * dcomps);

            auto g_z_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 598 * ccomps * dcomps);

            auto g_z_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_xxxx, g_z_x_y_xxxxy, g_z_x_y_xxxy, g_z_x_y_xxxyy, g_z_x_y_xxxyz, g_z_x_y_xxxz, g_z_x_y_xxyy, g_z_x_y_xxyyy, g_z_x_y_xxyyz, g_z_x_y_xxyz, g_z_x_y_xxyzz, g_z_x_y_xxzz, g_z_x_y_xyyy, g_z_x_y_xyyyy, g_z_x_y_xyyyz, g_z_x_y_xyyz, g_z_x_y_xyyzz, g_z_x_y_xyzz, g_z_x_y_xyzzz, g_z_x_y_xzzz, g_z_x_y_yyyy, g_z_x_y_yyyyy, g_z_x_y_yyyyz, g_z_x_y_yyyz, g_z_x_y_yyyzz, g_z_x_y_yyzz, g_z_x_y_yyzzz, g_z_x_y_yzzz, g_z_x_y_yzzzz, g_z_x_y_zzzz, g_z_x_yy_xxxx, g_z_x_yy_xxxy, g_z_x_yy_xxxz, g_z_x_yy_xxyy, g_z_x_yy_xxyz, g_z_x_yy_xxzz, g_z_x_yy_xyyy, g_z_x_yy_xyyz, g_z_x_yy_xyzz, g_z_x_yy_xzzz, g_z_x_yy_yyyy, g_z_x_yy_yyyz, g_z_x_yy_yyzz, g_z_x_yy_yzzz, g_z_x_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_xxxx[k] = -g_z_x_y_xxxx[k] * ab_y + g_z_x_y_xxxxy[k];

                g_z_x_yy_xxxy[k] = -g_z_x_y_xxxy[k] * ab_y + g_z_x_y_xxxyy[k];

                g_z_x_yy_xxxz[k] = -g_z_x_y_xxxz[k] * ab_y + g_z_x_y_xxxyz[k];

                g_z_x_yy_xxyy[k] = -g_z_x_y_xxyy[k] * ab_y + g_z_x_y_xxyyy[k];

                g_z_x_yy_xxyz[k] = -g_z_x_y_xxyz[k] * ab_y + g_z_x_y_xxyyz[k];

                g_z_x_yy_xxzz[k] = -g_z_x_y_xxzz[k] * ab_y + g_z_x_y_xxyzz[k];

                g_z_x_yy_xyyy[k] = -g_z_x_y_xyyy[k] * ab_y + g_z_x_y_xyyyy[k];

                g_z_x_yy_xyyz[k] = -g_z_x_y_xyyz[k] * ab_y + g_z_x_y_xyyyz[k];

                g_z_x_yy_xyzz[k] = -g_z_x_y_xyzz[k] * ab_y + g_z_x_y_xyyzz[k];

                g_z_x_yy_xzzz[k] = -g_z_x_y_xzzz[k] * ab_y + g_z_x_y_xyzzz[k];

                g_z_x_yy_yyyy[k] = -g_z_x_y_yyyy[k] * ab_y + g_z_x_y_yyyyy[k];

                g_z_x_yy_yyyz[k] = -g_z_x_y_yyyz[k] * ab_y + g_z_x_y_yyyyz[k];

                g_z_x_yy_yyzz[k] = -g_z_x_y_yyzz[k] * ab_y + g_z_x_y_yyyzz[k];

                g_z_x_yy_yzzz[k] = -g_z_x_y_yzzz[k] * ab_y + g_z_x_y_yyzzz[k];

                g_z_x_yy_zzzz[k] = -g_z_x_y_zzzz[k] * ab_y + g_z_x_y_yzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 605 * ccomps * dcomps);

            auto g_z_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 609 * ccomps * dcomps);

            auto g_z_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 611 * ccomps * dcomps);

            auto g_z_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_xxxx, g_z_x_yz_xxxy, g_z_x_yz_xxxz, g_z_x_yz_xxyy, g_z_x_yz_xxyz, g_z_x_yz_xxzz, g_z_x_yz_xyyy, g_z_x_yz_xyyz, g_z_x_yz_xyzz, g_z_x_yz_xzzz, g_z_x_yz_yyyy, g_z_x_yz_yyyz, g_z_x_yz_yyzz, g_z_x_yz_yzzz, g_z_x_yz_zzzz, g_z_x_z_xxxx, g_z_x_z_xxxxy, g_z_x_z_xxxy, g_z_x_z_xxxyy, g_z_x_z_xxxyz, g_z_x_z_xxxz, g_z_x_z_xxyy, g_z_x_z_xxyyy, g_z_x_z_xxyyz, g_z_x_z_xxyz, g_z_x_z_xxyzz, g_z_x_z_xxzz, g_z_x_z_xyyy, g_z_x_z_xyyyy, g_z_x_z_xyyyz, g_z_x_z_xyyz, g_z_x_z_xyyzz, g_z_x_z_xyzz, g_z_x_z_xyzzz, g_z_x_z_xzzz, g_z_x_z_yyyy, g_z_x_z_yyyyy, g_z_x_z_yyyyz, g_z_x_z_yyyz, g_z_x_z_yyyzz, g_z_x_z_yyzz, g_z_x_z_yyzzz, g_z_x_z_yzzz, g_z_x_z_yzzzz, g_z_x_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_xxxx[k] = -g_z_x_z_xxxx[k] * ab_y + g_z_x_z_xxxxy[k];

                g_z_x_yz_xxxy[k] = -g_z_x_z_xxxy[k] * ab_y + g_z_x_z_xxxyy[k];

                g_z_x_yz_xxxz[k] = -g_z_x_z_xxxz[k] * ab_y + g_z_x_z_xxxyz[k];

                g_z_x_yz_xxyy[k] = -g_z_x_z_xxyy[k] * ab_y + g_z_x_z_xxyyy[k];

                g_z_x_yz_xxyz[k] = -g_z_x_z_xxyz[k] * ab_y + g_z_x_z_xxyyz[k];

                g_z_x_yz_xxzz[k] = -g_z_x_z_xxzz[k] * ab_y + g_z_x_z_xxyzz[k];

                g_z_x_yz_xyyy[k] = -g_z_x_z_xyyy[k] * ab_y + g_z_x_z_xyyyy[k];

                g_z_x_yz_xyyz[k] = -g_z_x_z_xyyz[k] * ab_y + g_z_x_z_xyyyz[k];

                g_z_x_yz_xyzz[k] = -g_z_x_z_xyzz[k] * ab_y + g_z_x_z_xyyzz[k];

                g_z_x_yz_xzzz[k] = -g_z_x_z_xzzz[k] * ab_y + g_z_x_z_xyzzz[k];

                g_z_x_yz_yyyy[k] = -g_z_x_z_yyyy[k] * ab_y + g_z_x_z_yyyyy[k];

                g_z_x_yz_yyyz[k] = -g_z_x_z_yyyz[k] * ab_y + g_z_x_z_yyyyz[k];

                g_z_x_yz_yyzz[k] = -g_z_x_z_yyzz[k] * ab_y + g_z_x_z_yyyzz[k];

                g_z_x_yz_yzzz[k] = -g_z_x_z_yzzz[k] * ab_y + g_z_x_z_yyzzz[k];

                g_z_x_yz_zzzz[k] = -g_z_x_z_zzzz[k] * ab_y + g_z_x_z_yzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 615 * ccomps * dcomps);

            auto g_z_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 617 * ccomps * dcomps);

            auto g_z_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 619 * ccomps * dcomps);

            auto g_z_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 623 * ccomps * dcomps);

            auto g_z_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxxx, g_0_x_z_xxxy, g_0_x_z_xxxz, g_0_x_z_xxyy, g_0_x_z_xxyz, g_0_x_z_xxzz, g_0_x_z_xyyy, g_0_x_z_xyyz, g_0_x_z_xyzz, g_0_x_z_xzzz, g_0_x_z_yyyy, g_0_x_z_yyyz, g_0_x_z_yyzz, g_0_x_z_yzzz, g_0_x_z_zzzz, g_z_x_z_xxxx, g_z_x_z_xxxxz, g_z_x_z_xxxy, g_z_x_z_xxxyz, g_z_x_z_xxxz, g_z_x_z_xxxzz, g_z_x_z_xxyy, g_z_x_z_xxyyz, g_z_x_z_xxyz, g_z_x_z_xxyzz, g_z_x_z_xxzz, g_z_x_z_xxzzz, g_z_x_z_xyyy, g_z_x_z_xyyyz, g_z_x_z_xyyz, g_z_x_z_xyyzz, g_z_x_z_xyzz, g_z_x_z_xyzzz, g_z_x_z_xzzz, g_z_x_z_xzzzz, g_z_x_z_yyyy, g_z_x_z_yyyyz, g_z_x_z_yyyz, g_z_x_z_yyyzz, g_z_x_z_yyzz, g_z_x_z_yyzzz, g_z_x_z_yzzz, g_z_x_z_yzzzz, g_z_x_z_zzzz, g_z_x_z_zzzzz, g_z_x_zz_xxxx, g_z_x_zz_xxxy, g_z_x_zz_xxxz, g_z_x_zz_xxyy, g_z_x_zz_xxyz, g_z_x_zz_xxzz, g_z_x_zz_xyyy, g_z_x_zz_xyyz, g_z_x_zz_xyzz, g_z_x_zz_xzzz, g_z_x_zz_yyyy, g_z_x_zz_yyyz, g_z_x_zz_yyzz, g_z_x_zz_yzzz, g_z_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_xxxx[k] = -g_0_x_z_xxxx[k] - g_z_x_z_xxxx[k] * ab_z + g_z_x_z_xxxxz[k];

                g_z_x_zz_xxxy[k] = -g_0_x_z_xxxy[k] - g_z_x_z_xxxy[k] * ab_z + g_z_x_z_xxxyz[k];

                g_z_x_zz_xxxz[k] = -g_0_x_z_xxxz[k] - g_z_x_z_xxxz[k] * ab_z + g_z_x_z_xxxzz[k];

                g_z_x_zz_xxyy[k] = -g_0_x_z_xxyy[k] - g_z_x_z_xxyy[k] * ab_z + g_z_x_z_xxyyz[k];

                g_z_x_zz_xxyz[k] = -g_0_x_z_xxyz[k] - g_z_x_z_xxyz[k] * ab_z + g_z_x_z_xxyzz[k];

                g_z_x_zz_xxzz[k] = -g_0_x_z_xxzz[k] - g_z_x_z_xxzz[k] * ab_z + g_z_x_z_xxzzz[k];

                g_z_x_zz_xyyy[k] = -g_0_x_z_xyyy[k] - g_z_x_z_xyyy[k] * ab_z + g_z_x_z_xyyyz[k];

                g_z_x_zz_xyyz[k] = -g_0_x_z_xyyz[k] - g_z_x_z_xyyz[k] * ab_z + g_z_x_z_xyyzz[k];

                g_z_x_zz_xyzz[k] = -g_0_x_z_xyzz[k] - g_z_x_z_xyzz[k] * ab_z + g_z_x_z_xyzzz[k];

                g_z_x_zz_xzzz[k] = -g_0_x_z_xzzz[k] - g_z_x_z_xzzz[k] * ab_z + g_z_x_z_xzzzz[k];

                g_z_x_zz_yyyy[k] = -g_0_x_z_yyyy[k] - g_z_x_z_yyyy[k] * ab_z + g_z_x_z_yyyyz[k];

                g_z_x_zz_yyyz[k] = -g_0_x_z_yyyz[k] - g_z_x_z_yyyz[k] * ab_z + g_z_x_z_yyyzz[k];

                g_z_x_zz_yyzz[k] = -g_0_x_z_yyzz[k] - g_z_x_z_yyzz[k] * ab_z + g_z_x_z_yyzzz[k];

                g_z_x_zz_yzzz[k] = -g_0_x_z_yzzz[k] - g_z_x_z_yzzz[k] * ab_z + g_z_x_z_yzzzz[k];

                g_z_x_zz_zzzz[k] = -g_0_x_z_zzzz[k] - g_z_x_z_zzzz[k] * ab_z + g_z_x_z_zzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 635 * ccomps * dcomps);

            auto g_z_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 639 * ccomps * dcomps);

            auto g_z_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 641 * ccomps * dcomps);

            auto g_z_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 643 * ccomps * dcomps);

            auto g_z_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_xxxx, g_z_y_x_xxxxx, g_z_y_x_xxxxy, g_z_y_x_xxxxz, g_z_y_x_xxxy, g_z_y_x_xxxyy, g_z_y_x_xxxyz, g_z_y_x_xxxz, g_z_y_x_xxxzz, g_z_y_x_xxyy, g_z_y_x_xxyyy, g_z_y_x_xxyyz, g_z_y_x_xxyz, g_z_y_x_xxyzz, g_z_y_x_xxzz, g_z_y_x_xxzzz, g_z_y_x_xyyy, g_z_y_x_xyyyy, g_z_y_x_xyyyz, g_z_y_x_xyyz, g_z_y_x_xyyzz, g_z_y_x_xyzz, g_z_y_x_xyzzz, g_z_y_x_xzzz, g_z_y_x_xzzzz, g_z_y_x_yyyy, g_z_y_x_yyyz, g_z_y_x_yyzz, g_z_y_x_yzzz, g_z_y_x_zzzz, g_z_y_xx_xxxx, g_z_y_xx_xxxy, g_z_y_xx_xxxz, g_z_y_xx_xxyy, g_z_y_xx_xxyz, g_z_y_xx_xxzz, g_z_y_xx_xyyy, g_z_y_xx_xyyz, g_z_y_xx_xyzz, g_z_y_xx_xzzz, g_z_y_xx_yyyy, g_z_y_xx_yyyz, g_z_y_xx_yyzz, g_z_y_xx_yzzz, g_z_y_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_xxxx[k] = -g_z_y_x_xxxx[k] * ab_x + g_z_y_x_xxxxx[k];

                g_z_y_xx_xxxy[k] = -g_z_y_x_xxxy[k] * ab_x + g_z_y_x_xxxxy[k];

                g_z_y_xx_xxxz[k] = -g_z_y_x_xxxz[k] * ab_x + g_z_y_x_xxxxz[k];

                g_z_y_xx_xxyy[k] = -g_z_y_x_xxyy[k] * ab_x + g_z_y_x_xxxyy[k];

                g_z_y_xx_xxyz[k] = -g_z_y_x_xxyz[k] * ab_x + g_z_y_x_xxxyz[k];

                g_z_y_xx_xxzz[k] = -g_z_y_x_xxzz[k] * ab_x + g_z_y_x_xxxzz[k];

                g_z_y_xx_xyyy[k] = -g_z_y_x_xyyy[k] * ab_x + g_z_y_x_xxyyy[k];

                g_z_y_xx_xyyz[k] = -g_z_y_x_xyyz[k] * ab_x + g_z_y_x_xxyyz[k];

                g_z_y_xx_xyzz[k] = -g_z_y_x_xyzz[k] * ab_x + g_z_y_x_xxyzz[k];

                g_z_y_xx_xzzz[k] = -g_z_y_x_xzzz[k] * ab_x + g_z_y_x_xxzzz[k];

                g_z_y_xx_yyyy[k] = -g_z_y_x_yyyy[k] * ab_x + g_z_y_x_xyyyy[k];

                g_z_y_xx_yyyz[k] = -g_z_y_x_yyyz[k] * ab_x + g_z_y_x_xyyyz[k];

                g_z_y_xx_yyzz[k] = -g_z_y_x_yyzz[k] * ab_x + g_z_y_x_xyyzz[k];

                g_z_y_xx_yzzz[k] = -g_z_y_x_yzzz[k] * ab_x + g_z_y_x_xyzzz[k];

                g_z_y_xx_zzzz[k] = -g_z_y_x_zzzz[k] * ab_x + g_z_y_x_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_xxxx, g_z_y_xy_xxxy, g_z_y_xy_xxxz, g_z_y_xy_xxyy, g_z_y_xy_xxyz, g_z_y_xy_xxzz, g_z_y_xy_xyyy, g_z_y_xy_xyyz, g_z_y_xy_xyzz, g_z_y_xy_xzzz, g_z_y_xy_yyyy, g_z_y_xy_yyyz, g_z_y_xy_yyzz, g_z_y_xy_yzzz, g_z_y_xy_zzzz, g_z_y_y_xxxx, g_z_y_y_xxxxx, g_z_y_y_xxxxy, g_z_y_y_xxxxz, g_z_y_y_xxxy, g_z_y_y_xxxyy, g_z_y_y_xxxyz, g_z_y_y_xxxz, g_z_y_y_xxxzz, g_z_y_y_xxyy, g_z_y_y_xxyyy, g_z_y_y_xxyyz, g_z_y_y_xxyz, g_z_y_y_xxyzz, g_z_y_y_xxzz, g_z_y_y_xxzzz, g_z_y_y_xyyy, g_z_y_y_xyyyy, g_z_y_y_xyyyz, g_z_y_y_xyyz, g_z_y_y_xyyzz, g_z_y_y_xyzz, g_z_y_y_xyzzz, g_z_y_y_xzzz, g_z_y_y_xzzzz, g_z_y_y_yyyy, g_z_y_y_yyyz, g_z_y_y_yyzz, g_z_y_y_yzzz, g_z_y_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_xxxx[k] = -g_z_y_y_xxxx[k] * ab_x + g_z_y_y_xxxxx[k];

                g_z_y_xy_xxxy[k] = -g_z_y_y_xxxy[k] * ab_x + g_z_y_y_xxxxy[k];

                g_z_y_xy_xxxz[k] = -g_z_y_y_xxxz[k] * ab_x + g_z_y_y_xxxxz[k];

                g_z_y_xy_xxyy[k] = -g_z_y_y_xxyy[k] * ab_x + g_z_y_y_xxxyy[k];

                g_z_y_xy_xxyz[k] = -g_z_y_y_xxyz[k] * ab_x + g_z_y_y_xxxyz[k];

                g_z_y_xy_xxzz[k] = -g_z_y_y_xxzz[k] * ab_x + g_z_y_y_xxxzz[k];

                g_z_y_xy_xyyy[k] = -g_z_y_y_xyyy[k] * ab_x + g_z_y_y_xxyyy[k];

                g_z_y_xy_xyyz[k] = -g_z_y_y_xyyz[k] * ab_x + g_z_y_y_xxyyz[k];

                g_z_y_xy_xyzz[k] = -g_z_y_y_xyzz[k] * ab_x + g_z_y_y_xxyzz[k];

                g_z_y_xy_xzzz[k] = -g_z_y_y_xzzz[k] * ab_x + g_z_y_y_xxzzz[k];

                g_z_y_xy_yyyy[k] = -g_z_y_y_yyyy[k] * ab_x + g_z_y_y_xyyyy[k];

                g_z_y_xy_yyyz[k] = -g_z_y_y_yyyz[k] * ab_x + g_z_y_y_xyyyz[k];

                g_z_y_xy_yyzz[k] = -g_z_y_y_yyzz[k] * ab_x + g_z_y_y_xyyzz[k];

                g_z_y_xy_yzzz[k] = -g_z_y_y_yzzz[k] * ab_x + g_z_y_y_xyzzz[k];

                g_z_y_xy_zzzz[k] = -g_z_y_y_zzzz[k] * ab_x + g_z_y_y_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 671 * ccomps * dcomps);

            auto g_z_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_xxxx, g_z_y_xz_xxxy, g_z_y_xz_xxxz, g_z_y_xz_xxyy, g_z_y_xz_xxyz, g_z_y_xz_xxzz, g_z_y_xz_xyyy, g_z_y_xz_xyyz, g_z_y_xz_xyzz, g_z_y_xz_xzzz, g_z_y_xz_yyyy, g_z_y_xz_yyyz, g_z_y_xz_yyzz, g_z_y_xz_yzzz, g_z_y_xz_zzzz, g_z_y_z_xxxx, g_z_y_z_xxxxx, g_z_y_z_xxxxy, g_z_y_z_xxxxz, g_z_y_z_xxxy, g_z_y_z_xxxyy, g_z_y_z_xxxyz, g_z_y_z_xxxz, g_z_y_z_xxxzz, g_z_y_z_xxyy, g_z_y_z_xxyyy, g_z_y_z_xxyyz, g_z_y_z_xxyz, g_z_y_z_xxyzz, g_z_y_z_xxzz, g_z_y_z_xxzzz, g_z_y_z_xyyy, g_z_y_z_xyyyy, g_z_y_z_xyyyz, g_z_y_z_xyyz, g_z_y_z_xyyzz, g_z_y_z_xyzz, g_z_y_z_xyzzz, g_z_y_z_xzzz, g_z_y_z_xzzzz, g_z_y_z_yyyy, g_z_y_z_yyyz, g_z_y_z_yyzz, g_z_y_z_yzzz, g_z_y_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_xxxx[k] = -g_z_y_z_xxxx[k] * ab_x + g_z_y_z_xxxxx[k];

                g_z_y_xz_xxxy[k] = -g_z_y_z_xxxy[k] * ab_x + g_z_y_z_xxxxy[k];

                g_z_y_xz_xxxz[k] = -g_z_y_z_xxxz[k] * ab_x + g_z_y_z_xxxxz[k];

                g_z_y_xz_xxyy[k] = -g_z_y_z_xxyy[k] * ab_x + g_z_y_z_xxxyy[k];

                g_z_y_xz_xxyz[k] = -g_z_y_z_xxyz[k] * ab_x + g_z_y_z_xxxyz[k];

                g_z_y_xz_xxzz[k] = -g_z_y_z_xxzz[k] * ab_x + g_z_y_z_xxxzz[k];

                g_z_y_xz_xyyy[k] = -g_z_y_z_xyyy[k] * ab_x + g_z_y_z_xxyyy[k];

                g_z_y_xz_xyyz[k] = -g_z_y_z_xyyz[k] * ab_x + g_z_y_z_xxyyz[k];

                g_z_y_xz_xyzz[k] = -g_z_y_z_xyzz[k] * ab_x + g_z_y_z_xxyzz[k];

                g_z_y_xz_xzzz[k] = -g_z_y_z_xzzz[k] * ab_x + g_z_y_z_xxzzz[k];

                g_z_y_xz_yyyy[k] = -g_z_y_z_yyyy[k] * ab_x + g_z_y_z_xyyyy[k];

                g_z_y_xz_yyyz[k] = -g_z_y_z_yyyz[k] * ab_x + g_z_y_z_xyyyz[k];

                g_z_y_xz_yyzz[k] = -g_z_y_z_yyzz[k] * ab_x + g_z_y_z_xyyzz[k];

                g_z_y_xz_yzzz[k] = -g_z_y_z_yzzz[k] * ab_x + g_z_y_z_xyzzz[k];

                g_z_y_xz_zzzz[k] = -g_z_y_z_zzzz[k] * ab_x + g_z_y_z_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz, g_z_y_y_xxxx, g_z_y_y_xxxxy, g_z_y_y_xxxy, g_z_y_y_xxxyy, g_z_y_y_xxxyz, g_z_y_y_xxxz, g_z_y_y_xxyy, g_z_y_y_xxyyy, g_z_y_y_xxyyz, g_z_y_y_xxyz, g_z_y_y_xxyzz, g_z_y_y_xxzz, g_z_y_y_xyyy, g_z_y_y_xyyyy, g_z_y_y_xyyyz, g_z_y_y_xyyz, g_z_y_y_xyyzz, g_z_y_y_xyzz, g_z_y_y_xyzzz, g_z_y_y_xzzz, g_z_y_y_yyyy, g_z_y_y_yyyyy, g_z_y_y_yyyyz, g_z_y_y_yyyz, g_z_y_y_yyyzz, g_z_y_y_yyzz, g_z_y_y_yyzzz, g_z_y_y_yzzz, g_z_y_y_yzzzz, g_z_y_y_zzzz, g_z_y_yy_xxxx, g_z_y_yy_xxxy, g_z_y_yy_xxxz, g_z_y_yy_xxyy, g_z_y_yy_xxyz, g_z_y_yy_xxzz, g_z_y_yy_xyyy, g_z_y_yy_xyyz, g_z_y_yy_xyzz, g_z_y_yy_xzzz, g_z_y_yy_yyyy, g_z_y_yy_yyyz, g_z_y_yy_yyzz, g_z_y_yy_yzzz, g_z_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_xxxx[k] = g_z_0_y_xxxx[k] - g_z_y_y_xxxx[k] * ab_y + g_z_y_y_xxxxy[k];

                g_z_y_yy_xxxy[k] = g_z_0_y_xxxy[k] - g_z_y_y_xxxy[k] * ab_y + g_z_y_y_xxxyy[k];

                g_z_y_yy_xxxz[k] = g_z_0_y_xxxz[k] - g_z_y_y_xxxz[k] * ab_y + g_z_y_y_xxxyz[k];

                g_z_y_yy_xxyy[k] = g_z_0_y_xxyy[k] - g_z_y_y_xxyy[k] * ab_y + g_z_y_y_xxyyy[k];

                g_z_y_yy_xxyz[k] = g_z_0_y_xxyz[k] - g_z_y_y_xxyz[k] * ab_y + g_z_y_y_xxyyz[k];

                g_z_y_yy_xxzz[k] = g_z_0_y_xxzz[k] - g_z_y_y_xxzz[k] * ab_y + g_z_y_y_xxyzz[k];

                g_z_y_yy_xyyy[k] = g_z_0_y_xyyy[k] - g_z_y_y_xyyy[k] * ab_y + g_z_y_y_xyyyy[k];

                g_z_y_yy_xyyz[k] = g_z_0_y_xyyz[k] - g_z_y_y_xyyz[k] * ab_y + g_z_y_y_xyyyz[k];

                g_z_y_yy_xyzz[k] = g_z_0_y_xyzz[k] - g_z_y_y_xyzz[k] * ab_y + g_z_y_y_xyyzz[k];

                g_z_y_yy_xzzz[k] = g_z_0_y_xzzz[k] - g_z_y_y_xzzz[k] * ab_y + g_z_y_y_xyzzz[k];

                g_z_y_yy_yyyy[k] = g_z_0_y_yyyy[k] - g_z_y_y_yyyy[k] * ab_y + g_z_y_y_yyyyy[k];

                g_z_y_yy_yyyz[k] = g_z_0_y_yyyz[k] - g_z_y_y_yyyz[k] * ab_y + g_z_y_y_yyyyz[k];

                g_z_y_yy_yyzz[k] = g_z_0_y_yyzz[k] - g_z_y_y_yyzz[k] * ab_y + g_z_y_y_yyyzz[k];

                g_z_y_yy_yzzz[k] = g_z_0_y_yzzz[k] - g_z_y_y_yzzz[k] * ab_y + g_z_y_y_yyzzz[k];

                g_z_y_yy_zzzz[k] = g_z_0_y_zzzz[k] - g_z_y_y_zzzz[k] * ab_y + g_z_y_y_yzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 699 * ccomps * dcomps);

            auto g_z_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz, g_z_y_yz_xxxx, g_z_y_yz_xxxy, g_z_y_yz_xxxz, g_z_y_yz_xxyy, g_z_y_yz_xxyz, g_z_y_yz_xxzz, g_z_y_yz_xyyy, g_z_y_yz_xyyz, g_z_y_yz_xyzz, g_z_y_yz_xzzz, g_z_y_yz_yyyy, g_z_y_yz_yyyz, g_z_y_yz_yyzz, g_z_y_yz_yzzz, g_z_y_yz_zzzz, g_z_y_z_xxxx, g_z_y_z_xxxxy, g_z_y_z_xxxy, g_z_y_z_xxxyy, g_z_y_z_xxxyz, g_z_y_z_xxxz, g_z_y_z_xxyy, g_z_y_z_xxyyy, g_z_y_z_xxyyz, g_z_y_z_xxyz, g_z_y_z_xxyzz, g_z_y_z_xxzz, g_z_y_z_xyyy, g_z_y_z_xyyyy, g_z_y_z_xyyyz, g_z_y_z_xyyz, g_z_y_z_xyyzz, g_z_y_z_xyzz, g_z_y_z_xyzzz, g_z_y_z_xzzz, g_z_y_z_yyyy, g_z_y_z_yyyyy, g_z_y_z_yyyyz, g_z_y_z_yyyz, g_z_y_z_yyyzz, g_z_y_z_yyzz, g_z_y_z_yyzzz, g_z_y_z_yzzz, g_z_y_z_yzzzz, g_z_y_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_xxxx[k] = g_z_0_z_xxxx[k] - g_z_y_z_xxxx[k] * ab_y + g_z_y_z_xxxxy[k];

                g_z_y_yz_xxxy[k] = g_z_0_z_xxxy[k] - g_z_y_z_xxxy[k] * ab_y + g_z_y_z_xxxyy[k];

                g_z_y_yz_xxxz[k] = g_z_0_z_xxxz[k] - g_z_y_z_xxxz[k] * ab_y + g_z_y_z_xxxyz[k];

                g_z_y_yz_xxyy[k] = g_z_0_z_xxyy[k] - g_z_y_z_xxyy[k] * ab_y + g_z_y_z_xxyyy[k];

                g_z_y_yz_xxyz[k] = g_z_0_z_xxyz[k] - g_z_y_z_xxyz[k] * ab_y + g_z_y_z_xxyyz[k];

                g_z_y_yz_xxzz[k] = g_z_0_z_xxzz[k] - g_z_y_z_xxzz[k] * ab_y + g_z_y_z_xxyzz[k];

                g_z_y_yz_xyyy[k] = g_z_0_z_xyyy[k] - g_z_y_z_xyyy[k] * ab_y + g_z_y_z_xyyyy[k];

                g_z_y_yz_xyyz[k] = g_z_0_z_xyyz[k] - g_z_y_z_xyyz[k] * ab_y + g_z_y_z_xyyyz[k];

                g_z_y_yz_xyzz[k] = g_z_0_z_xyzz[k] - g_z_y_z_xyzz[k] * ab_y + g_z_y_z_xyyzz[k];

                g_z_y_yz_xzzz[k] = g_z_0_z_xzzz[k] - g_z_y_z_xzzz[k] * ab_y + g_z_y_z_xyzzz[k];

                g_z_y_yz_yyyy[k] = g_z_0_z_yyyy[k] - g_z_y_z_yyyy[k] * ab_y + g_z_y_z_yyyyy[k];

                g_z_y_yz_yyyz[k] = g_z_0_z_yyyz[k] - g_z_y_z_yyyz[k] * ab_y + g_z_y_z_yyyyz[k];

                g_z_y_yz_yyzz[k] = g_z_0_z_yyzz[k] - g_z_y_z_yyzz[k] * ab_y + g_z_y_z_yyyzz[k];

                g_z_y_yz_yzzz[k] = g_z_0_z_yzzz[k] - g_z_y_z_yzzz[k] * ab_y + g_z_y_z_yyzzz[k];

                g_z_y_yz_zzzz[k] = g_z_0_z_zzzz[k] - g_z_y_z_zzzz[k] * ab_y + g_z_y_z_yzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxxx, g_0_y_z_xxxy, g_0_y_z_xxxz, g_0_y_z_xxyy, g_0_y_z_xxyz, g_0_y_z_xxzz, g_0_y_z_xyyy, g_0_y_z_xyyz, g_0_y_z_xyzz, g_0_y_z_xzzz, g_0_y_z_yyyy, g_0_y_z_yyyz, g_0_y_z_yyzz, g_0_y_z_yzzz, g_0_y_z_zzzz, g_z_y_z_xxxx, g_z_y_z_xxxxz, g_z_y_z_xxxy, g_z_y_z_xxxyz, g_z_y_z_xxxz, g_z_y_z_xxxzz, g_z_y_z_xxyy, g_z_y_z_xxyyz, g_z_y_z_xxyz, g_z_y_z_xxyzz, g_z_y_z_xxzz, g_z_y_z_xxzzz, g_z_y_z_xyyy, g_z_y_z_xyyyz, g_z_y_z_xyyz, g_z_y_z_xyyzz, g_z_y_z_xyzz, g_z_y_z_xyzzz, g_z_y_z_xzzz, g_z_y_z_xzzzz, g_z_y_z_yyyy, g_z_y_z_yyyyz, g_z_y_z_yyyz, g_z_y_z_yyyzz, g_z_y_z_yyzz, g_z_y_z_yyzzz, g_z_y_z_yzzz, g_z_y_z_yzzzz, g_z_y_z_zzzz, g_z_y_z_zzzzz, g_z_y_zz_xxxx, g_z_y_zz_xxxy, g_z_y_zz_xxxz, g_z_y_zz_xxyy, g_z_y_zz_xxyz, g_z_y_zz_xxzz, g_z_y_zz_xyyy, g_z_y_zz_xyyz, g_z_y_zz_xyzz, g_z_y_zz_xzzz, g_z_y_zz_yyyy, g_z_y_zz_yyyz, g_z_y_zz_yyzz, g_z_y_zz_yzzz, g_z_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_xxxx[k] = -g_0_y_z_xxxx[k] - g_z_y_z_xxxx[k] * ab_z + g_z_y_z_xxxxz[k];

                g_z_y_zz_xxxy[k] = -g_0_y_z_xxxy[k] - g_z_y_z_xxxy[k] * ab_z + g_z_y_z_xxxyz[k];

                g_z_y_zz_xxxz[k] = -g_0_y_z_xxxz[k] - g_z_y_z_xxxz[k] * ab_z + g_z_y_z_xxxzz[k];

                g_z_y_zz_xxyy[k] = -g_0_y_z_xxyy[k] - g_z_y_z_xxyy[k] * ab_z + g_z_y_z_xxyyz[k];

                g_z_y_zz_xxyz[k] = -g_0_y_z_xxyz[k] - g_z_y_z_xxyz[k] * ab_z + g_z_y_z_xxyzz[k];

                g_z_y_zz_xxzz[k] = -g_0_y_z_xxzz[k] - g_z_y_z_xxzz[k] * ab_z + g_z_y_z_xxzzz[k];

                g_z_y_zz_xyyy[k] = -g_0_y_z_xyyy[k] - g_z_y_z_xyyy[k] * ab_z + g_z_y_z_xyyyz[k];

                g_z_y_zz_xyyz[k] = -g_0_y_z_xyyz[k] - g_z_y_z_xyyz[k] * ab_z + g_z_y_z_xyyzz[k];

                g_z_y_zz_xyzz[k] = -g_0_y_z_xyzz[k] - g_z_y_z_xyzz[k] * ab_z + g_z_y_z_xyzzz[k];

                g_z_y_zz_xzzz[k] = -g_0_y_z_xzzz[k] - g_z_y_z_xzzz[k] * ab_z + g_z_y_z_xzzzz[k];

                g_z_y_zz_yyyy[k] = -g_0_y_z_yyyy[k] - g_z_y_z_yyyy[k] * ab_z + g_z_y_z_yyyyz[k];

                g_z_y_zz_yyyz[k] = -g_0_y_z_yyyz[k] - g_z_y_z_yyyz[k] * ab_z + g_z_y_z_yyyzz[k];

                g_z_y_zz_yyzz[k] = -g_0_y_z_yyzz[k] - g_z_y_z_yyzz[k] * ab_z + g_z_y_z_yyzzz[k];

                g_z_y_zz_yzzz[k] = -g_0_y_z_yzzz[k] - g_z_y_z_yzzz[k] * ab_z + g_z_y_z_yzzzz[k];

                g_z_y_zz_zzzz[k] = -g_0_y_z_zzzz[k] - g_z_y_z_zzzz[k] * ab_z + g_z_y_z_zzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_xxxx, g_z_z_x_xxxxx, g_z_z_x_xxxxy, g_z_z_x_xxxxz, g_z_z_x_xxxy, g_z_z_x_xxxyy, g_z_z_x_xxxyz, g_z_z_x_xxxz, g_z_z_x_xxxzz, g_z_z_x_xxyy, g_z_z_x_xxyyy, g_z_z_x_xxyyz, g_z_z_x_xxyz, g_z_z_x_xxyzz, g_z_z_x_xxzz, g_z_z_x_xxzzz, g_z_z_x_xyyy, g_z_z_x_xyyyy, g_z_z_x_xyyyz, g_z_z_x_xyyz, g_z_z_x_xyyzz, g_z_z_x_xyzz, g_z_z_x_xyzzz, g_z_z_x_xzzz, g_z_z_x_xzzzz, g_z_z_x_yyyy, g_z_z_x_yyyz, g_z_z_x_yyzz, g_z_z_x_yzzz, g_z_z_x_zzzz, g_z_z_xx_xxxx, g_z_z_xx_xxxy, g_z_z_xx_xxxz, g_z_z_xx_xxyy, g_z_z_xx_xxyz, g_z_z_xx_xxzz, g_z_z_xx_xyyy, g_z_z_xx_xyyz, g_z_z_xx_xyzz, g_z_z_xx_xzzz, g_z_z_xx_yyyy, g_z_z_xx_yyyz, g_z_z_xx_yyzz, g_z_z_xx_yzzz, g_z_z_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_xxxx[k] = -g_z_z_x_xxxx[k] * ab_x + g_z_z_x_xxxxx[k];

                g_z_z_xx_xxxy[k] = -g_z_z_x_xxxy[k] * ab_x + g_z_z_x_xxxxy[k];

                g_z_z_xx_xxxz[k] = -g_z_z_x_xxxz[k] * ab_x + g_z_z_x_xxxxz[k];

                g_z_z_xx_xxyy[k] = -g_z_z_x_xxyy[k] * ab_x + g_z_z_x_xxxyy[k];

                g_z_z_xx_xxyz[k] = -g_z_z_x_xxyz[k] * ab_x + g_z_z_x_xxxyz[k];

                g_z_z_xx_xxzz[k] = -g_z_z_x_xxzz[k] * ab_x + g_z_z_x_xxxzz[k];

                g_z_z_xx_xyyy[k] = -g_z_z_x_xyyy[k] * ab_x + g_z_z_x_xxyyy[k];

                g_z_z_xx_xyyz[k] = -g_z_z_x_xyyz[k] * ab_x + g_z_z_x_xxyyz[k];

                g_z_z_xx_xyzz[k] = -g_z_z_x_xyzz[k] * ab_x + g_z_z_x_xxyzz[k];

                g_z_z_xx_xzzz[k] = -g_z_z_x_xzzz[k] * ab_x + g_z_z_x_xxzzz[k];

                g_z_z_xx_yyyy[k] = -g_z_z_x_yyyy[k] * ab_x + g_z_z_x_xyyyy[k];

                g_z_z_xx_yyyz[k] = -g_z_z_x_yyyz[k] * ab_x + g_z_z_x_xyyyz[k];

                g_z_z_xx_yyzz[k] = -g_z_z_x_yyzz[k] * ab_x + g_z_z_x_xyyzz[k];

                g_z_z_xx_yzzz[k] = -g_z_z_x_yzzz[k] * ab_x + g_z_z_x_xyzzz[k];

                g_z_z_xx_zzzz[k] = -g_z_z_x_zzzz[k] * ab_x + g_z_z_x_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_xxxx, g_z_z_xy_xxxy, g_z_z_xy_xxxz, g_z_z_xy_xxyy, g_z_z_xy_xxyz, g_z_z_xy_xxzz, g_z_z_xy_xyyy, g_z_z_xy_xyyz, g_z_z_xy_xyzz, g_z_z_xy_xzzz, g_z_z_xy_yyyy, g_z_z_xy_yyyz, g_z_z_xy_yyzz, g_z_z_xy_yzzz, g_z_z_xy_zzzz, g_z_z_y_xxxx, g_z_z_y_xxxxx, g_z_z_y_xxxxy, g_z_z_y_xxxxz, g_z_z_y_xxxy, g_z_z_y_xxxyy, g_z_z_y_xxxyz, g_z_z_y_xxxz, g_z_z_y_xxxzz, g_z_z_y_xxyy, g_z_z_y_xxyyy, g_z_z_y_xxyyz, g_z_z_y_xxyz, g_z_z_y_xxyzz, g_z_z_y_xxzz, g_z_z_y_xxzzz, g_z_z_y_xyyy, g_z_z_y_xyyyy, g_z_z_y_xyyyz, g_z_z_y_xyyz, g_z_z_y_xyyzz, g_z_z_y_xyzz, g_z_z_y_xyzzz, g_z_z_y_xzzz, g_z_z_y_xzzzz, g_z_z_y_yyyy, g_z_z_y_yyyz, g_z_z_y_yyzz, g_z_z_y_yzzz, g_z_z_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_xxxx[k] = -g_z_z_y_xxxx[k] * ab_x + g_z_z_y_xxxxx[k];

                g_z_z_xy_xxxy[k] = -g_z_z_y_xxxy[k] * ab_x + g_z_z_y_xxxxy[k];

                g_z_z_xy_xxxz[k] = -g_z_z_y_xxxz[k] * ab_x + g_z_z_y_xxxxz[k];

                g_z_z_xy_xxyy[k] = -g_z_z_y_xxyy[k] * ab_x + g_z_z_y_xxxyy[k];

                g_z_z_xy_xxyz[k] = -g_z_z_y_xxyz[k] * ab_x + g_z_z_y_xxxyz[k];

                g_z_z_xy_xxzz[k] = -g_z_z_y_xxzz[k] * ab_x + g_z_z_y_xxxzz[k];

                g_z_z_xy_xyyy[k] = -g_z_z_y_xyyy[k] * ab_x + g_z_z_y_xxyyy[k];

                g_z_z_xy_xyyz[k] = -g_z_z_y_xyyz[k] * ab_x + g_z_z_y_xxyyz[k];

                g_z_z_xy_xyzz[k] = -g_z_z_y_xyzz[k] * ab_x + g_z_z_y_xxyzz[k];

                g_z_z_xy_xzzz[k] = -g_z_z_y_xzzz[k] * ab_x + g_z_z_y_xxzzz[k];

                g_z_z_xy_yyyy[k] = -g_z_z_y_yyyy[k] * ab_x + g_z_z_y_xyyyy[k];

                g_z_z_xy_yyyz[k] = -g_z_z_y_yyyz[k] * ab_x + g_z_z_y_xyyyz[k];

                g_z_z_xy_yyzz[k] = -g_z_z_y_yyzz[k] * ab_x + g_z_z_y_xyyzz[k];

                g_z_z_xy_yzzz[k] = -g_z_z_y_yzzz[k] * ab_x + g_z_z_y_xyzzz[k];

                g_z_z_xy_zzzz[k] = -g_z_z_y_zzzz[k] * ab_x + g_z_z_y_xzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_xxxx, g_z_z_xz_xxxy, g_z_z_xz_xxxz, g_z_z_xz_xxyy, g_z_z_xz_xxyz, g_z_z_xz_xxzz, g_z_z_xz_xyyy, g_z_z_xz_xyyz, g_z_z_xz_xyzz, g_z_z_xz_xzzz, g_z_z_xz_yyyy, g_z_z_xz_yyyz, g_z_z_xz_yyzz, g_z_z_xz_yzzz, g_z_z_xz_zzzz, g_z_z_z_xxxx, g_z_z_z_xxxxx, g_z_z_z_xxxxy, g_z_z_z_xxxxz, g_z_z_z_xxxy, g_z_z_z_xxxyy, g_z_z_z_xxxyz, g_z_z_z_xxxz, g_z_z_z_xxxzz, g_z_z_z_xxyy, g_z_z_z_xxyyy, g_z_z_z_xxyyz, g_z_z_z_xxyz, g_z_z_z_xxyzz, g_z_z_z_xxzz, g_z_z_z_xxzzz, g_z_z_z_xyyy, g_z_z_z_xyyyy, g_z_z_z_xyyyz, g_z_z_z_xyyz, g_z_z_z_xyyzz, g_z_z_z_xyzz, g_z_z_z_xyzzz, g_z_z_z_xzzz, g_z_z_z_xzzzz, g_z_z_z_yyyy, g_z_z_z_yyyz, g_z_z_z_yyzz, g_z_z_z_yzzz, g_z_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_xxxx[k] = -g_z_z_z_xxxx[k] * ab_x + g_z_z_z_xxxxx[k];

                g_z_z_xz_xxxy[k] = -g_z_z_z_xxxy[k] * ab_x + g_z_z_z_xxxxy[k];

                g_z_z_xz_xxxz[k] = -g_z_z_z_xxxz[k] * ab_x + g_z_z_z_xxxxz[k];

                g_z_z_xz_xxyy[k] = -g_z_z_z_xxyy[k] * ab_x + g_z_z_z_xxxyy[k];

                g_z_z_xz_xxyz[k] = -g_z_z_z_xxyz[k] * ab_x + g_z_z_z_xxxyz[k];

                g_z_z_xz_xxzz[k] = -g_z_z_z_xxzz[k] * ab_x + g_z_z_z_xxxzz[k];

                g_z_z_xz_xyyy[k] = -g_z_z_z_xyyy[k] * ab_x + g_z_z_z_xxyyy[k];

                g_z_z_xz_xyyz[k] = -g_z_z_z_xyyz[k] * ab_x + g_z_z_z_xxyyz[k];

                g_z_z_xz_xyzz[k] = -g_z_z_z_xyzz[k] * ab_x + g_z_z_z_xxyzz[k];

                g_z_z_xz_xzzz[k] = -g_z_z_z_xzzz[k] * ab_x + g_z_z_z_xxzzz[k];

                g_z_z_xz_yyyy[k] = -g_z_z_z_yyyy[k] * ab_x + g_z_z_z_xyyyy[k];

                g_z_z_xz_yyyz[k] = -g_z_z_z_yyyz[k] * ab_x + g_z_z_z_xyyyz[k];

                g_z_z_xz_yyzz[k] = -g_z_z_z_yyzz[k] * ab_x + g_z_z_z_xyyzz[k];

                g_z_z_xz_yzzz[k] = -g_z_z_z_yzzz[k] * ab_x + g_z_z_z_xyzzz[k];

                g_z_z_xz_zzzz[k] = -g_z_z_z_zzzz[k] * ab_x + g_z_z_z_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_xxxx, g_z_z_y_xxxxy, g_z_z_y_xxxy, g_z_z_y_xxxyy, g_z_z_y_xxxyz, g_z_z_y_xxxz, g_z_z_y_xxyy, g_z_z_y_xxyyy, g_z_z_y_xxyyz, g_z_z_y_xxyz, g_z_z_y_xxyzz, g_z_z_y_xxzz, g_z_z_y_xyyy, g_z_z_y_xyyyy, g_z_z_y_xyyyz, g_z_z_y_xyyz, g_z_z_y_xyyzz, g_z_z_y_xyzz, g_z_z_y_xyzzz, g_z_z_y_xzzz, g_z_z_y_yyyy, g_z_z_y_yyyyy, g_z_z_y_yyyyz, g_z_z_y_yyyz, g_z_z_y_yyyzz, g_z_z_y_yyzz, g_z_z_y_yyzzz, g_z_z_y_yzzz, g_z_z_y_yzzzz, g_z_z_y_zzzz, g_z_z_yy_xxxx, g_z_z_yy_xxxy, g_z_z_yy_xxxz, g_z_z_yy_xxyy, g_z_z_yy_xxyz, g_z_z_yy_xxzz, g_z_z_yy_xyyy, g_z_z_yy_xyyz, g_z_z_yy_xyzz, g_z_z_yy_xzzz, g_z_z_yy_yyyy, g_z_z_yy_yyyz, g_z_z_yy_yyzz, g_z_z_yy_yzzz, g_z_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_xxxx[k] = -g_z_z_y_xxxx[k] * ab_y + g_z_z_y_xxxxy[k];

                g_z_z_yy_xxxy[k] = -g_z_z_y_xxxy[k] * ab_y + g_z_z_y_xxxyy[k];

                g_z_z_yy_xxxz[k] = -g_z_z_y_xxxz[k] * ab_y + g_z_z_y_xxxyz[k];

                g_z_z_yy_xxyy[k] = -g_z_z_y_xxyy[k] * ab_y + g_z_z_y_xxyyy[k];

                g_z_z_yy_xxyz[k] = -g_z_z_y_xxyz[k] * ab_y + g_z_z_y_xxyyz[k];

                g_z_z_yy_xxzz[k] = -g_z_z_y_xxzz[k] * ab_y + g_z_z_y_xxyzz[k];

                g_z_z_yy_xyyy[k] = -g_z_z_y_xyyy[k] * ab_y + g_z_z_y_xyyyy[k];

                g_z_z_yy_xyyz[k] = -g_z_z_y_xyyz[k] * ab_y + g_z_z_y_xyyyz[k];

                g_z_z_yy_xyzz[k] = -g_z_z_y_xyzz[k] * ab_y + g_z_z_y_xyyzz[k];

                g_z_z_yy_xzzz[k] = -g_z_z_y_xzzz[k] * ab_y + g_z_z_y_xyzzz[k];

                g_z_z_yy_yyyy[k] = -g_z_z_y_yyyy[k] * ab_y + g_z_z_y_yyyyy[k];

                g_z_z_yy_yyyz[k] = -g_z_z_y_yyyz[k] * ab_y + g_z_z_y_yyyyz[k];

                g_z_z_yy_yyzz[k] = -g_z_z_y_yyzz[k] * ab_y + g_z_z_y_yyyzz[k];

                g_z_z_yy_yzzz[k] = -g_z_z_y_yzzz[k] * ab_y + g_z_z_y_yyzzz[k];

                g_z_z_yy_zzzz[k] = -g_z_z_y_zzzz[k] * ab_y + g_z_z_y_yzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_xxxx, g_z_z_yz_xxxy, g_z_z_yz_xxxz, g_z_z_yz_xxyy, g_z_z_yz_xxyz, g_z_z_yz_xxzz, g_z_z_yz_xyyy, g_z_z_yz_xyyz, g_z_z_yz_xyzz, g_z_z_yz_xzzz, g_z_z_yz_yyyy, g_z_z_yz_yyyz, g_z_z_yz_yyzz, g_z_z_yz_yzzz, g_z_z_yz_zzzz, g_z_z_z_xxxx, g_z_z_z_xxxxy, g_z_z_z_xxxy, g_z_z_z_xxxyy, g_z_z_z_xxxyz, g_z_z_z_xxxz, g_z_z_z_xxyy, g_z_z_z_xxyyy, g_z_z_z_xxyyz, g_z_z_z_xxyz, g_z_z_z_xxyzz, g_z_z_z_xxzz, g_z_z_z_xyyy, g_z_z_z_xyyyy, g_z_z_z_xyyyz, g_z_z_z_xyyz, g_z_z_z_xyyzz, g_z_z_z_xyzz, g_z_z_z_xyzzz, g_z_z_z_xzzz, g_z_z_z_yyyy, g_z_z_z_yyyyy, g_z_z_z_yyyyz, g_z_z_z_yyyz, g_z_z_z_yyyzz, g_z_z_z_yyzz, g_z_z_z_yyzzz, g_z_z_z_yzzz, g_z_z_z_yzzzz, g_z_z_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_xxxx[k] = -g_z_z_z_xxxx[k] * ab_y + g_z_z_z_xxxxy[k];

                g_z_z_yz_xxxy[k] = -g_z_z_z_xxxy[k] * ab_y + g_z_z_z_xxxyy[k];

                g_z_z_yz_xxxz[k] = -g_z_z_z_xxxz[k] * ab_y + g_z_z_z_xxxyz[k];

                g_z_z_yz_xxyy[k] = -g_z_z_z_xxyy[k] * ab_y + g_z_z_z_xxyyy[k];

                g_z_z_yz_xxyz[k] = -g_z_z_z_xxyz[k] * ab_y + g_z_z_z_xxyyz[k];

                g_z_z_yz_xxzz[k] = -g_z_z_z_xxzz[k] * ab_y + g_z_z_z_xxyzz[k];

                g_z_z_yz_xyyy[k] = -g_z_z_z_xyyy[k] * ab_y + g_z_z_z_xyyyy[k];

                g_z_z_yz_xyyz[k] = -g_z_z_z_xyyz[k] * ab_y + g_z_z_z_xyyyz[k];

                g_z_z_yz_xyzz[k] = -g_z_z_z_xyzz[k] * ab_y + g_z_z_z_xyyzz[k];

                g_z_z_yz_xzzz[k] = -g_z_z_z_xzzz[k] * ab_y + g_z_z_z_xyzzz[k];

                g_z_z_yz_yyyy[k] = -g_z_z_z_yyyy[k] * ab_y + g_z_z_z_yyyyy[k];

                g_z_z_yz_yyyz[k] = -g_z_z_z_yyyz[k] * ab_y + g_z_z_z_yyyyz[k];

                g_z_z_yz_yyzz[k] = -g_z_z_z_yyzz[k] * ab_y + g_z_z_z_yyyzz[k];

                g_z_z_yz_yzzz[k] = -g_z_z_z_yzzz[k] * ab_y + g_z_z_z_yyzzz[k];

                g_z_z_yz_zzzz[k] = -g_z_z_z_zzzz[k] * ab_y + g_z_z_z_yzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxx, g_0_z_z_xxxy, g_0_z_z_xxxz, g_0_z_z_xxyy, g_0_z_z_xxyz, g_0_z_z_xxzz, g_0_z_z_xyyy, g_0_z_z_xyyz, g_0_z_z_xyzz, g_0_z_z_xzzz, g_0_z_z_yyyy, g_0_z_z_yyyz, g_0_z_z_yyzz, g_0_z_z_yzzz, g_0_z_z_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz, g_z_z_z_xxxx, g_z_z_z_xxxxz, g_z_z_z_xxxy, g_z_z_z_xxxyz, g_z_z_z_xxxz, g_z_z_z_xxxzz, g_z_z_z_xxyy, g_z_z_z_xxyyz, g_z_z_z_xxyz, g_z_z_z_xxyzz, g_z_z_z_xxzz, g_z_z_z_xxzzz, g_z_z_z_xyyy, g_z_z_z_xyyyz, g_z_z_z_xyyz, g_z_z_z_xyyzz, g_z_z_z_xyzz, g_z_z_z_xyzzz, g_z_z_z_xzzz, g_z_z_z_xzzzz, g_z_z_z_yyyy, g_z_z_z_yyyyz, g_z_z_z_yyyz, g_z_z_z_yyyzz, g_z_z_z_yyzz, g_z_z_z_yyzzz, g_z_z_z_yzzz, g_z_z_z_yzzzz, g_z_z_z_zzzz, g_z_z_z_zzzzz, g_z_z_zz_xxxx, g_z_z_zz_xxxy, g_z_z_zz_xxxz, g_z_z_zz_xxyy, g_z_z_zz_xxyz, g_z_z_zz_xxzz, g_z_z_zz_xyyy, g_z_z_zz_xyyz, g_z_z_zz_xyzz, g_z_z_zz_xzzz, g_z_z_zz_yyyy, g_z_z_zz_yyyz, g_z_z_zz_yyzz, g_z_z_zz_yzzz, g_z_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_xxxx[k] = -g_0_z_z_xxxx[k] + g_z_0_z_xxxx[k] - g_z_z_z_xxxx[k] * ab_z + g_z_z_z_xxxxz[k];

                g_z_z_zz_xxxy[k] = -g_0_z_z_xxxy[k] + g_z_0_z_xxxy[k] - g_z_z_z_xxxy[k] * ab_z + g_z_z_z_xxxyz[k];

                g_z_z_zz_xxxz[k] = -g_0_z_z_xxxz[k] + g_z_0_z_xxxz[k] - g_z_z_z_xxxz[k] * ab_z + g_z_z_z_xxxzz[k];

                g_z_z_zz_xxyy[k] = -g_0_z_z_xxyy[k] + g_z_0_z_xxyy[k] - g_z_z_z_xxyy[k] * ab_z + g_z_z_z_xxyyz[k];

                g_z_z_zz_xxyz[k] = -g_0_z_z_xxyz[k] + g_z_0_z_xxyz[k] - g_z_z_z_xxyz[k] * ab_z + g_z_z_z_xxyzz[k];

                g_z_z_zz_xxzz[k] = -g_0_z_z_xxzz[k] + g_z_0_z_xxzz[k] - g_z_z_z_xxzz[k] * ab_z + g_z_z_z_xxzzz[k];

                g_z_z_zz_xyyy[k] = -g_0_z_z_xyyy[k] + g_z_0_z_xyyy[k] - g_z_z_z_xyyy[k] * ab_z + g_z_z_z_xyyyz[k];

                g_z_z_zz_xyyz[k] = -g_0_z_z_xyyz[k] + g_z_0_z_xyyz[k] - g_z_z_z_xyyz[k] * ab_z + g_z_z_z_xyyzz[k];

                g_z_z_zz_xyzz[k] = -g_0_z_z_xyzz[k] + g_z_0_z_xyzz[k] - g_z_z_z_xyzz[k] * ab_z + g_z_z_z_xyzzz[k];

                g_z_z_zz_xzzz[k] = -g_0_z_z_xzzz[k] + g_z_0_z_xzzz[k] - g_z_z_z_xzzz[k] * ab_z + g_z_z_z_xzzzz[k];

                g_z_z_zz_yyyy[k] = -g_0_z_z_yyyy[k] + g_z_0_z_yyyy[k] - g_z_z_z_yyyy[k] * ab_z + g_z_z_z_yyyyz[k];

                g_z_z_zz_yyyz[k] = -g_0_z_z_yyyz[k] + g_z_0_z_yyyz[k] - g_z_z_z_yyyz[k] * ab_z + g_z_z_z_yyyzz[k];

                g_z_z_zz_yyzz[k] = -g_0_z_z_yyzz[k] + g_z_0_z_yyzz[k] - g_z_z_z_yyzz[k] * ab_z + g_z_z_z_yyzzz[k];

                g_z_z_zz_yzzz[k] = -g_0_z_z_yzzz[k] + g_z_0_z_yzzz[k] - g_z_z_z_yzzz[k] * ab_z + g_z_z_z_yzzzz[k];

                g_z_z_zz_zzzz[k] = -g_0_z_z_zzzz[k] + g_z_0_z_zzzz[k] - g_z_z_z_zzzz[k] * ab_z + g_z_z_z_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

