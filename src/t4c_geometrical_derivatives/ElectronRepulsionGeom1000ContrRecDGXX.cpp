#include "ElectronRepulsionGeom1000ContrRecDGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_dgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_dgxx,
                                            const size_t idx_pgxx,
                                            const size_t idx_geom_10_pgxx,
                                            const size_t idx_geom_10_phxx,
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

            const auto pg_off = idx_pgxx + i * dcomps + j;

            auto g_x_xxxx = cbuffer.data(pg_off + 0 * ccomps * dcomps);

            auto g_x_xxxy = cbuffer.data(pg_off + 1 * ccomps * dcomps);

            auto g_x_xxxz = cbuffer.data(pg_off + 2 * ccomps * dcomps);

            auto g_x_xxyy = cbuffer.data(pg_off + 3 * ccomps * dcomps);

            auto g_x_xxyz = cbuffer.data(pg_off + 4 * ccomps * dcomps);

            auto g_x_xxzz = cbuffer.data(pg_off + 5 * ccomps * dcomps);

            auto g_x_xyyy = cbuffer.data(pg_off + 6 * ccomps * dcomps);

            auto g_x_xyyz = cbuffer.data(pg_off + 7 * ccomps * dcomps);

            auto g_x_xyzz = cbuffer.data(pg_off + 8 * ccomps * dcomps);

            auto g_x_xzzz = cbuffer.data(pg_off + 9 * ccomps * dcomps);

            auto g_x_yyyy = cbuffer.data(pg_off + 10 * ccomps * dcomps);

            auto g_x_yyyz = cbuffer.data(pg_off + 11 * ccomps * dcomps);

            auto g_x_yyzz = cbuffer.data(pg_off + 12 * ccomps * dcomps);

            auto g_x_yzzz = cbuffer.data(pg_off + 13 * ccomps * dcomps);

            auto g_x_zzzz = cbuffer.data(pg_off + 14 * ccomps * dcomps);

            auto g_y_xxxx = cbuffer.data(pg_off + 15 * ccomps * dcomps);

            auto g_y_xxxy = cbuffer.data(pg_off + 16 * ccomps * dcomps);

            auto g_y_xxxz = cbuffer.data(pg_off + 17 * ccomps * dcomps);

            auto g_y_xxyy = cbuffer.data(pg_off + 18 * ccomps * dcomps);

            auto g_y_xxyz = cbuffer.data(pg_off + 19 * ccomps * dcomps);

            auto g_y_xxzz = cbuffer.data(pg_off + 20 * ccomps * dcomps);

            auto g_y_xyyy = cbuffer.data(pg_off + 21 * ccomps * dcomps);

            auto g_y_xyyz = cbuffer.data(pg_off + 22 * ccomps * dcomps);

            auto g_y_xyzz = cbuffer.data(pg_off + 23 * ccomps * dcomps);

            auto g_y_xzzz = cbuffer.data(pg_off + 24 * ccomps * dcomps);

            auto g_y_yyyy = cbuffer.data(pg_off + 25 * ccomps * dcomps);

            auto g_y_yyyz = cbuffer.data(pg_off + 26 * ccomps * dcomps);

            auto g_y_yyzz = cbuffer.data(pg_off + 27 * ccomps * dcomps);

            auto g_y_yzzz = cbuffer.data(pg_off + 28 * ccomps * dcomps);

            auto g_y_zzzz = cbuffer.data(pg_off + 29 * ccomps * dcomps);

            auto g_z_xxxx = cbuffer.data(pg_off + 30 * ccomps * dcomps);

            auto g_z_xxxy = cbuffer.data(pg_off + 31 * ccomps * dcomps);

            auto g_z_xxxz = cbuffer.data(pg_off + 32 * ccomps * dcomps);

            auto g_z_xxyy = cbuffer.data(pg_off + 33 * ccomps * dcomps);

            auto g_z_xxyz = cbuffer.data(pg_off + 34 * ccomps * dcomps);

            auto g_z_xxzz = cbuffer.data(pg_off + 35 * ccomps * dcomps);

            auto g_z_xyyy = cbuffer.data(pg_off + 36 * ccomps * dcomps);

            auto g_z_xyyz = cbuffer.data(pg_off + 37 * ccomps * dcomps);

            auto g_z_xyzz = cbuffer.data(pg_off + 38 * ccomps * dcomps);

            auto g_z_xzzz = cbuffer.data(pg_off + 39 * ccomps * dcomps);

            auto g_z_yyyy = cbuffer.data(pg_off + 40 * ccomps * dcomps);

            auto g_z_yyyz = cbuffer.data(pg_off + 41 * ccomps * dcomps);

            auto g_z_yyzz = cbuffer.data(pg_off + 42 * ccomps * dcomps);

            auto g_z_yzzz = cbuffer.data(pg_off + 43 * ccomps * dcomps);

            auto g_z_zzzz = cbuffer.data(pg_off + 44 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_10_off = idx_geom_10_phxx + i * dcomps + j;

            auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 188 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dgxx

            const auto dg_geom_10_off = idx_geom_10_dgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxx, g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxy, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxxzz, g_x_0_x_xxyy, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xxzzz, g_x_0_x_xyyy, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_xzzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz, g_x_0_xx_xxxx, g_x_0_xx_xxxy, g_x_0_xx_xxxz, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxzz, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyzz, g_x_0_xx_xzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyzz, g_x_0_xx_yzzz, g_x_0_xx_zzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz, g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy, g_x_xyyz, g_x_xyzz, g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_xxxx[k] = -g_x_xxxx[k] - g_x_0_x_xxxx[k] * ab_x + g_x_0_x_xxxxx[k];

                g_x_0_xx_xxxy[k] = -g_x_xxxy[k] - g_x_0_x_xxxy[k] * ab_x + g_x_0_x_xxxxy[k];

                g_x_0_xx_xxxz[k] = -g_x_xxxz[k] - g_x_0_x_xxxz[k] * ab_x + g_x_0_x_xxxxz[k];

                g_x_0_xx_xxyy[k] = -g_x_xxyy[k] - g_x_0_x_xxyy[k] * ab_x + g_x_0_x_xxxyy[k];

                g_x_0_xx_xxyz[k] = -g_x_xxyz[k] - g_x_0_x_xxyz[k] * ab_x + g_x_0_x_xxxyz[k];

                g_x_0_xx_xxzz[k] = -g_x_xxzz[k] - g_x_0_x_xxzz[k] * ab_x + g_x_0_x_xxxzz[k];

                g_x_0_xx_xyyy[k] = -g_x_xyyy[k] - g_x_0_x_xyyy[k] * ab_x + g_x_0_x_xxyyy[k];

                g_x_0_xx_xyyz[k] = -g_x_xyyz[k] - g_x_0_x_xyyz[k] * ab_x + g_x_0_x_xxyyz[k];

                g_x_0_xx_xyzz[k] = -g_x_xyzz[k] - g_x_0_x_xyzz[k] * ab_x + g_x_0_x_xxyzz[k];

                g_x_0_xx_xzzz[k] = -g_x_xzzz[k] - g_x_0_x_xzzz[k] * ab_x + g_x_0_x_xxzzz[k];

                g_x_0_xx_yyyy[k] = -g_x_yyyy[k] - g_x_0_x_yyyy[k] * ab_x + g_x_0_x_xyyyy[k];

                g_x_0_xx_yyyz[k] = -g_x_yyyz[k] - g_x_0_x_yyyz[k] * ab_x + g_x_0_x_xyyyz[k];

                g_x_0_xx_yyzz[k] = -g_x_yyzz[k] - g_x_0_x_yyzz[k] * ab_x + g_x_0_x_xyyzz[k];

                g_x_0_xx_yzzz[k] = -g_x_yzzz[k] - g_x_0_x_yzzz[k] * ab_x + g_x_0_x_xyzzz[k];

                g_x_0_xx_zzzz[k] = -g_x_zzzz[k] - g_x_0_x_zzzz[k] * ab_x + g_x_0_x_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxx, g_x_0_x_xxxxy, g_x_0_x_xxxy, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyz, g_x_0_x_yyyzz, g_x_0_x_yyzz, g_x_0_x_yyzzz, g_x_0_x_yzzz, g_x_0_x_yzzzz, g_x_0_x_zzzz, g_x_0_xy_xxxx, g_x_0_xy_xxxy, g_x_0_xy_xxxz, g_x_0_xy_xxyy, g_x_0_xy_xxyz, g_x_0_xy_xxzz, g_x_0_xy_xyyy, g_x_0_xy_xyyz, g_x_0_xy_xyzz, g_x_0_xy_xzzz, g_x_0_xy_yyyy, g_x_0_xy_yyyz, g_x_0_xy_yyzz, g_x_0_xy_yzzz, g_x_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_xxxx[k] = -g_x_0_x_xxxx[k] * ab_y + g_x_0_x_xxxxy[k];

                g_x_0_xy_xxxy[k] = -g_x_0_x_xxxy[k] * ab_y + g_x_0_x_xxxyy[k];

                g_x_0_xy_xxxz[k] = -g_x_0_x_xxxz[k] * ab_y + g_x_0_x_xxxyz[k];

                g_x_0_xy_xxyy[k] = -g_x_0_x_xxyy[k] * ab_y + g_x_0_x_xxyyy[k];

                g_x_0_xy_xxyz[k] = -g_x_0_x_xxyz[k] * ab_y + g_x_0_x_xxyyz[k];

                g_x_0_xy_xxzz[k] = -g_x_0_x_xxzz[k] * ab_y + g_x_0_x_xxyzz[k];

                g_x_0_xy_xyyy[k] = -g_x_0_x_xyyy[k] * ab_y + g_x_0_x_xyyyy[k];

                g_x_0_xy_xyyz[k] = -g_x_0_x_xyyz[k] * ab_y + g_x_0_x_xyyyz[k];

                g_x_0_xy_xyzz[k] = -g_x_0_x_xyzz[k] * ab_y + g_x_0_x_xyyzz[k];

                g_x_0_xy_xzzz[k] = -g_x_0_x_xzzz[k] * ab_y + g_x_0_x_xyzzz[k];

                g_x_0_xy_yyyy[k] = -g_x_0_x_yyyy[k] * ab_y + g_x_0_x_yyyyy[k];

                g_x_0_xy_yyyz[k] = -g_x_0_x_yyyz[k] * ab_y + g_x_0_x_yyyyz[k];

                g_x_0_xy_yyzz[k] = -g_x_0_x_yyzz[k] * ab_y + g_x_0_x_yyyzz[k];

                g_x_0_xy_yzzz[k] = -g_x_0_x_yzzz[k] * ab_y + g_x_0_x_yyzzz[k];

                g_x_0_xy_zzzz[k] = -g_x_0_x_zzzz[k] * ab_y + g_x_0_x_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxx, g_x_0_x_xxxxz, g_x_0_x_xxxy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxxzz, g_x_0_x_xxyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xxzzz, g_x_0_x_xyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_xzzzz, g_x_0_x_yyyy, g_x_0_x_yyyyz, g_x_0_x_yyyz, g_x_0_x_yyyzz, g_x_0_x_yyzz, g_x_0_x_yyzzz, g_x_0_x_yzzz, g_x_0_x_yzzzz, g_x_0_x_zzzz, g_x_0_x_zzzzz, g_x_0_xz_xxxx, g_x_0_xz_xxxy, g_x_0_xz_xxxz, g_x_0_xz_xxyy, g_x_0_xz_xxyz, g_x_0_xz_xxzz, g_x_0_xz_xyyy, g_x_0_xz_xyyz, g_x_0_xz_xyzz, g_x_0_xz_xzzz, g_x_0_xz_yyyy, g_x_0_xz_yyyz, g_x_0_xz_yyzz, g_x_0_xz_yzzz, g_x_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_xxxx[k] = -g_x_0_x_xxxx[k] * ab_z + g_x_0_x_xxxxz[k];

                g_x_0_xz_xxxy[k] = -g_x_0_x_xxxy[k] * ab_z + g_x_0_x_xxxyz[k];

                g_x_0_xz_xxxz[k] = -g_x_0_x_xxxz[k] * ab_z + g_x_0_x_xxxzz[k];

                g_x_0_xz_xxyy[k] = -g_x_0_x_xxyy[k] * ab_z + g_x_0_x_xxyyz[k];

                g_x_0_xz_xxyz[k] = -g_x_0_x_xxyz[k] * ab_z + g_x_0_x_xxyzz[k];

                g_x_0_xz_xxzz[k] = -g_x_0_x_xxzz[k] * ab_z + g_x_0_x_xxzzz[k];

                g_x_0_xz_xyyy[k] = -g_x_0_x_xyyy[k] * ab_z + g_x_0_x_xyyyz[k];

                g_x_0_xz_xyyz[k] = -g_x_0_x_xyyz[k] * ab_z + g_x_0_x_xyyzz[k];

                g_x_0_xz_xyzz[k] = -g_x_0_x_xyzz[k] * ab_z + g_x_0_x_xyzzz[k];

                g_x_0_xz_xzzz[k] = -g_x_0_x_xzzz[k] * ab_z + g_x_0_x_xzzzz[k];

                g_x_0_xz_yyyy[k] = -g_x_0_x_yyyy[k] * ab_z + g_x_0_x_yyyyz[k];

                g_x_0_xz_yyyz[k] = -g_x_0_x_yyyz[k] * ab_z + g_x_0_x_yyyzz[k];

                g_x_0_xz_yyzz[k] = -g_x_0_x_yyzz[k] * ab_z + g_x_0_x_yyzzz[k];

                g_x_0_xz_yzzz[k] = -g_x_0_x_yzzz[k] * ab_z + g_x_0_x_yzzzz[k];

                g_x_0_xz_zzzz[k] = -g_x_0_x_zzzz[k] * ab_z + g_x_0_x_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxx, g_x_0_y_xxxxy, g_x_0_y_xxxy, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyz, g_x_0_y_xxyzz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyz, g_x_0_y_xyyzz, g_x_0_y_xyzz, g_x_0_y_xyzzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyz, g_x_0_y_yyyzz, g_x_0_y_yyzz, g_x_0_y_yyzzz, g_x_0_y_yzzz, g_x_0_y_yzzzz, g_x_0_y_zzzz, g_x_0_yy_xxxx, g_x_0_yy_xxxy, g_x_0_yy_xxxz, g_x_0_yy_xxyy, g_x_0_yy_xxyz, g_x_0_yy_xxzz, g_x_0_yy_xyyy, g_x_0_yy_xyyz, g_x_0_yy_xyzz, g_x_0_yy_xzzz, g_x_0_yy_yyyy, g_x_0_yy_yyyz, g_x_0_yy_yyzz, g_x_0_yy_yzzz, g_x_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_xxxx[k] = -g_x_0_y_xxxx[k] * ab_y + g_x_0_y_xxxxy[k];

                g_x_0_yy_xxxy[k] = -g_x_0_y_xxxy[k] * ab_y + g_x_0_y_xxxyy[k];

                g_x_0_yy_xxxz[k] = -g_x_0_y_xxxz[k] * ab_y + g_x_0_y_xxxyz[k];

                g_x_0_yy_xxyy[k] = -g_x_0_y_xxyy[k] * ab_y + g_x_0_y_xxyyy[k];

                g_x_0_yy_xxyz[k] = -g_x_0_y_xxyz[k] * ab_y + g_x_0_y_xxyyz[k];

                g_x_0_yy_xxzz[k] = -g_x_0_y_xxzz[k] * ab_y + g_x_0_y_xxyzz[k];

                g_x_0_yy_xyyy[k] = -g_x_0_y_xyyy[k] * ab_y + g_x_0_y_xyyyy[k];

                g_x_0_yy_xyyz[k] = -g_x_0_y_xyyz[k] * ab_y + g_x_0_y_xyyyz[k];

                g_x_0_yy_xyzz[k] = -g_x_0_y_xyzz[k] * ab_y + g_x_0_y_xyyzz[k];

                g_x_0_yy_xzzz[k] = -g_x_0_y_xzzz[k] * ab_y + g_x_0_y_xyzzz[k];

                g_x_0_yy_yyyy[k] = -g_x_0_y_yyyy[k] * ab_y + g_x_0_y_yyyyy[k];

                g_x_0_yy_yyyz[k] = -g_x_0_y_yyyz[k] * ab_y + g_x_0_y_yyyyz[k];

                g_x_0_yy_yyzz[k] = -g_x_0_y_yyzz[k] * ab_y + g_x_0_y_yyyzz[k];

                g_x_0_yy_yzzz[k] = -g_x_0_y_yzzz[k] * ab_y + g_x_0_y_yyzzz[k];

                g_x_0_yy_zzzz[k] = -g_x_0_y_zzzz[k] * ab_y + g_x_0_y_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yz_xxxx, g_x_0_yz_xxxy, g_x_0_yz_xxxz, g_x_0_yz_xxyy, g_x_0_yz_xxyz, g_x_0_yz_xxzz, g_x_0_yz_xyyy, g_x_0_yz_xyyz, g_x_0_yz_xyzz, g_x_0_yz_xzzz, g_x_0_yz_yyyy, g_x_0_yz_yyyz, g_x_0_yz_yyzz, g_x_0_yz_yzzz, g_x_0_yz_zzzz, g_x_0_z_xxxx, g_x_0_z_xxxxy, g_x_0_z_xxxy, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyz, g_x_0_z_xxyzz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyz, g_x_0_z_xyyzz, g_x_0_z_xyzz, g_x_0_z_xyzzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyz, g_x_0_z_yyyzz, g_x_0_z_yyzz, g_x_0_z_yyzzz, g_x_0_z_yzzz, g_x_0_z_yzzzz, g_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_xxxx[k] = -g_x_0_z_xxxx[k] * ab_y + g_x_0_z_xxxxy[k];

                g_x_0_yz_xxxy[k] = -g_x_0_z_xxxy[k] * ab_y + g_x_0_z_xxxyy[k];

                g_x_0_yz_xxxz[k] = -g_x_0_z_xxxz[k] * ab_y + g_x_0_z_xxxyz[k];

                g_x_0_yz_xxyy[k] = -g_x_0_z_xxyy[k] * ab_y + g_x_0_z_xxyyy[k];

                g_x_0_yz_xxyz[k] = -g_x_0_z_xxyz[k] * ab_y + g_x_0_z_xxyyz[k];

                g_x_0_yz_xxzz[k] = -g_x_0_z_xxzz[k] * ab_y + g_x_0_z_xxyzz[k];

                g_x_0_yz_xyyy[k] = -g_x_0_z_xyyy[k] * ab_y + g_x_0_z_xyyyy[k];

                g_x_0_yz_xyyz[k] = -g_x_0_z_xyyz[k] * ab_y + g_x_0_z_xyyyz[k];

                g_x_0_yz_xyzz[k] = -g_x_0_z_xyzz[k] * ab_y + g_x_0_z_xyyzz[k];

                g_x_0_yz_xzzz[k] = -g_x_0_z_xzzz[k] * ab_y + g_x_0_z_xyzzz[k];

                g_x_0_yz_yyyy[k] = -g_x_0_z_yyyy[k] * ab_y + g_x_0_z_yyyyy[k];

                g_x_0_yz_yyyz[k] = -g_x_0_z_yyyz[k] * ab_y + g_x_0_z_yyyyz[k];

                g_x_0_yz_yyzz[k] = -g_x_0_z_yyzz[k] * ab_y + g_x_0_z_yyyzz[k];

                g_x_0_yz_yzzz[k] = -g_x_0_z_yzzz[k] * ab_y + g_x_0_z_yyzzz[k];

                g_x_0_yz_zzzz[k] = -g_x_0_z_zzzz[k] * ab_y + g_x_0_z_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxx, g_x_0_z_xxxxz, g_x_0_z_xxxy, g_x_0_z_xxxyz, g_x_0_z_xxxz, g_x_0_z_xxxzz, g_x_0_z_xxyy, g_x_0_z_xxyyz, g_x_0_z_xxyz, g_x_0_z_xxyzz, g_x_0_z_xxzz, g_x_0_z_xxzzz, g_x_0_z_xyyy, g_x_0_z_xyyyz, g_x_0_z_xyyz, g_x_0_z_xyyzz, g_x_0_z_xyzz, g_x_0_z_xyzzz, g_x_0_z_xzzz, g_x_0_z_xzzzz, g_x_0_z_yyyy, g_x_0_z_yyyyz, g_x_0_z_yyyz, g_x_0_z_yyyzz, g_x_0_z_yyzz, g_x_0_z_yyzzz, g_x_0_z_yzzz, g_x_0_z_yzzzz, g_x_0_z_zzzz, g_x_0_z_zzzzz, g_x_0_zz_xxxx, g_x_0_zz_xxxy, g_x_0_zz_xxxz, g_x_0_zz_xxyy, g_x_0_zz_xxyz, g_x_0_zz_xxzz, g_x_0_zz_xyyy, g_x_0_zz_xyyz, g_x_0_zz_xyzz, g_x_0_zz_xzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyz, g_x_0_zz_yyzz, g_x_0_zz_yzzz, g_x_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_xxxx[k] = -g_x_0_z_xxxx[k] * ab_z + g_x_0_z_xxxxz[k];

                g_x_0_zz_xxxy[k] = -g_x_0_z_xxxy[k] * ab_z + g_x_0_z_xxxyz[k];

                g_x_0_zz_xxxz[k] = -g_x_0_z_xxxz[k] * ab_z + g_x_0_z_xxxzz[k];

                g_x_0_zz_xxyy[k] = -g_x_0_z_xxyy[k] * ab_z + g_x_0_z_xxyyz[k];

                g_x_0_zz_xxyz[k] = -g_x_0_z_xxyz[k] * ab_z + g_x_0_z_xxyzz[k];

                g_x_0_zz_xxzz[k] = -g_x_0_z_xxzz[k] * ab_z + g_x_0_z_xxzzz[k];

                g_x_0_zz_xyyy[k] = -g_x_0_z_xyyy[k] * ab_z + g_x_0_z_xyyyz[k];

                g_x_0_zz_xyyz[k] = -g_x_0_z_xyyz[k] * ab_z + g_x_0_z_xyyzz[k];

                g_x_0_zz_xyzz[k] = -g_x_0_z_xyzz[k] * ab_z + g_x_0_z_xyzzz[k];

                g_x_0_zz_xzzz[k] = -g_x_0_z_xzzz[k] * ab_z + g_x_0_z_xzzzz[k];

                g_x_0_zz_yyyy[k] = -g_x_0_z_yyyy[k] * ab_z + g_x_0_z_yyyyz[k];

                g_x_0_zz_yyyz[k] = -g_x_0_z_yyyz[k] * ab_z + g_x_0_z_yyyzz[k];

                g_x_0_zz_yyzz[k] = -g_x_0_z_yyzz[k] * ab_z + g_x_0_z_yyzzz[k];

                g_x_0_zz_yzzz[k] = -g_x_0_z_yzzz[k] * ab_z + g_x_0_z_yzzzz[k];

                g_x_0_zz_zzzz[k] = -g_x_0_z_zzzz[k] * ab_z + g_x_0_z_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxxx, g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxy, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxz, g_y_0_x_xxxzz, g_y_0_x_xxyy, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyz, g_y_0_x_xxyzz, g_y_0_x_xxzz, g_y_0_x_xxzzz, g_y_0_x_xyyy, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyz, g_y_0_x_xyyzz, g_y_0_x_xyzz, g_y_0_x_xyzzz, g_y_0_x_xzzz, g_y_0_x_xzzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz, g_y_0_xx_xxxx, g_y_0_xx_xxxy, g_y_0_xx_xxxz, g_y_0_xx_xxyy, g_y_0_xx_xxyz, g_y_0_xx_xxzz, g_y_0_xx_xyyy, g_y_0_xx_xyyz, g_y_0_xx_xyzz, g_y_0_xx_xzzz, g_y_0_xx_yyyy, g_y_0_xx_yyyz, g_y_0_xx_yyzz, g_y_0_xx_yzzz, g_y_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_xxxx[k] = -g_y_0_x_xxxx[k] * ab_x + g_y_0_x_xxxxx[k];

                g_y_0_xx_xxxy[k] = -g_y_0_x_xxxy[k] * ab_x + g_y_0_x_xxxxy[k];

                g_y_0_xx_xxxz[k] = -g_y_0_x_xxxz[k] * ab_x + g_y_0_x_xxxxz[k];

                g_y_0_xx_xxyy[k] = -g_y_0_x_xxyy[k] * ab_x + g_y_0_x_xxxyy[k];

                g_y_0_xx_xxyz[k] = -g_y_0_x_xxyz[k] * ab_x + g_y_0_x_xxxyz[k];

                g_y_0_xx_xxzz[k] = -g_y_0_x_xxzz[k] * ab_x + g_y_0_x_xxxzz[k];

                g_y_0_xx_xyyy[k] = -g_y_0_x_xyyy[k] * ab_x + g_y_0_x_xxyyy[k];

                g_y_0_xx_xyyz[k] = -g_y_0_x_xyyz[k] * ab_x + g_y_0_x_xxyyz[k];

                g_y_0_xx_xyzz[k] = -g_y_0_x_xyzz[k] * ab_x + g_y_0_x_xxyzz[k];

                g_y_0_xx_xzzz[k] = -g_y_0_x_xzzz[k] * ab_x + g_y_0_x_xxzzz[k];

                g_y_0_xx_yyyy[k] = -g_y_0_x_yyyy[k] * ab_x + g_y_0_x_xyyyy[k];

                g_y_0_xx_yyyz[k] = -g_y_0_x_yyyz[k] * ab_x + g_y_0_x_xyyyz[k];

                g_y_0_xx_yyzz[k] = -g_y_0_x_yyzz[k] * ab_x + g_y_0_x_xyyzz[k];

                g_y_0_xx_yzzz[k] = -g_y_0_x_yzzz[k] * ab_x + g_y_0_x_xyzzz[k];

                g_y_0_xx_zzzz[k] = -g_y_0_x_zzzz[k] * ab_x + g_y_0_x_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xxxx, g_y_0_xy_xxxy, g_y_0_xy_xxxz, g_y_0_xy_xxyy, g_y_0_xy_xxyz, g_y_0_xy_xxzz, g_y_0_xy_xyyy, g_y_0_xy_xyyz, g_y_0_xy_xyzz, g_y_0_xy_xzzz, g_y_0_xy_yyyy, g_y_0_xy_yyyz, g_y_0_xy_yyzz, g_y_0_xy_yzzz, g_y_0_xy_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxy, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxxzz, g_y_0_y_xxyy, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xxzzz, g_y_0_y_xyyy, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_xzzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_xxxx[k] = -g_y_0_y_xxxx[k] * ab_x + g_y_0_y_xxxxx[k];

                g_y_0_xy_xxxy[k] = -g_y_0_y_xxxy[k] * ab_x + g_y_0_y_xxxxy[k];

                g_y_0_xy_xxxz[k] = -g_y_0_y_xxxz[k] * ab_x + g_y_0_y_xxxxz[k];

                g_y_0_xy_xxyy[k] = -g_y_0_y_xxyy[k] * ab_x + g_y_0_y_xxxyy[k];

                g_y_0_xy_xxyz[k] = -g_y_0_y_xxyz[k] * ab_x + g_y_0_y_xxxyz[k];

                g_y_0_xy_xxzz[k] = -g_y_0_y_xxzz[k] * ab_x + g_y_0_y_xxxzz[k];

                g_y_0_xy_xyyy[k] = -g_y_0_y_xyyy[k] * ab_x + g_y_0_y_xxyyy[k];

                g_y_0_xy_xyyz[k] = -g_y_0_y_xyyz[k] * ab_x + g_y_0_y_xxyyz[k];

                g_y_0_xy_xyzz[k] = -g_y_0_y_xyzz[k] * ab_x + g_y_0_y_xxyzz[k];

                g_y_0_xy_xzzz[k] = -g_y_0_y_xzzz[k] * ab_x + g_y_0_y_xxzzz[k];

                g_y_0_xy_yyyy[k] = -g_y_0_y_yyyy[k] * ab_x + g_y_0_y_xyyyy[k];

                g_y_0_xy_yyyz[k] = -g_y_0_y_yyyz[k] * ab_x + g_y_0_y_xyyyz[k];

                g_y_0_xy_yyzz[k] = -g_y_0_y_yyzz[k] * ab_x + g_y_0_y_xyyzz[k];

                g_y_0_xy_yzzz[k] = -g_y_0_y_yzzz[k] * ab_x + g_y_0_y_xyzzz[k];

                g_y_0_xy_zzzz[k] = -g_y_0_y_zzzz[k] * ab_x + g_y_0_y_xzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xz_xxxx, g_y_0_xz_xxxy, g_y_0_xz_xxxz, g_y_0_xz_xxyy, g_y_0_xz_xxyz, g_y_0_xz_xxzz, g_y_0_xz_xyyy, g_y_0_xz_xyyz, g_y_0_xz_xyzz, g_y_0_xz_xzzz, g_y_0_xz_yyyy, g_y_0_xz_yyyz, g_y_0_xz_yyzz, g_y_0_xz_yzzz, g_y_0_xz_zzzz, g_y_0_z_xxxx, g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxy, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxz, g_y_0_z_xxxzz, g_y_0_z_xxyy, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyz, g_y_0_z_xxyzz, g_y_0_z_xxzz, g_y_0_z_xxzzz, g_y_0_z_xyyy, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyz, g_y_0_z_xyyzz, g_y_0_z_xyzz, g_y_0_z_xyzzz, g_y_0_z_xzzz, g_y_0_z_xzzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_xxxx[k] = -g_y_0_z_xxxx[k] * ab_x + g_y_0_z_xxxxx[k];

                g_y_0_xz_xxxy[k] = -g_y_0_z_xxxy[k] * ab_x + g_y_0_z_xxxxy[k];

                g_y_0_xz_xxxz[k] = -g_y_0_z_xxxz[k] * ab_x + g_y_0_z_xxxxz[k];

                g_y_0_xz_xxyy[k] = -g_y_0_z_xxyy[k] * ab_x + g_y_0_z_xxxyy[k];

                g_y_0_xz_xxyz[k] = -g_y_0_z_xxyz[k] * ab_x + g_y_0_z_xxxyz[k];

                g_y_0_xz_xxzz[k] = -g_y_0_z_xxzz[k] * ab_x + g_y_0_z_xxxzz[k];

                g_y_0_xz_xyyy[k] = -g_y_0_z_xyyy[k] * ab_x + g_y_0_z_xxyyy[k];

                g_y_0_xz_xyyz[k] = -g_y_0_z_xyyz[k] * ab_x + g_y_0_z_xxyyz[k];

                g_y_0_xz_xyzz[k] = -g_y_0_z_xyzz[k] * ab_x + g_y_0_z_xxyzz[k];

                g_y_0_xz_xzzz[k] = -g_y_0_z_xzzz[k] * ab_x + g_y_0_z_xxzzz[k];

                g_y_0_xz_yyyy[k] = -g_y_0_z_yyyy[k] * ab_x + g_y_0_z_xyyyy[k];

                g_y_0_xz_yyyz[k] = -g_y_0_z_yyyz[k] * ab_x + g_y_0_z_xyyyz[k];

                g_y_0_xz_yyzz[k] = -g_y_0_z_yyzz[k] * ab_x + g_y_0_z_xyyzz[k];

                g_y_0_xz_yzzz[k] = -g_y_0_z_yzzz[k] * ab_x + g_y_0_z_xyzzz[k];

                g_y_0_xz_zzzz[k] = -g_y_0_z_zzzz[k] * ab_x + g_y_0_z_xzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxx, g_y_0_y_xxxxy, g_y_0_y_xxxy, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyz, g_y_0_y_yyyzz, g_y_0_y_yyzz, g_y_0_y_yyzzz, g_y_0_y_yzzz, g_y_0_y_yzzzz, g_y_0_y_zzzz, g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz, g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxyy, g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz, g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz, g_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_xxxx[k] = -g_y_xxxx[k] - g_y_0_y_xxxx[k] * ab_y + g_y_0_y_xxxxy[k];

                g_y_0_yy_xxxy[k] = -g_y_xxxy[k] - g_y_0_y_xxxy[k] * ab_y + g_y_0_y_xxxyy[k];

                g_y_0_yy_xxxz[k] = -g_y_xxxz[k] - g_y_0_y_xxxz[k] * ab_y + g_y_0_y_xxxyz[k];

                g_y_0_yy_xxyy[k] = -g_y_xxyy[k] - g_y_0_y_xxyy[k] * ab_y + g_y_0_y_xxyyy[k];

                g_y_0_yy_xxyz[k] = -g_y_xxyz[k] - g_y_0_y_xxyz[k] * ab_y + g_y_0_y_xxyyz[k];

                g_y_0_yy_xxzz[k] = -g_y_xxzz[k] - g_y_0_y_xxzz[k] * ab_y + g_y_0_y_xxyzz[k];

                g_y_0_yy_xyyy[k] = -g_y_xyyy[k] - g_y_0_y_xyyy[k] * ab_y + g_y_0_y_xyyyy[k];

                g_y_0_yy_xyyz[k] = -g_y_xyyz[k] - g_y_0_y_xyyz[k] * ab_y + g_y_0_y_xyyyz[k];

                g_y_0_yy_xyzz[k] = -g_y_xyzz[k] - g_y_0_y_xyzz[k] * ab_y + g_y_0_y_xyyzz[k];

                g_y_0_yy_xzzz[k] = -g_y_xzzz[k] - g_y_0_y_xzzz[k] * ab_y + g_y_0_y_xyzzz[k];

                g_y_0_yy_yyyy[k] = -g_y_yyyy[k] - g_y_0_y_yyyy[k] * ab_y + g_y_0_y_yyyyy[k];

                g_y_0_yy_yyyz[k] = -g_y_yyyz[k] - g_y_0_y_yyyz[k] * ab_y + g_y_0_y_yyyyz[k];

                g_y_0_yy_yyzz[k] = -g_y_yyzz[k] - g_y_0_y_yyzz[k] * ab_y + g_y_0_y_yyyzz[k];

                g_y_0_yy_yzzz[k] = -g_y_yzzz[k] - g_y_0_y_yzzz[k] * ab_y + g_y_0_y_yyzzz[k];

                g_y_0_yy_zzzz[k] = -g_y_zzzz[k] - g_y_0_y_zzzz[k] * ab_y + g_y_0_y_yzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxx, g_y_0_y_xxxxz, g_y_0_y_xxxy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxxzz, g_y_0_y_xxyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xxzzz, g_y_0_y_xyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_xzzzz, g_y_0_y_yyyy, g_y_0_y_yyyyz, g_y_0_y_yyyz, g_y_0_y_yyyzz, g_y_0_y_yyzz, g_y_0_y_yyzzz, g_y_0_y_yzzz, g_y_0_y_yzzzz, g_y_0_y_zzzz, g_y_0_y_zzzzz, g_y_0_yz_xxxx, g_y_0_yz_xxxy, g_y_0_yz_xxxz, g_y_0_yz_xxyy, g_y_0_yz_xxyz, g_y_0_yz_xxzz, g_y_0_yz_xyyy, g_y_0_yz_xyyz, g_y_0_yz_xyzz, g_y_0_yz_xzzz, g_y_0_yz_yyyy, g_y_0_yz_yyyz, g_y_0_yz_yyzz, g_y_0_yz_yzzz, g_y_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_xxxx[k] = -g_y_0_y_xxxx[k] * ab_z + g_y_0_y_xxxxz[k];

                g_y_0_yz_xxxy[k] = -g_y_0_y_xxxy[k] * ab_z + g_y_0_y_xxxyz[k];

                g_y_0_yz_xxxz[k] = -g_y_0_y_xxxz[k] * ab_z + g_y_0_y_xxxzz[k];

                g_y_0_yz_xxyy[k] = -g_y_0_y_xxyy[k] * ab_z + g_y_0_y_xxyyz[k];

                g_y_0_yz_xxyz[k] = -g_y_0_y_xxyz[k] * ab_z + g_y_0_y_xxyzz[k];

                g_y_0_yz_xxzz[k] = -g_y_0_y_xxzz[k] * ab_z + g_y_0_y_xxzzz[k];

                g_y_0_yz_xyyy[k] = -g_y_0_y_xyyy[k] * ab_z + g_y_0_y_xyyyz[k];

                g_y_0_yz_xyyz[k] = -g_y_0_y_xyyz[k] * ab_z + g_y_0_y_xyyzz[k];

                g_y_0_yz_xyzz[k] = -g_y_0_y_xyzz[k] * ab_z + g_y_0_y_xyzzz[k];

                g_y_0_yz_xzzz[k] = -g_y_0_y_xzzz[k] * ab_z + g_y_0_y_xzzzz[k];

                g_y_0_yz_yyyy[k] = -g_y_0_y_yyyy[k] * ab_z + g_y_0_y_yyyyz[k];

                g_y_0_yz_yyyz[k] = -g_y_0_y_yyyz[k] * ab_z + g_y_0_y_yyyzz[k];

                g_y_0_yz_yyzz[k] = -g_y_0_y_yyzz[k] * ab_z + g_y_0_y_yyzzz[k];

                g_y_0_yz_yzzz[k] = -g_y_0_y_yzzz[k] * ab_z + g_y_0_y_yzzzz[k];

                g_y_0_yz_zzzz[k] = -g_y_0_y_zzzz[k] * ab_z + g_y_0_y_zzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxx, g_y_0_z_xxxxz, g_y_0_z_xxxy, g_y_0_z_xxxyz, g_y_0_z_xxxz, g_y_0_z_xxxzz, g_y_0_z_xxyy, g_y_0_z_xxyyz, g_y_0_z_xxyz, g_y_0_z_xxyzz, g_y_0_z_xxzz, g_y_0_z_xxzzz, g_y_0_z_xyyy, g_y_0_z_xyyyz, g_y_0_z_xyyz, g_y_0_z_xyyzz, g_y_0_z_xyzz, g_y_0_z_xyzzz, g_y_0_z_xzzz, g_y_0_z_xzzzz, g_y_0_z_yyyy, g_y_0_z_yyyyz, g_y_0_z_yyyz, g_y_0_z_yyyzz, g_y_0_z_yyzz, g_y_0_z_yyzzz, g_y_0_z_yzzz, g_y_0_z_yzzzz, g_y_0_z_zzzz, g_y_0_z_zzzzz, g_y_0_zz_xxxx, g_y_0_zz_xxxy, g_y_0_zz_xxxz, g_y_0_zz_xxyy, g_y_0_zz_xxyz, g_y_0_zz_xxzz, g_y_0_zz_xyyy, g_y_0_zz_xyyz, g_y_0_zz_xyzz, g_y_0_zz_xzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyz, g_y_0_zz_yyzz, g_y_0_zz_yzzz, g_y_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_xxxx[k] = -g_y_0_z_xxxx[k] * ab_z + g_y_0_z_xxxxz[k];

                g_y_0_zz_xxxy[k] = -g_y_0_z_xxxy[k] * ab_z + g_y_0_z_xxxyz[k];

                g_y_0_zz_xxxz[k] = -g_y_0_z_xxxz[k] * ab_z + g_y_0_z_xxxzz[k];

                g_y_0_zz_xxyy[k] = -g_y_0_z_xxyy[k] * ab_z + g_y_0_z_xxyyz[k];

                g_y_0_zz_xxyz[k] = -g_y_0_z_xxyz[k] * ab_z + g_y_0_z_xxyzz[k];

                g_y_0_zz_xxzz[k] = -g_y_0_z_xxzz[k] * ab_z + g_y_0_z_xxzzz[k];

                g_y_0_zz_xyyy[k] = -g_y_0_z_xyyy[k] * ab_z + g_y_0_z_xyyyz[k];

                g_y_0_zz_xyyz[k] = -g_y_0_z_xyyz[k] * ab_z + g_y_0_z_xyyzz[k];

                g_y_0_zz_xyzz[k] = -g_y_0_z_xyzz[k] * ab_z + g_y_0_z_xyzzz[k];

                g_y_0_zz_xzzz[k] = -g_y_0_z_xzzz[k] * ab_z + g_y_0_z_xzzzz[k];

                g_y_0_zz_yyyy[k] = -g_y_0_z_yyyy[k] * ab_z + g_y_0_z_yyyyz[k];

                g_y_0_zz_yyyz[k] = -g_y_0_z_yyyz[k] * ab_z + g_y_0_z_yyyzz[k];

                g_y_0_zz_yyzz[k] = -g_y_0_z_yyzz[k] * ab_z + g_y_0_z_yyzzz[k];

                g_y_0_zz_yzzz[k] = -g_y_0_z_yzzz[k] * ab_z + g_y_0_z_yzzzz[k];

                g_y_0_zz_zzzz[k] = -g_y_0_z_zzzz[k] * ab_z + g_y_0_z_zzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 188 * ccomps * dcomps);

            auto g_z_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 189 * ccomps * dcomps);

            auto g_z_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 190 * ccomps * dcomps);

            auto g_z_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 191 * ccomps * dcomps);

            auto g_z_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 192 * ccomps * dcomps);

            auto g_z_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 193 * ccomps * dcomps);

            auto g_z_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxxx, g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxy, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxz, g_z_0_x_xxxzz, g_z_0_x_xxyy, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyz, g_z_0_x_xxyzz, g_z_0_x_xxzz, g_z_0_x_xxzzz, g_z_0_x_xyyy, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyz, g_z_0_x_xyyzz, g_z_0_x_xyzz, g_z_0_x_xyzzz, g_z_0_x_xzzz, g_z_0_x_xzzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz, g_z_0_xx_xxxx, g_z_0_xx_xxxy, g_z_0_xx_xxxz, g_z_0_xx_xxyy, g_z_0_xx_xxyz, g_z_0_xx_xxzz, g_z_0_xx_xyyy, g_z_0_xx_xyyz, g_z_0_xx_xyzz, g_z_0_xx_xzzz, g_z_0_xx_yyyy, g_z_0_xx_yyyz, g_z_0_xx_yyzz, g_z_0_xx_yzzz, g_z_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_xxxx[k] = -g_z_0_x_xxxx[k] * ab_x + g_z_0_x_xxxxx[k];

                g_z_0_xx_xxxy[k] = -g_z_0_x_xxxy[k] * ab_x + g_z_0_x_xxxxy[k];

                g_z_0_xx_xxxz[k] = -g_z_0_x_xxxz[k] * ab_x + g_z_0_x_xxxxz[k];

                g_z_0_xx_xxyy[k] = -g_z_0_x_xxyy[k] * ab_x + g_z_0_x_xxxyy[k];

                g_z_0_xx_xxyz[k] = -g_z_0_x_xxyz[k] * ab_x + g_z_0_x_xxxyz[k];

                g_z_0_xx_xxzz[k] = -g_z_0_x_xxzz[k] * ab_x + g_z_0_x_xxxzz[k];

                g_z_0_xx_xyyy[k] = -g_z_0_x_xyyy[k] * ab_x + g_z_0_x_xxyyy[k];

                g_z_0_xx_xyyz[k] = -g_z_0_x_xyyz[k] * ab_x + g_z_0_x_xxyyz[k];

                g_z_0_xx_xyzz[k] = -g_z_0_x_xyzz[k] * ab_x + g_z_0_x_xxyzz[k];

                g_z_0_xx_xzzz[k] = -g_z_0_x_xzzz[k] * ab_x + g_z_0_x_xxzzz[k];

                g_z_0_xx_yyyy[k] = -g_z_0_x_yyyy[k] * ab_x + g_z_0_x_xyyyy[k];

                g_z_0_xx_yyyz[k] = -g_z_0_x_yyyz[k] * ab_x + g_z_0_x_xyyyz[k];

                g_z_0_xx_yyzz[k] = -g_z_0_x_yyzz[k] * ab_x + g_z_0_x_xyyzz[k];

                g_z_0_xx_yzzz[k] = -g_z_0_x_yzzz[k] * ab_x + g_z_0_x_xyzzz[k];

                g_z_0_xx_zzzz[k] = -g_z_0_x_zzzz[k] * ab_x + g_z_0_x_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 195 * ccomps * dcomps);

            auto g_z_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 196 * ccomps * dcomps);

            auto g_z_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 197 * ccomps * dcomps);

            auto g_z_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 198 * ccomps * dcomps);

            auto g_z_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 199 * ccomps * dcomps);

            auto g_z_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 200 * ccomps * dcomps);

            auto g_z_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 201 * ccomps * dcomps);

            auto g_z_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 202 * ccomps * dcomps);

            auto g_z_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 203 * ccomps * dcomps);

            auto g_z_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 204 * ccomps * dcomps);

            auto g_z_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 205 * ccomps * dcomps);

            auto g_z_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 206 * ccomps * dcomps);

            auto g_z_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 207 * ccomps * dcomps);

            auto g_z_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 208 * ccomps * dcomps);

            auto g_z_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xy_xxxx, g_z_0_xy_xxxy, g_z_0_xy_xxxz, g_z_0_xy_xxyy, g_z_0_xy_xxyz, g_z_0_xy_xxzz, g_z_0_xy_xyyy, g_z_0_xy_xyyz, g_z_0_xy_xyzz, g_z_0_xy_xzzz, g_z_0_xy_yyyy, g_z_0_xy_yyyz, g_z_0_xy_yyzz, g_z_0_xy_yzzz, g_z_0_xy_zzzz, g_z_0_y_xxxx, g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxy, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxz, g_z_0_y_xxxzz, g_z_0_y_xxyy, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyz, g_z_0_y_xxyzz, g_z_0_y_xxzz, g_z_0_y_xxzzz, g_z_0_y_xyyy, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyz, g_z_0_y_xyyzz, g_z_0_y_xyzz, g_z_0_y_xyzzz, g_z_0_y_xzzz, g_z_0_y_xzzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_xxxx[k] = -g_z_0_y_xxxx[k] * ab_x + g_z_0_y_xxxxx[k];

                g_z_0_xy_xxxy[k] = -g_z_0_y_xxxy[k] * ab_x + g_z_0_y_xxxxy[k];

                g_z_0_xy_xxxz[k] = -g_z_0_y_xxxz[k] * ab_x + g_z_0_y_xxxxz[k];

                g_z_0_xy_xxyy[k] = -g_z_0_y_xxyy[k] * ab_x + g_z_0_y_xxxyy[k];

                g_z_0_xy_xxyz[k] = -g_z_0_y_xxyz[k] * ab_x + g_z_0_y_xxxyz[k];

                g_z_0_xy_xxzz[k] = -g_z_0_y_xxzz[k] * ab_x + g_z_0_y_xxxzz[k];

                g_z_0_xy_xyyy[k] = -g_z_0_y_xyyy[k] * ab_x + g_z_0_y_xxyyy[k];

                g_z_0_xy_xyyz[k] = -g_z_0_y_xyyz[k] * ab_x + g_z_0_y_xxyyz[k];

                g_z_0_xy_xyzz[k] = -g_z_0_y_xyzz[k] * ab_x + g_z_0_y_xxyzz[k];

                g_z_0_xy_xzzz[k] = -g_z_0_y_xzzz[k] * ab_x + g_z_0_y_xxzzz[k];

                g_z_0_xy_yyyy[k] = -g_z_0_y_yyyy[k] * ab_x + g_z_0_y_xyyyy[k];

                g_z_0_xy_yyyz[k] = -g_z_0_y_yyyz[k] * ab_x + g_z_0_y_xyyyz[k];

                g_z_0_xy_yyzz[k] = -g_z_0_y_yyzz[k] * ab_x + g_z_0_y_xyyzz[k];

                g_z_0_xy_yzzz[k] = -g_z_0_y_yzzz[k] * ab_x + g_z_0_y_xyzzz[k];

                g_z_0_xy_zzzz[k] = -g_z_0_y_zzzz[k] * ab_x + g_z_0_y_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 210 * ccomps * dcomps);

            auto g_z_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 211 * ccomps * dcomps);

            auto g_z_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 212 * ccomps * dcomps);

            auto g_z_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 213 * ccomps * dcomps);

            auto g_z_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 214 * ccomps * dcomps);

            auto g_z_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xxxx, g_z_0_xz_xxxy, g_z_0_xz_xxxz, g_z_0_xz_xxyy, g_z_0_xz_xxyz, g_z_0_xz_xxzz, g_z_0_xz_xyyy, g_z_0_xz_xyyz, g_z_0_xz_xyzz, g_z_0_xz_xzzz, g_z_0_xz_yyyy, g_z_0_xz_yyyz, g_z_0_xz_yyzz, g_z_0_xz_yzzz, g_z_0_xz_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxy, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxxzz, g_z_0_z_xxyy, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xxzzz, g_z_0_z_xyyy, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_xzzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_xxxx[k] = -g_z_0_z_xxxx[k] * ab_x + g_z_0_z_xxxxx[k];

                g_z_0_xz_xxxy[k] = -g_z_0_z_xxxy[k] * ab_x + g_z_0_z_xxxxy[k];

                g_z_0_xz_xxxz[k] = -g_z_0_z_xxxz[k] * ab_x + g_z_0_z_xxxxz[k];

                g_z_0_xz_xxyy[k] = -g_z_0_z_xxyy[k] * ab_x + g_z_0_z_xxxyy[k];

                g_z_0_xz_xxyz[k] = -g_z_0_z_xxyz[k] * ab_x + g_z_0_z_xxxyz[k];

                g_z_0_xz_xxzz[k] = -g_z_0_z_xxzz[k] * ab_x + g_z_0_z_xxxzz[k];

                g_z_0_xz_xyyy[k] = -g_z_0_z_xyyy[k] * ab_x + g_z_0_z_xxyyy[k];

                g_z_0_xz_xyyz[k] = -g_z_0_z_xyyz[k] * ab_x + g_z_0_z_xxyyz[k];

                g_z_0_xz_xyzz[k] = -g_z_0_z_xyzz[k] * ab_x + g_z_0_z_xxyzz[k];

                g_z_0_xz_xzzz[k] = -g_z_0_z_xzzz[k] * ab_x + g_z_0_z_xxzzz[k];

                g_z_0_xz_yyyy[k] = -g_z_0_z_yyyy[k] * ab_x + g_z_0_z_xyyyy[k];

                g_z_0_xz_yyyz[k] = -g_z_0_z_yyyz[k] * ab_x + g_z_0_z_xyyyz[k];

                g_z_0_xz_yyzz[k] = -g_z_0_z_yyzz[k] * ab_x + g_z_0_z_xyyzz[k];

                g_z_0_xz_yzzz[k] = -g_z_0_z_yzzz[k] * ab_x + g_z_0_z_xyzzz[k];

                g_z_0_xz_zzzz[k] = -g_z_0_z_zzzz[k] * ab_x + g_z_0_z_xzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxxx, g_z_0_y_xxxxy, g_z_0_y_xxxy, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyz, g_z_0_y_xxyzz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyz, g_z_0_y_xyyzz, g_z_0_y_xyzz, g_z_0_y_xyzzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyz, g_z_0_y_yyyzz, g_z_0_y_yyzz, g_z_0_y_yyzzz, g_z_0_y_yzzz, g_z_0_y_yzzzz, g_z_0_y_zzzz, g_z_0_yy_xxxx, g_z_0_yy_xxxy, g_z_0_yy_xxxz, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxzz, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyzz, g_z_0_yy_xzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyzz, g_z_0_yy_yzzz, g_z_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_xxxx[k] = -g_z_0_y_xxxx[k] * ab_y + g_z_0_y_xxxxy[k];

                g_z_0_yy_xxxy[k] = -g_z_0_y_xxxy[k] * ab_y + g_z_0_y_xxxyy[k];

                g_z_0_yy_xxxz[k] = -g_z_0_y_xxxz[k] * ab_y + g_z_0_y_xxxyz[k];

                g_z_0_yy_xxyy[k] = -g_z_0_y_xxyy[k] * ab_y + g_z_0_y_xxyyy[k];

                g_z_0_yy_xxyz[k] = -g_z_0_y_xxyz[k] * ab_y + g_z_0_y_xxyyz[k];

                g_z_0_yy_xxzz[k] = -g_z_0_y_xxzz[k] * ab_y + g_z_0_y_xxyzz[k];

                g_z_0_yy_xyyy[k] = -g_z_0_y_xyyy[k] * ab_y + g_z_0_y_xyyyy[k];

                g_z_0_yy_xyyz[k] = -g_z_0_y_xyyz[k] * ab_y + g_z_0_y_xyyyz[k];

                g_z_0_yy_xyzz[k] = -g_z_0_y_xyzz[k] * ab_y + g_z_0_y_xyyzz[k];

                g_z_0_yy_xzzz[k] = -g_z_0_y_xzzz[k] * ab_y + g_z_0_y_xyzzz[k];

                g_z_0_yy_yyyy[k] = -g_z_0_y_yyyy[k] * ab_y + g_z_0_y_yyyyy[k];

                g_z_0_yy_yyyz[k] = -g_z_0_y_yyyz[k] * ab_y + g_z_0_y_yyyyz[k];

                g_z_0_yy_yyzz[k] = -g_z_0_y_yyzz[k] * ab_y + g_z_0_y_yyyzz[k];

                g_z_0_yy_yzzz[k] = -g_z_0_y_yzzz[k] * ab_y + g_z_0_y_yyzzz[k];

                g_z_0_yy_zzzz[k] = -g_z_0_y_zzzz[k] * ab_y + g_z_0_y_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 244 * ccomps * dcomps);

            auto g_z_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 245 * ccomps * dcomps);

            auto g_z_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 246 * ccomps * dcomps);

            auto g_z_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 247 * ccomps * dcomps);

            auto g_z_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 248 * ccomps * dcomps);

            auto g_z_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 249 * ccomps * dcomps);

            auto g_z_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 250 * ccomps * dcomps);

            auto g_z_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xxxx, g_z_0_yz_xxxy, g_z_0_yz_xxxz, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxzz, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyzz, g_z_0_yz_xzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyzz, g_z_0_yz_yzzz, g_z_0_yz_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxxy, g_z_0_z_xxxy, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyz, g_z_0_z_yyyzz, g_z_0_z_yyzz, g_z_0_z_yyzzz, g_z_0_z_yzzz, g_z_0_z_yzzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_xxxx[k] = -g_z_0_z_xxxx[k] * ab_y + g_z_0_z_xxxxy[k];

                g_z_0_yz_xxxy[k] = -g_z_0_z_xxxy[k] * ab_y + g_z_0_z_xxxyy[k];

                g_z_0_yz_xxxz[k] = -g_z_0_z_xxxz[k] * ab_y + g_z_0_z_xxxyz[k];

                g_z_0_yz_xxyy[k] = -g_z_0_z_xxyy[k] * ab_y + g_z_0_z_xxyyy[k];

                g_z_0_yz_xxyz[k] = -g_z_0_z_xxyz[k] * ab_y + g_z_0_z_xxyyz[k];

                g_z_0_yz_xxzz[k] = -g_z_0_z_xxzz[k] * ab_y + g_z_0_z_xxyzz[k];

                g_z_0_yz_xyyy[k] = -g_z_0_z_xyyy[k] * ab_y + g_z_0_z_xyyyy[k];

                g_z_0_yz_xyyz[k] = -g_z_0_z_xyyz[k] * ab_y + g_z_0_z_xyyyz[k];

                g_z_0_yz_xyzz[k] = -g_z_0_z_xyzz[k] * ab_y + g_z_0_z_xyyzz[k];

                g_z_0_yz_xzzz[k] = -g_z_0_z_xzzz[k] * ab_y + g_z_0_z_xyzzz[k];

                g_z_0_yz_yyyy[k] = -g_z_0_z_yyyy[k] * ab_y + g_z_0_z_yyyyy[k];

                g_z_0_yz_yyyz[k] = -g_z_0_z_yyyz[k] * ab_y + g_z_0_z_yyyyz[k];

                g_z_0_yz_yyzz[k] = -g_z_0_z_yyzz[k] * ab_y + g_z_0_z_yyyzz[k];

                g_z_0_yz_yzzz[k] = -g_z_0_z_yzzz[k] * ab_y + g_z_0_z_yyzzz[k];

                g_z_0_yz_zzzz[k] = -g_z_0_z_zzzz[k] * ab_y + g_z_0_z_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxx, g_z_0_z_xxxxz, g_z_0_z_xxxy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxxzz, g_z_0_z_xxyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xxzzz, g_z_0_z_xyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_xzzzz, g_z_0_z_yyyy, g_z_0_z_yyyyz, g_z_0_z_yyyz, g_z_0_z_yyyzz, g_z_0_z_yyzz, g_z_0_z_yyzzz, g_z_0_z_yzzz, g_z_0_z_yzzzz, g_z_0_z_zzzz, g_z_0_z_zzzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz, g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy, g_z_xyyz, g_z_xyzz, g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_xxxx[k] = -g_z_xxxx[k] - g_z_0_z_xxxx[k] * ab_z + g_z_0_z_xxxxz[k];

                g_z_0_zz_xxxy[k] = -g_z_xxxy[k] - g_z_0_z_xxxy[k] * ab_z + g_z_0_z_xxxyz[k];

                g_z_0_zz_xxxz[k] = -g_z_xxxz[k] - g_z_0_z_xxxz[k] * ab_z + g_z_0_z_xxxzz[k];

                g_z_0_zz_xxyy[k] = -g_z_xxyy[k] - g_z_0_z_xxyy[k] * ab_z + g_z_0_z_xxyyz[k];

                g_z_0_zz_xxyz[k] = -g_z_xxyz[k] - g_z_0_z_xxyz[k] * ab_z + g_z_0_z_xxyzz[k];

                g_z_0_zz_xxzz[k] = -g_z_xxzz[k] - g_z_0_z_xxzz[k] * ab_z + g_z_0_z_xxzzz[k];

                g_z_0_zz_xyyy[k] = -g_z_xyyy[k] - g_z_0_z_xyyy[k] * ab_z + g_z_0_z_xyyyz[k];

                g_z_0_zz_xyyz[k] = -g_z_xyyz[k] - g_z_0_z_xyyz[k] * ab_z + g_z_0_z_xyyzz[k];

                g_z_0_zz_xyzz[k] = -g_z_xyzz[k] - g_z_0_z_xyzz[k] * ab_z + g_z_0_z_xyzzz[k];

                g_z_0_zz_xzzz[k] = -g_z_xzzz[k] - g_z_0_z_xzzz[k] * ab_z + g_z_0_z_xzzzz[k];

                g_z_0_zz_yyyy[k] = -g_z_yyyy[k] - g_z_0_z_yyyy[k] * ab_z + g_z_0_z_yyyyz[k];

                g_z_0_zz_yyyz[k] = -g_z_yyyz[k] - g_z_0_z_yyyz[k] * ab_z + g_z_0_z_yyyzz[k];

                g_z_0_zz_yyzz[k] = -g_z_yyzz[k] - g_z_0_z_yyzz[k] * ab_z + g_z_0_z_yyzzz[k];

                g_z_0_zz_yzzz[k] = -g_z_yzzz[k] - g_z_0_z_yzzz[k] * ab_z + g_z_0_z_yzzzz[k];

                g_z_0_zz_zzzz[k] = -g_z_zzzz[k] - g_z_0_z_zzzz[k] * ab_z + g_z_0_z_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

