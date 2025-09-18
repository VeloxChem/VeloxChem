#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_pgxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_pgxx,
                                              const size_t idx_geom_0010_sgxx,
                                              const size_t idx_geom_1010_sgxx,
                                              const size_t idx_geom_1010_shxx,
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
            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_0010_off = idx_geom_0010_sgxx + i * dcomps + j;

            auto g_0_0_x_0_0_xxxx = cbuffer.data(sg_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxy = cbuffer.data(sg_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxz = cbuffer.data(sg_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyy = cbuffer.data(sg_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyz = cbuffer.data(sg_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxzz = cbuffer.data(sg_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyy = cbuffer.data(sg_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyz = cbuffer.data(sg_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyzz = cbuffer.data(sg_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_0_xzzz = cbuffer.data(sg_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyy = cbuffer.data(sg_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyz = cbuffer.data(sg_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyzz = cbuffer.data(sg_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_0_yzzz = cbuffer.data(sg_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_0_zzzz = cbuffer.data(sg_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxx = cbuffer.data(sg_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxy = cbuffer.data(sg_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxz = cbuffer.data(sg_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyy = cbuffer.data(sg_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyz = cbuffer.data(sg_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxzz = cbuffer.data(sg_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyy = cbuffer.data(sg_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyz = cbuffer.data(sg_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyzz = cbuffer.data(sg_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_y_0_0_xzzz = cbuffer.data(sg_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyy = cbuffer.data(sg_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyz = cbuffer.data(sg_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyzz = cbuffer.data(sg_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_y_0_0_yzzz = cbuffer.data(sg_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_y_0_0_zzzz = cbuffer.data(sg_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxx = cbuffer.data(sg_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxy = cbuffer.data(sg_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxz = cbuffer.data(sg_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyy = cbuffer.data(sg_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyz = cbuffer.data(sg_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxzz = cbuffer.data(sg_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyy = cbuffer.data(sg_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyz = cbuffer.data(sg_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyzz = cbuffer.data(sg_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_z_0_0_xzzz = cbuffer.data(sg_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyy = cbuffer.data(sg_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyz = cbuffer.data(sg_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyzz = cbuffer.data(sg_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_z_0_0_yzzz = cbuffer.data(sg_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_z_0_0_zzzz = cbuffer.data(sg_geom_0010_off + 44 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_1010_off = idx_geom_1010_sgxx + i * dcomps + j;

            auto g_x_0_x_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_y_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_y_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_y_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_z_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_z_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_z_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 44 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 45 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 46 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 47 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 48 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 49 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 50 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 51 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 52 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 53 * ccomps * dcomps);

            auto g_y_0_x_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 54 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 55 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 56 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 57 * ccomps * dcomps);

            auto g_y_0_x_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 58 * ccomps * dcomps);

            auto g_y_0_x_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 59 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 60 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 61 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 62 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 63 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 64 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 65 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 66 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 67 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 68 * ccomps * dcomps);

            auto g_y_0_y_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 69 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 70 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 71 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 72 * ccomps * dcomps);

            auto g_y_0_y_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 73 * ccomps * dcomps);

            auto g_y_0_y_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 74 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 75 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 76 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 77 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 78 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 79 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 80 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 81 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 82 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 83 * ccomps * dcomps);

            auto g_y_0_z_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 84 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 85 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 86 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 87 * ccomps * dcomps);

            auto g_y_0_z_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 88 * ccomps * dcomps);

            auto g_y_0_z_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 89 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 90 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 91 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 92 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 93 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 94 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 95 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 96 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 97 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 98 * ccomps * dcomps);

            auto g_z_0_x_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 99 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 100 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 101 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 102 * ccomps * dcomps);

            auto g_z_0_x_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 103 * ccomps * dcomps);

            auto g_z_0_x_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 104 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 105 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 106 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 107 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 108 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 109 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 110 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 111 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 112 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 113 * ccomps * dcomps);

            auto g_z_0_y_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 114 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 115 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 116 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 117 * ccomps * dcomps);

            auto g_z_0_y_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 118 * ccomps * dcomps);

            auto g_z_0_y_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 119 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxx = cbuffer.data(sg_geom_1010_off + 120 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxy = cbuffer.data(sg_geom_1010_off + 121 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxz = cbuffer.data(sg_geom_1010_off + 122 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyy = cbuffer.data(sg_geom_1010_off + 123 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyz = cbuffer.data(sg_geom_1010_off + 124 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxzz = cbuffer.data(sg_geom_1010_off + 125 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyy = cbuffer.data(sg_geom_1010_off + 126 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyz = cbuffer.data(sg_geom_1010_off + 127 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyzz = cbuffer.data(sg_geom_1010_off + 128 * ccomps * dcomps);

            auto g_z_0_z_0_0_xzzz = cbuffer.data(sg_geom_1010_off + 129 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyy = cbuffer.data(sg_geom_1010_off + 130 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyz = cbuffer.data(sg_geom_1010_off + 131 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyzz = cbuffer.data(sg_geom_1010_off + 132 * ccomps * dcomps);

            auto g_z_0_z_0_0_yzzz = cbuffer.data(sg_geom_1010_off + 133 * ccomps * dcomps);

            auto g_z_0_z_0_0_zzzz = cbuffer.data(sg_geom_1010_off + 134 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_1010_off = idx_geom_1010_shxx + i * dcomps + j;

            auto g_x_0_x_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_y_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_y_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_y_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_z_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_z_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_z_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 62 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 63 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 64 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 65 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 66 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 67 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 68 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 69 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 70 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 71 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 72 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 73 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 74 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 75 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 76 * ccomps * dcomps);

            auto g_y_0_x_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 77 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 78 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 79 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 80 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 81 * ccomps * dcomps);

            auto g_y_0_x_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 82 * ccomps * dcomps);

            auto g_y_0_x_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 83 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 84 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 85 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 86 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 87 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 88 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 89 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 90 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 91 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 92 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 93 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 94 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 95 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 96 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 97 * ccomps * dcomps);

            auto g_y_0_y_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 98 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 99 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 100 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 101 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 102 * ccomps * dcomps);

            auto g_y_0_y_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 103 * ccomps * dcomps);

            auto g_y_0_y_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 104 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 105 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 106 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 107 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 108 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 109 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 110 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 111 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 112 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 113 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 114 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 115 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 116 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 117 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 118 * ccomps * dcomps);

            auto g_y_0_z_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 119 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 120 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 121 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 122 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 123 * ccomps * dcomps);

            auto g_y_0_z_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 124 * ccomps * dcomps);

            auto g_y_0_z_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 125 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 126 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 127 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 128 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 129 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 130 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 131 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 132 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 133 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 134 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 135 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 136 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 137 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 138 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 139 * ccomps * dcomps);

            auto g_z_0_x_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 140 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 141 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 142 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 143 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 144 * ccomps * dcomps);

            auto g_z_0_x_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 145 * ccomps * dcomps);

            auto g_z_0_x_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 146 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 147 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 148 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 149 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 150 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 151 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 152 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 153 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 154 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 155 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 156 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 157 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 158 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 159 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 160 * ccomps * dcomps);

            auto g_z_0_y_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 161 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 162 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 163 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 164 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 165 * ccomps * dcomps);

            auto g_z_0_y_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 166 * ccomps * dcomps);

            auto g_z_0_y_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 167 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxx = cbuffer.data(sh_geom_1010_off + 168 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxy = cbuffer.data(sh_geom_1010_off + 169 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxz = cbuffer.data(sh_geom_1010_off + 170 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyy = cbuffer.data(sh_geom_1010_off + 171 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyz = cbuffer.data(sh_geom_1010_off + 172 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxzz = cbuffer.data(sh_geom_1010_off + 173 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyy = cbuffer.data(sh_geom_1010_off + 174 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyz = cbuffer.data(sh_geom_1010_off + 175 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyzz = cbuffer.data(sh_geom_1010_off + 176 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxzzz = cbuffer.data(sh_geom_1010_off + 177 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyy = cbuffer.data(sh_geom_1010_off + 178 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyz = cbuffer.data(sh_geom_1010_off + 179 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyzz = cbuffer.data(sh_geom_1010_off + 180 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyzzz = cbuffer.data(sh_geom_1010_off + 181 * ccomps * dcomps);

            auto g_z_0_z_0_0_xzzzz = cbuffer.data(sh_geom_1010_off + 182 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyy = cbuffer.data(sh_geom_1010_off + 183 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyz = cbuffer.data(sh_geom_1010_off + 184 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyzz = cbuffer.data(sh_geom_1010_off + 185 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyzzz = cbuffer.data(sh_geom_1010_off + 186 * ccomps * dcomps);

            auto g_z_0_z_0_0_yzzzz = cbuffer.data(sh_geom_1010_off + 187 * ccomps * dcomps);

            auto g_z_0_z_0_0_zzzzz = cbuffer.data(sh_geom_1010_off + 188 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pgxx

            const auto pg_geom_1010_off = idx_geom_1010_pgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxx, g_0_0_x_0_0_xxxy, g_0_0_x_0_0_xxxz, g_0_0_x_0_0_xxyy, g_0_0_x_0_0_xxyz, g_0_0_x_0_0_xxzz, g_0_0_x_0_0_xyyy, g_0_0_x_0_0_xyyz, g_0_0_x_0_0_xyzz, g_0_0_x_0_0_xzzz, g_0_0_x_0_0_yyyy, g_0_0_x_0_0_yyyz, g_0_0_x_0_0_yyzz, g_0_0_x_0_0_yzzz, g_0_0_x_0_0_zzzz, g_x_0_x_0_0_xxxx, g_x_0_x_0_0_xxxxx, g_x_0_x_0_0_xxxxy, g_x_0_x_0_0_xxxxz, g_x_0_x_0_0_xxxy, g_x_0_x_0_0_xxxyy, g_x_0_x_0_0_xxxyz, g_x_0_x_0_0_xxxz, g_x_0_x_0_0_xxxzz, g_x_0_x_0_0_xxyy, g_x_0_x_0_0_xxyyy, g_x_0_x_0_0_xxyyz, g_x_0_x_0_0_xxyz, g_x_0_x_0_0_xxyzz, g_x_0_x_0_0_xxzz, g_x_0_x_0_0_xxzzz, g_x_0_x_0_0_xyyy, g_x_0_x_0_0_xyyyy, g_x_0_x_0_0_xyyyz, g_x_0_x_0_0_xyyz, g_x_0_x_0_0_xyyzz, g_x_0_x_0_0_xyzz, g_x_0_x_0_0_xyzzz, g_x_0_x_0_0_xzzz, g_x_0_x_0_0_xzzzz, g_x_0_x_0_0_yyyy, g_x_0_x_0_0_yyyz, g_x_0_x_0_0_yyzz, g_x_0_x_0_0_yzzz, g_x_0_x_0_0_zzzz, g_x_0_x_0_x_xxxx, g_x_0_x_0_x_xxxy, g_x_0_x_0_x_xxxz, g_x_0_x_0_x_xxyy, g_x_0_x_0_x_xxyz, g_x_0_x_0_x_xxzz, g_x_0_x_0_x_xyyy, g_x_0_x_0_x_xyyz, g_x_0_x_0_x_xyzz, g_x_0_x_0_x_xzzz, g_x_0_x_0_x_yyyy, g_x_0_x_0_x_yyyz, g_x_0_x_0_x_yyzz, g_x_0_x_0_x_yzzz, g_x_0_x_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_x_xxxx[k] = -g_0_0_x_0_0_xxxx[k] - g_x_0_x_0_0_xxxx[k] * ab_x + g_x_0_x_0_0_xxxxx[k];

                g_x_0_x_0_x_xxxy[k] = -g_0_0_x_0_0_xxxy[k] - g_x_0_x_0_0_xxxy[k] * ab_x + g_x_0_x_0_0_xxxxy[k];

                g_x_0_x_0_x_xxxz[k] = -g_0_0_x_0_0_xxxz[k] - g_x_0_x_0_0_xxxz[k] * ab_x + g_x_0_x_0_0_xxxxz[k];

                g_x_0_x_0_x_xxyy[k] = -g_0_0_x_0_0_xxyy[k] - g_x_0_x_0_0_xxyy[k] * ab_x + g_x_0_x_0_0_xxxyy[k];

                g_x_0_x_0_x_xxyz[k] = -g_0_0_x_0_0_xxyz[k] - g_x_0_x_0_0_xxyz[k] * ab_x + g_x_0_x_0_0_xxxyz[k];

                g_x_0_x_0_x_xxzz[k] = -g_0_0_x_0_0_xxzz[k] - g_x_0_x_0_0_xxzz[k] * ab_x + g_x_0_x_0_0_xxxzz[k];

                g_x_0_x_0_x_xyyy[k] = -g_0_0_x_0_0_xyyy[k] - g_x_0_x_0_0_xyyy[k] * ab_x + g_x_0_x_0_0_xxyyy[k];

                g_x_0_x_0_x_xyyz[k] = -g_0_0_x_0_0_xyyz[k] - g_x_0_x_0_0_xyyz[k] * ab_x + g_x_0_x_0_0_xxyyz[k];

                g_x_0_x_0_x_xyzz[k] = -g_0_0_x_0_0_xyzz[k] - g_x_0_x_0_0_xyzz[k] * ab_x + g_x_0_x_0_0_xxyzz[k];

                g_x_0_x_0_x_xzzz[k] = -g_0_0_x_0_0_xzzz[k] - g_x_0_x_0_0_xzzz[k] * ab_x + g_x_0_x_0_0_xxzzz[k];

                g_x_0_x_0_x_yyyy[k] = -g_0_0_x_0_0_yyyy[k] - g_x_0_x_0_0_yyyy[k] * ab_x + g_x_0_x_0_0_xyyyy[k];

                g_x_0_x_0_x_yyyz[k] = -g_0_0_x_0_0_yyyz[k] - g_x_0_x_0_0_yyyz[k] * ab_x + g_x_0_x_0_0_xyyyz[k];

                g_x_0_x_0_x_yyzz[k] = -g_0_0_x_0_0_yyzz[k] - g_x_0_x_0_0_yyzz[k] * ab_x + g_x_0_x_0_0_xyyzz[k];

                g_x_0_x_0_x_yzzz[k] = -g_0_0_x_0_0_yzzz[k] - g_x_0_x_0_0_yzzz[k] * ab_x + g_x_0_x_0_0_xyzzz[k];

                g_x_0_x_0_x_zzzz[k] = -g_0_0_x_0_0_zzzz[k] - g_x_0_x_0_0_zzzz[k] * ab_x + g_x_0_x_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xxxx, g_x_0_x_0_0_xxxxy, g_x_0_x_0_0_xxxy, g_x_0_x_0_0_xxxyy, g_x_0_x_0_0_xxxyz, g_x_0_x_0_0_xxxz, g_x_0_x_0_0_xxyy, g_x_0_x_0_0_xxyyy, g_x_0_x_0_0_xxyyz, g_x_0_x_0_0_xxyz, g_x_0_x_0_0_xxyzz, g_x_0_x_0_0_xxzz, g_x_0_x_0_0_xyyy, g_x_0_x_0_0_xyyyy, g_x_0_x_0_0_xyyyz, g_x_0_x_0_0_xyyz, g_x_0_x_0_0_xyyzz, g_x_0_x_0_0_xyzz, g_x_0_x_0_0_xyzzz, g_x_0_x_0_0_xzzz, g_x_0_x_0_0_yyyy, g_x_0_x_0_0_yyyyy, g_x_0_x_0_0_yyyyz, g_x_0_x_0_0_yyyz, g_x_0_x_0_0_yyyzz, g_x_0_x_0_0_yyzz, g_x_0_x_0_0_yyzzz, g_x_0_x_0_0_yzzz, g_x_0_x_0_0_yzzzz, g_x_0_x_0_0_zzzz, g_x_0_x_0_y_xxxx, g_x_0_x_0_y_xxxy, g_x_0_x_0_y_xxxz, g_x_0_x_0_y_xxyy, g_x_0_x_0_y_xxyz, g_x_0_x_0_y_xxzz, g_x_0_x_0_y_xyyy, g_x_0_x_0_y_xyyz, g_x_0_x_0_y_xyzz, g_x_0_x_0_y_xzzz, g_x_0_x_0_y_yyyy, g_x_0_x_0_y_yyyz, g_x_0_x_0_y_yyzz, g_x_0_x_0_y_yzzz, g_x_0_x_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_y_xxxx[k] = -g_x_0_x_0_0_xxxx[k] * ab_y + g_x_0_x_0_0_xxxxy[k];

                g_x_0_x_0_y_xxxy[k] = -g_x_0_x_0_0_xxxy[k] * ab_y + g_x_0_x_0_0_xxxyy[k];

                g_x_0_x_0_y_xxxz[k] = -g_x_0_x_0_0_xxxz[k] * ab_y + g_x_0_x_0_0_xxxyz[k];

                g_x_0_x_0_y_xxyy[k] = -g_x_0_x_0_0_xxyy[k] * ab_y + g_x_0_x_0_0_xxyyy[k];

                g_x_0_x_0_y_xxyz[k] = -g_x_0_x_0_0_xxyz[k] * ab_y + g_x_0_x_0_0_xxyyz[k];

                g_x_0_x_0_y_xxzz[k] = -g_x_0_x_0_0_xxzz[k] * ab_y + g_x_0_x_0_0_xxyzz[k];

                g_x_0_x_0_y_xyyy[k] = -g_x_0_x_0_0_xyyy[k] * ab_y + g_x_0_x_0_0_xyyyy[k];

                g_x_0_x_0_y_xyyz[k] = -g_x_0_x_0_0_xyyz[k] * ab_y + g_x_0_x_0_0_xyyyz[k];

                g_x_0_x_0_y_xyzz[k] = -g_x_0_x_0_0_xyzz[k] * ab_y + g_x_0_x_0_0_xyyzz[k];

                g_x_0_x_0_y_xzzz[k] = -g_x_0_x_0_0_xzzz[k] * ab_y + g_x_0_x_0_0_xyzzz[k];

                g_x_0_x_0_y_yyyy[k] = -g_x_0_x_0_0_yyyy[k] * ab_y + g_x_0_x_0_0_yyyyy[k];

                g_x_0_x_0_y_yyyz[k] = -g_x_0_x_0_0_yyyz[k] * ab_y + g_x_0_x_0_0_yyyyz[k];

                g_x_0_x_0_y_yyzz[k] = -g_x_0_x_0_0_yyzz[k] * ab_y + g_x_0_x_0_0_yyyzz[k];

                g_x_0_x_0_y_yzzz[k] = -g_x_0_x_0_0_yzzz[k] * ab_y + g_x_0_x_0_0_yyzzz[k];

                g_x_0_x_0_y_zzzz[k] = -g_x_0_x_0_0_zzzz[k] * ab_y + g_x_0_x_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xxxx, g_x_0_x_0_0_xxxxz, g_x_0_x_0_0_xxxy, g_x_0_x_0_0_xxxyz, g_x_0_x_0_0_xxxz, g_x_0_x_0_0_xxxzz, g_x_0_x_0_0_xxyy, g_x_0_x_0_0_xxyyz, g_x_0_x_0_0_xxyz, g_x_0_x_0_0_xxyzz, g_x_0_x_0_0_xxzz, g_x_0_x_0_0_xxzzz, g_x_0_x_0_0_xyyy, g_x_0_x_0_0_xyyyz, g_x_0_x_0_0_xyyz, g_x_0_x_0_0_xyyzz, g_x_0_x_0_0_xyzz, g_x_0_x_0_0_xyzzz, g_x_0_x_0_0_xzzz, g_x_0_x_0_0_xzzzz, g_x_0_x_0_0_yyyy, g_x_0_x_0_0_yyyyz, g_x_0_x_0_0_yyyz, g_x_0_x_0_0_yyyzz, g_x_0_x_0_0_yyzz, g_x_0_x_0_0_yyzzz, g_x_0_x_0_0_yzzz, g_x_0_x_0_0_yzzzz, g_x_0_x_0_0_zzzz, g_x_0_x_0_0_zzzzz, g_x_0_x_0_z_xxxx, g_x_0_x_0_z_xxxy, g_x_0_x_0_z_xxxz, g_x_0_x_0_z_xxyy, g_x_0_x_0_z_xxyz, g_x_0_x_0_z_xxzz, g_x_0_x_0_z_xyyy, g_x_0_x_0_z_xyyz, g_x_0_x_0_z_xyzz, g_x_0_x_0_z_xzzz, g_x_0_x_0_z_yyyy, g_x_0_x_0_z_yyyz, g_x_0_x_0_z_yyzz, g_x_0_x_0_z_yzzz, g_x_0_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_z_xxxx[k] = -g_x_0_x_0_0_xxxx[k] * ab_z + g_x_0_x_0_0_xxxxz[k];

                g_x_0_x_0_z_xxxy[k] = -g_x_0_x_0_0_xxxy[k] * ab_z + g_x_0_x_0_0_xxxyz[k];

                g_x_0_x_0_z_xxxz[k] = -g_x_0_x_0_0_xxxz[k] * ab_z + g_x_0_x_0_0_xxxzz[k];

                g_x_0_x_0_z_xxyy[k] = -g_x_0_x_0_0_xxyy[k] * ab_z + g_x_0_x_0_0_xxyyz[k];

                g_x_0_x_0_z_xxyz[k] = -g_x_0_x_0_0_xxyz[k] * ab_z + g_x_0_x_0_0_xxyzz[k];

                g_x_0_x_0_z_xxzz[k] = -g_x_0_x_0_0_xxzz[k] * ab_z + g_x_0_x_0_0_xxzzz[k];

                g_x_0_x_0_z_xyyy[k] = -g_x_0_x_0_0_xyyy[k] * ab_z + g_x_0_x_0_0_xyyyz[k];

                g_x_0_x_0_z_xyyz[k] = -g_x_0_x_0_0_xyyz[k] * ab_z + g_x_0_x_0_0_xyyzz[k];

                g_x_0_x_0_z_xyzz[k] = -g_x_0_x_0_0_xyzz[k] * ab_z + g_x_0_x_0_0_xyzzz[k];

                g_x_0_x_0_z_xzzz[k] = -g_x_0_x_0_0_xzzz[k] * ab_z + g_x_0_x_0_0_xzzzz[k];

                g_x_0_x_0_z_yyyy[k] = -g_x_0_x_0_0_yyyy[k] * ab_z + g_x_0_x_0_0_yyyyz[k];

                g_x_0_x_0_z_yyyz[k] = -g_x_0_x_0_0_yyyz[k] * ab_z + g_x_0_x_0_0_yyyzz[k];

                g_x_0_x_0_z_yyzz[k] = -g_x_0_x_0_0_yyzz[k] * ab_z + g_x_0_x_0_0_yyzzz[k];

                g_x_0_x_0_z_yzzz[k] = -g_x_0_x_0_0_yzzz[k] * ab_z + g_x_0_x_0_0_yzzzz[k];

                g_x_0_x_0_z_zzzz[k] = -g_x_0_x_0_0_zzzz[k] * ab_z + g_x_0_x_0_0_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_y_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_y_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_y_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxx, g_0_0_y_0_0_xxxy, g_0_0_y_0_0_xxxz, g_0_0_y_0_0_xxyy, g_0_0_y_0_0_xxyz, g_0_0_y_0_0_xxzz, g_0_0_y_0_0_xyyy, g_0_0_y_0_0_xyyz, g_0_0_y_0_0_xyzz, g_0_0_y_0_0_xzzz, g_0_0_y_0_0_yyyy, g_0_0_y_0_0_yyyz, g_0_0_y_0_0_yyzz, g_0_0_y_0_0_yzzz, g_0_0_y_0_0_zzzz, g_x_0_y_0_0_xxxx, g_x_0_y_0_0_xxxxx, g_x_0_y_0_0_xxxxy, g_x_0_y_0_0_xxxxz, g_x_0_y_0_0_xxxy, g_x_0_y_0_0_xxxyy, g_x_0_y_0_0_xxxyz, g_x_0_y_0_0_xxxz, g_x_0_y_0_0_xxxzz, g_x_0_y_0_0_xxyy, g_x_0_y_0_0_xxyyy, g_x_0_y_0_0_xxyyz, g_x_0_y_0_0_xxyz, g_x_0_y_0_0_xxyzz, g_x_0_y_0_0_xxzz, g_x_0_y_0_0_xxzzz, g_x_0_y_0_0_xyyy, g_x_0_y_0_0_xyyyy, g_x_0_y_0_0_xyyyz, g_x_0_y_0_0_xyyz, g_x_0_y_0_0_xyyzz, g_x_0_y_0_0_xyzz, g_x_0_y_0_0_xyzzz, g_x_0_y_0_0_xzzz, g_x_0_y_0_0_xzzzz, g_x_0_y_0_0_yyyy, g_x_0_y_0_0_yyyz, g_x_0_y_0_0_yyzz, g_x_0_y_0_0_yzzz, g_x_0_y_0_0_zzzz, g_x_0_y_0_x_xxxx, g_x_0_y_0_x_xxxy, g_x_0_y_0_x_xxxz, g_x_0_y_0_x_xxyy, g_x_0_y_0_x_xxyz, g_x_0_y_0_x_xxzz, g_x_0_y_0_x_xyyy, g_x_0_y_0_x_xyyz, g_x_0_y_0_x_xyzz, g_x_0_y_0_x_xzzz, g_x_0_y_0_x_yyyy, g_x_0_y_0_x_yyyz, g_x_0_y_0_x_yyzz, g_x_0_y_0_x_yzzz, g_x_0_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_x_xxxx[k] = -g_0_0_y_0_0_xxxx[k] - g_x_0_y_0_0_xxxx[k] * ab_x + g_x_0_y_0_0_xxxxx[k];

                g_x_0_y_0_x_xxxy[k] = -g_0_0_y_0_0_xxxy[k] - g_x_0_y_0_0_xxxy[k] * ab_x + g_x_0_y_0_0_xxxxy[k];

                g_x_0_y_0_x_xxxz[k] = -g_0_0_y_0_0_xxxz[k] - g_x_0_y_0_0_xxxz[k] * ab_x + g_x_0_y_0_0_xxxxz[k];

                g_x_0_y_0_x_xxyy[k] = -g_0_0_y_0_0_xxyy[k] - g_x_0_y_0_0_xxyy[k] * ab_x + g_x_0_y_0_0_xxxyy[k];

                g_x_0_y_0_x_xxyz[k] = -g_0_0_y_0_0_xxyz[k] - g_x_0_y_0_0_xxyz[k] * ab_x + g_x_0_y_0_0_xxxyz[k];

                g_x_0_y_0_x_xxzz[k] = -g_0_0_y_0_0_xxzz[k] - g_x_0_y_0_0_xxzz[k] * ab_x + g_x_0_y_0_0_xxxzz[k];

                g_x_0_y_0_x_xyyy[k] = -g_0_0_y_0_0_xyyy[k] - g_x_0_y_0_0_xyyy[k] * ab_x + g_x_0_y_0_0_xxyyy[k];

                g_x_0_y_0_x_xyyz[k] = -g_0_0_y_0_0_xyyz[k] - g_x_0_y_0_0_xyyz[k] * ab_x + g_x_0_y_0_0_xxyyz[k];

                g_x_0_y_0_x_xyzz[k] = -g_0_0_y_0_0_xyzz[k] - g_x_0_y_0_0_xyzz[k] * ab_x + g_x_0_y_0_0_xxyzz[k];

                g_x_0_y_0_x_xzzz[k] = -g_0_0_y_0_0_xzzz[k] - g_x_0_y_0_0_xzzz[k] * ab_x + g_x_0_y_0_0_xxzzz[k];

                g_x_0_y_0_x_yyyy[k] = -g_0_0_y_0_0_yyyy[k] - g_x_0_y_0_0_yyyy[k] * ab_x + g_x_0_y_0_0_xyyyy[k];

                g_x_0_y_0_x_yyyz[k] = -g_0_0_y_0_0_yyyz[k] - g_x_0_y_0_0_yyyz[k] * ab_x + g_x_0_y_0_0_xyyyz[k];

                g_x_0_y_0_x_yyzz[k] = -g_0_0_y_0_0_yyzz[k] - g_x_0_y_0_0_yyzz[k] * ab_x + g_x_0_y_0_0_xyyzz[k];

                g_x_0_y_0_x_yzzz[k] = -g_0_0_y_0_0_yzzz[k] - g_x_0_y_0_0_yzzz[k] * ab_x + g_x_0_y_0_0_xyzzz[k];

                g_x_0_y_0_x_zzzz[k] = -g_0_0_y_0_0_zzzz[k] - g_x_0_y_0_0_zzzz[k] * ab_x + g_x_0_y_0_0_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_y_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_y_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_y_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xxxx, g_x_0_y_0_0_xxxxy, g_x_0_y_0_0_xxxy, g_x_0_y_0_0_xxxyy, g_x_0_y_0_0_xxxyz, g_x_0_y_0_0_xxxz, g_x_0_y_0_0_xxyy, g_x_0_y_0_0_xxyyy, g_x_0_y_0_0_xxyyz, g_x_0_y_0_0_xxyz, g_x_0_y_0_0_xxyzz, g_x_0_y_0_0_xxzz, g_x_0_y_0_0_xyyy, g_x_0_y_0_0_xyyyy, g_x_0_y_0_0_xyyyz, g_x_0_y_0_0_xyyz, g_x_0_y_0_0_xyyzz, g_x_0_y_0_0_xyzz, g_x_0_y_0_0_xyzzz, g_x_0_y_0_0_xzzz, g_x_0_y_0_0_yyyy, g_x_0_y_0_0_yyyyy, g_x_0_y_0_0_yyyyz, g_x_0_y_0_0_yyyz, g_x_0_y_0_0_yyyzz, g_x_0_y_0_0_yyzz, g_x_0_y_0_0_yyzzz, g_x_0_y_0_0_yzzz, g_x_0_y_0_0_yzzzz, g_x_0_y_0_0_zzzz, g_x_0_y_0_y_xxxx, g_x_0_y_0_y_xxxy, g_x_0_y_0_y_xxxz, g_x_0_y_0_y_xxyy, g_x_0_y_0_y_xxyz, g_x_0_y_0_y_xxzz, g_x_0_y_0_y_xyyy, g_x_0_y_0_y_xyyz, g_x_0_y_0_y_xyzz, g_x_0_y_0_y_xzzz, g_x_0_y_0_y_yyyy, g_x_0_y_0_y_yyyz, g_x_0_y_0_y_yyzz, g_x_0_y_0_y_yzzz, g_x_0_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_y_xxxx[k] = -g_x_0_y_0_0_xxxx[k] * ab_y + g_x_0_y_0_0_xxxxy[k];

                g_x_0_y_0_y_xxxy[k] = -g_x_0_y_0_0_xxxy[k] * ab_y + g_x_0_y_0_0_xxxyy[k];

                g_x_0_y_0_y_xxxz[k] = -g_x_0_y_0_0_xxxz[k] * ab_y + g_x_0_y_0_0_xxxyz[k];

                g_x_0_y_0_y_xxyy[k] = -g_x_0_y_0_0_xxyy[k] * ab_y + g_x_0_y_0_0_xxyyy[k];

                g_x_0_y_0_y_xxyz[k] = -g_x_0_y_0_0_xxyz[k] * ab_y + g_x_0_y_0_0_xxyyz[k];

                g_x_0_y_0_y_xxzz[k] = -g_x_0_y_0_0_xxzz[k] * ab_y + g_x_0_y_0_0_xxyzz[k];

                g_x_0_y_0_y_xyyy[k] = -g_x_0_y_0_0_xyyy[k] * ab_y + g_x_0_y_0_0_xyyyy[k];

                g_x_0_y_0_y_xyyz[k] = -g_x_0_y_0_0_xyyz[k] * ab_y + g_x_0_y_0_0_xyyyz[k];

                g_x_0_y_0_y_xyzz[k] = -g_x_0_y_0_0_xyzz[k] * ab_y + g_x_0_y_0_0_xyyzz[k];

                g_x_0_y_0_y_xzzz[k] = -g_x_0_y_0_0_xzzz[k] * ab_y + g_x_0_y_0_0_xyzzz[k];

                g_x_0_y_0_y_yyyy[k] = -g_x_0_y_0_0_yyyy[k] * ab_y + g_x_0_y_0_0_yyyyy[k];

                g_x_0_y_0_y_yyyz[k] = -g_x_0_y_0_0_yyyz[k] * ab_y + g_x_0_y_0_0_yyyyz[k];

                g_x_0_y_0_y_yyzz[k] = -g_x_0_y_0_0_yyzz[k] * ab_y + g_x_0_y_0_0_yyyzz[k];

                g_x_0_y_0_y_yzzz[k] = -g_x_0_y_0_0_yzzz[k] * ab_y + g_x_0_y_0_0_yyzzz[k];

                g_x_0_y_0_y_zzzz[k] = -g_x_0_y_0_0_zzzz[k] * ab_y + g_x_0_y_0_0_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_y_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xxxx, g_x_0_y_0_0_xxxxz, g_x_0_y_0_0_xxxy, g_x_0_y_0_0_xxxyz, g_x_0_y_0_0_xxxz, g_x_0_y_0_0_xxxzz, g_x_0_y_0_0_xxyy, g_x_0_y_0_0_xxyyz, g_x_0_y_0_0_xxyz, g_x_0_y_0_0_xxyzz, g_x_0_y_0_0_xxzz, g_x_0_y_0_0_xxzzz, g_x_0_y_0_0_xyyy, g_x_0_y_0_0_xyyyz, g_x_0_y_0_0_xyyz, g_x_0_y_0_0_xyyzz, g_x_0_y_0_0_xyzz, g_x_0_y_0_0_xyzzz, g_x_0_y_0_0_xzzz, g_x_0_y_0_0_xzzzz, g_x_0_y_0_0_yyyy, g_x_0_y_0_0_yyyyz, g_x_0_y_0_0_yyyz, g_x_0_y_0_0_yyyzz, g_x_0_y_0_0_yyzz, g_x_0_y_0_0_yyzzz, g_x_0_y_0_0_yzzz, g_x_0_y_0_0_yzzzz, g_x_0_y_0_0_zzzz, g_x_0_y_0_0_zzzzz, g_x_0_y_0_z_xxxx, g_x_0_y_0_z_xxxy, g_x_0_y_0_z_xxxz, g_x_0_y_0_z_xxyy, g_x_0_y_0_z_xxyz, g_x_0_y_0_z_xxzz, g_x_0_y_0_z_xyyy, g_x_0_y_0_z_xyyz, g_x_0_y_0_z_xyzz, g_x_0_y_0_z_xzzz, g_x_0_y_0_z_yyyy, g_x_0_y_0_z_yyyz, g_x_0_y_0_z_yyzz, g_x_0_y_0_z_yzzz, g_x_0_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_z_xxxx[k] = -g_x_0_y_0_0_xxxx[k] * ab_z + g_x_0_y_0_0_xxxxz[k];

                g_x_0_y_0_z_xxxy[k] = -g_x_0_y_0_0_xxxy[k] * ab_z + g_x_0_y_0_0_xxxyz[k];

                g_x_0_y_0_z_xxxz[k] = -g_x_0_y_0_0_xxxz[k] * ab_z + g_x_0_y_0_0_xxxzz[k];

                g_x_0_y_0_z_xxyy[k] = -g_x_0_y_0_0_xxyy[k] * ab_z + g_x_0_y_0_0_xxyyz[k];

                g_x_0_y_0_z_xxyz[k] = -g_x_0_y_0_0_xxyz[k] * ab_z + g_x_0_y_0_0_xxyzz[k];

                g_x_0_y_0_z_xxzz[k] = -g_x_0_y_0_0_xxzz[k] * ab_z + g_x_0_y_0_0_xxzzz[k];

                g_x_0_y_0_z_xyyy[k] = -g_x_0_y_0_0_xyyy[k] * ab_z + g_x_0_y_0_0_xyyyz[k];

                g_x_0_y_0_z_xyyz[k] = -g_x_0_y_0_0_xyyz[k] * ab_z + g_x_0_y_0_0_xyyzz[k];

                g_x_0_y_0_z_xyzz[k] = -g_x_0_y_0_0_xyzz[k] * ab_z + g_x_0_y_0_0_xyzzz[k];

                g_x_0_y_0_z_xzzz[k] = -g_x_0_y_0_0_xzzz[k] * ab_z + g_x_0_y_0_0_xzzzz[k];

                g_x_0_y_0_z_yyyy[k] = -g_x_0_y_0_0_yyyy[k] * ab_z + g_x_0_y_0_0_yyyyz[k];

                g_x_0_y_0_z_yyyz[k] = -g_x_0_y_0_0_yyyz[k] * ab_z + g_x_0_y_0_0_yyyzz[k];

                g_x_0_y_0_z_yyzz[k] = -g_x_0_y_0_0_yyzz[k] * ab_z + g_x_0_y_0_0_yyzzz[k];

                g_x_0_y_0_z_yzzz[k] = -g_x_0_y_0_0_yzzz[k] * ab_z + g_x_0_y_0_0_yzzzz[k];

                g_x_0_y_0_z_zzzz[k] = -g_x_0_y_0_0_zzzz[k] * ab_z + g_x_0_y_0_0_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_z_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_z_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_z_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxx, g_0_0_z_0_0_xxxy, g_0_0_z_0_0_xxxz, g_0_0_z_0_0_xxyy, g_0_0_z_0_0_xxyz, g_0_0_z_0_0_xxzz, g_0_0_z_0_0_xyyy, g_0_0_z_0_0_xyyz, g_0_0_z_0_0_xyzz, g_0_0_z_0_0_xzzz, g_0_0_z_0_0_yyyy, g_0_0_z_0_0_yyyz, g_0_0_z_0_0_yyzz, g_0_0_z_0_0_yzzz, g_0_0_z_0_0_zzzz, g_x_0_z_0_0_xxxx, g_x_0_z_0_0_xxxxx, g_x_0_z_0_0_xxxxy, g_x_0_z_0_0_xxxxz, g_x_0_z_0_0_xxxy, g_x_0_z_0_0_xxxyy, g_x_0_z_0_0_xxxyz, g_x_0_z_0_0_xxxz, g_x_0_z_0_0_xxxzz, g_x_0_z_0_0_xxyy, g_x_0_z_0_0_xxyyy, g_x_0_z_0_0_xxyyz, g_x_0_z_0_0_xxyz, g_x_0_z_0_0_xxyzz, g_x_0_z_0_0_xxzz, g_x_0_z_0_0_xxzzz, g_x_0_z_0_0_xyyy, g_x_0_z_0_0_xyyyy, g_x_0_z_0_0_xyyyz, g_x_0_z_0_0_xyyz, g_x_0_z_0_0_xyyzz, g_x_0_z_0_0_xyzz, g_x_0_z_0_0_xyzzz, g_x_0_z_0_0_xzzz, g_x_0_z_0_0_xzzzz, g_x_0_z_0_0_yyyy, g_x_0_z_0_0_yyyz, g_x_0_z_0_0_yyzz, g_x_0_z_0_0_yzzz, g_x_0_z_0_0_zzzz, g_x_0_z_0_x_xxxx, g_x_0_z_0_x_xxxy, g_x_0_z_0_x_xxxz, g_x_0_z_0_x_xxyy, g_x_0_z_0_x_xxyz, g_x_0_z_0_x_xxzz, g_x_0_z_0_x_xyyy, g_x_0_z_0_x_xyyz, g_x_0_z_0_x_xyzz, g_x_0_z_0_x_xzzz, g_x_0_z_0_x_yyyy, g_x_0_z_0_x_yyyz, g_x_0_z_0_x_yyzz, g_x_0_z_0_x_yzzz, g_x_0_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_x_xxxx[k] = -g_0_0_z_0_0_xxxx[k] - g_x_0_z_0_0_xxxx[k] * ab_x + g_x_0_z_0_0_xxxxx[k];

                g_x_0_z_0_x_xxxy[k] = -g_0_0_z_0_0_xxxy[k] - g_x_0_z_0_0_xxxy[k] * ab_x + g_x_0_z_0_0_xxxxy[k];

                g_x_0_z_0_x_xxxz[k] = -g_0_0_z_0_0_xxxz[k] - g_x_0_z_0_0_xxxz[k] * ab_x + g_x_0_z_0_0_xxxxz[k];

                g_x_0_z_0_x_xxyy[k] = -g_0_0_z_0_0_xxyy[k] - g_x_0_z_0_0_xxyy[k] * ab_x + g_x_0_z_0_0_xxxyy[k];

                g_x_0_z_0_x_xxyz[k] = -g_0_0_z_0_0_xxyz[k] - g_x_0_z_0_0_xxyz[k] * ab_x + g_x_0_z_0_0_xxxyz[k];

                g_x_0_z_0_x_xxzz[k] = -g_0_0_z_0_0_xxzz[k] - g_x_0_z_0_0_xxzz[k] * ab_x + g_x_0_z_0_0_xxxzz[k];

                g_x_0_z_0_x_xyyy[k] = -g_0_0_z_0_0_xyyy[k] - g_x_0_z_0_0_xyyy[k] * ab_x + g_x_0_z_0_0_xxyyy[k];

                g_x_0_z_0_x_xyyz[k] = -g_0_0_z_0_0_xyyz[k] - g_x_0_z_0_0_xyyz[k] * ab_x + g_x_0_z_0_0_xxyyz[k];

                g_x_0_z_0_x_xyzz[k] = -g_0_0_z_0_0_xyzz[k] - g_x_0_z_0_0_xyzz[k] * ab_x + g_x_0_z_0_0_xxyzz[k];

                g_x_0_z_0_x_xzzz[k] = -g_0_0_z_0_0_xzzz[k] - g_x_0_z_0_0_xzzz[k] * ab_x + g_x_0_z_0_0_xxzzz[k];

                g_x_0_z_0_x_yyyy[k] = -g_0_0_z_0_0_yyyy[k] - g_x_0_z_0_0_yyyy[k] * ab_x + g_x_0_z_0_0_xyyyy[k];

                g_x_0_z_0_x_yyyz[k] = -g_0_0_z_0_0_yyyz[k] - g_x_0_z_0_0_yyyz[k] * ab_x + g_x_0_z_0_0_xyyyz[k];

                g_x_0_z_0_x_yyzz[k] = -g_0_0_z_0_0_yyzz[k] - g_x_0_z_0_0_yyzz[k] * ab_x + g_x_0_z_0_0_xyyzz[k];

                g_x_0_z_0_x_yzzz[k] = -g_0_0_z_0_0_yzzz[k] - g_x_0_z_0_0_yzzz[k] * ab_x + g_x_0_z_0_0_xyzzz[k];

                g_x_0_z_0_x_zzzz[k] = -g_0_0_z_0_0_zzzz[k] - g_x_0_z_0_0_zzzz[k] * ab_x + g_x_0_z_0_0_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_z_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_z_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_z_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xxxx, g_x_0_z_0_0_xxxxy, g_x_0_z_0_0_xxxy, g_x_0_z_0_0_xxxyy, g_x_0_z_0_0_xxxyz, g_x_0_z_0_0_xxxz, g_x_0_z_0_0_xxyy, g_x_0_z_0_0_xxyyy, g_x_0_z_0_0_xxyyz, g_x_0_z_0_0_xxyz, g_x_0_z_0_0_xxyzz, g_x_0_z_0_0_xxzz, g_x_0_z_0_0_xyyy, g_x_0_z_0_0_xyyyy, g_x_0_z_0_0_xyyyz, g_x_0_z_0_0_xyyz, g_x_0_z_0_0_xyyzz, g_x_0_z_0_0_xyzz, g_x_0_z_0_0_xyzzz, g_x_0_z_0_0_xzzz, g_x_0_z_0_0_yyyy, g_x_0_z_0_0_yyyyy, g_x_0_z_0_0_yyyyz, g_x_0_z_0_0_yyyz, g_x_0_z_0_0_yyyzz, g_x_0_z_0_0_yyzz, g_x_0_z_0_0_yyzzz, g_x_0_z_0_0_yzzz, g_x_0_z_0_0_yzzzz, g_x_0_z_0_0_zzzz, g_x_0_z_0_y_xxxx, g_x_0_z_0_y_xxxy, g_x_0_z_0_y_xxxz, g_x_0_z_0_y_xxyy, g_x_0_z_0_y_xxyz, g_x_0_z_0_y_xxzz, g_x_0_z_0_y_xyyy, g_x_0_z_0_y_xyyz, g_x_0_z_0_y_xyzz, g_x_0_z_0_y_xzzz, g_x_0_z_0_y_yyyy, g_x_0_z_0_y_yyyz, g_x_0_z_0_y_yyzz, g_x_0_z_0_y_yzzz, g_x_0_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_y_xxxx[k] = -g_x_0_z_0_0_xxxx[k] * ab_y + g_x_0_z_0_0_xxxxy[k];

                g_x_0_z_0_y_xxxy[k] = -g_x_0_z_0_0_xxxy[k] * ab_y + g_x_0_z_0_0_xxxyy[k];

                g_x_0_z_0_y_xxxz[k] = -g_x_0_z_0_0_xxxz[k] * ab_y + g_x_0_z_0_0_xxxyz[k];

                g_x_0_z_0_y_xxyy[k] = -g_x_0_z_0_0_xxyy[k] * ab_y + g_x_0_z_0_0_xxyyy[k];

                g_x_0_z_0_y_xxyz[k] = -g_x_0_z_0_0_xxyz[k] * ab_y + g_x_0_z_0_0_xxyyz[k];

                g_x_0_z_0_y_xxzz[k] = -g_x_0_z_0_0_xxzz[k] * ab_y + g_x_0_z_0_0_xxyzz[k];

                g_x_0_z_0_y_xyyy[k] = -g_x_0_z_0_0_xyyy[k] * ab_y + g_x_0_z_0_0_xyyyy[k];

                g_x_0_z_0_y_xyyz[k] = -g_x_0_z_0_0_xyyz[k] * ab_y + g_x_0_z_0_0_xyyyz[k];

                g_x_0_z_0_y_xyzz[k] = -g_x_0_z_0_0_xyzz[k] * ab_y + g_x_0_z_0_0_xyyzz[k];

                g_x_0_z_0_y_xzzz[k] = -g_x_0_z_0_0_xzzz[k] * ab_y + g_x_0_z_0_0_xyzzz[k];

                g_x_0_z_0_y_yyyy[k] = -g_x_0_z_0_0_yyyy[k] * ab_y + g_x_0_z_0_0_yyyyy[k];

                g_x_0_z_0_y_yyyz[k] = -g_x_0_z_0_0_yyyz[k] * ab_y + g_x_0_z_0_0_yyyyz[k];

                g_x_0_z_0_y_yyzz[k] = -g_x_0_z_0_0_yyzz[k] * ab_y + g_x_0_z_0_0_yyyzz[k];

                g_x_0_z_0_y_yzzz[k] = -g_x_0_z_0_0_yzzz[k] * ab_y + g_x_0_z_0_0_yyzzz[k];

                g_x_0_z_0_y_zzzz[k] = -g_x_0_z_0_0_zzzz[k] * ab_y + g_x_0_z_0_0_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_z_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_z_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_z_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xxxx, g_x_0_z_0_0_xxxxz, g_x_0_z_0_0_xxxy, g_x_0_z_0_0_xxxyz, g_x_0_z_0_0_xxxz, g_x_0_z_0_0_xxxzz, g_x_0_z_0_0_xxyy, g_x_0_z_0_0_xxyyz, g_x_0_z_0_0_xxyz, g_x_0_z_0_0_xxyzz, g_x_0_z_0_0_xxzz, g_x_0_z_0_0_xxzzz, g_x_0_z_0_0_xyyy, g_x_0_z_0_0_xyyyz, g_x_0_z_0_0_xyyz, g_x_0_z_0_0_xyyzz, g_x_0_z_0_0_xyzz, g_x_0_z_0_0_xyzzz, g_x_0_z_0_0_xzzz, g_x_0_z_0_0_xzzzz, g_x_0_z_0_0_yyyy, g_x_0_z_0_0_yyyyz, g_x_0_z_0_0_yyyz, g_x_0_z_0_0_yyyzz, g_x_0_z_0_0_yyzz, g_x_0_z_0_0_yyzzz, g_x_0_z_0_0_yzzz, g_x_0_z_0_0_yzzzz, g_x_0_z_0_0_zzzz, g_x_0_z_0_0_zzzzz, g_x_0_z_0_z_xxxx, g_x_0_z_0_z_xxxy, g_x_0_z_0_z_xxxz, g_x_0_z_0_z_xxyy, g_x_0_z_0_z_xxyz, g_x_0_z_0_z_xxzz, g_x_0_z_0_z_xyyy, g_x_0_z_0_z_xyyz, g_x_0_z_0_z_xyzz, g_x_0_z_0_z_xzzz, g_x_0_z_0_z_yyyy, g_x_0_z_0_z_yyyz, g_x_0_z_0_z_yyzz, g_x_0_z_0_z_yzzz, g_x_0_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_z_xxxx[k] = -g_x_0_z_0_0_xxxx[k] * ab_z + g_x_0_z_0_0_xxxxz[k];

                g_x_0_z_0_z_xxxy[k] = -g_x_0_z_0_0_xxxy[k] * ab_z + g_x_0_z_0_0_xxxyz[k];

                g_x_0_z_0_z_xxxz[k] = -g_x_0_z_0_0_xxxz[k] * ab_z + g_x_0_z_0_0_xxxzz[k];

                g_x_0_z_0_z_xxyy[k] = -g_x_0_z_0_0_xxyy[k] * ab_z + g_x_0_z_0_0_xxyyz[k];

                g_x_0_z_0_z_xxyz[k] = -g_x_0_z_0_0_xxyz[k] * ab_z + g_x_0_z_0_0_xxyzz[k];

                g_x_0_z_0_z_xxzz[k] = -g_x_0_z_0_0_xxzz[k] * ab_z + g_x_0_z_0_0_xxzzz[k];

                g_x_0_z_0_z_xyyy[k] = -g_x_0_z_0_0_xyyy[k] * ab_z + g_x_0_z_0_0_xyyyz[k];

                g_x_0_z_0_z_xyyz[k] = -g_x_0_z_0_0_xyyz[k] * ab_z + g_x_0_z_0_0_xyyzz[k];

                g_x_0_z_0_z_xyzz[k] = -g_x_0_z_0_0_xyzz[k] * ab_z + g_x_0_z_0_0_xyzzz[k];

                g_x_0_z_0_z_xzzz[k] = -g_x_0_z_0_0_xzzz[k] * ab_z + g_x_0_z_0_0_xzzzz[k];

                g_x_0_z_0_z_yyyy[k] = -g_x_0_z_0_0_yyyy[k] * ab_z + g_x_0_z_0_0_yyyyz[k];

                g_x_0_z_0_z_yyyz[k] = -g_x_0_z_0_0_yyyz[k] * ab_z + g_x_0_z_0_0_yyyzz[k];

                g_x_0_z_0_z_yyzz[k] = -g_x_0_z_0_0_yyzz[k] * ab_z + g_x_0_z_0_0_yyzzz[k];

                g_x_0_z_0_z_yzzz[k] = -g_x_0_z_0_0_yzzz[k] * ab_z + g_x_0_z_0_0_yzzzz[k];

                g_x_0_z_0_z_zzzz[k] = -g_x_0_z_0_0_zzzz[k] * ab_z + g_x_0_z_0_0_zzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 135 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 136 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 137 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 138 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 139 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 140 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 141 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 142 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 143 * ccomps * dcomps);

            auto g_y_0_x_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 144 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 145 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 146 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 147 * ccomps * dcomps);

            auto g_y_0_x_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 148 * ccomps * dcomps);

            auto g_y_0_x_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xxxx, g_y_0_x_0_0_xxxxx, g_y_0_x_0_0_xxxxy, g_y_0_x_0_0_xxxxz, g_y_0_x_0_0_xxxy, g_y_0_x_0_0_xxxyy, g_y_0_x_0_0_xxxyz, g_y_0_x_0_0_xxxz, g_y_0_x_0_0_xxxzz, g_y_0_x_0_0_xxyy, g_y_0_x_0_0_xxyyy, g_y_0_x_0_0_xxyyz, g_y_0_x_0_0_xxyz, g_y_0_x_0_0_xxyzz, g_y_0_x_0_0_xxzz, g_y_0_x_0_0_xxzzz, g_y_0_x_0_0_xyyy, g_y_0_x_0_0_xyyyy, g_y_0_x_0_0_xyyyz, g_y_0_x_0_0_xyyz, g_y_0_x_0_0_xyyzz, g_y_0_x_0_0_xyzz, g_y_0_x_0_0_xyzzz, g_y_0_x_0_0_xzzz, g_y_0_x_0_0_xzzzz, g_y_0_x_0_0_yyyy, g_y_0_x_0_0_yyyz, g_y_0_x_0_0_yyzz, g_y_0_x_0_0_yzzz, g_y_0_x_0_0_zzzz, g_y_0_x_0_x_xxxx, g_y_0_x_0_x_xxxy, g_y_0_x_0_x_xxxz, g_y_0_x_0_x_xxyy, g_y_0_x_0_x_xxyz, g_y_0_x_0_x_xxzz, g_y_0_x_0_x_xyyy, g_y_0_x_0_x_xyyz, g_y_0_x_0_x_xyzz, g_y_0_x_0_x_xzzz, g_y_0_x_0_x_yyyy, g_y_0_x_0_x_yyyz, g_y_0_x_0_x_yyzz, g_y_0_x_0_x_yzzz, g_y_0_x_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_x_xxxx[k] = -g_y_0_x_0_0_xxxx[k] * ab_x + g_y_0_x_0_0_xxxxx[k];

                g_y_0_x_0_x_xxxy[k] = -g_y_0_x_0_0_xxxy[k] * ab_x + g_y_0_x_0_0_xxxxy[k];

                g_y_0_x_0_x_xxxz[k] = -g_y_0_x_0_0_xxxz[k] * ab_x + g_y_0_x_0_0_xxxxz[k];

                g_y_0_x_0_x_xxyy[k] = -g_y_0_x_0_0_xxyy[k] * ab_x + g_y_0_x_0_0_xxxyy[k];

                g_y_0_x_0_x_xxyz[k] = -g_y_0_x_0_0_xxyz[k] * ab_x + g_y_0_x_0_0_xxxyz[k];

                g_y_0_x_0_x_xxzz[k] = -g_y_0_x_0_0_xxzz[k] * ab_x + g_y_0_x_0_0_xxxzz[k];

                g_y_0_x_0_x_xyyy[k] = -g_y_0_x_0_0_xyyy[k] * ab_x + g_y_0_x_0_0_xxyyy[k];

                g_y_0_x_0_x_xyyz[k] = -g_y_0_x_0_0_xyyz[k] * ab_x + g_y_0_x_0_0_xxyyz[k];

                g_y_0_x_0_x_xyzz[k] = -g_y_0_x_0_0_xyzz[k] * ab_x + g_y_0_x_0_0_xxyzz[k];

                g_y_0_x_0_x_xzzz[k] = -g_y_0_x_0_0_xzzz[k] * ab_x + g_y_0_x_0_0_xxzzz[k];

                g_y_0_x_0_x_yyyy[k] = -g_y_0_x_0_0_yyyy[k] * ab_x + g_y_0_x_0_0_xyyyy[k];

                g_y_0_x_0_x_yyyz[k] = -g_y_0_x_0_0_yyyz[k] * ab_x + g_y_0_x_0_0_xyyyz[k];

                g_y_0_x_0_x_yyzz[k] = -g_y_0_x_0_0_yyzz[k] * ab_x + g_y_0_x_0_0_xyyzz[k];

                g_y_0_x_0_x_yzzz[k] = -g_y_0_x_0_0_yzzz[k] * ab_x + g_y_0_x_0_0_xyzzz[k];

                g_y_0_x_0_x_zzzz[k] = -g_y_0_x_0_0_zzzz[k] * ab_x + g_y_0_x_0_0_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 150 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 151 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 152 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 153 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 154 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 155 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 156 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 157 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 158 * ccomps * dcomps);

            auto g_y_0_x_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 159 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 160 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 161 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 162 * ccomps * dcomps);

            auto g_y_0_x_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 163 * ccomps * dcomps);

            auto g_y_0_x_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxx, g_0_0_x_0_0_xxxy, g_0_0_x_0_0_xxxz, g_0_0_x_0_0_xxyy, g_0_0_x_0_0_xxyz, g_0_0_x_0_0_xxzz, g_0_0_x_0_0_xyyy, g_0_0_x_0_0_xyyz, g_0_0_x_0_0_xyzz, g_0_0_x_0_0_xzzz, g_0_0_x_0_0_yyyy, g_0_0_x_0_0_yyyz, g_0_0_x_0_0_yyzz, g_0_0_x_0_0_yzzz, g_0_0_x_0_0_zzzz, g_y_0_x_0_0_xxxx, g_y_0_x_0_0_xxxxy, g_y_0_x_0_0_xxxy, g_y_0_x_0_0_xxxyy, g_y_0_x_0_0_xxxyz, g_y_0_x_0_0_xxxz, g_y_0_x_0_0_xxyy, g_y_0_x_0_0_xxyyy, g_y_0_x_0_0_xxyyz, g_y_0_x_0_0_xxyz, g_y_0_x_0_0_xxyzz, g_y_0_x_0_0_xxzz, g_y_0_x_0_0_xyyy, g_y_0_x_0_0_xyyyy, g_y_0_x_0_0_xyyyz, g_y_0_x_0_0_xyyz, g_y_0_x_0_0_xyyzz, g_y_0_x_0_0_xyzz, g_y_0_x_0_0_xyzzz, g_y_0_x_0_0_xzzz, g_y_0_x_0_0_yyyy, g_y_0_x_0_0_yyyyy, g_y_0_x_0_0_yyyyz, g_y_0_x_0_0_yyyz, g_y_0_x_0_0_yyyzz, g_y_0_x_0_0_yyzz, g_y_0_x_0_0_yyzzz, g_y_0_x_0_0_yzzz, g_y_0_x_0_0_yzzzz, g_y_0_x_0_0_zzzz, g_y_0_x_0_y_xxxx, g_y_0_x_0_y_xxxy, g_y_0_x_0_y_xxxz, g_y_0_x_0_y_xxyy, g_y_0_x_0_y_xxyz, g_y_0_x_0_y_xxzz, g_y_0_x_0_y_xyyy, g_y_0_x_0_y_xyyz, g_y_0_x_0_y_xyzz, g_y_0_x_0_y_xzzz, g_y_0_x_0_y_yyyy, g_y_0_x_0_y_yyyz, g_y_0_x_0_y_yyzz, g_y_0_x_0_y_yzzz, g_y_0_x_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_y_xxxx[k] = -g_0_0_x_0_0_xxxx[k] - g_y_0_x_0_0_xxxx[k] * ab_y + g_y_0_x_0_0_xxxxy[k];

                g_y_0_x_0_y_xxxy[k] = -g_0_0_x_0_0_xxxy[k] - g_y_0_x_0_0_xxxy[k] * ab_y + g_y_0_x_0_0_xxxyy[k];

                g_y_0_x_0_y_xxxz[k] = -g_0_0_x_0_0_xxxz[k] - g_y_0_x_0_0_xxxz[k] * ab_y + g_y_0_x_0_0_xxxyz[k];

                g_y_0_x_0_y_xxyy[k] = -g_0_0_x_0_0_xxyy[k] - g_y_0_x_0_0_xxyy[k] * ab_y + g_y_0_x_0_0_xxyyy[k];

                g_y_0_x_0_y_xxyz[k] = -g_0_0_x_0_0_xxyz[k] - g_y_0_x_0_0_xxyz[k] * ab_y + g_y_0_x_0_0_xxyyz[k];

                g_y_0_x_0_y_xxzz[k] = -g_0_0_x_0_0_xxzz[k] - g_y_0_x_0_0_xxzz[k] * ab_y + g_y_0_x_0_0_xxyzz[k];

                g_y_0_x_0_y_xyyy[k] = -g_0_0_x_0_0_xyyy[k] - g_y_0_x_0_0_xyyy[k] * ab_y + g_y_0_x_0_0_xyyyy[k];

                g_y_0_x_0_y_xyyz[k] = -g_0_0_x_0_0_xyyz[k] - g_y_0_x_0_0_xyyz[k] * ab_y + g_y_0_x_0_0_xyyyz[k];

                g_y_0_x_0_y_xyzz[k] = -g_0_0_x_0_0_xyzz[k] - g_y_0_x_0_0_xyzz[k] * ab_y + g_y_0_x_0_0_xyyzz[k];

                g_y_0_x_0_y_xzzz[k] = -g_0_0_x_0_0_xzzz[k] - g_y_0_x_0_0_xzzz[k] * ab_y + g_y_0_x_0_0_xyzzz[k];

                g_y_0_x_0_y_yyyy[k] = -g_0_0_x_0_0_yyyy[k] - g_y_0_x_0_0_yyyy[k] * ab_y + g_y_0_x_0_0_yyyyy[k];

                g_y_0_x_0_y_yyyz[k] = -g_0_0_x_0_0_yyyz[k] - g_y_0_x_0_0_yyyz[k] * ab_y + g_y_0_x_0_0_yyyyz[k];

                g_y_0_x_0_y_yyzz[k] = -g_0_0_x_0_0_yyzz[k] - g_y_0_x_0_0_yyzz[k] * ab_y + g_y_0_x_0_0_yyyzz[k];

                g_y_0_x_0_y_yzzz[k] = -g_0_0_x_0_0_yzzz[k] - g_y_0_x_0_0_yzzz[k] * ab_y + g_y_0_x_0_0_yyzzz[k];

                g_y_0_x_0_y_zzzz[k] = -g_0_0_x_0_0_zzzz[k] - g_y_0_x_0_0_zzzz[k] * ab_y + g_y_0_x_0_0_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 165 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 166 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 167 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 168 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 169 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 170 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 171 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 172 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 173 * ccomps * dcomps);

            auto g_y_0_x_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 174 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 175 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 176 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 177 * ccomps * dcomps);

            auto g_y_0_x_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 178 * ccomps * dcomps);

            auto g_y_0_x_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xxxx, g_y_0_x_0_0_xxxxz, g_y_0_x_0_0_xxxy, g_y_0_x_0_0_xxxyz, g_y_0_x_0_0_xxxz, g_y_0_x_0_0_xxxzz, g_y_0_x_0_0_xxyy, g_y_0_x_0_0_xxyyz, g_y_0_x_0_0_xxyz, g_y_0_x_0_0_xxyzz, g_y_0_x_0_0_xxzz, g_y_0_x_0_0_xxzzz, g_y_0_x_0_0_xyyy, g_y_0_x_0_0_xyyyz, g_y_0_x_0_0_xyyz, g_y_0_x_0_0_xyyzz, g_y_0_x_0_0_xyzz, g_y_0_x_0_0_xyzzz, g_y_0_x_0_0_xzzz, g_y_0_x_0_0_xzzzz, g_y_0_x_0_0_yyyy, g_y_0_x_0_0_yyyyz, g_y_0_x_0_0_yyyz, g_y_0_x_0_0_yyyzz, g_y_0_x_0_0_yyzz, g_y_0_x_0_0_yyzzz, g_y_0_x_0_0_yzzz, g_y_0_x_0_0_yzzzz, g_y_0_x_0_0_zzzz, g_y_0_x_0_0_zzzzz, g_y_0_x_0_z_xxxx, g_y_0_x_0_z_xxxy, g_y_0_x_0_z_xxxz, g_y_0_x_0_z_xxyy, g_y_0_x_0_z_xxyz, g_y_0_x_0_z_xxzz, g_y_0_x_0_z_xyyy, g_y_0_x_0_z_xyyz, g_y_0_x_0_z_xyzz, g_y_0_x_0_z_xzzz, g_y_0_x_0_z_yyyy, g_y_0_x_0_z_yyyz, g_y_0_x_0_z_yyzz, g_y_0_x_0_z_yzzz, g_y_0_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_z_xxxx[k] = -g_y_0_x_0_0_xxxx[k] * ab_z + g_y_0_x_0_0_xxxxz[k];

                g_y_0_x_0_z_xxxy[k] = -g_y_0_x_0_0_xxxy[k] * ab_z + g_y_0_x_0_0_xxxyz[k];

                g_y_0_x_0_z_xxxz[k] = -g_y_0_x_0_0_xxxz[k] * ab_z + g_y_0_x_0_0_xxxzz[k];

                g_y_0_x_0_z_xxyy[k] = -g_y_0_x_0_0_xxyy[k] * ab_z + g_y_0_x_0_0_xxyyz[k];

                g_y_0_x_0_z_xxyz[k] = -g_y_0_x_0_0_xxyz[k] * ab_z + g_y_0_x_0_0_xxyzz[k];

                g_y_0_x_0_z_xxzz[k] = -g_y_0_x_0_0_xxzz[k] * ab_z + g_y_0_x_0_0_xxzzz[k];

                g_y_0_x_0_z_xyyy[k] = -g_y_0_x_0_0_xyyy[k] * ab_z + g_y_0_x_0_0_xyyyz[k];

                g_y_0_x_0_z_xyyz[k] = -g_y_0_x_0_0_xyyz[k] * ab_z + g_y_0_x_0_0_xyyzz[k];

                g_y_0_x_0_z_xyzz[k] = -g_y_0_x_0_0_xyzz[k] * ab_z + g_y_0_x_0_0_xyzzz[k];

                g_y_0_x_0_z_xzzz[k] = -g_y_0_x_0_0_xzzz[k] * ab_z + g_y_0_x_0_0_xzzzz[k];

                g_y_0_x_0_z_yyyy[k] = -g_y_0_x_0_0_yyyy[k] * ab_z + g_y_0_x_0_0_yyyyz[k];

                g_y_0_x_0_z_yyyz[k] = -g_y_0_x_0_0_yyyz[k] * ab_z + g_y_0_x_0_0_yyyzz[k];

                g_y_0_x_0_z_yyzz[k] = -g_y_0_x_0_0_yyzz[k] * ab_z + g_y_0_x_0_0_yyzzz[k];

                g_y_0_x_0_z_yzzz[k] = -g_y_0_x_0_0_yzzz[k] * ab_z + g_y_0_x_0_0_yzzzz[k];

                g_y_0_x_0_z_zzzz[k] = -g_y_0_x_0_0_zzzz[k] * ab_z + g_y_0_x_0_0_zzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 180 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 181 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 182 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 183 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 184 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 185 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 186 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 187 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 188 * ccomps * dcomps);

            auto g_y_0_y_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 189 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 190 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 191 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 192 * ccomps * dcomps);

            auto g_y_0_y_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 193 * ccomps * dcomps);

            auto g_y_0_y_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xxxx, g_y_0_y_0_0_xxxxx, g_y_0_y_0_0_xxxxy, g_y_0_y_0_0_xxxxz, g_y_0_y_0_0_xxxy, g_y_0_y_0_0_xxxyy, g_y_0_y_0_0_xxxyz, g_y_0_y_0_0_xxxz, g_y_0_y_0_0_xxxzz, g_y_0_y_0_0_xxyy, g_y_0_y_0_0_xxyyy, g_y_0_y_0_0_xxyyz, g_y_0_y_0_0_xxyz, g_y_0_y_0_0_xxyzz, g_y_0_y_0_0_xxzz, g_y_0_y_0_0_xxzzz, g_y_0_y_0_0_xyyy, g_y_0_y_0_0_xyyyy, g_y_0_y_0_0_xyyyz, g_y_0_y_0_0_xyyz, g_y_0_y_0_0_xyyzz, g_y_0_y_0_0_xyzz, g_y_0_y_0_0_xyzzz, g_y_0_y_0_0_xzzz, g_y_0_y_0_0_xzzzz, g_y_0_y_0_0_yyyy, g_y_0_y_0_0_yyyz, g_y_0_y_0_0_yyzz, g_y_0_y_0_0_yzzz, g_y_0_y_0_0_zzzz, g_y_0_y_0_x_xxxx, g_y_0_y_0_x_xxxy, g_y_0_y_0_x_xxxz, g_y_0_y_0_x_xxyy, g_y_0_y_0_x_xxyz, g_y_0_y_0_x_xxzz, g_y_0_y_0_x_xyyy, g_y_0_y_0_x_xyyz, g_y_0_y_0_x_xyzz, g_y_0_y_0_x_xzzz, g_y_0_y_0_x_yyyy, g_y_0_y_0_x_yyyz, g_y_0_y_0_x_yyzz, g_y_0_y_0_x_yzzz, g_y_0_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_x_xxxx[k] = -g_y_0_y_0_0_xxxx[k] * ab_x + g_y_0_y_0_0_xxxxx[k];

                g_y_0_y_0_x_xxxy[k] = -g_y_0_y_0_0_xxxy[k] * ab_x + g_y_0_y_0_0_xxxxy[k];

                g_y_0_y_0_x_xxxz[k] = -g_y_0_y_0_0_xxxz[k] * ab_x + g_y_0_y_0_0_xxxxz[k];

                g_y_0_y_0_x_xxyy[k] = -g_y_0_y_0_0_xxyy[k] * ab_x + g_y_0_y_0_0_xxxyy[k];

                g_y_0_y_0_x_xxyz[k] = -g_y_0_y_0_0_xxyz[k] * ab_x + g_y_0_y_0_0_xxxyz[k];

                g_y_0_y_0_x_xxzz[k] = -g_y_0_y_0_0_xxzz[k] * ab_x + g_y_0_y_0_0_xxxzz[k];

                g_y_0_y_0_x_xyyy[k] = -g_y_0_y_0_0_xyyy[k] * ab_x + g_y_0_y_0_0_xxyyy[k];

                g_y_0_y_0_x_xyyz[k] = -g_y_0_y_0_0_xyyz[k] * ab_x + g_y_0_y_0_0_xxyyz[k];

                g_y_0_y_0_x_xyzz[k] = -g_y_0_y_0_0_xyzz[k] * ab_x + g_y_0_y_0_0_xxyzz[k];

                g_y_0_y_0_x_xzzz[k] = -g_y_0_y_0_0_xzzz[k] * ab_x + g_y_0_y_0_0_xxzzz[k];

                g_y_0_y_0_x_yyyy[k] = -g_y_0_y_0_0_yyyy[k] * ab_x + g_y_0_y_0_0_xyyyy[k];

                g_y_0_y_0_x_yyyz[k] = -g_y_0_y_0_0_yyyz[k] * ab_x + g_y_0_y_0_0_xyyyz[k];

                g_y_0_y_0_x_yyzz[k] = -g_y_0_y_0_0_yyzz[k] * ab_x + g_y_0_y_0_0_xyyzz[k];

                g_y_0_y_0_x_yzzz[k] = -g_y_0_y_0_0_yzzz[k] * ab_x + g_y_0_y_0_0_xyzzz[k];

                g_y_0_y_0_x_zzzz[k] = -g_y_0_y_0_0_zzzz[k] * ab_x + g_y_0_y_0_0_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 195 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 196 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 197 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 198 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 199 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 200 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 201 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 202 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 203 * ccomps * dcomps);

            auto g_y_0_y_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 204 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 205 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 206 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 207 * ccomps * dcomps);

            auto g_y_0_y_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 208 * ccomps * dcomps);

            auto g_y_0_y_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxx, g_0_0_y_0_0_xxxy, g_0_0_y_0_0_xxxz, g_0_0_y_0_0_xxyy, g_0_0_y_0_0_xxyz, g_0_0_y_0_0_xxzz, g_0_0_y_0_0_xyyy, g_0_0_y_0_0_xyyz, g_0_0_y_0_0_xyzz, g_0_0_y_0_0_xzzz, g_0_0_y_0_0_yyyy, g_0_0_y_0_0_yyyz, g_0_0_y_0_0_yyzz, g_0_0_y_0_0_yzzz, g_0_0_y_0_0_zzzz, g_y_0_y_0_0_xxxx, g_y_0_y_0_0_xxxxy, g_y_0_y_0_0_xxxy, g_y_0_y_0_0_xxxyy, g_y_0_y_0_0_xxxyz, g_y_0_y_0_0_xxxz, g_y_0_y_0_0_xxyy, g_y_0_y_0_0_xxyyy, g_y_0_y_0_0_xxyyz, g_y_0_y_0_0_xxyz, g_y_0_y_0_0_xxyzz, g_y_0_y_0_0_xxzz, g_y_0_y_0_0_xyyy, g_y_0_y_0_0_xyyyy, g_y_0_y_0_0_xyyyz, g_y_0_y_0_0_xyyz, g_y_0_y_0_0_xyyzz, g_y_0_y_0_0_xyzz, g_y_0_y_0_0_xyzzz, g_y_0_y_0_0_xzzz, g_y_0_y_0_0_yyyy, g_y_0_y_0_0_yyyyy, g_y_0_y_0_0_yyyyz, g_y_0_y_0_0_yyyz, g_y_0_y_0_0_yyyzz, g_y_0_y_0_0_yyzz, g_y_0_y_0_0_yyzzz, g_y_0_y_0_0_yzzz, g_y_0_y_0_0_yzzzz, g_y_0_y_0_0_zzzz, g_y_0_y_0_y_xxxx, g_y_0_y_0_y_xxxy, g_y_0_y_0_y_xxxz, g_y_0_y_0_y_xxyy, g_y_0_y_0_y_xxyz, g_y_0_y_0_y_xxzz, g_y_0_y_0_y_xyyy, g_y_0_y_0_y_xyyz, g_y_0_y_0_y_xyzz, g_y_0_y_0_y_xzzz, g_y_0_y_0_y_yyyy, g_y_0_y_0_y_yyyz, g_y_0_y_0_y_yyzz, g_y_0_y_0_y_yzzz, g_y_0_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_y_xxxx[k] = -g_0_0_y_0_0_xxxx[k] - g_y_0_y_0_0_xxxx[k] * ab_y + g_y_0_y_0_0_xxxxy[k];

                g_y_0_y_0_y_xxxy[k] = -g_0_0_y_0_0_xxxy[k] - g_y_0_y_0_0_xxxy[k] * ab_y + g_y_0_y_0_0_xxxyy[k];

                g_y_0_y_0_y_xxxz[k] = -g_0_0_y_0_0_xxxz[k] - g_y_0_y_0_0_xxxz[k] * ab_y + g_y_0_y_0_0_xxxyz[k];

                g_y_0_y_0_y_xxyy[k] = -g_0_0_y_0_0_xxyy[k] - g_y_0_y_0_0_xxyy[k] * ab_y + g_y_0_y_0_0_xxyyy[k];

                g_y_0_y_0_y_xxyz[k] = -g_0_0_y_0_0_xxyz[k] - g_y_0_y_0_0_xxyz[k] * ab_y + g_y_0_y_0_0_xxyyz[k];

                g_y_0_y_0_y_xxzz[k] = -g_0_0_y_0_0_xxzz[k] - g_y_0_y_0_0_xxzz[k] * ab_y + g_y_0_y_0_0_xxyzz[k];

                g_y_0_y_0_y_xyyy[k] = -g_0_0_y_0_0_xyyy[k] - g_y_0_y_0_0_xyyy[k] * ab_y + g_y_0_y_0_0_xyyyy[k];

                g_y_0_y_0_y_xyyz[k] = -g_0_0_y_0_0_xyyz[k] - g_y_0_y_0_0_xyyz[k] * ab_y + g_y_0_y_0_0_xyyyz[k];

                g_y_0_y_0_y_xyzz[k] = -g_0_0_y_0_0_xyzz[k] - g_y_0_y_0_0_xyzz[k] * ab_y + g_y_0_y_0_0_xyyzz[k];

                g_y_0_y_0_y_xzzz[k] = -g_0_0_y_0_0_xzzz[k] - g_y_0_y_0_0_xzzz[k] * ab_y + g_y_0_y_0_0_xyzzz[k];

                g_y_0_y_0_y_yyyy[k] = -g_0_0_y_0_0_yyyy[k] - g_y_0_y_0_0_yyyy[k] * ab_y + g_y_0_y_0_0_yyyyy[k];

                g_y_0_y_0_y_yyyz[k] = -g_0_0_y_0_0_yyyz[k] - g_y_0_y_0_0_yyyz[k] * ab_y + g_y_0_y_0_0_yyyyz[k];

                g_y_0_y_0_y_yyzz[k] = -g_0_0_y_0_0_yyzz[k] - g_y_0_y_0_0_yyzz[k] * ab_y + g_y_0_y_0_0_yyyzz[k];

                g_y_0_y_0_y_yzzz[k] = -g_0_0_y_0_0_yzzz[k] - g_y_0_y_0_0_yzzz[k] * ab_y + g_y_0_y_0_0_yyzzz[k];

                g_y_0_y_0_y_zzzz[k] = -g_0_0_y_0_0_zzzz[k] - g_y_0_y_0_0_zzzz[k] * ab_y + g_y_0_y_0_0_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 210 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 211 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 212 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 213 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 214 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 215 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 216 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 217 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 218 * ccomps * dcomps);

            auto g_y_0_y_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 219 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 220 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 221 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 222 * ccomps * dcomps);

            auto g_y_0_y_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 223 * ccomps * dcomps);

            auto g_y_0_y_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xxxx, g_y_0_y_0_0_xxxxz, g_y_0_y_0_0_xxxy, g_y_0_y_0_0_xxxyz, g_y_0_y_0_0_xxxz, g_y_0_y_0_0_xxxzz, g_y_0_y_0_0_xxyy, g_y_0_y_0_0_xxyyz, g_y_0_y_0_0_xxyz, g_y_0_y_0_0_xxyzz, g_y_0_y_0_0_xxzz, g_y_0_y_0_0_xxzzz, g_y_0_y_0_0_xyyy, g_y_0_y_0_0_xyyyz, g_y_0_y_0_0_xyyz, g_y_0_y_0_0_xyyzz, g_y_0_y_0_0_xyzz, g_y_0_y_0_0_xyzzz, g_y_0_y_0_0_xzzz, g_y_0_y_0_0_xzzzz, g_y_0_y_0_0_yyyy, g_y_0_y_0_0_yyyyz, g_y_0_y_0_0_yyyz, g_y_0_y_0_0_yyyzz, g_y_0_y_0_0_yyzz, g_y_0_y_0_0_yyzzz, g_y_0_y_0_0_yzzz, g_y_0_y_0_0_yzzzz, g_y_0_y_0_0_zzzz, g_y_0_y_0_0_zzzzz, g_y_0_y_0_z_xxxx, g_y_0_y_0_z_xxxy, g_y_0_y_0_z_xxxz, g_y_0_y_0_z_xxyy, g_y_0_y_0_z_xxyz, g_y_0_y_0_z_xxzz, g_y_0_y_0_z_xyyy, g_y_0_y_0_z_xyyz, g_y_0_y_0_z_xyzz, g_y_0_y_0_z_xzzz, g_y_0_y_0_z_yyyy, g_y_0_y_0_z_yyyz, g_y_0_y_0_z_yyzz, g_y_0_y_0_z_yzzz, g_y_0_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_z_xxxx[k] = -g_y_0_y_0_0_xxxx[k] * ab_z + g_y_0_y_0_0_xxxxz[k];

                g_y_0_y_0_z_xxxy[k] = -g_y_0_y_0_0_xxxy[k] * ab_z + g_y_0_y_0_0_xxxyz[k];

                g_y_0_y_0_z_xxxz[k] = -g_y_0_y_0_0_xxxz[k] * ab_z + g_y_0_y_0_0_xxxzz[k];

                g_y_0_y_0_z_xxyy[k] = -g_y_0_y_0_0_xxyy[k] * ab_z + g_y_0_y_0_0_xxyyz[k];

                g_y_0_y_0_z_xxyz[k] = -g_y_0_y_0_0_xxyz[k] * ab_z + g_y_0_y_0_0_xxyzz[k];

                g_y_0_y_0_z_xxzz[k] = -g_y_0_y_0_0_xxzz[k] * ab_z + g_y_0_y_0_0_xxzzz[k];

                g_y_0_y_0_z_xyyy[k] = -g_y_0_y_0_0_xyyy[k] * ab_z + g_y_0_y_0_0_xyyyz[k];

                g_y_0_y_0_z_xyyz[k] = -g_y_0_y_0_0_xyyz[k] * ab_z + g_y_0_y_0_0_xyyzz[k];

                g_y_0_y_0_z_xyzz[k] = -g_y_0_y_0_0_xyzz[k] * ab_z + g_y_0_y_0_0_xyzzz[k];

                g_y_0_y_0_z_xzzz[k] = -g_y_0_y_0_0_xzzz[k] * ab_z + g_y_0_y_0_0_xzzzz[k];

                g_y_0_y_0_z_yyyy[k] = -g_y_0_y_0_0_yyyy[k] * ab_z + g_y_0_y_0_0_yyyyz[k];

                g_y_0_y_0_z_yyyz[k] = -g_y_0_y_0_0_yyyz[k] * ab_z + g_y_0_y_0_0_yyyzz[k];

                g_y_0_y_0_z_yyzz[k] = -g_y_0_y_0_0_yyzz[k] * ab_z + g_y_0_y_0_0_yyzzz[k];

                g_y_0_y_0_z_yzzz[k] = -g_y_0_y_0_0_yzzz[k] * ab_z + g_y_0_y_0_0_yzzzz[k];

                g_y_0_y_0_z_zzzz[k] = -g_y_0_y_0_0_zzzz[k] * ab_z + g_y_0_y_0_0_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 225 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 226 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 227 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 228 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 229 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 230 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 231 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 232 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 233 * ccomps * dcomps);

            auto g_y_0_z_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 234 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 235 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 236 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 237 * ccomps * dcomps);

            auto g_y_0_z_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 238 * ccomps * dcomps);

            auto g_y_0_z_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xxxx, g_y_0_z_0_0_xxxxx, g_y_0_z_0_0_xxxxy, g_y_0_z_0_0_xxxxz, g_y_0_z_0_0_xxxy, g_y_0_z_0_0_xxxyy, g_y_0_z_0_0_xxxyz, g_y_0_z_0_0_xxxz, g_y_0_z_0_0_xxxzz, g_y_0_z_0_0_xxyy, g_y_0_z_0_0_xxyyy, g_y_0_z_0_0_xxyyz, g_y_0_z_0_0_xxyz, g_y_0_z_0_0_xxyzz, g_y_0_z_0_0_xxzz, g_y_0_z_0_0_xxzzz, g_y_0_z_0_0_xyyy, g_y_0_z_0_0_xyyyy, g_y_0_z_0_0_xyyyz, g_y_0_z_0_0_xyyz, g_y_0_z_0_0_xyyzz, g_y_0_z_0_0_xyzz, g_y_0_z_0_0_xyzzz, g_y_0_z_0_0_xzzz, g_y_0_z_0_0_xzzzz, g_y_0_z_0_0_yyyy, g_y_0_z_0_0_yyyz, g_y_0_z_0_0_yyzz, g_y_0_z_0_0_yzzz, g_y_0_z_0_0_zzzz, g_y_0_z_0_x_xxxx, g_y_0_z_0_x_xxxy, g_y_0_z_0_x_xxxz, g_y_0_z_0_x_xxyy, g_y_0_z_0_x_xxyz, g_y_0_z_0_x_xxzz, g_y_0_z_0_x_xyyy, g_y_0_z_0_x_xyyz, g_y_0_z_0_x_xyzz, g_y_0_z_0_x_xzzz, g_y_0_z_0_x_yyyy, g_y_0_z_0_x_yyyz, g_y_0_z_0_x_yyzz, g_y_0_z_0_x_yzzz, g_y_0_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_x_xxxx[k] = -g_y_0_z_0_0_xxxx[k] * ab_x + g_y_0_z_0_0_xxxxx[k];

                g_y_0_z_0_x_xxxy[k] = -g_y_0_z_0_0_xxxy[k] * ab_x + g_y_0_z_0_0_xxxxy[k];

                g_y_0_z_0_x_xxxz[k] = -g_y_0_z_0_0_xxxz[k] * ab_x + g_y_0_z_0_0_xxxxz[k];

                g_y_0_z_0_x_xxyy[k] = -g_y_0_z_0_0_xxyy[k] * ab_x + g_y_0_z_0_0_xxxyy[k];

                g_y_0_z_0_x_xxyz[k] = -g_y_0_z_0_0_xxyz[k] * ab_x + g_y_0_z_0_0_xxxyz[k];

                g_y_0_z_0_x_xxzz[k] = -g_y_0_z_0_0_xxzz[k] * ab_x + g_y_0_z_0_0_xxxzz[k];

                g_y_0_z_0_x_xyyy[k] = -g_y_0_z_0_0_xyyy[k] * ab_x + g_y_0_z_0_0_xxyyy[k];

                g_y_0_z_0_x_xyyz[k] = -g_y_0_z_0_0_xyyz[k] * ab_x + g_y_0_z_0_0_xxyyz[k];

                g_y_0_z_0_x_xyzz[k] = -g_y_0_z_0_0_xyzz[k] * ab_x + g_y_0_z_0_0_xxyzz[k];

                g_y_0_z_0_x_xzzz[k] = -g_y_0_z_0_0_xzzz[k] * ab_x + g_y_0_z_0_0_xxzzz[k];

                g_y_0_z_0_x_yyyy[k] = -g_y_0_z_0_0_yyyy[k] * ab_x + g_y_0_z_0_0_xyyyy[k];

                g_y_0_z_0_x_yyyz[k] = -g_y_0_z_0_0_yyyz[k] * ab_x + g_y_0_z_0_0_xyyyz[k];

                g_y_0_z_0_x_yyzz[k] = -g_y_0_z_0_0_yyzz[k] * ab_x + g_y_0_z_0_0_xyyzz[k];

                g_y_0_z_0_x_yzzz[k] = -g_y_0_z_0_0_yzzz[k] * ab_x + g_y_0_z_0_0_xyzzz[k];

                g_y_0_z_0_x_zzzz[k] = -g_y_0_z_0_0_zzzz[k] * ab_x + g_y_0_z_0_0_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 240 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 241 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 242 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 243 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 244 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 245 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 246 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 247 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 248 * ccomps * dcomps);

            auto g_y_0_z_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 249 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 250 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 251 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_z_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_z_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxx, g_0_0_z_0_0_xxxy, g_0_0_z_0_0_xxxz, g_0_0_z_0_0_xxyy, g_0_0_z_0_0_xxyz, g_0_0_z_0_0_xxzz, g_0_0_z_0_0_xyyy, g_0_0_z_0_0_xyyz, g_0_0_z_0_0_xyzz, g_0_0_z_0_0_xzzz, g_0_0_z_0_0_yyyy, g_0_0_z_0_0_yyyz, g_0_0_z_0_0_yyzz, g_0_0_z_0_0_yzzz, g_0_0_z_0_0_zzzz, g_y_0_z_0_0_xxxx, g_y_0_z_0_0_xxxxy, g_y_0_z_0_0_xxxy, g_y_0_z_0_0_xxxyy, g_y_0_z_0_0_xxxyz, g_y_0_z_0_0_xxxz, g_y_0_z_0_0_xxyy, g_y_0_z_0_0_xxyyy, g_y_0_z_0_0_xxyyz, g_y_0_z_0_0_xxyz, g_y_0_z_0_0_xxyzz, g_y_0_z_0_0_xxzz, g_y_0_z_0_0_xyyy, g_y_0_z_0_0_xyyyy, g_y_0_z_0_0_xyyyz, g_y_0_z_0_0_xyyz, g_y_0_z_0_0_xyyzz, g_y_0_z_0_0_xyzz, g_y_0_z_0_0_xyzzz, g_y_0_z_0_0_xzzz, g_y_0_z_0_0_yyyy, g_y_0_z_0_0_yyyyy, g_y_0_z_0_0_yyyyz, g_y_0_z_0_0_yyyz, g_y_0_z_0_0_yyyzz, g_y_0_z_0_0_yyzz, g_y_0_z_0_0_yyzzz, g_y_0_z_0_0_yzzz, g_y_0_z_0_0_yzzzz, g_y_0_z_0_0_zzzz, g_y_0_z_0_y_xxxx, g_y_0_z_0_y_xxxy, g_y_0_z_0_y_xxxz, g_y_0_z_0_y_xxyy, g_y_0_z_0_y_xxyz, g_y_0_z_0_y_xxzz, g_y_0_z_0_y_xyyy, g_y_0_z_0_y_xyyz, g_y_0_z_0_y_xyzz, g_y_0_z_0_y_xzzz, g_y_0_z_0_y_yyyy, g_y_0_z_0_y_yyyz, g_y_0_z_0_y_yyzz, g_y_0_z_0_y_yzzz, g_y_0_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_y_xxxx[k] = -g_0_0_z_0_0_xxxx[k] - g_y_0_z_0_0_xxxx[k] * ab_y + g_y_0_z_0_0_xxxxy[k];

                g_y_0_z_0_y_xxxy[k] = -g_0_0_z_0_0_xxxy[k] - g_y_0_z_0_0_xxxy[k] * ab_y + g_y_0_z_0_0_xxxyy[k];

                g_y_0_z_0_y_xxxz[k] = -g_0_0_z_0_0_xxxz[k] - g_y_0_z_0_0_xxxz[k] * ab_y + g_y_0_z_0_0_xxxyz[k];

                g_y_0_z_0_y_xxyy[k] = -g_0_0_z_0_0_xxyy[k] - g_y_0_z_0_0_xxyy[k] * ab_y + g_y_0_z_0_0_xxyyy[k];

                g_y_0_z_0_y_xxyz[k] = -g_0_0_z_0_0_xxyz[k] - g_y_0_z_0_0_xxyz[k] * ab_y + g_y_0_z_0_0_xxyyz[k];

                g_y_0_z_0_y_xxzz[k] = -g_0_0_z_0_0_xxzz[k] - g_y_0_z_0_0_xxzz[k] * ab_y + g_y_0_z_0_0_xxyzz[k];

                g_y_0_z_0_y_xyyy[k] = -g_0_0_z_0_0_xyyy[k] - g_y_0_z_0_0_xyyy[k] * ab_y + g_y_0_z_0_0_xyyyy[k];

                g_y_0_z_0_y_xyyz[k] = -g_0_0_z_0_0_xyyz[k] - g_y_0_z_0_0_xyyz[k] * ab_y + g_y_0_z_0_0_xyyyz[k];

                g_y_0_z_0_y_xyzz[k] = -g_0_0_z_0_0_xyzz[k] - g_y_0_z_0_0_xyzz[k] * ab_y + g_y_0_z_0_0_xyyzz[k];

                g_y_0_z_0_y_xzzz[k] = -g_0_0_z_0_0_xzzz[k] - g_y_0_z_0_0_xzzz[k] * ab_y + g_y_0_z_0_0_xyzzz[k];

                g_y_0_z_0_y_yyyy[k] = -g_0_0_z_0_0_yyyy[k] - g_y_0_z_0_0_yyyy[k] * ab_y + g_y_0_z_0_0_yyyyy[k];

                g_y_0_z_0_y_yyyz[k] = -g_0_0_z_0_0_yyyz[k] - g_y_0_z_0_0_yyyz[k] * ab_y + g_y_0_z_0_0_yyyyz[k];

                g_y_0_z_0_y_yyzz[k] = -g_0_0_z_0_0_yyzz[k] - g_y_0_z_0_0_yyzz[k] * ab_y + g_y_0_z_0_0_yyyzz[k];

                g_y_0_z_0_y_yzzz[k] = -g_0_0_z_0_0_yzzz[k] - g_y_0_z_0_0_yzzz[k] * ab_y + g_y_0_z_0_0_yyzzz[k];

                g_y_0_z_0_y_zzzz[k] = -g_0_0_z_0_0_zzzz[k] - g_y_0_z_0_0_zzzz[k] * ab_y + g_y_0_z_0_0_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 257 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 260 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 263 * ccomps * dcomps);

            auto g_y_0_z_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 266 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_z_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_z_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xxxx, g_y_0_z_0_0_xxxxz, g_y_0_z_0_0_xxxy, g_y_0_z_0_0_xxxyz, g_y_0_z_0_0_xxxz, g_y_0_z_0_0_xxxzz, g_y_0_z_0_0_xxyy, g_y_0_z_0_0_xxyyz, g_y_0_z_0_0_xxyz, g_y_0_z_0_0_xxyzz, g_y_0_z_0_0_xxzz, g_y_0_z_0_0_xxzzz, g_y_0_z_0_0_xyyy, g_y_0_z_0_0_xyyyz, g_y_0_z_0_0_xyyz, g_y_0_z_0_0_xyyzz, g_y_0_z_0_0_xyzz, g_y_0_z_0_0_xyzzz, g_y_0_z_0_0_xzzz, g_y_0_z_0_0_xzzzz, g_y_0_z_0_0_yyyy, g_y_0_z_0_0_yyyyz, g_y_0_z_0_0_yyyz, g_y_0_z_0_0_yyyzz, g_y_0_z_0_0_yyzz, g_y_0_z_0_0_yyzzz, g_y_0_z_0_0_yzzz, g_y_0_z_0_0_yzzzz, g_y_0_z_0_0_zzzz, g_y_0_z_0_0_zzzzz, g_y_0_z_0_z_xxxx, g_y_0_z_0_z_xxxy, g_y_0_z_0_z_xxxz, g_y_0_z_0_z_xxyy, g_y_0_z_0_z_xxyz, g_y_0_z_0_z_xxzz, g_y_0_z_0_z_xyyy, g_y_0_z_0_z_xyyz, g_y_0_z_0_z_xyzz, g_y_0_z_0_z_xzzz, g_y_0_z_0_z_yyyy, g_y_0_z_0_z_yyyz, g_y_0_z_0_z_yyzz, g_y_0_z_0_z_yzzz, g_y_0_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_z_xxxx[k] = -g_y_0_z_0_0_xxxx[k] * ab_z + g_y_0_z_0_0_xxxxz[k];

                g_y_0_z_0_z_xxxy[k] = -g_y_0_z_0_0_xxxy[k] * ab_z + g_y_0_z_0_0_xxxyz[k];

                g_y_0_z_0_z_xxxz[k] = -g_y_0_z_0_0_xxxz[k] * ab_z + g_y_0_z_0_0_xxxzz[k];

                g_y_0_z_0_z_xxyy[k] = -g_y_0_z_0_0_xxyy[k] * ab_z + g_y_0_z_0_0_xxyyz[k];

                g_y_0_z_0_z_xxyz[k] = -g_y_0_z_0_0_xxyz[k] * ab_z + g_y_0_z_0_0_xxyzz[k];

                g_y_0_z_0_z_xxzz[k] = -g_y_0_z_0_0_xxzz[k] * ab_z + g_y_0_z_0_0_xxzzz[k];

                g_y_0_z_0_z_xyyy[k] = -g_y_0_z_0_0_xyyy[k] * ab_z + g_y_0_z_0_0_xyyyz[k];

                g_y_0_z_0_z_xyyz[k] = -g_y_0_z_0_0_xyyz[k] * ab_z + g_y_0_z_0_0_xyyzz[k];

                g_y_0_z_0_z_xyzz[k] = -g_y_0_z_0_0_xyzz[k] * ab_z + g_y_0_z_0_0_xyzzz[k];

                g_y_0_z_0_z_xzzz[k] = -g_y_0_z_0_0_xzzz[k] * ab_z + g_y_0_z_0_0_xzzzz[k];

                g_y_0_z_0_z_yyyy[k] = -g_y_0_z_0_0_yyyy[k] * ab_z + g_y_0_z_0_0_yyyyz[k];

                g_y_0_z_0_z_yyyz[k] = -g_y_0_z_0_0_yyyz[k] * ab_z + g_y_0_z_0_0_yyyzz[k];

                g_y_0_z_0_z_yyzz[k] = -g_y_0_z_0_0_yyzz[k] * ab_z + g_y_0_z_0_0_yyzzz[k];

                g_y_0_z_0_z_yzzz[k] = -g_y_0_z_0_0_yzzz[k] * ab_z + g_y_0_z_0_0_yzzzz[k];

                g_y_0_z_0_z_zzzz[k] = -g_y_0_z_0_0_zzzz[k] * ab_z + g_y_0_z_0_0_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 270 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 271 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 272 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 273 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 274 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 275 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 276 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 277 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 278 * ccomps * dcomps);

            auto g_z_0_x_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 279 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 280 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 281 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 282 * ccomps * dcomps);

            auto g_z_0_x_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 283 * ccomps * dcomps);

            auto g_z_0_x_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xxxx, g_z_0_x_0_0_xxxxx, g_z_0_x_0_0_xxxxy, g_z_0_x_0_0_xxxxz, g_z_0_x_0_0_xxxy, g_z_0_x_0_0_xxxyy, g_z_0_x_0_0_xxxyz, g_z_0_x_0_0_xxxz, g_z_0_x_0_0_xxxzz, g_z_0_x_0_0_xxyy, g_z_0_x_0_0_xxyyy, g_z_0_x_0_0_xxyyz, g_z_0_x_0_0_xxyz, g_z_0_x_0_0_xxyzz, g_z_0_x_0_0_xxzz, g_z_0_x_0_0_xxzzz, g_z_0_x_0_0_xyyy, g_z_0_x_0_0_xyyyy, g_z_0_x_0_0_xyyyz, g_z_0_x_0_0_xyyz, g_z_0_x_0_0_xyyzz, g_z_0_x_0_0_xyzz, g_z_0_x_0_0_xyzzz, g_z_0_x_0_0_xzzz, g_z_0_x_0_0_xzzzz, g_z_0_x_0_0_yyyy, g_z_0_x_0_0_yyyz, g_z_0_x_0_0_yyzz, g_z_0_x_0_0_yzzz, g_z_0_x_0_0_zzzz, g_z_0_x_0_x_xxxx, g_z_0_x_0_x_xxxy, g_z_0_x_0_x_xxxz, g_z_0_x_0_x_xxyy, g_z_0_x_0_x_xxyz, g_z_0_x_0_x_xxzz, g_z_0_x_0_x_xyyy, g_z_0_x_0_x_xyyz, g_z_0_x_0_x_xyzz, g_z_0_x_0_x_xzzz, g_z_0_x_0_x_yyyy, g_z_0_x_0_x_yyyz, g_z_0_x_0_x_yyzz, g_z_0_x_0_x_yzzz, g_z_0_x_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_x_xxxx[k] = -g_z_0_x_0_0_xxxx[k] * ab_x + g_z_0_x_0_0_xxxxx[k];

                g_z_0_x_0_x_xxxy[k] = -g_z_0_x_0_0_xxxy[k] * ab_x + g_z_0_x_0_0_xxxxy[k];

                g_z_0_x_0_x_xxxz[k] = -g_z_0_x_0_0_xxxz[k] * ab_x + g_z_0_x_0_0_xxxxz[k];

                g_z_0_x_0_x_xxyy[k] = -g_z_0_x_0_0_xxyy[k] * ab_x + g_z_0_x_0_0_xxxyy[k];

                g_z_0_x_0_x_xxyz[k] = -g_z_0_x_0_0_xxyz[k] * ab_x + g_z_0_x_0_0_xxxyz[k];

                g_z_0_x_0_x_xxzz[k] = -g_z_0_x_0_0_xxzz[k] * ab_x + g_z_0_x_0_0_xxxzz[k];

                g_z_0_x_0_x_xyyy[k] = -g_z_0_x_0_0_xyyy[k] * ab_x + g_z_0_x_0_0_xxyyy[k];

                g_z_0_x_0_x_xyyz[k] = -g_z_0_x_0_0_xyyz[k] * ab_x + g_z_0_x_0_0_xxyyz[k];

                g_z_0_x_0_x_xyzz[k] = -g_z_0_x_0_0_xyzz[k] * ab_x + g_z_0_x_0_0_xxyzz[k];

                g_z_0_x_0_x_xzzz[k] = -g_z_0_x_0_0_xzzz[k] * ab_x + g_z_0_x_0_0_xxzzz[k];

                g_z_0_x_0_x_yyyy[k] = -g_z_0_x_0_0_yyyy[k] * ab_x + g_z_0_x_0_0_xyyyy[k];

                g_z_0_x_0_x_yyyz[k] = -g_z_0_x_0_0_yyyz[k] * ab_x + g_z_0_x_0_0_xyyyz[k];

                g_z_0_x_0_x_yyzz[k] = -g_z_0_x_0_0_yyzz[k] * ab_x + g_z_0_x_0_0_xyyzz[k];

                g_z_0_x_0_x_yzzz[k] = -g_z_0_x_0_0_yzzz[k] * ab_x + g_z_0_x_0_0_xyzzz[k];

                g_z_0_x_0_x_zzzz[k] = -g_z_0_x_0_0_zzzz[k] * ab_x + g_z_0_x_0_0_xzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 285 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 286 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 287 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 288 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 289 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 290 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 291 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 292 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 293 * ccomps * dcomps);

            auto g_z_0_x_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 294 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 295 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 296 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 297 * ccomps * dcomps);

            auto g_z_0_x_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 298 * ccomps * dcomps);

            auto g_z_0_x_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xxxx, g_z_0_x_0_0_xxxxy, g_z_0_x_0_0_xxxy, g_z_0_x_0_0_xxxyy, g_z_0_x_0_0_xxxyz, g_z_0_x_0_0_xxxz, g_z_0_x_0_0_xxyy, g_z_0_x_0_0_xxyyy, g_z_0_x_0_0_xxyyz, g_z_0_x_0_0_xxyz, g_z_0_x_0_0_xxyzz, g_z_0_x_0_0_xxzz, g_z_0_x_0_0_xyyy, g_z_0_x_0_0_xyyyy, g_z_0_x_0_0_xyyyz, g_z_0_x_0_0_xyyz, g_z_0_x_0_0_xyyzz, g_z_0_x_0_0_xyzz, g_z_0_x_0_0_xyzzz, g_z_0_x_0_0_xzzz, g_z_0_x_0_0_yyyy, g_z_0_x_0_0_yyyyy, g_z_0_x_0_0_yyyyz, g_z_0_x_0_0_yyyz, g_z_0_x_0_0_yyyzz, g_z_0_x_0_0_yyzz, g_z_0_x_0_0_yyzzz, g_z_0_x_0_0_yzzz, g_z_0_x_0_0_yzzzz, g_z_0_x_0_0_zzzz, g_z_0_x_0_y_xxxx, g_z_0_x_0_y_xxxy, g_z_0_x_0_y_xxxz, g_z_0_x_0_y_xxyy, g_z_0_x_0_y_xxyz, g_z_0_x_0_y_xxzz, g_z_0_x_0_y_xyyy, g_z_0_x_0_y_xyyz, g_z_0_x_0_y_xyzz, g_z_0_x_0_y_xzzz, g_z_0_x_0_y_yyyy, g_z_0_x_0_y_yyyz, g_z_0_x_0_y_yyzz, g_z_0_x_0_y_yzzz, g_z_0_x_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_y_xxxx[k] = -g_z_0_x_0_0_xxxx[k] * ab_y + g_z_0_x_0_0_xxxxy[k];

                g_z_0_x_0_y_xxxy[k] = -g_z_0_x_0_0_xxxy[k] * ab_y + g_z_0_x_0_0_xxxyy[k];

                g_z_0_x_0_y_xxxz[k] = -g_z_0_x_0_0_xxxz[k] * ab_y + g_z_0_x_0_0_xxxyz[k];

                g_z_0_x_0_y_xxyy[k] = -g_z_0_x_0_0_xxyy[k] * ab_y + g_z_0_x_0_0_xxyyy[k];

                g_z_0_x_0_y_xxyz[k] = -g_z_0_x_0_0_xxyz[k] * ab_y + g_z_0_x_0_0_xxyyz[k];

                g_z_0_x_0_y_xxzz[k] = -g_z_0_x_0_0_xxzz[k] * ab_y + g_z_0_x_0_0_xxyzz[k];

                g_z_0_x_0_y_xyyy[k] = -g_z_0_x_0_0_xyyy[k] * ab_y + g_z_0_x_0_0_xyyyy[k];

                g_z_0_x_0_y_xyyz[k] = -g_z_0_x_0_0_xyyz[k] * ab_y + g_z_0_x_0_0_xyyyz[k];

                g_z_0_x_0_y_xyzz[k] = -g_z_0_x_0_0_xyzz[k] * ab_y + g_z_0_x_0_0_xyyzz[k];

                g_z_0_x_0_y_xzzz[k] = -g_z_0_x_0_0_xzzz[k] * ab_y + g_z_0_x_0_0_xyzzz[k];

                g_z_0_x_0_y_yyyy[k] = -g_z_0_x_0_0_yyyy[k] * ab_y + g_z_0_x_0_0_yyyyy[k];

                g_z_0_x_0_y_yyyz[k] = -g_z_0_x_0_0_yyyz[k] * ab_y + g_z_0_x_0_0_yyyyz[k];

                g_z_0_x_0_y_yyzz[k] = -g_z_0_x_0_0_yyzz[k] * ab_y + g_z_0_x_0_0_yyyzz[k];

                g_z_0_x_0_y_yzzz[k] = -g_z_0_x_0_0_yzzz[k] * ab_y + g_z_0_x_0_0_yyzzz[k];

                g_z_0_x_0_y_zzzz[k] = -g_z_0_x_0_0_zzzz[k] * ab_y + g_z_0_x_0_0_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 300 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 301 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 302 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 303 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 304 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 305 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 306 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 307 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 308 * ccomps * dcomps);

            auto g_z_0_x_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 309 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 310 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 311 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 312 * ccomps * dcomps);

            auto g_z_0_x_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 313 * ccomps * dcomps);

            auto g_z_0_x_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxx, g_0_0_x_0_0_xxxy, g_0_0_x_0_0_xxxz, g_0_0_x_0_0_xxyy, g_0_0_x_0_0_xxyz, g_0_0_x_0_0_xxzz, g_0_0_x_0_0_xyyy, g_0_0_x_0_0_xyyz, g_0_0_x_0_0_xyzz, g_0_0_x_0_0_xzzz, g_0_0_x_0_0_yyyy, g_0_0_x_0_0_yyyz, g_0_0_x_0_0_yyzz, g_0_0_x_0_0_yzzz, g_0_0_x_0_0_zzzz, g_z_0_x_0_0_xxxx, g_z_0_x_0_0_xxxxz, g_z_0_x_0_0_xxxy, g_z_0_x_0_0_xxxyz, g_z_0_x_0_0_xxxz, g_z_0_x_0_0_xxxzz, g_z_0_x_0_0_xxyy, g_z_0_x_0_0_xxyyz, g_z_0_x_0_0_xxyz, g_z_0_x_0_0_xxyzz, g_z_0_x_0_0_xxzz, g_z_0_x_0_0_xxzzz, g_z_0_x_0_0_xyyy, g_z_0_x_0_0_xyyyz, g_z_0_x_0_0_xyyz, g_z_0_x_0_0_xyyzz, g_z_0_x_0_0_xyzz, g_z_0_x_0_0_xyzzz, g_z_0_x_0_0_xzzz, g_z_0_x_0_0_xzzzz, g_z_0_x_0_0_yyyy, g_z_0_x_0_0_yyyyz, g_z_0_x_0_0_yyyz, g_z_0_x_0_0_yyyzz, g_z_0_x_0_0_yyzz, g_z_0_x_0_0_yyzzz, g_z_0_x_0_0_yzzz, g_z_0_x_0_0_yzzzz, g_z_0_x_0_0_zzzz, g_z_0_x_0_0_zzzzz, g_z_0_x_0_z_xxxx, g_z_0_x_0_z_xxxy, g_z_0_x_0_z_xxxz, g_z_0_x_0_z_xxyy, g_z_0_x_0_z_xxyz, g_z_0_x_0_z_xxzz, g_z_0_x_0_z_xyyy, g_z_0_x_0_z_xyyz, g_z_0_x_0_z_xyzz, g_z_0_x_0_z_xzzz, g_z_0_x_0_z_yyyy, g_z_0_x_0_z_yyyz, g_z_0_x_0_z_yyzz, g_z_0_x_0_z_yzzz, g_z_0_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_z_xxxx[k] = -g_0_0_x_0_0_xxxx[k] - g_z_0_x_0_0_xxxx[k] * ab_z + g_z_0_x_0_0_xxxxz[k];

                g_z_0_x_0_z_xxxy[k] = -g_0_0_x_0_0_xxxy[k] - g_z_0_x_0_0_xxxy[k] * ab_z + g_z_0_x_0_0_xxxyz[k];

                g_z_0_x_0_z_xxxz[k] = -g_0_0_x_0_0_xxxz[k] - g_z_0_x_0_0_xxxz[k] * ab_z + g_z_0_x_0_0_xxxzz[k];

                g_z_0_x_0_z_xxyy[k] = -g_0_0_x_0_0_xxyy[k] - g_z_0_x_0_0_xxyy[k] * ab_z + g_z_0_x_0_0_xxyyz[k];

                g_z_0_x_0_z_xxyz[k] = -g_0_0_x_0_0_xxyz[k] - g_z_0_x_0_0_xxyz[k] * ab_z + g_z_0_x_0_0_xxyzz[k];

                g_z_0_x_0_z_xxzz[k] = -g_0_0_x_0_0_xxzz[k] - g_z_0_x_0_0_xxzz[k] * ab_z + g_z_0_x_0_0_xxzzz[k];

                g_z_0_x_0_z_xyyy[k] = -g_0_0_x_0_0_xyyy[k] - g_z_0_x_0_0_xyyy[k] * ab_z + g_z_0_x_0_0_xyyyz[k];

                g_z_0_x_0_z_xyyz[k] = -g_0_0_x_0_0_xyyz[k] - g_z_0_x_0_0_xyyz[k] * ab_z + g_z_0_x_0_0_xyyzz[k];

                g_z_0_x_0_z_xyzz[k] = -g_0_0_x_0_0_xyzz[k] - g_z_0_x_0_0_xyzz[k] * ab_z + g_z_0_x_0_0_xyzzz[k];

                g_z_0_x_0_z_xzzz[k] = -g_0_0_x_0_0_xzzz[k] - g_z_0_x_0_0_xzzz[k] * ab_z + g_z_0_x_0_0_xzzzz[k];

                g_z_0_x_0_z_yyyy[k] = -g_0_0_x_0_0_yyyy[k] - g_z_0_x_0_0_yyyy[k] * ab_z + g_z_0_x_0_0_yyyyz[k];

                g_z_0_x_0_z_yyyz[k] = -g_0_0_x_0_0_yyyz[k] - g_z_0_x_0_0_yyyz[k] * ab_z + g_z_0_x_0_0_yyyzz[k];

                g_z_0_x_0_z_yyzz[k] = -g_0_0_x_0_0_yyzz[k] - g_z_0_x_0_0_yyzz[k] * ab_z + g_z_0_x_0_0_yyzzz[k];

                g_z_0_x_0_z_yzzz[k] = -g_0_0_x_0_0_yzzz[k] - g_z_0_x_0_0_yzzz[k] * ab_z + g_z_0_x_0_0_yzzzz[k];

                g_z_0_x_0_z_zzzz[k] = -g_0_0_x_0_0_zzzz[k] - g_z_0_x_0_0_zzzz[k] * ab_z + g_z_0_x_0_0_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 315 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 316 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 317 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 318 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 319 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 320 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 321 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 322 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 323 * ccomps * dcomps);

            auto g_z_0_y_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 324 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 325 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 326 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 327 * ccomps * dcomps);

            auto g_z_0_y_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 328 * ccomps * dcomps);

            auto g_z_0_y_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xxxx, g_z_0_y_0_0_xxxxx, g_z_0_y_0_0_xxxxy, g_z_0_y_0_0_xxxxz, g_z_0_y_0_0_xxxy, g_z_0_y_0_0_xxxyy, g_z_0_y_0_0_xxxyz, g_z_0_y_0_0_xxxz, g_z_0_y_0_0_xxxzz, g_z_0_y_0_0_xxyy, g_z_0_y_0_0_xxyyy, g_z_0_y_0_0_xxyyz, g_z_0_y_0_0_xxyz, g_z_0_y_0_0_xxyzz, g_z_0_y_0_0_xxzz, g_z_0_y_0_0_xxzzz, g_z_0_y_0_0_xyyy, g_z_0_y_0_0_xyyyy, g_z_0_y_0_0_xyyyz, g_z_0_y_0_0_xyyz, g_z_0_y_0_0_xyyzz, g_z_0_y_0_0_xyzz, g_z_0_y_0_0_xyzzz, g_z_0_y_0_0_xzzz, g_z_0_y_0_0_xzzzz, g_z_0_y_0_0_yyyy, g_z_0_y_0_0_yyyz, g_z_0_y_0_0_yyzz, g_z_0_y_0_0_yzzz, g_z_0_y_0_0_zzzz, g_z_0_y_0_x_xxxx, g_z_0_y_0_x_xxxy, g_z_0_y_0_x_xxxz, g_z_0_y_0_x_xxyy, g_z_0_y_0_x_xxyz, g_z_0_y_0_x_xxzz, g_z_0_y_0_x_xyyy, g_z_0_y_0_x_xyyz, g_z_0_y_0_x_xyzz, g_z_0_y_0_x_xzzz, g_z_0_y_0_x_yyyy, g_z_0_y_0_x_yyyz, g_z_0_y_0_x_yyzz, g_z_0_y_0_x_yzzz, g_z_0_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_x_xxxx[k] = -g_z_0_y_0_0_xxxx[k] * ab_x + g_z_0_y_0_0_xxxxx[k];

                g_z_0_y_0_x_xxxy[k] = -g_z_0_y_0_0_xxxy[k] * ab_x + g_z_0_y_0_0_xxxxy[k];

                g_z_0_y_0_x_xxxz[k] = -g_z_0_y_0_0_xxxz[k] * ab_x + g_z_0_y_0_0_xxxxz[k];

                g_z_0_y_0_x_xxyy[k] = -g_z_0_y_0_0_xxyy[k] * ab_x + g_z_0_y_0_0_xxxyy[k];

                g_z_0_y_0_x_xxyz[k] = -g_z_0_y_0_0_xxyz[k] * ab_x + g_z_0_y_0_0_xxxyz[k];

                g_z_0_y_0_x_xxzz[k] = -g_z_0_y_0_0_xxzz[k] * ab_x + g_z_0_y_0_0_xxxzz[k];

                g_z_0_y_0_x_xyyy[k] = -g_z_0_y_0_0_xyyy[k] * ab_x + g_z_0_y_0_0_xxyyy[k];

                g_z_0_y_0_x_xyyz[k] = -g_z_0_y_0_0_xyyz[k] * ab_x + g_z_0_y_0_0_xxyyz[k];

                g_z_0_y_0_x_xyzz[k] = -g_z_0_y_0_0_xyzz[k] * ab_x + g_z_0_y_0_0_xxyzz[k];

                g_z_0_y_0_x_xzzz[k] = -g_z_0_y_0_0_xzzz[k] * ab_x + g_z_0_y_0_0_xxzzz[k];

                g_z_0_y_0_x_yyyy[k] = -g_z_0_y_0_0_yyyy[k] * ab_x + g_z_0_y_0_0_xyyyy[k];

                g_z_0_y_0_x_yyyz[k] = -g_z_0_y_0_0_yyyz[k] * ab_x + g_z_0_y_0_0_xyyyz[k];

                g_z_0_y_0_x_yyzz[k] = -g_z_0_y_0_0_yyzz[k] * ab_x + g_z_0_y_0_0_xyyzz[k];

                g_z_0_y_0_x_yzzz[k] = -g_z_0_y_0_0_yzzz[k] * ab_x + g_z_0_y_0_0_xyzzz[k];

                g_z_0_y_0_x_zzzz[k] = -g_z_0_y_0_0_zzzz[k] * ab_x + g_z_0_y_0_0_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 330 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 331 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 332 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 333 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 334 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 335 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 336 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 337 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 338 * ccomps * dcomps);

            auto g_z_0_y_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 339 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 340 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 341 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 342 * ccomps * dcomps);

            auto g_z_0_y_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 343 * ccomps * dcomps);

            auto g_z_0_y_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xxxx, g_z_0_y_0_0_xxxxy, g_z_0_y_0_0_xxxy, g_z_0_y_0_0_xxxyy, g_z_0_y_0_0_xxxyz, g_z_0_y_0_0_xxxz, g_z_0_y_0_0_xxyy, g_z_0_y_0_0_xxyyy, g_z_0_y_0_0_xxyyz, g_z_0_y_0_0_xxyz, g_z_0_y_0_0_xxyzz, g_z_0_y_0_0_xxzz, g_z_0_y_0_0_xyyy, g_z_0_y_0_0_xyyyy, g_z_0_y_0_0_xyyyz, g_z_0_y_0_0_xyyz, g_z_0_y_0_0_xyyzz, g_z_0_y_0_0_xyzz, g_z_0_y_0_0_xyzzz, g_z_0_y_0_0_xzzz, g_z_0_y_0_0_yyyy, g_z_0_y_0_0_yyyyy, g_z_0_y_0_0_yyyyz, g_z_0_y_0_0_yyyz, g_z_0_y_0_0_yyyzz, g_z_0_y_0_0_yyzz, g_z_0_y_0_0_yyzzz, g_z_0_y_0_0_yzzz, g_z_0_y_0_0_yzzzz, g_z_0_y_0_0_zzzz, g_z_0_y_0_y_xxxx, g_z_0_y_0_y_xxxy, g_z_0_y_0_y_xxxz, g_z_0_y_0_y_xxyy, g_z_0_y_0_y_xxyz, g_z_0_y_0_y_xxzz, g_z_0_y_0_y_xyyy, g_z_0_y_0_y_xyyz, g_z_0_y_0_y_xyzz, g_z_0_y_0_y_xzzz, g_z_0_y_0_y_yyyy, g_z_0_y_0_y_yyyz, g_z_0_y_0_y_yyzz, g_z_0_y_0_y_yzzz, g_z_0_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_y_xxxx[k] = -g_z_0_y_0_0_xxxx[k] * ab_y + g_z_0_y_0_0_xxxxy[k];

                g_z_0_y_0_y_xxxy[k] = -g_z_0_y_0_0_xxxy[k] * ab_y + g_z_0_y_0_0_xxxyy[k];

                g_z_0_y_0_y_xxxz[k] = -g_z_0_y_0_0_xxxz[k] * ab_y + g_z_0_y_0_0_xxxyz[k];

                g_z_0_y_0_y_xxyy[k] = -g_z_0_y_0_0_xxyy[k] * ab_y + g_z_0_y_0_0_xxyyy[k];

                g_z_0_y_0_y_xxyz[k] = -g_z_0_y_0_0_xxyz[k] * ab_y + g_z_0_y_0_0_xxyyz[k];

                g_z_0_y_0_y_xxzz[k] = -g_z_0_y_0_0_xxzz[k] * ab_y + g_z_0_y_0_0_xxyzz[k];

                g_z_0_y_0_y_xyyy[k] = -g_z_0_y_0_0_xyyy[k] * ab_y + g_z_0_y_0_0_xyyyy[k];

                g_z_0_y_0_y_xyyz[k] = -g_z_0_y_0_0_xyyz[k] * ab_y + g_z_0_y_0_0_xyyyz[k];

                g_z_0_y_0_y_xyzz[k] = -g_z_0_y_0_0_xyzz[k] * ab_y + g_z_0_y_0_0_xyyzz[k];

                g_z_0_y_0_y_xzzz[k] = -g_z_0_y_0_0_xzzz[k] * ab_y + g_z_0_y_0_0_xyzzz[k];

                g_z_0_y_0_y_yyyy[k] = -g_z_0_y_0_0_yyyy[k] * ab_y + g_z_0_y_0_0_yyyyy[k];

                g_z_0_y_0_y_yyyz[k] = -g_z_0_y_0_0_yyyz[k] * ab_y + g_z_0_y_0_0_yyyyz[k];

                g_z_0_y_0_y_yyzz[k] = -g_z_0_y_0_0_yyzz[k] * ab_y + g_z_0_y_0_0_yyyzz[k];

                g_z_0_y_0_y_yzzz[k] = -g_z_0_y_0_0_yzzz[k] * ab_y + g_z_0_y_0_0_yyzzz[k];

                g_z_0_y_0_y_zzzz[k] = -g_z_0_y_0_0_zzzz[k] * ab_y + g_z_0_y_0_0_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 345 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 346 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 347 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 348 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 349 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 350 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 351 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 352 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 353 * ccomps * dcomps);

            auto g_z_0_y_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 354 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 355 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 356 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 357 * ccomps * dcomps);

            auto g_z_0_y_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 358 * ccomps * dcomps);

            auto g_z_0_y_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxx, g_0_0_y_0_0_xxxy, g_0_0_y_0_0_xxxz, g_0_0_y_0_0_xxyy, g_0_0_y_0_0_xxyz, g_0_0_y_0_0_xxzz, g_0_0_y_0_0_xyyy, g_0_0_y_0_0_xyyz, g_0_0_y_0_0_xyzz, g_0_0_y_0_0_xzzz, g_0_0_y_0_0_yyyy, g_0_0_y_0_0_yyyz, g_0_0_y_0_0_yyzz, g_0_0_y_0_0_yzzz, g_0_0_y_0_0_zzzz, g_z_0_y_0_0_xxxx, g_z_0_y_0_0_xxxxz, g_z_0_y_0_0_xxxy, g_z_0_y_0_0_xxxyz, g_z_0_y_0_0_xxxz, g_z_0_y_0_0_xxxzz, g_z_0_y_0_0_xxyy, g_z_0_y_0_0_xxyyz, g_z_0_y_0_0_xxyz, g_z_0_y_0_0_xxyzz, g_z_0_y_0_0_xxzz, g_z_0_y_0_0_xxzzz, g_z_0_y_0_0_xyyy, g_z_0_y_0_0_xyyyz, g_z_0_y_0_0_xyyz, g_z_0_y_0_0_xyyzz, g_z_0_y_0_0_xyzz, g_z_0_y_0_0_xyzzz, g_z_0_y_0_0_xzzz, g_z_0_y_0_0_xzzzz, g_z_0_y_0_0_yyyy, g_z_0_y_0_0_yyyyz, g_z_0_y_0_0_yyyz, g_z_0_y_0_0_yyyzz, g_z_0_y_0_0_yyzz, g_z_0_y_0_0_yyzzz, g_z_0_y_0_0_yzzz, g_z_0_y_0_0_yzzzz, g_z_0_y_0_0_zzzz, g_z_0_y_0_0_zzzzz, g_z_0_y_0_z_xxxx, g_z_0_y_0_z_xxxy, g_z_0_y_0_z_xxxz, g_z_0_y_0_z_xxyy, g_z_0_y_0_z_xxyz, g_z_0_y_0_z_xxzz, g_z_0_y_0_z_xyyy, g_z_0_y_0_z_xyyz, g_z_0_y_0_z_xyzz, g_z_0_y_0_z_xzzz, g_z_0_y_0_z_yyyy, g_z_0_y_0_z_yyyz, g_z_0_y_0_z_yyzz, g_z_0_y_0_z_yzzz, g_z_0_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_z_xxxx[k] = -g_0_0_y_0_0_xxxx[k] - g_z_0_y_0_0_xxxx[k] * ab_z + g_z_0_y_0_0_xxxxz[k];

                g_z_0_y_0_z_xxxy[k] = -g_0_0_y_0_0_xxxy[k] - g_z_0_y_0_0_xxxy[k] * ab_z + g_z_0_y_0_0_xxxyz[k];

                g_z_0_y_0_z_xxxz[k] = -g_0_0_y_0_0_xxxz[k] - g_z_0_y_0_0_xxxz[k] * ab_z + g_z_0_y_0_0_xxxzz[k];

                g_z_0_y_0_z_xxyy[k] = -g_0_0_y_0_0_xxyy[k] - g_z_0_y_0_0_xxyy[k] * ab_z + g_z_0_y_0_0_xxyyz[k];

                g_z_0_y_0_z_xxyz[k] = -g_0_0_y_0_0_xxyz[k] - g_z_0_y_0_0_xxyz[k] * ab_z + g_z_0_y_0_0_xxyzz[k];

                g_z_0_y_0_z_xxzz[k] = -g_0_0_y_0_0_xxzz[k] - g_z_0_y_0_0_xxzz[k] * ab_z + g_z_0_y_0_0_xxzzz[k];

                g_z_0_y_0_z_xyyy[k] = -g_0_0_y_0_0_xyyy[k] - g_z_0_y_0_0_xyyy[k] * ab_z + g_z_0_y_0_0_xyyyz[k];

                g_z_0_y_0_z_xyyz[k] = -g_0_0_y_0_0_xyyz[k] - g_z_0_y_0_0_xyyz[k] * ab_z + g_z_0_y_0_0_xyyzz[k];

                g_z_0_y_0_z_xyzz[k] = -g_0_0_y_0_0_xyzz[k] - g_z_0_y_0_0_xyzz[k] * ab_z + g_z_0_y_0_0_xyzzz[k];

                g_z_0_y_0_z_xzzz[k] = -g_0_0_y_0_0_xzzz[k] - g_z_0_y_0_0_xzzz[k] * ab_z + g_z_0_y_0_0_xzzzz[k];

                g_z_0_y_0_z_yyyy[k] = -g_0_0_y_0_0_yyyy[k] - g_z_0_y_0_0_yyyy[k] * ab_z + g_z_0_y_0_0_yyyyz[k];

                g_z_0_y_0_z_yyyz[k] = -g_0_0_y_0_0_yyyz[k] - g_z_0_y_0_0_yyyz[k] * ab_z + g_z_0_y_0_0_yyyzz[k];

                g_z_0_y_0_z_yyzz[k] = -g_0_0_y_0_0_yyzz[k] - g_z_0_y_0_0_yyzz[k] * ab_z + g_z_0_y_0_0_yyzzz[k];

                g_z_0_y_0_z_yzzz[k] = -g_0_0_y_0_0_yzzz[k] - g_z_0_y_0_0_yzzz[k] * ab_z + g_z_0_y_0_0_yzzzz[k];

                g_z_0_y_0_z_zzzz[k] = -g_0_0_y_0_0_zzzz[k] - g_z_0_y_0_0_zzzz[k] * ab_z + g_z_0_y_0_0_zzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_x_xxxx = cbuffer.data(pg_geom_1010_off + 360 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxy = cbuffer.data(pg_geom_1010_off + 361 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxz = cbuffer.data(pg_geom_1010_off + 362 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyy = cbuffer.data(pg_geom_1010_off + 363 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyz = cbuffer.data(pg_geom_1010_off + 364 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxzz = cbuffer.data(pg_geom_1010_off + 365 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyy = cbuffer.data(pg_geom_1010_off + 366 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyz = cbuffer.data(pg_geom_1010_off + 367 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyzz = cbuffer.data(pg_geom_1010_off + 368 * ccomps * dcomps);

            auto g_z_0_z_0_x_xzzz = cbuffer.data(pg_geom_1010_off + 369 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyy = cbuffer.data(pg_geom_1010_off + 370 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyz = cbuffer.data(pg_geom_1010_off + 371 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyzz = cbuffer.data(pg_geom_1010_off + 372 * ccomps * dcomps);

            auto g_z_0_z_0_x_yzzz = cbuffer.data(pg_geom_1010_off + 373 * ccomps * dcomps);

            auto g_z_0_z_0_x_zzzz = cbuffer.data(pg_geom_1010_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xxxx, g_z_0_z_0_0_xxxxx, g_z_0_z_0_0_xxxxy, g_z_0_z_0_0_xxxxz, g_z_0_z_0_0_xxxy, g_z_0_z_0_0_xxxyy, g_z_0_z_0_0_xxxyz, g_z_0_z_0_0_xxxz, g_z_0_z_0_0_xxxzz, g_z_0_z_0_0_xxyy, g_z_0_z_0_0_xxyyy, g_z_0_z_0_0_xxyyz, g_z_0_z_0_0_xxyz, g_z_0_z_0_0_xxyzz, g_z_0_z_0_0_xxzz, g_z_0_z_0_0_xxzzz, g_z_0_z_0_0_xyyy, g_z_0_z_0_0_xyyyy, g_z_0_z_0_0_xyyyz, g_z_0_z_0_0_xyyz, g_z_0_z_0_0_xyyzz, g_z_0_z_0_0_xyzz, g_z_0_z_0_0_xyzzz, g_z_0_z_0_0_xzzz, g_z_0_z_0_0_xzzzz, g_z_0_z_0_0_yyyy, g_z_0_z_0_0_yyyz, g_z_0_z_0_0_yyzz, g_z_0_z_0_0_yzzz, g_z_0_z_0_0_zzzz, g_z_0_z_0_x_xxxx, g_z_0_z_0_x_xxxy, g_z_0_z_0_x_xxxz, g_z_0_z_0_x_xxyy, g_z_0_z_0_x_xxyz, g_z_0_z_0_x_xxzz, g_z_0_z_0_x_xyyy, g_z_0_z_0_x_xyyz, g_z_0_z_0_x_xyzz, g_z_0_z_0_x_xzzz, g_z_0_z_0_x_yyyy, g_z_0_z_0_x_yyyz, g_z_0_z_0_x_yyzz, g_z_0_z_0_x_yzzz, g_z_0_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_x_xxxx[k] = -g_z_0_z_0_0_xxxx[k] * ab_x + g_z_0_z_0_0_xxxxx[k];

                g_z_0_z_0_x_xxxy[k] = -g_z_0_z_0_0_xxxy[k] * ab_x + g_z_0_z_0_0_xxxxy[k];

                g_z_0_z_0_x_xxxz[k] = -g_z_0_z_0_0_xxxz[k] * ab_x + g_z_0_z_0_0_xxxxz[k];

                g_z_0_z_0_x_xxyy[k] = -g_z_0_z_0_0_xxyy[k] * ab_x + g_z_0_z_0_0_xxxyy[k];

                g_z_0_z_0_x_xxyz[k] = -g_z_0_z_0_0_xxyz[k] * ab_x + g_z_0_z_0_0_xxxyz[k];

                g_z_0_z_0_x_xxzz[k] = -g_z_0_z_0_0_xxzz[k] * ab_x + g_z_0_z_0_0_xxxzz[k];

                g_z_0_z_0_x_xyyy[k] = -g_z_0_z_0_0_xyyy[k] * ab_x + g_z_0_z_0_0_xxyyy[k];

                g_z_0_z_0_x_xyyz[k] = -g_z_0_z_0_0_xyyz[k] * ab_x + g_z_0_z_0_0_xxyyz[k];

                g_z_0_z_0_x_xyzz[k] = -g_z_0_z_0_0_xyzz[k] * ab_x + g_z_0_z_0_0_xxyzz[k];

                g_z_0_z_0_x_xzzz[k] = -g_z_0_z_0_0_xzzz[k] * ab_x + g_z_0_z_0_0_xxzzz[k];

                g_z_0_z_0_x_yyyy[k] = -g_z_0_z_0_0_yyyy[k] * ab_x + g_z_0_z_0_0_xyyyy[k];

                g_z_0_z_0_x_yyyz[k] = -g_z_0_z_0_0_yyyz[k] * ab_x + g_z_0_z_0_0_xyyyz[k];

                g_z_0_z_0_x_yyzz[k] = -g_z_0_z_0_0_yyzz[k] * ab_x + g_z_0_z_0_0_xyyzz[k];

                g_z_0_z_0_x_yzzz[k] = -g_z_0_z_0_0_yzzz[k] * ab_x + g_z_0_z_0_0_xyzzz[k];

                g_z_0_z_0_x_zzzz[k] = -g_z_0_z_0_0_zzzz[k] * ab_x + g_z_0_z_0_0_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_y_xxxx = cbuffer.data(pg_geom_1010_off + 375 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxy = cbuffer.data(pg_geom_1010_off + 376 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxz = cbuffer.data(pg_geom_1010_off + 377 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyy = cbuffer.data(pg_geom_1010_off + 378 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyz = cbuffer.data(pg_geom_1010_off + 379 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxzz = cbuffer.data(pg_geom_1010_off + 380 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyy = cbuffer.data(pg_geom_1010_off + 381 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyz = cbuffer.data(pg_geom_1010_off + 382 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyzz = cbuffer.data(pg_geom_1010_off + 383 * ccomps * dcomps);

            auto g_z_0_z_0_y_xzzz = cbuffer.data(pg_geom_1010_off + 384 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyy = cbuffer.data(pg_geom_1010_off + 385 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyz = cbuffer.data(pg_geom_1010_off + 386 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyzz = cbuffer.data(pg_geom_1010_off + 387 * ccomps * dcomps);

            auto g_z_0_z_0_y_yzzz = cbuffer.data(pg_geom_1010_off + 388 * ccomps * dcomps);

            auto g_z_0_z_0_y_zzzz = cbuffer.data(pg_geom_1010_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xxxx, g_z_0_z_0_0_xxxxy, g_z_0_z_0_0_xxxy, g_z_0_z_0_0_xxxyy, g_z_0_z_0_0_xxxyz, g_z_0_z_0_0_xxxz, g_z_0_z_0_0_xxyy, g_z_0_z_0_0_xxyyy, g_z_0_z_0_0_xxyyz, g_z_0_z_0_0_xxyz, g_z_0_z_0_0_xxyzz, g_z_0_z_0_0_xxzz, g_z_0_z_0_0_xyyy, g_z_0_z_0_0_xyyyy, g_z_0_z_0_0_xyyyz, g_z_0_z_0_0_xyyz, g_z_0_z_0_0_xyyzz, g_z_0_z_0_0_xyzz, g_z_0_z_0_0_xyzzz, g_z_0_z_0_0_xzzz, g_z_0_z_0_0_yyyy, g_z_0_z_0_0_yyyyy, g_z_0_z_0_0_yyyyz, g_z_0_z_0_0_yyyz, g_z_0_z_0_0_yyyzz, g_z_0_z_0_0_yyzz, g_z_0_z_0_0_yyzzz, g_z_0_z_0_0_yzzz, g_z_0_z_0_0_yzzzz, g_z_0_z_0_0_zzzz, g_z_0_z_0_y_xxxx, g_z_0_z_0_y_xxxy, g_z_0_z_0_y_xxxz, g_z_0_z_0_y_xxyy, g_z_0_z_0_y_xxyz, g_z_0_z_0_y_xxzz, g_z_0_z_0_y_xyyy, g_z_0_z_0_y_xyyz, g_z_0_z_0_y_xyzz, g_z_0_z_0_y_xzzz, g_z_0_z_0_y_yyyy, g_z_0_z_0_y_yyyz, g_z_0_z_0_y_yyzz, g_z_0_z_0_y_yzzz, g_z_0_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_y_xxxx[k] = -g_z_0_z_0_0_xxxx[k] * ab_y + g_z_0_z_0_0_xxxxy[k];

                g_z_0_z_0_y_xxxy[k] = -g_z_0_z_0_0_xxxy[k] * ab_y + g_z_0_z_0_0_xxxyy[k];

                g_z_0_z_0_y_xxxz[k] = -g_z_0_z_0_0_xxxz[k] * ab_y + g_z_0_z_0_0_xxxyz[k];

                g_z_0_z_0_y_xxyy[k] = -g_z_0_z_0_0_xxyy[k] * ab_y + g_z_0_z_0_0_xxyyy[k];

                g_z_0_z_0_y_xxyz[k] = -g_z_0_z_0_0_xxyz[k] * ab_y + g_z_0_z_0_0_xxyyz[k];

                g_z_0_z_0_y_xxzz[k] = -g_z_0_z_0_0_xxzz[k] * ab_y + g_z_0_z_0_0_xxyzz[k];

                g_z_0_z_0_y_xyyy[k] = -g_z_0_z_0_0_xyyy[k] * ab_y + g_z_0_z_0_0_xyyyy[k];

                g_z_0_z_0_y_xyyz[k] = -g_z_0_z_0_0_xyyz[k] * ab_y + g_z_0_z_0_0_xyyyz[k];

                g_z_0_z_0_y_xyzz[k] = -g_z_0_z_0_0_xyzz[k] * ab_y + g_z_0_z_0_0_xyyzz[k];

                g_z_0_z_0_y_xzzz[k] = -g_z_0_z_0_0_xzzz[k] * ab_y + g_z_0_z_0_0_xyzzz[k];

                g_z_0_z_0_y_yyyy[k] = -g_z_0_z_0_0_yyyy[k] * ab_y + g_z_0_z_0_0_yyyyy[k];

                g_z_0_z_0_y_yyyz[k] = -g_z_0_z_0_0_yyyz[k] * ab_y + g_z_0_z_0_0_yyyyz[k];

                g_z_0_z_0_y_yyzz[k] = -g_z_0_z_0_0_yyzz[k] * ab_y + g_z_0_z_0_0_yyyzz[k];

                g_z_0_z_0_y_yzzz[k] = -g_z_0_z_0_0_yzzz[k] * ab_y + g_z_0_z_0_0_yyzzz[k];

                g_z_0_z_0_y_zzzz[k] = -g_z_0_z_0_0_zzzz[k] * ab_y + g_z_0_z_0_0_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_z_xxxx = cbuffer.data(pg_geom_1010_off + 390 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxy = cbuffer.data(pg_geom_1010_off + 391 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxz = cbuffer.data(pg_geom_1010_off + 392 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyy = cbuffer.data(pg_geom_1010_off + 393 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyz = cbuffer.data(pg_geom_1010_off + 394 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxzz = cbuffer.data(pg_geom_1010_off + 395 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyy = cbuffer.data(pg_geom_1010_off + 396 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyz = cbuffer.data(pg_geom_1010_off + 397 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyzz = cbuffer.data(pg_geom_1010_off + 398 * ccomps * dcomps);

            auto g_z_0_z_0_z_xzzz = cbuffer.data(pg_geom_1010_off + 399 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyy = cbuffer.data(pg_geom_1010_off + 400 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyz = cbuffer.data(pg_geom_1010_off + 401 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyzz = cbuffer.data(pg_geom_1010_off + 402 * ccomps * dcomps);

            auto g_z_0_z_0_z_yzzz = cbuffer.data(pg_geom_1010_off + 403 * ccomps * dcomps);

            auto g_z_0_z_0_z_zzzz = cbuffer.data(pg_geom_1010_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxx, g_0_0_z_0_0_xxxy, g_0_0_z_0_0_xxxz, g_0_0_z_0_0_xxyy, g_0_0_z_0_0_xxyz, g_0_0_z_0_0_xxzz, g_0_0_z_0_0_xyyy, g_0_0_z_0_0_xyyz, g_0_0_z_0_0_xyzz, g_0_0_z_0_0_xzzz, g_0_0_z_0_0_yyyy, g_0_0_z_0_0_yyyz, g_0_0_z_0_0_yyzz, g_0_0_z_0_0_yzzz, g_0_0_z_0_0_zzzz, g_z_0_z_0_0_xxxx, g_z_0_z_0_0_xxxxz, g_z_0_z_0_0_xxxy, g_z_0_z_0_0_xxxyz, g_z_0_z_0_0_xxxz, g_z_0_z_0_0_xxxzz, g_z_0_z_0_0_xxyy, g_z_0_z_0_0_xxyyz, g_z_0_z_0_0_xxyz, g_z_0_z_0_0_xxyzz, g_z_0_z_0_0_xxzz, g_z_0_z_0_0_xxzzz, g_z_0_z_0_0_xyyy, g_z_0_z_0_0_xyyyz, g_z_0_z_0_0_xyyz, g_z_0_z_0_0_xyyzz, g_z_0_z_0_0_xyzz, g_z_0_z_0_0_xyzzz, g_z_0_z_0_0_xzzz, g_z_0_z_0_0_xzzzz, g_z_0_z_0_0_yyyy, g_z_0_z_0_0_yyyyz, g_z_0_z_0_0_yyyz, g_z_0_z_0_0_yyyzz, g_z_0_z_0_0_yyzz, g_z_0_z_0_0_yyzzz, g_z_0_z_0_0_yzzz, g_z_0_z_0_0_yzzzz, g_z_0_z_0_0_zzzz, g_z_0_z_0_0_zzzzz, g_z_0_z_0_z_xxxx, g_z_0_z_0_z_xxxy, g_z_0_z_0_z_xxxz, g_z_0_z_0_z_xxyy, g_z_0_z_0_z_xxyz, g_z_0_z_0_z_xxzz, g_z_0_z_0_z_xyyy, g_z_0_z_0_z_xyyz, g_z_0_z_0_z_xyzz, g_z_0_z_0_z_xzzz, g_z_0_z_0_z_yyyy, g_z_0_z_0_z_yyyz, g_z_0_z_0_z_yyzz, g_z_0_z_0_z_yzzz, g_z_0_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_z_xxxx[k] = -g_0_0_z_0_0_xxxx[k] - g_z_0_z_0_0_xxxx[k] * ab_z + g_z_0_z_0_0_xxxxz[k];

                g_z_0_z_0_z_xxxy[k] = -g_0_0_z_0_0_xxxy[k] - g_z_0_z_0_0_xxxy[k] * ab_z + g_z_0_z_0_0_xxxyz[k];

                g_z_0_z_0_z_xxxz[k] = -g_0_0_z_0_0_xxxz[k] - g_z_0_z_0_0_xxxz[k] * ab_z + g_z_0_z_0_0_xxxzz[k];

                g_z_0_z_0_z_xxyy[k] = -g_0_0_z_0_0_xxyy[k] - g_z_0_z_0_0_xxyy[k] * ab_z + g_z_0_z_0_0_xxyyz[k];

                g_z_0_z_0_z_xxyz[k] = -g_0_0_z_0_0_xxyz[k] - g_z_0_z_0_0_xxyz[k] * ab_z + g_z_0_z_0_0_xxyzz[k];

                g_z_0_z_0_z_xxzz[k] = -g_0_0_z_0_0_xxzz[k] - g_z_0_z_0_0_xxzz[k] * ab_z + g_z_0_z_0_0_xxzzz[k];

                g_z_0_z_0_z_xyyy[k] = -g_0_0_z_0_0_xyyy[k] - g_z_0_z_0_0_xyyy[k] * ab_z + g_z_0_z_0_0_xyyyz[k];

                g_z_0_z_0_z_xyyz[k] = -g_0_0_z_0_0_xyyz[k] - g_z_0_z_0_0_xyyz[k] * ab_z + g_z_0_z_0_0_xyyzz[k];

                g_z_0_z_0_z_xyzz[k] = -g_0_0_z_0_0_xyzz[k] - g_z_0_z_0_0_xyzz[k] * ab_z + g_z_0_z_0_0_xyzzz[k];

                g_z_0_z_0_z_xzzz[k] = -g_0_0_z_0_0_xzzz[k] - g_z_0_z_0_0_xzzz[k] * ab_z + g_z_0_z_0_0_xzzzz[k];

                g_z_0_z_0_z_yyyy[k] = -g_0_0_z_0_0_yyyy[k] - g_z_0_z_0_0_yyyy[k] * ab_z + g_z_0_z_0_0_yyyyz[k];

                g_z_0_z_0_z_yyyz[k] = -g_0_0_z_0_0_yyyz[k] - g_z_0_z_0_0_yyyz[k] * ab_z + g_z_0_z_0_0_yyyzz[k];

                g_z_0_z_0_z_yyzz[k] = -g_0_0_z_0_0_yyzz[k] - g_z_0_z_0_0_yyzz[k] * ab_z + g_z_0_z_0_0_yyzzz[k];

                g_z_0_z_0_z_yzzz[k] = -g_0_0_z_0_0_yzzz[k] - g_z_0_z_0_0_yzzz[k] * ab_z + g_z_0_z_0_0_yzzzz[k];

                g_z_0_z_0_z_zzzz[k] = -g_0_0_z_0_0_zzzz[k] - g_z_0_z_0_0_zzzz[k] * ab_z + g_z_0_z_0_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

