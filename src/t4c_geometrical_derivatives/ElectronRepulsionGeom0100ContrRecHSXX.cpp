#include "ElectronRepulsionGeom0100ContrRecHSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_hsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hsxx,
                                            const size_t idx_gsxx,
                                            const size_t idx_geom_01_gsxx,
                                            const size_t idx_geom_01_gpxx,
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
            /// Set up components of auxilary buffer : GSSS

            const auto gs_off = idx_gsxx + i * dcomps + j;

            auto g_xxxx_0 = cbuffer.data(gs_off + 0 * ccomps * dcomps);

            auto g_xxxy_0 = cbuffer.data(gs_off + 1 * ccomps * dcomps);

            auto g_xxxz_0 = cbuffer.data(gs_off + 2 * ccomps * dcomps);

            auto g_xxyy_0 = cbuffer.data(gs_off + 3 * ccomps * dcomps);

            auto g_xxyz_0 = cbuffer.data(gs_off + 4 * ccomps * dcomps);

            auto g_xxzz_0 = cbuffer.data(gs_off + 5 * ccomps * dcomps);

            auto g_xyyy_0 = cbuffer.data(gs_off + 6 * ccomps * dcomps);

            auto g_xyyz_0 = cbuffer.data(gs_off + 7 * ccomps * dcomps);

            auto g_xyzz_0 = cbuffer.data(gs_off + 8 * ccomps * dcomps);

            auto g_xzzz_0 = cbuffer.data(gs_off + 9 * ccomps * dcomps);

            auto g_yyyy_0 = cbuffer.data(gs_off + 10 * ccomps * dcomps);

            auto g_yyyz_0 = cbuffer.data(gs_off + 11 * ccomps * dcomps);

            auto g_yyzz_0 = cbuffer.data(gs_off + 12 * ccomps * dcomps);

            auto g_yzzz_0 = cbuffer.data(gs_off + 13 * ccomps * dcomps);

            auto g_zzzz_0 = cbuffer.data(gs_off + 14 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GSSS

            const auto gs_geom_01_off = idx_geom_01_gsxx + i * dcomps + j;

            auto g_0_x_xxxx_0 = cbuffer.data(gs_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxy_0 = cbuffer.data(gs_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxz_0 = cbuffer.data(gs_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxyy_0 = cbuffer.data(gs_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxyz_0 = cbuffer.data(gs_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxzz_0 = cbuffer.data(gs_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xyyy_0 = cbuffer.data(gs_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xyyz_0 = cbuffer.data(gs_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xyzz_0 = cbuffer.data(gs_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xzzz_0 = cbuffer.data(gs_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_yyyy_0 = cbuffer.data(gs_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_yyyz_0 = cbuffer.data(gs_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_yyzz_0 = cbuffer.data(gs_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_yzzz_0 = cbuffer.data(gs_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_zzzz_0 = cbuffer.data(gs_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_xxxx_0 = cbuffer.data(gs_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_xxxy_0 = cbuffer.data(gs_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_xxxz_0 = cbuffer.data(gs_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_xxyy_0 = cbuffer.data(gs_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_xxyz_0 = cbuffer.data(gs_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_xxzz_0 = cbuffer.data(gs_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_xyyy_0 = cbuffer.data(gs_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_xyyz_0 = cbuffer.data(gs_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_xyzz_0 = cbuffer.data(gs_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_xzzz_0 = cbuffer.data(gs_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_yyyy_0 = cbuffer.data(gs_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_yyyz_0 = cbuffer.data(gs_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_yyzz_0 = cbuffer.data(gs_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_yzzz_0 = cbuffer.data(gs_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_zzzz_0 = cbuffer.data(gs_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_z_xxxx_0 = cbuffer.data(gs_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_z_xxxy_0 = cbuffer.data(gs_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_z_xxxz_0 = cbuffer.data(gs_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_z_xxyy_0 = cbuffer.data(gs_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_z_xxyz_0 = cbuffer.data(gs_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_z_xxzz_0 = cbuffer.data(gs_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_xyyy_0 = cbuffer.data(gs_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_xyyz_0 = cbuffer.data(gs_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_xyzz_0 = cbuffer.data(gs_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_xzzz_0 = cbuffer.data(gs_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_yyyy_0 = cbuffer.data(gs_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_yyyz_0 = cbuffer.data(gs_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_yyzz_0 = cbuffer.data(gs_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_yzzz_0 = cbuffer.data(gs_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_zzzz_0 = cbuffer.data(gs_geom_01_off + 44 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_hsxx

            const auto hs_geom_01_off = idx_geom_01_hsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxx_0 = cbuffer.data(hs_geom_01_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_0, g_0_x_xxxx_x, g_0_x_xxxxx_0, g_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxx_0[k] = g_xxxx_0[k] - g_0_x_xxxx_0[k] * ab_x + g_0_x_xxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxy_0 = cbuffer.data(hs_geom_01_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_0, g_0_x_xxxx_y, g_0_x_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxy_0[k] = -g_0_x_xxxx_0[k] * ab_y + g_0_x_xxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxz_0 = cbuffer.data(hs_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_0, g_0_x_xxxx_z, g_0_x_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxz_0[k] = -g_0_x_xxxx_0[k] * ab_z + g_0_x_xxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyy_0 = cbuffer.data(hs_geom_01_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxy_0, g_0_x_xxxy_y, g_0_x_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyy_0[k] = -g_0_x_xxxy_0[k] * ab_y + g_0_x_xxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyz_0 = cbuffer.data(hs_geom_01_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyz_0, g_0_x_xxxz_0, g_0_x_xxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyz_0[k] = -g_0_x_xxxz_0[k] * ab_y + g_0_x_xxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzz_0 = cbuffer.data(hs_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxz_0, g_0_x_xxxz_z, g_0_x_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzz_0[k] = -g_0_x_xxxz_0[k] * ab_z + g_0_x_xxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyy_0 = cbuffer.data(hs_geom_01_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyy_0, g_0_x_xxyy_y, g_0_x_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyy_0[k] = -g_0_x_xxyy_0[k] * ab_y + g_0_x_xxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyz_0 = cbuffer.data(hs_geom_01_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyz_0, g_0_x_xxyz_0, g_0_x_xxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyz_0[k] = -g_0_x_xxyz_0[k] * ab_y + g_0_x_xxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzz_0 = cbuffer.data(hs_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzz_0, g_0_x_xxzz_0, g_0_x_xxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzz_0[k] = -g_0_x_xxzz_0[k] * ab_y + g_0_x_xxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzz_0 = cbuffer.data(hs_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzz_0, g_0_x_xxzz_z, g_0_x_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzz_0[k] = -g_0_x_xxzz_0[k] * ab_z + g_0_x_xxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyy_0 = cbuffer.data(hs_geom_01_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyy_0, g_0_x_xyyy_y, g_0_x_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyy_0[k] = -g_0_x_xyyy_0[k] * ab_y + g_0_x_xyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyz_0 = cbuffer.data(hs_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyz_0, g_0_x_xyyz_0, g_0_x_xyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyz_0[k] = -g_0_x_xyyz_0[k] * ab_y + g_0_x_xyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzz_0 = cbuffer.data(hs_geom_01_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzz_0, g_0_x_xyzz_0, g_0_x_xyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzz_0[k] = -g_0_x_xyzz_0[k] * ab_y + g_0_x_xyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzz_0 = cbuffer.data(hs_geom_01_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzz_0, g_0_x_xzzz_0, g_0_x_xzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzz_0[k] = -g_0_x_xzzz_0[k] * ab_y + g_0_x_xzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzz_0 = cbuffer.data(hs_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzz_0, g_0_x_xzzz_z, g_0_x_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzz_0[k] = -g_0_x_xzzz_0[k] * ab_z + g_0_x_xzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyy_0 = cbuffer.data(hs_geom_01_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_0, g_0_x_yyyy_y, g_0_x_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyy_0[k] = -g_0_x_yyyy_0[k] * ab_y + g_0_x_yyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyz_0 = cbuffer.data(hs_geom_01_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyz_0, g_0_x_yyyz_0, g_0_x_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyz_0[k] = -g_0_x_yyyz_0[k] * ab_y + g_0_x_yyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzz_0 = cbuffer.data(hs_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzz_0, g_0_x_yyzz_0, g_0_x_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzz_0[k] = -g_0_x_yyzz_0[k] * ab_y + g_0_x_yyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzz_0 = cbuffer.data(hs_geom_01_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzz_0, g_0_x_yzzz_0, g_0_x_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzz_0[k] = -g_0_x_yzzz_0[k] * ab_y + g_0_x_yzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzz_0 = cbuffer.data(hs_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzz_0, g_0_x_zzzz_0, g_0_x_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzz_0[k] = -g_0_x_zzzz_0[k] * ab_y + g_0_x_zzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzz_0 = cbuffer.data(hs_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_0, g_0_x_zzzz_z, g_0_x_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzz_0[k] = -g_0_x_zzzz_0[k] * ab_z + g_0_x_zzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxx_0 = cbuffer.data(hs_geom_01_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_0, g_0_y_xxxx_x, g_0_y_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxx_0[k] = -g_0_y_xxxx_0[k] * ab_x + g_0_y_xxxx_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxy_0 = cbuffer.data(hs_geom_01_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxy_0, g_0_y_xxxy_0, g_0_y_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxy_0[k] = -g_0_y_xxxy_0[k] * ab_x + g_0_y_xxxy_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxz_0 = cbuffer.data(hs_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxz_0, g_0_y_xxxz_0, g_0_y_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxz_0[k] = -g_0_y_xxxz_0[k] * ab_x + g_0_y_xxxz_x[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyy_0 = cbuffer.data(hs_geom_01_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyy_0, g_0_y_xxyy_0, g_0_y_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyy_0[k] = -g_0_y_xxyy_0[k] * ab_x + g_0_y_xxyy_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyz_0 = cbuffer.data(hs_geom_01_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyz_0, g_0_y_xxyz_0, g_0_y_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyz_0[k] = -g_0_y_xxyz_0[k] * ab_x + g_0_y_xxyz_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzz_0 = cbuffer.data(hs_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzz_0, g_0_y_xxzz_0, g_0_y_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzz_0[k] = -g_0_y_xxzz_0[k] * ab_x + g_0_y_xxzz_x[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyy_0 = cbuffer.data(hs_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyy_0, g_0_y_xyyy_0, g_0_y_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyy_0[k] = -g_0_y_xyyy_0[k] * ab_x + g_0_y_xyyy_x[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyz_0 = cbuffer.data(hs_geom_01_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyz_0, g_0_y_xyyz_0, g_0_y_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyz_0[k] = -g_0_y_xyyz_0[k] * ab_x + g_0_y_xyyz_x[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzz_0 = cbuffer.data(hs_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzz_0, g_0_y_xyzz_0, g_0_y_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzz_0[k] = -g_0_y_xyzz_0[k] * ab_x + g_0_y_xyzz_x[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzz_0 = cbuffer.data(hs_geom_01_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzz_0, g_0_y_xzzz_0, g_0_y_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzz_0[k] = -g_0_y_xzzz_0[k] * ab_x + g_0_y_xzzz_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyy_0 = cbuffer.data(hs_geom_01_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyy_0, g_0_y_yyyy_0, g_0_y_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyy_0[k] = -g_0_y_yyyy_0[k] * ab_x + g_0_y_yyyy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyz_0 = cbuffer.data(hs_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyz_0, g_0_y_yyyz_0, g_0_y_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyz_0[k] = -g_0_y_yyyz_0[k] * ab_x + g_0_y_yyyz_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzz_0 = cbuffer.data(hs_geom_01_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzz_0, g_0_y_yyzz_0, g_0_y_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzz_0[k] = -g_0_y_yyzz_0[k] * ab_x + g_0_y_yyzz_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzz_0 = cbuffer.data(hs_geom_01_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzz_0, g_0_y_yzzz_0, g_0_y_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzz_0[k] = -g_0_y_yzzz_0[k] * ab_x + g_0_y_yzzz_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzz_0 = cbuffer.data(hs_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzz_0, g_0_y_zzzz_0, g_0_y_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzz_0[k] = -g_0_y_zzzz_0[k] * ab_x + g_0_y_zzzz_x[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyy_0 = cbuffer.data(hs_geom_01_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_0, g_0_y_yyyy_y, g_0_y_yyyyy_0, g_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyy_0[k] = g_yyyy_0[k] - g_0_y_yyyy_0[k] * ab_y + g_0_y_yyyy_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyz_0 = cbuffer.data(hs_geom_01_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_0, g_0_y_yyyy_z, g_0_y_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyz_0[k] = -g_0_y_yyyy_0[k] * ab_z + g_0_y_yyyy_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzz_0 = cbuffer.data(hs_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyz_0, g_0_y_yyyz_z, g_0_y_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzz_0[k] = -g_0_y_yyyz_0[k] * ab_z + g_0_y_yyyz_z[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzz_0 = cbuffer.data(hs_geom_01_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzz_0, g_0_y_yyzz_z, g_0_y_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzz_0[k] = -g_0_y_yyzz_0[k] * ab_z + g_0_y_yyzz_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzz_0 = cbuffer.data(hs_geom_01_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzz_0, g_0_y_yzzz_z, g_0_y_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzz_0[k] = -g_0_y_yzzz_0[k] * ab_z + g_0_y_yzzz_z[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzz_0 = cbuffer.data(hs_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_0, g_0_y_zzzz_z, g_0_y_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzz_0[k] = -g_0_y_zzzz_0[k] * ab_z + g_0_y_zzzz_z[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxx_0 = cbuffer.data(hs_geom_01_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_0, g_0_z_xxxx_x, g_0_z_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxx_0[k] = -g_0_z_xxxx_0[k] * ab_x + g_0_z_xxxx_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxy_0 = cbuffer.data(hs_geom_01_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxy_0, g_0_z_xxxy_0, g_0_z_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxy_0[k] = -g_0_z_xxxy_0[k] * ab_x + g_0_z_xxxy_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxz_0 = cbuffer.data(hs_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxz_0, g_0_z_xxxz_0, g_0_z_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxz_0[k] = -g_0_z_xxxz_0[k] * ab_x + g_0_z_xxxz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyy_0 = cbuffer.data(hs_geom_01_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyy_0, g_0_z_xxyy_0, g_0_z_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyy_0[k] = -g_0_z_xxyy_0[k] * ab_x + g_0_z_xxyy_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyz_0 = cbuffer.data(hs_geom_01_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyz_0, g_0_z_xxyz_0, g_0_z_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyz_0[k] = -g_0_z_xxyz_0[k] * ab_x + g_0_z_xxyz_x[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzz_0 = cbuffer.data(hs_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzz_0, g_0_z_xxzz_0, g_0_z_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzz_0[k] = -g_0_z_xxzz_0[k] * ab_x + g_0_z_xxzz_x[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyy_0 = cbuffer.data(hs_geom_01_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyy_0, g_0_z_xyyy_0, g_0_z_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyy_0[k] = -g_0_z_xyyy_0[k] * ab_x + g_0_z_xyyy_x[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyz_0 = cbuffer.data(hs_geom_01_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyz_0, g_0_z_xyyz_0, g_0_z_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyz_0[k] = -g_0_z_xyyz_0[k] * ab_x + g_0_z_xyyz_x[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzz_0 = cbuffer.data(hs_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzz_0, g_0_z_xyzz_0, g_0_z_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzz_0[k] = -g_0_z_xyzz_0[k] * ab_x + g_0_z_xyzz_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzz_0 = cbuffer.data(hs_geom_01_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzz_0, g_0_z_xzzz_0, g_0_z_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzz_0[k] = -g_0_z_xzzz_0[k] * ab_x + g_0_z_xzzz_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyy_0 = cbuffer.data(hs_geom_01_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyy_0, g_0_z_yyyy_0, g_0_z_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyy_0[k] = -g_0_z_yyyy_0[k] * ab_x + g_0_z_yyyy_x[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyz_0 = cbuffer.data(hs_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyz_0, g_0_z_yyyz_0, g_0_z_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyz_0[k] = -g_0_z_yyyz_0[k] * ab_x + g_0_z_yyyz_x[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzz_0 = cbuffer.data(hs_geom_01_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzz_0, g_0_z_yyzz_0, g_0_z_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzz_0[k] = -g_0_z_yyzz_0[k] * ab_x + g_0_z_yyzz_x[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzz_0 = cbuffer.data(hs_geom_01_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzz_0, g_0_z_yzzz_0, g_0_z_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzz_0[k] = -g_0_z_yzzz_0[k] * ab_x + g_0_z_yzzz_x[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzz_0 = cbuffer.data(hs_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzz_0, g_0_z_zzzz_0, g_0_z_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzz_0[k] = -g_0_z_zzzz_0[k] * ab_x + g_0_z_zzzz_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyy_0 = cbuffer.data(hs_geom_01_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_0, g_0_z_yyyy_y, g_0_z_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyy_0[k] = -g_0_z_yyyy_0[k] * ab_y + g_0_z_yyyy_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyz_0 = cbuffer.data(hs_geom_01_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyz_0, g_0_z_yyyz_0, g_0_z_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyz_0[k] = -g_0_z_yyyz_0[k] * ab_y + g_0_z_yyyz_y[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzz_0 = cbuffer.data(hs_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzz_0, g_0_z_yyzz_0, g_0_z_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzz_0[k] = -g_0_z_yyzz_0[k] * ab_y + g_0_z_yyzz_y[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzz_0 = cbuffer.data(hs_geom_01_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzz_0, g_0_z_yzzz_0, g_0_z_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzz_0[k] = -g_0_z_yzzz_0[k] * ab_y + g_0_z_yzzz_y[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzz_0 = cbuffer.data(hs_geom_01_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzz_0, g_0_z_zzzz_0, g_0_z_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzz_0[k] = -g_0_z_zzzz_0[k] * ab_y + g_0_z_zzzz_y[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzz_0 = cbuffer.data(hs_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_0, g_0_z_zzzz_z, g_0_z_zzzzz_0, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzz_0[k] = g_zzzz_0[k] - g_0_z_zzzz_0[k] * ab_z + g_0_z_zzzz_z[k];
            }
        }
    }
}

} // erirec namespace

