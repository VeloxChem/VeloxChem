#include "ElectronRepulsionGeom1100ContrRecHSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_hsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_hsxx,
                                            const size_t idx_geom_01_gsxx,
                                            const size_t idx_geom_10_gsxx,
                                            const size_t idx_geom_11_gsxx,
                                            const size_t idx_geom_11_gpxx,
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

            /// Set up components of auxilary buffer : GSSS

            const auto gs_geom_10_off = idx_geom_10_gsxx + i * dcomps + j;

            auto g_x_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 29 * ccomps * dcomps);

            auto g_z_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 30 * ccomps * dcomps);

            auto g_z_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 31 * ccomps * dcomps);

            auto g_z_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 32 * ccomps * dcomps);

            auto g_z_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 33 * ccomps * dcomps);

            auto g_z_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 34 * ccomps * dcomps);

            auto g_z_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 44 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GSSS

            const auto gs_geom_11_off = idx_geom_11_gsxx + i * dcomps + j;

            auto g_x_x_xxxx_0 = cbuffer.data(gs_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxy_0 = cbuffer.data(gs_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxz_0 = cbuffer.data(gs_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxyy_0 = cbuffer.data(gs_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxyz_0 = cbuffer.data(gs_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxzz_0 = cbuffer.data(gs_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xyyy_0 = cbuffer.data(gs_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xyyz_0 = cbuffer.data(gs_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xyzz_0 = cbuffer.data(gs_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xzzz_0 = cbuffer.data(gs_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_yyyy_0 = cbuffer.data(gs_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_yyyz_0 = cbuffer.data(gs_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_yyzz_0 = cbuffer.data(gs_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_yzzz_0 = cbuffer.data(gs_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_zzzz_0 = cbuffer.data(gs_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_y_xxxx_0 = cbuffer.data(gs_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_xxxy_0 = cbuffer.data(gs_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_xxxz_0 = cbuffer.data(gs_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_xxyy_0 = cbuffer.data(gs_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_xxyz_0 = cbuffer.data(gs_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_xxzz_0 = cbuffer.data(gs_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_xyyy_0 = cbuffer.data(gs_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_xyyz_0 = cbuffer.data(gs_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_xyzz_0 = cbuffer.data(gs_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_xzzz_0 = cbuffer.data(gs_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_yyyy_0 = cbuffer.data(gs_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_yyyz_0 = cbuffer.data(gs_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_yyzz_0 = cbuffer.data(gs_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_yzzz_0 = cbuffer.data(gs_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_zzzz_0 = cbuffer.data(gs_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_z_xxxx_0 = cbuffer.data(gs_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_z_xxxy_0 = cbuffer.data(gs_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_z_xxxz_0 = cbuffer.data(gs_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_z_xxyy_0 = cbuffer.data(gs_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_z_xxyz_0 = cbuffer.data(gs_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_z_xxzz_0 = cbuffer.data(gs_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_xyyy_0 = cbuffer.data(gs_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_xyyz_0 = cbuffer.data(gs_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_xyzz_0 = cbuffer.data(gs_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_xzzz_0 = cbuffer.data(gs_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_yyyy_0 = cbuffer.data(gs_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_yyyz_0 = cbuffer.data(gs_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_yyzz_0 = cbuffer.data(gs_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_yzzz_0 = cbuffer.data(gs_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_zzzz_0 = cbuffer.data(gs_geom_11_off + 44 * ccomps * dcomps);

            auto g_y_x_xxxx_0 = cbuffer.data(gs_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_x_xxxy_0 = cbuffer.data(gs_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_x_xxxz_0 = cbuffer.data(gs_geom_11_off + 47 * ccomps * dcomps);

            auto g_y_x_xxyy_0 = cbuffer.data(gs_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_x_xxyz_0 = cbuffer.data(gs_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_x_xxzz_0 = cbuffer.data(gs_geom_11_off + 50 * ccomps * dcomps);

            auto g_y_x_xyyy_0 = cbuffer.data(gs_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_x_xyyz_0 = cbuffer.data(gs_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_x_xyzz_0 = cbuffer.data(gs_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_xzzz_0 = cbuffer.data(gs_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_yyyy_0 = cbuffer.data(gs_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_yyyz_0 = cbuffer.data(gs_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_yyzz_0 = cbuffer.data(gs_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_yzzz_0 = cbuffer.data(gs_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_zzzz_0 = cbuffer.data(gs_geom_11_off + 59 * ccomps * dcomps);

            auto g_y_y_xxxx_0 = cbuffer.data(gs_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_y_xxxy_0 = cbuffer.data(gs_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_y_xxxz_0 = cbuffer.data(gs_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_y_xxyy_0 = cbuffer.data(gs_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_y_xxyz_0 = cbuffer.data(gs_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_y_xxzz_0 = cbuffer.data(gs_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_y_xyyy_0 = cbuffer.data(gs_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_y_xyyz_0 = cbuffer.data(gs_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_y_xyzz_0 = cbuffer.data(gs_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_y_xzzz_0 = cbuffer.data(gs_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_y_yyyy_0 = cbuffer.data(gs_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_y_yyyz_0 = cbuffer.data(gs_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_yyzz_0 = cbuffer.data(gs_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_yzzz_0 = cbuffer.data(gs_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_zzzz_0 = cbuffer.data(gs_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_z_xxxx_0 = cbuffer.data(gs_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_z_xxxy_0 = cbuffer.data(gs_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_z_xxxz_0 = cbuffer.data(gs_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_z_xxyy_0 = cbuffer.data(gs_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_z_xxyz_0 = cbuffer.data(gs_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_z_xxzz_0 = cbuffer.data(gs_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_z_xyyy_0 = cbuffer.data(gs_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_z_xyyz_0 = cbuffer.data(gs_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_z_xyzz_0 = cbuffer.data(gs_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_z_xzzz_0 = cbuffer.data(gs_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_z_yyyy_0 = cbuffer.data(gs_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_z_yyyz_0 = cbuffer.data(gs_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_z_yyzz_0 = cbuffer.data(gs_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_z_yzzz_0 = cbuffer.data(gs_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_z_zzzz_0 = cbuffer.data(gs_geom_11_off + 89 * ccomps * dcomps);

            auto g_z_x_xxxx_0 = cbuffer.data(gs_geom_11_off + 90 * ccomps * dcomps);

            auto g_z_x_xxxy_0 = cbuffer.data(gs_geom_11_off + 91 * ccomps * dcomps);

            auto g_z_x_xxxz_0 = cbuffer.data(gs_geom_11_off + 92 * ccomps * dcomps);

            auto g_z_x_xxyy_0 = cbuffer.data(gs_geom_11_off + 93 * ccomps * dcomps);

            auto g_z_x_xxyz_0 = cbuffer.data(gs_geom_11_off + 94 * ccomps * dcomps);

            auto g_z_x_xxzz_0 = cbuffer.data(gs_geom_11_off + 95 * ccomps * dcomps);

            auto g_z_x_xyyy_0 = cbuffer.data(gs_geom_11_off + 96 * ccomps * dcomps);

            auto g_z_x_xyyz_0 = cbuffer.data(gs_geom_11_off + 97 * ccomps * dcomps);

            auto g_z_x_xyzz_0 = cbuffer.data(gs_geom_11_off + 98 * ccomps * dcomps);

            auto g_z_x_xzzz_0 = cbuffer.data(gs_geom_11_off + 99 * ccomps * dcomps);

            auto g_z_x_yyyy_0 = cbuffer.data(gs_geom_11_off + 100 * ccomps * dcomps);

            auto g_z_x_yyyz_0 = cbuffer.data(gs_geom_11_off + 101 * ccomps * dcomps);

            auto g_z_x_yyzz_0 = cbuffer.data(gs_geom_11_off + 102 * ccomps * dcomps);

            auto g_z_x_yzzz_0 = cbuffer.data(gs_geom_11_off + 103 * ccomps * dcomps);

            auto g_z_x_zzzz_0 = cbuffer.data(gs_geom_11_off + 104 * ccomps * dcomps);

            auto g_z_y_xxxx_0 = cbuffer.data(gs_geom_11_off + 105 * ccomps * dcomps);

            auto g_z_y_xxxy_0 = cbuffer.data(gs_geom_11_off + 106 * ccomps * dcomps);

            auto g_z_y_xxxz_0 = cbuffer.data(gs_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_y_xxyy_0 = cbuffer.data(gs_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_y_xxyz_0 = cbuffer.data(gs_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_y_xxzz_0 = cbuffer.data(gs_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_y_xyyy_0 = cbuffer.data(gs_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_y_xyyz_0 = cbuffer.data(gs_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_y_xyzz_0 = cbuffer.data(gs_geom_11_off + 113 * ccomps * dcomps);

            auto g_z_y_xzzz_0 = cbuffer.data(gs_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_y_yyyy_0 = cbuffer.data(gs_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_y_yyyz_0 = cbuffer.data(gs_geom_11_off + 116 * ccomps * dcomps);

            auto g_z_y_yyzz_0 = cbuffer.data(gs_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_y_yzzz_0 = cbuffer.data(gs_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_y_zzzz_0 = cbuffer.data(gs_geom_11_off + 119 * ccomps * dcomps);

            auto g_z_z_xxxx_0 = cbuffer.data(gs_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_z_xxxy_0 = cbuffer.data(gs_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_z_xxxz_0 = cbuffer.data(gs_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_z_xxyy_0 = cbuffer.data(gs_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_z_xxyz_0 = cbuffer.data(gs_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_z_xxzz_0 = cbuffer.data(gs_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_z_xyyy_0 = cbuffer.data(gs_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_z_xyyz_0 = cbuffer.data(gs_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_z_xyzz_0 = cbuffer.data(gs_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_z_xzzz_0 = cbuffer.data(gs_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_z_yyyy_0 = cbuffer.data(gs_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_z_yyyz_0 = cbuffer.data(gs_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_z_yyzz_0 = cbuffer.data(gs_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_z_yzzz_0 = cbuffer.data(gs_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_z_zzzz_0 = cbuffer.data(gs_geom_11_off + 134 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GPSS

            const auto gp_geom_11_off = idx_geom_11_gpxx + i * dcomps + j;

            auto g_x_x_xxxx_x = cbuffer.data(gp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxx_y = cbuffer.data(gp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxx_z = cbuffer.data(gp_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxxy_x = cbuffer.data(gp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxy_y = cbuffer.data(gp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxy_z = cbuffer.data(gp_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxxz_x = cbuffer.data(gp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxz_y = cbuffer.data(gp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxz_z = cbuffer.data(gp_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxyy_x = cbuffer.data(gp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxyy_y = cbuffer.data(gp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxyy_z = cbuffer.data(gp_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxyz_x = cbuffer.data(gp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxyz_y = cbuffer.data(gp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxyz_z = cbuffer.data(gp_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxzz_x = cbuffer.data(gp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxzz_y = cbuffer.data(gp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxzz_z = cbuffer.data(gp_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xyyy_x = cbuffer.data(gp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xyyy_y = cbuffer.data(gp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xyyy_z = cbuffer.data(gp_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xyyz_x = cbuffer.data(gp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xyyz_y = cbuffer.data(gp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xyyz_z = cbuffer.data(gp_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xyzz_x = cbuffer.data(gp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xyzz_y = cbuffer.data(gp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xyzz_z = cbuffer.data(gp_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xzzz_x = cbuffer.data(gp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xzzz_y = cbuffer.data(gp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xzzz_z = cbuffer.data(gp_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_yyyy_x = cbuffer.data(gp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_yyyy_y = cbuffer.data(gp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_yyyy_z = cbuffer.data(gp_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_yyyz_x = cbuffer.data(gp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_yyyz_y = cbuffer.data(gp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_yyyz_z = cbuffer.data(gp_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_yyzz_x = cbuffer.data(gp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_yyzz_y = cbuffer.data(gp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_yyzz_z = cbuffer.data(gp_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_yzzz_x = cbuffer.data(gp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_yzzz_y = cbuffer.data(gp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_yzzz_z = cbuffer.data(gp_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_zzzz_x = cbuffer.data(gp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_zzzz_y = cbuffer.data(gp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_zzzz_z = cbuffer.data(gp_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_xxxx_x = cbuffer.data(gp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xxxx_y = cbuffer.data(gp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xxxx_z = cbuffer.data(gp_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_xxxy_x = cbuffer.data(gp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_xxxy_y = cbuffer.data(gp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_xxxy_z = cbuffer.data(gp_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_xxxz_x = cbuffer.data(gp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_xxxz_y = cbuffer.data(gp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_xxxz_z = cbuffer.data(gp_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_xxyy_x = cbuffer.data(gp_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_xxyy_y = cbuffer.data(gp_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_xxyy_z = cbuffer.data(gp_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_xxyz_x = cbuffer.data(gp_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_xxyz_y = cbuffer.data(gp_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_xxyz_z = cbuffer.data(gp_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_xxzz_x = cbuffer.data(gp_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_xxzz_y = cbuffer.data(gp_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_xxzz_z = cbuffer.data(gp_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_xyyy_x = cbuffer.data(gp_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_xyyy_y = cbuffer.data(gp_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_xyyy_z = cbuffer.data(gp_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_xyyz_x = cbuffer.data(gp_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_xyyz_y = cbuffer.data(gp_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_xyyz_z = cbuffer.data(gp_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_xyzz_x = cbuffer.data(gp_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_xyzz_y = cbuffer.data(gp_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_xyzz_z = cbuffer.data(gp_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_xzzz_x = cbuffer.data(gp_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_xzzz_y = cbuffer.data(gp_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_xzzz_z = cbuffer.data(gp_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_yyyy_x = cbuffer.data(gp_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_yyyy_y = cbuffer.data(gp_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_yyyy_z = cbuffer.data(gp_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_yyyz_x = cbuffer.data(gp_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_yyyz_y = cbuffer.data(gp_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_yyyz_z = cbuffer.data(gp_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_yyzz_x = cbuffer.data(gp_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_yyzz_y = cbuffer.data(gp_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_yyzz_z = cbuffer.data(gp_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_yzzz_x = cbuffer.data(gp_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_yzzz_y = cbuffer.data(gp_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_yzzz_z = cbuffer.data(gp_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_zzzz_x = cbuffer.data(gp_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_zzzz_y = cbuffer.data(gp_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_zzzz_z = cbuffer.data(gp_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_z_xxxx_x = cbuffer.data(gp_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_xxxx_y = cbuffer.data(gp_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_xxxx_z = cbuffer.data(gp_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_xxxy_x = cbuffer.data(gp_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_xxxy_y = cbuffer.data(gp_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_xxxy_z = cbuffer.data(gp_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_z_xxxz_x = cbuffer.data(gp_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_xxxz_y = cbuffer.data(gp_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_xxxz_z = cbuffer.data(gp_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_xxyy_x = cbuffer.data(gp_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_xxyy_y = cbuffer.data(gp_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_xxyy_z = cbuffer.data(gp_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_z_xxyz_x = cbuffer.data(gp_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_xxyz_y = cbuffer.data(gp_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_xxyz_z = cbuffer.data(gp_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_xxzz_x = cbuffer.data(gp_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_xxzz_y = cbuffer.data(gp_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_xxzz_z = cbuffer.data(gp_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_z_xyyy_x = cbuffer.data(gp_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_z_xyyy_y = cbuffer.data(gp_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_z_xyyy_z = cbuffer.data(gp_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_z_xyyz_x = cbuffer.data(gp_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_z_xyyz_y = cbuffer.data(gp_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_z_xyyz_z = cbuffer.data(gp_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_z_xyzz_x = cbuffer.data(gp_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_z_xyzz_y = cbuffer.data(gp_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_z_xyzz_z = cbuffer.data(gp_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_z_xzzz_x = cbuffer.data(gp_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_z_xzzz_y = cbuffer.data(gp_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_z_xzzz_z = cbuffer.data(gp_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_z_yyyy_x = cbuffer.data(gp_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_yyyy_y = cbuffer.data(gp_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_yyyy_z = cbuffer.data(gp_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_yyyz_x = cbuffer.data(gp_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_yyyz_y = cbuffer.data(gp_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_yyyz_z = cbuffer.data(gp_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_yyzz_x = cbuffer.data(gp_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_yyzz_y = cbuffer.data(gp_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_yyzz_z = cbuffer.data(gp_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_yzzz_x = cbuffer.data(gp_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_yzzz_y = cbuffer.data(gp_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_yzzz_z = cbuffer.data(gp_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_zzzz_x = cbuffer.data(gp_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_zzzz_y = cbuffer.data(gp_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_zzzz_z = cbuffer.data(gp_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_xxxx_x = cbuffer.data(gp_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_xxxx_y = cbuffer.data(gp_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_xxxx_z = cbuffer.data(gp_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_x_xxxy_x = cbuffer.data(gp_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_xxxy_y = cbuffer.data(gp_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_xxxy_z = cbuffer.data(gp_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_xxxz_x = cbuffer.data(gp_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_xxxz_y = cbuffer.data(gp_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_xxxz_z = cbuffer.data(gp_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_x_xxyy_x = cbuffer.data(gp_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_x_xxyy_y = cbuffer.data(gp_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_x_xxyy_z = cbuffer.data(gp_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_x_xxyz_x = cbuffer.data(gp_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_x_xxyz_y = cbuffer.data(gp_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_x_xxyz_z = cbuffer.data(gp_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_x_xxzz_x = cbuffer.data(gp_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_x_xxzz_y = cbuffer.data(gp_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_x_xxzz_z = cbuffer.data(gp_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_x_xyyy_x = cbuffer.data(gp_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_x_xyyy_y = cbuffer.data(gp_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_x_xyyy_z = cbuffer.data(gp_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_x_xyyz_x = cbuffer.data(gp_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_x_xyyz_y = cbuffer.data(gp_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_x_xyyz_z = cbuffer.data(gp_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_x_xyzz_x = cbuffer.data(gp_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_x_xyzz_y = cbuffer.data(gp_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_x_xyzz_z = cbuffer.data(gp_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_x_xzzz_x = cbuffer.data(gp_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_x_xzzz_y = cbuffer.data(gp_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_x_xzzz_z = cbuffer.data(gp_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_x_yyyy_x = cbuffer.data(gp_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_x_yyyy_y = cbuffer.data(gp_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_x_yyyy_z = cbuffer.data(gp_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_x_yyyz_x = cbuffer.data(gp_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_x_yyyz_y = cbuffer.data(gp_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_x_yyyz_z = cbuffer.data(gp_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_x_yyzz_x = cbuffer.data(gp_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_x_yyzz_y = cbuffer.data(gp_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_x_yyzz_z = cbuffer.data(gp_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_x_yzzz_x = cbuffer.data(gp_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_x_yzzz_y = cbuffer.data(gp_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_x_yzzz_z = cbuffer.data(gp_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_x_zzzz_x = cbuffer.data(gp_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_x_zzzz_y = cbuffer.data(gp_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_x_zzzz_z = cbuffer.data(gp_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_y_xxxx_x = cbuffer.data(gp_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_y_xxxx_y = cbuffer.data(gp_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_y_xxxx_z = cbuffer.data(gp_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_y_xxxy_x = cbuffer.data(gp_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_y_xxxy_y = cbuffer.data(gp_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_y_xxxy_z = cbuffer.data(gp_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_y_xxxz_x = cbuffer.data(gp_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_y_xxxz_y = cbuffer.data(gp_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_y_xxxz_z = cbuffer.data(gp_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_y_xxyy_x = cbuffer.data(gp_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_y_xxyy_y = cbuffer.data(gp_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_y_xxyy_z = cbuffer.data(gp_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_y_xxyz_x = cbuffer.data(gp_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_y_xxyz_y = cbuffer.data(gp_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_y_xxyz_z = cbuffer.data(gp_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_y_xxzz_x = cbuffer.data(gp_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_y_xxzz_y = cbuffer.data(gp_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_y_xxzz_z = cbuffer.data(gp_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_y_xyyy_x = cbuffer.data(gp_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_y_xyyy_y = cbuffer.data(gp_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_y_xyyy_z = cbuffer.data(gp_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_y_xyyz_x = cbuffer.data(gp_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_y_xyyz_y = cbuffer.data(gp_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_y_xyyz_z = cbuffer.data(gp_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_y_xyzz_x = cbuffer.data(gp_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_y_xyzz_y = cbuffer.data(gp_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_y_xyzz_z = cbuffer.data(gp_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_y_xzzz_x = cbuffer.data(gp_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_y_xzzz_y = cbuffer.data(gp_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_y_xzzz_z = cbuffer.data(gp_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_y_yyyy_x = cbuffer.data(gp_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_y_yyyy_y = cbuffer.data(gp_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_y_yyyy_z = cbuffer.data(gp_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_y_yyyz_x = cbuffer.data(gp_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_y_yyyz_y = cbuffer.data(gp_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_y_yyyz_z = cbuffer.data(gp_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_y_yyzz_x = cbuffer.data(gp_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_y_yyzz_y = cbuffer.data(gp_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_y_yyzz_z = cbuffer.data(gp_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_y_yzzz_x = cbuffer.data(gp_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_y_yzzz_y = cbuffer.data(gp_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_y_yzzz_z = cbuffer.data(gp_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_y_zzzz_x = cbuffer.data(gp_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_y_zzzz_y = cbuffer.data(gp_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_y_zzzz_z = cbuffer.data(gp_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_z_xxxx_x = cbuffer.data(gp_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_z_xxxx_y = cbuffer.data(gp_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_z_xxxx_z = cbuffer.data(gp_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_z_xxxy_x = cbuffer.data(gp_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_z_xxxy_y = cbuffer.data(gp_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_z_xxxy_z = cbuffer.data(gp_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_z_xxxz_x = cbuffer.data(gp_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_z_xxxz_y = cbuffer.data(gp_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_z_xxxz_z = cbuffer.data(gp_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_z_xxyy_x = cbuffer.data(gp_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_z_xxyy_y = cbuffer.data(gp_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_z_xxyy_z = cbuffer.data(gp_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_z_xxyz_x = cbuffer.data(gp_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_z_xxyz_y = cbuffer.data(gp_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_z_xxyz_z = cbuffer.data(gp_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_z_xxzz_x = cbuffer.data(gp_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_z_xxzz_y = cbuffer.data(gp_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_z_xxzz_z = cbuffer.data(gp_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_z_xyyy_x = cbuffer.data(gp_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_z_xyyy_y = cbuffer.data(gp_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_z_xyyy_z = cbuffer.data(gp_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_z_xyyz_x = cbuffer.data(gp_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_z_xyyz_y = cbuffer.data(gp_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_z_xyyz_z = cbuffer.data(gp_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_z_xyzz_x = cbuffer.data(gp_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_z_xyzz_y = cbuffer.data(gp_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_z_xyzz_z = cbuffer.data(gp_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_z_xzzz_x = cbuffer.data(gp_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_z_xzzz_y = cbuffer.data(gp_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_z_xzzz_z = cbuffer.data(gp_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_z_yyyy_x = cbuffer.data(gp_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_z_yyyy_y = cbuffer.data(gp_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_z_yyyy_z = cbuffer.data(gp_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_z_yyyz_x = cbuffer.data(gp_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_z_yyyz_y = cbuffer.data(gp_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_z_yyyz_z = cbuffer.data(gp_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_z_yyzz_x = cbuffer.data(gp_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_z_yyzz_y = cbuffer.data(gp_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_z_yyzz_z = cbuffer.data(gp_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_z_yzzz_x = cbuffer.data(gp_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_z_yzzz_y = cbuffer.data(gp_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_z_yzzz_z = cbuffer.data(gp_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_z_zzzz_x = cbuffer.data(gp_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_z_zzzz_y = cbuffer.data(gp_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_z_zzzz_z = cbuffer.data(gp_geom_11_off + 269 * ccomps * dcomps);

            auto g_z_x_xxxx_x = cbuffer.data(gp_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_x_xxxx_y = cbuffer.data(gp_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_x_xxxx_z = cbuffer.data(gp_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_x_xxxy_x = cbuffer.data(gp_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_x_xxxy_y = cbuffer.data(gp_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_x_xxxy_z = cbuffer.data(gp_geom_11_off + 275 * ccomps * dcomps);

            auto g_z_x_xxxz_x = cbuffer.data(gp_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_x_xxxz_y = cbuffer.data(gp_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_x_xxxz_z = cbuffer.data(gp_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_x_xxyy_x = cbuffer.data(gp_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_x_xxyy_y = cbuffer.data(gp_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_x_xxyy_z = cbuffer.data(gp_geom_11_off + 281 * ccomps * dcomps);

            auto g_z_x_xxyz_x = cbuffer.data(gp_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_x_xxyz_y = cbuffer.data(gp_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_x_xxyz_z = cbuffer.data(gp_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_x_xxzz_x = cbuffer.data(gp_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_x_xxzz_y = cbuffer.data(gp_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_x_xxzz_z = cbuffer.data(gp_geom_11_off + 287 * ccomps * dcomps);

            auto g_z_x_xyyy_x = cbuffer.data(gp_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_x_xyyy_y = cbuffer.data(gp_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_x_xyyy_z = cbuffer.data(gp_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_x_xyyz_x = cbuffer.data(gp_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_x_xyyz_y = cbuffer.data(gp_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_x_xyyz_z = cbuffer.data(gp_geom_11_off + 293 * ccomps * dcomps);

            auto g_z_x_xyzz_x = cbuffer.data(gp_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_x_xyzz_y = cbuffer.data(gp_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_x_xyzz_z = cbuffer.data(gp_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_x_xzzz_x = cbuffer.data(gp_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_x_xzzz_y = cbuffer.data(gp_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_x_xzzz_z = cbuffer.data(gp_geom_11_off + 299 * ccomps * dcomps);

            auto g_z_x_yyyy_x = cbuffer.data(gp_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_x_yyyy_y = cbuffer.data(gp_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_x_yyyy_z = cbuffer.data(gp_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_x_yyyz_x = cbuffer.data(gp_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_x_yyyz_y = cbuffer.data(gp_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_x_yyyz_z = cbuffer.data(gp_geom_11_off + 305 * ccomps * dcomps);

            auto g_z_x_yyzz_x = cbuffer.data(gp_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_x_yyzz_y = cbuffer.data(gp_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_x_yyzz_z = cbuffer.data(gp_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_x_yzzz_x = cbuffer.data(gp_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_x_yzzz_y = cbuffer.data(gp_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_x_yzzz_z = cbuffer.data(gp_geom_11_off + 311 * ccomps * dcomps);

            auto g_z_x_zzzz_x = cbuffer.data(gp_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_x_zzzz_y = cbuffer.data(gp_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_x_zzzz_z = cbuffer.data(gp_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_y_xxxx_x = cbuffer.data(gp_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_y_xxxx_y = cbuffer.data(gp_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_y_xxxx_z = cbuffer.data(gp_geom_11_off + 317 * ccomps * dcomps);

            auto g_z_y_xxxy_x = cbuffer.data(gp_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_y_xxxy_y = cbuffer.data(gp_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_y_xxxy_z = cbuffer.data(gp_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_y_xxxz_x = cbuffer.data(gp_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_y_xxxz_y = cbuffer.data(gp_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_y_xxxz_z = cbuffer.data(gp_geom_11_off + 323 * ccomps * dcomps);

            auto g_z_y_xxyy_x = cbuffer.data(gp_geom_11_off + 324 * ccomps * dcomps);

            auto g_z_y_xxyy_y = cbuffer.data(gp_geom_11_off + 325 * ccomps * dcomps);

            auto g_z_y_xxyy_z = cbuffer.data(gp_geom_11_off + 326 * ccomps * dcomps);

            auto g_z_y_xxyz_x = cbuffer.data(gp_geom_11_off + 327 * ccomps * dcomps);

            auto g_z_y_xxyz_y = cbuffer.data(gp_geom_11_off + 328 * ccomps * dcomps);

            auto g_z_y_xxyz_z = cbuffer.data(gp_geom_11_off + 329 * ccomps * dcomps);

            auto g_z_y_xxzz_x = cbuffer.data(gp_geom_11_off + 330 * ccomps * dcomps);

            auto g_z_y_xxzz_y = cbuffer.data(gp_geom_11_off + 331 * ccomps * dcomps);

            auto g_z_y_xxzz_z = cbuffer.data(gp_geom_11_off + 332 * ccomps * dcomps);

            auto g_z_y_xyyy_x = cbuffer.data(gp_geom_11_off + 333 * ccomps * dcomps);

            auto g_z_y_xyyy_y = cbuffer.data(gp_geom_11_off + 334 * ccomps * dcomps);

            auto g_z_y_xyyy_z = cbuffer.data(gp_geom_11_off + 335 * ccomps * dcomps);

            auto g_z_y_xyyz_x = cbuffer.data(gp_geom_11_off + 336 * ccomps * dcomps);

            auto g_z_y_xyyz_y = cbuffer.data(gp_geom_11_off + 337 * ccomps * dcomps);

            auto g_z_y_xyyz_z = cbuffer.data(gp_geom_11_off + 338 * ccomps * dcomps);

            auto g_z_y_xyzz_x = cbuffer.data(gp_geom_11_off + 339 * ccomps * dcomps);

            auto g_z_y_xyzz_y = cbuffer.data(gp_geom_11_off + 340 * ccomps * dcomps);

            auto g_z_y_xyzz_z = cbuffer.data(gp_geom_11_off + 341 * ccomps * dcomps);

            auto g_z_y_xzzz_x = cbuffer.data(gp_geom_11_off + 342 * ccomps * dcomps);

            auto g_z_y_xzzz_y = cbuffer.data(gp_geom_11_off + 343 * ccomps * dcomps);

            auto g_z_y_xzzz_z = cbuffer.data(gp_geom_11_off + 344 * ccomps * dcomps);

            auto g_z_y_yyyy_x = cbuffer.data(gp_geom_11_off + 345 * ccomps * dcomps);

            auto g_z_y_yyyy_y = cbuffer.data(gp_geom_11_off + 346 * ccomps * dcomps);

            auto g_z_y_yyyy_z = cbuffer.data(gp_geom_11_off + 347 * ccomps * dcomps);

            auto g_z_y_yyyz_x = cbuffer.data(gp_geom_11_off + 348 * ccomps * dcomps);

            auto g_z_y_yyyz_y = cbuffer.data(gp_geom_11_off + 349 * ccomps * dcomps);

            auto g_z_y_yyyz_z = cbuffer.data(gp_geom_11_off + 350 * ccomps * dcomps);

            auto g_z_y_yyzz_x = cbuffer.data(gp_geom_11_off + 351 * ccomps * dcomps);

            auto g_z_y_yyzz_y = cbuffer.data(gp_geom_11_off + 352 * ccomps * dcomps);

            auto g_z_y_yyzz_z = cbuffer.data(gp_geom_11_off + 353 * ccomps * dcomps);

            auto g_z_y_yzzz_x = cbuffer.data(gp_geom_11_off + 354 * ccomps * dcomps);

            auto g_z_y_yzzz_y = cbuffer.data(gp_geom_11_off + 355 * ccomps * dcomps);

            auto g_z_y_yzzz_z = cbuffer.data(gp_geom_11_off + 356 * ccomps * dcomps);

            auto g_z_y_zzzz_x = cbuffer.data(gp_geom_11_off + 357 * ccomps * dcomps);

            auto g_z_y_zzzz_y = cbuffer.data(gp_geom_11_off + 358 * ccomps * dcomps);

            auto g_z_y_zzzz_z = cbuffer.data(gp_geom_11_off + 359 * ccomps * dcomps);

            auto g_z_z_xxxx_x = cbuffer.data(gp_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_z_xxxx_y = cbuffer.data(gp_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_z_xxxx_z = cbuffer.data(gp_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_z_xxxy_x = cbuffer.data(gp_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_z_xxxy_y = cbuffer.data(gp_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_z_xxxy_z = cbuffer.data(gp_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_z_xxxz_x = cbuffer.data(gp_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_z_xxxz_y = cbuffer.data(gp_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_z_xxxz_z = cbuffer.data(gp_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_z_xxyy_x = cbuffer.data(gp_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_z_xxyy_y = cbuffer.data(gp_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_z_xxyy_z = cbuffer.data(gp_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_z_xxyz_x = cbuffer.data(gp_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_z_xxyz_y = cbuffer.data(gp_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_z_xxyz_z = cbuffer.data(gp_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_z_xxzz_x = cbuffer.data(gp_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_z_xxzz_y = cbuffer.data(gp_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_z_xxzz_z = cbuffer.data(gp_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_z_xyyy_x = cbuffer.data(gp_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_z_xyyy_y = cbuffer.data(gp_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_z_xyyy_z = cbuffer.data(gp_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_z_xyyz_x = cbuffer.data(gp_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_z_xyyz_y = cbuffer.data(gp_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_z_xyyz_z = cbuffer.data(gp_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_z_xyzz_x = cbuffer.data(gp_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_z_xyzz_y = cbuffer.data(gp_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_z_xyzz_z = cbuffer.data(gp_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_z_xzzz_x = cbuffer.data(gp_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_z_xzzz_y = cbuffer.data(gp_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_z_xzzz_z = cbuffer.data(gp_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_z_yyyy_x = cbuffer.data(gp_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_z_yyyy_y = cbuffer.data(gp_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_z_yyyy_z = cbuffer.data(gp_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_z_yyyz_x = cbuffer.data(gp_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_z_yyyz_y = cbuffer.data(gp_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_z_yyyz_z = cbuffer.data(gp_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_z_yyzz_x = cbuffer.data(gp_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_z_yyzz_y = cbuffer.data(gp_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_z_yyzz_z = cbuffer.data(gp_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_z_yzzz_x = cbuffer.data(gp_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_z_yzzz_y = cbuffer.data(gp_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_z_yzzz_z = cbuffer.data(gp_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_z_zzzz_x = cbuffer.data(gp_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_z_zzzz_y = cbuffer.data(gp_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_z_zzzz_z = cbuffer.data(gp_geom_11_off + 404 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hsxx

            const auto hs_geom_11_off = idx_geom_11_hsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxx_0 = cbuffer.data(hs_geom_11_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_0, g_x_0_xxxx_0, g_x_x_xxxx_0, g_x_x_xxxx_x, g_x_x_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxx_0[k] = -g_0_x_xxxx_0[k] + g_x_0_xxxx_0[k] - g_x_x_xxxx_0[k] * ab_x + g_x_x_xxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxy_0 = cbuffer.data(hs_geom_11_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxx_0, g_x_x_xxxx_y, g_x_x_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxy_0[k] = -g_x_x_xxxx_0[k] * ab_y + g_x_x_xxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxz_0 = cbuffer.data(hs_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxx_0, g_x_x_xxxx_z, g_x_x_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxz_0[k] = -g_x_x_xxxx_0[k] * ab_z + g_x_x_xxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxyy_0 = cbuffer.data(hs_geom_11_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxy_0, g_x_x_xxxy_y, g_x_x_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxyy_0[k] = -g_x_x_xxxy_0[k] * ab_y + g_x_x_xxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxyz_0 = cbuffer.data(hs_geom_11_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxyz_0, g_x_x_xxxz_0, g_x_x_xxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxyz_0[k] = -g_x_x_xxxz_0[k] * ab_y + g_x_x_xxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxzz_0 = cbuffer.data(hs_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxz_0, g_x_x_xxxz_z, g_x_x_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxzz_0[k] = -g_x_x_xxxz_0[k] * ab_z + g_x_x_xxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyyy_0 = cbuffer.data(hs_geom_11_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyy_0, g_x_x_xxyy_y, g_x_x_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyyy_0[k] = -g_x_x_xxyy_0[k] * ab_y + g_x_x_xxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyyz_0 = cbuffer.data(hs_geom_11_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyyz_0, g_x_x_xxyz_0, g_x_x_xxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyyz_0[k] = -g_x_x_xxyz_0[k] * ab_y + g_x_x_xxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyzz_0 = cbuffer.data(hs_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyzz_0, g_x_x_xxzz_0, g_x_x_xxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyzz_0[k] = -g_x_x_xxzz_0[k] * ab_y + g_x_x_xxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxzzz_0 = cbuffer.data(hs_geom_11_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxzz_0, g_x_x_xxzz_z, g_x_x_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxzzz_0[k] = -g_x_x_xxzz_0[k] * ab_z + g_x_x_xxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyyy_0 = cbuffer.data(hs_geom_11_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyy_0, g_x_x_xyyy_y, g_x_x_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyyy_0[k] = -g_x_x_xyyy_0[k] * ab_y + g_x_x_xyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyyz_0 = cbuffer.data(hs_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyyz_0, g_x_x_xyyz_0, g_x_x_xyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyyz_0[k] = -g_x_x_xyyz_0[k] * ab_y + g_x_x_xyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyzz_0 = cbuffer.data(hs_geom_11_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyzz_0, g_x_x_xyzz_0, g_x_x_xyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyzz_0[k] = -g_x_x_xyzz_0[k] * ab_y + g_x_x_xyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyzzz_0 = cbuffer.data(hs_geom_11_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyzzz_0, g_x_x_xzzz_0, g_x_x_xzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyzzz_0[k] = -g_x_x_xzzz_0[k] * ab_y + g_x_x_xzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzzzz_0 = cbuffer.data(hs_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xzzz_0, g_x_x_xzzz_z, g_x_x_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzzzz_0[k] = -g_x_x_xzzz_0[k] * ab_z + g_x_x_xzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyyy_0 = cbuffer.data(hs_geom_11_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyy_0, g_x_x_yyyy_y, g_x_x_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyyy_0[k] = -g_x_x_yyyy_0[k] * ab_y + g_x_x_yyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyyz_0 = cbuffer.data(hs_geom_11_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyyz_0, g_x_x_yyyz_0, g_x_x_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyyz_0[k] = -g_x_x_yyyz_0[k] * ab_y + g_x_x_yyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyzz_0 = cbuffer.data(hs_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyzz_0, g_x_x_yyzz_0, g_x_x_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyzz_0[k] = -g_x_x_yyzz_0[k] * ab_y + g_x_x_yyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyzzz_0 = cbuffer.data(hs_geom_11_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyzzz_0, g_x_x_yzzz_0, g_x_x_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyzzz_0[k] = -g_x_x_yzzz_0[k] * ab_y + g_x_x_yzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzzzz_0 = cbuffer.data(hs_geom_11_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzzzz_0, g_x_x_zzzz_0, g_x_x_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzzzz_0[k] = -g_x_x_zzzz_0[k] * ab_y + g_x_x_zzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzzzz_0 = cbuffer.data(hs_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zzzz_0, g_x_x_zzzz_z, g_x_x_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzzzz_0[k] = -g_x_x_zzzz_0[k] * ab_z + g_x_x_zzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxx_0 = cbuffer.data(hs_geom_11_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_0, g_x_y_xxxx_0, g_x_y_xxxx_x, g_x_y_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxx_0[k] = -g_0_y_xxxx_0[k] - g_x_y_xxxx_0[k] * ab_x + g_x_y_xxxx_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxy_0 = cbuffer.data(hs_geom_11_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_0, g_x_y_xxxxy_0, g_x_y_xxxy_0, g_x_y_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxy_0[k] = -g_0_y_xxxy_0[k] - g_x_y_xxxy_0[k] * ab_x + g_x_y_xxxy_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxz_0 = cbuffer.data(hs_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxx_0, g_x_y_xxxx_z, g_x_y_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxz_0[k] = -g_x_y_xxxx_0[k] * ab_z + g_x_y_xxxx_z[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxyy_0 = cbuffer.data(hs_geom_11_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_0, g_x_y_xxxyy_0, g_x_y_xxyy_0, g_x_y_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxyy_0[k] = -g_0_y_xxyy_0[k] - g_x_y_xxyy_0[k] * ab_x + g_x_y_xxyy_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxyz_0 = cbuffer.data(hs_geom_11_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxy_0, g_x_y_xxxy_z, g_x_y_xxxyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxyz_0[k] = -g_x_y_xxxy_0[k] * ab_z + g_x_y_xxxy_z[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxzz_0 = cbuffer.data(hs_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxz_0, g_x_y_xxxz_z, g_x_y_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxzz_0[k] = -g_x_y_xxxz_0[k] * ab_z + g_x_y_xxxz_z[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyyy_0 = cbuffer.data(hs_geom_11_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_0, g_x_y_xxyyy_0, g_x_y_xyyy_0, g_x_y_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyyy_0[k] = -g_0_y_xyyy_0[k] - g_x_y_xyyy_0[k] * ab_x + g_x_y_xyyy_x[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyyz_0 = cbuffer.data(hs_geom_11_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxyy_0, g_x_y_xxyy_z, g_x_y_xxyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyyz_0[k] = -g_x_y_xxyy_0[k] * ab_z + g_x_y_xxyy_z[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyzz_0 = cbuffer.data(hs_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxyz_0, g_x_y_xxyz_z, g_x_y_xxyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyzz_0[k] = -g_x_y_xxyz_0[k] * ab_z + g_x_y_xxyz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxzzz_0 = cbuffer.data(hs_geom_11_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxzz_0, g_x_y_xxzz_z, g_x_y_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxzzz_0[k] = -g_x_y_xxzz_0[k] * ab_z + g_x_y_xxzz_z[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyyy_0 = cbuffer.data(hs_geom_11_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_0, g_x_y_xyyyy_0, g_x_y_yyyy_0, g_x_y_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyyy_0[k] = -g_0_y_yyyy_0[k] - g_x_y_yyyy_0[k] * ab_x + g_x_y_yyyy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyyz_0 = cbuffer.data(hs_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyyy_0, g_x_y_xyyy_z, g_x_y_xyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyyz_0[k] = -g_x_y_xyyy_0[k] * ab_z + g_x_y_xyyy_z[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyzz_0 = cbuffer.data(hs_geom_11_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyyz_0, g_x_y_xyyz_z, g_x_y_xyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyzz_0[k] = -g_x_y_xyyz_0[k] * ab_z + g_x_y_xyyz_z[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyzzz_0 = cbuffer.data(hs_geom_11_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyzz_0, g_x_y_xyzz_z, g_x_y_xyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyzzz_0[k] = -g_x_y_xyzz_0[k] * ab_z + g_x_y_xyzz_z[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzzzz_0 = cbuffer.data(hs_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xzzz_0, g_x_y_xzzz_z, g_x_y_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzzzz_0[k] = -g_x_y_xzzz_0[k] * ab_z + g_x_y_xzzz_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyyy_0 = cbuffer.data(hs_geom_11_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_0, g_x_y_yyyy_0, g_x_y_yyyy_y, g_x_y_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyyy_0[k] = g_x_0_yyyy_0[k] - g_x_y_yyyy_0[k] * ab_y + g_x_y_yyyy_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyyz_0 = cbuffer.data(hs_geom_11_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyyy_0, g_x_y_yyyy_z, g_x_y_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyyz_0[k] = -g_x_y_yyyy_0[k] * ab_z + g_x_y_yyyy_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyzz_0 = cbuffer.data(hs_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyyz_0, g_x_y_yyyz_z, g_x_y_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyzz_0[k] = -g_x_y_yyyz_0[k] * ab_z + g_x_y_yyyz_z[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyzzz_0 = cbuffer.data(hs_geom_11_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyzz_0, g_x_y_yyzz_z, g_x_y_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyzzz_0[k] = -g_x_y_yyzz_0[k] * ab_z + g_x_y_yyzz_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzzzz_0 = cbuffer.data(hs_geom_11_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yzzz_0, g_x_y_yzzz_z, g_x_y_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzzzz_0[k] = -g_x_y_yzzz_0[k] * ab_z + g_x_y_yzzz_z[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzzzz_0 = cbuffer.data(hs_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zzzz_0, g_x_y_zzzz_z, g_x_y_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzzzz_0[k] = -g_x_y_zzzz_0[k] * ab_z + g_x_y_zzzz_z[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxx_0 = cbuffer.data(hs_geom_11_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_0, g_x_z_xxxx_0, g_x_z_xxxx_x, g_x_z_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxx_0[k] = -g_0_z_xxxx_0[k] - g_x_z_xxxx_0[k] * ab_x + g_x_z_xxxx_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxy_0 = cbuffer.data(hs_geom_11_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxx_0, g_x_z_xxxx_y, g_x_z_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxy_0[k] = -g_x_z_xxxx_0[k] * ab_y + g_x_z_xxxx_y[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxz_0 = cbuffer.data(hs_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_0, g_x_z_xxxxz_0, g_x_z_xxxz_0, g_x_z_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxz_0[k] = -g_0_z_xxxz_0[k] - g_x_z_xxxz_0[k] * ab_x + g_x_z_xxxz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxyy_0 = cbuffer.data(hs_geom_11_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxy_0, g_x_z_xxxy_y, g_x_z_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxyy_0[k] = -g_x_z_xxxy_0[k] * ab_y + g_x_z_xxxy_y[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxyz_0 = cbuffer.data(hs_geom_11_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxyz_0, g_x_z_xxxz_0, g_x_z_xxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxyz_0[k] = -g_x_z_xxxz_0[k] * ab_y + g_x_z_xxxz_y[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxzz_0 = cbuffer.data(hs_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_0, g_x_z_xxxzz_0, g_x_z_xxzz_0, g_x_z_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxzz_0[k] = -g_0_z_xxzz_0[k] - g_x_z_xxzz_0[k] * ab_x + g_x_z_xxzz_x[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyyy_0 = cbuffer.data(hs_geom_11_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyy_0, g_x_z_xxyy_y, g_x_z_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyyy_0[k] = -g_x_z_xxyy_0[k] * ab_y + g_x_z_xxyy_y[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyyz_0 = cbuffer.data(hs_geom_11_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyyz_0, g_x_z_xxyz_0, g_x_z_xxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyyz_0[k] = -g_x_z_xxyz_0[k] * ab_y + g_x_z_xxyz_y[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyzz_0 = cbuffer.data(hs_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyzz_0, g_x_z_xxzz_0, g_x_z_xxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyzz_0[k] = -g_x_z_xxzz_0[k] * ab_y + g_x_z_xxzz_y[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxzzz_0 = cbuffer.data(hs_geom_11_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_0, g_x_z_xxzzz_0, g_x_z_xzzz_0, g_x_z_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxzzz_0[k] = -g_0_z_xzzz_0[k] - g_x_z_xzzz_0[k] * ab_x + g_x_z_xzzz_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyyy_0 = cbuffer.data(hs_geom_11_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyy_0, g_x_z_xyyy_y, g_x_z_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyyy_0[k] = -g_x_z_xyyy_0[k] * ab_y + g_x_z_xyyy_y[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyyz_0 = cbuffer.data(hs_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyyz_0, g_x_z_xyyz_0, g_x_z_xyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyyz_0[k] = -g_x_z_xyyz_0[k] * ab_y + g_x_z_xyyz_y[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyzz_0 = cbuffer.data(hs_geom_11_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyzz_0, g_x_z_xyzz_0, g_x_z_xyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyzz_0[k] = -g_x_z_xyzz_0[k] * ab_y + g_x_z_xyzz_y[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyzzz_0 = cbuffer.data(hs_geom_11_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyzzz_0, g_x_z_xzzz_0, g_x_z_xzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyzzz_0[k] = -g_x_z_xzzz_0[k] * ab_y + g_x_z_xzzz_y[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzzzz_0 = cbuffer.data(hs_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_0, g_x_z_xzzzz_0, g_x_z_zzzz_0, g_x_z_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzzzz_0[k] = -g_0_z_zzzz_0[k] - g_x_z_zzzz_0[k] * ab_x + g_x_z_zzzz_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyyy_0 = cbuffer.data(hs_geom_11_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyy_0, g_x_z_yyyy_y, g_x_z_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyyy_0[k] = -g_x_z_yyyy_0[k] * ab_y + g_x_z_yyyy_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyyz_0 = cbuffer.data(hs_geom_11_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyyz_0, g_x_z_yyyz_0, g_x_z_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyyz_0[k] = -g_x_z_yyyz_0[k] * ab_y + g_x_z_yyyz_y[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyzz_0 = cbuffer.data(hs_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyzz_0, g_x_z_yyzz_0, g_x_z_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyzz_0[k] = -g_x_z_yyzz_0[k] * ab_y + g_x_z_yyzz_y[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyzzz_0 = cbuffer.data(hs_geom_11_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyzzz_0, g_x_z_yzzz_0, g_x_z_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyzzz_0[k] = -g_x_z_yzzz_0[k] * ab_y + g_x_z_yzzz_y[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzzzz_0 = cbuffer.data(hs_geom_11_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzzzz_0, g_x_z_zzzz_0, g_x_z_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzzzz_0[k] = -g_x_z_zzzz_0[k] * ab_y + g_x_z_zzzz_y[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzzzz_0 = cbuffer.data(hs_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_0, g_x_z_zzzz_0, g_x_z_zzzz_z, g_x_z_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzzzz_0[k] = g_x_0_zzzz_0[k] - g_x_z_zzzz_0[k] * ab_z + g_x_z_zzzz_z[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxx_0 = cbuffer.data(hs_geom_11_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxx_0, g_y_x_xxxx_0, g_y_x_xxxx_x, g_y_x_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxx_0[k] = g_y_0_xxxx_0[k] - g_y_x_xxxx_0[k] * ab_x + g_y_x_xxxx_x[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxy_0 = cbuffer.data(hs_geom_11_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_0, g_y_x_xxxxy_0, g_y_x_xxxy_0, g_y_x_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxy_0[k] = g_y_0_xxxy_0[k] - g_y_x_xxxy_0[k] * ab_x + g_y_x_xxxy_x[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxz_0 = cbuffer.data(hs_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxx_0, g_y_x_xxxx_z, g_y_x_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxz_0[k] = -g_y_x_xxxx_0[k] * ab_z + g_y_x_xxxx_z[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxyy_0 = cbuffer.data(hs_geom_11_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_0, g_y_x_xxxyy_0, g_y_x_xxyy_0, g_y_x_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxyy_0[k] = g_y_0_xxyy_0[k] - g_y_x_xxyy_0[k] * ab_x + g_y_x_xxyy_x[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxyz_0 = cbuffer.data(hs_geom_11_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxy_0, g_y_x_xxxy_z, g_y_x_xxxyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxyz_0[k] = -g_y_x_xxxy_0[k] * ab_z + g_y_x_xxxy_z[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxzz_0 = cbuffer.data(hs_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxz_0, g_y_x_xxxz_z, g_y_x_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxzz_0[k] = -g_y_x_xxxz_0[k] * ab_z + g_y_x_xxxz_z[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyyy_0 = cbuffer.data(hs_geom_11_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_0, g_y_x_xxyyy_0, g_y_x_xyyy_0, g_y_x_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyyy_0[k] = g_y_0_xyyy_0[k] - g_y_x_xyyy_0[k] * ab_x + g_y_x_xyyy_x[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyyz_0 = cbuffer.data(hs_geom_11_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxyy_0, g_y_x_xxyy_z, g_y_x_xxyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyyz_0[k] = -g_y_x_xxyy_0[k] * ab_z + g_y_x_xxyy_z[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyzz_0 = cbuffer.data(hs_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxyz_0, g_y_x_xxyz_z, g_y_x_xxyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyzz_0[k] = -g_y_x_xxyz_0[k] * ab_z + g_y_x_xxyz_z[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxzzz_0 = cbuffer.data(hs_geom_11_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxzz_0, g_y_x_xxzz_z, g_y_x_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxzzz_0[k] = -g_y_x_xxzz_0[k] * ab_z + g_y_x_xxzz_z[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyyy_0 = cbuffer.data(hs_geom_11_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_0, g_y_x_xyyyy_0, g_y_x_yyyy_0, g_y_x_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyyy_0[k] = g_y_0_yyyy_0[k] - g_y_x_yyyy_0[k] * ab_x + g_y_x_yyyy_x[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyyz_0 = cbuffer.data(hs_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyyy_0, g_y_x_xyyy_z, g_y_x_xyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyyz_0[k] = -g_y_x_xyyy_0[k] * ab_z + g_y_x_xyyy_z[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyzz_0 = cbuffer.data(hs_geom_11_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyyz_0, g_y_x_xyyz_z, g_y_x_xyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyzz_0[k] = -g_y_x_xyyz_0[k] * ab_z + g_y_x_xyyz_z[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyzzz_0 = cbuffer.data(hs_geom_11_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyzz_0, g_y_x_xyzz_z, g_y_x_xyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyzzz_0[k] = -g_y_x_xyzz_0[k] * ab_z + g_y_x_xyzz_z[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzzzz_0 = cbuffer.data(hs_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xzzz_0, g_y_x_xzzz_z, g_y_x_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzzzz_0[k] = -g_y_x_xzzz_0[k] * ab_z + g_y_x_xzzz_z[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyyy_0 = cbuffer.data(hs_geom_11_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_0, g_y_x_yyyy_0, g_y_x_yyyy_y, g_y_x_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyyy_0[k] = -g_0_x_yyyy_0[k] - g_y_x_yyyy_0[k] * ab_y + g_y_x_yyyy_y[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyyz_0 = cbuffer.data(hs_geom_11_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyyy_0, g_y_x_yyyy_z, g_y_x_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyyz_0[k] = -g_y_x_yyyy_0[k] * ab_z + g_y_x_yyyy_z[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyzz_0 = cbuffer.data(hs_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyyz_0, g_y_x_yyyz_z, g_y_x_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyzz_0[k] = -g_y_x_yyyz_0[k] * ab_z + g_y_x_yyyz_z[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyzzz_0 = cbuffer.data(hs_geom_11_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyzz_0, g_y_x_yyzz_z, g_y_x_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyzzz_0[k] = -g_y_x_yyzz_0[k] * ab_z + g_y_x_yyzz_z[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzzzz_0 = cbuffer.data(hs_geom_11_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yzzz_0, g_y_x_yzzz_z, g_y_x_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzzzz_0[k] = -g_y_x_yzzz_0[k] * ab_z + g_y_x_yzzz_z[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzzzz_0 = cbuffer.data(hs_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zzzz_0, g_y_x_zzzz_z, g_y_x_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzzzz_0[k] = -g_y_x_zzzz_0[k] * ab_z + g_y_x_zzzz_z[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxx_0 = cbuffer.data(hs_geom_11_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxx_0, g_y_y_xxxx_x, g_y_y_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxx_0[k] = -g_y_y_xxxx_0[k] * ab_x + g_y_y_xxxx_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxy_0 = cbuffer.data(hs_geom_11_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxxy_0, g_y_y_xxxy_0, g_y_y_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxy_0[k] = -g_y_y_xxxy_0[k] * ab_x + g_y_y_xxxy_x[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxz_0 = cbuffer.data(hs_geom_11_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxxz_0, g_y_y_xxxz_0, g_y_y_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxz_0[k] = -g_y_y_xxxz_0[k] * ab_x + g_y_y_xxxz_x[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxyy_0 = cbuffer.data(hs_geom_11_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxyy_0, g_y_y_xxyy_0, g_y_y_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxyy_0[k] = -g_y_y_xxyy_0[k] * ab_x + g_y_y_xxyy_x[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxyz_0 = cbuffer.data(hs_geom_11_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxyz_0, g_y_y_xxyz_0, g_y_y_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxyz_0[k] = -g_y_y_xxyz_0[k] * ab_x + g_y_y_xxyz_x[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxzz_0 = cbuffer.data(hs_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxzz_0, g_y_y_xxzz_0, g_y_y_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxzz_0[k] = -g_y_y_xxzz_0[k] * ab_x + g_y_y_xxzz_x[k];
            }

            /// Set up 90-91 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyyy_0 = cbuffer.data(hs_geom_11_off + 90 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyyy_0, g_y_y_xyyy_0, g_y_y_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyyy_0[k] = -g_y_y_xyyy_0[k] * ab_x + g_y_y_xyyy_x[k];
            }

            /// Set up 91-92 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyyz_0 = cbuffer.data(hs_geom_11_off + 91 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyyz_0, g_y_y_xyyz_0, g_y_y_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyyz_0[k] = -g_y_y_xyyz_0[k] * ab_x + g_y_y_xyyz_x[k];
            }

            /// Set up 92-93 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyzz_0 = cbuffer.data(hs_geom_11_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyzz_0, g_y_y_xyzz_0, g_y_y_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyzz_0[k] = -g_y_y_xyzz_0[k] * ab_x + g_y_y_xyzz_x[k];
            }

            /// Set up 93-94 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxzzz_0 = cbuffer.data(hs_geom_11_off + 93 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxzzz_0, g_y_y_xzzz_0, g_y_y_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxzzz_0[k] = -g_y_y_xzzz_0[k] * ab_x + g_y_y_xzzz_x[k];
            }

            /// Set up 94-95 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyyy_0 = cbuffer.data(hs_geom_11_off + 94 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyyy_0, g_y_y_yyyy_0, g_y_y_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyyy_0[k] = -g_y_y_yyyy_0[k] * ab_x + g_y_y_yyyy_x[k];
            }

            /// Set up 95-96 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyyz_0 = cbuffer.data(hs_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyyz_0, g_y_y_yyyz_0, g_y_y_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyyz_0[k] = -g_y_y_yyyz_0[k] * ab_x + g_y_y_yyyz_x[k];
            }

            /// Set up 96-97 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyzz_0 = cbuffer.data(hs_geom_11_off + 96 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyzz_0, g_y_y_yyzz_0, g_y_y_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyzz_0[k] = -g_y_y_yyzz_0[k] * ab_x + g_y_y_yyzz_x[k];
            }

            /// Set up 97-98 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyzzz_0 = cbuffer.data(hs_geom_11_off + 97 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyzzz_0, g_y_y_yzzz_0, g_y_y_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyzzz_0[k] = -g_y_y_yzzz_0[k] * ab_x + g_y_y_yzzz_x[k];
            }

            /// Set up 98-99 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzzzz_0 = cbuffer.data(hs_geom_11_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzzzz_0, g_y_y_zzzz_0, g_y_y_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzzzz_0[k] = -g_y_y_zzzz_0[k] * ab_x + g_y_y_zzzz_x[k];
            }

            /// Set up 99-100 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyyy_0 = cbuffer.data(hs_geom_11_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_0, g_y_0_yyyy_0, g_y_y_yyyy_0, g_y_y_yyyy_y, g_y_y_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyyy_0[k] = -g_0_y_yyyy_0[k] + g_y_0_yyyy_0[k] - g_y_y_yyyy_0[k] * ab_y + g_y_y_yyyy_y[k];
            }

            /// Set up 100-101 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyyz_0 = cbuffer.data(hs_geom_11_off + 100 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyyy_0, g_y_y_yyyy_z, g_y_y_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyyz_0[k] = -g_y_y_yyyy_0[k] * ab_z + g_y_y_yyyy_z[k];
            }

            /// Set up 101-102 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyzz_0 = cbuffer.data(hs_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyyz_0, g_y_y_yyyz_z, g_y_y_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyzz_0[k] = -g_y_y_yyyz_0[k] * ab_z + g_y_y_yyyz_z[k];
            }

            /// Set up 102-103 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyzzz_0 = cbuffer.data(hs_geom_11_off + 102 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyzz_0, g_y_y_yyzz_z, g_y_y_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyzzz_0[k] = -g_y_y_yyzz_0[k] * ab_z + g_y_y_yyzz_z[k];
            }

            /// Set up 103-104 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzzzz_0 = cbuffer.data(hs_geom_11_off + 103 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yzzz_0, g_y_y_yzzz_z, g_y_y_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzzzz_0[k] = -g_y_y_yzzz_0[k] * ab_z + g_y_y_yzzz_z[k];
            }

            /// Set up 104-105 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzzzz_0 = cbuffer.data(hs_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zzzz_0, g_y_y_zzzz_z, g_y_y_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzzzz_0[k] = -g_y_y_zzzz_0[k] * ab_z + g_y_y_zzzz_z[k];
            }

            /// Set up 105-106 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxx_0 = cbuffer.data(hs_geom_11_off + 105 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxx_0, g_y_z_xxxx_x, g_y_z_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxx_0[k] = -g_y_z_xxxx_0[k] * ab_x + g_y_z_xxxx_x[k];
            }

            /// Set up 106-107 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxy_0 = cbuffer.data(hs_geom_11_off + 106 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxxy_0, g_y_z_xxxy_0, g_y_z_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxy_0[k] = -g_y_z_xxxy_0[k] * ab_x + g_y_z_xxxy_x[k];
            }

            /// Set up 107-108 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxz_0 = cbuffer.data(hs_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxxz_0, g_y_z_xxxz_0, g_y_z_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxz_0[k] = -g_y_z_xxxz_0[k] * ab_x + g_y_z_xxxz_x[k];
            }

            /// Set up 108-109 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxyy_0 = cbuffer.data(hs_geom_11_off + 108 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxyy_0, g_y_z_xxyy_0, g_y_z_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxyy_0[k] = -g_y_z_xxyy_0[k] * ab_x + g_y_z_xxyy_x[k];
            }

            /// Set up 109-110 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxyz_0 = cbuffer.data(hs_geom_11_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxyz_0, g_y_z_xxyz_0, g_y_z_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxyz_0[k] = -g_y_z_xxyz_0[k] * ab_x + g_y_z_xxyz_x[k];
            }

            /// Set up 110-111 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxzz_0 = cbuffer.data(hs_geom_11_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxzz_0, g_y_z_xxzz_0, g_y_z_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxzz_0[k] = -g_y_z_xxzz_0[k] * ab_x + g_y_z_xxzz_x[k];
            }

            /// Set up 111-112 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyyy_0 = cbuffer.data(hs_geom_11_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyyy_0, g_y_z_xyyy_0, g_y_z_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyyy_0[k] = -g_y_z_xyyy_0[k] * ab_x + g_y_z_xyyy_x[k];
            }

            /// Set up 112-113 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyyz_0 = cbuffer.data(hs_geom_11_off + 112 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyyz_0, g_y_z_xyyz_0, g_y_z_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyyz_0[k] = -g_y_z_xyyz_0[k] * ab_x + g_y_z_xyyz_x[k];
            }

            /// Set up 113-114 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyzz_0 = cbuffer.data(hs_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyzz_0, g_y_z_xyzz_0, g_y_z_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyzz_0[k] = -g_y_z_xyzz_0[k] * ab_x + g_y_z_xyzz_x[k];
            }

            /// Set up 114-115 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxzzz_0 = cbuffer.data(hs_geom_11_off + 114 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxzzz_0, g_y_z_xzzz_0, g_y_z_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxzzz_0[k] = -g_y_z_xzzz_0[k] * ab_x + g_y_z_xzzz_x[k];
            }

            /// Set up 115-116 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyyy_0 = cbuffer.data(hs_geom_11_off + 115 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyyy_0, g_y_z_yyyy_0, g_y_z_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyyy_0[k] = -g_y_z_yyyy_0[k] * ab_x + g_y_z_yyyy_x[k];
            }

            /// Set up 116-117 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyyz_0 = cbuffer.data(hs_geom_11_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyyz_0, g_y_z_yyyz_0, g_y_z_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyyz_0[k] = -g_y_z_yyyz_0[k] * ab_x + g_y_z_yyyz_x[k];
            }

            /// Set up 117-118 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyzz_0 = cbuffer.data(hs_geom_11_off + 117 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyzz_0, g_y_z_yyzz_0, g_y_z_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyzz_0[k] = -g_y_z_yyzz_0[k] * ab_x + g_y_z_yyzz_x[k];
            }

            /// Set up 118-119 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyzzz_0 = cbuffer.data(hs_geom_11_off + 118 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyzzz_0, g_y_z_yzzz_0, g_y_z_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyzzz_0[k] = -g_y_z_yzzz_0[k] * ab_x + g_y_z_yzzz_x[k];
            }

            /// Set up 119-120 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzzzz_0 = cbuffer.data(hs_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzzzz_0, g_y_z_zzzz_0, g_y_z_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzzzz_0[k] = -g_y_z_zzzz_0[k] * ab_x + g_y_z_zzzz_x[k];
            }

            /// Set up 120-121 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyyy_0 = cbuffer.data(hs_geom_11_off + 120 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_0, g_y_z_yyyy_0, g_y_z_yyyy_y, g_y_z_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyyy_0[k] = -g_0_z_yyyy_0[k] - g_y_z_yyyy_0[k] * ab_y + g_y_z_yyyy_y[k];
            }

            /// Set up 121-122 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyyz_0 = cbuffer.data(hs_geom_11_off + 121 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_0, g_y_z_yyyyz_0, g_y_z_yyyz_0, g_y_z_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyyz_0[k] = -g_0_z_yyyz_0[k] - g_y_z_yyyz_0[k] * ab_y + g_y_z_yyyz_y[k];
            }

            /// Set up 122-123 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyzz_0 = cbuffer.data(hs_geom_11_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_0, g_y_z_yyyzz_0, g_y_z_yyzz_0, g_y_z_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyzz_0[k] = -g_0_z_yyzz_0[k] - g_y_z_yyzz_0[k] * ab_y + g_y_z_yyzz_y[k];
            }

            /// Set up 123-124 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyzzz_0 = cbuffer.data(hs_geom_11_off + 123 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_0, g_y_z_yyzzz_0, g_y_z_yzzz_0, g_y_z_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyzzz_0[k] = -g_0_z_yzzz_0[k] - g_y_z_yzzz_0[k] * ab_y + g_y_z_yzzz_y[k];
            }

            /// Set up 124-125 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzzzz_0 = cbuffer.data(hs_geom_11_off + 124 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_0, g_y_z_yzzzz_0, g_y_z_zzzz_0, g_y_z_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzzzz_0[k] = -g_0_z_zzzz_0[k] - g_y_z_zzzz_0[k] * ab_y + g_y_z_zzzz_y[k];
            }

            /// Set up 125-126 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzzzz_0 = cbuffer.data(hs_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_0, g_y_z_zzzz_0, g_y_z_zzzz_z, g_y_z_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzzzz_0[k] = g_y_0_zzzz_0[k] - g_y_z_zzzz_0[k] * ab_z + g_y_z_zzzz_z[k];
            }

            /// Set up 126-127 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxx_0 = cbuffer.data(hs_geom_11_off + 126 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxx_0, g_z_x_xxxx_0, g_z_x_xxxx_x, g_z_x_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxx_0[k] = g_z_0_xxxx_0[k] - g_z_x_xxxx_0[k] * ab_x + g_z_x_xxxx_x[k];
            }

            /// Set up 127-128 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxy_0 = cbuffer.data(hs_geom_11_off + 127 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxx_0, g_z_x_xxxx_y, g_z_x_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxy_0[k] = -g_z_x_xxxx_0[k] * ab_y + g_z_x_xxxx_y[k];
            }

            /// Set up 128-129 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxz_0 = cbuffer.data(hs_geom_11_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_0, g_z_x_xxxxz_0, g_z_x_xxxz_0, g_z_x_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxz_0[k] = g_z_0_xxxz_0[k] - g_z_x_xxxz_0[k] * ab_x + g_z_x_xxxz_x[k];
            }

            /// Set up 129-130 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxyy_0 = cbuffer.data(hs_geom_11_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxy_0, g_z_x_xxxy_y, g_z_x_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxyy_0[k] = -g_z_x_xxxy_0[k] * ab_y + g_z_x_xxxy_y[k];
            }

            /// Set up 130-131 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxyz_0 = cbuffer.data(hs_geom_11_off + 130 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxyz_0, g_z_x_xxxz_0, g_z_x_xxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxyz_0[k] = -g_z_x_xxxz_0[k] * ab_y + g_z_x_xxxz_y[k];
            }

            /// Set up 131-132 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxzz_0 = cbuffer.data(hs_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_0, g_z_x_xxxzz_0, g_z_x_xxzz_0, g_z_x_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxzz_0[k] = g_z_0_xxzz_0[k] - g_z_x_xxzz_0[k] * ab_x + g_z_x_xxzz_x[k];
            }

            /// Set up 132-133 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyyy_0 = cbuffer.data(hs_geom_11_off + 132 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyy_0, g_z_x_xxyy_y, g_z_x_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyyy_0[k] = -g_z_x_xxyy_0[k] * ab_y + g_z_x_xxyy_y[k];
            }

            /// Set up 133-134 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyyz_0 = cbuffer.data(hs_geom_11_off + 133 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyyz_0, g_z_x_xxyz_0, g_z_x_xxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyyz_0[k] = -g_z_x_xxyz_0[k] * ab_y + g_z_x_xxyz_y[k];
            }

            /// Set up 134-135 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyzz_0 = cbuffer.data(hs_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyzz_0, g_z_x_xxzz_0, g_z_x_xxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyzz_0[k] = -g_z_x_xxzz_0[k] * ab_y + g_z_x_xxzz_y[k];
            }

            /// Set up 135-136 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxzzz_0 = cbuffer.data(hs_geom_11_off + 135 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_0, g_z_x_xxzzz_0, g_z_x_xzzz_0, g_z_x_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxzzz_0[k] = g_z_0_xzzz_0[k] - g_z_x_xzzz_0[k] * ab_x + g_z_x_xzzz_x[k];
            }

            /// Set up 136-137 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyyy_0 = cbuffer.data(hs_geom_11_off + 136 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyy_0, g_z_x_xyyy_y, g_z_x_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyyy_0[k] = -g_z_x_xyyy_0[k] * ab_y + g_z_x_xyyy_y[k];
            }

            /// Set up 137-138 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyyz_0 = cbuffer.data(hs_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyyz_0, g_z_x_xyyz_0, g_z_x_xyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyyz_0[k] = -g_z_x_xyyz_0[k] * ab_y + g_z_x_xyyz_y[k];
            }

            /// Set up 138-139 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyzz_0 = cbuffer.data(hs_geom_11_off + 138 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyzz_0, g_z_x_xyzz_0, g_z_x_xyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyzz_0[k] = -g_z_x_xyzz_0[k] * ab_y + g_z_x_xyzz_y[k];
            }

            /// Set up 139-140 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyzzz_0 = cbuffer.data(hs_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyzzz_0, g_z_x_xzzz_0, g_z_x_xzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyzzz_0[k] = -g_z_x_xzzz_0[k] * ab_y + g_z_x_xzzz_y[k];
            }

            /// Set up 140-141 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzzzz_0 = cbuffer.data(hs_geom_11_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_0, g_z_x_xzzzz_0, g_z_x_zzzz_0, g_z_x_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzzzz_0[k] = g_z_0_zzzz_0[k] - g_z_x_zzzz_0[k] * ab_x + g_z_x_zzzz_x[k];
            }

            /// Set up 141-142 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyyy_0 = cbuffer.data(hs_geom_11_off + 141 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyy_0, g_z_x_yyyy_y, g_z_x_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyyy_0[k] = -g_z_x_yyyy_0[k] * ab_y + g_z_x_yyyy_y[k];
            }

            /// Set up 142-143 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyyz_0 = cbuffer.data(hs_geom_11_off + 142 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyyz_0, g_z_x_yyyz_0, g_z_x_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyyz_0[k] = -g_z_x_yyyz_0[k] * ab_y + g_z_x_yyyz_y[k];
            }

            /// Set up 143-144 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyzz_0 = cbuffer.data(hs_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyzz_0, g_z_x_yyzz_0, g_z_x_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyzz_0[k] = -g_z_x_yyzz_0[k] * ab_y + g_z_x_yyzz_y[k];
            }

            /// Set up 144-145 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyzzz_0 = cbuffer.data(hs_geom_11_off + 144 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyzzz_0, g_z_x_yzzz_0, g_z_x_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyzzz_0[k] = -g_z_x_yzzz_0[k] * ab_y + g_z_x_yzzz_y[k];
            }

            /// Set up 145-146 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzzzz_0 = cbuffer.data(hs_geom_11_off + 145 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzzzz_0, g_z_x_zzzz_0, g_z_x_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzzzz_0[k] = -g_z_x_zzzz_0[k] * ab_y + g_z_x_zzzz_y[k];
            }

            /// Set up 146-147 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzzzz_0 = cbuffer.data(hs_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_0, g_z_x_zzzz_0, g_z_x_zzzz_z, g_z_x_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzzzz_0[k] = -g_0_x_zzzz_0[k] - g_z_x_zzzz_0[k] * ab_z + g_z_x_zzzz_z[k];
            }

            /// Set up 147-148 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxx_0 = cbuffer.data(hs_geom_11_off + 147 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxx_0, g_z_y_xxxx_x, g_z_y_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxx_0[k] = -g_z_y_xxxx_0[k] * ab_x + g_z_y_xxxx_x[k];
            }

            /// Set up 148-149 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxy_0 = cbuffer.data(hs_geom_11_off + 148 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxxy_0, g_z_y_xxxy_0, g_z_y_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxy_0[k] = -g_z_y_xxxy_0[k] * ab_x + g_z_y_xxxy_x[k];
            }

            /// Set up 149-150 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxz_0 = cbuffer.data(hs_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxxz_0, g_z_y_xxxz_0, g_z_y_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxz_0[k] = -g_z_y_xxxz_0[k] * ab_x + g_z_y_xxxz_x[k];
            }

            /// Set up 150-151 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxyy_0 = cbuffer.data(hs_geom_11_off + 150 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxyy_0, g_z_y_xxyy_0, g_z_y_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxyy_0[k] = -g_z_y_xxyy_0[k] * ab_x + g_z_y_xxyy_x[k];
            }

            /// Set up 151-152 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxyz_0 = cbuffer.data(hs_geom_11_off + 151 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxyz_0, g_z_y_xxyz_0, g_z_y_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxyz_0[k] = -g_z_y_xxyz_0[k] * ab_x + g_z_y_xxyz_x[k];
            }

            /// Set up 152-153 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxzz_0 = cbuffer.data(hs_geom_11_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxzz_0, g_z_y_xxzz_0, g_z_y_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxzz_0[k] = -g_z_y_xxzz_0[k] * ab_x + g_z_y_xxzz_x[k];
            }

            /// Set up 153-154 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyyy_0 = cbuffer.data(hs_geom_11_off + 153 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyyy_0, g_z_y_xyyy_0, g_z_y_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyyy_0[k] = -g_z_y_xyyy_0[k] * ab_x + g_z_y_xyyy_x[k];
            }

            /// Set up 154-155 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyyz_0 = cbuffer.data(hs_geom_11_off + 154 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyyz_0, g_z_y_xyyz_0, g_z_y_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyyz_0[k] = -g_z_y_xyyz_0[k] * ab_x + g_z_y_xyyz_x[k];
            }

            /// Set up 155-156 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyzz_0 = cbuffer.data(hs_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyzz_0, g_z_y_xyzz_0, g_z_y_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyzz_0[k] = -g_z_y_xyzz_0[k] * ab_x + g_z_y_xyzz_x[k];
            }

            /// Set up 156-157 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxzzz_0 = cbuffer.data(hs_geom_11_off + 156 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxzzz_0, g_z_y_xzzz_0, g_z_y_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxzzz_0[k] = -g_z_y_xzzz_0[k] * ab_x + g_z_y_xzzz_x[k];
            }

            /// Set up 157-158 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyyy_0 = cbuffer.data(hs_geom_11_off + 157 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyyy_0, g_z_y_yyyy_0, g_z_y_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyyy_0[k] = -g_z_y_yyyy_0[k] * ab_x + g_z_y_yyyy_x[k];
            }

            /// Set up 158-159 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyyz_0 = cbuffer.data(hs_geom_11_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyyz_0, g_z_y_yyyz_0, g_z_y_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyyz_0[k] = -g_z_y_yyyz_0[k] * ab_x + g_z_y_yyyz_x[k];
            }

            /// Set up 159-160 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyzz_0 = cbuffer.data(hs_geom_11_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyzz_0, g_z_y_yyzz_0, g_z_y_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyzz_0[k] = -g_z_y_yyzz_0[k] * ab_x + g_z_y_yyzz_x[k];
            }

            /// Set up 160-161 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyzzz_0 = cbuffer.data(hs_geom_11_off + 160 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyzzz_0, g_z_y_yzzz_0, g_z_y_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyzzz_0[k] = -g_z_y_yzzz_0[k] * ab_x + g_z_y_yzzz_x[k];
            }

            /// Set up 161-162 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzzzz_0 = cbuffer.data(hs_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzzzz_0, g_z_y_zzzz_0, g_z_y_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzzzz_0[k] = -g_z_y_zzzz_0[k] * ab_x + g_z_y_zzzz_x[k];
            }

            /// Set up 162-163 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyyy_0 = cbuffer.data(hs_geom_11_off + 162 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyy_0, g_z_y_yyyy_0, g_z_y_yyyy_y, g_z_y_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyyy_0[k] = g_z_0_yyyy_0[k] - g_z_y_yyyy_0[k] * ab_y + g_z_y_yyyy_y[k];
            }

            /// Set up 163-164 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyyz_0 = cbuffer.data(hs_geom_11_off + 163 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_0, g_z_y_yyyyz_0, g_z_y_yyyz_0, g_z_y_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyyz_0[k] = g_z_0_yyyz_0[k] - g_z_y_yyyz_0[k] * ab_y + g_z_y_yyyz_y[k];
            }

            /// Set up 164-165 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyzz_0 = cbuffer.data(hs_geom_11_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_0, g_z_y_yyyzz_0, g_z_y_yyzz_0, g_z_y_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyzz_0[k] = g_z_0_yyzz_0[k] - g_z_y_yyzz_0[k] * ab_y + g_z_y_yyzz_y[k];
            }

            /// Set up 165-166 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyzzz_0 = cbuffer.data(hs_geom_11_off + 165 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_0, g_z_y_yyzzz_0, g_z_y_yzzz_0, g_z_y_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyzzz_0[k] = g_z_0_yzzz_0[k] - g_z_y_yzzz_0[k] * ab_y + g_z_y_yzzz_y[k];
            }

            /// Set up 166-167 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzzzz_0 = cbuffer.data(hs_geom_11_off + 166 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_0, g_z_y_yzzzz_0, g_z_y_zzzz_0, g_z_y_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzzzz_0[k] = g_z_0_zzzz_0[k] - g_z_y_zzzz_0[k] * ab_y + g_z_y_zzzz_y[k];
            }

            /// Set up 167-168 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzzzz_0 = cbuffer.data(hs_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_0, g_z_y_zzzz_0, g_z_y_zzzz_z, g_z_y_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzzzz_0[k] = -g_0_y_zzzz_0[k] - g_z_y_zzzz_0[k] * ab_z + g_z_y_zzzz_z[k];
            }

            /// Set up 168-169 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxx_0 = cbuffer.data(hs_geom_11_off + 168 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxx_0, g_z_z_xxxx_x, g_z_z_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxx_0[k] = -g_z_z_xxxx_0[k] * ab_x + g_z_z_xxxx_x[k];
            }

            /// Set up 169-170 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxy_0 = cbuffer.data(hs_geom_11_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxxy_0, g_z_z_xxxy_0, g_z_z_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxy_0[k] = -g_z_z_xxxy_0[k] * ab_x + g_z_z_xxxy_x[k];
            }

            /// Set up 170-171 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxz_0 = cbuffer.data(hs_geom_11_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxxz_0, g_z_z_xxxz_0, g_z_z_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxz_0[k] = -g_z_z_xxxz_0[k] * ab_x + g_z_z_xxxz_x[k];
            }

            /// Set up 171-172 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxyy_0 = cbuffer.data(hs_geom_11_off + 171 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxyy_0, g_z_z_xxyy_0, g_z_z_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxyy_0[k] = -g_z_z_xxyy_0[k] * ab_x + g_z_z_xxyy_x[k];
            }

            /// Set up 172-173 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxyz_0 = cbuffer.data(hs_geom_11_off + 172 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxyz_0, g_z_z_xxyz_0, g_z_z_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxyz_0[k] = -g_z_z_xxyz_0[k] * ab_x + g_z_z_xxyz_x[k];
            }

            /// Set up 173-174 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxzz_0 = cbuffer.data(hs_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxzz_0, g_z_z_xxzz_0, g_z_z_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxzz_0[k] = -g_z_z_xxzz_0[k] * ab_x + g_z_z_xxzz_x[k];
            }

            /// Set up 174-175 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyyy_0 = cbuffer.data(hs_geom_11_off + 174 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyyy_0, g_z_z_xyyy_0, g_z_z_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyyy_0[k] = -g_z_z_xyyy_0[k] * ab_x + g_z_z_xyyy_x[k];
            }

            /// Set up 175-176 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyyz_0 = cbuffer.data(hs_geom_11_off + 175 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyyz_0, g_z_z_xyyz_0, g_z_z_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyyz_0[k] = -g_z_z_xyyz_0[k] * ab_x + g_z_z_xyyz_x[k];
            }

            /// Set up 176-177 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyzz_0 = cbuffer.data(hs_geom_11_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyzz_0, g_z_z_xyzz_0, g_z_z_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyzz_0[k] = -g_z_z_xyzz_0[k] * ab_x + g_z_z_xyzz_x[k];
            }

            /// Set up 177-178 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxzzz_0 = cbuffer.data(hs_geom_11_off + 177 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxzzz_0, g_z_z_xzzz_0, g_z_z_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxzzz_0[k] = -g_z_z_xzzz_0[k] * ab_x + g_z_z_xzzz_x[k];
            }

            /// Set up 178-179 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyyy_0 = cbuffer.data(hs_geom_11_off + 178 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyyy_0, g_z_z_yyyy_0, g_z_z_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyyy_0[k] = -g_z_z_yyyy_0[k] * ab_x + g_z_z_yyyy_x[k];
            }

            /// Set up 179-180 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyyz_0 = cbuffer.data(hs_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyyz_0, g_z_z_yyyz_0, g_z_z_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyyz_0[k] = -g_z_z_yyyz_0[k] * ab_x + g_z_z_yyyz_x[k];
            }

            /// Set up 180-181 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyzz_0 = cbuffer.data(hs_geom_11_off + 180 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyzz_0, g_z_z_yyzz_0, g_z_z_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyzz_0[k] = -g_z_z_yyzz_0[k] * ab_x + g_z_z_yyzz_x[k];
            }

            /// Set up 181-182 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyzzz_0 = cbuffer.data(hs_geom_11_off + 181 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyzzz_0, g_z_z_yzzz_0, g_z_z_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyzzz_0[k] = -g_z_z_yzzz_0[k] * ab_x + g_z_z_yzzz_x[k];
            }

            /// Set up 182-183 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzzzz_0 = cbuffer.data(hs_geom_11_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzzzz_0, g_z_z_zzzz_0, g_z_z_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzzzz_0[k] = -g_z_z_zzzz_0[k] * ab_x + g_z_z_zzzz_x[k];
            }

            /// Set up 183-184 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyyy_0 = cbuffer.data(hs_geom_11_off + 183 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyy_0, g_z_z_yyyy_y, g_z_z_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyyy_0[k] = -g_z_z_yyyy_0[k] * ab_y + g_z_z_yyyy_y[k];
            }

            /// Set up 184-185 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyyz_0 = cbuffer.data(hs_geom_11_off + 184 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyyz_0, g_z_z_yyyz_0, g_z_z_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyyz_0[k] = -g_z_z_yyyz_0[k] * ab_y + g_z_z_yyyz_y[k];
            }

            /// Set up 185-186 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyzz_0 = cbuffer.data(hs_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyzz_0, g_z_z_yyzz_0, g_z_z_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyzz_0[k] = -g_z_z_yyzz_0[k] * ab_y + g_z_z_yyzz_y[k];
            }

            /// Set up 186-187 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyzzz_0 = cbuffer.data(hs_geom_11_off + 186 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyzzz_0, g_z_z_yzzz_0, g_z_z_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyzzz_0[k] = -g_z_z_yzzz_0[k] * ab_y + g_z_z_yzzz_y[k];
            }

            /// Set up 187-188 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzzzz_0 = cbuffer.data(hs_geom_11_off + 187 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzzzz_0, g_z_z_zzzz_0, g_z_z_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzzzz_0[k] = -g_z_z_zzzz_0[k] * ab_y + g_z_z_zzzz_y[k];
            }

            /// Set up 188-189 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzzzz_0 = cbuffer.data(hs_geom_11_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_0, g_z_0_zzzz_0, g_z_z_zzzz_0, g_z_z_zzzz_z, g_z_z_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzzzz_0[k] = -g_0_z_zzzz_0[k] + g_z_0_zzzz_0[k] - g_z_z_zzzz_0[k] * ab_z + g_z_z_zzzz_z[k];
            }
        }
    }
}

} // erirec namespace

