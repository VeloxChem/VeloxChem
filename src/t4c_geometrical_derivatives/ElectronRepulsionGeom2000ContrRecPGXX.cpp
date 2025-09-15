#include "ElectronRepulsionGeom2000ContrRecPGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_pgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_pgxx,
                                            const size_t idx_geom_10_sgxx,
                                            const size_t idx_geom_20_sgxx,
                                            const size_t idx_geom_20_shxx,
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

            const auto sg_geom_10_off = idx_geom_10_sgxx + i * dcomps + j;

            auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 29 * ccomps * dcomps);

            auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * ccomps * dcomps);

            auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 31 * ccomps * dcomps);

            auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 32 * ccomps * dcomps);

            auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 33 * ccomps * dcomps);

            auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 34 * ccomps * dcomps);

            auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 44 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_20_off = idx_geom_20_sgxx + i * dcomps + j;

            auto g_xx_0_0_xxxx = cbuffer.data(sg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxy = cbuffer.data(sg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxz = cbuffer.data(sg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxyy = cbuffer.data(sg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxyz = cbuffer.data(sg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxzz = cbuffer.data(sg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xyyy = cbuffer.data(sg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xyyz = cbuffer.data(sg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xyzz = cbuffer.data(sg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xzzz = cbuffer.data(sg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_yyyy = cbuffer.data(sg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_yyyz = cbuffer.data(sg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_yyzz = cbuffer.data(sg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_yzzz = cbuffer.data(sg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_zzzz = cbuffer.data(sg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xy_0_0_xxxx = cbuffer.data(sg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xy_0_0_xxxy = cbuffer.data(sg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xy_0_0_xxxz = cbuffer.data(sg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_0_xxyy = cbuffer.data(sg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_0_xxyz = cbuffer.data(sg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xy_0_0_xxzz = cbuffer.data(sg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_0_xyyy = cbuffer.data(sg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_0_xyyz = cbuffer.data(sg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_0_xyzz = cbuffer.data(sg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_0_xzzz = cbuffer.data(sg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_0_yyyy = cbuffer.data(sg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_0_yyyz = cbuffer.data(sg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_0_yyzz = cbuffer.data(sg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_0_yzzz = cbuffer.data(sg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_0_zzzz = cbuffer.data(sg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 44 * ccomps * dcomps);

            auto g_yy_0_0_xxxx = cbuffer.data(sg_geom_20_off + 45 * ccomps * dcomps);

            auto g_yy_0_0_xxxy = cbuffer.data(sg_geom_20_off + 46 * ccomps * dcomps);

            auto g_yy_0_0_xxxz = cbuffer.data(sg_geom_20_off + 47 * ccomps * dcomps);

            auto g_yy_0_0_xxyy = cbuffer.data(sg_geom_20_off + 48 * ccomps * dcomps);

            auto g_yy_0_0_xxyz = cbuffer.data(sg_geom_20_off + 49 * ccomps * dcomps);

            auto g_yy_0_0_xxzz = cbuffer.data(sg_geom_20_off + 50 * ccomps * dcomps);

            auto g_yy_0_0_xyyy = cbuffer.data(sg_geom_20_off + 51 * ccomps * dcomps);

            auto g_yy_0_0_xyyz = cbuffer.data(sg_geom_20_off + 52 * ccomps * dcomps);

            auto g_yy_0_0_xyzz = cbuffer.data(sg_geom_20_off + 53 * ccomps * dcomps);

            auto g_yy_0_0_xzzz = cbuffer.data(sg_geom_20_off + 54 * ccomps * dcomps);

            auto g_yy_0_0_yyyy = cbuffer.data(sg_geom_20_off + 55 * ccomps * dcomps);

            auto g_yy_0_0_yyyz = cbuffer.data(sg_geom_20_off + 56 * ccomps * dcomps);

            auto g_yy_0_0_yyzz = cbuffer.data(sg_geom_20_off + 57 * ccomps * dcomps);

            auto g_yy_0_0_yzzz = cbuffer.data(sg_geom_20_off + 58 * ccomps * dcomps);

            auto g_yy_0_0_zzzz = cbuffer.data(sg_geom_20_off + 59 * ccomps * dcomps);

            auto g_yz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 60 * ccomps * dcomps);

            auto g_yz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 61 * ccomps * dcomps);

            auto g_yz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 62 * ccomps * dcomps);

            auto g_yz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 63 * ccomps * dcomps);

            auto g_yz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 64 * ccomps * dcomps);

            auto g_yz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 65 * ccomps * dcomps);

            auto g_yz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 66 * ccomps * dcomps);

            auto g_yz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 67 * ccomps * dcomps);

            auto g_yz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 68 * ccomps * dcomps);

            auto g_yz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 69 * ccomps * dcomps);

            auto g_yz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 70 * ccomps * dcomps);

            auto g_yz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 71 * ccomps * dcomps);

            auto g_yz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 72 * ccomps * dcomps);

            auto g_yz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 73 * ccomps * dcomps);

            auto g_yz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 74 * ccomps * dcomps);

            auto g_zz_0_0_xxxx = cbuffer.data(sg_geom_20_off + 75 * ccomps * dcomps);

            auto g_zz_0_0_xxxy = cbuffer.data(sg_geom_20_off + 76 * ccomps * dcomps);

            auto g_zz_0_0_xxxz = cbuffer.data(sg_geom_20_off + 77 * ccomps * dcomps);

            auto g_zz_0_0_xxyy = cbuffer.data(sg_geom_20_off + 78 * ccomps * dcomps);

            auto g_zz_0_0_xxyz = cbuffer.data(sg_geom_20_off + 79 * ccomps * dcomps);

            auto g_zz_0_0_xxzz = cbuffer.data(sg_geom_20_off + 80 * ccomps * dcomps);

            auto g_zz_0_0_xyyy = cbuffer.data(sg_geom_20_off + 81 * ccomps * dcomps);

            auto g_zz_0_0_xyyz = cbuffer.data(sg_geom_20_off + 82 * ccomps * dcomps);

            auto g_zz_0_0_xyzz = cbuffer.data(sg_geom_20_off + 83 * ccomps * dcomps);

            auto g_zz_0_0_xzzz = cbuffer.data(sg_geom_20_off + 84 * ccomps * dcomps);

            auto g_zz_0_0_yyyy = cbuffer.data(sg_geom_20_off + 85 * ccomps * dcomps);

            auto g_zz_0_0_yyyz = cbuffer.data(sg_geom_20_off + 86 * ccomps * dcomps);

            auto g_zz_0_0_yyzz = cbuffer.data(sg_geom_20_off + 87 * ccomps * dcomps);

            auto g_zz_0_0_yzzz = cbuffer.data(sg_geom_20_off + 88 * ccomps * dcomps);

            auto g_zz_0_0_zzzz = cbuffer.data(sg_geom_20_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_20_off = idx_geom_20_shxx + i * dcomps + j;

            auto g_xx_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 44 * ccomps * dcomps);

            auto g_xz_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 45 * ccomps * dcomps);

            auto g_xz_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 46 * ccomps * dcomps);

            auto g_xz_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 47 * ccomps * dcomps);

            auto g_xz_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 48 * ccomps * dcomps);

            auto g_xz_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 49 * ccomps * dcomps);

            auto g_xz_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 50 * ccomps * dcomps);

            auto g_xz_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 51 * ccomps * dcomps);

            auto g_xz_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 52 * ccomps * dcomps);

            auto g_xz_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 53 * ccomps * dcomps);

            auto g_xz_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 54 * ccomps * dcomps);

            auto g_xz_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 55 * ccomps * dcomps);

            auto g_xz_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 56 * ccomps * dcomps);

            auto g_xz_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 57 * ccomps * dcomps);

            auto g_xz_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 58 * ccomps * dcomps);

            auto g_xz_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 62 * ccomps * dcomps);

            auto g_yy_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 63 * ccomps * dcomps);

            auto g_yy_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 64 * ccomps * dcomps);

            auto g_yy_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 65 * ccomps * dcomps);

            auto g_yy_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 66 * ccomps * dcomps);

            auto g_yy_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 67 * ccomps * dcomps);

            auto g_yy_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 68 * ccomps * dcomps);

            auto g_yy_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 69 * ccomps * dcomps);

            auto g_yy_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 70 * ccomps * dcomps);

            auto g_yy_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 71 * ccomps * dcomps);

            auto g_yy_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 72 * ccomps * dcomps);

            auto g_yy_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 73 * ccomps * dcomps);

            auto g_yy_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 74 * ccomps * dcomps);

            auto g_yy_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 75 * ccomps * dcomps);

            auto g_yy_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 76 * ccomps * dcomps);

            auto g_yy_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 77 * ccomps * dcomps);

            auto g_yy_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 78 * ccomps * dcomps);

            auto g_yy_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 79 * ccomps * dcomps);

            auto g_yy_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 80 * ccomps * dcomps);

            auto g_yy_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 81 * ccomps * dcomps);

            auto g_yy_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 82 * ccomps * dcomps);

            auto g_yy_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 83 * ccomps * dcomps);

            auto g_yz_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 84 * ccomps * dcomps);

            auto g_yz_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 85 * ccomps * dcomps);

            auto g_yz_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 86 * ccomps * dcomps);

            auto g_yz_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 87 * ccomps * dcomps);

            auto g_yz_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 88 * ccomps * dcomps);

            auto g_yz_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 89 * ccomps * dcomps);

            auto g_yz_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 90 * ccomps * dcomps);

            auto g_yz_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 91 * ccomps * dcomps);

            auto g_yz_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 92 * ccomps * dcomps);

            auto g_yz_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 93 * ccomps * dcomps);

            auto g_yz_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 94 * ccomps * dcomps);

            auto g_yz_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 95 * ccomps * dcomps);

            auto g_yz_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 96 * ccomps * dcomps);

            auto g_yz_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 97 * ccomps * dcomps);

            auto g_yz_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 98 * ccomps * dcomps);

            auto g_yz_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 99 * ccomps * dcomps);

            auto g_yz_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 100 * ccomps * dcomps);

            auto g_yz_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 101 * ccomps * dcomps);

            auto g_yz_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 102 * ccomps * dcomps);

            auto g_yz_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 103 * ccomps * dcomps);

            auto g_yz_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 104 * ccomps * dcomps);

            auto g_zz_0_0_xxxxx = cbuffer.data(sh_geom_20_off + 105 * ccomps * dcomps);

            auto g_zz_0_0_xxxxy = cbuffer.data(sh_geom_20_off + 106 * ccomps * dcomps);

            auto g_zz_0_0_xxxxz = cbuffer.data(sh_geom_20_off + 107 * ccomps * dcomps);

            auto g_zz_0_0_xxxyy = cbuffer.data(sh_geom_20_off + 108 * ccomps * dcomps);

            auto g_zz_0_0_xxxyz = cbuffer.data(sh_geom_20_off + 109 * ccomps * dcomps);

            auto g_zz_0_0_xxxzz = cbuffer.data(sh_geom_20_off + 110 * ccomps * dcomps);

            auto g_zz_0_0_xxyyy = cbuffer.data(sh_geom_20_off + 111 * ccomps * dcomps);

            auto g_zz_0_0_xxyyz = cbuffer.data(sh_geom_20_off + 112 * ccomps * dcomps);

            auto g_zz_0_0_xxyzz = cbuffer.data(sh_geom_20_off + 113 * ccomps * dcomps);

            auto g_zz_0_0_xxzzz = cbuffer.data(sh_geom_20_off + 114 * ccomps * dcomps);

            auto g_zz_0_0_xyyyy = cbuffer.data(sh_geom_20_off + 115 * ccomps * dcomps);

            auto g_zz_0_0_xyyyz = cbuffer.data(sh_geom_20_off + 116 * ccomps * dcomps);

            auto g_zz_0_0_xyyzz = cbuffer.data(sh_geom_20_off + 117 * ccomps * dcomps);

            auto g_zz_0_0_xyzzz = cbuffer.data(sh_geom_20_off + 118 * ccomps * dcomps);

            auto g_zz_0_0_xzzzz = cbuffer.data(sh_geom_20_off + 119 * ccomps * dcomps);

            auto g_zz_0_0_yyyyy = cbuffer.data(sh_geom_20_off + 120 * ccomps * dcomps);

            auto g_zz_0_0_yyyyz = cbuffer.data(sh_geom_20_off + 121 * ccomps * dcomps);

            auto g_zz_0_0_yyyzz = cbuffer.data(sh_geom_20_off + 122 * ccomps * dcomps);

            auto g_zz_0_0_yyzzz = cbuffer.data(sh_geom_20_off + 123 * ccomps * dcomps);

            auto g_zz_0_0_yzzzz = cbuffer.data(sh_geom_20_off + 124 * ccomps * dcomps);

            auto g_zz_0_0_zzzzz = cbuffer.data(sh_geom_20_off + 125 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pgxx

            const auto pg_geom_20_off = idx_geom_20_pgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz, g_xx_0_0_xxxx, g_xx_0_0_xxxxx, g_xx_0_0_xxxxy, g_xx_0_0_xxxxz, g_xx_0_0_xxxy, g_xx_0_0_xxxyy, g_xx_0_0_xxxyz, g_xx_0_0_xxxz, g_xx_0_0_xxxzz, g_xx_0_0_xxyy, g_xx_0_0_xxyyy, g_xx_0_0_xxyyz, g_xx_0_0_xxyz, g_xx_0_0_xxyzz, g_xx_0_0_xxzz, g_xx_0_0_xxzzz, g_xx_0_0_xyyy, g_xx_0_0_xyyyy, g_xx_0_0_xyyyz, g_xx_0_0_xyyz, g_xx_0_0_xyyzz, g_xx_0_0_xyzz, g_xx_0_0_xyzzz, g_xx_0_0_xzzz, g_xx_0_0_xzzzz, g_xx_0_0_yyyy, g_xx_0_0_yyyz, g_xx_0_0_yyzz, g_xx_0_0_yzzz, g_xx_0_0_zzzz, g_xx_0_x_xxxx, g_xx_0_x_xxxy, g_xx_0_x_xxxz, g_xx_0_x_xxyy, g_xx_0_x_xxyz, g_xx_0_x_xxzz, g_xx_0_x_xyyy, g_xx_0_x_xyyz, g_xx_0_x_xyzz, g_xx_0_x_xzzz, g_xx_0_x_yyyy, g_xx_0_x_yyyz, g_xx_0_x_yyzz, g_xx_0_x_yzzz, g_xx_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_x_xxxx[k] = -2.0 * g_x_0_0_xxxx[k] - g_xx_0_0_xxxx[k] * ab_x + g_xx_0_0_xxxxx[k];

                g_xx_0_x_xxxy[k] = -2.0 * g_x_0_0_xxxy[k] - g_xx_0_0_xxxy[k] * ab_x + g_xx_0_0_xxxxy[k];

                g_xx_0_x_xxxz[k] = -2.0 * g_x_0_0_xxxz[k] - g_xx_0_0_xxxz[k] * ab_x + g_xx_0_0_xxxxz[k];

                g_xx_0_x_xxyy[k] = -2.0 * g_x_0_0_xxyy[k] - g_xx_0_0_xxyy[k] * ab_x + g_xx_0_0_xxxyy[k];

                g_xx_0_x_xxyz[k] = -2.0 * g_x_0_0_xxyz[k] - g_xx_0_0_xxyz[k] * ab_x + g_xx_0_0_xxxyz[k];

                g_xx_0_x_xxzz[k] = -2.0 * g_x_0_0_xxzz[k] - g_xx_0_0_xxzz[k] * ab_x + g_xx_0_0_xxxzz[k];

                g_xx_0_x_xyyy[k] = -2.0 * g_x_0_0_xyyy[k] - g_xx_0_0_xyyy[k] * ab_x + g_xx_0_0_xxyyy[k];

                g_xx_0_x_xyyz[k] = -2.0 * g_x_0_0_xyyz[k] - g_xx_0_0_xyyz[k] * ab_x + g_xx_0_0_xxyyz[k];

                g_xx_0_x_xyzz[k] = -2.0 * g_x_0_0_xyzz[k] - g_xx_0_0_xyzz[k] * ab_x + g_xx_0_0_xxyzz[k];

                g_xx_0_x_xzzz[k] = -2.0 * g_x_0_0_xzzz[k] - g_xx_0_0_xzzz[k] * ab_x + g_xx_0_0_xxzzz[k];

                g_xx_0_x_yyyy[k] = -2.0 * g_x_0_0_yyyy[k] - g_xx_0_0_yyyy[k] * ab_x + g_xx_0_0_xyyyy[k];

                g_xx_0_x_yyyz[k] = -2.0 * g_x_0_0_yyyz[k] - g_xx_0_0_yyyz[k] * ab_x + g_xx_0_0_xyyyz[k];

                g_xx_0_x_yyzz[k] = -2.0 * g_x_0_0_yyzz[k] - g_xx_0_0_yyzz[k] * ab_x + g_xx_0_0_xyyzz[k];

                g_xx_0_x_yzzz[k] = -2.0 * g_x_0_0_yzzz[k] - g_xx_0_0_yzzz[k] * ab_x + g_xx_0_0_xyzzz[k];

                g_xx_0_x_zzzz[k] = -2.0 * g_x_0_0_zzzz[k] - g_xx_0_0_zzzz[k] * ab_x + g_xx_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xx_0_0_xxxx, g_xx_0_0_xxxxy, g_xx_0_0_xxxy, g_xx_0_0_xxxyy, g_xx_0_0_xxxyz, g_xx_0_0_xxxz, g_xx_0_0_xxyy, g_xx_0_0_xxyyy, g_xx_0_0_xxyyz, g_xx_0_0_xxyz, g_xx_0_0_xxyzz, g_xx_0_0_xxzz, g_xx_0_0_xyyy, g_xx_0_0_xyyyy, g_xx_0_0_xyyyz, g_xx_0_0_xyyz, g_xx_0_0_xyyzz, g_xx_0_0_xyzz, g_xx_0_0_xyzzz, g_xx_0_0_xzzz, g_xx_0_0_yyyy, g_xx_0_0_yyyyy, g_xx_0_0_yyyyz, g_xx_0_0_yyyz, g_xx_0_0_yyyzz, g_xx_0_0_yyzz, g_xx_0_0_yyzzz, g_xx_0_0_yzzz, g_xx_0_0_yzzzz, g_xx_0_0_zzzz, g_xx_0_y_xxxx, g_xx_0_y_xxxy, g_xx_0_y_xxxz, g_xx_0_y_xxyy, g_xx_0_y_xxyz, g_xx_0_y_xxzz, g_xx_0_y_xyyy, g_xx_0_y_xyyz, g_xx_0_y_xyzz, g_xx_0_y_xzzz, g_xx_0_y_yyyy, g_xx_0_y_yyyz, g_xx_0_y_yyzz, g_xx_0_y_yzzz, g_xx_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_y_xxxx[k] = -g_xx_0_0_xxxx[k] * ab_y + g_xx_0_0_xxxxy[k];

                g_xx_0_y_xxxy[k] = -g_xx_0_0_xxxy[k] * ab_y + g_xx_0_0_xxxyy[k];

                g_xx_0_y_xxxz[k] = -g_xx_0_0_xxxz[k] * ab_y + g_xx_0_0_xxxyz[k];

                g_xx_0_y_xxyy[k] = -g_xx_0_0_xxyy[k] * ab_y + g_xx_0_0_xxyyy[k];

                g_xx_0_y_xxyz[k] = -g_xx_0_0_xxyz[k] * ab_y + g_xx_0_0_xxyyz[k];

                g_xx_0_y_xxzz[k] = -g_xx_0_0_xxzz[k] * ab_y + g_xx_0_0_xxyzz[k];

                g_xx_0_y_xyyy[k] = -g_xx_0_0_xyyy[k] * ab_y + g_xx_0_0_xyyyy[k];

                g_xx_0_y_xyyz[k] = -g_xx_0_0_xyyz[k] * ab_y + g_xx_0_0_xyyyz[k];

                g_xx_0_y_xyzz[k] = -g_xx_0_0_xyzz[k] * ab_y + g_xx_0_0_xyyzz[k];

                g_xx_0_y_xzzz[k] = -g_xx_0_0_xzzz[k] * ab_y + g_xx_0_0_xyzzz[k];

                g_xx_0_y_yyyy[k] = -g_xx_0_0_yyyy[k] * ab_y + g_xx_0_0_yyyyy[k];

                g_xx_0_y_yyyz[k] = -g_xx_0_0_yyyz[k] * ab_y + g_xx_0_0_yyyyz[k];

                g_xx_0_y_yyzz[k] = -g_xx_0_0_yyzz[k] * ab_y + g_xx_0_0_yyyzz[k];

                g_xx_0_y_yzzz[k] = -g_xx_0_0_yzzz[k] * ab_y + g_xx_0_0_yyzzz[k];

                g_xx_0_y_zzzz[k] = -g_xx_0_0_zzzz[k] * ab_y + g_xx_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xx_0_0_xxxx, g_xx_0_0_xxxxz, g_xx_0_0_xxxy, g_xx_0_0_xxxyz, g_xx_0_0_xxxz, g_xx_0_0_xxxzz, g_xx_0_0_xxyy, g_xx_0_0_xxyyz, g_xx_0_0_xxyz, g_xx_0_0_xxyzz, g_xx_0_0_xxzz, g_xx_0_0_xxzzz, g_xx_0_0_xyyy, g_xx_0_0_xyyyz, g_xx_0_0_xyyz, g_xx_0_0_xyyzz, g_xx_0_0_xyzz, g_xx_0_0_xyzzz, g_xx_0_0_xzzz, g_xx_0_0_xzzzz, g_xx_0_0_yyyy, g_xx_0_0_yyyyz, g_xx_0_0_yyyz, g_xx_0_0_yyyzz, g_xx_0_0_yyzz, g_xx_0_0_yyzzz, g_xx_0_0_yzzz, g_xx_0_0_yzzzz, g_xx_0_0_zzzz, g_xx_0_0_zzzzz, g_xx_0_z_xxxx, g_xx_0_z_xxxy, g_xx_0_z_xxxz, g_xx_0_z_xxyy, g_xx_0_z_xxyz, g_xx_0_z_xxzz, g_xx_0_z_xyyy, g_xx_0_z_xyyz, g_xx_0_z_xyzz, g_xx_0_z_xzzz, g_xx_0_z_yyyy, g_xx_0_z_yyyz, g_xx_0_z_yyzz, g_xx_0_z_yzzz, g_xx_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_z_xxxx[k] = -g_xx_0_0_xxxx[k] * ab_z + g_xx_0_0_xxxxz[k];

                g_xx_0_z_xxxy[k] = -g_xx_0_0_xxxy[k] * ab_z + g_xx_0_0_xxxyz[k];

                g_xx_0_z_xxxz[k] = -g_xx_0_0_xxxz[k] * ab_z + g_xx_0_0_xxxzz[k];

                g_xx_0_z_xxyy[k] = -g_xx_0_0_xxyy[k] * ab_z + g_xx_0_0_xxyyz[k];

                g_xx_0_z_xxyz[k] = -g_xx_0_0_xxyz[k] * ab_z + g_xx_0_0_xxyzz[k];

                g_xx_0_z_xxzz[k] = -g_xx_0_0_xxzz[k] * ab_z + g_xx_0_0_xxzzz[k];

                g_xx_0_z_xyyy[k] = -g_xx_0_0_xyyy[k] * ab_z + g_xx_0_0_xyyyz[k];

                g_xx_0_z_xyyz[k] = -g_xx_0_0_xyyz[k] * ab_z + g_xx_0_0_xyyzz[k];

                g_xx_0_z_xyzz[k] = -g_xx_0_0_xyzz[k] * ab_z + g_xx_0_0_xyzzz[k];

                g_xx_0_z_xzzz[k] = -g_xx_0_0_xzzz[k] * ab_z + g_xx_0_0_xzzzz[k];

                g_xx_0_z_yyyy[k] = -g_xx_0_0_yyyy[k] * ab_z + g_xx_0_0_yyyyz[k];

                g_xx_0_z_yyyz[k] = -g_xx_0_0_yyyz[k] * ab_z + g_xx_0_0_yyyzz[k];

                g_xx_0_z_yyzz[k] = -g_xx_0_0_yyzz[k] * ab_z + g_xx_0_0_yyzzz[k];

                g_xx_0_z_yzzz[k] = -g_xx_0_0_yzzz[k] * ab_z + g_xx_0_0_yzzzz[k];

                g_xx_0_z_zzzz[k] = -g_xx_0_0_zzzz[k] * ab_z + g_xx_0_0_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xy_0_0_xxxx, g_xy_0_0_xxxxx, g_xy_0_0_xxxxy, g_xy_0_0_xxxxz, g_xy_0_0_xxxy, g_xy_0_0_xxxyy, g_xy_0_0_xxxyz, g_xy_0_0_xxxz, g_xy_0_0_xxxzz, g_xy_0_0_xxyy, g_xy_0_0_xxyyy, g_xy_0_0_xxyyz, g_xy_0_0_xxyz, g_xy_0_0_xxyzz, g_xy_0_0_xxzz, g_xy_0_0_xxzzz, g_xy_0_0_xyyy, g_xy_0_0_xyyyy, g_xy_0_0_xyyyz, g_xy_0_0_xyyz, g_xy_0_0_xyyzz, g_xy_0_0_xyzz, g_xy_0_0_xyzzz, g_xy_0_0_xzzz, g_xy_0_0_xzzzz, g_xy_0_0_yyyy, g_xy_0_0_yyyz, g_xy_0_0_yyzz, g_xy_0_0_yzzz, g_xy_0_0_zzzz, g_xy_0_x_xxxx, g_xy_0_x_xxxy, g_xy_0_x_xxxz, g_xy_0_x_xxyy, g_xy_0_x_xxyz, g_xy_0_x_xxzz, g_xy_0_x_xyyy, g_xy_0_x_xyyz, g_xy_0_x_xyzz, g_xy_0_x_xzzz, g_xy_0_x_yyyy, g_xy_0_x_yyyz, g_xy_0_x_yyzz, g_xy_0_x_yzzz, g_xy_0_x_zzzz, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_x_xxxx[k] = -g_y_0_0_xxxx[k] - g_xy_0_0_xxxx[k] * ab_x + g_xy_0_0_xxxxx[k];

                g_xy_0_x_xxxy[k] = -g_y_0_0_xxxy[k] - g_xy_0_0_xxxy[k] * ab_x + g_xy_0_0_xxxxy[k];

                g_xy_0_x_xxxz[k] = -g_y_0_0_xxxz[k] - g_xy_0_0_xxxz[k] * ab_x + g_xy_0_0_xxxxz[k];

                g_xy_0_x_xxyy[k] = -g_y_0_0_xxyy[k] - g_xy_0_0_xxyy[k] * ab_x + g_xy_0_0_xxxyy[k];

                g_xy_0_x_xxyz[k] = -g_y_0_0_xxyz[k] - g_xy_0_0_xxyz[k] * ab_x + g_xy_0_0_xxxyz[k];

                g_xy_0_x_xxzz[k] = -g_y_0_0_xxzz[k] - g_xy_0_0_xxzz[k] * ab_x + g_xy_0_0_xxxzz[k];

                g_xy_0_x_xyyy[k] = -g_y_0_0_xyyy[k] - g_xy_0_0_xyyy[k] * ab_x + g_xy_0_0_xxyyy[k];

                g_xy_0_x_xyyz[k] = -g_y_0_0_xyyz[k] - g_xy_0_0_xyyz[k] * ab_x + g_xy_0_0_xxyyz[k];

                g_xy_0_x_xyzz[k] = -g_y_0_0_xyzz[k] - g_xy_0_0_xyzz[k] * ab_x + g_xy_0_0_xxyzz[k];

                g_xy_0_x_xzzz[k] = -g_y_0_0_xzzz[k] - g_xy_0_0_xzzz[k] * ab_x + g_xy_0_0_xxzzz[k];

                g_xy_0_x_yyyy[k] = -g_y_0_0_yyyy[k] - g_xy_0_0_yyyy[k] * ab_x + g_xy_0_0_xyyyy[k];

                g_xy_0_x_yyyz[k] = -g_y_0_0_yyyz[k] - g_xy_0_0_yyyz[k] * ab_x + g_xy_0_0_xyyyz[k];

                g_xy_0_x_yyzz[k] = -g_y_0_0_yyzz[k] - g_xy_0_0_yyzz[k] * ab_x + g_xy_0_0_xyyzz[k];

                g_xy_0_x_yzzz[k] = -g_y_0_0_yzzz[k] - g_xy_0_0_yzzz[k] * ab_x + g_xy_0_0_xyzzz[k];

                g_xy_0_x_zzzz[k] = -g_y_0_0_zzzz[k] - g_xy_0_0_zzzz[k] * ab_x + g_xy_0_0_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz, g_xy_0_0_xxxx, g_xy_0_0_xxxxy, g_xy_0_0_xxxy, g_xy_0_0_xxxyy, g_xy_0_0_xxxyz, g_xy_0_0_xxxz, g_xy_0_0_xxyy, g_xy_0_0_xxyyy, g_xy_0_0_xxyyz, g_xy_0_0_xxyz, g_xy_0_0_xxyzz, g_xy_0_0_xxzz, g_xy_0_0_xyyy, g_xy_0_0_xyyyy, g_xy_0_0_xyyyz, g_xy_0_0_xyyz, g_xy_0_0_xyyzz, g_xy_0_0_xyzz, g_xy_0_0_xyzzz, g_xy_0_0_xzzz, g_xy_0_0_yyyy, g_xy_0_0_yyyyy, g_xy_0_0_yyyyz, g_xy_0_0_yyyz, g_xy_0_0_yyyzz, g_xy_0_0_yyzz, g_xy_0_0_yyzzz, g_xy_0_0_yzzz, g_xy_0_0_yzzzz, g_xy_0_0_zzzz, g_xy_0_y_xxxx, g_xy_0_y_xxxy, g_xy_0_y_xxxz, g_xy_0_y_xxyy, g_xy_0_y_xxyz, g_xy_0_y_xxzz, g_xy_0_y_xyyy, g_xy_0_y_xyyz, g_xy_0_y_xyzz, g_xy_0_y_xzzz, g_xy_0_y_yyyy, g_xy_0_y_yyyz, g_xy_0_y_yyzz, g_xy_0_y_yzzz, g_xy_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_y_xxxx[k] = -g_x_0_0_xxxx[k] - g_xy_0_0_xxxx[k] * ab_y + g_xy_0_0_xxxxy[k];

                g_xy_0_y_xxxy[k] = -g_x_0_0_xxxy[k] - g_xy_0_0_xxxy[k] * ab_y + g_xy_0_0_xxxyy[k];

                g_xy_0_y_xxxz[k] = -g_x_0_0_xxxz[k] - g_xy_0_0_xxxz[k] * ab_y + g_xy_0_0_xxxyz[k];

                g_xy_0_y_xxyy[k] = -g_x_0_0_xxyy[k] - g_xy_0_0_xxyy[k] * ab_y + g_xy_0_0_xxyyy[k];

                g_xy_0_y_xxyz[k] = -g_x_0_0_xxyz[k] - g_xy_0_0_xxyz[k] * ab_y + g_xy_0_0_xxyyz[k];

                g_xy_0_y_xxzz[k] = -g_x_0_0_xxzz[k] - g_xy_0_0_xxzz[k] * ab_y + g_xy_0_0_xxyzz[k];

                g_xy_0_y_xyyy[k] = -g_x_0_0_xyyy[k] - g_xy_0_0_xyyy[k] * ab_y + g_xy_0_0_xyyyy[k];

                g_xy_0_y_xyyz[k] = -g_x_0_0_xyyz[k] - g_xy_0_0_xyyz[k] * ab_y + g_xy_0_0_xyyyz[k];

                g_xy_0_y_xyzz[k] = -g_x_0_0_xyzz[k] - g_xy_0_0_xyzz[k] * ab_y + g_xy_0_0_xyyzz[k];

                g_xy_0_y_xzzz[k] = -g_x_0_0_xzzz[k] - g_xy_0_0_xzzz[k] * ab_y + g_xy_0_0_xyzzz[k];

                g_xy_0_y_yyyy[k] = -g_x_0_0_yyyy[k] - g_xy_0_0_yyyy[k] * ab_y + g_xy_0_0_yyyyy[k];

                g_xy_0_y_yyyz[k] = -g_x_0_0_yyyz[k] - g_xy_0_0_yyyz[k] * ab_y + g_xy_0_0_yyyyz[k];

                g_xy_0_y_yyzz[k] = -g_x_0_0_yyzz[k] - g_xy_0_0_yyzz[k] * ab_y + g_xy_0_0_yyyzz[k];

                g_xy_0_y_yzzz[k] = -g_x_0_0_yzzz[k] - g_xy_0_0_yzzz[k] * ab_y + g_xy_0_0_yyzzz[k];

                g_xy_0_y_zzzz[k] = -g_x_0_0_zzzz[k] - g_xy_0_0_zzzz[k] * ab_y + g_xy_0_0_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xy_0_0_xxxx, g_xy_0_0_xxxxz, g_xy_0_0_xxxy, g_xy_0_0_xxxyz, g_xy_0_0_xxxz, g_xy_0_0_xxxzz, g_xy_0_0_xxyy, g_xy_0_0_xxyyz, g_xy_0_0_xxyz, g_xy_0_0_xxyzz, g_xy_0_0_xxzz, g_xy_0_0_xxzzz, g_xy_0_0_xyyy, g_xy_0_0_xyyyz, g_xy_0_0_xyyz, g_xy_0_0_xyyzz, g_xy_0_0_xyzz, g_xy_0_0_xyzzz, g_xy_0_0_xzzz, g_xy_0_0_xzzzz, g_xy_0_0_yyyy, g_xy_0_0_yyyyz, g_xy_0_0_yyyz, g_xy_0_0_yyyzz, g_xy_0_0_yyzz, g_xy_0_0_yyzzz, g_xy_0_0_yzzz, g_xy_0_0_yzzzz, g_xy_0_0_zzzz, g_xy_0_0_zzzzz, g_xy_0_z_xxxx, g_xy_0_z_xxxy, g_xy_0_z_xxxz, g_xy_0_z_xxyy, g_xy_0_z_xxyz, g_xy_0_z_xxzz, g_xy_0_z_xyyy, g_xy_0_z_xyyz, g_xy_0_z_xyzz, g_xy_0_z_xzzz, g_xy_0_z_yyyy, g_xy_0_z_yyyz, g_xy_0_z_yyzz, g_xy_0_z_yzzz, g_xy_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_z_xxxx[k] = -g_xy_0_0_xxxx[k] * ab_z + g_xy_0_0_xxxxz[k];

                g_xy_0_z_xxxy[k] = -g_xy_0_0_xxxy[k] * ab_z + g_xy_0_0_xxxyz[k];

                g_xy_0_z_xxxz[k] = -g_xy_0_0_xxxz[k] * ab_z + g_xy_0_0_xxxzz[k];

                g_xy_0_z_xxyy[k] = -g_xy_0_0_xxyy[k] * ab_z + g_xy_0_0_xxyyz[k];

                g_xy_0_z_xxyz[k] = -g_xy_0_0_xxyz[k] * ab_z + g_xy_0_0_xxyzz[k];

                g_xy_0_z_xxzz[k] = -g_xy_0_0_xxzz[k] * ab_z + g_xy_0_0_xxzzz[k];

                g_xy_0_z_xyyy[k] = -g_xy_0_0_xyyy[k] * ab_z + g_xy_0_0_xyyyz[k];

                g_xy_0_z_xyyz[k] = -g_xy_0_0_xyyz[k] * ab_z + g_xy_0_0_xyyzz[k];

                g_xy_0_z_xyzz[k] = -g_xy_0_0_xyzz[k] * ab_z + g_xy_0_0_xyzzz[k];

                g_xy_0_z_xzzz[k] = -g_xy_0_0_xzzz[k] * ab_z + g_xy_0_0_xzzzz[k];

                g_xy_0_z_yyyy[k] = -g_xy_0_0_yyyy[k] * ab_z + g_xy_0_0_yyyyz[k];

                g_xy_0_z_yyyz[k] = -g_xy_0_0_yyyz[k] * ab_z + g_xy_0_0_yyyzz[k];

                g_xy_0_z_yyzz[k] = -g_xy_0_0_yyzz[k] * ab_z + g_xy_0_0_yyzzz[k];

                g_xy_0_z_yzzz[k] = -g_xy_0_0_yzzz[k] * ab_z + g_xy_0_0_yzzzz[k];

                g_xy_0_z_zzzz[k] = -g_xy_0_0_zzzz[k] * ab_z + g_xy_0_0_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xz_0_0_xxxx, g_xz_0_0_xxxxx, g_xz_0_0_xxxxy, g_xz_0_0_xxxxz, g_xz_0_0_xxxy, g_xz_0_0_xxxyy, g_xz_0_0_xxxyz, g_xz_0_0_xxxz, g_xz_0_0_xxxzz, g_xz_0_0_xxyy, g_xz_0_0_xxyyy, g_xz_0_0_xxyyz, g_xz_0_0_xxyz, g_xz_0_0_xxyzz, g_xz_0_0_xxzz, g_xz_0_0_xxzzz, g_xz_0_0_xyyy, g_xz_0_0_xyyyy, g_xz_0_0_xyyyz, g_xz_0_0_xyyz, g_xz_0_0_xyyzz, g_xz_0_0_xyzz, g_xz_0_0_xyzzz, g_xz_0_0_xzzz, g_xz_0_0_xzzzz, g_xz_0_0_yyyy, g_xz_0_0_yyyz, g_xz_0_0_yyzz, g_xz_0_0_yzzz, g_xz_0_0_zzzz, g_xz_0_x_xxxx, g_xz_0_x_xxxy, g_xz_0_x_xxxz, g_xz_0_x_xxyy, g_xz_0_x_xxyz, g_xz_0_x_xxzz, g_xz_0_x_xyyy, g_xz_0_x_xyyz, g_xz_0_x_xyzz, g_xz_0_x_xzzz, g_xz_0_x_yyyy, g_xz_0_x_yyyz, g_xz_0_x_yyzz, g_xz_0_x_yzzz, g_xz_0_x_zzzz, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_x_xxxx[k] = -g_z_0_0_xxxx[k] - g_xz_0_0_xxxx[k] * ab_x + g_xz_0_0_xxxxx[k];

                g_xz_0_x_xxxy[k] = -g_z_0_0_xxxy[k] - g_xz_0_0_xxxy[k] * ab_x + g_xz_0_0_xxxxy[k];

                g_xz_0_x_xxxz[k] = -g_z_0_0_xxxz[k] - g_xz_0_0_xxxz[k] * ab_x + g_xz_0_0_xxxxz[k];

                g_xz_0_x_xxyy[k] = -g_z_0_0_xxyy[k] - g_xz_0_0_xxyy[k] * ab_x + g_xz_0_0_xxxyy[k];

                g_xz_0_x_xxyz[k] = -g_z_0_0_xxyz[k] - g_xz_0_0_xxyz[k] * ab_x + g_xz_0_0_xxxyz[k];

                g_xz_0_x_xxzz[k] = -g_z_0_0_xxzz[k] - g_xz_0_0_xxzz[k] * ab_x + g_xz_0_0_xxxzz[k];

                g_xz_0_x_xyyy[k] = -g_z_0_0_xyyy[k] - g_xz_0_0_xyyy[k] * ab_x + g_xz_0_0_xxyyy[k];

                g_xz_0_x_xyyz[k] = -g_z_0_0_xyyz[k] - g_xz_0_0_xyyz[k] * ab_x + g_xz_0_0_xxyyz[k];

                g_xz_0_x_xyzz[k] = -g_z_0_0_xyzz[k] - g_xz_0_0_xyzz[k] * ab_x + g_xz_0_0_xxyzz[k];

                g_xz_0_x_xzzz[k] = -g_z_0_0_xzzz[k] - g_xz_0_0_xzzz[k] * ab_x + g_xz_0_0_xxzzz[k];

                g_xz_0_x_yyyy[k] = -g_z_0_0_yyyy[k] - g_xz_0_0_yyyy[k] * ab_x + g_xz_0_0_xyyyy[k];

                g_xz_0_x_yyyz[k] = -g_z_0_0_yyyz[k] - g_xz_0_0_yyyz[k] * ab_x + g_xz_0_0_xyyyz[k];

                g_xz_0_x_yyzz[k] = -g_z_0_0_yyzz[k] - g_xz_0_0_yyzz[k] * ab_x + g_xz_0_0_xyyzz[k];

                g_xz_0_x_yzzz[k] = -g_z_0_0_yzzz[k] - g_xz_0_0_yzzz[k] * ab_x + g_xz_0_0_xyzzz[k];

                g_xz_0_x_zzzz[k] = -g_z_0_0_zzzz[k] - g_xz_0_0_zzzz[k] * ab_x + g_xz_0_0_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xz_0_0_xxxx, g_xz_0_0_xxxxy, g_xz_0_0_xxxy, g_xz_0_0_xxxyy, g_xz_0_0_xxxyz, g_xz_0_0_xxxz, g_xz_0_0_xxyy, g_xz_0_0_xxyyy, g_xz_0_0_xxyyz, g_xz_0_0_xxyz, g_xz_0_0_xxyzz, g_xz_0_0_xxzz, g_xz_0_0_xyyy, g_xz_0_0_xyyyy, g_xz_0_0_xyyyz, g_xz_0_0_xyyz, g_xz_0_0_xyyzz, g_xz_0_0_xyzz, g_xz_0_0_xyzzz, g_xz_0_0_xzzz, g_xz_0_0_yyyy, g_xz_0_0_yyyyy, g_xz_0_0_yyyyz, g_xz_0_0_yyyz, g_xz_0_0_yyyzz, g_xz_0_0_yyzz, g_xz_0_0_yyzzz, g_xz_0_0_yzzz, g_xz_0_0_yzzzz, g_xz_0_0_zzzz, g_xz_0_y_xxxx, g_xz_0_y_xxxy, g_xz_0_y_xxxz, g_xz_0_y_xxyy, g_xz_0_y_xxyz, g_xz_0_y_xxzz, g_xz_0_y_xyyy, g_xz_0_y_xyyz, g_xz_0_y_xyzz, g_xz_0_y_xzzz, g_xz_0_y_yyyy, g_xz_0_y_yyyz, g_xz_0_y_yyzz, g_xz_0_y_yzzz, g_xz_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_y_xxxx[k] = -g_xz_0_0_xxxx[k] * ab_y + g_xz_0_0_xxxxy[k];

                g_xz_0_y_xxxy[k] = -g_xz_0_0_xxxy[k] * ab_y + g_xz_0_0_xxxyy[k];

                g_xz_0_y_xxxz[k] = -g_xz_0_0_xxxz[k] * ab_y + g_xz_0_0_xxxyz[k];

                g_xz_0_y_xxyy[k] = -g_xz_0_0_xxyy[k] * ab_y + g_xz_0_0_xxyyy[k];

                g_xz_0_y_xxyz[k] = -g_xz_0_0_xxyz[k] * ab_y + g_xz_0_0_xxyyz[k];

                g_xz_0_y_xxzz[k] = -g_xz_0_0_xxzz[k] * ab_y + g_xz_0_0_xxyzz[k];

                g_xz_0_y_xyyy[k] = -g_xz_0_0_xyyy[k] * ab_y + g_xz_0_0_xyyyy[k];

                g_xz_0_y_xyyz[k] = -g_xz_0_0_xyyz[k] * ab_y + g_xz_0_0_xyyyz[k];

                g_xz_0_y_xyzz[k] = -g_xz_0_0_xyzz[k] * ab_y + g_xz_0_0_xyyzz[k];

                g_xz_0_y_xzzz[k] = -g_xz_0_0_xzzz[k] * ab_y + g_xz_0_0_xyzzz[k];

                g_xz_0_y_yyyy[k] = -g_xz_0_0_yyyy[k] * ab_y + g_xz_0_0_yyyyy[k];

                g_xz_0_y_yyyz[k] = -g_xz_0_0_yyyz[k] * ab_y + g_xz_0_0_yyyyz[k];

                g_xz_0_y_yyzz[k] = -g_xz_0_0_yyzz[k] * ab_y + g_xz_0_0_yyyzz[k];

                g_xz_0_y_yzzz[k] = -g_xz_0_0_yzzz[k] * ab_y + g_xz_0_0_yyzzz[k];

                g_xz_0_y_zzzz[k] = -g_xz_0_0_zzzz[k] * ab_y + g_xz_0_0_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz, g_xz_0_0_xxxx, g_xz_0_0_xxxxz, g_xz_0_0_xxxy, g_xz_0_0_xxxyz, g_xz_0_0_xxxz, g_xz_0_0_xxxzz, g_xz_0_0_xxyy, g_xz_0_0_xxyyz, g_xz_0_0_xxyz, g_xz_0_0_xxyzz, g_xz_0_0_xxzz, g_xz_0_0_xxzzz, g_xz_0_0_xyyy, g_xz_0_0_xyyyz, g_xz_0_0_xyyz, g_xz_0_0_xyyzz, g_xz_0_0_xyzz, g_xz_0_0_xyzzz, g_xz_0_0_xzzz, g_xz_0_0_xzzzz, g_xz_0_0_yyyy, g_xz_0_0_yyyyz, g_xz_0_0_yyyz, g_xz_0_0_yyyzz, g_xz_0_0_yyzz, g_xz_0_0_yyzzz, g_xz_0_0_yzzz, g_xz_0_0_yzzzz, g_xz_0_0_zzzz, g_xz_0_0_zzzzz, g_xz_0_z_xxxx, g_xz_0_z_xxxy, g_xz_0_z_xxxz, g_xz_0_z_xxyy, g_xz_0_z_xxyz, g_xz_0_z_xxzz, g_xz_0_z_xyyy, g_xz_0_z_xyyz, g_xz_0_z_xyzz, g_xz_0_z_xzzz, g_xz_0_z_yyyy, g_xz_0_z_yyyz, g_xz_0_z_yyzz, g_xz_0_z_yzzz, g_xz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_z_xxxx[k] = -g_x_0_0_xxxx[k] - g_xz_0_0_xxxx[k] * ab_z + g_xz_0_0_xxxxz[k];

                g_xz_0_z_xxxy[k] = -g_x_0_0_xxxy[k] - g_xz_0_0_xxxy[k] * ab_z + g_xz_0_0_xxxyz[k];

                g_xz_0_z_xxxz[k] = -g_x_0_0_xxxz[k] - g_xz_0_0_xxxz[k] * ab_z + g_xz_0_0_xxxzz[k];

                g_xz_0_z_xxyy[k] = -g_x_0_0_xxyy[k] - g_xz_0_0_xxyy[k] * ab_z + g_xz_0_0_xxyyz[k];

                g_xz_0_z_xxyz[k] = -g_x_0_0_xxyz[k] - g_xz_0_0_xxyz[k] * ab_z + g_xz_0_0_xxyzz[k];

                g_xz_0_z_xxzz[k] = -g_x_0_0_xxzz[k] - g_xz_0_0_xxzz[k] * ab_z + g_xz_0_0_xxzzz[k];

                g_xz_0_z_xyyy[k] = -g_x_0_0_xyyy[k] - g_xz_0_0_xyyy[k] * ab_z + g_xz_0_0_xyyyz[k];

                g_xz_0_z_xyyz[k] = -g_x_0_0_xyyz[k] - g_xz_0_0_xyyz[k] * ab_z + g_xz_0_0_xyyzz[k];

                g_xz_0_z_xyzz[k] = -g_x_0_0_xyzz[k] - g_xz_0_0_xyzz[k] * ab_z + g_xz_0_0_xyzzz[k];

                g_xz_0_z_xzzz[k] = -g_x_0_0_xzzz[k] - g_xz_0_0_xzzz[k] * ab_z + g_xz_0_0_xzzzz[k];

                g_xz_0_z_yyyy[k] = -g_x_0_0_yyyy[k] - g_xz_0_0_yyyy[k] * ab_z + g_xz_0_0_yyyyz[k];

                g_xz_0_z_yyyz[k] = -g_x_0_0_yyyz[k] - g_xz_0_0_yyyz[k] * ab_z + g_xz_0_0_yyyzz[k];

                g_xz_0_z_yyzz[k] = -g_x_0_0_yyzz[k] - g_xz_0_0_yyzz[k] * ab_z + g_xz_0_0_yyzzz[k];

                g_xz_0_z_yzzz[k] = -g_x_0_0_yzzz[k] - g_xz_0_0_yzzz[k] * ab_z + g_xz_0_0_yzzzz[k];

                g_xz_0_z_zzzz[k] = -g_x_0_0_zzzz[k] - g_xz_0_0_zzzz[k] * ab_z + g_xz_0_0_zzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yy_0_0_xxxx, g_yy_0_0_xxxxx, g_yy_0_0_xxxxy, g_yy_0_0_xxxxz, g_yy_0_0_xxxy, g_yy_0_0_xxxyy, g_yy_0_0_xxxyz, g_yy_0_0_xxxz, g_yy_0_0_xxxzz, g_yy_0_0_xxyy, g_yy_0_0_xxyyy, g_yy_0_0_xxyyz, g_yy_0_0_xxyz, g_yy_0_0_xxyzz, g_yy_0_0_xxzz, g_yy_0_0_xxzzz, g_yy_0_0_xyyy, g_yy_0_0_xyyyy, g_yy_0_0_xyyyz, g_yy_0_0_xyyz, g_yy_0_0_xyyzz, g_yy_0_0_xyzz, g_yy_0_0_xyzzz, g_yy_0_0_xzzz, g_yy_0_0_xzzzz, g_yy_0_0_yyyy, g_yy_0_0_yyyz, g_yy_0_0_yyzz, g_yy_0_0_yzzz, g_yy_0_0_zzzz, g_yy_0_x_xxxx, g_yy_0_x_xxxy, g_yy_0_x_xxxz, g_yy_0_x_xxyy, g_yy_0_x_xxyz, g_yy_0_x_xxzz, g_yy_0_x_xyyy, g_yy_0_x_xyyz, g_yy_0_x_xyzz, g_yy_0_x_xzzz, g_yy_0_x_yyyy, g_yy_0_x_yyyz, g_yy_0_x_yyzz, g_yy_0_x_yzzz, g_yy_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_x_xxxx[k] = -g_yy_0_0_xxxx[k] * ab_x + g_yy_0_0_xxxxx[k];

                g_yy_0_x_xxxy[k] = -g_yy_0_0_xxxy[k] * ab_x + g_yy_0_0_xxxxy[k];

                g_yy_0_x_xxxz[k] = -g_yy_0_0_xxxz[k] * ab_x + g_yy_0_0_xxxxz[k];

                g_yy_0_x_xxyy[k] = -g_yy_0_0_xxyy[k] * ab_x + g_yy_0_0_xxxyy[k];

                g_yy_0_x_xxyz[k] = -g_yy_0_0_xxyz[k] * ab_x + g_yy_0_0_xxxyz[k];

                g_yy_0_x_xxzz[k] = -g_yy_0_0_xxzz[k] * ab_x + g_yy_0_0_xxxzz[k];

                g_yy_0_x_xyyy[k] = -g_yy_0_0_xyyy[k] * ab_x + g_yy_0_0_xxyyy[k];

                g_yy_0_x_xyyz[k] = -g_yy_0_0_xyyz[k] * ab_x + g_yy_0_0_xxyyz[k];

                g_yy_0_x_xyzz[k] = -g_yy_0_0_xyzz[k] * ab_x + g_yy_0_0_xxyzz[k];

                g_yy_0_x_xzzz[k] = -g_yy_0_0_xzzz[k] * ab_x + g_yy_0_0_xxzzz[k];

                g_yy_0_x_yyyy[k] = -g_yy_0_0_yyyy[k] * ab_x + g_yy_0_0_xyyyy[k];

                g_yy_0_x_yyyz[k] = -g_yy_0_0_yyyz[k] * ab_x + g_yy_0_0_xyyyz[k];

                g_yy_0_x_yyzz[k] = -g_yy_0_0_yyzz[k] * ab_x + g_yy_0_0_xyyzz[k];

                g_yy_0_x_yzzz[k] = -g_yy_0_0_yzzz[k] * ab_x + g_yy_0_0_xyzzz[k];

                g_yy_0_x_zzzz[k] = -g_yy_0_0_zzzz[k] * ab_x + g_yy_0_0_xzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz, g_yy_0_0_xxxx, g_yy_0_0_xxxxy, g_yy_0_0_xxxy, g_yy_0_0_xxxyy, g_yy_0_0_xxxyz, g_yy_0_0_xxxz, g_yy_0_0_xxyy, g_yy_0_0_xxyyy, g_yy_0_0_xxyyz, g_yy_0_0_xxyz, g_yy_0_0_xxyzz, g_yy_0_0_xxzz, g_yy_0_0_xyyy, g_yy_0_0_xyyyy, g_yy_0_0_xyyyz, g_yy_0_0_xyyz, g_yy_0_0_xyyzz, g_yy_0_0_xyzz, g_yy_0_0_xyzzz, g_yy_0_0_xzzz, g_yy_0_0_yyyy, g_yy_0_0_yyyyy, g_yy_0_0_yyyyz, g_yy_0_0_yyyz, g_yy_0_0_yyyzz, g_yy_0_0_yyzz, g_yy_0_0_yyzzz, g_yy_0_0_yzzz, g_yy_0_0_yzzzz, g_yy_0_0_zzzz, g_yy_0_y_xxxx, g_yy_0_y_xxxy, g_yy_0_y_xxxz, g_yy_0_y_xxyy, g_yy_0_y_xxyz, g_yy_0_y_xxzz, g_yy_0_y_xyyy, g_yy_0_y_xyyz, g_yy_0_y_xyzz, g_yy_0_y_xzzz, g_yy_0_y_yyyy, g_yy_0_y_yyyz, g_yy_0_y_yyzz, g_yy_0_y_yzzz, g_yy_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_y_xxxx[k] = -2.0 * g_y_0_0_xxxx[k] - g_yy_0_0_xxxx[k] * ab_y + g_yy_0_0_xxxxy[k];

                g_yy_0_y_xxxy[k] = -2.0 * g_y_0_0_xxxy[k] - g_yy_0_0_xxxy[k] * ab_y + g_yy_0_0_xxxyy[k];

                g_yy_0_y_xxxz[k] = -2.0 * g_y_0_0_xxxz[k] - g_yy_0_0_xxxz[k] * ab_y + g_yy_0_0_xxxyz[k];

                g_yy_0_y_xxyy[k] = -2.0 * g_y_0_0_xxyy[k] - g_yy_0_0_xxyy[k] * ab_y + g_yy_0_0_xxyyy[k];

                g_yy_0_y_xxyz[k] = -2.0 * g_y_0_0_xxyz[k] - g_yy_0_0_xxyz[k] * ab_y + g_yy_0_0_xxyyz[k];

                g_yy_0_y_xxzz[k] = -2.0 * g_y_0_0_xxzz[k] - g_yy_0_0_xxzz[k] * ab_y + g_yy_0_0_xxyzz[k];

                g_yy_0_y_xyyy[k] = -2.0 * g_y_0_0_xyyy[k] - g_yy_0_0_xyyy[k] * ab_y + g_yy_0_0_xyyyy[k];

                g_yy_0_y_xyyz[k] = -2.0 * g_y_0_0_xyyz[k] - g_yy_0_0_xyyz[k] * ab_y + g_yy_0_0_xyyyz[k];

                g_yy_0_y_xyzz[k] = -2.0 * g_y_0_0_xyzz[k] - g_yy_0_0_xyzz[k] * ab_y + g_yy_0_0_xyyzz[k];

                g_yy_0_y_xzzz[k] = -2.0 * g_y_0_0_xzzz[k] - g_yy_0_0_xzzz[k] * ab_y + g_yy_0_0_xyzzz[k];

                g_yy_0_y_yyyy[k] = -2.0 * g_y_0_0_yyyy[k] - g_yy_0_0_yyyy[k] * ab_y + g_yy_0_0_yyyyy[k];

                g_yy_0_y_yyyz[k] = -2.0 * g_y_0_0_yyyz[k] - g_yy_0_0_yyyz[k] * ab_y + g_yy_0_0_yyyyz[k];

                g_yy_0_y_yyzz[k] = -2.0 * g_y_0_0_yyzz[k] - g_yy_0_0_yyzz[k] * ab_y + g_yy_0_0_yyyzz[k];

                g_yy_0_y_yzzz[k] = -2.0 * g_y_0_0_yzzz[k] - g_yy_0_0_yzzz[k] * ab_y + g_yy_0_0_yyzzz[k];

                g_yy_0_y_zzzz[k] = -2.0 * g_y_0_0_zzzz[k] - g_yy_0_0_zzzz[k] * ab_y + g_yy_0_0_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yy_0_0_xxxx, g_yy_0_0_xxxxz, g_yy_0_0_xxxy, g_yy_0_0_xxxyz, g_yy_0_0_xxxz, g_yy_0_0_xxxzz, g_yy_0_0_xxyy, g_yy_0_0_xxyyz, g_yy_0_0_xxyz, g_yy_0_0_xxyzz, g_yy_0_0_xxzz, g_yy_0_0_xxzzz, g_yy_0_0_xyyy, g_yy_0_0_xyyyz, g_yy_0_0_xyyz, g_yy_0_0_xyyzz, g_yy_0_0_xyzz, g_yy_0_0_xyzzz, g_yy_0_0_xzzz, g_yy_0_0_xzzzz, g_yy_0_0_yyyy, g_yy_0_0_yyyyz, g_yy_0_0_yyyz, g_yy_0_0_yyyzz, g_yy_0_0_yyzz, g_yy_0_0_yyzzz, g_yy_0_0_yzzz, g_yy_0_0_yzzzz, g_yy_0_0_zzzz, g_yy_0_0_zzzzz, g_yy_0_z_xxxx, g_yy_0_z_xxxy, g_yy_0_z_xxxz, g_yy_0_z_xxyy, g_yy_0_z_xxyz, g_yy_0_z_xxzz, g_yy_0_z_xyyy, g_yy_0_z_xyyz, g_yy_0_z_xyzz, g_yy_0_z_xzzz, g_yy_0_z_yyyy, g_yy_0_z_yyyz, g_yy_0_z_yyzz, g_yy_0_z_yzzz, g_yy_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_z_xxxx[k] = -g_yy_0_0_xxxx[k] * ab_z + g_yy_0_0_xxxxz[k];

                g_yy_0_z_xxxy[k] = -g_yy_0_0_xxxy[k] * ab_z + g_yy_0_0_xxxyz[k];

                g_yy_0_z_xxxz[k] = -g_yy_0_0_xxxz[k] * ab_z + g_yy_0_0_xxxzz[k];

                g_yy_0_z_xxyy[k] = -g_yy_0_0_xxyy[k] * ab_z + g_yy_0_0_xxyyz[k];

                g_yy_0_z_xxyz[k] = -g_yy_0_0_xxyz[k] * ab_z + g_yy_0_0_xxyzz[k];

                g_yy_0_z_xxzz[k] = -g_yy_0_0_xxzz[k] * ab_z + g_yy_0_0_xxzzz[k];

                g_yy_0_z_xyyy[k] = -g_yy_0_0_xyyy[k] * ab_z + g_yy_0_0_xyyyz[k];

                g_yy_0_z_xyyz[k] = -g_yy_0_0_xyyz[k] * ab_z + g_yy_0_0_xyyzz[k];

                g_yy_0_z_xyzz[k] = -g_yy_0_0_xyzz[k] * ab_z + g_yy_0_0_xyzzz[k];

                g_yy_0_z_xzzz[k] = -g_yy_0_0_xzzz[k] * ab_z + g_yy_0_0_xzzzz[k];

                g_yy_0_z_yyyy[k] = -g_yy_0_0_yyyy[k] * ab_z + g_yy_0_0_yyyyz[k];

                g_yy_0_z_yyyz[k] = -g_yy_0_0_yyyz[k] * ab_z + g_yy_0_0_yyyzz[k];

                g_yy_0_z_yyzz[k] = -g_yy_0_0_yyzz[k] * ab_z + g_yy_0_0_yyzzz[k];

                g_yy_0_z_yzzz[k] = -g_yy_0_0_yzzz[k] * ab_z + g_yy_0_0_yzzzz[k];

                g_yy_0_z_zzzz[k] = -g_yy_0_0_zzzz[k] * ab_z + g_yy_0_0_zzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yz_0_0_xxxx, g_yz_0_0_xxxxx, g_yz_0_0_xxxxy, g_yz_0_0_xxxxz, g_yz_0_0_xxxy, g_yz_0_0_xxxyy, g_yz_0_0_xxxyz, g_yz_0_0_xxxz, g_yz_0_0_xxxzz, g_yz_0_0_xxyy, g_yz_0_0_xxyyy, g_yz_0_0_xxyyz, g_yz_0_0_xxyz, g_yz_0_0_xxyzz, g_yz_0_0_xxzz, g_yz_0_0_xxzzz, g_yz_0_0_xyyy, g_yz_0_0_xyyyy, g_yz_0_0_xyyyz, g_yz_0_0_xyyz, g_yz_0_0_xyyzz, g_yz_0_0_xyzz, g_yz_0_0_xyzzz, g_yz_0_0_xzzz, g_yz_0_0_xzzzz, g_yz_0_0_yyyy, g_yz_0_0_yyyz, g_yz_0_0_yyzz, g_yz_0_0_yzzz, g_yz_0_0_zzzz, g_yz_0_x_xxxx, g_yz_0_x_xxxy, g_yz_0_x_xxxz, g_yz_0_x_xxyy, g_yz_0_x_xxyz, g_yz_0_x_xxzz, g_yz_0_x_xyyy, g_yz_0_x_xyyz, g_yz_0_x_xyzz, g_yz_0_x_xzzz, g_yz_0_x_yyyy, g_yz_0_x_yyyz, g_yz_0_x_yyzz, g_yz_0_x_yzzz, g_yz_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_x_xxxx[k] = -g_yz_0_0_xxxx[k] * ab_x + g_yz_0_0_xxxxx[k];

                g_yz_0_x_xxxy[k] = -g_yz_0_0_xxxy[k] * ab_x + g_yz_0_0_xxxxy[k];

                g_yz_0_x_xxxz[k] = -g_yz_0_0_xxxz[k] * ab_x + g_yz_0_0_xxxxz[k];

                g_yz_0_x_xxyy[k] = -g_yz_0_0_xxyy[k] * ab_x + g_yz_0_0_xxxyy[k];

                g_yz_0_x_xxyz[k] = -g_yz_0_0_xxyz[k] * ab_x + g_yz_0_0_xxxyz[k];

                g_yz_0_x_xxzz[k] = -g_yz_0_0_xxzz[k] * ab_x + g_yz_0_0_xxxzz[k];

                g_yz_0_x_xyyy[k] = -g_yz_0_0_xyyy[k] * ab_x + g_yz_0_0_xxyyy[k];

                g_yz_0_x_xyyz[k] = -g_yz_0_0_xyyz[k] * ab_x + g_yz_0_0_xxyyz[k];

                g_yz_0_x_xyzz[k] = -g_yz_0_0_xyzz[k] * ab_x + g_yz_0_0_xxyzz[k];

                g_yz_0_x_xzzz[k] = -g_yz_0_0_xzzz[k] * ab_x + g_yz_0_0_xxzzz[k];

                g_yz_0_x_yyyy[k] = -g_yz_0_0_yyyy[k] * ab_x + g_yz_0_0_xyyyy[k];

                g_yz_0_x_yyyz[k] = -g_yz_0_0_yyyz[k] * ab_x + g_yz_0_0_xyyyz[k];

                g_yz_0_x_yyzz[k] = -g_yz_0_0_yyzz[k] * ab_x + g_yz_0_0_xyyzz[k];

                g_yz_0_x_yzzz[k] = -g_yz_0_0_yzzz[k] * ab_x + g_yz_0_0_xyzzz[k];

                g_yz_0_x_zzzz[k] = -g_yz_0_0_zzzz[k] * ab_x + g_yz_0_0_xzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yz_0_0_xxxx, g_yz_0_0_xxxxy, g_yz_0_0_xxxy, g_yz_0_0_xxxyy, g_yz_0_0_xxxyz, g_yz_0_0_xxxz, g_yz_0_0_xxyy, g_yz_0_0_xxyyy, g_yz_0_0_xxyyz, g_yz_0_0_xxyz, g_yz_0_0_xxyzz, g_yz_0_0_xxzz, g_yz_0_0_xyyy, g_yz_0_0_xyyyy, g_yz_0_0_xyyyz, g_yz_0_0_xyyz, g_yz_0_0_xyyzz, g_yz_0_0_xyzz, g_yz_0_0_xyzzz, g_yz_0_0_xzzz, g_yz_0_0_yyyy, g_yz_0_0_yyyyy, g_yz_0_0_yyyyz, g_yz_0_0_yyyz, g_yz_0_0_yyyzz, g_yz_0_0_yyzz, g_yz_0_0_yyzzz, g_yz_0_0_yzzz, g_yz_0_0_yzzzz, g_yz_0_0_zzzz, g_yz_0_y_xxxx, g_yz_0_y_xxxy, g_yz_0_y_xxxz, g_yz_0_y_xxyy, g_yz_0_y_xxyz, g_yz_0_y_xxzz, g_yz_0_y_xyyy, g_yz_0_y_xyyz, g_yz_0_y_xyzz, g_yz_0_y_xzzz, g_yz_0_y_yyyy, g_yz_0_y_yyyz, g_yz_0_y_yyzz, g_yz_0_y_yzzz, g_yz_0_y_zzzz, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_y_xxxx[k] = -g_z_0_0_xxxx[k] - g_yz_0_0_xxxx[k] * ab_y + g_yz_0_0_xxxxy[k];

                g_yz_0_y_xxxy[k] = -g_z_0_0_xxxy[k] - g_yz_0_0_xxxy[k] * ab_y + g_yz_0_0_xxxyy[k];

                g_yz_0_y_xxxz[k] = -g_z_0_0_xxxz[k] - g_yz_0_0_xxxz[k] * ab_y + g_yz_0_0_xxxyz[k];

                g_yz_0_y_xxyy[k] = -g_z_0_0_xxyy[k] - g_yz_0_0_xxyy[k] * ab_y + g_yz_0_0_xxyyy[k];

                g_yz_0_y_xxyz[k] = -g_z_0_0_xxyz[k] - g_yz_0_0_xxyz[k] * ab_y + g_yz_0_0_xxyyz[k];

                g_yz_0_y_xxzz[k] = -g_z_0_0_xxzz[k] - g_yz_0_0_xxzz[k] * ab_y + g_yz_0_0_xxyzz[k];

                g_yz_0_y_xyyy[k] = -g_z_0_0_xyyy[k] - g_yz_0_0_xyyy[k] * ab_y + g_yz_0_0_xyyyy[k];

                g_yz_0_y_xyyz[k] = -g_z_0_0_xyyz[k] - g_yz_0_0_xyyz[k] * ab_y + g_yz_0_0_xyyyz[k];

                g_yz_0_y_xyzz[k] = -g_z_0_0_xyzz[k] - g_yz_0_0_xyzz[k] * ab_y + g_yz_0_0_xyyzz[k];

                g_yz_0_y_xzzz[k] = -g_z_0_0_xzzz[k] - g_yz_0_0_xzzz[k] * ab_y + g_yz_0_0_xyzzz[k];

                g_yz_0_y_yyyy[k] = -g_z_0_0_yyyy[k] - g_yz_0_0_yyyy[k] * ab_y + g_yz_0_0_yyyyy[k];

                g_yz_0_y_yyyz[k] = -g_z_0_0_yyyz[k] - g_yz_0_0_yyyz[k] * ab_y + g_yz_0_0_yyyyz[k];

                g_yz_0_y_yyzz[k] = -g_z_0_0_yyzz[k] - g_yz_0_0_yyzz[k] * ab_y + g_yz_0_0_yyyzz[k];

                g_yz_0_y_yzzz[k] = -g_z_0_0_yzzz[k] - g_yz_0_0_yzzz[k] * ab_y + g_yz_0_0_yyzzz[k];

                g_yz_0_y_zzzz[k] = -g_z_0_0_zzzz[k] - g_yz_0_0_zzzz[k] * ab_y + g_yz_0_0_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz, g_yz_0_0_xxxx, g_yz_0_0_xxxxz, g_yz_0_0_xxxy, g_yz_0_0_xxxyz, g_yz_0_0_xxxz, g_yz_0_0_xxxzz, g_yz_0_0_xxyy, g_yz_0_0_xxyyz, g_yz_0_0_xxyz, g_yz_0_0_xxyzz, g_yz_0_0_xxzz, g_yz_0_0_xxzzz, g_yz_0_0_xyyy, g_yz_0_0_xyyyz, g_yz_0_0_xyyz, g_yz_0_0_xyyzz, g_yz_0_0_xyzz, g_yz_0_0_xyzzz, g_yz_0_0_xzzz, g_yz_0_0_xzzzz, g_yz_0_0_yyyy, g_yz_0_0_yyyyz, g_yz_0_0_yyyz, g_yz_0_0_yyyzz, g_yz_0_0_yyzz, g_yz_0_0_yyzzz, g_yz_0_0_yzzz, g_yz_0_0_yzzzz, g_yz_0_0_zzzz, g_yz_0_0_zzzzz, g_yz_0_z_xxxx, g_yz_0_z_xxxy, g_yz_0_z_xxxz, g_yz_0_z_xxyy, g_yz_0_z_xxyz, g_yz_0_z_xxzz, g_yz_0_z_xyyy, g_yz_0_z_xyyz, g_yz_0_z_xyzz, g_yz_0_z_xzzz, g_yz_0_z_yyyy, g_yz_0_z_yyyz, g_yz_0_z_yyzz, g_yz_0_z_yzzz, g_yz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_z_xxxx[k] = -g_y_0_0_xxxx[k] - g_yz_0_0_xxxx[k] * ab_z + g_yz_0_0_xxxxz[k];

                g_yz_0_z_xxxy[k] = -g_y_0_0_xxxy[k] - g_yz_0_0_xxxy[k] * ab_z + g_yz_0_0_xxxyz[k];

                g_yz_0_z_xxxz[k] = -g_y_0_0_xxxz[k] - g_yz_0_0_xxxz[k] * ab_z + g_yz_0_0_xxxzz[k];

                g_yz_0_z_xxyy[k] = -g_y_0_0_xxyy[k] - g_yz_0_0_xxyy[k] * ab_z + g_yz_0_0_xxyyz[k];

                g_yz_0_z_xxyz[k] = -g_y_0_0_xxyz[k] - g_yz_0_0_xxyz[k] * ab_z + g_yz_0_0_xxyzz[k];

                g_yz_0_z_xxzz[k] = -g_y_0_0_xxzz[k] - g_yz_0_0_xxzz[k] * ab_z + g_yz_0_0_xxzzz[k];

                g_yz_0_z_xyyy[k] = -g_y_0_0_xyyy[k] - g_yz_0_0_xyyy[k] * ab_z + g_yz_0_0_xyyyz[k];

                g_yz_0_z_xyyz[k] = -g_y_0_0_xyyz[k] - g_yz_0_0_xyyz[k] * ab_z + g_yz_0_0_xyyzz[k];

                g_yz_0_z_xyzz[k] = -g_y_0_0_xyzz[k] - g_yz_0_0_xyzz[k] * ab_z + g_yz_0_0_xyzzz[k];

                g_yz_0_z_xzzz[k] = -g_y_0_0_xzzz[k] - g_yz_0_0_xzzz[k] * ab_z + g_yz_0_0_xzzzz[k];

                g_yz_0_z_yyyy[k] = -g_y_0_0_yyyy[k] - g_yz_0_0_yyyy[k] * ab_z + g_yz_0_0_yyyyz[k];

                g_yz_0_z_yyyz[k] = -g_y_0_0_yyyz[k] - g_yz_0_0_yyyz[k] * ab_z + g_yz_0_0_yyyzz[k];

                g_yz_0_z_yyzz[k] = -g_y_0_0_yyzz[k] - g_yz_0_0_yyzz[k] * ab_z + g_yz_0_0_yyzzz[k];

                g_yz_0_z_yzzz[k] = -g_y_0_0_yzzz[k] - g_yz_0_0_yzzz[k] * ab_z + g_yz_0_0_yzzzz[k];

                g_yz_0_z_zzzz[k] = -g_y_0_0_zzzz[k] - g_yz_0_0_zzzz[k] * ab_z + g_yz_0_0_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_zz_0_0_xxxx, g_zz_0_0_xxxxx, g_zz_0_0_xxxxy, g_zz_0_0_xxxxz, g_zz_0_0_xxxy, g_zz_0_0_xxxyy, g_zz_0_0_xxxyz, g_zz_0_0_xxxz, g_zz_0_0_xxxzz, g_zz_0_0_xxyy, g_zz_0_0_xxyyy, g_zz_0_0_xxyyz, g_zz_0_0_xxyz, g_zz_0_0_xxyzz, g_zz_0_0_xxzz, g_zz_0_0_xxzzz, g_zz_0_0_xyyy, g_zz_0_0_xyyyy, g_zz_0_0_xyyyz, g_zz_0_0_xyyz, g_zz_0_0_xyyzz, g_zz_0_0_xyzz, g_zz_0_0_xyzzz, g_zz_0_0_xzzz, g_zz_0_0_xzzzz, g_zz_0_0_yyyy, g_zz_0_0_yyyz, g_zz_0_0_yyzz, g_zz_0_0_yzzz, g_zz_0_0_zzzz, g_zz_0_x_xxxx, g_zz_0_x_xxxy, g_zz_0_x_xxxz, g_zz_0_x_xxyy, g_zz_0_x_xxyz, g_zz_0_x_xxzz, g_zz_0_x_xyyy, g_zz_0_x_xyyz, g_zz_0_x_xyzz, g_zz_0_x_xzzz, g_zz_0_x_yyyy, g_zz_0_x_yyyz, g_zz_0_x_yyzz, g_zz_0_x_yzzz, g_zz_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_x_xxxx[k] = -g_zz_0_0_xxxx[k] * ab_x + g_zz_0_0_xxxxx[k];

                g_zz_0_x_xxxy[k] = -g_zz_0_0_xxxy[k] * ab_x + g_zz_0_0_xxxxy[k];

                g_zz_0_x_xxxz[k] = -g_zz_0_0_xxxz[k] * ab_x + g_zz_0_0_xxxxz[k];

                g_zz_0_x_xxyy[k] = -g_zz_0_0_xxyy[k] * ab_x + g_zz_0_0_xxxyy[k];

                g_zz_0_x_xxyz[k] = -g_zz_0_0_xxyz[k] * ab_x + g_zz_0_0_xxxyz[k];

                g_zz_0_x_xxzz[k] = -g_zz_0_0_xxzz[k] * ab_x + g_zz_0_0_xxxzz[k];

                g_zz_0_x_xyyy[k] = -g_zz_0_0_xyyy[k] * ab_x + g_zz_0_0_xxyyy[k];

                g_zz_0_x_xyyz[k] = -g_zz_0_0_xyyz[k] * ab_x + g_zz_0_0_xxyyz[k];

                g_zz_0_x_xyzz[k] = -g_zz_0_0_xyzz[k] * ab_x + g_zz_0_0_xxyzz[k];

                g_zz_0_x_xzzz[k] = -g_zz_0_0_xzzz[k] * ab_x + g_zz_0_0_xxzzz[k];

                g_zz_0_x_yyyy[k] = -g_zz_0_0_yyyy[k] * ab_x + g_zz_0_0_xyyyy[k];

                g_zz_0_x_yyyz[k] = -g_zz_0_0_yyyz[k] * ab_x + g_zz_0_0_xyyyz[k];

                g_zz_0_x_yyzz[k] = -g_zz_0_0_yyzz[k] * ab_x + g_zz_0_0_xyyzz[k];

                g_zz_0_x_yzzz[k] = -g_zz_0_0_yzzz[k] * ab_x + g_zz_0_0_xyzzz[k];

                g_zz_0_x_zzzz[k] = -g_zz_0_0_zzzz[k] * ab_x + g_zz_0_0_xzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_zz_0_0_xxxx, g_zz_0_0_xxxxy, g_zz_0_0_xxxy, g_zz_0_0_xxxyy, g_zz_0_0_xxxyz, g_zz_0_0_xxxz, g_zz_0_0_xxyy, g_zz_0_0_xxyyy, g_zz_0_0_xxyyz, g_zz_0_0_xxyz, g_zz_0_0_xxyzz, g_zz_0_0_xxzz, g_zz_0_0_xyyy, g_zz_0_0_xyyyy, g_zz_0_0_xyyyz, g_zz_0_0_xyyz, g_zz_0_0_xyyzz, g_zz_0_0_xyzz, g_zz_0_0_xyzzz, g_zz_0_0_xzzz, g_zz_0_0_yyyy, g_zz_0_0_yyyyy, g_zz_0_0_yyyyz, g_zz_0_0_yyyz, g_zz_0_0_yyyzz, g_zz_0_0_yyzz, g_zz_0_0_yyzzz, g_zz_0_0_yzzz, g_zz_0_0_yzzzz, g_zz_0_0_zzzz, g_zz_0_y_xxxx, g_zz_0_y_xxxy, g_zz_0_y_xxxz, g_zz_0_y_xxyy, g_zz_0_y_xxyz, g_zz_0_y_xxzz, g_zz_0_y_xyyy, g_zz_0_y_xyyz, g_zz_0_y_xyzz, g_zz_0_y_xzzz, g_zz_0_y_yyyy, g_zz_0_y_yyyz, g_zz_0_y_yyzz, g_zz_0_y_yzzz, g_zz_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_y_xxxx[k] = -g_zz_0_0_xxxx[k] * ab_y + g_zz_0_0_xxxxy[k];

                g_zz_0_y_xxxy[k] = -g_zz_0_0_xxxy[k] * ab_y + g_zz_0_0_xxxyy[k];

                g_zz_0_y_xxxz[k] = -g_zz_0_0_xxxz[k] * ab_y + g_zz_0_0_xxxyz[k];

                g_zz_0_y_xxyy[k] = -g_zz_0_0_xxyy[k] * ab_y + g_zz_0_0_xxyyy[k];

                g_zz_0_y_xxyz[k] = -g_zz_0_0_xxyz[k] * ab_y + g_zz_0_0_xxyyz[k];

                g_zz_0_y_xxzz[k] = -g_zz_0_0_xxzz[k] * ab_y + g_zz_0_0_xxyzz[k];

                g_zz_0_y_xyyy[k] = -g_zz_0_0_xyyy[k] * ab_y + g_zz_0_0_xyyyy[k];

                g_zz_0_y_xyyz[k] = -g_zz_0_0_xyyz[k] * ab_y + g_zz_0_0_xyyyz[k];

                g_zz_0_y_xyzz[k] = -g_zz_0_0_xyzz[k] * ab_y + g_zz_0_0_xyyzz[k];

                g_zz_0_y_xzzz[k] = -g_zz_0_0_xzzz[k] * ab_y + g_zz_0_0_xyzzz[k];

                g_zz_0_y_yyyy[k] = -g_zz_0_0_yyyy[k] * ab_y + g_zz_0_0_yyyyy[k];

                g_zz_0_y_yyyz[k] = -g_zz_0_0_yyyz[k] * ab_y + g_zz_0_0_yyyyz[k];

                g_zz_0_y_yyzz[k] = -g_zz_0_0_yyzz[k] * ab_y + g_zz_0_0_yyyzz[k];

                g_zz_0_y_yzzz[k] = -g_zz_0_0_yzzz[k] * ab_y + g_zz_0_0_yyzzz[k];

                g_zz_0_y_zzzz[k] = -g_zz_0_0_zzzz[k] * ab_y + g_zz_0_0_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz, g_zz_0_0_xxxx, g_zz_0_0_xxxxz, g_zz_0_0_xxxy, g_zz_0_0_xxxyz, g_zz_0_0_xxxz, g_zz_0_0_xxxzz, g_zz_0_0_xxyy, g_zz_0_0_xxyyz, g_zz_0_0_xxyz, g_zz_0_0_xxyzz, g_zz_0_0_xxzz, g_zz_0_0_xxzzz, g_zz_0_0_xyyy, g_zz_0_0_xyyyz, g_zz_0_0_xyyz, g_zz_0_0_xyyzz, g_zz_0_0_xyzz, g_zz_0_0_xyzzz, g_zz_0_0_xzzz, g_zz_0_0_xzzzz, g_zz_0_0_yyyy, g_zz_0_0_yyyyz, g_zz_0_0_yyyz, g_zz_0_0_yyyzz, g_zz_0_0_yyzz, g_zz_0_0_yyzzz, g_zz_0_0_yzzz, g_zz_0_0_yzzzz, g_zz_0_0_zzzz, g_zz_0_0_zzzzz, g_zz_0_z_xxxx, g_zz_0_z_xxxy, g_zz_0_z_xxxz, g_zz_0_z_xxyy, g_zz_0_z_xxyz, g_zz_0_z_xxzz, g_zz_0_z_xyyy, g_zz_0_z_xyyz, g_zz_0_z_xyzz, g_zz_0_z_xzzz, g_zz_0_z_yyyy, g_zz_0_z_yyyz, g_zz_0_z_yyzz, g_zz_0_z_yzzz, g_zz_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_z_xxxx[k] = -2.0 * g_z_0_0_xxxx[k] - g_zz_0_0_xxxx[k] * ab_z + g_zz_0_0_xxxxz[k];

                g_zz_0_z_xxxy[k] = -2.0 * g_z_0_0_xxxy[k] - g_zz_0_0_xxxy[k] * ab_z + g_zz_0_0_xxxyz[k];

                g_zz_0_z_xxxz[k] = -2.0 * g_z_0_0_xxxz[k] - g_zz_0_0_xxxz[k] * ab_z + g_zz_0_0_xxxzz[k];

                g_zz_0_z_xxyy[k] = -2.0 * g_z_0_0_xxyy[k] - g_zz_0_0_xxyy[k] * ab_z + g_zz_0_0_xxyyz[k];

                g_zz_0_z_xxyz[k] = -2.0 * g_z_0_0_xxyz[k] - g_zz_0_0_xxyz[k] * ab_z + g_zz_0_0_xxyzz[k];

                g_zz_0_z_xxzz[k] = -2.0 * g_z_0_0_xxzz[k] - g_zz_0_0_xxzz[k] * ab_z + g_zz_0_0_xxzzz[k];

                g_zz_0_z_xyyy[k] = -2.0 * g_z_0_0_xyyy[k] - g_zz_0_0_xyyy[k] * ab_z + g_zz_0_0_xyyyz[k];

                g_zz_0_z_xyyz[k] = -2.0 * g_z_0_0_xyyz[k] - g_zz_0_0_xyyz[k] * ab_z + g_zz_0_0_xyyzz[k];

                g_zz_0_z_xyzz[k] = -2.0 * g_z_0_0_xyzz[k] - g_zz_0_0_xyzz[k] * ab_z + g_zz_0_0_xyzzz[k];

                g_zz_0_z_xzzz[k] = -2.0 * g_z_0_0_xzzz[k] - g_zz_0_0_xzzz[k] * ab_z + g_zz_0_0_xzzzz[k];

                g_zz_0_z_yyyy[k] = -2.0 * g_z_0_0_yyyy[k] - g_zz_0_0_yyyy[k] * ab_z + g_zz_0_0_yyyyz[k];

                g_zz_0_z_yyyz[k] = -2.0 * g_z_0_0_yyyz[k] - g_zz_0_0_yyyz[k] * ab_z + g_zz_0_0_yyyzz[k];

                g_zz_0_z_yyzz[k] = -2.0 * g_z_0_0_yyzz[k] - g_zz_0_0_yyzz[k] * ab_z + g_zz_0_0_yyzzz[k];

                g_zz_0_z_yzzz[k] = -2.0 * g_z_0_0_yzzz[k] - g_zz_0_0_yzzz[k] * ab_z + g_zz_0_0_yzzzz[k];

                g_zz_0_z_zzzz[k] = -2.0 * g_z_0_0_zzzz[k] - g_zz_0_0_zzzz[k] * ab_z + g_zz_0_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

