#include "ElectronRepulsionGeom2000ContrRecPHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_phxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_phxx,
                                            const size_t idx_geom_10_shxx,
                                            const size_t idx_geom_20_shxx,
                                            const size_t idx_geom_20_sixx,
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
            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_10_off = idx_geom_10_shxx + i * dcomps + j;

            auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 53 * ccomps * dcomps);

            auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 54 * ccomps * dcomps);

            auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 62 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_20_off = idx_geom_20_sixx + i * dcomps + j;

            auto g_xx_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 55 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 56 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 57 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 58 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 62 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 65 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 68 * ccomps * dcomps);

            auto g_xz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 69 * ccomps * dcomps);

            auto g_xz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 83 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 84 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 85 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 86 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 87 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 88 * ccomps * dcomps);

            auto g_yy_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 89 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 92 * ccomps * dcomps);

            auto g_yy_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 95 * ccomps * dcomps);

            auto g_yy_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 98 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 99 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 101 * ccomps * dcomps);

            auto g_yy_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 104 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 111 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 112 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 113 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 114 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 115 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 116 * ccomps * dcomps);

            auto g_yz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 117 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 118 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 119 * ccomps * dcomps);

            auto g_yz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 122 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 125 * ccomps * dcomps);

            auto g_yz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 128 * ccomps * dcomps);

            auto g_yz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 129 * ccomps * dcomps);

            auto g_yz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 131 * ccomps * dcomps);

            auto g_yz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 134 * ccomps * dcomps);

            auto g_yz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 137 * ccomps * dcomps);

            auto g_yz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 139 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 140 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 141 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 142 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 143 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 144 * ccomps * dcomps);

            auto g_zz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 145 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 146 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 147 * ccomps * dcomps);

            auto g_zz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 148 * ccomps * dcomps);

            auto g_zz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 149 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 152 * ccomps * dcomps);

            auto g_zz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 155 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 158 * ccomps * dcomps);

            auto g_zz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 159 * ccomps * dcomps);

            auto g_zz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 161 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 164 * ccomps * dcomps);

            auto g_zz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 167 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_phxx

            const auto ph_geom_20_off = idx_geom_20_phxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_xx_0_0_xxxxx, g_xx_0_0_xxxxxx, g_xx_0_0_xxxxxy, g_xx_0_0_xxxxxz, g_xx_0_0_xxxxy, g_xx_0_0_xxxxyy, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxz, g_xx_0_0_xxxxzz, g_xx_0_0_xxxyy, g_xx_0_0_xxxyyy, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxzz, g_xx_0_0_xxxzzz, g_xx_0_0_xxyyy, g_xx_0_0_xxyyyy, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxzzz, g_xx_0_0_xxzzzz, g_xx_0_0_xyyyy, g_xx_0_0_xyyyyy, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xzzzz, g_xx_0_0_xzzzzz, g_xx_0_0_yyyyy, g_xx_0_0_yyyyz, g_xx_0_0_yyyzz, g_xx_0_0_yyzzz, g_xx_0_0_yzzzz, g_xx_0_0_zzzzz, g_xx_0_x_xxxxx, g_xx_0_x_xxxxy, g_xx_0_x_xxxxz, g_xx_0_x_xxxyy, g_xx_0_x_xxxyz, g_xx_0_x_xxxzz, g_xx_0_x_xxyyy, g_xx_0_x_xxyyz, g_xx_0_x_xxyzz, g_xx_0_x_xxzzz, g_xx_0_x_xyyyy, g_xx_0_x_xyyyz, g_xx_0_x_xyyzz, g_xx_0_x_xyzzz, g_xx_0_x_xzzzz, g_xx_0_x_yyyyy, g_xx_0_x_yyyyz, g_xx_0_x_yyyzz, g_xx_0_x_yyzzz, g_xx_0_x_yzzzz, g_xx_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_x_xxxxx[k] = -2.0 * g_x_0_0_xxxxx[k] - g_xx_0_0_xxxxx[k] * ab_x + g_xx_0_0_xxxxxx[k];

                g_xx_0_x_xxxxy[k] = -2.0 * g_x_0_0_xxxxy[k] - g_xx_0_0_xxxxy[k] * ab_x + g_xx_0_0_xxxxxy[k];

                g_xx_0_x_xxxxz[k] = -2.0 * g_x_0_0_xxxxz[k] - g_xx_0_0_xxxxz[k] * ab_x + g_xx_0_0_xxxxxz[k];

                g_xx_0_x_xxxyy[k] = -2.0 * g_x_0_0_xxxyy[k] - g_xx_0_0_xxxyy[k] * ab_x + g_xx_0_0_xxxxyy[k];

                g_xx_0_x_xxxyz[k] = -2.0 * g_x_0_0_xxxyz[k] - g_xx_0_0_xxxyz[k] * ab_x + g_xx_0_0_xxxxyz[k];

                g_xx_0_x_xxxzz[k] = -2.0 * g_x_0_0_xxxzz[k] - g_xx_0_0_xxxzz[k] * ab_x + g_xx_0_0_xxxxzz[k];

                g_xx_0_x_xxyyy[k] = -2.0 * g_x_0_0_xxyyy[k] - g_xx_0_0_xxyyy[k] * ab_x + g_xx_0_0_xxxyyy[k];

                g_xx_0_x_xxyyz[k] = -2.0 * g_x_0_0_xxyyz[k] - g_xx_0_0_xxyyz[k] * ab_x + g_xx_0_0_xxxyyz[k];

                g_xx_0_x_xxyzz[k] = -2.0 * g_x_0_0_xxyzz[k] - g_xx_0_0_xxyzz[k] * ab_x + g_xx_0_0_xxxyzz[k];

                g_xx_0_x_xxzzz[k] = -2.0 * g_x_0_0_xxzzz[k] - g_xx_0_0_xxzzz[k] * ab_x + g_xx_0_0_xxxzzz[k];

                g_xx_0_x_xyyyy[k] = -2.0 * g_x_0_0_xyyyy[k] - g_xx_0_0_xyyyy[k] * ab_x + g_xx_0_0_xxyyyy[k];

                g_xx_0_x_xyyyz[k] = -2.0 * g_x_0_0_xyyyz[k] - g_xx_0_0_xyyyz[k] * ab_x + g_xx_0_0_xxyyyz[k];

                g_xx_0_x_xyyzz[k] = -2.0 * g_x_0_0_xyyzz[k] - g_xx_0_0_xyyzz[k] * ab_x + g_xx_0_0_xxyyzz[k];

                g_xx_0_x_xyzzz[k] = -2.0 * g_x_0_0_xyzzz[k] - g_xx_0_0_xyzzz[k] * ab_x + g_xx_0_0_xxyzzz[k];

                g_xx_0_x_xzzzz[k] = -2.0 * g_x_0_0_xzzzz[k] - g_xx_0_0_xzzzz[k] * ab_x + g_xx_0_0_xxzzzz[k];

                g_xx_0_x_yyyyy[k] = -2.0 * g_x_0_0_yyyyy[k] - g_xx_0_0_yyyyy[k] * ab_x + g_xx_0_0_xyyyyy[k];

                g_xx_0_x_yyyyz[k] = -2.0 * g_x_0_0_yyyyz[k] - g_xx_0_0_yyyyz[k] * ab_x + g_xx_0_0_xyyyyz[k];

                g_xx_0_x_yyyzz[k] = -2.0 * g_x_0_0_yyyzz[k] - g_xx_0_0_yyyzz[k] * ab_x + g_xx_0_0_xyyyzz[k];

                g_xx_0_x_yyzzz[k] = -2.0 * g_x_0_0_yyzzz[k] - g_xx_0_0_yyzzz[k] * ab_x + g_xx_0_0_xyyzzz[k];

                g_xx_0_x_yzzzz[k] = -2.0 * g_x_0_0_yzzzz[k] - g_xx_0_0_yzzzz[k] * ab_x + g_xx_0_0_xyzzzz[k];

                g_xx_0_x_zzzzz[k] = -2.0 * g_x_0_0_zzzzz[k] - g_xx_0_0_zzzzz[k] * ab_x + g_xx_0_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xx_0_0_xxxxx, g_xx_0_0_xxxxxy, g_xx_0_0_xxxxy, g_xx_0_0_xxxxyy, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxz, g_xx_0_0_xxxyy, g_xx_0_0_xxxyyy, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxzz, g_xx_0_0_xxyyy, g_xx_0_0_xxyyyy, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxzzz, g_xx_0_0_xyyyy, g_xx_0_0_xyyyyy, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xzzzz, g_xx_0_0_yyyyy, g_xx_0_0_yyyyyy, g_xx_0_0_yyyyyz, g_xx_0_0_yyyyz, g_xx_0_0_yyyyzz, g_xx_0_0_yyyzz, g_xx_0_0_yyyzzz, g_xx_0_0_yyzzz, g_xx_0_0_yyzzzz, g_xx_0_0_yzzzz, g_xx_0_0_yzzzzz, g_xx_0_0_zzzzz, g_xx_0_y_xxxxx, g_xx_0_y_xxxxy, g_xx_0_y_xxxxz, g_xx_0_y_xxxyy, g_xx_0_y_xxxyz, g_xx_0_y_xxxzz, g_xx_0_y_xxyyy, g_xx_0_y_xxyyz, g_xx_0_y_xxyzz, g_xx_0_y_xxzzz, g_xx_0_y_xyyyy, g_xx_0_y_xyyyz, g_xx_0_y_xyyzz, g_xx_0_y_xyzzz, g_xx_0_y_xzzzz, g_xx_0_y_yyyyy, g_xx_0_y_yyyyz, g_xx_0_y_yyyzz, g_xx_0_y_yyzzz, g_xx_0_y_yzzzz, g_xx_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_y_xxxxx[k] = -g_xx_0_0_xxxxx[k] * ab_y + g_xx_0_0_xxxxxy[k];

                g_xx_0_y_xxxxy[k] = -g_xx_0_0_xxxxy[k] * ab_y + g_xx_0_0_xxxxyy[k];

                g_xx_0_y_xxxxz[k] = -g_xx_0_0_xxxxz[k] * ab_y + g_xx_0_0_xxxxyz[k];

                g_xx_0_y_xxxyy[k] = -g_xx_0_0_xxxyy[k] * ab_y + g_xx_0_0_xxxyyy[k];

                g_xx_0_y_xxxyz[k] = -g_xx_0_0_xxxyz[k] * ab_y + g_xx_0_0_xxxyyz[k];

                g_xx_0_y_xxxzz[k] = -g_xx_0_0_xxxzz[k] * ab_y + g_xx_0_0_xxxyzz[k];

                g_xx_0_y_xxyyy[k] = -g_xx_0_0_xxyyy[k] * ab_y + g_xx_0_0_xxyyyy[k];

                g_xx_0_y_xxyyz[k] = -g_xx_0_0_xxyyz[k] * ab_y + g_xx_0_0_xxyyyz[k];

                g_xx_0_y_xxyzz[k] = -g_xx_0_0_xxyzz[k] * ab_y + g_xx_0_0_xxyyzz[k];

                g_xx_0_y_xxzzz[k] = -g_xx_0_0_xxzzz[k] * ab_y + g_xx_0_0_xxyzzz[k];

                g_xx_0_y_xyyyy[k] = -g_xx_0_0_xyyyy[k] * ab_y + g_xx_0_0_xyyyyy[k];

                g_xx_0_y_xyyyz[k] = -g_xx_0_0_xyyyz[k] * ab_y + g_xx_0_0_xyyyyz[k];

                g_xx_0_y_xyyzz[k] = -g_xx_0_0_xyyzz[k] * ab_y + g_xx_0_0_xyyyzz[k];

                g_xx_0_y_xyzzz[k] = -g_xx_0_0_xyzzz[k] * ab_y + g_xx_0_0_xyyzzz[k];

                g_xx_0_y_xzzzz[k] = -g_xx_0_0_xzzzz[k] * ab_y + g_xx_0_0_xyzzzz[k];

                g_xx_0_y_yyyyy[k] = -g_xx_0_0_yyyyy[k] * ab_y + g_xx_0_0_yyyyyy[k];

                g_xx_0_y_yyyyz[k] = -g_xx_0_0_yyyyz[k] * ab_y + g_xx_0_0_yyyyyz[k];

                g_xx_0_y_yyyzz[k] = -g_xx_0_0_yyyzz[k] * ab_y + g_xx_0_0_yyyyzz[k];

                g_xx_0_y_yyzzz[k] = -g_xx_0_0_yyzzz[k] * ab_y + g_xx_0_0_yyyzzz[k];

                g_xx_0_y_yzzzz[k] = -g_xx_0_0_yzzzz[k] * ab_y + g_xx_0_0_yyzzzz[k];

                g_xx_0_y_zzzzz[k] = -g_xx_0_0_zzzzz[k] * ab_y + g_xx_0_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xx_0_0_xxxxx, g_xx_0_0_xxxxxz, g_xx_0_0_xxxxy, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxz, g_xx_0_0_xxxxzz, g_xx_0_0_xxxyy, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxzz, g_xx_0_0_xxxzzz, g_xx_0_0_xxyyy, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxzzz, g_xx_0_0_xxzzzz, g_xx_0_0_xyyyy, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xzzzz, g_xx_0_0_xzzzzz, g_xx_0_0_yyyyy, g_xx_0_0_yyyyyz, g_xx_0_0_yyyyz, g_xx_0_0_yyyyzz, g_xx_0_0_yyyzz, g_xx_0_0_yyyzzz, g_xx_0_0_yyzzz, g_xx_0_0_yyzzzz, g_xx_0_0_yzzzz, g_xx_0_0_yzzzzz, g_xx_0_0_zzzzz, g_xx_0_0_zzzzzz, g_xx_0_z_xxxxx, g_xx_0_z_xxxxy, g_xx_0_z_xxxxz, g_xx_0_z_xxxyy, g_xx_0_z_xxxyz, g_xx_0_z_xxxzz, g_xx_0_z_xxyyy, g_xx_0_z_xxyyz, g_xx_0_z_xxyzz, g_xx_0_z_xxzzz, g_xx_0_z_xyyyy, g_xx_0_z_xyyyz, g_xx_0_z_xyyzz, g_xx_0_z_xyzzz, g_xx_0_z_xzzzz, g_xx_0_z_yyyyy, g_xx_0_z_yyyyz, g_xx_0_z_yyyzz, g_xx_0_z_yyzzz, g_xx_0_z_yzzzz, g_xx_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_z_xxxxx[k] = -g_xx_0_0_xxxxx[k] * ab_z + g_xx_0_0_xxxxxz[k];

                g_xx_0_z_xxxxy[k] = -g_xx_0_0_xxxxy[k] * ab_z + g_xx_0_0_xxxxyz[k];

                g_xx_0_z_xxxxz[k] = -g_xx_0_0_xxxxz[k] * ab_z + g_xx_0_0_xxxxzz[k];

                g_xx_0_z_xxxyy[k] = -g_xx_0_0_xxxyy[k] * ab_z + g_xx_0_0_xxxyyz[k];

                g_xx_0_z_xxxyz[k] = -g_xx_0_0_xxxyz[k] * ab_z + g_xx_0_0_xxxyzz[k];

                g_xx_0_z_xxxzz[k] = -g_xx_0_0_xxxzz[k] * ab_z + g_xx_0_0_xxxzzz[k];

                g_xx_0_z_xxyyy[k] = -g_xx_0_0_xxyyy[k] * ab_z + g_xx_0_0_xxyyyz[k];

                g_xx_0_z_xxyyz[k] = -g_xx_0_0_xxyyz[k] * ab_z + g_xx_0_0_xxyyzz[k];

                g_xx_0_z_xxyzz[k] = -g_xx_0_0_xxyzz[k] * ab_z + g_xx_0_0_xxyzzz[k];

                g_xx_0_z_xxzzz[k] = -g_xx_0_0_xxzzz[k] * ab_z + g_xx_0_0_xxzzzz[k];

                g_xx_0_z_xyyyy[k] = -g_xx_0_0_xyyyy[k] * ab_z + g_xx_0_0_xyyyyz[k];

                g_xx_0_z_xyyyz[k] = -g_xx_0_0_xyyyz[k] * ab_z + g_xx_0_0_xyyyzz[k];

                g_xx_0_z_xyyzz[k] = -g_xx_0_0_xyyzz[k] * ab_z + g_xx_0_0_xyyzzz[k];

                g_xx_0_z_xyzzz[k] = -g_xx_0_0_xyzzz[k] * ab_z + g_xx_0_0_xyzzzz[k];

                g_xx_0_z_xzzzz[k] = -g_xx_0_0_xzzzz[k] * ab_z + g_xx_0_0_xzzzzz[k];

                g_xx_0_z_yyyyy[k] = -g_xx_0_0_yyyyy[k] * ab_z + g_xx_0_0_yyyyyz[k];

                g_xx_0_z_yyyyz[k] = -g_xx_0_0_yyyyz[k] * ab_z + g_xx_0_0_yyyyzz[k];

                g_xx_0_z_yyyzz[k] = -g_xx_0_0_yyyzz[k] * ab_z + g_xx_0_0_yyyzzz[k];

                g_xx_0_z_yyzzz[k] = -g_xx_0_0_yyzzz[k] * ab_z + g_xx_0_0_yyzzzz[k];

                g_xx_0_z_yzzzz[k] = -g_xx_0_0_yzzzz[k] * ab_z + g_xx_0_0_yzzzzz[k];

                g_xx_0_z_zzzzz[k] = -g_xx_0_0_zzzzz[k] * ab_z + g_xx_0_0_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xy_0_0_xxxxx, g_xy_0_0_xxxxxx, g_xy_0_0_xxxxxy, g_xy_0_0_xxxxxz, g_xy_0_0_xxxxy, g_xy_0_0_xxxxyy, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxz, g_xy_0_0_xxxxzz, g_xy_0_0_xxxyy, g_xy_0_0_xxxyyy, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxzz, g_xy_0_0_xxxzzz, g_xy_0_0_xxyyy, g_xy_0_0_xxyyyy, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxzzz, g_xy_0_0_xxzzzz, g_xy_0_0_xyyyy, g_xy_0_0_xyyyyy, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xzzzz, g_xy_0_0_xzzzzz, g_xy_0_0_yyyyy, g_xy_0_0_yyyyz, g_xy_0_0_yyyzz, g_xy_0_0_yyzzz, g_xy_0_0_yzzzz, g_xy_0_0_zzzzz, g_xy_0_x_xxxxx, g_xy_0_x_xxxxy, g_xy_0_x_xxxxz, g_xy_0_x_xxxyy, g_xy_0_x_xxxyz, g_xy_0_x_xxxzz, g_xy_0_x_xxyyy, g_xy_0_x_xxyyz, g_xy_0_x_xxyzz, g_xy_0_x_xxzzz, g_xy_0_x_xyyyy, g_xy_0_x_xyyyz, g_xy_0_x_xyyzz, g_xy_0_x_xyzzz, g_xy_0_x_xzzzz, g_xy_0_x_yyyyy, g_xy_0_x_yyyyz, g_xy_0_x_yyyzz, g_xy_0_x_yyzzz, g_xy_0_x_yzzzz, g_xy_0_x_zzzzz, g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_x_xxxxx[k] = -g_y_0_0_xxxxx[k] - g_xy_0_0_xxxxx[k] * ab_x + g_xy_0_0_xxxxxx[k];

                g_xy_0_x_xxxxy[k] = -g_y_0_0_xxxxy[k] - g_xy_0_0_xxxxy[k] * ab_x + g_xy_0_0_xxxxxy[k];

                g_xy_0_x_xxxxz[k] = -g_y_0_0_xxxxz[k] - g_xy_0_0_xxxxz[k] * ab_x + g_xy_0_0_xxxxxz[k];

                g_xy_0_x_xxxyy[k] = -g_y_0_0_xxxyy[k] - g_xy_0_0_xxxyy[k] * ab_x + g_xy_0_0_xxxxyy[k];

                g_xy_0_x_xxxyz[k] = -g_y_0_0_xxxyz[k] - g_xy_0_0_xxxyz[k] * ab_x + g_xy_0_0_xxxxyz[k];

                g_xy_0_x_xxxzz[k] = -g_y_0_0_xxxzz[k] - g_xy_0_0_xxxzz[k] * ab_x + g_xy_0_0_xxxxzz[k];

                g_xy_0_x_xxyyy[k] = -g_y_0_0_xxyyy[k] - g_xy_0_0_xxyyy[k] * ab_x + g_xy_0_0_xxxyyy[k];

                g_xy_0_x_xxyyz[k] = -g_y_0_0_xxyyz[k] - g_xy_0_0_xxyyz[k] * ab_x + g_xy_0_0_xxxyyz[k];

                g_xy_0_x_xxyzz[k] = -g_y_0_0_xxyzz[k] - g_xy_0_0_xxyzz[k] * ab_x + g_xy_0_0_xxxyzz[k];

                g_xy_0_x_xxzzz[k] = -g_y_0_0_xxzzz[k] - g_xy_0_0_xxzzz[k] * ab_x + g_xy_0_0_xxxzzz[k];

                g_xy_0_x_xyyyy[k] = -g_y_0_0_xyyyy[k] - g_xy_0_0_xyyyy[k] * ab_x + g_xy_0_0_xxyyyy[k];

                g_xy_0_x_xyyyz[k] = -g_y_0_0_xyyyz[k] - g_xy_0_0_xyyyz[k] * ab_x + g_xy_0_0_xxyyyz[k];

                g_xy_0_x_xyyzz[k] = -g_y_0_0_xyyzz[k] - g_xy_0_0_xyyzz[k] * ab_x + g_xy_0_0_xxyyzz[k];

                g_xy_0_x_xyzzz[k] = -g_y_0_0_xyzzz[k] - g_xy_0_0_xyzzz[k] * ab_x + g_xy_0_0_xxyzzz[k];

                g_xy_0_x_xzzzz[k] = -g_y_0_0_xzzzz[k] - g_xy_0_0_xzzzz[k] * ab_x + g_xy_0_0_xxzzzz[k];

                g_xy_0_x_yyyyy[k] = -g_y_0_0_yyyyy[k] - g_xy_0_0_yyyyy[k] * ab_x + g_xy_0_0_xyyyyy[k];

                g_xy_0_x_yyyyz[k] = -g_y_0_0_yyyyz[k] - g_xy_0_0_yyyyz[k] * ab_x + g_xy_0_0_xyyyyz[k];

                g_xy_0_x_yyyzz[k] = -g_y_0_0_yyyzz[k] - g_xy_0_0_yyyzz[k] * ab_x + g_xy_0_0_xyyyzz[k];

                g_xy_0_x_yyzzz[k] = -g_y_0_0_yyzzz[k] - g_xy_0_0_yyzzz[k] * ab_x + g_xy_0_0_xyyzzz[k];

                g_xy_0_x_yzzzz[k] = -g_y_0_0_yzzzz[k] - g_xy_0_0_yzzzz[k] * ab_x + g_xy_0_0_xyzzzz[k];

                g_xy_0_x_zzzzz[k] = -g_y_0_0_zzzzz[k] - g_xy_0_0_zzzzz[k] * ab_x + g_xy_0_0_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_xy_0_0_xxxxx, g_xy_0_0_xxxxxy, g_xy_0_0_xxxxy, g_xy_0_0_xxxxyy, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxz, g_xy_0_0_xxxyy, g_xy_0_0_xxxyyy, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxzz, g_xy_0_0_xxyyy, g_xy_0_0_xxyyyy, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxzzz, g_xy_0_0_xyyyy, g_xy_0_0_xyyyyy, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xzzzz, g_xy_0_0_yyyyy, g_xy_0_0_yyyyyy, g_xy_0_0_yyyyyz, g_xy_0_0_yyyyz, g_xy_0_0_yyyyzz, g_xy_0_0_yyyzz, g_xy_0_0_yyyzzz, g_xy_0_0_yyzzz, g_xy_0_0_yyzzzz, g_xy_0_0_yzzzz, g_xy_0_0_yzzzzz, g_xy_0_0_zzzzz, g_xy_0_y_xxxxx, g_xy_0_y_xxxxy, g_xy_0_y_xxxxz, g_xy_0_y_xxxyy, g_xy_0_y_xxxyz, g_xy_0_y_xxxzz, g_xy_0_y_xxyyy, g_xy_0_y_xxyyz, g_xy_0_y_xxyzz, g_xy_0_y_xxzzz, g_xy_0_y_xyyyy, g_xy_0_y_xyyyz, g_xy_0_y_xyyzz, g_xy_0_y_xyzzz, g_xy_0_y_xzzzz, g_xy_0_y_yyyyy, g_xy_0_y_yyyyz, g_xy_0_y_yyyzz, g_xy_0_y_yyzzz, g_xy_0_y_yzzzz, g_xy_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_y_xxxxx[k] = -g_x_0_0_xxxxx[k] - g_xy_0_0_xxxxx[k] * ab_y + g_xy_0_0_xxxxxy[k];

                g_xy_0_y_xxxxy[k] = -g_x_0_0_xxxxy[k] - g_xy_0_0_xxxxy[k] * ab_y + g_xy_0_0_xxxxyy[k];

                g_xy_0_y_xxxxz[k] = -g_x_0_0_xxxxz[k] - g_xy_0_0_xxxxz[k] * ab_y + g_xy_0_0_xxxxyz[k];

                g_xy_0_y_xxxyy[k] = -g_x_0_0_xxxyy[k] - g_xy_0_0_xxxyy[k] * ab_y + g_xy_0_0_xxxyyy[k];

                g_xy_0_y_xxxyz[k] = -g_x_0_0_xxxyz[k] - g_xy_0_0_xxxyz[k] * ab_y + g_xy_0_0_xxxyyz[k];

                g_xy_0_y_xxxzz[k] = -g_x_0_0_xxxzz[k] - g_xy_0_0_xxxzz[k] * ab_y + g_xy_0_0_xxxyzz[k];

                g_xy_0_y_xxyyy[k] = -g_x_0_0_xxyyy[k] - g_xy_0_0_xxyyy[k] * ab_y + g_xy_0_0_xxyyyy[k];

                g_xy_0_y_xxyyz[k] = -g_x_0_0_xxyyz[k] - g_xy_0_0_xxyyz[k] * ab_y + g_xy_0_0_xxyyyz[k];

                g_xy_0_y_xxyzz[k] = -g_x_0_0_xxyzz[k] - g_xy_0_0_xxyzz[k] * ab_y + g_xy_0_0_xxyyzz[k];

                g_xy_0_y_xxzzz[k] = -g_x_0_0_xxzzz[k] - g_xy_0_0_xxzzz[k] * ab_y + g_xy_0_0_xxyzzz[k];

                g_xy_0_y_xyyyy[k] = -g_x_0_0_xyyyy[k] - g_xy_0_0_xyyyy[k] * ab_y + g_xy_0_0_xyyyyy[k];

                g_xy_0_y_xyyyz[k] = -g_x_0_0_xyyyz[k] - g_xy_0_0_xyyyz[k] * ab_y + g_xy_0_0_xyyyyz[k];

                g_xy_0_y_xyyzz[k] = -g_x_0_0_xyyzz[k] - g_xy_0_0_xyyzz[k] * ab_y + g_xy_0_0_xyyyzz[k];

                g_xy_0_y_xyzzz[k] = -g_x_0_0_xyzzz[k] - g_xy_0_0_xyzzz[k] * ab_y + g_xy_0_0_xyyzzz[k];

                g_xy_0_y_xzzzz[k] = -g_x_0_0_xzzzz[k] - g_xy_0_0_xzzzz[k] * ab_y + g_xy_0_0_xyzzzz[k];

                g_xy_0_y_yyyyy[k] = -g_x_0_0_yyyyy[k] - g_xy_0_0_yyyyy[k] * ab_y + g_xy_0_0_yyyyyy[k];

                g_xy_0_y_yyyyz[k] = -g_x_0_0_yyyyz[k] - g_xy_0_0_yyyyz[k] * ab_y + g_xy_0_0_yyyyyz[k];

                g_xy_0_y_yyyzz[k] = -g_x_0_0_yyyzz[k] - g_xy_0_0_yyyzz[k] * ab_y + g_xy_0_0_yyyyzz[k];

                g_xy_0_y_yyzzz[k] = -g_x_0_0_yyzzz[k] - g_xy_0_0_yyzzz[k] * ab_y + g_xy_0_0_yyyzzz[k];

                g_xy_0_y_yzzzz[k] = -g_x_0_0_yzzzz[k] - g_xy_0_0_yzzzz[k] * ab_y + g_xy_0_0_yyzzzz[k];

                g_xy_0_y_zzzzz[k] = -g_x_0_0_zzzzz[k] - g_xy_0_0_zzzzz[k] * ab_y + g_xy_0_0_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xy_0_0_xxxxx, g_xy_0_0_xxxxxz, g_xy_0_0_xxxxy, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxz, g_xy_0_0_xxxxzz, g_xy_0_0_xxxyy, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxzz, g_xy_0_0_xxxzzz, g_xy_0_0_xxyyy, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxzzz, g_xy_0_0_xxzzzz, g_xy_0_0_xyyyy, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xzzzz, g_xy_0_0_xzzzzz, g_xy_0_0_yyyyy, g_xy_0_0_yyyyyz, g_xy_0_0_yyyyz, g_xy_0_0_yyyyzz, g_xy_0_0_yyyzz, g_xy_0_0_yyyzzz, g_xy_0_0_yyzzz, g_xy_0_0_yyzzzz, g_xy_0_0_yzzzz, g_xy_0_0_yzzzzz, g_xy_0_0_zzzzz, g_xy_0_0_zzzzzz, g_xy_0_z_xxxxx, g_xy_0_z_xxxxy, g_xy_0_z_xxxxz, g_xy_0_z_xxxyy, g_xy_0_z_xxxyz, g_xy_0_z_xxxzz, g_xy_0_z_xxyyy, g_xy_0_z_xxyyz, g_xy_0_z_xxyzz, g_xy_0_z_xxzzz, g_xy_0_z_xyyyy, g_xy_0_z_xyyyz, g_xy_0_z_xyyzz, g_xy_0_z_xyzzz, g_xy_0_z_xzzzz, g_xy_0_z_yyyyy, g_xy_0_z_yyyyz, g_xy_0_z_yyyzz, g_xy_0_z_yyzzz, g_xy_0_z_yzzzz, g_xy_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_z_xxxxx[k] = -g_xy_0_0_xxxxx[k] * ab_z + g_xy_0_0_xxxxxz[k];

                g_xy_0_z_xxxxy[k] = -g_xy_0_0_xxxxy[k] * ab_z + g_xy_0_0_xxxxyz[k];

                g_xy_0_z_xxxxz[k] = -g_xy_0_0_xxxxz[k] * ab_z + g_xy_0_0_xxxxzz[k];

                g_xy_0_z_xxxyy[k] = -g_xy_0_0_xxxyy[k] * ab_z + g_xy_0_0_xxxyyz[k];

                g_xy_0_z_xxxyz[k] = -g_xy_0_0_xxxyz[k] * ab_z + g_xy_0_0_xxxyzz[k];

                g_xy_0_z_xxxzz[k] = -g_xy_0_0_xxxzz[k] * ab_z + g_xy_0_0_xxxzzz[k];

                g_xy_0_z_xxyyy[k] = -g_xy_0_0_xxyyy[k] * ab_z + g_xy_0_0_xxyyyz[k];

                g_xy_0_z_xxyyz[k] = -g_xy_0_0_xxyyz[k] * ab_z + g_xy_0_0_xxyyzz[k];

                g_xy_0_z_xxyzz[k] = -g_xy_0_0_xxyzz[k] * ab_z + g_xy_0_0_xxyzzz[k];

                g_xy_0_z_xxzzz[k] = -g_xy_0_0_xxzzz[k] * ab_z + g_xy_0_0_xxzzzz[k];

                g_xy_0_z_xyyyy[k] = -g_xy_0_0_xyyyy[k] * ab_z + g_xy_0_0_xyyyyz[k];

                g_xy_0_z_xyyyz[k] = -g_xy_0_0_xyyyz[k] * ab_z + g_xy_0_0_xyyyzz[k];

                g_xy_0_z_xyyzz[k] = -g_xy_0_0_xyyzz[k] * ab_z + g_xy_0_0_xyyzzz[k];

                g_xy_0_z_xyzzz[k] = -g_xy_0_0_xyzzz[k] * ab_z + g_xy_0_0_xyzzzz[k];

                g_xy_0_z_xzzzz[k] = -g_xy_0_0_xzzzz[k] * ab_z + g_xy_0_0_xzzzzz[k];

                g_xy_0_z_yyyyy[k] = -g_xy_0_0_yyyyy[k] * ab_z + g_xy_0_0_yyyyyz[k];

                g_xy_0_z_yyyyz[k] = -g_xy_0_0_yyyyz[k] * ab_z + g_xy_0_0_yyyyzz[k];

                g_xy_0_z_yyyzz[k] = -g_xy_0_0_yyyzz[k] * ab_z + g_xy_0_0_yyyzzz[k];

                g_xy_0_z_yyzzz[k] = -g_xy_0_0_yyzzz[k] * ab_z + g_xy_0_0_yyzzzz[k];

                g_xy_0_z_yzzzz[k] = -g_xy_0_0_yzzzz[k] * ab_z + g_xy_0_0_yzzzzz[k];

                g_xy_0_z_zzzzz[k] = -g_xy_0_0_zzzzz[k] * ab_z + g_xy_0_0_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xz_0_0_xxxxx, g_xz_0_0_xxxxxx, g_xz_0_0_xxxxxy, g_xz_0_0_xxxxxz, g_xz_0_0_xxxxy, g_xz_0_0_xxxxyy, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxz, g_xz_0_0_xxxxzz, g_xz_0_0_xxxyy, g_xz_0_0_xxxyyy, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxzz, g_xz_0_0_xxxzzz, g_xz_0_0_xxyyy, g_xz_0_0_xxyyyy, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxzzz, g_xz_0_0_xxzzzz, g_xz_0_0_xyyyy, g_xz_0_0_xyyyyy, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xzzzz, g_xz_0_0_xzzzzz, g_xz_0_0_yyyyy, g_xz_0_0_yyyyz, g_xz_0_0_yyyzz, g_xz_0_0_yyzzz, g_xz_0_0_yzzzz, g_xz_0_0_zzzzz, g_xz_0_x_xxxxx, g_xz_0_x_xxxxy, g_xz_0_x_xxxxz, g_xz_0_x_xxxyy, g_xz_0_x_xxxyz, g_xz_0_x_xxxzz, g_xz_0_x_xxyyy, g_xz_0_x_xxyyz, g_xz_0_x_xxyzz, g_xz_0_x_xxzzz, g_xz_0_x_xyyyy, g_xz_0_x_xyyyz, g_xz_0_x_xyyzz, g_xz_0_x_xyzzz, g_xz_0_x_xzzzz, g_xz_0_x_yyyyy, g_xz_0_x_yyyyz, g_xz_0_x_yyyzz, g_xz_0_x_yyzzz, g_xz_0_x_yzzzz, g_xz_0_x_zzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_x_xxxxx[k] = -g_z_0_0_xxxxx[k] - g_xz_0_0_xxxxx[k] * ab_x + g_xz_0_0_xxxxxx[k];

                g_xz_0_x_xxxxy[k] = -g_z_0_0_xxxxy[k] - g_xz_0_0_xxxxy[k] * ab_x + g_xz_0_0_xxxxxy[k];

                g_xz_0_x_xxxxz[k] = -g_z_0_0_xxxxz[k] - g_xz_0_0_xxxxz[k] * ab_x + g_xz_0_0_xxxxxz[k];

                g_xz_0_x_xxxyy[k] = -g_z_0_0_xxxyy[k] - g_xz_0_0_xxxyy[k] * ab_x + g_xz_0_0_xxxxyy[k];

                g_xz_0_x_xxxyz[k] = -g_z_0_0_xxxyz[k] - g_xz_0_0_xxxyz[k] * ab_x + g_xz_0_0_xxxxyz[k];

                g_xz_0_x_xxxzz[k] = -g_z_0_0_xxxzz[k] - g_xz_0_0_xxxzz[k] * ab_x + g_xz_0_0_xxxxzz[k];

                g_xz_0_x_xxyyy[k] = -g_z_0_0_xxyyy[k] - g_xz_0_0_xxyyy[k] * ab_x + g_xz_0_0_xxxyyy[k];

                g_xz_0_x_xxyyz[k] = -g_z_0_0_xxyyz[k] - g_xz_0_0_xxyyz[k] * ab_x + g_xz_0_0_xxxyyz[k];

                g_xz_0_x_xxyzz[k] = -g_z_0_0_xxyzz[k] - g_xz_0_0_xxyzz[k] * ab_x + g_xz_0_0_xxxyzz[k];

                g_xz_0_x_xxzzz[k] = -g_z_0_0_xxzzz[k] - g_xz_0_0_xxzzz[k] * ab_x + g_xz_0_0_xxxzzz[k];

                g_xz_0_x_xyyyy[k] = -g_z_0_0_xyyyy[k] - g_xz_0_0_xyyyy[k] * ab_x + g_xz_0_0_xxyyyy[k];

                g_xz_0_x_xyyyz[k] = -g_z_0_0_xyyyz[k] - g_xz_0_0_xyyyz[k] * ab_x + g_xz_0_0_xxyyyz[k];

                g_xz_0_x_xyyzz[k] = -g_z_0_0_xyyzz[k] - g_xz_0_0_xyyzz[k] * ab_x + g_xz_0_0_xxyyzz[k];

                g_xz_0_x_xyzzz[k] = -g_z_0_0_xyzzz[k] - g_xz_0_0_xyzzz[k] * ab_x + g_xz_0_0_xxyzzz[k];

                g_xz_0_x_xzzzz[k] = -g_z_0_0_xzzzz[k] - g_xz_0_0_xzzzz[k] * ab_x + g_xz_0_0_xxzzzz[k];

                g_xz_0_x_yyyyy[k] = -g_z_0_0_yyyyy[k] - g_xz_0_0_yyyyy[k] * ab_x + g_xz_0_0_xyyyyy[k];

                g_xz_0_x_yyyyz[k] = -g_z_0_0_yyyyz[k] - g_xz_0_0_yyyyz[k] * ab_x + g_xz_0_0_xyyyyz[k];

                g_xz_0_x_yyyzz[k] = -g_z_0_0_yyyzz[k] - g_xz_0_0_yyyzz[k] * ab_x + g_xz_0_0_xyyyzz[k];

                g_xz_0_x_yyzzz[k] = -g_z_0_0_yyzzz[k] - g_xz_0_0_yyzzz[k] * ab_x + g_xz_0_0_xyyzzz[k];

                g_xz_0_x_yzzzz[k] = -g_z_0_0_yzzzz[k] - g_xz_0_0_yzzzz[k] * ab_x + g_xz_0_0_xyzzzz[k];

                g_xz_0_x_zzzzz[k] = -g_z_0_0_zzzzz[k] - g_xz_0_0_zzzzz[k] * ab_x + g_xz_0_0_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_xz_0_0_xxxxx, g_xz_0_0_xxxxxy, g_xz_0_0_xxxxy, g_xz_0_0_xxxxyy, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxz, g_xz_0_0_xxxyy, g_xz_0_0_xxxyyy, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxzz, g_xz_0_0_xxyyy, g_xz_0_0_xxyyyy, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxzzz, g_xz_0_0_xyyyy, g_xz_0_0_xyyyyy, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xzzzz, g_xz_0_0_yyyyy, g_xz_0_0_yyyyyy, g_xz_0_0_yyyyyz, g_xz_0_0_yyyyz, g_xz_0_0_yyyyzz, g_xz_0_0_yyyzz, g_xz_0_0_yyyzzz, g_xz_0_0_yyzzz, g_xz_0_0_yyzzzz, g_xz_0_0_yzzzz, g_xz_0_0_yzzzzz, g_xz_0_0_zzzzz, g_xz_0_y_xxxxx, g_xz_0_y_xxxxy, g_xz_0_y_xxxxz, g_xz_0_y_xxxyy, g_xz_0_y_xxxyz, g_xz_0_y_xxxzz, g_xz_0_y_xxyyy, g_xz_0_y_xxyyz, g_xz_0_y_xxyzz, g_xz_0_y_xxzzz, g_xz_0_y_xyyyy, g_xz_0_y_xyyyz, g_xz_0_y_xyyzz, g_xz_0_y_xyzzz, g_xz_0_y_xzzzz, g_xz_0_y_yyyyy, g_xz_0_y_yyyyz, g_xz_0_y_yyyzz, g_xz_0_y_yyzzz, g_xz_0_y_yzzzz, g_xz_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_y_xxxxx[k] = -g_xz_0_0_xxxxx[k] * ab_y + g_xz_0_0_xxxxxy[k];

                g_xz_0_y_xxxxy[k] = -g_xz_0_0_xxxxy[k] * ab_y + g_xz_0_0_xxxxyy[k];

                g_xz_0_y_xxxxz[k] = -g_xz_0_0_xxxxz[k] * ab_y + g_xz_0_0_xxxxyz[k];

                g_xz_0_y_xxxyy[k] = -g_xz_0_0_xxxyy[k] * ab_y + g_xz_0_0_xxxyyy[k];

                g_xz_0_y_xxxyz[k] = -g_xz_0_0_xxxyz[k] * ab_y + g_xz_0_0_xxxyyz[k];

                g_xz_0_y_xxxzz[k] = -g_xz_0_0_xxxzz[k] * ab_y + g_xz_0_0_xxxyzz[k];

                g_xz_0_y_xxyyy[k] = -g_xz_0_0_xxyyy[k] * ab_y + g_xz_0_0_xxyyyy[k];

                g_xz_0_y_xxyyz[k] = -g_xz_0_0_xxyyz[k] * ab_y + g_xz_0_0_xxyyyz[k];

                g_xz_0_y_xxyzz[k] = -g_xz_0_0_xxyzz[k] * ab_y + g_xz_0_0_xxyyzz[k];

                g_xz_0_y_xxzzz[k] = -g_xz_0_0_xxzzz[k] * ab_y + g_xz_0_0_xxyzzz[k];

                g_xz_0_y_xyyyy[k] = -g_xz_0_0_xyyyy[k] * ab_y + g_xz_0_0_xyyyyy[k];

                g_xz_0_y_xyyyz[k] = -g_xz_0_0_xyyyz[k] * ab_y + g_xz_0_0_xyyyyz[k];

                g_xz_0_y_xyyzz[k] = -g_xz_0_0_xyyzz[k] * ab_y + g_xz_0_0_xyyyzz[k];

                g_xz_0_y_xyzzz[k] = -g_xz_0_0_xyzzz[k] * ab_y + g_xz_0_0_xyyzzz[k];

                g_xz_0_y_xzzzz[k] = -g_xz_0_0_xzzzz[k] * ab_y + g_xz_0_0_xyzzzz[k];

                g_xz_0_y_yyyyy[k] = -g_xz_0_0_yyyyy[k] * ab_y + g_xz_0_0_yyyyyy[k];

                g_xz_0_y_yyyyz[k] = -g_xz_0_0_yyyyz[k] * ab_y + g_xz_0_0_yyyyyz[k];

                g_xz_0_y_yyyzz[k] = -g_xz_0_0_yyyzz[k] * ab_y + g_xz_0_0_yyyyzz[k];

                g_xz_0_y_yyzzz[k] = -g_xz_0_0_yyzzz[k] * ab_y + g_xz_0_0_yyyzzz[k];

                g_xz_0_y_yzzzz[k] = -g_xz_0_0_yzzzz[k] * ab_y + g_xz_0_0_yyzzzz[k];

                g_xz_0_y_zzzzz[k] = -g_xz_0_0_zzzzz[k] * ab_y + g_xz_0_0_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_xz_0_0_xxxxx, g_xz_0_0_xxxxxz, g_xz_0_0_xxxxy, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxz, g_xz_0_0_xxxxzz, g_xz_0_0_xxxyy, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxzz, g_xz_0_0_xxxzzz, g_xz_0_0_xxyyy, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxzzz, g_xz_0_0_xxzzzz, g_xz_0_0_xyyyy, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xzzzz, g_xz_0_0_xzzzzz, g_xz_0_0_yyyyy, g_xz_0_0_yyyyyz, g_xz_0_0_yyyyz, g_xz_0_0_yyyyzz, g_xz_0_0_yyyzz, g_xz_0_0_yyyzzz, g_xz_0_0_yyzzz, g_xz_0_0_yyzzzz, g_xz_0_0_yzzzz, g_xz_0_0_yzzzzz, g_xz_0_0_zzzzz, g_xz_0_0_zzzzzz, g_xz_0_z_xxxxx, g_xz_0_z_xxxxy, g_xz_0_z_xxxxz, g_xz_0_z_xxxyy, g_xz_0_z_xxxyz, g_xz_0_z_xxxzz, g_xz_0_z_xxyyy, g_xz_0_z_xxyyz, g_xz_0_z_xxyzz, g_xz_0_z_xxzzz, g_xz_0_z_xyyyy, g_xz_0_z_xyyyz, g_xz_0_z_xyyzz, g_xz_0_z_xyzzz, g_xz_0_z_xzzzz, g_xz_0_z_yyyyy, g_xz_0_z_yyyyz, g_xz_0_z_yyyzz, g_xz_0_z_yyzzz, g_xz_0_z_yzzzz, g_xz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_z_xxxxx[k] = -g_x_0_0_xxxxx[k] - g_xz_0_0_xxxxx[k] * ab_z + g_xz_0_0_xxxxxz[k];

                g_xz_0_z_xxxxy[k] = -g_x_0_0_xxxxy[k] - g_xz_0_0_xxxxy[k] * ab_z + g_xz_0_0_xxxxyz[k];

                g_xz_0_z_xxxxz[k] = -g_x_0_0_xxxxz[k] - g_xz_0_0_xxxxz[k] * ab_z + g_xz_0_0_xxxxzz[k];

                g_xz_0_z_xxxyy[k] = -g_x_0_0_xxxyy[k] - g_xz_0_0_xxxyy[k] * ab_z + g_xz_0_0_xxxyyz[k];

                g_xz_0_z_xxxyz[k] = -g_x_0_0_xxxyz[k] - g_xz_0_0_xxxyz[k] * ab_z + g_xz_0_0_xxxyzz[k];

                g_xz_0_z_xxxzz[k] = -g_x_0_0_xxxzz[k] - g_xz_0_0_xxxzz[k] * ab_z + g_xz_0_0_xxxzzz[k];

                g_xz_0_z_xxyyy[k] = -g_x_0_0_xxyyy[k] - g_xz_0_0_xxyyy[k] * ab_z + g_xz_0_0_xxyyyz[k];

                g_xz_0_z_xxyyz[k] = -g_x_0_0_xxyyz[k] - g_xz_0_0_xxyyz[k] * ab_z + g_xz_0_0_xxyyzz[k];

                g_xz_0_z_xxyzz[k] = -g_x_0_0_xxyzz[k] - g_xz_0_0_xxyzz[k] * ab_z + g_xz_0_0_xxyzzz[k];

                g_xz_0_z_xxzzz[k] = -g_x_0_0_xxzzz[k] - g_xz_0_0_xxzzz[k] * ab_z + g_xz_0_0_xxzzzz[k];

                g_xz_0_z_xyyyy[k] = -g_x_0_0_xyyyy[k] - g_xz_0_0_xyyyy[k] * ab_z + g_xz_0_0_xyyyyz[k];

                g_xz_0_z_xyyyz[k] = -g_x_0_0_xyyyz[k] - g_xz_0_0_xyyyz[k] * ab_z + g_xz_0_0_xyyyzz[k];

                g_xz_0_z_xyyzz[k] = -g_x_0_0_xyyzz[k] - g_xz_0_0_xyyzz[k] * ab_z + g_xz_0_0_xyyzzz[k];

                g_xz_0_z_xyzzz[k] = -g_x_0_0_xyzzz[k] - g_xz_0_0_xyzzz[k] * ab_z + g_xz_0_0_xyzzzz[k];

                g_xz_0_z_xzzzz[k] = -g_x_0_0_xzzzz[k] - g_xz_0_0_xzzzz[k] * ab_z + g_xz_0_0_xzzzzz[k];

                g_xz_0_z_yyyyy[k] = -g_x_0_0_yyyyy[k] - g_xz_0_0_yyyyy[k] * ab_z + g_xz_0_0_yyyyyz[k];

                g_xz_0_z_yyyyz[k] = -g_x_0_0_yyyyz[k] - g_xz_0_0_yyyyz[k] * ab_z + g_xz_0_0_yyyyzz[k];

                g_xz_0_z_yyyzz[k] = -g_x_0_0_yyyzz[k] - g_xz_0_0_yyyzz[k] * ab_z + g_xz_0_0_yyyzzz[k];

                g_xz_0_z_yyzzz[k] = -g_x_0_0_yyzzz[k] - g_xz_0_0_yyzzz[k] * ab_z + g_xz_0_0_yyzzzz[k];

                g_xz_0_z_yzzzz[k] = -g_x_0_0_yzzzz[k] - g_xz_0_0_yzzzz[k] * ab_z + g_xz_0_0_yzzzzz[k];

                g_xz_0_z_zzzzz[k] = -g_x_0_0_zzzzz[k] - g_xz_0_0_zzzzz[k] * ab_z + g_xz_0_0_zzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yy_0_0_xxxxx, g_yy_0_0_xxxxxx, g_yy_0_0_xxxxxy, g_yy_0_0_xxxxxz, g_yy_0_0_xxxxy, g_yy_0_0_xxxxyy, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxz, g_yy_0_0_xxxxzz, g_yy_0_0_xxxyy, g_yy_0_0_xxxyyy, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxzz, g_yy_0_0_xxxzzz, g_yy_0_0_xxyyy, g_yy_0_0_xxyyyy, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxzzz, g_yy_0_0_xxzzzz, g_yy_0_0_xyyyy, g_yy_0_0_xyyyyy, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xzzzz, g_yy_0_0_xzzzzz, g_yy_0_0_yyyyy, g_yy_0_0_yyyyz, g_yy_0_0_yyyzz, g_yy_0_0_yyzzz, g_yy_0_0_yzzzz, g_yy_0_0_zzzzz, g_yy_0_x_xxxxx, g_yy_0_x_xxxxy, g_yy_0_x_xxxxz, g_yy_0_x_xxxyy, g_yy_0_x_xxxyz, g_yy_0_x_xxxzz, g_yy_0_x_xxyyy, g_yy_0_x_xxyyz, g_yy_0_x_xxyzz, g_yy_0_x_xxzzz, g_yy_0_x_xyyyy, g_yy_0_x_xyyyz, g_yy_0_x_xyyzz, g_yy_0_x_xyzzz, g_yy_0_x_xzzzz, g_yy_0_x_yyyyy, g_yy_0_x_yyyyz, g_yy_0_x_yyyzz, g_yy_0_x_yyzzz, g_yy_0_x_yzzzz, g_yy_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_x_xxxxx[k] = -g_yy_0_0_xxxxx[k] * ab_x + g_yy_0_0_xxxxxx[k];

                g_yy_0_x_xxxxy[k] = -g_yy_0_0_xxxxy[k] * ab_x + g_yy_0_0_xxxxxy[k];

                g_yy_0_x_xxxxz[k] = -g_yy_0_0_xxxxz[k] * ab_x + g_yy_0_0_xxxxxz[k];

                g_yy_0_x_xxxyy[k] = -g_yy_0_0_xxxyy[k] * ab_x + g_yy_0_0_xxxxyy[k];

                g_yy_0_x_xxxyz[k] = -g_yy_0_0_xxxyz[k] * ab_x + g_yy_0_0_xxxxyz[k];

                g_yy_0_x_xxxzz[k] = -g_yy_0_0_xxxzz[k] * ab_x + g_yy_0_0_xxxxzz[k];

                g_yy_0_x_xxyyy[k] = -g_yy_0_0_xxyyy[k] * ab_x + g_yy_0_0_xxxyyy[k];

                g_yy_0_x_xxyyz[k] = -g_yy_0_0_xxyyz[k] * ab_x + g_yy_0_0_xxxyyz[k];

                g_yy_0_x_xxyzz[k] = -g_yy_0_0_xxyzz[k] * ab_x + g_yy_0_0_xxxyzz[k];

                g_yy_0_x_xxzzz[k] = -g_yy_0_0_xxzzz[k] * ab_x + g_yy_0_0_xxxzzz[k];

                g_yy_0_x_xyyyy[k] = -g_yy_0_0_xyyyy[k] * ab_x + g_yy_0_0_xxyyyy[k];

                g_yy_0_x_xyyyz[k] = -g_yy_0_0_xyyyz[k] * ab_x + g_yy_0_0_xxyyyz[k];

                g_yy_0_x_xyyzz[k] = -g_yy_0_0_xyyzz[k] * ab_x + g_yy_0_0_xxyyzz[k];

                g_yy_0_x_xyzzz[k] = -g_yy_0_0_xyzzz[k] * ab_x + g_yy_0_0_xxyzzz[k];

                g_yy_0_x_xzzzz[k] = -g_yy_0_0_xzzzz[k] * ab_x + g_yy_0_0_xxzzzz[k];

                g_yy_0_x_yyyyy[k] = -g_yy_0_0_yyyyy[k] * ab_x + g_yy_0_0_xyyyyy[k];

                g_yy_0_x_yyyyz[k] = -g_yy_0_0_yyyyz[k] * ab_x + g_yy_0_0_xyyyyz[k];

                g_yy_0_x_yyyzz[k] = -g_yy_0_0_yyyzz[k] * ab_x + g_yy_0_0_xyyyzz[k];

                g_yy_0_x_yyzzz[k] = -g_yy_0_0_yyzzz[k] * ab_x + g_yy_0_0_xyyzzz[k];

                g_yy_0_x_yzzzz[k] = -g_yy_0_0_yzzzz[k] * ab_x + g_yy_0_0_xyzzzz[k];

                g_yy_0_x_zzzzz[k] = -g_yy_0_0_zzzzz[k] * ab_x + g_yy_0_0_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_yy_0_0_xxxxx, g_yy_0_0_xxxxxy, g_yy_0_0_xxxxy, g_yy_0_0_xxxxyy, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxz, g_yy_0_0_xxxyy, g_yy_0_0_xxxyyy, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxzz, g_yy_0_0_xxyyy, g_yy_0_0_xxyyyy, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxzzz, g_yy_0_0_xyyyy, g_yy_0_0_xyyyyy, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xzzzz, g_yy_0_0_yyyyy, g_yy_0_0_yyyyyy, g_yy_0_0_yyyyyz, g_yy_0_0_yyyyz, g_yy_0_0_yyyyzz, g_yy_0_0_yyyzz, g_yy_0_0_yyyzzz, g_yy_0_0_yyzzz, g_yy_0_0_yyzzzz, g_yy_0_0_yzzzz, g_yy_0_0_yzzzzz, g_yy_0_0_zzzzz, g_yy_0_y_xxxxx, g_yy_0_y_xxxxy, g_yy_0_y_xxxxz, g_yy_0_y_xxxyy, g_yy_0_y_xxxyz, g_yy_0_y_xxxzz, g_yy_0_y_xxyyy, g_yy_0_y_xxyyz, g_yy_0_y_xxyzz, g_yy_0_y_xxzzz, g_yy_0_y_xyyyy, g_yy_0_y_xyyyz, g_yy_0_y_xyyzz, g_yy_0_y_xyzzz, g_yy_0_y_xzzzz, g_yy_0_y_yyyyy, g_yy_0_y_yyyyz, g_yy_0_y_yyyzz, g_yy_0_y_yyzzz, g_yy_0_y_yzzzz, g_yy_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_y_xxxxx[k] = -2.0 * g_y_0_0_xxxxx[k] - g_yy_0_0_xxxxx[k] * ab_y + g_yy_0_0_xxxxxy[k];

                g_yy_0_y_xxxxy[k] = -2.0 * g_y_0_0_xxxxy[k] - g_yy_0_0_xxxxy[k] * ab_y + g_yy_0_0_xxxxyy[k];

                g_yy_0_y_xxxxz[k] = -2.0 * g_y_0_0_xxxxz[k] - g_yy_0_0_xxxxz[k] * ab_y + g_yy_0_0_xxxxyz[k];

                g_yy_0_y_xxxyy[k] = -2.0 * g_y_0_0_xxxyy[k] - g_yy_0_0_xxxyy[k] * ab_y + g_yy_0_0_xxxyyy[k];

                g_yy_0_y_xxxyz[k] = -2.0 * g_y_0_0_xxxyz[k] - g_yy_0_0_xxxyz[k] * ab_y + g_yy_0_0_xxxyyz[k];

                g_yy_0_y_xxxzz[k] = -2.0 * g_y_0_0_xxxzz[k] - g_yy_0_0_xxxzz[k] * ab_y + g_yy_0_0_xxxyzz[k];

                g_yy_0_y_xxyyy[k] = -2.0 * g_y_0_0_xxyyy[k] - g_yy_0_0_xxyyy[k] * ab_y + g_yy_0_0_xxyyyy[k];

                g_yy_0_y_xxyyz[k] = -2.0 * g_y_0_0_xxyyz[k] - g_yy_0_0_xxyyz[k] * ab_y + g_yy_0_0_xxyyyz[k];

                g_yy_0_y_xxyzz[k] = -2.0 * g_y_0_0_xxyzz[k] - g_yy_0_0_xxyzz[k] * ab_y + g_yy_0_0_xxyyzz[k];

                g_yy_0_y_xxzzz[k] = -2.0 * g_y_0_0_xxzzz[k] - g_yy_0_0_xxzzz[k] * ab_y + g_yy_0_0_xxyzzz[k];

                g_yy_0_y_xyyyy[k] = -2.0 * g_y_0_0_xyyyy[k] - g_yy_0_0_xyyyy[k] * ab_y + g_yy_0_0_xyyyyy[k];

                g_yy_0_y_xyyyz[k] = -2.0 * g_y_0_0_xyyyz[k] - g_yy_0_0_xyyyz[k] * ab_y + g_yy_0_0_xyyyyz[k];

                g_yy_0_y_xyyzz[k] = -2.0 * g_y_0_0_xyyzz[k] - g_yy_0_0_xyyzz[k] * ab_y + g_yy_0_0_xyyyzz[k];

                g_yy_0_y_xyzzz[k] = -2.0 * g_y_0_0_xyzzz[k] - g_yy_0_0_xyzzz[k] * ab_y + g_yy_0_0_xyyzzz[k];

                g_yy_0_y_xzzzz[k] = -2.0 * g_y_0_0_xzzzz[k] - g_yy_0_0_xzzzz[k] * ab_y + g_yy_0_0_xyzzzz[k];

                g_yy_0_y_yyyyy[k] = -2.0 * g_y_0_0_yyyyy[k] - g_yy_0_0_yyyyy[k] * ab_y + g_yy_0_0_yyyyyy[k];

                g_yy_0_y_yyyyz[k] = -2.0 * g_y_0_0_yyyyz[k] - g_yy_0_0_yyyyz[k] * ab_y + g_yy_0_0_yyyyyz[k];

                g_yy_0_y_yyyzz[k] = -2.0 * g_y_0_0_yyyzz[k] - g_yy_0_0_yyyzz[k] * ab_y + g_yy_0_0_yyyyzz[k];

                g_yy_0_y_yyzzz[k] = -2.0 * g_y_0_0_yyzzz[k] - g_yy_0_0_yyzzz[k] * ab_y + g_yy_0_0_yyyzzz[k];

                g_yy_0_y_yzzzz[k] = -2.0 * g_y_0_0_yzzzz[k] - g_yy_0_0_yzzzz[k] * ab_y + g_yy_0_0_yyzzzz[k];

                g_yy_0_y_zzzzz[k] = -2.0 * g_y_0_0_zzzzz[k] - g_yy_0_0_zzzzz[k] * ab_y + g_yy_0_0_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yy_0_0_xxxxx, g_yy_0_0_xxxxxz, g_yy_0_0_xxxxy, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxz, g_yy_0_0_xxxxzz, g_yy_0_0_xxxyy, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxzz, g_yy_0_0_xxxzzz, g_yy_0_0_xxyyy, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxzzz, g_yy_0_0_xxzzzz, g_yy_0_0_xyyyy, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xzzzz, g_yy_0_0_xzzzzz, g_yy_0_0_yyyyy, g_yy_0_0_yyyyyz, g_yy_0_0_yyyyz, g_yy_0_0_yyyyzz, g_yy_0_0_yyyzz, g_yy_0_0_yyyzzz, g_yy_0_0_yyzzz, g_yy_0_0_yyzzzz, g_yy_0_0_yzzzz, g_yy_0_0_yzzzzz, g_yy_0_0_zzzzz, g_yy_0_0_zzzzzz, g_yy_0_z_xxxxx, g_yy_0_z_xxxxy, g_yy_0_z_xxxxz, g_yy_0_z_xxxyy, g_yy_0_z_xxxyz, g_yy_0_z_xxxzz, g_yy_0_z_xxyyy, g_yy_0_z_xxyyz, g_yy_0_z_xxyzz, g_yy_0_z_xxzzz, g_yy_0_z_xyyyy, g_yy_0_z_xyyyz, g_yy_0_z_xyyzz, g_yy_0_z_xyzzz, g_yy_0_z_xzzzz, g_yy_0_z_yyyyy, g_yy_0_z_yyyyz, g_yy_0_z_yyyzz, g_yy_0_z_yyzzz, g_yy_0_z_yzzzz, g_yy_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_z_xxxxx[k] = -g_yy_0_0_xxxxx[k] * ab_z + g_yy_0_0_xxxxxz[k];

                g_yy_0_z_xxxxy[k] = -g_yy_0_0_xxxxy[k] * ab_z + g_yy_0_0_xxxxyz[k];

                g_yy_0_z_xxxxz[k] = -g_yy_0_0_xxxxz[k] * ab_z + g_yy_0_0_xxxxzz[k];

                g_yy_0_z_xxxyy[k] = -g_yy_0_0_xxxyy[k] * ab_z + g_yy_0_0_xxxyyz[k];

                g_yy_0_z_xxxyz[k] = -g_yy_0_0_xxxyz[k] * ab_z + g_yy_0_0_xxxyzz[k];

                g_yy_0_z_xxxzz[k] = -g_yy_0_0_xxxzz[k] * ab_z + g_yy_0_0_xxxzzz[k];

                g_yy_0_z_xxyyy[k] = -g_yy_0_0_xxyyy[k] * ab_z + g_yy_0_0_xxyyyz[k];

                g_yy_0_z_xxyyz[k] = -g_yy_0_0_xxyyz[k] * ab_z + g_yy_0_0_xxyyzz[k];

                g_yy_0_z_xxyzz[k] = -g_yy_0_0_xxyzz[k] * ab_z + g_yy_0_0_xxyzzz[k];

                g_yy_0_z_xxzzz[k] = -g_yy_0_0_xxzzz[k] * ab_z + g_yy_0_0_xxzzzz[k];

                g_yy_0_z_xyyyy[k] = -g_yy_0_0_xyyyy[k] * ab_z + g_yy_0_0_xyyyyz[k];

                g_yy_0_z_xyyyz[k] = -g_yy_0_0_xyyyz[k] * ab_z + g_yy_0_0_xyyyzz[k];

                g_yy_0_z_xyyzz[k] = -g_yy_0_0_xyyzz[k] * ab_z + g_yy_0_0_xyyzzz[k];

                g_yy_0_z_xyzzz[k] = -g_yy_0_0_xyzzz[k] * ab_z + g_yy_0_0_xyzzzz[k];

                g_yy_0_z_xzzzz[k] = -g_yy_0_0_xzzzz[k] * ab_z + g_yy_0_0_xzzzzz[k];

                g_yy_0_z_yyyyy[k] = -g_yy_0_0_yyyyy[k] * ab_z + g_yy_0_0_yyyyyz[k];

                g_yy_0_z_yyyyz[k] = -g_yy_0_0_yyyyz[k] * ab_z + g_yy_0_0_yyyyzz[k];

                g_yy_0_z_yyyzz[k] = -g_yy_0_0_yyyzz[k] * ab_z + g_yy_0_0_yyyzzz[k];

                g_yy_0_z_yyzzz[k] = -g_yy_0_0_yyzzz[k] * ab_z + g_yy_0_0_yyzzzz[k];

                g_yy_0_z_yzzzz[k] = -g_yy_0_0_yzzzz[k] * ab_z + g_yy_0_0_yzzzzz[k];

                g_yy_0_z_zzzzz[k] = -g_yy_0_0_zzzzz[k] * ab_z + g_yy_0_0_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yz_0_0_xxxxx, g_yz_0_0_xxxxxx, g_yz_0_0_xxxxxy, g_yz_0_0_xxxxxz, g_yz_0_0_xxxxy, g_yz_0_0_xxxxyy, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxz, g_yz_0_0_xxxxzz, g_yz_0_0_xxxyy, g_yz_0_0_xxxyyy, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxzz, g_yz_0_0_xxxzzz, g_yz_0_0_xxyyy, g_yz_0_0_xxyyyy, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxzzz, g_yz_0_0_xxzzzz, g_yz_0_0_xyyyy, g_yz_0_0_xyyyyy, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xzzzz, g_yz_0_0_xzzzzz, g_yz_0_0_yyyyy, g_yz_0_0_yyyyz, g_yz_0_0_yyyzz, g_yz_0_0_yyzzz, g_yz_0_0_yzzzz, g_yz_0_0_zzzzz, g_yz_0_x_xxxxx, g_yz_0_x_xxxxy, g_yz_0_x_xxxxz, g_yz_0_x_xxxyy, g_yz_0_x_xxxyz, g_yz_0_x_xxxzz, g_yz_0_x_xxyyy, g_yz_0_x_xxyyz, g_yz_0_x_xxyzz, g_yz_0_x_xxzzz, g_yz_0_x_xyyyy, g_yz_0_x_xyyyz, g_yz_0_x_xyyzz, g_yz_0_x_xyzzz, g_yz_0_x_xzzzz, g_yz_0_x_yyyyy, g_yz_0_x_yyyyz, g_yz_0_x_yyyzz, g_yz_0_x_yyzzz, g_yz_0_x_yzzzz, g_yz_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_x_xxxxx[k] = -g_yz_0_0_xxxxx[k] * ab_x + g_yz_0_0_xxxxxx[k];

                g_yz_0_x_xxxxy[k] = -g_yz_0_0_xxxxy[k] * ab_x + g_yz_0_0_xxxxxy[k];

                g_yz_0_x_xxxxz[k] = -g_yz_0_0_xxxxz[k] * ab_x + g_yz_0_0_xxxxxz[k];

                g_yz_0_x_xxxyy[k] = -g_yz_0_0_xxxyy[k] * ab_x + g_yz_0_0_xxxxyy[k];

                g_yz_0_x_xxxyz[k] = -g_yz_0_0_xxxyz[k] * ab_x + g_yz_0_0_xxxxyz[k];

                g_yz_0_x_xxxzz[k] = -g_yz_0_0_xxxzz[k] * ab_x + g_yz_0_0_xxxxzz[k];

                g_yz_0_x_xxyyy[k] = -g_yz_0_0_xxyyy[k] * ab_x + g_yz_0_0_xxxyyy[k];

                g_yz_0_x_xxyyz[k] = -g_yz_0_0_xxyyz[k] * ab_x + g_yz_0_0_xxxyyz[k];

                g_yz_0_x_xxyzz[k] = -g_yz_0_0_xxyzz[k] * ab_x + g_yz_0_0_xxxyzz[k];

                g_yz_0_x_xxzzz[k] = -g_yz_0_0_xxzzz[k] * ab_x + g_yz_0_0_xxxzzz[k];

                g_yz_0_x_xyyyy[k] = -g_yz_0_0_xyyyy[k] * ab_x + g_yz_0_0_xxyyyy[k];

                g_yz_0_x_xyyyz[k] = -g_yz_0_0_xyyyz[k] * ab_x + g_yz_0_0_xxyyyz[k];

                g_yz_0_x_xyyzz[k] = -g_yz_0_0_xyyzz[k] * ab_x + g_yz_0_0_xxyyzz[k];

                g_yz_0_x_xyzzz[k] = -g_yz_0_0_xyzzz[k] * ab_x + g_yz_0_0_xxyzzz[k];

                g_yz_0_x_xzzzz[k] = -g_yz_0_0_xzzzz[k] * ab_x + g_yz_0_0_xxzzzz[k];

                g_yz_0_x_yyyyy[k] = -g_yz_0_0_yyyyy[k] * ab_x + g_yz_0_0_xyyyyy[k];

                g_yz_0_x_yyyyz[k] = -g_yz_0_0_yyyyz[k] * ab_x + g_yz_0_0_xyyyyz[k];

                g_yz_0_x_yyyzz[k] = -g_yz_0_0_yyyzz[k] * ab_x + g_yz_0_0_xyyyzz[k];

                g_yz_0_x_yyzzz[k] = -g_yz_0_0_yyzzz[k] * ab_x + g_yz_0_0_xyyzzz[k];

                g_yz_0_x_yzzzz[k] = -g_yz_0_0_yzzzz[k] * ab_x + g_yz_0_0_xyzzzz[k];

                g_yz_0_x_zzzzz[k] = -g_yz_0_0_zzzzz[k] * ab_x + g_yz_0_0_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_yz_0_0_xxxxx, g_yz_0_0_xxxxxy, g_yz_0_0_xxxxy, g_yz_0_0_xxxxyy, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxz, g_yz_0_0_xxxyy, g_yz_0_0_xxxyyy, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxzz, g_yz_0_0_xxyyy, g_yz_0_0_xxyyyy, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxzzz, g_yz_0_0_xyyyy, g_yz_0_0_xyyyyy, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xzzzz, g_yz_0_0_yyyyy, g_yz_0_0_yyyyyy, g_yz_0_0_yyyyyz, g_yz_0_0_yyyyz, g_yz_0_0_yyyyzz, g_yz_0_0_yyyzz, g_yz_0_0_yyyzzz, g_yz_0_0_yyzzz, g_yz_0_0_yyzzzz, g_yz_0_0_yzzzz, g_yz_0_0_yzzzzz, g_yz_0_0_zzzzz, g_yz_0_y_xxxxx, g_yz_0_y_xxxxy, g_yz_0_y_xxxxz, g_yz_0_y_xxxyy, g_yz_0_y_xxxyz, g_yz_0_y_xxxzz, g_yz_0_y_xxyyy, g_yz_0_y_xxyyz, g_yz_0_y_xxyzz, g_yz_0_y_xxzzz, g_yz_0_y_xyyyy, g_yz_0_y_xyyyz, g_yz_0_y_xyyzz, g_yz_0_y_xyzzz, g_yz_0_y_xzzzz, g_yz_0_y_yyyyy, g_yz_0_y_yyyyz, g_yz_0_y_yyyzz, g_yz_0_y_yyzzz, g_yz_0_y_yzzzz, g_yz_0_y_zzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_y_xxxxx[k] = -g_z_0_0_xxxxx[k] - g_yz_0_0_xxxxx[k] * ab_y + g_yz_0_0_xxxxxy[k];

                g_yz_0_y_xxxxy[k] = -g_z_0_0_xxxxy[k] - g_yz_0_0_xxxxy[k] * ab_y + g_yz_0_0_xxxxyy[k];

                g_yz_0_y_xxxxz[k] = -g_z_0_0_xxxxz[k] - g_yz_0_0_xxxxz[k] * ab_y + g_yz_0_0_xxxxyz[k];

                g_yz_0_y_xxxyy[k] = -g_z_0_0_xxxyy[k] - g_yz_0_0_xxxyy[k] * ab_y + g_yz_0_0_xxxyyy[k];

                g_yz_0_y_xxxyz[k] = -g_z_0_0_xxxyz[k] - g_yz_0_0_xxxyz[k] * ab_y + g_yz_0_0_xxxyyz[k];

                g_yz_0_y_xxxzz[k] = -g_z_0_0_xxxzz[k] - g_yz_0_0_xxxzz[k] * ab_y + g_yz_0_0_xxxyzz[k];

                g_yz_0_y_xxyyy[k] = -g_z_0_0_xxyyy[k] - g_yz_0_0_xxyyy[k] * ab_y + g_yz_0_0_xxyyyy[k];

                g_yz_0_y_xxyyz[k] = -g_z_0_0_xxyyz[k] - g_yz_0_0_xxyyz[k] * ab_y + g_yz_0_0_xxyyyz[k];

                g_yz_0_y_xxyzz[k] = -g_z_0_0_xxyzz[k] - g_yz_0_0_xxyzz[k] * ab_y + g_yz_0_0_xxyyzz[k];

                g_yz_0_y_xxzzz[k] = -g_z_0_0_xxzzz[k] - g_yz_0_0_xxzzz[k] * ab_y + g_yz_0_0_xxyzzz[k];

                g_yz_0_y_xyyyy[k] = -g_z_0_0_xyyyy[k] - g_yz_0_0_xyyyy[k] * ab_y + g_yz_0_0_xyyyyy[k];

                g_yz_0_y_xyyyz[k] = -g_z_0_0_xyyyz[k] - g_yz_0_0_xyyyz[k] * ab_y + g_yz_0_0_xyyyyz[k];

                g_yz_0_y_xyyzz[k] = -g_z_0_0_xyyzz[k] - g_yz_0_0_xyyzz[k] * ab_y + g_yz_0_0_xyyyzz[k];

                g_yz_0_y_xyzzz[k] = -g_z_0_0_xyzzz[k] - g_yz_0_0_xyzzz[k] * ab_y + g_yz_0_0_xyyzzz[k];

                g_yz_0_y_xzzzz[k] = -g_z_0_0_xzzzz[k] - g_yz_0_0_xzzzz[k] * ab_y + g_yz_0_0_xyzzzz[k];

                g_yz_0_y_yyyyy[k] = -g_z_0_0_yyyyy[k] - g_yz_0_0_yyyyy[k] * ab_y + g_yz_0_0_yyyyyy[k];

                g_yz_0_y_yyyyz[k] = -g_z_0_0_yyyyz[k] - g_yz_0_0_yyyyz[k] * ab_y + g_yz_0_0_yyyyyz[k];

                g_yz_0_y_yyyzz[k] = -g_z_0_0_yyyzz[k] - g_yz_0_0_yyyzz[k] * ab_y + g_yz_0_0_yyyyzz[k];

                g_yz_0_y_yyzzz[k] = -g_z_0_0_yyzzz[k] - g_yz_0_0_yyzzz[k] * ab_y + g_yz_0_0_yyyzzz[k];

                g_yz_0_y_yzzzz[k] = -g_z_0_0_yzzzz[k] - g_yz_0_0_yzzzz[k] * ab_y + g_yz_0_0_yyzzzz[k];

                g_yz_0_y_zzzzz[k] = -g_z_0_0_zzzzz[k] - g_yz_0_0_zzzzz[k] * ab_y + g_yz_0_0_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_yz_0_0_xxxxx, g_yz_0_0_xxxxxz, g_yz_0_0_xxxxy, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxz, g_yz_0_0_xxxxzz, g_yz_0_0_xxxyy, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxzz, g_yz_0_0_xxxzzz, g_yz_0_0_xxyyy, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxzzz, g_yz_0_0_xxzzzz, g_yz_0_0_xyyyy, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xzzzz, g_yz_0_0_xzzzzz, g_yz_0_0_yyyyy, g_yz_0_0_yyyyyz, g_yz_0_0_yyyyz, g_yz_0_0_yyyyzz, g_yz_0_0_yyyzz, g_yz_0_0_yyyzzz, g_yz_0_0_yyzzz, g_yz_0_0_yyzzzz, g_yz_0_0_yzzzz, g_yz_0_0_yzzzzz, g_yz_0_0_zzzzz, g_yz_0_0_zzzzzz, g_yz_0_z_xxxxx, g_yz_0_z_xxxxy, g_yz_0_z_xxxxz, g_yz_0_z_xxxyy, g_yz_0_z_xxxyz, g_yz_0_z_xxxzz, g_yz_0_z_xxyyy, g_yz_0_z_xxyyz, g_yz_0_z_xxyzz, g_yz_0_z_xxzzz, g_yz_0_z_xyyyy, g_yz_0_z_xyyyz, g_yz_0_z_xyyzz, g_yz_0_z_xyzzz, g_yz_0_z_xzzzz, g_yz_0_z_yyyyy, g_yz_0_z_yyyyz, g_yz_0_z_yyyzz, g_yz_0_z_yyzzz, g_yz_0_z_yzzzz, g_yz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_z_xxxxx[k] = -g_y_0_0_xxxxx[k] - g_yz_0_0_xxxxx[k] * ab_z + g_yz_0_0_xxxxxz[k];

                g_yz_0_z_xxxxy[k] = -g_y_0_0_xxxxy[k] - g_yz_0_0_xxxxy[k] * ab_z + g_yz_0_0_xxxxyz[k];

                g_yz_0_z_xxxxz[k] = -g_y_0_0_xxxxz[k] - g_yz_0_0_xxxxz[k] * ab_z + g_yz_0_0_xxxxzz[k];

                g_yz_0_z_xxxyy[k] = -g_y_0_0_xxxyy[k] - g_yz_0_0_xxxyy[k] * ab_z + g_yz_0_0_xxxyyz[k];

                g_yz_0_z_xxxyz[k] = -g_y_0_0_xxxyz[k] - g_yz_0_0_xxxyz[k] * ab_z + g_yz_0_0_xxxyzz[k];

                g_yz_0_z_xxxzz[k] = -g_y_0_0_xxxzz[k] - g_yz_0_0_xxxzz[k] * ab_z + g_yz_0_0_xxxzzz[k];

                g_yz_0_z_xxyyy[k] = -g_y_0_0_xxyyy[k] - g_yz_0_0_xxyyy[k] * ab_z + g_yz_0_0_xxyyyz[k];

                g_yz_0_z_xxyyz[k] = -g_y_0_0_xxyyz[k] - g_yz_0_0_xxyyz[k] * ab_z + g_yz_0_0_xxyyzz[k];

                g_yz_0_z_xxyzz[k] = -g_y_0_0_xxyzz[k] - g_yz_0_0_xxyzz[k] * ab_z + g_yz_0_0_xxyzzz[k];

                g_yz_0_z_xxzzz[k] = -g_y_0_0_xxzzz[k] - g_yz_0_0_xxzzz[k] * ab_z + g_yz_0_0_xxzzzz[k];

                g_yz_0_z_xyyyy[k] = -g_y_0_0_xyyyy[k] - g_yz_0_0_xyyyy[k] * ab_z + g_yz_0_0_xyyyyz[k];

                g_yz_0_z_xyyyz[k] = -g_y_0_0_xyyyz[k] - g_yz_0_0_xyyyz[k] * ab_z + g_yz_0_0_xyyyzz[k];

                g_yz_0_z_xyyzz[k] = -g_y_0_0_xyyzz[k] - g_yz_0_0_xyyzz[k] * ab_z + g_yz_0_0_xyyzzz[k];

                g_yz_0_z_xyzzz[k] = -g_y_0_0_xyzzz[k] - g_yz_0_0_xyzzz[k] * ab_z + g_yz_0_0_xyzzzz[k];

                g_yz_0_z_xzzzz[k] = -g_y_0_0_xzzzz[k] - g_yz_0_0_xzzzz[k] * ab_z + g_yz_0_0_xzzzzz[k];

                g_yz_0_z_yyyyy[k] = -g_y_0_0_yyyyy[k] - g_yz_0_0_yyyyy[k] * ab_z + g_yz_0_0_yyyyyz[k];

                g_yz_0_z_yyyyz[k] = -g_y_0_0_yyyyz[k] - g_yz_0_0_yyyyz[k] * ab_z + g_yz_0_0_yyyyzz[k];

                g_yz_0_z_yyyzz[k] = -g_y_0_0_yyyzz[k] - g_yz_0_0_yyyzz[k] * ab_z + g_yz_0_0_yyyzzz[k];

                g_yz_0_z_yyzzz[k] = -g_y_0_0_yyzzz[k] - g_yz_0_0_yyzzz[k] * ab_z + g_yz_0_0_yyzzzz[k];

                g_yz_0_z_yzzzz[k] = -g_y_0_0_yzzzz[k] - g_yz_0_0_yzzzz[k] * ab_z + g_yz_0_0_yzzzzz[k];

                g_yz_0_z_zzzzz[k] = -g_y_0_0_zzzzz[k] - g_yz_0_0_zzzzz[k] * ab_z + g_yz_0_0_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_zz_0_0_xxxxx, g_zz_0_0_xxxxxx, g_zz_0_0_xxxxxy, g_zz_0_0_xxxxxz, g_zz_0_0_xxxxy, g_zz_0_0_xxxxyy, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxz, g_zz_0_0_xxxxzz, g_zz_0_0_xxxyy, g_zz_0_0_xxxyyy, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxzz, g_zz_0_0_xxxzzz, g_zz_0_0_xxyyy, g_zz_0_0_xxyyyy, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxzzz, g_zz_0_0_xxzzzz, g_zz_0_0_xyyyy, g_zz_0_0_xyyyyy, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xzzzz, g_zz_0_0_xzzzzz, g_zz_0_0_yyyyy, g_zz_0_0_yyyyz, g_zz_0_0_yyyzz, g_zz_0_0_yyzzz, g_zz_0_0_yzzzz, g_zz_0_0_zzzzz, g_zz_0_x_xxxxx, g_zz_0_x_xxxxy, g_zz_0_x_xxxxz, g_zz_0_x_xxxyy, g_zz_0_x_xxxyz, g_zz_0_x_xxxzz, g_zz_0_x_xxyyy, g_zz_0_x_xxyyz, g_zz_0_x_xxyzz, g_zz_0_x_xxzzz, g_zz_0_x_xyyyy, g_zz_0_x_xyyyz, g_zz_0_x_xyyzz, g_zz_0_x_xyzzz, g_zz_0_x_xzzzz, g_zz_0_x_yyyyy, g_zz_0_x_yyyyz, g_zz_0_x_yyyzz, g_zz_0_x_yyzzz, g_zz_0_x_yzzzz, g_zz_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_x_xxxxx[k] = -g_zz_0_0_xxxxx[k] * ab_x + g_zz_0_0_xxxxxx[k];

                g_zz_0_x_xxxxy[k] = -g_zz_0_0_xxxxy[k] * ab_x + g_zz_0_0_xxxxxy[k];

                g_zz_0_x_xxxxz[k] = -g_zz_0_0_xxxxz[k] * ab_x + g_zz_0_0_xxxxxz[k];

                g_zz_0_x_xxxyy[k] = -g_zz_0_0_xxxyy[k] * ab_x + g_zz_0_0_xxxxyy[k];

                g_zz_0_x_xxxyz[k] = -g_zz_0_0_xxxyz[k] * ab_x + g_zz_0_0_xxxxyz[k];

                g_zz_0_x_xxxzz[k] = -g_zz_0_0_xxxzz[k] * ab_x + g_zz_0_0_xxxxzz[k];

                g_zz_0_x_xxyyy[k] = -g_zz_0_0_xxyyy[k] * ab_x + g_zz_0_0_xxxyyy[k];

                g_zz_0_x_xxyyz[k] = -g_zz_0_0_xxyyz[k] * ab_x + g_zz_0_0_xxxyyz[k];

                g_zz_0_x_xxyzz[k] = -g_zz_0_0_xxyzz[k] * ab_x + g_zz_0_0_xxxyzz[k];

                g_zz_0_x_xxzzz[k] = -g_zz_0_0_xxzzz[k] * ab_x + g_zz_0_0_xxxzzz[k];

                g_zz_0_x_xyyyy[k] = -g_zz_0_0_xyyyy[k] * ab_x + g_zz_0_0_xxyyyy[k];

                g_zz_0_x_xyyyz[k] = -g_zz_0_0_xyyyz[k] * ab_x + g_zz_0_0_xxyyyz[k];

                g_zz_0_x_xyyzz[k] = -g_zz_0_0_xyyzz[k] * ab_x + g_zz_0_0_xxyyzz[k];

                g_zz_0_x_xyzzz[k] = -g_zz_0_0_xyzzz[k] * ab_x + g_zz_0_0_xxyzzz[k];

                g_zz_0_x_xzzzz[k] = -g_zz_0_0_xzzzz[k] * ab_x + g_zz_0_0_xxzzzz[k];

                g_zz_0_x_yyyyy[k] = -g_zz_0_0_yyyyy[k] * ab_x + g_zz_0_0_xyyyyy[k];

                g_zz_0_x_yyyyz[k] = -g_zz_0_0_yyyyz[k] * ab_x + g_zz_0_0_xyyyyz[k];

                g_zz_0_x_yyyzz[k] = -g_zz_0_0_yyyzz[k] * ab_x + g_zz_0_0_xyyyzz[k];

                g_zz_0_x_yyzzz[k] = -g_zz_0_0_yyzzz[k] * ab_x + g_zz_0_0_xyyzzz[k];

                g_zz_0_x_yzzzz[k] = -g_zz_0_0_yzzzz[k] * ab_x + g_zz_0_0_xyzzzz[k];

                g_zz_0_x_zzzzz[k] = -g_zz_0_0_zzzzz[k] * ab_x + g_zz_0_0_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_zz_0_0_xxxxx, g_zz_0_0_xxxxxy, g_zz_0_0_xxxxy, g_zz_0_0_xxxxyy, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxz, g_zz_0_0_xxxyy, g_zz_0_0_xxxyyy, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxzz, g_zz_0_0_xxyyy, g_zz_0_0_xxyyyy, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxzzz, g_zz_0_0_xyyyy, g_zz_0_0_xyyyyy, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xzzzz, g_zz_0_0_yyyyy, g_zz_0_0_yyyyyy, g_zz_0_0_yyyyyz, g_zz_0_0_yyyyz, g_zz_0_0_yyyyzz, g_zz_0_0_yyyzz, g_zz_0_0_yyyzzz, g_zz_0_0_yyzzz, g_zz_0_0_yyzzzz, g_zz_0_0_yzzzz, g_zz_0_0_yzzzzz, g_zz_0_0_zzzzz, g_zz_0_y_xxxxx, g_zz_0_y_xxxxy, g_zz_0_y_xxxxz, g_zz_0_y_xxxyy, g_zz_0_y_xxxyz, g_zz_0_y_xxxzz, g_zz_0_y_xxyyy, g_zz_0_y_xxyyz, g_zz_0_y_xxyzz, g_zz_0_y_xxzzz, g_zz_0_y_xyyyy, g_zz_0_y_xyyyz, g_zz_0_y_xyyzz, g_zz_0_y_xyzzz, g_zz_0_y_xzzzz, g_zz_0_y_yyyyy, g_zz_0_y_yyyyz, g_zz_0_y_yyyzz, g_zz_0_y_yyzzz, g_zz_0_y_yzzzz, g_zz_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_y_xxxxx[k] = -g_zz_0_0_xxxxx[k] * ab_y + g_zz_0_0_xxxxxy[k];

                g_zz_0_y_xxxxy[k] = -g_zz_0_0_xxxxy[k] * ab_y + g_zz_0_0_xxxxyy[k];

                g_zz_0_y_xxxxz[k] = -g_zz_0_0_xxxxz[k] * ab_y + g_zz_0_0_xxxxyz[k];

                g_zz_0_y_xxxyy[k] = -g_zz_0_0_xxxyy[k] * ab_y + g_zz_0_0_xxxyyy[k];

                g_zz_0_y_xxxyz[k] = -g_zz_0_0_xxxyz[k] * ab_y + g_zz_0_0_xxxyyz[k];

                g_zz_0_y_xxxzz[k] = -g_zz_0_0_xxxzz[k] * ab_y + g_zz_0_0_xxxyzz[k];

                g_zz_0_y_xxyyy[k] = -g_zz_0_0_xxyyy[k] * ab_y + g_zz_0_0_xxyyyy[k];

                g_zz_0_y_xxyyz[k] = -g_zz_0_0_xxyyz[k] * ab_y + g_zz_0_0_xxyyyz[k];

                g_zz_0_y_xxyzz[k] = -g_zz_0_0_xxyzz[k] * ab_y + g_zz_0_0_xxyyzz[k];

                g_zz_0_y_xxzzz[k] = -g_zz_0_0_xxzzz[k] * ab_y + g_zz_0_0_xxyzzz[k];

                g_zz_0_y_xyyyy[k] = -g_zz_0_0_xyyyy[k] * ab_y + g_zz_0_0_xyyyyy[k];

                g_zz_0_y_xyyyz[k] = -g_zz_0_0_xyyyz[k] * ab_y + g_zz_0_0_xyyyyz[k];

                g_zz_0_y_xyyzz[k] = -g_zz_0_0_xyyzz[k] * ab_y + g_zz_0_0_xyyyzz[k];

                g_zz_0_y_xyzzz[k] = -g_zz_0_0_xyzzz[k] * ab_y + g_zz_0_0_xyyzzz[k];

                g_zz_0_y_xzzzz[k] = -g_zz_0_0_xzzzz[k] * ab_y + g_zz_0_0_xyzzzz[k];

                g_zz_0_y_yyyyy[k] = -g_zz_0_0_yyyyy[k] * ab_y + g_zz_0_0_yyyyyy[k];

                g_zz_0_y_yyyyz[k] = -g_zz_0_0_yyyyz[k] * ab_y + g_zz_0_0_yyyyyz[k];

                g_zz_0_y_yyyzz[k] = -g_zz_0_0_yyyzz[k] * ab_y + g_zz_0_0_yyyyzz[k];

                g_zz_0_y_yyzzz[k] = -g_zz_0_0_yyzzz[k] * ab_y + g_zz_0_0_yyyzzz[k];

                g_zz_0_y_yzzzz[k] = -g_zz_0_0_yzzzz[k] * ab_y + g_zz_0_0_yyzzzz[k];

                g_zz_0_y_zzzzz[k] = -g_zz_0_0_zzzzz[k] * ab_y + g_zz_0_0_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_zz_0_0_xxxxx, g_zz_0_0_xxxxxz, g_zz_0_0_xxxxy, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxz, g_zz_0_0_xxxxzz, g_zz_0_0_xxxyy, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxzz, g_zz_0_0_xxxzzz, g_zz_0_0_xxyyy, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxzzz, g_zz_0_0_xxzzzz, g_zz_0_0_xyyyy, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xzzzz, g_zz_0_0_xzzzzz, g_zz_0_0_yyyyy, g_zz_0_0_yyyyyz, g_zz_0_0_yyyyz, g_zz_0_0_yyyyzz, g_zz_0_0_yyyzz, g_zz_0_0_yyyzzz, g_zz_0_0_yyzzz, g_zz_0_0_yyzzzz, g_zz_0_0_yzzzz, g_zz_0_0_yzzzzz, g_zz_0_0_zzzzz, g_zz_0_0_zzzzzz, g_zz_0_z_xxxxx, g_zz_0_z_xxxxy, g_zz_0_z_xxxxz, g_zz_0_z_xxxyy, g_zz_0_z_xxxyz, g_zz_0_z_xxxzz, g_zz_0_z_xxyyy, g_zz_0_z_xxyyz, g_zz_0_z_xxyzz, g_zz_0_z_xxzzz, g_zz_0_z_xyyyy, g_zz_0_z_xyyyz, g_zz_0_z_xyyzz, g_zz_0_z_xyzzz, g_zz_0_z_xzzzz, g_zz_0_z_yyyyy, g_zz_0_z_yyyyz, g_zz_0_z_yyyzz, g_zz_0_z_yyzzz, g_zz_0_z_yzzzz, g_zz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_z_xxxxx[k] = -2.0 * g_z_0_0_xxxxx[k] - g_zz_0_0_xxxxx[k] * ab_z + g_zz_0_0_xxxxxz[k];

                g_zz_0_z_xxxxy[k] = -2.0 * g_z_0_0_xxxxy[k] - g_zz_0_0_xxxxy[k] * ab_z + g_zz_0_0_xxxxyz[k];

                g_zz_0_z_xxxxz[k] = -2.0 * g_z_0_0_xxxxz[k] - g_zz_0_0_xxxxz[k] * ab_z + g_zz_0_0_xxxxzz[k];

                g_zz_0_z_xxxyy[k] = -2.0 * g_z_0_0_xxxyy[k] - g_zz_0_0_xxxyy[k] * ab_z + g_zz_0_0_xxxyyz[k];

                g_zz_0_z_xxxyz[k] = -2.0 * g_z_0_0_xxxyz[k] - g_zz_0_0_xxxyz[k] * ab_z + g_zz_0_0_xxxyzz[k];

                g_zz_0_z_xxxzz[k] = -2.0 * g_z_0_0_xxxzz[k] - g_zz_0_0_xxxzz[k] * ab_z + g_zz_0_0_xxxzzz[k];

                g_zz_0_z_xxyyy[k] = -2.0 * g_z_0_0_xxyyy[k] - g_zz_0_0_xxyyy[k] * ab_z + g_zz_0_0_xxyyyz[k];

                g_zz_0_z_xxyyz[k] = -2.0 * g_z_0_0_xxyyz[k] - g_zz_0_0_xxyyz[k] * ab_z + g_zz_0_0_xxyyzz[k];

                g_zz_0_z_xxyzz[k] = -2.0 * g_z_0_0_xxyzz[k] - g_zz_0_0_xxyzz[k] * ab_z + g_zz_0_0_xxyzzz[k];

                g_zz_0_z_xxzzz[k] = -2.0 * g_z_0_0_xxzzz[k] - g_zz_0_0_xxzzz[k] * ab_z + g_zz_0_0_xxzzzz[k];

                g_zz_0_z_xyyyy[k] = -2.0 * g_z_0_0_xyyyy[k] - g_zz_0_0_xyyyy[k] * ab_z + g_zz_0_0_xyyyyz[k];

                g_zz_0_z_xyyyz[k] = -2.0 * g_z_0_0_xyyyz[k] - g_zz_0_0_xyyyz[k] * ab_z + g_zz_0_0_xyyyzz[k];

                g_zz_0_z_xyyzz[k] = -2.0 * g_z_0_0_xyyzz[k] - g_zz_0_0_xyyzz[k] * ab_z + g_zz_0_0_xyyzzz[k];

                g_zz_0_z_xyzzz[k] = -2.0 * g_z_0_0_xyzzz[k] - g_zz_0_0_xyzzz[k] * ab_z + g_zz_0_0_xyzzzz[k];

                g_zz_0_z_xzzzz[k] = -2.0 * g_z_0_0_xzzzz[k] - g_zz_0_0_xzzzz[k] * ab_z + g_zz_0_0_xzzzzz[k];

                g_zz_0_z_yyyyy[k] = -2.0 * g_z_0_0_yyyyy[k] - g_zz_0_0_yyyyy[k] * ab_z + g_zz_0_0_yyyyyz[k];

                g_zz_0_z_yyyyz[k] = -2.0 * g_z_0_0_yyyyz[k] - g_zz_0_0_yyyyz[k] * ab_z + g_zz_0_0_yyyyzz[k];

                g_zz_0_z_yyyzz[k] = -2.0 * g_z_0_0_yyyzz[k] - g_zz_0_0_yyyzz[k] * ab_z + g_zz_0_0_yyyzzz[k];

                g_zz_0_z_yyzzz[k] = -2.0 * g_z_0_0_yyzzz[k] - g_zz_0_0_yyzzz[k] * ab_z + g_zz_0_0_yyzzzz[k];

                g_zz_0_z_yzzzz[k] = -2.0 * g_z_0_0_yzzzz[k] - g_zz_0_0_yzzzz[k] * ab_z + g_zz_0_0_yzzzzz[k];

                g_zz_0_z_zzzzz[k] = -2.0 * g_z_0_0_zzzzz[k] - g_zz_0_0_zzzzz[k] * ab_z + g_zz_0_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

