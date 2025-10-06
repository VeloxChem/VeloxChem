#include "ElectronRepulsionGeom2000ContrRecISXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_isxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_isxx,
                                            const size_t idx_geom_10_hsxx,
                                            const size_t idx_geom_20_hsxx,
                                            const size_t idx_geom_20_hpxx,
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
            /// Set up components of auxilary buffer : HSSS

            const auto hs_geom_10_off = idx_geom_10_hsxx + i * dcomps + j;

            auto g_x_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 53 * ccomps * dcomps);

            auto g_z_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 54 * ccomps * dcomps);

            auto g_z_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 62 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HSSS

            const auto hs_geom_20_off = idx_geom_20_hsxx + i * dcomps + j;

            auto g_xx_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 44 * ccomps * dcomps);

            auto g_xz_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 45 * ccomps * dcomps);

            auto g_xz_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 46 * ccomps * dcomps);

            auto g_xz_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 47 * ccomps * dcomps);

            auto g_xz_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 48 * ccomps * dcomps);

            auto g_xz_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 49 * ccomps * dcomps);

            auto g_xz_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 50 * ccomps * dcomps);

            auto g_xz_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 51 * ccomps * dcomps);

            auto g_xz_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 52 * ccomps * dcomps);

            auto g_xz_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 53 * ccomps * dcomps);

            auto g_xz_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 54 * ccomps * dcomps);

            auto g_xz_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 55 * ccomps * dcomps);

            auto g_xz_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 56 * ccomps * dcomps);

            auto g_xz_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 57 * ccomps * dcomps);

            auto g_xz_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 58 * ccomps * dcomps);

            auto g_xz_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 62 * ccomps * dcomps);

            auto g_yy_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 63 * ccomps * dcomps);

            auto g_yy_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 64 * ccomps * dcomps);

            auto g_yy_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 65 * ccomps * dcomps);

            auto g_yy_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 66 * ccomps * dcomps);

            auto g_yy_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 67 * ccomps * dcomps);

            auto g_yy_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 68 * ccomps * dcomps);

            auto g_yy_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 69 * ccomps * dcomps);

            auto g_yy_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 70 * ccomps * dcomps);

            auto g_yy_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 71 * ccomps * dcomps);

            auto g_yy_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 72 * ccomps * dcomps);

            auto g_yy_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 73 * ccomps * dcomps);

            auto g_yy_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 74 * ccomps * dcomps);

            auto g_yy_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 75 * ccomps * dcomps);

            auto g_yy_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 76 * ccomps * dcomps);

            auto g_yy_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 77 * ccomps * dcomps);

            auto g_yy_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 78 * ccomps * dcomps);

            auto g_yy_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 79 * ccomps * dcomps);

            auto g_yy_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 80 * ccomps * dcomps);

            auto g_yy_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 81 * ccomps * dcomps);

            auto g_yy_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 82 * ccomps * dcomps);

            auto g_yy_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 83 * ccomps * dcomps);

            auto g_yz_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 84 * ccomps * dcomps);

            auto g_yz_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 85 * ccomps * dcomps);

            auto g_yz_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 86 * ccomps * dcomps);

            auto g_yz_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 87 * ccomps * dcomps);

            auto g_yz_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 88 * ccomps * dcomps);

            auto g_yz_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 89 * ccomps * dcomps);

            auto g_yz_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 90 * ccomps * dcomps);

            auto g_yz_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 91 * ccomps * dcomps);

            auto g_yz_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 92 * ccomps * dcomps);

            auto g_yz_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 93 * ccomps * dcomps);

            auto g_yz_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 94 * ccomps * dcomps);

            auto g_yz_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 95 * ccomps * dcomps);

            auto g_yz_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 96 * ccomps * dcomps);

            auto g_yz_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 97 * ccomps * dcomps);

            auto g_yz_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 98 * ccomps * dcomps);

            auto g_yz_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 99 * ccomps * dcomps);

            auto g_yz_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 100 * ccomps * dcomps);

            auto g_yz_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 101 * ccomps * dcomps);

            auto g_yz_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 102 * ccomps * dcomps);

            auto g_yz_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 103 * ccomps * dcomps);

            auto g_yz_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 104 * ccomps * dcomps);

            auto g_zz_0_xxxxx_0 = cbuffer.data(hs_geom_20_off + 105 * ccomps * dcomps);

            auto g_zz_0_xxxxy_0 = cbuffer.data(hs_geom_20_off + 106 * ccomps * dcomps);

            auto g_zz_0_xxxxz_0 = cbuffer.data(hs_geom_20_off + 107 * ccomps * dcomps);

            auto g_zz_0_xxxyy_0 = cbuffer.data(hs_geom_20_off + 108 * ccomps * dcomps);

            auto g_zz_0_xxxyz_0 = cbuffer.data(hs_geom_20_off + 109 * ccomps * dcomps);

            auto g_zz_0_xxxzz_0 = cbuffer.data(hs_geom_20_off + 110 * ccomps * dcomps);

            auto g_zz_0_xxyyy_0 = cbuffer.data(hs_geom_20_off + 111 * ccomps * dcomps);

            auto g_zz_0_xxyyz_0 = cbuffer.data(hs_geom_20_off + 112 * ccomps * dcomps);

            auto g_zz_0_xxyzz_0 = cbuffer.data(hs_geom_20_off + 113 * ccomps * dcomps);

            auto g_zz_0_xxzzz_0 = cbuffer.data(hs_geom_20_off + 114 * ccomps * dcomps);

            auto g_zz_0_xyyyy_0 = cbuffer.data(hs_geom_20_off + 115 * ccomps * dcomps);

            auto g_zz_0_xyyyz_0 = cbuffer.data(hs_geom_20_off + 116 * ccomps * dcomps);

            auto g_zz_0_xyyzz_0 = cbuffer.data(hs_geom_20_off + 117 * ccomps * dcomps);

            auto g_zz_0_xyzzz_0 = cbuffer.data(hs_geom_20_off + 118 * ccomps * dcomps);

            auto g_zz_0_xzzzz_0 = cbuffer.data(hs_geom_20_off + 119 * ccomps * dcomps);

            auto g_zz_0_yyyyy_0 = cbuffer.data(hs_geom_20_off + 120 * ccomps * dcomps);

            auto g_zz_0_yyyyz_0 = cbuffer.data(hs_geom_20_off + 121 * ccomps * dcomps);

            auto g_zz_0_yyyzz_0 = cbuffer.data(hs_geom_20_off + 122 * ccomps * dcomps);

            auto g_zz_0_yyzzz_0 = cbuffer.data(hs_geom_20_off + 123 * ccomps * dcomps);

            auto g_zz_0_yzzzz_0 = cbuffer.data(hs_geom_20_off + 124 * ccomps * dcomps);

            auto g_zz_0_zzzzz_0 = cbuffer.data(hs_geom_20_off + 125 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_20_off = idx_geom_20_hpxx + i * dcomps + j;

            auto g_xx_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 149 * ccomps * dcomps);

            auto g_xz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 209 * ccomps * dcomps);

            auto g_yy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 239 * ccomps * dcomps);

            auto g_yy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 240 * ccomps * dcomps);

            auto g_yy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 241 * ccomps * dcomps);

            auto g_yy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 242 * ccomps * dcomps);

            auto g_yy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 243 * ccomps * dcomps);

            auto g_yy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 244 * ccomps * dcomps);

            auto g_yy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 245 * ccomps * dcomps);

            auto g_yy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 246 * ccomps * dcomps);

            auto g_yy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 247 * ccomps * dcomps);

            auto g_yy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 248 * ccomps * dcomps);

            auto g_yy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 249 * ccomps * dcomps);

            auto g_yy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 250 * ccomps * dcomps);

            auto g_yy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 269 * ccomps * dcomps);

            auto g_yz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 299 * ccomps * dcomps);

            auto g_yz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 300 * ccomps * dcomps);

            auto g_yz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 301 * ccomps * dcomps);

            auto g_yz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 302 * ccomps * dcomps);

            auto g_yz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 303 * ccomps * dcomps);

            auto g_yz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 304 * ccomps * dcomps);

            auto g_yz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 305 * ccomps * dcomps);

            auto g_yz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 306 * ccomps * dcomps);

            auto g_yz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 307 * ccomps * dcomps);

            auto g_yz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 308 * ccomps * dcomps);

            auto g_yz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 309 * ccomps * dcomps);

            auto g_yz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 310 * ccomps * dcomps);

            auto g_yz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 311 * ccomps * dcomps);

            auto g_yz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 312 * ccomps * dcomps);

            auto g_yz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 313 * ccomps * dcomps);

            auto g_yz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 329 * ccomps * dcomps);

            auto g_zz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 359 * ccomps * dcomps);

            auto g_zz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 360 * ccomps * dcomps);

            auto g_zz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 361 * ccomps * dcomps);

            auto g_zz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 362 * ccomps * dcomps);

            auto g_zz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 363 * ccomps * dcomps);

            auto g_zz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 364 * ccomps * dcomps);

            auto g_zz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 365 * ccomps * dcomps);

            auto g_zz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 366 * ccomps * dcomps);

            auto g_zz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 367 * ccomps * dcomps);

            auto g_zz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 368 * ccomps * dcomps);

            auto g_zz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 369 * ccomps * dcomps);

            auto g_zz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 370 * ccomps * dcomps);

            auto g_zz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 371 * ccomps * dcomps);

            auto g_zz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 372 * ccomps * dcomps);

            auto g_zz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 373 * ccomps * dcomps);

            auto g_zz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 374 * ccomps * dcomps);

            auto g_zz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 375 * ccomps * dcomps);

            auto g_zz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 376 * ccomps * dcomps);

            auto g_zz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 377 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_isxx

            const auto is_geom_20_off = idx_geom_20_isxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_0, g_xx_0_xxxxx_0, g_xx_0_xxxxx_x, g_xx_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxx_0[k] = -2.0 * g_x_0_xxxxx_0[k] - g_xx_0_xxxxx_0[k] * ab_x + g_xx_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_0, g_xx_0_xxxxx_y, g_xx_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxy_0[k] = -g_xx_0_xxxxx_0[k] * ab_y + g_xx_0_xxxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_0, g_xx_0_xxxxx_z, g_xx_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxz_0[k] = -g_xx_0_xxxxx_0[k] * ab_z + g_xx_0_xxxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxy_0, g_xx_0_xxxxy_y, g_xx_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyy_0[k] = -g_xx_0_xxxxy_0[k] * ab_y + g_xx_0_xxxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxyz_0, g_xx_0_xxxxz_0, g_xx_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyz_0[k] = -g_xx_0_xxxxz_0[k] * ab_y + g_xx_0_xxxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxz_0, g_xx_0_xxxxz_z, g_xx_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxzz_0[k] = -g_xx_0_xxxxz_0[k] * ab_z + g_xx_0_xxxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyy_0, g_xx_0_xxxyy_y, g_xx_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyy_0[k] = -g_xx_0_xxxyy_0[k] * ab_y + g_xx_0_xxxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyyz_0, g_xx_0_xxxyz_0, g_xx_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyz_0[k] = -g_xx_0_xxxyz_0[k] * ab_y + g_xx_0_xxxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyzz_0, g_xx_0_xxxzz_0, g_xx_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyzz_0[k] = -g_xx_0_xxxzz_0[k] * ab_y + g_xx_0_xxxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxzz_0, g_xx_0_xxxzz_z, g_xx_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzzz_0[k] = -g_xx_0_xxxzz_0[k] * ab_z + g_xx_0_xxxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyy_0, g_xx_0_xxyyy_y, g_xx_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyy_0[k] = -g_xx_0_xxyyy_0[k] * ab_y + g_xx_0_xxyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyyz_0, g_xx_0_xxyyz_0, g_xx_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyz_0[k] = -g_xx_0_xxyyz_0[k] * ab_y + g_xx_0_xxyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyzz_0, g_xx_0_xxyzz_0, g_xx_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyzz_0[k] = -g_xx_0_xxyzz_0[k] * ab_y + g_xx_0_xxyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzzz_0, g_xx_0_xxzzz_0, g_xx_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzzz_0[k] = -g_xx_0_xxzzz_0[k] * ab_y + g_xx_0_xxzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzzz_0, g_xx_0_xxzzz_z, g_xx_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzzz_0[k] = -g_xx_0_xxzzz_0[k] * ab_z + g_xx_0_xxzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyy_0, g_xx_0_xyyyy_y, g_xx_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyy_0[k] = -g_xx_0_xyyyy_0[k] * ab_y + g_xx_0_xyyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyyz_0, g_xx_0_xyyyz_0, g_xx_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyz_0[k] = -g_xx_0_xyyyz_0[k] * ab_y + g_xx_0_xyyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyzz_0, g_xx_0_xyyzz_0, g_xx_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyzz_0[k] = -g_xx_0_xyyzz_0[k] * ab_y + g_xx_0_xyyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzzz_0, g_xx_0_xyzzz_0, g_xx_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzzz_0[k] = -g_xx_0_xyzzz_0[k] * ab_y + g_xx_0_xyzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzzz_0, g_xx_0_xzzzz_0, g_xx_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzzz_0[k] = -g_xx_0_xzzzz_0[k] * ab_y + g_xx_0_xzzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzzz_0, g_xx_0_xzzzz_z, g_xx_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzzz_0[k] = -g_xx_0_xzzzz_0[k] * ab_z + g_xx_0_xzzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyy_0, g_xx_0_yyyyy_y, g_xx_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyy_0[k] = -g_xx_0_yyyyy_0[k] * ab_y + g_xx_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyyz_0, g_xx_0_yyyyz_0, g_xx_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyz_0[k] = -g_xx_0_yyyyz_0[k] * ab_y + g_xx_0_yyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyzz_0, g_xx_0_yyyzz_0, g_xx_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyzz_0[k] = -g_xx_0_yyyzz_0[k] * ab_y + g_xx_0_yyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzzz_0, g_xx_0_yyzzz_0, g_xx_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzzz_0[k] = -g_xx_0_yyzzz_0[k] * ab_y + g_xx_0_yyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzzz_0, g_xx_0_yzzzz_0, g_xx_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzzz_0[k] = -g_xx_0_yzzzz_0[k] * ab_y + g_xx_0_yzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzzz_0, g_xx_0_zzzzz_0, g_xx_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzzz_0[k] = -g_xx_0_zzzzz_0[k] * ab_y + g_xx_0_zzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzzz_0, g_xx_0_zzzzz_z, g_xx_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzzz_0[k] = -g_xx_0_zzzzz_0[k] * ab_z + g_xx_0_zzzzz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_0, g_xy_0_xxxxx_x, g_xy_0_xxxxxx_0, g_y_0_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxx_0[k] = -g_y_0_xxxxx_0[k] - g_xy_0_xxxxx_0[k] * ab_x + g_xy_0_xxxxx_x[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxxy_0, g_xy_0_xxxxy_0, g_xy_0_xxxxy_x, g_y_0_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxy_0[k] = -g_y_0_xxxxy_0[k] - g_xy_0_xxxxy_0[k] * ab_x + g_xy_0_xxxxy_x[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_0, g_xy_0_xxxxx_z, g_xy_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxz_0[k] = -g_xy_0_xxxxx_0[k] * ab_z + g_xy_0_xxxxx_z[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxyy_0, g_xy_0_xxxyy_0, g_xy_0_xxxyy_x, g_y_0_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyy_0[k] = -g_y_0_xxxyy_0[k] - g_xy_0_xxxyy_0[k] * ab_x + g_xy_0_xxxyy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_0, g_xy_0_xxxxy_z, g_xy_0_xxxxyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyz_0[k] = -g_xy_0_xxxxy_0[k] * ab_z + g_xy_0_xxxxy_z[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxz_0, g_xy_0_xxxxz_z, g_xy_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxzz_0[k] = -g_xy_0_xxxxz_0[k] * ab_z + g_xy_0_xxxxz_z[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyyy_0, g_xy_0_xxyyy_0, g_xy_0_xxyyy_x, g_y_0_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyy_0[k] = -g_y_0_xxyyy_0[k] - g_xy_0_xxyyy_0[k] * ab_x + g_xy_0_xxyyy_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_0, g_xy_0_xxxyy_z, g_xy_0_xxxyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyz_0[k] = -g_xy_0_xxxyy_0[k] * ab_z + g_xy_0_xxxyy_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyz_0, g_xy_0_xxxyz_z, g_xy_0_xxxyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyzz_0[k] = -g_xy_0_xxxyz_0[k] * ab_z + g_xy_0_xxxyz_z[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxzz_0, g_xy_0_xxxzz_z, g_xy_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzzz_0[k] = -g_xy_0_xxxzz_0[k] * ab_z + g_xy_0_xxxzz_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyyy_0, g_xy_0_xyyyy_0, g_xy_0_xyyyy_x, g_y_0_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyy_0[k] = -g_y_0_xyyyy_0[k] - g_xy_0_xyyyy_0[k] * ab_x + g_xy_0_xyyyy_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_0, g_xy_0_xxyyy_z, g_xy_0_xxyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyz_0[k] = -g_xy_0_xxyyy_0[k] * ab_z + g_xy_0_xxyyy_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyz_0, g_xy_0_xxyyz_z, g_xy_0_xxyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyzz_0[k] = -g_xy_0_xxyyz_0[k] * ab_z + g_xy_0_xxyyz_z[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyzz_0, g_xy_0_xxyzz_z, g_xy_0_xxyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzzz_0[k] = -g_xy_0_xxyzz_0[k] * ab_z + g_xy_0_xxyzz_z[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzzz_0, g_xy_0_xxzzz_z, g_xy_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzzz_0[k] = -g_xy_0_xxzzz_0[k] * ab_z + g_xy_0_xxzzz_z[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyyy_0, g_xy_0_yyyyy_0, g_xy_0_yyyyy_x, g_y_0_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyy_0[k] = -g_y_0_yyyyy_0[k] - g_xy_0_yyyyy_0[k] * ab_x + g_xy_0_yyyyy_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_0, g_xy_0_xyyyy_z, g_xy_0_xyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyz_0[k] = -g_xy_0_xyyyy_0[k] * ab_z + g_xy_0_xyyyy_z[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyz_0, g_xy_0_xyyyz_z, g_xy_0_xyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyzz_0[k] = -g_xy_0_xyyyz_0[k] * ab_z + g_xy_0_xyyyz_z[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyzz_0, g_xy_0_xyyzz_z, g_xy_0_xyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzzz_0[k] = -g_xy_0_xyyzz_0[k] * ab_z + g_xy_0_xyyzz_z[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzzz_0, g_xy_0_xyzzz_z, g_xy_0_xyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzzz_0[k] = -g_xy_0_xyzzz_0[k] * ab_z + g_xy_0_xyzzz_z[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzzz_0, g_xy_0_xzzzz_z, g_xy_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzzz_0[k] = -g_xy_0_xzzzz_0[k] * ab_z + g_xy_0_xzzzz_z[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_0, g_xy_0_yyyyy_0, g_xy_0_yyyyy_y, g_xy_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyy_0[k] = -g_x_0_yyyyy_0[k] - g_xy_0_yyyyy_0[k] * ab_y + g_xy_0_yyyyy_y[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyy_0, g_xy_0_yyyyy_z, g_xy_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyz_0[k] = -g_xy_0_yyyyy_0[k] * ab_z + g_xy_0_yyyyy_z[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyz_0, g_xy_0_yyyyz_z, g_xy_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyzz_0[k] = -g_xy_0_yyyyz_0[k] * ab_z + g_xy_0_yyyyz_z[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyzz_0, g_xy_0_yyyzz_z, g_xy_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzzz_0[k] = -g_xy_0_yyyzz_0[k] * ab_z + g_xy_0_yyyzz_z[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzzz_0, g_xy_0_yyzzz_z, g_xy_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzzz_0[k] = -g_xy_0_yyzzz_0[k] * ab_z + g_xy_0_yyzzz_z[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzzz_0, g_xy_0_yzzzz_z, g_xy_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzzz_0[k] = -g_xy_0_yzzzz_0[k] * ab_z + g_xy_0_yzzzz_z[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzzz_0, g_xy_0_zzzzz_z, g_xy_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzzz_0[k] = -g_xy_0_zzzzz_0[k] * ab_z + g_xy_0_zzzzz_z[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_0, g_xz_0_xxxxx_x, g_xz_0_xxxxxx_0, g_z_0_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxx_0[k] = -g_z_0_xxxxx_0[k] - g_xz_0_xxxxx_0[k] * ab_x + g_xz_0_xxxxx_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_0, g_xz_0_xxxxx_y, g_xz_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxy_0[k] = -g_xz_0_xxxxx_0[k] * ab_y + g_xz_0_xxxxx_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxxz_0, g_xz_0_xxxxz_0, g_xz_0_xxxxz_x, g_z_0_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxz_0[k] = -g_z_0_xxxxz_0[k] - g_xz_0_xxxxz_0[k] * ab_x + g_xz_0_xxxxz_x[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxy_0, g_xz_0_xxxxy_y, g_xz_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyy_0[k] = -g_xz_0_xxxxy_0[k] * ab_y + g_xz_0_xxxxy_y[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxyz_0, g_xz_0_xxxxz_0, g_xz_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyz_0[k] = -g_xz_0_xxxxz_0[k] * ab_y + g_xz_0_xxxxz_y[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxzz_0, g_xz_0_xxxzz_0, g_xz_0_xxxzz_x, g_z_0_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxzz_0[k] = -g_z_0_xxxzz_0[k] - g_xz_0_xxxzz_0[k] * ab_x + g_xz_0_xxxzz_x[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyy_0, g_xz_0_xxxyy_y, g_xz_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyy_0[k] = -g_xz_0_xxxyy_0[k] * ab_y + g_xz_0_xxxyy_y[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyyz_0, g_xz_0_xxxyz_0, g_xz_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyz_0[k] = -g_xz_0_xxxyz_0[k] * ab_y + g_xz_0_xxxyz_y[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyzz_0, g_xz_0_xxxzz_0, g_xz_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyzz_0[k] = -g_xz_0_xxxzz_0[k] * ab_y + g_xz_0_xxxzz_y[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzzz_0, g_xz_0_xxzzz_0, g_xz_0_xxzzz_x, g_z_0_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzzz_0[k] = -g_z_0_xxzzz_0[k] - g_xz_0_xxzzz_0[k] * ab_x + g_xz_0_xxzzz_x[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyy_0, g_xz_0_xxyyy_y, g_xz_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyy_0[k] = -g_xz_0_xxyyy_0[k] * ab_y + g_xz_0_xxyyy_y[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyyz_0, g_xz_0_xxyyz_0, g_xz_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyz_0[k] = -g_xz_0_xxyyz_0[k] * ab_y + g_xz_0_xxyyz_y[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyzz_0, g_xz_0_xxyzz_0, g_xz_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyzz_0[k] = -g_xz_0_xxyzz_0[k] * ab_y + g_xz_0_xxyzz_y[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzzz_0, g_xz_0_xxzzz_0, g_xz_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzzz_0[k] = -g_xz_0_xxzzz_0[k] * ab_y + g_xz_0_xxzzz_y[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzzz_0, g_xz_0_xzzzz_0, g_xz_0_xzzzz_x, g_z_0_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzzz_0[k] = -g_z_0_xzzzz_0[k] - g_xz_0_xzzzz_0[k] * ab_x + g_xz_0_xzzzz_x[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyy_0, g_xz_0_xyyyy_y, g_xz_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyy_0[k] = -g_xz_0_xyyyy_0[k] * ab_y + g_xz_0_xyyyy_y[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyyz_0, g_xz_0_xyyyz_0, g_xz_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyz_0[k] = -g_xz_0_xyyyz_0[k] * ab_y + g_xz_0_xyyyz_y[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyzz_0, g_xz_0_xyyzz_0, g_xz_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyzz_0[k] = -g_xz_0_xyyzz_0[k] * ab_y + g_xz_0_xyyzz_y[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzzz_0, g_xz_0_xyzzz_0, g_xz_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzzz_0[k] = -g_xz_0_xyzzz_0[k] * ab_y + g_xz_0_xyzzz_y[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzzz_0, g_xz_0_xzzzz_0, g_xz_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzzz_0[k] = -g_xz_0_xzzzz_0[k] * ab_y + g_xz_0_xzzzz_y[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzzz_0, g_xz_0_zzzzz_0, g_xz_0_zzzzz_x, g_z_0_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzzz_0[k] = -g_z_0_zzzzz_0[k] - g_xz_0_zzzzz_0[k] * ab_x + g_xz_0_zzzzz_x[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyy_0, g_xz_0_yyyyy_y, g_xz_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyy_0[k] = -g_xz_0_yyyyy_0[k] * ab_y + g_xz_0_yyyyy_y[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyyz_0, g_xz_0_yyyyz_0, g_xz_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyz_0[k] = -g_xz_0_yyyyz_0[k] * ab_y + g_xz_0_yyyyz_y[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyzz_0, g_xz_0_yyyzz_0, g_xz_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyzz_0[k] = -g_xz_0_yyyzz_0[k] * ab_y + g_xz_0_yyyzz_y[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzzz_0, g_xz_0_yyzzz_0, g_xz_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzzz_0[k] = -g_xz_0_yyzzz_0[k] * ab_y + g_xz_0_yyzzz_y[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzzz_0, g_xz_0_yzzzz_0, g_xz_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzzz_0[k] = -g_xz_0_yzzzz_0[k] * ab_y + g_xz_0_yzzzz_y[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzzz_0, g_xz_0_zzzzz_0, g_xz_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzzz_0[k] = -g_xz_0_zzzzz_0[k] * ab_y + g_xz_0_zzzzz_y[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_0, g_xz_0_zzzzz_0, g_xz_0_zzzzz_z, g_xz_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzzz_0[k] = -g_x_0_zzzzz_0[k] - g_xz_0_zzzzz_0[k] * ab_z + g_xz_0_zzzzz_z[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxx_0, g_yy_0_xxxxx_x, g_yy_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxx_0[k] = -g_yy_0_xxxxx_0[k] * ab_x + g_yy_0_xxxxx_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxy_0, g_yy_0_xxxxy_0, g_yy_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxy_0[k] = -g_yy_0_xxxxy_0[k] * ab_x + g_yy_0_xxxxy_x[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxz_0, g_yy_0_xxxxz_0, g_yy_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxz_0[k] = -g_yy_0_xxxxz_0[k] * ab_x + g_yy_0_xxxxz_x[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyy_0, g_yy_0_xxxyy_0, g_yy_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyy_0[k] = -g_yy_0_xxxyy_0[k] * ab_x + g_yy_0_xxxyy_x[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyz_0, g_yy_0_xxxyz_0, g_yy_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyz_0[k] = -g_yy_0_xxxyz_0[k] * ab_x + g_yy_0_xxxyz_x[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxzz_0, g_yy_0_xxxzz_0, g_yy_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxzz_0[k] = -g_yy_0_xxxzz_0[k] * ab_x + g_yy_0_xxxzz_x[k];
            }

            /// Set up 90-91 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 90 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyy_0, g_yy_0_xxyyy_0, g_yy_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyy_0[k] = -g_yy_0_xxyyy_0[k] * ab_x + g_yy_0_xxyyy_x[k];
            }

            /// Set up 91-92 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 91 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyz_0, g_yy_0_xxyyz_0, g_yy_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyz_0[k] = -g_yy_0_xxyyz_0[k] * ab_x + g_yy_0_xxyyz_x[k];
            }

            /// Set up 92-93 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyzz_0, g_yy_0_xxyzz_0, g_yy_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyzz_0[k] = -g_yy_0_xxyzz_0[k] * ab_x + g_yy_0_xxyzz_x[k];
            }

            /// Set up 93-94 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 93 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzzz_0, g_yy_0_xxzzz_0, g_yy_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzzz_0[k] = -g_yy_0_xxzzz_0[k] * ab_x + g_yy_0_xxzzz_x[k];
            }

            /// Set up 94-95 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 94 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyy_0, g_yy_0_xyyyy_0, g_yy_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyy_0[k] = -g_yy_0_xyyyy_0[k] * ab_x + g_yy_0_xyyyy_x[k];
            }

            /// Set up 95-96 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyz_0, g_yy_0_xyyyz_0, g_yy_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyz_0[k] = -g_yy_0_xyyyz_0[k] * ab_x + g_yy_0_xyyyz_x[k];
            }

            /// Set up 96-97 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 96 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyzz_0, g_yy_0_xyyzz_0, g_yy_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyzz_0[k] = -g_yy_0_xyyzz_0[k] * ab_x + g_yy_0_xyyzz_x[k];
            }

            /// Set up 97-98 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 97 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzzz_0, g_yy_0_xyzzz_0, g_yy_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzzz_0[k] = -g_yy_0_xyzzz_0[k] * ab_x + g_yy_0_xyzzz_x[k];
            }

            /// Set up 98-99 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzzz_0, g_yy_0_xzzzz_0, g_yy_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzzz_0[k] = -g_yy_0_xzzzz_0[k] * ab_x + g_yy_0_xzzzz_x[k];
            }

            /// Set up 99-100 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyy_0, g_yy_0_yyyyy_0, g_yy_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyy_0[k] = -g_yy_0_yyyyy_0[k] * ab_x + g_yy_0_yyyyy_x[k];
            }

            /// Set up 100-101 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 100 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyz_0, g_yy_0_yyyyz_0, g_yy_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyz_0[k] = -g_yy_0_yyyyz_0[k] * ab_x + g_yy_0_yyyyz_x[k];
            }

            /// Set up 101-102 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyzz_0, g_yy_0_yyyzz_0, g_yy_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyzz_0[k] = -g_yy_0_yyyzz_0[k] * ab_x + g_yy_0_yyyzz_x[k];
            }

            /// Set up 102-103 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 102 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzzz_0, g_yy_0_yyzzz_0, g_yy_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzzz_0[k] = -g_yy_0_yyzzz_0[k] * ab_x + g_yy_0_yyzzz_x[k];
            }

            /// Set up 103-104 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 103 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzzz_0, g_yy_0_yzzzz_0, g_yy_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzzz_0[k] = -g_yy_0_yzzzz_0[k] * ab_x + g_yy_0_yzzzz_x[k];
            }

            /// Set up 104-105 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzzz_0, g_yy_0_zzzzz_0, g_yy_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzzz_0[k] = -g_yy_0_zzzzz_0[k] * ab_x + g_yy_0_zzzzz_x[k];
            }

            /// Set up 105-106 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 105 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_0, g_yy_0_yyyyy_0, g_yy_0_yyyyy_y, g_yy_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyy_0[k] = -2.0 * g_y_0_yyyyy_0[k] - g_yy_0_yyyyy_0[k] * ab_y + g_yy_0_yyyyy_y[k];
            }

            /// Set up 106-107 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 106 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyy_0, g_yy_0_yyyyy_z, g_yy_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyz_0[k] = -g_yy_0_yyyyy_0[k] * ab_z + g_yy_0_yyyyy_z[k];
            }

            /// Set up 107-108 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyz_0, g_yy_0_yyyyz_z, g_yy_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyzz_0[k] = -g_yy_0_yyyyz_0[k] * ab_z + g_yy_0_yyyyz_z[k];
            }

            /// Set up 108-109 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 108 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyzz_0, g_yy_0_yyyzz_z, g_yy_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzzz_0[k] = -g_yy_0_yyyzz_0[k] * ab_z + g_yy_0_yyyzz_z[k];
            }

            /// Set up 109-110 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzzz_0, g_yy_0_yyzzz_z, g_yy_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzzz_0[k] = -g_yy_0_yyzzz_0[k] * ab_z + g_yy_0_yyzzz_z[k];
            }

            /// Set up 110-111 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzzz_0, g_yy_0_yzzzz_z, g_yy_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzzz_0[k] = -g_yy_0_yzzzz_0[k] * ab_z + g_yy_0_yzzzz_z[k];
            }

            /// Set up 111-112 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzzz_0, g_yy_0_zzzzz_z, g_yy_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzzz_0[k] = -g_yy_0_zzzzz_0[k] * ab_z + g_yy_0_zzzzz_z[k];
            }

            /// Set up 112-113 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 112 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxx_0, g_yz_0_xxxxx_x, g_yz_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxx_0[k] = -g_yz_0_xxxxx_0[k] * ab_x + g_yz_0_xxxxx_x[k];
            }

            /// Set up 113-114 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxy_0, g_yz_0_xxxxy_0, g_yz_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxy_0[k] = -g_yz_0_xxxxy_0[k] * ab_x + g_yz_0_xxxxy_x[k];
            }

            /// Set up 114-115 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 114 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxz_0, g_yz_0_xxxxz_0, g_yz_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxz_0[k] = -g_yz_0_xxxxz_0[k] * ab_x + g_yz_0_xxxxz_x[k];
            }

            /// Set up 115-116 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 115 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyy_0, g_yz_0_xxxyy_0, g_yz_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyy_0[k] = -g_yz_0_xxxyy_0[k] * ab_x + g_yz_0_xxxyy_x[k];
            }

            /// Set up 116-117 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyz_0, g_yz_0_xxxyz_0, g_yz_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyz_0[k] = -g_yz_0_xxxyz_0[k] * ab_x + g_yz_0_xxxyz_x[k];
            }

            /// Set up 117-118 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 117 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxzz_0, g_yz_0_xxxzz_0, g_yz_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxzz_0[k] = -g_yz_0_xxxzz_0[k] * ab_x + g_yz_0_xxxzz_x[k];
            }

            /// Set up 118-119 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 118 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyy_0, g_yz_0_xxyyy_0, g_yz_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyy_0[k] = -g_yz_0_xxyyy_0[k] * ab_x + g_yz_0_xxyyy_x[k];
            }

            /// Set up 119-120 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyz_0, g_yz_0_xxyyz_0, g_yz_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyz_0[k] = -g_yz_0_xxyyz_0[k] * ab_x + g_yz_0_xxyyz_x[k];
            }

            /// Set up 120-121 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 120 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyzz_0, g_yz_0_xxyzz_0, g_yz_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyzz_0[k] = -g_yz_0_xxyzz_0[k] * ab_x + g_yz_0_xxyzz_x[k];
            }

            /// Set up 121-122 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 121 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzzz_0, g_yz_0_xxzzz_0, g_yz_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzzz_0[k] = -g_yz_0_xxzzz_0[k] * ab_x + g_yz_0_xxzzz_x[k];
            }

            /// Set up 122-123 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyy_0, g_yz_0_xyyyy_0, g_yz_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyy_0[k] = -g_yz_0_xyyyy_0[k] * ab_x + g_yz_0_xyyyy_x[k];
            }

            /// Set up 123-124 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 123 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyz_0, g_yz_0_xyyyz_0, g_yz_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyz_0[k] = -g_yz_0_xyyyz_0[k] * ab_x + g_yz_0_xyyyz_x[k];
            }

            /// Set up 124-125 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 124 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyzz_0, g_yz_0_xyyzz_0, g_yz_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyzz_0[k] = -g_yz_0_xyyzz_0[k] * ab_x + g_yz_0_xyyzz_x[k];
            }

            /// Set up 125-126 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzzz_0, g_yz_0_xyzzz_0, g_yz_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzzz_0[k] = -g_yz_0_xyzzz_0[k] * ab_x + g_yz_0_xyzzz_x[k];
            }

            /// Set up 126-127 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 126 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzzz_0, g_yz_0_xzzzz_0, g_yz_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzzz_0[k] = -g_yz_0_xzzzz_0[k] * ab_x + g_yz_0_xzzzz_x[k];
            }

            /// Set up 127-128 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 127 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyy_0, g_yz_0_yyyyy_0, g_yz_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyy_0[k] = -g_yz_0_yyyyy_0[k] * ab_x + g_yz_0_yyyyy_x[k];
            }

            /// Set up 128-129 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyz_0, g_yz_0_yyyyz_0, g_yz_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyz_0[k] = -g_yz_0_yyyyz_0[k] * ab_x + g_yz_0_yyyyz_x[k];
            }

            /// Set up 129-130 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyzz_0, g_yz_0_yyyzz_0, g_yz_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyzz_0[k] = -g_yz_0_yyyzz_0[k] * ab_x + g_yz_0_yyyzz_x[k];
            }

            /// Set up 130-131 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 130 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzzz_0, g_yz_0_yyzzz_0, g_yz_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzzz_0[k] = -g_yz_0_yyzzz_0[k] * ab_x + g_yz_0_yyzzz_x[k];
            }

            /// Set up 131-132 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzzz_0, g_yz_0_yzzzz_0, g_yz_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzzz_0[k] = -g_yz_0_yzzzz_0[k] * ab_x + g_yz_0_yzzzz_x[k];
            }

            /// Set up 132-133 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 132 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzzz_0, g_yz_0_zzzzz_0, g_yz_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzzz_0[k] = -g_yz_0_zzzzz_0[k] * ab_x + g_yz_0_zzzzz_x[k];
            }

            /// Set up 133-134 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 133 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyy_0, g_yz_0_yyyyy_y, g_yz_0_yyyyyy_0, g_z_0_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyy_0[k] = -g_z_0_yyyyy_0[k] - g_yz_0_yyyyy_0[k] * ab_y + g_yz_0_yyyyy_y[k];
            }

            /// Set up 134-135 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyyz_0, g_yz_0_yyyyz_0, g_yz_0_yyyyz_y, g_z_0_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyz_0[k] = -g_z_0_yyyyz_0[k] - g_yz_0_yyyyz_0[k] * ab_y + g_yz_0_yyyyz_y[k];
            }

            /// Set up 135-136 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 135 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyzz_0, g_yz_0_yyyzz_0, g_yz_0_yyyzz_y, g_z_0_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyzz_0[k] = -g_z_0_yyyzz_0[k] - g_yz_0_yyyzz_0[k] * ab_y + g_yz_0_yyyzz_y[k];
            }

            /// Set up 136-137 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 136 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzzz_0, g_yz_0_yyzzz_0, g_yz_0_yyzzz_y, g_z_0_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzzz_0[k] = -g_z_0_yyzzz_0[k] - g_yz_0_yyzzz_0[k] * ab_y + g_yz_0_yyzzz_y[k];
            }

            /// Set up 137-138 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzzz_0, g_yz_0_yzzzz_0, g_yz_0_yzzzz_y, g_z_0_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzzz_0[k] = -g_z_0_yzzzz_0[k] - g_yz_0_yzzzz_0[k] * ab_y + g_yz_0_yzzzz_y[k];
            }

            /// Set up 138-139 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 138 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzzz_0, g_yz_0_zzzzz_0, g_yz_0_zzzzz_y, g_z_0_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzzz_0[k] = -g_z_0_zzzzz_0[k] - g_yz_0_zzzzz_0[k] * ab_y + g_yz_0_zzzzz_y[k];
            }

            /// Set up 139-140 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_0, g_yz_0_zzzzz_0, g_yz_0_zzzzz_z, g_yz_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzzz_0[k] = -g_y_0_zzzzz_0[k] - g_yz_0_zzzzz_0[k] * ab_z + g_yz_0_zzzzz_z[k];
            }

            /// Set up 140-141 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxx_0 = cbuffer.data(is_geom_20_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxx_0, g_zz_0_xxxxx_x, g_zz_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxx_0[k] = -g_zz_0_xxxxx_0[k] * ab_x + g_zz_0_xxxxx_x[k];
            }

            /// Set up 141-142 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxy_0 = cbuffer.data(is_geom_20_off + 141 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxy_0, g_zz_0_xxxxy_0, g_zz_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxy_0[k] = -g_zz_0_xxxxy_0[k] * ab_x + g_zz_0_xxxxy_x[k];
            }

            /// Set up 142-143 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxz_0 = cbuffer.data(is_geom_20_off + 142 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxz_0, g_zz_0_xxxxz_0, g_zz_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxz_0[k] = -g_zz_0_xxxxz_0[k] * ab_x + g_zz_0_xxxxz_x[k];
            }

            /// Set up 143-144 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyy_0 = cbuffer.data(is_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyy_0, g_zz_0_xxxyy_0, g_zz_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyy_0[k] = -g_zz_0_xxxyy_0[k] * ab_x + g_zz_0_xxxyy_x[k];
            }

            /// Set up 144-145 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyz_0 = cbuffer.data(is_geom_20_off + 144 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyz_0, g_zz_0_xxxyz_0, g_zz_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyz_0[k] = -g_zz_0_xxxyz_0[k] * ab_x + g_zz_0_xxxyz_x[k];
            }

            /// Set up 145-146 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxzz_0 = cbuffer.data(is_geom_20_off + 145 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxzz_0, g_zz_0_xxxzz_0, g_zz_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxzz_0[k] = -g_zz_0_xxxzz_0[k] * ab_x + g_zz_0_xxxzz_x[k];
            }

            /// Set up 146-147 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyy_0 = cbuffer.data(is_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyy_0, g_zz_0_xxyyy_0, g_zz_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyy_0[k] = -g_zz_0_xxyyy_0[k] * ab_x + g_zz_0_xxyyy_x[k];
            }

            /// Set up 147-148 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyz_0 = cbuffer.data(is_geom_20_off + 147 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyz_0, g_zz_0_xxyyz_0, g_zz_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyz_0[k] = -g_zz_0_xxyyz_0[k] * ab_x + g_zz_0_xxyyz_x[k];
            }

            /// Set up 148-149 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyzz_0 = cbuffer.data(is_geom_20_off + 148 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyzz_0, g_zz_0_xxyzz_0, g_zz_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyzz_0[k] = -g_zz_0_xxyzz_0[k] * ab_x + g_zz_0_xxyzz_x[k];
            }

            /// Set up 149-150 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzzz_0 = cbuffer.data(is_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzzz_0, g_zz_0_xxzzz_0, g_zz_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzzz_0[k] = -g_zz_0_xxzzz_0[k] * ab_x + g_zz_0_xxzzz_x[k];
            }

            /// Set up 150-151 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyy_0 = cbuffer.data(is_geom_20_off + 150 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyy_0, g_zz_0_xyyyy_0, g_zz_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyy_0[k] = -g_zz_0_xyyyy_0[k] * ab_x + g_zz_0_xyyyy_x[k];
            }

            /// Set up 151-152 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyz_0 = cbuffer.data(is_geom_20_off + 151 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyz_0, g_zz_0_xyyyz_0, g_zz_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyz_0[k] = -g_zz_0_xyyyz_0[k] * ab_x + g_zz_0_xyyyz_x[k];
            }

            /// Set up 152-153 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyzz_0 = cbuffer.data(is_geom_20_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyzz_0, g_zz_0_xyyzz_0, g_zz_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyzz_0[k] = -g_zz_0_xyyzz_0[k] * ab_x + g_zz_0_xyyzz_x[k];
            }

            /// Set up 153-154 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzzz_0 = cbuffer.data(is_geom_20_off + 153 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzzz_0, g_zz_0_xyzzz_0, g_zz_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzzz_0[k] = -g_zz_0_xyzzz_0[k] * ab_x + g_zz_0_xyzzz_x[k];
            }

            /// Set up 154-155 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzzz_0 = cbuffer.data(is_geom_20_off + 154 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzzz_0, g_zz_0_xzzzz_0, g_zz_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzzz_0[k] = -g_zz_0_xzzzz_0[k] * ab_x + g_zz_0_xzzzz_x[k];
            }

            /// Set up 155-156 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyy_0 = cbuffer.data(is_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyy_0, g_zz_0_yyyyy_0, g_zz_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyy_0[k] = -g_zz_0_yyyyy_0[k] * ab_x + g_zz_0_yyyyy_x[k];
            }

            /// Set up 156-157 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyz_0 = cbuffer.data(is_geom_20_off + 156 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyz_0, g_zz_0_yyyyz_0, g_zz_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyz_0[k] = -g_zz_0_yyyyz_0[k] * ab_x + g_zz_0_yyyyz_x[k];
            }

            /// Set up 157-158 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyzz_0 = cbuffer.data(is_geom_20_off + 157 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyzz_0, g_zz_0_yyyzz_0, g_zz_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyzz_0[k] = -g_zz_0_yyyzz_0[k] * ab_x + g_zz_0_yyyzz_x[k];
            }

            /// Set up 158-159 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzzz_0 = cbuffer.data(is_geom_20_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzzz_0, g_zz_0_yyzzz_0, g_zz_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzzz_0[k] = -g_zz_0_yyzzz_0[k] * ab_x + g_zz_0_yyzzz_x[k];
            }

            /// Set up 159-160 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzzz_0 = cbuffer.data(is_geom_20_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzzz_0, g_zz_0_yzzzz_0, g_zz_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzzz_0[k] = -g_zz_0_yzzzz_0[k] * ab_x + g_zz_0_yzzzz_x[k];
            }

            /// Set up 160-161 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzzz_0 = cbuffer.data(is_geom_20_off + 160 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzzz_0, g_zz_0_zzzzz_0, g_zz_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzzz_0[k] = -g_zz_0_zzzzz_0[k] * ab_x + g_zz_0_zzzzz_x[k];
            }

            /// Set up 161-162 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyy_0 = cbuffer.data(is_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyy_0, g_zz_0_yyyyy_y, g_zz_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyy_0[k] = -g_zz_0_yyyyy_0[k] * ab_y + g_zz_0_yyyyy_y[k];
            }

            /// Set up 162-163 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyz_0 = cbuffer.data(is_geom_20_off + 162 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyyz_0, g_zz_0_yyyyz_0, g_zz_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyz_0[k] = -g_zz_0_yyyyz_0[k] * ab_y + g_zz_0_yyyyz_y[k];
            }

            /// Set up 163-164 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyzz_0 = cbuffer.data(is_geom_20_off + 163 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyzz_0, g_zz_0_yyyzz_0, g_zz_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyzz_0[k] = -g_zz_0_yyyzz_0[k] * ab_y + g_zz_0_yyyzz_y[k];
            }

            /// Set up 164-165 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzzz_0 = cbuffer.data(is_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzzz_0, g_zz_0_yyzzz_0, g_zz_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzzz_0[k] = -g_zz_0_yyzzz_0[k] * ab_y + g_zz_0_yyzzz_y[k];
            }

            /// Set up 165-166 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzzz_0 = cbuffer.data(is_geom_20_off + 165 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzzz_0, g_zz_0_yzzzz_0, g_zz_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzzz_0[k] = -g_zz_0_yzzzz_0[k] * ab_y + g_zz_0_yzzzz_y[k];
            }

            /// Set up 166-167 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzzz_0 = cbuffer.data(is_geom_20_off + 166 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzzz_0, g_zz_0_zzzzz_0, g_zz_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzzz_0[k] = -g_zz_0_zzzzz_0[k] * ab_y + g_zz_0_zzzzz_y[k];
            }

            /// Set up 167-168 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzzz_0 = cbuffer.data(is_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_0, g_zz_0_zzzzz_0, g_zz_0_zzzzz_z, g_zz_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzzz_0[k] = -2.0 * g_z_0_zzzzz_0[k] - g_zz_0_zzzzz_0[k] * ab_z + g_zz_0_zzzzz_z[k];
            }
        }
    }
}

} // erirec namespace

