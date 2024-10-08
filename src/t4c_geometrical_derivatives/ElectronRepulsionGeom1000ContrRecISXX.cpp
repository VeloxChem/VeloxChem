#include "ElectronRepulsionGeom1000ContrRecISXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_isxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_isxx,
                                            const size_t idx_hsxx,
                                            const size_t idx_geom_10_hsxx,
                                            const size_t idx_geom_10_hpxx,
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

            const auto hs_off = idx_hsxx + i * dcomps + j;

            auto g_xxxxx_0 = cbuffer.data(hs_off + 0 * ccomps * dcomps);

            auto g_xxxxy_0 = cbuffer.data(hs_off + 1 * ccomps * dcomps);

            auto g_xxxxz_0 = cbuffer.data(hs_off + 2 * ccomps * dcomps);

            auto g_xxxyy_0 = cbuffer.data(hs_off + 3 * ccomps * dcomps);

            auto g_xxxyz_0 = cbuffer.data(hs_off + 4 * ccomps * dcomps);

            auto g_xxxzz_0 = cbuffer.data(hs_off + 5 * ccomps * dcomps);

            auto g_xxyyy_0 = cbuffer.data(hs_off + 6 * ccomps * dcomps);

            auto g_xxyyz_0 = cbuffer.data(hs_off + 7 * ccomps * dcomps);

            auto g_xxyzz_0 = cbuffer.data(hs_off + 8 * ccomps * dcomps);

            auto g_xxzzz_0 = cbuffer.data(hs_off + 9 * ccomps * dcomps);

            auto g_xyyyy_0 = cbuffer.data(hs_off + 10 * ccomps * dcomps);

            auto g_xyyyz_0 = cbuffer.data(hs_off + 11 * ccomps * dcomps);

            auto g_xyyzz_0 = cbuffer.data(hs_off + 12 * ccomps * dcomps);

            auto g_xyzzz_0 = cbuffer.data(hs_off + 13 * ccomps * dcomps);

            auto g_xzzzz_0 = cbuffer.data(hs_off + 14 * ccomps * dcomps);

            auto g_yyyyy_0 = cbuffer.data(hs_off + 15 * ccomps * dcomps);

            auto g_yyyyz_0 = cbuffer.data(hs_off + 16 * ccomps * dcomps);

            auto g_yyyzz_0 = cbuffer.data(hs_off + 17 * ccomps * dcomps);

            auto g_yyzzz_0 = cbuffer.data(hs_off + 18 * ccomps * dcomps);

            auto g_yzzzz_0 = cbuffer.data(hs_off + 19 * ccomps * dcomps);

            auto g_zzzzz_0 = cbuffer.data(hs_off + 20 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_10_off = idx_geom_10_hpxx + i * dcomps + j;

            auto g_x_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 188 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_isxx

            const auto is_geom_10_off = idx_geom_10_isxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_0, g_x_0_xxxxx_x, g_x_0_xxxxxx_0, g_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_0[k] = -g_xxxxx_0[k] - g_x_0_xxxxx_0[k] * ab_x + g_x_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_0, g_x_0_xxxxx_y, g_x_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_0[k] = -g_x_0_xxxxx_0[k] * ab_y + g_x_0_xxxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_0, g_x_0_xxxxx_z, g_x_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_0[k] = -g_x_0_xxxxx_0[k] * ab_z + g_x_0_xxxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxy_0, g_x_0_xxxxy_y, g_x_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_0[k] = -g_x_0_xxxxy_0[k] * ab_y + g_x_0_xxxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxyz_0, g_x_0_xxxxz_0, g_x_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_0[k] = -g_x_0_xxxxz_0[k] * ab_y + g_x_0_xxxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxz_0, g_x_0_xxxxz_z, g_x_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_0[k] = -g_x_0_xxxxz_0[k] * ab_z + g_x_0_xxxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyy_0, g_x_0_xxxyy_y, g_x_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_0[k] = -g_x_0_xxxyy_0[k] * ab_y + g_x_0_xxxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyyz_0, g_x_0_xxxyz_0, g_x_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_0[k] = -g_x_0_xxxyz_0[k] * ab_y + g_x_0_xxxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyzz_0, g_x_0_xxxzz_0, g_x_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_0[k] = -g_x_0_xxxzz_0[k] * ab_y + g_x_0_xxxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxzz_0, g_x_0_xxxzz_z, g_x_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_0[k] = -g_x_0_xxxzz_0[k] * ab_z + g_x_0_xxxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyy_0, g_x_0_xxyyy_y, g_x_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_0[k] = -g_x_0_xxyyy_0[k] * ab_y + g_x_0_xxyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyyz_0, g_x_0_xxyyz_0, g_x_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_0[k] = -g_x_0_xxyyz_0[k] * ab_y + g_x_0_xxyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyzz_0, g_x_0_xxyzz_0, g_x_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_0[k] = -g_x_0_xxyzz_0[k] * ab_y + g_x_0_xxyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzzz_0, g_x_0_xxzzz_0, g_x_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_0[k] = -g_x_0_xxzzz_0[k] * ab_y + g_x_0_xxzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzzz_0, g_x_0_xxzzz_z, g_x_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_0[k] = -g_x_0_xxzzz_0[k] * ab_z + g_x_0_xxzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyy_0, g_x_0_xyyyy_y, g_x_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_0[k] = -g_x_0_xyyyy_0[k] * ab_y + g_x_0_xyyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyyz_0, g_x_0_xyyyz_0, g_x_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_0[k] = -g_x_0_xyyyz_0[k] * ab_y + g_x_0_xyyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyzz_0, g_x_0_xyyzz_0, g_x_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_0[k] = -g_x_0_xyyzz_0[k] * ab_y + g_x_0_xyyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzzz_0, g_x_0_xyzzz_0, g_x_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_0[k] = -g_x_0_xyzzz_0[k] * ab_y + g_x_0_xyzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzzz_0, g_x_0_xzzzz_0, g_x_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_0[k] = -g_x_0_xzzzz_0[k] * ab_y + g_x_0_xzzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzzz_0, g_x_0_xzzzz_z, g_x_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_0[k] = -g_x_0_xzzzz_0[k] * ab_z + g_x_0_xzzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_0, g_x_0_yyyyy_y, g_x_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_0[k] = -g_x_0_yyyyy_0[k] * ab_y + g_x_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyyz_0, g_x_0_yyyyz_0, g_x_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_0[k] = -g_x_0_yyyyz_0[k] * ab_y + g_x_0_yyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyzz_0, g_x_0_yyyzz_0, g_x_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_0[k] = -g_x_0_yyyzz_0[k] * ab_y + g_x_0_yyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzzz_0, g_x_0_yyzzz_0, g_x_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_0[k] = -g_x_0_yyzzz_0[k] * ab_y + g_x_0_yyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzzz_0, g_x_0_yzzzz_0, g_x_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_0[k] = -g_x_0_yzzzz_0[k] * ab_y + g_x_0_yzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzzz_0, g_x_0_zzzzz_0, g_x_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_0[k] = -g_x_0_zzzzz_0[k] * ab_y + g_x_0_zzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_0, g_x_0_zzzzz_z, g_x_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_0[k] = -g_x_0_zzzzz_0[k] * ab_z + g_x_0_zzzzz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxx_0, g_y_0_xxxxx_x, g_y_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_0[k] = -g_y_0_xxxxx_0[k] * ab_x + g_y_0_xxxxx_x[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxy_0, g_y_0_xxxxy_0, g_y_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_0[k] = -g_y_0_xxxxy_0[k] * ab_x + g_y_0_xxxxy_x[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxz_0, g_y_0_xxxxz_0, g_y_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_0[k] = -g_y_0_xxxxz_0[k] * ab_x + g_y_0_xxxxz_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyy_0, g_y_0_xxxyy_0, g_y_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_0[k] = -g_y_0_xxxyy_0[k] * ab_x + g_y_0_xxxyy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyz_0, g_y_0_xxxyz_0, g_y_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_0[k] = -g_y_0_xxxyz_0[k] * ab_x + g_y_0_xxxyz_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxzz_0, g_y_0_xxxzz_0, g_y_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_0[k] = -g_y_0_xxxzz_0[k] * ab_x + g_y_0_xxxzz_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyy_0, g_y_0_xxyyy_0, g_y_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_0[k] = -g_y_0_xxyyy_0[k] * ab_x + g_y_0_xxyyy_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyz_0, g_y_0_xxyyz_0, g_y_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_0[k] = -g_y_0_xxyyz_0[k] * ab_x + g_y_0_xxyyz_x[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyzz_0, g_y_0_xxyzz_0, g_y_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_0[k] = -g_y_0_xxyzz_0[k] * ab_x + g_y_0_xxyzz_x[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzzz_0, g_y_0_xxzzz_0, g_y_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_0[k] = -g_y_0_xxzzz_0[k] * ab_x + g_y_0_xxzzz_x[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyy_0, g_y_0_xyyyy_0, g_y_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_0[k] = -g_y_0_xyyyy_0[k] * ab_x + g_y_0_xyyyy_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyz_0, g_y_0_xyyyz_0, g_y_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_0[k] = -g_y_0_xyyyz_0[k] * ab_x + g_y_0_xyyyz_x[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyzz_0, g_y_0_xyyzz_0, g_y_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_0[k] = -g_y_0_xyyzz_0[k] * ab_x + g_y_0_xyyzz_x[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzzz_0, g_y_0_xyzzz_0, g_y_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_0[k] = -g_y_0_xyzzz_0[k] * ab_x + g_y_0_xyzzz_x[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzzz_0, g_y_0_xzzzz_0, g_y_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_0[k] = -g_y_0_xzzzz_0[k] * ab_x + g_y_0_xzzzz_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyy_0, g_y_0_yyyyy_0, g_y_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_0[k] = -g_y_0_yyyyy_0[k] * ab_x + g_y_0_yyyyy_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyz_0, g_y_0_yyyyz_0, g_y_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_0[k] = -g_y_0_yyyyz_0[k] * ab_x + g_y_0_yyyyz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyzz_0, g_y_0_yyyzz_0, g_y_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_0[k] = -g_y_0_yyyzz_0[k] * ab_x + g_y_0_yyyzz_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzzz_0, g_y_0_yyzzz_0, g_y_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_0[k] = -g_y_0_yyzzz_0[k] * ab_x + g_y_0_yyzzz_x[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzzz_0, g_y_0_yzzzz_0, g_y_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_0[k] = -g_y_0_yzzzz_0[k] * ab_x + g_y_0_yzzzz_x[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzzz_0, g_y_0_zzzzz_0, g_y_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_0[k] = -g_y_0_zzzzz_0[k] * ab_x + g_y_0_zzzzz_x[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_0, g_y_0_yyyyy_y, g_y_0_yyyyyy_0, g_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_0[k] = -g_yyyyy_0[k] - g_y_0_yyyyy_0[k] * ab_y + g_y_0_yyyyy_y[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_0, g_y_0_yyyyy_z, g_y_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_0[k] = -g_y_0_yyyyy_0[k] * ab_z + g_y_0_yyyyy_z[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyz_0, g_y_0_yyyyz_z, g_y_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_0[k] = -g_y_0_yyyyz_0[k] * ab_z + g_y_0_yyyyz_z[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyzz_0, g_y_0_yyyzz_z, g_y_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_0[k] = -g_y_0_yyyzz_0[k] * ab_z + g_y_0_yyyzz_z[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzzz_0, g_y_0_yyzzz_z, g_y_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_0[k] = -g_y_0_yyzzz_0[k] * ab_z + g_y_0_yyzzz_z[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzzz_0, g_y_0_yzzzz_z, g_y_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_0[k] = -g_y_0_yzzzz_0[k] * ab_z + g_y_0_yzzzz_z[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_0, g_y_0_zzzzz_z, g_y_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_0[k] = -g_y_0_zzzzz_0[k] * ab_z + g_y_0_zzzzz_z[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxx_0, g_z_0_xxxxx_x, g_z_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_0[k] = -g_z_0_xxxxx_0[k] * ab_x + g_z_0_xxxxx_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxy_0, g_z_0_xxxxy_0, g_z_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_0[k] = -g_z_0_xxxxy_0[k] * ab_x + g_z_0_xxxxy_x[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxz_0, g_z_0_xxxxz_0, g_z_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_0[k] = -g_z_0_xxxxz_0[k] * ab_x + g_z_0_xxxxz_x[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyy_0, g_z_0_xxxyy_0, g_z_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_0[k] = -g_z_0_xxxyy_0[k] * ab_x + g_z_0_xxxyy_x[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyz_0, g_z_0_xxxyz_0, g_z_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_0[k] = -g_z_0_xxxyz_0[k] * ab_x + g_z_0_xxxyz_x[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxzz_0, g_z_0_xxxzz_0, g_z_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_0[k] = -g_z_0_xxxzz_0[k] * ab_x + g_z_0_xxxzz_x[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyy_0, g_z_0_xxyyy_0, g_z_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_0[k] = -g_z_0_xxyyy_0[k] * ab_x + g_z_0_xxyyy_x[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyz_0, g_z_0_xxyyz_0, g_z_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_0[k] = -g_z_0_xxyyz_0[k] * ab_x + g_z_0_xxyyz_x[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyzz_0, g_z_0_xxyzz_0, g_z_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_0[k] = -g_z_0_xxyzz_0[k] * ab_x + g_z_0_xxyzz_x[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzzz_0, g_z_0_xxzzz_0, g_z_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_0[k] = -g_z_0_xxzzz_0[k] * ab_x + g_z_0_xxzzz_x[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyy_0, g_z_0_xyyyy_0, g_z_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_0[k] = -g_z_0_xyyyy_0[k] * ab_x + g_z_0_xyyyy_x[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyz_0, g_z_0_xyyyz_0, g_z_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_0[k] = -g_z_0_xyyyz_0[k] * ab_x + g_z_0_xyyyz_x[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyzz_0, g_z_0_xyyzz_0, g_z_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_0[k] = -g_z_0_xyyzz_0[k] * ab_x + g_z_0_xyyzz_x[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzzz_0, g_z_0_xyzzz_0, g_z_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_0[k] = -g_z_0_xyzzz_0[k] * ab_x + g_z_0_xyzzz_x[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzzz_0, g_z_0_xzzzz_0, g_z_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_0[k] = -g_z_0_xzzzz_0[k] * ab_x + g_z_0_xzzzz_x[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyy_0, g_z_0_yyyyy_0, g_z_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_0[k] = -g_z_0_yyyyy_0[k] * ab_x + g_z_0_yyyyy_x[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyz_0, g_z_0_yyyyz_0, g_z_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_0[k] = -g_z_0_yyyyz_0[k] * ab_x + g_z_0_yyyyz_x[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyzz_0, g_z_0_yyyzz_0, g_z_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_0[k] = -g_z_0_yyyzz_0[k] * ab_x + g_z_0_yyyzz_x[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzzz_0, g_z_0_yyzzz_0, g_z_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_0[k] = -g_z_0_yyzzz_0[k] * ab_x + g_z_0_yyzzz_x[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzzz_0, g_z_0_yzzzz_0, g_z_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_0[k] = -g_z_0_yzzzz_0[k] * ab_x + g_z_0_yzzzz_x[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzzz_0, g_z_0_zzzzz_0, g_z_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_0[k] = -g_z_0_zzzzz_0[k] * ab_x + g_z_0_zzzzz_x[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyy_0, g_z_0_yyyyy_y, g_z_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_0[k] = -g_z_0_yyyyy_0[k] * ab_y + g_z_0_yyyyy_y[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyyz_0, g_z_0_yyyyz_0, g_z_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_0[k] = -g_z_0_yyyyz_0[k] * ab_y + g_z_0_yyyyz_y[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyzz_0, g_z_0_yyyzz_0, g_z_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_0[k] = -g_z_0_yyyzz_0[k] * ab_y + g_z_0_yyyzz_y[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzzz_0, g_z_0_yyzzz_0, g_z_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_0[k] = -g_z_0_yyzzz_0[k] * ab_y + g_z_0_yyzzz_y[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzzz_0, g_z_0_yzzzz_0, g_z_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_0[k] = -g_z_0_yzzzz_0[k] * ab_y + g_z_0_yzzzz_y[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzzz_0, g_z_0_zzzzz_0, g_z_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_0[k] = -g_z_0_zzzzz_0[k] * ab_y + g_z_0_zzzzz_y[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_0, g_z_0_zzzzz_z, g_z_0_zzzzzz_0, g_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_0[k] = -g_zzzzz_0[k] - g_z_0_zzzzz_0[k] * ab_z + g_z_0_zzzzz_z[k];
            }
        }
    }
}

} // erirec namespace

