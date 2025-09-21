#include "ElectronRepulsionGeom1010ContrRecISXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_isxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_isxx,
                                              const size_t idx_geom_0010_hsxx,
                                              const size_t idx_geom_1010_hsxx,
                                              const size_t idx_geom_1010_hpxx,
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

            const auto hs_geom_0010_off = idx_geom_0010_hsxx + i * dcomps + j;

            auto g_0_0_x_0_xxxxx_0 = cbuffer.data(hs_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_0 = cbuffer.data(hs_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_0 = cbuffer.data(hs_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_0 = cbuffer.data(hs_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_0 = cbuffer.data(hs_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_0 = cbuffer.data(hs_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_0 = cbuffer.data(hs_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_0 = cbuffer.data(hs_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_0 = cbuffer.data(hs_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_0 = cbuffer.data(hs_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_0 = cbuffer.data(hs_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_0 = cbuffer.data(hs_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_0 = cbuffer.data(hs_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_0 = cbuffer.data(hs_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_0 = cbuffer.data(hs_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_0 = cbuffer.data(hs_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_0 = cbuffer.data(hs_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_0 = cbuffer.data(hs_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_0 = cbuffer.data(hs_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_0 = cbuffer.data(hs_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_0 = cbuffer.data(hs_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_0 = cbuffer.data(hs_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_0 = cbuffer.data(hs_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_0 = cbuffer.data(hs_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_0 = cbuffer.data(hs_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_0 = cbuffer.data(hs_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_0 = cbuffer.data(hs_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_0 = cbuffer.data(hs_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_0 = cbuffer.data(hs_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_0 = cbuffer.data(hs_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_0 = cbuffer.data(hs_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_0 = cbuffer.data(hs_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_0 = cbuffer.data(hs_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_0 = cbuffer.data(hs_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_0 = cbuffer.data(hs_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_0 = cbuffer.data(hs_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_0 = cbuffer.data(hs_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_0 = cbuffer.data(hs_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_0 = cbuffer.data(hs_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_0 = cbuffer.data(hs_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_0 = cbuffer.data(hs_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_0 = cbuffer.data(hs_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_0 = cbuffer.data(hs_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_0 = cbuffer.data(hs_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_0 = cbuffer.data(hs_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_0 = cbuffer.data(hs_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_0 = cbuffer.data(hs_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_0 = cbuffer.data(hs_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_0 = cbuffer.data(hs_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_0 = cbuffer.data(hs_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_0 = cbuffer.data(hs_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_0 = cbuffer.data(hs_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_0 = cbuffer.data(hs_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_0 = cbuffer.data(hs_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_0 = cbuffer.data(hs_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_0 = cbuffer.data(hs_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_0 = cbuffer.data(hs_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_0 = cbuffer.data(hs_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_0 = cbuffer.data(hs_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_0 = cbuffer.data(hs_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_0 = cbuffer.data(hs_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_0 = cbuffer.data(hs_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_0 = cbuffer.data(hs_geom_0010_off + 62 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HSSS

            const auto hs_geom_1010_off = idx_geom_1010_hsxx + i * dcomps + j;

            auto g_x_0_x_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 62 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 63 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 64 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 65 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 66 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 67 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 68 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 69 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 70 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 71 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 72 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 73 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 74 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 75 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 76 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 77 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 78 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 79 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 80 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 81 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 82 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 83 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 84 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 85 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 86 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 87 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 88 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 89 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 90 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 91 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 92 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 93 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 94 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 95 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 96 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 97 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 98 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 99 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 100 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 101 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 102 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 103 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 104 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 105 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 106 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 107 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 108 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 109 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 110 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 111 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 112 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 113 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 114 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 115 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 116 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 117 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 118 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 119 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 120 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 121 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 122 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 123 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 124 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 125 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 126 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 127 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 128 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 129 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 130 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 131 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 132 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 133 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 134 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 135 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 136 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 137 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 138 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 139 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 140 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 141 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 142 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 143 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 144 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 145 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 146 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 147 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 148 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 149 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 150 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 151 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 152 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 153 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 154 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 155 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 156 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 157 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 158 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 159 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 160 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 161 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 162 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 163 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 164 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 165 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 166 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 167 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_0 = cbuffer.data(hs_geom_1010_off + 168 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_0 = cbuffer.data(hs_geom_1010_off + 169 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_0 = cbuffer.data(hs_geom_1010_off + 170 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_0 = cbuffer.data(hs_geom_1010_off + 171 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_0 = cbuffer.data(hs_geom_1010_off + 172 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_0 = cbuffer.data(hs_geom_1010_off + 173 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_0 = cbuffer.data(hs_geom_1010_off + 174 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_0 = cbuffer.data(hs_geom_1010_off + 175 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_0 = cbuffer.data(hs_geom_1010_off + 176 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_0 = cbuffer.data(hs_geom_1010_off + 177 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_0 = cbuffer.data(hs_geom_1010_off + 178 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_0 = cbuffer.data(hs_geom_1010_off + 179 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_0 = cbuffer.data(hs_geom_1010_off + 180 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_0 = cbuffer.data(hs_geom_1010_off + 181 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_0 = cbuffer.data(hs_geom_1010_off + 182 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_0 = cbuffer.data(hs_geom_1010_off + 183 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_0 = cbuffer.data(hs_geom_1010_off + 184 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_0 = cbuffer.data(hs_geom_1010_off + 185 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_0 = cbuffer.data(hs_geom_1010_off + 186 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_0 = cbuffer.data(hs_geom_1010_off + 187 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_0 = cbuffer.data(hs_geom_1010_off + 188 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_1010_off = idx_geom_1010_hpxx + i * dcomps + j;

            auto g_x_0_x_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 188 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 189 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 190 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 191 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 192 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 193 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 194 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 195 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 196 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 197 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 198 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 199 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 200 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 201 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 202 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 203 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 204 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 205 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 206 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 207 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 208 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 209 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 210 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 211 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 212 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 213 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 214 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 215 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 216 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 217 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 218 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 219 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 220 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 221 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 222 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 223 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 224 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 225 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 226 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 227 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 228 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 229 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 230 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 231 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 232 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 233 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 234 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 235 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 236 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 237 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 238 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 239 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 240 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 241 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 242 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 243 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 244 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 245 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 246 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 247 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 248 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 249 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 250 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 251 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 254 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 257 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 260 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 263 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 266 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 269 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 270 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 271 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 272 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 273 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 274 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 275 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 276 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 277 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 278 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 279 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 280 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 281 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 282 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 283 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 284 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 285 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 286 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 287 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 288 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 289 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 290 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 291 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 292 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 293 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 294 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 295 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 296 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 297 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 298 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 299 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 300 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 301 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 302 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 303 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 304 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 305 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 306 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 307 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 308 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 309 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 310 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 311 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 312 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 313 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 314 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 315 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 316 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 317 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 318 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 319 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 320 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 321 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 322 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 323 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 324 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 325 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 326 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 327 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 328 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 329 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 330 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 331 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 332 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 333 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 334 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 335 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 336 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 337 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 338 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 339 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 340 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 341 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 342 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 343 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 344 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 345 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 346 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 347 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 348 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 349 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 350 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 351 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 352 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 353 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 354 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 355 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 356 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 357 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 358 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 359 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 360 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 361 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 362 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 363 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 364 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 365 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 366 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 367 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 368 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 369 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 370 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 371 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 372 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 373 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 374 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 375 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 376 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 377 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 378 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 379 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 380 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 381 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 382 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 383 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 384 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 385 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 386 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 387 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 388 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 389 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 390 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 391 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 392 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 393 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 394 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 395 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 396 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 397 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 398 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 399 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 400 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 401 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 402 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 403 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 404 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 405 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 406 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 407 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 408 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 409 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 410 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 411 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 412 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 413 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 414 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 415 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 416 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 417 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 418 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 419 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 420 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 421 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 422 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 423 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 424 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 425 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 426 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 427 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 428 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 429 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 430 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 431 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 432 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 433 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 434 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 435 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 436 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 437 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 438 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 439 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 440 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 441 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 442 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 443 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 444 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 445 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 446 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 447 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 448 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 449 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 450 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 451 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 452 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 453 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 454 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 455 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 456 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 457 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 458 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 459 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 460 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 461 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 462 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 463 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 464 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 465 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 466 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 467 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 468 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 469 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 470 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 471 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 472 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 473 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 474 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 475 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 476 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 477 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 478 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 479 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 480 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 481 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 482 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 483 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 484 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 485 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 486 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 487 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 488 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 489 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 490 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 491 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 492 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 493 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 494 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 495 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 496 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 497 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 498 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 499 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 500 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 501 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 502 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 503 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_x = cbuffer.data(hp_geom_1010_off + 504 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_y = cbuffer.data(hp_geom_1010_off + 505 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_z = cbuffer.data(hp_geom_1010_off + 506 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_x = cbuffer.data(hp_geom_1010_off + 507 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_y = cbuffer.data(hp_geom_1010_off + 508 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_z = cbuffer.data(hp_geom_1010_off + 509 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_x = cbuffer.data(hp_geom_1010_off + 510 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_y = cbuffer.data(hp_geom_1010_off + 511 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_z = cbuffer.data(hp_geom_1010_off + 512 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_x = cbuffer.data(hp_geom_1010_off + 513 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_y = cbuffer.data(hp_geom_1010_off + 514 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_z = cbuffer.data(hp_geom_1010_off + 515 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_x = cbuffer.data(hp_geom_1010_off + 516 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_y = cbuffer.data(hp_geom_1010_off + 517 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_z = cbuffer.data(hp_geom_1010_off + 518 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_x = cbuffer.data(hp_geom_1010_off + 519 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_y = cbuffer.data(hp_geom_1010_off + 520 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_z = cbuffer.data(hp_geom_1010_off + 521 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_x = cbuffer.data(hp_geom_1010_off + 522 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_y = cbuffer.data(hp_geom_1010_off + 523 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_z = cbuffer.data(hp_geom_1010_off + 524 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_x = cbuffer.data(hp_geom_1010_off + 525 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_y = cbuffer.data(hp_geom_1010_off + 526 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_z = cbuffer.data(hp_geom_1010_off + 527 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_x = cbuffer.data(hp_geom_1010_off + 528 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_y = cbuffer.data(hp_geom_1010_off + 529 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_z = cbuffer.data(hp_geom_1010_off + 530 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_x = cbuffer.data(hp_geom_1010_off + 531 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_y = cbuffer.data(hp_geom_1010_off + 532 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_z = cbuffer.data(hp_geom_1010_off + 533 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_x = cbuffer.data(hp_geom_1010_off + 534 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_y = cbuffer.data(hp_geom_1010_off + 535 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_z = cbuffer.data(hp_geom_1010_off + 536 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_x = cbuffer.data(hp_geom_1010_off + 537 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_y = cbuffer.data(hp_geom_1010_off + 538 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_z = cbuffer.data(hp_geom_1010_off + 539 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_x = cbuffer.data(hp_geom_1010_off + 540 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_y = cbuffer.data(hp_geom_1010_off + 541 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_z = cbuffer.data(hp_geom_1010_off + 542 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_x = cbuffer.data(hp_geom_1010_off + 543 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_y = cbuffer.data(hp_geom_1010_off + 544 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_z = cbuffer.data(hp_geom_1010_off + 545 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_x = cbuffer.data(hp_geom_1010_off + 546 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_y = cbuffer.data(hp_geom_1010_off + 547 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_z = cbuffer.data(hp_geom_1010_off + 548 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_x = cbuffer.data(hp_geom_1010_off + 549 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_y = cbuffer.data(hp_geom_1010_off + 550 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_z = cbuffer.data(hp_geom_1010_off + 551 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_x = cbuffer.data(hp_geom_1010_off + 552 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_y = cbuffer.data(hp_geom_1010_off + 553 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_z = cbuffer.data(hp_geom_1010_off + 554 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_x = cbuffer.data(hp_geom_1010_off + 555 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_y = cbuffer.data(hp_geom_1010_off + 556 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_z = cbuffer.data(hp_geom_1010_off + 557 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_x = cbuffer.data(hp_geom_1010_off + 558 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_y = cbuffer.data(hp_geom_1010_off + 559 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_z = cbuffer.data(hp_geom_1010_off + 560 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_x = cbuffer.data(hp_geom_1010_off + 561 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_y = cbuffer.data(hp_geom_1010_off + 562 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_z = cbuffer.data(hp_geom_1010_off + 563 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_x = cbuffer.data(hp_geom_1010_off + 564 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_y = cbuffer.data(hp_geom_1010_off + 565 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_z = cbuffer.data(hp_geom_1010_off + 566 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_isxx

            const auto is_geom_1010_off = idx_geom_1010_isxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxxxx_0, g_x_0_x_0_xxxxx_0, g_x_0_x_0_xxxxx_x, g_x_0_x_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxx_0[k] = -g_0_0_x_0_xxxxx_0[k] - g_x_0_x_0_xxxxx_0[k] * ab_x + g_x_0_x_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_0, g_x_0_x_0_xxxxx_y, g_x_0_x_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxy_0[k] = -g_x_0_x_0_xxxxx_0[k] * ab_y + g_x_0_x_0_xxxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_0, g_x_0_x_0_xxxxx_z, g_x_0_x_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxz_0[k] = -g_x_0_x_0_xxxxx_0[k] * ab_z + g_x_0_x_0_xxxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxy_0, g_x_0_x_0_xxxxy_y, g_x_0_x_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyy_0[k] = -g_x_0_x_0_xxxxy_0[k] * ab_y + g_x_0_x_0_xxxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxyz_0, g_x_0_x_0_xxxxz_0, g_x_0_x_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyz_0[k] = -g_x_0_x_0_xxxxz_0[k] * ab_y + g_x_0_x_0_xxxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxz_0, g_x_0_x_0_xxxxz_z, g_x_0_x_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxzz_0[k] = -g_x_0_x_0_xxxxz_0[k] * ab_z + g_x_0_x_0_xxxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyy_0, g_x_0_x_0_xxxyy_y, g_x_0_x_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyy_0[k] = -g_x_0_x_0_xxxyy_0[k] * ab_y + g_x_0_x_0_xxxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyyz_0, g_x_0_x_0_xxxyz_0, g_x_0_x_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyz_0[k] = -g_x_0_x_0_xxxyz_0[k] * ab_y + g_x_0_x_0_xxxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyzz_0, g_x_0_x_0_xxxzz_0, g_x_0_x_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyzz_0[k] = -g_x_0_x_0_xxxzz_0[k] * ab_y + g_x_0_x_0_xxxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxzz_0, g_x_0_x_0_xxxzz_z, g_x_0_x_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxzzz_0[k] = -g_x_0_x_0_xxxzz_0[k] * ab_z + g_x_0_x_0_xxxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyy_0, g_x_0_x_0_xxyyy_y, g_x_0_x_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyy_0[k] = -g_x_0_x_0_xxyyy_0[k] * ab_y + g_x_0_x_0_xxyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyyz_0, g_x_0_x_0_xxyyz_0, g_x_0_x_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyz_0[k] = -g_x_0_x_0_xxyyz_0[k] * ab_y + g_x_0_x_0_xxyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyzz_0, g_x_0_x_0_xxyzz_0, g_x_0_x_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyzz_0[k] = -g_x_0_x_0_xxyzz_0[k] * ab_y + g_x_0_x_0_xxyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyzzz_0, g_x_0_x_0_xxzzz_0, g_x_0_x_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyzzz_0[k] = -g_x_0_x_0_xxzzz_0[k] * ab_y + g_x_0_x_0_xxzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxzzz_0, g_x_0_x_0_xxzzz_z, g_x_0_x_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzzzz_0[k] = -g_x_0_x_0_xxzzz_0[k] * ab_z + g_x_0_x_0_xxzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyy_0, g_x_0_x_0_xyyyy_y, g_x_0_x_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyy_0[k] = -g_x_0_x_0_xyyyy_0[k] * ab_y + g_x_0_x_0_xyyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyyz_0, g_x_0_x_0_xyyyz_0, g_x_0_x_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyz_0[k] = -g_x_0_x_0_xyyyz_0[k] * ab_y + g_x_0_x_0_xyyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyzz_0, g_x_0_x_0_xyyzz_0, g_x_0_x_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyzz_0[k] = -g_x_0_x_0_xyyzz_0[k] * ab_y + g_x_0_x_0_xyyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyzzz_0, g_x_0_x_0_xyzzz_0, g_x_0_x_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyzzz_0[k] = -g_x_0_x_0_xyzzz_0[k] * ab_y + g_x_0_x_0_xyzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzzzz_0, g_x_0_x_0_xzzzz_0, g_x_0_x_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzzzz_0[k] = -g_x_0_x_0_xzzzz_0[k] * ab_y + g_x_0_x_0_xzzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzzzz_0, g_x_0_x_0_xzzzz_z, g_x_0_x_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzzzz_0[k] = -g_x_0_x_0_xzzzz_0[k] * ab_z + g_x_0_x_0_xzzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyy_0, g_x_0_x_0_yyyyy_y, g_x_0_x_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyy_0[k] = -g_x_0_x_0_yyyyy_0[k] * ab_y + g_x_0_x_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyyz_0, g_x_0_x_0_yyyyz_0, g_x_0_x_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyz_0[k] = -g_x_0_x_0_yyyyz_0[k] * ab_y + g_x_0_x_0_yyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyzz_0, g_x_0_x_0_yyyzz_0, g_x_0_x_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyzz_0[k] = -g_x_0_x_0_yyyzz_0[k] * ab_y + g_x_0_x_0_yyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyzzz_0, g_x_0_x_0_yyzzz_0, g_x_0_x_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyzzz_0[k] = -g_x_0_x_0_yyzzz_0[k] * ab_y + g_x_0_x_0_yyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzzzz_0, g_x_0_x_0_yzzzz_0, g_x_0_x_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzzzz_0[k] = -g_x_0_x_0_yzzzz_0[k] * ab_y + g_x_0_x_0_yzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzzzz_0, g_x_0_x_0_zzzzz_0, g_x_0_x_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzzzz_0[k] = -g_x_0_x_0_zzzzz_0[k] * ab_y + g_x_0_x_0_zzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzzzz_0, g_x_0_x_0_zzzzz_z, g_x_0_x_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzzzz_0[k] = -g_x_0_x_0_zzzzz_0[k] * ab_z + g_x_0_x_0_zzzzz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxxxx_0, g_x_0_y_0_xxxxx_0, g_x_0_y_0_xxxxx_x, g_x_0_y_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxx_0[k] = -g_0_0_y_0_xxxxx_0[k] - g_x_0_y_0_xxxxx_0[k] * ab_x + g_x_0_y_0_xxxxx_x[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_0, g_x_0_y_0_xxxxx_y, g_x_0_y_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxy_0[k] = -g_x_0_y_0_xxxxx_0[k] * ab_y + g_x_0_y_0_xxxxx_y[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_0, g_x_0_y_0_xxxxx_z, g_x_0_y_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxz_0[k] = -g_x_0_y_0_xxxxx_0[k] * ab_z + g_x_0_y_0_xxxxx_z[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxy_0, g_x_0_y_0_xxxxy_y, g_x_0_y_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyy_0[k] = -g_x_0_y_0_xxxxy_0[k] * ab_y + g_x_0_y_0_xxxxy_y[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxyz_0, g_x_0_y_0_xxxxz_0, g_x_0_y_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyz_0[k] = -g_x_0_y_0_xxxxz_0[k] * ab_y + g_x_0_y_0_xxxxz_y[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxz_0, g_x_0_y_0_xxxxz_z, g_x_0_y_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxzz_0[k] = -g_x_0_y_0_xxxxz_0[k] * ab_z + g_x_0_y_0_xxxxz_z[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyy_0, g_x_0_y_0_xxxyy_y, g_x_0_y_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyy_0[k] = -g_x_0_y_0_xxxyy_0[k] * ab_y + g_x_0_y_0_xxxyy_y[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyyz_0, g_x_0_y_0_xxxyz_0, g_x_0_y_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyz_0[k] = -g_x_0_y_0_xxxyz_0[k] * ab_y + g_x_0_y_0_xxxyz_y[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyzz_0, g_x_0_y_0_xxxzz_0, g_x_0_y_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyzz_0[k] = -g_x_0_y_0_xxxzz_0[k] * ab_y + g_x_0_y_0_xxxzz_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxzz_0, g_x_0_y_0_xxxzz_z, g_x_0_y_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxzzz_0[k] = -g_x_0_y_0_xxxzz_0[k] * ab_z + g_x_0_y_0_xxxzz_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyy_0, g_x_0_y_0_xxyyy_y, g_x_0_y_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyy_0[k] = -g_x_0_y_0_xxyyy_0[k] * ab_y + g_x_0_y_0_xxyyy_y[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyyz_0, g_x_0_y_0_xxyyz_0, g_x_0_y_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyz_0[k] = -g_x_0_y_0_xxyyz_0[k] * ab_y + g_x_0_y_0_xxyyz_y[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyzz_0, g_x_0_y_0_xxyzz_0, g_x_0_y_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyzz_0[k] = -g_x_0_y_0_xxyzz_0[k] * ab_y + g_x_0_y_0_xxyzz_y[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyzzz_0, g_x_0_y_0_xxzzz_0, g_x_0_y_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyzzz_0[k] = -g_x_0_y_0_xxzzz_0[k] * ab_y + g_x_0_y_0_xxzzz_y[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxzzz_0, g_x_0_y_0_xxzzz_z, g_x_0_y_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzzzz_0[k] = -g_x_0_y_0_xxzzz_0[k] * ab_z + g_x_0_y_0_xxzzz_z[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyy_0, g_x_0_y_0_xyyyy_y, g_x_0_y_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyy_0[k] = -g_x_0_y_0_xyyyy_0[k] * ab_y + g_x_0_y_0_xyyyy_y[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyyz_0, g_x_0_y_0_xyyyz_0, g_x_0_y_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyz_0[k] = -g_x_0_y_0_xyyyz_0[k] * ab_y + g_x_0_y_0_xyyyz_y[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyzz_0, g_x_0_y_0_xyyzz_0, g_x_0_y_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyzz_0[k] = -g_x_0_y_0_xyyzz_0[k] * ab_y + g_x_0_y_0_xyyzz_y[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyzzz_0, g_x_0_y_0_xyzzz_0, g_x_0_y_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyzzz_0[k] = -g_x_0_y_0_xyzzz_0[k] * ab_y + g_x_0_y_0_xyzzz_y[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzzzz_0, g_x_0_y_0_xzzzz_0, g_x_0_y_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzzzz_0[k] = -g_x_0_y_0_xzzzz_0[k] * ab_y + g_x_0_y_0_xzzzz_y[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzzzz_0, g_x_0_y_0_xzzzz_z, g_x_0_y_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzzzz_0[k] = -g_x_0_y_0_xzzzz_0[k] * ab_z + g_x_0_y_0_xzzzz_z[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyy_0, g_x_0_y_0_yyyyy_y, g_x_0_y_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyy_0[k] = -g_x_0_y_0_yyyyy_0[k] * ab_y + g_x_0_y_0_yyyyy_y[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyyz_0, g_x_0_y_0_yyyyz_0, g_x_0_y_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyz_0[k] = -g_x_0_y_0_yyyyz_0[k] * ab_y + g_x_0_y_0_yyyyz_y[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyzz_0, g_x_0_y_0_yyyzz_0, g_x_0_y_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyzz_0[k] = -g_x_0_y_0_yyyzz_0[k] * ab_y + g_x_0_y_0_yyyzz_y[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyzzz_0, g_x_0_y_0_yyzzz_0, g_x_0_y_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyzzz_0[k] = -g_x_0_y_0_yyzzz_0[k] * ab_y + g_x_0_y_0_yyzzz_y[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzzzz_0, g_x_0_y_0_yzzzz_0, g_x_0_y_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzzzz_0[k] = -g_x_0_y_0_yzzzz_0[k] * ab_y + g_x_0_y_0_yzzzz_y[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzzzz_0, g_x_0_y_0_zzzzz_0, g_x_0_y_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzzzz_0[k] = -g_x_0_y_0_zzzzz_0[k] * ab_y + g_x_0_y_0_zzzzz_y[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzzzz_0, g_x_0_y_0_zzzzz_z, g_x_0_y_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzzzz_0[k] = -g_x_0_y_0_zzzzz_0[k] * ab_z + g_x_0_y_0_zzzzz_z[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxxxx_0, g_x_0_z_0_xxxxx_0, g_x_0_z_0_xxxxx_x, g_x_0_z_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxx_0[k] = -g_0_0_z_0_xxxxx_0[k] - g_x_0_z_0_xxxxx_0[k] * ab_x + g_x_0_z_0_xxxxx_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_0, g_x_0_z_0_xxxxx_y, g_x_0_z_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxy_0[k] = -g_x_0_z_0_xxxxx_0[k] * ab_y + g_x_0_z_0_xxxxx_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_0, g_x_0_z_0_xxxxx_z, g_x_0_z_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxz_0[k] = -g_x_0_z_0_xxxxx_0[k] * ab_z + g_x_0_z_0_xxxxx_z[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxy_0, g_x_0_z_0_xxxxy_y, g_x_0_z_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyy_0[k] = -g_x_0_z_0_xxxxy_0[k] * ab_y + g_x_0_z_0_xxxxy_y[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxyz_0, g_x_0_z_0_xxxxz_0, g_x_0_z_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyz_0[k] = -g_x_0_z_0_xxxxz_0[k] * ab_y + g_x_0_z_0_xxxxz_y[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxz_0, g_x_0_z_0_xxxxz_z, g_x_0_z_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxzz_0[k] = -g_x_0_z_0_xxxxz_0[k] * ab_z + g_x_0_z_0_xxxxz_z[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyy_0, g_x_0_z_0_xxxyy_y, g_x_0_z_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyy_0[k] = -g_x_0_z_0_xxxyy_0[k] * ab_y + g_x_0_z_0_xxxyy_y[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyyz_0, g_x_0_z_0_xxxyz_0, g_x_0_z_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyz_0[k] = -g_x_0_z_0_xxxyz_0[k] * ab_y + g_x_0_z_0_xxxyz_y[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyzz_0, g_x_0_z_0_xxxzz_0, g_x_0_z_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyzz_0[k] = -g_x_0_z_0_xxxzz_0[k] * ab_y + g_x_0_z_0_xxxzz_y[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxzz_0, g_x_0_z_0_xxxzz_z, g_x_0_z_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxzzz_0[k] = -g_x_0_z_0_xxxzz_0[k] * ab_z + g_x_0_z_0_xxxzz_z[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyy_0, g_x_0_z_0_xxyyy_y, g_x_0_z_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyy_0[k] = -g_x_0_z_0_xxyyy_0[k] * ab_y + g_x_0_z_0_xxyyy_y[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyyz_0, g_x_0_z_0_xxyyz_0, g_x_0_z_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyz_0[k] = -g_x_0_z_0_xxyyz_0[k] * ab_y + g_x_0_z_0_xxyyz_y[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyzz_0, g_x_0_z_0_xxyzz_0, g_x_0_z_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyzz_0[k] = -g_x_0_z_0_xxyzz_0[k] * ab_y + g_x_0_z_0_xxyzz_y[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyzzz_0, g_x_0_z_0_xxzzz_0, g_x_0_z_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyzzz_0[k] = -g_x_0_z_0_xxzzz_0[k] * ab_y + g_x_0_z_0_xxzzz_y[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxzzz_0, g_x_0_z_0_xxzzz_z, g_x_0_z_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzzzz_0[k] = -g_x_0_z_0_xxzzz_0[k] * ab_z + g_x_0_z_0_xxzzz_z[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyy_0, g_x_0_z_0_xyyyy_y, g_x_0_z_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyy_0[k] = -g_x_0_z_0_xyyyy_0[k] * ab_y + g_x_0_z_0_xyyyy_y[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyyz_0, g_x_0_z_0_xyyyz_0, g_x_0_z_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyz_0[k] = -g_x_0_z_0_xyyyz_0[k] * ab_y + g_x_0_z_0_xyyyz_y[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyzz_0, g_x_0_z_0_xyyzz_0, g_x_0_z_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyzz_0[k] = -g_x_0_z_0_xyyzz_0[k] * ab_y + g_x_0_z_0_xyyzz_y[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyzzz_0, g_x_0_z_0_xyzzz_0, g_x_0_z_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyzzz_0[k] = -g_x_0_z_0_xyzzz_0[k] * ab_y + g_x_0_z_0_xyzzz_y[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzzzz_0, g_x_0_z_0_xzzzz_0, g_x_0_z_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzzzz_0[k] = -g_x_0_z_0_xzzzz_0[k] * ab_y + g_x_0_z_0_xzzzz_y[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzzzz_0, g_x_0_z_0_xzzzz_z, g_x_0_z_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzzzz_0[k] = -g_x_0_z_0_xzzzz_0[k] * ab_z + g_x_0_z_0_xzzzz_z[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyy_0, g_x_0_z_0_yyyyy_y, g_x_0_z_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyy_0[k] = -g_x_0_z_0_yyyyy_0[k] * ab_y + g_x_0_z_0_yyyyy_y[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyyz_0, g_x_0_z_0_yyyyz_0, g_x_0_z_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyz_0[k] = -g_x_0_z_0_yyyyz_0[k] * ab_y + g_x_0_z_0_yyyyz_y[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyzz_0, g_x_0_z_0_yyyzz_0, g_x_0_z_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyzz_0[k] = -g_x_0_z_0_yyyzz_0[k] * ab_y + g_x_0_z_0_yyyzz_y[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyzzz_0, g_x_0_z_0_yyzzz_0, g_x_0_z_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyzzz_0[k] = -g_x_0_z_0_yyzzz_0[k] * ab_y + g_x_0_z_0_yyzzz_y[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzzzz_0, g_x_0_z_0_yzzzz_0, g_x_0_z_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzzzz_0[k] = -g_x_0_z_0_yzzzz_0[k] * ab_y + g_x_0_z_0_yzzzz_y[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzzzz_0, g_x_0_z_0_zzzzz_0, g_x_0_z_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzzzz_0[k] = -g_x_0_z_0_zzzzz_0[k] * ab_y + g_x_0_z_0_zzzzz_y[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzzzz_0, g_x_0_z_0_zzzzz_z, g_x_0_z_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzzzz_0[k] = -g_x_0_z_0_zzzzz_0[k] * ab_z + g_x_0_z_0_zzzzz_z[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxx_0, g_y_0_x_0_xxxxx_x, g_y_0_x_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxx_0[k] = -g_y_0_x_0_xxxxx_0[k] * ab_x + g_y_0_x_0_xxxxx_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxy_0, g_y_0_x_0_xxxxy_0, g_y_0_x_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxy_0[k] = -g_y_0_x_0_xxxxy_0[k] * ab_x + g_y_0_x_0_xxxxy_x[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxz_0, g_y_0_x_0_xxxxz_0, g_y_0_x_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxz_0[k] = -g_y_0_x_0_xxxxz_0[k] * ab_x + g_y_0_x_0_xxxxz_x[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyy_0, g_y_0_x_0_xxxyy_0, g_y_0_x_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyy_0[k] = -g_y_0_x_0_xxxyy_0[k] * ab_x + g_y_0_x_0_xxxyy_x[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyz_0, g_y_0_x_0_xxxyz_0, g_y_0_x_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyz_0[k] = -g_y_0_x_0_xxxyz_0[k] * ab_x + g_y_0_x_0_xxxyz_x[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxzz_0, g_y_0_x_0_xxxzz_0, g_y_0_x_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxzz_0[k] = -g_y_0_x_0_xxxzz_0[k] * ab_x + g_y_0_x_0_xxxzz_x[k];
            }

            /// Set up 90-91 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 90 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyy_0, g_y_0_x_0_xxyyy_0, g_y_0_x_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyy_0[k] = -g_y_0_x_0_xxyyy_0[k] * ab_x + g_y_0_x_0_xxyyy_x[k];
            }

            /// Set up 91-92 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 91 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyz_0, g_y_0_x_0_xxyyz_0, g_y_0_x_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyz_0[k] = -g_y_0_x_0_xxyyz_0[k] * ab_x + g_y_0_x_0_xxyyz_x[k];
            }

            /// Set up 92-93 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyzz_0, g_y_0_x_0_xxyzz_0, g_y_0_x_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyzz_0[k] = -g_y_0_x_0_xxyzz_0[k] * ab_x + g_y_0_x_0_xxyzz_x[k];
            }

            /// Set up 93-94 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 93 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxzzz_0, g_y_0_x_0_xxzzz_0, g_y_0_x_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxzzz_0[k] = -g_y_0_x_0_xxzzz_0[k] * ab_x + g_y_0_x_0_xxzzz_x[k];
            }

            /// Set up 94-95 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 94 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyy_0, g_y_0_x_0_xyyyy_0, g_y_0_x_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyy_0[k] = -g_y_0_x_0_xyyyy_0[k] * ab_x + g_y_0_x_0_xyyyy_x[k];
            }

            /// Set up 95-96 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyz_0, g_y_0_x_0_xyyyz_0, g_y_0_x_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyz_0[k] = -g_y_0_x_0_xyyyz_0[k] * ab_x + g_y_0_x_0_xyyyz_x[k];
            }

            /// Set up 96-97 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 96 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyzz_0, g_y_0_x_0_xyyzz_0, g_y_0_x_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyzz_0[k] = -g_y_0_x_0_xyyzz_0[k] * ab_x + g_y_0_x_0_xyyzz_x[k];
            }

            /// Set up 97-98 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 97 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyzzz_0, g_y_0_x_0_xyzzz_0, g_y_0_x_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyzzz_0[k] = -g_y_0_x_0_xyzzz_0[k] * ab_x + g_y_0_x_0_xyzzz_x[k];
            }

            /// Set up 98-99 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzzzz_0, g_y_0_x_0_xzzzz_0, g_y_0_x_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzzzz_0[k] = -g_y_0_x_0_xzzzz_0[k] * ab_x + g_y_0_x_0_xzzzz_x[k];
            }

            /// Set up 99-100 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyy_0, g_y_0_x_0_yyyyy_0, g_y_0_x_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyy_0[k] = -g_y_0_x_0_yyyyy_0[k] * ab_x + g_y_0_x_0_yyyyy_x[k];
            }

            /// Set up 100-101 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 100 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyz_0, g_y_0_x_0_yyyyz_0, g_y_0_x_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyz_0[k] = -g_y_0_x_0_yyyyz_0[k] * ab_x + g_y_0_x_0_yyyyz_x[k];
            }

            /// Set up 101-102 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyzz_0, g_y_0_x_0_yyyzz_0, g_y_0_x_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyzz_0[k] = -g_y_0_x_0_yyyzz_0[k] * ab_x + g_y_0_x_0_yyyzz_x[k];
            }

            /// Set up 102-103 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 102 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyzzz_0, g_y_0_x_0_yyzzz_0, g_y_0_x_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyzzz_0[k] = -g_y_0_x_0_yyzzz_0[k] * ab_x + g_y_0_x_0_yyzzz_x[k];
            }

            /// Set up 103-104 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 103 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzzzz_0, g_y_0_x_0_yzzzz_0, g_y_0_x_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzzzz_0[k] = -g_y_0_x_0_yzzzz_0[k] * ab_x + g_y_0_x_0_yzzzz_x[k];
            }

            /// Set up 104-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzzzz_0, g_y_0_x_0_zzzzz_0, g_y_0_x_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzzzz_0[k] = -g_y_0_x_0_zzzzz_0[k] * ab_x + g_y_0_x_0_zzzzz_x[k];
            }

            /// Set up 105-106 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 105 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyyyy_0, g_y_0_x_0_yyyyy_0, g_y_0_x_0_yyyyy_y, g_y_0_x_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyy_0[k] = -g_0_0_x_0_yyyyy_0[k] - g_y_0_x_0_yyyyy_0[k] * ab_y + g_y_0_x_0_yyyyy_y[k];
            }

            /// Set up 106-107 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 106 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyy_0, g_y_0_x_0_yyyyy_z, g_y_0_x_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyz_0[k] = -g_y_0_x_0_yyyyy_0[k] * ab_z + g_y_0_x_0_yyyyy_z[k];
            }

            /// Set up 107-108 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyz_0, g_y_0_x_0_yyyyz_z, g_y_0_x_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyzz_0[k] = -g_y_0_x_0_yyyyz_0[k] * ab_z + g_y_0_x_0_yyyyz_z[k];
            }

            /// Set up 108-109 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 108 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyzz_0, g_y_0_x_0_yyyzz_z, g_y_0_x_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyzzz_0[k] = -g_y_0_x_0_yyyzz_0[k] * ab_z + g_y_0_x_0_yyyzz_z[k];
            }

            /// Set up 109-110 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyzzz_0, g_y_0_x_0_yyzzz_z, g_y_0_x_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzzzz_0[k] = -g_y_0_x_0_yyzzz_0[k] * ab_z + g_y_0_x_0_yyzzz_z[k];
            }

            /// Set up 110-111 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzzzz_0, g_y_0_x_0_yzzzz_z, g_y_0_x_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzzzz_0[k] = -g_y_0_x_0_yzzzz_0[k] * ab_z + g_y_0_x_0_yzzzz_z[k];
            }

            /// Set up 111-112 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzzzz_0, g_y_0_x_0_zzzzz_z, g_y_0_x_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzzzz_0[k] = -g_y_0_x_0_zzzzz_0[k] * ab_z + g_y_0_x_0_zzzzz_z[k];
            }

            /// Set up 112-113 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 112 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxx_0, g_y_0_y_0_xxxxx_x, g_y_0_y_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxx_0[k] = -g_y_0_y_0_xxxxx_0[k] * ab_x + g_y_0_y_0_xxxxx_x[k];
            }

            /// Set up 113-114 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxy_0, g_y_0_y_0_xxxxy_0, g_y_0_y_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxy_0[k] = -g_y_0_y_0_xxxxy_0[k] * ab_x + g_y_0_y_0_xxxxy_x[k];
            }

            /// Set up 114-115 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 114 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxz_0, g_y_0_y_0_xxxxz_0, g_y_0_y_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxz_0[k] = -g_y_0_y_0_xxxxz_0[k] * ab_x + g_y_0_y_0_xxxxz_x[k];
            }

            /// Set up 115-116 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 115 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyy_0, g_y_0_y_0_xxxyy_0, g_y_0_y_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyy_0[k] = -g_y_0_y_0_xxxyy_0[k] * ab_x + g_y_0_y_0_xxxyy_x[k];
            }

            /// Set up 116-117 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyz_0, g_y_0_y_0_xxxyz_0, g_y_0_y_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyz_0[k] = -g_y_0_y_0_xxxyz_0[k] * ab_x + g_y_0_y_0_xxxyz_x[k];
            }

            /// Set up 117-118 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 117 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxzz_0, g_y_0_y_0_xxxzz_0, g_y_0_y_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxzz_0[k] = -g_y_0_y_0_xxxzz_0[k] * ab_x + g_y_0_y_0_xxxzz_x[k];
            }

            /// Set up 118-119 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 118 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyy_0, g_y_0_y_0_xxyyy_0, g_y_0_y_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyy_0[k] = -g_y_0_y_0_xxyyy_0[k] * ab_x + g_y_0_y_0_xxyyy_x[k];
            }

            /// Set up 119-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyz_0, g_y_0_y_0_xxyyz_0, g_y_0_y_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyz_0[k] = -g_y_0_y_0_xxyyz_0[k] * ab_x + g_y_0_y_0_xxyyz_x[k];
            }

            /// Set up 120-121 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 120 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyzz_0, g_y_0_y_0_xxyzz_0, g_y_0_y_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyzz_0[k] = -g_y_0_y_0_xxyzz_0[k] * ab_x + g_y_0_y_0_xxyzz_x[k];
            }

            /// Set up 121-122 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 121 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxzzz_0, g_y_0_y_0_xxzzz_0, g_y_0_y_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxzzz_0[k] = -g_y_0_y_0_xxzzz_0[k] * ab_x + g_y_0_y_0_xxzzz_x[k];
            }

            /// Set up 122-123 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyy_0, g_y_0_y_0_xyyyy_0, g_y_0_y_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyy_0[k] = -g_y_0_y_0_xyyyy_0[k] * ab_x + g_y_0_y_0_xyyyy_x[k];
            }

            /// Set up 123-124 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 123 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyz_0, g_y_0_y_0_xyyyz_0, g_y_0_y_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyz_0[k] = -g_y_0_y_0_xyyyz_0[k] * ab_x + g_y_0_y_0_xyyyz_x[k];
            }

            /// Set up 124-125 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 124 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyzz_0, g_y_0_y_0_xyyzz_0, g_y_0_y_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyzz_0[k] = -g_y_0_y_0_xyyzz_0[k] * ab_x + g_y_0_y_0_xyyzz_x[k];
            }

            /// Set up 125-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyzzz_0, g_y_0_y_0_xyzzz_0, g_y_0_y_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyzzz_0[k] = -g_y_0_y_0_xyzzz_0[k] * ab_x + g_y_0_y_0_xyzzz_x[k];
            }

            /// Set up 126-127 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 126 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzzzz_0, g_y_0_y_0_xzzzz_0, g_y_0_y_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzzzz_0[k] = -g_y_0_y_0_xzzzz_0[k] * ab_x + g_y_0_y_0_xzzzz_x[k];
            }

            /// Set up 127-128 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 127 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyy_0, g_y_0_y_0_yyyyy_0, g_y_0_y_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyy_0[k] = -g_y_0_y_0_yyyyy_0[k] * ab_x + g_y_0_y_0_yyyyy_x[k];
            }

            /// Set up 128-129 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyz_0, g_y_0_y_0_yyyyz_0, g_y_0_y_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyz_0[k] = -g_y_0_y_0_yyyyz_0[k] * ab_x + g_y_0_y_0_yyyyz_x[k];
            }

            /// Set up 129-130 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyzz_0, g_y_0_y_0_yyyzz_0, g_y_0_y_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyzz_0[k] = -g_y_0_y_0_yyyzz_0[k] * ab_x + g_y_0_y_0_yyyzz_x[k];
            }

            /// Set up 130-131 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 130 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyzzz_0, g_y_0_y_0_yyzzz_0, g_y_0_y_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyzzz_0[k] = -g_y_0_y_0_yyzzz_0[k] * ab_x + g_y_0_y_0_yyzzz_x[k];
            }

            /// Set up 131-132 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzzzz_0, g_y_0_y_0_yzzzz_0, g_y_0_y_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzzzz_0[k] = -g_y_0_y_0_yzzzz_0[k] * ab_x + g_y_0_y_0_yzzzz_x[k];
            }

            /// Set up 132-133 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 132 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzzzz_0, g_y_0_y_0_zzzzz_0, g_y_0_y_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzzzz_0[k] = -g_y_0_y_0_zzzzz_0[k] * ab_x + g_y_0_y_0_zzzzz_x[k];
            }

            /// Set up 133-134 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 133 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyyyy_0, g_y_0_y_0_yyyyy_0, g_y_0_y_0_yyyyy_y, g_y_0_y_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyy_0[k] = -g_0_0_y_0_yyyyy_0[k] - g_y_0_y_0_yyyyy_0[k] * ab_y + g_y_0_y_0_yyyyy_y[k];
            }

            /// Set up 134-135 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyy_0, g_y_0_y_0_yyyyy_z, g_y_0_y_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyz_0[k] = -g_y_0_y_0_yyyyy_0[k] * ab_z + g_y_0_y_0_yyyyy_z[k];
            }

            /// Set up 135-136 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 135 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyz_0, g_y_0_y_0_yyyyz_z, g_y_0_y_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyzz_0[k] = -g_y_0_y_0_yyyyz_0[k] * ab_z + g_y_0_y_0_yyyyz_z[k];
            }

            /// Set up 136-137 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 136 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyzz_0, g_y_0_y_0_yyyzz_z, g_y_0_y_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyzzz_0[k] = -g_y_0_y_0_yyyzz_0[k] * ab_z + g_y_0_y_0_yyyzz_z[k];
            }

            /// Set up 137-138 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyzzz_0, g_y_0_y_0_yyzzz_z, g_y_0_y_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzzzz_0[k] = -g_y_0_y_0_yyzzz_0[k] * ab_z + g_y_0_y_0_yyzzz_z[k];
            }

            /// Set up 138-139 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 138 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzzzz_0, g_y_0_y_0_yzzzz_z, g_y_0_y_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzzzz_0[k] = -g_y_0_y_0_yzzzz_0[k] * ab_z + g_y_0_y_0_yzzzz_z[k];
            }

            /// Set up 139-140 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzzzz_0, g_y_0_y_0_zzzzz_z, g_y_0_y_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzzzz_0[k] = -g_y_0_y_0_zzzzz_0[k] * ab_z + g_y_0_y_0_zzzzz_z[k];
            }

            /// Set up 140-141 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxx_0, g_y_0_z_0_xxxxx_x, g_y_0_z_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxx_0[k] = -g_y_0_z_0_xxxxx_0[k] * ab_x + g_y_0_z_0_xxxxx_x[k];
            }

            /// Set up 141-142 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 141 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxy_0, g_y_0_z_0_xxxxy_0, g_y_0_z_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxy_0[k] = -g_y_0_z_0_xxxxy_0[k] * ab_x + g_y_0_z_0_xxxxy_x[k];
            }

            /// Set up 142-143 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 142 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxz_0, g_y_0_z_0_xxxxz_0, g_y_0_z_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxz_0[k] = -g_y_0_z_0_xxxxz_0[k] * ab_x + g_y_0_z_0_xxxxz_x[k];
            }

            /// Set up 143-144 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyy_0, g_y_0_z_0_xxxyy_0, g_y_0_z_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyy_0[k] = -g_y_0_z_0_xxxyy_0[k] * ab_x + g_y_0_z_0_xxxyy_x[k];
            }

            /// Set up 144-145 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 144 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyz_0, g_y_0_z_0_xxxyz_0, g_y_0_z_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyz_0[k] = -g_y_0_z_0_xxxyz_0[k] * ab_x + g_y_0_z_0_xxxyz_x[k];
            }

            /// Set up 145-146 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 145 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxzz_0, g_y_0_z_0_xxxzz_0, g_y_0_z_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxzz_0[k] = -g_y_0_z_0_xxxzz_0[k] * ab_x + g_y_0_z_0_xxxzz_x[k];
            }

            /// Set up 146-147 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyy_0, g_y_0_z_0_xxyyy_0, g_y_0_z_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyy_0[k] = -g_y_0_z_0_xxyyy_0[k] * ab_x + g_y_0_z_0_xxyyy_x[k];
            }

            /// Set up 147-148 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 147 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyz_0, g_y_0_z_0_xxyyz_0, g_y_0_z_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyz_0[k] = -g_y_0_z_0_xxyyz_0[k] * ab_x + g_y_0_z_0_xxyyz_x[k];
            }

            /// Set up 148-149 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 148 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyzz_0, g_y_0_z_0_xxyzz_0, g_y_0_z_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyzz_0[k] = -g_y_0_z_0_xxyzz_0[k] * ab_x + g_y_0_z_0_xxyzz_x[k];
            }

            /// Set up 149-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxzzz_0, g_y_0_z_0_xxzzz_0, g_y_0_z_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxzzz_0[k] = -g_y_0_z_0_xxzzz_0[k] * ab_x + g_y_0_z_0_xxzzz_x[k];
            }

            /// Set up 150-151 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 150 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyy_0, g_y_0_z_0_xyyyy_0, g_y_0_z_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyy_0[k] = -g_y_0_z_0_xyyyy_0[k] * ab_x + g_y_0_z_0_xyyyy_x[k];
            }

            /// Set up 151-152 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 151 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyz_0, g_y_0_z_0_xyyyz_0, g_y_0_z_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyz_0[k] = -g_y_0_z_0_xyyyz_0[k] * ab_x + g_y_0_z_0_xyyyz_x[k];
            }

            /// Set up 152-153 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyzz_0, g_y_0_z_0_xyyzz_0, g_y_0_z_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyzz_0[k] = -g_y_0_z_0_xyyzz_0[k] * ab_x + g_y_0_z_0_xyyzz_x[k];
            }

            /// Set up 153-154 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 153 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyzzz_0, g_y_0_z_0_xyzzz_0, g_y_0_z_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyzzz_0[k] = -g_y_0_z_0_xyzzz_0[k] * ab_x + g_y_0_z_0_xyzzz_x[k];
            }

            /// Set up 154-155 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 154 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzzzz_0, g_y_0_z_0_xzzzz_0, g_y_0_z_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzzzz_0[k] = -g_y_0_z_0_xzzzz_0[k] * ab_x + g_y_0_z_0_xzzzz_x[k];
            }

            /// Set up 155-156 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyy_0, g_y_0_z_0_yyyyy_0, g_y_0_z_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyy_0[k] = -g_y_0_z_0_yyyyy_0[k] * ab_x + g_y_0_z_0_yyyyy_x[k];
            }

            /// Set up 156-157 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 156 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyz_0, g_y_0_z_0_yyyyz_0, g_y_0_z_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyz_0[k] = -g_y_0_z_0_yyyyz_0[k] * ab_x + g_y_0_z_0_yyyyz_x[k];
            }

            /// Set up 157-158 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 157 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyzz_0, g_y_0_z_0_yyyzz_0, g_y_0_z_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyzz_0[k] = -g_y_0_z_0_yyyzz_0[k] * ab_x + g_y_0_z_0_yyyzz_x[k];
            }

            /// Set up 158-159 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyzzz_0, g_y_0_z_0_yyzzz_0, g_y_0_z_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyzzz_0[k] = -g_y_0_z_0_yyzzz_0[k] * ab_x + g_y_0_z_0_yyzzz_x[k];
            }

            /// Set up 159-160 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzzzz_0, g_y_0_z_0_yzzzz_0, g_y_0_z_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzzzz_0[k] = -g_y_0_z_0_yzzzz_0[k] * ab_x + g_y_0_z_0_yzzzz_x[k];
            }

            /// Set up 160-161 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 160 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzzzz_0, g_y_0_z_0_zzzzz_0, g_y_0_z_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzzzz_0[k] = -g_y_0_z_0_zzzzz_0[k] * ab_x + g_y_0_z_0_zzzzz_x[k];
            }

            /// Set up 161-162 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyyyy_0, g_y_0_z_0_yyyyy_0, g_y_0_z_0_yyyyy_y, g_y_0_z_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyy_0[k] = -g_0_0_z_0_yyyyy_0[k] - g_y_0_z_0_yyyyy_0[k] * ab_y + g_y_0_z_0_yyyyy_y[k];
            }

            /// Set up 162-163 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 162 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyy_0, g_y_0_z_0_yyyyy_z, g_y_0_z_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyz_0[k] = -g_y_0_z_0_yyyyy_0[k] * ab_z + g_y_0_z_0_yyyyy_z[k];
            }

            /// Set up 163-164 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 163 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyz_0, g_y_0_z_0_yyyyz_z, g_y_0_z_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyzz_0[k] = -g_y_0_z_0_yyyyz_0[k] * ab_z + g_y_0_z_0_yyyyz_z[k];
            }

            /// Set up 164-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyzz_0, g_y_0_z_0_yyyzz_z, g_y_0_z_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyzzz_0[k] = -g_y_0_z_0_yyyzz_0[k] * ab_z + g_y_0_z_0_yyyzz_z[k];
            }

            /// Set up 165-166 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 165 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyzzz_0, g_y_0_z_0_yyzzz_z, g_y_0_z_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzzzz_0[k] = -g_y_0_z_0_yyzzz_0[k] * ab_z + g_y_0_z_0_yyzzz_z[k];
            }

            /// Set up 166-167 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 166 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzzzz_0, g_y_0_z_0_yzzzz_z, g_y_0_z_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzzzz_0[k] = -g_y_0_z_0_yzzzz_0[k] * ab_z + g_y_0_z_0_yzzzz_z[k];
            }

            /// Set up 167-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzzzz_0, g_y_0_z_0_zzzzz_z, g_y_0_z_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzzzz_0[k] = -g_y_0_z_0_zzzzz_0[k] * ab_z + g_y_0_z_0_zzzzz_z[k];
            }

            /// Set up 168-169 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 168 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxx_0, g_z_0_x_0_xxxxx_x, g_z_0_x_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxx_0[k] = -g_z_0_x_0_xxxxx_0[k] * ab_x + g_z_0_x_0_xxxxx_x[k];
            }

            /// Set up 169-170 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxy_0, g_z_0_x_0_xxxxy_0, g_z_0_x_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxy_0[k] = -g_z_0_x_0_xxxxy_0[k] * ab_x + g_z_0_x_0_xxxxy_x[k];
            }

            /// Set up 170-171 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxz_0, g_z_0_x_0_xxxxz_0, g_z_0_x_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxz_0[k] = -g_z_0_x_0_xxxxz_0[k] * ab_x + g_z_0_x_0_xxxxz_x[k];
            }

            /// Set up 171-172 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 171 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyy_0, g_z_0_x_0_xxxyy_0, g_z_0_x_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyy_0[k] = -g_z_0_x_0_xxxyy_0[k] * ab_x + g_z_0_x_0_xxxyy_x[k];
            }

            /// Set up 172-173 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 172 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyz_0, g_z_0_x_0_xxxyz_0, g_z_0_x_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyz_0[k] = -g_z_0_x_0_xxxyz_0[k] * ab_x + g_z_0_x_0_xxxyz_x[k];
            }

            /// Set up 173-174 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxzz_0, g_z_0_x_0_xxxzz_0, g_z_0_x_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxzz_0[k] = -g_z_0_x_0_xxxzz_0[k] * ab_x + g_z_0_x_0_xxxzz_x[k];
            }

            /// Set up 174-175 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 174 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyy_0, g_z_0_x_0_xxyyy_0, g_z_0_x_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyy_0[k] = -g_z_0_x_0_xxyyy_0[k] * ab_x + g_z_0_x_0_xxyyy_x[k];
            }

            /// Set up 175-176 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 175 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyz_0, g_z_0_x_0_xxyyz_0, g_z_0_x_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyz_0[k] = -g_z_0_x_0_xxyyz_0[k] * ab_x + g_z_0_x_0_xxyyz_x[k];
            }

            /// Set up 176-177 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyzz_0, g_z_0_x_0_xxyzz_0, g_z_0_x_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyzz_0[k] = -g_z_0_x_0_xxyzz_0[k] * ab_x + g_z_0_x_0_xxyzz_x[k];
            }

            /// Set up 177-178 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 177 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxzzz_0, g_z_0_x_0_xxzzz_0, g_z_0_x_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxzzz_0[k] = -g_z_0_x_0_xxzzz_0[k] * ab_x + g_z_0_x_0_xxzzz_x[k];
            }

            /// Set up 178-179 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 178 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyy_0, g_z_0_x_0_xyyyy_0, g_z_0_x_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyy_0[k] = -g_z_0_x_0_xyyyy_0[k] * ab_x + g_z_0_x_0_xyyyy_x[k];
            }

            /// Set up 179-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyz_0, g_z_0_x_0_xyyyz_0, g_z_0_x_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyz_0[k] = -g_z_0_x_0_xyyyz_0[k] * ab_x + g_z_0_x_0_xyyyz_x[k];
            }

            /// Set up 180-181 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 180 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyzz_0, g_z_0_x_0_xyyzz_0, g_z_0_x_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyzz_0[k] = -g_z_0_x_0_xyyzz_0[k] * ab_x + g_z_0_x_0_xyyzz_x[k];
            }

            /// Set up 181-182 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 181 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyzzz_0, g_z_0_x_0_xyzzz_0, g_z_0_x_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyzzz_0[k] = -g_z_0_x_0_xyzzz_0[k] * ab_x + g_z_0_x_0_xyzzz_x[k];
            }

            /// Set up 182-183 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzzzz_0, g_z_0_x_0_xzzzz_0, g_z_0_x_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzzzz_0[k] = -g_z_0_x_0_xzzzz_0[k] * ab_x + g_z_0_x_0_xzzzz_x[k];
            }

            /// Set up 183-184 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 183 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyy_0, g_z_0_x_0_yyyyy_0, g_z_0_x_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyy_0[k] = -g_z_0_x_0_yyyyy_0[k] * ab_x + g_z_0_x_0_yyyyy_x[k];
            }

            /// Set up 184-185 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 184 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyz_0, g_z_0_x_0_yyyyz_0, g_z_0_x_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyz_0[k] = -g_z_0_x_0_yyyyz_0[k] * ab_x + g_z_0_x_0_yyyyz_x[k];
            }

            /// Set up 185-186 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyzz_0, g_z_0_x_0_yyyzz_0, g_z_0_x_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyzz_0[k] = -g_z_0_x_0_yyyzz_0[k] * ab_x + g_z_0_x_0_yyyzz_x[k];
            }

            /// Set up 186-187 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 186 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyzzz_0, g_z_0_x_0_yyzzz_0, g_z_0_x_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyzzz_0[k] = -g_z_0_x_0_yyzzz_0[k] * ab_x + g_z_0_x_0_yyzzz_x[k];
            }

            /// Set up 187-188 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 187 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzzzz_0, g_z_0_x_0_yzzzz_0, g_z_0_x_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzzzz_0[k] = -g_z_0_x_0_yzzzz_0[k] * ab_x + g_z_0_x_0_yzzzz_x[k];
            }

            /// Set up 188-189 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzzzz_0, g_z_0_x_0_zzzzz_0, g_z_0_x_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzzzz_0[k] = -g_z_0_x_0_zzzzz_0[k] * ab_x + g_z_0_x_0_zzzzz_x[k];
            }

            /// Set up 189-190 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyy_0, g_z_0_x_0_yyyyy_y, g_z_0_x_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyy_0[k] = -g_z_0_x_0_yyyyy_0[k] * ab_y + g_z_0_x_0_yyyyy_y[k];
            }

            /// Set up 190-191 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 190 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyyz_0, g_z_0_x_0_yyyyz_0, g_z_0_x_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyz_0[k] = -g_z_0_x_0_yyyyz_0[k] * ab_y + g_z_0_x_0_yyyyz_y[k];
            }

            /// Set up 191-192 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyzz_0, g_z_0_x_0_yyyzz_0, g_z_0_x_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyzz_0[k] = -g_z_0_x_0_yyyzz_0[k] * ab_y + g_z_0_x_0_yyyzz_y[k];
            }

            /// Set up 192-193 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 192 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyzzz_0, g_z_0_x_0_yyzzz_0, g_z_0_x_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyzzz_0[k] = -g_z_0_x_0_yyzzz_0[k] * ab_y + g_z_0_x_0_yyzzz_y[k];
            }

            /// Set up 193-194 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 193 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzzzz_0, g_z_0_x_0_yzzzz_0, g_z_0_x_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzzzz_0[k] = -g_z_0_x_0_yzzzz_0[k] * ab_y + g_z_0_x_0_yzzzz_y[k];
            }

            /// Set up 194-195 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzzzz_0, g_z_0_x_0_zzzzz_0, g_z_0_x_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzzzz_0[k] = -g_z_0_x_0_zzzzz_0[k] * ab_y + g_z_0_x_0_zzzzz_y[k];
            }

            /// Set up 195-196 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzzzz_0, g_z_0_x_0_zzzzz_0, g_z_0_x_0_zzzzz_z, g_z_0_x_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzzzz_0[k] = -g_0_0_x_0_zzzzz_0[k] - g_z_0_x_0_zzzzz_0[k] * ab_z + g_z_0_x_0_zzzzz_z[k];
            }

            /// Set up 196-197 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 196 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxx_0, g_z_0_y_0_xxxxx_x, g_z_0_y_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxx_0[k] = -g_z_0_y_0_xxxxx_0[k] * ab_x + g_z_0_y_0_xxxxx_x[k];
            }

            /// Set up 197-198 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxy_0, g_z_0_y_0_xxxxy_0, g_z_0_y_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxy_0[k] = -g_z_0_y_0_xxxxy_0[k] * ab_x + g_z_0_y_0_xxxxy_x[k];
            }

            /// Set up 198-199 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 198 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxz_0, g_z_0_y_0_xxxxz_0, g_z_0_y_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxz_0[k] = -g_z_0_y_0_xxxxz_0[k] * ab_x + g_z_0_y_0_xxxxz_x[k];
            }

            /// Set up 199-200 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyy_0, g_z_0_y_0_xxxyy_0, g_z_0_y_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyy_0[k] = -g_z_0_y_0_xxxyy_0[k] * ab_x + g_z_0_y_0_xxxyy_x[k];
            }

            /// Set up 200-201 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyz_0, g_z_0_y_0_xxxyz_0, g_z_0_y_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyz_0[k] = -g_z_0_y_0_xxxyz_0[k] * ab_x + g_z_0_y_0_xxxyz_x[k];
            }

            /// Set up 201-202 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 201 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxzz_0, g_z_0_y_0_xxxzz_0, g_z_0_y_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxzz_0[k] = -g_z_0_y_0_xxxzz_0[k] * ab_x + g_z_0_y_0_xxxzz_x[k];
            }

            /// Set up 202-203 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 202 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyy_0, g_z_0_y_0_xxyyy_0, g_z_0_y_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyy_0[k] = -g_z_0_y_0_xxyyy_0[k] * ab_x + g_z_0_y_0_xxyyy_x[k];
            }

            /// Set up 203-204 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyz_0, g_z_0_y_0_xxyyz_0, g_z_0_y_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyz_0[k] = -g_z_0_y_0_xxyyz_0[k] * ab_x + g_z_0_y_0_xxyyz_x[k];
            }

            /// Set up 204-205 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 204 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyzz_0, g_z_0_y_0_xxyzz_0, g_z_0_y_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyzz_0[k] = -g_z_0_y_0_xxyzz_0[k] * ab_x + g_z_0_y_0_xxyzz_x[k];
            }

            /// Set up 205-206 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 205 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxzzz_0, g_z_0_y_0_xxzzz_0, g_z_0_y_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxzzz_0[k] = -g_z_0_y_0_xxzzz_0[k] * ab_x + g_z_0_y_0_xxzzz_x[k];
            }

            /// Set up 206-207 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyy_0, g_z_0_y_0_xyyyy_0, g_z_0_y_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyy_0[k] = -g_z_0_y_0_xyyyy_0[k] * ab_x + g_z_0_y_0_xyyyy_x[k];
            }

            /// Set up 207-208 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 207 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyz_0, g_z_0_y_0_xyyyz_0, g_z_0_y_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyz_0[k] = -g_z_0_y_0_xyyyz_0[k] * ab_x + g_z_0_y_0_xyyyz_x[k];
            }

            /// Set up 208-209 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 208 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyzz_0, g_z_0_y_0_xyyzz_0, g_z_0_y_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyzz_0[k] = -g_z_0_y_0_xyyzz_0[k] * ab_x + g_z_0_y_0_xyyzz_x[k];
            }

            /// Set up 209-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyzzz_0, g_z_0_y_0_xyzzz_0, g_z_0_y_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyzzz_0[k] = -g_z_0_y_0_xyzzz_0[k] * ab_x + g_z_0_y_0_xyzzz_x[k];
            }

            /// Set up 210-211 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 210 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzzzz_0, g_z_0_y_0_xzzzz_0, g_z_0_y_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzzzz_0[k] = -g_z_0_y_0_xzzzz_0[k] * ab_x + g_z_0_y_0_xzzzz_x[k];
            }

            /// Set up 211-212 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 211 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyy_0, g_z_0_y_0_yyyyy_0, g_z_0_y_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyy_0[k] = -g_z_0_y_0_yyyyy_0[k] * ab_x + g_z_0_y_0_yyyyy_x[k];
            }

            /// Set up 212-213 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyz_0, g_z_0_y_0_yyyyz_0, g_z_0_y_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyz_0[k] = -g_z_0_y_0_yyyyz_0[k] * ab_x + g_z_0_y_0_yyyyz_x[k];
            }

            /// Set up 213-214 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 213 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyzz_0, g_z_0_y_0_yyyzz_0, g_z_0_y_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyzz_0[k] = -g_z_0_y_0_yyyzz_0[k] * ab_x + g_z_0_y_0_yyyzz_x[k];
            }

            /// Set up 214-215 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 214 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyzzz_0, g_z_0_y_0_yyzzz_0, g_z_0_y_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyzzz_0[k] = -g_z_0_y_0_yyzzz_0[k] * ab_x + g_z_0_y_0_yyzzz_x[k];
            }

            /// Set up 215-216 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzzzz_0, g_z_0_y_0_yzzzz_0, g_z_0_y_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzzzz_0[k] = -g_z_0_y_0_yzzzz_0[k] * ab_x + g_z_0_y_0_yzzzz_x[k];
            }

            /// Set up 216-217 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 216 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzzzz_0, g_z_0_y_0_zzzzz_0, g_z_0_y_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzzzz_0[k] = -g_z_0_y_0_zzzzz_0[k] * ab_x + g_z_0_y_0_zzzzz_x[k];
            }

            /// Set up 217-218 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 217 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyy_0, g_z_0_y_0_yyyyy_y, g_z_0_y_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyy_0[k] = -g_z_0_y_0_yyyyy_0[k] * ab_y + g_z_0_y_0_yyyyy_y[k];
            }

            /// Set up 218-219 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyyz_0, g_z_0_y_0_yyyyz_0, g_z_0_y_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyz_0[k] = -g_z_0_y_0_yyyyz_0[k] * ab_y + g_z_0_y_0_yyyyz_y[k];
            }

            /// Set up 219-220 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyzz_0, g_z_0_y_0_yyyzz_0, g_z_0_y_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyzz_0[k] = -g_z_0_y_0_yyyzz_0[k] * ab_y + g_z_0_y_0_yyyzz_y[k];
            }

            /// Set up 220-221 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 220 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyzzz_0, g_z_0_y_0_yyzzz_0, g_z_0_y_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyzzz_0[k] = -g_z_0_y_0_yyzzz_0[k] * ab_y + g_z_0_y_0_yyzzz_y[k];
            }

            /// Set up 221-222 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzzzz_0, g_z_0_y_0_yzzzz_0, g_z_0_y_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzzzz_0[k] = -g_z_0_y_0_yzzzz_0[k] * ab_y + g_z_0_y_0_yzzzz_y[k];
            }

            /// Set up 222-223 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 222 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzzzz_0, g_z_0_y_0_zzzzz_0, g_z_0_y_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzzzz_0[k] = -g_z_0_y_0_zzzzz_0[k] * ab_y + g_z_0_y_0_zzzzz_y[k];
            }

            /// Set up 223-224 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzzzz_0, g_z_0_y_0_zzzzz_0, g_z_0_y_0_zzzzz_z, g_z_0_y_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzzzz_0[k] = -g_0_0_y_0_zzzzz_0[k] - g_z_0_y_0_zzzzz_0[k] * ab_z + g_z_0_y_0_zzzzz_z[k];
            }

            /// Set up 224-225 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxx_0 = cbuffer.data(is_geom_1010_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxx_0, g_z_0_z_0_xxxxx_x, g_z_0_z_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxx_0[k] = -g_z_0_z_0_xxxxx_0[k] * ab_x + g_z_0_z_0_xxxxx_x[k];
            }

            /// Set up 225-226 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxy_0 = cbuffer.data(is_geom_1010_off + 225 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxy_0, g_z_0_z_0_xxxxy_0, g_z_0_z_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxy_0[k] = -g_z_0_z_0_xxxxy_0[k] * ab_x + g_z_0_z_0_xxxxy_x[k];
            }

            /// Set up 226-227 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxz_0 = cbuffer.data(is_geom_1010_off + 226 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxz_0, g_z_0_z_0_xxxxz_0, g_z_0_z_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxz_0[k] = -g_z_0_z_0_xxxxz_0[k] * ab_x + g_z_0_z_0_xxxxz_x[k];
            }

            /// Set up 227-228 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyy_0 = cbuffer.data(is_geom_1010_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyy_0, g_z_0_z_0_xxxyy_0, g_z_0_z_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyy_0[k] = -g_z_0_z_0_xxxyy_0[k] * ab_x + g_z_0_z_0_xxxyy_x[k];
            }

            /// Set up 228-229 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyz_0 = cbuffer.data(is_geom_1010_off + 228 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyz_0, g_z_0_z_0_xxxyz_0, g_z_0_z_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyz_0[k] = -g_z_0_z_0_xxxyz_0[k] * ab_x + g_z_0_z_0_xxxyz_x[k];
            }

            /// Set up 229-230 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxzz_0 = cbuffer.data(is_geom_1010_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxzz_0, g_z_0_z_0_xxxzz_0, g_z_0_z_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxzz_0[k] = -g_z_0_z_0_xxxzz_0[k] * ab_x + g_z_0_z_0_xxxzz_x[k];
            }

            /// Set up 230-231 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyy_0 = cbuffer.data(is_geom_1010_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyy_0, g_z_0_z_0_xxyyy_0, g_z_0_z_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyy_0[k] = -g_z_0_z_0_xxyyy_0[k] * ab_x + g_z_0_z_0_xxyyy_x[k];
            }

            /// Set up 231-232 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyz_0 = cbuffer.data(is_geom_1010_off + 231 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyz_0, g_z_0_z_0_xxyyz_0, g_z_0_z_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyz_0[k] = -g_z_0_z_0_xxyyz_0[k] * ab_x + g_z_0_z_0_xxyyz_x[k];
            }

            /// Set up 232-233 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyzz_0 = cbuffer.data(is_geom_1010_off + 232 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyzz_0, g_z_0_z_0_xxyzz_0, g_z_0_z_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyzz_0[k] = -g_z_0_z_0_xxyzz_0[k] * ab_x + g_z_0_z_0_xxyzz_x[k];
            }

            /// Set up 233-234 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxzzz_0 = cbuffer.data(is_geom_1010_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxzzz_0, g_z_0_z_0_xxzzz_0, g_z_0_z_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxzzz_0[k] = -g_z_0_z_0_xxzzz_0[k] * ab_x + g_z_0_z_0_xxzzz_x[k];
            }

            /// Set up 234-235 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyy_0 = cbuffer.data(is_geom_1010_off + 234 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyy_0, g_z_0_z_0_xyyyy_0, g_z_0_z_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyy_0[k] = -g_z_0_z_0_xyyyy_0[k] * ab_x + g_z_0_z_0_xyyyy_x[k];
            }

            /// Set up 235-236 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyz_0 = cbuffer.data(is_geom_1010_off + 235 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyz_0, g_z_0_z_0_xyyyz_0, g_z_0_z_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyz_0[k] = -g_z_0_z_0_xyyyz_0[k] * ab_x + g_z_0_z_0_xyyyz_x[k];
            }

            /// Set up 236-237 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyzz_0 = cbuffer.data(is_geom_1010_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyzz_0, g_z_0_z_0_xyyzz_0, g_z_0_z_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyzz_0[k] = -g_z_0_z_0_xyyzz_0[k] * ab_x + g_z_0_z_0_xyyzz_x[k];
            }

            /// Set up 237-238 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyzzz_0 = cbuffer.data(is_geom_1010_off + 237 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyzzz_0, g_z_0_z_0_xyzzz_0, g_z_0_z_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyzzz_0[k] = -g_z_0_z_0_xyzzz_0[k] * ab_x + g_z_0_z_0_xyzzz_x[k];
            }

            /// Set up 238-239 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzzzz_0 = cbuffer.data(is_geom_1010_off + 238 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzzzz_0, g_z_0_z_0_xzzzz_0, g_z_0_z_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzzzz_0[k] = -g_z_0_z_0_xzzzz_0[k] * ab_x + g_z_0_z_0_xzzzz_x[k];
            }

            /// Set up 239-240 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyy_0 = cbuffer.data(is_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyy_0, g_z_0_z_0_yyyyy_0, g_z_0_z_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyy_0[k] = -g_z_0_z_0_yyyyy_0[k] * ab_x + g_z_0_z_0_yyyyy_x[k];
            }

            /// Set up 240-241 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyz_0 = cbuffer.data(is_geom_1010_off + 240 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyz_0, g_z_0_z_0_yyyyz_0, g_z_0_z_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyz_0[k] = -g_z_0_z_0_yyyyz_0[k] * ab_x + g_z_0_z_0_yyyyz_x[k];
            }

            /// Set up 241-242 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyzz_0 = cbuffer.data(is_geom_1010_off + 241 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyzz_0, g_z_0_z_0_yyyzz_0, g_z_0_z_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyzz_0[k] = -g_z_0_z_0_yyyzz_0[k] * ab_x + g_z_0_z_0_yyyzz_x[k];
            }

            /// Set up 242-243 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyzzz_0 = cbuffer.data(is_geom_1010_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyzzz_0, g_z_0_z_0_yyzzz_0, g_z_0_z_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyzzz_0[k] = -g_z_0_z_0_yyzzz_0[k] * ab_x + g_z_0_z_0_yyzzz_x[k];
            }

            /// Set up 243-244 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzzzz_0 = cbuffer.data(is_geom_1010_off + 243 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzzzz_0, g_z_0_z_0_yzzzz_0, g_z_0_z_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzzzz_0[k] = -g_z_0_z_0_yzzzz_0[k] * ab_x + g_z_0_z_0_yzzzz_x[k];
            }

            /// Set up 244-245 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzzzz_0 = cbuffer.data(is_geom_1010_off + 244 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzzzz_0, g_z_0_z_0_zzzzz_0, g_z_0_z_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzzzz_0[k] = -g_z_0_z_0_zzzzz_0[k] * ab_x + g_z_0_z_0_zzzzz_x[k];
            }

            /// Set up 245-246 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyy_0 = cbuffer.data(is_geom_1010_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyy_0, g_z_0_z_0_yyyyy_y, g_z_0_z_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyy_0[k] = -g_z_0_z_0_yyyyy_0[k] * ab_y + g_z_0_z_0_yyyyy_y[k];
            }

            /// Set up 246-247 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyz_0 = cbuffer.data(is_geom_1010_off + 246 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyyz_0, g_z_0_z_0_yyyyz_0, g_z_0_z_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyz_0[k] = -g_z_0_z_0_yyyyz_0[k] * ab_y + g_z_0_z_0_yyyyz_y[k];
            }

            /// Set up 247-248 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyzz_0 = cbuffer.data(is_geom_1010_off + 247 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyzz_0, g_z_0_z_0_yyyzz_0, g_z_0_z_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyzz_0[k] = -g_z_0_z_0_yyyzz_0[k] * ab_y + g_z_0_z_0_yyyzz_y[k];
            }

            /// Set up 248-249 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyzzz_0 = cbuffer.data(is_geom_1010_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyzzz_0, g_z_0_z_0_yyzzz_0, g_z_0_z_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyzzz_0[k] = -g_z_0_z_0_yyzzz_0[k] * ab_y + g_z_0_z_0_yyzzz_y[k];
            }

            /// Set up 249-250 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzzzz_0 = cbuffer.data(is_geom_1010_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzzzz_0, g_z_0_z_0_yzzzz_0, g_z_0_z_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzzzz_0[k] = -g_z_0_z_0_yzzzz_0[k] * ab_y + g_z_0_z_0_yzzzz_y[k];
            }

            /// Set up 250-251 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzzzz_0 = cbuffer.data(is_geom_1010_off + 250 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzzzz_0, g_z_0_z_0_zzzzz_0, g_z_0_z_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzzzz_0[k] = -g_z_0_z_0_zzzzz_0[k] * ab_y + g_z_0_z_0_zzzzz_y[k];
            }

            /// Set up 251-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzzzz_0 = cbuffer.data(is_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzzzz_0, g_z_0_z_0_zzzzz_0, g_z_0_z_0_zzzzz_z, g_z_0_z_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzzzz_0[k] = -g_0_0_z_0_zzzzz_0[k] - g_z_0_z_0_zzzzz_0[k] * ab_z + g_z_0_z_0_zzzzz_z[k];
            }
        }
    }
}

} // erirec namespace

