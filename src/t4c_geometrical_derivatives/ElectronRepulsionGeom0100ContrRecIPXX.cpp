#include "ElectronRepulsionGeom0100ContrRecIPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_ipxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ipxx,
                                            const size_t idx_hpxx,
                                            const size_t idx_geom_01_hpxx,
                                            const size_t idx_geom_01_hdxx,
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
            /// Set up components of auxilary buffer : HPSS

            const auto hp_off = idx_hpxx + i * dcomps + j;

            auto g_xxxxx_x = cbuffer.data(hp_off + 0 * ccomps * dcomps);

            auto g_xxxxx_y = cbuffer.data(hp_off + 1 * ccomps * dcomps);

            auto g_xxxxx_z = cbuffer.data(hp_off + 2 * ccomps * dcomps);

            auto g_xxxxy_x = cbuffer.data(hp_off + 3 * ccomps * dcomps);

            auto g_xxxxy_y = cbuffer.data(hp_off + 4 * ccomps * dcomps);

            auto g_xxxxy_z = cbuffer.data(hp_off + 5 * ccomps * dcomps);

            auto g_xxxxz_x = cbuffer.data(hp_off + 6 * ccomps * dcomps);

            auto g_xxxxz_y = cbuffer.data(hp_off + 7 * ccomps * dcomps);

            auto g_xxxxz_z = cbuffer.data(hp_off + 8 * ccomps * dcomps);

            auto g_xxxyy_x = cbuffer.data(hp_off + 9 * ccomps * dcomps);

            auto g_xxxyy_y = cbuffer.data(hp_off + 10 * ccomps * dcomps);

            auto g_xxxyy_z = cbuffer.data(hp_off + 11 * ccomps * dcomps);

            auto g_xxxyz_x = cbuffer.data(hp_off + 12 * ccomps * dcomps);

            auto g_xxxyz_y = cbuffer.data(hp_off + 13 * ccomps * dcomps);

            auto g_xxxyz_z = cbuffer.data(hp_off + 14 * ccomps * dcomps);

            auto g_xxxzz_x = cbuffer.data(hp_off + 15 * ccomps * dcomps);

            auto g_xxxzz_y = cbuffer.data(hp_off + 16 * ccomps * dcomps);

            auto g_xxxzz_z = cbuffer.data(hp_off + 17 * ccomps * dcomps);

            auto g_xxyyy_x = cbuffer.data(hp_off + 18 * ccomps * dcomps);

            auto g_xxyyy_y = cbuffer.data(hp_off + 19 * ccomps * dcomps);

            auto g_xxyyy_z = cbuffer.data(hp_off + 20 * ccomps * dcomps);

            auto g_xxyyz_x = cbuffer.data(hp_off + 21 * ccomps * dcomps);

            auto g_xxyyz_y = cbuffer.data(hp_off + 22 * ccomps * dcomps);

            auto g_xxyyz_z = cbuffer.data(hp_off + 23 * ccomps * dcomps);

            auto g_xxyzz_x = cbuffer.data(hp_off + 24 * ccomps * dcomps);

            auto g_xxyzz_y = cbuffer.data(hp_off + 25 * ccomps * dcomps);

            auto g_xxyzz_z = cbuffer.data(hp_off + 26 * ccomps * dcomps);

            auto g_xxzzz_x = cbuffer.data(hp_off + 27 * ccomps * dcomps);

            auto g_xxzzz_y = cbuffer.data(hp_off + 28 * ccomps * dcomps);

            auto g_xxzzz_z = cbuffer.data(hp_off + 29 * ccomps * dcomps);

            auto g_xyyyy_x = cbuffer.data(hp_off + 30 * ccomps * dcomps);

            auto g_xyyyy_y = cbuffer.data(hp_off + 31 * ccomps * dcomps);

            auto g_xyyyy_z = cbuffer.data(hp_off + 32 * ccomps * dcomps);

            auto g_xyyyz_x = cbuffer.data(hp_off + 33 * ccomps * dcomps);

            auto g_xyyyz_y = cbuffer.data(hp_off + 34 * ccomps * dcomps);

            auto g_xyyyz_z = cbuffer.data(hp_off + 35 * ccomps * dcomps);

            auto g_xyyzz_x = cbuffer.data(hp_off + 36 * ccomps * dcomps);

            auto g_xyyzz_y = cbuffer.data(hp_off + 37 * ccomps * dcomps);

            auto g_xyyzz_z = cbuffer.data(hp_off + 38 * ccomps * dcomps);

            auto g_xyzzz_x = cbuffer.data(hp_off + 39 * ccomps * dcomps);

            auto g_xyzzz_y = cbuffer.data(hp_off + 40 * ccomps * dcomps);

            auto g_xyzzz_z = cbuffer.data(hp_off + 41 * ccomps * dcomps);

            auto g_xzzzz_x = cbuffer.data(hp_off + 42 * ccomps * dcomps);

            auto g_xzzzz_y = cbuffer.data(hp_off + 43 * ccomps * dcomps);

            auto g_xzzzz_z = cbuffer.data(hp_off + 44 * ccomps * dcomps);

            auto g_yyyyy_x = cbuffer.data(hp_off + 45 * ccomps * dcomps);

            auto g_yyyyy_y = cbuffer.data(hp_off + 46 * ccomps * dcomps);

            auto g_yyyyy_z = cbuffer.data(hp_off + 47 * ccomps * dcomps);

            auto g_yyyyz_x = cbuffer.data(hp_off + 48 * ccomps * dcomps);

            auto g_yyyyz_y = cbuffer.data(hp_off + 49 * ccomps * dcomps);

            auto g_yyyyz_z = cbuffer.data(hp_off + 50 * ccomps * dcomps);

            auto g_yyyzz_x = cbuffer.data(hp_off + 51 * ccomps * dcomps);

            auto g_yyyzz_y = cbuffer.data(hp_off + 52 * ccomps * dcomps);

            auto g_yyyzz_z = cbuffer.data(hp_off + 53 * ccomps * dcomps);

            auto g_yyzzz_x = cbuffer.data(hp_off + 54 * ccomps * dcomps);

            auto g_yyzzz_y = cbuffer.data(hp_off + 55 * ccomps * dcomps);

            auto g_yyzzz_z = cbuffer.data(hp_off + 56 * ccomps * dcomps);

            auto g_yzzzz_x = cbuffer.data(hp_off + 57 * ccomps * dcomps);

            auto g_yzzzz_y = cbuffer.data(hp_off + 58 * ccomps * dcomps);

            auto g_yzzzz_z = cbuffer.data(hp_off + 59 * ccomps * dcomps);

            auto g_zzzzz_x = cbuffer.data(hp_off + 60 * ccomps * dcomps);

            auto g_zzzzz_y = cbuffer.data(hp_off + 61 * ccomps * dcomps);

            auto g_zzzzz_z = cbuffer.data(hp_off + 62 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_01_off = idx_geom_01_hpxx + i * dcomps + j;

            auto g_0_x_xxxxx_x = cbuffer.data(hp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_y = cbuffer.data(hp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_z = cbuffer.data(hp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxy_x = cbuffer.data(hp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxy_y = cbuffer.data(hp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxy_z = cbuffer.data(hp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxz_x = cbuffer.data(hp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxz_y = cbuffer.data(hp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxz_z = cbuffer.data(hp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxyy_x = cbuffer.data(hp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxyy_y = cbuffer.data(hp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxyy_z = cbuffer.data(hp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxyz_x = cbuffer.data(hp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxyz_y = cbuffer.data(hp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxyz_z = cbuffer.data(hp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxzz_x = cbuffer.data(hp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxzz_y = cbuffer.data(hp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxzz_z = cbuffer.data(hp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxyyy_x = cbuffer.data(hp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxyyy_y = cbuffer.data(hp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxyyy_z = cbuffer.data(hp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxyyz_x = cbuffer.data(hp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxyyz_y = cbuffer.data(hp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxyyz_z = cbuffer.data(hp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxyzz_x = cbuffer.data(hp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxyzz_y = cbuffer.data(hp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxyzz_z = cbuffer.data(hp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxzzz_x = cbuffer.data(hp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxzzz_y = cbuffer.data(hp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxzzz_z = cbuffer.data(hp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xyyyy_x = cbuffer.data(hp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xyyyy_y = cbuffer.data(hp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xyyyy_z = cbuffer.data(hp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xyyyz_x = cbuffer.data(hp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xyyyz_y = cbuffer.data(hp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xyyyz_z = cbuffer.data(hp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xyyzz_x = cbuffer.data(hp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyyzz_y = cbuffer.data(hp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyyzz_z = cbuffer.data(hp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyzzz_x = cbuffer.data(hp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyzzz_y = cbuffer.data(hp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyzzz_z = cbuffer.data(hp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xzzzz_x = cbuffer.data(hp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xzzzz_y = cbuffer.data(hp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xzzzz_z = cbuffer.data(hp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yyyyy_x = cbuffer.data(hp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yyyyy_y = cbuffer.data(hp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yyyyy_z = cbuffer.data(hp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yyyyz_x = cbuffer.data(hp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yyyyz_y = cbuffer.data(hp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yyyyz_z = cbuffer.data(hp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_yyyzz_x = cbuffer.data(hp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yyyzz_y = cbuffer.data(hp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yyyzz_z = cbuffer.data(hp_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_yyzzz_x = cbuffer.data(hp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_yyzzz_y = cbuffer.data(hp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_yyzzz_z = cbuffer.data(hp_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_yzzzz_x = cbuffer.data(hp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_yzzzz_y = cbuffer.data(hp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_yzzzz_z = cbuffer.data(hp_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_zzzzz_x = cbuffer.data(hp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_zzzzz_y = cbuffer.data(hp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_zzzzz_z = cbuffer.data(hp_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_xxxxx_x = cbuffer.data(hp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xxxxx_y = cbuffer.data(hp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xxxxx_z = cbuffer.data(hp_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_xxxxy_x = cbuffer.data(hp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xxxxy_y = cbuffer.data(hp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xxxxy_z = cbuffer.data(hp_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_xxxxz_x = cbuffer.data(hp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xxxxz_y = cbuffer.data(hp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xxxxz_z = cbuffer.data(hp_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_xxxyy_x = cbuffer.data(hp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xxxyy_y = cbuffer.data(hp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xxxyy_z = cbuffer.data(hp_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_xxxyz_x = cbuffer.data(hp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_xxxyz_y = cbuffer.data(hp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_xxxyz_z = cbuffer.data(hp_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_xxxzz_x = cbuffer.data(hp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_xxxzz_y = cbuffer.data(hp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_xxxzz_z = cbuffer.data(hp_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_xxyyy_x = cbuffer.data(hp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_xxyyy_y = cbuffer.data(hp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_xxyyy_z = cbuffer.data(hp_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_xxyyz_x = cbuffer.data(hp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xxyyz_y = cbuffer.data(hp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xxyyz_z = cbuffer.data(hp_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_xxyzz_x = cbuffer.data(hp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xxyzz_y = cbuffer.data(hp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xxyzz_z = cbuffer.data(hp_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xxzzz_x = cbuffer.data(hp_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxzzz_y = cbuffer.data(hp_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxzzz_z = cbuffer.data(hp_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xyyyy_x = cbuffer.data(hp_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xyyyy_y = cbuffer.data(hp_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xyyyy_z = cbuffer.data(hp_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xyyyz_x = cbuffer.data(hp_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xyyyz_y = cbuffer.data(hp_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xyyyz_z = cbuffer.data(hp_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xyyzz_x = cbuffer.data(hp_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xyyzz_y = cbuffer.data(hp_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xyyzz_z = cbuffer.data(hp_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xyzzz_x = cbuffer.data(hp_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xyzzz_y = cbuffer.data(hp_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xyzzz_z = cbuffer.data(hp_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xzzzz_x = cbuffer.data(hp_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xzzzz_y = cbuffer.data(hp_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xzzzz_z = cbuffer.data(hp_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_yyyyy_x = cbuffer.data(hp_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_yyyyy_y = cbuffer.data(hp_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_yyyyy_z = cbuffer.data(hp_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_yyyyz_x = cbuffer.data(hp_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_yyyyz_y = cbuffer.data(hp_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_yyyyz_z = cbuffer.data(hp_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_yyyzz_x = cbuffer.data(hp_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_yyyzz_y = cbuffer.data(hp_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_yyyzz_z = cbuffer.data(hp_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_yyzzz_x = cbuffer.data(hp_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_yyzzz_y = cbuffer.data(hp_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_yyzzz_z = cbuffer.data(hp_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_yzzzz_x = cbuffer.data(hp_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_yzzzz_y = cbuffer.data(hp_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_yzzzz_z = cbuffer.data(hp_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_zzzzz_x = cbuffer.data(hp_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_zzzzz_y = cbuffer.data(hp_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_zzzzz_z = cbuffer.data(hp_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_xxxxx_x = cbuffer.data(hp_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_xxxxx_y = cbuffer.data(hp_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_xxxxx_z = cbuffer.data(hp_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_xxxxy_x = cbuffer.data(hp_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_xxxxy_y = cbuffer.data(hp_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_xxxxy_z = cbuffer.data(hp_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_xxxxz_x = cbuffer.data(hp_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_xxxxz_y = cbuffer.data(hp_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_xxxxz_z = cbuffer.data(hp_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_z_xxxyy_x = cbuffer.data(hp_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_xxxyy_y = cbuffer.data(hp_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_xxxyy_z = cbuffer.data(hp_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_z_xxxyz_x = cbuffer.data(hp_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_xxxyz_y = cbuffer.data(hp_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_xxxyz_z = cbuffer.data(hp_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_z_xxxzz_x = cbuffer.data(hp_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_xxxzz_y = cbuffer.data(hp_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_xxxzz_z = cbuffer.data(hp_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_z_xxyyy_x = cbuffer.data(hp_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_xxyyy_y = cbuffer.data(hp_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_xxyyy_z = cbuffer.data(hp_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_z_xxyyz_x = cbuffer.data(hp_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_xxyyz_y = cbuffer.data(hp_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_xxyyz_z = cbuffer.data(hp_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_z_xxyzz_x = cbuffer.data(hp_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_xxyzz_y = cbuffer.data(hp_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_xxyzz_z = cbuffer.data(hp_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_z_xxzzz_x = cbuffer.data(hp_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_xxzzz_y = cbuffer.data(hp_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_xxzzz_z = cbuffer.data(hp_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_z_xyyyy_x = cbuffer.data(hp_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_xyyyy_y = cbuffer.data(hp_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_xyyyy_z = cbuffer.data(hp_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_z_xyyyz_x = cbuffer.data(hp_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_xyyyz_y = cbuffer.data(hp_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_xyyyz_z = cbuffer.data(hp_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_z_xyyzz_x = cbuffer.data(hp_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_xyyzz_y = cbuffer.data(hp_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_xyyzz_z = cbuffer.data(hp_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_z_xyzzz_x = cbuffer.data(hp_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_xyzzz_y = cbuffer.data(hp_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_xyzzz_z = cbuffer.data(hp_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_xzzzz_x = cbuffer.data(hp_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_xzzzz_y = cbuffer.data(hp_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_xzzzz_z = cbuffer.data(hp_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_yyyyy_x = cbuffer.data(hp_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_yyyyy_y = cbuffer.data(hp_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_yyyyy_z = cbuffer.data(hp_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_yyyyz_x = cbuffer.data(hp_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_yyyyz_y = cbuffer.data(hp_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_yyyyz_z = cbuffer.data(hp_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_yyyzz_x = cbuffer.data(hp_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_yyyzz_y = cbuffer.data(hp_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_yyyzz_z = cbuffer.data(hp_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_yyzzz_x = cbuffer.data(hp_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_yyzzz_y = cbuffer.data(hp_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_yyzzz_z = cbuffer.data(hp_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_yzzzz_x = cbuffer.data(hp_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_yzzzz_y = cbuffer.data(hp_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_yzzzz_z = cbuffer.data(hp_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_zzzzz_x = cbuffer.data(hp_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_zzzzz_y = cbuffer.data(hp_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_zzzzz_z = cbuffer.data(hp_geom_01_off + 188 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_01_off = idx_geom_01_hdxx + i * dcomps + j;

            auto g_0_x_xxxxx_xx = cbuffer.data(hd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xy = cbuffer.data(hd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xz = cbuffer.data(hd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_yy = cbuffer.data(hd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_yz = cbuffer.data(hd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_zz = cbuffer.data(hd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxy_xx = cbuffer.data(hd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxy_xy = cbuffer.data(hd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxy_xz = cbuffer.data(hd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxy_yy = cbuffer.data(hd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxy_yz = cbuffer.data(hd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxy_zz = cbuffer.data(hd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxz_xx = cbuffer.data(hd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxz_xy = cbuffer.data(hd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxz_xz = cbuffer.data(hd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxz_yy = cbuffer.data(hd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxz_yz = cbuffer.data(hd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxz_zz = cbuffer.data(hd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxyy_xx = cbuffer.data(hd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxyy_xy = cbuffer.data(hd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxyy_xz = cbuffer.data(hd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxyy_yy = cbuffer.data(hd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxyy_yz = cbuffer.data(hd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxyy_zz = cbuffer.data(hd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxyz_xx = cbuffer.data(hd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxyz_xy = cbuffer.data(hd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxyz_xz = cbuffer.data(hd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxyz_yy = cbuffer.data(hd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxyz_yz = cbuffer.data(hd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxyz_zz = cbuffer.data(hd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxzz_xx = cbuffer.data(hd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxzz_xy = cbuffer.data(hd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxzz_xz = cbuffer.data(hd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxzz_yy = cbuffer.data(hd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxzz_yz = cbuffer.data(hd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxzz_zz = cbuffer.data(hd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxyyy_xx = cbuffer.data(hd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyyy_xy = cbuffer.data(hd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyyy_xz = cbuffer.data(hd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyyy_yy = cbuffer.data(hd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyyy_yz = cbuffer.data(hd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyyy_zz = cbuffer.data(hd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxyyz_xx = cbuffer.data(hd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxyyz_xy = cbuffer.data(hd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxyyz_xz = cbuffer.data(hd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxyyz_yy = cbuffer.data(hd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyyz_yz = cbuffer.data(hd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyyz_zz = cbuffer.data(hd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxyzz_xx = cbuffer.data(hd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyzz_xy = cbuffer.data(hd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxyzz_xz = cbuffer.data(hd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxyzz_yy = cbuffer.data(hd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxyzz_yz = cbuffer.data(hd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxyzz_zz = cbuffer.data(hd_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxzzz_xx = cbuffer.data(hd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxzzz_xy = cbuffer.data(hd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxzzz_xz = cbuffer.data(hd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxzzz_yy = cbuffer.data(hd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxzzz_yz = cbuffer.data(hd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxzzz_zz = cbuffer.data(hd_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xyyyy_xx = cbuffer.data(hd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyyyy_xy = cbuffer.data(hd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyyyy_xz = cbuffer.data(hd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyyyy_yy = cbuffer.data(hd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyyy_yz = cbuffer.data(hd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyyy_zz = cbuffer.data(hd_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyyyz_xx = cbuffer.data(hd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyyz_xy = cbuffer.data(hd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyyz_xz = cbuffer.data(hd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyyyz_yy = cbuffer.data(hd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyyyz_yz = cbuffer.data(hd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyyz_zz = cbuffer.data(hd_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyyzz_xx = cbuffer.data(hd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyzz_xy = cbuffer.data(hd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyzz_xz = cbuffer.data(hd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyyzz_yy = cbuffer.data(hd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyzz_yz = cbuffer.data(hd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyzz_zz = cbuffer.data(hd_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyzzz_xx = cbuffer.data(hd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyzzz_xy = cbuffer.data(hd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyzzz_xz = cbuffer.data(hd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyzzz_yy = cbuffer.data(hd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyzzz_yz = cbuffer.data(hd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyzzz_zz = cbuffer.data(hd_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xzzzz_xx = cbuffer.data(hd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzzzz_xy = cbuffer.data(hd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzzzz_xz = cbuffer.data(hd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzzzz_yy = cbuffer.data(hd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzzzz_yz = cbuffer.data(hd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzzzz_zz = cbuffer.data(hd_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yyyyy_xx = cbuffer.data(hd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyyyy_xy = cbuffer.data(hd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyyyy_xz = cbuffer.data(hd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyyyy_yy = cbuffer.data(hd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyyyy_yz = cbuffer.data(hd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyyyy_zz = cbuffer.data(hd_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yyyyz_xx = cbuffer.data(hd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyyyz_xy = cbuffer.data(hd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyyyz_xz = cbuffer.data(hd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyyyz_yy = cbuffer.data(hd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyyyz_yz = cbuffer.data(hd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyyyz_zz = cbuffer.data(hd_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyyzz_xx = cbuffer.data(hd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyyzz_xy = cbuffer.data(hd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyyzz_xz = cbuffer.data(hd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyyzz_yy = cbuffer.data(hd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyyzz_yz = cbuffer.data(hd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyyzz_zz = cbuffer.data(hd_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyzzz_xx = cbuffer.data(hd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyzzz_xy = cbuffer.data(hd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyzzz_xz = cbuffer.data(hd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyzzz_yy = cbuffer.data(hd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyzzz_yz = cbuffer.data(hd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyzzz_zz = cbuffer.data(hd_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yzzzz_xx = cbuffer.data(hd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yzzzz_xy = cbuffer.data(hd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yzzzz_xz = cbuffer.data(hd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yzzzz_yy = cbuffer.data(hd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yzzzz_yz = cbuffer.data(hd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yzzzz_zz = cbuffer.data(hd_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_zzzzz_xx = cbuffer.data(hd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_zzzzz_xy = cbuffer.data(hd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_zzzzz_xz = cbuffer.data(hd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_zzzzz_yy = cbuffer.data(hd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_zzzzz_yz = cbuffer.data(hd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_zzzzz_zz = cbuffer.data(hd_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xxxxx_xx = cbuffer.data(hd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxxxx_xy = cbuffer.data(hd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxxxx_xz = cbuffer.data(hd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xxxxx_yy = cbuffer.data(hd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xxxxx_yz = cbuffer.data(hd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xxxxx_zz = cbuffer.data(hd_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xxxxy_xx = cbuffer.data(hd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xxxxy_xy = cbuffer.data(hd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xxxxy_xz = cbuffer.data(hd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xxxxy_yy = cbuffer.data(hd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xxxxy_yz = cbuffer.data(hd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xxxxy_zz = cbuffer.data(hd_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xxxxz_xx = cbuffer.data(hd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xxxxz_xy = cbuffer.data(hd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xxxxz_xz = cbuffer.data(hd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xxxxz_yy = cbuffer.data(hd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xxxxz_yz = cbuffer.data(hd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xxxxz_zz = cbuffer.data(hd_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xxxyy_xx = cbuffer.data(hd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xxxyy_xy = cbuffer.data(hd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xxxyy_xz = cbuffer.data(hd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xxxyy_yy = cbuffer.data(hd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xxxyy_yz = cbuffer.data(hd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xxxyy_zz = cbuffer.data(hd_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xxxyz_xx = cbuffer.data(hd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxyz_xy = cbuffer.data(hd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxyz_xz = cbuffer.data(hd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxxyz_yy = cbuffer.data(hd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxxyz_yz = cbuffer.data(hd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxxyz_zz = cbuffer.data(hd_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxxzz_xx = cbuffer.data(hd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxxzz_xy = cbuffer.data(hd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxxzz_xz = cbuffer.data(hd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxxzz_yy = cbuffer.data(hd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxxzz_yz = cbuffer.data(hd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxxzz_zz = cbuffer.data(hd_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxyyy_xx = cbuffer.data(hd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxyyy_xy = cbuffer.data(hd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxyyy_xz = cbuffer.data(hd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxyyy_yy = cbuffer.data(hd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxyyy_yz = cbuffer.data(hd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxyyy_zz = cbuffer.data(hd_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxyyz_xx = cbuffer.data(hd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxyyz_xy = cbuffer.data(hd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxyyz_xz = cbuffer.data(hd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxyyz_yy = cbuffer.data(hd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxyyz_yz = cbuffer.data(hd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxyyz_zz = cbuffer.data(hd_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxyzz_xx = cbuffer.data(hd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxyzz_xy = cbuffer.data(hd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxyzz_xz = cbuffer.data(hd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxyzz_yy = cbuffer.data(hd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxyzz_yz = cbuffer.data(hd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxyzz_zz = cbuffer.data(hd_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxzzz_xx = cbuffer.data(hd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxzzz_xy = cbuffer.data(hd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxzzz_xz = cbuffer.data(hd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxzzz_yy = cbuffer.data(hd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxzzz_yz = cbuffer.data(hd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxzzz_zz = cbuffer.data(hd_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xyyyy_xx = cbuffer.data(hd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xyyyy_xy = cbuffer.data(hd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xyyyy_xz = cbuffer.data(hd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xyyyy_yy = cbuffer.data(hd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xyyyy_yz = cbuffer.data(hd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xyyyy_zz = cbuffer.data(hd_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xyyyz_xx = cbuffer.data(hd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xyyyz_xy = cbuffer.data(hd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xyyyz_xz = cbuffer.data(hd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xyyyz_yy = cbuffer.data(hd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyyyz_yz = cbuffer.data(hd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyyyz_zz = cbuffer.data(hd_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xyyzz_xx = cbuffer.data(hd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyyzz_xy = cbuffer.data(hd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyyzz_xz = cbuffer.data(hd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyyzz_yy = cbuffer.data(hd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyyzz_yz = cbuffer.data(hd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyyzz_zz = cbuffer.data(hd_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xyzzz_xx = cbuffer.data(hd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyzzz_xy = cbuffer.data(hd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyzzz_xz = cbuffer.data(hd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyzzz_yy = cbuffer.data(hd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyzzz_yz = cbuffer.data(hd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyzzz_zz = cbuffer.data(hd_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xzzzz_xx = cbuffer.data(hd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xzzzz_xy = cbuffer.data(hd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xzzzz_xz = cbuffer.data(hd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xzzzz_yy = cbuffer.data(hd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xzzzz_yz = cbuffer.data(hd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xzzzz_zz = cbuffer.data(hd_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_yyyyy_xx = cbuffer.data(hd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_yyyyy_xy = cbuffer.data(hd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_yyyyy_xz = cbuffer.data(hd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_yyyyy_yy = cbuffer.data(hd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_yyyyy_yz = cbuffer.data(hd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_yyyyy_zz = cbuffer.data(hd_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_yyyyz_xx = cbuffer.data(hd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_yyyyz_xy = cbuffer.data(hd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_yyyyz_xz = cbuffer.data(hd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_yyyyz_yy = cbuffer.data(hd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_yyyyz_yz = cbuffer.data(hd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_yyyyz_zz = cbuffer.data(hd_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_yyyzz_xx = cbuffer.data(hd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_yyyzz_xy = cbuffer.data(hd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_yyyzz_xz = cbuffer.data(hd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_yyyzz_yy = cbuffer.data(hd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_yyyzz_yz = cbuffer.data(hd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_yyyzz_zz = cbuffer.data(hd_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_yyzzz_xx = cbuffer.data(hd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_yyzzz_xy = cbuffer.data(hd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_yyzzz_xz = cbuffer.data(hd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_yyzzz_yy = cbuffer.data(hd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_yyzzz_yz = cbuffer.data(hd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_yyzzz_zz = cbuffer.data(hd_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_yzzzz_xx = cbuffer.data(hd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yzzzz_xy = cbuffer.data(hd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yzzzz_xz = cbuffer.data(hd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yzzzz_yy = cbuffer.data(hd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yzzzz_yz = cbuffer.data(hd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yzzzz_zz = cbuffer.data(hd_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_zzzzz_xx = cbuffer.data(hd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_zzzzz_xy = cbuffer.data(hd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_zzzzz_xz = cbuffer.data(hd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_zzzzz_yy = cbuffer.data(hd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_zzzzz_yz = cbuffer.data(hd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_zzzzz_zz = cbuffer.data(hd_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_xxxxx_xx = cbuffer.data(hd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xxxxx_xy = cbuffer.data(hd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xxxxx_xz = cbuffer.data(hd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xxxxx_yy = cbuffer.data(hd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xxxxx_yz = cbuffer.data(hd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xxxxx_zz = cbuffer.data(hd_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xxxxy_xx = cbuffer.data(hd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xxxxy_xy = cbuffer.data(hd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xxxxy_xz = cbuffer.data(hd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_xxxxy_yy = cbuffer.data(hd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xxxxy_yz = cbuffer.data(hd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xxxxy_zz = cbuffer.data(hd_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_xxxxz_xx = cbuffer.data(hd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xxxxz_xy = cbuffer.data(hd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xxxxz_xz = cbuffer.data(hd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_xxxxz_yy = cbuffer.data(hd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xxxxz_yz = cbuffer.data(hd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xxxxz_zz = cbuffer.data(hd_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_xxxyy_xx = cbuffer.data(hd_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xxxyy_xy = cbuffer.data(hd_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xxxyy_xz = cbuffer.data(hd_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_xxxyy_yy = cbuffer.data(hd_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xxxyy_yz = cbuffer.data(hd_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xxxyy_zz = cbuffer.data(hd_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_xxxyz_xx = cbuffer.data(hd_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xxxyz_xy = cbuffer.data(hd_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xxxyz_xz = cbuffer.data(hd_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_xxxyz_yy = cbuffer.data(hd_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xxxyz_yz = cbuffer.data(hd_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xxxyz_zz = cbuffer.data(hd_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_xxxzz_xx = cbuffer.data(hd_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xxxzz_xy = cbuffer.data(hd_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xxxzz_xz = cbuffer.data(hd_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_xxxzz_yy = cbuffer.data(hd_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xxxzz_yz = cbuffer.data(hd_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xxxzz_zz = cbuffer.data(hd_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_xxyyy_xx = cbuffer.data(hd_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xxyyy_xy = cbuffer.data(hd_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xxyyy_xz = cbuffer.data(hd_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_xxyyy_yy = cbuffer.data(hd_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xxyyy_yz = cbuffer.data(hd_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xxyyy_zz = cbuffer.data(hd_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_xxyyz_xx = cbuffer.data(hd_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xxyyz_xy = cbuffer.data(hd_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xxyyz_xz = cbuffer.data(hd_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_xxyyz_yy = cbuffer.data(hd_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xxyyz_yz = cbuffer.data(hd_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xxyyz_zz = cbuffer.data(hd_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xxyzz_xx = cbuffer.data(hd_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxyzz_xy = cbuffer.data(hd_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxyzz_xz = cbuffer.data(hd_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxyzz_yy = cbuffer.data(hd_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxyzz_yz = cbuffer.data(hd_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxyzz_zz = cbuffer.data(hd_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxzzz_xx = cbuffer.data(hd_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxzzz_xy = cbuffer.data(hd_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxzzz_xz = cbuffer.data(hd_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxzzz_yy = cbuffer.data(hd_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxzzz_yz = cbuffer.data(hd_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxzzz_zz = cbuffer.data(hd_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xyyyy_xx = cbuffer.data(hd_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xyyyy_xy = cbuffer.data(hd_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xyyyy_xz = cbuffer.data(hd_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xyyyy_yy = cbuffer.data(hd_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xyyyy_yz = cbuffer.data(hd_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xyyyy_zz = cbuffer.data(hd_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xyyyz_xx = cbuffer.data(hd_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xyyyz_xy = cbuffer.data(hd_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xyyyz_xz = cbuffer.data(hd_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xyyyz_yy = cbuffer.data(hd_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xyyyz_yz = cbuffer.data(hd_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xyyyz_zz = cbuffer.data(hd_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xyyzz_xx = cbuffer.data(hd_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xyyzz_xy = cbuffer.data(hd_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xyyzz_xz = cbuffer.data(hd_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xyyzz_yy = cbuffer.data(hd_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xyyzz_yz = cbuffer.data(hd_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xyyzz_zz = cbuffer.data(hd_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_xyzzz_xx = cbuffer.data(hd_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xyzzz_xy = cbuffer.data(hd_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xyzzz_xz = cbuffer.data(hd_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xyzzz_yy = cbuffer.data(hd_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xyzzz_yz = cbuffer.data(hd_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xyzzz_zz = cbuffer.data(hd_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xzzzz_xx = cbuffer.data(hd_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xzzzz_xy = cbuffer.data(hd_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xzzzz_xz = cbuffer.data(hd_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xzzzz_yy = cbuffer.data(hd_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xzzzz_yz = cbuffer.data(hd_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xzzzz_zz = cbuffer.data(hd_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_yyyyy_xx = cbuffer.data(hd_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_yyyyy_xy = cbuffer.data(hd_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_yyyyy_xz = cbuffer.data(hd_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_yyyyy_yy = cbuffer.data(hd_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_yyyyy_yz = cbuffer.data(hd_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_yyyyy_zz = cbuffer.data(hd_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_yyyyz_xx = cbuffer.data(hd_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_yyyyz_xy = cbuffer.data(hd_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_yyyyz_xz = cbuffer.data(hd_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_yyyyz_yy = cbuffer.data(hd_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_yyyyz_yz = cbuffer.data(hd_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_yyyyz_zz = cbuffer.data(hd_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_yyyzz_xx = cbuffer.data(hd_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_yyyzz_xy = cbuffer.data(hd_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_yyyzz_xz = cbuffer.data(hd_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_yyyzz_yy = cbuffer.data(hd_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_yyyzz_yz = cbuffer.data(hd_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_yyyzz_zz = cbuffer.data(hd_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_yyzzz_xx = cbuffer.data(hd_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_yyzzz_xy = cbuffer.data(hd_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_yyzzz_xz = cbuffer.data(hd_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_yyzzz_yy = cbuffer.data(hd_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_yyzzz_yz = cbuffer.data(hd_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_yyzzz_zz = cbuffer.data(hd_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_yzzzz_xx = cbuffer.data(hd_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_yzzzz_xy = cbuffer.data(hd_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_yzzzz_xz = cbuffer.data(hd_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_yzzzz_yy = cbuffer.data(hd_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_yzzzz_yz = cbuffer.data(hd_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_yzzzz_zz = cbuffer.data(hd_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_zzzzz_xx = cbuffer.data(hd_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_zzzzz_xy = cbuffer.data(hd_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_zzzzz_xz = cbuffer.data(hd_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_zzzzz_yy = cbuffer.data(hd_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_zzzzz_yz = cbuffer.data(hd_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_zzzzz_zz = cbuffer.data(hd_geom_01_off + 377 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ipxx

            const auto ip_geom_01_off = idx_geom_01_ipxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxx_x = cbuffer.data(ip_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_y = cbuffer.data(ip_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_z = cbuffer.data(ip_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_x, g_0_x_xxxxx_xx, g_0_x_xxxxx_xy, g_0_x_xxxxx_xz, g_0_x_xxxxx_y, g_0_x_xxxxx_z, g_0_x_xxxxxx_x, g_0_x_xxxxxx_y, g_0_x_xxxxxx_z, g_xxxxx_x, g_xxxxx_y, g_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxx_x[k] = g_xxxxx_x[k] - g_0_x_xxxxx_x[k] * ab_x + g_0_x_xxxxx_xx[k];

                g_0_x_xxxxxx_y[k] = g_xxxxx_y[k] - g_0_x_xxxxx_y[k] * ab_x + g_0_x_xxxxx_xy[k];

                g_0_x_xxxxxx_z[k] = g_xxxxx_z[k] - g_0_x_xxxxx_z[k] * ab_x + g_0_x_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxy_x = cbuffer.data(ip_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxy_y = cbuffer.data(ip_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxy_z = cbuffer.data(ip_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_x, g_0_x_xxxxx_xy, g_0_x_xxxxx_y, g_0_x_xxxxx_yy, g_0_x_xxxxx_yz, g_0_x_xxxxx_z, g_0_x_xxxxxy_x, g_0_x_xxxxxy_y, g_0_x_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxy_x[k] = -g_0_x_xxxxx_x[k] * ab_y + g_0_x_xxxxx_xy[k];

                g_0_x_xxxxxy_y[k] = -g_0_x_xxxxx_y[k] * ab_y + g_0_x_xxxxx_yy[k];

                g_0_x_xxxxxy_z[k] = -g_0_x_xxxxx_z[k] * ab_y + g_0_x_xxxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxz_x = cbuffer.data(ip_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxz_y = cbuffer.data(ip_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxz_z = cbuffer.data(ip_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_x, g_0_x_xxxxx_xz, g_0_x_xxxxx_y, g_0_x_xxxxx_yz, g_0_x_xxxxx_z, g_0_x_xxxxx_zz, g_0_x_xxxxxz_x, g_0_x_xxxxxz_y, g_0_x_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxz_x[k] = -g_0_x_xxxxx_x[k] * ab_z + g_0_x_xxxxx_xz[k];

                g_0_x_xxxxxz_y[k] = -g_0_x_xxxxx_y[k] * ab_z + g_0_x_xxxxx_yz[k];

                g_0_x_xxxxxz_z[k] = -g_0_x_xxxxx_z[k] * ab_z + g_0_x_xxxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyy_x = cbuffer.data(ip_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxyy_y = cbuffer.data(ip_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxyy_z = cbuffer.data(ip_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxy_x, g_0_x_xxxxy_xy, g_0_x_xxxxy_y, g_0_x_xxxxy_yy, g_0_x_xxxxy_yz, g_0_x_xxxxy_z, g_0_x_xxxxyy_x, g_0_x_xxxxyy_y, g_0_x_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyy_x[k] = -g_0_x_xxxxy_x[k] * ab_y + g_0_x_xxxxy_xy[k];

                g_0_x_xxxxyy_y[k] = -g_0_x_xxxxy_y[k] * ab_y + g_0_x_xxxxy_yy[k];

                g_0_x_xxxxyy_z[k] = -g_0_x_xxxxy_z[k] * ab_y + g_0_x_xxxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyz_x = cbuffer.data(ip_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxyz_y = cbuffer.data(ip_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxyz_z = cbuffer.data(ip_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyz_x, g_0_x_xxxxyz_y, g_0_x_xxxxyz_z, g_0_x_xxxxz_x, g_0_x_xxxxz_xy, g_0_x_xxxxz_y, g_0_x_xxxxz_yy, g_0_x_xxxxz_yz, g_0_x_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyz_x[k] = -g_0_x_xxxxz_x[k] * ab_y + g_0_x_xxxxz_xy[k];

                g_0_x_xxxxyz_y[k] = -g_0_x_xxxxz_y[k] * ab_y + g_0_x_xxxxz_yy[k];

                g_0_x_xxxxyz_z[k] = -g_0_x_xxxxz_z[k] * ab_y + g_0_x_xxxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzz_x = cbuffer.data(ip_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxzz_y = cbuffer.data(ip_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxzz_z = cbuffer.data(ip_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxz_x, g_0_x_xxxxz_xz, g_0_x_xxxxz_y, g_0_x_xxxxz_yz, g_0_x_xxxxz_z, g_0_x_xxxxz_zz, g_0_x_xxxxzz_x, g_0_x_xxxxzz_y, g_0_x_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzz_x[k] = -g_0_x_xxxxz_x[k] * ab_z + g_0_x_xxxxz_xz[k];

                g_0_x_xxxxzz_y[k] = -g_0_x_xxxxz_y[k] * ab_z + g_0_x_xxxxz_yz[k];

                g_0_x_xxxxzz_z[k] = -g_0_x_xxxxz_z[k] * ab_z + g_0_x_xxxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyy_x = cbuffer.data(ip_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxyyy_y = cbuffer.data(ip_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxyyy_z = cbuffer.data(ip_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyy_x, g_0_x_xxxyy_xy, g_0_x_xxxyy_y, g_0_x_xxxyy_yy, g_0_x_xxxyy_yz, g_0_x_xxxyy_z, g_0_x_xxxyyy_x, g_0_x_xxxyyy_y, g_0_x_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyy_x[k] = -g_0_x_xxxyy_x[k] * ab_y + g_0_x_xxxyy_xy[k];

                g_0_x_xxxyyy_y[k] = -g_0_x_xxxyy_y[k] * ab_y + g_0_x_xxxyy_yy[k];

                g_0_x_xxxyyy_z[k] = -g_0_x_xxxyy_z[k] * ab_y + g_0_x_xxxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyz_x = cbuffer.data(ip_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxyyz_y = cbuffer.data(ip_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxyyz_z = cbuffer.data(ip_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyz_x, g_0_x_xxxyyz_y, g_0_x_xxxyyz_z, g_0_x_xxxyz_x, g_0_x_xxxyz_xy, g_0_x_xxxyz_y, g_0_x_xxxyz_yy, g_0_x_xxxyz_yz, g_0_x_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyz_x[k] = -g_0_x_xxxyz_x[k] * ab_y + g_0_x_xxxyz_xy[k];

                g_0_x_xxxyyz_y[k] = -g_0_x_xxxyz_y[k] * ab_y + g_0_x_xxxyz_yy[k];

                g_0_x_xxxyyz_z[k] = -g_0_x_xxxyz_z[k] * ab_y + g_0_x_xxxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzz_x = cbuffer.data(ip_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxyzz_y = cbuffer.data(ip_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxyzz_z = cbuffer.data(ip_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzz_x, g_0_x_xxxyzz_y, g_0_x_xxxyzz_z, g_0_x_xxxzz_x, g_0_x_xxxzz_xy, g_0_x_xxxzz_y, g_0_x_xxxzz_yy, g_0_x_xxxzz_yz, g_0_x_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzz_x[k] = -g_0_x_xxxzz_x[k] * ab_y + g_0_x_xxxzz_xy[k];

                g_0_x_xxxyzz_y[k] = -g_0_x_xxxzz_y[k] * ab_y + g_0_x_xxxzz_yy[k];

                g_0_x_xxxyzz_z[k] = -g_0_x_xxxzz_z[k] * ab_y + g_0_x_xxxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzz_x = cbuffer.data(ip_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxzzz_y = cbuffer.data(ip_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxzzz_z = cbuffer.data(ip_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzz_x, g_0_x_xxxzz_xz, g_0_x_xxxzz_y, g_0_x_xxxzz_yz, g_0_x_xxxzz_z, g_0_x_xxxzz_zz, g_0_x_xxxzzz_x, g_0_x_xxxzzz_y, g_0_x_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzz_x[k] = -g_0_x_xxxzz_x[k] * ab_z + g_0_x_xxxzz_xz[k];

                g_0_x_xxxzzz_y[k] = -g_0_x_xxxzz_y[k] * ab_z + g_0_x_xxxzz_yz[k];

                g_0_x_xxxzzz_z[k] = -g_0_x_xxxzz_z[k] * ab_z + g_0_x_xxxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyy_x = cbuffer.data(ip_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxyyyy_y = cbuffer.data(ip_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxyyyy_z = cbuffer.data(ip_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyy_x, g_0_x_xxyyy_xy, g_0_x_xxyyy_y, g_0_x_xxyyy_yy, g_0_x_xxyyy_yz, g_0_x_xxyyy_z, g_0_x_xxyyyy_x, g_0_x_xxyyyy_y, g_0_x_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyy_x[k] = -g_0_x_xxyyy_x[k] * ab_y + g_0_x_xxyyy_xy[k];

                g_0_x_xxyyyy_y[k] = -g_0_x_xxyyy_y[k] * ab_y + g_0_x_xxyyy_yy[k];

                g_0_x_xxyyyy_z[k] = -g_0_x_xxyyy_z[k] * ab_y + g_0_x_xxyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyz_x = cbuffer.data(ip_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxyyyz_y = cbuffer.data(ip_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxyyyz_z = cbuffer.data(ip_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyz_x, g_0_x_xxyyyz_y, g_0_x_xxyyyz_z, g_0_x_xxyyz_x, g_0_x_xxyyz_xy, g_0_x_xxyyz_y, g_0_x_xxyyz_yy, g_0_x_xxyyz_yz, g_0_x_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyz_x[k] = -g_0_x_xxyyz_x[k] * ab_y + g_0_x_xxyyz_xy[k];

                g_0_x_xxyyyz_y[k] = -g_0_x_xxyyz_y[k] * ab_y + g_0_x_xxyyz_yy[k];

                g_0_x_xxyyyz_z[k] = -g_0_x_xxyyz_z[k] * ab_y + g_0_x_xxyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzz_x = cbuffer.data(ip_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyyzz_y = cbuffer.data(ip_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyyzz_z = cbuffer.data(ip_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzz_x, g_0_x_xxyyzz_y, g_0_x_xxyyzz_z, g_0_x_xxyzz_x, g_0_x_xxyzz_xy, g_0_x_xxyzz_y, g_0_x_xxyzz_yy, g_0_x_xxyzz_yz, g_0_x_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzz_x[k] = -g_0_x_xxyzz_x[k] * ab_y + g_0_x_xxyzz_xy[k];

                g_0_x_xxyyzz_y[k] = -g_0_x_xxyzz_y[k] * ab_y + g_0_x_xxyzz_yy[k];

                g_0_x_xxyyzz_z[k] = -g_0_x_xxyzz_z[k] * ab_y + g_0_x_xxyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzz_x = cbuffer.data(ip_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyzzz_y = cbuffer.data(ip_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyzzz_z = cbuffer.data(ip_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzz_x, g_0_x_xxyzzz_y, g_0_x_xxyzzz_z, g_0_x_xxzzz_x, g_0_x_xxzzz_xy, g_0_x_xxzzz_y, g_0_x_xxzzz_yy, g_0_x_xxzzz_yz, g_0_x_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzz_x[k] = -g_0_x_xxzzz_x[k] * ab_y + g_0_x_xxzzz_xy[k];

                g_0_x_xxyzzz_y[k] = -g_0_x_xxzzz_y[k] * ab_y + g_0_x_xxzzz_yy[k];

                g_0_x_xxyzzz_z[k] = -g_0_x_xxzzz_z[k] * ab_y + g_0_x_xxzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzz_x = cbuffer.data(ip_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxzzzz_y = cbuffer.data(ip_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxzzzz_z = cbuffer.data(ip_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzz_x, g_0_x_xxzzz_xz, g_0_x_xxzzz_y, g_0_x_xxzzz_yz, g_0_x_xxzzz_z, g_0_x_xxzzz_zz, g_0_x_xxzzzz_x, g_0_x_xxzzzz_y, g_0_x_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzz_x[k] = -g_0_x_xxzzz_x[k] * ab_z + g_0_x_xxzzz_xz[k];

                g_0_x_xxzzzz_y[k] = -g_0_x_xxzzz_y[k] * ab_z + g_0_x_xxzzz_yz[k];

                g_0_x_xxzzzz_z[k] = -g_0_x_xxzzz_z[k] * ab_z + g_0_x_xxzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyy_x = cbuffer.data(ip_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyyyyy_y = cbuffer.data(ip_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyyyyy_z = cbuffer.data(ip_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyy_x, g_0_x_xyyyy_xy, g_0_x_xyyyy_y, g_0_x_xyyyy_yy, g_0_x_xyyyy_yz, g_0_x_xyyyy_z, g_0_x_xyyyyy_x, g_0_x_xyyyyy_y, g_0_x_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyy_x[k] = -g_0_x_xyyyy_x[k] * ab_y + g_0_x_xyyyy_xy[k];

                g_0_x_xyyyyy_y[k] = -g_0_x_xyyyy_y[k] * ab_y + g_0_x_xyyyy_yy[k];

                g_0_x_xyyyyy_z[k] = -g_0_x_xyyyy_z[k] * ab_y + g_0_x_xyyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyz_x = cbuffer.data(ip_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyyyyz_y = cbuffer.data(ip_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyyyyz_z = cbuffer.data(ip_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyz_x, g_0_x_xyyyyz_y, g_0_x_xyyyyz_z, g_0_x_xyyyz_x, g_0_x_xyyyz_xy, g_0_x_xyyyz_y, g_0_x_xyyyz_yy, g_0_x_xyyyz_yz, g_0_x_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyz_x[k] = -g_0_x_xyyyz_x[k] * ab_y + g_0_x_xyyyz_xy[k];

                g_0_x_xyyyyz_y[k] = -g_0_x_xyyyz_y[k] * ab_y + g_0_x_xyyyz_yy[k];

                g_0_x_xyyyyz_z[k] = -g_0_x_xyyyz_z[k] * ab_y + g_0_x_xyyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzz_x = cbuffer.data(ip_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyyyzz_y = cbuffer.data(ip_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyyyzz_z = cbuffer.data(ip_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzz_x, g_0_x_xyyyzz_y, g_0_x_xyyyzz_z, g_0_x_xyyzz_x, g_0_x_xyyzz_xy, g_0_x_xyyzz_y, g_0_x_xyyzz_yy, g_0_x_xyyzz_yz, g_0_x_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzz_x[k] = -g_0_x_xyyzz_x[k] * ab_y + g_0_x_xyyzz_xy[k];

                g_0_x_xyyyzz_y[k] = -g_0_x_xyyzz_y[k] * ab_y + g_0_x_xyyzz_yy[k];

                g_0_x_xyyyzz_z[k] = -g_0_x_xyyzz_z[k] * ab_y + g_0_x_xyyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzz_x = cbuffer.data(ip_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xyyzzz_y = cbuffer.data(ip_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xyyzzz_z = cbuffer.data(ip_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzz_x, g_0_x_xyyzzz_y, g_0_x_xyyzzz_z, g_0_x_xyzzz_x, g_0_x_xyzzz_xy, g_0_x_xyzzz_y, g_0_x_xyzzz_yy, g_0_x_xyzzz_yz, g_0_x_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzz_x[k] = -g_0_x_xyzzz_x[k] * ab_y + g_0_x_xyzzz_xy[k];

                g_0_x_xyyzzz_y[k] = -g_0_x_xyzzz_y[k] * ab_y + g_0_x_xyzzz_yy[k];

                g_0_x_xyyzzz_z[k] = -g_0_x_xyzzz_z[k] * ab_y + g_0_x_xyzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzz_x = cbuffer.data(ip_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xyzzzz_y = cbuffer.data(ip_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xyzzzz_z = cbuffer.data(ip_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzz_x, g_0_x_xyzzzz_y, g_0_x_xyzzzz_z, g_0_x_xzzzz_x, g_0_x_xzzzz_xy, g_0_x_xzzzz_y, g_0_x_xzzzz_yy, g_0_x_xzzzz_yz, g_0_x_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzz_x[k] = -g_0_x_xzzzz_x[k] * ab_y + g_0_x_xzzzz_xy[k];

                g_0_x_xyzzzz_y[k] = -g_0_x_xzzzz_y[k] * ab_y + g_0_x_xzzzz_yy[k];

                g_0_x_xyzzzz_z[k] = -g_0_x_xzzzz_z[k] * ab_y + g_0_x_xzzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzz_x = cbuffer.data(ip_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xzzzzz_y = cbuffer.data(ip_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xzzzzz_z = cbuffer.data(ip_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzz_x, g_0_x_xzzzz_xz, g_0_x_xzzzz_y, g_0_x_xzzzz_yz, g_0_x_xzzzz_z, g_0_x_xzzzz_zz, g_0_x_xzzzzz_x, g_0_x_xzzzzz_y, g_0_x_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzz_x[k] = -g_0_x_xzzzz_x[k] * ab_z + g_0_x_xzzzz_xz[k];

                g_0_x_xzzzzz_y[k] = -g_0_x_xzzzz_y[k] * ab_z + g_0_x_xzzzz_yz[k];

                g_0_x_xzzzzz_z[k] = -g_0_x_xzzzz_z[k] * ab_z + g_0_x_xzzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyy_x = cbuffer.data(ip_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyyyyy_y = cbuffer.data(ip_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyyyyy_z = cbuffer.data(ip_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyy_x, g_0_x_yyyyy_xy, g_0_x_yyyyy_y, g_0_x_yyyyy_yy, g_0_x_yyyyy_yz, g_0_x_yyyyy_z, g_0_x_yyyyyy_x, g_0_x_yyyyyy_y, g_0_x_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyy_x[k] = -g_0_x_yyyyy_x[k] * ab_y + g_0_x_yyyyy_xy[k];

                g_0_x_yyyyyy_y[k] = -g_0_x_yyyyy_y[k] * ab_y + g_0_x_yyyyy_yy[k];

                g_0_x_yyyyyy_z[k] = -g_0_x_yyyyy_z[k] * ab_y + g_0_x_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyz_x = cbuffer.data(ip_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyyyyz_y = cbuffer.data(ip_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyyyyz_z = cbuffer.data(ip_geom_01_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyz_x, g_0_x_yyyyyz_y, g_0_x_yyyyyz_z, g_0_x_yyyyz_x, g_0_x_yyyyz_xy, g_0_x_yyyyz_y, g_0_x_yyyyz_yy, g_0_x_yyyyz_yz, g_0_x_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyz_x[k] = -g_0_x_yyyyz_x[k] * ab_y + g_0_x_yyyyz_xy[k];

                g_0_x_yyyyyz_y[k] = -g_0_x_yyyyz_y[k] * ab_y + g_0_x_yyyyz_yy[k];

                g_0_x_yyyyyz_z[k] = -g_0_x_yyyyz_z[k] * ab_y + g_0_x_yyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzz_x = cbuffer.data(ip_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyyyzz_y = cbuffer.data(ip_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyyyzz_z = cbuffer.data(ip_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzz_x, g_0_x_yyyyzz_y, g_0_x_yyyyzz_z, g_0_x_yyyzz_x, g_0_x_yyyzz_xy, g_0_x_yyyzz_y, g_0_x_yyyzz_yy, g_0_x_yyyzz_yz, g_0_x_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzz_x[k] = -g_0_x_yyyzz_x[k] * ab_y + g_0_x_yyyzz_xy[k];

                g_0_x_yyyyzz_y[k] = -g_0_x_yyyzz_y[k] * ab_y + g_0_x_yyyzz_yy[k];

                g_0_x_yyyyzz_z[k] = -g_0_x_yyyzz_z[k] * ab_y + g_0_x_yyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzz_x = cbuffer.data(ip_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyyzzz_y = cbuffer.data(ip_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyyzzz_z = cbuffer.data(ip_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzz_x, g_0_x_yyyzzz_y, g_0_x_yyyzzz_z, g_0_x_yyzzz_x, g_0_x_yyzzz_xy, g_0_x_yyzzz_y, g_0_x_yyzzz_yy, g_0_x_yyzzz_yz, g_0_x_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzz_x[k] = -g_0_x_yyzzz_x[k] * ab_y + g_0_x_yyzzz_xy[k];

                g_0_x_yyyzzz_y[k] = -g_0_x_yyzzz_y[k] * ab_y + g_0_x_yyzzz_yy[k];

                g_0_x_yyyzzz_z[k] = -g_0_x_yyzzz_z[k] * ab_y + g_0_x_yyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzz_x = cbuffer.data(ip_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyzzzz_y = cbuffer.data(ip_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyzzzz_z = cbuffer.data(ip_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzz_x, g_0_x_yyzzzz_y, g_0_x_yyzzzz_z, g_0_x_yzzzz_x, g_0_x_yzzzz_xy, g_0_x_yzzzz_y, g_0_x_yzzzz_yy, g_0_x_yzzzz_yz, g_0_x_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzz_x[k] = -g_0_x_yzzzz_x[k] * ab_y + g_0_x_yzzzz_xy[k];

                g_0_x_yyzzzz_y[k] = -g_0_x_yzzzz_y[k] * ab_y + g_0_x_yzzzz_yy[k];

                g_0_x_yyzzzz_z[k] = -g_0_x_yzzzz_z[k] * ab_y + g_0_x_yzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzz_x = cbuffer.data(ip_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yzzzzz_y = cbuffer.data(ip_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzzzzz_z = cbuffer.data(ip_geom_01_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzz_x, g_0_x_yzzzzz_y, g_0_x_yzzzzz_z, g_0_x_zzzzz_x, g_0_x_zzzzz_xy, g_0_x_zzzzz_y, g_0_x_zzzzz_yy, g_0_x_zzzzz_yz, g_0_x_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzz_x[k] = -g_0_x_zzzzz_x[k] * ab_y + g_0_x_zzzzz_xy[k];

                g_0_x_yzzzzz_y[k] = -g_0_x_zzzzz_y[k] * ab_y + g_0_x_zzzzz_yy[k];

                g_0_x_yzzzzz_z[k] = -g_0_x_zzzzz_z[k] * ab_y + g_0_x_zzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzz_x = cbuffer.data(ip_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_zzzzzz_y = cbuffer.data(ip_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_zzzzzz_z = cbuffer.data(ip_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzz_x, g_0_x_zzzzz_xz, g_0_x_zzzzz_y, g_0_x_zzzzz_yz, g_0_x_zzzzz_z, g_0_x_zzzzz_zz, g_0_x_zzzzzz_x, g_0_x_zzzzzz_y, g_0_x_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzz_x[k] = -g_0_x_zzzzz_x[k] * ab_z + g_0_x_zzzzz_xz[k];

                g_0_x_zzzzzz_y[k] = -g_0_x_zzzzz_y[k] * ab_z + g_0_x_zzzzz_yz[k];

                g_0_x_zzzzzz_z[k] = -g_0_x_zzzzz_z[k] * ab_z + g_0_x_zzzzz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxx_x = cbuffer.data(ip_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xxxxxx_y = cbuffer.data(ip_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xxxxxx_z = cbuffer.data(ip_geom_01_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxx_x, g_0_y_xxxxx_xx, g_0_y_xxxxx_xy, g_0_y_xxxxx_xz, g_0_y_xxxxx_y, g_0_y_xxxxx_z, g_0_y_xxxxxx_x, g_0_y_xxxxxx_y, g_0_y_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxx_x[k] = -g_0_y_xxxxx_x[k] * ab_x + g_0_y_xxxxx_xx[k];

                g_0_y_xxxxxx_y[k] = -g_0_y_xxxxx_y[k] * ab_x + g_0_y_xxxxx_xy[k];

                g_0_y_xxxxxx_z[k] = -g_0_y_xxxxx_z[k] * ab_x + g_0_y_xxxxx_xz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxy_x = cbuffer.data(ip_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xxxxxy_y = cbuffer.data(ip_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xxxxxy_z = cbuffer.data(ip_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxy_x, g_0_y_xxxxxy_y, g_0_y_xxxxxy_z, g_0_y_xxxxy_x, g_0_y_xxxxy_xx, g_0_y_xxxxy_xy, g_0_y_xxxxy_xz, g_0_y_xxxxy_y, g_0_y_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxy_x[k] = -g_0_y_xxxxy_x[k] * ab_x + g_0_y_xxxxy_xx[k];

                g_0_y_xxxxxy_y[k] = -g_0_y_xxxxy_y[k] * ab_x + g_0_y_xxxxy_xy[k];

                g_0_y_xxxxxy_z[k] = -g_0_y_xxxxy_z[k] * ab_x + g_0_y_xxxxy_xz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxz_x = cbuffer.data(ip_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxxxxz_y = cbuffer.data(ip_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxxxxz_z = cbuffer.data(ip_geom_01_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxz_x, g_0_y_xxxxxz_y, g_0_y_xxxxxz_z, g_0_y_xxxxz_x, g_0_y_xxxxz_xx, g_0_y_xxxxz_xy, g_0_y_xxxxz_xz, g_0_y_xxxxz_y, g_0_y_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxz_x[k] = -g_0_y_xxxxz_x[k] * ab_x + g_0_y_xxxxz_xx[k];

                g_0_y_xxxxxz_y[k] = -g_0_y_xxxxz_y[k] * ab_x + g_0_y_xxxxz_xy[k];

                g_0_y_xxxxxz_z[k] = -g_0_y_xxxxz_z[k] * ab_x + g_0_y_xxxxz_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyy_x = cbuffer.data(ip_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xxxxyy_y = cbuffer.data(ip_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xxxxyy_z = cbuffer.data(ip_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyy_x, g_0_y_xxxxyy_y, g_0_y_xxxxyy_z, g_0_y_xxxyy_x, g_0_y_xxxyy_xx, g_0_y_xxxyy_xy, g_0_y_xxxyy_xz, g_0_y_xxxyy_y, g_0_y_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyy_x[k] = -g_0_y_xxxyy_x[k] * ab_x + g_0_y_xxxyy_xx[k];

                g_0_y_xxxxyy_y[k] = -g_0_y_xxxyy_y[k] * ab_x + g_0_y_xxxyy_xy[k];

                g_0_y_xxxxyy_z[k] = -g_0_y_xxxyy_z[k] * ab_x + g_0_y_xxxyy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyz_x = cbuffer.data(ip_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xxxxyz_y = cbuffer.data(ip_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xxxxyz_z = cbuffer.data(ip_geom_01_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyz_x, g_0_y_xxxxyz_y, g_0_y_xxxxyz_z, g_0_y_xxxyz_x, g_0_y_xxxyz_xx, g_0_y_xxxyz_xy, g_0_y_xxxyz_xz, g_0_y_xxxyz_y, g_0_y_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyz_x[k] = -g_0_y_xxxyz_x[k] * ab_x + g_0_y_xxxyz_xx[k];

                g_0_y_xxxxyz_y[k] = -g_0_y_xxxyz_y[k] * ab_x + g_0_y_xxxyz_xy[k];

                g_0_y_xxxxyz_z[k] = -g_0_y_xxxyz_z[k] * ab_x + g_0_y_xxxyz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzz_x = cbuffer.data(ip_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxxxzz_y = cbuffer.data(ip_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxxxzz_z = cbuffer.data(ip_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzz_x, g_0_y_xxxxzz_y, g_0_y_xxxxzz_z, g_0_y_xxxzz_x, g_0_y_xxxzz_xx, g_0_y_xxxzz_xy, g_0_y_xxxzz_xz, g_0_y_xxxzz_y, g_0_y_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzz_x[k] = -g_0_y_xxxzz_x[k] * ab_x + g_0_y_xxxzz_xx[k];

                g_0_y_xxxxzz_y[k] = -g_0_y_xxxzz_y[k] * ab_x + g_0_y_xxxzz_xy[k];

                g_0_y_xxxxzz_z[k] = -g_0_y_xxxzz_z[k] * ab_x + g_0_y_xxxzz_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyy_x = cbuffer.data(ip_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxxyyy_y = cbuffer.data(ip_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxxyyy_z = cbuffer.data(ip_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyy_x, g_0_y_xxxyyy_y, g_0_y_xxxyyy_z, g_0_y_xxyyy_x, g_0_y_xxyyy_xx, g_0_y_xxyyy_xy, g_0_y_xxyyy_xz, g_0_y_xxyyy_y, g_0_y_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyy_x[k] = -g_0_y_xxyyy_x[k] * ab_x + g_0_y_xxyyy_xx[k];

                g_0_y_xxxyyy_y[k] = -g_0_y_xxyyy_y[k] * ab_x + g_0_y_xxyyy_xy[k];

                g_0_y_xxxyyy_z[k] = -g_0_y_xxyyy_z[k] * ab_x + g_0_y_xxyyy_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyz_x = cbuffer.data(ip_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxxyyz_y = cbuffer.data(ip_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxxyyz_z = cbuffer.data(ip_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyz_x, g_0_y_xxxyyz_y, g_0_y_xxxyyz_z, g_0_y_xxyyz_x, g_0_y_xxyyz_xx, g_0_y_xxyyz_xy, g_0_y_xxyyz_xz, g_0_y_xxyyz_y, g_0_y_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyz_x[k] = -g_0_y_xxyyz_x[k] * ab_x + g_0_y_xxyyz_xx[k];

                g_0_y_xxxyyz_y[k] = -g_0_y_xxyyz_y[k] * ab_x + g_0_y_xxyyz_xy[k];

                g_0_y_xxxyyz_z[k] = -g_0_y_xxyyz_z[k] * ab_x + g_0_y_xxyyz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzz_x = cbuffer.data(ip_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxxyzz_y = cbuffer.data(ip_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxxyzz_z = cbuffer.data(ip_geom_01_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzz_x, g_0_y_xxxyzz_y, g_0_y_xxxyzz_z, g_0_y_xxyzz_x, g_0_y_xxyzz_xx, g_0_y_xxyzz_xy, g_0_y_xxyzz_xz, g_0_y_xxyzz_y, g_0_y_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzz_x[k] = -g_0_y_xxyzz_x[k] * ab_x + g_0_y_xxyzz_xx[k];

                g_0_y_xxxyzz_y[k] = -g_0_y_xxyzz_y[k] * ab_x + g_0_y_xxyzz_xy[k];

                g_0_y_xxxyzz_z[k] = -g_0_y_xxyzz_z[k] * ab_x + g_0_y_xxyzz_xz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzz_x = cbuffer.data(ip_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxxzzz_y = cbuffer.data(ip_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxxzzz_z = cbuffer.data(ip_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzz_x, g_0_y_xxxzzz_y, g_0_y_xxxzzz_z, g_0_y_xxzzz_x, g_0_y_xxzzz_xx, g_0_y_xxzzz_xy, g_0_y_xxzzz_xz, g_0_y_xxzzz_y, g_0_y_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzz_x[k] = -g_0_y_xxzzz_x[k] * ab_x + g_0_y_xxzzz_xx[k];

                g_0_y_xxxzzz_y[k] = -g_0_y_xxzzz_y[k] * ab_x + g_0_y_xxzzz_xy[k];

                g_0_y_xxxzzz_z[k] = -g_0_y_xxzzz_z[k] * ab_x + g_0_y_xxzzz_xz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyy_x = cbuffer.data(ip_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxyyyy_y = cbuffer.data(ip_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxyyyy_z = cbuffer.data(ip_geom_01_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyy_x, g_0_y_xxyyyy_y, g_0_y_xxyyyy_z, g_0_y_xyyyy_x, g_0_y_xyyyy_xx, g_0_y_xyyyy_xy, g_0_y_xyyyy_xz, g_0_y_xyyyy_y, g_0_y_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyy_x[k] = -g_0_y_xyyyy_x[k] * ab_x + g_0_y_xyyyy_xx[k];

                g_0_y_xxyyyy_y[k] = -g_0_y_xyyyy_y[k] * ab_x + g_0_y_xyyyy_xy[k];

                g_0_y_xxyyyy_z[k] = -g_0_y_xyyyy_z[k] * ab_x + g_0_y_xyyyy_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyz_x = cbuffer.data(ip_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxyyyz_y = cbuffer.data(ip_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxyyyz_z = cbuffer.data(ip_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyz_x, g_0_y_xxyyyz_y, g_0_y_xxyyyz_z, g_0_y_xyyyz_x, g_0_y_xyyyz_xx, g_0_y_xyyyz_xy, g_0_y_xyyyz_xz, g_0_y_xyyyz_y, g_0_y_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyz_x[k] = -g_0_y_xyyyz_x[k] * ab_x + g_0_y_xyyyz_xx[k];

                g_0_y_xxyyyz_y[k] = -g_0_y_xyyyz_y[k] * ab_x + g_0_y_xyyyz_xy[k];

                g_0_y_xxyyyz_z[k] = -g_0_y_xyyyz_z[k] * ab_x + g_0_y_xyyyz_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzz_x = cbuffer.data(ip_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxyyzz_y = cbuffer.data(ip_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxyyzz_z = cbuffer.data(ip_geom_01_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzz_x, g_0_y_xxyyzz_y, g_0_y_xxyyzz_z, g_0_y_xyyzz_x, g_0_y_xyyzz_xx, g_0_y_xyyzz_xy, g_0_y_xyyzz_xz, g_0_y_xyyzz_y, g_0_y_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzz_x[k] = -g_0_y_xyyzz_x[k] * ab_x + g_0_y_xyyzz_xx[k];

                g_0_y_xxyyzz_y[k] = -g_0_y_xyyzz_y[k] * ab_x + g_0_y_xyyzz_xy[k];

                g_0_y_xxyyzz_z[k] = -g_0_y_xyyzz_z[k] * ab_x + g_0_y_xyyzz_xz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzz_x = cbuffer.data(ip_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxyzzz_y = cbuffer.data(ip_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxyzzz_z = cbuffer.data(ip_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzz_x, g_0_y_xxyzzz_y, g_0_y_xxyzzz_z, g_0_y_xyzzz_x, g_0_y_xyzzz_xx, g_0_y_xyzzz_xy, g_0_y_xyzzz_xz, g_0_y_xyzzz_y, g_0_y_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzz_x[k] = -g_0_y_xyzzz_x[k] * ab_x + g_0_y_xyzzz_xx[k];

                g_0_y_xxyzzz_y[k] = -g_0_y_xyzzz_y[k] * ab_x + g_0_y_xyzzz_xy[k];

                g_0_y_xxyzzz_z[k] = -g_0_y_xyzzz_z[k] * ab_x + g_0_y_xyzzz_xz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzz_x = cbuffer.data(ip_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxzzzz_y = cbuffer.data(ip_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxzzzz_z = cbuffer.data(ip_geom_01_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzz_x, g_0_y_xxzzzz_y, g_0_y_xxzzzz_z, g_0_y_xzzzz_x, g_0_y_xzzzz_xx, g_0_y_xzzzz_xy, g_0_y_xzzzz_xz, g_0_y_xzzzz_y, g_0_y_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzz_x[k] = -g_0_y_xzzzz_x[k] * ab_x + g_0_y_xzzzz_xx[k];

                g_0_y_xxzzzz_y[k] = -g_0_y_xzzzz_y[k] * ab_x + g_0_y_xzzzz_xy[k];

                g_0_y_xxzzzz_z[k] = -g_0_y_xzzzz_z[k] * ab_x + g_0_y_xzzzz_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyy_x = cbuffer.data(ip_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyyyyy_y = cbuffer.data(ip_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyyyyy_z = cbuffer.data(ip_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyy_x, g_0_y_xyyyyy_y, g_0_y_xyyyyy_z, g_0_y_yyyyy_x, g_0_y_yyyyy_xx, g_0_y_yyyyy_xy, g_0_y_yyyyy_xz, g_0_y_yyyyy_y, g_0_y_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyy_x[k] = -g_0_y_yyyyy_x[k] * ab_x + g_0_y_yyyyy_xx[k];

                g_0_y_xyyyyy_y[k] = -g_0_y_yyyyy_y[k] * ab_x + g_0_y_yyyyy_xy[k];

                g_0_y_xyyyyy_z[k] = -g_0_y_yyyyy_z[k] * ab_x + g_0_y_yyyyy_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyz_x = cbuffer.data(ip_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyyyyz_y = cbuffer.data(ip_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyyyyz_z = cbuffer.data(ip_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyz_x, g_0_y_xyyyyz_y, g_0_y_xyyyyz_z, g_0_y_yyyyz_x, g_0_y_yyyyz_xx, g_0_y_yyyyz_xy, g_0_y_yyyyz_xz, g_0_y_yyyyz_y, g_0_y_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyz_x[k] = -g_0_y_yyyyz_x[k] * ab_x + g_0_y_yyyyz_xx[k];

                g_0_y_xyyyyz_y[k] = -g_0_y_yyyyz_y[k] * ab_x + g_0_y_yyyyz_xy[k];

                g_0_y_xyyyyz_z[k] = -g_0_y_yyyyz_z[k] * ab_x + g_0_y_yyyyz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzz_x = cbuffer.data(ip_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyyyzz_y = cbuffer.data(ip_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyyyzz_z = cbuffer.data(ip_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzz_x, g_0_y_xyyyzz_y, g_0_y_xyyyzz_z, g_0_y_yyyzz_x, g_0_y_yyyzz_xx, g_0_y_yyyzz_xy, g_0_y_yyyzz_xz, g_0_y_yyyzz_y, g_0_y_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzz_x[k] = -g_0_y_yyyzz_x[k] * ab_x + g_0_y_yyyzz_xx[k];

                g_0_y_xyyyzz_y[k] = -g_0_y_yyyzz_y[k] * ab_x + g_0_y_yyyzz_xy[k];

                g_0_y_xyyyzz_z[k] = -g_0_y_yyyzz_z[k] * ab_x + g_0_y_yyyzz_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzz_x = cbuffer.data(ip_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyyzzz_y = cbuffer.data(ip_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyyzzz_z = cbuffer.data(ip_geom_01_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzz_x, g_0_y_xyyzzz_y, g_0_y_xyyzzz_z, g_0_y_yyzzz_x, g_0_y_yyzzz_xx, g_0_y_yyzzz_xy, g_0_y_yyzzz_xz, g_0_y_yyzzz_y, g_0_y_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzz_x[k] = -g_0_y_yyzzz_x[k] * ab_x + g_0_y_yyzzz_xx[k];

                g_0_y_xyyzzz_y[k] = -g_0_y_yyzzz_y[k] * ab_x + g_0_y_yyzzz_xy[k];

                g_0_y_xyyzzz_z[k] = -g_0_y_yyzzz_z[k] * ab_x + g_0_y_yyzzz_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzz_x = cbuffer.data(ip_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyzzzz_y = cbuffer.data(ip_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyzzzz_z = cbuffer.data(ip_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzz_x, g_0_y_xyzzzz_y, g_0_y_xyzzzz_z, g_0_y_yzzzz_x, g_0_y_yzzzz_xx, g_0_y_yzzzz_xy, g_0_y_yzzzz_xz, g_0_y_yzzzz_y, g_0_y_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzz_x[k] = -g_0_y_yzzzz_x[k] * ab_x + g_0_y_yzzzz_xx[k];

                g_0_y_xyzzzz_y[k] = -g_0_y_yzzzz_y[k] * ab_x + g_0_y_yzzzz_xy[k];

                g_0_y_xyzzzz_z[k] = -g_0_y_yzzzz_z[k] * ab_x + g_0_y_yzzzz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzz_x = cbuffer.data(ip_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xzzzzz_y = cbuffer.data(ip_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xzzzzz_z = cbuffer.data(ip_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzz_x, g_0_y_xzzzzz_y, g_0_y_xzzzzz_z, g_0_y_zzzzz_x, g_0_y_zzzzz_xx, g_0_y_zzzzz_xy, g_0_y_zzzzz_xz, g_0_y_zzzzz_y, g_0_y_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzz_x[k] = -g_0_y_zzzzz_x[k] * ab_x + g_0_y_zzzzz_xx[k];

                g_0_y_xzzzzz_y[k] = -g_0_y_zzzzz_y[k] * ab_x + g_0_y_zzzzz_xy[k];

                g_0_y_xzzzzz_z[k] = -g_0_y_zzzzz_z[k] * ab_x + g_0_y_zzzzz_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyy_x = cbuffer.data(ip_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_yyyyyy_y = cbuffer.data(ip_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_yyyyyy_z = cbuffer.data(ip_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_x, g_0_y_yyyyy_xy, g_0_y_yyyyy_y, g_0_y_yyyyy_yy, g_0_y_yyyyy_yz, g_0_y_yyyyy_z, g_0_y_yyyyyy_x, g_0_y_yyyyyy_y, g_0_y_yyyyyy_z, g_yyyyy_x, g_yyyyy_y, g_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyy_x[k] = g_yyyyy_x[k] - g_0_y_yyyyy_x[k] * ab_y + g_0_y_yyyyy_xy[k];

                g_0_y_yyyyyy_y[k] = g_yyyyy_y[k] - g_0_y_yyyyy_y[k] * ab_y + g_0_y_yyyyy_yy[k];

                g_0_y_yyyyyy_z[k] = g_yyyyy_z[k] - g_0_y_yyyyy_z[k] * ab_y + g_0_y_yyyyy_yz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyz_x = cbuffer.data(ip_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yyyyyz_y = cbuffer.data(ip_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yyyyyz_z = cbuffer.data(ip_geom_01_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_x, g_0_y_yyyyy_xz, g_0_y_yyyyy_y, g_0_y_yyyyy_yz, g_0_y_yyyyy_z, g_0_y_yyyyy_zz, g_0_y_yyyyyz_x, g_0_y_yyyyyz_y, g_0_y_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyz_x[k] = -g_0_y_yyyyy_x[k] * ab_z + g_0_y_yyyyy_xz[k];

                g_0_y_yyyyyz_y[k] = -g_0_y_yyyyy_y[k] * ab_z + g_0_y_yyyyy_yz[k];

                g_0_y_yyyyyz_z[k] = -g_0_y_yyyyy_z[k] * ab_z + g_0_y_yyyyy_zz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzz_x = cbuffer.data(ip_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yyyyzz_y = cbuffer.data(ip_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yyyyzz_z = cbuffer.data(ip_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyz_x, g_0_y_yyyyz_xz, g_0_y_yyyyz_y, g_0_y_yyyyz_yz, g_0_y_yyyyz_z, g_0_y_yyyyz_zz, g_0_y_yyyyzz_x, g_0_y_yyyyzz_y, g_0_y_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzz_x[k] = -g_0_y_yyyyz_x[k] * ab_z + g_0_y_yyyyz_xz[k];

                g_0_y_yyyyzz_y[k] = -g_0_y_yyyyz_y[k] * ab_z + g_0_y_yyyyz_yz[k];

                g_0_y_yyyyzz_z[k] = -g_0_y_yyyyz_z[k] * ab_z + g_0_y_yyyyz_zz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzz_x = cbuffer.data(ip_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yyyzzz_y = cbuffer.data(ip_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yyyzzz_z = cbuffer.data(ip_geom_01_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzz_x, g_0_y_yyyzz_xz, g_0_y_yyyzz_y, g_0_y_yyyzz_yz, g_0_y_yyyzz_z, g_0_y_yyyzz_zz, g_0_y_yyyzzz_x, g_0_y_yyyzzz_y, g_0_y_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzz_x[k] = -g_0_y_yyyzz_x[k] * ab_z + g_0_y_yyyzz_xz[k];

                g_0_y_yyyzzz_y[k] = -g_0_y_yyyzz_y[k] * ab_z + g_0_y_yyyzz_yz[k];

                g_0_y_yyyzzz_z[k] = -g_0_y_yyyzz_z[k] * ab_z + g_0_y_yyyzz_zz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzz_x = cbuffer.data(ip_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyzzzz_y = cbuffer.data(ip_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyzzzz_z = cbuffer.data(ip_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzz_x, g_0_y_yyzzz_xz, g_0_y_yyzzz_y, g_0_y_yyzzz_yz, g_0_y_yyzzz_z, g_0_y_yyzzz_zz, g_0_y_yyzzzz_x, g_0_y_yyzzzz_y, g_0_y_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzz_x[k] = -g_0_y_yyzzz_x[k] * ab_z + g_0_y_yyzzz_xz[k];

                g_0_y_yyzzzz_y[k] = -g_0_y_yyzzz_y[k] * ab_z + g_0_y_yyzzz_yz[k];

                g_0_y_yyzzzz_z[k] = -g_0_y_yyzzz_z[k] * ab_z + g_0_y_yyzzz_zz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzz_x = cbuffer.data(ip_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yzzzzz_y = cbuffer.data(ip_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yzzzzz_z = cbuffer.data(ip_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzz_x, g_0_y_yzzzz_xz, g_0_y_yzzzz_y, g_0_y_yzzzz_yz, g_0_y_yzzzz_z, g_0_y_yzzzz_zz, g_0_y_yzzzzz_x, g_0_y_yzzzzz_y, g_0_y_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzz_x[k] = -g_0_y_yzzzz_x[k] * ab_z + g_0_y_yzzzz_xz[k];

                g_0_y_yzzzzz_y[k] = -g_0_y_yzzzz_y[k] * ab_z + g_0_y_yzzzz_yz[k];

                g_0_y_yzzzzz_z[k] = -g_0_y_yzzzz_z[k] * ab_z + g_0_y_yzzzz_zz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzz_x = cbuffer.data(ip_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_zzzzzz_y = cbuffer.data(ip_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_zzzzzz_z = cbuffer.data(ip_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzz_x, g_0_y_zzzzz_xz, g_0_y_zzzzz_y, g_0_y_zzzzz_yz, g_0_y_zzzzz_z, g_0_y_zzzzz_zz, g_0_y_zzzzzz_x, g_0_y_zzzzzz_y, g_0_y_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzz_x[k] = -g_0_y_zzzzz_x[k] * ab_z + g_0_y_zzzzz_xz[k];

                g_0_y_zzzzzz_y[k] = -g_0_y_zzzzz_y[k] * ab_z + g_0_y_zzzzz_yz[k];

                g_0_y_zzzzzz_z[k] = -g_0_y_zzzzz_z[k] * ab_z + g_0_y_zzzzz_zz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxx_x = cbuffer.data(ip_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_xxxxxx_y = cbuffer.data(ip_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_xxxxxx_z = cbuffer.data(ip_geom_01_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxx_x, g_0_z_xxxxx_xx, g_0_z_xxxxx_xy, g_0_z_xxxxx_xz, g_0_z_xxxxx_y, g_0_z_xxxxx_z, g_0_z_xxxxxx_x, g_0_z_xxxxxx_y, g_0_z_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxx_x[k] = -g_0_z_xxxxx_x[k] * ab_x + g_0_z_xxxxx_xx[k];

                g_0_z_xxxxxx_y[k] = -g_0_z_xxxxx_y[k] * ab_x + g_0_z_xxxxx_xy[k];

                g_0_z_xxxxxx_z[k] = -g_0_z_xxxxx_z[k] * ab_x + g_0_z_xxxxx_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxy_x = cbuffer.data(ip_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_xxxxxy_y = cbuffer.data(ip_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_xxxxxy_z = cbuffer.data(ip_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxy_x, g_0_z_xxxxxy_y, g_0_z_xxxxxy_z, g_0_z_xxxxy_x, g_0_z_xxxxy_xx, g_0_z_xxxxy_xy, g_0_z_xxxxy_xz, g_0_z_xxxxy_y, g_0_z_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxy_x[k] = -g_0_z_xxxxy_x[k] * ab_x + g_0_z_xxxxy_xx[k];

                g_0_z_xxxxxy_y[k] = -g_0_z_xxxxy_y[k] * ab_x + g_0_z_xxxxy_xy[k];

                g_0_z_xxxxxy_z[k] = -g_0_z_xxxxy_z[k] * ab_x + g_0_z_xxxxy_xz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxz_x = cbuffer.data(ip_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_xxxxxz_y = cbuffer.data(ip_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_xxxxxz_z = cbuffer.data(ip_geom_01_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxz_x, g_0_z_xxxxxz_y, g_0_z_xxxxxz_z, g_0_z_xxxxz_x, g_0_z_xxxxz_xx, g_0_z_xxxxz_xy, g_0_z_xxxxz_xz, g_0_z_xxxxz_y, g_0_z_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxz_x[k] = -g_0_z_xxxxz_x[k] * ab_x + g_0_z_xxxxz_xx[k];

                g_0_z_xxxxxz_y[k] = -g_0_z_xxxxz_y[k] * ab_x + g_0_z_xxxxz_xy[k];

                g_0_z_xxxxxz_z[k] = -g_0_z_xxxxz_z[k] * ab_x + g_0_z_xxxxz_xz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyy_x = cbuffer.data(ip_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_xxxxyy_y = cbuffer.data(ip_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_xxxxyy_z = cbuffer.data(ip_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyy_x, g_0_z_xxxxyy_y, g_0_z_xxxxyy_z, g_0_z_xxxyy_x, g_0_z_xxxyy_xx, g_0_z_xxxyy_xy, g_0_z_xxxyy_xz, g_0_z_xxxyy_y, g_0_z_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyy_x[k] = -g_0_z_xxxyy_x[k] * ab_x + g_0_z_xxxyy_xx[k];

                g_0_z_xxxxyy_y[k] = -g_0_z_xxxyy_y[k] * ab_x + g_0_z_xxxyy_xy[k];

                g_0_z_xxxxyy_z[k] = -g_0_z_xxxyy_z[k] * ab_x + g_0_z_xxxyy_xz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyz_x = cbuffer.data(ip_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xxxxyz_y = cbuffer.data(ip_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xxxxyz_z = cbuffer.data(ip_geom_01_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyz_x, g_0_z_xxxxyz_y, g_0_z_xxxxyz_z, g_0_z_xxxyz_x, g_0_z_xxxyz_xx, g_0_z_xxxyz_xy, g_0_z_xxxyz_xz, g_0_z_xxxyz_y, g_0_z_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyz_x[k] = -g_0_z_xxxyz_x[k] * ab_x + g_0_z_xxxyz_xx[k];

                g_0_z_xxxxyz_y[k] = -g_0_z_xxxyz_y[k] * ab_x + g_0_z_xxxyz_xy[k];

                g_0_z_xxxxyz_z[k] = -g_0_z_xxxyz_z[k] * ab_x + g_0_z_xxxyz_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzz_x = cbuffer.data(ip_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xxxxzz_y = cbuffer.data(ip_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xxxxzz_z = cbuffer.data(ip_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzz_x, g_0_z_xxxxzz_y, g_0_z_xxxxzz_z, g_0_z_xxxzz_x, g_0_z_xxxzz_xx, g_0_z_xxxzz_xy, g_0_z_xxxzz_xz, g_0_z_xxxzz_y, g_0_z_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzz_x[k] = -g_0_z_xxxzz_x[k] * ab_x + g_0_z_xxxzz_xx[k];

                g_0_z_xxxxzz_y[k] = -g_0_z_xxxzz_y[k] * ab_x + g_0_z_xxxzz_xy[k];

                g_0_z_xxxxzz_z[k] = -g_0_z_xxxzz_z[k] * ab_x + g_0_z_xxxzz_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyy_x = cbuffer.data(ip_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xxxyyy_y = cbuffer.data(ip_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xxxyyy_z = cbuffer.data(ip_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyy_x, g_0_z_xxxyyy_y, g_0_z_xxxyyy_z, g_0_z_xxyyy_x, g_0_z_xxyyy_xx, g_0_z_xxyyy_xy, g_0_z_xxyyy_xz, g_0_z_xxyyy_y, g_0_z_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyy_x[k] = -g_0_z_xxyyy_x[k] * ab_x + g_0_z_xxyyy_xx[k];

                g_0_z_xxxyyy_y[k] = -g_0_z_xxyyy_y[k] * ab_x + g_0_z_xxyyy_xy[k];

                g_0_z_xxxyyy_z[k] = -g_0_z_xxyyy_z[k] * ab_x + g_0_z_xxyyy_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyz_x = cbuffer.data(ip_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xxxyyz_y = cbuffer.data(ip_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xxxyyz_z = cbuffer.data(ip_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyz_x, g_0_z_xxxyyz_y, g_0_z_xxxyyz_z, g_0_z_xxyyz_x, g_0_z_xxyyz_xx, g_0_z_xxyyz_xy, g_0_z_xxyyz_xz, g_0_z_xxyyz_y, g_0_z_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyz_x[k] = -g_0_z_xxyyz_x[k] * ab_x + g_0_z_xxyyz_xx[k];

                g_0_z_xxxyyz_y[k] = -g_0_z_xxyyz_y[k] * ab_x + g_0_z_xxyyz_xy[k];

                g_0_z_xxxyyz_z[k] = -g_0_z_xxyyz_z[k] * ab_x + g_0_z_xxyyz_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzz_x = cbuffer.data(ip_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xxxyzz_y = cbuffer.data(ip_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xxxyzz_z = cbuffer.data(ip_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzz_x, g_0_z_xxxyzz_y, g_0_z_xxxyzz_z, g_0_z_xxyzz_x, g_0_z_xxyzz_xx, g_0_z_xxyzz_xy, g_0_z_xxyzz_xz, g_0_z_xxyzz_y, g_0_z_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzz_x[k] = -g_0_z_xxyzz_x[k] * ab_x + g_0_z_xxyzz_xx[k];

                g_0_z_xxxyzz_y[k] = -g_0_z_xxyzz_y[k] * ab_x + g_0_z_xxyzz_xy[k];

                g_0_z_xxxyzz_z[k] = -g_0_z_xxyzz_z[k] * ab_x + g_0_z_xxyzz_xz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzz_x = cbuffer.data(ip_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xxxzzz_y = cbuffer.data(ip_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xxxzzz_z = cbuffer.data(ip_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzz_x, g_0_z_xxxzzz_y, g_0_z_xxxzzz_z, g_0_z_xxzzz_x, g_0_z_xxzzz_xx, g_0_z_xxzzz_xy, g_0_z_xxzzz_xz, g_0_z_xxzzz_y, g_0_z_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzz_x[k] = -g_0_z_xxzzz_x[k] * ab_x + g_0_z_xxzzz_xx[k];

                g_0_z_xxxzzz_y[k] = -g_0_z_xxzzz_y[k] * ab_x + g_0_z_xxzzz_xy[k];

                g_0_z_xxxzzz_z[k] = -g_0_z_xxzzz_z[k] * ab_x + g_0_z_xxzzz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyy_x = cbuffer.data(ip_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xxyyyy_y = cbuffer.data(ip_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxyyyy_z = cbuffer.data(ip_geom_01_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyy_x, g_0_z_xxyyyy_y, g_0_z_xxyyyy_z, g_0_z_xyyyy_x, g_0_z_xyyyy_xx, g_0_z_xyyyy_xy, g_0_z_xyyyy_xz, g_0_z_xyyyy_y, g_0_z_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyy_x[k] = -g_0_z_xyyyy_x[k] * ab_x + g_0_z_xyyyy_xx[k];

                g_0_z_xxyyyy_y[k] = -g_0_z_xyyyy_y[k] * ab_x + g_0_z_xyyyy_xy[k];

                g_0_z_xxyyyy_z[k] = -g_0_z_xyyyy_z[k] * ab_x + g_0_z_xyyyy_xz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyz_x = cbuffer.data(ip_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxyyyz_y = cbuffer.data(ip_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxyyyz_z = cbuffer.data(ip_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyz_x, g_0_z_xxyyyz_y, g_0_z_xxyyyz_z, g_0_z_xyyyz_x, g_0_z_xyyyz_xx, g_0_z_xyyyz_xy, g_0_z_xyyyz_xz, g_0_z_xyyyz_y, g_0_z_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyz_x[k] = -g_0_z_xyyyz_x[k] * ab_x + g_0_z_xyyyz_xx[k];

                g_0_z_xxyyyz_y[k] = -g_0_z_xyyyz_y[k] * ab_x + g_0_z_xyyyz_xy[k];

                g_0_z_xxyyyz_z[k] = -g_0_z_xyyyz_z[k] * ab_x + g_0_z_xyyyz_xz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzz_x = cbuffer.data(ip_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxyyzz_y = cbuffer.data(ip_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxyyzz_z = cbuffer.data(ip_geom_01_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzz_x, g_0_z_xxyyzz_y, g_0_z_xxyyzz_z, g_0_z_xyyzz_x, g_0_z_xyyzz_xx, g_0_z_xyyzz_xy, g_0_z_xyyzz_xz, g_0_z_xyyzz_y, g_0_z_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzz_x[k] = -g_0_z_xyyzz_x[k] * ab_x + g_0_z_xyyzz_xx[k];

                g_0_z_xxyyzz_y[k] = -g_0_z_xyyzz_y[k] * ab_x + g_0_z_xyyzz_xy[k];

                g_0_z_xxyyzz_z[k] = -g_0_z_xyyzz_z[k] * ab_x + g_0_z_xyyzz_xz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzz_x = cbuffer.data(ip_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxyzzz_y = cbuffer.data(ip_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxyzzz_z = cbuffer.data(ip_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzz_x, g_0_z_xxyzzz_y, g_0_z_xxyzzz_z, g_0_z_xyzzz_x, g_0_z_xyzzz_xx, g_0_z_xyzzz_xy, g_0_z_xyzzz_xz, g_0_z_xyzzz_y, g_0_z_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzz_x[k] = -g_0_z_xyzzz_x[k] * ab_x + g_0_z_xyzzz_xx[k];

                g_0_z_xxyzzz_y[k] = -g_0_z_xyzzz_y[k] * ab_x + g_0_z_xyzzz_xy[k];

                g_0_z_xxyzzz_z[k] = -g_0_z_xyzzz_z[k] * ab_x + g_0_z_xyzzz_xz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzz_x = cbuffer.data(ip_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxzzzz_y = cbuffer.data(ip_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxzzzz_z = cbuffer.data(ip_geom_01_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzz_x, g_0_z_xxzzzz_y, g_0_z_xxzzzz_z, g_0_z_xzzzz_x, g_0_z_xzzzz_xx, g_0_z_xzzzz_xy, g_0_z_xzzzz_xz, g_0_z_xzzzz_y, g_0_z_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzz_x[k] = -g_0_z_xzzzz_x[k] * ab_x + g_0_z_xzzzz_xx[k];

                g_0_z_xxzzzz_y[k] = -g_0_z_xzzzz_y[k] * ab_x + g_0_z_xzzzz_xy[k];

                g_0_z_xxzzzz_z[k] = -g_0_z_xzzzz_z[k] * ab_x + g_0_z_xzzzz_xz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyy_x = cbuffer.data(ip_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xyyyyy_y = cbuffer.data(ip_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xyyyyy_z = cbuffer.data(ip_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyy_x, g_0_z_xyyyyy_y, g_0_z_xyyyyy_z, g_0_z_yyyyy_x, g_0_z_yyyyy_xx, g_0_z_yyyyy_xy, g_0_z_yyyyy_xz, g_0_z_yyyyy_y, g_0_z_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyy_x[k] = -g_0_z_yyyyy_x[k] * ab_x + g_0_z_yyyyy_xx[k];

                g_0_z_xyyyyy_y[k] = -g_0_z_yyyyy_y[k] * ab_x + g_0_z_yyyyy_xy[k];

                g_0_z_xyyyyy_z[k] = -g_0_z_yyyyy_z[k] * ab_x + g_0_z_yyyyy_xz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyz_x = cbuffer.data(ip_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xyyyyz_y = cbuffer.data(ip_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xyyyyz_z = cbuffer.data(ip_geom_01_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyz_x, g_0_z_xyyyyz_y, g_0_z_xyyyyz_z, g_0_z_yyyyz_x, g_0_z_yyyyz_xx, g_0_z_yyyyz_xy, g_0_z_yyyyz_xz, g_0_z_yyyyz_y, g_0_z_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyz_x[k] = -g_0_z_yyyyz_x[k] * ab_x + g_0_z_yyyyz_xx[k];

                g_0_z_xyyyyz_y[k] = -g_0_z_yyyyz_y[k] * ab_x + g_0_z_yyyyz_xy[k];

                g_0_z_xyyyyz_z[k] = -g_0_z_yyyyz_z[k] * ab_x + g_0_z_yyyyz_xz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzz_x = cbuffer.data(ip_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xyyyzz_y = cbuffer.data(ip_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xyyyzz_z = cbuffer.data(ip_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzz_x, g_0_z_xyyyzz_y, g_0_z_xyyyzz_z, g_0_z_yyyzz_x, g_0_z_yyyzz_xx, g_0_z_yyyzz_xy, g_0_z_yyyzz_xz, g_0_z_yyyzz_y, g_0_z_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzz_x[k] = -g_0_z_yyyzz_x[k] * ab_x + g_0_z_yyyzz_xx[k];

                g_0_z_xyyyzz_y[k] = -g_0_z_yyyzz_y[k] * ab_x + g_0_z_yyyzz_xy[k];

                g_0_z_xyyyzz_z[k] = -g_0_z_yyyzz_z[k] * ab_x + g_0_z_yyyzz_xz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzz_x = cbuffer.data(ip_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xyyzzz_y = cbuffer.data(ip_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xyyzzz_z = cbuffer.data(ip_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzz_x, g_0_z_xyyzzz_y, g_0_z_xyyzzz_z, g_0_z_yyzzz_x, g_0_z_yyzzz_xx, g_0_z_yyzzz_xy, g_0_z_yyzzz_xz, g_0_z_yyzzz_y, g_0_z_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzz_x[k] = -g_0_z_yyzzz_x[k] * ab_x + g_0_z_yyzzz_xx[k];

                g_0_z_xyyzzz_y[k] = -g_0_z_yyzzz_y[k] * ab_x + g_0_z_yyzzz_xy[k];

                g_0_z_xyyzzz_z[k] = -g_0_z_yyzzz_z[k] * ab_x + g_0_z_yyzzz_xz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzz_x = cbuffer.data(ip_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xyzzzz_y = cbuffer.data(ip_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xyzzzz_z = cbuffer.data(ip_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzz_x, g_0_z_xyzzzz_y, g_0_z_xyzzzz_z, g_0_z_yzzzz_x, g_0_z_yzzzz_xx, g_0_z_yzzzz_xy, g_0_z_yzzzz_xz, g_0_z_yzzzz_y, g_0_z_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzz_x[k] = -g_0_z_yzzzz_x[k] * ab_x + g_0_z_yzzzz_xx[k];

                g_0_z_xyzzzz_y[k] = -g_0_z_yzzzz_y[k] * ab_x + g_0_z_yzzzz_xy[k];

                g_0_z_xyzzzz_z[k] = -g_0_z_yzzzz_z[k] * ab_x + g_0_z_yzzzz_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzz_x = cbuffer.data(ip_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xzzzzz_y = cbuffer.data(ip_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xzzzzz_z = cbuffer.data(ip_geom_01_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzz_x, g_0_z_xzzzzz_y, g_0_z_xzzzzz_z, g_0_z_zzzzz_x, g_0_z_zzzzz_xx, g_0_z_zzzzz_xy, g_0_z_zzzzz_xz, g_0_z_zzzzz_y, g_0_z_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzz_x[k] = -g_0_z_zzzzz_x[k] * ab_x + g_0_z_zzzzz_xx[k];

                g_0_z_xzzzzz_y[k] = -g_0_z_zzzzz_y[k] * ab_x + g_0_z_zzzzz_xy[k];

                g_0_z_xzzzzz_z[k] = -g_0_z_zzzzz_z[k] * ab_x + g_0_z_zzzzz_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyy_x = cbuffer.data(ip_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_yyyyyy_y = cbuffer.data(ip_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_yyyyyy_z = cbuffer.data(ip_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyy_x, g_0_z_yyyyy_xy, g_0_z_yyyyy_y, g_0_z_yyyyy_yy, g_0_z_yyyyy_yz, g_0_z_yyyyy_z, g_0_z_yyyyyy_x, g_0_z_yyyyyy_y, g_0_z_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyy_x[k] = -g_0_z_yyyyy_x[k] * ab_y + g_0_z_yyyyy_xy[k];

                g_0_z_yyyyyy_y[k] = -g_0_z_yyyyy_y[k] * ab_y + g_0_z_yyyyy_yy[k];

                g_0_z_yyyyyy_z[k] = -g_0_z_yyyyy_z[k] * ab_y + g_0_z_yyyyy_yz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyz_x = cbuffer.data(ip_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_yyyyyz_y = cbuffer.data(ip_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_yyyyyz_z = cbuffer.data(ip_geom_01_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyz_x, g_0_z_yyyyyz_y, g_0_z_yyyyyz_z, g_0_z_yyyyz_x, g_0_z_yyyyz_xy, g_0_z_yyyyz_y, g_0_z_yyyyz_yy, g_0_z_yyyyz_yz, g_0_z_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyz_x[k] = -g_0_z_yyyyz_x[k] * ab_y + g_0_z_yyyyz_xy[k];

                g_0_z_yyyyyz_y[k] = -g_0_z_yyyyz_y[k] * ab_y + g_0_z_yyyyz_yy[k];

                g_0_z_yyyyyz_z[k] = -g_0_z_yyyyz_z[k] * ab_y + g_0_z_yyyyz_yz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzz_x = cbuffer.data(ip_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_yyyyzz_y = cbuffer.data(ip_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_yyyyzz_z = cbuffer.data(ip_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzz_x, g_0_z_yyyyzz_y, g_0_z_yyyyzz_z, g_0_z_yyyzz_x, g_0_z_yyyzz_xy, g_0_z_yyyzz_y, g_0_z_yyyzz_yy, g_0_z_yyyzz_yz, g_0_z_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzz_x[k] = -g_0_z_yyyzz_x[k] * ab_y + g_0_z_yyyzz_xy[k];

                g_0_z_yyyyzz_y[k] = -g_0_z_yyyzz_y[k] * ab_y + g_0_z_yyyzz_yy[k];

                g_0_z_yyyyzz_z[k] = -g_0_z_yyyzz_z[k] * ab_y + g_0_z_yyyzz_yz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzz_x = cbuffer.data(ip_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yyyzzz_y = cbuffer.data(ip_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yyyzzz_z = cbuffer.data(ip_geom_01_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzz_x, g_0_z_yyyzzz_y, g_0_z_yyyzzz_z, g_0_z_yyzzz_x, g_0_z_yyzzz_xy, g_0_z_yyzzz_y, g_0_z_yyzzz_yy, g_0_z_yyzzz_yz, g_0_z_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzz_x[k] = -g_0_z_yyzzz_x[k] * ab_y + g_0_z_yyzzz_xy[k];

                g_0_z_yyyzzz_y[k] = -g_0_z_yyzzz_y[k] * ab_y + g_0_z_yyzzz_yy[k];

                g_0_z_yyyzzz_z[k] = -g_0_z_yyzzz_z[k] * ab_y + g_0_z_yyzzz_yz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzz_x = cbuffer.data(ip_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yyzzzz_y = cbuffer.data(ip_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yyzzzz_z = cbuffer.data(ip_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzz_x, g_0_z_yyzzzz_y, g_0_z_yyzzzz_z, g_0_z_yzzzz_x, g_0_z_yzzzz_xy, g_0_z_yzzzz_y, g_0_z_yzzzz_yy, g_0_z_yzzzz_yz, g_0_z_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzz_x[k] = -g_0_z_yzzzz_x[k] * ab_y + g_0_z_yzzzz_xy[k];

                g_0_z_yyzzzz_y[k] = -g_0_z_yzzzz_y[k] * ab_y + g_0_z_yzzzz_yy[k];

                g_0_z_yyzzzz_z[k] = -g_0_z_yzzzz_z[k] * ab_y + g_0_z_yzzzz_yz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzz_x = cbuffer.data(ip_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yzzzzz_y = cbuffer.data(ip_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yzzzzz_z = cbuffer.data(ip_geom_01_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzz_x, g_0_z_yzzzzz_y, g_0_z_yzzzzz_z, g_0_z_zzzzz_x, g_0_z_zzzzz_xy, g_0_z_zzzzz_y, g_0_z_zzzzz_yy, g_0_z_zzzzz_yz, g_0_z_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzz_x[k] = -g_0_z_zzzzz_x[k] * ab_y + g_0_z_zzzzz_xy[k];

                g_0_z_yzzzzz_y[k] = -g_0_z_zzzzz_y[k] * ab_y + g_0_z_zzzzz_yy[k];

                g_0_z_yzzzzz_z[k] = -g_0_z_zzzzz_z[k] * ab_y + g_0_z_zzzzz_yz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzz_x = cbuffer.data(ip_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_zzzzzz_y = cbuffer.data(ip_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_zzzzzz_z = cbuffer.data(ip_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzz_x, g_0_z_zzzzz_xz, g_0_z_zzzzz_y, g_0_z_zzzzz_yz, g_0_z_zzzzz_z, g_0_z_zzzzz_zz, g_0_z_zzzzzz_x, g_0_z_zzzzzz_y, g_0_z_zzzzzz_z, g_zzzzz_x, g_zzzzz_y, g_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzz_x[k] = g_zzzzz_x[k] - g_0_z_zzzzz_x[k] * ab_z + g_0_z_zzzzz_xz[k];

                g_0_z_zzzzzz_y[k] = g_zzzzz_y[k] - g_0_z_zzzzz_y[k] * ab_z + g_0_z_zzzzz_yz[k];

                g_0_z_zzzzzz_z[k] = g_zzzzz_z[k] - g_0_z_zzzzz_z[k] * ab_z + g_0_z_zzzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

