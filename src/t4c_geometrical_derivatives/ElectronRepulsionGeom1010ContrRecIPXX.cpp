#include "ElectronRepulsionGeom1010ContrRecIPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_ipxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_ipxx,
                                              const size_t idx_geom_0010_hpxx,
                                              const size_t idx_geom_1010_hpxx,
                                              const size_t idx_geom_1010_hdxx,
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

            const auto hp_geom_0010_off = idx_geom_0010_hpxx + i * dcomps + j;

            auto g_0_0_x_0_xxxxx_x = cbuffer.data(hp_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_y = cbuffer.data(hp_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxx_z = cbuffer.data(hp_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_x = cbuffer.data(hp_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_y = cbuffer.data(hp_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxy_z = cbuffer.data(hp_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_x = cbuffer.data(hp_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_y = cbuffer.data(hp_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xxxxz_z = cbuffer.data(hp_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_x = cbuffer.data(hp_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_y = cbuffer.data(hp_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyy_z = cbuffer.data(hp_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_x = cbuffer.data(hp_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_y = cbuffer.data(hp_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xxxyz_z = cbuffer.data(hp_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_x = cbuffer.data(hp_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_y = cbuffer.data(hp_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_xxxzz_z = cbuffer.data(hp_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_x = cbuffer.data(hp_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_y = cbuffer.data(hp_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyy_z = cbuffer.data(hp_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_x = cbuffer.data(hp_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_y = cbuffer.data(hp_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_xxyyz_z = cbuffer.data(hp_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_x = cbuffer.data(hp_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_y = cbuffer.data(hp_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_xxyzz_z = cbuffer.data(hp_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_x = cbuffer.data(hp_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_y = cbuffer.data(hp_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_x_0_xxzzz_z = cbuffer.data(hp_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_x = cbuffer.data(hp_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_y = cbuffer.data(hp_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyy_z = cbuffer.data(hp_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_x = cbuffer.data(hp_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_y = cbuffer.data(hp_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_x_0_xyyyz_z = cbuffer.data(hp_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_x = cbuffer.data(hp_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_y = cbuffer.data(hp_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_x_0_xyyzz_z = cbuffer.data(hp_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_x = cbuffer.data(hp_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_y = cbuffer.data(hp_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_x_0_xyzzz_z = cbuffer.data(hp_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_x = cbuffer.data(hp_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_y = cbuffer.data(hp_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_x_0_xzzzz_z = cbuffer.data(hp_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_x = cbuffer.data(hp_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_y = cbuffer.data(hp_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyy_z = cbuffer.data(hp_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_x = cbuffer.data(hp_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_y = cbuffer.data(hp_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_x_0_yyyyz_z = cbuffer.data(hp_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_x = cbuffer.data(hp_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_y = cbuffer.data(hp_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_x_0_yyyzz_z = cbuffer.data(hp_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_x = cbuffer.data(hp_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_y = cbuffer.data(hp_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_x_0_yyzzz_z = cbuffer.data(hp_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_x = cbuffer.data(hp_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_y = cbuffer.data(hp_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_x_0_yzzzz_z = cbuffer.data(hp_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_x = cbuffer.data(hp_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_y = cbuffer.data(hp_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_x_0_zzzzz_z = cbuffer.data(hp_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_x = cbuffer.data(hp_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_y = cbuffer.data(hp_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxx_z = cbuffer.data(hp_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_x = cbuffer.data(hp_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_y = cbuffer.data(hp_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxy_z = cbuffer.data(hp_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_x = cbuffer.data(hp_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_y = cbuffer.data(hp_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_y_0_xxxxz_z = cbuffer.data(hp_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_x = cbuffer.data(hp_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_y = cbuffer.data(hp_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyy_z = cbuffer.data(hp_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_x = cbuffer.data(hp_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_y = cbuffer.data(hp_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_y_0_xxxyz_z = cbuffer.data(hp_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_x = cbuffer.data(hp_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_y = cbuffer.data(hp_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_y_0_xxxzz_z = cbuffer.data(hp_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_x = cbuffer.data(hp_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_y = cbuffer.data(hp_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyy_z = cbuffer.data(hp_geom_0010_off + 83 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_x = cbuffer.data(hp_geom_0010_off + 84 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_y = cbuffer.data(hp_geom_0010_off + 85 * ccomps * dcomps);

            auto g_0_0_y_0_xxyyz_z = cbuffer.data(hp_geom_0010_off + 86 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_x = cbuffer.data(hp_geom_0010_off + 87 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_y = cbuffer.data(hp_geom_0010_off + 88 * ccomps * dcomps);

            auto g_0_0_y_0_xxyzz_z = cbuffer.data(hp_geom_0010_off + 89 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_x = cbuffer.data(hp_geom_0010_off + 90 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_y = cbuffer.data(hp_geom_0010_off + 91 * ccomps * dcomps);

            auto g_0_0_y_0_xxzzz_z = cbuffer.data(hp_geom_0010_off + 92 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_x = cbuffer.data(hp_geom_0010_off + 93 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_y = cbuffer.data(hp_geom_0010_off + 94 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyy_z = cbuffer.data(hp_geom_0010_off + 95 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_x = cbuffer.data(hp_geom_0010_off + 96 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_y = cbuffer.data(hp_geom_0010_off + 97 * ccomps * dcomps);

            auto g_0_0_y_0_xyyyz_z = cbuffer.data(hp_geom_0010_off + 98 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_x = cbuffer.data(hp_geom_0010_off + 99 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_y = cbuffer.data(hp_geom_0010_off + 100 * ccomps * dcomps);

            auto g_0_0_y_0_xyyzz_z = cbuffer.data(hp_geom_0010_off + 101 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_x = cbuffer.data(hp_geom_0010_off + 102 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_y = cbuffer.data(hp_geom_0010_off + 103 * ccomps * dcomps);

            auto g_0_0_y_0_xyzzz_z = cbuffer.data(hp_geom_0010_off + 104 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_x = cbuffer.data(hp_geom_0010_off + 105 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_y = cbuffer.data(hp_geom_0010_off + 106 * ccomps * dcomps);

            auto g_0_0_y_0_xzzzz_z = cbuffer.data(hp_geom_0010_off + 107 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_x = cbuffer.data(hp_geom_0010_off + 108 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_y = cbuffer.data(hp_geom_0010_off + 109 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyy_z = cbuffer.data(hp_geom_0010_off + 110 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_x = cbuffer.data(hp_geom_0010_off + 111 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_y = cbuffer.data(hp_geom_0010_off + 112 * ccomps * dcomps);

            auto g_0_0_y_0_yyyyz_z = cbuffer.data(hp_geom_0010_off + 113 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_x = cbuffer.data(hp_geom_0010_off + 114 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_y = cbuffer.data(hp_geom_0010_off + 115 * ccomps * dcomps);

            auto g_0_0_y_0_yyyzz_z = cbuffer.data(hp_geom_0010_off + 116 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_x = cbuffer.data(hp_geom_0010_off + 117 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_y = cbuffer.data(hp_geom_0010_off + 118 * ccomps * dcomps);

            auto g_0_0_y_0_yyzzz_z = cbuffer.data(hp_geom_0010_off + 119 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_x = cbuffer.data(hp_geom_0010_off + 120 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_y = cbuffer.data(hp_geom_0010_off + 121 * ccomps * dcomps);

            auto g_0_0_y_0_yzzzz_z = cbuffer.data(hp_geom_0010_off + 122 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_x = cbuffer.data(hp_geom_0010_off + 123 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_y = cbuffer.data(hp_geom_0010_off + 124 * ccomps * dcomps);

            auto g_0_0_y_0_zzzzz_z = cbuffer.data(hp_geom_0010_off + 125 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_x = cbuffer.data(hp_geom_0010_off + 126 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_y = cbuffer.data(hp_geom_0010_off + 127 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxx_z = cbuffer.data(hp_geom_0010_off + 128 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_x = cbuffer.data(hp_geom_0010_off + 129 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_y = cbuffer.data(hp_geom_0010_off + 130 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxy_z = cbuffer.data(hp_geom_0010_off + 131 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_x = cbuffer.data(hp_geom_0010_off + 132 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_y = cbuffer.data(hp_geom_0010_off + 133 * ccomps * dcomps);

            auto g_0_0_z_0_xxxxz_z = cbuffer.data(hp_geom_0010_off + 134 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_x = cbuffer.data(hp_geom_0010_off + 135 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_y = cbuffer.data(hp_geom_0010_off + 136 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyy_z = cbuffer.data(hp_geom_0010_off + 137 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_x = cbuffer.data(hp_geom_0010_off + 138 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_y = cbuffer.data(hp_geom_0010_off + 139 * ccomps * dcomps);

            auto g_0_0_z_0_xxxyz_z = cbuffer.data(hp_geom_0010_off + 140 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_x = cbuffer.data(hp_geom_0010_off + 141 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_y = cbuffer.data(hp_geom_0010_off + 142 * ccomps * dcomps);

            auto g_0_0_z_0_xxxzz_z = cbuffer.data(hp_geom_0010_off + 143 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_x = cbuffer.data(hp_geom_0010_off + 144 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_y = cbuffer.data(hp_geom_0010_off + 145 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyy_z = cbuffer.data(hp_geom_0010_off + 146 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_x = cbuffer.data(hp_geom_0010_off + 147 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_y = cbuffer.data(hp_geom_0010_off + 148 * ccomps * dcomps);

            auto g_0_0_z_0_xxyyz_z = cbuffer.data(hp_geom_0010_off + 149 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_x = cbuffer.data(hp_geom_0010_off + 150 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_y = cbuffer.data(hp_geom_0010_off + 151 * ccomps * dcomps);

            auto g_0_0_z_0_xxyzz_z = cbuffer.data(hp_geom_0010_off + 152 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_x = cbuffer.data(hp_geom_0010_off + 153 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_y = cbuffer.data(hp_geom_0010_off + 154 * ccomps * dcomps);

            auto g_0_0_z_0_xxzzz_z = cbuffer.data(hp_geom_0010_off + 155 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_x = cbuffer.data(hp_geom_0010_off + 156 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_y = cbuffer.data(hp_geom_0010_off + 157 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyy_z = cbuffer.data(hp_geom_0010_off + 158 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_x = cbuffer.data(hp_geom_0010_off + 159 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_y = cbuffer.data(hp_geom_0010_off + 160 * ccomps * dcomps);

            auto g_0_0_z_0_xyyyz_z = cbuffer.data(hp_geom_0010_off + 161 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_x = cbuffer.data(hp_geom_0010_off + 162 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_y = cbuffer.data(hp_geom_0010_off + 163 * ccomps * dcomps);

            auto g_0_0_z_0_xyyzz_z = cbuffer.data(hp_geom_0010_off + 164 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_x = cbuffer.data(hp_geom_0010_off + 165 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_y = cbuffer.data(hp_geom_0010_off + 166 * ccomps * dcomps);

            auto g_0_0_z_0_xyzzz_z = cbuffer.data(hp_geom_0010_off + 167 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_x = cbuffer.data(hp_geom_0010_off + 168 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_y = cbuffer.data(hp_geom_0010_off + 169 * ccomps * dcomps);

            auto g_0_0_z_0_xzzzz_z = cbuffer.data(hp_geom_0010_off + 170 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_x = cbuffer.data(hp_geom_0010_off + 171 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_y = cbuffer.data(hp_geom_0010_off + 172 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyy_z = cbuffer.data(hp_geom_0010_off + 173 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_x = cbuffer.data(hp_geom_0010_off + 174 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_y = cbuffer.data(hp_geom_0010_off + 175 * ccomps * dcomps);

            auto g_0_0_z_0_yyyyz_z = cbuffer.data(hp_geom_0010_off + 176 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_x = cbuffer.data(hp_geom_0010_off + 177 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_y = cbuffer.data(hp_geom_0010_off + 178 * ccomps * dcomps);

            auto g_0_0_z_0_yyyzz_z = cbuffer.data(hp_geom_0010_off + 179 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_x = cbuffer.data(hp_geom_0010_off + 180 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_y = cbuffer.data(hp_geom_0010_off + 181 * ccomps * dcomps);

            auto g_0_0_z_0_yyzzz_z = cbuffer.data(hp_geom_0010_off + 182 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_x = cbuffer.data(hp_geom_0010_off + 183 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_y = cbuffer.data(hp_geom_0010_off + 184 * ccomps * dcomps);

            auto g_0_0_z_0_yzzzz_z = cbuffer.data(hp_geom_0010_off + 185 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_x = cbuffer.data(hp_geom_0010_off + 186 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_y = cbuffer.data(hp_geom_0010_off + 187 * ccomps * dcomps);

            auto g_0_0_z_0_zzzzz_z = cbuffer.data(hp_geom_0010_off + 188 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_1010_off = idx_geom_1010_hdxx + i * dcomps + j;

            auto g_x_0_x_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 377 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 378 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 379 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 380 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 381 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 382 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 383 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 384 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 385 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 386 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 387 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 388 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 389 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 390 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 391 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 392 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 393 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 394 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 395 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 396 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 397 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 398 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 399 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 400 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 401 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 402 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 403 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 404 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 405 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 406 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 407 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 408 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 409 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 410 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 411 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 412 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 413 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 414 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 415 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 416 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 417 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 418 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 419 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 420 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 421 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 422 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 423 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 424 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 425 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 426 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 427 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 428 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 429 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 430 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 431 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 432 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 433 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 434 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 435 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 436 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 437 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 438 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 439 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 440 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 441 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 442 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 443 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 444 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 445 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 446 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 447 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 448 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 449 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 450 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 451 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 452 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 453 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 454 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 455 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 456 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 457 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 458 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 459 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 460 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 461 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 462 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 463 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 464 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 465 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 466 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 467 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 468 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 469 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 470 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 471 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 472 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 473 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 474 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 475 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 476 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 477 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 478 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 479 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 480 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 481 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 482 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 483 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 484 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 485 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 486 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 487 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 488 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 489 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 490 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 491 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 492 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 493 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 494 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 495 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 496 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 497 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 498 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 499 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 500 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 501 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 502 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 503 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 504 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 505 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 506 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 507 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 508 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 509 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 510 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 511 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 512 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 513 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 514 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 515 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 516 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 517 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 518 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 519 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 520 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 521 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 522 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 523 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 524 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 525 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 526 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 527 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 528 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 529 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 530 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 531 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 532 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 533 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 534 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 535 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 536 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 537 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 538 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 539 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 540 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 541 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 542 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 543 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 544 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 545 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 546 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 547 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 548 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 549 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 550 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 551 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 552 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 553 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 554 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 555 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 556 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 557 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 558 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 559 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 560 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 561 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 562 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 563 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 564 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 565 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 566 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 567 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 568 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 569 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 570 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 571 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 572 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 573 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 574 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 575 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 576 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 577 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 578 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 579 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 580 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 581 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 582 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 583 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 584 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 585 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 586 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 587 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 588 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 589 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 590 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 591 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 592 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 593 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 594 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 595 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 596 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 597 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 598 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 599 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 600 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 601 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 602 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 603 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 604 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 605 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 606 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 607 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 608 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 609 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 610 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 611 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 612 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 613 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 614 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 615 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 616 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 617 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 618 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 619 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 620 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 621 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 622 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 623 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 624 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 625 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 626 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 627 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 628 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 629 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 671 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 755 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 756 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 757 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 758 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 759 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 760 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 761 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 762 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 763 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 764 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 765 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 766 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 767 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 768 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 769 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 770 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 771 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 772 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 773 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 774 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 775 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 776 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 777 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 778 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 779 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 780 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 781 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 782 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 783 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 784 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 785 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 786 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 787 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 788 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 789 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 790 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 791 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 792 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 793 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 794 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 795 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 796 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 797 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 798 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 799 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 800 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 801 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 802 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 803 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 804 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 805 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 806 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 807 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 808 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 809 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 810 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 811 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 812 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 813 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 814 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 815 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 816 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 817 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 818 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 819 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 820 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 821 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 822 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 823 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 824 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 825 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 826 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 827 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 828 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 829 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 830 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 831 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 832 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 833 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 834 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 835 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 836 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 837 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 838 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 839 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 840 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 841 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 842 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 843 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 844 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 845 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 846 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 847 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 848 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 849 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 850 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 851 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 852 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 853 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 854 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 855 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 856 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 857 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 858 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 859 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 860 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 861 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 862 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 863 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 864 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 865 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 866 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 867 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 868 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 869 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 870 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 871 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 872 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 873 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 874 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 875 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 876 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 877 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 878 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 879 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 880 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 881 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 882 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 883 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 884 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 885 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 886 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 887 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 888 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 889 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 890 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 891 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 892 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 893 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 894 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 895 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 896 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 897 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 898 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 899 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 900 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 901 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 902 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 903 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 904 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 905 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 906 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 907 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 908 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 909 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 910 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 911 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 912 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 913 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 914 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 915 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 916 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 917 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 918 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 919 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 920 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 921 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 922 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 923 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 924 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 925 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 926 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 927 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 928 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 929 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 930 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 931 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 932 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 933 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 934 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 935 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 936 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 937 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 938 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 939 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 940 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 941 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 942 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 943 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 944 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 945 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 946 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 947 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 948 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 949 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 950 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 951 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 952 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 953 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 954 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 955 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 956 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 957 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 958 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 959 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 960 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 961 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 962 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 963 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 964 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 965 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 966 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 967 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 968 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 969 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 970 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 971 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 972 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 973 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 974 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 975 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 976 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 977 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 978 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 979 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 980 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 981 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 982 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 983 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 984 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 985 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 986 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 987 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 988 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 989 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 990 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 991 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 992 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 993 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 994 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 995 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 996 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 997 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 998 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 999 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xx = cbuffer.data(hd_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xy = cbuffer.data(hd_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_xz = cbuffer.data(hd_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_yy = cbuffer.data(hd_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_yz = cbuffer.data(hd_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxx_zz = cbuffer.data(hd_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xx = cbuffer.data(hd_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xy = cbuffer.data(hd_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_xz = cbuffer.data(hd_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_yy = cbuffer.data(hd_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_yz = cbuffer.data(hd_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxy_zz = cbuffer.data(hd_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xx = cbuffer.data(hd_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xy = cbuffer.data(hd_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_xz = cbuffer.data(hd_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_yy = cbuffer.data(hd_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_yz = cbuffer.data(hd_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxz_zz = cbuffer.data(hd_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xx = cbuffer.data(hd_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xy = cbuffer.data(hd_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_xz = cbuffer.data(hd_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_yy = cbuffer.data(hd_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_yz = cbuffer.data(hd_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyy_zz = cbuffer.data(hd_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xx = cbuffer.data(hd_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xy = cbuffer.data(hd_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_xz = cbuffer.data(hd_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_yy = cbuffer.data(hd_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_yz = cbuffer.data(hd_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyz_zz = cbuffer.data(hd_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xx = cbuffer.data(hd_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xy = cbuffer.data(hd_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_xz = cbuffer.data(hd_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_yy = cbuffer.data(hd_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_yz = cbuffer.data(hd_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzz_zz = cbuffer.data(hd_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xx = cbuffer.data(hd_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xy = cbuffer.data(hd_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_xz = cbuffer.data(hd_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_yy = cbuffer.data(hd_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_yz = cbuffer.data(hd_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyy_zz = cbuffer.data(hd_geom_1010_off + 1049 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xx = cbuffer.data(hd_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xy = cbuffer.data(hd_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_xz = cbuffer.data(hd_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_yy = cbuffer.data(hd_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_yz = cbuffer.data(hd_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyz_zz = cbuffer.data(hd_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xx = cbuffer.data(hd_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xy = cbuffer.data(hd_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_xz = cbuffer.data(hd_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_yy = cbuffer.data(hd_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_yz = cbuffer.data(hd_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzz_zz = cbuffer.data(hd_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xx = cbuffer.data(hd_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xy = cbuffer.data(hd_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_xz = cbuffer.data(hd_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_yy = cbuffer.data(hd_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_yz = cbuffer.data(hd_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzz_zz = cbuffer.data(hd_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xx = cbuffer.data(hd_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xy = cbuffer.data(hd_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_xz = cbuffer.data(hd_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_yy = cbuffer.data(hd_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_yz = cbuffer.data(hd_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyy_zz = cbuffer.data(hd_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xx = cbuffer.data(hd_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xy = cbuffer.data(hd_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_xz = cbuffer.data(hd_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_yy = cbuffer.data(hd_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_yz = cbuffer.data(hd_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyz_zz = cbuffer.data(hd_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xx = cbuffer.data(hd_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xy = cbuffer.data(hd_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_xz = cbuffer.data(hd_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_yy = cbuffer.data(hd_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_yz = cbuffer.data(hd_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzz_zz = cbuffer.data(hd_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xx = cbuffer.data(hd_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xy = cbuffer.data(hd_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_xz = cbuffer.data(hd_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_yy = cbuffer.data(hd_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_yz = cbuffer.data(hd_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzz_zz = cbuffer.data(hd_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xx = cbuffer.data(hd_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xy = cbuffer.data(hd_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_xz = cbuffer.data(hd_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_yy = cbuffer.data(hd_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_yz = cbuffer.data(hd_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzz_zz = cbuffer.data(hd_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xx = cbuffer.data(hd_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xy = cbuffer.data(hd_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_xz = cbuffer.data(hd_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_yy = cbuffer.data(hd_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_yz = cbuffer.data(hd_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyy_zz = cbuffer.data(hd_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xx = cbuffer.data(hd_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xy = cbuffer.data(hd_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_xz = cbuffer.data(hd_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_yy = cbuffer.data(hd_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_yz = cbuffer.data(hd_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyz_zz = cbuffer.data(hd_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xx = cbuffer.data(hd_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xy = cbuffer.data(hd_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_xz = cbuffer.data(hd_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_yy = cbuffer.data(hd_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_yz = cbuffer.data(hd_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzz_zz = cbuffer.data(hd_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xx = cbuffer.data(hd_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xy = cbuffer.data(hd_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_xz = cbuffer.data(hd_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_yy = cbuffer.data(hd_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_yz = cbuffer.data(hd_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzz_zz = cbuffer.data(hd_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xx = cbuffer.data(hd_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xy = cbuffer.data(hd_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_xz = cbuffer.data(hd_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_yy = cbuffer.data(hd_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_yz = cbuffer.data(hd_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzz_zz = cbuffer.data(hd_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xx = cbuffer.data(hd_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xy = cbuffer.data(hd_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_xz = cbuffer.data(hd_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_yy = cbuffer.data(hd_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_yz = cbuffer.data(hd_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzz_zz = cbuffer.data(hd_geom_1010_off + 1133 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ipxx

            const auto ip_geom_1010_off = idx_geom_1010_ipxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxxxx_x, g_0_0_x_0_xxxxx_y, g_0_0_x_0_xxxxx_z, g_x_0_x_0_xxxxx_x, g_x_0_x_0_xxxxx_xx, g_x_0_x_0_xxxxx_xy, g_x_0_x_0_xxxxx_xz, g_x_0_x_0_xxxxx_y, g_x_0_x_0_xxxxx_z, g_x_0_x_0_xxxxxx_x, g_x_0_x_0_xxxxxx_y, g_x_0_x_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxx_x[k] = -g_0_0_x_0_xxxxx_x[k] - g_x_0_x_0_xxxxx_x[k] * ab_x + g_x_0_x_0_xxxxx_xx[k];

                g_x_0_x_0_xxxxxx_y[k] = -g_0_0_x_0_xxxxx_y[k] - g_x_0_x_0_xxxxx_y[k] * ab_x + g_x_0_x_0_xxxxx_xy[k];

                g_x_0_x_0_xxxxxx_z[k] = -g_0_0_x_0_xxxxx_z[k] - g_x_0_x_0_xxxxx_z[k] * ab_x + g_x_0_x_0_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_x, g_x_0_x_0_xxxxx_xy, g_x_0_x_0_xxxxx_y, g_x_0_x_0_xxxxx_yy, g_x_0_x_0_xxxxx_yz, g_x_0_x_0_xxxxx_z, g_x_0_x_0_xxxxxy_x, g_x_0_x_0_xxxxxy_y, g_x_0_x_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxy_x[k] = -g_x_0_x_0_xxxxx_x[k] * ab_y + g_x_0_x_0_xxxxx_xy[k];

                g_x_0_x_0_xxxxxy_y[k] = -g_x_0_x_0_xxxxx_y[k] * ab_y + g_x_0_x_0_xxxxx_yy[k];

                g_x_0_x_0_xxxxxy_z[k] = -g_x_0_x_0_xxxxx_z[k] * ab_y + g_x_0_x_0_xxxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxx_x, g_x_0_x_0_xxxxx_xz, g_x_0_x_0_xxxxx_y, g_x_0_x_0_xxxxx_yz, g_x_0_x_0_xxxxx_z, g_x_0_x_0_xxxxx_zz, g_x_0_x_0_xxxxxz_x, g_x_0_x_0_xxxxxz_y, g_x_0_x_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxxz_x[k] = -g_x_0_x_0_xxxxx_x[k] * ab_z + g_x_0_x_0_xxxxx_xz[k];

                g_x_0_x_0_xxxxxz_y[k] = -g_x_0_x_0_xxxxx_y[k] * ab_z + g_x_0_x_0_xxxxx_yz[k];

                g_x_0_x_0_xxxxxz_z[k] = -g_x_0_x_0_xxxxx_z[k] * ab_z + g_x_0_x_0_xxxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxy_x, g_x_0_x_0_xxxxy_xy, g_x_0_x_0_xxxxy_y, g_x_0_x_0_xxxxy_yy, g_x_0_x_0_xxxxy_yz, g_x_0_x_0_xxxxy_z, g_x_0_x_0_xxxxyy_x, g_x_0_x_0_xxxxyy_y, g_x_0_x_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyy_x[k] = -g_x_0_x_0_xxxxy_x[k] * ab_y + g_x_0_x_0_xxxxy_xy[k];

                g_x_0_x_0_xxxxyy_y[k] = -g_x_0_x_0_xxxxy_y[k] * ab_y + g_x_0_x_0_xxxxy_yy[k];

                g_x_0_x_0_xxxxyy_z[k] = -g_x_0_x_0_xxxxy_z[k] * ab_y + g_x_0_x_0_xxxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxyz_x, g_x_0_x_0_xxxxyz_y, g_x_0_x_0_xxxxyz_z, g_x_0_x_0_xxxxz_x, g_x_0_x_0_xxxxz_xy, g_x_0_x_0_xxxxz_y, g_x_0_x_0_xxxxz_yy, g_x_0_x_0_xxxxz_yz, g_x_0_x_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxyz_x[k] = -g_x_0_x_0_xxxxz_x[k] * ab_y + g_x_0_x_0_xxxxz_xy[k];

                g_x_0_x_0_xxxxyz_y[k] = -g_x_0_x_0_xxxxz_y[k] * ab_y + g_x_0_x_0_xxxxz_yy[k];

                g_x_0_x_0_xxxxyz_z[k] = -g_x_0_x_0_xxxxz_z[k] * ab_y + g_x_0_x_0_xxxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxxz_x, g_x_0_x_0_xxxxz_xz, g_x_0_x_0_xxxxz_y, g_x_0_x_0_xxxxz_yz, g_x_0_x_0_xxxxz_z, g_x_0_x_0_xxxxz_zz, g_x_0_x_0_xxxxzz_x, g_x_0_x_0_xxxxzz_y, g_x_0_x_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxxzz_x[k] = -g_x_0_x_0_xxxxz_x[k] * ab_z + g_x_0_x_0_xxxxz_xz[k];

                g_x_0_x_0_xxxxzz_y[k] = -g_x_0_x_0_xxxxz_y[k] * ab_z + g_x_0_x_0_xxxxz_yz[k];

                g_x_0_x_0_xxxxzz_z[k] = -g_x_0_x_0_xxxxz_z[k] * ab_z + g_x_0_x_0_xxxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyy_x, g_x_0_x_0_xxxyy_xy, g_x_0_x_0_xxxyy_y, g_x_0_x_0_xxxyy_yy, g_x_0_x_0_xxxyy_yz, g_x_0_x_0_xxxyy_z, g_x_0_x_0_xxxyyy_x, g_x_0_x_0_xxxyyy_y, g_x_0_x_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyy_x[k] = -g_x_0_x_0_xxxyy_x[k] * ab_y + g_x_0_x_0_xxxyy_xy[k];

                g_x_0_x_0_xxxyyy_y[k] = -g_x_0_x_0_xxxyy_y[k] * ab_y + g_x_0_x_0_xxxyy_yy[k];

                g_x_0_x_0_xxxyyy_z[k] = -g_x_0_x_0_xxxyy_z[k] * ab_y + g_x_0_x_0_xxxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyyz_x, g_x_0_x_0_xxxyyz_y, g_x_0_x_0_xxxyyz_z, g_x_0_x_0_xxxyz_x, g_x_0_x_0_xxxyz_xy, g_x_0_x_0_xxxyz_y, g_x_0_x_0_xxxyz_yy, g_x_0_x_0_xxxyz_yz, g_x_0_x_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyyz_x[k] = -g_x_0_x_0_xxxyz_x[k] * ab_y + g_x_0_x_0_xxxyz_xy[k];

                g_x_0_x_0_xxxyyz_y[k] = -g_x_0_x_0_xxxyz_y[k] * ab_y + g_x_0_x_0_xxxyz_yy[k];

                g_x_0_x_0_xxxyyz_z[k] = -g_x_0_x_0_xxxyz_z[k] * ab_y + g_x_0_x_0_xxxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxyzz_x, g_x_0_x_0_xxxyzz_y, g_x_0_x_0_xxxyzz_z, g_x_0_x_0_xxxzz_x, g_x_0_x_0_xxxzz_xy, g_x_0_x_0_xxxzz_y, g_x_0_x_0_xxxzz_yy, g_x_0_x_0_xxxzz_yz, g_x_0_x_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxyzz_x[k] = -g_x_0_x_0_xxxzz_x[k] * ab_y + g_x_0_x_0_xxxzz_xy[k];

                g_x_0_x_0_xxxyzz_y[k] = -g_x_0_x_0_xxxzz_y[k] * ab_y + g_x_0_x_0_xxxzz_yy[k];

                g_x_0_x_0_xxxyzz_z[k] = -g_x_0_x_0_xxxzz_z[k] * ab_y + g_x_0_x_0_xxxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxxzz_x, g_x_0_x_0_xxxzz_xz, g_x_0_x_0_xxxzz_y, g_x_0_x_0_xxxzz_yz, g_x_0_x_0_xxxzz_z, g_x_0_x_0_xxxzz_zz, g_x_0_x_0_xxxzzz_x, g_x_0_x_0_xxxzzz_y, g_x_0_x_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxzzz_x[k] = -g_x_0_x_0_xxxzz_x[k] * ab_z + g_x_0_x_0_xxxzz_xz[k];

                g_x_0_x_0_xxxzzz_y[k] = -g_x_0_x_0_xxxzz_y[k] * ab_z + g_x_0_x_0_xxxzz_yz[k];

                g_x_0_x_0_xxxzzz_z[k] = -g_x_0_x_0_xxxzz_z[k] * ab_z + g_x_0_x_0_xxxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyy_x, g_x_0_x_0_xxyyy_xy, g_x_0_x_0_xxyyy_y, g_x_0_x_0_xxyyy_yy, g_x_0_x_0_xxyyy_yz, g_x_0_x_0_xxyyy_z, g_x_0_x_0_xxyyyy_x, g_x_0_x_0_xxyyyy_y, g_x_0_x_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyy_x[k] = -g_x_0_x_0_xxyyy_x[k] * ab_y + g_x_0_x_0_xxyyy_xy[k];

                g_x_0_x_0_xxyyyy_y[k] = -g_x_0_x_0_xxyyy_y[k] * ab_y + g_x_0_x_0_xxyyy_yy[k];

                g_x_0_x_0_xxyyyy_z[k] = -g_x_0_x_0_xxyyy_z[k] * ab_y + g_x_0_x_0_xxyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyyz_x, g_x_0_x_0_xxyyyz_y, g_x_0_x_0_xxyyyz_z, g_x_0_x_0_xxyyz_x, g_x_0_x_0_xxyyz_xy, g_x_0_x_0_xxyyz_y, g_x_0_x_0_xxyyz_yy, g_x_0_x_0_xxyyz_yz, g_x_0_x_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyyz_x[k] = -g_x_0_x_0_xxyyz_x[k] * ab_y + g_x_0_x_0_xxyyz_xy[k];

                g_x_0_x_0_xxyyyz_y[k] = -g_x_0_x_0_xxyyz_y[k] * ab_y + g_x_0_x_0_xxyyz_yy[k];

                g_x_0_x_0_xxyyyz_z[k] = -g_x_0_x_0_xxyyz_z[k] * ab_y + g_x_0_x_0_xxyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyyzz_x, g_x_0_x_0_xxyyzz_y, g_x_0_x_0_xxyyzz_z, g_x_0_x_0_xxyzz_x, g_x_0_x_0_xxyzz_xy, g_x_0_x_0_xxyzz_y, g_x_0_x_0_xxyzz_yy, g_x_0_x_0_xxyzz_yz, g_x_0_x_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyyzz_x[k] = -g_x_0_x_0_xxyzz_x[k] * ab_y + g_x_0_x_0_xxyzz_xy[k];

                g_x_0_x_0_xxyyzz_y[k] = -g_x_0_x_0_xxyzz_y[k] * ab_y + g_x_0_x_0_xxyzz_yy[k];

                g_x_0_x_0_xxyyzz_z[k] = -g_x_0_x_0_xxyzz_z[k] * ab_y + g_x_0_x_0_xxyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyzzz_x, g_x_0_x_0_xxyzzz_y, g_x_0_x_0_xxyzzz_z, g_x_0_x_0_xxzzz_x, g_x_0_x_0_xxzzz_xy, g_x_0_x_0_xxzzz_y, g_x_0_x_0_xxzzz_yy, g_x_0_x_0_xxzzz_yz, g_x_0_x_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyzzz_x[k] = -g_x_0_x_0_xxzzz_x[k] * ab_y + g_x_0_x_0_xxzzz_xy[k];

                g_x_0_x_0_xxyzzz_y[k] = -g_x_0_x_0_xxzzz_y[k] * ab_y + g_x_0_x_0_xxzzz_yy[k];

                g_x_0_x_0_xxyzzz_z[k] = -g_x_0_x_0_xxzzz_z[k] * ab_y + g_x_0_x_0_xxzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxzzz_x, g_x_0_x_0_xxzzz_xz, g_x_0_x_0_xxzzz_y, g_x_0_x_0_xxzzz_yz, g_x_0_x_0_xxzzz_z, g_x_0_x_0_xxzzz_zz, g_x_0_x_0_xxzzzz_x, g_x_0_x_0_xxzzzz_y, g_x_0_x_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzzzz_x[k] = -g_x_0_x_0_xxzzz_x[k] * ab_z + g_x_0_x_0_xxzzz_xz[k];

                g_x_0_x_0_xxzzzz_y[k] = -g_x_0_x_0_xxzzz_y[k] * ab_z + g_x_0_x_0_xxzzz_yz[k];

                g_x_0_x_0_xxzzzz_z[k] = -g_x_0_x_0_xxzzz_z[k] * ab_z + g_x_0_x_0_xxzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyy_x, g_x_0_x_0_xyyyy_xy, g_x_0_x_0_xyyyy_y, g_x_0_x_0_xyyyy_yy, g_x_0_x_0_xyyyy_yz, g_x_0_x_0_xyyyy_z, g_x_0_x_0_xyyyyy_x, g_x_0_x_0_xyyyyy_y, g_x_0_x_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyy_x[k] = -g_x_0_x_0_xyyyy_x[k] * ab_y + g_x_0_x_0_xyyyy_xy[k];

                g_x_0_x_0_xyyyyy_y[k] = -g_x_0_x_0_xyyyy_y[k] * ab_y + g_x_0_x_0_xyyyy_yy[k];

                g_x_0_x_0_xyyyyy_z[k] = -g_x_0_x_0_xyyyy_z[k] * ab_y + g_x_0_x_0_xyyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyyz_x, g_x_0_x_0_xyyyyz_y, g_x_0_x_0_xyyyyz_z, g_x_0_x_0_xyyyz_x, g_x_0_x_0_xyyyz_xy, g_x_0_x_0_xyyyz_y, g_x_0_x_0_xyyyz_yy, g_x_0_x_0_xyyyz_yz, g_x_0_x_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyyz_x[k] = -g_x_0_x_0_xyyyz_x[k] * ab_y + g_x_0_x_0_xyyyz_xy[k];

                g_x_0_x_0_xyyyyz_y[k] = -g_x_0_x_0_xyyyz_y[k] * ab_y + g_x_0_x_0_xyyyz_yy[k];

                g_x_0_x_0_xyyyyz_z[k] = -g_x_0_x_0_xyyyz_z[k] * ab_y + g_x_0_x_0_xyyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyyzz_x, g_x_0_x_0_xyyyzz_y, g_x_0_x_0_xyyyzz_z, g_x_0_x_0_xyyzz_x, g_x_0_x_0_xyyzz_xy, g_x_0_x_0_xyyzz_y, g_x_0_x_0_xyyzz_yy, g_x_0_x_0_xyyzz_yz, g_x_0_x_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyyzz_x[k] = -g_x_0_x_0_xyyzz_x[k] * ab_y + g_x_0_x_0_xyyzz_xy[k];

                g_x_0_x_0_xyyyzz_y[k] = -g_x_0_x_0_xyyzz_y[k] * ab_y + g_x_0_x_0_xyyzz_yy[k];

                g_x_0_x_0_xyyyzz_z[k] = -g_x_0_x_0_xyyzz_z[k] * ab_y + g_x_0_x_0_xyyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyzzz_x, g_x_0_x_0_xyyzzz_y, g_x_0_x_0_xyyzzz_z, g_x_0_x_0_xyzzz_x, g_x_0_x_0_xyzzz_xy, g_x_0_x_0_xyzzz_y, g_x_0_x_0_xyzzz_yy, g_x_0_x_0_xyzzz_yz, g_x_0_x_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyzzz_x[k] = -g_x_0_x_0_xyzzz_x[k] * ab_y + g_x_0_x_0_xyzzz_xy[k];

                g_x_0_x_0_xyyzzz_y[k] = -g_x_0_x_0_xyzzz_y[k] * ab_y + g_x_0_x_0_xyzzz_yy[k];

                g_x_0_x_0_xyyzzz_z[k] = -g_x_0_x_0_xyzzz_z[k] * ab_y + g_x_0_x_0_xyzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzzzz_x, g_x_0_x_0_xyzzzz_y, g_x_0_x_0_xyzzzz_z, g_x_0_x_0_xzzzz_x, g_x_0_x_0_xzzzz_xy, g_x_0_x_0_xzzzz_y, g_x_0_x_0_xzzzz_yy, g_x_0_x_0_xzzzz_yz, g_x_0_x_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzzzz_x[k] = -g_x_0_x_0_xzzzz_x[k] * ab_y + g_x_0_x_0_xzzzz_xy[k];

                g_x_0_x_0_xyzzzz_y[k] = -g_x_0_x_0_xzzzz_y[k] * ab_y + g_x_0_x_0_xzzzz_yy[k];

                g_x_0_x_0_xyzzzz_z[k] = -g_x_0_x_0_xzzzz_z[k] * ab_y + g_x_0_x_0_xzzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzzzz_x, g_x_0_x_0_xzzzz_xz, g_x_0_x_0_xzzzz_y, g_x_0_x_0_xzzzz_yz, g_x_0_x_0_xzzzz_z, g_x_0_x_0_xzzzz_zz, g_x_0_x_0_xzzzzz_x, g_x_0_x_0_xzzzzz_y, g_x_0_x_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzzzz_x[k] = -g_x_0_x_0_xzzzz_x[k] * ab_z + g_x_0_x_0_xzzzz_xz[k];

                g_x_0_x_0_xzzzzz_y[k] = -g_x_0_x_0_xzzzz_y[k] * ab_z + g_x_0_x_0_xzzzz_yz[k];

                g_x_0_x_0_xzzzzz_z[k] = -g_x_0_x_0_xzzzz_z[k] * ab_z + g_x_0_x_0_xzzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyy_x, g_x_0_x_0_yyyyy_xy, g_x_0_x_0_yyyyy_y, g_x_0_x_0_yyyyy_yy, g_x_0_x_0_yyyyy_yz, g_x_0_x_0_yyyyy_z, g_x_0_x_0_yyyyyy_x, g_x_0_x_0_yyyyyy_y, g_x_0_x_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyy_x[k] = -g_x_0_x_0_yyyyy_x[k] * ab_y + g_x_0_x_0_yyyyy_xy[k];

                g_x_0_x_0_yyyyyy_y[k] = -g_x_0_x_0_yyyyy_y[k] * ab_y + g_x_0_x_0_yyyyy_yy[k];

                g_x_0_x_0_yyyyyy_z[k] = -g_x_0_x_0_yyyyy_z[k] * ab_y + g_x_0_x_0_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyyz_x, g_x_0_x_0_yyyyyz_y, g_x_0_x_0_yyyyyz_z, g_x_0_x_0_yyyyz_x, g_x_0_x_0_yyyyz_xy, g_x_0_x_0_yyyyz_y, g_x_0_x_0_yyyyz_yy, g_x_0_x_0_yyyyz_yz, g_x_0_x_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyyz_x[k] = -g_x_0_x_0_yyyyz_x[k] * ab_y + g_x_0_x_0_yyyyz_xy[k];

                g_x_0_x_0_yyyyyz_y[k] = -g_x_0_x_0_yyyyz_y[k] * ab_y + g_x_0_x_0_yyyyz_yy[k];

                g_x_0_x_0_yyyyyz_z[k] = -g_x_0_x_0_yyyyz_z[k] * ab_y + g_x_0_x_0_yyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyyzz_x, g_x_0_x_0_yyyyzz_y, g_x_0_x_0_yyyyzz_z, g_x_0_x_0_yyyzz_x, g_x_0_x_0_yyyzz_xy, g_x_0_x_0_yyyzz_y, g_x_0_x_0_yyyzz_yy, g_x_0_x_0_yyyzz_yz, g_x_0_x_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyyzz_x[k] = -g_x_0_x_0_yyyzz_x[k] * ab_y + g_x_0_x_0_yyyzz_xy[k];

                g_x_0_x_0_yyyyzz_y[k] = -g_x_0_x_0_yyyzz_y[k] * ab_y + g_x_0_x_0_yyyzz_yy[k];

                g_x_0_x_0_yyyyzz_z[k] = -g_x_0_x_0_yyyzz_z[k] * ab_y + g_x_0_x_0_yyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyzzz_x, g_x_0_x_0_yyyzzz_y, g_x_0_x_0_yyyzzz_z, g_x_0_x_0_yyzzz_x, g_x_0_x_0_yyzzz_xy, g_x_0_x_0_yyzzz_y, g_x_0_x_0_yyzzz_yy, g_x_0_x_0_yyzzz_yz, g_x_0_x_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyzzz_x[k] = -g_x_0_x_0_yyzzz_x[k] * ab_y + g_x_0_x_0_yyzzz_xy[k];

                g_x_0_x_0_yyyzzz_y[k] = -g_x_0_x_0_yyzzz_y[k] * ab_y + g_x_0_x_0_yyzzz_yy[k];

                g_x_0_x_0_yyyzzz_z[k] = -g_x_0_x_0_yyzzz_z[k] * ab_y + g_x_0_x_0_yyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzzzz_x, g_x_0_x_0_yyzzzz_y, g_x_0_x_0_yyzzzz_z, g_x_0_x_0_yzzzz_x, g_x_0_x_0_yzzzz_xy, g_x_0_x_0_yzzzz_y, g_x_0_x_0_yzzzz_yy, g_x_0_x_0_yzzzz_yz, g_x_0_x_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzzzz_x[k] = -g_x_0_x_0_yzzzz_x[k] * ab_y + g_x_0_x_0_yzzzz_xy[k];

                g_x_0_x_0_yyzzzz_y[k] = -g_x_0_x_0_yzzzz_y[k] * ab_y + g_x_0_x_0_yzzzz_yy[k];

                g_x_0_x_0_yyzzzz_z[k] = -g_x_0_x_0_yzzzz_z[k] * ab_y + g_x_0_x_0_yzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzzzz_x, g_x_0_x_0_yzzzzz_y, g_x_0_x_0_yzzzzz_z, g_x_0_x_0_zzzzz_x, g_x_0_x_0_zzzzz_xy, g_x_0_x_0_zzzzz_y, g_x_0_x_0_zzzzz_yy, g_x_0_x_0_zzzzz_yz, g_x_0_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzzzz_x[k] = -g_x_0_x_0_zzzzz_x[k] * ab_y + g_x_0_x_0_zzzzz_xy[k];

                g_x_0_x_0_yzzzzz_y[k] = -g_x_0_x_0_zzzzz_y[k] * ab_y + g_x_0_x_0_zzzzz_yy[k];

                g_x_0_x_0_yzzzzz_z[k] = -g_x_0_x_0_zzzzz_z[k] * ab_y + g_x_0_x_0_zzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzzzz_x, g_x_0_x_0_zzzzz_xz, g_x_0_x_0_zzzzz_y, g_x_0_x_0_zzzzz_yz, g_x_0_x_0_zzzzz_z, g_x_0_x_0_zzzzz_zz, g_x_0_x_0_zzzzzz_x, g_x_0_x_0_zzzzzz_y, g_x_0_x_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzzzz_x[k] = -g_x_0_x_0_zzzzz_x[k] * ab_z + g_x_0_x_0_zzzzz_xz[k];

                g_x_0_x_0_zzzzzz_y[k] = -g_x_0_x_0_zzzzz_y[k] * ab_z + g_x_0_x_0_zzzzz_yz[k];

                g_x_0_x_0_zzzzzz_z[k] = -g_x_0_x_0_zzzzz_z[k] * ab_z + g_x_0_x_0_zzzzz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxxxx_x, g_0_0_y_0_xxxxx_y, g_0_0_y_0_xxxxx_z, g_x_0_y_0_xxxxx_x, g_x_0_y_0_xxxxx_xx, g_x_0_y_0_xxxxx_xy, g_x_0_y_0_xxxxx_xz, g_x_0_y_0_xxxxx_y, g_x_0_y_0_xxxxx_z, g_x_0_y_0_xxxxxx_x, g_x_0_y_0_xxxxxx_y, g_x_0_y_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxx_x[k] = -g_0_0_y_0_xxxxx_x[k] - g_x_0_y_0_xxxxx_x[k] * ab_x + g_x_0_y_0_xxxxx_xx[k];

                g_x_0_y_0_xxxxxx_y[k] = -g_0_0_y_0_xxxxx_y[k] - g_x_0_y_0_xxxxx_y[k] * ab_x + g_x_0_y_0_xxxxx_xy[k];

                g_x_0_y_0_xxxxxx_z[k] = -g_0_0_y_0_xxxxx_z[k] - g_x_0_y_0_xxxxx_z[k] * ab_x + g_x_0_y_0_xxxxx_xz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_x, g_x_0_y_0_xxxxx_xy, g_x_0_y_0_xxxxx_y, g_x_0_y_0_xxxxx_yy, g_x_0_y_0_xxxxx_yz, g_x_0_y_0_xxxxx_z, g_x_0_y_0_xxxxxy_x, g_x_0_y_0_xxxxxy_y, g_x_0_y_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxy_x[k] = -g_x_0_y_0_xxxxx_x[k] * ab_y + g_x_0_y_0_xxxxx_xy[k];

                g_x_0_y_0_xxxxxy_y[k] = -g_x_0_y_0_xxxxx_y[k] * ab_y + g_x_0_y_0_xxxxx_yy[k];

                g_x_0_y_0_xxxxxy_z[k] = -g_x_0_y_0_xxxxx_z[k] * ab_y + g_x_0_y_0_xxxxx_yz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxx_x, g_x_0_y_0_xxxxx_xz, g_x_0_y_0_xxxxx_y, g_x_0_y_0_xxxxx_yz, g_x_0_y_0_xxxxx_z, g_x_0_y_0_xxxxx_zz, g_x_0_y_0_xxxxxz_x, g_x_0_y_0_xxxxxz_y, g_x_0_y_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxxz_x[k] = -g_x_0_y_0_xxxxx_x[k] * ab_z + g_x_0_y_0_xxxxx_xz[k];

                g_x_0_y_0_xxxxxz_y[k] = -g_x_0_y_0_xxxxx_y[k] * ab_z + g_x_0_y_0_xxxxx_yz[k];

                g_x_0_y_0_xxxxxz_z[k] = -g_x_0_y_0_xxxxx_z[k] * ab_z + g_x_0_y_0_xxxxx_zz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxy_x, g_x_0_y_0_xxxxy_xy, g_x_0_y_0_xxxxy_y, g_x_0_y_0_xxxxy_yy, g_x_0_y_0_xxxxy_yz, g_x_0_y_0_xxxxy_z, g_x_0_y_0_xxxxyy_x, g_x_0_y_0_xxxxyy_y, g_x_0_y_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyy_x[k] = -g_x_0_y_0_xxxxy_x[k] * ab_y + g_x_0_y_0_xxxxy_xy[k];

                g_x_0_y_0_xxxxyy_y[k] = -g_x_0_y_0_xxxxy_y[k] * ab_y + g_x_0_y_0_xxxxy_yy[k];

                g_x_0_y_0_xxxxyy_z[k] = -g_x_0_y_0_xxxxy_z[k] * ab_y + g_x_0_y_0_xxxxy_yz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxyz_x, g_x_0_y_0_xxxxyz_y, g_x_0_y_0_xxxxyz_z, g_x_0_y_0_xxxxz_x, g_x_0_y_0_xxxxz_xy, g_x_0_y_0_xxxxz_y, g_x_0_y_0_xxxxz_yy, g_x_0_y_0_xxxxz_yz, g_x_0_y_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxyz_x[k] = -g_x_0_y_0_xxxxz_x[k] * ab_y + g_x_0_y_0_xxxxz_xy[k];

                g_x_0_y_0_xxxxyz_y[k] = -g_x_0_y_0_xxxxz_y[k] * ab_y + g_x_0_y_0_xxxxz_yy[k];

                g_x_0_y_0_xxxxyz_z[k] = -g_x_0_y_0_xxxxz_z[k] * ab_y + g_x_0_y_0_xxxxz_yz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_y_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxxz_x, g_x_0_y_0_xxxxz_xz, g_x_0_y_0_xxxxz_y, g_x_0_y_0_xxxxz_yz, g_x_0_y_0_xxxxz_z, g_x_0_y_0_xxxxz_zz, g_x_0_y_0_xxxxzz_x, g_x_0_y_0_xxxxzz_y, g_x_0_y_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxxzz_x[k] = -g_x_0_y_0_xxxxz_x[k] * ab_z + g_x_0_y_0_xxxxz_xz[k];

                g_x_0_y_0_xxxxzz_y[k] = -g_x_0_y_0_xxxxz_y[k] * ab_z + g_x_0_y_0_xxxxz_yz[k];

                g_x_0_y_0_xxxxzz_z[k] = -g_x_0_y_0_xxxxz_z[k] * ab_z + g_x_0_y_0_xxxxz_zz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyy_x, g_x_0_y_0_xxxyy_xy, g_x_0_y_0_xxxyy_y, g_x_0_y_0_xxxyy_yy, g_x_0_y_0_xxxyy_yz, g_x_0_y_0_xxxyy_z, g_x_0_y_0_xxxyyy_x, g_x_0_y_0_xxxyyy_y, g_x_0_y_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyy_x[k] = -g_x_0_y_0_xxxyy_x[k] * ab_y + g_x_0_y_0_xxxyy_xy[k];

                g_x_0_y_0_xxxyyy_y[k] = -g_x_0_y_0_xxxyy_y[k] * ab_y + g_x_0_y_0_xxxyy_yy[k];

                g_x_0_y_0_xxxyyy_z[k] = -g_x_0_y_0_xxxyy_z[k] * ab_y + g_x_0_y_0_xxxyy_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyyz_x, g_x_0_y_0_xxxyyz_y, g_x_0_y_0_xxxyyz_z, g_x_0_y_0_xxxyz_x, g_x_0_y_0_xxxyz_xy, g_x_0_y_0_xxxyz_y, g_x_0_y_0_xxxyz_yy, g_x_0_y_0_xxxyz_yz, g_x_0_y_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyyz_x[k] = -g_x_0_y_0_xxxyz_x[k] * ab_y + g_x_0_y_0_xxxyz_xy[k];

                g_x_0_y_0_xxxyyz_y[k] = -g_x_0_y_0_xxxyz_y[k] * ab_y + g_x_0_y_0_xxxyz_yy[k];

                g_x_0_y_0_xxxyyz_z[k] = -g_x_0_y_0_xxxyz_z[k] * ab_y + g_x_0_y_0_xxxyz_yz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_y_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxyzz_x, g_x_0_y_0_xxxyzz_y, g_x_0_y_0_xxxyzz_z, g_x_0_y_0_xxxzz_x, g_x_0_y_0_xxxzz_xy, g_x_0_y_0_xxxzz_y, g_x_0_y_0_xxxzz_yy, g_x_0_y_0_xxxzz_yz, g_x_0_y_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxyzz_x[k] = -g_x_0_y_0_xxxzz_x[k] * ab_y + g_x_0_y_0_xxxzz_xy[k];

                g_x_0_y_0_xxxyzz_y[k] = -g_x_0_y_0_xxxzz_y[k] * ab_y + g_x_0_y_0_xxxzz_yy[k];

                g_x_0_y_0_xxxyzz_z[k] = -g_x_0_y_0_xxxzz_z[k] * ab_y + g_x_0_y_0_xxxzz_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_y_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxxzz_x, g_x_0_y_0_xxxzz_xz, g_x_0_y_0_xxxzz_y, g_x_0_y_0_xxxzz_yz, g_x_0_y_0_xxxzz_z, g_x_0_y_0_xxxzz_zz, g_x_0_y_0_xxxzzz_x, g_x_0_y_0_xxxzzz_y, g_x_0_y_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxzzz_x[k] = -g_x_0_y_0_xxxzz_x[k] * ab_z + g_x_0_y_0_xxxzz_xz[k];

                g_x_0_y_0_xxxzzz_y[k] = -g_x_0_y_0_xxxzz_y[k] * ab_z + g_x_0_y_0_xxxzz_yz[k];

                g_x_0_y_0_xxxzzz_z[k] = -g_x_0_y_0_xxxzz_z[k] * ab_z + g_x_0_y_0_xxxzz_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyy_x, g_x_0_y_0_xxyyy_xy, g_x_0_y_0_xxyyy_y, g_x_0_y_0_xxyyy_yy, g_x_0_y_0_xxyyy_yz, g_x_0_y_0_xxyyy_z, g_x_0_y_0_xxyyyy_x, g_x_0_y_0_xxyyyy_y, g_x_0_y_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyy_x[k] = -g_x_0_y_0_xxyyy_x[k] * ab_y + g_x_0_y_0_xxyyy_xy[k];

                g_x_0_y_0_xxyyyy_y[k] = -g_x_0_y_0_xxyyy_y[k] * ab_y + g_x_0_y_0_xxyyy_yy[k];

                g_x_0_y_0_xxyyyy_z[k] = -g_x_0_y_0_xxyyy_z[k] * ab_y + g_x_0_y_0_xxyyy_yz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyyz_x, g_x_0_y_0_xxyyyz_y, g_x_0_y_0_xxyyyz_z, g_x_0_y_0_xxyyz_x, g_x_0_y_0_xxyyz_xy, g_x_0_y_0_xxyyz_y, g_x_0_y_0_xxyyz_yy, g_x_0_y_0_xxyyz_yz, g_x_0_y_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyyz_x[k] = -g_x_0_y_0_xxyyz_x[k] * ab_y + g_x_0_y_0_xxyyz_xy[k];

                g_x_0_y_0_xxyyyz_y[k] = -g_x_0_y_0_xxyyz_y[k] * ab_y + g_x_0_y_0_xxyyz_yy[k];

                g_x_0_y_0_xxyyyz_z[k] = -g_x_0_y_0_xxyyz_z[k] * ab_y + g_x_0_y_0_xxyyz_yz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_y_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyyzz_x, g_x_0_y_0_xxyyzz_y, g_x_0_y_0_xxyyzz_z, g_x_0_y_0_xxyzz_x, g_x_0_y_0_xxyzz_xy, g_x_0_y_0_xxyzz_y, g_x_0_y_0_xxyzz_yy, g_x_0_y_0_xxyzz_yz, g_x_0_y_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyyzz_x[k] = -g_x_0_y_0_xxyzz_x[k] * ab_y + g_x_0_y_0_xxyzz_xy[k];

                g_x_0_y_0_xxyyzz_y[k] = -g_x_0_y_0_xxyzz_y[k] * ab_y + g_x_0_y_0_xxyzz_yy[k];

                g_x_0_y_0_xxyyzz_z[k] = -g_x_0_y_0_xxyzz_z[k] * ab_y + g_x_0_y_0_xxyzz_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_y_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyzzz_x, g_x_0_y_0_xxyzzz_y, g_x_0_y_0_xxyzzz_z, g_x_0_y_0_xxzzz_x, g_x_0_y_0_xxzzz_xy, g_x_0_y_0_xxzzz_y, g_x_0_y_0_xxzzz_yy, g_x_0_y_0_xxzzz_yz, g_x_0_y_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyzzz_x[k] = -g_x_0_y_0_xxzzz_x[k] * ab_y + g_x_0_y_0_xxzzz_xy[k];

                g_x_0_y_0_xxyzzz_y[k] = -g_x_0_y_0_xxzzz_y[k] * ab_y + g_x_0_y_0_xxzzz_yy[k];

                g_x_0_y_0_xxyzzz_z[k] = -g_x_0_y_0_xxzzz_z[k] * ab_y + g_x_0_y_0_xxzzz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_y_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxzzz_x, g_x_0_y_0_xxzzz_xz, g_x_0_y_0_xxzzz_y, g_x_0_y_0_xxzzz_yz, g_x_0_y_0_xxzzz_z, g_x_0_y_0_xxzzz_zz, g_x_0_y_0_xxzzzz_x, g_x_0_y_0_xxzzzz_y, g_x_0_y_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzzzz_x[k] = -g_x_0_y_0_xxzzz_x[k] * ab_z + g_x_0_y_0_xxzzz_xz[k];

                g_x_0_y_0_xxzzzz_y[k] = -g_x_0_y_0_xxzzz_y[k] * ab_z + g_x_0_y_0_xxzzz_yz[k];

                g_x_0_y_0_xxzzzz_z[k] = -g_x_0_y_0_xxzzz_z[k] * ab_z + g_x_0_y_0_xxzzz_zz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyy_x, g_x_0_y_0_xyyyy_xy, g_x_0_y_0_xyyyy_y, g_x_0_y_0_xyyyy_yy, g_x_0_y_0_xyyyy_yz, g_x_0_y_0_xyyyy_z, g_x_0_y_0_xyyyyy_x, g_x_0_y_0_xyyyyy_y, g_x_0_y_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyy_x[k] = -g_x_0_y_0_xyyyy_x[k] * ab_y + g_x_0_y_0_xyyyy_xy[k];

                g_x_0_y_0_xyyyyy_y[k] = -g_x_0_y_0_xyyyy_y[k] * ab_y + g_x_0_y_0_xyyyy_yy[k];

                g_x_0_y_0_xyyyyy_z[k] = -g_x_0_y_0_xyyyy_z[k] * ab_y + g_x_0_y_0_xyyyy_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyyz_x, g_x_0_y_0_xyyyyz_y, g_x_0_y_0_xyyyyz_z, g_x_0_y_0_xyyyz_x, g_x_0_y_0_xyyyz_xy, g_x_0_y_0_xyyyz_y, g_x_0_y_0_xyyyz_yy, g_x_0_y_0_xyyyz_yz, g_x_0_y_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyyz_x[k] = -g_x_0_y_0_xyyyz_x[k] * ab_y + g_x_0_y_0_xyyyz_xy[k];

                g_x_0_y_0_xyyyyz_y[k] = -g_x_0_y_0_xyyyz_y[k] * ab_y + g_x_0_y_0_xyyyz_yy[k];

                g_x_0_y_0_xyyyyz_z[k] = -g_x_0_y_0_xyyyz_z[k] * ab_y + g_x_0_y_0_xyyyz_yz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_y_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyyzz_x, g_x_0_y_0_xyyyzz_y, g_x_0_y_0_xyyyzz_z, g_x_0_y_0_xyyzz_x, g_x_0_y_0_xyyzz_xy, g_x_0_y_0_xyyzz_y, g_x_0_y_0_xyyzz_yy, g_x_0_y_0_xyyzz_yz, g_x_0_y_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyyzz_x[k] = -g_x_0_y_0_xyyzz_x[k] * ab_y + g_x_0_y_0_xyyzz_xy[k];

                g_x_0_y_0_xyyyzz_y[k] = -g_x_0_y_0_xyyzz_y[k] * ab_y + g_x_0_y_0_xyyzz_yy[k];

                g_x_0_y_0_xyyyzz_z[k] = -g_x_0_y_0_xyyzz_z[k] * ab_y + g_x_0_y_0_xyyzz_yz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_y_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyzzz_x, g_x_0_y_0_xyyzzz_y, g_x_0_y_0_xyyzzz_z, g_x_0_y_0_xyzzz_x, g_x_0_y_0_xyzzz_xy, g_x_0_y_0_xyzzz_y, g_x_0_y_0_xyzzz_yy, g_x_0_y_0_xyzzz_yz, g_x_0_y_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyzzz_x[k] = -g_x_0_y_0_xyzzz_x[k] * ab_y + g_x_0_y_0_xyzzz_xy[k];

                g_x_0_y_0_xyyzzz_y[k] = -g_x_0_y_0_xyzzz_y[k] * ab_y + g_x_0_y_0_xyzzz_yy[k];

                g_x_0_y_0_xyyzzz_z[k] = -g_x_0_y_0_xyzzz_z[k] * ab_y + g_x_0_y_0_xyzzz_yz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_y_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzzzz_x, g_x_0_y_0_xyzzzz_y, g_x_0_y_0_xyzzzz_z, g_x_0_y_0_xzzzz_x, g_x_0_y_0_xzzzz_xy, g_x_0_y_0_xzzzz_y, g_x_0_y_0_xzzzz_yy, g_x_0_y_0_xzzzz_yz, g_x_0_y_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzzzz_x[k] = -g_x_0_y_0_xzzzz_x[k] * ab_y + g_x_0_y_0_xzzzz_xy[k];

                g_x_0_y_0_xyzzzz_y[k] = -g_x_0_y_0_xzzzz_y[k] * ab_y + g_x_0_y_0_xzzzz_yy[k];

                g_x_0_y_0_xyzzzz_z[k] = -g_x_0_y_0_xzzzz_z[k] * ab_y + g_x_0_y_0_xzzzz_yz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_y_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzzzz_x, g_x_0_y_0_xzzzz_xz, g_x_0_y_0_xzzzz_y, g_x_0_y_0_xzzzz_yz, g_x_0_y_0_xzzzz_z, g_x_0_y_0_xzzzz_zz, g_x_0_y_0_xzzzzz_x, g_x_0_y_0_xzzzzz_y, g_x_0_y_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzzzz_x[k] = -g_x_0_y_0_xzzzz_x[k] * ab_z + g_x_0_y_0_xzzzz_xz[k];

                g_x_0_y_0_xzzzzz_y[k] = -g_x_0_y_0_xzzzz_y[k] * ab_z + g_x_0_y_0_xzzzz_yz[k];

                g_x_0_y_0_xzzzzz_z[k] = -g_x_0_y_0_xzzzz_z[k] * ab_z + g_x_0_y_0_xzzzz_zz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyy_x, g_x_0_y_0_yyyyy_xy, g_x_0_y_0_yyyyy_y, g_x_0_y_0_yyyyy_yy, g_x_0_y_0_yyyyy_yz, g_x_0_y_0_yyyyy_z, g_x_0_y_0_yyyyyy_x, g_x_0_y_0_yyyyyy_y, g_x_0_y_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyy_x[k] = -g_x_0_y_0_yyyyy_x[k] * ab_y + g_x_0_y_0_yyyyy_xy[k];

                g_x_0_y_0_yyyyyy_y[k] = -g_x_0_y_0_yyyyy_y[k] * ab_y + g_x_0_y_0_yyyyy_yy[k];

                g_x_0_y_0_yyyyyy_z[k] = -g_x_0_y_0_yyyyy_z[k] * ab_y + g_x_0_y_0_yyyyy_yz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyyz_x, g_x_0_y_0_yyyyyz_y, g_x_0_y_0_yyyyyz_z, g_x_0_y_0_yyyyz_x, g_x_0_y_0_yyyyz_xy, g_x_0_y_0_yyyyz_y, g_x_0_y_0_yyyyz_yy, g_x_0_y_0_yyyyz_yz, g_x_0_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyyz_x[k] = -g_x_0_y_0_yyyyz_x[k] * ab_y + g_x_0_y_0_yyyyz_xy[k];

                g_x_0_y_0_yyyyyz_y[k] = -g_x_0_y_0_yyyyz_y[k] * ab_y + g_x_0_y_0_yyyyz_yy[k];

                g_x_0_y_0_yyyyyz_z[k] = -g_x_0_y_0_yyyyz_z[k] * ab_y + g_x_0_y_0_yyyyz_yz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_y_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyyzz_x, g_x_0_y_0_yyyyzz_y, g_x_0_y_0_yyyyzz_z, g_x_0_y_0_yyyzz_x, g_x_0_y_0_yyyzz_xy, g_x_0_y_0_yyyzz_y, g_x_0_y_0_yyyzz_yy, g_x_0_y_0_yyyzz_yz, g_x_0_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyyzz_x[k] = -g_x_0_y_0_yyyzz_x[k] * ab_y + g_x_0_y_0_yyyzz_xy[k];

                g_x_0_y_0_yyyyzz_y[k] = -g_x_0_y_0_yyyzz_y[k] * ab_y + g_x_0_y_0_yyyzz_yy[k];

                g_x_0_y_0_yyyyzz_z[k] = -g_x_0_y_0_yyyzz_z[k] * ab_y + g_x_0_y_0_yyyzz_yz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_y_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyzzz_x, g_x_0_y_0_yyyzzz_y, g_x_0_y_0_yyyzzz_z, g_x_0_y_0_yyzzz_x, g_x_0_y_0_yyzzz_xy, g_x_0_y_0_yyzzz_y, g_x_0_y_0_yyzzz_yy, g_x_0_y_0_yyzzz_yz, g_x_0_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyzzz_x[k] = -g_x_0_y_0_yyzzz_x[k] * ab_y + g_x_0_y_0_yyzzz_xy[k];

                g_x_0_y_0_yyyzzz_y[k] = -g_x_0_y_0_yyzzz_y[k] * ab_y + g_x_0_y_0_yyzzz_yy[k];

                g_x_0_y_0_yyyzzz_z[k] = -g_x_0_y_0_yyzzz_z[k] * ab_y + g_x_0_y_0_yyzzz_yz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_y_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzzzz_x, g_x_0_y_0_yyzzzz_y, g_x_0_y_0_yyzzzz_z, g_x_0_y_0_yzzzz_x, g_x_0_y_0_yzzzz_xy, g_x_0_y_0_yzzzz_y, g_x_0_y_0_yzzzz_yy, g_x_0_y_0_yzzzz_yz, g_x_0_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzzzz_x[k] = -g_x_0_y_0_yzzzz_x[k] * ab_y + g_x_0_y_0_yzzzz_xy[k];

                g_x_0_y_0_yyzzzz_y[k] = -g_x_0_y_0_yzzzz_y[k] * ab_y + g_x_0_y_0_yzzzz_yy[k];

                g_x_0_y_0_yyzzzz_z[k] = -g_x_0_y_0_yzzzz_z[k] * ab_y + g_x_0_y_0_yzzzz_yz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_y_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzzzz_x, g_x_0_y_0_yzzzzz_y, g_x_0_y_0_yzzzzz_z, g_x_0_y_0_zzzzz_x, g_x_0_y_0_zzzzz_xy, g_x_0_y_0_zzzzz_y, g_x_0_y_0_zzzzz_yy, g_x_0_y_0_zzzzz_yz, g_x_0_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzzzz_x[k] = -g_x_0_y_0_zzzzz_x[k] * ab_y + g_x_0_y_0_zzzzz_xy[k];

                g_x_0_y_0_yzzzzz_y[k] = -g_x_0_y_0_zzzzz_y[k] * ab_y + g_x_0_y_0_zzzzz_yy[k];

                g_x_0_y_0_yzzzzz_z[k] = -g_x_0_y_0_zzzzz_z[k] * ab_y + g_x_0_y_0_zzzzz_yz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_y_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzzzz_x, g_x_0_y_0_zzzzz_xz, g_x_0_y_0_zzzzz_y, g_x_0_y_0_zzzzz_yz, g_x_0_y_0_zzzzz_z, g_x_0_y_0_zzzzz_zz, g_x_0_y_0_zzzzzz_x, g_x_0_y_0_zzzzzz_y, g_x_0_y_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzzzz_x[k] = -g_x_0_y_0_zzzzz_x[k] * ab_z + g_x_0_y_0_zzzzz_xz[k];

                g_x_0_y_0_zzzzzz_y[k] = -g_x_0_y_0_zzzzz_y[k] * ab_z + g_x_0_y_0_zzzzz_yz[k];

                g_x_0_y_0_zzzzzz_z[k] = -g_x_0_y_0_zzzzz_z[k] * ab_z + g_x_0_y_0_zzzzz_zz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxxxx_x, g_0_0_z_0_xxxxx_y, g_0_0_z_0_xxxxx_z, g_x_0_z_0_xxxxx_x, g_x_0_z_0_xxxxx_xx, g_x_0_z_0_xxxxx_xy, g_x_0_z_0_xxxxx_xz, g_x_0_z_0_xxxxx_y, g_x_0_z_0_xxxxx_z, g_x_0_z_0_xxxxxx_x, g_x_0_z_0_xxxxxx_y, g_x_0_z_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxx_x[k] = -g_0_0_z_0_xxxxx_x[k] - g_x_0_z_0_xxxxx_x[k] * ab_x + g_x_0_z_0_xxxxx_xx[k];

                g_x_0_z_0_xxxxxx_y[k] = -g_0_0_z_0_xxxxx_y[k] - g_x_0_z_0_xxxxx_y[k] * ab_x + g_x_0_z_0_xxxxx_xy[k];

                g_x_0_z_0_xxxxxx_z[k] = -g_0_0_z_0_xxxxx_z[k] - g_x_0_z_0_xxxxx_z[k] * ab_x + g_x_0_z_0_xxxxx_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_x, g_x_0_z_0_xxxxx_xy, g_x_0_z_0_xxxxx_y, g_x_0_z_0_xxxxx_yy, g_x_0_z_0_xxxxx_yz, g_x_0_z_0_xxxxx_z, g_x_0_z_0_xxxxxy_x, g_x_0_z_0_xxxxxy_y, g_x_0_z_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxy_x[k] = -g_x_0_z_0_xxxxx_x[k] * ab_y + g_x_0_z_0_xxxxx_xy[k];

                g_x_0_z_0_xxxxxy_y[k] = -g_x_0_z_0_xxxxx_y[k] * ab_y + g_x_0_z_0_xxxxx_yy[k];

                g_x_0_z_0_xxxxxy_z[k] = -g_x_0_z_0_xxxxx_z[k] * ab_y + g_x_0_z_0_xxxxx_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxx_x, g_x_0_z_0_xxxxx_xz, g_x_0_z_0_xxxxx_y, g_x_0_z_0_xxxxx_yz, g_x_0_z_0_xxxxx_z, g_x_0_z_0_xxxxx_zz, g_x_0_z_0_xxxxxz_x, g_x_0_z_0_xxxxxz_y, g_x_0_z_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxxz_x[k] = -g_x_0_z_0_xxxxx_x[k] * ab_z + g_x_0_z_0_xxxxx_xz[k];

                g_x_0_z_0_xxxxxz_y[k] = -g_x_0_z_0_xxxxx_y[k] * ab_z + g_x_0_z_0_xxxxx_yz[k];

                g_x_0_z_0_xxxxxz_z[k] = -g_x_0_z_0_xxxxx_z[k] * ab_z + g_x_0_z_0_xxxxx_zz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxy_x, g_x_0_z_0_xxxxy_xy, g_x_0_z_0_xxxxy_y, g_x_0_z_0_xxxxy_yy, g_x_0_z_0_xxxxy_yz, g_x_0_z_0_xxxxy_z, g_x_0_z_0_xxxxyy_x, g_x_0_z_0_xxxxyy_y, g_x_0_z_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyy_x[k] = -g_x_0_z_0_xxxxy_x[k] * ab_y + g_x_0_z_0_xxxxy_xy[k];

                g_x_0_z_0_xxxxyy_y[k] = -g_x_0_z_0_xxxxy_y[k] * ab_y + g_x_0_z_0_xxxxy_yy[k];

                g_x_0_z_0_xxxxyy_z[k] = -g_x_0_z_0_xxxxy_z[k] * ab_y + g_x_0_z_0_xxxxy_yz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxyz_x, g_x_0_z_0_xxxxyz_y, g_x_0_z_0_xxxxyz_z, g_x_0_z_0_xxxxz_x, g_x_0_z_0_xxxxz_xy, g_x_0_z_0_xxxxz_y, g_x_0_z_0_xxxxz_yy, g_x_0_z_0_xxxxz_yz, g_x_0_z_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxyz_x[k] = -g_x_0_z_0_xxxxz_x[k] * ab_y + g_x_0_z_0_xxxxz_xy[k];

                g_x_0_z_0_xxxxyz_y[k] = -g_x_0_z_0_xxxxz_y[k] * ab_y + g_x_0_z_0_xxxxz_yy[k];

                g_x_0_z_0_xxxxyz_z[k] = -g_x_0_z_0_xxxxz_z[k] * ab_y + g_x_0_z_0_xxxxz_yz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_z_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxxz_x, g_x_0_z_0_xxxxz_xz, g_x_0_z_0_xxxxz_y, g_x_0_z_0_xxxxz_yz, g_x_0_z_0_xxxxz_z, g_x_0_z_0_xxxxz_zz, g_x_0_z_0_xxxxzz_x, g_x_0_z_0_xxxxzz_y, g_x_0_z_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxxzz_x[k] = -g_x_0_z_0_xxxxz_x[k] * ab_z + g_x_0_z_0_xxxxz_xz[k];

                g_x_0_z_0_xxxxzz_y[k] = -g_x_0_z_0_xxxxz_y[k] * ab_z + g_x_0_z_0_xxxxz_yz[k];

                g_x_0_z_0_xxxxzz_z[k] = -g_x_0_z_0_xxxxz_z[k] * ab_z + g_x_0_z_0_xxxxz_zz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyy_x, g_x_0_z_0_xxxyy_xy, g_x_0_z_0_xxxyy_y, g_x_0_z_0_xxxyy_yy, g_x_0_z_0_xxxyy_yz, g_x_0_z_0_xxxyy_z, g_x_0_z_0_xxxyyy_x, g_x_0_z_0_xxxyyy_y, g_x_0_z_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyy_x[k] = -g_x_0_z_0_xxxyy_x[k] * ab_y + g_x_0_z_0_xxxyy_xy[k];

                g_x_0_z_0_xxxyyy_y[k] = -g_x_0_z_0_xxxyy_y[k] * ab_y + g_x_0_z_0_xxxyy_yy[k];

                g_x_0_z_0_xxxyyy_z[k] = -g_x_0_z_0_xxxyy_z[k] * ab_y + g_x_0_z_0_xxxyy_yz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyyz_x, g_x_0_z_0_xxxyyz_y, g_x_0_z_0_xxxyyz_z, g_x_0_z_0_xxxyz_x, g_x_0_z_0_xxxyz_xy, g_x_0_z_0_xxxyz_y, g_x_0_z_0_xxxyz_yy, g_x_0_z_0_xxxyz_yz, g_x_0_z_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyyz_x[k] = -g_x_0_z_0_xxxyz_x[k] * ab_y + g_x_0_z_0_xxxyz_xy[k];

                g_x_0_z_0_xxxyyz_y[k] = -g_x_0_z_0_xxxyz_y[k] * ab_y + g_x_0_z_0_xxxyz_yy[k];

                g_x_0_z_0_xxxyyz_z[k] = -g_x_0_z_0_xxxyz_z[k] * ab_y + g_x_0_z_0_xxxyz_yz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_z_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxyzz_x, g_x_0_z_0_xxxyzz_y, g_x_0_z_0_xxxyzz_z, g_x_0_z_0_xxxzz_x, g_x_0_z_0_xxxzz_xy, g_x_0_z_0_xxxzz_y, g_x_0_z_0_xxxzz_yy, g_x_0_z_0_xxxzz_yz, g_x_0_z_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxyzz_x[k] = -g_x_0_z_0_xxxzz_x[k] * ab_y + g_x_0_z_0_xxxzz_xy[k];

                g_x_0_z_0_xxxyzz_y[k] = -g_x_0_z_0_xxxzz_y[k] * ab_y + g_x_0_z_0_xxxzz_yy[k];

                g_x_0_z_0_xxxyzz_z[k] = -g_x_0_z_0_xxxzz_z[k] * ab_y + g_x_0_z_0_xxxzz_yz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_z_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxxzz_x, g_x_0_z_0_xxxzz_xz, g_x_0_z_0_xxxzz_y, g_x_0_z_0_xxxzz_yz, g_x_0_z_0_xxxzz_z, g_x_0_z_0_xxxzz_zz, g_x_0_z_0_xxxzzz_x, g_x_0_z_0_xxxzzz_y, g_x_0_z_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxzzz_x[k] = -g_x_0_z_0_xxxzz_x[k] * ab_z + g_x_0_z_0_xxxzz_xz[k];

                g_x_0_z_0_xxxzzz_y[k] = -g_x_0_z_0_xxxzz_y[k] * ab_z + g_x_0_z_0_xxxzz_yz[k];

                g_x_0_z_0_xxxzzz_z[k] = -g_x_0_z_0_xxxzz_z[k] * ab_z + g_x_0_z_0_xxxzz_zz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyy_x, g_x_0_z_0_xxyyy_xy, g_x_0_z_0_xxyyy_y, g_x_0_z_0_xxyyy_yy, g_x_0_z_0_xxyyy_yz, g_x_0_z_0_xxyyy_z, g_x_0_z_0_xxyyyy_x, g_x_0_z_0_xxyyyy_y, g_x_0_z_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyy_x[k] = -g_x_0_z_0_xxyyy_x[k] * ab_y + g_x_0_z_0_xxyyy_xy[k];

                g_x_0_z_0_xxyyyy_y[k] = -g_x_0_z_0_xxyyy_y[k] * ab_y + g_x_0_z_0_xxyyy_yy[k];

                g_x_0_z_0_xxyyyy_z[k] = -g_x_0_z_0_xxyyy_z[k] * ab_y + g_x_0_z_0_xxyyy_yz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyyz_x, g_x_0_z_0_xxyyyz_y, g_x_0_z_0_xxyyyz_z, g_x_0_z_0_xxyyz_x, g_x_0_z_0_xxyyz_xy, g_x_0_z_0_xxyyz_y, g_x_0_z_0_xxyyz_yy, g_x_0_z_0_xxyyz_yz, g_x_0_z_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyyz_x[k] = -g_x_0_z_0_xxyyz_x[k] * ab_y + g_x_0_z_0_xxyyz_xy[k];

                g_x_0_z_0_xxyyyz_y[k] = -g_x_0_z_0_xxyyz_y[k] * ab_y + g_x_0_z_0_xxyyz_yy[k];

                g_x_0_z_0_xxyyyz_z[k] = -g_x_0_z_0_xxyyz_z[k] * ab_y + g_x_0_z_0_xxyyz_yz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_z_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyyzz_x, g_x_0_z_0_xxyyzz_y, g_x_0_z_0_xxyyzz_z, g_x_0_z_0_xxyzz_x, g_x_0_z_0_xxyzz_xy, g_x_0_z_0_xxyzz_y, g_x_0_z_0_xxyzz_yy, g_x_0_z_0_xxyzz_yz, g_x_0_z_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyyzz_x[k] = -g_x_0_z_0_xxyzz_x[k] * ab_y + g_x_0_z_0_xxyzz_xy[k];

                g_x_0_z_0_xxyyzz_y[k] = -g_x_0_z_0_xxyzz_y[k] * ab_y + g_x_0_z_0_xxyzz_yy[k];

                g_x_0_z_0_xxyyzz_z[k] = -g_x_0_z_0_xxyzz_z[k] * ab_y + g_x_0_z_0_xxyzz_yz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_z_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyzzz_x, g_x_0_z_0_xxyzzz_y, g_x_0_z_0_xxyzzz_z, g_x_0_z_0_xxzzz_x, g_x_0_z_0_xxzzz_xy, g_x_0_z_0_xxzzz_y, g_x_0_z_0_xxzzz_yy, g_x_0_z_0_xxzzz_yz, g_x_0_z_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyzzz_x[k] = -g_x_0_z_0_xxzzz_x[k] * ab_y + g_x_0_z_0_xxzzz_xy[k];

                g_x_0_z_0_xxyzzz_y[k] = -g_x_0_z_0_xxzzz_y[k] * ab_y + g_x_0_z_0_xxzzz_yy[k];

                g_x_0_z_0_xxyzzz_z[k] = -g_x_0_z_0_xxzzz_z[k] * ab_y + g_x_0_z_0_xxzzz_yz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_z_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxzzz_x, g_x_0_z_0_xxzzz_xz, g_x_0_z_0_xxzzz_y, g_x_0_z_0_xxzzz_yz, g_x_0_z_0_xxzzz_z, g_x_0_z_0_xxzzz_zz, g_x_0_z_0_xxzzzz_x, g_x_0_z_0_xxzzzz_y, g_x_0_z_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzzzz_x[k] = -g_x_0_z_0_xxzzz_x[k] * ab_z + g_x_0_z_0_xxzzz_xz[k];

                g_x_0_z_0_xxzzzz_y[k] = -g_x_0_z_0_xxzzz_y[k] * ab_z + g_x_0_z_0_xxzzz_yz[k];

                g_x_0_z_0_xxzzzz_z[k] = -g_x_0_z_0_xxzzz_z[k] * ab_z + g_x_0_z_0_xxzzz_zz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyy_x, g_x_0_z_0_xyyyy_xy, g_x_0_z_0_xyyyy_y, g_x_0_z_0_xyyyy_yy, g_x_0_z_0_xyyyy_yz, g_x_0_z_0_xyyyy_z, g_x_0_z_0_xyyyyy_x, g_x_0_z_0_xyyyyy_y, g_x_0_z_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyy_x[k] = -g_x_0_z_0_xyyyy_x[k] * ab_y + g_x_0_z_0_xyyyy_xy[k];

                g_x_0_z_0_xyyyyy_y[k] = -g_x_0_z_0_xyyyy_y[k] * ab_y + g_x_0_z_0_xyyyy_yy[k];

                g_x_0_z_0_xyyyyy_z[k] = -g_x_0_z_0_xyyyy_z[k] * ab_y + g_x_0_z_0_xyyyy_yz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyyz_x, g_x_0_z_0_xyyyyz_y, g_x_0_z_0_xyyyyz_z, g_x_0_z_0_xyyyz_x, g_x_0_z_0_xyyyz_xy, g_x_0_z_0_xyyyz_y, g_x_0_z_0_xyyyz_yy, g_x_0_z_0_xyyyz_yz, g_x_0_z_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyyz_x[k] = -g_x_0_z_0_xyyyz_x[k] * ab_y + g_x_0_z_0_xyyyz_xy[k];

                g_x_0_z_0_xyyyyz_y[k] = -g_x_0_z_0_xyyyz_y[k] * ab_y + g_x_0_z_0_xyyyz_yy[k];

                g_x_0_z_0_xyyyyz_z[k] = -g_x_0_z_0_xyyyz_z[k] * ab_y + g_x_0_z_0_xyyyz_yz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_z_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyyzz_x, g_x_0_z_0_xyyyzz_y, g_x_0_z_0_xyyyzz_z, g_x_0_z_0_xyyzz_x, g_x_0_z_0_xyyzz_xy, g_x_0_z_0_xyyzz_y, g_x_0_z_0_xyyzz_yy, g_x_0_z_0_xyyzz_yz, g_x_0_z_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyyzz_x[k] = -g_x_0_z_0_xyyzz_x[k] * ab_y + g_x_0_z_0_xyyzz_xy[k];

                g_x_0_z_0_xyyyzz_y[k] = -g_x_0_z_0_xyyzz_y[k] * ab_y + g_x_0_z_0_xyyzz_yy[k];

                g_x_0_z_0_xyyyzz_z[k] = -g_x_0_z_0_xyyzz_z[k] * ab_y + g_x_0_z_0_xyyzz_yz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_z_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyzzz_x, g_x_0_z_0_xyyzzz_y, g_x_0_z_0_xyyzzz_z, g_x_0_z_0_xyzzz_x, g_x_0_z_0_xyzzz_xy, g_x_0_z_0_xyzzz_y, g_x_0_z_0_xyzzz_yy, g_x_0_z_0_xyzzz_yz, g_x_0_z_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyzzz_x[k] = -g_x_0_z_0_xyzzz_x[k] * ab_y + g_x_0_z_0_xyzzz_xy[k];

                g_x_0_z_0_xyyzzz_y[k] = -g_x_0_z_0_xyzzz_y[k] * ab_y + g_x_0_z_0_xyzzz_yy[k];

                g_x_0_z_0_xyyzzz_z[k] = -g_x_0_z_0_xyzzz_z[k] * ab_y + g_x_0_z_0_xyzzz_yz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_z_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzzzz_x, g_x_0_z_0_xyzzzz_y, g_x_0_z_0_xyzzzz_z, g_x_0_z_0_xzzzz_x, g_x_0_z_0_xzzzz_xy, g_x_0_z_0_xzzzz_y, g_x_0_z_0_xzzzz_yy, g_x_0_z_0_xzzzz_yz, g_x_0_z_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzzzz_x[k] = -g_x_0_z_0_xzzzz_x[k] * ab_y + g_x_0_z_0_xzzzz_xy[k];

                g_x_0_z_0_xyzzzz_y[k] = -g_x_0_z_0_xzzzz_y[k] * ab_y + g_x_0_z_0_xzzzz_yy[k];

                g_x_0_z_0_xyzzzz_z[k] = -g_x_0_z_0_xzzzz_z[k] * ab_y + g_x_0_z_0_xzzzz_yz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_z_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzzzz_x, g_x_0_z_0_xzzzz_xz, g_x_0_z_0_xzzzz_y, g_x_0_z_0_xzzzz_yz, g_x_0_z_0_xzzzz_z, g_x_0_z_0_xzzzz_zz, g_x_0_z_0_xzzzzz_x, g_x_0_z_0_xzzzzz_y, g_x_0_z_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzzzz_x[k] = -g_x_0_z_0_xzzzz_x[k] * ab_z + g_x_0_z_0_xzzzz_xz[k];

                g_x_0_z_0_xzzzzz_y[k] = -g_x_0_z_0_xzzzz_y[k] * ab_z + g_x_0_z_0_xzzzz_yz[k];

                g_x_0_z_0_xzzzzz_z[k] = -g_x_0_z_0_xzzzz_z[k] * ab_z + g_x_0_z_0_xzzzz_zz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyy_x, g_x_0_z_0_yyyyy_xy, g_x_0_z_0_yyyyy_y, g_x_0_z_0_yyyyy_yy, g_x_0_z_0_yyyyy_yz, g_x_0_z_0_yyyyy_z, g_x_0_z_0_yyyyyy_x, g_x_0_z_0_yyyyyy_y, g_x_0_z_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyy_x[k] = -g_x_0_z_0_yyyyy_x[k] * ab_y + g_x_0_z_0_yyyyy_xy[k];

                g_x_0_z_0_yyyyyy_y[k] = -g_x_0_z_0_yyyyy_y[k] * ab_y + g_x_0_z_0_yyyyy_yy[k];

                g_x_0_z_0_yyyyyy_z[k] = -g_x_0_z_0_yyyyy_z[k] * ab_y + g_x_0_z_0_yyyyy_yz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyyz_x, g_x_0_z_0_yyyyyz_y, g_x_0_z_0_yyyyyz_z, g_x_0_z_0_yyyyz_x, g_x_0_z_0_yyyyz_xy, g_x_0_z_0_yyyyz_y, g_x_0_z_0_yyyyz_yy, g_x_0_z_0_yyyyz_yz, g_x_0_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyyz_x[k] = -g_x_0_z_0_yyyyz_x[k] * ab_y + g_x_0_z_0_yyyyz_xy[k];

                g_x_0_z_0_yyyyyz_y[k] = -g_x_0_z_0_yyyyz_y[k] * ab_y + g_x_0_z_0_yyyyz_yy[k];

                g_x_0_z_0_yyyyyz_z[k] = -g_x_0_z_0_yyyyz_z[k] * ab_y + g_x_0_z_0_yyyyz_yz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_z_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyyzz_x, g_x_0_z_0_yyyyzz_y, g_x_0_z_0_yyyyzz_z, g_x_0_z_0_yyyzz_x, g_x_0_z_0_yyyzz_xy, g_x_0_z_0_yyyzz_y, g_x_0_z_0_yyyzz_yy, g_x_0_z_0_yyyzz_yz, g_x_0_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyyzz_x[k] = -g_x_0_z_0_yyyzz_x[k] * ab_y + g_x_0_z_0_yyyzz_xy[k];

                g_x_0_z_0_yyyyzz_y[k] = -g_x_0_z_0_yyyzz_y[k] * ab_y + g_x_0_z_0_yyyzz_yy[k];

                g_x_0_z_0_yyyyzz_z[k] = -g_x_0_z_0_yyyzz_z[k] * ab_y + g_x_0_z_0_yyyzz_yz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_z_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyzzz_x, g_x_0_z_0_yyyzzz_y, g_x_0_z_0_yyyzzz_z, g_x_0_z_0_yyzzz_x, g_x_0_z_0_yyzzz_xy, g_x_0_z_0_yyzzz_y, g_x_0_z_0_yyzzz_yy, g_x_0_z_0_yyzzz_yz, g_x_0_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyzzz_x[k] = -g_x_0_z_0_yyzzz_x[k] * ab_y + g_x_0_z_0_yyzzz_xy[k];

                g_x_0_z_0_yyyzzz_y[k] = -g_x_0_z_0_yyzzz_y[k] * ab_y + g_x_0_z_0_yyzzz_yy[k];

                g_x_0_z_0_yyyzzz_z[k] = -g_x_0_z_0_yyzzz_z[k] * ab_y + g_x_0_z_0_yyzzz_yz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_z_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzzzz_x, g_x_0_z_0_yyzzzz_y, g_x_0_z_0_yyzzzz_z, g_x_0_z_0_yzzzz_x, g_x_0_z_0_yzzzz_xy, g_x_0_z_0_yzzzz_y, g_x_0_z_0_yzzzz_yy, g_x_0_z_0_yzzzz_yz, g_x_0_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzzzz_x[k] = -g_x_0_z_0_yzzzz_x[k] * ab_y + g_x_0_z_0_yzzzz_xy[k];

                g_x_0_z_0_yyzzzz_y[k] = -g_x_0_z_0_yzzzz_y[k] * ab_y + g_x_0_z_0_yzzzz_yy[k];

                g_x_0_z_0_yyzzzz_z[k] = -g_x_0_z_0_yzzzz_z[k] * ab_y + g_x_0_z_0_yzzzz_yz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_z_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzzzz_x, g_x_0_z_0_yzzzzz_y, g_x_0_z_0_yzzzzz_z, g_x_0_z_0_zzzzz_x, g_x_0_z_0_zzzzz_xy, g_x_0_z_0_zzzzz_y, g_x_0_z_0_zzzzz_yy, g_x_0_z_0_zzzzz_yz, g_x_0_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzzzz_x[k] = -g_x_0_z_0_zzzzz_x[k] * ab_y + g_x_0_z_0_zzzzz_xy[k];

                g_x_0_z_0_yzzzzz_y[k] = -g_x_0_z_0_zzzzz_y[k] * ab_y + g_x_0_z_0_zzzzz_yy[k];

                g_x_0_z_0_yzzzzz_z[k] = -g_x_0_z_0_zzzzz_z[k] * ab_y + g_x_0_z_0_zzzzz_yz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_z_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzzzz_x, g_x_0_z_0_zzzzz_xz, g_x_0_z_0_zzzzz_y, g_x_0_z_0_zzzzz_yz, g_x_0_z_0_zzzzz_z, g_x_0_z_0_zzzzz_zz, g_x_0_z_0_zzzzzz_x, g_x_0_z_0_zzzzzz_y, g_x_0_z_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzzzz_x[k] = -g_x_0_z_0_zzzzz_x[k] * ab_z + g_x_0_z_0_zzzzz_xz[k];

                g_x_0_z_0_zzzzzz_y[k] = -g_x_0_z_0_zzzzz_y[k] * ab_z + g_x_0_z_0_zzzzz_yz[k];

                g_x_0_z_0_zzzzzz_z[k] = -g_x_0_z_0_zzzzz_z[k] * ab_z + g_x_0_z_0_zzzzz_zz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxx_x, g_y_0_x_0_xxxxx_xx, g_y_0_x_0_xxxxx_xy, g_y_0_x_0_xxxxx_xz, g_y_0_x_0_xxxxx_y, g_y_0_x_0_xxxxx_z, g_y_0_x_0_xxxxxx_x, g_y_0_x_0_xxxxxx_y, g_y_0_x_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxx_x[k] = -g_y_0_x_0_xxxxx_x[k] * ab_x + g_y_0_x_0_xxxxx_xx[k];

                g_y_0_x_0_xxxxxx_y[k] = -g_y_0_x_0_xxxxx_y[k] * ab_x + g_y_0_x_0_xxxxx_xy[k];

                g_y_0_x_0_xxxxxx_z[k] = -g_y_0_x_0_xxxxx_z[k] * ab_x + g_y_0_x_0_xxxxx_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxy_x, g_y_0_x_0_xxxxxy_y, g_y_0_x_0_xxxxxy_z, g_y_0_x_0_xxxxy_x, g_y_0_x_0_xxxxy_xx, g_y_0_x_0_xxxxy_xy, g_y_0_x_0_xxxxy_xz, g_y_0_x_0_xxxxy_y, g_y_0_x_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxy_x[k] = -g_y_0_x_0_xxxxy_x[k] * ab_x + g_y_0_x_0_xxxxy_xx[k];

                g_y_0_x_0_xxxxxy_y[k] = -g_y_0_x_0_xxxxy_y[k] * ab_x + g_y_0_x_0_xxxxy_xy[k];

                g_y_0_x_0_xxxxxy_z[k] = -g_y_0_x_0_xxxxy_z[k] * ab_x + g_y_0_x_0_xxxxy_xz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxxz_x, g_y_0_x_0_xxxxxz_y, g_y_0_x_0_xxxxxz_z, g_y_0_x_0_xxxxz_x, g_y_0_x_0_xxxxz_xx, g_y_0_x_0_xxxxz_xy, g_y_0_x_0_xxxxz_xz, g_y_0_x_0_xxxxz_y, g_y_0_x_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxxz_x[k] = -g_y_0_x_0_xxxxz_x[k] * ab_x + g_y_0_x_0_xxxxz_xx[k];

                g_y_0_x_0_xxxxxz_y[k] = -g_y_0_x_0_xxxxz_y[k] * ab_x + g_y_0_x_0_xxxxz_xy[k];

                g_y_0_x_0_xxxxxz_z[k] = -g_y_0_x_0_xxxxz_z[k] * ab_x + g_y_0_x_0_xxxxz_xz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyy_x, g_y_0_x_0_xxxxyy_y, g_y_0_x_0_xxxxyy_z, g_y_0_x_0_xxxyy_x, g_y_0_x_0_xxxyy_xx, g_y_0_x_0_xxxyy_xy, g_y_0_x_0_xxxyy_xz, g_y_0_x_0_xxxyy_y, g_y_0_x_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyy_x[k] = -g_y_0_x_0_xxxyy_x[k] * ab_x + g_y_0_x_0_xxxyy_xx[k];

                g_y_0_x_0_xxxxyy_y[k] = -g_y_0_x_0_xxxyy_y[k] * ab_x + g_y_0_x_0_xxxyy_xy[k];

                g_y_0_x_0_xxxxyy_z[k] = -g_y_0_x_0_xxxyy_z[k] * ab_x + g_y_0_x_0_xxxyy_xz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxyz_x, g_y_0_x_0_xxxxyz_y, g_y_0_x_0_xxxxyz_z, g_y_0_x_0_xxxyz_x, g_y_0_x_0_xxxyz_xx, g_y_0_x_0_xxxyz_xy, g_y_0_x_0_xxxyz_xz, g_y_0_x_0_xxxyz_y, g_y_0_x_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxyz_x[k] = -g_y_0_x_0_xxxyz_x[k] * ab_x + g_y_0_x_0_xxxyz_xx[k];

                g_y_0_x_0_xxxxyz_y[k] = -g_y_0_x_0_xxxyz_y[k] * ab_x + g_y_0_x_0_xxxyz_xy[k];

                g_y_0_x_0_xxxxyz_z[k] = -g_y_0_x_0_xxxyz_z[k] * ab_x + g_y_0_x_0_xxxyz_xz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_x_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxxzz_x, g_y_0_x_0_xxxxzz_y, g_y_0_x_0_xxxxzz_z, g_y_0_x_0_xxxzz_x, g_y_0_x_0_xxxzz_xx, g_y_0_x_0_xxxzz_xy, g_y_0_x_0_xxxzz_xz, g_y_0_x_0_xxxzz_y, g_y_0_x_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxxzz_x[k] = -g_y_0_x_0_xxxzz_x[k] * ab_x + g_y_0_x_0_xxxzz_xx[k];

                g_y_0_x_0_xxxxzz_y[k] = -g_y_0_x_0_xxxzz_y[k] * ab_x + g_y_0_x_0_xxxzz_xy[k];

                g_y_0_x_0_xxxxzz_z[k] = -g_y_0_x_0_xxxzz_z[k] * ab_x + g_y_0_x_0_xxxzz_xz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 270 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 271 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyy_x, g_y_0_x_0_xxxyyy_y, g_y_0_x_0_xxxyyy_z, g_y_0_x_0_xxyyy_x, g_y_0_x_0_xxyyy_xx, g_y_0_x_0_xxyyy_xy, g_y_0_x_0_xxyyy_xz, g_y_0_x_0_xxyyy_y, g_y_0_x_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyy_x[k] = -g_y_0_x_0_xxyyy_x[k] * ab_x + g_y_0_x_0_xxyyy_xx[k];

                g_y_0_x_0_xxxyyy_y[k] = -g_y_0_x_0_xxyyy_y[k] * ab_x + g_y_0_x_0_xxyyy_xy[k];

                g_y_0_x_0_xxxyyy_z[k] = -g_y_0_x_0_xxyyy_z[k] * ab_x + g_y_0_x_0_xxyyy_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 273 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 274 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyyz_x, g_y_0_x_0_xxxyyz_y, g_y_0_x_0_xxxyyz_z, g_y_0_x_0_xxyyz_x, g_y_0_x_0_xxyyz_xx, g_y_0_x_0_xxyyz_xy, g_y_0_x_0_xxyyz_xz, g_y_0_x_0_xxyyz_y, g_y_0_x_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyyz_x[k] = -g_y_0_x_0_xxyyz_x[k] * ab_x + g_y_0_x_0_xxyyz_xx[k];

                g_y_0_x_0_xxxyyz_y[k] = -g_y_0_x_0_xxyyz_y[k] * ab_x + g_y_0_x_0_xxyyz_xy[k];

                g_y_0_x_0_xxxyyz_z[k] = -g_y_0_x_0_xxyyz_z[k] * ab_x + g_y_0_x_0_xxyyz_xz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 276 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 277 * ccomps * dcomps);

            auto g_y_0_x_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxyzz_x, g_y_0_x_0_xxxyzz_y, g_y_0_x_0_xxxyzz_z, g_y_0_x_0_xxyzz_x, g_y_0_x_0_xxyzz_xx, g_y_0_x_0_xxyzz_xy, g_y_0_x_0_xxyzz_xz, g_y_0_x_0_xxyzz_y, g_y_0_x_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxyzz_x[k] = -g_y_0_x_0_xxyzz_x[k] * ab_x + g_y_0_x_0_xxyzz_xx[k];

                g_y_0_x_0_xxxyzz_y[k] = -g_y_0_x_0_xxyzz_y[k] * ab_x + g_y_0_x_0_xxyzz_xy[k];

                g_y_0_x_0_xxxyzz_z[k] = -g_y_0_x_0_xxyzz_z[k] * ab_x + g_y_0_x_0_xxyzz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 279 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 280 * ccomps * dcomps);

            auto g_y_0_x_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxzzz_x, g_y_0_x_0_xxxzzz_y, g_y_0_x_0_xxxzzz_z, g_y_0_x_0_xxzzz_x, g_y_0_x_0_xxzzz_xx, g_y_0_x_0_xxzzz_xy, g_y_0_x_0_xxzzz_xz, g_y_0_x_0_xxzzz_y, g_y_0_x_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxzzz_x[k] = -g_y_0_x_0_xxzzz_x[k] * ab_x + g_y_0_x_0_xxzzz_xx[k];

                g_y_0_x_0_xxxzzz_y[k] = -g_y_0_x_0_xxzzz_y[k] * ab_x + g_y_0_x_0_xxzzz_xy[k];

                g_y_0_x_0_xxxzzz_z[k] = -g_y_0_x_0_xxzzz_z[k] * ab_x + g_y_0_x_0_xxzzz_xz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 282 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 283 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyy_x, g_y_0_x_0_xxyyyy_y, g_y_0_x_0_xxyyyy_z, g_y_0_x_0_xyyyy_x, g_y_0_x_0_xyyyy_xx, g_y_0_x_0_xyyyy_xy, g_y_0_x_0_xyyyy_xz, g_y_0_x_0_xyyyy_y, g_y_0_x_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyy_x[k] = -g_y_0_x_0_xyyyy_x[k] * ab_x + g_y_0_x_0_xyyyy_xx[k];

                g_y_0_x_0_xxyyyy_y[k] = -g_y_0_x_0_xyyyy_y[k] * ab_x + g_y_0_x_0_xyyyy_xy[k];

                g_y_0_x_0_xxyyyy_z[k] = -g_y_0_x_0_xyyyy_z[k] * ab_x + g_y_0_x_0_xyyyy_xz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 285 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 286 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyyz_x, g_y_0_x_0_xxyyyz_y, g_y_0_x_0_xxyyyz_z, g_y_0_x_0_xyyyz_x, g_y_0_x_0_xyyyz_xx, g_y_0_x_0_xyyyz_xy, g_y_0_x_0_xyyyz_xz, g_y_0_x_0_xyyyz_y, g_y_0_x_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyyz_x[k] = -g_y_0_x_0_xyyyz_x[k] * ab_x + g_y_0_x_0_xyyyz_xx[k];

                g_y_0_x_0_xxyyyz_y[k] = -g_y_0_x_0_xyyyz_y[k] * ab_x + g_y_0_x_0_xyyyz_xy[k];

                g_y_0_x_0_xxyyyz_z[k] = -g_y_0_x_0_xyyyz_z[k] * ab_x + g_y_0_x_0_xyyyz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 288 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 289 * ccomps * dcomps);

            auto g_y_0_x_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyyzz_x, g_y_0_x_0_xxyyzz_y, g_y_0_x_0_xxyyzz_z, g_y_0_x_0_xyyzz_x, g_y_0_x_0_xyyzz_xx, g_y_0_x_0_xyyzz_xy, g_y_0_x_0_xyyzz_xz, g_y_0_x_0_xyyzz_y, g_y_0_x_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyyzz_x[k] = -g_y_0_x_0_xyyzz_x[k] * ab_x + g_y_0_x_0_xyyzz_xx[k];

                g_y_0_x_0_xxyyzz_y[k] = -g_y_0_x_0_xyyzz_y[k] * ab_x + g_y_0_x_0_xyyzz_xy[k];

                g_y_0_x_0_xxyyzz_z[k] = -g_y_0_x_0_xyyzz_z[k] * ab_x + g_y_0_x_0_xyyzz_xz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 291 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 292 * ccomps * dcomps);

            auto g_y_0_x_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyzzz_x, g_y_0_x_0_xxyzzz_y, g_y_0_x_0_xxyzzz_z, g_y_0_x_0_xyzzz_x, g_y_0_x_0_xyzzz_xx, g_y_0_x_0_xyzzz_xy, g_y_0_x_0_xyzzz_xz, g_y_0_x_0_xyzzz_y, g_y_0_x_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyzzz_x[k] = -g_y_0_x_0_xyzzz_x[k] * ab_x + g_y_0_x_0_xyzzz_xx[k];

                g_y_0_x_0_xxyzzz_y[k] = -g_y_0_x_0_xyzzz_y[k] * ab_x + g_y_0_x_0_xyzzz_xy[k];

                g_y_0_x_0_xxyzzz_z[k] = -g_y_0_x_0_xyzzz_z[k] * ab_x + g_y_0_x_0_xyzzz_xz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 294 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 295 * ccomps * dcomps);

            auto g_y_0_x_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzzzz_x, g_y_0_x_0_xxzzzz_y, g_y_0_x_0_xxzzzz_z, g_y_0_x_0_xzzzz_x, g_y_0_x_0_xzzzz_xx, g_y_0_x_0_xzzzz_xy, g_y_0_x_0_xzzzz_xz, g_y_0_x_0_xzzzz_y, g_y_0_x_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzzzz_x[k] = -g_y_0_x_0_xzzzz_x[k] * ab_x + g_y_0_x_0_xzzzz_xx[k];

                g_y_0_x_0_xxzzzz_y[k] = -g_y_0_x_0_xzzzz_y[k] * ab_x + g_y_0_x_0_xzzzz_xy[k];

                g_y_0_x_0_xxzzzz_z[k] = -g_y_0_x_0_xzzzz_z[k] * ab_x + g_y_0_x_0_xzzzz_xz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 297 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 298 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyy_x, g_y_0_x_0_xyyyyy_y, g_y_0_x_0_xyyyyy_z, g_y_0_x_0_yyyyy_x, g_y_0_x_0_yyyyy_xx, g_y_0_x_0_yyyyy_xy, g_y_0_x_0_yyyyy_xz, g_y_0_x_0_yyyyy_y, g_y_0_x_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyy_x[k] = -g_y_0_x_0_yyyyy_x[k] * ab_x + g_y_0_x_0_yyyyy_xx[k];

                g_y_0_x_0_xyyyyy_y[k] = -g_y_0_x_0_yyyyy_y[k] * ab_x + g_y_0_x_0_yyyyy_xy[k];

                g_y_0_x_0_xyyyyy_z[k] = -g_y_0_x_0_yyyyy_z[k] * ab_x + g_y_0_x_0_yyyyy_xz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 300 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 301 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyyz_x, g_y_0_x_0_xyyyyz_y, g_y_0_x_0_xyyyyz_z, g_y_0_x_0_yyyyz_x, g_y_0_x_0_yyyyz_xx, g_y_0_x_0_yyyyz_xy, g_y_0_x_0_yyyyz_xz, g_y_0_x_0_yyyyz_y, g_y_0_x_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyyz_x[k] = -g_y_0_x_0_yyyyz_x[k] * ab_x + g_y_0_x_0_yyyyz_xx[k];

                g_y_0_x_0_xyyyyz_y[k] = -g_y_0_x_0_yyyyz_y[k] * ab_x + g_y_0_x_0_yyyyz_xy[k];

                g_y_0_x_0_xyyyyz_z[k] = -g_y_0_x_0_yyyyz_z[k] * ab_x + g_y_0_x_0_yyyyz_xz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 303 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 304 * ccomps * dcomps);

            auto g_y_0_x_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyyzz_x, g_y_0_x_0_xyyyzz_y, g_y_0_x_0_xyyyzz_z, g_y_0_x_0_yyyzz_x, g_y_0_x_0_yyyzz_xx, g_y_0_x_0_yyyzz_xy, g_y_0_x_0_yyyzz_xz, g_y_0_x_0_yyyzz_y, g_y_0_x_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyyzz_x[k] = -g_y_0_x_0_yyyzz_x[k] * ab_x + g_y_0_x_0_yyyzz_xx[k];

                g_y_0_x_0_xyyyzz_y[k] = -g_y_0_x_0_yyyzz_y[k] * ab_x + g_y_0_x_0_yyyzz_xy[k];

                g_y_0_x_0_xyyyzz_z[k] = -g_y_0_x_0_yyyzz_z[k] * ab_x + g_y_0_x_0_yyyzz_xz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 306 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 307 * ccomps * dcomps);

            auto g_y_0_x_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyzzz_x, g_y_0_x_0_xyyzzz_y, g_y_0_x_0_xyyzzz_z, g_y_0_x_0_yyzzz_x, g_y_0_x_0_yyzzz_xx, g_y_0_x_0_yyzzz_xy, g_y_0_x_0_yyzzz_xz, g_y_0_x_0_yyzzz_y, g_y_0_x_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyzzz_x[k] = -g_y_0_x_0_yyzzz_x[k] * ab_x + g_y_0_x_0_yyzzz_xx[k];

                g_y_0_x_0_xyyzzz_y[k] = -g_y_0_x_0_yyzzz_y[k] * ab_x + g_y_0_x_0_yyzzz_xy[k];

                g_y_0_x_0_xyyzzz_z[k] = -g_y_0_x_0_yyzzz_z[k] * ab_x + g_y_0_x_0_yyzzz_xz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 309 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 310 * ccomps * dcomps);

            auto g_y_0_x_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzzzz_x, g_y_0_x_0_xyzzzz_y, g_y_0_x_0_xyzzzz_z, g_y_0_x_0_yzzzz_x, g_y_0_x_0_yzzzz_xx, g_y_0_x_0_yzzzz_xy, g_y_0_x_0_yzzzz_xz, g_y_0_x_0_yzzzz_y, g_y_0_x_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzzzz_x[k] = -g_y_0_x_0_yzzzz_x[k] * ab_x + g_y_0_x_0_yzzzz_xx[k];

                g_y_0_x_0_xyzzzz_y[k] = -g_y_0_x_0_yzzzz_y[k] * ab_x + g_y_0_x_0_yzzzz_xy[k];

                g_y_0_x_0_xyzzzz_z[k] = -g_y_0_x_0_yzzzz_z[k] * ab_x + g_y_0_x_0_yzzzz_xz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 312 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 313 * ccomps * dcomps);

            auto g_y_0_x_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzzzz_x, g_y_0_x_0_xzzzzz_y, g_y_0_x_0_xzzzzz_z, g_y_0_x_0_zzzzz_x, g_y_0_x_0_zzzzz_xx, g_y_0_x_0_zzzzz_xy, g_y_0_x_0_zzzzz_xz, g_y_0_x_0_zzzzz_y, g_y_0_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzzzz_x[k] = -g_y_0_x_0_zzzzz_x[k] * ab_x + g_y_0_x_0_zzzzz_xx[k];

                g_y_0_x_0_xzzzzz_y[k] = -g_y_0_x_0_zzzzz_y[k] * ab_x + g_y_0_x_0_zzzzz_xy[k];

                g_y_0_x_0_xzzzzz_z[k] = -g_y_0_x_0_zzzzz_z[k] * ab_x + g_y_0_x_0_zzzzz_xz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 315 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 316 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyyyy_x, g_0_0_x_0_yyyyy_y, g_0_0_x_0_yyyyy_z, g_y_0_x_0_yyyyy_x, g_y_0_x_0_yyyyy_xy, g_y_0_x_0_yyyyy_y, g_y_0_x_0_yyyyy_yy, g_y_0_x_0_yyyyy_yz, g_y_0_x_0_yyyyy_z, g_y_0_x_0_yyyyyy_x, g_y_0_x_0_yyyyyy_y, g_y_0_x_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyy_x[k] = -g_0_0_x_0_yyyyy_x[k] - g_y_0_x_0_yyyyy_x[k] * ab_y + g_y_0_x_0_yyyyy_xy[k];

                g_y_0_x_0_yyyyyy_y[k] = -g_0_0_x_0_yyyyy_y[k] - g_y_0_x_0_yyyyy_y[k] * ab_y + g_y_0_x_0_yyyyy_yy[k];

                g_y_0_x_0_yyyyyy_z[k] = -g_0_0_x_0_yyyyy_z[k] - g_y_0_x_0_yyyyy_z[k] * ab_y + g_y_0_x_0_yyyyy_yz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 318 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 319 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyy_x, g_y_0_x_0_yyyyy_xz, g_y_0_x_0_yyyyy_y, g_y_0_x_0_yyyyy_yz, g_y_0_x_0_yyyyy_z, g_y_0_x_0_yyyyy_zz, g_y_0_x_0_yyyyyz_x, g_y_0_x_0_yyyyyz_y, g_y_0_x_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyyz_x[k] = -g_y_0_x_0_yyyyy_x[k] * ab_z + g_y_0_x_0_yyyyy_xz[k];

                g_y_0_x_0_yyyyyz_y[k] = -g_y_0_x_0_yyyyy_y[k] * ab_z + g_y_0_x_0_yyyyy_yz[k];

                g_y_0_x_0_yyyyyz_z[k] = -g_y_0_x_0_yyyyy_z[k] * ab_z + g_y_0_x_0_yyyyy_zz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 321 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 322 * ccomps * dcomps);

            auto g_y_0_x_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyyz_x, g_y_0_x_0_yyyyz_xz, g_y_0_x_0_yyyyz_y, g_y_0_x_0_yyyyz_yz, g_y_0_x_0_yyyyz_z, g_y_0_x_0_yyyyz_zz, g_y_0_x_0_yyyyzz_x, g_y_0_x_0_yyyyzz_y, g_y_0_x_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyyzz_x[k] = -g_y_0_x_0_yyyyz_x[k] * ab_z + g_y_0_x_0_yyyyz_xz[k];

                g_y_0_x_0_yyyyzz_y[k] = -g_y_0_x_0_yyyyz_y[k] * ab_z + g_y_0_x_0_yyyyz_yz[k];

                g_y_0_x_0_yyyyzz_z[k] = -g_y_0_x_0_yyyyz_z[k] * ab_z + g_y_0_x_0_yyyyz_zz[k];
            }

            /// Set up 324-327 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 324 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 325 * ccomps * dcomps);

            auto g_y_0_x_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 326 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyyzz_x, g_y_0_x_0_yyyzz_xz, g_y_0_x_0_yyyzz_y, g_y_0_x_0_yyyzz_yz, g_y_0_x_0_yyyzz_z, g_y_0_x_0_yyyzz_zz, g_y_0_x_0_yyyzzz_x, g_y_0_x_0_yyyzzz_y, g_y_0_x_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyzzz_x[k] = -g_y_0_x_0_yyyzz_x[k] * ab_z + g_y_0_x_0_yyyzz_xz[k];

                g_y_0_x_0_yyyzzz_y[k] = -g_y_0_x_0_yyyzz_y[k] * ab_z + g_y_0_x_0_yyyzz_yz[k];

                g_y_0_x_0_yyyzzz_z[k] = -g_y_0_x_0_yyyzz_z[k] * ab_z + g_y_0_x_0_yyyzz_zz[k];
            }

            /// Set up 327-330 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 327 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 328 * ccomps * dcomps);

            auto g_y_0_x_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyzzz_x, g_y_0_x_0_yyzzz_xz, g_y_0_x_0_yyzzz_y, g_y_0_x_0_yyzzz_yz, g_y_0_x_0_yyzzz_z, g_y_0_x_0_yyzzz_zz, g_y_0_x_0_yyzzzz_x, g_y_0_x_0_yyzzzz_y, g_y_0_x_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzzzz_x[k] = -g_y_0_x_0_yyzzz_x[k] * ab_z + g_y_0_x_0_yyzzz_xz[k];

                g_y_0_x_0_yyzzzz_y[k] = -g_y_0_x_0_yyzzz_y[k] * ab_z + g_y_0_x_0_yyzzz_yz[k];

                g_y_0_x_0_yyzzzz_z[k] = -g_y_0_x_0_yyzzz_z[k] * ab_z + g_y_0_x_0_yyzzz_zz[k];
            }

            /// Set up 330-333 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 330 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 331 * ccomps * dcomps);

            auto g_y_0_x_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 332 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzzzz_x, g_y_0_x_0_yzzzz_xz, g_y_0_x_0_yzzzz_y, g_y_0_x_0_yzzzz_yz, g_y_0_x_0_yzzzz_z, g_y_0_x_0_yzzzz_zz, g_y_0_x_0_yzzzzz_x, g_y_0_x_0_yzzzzz_y, g_y_0_x_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzzzz_x[k] = -g_y_0_x_0_yzzzz_x[k] * ab_z + g_y_0_x_0_yzzzz_xz[k];

                g_y_0_x_0_yzzzzz_y[k] = -g_y_0_x_0_yzzzz_y[k] * ab_z + g_y_0_x_0_yzzzz_yz[k];

                g_y_0_x_0_yzzzzz_z[k] = -g_y_0_x_0_yzzzz_z[k] * ab_z + g_y_0_x_0_yzzzz_zz[k];
            }

            /// Set up 333-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 333 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 334 * ccomps * dcomps);

            auto g_y_0_x_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzzzz_x, g_y_0_x_0_zzzzz_xz, g_y_0_x_0_zzzzz_y, g_y_0_x_0_zzzzz_yz, g_y_0_x_0_zzzzz_z, g_y_0_x_0_zzzzz_zz, g_y_0_x_0_zzzzzz_x, g_y_0_x_0_zzzzzz_y, g_y_0_x_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzzzz_x[k] = -g_y_0_x_0_zzzzz_x[k] * ab_z + g_y_0_x_0_zzzzz_xz[k];

                g_y_0_x_0_zzzzzz_y[k] = -g_y_0_x_0_zzzzz_y[k] * ab_z + g_y_0_x_0_zzzzz_yz[k];

                g_y_0_x_0_zzzzzz_z[k] = -g_y_0_x_0_zzzzz_z[k] * ab_z + g_y_0_x_0_zzzzz_zz[k];
            }

            /// Set up 336-339 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 336 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 337 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 338 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxx_x, g_y_0_y_0_xxxxx_xx, g_y_0_y_0_xxxxx_xy, g_y_0_y_0_xxxxx_xz, g_y_0_y_0_xxxxx_y, g_y_0_y_0_xxxxx_z, g_y_0_y_0_xxxxxx_x, g_y_0_y_0_xxxxxx_y, g_y_0_y_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxx_x[k] = -g_y_0_y_0_xxxxx_x[k] * ab_x + g_y_0_y_0_xxxxx_xx[k];

                g_y_0_y_0_xxxxxx_y[k] = -g_y_0_y_0_xxxxx_y[k] * ab_x + g_y_0_y_0_xxxxx_xy[k];

                g_y_0_y_0_xxxxxx_z[k] = -g_y_0_y_0_xxxxx_z[k] * ab_x + g_y_0_y_0_xxxxx_xz[k];
            }

            /// Set up 339-342 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 339 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 340 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxy_x, g_y_0_y_0_xxxxxy_y, g_y_0_y_0_xxxxxy_z, g_y_0_y_0_xxxxy_x, g_y_0_y_0_xxxxy_xx, g_y_0_y_0_xxxxy_xy, g_y_0_y_0_xxxxy_xz, g_y_0_y_0_xxxxy_y, g_y_0_y_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxy_x[k] = -g_y_0_y_0_xxxxy_x[k] * ab_x + g_y_0_y_0_xxxxy_xx[k];

                g_y_0_y_0_xxxxxy_y[k] = -g_y_0_y_0_xxxxy_y[k] * ab_x + g_y_0_y_0_xxxxy_xy[k];

                g_y_0_y_0_xxxxxy_z[k] = -g_y_0_y_0_xxxxy_z[k] * ab_x + g_y_0_y_0_xxxxy_xz[k];
            }

            /// Set up 342-345 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 342 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 343 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxxz_x, g_y_0_y_0_xxxxxz_y, g_y_0_y_0_xxxxxz_z, g_y_0_y_0_xxxxz_x, g_y_0_y_0_xxxxz_xx, g_y_0_y_0_xxxxz_xy, g_y_0_y_0_xxxxz_xz, g_y_0_y_0_xxxxz_y, g_y_0_y_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxxz_x[k] = -g_y_0_y_0_xxxxz_x[k] * ab_x + g_y_0_y_0_xxxxz_xx[k];

                g_y_0_y_0_xxxxxz_y[k] = -g_y_0_y_0_xxxxz_y[k] * ab_x + g_y_0_y_0_xxxxz_xy[k];

                g_y_0_y_0_xxxxxz_z[k] = -g_y_0_y_0_xxxxz_z[k] * ab_x + g_y_0_y_0_xxxxz_xz[k];
            }

            /// Set up 345-348 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 345 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 346 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyy_x, g_y_0_y_0_xxxxyy_y, g_y_0_y_0_xxxxyy_z, g_y_0_y_0_xxxyy_x, g_y_0_y_0_xxxyy_xx, g_y_0_y_0_xxxyy_xy, g_y_0_y_0_xxxyy_xz, g_y_0_y_0_xxxyy_y, g_y_0_y_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyy_x[k] = -g_y_0_y_0_xxxyy_x[k] * ab_x + g_y_0_y_0_xxxyy_xx[k];

                g_y_0_y_0_xxxxyy_y[k] = -g_y_0_y_0_xxxyy_y[k] * ab_x + g_y_0_y_0_xxxyy_xy[k];

                g_y_0_y_0_xxxxyy_z[k] = -g_y_0_y_0_xxxyy_z[k] * ab_x + g_y_0_y_0_xxxyy_xz[k];
            }

            /// Set up 348-351 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 348 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 349 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 350 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxyz_x, g_y_0_y_0_xxxxyz_y, g_y_0_y_0_xxxxyz_z, g_y_0_y_0_xxxyz_x, g_y_0_y_0_xxxyz_xx, g_y_0_y_0_xxxyz_xy, g_y_0_y_0_xxxyz_xz, g_y_0_y_0_xxxyz_y, g_y_0_y_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxyz_x[k] = -g_y_0_y_0_xxxyz_x[k] * ab_x + g_y_0_y_0_xxxyz_xx[k];

                g_y_0_y_0_xxxxyz_y[k] = -g_y_0_y_0_xxxyz_y[k] * ab_x + g_y_0_y_0_xxxyz_xy[k];

                g_y_0_y_0_xxxxyz_z[k] = -g_y_0_y_0_xxxyz_z[k] * ab_x + g_y_0_y_0_xxxyz_xz[k];
            }

            /// Set up 351-354 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 351 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 352 * ccomps * dcomps);

            auto g_y_0_y_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxxzz_x, g_y_0_y_0_xxxxzz_y, g_y_0_y_0_xxxxzz_z, g_y_0_y_0_xxxzz_x, g_y_0_y_0_xxxzz_xx, g_y_0_y_0_xxxzz_xy, g_y_0_y_0_xxxzz_xz, g_y_0_y_0_xxxzz_y, g_y_0_y_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxxzz_x[k] = -g_y_0_y_0_xxxzz_x[k] * ab_x + g_y_0_y_0_xxxzz_xx[k];

                g_y_0_y_0_xxxxzz_y[k] = -g_y_0_y_0_xxxzz_y[k] * ab_x + g_y_0_y_0_xxxzz_xy[k];

                g_y_0_y_0_xxxxzz_z[k] = -g_y_0_y_0_xxxzz_z[k] * ab_x + g_y_0_y_0_xxxzz_xz[k];
            }

            /// Set up 354-357 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 354 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 355 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyy_x, g_y_0_y_0_xxxyyy_y, g_y_0_y_0_xxxyyy_z, g_y_0_y_0_xxyyy_x, g_y_0_y_0_xxyyy_xx, g_y_0_y_0_xxyyy_xy, g_y_0_y_0_xxyyy_xz, g_y_0_y_0_xxyyy_y, g_y_0_y_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyy_x[k] = -g_y_0_y_0_xxyyy_x[k] * ab_x + g_y_0_y_0_xxyyy_xx[k];

                g_y_0_y_0_xxxyyy_y[k] = -g_y_0_y_0_xxyyy_y[k] * ab_x + g_y_0_y_0_xxyyy_xy[k];

                g_y_0_y_0_xxxyyy_z[k] = -g_y_0_y_0_xxyyy_z[k] * ab_x + g_y_0_y_0_xxyyy_xz[k];
            }

            /// Set up 357-360 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 357 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 358 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyyz_x, g_y_0_y_0_xxxyyz_y, g_y_0_y_0_xxxyyz_z, g_y_0_y_0_xxyyz_x, g_y_0_y_0_xxyyz_xx, g_y_0_y_0_xxyyz_xy, g_y_0_y_0_xxyyz_xz, g_y_0_y_0_xxyyz_y, g_y_0_y_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyyz_x[k] = -g_y_0_y_0_xxyyz_x[k] * ab_x + g_y_0_y_0_xxyyz_xx[k];

                g_y_0_y_0_xxxyyz_y[k] = -g_y_0_y_0_xxyyz_y[k] * ab_x + g_y_0_y_0_xxyyz_xy[k];

                g_y_0_y_0_xxxyyz_z[k] = -g_y_0_y_0_xxyyz_z[k] * ab_x + g_y_0_y_0_xxyyz_xz[k];
            }

            /// Set up 360-363 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 360 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 361 * ccomps * dcomps);

            auto g_y_0_y_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 362 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxyzz_x, g_y_0_y_0_xxxyzz_y, g_y_0_y_0_xxxyzz_z, g_y_0_y_0_xxyzz_x, g_y_0_y_0_xxyzz_xx, g_y_0_y_0_xxyzz_xy, g_y_0_y_0_xxyzz_xz, g_y_0_y_0_xxyzz_y, g_y_0_y_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxyzz_x[k] = -g_y_0_y_0_xxyzz_x[k] * ab_x + g_y_0_y_0_xxyzz_xx[k];

                g_y_0_y_0_xxxyzz_y[k] = -g_y_0_y_0_xxyzz_y[k] * ab_x + g_y_0_y_0_xxyzz_xy[k];

                g_y_0_y_0_xxxyzz_z[k] = -g_y_0_y_0_xxyzz_z[k] * ab_x + g_y_0_y_0_xxyzz_xz[k];
            }

            /// Set up 363-366 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 363 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 364 * ccomps * dcomps);

            auto g_y_0_y_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxzzz_x, g_y_0_y_0_xxxzzz_y, g_y_0_y_0_xxxzzz_z, g_y_0_y_0_xxzzz_x, g_y_0_y_0_xxzzz_xx, g_y_0_y_0_xxzzz_xy, g_y_0_y_0_xxzzz_xz, g_y_0_y_0_xxzzz_y, g_y_0_y_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxzzz_x[k] = -g_y_0_y_0_xxzzz_x[k] * ab_x + g_y_0_y_0_xxzzz_xx[k];

                g_y_0_y_0_xxxzzz_y[k] = -g_y_0_y_0_xxzzz_y[k] * ab_x + g_y_0_y_0_xxzzz_xy[k];

                g_y_0_y_0_xxxzzz_z[k] = -g_y_0_y_0_xxzzz_z[k] * ab_x + g_y_0_y_0_xxzzz_xz[k];
            }

            /// Set up 366-369 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 366 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 367 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 368 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyy_x, g_y_0_y_0_xxyyyy_y, g_y_0_y_0_xxyyyy_z, g_y_0_y_0_xyyyy_x, g_y_0_y_0_xyyyy_xx, g_y_0_y_0_xyyyy_xy, g_y_0_y_0_xyyyy_xz, g_y_0_y_0_xyyyy_y, g_y_0_y_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyy_x[k] = -g_y_0_y_0_xyyyy_x[k] * ab_x + g_y_0_y_0_xyyyy_xx[k];

                g_y_0_y_0_xxyyyy_y[k] = -g_y_0_y_0_xyyyy_y[k] * ab_x + g_y_0_y_0_xyyyy_xy[k];

                g_y_0_y_0_xxyyyy_z[k] = -g_y_0_y_0_xyyyy_z[k] * ab_x + g_y_0_y_0_xyyyy_xz[k];
            }

            /// Set up 369-372 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 369 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 370 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyyz_x, g_y_0_y_0_xxyyyz_y, g_y_0_y_0_xxyyyz_z, g_y_0_y_0_xyyyz_x, g_y_0_y_0_xyyyz_xx, g_y_0_y_0_xyyyz_xy, g_y_0_y_0_xyyyz_xz, g_y_0_y_0_xyyyz_y, g_y_0_y_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyyz_x[k] = -g_y_0_y_0_xyyyz_x[k] * ab_x + g_y_0_y_0_xyyyz_xx[k];

                g_y_0_y_0_xxyyyz_y[k] = -g_y_0_y_0_xyyyz_y[k] * ab_x + g_y_0_y_0_xyyyz_xy[k];

                g_y_0_y_0_xxyyyz_z[k] = -g_y_0_y_0_xyyyz_z[k] * ab_x + g_y_0_y_0_xyyyz_xz[k];
            }

            /// Set up 372-375 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 372 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 373 * ccomps * dcomps);

            auto g_y_0_y_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyyzz_x, g_y_0_y_0_xxyyzz_y, g_y_0_y_0_xxyyzz_z, g_y_0_y_0_xyyzz_x, g_y_0_y_0_xyyzz_xx, g_y_0_y_0_xyyzz_xy, g_y_0_y_0_xyyzz_xz, g_y_0_y_0_xyyzz_y, g_y_0_y_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyyzz_x[k] = -g_y_0_y_0_xyyzz_x[k] * ab_x + g_y_0_y_0_xyyzz_xx[k];

                g_y_0_y_0_xxyyzz_y[k] = -g_y_0_y_0_xyyzz_y[k] * ab_x + g_y_0_y_0_xyyzz_xy[k];

                g_y_0_y_0_xxyyzz_z[k] = -g_y_0_y_0_xyyzz_z[k] * ab_x + g_y_0_y_0_xyyzz_xz[k];
            }

            /// Set up 375-378 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 375 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 376 * ccomps * dcomps);

            auto g_y_0_y_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyzzz_x, g_y_0_y_0_xxyzzz_y, g_y_0_y_0_xxyzzz_z, g_y_0_y_0_xyzzz_x, g_y_0_y_0_xyzzz_xx, g_y_0_y_0_xyzzz_xy, g_y_0_y_0_xyzzz_xz, g_y_0_y_0_xyzzz_y, g_y_0_y_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyzzz_x[k] = -g_y_0_y_0_xyzzz_x[k] * ab_x + g_y_0_y_0_xyzzz_xx[k];

                g_y_0_y_0_xxyzzz_y[k] = -g_y_0_y_0_xyzzz_y[k] * ab_x + g_y_0_y_0_xyzzz_xy[k];

                g_y_0_y_0_xxyzzz_z[k] = -g_y_0_y_0_xyzzz_z[k] * ab_x + g_y_0_y_0_xyzzz_xz[k];
            }

            /// Set up 378-381 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 378 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 379 * ccomps * dcomps);

            auto g_y_0_y_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 380 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzzzz_x, g_y_0_y_0_xxzzzz_y, g_y_0_y_0_xxzzzz_z, g_y_0_y_0_xzzzz_x, g_y_0_y_0_xzzzz_xx, g_y_0_y_0_xzzzz_xy, g_y_0_y_0_xzzzz_xz, g_y_0_y_0_xzzzz_y, g_y_0_y_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzzzz_x[k] = -g_y_0_y_0_xzzzz_x[k] * ab_x + g_y_0_y_0_xzzzz_xx[k];

                g_y_0_y_0_xxzzzz_y[k] = -g_y_0_y_0_xzzzz_y[k] * ab_x + g_y_0_y_0_xzzzz_xy[k];

                g_y_0_y_0_xxzzzz_z[k] = -g_y_0_y_0_xzzzz_z[k] * ab_x + g_y_0_y_0_xzzzz_xz[k];
            }

            /// Set up 381-384 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 381 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 382 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyy_x, g_y_0_y_0_xyyyyy_y, g_y_0_y_0_xyyyyy_z, g_y_0_y_0_yyyyy_x, g_y_0_y_0_yyyyy_xx, g_y_0_y_0_yyyyy_xy, g_y_0_y_0_yyyyy_xz, g_y_0_y_0_yyyyy_y, g_y_0_y_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyy_x[k] = -g_y_0_y_0_yyyyy_x[k] * ab_x + g_y_0_y_0_yyyyy_xx[k];

                g_y_0_y_0_xyyyyy_y[k] = -g_y_0_y_0_yyyyy_y[k] * ab_x + g_y_0_y_0_yyyyy_xy[k];

                g_y_0_y_0_xyyyyy_z[k] = -g_y_0_y_0_yyyyy_z[k] * ab_x + g_y_0_y_0_yyyyy_xz[k];
            }

            /// Set up 384-387 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 384 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 385 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 386 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyyz_x, g_y_0_y_0_xyyyyz_y, g_y_0_y_0_xyyyyz_z, g_y_0_y_0_yyyyz_x, g_y_0_y_0_yyyyz_xx, g_y_0_y_0_yyyyz_xy, g_y_0_y_0_yyyyz_xz, g_y_0_y_0_yyyyz_y, g_y_0_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyyz_x[k] = -g_y_0_y_0_yyyyz_x[k] * ab_x + g_y_0_y_0_yyyyz_xx[k];

                g_y_0_y_0_xyyyyz_y[k] = -g_y_0_y_0_yyyyz_y[k] * ab_x + g_y_0_y_0_yyyyz_xy[k];

                g_y_0_y_0_xyyyyz_z[k] = -g_y_0_y_0_yyyyz_z[k] * ab_x + g_y_0_y_0_yyyyz_xz[k];
            }

            /// Set up 387-390 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 387 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 388 * ccomps * dcomps);

            auto g_y_0_y_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyyzz_x, g_y_0_y_0_xyyyzz_y, g_y_0_y_0_xyyyzz_z, g_y_0_y_0_yyyzz_x, g_y_0_y_0_yyyzz_xx, g_y_0_y_0_yyyzz_xy, g_y_0_y_0_yyyzz_xz, g_y_0_y_0_yyyzz_y, g_y_0_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyyzz_x[k] = -g_y_0_y_0_yyyzz_x[k] * ab_x + g_y_0_y_0_yyyzz_xx[k];

                g_y_0_y_0_xyyyzz_y[k] = -g_y_0_y_0_yyyzz_y[k] * ab_x + g_y_0_y_0_yyyzz_xy[k];

                g_y_0_y_0_xyyyzz_z[k] = -g_y_0_y_0_yyyzz_z[k] * ab_x + g_y_0_y_0_yyyzz_xz[k];
            }

            /// Set up 390-393 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 390 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 391 * ccomps * dcomps);

            auto g_y_0_y_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 392 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyzzz_x, g_y_0_y_0_xyyzzz_y, g_y_0_y_0_xyyzzz_z, g_y_0_y_0_yyzzz_x, g_y_0_y_0_yyzzz_xx, g_y_0_y_0_yyzzz_xy, g_y_0_y_0_yyzzz_xz, g_y_0_y_0_yyzzz_y, g_y_0_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyzzz_x[k] = -g_y_0_y_0_yyzzz_x[k] * ab_x + g_y_0_y_0_yyzzz_xx[k];

                g_y_0_y_0_xyyzzz_y[k] = -g_y_0_y_0_yyzzz_y[k] * ab_x + g_y_0_y_0_yyzzz_xy[k];

                g_y_0_y_0_xyyzzz_z[k] = -g_y_0_y_0_yyzzz_z[k] * ab_x + g_y_0_y_0_yyzzz_xz[k];
            }

            /// Set up 393-396 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 393 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 394 * ccomps * dcomps);

            auto g_y_0_y_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzzzz_x, g_y_0_y_0_xyzzzz_y, g_y_0_y_0_xyzzzz_z, g_y_0_y_0_yzzzz_x, g_y_0_y_0_yzzzz_xx, g_y_0_y_0_yzzzz_xy, g_y_0_y_0_yzzzz_xz, g_y_0_y_0_yzzzz_y, g_y_0_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzzzz_x[k] = -g_y_0_y_0_yzzzz_x[k] * ab_x + g_y_0_y_0_yzzzz_xx[k];

                g_y_0_y_0_xyzzzz_y[k] = -g_y_0_y_0_yzzzz_y[k] * ab_x + g_y_0_y_0_yzzzz_xy[k];

                g_y_0_y_0_xyzzzz_z[k] = -g_y_0_y_0_yzzzz_z[k] * ab_x + g_y_0_y_0_yzzzz_xz[k];
            }

            /// Set up 396-399 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 396 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 397 * ccomps * dcomps);

            auto g_y_0_y_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzzzz_x, g_y_0_y_0_xzzzzz_y, g_y_0_y_0_xzzzzz_z, g_y_0_y_0_zzzzz_x, g_y_0_y_0_zzzzz_xx, g_y_0_y_0_zzzzz_xy, g_y_0_y_0_zzzzz_xz, g_y_0_y_0_zzzzz_y, g_y_0_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzzzz_x[k] = -g_y_0_y_0_zzzzz_x[k] * ab_x + g_y_0_y_0_zzzzz_xx[k];

                g_y_0_y_0_xzzzzz_y[k] = -g_y_0_y_0_zzzzz_y[k] * ab_x + g_y_0_y_0_zzzzz_xy[k];

                g_y_0_y_0_xzzzzz_z[k] = -g_y_0_y_0_zzzzz_z[k] * ab_x + g_y_0_y_0_zzzzz_xz[k];
            }

            /// Set up 399-402 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 399 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 400 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyyyy_x, g_0_0_y_0_yyyyy_y, g_0_0_y_0_yyyyy_z, g_y_0_y_0_yyyyy_x, g_y_0_y_0_yyyyy_xy, g_y_0_y_0_yyyyy_y, g_y_0_y_0_yyyyy_yy, g_y_0_y_0_yyyyy_yz, g_y_0_y_0_yyyyy_z, g_y_0_y_0_yyyyyy_x, g_y_0_y_0_yyyyyy_y, g_y_0_y_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyy_x[k] = -g_0_0_y_0_yyyyy_x[k] - g_y_0_y_0_yyyyy_x[k] * ab_y + g_y_0_y_0_yyyyy_xy[k];

                g_y_0_y_0_yyyyyy_y[k] = -g_0_0_y_0_yyyyy_y[k] - g_y_0_y_0_yyyyy_y[k] * ab_y + g_y_0_y_0_yyyyy_yy[k];

                g_y_0_y_0_yyyyyy_z[k] = -g_0_0_y_0_yyyyy_z[k] - g_y_0_y_0_yyyyy_z[k] * ab_y + g_y_0_y_0_yyyyy_yz[k];
            }

            /// Set up 402-405 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 402 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 403 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyy_x, g_y_0_y_0_yyyyy_xz, g_y_0_y_0_yyyyy_y, g_y_0_y_0_yyyyy_yz, g_y_0_y_0_yyyyy_z, g_y_0_y_0_yyyyy_zz, g_y_0_y_0_yyyyyz_x, g_y_0_y_0_yyyyyz_y, g_y_0_y_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyyz_x[k] = -g_y_0_y_0_yyyyy_x[k] * ab_z + g_y_0_y_0_yyyyy_xz[k];

                g_y_0_y_0_yyyyyz_y[k] = -g_y_0_y_0_yyyyy_y[k] * ab_z + g_y_0_y_0_yyyyy_yz[k];

                g_y_0_y_0_yyyyyz_z[k] = -g_y_0_y_0_yyyyy_z[k] * ab_z + g_y_0_y_0_yyyyy_zz[k];
            }

            /// Set up 405-408 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 405 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 406 * ccomps * dcomps);

            auto g_y_0_y_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyyz_x, g_y_0_y_0_yyyyz_xz, g_y_0_y_0_yyyyz_y, g_y_0_y_0_yyyyz_yz, g_y_0_y_0_yyyyz_z, g_y_0_y_0_yyyyz_zz, g_y_0_y_0_yyyyzz_x, g_y_0_y_0_yyyyzz_y, g_y_0_y_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyyzz_x[k] = -g_y_0_y_0_yyyyz_x[k] * ab_z + g_y_0_y_0_yyyyz_xz[k];

                g_y_0_y_0_yyyyzz_y[k] = -g_y_0_y_0_yyyyz_y[k] * ab_z + g_y_0_y_0_yyyyz_yz[k];

                g_y_0_y_0_yyyyzz_z[k] = -g_y_0_y_0_yyyyz_z[k] * ab_z + g_y_0_y_0_yyyyz_zz[k];
            }

            /// Set up 408-411 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 408 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 409 * ccomps * dcomps);

            auto g_y_0_y_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 410 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyyzz_x, g_y_0_y_0_yyyzz_xz, g_y_0_y_0_yyyzz_y, g_y_0_y_0_yyyzz_yz, g_y_0_y_0_yyyzz_z, g_y_0_y_0_yyyzz_zz, g_y_0_y_0_yyyzzz_x, g_y_0_y_0_yyyzzz_y, g_y_0_y_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyzzz_x[k] = -g_y_0_y_0_yyyzz_x[k] * ab_z + g_y_0_y_0_yyyzz_xz[k];

                g_y_0_y_0_yyyzzz_y[k] = -g_y_0_y_0_yyyzz_y[k] * ab_z + g_y_0_y_0_yyyzz_yz[k];

                g_y_0_y_0_yyyzzz_z[k] = -g_y_0_y_0_yyyzz_z[k] * ab_z + g_y_0_y_0_yyyzz_zz[k];
            }

            /// Set up 411-414 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 411 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 412 * ccomps * dcomps);

            auto g_y_0_y_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyzzz_x, g_y_0_y_0_yyzzz_xz, g_y_0_y_0_yyzzz_y, g_y_0_y_0_yyzzz_yz, g_y_0_y_0_yyzzz_z, g_y_0_y_0_yyzzz_zz, g_y_0_y_0_yyzzzz_x, g_y_0_y_0_yyzzzz_y, g_y_0_y_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzzzz_x[k] = -g_y_0_y_0_yyzzz_x[k] * ab_z + g_y_0_y_0_yyzzz_xz[k];

                g_y_0_y_0_yyzzzz_y[k] = -g_y_0_y_0_yyzzz_y[k] * ab_z + g_y_0_y_0_yyzzz_yz[k];

                g_y_0_y_0_yyzzzz_z[k] = -g_y_0_y_0_yyzzz_z[k] * ab_z + g_y_0_y_0_yyzzz_zz[k];
            }

            /// Set up 414-417 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 414 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 415 * ccomps * dcomps);

            auto g_y_0_y_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 416 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzzzz_x, g_y_0_y_0_yzzzz_xz, g_y_0_y_0_yzzzz_y, g_y_0_y_0_yzzzz_yz, g_y_0_y_0_yzzzz_z, g_y_0_y_0_yzzzz_zz, g_y_0_y_0_yzzzzz_x, g_y_0_y_0_yzzzzz_y, g_y_0_y_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzzzz_x[k] = -g_y_0_y_0_yzzzz_x[k] * ab_z + g_y_0_y_0_yzzzz_xz[k];

                g_y_0_y_0_yzzzzz_y[k] = -g_y_0_y_0_yzzzz_y[k] * ab_z + g_y_0_y_0_yzzzz_yz[k];

                g_y_0_y_0_yzzzzz_z[k] = -g_y_0_y_0_yzzzz_z[k] * ab_z + g_y_0_y_0_yzzzz_zz[k];
            }

            /// Set up 417-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 417 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 418 * ccomps * dcomps);

            auto g_y_0_y_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzzzz_x, g_y_0_y_0_zzzzz_xz, g_y_0_y_0_zzzzz_y, g_y_0_y_0_zzzzz_yz, g_y_0_y_0_zzzzz_z, g_y_0_y_0_zzzzz_zz, g_y_0_y_0_zzzzzz_x, g_y_0_y_0_zzzzzz_y, g_y_0_y_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzzzz_x[k] = -g_y_0_y_0_zzzzz_x[k] * ab_z + g_y_0_y_0_zzzzz_xz[k];

                g_y_0_y_0_zzzzzz_y[k] = -g_y_0_y_0_zzzzz_y[k] * ab_z + g_y_0_y_0_zzzzz_yz[k];

                g_y_0_y_0_zzzzzz_z[k] = -g_y_0_y_0_zzzzz_z[k] * ab_z + g_y_0_y_0_zzzzz_zz[k];
            }

            /// Set up 420-423 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 420 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 421 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 422 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxx_x, g_y_0_z_0_xxxxx_xx, g_y_0_z_0_xxxxx_xy, g_y_0_z_0_xxxxx_xz, g_y_0_z_0_xxxxx_y, g_y_0_z_0_xxxxx_z, g_y_0_z_0_xxxxxx_x, g_y_0_z_0_xxxxxx_y, g_y_0_z_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxx_x[k] = -g_y_0_z_0_xxxxx_x[k] * ab_x + g_y_0_z_0_xxxxx_xx[k];

                g_y_0_z_0_xxxxxx_y[k] = -g_y_0_z_0_xxxxx_y[k] * ab_x + g_y_0_z_0_xxxxx_xy[k];

                g_y_0_z_0_xxxxxx_z[k] = -g_y_0_z_0_xxxxx_z[k] * ab_x + g_y_0_z_0_xxxxx_xz[k];
            }

            /// Set up 423-426 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 423 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 424 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxy_x, g_y_0_z_0_xxxxxy_y, g_y_0_z_0_xxxxxy_z, g_y_0_z_0_xxxxy_x, g_y_0_z_0_xxxxy_xx, g_y_0_z_0_xxxxy_xy, g_y_0_z_0_xxxxy_xz, g_y_0_z_0_xxxxy_y, g_y_0_z_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxy_x[k] = -g_y_0_z_0_xxxxy_x[k] * ab_x + g_y_0_z_0_xxxxy_xx[k];

                g_y_0_z_0_xxxxxy_y[k] = -g_y_0_z_0_xxxxy_y[k] * ab_x + g_y_0_z_0_xxxxy_xy[k];

                g_y_0_z_0_xxxxxy_z[k] = -g_y_0_z_0_xxxxy_z[k] * ab_x + g_y_0_z_0_xxxxy_xz[k];
            }

            /// Set up 426-429 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 426 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 427 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 428 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxxz_x, g_y_0_z_0_xxxxxz_y, g_y_0_z_0_xxxxxz_z, g_y_0_z_0_xxxxz_x, g_y_0_z_0_xxxxz_xx, g_y_0_z_0_xxxxz_xy, g_y_0_z_0_xxxxz_xz, g_y_0_z_0_xxxxz_y, g_y_0_z_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxxz_x[k] = -g_y_0_z_0_xxxxz_x[k] * ab_x + g_y_0_z_0_xxxxz_xx[k];

                g_y_0_z_0_xxxxxz_y[k] = -g_y_0_z_0_xxxxz_y[k] * ab_x + g_y_0_z_0_xxxxz_xy[k];

                g_y_0_z_0_xxxxxz_z[k] = -g_y_0_z_0_xxxxz_z[k] * ab_x + g_y_0_z_0_xxxxz_xz[k];
            }

            /// Set up 429-432 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 429 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 430 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyy_x, g_y_0_z_0_xxxxyy_y, g_y_0_z_0_xxxxyy_z, g_y_0_z_0_xxxyy_x, g_y_0_z_0_xxxyy_xx, g_y_0_z_0_xxxyy_xy, g_y_0_z_0_xxxyy_xz, g_y_0_z_0_xxxyy_y, g_y_0_z_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyy_x[k] = -g_y_0_z_0_xxxyy_x[k] * ab_x + g_y_0_z_0_xxxyy_xx[k];

                g_y_0_z_0_xxxxyy_y[k] = -g_y_0_z_0_xxxyy_y[k] * ab_x + g_y_0_z_0_xxxyy_xy[k];

                g_y_0_z_0_xxxxyy_z[k] = -g_y_0_z_0_xxxyy_z[k] * ab_x + g_y_0_z_0_xxxyy_xz[k];
            }

            /// Set up 432-435 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 432 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 433 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxyz_x, g_y_0_z_0_xxxxyz_y, g_y_0_z_0_xxxxyz_z, g_y_0_z_0_xxxyz_x, g_y_0_z_0_xxxyz_xx, g_y_0_z_0_xxxyz_xy, g_y_0_z_0_xxxyz_xz, g_y_0_z_0_xxxyz_y, g_y_0_z_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxyz_x[k] = -g_y_0_z_0_xxxyz_x[k] * ab_x + g_y_0_z_0_xxxyz_xx[k];

                g_y_0_z_0_xxxxyz_y[k] = -g_y_0_z_0_xxxyz_y[k] * ab_x + g_y_0_z_0_xxxyz_xy[k];

                g_y_0_z_0_xxxxyz_z[k] = -g_y_0_z_0_xxxyz_z[k] * ab_x + g_y_0_z_0_xxxyz_xz[k];
            }

            /// Set up 435-438 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 435 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 436 * ccomps * dcomps);

            auto g_y_0_z_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxxzz_x, g_y_0_z_0_xxxxzz_y, g_y_0_z_0_xxxxzz_z, g_y_0_z_0_xxxzz_x, g_y_0_z_0_xxxzz_xx, g_y_0_z_0_xxxzz_xy, g_y_0_z_0_xxxzz_xz, g_y_0_z_0_xxxzz_y, g_y_0_z_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxxzz_x[k] = -g_y_0_z_0_xxxzz_x[k] * ab_x + g_y_0_z_0_xxxzz_xx[k];

                g_y_0_z_0_xxxxzz_y[k] = -g_y_0_z_0_xxxzz_y[k] * ab_x + g_y_0_z_0_xxxzz_xy[k];

                g_y_0_z_0_xxxxzz_z[k] = -g_y_0_z_0_xxxzz_z[k] * ab_x + g_y_0_z_0_xxxzz_xz[k];
            }

            /// Set up 438-441 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 438 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 439 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyy_x, g_y_0_z_0_xxxyyy_y, g_y_0_z_0_xxxyyy_z, g_y_0_z_0_xxyyy_x, g_y_0_z_0_xxyyy_xx, g_y_0_z_0_xxyyy_xy, g_y_0_z_0_xxyyy_xz, g_y_0_z_0_xxyyy_y, g_y_0_z_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyy_x[k] = -g_y_0_z_0_xxyyy_x[k] * ab_x + g_y_0_z_0_xxyyy_xx[k];

                g_y_0_z_0_xxxyyy_y[k] = -g_y_0_z_0_xxyyy_y[k] * ab_x + g_y_0_z_0_xxyyy_xy[k];

                g_y_0_z_0_xxxyyy_z[k] = -g_y_0_z_0_xxyyy_z[k] * ab_x + g_y_0_z_0_xxyyy_xz[k];
            }

            /// Set up 441-444 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 441 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 442 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyyz_x, g_y_0_z_0_xxxyyz_y, g_y_0_z_0_xxxyyz_z, g_y_0_z_0_xxyyz_x, g_y_0_z_0_xxyyz_xx, g_y_0_z_0_xxyyz_xy, g_y_0_z_0_xxyyz_xz, g_y_0_z_0_xxyyz_y, g_y_0_z_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyyz_x[k] = -g_y_0_z_0_xxyyz_x[k] * ab_x + g_y_0_z_0_xxyyz_xx[k];

                g_y_0_z_0_xxxyyz_y[k] = -g_y_0_z_0_xxyyz_y[k] * ab_x + g_y_0_z_0_xxyyz_xy[k];

                g_y_0_z_0_xxxyyz_z[k] = -g_y_0_z_0_xxyyz_z[k] * ab_x + g_y_0_z_0_xxyyz_xz[k];
            }

            /// Set up 444-447 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 444 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 445 * ccomps * dcomps);

            auto g_y_0_z_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 446 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxyzz_x, g_y_0_z_0_xxxyzz_y, g_y_0_z_0_xxxyzz_z, g_y_0_z_0_xxyzz_x, g_y_0_z_0_xxyzz_xx, g_y_0_z_0_xxyzz_xy, g_y_0_z_0_xxyzz_xz, g_y_0_z_0_xxyzz_y, g_y_0_z_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxyzz_x[k] = -g_y_0_z_0_xxyzz_x[k] * ab_x + g_y_0_z_0_xxyzz_xx[k];

                g_y_0_z_0_xxxyzz_y[k] = -g_y_0_z_0_xxyzz_y[k] * ab_x + g_y_0_z_0_xxyzz_xy[k];

                g_y_0_z_0_xxxyzz_z[k] = -g_y_0_z_0_xxyzz_z[k] * ab_x + g_y_0_z_0_xxyzz_xz[k];
            }

            /// Set up 447-450 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 447 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 448 * ccomps * dcomps);

            auto g_y_0_z_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxzzz_x, g_y_0_z_0_xxxzzz_y, g_y_0_z_0_xxxzzz_z, g_y_0_z_0_xxzzz_x, g_y_0_z_0_xxzzz_xx, g_y_0_z_0_xxzzz_xy, g_y_0_z_0_xxzzz_xz, g_y_0_z_0_xxzzz_y, g_y_0_z_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxzzz_x[k] = -g_y_0_z_0_xxzzz_x[k] * ab_x + g_y_0_z_0_xxzzz_xx[k];

                g_y_0_z_0_xxxzzz_y[k] = -g_y_0_z_0_xxzzz_y[k] * ab_x + g_y_0_z_0_xxzzz_xy[k];

                g_y_0_z_0_xxxzzz_z[k] = -g_y_0_z_0_xxzzz_z[k] * ab_x + g_y_0_z_0_xxzzz_xz[k];
            }

            /// Set up 450-453 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 450 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 451 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 452 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyy_x, g_y_0_z_0_xxyyyy_y, g_y_0_z_0_xxyyyy_z, g_y_0_z_0_xyyyy_x, g_y_0_z_0_xyyyy_xx, g_y_0_z_0_xyyyy_xy, g_y_0_z_0_xyyyy_xz, g_y_0_z_0_xyyyy_y, g_y_0_z_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyy_x[k] = -g_y_0_z_0_xyyyy_x[k] * ab_x + g_y_0_z_0_xyyyy_xx[k];

                g_y_0_z_0_xxyyyy_y[k] = -g_y_0_z_0_xyyyy_y[k] * ab_x + g_y_0_z_0_xyyyy_xy[k];

                g_y_0_z_0_xxyyyy_z[k] = -g_y_0_z_0_xyyyy_z[k] * ab_x + g_y_0_z_0_xyyyy_xz[k];
            }

            /// Set up 453-456 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 453 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 454 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyyz_x, g_y_0_z_0_xxyyyz_y, g_y_0_z_0_xxyyyz_z, g_y_0_z_0_xyyyz_x, g_y_0_z_0_xyyyz_xx, g_y_0_z_0_xyyyz_xy, g_y_0_z_0_xyyyz_xz, g_y_0_z_0_xyyyz_y, g_y_0_z_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyyz_x[k] = -g_y_0_z_0_xyyyz_x[k] * ab_x + g_y_0_z_0_xyyyz_xx[k];

                g_y_0_z_0_xxyyyz_y[k] = -g_y_0_z_0_xyyyz_y[k] * ab_x + g_y_0_z_0_xyyyz_xy[k];

                g_y_0_z_0_xxyyyz_z[k] = -g_y_0_z_0_xyyyz_z[k] * ab_x + g_y_0_z_0_xyyyz_xz[k];
            }

            /// Set up 456-459 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 456 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 457 * ccomps * dcomps);

            auto g_y_0_z_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 458 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyyzz_x, g_y_0_z_0_xxyyzz_y, g_y_0_z_0_xxyyzz_z, g_y_0_z_0_xyyzz_x, g_y_0_z_0_xyyzz_xx, g_y_0_z_0_xyyzz_xy, g_y_0_z_0_xyyzz_xz, g_y_0_z_0_xyyzz_y, g_y_0_z_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyyzz_x[k] = -g_y_0_z_0_xyyzz_x[k] * ab_x + g_y_0_z_0_xyyzz_xx[k];

                g_y_0_z_0_xxyyzz_y[k] = -g_y_0_z_0_xyyzz_y[k] * ab_x + g_y_0_z_0_xyyzz_xy[k];

                g_y_0_z_0_xxyyzz_z[k] = -g_y_0_z_0_xyyzz_z[k] * ab_x + g_y_0_z_0_xyyzz_xz[k];
            }

            /// Set up 459-462 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 459 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 460 * ccomps * dcomps);

            auto g_y_0_z_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyzzz_x, g_y_0_z_0_xxyzzz_y, g_y_0_z_0_xxyzzz_z, g_y_0_z_0_xyzzz_x, g_y_0_z_0_xyzzz_xx, g_y_0_z_0_xyzzz_xy, g_y_0_z_0_xyzzz_xz, g_y_0_z_0_xyzzz_y, g_y_0_z_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyzzz_x[k] = -g_y_0_z_0_xyzzz_x[k] * ab_x + g_y_0_z_0_xyzzz_xx[k];

                g_y_0_z_0_xxyzzz_y[k] = -g_y_0_z_0_xyzzz_y[k] * ab_x + g_y_0_z_0_xyzzz_xy[k];

                g_y_0_z_0_xxyzzz_z[k] = -g_y_0_z_0_xyzzz_z[k] * ab_x + g_y_0_z_0_xyzzz_xz[k];
            }

            /// Set up 462-465 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 462 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 463 * ccomps * dcomps);

            auto g_y_0_z_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzzzz_x, g_y_0_z_0_xxzzzz_y, g_y_0_z_0_xxzzzz_z, g_y_0_z_0_xzzzz_x, g_y_0_z_0_xzzzz_xx, g_y_0_z_0_xzzzz_xy, g_y_0_z_0_xzzzz_xz, g_y_0_z_0_xzzzz_y, g_y_0_z_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzzzz_x[k] = -g_y_0_z_0_xzzzz_x[k] * ab_x + g_y_0_z_0_xzzzz_xx[k];

                g_y_0_z_0_xxzzzz_y[k] = -g_y_0_z_0_xzzzz_y[k] * ab_x + g_y_0_z_0_xzzzz_xy[k];

                g_y_0_z_0_xxzzzz_z[k] = -g_y_0_z_0_xzzzz_z[k] * ab_x + g_y_0_z_0_xzzzz_xz[k];
            }

            /// Set up 465-468 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 465 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 466 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyy_x, g_y_0_z_0_xyyyyy_y, g_y_0_z_0_xyyyyy_z, g_y_0_z_0_yyyyy_x, g_y_0_z_0_yyyyy_xx, g_y_0_z_0_yyyyy_xy, g_y_0_z_0_yyyyy_xz, g_y_0_z_0_yyyyy_y, g_y_0_z_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyy_x[k] = -g_y_0_z_0_yyyyy_x[k] * ab_x + g_y_0_z_0_yyyyy_xx[k];

                g_y_0_z_0_xyyyyy_y[k] = -g_y_0_z_0_yyyyy_y[k] * ab_x + g_y_0_z_0_yyyyy_xy[k];

                g_y_0_z_0_xyyyyy_z[k] = -g_y_0_z_0_yyyyy_z[k] * ab_x + g_y_0_z_0_yyyyy_xz[k];
            }

            /// Set up 468-471 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 468 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 469 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 470 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyyz_x, g_y_0_z_0_xyyyyz_y, g_y_0_z_0_xyyyyz_z, g_y_0_z_0_yyyyz_x, g_y_0_z_0_yyyyz_xx, g_y_0_z_0_yyyyz_xy, g_y_0_z_0_yyyyz_xz, g_y_0_z_0_yyyyz_y, g_y_0_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyyz_x[k] = -g_y_0_z_0_yyyyz_x[k] * ab_x + g_y_0_z_0_yyyyz_xx[k];

                g_y_0_z_0_xyyyyz_y[k] = -g_y_0_z_0_yyyyz_y[k] * ab_x + g_y_0_z_0_yyyyz_xy[k];

                g_y_0_z_0_xyyyyz_z[k] = -g_y_0_z_0_yyyyz_z[k] * ab_x + g_y_0_z_0_yyyyz_xz[k];
            }

            /// Set up 471-474 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 471 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 472 * ccomps * dcomps);

            auto g_y_0_z_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyyzz_x, g_y_0_z_0_xyyyzz_y, g_y_0_z_0_xyyyzz_z, g_y_0_z_0_yyyzz_x, g_y_0_z_0_yyyzz_xx, g_y_0_z_0_yyyzz_xy, g_y_0_z_0_yyyzz_xz, g_y_0_z_0_yyyzz_y, g_y_0_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyyzz_x[k] = -g_y_0_z_0_yyyzz_x[k] * ab_x + g_y_0_z_0_yyyzz_xx[k];

                g_y_0_z_0_xyyyzz_y[k] = -g_y_0_z_0_yyyzz_y[k] * ab_x + g_y_0_z_0_yyyzz_xy[k];

                g_y_0_z_0_xyyyzz_z[k] = -g_y_0_z_0_yyyzz_z[k] * ab_x + g_y_0_z_0_yyyzz_xz[k];
            }

            /// Set up 474-477 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 474 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 475 * ccomps * dcomps);

            auto g_y_0_z_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 476 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyzzz_x, g_y_0_z_0_xyyzzz_y, g_y_0_z_0_xyyzzz_z, g_y_0_z_0_yyzzz_x, g_y_0_z_0_yyzzz_xx, g_y_0_z_0_yyzzz_xy, g_y_0_z_0_yyzzz_xz, g_y_0_z_0_yyzzz_y, g_y_0_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyzzz_x[k] = -g_y_0_z_0_yyzzz_x[k] * ab_x + g_y_0_z_0_yyzzz_xx[k];

                g_y_0_z_0_xyyzzz_y[k] = -g_y_0_z_0_yyzzz_y[k] * ab_x + g_y_0_z_0_yyzzz_xy[k];

                g_y_0_z_0_xyyzzz_z[k] = -g_y_0_z_0_yyzzz_z[k] * ab_x + g_y_0_z_0_yyzzz_xz[k];
            }

            /// Set up 477-480 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 477 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 478 * ccomps * dcomps);

            auto g_y_0_z_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzzzz_x, g_y_0_z_0_xyzzzz_y, g_y_0_z_0_xyzzzz_z, g_y_0_z_0_yzzzz_x, g_y_0_z_0_yzzzz_xx, g_y_0_z_0_yzzzz_xy, g_y_0_z_0_yzzzz_xz, g_y_0_z_0_yzzzz_y, g_y_0_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzzzz_x[k] = -g_y_0_z_0_yzzzz_x[k] * ab_x + g_y_0_z_0_yzzzz_xx[k];

                g_y_0_z_0_xyzzzz_y[k] = -g_y_0_z_0_yzzzz_y[k] * ab_x + g_y_0_z_0_yzzzz_xy[k];

                g_y_0_z_0_xyzzzz_z[k] = -g_y_0_z_0_yzzzz_z[k] * ab_x + g_y_0_z_0_yzzzz_xz[k];
            }

            /// Set up 480-483 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 480 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 481 * ccomps * dcomps);

            auto g_y_0_z_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzzzz_x, g_y_0_z_0_xzzzzz_y, g_y_0_z_0_xzzzzz_z, g_y_0_z_0_zzzzz_x, g_y_0_z_0_zzzzz_xx, g_y_0_z_0_zzzzz_xy, g_y_0_z_0_zzzzz_xz, g_y_0_z_0_zzzzz_y, g_y_0_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzzzz_x[k] = -g_y_0_z_0_zzzzz_x[k] * ab_x + g_y_0_z_0_zzzzz_xx[k];

                g_y_0_z_0_xzzzzz_y[k] = -g_y_0_z_0_zzzzz_y[k] * ab_x + g_y_0_z_0_zzzzz_xy[k];

                g_y_0_z_0_xzzzzz_z[k] = -g_y_0_z_0_zzzzz_z[k] * ab_x + g_y_0_z_0_zzzzz_xz[k];
            }

            /// Set up 483-486 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 483 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 484 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyyyy_x, g_0_0_z_0_yyyyy_y, g_0_0_z_0_yyyyy_z, g_y_0_z_0_yyyyy_x, g_y_0_z_0_yyyyy_xy, g_y_0_z_0_yyyyy_y, g_y_0_z_0_yyyyy_yy, g_y_0_z_0_yyyyy_yz, g_y_0_z_0_yyyyy_z, g_y_0_z_0_yyyyyy_x, g_y_0_z_0_yyyyyy_y, g_y_0_z_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyy_x[k] = -g_0_0_z_0_yyyyy_x[k] - g_y_0_z_0_yyyyy_x[k] * ab_y + g_y_0_z_0_yyyyy_xy[k];

                g_y_0_z_0_yyyyyy_y[k] = -g_0_0_z_0_yyyyy_y[k] - g_y_0_z_0_yyyyy_y[k] * ab_y + g_y_0_z_0_yyyyy_yy[k];

                g_y_0_z_0_yyyyyy_z[k] = -g_0_0_z_0_yyyyy_z[k] - g_y_0_z_0_yyyyy_z[k] * ab_y + g_y_0_z_0_yyyyy_yz[k];
            }

            /// Set up 486-489 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 486 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 487 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 488 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyy_x, g_y_0_z_0_yyyyy_xz, g_y_0_z_0_yyyyy_y, g_y_0_z_0_yyyyy_yz, g_y_0_z_0_yyyyy_z, g_y_0_z_0_yyyyy_zz, g_y_0_z_0_yyyyyz_x, g_y_0_z_0_yyyyyz_y, g_y_0_z_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyyz_x[k] = -g_y_0_z_0_yyyyy_x[k] * ab_z + g_y_0_z_0_yyyyy_xz[k];

                g_y_0_z_0_yyyyyz_y[k] = -g_y_0_z_0_yyyyy_y[k] * ab_z + g_y_0_z_0_yyyyy_yz[k];

                g_y_0_z_0_yyyyyz_z[k] = -g_y_0_z_0_yyyyy_z[k] * ab_z + g_y_0_z_0_yyyyy_zz[k];
            }

            /// Set up 489-492 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 489 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 490 * ccomps * dcomps);

            auto g_y_0_z_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyyz_x, g_y_0_z_0_yyyyz_xz, g_y_0_z_0_yyyyz_y, g_y_0_z_0_yyyyz_yz, g_y_0_z_0_yyyyz_z, g_y_0_z_0_yyyyz_zz, g_y_0_z_0_yyyyzz_x, g_y_0_z_0_yyyyzz_y, g_y_0_z_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyyzz_x[k] = -g_y_0_z_0_yyyyz_x[k] * ab_z + g_y_0_z_0_yyyyz_xz[k];

                g_y_0_z_0_yyyyzz_y[k] = -g_y_0_z_0_yyyyz_y[k] * ab_z + g_y_0_z_0_yyyyz_yz[k];

                g_y_0_z_0_yyyyzz_z[k] = -g_y_0_z_0_yyyyz_z[k] * ab_z + g_y_0_z_0_yyyyz_zz[k];
            }

            /// Set up 492-495 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 492 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 493 * ccomps * dcomps);

            auto g_y_0_z_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyyzz_x, g_y_0_z_0_yyyzz_xz, g_y_0_z_0_yyyzz_y, g_y_0_z_0_yyyzz_yz, g_y_0_z_0_yyyzz_z, g_y_0_z_0_yyyzz_zz, g_y_0_z_0_yyyzzz_x, g_y_0_z_0_yyyzzz_y, g_y_0_z_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyzzz_x[k] = -g_y_0_z_0_yyyzz_x[k] * ab_z + g_y_0_z_0_yyyzz_xz[k];

                g_y_0_z_0_yyyzzz_y[k] = -g_y_0_z_0_yyyzz_y[k] * ab_z + g_y_0_z_0_yyyzz_yz[k];

                g_y_0_z_0_yyyzzz_z[k] = -g_y_0_z_0_yyyzz_z[k] * ab_z + g_y_0_z_0_yyyzz_zz[k];
            }

            /// Set up 495-498 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 495 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 496 * ccomps * dcomps);

            auto g_y_0_z_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyzzz_x, g_y_0_z_0_yyzzz_xz, g_y_0_z_0_yyzzz_y, g_y_0_z_0_yyzzz_yz, g_y_0_z_0_yyzzz_z, g_y_0_z_0_yyzzz_zz, g_y_0_z_0_yyzzzz_x, g_y_0_z_0_yyzzzz_y, g_y_0_z_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzzzz_x[k] = -g_y_0_z_0_yyzzz_x[k] * ab_z + g_y_0_z_0_yyzzz_xz[k];

                g_y_0_z_0_yyzzzz_y[k] = -g_y_0_z_0_yyzzz_y[k] * ab_z + g_y_0_z_0_yyzzz_yz[k];

                g_y_0_z_0_yyzzzz_z[k] = -g_y_0_z_0_yyzzz_z[k] * ab_z + g_y_0_z_0_yyzzz_zz[k];
            }

            /// Set up 498-501 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 498 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 499 * ccomps * dcomps);

            auto g_y_0_z_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 500 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzzzz_x, g_y_0_z_0_yzzzz_xz, g_y_0_z_0_yzzzz_y, g_y_0_z_0_yzzzz_yz, g_y_0_z_0_yzzzz_z, g_y_0_z_0_yzzzz_zz, g_y_0_z_0_yzzzzz_x, g_y_0_z_0_yzzzzz_y, g_y_0_z_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzzzz_x[k] = -g_y_0_z_0_yzzzz_x[k] * ab_z + g_y_0_z_0_yzzzz_xz[k];

                g_y_0_z_0_yzzzzz_y[k] = -g_y_0_z_0_yzzzz_y[k] * ab_z + g_y_0_z_0_yzzzz_yz[k];

                g_y_0_z_0_yzzzzz_z[k] = -g_y_0_z_0_yzzzz_z[k] * ab_z + g_y_0_z_0_yzzzz_zz[k];
            }

            /// Set up 501-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 501 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 502 * ccomps * dcomps);

            auto g_y_0_z_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzzzz_x, g_y_0_z_0_zzzzz_xz, g_y_0_z_0_zzzzz_y, g_y_0_z_0_zzzzz_yz, g_y_0_z_0_zzzzz_z, g_y_0_z_0_zzzzz_zz, g_y_0_z_0_zzzzzz_x, g_y_0_z_0_zzzzzz_y, g_y_0_z_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzzzz_x[k] = -g_y_0_z_0_zzzzz_x[k] * ab_z + g_y_0_z_0_zzzzz_xz[k];

                g_y_0_z_0_zzzzzz_y[k] = -g_y_0_z_0_zzzzz_y[k] * ab_z + g_y_0_z_0_zzzzz_yz[k];

                g_y_0_z_0_zzzzzz_z[k] = -g_y_0_z_0_zzzzz_z[k] * ab_z + g_y_0_z_0_zzzzz_zz[k];
            }

            /// Set up 504-507 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 504 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 505 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 506 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxx_x, g_z_0_x_0_xxxxx_xx, g_z_0_x_0_xxxxx_xy, g_z_0_x_0_xxxxx_xz, g_z_0_x_0_xxxxx_y, g_z_0_x_0_xxxxx_z, g_z_0_x_0_xxxxxx_x, g_z_0_x_0_xxxxxx_y, g_z_0_x_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxx_x[k] = -g_z_0_x_0_xxxxx_x[k] * ab_x + g_z_0_x_0_xxxxx_xx[k];

                g_z_0_x_0_xxxxxx_y[k] = -g_z_0_x_0_xxxxx_y[k] * ab_x + g_z_0_x_0_xxxxx_xy[k];

                g_z_0_x_0_xxxxxx_z[k] = -g_z_0_x_0_xxxxx_z[k] * ab_x + g_z_0_x_0_xxxxx_xz[k];
            }

            /// Set up 507-510 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 507 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 508 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxy_x, g_z_0_x_0_xxxxxy_y, g_z_0_x_0_xxxxxy_z, g_z_0_x_0_xxxxy_x, g_z_0_x_0_xxxxy_xx, g_z_0_x_0_xxxxy_xy, g_z_0_x_0_xxxxy_xz, g_z_0_x_0_xxxxy_y, g_z_0_x_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxy_x[k] = -g_z_0_x_0_xxxxy_x[k] * ab_x + g_z_0_x_0_xxxxy_xx[k];

                g_z_0_x_0_xxxxxy_y[k] = -g_z_0_x_0_xxxxy_y[k] * ab_x + g_z_0_x_0_xxxxy_xy[k];

                g_z_0_x_0_xxxxxy_z[k] = -g_z_0_x_0_xxxxy_z[k] * ab_x + g_z_0_x_0_xxxxy_xz[k];
            }

            /// Set up 510-513 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 510 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 511 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 512 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxxz_x, g_z_0_x_0_xxxxxz_y, g_z_0_x_0_xxxxxz_z, g_z_0_x_0_xxxxz_x, g_z_0_x_0_xxxxz_xx, g_z_0_x_0_xxxxz_xy, g_z_0_x_0_xxxxz_xz, g_z_0_x_0_xxxxz_y, g_z_0_x_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxxz_x[k] = -g_z_0_x_0_xxxxz_x[k] * ab_x + g_z_0_x_0_xxxxz_xx[k];

                g_z_0_x_0_xxxxxz_y[k] = -g_z_0_x_0_xxxxz_y[k] * ab_x + g_z_0_x_0_xxxxz_xy[k];

                g_z_0_x_0_xxxxxz_z[k] = -g_z_0_x_0_xxxxz_z[k] * ab_x + g_z_0_x_0_xxxxz_xz[k];
            }

            /// Set up 513-516 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 513 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 514 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyy_x, g_z_0_x_0_xxxxyy_y, g_z_0_x_0_xxxxyy_z, g_z_0_x_0_xxxyy_x, g_z_0_x_0_xxxyy_xx, g_z_0_x_0_xxxyy_xy, g_z_0_x_0_xxxyy_xz, g_z_0_x_0_xxxyy_y, g_z_0_x_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyy_x[k] = -g_z_0_x_0_xxxyy_x[k] * ab_x + g_z_0_x_0_xxxyy_xx[k];

                g_z_0_x_0_xxxxyy_y[k] = -g_z_0_x_0_xxxyy_y[k] * ab_x + g_z_0_x_0_xxxyy_xy[k];

                g_z_0_x_0_xxxxyy_z[k] = -g_z_0_x_0_xxxyy_z[k] * ab_x + g_z_0_x_0_xxxyy_xz[k];
            }

            /// Set up 516-519 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 516 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 517 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 518 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxyz_x, g_z_0_x_0_xxxxyz_y, g_z_0_x_0_xxxxyz_z, g_z_0_x_0_xxxyz_x, g_z_0_x_0_xxxyz_xx, g_z_0_x_0_xxxyz_xy, g_z_0_x_0_xxxyz_xz, g_z_0_x_0_xxxyz_y, g_z_0_x_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxyz_x[k] = -g_z_0_x_0_xxxyz_x[k] * ab_x + g_z_0_x_0_xxxyz_xx[k];

                g_z_0_x_0_xxxxyz_y[k] = -g_z_0_x_0_xxxyz_y[k] * ab_x + g_z_0_x_0_xxxyz_xy[k];

                g_z_0_x_0_xxxxyz_z[k] = -g_z_0_x_0_xxxyz_z[k] * ab_x + g_z_0_x_0_xxxyz_xz[k];
            }

            /// Set up 519-522 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 519 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 520 * ccomps * dcomps);

            auto g_z_0_x_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxxzz_x, g_z_0_x_0_xxxxzz_y, g_z_0_x_0_xxxxzz_z, g_z_0_x_0_xxxzz_x, g_z_0_x_0_xxxzz_xx, g_z_0_x_0_xxxzz_xy, g_z_0_x_0_xxxzz_xz, g_z_0_x_0_xxxzz_y, g_z_0_x_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxxzz_x[k] = -g_z_0_x_0_xxxzz_x[k] * ab_x + g_z_0_x_0_xxxzz_xx[k];

                g_z_0_x_0_xxxxzz_y[k] = -g_z_0_x_0_xxxzz_y[k] * ab_x + g_z_0_x_0_xxxzz_xy[k];

                g_z_0_x_0_xxxxzz_z[k] = -g_z_0_x_0_xxxzz_z[k] * ab_x + g_z_0_x_0_xxxzz_xz[k];
            }

            /// Set up 522-525 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 522 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 523 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyy_x, g_z_0_x_0_xxxyyy_y, g_z_0_x_0_xxxyyy_z, g_z_0_x_0_xxyyy_x, g_z_0_x_0_xxyyy_xx, g_z_0_x_0_xxyyy_xy, g_z_0_x_0_xxyyy_xz, g_z_0_x_0_xxyyy_y, g_z_0_x_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyy_x[k] = -g_z_0_x_0_xxyyy_x[k] * ab_x + g_z_0_x_0_xxyyy_xx[k];

                g_z_0_x_0_xxxyyy_y[k] = -g_z_0_x_0_xxyyy_y[k] * ab_x + g_z_0_x_0_xxyyy_xy[k];

                g_z_0_x_0_xxxyyy_z[k] = -g_z_0_x_0_xxyyy_z[k] * ab_x + g_z_0_x_0_xxyyy_xz[k];
            }

            /// Set up 525-528 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 525 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 526 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyyz_x, g_z_0_x_0_xxxyyz_y, g_z_0_x_0_xxxyyz_z, g_z_0_x_0_xxyyz_x, g_z_0_x_0_xxyyz_xx, g_z_0_x_0_xxyyz_xy, g_z_0_x_0_xxyyz_xz, g_z_0_x_0_xxyyz_y, g_z_0_x_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyyz_x[k] = -g_z_0_x_0_xxyyz_x[k] * ab_x + g_z_0_x_0_xxyyz_xx[k];

                g_z_0_x_0_xxxyyz_y[k] = -g_z_0_x_0_xxyyz_y[k] * ab_x + g_z_0_x_0_xxyyz_xy[k];

                g_z_0_x_0_xxxyyz_z[k] = -g_z_0_x_0_xxyyz_z[k] * ab_x + g_z_0_x_0_xxyyz_xz[k];
            }

            /// Set up 528-531 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 528 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 529 * ccomps * dcomps);

            auto g_z_0_x_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 530 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxyzz_x, g_z_0_x_0_xxxyzz_y, g_z_0_x_0_xxxyzz_z, g_z_0_x_0_xxyzz_x, g_z_0_x_0_xxyzz_xx, g_z_0_x_0_xxyzz_xy, g_z_0_x_0_xxyzz_xz, g_z_0_x_0_xxyzz_y, g_z_0_x_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxyzz_x[k] = -g_z_0_x_0_xxyzz_x[k] * ab_x + g_z_0_x_0_xxyzz_xx[k];

                g_z_0_x_0_xxxyzz_y[k] = -g_z_0_x_0_xxyzz_y[k] * ab_x + g_z_0_x_0_xxyzz_xy[k];

                g_z_0_x_0_xxxyzz_z[k] = -g_z_0_x_0_xxyzz_z[k] * ab_x + g_z_0_x_0_xxyzz_xz[k];
            }

            /// Set up 531-534 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 531 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 532 * ccomps * dcomps);

            auto g_z_0_x_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxzzz_x, g_z_0_x_0_xxxzzz_y, g_z_0_x_0_xxxzzz_z, g_z_0_x_0_xxzzz_x, g_z_0_x_0_xxzzz_xx, g_z_0_x_0_xxzzz_xy, g_z_0_x_0_xxzzz_xz, g_z_0_x_0_xxzzz_y, g_z_0_x_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxzzz_x[k] = -g_z_0_x_0_xxzzz_x[k] * ab_x + g_z_0_x_0_xxzzz_xx[k];

                g_z_0_x_0_xxxzzz_y[k] = -g_z_0_x_0_xxzzz_y[k] * ab_x + g_z_0_x_0_xxzzz_xy[k];

                g_z_0_x_0_xxxzzz_z[k] = -g_z_0_x_0_xxzzz_z[k] * ab_x + g_z_0_x_0_xxzzz_xz[k];
            }

            /// Set up 534-537 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 534 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 535 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 536 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyy_x, g_z_0_x_0_xxyyyy_y, g_z_0_x_0_xxyyyy_z, g_z_0_x_0_xyyyy_x, g_z_0_x_0_xyyyy_xx, g_z_0_x_0_xyyyy_xy, g_z_0_x_0_xyyyy_xz, g_z_0_x_0_xyyyy_y, g_z_0_x_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyy_x[k] = -g_z_0_x_0_xyyyy_x[k] * ab_x + g_z_0_x_0_xyyyy_xx[k];

                g_z_0_x_0_xxyyyy_y[k] = -g_z_0_x_0_xyyyy_y[k] * ab_x + g_z_0_x_0_xyyyy_xy[k];

                g_z_0_x_0_xxyyyy_z[k] = -g_z_0_x_0_xyyyy_z[k] * ab_x + g_z_0_x_0_xyyyy_xz[k];
            }

            /// Set up 537-540 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 537 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 538 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyyz_x, g_z_0_x_0_xxyyyz_y, g_z_0_x_0_xxyyyz_z, g_z_0_x_0_xyyyz_x, g_z_0_x_0_xyyyz_xx, g_z_0_x_0_xyyyz_xy, g_z_0_x_0_xyyyz_xz, g_z_0_x_0_xyyyz_y, g_z_0_x_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyyz_x[k] = -g_z_0_x_0_xyyyz_x[k] * ab_x + g_z_0_x_0_xyyyz_xx[k];

                g_z_0_x_0_xxyyyz_y[k] = -g_z_0_x_0_xyyyz_y[k] * ab_x + g_z_0_x_0_xyyyz_xy[k];

                g_z_0_x_0_xxyyyz_z[k] = -g_z_0_x_0_xyyyz_z[k] * ab_x + g_z_0_x_0_xyyyz_xz[k];
            }

            /// Set up 540-543 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 540 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 541 * ccomps * dcomps);

            auto g_z_0_x_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 542 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyyzz_x, g_z_0_x_0_xxyyzz_y, g_z_0_x_0_xxyyzz_z, g_z_0_x_0_xyyzz_x, g_z_0_x_0_xyyzz_xx, g_z_0_x_0_xyyzz_xy, g_z_0_x_0_xyyzz_xz, g_z_0_x_0_xyyzz_y, g_z_0_x_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyyzz_x[k] = -g_z_0_x_0_xyyzz_x[k] * ab_x + g_z_0_x_0_xyyzz_xx[k];

                g_z_0_x_0_xxyyzz_y[k] = -g_z_0_x_0_xyyzz_y[k] * ab_x + g_z_0_x_0_xyyzz_xy[k];

                g_z_0_x_0_xxyyzz_z[k] = -g_z_0_x_0_xyyzz_z[k] * ab_x + g_z_0_x_0_xyyzz_xz[k];
            }

            /// Set up 543-546 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 543 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 544 * ccomps * dcomps);

            auto g_z_0_x_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyzzz_x, g_z_0_x_0_xxyzzz_y, g_z_0_x_0_xxyzzz_z, g_z_0_x_0_xyzzz_x, g_z_0_x_0_xyzzz_xx, g_z_0_x_0_xyzzz_xy, g_z_0_x_0_xyzzz_xz, g_z_0_x_0_xyzzz_y, g_z_0_x_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyzzz_x[k] = -g_z_0_x_0_xyzzz_x[k] * ab_x + g_z_0_x_0_xyzzz_xx[k];

                g_z_0_x_0_xxyzzz_y[k] = -g_z_0_x_0_xyzzz_y[k] * ab_x + g_z_0_x_0_xyzzz_xy[k];

                g_z_0_x_0_xxyzzz_z[k] = -g_z_0_x_0_xyzzz_z[k] * ab_x + g_z_0_x_0_xyzzz_xz[k];
            }

            /// Set up 546-549 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 546 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 547 * ccomps * dcomps);

            auto g_z_0_x_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 548 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzzzz_x, g_z_0_x_0_xxzzzz_y, g_z_0_x_0_xxzzzz_z, g_z_0_x_0_xzzzz_x, g_z_0_x_0_xzzzz_xx, g_z_0_x_0_xzzzz_xy, g_z_0_x_0_xzzzz_xz, g_z_0_x_0_xzzzz_y, g_z_0_x_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzzzz_x[k] = -g_z_0_x_0_xzzzz_x[k] * ab_x + g_z_0_x_0_xzzzz_xx[k];

                g_z_0_x_0_xxzzzz_y[k] = -g_z_0_x_0_xzzzz_y[k] * ab_x + g_z_0_x_0_xzzzz_xy[k];

                g_z_0_x_0_xxzzzz_z[k] = -g_z_0_x_0_xzzzz_z[k] * ab_x + g_z_0_x_0_xzzzz_xz[k];
            }

            /// Set up 549-552 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 549 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 550 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyy_x, g_z_0_x_0_xyyyyy_y, g_z_0_x_0_xyyyyy_z, g_z_0_x_0_yyyyy_x, g_z_0_x_0_yyyyy_xx, g_z_0_x_0_yyyyy_xy, g_z_0_x_0_yyyyy_xz, g_z_0_x_0_yyyyy_y, g_z_0_x_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyy_x[k] = -g_z_0_x_0_yyyyy_x[k] * ab_x + g_z_0_x_0_yyyyy_xx[k];

                g_z_0_x_0_xyyyyy_y[k] = -g_z_0_x_0_yyyyy_y[k] * ab_x + g_z_0_x_0_yyyyy_xy[k];

                g_z_0_x_0_xyyyyy_z[k] = -g_z_0_x_0_yyyyy_z[k] * ab_x + g_z_0_x_0_yyyyy_xz[k];
            }

            /// Set up 552-555 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 552 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 553 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyyz_x, g_z_0_x_0_xyyyyz_y, g_z_0_x_0_xyyyyz_z, g_z_0_x_0_yyyyz_x, g_z_0_x_0_yyyyz_xx, g_z_0_x_0_yyyyz_xy, g_z_0_x_0_yyyyz_xz, g_z_0_x_0_yyyyz_y, g_z_0_x_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyyz_x[k] = -g_z_0_x_0_yyyyz_x[k] * ab_x + g_z_0_x_0_yyyyz_xx[k];

                g_z_0_x_0_xyyyyz_y[k] = -g_z_0_x_0_yyyyz_y[k] * ab_x + g_z_0_x_0_yyyyz_xy[k];

                g_z_0_x_0_xyyyyz_z[k] = -g_z_0_x_0_yyyyz_z[k] * ab_x + g_z_0_x_0_yyyyz_xz[k];
            }

            /// Set up 555-558 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 555 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 556 * ccomps * dcomps);

            auto g_z_0_x_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyyzz_x, g_z_0_x_0_xyyyzz_y, g_z_0_x_0_xyyyzz_z, g_z_0_x_0_yyyzz_x, g_z_0_x_0_yyyzz_xx, g_z_0_x_0_yyyzz_xy, g_z_0_x_0_yyyzz_xz, g_z_0_x_0_yyyzz_y, g_z_0_x_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyyzz_x[k] = -g_z_0_x_0_yyyzz_x[k] * ab_x + g_z_0_x_0_yyyzz_xx[k];

                g_z_0_x_0_xyyyzz_y[k] = -g_z_0_x_0_yyyzz_y[k] * ab_x + g_z_0_x_0_yyyzz_xy[k];

                g_z_0_x_0_xyyyzz_z[k] = -g_z_0_x_0_yyyzz_z[k] * ab_x + g_z_0_x_0_yyyzz_xz[k];
            }

            /// Set up 558-561 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 558 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 559 * ccomps * dcomps);

            auto g_z_0_x_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 560 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyzzz_x, g_z_0_x_0_xyyzzz_y, g_z_0_x_0_xyyzzz_z, g_z_0_x_0_yyzzz_x, g_z_0_x_0_yyzzz_xx, g_z_0_x_0_yyzzz_xy, g_z_0_x_0_yyzzz_xz, g_z_0_x_0_yyzzz_y, g_z_0_x_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyzzz_x[k] = -g_z_0_x_0_yyzzz_x[k] * ab_x + g_z_0_x_0_yyzzz_xx[k];

                g_z_0_x_0_xyyzzz_y[k] = -g_z_0_x_0_yyzzz_y[k] * ab_x + g_z_0_x_0_yyzzz_xy[k];

                g_z_0_x_0_xyyzzz_z[k] = -g_z_0_x_0_yyzzz_z[k] * ab_x + g_z_0_x_0_yyzzz_xz[k];
            }

            /// Set up 561-564 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 561 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 562 * ccomps * dcomps);

            auto g_z_0_x_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzzzz_x, g_z_0_x_0_xyzzzz_y, g_z_0_x_0_xyzzzz_z, g_z_0_x_0_yzzzz_x, g_z_0_x_0_yzzzz_xx, g_z_0_x_0_yzzzz_xy, g_z_0_x_0_yzzzz_xz, g_z_0_x_0_yzzzz_y, g_z_0_x_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzzzz_x[k] = -g_z_0_x_0_yzzzz_x[k] * ab_x + g_z_0_x_0_yzzzz_xx[k];

                g_z_0_x_0_xyzzzz_y[k] = -g_z_0_x_0_yzzzz_y[k] * ab_x + g_z_0_x_0_yzzzz_xy[k];

                g_z_0_x_0_xyzzzz_z[k] = -g_z_0_x_0_yzzzz_z[k] * ab_x + g_z_0_x_0_yzzzz_xz[k];
            }

            /// Set up 564-567 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 564 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 565 * ccomps * dcomps);

            auto g_z_0_x_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzzzz_x, g_z_0_x_0_xzzzzz_y, g_z_0_x_0_xzzzzz_z, g_z_0_x_0_zzzzz_x, g_z_0_x_0_zzzzz_xx, g_z_0_x_0_zzzzz_xy, g_z_0_x_0_zzzzz_xz, g_z_0_x_0_zzzzz_y, g_z_0_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzzzz_x[k] = -g_z_0_x_0_zzzzz_x[k] * ab_x + g_z_0_x_0_zzzzz_xx[k];

                g_z_0_x_0_xzzzzz_y[k] = -g_z_0_x_0_zzzzz_y[k] * ab_x + g_z_0_x_0_zzzzz_xy[k];

                g_z_0_x_0_xzzzzz_z[k] = -g_z_0_x_0_zzzzz_z[k] * ab_x + g_z_0_x_0_zzzzz_xz[k];
            }

            /// Set up 567-570 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 567 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 568 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyy_x, g_z_0_x_0_yyyyy_xy, g_z_0_x_0_yyyyy_y, g_z_0_x_0_yyyyy_yy, g_z_0_x_0_yyyyy_yz, g_z_0_x_0_yyyyy_z, g_z_0_x_0_yyyyyy_x, g_z_0_x_0_yyyyyy_y, g_z_0_x_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyy_x[k] = -g_z_0_x_0_yyyyy_x[k] * ab_y + g_z_0_x_0_yyyyy_xy[k];

                g_z_0_x_0_yyyyyy_y[k] = -g_z_0_x_0_yyyyy_y[k] * ab_y + g_z_0_x_0_yyyyy_yy[k];

                g_z_0_x_0_yyyyyy_z[k] = -g_z_0_x_0_yyyyy_z[k] * ab_y + g_z_0_x_0_yyyyy_yz[k];
            }

            /// Set up 570-573 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 570 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 571 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 572 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyyz_x, g_z_0_x_0_yyyyyz_y, g_z_0_x_0_yyyyyz_z, g_z_0_x_0_yyyyz_x, g_z_0_x_0_yyyyz_xy, g_z_0_x_0_yyyyz_y, g_z_0_x_0_yyyyz_yy, g_z_0_x_0_yyyyz_yz, g_z_0_x_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyyz_x[k] = -g_z_0_x_0_yyyyz_x[k] * ab_y + g_z_0_x_0_yyyyz_xy[k];

                g_z_0_x_0_yyyyyz_y[k] = -g_z_0_x_0_yyyyz_y[k] * ab_y + g_z_0_x_0_yyyyz_yy[k];

                g_z_0_x_0_yyyyyz_z[k] = -g_z_0_x_0_yyyyz_z[k] * ab_y + g_z_0_x_0_yyyyz_yz[k];
            }

            /// Set up 573-576 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 573 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 574 * ccomps * dcomps);

            auto g_z_0_x_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyyzz_x, g_z_0_x_0_yyyyzz_y, g_z_0_x_0_yyyyzz_z, g_z_0_x_0_yyyzz_x, g_z_0_x_0_yyyzz_xy, g_z_0_x_0_yyyzz_y, g_z_0_x_0_yyyzz_yy, g_z_0_x_0_yyyzz_yz, g_z_0_x_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyyzz_x[k] = -g_z_0_x_0_yyyzz_x[k] * ab_y + g_z_0_x_0_yyyzz_xy[k];

                g_z_0_x_0_yyyyzz_y[k] = -g_z_0_x_0_yyyzz_y[k] * ab_y + g_z_0_x_0_yyyzz_yy[k];

                g_z_0_x_0_yyyyzz_z[k] = -g_z_0_x_0_yyyzz_z[k] * ab_y + g_z_0_x_0_yyyzz_yz[k];
            }

            /// Set up 576-579 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 576 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 577 * ccomps * dcomps);

            auto g_z_0_x_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 578 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyzzz_x, g_z_0_x_0_yyyzzz_y, g_z_0_x_0_yyyzzz_z, g_z_0_x_0_yyzzz_x, g_z_0_x_0_yyzzz_xy, g_z_0_x_0_yyzzz_y, g_z_0_x_0_yyzzz_yy, g_z_0_x_0_yyzzz_yz, g_z_0_x_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyzzz_x[k] = -g_z_0_x_0_yyzzz_x[k] * ab_y + g_z_0_x_0_yyzzz_xy[k];

                g_z_0_x_0_yyyzzz_y[k] = -g_z_0_x_0_yyzzz_y[k] * ab_y + g_z_0_x_0_yyzzz_yy[k];

                g_z_0_x_0_yyyzzz_z[k] = -g_z_0_x_0_yyzzz_z[k] * ab_y + g_z_0_x_0_yyzzz_yz[k];
            }

            /// Set up 579-582 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 579 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 580 * ccomps * dcomps);

            auto g_z_0_x_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzzzz_x, g_z_0_x_0_yyzzzz_y, g_z_0_x_0_yyzzzz_z, g_z_0_x_0_yzzzz_x, g_z_0_x_0_yzzzz_xy, g_z_0_x_0_yzzzz_y, g_z_0_x_0_yzzzz_yy, g_z_0_x_0_yzzzz_yz, g_z_0_x_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzzzz_x[k] = -g_z_0_x_0_yzzzz_x[k] * ab_y + g_z_0_x_0_yzzzz_xy[k];

                g_z_0_x_0_yyzzzz_y[k] = -g_z_0_x_0_yzzzz_y[k] * ab_y + g_z_0_x_0_yzzzz_yy[k];

                g_z_0_x_0_yyzzzz_z[k] = -g_z_0_x_0_yzzzz_z[k] * ab_y + g_z_0_x_0_yzzzz_yz[k];
            }

            /// Set up 582-585 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 582 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 583 * ccomps * dcomps);

            auto g_z_0_x_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzzzz_x, g_z_0_x_0_yzzzzz_y, g_z_0_x_0_yzzzzz_z, g_z_0_x_0_zzzzz_x, g_z_0_x_0_zzzzz_xy, g_z_0_x_0_zzzzz_y, g_z_0_x_0_zzzzz_yy, g_z_0_x_0_zzzzz_yz, g_z_0_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzzzz_x[k] = -g_z_0_x_0_zzzzz_x[k] * ab_y + g_z_0_x_0_zzzzz_xy[k];

                g_z_0_x_0_yzzzzz_y[k] = -g_z_0_x_0_zzzzz_y[k] * ab_y + g_z_0_x_0_zzzzz_yy[k];

                g_z_0_x_0_yzzzzz_z[k] = -g_z_0_x_0_zzzzz_z[k] * ab_y + g_z_0_x_0_zzzzz_yz[k];
            }

            /// Set up 585-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 585 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 586 * ccomps * dcomps);

            auto g_z_0_x_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzzzz_x, g_0_0_x_0_zzzzz_y, g_0_0_x_0_zzzzz_z, g_z_0_x_0_zzzzz_x, g_z_0_x_0_zzzzz_xz, g_z_0_x_0_zzzzz_y, g_z_0_x_0_zzzzz_yz, g_z_0_x_0_zzzzz_z, g_z_0_x_0_zzzzz_zz, g_z_0_x_0_zzzzzz_x, g_z_0_x_0_zzzzzz_y, g_z_0_x_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzzzz_x[k] = -g_0_0_x_0_zzzzz_x[k] - g_z_0_x_0_zzzzz_x[k] * ab_z + g_z_0_x_0_zzzzz_xz[k];

                g_z_0_x_0_zzzzzz_y[k] = -g_0_0_x_0_zzzzz_y[k] - g_z_0_x_0_zzzzz_y[k] * ab_z + g_z_0_x_0_zzzzz_yz[k];

                g_z_0_x_0_zzzzzz_z[k] = -g_0_0_x_0_zzzzz_z[k] - g_z_0_x_0_zzzzz_z[k] * ab_z + g_z_0_x_0_zzzzz_zz[k];
            }

            /// Set up 588-591 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 588 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 589 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 590 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxx_x, g_z_0_y_0_xxxxx_xx, g_z_0_y_0_xxxxx_xy, g_z_0_y_0_xxxxx_xz, g_z_0_y_0_xxxxx_y, g_z_0_y_0_xxxxx_z, g_z_0_y_0_xxxxxx_x, g_z_0_y_0_xxxxxx_y, g_z_0_y_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxx_x[k] = -g_z_0_y_0_xxxxx_x[k] * ab_x + g_z_0_y_0_xxxxx_xx[k];

                g_z_0_y_0_xxxxxx_y[k] = -g_z_0_y_0_xxxxx_y[k] * ab_x + g_z_0_y_0_xxxxx_xy[k];

                g_z_0_y_0_xxxxxx_z[k] = -g_z_0_y_0_xxxxx_z[k] * ab_x + g_z_0_y_0_xxxxx_xz[k];
            }

            /// Set up 591-594 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 591 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 592 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxy_x, g_z_0_y_0_xxxxxy_y, g_z_0_y_0_xxxxxy_z, g_z_0_y_0_xxxxy_x, g_z_0_y_0_xxxxy_xx, g_z_0_y_0_xxxxy_xy, g_z_0_y_0_xxxxy_xz, g_z_0_y_0_xxxxy_y, g_z_0_y_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxy_x[k] = -g_z_0_y_0_xxxxy_x[k] * ab_x + g_z_0_y_0_xxxxy_xx[k];

                g_z_0_y_0_xxxxxy_y[k] = -g_z_0_y_0_xxxxy_y[k] * ab_x + g_z_0_y_0_xxxxy_xy[k];

                g_z_0_y_0_xxxxxy_z[k] = -g_z_0_y_0_xxxxy_z[k] * ab_x + g_z_0_y_0_xxxxy_xz[k];
            }

            /// Set up 594-597 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 594 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 595 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 596 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxxz_x, g_z_0_y_0_xxxxxz_y, g_z_0_y_0_xxxxxz_z, g_z_0_y_0_xxxxz_x, g_z_0_y_0_xxxxz_xx, g_z_0_y_0_xxxxz_xy, g_z_0_y_0_xxxxz_xz, g_z_0_y_0_xxxxz_y, g_z_0_y_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxxz_x[k] = -g_z_0_y_0_xxxxz_x[k] * ab_x + g_z_0_y_0_xxxxz_xx[k];

                g_z_0_y_0_xxxxxz_y[k] = -g_z_0_y_0_xxxxz_y[k] * ab_x + g_z_0_y_0_xxxxz_xy[k];

                g_z_0_y_0_xxxxxz_z[k] = -g_z_0_y_0_xxxxz_z[k] * ab_x + g_z_0_y_0_xxxxz_xz[k];
            }

            /// Set up 597-600 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 597 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 598 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyy_x, g_z_0_y_0_xxxxyy_y, g_z_0_y_0_xxxxyy_z, g_z_0_y_0_xxxyy_x, g_z_0_y_0_xxxyy_xx, g_z_0_y_0_xxxyy_xy, g_z_0_y_0_xxxyy_xz, g_z_0_y_0_xxxyy_y, g_z_0_y_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyy_x[k] = -g_z_0_y_0_xxxyy_x[k] * ab_x + g_z_0_y_0_xxxyy_xx[k];

                g_z_0_y_0_xxxxyy_y[k] = -g_z_0_y_0_xxxyy_y[k] * ab_x + g_z_0_y_0_xxxyy_xy[k];

                g_z_0_y_0_xxxxyy_z[k] = -g_z_0_y_0_xxxyy_z[k] * ab_x + g_z_0_y_0_xxxyy_xz[k];
            }

            /// Set up 600-603 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 600 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 601 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 602 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxyz_x, g_z_0_y_0_xxxxyz_y, g_z_0_y_0_xxxxyz_z, g_z_0_y_0_xxxyz_x, g_z_0_y_0_xxxyz_xx, g_z_0_y_0_xxxyz_xy, g_z_0_y_0_xxxyz_xz, g_z_0_y_0_xxxyz_y, g_z_0_y_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxyz_x[k] = -g_z_0_y_0_xxxyz_x[k] * ab_x + g_z_0_y_0_xxxyz_xx[k];

                g_z_0_y_0_xxxxyz_y[k] = -g_z_0_y_0_xxxyz_y[k] * ab_x + g_z_0_y_0_xxxyz_xy[k];

                g_z_0_y_0_xxxxyz_z[k] = -g_z_0_y_0_xxxyz_z[k] * ab_x + g_z_0_y_0_xxxyz_xz[k];
            }

            /// Set up 603-606 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 603 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 604 * ccomps * dcomps);

            auto g_z_0_y_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxxzz_x, g_z_0_y_0_xxxxzz_y, g_z_0_y_0_xxxxzz_z, g_z_0_y_0_xxxzz_x, g_z_0_y_0_xxxzz_xx, g_z_0_y_0_xxxzz_xy, g_z_0_y_0_xxxzz_xz, g_z_0_y_0_xxxzz_y, g_z_0_y_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxxzz_x[k] = -g_z_0_y_0_xxxzz_x[k] * ab_x + g_z_0_y_0_xxxzz_xx[k];

                g_z_0_y_0_xxxxzz_y[k] = -g_z_0_y_0_xxxzz_y[k] * ab_x + g_z_0_y_0_xxxzz_xy[k];

                g_z_0_y_0_xxxxzz_z[k] = -g_z_0_y_0_xxxzz_z[k] * ab_x + g_z_0_y_0_xxxzz_xz[k];
            }

            /// Set up 606-609 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 606 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 607 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyy_x, g_z_0_y_0_xxxyyy_y, g_z_0_y_0_xxxyyy_z, g_z_0_y_0_xxyyy_x, g_z_0_y_0_xxyyy_xx, g_z_0_y_0_xxyyy_xy, g_z_0_y_0_xxyyy_xz, g_z_0_y_0_xxyyy_y, g_z_0_y_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyy_x[k] = -g_z_0_y_0_xxyyy_x[k] * ab_x + g_z_0_y_0_xxyyy_xx[k];

                g_z_0_y_0_xxxyyy_y[k] = -g_z_0_y_0_xxyyy_y[k] * ab_x + g_z_0_y_0_xxyyy_xy[k];

                g_z_0_y_0_xxxyyy_z[k] = -g_z_0_y_0_xxyyy_z[k] * ab_x + g_z_0_y_0_xxyyy_xz[k];
            }

            /// Set up 609-612 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 609 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 610 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyyz_x, g_z_0_y_0_xxxyyz_y, g_z_0_y_0_xxxyyz_z, g_z_0_y_0_xxyyz_x, g_z_0_y_0_xxyyz_xx, g_z_0_y_0_xxyyz_xy, g_z_0_y_0_xxyyz_xz, g_z_0_y_0_xxyyz_y, g_z_0_y_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyyz_x[k] = -g_z_0_y_0_xxyyz_x[k] * ab_x + g_z_0_y_0_xxyyz_xx[k];

                g_z_0_y_0_xxxyyz_y[k] = -g_z_0_y_0_xxyyz_y[k] * ab_x + g_z_0_y_0_xxyyz_xy[k];

                g_z_0_y_0_xxxyyz_z[k] = -g_z_0_y_0_xxyyz_z[k] * ab_x + g_z_0_y_0_xxyyz_xz[k];
            }

            /// Set up 612-615 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 612 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 613 * ccomps * dcomps);

            auto g_z_0_y_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxyzz_x, g_z_0_y_0_xxxyzz_y, g_z_0_y_0_xxxyzz_z, g_z_0_y_0_xxyzz_x, g_z_0_y_0_xxyzz_xx, g_z_0_y_0_xxyzz_xy, g_z_0_y_0_xxyzz_xz, g_z_0_y_0_xxyzz_y, g_z_0_y_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxyzz_x[k] = -g_z_0_y_0_xxyzz_x[k] * ab_x + g_z_0_y_0_xxyzz_xx[k];

                g_z_0_y_0_xxxyzz_y[k] = -g_z_0_y_0_xxyzz_y[k] * ab_x + g_z_0_y_0_xxyzz_xy[k];

                g_z_0_y_0_xxxyzz_z[k] = -g_z_0_y_0_xxyzz_z[k] * ab_x + g_z_0_y_0_xxyzz_xz[k];
            }

            /// Set up 615-618 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 615 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 616 * ccomps * dcomps);

            auto g_z_0_y_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxzzz_x, g_z_0_y_0_xxxzzz_y, g_z_0_y_0_xxxzzz_z, g_z_0_y_0_xxzzz_x, g_z_0_y_0_xxzzz_xx, g_z_0_y_0_xxzzz_xy, g_z_0_y_0_xxzzz_xz, g_z_0_y_0_xxzzz_y, g_z_0_y_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxzzz_x[k] = -g_z_0_y_0_xxzzz_x[k] * ab_x + g_z_0_y_0_xxzzz_xx[k];

                g_z_0_y_0_xxxzzz_y[k] = -g_z_0_y_0_xxzzz_y[k] * ab_x + g_z_0_y_0_xxzzz_xy[k];

                g_z_0_y_0_xxxzzz_z[k] = -g_z_0_y_0_xxzzz_z[k] * ab_x + g_z_0_y_0_xxzzz_xz[k];
            }

            /// Set up 618-621 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 618 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 619 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 620 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyy_x, g_z_0_y_0_xxyyyy_y, g_z_0_y_0_xxyyyy_z, g_z_0_y_0_xyyyy_x, g_z_0_y_0_xyyyy_xx, g_z_0_y_0_xyyyy_xy, g_z_0_y_0_xyyyy_xz, g_z_0_y_0_xyyyy_y, g_z_0_y_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyy_x[k] = -g_z_0_y_0_xyyyy_x[k] * ab_x + g_z_0_y_0_xyyyy_xx[k];

                g_z_0_y_0_xxyyyy_y[k] = -g_z_0_y_0_xyyyy_y[k] * ab_x + g_z_0_y_0_xyyyy_xy[k];

                g_z_0_y_0_xxyyyy_z[k] = -g_z_0_y_0_xyyyy_z[k] * ab_x + g_z_0_y_0_xyyyy_xz[k];
            }

            /// Set up 621-624 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 621 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 622 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyyz_x, g_z_0_y_0_xxyyyz_y, g_z_0_y_0_xxyyyz_z, g_z_0_y_0_xyyyz_x, g_z_0_y_0_xyyyz_xx, g_z_0_y_0_xyyyz_xy, g_z_0_y_0_xyyyz_xz, g_z_0_y_0_xyyyz_y, g_z_0_y_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyyz_x[k] = -g_z_0_y_0_xyyyz_x[k] * ab_x + g_z_0_y_0_xyyyz_xx[k];

                g_z_0_y_0_xxyyyz_y[k] = -g_z_0_y_0_xyyyz_y[k] * ab_x + g_z_0_y_0_xyyyz_xy[k];

                g_z_0_y_0_xxyyyz_z[k] = -g_z_0_y_0_xyyyz_z[k] * ab_x + g_z_0_y_0_xyyyz_xz[k];
            }

            /// Set up 624-627 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 624 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 625 * ccomps * dcomps);

            auto g_z_0_y_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 626 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyyzz_x, g_z_0_y_0_xxyyzz_y, g_z_0_y_0_xxyyzz_z, g_z_0_y_0_xyyzz_x, g_z_0_y_0_xyyzz_xx, g_z_0_y_0_xyyzz_xy, g_z_0_y_0_xyyzz_xz, g_z_0_y_0_xyyzz_y, g_z_0_y_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyyzz_x[k] = -g_z_0_y_0_xyyzz_x[k] * ab_x + g_z_0_y_0_xyyzz_xx[k];

                g_z_0_y_0_xxyyzz_y[k] = -g_z_0_y_0_xyyzz_y[k] * ab_x + g_z_0_y_0_xyyzz_xy[k];

                g_z_0_y_0_xxyyzz_z[k] = -g_z_0_y_0_xyyzz_z[k] * ab_x + g_z_0_y_0_xyyzz_xz[k];
            }

            /// Set up 627-630 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 627 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 628 * ccomps * dcomps);

            auto g_z_0_y_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyzzz_x, g_z_0_y_0_xxyzzz_y, g_z_0_y_0_xxyzzz_z, g_z_0_y_0_xyzzz_x, g_z_0_y_0_xyzzz_xx, g_z_0_y_0_xyzzz_xy, g_z_0_y_0_xyzzz_xz, g_z_0_y_0_xyzzz_y, g_z_0_y_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyzzz_x[k] = -g_z_0_y_0_xyzzz_x[k] * ab_x + g_z_0_y_0_xyzzz_xx[k];

                g_z_0_y_0_xxyzzz_y[k] = -g_z_0_y_0_xyzzz_y[k] * ab_x + g_z_0_y_0_xyzzz_xy[k];

                g_z_0_y_0_xxyzzz_z[k] = -g_z_0_y_0_xyzzz_z[k] * ab_x + g_z_0_y_0_xyzzz_xz[k];
            }

            /// Set up 630-633 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 630 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 631 * ccomps * dcomps);

            auto g_z_0_y_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 632 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzzzz_x, g_z_0_y_0_xxzzzz_y, g_z_0_y_0_xxzzzz_z, g_z_0_y_0_xzzzz_x, g_z_0_y_0_xzzzz_xx, g_z_0_y_0_xzzzz_xy, g_z_0_y_0_xzzzz_xz, g_z_0_y_0_xzzzz_y, g_z_0_y_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzzzz_x[k] = -g_z_0_y_0_xzzzz_x[k] * ab_x + g_z_0_y_0_xzzzz_xx[k];

                g_z_0_y_0_xxzzzz_y[k] = -g_z_0_y_0_xzzzz_y[k] * ab_x + g_z_0_y_0_xzzzz_xy[k];

                g_z_0_y_0_xxzzzz_z[k] = -g_z_0_y_0_xzzzz_z[k] * ab_x + g_z_0_y_0_xzzzz_xz[k];
            }

            /// Set up 633-636 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 633 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 634 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyy_x, g_z_0_y_0_xyyyyy_y, g_z_0_y_0_xyyyyy_z, g_z_0_y_0_yyyyy_x, g_z_0_y_0_yyyyy_xx, g_z_0_y_0_yyyyy_xy, g_z_0_y_0_yyyyy_xz, g_z_0_y_0_yyyyy_y, g_z_0_y_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyy_x[k] = -g_z_0_y_0_yyyyy_x[k] * ab_x + g_z_0_y_0_yyyyy_xx[k];

                g_z_0_y_0_xyyyyy_y[k] = -g_z_0_y_0_yyyyy_y[k] * ab_x + g_z_0_y_0_yyyyy_xy[k];

                g_z_0_y_0_xyyyyy_z[k] = -g_z_0_y_0_yyyyy_z[k] * ab_x + g_z_0_y_0_yyyyy_xz[k];
            }

            /// Set up 636-639 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 636 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 637 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 638 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyyz_x, g_z_0_y_0_xyyyyz_y, g_z_0_y_0_xyyyyz_z, g_z_0_y_0_yyyyz_x, g_z_0_y_0_yyyyz_xx, g_z_0_y_0_yyyyz_xy, g_z_0_y_0_yyyyz_xz, g_z_0_y_0_yyyyz_y, g_z_0_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyyz_x[k] = -g_z_0_y_0_yyyyz_x[k] * ab_x + g_z_0_y_0_yyyyz_xx[k];

                g_z_0_y_0_xyyyyz_y[k] = -g_z_0_y_0_yyyyz_y[k] * ab_x + g_z_0_y_0_yyyyz_xy[k];

                g_z_0_y_0_xyyyyz_z[k] = -g_z_0_y_0_yyyyz_z[k] * ab_x + g_z_0_y_0_yyyyz_xz[k];
            }

            /// Set up 639-642 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 639 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 640 * ccomps * dcomps);

            auto g_z_0_y_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyyzz_x, g_z_0_y_0_xyyyzz_y, g_z_0_y_0_xyyyzz_z, g_z_0_y_0_yyyzz_x, g_z_0_y_0_yyyzz_xx, g_z_0_y_0_yyyzz_xy, g_z_0_y_0_yyyzz_xz, g_z_0_y_0_yyyzz_y, g_z_0_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyyzz_x[k] = -g_z_0_y_0_yyyzz_x[k] * ab_x + g_z_0_y_0_yyyzz_xx[k];

                g_z_0_y_0_xyyyzz_y[k] = -g_z_0_y_0_yyyzz_y[k] * ab_x + g_z_0_y_0_yyyzz_xy[k];

                g_z_0_y_0_xyyyzz_z[k] = -g_z_0_y_0_yyyzz_z[k] * ab_x + g_z_0_y_0_yyyzz_xz[k];
            }

            /// Set up 642-645 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 642 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 643 * ccomps * dcomps);

            auto g_z_0_y_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyzzz_x, g_z_0_y_0_xyyzzz_y, g_z_0_y_0_xyyzzz_z, g_z_0_y_0_yyzzz_x, g_z_0_y_0_yyzzz_xx, g_z_0_y_0_yyzzz_xy, g_z_0_y_0_yyzzz_xz, g_z_0_y_0_yyzzz_y, g_z_0_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyzzz_x[k] = -g_z_0_y_0_yyzzz_x[k] * ab_x + g_z_0_y_0_yyzzz_xx[k];

                g_z_0_y_0_xyyzzz_y[k] = -g_z_0_y_0_yyzzz_y[k] * ab_x + g_z_0_y_0_yyzzz_xy[k];

                g_z_0_y_0_xyyzzz_z[k] = -g_z_0_y_0_yyzzz_z[k] * ab_x + g_z_0_y_0_yyzzz_xz[k];
            }

            /// Set up 645-648 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 645 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 646 * ccomps * dcomps);

            auto g_z_0_y_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzzzz_x, g_z_0_y_0_xyzzzz_y, g_z_0_y_0_xyzzzz_z, g_z_0_y_0_yzzzz_x, g_z_0_y_0_yzzzz_xx, g_z_0_y_0_yzzzz_xy, g_z_0_y_0_yzzzz_xz, g_z_0_y_0_yzzzz_y, g_z_0_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzzzz_x[k] = -g_z_0_y_0_yzzzz_x[k] * ab_x + g_z_0_y_0_yzzzz_xx[k];

                g_z_0_y_0_xyzzzz_y[k] = -g_z_0_y_0_yzzzz_y[k] * ab_x + g_z_0_y_0_yzzzz_xy[k];

                g_z_0_y_0_xyzzzz_z[k] = -g_z_0_y_0_yzzzz_z[k] * ab_x + g_z_0_y_0_yzzzz_xz[k];
            }

            /// Set up 648-651 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 648 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 649 * ccomps * dcomps);

            auto g_z_0_y_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzzzz_x, g_z_0_y_0_xzzzzz_y, g_z_0_y_0_xzzzzz_z, g_z_0_y_0_zzzzz_x, g_z_0_y_0_zzzzz_xx, g_z_0_y_0_zzzzz_xy, g_z_0_y_0_zzzzz_xz, g_z_0_y_0_zzzzz_y, g_z_0_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzzzz_x[k] = -g_z_0_y_0_zzzzz_x[k] * ab_x + g_z_0_y_0_zzzzz_xx[k];

                g_z_0_y_0_xzzzzz_y[k] = -g_z_0_y_0_zzzzz_y[k] * ab_x + g_z_0_y_0_zzzzz_xy[k];

                g_z_0_y_0_xzzzzz_z[k] = -g_z_0_y_0_zzzzz_z[k] * ab_x + g_z_0_y_0_zzzzz_xz[k];
            }

            /// Set up 651-654 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 651 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 652 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyy_x, g_z_0_y_0_yyyyy_xy, g_z_0_y_0_yyyyy_y, g_z_0_y_0_yyyyy_yy, g_z_0_y_0_yyyyy_yz, g_z_0_y_0_yyyyy_z, g_z_0_y_0_yyyyyy_x, g_z_0_y_0_yyyyyy_y, g_z_0_y_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyy_x[k] = -g_z_0_y_0_yyyyy_x[k] * ab_y + g_z_0_y_0_yyyyy_xy[k];

                g_z_0_y_0_yyyyyy_y[k] = -g_z_0_y_0_yyyyy_y[k] * ab_y + g_z_0_y_0_yyyyy_yy[k];

                g_z_0_y_0_yyyyyy_z[k] = -g_z_0_y_0_yyyyy_z[k] * ab_y + g_z_0_y_0_yyyyy_yz[k];
            }

            /// Set up 654-657 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 654 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 655 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 656 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyyz_x, g_z_0_y_0_yyyyyz_y, g_z_0_y_0_yyyyyz_z, g_z_0_y_0_yyyyz_x, g_z_0_y_0_yyyyz_xy, g_z_0_y_0_yyyyz_y, g_z_0_y_0_yyyyz_yy, g_z_0_y_0_yyyyz_yz, g_z_0_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyyz_x[k] = -g_z_0_y_0_yyyyz_x[k] * ab_y + g_z_0_y_0_yyyyz_xy[k];

                g_z_0_y_0_yyyyyz_y[k] = -g_z_0_y_0_yyyyz_y[k] * ab_y + g_z_0_y_0_yyyyz_yy[k];

                g_z_0_y_0_yyyyyz_z[k] = -g_z_0_y_0_yyyyz_z[k] * ab_y + g_z_0_y_0_yyyyz_yz[k];
            }

            /// Set up 657-660 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 657 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 658 * ccomps * dcomps);

            auto g_z_0_y_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyyzz_x, g_z_0_y_0_yyyyzz_y, g_z_0_y_0_yyyyzz_z, g_z_0_y_0_yyyzz_x, g_z_0_y_0_yyyzz_xy, g_z_0_y_0_yyyzz_y, g_z_0_y_0_yyyzz_yy, g_z_0_y_0_yyyzz_yz, g_z_0_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyyzz_x[k] = -g_z_0_y_0_yyyzz_x[k] * ab_y + g_z_0_y_0_yyyzz_xy[k];

                g_z_0_y_0_yyyyzz_y[k] = -g_z_0_y_0_yyyzz_y[k] * ab_y + g_z_0_y_0_yyyzz_yy[k];

                g_z_0_y_0_yyyyzz_z[k] = -g_z_0_y_0_yyyzz_z[k] * ab_y + g_z_0_y_0_yyyzz_yz[k];
            }

            /// Set up 660-663 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 660 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 661 * ccomps * dcomps);

            auto g_z_0_y_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 662 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyzzz_x, g_z_0_y_0_yyyzzz_y, g_z_0_y_0_yyyzzz_z, g_z_0_y_0_yyzzz_x, g_z_0_y_0_yyzzz_xy, g_z_0_y_0_yyzzz_y, g_z_0_y_0_yyzzz_yy, g_z_0_y_0_yyzzz_yz, g_z_0_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyzzz_x[k] = -g_z_0_y_0_yyzzz_x[k] * ab_y + g_z_0_y_0_yyzzz_xy[k];

                g_z_0_y_0_yyyzzz_y[k] = -g_z_0_y_0_yyzzz_y[k] * ab_y + g_z_0_y_0_yyzzz_yy[k];

                g_z_0_y_0_yyyzzz_z[k] = -g_z_0_y_0_yyzzz_z[k] * ab_y + g_z_0_y_0_yyzzz_yz[k];
            }

            /// Set up 663-666 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 663 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 664 * ccomps * dcomps);

            auto g_z_0_y_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzzzz_x, g_z_0_y_0_yyzzzz_y, g_z_0_y_0_yyzzzz_z, g_z_0_y_0_yzzzz_x, g_z_0_y_0_yzzzz_xy, g_z_0_y_0_yzzzz_y, g_z_0_y_0_yzzzz_yy, g_z_0_y_0_yzzzz_yz, g_z_0_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzzzz_x[k] = -g_z_0_y_0_yzzzz_x[k] * ab_y + g_z_0_y_0_yzzzz_xy[k];

                g_z_0_y_0_yyzzzz_y[k] = -g_z_0_y_0_yzzzz_y[k] * ab_y + g_z_0_y_0_yzzzz_yy[k];

                g_z_0_y_0_yyzzzz_z[k] = -g_z_0_y_0_yzzzz_z[k] * ab_y + g_z_0_y_0_yzzzz_yz[k];
            }

            /// Set up 666-669 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 666 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 667 * ccomps * dcomps);

            auto g_z_0_y_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 668 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzzzz_x, g_z_0_y_0_yzzzzz_y, g_z_0_y_0_yzzzzz_z, g_z_0_y_0_zzzzz_x, g_z_0_y_0_zzzzz_xy, g_z_0_y_0_zzzzz_y, g_z_0_y_0_zzzzz_yy, g_z_0_y_0_zzzzz_yz, g_z_0_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzzzz_x[k] = -g_z_0_y_0_zzzzz_x[k] * ab_y + g_z_0_y_0_zzzzz_xy[k];

                g_z_0_y_0_yzzzzz_y[k] = -g_z_0_y_0_zzzzz_y[k] * ab_y + g_z_0_y_0_zzzzz_yy[k];

                g_z_0_y_0_yzzzzz_z[k] = -g_z_0_y_0_zzzzz_z[k] * ab_y + g_z_0_y_0_zzzzz_yz[k];
            }

            /// Set up 669-672 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 669 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 670 * ccomps * dcomps);

            auto g_z_0_y_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzzzz_x, g_0_0_y_0_zzzzz_y, g_0_0_y_0_zzzzz_z, g_z_0_y_0_zzzzz_x, g_z_0_y_0_zzzzz_xz, g_z_0_y_0_zzzzz_y, g_z_0_y_0_zzzzz_yz, g_z_0_y_0_zzzzz_z, g_z_0_y_0_zzzzz_zz, g_z_0_y_0_zzzzzz_x, g_z_0_y_0_zzzzzz_y, g_z_0_y_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzzzz_x[k] = -g_0_0_y_0_zzzzz_x[k] - g_z_0_y_0_zzzzz_x[k] * ab_z + g_z_0_y_0_zzzzz_xz[k];

                g_z_0_y_0_zzzzzz_y[k] = -g_0_0_y_0_zzzzz_y[k] - g_z_0_y_0_zzzzz_y[k] * ab_z + g_z_0_y_0_zzzzz_yz[k];

                g_z_0_y_0_zzzzzz_z[k] = -g_0_0_y_0_zzzzz_z[k] - g_z_0_y_0_zzzzz_z[k] * ab_z + g_z_0_y_0_zzzzz_zz[k];
            }

            /// Set up 672-675 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxx_x = cbuffer.data(ip_geom_1010_off + 672 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_y = cbuffer.data(ip_geom_1010_off + 673 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxx_z = cbuffer.data(ip_geom_1010_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxx_x, g_z_0_z_0_xxxxx_xx, g_z_0_z_0_xxxxx_xy, g_z_0_z_0_xxxxx_xz, g_z_0_z_0_xxxxx_y, g_z_0_z_0_xxxxx_z, g_z_0_z_0_xxxxxx_x, g_z_0_z_0_xxxxxx_y, g_z_0_z_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxx_x[k] = -g_z_0_z_0_xxxxx_x[k] * ab_x + g_z_0_z_0_xxxxx_xx[k];

                g_z_0_z_0_xxxxxx_y[k] = -g_z_0_z_0_xxxxx_y[k] * ab_x + g_z_0_z_0_xxxxx_xy[k];

                g_z_0_z_0_xxxxxx_z[k] = -g_z_0_z_0_xxxxx_z[k] * ab_x + g_z_0_z_0_xxxxx_xz[k];
            }

            /// Set up 675-678 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxy_x = cbuffer.data(ip_geom_1010_off + 675 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_y = cbuffer.data(ip_geom_1010_off + 676 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxy_z = cbuffer.data(ip_geom_1010_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxy_x, g_z_0_z_0_xxxxxy_y, g_z_0_z_0_xxxxxy_z, g_z_0_z_0_xxxxy_x, g_z_0_z_0_xxxxy_xx, g_z_0_z_0_xxxxy_xy, g_z_0_z_0_xxxxy_xz, g_z_0_z_0_xxxxy_y, g_z_0_z_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxy_x[k] = -g_z_0_z_0_xxxxy_x[k] * ab_x + g_z_0_z_0_xxxxy_xx[k];

                g_z_0_z_0_xxxxxy_y[k] = -g_z_0_z_0_xxxxy_y[k] * ab_x + g_z_0_z_0_xxxxy_xy[k];

                g_z_0_z_0_xxxxxy_z[k] = -g_z_0_z_0_xxxxy_z[k] * ab_x + g_z_0_z_0_xxxxy_xz[k];
            }

            /// Set up 678-681 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxxz_x = cbuffer.data(ip_geom_1010_off + 678 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_y = cbuffer.data(ip_geom_1010_off + 679 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxxz_z = cbuffer.data(ip_geom_1010_off + 680 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxxz_x, g_z_0_z_0_xxxxxz_y, g_z_0_z_0_xxxxxz_z, g_z_0_z_0_xxxxz_x, g_z_0_z_0_xxxxz_xx, g_z_0_z_0_xxxxz_xy, g_z_0_z_0_xxxxz_xz, g_z_0_z_0_xxxxz_y, g_z_0_z_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxxz_x[k] = -g_z_0_z_0_xxxxz_x[k] * ab_x + g_z_0_z_0_xxxxz_xx[k];

                g_z_0_z_0_xxxxxz_y[k] = -g_z_0_z_0_xxxxz_y[k] * ab_x + g_z_0_z_0_xxxxz_xy[k];

                g_z_0_z_0_xxxxxz_z[k] = -g_z_0_z_0_xxxxz_z[k] * ab_x + g_z_0_z_0_xxxxz_xz[k];
            }

            /// Set up 681-684 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyy_x = cbuffer.data(ip_geom_1010_off + 681 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_y = cbuffer.data(ip_geom_1010_off + 682 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyy_z = cbuffer.data(ip_geom_1010_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyy_x, g_z_0_z_0_xxxxyy_y, g_z_0_z_0_xxxxyy_z, g_z_0_z_0_xxxyy_x, g_z_0_z_0_xxxyy_xx, g_z_0_z_0_xxxyy_xy, g_z_0_z_0_xxxyy_xz, g_z_0_z_0_xxxyy_y, g_z_0_z_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyy_x[k] = -g_z_0_z_0_xxxyy_x[k] * ab_x + g_z_0_z_0_xxxyy_xx[k];

                g_z_0_z_0_xxxxyy_y[k] = -g_z_0_z_0_xxxyy_y[k] * ab_x + g_z_0_z_0_xxxyy_xy[k];

                g_z_0_z_0_xxxxyy_z[k] = -g_z_0_z_0_xxxyy_z[k] * ab_x + g_z_0_z_0_xxxyy_xz[k];
            }

            /// Set up 684-687 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxyz_x = cbuffer.data(ip_geom_1010_off + 684 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_y = cbuffer.data(ip_geom_1010_off + 685 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxyz_z = cbuffer.data(ip_geom_1010_off + 686 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxyz_x, g_z_0_z_0_xxxxyz_y, g_z_0_z_0_xxxxyz_z, g_z_0_z_0_xxxyz_x, g_z_0_z_0_xxxyz_xx, g_z_0_z_0_xxxyz_xy, g_z_0_z_0_xxxyz_xz, g_z_0_z_0_xxxyz_y, g_z_0_z_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxyz_x[k] = -g_z_0_z_0_xxxyz_x[k] * ab_x + g_z_0_z_0_xxxyz_xx[k];

                g_z_0_z_0_xxxxyz_y[k] = -g_z_0_z_0_xxxyz_y[k] * ab_x + g_z_0_z_0_xxxyz_xy[k];

                g_z_0_z_0_xxxxyz_z[k] = -g_z_0_z_0_xxxyz_z[k] * ab_x + g_z_0_z_0_xxxyz_xz[k];
            }

            /// Set up 687-690 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxxzz_x = cbuffer.data(ip_geom_1010_off + 687 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_y = cbuffer.data(ip_geom_1010_off + 688 * ccomps * dcomps);

            auto g_z_0_z_0_xxxxzz_z = cbuffer.data(ip_geom_1010_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxxzz_x, g_z_0_z_0_xxxxzz_y, g_z_0_z_0_xxxxzz_z, g_z_0_z_0_xxxzz_x, g_z_0_z_0_xxxzz_xx, g_z_0_z_0_xxxzz_xy, g_z_0_z_0_xxxzz_xz, g_z_0_z_0_xxxzz_y, g_z_0_z_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxxzz_x[k] = -g_z_0_z_0_xxxzz_x[k] * ab_x + g_z_0_z_0_xxxzz_xx[k];

                g_z_0_z_0_xxxxzz_y[k] = -g_z_0_z_0_xxxzz_y[k] * ab_x + g_z_0_z_0_xxxzz_xy[k];

                g_z_0_z_0_xxxxzz_z[k] = -g_z_0_z_0_xxxzz_z[k] * ab_x + g_z_0_z_0_xxxzz_xz[k];
            }

            /// Set up 690-693 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyy_x = cbuffer.data(ip_geom_1010_off + 690 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_y = cbuffer.data(ip_geom_1010_off + 691 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyy_z = cbuffer.data(ip_geom_1010_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyy_x, g_z_0_z_0_xxxyyy_y, g_z_0_z_0_xxxyyy_z, g_z_0_z_0_xxyyy_x, g_z_0_z_0_xxyyy_xx, g_z_0_z_0_xxyyy_xy, g_z_0_z_0_xxyyy_xz, g_z_0_z_0_xxyyy_y, g_z_0_z_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyy_x[k] = -g_z_0_z_0_xxyyy_x[k] * ab_x + g_z_0_z_0_xxyyy_xx[k];

                g_z_0_z_0_xxxyyy_y[k] = -g_z_0_z_0_xxyyy_y[k] * ab_x + g_z_0_z_0_xxyyy_xy[k];

                g_z_0_z_0_xxxyyy_z[k] = -g_z_0_z_0_xxyyy_z[k] * ab_x + g_z_0_z_0_xxyyy_xz[k];
            }

            /// Set up 693-696 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyyz_x = cbuffer.data(ip_geom_1010_off + 693 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_y = cbuffer.data(ip_geom_1010_off + 694 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyyz_z = cbuffer.data(ip_geom_1010_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyyz_x, g_z_0_z_0_xxxyyz_y, g_z_0_z_0_xxxyyz_z, g_z_0_z_0_xxyyz_x, g_z_0_z_0_xxyyz_xx, g_z_0_z_0_xxyyz_xy, g_z_0_z_0_xxyyz_xz, g_z_0_z_0_xxyyz_y, g_z_0_z_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyyz_x[k] = -g_z_0_z_0_xxyyz_x[k] * ab_x + g_z_0_z_0_xxyyz_xx[k];

                g_z_0_z_0_xxxyyz_y[k] = -g_z_0_z_0_xxyyz_y[k] * ab_x + g_z_0_z_0_xxyyz_xy[k];

                g_z_0_z_0_xxxyyz_z[k] = -g_z_0_z_0_xxyyz_z[k] * ab_x + g_z_0_z_0_xxyyz_xz[k];
            }

            /// Set up 696-699 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxyzz_x = cbuffer.data(ip_geom_1010_off + 696 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_y = cbuffer.data(ip_geom_1010_off + 697 * ccomps * dcomps);

            auto g_z_0_z_0_xxxyzz_z = cbuffer.data(ip_geom_1010_off + 698 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxyzz_x, g_z_0_z_0_xxxyzz_y, g_z_0_z_0_xxxyzz_z, g_z_0_z_0_xxyzz_x, g_z_0_z_0_xxyzz_xx, g_z_0_z_0_xxyzz_xy, g_z_0_z_0_xxyzz_xz, g_z_0_z_0_xxyzz_y, g_z_0_z_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxyzz_x[k] = -g_z_0_z_0_xxyzz_x[k] * ab_x + g_z_0_z_0_xxyzz_xx[k];

                g_z_0_z_0_xxxyzz_y[k] = -g_z_0_z_0_xxyzz_y[k] * ab_x + g_z_0_z_0_xxyzz_xy[k];

                g_z_0_z_0_xxxyzz_z[k] = -g_z_0_z_0_xxyzz_z[k] * ab_x + g_z_0_z_0_xxyzz_xz[k];
            }

            /// Set up 699-702 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxzzz_x = cbuffer.data(ip_geom_1010_off + 699 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_y = cbuffer.data(ip_geom_1010_off + 700 * ccomps * dcomps);

            auto g_z_0_z_0_xxxzzz_z = cbuffer.data(ip_geom_1010_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxzzz_x, g_z_0_z_0_xxxzzz_y, g_z_0_z_0_xxxzzz_z, g_z_0_z_0_xxzzz_x, g_z_0_z_0_xxzzz_xx, g_z_0_z_0_xxzzz_xy, g_z_0_z_0_xxzzz_xz, g_z_0_z_0_xxzzz_y, g_z_0_z_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxzzz_x[k] = -g_z_0_z_0_xxzzz_x[k] * ab_x + g_z_0_z_0_xxzzz_xx[k];

                g_z_0_z_0_xxxzzz_y[k] = -g_z_0_z_0_xxzzz_y[k] * ab_x + g_z_0_z_0_xxzzz_xy[k];

                g_z_0_z_0_xxxzzz_z[k] = -g_z_0_z_0_xxzzz_z[k] * ab_x + g_z_0_z_0_xxzzz_xz[k];
            }

            /// Set up 702-705 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyy_x = cbuffer.data(ip_geom_1010_off + 702 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_y = cbuffer.data(ip_geom_1010_off + 703 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyy_z = cbuffer.data(ip_geom_1010_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyy_x, g_z_0_z_0_xxyyyy_y, g_z_0_z_0_xxyyyy_z, g_z_0_z_0_xyyyy_x, g_z_0_z_0_xyyyy_xx, g_z_0_z_0_xyyyy_xy, g_z_0_z_0_xyyyy_xz, g_z_0_z_0_xyyyy_y, g_z_0_z_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyy_x[k] = -g_z_0_z_0_xyyyy_x[k] * ab_x + g_z_0_z_0_xyyyy_xx[k];

                g_z_0_z_0_xxyyyy_y[k] = -g_z_0_z_0_xyyyy_y[k] * ab_x + g_z_0_z_0_xyyyy_xy[k];

                g_z_0_z_0_xxyyyy_z[k] = -g_z_0_z_0_xyyyy_z[k] * ab_x + g_z_0_z_0_xyyyy_xz[k];
            }

            /// Set up 705-708 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyyz_x = cbuffer.data(ip_geom_1010_off + 705 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_y = cbuffer.data(ip_geom_1010_off + 706 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyyz_z = cbuffer.data(ip_geom_1010_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyyz_x, g_z_0_z_0_xxyyyz_y, g_z_0_z_0_xxyyyz_z, g_z_0_z_0_xyyyz_x, g_z_0_z_0_xyyyz_xx, g_z_0_z_0_xyyyz_xy, g_z_0_z_0_xyyyz_xz, g_z_0_z_0_xyyyz_y, g_z_0_z_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyyz_x[k] = -g_z_0_z_0_xyyyz_x[k] * ab_x + g_z_0_z_0_xyyyz_xx[k];

                g_z_0_z_0_xxyyyz_y[k] = -g_z_0_z_0_xyyyz_y[k] * ab_x + g_z_0_z_0_xyyyz_xy[k];

                g_z_0_z_0_xxyyyz_z[k] = -g_z_0_z_0_xyyyz_z[k] * ab_x + g_z_0_z_0_xyyyz_xz[k];
            }

            /// Set up 708-711 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyyzz_x = cbuffer.data(ip_geom_1010_off + 708 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_y = cbuffer.data(ip_geom_1010_off + 709 * ccomps * dcomps);

            auto g_z_0_z_0_xxyyzz_z = cbuffer.data(ip_geom_1010_off + 710 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyyzz_x, g_z_0_z_0_xxyyzz_y, g_z_0_z_0_xxyyzz_z, g_z_0_z_0_xyyzz_x, g_z_0_z_0_xyyzz_xx, g_z_0_z_0_xyyzz_xy, g_z_0_z_0_xyyzz_xz, g_z_0_z_0_xyyzz_y, g_z_0_z_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyyzz_x[k] = -g_z_0_z_0_xyyzz_x[k] * ab_x + g_z_0_z_0_xyyzz_xx[k];

                g_z_0_z_0_xxyyzz_y[k] = -g_z_0_z_0_xyyzz_y[k] * ab_x + g_z_0_z_0_xyyzz_xy[k];

                g_z_0_z_0_xxyyzz_z[k] = -g_z_0_z_0_xyyzz_z[k] * ab_x + g_z_0_z_0_xyyzz_xz[k];
            }

            /// Set up 711-714 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyzzz_x = cbuffer.data(ip_geom_1010_off + 711 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_y = cbuffer.data(ip_geom_1010_off + 712 * ccomps * dcomps);

            auto g_z_0_z_0_xxyzzz_z = cbuffer.data(ip_geom_1010_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyzzz_x, g_z_0_z_0_xxyzzz_y, g_z_0_z_0_xxyzzz_z, g_z_0_z_0_xyzzz_x, g_z_0_z_0_xyzzz_xx, g_z_0_z_0_xyzzz_xy, g_z_0_z_0_xyzzz_xz, g_z_0_z_0_xyzzz_y, g_z_0_z_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyzzz_x[k] = -g_z_0_z_0_xyzzz_x[k] * ab_x + g_z_0_z_0_xyzzz_xx[k];

                g_z_0_z_0_xxyzzz_y[k] = -g_z_0_z_0_xyzzz_y[k] * ab_x + g_z_0_z_0_xyzzz_xy[k];

                g_z_0_z_0_xxyzzz_z[k] = -g_z_0_z_0_xyzzz_z[k] * ab_x + g_z_0_z_0_xyzzz_xz[k];
            }

            /// Set up 714-717 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzzzz_x = cbuffer.data(ip_geom_1010_off + 714 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_y = cbuffer.data(ip_geom_1010_off + 715 * ccomps * dcomps);

            auto g_z_0_z_0_xxzzzz_z = cbuffer.data(ip_geom_1010_off + 716 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzzzz_x, g_z_0_z_0_xxzzzz_y, g_z_0_z_0_xxzzzz_z, g_z_0_z_0_xzzzz_x, g_z_0_z_0_xzzzz_xx, g_z_0_z_0_xzzzz_xy, g_z_0_z_0_xzzzz_xz, g_z_0_z_0_xzzzz_y, g_z_0_z_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzzzz_x[k] = -g_z_0_z_0_xzzzz_x[k] * ab_x + g_z_0_z_0_xzzzz_xx[k];

                g_z_0_z_0_xxzzzz_y[k] = -g_z_0_z_0_xzzzz_y[k] * ab_x + g_z_0_z_0_xzzzz_xy[k];

                g_z_0_z_0_xxzzzz_z[k] = -g_z_0_z_0_xzzzz_z[k] * ab_x + g_z_0_z_0_xzzzz_xz[k];
            }

            /// Set up 717-720 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyy_x = cbuffer.data(ip_geom_1010_off + 717 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_y = cbuffer.data(ip_geom_1010_off + 718 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyy_z = cbuffer.data(ip_geom_1010_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyy_x, g_z_0_z_0_xyyyyy_y, g_z_0_z_0_xyyyyy_z, g_z_0_z_0_yyyyy_x, g_z_0_z_0_yyyyy_xx, g_z_0_z_0_yyyyy_xy, g_z_0_z_0_yyyyy_xz, g_z_0_z_0_yyyyy_y, g_z_0_z_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyy_x[k] = -g_z_0_z_0_yyyyy_x[k] * ab_x + g_z_0_z_0_yyyyy_xx[k];

                g_z_0_z_0_xyyyyy_y[k] = -g_z_0_z_0_yyyyy_y[k] * ab_x + g_z_0_z_0_yyyyy_xy[k];

                g_z_0_z_0_xyyyyy_z[k] = -g_z_0_z_0_yyyyy_z[k] * ab_x + g_z_0_z_0_yyyyy_xz[k];
            }

            /// Set up 720-723 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyyz_x = cbuffer.data(ip_geom_1010_off + 720 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_y = cbuffer.data(ip_geom_1010_off + 721 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyyz_z = cbuffer.data(ip_geom_1010_off + 722 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyyz_x, g_z_0_z_0_xyyyyz_y, g_z_0_z_0_xyyyyz_z, g_z_0_z_0_yyyyz_x, g_z_0_z_0_yyyyz_xx, g_z_0_z_0_yyyyz_xy, g_z_0_z_0_yyyyz_xz, g_z_0_z_0_yyyyz_y, g_z_0_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyyz_x[k] = -g_z_0_z_0_yyyyz_x[k] * ab_x + g_z_0_z_0_yyyyz_xx[k];

                g_z_0_z_0_xyyyyz_y[k] = -g_z_0_z_0_yyyyz_y[k] * ab_x + g_z_0_z_0_yyyyz_xy[k];

                g_z_0_z_0_xyyyyz_z[k] = -g_z_0_z_0_yyyyz_z[k] * ab_x + g_z_0_z_0_yyyyz_xz[k];
            }

            /// Set up 723-726 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyyzz_x = cbuffer.data(ip_geom_1010_off + 723 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_y = cbuffer.data(ip_geom_1010_off + 724 * ccomps * dcomps);

            auto g_z_0_z_0_xyyyzz_z = cbuffer.data(ip_geom_1010_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyyzz_x, g_z_0_z_0_xyyyzz_y, g_z_0_z_0_xyyyzz_z, g_z_0_z_0_yyyzz_x, g_z_0_z_0_yyyzz_xx, g_z_0_z_0_yyyzz_xy, g_z_0_z_0_yyyzz_xz, g_z_0_z_0_yyyzz_y, g_z_0_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyyzz_x[k] = -g_z_0_z_0_yyyzz_x[k] * ab_x + g_z_0_z_0_yyyzz_xx[k];

                g_z_0_z_0_xyyyzz_y[k] = -g_z_0_z_0_yyyzz_y[k] * ab_x + g_z_0_z_0_yyyzz_xy[k];

                g_z_0_z_0_xyyyzz_z[k] = -g_z_0_z_0_yyyzz_z[k] * ab_x + g_z_0_z_0_yyyzz_xz[k];
            }

            /// Set up 726-729 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyzzz_x = cbuffer.data(ip_geom_1010_off + 726 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_y = cbuffer.data(ip_geom_1010_off + 727 * ccomps * dcomps);

            auto g_z_0_z_0_xyyzzz_z = cbuffer.data(ip_geom_1010_off + 728 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyzzz_x, g_z_0_z_0_xyyzzz_y, g_z_0_z_0_xyyzzz_z, g_z_0_z_0_yyzzz_x, g_z_0_z_0_yyzzz_xx, g_z_0_z_0_yyzzz_xy, g_z_0_z_0_yyzzz_xz, g_z_0_z_0_yyzzz_y, g_z_0_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyzzz_x[k] = -g_z_0_z_0_yyzzz_x[k] * ab_x + g_z_0_z_0_yyzzz_xx[k];

                g_z_0_z_0_xyyzzz_y[k] = -g_z_0_z_0_yyzzz_y[k] * ab_x + g_z_0_z_0_yyzzz_xy[k];

                g_z_0_z_0_xyyzzz_z[k] = -g_z_0_z_0_yyzzz_z[k] * ab_x + g_z_0_z_0_yyzzz_xz[k];
            }

            /// Set up 729-732 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzzzz_x = cbuffer.data(ip_geom_1010_off + 729 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_y = cbuffer.data(ip_geom_1010_off + 730 * ccomps * dcomps);

            auto g_z_0_z_0_xyzzzz_z = cbuffer.data(ip_geom_1010_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzzzz_x, g_z_0_z_0_xyzzzz_y, g_z_0_z_0_xyzzzz_z, g_z_0_z_0_yzzzz_x, g_z_0_z_0_yzzzz_xx, g_z_0_z_0_yzzzz_xy, g_z_0_z_0_yzzzz_xz, g_z_0_z_0_yzzzz_y, g_z_0_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzzzz_x[k] = -g_z_0_z_0_yzzzz_x[k] * ab_x + g_z_0_z_0_yzzzz_xx[k];

                g_z_0_z_0_xyzzzz_y[k] = -g_z_0_z_0_yzzzz_y[k] * ab_x + g_z_0_z_0_yzzzz_xy[k];

                g_z_0_z_0_xyzzzz_z[k] = -g_z_0_z_0_yzzzz_z[k] * ab_x + g_z_0_z_0_yzzzz_xz[k];
            }

            /// Set up 732-735 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzzzz_x = cbuffer.data(ip_geom_1010_off + 732 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_y = cbuffer.data(ip_geom_1010_off + 733 * ccomps * dcomps);

            auto g_z_0_z_0_xzzzzz_z = cbuffer.data(ip_geom_1010_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzzzz_x, g_z_0_z_0_xzzzzz_y, g_z_0_z_0_xzzzzz_z, g_z_0_z_0_zzzzz_x, g_z_0_z_0_zzzzz_xx, g_z_0_z_0_zzzzz_xy, g_z_0_z_0_zzzzz_xz, g_z_0_z_0_zzzzz_y, g_z_0_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzzzz_x[k] = -g_z_0_z_0_zzzzz_x[k] * ab_x + g_z_0_z_0_zzzzz_xx[k];

                g_z_0_z_0_xzzzzz_y[k] = -g_z_0_z_0_zzzzz_y[k] * ab_x + g_z_0_z_0_zzzzz_xy[k];

                g_z_0_z_0_xzzzzz_z[k] = -g_z_0_z_0_zzzzz_z[k] * ab_x + g_z_0_z_0_zzzzz_xz[k];
            }

            /// Set up 735-738 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyy_x = cbuffer.data(ip_geom_1010_off + 735 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_y = cbuffer.data(ip_geom_1010_off + 736 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyy_z = cbuffer.data(ip_geom_1010_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyy_x, g_z_0_z_0_yyyyy_xy, g_z_0_z_0_yyyyy_y, g_z_0_z_0_yyyyy_yy, g_z_0_z_0_yyyyy_yz, g_z_0_z_0_yyyyy_z, g_z_0_z_0_yyyyyy_x, g_z_0_z_0_yyyyyy_y, g_z_0_z_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyy_x[k] = -g_z_0_z_0_yyyyy_x[k] * ab_y + g_z_0_z_0_yyyyy_xy[k];

                g_z_0_z_0_yyyyyy_y[k] = -g_z_0_z_0_yyyyy_y[k] * ab_y + g_z_0_z_0_yyyyy_yy[k];

                g_z_0_z_0_yyyyyy_z[k] = -g_z_0_z_0_yyyyy_z[k] * ab_y + g_z_0_z_0_yyyyy_yz[k];
            }

            /// Set up 738-741 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyyz_x = cbuffer.data(ip_geom_1010_off + 738 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_y = cbuffer.data(ip_geom_1010_off + 739 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyyz_z = cbuffer.data(ip_geom_1010_off + 740 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyyz_x, g_z_0_z_0_yyyyyz_y, g_z_0_z_0_yyyyyz_z, g_z_0_z_0_yyyyz_x, g_z_0_z_0_yyyyz_xy, g_z_0_z_0_yyyyz_y, g_z_0_z_0_yyyyz_yy, g_z_0_z_0_yyyyz_yz, g_z_0_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyyz_x[k] = -g_z_0_z_0_yyyyz_x[k] * ab_y + g_z_0_z_0_yyyyz_xy[k];

                g_z_0_z_0_yyyyyz_y[k] = -g_z_0_z_0_yyyyz_y[k] * ab_y + g_z_0_z_0_yyyyz_yy[k];

                g_z_0_z_0_yyyyyz_z[k] = -g_z_0_z_0_yyyyz_z[k] * ab_y + g_z_0_z_0_yyyyz_yz[k];
            }

            /// Set up 741-744 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyyzz_x = cbuffer.data(ip_geom_1010_off + 741 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_y = cbuffer.data(ip_geom_1010_off + 742 * ccomps * dcomps);

            auto g_z_0_z_0_yyyyzz_z = cbuffer.data(ip_geom_1010_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyyzz_x, g_z_0_z_0_yyyyzz_y, g_z_0_z_0_yyyyzz_z, g_z_0_z_0_yyyzz_x, g_z_0_z_0_yyyzz_xy, g_z_0_z_0_yyyzz_y, g_z_0_z_0_yyyzz_yy, g_z_0_z_0_yyyzz_yz, g_z_0_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyyzz_x[k] = -g_z_0_z_0_yyyzz_x[k] * ab_y + g_z_0_z_0_yyyzz_xy[k];

                g_z_0_z_0_yyyyzz_y[k] = -g_z_0_z_0_yyyzz_y[k] * ab_y + g_z_0_z_0_yyyzz_yy[k];

                g_z_0_z_0_yyyyzz_z[k] = -g_z_0_z_0_yyyzz_z[k] * ab_y + g_z_0_z_0_yyyzz_yz[k];
            }

            /// Set up 744-747 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyzzz_x = cbuffer.data(ip_geom_1010_off + 744 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_y = cbuffer.data(ip_geom_1010_off + 745 * ccomps * dcomps);

            auto g_z_0_z_0_yyyzzz_z = cbuffer.data(ip_geom_1010_off + 746 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyzzz_x, g_z_0_z_0_yyyzzz_y, g_z_0_z_0_yyyzzz_z, g_z_0_z_0_yyzzz_x, g_z_0_z_0_yyzzz_xy, g_z_0_z_0_yyzzz_y, g_z_0_z_0_yyzzz_yy, g_z_0_z_0_yyzzz_yz, g_z_0_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyzzz_x[k] = -g_z_0_z_0_yyzzz_x[k] * ab_y + g_z_0_z_0_yyzzz_xy[k];

                g_z_0_z_0_yyyzzz_y[k] = -g_z_0_z_0_yyzzz_y[k] * ab_y + g_z_0_z_0_yyzzz_yy[k];

                g_z_0_z_0_yyyzzz_z[k] = -g_z_0_z_0_yyzzz_z[k] * ab_y + g_z_0_z_0_yyzzz_yz[k];
            }

            /// Set up 747-750 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzzzz_x = cbuffer.data(ip_geom_1010_off + 747 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_y = cbuffer.data(ip_geom_1010_off + 748 * ccomps * dcomps);

            auto g_z_0_z_0_yyzzzz_z = cbuffer.data(ip_geom_1010_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzzzz_x, g_z_0_z_0_yyzzzz_y, g_z_0_z_0_yyzzzz_z, g_z_0_z_0_yzzzz_x, g_z_0_z_0_yzzzz_xy, g_z_0_z_0_yzzzz_y, g_z_0_z_0_yzzzz_yy, g_z_0_z_0_yzzzz_yz, g_z_0_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzzzz_x[k] = -g_z_0_z_0_yzzzz_x[k] * ab_y + g_z_0_z_0_yzzzz_xy[k];

                g_z_0_z_0_yyzzzz_y[k] = -g_z_0_z_0_yzzzz_y[k] * ab_y + g_z_0_z_0_yzzzz_yy[k];

                g_z_0_z_0_yyzzzz_z[k] = -g_z_0_z_0_yzzzz_z[k] * ab_y + g_z_0_z_0_yzzzz_yz[k];
            }

            /// Set up 750-753 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzzzz_x = cbuffer.data(ip_geom_1010_off + 750 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_y = cbuffer.data(ip_geom_1010_off + 751 * ccomps * dcomps);

            auto g_z_0_z_0_yzzzzz_z = cbuffer.data(ip_geom_1010_off + 752 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzzzz_x, g_z_0_z_0_yzzzzz_y, g_z_0_z_0_yzzzzz_z, g_z_0_z_0_zzzzz_x, g_z_0_z_0_zzzzz_xy, g_z_0_z_0_zzzzz_y, g_z_0_z_0_zzzzz_yy, g_z_0_z_0_zzzzz_yz, g_z_0_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzzzz_x[k] = -g_z_0_z_0_zzzzz_x[k] * ab_y + g_z_0_z_0_zzzzz_xy[k];

                g_z_0_z_0_yzzzzz_y[k] = -g_z_0_z_0_zzzzz_y[k] * ab_y + g_z_0_z_0_zzzzz_yy[k];

                g_z_0_z_0_yzzzzz_z[k] = -g_z_0_z_0_zzzzz_z[k] * ab_y + g_z_0_z_0_zzzzz_yz[k];
            }

            /// Set up 753-756 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzzzz_x = cbuffer.data(ip_geom_1010_off + 753 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_y = cbuffer.data(ip_geom_1010_off + 754 * ccomps * dcomps);

            auto g_z_0_z_0_zzzzzz_z = cbuffer.data(ip_geom_1010_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzzzz_x, g_0_0_z_0_zzzzz_y, g_0_0_z_0_zzzzz_z, g_z_0_z_0_zzzzz_x, g_z_0_z_0_zzzzz_xz, g_z_0_z_0_zzzzz_y, g_z_0_z_0_zzzzz_yz, g_z_0_z_0_zzzzz_z, g_z_0_z_0_zzzzz_zz, g_z_0_z_0_zzzzzz_x, g_z_0_z_0_zzzzzz_y, g_z_0_z_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzzzz_x[k] = -g_0_0_z_0_zzzzz_x[k] - g_z_0_z_0_zzzzz_x[k] * ab_z + g_z_0_z_0_zzzzz_xz[k];

                g_z_0_z_0_zzzzzz_y[k] = -g_0_0_z_0_zzzzz_y[k] - g_z_0_z_0_zzzzz_y[k] * ab_z + g_z_0_z_0_zzzzz_yz[k];

                g_z_0_z_0_zzzzzz_z[k] = -g_0_0_z_0_zzzzz_z[k] - g_z_0_z_0_zzzzz_z[k] * ab_z + g_z_0_z_0_zzzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

