#include "ElectronRepulsionGeom1010ContrRecFHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_fhxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_fhxx,
                                              const size_t idx_geom_0010_dhxx,
                                              const size_t idx_geom_1010_dhxx,
                                              const size_t idx_geom_1010_dixx,
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
            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_0010_off = idx_geom_0010_dhxx + i * dcomps + j;

            auto g_0_0_x_0_xx_xxxxx = cbuffer.data(dh_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxxxy = cbuffer.data(dh_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxxxz = cbuffer.data(dh_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxxyy = cbuffer.data(dh_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxxyz = cbuffer.data(dh_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxxzz = cbuffer.data(dh_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxyyy = cbuffer.data(dh_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxyyz = cbuffer.data(dh_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxyzz = cbuffer.data(dh_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xxzzz = cbuffer.data(dh_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xyyyy = cbuffer.data(dh_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xyyyz = cbuffer.data(dh_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xyyzz = cbuffer.data(dh_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xyzzz = cbuffer.data(dh_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xx_xzzzz = cbuffer.data(dh_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_xx_yyyyy = cbuffer.data(dh_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_xx_yyyyz = cbuffer.data(dh_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_xx_yyyzz = cbuffer.data(dh_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_xx_yyzzz = cbuffer.data(dh_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_xx_yzzzz = cbuffer.data(dh_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_xx_zzzzz = cbuffer.data(dh_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxxx = cbuffer.data(dh_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxxy = cbuffer.data(dh_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxxz = cbuffer.data(dh_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxyy = cbuffer.data(dh_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxyz = cbuffer.data(dh_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxxzz = cbuffer.data(dh_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxyyy = cbuffer.data(dh_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxyyz = cbuffer.data(dh_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxyzz = cbuffer.data(dh_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xxzzz = cbuffer.data(dh_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xyyyy = cbuffer.data(dh_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xyyyz = cbuffer.data(dh_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xyyzz = cbuffer.data(dh_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xyzzz = cbuffer.data(dh_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_x_0_xy_xzzzz = cbuffer.data(dh_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_x_0_xy_yyyyy = cbuffer.data(dh_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_x_0_xy_yyyyz = cbuffer.data(dh_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_x_0_xy_yyyzz = cbuffer.data(dh_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_x_0_xy_yyzzz = cbuffer.data(dh_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_x_0_xy_yzzzz = cbuffer.data(dh_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_x_0_xy_zzzzz = cbuffer.data(dh_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxxx = cbuffer.data(dh_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxxy = cbuffer.data(dh_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxxz = cbuffer.data(dh_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxyy = cbuffer.data(dh_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxyz = cbuffer.data(dh_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxxzz = cbuffer.data(dh_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxyyy = cbuffer.data(dh_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxyyz = cbuffer.data(dh_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxyzz = cbuffer.data(dh_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xxzzz = cbuffer.data(dh_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xyyyy = cbuffer.data(dh_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xyyyz = cbuffer.data(dh_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xyyzz = cbuffer.data(dh_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xyzzz = cbuffer.data(dh_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_x_0_xz_xzzzz = cbuffer.data(dh_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_x_0_xz_yyyyy = cbuffer.data(dh_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_x_0_xz_yyyyz = cbuffer.data(dh_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_x_0_xz_yyyzz = cbuffer.data(dh_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_x_0_xz_yyzzz = cbuffer.data(dh_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_x_0_xz_yzzzz = cbuffer.data(dh_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_x_0_xz_zzzzz = cbuffer.data(dh_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxxx = cbuffer.data(dh_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxxy = cbuffer.data(dh_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxxz = cbuffer.data(dh_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxyy = cbuffer.data(dh_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxyz = cbuffer.data(dh_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxxzz = cbuffer.data(dh_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxyyy = cbuffer.data(dh_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxyyz = cbuffer.data(dh_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxyzz = cbuffer.data(dh_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xxzzz = cbuffer.data(dh_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xyyyy = cbuffer.data(dh_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xyyyz = cbuffer.data(dh_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xyyzz = cbuffer.data(dh_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xyzzz = cbuffer.data(dh_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_x_0_yy_xzzzz = cbuffer.data(dh_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_x_0_yy_yyyyy = cbuffer.data(dh_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_x_0_yy_yyyyz = cbuffer.data(dh_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_x_0_yy_yyyzz = cbuffer.data(dh_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_x_0_yy_yyzzz = cbuffer.data(dh_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_x_0_yy_yzzzz = cbuffer.data(dh_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_x_0_yy_zzzzz = cbuffer.data(dh_geom_0010_off + 83 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxxx = cbuffer.data(dh_geom_0010_off + 84 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxxy = cbuffer.data(dh_geom_0010_off + 85 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxxz = cbuffer.data(dh_geom_0010_off + 86 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxyy = cbuffer.data(dh_geom_0010_off + 87 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxyz = cbuffer.data(dh_geom_0010_off + 88 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxxzz = cbuffer.data(dh_geom_0010_off + 89 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxyyy = cbuffer.data(dh_geom_0010_off + 90 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxyyz = cbuffer.data(dh_geom_0010_off + 91 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxyzz = cbuffer.data(dh_geom_0010_off + 92 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xxzzz = cbuffer.data(dh_geom_0010_off + 93 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xyyyy = cbuffer.data(dh_geom_0010_off + 94 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xyyyz = cbuffer.data(dh_geom_0010_off + 95 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xyyzz = cbuffer.data(dh_geom_0010_off + 96 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xyzzz = cbuffer.data(dh_geom_0010_off + 97 * ccomps * dcomps);

            auto g_0_0_x_0_yz_xzzzz = cbuffer.data(dh_geom_0010_off + 98 * ccomps * dcomps);

            auto g_0_0_x_0_yz_yyyyy = cbuffer.data(dh_geom_0010_off + 99 * ccomps * dcomps);

            auto g_0_0_x_0_yz_yyyyz = cbuffer.data(dh_geom_0010_off + 100 * ccomps * dcomps);

            auto g_0_0_x_0_yz_yyyzz = cbuffer.data(dh_geom_0010_off + 101 * ccomps * dcomps);

            auto g_0_0_x_0_yz_yyzzz = cbuffer.data(dh_geom_0010_off + 102 * ccomps * dcomps);

            auto g_0_0_x_0_yz_yzzzz = cbuffer.data(dh_geom_0010_off + 103 * ccomps * dcomps);

            auto g_0_0_x_0_yz_zzzzz = cbuffer.data(dh_geom_0010_off + 104 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxxx = cbuffer.data(dh_geom_0010_off + 105 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxxy = cbuffer.data(dh_geom_0010_off + 106 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxxz = cbuffer.data(dh_geom_0010_off + 107 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxyy = cbuffer.data(dh_geom_0010_off + 108 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxyz = cbuffer.data(dh_geom_0010_off + 109 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxxzz = cbuffer.data(dh_geom_0010_off + 110 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxyyy = cbuffer.data(dh_geom_0010_off + 111 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxyyz = cbuffer.data(dh_geom_0010_off + 112 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxyzz = cbuffer.data(dh_geom_0010_off + 113 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xxzzz = cbuffer.data(dh_geom_0010_off + 114 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xyyyy = cbuffer.data(dh_geom_0010_off + 115 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xyyyz = cbuffer.data(dh_geom_0010_off + 116 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xyyzz = cbuffer.data(dh_geom_0010_off + 117 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xyzzz = cbuffer.data(dh_geom_0010_off + 118 * ccomps * dcomps);

            auto g_0_0_x_0_zz_xzzzz = cbuffer.data(dh_geom_0010_off + 119 * ccomps * dcomps);

            auto g_0_0_x_0_zz_yyyyy = cbuffer.data(dh_geom_0010_off + 120 * ccomps * dcomps);

            auto g_0_0_x_0_zz_yyyyz = cbuffer.data(dh_geom_0010_off + 121 * ccomps * dcomps);

            auto g_0_0_x_0_zz_yyyzz = cbuffer.data(dh_geom_0010_off + 122 * ccomps * dcomps);

            auto g_0_0_x_0_zz_yyzzz = cbuffer.data(dh_geom_0010_off + 123 * ccomps * dcomps);

            auto g_0_0_x_0_zz_yzzzz = cbuffer.data(dh_geom_0010_off + 124 * ccomps * dcomps);

            auto g_0_0_x_0_zz_zzzzz = cbuffer.data(dh_geom_0010_off + 125 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxxx = cbuffer.data(dh_geom_0010_off + 126 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxxy = cbuffer.data(dh_geom_0010_off + 127 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxxz = cbuffer.data(dh_geom_0010_off + 128 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxyy = cbuffer.data(dh_geom_0010_off + 129 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxyz = cbuffer.data(dh_geom_0010_off + 130 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxxzz = cbuffer.data(dh_geom_0010_off + 131 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxyyy = cbuffer.data(dh_geom_0010_off + 132 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxyyz = cbuffer.data(dh_geom_0010_off + 133 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxyzz = cbuffer.data(dh_geom_0010_off + 134 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xxzzz = cbuffer.data(dh_geom_0010_off + 135 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xyyyy = cbuffer.data(dh_geom_0010_off + 136 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xyyyz = cbuffer.data(dh_geom_0010_off + 137 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xyyzz = cbuffer.data(dh_geom_0010_off + 138 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xyzzz = cbuffer.data(dh_geom_0010_off + 139 * ccomps * dcomps);

            auto g_0_0_y_0_xx_xzzzz = cbuffer.data(dh_geom_0010_off + 140 * ccomps * dcomps);

            auto g_0_0_y_0_xx_yyyyy = cbuffer.data(dh_geom_0010_off + 141 * ccomps * dcomps);

            auto g_0_0_y_0_xx_yyyyz = cbuffer.data(dh_geom_0010_off + 142 * ccomps * dcomps);

            auto g_0_0_y_0_xx_yyyzz = cbuffer.data(dh_geom_0010_off + 143 * ccomps * dcomps);

            auto g_0_0_y_0_xx_yyzzz = cbuffer.data(dh_geom_0010_off + 144 * ccomps * dcomps);

            auto g_0_0_y_0_xx_yzzzz = cbuffer.data(dh_geom_0010_off + 145 * ccomps * dcomps);

            auto g_0_0_y_0_xx_zzzzz = cbuffer.data(dh_geom_0010_off + 146 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxxx = cbuffer.data(dh_geom_0010_off + 147 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxxy = cbuffer.data(dh_geom_0010_off + 148 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxxz = cbuffer.data(dh_geom_0010_off + 149 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxyy = cbuffer.data(dh_geom_0010_off + 150 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxyz = cbuffer.data(dh_geom_0010_off + 151 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxxzz = cbuffer.data(dh_geom_0010_off + 152 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxyyy = cbuffer.data(dh_geom_0010_off + 153 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxyyz = cbuffer.data(dh_geom_0010_off + 154 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxyzz = cbuffer.data(dh_geom_0010_off + 155 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xxzzz = cbuffer.data(dh_geom_0010_off + 156 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xyyyy = cbuffer.data(dh_geom_0010_off + 157 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xyyyz = cbuffer.data(dh_geom_0010_off + 158 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xyyzz = cbuffer.data(dh_geom_0010_off + 159 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xyzzz = cbuffer.data(dh_geom_0010_off + 160 * ccomps * dcomps);

            auto g_0_0_y_0_xy_xzzzz = cbuffer.data(dh_geom_0010_off + 161 * ccomps * dcomps);

            auto g_0_0_y_0_xy_yyyyy = cbuffer.data(dh_geom_0010_off + 162 * ccomps * dcomps);

            auto g_0_0_y_0_xy_yyyyz = cbuffer.data(dh_geom_0010_off + 163 * ccomps * dcomps);

            auto g_0_0_y_0_xy_yyyzz = cbuffer.data(dh_geom_0010_off + 164 * ccomps * dcomps);

            auto g_0_0_y_0_xy_yyzzz = cbuffer.data(dh_geom_0010_off + 165 * ccomps * dcomps);

            auto g_0_0_y_0_xy_yzzzz = cbuffer.data(dh_geom_0010_off + 166 * ccomps * dcomps);

            auto g_0_0_y_0_xy_zzzzz = cbuffer.data(dh_geom_0010_off + 167 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxxx = cbuffer.data(dh_geom_0010_off + 168 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxxy = cbuffer.data(dh_geom_0010_off + 169 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxxz = cbuffer.data(dh_geom_0010_off + 170 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxyy = cbuffer.data(dh_geom_0010_off + 171 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxyz = cbuffer.data(dh_geom_0010_off + 172 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxxzz = cbuffer.data(dh_geom_0010_off + 173 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxyyy = cbuffer.data(dh_geom_0010_off + 174 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxyyz = cbuffer.data(dh_geom_0010_off + 175 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxyzz = cbuffer.data(dh_geom_0010_off + 176 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xxzzz = cbuffer.data(dh_geom_0010_off + 177 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xyyyy = cbuffer.data(dh_geom_0010_off + 178 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xyyyz = cbuffer.data(dh_geom_0010_off + 179 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xyyzz = cbuffer.data(dh_geom_0010_off + 180 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xyzzz = cbuffer.data(dh_geom_0010_off + 181 * ccomps * dcomps);

            auto g_0_0_y_0_xz_xzzzz = cbuffer.data(dh_geom_0010_off + 182 * ccomps * dcomps);

            auto g_0_0_y_0_xz_yyyyy = cbuffer.data(dh_geom_0010_off + 183 * ccomps * dcomps);

            auto g_0_0_y_0_xz_yyyyz = cbuffer.data(dh_geom_0010_off + 184 * ccomps * dcomps);

            auto g_0_0_y_0_xz_yyyzz = cbuffer.data(dh_geom_0010_off + 185 * ccomps * dcomps);

            auto g_0_0_y_0_xz_yyzzz = cbuffer.data(dh_geom_0010_off + 186 * ccomps * dcomps);

            auto g_0_0_y_0_xz_yzzzz = cbuffer.data(dh_geom_0010_off + 187 * ccomps * dcomps);

            auto g_0_0_y_0_xz_zzzzz = cbuffer.data(dh_geom_0010_off + 188 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxxx = cbuffer.data(dh_geom_0010_off + 189 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxxy = cbuffer.data(dh_geom_0010_off + 190 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxxz = cbuffer.data(dh_geom_0010_off + 191 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxyy = cbuffer.data(dh_geom_0010_off + 192 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxyz = cbuffer.data(dh_geom_0010_off + 193 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxxzz = cbuffer.data(dh_geom_0010_off + 194 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxyyy = cbuffer.data(dh_geom_0010_off + 195 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxyyz = cbuffer.data(dh_geom_0010_off + 196 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxyzz = cbuffer.data(dh_geom_0010_off + 197 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xxzzz = cbuffer.data(dh_geom_0010_off + 198 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xyyyy = cbuffer.data(dh_geom_0010_off + 199 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xyyyz = cbuffer.data(dh_geom_0010_off + 200 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xyyzz = cbuffer.data(dh_geom_0010_off + 201 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xyzzz = cbuffer.data(dh_geom_0010_off + 202 * ccomps * dcomps);

            auto g_0_0_y_0_yy_xzzzz = cbuffer.data(dh_geom_0010_off + 203 * ccomps * dcomps);

            auto g_0_0_y_0_yy_yyyyy = cbuffer.data(dh_geom_0010_off + 204 * ccomps * dcomps);

            auto g_0_0_y_0_yy_yyyyz = cbuffer.data(dh_geom_0010_off + 205 * ccomps * dcomps);

            auto g_0_0_y_0_yy_yyyzz = cbuffer.data(dh_geom_0010_off + 206 * ccomps * dcomps);

            auto g_0_0_y_0_yy_yyzzz = cbuffer.data(dh_geom_0010_off + 207 * ccomps * dcomps);

            auto g_0_0_y_0_yy_yzzzz = cbuffer.data(dh_geom_0010_off + 208 * ccomps * dcomps);

            auto g_0_0_y_0_yy_zzzzz = cbuffer.data(dh_geom_0010_off + 209 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxxx = cbuffer.data(dh_geom_0010_off + 210 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxxy = cbuffer.data(dh_geom_0010_off + 211 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxxz = cbuffer.data(dh_geom_0010_off + 212 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxyy = cbuffer.data(dh_geom_0010_off + 213 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxyz = cbuffer.data(dh_geom_0010_off + 214 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxxzz = cbuffer.data(dh_geom_0010_off + 215 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxyyy = cbuffer.data(dh_geom_0010_off + 216 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxyyz = cbuffer.data(dh_geom_0010_off + 217 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxyzz = cbuffer.data(dh_geom_0010_off + 218 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xxzzz = cbuffer.data(dh_geom_0010_off + 219 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xyyyy = cbuffer.data(dh_geom_0010_off + 220 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xyyyz = cbuffer.data(dh_geom_0010_off + 221 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xyyzz = cbuffer.data(dh_geom_0010_off + 222 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xyzzz = cbuffer.data(dh_geom_0010_off + 223 * ccomps * dcomps);

            auto g_0_0_y_0_yz_xzzzz = cbuffer.data(dh_geom_0010_off + 224 * ccomps * dcomps);

            auto g_0_0_y_0_yz_yyyyy = cbuffer.data(dh_geom_0010_off + 225 * ccomps * dcomps);

            auto g_0_0_y_0_yz_yyyyz = cbuffer.data(dh_geom_0010_off + 226 * ccomps * dcomps);

            auto g_0_0_y_0_yz_yyyzz = cbuffer.data(dh_geom_0010_off + 227 * ccomps * dcomps);

            auto g_0_0_y_0_yz_yyzzz = cbuffer.data(dh_geom_0010_off + 228 * ccomps * dcomps);

            auto g_0_0_y_0_yz_yzzzz = cbuffer.data(dh_geom_0010_off + 229 * ccomps * dcomps);

            auto g_0_0_y_0_yz_zzzzz = cbuffer.data(dh_geom_0010_off + 230 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxxx = cbuffer.data(dh_geom_0010_off + 231 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxxy = cbuffer.data(dh_geom_0010_off + 232 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxxz = cbuffer.data(dh_geom_0010_off + 233 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxyy = cbuffer.data(dh_geom_0010_off + 234 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxyz = cbuffer.data(dh_geom_0010_off + 235 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxxzz = cbuffer.data(dh_geom_0010_off + 236 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxyyy = cbuffer.data(dh_geom_0010_off + 237 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxyyz = cbuffer.data(dh_geom_0010_off + 238 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxyzz = cbuffer.data(dh_geom_0010_off + 239 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xxzzz = cbuffer.data(dh_geom_0010_off + 240 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xyyyy = cbuffer.data(dh_geom_0010_off + 241 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xyyyz = cbuffer.data(dh_geom_0010_off + 242 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xyyzz = cbuffer.data(dh_geom_0010_off + 243 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xyzzz = cbuffer.data(dh_geom_0010_off + 244 * ccomps * dcomps);

            auto g_0_0_y_0_zz_xzzzz = cbuffer.data(dh_geom_0010_off + 245 * ccomps * dcomps);

            auto g_0_0_y_0_zz_yyyyy = cbuffer.data(dh_geom_0010_off + 246 * ccomps * dcomps);

            auto g_0_0_y_0_zz_yyyyz = cbuffer.data(dh_geom_0010_off + 247 * ccomps * dcomps);

            auto g_0_0_y_0_zz_yyyzz = cbuffer.data(dh_geom_0010_off + 248 * ccomps * dcomps);

            auto g_0_0_y_0_zz_yyzzz = cbuffer.data(dh_geom_0010_off + 249 * ccomps * dcomps);

            auto g_0_0_y_0_zz_yzzzz = cbuffer.data(dh_geom_0010_off + 250 * ccomps * dcomps);

            auto g_0_0_y_0_zz_zzzzz = cbuffer.data(dh_geom_0010_off + 251 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxxx = cbuffer.data(dh_geom_0010_off + 252 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxxy = cbuffer.data(dh_geom_0010_off + 253 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxxz = cbuffer.data(dh_geom_0010_off + 254 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxyy = cbuffer.data(dh_geom_0010_off + 255 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxyz = cbuffer.data(dh_geom_0010_off + 256 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxxzz = cbuffer.data(dh_geom_0010_off + 257 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxyyy = cbuffer.data(dh_geom_0010_off + 258 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxyyz = cbuffer.data(dh_geom_0010_off + 259 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxyzz = cbuffer.data(dh_geom_0010_off + 260 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xxzzz = cbuffer.data(dh_geom_0010_off + 261 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xyyyy = cbuffer.data(dh_geom_0010_off + 262 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xyyyz = cbuffer.data(dh_geom_0010_off + 263 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xyyzz = cbuffer.data(dh_geom_0010_off + 264 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xyzzz = cbuffer.data(dh_geom_0010_off + 265 * ccomps * dcomps);

            auto g_0_0_z_0_xx_xzzzz = cbuffer.data(dh_geom_0010_off + 266 * ccomps * dcomps);

            auto g_0_0_z_0_xx_yyyyy = cbuffer.data(dh_geom_0010_off + 267 * ccomps * dcomps);

            auto g_0_0_z_0_xx_yyyyz = cbuffer.data(dh_geom_0010_off + 268 * ccomps * dcomps);

            auto g_0_0_z_0_xx_yyyzz = cbuffer.data(dh_geom_0010_off + 269 * ccomps * dcomps);

            auto g_0_0_z_0_xx_yyzzz = cbuffer.data(dh_geom_0010_off + 270 * ccomps * dcomps);

            auto g_0_0_z_0_xx_yzzzz = cbuffer.data(dh_geom_0010_off + 271 * ccomps * dcomps);

            auto g_0_0_z_0_xx_zzzzz = cbuffer.data(dh_geom_0010_off + 272 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxxx = cbuffer.data(dh_geom_0010_off + 273 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxxy = cbuffer.data(dh_geom_0010_off + 274 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxxz = cbuffer.data(dh_geom_0010_off + 275 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxyy = cbuffer.data(dh_geom_0010_off + 276 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxyz = cbuffer.data(dh_geom_0010_off + 277 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxxzz = cbuffer.data(dh_geom_0010_off + 278 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxyyy = cbuffer.data(dh_geom_0010_off + 279 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxyyz = cbuffer.data(dh_geom_0010_off + 280 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxyzz = cbuffer.data(dh_geom_0010_off + 281 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xxzzz = cbuffer.data(dh_geom_0010_off + 282 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xyyyy = cbuffer.data(dh_geom_0010_off + 283 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xyyyz = cbuffer.data(dh_geom_0010_off + 284 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xyyzz = cbuffer.data(dh_geom_0010_off + 285 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xyzzz = cbuffer.data(dh_geom_0010_off + 286 * ccomps * dcomps);

            auto g_0_0_z_0_xy_xzzzz = cbuffer.data(dh_geom_0010_off + 287 * ccomps * dcomps);

            auto g_0_0_z_0_xy_yyyyy = cbuffer.data(dh_geom_0010_off + 288 * ccomps * dcomps);

            auto g_0_0_z_0_xy_yyyyz = cbuffer.data(dh_geom_0010_off + 289 * ccomps * dcomps);

            auto g_0_0_z_0_xy_yyyzz = cbuffer.data(dh_geom_0010_off + 290 * ccomps * dcomps);

            auto g_0_0_z_0_xy_yyzzz = cbuffer.data(dh_geom_0010_off + 291 * ccomps * dcomps);

            auto g_0_0_z_0_xy_yzzzz = cbuffer.data(dh_geom_0010_off + 292 * ccomps * dcomps);

            auto g_0_0_z_0_xy_zzzzz = cbuffer.data(dh_geom_0010_off + 293 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxxx = cbuffer.data(dh_geom_0010_off + 294 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxxy = cbuffer.data(dh_geom_0010_off + 295 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxxz = cbuffer.data(dh_geom_0010_off + 296 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxyy = cbuffer.data(dh_geom_0010_off + 297 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxyz = cbuffer.data(dh_geom_0010_off + 298 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxxzz = cbuffer.data(dh_geom_0010_off + 299 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxyyy = cbuffer.data(dh_geom_0010_off + 300 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxyyz = cbuffer.data(dh_geom_0010_off + 301 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxyzz = cbuffer.data(dh_geom_0010_off + 302 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xxzzz = cbuffer.data(dh_geom_0010_off + 303 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xyyyy = cbuffer.data(dh_geom_0010_off + 304 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xyyyz = cbuffer.data(dh_geom_0010_off + 305 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xyyzz = cbuffer.data(dh_geom_0010_off + 306 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xyzzz = cbuffer.data(dh_geom_0010_off + 307 * ccomps * dcomps);

            auto g_0_0_z_0_xz_xzzzz = cbuffer.data(dh_geom_0010_off + 308 * ccomps * dcomps);

            auto g_0_0_z_0_xz_yyyyy = cbuffer.data(dh_geom_0010_off + 309 * ccomps * dcomps);

            auto g_0_0_z_0_xz_yyyyz = cbuffer.data(dh_geom_0010_off + 310 * ccomps * dcomps);

            auto g_0_0_z_0_xz_yyyzz = cbuffer.data(dh_geom_0010_off + 311 * ccomps * dcomps);

            auto g_0_0_z_0_xz_yyzzz = cbuffer.data(dh_geom_0010_off + 312 * ccomps * dcomps);

            auto g_0_0_z_0_xz_yzzzz = cbuffer.data(dh_geom_0010_off + 313 * ccomps * dcomps);

            auto g_0_0_z_0_xz_zzzzz = cbuffer.data(dh_geom_0010_off + 314 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxxx = cbuffer.data(dh_geom_0010_off + 315 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxxy = cbuffer.data(dh_geom_0010_off + 316 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxxz = cbuffer.data(dh_geom_0010_off + 317 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxyy = cbuffer.data(dh_geom_0010_off + 318 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxyz = cbuffer.data(dh_geom_0010_off + 319 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxxzz = cbuffer.data(dh_geom_0010_off + 320 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxyyy = cbuffer.data(dh_geom_0010_off + 321 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxyyz = cbuffer.data(dh_geom_0010_off + 322 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxyzz = cbuffer.data(dh_geom_0010_off + 323 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xxzzz = cbuffer.data(dh_geom_0010_off + 324 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xyyyy = cbuffer.data(dh_geom_0010_off + 325 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xyyyz = cbuffer.data(dh_geom_0010_off + 326 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xyyzz = cbuffer.data(dh_geom_0010_off + 327 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xyzzz = cbuffer.data(dh_geom_0010_off + 328 * ccomps * dcomps);

            auto g_0_0_z_0_yy_xzzzz = cbuffer.data(dh_geom_0010_off + 329 * ccomps * dcomps);

            auto g_0_0_z_0_yy_yyyyy = cbuffer.data(dh_geom_0010_off + 330 * ccomps * dcomps);

            auto g_0_0_z_0_yy_yyyyz = cbuffer.data(dh_geom_0010_off + 331 * ccomps * dcomps);

            auto g_0_0_z_0_yy_yyyzz = cbuffer.data(dh_geom_0010_off + 332 * ccomps * dcomps);

            auto g_0_0_z_0_yy_yyzzz = cbuffer.data(dh_geom_0010_off + 333 * ccomps * dcomps);

            auto g_0_0_z_0_yy_yzzzz = cbuffer.data(dh_geom_0010_off + 334 * ccomps * dcomps);

            auto g_0_0_z_0_yy_zzzzz = cbuffer.data(dh_geom_0010_off + 335 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxxx = cbuffer.data(dh_geom_0010_off + 336 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxxy = cbuffer.data(dh_geom_0010_off + 337 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxxz = cbuffer.data(dh_geom_0010_off + 338 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxyy = cbuffer.data(dh_geom_0010_off + 339 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxyz = cbuffer.data(dh_geom_0010_off + 340 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxxzz = cbuffer.data(dh_geom_0010_off + 341 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxyyy = cbuffer.data(dh_geom_0010_off + 342 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxyyz = cbuffer.data(dh_geom_0010_off + 343 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxyzz = cbuffer.data(dh_geom_0010_off + 344 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xxzzz = cbuffer.data(dh_geom_0010_off + 345 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xyyyy = cbuffer.data(dh_geom_0010_off + 346 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xyyyz = cbuffer.data(dh_geom_0010_off + 347 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xyyzz = cbuffer.data(dh_geom_0010_off + 348 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xyzzz = cbuffer.data(dh_geom_0010_off + 349 * ccomps * dcomps);

            auto g_0_0_z_0_yz_xzzzz = cbuffer.data(dh_geom_0010_off + 350 * ccomps * dcomps);

            auto g_0_0_z_0_yz_yyyyy = cbuffer.data(dh_geom_0010_off + 351 * ccomps * dcomps);

            auto g_0_0_z_0_yz_yyyyz = cbuffer.data(dh_geom_0010_off + 352 * ccomps * dcomps);

            auto g_0_0_z_0_yz_yyyzz = cbuffer.data(dh_geom_0010_off + 353 * ccomps * dcomps);

            auto g_0_0_z_0_yz_yyzzz = cbuffer.data(dh_geom_0010_off + 354 * ccomps * dcomps);

            auto g_0_0_z_0_yz_yzzzz = cbuffer.data(dh_geom_0010_off + 355 * ccomps * dcomps);

            auto g_0_0_z_0_yz_zzzzz = cbuffer.data(dh_geom_0010_off + 356 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxxx = cbuffer.data(dh_geom_0010_off + 357 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxxy = cbuffer.data(dh_geom_0010_off + 358 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxxz = cbuffer.data(dh_geom_0010_off + 359 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxyy = cbuffer.data(dh_geom_0010_off + 360 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxyz = cbuffer.data(dh_geom_0010_off + 361 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxxzz = cbuffer.data(dh_geom_0010_off + 362 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxyyy = cbuffer.data(dh_geom_0010_off + 363 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxyyz = cbuffer.data(dh_geom_0010_off + 364 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxyzz = cbuffer.data(dh_geom_0010_off + 365 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xxzzz = cbuffer.data(dh_geom_0010_off + 366 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xyyyy = cbuffer.data(dh_geom_0010_off + 367 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xyyyz = cbuffer.data(dh_geom_0010_off + 368 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xyyzz = cbuffer.data(dh_geom_0010_off + 369 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xyzzz = cbuffer.data(dh_geom_0010_off + 370 * ccomps * dcomps);

            auto g_0_0_z_0_zz_xzzzz = cbuffer.data(dh_geom_0010_off + 371 * ccomps * dcomps);

            auto g_0_0_z_0_zz_yyyyy = cbuffer.data(dh_geom_0010_off + 372 * ccomps * dcomps);

            auto g_0_0_z_0_zz_yyyyz = cbuffer.data(dh_geom_0010_off + 373 * ccomps * dcomps);

            auto g_0_0_z_0_zz_yyyzz = cbuffer.data(dh_geom_0010_off + 374 * ccomps * dcomps);

            auto g_0_0_z_0_zz_yyzzz = cbuffer.data(dh_geom_0010_off + 375 * ccomps * dcomps);

            auto g_0_0_z_0_zz_yzzzz = cbuffer.data(dh_geom_0010_off + 376 * ccomps * dcomps);

            auto g_0_0_z_0_zz_zzzzz = cbuffer.data(dh_geom_0010_off + 377 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_1010_off = idx_geom_1010_dhxx + i * dcomps + j;

            auto g_x_0_x_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_y_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_y_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_y_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_y_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_z_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_z_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_z_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_z_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_z_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_z_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 377 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 378 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 379 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 380 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 381 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 382 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 383 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 384 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 385 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 386 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 387 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 388 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 389 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 390 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 391 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 392 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 393 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 394 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 395 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 396 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 397 * ccomps * dcomps);

            auto g_y_0_x_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 398 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 399 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 400 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 401 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 402 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 403 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 404 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 405 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 406 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 407 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 408 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 409 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 410 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 411 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 412 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 413 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 414 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 415 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 416 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 417 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 418 * ccomps * dcomps);

            auto g_y_0_x_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 419 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 420 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 421 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 422 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 423 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 424 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 425 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 426 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 427 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 428 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 429 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 430 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 431 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 432 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 433 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 434 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 435 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 436 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 437 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 438 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 439 * ccomps * dcomps);

            auto g_y_0_x_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 440 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 441 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 442 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 443 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 444 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 445 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 446 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 447 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 448 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 449 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 450 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 451 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 452 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 453 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 454 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 455 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 456 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 457 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 458 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 459 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 460 * ccomps * dcomps);

            auto g_y_0_x_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 461 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 462 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 463 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 464 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 465 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 466 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 467 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 468 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 469 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 470 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 471 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 472 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 473 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 474 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 475 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 476 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 477 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 478 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 479 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 480 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 481 * ccomps * dcomps);

            auto g_y_0_x_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 482 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 483 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 484 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 485 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 486 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 487 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 488 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 489 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 490 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 491 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 492 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 493 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 494 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 495 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 496 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 497 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 498 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 499 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 500 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 501 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 502 * ccomps * dcomps);

            auto g_y_0_x_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 503 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 504 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 505 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 506 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 507 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 508 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 509 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 510 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 511 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 512 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 513 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 514 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 515 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 516 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 517 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 518 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 519 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 520 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 521 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 522 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 523 * ccomps * dcomps);

            auto g_y_0_y_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 524 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 525 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 526 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 527 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 528 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 529 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 530 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 531 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 532 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 533 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 534 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 535 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 536 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 537 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 538 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 539 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 540 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 541 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 542 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 543 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 544 * ccomps * dcomps);

            auto g_y_0_y_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 545 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 546 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 547 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 548 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 549 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 550 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 551 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 552 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 553 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 554 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 555 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 556 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 557 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 558 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 559 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 560 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 561 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 562 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 563 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 564 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 565 * ccomps * dcomps);

            auto g_y_0_y_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 566 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 567 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 568 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 569 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 570 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 571 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 572 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 573 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 574 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 575 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 576 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 577 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 578 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 579 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 580 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 581 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 582 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 583 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 584 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 585 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 586 * ccomps * dcomps);

            auto g_y_0_y_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 587 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 588 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 589 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 590 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 591 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 592 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 593 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 594 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 595 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 596 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 597 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 598 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 599 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 600 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 601 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 602 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 603 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 604 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 605 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 606 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 607 * ccomps * dcomps);

            auto g_y_0_y_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 608 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 609 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 610 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 611 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 612 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 613 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 614 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 615 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 616 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 617 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 618 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 619 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 620 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 621 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 622 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 623 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 624 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 625 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 626 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 627 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 628 * ccomps * dcomps);

            auto g_y_0_y_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 629 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_z_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_z_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 671 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_z_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_z_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_z_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_z_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 755 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 756 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 757 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 758 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 759 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 760 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 761 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 762 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 763 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 764 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 765 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 766 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 767 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 768 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 769 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 770 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 771 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 772 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 773 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 774 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 775 * ccomps * dcomps);

            auto g_z_0_x_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 776 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 777 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 778 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 779 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 780 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 781 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 782 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 783 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 784 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 785 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 786 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 787 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 788 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 789 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 790 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 791 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 792 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 793 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 794 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 795 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 796 * ccomps * dcomps);

            auto g_z_0_x_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 797 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 798 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 799 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 800 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 801 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 802 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 803 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 804 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 805 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 806 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 807 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 808 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 809 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 810 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 811 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 812 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 813 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 814 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 815 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 816 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 817 * ccomps * dcomps);

            auto g_z_0_x_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 818 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 819 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 820 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 821 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 822 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 823 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 824 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 825 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 826 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 827 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 828 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 829 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 830 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 831 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 832 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 833 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 834 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 835 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 836 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 837 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 838 * ccomps * dcomps);

            auto g_z_0_x_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 839 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 840 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 841 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 842 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 843 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 844 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 845 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 846 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 847 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 848 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 849 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 850 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 851 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 852 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 853 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 854 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 855 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 856 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 857 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 858 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 859 * ccomps * dcomps);

            auto g_z_0_x_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 860 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 861 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 862 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 863 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 864 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 865 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 866 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 867 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 868 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 869 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 870 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 871 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 872 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 873 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 874 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 875 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 876 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 877 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 878 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 879 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 880 * ccomps * dcomps);

            auto g_z_0_x_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 881 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 882 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 883 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 884 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 885 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 886 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 887 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 888 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 889 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 890 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 891 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 892 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 893 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 894 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 895 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 896 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 897 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 898 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 899 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 900 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 901 * ccomps * dcomps);

            auto g_z_0_y_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 902 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 903 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 904 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 905 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 906 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 907 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 908 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 909 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 910 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 911 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 912 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 913 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 914 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 915 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 916 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 917 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 918 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 919 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 920 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 921 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 922 * ccomps * dcomps);

            auto g_z_0_y_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 923 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 924 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 925 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 926 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 927 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 928 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 929 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 930 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 931 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 932 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 933 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 934 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 935 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 936 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 937 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 938 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 939 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 940 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 941 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 942 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 943 * ccomps * dcomps);

            auto g_z_0_y_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 944 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 945 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 946 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 947 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 948 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 949 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 950 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 951 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 952 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 953 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 954 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 955 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 956 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 957 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 958 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 959 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 960 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 961 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 962 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 963 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 964 * ccomps * dcomps);

            auto g_z_0_y_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 965 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 966 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 967 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 968 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 969 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 970 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 971 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 972 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 973 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 974 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 975 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 976 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 977 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 978 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 979 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 980 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 981 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 982 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 983 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 984 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 985 * ccomps * dcomps);

            auto g_z_0_y_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 986 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 987 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 988 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 989 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 990 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 991 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 992 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 993 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 994 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 995 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 996 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 997 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 998 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 999 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_z_0_y_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxx = cbuffer.data(dh_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxy = cbuffer.data(dh_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxz = cbuffer.data(dh_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxyy = cbuffer.data(dh_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxyz = cbuffer.data(dh_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxzz = cbuffer.data(dh_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyyy = cbuffer.data(dh_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyyz = cbuffer.data(dh_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyzz = cbuffer.data(dh_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxzzz = cbuffer.data(dh_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyyy = cbuffer.data(dh_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyyz = cbuffer.data(dh_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyzz = cbuffer.data(dh_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyzzz = cbuffer.data(dh_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xzzzz = cbuffer.data(dh_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyyy = cbuffer.data(dh_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyyz = cbuffer.data(dh_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyzz = cbuffer.data(dh_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyzzz = cbuffer.data(dh_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yzzzz = cbuffer.data(dh_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_z_0_z_0_xx_zzzzz = cbuffer.data(dh_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxx = cbuffer.data(dh_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxy = cbuffer.data(dh_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxz = cbuffer.data(dh_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxyy = cbuffer.data(dh_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxyz = cbuffer.data(dh_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxzz = cbuffer.data(dh_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyyy = cbuffer.data(dh_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyyz = cbuffer.data(dh_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyzz = cbuffer.data(dh_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxzzz = cbuffer.data(dh_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyyy = cbuffer.data(dh_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyyz = cbuffer.data(dh_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyzz = cbuffer.data(dh_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyzzz = cbuffer.data(dh_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xzzzz = cbuffer.data(dh_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyyy = cbuffer.data(dh_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyyz = cbuffer.data(dh_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyzz = cbuffer.data(dh_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyzzz = cbuffer.data(dh_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yzzzz = cbuffer.data(dh_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_z_0_z_0_xy_zzzzz = cbuffer.data(dh_geom_1010_off + 1049 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxx = cbuffer.data(dh_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxy = cbuffer.data(dh_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxz = cbuffer.data(dh_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxyy = cbuffer.data(dh_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxyz = cbuffer.data(dh_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxzz = cbuffer.data(dh_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyyy = cbuffer.data(dh_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyyz = cbuffer.data(dh_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyzz = cbuffer.data(dh_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxzzz = cbuffer.data(dh_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyyy = cbuffer.data(dh_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyyz = cbuffer.data(dh_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyzz = cbuffer.data(dh_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyzzz = cbuffer.data(dh_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xzzzz = cbuffer.data(dh_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyyy = cbuffer.data(dh_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyyz = cbuffer.data(dh_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyzz = cbuffer.data(dh_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyzzz = cbuffer.data(dh_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yzzzz = cbuffer.data(dh_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_z_0_z_0_xz_zzzzz = cbuffer.data(dh_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxx = cbuffer.data(dh_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxy = cbuffer.data(dh_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxz = cbuffer.data(dh_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxyy = cbuffer.data(dh_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxyz = cbuffer.data(dh_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxzz = cbuffer.data(dh_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyyy = cbuffer.data(dh_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyyz = cbuffer.data(dh_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyzz = cbuffer.data(dh_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxzzz = cbuffer.data(dh_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyyy = cbuffer.data(dh_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyyz = cbuffer.data(dh_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyzz = cbuffer.data(dh_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyzzz = cbuffer.data(dh_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xzzzz = cbuffer.data(dh_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyyy = cbuffer.data(dh_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyyz = cbuffer.data(dh_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyzz = cbuffer.data(dh_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyzzz = cbuffer.data(dh_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yzzzz = cbuffer.data(dh_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_z_0_z_0_yy_zzzzz = cbuffer.data(dh_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxx = cbuffer.data(dh_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxy = cbuffer.data(dh_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxz = cbuffer.data(dh_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxyy = cbuffer.data(dh_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxyz = cbuffer.data(dh_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxzz = cbuffer.data(dh_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyyy = cbuffer.data(dh_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyyz = cbuffer.data(dh_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyzz = cbuffer.data(dh_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxzzz = cbuffer.data(dh_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyyy = cbuffer.data(dh_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyyz = cbuffer.data(dh_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyzz = cbuffer.data(dh_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyzzz = cbuffer.data(dh_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xzzzz = cbuffer.data(dh_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyyy = cbuffer.data(dh_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyyz = cbuffer.data(dh_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyzz = cbuffer.data(dh_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyzzz = cbuffer.data(dh_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yzzzz = cbuffer.data(dh_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_z_0_z_0_yz_zzzzz = cbuffer.data(dh_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxx = cbuffer.data(dh_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxy = cbuffer.data(dh_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxz = cbuffer.data(dh_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxyy = cbuffer.data(dh_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxyz = cbuffer.data(dh_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxzz = cbuffer.data(dh_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyyy = cbuffer.data(dh_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyyz = cbuffer.data(dh_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyzz = cbuffer.data(dh_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxzzz = cbuffer.data(dh_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyyy = cbuffer.data(dh_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyyz = cbuffer.data(dh_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyzz = cbuffer.data(dh_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyzzz = cbuffer.data(dh_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xzzzz = cbuffer.data(dh_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyyy = cbuffer.data(dh_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyyz = cbuffer.data(dh_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyzz = cbuffer.data(dh_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyzzz = cbuffer.data(dh_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yzzzz = cbuffer.data(dh_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_z_0_z_0_zz_zzzzz = cbuffer.data(dh_geom_1010_off + 1133 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DISS

            const auto di_geom_1010_off = idx_geom_1010_dixx + i * dcomps + j;

            auto g_x_0_x_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_x_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_x_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_x_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_x_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_x_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_y_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_y_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_y_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 251 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 272 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 293 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_y_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_y_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_y_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 314 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_y_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_y_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_y_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 335 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_z_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 356 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_z_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_z_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 377 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 383 * ccomps * dcomps);

            auto g_x_0_z_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 389 * ccomps * dcomps);

            auto g_x_0_z_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_z_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 395 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 398 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 401 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 404 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 407 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_z_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 413 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_z_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_z_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 419 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 425 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 431 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 434 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 437 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_z_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 440 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 443 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_z_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_z_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 449 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 450 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 451 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 452 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 453 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 454 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 455 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 456 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 457 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 458 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 459 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 460 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 461 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 462 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 463 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 464 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 465 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 466 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 467 * ccomps * dcomps);

            auto g_x_0_z_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 468 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 469 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 470 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 471 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 472 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 473 * ccomps * dcomps);

            auto g_x_0_z_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 474 * ccomps * dcomps);

            auto g_x_0_z_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 475 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 476 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 477 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 478 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 479 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 480 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 481 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 482 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 483 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 484 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 485 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 486 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 487 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 488 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 489 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 490 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 491 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 492 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 493 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 494 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 495 * ccomps * dcomps);

            auto g_x_0_z_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 496 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 497 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 498 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 499 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 500 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 501 * ccomps * dcomps);

            auto g_x_0_z_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 502 * ccomps * dcomps);

            auto g_x_0_z_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 503 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 504 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 505 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 506 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 507 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 508 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 509 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 510 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 511 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 512 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 513 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 514 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 515 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 516 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 517 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 518 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 519 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 520 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 521 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 522 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 523 * ccomps * dcomps);

            auto g_y_0_x_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 524 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 525 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 526 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 527 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 528 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 529 * ccomps * dcomps);

            auto g_y_0_x_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 530 * ccomps * dcomps);

            auto g_y_0_x_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 531 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 532 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 533 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 534 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 535 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 536 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 537 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 538 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 539 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 540 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 541 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 542 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 543 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 544 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 545 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 546 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 547 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 548 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 549 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 550 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 551 * ccomps * dcomps);

            auto g_y_0_x_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 552 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 553 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 554 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 555 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 556 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 557 * ccomps * dcomps);

            auto g_y_0_x_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 558 * ccomps * dcomps);

            auto g_y_0_x_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 559 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 560 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 561 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 562 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 563 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 564 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 565 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 566 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 567 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 568 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 569 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 570 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 571 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 572 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 573 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 574 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 575 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 576 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 577 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 578 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 579 * ccomps * dcomps);

            auto g_y_0_x_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 580 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 581 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 582 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 583 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 584 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 585 * ccomps * dcomps);

            auto g_y_0_x_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 586 * ccomps * dcomps);

            auto g_y_0_x_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 587 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 588 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 589 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 590 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 591 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 592 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 593 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 594 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 595 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 596 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 597 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 598 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 599 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 600 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 601 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 602 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 603 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 604 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 605 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 606 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 607 * ccomps * dcomps);

            auto g_y_0_x_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 608 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 609 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 610 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 611 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 612 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 613 * ccomps * dcomps);

            auto g_y_0_x_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 614 * ccomps * dcomps);

            auto g_y_0_x_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 615 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 616 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 617 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 618 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 619 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 620 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 621 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 622 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 623 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 624 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 625 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 626 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 627 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 628 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 629 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_x_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_x_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_x_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 650 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_x_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_x_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_x_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 671 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_y_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 692 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_y_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_y_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 713 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_y_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_y_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_y_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 734 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_y_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_y_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_y_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 755 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 761 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 764 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 767 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 773 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_y_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 776 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 779 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_y_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_y_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 785 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 791 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 794 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 797 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 803 * ccomps * dcomps);

            auto g_y_0_y_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 809 * ccomps * dcomps);

            auto g_y_0_y_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_y_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 815 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 818 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 821 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 824 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 827 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_y_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 833 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_y_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_y_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 839 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 845 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 851 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 854 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 857 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_z_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 860 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 863 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_z_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_z_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 869 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 875 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 881 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 884 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 887 * ccomps * dcomps);

            auto g_y_0_z_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 893 * ccomps * dcomps);

            auto g_y_0_z_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_z_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 899 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 900 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 901 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 902 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 903 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 904 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 905 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 906 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 907 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 908 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 909 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 910 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 911 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 912 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 913 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 914 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 915 * ccomps * dcomps);

            auto g_y_0_z_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 916 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 917 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 918 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 919 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 920 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 921 * ccomps * dcomps);

            auto g_y_0_z_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 922 * ccomps * dcomps);

            auto g_y_0_z_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 923 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 924 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 925 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 926 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 927 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 928 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 929 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 930 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 931 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 932 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 933 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 934 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 935 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 936 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 937 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 938 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 939 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 940 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 941 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 942 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 943 * ccomps * dcomps);

            auto g_y_0_z_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 944 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 945 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 946 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 947 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 948 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 949 * ccomps * dcomps);

            auto g_y_0_z_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 950 * ccomps * dcomps);

            auto g_y_0_z_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 951 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 952 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 953 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 954 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 955 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 956 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 957 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 958 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 959 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 960 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 961 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 962 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 963 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 964 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 965 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 966 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 967 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 968 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 969 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 970 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 971 * ccomps * dcomps);

            auto g_y_0_z_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 972 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 973 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 974 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 975 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 976 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 977 * ccomps * dcomps);

            auto g_y_0_z_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 978 * ccomps * dcomps);

            auto g_y_0_z_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 979 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 980 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 981 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 982 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 983 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 984 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 985 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 986 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 987 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 988 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 989 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 990 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 991 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 992 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 993 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 994 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 995 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 996 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 997 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 998 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 999 * ccomps * dcomps);

            auto g_y_0_z_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_y_0_z_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_y_0_z_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 1007 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_z_0_x_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 1028 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_z_0_x_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_z_0_x_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 1049 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_z_0_x_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_z_0_x_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_z_0_x_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 1070 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_z_0_x_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_z_0_x_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_z_0_x_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 1091 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_z_0_x_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 1112 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_z_0_x_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_z_0_x_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 1133 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 1139 * ccomps * dcomps);

            auto g_z_0_x_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 1145 * ccomps * dcomps);

            auto g_z_0_x_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_z_0_x_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 1151 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 1154 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 1157 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 1163 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_z_0_x_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 1169 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_z_0_x_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_z_0_x_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 1175 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 1181 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 1184 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 1187 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 1193 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_z_0_y_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 1196 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 1199 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_z_0_y_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_z_0_y_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 1205 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 1211 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 1214 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 1217 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 1223 * ccomps * dcomps);

            auto g_z_0_y_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 1229 * ccomps * dcomps);

            auto g_z_0_y_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_z_0_y_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 1235 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 1238 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 1241 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 1244 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 1247 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_z_0_y_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 1253 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_z_0_y_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_z_0_y_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 1259 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 1265 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 1271 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 1274 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 1277 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_z_0_y_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 1280 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 1283 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_z_0_y_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_z_0_y_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 1289 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 1295 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 1301 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 1304 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 1307 * ccomps * dcomps);

            auto g_z_0_y_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 1313 * ccomps * dcomps);

            auto g_z_0_y_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_z_0_y_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 1319 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 1322 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 1325 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 1331 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 1334 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_z_0_y_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 1337 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_z_0_y_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_z_0_y_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 1343 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxxx = cbuffer.data(di_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxxy = cbuffer.data(di_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxxz = cbuffer.data(di_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxyy = cbuffer.data(di_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxyz = cbuffer.data(di_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxxzz = cbuffer.data(di_geom_1010_off + 1349 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxyyy = cbuffer.data(di_geom_1010_off + 1350 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxyyz = cbuffer.data(di_geom_1010_off + 1351 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxyzz = cbuffer.data(di_geom_1010_off + 1352 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxxzzz = cbuffer.data(di_geom_1010_off + 1353 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyyyy = cbuffer.data(di_geom_1010_off + 1354 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyyyz = cbuffer.data(di_geom_1010_off + 1355 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyyzz = cbuffer.data(di_geom_1010_off + 1356 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxyzzz = cbuffer.data(di_geom_1010_off + 1357 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xxzzzz = cbuffer.data(di_geom_1010_off + 1358 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyyyy = cbuffer.data(di_geom_1010_off + 1359 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyyyz = cbuffer.data(di_geom_1010_off + 1360 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyyzz = cbuffer.data(di_geom_1010_off + 1361 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyyzzz = cbuffer.data(di_geom_1010_off + 1362 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xyzzzz = cbuffer.data(di_geom_1010_off + 1363 * ccomps * dcomps);

            auto g_z_0_z_0_xx_xzzzzz = cbuffer.data(di_geom_1010_off + 1364 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyyyy = cbuffer.data(di_geom_1010_off + 1365 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyyyz = cbuffer.data(di_geom_1010_off + 1366 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyyzz = cbuffer.data(di_geom_1010_off + 1367 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyyzzz = cbuffer.data(di_geom_1010_off + 1368 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yyzzzz = cbuffer.data(di_geom_1010_off + 1369 * ccomps * dcomps);

            auto g_z_0_z_0_xx_yzzzzz = cbuffer.data(di_geom_1010_off + 1370 * ccomps * dcomps);

            auto g_z_0_z_0_xx_zzzzzz = cbuffer.data(di_geom_1010_off + 1371 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxxx = cbuffer.data(di_geom_1010_off + 1372 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxxy = cbuffer.data(di_geom_1010_off + 1373 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxxz = cbuffer.data(di_geom_1010_off + 1374 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxyy = cbuffer.data(di_geom_1010_off + 1375 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxyz = cbuffer.data(di_geom_1010_off + 1376 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxxzz = cbuffer.data(di_geom_1010_off + 1377 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxyyy = cbuffer.data(di_geom_1010_off + 1378 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxyyz = cbuffer.data(di_geom_1010_off + 1379 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxyzz = cbuffer.data(di_geom_1010_off + 1380 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxxzzz = cbuffer.data(di_geom_1010_off + 1381 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyyyy = cbuffer.data(di_geom_1010_off + 1382 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyyyz = cbuffer.data(di_geom_1010_off + 1383 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyyzz = cbuffer.data(di_geom_1010_off + 1384 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxyzzz = cbuffer.data(di_geom_1010_off + 1385 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xxzzzz = cbuffer.data(di_geom_1010_off + 1386 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyyyy = cbuffer.data(di_geom_1010_off + 1387 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyyyz = cbuffer.data(di_geom_1010_off + 1388 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyyzz = cbuffer.data(di_geom_1010_off + 1389 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyyzzz = cbuffer.data(di_geom_1010_off + 1390 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xyzzzz = cbuffer.data(di_geom_1010_off + 1391 * ccomps * dcomps);

            auto g_z_0_z_0_xy_xzzzzz = cbuffer.data(di_geom_1010_off + 1392 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyyyy = cbuffer.data(di_geom_1010_off + 1393 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyyyz = cbuffer.data(di_geom_1010_off + 1394 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyyzz = cbuffer.data(di_geom_1010_off + 1395 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyyzzz = cbuffer.data(di_geom_1010_off + 1396 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yyzzzz = cbuffer.data(di_geom_1010_off + 1397 * ccomps * dcomps);

            auto g_z_0_z_0_xy_yzzzzz = cbuffer.data(di_geom_1010_off + 1398 * ccomps * dcomps);

            auto g_z_0_z_0_xy_zzzzzz = cbuffer.data(di_geom_1010_off + 1399 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxxx = cbuffer.data(di_geom_1010_off + 1400 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxxy = cbuffer.data(di_geom_1010_off + 1401 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxxz = cbuffer.data(di_geom_1010_off + 1402 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxyy = cbuffer.data(di_geom_1010_off + 1403 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxyz = cbuffer.data(di_geom_1010_off + 1404 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxxzz = cbuffer.data(di_geom_1010_off + 1405 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxyyy = cbuffer.data(di_geom_1010_off + 1406 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxyyz = cbuffer.data(di_geom_1010_off + 1407 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxyzz = cbuffer.data(di_geom_1010_off + 1408 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxxzzz = cbuffer.data(di_geom_1010_off + 1409 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyyyy = cbuffer.data(di_geom_1010_off + 1410 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyyyz = cbuffer.data(di_geom_1010_off + 1411 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyyzz = cbuffer.data(di_geom_1010_off + 1412 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxyzzz = cbuffer.data(di_geom_1010_off + 1413 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xxzzzz = cbuffer.data(di_geom_1010_off + 1414 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyyyy = cbuffer.data(di_geom_1010_off + 1415 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyyyz = cbuffer.data(di_geom_1010_off + 1416 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyyzz = cbuffer.data(di_geom_1010_off + 1417 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyyzzz = cbuffer.data(di_geom_1010_off + 1418 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xyzzzz = cbuffer.data(di_geom_1010_off + 1419 * ccomps * dcomps);

            auto g_z_0_z_0_xz_xzzzzz = cbuffer.data(di_geom_1010_off + 1420 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyyyy = cbuffer.data(di_geom_1010_off + 1421 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyyyz = cbuffer.data(di_geom_1010_off + 1422 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyyzz = cbuffer.data(di_geom_1010_off + 1423 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyyzzz = cbuffer.data(di_geom_1010_off + 1424 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yyzzzz = cbuffer.data(di_geom_1010_off + 1425 * ccomps * dcomps);

            auto g_z_0_z_0_xz_yzzzzz = cbuffer.data(di_geom_1010_off + 1426 * ccomps * dcomps);

            auto g_z_0_z_0_xz_zzzzzz = cbuffer.data(di_geom_1010_off + 1427 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxxx = cbuffer.data(di_geom_1010_off + 1428 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxxy = cbuffer.data(di_geom_1010_off + 1429 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxxz = cbuffer.data(di_geom_1010_off + 1430 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxyy = cbuffer.data(di_geom_1010_off + 1431 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxyz = cbuffer.data(di_geom_1010_off + 1432 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxxzz = cbuffer.data(di_geom_1010_off + 1433 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxyyy = cbuffer.data(di_geom_1010_off + 1434 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxyyz = cbuffer.data(di_geom_1010_off + 1435 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxyzz = cbuffer.data(di_geom_1010_off + 1436 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxxzzz = cbuffer.data(di_geom_1010_off + 1437 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyyyy = cbuffer.data(di_geom_1010_off + 1438 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyyyz = cbuffer.data(di_geom_1010_off + 1439 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyyzz = cbuffer.data(di_geom_1010_off + 1440 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxyzzz = cbuffer.data(di_geom_1010_off + 1441 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xxzzzz = cbuffer.data(di_geom_1010_off + 1442 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyyyy = cbuffer.data(di_geom_1010_off + 1443 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyyyz = cbuffer.data(di_geom_1010_off + 1444 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyyzz = cbuffer.data(di_geom_1010_off + 1445 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyyzzz = cbuffer.data(di_geom_1010_off + 1446 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xyzzzz = cbuffer.data(di_geom_1010_off + 1447 * ccomps * dcomps);

            auto g_z_0_z_0_yy_xzzzzz = cbuffer.data(di_geom_1010_off + 1448 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyyyy = cbuffer.data(di_geom_1010_off + 1449 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyyyz = cbuffer.data(di_geom_1010_off + 1450 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyyzz = cbuffer.data(di_geom_1010_off + 1451 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyyzzz = cbuffer.data(di_geom_1010_off + 1452 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yyzzzz = cbuffer.data(di_geom_1010_off + 1453 * ccomps * dcomps);

            auto g_z_0_z_0_yy_yzzzzz = cbuffer.data(di_geom_1010_off + 1454 * ccomps * dcomps);

            auto g_z_0_z_0_yy_zzzzzz = cbuffer.data(di_geom_1010_off + 1455 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxxx = cbuffer.data(di_geom_1010_off + 1456 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxxy = cbuffer.data(di_geom_1010_off + 1457 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxxz = cbuffer.data(di_geom_1010_off + 1458 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxyy = cbuffer.data(di_geom_1010_off + 1459 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxyz = cbuffer.data(di_geom_1010_off + 1460 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxxzz = cbuffer.data(di_geom_1010_off + 1461 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxyyy = cbuffer.data(di_geom_1010_off + 1462 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxyyz = cbuffer.data(di_geom_1010_off + 1463 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxyzz = cbuffer.data(di_geom_1010_off + 1464 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxxzzz = cbuffer.data(di_geom_1010_off + 1465 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyyyy = cbuffer.data(di_geom_1010_off + 1466 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyyyz = cbuffer.data(di_geom_1010_off + 1467 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyyzz = cbuffer.data(di_geom_1010_off + 1468 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxyzzz = cbuffer.data(di_geom_1010_off + 1469 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xxzzzz = cbuffer.data(di_geom_1010_off + 1470 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyyyy = cbuffer.data(di_geom_1010_off + 1471 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyyyz = cbuffer.data(di_geom_1010_off + 1472 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyyzz = cbuffer.data(di_geom_1010_off + 1473 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyyzzz = cbuffer.data(di_geom_1010_off + 1474 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xyzzzz = cbuffer.data(di_geom_1010_off + 1475 * ccomps * dcomps);

            auto g_z_0_z_0_yz_xzzzzz = cbuffer.data(di_geom_1010_off + 1476 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyyyy = cbuffer.data(di_geom_1010_off + 1477 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyyyz = cbuffer.data(di_geom_1010_off + 1478 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyyzz = cbuffer.data(di_geom_1010_off + 1479 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyyzzz = cbuffer.data(di_geom_1010_off + 1480 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yyzzzz = cbuffer.data(di_geom_1010_off + 1481 * ccomps * dcomps);

            auto g_z_0_z_0_yz_yzzzzz = cbuffer.data(di_geom_1010_off + 1482 * ccomps * dcomps);

            auto g_z_0_z_0_yz_zzzzzz = cbuffer.data(di_geom_1010_off + 1483 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxxx = cbuffer.data(di_geom_1010_off + 1484 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxxy = cbuffer.data(di_geom_1010_off + 1485 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxxz = cbuffer.data(di_geom_1010_off + 1486 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxyy = cbuffer.data(di_geom_1010_off + 1487 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxyz = cbuffer.data(di_geom_1010_off + 1488 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxxzz = cbuffer.data(di_geom_1010_off + 1489 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxyyy = cbuffer.data(di_geom_1010_off + 1490 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxyyz = cbuffer.data(di_geom_1010_off + 1491 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxyzz = cbuffer.data(di_geom_1010_off + 1492 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxxzzz = cbuffer.data(di_geom_1010_off + 1493 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyyyy = cbuffer.data(di_geom_1010_off + 1494 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyyyz = cbuffer.data(di_geom_1010_off + 1495 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyyzz = cbuffer.data(di_geom_1010_off + 1496 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxyzzz = cbuffer.data(di_geom_1010_off + 1497 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xxzzzz = cbuffer.data(di_geom_1010_off + 1498 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyyyy = cbuffer.data(di_geom_1010_off + 1499 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyyyz = cbuffer.data(di_geom_1010_off + 1500 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyyzz = cbuffer.data(di_geom_1010_off + 1501 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyyzzz = cbuffer.data(di_geom_1010_off + 1502 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xyzzzz = cbuffer.data(di_geom_1010_off + 1503 * ccomps * dcomps);

            auto g_z_0_z_0_zz_xzzzzz = cbuffer.data(di_geom_1010_off + 1504 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyyyy = cbuffer.data(di_geom_1010_off + 1505 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyyyz = cbuffer.data(di_geom_1010_off + 1506 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyyzz = cbuffer.data(di_geom_1010_off + 1507 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyyzzz = cbuffer.data(di_geom_1010_off + 1508 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yyzzzz = cbuffer.data(di_geom_1010_off + 1509 * ccomps * dcomps);

            auto g_z_0_z_0_zz_yzzzzz = cbuffer.data(di_geom_1010_off + 1510 * ccomps * dcomps);

            auto g_z_0_z_0_zz_zzzzzz = cbuffer.data(di_geom_1010_off + 1511 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fhxx

            const auto fh_geom_1010_off = idx_geom_1010_fhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xx_xxxxx, g_0_0_x_0_xx_xxxxy, g_0_0_x_0_xx_xxxxz, g_0_0_x_0_xx_xxxyy, g_0_0_x_0_xx_xxxyz, g_0_0_x_0_xx_xxxzz, g_0_0_x_0_xx_xxyyy, g_0_0_x_0_xx_xxyyz, g_0_0_x_0_xx_xxyzz, g_0_0_x_0_xx_xxzzz, g_0_0_x_0_xx_xyyyy, g_0_0_x_0_xx_xyyyz, g_0_0_x_0_xx_xyyzz, g_0_0_x_0_xx_xyzzz, g_0_0_x_0_xx_xzzzz, g_0_0_x_0_xx_yyyyy, g_0_0_x_0_xx_yyyyz, g_0_0_x_0_xx_yyyzz, g_0_0_x_0_xx_yyzzz, g_0_0_x_0_xx_yzzzz, g_0_0_x_0_xx_zzzzz, g_x_0_x_0_xx_xxxxx, g_x_0_x_0_xx_xxxxxx, g_x_0_x_0_xx_xxxxxy, g_x_0_x_0_xx_xxxxxz, g_x_0_x_0_xx_xxxxy, g_x_0_x_0_xx_xxxxyy, g_x_0_x_0_xx_xxxxyz, g_x_0_x_0_xx_xxxxz, g_x_0_x_0_xx_xxxxzz, g_x_0_x_0_xx_xxxyy, g_x_0_x_0_xx_xxxyyy, g_x_0_x_0_xx_xxxyyz, g_x_0_x_0_xx_xxxyz, g_x_0_x_0_xx_xxxyzz, g_x_0_x_0_xx_xxxzz, g_x_0_x_0_xx_xxxzzz, g_x_0_x_0_xx_xxyyy, g_x_0_x_0_xx_xxyyyy, g_x_0_x_0_xx_xxyyyz, g_x_0_x_0_xx_xxyyz, g_x_0_x_0_xx_xxyyzz, g_x_0_x_0_xx_xxyzz, g_x_0_x_0_xx_xxyzzz, g_x_0_x_0_xx_xxzzz, g_x_0_x_0_xx_xxzzzz, g_x_0_x_0_xx_xyyyy, g_x_0_x_0_xx_xyyyyy, g_x_0_x_0_xx_xyyyyz, g_x_0_x_0_xx_xyyyz, g_x_0_x_0_xx_xyyyzz, g_x_0_x_0_xx_xyyzz, g_x_0_x_0_xx_xyyzzz, g_x_0_x_0_xx_xyzzz, g_x_0_x_0_xx_xyzzzz, g_x_0_x_0_xx_xzzzz, g_x_0_x_0_xx_xzzzzz, g_x_0_x_0_xx_yyyyy, g_x_0_x_0_xx_yyyyz, g_x_0_x_0_xx_yyyzz, g_x_0_x_0_xx_yyzzz, g_x_0_x_0_xx_yzzzz, g_x_0_x_0_xx_zzzzz, g_x_0_x_0_xxx_xxxxx, g_x_0_x_0_xxx_xxxxy, g_x_0_x_0_xxx_xxxxz, g_x_0_x_0_xxx_xxxyy, g_x_0_x_0_xxx_xxxyz, g_x_0_x_0_xxx_xxxzz, g_x_0_x_0_xxx_xxyyy, g_x_0_x_0_xxx_xxyyz, g_x_0_x_0_xxx_xxyzz, g_x_0_x_0_xxx_xxzzz, g_x_0_x_0_xxx_xyyyy, g_x_0_x_0_xxx_xyyyz, g_x_0_x_0_xxx_xyyzz, g_x_0_x_0_xxx_xyzzz, g_x_0_x_0_xxx_xzzzz, g_x_0_x_0_xxx_yyyyy, g_x_0_x_0_xxx_yyyyz, g_x_0_x_0_xxx_yyyzz, g_x_0_x_0_xxx_yyzzz, g_x_0_x_0_xxx_yzzzz, g_x_0_x_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxx_xxxxx[k] = -g_0_0_x_0_xx_xxxxx[k] - g_x_0_x_0_xx_xxxxx[k] * ab_x + g_x_0_x_0_xx_xxxxxx[k];

                g_x_0_x_0_xxx_xxxxy[k] = -g_0_0_x_0_xx_xxxxy[k] - g_x_0_x_0_xx_xxxxy[k] * ab_x + g_x_0_x_0_xx_xxxxxy[k];

                g_x_0_x_0_xxx_xxxxz[k] = -g_0_0_x_0_xx_xxxxz[k] - g_x_0_x_0_xx_xxxxz[k] * ab_x + g_x_0_x_0_xx_xxxxxz[k];

                g_x_0_x_0_xxx_xxxyy[k] = -g_0_0_x_0_xx_xxxyy[k] - g_x_0_x_0_xx_xxxyy[k] * ab_x + g_x_0_x_0_xx_xxxxyy[k];

                g_x_0_x_0_xxx_xxxyz[k] = -g_0_0_x_0_xx_xxxyz[k] - g_x_0_x_0_xx_xxxyz[k] * ab_x + g_x_0_x_0_xx_xxxxyz[k];

                g_x_0_x_0_xxx_xxxzz[k] = -g_0_0_x_0_xx_xxxzz[k] - g_x_0_x_0_xx_xxxzz[k] * ab_x + g_x_0_x_0_xx_xxxxzz[k];

                g_x_0_x_0_xxx_xxyyy[k] = -g_0_0_x_0_xx_xxyyy[k] - g_x_0_x_0_xx_xxyyy[k] * ab_x + g_x_0_x_0_xx_xxxyyy[k];

                g_x_0_x_0_xxx_xxyyz[k] = -g_0_0_x_0_xx_xxyyz[k] - g_x_0_x_0_xx_xxyyz[k] * ab_x + g_x_0_x_0_xx_xxxyyz[k];

                g_x_0_x_0_xxx_xxyzz[k] = -g_0_0_x_0_xx_xxyzz[k] - g_x_0_x_0_xx_xxyzz[k] * ab_x + g_x_0_x_0_xx_xxxyzz[k];

                g_x_0_x_0_xxx_xxzzz[k] = -g_0_0_x_0_xx_xxzzz[k] - g_x_0_x_0_xx_xxzzz[k] * ab_x + g_x_0_x_0_xx_xxxzzz[k];

                g_x_0_x_0_xxx_xyyyy[k] = -g_0_0_x_0_xx_xyyyy[k] - g_x_0_x_0_xx_xyyyy[k] * ab_x + g_x_0_x_0_xx_xxyyyy[k];

                g_x_0_x_0_xxx_xyyyz[k] = -g_0_0_x_0_xx_xyyyz[k] - g_x_0_x_0_xx_xyyyz[k] * ab_x + g_x_0_x_0_xx_xxyyyz[k];

                g_x_0_x_0_xxx_xyyzz[k] = -g_0_0_x_0_xx_xyyzz[k] - g_x_0_x_0_xx_xyyzz[k] * ab_x + g_x_0_x_0_xx_xxyyzz[k];

                g_x_0_x_0_xxx_xyzzz[k] = -g_0_0_x_0_xx_xyzzz[k] - g_x_0_x_0_xx_xyzzz[k] * ab_x + g_x_0_x_0_xx_xxyzzz[k];

                g_x_0_x_0_xxx_xzzzz[k] = -g_0_0_x_0_xx_xzzzz[k] - g_x_0_x_0_xx_xzzzz[k] * ab_x + g_x_0_x_0_xx_xxzzzz[k];

                g_x_0_x_0_xxx_yyyyy[k] = -g_0_0_x_0_xx_yyyyy[k] - g_x_0_x_0_xx_yyyyy[k] * ab_x + g_x_0_x_0_xx_xyyyyy[k];

                g_x_0_x_0_xxx_yyyyz[k] = -g_0_0_x_0_xx_yyyyz[k] - g_x_0_x_0_xx_yyyyz[k] * ab_x + g_x_0_x_0_xx_xyyyyz[k];

                g_x_0_x_0_xxx_yyyzz[k] = -g_0_0_x_0_xx_yyyzz[k] - g_x_0_x_0_xx_yyyzz[k] * ab_x + g_x_0_x_0_xx_xyyyzz[k];

                g_x_0_x_0_xxx_yyzzz[k] = -g_0_0_x_0_xx_yyzzz[k] - g_x_0_x_0_xx_yyzzz[k] * ab_x + g_x_0_x_0_xx_xyyzzz[k];

                g_x_0_x_0_xxx_yzzzz[k] = -g_0_0_x_0_xx_yzzzz[k] - g_x_0_x_0_xx_yzzzz[k] * ab_x + g_x_0_x_0_xx_xyzzzz[k];

                g_x_0_x_0_xxx_zzzzz[k] = -g_0_0_x_0_xx_zzzzz[k] - g_x_0_x_0_xx_zzzzz[k] * ab_x + g_x_0_x_0_xx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xx_xxxxx, g_x_0_x_0_xx_xxxxxy, g_x_0_x_0_xx_xxxxy, g_x_0_x_0_xx_xxxxyy, g_x_0_x_0_xx_xxxxyz, g_x_0_x_0_xx_xxxxz, g_x_0_x_0_xx_xxxyy, g_x_0_x_0_xx_xxxyyy, g_x_0_x_0_xx_xxxyyz, g_x_0_x_0_xx_xxxyz, g_x_0_x_0_xx_xxxyzz, g_x_0_x_0_xx_xxxzz, g_x_0_x_0_xx_xxyyy, g_x_0_x_0_xx_xxyyyy, g_x_0_x_0_xx_xxyyyz, g_x_0_x_0_xx_xxyyz, g_x_0_x_0_xx_xxyyzz, g_x_0_x_0_xx_xxyzz, g_x_0_x_0_xx_xxyzzz, g_x_0_x_0_xx_xxzzz, g_x_0_x_0_xx_xyyyy, g_x_0_x_0_xx_xyyyyy, g_x_0_x_0_xx_xyyyyz, g_x_0_x_0_xx_xyyyz, g_x_0_x_0_xx_xyyyzz, g_x_0_x_0_xx_xyyzz, g_x_0_x_0_xx_xyyzzz, g_x_0_x_0_xx_xyzzz, g_x_0_x_0_xx_xyzzzz, g_x_0_x_0_xx_xzzzz, g_x_0_x_0_xx_yyyyy, g_x_0_x_0_xx_yyyyyy, g_x_0_x_0_xx_yyyyyz, g_x_0_x_0_xx_yyyyz, g_x_0_x_0_xx_yyyyzz, g_x_0_x_0_xx_yyyzz, g_x_0_x_0_xx_yyyzzz, g_x_0_x_0_xx_yyzzz, g_x_0_x_0_xx_yyzzzz, g_x_0_x_0_xx_yzzzz, g_x_0_x_0_xx_yzzzzz, g_x_0_x_0_xx_zzzzz, g_x_0_x_0_xxy_xxxxx, g_x_0_x_0_xxy_xxxxy, g_x_0_x_0_xxy_xxxxz, g_x_0_x_0_xxy_xxxyy, g_x_0_x_0_xxy_xxxyz, g_x_0_x_0_xxy_xxxzz, g_x_0_x_0_xxy_xxyyy, g_x_0_x_0_xxy_xxyyz, g_x_0_x_0_xxy_xxyzz, g_x_0_x_0_xxy_xxzzz, g_x_0_x_0_xxy_xyyyy, g_x_0_x_0_xxy_xyyyz, g_x_0_x_0_xxy_xyyzz, g_x_0_x_0_xxy_xyzzz, g_x_0_x_0_xxy_xzzzz, g_x_0_x_0_xxy_yyyyy, g_x_0_x_0_xxy_yyyyz, g_x_0_x_0_xxy_yyyzz, g_x_0_x_0_xxy_yyzzz, g_x_0_x_0_xxy_yzzzz, g_x_0_x_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxy_xxxxx[k] = -g_x_0_x_0_xx_xxxxx[k] * ab_y + g_x_0_x_0_xx_xxxxxy[k];

                g_x_0_x_0_xxy_xxxxy[k] = -g_x_0_x_0_xx_xxxxy[k] * ab_y + g_x_0_x_0_xx_xxxxyy[k];

                g_x_0_x_0_xxy_xxxxz[k] = -g_x_0_x_0_xx_xxxxz[k] * ab_y + g_x_0_x_0_xx_xxxxyz[k];

                g_x_0_x_0_xxy_xxxyy[k] = -g_x_0_x_0_xx_xxxyy[k] * ab_y + g_x_0_x_0_xx_xxxyyy[k];

                g_x_0_x_0_xxy_xxxyz[k] = -g_x_0_x_0_xx_xxxyz[k] * ab_y + g_x_0_x_0_xx_xxxyyz[k];

                g_x_0_x_0_xxy_xxxzz[k] = -g_x_0_x_0_xx_xxxzz[k] * ab_y + g_x_0_x_0_xx_xxxyzz[k];

                g_x_0_x_0_xxy_xxyyy[k] = -g_x_0_x_0_xx_xxyyy[k] * ab_y + g_x_0_x_0_xx_xxyyyy[k];

                g_x_0_x_0_xxy_xxyyz[k] = -g_x_0_x_0_xx_xxyyz[k] * ab_y + g_x_0_x_0_xx_xxyyyz[k];

                g_x_0_x_0_xxy_xxyzz[k] = -g_x_0_x_0_xx_xxyzz[k] * ab_y + g_x_0_x_0_xx_xxyyzz[k];

                g_x_0_x_0_xxy_xxzzz[k] = -g_x_0_x_0_xx_xxzzz[k] * ab_y + g_x_0_x_0_xx_xxyzzz[k];

                g_x_0_x_0_xxy_xyyyy[k] = -g_x_0_x_0_xx_xyyyy[k] * ab_y + g_x_0_x_0_xx_xyyyyy[k];

                g_x_0_x_0_xxy_xyyyz[k] = -g_x_0_x_0_xx_xyyyz[k] * ab_y + g_x_0_x_0_xx_xyyyyz[k];

                g_x_0_x_0_xxy_xyyzz[k] = -g_x_0_x_0_xx_xyyzz[k] * ab_y + g_x_0_x_0_xx_xyyyzz[k];

                g_x_0_x_0_xxy_xyzzz[k] = -g_x_0_x_0_xx_xyzzz[k] * ab_y + g_x_0_x_0_xx_xyyzzz[k];

                g_x_0_x_0_xxy_xzzzz[k] = -g_x_0_x_0_xx_xzzzz[k] * ab_y + g_x_0_x_0_xx_xyzzzz[k];

                g_x_0_x_0_xxy_yyyyy[k] = -g_x_0_x_0_xx_yyyyy[k] * ab_y + g_x_0_x_0_xx_yyyyyy[k];

                g_x_0_x_0_xxy_yyyyz[k] = -g_x_0_x_0_xx_yyyyz[k] * ab_y + g_x_0_x_0_xx_yyyyyz[k];

                g_x_0_x_0_xxy_yyyzz[k] = -g_x_0_x_0_xx_yyyzz[k] * ab_y + g_x_0_x_0_xx_yyyyzz[k];

                g_x_0_x_0_xxy_yyzzz[k] = -g_x_0_x_0_xx_yyzzz[k] * ab_y + g_x_0_x_0_xx_yyyzzz[k];

                g_x_0_x_0_xxy_yzzzz[k] = -g_x_0_x_0_xx_yzzzz[k] * ab_y + g_x_0_x_0_xx_yyzzzz[k];

                g_x_0_x_0_xxy_zzzzz[k] = -g_x_0_x_0_xx_zzzzz[k] * ab_y + g_x_0_x_0_xx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xx_xxxxx, g_x_0_x_0_xx_xxxxxz, g_x_0_x_0_xx_xxxxy, g_x_0_x_0_xx_xxxxyz, g_x_0_x_0_xx_xxxxz, g_x_0_x_0_xx_xxxxzz, g_x_0_x_0_xx_xxxyy, g_x_0_x_0_xx_xxxyyz, g_x_0_x_0_xx_xxxyz, g_x_0_x_0_xx_xxxyzz, g_x_0_x_0_xx_xxxzz, g_x_0_x_0_xx_xxxzzz, g_x_0_x_0_xx_xxyyy, g_x_0_x_0_xx_xxyyyz, g_x_0_x_0_xx_xxyyz, g_x_0_x_0_xx_xxyyzz, g_x_0_x_0_xx_xxyzz, g_x_0_x_0_xx_xxyzzz, g_x_0_x_0_xx_xxzzz, g_x_0_x_0_xx_xxzzzz, g_x_0_x_0_xx_xyyyy, g_x_0_x_0_xx_xyyyyz, g_x_0_x_0_xx_xyyyz, g_x_0_x_0_xx_xyyyzz, g_x_0_x_0_xx_xyyzz, g_x_0_x_0_xx_xyyzzz, g_x_0_x_0_xx_xyzzz, g_x_0_x_0_xx_xyzzzz, g_x_0_x_0_xx_xzzzz, g_x_0_x_0_xx_xzzzzz, g_x_0_x_0_xx_yyyyy, g_x_0_x_0_xx_yyyyyz, g_x_0_x_0_xx_yyyyz, g_x_0_x_0_xx_yyyyzz, g_x_0_x_0_xx_yyyzz, g_x_0_x_0_xx_yyyzzz, g_x_0_x_0_xx_yyzzz, g_x_0_x_0_xx_yyzzzz, g_x_0_x_0_xx_yzzzz, g_x_0_x_0_xx_yzzzzz, g_x_0_x_0_xx_zzzzz, g_x_0_x_0_xx_zzzzzz, g_x_0_x_0_xxz_xxxxx, g_x_0_x_0_xxz_xxxxy, g_x_0_x_0_xxz_xxxxz, g_x_0_x_0_xxz_xxxyy, g_x_0_x_0_xxz_xxxyz, g_x_0_x_0_xxz_xxxzz, g_x_0_x_0_xxz_xxyyy, g_x_0_x_0_xxz_xxyyz, g_x_0_x_0_xxz_xxyzz, g_x_0_x_0_xxz_xxzzz, g_x_0_x_0_xxz_xyyyy, g_x_0_x_0_xxz_xyyyz, g_x_0_x_0_xxz_xyyzz, g_x_0_x_0_xxz_xyzzz, g_x_0_x_0_xxz_xzzzz, g_x_0_x_0_xxz_yyyyy, g_x_0_x_0_xxz_yyyyz, g_x_0_x_0_xxz_yyyzz, g_x_0_x_0_xxz_yyzzz, g_x_0_x_0_xxz_yzzzz, g_x_0_x_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxz_xxxxx[k] = -g_x_0_x_0_xx_xxxxx[k] * ab_z + g_x_0_x_0_xx_xxxxxz[k];

                g_x_0_x_0_xxz_xxxxy[k] = -g_x_0_x_0_xx_xxxxy[k] * ab_z + g_x_0_x_0_xx_xxxxyz[k];

                g_x_0_x_0_xxz_xxxxz[k] = -g_x_0_x_0_xx_xxxxz[k] * ab_z + g_x_0_x_0_xx_xxxxzz[k];

                g_x_0_x_0_xxz_xxxyy[k] = -g_x_0_x_0_xx_xxxyy[k] * ab_z + g_x_0_x_0_xx_xxxyyz[k];

                g_x_0_x_0_xxz_xxxyz[k] = -g_x_0_x_0_xx_xxxyz[k] * ab_z + g_x_0_x_0_xx_xxxyzz[k];

                g_x_0_x_0_xxz_xxxzz[k] = -g_x_0_x_0_xx_xxxzz[k] * ab_z + g_x_0_x_0_xx_xxxzzz[k];

                g_x_0_x_0_xxz_xxyyy[k] = -g_x_0_x_0_xx_xxyyy[k] * ab_z + g_x_0_x_0_xx_xxyyyz[k];

                g_x_0_x_0_xxz_xxyyz[k] = -g_x_0_x_0_xx_xxyyz[k] * ab_z + g_x_0_x_0_xx_xxyyzz[k];

                g_x_0_x_0_xxz_xxyzz[k] = -g_x_0_x_0_xx_xxyzz[k] * ab_z + g_x_0_x_0_xx_xxyzzz[k];

                g_x_0_x_0_xxz_xxzzz[k] = -g_x_0_x_0_xx_xxzzz[k] * ab_z + g_x_0_x_0_xx_xxzzzz[k];

                g_x_0_x_0_xxz_xyyyy[k] = -g_x_0_x_0_xx_xyyyy[k] * ab_z + g_x_0_x_0_xx_xyyyyz[k];

                g_x_0_x_0_xxz_xyyyz[k] = -g_x_0_x_0_xx_xyyyz[k] * ab_z + g_x_0_x_0_xx_xyyyzz[k];

                g_x_0_x_0_xxz_xyyzz[k] = -g_x_0_x_0_xx_xyyzz[k] * ab_z + g_x_0_x_0_xx_xyyzzz[k];

                g_x_0_x_0_xxz_xyzzz[k] = -g_x_0_x_0_xx_xyzzz[k] * ab_z + g_x_0_x_0_xx_xyzzzz[k];

                g_x_0_x_0_xxz_xzzzz[k] = -g_x_0_x_0_xx_xzzzz[k] * ab_z + g_x_0_x_0_xx_xzzzzz[k];

                g_x_0_x_0_xxz_yyyyy[k] = -g_x_0_x_0_xx_yyyyy[k] * ab_z + g_x_0_x_0_xx_yyyyyz[k];

                g_x_0_x_0_xxz_yyyyz[k] = -g_x_0_x_0_xx_yyyyz[k] * ab_z + g_x_0_x_0_xx_yyyyzz[k];

                g_x_0_x_0_xxz_yyyzz[k] = -g_x_0_x_0_xx_yyyzz[k] * ab_z + g_x_0_x_0_xx_yyyzzz[k];

                g_x_0_x_0_xxz_yyzzz[k] = -g_x_0_x_0_xx_yyzzz[k] * ab_z + g_x_0_x_0_xx_yyzzzz[k];

                g_x_0_x_0_xxz_yzzzz[k] = -g_x_0_x_0_xx_yzzzz[k] * ab_z + g_x_0_x_0_xx_yzzzzz[k];

                g_x_0_x_0_xxz_zzzzz[k] = -g_x_0_x_0_xx_zzzzz[k] * ab_z + g_x_0_x_0_xx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xy_xxxxx, g_x_0_x_0_xy_xxxxxy, g_x_0_x_0_xy_xxxxy, g_x_0_x_0_xy_xxxxyy, g_x_0_x_0_xy_xxxxyz, g_x_0_x_0_xy_xxxxz, g_x_0_x_0_xy_xxxyy, g_x_0_x_0_xy_xxxyyy, g_x_0_x_0_xy_xxxyyz, g_x_0_x_0_xy_xxxyz, g_x_0_x_0_xy_xxxyzz, g_x_0_x_0_xy_xxxzz, g_x_0_x_0_xy_xxyyy, g_x_0_x_0_xy_xxyyyy, g_x_0_x_0_xy_xxyyyz, g_x_0_x_0_xy_xxyyz, g_x_0_x_0_xy_xxyyzz, g_x_0_x_0_xy_xxyzz, g_x_0_x_0_xy_xxyzzz, g_x_0_x_0_xy_xxzzz, g_x_0_x_0_xy_xyyyy, g_x_0_x_0_xy_xyyyyy, g_x_0_x_0_xy_xyyyyz, g_x_0_x_0_xy_xyyyz, g_x_0_x_0_xy_xyyyzz, g_x_0_x_0_xy_xyyzz, g_x_0_x_0_xy_xyyzzz, g_x_0_x_0_xy_xyzzz, g_x_0_x_0_xy_xyzzzz, g_x_0_x_0_xy_xzzzz, g_x_0_x_0_xy_yyyyy, g_x_0_x_0_xy_yyyyyy, g_x_0_x_0_xy_yyyyyz, g_x_0_x_0_xy_yyyyz, g_x_0_x_0_xy_yyyyzz, g_x_0_x_0_xy_yyyzz, g_x_0_x_0_xy_yyyzzz, g_x_0_x_0_xy_yyzzz, g_x_0_x_0_xy_yyzzzz, g_x_0_x_0_xy_yzzzz, g_x_0_x_0_xy_yzzzzz, g_x_0_x_0_xy_zzzzz, g_x_0_x_0_xyy_xxxxx, g_x_0_x_0_xyy_xxxxy, g_x_0_x_0_xyy_xxxxz, g_x_0_x_0_xyy_xxxyy, g_x_0_x_0_xyy_xxxyz, g_x_0_x_0_xyy_xxxzz, g_x_0_x_0_xyy_xxyyy, g_x_0_x_0_xyy_xxyyz, g_x_0_x_0_xyy_xxyzz, g_x_0_x_0_xyy_xxzzz, g_x_0_x_0_xyy_xyyyy, g_x_0_x_0_xyy_xyyyz, g_x_0_x_0_xyy_xyyzz, g_x_0_x_0_xyy_xyzzz, g_x_0_x_0_xyy_xzzzz, g_x_0_x_0_xyy_yyyyy, g_x_0_x_0_xyy_yyyyz, g_x_0_x_0_xyy_yyyzz, g_x_0_x_0_xyy_yyzzz, g_x_0_x_0_xyy_yzzzz, g_x_0_x_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyy_xxxxx[k] = -g_x_0_x_0_xy_xxxxx[k] * ab_y + g_x_0_x_0_xy_xxxxxy[k];

                g_x_0_x_0_xyy_xxxxy[k] = -g_x_0_x_0_xy_xxxxy[k] * ab_y + g_x_0_x_0_xy_xxxxyy[k];

                g_x_0_x_0_xyy_xxxxz[k] = -g_x_0_x_0_xy_xxxxz[k] * ab_y + g_x_0_x_0_xy_xxxxyz[k];

                g_x_0_x_0_xyy_xxxyy[k] = -g_x_0_x_0_xy_xxxyy[k] * ab_y + g_x_0_x_0_xy_xxxyyy[k];

                g_x_0_x_0_xyy_xxxyz[k] = -g_x_0_x_0_xy_xxxyz[k] * ab_y + g_x_0_x_0_xy_xxxyyz[k];

                g_x_0_x_0_xyy_xxxzz[k] = -g_x_0_x_0_xy_xxxzz[k] * ab_y + g_x_0_x_0_xy_xxxyzz[k];

                g_x_0_x_0_xyy_xxyyy[k] = -g_x_0_x_0_xy_xxyyy[k] * ab_y + g_x_0_x_0_xy_xxyyyy[k];

                g_x_0_x_0_xyy_xxyyz[k] = -g_x_0_x_0_xy_xxyyz[k] * ab_y + g_x_0_x_0_xy_xxyyyz[k];

                g_x_0_x_0_xyy_xxyzz[k] = -g_x_0_x_0_xy_xxyzz[k] * ab_y + g_x_0_x_0_xy_xxyyzz[k];

                g_x_0_x_0_xyy_xxzzz[k] = -g_x_0_x_0_xy_xxzzz[k] * ab_y + g_x_0_x_0_xy_xxyzzz[k];

                g_x_0_x_0_xyy_xyyyy[k] = -g_x_0_x_0_xy_xyyyy[k] * ab_y + g_x_0_x_0_xy_xyyyyy[k];

                g_x_0_x_0_xyy_xyyyz[k] = -g_x_0_x_0_xy_xyyyz[k] * ab_y + g_x_0_x_0_xy_xyyyyz[k];

                g_x_0_x_0_xyy_xyyzz[k] = -g_x_0_x_0_xy_xyyzz[k] * ab_y + g_x_0_x_0_xy_xyyyzz[k];

                g_x_0_x_0_xyy_xyzzz[k] = -g_x_0_x_0_xy_xyzzz[k] * ab_y + g_x_0_x_0_xy_xyyzzz[k];

                g_x_0_x_0_xyy_xzzzz[k] = -g_x_0_x_0_xy_xzzzz[k] * ab_y + g_x_0_x_0_xy_xyzzzz[k];

                g_x_0_x_0_xyy_yyyyy[k] = -g_x_0_x_0_xy_yyyyy[k] * ab_y + g_x_0_x_0_xy_yyyyyy[k];

                g_x_0_x_0_xyy_yyyyz[k] = -g_x_0_x_0_xy_yyyyz[k] * ab_y + g_x_0_x_0_xy_yyyyyz[k];

                g_x_0_x_0_xyy_yyyzz[k] = -g_x_0_x_0_xy_yyyzz[k] * ab_y + g_x_0_x_0_xy_yyyyzz[k];

                g_x_0_x_0_xyy_yyzzz[k] = -g_x_0_x_0_xy_yyzzz[k] * ab_y + g_x_0_x_0_xy_yyyzzz[k];

                g_x_0_x_0_xyy_yzzzz[k] = -g_x_0_x_0_xy_yzzzz[k] * ab_y + g_x_0_x_0_xy_yyzzzz[k];

                g_x_0_x_0_xyy_zzzzz[k] = -g_x_0_x_0_xy_zzzzz[k] * ab_y + g_x_0_x_0_xy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyz_xxxxx, g_x_0_x_0_xyz_xxxxy, g_x_0_x_0_xyz_xxxxz, g_x_0_x_0_xyz_xxxyy, g_x_0_x_0_xyz_xxxyz, g_x_0_x_0_xyz_xxxzz, g_x_0_x_0_xyz_xxyyy, g_x_0_x_0_xyz_xxyyz, g_x_0_x_0_xyz_xxyzz, g_x_0_x_0_xyz_xxzzz, g_x_0_x_0_xyz_xyyyy, g_x_0_x_0_xyz_xyyyz, g_x_0_x_0_xyz_xyyzz, g_x_0_x_0_xyz_xyzzz, g_x_0_x_0_xyz_xzzzz, g_x_0_x_0_xyz_yyyyy, g_x_0_x_0_xyz_yyyyz, g_x_0_x_0_xyz_yyyzz, g_x_0_x_0_xyz_yyzzz, g_x_0_x_0_xyz_yzzzz, g_x_0_x_0_xyz_zzzzz, g_x_0_x_0_xz_xxxxx, g_x_0_x_0_xz_xxxxxy, g_x_0_x_0_xz_xxxxy, g_x_0_x_0_xz_xxxxyy, g_x_0_x_0_xz_xxxxyz, g_x_0_x_0_xz_xxxxz, g_x_0_x_0_xz_xxxyy, g_x_0_x_0_xz_xxxyyy, g_x_0_x_0_xz_xxxyyz, g_x_0_x_0_xz_xxxyz, g_x_0_x_0_xz_xxxyzz, g_x_0_x_0_xz_xxxzz, g_x_0_x_0_xz_xxyyy, g_x_0_x_0_xz_xxyyyy, g_x_0_x_0_xz_xxyyyz, g_x_0_x_0_xz_xxyyz, g_x_0_x_0_xz_xxyyzz, g_x_0_x_0_xz_xxyzz, g_x_0_x_0_xz_xxyzzz, g_x_0_x_0_xz_xxzzz, g_x_0_x_0_xz_xyyyy, g_x_0_x_0_xz_xyyyyy, g_x_0_x_0_xz_xyyyyz, g_x_0_x_0_xz_xyyyz, g_x_0_x_0_xz_xyyyzz, g_x_0_x_0_xz_xyyzz, g_x_0_x_0_xz_xyyzzz, g_x_0_x_0_xz_xyzzz, g_x_0_x_0_xz_xyzzzz, g_x_0_x_0_xz_xzzzz, g_x_0_x_0_xz_yyyyy, g_x_0_x_0_xz_yyyyyy, g_x_0_x_0_xz_yyyyyz, g_x_0_x_0_xz_yyyyz, g_x_0_x_0_xz_yyyyzz, g_x_0_x_0_xz_yyyzz, g_x_0_x_0_xz_yyyzzz, g_x_0_x_0_xz_yyzzz, g_x_0_x_0_xz_yyzzzz, g_x_0_x_0_xz_yzzzz, g_x_0_x_0_xz_yzzzzz, g_x_0_x_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyz_xxxxx[k] = -g_x_0_x_0_xz_xxxxx[k] * ab_y + g_x_0_x_0_xz_xxxxxy[k];

                g_x_0_x_0_xyz_xxxxy[k] = -g_x_0_x_0_xz_xxxxy[k] * ab_y + g_x_0_x_0_xz_xxxxyy[k];

                g_x_0_x_0_xyz_xxxxz[k] = -g_x_0_x_0_xz_xxxxz[k] * ab_y + g_x_0_x_0_xz_xxxxyz[k];

                g_x_0_x_0_xyz_xxxyy[k] = -g_x_0_x_0_xz_xxxyy[k] * ab_y + g_x_0_x_0_xz_xxxyyy[k];

                g_x_0_x_0_xyz_xxxyz[k] = -g_x_0_x_0_xz_xxxyz[k] * ab_y + g_x_0_x_0_xz_xxxyyz[k];

                g_x_0_x_0_xyz_xxxzz[k] = -g_x_0_x_0_xz_xxxzz[k] * ab_y + g_x_0_x_0_xz_xxxyzz[k];

                g_x_0_x_0_xyz_xxyyy[k] = -g_x_0_x_0_xz_xxyyy[k] * ab_y + g_x_0_x_0_xz_xxyyyy[k];

                g_x_0_x_0_xyz_xxyyz[k] = -g_x_0_x_0_xz_xxyyz[k] * ab_y + g_x_0_x_0_xz_xxyyyz[k];

                g_x_0_x_0_xyz_xxyzz[k] = -g_x_0_x_0_xz_xxyzz[k] * ab_y + g_x_0_x_0_xz_xxyyzz[k];

                g_x_0_x_0_xyz_xxzzz[k] = -g_x_0_x_0_xz_xxzzz[k] * ab_y + g_x_0_x_0_xz_xxyzzz[k];

                g_x_0_x_0_xyz_xyyyy[k] = -g_x_0_x_0_xz_xyyyy[k] * ab_y + g_x_0_x_0_xz_xyyyyy[k];

                g_x_0_x_0_xyz_xyyyz[k] = -g_x_0_x_0_xz_xyyyz[k] * ab_y + g_x_0_x_0_xz_xyyyyz[k];

                g_x_0_x_0_xyz_xyyzz[k] = -g_x_0_x_0_xz_xyyzz[k] * ab_y + g_x_0_x_0_xz_xyyyzz[k];

                g_x_0_x_0_xyz_xyzzz[k] = -g_x_0_x_0_xz_xyzzz[k] * ab_y + g_x_0_x_0_xz_xyyzzz[k];

                g_x_0_x_0_xyz_xzzzz[k] = -g_x_0_x_0_xz_xzzzz[k] * ab_y + g_x_0_x_0_xz_xyzzzz[k];

                g_x_0_x_0_xyz_yyyyy[k] = -g_x_0_x_0_xz_yyyyy[k] * ab_y + g_x_0_x_0_xz_yyyyyy[k];

                g_x_0_x_0_xyz_yyyyz[k] = -g_x_0_x_0_xz_yyyyz[k] * ab_y + g_x_0_x_0_xz_yyyyyz[k];

                g_x_0_x_0_xyz_yyyzz[k] = -g_x_0_x_0_xz_yyyzz[k] * ab_y + g_x_0_x_0_xz_yyyyzz[k];

                g_x_0_x_0_xyz_yyzzz[k] = -g_x_0_x_0_xz_yyzzz[k] * ab_y + g_x_0_x_0_xz_yyyzzz[k];

                g_x_0_x_0_xyz_yzzzz[k] = -g_x_0_x_0_xz_yzzzz[k] * ab_y + g_x_0_x_0_xz_yyzzzz[k];

                g_x_0_x_0_xyz_zzzzz[k] = -g_x_0_x_0_xz_zzzzz[k] * ab_y + g_x_0_x_0_xz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xz_xxxxx, g_x_0_x_0_xz_xxxxxz, g_x_0_x_0_xz_xxxxy, g_x_0_x_0_xz_xxxxyz, g_x_0_x_0_xz_xxxxz, g_x_0_x_0_xz_xxxxzz, g_x_0_x_0_xz_xxxyy, g_x_0_x_0_xz_xxxyyz, g_x_0_x_0_xz_xxxyz, g_x_0_x_0_xz_xxxyzz, g_x_0_x_0_xz_xxxzz, g_x_0_x_0_xz_xxxzzz, g_x_0_x_0_xz_xxyyy, g_x_0_x_0_xz_xxyyyz, g_x_0_x_0_xz_xxyyz, g_x_0_x_0_xz_xxyyzz, g_x_0_x_0_xz_xxyzz, g_x_0_x_0_xz_xxyzzz, g_x_0_x_0_xz_xxzzz, g_x_0_x_0_xz_xxzzzz, g_x_0_x_0_xz_xyyyy, g_x_0_x_0_xz_xyyyyz, g_x_0_x_0_xz_xyyyz, g_x_0_x_0_xz_xyyyzz, g_x_0_x_0_xz_xyyzz, g_x_0_x_0_xz_xyyzzz, g_x_0_x_0_xz_xyzzz, g_x_0_x_0_xz_xyzzzz, g_x_0_x_0_xz_xzzzz, g_x_0_x_0_xz_xzzzzz, g_x_0_x_0_xz_yyyyy, g_x_0_x_0_xz_yyyyyz, g_x_0_x_0_xz_yyyyz, g_x_0_x_0_xz_yyyyzz, g_x_0_x_0_xz_yyyzz, g_x_0_x_0_xz_yyyzzz, g_x_0_x_0_xz_yyzzz, g_x_0_x_0_xz_yyzzzz, g_x_0_x_0_xz_yzzzz, g_x_0_x_0_xz_yzzzzz, g_x_0_x_0_xz_zzzzz, g_x_0_x_0_xz_zzzzzz, g_x_0_x_0_xzz_xxxxx, g_x_0_x_0_xzz_xxxxy, g_x_0_x_0_xzz_xxxxz, g_x_0_x_0_xzz_xxxyy, g_x_0_x_0_xzz_xxxyz, g_x_0_x_0_xzz_xxxzz, g_x_0_x_0_xzz_xxyyy, g_x_0_x_0_xzz_xxyyz, g_x_0_x_0_xzz_xxyzz, g_x_0_x_0_xzz_xxzzz, g_x_0_x_0_xzz_xyyyy, g_x_0_x_0_xzz_xyyyz, g_x_0_x_0_xzz_xyyzz, g_x_0_x_0_xzz_xyzzz, g_x_0_x_0_xzz_xzzzz, g_x_0_x_0_xzz_yyyyy, g_x_0_x_0_xzz_yyyyz, g_x_0_x_0_xzz_yyyzz, g_x_0_x_0_xzz_yyzzz, g_x_0_x_0_xzz_yzzzz, g_x_0_x_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzz_xxxxx[k] = -g_x_0_x_0_xz_xxxxx[k] * ab_z + g_x_0_x_0_xz_xxxxxz[k];

                g_x_0_x_0_xzz_xxxxy[k] = -g_x_0_x_0_xz_xxxxy[k] * ab_z + g_x_0_x_0_xz_xxxxyz[k];

                g_x_0_x_0_xzz_xxxxz[k] = -g_x_0_x_0_xz_xxxxz[k] * ab_z + g_x_0_x_0_xz_xxxxzz[k];

                g_x_0_x_0_xzz_xxxyy[k] = -g_x_0_x_0_xz_xxxyy[k] * ab_z + g_x_0_x_0_xz_xxxyyz[k];

                g_x_0_x_0_xzz_xxxyz[k] = -g_x_0_x_0_xz_xxxyz[k] * ab_z + g_x_0_x_0_xz_xxxyzz[k];

                g_x_0_x_0_xzz_xxxzz[k] = -g_x_0_x_0_xz_xxxzz[k] * ab_z + g_x_0_x_0_xz_xxxzzz[k];

                g_x_0_x_0_xzz_xxyyy[k] = -g_x_0_x_0_xz_xxyyy[k] * ab_z + g_x_0_x_0_xz_xxyyyz[k];

                g_x_0_x_0_xzz_xxyyz[k] = -g_x_0_x_0_xz_xxyyz[k] * ab_z + g_x_0_x_0_xz_xxyyzz[k];

                g_x_0_x_0_xzz_xxyzz[k] = -g_x_0_x_0_xz_xxyzz[k] * ab_z + g_x_0_x_0_xz_xxyzzz[k];

                g_x_0_x_0_xzz_xxzzz[k] = -g_x_0_x_0_xz_xxzzz[k] * ab_z + g_x_0_x_0_xz_xxzzzz[k];

                g_x_0_x_0_xzz_xyyyy[k] = -g_x_0_x_0_xz_xyyyy[k] * ab_z + g_x_0_x_0_xz_xyyyyz[k];

                g_x_0_x_0_xzz_xyyyz[k] = -g_x_0_x_0_xz_xyyyz[k] * ab_z + g_x_0_x_0_xz_xyyyzz[k];

                g_x_0_x_0_xzz_xyyzz[k] = -g_x_0_x_0_xz_xyyzz[k] * ab_z + g_x_0_x_0_xz_xyyzzz[k];

                g_x_0_x_0_xzz_xyzzz[k] = -g_x_0_x_0_xz_xyzzz[k] * ab_z + g_x_0_x_0_xz_xyzzzz[k];

                g_x_0_x_0_xzz_xzzzz[k] = -g_x_0_x_0_xz_xzzzz[k] * ab_z + g_x_0_x_0_xz_xzzzzz[k];

                g_x_0_x_0_xzz_yyyyy[k] = -g_x_0_x_0_xz_yyyyy[k] * ab_z + g_x_0_x_0_xz_yyyyyz[k];

                g_x_0_x_0_xzz_yyyyz[k] = -g_x_0_x_0_xz_yyyyz[k] * ab_z + g_x_0_x_0_xz_yyyyzz[k];

                g_x_0_x_0_xzz_yyyzz[k] = -g_x_0_x_0_xz_yyyzz[k] * ab_z + g_x_0_x_0_xz_yyyzzz[k];

                g_x_0_x_0_xzz_yyzzz[k] = -g_x_0_x_0_xz_yyzzz[k] * ab_z + g_x_0_x_0_xz_yyzzzz[k];

                g_x_0_x_0_xzz_yzzzz[k] = -g_x_0_x_0_xz_yzzzz[k] * ab_z + g_x_0_x_0_xz_yzzzzz[k];

                g_x_0_x_0_xzz_zzzzz[k] = -g_x_0_x_0_xz_zzzzz[k] * ab_z + g_x_0_x_0_xz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yy_xxxxx, g_x_0_x_0_yy_xxxxxy, g_x_0_x_0_yy_xxxxy, g_x_0_x_0_yy_xxxxyy, g_x_0_x_0_yy_xxxxyz, g_x_0_x_0_yy_xxxxz, g_x_0_x_0_yy_xxxyy, g_x_0_x_0_yy_xxxyyy, g_x_0_x_0_yy_xxxyyz, g_x_0_x_0_yy_xxxyz, g_x_0_x_0_yy_xxxyzz, g_x_0_x_0_yy_xxxzz, g_x_0_x_0_yy_xxyyy, g_x_0_x_0_yy_xxyyyy, g_x_0_x_0_yy_xxyyyz, g_x_0_x_0_yy_xxyyz, g_x_0_x_0_yy_xxyyzz, g_x_0_x_0_yy_xxyzz, g_x_0_x_0_yy_xxyzzz, g_x_0_x_0_yy_xxzzz, g_x_0_x_0_yy_xyyyy, g_x_0_x_0_yy_xyyyyy, g_x_0_x_0_yy_xyyyyz, g_x_0_x_0_yy_xyyyz, g_x_0_x_0_yy_xyyyzz, g_x_0_x_0_yy_xyyzz, g_x_0_x_0_yy_xyyzzz, g_x_0_x_0_yy_xyzzz, g_x_0_x_0_yy_xyzzzz, g_x_0_x_0_yy_xzzzz, g_x_0_x_0_yy_yyyyy, g_x_0_x_0_yy_yyyyyy, g_x_0_x_0_yy_yyyyyz, g_x_0_x_0_yy_yyyyz, g_x_0_x_0_yy_yyyyzz, g_x_0_x_0_yy_yyyzz, g_x_0_x_0_yy_yyyzzz, g_x_0_x_0_yy_yyzzz, g_x_0_x_0_yy_yyzzzz, g_x_0_x_0_yy_yzzzz, g_x_0_x_0_yy_yzzzzz, g_x_0_x_0_yy_zzzzz, g_x_0_x_0_yyy_xxxxx, g_x_0_x_0_yyy_xxxxy, g_x_0_x_0_yyy_xxxxz, g_x_0_x_0_yyy_xxxyy, g_x_0_x_0_yyy_xxxyz, g_x_0_x_0_yyy_xxxzz, g_x_0_x_0_yyy_xxyyy, g_x_0_x_0_yyy_xxyyz, g_x_0_x_0_yyy_xxyzz, g_x_0_x_0_yyy_xxzzz, g_x_0_x_0_yyy_xyyyy, g_x_0_x_0_yyy_xyyyz, g_x_0_x_0_yyy_xyyzz, g_x_0_x_0_yyy_xyzzz, g_x_0_x_0_yyy_xzzzz, g_x_0_x_0_yyy_yyyyy, g_x_0_x_0_yyy_yyyyz, g_x_0_x_0_yyy_yyyzz, g_x_0_x_0_yyy_yyzzz, g_x_0_x_0_yyy_yzzzz, g_x_0_x_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyy_xxxxx[k] = -g_x_0_x_0_yy_xxxxx[k] * ab_y + g_x_0_x_0_yy_xxxxxy[k];

                g_x_0_x_0_yyy_xxxxy[k] = -g_x_0_x_0_yy_xxxxy[k] * ab_y + g_x_0_x_0_yy_xxxxyy[k];

                g_x_0_x_0_yyy_xxxxz[k] = -g_x_0_x_0_yy_xxxxz[k] * ab_y + g_x_0_x_0_yy_xxxxyz[k];

                g_x_0_x_0_yyy_xxxyy[k] = -g_x_0_x_0_yy_xxxyy[k] * ab_y + g_x_0_x_0_yy_xxxyyy[k];

                g_x_0_x_0_yyy_xxxyz[k] = -g_x_0_x_0_yy_xxxyz[k] * ab_y + g_x_0_x_0_yy_xxxyyz[k];

                g_x_0_x_0_yyy_xxxzz[k] = -g_x_0_x_0_yy_xxxzz[k] * ab_y + g_x_0_x_0_yy_xxxyzz[k];

                g_x_0_x_0_yyy_xxyyy[k] = -g_x_0_x_0_yy_xxyyy[k] * ab_y + g_x_0_x_0_yy_xxyyyy[k];

                g_x_0_x_0_yyy_xxyyz[k] = -g_x_0_x_0_yy_xxyyz[k] * ab_y + g_x_0_x_0_yy_xxyyyz[k];

                g_x_0_x_0_yyy_xxyzz[k] = -g_x_0_x_0_yy_xxyzz[k] * ab_y + g_x_0_x_0_yy_xxyyzz[k];

                g_x_0_x_0_yyy_xxzzz[k] = -g_x_0_x_0_yy_xxzzz[k] * ab_y + g_x_0_x_0_yy_xxyzzz[k];

                g_x_0_x_0_yyy_xyyyy[k] = -g_x_0_x_0_yy_xyyyy[k] * ab_y + g_x_0_x_0_yy_xyyyyy[k];

                g_x_0_x_0_yyy_xyyyz[k] = -g_x_0_x_0_yy_xyyyz[k] * ab_y + g_x_0_x_0_yy_xyyyyz[k];

                g_x_0_x_0_yyy_xyyzz[k] = -g_x_0_x_0_yy_xyyzz[k] * ab_y + g_x_0_x_0_yy_xyyyzz[k];

                g_x_0_x_0_yyy_xyzzz[k] = -g_x_0_x_0_yy_xyzzz[k] * ab_y + g_x_0_x_0_yy_xyyzzz[k];

                g_x_0_x_0_yyy_xzzzz[k] = -g_x_0_x_0_yy_xzzzz[k] * ab_y + g_x_0_x_0_yy_xyzzzz[k];

                g_x_0_x_0_yyy_yyyyy[k] = -g_x_0_x_0_yy_yyyyy[k] * ab_y + g_x_0_x_0_yy_yyyyyy[k];

                g_x_0_x_0_yyy_yyyyz[k] = -g_x_0_x_0_yy_yyyyz[k] * ab_y + g_x_0_x_0_yy_yyyyyz[k];

                g_x_0_x_0_yyy_yyyzz[k] = -g_x_0_x_0_yy_yyyzz[k] * ab_y + g_x_0_x_0_yy_yyyyzz[k];

                g_x_0_x_0_yyy_yyzzz[k] = -g_x_0_x_0_yy_yyzzz[k] * ab_y + g_x_0_x_0_yy_yyyzzz[k];

                g_x_0_x_0_yyy_yzzzz[k] = -g_x_0_x_0_yy_yzzzz[k] * ab_y + g_x_0_x_0_yy_yyzzzz[k];

                g_x_0_x_0_yyy_zzzzz[k] = -g_x_0_x_0_yy_zzzzz[k] * ab_y + g_x_0_x_0_yy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyz_xxxxx, g_x_0_x_0_yyz_xxxxy, g_x_0_x_0_yyz_xxxxz, g_x_0_x_0_yyz_xxxyy, g_x_0_x_0_yyz_xxxyz, g_x_0_x_0_yyz_xxxzz, g_x_0_x_0_yyz_xxyyy, g_x_0_x_0_yyz_xxyyz, g_x_0_x_0_yyz_xxyzz, g_x_0_x_0_yyz_xxzzz, g_x_0_x_0_yyz_xyyyy, g_x_0_x_0_yyz_xyyyz, g_x_0_x_0_yyz_xyyzz, g_x_0_x_0_yyz_xyzzz, g_x_0_x_0_yyz_xzzzz, g_x_0_x_0_yyz_yyyyy, g_x_0_x_0_yyz_yyyyz, g_x_0_x_0_yyz_yyyzz, g_x_0_x_0_yyz_yyzzz, g_x_0_x_0_yyz_yzzzz, g_x_0_x_0_yyz_zzzzz, g_x_0_x_0_yz_xxxxx, g_x_0_x_0_yz_xxxxxy, g_x_0_x_0_yz_xxxxy, g_x_0_x_0_yz_xxxxyy, g_x_0_x_0_yz_xxxxyz, g_x_0_x_0_yz_xxxxz, g_x_0_x_0_yz_xxxyy, g_x_0_x_0_yz_xxxyyy, g_x_0_x_0_yz_xxxyyz, g_x_0_x_0_yz_xxxyz, g_x_0_x_0_yz_xxxyzz, g_x_0_x_0_yz_xxxzz, g_x_0_x_0_yz_xxyyy, g_x_0_x_0_yz_xxyyyy, g_x_0_x_0_yz_xxyyyz, g_x_0_x_0_yz_xxyyz, g_x_0_x_0_yz_xxyyzz, g_x_0_x_0_yz_xxyzz, g_x_0_x_0_yz_xxyzzz, g_x_0_x_0_yz_xxzzz, g_x_0_x_0_yz_xyyyy, g_x_0_x_0_yz_xyyyyy, g_x_0_x_0_yz_xyyyyz, g_x_0_x_0_yz_xyyyz, g_x_0_x_0_yz_xyyyzz, g_x_0_x_0_yz_xyyzz, g_x_0_x_0_yz_xyyzzz, g_x_0_x_0_yz_xyzzz, g_x_0_x_0_yz_xyzzzz, g_x_0_x_0_yz_xzzzz, g_x_0_x_0_yz_yyyyy, g_x_0_x_0_yz_yyyyyy, g_x_0_x_0_yz_yyyyyz, g_x_0_x_0_yz_yyyyz, g_x_0_x_0_yz_yyyyzz, g_x_0_x_0_yz_yyyzz, g_x_0_x_0_yz_yyyzzz, g_x_0_x_0_yz_yyzzz, g_x_0_x_0_yz_yyzzzz, g_x_0_x_0_yz_yzzzz, g_x_0_x_0_yz_yzzzzz, g_x_0_x_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyz_xxxxx[k] = -g_x_0_x_0_yz_xxxxx[k] * ab_y + g_x_0_x_0_yz_xxxxxy[k];

                g_x_0_x_0_yyz_xxxxy[k] = -g_x_0_x_0_yz_xxxxy[k] * ab_y + g_x_0_x_0_yz_xxxxyy[k];

                g_x_0_x_0_yyz_xxxxz[k] = -g_x_0_x_0_yz_xxxxz[k] * ab_y + g_x_0_x_0_yz_xxxxyz[k];

                g_x_0_x_0_yyz_xxxyy[k] = -g_x_0_x_0_yz_xxxyy[k] * ab_y + g_x_0_x_0_yz_xxxyyy[k];

                g_x_0_x_0_yyz_xxxyz[k] = -g_x_0_x_0_yz_xxxyz[k] * ab_y + g_x_0_x_0_yz_xxxyyz[k];

                g_x_0_x_0_yyz_xxxzz[k] = -g_x_0_x_0_yz_xxxzz[k] * ab_y + g_x_0_x_0_yz_xxxyzz[k];

                g_x_0_x_0_yyz_xxyyy[k] = -g_x_0_x_0_yz_xxyyy[k] * ab_y + g_x_0_x_0_yz_xxyyyy[k];

                g_x_0_x_0_yyz_xxyyz[k] = -g_x_0_x_0_yz_xxyyz[k] * ab_y + g_x_0_x_0_yz_xxyyyz[k];

                g_x_0_x_0_yyz_xxyzz[k] = -g_x_0_x_0_yz_xxyzz[k] * ab_y + g_x_0_x_0_yz_xxyyzz[k];

                g_x_0_x_0_yyz_xxzzz[k] = -g_x_0_x_0_yz_xxzzz[k] * ab_y + g_x_0_x_0_yz_xxyzzz[k];

                g_x_0_x_0_yyz_xyyyy[k] = -g_x_0_x_0_yz_xyyyy[k] * ab_y + g_x_0_x_0_yz_xyyyyy[k];

                g_x_0_x_0_yyz_xyyyz[k] = -g_x_0_x_0_yz_xyyyz[k] * ab_y + g_x_0_x_0_yz_xyyyyz[k];

                g_x_0_x_0_yyz_xyyzz[k] = -g_x_0_x_0_yz_xyyzz[k] * ab_y + g_x_0_x_0_yz_xyyyzz[k];

                g_x_0_x_0_yyz_xyzzz[k] = -g_x_0_x_0_yz_xyzzz[k] * ab_y + g_x_0_x_0_yz_xyyzzz[k];

                g_x_0_x_0_yyz_xzzzz[k] = -g_x_0_x_0_yz_xzzzz[k] * ab_y + g_x_0_x_0_yz_xyzzzz[k];

                g_x_0_x_0_yyz_yyyyy[k] = -g_x_0_x_0_yz_yyyyy[k] * ab_y + g_x_0_x_0_yz_yyyyyy[k];

                g_x_0_x_0_yyz_yyyyz[k] = -g_x_0_x_0_yz_yyyyz[k] * ab_y + g_x_0_x_0_yz_yyyyyz[k];

                g_x_0_x_0_yyz_yyyzz[k] = -g_x_0_x_0_yz_yyyzz[k] * ab_y + g_x_0_x_0_yz_yyyyzz[k];

                g_x_0_x_0_yyz_yyzzz[k] = -g_x_0_x_0_yz_yyzzz[k] * ab_y + g_x_0_x_0_yz_yyyzzz[k];

                g_x_0_x_0_yyz_yzzzz[k] = -g_x_0_x_0_yz_yzzzz[k] * ab_y + g_x_0_x_0_yz_yyzzzz[k];

                g_x_0_x_0_yyz_zzzzz[k] = -g_x_0_x_0_yz_zzzzz[k] * ab_y + g_x_0_x_0_yz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzz_xxxxx, g_x_0_x_0_yzz_xxxxy, g_x_0_x_0_yzz_xxxxz, g_x_0_x_0_yzz_xxxyy, g_x_0_x_0_yzz_xxxyz, g_x_0_x_0_yzz_xxxzz, g_x_0_x_0_yzz_xxyyy, g_x_0_x_0_yzz_xxyyz, g_x_0_x_0_yzz_xxyzz, g_x_0_x_0_yzz_xxzzz, g_x_0_x_0_yzz_xyyyy, g_x_0_x_0_yzz_xyyyz, g_x_0_x_0_yzz_xyyzz, g_x_0_x_0_yzz_xyzzz, g_x_0_x_0_yzz_xzzzz, g_x_0_x_0_yzz_yyyyy, g_x_0_x_0_yzz_yyyyz, g_x_0_x_0_yzz_yyyzz, g_x_0_x_0_yzz_yyzzz, g_x_0_x_0_yzz_yzzzz, g_x_0_x_0_yzz_zzzzz, g_x_0_x_0_zz_xxxxx, g_x_0_x_0_zz_xxxxxy, g_x_0_x_0_zz_xxxxy, g_x_0_x_0_zz_xxxxyy, g_x_0_x_0_zz_xxxxyz, g_x_0_x_0_zz_xxxxz, g_x_0_x_0_zz_xxxyy, g_x_0_x_0_zz_xxxyyy, g_x_0_x_0_zz_xxxyyz, g_x_0_x_0_zz_xxxyz, g_x_0_x_0_zz_xxxyzz, g_x_0_x_0_zz_xxxzz, g_x_0_x_0_zz_xxyyy, g_x_0_x_0_zz_xxyyyy, g_x_0_x_0_zz_xxyyyz, g_x_0_x_0_zz_xxyyz, g_x_0_x_0_zz_xxyyzz, g_x_0_x_0_zz_xxyzz, g_x_0_x_0_zz_xxyzzz, g_x_0_x_0_zz_xxzzz, g_x_0_x_0_zz_xyyyy, g_x_0_x_0_zz_xyyyyy, g_x_0_x_0_zz_xyyyyz, g_x_0_x_0_zz_xyyyz, g_x_0_x_0_zz_xyyyzz, g_x_0_x_0_zz_xyyzz, g_x_0_x_0_zz_xyyzzz, g_x_0_x_0_zz_xyzzz, g_x_0_x_0_zz_xyzzzz, g_x_0_x_0_zz_xzzzz, g_x_0_x_0_zz_yyyyy, g_x_0_x_0_zz_yyyyyy, g_x_0_x_0_zz_yyyyyz, g_x_0_x_0_zz_yyyyz, g_x_0_x_0_zz_yyyyzz, g_x_0_x_0_zz_yyyzz, g_x_0_x_0_zz_yyyzzz, g_x_0_x_0_zz_yyzzz, g_x_0_x_0_zz_yyzzzz, g_x_0_x_0_zz_yzzzz, g_x_0_x_0_zz_yzzzzz, g_x_0_x_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzz_xxxxx[k] = -g_x_0_x_0_zz_xxxxx[k] * ab_y + g_x_0_x_0_zz_xxxxxy[k];

                g_x_0_x_0_yzz_xxxxy[k] = -g_x_0_x_0_zz_xxxxy[k] * ab_y + g_x_0_x_0_zz_xxxxyy[k];

                g_x_0_x_0_yzz_xxxxz[k] = -g_x_0_x_0_zz_xxxxz[k] * ab_y + g_x_0_x_0_zz_xxxxyz[k];

                g_x_0_x_0_yzz_xxxyy[k] = -g_x_0_x_0_zz_xxxyy[k] * ab_y + g_x_0_x_0_zz_xxxyyy[k];

                g_x_0_x_0_yzz_xxxyz[k] = -g_x_0_x_0_zz_xxxyz[k] * ab_y + g_x_0_x_0_zz_xxxyyz[k];

                g_x_0_x_0_yzz_xxxzz[k] = -g_x_0_x_0_zz_xxxzz[k] * ab_y + g_x_0_x_0_zz_xxxyzz[k];

                g_x_0_x_0_yzz_xxyyy[k] = -g_x_0_x_0_zz_xxyyy[k] * ab_y + g_x_0_x_0_zz_xxyyyy[k];

                g_x_0_x_0_yzz_xxyyz[k] = -g_x_0_x_0_zz_xxyyz[k] * ab_y + g_x_0_x_0_zz_xxyyyz[k];

                g_x_0_x_0_yzz_xxyzz[k] = -g_x_0_x_0_zz_xxyzz[k] * ab_y + g_x_0_x_0_zz_xxyyzz[k];

                g_x_0_x_0_yzz_xxzzz[k] = -g_x_0_x_0_zz_xxzzz[k] * ab_y + g_x_0_x_0_zz_xxyzzz[k];

                g_x_0_x_0_yzz_xyyyy[k] = -g_x_0_x_0_zz_xyyyy[k] * ab_y + g_x_0_x_0_zz_xyyyyy[k];

                g_x_0_x_0_yzz_xyyyz[k] = -g_x_0_x_0_zz_xyyyz[k] * ab_y + g_x_0_x_0_zz_xyyyyz[k];

                g_x_0_x_0_yzz_xyyzz[k] = -g_x_0_x_0_zz_xyyzz[k] * ab_y + g_x_0_x_0_zz_xyyyzz[k];

                g_x_0_x_0_yzz_xyzzz[k] = -g_x_0_x_0_zz_xyzzz[k] * ab_y + g_x_0_x_0_zz_xyyzzz[k];

                g_x_0_x_0_yzz_xzzzz[k] = -g_x_0_x_0_zz_xzzzz[k] * ab_y + g_x_0_x_0_zz_xyzzzz[k];

                g_x_0_x_0_yzz_yyyyy[k] = -g_x_0_x_0_zz_yyyyy[k] * ab_y + g_x_0_x_0_zz_yyyyyy[k];

                g_x_0_x_0_yzz_yyyyz[k] = -g_x_0_x_0_zz_yyyyz[k] * ab_y + g_x_0_x_0_zz_yyyyyz[k];

                g_x_0_x_0_yzz_yyyzz[k] = -g_x_0_x_0_zz_yyyzz[k] * ab_y + g_x_0_x_0_zz_yyyyzz[k];

                g_x_0_x_0_yzz_yyzzz[k] = -g_x_0_x_0_zz_yyzzz[k] * ab_y + g_x_0_x_0_zz_yyyzzz[k];

                g_x_0_x_0_yzz_yzzzz[k] = -g_x_0_x_0_zz_yzzzz[k] * ab_y + g_x_0_x_0_zz_yyzzzz[k];

                g_x_0_x_0_yzz_zzzzz[k] = -g_x_0_x_0_zz_zzzzz[k] * ab_y + g_x_0_x_0_zz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 195 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zz_xxxxx, g_x_0_x_0_zz_xxxxxz, g_x_0_x_0_zz_xxxxy, g_x_0_x_0_zz_xxxxyz, g_x_0_x_0_zz_xxxxz, g_x_0_x_0_zz_xxxxzz, g_x_0_x_0_zz_xxxyy, g_x_0_x_0_zz_xxxyyz, g_x_0_x_0_zz_xxxyz, g_x_0_x_0_zz_xxxyzz, g_x_0_x_0_zz_xxxzz, g_x_0_x_0_zz_xxxzzz, g_x_0_x_0_zz_xxyyy, g_x_0_x_0_zz_xxyyyz, g_x_0_x_0_zz_xxyyz, g_x_0_x_0_zz_xxyyzz, g_x_0_x_0_zz_xxyzz, g_x_0_x_0_zz_xxyzzz, g_x_0_x_0_zz_xxzzz, g_x_0_x_0_zz_xxzzzz, g_x_0_x_0_zz_xyyyy, g_x_0_x_0_zz_xyyyyz, g_x_0_x_0_zz_xyyyz, g_x_0_x_0_zz_xyyyzz, g_x_0_x_0_zz_xyyzz, g_x_0_x_0_zz_xyyzzz, g_x_0_x_0_zz_xyzzz, g_x_0_x_0_zz_xyzzzz, g_x_0_x_0_zz_xzzzz, g_x_0_x_0_zz_xzzzzz, g_x_0_x_0_zz_yyyyy, g_x_0_x_0_zz_yyyyyz, g_x_0_x_0_zz_yyyyz, g_x_0_x_0_zz_yyyyzz, g_x_0_x_0_zz_yyyzz, g_x_0_x_0_zz_yyyzzz, g_x_0_x_0_zz_yyzzz, g_x_0_x_0_zz_yyzzzz, g_x_0_x_0_zz_yzzzz, g_x_0_x_0_zz_yzzzzz, g_x_0_x_0_zz_zzzzz, g_x_0_x_0_zz_zzzzzz, g_x_0_x_0_zzz_xxxxx, g_x_0_x_0_zzz_xxxxy, g_x_0_x_0_zzz_xxxxz, g_x_0_x_0_zzz_xxxyy, g_x_0_x_0_zzz_xxxyz, g_x_0_x_0_zzz_xxxzz, g_x_0_x_0_zzz_xxyyy, g_x_0_x_0_zzz_xxyyz, g_x_0_x_0_zzz_xxyzz, g_x_0_x_0_zzz_xxzzz, g_x_0_x_0_zzz_xyyyy, g_x_0_x_0_zzz_xyyyz, g_x_0_x_0_zzz_xyyzz, g_x_0_x_0_zzz_xyzzz, g_x_0_x_0_zzz_xzzzz, g_x_0_x_0_zzz_yyyyy, g_x_0_x_0_zzz_yyyyz, g_x_0_x_0_zzz_yyyzz, g_x_0_x_0_zzz_yyzzz, g_x_0_x_0_zzz_yzzzz, g_x_0_x_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzz_xxxxx[k] = -g_x_0_x_0_zz_xxxxx[k] * ab_z + g_x_0_x_0_zz_xxxxxz[k];

                g_x_0_x_0_zzz_xxxxy[k] = -g_x_0_x_0_zz_xxxxy[k] * ab_z + g_x_0_x_0_zz_xxxxyz[k];

                g_x_0_x_0_zzz_xxxxz[k] = -g_x_0_x_0_zz_xxxxz[k] * ab_z + g_x_0_x_0_zz_xxxxzz[k];

                g_x_0_x_0_zzz_xxxyy[k] = -g_x_0_x_0_zz_xxxyy[k] * ab_z + g_x_0_x_0_zz_xxxyyz[k];

                g_x_0_x_0_zzz_xxxyz[k] = -g_x_0_x_0_zz_xxxyz[k] * ab_z + g_x_0_x_0_zz_xxxyzz[k];

                g_x_0_x_0_zzz_xxxzz[k] = -g_x_0_x_0_zz_xxxzz[k] * ab_z + g_x_0_x_0_zz_xxxzzz[k];

                g_x_0_x_0_zzz_xxyyy[k] = -g_x_0_x_0_zz_xxyyy[k] * ab_z + g_x_0_x_0_zz_xxyyyz[k];

                g_x_0_x_0_zzz_xxyyz[k] = -g_x_0_x_0_zz_xxyyz[k] * ab_z + g_x_0_x_0_zz_xxyyzz[k];

                g_x_0_x_0_zzz_xxyzz[k] = -g_x_0_x_0_zz_xxyzz[k] * ab_z + g_x_0_x_0_zz_xxyzzz[k];

                g_x_0_x_0_zzz_xxzzz[k] = -g_x_0_x_0_zz_xxzzz[k] * ab_z + g_x_0_x_0_zz_xxzzzz[k];

                g_x_0_x_0_zzz_xyyyy[k] = -g_x_0_x_0_zz_xyyyy[k] * ab_z + g_x_0_x_0_zz_xyyyyz[k];

                g_x_0_x_0_zzz_xyyyz[k] = -g_x_0_x_0_zz_xyyyz[k] * ab_z + g_x_0_x_0_zz_xyyyzz[k];

                g_x_0_x_0_zzz_xyyzz[k] = -g_x_0_x_0_zz_xyyzz[k] * ab_z + g_x_0_x_0_zz_xyyzzz[k];

                g_x_0_x_0_zzz_xyzzz[k] = -g_x_0_x_0_zz_xyzzz[k] * ab_z + g_x_0_x_0_zz_xyzzzz[k];

                g_x_0_x_0_zzz_xzzzz[k] = -g_x_0_x_0_zz_xzzzz[k] * ab_z + g_x_0_x_0_zz_xzzzzz[k];

                g_x_0_x_0_zzz_yyyyy[k] = -g_x_0_x_0_zz_yyyyy[k] * ab_z + g_x_0_x_0_zz_yyyyyz[k];

                g_x_0_x_0_zzz_yyyyz[k] = -g_x_0_x_0_zz_yyyyz[k] * ab_z + g_x_0_x_0_zz_yyyyzz[k];

                g_x_0_x_0_zzz_yyyzz[k] = -g_x_0_x_0_zz_yyyzz[k] * ab_z + g_x_0_x_0_zz_yyyzzz[k];

                g_x_0_x_0_zzz_yyzzz[k] = -g_x_0_x_0_zz_yyzzz[k] * ab_z + g_x_0_x_0_zz_yyzzzz[k];

                g_x_0_x_0_zzz_yzzzz[k] = -g_x_0_x_0_zz_yzzzz[k] * ab_z + g_x_0_x_0_zz_yzzzzz[k];

                g_x_0_x_0_zzz_zzzzz[k] = -g_x_0_x_0_zz_zzzzz[k] * ab_z + g_x_0_x_0_zz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 223 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xx_xxxxx, g_0_0_y_0_xx_xxxxy, g_0_0_y_0_xx_xxxxz, g_0_0_y_0_xx_xxxyy, g_0_0_y_0_xx_xxxyz, g_0_0_y_0_xx_xxxzz, g_0_0_y_0_xx_xxyyy, g_0_0_y_0_xx_xxyyz, g_0_0_y_0_xx_xxyzz, g_0_0_y_0_xx_xxzzz, g_0_0_y_0_xx_xyyyy, g_0_0_y_0_xx_xyyyz, g_0_0_y_0_xx_xyyzz, g_0_0_y_0_xx_xyzzz, g_0_0_y_0_xx_xzzzz, g_0_0_y_0_xx_yyyyy, g_0_0_y_0_xx_yyyyz, g_0_0_y_0_xx_yyyzz, g_0_0_y_0_xx_yyzzz, g_0_0_y_0_xx_yzzzz, g_0_0_y_0_xx_zzzzz, g_x_0_y_0_xx_xxxxx, g_x_0_y_0_xx_xxxxxx, g_x_0_y_0_xx_xxxxxy, g_x_0_y_0_xx_xxxxxz, g_x_0_y_0_xx_xxxxy, g_x_0_y_0_xx_xxxxyy, g_x_0_y_0_xx_xxxxyz, g_x_0_y_0_xx_xxxxz, g_x_0_y_0_xx_xxxxzz, g_x_0_y_0_xx_xxxyy, g_x_0_y_0_xx_xxxyyy, g_x_0_y_0_xx_xxxyyz, g_x_0_y_0_xx_xxxyz, g_x_0_y_0_xx_xxxyzz, g_x_0_y_0_xx_xxxzz, g_x_0_y_0_xx_xxxzzz, g_x_0_y_0_xx_xxyyy, g_x_0_y_0_xx_xxyyyy, g_x_0_y_0_xx_xxyyyz, g_x_0_y_0_xx_xxyyz, g_x_0_y_0_xx_xxyyzz, g_x_0_y_0_xx_xxyzz, g_x_0_y_0_xx_xxyzzz, g_x_0_y_0_xx_xxzzz, g_x_0_y_0_xx_xxzzzz, g_x_0_y_0_xx_xyyyy, g_x_0_y_0_xx_xyyyyy, g_x_0_y_0_xx_xyyyyz, g_x_0_y_0_xx_xyyyz, g_x_0_y_0_xx_xyyyzz, g_x_0_y_0_xx_xyyzz, g_x_0_y_0_xx_xyyzzz, g_x_0_y_0_xx_xyzzz, g_x_0_y_0_xx_xyzzzz, g_x_0_y_0_xx_xzzzz, g_x_0_y_0_xx_xzzzzz, g_x_0_y_0_xx_yyyyy, g_x_0_y_0_xx_yyyyz, g_x_0_y_0_xx_yyyzz, g_x_0_y_0_xx_yyzzz, g_x_0_y_0_xx_yzzzz, g_x_0_y_0_xx_zzzzz, g_x_0_y_0_xxx_xxxxx, g_x_0_y_0_xxx_xxxxy, g_x_0_y_0_xxx_xxxxz, g_x_0_y_0_xxx_xxxyy, g_x_0_y_0_xxx_xxxyz, g_x_0_y_0_xxx_xxxzz, g_x_0_y_0_xxx_xxyyy, g_x_0_y_0_xxx_xxyyz, g_x_0_y_0_xxx_xxyzz, g_x_0_y_0_xxx_xxzzz, g_x_0_y_0_xxx_xyyyy, g_x_0_y_0_xxx_xyyyz, g_x_0_y_0_xxx_xyyzz, g_x_0_y_0_xxx_xyzzz, g_x_0_y_0_xxx_xzzzz, g_x_0_y_0_xxx_yyyyy, g_x_0_y_0_xxx_yyyyz, g_x_0_y_0_xxx_yyyzz, g_x_0_y_0_xxx_yyzzz, g_x_0_y_0_xxx_yzzzz, g_x_0_y_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxx_xxxxx[k] = -g_0_0_y_0_xx_xxxxx[k] - g_x_0_y_0_xx_xxxxx[k] * ab_x + g_x_0_y_0_xx_xxxxxx[k];

                g_x_0_y_0_xxx_xxxxy[k] = -g_0_0_y_0_xx_xxxxy[k] - g_x_0_y_0_xx_xxxxy[k] * ab_x + g_x_0_y_0_xx_xxxxxy[k];

                g_x_0_y_0_xxx_xxxxz[k] = -g_0_0_y_0_xx_xxxxz[k] - g_x_0_y_0_xx_xxxxz[k] * ab_x + g_x_0_y_0_xx_xxxxxz[k];

                g_x_0_y_0_xxx_xxxyy[k] = -g_0_0_y_0_xx_xxxyy[k] - g_x_0_y_0_xx_xxxyy[k] * ab_x + g_x_0_y_0_xx_xxxxyy[k];

                g_x_0_y_0_xxx_xxxyz[k] = -g_0_0_y_0_xx_xxxyz[k] - g_x_0_y_0_xx_xxxyz[k] * ab_x + g_x_0_y_0_xx_xxxxyz[k];

                g_x_0_y_0_xxx_xxxzz[k] = -g_0_0_y_0_xx_xxxzz[k] - g_x_0_y_0_xx_xxxzz[k] * ab_x + g_x_0_y_0_xx_xxxxzz[k];

                g_x_0_y_0_xxx_xxyyy[k] = -g_0_0_y_0_xx_xxyyy[k] - g_x_0_y_0_xx_xxyyy[k] * ab_x + g_x_0_y_0_xx_xxxyyy[k];

                g_x_0_y_0_xxx_xxyyz[k] = -g_0_0_y_0_xx_xxyyz[k] - g_x_0_y_0_xx_xxyyz[k] * ab_x + g_x_0_y_0_xx_xxxyyz[k];

                g_x_0_y_0_xxx_xxyzz[k] = -g_0_0_y_0_xx_xxyzz[k] - g_x_0_y_0_xx_xxyzz[k] * ab_x + g_x_0_y_0_xx_xxxyzz[k];

                g_x_0_y_0_xxx_xxzzz[k] = -g_0_0_y_0_xx_xxzzz[k] - g_x_0_y_0_xx_xxzzz[k] * ab_x + g_x_0_y_0_xx_xxxzzz[k];

                g_x_0_y_0_xxx_xyyyy[k] = -g_0_0_y_0_xx_xyyyy[k] - g_x_0_y_0_xx_xyyyy[k] * ab_x + g_x_0_y_0_xx_xxyyyy[k];

                g_x_0_y_0_xxx_xyyyz[k] = -g_0_0_y_0_xx_xyyyz[k] - g_x_0_y_0_xx_xyyyz[k] * ab_x + g_x_0_y_0_xx_xxyyyz[k];

                g_x_0_y_0_xxx_xyyzz[k] = -g_0_0_y_0_xx_xyyzz[k] - g_x_0_y_0_xx_xyyzz[k] * ab_x + g_x_0_y_0_xx_xxyyzz[k];

                g_x_0_y_0_xxx_xyzzz[k] = -g_0_0_y_0_xx_xyzzz[k] - g_x_0_y_0_xx_xyzzz[k] * ab_x + g_x_0_y_0_xx_xxyzzz[k];

                g_x_0_y_0_xxx_xzzzz[k] = -g_0_0_y_0_xx_xzzzz[k] - g_x_0_y_0_xx_xzzzz[k] * ab_x + g_x_0_y_0_xx_xxzzzz[k];

                g_x_0_y_0_xxx_yyyyy[k] = -g_0_0_y_0_xx_yyyyy[k] - g_x_0_y_0_xx_yyyyy[k] * ab_x + g_x_0_y_0_xx_xyyyyy[k];

                g_x_0_y_0_xxx_yyyyz[k] = -g_0_0_y_0_xx_yyyyz[k] - g_x_0_y_0_xx_yyyyz[k] * ab_x + g_x_0_y_0_xx_xyyyyz[k];

                g_x_0_y_0_xxx_yyyzz[k] = -g_0_0_y_0_xx_yyyzz[k] - g_x_0_y_0_xx_yyyzz[k] * ab_x + g_x_0_y_0_xx_xyyyzz[k];

                g_x_0_y_0_xxx_yyzzz[k] = -g_0_0_y_0_xx_yyzzz[k] - g_x_0_y_0_xx_yyzzz[k] * ab_x + g_x_0_y_0_xx_xyyzzz[k];

                g_x_0_y_0_xxx_yzzzz[k] = -g_0_0_y_0_xx_yzzzz[k] - g_x_0_y_0_xx_yzzzz[k] * ab_x + g_x_0_y_0_xx_xyzzzz[k];

                g_x_0_y_0_xxx_zzzzz[k] = -g_0_0_y_0_xx_zzzzz[k] - g_x_0_y_0_xx_zzzzz[k] * ab_x + g_x_0_y_0_xx_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xx_xxxxx, g_x_0_y_0_xx_xxxxxy, g_x_0_y_0_xx_xxxxy, g_x_0_y_0_xx_xxxxyy, g_x_0_y_0_xx_xxxxyz, g_x_0_y_0_xx_xxxxz, g_x_0_y_0_xx_xxxyy, g_x_0_y_0_xx_xxxyyy, g_x_0_y_0_xx_xxxyyz, g_x_0_y_0_xx_xxxyz, g_x_0_y_0_xx_xxxyzz, g_x_0_y_0_xx_xxxzz, g_x_0_y_0_xx_xxyyy, g_x_0_y_0_xx_xxyyyy, g_x_0_y_0_xx_xxyyyz, g_x_0_y_0_xx_xxyyz, g_x_0_y_0_xx_xxyyzz, g_x_0_y_0_xx_xxyzz, g_x_0_y_0_xx_xxyzzz, g_x_0_y_0_xx_xxzzz, g_x_0_y_0_xx_xyyyy, g_x_0_y_0_xx_xyyyyy, g_x_0_y_0_xx_xyyyyz, g_x_0_y_0_xx_xyyyz, g_x_0_y_0_xx_xyyyzz, g_x_0_y_0_xx_xyyzz, g_x_0_y_0_xx_xyyzzz, g_x_0_y_0_xx_xyzzz, g_x_0_y_0_xx_xyzzzz, g_x_0_y_0_xx_xzzzz, g_x_0_y_0_xx_yyyyy, g_x_0_y_0_xx_yyyyyy, g_x_0_y_0_xx_yyyyyz, g_x_0_y_0_xx_yyyyz, g_x_0_y_0_xx_yyyyzz, g_x_0_y_0_xx_yyyzz, g_x_0_y_0_xx_yyyzzz, g_x_0_y_0_xx_yyzzz, g_x_0_y_0_xx_yyzzzz, g_x_0_y_0_xx_yzzzz, g_x_0_y_0_xx_yzzzzz, g_x_0_y_0_xx_zzzzz, g_x_0_y_0_xxy_xxxxx, g_x_0_y_0_xxy_xxxxy, g_x_0_y_0_xxy_xxxxz, g_x_0_y_0_xxy_xxxyy, g_x_0_y_0_xxy_xxxyz, g_x_0_y_0_xxy_xxxzz, g_x_0_y_0_xxy_xxyyy, g_x_0_y_0_xxy_xxyyz, g_x_0_y_0_xxy_xxyzz, g_x_0_y_0_xxy_xxzzz, g_x_0_y_0_xxy_xyyyy, g_x_0_y_0_xxy_xyyyz, g_x_0_y_0_xxy_xyyzz, g_x_0_y_0_xxy_xyzzz, g_x_0_y_0_xxy_xzzzz, g_x_0_y_0_xxy_yyyyy, g_x_0_y_0_xxy_yyyyz, g_x_0_y_0_xxy_yyyzz, g_x_0_y_0_xxy_yyzzz, g_x_0_y_0_xxy_yzzzz, g_x_0_y_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxy_xxxxx[k] = -g_x_0_y_0_xx_xxxxx[k] * ab_y + g_x_0_y_0_xx_xxxxxy[k];

                g_x_0_y_0_xxy_xxxxy[k] = -g_x_0_y_0_xx_xxxxy[k] * ab_y + g_x_0_y_0_xx_xxxxyy[k];

                g_x_0_y_0_xxy_xxxxz[k] = -g_x_0_y_0_xx_xxxxz[k] * ab_y + g_x_0_y_0_xx_xxxxyz[k];

                g_x_0_y_0_xxy_xxxyy[k] = -g_x_0_y_0_xx_xxxyy[k] * ab_y + g_x_0_y_0_xx_xxxyyy[k];

                g_x_0_y_0_xxy_xxxyz[k] = -g_x_0_y_0_xx_xxxyz[k] * ab_y + g_x_0_y_0_xx_xxxyyz[k];

                g_x_0_y_0_xxy_xxxzz[k] = -g_x_0_y_0_xx_xxxzz[k] * ab_y + g_x_0_y_0_xx_xxxyzz[k];

                g_x_0_y_0_xxy_xxyyy[k] = -g_x_0_y_0_xx_xxyyy[k] * ab_y + g_x_0_y_0_xx_xxyyyy[k];

                g_x_0_y_0_xxy_xxyyz[k] = -g_x_0_y_0_xx_xxyyz[k] * ab_y + g_x_0_y_0_xx_xxyyyz[k];

                g_x_0_y_0_xxy_xxyzz[k] = -g_x_0_y_0_xx_xxyzz[k] * ab_y + g_x_0_y_0_xx_xxyyzz[k];

                g_x_0_y_0_xxy_xxzzz[k] = -g_x_0_y_0_xx_xxzzz[k] * ab_y + g_x_0_y_0_xx_xxyzzz[k];

                g_x_0_y_0_xxy_xyyyy[k] = -g_x_0_y_0_xx_xyyyy[k] * ab_y + g_x_0_y_0_xx_xyyyyy[k];

                g_x_0_y_0_xxy_xyyyz[k] = -g_x_0_y_0_xx_xyyyz[k] * ab_y + g_x_0_y_0_xx_xyyyyz[k];

                g_x_0_y_0_xxy_xyyzz[k] = -g_x_0_y_0_xx_xyyzz[k] * ab_y + g_x_0_y_0_xx_xyyyzz[k];

                g_x_0_y_0_xxy_xyzzz[k] = -g_x_0_y_0_xx_xyzzz[k] * ab_y + g_x_0_y_0_xx_xyyzzz[k];

                g_x_0_y_0_xxy_xzzzz[k] = -g_x_0_y_0_xx_xzzzz[k] * ab_y + g_x_0_y_0_xx_xyzzzz[k];

                g_x_0_y_0_xxy_yyyyy[k] = -g_x_0_y_0_xx_yyyyy[k] * ab_y + g_x_0_y_0_xx_yyyyyy[k];

                g_x_0_y_0_xxy_yyyyz[k] = -g_x_0_y_0_xx_yyyyz[k] * ab_y + g_x_0_y_0_xx_yyyyyz[k];

                g_x_0_y_0_xxy_yyyzz[k] = -g_x_0_y_0_xx_yyyzz[k] * ab_y + g_x_0_y_0_xx_yyyyzz[k];

                g_x_0_y_0_xxy_yyzzz[k] = -g_x_0_y_0_xx_yyzzz[k] * ab_y + g_x_0_y_0_xx_yyyzzz[k];

                g_x_0_y_0_xxy_yzzzz[k] = -g_x_0_y_0_xx_yzzzz[k] * ab_y + g_x_0_y_0_xx_yyzzzz[k];

                g_x_0_y_0_xxy_zzzzz[k] = -g_x_0_y_0_xx_zzzzz[k] * ab_y + g_x_0_y_0_xx_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 252 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 253 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 254 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 255 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 256 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 257 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 258 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 259 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 260 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 261 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 262 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 263 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 264 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 265 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 266 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 267 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 268 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 269 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 270 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 271 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xx_xxxxx, g_x_0_y_0_xx_xxxxxz, g_x_0_y_0_xx_xxxxy, g_x_0_y_0_xx_xxxxyz, g_x_0_y_0_xx_xxxxz, g_x_0_y_0_xx_xxxxzz, g_x_0_y_0_xx_xxxyy, g_x_0_y_0_xx_xxxyyz, g_x_0_y_0_xx_xxxyz, g_x_0_y_0_xx_xxxyzz, g_x_0_y_0_xx_xxxzz, g_x_0_y_0_xx_xxxzzz, g_x_0_y_0_xx_xxyyy, g_x_0_y_0_xx_xxyyyz, g_x_0_y_0_xx_xxyyz, g_x_0_y_0_xx_xxyyzz, g_x_0_y_0_xx_xxyzz, g_x_0_y_0_xx_xxyzzz, g_x_0_y_0_xx_xxzzz, g_x_0_y_0_xx_xxzzzz, g_x_0_y_0_xx_xyyyy, g_x_0_y_0_xx_xyyyyz, g_x_0_y_0_xx_xyyyz, g_x_0_y_0_xx_xyyyzz, g_x_0_y_0_xx_xyyzz, g_x_0_y_0_xx_xyyzzz, g_x_0_y_0_xx_xyzzz, g_x_0_y_0_xx_xyzzzz, g_x_0_y_0_xx_xzzzz, g_x_0_y_0_xx_xzzzzz, g_x_0_y_0_xx_yyyyy, g_x_0_y_0_xx_yyyyyz, g_x_0_y_0_xx_yyyyz, g_x_0_y_0_xx_yyyyzz, g_x_0_y_0_xx_yyyzz, g_x_0_y_0_xx_yyyzzz, g_x_0_y_0_xx_yyzzz, g_x_0_y_0_xx_yyzzzz, g_x_0_y_0_xx_yzzzz, g_x_0_y_0_xx_yzzzzz, g_x_0_y_0_xx_zzzzz, g_x_0_y_0_xx_zzzzzz, g_x_0_y_0_xxz_xxxxx, g_x_0_y_0_xxz_xxxxy, g_x_0_y_0_xxz_xxxxz, g_x_0_y_0_xxz_xxxyy, g_x_0_y_0_xxz_xxxyz, g_x_0_y_0_xxz_xxxzz, g_x_0_y_0_xxz_xxyyy, g_x_0_y_0_xxz_xxyyz, g_x_0_y_0_xxz_xxyzz, g_x_0_y_0_xxz_xxzzz, g_x_0_y_0_xxz_xyyyy, g_x_0_y_0_xxz_xyyyz, g_x_0_y_0_xxz_xyyzz, g_x_0_y_0_xxz_xyzzz, g_x_0_y_0_xxz_xzzzz, g_x_0_y_0_xxz_yyyyy, g_x_0_y_0_xxz_yyyyz, g_x_0_y_0_xxz_yyyzz, g_x_0_y_0_xxz_yyzzz, g_x_0_y_0_xxz_yzzzz, g_x_0_y_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxz_xxxxx[k] = -g_x_0_y_0_xx_xxxxx[k] * ab_z + g_x_0_y_0_xx_xxxxxz[k];

                g_x_0_y_0_xxz_xxxxy[k] = -g_x_0_y_0_xx_xxxxy[k] * ab_z + g_x_0_y_0_xx_xxxxyz[k];

                g_x_0_y_0_xxz_xxxxz[k] = -g_x_0_y_0_xx_xxxxz[k] * ab_z + g_x_0_y_0_xx_xxxxzz[k];

                g_x_0_y_0_xxz_xxxyy[k] = -g_x_0_y_0_xx_xxxyy[k] * ab_z + g_x_0_y_0_xx_xxxyyz[k];

                g_x_0_y_0_xxz_xxxyz[k] = -g_x_0_y_0_xx_xxxyz[k] * ab_z + g_x_0_y_0_xx_xxxyzz[k];

                g_x_0_y_0_xxz_xxxzz[k] = -g_x_0_y_0_xx_xxxzz[k] * ab_z + g_x_0_y_0_xx_xxxzzz[k];

                g_x_0_y_0_xxz_xxyyy[k] = -g_x_0_y_0_xx_xxyyy[k] * ab_z + g_x_0_y_0_xx_xxyyyz[k];

                g_x_0_y_0_xxz_xxyyz[k] = -g_x_0_y_0_xx_xxyyz[k] * ab_z + g_x_0_y_0_xx_xxyyzz[k];

                g_x_0_y_0_xxz_xxyzz[k] = -g_x_0_y_0_xx_xxyzz[k] * ab_z + g_x_0_y_0_xx_xxyzzz[k];

                g_x_0_y_0_xxz_xxzzz[k] = -g_x_0_y_0_xx_xxzzz[k] * ab_z + g_x_0_y_0_xx_xxzzzz[k];

                g_x_0_y_0_xxz_xyyyy[k] = -g_x_0_y_0_xx_xyyyy[k] * ab_z + g_x_0_y_0_xx_xyyyyz[k];

                g_x_0_y_0_xxz_xyyyz[k] = -g_x_0_y_0_xx_xyyyz[k] * ab_z + g_x_0_y_0_xx_xyyyzz[k];

                g_x_0_y_0_xxz_xyyzz[k] = -g_x_0_y_0_xx_xyyzz[k] * ab_z + g_x_0_y_0_xx_xyyzzz[k];

                g_x_0_y_0_xxz_xyzzz[k] = -g_x_0_y_0_xx_xyzzz[k] * ab_z + g_x_0_y_0_xx_xyzzzz[k];

                g_x_0_y_0_xxz_xzzzz[k] = -g_x_0_y_0_xx_xzzzz[k] * ab_z + g_x_0_y_0_xx_xzzzzz[k];

                g_x_0_y_0_xxz_yyyyy[k] = -g_x_0_y_0_xx_yyyyy[k] * ab_z + g_x_0_y_0_xx_yyyyyz[k];

                g_x_0_y_0_xxz_yyyyz[k] = -g_x_0_y_0_xx_yyyyz[k] * ab_z + g_x_0_y_0_xx_yyyyzz[k];

                g_x_0_y_0_xxz_yyyzz[k] = -g_x_0_y_0_xx_yyyzz[k] * ab_z + g_x_0_y_0_xx_yyyzzz[k];

                g_x_0_y_0_xxz_yyzzz[k] = -g_x_0_y_0_xx_yyzzz[k] * ab_z + g_x_0_y_0_xx_yyzzzz[k];

                g_x_0_y_0_xxz_yzzzz[k] = -g_x_0_y_0_xx_yzzzz[k] * ab_z + g_x_0_y_0_xx_yzzzzz[k];

                g_x_0_y_0_xxz_zzzzz[k] = -g_x_0_y_0_xx_zzzzz[k] * ab_z + g_x_0_y_0_xx_zzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 273 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 274 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 275 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 276 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 277 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 278 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 279 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 280 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 281 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 282 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 283 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 284 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 285 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 286 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 287 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 288 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 289 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 290 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 291 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 292 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xy_xxxxx, g_x_0_y_0_xy_xxxxxy, g_x_0_y_0_xy_xxxxy, g_x_0_y_0_xy_xxxxyy, g_x_0_y_0_xy_xxxxyz, g_x_0_y_0_xy_xxxxz, g_x_0_y_0_xy_xxxyy, g_x_0_y_0_xy_xxxyyy, g_x_0_y_0_xy_xxxyyz, g_x_0_y_0_xy_xxxyz, g_x_0_y_0_xy_xxxyzz, g_x_0_y_0_xy_xxxzz, g_x_0_y_0_xy_xxyyy, g_x_0_y_0_xy_xxyyyy, g_x_0_y_0_xy_xxyyyz, g_x_0_y_0_xy_xxyyz, g_x_0_y_0_xy_xxyyzz, g_x_0_y_0_xy_xxyzz, g_x_0_y_0_xy_xxyzzz, g_x_0_y_0_xy_xxzzz, g_x_0_y_0_xy_xyyyy, g_x_0_y_0_xy_xyyyyy, g_x_0_y_0_xy_xyyyyz, g_x_0_y_0_xy_xyyyz, g_x_0_y_0_xy_xyyyzz, g_x_0_y_0_xy_xyyzz, g_x_0_y_0_xy_xyyzzz, g_x_0_y_0_xy_xyzzz, g_x_0_y_0_xy_xyzzzz, g_x_0_y_0_xy_xzzzz, g_x_0_y_0_xy_yyyyy, g_x_0_y_0_xy_yyyyyy, g_x_0_y_0_xy_yyyyyz, g_x_0_y_0_xy_yyyyz, g_x_0_y_0_xy_yyyyzz, g_x_0_y_0_xy_yyyzz, g_x_0_y_0_xy_yyyzzz, g_x_0_y_0_xy_yyzzz, g_x_0_y_0_xy_yyzzzz, g_x_0_y_0_xy_yzzzz, g_x_0_y_0_xy_yzzzzz, g_x_0_y_0_xy_zzzzz, g_x_0_y_0_xyy_xxxxx, g_x_0_y_0_xyy_xxxxy, g_x_0_y_0_xyy_xxxxz, g_x_0_y_0_xyy_xxxyy, g_x_0_y_0_xyy_xxxyz, g_x_0_y_0_xyy_xxxzz, g_x_0_y_0_xyy_xxyyy, g_x_0_y_0_xyy_xxyyz, g_x_0_y_0_xyy_xxyzz, g_x_0_y_0_xyy_xxzzz, g_x_0_y_0_xyy_xyyyy, g_x_0_y_0_xyy_xyyyz, g_x_0_y_0_xyy_xyyzz, g_x_0_y_0_xyy_xyzzz, g_x_0_y_0_xyy_xzzzz, g_x_0_y_0_xyy_yyyyy, g_x_0_y_0_xyy_yyyyz, g_x_0_y_0_xyy_yyyzz, g_x_0_y_0_xyy_yyzzz, g_x_0_y_0_xyy_yzzzz, g_x_0_y_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyy_xxxxx[k] = -g_x_0_y_0_xy_xxxxx[k] * ab_y + g_x_0_y_0_xy_xxxxxy[k];

                g_x_0_y_0_xyy_xxxxy[k] = -g_x_0_y_0_xy_xxxxy[k] * ab_y + g_x_0_y_0_xy_xxxxyy[k];

                g_x_0_y_0_xyy_xxxxz[k] = -g_x_0_y_0_xy_xxxxz[k] * ab_y + g_x_0_y_0_xy_xxxxyz[k];

                g_x_0_y_0_xyy_xxxyy[k] = -g_x_0_y_0_xy_xxxyy[k] * ab_y + g_x_0_y_0_xy_xxxyyy[k];

                g_x_0_y_0_xyy_xxxyz[k] = -g_x_0_y_0_xy_xxxyz[k] * ab_y + g_x_0_y_0_xy_xxxyyz[k];

                g_x_0_y_0_xyy_xxxzz[k] = -g_x_0_y_0_xy_xxxzz[k] * ab_y + g_x_0_y_0_xy_xxxyzz[k];

                g_x_0_y_0_xyy_xxyyy[k] = -g_x_0_y_0_xy_xxyyy[k] * ab_y + g_x_0_y_0_xy_xxyyyy[k];

                g_x_0_y_0_xyy_xxyyz[k] = -g_x_0_y_0_xy_xxyyz[k] * ab_y + g_x_0_y_0_xy_xxyyyz[k];

                g_x_0_y_0_xyy_xxyzz[k] = -g_x_0_y_0_xy_xxyzz[k] * ab_y + g_x_0_y_0_xy_xxyyzz[k];

                g_x_0_y_0_xyy_xxzzz[k] = -g_x_0_y_0_xy_xxzzz[k] * ab_y + g_x_0_y_0_xy_xxyzzz[k];

                g_x_0_y_0_xyy_xyyyy[k] = -g_x_0_y_0_xy_xyyyy[k] * ab_y + g_x_0_y_0_xy_xyyyyy[k];

                g_x_0_y_0_xyy_xyyyz[k] = -g_x_0_y_0_xy_xyyyz[k] * ab_y + g_x_0_y_0_xy_xyyyyz[k];

                g_x_0_y_0_xyy_xyyzz[k] = -g_x_0_y_0_xy_xyyzz[k] * ab_y + g_x_0_y_0_xy_xyyyzz[k];

                g_x_0_y_0_xyy_xyzzz[k] = -g_x_0_y_0_xy_xyzzz[k] * ab_y + g_x_0_y_0_xy_xyyzzz[k];

                g_x_0_y_0_xyy_xzzzz[k] = -g_x_0_y_0_xy_xzzzz[k] * ab_y + g_x_0_y_0_xy_xyzzzz[k];

                g_x_0_y_0_xyy_yyyyy[k] = -g_x_0_y_0_xy_yyyyy[k] * ab_y + g_x_0_y_0_xy_yyyyyy[k];

                g_x_0_y_0_xyy_yyyyz[k] = -g_x_0_y_0_xy_yyyyz[k] * ab_y + g_x_0_y_0_xy_yyyyyz[k];

                g_x_0_y_0_xyy_yyyzz[k] = -g_x_0_y_0_xy_yyyzz[k] * ab_y + g_x_0_y_0_xy_yyyyzz[k];

                g_x_0_y_0_xyy_yyzzz[k] = -g_x_0_y_0_xy_yyzzz[k] * ab_y + g_x_0_y_0_xy_yyyzzz[k];

                g_x_0_y_0_xyy_yzzzz[k] = -g_x_0_y_0_xy_yzzzz[k] * ab_y + g_x_0_y_0_xy_yyzzzz[k];

                g_x_0_y_0_xyy_zzzzz[k] = -g_x_0_y_0_xy_zzzzz[k] * ab_y + g_x_0_y_0_xy_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 294 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 295 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 296 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 297 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 298 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 299 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 300 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 301 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 302 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 303 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 304 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 305 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 306 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 307 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 308 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 309 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 310 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 311 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 312 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 313 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyz_xxxxx, g_x_0_y_0_xyz_xxxxy, g_x_0_y_0_xyz_xxxxz, g_x_0_y_0_xyz_xxxyy, g_x_0_y_0_xyz_xxxyz, g_x_0_y_0_xyz_xxxzz, g_x_0_y_0_xyz_xxyyy, g_x_0_y_0_xyz_xxyyz, g_x_0_y_0_xyz_xxyzz, g_x_0_y_0_xyz_xxzzz, g_x_0_y_0_xyz_xyyyy, g_x_0_y_0_xyz_xyyyz, g_x_0_y_0_xyz_xyyzz, g_x_0_y_0_xyz_xyzzz, g_x_0_y_0_xyz_xzzzz, g_x_0_y_0_xyz_yyyyy, g_x_0_y_0_xyz_yyyyz, g_x_0_y_0_xyz_yyyzz, g_x_0_y_0_xyz_yyzzz, g_x_0_y_0_xyz_yzzzz, g_x_0_y_0_xyz_zzzzz, g_x_0_y_0_xz_xxxxx, g_x_0_y_0_xz_xxxxxy, g_x_0_y_0_xz_xxxxy, g_x_0_y_0_xz_xxxxyy, g_x_0_y_0_xz_xxxxyz, g_x_0_y_0_xz_xxxxz, g_x_0_y_0_xz_xxxyy, g_x_0_y_0_xz_xxxyyy, g_x_0_y_0_xz_xxxyyz, g_x_0_y_0_xz_xxxyz, g_x_0_y_0_xz_xxxyzz, g_x_0_y_0_xz_xxxzz, g_x_0_y_0_xz_xxyyy, g_x_0_y_0_xz_xxyyyy, g_x_0_y_0_xz_xxyyyz, g_x_0_y_0_xz_xxyyz, g_x_0_y_0_xz_xxyyzz, g_x_0_y_0_xz_xxyzz, g_x_0_y_0_xz_xxyzzz, g_x_0_y_0_xz_xxzzz, g_x_0_y_0_xz_xyyyy, g_x_0_y_0_xz_xyyyyy, g_x_0_y_0_xz_xyyyyz, g_x_0_y_0_xz_xyyyz, g_x_0_y_0_xz_xyyyzz, g_x_0_y_0_xz_xyyzz, g_x_0_y_0_xz_xyyzzz, g_x_0_y_0_xz_xyzzz, g_x_0_y_0_xz_xyzzzz, g_x_0_y_0_xz_xzzzz, g_x_0_y_0_xz_yyyyy, g_x_0_y_0_xz_yyyyyy, g_x_0_y_0_xz_yyyyyz, g_x_0_y_0_xz_yyyyz, g_x_0_y_0_xz_yyyyzz, g_x_0_y_0_xz_yyyzz, g_x_0_y_0_xz_yyyzzz, g_x_0_y_0_xz_yyzzz, g_x_0_y_0_xz_yyzzzz, g_x_0_y_0_xz_yzzzz, g_x_0_y_0_xz_yzzzzz, g_x_0_y_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyz_xxxxx[k] = -g_x_0_y_0_xz_xxxxx[k] * ab_y + g_x_0_y_0_xz_xxxxxy[k];

                g_x_0_y_0_xyz_xxxxy[k] = -g_x_0_y_0_xz_xxxxy[k] * ab_y + g_x_0_y_0_xz_xxxxyy[k];

                g_x_0_y_0_xyz_xxxxz[k] = -g_x_0_y_0_xz_xxxxz[k] * ab_y + g_x_0_y_0_xz_xxxxyz[k];

                g_x_0_y_0_xyz_xxxyy[k] = -g_x_0_y_0_xz_xxxyy[k] * ab_y + g_x_0_y_0_xz_xxxyyy[k];

                g_x_0_y_0_xyz_xxxyz[k] = -g_x_0_y_0_xz_xxxyz[k] * ab_y + g_x_0_y_0_xz_xxxyyz[k];

                g_x_0_y_0_xyz_xxxzz[k] = -g_x_0_y_0_xz_xxxzz[k] * ab_y + g_x_0_y_0_xz_xxxyzz[k];

                g_x_0_y_0_xyz_xxyyy[k] = -g_x_0_y_0_xz_xxyyy[k] * ab_y + g_x_0_y_0_xz_xxyyyy[k];

                g_x_0_y_0_xyz_xxyyz[k] = -g_x_0_y_0_xz_xxyyz[k] * ab_y + g_x_0_y_0_xz_xxyyyz[k];

                g_x_0_y_0_xyz_xxyzz[k] = -g_x_0_y_0_xz_xxyzz[k] * ab_y + g_x_0_y_0_xz_xxyyzz[k];

                g_x_0_y_0_xyz_xxzzz[k] = -g_x_0_y_0_xz_xxzzz[k] * ab_y + g_x_0_y_0_xz_xxyzzz[k];

                g_x_0_y_0_xyz_xyyyy[k] = -g_x_0_y_0_xz_xyyyy[k] * ab_y + g_x_0_y_0_xz_xyyyyy[k];

                g_x_0_y_0_xyz_xyyyz[k] = -g_x_0_y_0_xz_xyyyz[k] * ab_y + g_x_0_y_0_xz_xyyyyz[k];

                g_x_0_y_0_xyz_xyyzz[k] = -g_x_0_y_0_xz_xyyzz[k] * ab_y + g_x_0_y_0_xz_xyyyzz[k];

                g_x_0_y_0_xyz_xyzzz[k] = -g_x_0_y_0_xz_xyzzz[k] * ab_y + g_x_0_y_0_xz_xyyzzz[k];

                g_x_0_y_0_xyz_xzzzz[k] = -g_x_0_y_0_xz_xzzzz[k] * ab_y + g_x_0_y_0_xz_xyzzzz[k];

                g_x_0_y_0_xyz_yyyyy[k] = -g_x_0_y_0_xz_yyyyy[k] * ab_y + g_x_0_y_0_xz_yyyyyy[k];

                g_x_0_y_0_xyz_yyyyz[k] = -g_x_0_y_0_xz_yyyyz[k] * ab_y + g_x_0_y_0_xz_yyyyyz[k];

                g_x_0_y_0_xyz_yyyzz[k] = -g_x_0_y_0_xz_yyyzz[k] * ab_y + g_x_0_y_0_xz_yyyyzz[k];

                g_x_0_y_0_xyz_yyzzz[k] = -g_x_0_y_0_xz_yyzzz[k] * ab_y + g_x_0_y_0_xz_yyyzzz[k];

                g_x_0_y_0_xyz_yzzzz[k] = -g_x_0_y_0_xz_yzzzz[k] * ab_y + g_x_0_y_0_xz_yyzzzz[k];

                g_x_0_y_0_xyz_zzzzz[k] = -g_x_0_y_0_xz_zzzzz[k] * ab_y + g_x_0_y_0_xz_yzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 315 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 316 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 317 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 318 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 319 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 320 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 321 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 322 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 323 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 324 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 325 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 326 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 327 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 328 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 329 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 330 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 331 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 332 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 333 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 334 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xz_xxxxx, g_x_0_y_0_xz_xxxxxz, g_x_0_y_0_xz_xxxxy, g_x_0_y_0_xz_xxxxyz, g_x_0_y_0_xz_xxxxz, g_x_0_y_0_xz_xxxxzz, g_x_0_y_0_xz_xxxyy, g_x_0_y_0_xz_xxxyyz, g_x_0_y_0_xz_xxxyz, g_x_0_y_0_xz_xxxyzz, g_x_0_y_0_xz_xxxzz, g_x_0_y_0_xz_xxxzzz, g_x_0_y_0_xz_xxyyy, g_x_0_y_0_xz_xxyyyz, g_x_0_y_0_xz_xxyyz, g_x_0_y_0_xz_xxyyzz, g_x_0_y_0_xz_xxyzz, g_x_0_y_0_xz_xxyzzz, g_x_0_y_0_xz_xxzzz, g_x_0_y_0_xz_xxzzzz, g_x_0_y_0_xz_xyyyy, g_x_0_y_0_xz_xyyyyz, g_x_0_y_0_xz_xyyyz, g_x_0_y_0_xz_xyyyzz, g_x_0_y_0_xz_xyyzz, g_x_0_y_0_xz_xyyzzz, g_x_0_y_0_xz_xyzzz, g_x_0_y_0_xz_xyzzzz, g_x_0_y_0_xz_xzzzz, g_x_0_y_0_xz_xzzzzz, g_x_0_y_0_xz_yyyyy, g_x_0_y_0_xz_yyyyyz, g_x_0_y_0_xz_yyyyz, g_x_0_y_0_xz_yyyyzz, g_x_0_y_0_xz_yyyzz, g_x_0_y_0_xz_yyyzzz, g_x_0_y_0_xz_yyzzz, g_x_0_y_0_xz_yyzzzz, g_x_0_y_0_xz_yzzzz, g_x_0_y_0_xz_yzzzzz, g_x_0_y_0_xz_zzzzz, g_x_0_y_0_xz_zzzzzz, g_x_0_y_0_xzz_xxxxx, g_x_0_y_0_xzz_xxxxy, g_x_0_y_0_xzz_xxxxz, g_x_0_y_0_xzz_xxxyy, g_x_0_y_0_xzz_xxxyz, g_x_0_y_0_xzz_xxxzz, g_x_0_y_0_xzz_xxyyy, g_x_0_y_0_xzz_xxyyz, g_x_0_y_0_xzz_xxyzz, g_x_0_y_0_xzz_xxzzz, g_x_0_y_0_xzz_xyyyy, g_x_0_y_0_xzz_xyyyz, g_x_0_y_0_xzz_xyyzz, g_x_0_y_0_xzz_xyzzz, g_x_0_y_0_xzz_xzzzz, g_x_0_y_0_xzz_yyyyy, g_x_0_y_0_xzz_yyyyz, g_x_0_y_0_xzz_yyyzz, g_x_0_y_0_xzz_yyzzz, g_x_0_y_0_xzz_yzzzz, g_x_0_y_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzz_xxxxx[k] = -g_x_0_y_0_xz_xxxxx[k] * ab_z + g_x_0_y_0_xz_xxxxxz[k];

                g_x_0_y_0_xzz_xxxxy[k] = -g_x_0_y_0_xz_xxxxy[k] * ab_z + g_x_0_y_0_xz_xxxxyz[k];

                g_x_0_y_0_xzz_xxxxz[k] = -g_x_0_y_0_xz_xxxxz[k] * ab_z + g_x_0_y_0_xz_xxxxzz[k];

                g_x_0_y_0_xzz_xxxyy[k] = -g_x_0_y_0_xz_xxxyy[k] * ab_z + g_x_0_y_0_xz_xxxyyz[k];

                g_x_0_y_0_xzz_xxxyz[k] = -g_x_0_y_0_xz_xxxyz[k] * ab_z + g_x_0_y_0_xz_xxxyzz[k];

                g_x_0_y_0_xzz_xxxzz[k] = -g_x_0_y_0_xz_xxxzz[k] * ab_z + g_x_0_y_0_xz_xxxzzz[k];

                g_x_0_y_0_xzz_xxyyy[k] = -g_x_0_y_0_xz_xxyyy[k] * ab_z + g_x_0_y_0_xz_xxyyyz[k];

                g_x_0_y_0_xzz_xxyyz[k] = -g_x_0_y_0_xz_xxyyz[k] * ab_z + g_x_0_y_0_xz_xxyyzz[k];

                g_x_0_y_0_xzz_xxyzz[k] = -g_x_0_y_0_xz_xxyzz[k] * ab_z + g_x_0_y_0_xz_xxyzzz[k];

                g_x_0_y_0_xzz_xxzzz[k] = -g_x_0_y_0_xz_xxzzz[k] * ab_z + g_x_0_y_0_xz_xxzzzz[k];

                g_x_0_y_0_xzz_xyyyy[k] = -g_x_0_y_0_xz_xyyyy[k] * ab_z + g_x_0_y_0_xz_xyyyyz[k];

                g_x_0_y_0_xzz_xyyyz[k] = -g_x_0_y_0_xz_xyyyz[k] * ab_z + g_x_0_y_0_xz_xyyyzz[k];

                g_x_0_y_0_xzz_xyyzz[k] = -g_x_0_y_0_xz_xyyzz[k] * ab_z + g_x_0_y_0_xz_xyyzzz[k];

                g_x_0_y_0_xzz_xyzzz[k] = -g_x_0_y_0_xz_xyzzz[k] * ab_z + g_x_0_y_0_xz_xyzzzz[k];

                g_x_0_y_0_xzz_xzzzz[k] = -g_x_0_y_0_xz_xzzzz[k] * ab_z + g_x_0_y_0_xz_xzzzzz[k];

                g_x_0_y_0_xzz_yyyyy[k] = -g_x_0_y_0_xz_yyyyy[k] * ab_z + g_x_0_y_0_xz_yyyyyz[k];

                g_x_0_y_0_xzz_yyyyz[k] = -g_x_0_y_0_xz_yyyyz[k] * ab_z + g_x_0_y_0_xz_yyyyzz[k];

                g_x_0_y_0_xzz_yyyzz[k] = -g_x_0_y_0_xz_yyyzz[k] * ab_z + g_x_0_y_0_xz_yyyzzz[k];

                g_x_0_y_0_xzz_yyzzz[k] = -g_x_0_y_0_xz_yyzzz[k] * ab_z + g_x_0_y_0_xz_yyzzzz[k];

                g_x_0_y_0_xzz_yzzzz[k] = -g_x_0_y_0_xz_yzzzz[k] * ab_z + g_x_0_y_0_xz_yzzzzz[k];

                g_x_0_y_0_xzz_zzzzz[k] = -g_x_0_y_0_xz_zzzzz[k] * ab_z + g_x_0_y_0_xz_zzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 336 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 337 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 338 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 339 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 340 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 341 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 342 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 343 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 344 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 345 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 346 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 347 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 348 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 349 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 350 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 351 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 352 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 353 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 354 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 355 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yy_xxxxx, g_x_0_y_0_yy_xxxxxy, g_x_0_y_0_yy_xxxxy, g_x_0_y_0_yy_xxxxyy, g_x_0_y_0_yy_xxxxyz, g_x_0_y_0_yy_xxxxz, g_x_0_y_0_yy_xxxyy, g_x_0_y_0_yy_xxxyyy, g_x_0_y_0_yy_xxxyyz, g_x_0_y_0_yy_xxxyz, g_x_0_y_0_yy_xxxyzz, g_x_0_y_0_yy_xxxzz, g_x_0_y_0_yy_xxyyy, g_x_0_y_0_yy_xxyyyy, g_x_0_y_0_yy_xxyyyz, g_x_0_y_0_yy_xxyyz, g_x_0_y_0_yy_xxyyzz, g_x_0_y_0_yy_xxyzz, g_x_0_y_0_yy_xxyzzz, g_x_0_y_0_yy_xxzzz, g_x_0_y_0_yy_xyyyy, g_x_0_y_0_yy_xyyyyy, g_x_0_y_0_yy_xyyyyz, g_x_0_y_0_yy_xyyyz, g_x_0_y_0_yy_xyyyzz, g_x_0_y_0_yy_xyyzz, g_x_0_y_0_yy_xyyzzz, g_x_0_y_0_yy_xyzzz, g_x_0_y_0_yy_xyzzzz, g_x_0_y_0_yy_xzzzz, g_x_0_y_0_yy_yyyyy, g_x_0_y_0_yy_yyyyyy, g_x_0_y_0_yy_yyyyyz, g_x_0_y_0_yy_yyyyz, g_x_0_y_0_yy_yyyyzz, g_x_0_y_0_yy_yyyzz, g_x_0_y_0_yy_yyyzzz, g_x_0_y_0_yy_yyzzz, g_x_0_y_0_yy_yyzzzz, g_x_0_y_0_yy_yzzzz, g_x_0_y_0_yy_yzzzzz, g_x_0_y_0_yy_zzzzz, g_x_0_y_0_yyy_xxxxx, g_x_0_y_0_yyy_xxxxy, g_x_0_y_0_yyy_xxxxz, g_x_0_y_0_yyy_xxxyy, g_x_0_y_0_yyy_xxxyz, g_x_0_y_0_yyy_xxxzz, g_x_0_y_0_yyy_xxyyy, g_x_0_y_0_yyy_xxyyz, g_x_0_y_0_yyy_xxyzz, g_x_0_y_0_yyy_xxzzz, g_x_0_y_0_yyy_xyyyy, g_x_0_y_0_yyy_xyyyz, g_x_0_y_0_yyy_xyyzz, g_x_0_y_0_yyy_xyzzz, g_x_0_y_0_yyy_xzzzz, g_x_0_y_0_yyy_yyyyy, g_x_0_y_0_yyy_yyyyz, g_x_0_y_0_yyy_yyyzz, g_x_0_y_0_yyy_yyzzz, g_x_0_y_0_yyy_yzzzz, g_x_0_y_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyy_xxxxx[k] = -g_x_0_y_0_yy_xxxxx[k] * ab_y + g_x_0_y_0_yy_xxxxxy[k];

                g_x_0_y_0_yyy_xxxxy[k] = -g_x_0_y_0_yy_xxxxy[k] * ab_y + g_x_0_y_0_yy_xxxxyy[k];

                g_x_0_y_0_yyy_xxxxz[k] = -g_x_0_y_0_yy_xxxxz[k] * ab_y + g_x_0_y_0_yy_xxxxyz[k];

                g_x_0_y_0_yyy_xxxyy[k] = -g_x_0_y_0_yy_xxxyy[k] * ab_y + g_x_0_y_0_yy_xxxyyy[k];

                g_x_0_y_0_yyy_xxxyz[k] = -g_x_0_y_0_yy_xxxyz[k] * ab_y + g_x_0_y_0_yy_xxxyyz[k];

                g_x_0_y_0_yyy_xxxzz[k] = -g_x_0_y_0_yy_xxxzz[k] * ab_y + g_x_0_y_0_yy_xxxyzz[k];

                g_x_0_y_0_yyy_xxyyy[k] = -g_x_0_y_0_yy_xxyyy[k] * ab_y + g_x_0_y_0_yy_xxyyyy[k];

                g_x_0_y_0_yyy_xxyyz[k] = -g_x_0_y_0_yy_xxyyz[k] * ab_y + g_x_0_y_0_yy_xxyyyz[k];

                g_x_0_y_0_yyy_xxyzz[k] = -g_x_0_y_0_yy_xxyzz[k] * ab_y + g_x_0_y_0_yy_xxyyzz[k];

                g_x_0_y_0_yyy_xxzzz[k] = -g_x_0_y_0_yy_xxzzz[k] * ab_y + g_x_0_y_0_yy_xxyzzz[k];

                g_x_0_y_0_yyy_xyyyy[k] = -g_x_0_y_0_yy_xyyyy[k] * ab_y + g_x_0_y_0_yy_xyyyyy[k];

                g_x_0_y_0_yyy_xyyyz[k] = -g_x_0_y_0_yy_xyyyz[k] * ab_y + g_x_0_y_0_yy_xyyyyz[k];

                g_x_0_y_0_yyy_xyyzz[k] = -g_x_0_y_0_yy_xyyzz[k] * ab_y + g_x_0_y_0_yy_xyyyzz[k];

                g_x_0_y_0_yyy_xyzzz[k] = -g_x_0_y_0_yy_xyzzz[k] * ab_y + g_x_0_y_0_yy_xyyzzz[k];

                g_x_0_y_0_yyy_xzzzz[k] = -g_x_0_y_0_yy_xzzzz[k] * ab_y + g_x_0_y_0_yy_xyzzzz[k];

                g_x_0_y_0_yyy_yyyyy[k] = -g_x_0_y_0_yy_yyyyy[k] * ab_y + g_x_0_y_0_yy_yyyyyy[k];

                g_x_0_y_0_yyy_yyyyz[k] = -g_x_0_y_0_yy_yyyyz[k] * ab_y + g_x_0_y_0_yy_yyyyyz[k];

                g_x_0_y_0_yyy_yyyzz[k] = -g_x_0_y_0_yy_yyyzz[k] * ab_y + g_x_0_y_0_yy_yyyyzz[k];

                g_x_0_y_0_yyy_yyzzz[k] = -g_x_0_y_0_yy_yyzzz[k] * ab_y + g_x_0_y_0_yy_yyyzzz[k];

                g_x_0_y_0_yyy_yzzzz[k] = -g_x_0_y_0_yy_yzzzz[k] * ab_y + g_x_0_y_0_yy_yyzzzz[k];

                g_x_0_y_0_yyy_zzzzz[k] = -g_x_0_y_0_yy_zzzzz[k] * ab_y + g_x_0_y_0_yy_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 357 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 358 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 359 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 360 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 361 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 362 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 363 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 364 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 365 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 366 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 367 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 368 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 369 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 370 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 371 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 372 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 373 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 374 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 375 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 376 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyz_xxxxx, g_x_0_y_0_yyz_xxxxy, g_x_0_y_0_yyz_xxxxz, g_x_0_y_0_yyz_xxxyy, g_x_0_y_0_yyz_xxxyz, g_x_0_y_0_yyz_xxxzz, g_x_0_y_0_yyz_xxyyy, g_x_0_y_0_yyz_xxyyz, g_x_0_y_0_yyz_xxyzz, g_x_0_y_0_yyz_xxzzz, g_x_0_y_0_yyz_xyyyy, g_x_0_y_0_yyz_xyyyz, g_x_0_y_0_yyz_xyyzz, g_x_0_y_0_yyz_xyzzz, g_x_0_y_0_yyz_xzzzz, g_x_0_y_0_yyz_yyyyy, g_x_0_y_0_yyz_yyyyz, g_x_0_y_0_yyz_yyyzz, g_x_0_y_0_yyz_yyzzz, g_x_0_y_0_yyz_yzzzz, g_x_0_y_0_yyz_zzzzz, g_x_0_y_0_yz_xxxxx, g_x_0_y_0_yz_xxxxxy, g_x_0_y_0_yz_xxxxy, g_x_0_y_0_yz_xxxxyy, g_x_0_y_0_yz_xxxxyz, g_x_0_y_0_yz_xxxxz, g_x_0_y_0_yz_xxxyy, g_x_0_y_0_yz_xxxyyy, g_x_0_y_0_yz_xxxyyz, g_x_0_y_0_yz_xxxyz, g_x_0_y_0_yz_xxxyzz, g_x_0_y_0_yz_xxxzz, g_x_0_y_0_yz_xxyyy, g_x_0_y_0_yz_xxyyyy, g_x_0_y_0_yz_xxyyyz, g_x_0_y_0_yz_xxyyz, g_x_0_y_0_yz_xxyyzz, g_x_0_y_0_yz_xxyzz, g_x_0_y_0_yz_xxyzzz, g_x_0_y_0_yz_xxzzz, g_x_0_y_0_yz_xyyyy, g_x_0_y_0_yz_xyyyyy, g_x_0_y_0_yz_xyyyyz, g_x_0_y_0_yz_xyyyz, g_x_0_y_0_yz_xyyyzz, g_x_0_y_0_yz_xyyzz, g_x_0_y_0_yz_xyyzzz, g_x_0_y_0_yz_xyzzz, g_x_0_y_0_yz_xyzzzz, g_x_0_y_0_yz_xzzzz, g_x_0_y_0_yz_yyyyy, g_x_0_y_0_yz_yyyyyy, g_x_0_y_0_yz_yyyyyz, g_x_0_y_0_yz_yyyyz, g_x_0_y_0_yz_yyyyzz, g_x_0_y_0_yz_yyyzz, g_x_0_y_0_yz_yyyzzz, g_x_0_y_0_yz_yyzzz, g_x_0_y_0_yz_yyzzzz, g_x_0_y_0_yz_yzzzz, g_x_0_y_0_yz_yzzzzz, g_x_0_y_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyz_xxxxx[k] = -g_x_0_y_0_yz_xxxxx[k] * ab_y + g_x_0_y_0_yz_xxxxxy[k];

                g_x_0_y_0_yyz_xxxxy[k] = -g_x_0_y_0_yz_xxxxy[k] * ab_y + g_x_0_y_0_yz_xxxxyy[k];

                g_x_0_y_0_yyz_xxxxz[k] = -g_x_0_y_0_yz_xxxxz[k] * ab_y + g_x_0_y_0_yz_xxxxyz[k];

                g_x_0_y_0_yyz_xxxyy[k] = -g_x_0_y_0_yz_xxxyy[k] * ab_y + g_x_0_y_0_yz_xxxyyy[k];

                g_x_0_y_0_yyz_xxxyz[k] = -g_x_0_y_0_yz_xxxyz[k] * ab_y + g_x_0_y_0_yz_xxxyyz[k];

                g_x_0_y_0_yyz_xxxzz[k] = -g_x_0_y_0_yz_xxxzz[k] * ab_y + g_x_0_y_0_yz_xxxyzz[k];

                g_x_0_y_0_yyz_xxyyy[k] = -g_x_0_y_0_yz_xxyyy[k] * ab_y + g_x_0_y_0_yz_xxyyyy[k];

                g_x_0_y_0_yyz_xxyyz[k] = -g_x_0_y_0_yz_xxyyz[k] * ab_y + g_x_0_y_0_yz_xxyyyz[k];

                g_x_0_y_0_yyz_xxyzz[k] = -g_x_0_y_0_yz_xxyzz[k] * ab_y + g_x_0_y_0_yz_xxyyzz[k];

                g_x_0_y_0_yyz_xxzzz[k] = -g_x_0_y_0_yz_xxzzz[k] * ab_y + g_x_0_y_0_yz_xxyzzz[k];

                g_x_0_y_0_yyz_xyyyy[k] = -g_x_0_y_0_yz_xyyyy[k] * ab_y + g_x_0_y_0_yz_xyyyyy[k];

                g_x_0_y_0_yyz_xyyyz[k] = -g_x_0_y_0_yz_xyyyz[k] * ab_y + g_x_0_y_0_yz_xyyyyz[k];

                g_x_0_y_0_yyz_xyyzz[k] = -g_x_0_y_0_yz_xyyzz[k] * ab_y + g_x_0_y_0_yz_xyyyzz[k];

                g_x_0_y_0_yyz_xyzzz[k] = -g_x_0_y_0_yz_xyzzz[k] * ab_y + g_x_0_y_0_yz_xyyzzz[k];

                g_x_0_y_0_yyz_xzzzz[k] = -g_x_0_y_0_yz_xzzzz[k] * ab_y + g_x_0_y_0_yz_xyzzzz[k];

                g_x_0_y_0_yyz_yyyyy[k] = -g_x_0_y_0_yz_yyyyy[k] * ab_y + g_x_0_y_0_yz_yyyyyy[k];

                g_x_0_y_0_yyz_yyyyz[k] = -g_x_0_y_0_yz_yyyyz[k] * ab_y + g_x_0_y_0_yz_yyyyyz[k];

                g_x_0_y_0_yyz_yyyzz[k] = -g_x_0_y_0_yz_yyyzz[k] * ab_y + g_x_0_y_0_yz_yyyyzz[k];

                g_x_0_y_0_yyz_yyzzz[k] = -g_x_0_y_0_yz_yyzzz[k] * ab_y + g_x_0_y_0_yz_yyyzzz[k];

                g_x_0_y_0_yyz_yzzzz[k] = -g_x_0_y_0_yz_yzzzz[k] * ab_y + g_x_0_y_0_yz_yyzzzz[k];

                g_x_0_y_0_yyz_zzzzz[k] = -g_x_0_y_0_yz_zzzzz[k] * ab_y + g_x_0_y_0_yz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 378 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 379 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 380 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 381 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 382 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 383 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 384 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 385 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 386 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 387 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 388 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 389 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 390 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 391 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 392 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 393 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 394 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 395 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 396 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 397 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzz_xxxxx, g_x_0_y_0_yzz_xxxxy, g_x_0_y_0_yzz_xxxxz, g_x_0_y_0_yzz_xxxyy, g_x_0_y_0_yzz_xxxyz, g_x_0_y_0_yzz_xxxzz, g_x_0_y_0_yzz_xxyyy, g_x_0_y_0_yzz_xxyyz, g_x_0_y_0_yzz_xxyzz, g_x_0_y_0_yzz_xxzzz, g_x_0_y_0_yzz_xyyyy, g_x_0_y_0_yzz_xyyyz, g_x_0_y_0_yzz_xyyzz, g_x_0_y_0_yzz_xyzzz, g_x_0_y_0_yzz_xzzzz, g_x_0_y_0_yzz_yyyyy, g_x_0_y_0_yzz_yyyyz, g_x_0_y_0_yzz_yyyzz, g_x_0_y_0_yzz_yyzzz, g_x_0_y_0_yzz_yzzzz, g_x_0_y_0_yzz_zzzzz, g_x_0_y_0_zz_xxxxx, g_x_0_y_0_zz_xxxxxy, g_x_0_y_0_zz_xxxxy, g_x_0_y_0_zz_xxxxyy, g_x_0_y_0_zz_xxxxyz, g_x_0_y_0_zz_xxxxz, g_x_0_y_0_zz_xxxyy, g_x_0_y_0_zz_xxxyyy, g_x_0_y_0_zz_xxxyyz, g_x_0_y_0_zz_xxxyz, g_x_0_y_0_zz_xxxyzz, g_x_0_y_0_zz_xxxzz, g_x_0_y_0_zz_xxyyy, g_x_0_y_0_zz_xxyyyy, g_x_0_y_0_zz_xxyyyz, g_x_0_y_0_zz_xxyyz, g_x_0_y_0_zz_xxyyzz, g_x_0_y_0_zz_xxyzz, g_x_0_y_0_zz_xxyzzz, g_x_0_y_0_zz_xxzzz, g_x_0_y_0_zz_xyyyy, g_x_0_y_0_zz_xyyyyy, g_x_0_y_0_zz_xyyyyz, g_x_0_y_0_zz_xyyyz, g_x_0_y_0_zz_xyyyzz, g_x_0_y_0_zz_xyyzz, g_x_0_y_0_zz_xyyzzz, g_x_0_y_0_zz_xyzzz, g_x_0_y_0_zz_xyzzzz, g_x_0_y_0_zz_xzzzz, g_x_0_y_0_zz_yyyyy, g_x_0_y_0_zz_yyyyyy, g_x_0_y_0_zz_yyyyyz, g_x_0_y_0_zz_yyyyz, g_x_0_y_0_zz_yyyyzz, g_x_0_y_0_zz_yyyzz, g_x_0_y_0_zz_yyyzzz, g_x_0_y_0_zz_yyzzz, g_x_0_y_0_zz_yyzzzz, g_x_0_y_0_zz_yzzzz, g_x_0_y_0_zz_yzzzzz, g_x_0_y_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzz_xxxxx[k] = -g_x_0_y_0_zz_xxxxx[k] * ab_y + g_x_0_y_0_zz_xxxxxy[k];

                g_x_0_y_0_yzz_xxxxy[k] = -g_x_0_y_0_zz_xxxxy[k] * ab_y + g_x_0_y_0_zz_xxxxyy[k];

                g_x_0_y_0_yzz_xxxxz[k] = -g_x_0_y_0_zz_xxxxz[k] * ab_y + g_x_0_y_0_zz_xxxxyz[k];

                g_x_0_y_0_yzz_xxxyy[k] = -g_x_0_y_0_zz_xxxyy[k] * ab_y + g_x_0_y_0_zz_xxxyyy[k];

                g_x_0_y_0_yzz_xxxyz[k] = -g_x_0_y_0_zz_xxxyz[k] * ab_y + g_x_0_y_0_zz_xxxyyz[k];

                g_x_0_y_0_yzz_xxxzz[k] = -g_x_0_y_0_zz_xxxzz[k] * ab_y + g_x_0_y_0_zz_xxxyzz[k];

                g_x_0_y_0_yzz_xxyyy[k] = -g_x_0_y_0_zz_xxyyy[k] * ab_y + g_x_0_y_0_zz_xxyyyy[k];

                g_x_0_y_0_yzz_xxyyz[k] = -g_x_0_y_0_zz_xxyyz[k] * ab_y + g_x_0_y_0_zz_xxyyyz[k];

                g_x_0_y_0_yzz_xxyzz[k] = -g_x_0_y_0_zz_xxyzz[k] * ab_y + g_x_0_y_0_zz_xxyyzz[k];

                g_x_0_y_0_yzz_xxzzz[k] = -g_x_0_y_0_zz_xxzzz[k] * ab_y + g_x_0_y_0_zz_xxyzzz[k];

                g_x_0_y_0_yzz_xyyyy[k] = -g_x_0_y_0_zz_xyyyy[k] * ab_y + g_x_0_y_0_zz_xyyyyy[k];

                g_x_0_y_0_yzz_xyyyz[k] = -g_x_0_y_0_zz_xyyyz[k] * ab_y + g_x_0_y_0_zz_xyyyyz[k];

                g_x_0_y_0_yzz_xyyzz[k] = -g_x_0_y_0_zz_xyyzz[k] * ab_y + g_x_0_y_0_zz_xyyyzz[k];

                g_x_0_y_0_yzz_xyzzz[k] = -g_x_0_y_0_zz_xyzzz[k] * ab_y + g_x_0_y_0_zz_xyyzzz[k];

                g_x_0_y_0_yzz_xzzzz[k] = -g_x_0_y_0_zz_xzzzz[k] * ab_y + g_x_0_y_0_zz_xyzzzz[k];

                g_x_0_y_0_yzz_yyyyy[k] = -g_x_0_y_0_zz_yyyyy[k] * ab_y + g_x_0_y_0_zz_yyyyyy[k];

                g_x_0_y_0_yzz_yyyyz[k] = -g_x_0_y_0_zz_yyyyz[k] * ab_y + g_x_0_y_0_zz_yyyyyz[k];

                g_x_0_y_0_yzz_yyyzz[k] = -g_x_0_y_0_zz_yyyzz[k] * ab_y + g_x_0_y_0_zz_yyyyzz[k];

                g_x_0_y_0_yzz_yyzzz[k] = -g_x_0_y_0_zz_yyzzz[k] * ab_y + g_x_0_y_0_zz_yyyzzz[k];

                g_x_0_y_0_yzz_yzzzz[k] = -g_x_0_y_0_zz_yzzzz[k] * ab_y + g_x_0_y_0_zz_yyzzzz[k];

                g_x_0_y_0_yzz_zzzzz[k] = -g_x_0_y_0_zz_zzzzz[k] * ab_y + g_x_0_y_0_zz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 399 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 400 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 401 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 402 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 403 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 404 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 405 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 406 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 407 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 408 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 409 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 410 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 411 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 412 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 413 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 414 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 415 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 416 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 417 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 418 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zz_xxxxx, g_x_0_y_0_zz_xxxxxz, g_x_0_y_0_zz_xxxxy, g_x_0_y_0_zz_xxxxyz, g_x_0_y_0_zz_xxxxz, g_x_0_y_0_zz_xxxxzz, g_x_0_y_0_zz_xxxyy, g_x_0_y_0_zz_xxxyyz, g_x_0_y_0_zz_xxxyz, g_x_0_y_0_zz_xxxyzz, g_x_0_y_0_zz_xxxzz, g_x_0_y_0_zz_xxxzzz, g_x_0_y_0_zz_xxyyy, g_x_0_y_0_zz_xxyyyz, g_x_0_y_0_zz_xxyyz, g_x_0_y_0_zz_xxyyzz, g_x_0_y_0_zz_xxyzz, g_x_0_y_0_zz_xxyzzz, g_x_0_y_0_zz_xxzzz, g_x_0_y_0_zz_xxzzzz, g_x_0_y_0_zz_xyyyy, g_x_0_y_0_zz_xyyyyz, g_x_0_y_0_zz_xyyyz, g_x_0_y_0_zz_xyyyzz, g_x_0_y_0_zz_xyyzz, g_x_0_y_0_zz_xyyzzz, g_x_0_y_0_zz_xyzzz, g_x_0_y_0_zz_xyzzzz, g_x_0_y_0_zz_xzzzz, g_x_0_y_0_zz_xzzzzz, g_x_0_y_0_zz_yyyyy, g_x_0_y_0_zz_yyyyyz, g_x_0_y_0_zz_yyyyz, g_x_0_y_0_zz_yyyyzz, g_x_0_y_0_zz_yyyzz, g_x_0_y_0_zz_yyyzzz, g_x_0_y_0_zz_yyzzz, g_x_0_y_0_zz_yyzzzz, g_x_0_y_0_zz_yzzzz, g_x_0_y_0_zz_yzzzzz, g_x_0_y_0_zz_zzzzz, g_x_0_y_0_zz_zzzzzz, g_x_0_y_0_zzz_xxxxx, g_x_0_y_0_zzz_xxxxy, g_x_0_y_0_zzz_xxxxz, g_x_0_y_0_zzz_xxxyy, g_x_0_y_0_zzz_xxxyz, g_x_0_y_0_zzz_xxxzz, g_x_0_y_0_zzz_xxyyy, g_x_0_y_0_zzz_xxyyz, g_x_0_y_0_zzz_xxyzz, g_x_0_y_0_zzz_xxzzz, g_x_0_y_0_zzz_xyyyy, g_x_0_y_0_zzz_xyyyz, g_x_0_y_0_zzz_xyyzz, g_x_0_y_0_zzz_xyzzz, g_x_0_y_0_zzz_xzzzz, g_x_0_y_0_zzz_yyyyy, g_x_0_y_0_zzz_yyyyz, g_x_0_y_0_zzz_yyyzz, g_x_0_y_0_zzz_yyzzz, g_x_0_y_0_zzz_yzzzz, g_x_0_y_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzz_xxxxx[k] = -g_x_0_y_0_zz_xxxxx[k] * ab_z + g_x_0_y_0_zz_xxxxxz[k];

                g_x_0_y_0_zzz_xxxxy[k] = -g_x_0_y_0_zz_xxxxy[k] * ab_z + g_x_0_y_0_zz_xxxxyz[k];

                g_x_0_y_0_zzz_xxxxz[k] = -g_x_0_y_0_zz_xxxxz[k] * ab_z + g_x_0_y_0_zz_xxxxzz[k];

                g_x_0_y_0_zzz_xxxyy[k] = -g_x_0_y_0_zz_xxxyy[k] * ab_z + g_x_0_y_0_zz_xxxyyz[k];

                g_x_0_y_0_zzz_xxxyz[k] = -g_x_0_y_0_zz_xxxyz[k] * ab_z + g_x_0_y_0_zz_xxxyzz[k];

                g_x_0_y_0_zzz_xxxzz[k] = -g_x_0_y_0_zz_xxxzz[k] * ab_z + g_x_0_y_0_zz_xxxzzz[k];

                g_x_0_y_0_zzz_xxyyy[k] = -g_x_0_y_0_zz_xxyyy[k] * ab_z + g_x_0_y_0_zz_xxyyyz[k];

                g_x_0_y_0_zzz_xxyyz[k] = -g_x_0_y_0_zz_xxyyz[k] * ab_z + g_x_0_y_0_zz_xxyyzz[k];

                g_x_0_y_0_zzz_xxyzz[k] = -g_x_0_y_0_zz_xxyzz[k] * ab_z + g_x_0_y_0_zz_xxyzzz[k];

                g_x_0_y_0_zzz_xxzzz[k] = -g_x_0_y_0_zz_xxzzz[k] * ab_z + g_x_0_y_0_zz_xxzzzz[k];

                g_x_0_y_0_zzz_xyyyy[k] = -g_x_0_y_0_zz_xyyyy[k] * ab_z + g_x_0_y_0_zz_xyyyyz[k];

                g_x_0_y_0_zzz_xyyyz[k] = -g_x_0_y_0_zz_xyyyz[k] * ab_z + g_x_0_y_0_zz_xyyyzz[k];

                g_x_0_y_0_zzz_xyyzz[k] = -g_x_0_y_0_zz_xyyzz[k] * ab_z + g_x_0_y_0_zz_xyyzzz[k];

                g_x_0_y_0_zzz_xyzzz[k] = -g_x_0_y_0_zz_xyzzz[k] * ab_z + g_x_0_y_0_zz_xyzzzz[k];

                g_x_0_y_0_zzz_xzzzz[k] = -g_x_0_y_0_zz_xzzzz[k] * ab_z + g_x_0_y_0_zz_xzzzzz[k];

                g_x_0_y_0_zzz_yyyyy[k] = -g_x_0_y_0_zz_yyyyy[k] * ab_z + g_x_0_y_0_zz_yyyyyz[k];

                g_x_0_y_0_zzz_yyyyz[k] = -g_x_0_y_0_zz_yyyyz[k] * ab_z + g_x_0_y_0_zz_yyyyzz[k];

                g_x_0_y_0_zzz_yyyzz[k] = -g_x_0_y_0_zz_yyyzz[k] * ab_z + g_x_0_y_0_zz_yyyzzz[k];

                g_x_0_y_0_zzz_yyzzz[k] = -g_x_0_y_0_zz_yyzzz[k] * ab_z + g_x_0_y_0_zz_yyzzzz[k];

                g_x_0_y_0_zzz_yzzzz[k] = -g_x_0_y_0_zz_yzzzz[k] * ab_z + g_x_0_y_0_zz_yzzzzz[k];

                g_x_0_y_0_zzz_zzzzz[k] = -g_x_0_y_0_zz_zzzzz[k] * ab_z + g_x_0_y_0_zz_zzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 420 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 421 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 422 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 423 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 424 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 425 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 426 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 427 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 428 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 429 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 430 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 431 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 432 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 433 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 434 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 435 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 436 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 437 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 438 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 439 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xx_xxxxx, g_0_0_z_0_xx_xxxxy, g_0_0_z_0_xx_xxxxz, g_0_0_z_0_xx_xxxyy, g_0_0_z_0_xx_xxxyz, g_0_0_z_0_xx_xxxzz, g_0_0_z_0_xx_xxyyy, g_0_0_z_0_xx_xxyyz, g_0_0_z_0_xx_xxyzz, g_0_0_z_0_xx_xxzzz, g_0_0_z_0_xx_xyyyy, g_0_0_z_0_xx_xyyyz, g_0_0_z_0_xx_xyyzz, g_0_0_z_0_xx_xyzzz, g_0_0_z_0_xx_xzzzz, g_0_0_z_0_xx_yyyyy, g_0_0_z_0_xx_yyyyz, g_0_0_z_0_xx_yyyzz, g_0_0_z_0_xx_yyzzz, g_0_0_z_0_xx_yzzzz, g_0_0_z_0_xx_zzzzz, g_x_0_z_0_xx_xxxxx, g_x_0_z_0_xx_xxxxxx, g_x_0_z_0_xx_xxxxxy, g_x_0_z_0_xx_xxxxxz, g_x_0_z_0_xx_xxxxy, g_x_0_z_0_xx_xxxxyy, g_x_0_z_0_xx_xxxxyz, g_x_0_z_0_xx_xxxxz, g_x_0_z_0_xx_xxxxzz, g_x_0_z_0_xx_xxxyy, g_x_0_z_0_xx_xxxyyy, g_x_0_z_0_xx_xxxyyz, g_x_0_z_0_xx_xxxyz, g_x_0_z_0_xx_xxxyzz, g_x_0_z_0_xx_xxxzz, g_x_0_z_0_xx_xxxzzz, g_x_0_z_0_xx_xxyyy, g_x_0_z_0_xx_xxyyyy, g_x_0_z_0_xx_xxyyyz, g_x_0_z_0_xx_xxyyz, g_x_0_z_0_xx_xxyyzz, g_x_0_z_0_xx_xxyzz, g_x_0_z_0_xx_xxyzzz, g_x_0_z_0_xx_xxzzz, g_x_0_z_0_xx_xxzzzz, g_x_0_z_0_xx_xyyyy, g_x_0_z_0_xx_xyyyyy, g_x_0_z_0_xx_xyyyyz, g_x_0_z_0_xx_xyyyz, g_x_0_z_0_xx_xyyyzz, g_x_0_z_0_xx_xyyzz, g_x_0_z_0_xx_xyyzzz, g_x_0_z_0_xx_xyzzz, g_x_0_z_0_xx_xyzzzz, g_x_0_z_0_xx_xzzzz, g_x_0_z_0_xx_xzzzzz, g_x_0_z_0_xx_yyyyy, g_x_0_z_0_xx_yyyyz, g_x_0_z_0_xx_yyyzz, g_x_0_z_0_xx_yyzzz, g_x_0_z_0_xx_yzzzz, g_x_0_z_0_xx_zzzzz, g_x_0_z_0_xxx_xxxxx, g_x_0_z_0_xxx_xxxxy, g_x_0_z_0_xxx_xxxxz, g_x_0_z_0_xxx_xxxyy, g_x_0_z_0_xxx_xxxyz, g_x_0_z_0_xxx_xxxzz, g_x_0_z_0_xxx_xxyyy, g_x_0_z_0_xxx_xxyyz, g_x_0_z_0_xxx_xxyzz, g_x_0_z_0_xxx_xxzzz, g_x_0_z_0_xxx_xyyyy, g_x_0_z_0_xxx_xyyyz, g_x_0_z_0_xxx_xyyzz, g_x_0_z_0_xxx_xyzzz, g_x_0_z_0_xxx_xzzzz, g_x_0_z_0_xxx_yyyyy, g_x_0_z_0_xxx_yyyyz, g_x_0_z_0_xxx_yyyzz, g_x_0_z_0_xxx_yyzzz, g_x_0_z_0_xxx_yzzzz, g_x_0_z_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxx_xxxxx[k] = -g_0_0_z_0_xx_xxxxx[k] - g_x_0_z_0_xx_xxxxx[k] * ab_x + g_x_0_z_0_xx_xxxxxx[k];

                g_x_0_z_0_xxx_xxxxy[k] = -g_0_0_z_0_xx_xxxxy[k] - g_x_0_z_0_xx_xxxxy[k] * ab_x + g_x_0_z_0_xx_xxxxxy[k];

                g_x_0_z_0_xxx_xxxxz[k] = -g_0_0_z_0_xx_xxxxz[k] - g_x_0_z_0_xx_xxxxz[k] * ab_x + g_x_0_z_0_xx_xxxxxz[k];

                g_x_0_z_0_xxx_xxxyy[k] = -g_0_0_z_0_xx_xxxyy[k] - g_x_0_z_0_xx_xxxyy[k] * ab_x + g_x_0_z_0_xx_xxxxyy[k];

                g_x_0_z_0_xxx_xxxyz[k] = -g_0_0_z_0_xx_xxxyz[k] - g_x_0_z_0_xx_xxxyz[k] * ab_x + g_x_0_z_0_xx_xxxxyz[k];

                g_x_0_z_0_xxx_xxxzz[k] = -g_0_0_z_0_xx_xxxzz[k] - g_x_0_z_0_xx_xxxzz[k] * ab_x + g_x_0_z_0_xx_xxxxzz[k];

                g_x_0_z_0_xxx_xxyyy[k] = -g_0_0_z_0_xx_xxyyy[k] - g_x_0_z_0_xx_xxyyy[k] * ab_x + g_x_0_z_0_xx_xxxyyy[k];

                g_x_0_z_0_xxx_xxyyz[k] = -g_0_0_z_0_xx_xxyyz[k] - g_x_0_z_0_xx_xxyyz[k] * ab_x + g_x_0_z_0_xx_xxxyyz[k];

                g_x_0_z_0_xxx_xxyzz[k] = -g_0_0_z_0_xx_xxyzz[k] - g_x_0_z_0_xx_xxyzz[k] * ab_x + g_x_0_z_0_xx_xxxyzz[k];

                g_x_0_z_0_xxx_xxzzz[k] = -g_0_0_z_0_xx_xxzzz[k] - g_x_0_z_0_xx_xxzzz[k] * ab_x + g_x_0_z_0_xx_xxxzzz[k];

                g_x_0_z_0_xxx_xyyyy[k] = -g_0_0_z_0_xx_xyyyy[k] - g_x_0_z_0_xx_xyyyy[k] * ab_x + g_x_0_z_0_xx_xxyyyy[k];

                g_x_0_z_0_xxx_xyyyz[k] = -g_0_0_z_0_xx_xyyyz[k] - g_x_0_z_0_xx_xyyyz[k] * ab_x + g_x_0_z_0_xx_xxyyyz[k];

                g_x_0_z_0_xxx_xyyzz[k] = -g_0_0_z_0_xx_xyyzz[k] - g_x_0_z_0_xx_xyyzz[k] * ab_x + g_x_0_z_0_xx_xxyyzz[k];

                g_x_0_z_0_xxx_xyzzz[k] = -g_0_0_z_0_xx_xyzzz[k] - g_x_0_z_0_xx_xyzzz[k] * ab_x + g_x_0_z_0_xx_xxyzzz[k];

                g_x_0_z_0_xxx_xzzzz[k] = -g_0_0_z_0_xx_xzzzz[k] - g_x_0_z_0_xx_xzzzz[k] * ab_x + g_x_0_z_0_xx_xxzzzz[k];

                g_x_0_z_0_xxx_yyyyy[k] = -g_0_0_z_0_xx_yyyyy[k] - g_x_0_z_0_xx_yyyyy[k] * ab_x + g_x_0_z_0_xx_xyyyyy[k];

                g_x_0_z_0_xxx_yyyyz[k] = -g_0_0_z_0_xx_yyyyz[k] - g_x_0_z_0_xx_yyyyz[k] * ab_x + g_x_0_z_0_xx_xyyyyz[k];

                g_x_0_z_0_xxx_yyyzz[k] = -g_0_0_z_0_xx_yyyzz[k] - g_x_0_z_0_xx_yyyzz[k] * ab_x + g_x_0_z_0_xx_xyyyzz[k];

                g_x_0_z_0_xxx_yyzzz[k] = -g_0_0_z_0_xx_yyzzz[k] - g_x_0_z_0_xx_yyzzz[k] * ab_x + g_x_0_z_0_xx_xyyzzz[k];

                g_x_0_z_0_xxx_yzzzz[k] = -g_0_0_z_0_xx_yzzzz[k] - g_x_0_z_0_xx_yzzzz[k] * ab_x + g_x_0_z_0_xx_xyzzzz[k];

                g_x_0_z_0_xxx_zzzzz[k] = -g_0_0_z_0_xx_zzzzz[k] - g_x_0_z_0_xx_zzzzz[k] * ab_x + g_x_0_z_0_xx_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 441 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 442 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 443 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 444 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 445 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 446 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 447 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 448 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 449 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 450 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 451 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 452 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 453 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 454 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 455 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 456 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 457 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 458 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 459 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 460 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xx_xxxxx, g_x_0_z_0_xx_xxxxxy, g_x_0_z_0_xx_xxxxy, g_x_0_z_0_xx_xxxxyy, g_x_0_z_0_xx_xxxxyz, g_x_0_z_0_xx_xxxxz, g_x_0_z_0_xx_xxxyy, g_x_0_z_0_xx_xxxyyy, g_x_0_z_0_xx_xxxyyz, g_x_0_z_0_xx_xxxyz, g_x_0_z_0_xx_xxxyzz, g_x_0_z_0_xx_xxxzz, g_x_0_z_0_xx_xxyyy, g_x_0_z_0_xx_xxyyyy, g_x_0_z_0_xx_xxyyyz, g_x_0_z_0_xx_xxyyz, g_x_0_z_0_xx_xxyyzz, g_x_0_z_0_xx_xxyzz, g_x_0_z_0_xx_xxyzzz, g_x_0_z_0_xx_xxzzz, g_x_0_z_0_xx_xyyyy, g_x_0_z_0_xx_xyyyyy, g_x_0_z_0_xx_xyyyyz, g_x_0_z_0_xx_xyyyz, g_x_0_z_0_xx_xyyyzz, g_x_0_z_0_xx_xyyzz, g_x_0_z_0_xx_xyyzzz, g_x_0_z_0_xx_xyzzz, g_x_0_z_0_xx_xyzzzz, g_x_0_z_0_xx_xzzzz, g_x_0_z_0_xx_yyyyy, g_x_0_z_0_xx_yyyyyy, g_x_0_z_0_xx_yyyyyz, g_x_0_z_0_xx_yyyyz, g_x_0_z_0_xx_yyyyzz, g_x_0_z_0_xx_yyyzz, g_x_0_z_0_xx_yyyzzz, g_x_0_z_0_xx_yyzzz, g_x_0_z_0_xx_yyzzzz, g_x_0_z_0_xx_yzzzz, g_x_0_z_0_xx_yzzzzz, g_x_0_z_0_xx_zzzzz, g_x_0_z_0_xxy_xxxxx, g_x_0_z_0_xxy_xxxxy, g_x_0_z_0_xxy_xxxxz, g_x_0_z_0_xxy_xxxyy, g_x_0_z_0_xxy_xxxyz, g_x_0_z_0_xxy_xxxzz, g_x_0_z_0_xxy_xxyyy, g_x_0_z_0_xxy_xxyyz, g_x_0_z_0_xxy_xxyzz, g_x_0_z_0_xxy_xxzzz, g_x_0_z_0_xxy_xyyyy, g_x_0_z_0_xxy_xyyyz, g_x_0_z_0_xxy_xyyzz, g_x_0_z_0_xxy_xyzzz, g_x_0_z_0_xxy_xzzzz, g_x_0_z_0_xxy_yyyyy, g_x_0_z_0_xxy_yyyyz, g_x_0_z_0_xxy_yyyzz, g_x_0_z_0_xxy_yyzzz, g_x_0_z_0_xxy_yzzzz, g_x_0_z_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxy_xxxxx[k] = -g_x_0_z_0_xx_xxxxx[k] * ab_y + g_x_0_z_0_xx_xxxxxy[k];

                g_x_0_z_0_xxy_xxxxy[k] = -g_x_0_z_0_xx_xxxxy[k] * ab_y + g_x_0_z_0_xx_xxxxyy[k];

                g_x_0_z_0_xxy_xxxxz[k] = -g_x_0_z_0_xx_xxxxz[k] * ab_y + g_x_0_z_0_xx_xxxxyz[k];

                g_x_0_z_0_xxy_xxxyy[k] = -g_x_0_z_0_xx_xxxyy[k] * ab_y + g_x_0_z_0_xx_xxxyyy[k];

                g_x_0_z_0_xxy_xxxyz[k] = -g_x_0_z_0_xx_xxxyz[k] * ab_y + g_x_0_z_0_xx_xxxyyz[k];

                g_x_0_z_0_xxy_xxxzz[k] = -g_x_0_z_0_xx_xxxzz[k] * ab_y + g_x_0_z_0_xx_xxxyzz[k];

                g_x_0_z_0_xxy_xxyyy[k] = -g_x_0_z_0_xx_xxyyy[k] * ab_y + g_x_0_z_0_xx_xxyyyy[k];

                g_x_0_z_0_xxy_xxyyz[k] = -g_x_0_z_0_xx_xxyyz[k] * ab_y + g_x_0_z_0_xx_xxyyyz[k];

                g_x_0_z_0_xxy_xxyzz[k] = -g_x_0_z_0_xx_xxyzz[k] * ab_y + g_x_0_z_0_xx_xxyyzz[k];

                g_x_0_z_0_xxy_xxzzz[k] = -g_x_0_z_0_xx_xxzzz[k] * ab_y + g_x_0_z_0_xx_xxyzzz[k];

                g_x_0_z_0_xxy_xyyyy[k] = -g_x_0_z_0_xx_xyyyy[k] * ab_y + g_x_0_z_0_xx_xyyyyy[k];

                g_x_0_z_0_xxy_xyyyz[k] = -g_x_0_z_0_xx_xyyyz[k] * ab_y + g_x_0_z_0_xx_xyyyyz[k];

                g_x_0_z_0_xxy_xyyzz[k] = -g_x_0_z_0_xx_xyyzz[k] * ab_y + g_x_0_z_0_xx_xyyyzz[k];

                g_x_0_z_0_xxy_xyzzz[k] = -g_x_0_z_0_xx_xyzzz[k] * ab_y + g_x_0_z_0_xx_xyyzzz[k];

                g_x_0_z_0_xxy_xzzzz[k] = -g_x_0_z_0_xx_xzzzz[k] * ab_y + g_x_0_z_0_xx_xyzzzz[k];

                g_x_0_z_0_xxy_yyyyy[k] = -g_x_0_z_0_xx_yyyyy[k] * ab_y + g_x_0_z_0_xx_yyyyyy[k];

                g_x_0_z_0_xxy_yyyyz[k] = -g_x_0_z_0_xx_yyyyz[k] * ab_y + g_x_0_z_0_xx_yyyyyz[k];

                g_x_0_z_0_xxy_yyyzz[k] = -g_x_0_z_0_xx_yyyzz[k] * ab_y + g_x_0_z_0_xx_yyyyzz[k];

                g_x_0_z_0_xxy_yyzzz[k] = -g_x_0_z_0_xx_yyzzz[k] * ab_y + g_x_0_z_0_xx_yyyzzz[k];

                g_x_0_z_0_xxy_yzzzz[k] = -g_x_0_z_0_xx_yzzzz[k] * ab_y + g_x_0_z_0_xx_yyzzzz[k];

                g_x_0_z_0_xxy_zzzzz[k] = -g_x_0_z_0_xx_zzzzz[k] * ab_y + g_x_0_z_0_xx_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 462 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 463 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 464 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 465 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 466 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 467 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 468 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 469 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 470 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 471 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 472 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 473 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 474 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 475 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 476 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 477 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 478 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 479 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 480 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 481 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xx_xxxxx, g_x_0_z_0_xx_xxxxxz, g_x_0_z_0_xx_xxxxy, g_x_0_z_0_xx_xxxxyz, g_x_0_z_0_xx_xxxxz, g_x_0_z_0_xx_xxxxzz, g_x_0_z_0_xx_xxxyy, g_x_0_z_0_xx_xxxyyz, g_x_0_z_0_xx_xxxyz, g_x_0_z_0_xx_xxxyzz, g_x_0_z_0_xx_xxxzz, g_x_0_z_0_xx_xxxzzz, g_x_0_z_0_xx_xxyyy, g_x_0_z_0_xx_xxyyyz, g_x_0_z_0_xx_xxyyz, g_x_0_z_0_xx_xxyyzz, g_x_0_z_0_xx_xxyzz, g_x_0_z_0_xx_xxyzzz, g_x_0_z_0_xx_xxzzz, g_x_0_z_0_xx_xxzzzz, g_x_0_z_0_xx_xyyyy, g_x_0_z_0_xx_xyyyyz, g_x_0_z_0_xx_xyyyz, g_x_0_z_0_xx_xyyyzz, g_x_0_z_0_xx_xyyzz, g_x_0_z_0_xx_xyyzzz, g_x_0_z_0_xx_xyzzz, g_x_0_z_0_xx_xyzzzz, g_x_0_z_0_xx_xzzzz, g_x_0_z_0_xx_xzzzzz, g_x_0_z_0_xx_yyyyy, g_x_0_z_0_xx_yyyyyz, g_x_0_z_0_xx_yyyyz, g_x_0_z_0_xx_yyyyzz, g_x_0_z_0_xx_yyyzz, g_x_0_z_0_xx_yyyzzz, g_x_0_z_0_xx_yyzzz, g_x_0_z_0_xx_yyzzzz, g_x_0_z_0_xx_yzzzz, g_x_0_z_0_xx_yzzzzz, g_x_0_z_0_xx_zzzzz, g_x_0_z_0_xx_zzzzzz, g_x_0_z_0_xxz_xxxxx, g_x_0_z_0_xxz_xxxxy, g_x_0_z_0_xxz_xxxxz, g_x_0_z_0_xxz_xxxyy, g_x_0_z_0_xxz_xxxyz, g_x_0_z_0_xxz_xxxzz, g_x_0_z_0_xxz_xxyyy, g_x_0_z_0_xxz_xxyyz, g_x_0_z_0_xxz_xxyzz, g_x_0_z_0_xxz_xxzzz, g_x_0_z_0_xxz_xyyyy, g_x_0_z_0_xxz_xyyyz, g_x_0_z_0_xxz_xyyzz, g_x_0_z_0_xxz_xyzzz, g_x_0_z_0_xxz_xzzzz, g_x_0_z_0_xxz_yyyyy, g_x_0_z_0_xxz_yyyyz, g_x_0_z_0_xxz_yyyzz, g_x_0_z_0_xxz_yyzzz, g_x_0_z_0_xxz_yzzzz, g_x_0_z_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxz_xxxxx[k] = -g_x_0_z_0_xx_xxxxx[k] * ab_z + g_x_0_z_0_xx_xxxxxz[k];

                g_x_0_z_0_xxz_xxxxy[k] = -g_x_0_z_0_xx_xxxxy[k] * ab_z + g_x_0_z_0_xx_xxxxyz[k];

                g_x_0_z_0_xxz_xxxxz[k] = -g_x_0_z_0_xx_xxxxz[k] * ab_z + g_x_0_z_0_xx_xxxxzz[k];

                g_x_0_z_0_xxz_xxxyy[k] = -g_x_0_z_0_xx_xxxyy[k] * ab_z + g_x_0_z_0_xx_xxxyyz[k];

                g_x_0_z_0_xxz_xxxyz[k] = -g_x_0_z_0_xx_xxxyz[k] * ab_z + g_x_0_z_0_xx_xxxyzz[k];

                g_x_0_z_0_xxz_xxxzz[k] = -g_x_0_z_0_xx_xxxzz[k] * ab_z + g_x_0_z_0_xx_xxxzzz[k];

                g_x_0_z_0_xxz_xxyyy[k] = -g_x_0_z_0_xx_xxyyy[k] * ab_z + g_x_0_z_0_xx_xxyyyz[k];

                g_x_0_z_0_xxz_xxyyz[k] = -g_x_0_z_0_xx_xxyyz[k] * ab_z + g_x_0_z_0_xx_xxyyzz[k];

                g_x_0_z_0_xxz_xxyzz[k] = -g_x_0_z_0_xx_xxyzz[k] * ab_z + g_x_0_z_0_xx_xxyzzz[k];

                g_x_0_z_0_xxz_xxzzz[k] = -g_x_0_z_0_xx_xxzzz[k] * ab_z + g_x_0_z_0_xx_xxzzzz[k];

                g_x_0_z_0_xxz_xyyyy[k] = -g_x_0_z_0_xx_xyyyy[k] * ab_z + g_x_0_z_0_xx_xyyyyz[k];

                g_x_0_z_0_xxz_xyyyz[k] = -g_x_0_z_0_xx_xyyyz[k] * ab_z + g_x_0_z_0_xx_xyyyzz[k];

                g_x_0_z_0_xxz_xyyzz[k] = -g_x_0_z_0_xx_xyyzz[k] * ab_z + g_x_0_z_0_xx_xyyzzz[k];

                g_x_0_z_0_xxz_xyzzz[k] = -g_x_0_z_0_xx_xyzzz[k] * ab_z + g_x_0_z_0_xx_xyzzzz[k];

                g_x_0_z_0_xxz_xzzzz[k] = -g_x_0_z_0_xx_xzzzz[k] * ab_z + g_x_0_z_0_xx_xzzzzz[k];

                g_x_0_z_0_xxz_yyyyy[k] = -g_x_0_z_0_xx_yyyyy[k] * ab_z + g_x_0_z_0_xx_yyyyyz[k];

                g_x_0_z_0_xxz_yyyyz[k] = -g_x_0_z_0_xx_yyyyz[k] * ab_z + g_x_0_z_0_xx_yyyyzz[k];

                g_x_0_z_0_xxz_yyyzz[k] = -g_x_0_z_0_xx_yyyzz[k] * ab_z + g_x_0_z_0_xx_yyyzzz[k];

                g_x_0_z_0_xxz_yyzzz[k] = -g_x_0_z_0_xx_yyzzz[k] * ab_z + g_x_0_z_0_xx_yyzzzz[k];

                g_x_0_z_0_xxz_yzzzz[k] = -g_x_0_z_0_xx_yzzzz[k] * ab_z + g_x_0_z_0_xx_yzzzzz[k];

                g_x_0_z_0_xxz_zzzzz[k] = -g_x_0_z_0_xx_zzzzz[k] * ab_z + g_x_0_z_0_xx_zzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 483 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 484 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 485 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 486 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 487 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 488 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 489 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 490 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 491 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 492 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 493 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 494 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 495 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 496 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 497 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 498 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 499 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 500 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 501 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 502 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xy_xxxxx, g_x_0_z_0_xy_xxxxxy, g_x_0_z_0_xy_xxxxy, g_x_0_z_0_xy_xxxxyy, g_x_0_z_0_xy_xxxxyz, g_x_0_z_0_xy_xxxxz, g_x_0_z_0_xy_xxxyy, g_x_0_z_0_xy_xxxyyy, g_x_0_z_0_xy_xxxyyz, g_x_0_z_0_xy_xxxyz, g_x_0_z_0_xy_xxxyzz, g_x_0_z_0_xy_xxxzz, g_x_0_z_0_xy_xxyyy, g_x_0_z_0_xy_xxyyyy, g_x_0_z_0_xy_xxyyyz, g_x_0_z_0_xy_xxyyz, g_x_0_z_0_xy_xxyyzz, g_x_0_z_0_xy_xxyzz, g_x_0_z_0_xy_xxyzzz, g_x_0_z_0_xy_xxzzz, g_x_0_z_0_xy_xyyyy, g_x_0_z_0_xy_xyyyyy, g_x_0_z_0_xy_xyyyyz, g_x_0_z_0_xy_xyyyz, g_x_0_z_0_xy_xyyyzz, g_x_0_z_0_xy_xyyzz, g_x_0_z_0_xy_xyyzzz, g_x_0_z_0_xy_xyzzz, g_x_0_z_0_xy_xyzzzz, g_x_0_z_0_xy_xzzzz, g_x_0_z_0_xy_yyyyy, g_x_0_z_0_xy_yyyyyy, g_x_0_z_0_xy_yyyyyz, g_x_0_z_0_xy_yyyyz, g_x_0_z_0_xy_yyyyzz, g_x_0_z_0_xy_yyyzz, g_x_0_z_0_xy_yyyzzz, g_x_0_z_0_xy_yyzzz, g_x_0_z_0_xy_yyzzzz, g_x_0_z_0_xy_yzzzz, g_x_0_z_0_xy_yzzzzz, g_x_0_z_0_xy_zzzzz, g_x_0_z_0_xyy_xxxxx, g_x_0_z_0_xyy_xxxxy, g_x_0_z_0_xyy_xxxxz, g_x_0_z_0_xyy_xxxyy, g_x_0_z_0_xyy_xxxyz, g_x_0_z_0_xyy_xxxzz, g_x_0_z_0_xyy_xxyyy, g_x_0_z_0_xyy_xxyyz, g_x_0_z_0_xyy_xxyzz, g_x_0_z_0_xyy_xxzzz, g_x_0_z_0_xyy_xyyyy, g_x_0_z_0_xyy_xyyyz, g_x_0_z_0_xyy_xyyzz, g_x_0_z_0_xyy_xyzzz, g_x_0_z_0_xyy_xzzzz, g_x_0_z_0_xyy_yyyyy, g_x_0_z_0_xyy_yyyyz, g_x_0_z_0_xyy_yyyzz, g_x_0_z_0_xyy_yyzzz, g_x_0_z_0_xyy_yzzzz, g_x_0_z_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyy_xxxxx[k] = -g_x_0_z_0_xy_xxxxx[k] * ab_y + g_x_0_z_0_xy_xxxxxy[k];

                g_x_0_z_0_xyy_xxxxy[k] = -g_x_0_z_0_xy_xxxxy[k] * ab_y + g_x_0_z_0_xy_xxxxyy[k];

                g_x_0_z_0_xyy_xxxxz[k] = -g_x_0_z_0_xy_xxxxz[k] * ab_y + g_x_0_z_0_xy_xxxxyz[k];

                g_x_0_z_0_xyy_xxxyy[k] = -g_x_0_z_0_xy_xxxyy[k] * ab_y + g_x_0_z_0_xy_xxxyyy[k];

                g_x_0_z_0_xyy_xxxyz[k] = -g_x_0_z_0_xy_xxxyz[k] * ab_y + g_x_0_z_0_xy_xxxyyz[k];

                g_x_0_z_0_xyy_xxxzz[k] = -g_x_0_z_0_xy_xxxzz[k] * ab_y + g_x_0_z_0_xy_xxxyzz[k];

                g_x_0_z_0_xyy_xxyyy[k] = -g_x_0_z_0_xy_xxyyy[k] * ab_y + g_x_0_z_0_xy_xxyyyy[k];

                g_x_0_z_0_xyy_xxyyz[k] = -g_x_0_z_0_xy_xxyyz[k] * ab_y + g_x_0_z_0_xy_xxyyyz[k];

                g_x_0_z_0_xyy_xxyzz[k] = -g_x_0_z_0_xy_xxyzz[k] * ab_y + g_x_0_z_0_xy_xxyyzz[k];

                g_x_0_z_0_xyy_xxzzz[k] = -g_x_0_z_0_xy_xxzzz[k] * ab_y + g_x_0_z_0_xy_xxyzzz[k];

                g_x_0_z_0_xyy_xyyyy[k] = -g_x_0_z_0_xy_xyyyy[k] * ab_y + g_x_0_z_0_xy_xyyyyy[k];

                g_x_0_z_0_xyy_xyyyz[k] = -g_x_0_z_0_xy_xyyyz[k] * ab_y + g_x_0_z_0_xy_xyyyyz[k];

                g_x_0_z_0_xyy_xyyzz[k] = -g_x_0_z_0_xy_xyyzz[k] * ab_y + g_x_0_z_0_xy_xyyyzz[k];

                g_x_0_z_0_xyy_xyzzz[k] = -g_x_0_z_0_xy_xyzzz[k] * ab_y + g_x_0_z_0_xy_xyyzzz[k];

                g_x_0_z_0_xyy_xzzzz[k] = -g_x_0_z_0_xy_xzzzz[k] * ab_y + g_x_0_z_0_xy_xyzzzz[k];

                g_x_0_z_0_xyy_yyyyy[k] = -g_x_0_z_0_xy_yyyyy[k] * ab_y + g_x_0_z_0_xy_yyyyyy[k];

                g_x_0_z_0_xyy_yyyyz[k] = -g_x_0_z_0_xy_yyyyz[k] * ab_y + g_x_0_z_0_xy_yyyyyz[k];

                g_x_0_z_0_xyy_yyyzz[k] = -g_x_0_z_0_xy_yyyzz[k] * ab_y + g_x_0_z_0_xy_yyyyzz[k];

                g_x_0_z_0_xyy_yyzzz[k] = -g_x_0_z_0_xy_yyzzz[k] * ab_y + g_x_0_z_0_xy_yyyzzz[k];

                g_x_0_z_0_xyy_yzzzz[k] = -g_x_0_z_0_xy_yzzzz[k] * ab_y + g_x_0_z_0_xy_yyzzzz[k];

                g_x_0_z_0_xyy_zzzzz[k] = -g_x_0_z_0_xy_zzzzz[k] * ab_y + g_x_0_z_0_xy_yzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 504 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 505 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 506 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 507 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 508 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 509 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 510 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 511 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 512 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 513 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 514 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 515 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 516 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 517 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 518 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 519 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 520 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 521 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 522 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 523 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyz_xxxxx, g_x_0_z_0_xyz_xxxxy, g_x_0_z_0_xyz_xxxxz, g_x_0_z_0_xyz_xxxyy, g_x_0_z_0_xyz_xxxyz, g_x_0_z_0_xyz_xxxzz, g_x_0_z_0_xyz_xxyyy, g_x_0_z_0_xyz_xxyyz, g_x_0_z_0_xyz_xxyzz, g_x_0_z_0_xyz_xxzzz, g_x_0_z_0_xyz_xyyyy, g_x_0_z_0_xyz_xyyyz, g_x_0_z_0_xyz_xyyzz, g_x_0_z_0_xyz_xyzzz, g_x_0_z_0_xyz_xzzzz, g_x_0_z_0_xyz_yyyyy, g_x_0_z_0_xyz_yyyyz, g_x_0_z_0_xyz_yyyzz, g_x_0_z_0_xyz_yyzzz, g_x_0_z_0_xyz_yzzzz, g_x_0_z_0_xyz_zzzzz, g_x_0_z_0_xz_xxxxx, g_x_0_z_0_xz_xxxxxy, g_x_0_z_0_xz_xxxxy, g_x_0_z_0_xz_xxxxyy, g_x_0_z_0_xz_xxxxyz, g_x_0_z_0_xz_xxxxz, g_x_0_z_0_xz_xxxyy, g_x_0_z_0_xz_xxxyyy, g_x_0_z_0_xz_xxxyyz, g_x_0_z_0_xz_xxxyz, g_x_0_z_0_xz_xxxyzz, g_x_0_z_0_xz_xxxzz, g_x_0_z_0_xz_xxyyy, g_x_0_z_0_xz_xxyyyy, g_x_0_z_0_xz_xxyyyz, g_x_0_z_0_xz_xxyyz, g_x_0_z_0_xz_xxyyzz, g_x_0_z_0_xz_xxyzz, g_x_0_z_0_xz_xxyzzz, g_x_0_z_0_xz_xxzzz, g_x_0_z_0_xz_xyyyy, g_x_0_z_0_xz_xyyyyy, g_x_0_z_0_xz_xyyyyz, g_x_0_z_0_xz_xyyyz, g_x_0_z_0_xz_xyyyzz, g_x_0_z_0_xz_xyyzz, g_x_0_z_0_xz_xyyzzz, g_x_0_z_0_xz_xyzzz, g_x_0_z_0_xz_xyzzzz, g_x_0_z_0_xz_xzzzz, g_x_0_z_0_xz_yyyyy, g_x_0_z_0_xz_yyyyyy, g_x_0_z_0_xz_yyyyyz, g_x_0_z_0_xz_yyyyz, g_x_0_z_0_xz_yyyyzz, g_x_0_z_0_xz_yyyzz, g_x_0_z_0_xz_yyyzzz, g_x_0_z_0_xz_yyzzz, g_x_0_z_0_xz_yyzzzz, g_x_0_z_0_xz_yzzzz, g_x_0_z_0_xz_yzzzzz, g_x_0_z_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyz_xxxxx[k] = -g_x_0_z_0_xz_xxxxx[k] * ab_y + g_x_0_z_0_xz_xxxxxy[k];

                g_x_0_z_0_xyz_xxxxy[k] = -g_x_0_z_0_xz_xxxxy[k] * ab_y + g_x_0_z_0_xz_xxxxyy[k];

                g_x_0_z_0_xyz_xxxxz[k] = -g_x_0_z_0_xz_xxxxz[k] * ab_y + g_x_0_z_0_xz_xxxxyz[k];

                g_x_0_z_0_xyz_xxxyy[k] = -g_x_0_z_0_xz_xxxyy[k] * ab_y + g_x_0_z_0_xz_xxxyyy[k];

                g_x_0_z_0_xyz_xxxyz[k] = -g_x_0_z_0_xz_xxxyz[k] * ab_y + g_x_0_z_0_xz_xxxyyz[k];

                g_x_0_z_0_xyz_xxxzz[k] = -g_x_0_z_0_xz_xxxzz[k] * ab_y + g_x_0_z_0_xz_xxxyzz[k];

                g_x_0_z_0_xyz_xxyyy[k] = -g_x_0_z_0_xz_xxyyy[k] * ab_y + g_x_0_z_0_xz_xxyyyy[k];

                g_x_0_z_0_xyz_xxyyz[k] = -g_x_0_z_0_xz_xxyyz[k] * ab_y + g_x_0_z_0_xz_xxyyyz[k];

                g_x_0_z_0_xyz_xxyzz[k] = -g_x_0_z_0_xz_xxyzz[k] * ab_y + g_x_0_z_0_xz_xxyyzz[k];

                g_x_0_z_0_xyz_xxzzz[k] = -g_x_0_z_0_xz_xxzzz[k] * ab_y + g_x_0_z_0_xz_xxyzzz[k];

                g_x_0_z_0_xyz_xyyyy[k] = -g_x_0_z_0_xz_xyyyy[k] * ab_y + g_x_0_z_0_xz_xyyyyy[k];

                g_x_0_z_0_xyz_xyyyz[k] = -g_x_0_z_0_xz_xyyyz[k] * ab_y + g_x_0_z_0_xz_xyyyyz[k];

                g_x_0_z_0_xyz_xyyzz[k] = -g_x_0_z_0_xz_xyyzz[k] * ab_y + g_x_0_z_0_xz_xyyyzz[k];

                g_x_0_z_0_xyz_xyzzz[k] = -g_x_0_z_0_xz_xyzzz[k] * ab_y + g_x_0_z_0_xz_xyyzzz[k];

                g_x_0_z_0_xyz_xzzzz[k] = -g_x_0_z_0_xz_xzzzz[k] * ab_y + g_x_0_z_0_xz_xyzzzz[k];

                g_x_0_z_0_xyz_yyyyy[k] = -g_x_0_z_0_xz_yyyyy[k] * ab_y + g_x_0_z_0_xz_yyyyyy[k];

                g_x_0_z_0_xyz_yyyyz[k] = -g_x_0_z_0_xz_yyyyz[k] * ab_y + g_x_0_z_0_xz_yyyyyz[k];

                g_x_0_z_0_xyz_yyyzz[k] = -g_x_0_z_0_xz_yyyzz[k] * ab_y + g_x_0_z_0_xz_yyyyzz[k];

                g_x_0_z_0_xyz_yyzzz[k] = -g_x_0_z_0_xz_yyzzz[k] * ab_y + g_x_0_z_0_xz_yyyzzz[k];

                g_x_0_z_0_xyz_yzzzz[k] = -g_x_0_z_0_xz_yzzzz[k] * ab_y + g_x_0_z_0_xz_yyzzzz[k];

                g_x_0_z_0_xyz_zzzzz[k] = -g_x_0_z_0_xz_zzzzz[k] * ab_y + g_x_0_z_0_xz_yzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 525 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 526 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 527 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 528 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 529 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 530 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 531 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 532 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 533 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 534 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 535 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 536 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 537 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 538 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 539 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 540 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 541 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 542 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 543 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 544 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xz_xxxxx, g_x_0_z_0_xz_xxxxxz, g_x_0_z_0_xz_xxxxy, g_x_0_z_0_xz_xxxxyz, g_x_0_z_0_xz_xxxxz, g_x_0_z_0_xz_xxxxzz, g_x_0_z_0_xz_xxxyy, g_x_0_z_0_xz_xxxyyz, g_x_0_z_0_xz_xxxyz, g_x_0_z_0_xz_xxxyzz, g_x_0_z_0_xz_xxxzz, g_x_0_z_0_xz_xxxzzz, g_x_0_z_0_xz_xxyyy, g_x_0_z_0_xz_xxyyyz, g_x_0_z_0_xz_xxyyz, g_x_0_z_0_xz_xxyyzz, g_x_0_z_0_xz_xxyzz, g_x_0_z_0_xz_xxyzzz, g_x_0_z_0_xz_xxzzz, g_x_0_z_0_xz_xxzzzz, g_x_0_z_0_xz_xyyyy, g_x_0_z_0_xz_xyyyyz, g_x_0_z_0_xz_xyyyz, g_x_0_z_0_xz_xyyyzz, g_x_0_z_0_xz_xyyzz, g_x_0_z_0_xz_xyyzzz, g_x_0_z_0_xz_xyzzz, g_x_0_z_0_xz_xyzzzz, g_x_0_z_0_xz_xzzzz, g_x_0_z_0_xz_xzzzzz, g_x_0_z_0_xz_yyyyy, g_x_0_z_0_xz_yyyyyz, g_x_0_z_0_xz_yyyyz, g_x_0_z_0_xz_yyyyzz, g_x_0_z_0_xz_yyyzz, g_x_0_z_0_xz_yyyzzz, g_x_0_z_0_xz_yyzzz, g_x_0_z_0_xz_yyzzzz, g_x_0_z_0_xz_yzzzz, g_x_0_z_0_xz_yzzzzz, g_x_0_z_0_xz_zzzzz, g_x_0_z_0_xz_zzzzzz, g_x_0_z_0_xzz_xxxxx, g_x_0_z_0_xzz_xxxxy, g_x_0_z_0_xzz_xxxxz, g_x_0_z_0_xzz_xxxyy, g_x_0_z_0_xzz_xxxyz, g_x_0_z_0_xzz_xxxzz, g_x_0_z_0_xzz_xxyyy, g_x_0_z_0_xzz_xxyyz, g_x_0_z_0_xzz_xxyzz, g_x_0_z_0_xzz_xxzzz, g_x_0_z_0_xzz_xyyyy, g_x_0_z_0_xzz_xyyyz, g_x_0_z_0_xzz_xyyzz, g_x_0_z_0_xzz_xyzzz, g_x_0_z_0_xzz_xzzzz, g_x_0_z_0_xzz_yyyyy, g_x_0_z_0_xzz_yyyyz, g_x_0_z_0_xzz_yyyzz, g_x_0_z_0_xzz_yyzzz, g_x_0_z_0_xzz_yzzzz, g_x_0_z_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzz_xxxxx[k] = -g_x_0_z_0_xz_xxxxx[k] * ab_z + g_x_0_z_0_xz_xxxxxz[k];

                g_x_0_z_0_xzz_xxxxy[k] = -g_x_0_z_0_xz_xxxxy[k] * ab_z + g_x_0_z_0_xz_xxxxyz[k];

                g_x_0_z_0_xzz_xxxxz[k] = -g_x_0_z_0_xz_xxxxz[k] * ab_z + g_x_0_z_0_xz_xxxxzz[k];

                g_x_0_z_0_xzz_xxxyy[k] = -g_x_0_z_0_xz_xxxyy[k] * ab_z + g_x_0_z_0_xz_xxxyyz[k];

                g_x_0_z_0_xzz_xxxyz[k] = -g_x_0_z_0_xz_xxxyz[k] * ab_z + g_x_0_z_0_xz_xxxyzz[k];

                g_x_0_z_0_xzz_xxxzz[k] = -g_x_0_z_0_xz_xxxzz[k] * ab_z + g_x_0_z_0_xz_xxxzzz[k];

                g_x_0_z_0_xzz_xxyyy[k] = -g_x_0_z_0_xz_xxyyy[k] * ab_z + g_x_0_z_0_xz_xxyyyz[k];

                g_x_0_z_0_xzz_xxyyz[k] = -g_x_0_z_0_xz_xxyyz[k] * ab_z + g_x_0_z_0_xz_xxyyzz[k];

                g_x_0_z_0_xzz_xxyzz[k] = -g_x_0_z_0_xz_xxyzz[k] * ab_z + g_x_0_z_0_xz_xxyzzz[k];

                g_x_0_z_0_xzz_xxzzz[k] = -g_x_0_z_0_xz_xxzzz[k] * ab_z + g_x_0_z_0_xz_xxzzzz[k];

                g_x_0_z_0_xzz_xyyyy[k] = -g_x_0_z_0_xz_xyyyy[k] * ab_z + g_x_0_z_0_xz_xyyyyz[k];

                g_x_0_z_0_xzz_xyyyz[k] = -g_x_0_z_0_xz_xyyyz[k] * ab_z + g_x_0_z_0_xz_xyyyzz[k];

                g_x_0_z_0_xzz_xyyzz[k] = -g_x_0_z_0_xz_xyyzz[k] * ab_z + g_x_0_z_0_xz_xyyzzz[k];

                g_x_0_z_0_xzz_xyzzz[k] = -g_x_0_z_0_xz_xyzzz[k] * ab_z + g_x_0_z_0_xz_xyzzzz[k];

                g_x_0_z_0_xzz_xzzzz[k] = -g_x_0_z_0_xz_xzzzz[k] * ab_z + g_x_0_z_0_xz_xzzzzz[k];

                g_x_0_z_0_xzz_yyyyy[k] = -g_x_0_z_0_xz_yyyyy[k] * ab_z + g_x_0_z_0_xz_yyyyyz[k];

                g_x_0_z_0_xzz_yyyyz[k] = -g_x_0_z_0_xz_yyyyz[k] * ab_z + g_x_0_z_0_xz_yyyyzz[k];

                g_x_0_z_0_xzz_yyyzz[k] = -g_x_0_z_0_xz_yyyzz[k] * ab_z + g_x_0_z_0_xz_yyyzzz[k];

                g_x_0_z_0_xzz_yyzzz[k] = -g_x_0_z_0_xz_yyzzz[k] * ab_z + g_x_0_z_0_xz_yyzzzz[k];

                g_x_0_z_0_xzz_yzzzz[k] = -g_x_0_z_0_xz_yzzzz[k] * ab_z + g_x_0_z_0_xz_yzzzzz[k];

                g_x_0_z_0_xzz_zzzzz[k] = -g_x_0_z_0_xz_zzzzz[k] * ab_z + g_x_0_z_0_xz_zzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 546 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 547 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 548 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 549 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 550 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 551 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 552 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 553 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 554 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 555 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 556 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 557 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 558 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 559 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 560 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 561 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 562 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 563 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 564 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 565 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yy_xxxxx, g_x_0_z_0_yy_xxxxxy, g_x_0_z_0_yy_xxxxy, g_x_0_z_0_yy_xxxxyy, g_x_0_z_0_yy_xxxxyz, g_x_0_z_0_yy_xxxxz, g_x_0_z_0_yy_xxxyy, g_x_0_z_0_yy_xxxyyy, g_x_0_z_0_yy_xxxyyz, g_x_0_z_0_yy_xxxyz, g_x_0_z_0_yy_xxxyzz, g_x_0_z_0_yy_xxxzz, g_x_0_z_0_yy_xxyyy, g_x_0_z_0_yy_xxyyyy, g_x_0_z_0_yy_xxyyyz, g_x_0_z_0_yy_xxyyz, g_x_0_z_0_yy_xxyyzz, g_x_0_z_0_yy_xxyzz, g_x_0_z_0_yy_xxyzzz, g_x_0_z_0_yy_xxzzz, g_x_0_z_0_yy_xyyyy, g_x_0_z_0_yy_xyyyyy, g_x_0_z_0_yy_xyyyyz, g_x_0_z_0_yy_xyyyz, g_x_0_z_0_yy_xyyyzz, g_x_0_z_0_yy_xyyzz, g_x_0_z_0_yy_xyyzzz, g_x_0_z_0_yy_xyzzz, g_x_0_z_0_yy_xyzzzz, g_x_0_z_0_yy_xzzzz, g_x_0_z_0_yy_yyyyy, g_x_0_z_0_yy_yyyyyy, g_x_0_z_0_yy_yyyyyz, g_x_0_z_0_yy_yyyyz, g_x_0_z_0_yy_yyyyzz, g_x_0_z_0_yy_yyyzz, g_x_0_z_0_yy_yyyzzz, g_x_0_z_0_yy_yyzzz, g_x_0_z_0_yy_yyzzzz, g_x_0_z_0_yy_yzzzz, g_x_0_z_0_yy_yzzzzz, g_x_0_z_0_yy_zzzzz, g_x_0_z_0_yyy_xxxxx, g_x_0_z_0_yyy_xxxxy, g_x_0_z_0_yyy_xxxxz, g_x_0_z_0_yyy_xxxyy, g_x_0_z_0_yyy_xxxyz, g_x_0_z_0_yyy_xxxzz, g_x_0_z_0_yyy_xxyyy, g_x_0_z_0_yyy_xxyyz, g_x_0_z_0_yyy_xxyzz, g_x_0_z_0_yyy_xxzzz, g_x_0_z_0_yyy_xyyyy, g_x_0_z_0_yyy_xyyyz, g_x_0_z_0_yyy_xyyzz, g_x_0_z_0_yyy_xyzzz, g_x_0_z_0_yyy_xzzzz, g_x_0_z_0_yyy_yyyyy, g_x_0_z_0_yyy_yyyyz, g_x_0_z_0_yyy_yyyzz, g_x_0_z_0_yyy_yyzzz, g_x_0_z_0_yyy_yzzzz, g_x_0_z_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyy_xxxxx[k] = -g_x_0_z_0_yy_xxxxx[k] * ab_y + g_x_0_z_0_yy_xxxxxy[k];

                g_x_0_z_0_yyy_xxxxy[k] = -g_x_0_z_0_yy_xxxxy[k] * ab_y + g_x_0_z_0_yy_xxxxyy[k];

                g_x_0_z_0_yyy_xxxxz[k] = -g_x_0_z_0_yy_xxxxz[k] * ab_y + g_x_0_z_0_yy_xxxxyz[k];

                g_x_0_z_0_yyy_xxxyy[k] = -g_x_0_z_0_yy_xxxyy[k] * ab_y + g_x_0_z_0_yy_xxxyyy[k];

                g_x_0_z_0_yyy_xxxyz[k] = -g_x_0_z_0_yy_xxxyz[k] * ab_y + g_x_0_z_0_yy_xxxyyz[k];

                g_x_0_z_0_yyy_xxxzz[k] = -g_x_0_z_0_yy_xxxzz[k] * ab_y + g_x_0_z_0_yy_xxxyzz[k];

                g_x_0_z_0_yyy_xxyyy[k] = -g_x_0_z_0_yy_xxyyy[k] * ab_y + g_x_0_z_0_yy_xxyyyy[k];

                g_x_0_z_0_yyy_xxyyz[k] = -g_x_0_z_0_yy_xxyyz[k] * ab_y + g_x_0_z_0_yy_xxyyyz[k];

                g_x_0_z_0_yyy_xxyzz[k] = -g_x_0_z_0_yy_xxyzz[k] * ab_y + g_x_0_z_0_yy_xxyyzz[k];

                g_x_0_z_0_yyy_xxzzz[k] = -g_x_0_z_0_yy_xxzzz[k] * ab_y + g_x_0_z_0_yy_xxyzzz[k];

                g_x_0_z_0_yyy_xyyyy[k] = -g_x_0_z_0_yy_xyyyy[k] * ab_y + g_x_0_z_0_yy_xyyyyy[k];

                g_x_0_z_0_yyy_xyyyz[k] = -g_x_0_z_0_yy_xyyyz[k] * ab_y + g_x_0_z_0_yy_xyyyyz[k];

                g_x_0_z_0_yyy_xyyzz[k] = -g_x_0_z_0_yy_xyyzz[k] * ab_y + g_x_0_z_0_yy_xyyyzz[k];

                g_x_0_z_0_yyy_xyzzz[k] = -g_x_0_z_0_yy_xyzzz[k] * ab_y + g_x_0_z_0_yy_xyyzzz[k];

                g_x_0_z_0_yyy_xzzzz[k] = -g_x_0_z_0_yy_xzzzz[k] * ab_y + g_x_0_z_0_yy_xyzzzz[k];

                g_x_0_z_0_yyy_yyyyy[k] = -g_x_0_z_0_yy_yyyyy[k] * ab_y + g_x_0_z_0_yy_yyyyyy[k];

                g_x_0_z_0_yyy_yyyyz[k] = -g_x_0_z_0_yy_yyyyz[k] * ab_y + g_x_0_z_0_yy_yyyyyz[k];

                g_x_0_z_0_yyy_yyyzz[k] = -g_x_0_z_0_yy_yyyzz[k] * ab_y + g_x_0_z_0_yy_yyyyzz[k];

                g_x_0_z_0_yyy_yyzzz[k] = -g_x_0_z_0_yy_yyzzz[k] * ab_y + g_x_0_z_0_yy_yyyzzz[k];

                g_x_0_z_0_yyy_yzzzz[k] = -g_x_0_z_0_yy_yzzzz[k] * ab_y + g_x_0_z_0_yy_yyzzzz[k];

                g_x_0_z_0_yyy_zzzzz[k] = -g_x_0_z_0_yy_zzzzz[k] * ab_y + g_x_0_z_0_yy_yzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 567 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 568 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 569 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 570 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 571 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 572 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 573 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 574 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 575 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 576 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 577 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 578 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 579 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 580 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 581 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 582 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 583 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 584 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 585 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 586 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyz_xxxxx, g_x_0_z_0_yyz_xxxxy, g_x_0_z_0_yyz_xxxxz, g_x_0_z_0_yyz_xxxyy, g_x_0_z_0_yyz_xxxyz, g_x_0_z_0_yyz_xxxzz, g_x_0_z_0_yyz_xxyyy, g_x_0_z_0_yyz_xxyyz, g_x_0_z_0_yyz_xxyzz, g_x_0_z_0_yyz_xxzzz, g_x_0_z_0_yyz_xyyyy, g_x_0_z_0_yyz_xyyyz, g_x_0_z_0_yyz_xyyzz, g_x_0_z_0_yyz_xyzzz, g_x_0_z_0_yyz_xzzzz, g_x_0_z_0_yyz_yyyyy, g_x_0_z_0_yyz_yyyyz, g_x_0_z_0_yyz_yyyzz, g_x_0_z_0_yyz_yyzzz, g_x_0_z_0_yyz_yzzzz, g_x_0_z_0_yyz_zzzzz, g_x_0_z_0_yz_xxxxx, g_x_0_z_0_yz_xxxxxy, g_x_0_z_0_yz_xxxxy, g_x_0_z_0_yz_xxxxyy, g_x_0_z_0_yz_xxxxyz, g_x_0_z_0_yz_xxxxz, g_x_0_z_0_yz_xxxyy, g_x_0_z_0_yz_xxxyyy, g_x_0_z_0_yz_xxxyyz, g_x_0_z_0_yz_xxxyz, g_x_0_z_0_yz_xxxyzz, g_x_0_z_0_yz_xxxzz, g_x_0_z_0_yz_xxyyy, g_x_0_z_0_yz_xxyyyy, g_x_0_z_0_yz_xxyyyz, g_x_0_z_0_yz_xxyyz, g_x_0_z_0_yz_xxyyzz, g_x_0_z_0_yz_xxyzz, g_x_0_z_0_yz_xxyzzz, g_x_0_z_0_yz_xxzzz, g_x_0_z_0_yz_xyyyy, g_x_0_z_0_yz_xyyyyy, g_x_0_z_0_yz_xyyyyz, g_x_0_z_0_yz_xyyyz, g_x_0_z_0_yz_xyyyzz, g_x_0_z_0_yz_xyyzz, g_x_0_z_0_yz_xyyzzz, g_x_0_z_0_yz_xyzzz, g_x_0_z_0_yz_xyzzzz, g_x_0_z_0_yz_xzzzz, g_x_0_z_0_yz_yyyyy, g_x_0_z_0_yz_yyyyyy, g_x_0_z_0_yz_yyyyyz, g_x_0_z_0_yz_yyyyz, g_x_0_z_0_yz_yyyyzz, g_x_0_z_0_yz_yyyzz, g_x_0_z_0_yz_yyyzzz, g_x_0_z_0_yz_yyzzz, g_x_0_z_0_yz_yyzzzz, g_x_0_z_0_yz_yzzzz, g_x_0_z_0_yz_yzzzzz, g_x_0_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyz_xxxxx[k] = -g_x_0_z_0_yz_xxxxx[k] * ab_y + g_x_0_z_0_yz_xxxxxy[k];

                g_x_0_z_0_yyz_xxxxy[k] = -g_x_0_z_0_yz_xxxxy[k] * ab_y + g_x_0_z_0_yz_xxxxyy[k];

                g_x_0_z_0_yyz_xxxxz[k] = -g_x_0_z_0_yz_xxxxz[k] * ab_y + g_x_0_z_0_yz_xxxxyz[k];

                g_x_0_z_0_yyz_xxxyy[k] = -g_x_0_z_0_yz_xxxyy[k] * ab_y + g_x_0_z_0_yz_xxxyyy[k];

                g_x_0_z_0_yyz_xxxyz[k] = -g_x_0_z_0_yz_xxxyz[k] * ab_y + g_x_0_z_0_yz_xxxyyz[k];

                g_x_0_z_0_yyz_xxxzz[k] = -g_x_0_z_0_yz_xxxzz[k] * ab_y + g_x_0_z_0_yz_xxxyzz[k];

                g_x_0_z_0_yyz_xxyyy[k] = -g_x_0_z_0_yz_xxyyy[k] * ab_y + g_x_0_z_0_yz_xxyyyy[k];

                g_x_0_z_0_yyz_xxyyz[k] = -g_x_0_z_0_yz_xxyyz[k] * ab_y + g_x_0_z_0_yz_xxyyyz[k];

                g_x_0_z_0_yyz_xxyzz[k] = -g_x_0_z_0_yz_xxyzz[k] * ab_y + g_x_0_z_0_yz_xxyyzz[k];

                g_x_0_z_0_yyz_xxzzz[k] = -g_x_0_z_0_yz_xxzzz[k] * ab_y + g_x_0_z_0_yz_xxyzzz[k];

                g_x_0_z_0_yyz_xyyyy[k] = -g_x_0_z_0_yz_xyyyy[k] * ab_y + g_x_0_z_0_yz_xyyyyy[k];

                g_x_0_z_0_yyz_xyyyz[k] = -g_x_0_z_0_yz_xyyyz[k] * ab_y + g_x_0_z_0_yz_xyyyyz[k];

                g_x_0_z_0_yyz_xyyzz[k] = -g_x_0_z_0_yz_xyyzz[k] * ab_y + g_x_0_z_0_yz_xyyyzz[k];

                g_x_0_z_0_yyz_xyzzz[k] = -g_x_0_z_0_yz_xyzzz[k] * ab_y + g_x_0_z_0_yz_xyyzzz[k];

                g_x_0_z_0_yyz_xzzzz[k] = -g_x_0_z_0_yz_xzzzz[k] * ab_y + g_x_0_z_0_yz_xyzzzz[k];

                g_x_0_z_0_yyz_yyyyy[k] = -g_x_0_z_0_yz_yyyyy[k] * ab_y + g_x_0_z_0_yz_yyyyyy[k];

                g_x_0_z_0_yyz_yyyyz[k] = -g_x_0_z_0_yz_yyyyz[k] * ab_y + g_x_0_z_0_yz_yyyyyz[k];

                g_x_0_z_0_yyz_yyyzz[k] = -g_x_0_z_0_yz_yyyzz[k] * ab_y + g_x_0_z_0_yz_yyyyzz[k];

                g_x_0_z_0_yyz_yyzzz[k] = -g_x_0_z_0_yz_yyzzz[k] * ab_y + g_x_0_z_0_yz_yyyzzz[k];

                g_x_0_z_0_yyz_yzzzz[k] = -g_x_0_z_0_yz_yzzzz[k] * ab_y + g_x_0_z_0_yz_yyzzzz[k];

                g_x_0_z_0_yyz_zzzzz[k] = -g_x_0_z_0_yz_zzzzz[k] * ab_y + g_x_0_z_0_yz_yzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 588 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 589 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 590 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 591 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 592 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 593 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 594 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 595 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 596 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 597 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 598 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 599 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 600 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 601 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 602 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 603 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 604 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 605 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 606 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 607 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzz_xxxxx, g_x_0_z_0_yzz_xxxxy, g_x_0_z_0_yzz_xxxxz, g_x_0_z_0_yzz_xxxyy, g_x_0_z_0_yzz_xxxyz, g_x_0_z_0_yzz_xxxzz, g_x_0_z_0_yzz_xxyyy, g_x_0_z_0_yzz_xxyyz, g_x_0_z_0_yzz_xxyzz, g_x_0_z_0_yzz_xxzzz, g_x_0_z_0_yzz_xyyyy, g_x_0_z_0_yzz_xyyyz, g_x_0_z_0_yzz_xyyzz, g_x_0_z_0_yzz_xyzzz, g_x_0_z_0_yzz_xzzzz, g_x_0_z_0_yzz_yyyyy, g_x_0_z_0_yzz_yyyyz, g_x_0_z_0_yzz_yyyzz, g_x_0_z_0_yzz_yyzzz, g_x_0_z_0_yzz_yzzzz, g_x_0_z_0_yzz_zzzzz, g_x_0_z_0_zz_xxxxx, g_x_0_z_0_zz_xxxxxy, g_x_0_z_0_zz_xxxxy, g_x_0_z_0_zz_xxxxyy, g_x_0_z_0_zz_xxxxyz, g_x_0_z_0_zz_xxxxz, g_x_0_z_0_zz_xxxyy, g_x_0_z_0_zz_xxxyyy, g_x_0_z_0_zz_xxxyyz, g_x_0_z_0_zz_xxxyz, g_x_0_z_0_zz_xxxyzz, g_x_0_z_0_zz_xxxzz, g_x_0_z_0_zz_xxyyy, g_x_0_z_0_zz_xxyyyy, g_x_0_z_0_zz_xxyyyz, g_x_0_z_0_zz_xxyyz, g_x_0_z_0_zz_xxyyzz, g_x_0_z_0_zz_xxyzz, g_x_0_z_0_zz_xxyzzz, g_x_0_z_0_zz_xxzzz, g_x_0_z_0_zz_xyyyy, g_x_0_z_0_zz_xyyyyy, g_x_0_z_0_zz_xyyyyz, g_x_0_z_0_zz_xyyyz, g_x_0_z_0_zz_xyyyzz, g_x_0_z_0_zz_xyyzz, g_x_0_z_0_zz_xyyzzz, g_x_0_z_0_zz_xyzzz, g_x_0_z_0_zz_xyzzzz, g_x_0_z_0_zz_xzzzz, g_x_0_z_0_zz_yyyyy, g_x_0_z_0_zz_yyyyyy, g_x_0_z_0_zz_yyyyyz, g_x_0_z_0_zz_yyyyz, g_x_0_z_0_zz_yyyyzz, g_x_0_z_0_zz_yyyzz, g_x_0_z_0_zz_yyyzzz, g_x_0_z_0_zz_yyzzz, g_x_0_z_0_zz_yyzzzz, g_x_0_z_0_zz_yzzzz, g_x_0_z_0_zz_yzzzzz, g_x_0_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzz_xxxxx[k] = -g_x_0_z_0_zz_xxxxx[k] * ab_y + g_x_0_z_0_zz_xxxxxy[k];

                g_x_0_z_0_yzz_xxxxy[k] = -g_x_0_z_0_zz_xxxxy[k] * ab_y + g_x_0_z_0_zz_xxxxyy[k];

                g_x_0_z_0_yzz_xxxxz[k] = -g_x_0_z_0_zz_xxxxz[k] * ab_y + g_x_0_z_0_zz_xxxxyz[k];

                g_x_0_z_0_yzz_xxxyy[k] = -g_x_0_z_0_zz_xxxyy[k] * ab_y + g_x_0_z_0_zz_xxxyyy[k];

                g_x_0_z_0_yzz_xxxyz[k] = -g_x_0_z_0_zz_xxxyz[k] * ab_y + g_x_0_z_0_zz_xxxyyz[k];

                g_x_0_z_0_yzz_xxxzz[k] = -g_x_0_z_0_zz_xxxzz[k] * ab_y + g_x_0_z_0_zz_xxxyzz[k];

                g_x_0_z_0_yzz_xxyyy[k] = -g_x_0_z_0_zz_xxyyy[k] * ab_y + g_x_0_z_0_zz_xxyyyy[k];

                g_x_0_z_0_yzz_xxyyz[k] = -g_x_0_z_0_zz_xxyyz[k] * ab_y + g_x_0_z_0_zz_xxyyyz[k];

                g_x_0_z_0_yzz_xxyzz[k] = -g_x_0_z_0_zz_xxyzz[k] * ab_y + g_x_0_z_0_zz_xxyyzz[k];

                g_x_0_z_0_yzz_xxzzz[k] = -g_x_0_z_0_zz_xxzzz[k] * ab_y + g_x_0_z_0_zz_xxyzzz[k];

                g_x_0_z_0_yzz_xyyyy[k] = -g_x_0_z_0_zz_xyyyy[k] * ab_y + g_x_0_z_0_zz_xyyyyy[k];

                g_x_0_z_0_yzz_xyyyz[k] = -g_x_0_z_0_zz_xyyyz[k] * ab_y + g_x_0_z_0_zz_xyyyyz[k];

                g_x_0_z_0_yzz_xyyzz[k] = -g_x_0_z_0_zz_xyyzz[k] * ab_y + g_x_0_z_0_zz_xyyyzz[k];

                g_x_0_z_0_yzz_xyzzz[k] = -g_x_0_z_0_zz_xyzzz[k] * ab_y + g_x_0_z_0_zz_xyyzzz[k];

                g_x_0_z_0_yzz_xzzzz[k] = -g_x_0_z_0_zz_xzzzz[k] * ab_y + g_x_0_z_0_zz_xyzzzz[k];

                g_x_0_z_0_yzz_yyyyy[k] = -g_x_0_z_0_zz_yyyyy[k] * ab_y + g_x_0_z_0_zz_yyyyyy[k];

                g_x_0_z_0_yzz_yyyyz[k] = -g_x_0_z_0_zz_yyyyz[k] * ab_y + g_x_0_z_0_zz_yyyyyz[k];

                g_x_0_z_0_yzz_yyyzz[k] = -g_x_0_z_0_zz_yyyzz[k] * ab_y + g_x_0_z_0_zz_yyyyzz[k];

                g_x_0_z_0_yzz_yyzzz[k] = -g_x_0_z_0_zz_yyzzz[k] * ab_y + g_x_0_z_0_zz_yyyzzz[k];

                g_x_0_z_0_yzz_yzzzz[k] = -g_x_0_z_0_zz_yzzzz[k] * ab_y + g_x_0_z_0_zz_yyzzzz[k];

                g_x_0_z_0_yzz_zzzzz[k] = -g_x_0_z_0_zz_zzzzz[k] * ab_y + g_x_0_z_0_zz_yzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 609 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 610 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 611 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 612 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 613 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 614 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 615 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 616 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 617 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 618 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 619 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 620 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 621 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 622 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 623 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 624 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 625 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 626 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 627 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 628 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zz_xxxxx, g_x_0_z_0_zz_xxxxxz, g_x_0_z_0_zz_xxxxy, g_x_0_z_0_zz_xxxxyz, g_x_0_z_0_zz_xxxxz, g_x_0_z_0_zz_xxxxzz, g_x_0_z_0_zz_xxxyy, g_x_0_z_0_zz_xxxyyz, g_x_0_z_0_zz_xxxyz, g_x_0_z_0_zz_xxxyzz, g_x_0_z_0_zz_xxxzz, g_x_0_z_0_zz_xxxzzz, g_x_0_z_0_zz_xxyyy, g_x_0_z_0_zz_xxyyyz, g_x_0_z_0_zz_xxyyz, g_x_0_z_0_zz_xxyyzz, g_x_0_z_0_zz_xxyzz, g_x_0_z_0_zz_xxyzzz, g_x_0_z_0_zz_xxzzz, g_x_0_z_0_zz_xxzzzz, g_x_0_z_0_zz_xyyyy, g_x_0_z_0_zz_xyyyyz, g_x_0_z_0_zz_xyyyz, g_x_0_z_0_zz_xyyyzz, g_x_0_z_0_zz_xyyzz, g_x_0_z_0_zz_xyyzzz, g_x_0_z_0_zz_xyzzz, g_x_0_z_0_zz_xyzzzz, g_x_0_z_0_zz_xzzzz, g_x_0_z_0_zz_xzzzzz, g_x_0_z_0_zz_yyyyy, g_x_0_z_0_zz_yyyyyz, g_x_0_z_0_zz_yyyyz, g_x_0_z_0_zz_yyyyzz, g_x_0_z_0_zz_yyyzz, g_x_0_z_0_zz_yyyzzz, g_x_0_z_0_zz_yyzzz, g_x_0_z_0_zz_yyzzzz, g_x_0_z_0_zz_yzzzz, g_x_0_z_0_zz_yzzzzz, g_x_0_z_0_zz_zzzzz, g_x_0_z_0_zz_zzzzzz, g_x_0_z_0_zzz_xxxxx, g_x_0_z_0_zzz_xxxxy, g_x_0_z_0_zzz_xxxxz, g_x_0_z_0_zzz_xxxyy, g_x_0_z_0_zzz_xxxyz, g_x_0_z_0_zzz_xxxzz, g_x_0_z_0_zzz_xxyyy, g_x_0_z_0_zzz_xxyyz, g_x_0_z_0_zzz_xxyzz, g_x_0_z_0_zzz_xxzzz, g_x_0_z_0_zzz_xyyyy, g_x_0_z_0_zzz_xyyyz, g_x_0_z_0_zzz_xyyzz, g_x_0_z_0_zzz_xyzzz, g_x_0_z_0_zzz_xzzzz, g_x_0_z_0_zzz_yyyyy, g_x_0_z_0_zzz_yyyyz, g_x_0_z_0_zzz_yyyzz, g_x_0_z_0_zzz_yyzzz, g_x_0_z_0_zzz_yzzzz, g_x_0_z_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzz_xxxxx[k] = -g_x_0_z_0_zz_xxxxx[k] * ab_z + g_x_0_z_0_zz_xxxxxz[k];

                g_x_0_z_0_zzz_xxxxy[k] = -g_x_0_z_0_zz_xxxxy[k] * ab_z + g_x_0_z_0_zz_xxxxyz[k];

                g_x_0_z_0_zzz_xxxxz[k] = -g_x_0_z_0_zz_xxxxz[k] * ab_z + g_x_0_z_0_zz_xxxxzz[k];

                g_x_0_z_0_zzz_xxxyy[k] = -g_x_0_z_0_zz_xxxyy[k] * ab_z + g_x_0_z_0_zz_xxxyyz[k];

                g_x_0_z_0_zzz_xxxyz[k] = -g_x_0_z_0_zz_xxxyz[k] * ab_z + g_x_0_z_0_zz_xxxyzz[k];

                g_x_0_z_0_zzz_xxxzz[k] = -g_x_0_z_0_zz_xxxzz[k] * ab_z + g_x_0_z_0_zz_xxxzzz[k];

                g_x_0_z_0_zzz_xxyyy[k] = -g_x_0_z_0_zz_xxyyy[k] * ab_z + g_x_0_z_0_zz_xxyyyz[k];

                g_x_0_z_0_zzz_xxyyz[k] = -g_x_0_z_0_zz_xxyyz[k] * ab_z + g_x_0_z_0_zz_xxyyzz[k];

                g_x_0_z_0_zzz_xxyzz[k] = -g_x_0_z_0_zz_xxyzz[k] * ab_z + g_x_0_z_0_zz_xxyzzz[k];

                g_x_0_z_0_zzz_xxzzz[k] = -g_x_0_z_0_zz_xxzzz[k] * ab_z + g_x_0_z_0_zz_xxzzzz[k];

                g_x_0_z_0_zzz_xyyyy[k] = -g_x_0_z_0_zz_xyyyy[k] * ab_z + g_x_0_z_0_zz_xyyyyz[k];

                g_x_0_z_0_zzz_xyyyz[k] = -g_x_0_z_0_zz_xyyyz[k] * ab_z + g_x_0_z_0_zz_xyyyzz[k];

                g_x_0_z_0_zzz_xyyzz[k] = -g_x_0_z_0_zz_xyyzz[k] * ab_z + g_x_0_z_0_zz_xyyzzz[k];

                g_x_0_z_0_zzz_xyzzz[k] = -g_x_0_z_0_zz_xyzzz[k] * ab_z + g_x_0_z_0_zz_xyzzzz[k];

                g_x_0_z_0_zzz_xzzzz[k] = -g_x_0_z_0_zz_xzzzz[k] * ab_z + g_x_0_z_0_zz_xzzzzz[k];

                g_x_0_z_0_zzz_yyyyy[k] = -g_x_0_z_0_zz_yyyyy[k] * ab_z + g_x_0_z_0_zz_yyyyyz[k];

                g_x_0_z_0_zzz_yyyyz[k] = -g_x_0_z_0_zz_yyyyz[k] * ab_z + g_x_0_z_0_zz_yyyyzz[k];

                g_x_0_z_0_zzz_yyyzz[k] = -g_x_0_z_0_zz_yyyzz[k] * ab_z + g_x_0_z_0_zz_yyyzzz[k];

                g_x_0_z_0_zzz_yyzzz[k] = -g_x_0_z_0_zz_yyzzz[k] * ab_z + g_x_0_z_0_zz_yyzzzz[k];

                g_x_0_z_0_zzz_yzzzz[k] = -g_x_0_z_0_zz_yzzzz[k] * ab_z + g_x_0_z_0_zz_yzzzzz[k];

                g_x_0_z_0_zzz_zzzzz[k] = -g_x_0_z_0_zz_zzzzz[k] * ab_z + g_x_0_z_0_zz_zzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 630 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 631 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 632 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 633 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 634 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 635 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 636 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 637 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 638 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 639 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 640 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 641 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 642 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 643 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 644 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 645 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 646 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 647 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 648 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 649 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xx_xxxxx, g_y_0_x_0_xx_xxxxxx, g_y_0_x_0_xx_xxxxxy, g_y_0_x_0_xx_xxxxxz, g_y_0_x_0_xx_xxxxy, g_y_0_x_0_xx_xxxxyy, g_y_0_x_0_xx_xxxxyz, g_y_0_x_0_xx_xxxxz, g_y_0_x_0_xx_xxxxzz, g_y_0_x_0_xx_xxxyy, g_y_0_x_0_xx_xxxyyy, g_y_0_x_0_xx_xxxyyz, g_y_0_x_0_xx_xxxyz, g_y_0_x_0_xx_xxxyzz, g_y_0_x_0_xx_xxxzz, g_y_0_x_0_xx_xxxzzz, g_y_0_x_0_xx_xxyyy, g_y_0_x_0_xx_xxyyyy, g_y_0_x_0_xx_xxyyyz, g_y_0_x_0_xx_xxyyz, g_y_0_x_0_xx_xxyyzz, g_y_0_x_0_xx_xxyzz, g_y_0_x_0_xx_xxyzzz, g_y_0_x_0_xx_xxzzz, g_y_0_x_0_xx_xxzzzz, g_y_0_x_0_xx_xyyyy, g_y_0_x_0_xx_xyyyyy, g_y_0_x_0_xx_xyyyyz, g_y_0_x_0_xx_xyyyz, g_y_0_x_0_xx_xyyyzz, g_y_0_x_0_xx_xyyzz, g_y_0_x_0_xx_xyyzzz, g_y_0_x_0_xx_xyzzz, g_y_0_x_0_xx_xyzzzz, g_y_0_x_0_xx_xzzzz, g_y_0_x_0_xx_xzzzzz, g_y_0_x_0_xx_yyyyy, g_y_0_x_0_xx_yyyyz, g_y_0_x_0_xx_yyyzz, g_y_0_x_0_xx_yyzzz, g_y_0_x_0_xx_yzzzz, g_y_0_x_0_xx_zzzzz, g_y_0_x_0_xxx_xxxxx, g_y_0_x_0_xxx_xxxxy, g_y_0_x_0_xxx_xxxxz, g_y_0_x_0_xxx_xxxyy, g_y_0_x_0_xxx_xxxyz, g_y_0_x_0_xxx_xxxzz, g_y_0_x_0_xxx_xxyyy, g_y_0_x_0_xxx_xxyyz, g_y_0_x_0_xxx_xxyzz, g_y_0_x_0_xxx_xxzzz, g_y_0_x_0_xxx_xyyyy, g_y_0_x_0_xxx_xyyyz, g_y_0_x_0_xxx_xyyzz, g_y_0_x_0_xxx_xyzzz, g_y_0_x_0_xxx_xzzzz, g_y_0_x_0_xxx_yyyyy, g_y_0_x_0_xxx_yyyyz, g_y_0_x_0_xxx_yyyzz, g_y_0_x_0_xxx_yyzzz, g_y_0_x_0_xxx_yzzzz, g_y_0_x_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxx_xxxxx[k] = -g_y_0_x_0_xx_xxxxx[k] * ab_x + g_y_0_x_0_xx_xxxxxx[k];

                g_y_0_x_0_xxx_xxxxy[k] = -g_y_0_x_0_xx_xxxxy[k] * ab_x + g_y_0_x_0_xx_xxxxxy[k];

                g_y_0_x_0_xxx_xxxxz[k] = -g_y_0_x_0_xx_xxxxz[k] * ab_x + g_y_0_x_0_xx_xxxxxz[k];

                g_y_0_x_0_xxx_xxxyy[k] = -g_y_0_x_0_xx_xxxyy[k] * ab_x + g_y_0_x_0_xx_xxxxyy[k];

                g_y_0_x_0_xxx_xxxyz[k] = -g_y_0_x_0_xx_xxxyz[k] * ab_x + g_y_0_x_0_xx_xxxxyz[k];

                g_y_0_x_0_xxx_xxxzz[k] = -g_y_0_x_0_xx_xxxzz[k] * ab_x + g_y_0_x_0_xx_xxxxzz[k];

                g_y_0_x_0_xxx_xxyyy[k] = -g_y_0_x_0_xx_xxyyy[k] * ab_x + g_y_0_x_0_xx_xxxyyy[k];

                g_y_0_x_0_xxx_xxyyz[k] = -g_y_0_x_0_xx_xxyyz[k] * ab_x + g_y_0_x_0_xx_xxxyyz[k];

                g_y_0_x_0_xxx_xxyzz[k] = -g_y_0_x_0_xx_xxyzz[k] * ab_x + g_y_0_x_0_xx_xxxyzz[k];

                g_y_0_x_0_xxx_xxzzz[k] = -g_y_0_x_0_xx_xxzzz[k] * ab_x + g_y_0_x_0_xx_xxxzzz[k];

                g_y_0_x_0_xxx_xyyyy[k] = -g_y_0_x_0_xx_xyyyy[k] * ab_x + g_y_0_x_0_xx_xxyyyy[k];

                g_y_0_x_0_xxx_xyyyz[k] = -g_y_0_x_0_xx_xyyyz[k] * ab_x + g_y_0_x_0_xx_xxyyyz[k];

                g_y_0_x_0_xxx_xyyzz[k] = -g_y_0_x_0_xx_xyyzz[k] * ab_x + g_y_0_x_0_xx_xxyyzz[k];

                g_y_0_x_0_xxx_xyzzz[k] = -g_y_0_x_0_xx_xyzzz[k] * ab_x + g_y_0_x_0_xx_xxyzzz[k];

                g_y_0_x_0_xxx_xzzzz[k] = -g_y_0_x_0_xx_xzzzz[k] * ab_x + g_y_0_x_0_xx_xxzzzz[k];

                g_y_0_x_0_xxx_yyyyy[k] = -g_y_0_x_0_xx_yyyyy[k] * ab_x + g_y_0_x_0_xx_xyyyyy[k];

                g_y_0_x_0_xxx_yyyyz[k] = -g_y_0_x_0_xx_yyyyz[k] * ab_x + g_y_0_x_0_xx_xyyyyz[k];

                g_y_0_x_0_xxx_yyyzz[k] = -g_y_0_x_0_xx_yyyzz[k] * ab_x + g_y_0_x_0_xx_xyyyzz[k];

                g_y_0_x_0_xxx_yyzzz[k] = -g_y_0_x_0_xx_yyzzz[k] * ab_x + g_y_0_x_0_xx_xyyzzz[k];

                g_y_0_x_0_xxx_yzzzz[k] = -g_y_0_x_0_xx_yzzzz[k] * ab_x + g_y_0_x_0_xx_xyzzzz[k];

                g_y_0_x_0_xxx_zzzzz[k] = -g_y_0_x_0_xx_zzzzz[k] * ab_x + g_y_0_x_0_xx_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 651 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 652 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 653 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 654 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 655 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 656 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 657 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 658 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 659 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 660 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 661 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 662 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 663 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 664 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 665 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 666 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 667 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 668 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 669 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 670 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxy_xxxxx, g_y_0_x_0_xxy_xxxxy, g_y_0_x_0_xxy_xxxxz, g_y_0_x_0_xxy_xxxyy, g_y_0_x_0_xxy_xxxyz, g_y_0_x_0_xxy_xxxzz, g_y_0_x_0_xxy_xxyyy, g_y_0_x_0_xxy_xxyyz, g_y_0_x_0_xxy_xxyzz, g_y_0_x_0_xxy_xxzzz, g_y_0_x_0_xxy_xyyyy, g_y_0_x_0_xxy_xyyyz, g_y_0_x_0_xxy_xyyzz, g_y_0_x_0_xxy_xyzzz, g_y_0_x_0_xxy_xzzzz, g_y_0_x_0_xxy_yyyyy, g_y_0_x_0_xxy_yyyyz, g_y_0_x_0_xxy_yyyzz, g_y_0_x_0_xxy_yyzzz, g_y_0_x_0_xxy_yzzzz, g_y_0_x_0_xxy_zzzzz, g_y_0_x_0_xy_xxxxx, g_y_0_x_0_xy_xxxxxx, g_y_0_x_0_xy_xxxxxy, g_y_0_x_0_xy_xxxxxz, g_y_0_x_0_xy_xxxxy, g_y_0_x_0_xy_xxxxyy, g_y_0_x_0_xy_xxxxyz, g_y_0_x_0_xy_xxxxz, g_y_0_x_0_xy_xxxxzz, g_y_0_x_0_xy_xxxyy, g_y_0_x_0_xy_xxxyyy, g_y_0_x_0_xy_xxxyyz, g_y_0_x_0_xy_xxxyz, g_y_0_x_0_xy_xxxyzz, g_y_0_x_0_xy_xxxzz, g_y_0_x_0_xy_xxxzzz, g_y_0_x_0_xy_xxyyy, g_y_0_x_0_xy_xxyyyy, g_y_0_x_0_xy_xxyyyz, g_y_0_x_0_xy_xxyyz, g_y_0_x_0_xy_xxyyzz, g_y_0_x_0_xy_xxyzz, g_y_0_x_0_xy_xxyzzz, g_y_0_x_0_xy_xxzzz, g_y_0_x_0_xy_xxzzzz, g_y_0_x_0_xy_xyyyy, g_y_0_x_0_xy_xyyyyy, g_y_0_x_0_xy_xyyyyz, g_y_0_x_0_xy_xyyyz, g_y_0_x_0_xy_xyyyzz, g_y_0_x_0_xy_xyyzz, g_y_0_x_0_xy_xyyzzz, g_y_0_x_0_xy_xyzzz, g_y_0_x_0_xy_xyzzzz, g_y_0_x_0_xy_xzzzz, g_y_0_x_0_xy_xzzzzz, g_y_0_x_0_xy_yyyyy, g_y_0_x_0_xy_yyyyz, g_y_0_x_0_xy_yyyzz, g_y_0_x_0_xy_yyzzz, g_y_0_x_0_xy_yzzzz, g_y_0_x_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxy_xxxxx[k] = -g_y_0_x_0_xy_xxxxx[k] * ab_x + g_y_0_x_0_xy_xxxxxx[k];

                g_y_0_x_0_xxy_xxxxy[k] = -g_y_0_x_0_xy_xxxxy[k] * ab_x + g_y_0_x_0_xy_xxxxxy[k];

                g_y_0_x_0_xxy_xxxxz[k] = -g_y_0_x_0_xy_xxxxz[k] * ab_x + g_y_0_x_0_xy_xxxxxz[k];

                g_y_0_x_0_xxy_xxxyy[k] = -g_y_0_x_0_xy_xxxyy[k] * ab_x + g_y_0_x_0_xy_xxxxyy[k];

                g_y_0_x_0_xxy_xxxyz[k] = -g_y_0_x_0_xy_xxxyz[k] * ab_x + g_y_0_x_0_xy_xxxxyz[k];

                g_y_0_x_0_xxy_xxxzz[k] = -g_y_0_x_0_xy_xxxzz[k] * ab_x + g_y_0_x_0_xy_xxxxzz[k];

                g_y_0_x_0_xxy_xxyyy[k] = -g_y_0_x_0_xy_xxyyy[k] * ab_x + g_y_0_x_0_xy_xxxyyy[k];

                g_y_0_x_0_xxy_xxyyz[k] = -g_y_0_x_0_xy_xxyyz[k] * ab_x + g_y_0_x_0_xy_xxxyyz[k];

                g_y_0_x_0_xxy_xxyzz[k] = -g_y_0_x_0_xy_xxyzz[k] * ab_x + g_y_0_x_0_xy_xxxyzz[k];

                g_y_0_x_0_xxy_xxzzz[k] = -g_y_0_x_0_xy_xxzzz[k] * ab_x + g_y_0_x_0_xy_xxxzzz[k];

                g_y_0_x_0_xxy_xyyyy[k] = -g_y_0_x_0_xy_xyyyy[k] * ab_x + g_y_0_x_0_xy_xxyyyy[k];

                g_y_0_x_0_xxy_xyyyz[k] = -g_y_0_x_0_xy_xyyyz[k] * ab_x + g_y_0_x_0_xy_xxyyyz[k];

                g_y_0_x_0_xxy_xyyzz[k] = -g_y_0_x_0_xy_xyyzz[k] * ab_x + g_y_0_x_0_xy_xxyyzz[k];

                g_y_0_x_0_xxy_xyzzz[k] = -g_y_0_x_0_xy_xyzzz[k] * ab_x + g_y_0_x_0_xy_xxyzzz[k];

                g_y_0_x_0_xxy_xzzzz[k] = -g_y_0_x_0_xy_xzzzz[k] * ab_x + g_y_0_x_0_xy_xxzzzz[k];

                g_y_0_x_0_xxy_yyyyy[k] = -g_y_0_x_0_xy_yyyyy[k] * ab_x + g_y_0_x_0_xy_xyyyyy[k];

                g_y_0_x_0_xxy_yyyyz[k] = -g_y_0_x_0_xy_yyyyz[k] * ab_x + g_y_0_x_0_xy_xyyyyz[k];

                g_y_0_x_0_xxy_yyyzz[k] = -g_y_0_x_0_xy_yyyzz[k] * ab_x + g_y_0_x_0_xy_xyyyzz[k];

                g_y_0_x_0_xxy_yyzzz[k] = -g_y_0_x_0_xy_yyzzz[k] * ab_x + g_y_0_x_0_xy_xyyzzz[k];

                g_y_0_x_0_xxy_yzzzz[k] = -g_y_0_x_0_xy_yzzzz[k] * ab_x + g_y_0_x_0_xy_xyzzzz[k];

                g_y_0_x_0_xxy_zzzzz[k] = -g_y_0_x_0_xy_zzzzz[k] * ab_x + g_y_0_x_0_xy_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 672 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 673 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 674 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 675 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 676 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 677 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 678 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 679 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 680 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 681 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 682 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 683 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 684 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 685 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 686 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 687 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 688 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 689 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 690 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 691 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxz_xxxxx, g_y_0_x_0_xxz_xxxxy, g_y_0_x_0_xxz_xxxxz, g_y_0_x_0_xxz_xxxyy, g_y_0_x_0_xxz_xxxyz, g_y_0_x_0_xxz_xxxzz, g_y_0_x_0_xxz_xxyyy, g_y_0_x_0_xxz_xxyyz, g_y_0_x_0_xxz_xxyzz, g_y_0_x_0_xxz_xxzzz, g_y_0_x_0_xxz_xyyyy, g_y_0_x_0_xxz_xyyyz, g_y_0_x_0_xxz_xyyzz, g_y_0_x_0_xxz_xyzzz, g_y_0_x_0_xxz_xzzzz, g_y_0_x_0_xxz_yyyyy, g_y_0_x_0_xxz_yyyyz, g_y_0_x_0_xxz_yyyzz, g_y_0_x_0_xxz_yyzzz, g_y_0_x_0_xxz_yzzzz, g_y_0_x_0_xxz_zzzzz, g_y_0_x_0_xz_xxxxx, g_y_0_x_0_xz_xxxxxx, g_y_0_x_0_xz_xxxxxy, g_y_0_x_0_xz_xxxxxz, g_y_0_x_0_xz_xxxxy, g_y_0_x_0_xz_xxxxyy, g_y_0_x_0_xz_xxxxyz, g_y_0_x_0_xz_xxxxz, g_y_0_x_0_xz_xxxxzz, g_y_0_x_0_xz_xxxyy, g_y_0_x_0_xz_xxxyyy, g_y_0_x_0_xz_xxxyyz, g_y_0_x_0_xz_xxxyz, g_y_0_x_0_xz_xxxyzz, g_y_0_x_0_xz_xxxzz, g_y_0_x_0_xz_xxxzzz, g_y_0_x_0_xz_xxyyy, g_y_0_x_0_xz_xxyyyy, g_y_0_x_0_xz_xxyyyz, g_y_0_x_0_xz_xxyyz, g_y_0_x_0_xz_xxyyzz, g_y_0_x_0_xz_xxyzz, g_y_0_x_0_xz_xxyzzz, g_y_0_x_0_xz_xxzzz, g_y_0_x_0_xz_xxzzzz, g_y_0_x_0_xz_xyyyy, g_y_0_x_0_xz_xyyyyy, g_y_0_x_0_xz_xyyyyz, g_y_0_x_0_xz_xyyyz, g_y_0_x_0_xz_xyyyzz, g_y_0_x_0_xz_xyyzz, g_y_0_x_0_xz_xyyzzz, g_y_0_x_0_xz_xyzzz, g_y_0_x_0_xz_xyzzzz, g_y_0_x_0_xz_xzzzz, g_y_0_x_0_xz_xzzzzz, g_y_0_x_0_xz_yyyyy, g_y_0_x_0_xz_yyyyz, g_y_0_x_0_xz_yyyzz, g_y_0_x_0_xz_yyzzz, g_y_0_x_0_xz_yzzzz, g_y_0_x_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxz_xxxxx[k] = -g_y_0_x_0_xz_xxxxx[k] * ab_x + g_y_0_x_0_xz_xxxxxx[k];

                g_y_0_x_0_xxz_xxxxy[k] = -g_y_0_x_0_xz_xxxxy[k] * ab_x + g_y_0_x_0_xz_xxxxxy[k];

                g_y_0_x_0_xxz_xxxxz[k] = -g_y_0_x_0_xz_xxxxz[k] * ab_x + g_y_0_x_0_xz_xxxxxz[k];

                g_y_0_x_0_xxz_xxxyy[k] = -g_y_0_x_0_xz_xxxyy[k] * ab_x + g_y_0_x_0_xz_xxxxyy[k];

                g_y_0_x_0_xxz_xxxyz[k] = -g_y_0_x_0_xz_xxxyz[k] * ab_x + g_y_0_x_0_xz_xxxxyz[k];

                g_y_0_x_0_xxz_xxxzz[k] = -g_y_0_x_0_xz_xxxzz[k] * ab_x + g_y_0_x_0_xz_xxxxzz[k];

                g_y_0_x_0_xxz_xxyyy[k] = -g_y_0_x_0_xz_xxyyy[k] * ab_x + g_y_0_x_0_xz_xxxyyy[k];

                g_y_0_x_0_xxz_xxyyz[k] = -g_y_0_x_0_xz_xxyyz[k] * ab_x + g_y_0_x_0_xz_xxxyyz[k];

                g_y_0_x_0_xxz_xxyzz[k] = -g_y_0_x_0_xz_xxyzz[k] * ab_x + g_y_0_x_0_xz_xxxyzz[k];

                g_y_0_x_0_xxz_xxzzz[k] = -g_y_0_x_0_xz_xxzzz[k] * ab_x + g_y_0_x_0_xz_xxxzzz[k];

                g_y_0_x_0_xxz_xyyyy[k] = -g_y_0_x_0_xz_xyyyy[k] * ab_x + g_y_0_x_0_xz_xxyyyy[k];

                g_y_0_x_0_xxz_xyyyz[k] = -g_y_0_x_0_xz_xyyyz[k] * ab_x + g_y_0_x_0_xz_xxyyyz[k];

                g_y_0_x_0_xxz_xyyzz[k] = -g_y_0_x_0_xz_xyyzz[k] * ab_x + g_y_0_x_0_xz_xxyyzz[k];

                g_y_0_x_0_xxz_xyzzz[k] = -g_y_0_x_0_xz_xyzzz[k] * ab_x + g_y_0_x_0_xz_xxyzzz[k];

                g_y_0_x_0_xxz_xzzzz[k] = -g_y_0_x_0_xz_xzzzz[k] * ab_x + g_y_0_x_0_xz_xxzzzz[k];

                g_y_0_x_0_xxz_yyyyy[k] = -g_y_0_x_0_xz_yyyyy[k] * ab_x + g_y_0_x_0_xz_xyyyyy[k];

                g_y_0_x_0_xxz_yyyyz[k] = -g_y_0_x_0_xz_yyyyz[k] * ab_x + g_y_0_x_0_xz_xyyyyz[k];

                g_y_0_x_0_xxz_yyyzz[k] = -g_y_0_x_0_xz_yyyzz[k] * ab_x + g_y_0_x_0_xz_xyyyzz[k];

                g_y_0_x_0_xxz_yyzzz[k] = -g_y_0_x_0_xz_yyzzz[k] * ab_x + g_y_0_x_0_xz_xyyzzz[k];

                g_y_0_x_0_xxz_yzzzz[k] = -g_y_0_x_0_xz_yzzzz[k] * ab_x + g_y_0_x_0_xz_xyzzzz[k];

                g_y_0_x_0_xxz_zzzzz[k] = -g_y_0_x_0_xz_zzzzz[k] * ab_x + g_y_0_x_0_xz_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 693 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 694 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 695 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 696 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 697 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 698 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 699 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 700 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 701 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 702 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 703 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 704 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 705 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 706 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 707 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 708 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 709 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 710 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 711 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 712 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyy_xxxxx, g_y_0_x_0_xyy_xxxxy, g_y_0_x_0_xyy_xxxxz, g_y_0_x_0_xyy_xxxyy, g_y_0_x_0_xyy_xxxyz, g_y_0_x_0_xyy_xxxzz, g_y_0_x_0_xyy_xxyyy, g_y_0_x_0_xyy_xxyyz, g_y_0_x_0_xyy_xxyzz, g_y_0_x_0_xyy_xxzzz, g_y_0_x_0_xyy_xyyyy, g_y_0_x_0_xyy_xyyyz, g_y_0_x_0_xyy_xyyzz, g_y_0_x_0_xyy_xyzzz, g_y_0_x_0_xyy_xzzzz, g_y_0_x_0_xyy_yyyyy, g_y_0_x_0_xyy_yyyyz, g_y_0_x_0_xyy_yyyzz, g_y_0_x_0_xyy_yyzzz, g_y_0_x_0_xyy_yzzzz, g_y_0_x_0_xyy_zzzzz, g_y_0_x_0_yy_xxxxx, g_y_0_x_0_yy_xxxxxx, g_y_0_x_0_yy_xxxxxy, g_y_0_x_0_yy_xxxxxz, g_y_0_x_0_yy_xxxxy, g_y_0_x_0_yy_xxxxyy, g_y_0_x_0_yy_xxxxyz, g_y_0_x_0_yy_xxxxz, g_y_0_x_0_yy_xxxxzz, g_y_0_x_0_yy_xxxyy, g_y_0_x_0_yy_xxxyyy, g_y_0_x_0_yy_xxxyyz, g_y_0_x_0_yy_xxxyz, g_y_0_x_0_yy_xxxyzz, g_y_0_x_0_yy_xxxzz, g_y_0_x_0_yy_xxxzzz, g_y_0_x_0_yy_xxyyy, g_y_0_x_0_yy_xxyyyy, g_y_0_x_0_yy_xxyyyz, g_y_0_x_0_yy_xxyyz, g_y_0_x_0_yy_xxyyzz, g_y_0_x_0_yy_xxyzz, g_y_0_x_0_yy_xxyzzz, g_y_0_x_0_yy_xxzzz, g_y_0_x_0_yy_xxzzzz, g_y_0_x_0_yy_xyyyy, g_y_0_x_0_yy_xyyyyy, g_y_0_x_0_yy_xyyyyz, g_y_0_x_0_yy_xyyyz, g_y_0_x_0_yy_xyyyzz, g_y_0_x_0_yy_xyyzz, g_y_0_x_0_yy_xyyzzz, g_y_0_x_0_yy_xyzzz, g_y_0_x_0_yy_xyzzzz, g_y_0_x_0_yy_xzzzz, g_y_0_x_0_yy_xzzzzz, g_y_0_x_0_yy_yyyyy, g_y_0_x_0_yy_yyyyz, g_y_0_x_0_yy_yyyzz, g_y_0_x_0_yy_yyzzz, g_y_0_x_0_yy_yzzzz, g_y_0_x_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyy_xxxxx[k] = -g_y_0_x_0_yy_xxxxx[k] * ab_x + g_y_0_x_0_yy_xxxxxx[k];

                g_y_0_x_0_xyy_xxxxy[k] = -g_y_0_x_0_yy_xxxxy[k] * ab_x + g_y_0_x_0_yy_xxxxxy[k];

                g_y_0_x_0_xyy_xxxxz[k] = -g_y_0_x_0_yy_xxxxz[k] * ab_x + g_y_0_x_0_yy_xxxxxz[k];

                g_y_0_x_0_xyy_xxxyy[k] = -g_y_0_x_0_yy_xxxyy[k] * ab_x + g_y_0_x_0_yy_xxxxyy[k];

                g_y_0_x_0_xyy_xxxyz[k] = -g_y_0_x_0_yy_xxxyz[k] * ab_x + g_y_0_x_0_yy_xxxxyz[k];

                g_y_0_x_0_xyy_xxxzz[k] = -g_y_0_x_0_yy_xxxzz[k] * ab_x + g_y_0_x_0_yy_xxxxzz[k];

                g_y_0_x_0_xyy_xxyyy[k] = -g_y_0_x_0_yy_xxyyy[k] * ab_x + g_y_0_x_0_yy_xxxyyy[k];

                g_y_0_x_0_xyy_xxyyz[k] = -g_y_0_x_0_yy_xxyyz[k] * ab_x + g_y_0_x_0_yy_xxxyyz[k];

                g_y_0_x_0_xyy_xxyzz[k] = -g_y_0_x_0_yy_xxyzz[k] * ab_x + g_y_0_x_0_yy_xxxyzz[k];

                g_y_0_x_0_xyy_xxzzz[k] = -g_y_0_x_0_yy_xxzzz[k] * ab_x + g_y_0_x_0_yy_xxxzzz[k];

                g_y_0_x_0_xyy_xyyyy[k] = -g_y_0_x_0_yy_xyyyy[k] * ab_x + g_y_0_x_0_yy_xxyyyy[k];

                g_y_0_x_0_xyy_xyyyz[k] = -g_y_0_x_0_yy_xyyyz[k] * ab_x + g_y_0_x_0_yy_xxyyyz[k];

                g_y_0_x_0_xyy_xyyzz[k] = -g_y_0_x_0_yy_xyyzz[k] * ab_x + g_y_0_x_0_yy_xxyyzz[k];

                g_y_0_x_0_xyy_xyzzz[k] = -g_y_0_x_0_yy_xyzzz[k] * ab_x + g_y_0_x_0_yy_xxyzzz[k];

                g_y_0_x_0_xyy_xzzzz[k] = -g_y_0_x_0_yy_xzzzz[k] * ab_x + g_y_0_x_0_yy_xxzzzz[k];

                g_y_0_x_0_xyy_yyyyy[k] = -g_y_0_x_0_yy_yyyyy[k] * ab_x + g_y_0_x_0_yy_xyyyyy[k];

                g_y_0_x_0_xyy_yyyyz[k] = -g_y_0_x_0_yy_yyyyz[k] * ab_x + g_y_0_x_0_yy_xyyyyz[k];

                g_y_0_x_0_xyy_yyyzz[k] = -g_y_0_x_0_yy_yyyzz[k] * ab_x + g_y_0_x_0_yy_xyyyzz[k];

                g_y_0_x_0_xyy_yyzzz[k] = -g_y_0_x_0_yy_yyzzz[k] * ab_x + g_y_0_x_0_yy_xyyzzz[k];

                g_y_0_x_0_xyy_yzzzz[k] = -g_y_0_x_0_yy_yzzzz[k] * ab_x + g_y_0_x_0_yy_xyzzzz[k];

                g_y_0_x_0_xyy_zzzzz[k] = -g_y_0_x_0_yy_zzzzz[k] * ab_x + g_y_0_x_0_yy_xzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 714 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 715 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 716 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 717 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 718 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 719 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 720 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 721 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 722 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 723 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 724 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 725 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 726 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 727 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 728 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 729 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 730 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 731 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 732 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 733 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyz_xxxxx, g_y_0_x_0_xyz_xxxxy, g_y_0_x_0_xyz_xxxxz, g_y_0_x_0_xyz_xxxyy, g_y_0_x_0_xyz_xxxyz, g_y_0_x_0_xyz_xxxzz, g_y_0_x_0_xyz_xxyyy, g_y_0_x_0_xyz_xxyyz, g_y_0_x_0_xyz_xxyzz, g_y_0_x_0_xyz_xxzzz, g_y_0_x_0_xyz_xyyyy, g_y_0_x_0_xyz_xyyyz, g_y_0_x_0_xyz_xyyzz, g_y_0_x_0_xyz_xyzzz, g_y_0_x_0_xyz_xzzzz, g_y_0_x_0_xyz_yyyyy, g_y_0_x_0_xyz_yyyyz, g_y_0_x_0_xyz_yyyzz, g_y_0_x_0_xyz_yyzzz, g_y_0_x_0_xyz_yzzzz, g_y_0_x_0_xyz_zzzzz, g_y_0_x_0_yz_xxxxx, g_y_0_x_0_yz_xxxxxx, g_y_0_x_0_yz_xxxxxy, g_y_0_x_0_yz_xxxxxz, g_y_0_x_0_yz_xxxxy, g_y_0_x_0_yz_xxxxyy, g_y_0_x_0_yz_xxxxyz, g_y_0_x_0_yz_xxxxz, g_y_0_x_0_yz_xxxxzz, g_y_0_x_0_yz_xxxyy, g_y_0_x_0_yz_xxxyyy, g_y_0_x_0_yz_xxxyyz, g_y_0_x_0_yz_xxxyz, g_y_0_x_0_yz_xxxyzz, g_y_0_x_0_yz_xxxzz, g_y_0_x_0_yz_xxxzzz, g_y_0_x_0_yz_xxyyy, g_y_0_x_0_yz_xxyyyy, g_y_0_x_0_yz_xxyyyz, g_y_0_x_0_yz_xxyyz, g_y_0_x_0_yz_xxyyzz, g_y_0_x_0_yz_xxyzz, g_y_0_x_0_yz_xxyzzz, g_y_0_x_0_yz_xxzzz, g_y_0_x_0_yz_xxzzzz, g_y_0_x_0_yz_xyyyy, g_y_0_x_0_yz_xyyyyy, g_y_0_x_0_yz_xyyyyz, g_y_0_x_0_yz_xyyyz, g_y_0_x_0_yz_xyyyzz, g_y_0_x_0_yz_xyyzz, g_y_0_x_0_yz_xyyzzz, g_y_0_x_0_yz_xyzzz, g_y_0_x_0_yz_xyzzzz, g_y_0_x_0_yz_xzzzz, g_y_0_x_0_yz_xzzzzz, g_y_0_x_0_yz_yyyyy, g_y_0_x_0_yz_yyyyz, g_y_0_x_0_yz_yyyzz, g_y_0_x_0_yz_yyzzz, g_y_0_x_0_yz_yzzzz, g_y_0_x_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyz_xxxxx[k] = -g_y_0_x_0_yz_xxxxx[k] * ab_x + g_y_0_x_0_yz_xxxxxx[k];

                g_y_0_x_0_xyz_xxxxy[k] = -g_y_0_x_0_yz_xxxxy[k] * ab_x + g_y_0_x_0_yz_xxxxxy[k];

                g_y_0_x_0_xyz_xxxxz[k] = -g_y_0_x_0_yz_xxxxz[k] * ab_x + g_y_0_x_0_yz_xxxxxz[k];

                g_y_0_x_0_xyz_xxxyy[k] = -g_y_0_x_0_yz_xxxyy[k] * ab_x + g_y_0_x_0_yz_xxxxyy[k];

                g_y_0_x_0_xyz_xxxyz[k] = -g_y_0_x_0_yz_xxxyz[k] * ab_x + g_y_0_x_0_yz_xxxxyz[k];

                g_y_0_x_0_xyz_xxxzz[k] = -g_y_0_x_0_yz_xxxzz[k] * ab_x + g_y_0_x_0_yz_xxxxzz[k];

                g_y_0_x_0_xyz_xxyyy[k] = -g_y_0_x_0_yz_xxyyy[k] * ab_x + g_y_0_x_0_yz_xxxyyy[k];

                g_y_0_x_0_xyz_xxyyz[k] = -g_y_0_x_0_yz_xxyyz[k] * ab_x + g_y_0_x_0_yz_xxxyyz[k];

                g_y_0_x_0_xyz_xxyzz[k] = -g_y_0_x_0_yz_xxyzz[k] * ab_x + g_y_0_x_0_yz_xxxyzz[k];

                g_y_0_x_0_xyz_xxzzz[k] = -g_y_0_x_0_yz_xxzzz[k] * ab_x + g_y_0_x_0_yz_xxxzzz[k];

                g_y_0_x_0_xyz_xyyyy[k] = -g_y_0_x_0_yz_xyyyy[k] * ab_x + g_y_0_x_0_yz_xxyyyy[k];

                g_y_0_x_0_xyz_xyyyz[k] = -g_y_0_x_0_yz_xyyyz[k] * ab_x + g_y_0_x_0_yz_xxyyyz[k];

                g_y_0_x_0_xyz_xyyzz[k] = -g_y_0_x_0_yz_xyyzz[k] * ab_x + g_y_0_x_0_yz_xxyyzz[k];

                g_y_0_x_0_xyz_xyzzz[k] = -g_y_0_x_0_yz_xyzzz[k] * ab_x + g_y_0_x_0_yz_xxyzzz[k];

                g_y_0_x_0_xyz_xzzzz[k] = -g_y_0_x_0_yz_xzzzz[k] * ab_x + g_y_0_x_0_yz_xxzzzz[k];

                g_y_0_x_0_xyz_yyyyy[k] = -g_y_0_x_0_yz_yyyyy[k] * ab_x + g_y_0_x_0_yz_xyyyyy[k];

                g_y_0_x_0_xyz_yyyyz[k] = -g_y_0_x_0_yz_yyyyz[k] * ab_x + g_y_0_x_0_yz_xyyyyz[k];

                g_y_0_x_0_xyz_yyyzz[k] = -g_y_0_x_0_yz_yyyzz[k] * ab_x + g_y_0_x_0_yz_xyyyzz[k];

                g_y_0_x_0_xyz_yyzzz[k] = -g_y_0_x_0_yz_yyzzz[k] * ab_x + g_y_0_x_0_yz_xyyzzz[k];

                g_y_0_x_0_xyz_yzzzz[k] = -g_y_0_x_0_yz_yzzzz[k] * ab_x + g_y_0_x_0_yz_xyzzzz[k];

                g_y_0_x_0_xyz_zzzzz[k] = -g_y_0_x_0_yz_zzzzz[k] * ab_x + g_y_0_x_0_yz_xzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 735 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 736 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 737 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 738 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 739 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 740 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 741 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 742 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 743 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 744 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 745 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 746 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 747 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 748 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 749 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 750 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 751 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 752 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 753 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 754 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzz_xxxxx, g_y_0_x_0_xzz_xxxxy, g_y_0_x_0_xzz_xxxxz, g_y_0_x_0_xzz_xxxyy, g_y_0_x_0_xzz_xxxyz, g_y_0_x_0_xzz_xxxzz, g_y_0_x_0_xzz_xxyyy, g_y_0_x_0_xzz_xxyyz, g_y_0_x_0_xzz_xxyzz, g_y_0_x_0_xzz_xxzzz, g_y_0_x_0_xzz_xyyyy, g_y_0_x_0_xzz_xyyyz, g_y_0_x_0_xzz_xyyzz, g_y_0_x_0_xzz_xyzzz, g_y_0_x_0_xzz_xzzzz, g_y_0_x_0_xzz_yyyyy, g_y_0_x_0_xzz_yyyyz, g_y_0_x_0_xzz_yyyzz, g_y_0_x_0_xzz_yyzzz, g_y_0_x_0_xzz_yzzzz, g_y_0_x_0_xzz_zzzzz, g_y_0_x_0_zz_xxxxx, g_y_0_x_0_zz_xxxxxx, g_y_0_x_0_zz_xxxxxy, g_y_0_x_0_zz_xxxxxz, g_y_0_x_0_zz_xxxxy, g_y_0_x_0_zz_xxxxyy, g_y_0_x_0_zz_xxxxyz, g_y_0_x_0_zz_xxxxz, g_y_0_x_0_zz_xxxxzz, g_y_0_x_0_zz_xxxyy, g_y_0_x_0_zz_xxxyyy, g_y_0_x_0_zz_xxxyyz, g_y_0_x_0_zz_xxxyz, g_y_0_x_0_zz_xxxyzz, g_y_0_x_0_zz_xxxzz, g_y_0_x_0_zz_xxxzzz, g_y_0_x_0_zz_xxyyy, g_y_0_x_0_zz_xxyyyy, g_y_0_x_0_zz_xxyyyz, g_y_0_x_0_zz_xxyyz, g_y_0_x_0_zz_xxyyzz, g_y_0_x_0_zz_xxyzz, g_y_0_x_0_zz_xxyzzz, g_y_0_x_0_zz_xxzzz, g_y_0_x_0_zz_xxzzzz, g_y_0_x_0_zz_xyyyy, g_y_0_x_0_zz_xyyyyy, g_y_0_x_0_zz_xyyyyz, g_y_0_x_0_zz_xyyyz, g_y_0_x_0_zz_xyyyzz, g_y_0_x_0_zz_xyyzz, g_y_0_x_0_zz_xyyzzz, g_y_0_x_0_zz_xyzzz, g_y_0_x_0_zz_xyzzzz, g_y_0_x_0_zz_xzzzz, g_y_0_x_0_zz_xzzzzz, g_y_0_x_0_zz_yyyyy, g_y_0_x_0_zz_yyyyz, g_y_0_x_0_zz_yyyzz, g_y_0_x_0_zz_yyzzz, g_y_0_x_0_zz_yzzzz, g_y_0_x_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzz_xxxxx[k] = -g_y_0_x_0_zz_xxxxx[k] * ab_x + g_y_0_x_0_zz_xxxxxx[k];

                g_y_0_x_0_xzz_xxxxy[k] = -g_y_0_x_0_zz_xxxxy[k] * ab_x + g_y_0_x_0_zz_xxxxxy[k];

                g_y_0_x_0_xzz_xxxxz[k] = -g_y_0_x_0_zz_xxxxz[k] * ab_x + g_y_0_x_0_zz_xxxxxz[k];

                g_y_0_x_0_xzz_xxxyy[k] = -g_y_0_x_0_zz_xxxyy[k] * ab_x + g_y_0_x_0_zz_xxxxyy[k];

                g_y_0_x_0_xzz_xxxyz[k] = -g_y_0_x_0_zz_xxxyz[k] * ab_x + g_y_0_x_0_zz_xxxxyz[k];

                g_y_0_x_0_xzz_xxxzz[k] = -g_y_0_x_0_zz_xxxzz[k] * ab_x + g_y_0_x_0_zz_xxxxzz[k];

                g_y_0_x_0_xzz_xxyyy[k] = -g_y_0_x_0_zz_xxyyy[k] * ab_x + g_y_0_x_0_zz_xxxyyy[k];

                g_y_0_x_0_xzz_xxyyz[k] = -g_y_0_x_0_zz_xxyyz[k] * ab_x + g_y_0_x_0_zz_xxxyyz[k];

                g_y_0_x_0_xzz_xxyzz[k] = -g_y_0_x_0_zz_xxyzz[k] * ab_x + g_y_0_x_0_zz_xxxyzz[k];

                g_y_0_x_0_xzz_xxzzz[k] = -g_y_0_x_0_zz_xxzzz[k] * ab_x + g_y_0_x_0_zz_xxxzzz[k];

                g_y_0_x_0_xzz_xyyyy[k] = -g_y_0_x_0_zz_xyyyy[k] * ab_x + g_y_0_x_0_zz_xxyyyy[k];

                g_y_0_x_0_xzz_xyyyz[k] = -g_y_0_x_0_zz_xyyyz[k] * ab_x + g_y_0_x_0_zz_xxyyyz[k];

                g_y_0_x_0_xzz_xyyzz[k] = -g_y_0_x_0_zz_xyyzz[k] * ab_x + g_y_0_x_0_zz_xxyyzz[k];

                g_y_0_x_0_xzz_xyzzz[k] = -g_y_0_x_0_zz_xyzzz[k] * ab_x + g_y_0_x_0_zz_xxyzzz[k];

                g_y_0_x_0_xzz_xzzzz[k] = -g_y_0_x_0_zz_xzzzz[k] * ab_x + g_y_0_x_0_zz_xxzzzz[k];

                g_y_0_x_0_xzz_yyyyy[k] = -g_y_0_x_0_zz_yyyyy[k] * ab_x + g_y_0_x_0_zz_xyyyyy[k];

                g_y_0_x_0_xzz_yyyyz[k] = -g_y_0_x_0_zz_yyyyz[k] * ab_x + g_y_0_x_0_zz_xyyyyz[k];

                g_y_0_x_0_xzz_yyyzz[k] = -g_y_0_x_0_zz_yyyzz[k] * ab_x + g_y_0_x_0_zz_xyyyzz[k];

                g_y_0_x_0_xzz_yyzzz[k] = -g_y_0_x_0_zz_yyzzz[k] * ab_x + g_y_0_x_0_zz_xyyzzz[k];

                g_y_0_x_0_xzz_yzzzz[k] = -g_y_0_x_0_zz_yzzzz[k] * ab_x + g_y_0_x_0_zz_xyzzzz[k];

                g_y_0_x_0_xzz_zzzzz[k] = -g_y_0_x_0_zz_zzzzz[k] * ab_x + g_y_0_x_0_zz_xzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 756 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 757 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 758 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 759 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 760 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 761 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 762 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 763 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 764 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 765 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 766 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 767 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 768 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 769 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 770 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 771 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 772 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 773 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 774 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 775 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 776 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yy_xxxxx, g_0_0_x_0_yy_xxxxy, g_0_0_x_0_yy_xxxxz, g_0_0_x_0_yy_xxxyy, g_0_0_x_0_yy_xxxyz, g_0_0_x_0_yy_xxxzz, g_0_0_x_0_yy_xxyyy, g_0_0_x_0_yy_xxyyz, g_0_0_x_0_yy_xxyzz, g_0_0_x_0_yy_xxzzz, g_0_0_x_0_yy_xyyyy, g_0_0_x_0_yy_xyyyz, g_0_0_x_0_yy_xyyzz, g_0_0_x_0_yy_xyzzz, g_0_0_x_0_yy_xzzzz, g_0_0_x_0_yy_yyyyy, g_0_0_x_0_yy_yyyyz, g_0_0_x_0_yy_yyyzz, g_0_0_x_0_yy_yyzzz, g_0_0_x_0_yy_yzzzz, g_0_0_x_0_yy_zzzzz, g_y_0_x_0_yy_xxxxx, g_y_0_x_0_yy_xxxxxy, g_y_0_x_0_yy_xxxxy, g_y_0_x_0_yy_xxxxyy, g_y_0_x_0_yy_xxxxyz, g_y_0_x_0_yy_xxxxz, g_y_0_x_0_yy_xxxyy, g_y_0_x_0_yy_xxxyyy, g_y_0_x_0_yy_xxxyyz, g_y_0_x_0_yy_xxxyz, g_y_0_x_0_yy_xxxyzz, g_y_0_x_0_yy_xxxzz, g_y_0_x_0_yy_xxyyy, g_y_0_x_0_yy_xxyyyy, g_y_0_x_0_yy_xxyyyz, g_y_0_x_0_yy_xxyyz, g_y_0_x_0_yy_xxyyzz, g_y_0_x_0_yy_xxyzz, g_y_0_x_0_yy_xxyzzz, g_y_0_x_0_yy_xxzzz, g_y_0_x_0_yy_xyyyy, g_y_0_x_0_yy_xyyyyy, g_y_0_x_0_yy_xyyyyz, g_y_0_x_0_yy_xyyyz, g_y_0_x_0_yy_xyyyzz, g_y_0_x_0_yy_xyyzz, g_y_0_x_0_yy_xyyzzz, g_y_0_x_0_yy_xyzzz, g_y_0_x_0_yy_xyzzzz, g_y_0_x_0_yy_xzzzz, g_y_0_x_0_yy_yyyyy, g_y_0_x_0_yy_yyyyyy, g_y_0_x_0_yy_yyyyyz, g_y_0_x_0_yy_yyyyz, g_y_0_x_0_yy_yyyyzz, g_y_0_x_0_yy_yyyzz, g_y_0_x_0_yy_yyyzzz, g_y_0_x_0_yy_yyzzz, g_y_0_x_0_yy_yyzzzz, g_y_0_x_0_yy_yzzzz, g_y_0_x_0_yy_yzzzzz, g_y_0_x_0_yy_zzzzz, g_y_0_x_0_yyy_xxxxx, g_y_0_x_0_yyy_xxxxy, g_y_0_x_0_yyy_xxxxz, g_y_0_x_0_yyy_xxxyy, g_y_0_x_0_yyy_xxxyz, g_y_0_x_0_yyy_xxxzz, g_y_0_x_0_yyy_xxyyy, g_y_0_x_0_yyy_xxyyz, g_y_0_x_0_yyy_xxyzz, g_y_0_x_0_yyy_xxzzz, g_y_0_x_0_yyy_xyyyy, g_y_0_x_0_yyy_xyyyz, g_y_0_x_0_yyy_xyyzz, g_y_0_x_0_yyy_xyzzz, g_y_0_x_0_yyy_xzzzz, g_y_0_x_0_yyy_yyyyy, g_y_0_x_0_yyy_yyyyz, g_y_0_x_0_yyy_yyyzz, g_y_0_x_0_yyy_yyzzz, g_y_0_x_0_yyy_yzzzz, g_y_0_x_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyy_xxxxx[k] = -g_0_0_x_0_yy_xxxxx[k] - g_y_0_x_0_yy_xxxxx[k] * ab_y + g_y_0_x_0_yy_xxxxxy[k];

                g_y_0_x_0_yyy_xxxxy[k] = -g_0_0_x_0_yy_xxxxy[k] - g_y_0_x_0_yy_xxxxy[k] * ab_y + g_y_0_x_0_yy_xxxxyy[k];

                g_y_0_x_0_yyy_xxxxz[k] = -g_0_0_x_0_yy_xxxxz[k] - g_y_0_x_0_yy_xxxxz[k] * ab_y + g_y_0_x_0_yy_xxxxyz[k];

                g_y_0_x_0_yyy_xxxyy[k] = -g_0_0_x_0_yy_xxxyy[k] - g_y_0_x_0_yy_xxxyy[k] * ab_y + g_y_0_x_0_yy_xxxyyy[k];

                g_y_0_x_0_yyy_xxxyz[k] = -g_0_0_x_0_yy_xxxyz[k] - g_y_0_x_0_yy_xxxyz[k] * ab_y + g_y_0_x_0_yy_xxxyyz[k];

                g_y_0_x_0_yyy_xxxzz[k] = -g_0_0_x_0_yy_xxxzz[k] - g_y_0_x_0_yy_xxxzz[k] * ab_y + g_y_0_x_0_yy_xxxyzz[k];

                g_y_0_x_0_yyy_xxyyy[k] = -g_0_0_x_0_yy_xxyyy[k] - g_y_0_x_0_yy_xxyyy[k] * ab_y + g_y_0_x_0_yy_xxyyyy[k];

                g_y_0_x_0_yyy_xxyyz[k] = -g_0_0_x_0_yy_xxyyz[k] - g_y_0_x_0_yy_xxyyz[k] * ab_y + g_y_0_x_0_yy_xxyyyz[k];

                g_y_0_x_0_yyy_xxyzz[k] = -g_0_0_x_0_yy_xxyzz[k] - g_y_0_x_0_yy_xxyzz[k] * ab_y + g_y_0_x_0_yy_xxyyzz[k];

                g_y_0_x_0_yyy_xxzzz[k] = -g_0_0_x_0_yy_xxzzz[k] - g_y_0_x_0_yy_xxzzz[k] * ab_y + g_y_0_x_0_yy_xxyzzz[k];

                g_y_0_x_0_yyy_xyyyy[k] = -g_0_0_x_0_yy_xyyyy[k] - g_y_0_x_0_yy_xyyyy[k] * ab_y + g_y_0_x_0_yy_xyyyyy[k];

                g_y_0_x_0_yyy_xyyyz[k] = -g_0_0_x_0_yy_xyyyz[k] - g_y_0_x_0_yy_xyyyz[k] * ab_y + g_y_0_x_0_yy_xyyyyz[k];

                g_y_0_x_0_yyy_xyyzz[k] = -g_0_0_x_0_yy_xyyzz[k] - g_y_0_x_0_yy_xyyzz[k] * ab_y + g_y_0_x_0_yy_xyyyzz[k];

                g_y_0_x_0_yyy_xyzzz[k] = -g_0_0_x_0_yy_xyzzz[k] - g_y_0_x_0_yy_xyzzz[k] * ab_y + g_y_0_x_0_yy_xyyzzz[k];

                g_y_0_x_0_yyy_xzzzz[k] = -g_0_0_x_0_yy_xzzzz[k] - g_y_0_x_0_yy_xzzzz[k] * ab_y + g_y_0_x_0_yy_xyzzzz[k];

                g_y_0_x_0_yyy_yyyyy[k] = -g_0_0_x_0_yy_yyyyy[k] - g_y_0_x_0_yy_yyyyy[k] * ab_y + g_y_0_x_0_yy_yyyyyy[k];

                g_y_0_x_0_yyy_yyyyz[k] = -g_0_0_x_0_yy_yyyyz[k] - g_y_0_x_0_yy_yyyyz[k] * ab_y + g_y_0_x_0_yy_yyyyyz[k];

                g_y_0_x_0_yyy_yyyzz[k] = -g_0_0_x_0_yy_yyyzz[k] - g_y_0_x_0_yy_yyyzz[k] * ab_y + g_y_0_x_0_yy_yyyyzz[k];

                g_y_0_x_0_yyy_yyzzz[k] = -g_0_0_x_0_yy_yyzzz[k] - g_y_0_x_0_yy_yyzzz[k] * ab_y + g_y_0_x_0_yy_yyyzzz[k];

                g_y_0_x_0_yyy_yzzzz[k] = -g_0_0_x_0_yy_yzzzz[k] - g_y_0_x_0_yy_yzzzz[k] * ab_y + g_y_0_x_0_yy_yyzzzz[k];

                g_y_0_x_0_yyy_zzzzz[k] = -g_0_0_x_0_yy_zzzzz[k] - g_y_0_x_0_yy_zzzzz[k] * ab_y + g_y_0_x_0_yy_yzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 777 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 778 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 779 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 780 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 781 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 782 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 783 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 784 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 785 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 786 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 787 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 788 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 789 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 790 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 791 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 792 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 793 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 794 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 795 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 796 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yy_xxxxx, g_y_0_x_0_yy_xxxxxz, g_y_0_x_0_yy_xxxxy, g_y_0_x_0_yy_xxxxyz, g_y_0_x_0_yy_xxxxz, g_y_0_x_0_yy_xxxxzz, g_y_0_x_0_yy_xxxyy, g_y_0_x_0_yy_xxxyyz, g_y_0_x_0_yy_xxxyz, g_y_0_x_0_yy_xxxyzz, g_y_0_x_0_yy_xxxzz, g_y_0_x_0_yy_xxxzzz, g_y_0_x_0_yy_xxyyy, g_y_0_x_0_yy_xxyyyz, g_y_0_x_0_yy_xxyyz, g_y_0_x_0_yy_xxyyzz, g_y_0_x_0_yy_xxyzz, g_y_0_x_0_yy_xxyzzz, g_y_0_x_0_yy_xxzzz, g_y_0_x_0_yy_xxzzzz, g_y_0_x_0_yy_xyyyy, g_y_0_x_0_yy_xyyyyz, g_y_0_x_0_yy_xyyyz, g_y_0_x_0_yy_xyyyzz, g_y_0_x_0_yy_xyyzz, g_y_0_x_0_yy_xyyzzz, g_y_0_x_0_yy_xyzzz, g_y_0_x_0_yy_xyzzzz, g_y_0_x_0_yy_xzzzz, g_y_0_x_0_yy_xzzzzz, g_y_0_x_0_yy_yyyyy, g_y_0_x_0_yy_yyyyyz, g_y_0_x_0_yy_yyyyz, g_y_0_x_0_yy_yyyyzz, g_y_0_x_0_yy_yyyzz, g_y_0_x_0_yy_yyyzzz, g_y_0_x_0_yy_yyzzz, g_y_0_x_0_yy_yyzzzz, g_y_0_x_0_yy_yzzzz, g_y_0_x_0_yy_yzzzzz, g_y_0_x_0_yy_zzzzz, g_y_0_x_0_yy_zzzzzz, g_y_0_x_0_yyz_xxxxx, g_y_0_x_0_yyz_xxxxy, g_y_0_x_0_yyz_xxxxz, g_y_0_x_0_yyz_xxxyy, g_y_0_x_0_yyz_xxxyz, g_y_0_x_0_yyz_xxxzz, g_y_0_x_0_yyz_xxyyy, g_y_0_x_0_yyz_xxyyz, g_y_0_x_0_yyz_xxyzz, g_y_0_x_0_yyz_xxzzz, g_y_0_x_0_yyz_xyyyy, g_y_0_x_0_yyz_xyyyz, g_y_0_x_0_yyz_xyyzz, g_y_0_x_0_yyz_xyzzz, g_y_0_x_0_yyz_xzzzz, g_y_0_x_0_yyz_yyyyy, g_y_0_x_0_yyz_yyyyz, g_y_0_x_0_yyz_yyyzz, g_y_0_x_0_yyz_yyzzz, g_y_0_x_0_yyz_yzzzz, g_y_0_x_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyz_xxxxx[k] = -g_y_0_x_0_yy_xxxxx[k] * ab_z + g_y_0_x_0_yy_xxxxxz[k];

                g_y_0_x_0_yyz_xxxxy[k] = -g_y_0_x_0_yy_xxxxy[k] * ab_z + g_y_0_x_0_yy_xxxxyz[k];

                g_y_0_x_0_yyz_xxxxz[k] = -g_y_0_x_0_yy_xxxxz[k] * ab_z + g_y_0_x_0_yy_xxxxzz[k];

                g_y_0_x_0_yyz_xxxyy[k] = -g_y_0_x_0_yy_xxxyy[k] * ab_z + g_y_0_x_0_yy_xxxyyz[k];

                g_y_0_x_0_yyz_xxxyz[k] = -g_y_0_x_0_yy_xxxyz[k] * ab_z + g_y_0_x_0_yy_xxxyzz[k];

                g_y_0_x_0_yyz_xxxzz[k] = -g_y_0_x_0_yy_xxxzz[k] * ab_z + g_y_0_x_0_yy_xxxzzz[k];

                g_y_0_x_0_yyz_xxyyy[k] = -g_y_0_x_0_yy_xxyyy[k] * ab_z + g_y_0_x_0_yy_xxyyyz[k];

                g_y_0_x_0_yyz_xxyyz[k] = -g_y_0_x_0_yy_xxyyz[k] * ab_z + g_y_0_x_0_yy_xxyyzz[k];

                g_y_0_x_0_yyz_xxyzz[k] = -g_y_0_x_0_yy_xxyzz[k] * ab_z + g_y_0_x_0_yy_xxyzzz[k];

                g_y_0_x_0_yyz_xxzzz[k] = -g_y_0_x_0_yy_xxzzz[k] * ab_z + g_y_0_x_0_yy_xxzzzz[k];

                g_y_0_x_0_yyz_xyyyy[k] = -g_y_0_x_0_yy_xyyyy[k] * ab_z + g_y_0_x_0_yy_xyyyyz[k];

                g_y_0_x_0_yyz_xyyyz[k] = -g_y_0_x_0_yy_xyyyz[k] * ab_z + g_y_0_x_0_yy_xyyyzz[k];

                g_y_0_x_0_yyz_xyyzz[k] = -g_y_0_x_0_yy_xyyzz[k] * ab_z + g_y_0_x_0_yy_xyyzzz[k];

                g_y_0_x_0_yyz_xyzzz[k] = -g_y_0_x_0_yy_xyzzz[k] * ab_z + g_y_0_x_0_yy_xyzzzz[k];

                g_y_0_x_0_yyz_xzzzz[k] = -g_y_0_x_0_yy_xzzzz[k] * ab_z + g_y_0_x_0_yy_xzzzzz[k];

                g_y_0_x_0_yyz_yyyyy[k] = -g_y_0_x_0_yy_yyyyy[k] * ab_z + g_y_0_x_0_yy_yyyyyz[k];

                g_y_0_x_0_yyz_yyyyz[k] = -g_y_0_x_0_yy_yyyyz[k] * ab_z + g_y_0_x_0_yy_yyyyzz[k];

                g_y_0_x_0_yyz_yyyzz[k] = -g_y_0_x_0_yy_yyyzz[k] * ab_z + g_y_0_x_0_yy_yyyzzz[k];

                g_y_0_x_0_yyz_yyzzz[k] = -g_y_0_x_0_yy_yyzzz[k] * ab_z + g_y_0_x_0_yy_yyzzzz[k];

                g_y_0_x_0_yyz_yzzzz[k] = -g_y_0_x_0_yy_yzzzz[k] * ab_z + g_y_0_x_0_yy_yzzzzz[k];

                g_y_0_x_0_yyz_zzzzz[k] = -g_y_0_x_0_yy_zzzzz[k] * ab_z + g_y_0_x_0_yy_zzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 798 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 799 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 800 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 801 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 802 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 803 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 804 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 805 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 806 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 807 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 808 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 809 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 810 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 811 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 812 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 813 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 814 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 815 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 816 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 817 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 818 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yz_xxxxx, g_y_0_x_0_yz_xxxxxz, g_y_0_x_0_yz_xxxxy, g_y_0_x_0_yz_xxxxyz, g_y_0_x_0_yz_xxxxz, g_y_0_x_0_yz_xxxxzz, g_y_0_x_0_yz_xxxyy, g_y_0_x_0_yz_xxxyyz, g_y_0_x_0_yz_xxxyz, g_y_0_x_0_yz_xxxyzz, g_y_0_x_0_yz_xxxzz, g_y_0_x_0_yz_xxxzzz, g_y_0_x_0_yz_xxyyy, g_y_0_x_0_yz_xxyyyz, g_y_0_x_0_yz_xxyyz, g_y_0_x_0_yz_xxyyzz, g_y_0_x_0_yz_xxyzz, g_y_0_x_0_yz_xxyzzz, g_y_0_x_0_yz_xxzzz, g_y_0_x_0_yz_xxzzzz, g_y_0_x_0_yz_xyyyy, g_y_0_x_0_yz_xyyyyz, g_y_0_x_0_yz_xyyyz, g_y_0_x_0_yz_xyyyzz, g_y_0_x_0_yz_xyyzz, g_y_0_x_0_yz_xyyzzz, g_y_0_x_0_yz_xyzzz, g_y_0_x_0_yz_xyzzzz, g_y_0_x_0_yz_xzzzz, g_y_0_x_0_yz_xzzzzz, g_y_0_x_0_yz_yyyyy, g_y_0_x_0_yz_yyyyyz, g_y_0_x_0_yz_yyyyz, g_y_0_x_0_yz_yyyyzz, g_y_0_x_0_yz_yyyzz, g_y_0_x_0_yz_yyyzzz, g_y_0_x_0_yz_yyzzz, g_y_0_x_0_yz_yyzzzz, g_y_0_x_0_yz_yzzzz, g_y_0_x_0_yz_yzzzzz, g_y_0_x_0_yz_zzzzz, g_y_0_x_0_yz_zzzzzz, g_y_0_x_0_yzz_xxxxx, g_y_0_x_0_yzz_xxxxy, g_y_0_x_0_yzz_xxxxz, g_y_0_x_0_yzz_xxxyy, g_y_0_x_0_yzz_xxxyz, g_y_0_x_0_yzz_xxxzz, g_y_0_x_0_yzz_xxyyy, g_y_0_x_0_yzz_xxyyz, g_y_0_x_0_yzz_xxyzz, g_y_0_x_0_yzz_xxzzz, g_y_0_x_0_yzz_xyyyy, g_y_0_x_0_yzz_xyyyz, g_y_0_x_0_yzz_xyyzz, g_y_0_x_0_yzz_xyzzz, g_y_0_x_0_yzz_xzzzz, g_y_0_x_0_yzz_yyyyy, g_y_0_x_0_yzz_yyyyz, g_y_0_x_0_yzz_yyyzz, g_y_0_x_0_yzz_yyzzz, g_y_0_x_0_yzz_yzzzz, g_y_0_x_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzz_xxxxx[k] = -g_y_0_x_0_yz_xxxxx[k] * ab_z + g_y_0_x_0_yz_xxxxxz[k];

                g_y_0_x_0_yzz_xxxxy[k] = -g_y_0_x_0_yz_xxxxy[k] * ab_z + g_y_0_x_0_yz_xxxxyz[k];

                g_y_0_x_0_yzz_xxxxz[k] = -g_y_0_x_0_yz_xxxxz[k] * ab_z + g_y_0_x_0_yz_xxxxzz[k];

                g_y_0_x_0_yzz_xxxyy[k] = -g_y_0_x_0_yz_xxxyy[k] * ab_z + g_y_0_x_0_yz_xxxyyz[k];

                g_y_0_x_0_yzz_xxxyz[k] = -g_y_0_x_0_yz_xxxyz[k] * ab_z + g_y_0_x_0_yz_xxxyzz[k];

                g_y_0_x_0_yzz_xxxzz[k] = -g_y_0_x_0_yz_xxxzz[k] * ab_z + g_y_0_x_0_yz_xxxzzz[k];

                g_y_0_x_0_yzz_xxyyy[k] = -g_y_0_x_0_yz_xxyyy[k] * ab_z + g_y_0_x_0_yz_xxyyyz[k];

                g_y_0_x_0_yzz_xxyyz[k] = -g_y_0_x_0_yz_xxyyz[k] * ab_z + g_y_0_x_0_yz_xxyyzz[k];

                g_y_0_x_0_yzz_xxyzz[k] = -g_y_0_x_0_yz_xxyzz[k] * ab_z + g_y_0_x_0_yz_xxyzzz[k];

                g_y_0_x_0_yzz_xxzzz[k] = -g_y_0_x_0_yz_xxzzz[k] * ab_z + g_y_0_x_0_yz_xxzzzz[k];

                g_y_0_x_0_yzz_xyyyy[k] = -g_y_0_x_0_yz_xyyyy[k] * ab_z + g_y_0_x_0_yz_xyyyyz[k];

                g_y_0_x_0_yzz_xyyyz[k] = -g_y_0_x_0_yz_xyyyz[k] * ab_z + g_y_0_x_0_yz_xyyyzz[k];

                g_y_0_x_0_yzz_xyyzz[k] = -g_y_0_x_0_yz_xyyzz[k] * ab_z + g_y_0_x_0_yz_xyyzzz[k];

                g_y_0_x_0_yzz_xyzzz[k] = -g_y_0_x_0_yz_xyzzz[k] * ab_z + g_y_0_x_0_yz_xyzzzz[k];

                g_y_0_x_0_yzz_xzzzz[k] = -g_y_0_x_0_yz_xzzzz[k] * ab_z + g_y_0_x_0_yz_xzzzzz[k];

                g_y_0_x_0_yzz_yyyyy[k] = -g_y_0_x_0_yz_yyyyy[k] * ab_z + g_y_0_x_0_yz_yyyyyz[k];

                g_y_0_x_0_yzz_yyyyz[k] = -g_y_0_x_0_yz_yyyyz[k] * ab_z + g_y_0_x_0_yz_yyyyzz[k];

                g_y_0_x_0_yzz_yyyzz[k] = -g_y_0_x_0_yz_yyyzz[k] * ab_z + g_y_0_x_0_yz_yyyzzz[k];

                g_y_0_x_0_yzz_yyzzz[k] = -g_y_0_x_0_yz_yyzzz[k] * ab_z + g_y_0_x_0_yz_yyzzzz[k];

                g_y_0_x_0_yzz_yzzzz[k] = -g_y_0_x_0_yz_yzzzz[k] * ab_z + g_y_0_x_0_yz_yzzzzz[k];

                g_y_0_x_0_yzz_zzzzz[k] = -g_y_0_x_0_yz_zzzzz[k] * ab_z + g_y_0_x_0_yz_zzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 819 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 820 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 821 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 822 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 823 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 824 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 825 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 826 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 827 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 828 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 829 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 830 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 831 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 832 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 833 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 834 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 835 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 836 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 837 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 838 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zz_xxxxx, g_y_0_x_0_zz_xxxxxz, g_y_0_x_0_zz_xxxxy, g_y_0_x_0_zz_xxxxyz, g_y_0_x_0_zz_xxxxz, g_y_0_x_0_zz_xxxxzz, g_y_0_x_0_zz_xxxyy, g_y_0_x_0_zz_xxxyyz, g_y_0_x_0_zz_xxxyz, g_y_0_x_0_zz_xxxyzz, g_y_0_x_0_zz_xxxzz, g_y_0_x_0_zz_xxxzzz, g_y_0_x_0_zz_xxyyy, g_y_0_x_0_zz_xxyyyz, g_y_0_x_0_zz_xxyyz, g_y_0_x_0_zz_xxyyzz, g_y_0_x_0_zz_xxyzz, g_y_0_x_0_zz_xxyzzz, g_y_0_x_0_zz_xxzzz, g_y_0_x_0_zz_xxzzzz, g_y_0_x_0_zz_xyyyy, g_y_0_x_0_zz_xyyyyz, g_y_0_x_0_zz_xyyyz, g_y_0_x_0_zz_xyyyzz, g_y_0_x_0_zz_xyyzz, g_y_0_x_0_zz_xyyzzz, g_y_0_x_0_zz_xyzzz, g_y_0_x_0_zz_xyzzzz, g_y_0_x_0_zz_xzzzz, g_y_0_x_0_zz_xzzzzz, g_y_0_x_0_zz_yyyyy, g_y_0_x_0_zz_yyyyyz, g_y_0_x_0_zz_yyyyz, g_y_0_x_0_zz_yyyyzz, g_y_0_x_0_zz_yyyzz, g_y_0_x_0_zz_yyyzzz, g_y_0_x_0_zz_yyzzz, g_y_0_x_0_zz_yyzzzz, g_y_0_x_0_zz_yzzzz, g_y_0_x_0_zz_yzzzzz, g_y_0_x_0_zz_zzzzz, g_y_0_x_0_zz_zzzzzz, g_y_0_x_0_zzz_xxxxx, g_y_0_x_0_zzz_xxxxy, g_y_0_x_0_zzz_xxxxz, g_y_0_x_0_zzz_xxxyy, g_y_0_x_0_zzz_xxxyz, g_y_0_x_0_zzz_xxxzz, g_y_0_x_0_zzz_xxyyy, g_y_0_x_0_zzz_xxyyz, g_y_0_x_0_zzz_xxyzz, g_y_0_x_0_zzz_xxzzz, g_y_0_x_0_zzz_xyyyy, g_y_0_x_0_zzz_xyyyz, g_y_0_x_0_zzz_xyyzz, g_y_0_x_0_zzz_xyzzz, g_y_0_x_0_zzz_xzzzz, g_y_0_x_0_zzz_yyyyy, g_y_0_x_0_zzz_yyyyz, g_y_0_x_0_zzz_yyyzz, g_y_0_x_0_zzz_yyzzz, g_y_0_x_0_zzz_yzzzz, g_y_0_x_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzz_xxxxx[k] = -g_y_0_x_0_zz_xxxxx[k] * ab_z + g_y_0_x_0_zz_xxxxxz[k];

                g_y_0_x_0_zzz_xxxxy[k] = -g_y_0_x_0_zz_xxxxy[k] * ab_z + g_y_0_x_0_zz_xxxxyz[k];

                g_y_0_x_0_zzz_xxxxz[k] = -g_y_0_x_0_zz_xxxxz[k] * ab_z + g_y_0_x_0_zz_xxxxzz[k];

                g_y_0_x_0_zzz_xxxyy[k] = -g_y_0_x_0_zz_xxxyy[k] * ab_z + g_y_0_x_0_zz_xxxyyz[k];

                g_y_0_x_0_zzz_xxxyz[k] = -g_y_0_x_0_zz_xxxyz[k] * ab_z + g_y_0_x_0_zz_xxxyzz[k];

                g_y_0_x_0_zzz_xxxzz[k] = -g_y_0_x_0_zz_xxxzz[k] * ab_z + g_y_0_x_0_zz_xxxzzz[k];

                g_y_0_x_0_zzz_xxyyy[k] = -g_y_0_x_0_zz_xxyyy[k] * ab_z + g_y_0_x_0_zz_xxyyyz[k];

                g_y_0_x_0_zzz_xxyyz[k] = -g_y_0_x_0_zz_xxyyz[k] * ab_z + g_y_0_x_0_zz_xxyyzz[k];

                g_y_0_x_0_zzz_xxyzz[k] = -g_y_0_x_0_zz_xxyzz[k] * ab_z + g_y_0_x_0_zz_xxyzzz[k];

                g_y_0_x_0_zzz_xxzzz[k] = -g_y_0_x_0_zz_xxzzz[k] * ab_z + g_y_0_x_0_zz_xxzzzz[k];

                g_y_0_x_0_zzz_xyyyy[k] = -g_y_0_x_0_zz_xyyyy[k] * ab_z + g_y_0_x_0_zz_xyyyyz[k];

                g_y_0_x_0_zzz_xyyyz[k] = -g_y_0_x_0_zz_xyyyz[k] * ab_z + g_y_0_x_0_zz_xyyyzz[k];

                g_y_0_x_0_zzz_xyyzz[k] = -g_y_0_x_0_zz_xyyzz[k] * ab_z + g_y_0_x_0_zz_xyyzzz[k];

                g_y_0_x_0_zzz_xyzzz[k] = -g_y_0_x_0_zz_xyzzz[k] * ab_z + g_y_0_x_0_zz_xyzzzz[k];

                g_y_0_x_0_zzz_xzzzz[k] = -g_y_0_x_0_zz_xzzzz[k] * ab_z + g_y_0_x_0_zz_xzzzzz[k];

                g_y_0_x_0_zzz_yyyyy[k] = -g_y_0_x_0_zz_yyyyy[k] * ab_z + g_y_0_x_0_zz_yyyyyz[k];

                g_y_0_x_0_zzz_yyyyz[k] = -g_y_0_x_0_zz_yyyyz[k] * ab_z + g_y_0_x_0_zz_yyyyzz[k];

                g_y_0_x_0_zzz_yyyzz[k] = -g_y_0_x_0_zz_yyyzz[k] * ab_z + g_y_0_x_0_zz_yyyzzz[k];

                g_y_0_x_0_zzz_yyzzz[k] = -g_y_0_x_0_zz_yyzzz[k] * ab_z + g_y_0_x_0_zz_yyzzzz[k];

                g_y_0_x_0_zzz_yzzzz[k] = -g_y_0_x_0_zz_yzzzz[k] * ab_z + g_y_0_x_0_zz_yzzzzz[k];

                g_y_0_x_0_zzz_zzzzz[k] = -g_y_0_x_0_zz_zzzzz[k] * ab_z + g_y_0_x_0_zz_zzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 840 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 841 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 842 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 843 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 844 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 845 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 846 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 847 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 848 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 849 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 850 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 851 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 852 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 853 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 854 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 855 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 856 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 857 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 858 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 859 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 860 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xx_xxxxx, g_y_0_y_0_xx_xxxxxx, g_y_0_y_0_xx_xxxxxy, g_y_0_y_0_xx_xxxxxz, g_y_0_y_0_xx_xxxxy, g_y_0_y_0_xx_xxxxyy, g_y_0_y_0_xx_xxxxyz, g_y_0_y_0_xx_xxxxz, g_y_0_y_0_xx_xxxxzz, g_y_0_y_0_xx_xxxyy, g_y_0_y_0_xx_xxxyyy, g_y_0_y_0_xx_xxxyyz, g_y_0_y_0_xx_xxxyz, g_y_0_y_0_xx_xxxyzz, g_y_0_y_0_xx_xxxzz, g_y_0_y_0_xx_xxxzzz, g_y_0_y_0_xx_xxyyy, g_y_0_y_0_xx_xxyyyy, g_y_0_y_0_xx_xxyyyz, g_y_0_y_0_xx_xxyyz, g_y_0_y_0_xx_xxyyzz, g_y_0_y_0_xx_xxyzz, g_y_0_y_0_xx_xxyzzz, g_y_0_y_0_xx_xxzzz, g_y_0_y_0_xx_xxzzzz, g_y_0_y_0_xx_xyyyy, g_y_0_y_0_xx_xyyyyy, g_y_0_y_0_xx_xyyyyz, g_y_0_y_0_xx_xyyyz, g_y_0_y_0_xx_xyyyzz, g_y_0_y_0_xx_xyyzz, g_y_0_y_0_xx_xyyzzz, g_y_0_y_0_xx_xyzzz, g_y_0_y_0_xx_xyzzzz, g_y_0_y_0_xx_xzzzz, g_y_0_y_0_xx_xzzzzz, g_y_0_y_0_xx_yyyyy, g_y_0_y_0_xx_yyyyz, g_y_0_y_0_xx_yyyzz, g_y_0_y_0_xx_yyzzz, g_y_0_y_0_xx_yzzzz, g_y_0_y_0_xx_zzzzz, g_y_0_y_0_xxx_xxxxx, g_y_0_y_0_xxx_xxxxy, g_y_0_y_0_xxx_xxxxz, g_y_0_y_0_xxx_xxxyy, g_y_0_y_0_xxx_xxxyz, g_y_0_y_0_xxx_xxxzz, g_y_0_y_0_xxx_xxyyy, g_y_0_y_0_xxx_xxyyz, g_y_0_y_0_xxx_xxyzz, g_y_0_y_0_xxx_xxzzz, g_y_0_y_0_xxx_xyyyy, g_y_0_y_0_xxx_xyyyz, g_y_0_y_0_xxx_xyyzz, g_y_0_y_0_xxx_xyzzz, g_y_0_y_0_xxx_xzzzz, g_y_0_y_0_xxx_yyyyy, g_y_0_y_0_xxx_yyyyz, g_y_0_y_0_xxx_yyyzz, g_y_0_y_0_xxx_yyzzz, g_y_0_y_0_xxx_yzzzz, g_y_0_y_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxx_xxxxx[k] = -g_y_0_y_0_xx_xxxxx[k] * ab_x + g_y_0_y_0_xx_xxxxxx[k];

                g_y_0_y_0_xxx_xxxxy[k] = -g_y_0_y_0_xx_xxxxy[k] * ab_x + g_y_0_y_0_xx_xxxxxy[k];

                g_y_0_y_0_xxx_xxxxz[k] = -g_y_0_y_0_xx_xxxxz[k] * ab_x + g_y_0_y_0_xx_xxxxxz[k];

                g_y_0_y_0_xxx_xxxyy[k] = -g_y_0_y_0_xx_xxxyy[k] * ab_x + g_y_0_y_0_xx_xxxxyy[k];

                g_y_0_y_0_xxx_xxxyz[k] = -g_y_0_y_0_xx_xxxyz[k] * ab_x + g_y_0_y_0_xx_xxxxyz[k];

                g_y_0_y_0_xxx_xxxzz[k] = -g_y_0_y_0_xx_xxxzz[k] * ab_x + g_y_0_y_0_xx_xxxxzz[k];

                g_y_0_y_0_xxx_xxyyy[k] = -g_y_0_y_0_xx_xxyyy[k] * ab_x + g_y_0_y_0_xx_xxxyyy[k];

                g_y_0_y_0_xxx_xxyyz[k] = -g_y_0_y_0_xx_xxyyz[k] * ab_x + g_y_0_y_0_xx_xxxyyz[k];

                g_y_0_y_0_xxx_xxyzz[k] = -g_y_0_y_0_xx_xxyzz[k] * ab_x + g_y_0_y_0_xx_xxxyzz[k];

                g_y_0_y_0_xxx_xxzzz[k] = -g_y_0_y_0_xx_xxzzz[k] * ab_x + g_y_0_y_0_xx_xxxzzz[k];

                g_y_0_y_0_xxx_xyyyy[k] = -g_y_0_y_0_xx_xyyyy[k] * ab_x + g_y_0_y_0_xx_xxyyyy[k];

                g_y_0_y_0_xxx_xyyyz[k] = -g_y_0_y_0_xx_xyyyz[k] * ab_x + g_y_0_y_0_xx_xxyyyz[k];

                g_y_0_y_0_xxx_xyyzz[k] = -g_y_0_y_0_xx_xyyzz[k] * ab_x + g_y_0_y_0_xx_xxyyzz[k];

                g_y_0_y_0_xxx_xyzzz[k] = -g_y_0_y_0_xx_xyzzz[k] * ab_x + g_y_0_y_0_xx_xxyzzz[k];

                g_y_0_y_0_xxx_xzzzz[k] = -g_y_0_y_0_xx_xzzzz[k] * ab_x + g_y_0_y_0_xx_xxzzzz[k];

                g_y_0_y_0_xxx_yyyyy[k] = -g_y_0_y_0_xx_yyyyy[k] * ab_x + g_y_0_y_0_xx_xyyyyy[k];

                g_y_0_y_0_xxx_yyyyz[k] = -g_y_0_y_0_xx_yyyyz[k] * ab_x + g_y_0_y_0_xx_xyyyyz[k];

                g_y_0_y_0_xxx_yyyzz[k] = -g_y_0_y_0_xx_yyyzz[k] * ab_x + g_y_0_y_0_xx_xyyyzz[k];

                g_y_0_y_0_xxx_yyzzz[k] = -g_y_0_y_0_xx_yyzzz[k] * ab_x + g_y_0_y_0_xx_xyyzzz[k];

                g_y_0_y_0_xxx_yzzzz[k] = -g_y_0_y_0_xx_yzzzz[k] * ab_x + g_y_0_y_0_xx_xyzzzz[k];

                g_y_0_y_0_xxx_zzzzz[k] = -g_y_0_y_0_xx_zzzzz[k] * ab_x + g_y_0_y_0_xx_xzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 861 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 862 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 863 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 864 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 865 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 866 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 867 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 868 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 869 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 870 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 871 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 872 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 873 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 874 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 875 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 876 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 877 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 878 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 879 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 880 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxy_xxxxx, g_y_0_y_0_xxy_xxxxy, g_y_0_y_0_xxy_xxxxz, g_y_0_y_0_xxy_xxxyy, g_y_0_y_0_xxy_xxxyz, g_y_0_y_0_xxy_xxxzz, g_y_0_y_0_xxy_xxyyy, g_y_0_y_0_xxy_xxyyz, g_y_0_y_0_xxy_xxyzz, g_y_0_y_0_xxy_xxzzz, g_y_0_y_0_xxy_xyyyy, g_y_0_y_0_xxy_xyyyz, g_y_0_y_0_xxy_xyyzz, g_y_0_y_0_xxy_xyzzz, g_y_0_y_0_xxy_xzzzz, g_y_0_y_0_xxy_yyyyy, g_y_0_y_0_xxy_yyyyz, g_y_0_y_0_xxy_yyyzz, g_y_0_y_0_xxy_yyzzz, g_y_0_y_0_xxy_yzzzz, g_y_0_y_0_xxy_zzzzz, g_y_0_y_0_xy_xxxxx, g_y_0_y_0_xy_xxxxxx, g_y_0_y_0_xy_xxxxxy, g_y_0_y_0_xy_xxxxxz, g_y_0_y_0_xy_xxxxy, g_y_0_y_0_xy_xxxxyy, g_y_0_y_0_xy_xxxxyz, g_y_0_y_0_xy_xxxxz, g_y_0_y_0_xy_xxxxzz, g_y_0_y_0_xy_xxxyy, g_y_0_y_0_xy_xxxyyy, g_y_0_y_0_xy_xxxyyz, g_y_0_y_0_xy_xxxyz, g_y_0_y_0_xy_xxxyzz, g_y_0_y_0_xy_xxxzz, g_y_0_y_0_xy_xxxzzz, g_y_0_y_0_xy_xxyyy, g_y_0_y_0_xy_xxyyyy, g_y_0_y_0_xy_xxyyyz, g_y_0_y_0_xy_xxyyz, g_y_0_y_0_xy_xxyyzz, g_y_0_y_0_xy_xxyzz, g_y_0_y_0_xy_xxyzzz, g_y_0_y_0_xy_xxzzz, g_y_0_y_0_xy_xxzzzz, g_y_0_y_0_xy_xyyyy, g_y_0_y_0_xy_xyyyyy, g_y_0_y_0_xy_xyyyyz, g_y_0_y_0_xy_xyyyz, g_y_0_y_0_xy_xyyyzz, g_y_0_y_0_xy_xyyzz, g_y_0_y_0_xy_xyyzzz, g_y_0_y_0_xy_xyzzz, g_y_0_y_0_xy_xyzzzz, g_y_0_y_0_xy_xzzzz, g_y_0_y_0_xy_xzzzzz, g_y_0_y_0_xy_yyyyy, g_y_0_y_0_xy_yyyyz, g_y_0_y_0_xy_yyyzz, g_y_0_y_0_xy_yyzzz, g_y_0_y_0_xy_yzzzz, g_y_0_y_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxy_xxxxx[k] = -g_y_0_y_0_xy_xxxxx[k] * ab_x + g_y_0_y_0_xy_xxxxxx[k];

                g_y_0_y_0_xxy_xxxxy[k] = -g_y_0_y_0_xy_xxxxy[k] * ab_x + g_y_0_y_0_xy_xxxxxy[k];

                g_y_0_y_0_xxy_xxxxz[k] = -g_y_0_y_0_xy_xxxxz[k] * ab_x + g_y_0_y_0_xy_xxxxxz[k];

                g_y_0_y_0_xxy_xxxyy[k] = -g_y_0_y_0_xy_xxxyy[k] * ab_x + g_y_0_y_0_xy_xxxxyy[k];

                g_y_0_y_0_xxy_xxxyz[k] = -g_y_0_y_0_xy_xxxyz[k] * ab_x + g_y_0_y_0_xy_xxxxyz[k];

                g_y_0_y_0_xxy_xxxzz[k] = -g_y_0_y_0_xy_xxxzz[k] * ab_x + g_y_0_y_0_xy_xxxxzz[k];

                g_y_0_y_0_xxy_xxyyy[k] = -g_y_0_y_0_xy_xxyyy[k] * ab_x + g_y_0_y_0_xy_xxxyyy[k];

                g_y_0_y_0_xxy_xxyyz[k] = -g_y_0_y_0_xy_xxyyz[k] * ab_x + g_y_0_y_0_xy_xxxyyz[k];

                g_y_0_y_0_xxy_xxyzz[k] = -g_y_0_y_0_xy_xxyzz[k] * ab_x + g_y_0_y_0_xy_xxxyzz[k];

                g_y_0_y_0_xxy_xxzzz[k] = -g_y_0_y_0_xy_xxzzz[k] * ab_x + g_y_0_y_0_xy_xxxzzz[k];

                g_y_0_y_0_xxy_xyyyy[k] = -g_y_0_y_0_xy_xyyyy[k] * ab_x + g_y_0_y_0_xy_xxyyyy[k];

                g_y_0_y_0_xxy_xyyyz[k] = -g_y_0_y_0_xy_xyyyz[k] * ab_x + g_y_0_y_0_xy_xxyyyz[k];

                g_y_0_y_0_xxy_xyyzz[k] = -g_y_0_y_0_xy_xyyzz[k] * ab_x + g_y_0_y_0_xy_xxyyzz[k];

                g_y_0_y_0_xxy_xyzzz[k] = -g_y_0_y_0_xy_xyzzz[k] * ab_x + g_y_0_y_0_xy_xxyzzz[k];

                g_y_0_y_0_xxy_xzzzz[k] = -g_y_0_y_0_xy_xzzzz[k] * ab_x + g_y_0_y_0_xy_xxzzzz[k];

                g_y_0_y_0_xxy_yyyyy[k] = -g_y_0_y_0_xy_yyyyy[k] * ab_x + g_y_0_y_0_xy_xyyyyy[k];

                g_y_0_y_0_xxy_yyyyz[k] = -g_y_0_y_0_xy_yyyyz[k] * ab_x + g_y_0_y_0_xy_xyyyyz[k];

                g_y_0_y_0_xxy_yyyzz[k] = -g_y_0_y_0_xy_yyyzz[k] * ab_x + g_y_0_y_0_xy_xyyyzz[k];

                g_y_0_y_0_xxy_yyzzz[k] = -g_y_0_y_0_xy_yyzzz[k] * ab_x + g_y_0_y_0_xy_xyyzzz[k];

                g_y_0_y_0_xxy_yzzzz[k] = -g_y_0_y_0_xy_yzzzz[k] * ab_x + g_y_0_y_0_xy_xyzzzz[k];

                g_y_0_y_0_xxy_zzzzz[k] = -g_y_0_y_0_xy_zzzzz[k] * ab_x + g_y_0_y_0_xy_xzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 882 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 883 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 884 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 885 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 886 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 887 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 888 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 889 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 890 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 891 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 892 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 893 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 894 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 895 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 896 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 897 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 898 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 899 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 900 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 901 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 902 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxz_xxxxx, g_y_0_y_0_xxz_xxxxy, g_y_0_y_0_xxz_xxxxz, g_y_0_y_0_xxz_xxxyy, g_y_0_y_0_xxz_xxxyz, g_y_0_y_0_xxz_xxxzz, g_y_0_y_0_xxz_xxyyy, g_y_0_y_0_xxz_xxyyz, g_y_0_y_0_xxz_xxyzz, g_y_0_y_0_xxz_xxzzz, g_y_0_y_0_xxz_xyyyy, g_y_0_y_0_xxz_xyyyz, g_y_0_y_0_xxz_xyyzz, g_y_0_y_0_xxz_xyzzz, g_y_0_y_0_xxz_xzzzz, g_y_0_y_0_xxz_yyyyy, g_y_0_y_0_xxz_yyyyz, g_y_0_y_0_xxz_yyyzz, g_y_0_y_0_xxz_yyzzz, g_y_0_y_0_xxz_yzzzz, g_y_0_y_0_xxz_zzzzz, g_y_0_y_0_xz_xxxxx, g_y_0_y_0_xz_xxxxxx, g_y_0_y_0_xz_xxxxxy, g_y_0_y_0_xz_xxxxxz, g_y_0_y_0_xz_xxxxy, g_y_0_y_0_xz_xxxxyy, g_y_0_y_0_xz_xxxxyz, g_y_0_y_0_xz_xxxxz, g_y_0_y_0_xz_xxxxzz, g_y_0_y_0_xz_xxxyy, g_y_0_y_0_xz_xxxyyy, g_y_0_y_0_xz_xxxyyz, g_y_0_y_0_xz_xxxyz, g_y_0_y_0_xz_xxxyzz, g_y_0_y_0_xz_xxxzz, g_y_0_y_0_xz_xxxzzz, g_y_0_y_0_xz_xxyyy, g_y_0_y_0_xz_xxyyyy, g_y_0_y_0_xz_xxyyyz, g_y_0_y_0_xz_xxyyz, g_y_0_y_0_xz_xxyyzz, g_y_0_y_0_xz_xxyzz, g_y_0_y_0_xz_xxyzzz, g_y_0_y_0_xz_xxzzz, g_y_0_y_0_xz_xxzzzz, g_y_0_y_0_xz_xyyyy, g_y_0_y_0_xz_xyyyyy, g_y_0_y_0_xz_xyyyyz, g_y_0_y_0_xz_xyyyz, g_y_0_y_0_xz_xyyyzz, g_y_0_y_0_xz_xyyzz, g_y_0_y_0_xz_xyyzzz, g_y_0_y_0_xz_xyzzz, g_y_0_y_0_xz_xyzzzz, g_y_0_y_0_xz_xzzzz, g_y_0_y_0_xz_xzzzzz, g_y_0_y_0_xz_yyyyy, g_y_0_y_0_xz_yyyyz, g_y_0_y_0_xz_yyyzz, g_y_0_y_0_xz_yyzzz, g_y_0_y_0_xz_yzzzz, g_y_0_y_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxz_xxxxx[k] = -g_y_0_y_0_xz_xxxxx[k] * ab_x + g_y_0_y_0_xz_xxxxxx[k];

                g_y_0_y_0_xxz_xxxxy[k] = -g_y_0_y_0_xz_xxxxy[k] * ab_x + g_y_0_y_0_xz_xxxxxy[k];

                g_y_0_y_0_xxz_xxxxz[k] = -g_y_0_y_0_xz_xxxxz[k] * ab_x + g_y_0_y_0_xz_xxxxxz[k];

                g_y_0_y_0_xxz_xxxyy[k] = -g_y_0_y_0_xz_xxxyy[k] * ab_x + g_y_0_y_0_xz_xxxxyy[k];

                g_y_0_y_0_xxz_xxxyz[k] = -g_y_0_y_0_xz_xxxyz[k] * ab_x + g_y_0_y_0_xz_xxxxyz[k];

                g_y_0_y_0_xxz_xxxzz[k] = -g_y_0_y_0_xz_xxxzz[k] * ab_x + g_y_0_y_0_xz_xxxxzz[k];

                g_y_0_y_0_xxz_xxyyy[k] = -g_y_0_y_0_xz_xxyyy[k] * ab_x + g_y_0_y_0_xz_xxxyyy[k];

                g_y_0_y_0_xxz_xxyyz[k] = -g_y_0_y_0_xz_xxyyz[k] * ab_x + g_y_0_y_0_xz_xxxyyz[k];

                g_y_0_y_0_xxz_xxyzz[k] = -g_y_0_y_0_xz_xxyzz[k] * ab_x + g_y_0_y_0_xz_xxxyzz[k];

                g_y_0_y_0_xxz_xxzzz[k] = -g_y_0_y_0_xz_xxzzz[k] * ab_x + g_y_0_y_0_xz_xxxzzz[k];

                g_y_0_y_0_xxz_xyyyy[k] = -g_y_0_y_0_xz_xyyyy[k] * ab_x + g_y_0_y_0_xz_xxyyyy[k];

                g_y_0_y_0_xxz_xyyyz[k] = -g_y_0_y_0_xz_xyyyz[k] * ab_x + g_y_0_y_0_xz_xxyyyz[k];

                g_y_0_y_0_xxz_xyyzz[k] = -g_y_0_y_0_xz_xyyzz[k] * ab_x + g_y_0_y_0_xz_xxyyzz[k];

                g_y_0_y_0_xxz_xyzzz[k] = -g_y_0_y_0_xz_xyzzz[k] * ab_x + g_y_0_y_0_xz_xxyzzz[k];

                g_y_0_y_0_xxz_xzzzz[k] = -g_y_0_y_0_xz_xzzzz[k] * ab_x + g_y_0_y_0_xz_xxzzzz[k];

                g_y_0_y_0_xxz_yyyyy[k] = -g_y_0_y_0_xz_yyyyy[k] * ab_x + g_y_0_y_0_xz_xyyyyy[k];

                g_y_0_y_0_xxz_yyyyz[k] = -g_y_0_y_0_xz_yyyyz[k] * ab_x + g_y_0_y_0_xz_xyyyyz[k];

                g_y_0_y_0_xxz_yyyzz[k] = -g_y_0_y_0_xz_yyyzz[k] * ab_x + g_y_0_y_0_xz_xyyyzz[k];

                g_y_0_y_0_xxz_yyzzz[k] = -g_y_0_y_0_xz_yyzzz[k] * ab_x + g_y_0_y_0_xz_xyyzzz[k];

                g_y_0_y_0_xxz_yzzzz[k] = -g_y_0_y_0_xz_yzzzz[k] * ab_x + g_y_0_y_0_xz_xyzzzz[k];

                g_y_0_y_0_xxz_zzzzz[k] = -g_y_0_y_0_xz_zzzzz[k] * ab_x + g_y_0_y_0_xz_xzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 903 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 904 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 905 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 906 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 907 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 908 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 909 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 910 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 911 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 912 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 913 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 914 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 915 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 916 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 917 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 918 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 919 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 920 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 921 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 922 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyy_xxxxx, g_y_0_y_0_xyy_xxxxy, g_y_0_y_0_xyy_xxxxz, g_y_0_y_0_xyy_xxxyy, g_y_0_y_0_xyy_xxxyz, g_y_0_y_0_xyy_xxxzz, g_y_0_y_0_xyy_xxyyy, g_y_0_y_0_xyy_xxyyz, g_y_0_y_0_xyy_xxyzz, g_y_0_y_0_xyy_xxzzz, g_y_0_y_0_xyy_xyyyy, g_y_0_y_0_xyy_xyyyz, g_y_0_y_0_xyy_xyyzz, g_y_0_y_0_xyy_xyzzz, g_y_0_y_0_xyy_xzzzz, g_y_0_y_0_xyy_yyyyy, g_y_0_y_0_xyy_yyyyz, g_y_0_y_0_xyy_yyyzz, g_y_0_y_0_xyy_yyzzz, g_y_0_y_0_xyy_yzzzz, g_y_0_y_0_xyy_zzzzz, g_y_0_y_0_yy_xxxxx, g_y_0_y_0_yy_xxxxxx, g_y_0_y_0_yy_xxxxxy, g_y_0_y_0_yy_xxxxxz, g_y_0_y_0_yy_xxxxy, g_y_0_y_0_yy_xxxxyy, g_y_0_y_0_yy_xxxxyz, g_y_0_y_0_yy_xxxxz, g_y_0_y_0_yy_xxxxzz, g_y_0_y_0_yy_xxxyy, g_y_0_y_0_yy_xxxyyy, g_y_0_y_0_yy_xxxyyz, g_y_0_y_0_yy_xxxyz, g_y_0_y_0_yy_xxxyzz, g_y_0_y_0_yy_xxxzz, g_y_0_y_0_yy_xxxzzz, g_y_0_y_0_yy_xxyyy, g_y_0_y_0_yy_xxyyyy, g_y_0_y_0_yy_xxyyyz, g_y_0_y_0_yy_xxyyz, g_y_0_y_0_yy_xxyyzz, g_y_0_y_0_yy_xxyzz, g_y_0_y_0_yy_xxyzzz, g_y_0_y_0_yy_xxzzz, g_y_0_y_0_yy_xxzzzz, g_y_0_y_0_yy_xyyyy, g_y_0_y_0_yy_xyyyyy, g_y_0_y_0_yy_xyyyyz, g_y_0_y_0_yy_xyyyz, g_y_0_y_0_yy_xyyyzz, g_y_0_y_0_yy_xyyzz, g_y_0_y_0_yy_xyyzzz, g_y_0_y_0_yy_xyzzz, g_y_0_y_0_yy_xyzzzz, g_y_0_y_0_yy_xzzzz, g_y_0_y_0_yy_xzzzzz, g_y_0_y_0_yy_yyyyy, g_y_0_y_0_yy_yyyyz, g_y_0_y_0_yy_yyyzz, g_y_0_y_0_yy_yyzzz, g_y_0_y_0_yy_yzzzz, g_y_0_y_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyy_xxxxx[k] = -g_y_0_y_0_yy_xxxxx[k] * ab_x + g_y_0_y_0_yy_xxxxxx[k];

                g_y_0_y_0_xyy_xxxxy[k] = -g_y_0_y_0_yy_xxxxy[k] * ab_x + g_y_0_y_0_yy_xxxxxy[k];

                g_y_0_y_0_xyy_xxxxz[k] = -g_y_0_y_0_yy_xxxxz[k] * ab_x + g_y_0_y_0_yy_xxxxxz[k];

                g_y_0_y_0_xyy_xxxyy[k] = -g_y_0_y_0_yy_xxxyy[k] * ab_x + g_y_0_y_0_yy_xxxxyy[k];

                g_y_0_y_0_xyy_xxxyz[k] = -g_y_0_y_0_yy_xxxyz[k] * ab_x + g_y_0_y_0_yy_xxxxyz[k];

                g_y_0_y_0_xyy_xxxzz[k] = -g_y_0_y_0_yy_xxxzz[k] * ab_x + g_y_0_y_0_yy_xxxxzz[k];

                g_y_0_y_0_xyy_xxyyy[k] = -g_y_0_y_0_yy_xxyyy[k] * ab_x + g_y_0_y_0_yy_xxxyyy[k];

                g_y_0_y_0_xyy_xxyyz[k] = -g_y_0_y_0_yy_xxyyz[k] * ab_x + g_y_0_y_0_yy_xxxyyz[k];

                g_y_0_y_0_xyy_xxyzz[k] = -g_y_0_y_0_yy_xxyzz[k] * ab_x + g_y_0_y_0_yy_xxxyzz[k];

                g_y_0_y_0_xyy_xxzzz[k] = -g_y_0_y_0_yy_xxzzz[k] * ab_x + g_y_0_y_0_yy_xxxzzz[k];

                g_y_0_y_0_xyy_xyyyy[k] = -g_y_0_y_0_yy_xyyyy[k] * ab_x + g_y_0_y_0_yy_xxyyyy[k];

                g_y_0_y_0_xyy_xyyyz[k] = -g_y_0_y_0_yy_xyyyz[k] * ab_x + g_y_0_y_0_yy_xxyyyz[k];

                g_y_0_y_0_xyy_xyyzz[k] = -g_y_0_y_0_yy_xyyzz[k] * ab_x + g_y_0_y_0_yy_xxyyzz[k];

                g_y_0_y_0_xyy_xyzzz[k] = -g_y_0_y_0_yy_xyzzz[k] * ab_x + g_y_0_y_0_yy_xxyzzz[k];

                g_y_0_y_0_xyy_xzzzz[k] = -g_y_0_y_0_yy_xzzzz[k] * ab_x + g_y_0_y_0_yy_xxzzzz[k];

                g_y_0_y_0_xyy_yyyyy[k] = -g_y_0_y_0_yy_yyyyy[k] * ab_x + g_y_0_y_0_yy_xyyyyy[k];

                g_y_0_y_0_xyy_yyyyz[k] = -g_y_0_y_0_yy_yyyyz[k] * ab_x + g_y_0_y_0_yy_xyyyyz[k];

                g_y_0_y_0_xyy_yyyzz[k] = -g_y_0_y_0_yy_yyyzz[k] * ab_x + g_y_0_y_0_yy_xyyyzz[k];

                g_y_0_y_0_xyy_yyzzz[k] = -g_y_0_y_0_yy_yyzzz[k] * ab_x + g_y_0_y_0_yy_xyyzzz[k];

                g_y_0_y_0_xyy_yzzzz[k] = -g_y_0_y_0_yy_yzzzz[k] * ab_x + g_y_0_y_0_yy_xyzzzz[k];

                g_y_0_y_0_xyy_zzzzz[k] = -g_y_0_y_0_yy_zzzzz[k] * ab_x + g_y_0_y_0_yy_xzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 924 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 925 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 926 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 927 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 928 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 929 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 930 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 931 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 932 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 933 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 934 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 935 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 936 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 937 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 938 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 939 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 940 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 941 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 942 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 943 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyz_xxxxx, g_y_0_y_0_xyz_xxxxy, g_y_0_y_0_xyz_xxxxz, g_y_0_y_0_xyz_xxxyy, g_y_0_y_0_xyz_xxxyz, g_y_0_y_0_xyz_xxxzz, g_y_0_y_0_xyz_xxyyy, g_y_0_y_0_xyz_xxyyz, g_y_0_y_0_xyz_xxyzz, g_y_0_y_0_xyz_xxzzz, g_y_0_y_0_xyz_xyyyy, g_y_0_y_0_xyz_xyyyz, g_y_0_y_0_xyz_xyyzz, g_y_0_y_0_xyz_xyzzz, g_y_0_y_0_xyz_xzzzz, g_y_0_y_0_xyz_yyyyy, g_y_0_y_0_xyz_yyyyz, g_y_0_y_0_xyz_yyyzz, g_y_0_y_0_xyz_yyzzz, g_y_0_y_0_xyz_yzzzz, g_y_0_y_0_xyz_zzzzz, g_y_0_y_0_yz_xxxxx, g_y_0_y_0_yz_xxxxxx, g_y_0_y_0_yz_xxxxxy, g_y_0_y_0_yz_xxxxxz, g_y_0_y_0_yz_xxxxy, g_y_0_y_0_yz_xxxxyy, g_y_0_y_0_yz_xxxxyz, g_y_0_y_0_yz_xxxxz, g_y_0_y_0_yz_xxxxzz, g_y_0_y_0_yz_xxxyy, g_y_0_y_0_yz_xxxyyy, g_y_0_y_0_yz_xxxyyz, g_y_0_y_0_yz_xxxyz, g_y_0_y_0_yz_xxxyzz, g_y_0_y_0_yz_xxxzz, g_y_0_y_0_yz_xxxzzz, g_y_0_y_0_yz_xxyyy, g_y_0_y_0_yz_xxyyyy, g_y_0_y_0_yz_xxyyyz, g_y_0_y_0_yz_xxyyz, g_y_0_y_0_yz_xxyyzz, g_y_0_y_0_yz_xxyzz, g_y_0_y_0_yz_xxyzzz, g_y_0_y_0_yz_xxzzz, g_y_0_y_0_yz_xxzzzz, g_y_0_y_0_yz_xyyyy, g_y_0_y_0_yz_xyyyyy, g_y_0_y_0_yz_xyyyyz, g_y_0_y_0_yz_xyyyz, g_y_0_y_0_yz_xyyyzz, g_y_0_y_0_yz_xyyzz, g_y_0_y_0_yz_xyyzzz, g_y_0_y_0_yz_xyzzz, g_y_0_y_0_yz_xyzzzz, g_y_0_y_0_yz_xzzzz, g_y_0_y_0_yz_xzzzzz, g_y_0_y_0_yz_yyyyy, g_y_0_y_0_yz_yyyyz, g_y_0_y_0_yz_yyyzz, g_y_0_y_0_yz_yyzzz, g_y_0_y_0_yz_yzzzz, g_y_0_y_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyz_xxxxx[k] = -g_y_0_y_0_yz_xxxxx[k] * ab_x + g_y_0_y_0_yz_xxxxxx[k];

                g_y_0_y_0_xyz_xxxxy[k] = -g_y_0_y_0_yz_xxxxy[k] * ab_x + g_y_0_y_0_yz_xxxxxy[k];

                g_y_0_y_0_xyz_xxxxz[k] = -g_y_0_y_0_yz_xxxxz[k] * ab_x + g_y_0_y_0_yz_xxxxxz[k];

                g_y_0_y_0_xyz_xxxyy[k] = -g_y_0_y_0_yz_xxxyy[k] * ab_x + g_y_0_y_0_yz_xxxxyy[k];

                g_y_0_y_0_xyz_xxxyz[k] = -g_y_0_y_0_yz_xxxyz[k] * ab_x + g_y_0_y_0_yz_xxxxyz[k];

                g_y_0_y_0_xyz_xxxzz[k] = -g_y_0_y_0_yz_xxxzz[k] * ab_x + g_y_0_y_0_yz_xxxxzz[k];

                g_y_0_y_0_xyz_xxyyy[k] = -g_y_0_y_0_yz_xxyyy[k] * ab_x + g_y_0_y_0_yz_xxxyyy[k];

                g_y_0_y_0_xyz_xxyyz[k] = -g_y_0_y_0_yz_xxyyz[k] * ab_x + g_y_0_y_0_yz_xxxyyz[k];

                g_y_0_y_0_xyz_xxyzz[k] = -g_y_0_y_0_yz_xxyzz[k] * ab_x + g_y_0_y_0_yz_xxxyzz[k];

                g_y_0_y_0_xyz_xxzzz[k] = -g_y_0_y_0_yz_xxzzz[k] * ab_x + g_y_0_y_0_yz_xxxzzz[k];

                g_y_0_y_0_xyz_xyyyy[k] = -g_y_0_y_0_yz_xyyyy[k] * ab_x + g_y_0_y_0_yz_xxyyyy[k];

                g_y_0_y_0_xyz_xyyyz[k] = -g_y_0_y_0_yz_xyyyz[k] * ab_x + g_y_0_y_0_yz_xxyyyz[k];

                g_y_0_y_0_xyz_xyyzz[k] = -g_y_0_y_0_yz_xyyzz[k] * ab_x + g_y_0_y_0_yz_xxyyzz[k];

                g_y_0_y_0_xyz_xyzzz[k] = -g_y_0_y_0_yz_xyzzz[k] * ab_x + g_y_0_y_0_yz_xxyzzz[k];

                g_y_0_y_0_xyz_xzzzz[k] = -g_y_0_y_0_yz_xzzzz[k] * ab_x + g_y_0_y_0_yz_xxzzzz[k];

                g_y_0_y_0_xyz_yyyyy[k] = -g_y_0_y_0_yz_yyyyy[k] * ab_x + g_y_0_y_0_yz_xyyyyy[k];

                g_y_0_y_0_xyz_yyyyz[k] = -g_y_0_y_0_yz_yyyyz[k] * ab_x + g_y_0_y_0_yz_xyyyyz[k];

                g_y_0_y_0_xyz_yyyzz[k] = -g_y_0_y_0_yz_yyyzz[k] * ab_x + g_y_0_y_0_yz_xyyyzz[k];

                g_y_0_y_0_xyz_yyzzz[k] = -g_y_0_y_0_yz_yyzzz[k] * ab_x + g_y_0_y_0_yz_xyyzzz[k];

                g_y_0_y_0_xyz_yzzzz[k] = -g_y_0_y_0_yz_yzzzz[k] * ab_x + g_y_0_y_0_yz_xyzzzz[k];

                g_y_0_y_0_xyz_zzzzz[k] = -g_y_0_y_0_yz_zzzzz[k] * ab_x + g_y_0_y_0_yz_xzzzzz[k];
            }

            /// Set up 945-966 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 945 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 946 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 947 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 948 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 949 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 950 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 951 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 952 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 953 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 954 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 955 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 956 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 957 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 958 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 959 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 960 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 961 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 962 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 963 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 964 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzz_xxxxx, g_y_0_y_0_xzz_xxxxy, g_y_0_y_0_xzz_xxxxz, g_y_0_y_0_xzz_xxxyy, g_y_0_y_0_xzz_xxxyz, g_y_0_y_0_xzz_xxxzz, g_y_0_y_0_xzz_xxyyy, g_y_0_y_0_xzz_xxyyz, g_y_0_y_0_xzz_xxyzz, g_y_0_y_0_xzz_xxzzz, g_y_0_y_0_xzz_xyyyy, g_y_0_y_0_xzz_xyyyz, g_y_0_y_0_xzz_xyyzz, g_y_0_y_0_xzz_xyzzz, g_y_0_y_0_xzz_xzzzz, g_y_0_y_0_xzz_yyyyy, g_y_0_y_0_xzz_yyyyz, g_y_0_y_0_xzz_yyyzz, g_y_0_y_0_xzz_yyzzz, g_y_0_y_0_xzz_yzzzz, g_y_0_y_0_xzz_zzzzz, g_y_0_y_0_zz_xxxxx, g_y_0_y_0_zz_xxxxxx, g_y_0_y_0_zz_xxxxxy, g_y_0_y_0_zz_xxxxxz, g_y_0_y_0_zz_xxxxy, g_y_0_y_0_zz_xxxxyy, g_y_0_y_0_zz_xxxxyz, g_y_0_y_0_zz_xxxxz, g_y_0_y_0_zz_xxxxzz, g_y_0_y_0_zz_xxxyy, g_y_0_y_0_zz_xxxyyy, g_y_0_y_0_zz_xxxyyz, g_y_0_y_0_zz_xxxyz, g_y_0_y_0_zz_xxxyzz, g_y_0_y_0_zz_xxxzz, g_y_0_y_0_zz_xxxzzz, g_y_0_y_0_zz_xxyyy, g_y_0_y_0_zz_xxyyyy, g_y_0_y_0_zz_xxyyyz, g_y_0_y_0_zz_xxyyz, g_y_0_y_0_zz_xxyyzz, g_y_0_y_0_zz_xxyzz, g_y_0_y_0_zz_xxyzzz, g_y_0_y_0_zz_xxzzz, g_y_0_y_0_zz_xxzzzz, g_y_0_y_0_zz_xyyyy, g_y_0_y_0_zz_xyyyyy, g_y_0_y_0_zz_xyyyyz, g_y_0_y_0_zz_xyyyz, g_y_0_y_0_zz_xyyyzz, g_y_0_y_0_zz_xyyzz, g_y_0_y_0_zz_xyyzzz, g_y_0_y_0_zz_xyzzz, g_y_0_y_0_zz_xyzzzz, g_y_0_y_0_zz_xzzzz, g_y_0_y_0_zz_xzzzzz, g_y_0_y_0_zz_yyyyy, g_y_0_y_0_zz_yyyyz, g_y_0_y_0_zz_yyyzz, g_y_0_y_0_zz_yyzzz, g_y_0_y_0_zz_yzzzz, g_y_0_y_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzz_xxxxx[k] = -g_y_0_y_0_zz_xxxxx[k] * ab_x + g_y_0_y_0_zz_xxxxxx[k];

                g_y_0_y_0_xzz_xxxxy[k] = -g_y_0_y_0_zz_xxxxy[k] * ab_x + g_y_0_y_0_zz_xxxxxy[k];

                g_y_0_y_0_xzz_xxxxz[k] = -g_y_0_y_0_zz_xxxxz[k] * ab_x + g_y_0_y_0_zz_xxxxxz[k];

                g_y_0_y_0_xzz_xxxyy[k] = -g_y_0_y_0_zz_xxxyy[k] * ab_x + g_y_0_y_0_zz_xxxxyy[k];

                g_y_0_y_0_xzz_xxxyz[k] = -g_y_0_y_0_zz_xxxyz[k] * ab_x + g_y_0_y_0_zz_xxxxyz[k];

                g_y_0_y_0_xzz_xxxzz[k] = -g_y_0_y_0_zz_xxxzz[k] * ab_x + g_y_0_y_0_zz_xxxxzz[k];

                g_y_0_y_0_xzz_xxyyy[k] = -g_y_0_y_0_zz_xxyyy[k] * ab_x + g_y_0_y_0_zz_xxxyyy[k];

                g_y_0_y_0_xzz_xxyyz[k] = -g_y_0_y_0_zz_xxyyz[k] * ab_x + g_y_0_y_0_zz_xxxyyz[k];

                g_y_0_y_0_xzz_xxyzz[k] = -g_y_0_y_0_zz_xxyzz[k] * ab_x + g_y_0_y_0_zz_xxxyzz[k];

                g_y_0_y_0_xzz_xxzzz[k] = -g_y_0_y_0_zz_xxzzz[k] * ab_x + g_y_0_y_0_zz_xxxzzz[k];

                g_y_0_y_0_xzz_xyyyy[k] = -g_y_0_y_0_zz_xyyyy[k] * ab_x + g_y_0_y_0_zz_xxyyyy[k];

                g_y_0_y_0_xzz_xyyyz[k] = -g_y_0_y_0_zz_xyyyz[k] * ab_x + g_y_0_y_0_zz_xxyyyz[k];

                g_y_0_y_0_xzz_xyyzz[k] = -g_y_0_y_0_zz_xyyzz[k] * ab_x + g_y_0_y_0_zz_xxyyzz[k];

                g_y_0_y_0_xzz_xyzzz[k] = -g_y_0_y_0_zz_xyzzz[k] * ab_x + g_y_0_y_0_zz_xxyzzz[k];

                g_y_0_y_0_xzz_xzzzz[k] = -g_y_0_y_0_zz_xzzzz[k] * ab_x + g_y_0_y_0_zz_xxzzzz[k];

                g_y_0_y_0_xzz_yyyyy[k] = -g_y_0_y_0_zz_yyyyy[k] * ab_x + g_y_0_y_0_zz_xyyyyy[k];

                g_y_0_y_0_xzz_yyyyz[k] = -g_y_0_y_0_zz_yyyyz[k] * ab_x + g_y_0_y_0_zz_xyyyyz[k];

                g_y_0_y_0_xzz_yyyzz[k] = -g_y_0_y_0_zz_yyyzz[k] * ab_x + g_y_0_y_0_zz_xyyyzz[k];

                g_y_0_y_0_xzz_yyzzz[k] = -g_y_0_y_0_zz_yyzzz[k] * ab_x + g_y_0_y_0_zz_xyyzzz[k];

                g_y_0_y_0_xzz_yzzzz[k] = -g_y_0_y_0_zz_yzzzz[k] * ab_x + g_y_0_y_0_zz_xyzzzz[k];

                g_y_0_y_0_xzz_zzzzz[k] = -g_y_0_y_0_zz_zzzzz[k] * ab_x + g_y_0_y_0_zz_xzzzzz[k];
            }

            /// Set up 966-987 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 966 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 967 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 968 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 969 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 970 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 971 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 972 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 973 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 974 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 975 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 976 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 977 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 978 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 979 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 980 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 981 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 982 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 983 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 984 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 985 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 986 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yy_xxxxx, g_0_0_y_0_yy_xxxxy, g_0_0_y_0_yy_xxxxz, g_0_0_y_0_yy_xxxyy, g_0_0_y_0_yy_xxxyz, g_0_0_y_0_yy_xxxzz, g_0_0_y_0_yy_xxyyy, g_0_0_y_0_yy_xxyyz, g_0_0_y_0_yy_xxyzz, g_0_0_y_0_yy_xxzzz, g_0_0_y_0_yy_xyyyy, g_0_0_y_0_yy_xyyyz, g_0_0_y_0_yy_xyyzz, g_0_0_y_0_yy_xyzzz, g_0_0_y_0_yy_xzzzz, g_0_0_y_0_yy_yyyyy, g_0_0_y_0_yy_yyyyz, g_0_0_y_0_yy_yyyzz, g_0_0_y_0_yy_yyzzz, g_0_0_y_0_yy_yzzzz, g_0_0_y_0_yy_zzzzz, g_y_0_y_0_yy_xxxxx, g_y_0_y_0_yy_xxxxxy, g_y_0_y_0_yy_xxxxy, g_y_0_y_0_yy_xxxxyy, g_y_0_y_0_yy_xxxxyz, g_y_0_y_0_yy_xxxxz, g_y_0_y_0_yy_xxxyy, g_y_0_y_0_yy_xxxyyy, g_y_0_y_0_yy_xxxyyz, g_y_0_y_0_yy_xxxyz, g_y_0_y_0_yy_xxxyzz, g_y_0_y_0_yy_xxxzz, g_y_0_y_0_yy_xxyyy, g_y_0_y_0_yy_xxyyyy, g_y_0_y_0_yy_xxyyyz, g_y_0_y_0_yy_xxyyz, g_y_0_y_0_yy_xxyyzz, g_y_0_y_0_yy_xxyzz, g_y_0_y_0_yy_xxyzzz, g_y_0_y_0_yy_xxzzz, g_y_0_y_0_yy_xyyyy, g_y_0_y_0_yy_xyyyyy, g_y_0_y_0_yy_xyyyyz, g_y_0_y_0_yy_xyyyz, g_y_0_y_0_yy_xyyyzz, g_y_0_y_0_yy_xyyzz, g_y_0_y_0_yy_xyyzzz, g_y_0_y_0_yy_xyzzz, g_y_0_y_0_yy_xyzzzz, g_y_0_y_0_yy_xzzzz, g_y_0_y_0_yy_yyyyy, g_y_0_y_0_yy_yyyyyy, g_y_0_y_0_yy_yyyyyz, g_y_0_y_0_yy_yyyyz, g_y_0_y_0_yy_yyyyzz, g_y_0_y_0_yy_yyyzz, g_y_0_y_0_yy_yyyzzz, g_y_0_y_0_yy_yyzzz, g_y_0_y_0_yy_yyzzzz, g_y_0_y_0_yy_yzzzz, g_y_0_y_0_yy_yzzzzz, g_y_0_y_0_yy_zzzzz, g_y_0_y_0_yyy_xxxxx, g_y_0_y_0_yyy_xxxxy, g_y_0_y_0_yyy_xxxxz, g_y_0_y_0_yyy_xxxyy, g_y_0_y_0_yyy_xxxyz, g_y_0_y_0_yyy_xxxzz, g_y_0_y_0_yyy_xxyyy, g_y_0_y_0_yyy_xxyyz, g_y_0_y_0_yyy_xxyzz, g_y_0_y_0_yyy_xxzzz, g_y_0_y_0_yyy_xyyyy, g_y_0_y_0_yyy_xyyyz, g_y_0_y_0_yyy_xyyzz, g_y_0_y_0_yyy_xyzzz, g_y_0_y_0_yyy_xzzzz, g_y_0_y_0_yyy_yyyyy, g_y_0_y_0_yyy_yyyyz, g_y_0_y_0_yyy_yyyzz, g_y_0_y_0_yyy_yyzzz, g_y_0_y_0_yyy_yzzzz, g_y_0_y_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyy_xxxxx[k] = -g_0_0_y_0_yy_xxxxx[k] - g_y_0_y_0_yy_xxxxx[k] * ab_y + g_y_0_y_0_yy_xxxxxy[k];

                g_y_0_y_0_yyy_xxxxy[k] = -g_0_0_y_0_yy_xxxxy[k] - g_y_0_y_0_yy_xxxxy[k] * ab_y + g_y_0_y_0_yy_xxxxyy[k];

                g_y_0_y_0_yyy_xxxxz[k] = -g_0_0_y_0_yy_xxxxz[k] - g_y_0_y_0_yy_xxxxz[k] * ab_y + g_y_0_y_0_yy_xxxxyz[k];

                g_y_0_y_0_yyy_xxxyy[k] = -g_0_0_y_0_yy_xxxyy[k] - g_y_0_y_0_yy_xxxyy[k] * ab_y + g_y_0_y_0_yy_xxxyyy[k];

                g_y_0_y_0_yyy_xxxyz[k] = -g_0_0_y_0_yy_xxxyz[k] - g_y_0_y_0_yy_xxxyz[k] * ab_y + g_y_0_y_0_yy_xxxyyz[k];

                g_y_0_y_0_yyy_xxxzz[k] = -g_0_0_y_0_yy_xxxzz[k] - g_y_0_y_0_yy_xxxzz[k] * ab_y + g_y_0_y_0_yy_xxxyzz[k];

                g_y_0_y_0_yyy_xxyyy[k] = -g_0_0_y_0_yy_xxyyy[k] - g_y_0_y_0_yy_xxyyy[k] * ab_y + g_y_0_y_0_yy_xxyyyy[k];

                g_y_0_y_0_yyy_xxyyz[k] = -g_0_0_y_0_yy_xxyyz[k] - g_y_0_y_0_yy_xxyyz[k] * ab_y + g_y_0_y_0_yy_xxyyyz[k];

                g_y_0_y_0_yyy_xxyzz[k] = -g_0_0_y_0_yy_xxyzz[k] - g_y_0_y_0_yy_xxyzz[k] * ab_y + g_y_0_y_0_yy_xxyyzz[k];

                g_y_0_y_0_yyy_xxzzz[k] = -g_0_0_y_0_yy_xxzzz[k] - g_y_0_y_0_yy_xxzzz[k] * ab_y + g_y_0_y_0_yy_xxyzzz[k];

                g_y_0_y_0_yyy_xyyyy[k] = -g_0_0_y_0_yy_xyyyy[k] - g_y_0_y_0_yy_xyyyy[k] * ab_y + g_y_0_y_0_yy_xyyyyy[k];

                g_y_0_y_0_yyy_xyyyz[k] = -g_0_0_y_0_yy_xyyyz[k] - g_y_0_y_0_yy_xyyyz[k] * ab_y + g_y_0_y_0_yy_xyyyyz[k];

                g_y_0_y_0_yyy_xyyzz[k] = -g_0_0_y_0_yy_xyyzz[k] - g_y_0_y_0_yy_xyyzz[k] * ab_y + g_y_0_y_0_yy_xyyyzz[k];

                g_y_0_y_0_yyy_xyzzz[k] = -g_0_0_y_0_yy_xyzzz[k] - g_y_0_y_0_yy_xyzzz[k] * ab_y + g_y_0_y_0_yy_xyyzzz[k];

                g_y_0_y_0_yyy_xzzzz[k] = -g_0_0_y_0_yy_xzzzz[k] - g_y_0_y_0_yy_xzzzz[k] * ab_y + g_y_0_y_0_yy_xyzzzz[k];

                g_y_0_y_0_yyy_yyyyy[k] = -g_0_0_y_0_yy_yyyyy[k] - g_y_0_y_0_yy_yyyyy[k] * ab_y + g_y_0_y_0_yy_yyyyyy[k];

                g_y_0_y_0_yyy_yyyyz[k] = -g_0_0_y_0_yy_yyyyz[k] - g_y_0_y_0_yy_yyyyz[k] * ab_y + g_y_0_y_0_yy_yyyyyz[k];

                g_y_0_y_0_yyy_yyyzz[k] = -g_0_0_y_0_yy_yyyzz[k] - g_y_0_y_0_yy_yyyzz[k] * ab_y + g_y_0_y_0_yy_yyyyzz[k];

                g_y_0_y_0_yyy_yyzzz[k] = -g_0_0_y_0_yy_yyzzz[k] - g_y_0_y_0_yy_yyzzz[k] * ab_y + g_y_0_y_0_yy_yyyzzz[k];

                g_y_0_y_0_yyy_yzzzz[k] = -g_0_0_y_0_yy_yzzzz[k] - g_y_0_y_0_yy_yzzzz[k] * ab_y + g_y_0_y_0_yy_yyzzzz[k];

                g_y_0_y_0_yyy_zzzzz[k] = -g_0_0_y_0_yy_zzzzz[k] - g_y_0_y_0_yy_zzzzz[k] * ab_y + g_y_0_y_0_yy_yzzzzz[k];
            }

            /// Set up 987-1008 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 987 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 988 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 989 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 990 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 991 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 992 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 993 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 994 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 995 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 996 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 997 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 998 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 999 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1000 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1001 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1002 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1003 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1004 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1005 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1006 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yy_xxxxx, g_y_0_y_0_yy_xxxxxz, g_y_0_y_0_yy_xxxxy, g_y_0_y_0_yy_xxxxyz, g_y_0_y_0_yy_xxxxz, g_y_0_y_0_yy_xxxxzz, g_y_0_y_0_yy_xxxyy, g_y_0_y_0_yy_xxxyyz, g_y_0_y_0_yy_xxxyz, g_y_0_y_0_yy_xxxyzz, g_y_0_y_0_yy_xxxzz, g_y_0_y_0_yy_xxxzzz, g_y_0_y_0_yy_xxyyy, g_y_0_y_0_yy_xxyyyz, g_y_0_y_0_yy_xxyyz, g_y_0_y_0_yy_xxyyzz, g_y_0_y_0_yy_xxyzz, g_y_0_y_0_yy_xxyzzz, g_y_0_y_0_yy_xxzzz, g_y_0_y_0_yy_xxzzzz, g_y_0_y_0_yy_xyyyy, g_y_0_y_0_yy_xyyyyz, g_y_0_y_0_yy_xyyyz, g_y_0_y_0_yy_xyyyzz, g_y_0_y_0_yy_xyyzz, g_y_0_y_0_yy_xyyzzz, g_y_0_y_0_yy_xyzzz, g_y_0_y_0_yy_xyzzzz, g_y_0_y_0_yy_xzzzz, g_y_0_y_0_yy_xzzzzz, g_y_0_y_0_yy_yyyyy, g_y_0_y_0_yy_yyyyyz, g_y_0_y_0_yy_yyyyz, g_y_0_y_0_yy_yyyyzz, g_y_0_y_0_yy_yyyzz, g_y_0_y_0_yy_yyyzzz, g_y_0_y_0_yy_yyzzz, g_y_0_y_0_yy_yyzzzz, g_y_0_y_0_yy_yzzzz, g_y_0_y_0_yy_yzzzzz, g_y_0_y_0_yy_zzzzz, g_y_0_y_0_yy_zzzzzz, g_y_0_y_0_yyz_xxxxx, g_y_0_y_0_yyz_xxxxy, g_y_0_y_0_yyz_xxxxz, g_y_0_y_0_yyz_xxxyy, g_y_0_y_0_yyz_xxxyz, g_y_0_y_0_yyz_xxxzz, g_y_0_y_0_yyz_xxyyy, g_y_0_y_0_yyz_xxyyz, g_y_0_y_0_yyz_xxyzz, g_y_0_y_0_yyz_xxzzz, g_y_0_y_0_yyz_xyyyy, g_y_0_y_0_yyz_xyyyz, g_y_0_y_0_yyz_xyyzz, g_y_0_y_0_yyz_xyzzz, g_y_0_y_0_yyz_xzzzz, g_y_0_y_0_yyz_yyyyy, g_y_0_y_0_yyz_yyyyz, g_y_0_y_0_yyz_yyyzz, g_y_0_y_0_yyz_yyzzz, g_y_0_y_0_yyz_yzzzz, g_y_0_y_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyz_xxxxx[k] = -g_y_0_y_0_yy_xxxxx[k] * ab_z + g_y_0_y_0_yy_xxxxxz[k];

                g_y_0_y_0_yyz_xxxxy[k] = -g_y_0_y_0_yy_xxxxy[k] * ab_z + g_y_0_y_0_yy_xxxxyz[k];

                g_y_0_y_0_yyz_xxxxz[k] = -g_y_0_y_0_yy_xxxxz[k] * ab_z + g_y_0_y_0_yy_xxxxzz[k];

                g_y_0_y_0_yyz_xxxyy[k] = -g_y_0_y_0_yy_xxxyy[k] * ab_z + g_y_0_y_0_yy_xxxyyz[k];

                g_y_0_y_0_yyz_xxxyz[k] = -g_y_0_y_0_yy_xxxyz[k] * ab_z + g_y_0_y_0_yy_xxxyzz[k];

                g_y_0_y_0_yyz_xxxzz[k] = -g_y_0_y_0_yy_xxxzz[k] * ab_z + g_y_0_y_0_yy_xxxzzz[k];

                g_y_0_y_0_yyz_xxyyy[k] = -g_y_0_y_0_yy_xxyyy[k] * ab_z + g_y_0_y_0_yy_xxyyyz[k];

                g_y_0_y_0_yyz_xxyyz[k] = -g_y_0_y_0_yy_xxyyz[k] * ab_z + g_y_0_y_0_yy_xxyyzz[k];

                g_y_0_y_0_yyz_xxyzz[k] = -g_y_0_y_0_yy_xxyzz[k] * ab_z + g_y_0_y_0_yy_xxyzzz[k];

                g_y_0_y_0_yyz_xxzzz[k] = -g_y_0_y_0_yy_xxzzz[k] * ab_z + g_y_0_y_0_yy_xxzzzz[k];

                g_y_0_y_0_yyz_xyyyy[k] = -g_y_0_y_0_yy_xyyyy[k] * ab_z + g_y_0_y_0_yy_xyyyyz[k];

                g_y_0_y_0_yyz_xyyyz[k] = -g_y_0_y_0_yy_xyyyz[k] * ab_z + g_y_0_y_0_yy_xyyyzz[k];

                g_y_0_y_0_yyz_xyyzz[k] = -g_y_0_y_0_yy_xyyzz[k] * ab_z + g_y_0_y_0_yy_xyyzzz[k];

                g_y_0_y_0_yyz_xyzzz[k] = -g_y_0_y_0_yy_xyzzz[k] * ab_z + g_y_0_y_0_yy_xyzzzz[k];

                g_y_0_y_0_yyz_xzzzz[k] = -g_y_0_y_0_yy_xzzzz[k] * ab_z + g_y_0_y_0_yy_xzzzzz[k];

                g_y_0_y_0_yyz_yyyyy[k] = -g_y_0_y_0_yy_yyyyy[k] * ab_z + g_y_0_y_0_yy_yyyyyz[k];

                g_y_0_y_0_yyz_yyyyz[k] = -g_y_0_y_0_yy_yyyyz[k] * ab_z + g_y_0_y_0_yy_yyyyzz[k];

                g_y_0_y_0_yyz_yyyzz[k] = -g_y_0_y_0_yy_yyyzz[k] * ab_z + g_y_0_y_0_yy_yyyzzz[k];

                g_y_0_y_0_yyz_yyzzz[k] = -g_y_0_y_0_yy_yyzzz[k] * ab_z + g_y_0_y_0_yy_yyzzzz[k];

                g_y_0_y_0_yyz_yzzzz[k] = -g_y_0_y_0_yy_yzzzz[k] * ab_z + g_y_0_y_0_yy_yzzzzz[k];

                g_y_0_y_0_yyz_zzzzz[k] = -g_y_0_y_0_yy_zzzzz[k] * ab_z + g_y_0_y_0_yy_zzzzzz[k];
            }

            /// Set up 1008-1029 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1008 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1009 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1010 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1011 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1012 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1013 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1014 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1015 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1016 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1017 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1018 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1019 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1020 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1021 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1022 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1023 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1024 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1025 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1026 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1027 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1028 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yz_xxxxx, g_y_0_y_0_yz_xxxxxz, g_y_0_y_0_yz_xxxxy, g_y_0_y_0_yz_xxxxyz, g_y_0_y_0_yz_xxxxz, g_y_0_y_0_yz_xxxxzz, g_y_0_y_0_yz_xxxyy, g_y_0_y_0_yz_xxxyyz, g_y_0_y_0_yz_xxxyz, g_y_0_y_0_yz_xxxyzz, g_y_0_y_0_yz_xxxzz, g_y_0_y_0_yz_xxxzzz, g_y_0_y_0_yz_xxyyy, g_y_0_y_0_yz_xxyyyz, g_y_0_y_0_yz_xxyyz, g_y_0_y_0_yz_xxyyzz, g_y_0_y_0_yz_xxyzz, g_y_0_y_0_yz_xxyzzz, g_y_0_y_0_yz_xxzzz, g_y_0_y_0_yz_xxzzzz, g_y_0_y_0_yz_xyyyy, g_y_0_y_0_yz_xyyyyz, g_y_0_y_0_yz_xyyyz, g_y_0_y_0_yz_xyyyzz, g_y_0_y_0_yz_xyyzz, g_y_0_y_0_yz_xyyzzz, g_y_0_y_0_yz_xyzzz, g_y_0_y_0_yz_xyzzzz, g_y_0_y_0_yz_xzzzz, g_y_0_y_0_yz_xzzzzz, g_y_0_y_0_yz_yyyyy, g_y_0_y_0_yz_yyyyyz, g_y_0_y_0_yz_yyyyz, g_y_0_y_0_yz_yyyyzz, g_y_0_y_0_yz_yyyzz, g_y_0_y_0_yz_yyyzzz, g_y_0_y_0_yz_yyzzz, g_y_0_y_0_yz_yyzzzz, g_y_0_y_0_yz_yzzzz, g_y_0_y_0_yz_yzzzzz, g_y_0_y_0_yz_zzzzz, g_y_0_y_0_yz_zzzzzz, g_y_0_y_0_yzz_xxxxx, g_y_0_y_0_yzz_xxxxy, g_y_0_y_0_yzz_xxxxz, g_y_0_y_0_yzz_xxxyy, g_y_0_y_0_yzz_xxxyz, g_y_0_y_0_yzz_xxxzz, g_y_0_y_0_yzz_xxyyy, g_y_0_y_0_yzz_xxyyz, g_y_0_y_0_yzz_xxyzz, g_y_0_y_0_yzz_xxzzz, g_y_0_y_0_yzz_xyyyy, g_y_0_y_0_yzz_xyyyz, g_y_0_y_0_yzz_xyyzz, g_y_0_y_0_yzz_xyzzz, g_y_0_y_0_yzz_xzzzz, g_y_0_y_0_yzz_yyyyy, g_y_0_y_0_yzz_yyyyz, g_y_0_y_0_yzz_yyyzz, g_y_0_y_0_yzz_yyzzz, g_y_0_y_0_yzz_yzzzz, g_y_0_y_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzz_xxxxx[k] = -g_y_0_y_0_yz_xxxxx[k] * ab_z + g_y_0_y_0_yz_xxxxxz[k];

                g_y_0_y_0_yzz_xxxxy[k] = -g_y_0_y_0_yz_xxxxy[k] * ab_z + g_y_0_y_0_yz_xxxxyz[k];

                g_y_0_y_0_yzz_xxxxz[k] = -g_y_0_y_0_yz_xxxxz[k] * ab_z + g_y_0_y_0_yz_xxxxzz[k];

                g_y_0_y_0_yzz_xxxyy[k] = -g_y_0_y_0_yz_xxxyy[k] * ab_z + g_y_0_y_0_yz_xxxyyz[k];

                g_y_0_y_0_yzz_xxxyz[k] = -g_y_0_y_0_yz_xxxyz[k] * ab_z + g_y_0_y_0_yz_xxxyzz[k];

                g_y_0_y_0_yzz_xxxzz[k] = -g_y_0_y_0_yz_xxxzz[k] * ab_z + g_y_0_y_0_yz_xxxzzz[k];

                g_y_0_y_0_yzz_xxyyy[k] = -g_y_0_y_0_yz_xxyyy[k] * ab_z + g_y_0_y_0_yz_xxyyyz[k];

                g_y_0_y_0_yzz_xxyyz[k] = -g_y_0_y_0_yz_xxyyz[k] * ab_z + g_y_0_y_0_yz_xxyyzz[k];

                g_y_0_y_0_yzz_xxyzz[k] = -g_y_0_y_0_yz_xxyzz[k] * ab_z + g_y_0_y_0_yz_xxyzzz[k];

                g_y_0_y_0_yzz_xxzzz[k] = -g_y_0_y_0_yz_xxzzz[k] * ab_z + g_y_0_y_0_yz_xxzzzz[k];

                g_y_0_y_0_yzz_xyyyy[k] = -g_y_0_y_0_yz_xyyyy[k] * ab_z + g_y_0_y_0_yz_xyyyyz[k];

                g_y_0_y_0_yzz_xyyyz[k] = -g_y_0_y_0_yz_xyyyz[k] * ab_z + g_y_0_y_0_yz_xyyyzz[k];

                g_y_0_y_0_yzz_xyyzz[k] = -g_y_0_y_0_yz_xyyzz[k] * ab_z + g_y_0_y_0_yz_xyyzzz[k];

                g_y_0_y_0_yzz_xyzzz[k] = -g_y_0_y_0_yz_xyzzz[k] * ab_z + g_y_0_y_0_yz_xyzzzz[k];

                g_y_0_y_0_yzz_xzzzz[k] = -g_y_0_y_0_yz_xzzzz[k] * ab_z + g_y_0_y_0_yz_xzzzzz[k];

                g_y_0_y_0_yzz_yyyyy[k] = -g_y_0_y_0_yz_yyyyy[k] * ab_z + g_y_0_y_0_yz_yyyyyz[k];

                g_y_0_y_0_yzz_yyyyz[k] = -g_y_0_y_0_yz_yyyyz[k] * ab_z + g_y_0_y_0_yz_yyyyzz[k];

                g_y_0_y_0_yzz_yyyzz[k] = -g_y_0_y_0_yz_yyyzz[k] * ab_z + g_y_0_y_0_yz_yyyzzz[k];

                g_y_0_y_0_yzz_yyzzz[k] = -g_y_0_y_0_yz_yyzzz[k] * ab_z + g_y_0_y_0_yz_yyzzzz[k];

                g_y_0_y_0_yzz_yzzzz[k] = -g_y_0_y_0_yz_yzzzz[k] * ab_z + g_y_0_y_0_yz_yzzzzz[k];

                g_y_0_y_0_yzz_zzzzz[k] = -g_y_0_y_0_yz_zzzzz[k] * ab_z + g_y_0_y_0_yz_zzzzzz[k];
            }

            /// Set up 1029-1050 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1029 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1030 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1031 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1032 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1033 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1034 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1035 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1036 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1037 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1038 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1039 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1040 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1041 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1042 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1043 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1044 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1045 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1046 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1047 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1048 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zz_xxxxx, g_y_0_y_0_zz_xxxxxz, g_y_0_y_0_zz_xxxxy, g_y_0_y_0_zz_xxxxyz, g_y_0_y_0_zz_xxxxz, g_y_0_y_0_zz_xxxxzz, g_y_0_y_0_zz_xxxyy, g_y_0_y_0_zz_xxxyyz, g_y_0_y_0_zz_xxxyz, g_y_0_y_0_zz_xxxyzz, g_y_0_y_0_zz_xxxzz, g_y_0_y_0_zz_xxxzzz, g_y_0_y_0_zz_xxyyy, g_y_0_y_0_zz_xxyyyz, g_y_0_y_0_zz_xxyyz, g_y_0_y_0_zz_xxyyzz, g_y_0_y_0_zz_xxyzz, g_y_0_y_0_zz_xxyzzz, g_y_0_y_0_zz_xxzzz, g_y_0_y_0_zz_xxzzzz, g_y_0_y_0_zz_xyyyy, g_y_0_y_0_zz_xyyyyz, g_y_0_y_0_zz_xyyyz, g_y_0_y_0_zz_xyyyzz, g_y_0_y_0_zz_xyyzz, g_y_0_y_0_zz_xyyzzz, g_y_0_y_0_zz_xyzzz, g_y_0_y_0_zz_xyzzzz, g_y_0_y_0_zz_xzzzz, g_y_0_y_0_zz_xzzzzz, g_y_0_y_0_zz_yyyyy, g_y_0_y_0_zz_yyyyyz, g_y_0_y_0_zz_yyyyz, g_y_0_y_0_zz_yyyyzz, g_y_0_y_0_zz_yyyzz, g_y_0_y_0_zz_yyyzzz, g_y_0_y_0_zz_yyzzz, g_y_0_y_0_zz_yyzzzz, g_y_0_y_0_zz_yzzzz, g_y_0_y_0_zz_yzzzzz, g_y_0_y_0_zz_zzzzz, g_y_0_y_0_zz_zzzzzz, g_y_0_y_0_zzz_xxxxx, g_y_0_y_0_zzz_xxxxy, g_y_0_y_0_zzz_xxxxz, g_y_0_y_0_zzz_xxxyy, g_y_0_y_0_zzz_xxxyz, g_y_0_y_0_zzz_xxxzz, g_y_0_y_0_zzz_xxyyy, g_y_0_y_0_zzz_xxyyz, g_y_0_y_0_zzz_xxyzz, g_y_0_y_0_zzz_xxzzz, g_y_0_y_0_zzz_xyyyy, g_y_0_y_0_zzz_xyyyz, g_y_0_y_0_zzz_xyyzz, g_y_0_y_0_zzz_xyzzz, g_y_0_y_0_zzz_xzzzz, g_y_0_y_0_zzz_yyyyy, g_y_0_y_0_zzz_yyyyz, g_y_0_y_0_zzz_yyyzz, g_y_0_y_0_zzz_yyzzz, g_y_0_y_0_zzz_yzzzz, g_y_0_y_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzz_xxxxx[k] = -g_y_0_y_0_zz_xxxxx[k] * ab_z + g_y_0_y_0_zz_xxxxxz[k];

                g_y_0_y_0_zzz_xxxxy[k] = -g_y_0_y_0_zz_xxxxy[k] * ab_z + g_y_0_y_0_zz_xxxxyz[k];

                g_y_0_y_0_zzz_xxxxz[k] = -g_y_0_y_0_zz_xxxxz[k] * ab_z + g_y_0_y_0_zz_xxxxzz[k];

                g_y_0_y_0_zzz_xxxyy[k] = -g_y_0_y_0_zz_xxxyy[k] * ab_z + g_y_0_y_0_zz_xxxyyz[k];

                g_y_0_y_0_zzz_xxxyz[k] = -g_y_0_y_0_zz_xxxyz[k] * ab_z + g_y_0_y_0_zz_xxxyzz[k];

                g_y_0_y_0_zzz_xxxzz[k] = -g_y_0_y_0_zz_xxxzz[k] * ab_z + g_y_0_y_0_zz_xxxzzz[k];

                g_y_0_y_0_zzz_xxyyy[k] = -g_y_0_y_0_zz_xxyyy[k] * ab_z + g_y_0_y_0_zz_xxyyyz[k];

                g_y_0_y_0_zzz_xxyyz[k] = -g_y_0_y_0_zz_xxyyz[k] * ab_z + g_y_0_y_0_zz_xxyyzz[k];

                g_y_0_y_0_zzz_xxyzz[k] = -g_y_0_y_0_zz_xxyzz[k] * ab_z + g_y_0_y_0_zz_xxyzzz[k];

                g_y_0_y_0_zzz_xxzzz[k] = -g_y_0_y_0_zz_xxzzz[k] * ab_z + g_y_0_y_0_zz_xxzzzz[k];

                g_y_0_y_0_zzz_xyyyy[k] = -g_y_0_y_0_zz_xyyyy[k] * ab_z + g_y_0_y_0_zz_xyyyyz[k];

                g_y_0_y_0_zzz_xyyyz[k] = -g_y_0_y_0_zz_xyyyz[k] * ab_z + g_y_0_y_0_zz_xyyyzz[k];

                g_y_0_y_0_zzz_xyyzz[k] = -g_y_0_y_0_zz_xyyzz[k] * ab_z + g_y_0_y_0_zz_xyyzzz[k];

                g_y_0_y_0_zzz_xyzzz[k] = -g_y_0_y_0_zz_xyzzz[k] * ab_z + g_y_0_y_0_zz_xyzzzz[k];

                g_y_0_y_0_zzz_xzzzz[k] = -g_y_0_y_0_zz_xzzzz[k] * ab_z + g_y_0_y_0_zz_xzzzzz[k];

                g_y_0_y_0_zzz_yyyyy[k] = -g_y_0_y_0_zz_yyyyy[k] * ab_z + g_y_0_y_0_zz_yyyyyz[k];

                g_y_0_y_0_zzz_yyyyz[k] = -g_y_0_y_0_zz_yyyyz[k] * ab_z + g_y_0_y_0_zz_yyyyzz[k];

                g_y_0_y_0_zzz_yyyzz[k] = -g_y_0_y_0_zz_yyyzz[k] * ab_z + g_y_0_y_0_zz_yyyzzz[k];

                g_y_0_y_0_zzz_yyzzz[k] = -g_y_0_y_0_zz_yyzzz[k] * ab_z + g_y_0_y_0_zz_yyzzzz[k];

                g_y_0_y_0_zzz_yzzzz[k] = -g_y_0_y_0_zz_yzzzz[k] * ab_z + g_y_0_y_0_zz_yzzzzz[k];

                g_y_0_y_0_zzz_zzzzz[k] = -g_y_0_y_0_zz_zzzzz[k] * ab_z + g_y_0_y_0_zz_zzzzzz[k];
            }

            /// Set up 1050-1071 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 1050 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 1051 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 1052 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 1053 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 1054 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 1055 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 1056 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 1057 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 1058 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 1059 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 1060 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 1061 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 1062 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 1063 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 1064 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 1065 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 1066 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 1067 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 1068 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 1069 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 1070 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xx_xxxxx, g_y_0_z_0_xx_xxxxxx, g_y_0_z_0_xx_xxxxxy, g_y_0_z_0_xx_xxxxxz, g_y_0_z_0_xx_xxxxy, g_y_0_z_0_xx_xxxxyy, g_y_0_z_0_xx_xxxxyz, g_y_0_z_0_xx_xxxxz, g_y_0_z_0_xx_xxxxzz, g_y_0_z_0_xx_xxxyy, g_y_0_z_0_xx_xxxyyy, g_y_0_z_0_xx_xxxyyz, g_y_0_z_0_xx_xxxyz, g_y_0_z_0_xx_xxxyzz, g_y_0_z_0_xx_xxxzz, g_y_0_z_0_xx_xxxzzz, g_y_0_z_0_xx_xxyyy, g_y_0_z_0_xx_xxyyyy, g_y_0_z_0_xx_xxyyyz, g_y_0_z_0_xx_xxyyz, g_y_0_z_0_xx_xxyyzz, g_y_0_z_0_xx_xxyzz, g_y_0_z_0_xx_xxyzzz, g_y_0_z_0_xx_xxzzz, g_y_0_z_0_xx_xxzzzz, g_y_0_z_0_xx_xyyyy, g_y_0_z_0_xx_xyyyyy, g_y_0_z_0_xx_xyyyyz, g_y_0_z_0_xx_xyyyz, g_y_0_z_0_xx_xyyyzz, g_y_0_z_0_xx_xyyzz, g_y_0_z_0_xx_xyyzzz, g_y_0_z_0_xx_xyzzz, g_y_0_z_0_xx_xyzzzz, g_y_0_z_0_xx_xzzzz, g_y_0_z_0_xx_xzzzzz, g_y_0_z_0_xx_yyyyy, g_y_0_z_0_xx_yyyyz, g_y_0_z_0_xx_yyyzz, g_y_0_z_0_xx_yyzzz, g_y_0_z_0_xx_yzzzz, g_y_0_z_0_xx_zzzzz, g_y_0_z_0_xxx_xxxxx, g_y_0_z_0_xxx_xxxxy, g_y_0_z_0_xxx_xxxxz, g_y_0_z_0_xxx_xxxyy, g_y_0_z_0_xxx_xxxyz, g_y_0_z_0_xxx_xxxzz, g_y_0_z_0_xxx_xxyyy, g_y_0_z_0_xxx_xxyyz, g_y_0_z_0_xxx_xxyzz, g_y_0_z_0_xxx_xxzzz, g_y_0_z_0_xxx_xyyyy, g_y_0_z_0_xxx_xyyyz, g_y_0_z_0_xxx_xyyzz, g_y_0_z_0_xxx_xyzzz, g_y_0_z_0_xxx_xzzzz, g_y_0_z_0_xxx_yyyyy, g_y_0_z_0_xxx_yyyyz, g_y_0_z_0_xxx_yyyzz, g_y_0_z_0_xxx_yyzzz, g_y_0_z_0_xxx_yzzzz, g_y_0_z_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxx_xxxxx[k] = -g_y_0_z_0_xx_xxxxx[k] * ab_x + g_y_0_z_0_xx_xxxxxx[k];

                g_y_0_z_0_xxx_xxxxy[k] = -g_y_0_z_0_xx_xxxxy[k] * ab_x + g_y_0_z_0_xx_xxxxxy[k];

                g_y_0_z_0_xxx_xxxxz[k] = -g_y_0_z_0_xx_xxxxz[k] * ab_x + g_y_0_z_0_xx_xxxxxz[k];

                g_y_0_z_0_xxx_xxxyy[k] = -g_y_0_z_0_xx_xxxyy[k] * ab_x + g_y_0_z_0_xx_xxxxyy[k];

                g_y_0_z_0_xxx_xxxyz[k] = -g_y_0_z_0_xx_xxxyz[k] * ab_x + g_y_0_z_0_xx_xxxxyz[k];

                g_y_0_z_0_xxx_xxxzz[k] = -g_y_0_z_0_xx_xxxzz[k] * ab_x + g_y_0_z_0_xx_xxxxzz[k];

                g_y_0_z_0_xxx_xxyyy[k] = -g_y_0_z_0_xx_xxyyy[k] * ab_x + g_y_0_z_0_xx_xxxyyy[k];

                g_y_0_z_0_xxx_xxyyz[k] = -g_y_0_z_0_xx_xxyyz[k] * ab_x + g_y_0_z_0_xx_xxxyyz[k];

                g_y_0_z_0_xxx_xxyzz[k] = -g_y_0_z_0_xx_xxyzz[k] * ab_x + g_y_0_z_0_xx_xxxyzz[k];

                g_y_0_z_0_xxx_xxzzz[k] = -g_y_0_z_0_xx_xxzzz[k] * ab_x + g_y_0_z_0_xx_xxxzzz[k];

                g_y_0_z_0_xxx_xyyyy[k] = -g_y_0_z_0_xx_xyyyy[k] * ab_x + g_y_0_z_0_xx_xxyyyy[k];

                g_y_0_z_0_xxx_xyyyz[k] = -g_y_0_z_0_xx_xyyyz[k] * ab_x + g_y_0_z_0_xx_xxyyyz[k];

                g_y_0_z_0_xxx_xyyzz[k] = -g_y_0_z_0_xx_xyyzz[k] * ab_x + g_y_0_z_0_xx_xxyyzz[k];

                g_y_0_z_0_xxx_xyzzz[k] = -g_y_0_z_0_xx_xyzzz[k] * ab_x + g_y_0_z_0_xx_xxyzzz[k];

                g_y_0_z_0_xxx_xzzzz[k] = -g_y_0_z_0_xx_xzzzz[k] * ab_x + g_y_0_z_0_xx_xxzzzz[k];

                g_y_0_z_0_xxx_yyyyy[k] = -g_y_0_z_0_xx_yyyyy[k] * ab_x + g_y_0_z_0_xx_xyyyyy[k];

                g_y_0_z_0_xxx_yyyyz[k] = -g_y_0_z_0_xx_yyyyz[k] * ab_x + g_y_0_z_0_xx_xyyyyz[k];

                g_y_0_z_0_xxx_yyyzz[k] = -g_y_0_z_0_xx_yyyzz[k] * ab_x + g_y_0_z_0_xx_xyyyzz[k];

                g_y_0_z_0_xxx_yyzzz[k] = -g_y_0_z_0_xx_yyzzz[k] * ab_x + g_y_0_z_0_xx_xyyzzz[k];

                g_y_0_z_0_xxx_yzzzz[k] = -g_y_0_z_0_xx_yzzzz[k] * ab_x + g_y_0_z_0_xx_xyzzzz[k];

                g_y_0_z_0_xxx_zzzzz[k] = -g_y_0_z_0_xx_zzzzz[k] * ab_x + g_y_0_z_0_xx_xzzzzz[k];
            }

            /// Set up 1071-1092 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 1071 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 1072 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 1073 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 1074 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 1075 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 1076 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 1077 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 1078 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 1079 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 1080 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 1081 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 1082 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 1083 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 1084 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 1085 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 1086 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 1087 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 1088 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 1089 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 1090 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxy_xxxxx, g_y_0_z_0_xxy_xxxxy, g_y_0_z_0_xxy_xxxxz, g_y_0_z_0_xxy_xxxyy, g_y_0_z_0_xxy_xxxyz, g_y_0_z_0_xxy_xxxzz, g_y_0_z_0_xxy_xxyyy, g_y_0_z_0_xxy_xxyyz, g_y_0_z_0_xxy_xxyzz, g_y_0_z_0_xxy_xxzzz, g_y_0_z_0_xxy_xyyyy, g_y_0_z_0_xxy_xyyyz, g_y_0_z_0_xxy_xyyzz, g_y_0_z_0_xxy_xyzzz, g_y_0_z_0_xxy_xzzzz, g_y_0_z_0_xxy_yyyyy, g_y_0_z_0_xxy_yyyyz, g_y_0_z_0_xxy_yyyzz, g_y_0_z_0_xxy_yyzzz, g_y_0_z_0_xxy_yzzzz, g_y_0_z_0_xxy_zzzzz, g_y_0_z_0_xy_xxxxx, g_y_0_z_0_xy_xxxxxx, g_y_0_z_0_xy_xxxxxy, g_y_0_z_0_xy_xxxxxz, g_y_0_z_0_xy_xxxxy, g_y_0_z_0_xy_xxxxyy, g_y_0_z_0_xy_xxxxyz, g_y_0_z_0_xy_xxxxz, g_y_0_z_0_xy_xxxxzz, g_y_0_z_0_xy_xxxyy, g_y_0_z_0_xy_xxxyyy, g_y_0_z_0_xy_xxxyyz, g_y_0_z_0_xy_xxxyz, g_y_0_z_0_xy_xxxyzz, g_y_0_z_0_xy_xxxzz, g_y_0_z_0_xy_xxxzzz, g_y_0_z_0_xy_xxyyy, g_y_0_z_0_xy_xxyyyy, g_y_0_z_0_xy_xxyyyz, g_y_0_z_0_xy_xxyyz, g_y_0_z_0_xy_xxyyzz, g_y_0_z_0_xy_xxyzz, g_y_0_z_0_xy_xxyzzz, g_y_0_z_0_xy_xxzzz, g_y_0_z_0_xy_xxzzzz, g_y_0_z_0_xy_xyyyy, g_y_0_z_0_xy_xyyyyy, g_y_0_z_0_xy_xyyyyz, g_y_0_z_0_xy_xyyyz, g_y_0_z_0_xy_xyyyzz, g_y_0_z_0_xy_xyyzz, g_y_0_z_0_xy_xyyzzz, g_y_0_z_0_xy_xyzzz, g_y_0_z_0_xy_xyzzzz, g_y_0_z_0_xy_xzzzz, g_y_0_z_0_xy_xzzzzz, g_y_0_z_0_xy_yyyyy, g_y_0_z_0_xy_yyyyz, g_y_0_z_0_xy_yyyzz, g_y_0_z_0_xy_yyzzz, g_y_0_z_0_xy_yzzzz, g_y_0_z_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxy_xxxxx[k] = -g_y_0_z_0_xy_xxxxx[k] * ab_x + g_y_0_z_0_xy_xxxxxx[k];

                g_y_0_z_0_xxy_xxxxy[k] = -g_y_0_z_0_xy_xxxxy[k] * ab_x + g_y_0_z_0_xy_xxxxxy[k];

                g_y_0_z_0_xxy_xxxxz[k] = -g_y_0_z_0_xy_xxxxz[k] * ab_x + g_y_0_z_0_xy_xxxxxz[k];

                g_y_0_z_0_xxy_xxxyy[k] = -g_y_0_z_0_xy_xxxyy[k] * ab_x + g_y_0_z_0_xy_xxxxyy[k];

                g_y_0_z_0_xxy_xxxyz[k] = -g_y_0_z_0_xy_xxxyz[k] * ab_x + g_y_0_z_0_xy_xxxxyz[k];

                g_y_0_z_0_xxy_xxxzz[k] = -g_y_0_z_0_xy_xxxzz[k] * ab_x + g_y_0_z_0_xy_xxxxzz[k];

                g_y_0_z_0_xxy_xxyyy[k] = -g_y_0_z_0_xy_xxyyy[k] * ab_x + g_y_0_z_0_xy_xxxyyy[k];

                g_y_0_z_0_xxy_xxyyz[k] = -g_y_0_z_0_xy_xxyyz[k] * ab_x + g_y_0_z_0_xy_xxxyyz[k];

                g_y_0_z_0_xxy_xxyzz[k] = -g_y_0_z_0_xy_xxyzz[k] * ab_x + g_y_0_z_0_xy_xxxyzz[k];

                g_y_0_z_0_xxy_xxzzz[k] = -g_y_0_z_0_xy_xxzzz[k] * ab_x + g_y_0_z_0_xy_xxxzzz[k];

                g_y_0_z_0_xxy_xyyyy[k] = -g_y_0_z_0_xy_xyyyy[k] * ab_x + g_y_0_z_0_xy_xxyyyy[k];

                g_y_0_z_0_xxy_xyyyz[k] = -g_y_0_z_0_xy_xyyyz[k] * ab_x + g_y_0_z_0_xy_xxyyyz[k];

                g_y_0_z_0_xxy_xyyzz[k] = -g_y_0_z_0_xy_xyyzz[k] * ab_x + g_y_0_z_0_xy_xxyyzz[k];

                g_y_0_z_0_xxy_xyzzz[k] = -g_y_0_z_0_xy_xyzzz[k] * ab_x + g_y_0_z_0_xy_xxyzzz[k];

                g_y_0_z_0_xxy_xzzzz[k] = -g_y_0_z_0_xy_xzzzz[k] * ab_x + g_y_0_z_0_xy_xxzzzz[k];

                g_y_0_z_0_xxy_yyyyy[k] = -g_y_0_z_0_xy_yyyyy[k] * ab_x + g_y_0_z_0_xy_xyyyyy[k];

                g_y_0_z_0_xxy_yyyyz[k] = -g_y_0_z_0_xy_yyyyz[k] * ab_x + g_y_0_z_0_xy_xyyyyz[k];

                g_y_0_z_0_xxy_yyyzz[k] = -g_y_0_z_0_xy_yyyzz[k] * ab_x + g_y_0_z_0_xy_xyyyzz[k];

                g_y_0_z_0_xxy_yyzzz[k] = -g_y_0_z_0_xy_yyzzz[k] * ab_x + g_y_0_z_0_xy_xyyzzz[k];

                g_y_0_z_0_xxy_yzzzz[k] = -g_y_0_z_0_xy_yzzzz[k] * ab_x + g_y_0_z_0_xy_xyzzzz[k];

                g_y_0_z_0_xxy_zzzzz[k] = -g_y_0_z_0_xy_zzzzz[k] * ab_x + g_y_0_z_0_xy_xzzzzz[k];
            }

            /// Set up 1092-1113 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 1092 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 1093 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 1094 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 1095 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 1096 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 1097 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 1098 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 1099 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 1100 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 1101 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 1102 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 1103 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 1104 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 1105 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 1106 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 1107 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 1108 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 1109 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 1110 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 1111 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 1112 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxz_xxxxx, g_y_0_z_0_xxz_xxxxy, g_y_0_z_0_xxz_xxxxz, g_y_0_z_0_xxz_xxxyy, g_y_0_z_0_xxz_xxxyz, g_y_0_z_0_xxz_xxxzz, g_y_0_z_0_xxz_xxyyy, g_y_0_z_0_xxz_xxyyz, g_y_0_z_0_xxz_xxyzz, g_y_0_z_0_xxz_xxzzz, g_y_0_z_0_xxz_xyyyy, g_y_0_z_0_xxz_xyyyz, g_y_0_z_0_xxz_xyyzz, g_y_0_z_0_xxz_xyzzz, g_y_0_z_0_xxz_xzzzz, g_y_0_z_0_xxz_yyyyy, g_y_0_z_0_xxz_yyyyz, g_y_0_z_0_xxz_yyyzz, g_y_0_z_0_xxz_yyzzz, g_y_0_z_0_xxz_yzzzz, g_y_0_z_0_xxz_zzzzz, g_y_0_z_0_xz_xxxxx, g_y_0_z_0_xz_xxxxxx, g_y_0_z_0_xz_xxxxxy, g_y_0_z_0_xz_xxxxxz, g_y_0_z_0_xz_xxxxy, g_y_0_z_0_xz_xxxxyy, g_y_0_z_0_xz_xxxxyz, g_y_0_z_0_xz_xxxxz, g_y_0_z_0_xz_xxxxzz, g_y_0_z_0_xz_xxxyy, g_y_0_z_0_xz_xxxyyy, g_y_0_z_0_xz_xxxyyz, g_y_0_z_0_xz_xxxyz, g_y_0_z_0_xz_xxxyzz, g_y_0_z_0_xz_xxxzz, g_y_0_z_0_xz_xxxzzz, g_y_0_z_0_xz_xxyyy, g_y_0_z_0_xz_xxyyyy, g_y_0_z_0_xz_xxyyyz, g_y_0_z_0_xz_xxyyz, g_y_0_z_0_xz_xxyyzz, g_y_0_z_0_xz_xxyzz, g_y_0_z_0_xz_xxyzzz, g_y_0_z_0_xz_xxzzz, g_y_0_z_0_xz_xxzzzz, g_y_0_z_0_xz_xyyyy, g_y_0_z_0_xz_xyyyyy, g_y_0_z_0_xz_xyyyyz, g_y_0_z_0_xz_xyyyz, g_y_0_z_0_xz_xyyyzz, g_y_0_z_0_xz_xyyzz, g_y_0_z_0_xz_xyyzzz, g_y_0_z_0_xz_xyzzz, g_y_0_z_0_xz_xyzzzz, g_y_0_z_0_xz_xzzzz, g_y_0_z_0_xz_xzzzzz, g_y_0_z_0_xz_yyyyy, g_y_0_z_0_xz_yyyyz, g_y_0_z_0_xz_yyyzz, g_y_0_z_0_xz_yyzzz, g_y_0_z_0_xz_yzzzz, g_y_0_z_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxz_xxxxx[k] = -g_y_0_z_0_xz_xxxxx[k] * ab_x + g_y_0_z_0_xz_xxxxxx[k];

                g_y_0_z_0_xxz_xxxxy[k] = -g_y_0_z_0_xz_xxxxy[k] * ab_x + g_y_0_z_0_xz_xxxxxy[k];

                g_y_0_z_0_xxz_xxxxz[k] = -g_y_0_z_0_xz_xxxxz[k] * ab_x + g_y_0_z_0_xz_xxxxxz[k];

                g_y_0_z_0_xxz_xxxyy[k] = -g_y_0_z_0_xz_xxxyy[k] * ab_x + g_y_0_z_0_xz_xxxxyy[k];

                g_y_0_z_0_xxz_xxxyz[k] = -g_y_0_z_0_xz_xxxyz[k] * ab_x + g_y_0_z_0_xz_xxxxyz[k];

                g_y_0_z_0_xxz_xxxzz[k] = -g_y_0_z_0_xz_xxxzz[k] * ab_x + g_y_0_z_0_xz_xxxxzz[k];

                g_y_0_z_0_xxz_xxyyy[k] = -g_y_0_z_0_xz_xxyyy[k] * ab_x + g_y_0_z_0_xz_xxxyyy[k];

                g_y_0_z_0_xxz_xxyyz[k] = -g_y_0_z_0_xz_xxyyz[k] * ab_x + g_y_0_z_0_xz_xxxyyz[k];

                g_y_0_z_0_xxz_xxyzz[k] = -g_y_0_z_0_xz_xxyzz[k] * ab_x + g_y_0_z_0_xz_xxxyzz[k];

                g_y_0_z_0_xxz_xxzzz[k] = -g_y_0_z_0_xz_xxzzz[k] * ab_x + g_y_0_z_0_xz_xxxzzz[k];

                g_y_0_z_0_xxz_xyyyy[k] = -g_y_0_z_0_xz_xyyyy[k] * ab_x + g_y_0_z_0_xz_xxyyyy[k];

                g_y_0_z_0_xxz_xyyyz[k] = -g_y_0_z_0_xz_xyyyz[k] * ab_x + g_y_0_z_0_xz_xxyyyz[k];

                g_y_0_z_0_xxz_xyyzz[k] = -g_y_0_z_0_xz_xyyzz[k] * ab_x + g_y_0_z_0_xz_xxyyzz[k];

                g_y_0_z_0_xxz_xyzzz[k] = -g_y_0_z_0_xz_xyzzz[k] * ab_x + g_y_0_z_0_xz_xxyzzz[k];

                g_y_0_z_0_xxz_xzzzz[k] = -g_y_0_z_0_xz_xzzzz[k] * ab_x + g_y_0_z_0_xz_xxzzzz[k];

                g_y_0_z_0_xxz_yyyyy[k] = -g_y_0_z_0_xz_yyyyy[k] * ab_x + g_y_0_z_0_xz_xyyyyy[k];

                g_y_0_z_0_xxz_yyyyz[k] = -g_y_0_z_0_xz_yyyyz[k] * ab_x + g_y_0_z_0_xz_xyyyyz[k];

                g_y_0_z_0_xxz_yyyzz[k] = -g_y_0_z_0_xz_yyyzz[k] * ab_x + g_y_0_z_0_xz_xyyyzz[k];

                g_y_0_z_0_xxz_yyzzz[k] = -g_y_0_z_0_xz_yyzzz[k] * ab_x + g_y_0_z_0_xz_xyyzzz[k];

                g_y_0_z_0_xxz_yzzzz[k] = -g_y_0_z_0_xz_yzzzz[k] * ab_x + g_y_0_z_0_xz_xyzzzz[k];

                g_y_0_z_0_xxz_zzzzz[k] = -g_y_0_z_0_xz_zzzzz[k] * ab_x + g_y_0_z_0_xz_xzzzzz[k];
            }

            /// Set up 1113-1134 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1113 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1114 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1115 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1116 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1117 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1118 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1119 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1120 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1121 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1122 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1123 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1124 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1125 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1126 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1127 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1128 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1129 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1130 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1131 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1132 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyy_xxxxx, g_y_0_z_0_xyy_xxxxy, g_y_0_z_0_xyy_xxxxz, g_y_0_z_0_xyy_xxxyy, g_y_0_z_0_xyy_xxxyz, g_y_0_z_0_xyy_xxxzz, g_y_0_z_0_xyy_xxyyy, g_y_0_z_0_xyy_xxyyz, g_y_0_z_0_xyy_xxyzz, g_y_0_z_0_xyy_xxzzz, g_y_0_z_0_xyy_xyyyy, g_y_0_z_0_xyy_xyyyz, g_y_0_z_0_xyy_xyyzz, g_y_0_z_0_xyy_xyzzz, g_y_0_z_0_xyy_xzzzz, g_y_0_z_0_xyy_yyyyy, g_y_0_z_0_xyy_yyyyz, g_y_0_z_0_xyy_yyyzz, g_y_0_z_0_xyy_yyzzz, g_y_0_z_0_xyy_yzzzz, g_y_0_z_0_xyy_zzzzz, g_y_0_z_0_yy_xxxxx, g_y_0_z_0_yy_xxxxxx, g_y_0_z_0_yy_xxxxxy, g_y_0_z_0_yy_xxxxxz, g_y_0_z_0_yy_xxxxy, g_y_0_z_0_yy_xxxxyy, g_y_0_z_0_yy_xxxxyz, g_y_0_z_0_yy_xxxxz, g_y_0_z_0_yy_xxxxzz, g_y_0_z_0_yy_xxxyy, g_y_0_z_0_yy_xxxyyy, g_y_0_z_0_yy_xxxyyz, g_y_0_z_0_yy_xxxyz, g_y_0_z_0_yy_xxxyzz, g_y_0_z_0_yy_xxxzz, g_y_0_z_0_yy_xxxzzz, g_y_0_z_0_yy_xxyyy, g_y_0_z_0_yy_xxyyyy, g_y_0_z_0_yy_xxyyyz, g_y_0_z_0_yy_xxyyz, g_y_0_z_0_yy_xxyyzz, g_y_0_z_0_yy_xxyzz, g_y_0_z_0_yy_xxyzzz, g_y_0_z_0_yy_xxzzz, g_y_0_z_0_yy_xxzzzz, g_y_0_z_0_yy_xyyyy, g_y_0_z_0_yy_xyyyyy, g_y_0_z_0_yy_xyyyyz, g_y_0_z_0_yy_xyyyz, g_y_0_z_0_yy_xyyyzz, g_y_0_z_0_yy_xyyzz, g_y_0_z_0_yy_xyyzzz, g_y_0_z_0_yy_xyzzz, g_y_0_z_0_yy_xyzzzz, g_y_0_z_0_yy_xzzzz, g_y_0_z_0_yy_xzzzzz, g_y_0_z_0_yy_yyyyy, g_y_0_z_0_yy_yyyyz, g_y_0_z_0_yy_yyyzz, g_y_0_z_0_yy_yyzzz, g_y_0_z_0_yy_yzzzz, g_y_0_z_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyy_xxxxx[k] = -g_y_0_z_0_yy_xxxxx[k] * ab_x + g_y_0_z_0_yy_xxxxxx[k];

                g_y_0_z_0_xyy_xxxxy[k] = -g_y_0_z_0_yy_xxxxy[k] * ab_x + g_y_0_z_0_yy_xxxxxy[k];

                g_y_0_z_0_xyy_xxxxz[k] = -g_y_0_z_0_yy_xxxxz[k] * ab_x + g_y_0_z_0_yy_xxxxxz[k];

                g_y_0_z_0_xyy_xxxyy[k] = -g_y_0_z_0_yy_xxxyy[k] * ab_x + g_y_0_z_0_yy_xxxxyy[k];

                g_y_0_z_0_xyy_xxxyz[k] = -g_y_0_z_0_yy_xxxyz[k] * ab_x + g_y_0_z_0_yy_xxxxyz[k];

                g_y_0_z_0_xyy_xxxzz[k] = -g_y_0_z_0_yy_xxxzz[k] * ab_x + g_y_0_z_0_yy_xxxxzz[k];

                g_y_0_z_0_xyy_xxyyy[k] = -g_y_0_z_0_yy_xxyyy[k] * ab_x + g_y_0_z_0_yy_xxxyyy[k];

                g_y_0_z_0_xyy_xxyyz[k] = -g_y_0_z_0_yy_xxyyz[k] * ab_x + g_y_0_z_0_yy_xxxyyz[k];

                g_y_0_z_0_xyy_xxyzz[k] = -g_y_0_z_0_yy_xxyzz[k] * ab_x + g_y_0_z_0_yy_xxxyzz[k];

                g_y_0_z_0_xyy_xxzzz[k] = -g_y_0_z_0_yy_xxzzz[k] * ab_x + g_y_0_z_0_yy_xxxzzz[k];

                g_y_0_z_0_xyy_xyyyy[k] = -g_y_0_z_0_yy_xyyyy[k] * ab_x + g_y_0_z_0_yy_xxyyyy[k];

                g_y_0_z_0_xyy_xyyyz[k] = -g_y_0_z_0_yy_xyyyz[k] * ab_x + g_y_0_z_0_yy_xxyyyz[k];

                g_y_0_z_0_xyy_xyyzz[k] = -g_y_0_z_0_yy_xyyzz[k] * ab_x + g_y_0_z_0_yy_xxyyzz[k];

                g_y_0_z_0_xyy_xyzzz[k] = -g_y_0_z_0_yy_xyzzz[k] * ab_x + g_y_0_z_0_yy_xxyzzz[k];

                g_y_0_z_0_xyy_xzzzz[k] = -g_y_0_z_0_yy_xzzzz[k] * ab_x + g_y_0_z_0_yy_xxzzzz[k];

                g_y_0_z_0_xyy_yyyyy[k] = -g_y_0_z_0_yy_yyyyy[k] * ab_x + g_y_0_z_0_yy_xyyyyy[k];

                g_y_0_z_0_xyy_yyyyz[k] = -g_y_0_z_0_yy_yyyyz[k] * ab_x + g_y_0_z_0_yy_xyyyyz[k];

                g_y_0_z_0_xyy_yyyzz[k] = -g_y_0_z_0_yy_yyyzz[k] * ab_x + g_y_0_z_0_yy_xyyyzz[k];

                g_y_0_z_0_xyy_yyzzz[k] = -g_y_0_z_0_yy_yyzzz[k] * ab_x + g_y_0_z_0_yy_xyyzzz[k];

                g_y_0_z_0_xyy_yzzzz[k] = -g_y_0_z_0_yy_yzzzz[k] * ab_x + g_y_0_z_0_yy_xyzzzz[k];

                g_y_0_z_0_xyy_zzzzz[k] = -g_y_0_z_0_yy_zzzzz[k] * ab_x + g_y_0_z_0_yy_xzzzzz[k];
            }

            /// Set up 1134-1155 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1134 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1135 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1136 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1137 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1138 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1139 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1140 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1141 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1142 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1143 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1144 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1145 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1146 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1147 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1148 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1149 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1150 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1151 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1152 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1153 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyz_xxxxx, g_y_0_z_0_xyz_xxxxy, g_y_0_z_0_xyz_xxxxz, g_y_0_z_0_xyz_xxxyy, g_y_0_z_0_xyz_xxxyz, g_y_0_z_0_xyz_xxxzz, g_y_0_z_0_xyz_xxyyy, g_y_0_z_0_xyz_xxyyz, g_y_0_z_0_xyz_xxyzz, g_y_0_z_0_xyz_xxzzz, g_y_0_z_0_xyz_xyyyy, g_y_0_z_0_xyz_xyyyz, g_y_0_z_0_xyz_xyyzz, g_y_0_z_0_xyz_xyzzz, g_y_0_z_0_xyz_xzzzz, g_y_0_z_0_xyz_yyyyy, g_y_0_z_0_xyz_yyyyz, g_y_0_z_0_xyz_yyyzz, g_y_0_z_0_xyz_yyzzz, g_y_0_z_0_xyz_yzzzz, g_y_0_z_0_xyz_zzzzz, g_y_0_z_0_yz_xxxxx, g_y_0_z_0_yz_xxxxxx, g_y_0_z_0_yz_xxxxxy, g_y_0_z_0_yz_xxxxxz, g_y_0_z_0_yz_xxxxy, g_y_0_z_0_yz_xxxxyy, g_y_0_z_0_yz_xxxxyz, g_y_0_z_0_yz_xxxxz, g_y_0_z_0_yz_xxxxzz, g_y_0_z_0_yz_xxxyy, g_y_0_z_0_yz_xxxyyy, g_y_0_z_0_yz_xxxyyz, g_y_0_z_0_yz_xxxyz, g_y_0_z_0_yz_xxxyzz, g_y_0_z_0_yz_xxxzz, g_y_0_z_0_yz_xxxzzz, g_y_0_z_0_yz_xxyyy, g_y_0_z_0_yz_xxyyyy, g_y_0_z_0_yz_xxyyyz, g_y_0_z_0_yz_xxyyz, g_y_0_z_0_yz_xxyyzz, g_y_0_z_0_yz_xxyzz, g_y_0_z_0_yz_xxyzzz, g_y_0_z_0_yz_xxzzz, g_y_0_z_0_yz_xxzzzz, g_y_0_z_0_yz_xyyyy, g_y_0_z_0_yz_xyyyyy, g_y_0_z_0_yz_xyyyyz, g_y_0_z_0_yz_xyyyz, g_y_0_z_0_yz_xyyyzz, g_y_0_z_0_yz_xyyzz, g_y_0_z_0_yz_xyyzzz, g_y_0_z_0_yz_xyzzz, g_y_0_z_0_yz_xyzzzz, g_y_0_z_0_yz_xzzzz, g_y_0_z_0_yz_xzzzzz, g_y_0_z_0_yz_yyyyy, g_y_0_z_0_yz_yyyyz, g_y_0_z_0_yz_yyyzz, g_y_0_z_0_yz_yyzzz, g_y_0_z_0_yz_yzzzz, g_y_0_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyz_xxxxx[k] = -g_y_0_z_0_yz_xxxxx[k] * ab_x + g_y_0_z_0_yz_xxxxxx[k];

                g_y_0_z_0_xyz_xxxxy[k] = -g_y_0_z_0_yz_xxxxy[k] * ab_x + g_y_0_z_0_yz_xxxxxy[k];

                g_y_0_z_0_xyz_xxxxz[k] = -g_y_0_z_0_yz_xxxxz[k] * ab_x + g_y_0_z_0_yz_xxxxxz[k];

                g_y_0_z_0_xyz_xxxyy[k] = -g_y_0_z_0_yz_xxxyy[k] * ab_x + g_y_0_z_0_yz_xxxxyy[k];

                g_y_0_z_0_xyz_xxxyz[k] = -g_y_0_z_0_yz_xxxyz[k] * ab_x + g_y_0_z_0_yz_xxxxyz[k];

                g_y_0_z_0_xyz_xxxzz[k] = -g_y_0_z_0_yz_xxxzz[k] * ab_x + g_y_0_z_0_yz_xxxxzz[k];

                g_y_0_z_0_xyz_xxyyy[k] = -g_y_0_z_0_yz_xxyyy[k] * ab_x + g_y_0_z_0_yz_xxxyyy[k];

                g_y_0_z_0_xyz_xxyyz[k] = -g_y_0_z_0_yz_xxyyz[k] * ab_x + g_y_0_z_0_yz_xxxyyz[k];

                g_y_0_z_0_xyz_xxyzz[k] = -g_y_0_z_0_yz_xxyzz[k] * ab_x + g_y_0_z_0_yz_xxxyzz[k];

                g_y_0_z_0_xyz_xxzzz[k] = -g_y_0_z_0_yz_xxzzz[k] * ab_x + g_y_0_z_0_yz_xxxzzz[k];

                g_y_0_z_0_xyz_xyyyy[k] = -g_y_0_z_0_yz_xyyyy[k] * ab_x + g_y_0_z_0_yz_xxyyyy[k];

                g_y_0_z_0_xyz_xyyyz[k] = -g_y_0_z_0_yz_xyyyz[k] * ab_x + g_y_0_z_0_yz_xxyyyz[k];

                g_y_0_z_0_xyz_xyyzz[k] = -g_y_0_z_0_yz_xyyzz[k] * ab_x + g_y_0_z_0_yz_xxyyzz[k];

                g_y_0_z_0_xyz_xyzzz[k] = -g_y_0_z_0_yz_xyzzz[k] * ab_x + g_y_0_z_0_yz_xxyzzz[k];

                g_y_0_z_0_xyz_xzzzz[k] = -g_y_0_z_0_yz_xzzzz[k] * ab_x + g_y_0_z_0_yz_xxzzzz[k];

                g_y_0_z_0_xyz_yyyyy[k] = -g_y_0_z_0_yz_yyyyy[k] * ab_x + g_y_0_z_0_yz_xyyyyy[k];

                g_y_0_z_0_xyz_yyyyz[k] = -g_y_0_z_0_yz_yyyyz[k] * ab_x + g_y_0_z_0_yz_xyyyyz[k];

                g_y_0_z_0_xyz_yyyzz[k] = -g_y_0_z_0_yz_yyyzz[k] * ab_x + g_y_0_z_0_yz_xyyyzz[k];

                g_y_0_z_0_xyz_yyzzz[k] = -g_y_0_z_0_yz_yyzzz[k] * ab_x + g_y_0_z_0_yz_xyyzzz[k];

                g_y_0_z_0_xyz_yzzzz[k] = -g_y_0_z_0_yz_yzzzz[k] * ab_x + g_y_0_z_0_yz_xyzzzz[k];

                g_y_0_z_0_xyz_zzzzz[k] = -g_y_0_z_0_yz_zzzzz[k] * ab_x + g_y_0_z_0_yz_xzzzzz[k];
            }

            /// Set up 1155-1176 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1155 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1156 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1157 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1158 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1159 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1160 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1161 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1162 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1163 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1164 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1165 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1166 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1167 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1168 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1169 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1170 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1171 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1172 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1173 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1174 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzz_xxxxx, g_y_0_z_0_xzz_xxxxy, g_y_0_z_0_xzz_xxxxz, g_y_0_z_0_xzz_xxxyy, g_y_0_z_0_xzz_xxxyz, g_y_0_z_0_xzz_xxxzz, g_y_0_z_0_xzz_xxyyy, g_y_0_z_0_xzz_xxyyz, g_y_0_z_0_xzz_xxyzz, g_y_0_z_0_xzz_xxzzz, g_y_0_z_0_xzz_xyyyy, g_y_0_z_0_xzz_xyyyz, g_y_0_z_0_xzz_xyyzz, g_y_0_z_0_xzz_xyzzz, g_y_0_z_0_xzz_xzzzz, g_y_0_z_0_xzz_yyyyy, g_y_0_z_0_xzz_yyyyz, g_y_0_z_0_xzz_yyyzz, g_y_0_z_0_xzz_yyzzz, g_y_0_z_0_xzz_yzzzz, g_y_0_z_0_xzz_zzzzz, g_y_0_z_0_zz_xxxxx, g_y_0_z_0_zz_xxxxxx, g_y_0_z_0_zz_xxxxxy, g_y_0_z_0_zz_xxxxxz, g_y_0_z_0_zz_xxxxy, g_y_0_z_0_zz_xxxxyy, g_y_0_z_0_zz_xxxxyz, g_y_0_z_0_zz_xxxxz, g_y_0_z_0_zz_xxxxzz, g_y_0_z_0_zz_xxxyy, g_y_0_z_0_zz_xxxyyy, g_y_0_z_0_zz_xxxyyz, g_y_0_z_0_zz_xxxyz, g_y_0_z_0_zz_xxxyzz, g_y_0_z_0_zz_xxxzz, g_y_0_z_0_zz_xxxzzz, g_y_0_z_0_zz_xxyyy, g_y_0_z_0_zz_xxyyyy, g_y_0_z_0_zz_xxyyyz, g_y_0_z_0_zz_xxyyz, g_y_0_z_0_zz_xxyyzz, g_y_0_z_0_zz_xxyzz, g_y_0_z_0_zz_xxyzzz, g_y_0_z_0_zz_xxzzz, g_y_0_z_0_zz_xxzzzz, g_y_0_z_0_zz_xyyyy, g_y_0_z_0_zz_xyyyyy, g_y_0_z_0_zz_xyyyyz, g_y_0_z_0_zz_xyyyz, g_y_0_z_0_zz_xyyyzz, g_y_0_z_0_zz_xyyzz, g_y_0_z_0_zz_xyyzzz, g_y_0_z_0_zz_xyzzz, g_y_0_z_0_zz_xyzzzz, g_y_0_z_0_zz_xzzzz, g_y_0_z_0_zz_xzzzzz, g_y_0_z_0_zz_yyyyy, g_y_0_z_0_zz_yyyyz, g_y_0_z_0_zz_yyyzz, g_y_0_z_0_zz_yyzzz, g_y_0_z_0_zz_yzzzz, g_y_0_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzz_xxxxx[k] = -g_y_0_z_0_zz_xxxxx[k] * ab_x + g_y_0_z_0_zz_xxxxxx[k];

                g_y_0_z_0_xzz_xxxxy[k] = -g_y_0_z_0_zz_xxxxy[k] * ab_x + g_y_0_z_0_zz_xxxxxy[k];

                g_y_0_z_0_xzz_xxxxz[k] = -g_y_0_z_0_zz_xxxxz[k] * ab_x + g_y_0_z_0_zz_xxxxxz[k];

                g_y_0_z_0_xzz_xxxyy[k] = -g_y_0_z_0_zz_xxxyy[k] * ab_x + g_y_0_z_0_zz_xxxxyy[k];

                g_y_0_z_0_xzz_xxxyz[k] = -g_y_0_z_0_zz_xxxyz[k] * ab_x + g_y_0_z_0_zz_xxxxyz[k];

                g_y_0_z_0_xzz_xxxzz[k] = -g_y_0_z_0_zz_xxxzz[k] * ab_x + g_y_0_z_0_zz_xxxxzz[k];

                g_y_0_z_0_xzz_xxyyy[k] = -g_y_0_z_0_zz_xxyyy[k] * ab_x + g_y_0_z_0_zz_xxxyyy[k];

                g_y_0_z_0_xzz_xxyyz[k] = -g_y_0_z_0_zz_xxyyz[k] * ab_x + g_y_0_z_0_zz_xxxyyz[k];

                g_y_0_z_0_xzz_xxyzz[k] = -g_y_0_z_0_zz_xxyzz[k] * ab_x + g_y_0_z_0_zz_xxxyzz[k];

                g_y_0_z_0_xzz_xxzzz[k] = -g_y_0_z_0_zz_xxzzz[k] * ab_x + g_y_0_z_0_zz_xxxzzz[k];

                g_y_0_z_0_xzz_xyyyy[k] = -g_y_0_z_0_zz_xyyyy[k] * ab_x + g_y_0_z_0_zz_xxyyyy[k];

                g_y_0_z_0_xzz_xyyyz[k] = -g_y_0_z_0_zz_xyyyz[k] * ab_x + g_y_0_z_0_zz_xxyyyz[k];

                g_y_0_z_0_xzz_xyyzz[k] = -g_y_0_z_0_zz_xyyzz[k] * ab_x + g_y_0_z_0_zz_xxyyzz[k];

                g_y_0_z_0_xzz_xyzzz[k] = -g_y_0_z_0_zz_xyzzz[k] * ab_x + g_y_0_z_0_zz_xxyzzz[k];

                g_y_0_z_0_xzz_xzzzz[k] = -g_y_0_z_0_zz_xzzzz[k] * ab_x + g_y_0_z_0_zz_xxzzzz[k];

                g_y_0_z_0_xzz_yyyyy[k] = -g_y_0_z_0_zz_yyyyy[k] * ab_x + g_y_0_z_0_zz_xyyyyy[k];

                g_y_0_z_0_xzz_yyyyz[k] = -g_y_0_z_0_zz_yyyyz[k] * ab_x + g_y_0_z_0_zz_xyyyyz[k];

                g_y_0_z_0_xzz_yyyzz[k] = -g_y_0_z_0_zz_yyyzz[k] * ab_x + g_y_0_z_0_zz_xyyyzz[k];

                g_y_0_z_0_xzz_yyzzz[k] = -g_y_0_z_0_zz_yyzzz[k] * ab_x + g_y_0_z_0_zz_xyyzzz[k];

                g_y_0_z_0_xzz_yzzzz[k] = -g_y_0_z_0_zz_yzzzz[k] * ab_x + g_y_0_z_0_zz_xyzzzz[k];

                g_y_0_z_0_xzz_zzzzz[k] = -g_y_0_z_0_zz_zzzzz[k] * ab_x + g_y_0_z_0_zz_xzzzzz[k];
            }

            /// Set up 1176-1197 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1176 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1177 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1178 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1179 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1180 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1181 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1182 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1183 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1184 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1185 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1186 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1187 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1188 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1189 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1190 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1191 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1192 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1193 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1194 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1195 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1196 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yy_xxxxx, g_0_0_z_0_yy_xxxxy, g_0_0_z_0_yy_xxxxz, g_0_0_z_0_yy_xxxyy, g_0_0_z_0_yy_xxxyz, g_0_0_z_0_yy_xxxzz, g_0_0_z_0_yy_xxyyy, g_0_0_z_0_yy_xxyyz, g_0_0_z_0_yy_xxyzz, g_0_0_z_0_yy_xxzzz, g_0_0_z_0_yy_xyyyy, g_0_0_z_0_yy_xyyyz, g_0_0_z_0_yy_xyyzz, g_0_0_z_0_yy_xyzzz, g_0_0_z_0_yy_xzzzz, g_0_0_z_0_yy_yyyyy, g_0_0_z_0_yy_yyyyz, g_0_0_z_0_yy_yyyzz, g_0_0_z_0_yy_yyzzz, g_0_0_z_0_yy_yzzzz, g_0_0_z_0_yy_zzzzz, g_y_0_z_0_yy_xxxxx, g_y_0_z_0_yy_xxxxxy, g_y_0_z_0_yy_xxxxy, g_y_0_z_0_yy_xxxxyy, g_y_0_z_0_yy_xxxxyz, g_y_0_z_0_yy_xxxxz, g_y_0_z_0_yy_xxxyy, g_y_0_z_0_yy_xxxyyy, g_y_0_z_0_yy_xxxyyz, g_y_0_z_0_yy_xxxyz, g_y_0_z_0_yy_xxxyzz, g_y_0_z_0_yy_xxxzz, g_y_0_z_0_yy_xxyyy, g_y_0_z_0_yy_xxyyyy, g_y_0_z_0_yy_xxyyyz, g_y_0_z_0_yy_xxyyz, g_y_0_z_0_yy_xxyyzz, g_y_0_z_0_yy_xxyzz, g_y_0_z_0_yy_xxyzzz, g_y_0_z_0_yy_xxzzz, g_y_0_z_0_yy_xyyyy, g_y_0_z_0_yy_xyyyyy, g_y_0_z_0_yy_xyyyyz, g_y_0_z_0_yy_xyyyz, g_y_0_z_0_yy_xyyyzz, g_y_0_z_0_yy_xyyzz, g_y_0_z_0_yy_xyyzzz, g_y_0_z_0_yy_xyzzz, g_y_0_z_0_yy_xyzzzz, g_y_0_z_0_yy_xzzzz, g_y_0_z_0_yy_yyyyy, g_y_0_z_0_yy_yyyyyy, g_y_0_z_0_yy_yyyyyz, g_y_0_z_0_yy_yyyyz, g_y_0_z_0_yy_yyyyzz, g_y_0_z_0_yy_yyyzz, g_y_0_z_0_yy_yyyzzz, g_y_0_z_0_yy_yyzzz, g_y_0_z_0_yy_yyzzzz, g_y_0_z_0_yy_yzzzz, g_y_0_z_0_yy_yzzzzz, g_y_0_z_0_yy_zzzzz, g_y_0_z_0_yyy_xxxxx, g_y_0_z_0_yyy_xxxxy, g_y_0_z_0_yyy_xxxxz, g_y_0_z_0_yyy_xxxyy, g_y_0_z_0_yyy_xxxyz, g_y_0_z_0_yyy_xxxzz, g_y_0_z_0_yyy_xxyyy, g_y_0_z_0_yyy_xxyyz, g_y_0_z_0_yyy_xxyzz, g_y_0_z_0_yyy_xxzzz, g_y_0_z_0_yyy_xyyyy, g_y_0_z_0_yyy_xyyyz, g_y_0_z_0_yyy_xyyzz, g_y_0_z_0_yyy_xyzzz, g_y_0_z_0_yyy_xzzzz, g_y_0_z_0_yyy_yyyyy, g_y_0_z_0_yyy_yyyyz, g_y_0_z_0_yyy_yyyzz, g_y_0_z_0_yyy_yyzzz, g_y_0_z_0_yyy_yzzzz, g_y_0_z_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyy_xxxxx[k] = -g_0_0_z_0_yy_xxxxx[k] - g_y_0_z_0_yy_xxxxx[k] * ab_y + g_y_0_z_0_yy_xxxxxy[k];

                g_y_0_z_0_yyy_xxxxy[k] = -g_0_0_z_0_yy_xxxxy[k] - g_y_0_z_0_yy_xxxxy[k] * ab_y + g_y_0_z_0_yy_xxxxyy[k];

                g_y_0_z_0_yyy_xxxxz[k] = -g_0_0_z_0_yy_xxxxz[k] - g_y_0_z_0_yy_xxxxz[k] * ab_y + g_y_0_z_0_yy_xxxxyz[k];

                g_y_0_z_0_yyy_xxxyy[k] = -g_0_0_z_0_yy_xxxyy[k] - g_y_0_z_0_yy_xxxyy[k] * ab_y + g_y_0_z_0_yy_xxxyyy[k];

                g_y_0_z_0_yyy_xxxyz[k] = -g_0_0_z_0_yy_xxxyz[k] - g_y_0_z_0_yy_xxxyz[k] * ab_y + g_y_0_z_0_yy_xxxyyz[k];

                g_y_0_z_0_yyy_xxxzz[k] = -g_0_0_z_0_yy_xxxzz[k] - g_y_0_z_0_yy_xxxzz[k] * ab_y + g_y_0_z_0_yy_xxxyzz[k];

                g_y_0_z_0_yyy_xxyyy[k] = -g_0_0_z_0_yy_xxyyy[k] - g_y_0_z_0_yy_xxyyy[k] * ab_y + g_y_0_z_0_yy_xxyyyy[k];

                g_y_0_z_0_yyy_xxyyz[k] = -g_0_0_z_0_yy_xxyyz[k] - g_y_0_z_0_yy_xxyyz[k] * ab_y + g_y_0_z_0_yy_xxyyyz[k];

                g_y_0_z_0_yyy_xxyzz[k] = -g_0_0_z_0_yy_xxyzz[k] - g_y_0_z_0_yy_xxyzz[k] * ab_y + g_y_0_z_0_yy_xxyyzz[k];

                g_y_0_z_0_yyy_xxzzz[k] = -g_0_0_z_0_yy_xxzzz[k] - g_y_0_z_0_yy_xxzzz[k] * ab_y + g_y_0_z_0_yy_xxyzzz[k];

                g_y_0_z_0_yyy_xyyyy[k] = -g_0_0_z_0_yy_xyyyy[k] - g_y_0_z_0_yy_xyyyy[k] * ab_y + g_y_0_z_0_yy_xyyyyy[k];

                g_y_0_z_0_yyy_xyyyz[k] = -g_0_0_z_0_yy_xyyyz[k] - g_y_0_z_0_yy_xyyyz[k] * ab_y + g_y_0_z_0_yy_xyyyyz[k];

                g_y_0_z_0_yyy_xyyzz[k] = -g_0_0_z_0_yy_xyyzz[k] - g_y_0_z_0_yy_xyyzz[k] * ab_y + g_y_0_z_0_yy_xyyyzz[k];

                g_y_0_z_0_yyy_xyzzz[k] = -g_0_0_z_0_yy_xyzzz[k] - g_y_0_z_0_yy_xyzzz[k] * ab_y + g_y_0_z_0_yy_xyyzzz[k];

                g_y_0_z_0_yyy_xzzzz[k] = -g_0_0_z_0_yy_xzzzz[k] - g_y_0_z_0_yy_xzzzz[k] * ab_y + g_y_0_z_0_yy_xyzzzz[k];

                g_y_0_z_0_yyy_yyyyy[k] = -g_0_0_z_0_yy_yyyyy[k] - g_y_0_z_0_yy_yyyyy[k] * ab_y + g_y_0_z_0_yy_yyyyyy[k];

                g_y_0_z_0_yyy_yyyyz[k] = -g_0_0_z_0_yy_yyyyz[k] - g_y_0_z_0_yy_yyyyz[k] * ab_y + g_y_0_z_0_yy_yyyyyz[k];

                g_y_0_z_0_yyy_yyyzz[k] = -g_0_0_z_0_yy_yyyzz[k] - g_y_0_z_0_yy_yyyzz[k] * ab_y + g_y_0_z_0_yy_yyyyzz[k];

                g_y_0_z_0_yyy_yyzzz[k] = -g_0_0_z_0_yy_yyzzz[k] - g_y_0_z_0_yy_yyzzz[k] * ab_y + g_y_0_z_0_yy_yyyzzz[k];

                g_y_0_z_0_yyy_yzzzz[k] = -g_0_0_z_0_yy_yzzzz[k] - g_y_0_z_0_yy_yzzzz[k] * ab_y + g_y_0_z_0_yy_yyzzzz[k];

                g_y_0_z_0_yyy_zzzzz[k] = -g_0_0_z_0_yy_zzzzz[k] - g_y_0_z_0_yy_zzzzz[k] * ab_y + g_y_0_z_0_yy_yzzzzz[k];
            }

            /// Set up 1197-1218 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1197 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1198 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1199 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1200 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1201 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1202 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1203 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1204 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1205 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1206 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1207 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1208 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1209 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1210 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1211 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1212 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1213 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1214 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1215 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1216 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1217 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yy_xxxxx, g_y_0_z_0_yy_xxxxxz, g_y_0_z_0_yy_xxxxy, g_y_0_z_0_yy_xxxxyz, g_y_0_z_0_yy_xxxxz, g_y_0_z_0_yy_xxxxzz, g_y_0_z_0_yy_xxxyy, g_y_0_z_0_yy_xxxyyz, g_y_0_z_0_yy_xxxyz, g_y_0_z_0_yy_xxxyzz, g_y_0_z_0_yy_xxxzz, g_y_0_z_0_yy_xxxzzz, g_y_0_z_0_yy_xxyyy, g_y_0_z_0_yy_xxyyyz, g_y_0_z_0_yy_xxyyz, g_y_0_z_0_yy_xxyyzz, g_y_0_z_0_yy_xxyzz, g_y_0_z_0_yy_xxyzzz, g_y_0_z_0_yy_xxzzz, g_y_0_z_0_yy_xxzzzz, g_y_0_z_0_yy_xyyyy, g_y_0_z_0_yy_xyyyyz, g_y_0_z_0_yy_xyyyz, g_y_0_z_0_yy_xyyyzz, g_y_0_z_0_yy_xyyzz, g_y_0_z_0_yy_xyyzzz, g_y_0_z_0_yy_xyzzz, g_y_0_z_0_yy_xyzzzz, g_y_0_z_0_yy_xzzzz, g_y_0_z_0_yy_xzzzzz, g_y_0_z_0_yy_yyyyy, g_y_0_z_0_yy_yyyyyz, g_y_0_z_0_yy_yyyyz, g_y_0_z_0_yy_yyyyzz, g_y_0_z_0_yy_yyyzz, g_y_0_z_0_yy_yyyzzz, g_y_0_z_0_yy_yyzzz, g_y_0_z_0_yy_yyzzzz, g_y_0_z_0_yy_yzzzz, g_y_0_z_0_yy_yzzzzz, g_y_0_z_0_yy_zzzzz, g_y_0_z_0_yy_zzzzzz, g_y_0_z_0_yyz_xxxxx, g_y_0_z_0_yyz_xxxxy, g_y_0_z_0_yyz_xxxxz, g_y_0_z_0_yyz_xxxyy, g_y_0_z_0_yyz_xxxyz, g_y_0_z_0_yyz_xxxzz, g_y_0_z_0_yyz_xxyyy, g_y_0_z_0_yyz_xxyyz, g_y_0_z_0_yyz_xxyzz, g_y_0_z_0_yyz_xxzzz, g_y_0_z_0_yyz_xyyyy, g_y_0_z_0_yyz_xyyyz, g_y_0_z_0_yyz_xyyzz, g_y_0_z_0_yyz_xyzzz, g_y_0_z_0_yyz_xzzzz, g_y_0_z_0_yyz_yyyyy, g_y_0_z_0_yyz_yyyyz, g_y_0_z_0_yyz_yyyzz, g_y_0_z_0_yyz_yyzzz, g_y_0_z_0_yyz_yzzzz, g_y_0_z_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyz_xxxxx[k] = -g_y_0_z_0_yy_xxxxx[k] * ab_z + g_y_0_z_0_yy_xxxxxz[k];

                g_y_0_z_0_yyz_xxxxy[k] = -g_y_0_z_0_yy_xxxxy[k] * ab_z + g_y_0_z_0_yy_xxxxyz[k];

                g_y_0_z_0_yyz_xxxxz[k] = -g_y_0_z_0_yy_xxxxz[k] * ab_z + g_y_0_z_0_yy_xxxxzz[k];

                g_y_0_z_0_yyz_xxxyy[k] = -g_y_0_z_0_yy_xxxyy[k] * ab_z + g_y_0_z_0_yy_xxxyyz[k];

                g_y_0_z_0_yyz_xxxyz[k] = -g_y_0_z_0_yy_xxxyz[k] * ab_z + g_y_0_z_0_yy_xxxyzz[k];

                g_y_0_z_0_yyz_xxxzz[k] = -g_y_0_z_0_yy_xxxzz[k] * ab_z + g_y_0_z_0_yy_xxxzzz[k];

                g_y_0_z_0_yyz_xxyyy[k] = -g_y_0_z_0_yy_xxyyy[k] * ab_z + g_y_0_z_0_yy_xxyyyz[k];

                g_y_0_z_0_yyz_xxyyz[k] = -g_y_0_z_0_yy_xxyyz[k] * ab_z + g_y_0_z_0_yy_xxyyzz[k];

                g_y_0_z_0_yyz_xxyzz[k] = -g_y_0_z_0_yy_xxyzz[k] * ab_z + g_y_0_z_0_yy_xxyzzz[k];

                g_y_0_z_0_yyz_xxzzz[k] = -g_y_0_z_0_yy_xxzzz[k] * ab_z + g_y_0_z_0_yy_xxzzzz[k];

                g_y_0_z_0_yyz_xyyyy[k] = -g_y_0_z_0_yy_xyyyy[k] * ab_z + g_y_0_z_0_yy_xyyyyz[k];

                g_y_0_z_0_yyz_xyyyz[k] = -g_y_0_z_0_yy_xyyyz[k] * ab_z + g_y_0_z_0_yy_xyyyzz[k];

                g_y_0_z_0_yyz_xyyzz[k] = -g_y_0_z_0_yy_xyyzz[k] * ab_z + g_y_0_z_0_yy_xyyzzz[k];

                g_y_0_z_0_yyz_xyzzz[k] = -g_y_0_z_0_yy_xyzzz[k] * ab_z + g_y_0_z_0_yy_xyzzzz[k];

                g_y_0_z_0_yyz_xzzzz[k] = -g_y_0_z_0_yy_xzzzz[k] * ab_z + g_y_0_z_0_yy_xzzzzz[k];

                g_y_0_z_0_yyz_yyyyy[k] = -g_y_0_z_0_yy_yyyyy[k] * ab_z + g_y_0_z_0_yy_yyyyyz[k];

                g_y_0_z_0_yyz_yyyyz[k] = -g_y_0_z_0_yy_yyyyz[k] * ab_z + g_y_0_z_0_yy_yyyyzz[k];

                g_y_0_z_0_yyz_yyyzz[k] = -g_y_0_z_0_yy_yyyzz[k] * ab_z + g_y_0_z_0_yy_yyyzzz[k];

                g_y_0_z_0_yyz_yyzzz[k] = -g_y_0_z_0_yy_yyzzz[k] * ab_z + g_y_0_z_0_yy_yyzzzz[k];

                g_y_0_z_0_yyz_yzzzz[k] = -g_y_0_z_0_yy_yzzzz[k] * ab_z + g_y_0_z_0_yy_yzzzzz[k];

                g_y_0_z_0_yyz_zzzzz[k] = -g_y_0_z_0_yy_zzzzz[k] * ab_z + g_y_0_z_0_yy_zzzzzz[k];
            }

            /// Set up 1218-1239 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1218 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1219 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1220 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1221 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1222 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1223 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1224 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1225 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1226 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1227 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1228 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1229 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1230 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1231 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1232 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1233 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1234 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1235 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1236 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1237 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1238 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yz_xxxxx, g_y_0_z_0_yz_xxxxxz, g_y_0_z_0_yz_xxxxy, g_y_0_z_0_yz_xxxxyz, g_y_0_z_0_yz_xxxxz, g_y_0_z_0_yz_xxxxzz, g_y_0_z_0_yz_xxxyy, g_y_0_z_0_yz_xxxyyz, g_y_0_z_0_yz_xxxyz, g_y_0_z_0_yz_xxxyzz, g_y_0_z_0_yz_xxxzz, g_y_0_z_0_yz_xxxzzz, g_y_0_z_0_yz_xxyyy, g_y_0_z_0_yz_xxyyyz, g_y_0_z_0_yz_xxyyz, g_y_0_z_0_yz_xxyyzz, g_y_0_z_0_yz_xxyzz, g_y_0_z_0_yz_xxyzzz, g_y_0_z_0_yz_xxzzz, g_y_0_z_0_yz_xxzzzz, g_y_0_z_0_yz_xyyyy, g_y_0_z_0_yz_xyyyyz, g_y_0_z_0_yz_xyyyz, g_y_0_z_0_yz_xyyyzz, g_y_0_z_0_yz_xyyzz, g_y_0_z_0_yz_xyyzzz, g_y_0_z_0_yz_xyzzz, g_y_0_z_0_yz_xyzzzz, g_y_0_z_0_yz_xzzzz, g_y_0_z_0_yz_xzzzzz, g_y_0_z_0_yz_yyyyy, g_y_0_z_0_yz_yyyyyz, g_y_0_z_0_yz_yyyyz, g_y_0_z_0_yz_yyyyzz, g_y_0_z_0_yz_yyyzz, g_y_0_z_0_yz_yyyzzz, g_y_0_z_0_yz_yyzzz, g_y_0_z_0_yz_yyzzzz, g_y_0_z_0_yz_yzzzz, g_y_0_z_0_yz_yzzzzz, g_y_0_z_0_yz_zzzzz, g_y_0_z_0_yz_zzzzzz, g_y_0_z_0_yzz_xxxxx, g_y_0_z_0_yzz_xxxxy, g_y_0_z_0_yzz_xxxxz, g_y_0_z_0_yzz_xxxyy, g_y_0_z_0_yzz_xxxyz, g_y_0_z_0_yzz_xxxzz, g_y_0_z_0_yzz_xxyyy, g_y_0_z_0_yzz_xxyyz, g_y_0_z_0_yzz_xxyzz, g_y_0_z_0_yzz_xxzzz, g_y_0_z_0_yzz_xyyyy, g_y_0_z_0_yzz_xyyyz, g_y_0_z_0_yzz_xyyzz, g_y_0_z_0_yzz_xyzzz, g_y_0_z_0_yzz_xzzzz, g_y_0_z_0_yzz_yyyyy, g_y_0_z_0_yzz_yyyyz, g_y_0_z_0_yzz_yyyzz, g_y_0_z_0_yzz_yyzzz, g_y_0_z_0_yzz_yzzzz, g_y_0_z_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzz_xxxxx[k] = -g_y_0_z_0_yz_xxxxx[k] * ab_z + g_y_0_z_0_yz_xxxxxz[k];

                g_y_0_z_0_yzz_xxxxy[k] = -g_y_0_z_0_yz_xxxxy[k] * ab_z + g_y_0_z_0_yz_xxxxyz[k];

                g_y_0_z_0_yzz_xxxxz[k] = -g_y_0_z_0_yz_xxxxz[k] * ab_z + g_y_0_z_0_yz_xxxxzz[k];

                g_y_0_z_0_yzz_xxxyy[k] = -g_y_0_z_0_yz_xxxyy[k] * ab_z + g_y_0_z_0_yz_xxxyyz[k];

                g_y_0_z_0_yzz_xxxyz[k] = -g_y_0_z_0_yz_xxxyz[k] * ab_z + g_y_0_z_0_yz_xxxyzz[k];

                g_y_0_z_0_yzz_xxxzz[k] = -g_y_0_z_0_yz_xxxzz[k] * ab_z + g_y_0_z_0_yz_xxxzzz[k];

                g_y_0_z_0_yzz_xxyyy[k] = -g_y_0_z_0_yz_xxyyy[k] * ab_z + g_y_0_z_0_yz_xxyyyz[k];

                g_y_0_z_0_yzz_xxyyz[k] = -g_y_0_z_0_yz_xxyyz[k] * ab_z + g_y_0_z_0_yz_xxyyzz[k];

                g_y_0_z_0_yzz_xxyzz[k] = -g_y_0_z_0_yz_xxyzz[k] * ab_z + g_y_0_z_0_yz_xxyzzz[k];

                g_y_0_z_0_yzz_xxzzz[k] = -g_y_0_z_0_yz_xxzzz[k] * ab_z + g_y_0_z_0_yz_xxzzzz[k];

                g_y_0_z_0_yzz_xyyyy[k] = -g_y_0_z_0_yz_xyyyy[k] * ab_z + g_y_0_z_0_yz_xyyyyz[k];

                g_y_0_z_0_yzz_xyyyz[k] = -g_y_0_z_0_yz_xyyyz[k] * ab_z + g_y_0_z_0_yz_xyyyzz[k];

                g_y_0_z_0_yzz_xyyzz[k] = -g_y_0_z_0_yz_xyyzz[k] * ab_z + g_y_0_z_0_yz_xyyzzz[k];

                g_y_0_z_0_yzz_xyzzz[k] = -g_y_0_z_0_yz_xyzzz[k] * ab_z + g_y_0_z_0_yz_xyzzzz[k];

                g_y_0_z_0_yzz_xzzzz[k] = -g_y_0_z_0_yz_xzzzz[k] * ab_z + g_y_0_z_0_yz_xzzzzz[k];

                g_y_0_z_0_yzz_yyyyy[k] = -g_y_0_z_0_yz_yyyyy[k] * ab_z + g_y_0_z_0_yz_yyyyyz[k];

                g_y_0_z_0_yzz_yyyyz[k] = -g_y_0_z_0_yz_yyyyz[k] * ab_z + g_y_0_z_0_yz_yyyyzz[k];

                g_y_0_z_0_yzz_yyyzz[k] = -g_y_0_z_0_yz_yyyzz[k] * ab_z + g_y_0_z_0_yz_yyyzzz[k];

                g_y_0_z_0_yzz_yyzzz[k] = -g_y_0_z_0_yz_yyzzz[k] * ab_z + g_y_0_z_0_yz_yyzzzz[k];

                g_y_0_z_0_yzz_yzzzz[k] = -g_y_0_z_0_yz_yzzzz[k] * ab_z + g_y_0_z_0_yz_yzzzzz[k];

                g_y_0_z_0_yzz_zzzzz[k] = -g_y_0_z_0_yz_zzzzz[k] * ab_z + g_y_0_z_0_yz_zzzzzz[k];
            }

            /// Set up 1239-1260 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1239 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1240 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1241 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1242 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1243 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1244 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1245 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1246 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1247 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1248 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1249 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1250 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1251 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1252 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1253 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1254 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1255 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1256 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1257 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1258 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zz_xxxxx, g_y_0_z_0_zz_xxxxxz, g_y_0_z_0_zz_xxxxy, g_y_0_z_0_zz_xxxxyz, g_y_0_z_0_zz_xxxxz, g_y_0_z_0_zz_xxxxzz, g_y_0_z_0_zz_xxxyy, g_y_0_z_0_zz_xxxyyz, g_y_0_z_0_zz_xxxyz, g_y_0_z_0_zz_xxxyzz, g_y_0_z_0_zz_xxxzz, g_y_0_z_0_zz_xxxzzz, g_y_0_z_0_zz_xxyyy, g_y_0_z_0_zz_xxyyyz, g_y_0_z_0_zz_xxyyz, g_y_0_z_0_zz_xxyyzz, g_y_0_z_0_zz_xxyzz, g_y_0_z_0_zz_xxyzzz, g_y_0_z_0_zz_xxzzz, g_y_0_z_0_zz_xxzzzz, g_y_0_z_0_zz_xyyyy, g_y_0_z_0_zz_xyyyyz, g_y_0_z_0_zz_xyyyz, g_y_0_z_0_zz_xyyyzz, g_y_0_z_0_zz_xyyzz, g_y_0_z_0_zz_xyyzzz, g_y_0_z_0_zz_xyzzz, g_y_0_z_0_zz_xyzzzz, g_y_0_z_0_zz_xzzzz, g_y_0_z_0_zz_xzzzzz, g_y_0_z_0_zz_yyyyy, g_y_0_z_0_zz_yyyyyz, g_y_0_z_0_zz_yyyyz, g_y_0_z_0_zz_yyyyzz, g_y_0_z_0_zz_yyyzz, g_y_0_z_0_zz_yyyzzz, g_y_0_z_0_zz_yyzzz, g_y_0_z_0_zz_yyzzzz, g_y_0_z_0_zz_yzzzz, g_y_0_z_0_zz_yzzzzz, g_y_0_z_0_zz_zzzzz, g_y_0_z_0_zz_zzzzzz, g_y_0_z_0_zzz_xxxxx, g_y_0_z_0_zzz_xxxxy, g_y_0_z_0_zzz_xxxxz, g_y_0_z_0_zzz_xxxyy, g_y_0_z_0_zzz_xxxyz, g_y_0_z_0_zzz_xxxzz, g_y_0_z_0_zzz_xxyyy, g_y_0_z_0_zzz_xxyyz, g_y_0_z_0_zzz_xxyzz, g_y_0_z_0_zzz_xxzzz, g_y_0_z_0_zzz_xyyyy, g_y_0_z_0_zzz_xyyyz, g_y_0_z_0_zzz_xyyzz, g_y_0_z_0_zzz_xyzzz, g_y_0_z_0_zzz_xzzzz, g_y_0_z_0_zzz_yyyyy, g_y_0_z_0_zzz_yyyyz, g_y_0_z_0_zzz_yyyzz, g_y_0_z_0_zzz_yyzzz, g_y_0_z_0_zzz_yzzzz, g_y_0_z_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzz_xxxxx[k] = -g_y_0_z_0_zz_xxxxx[k] * ab_z + g_y_0_z_0_zz_xxxxxz[k];

                g_y_0_z_0_zzz_xxxxy[k] = -g_y_0_z_0_zz_xxxxy[k] * ab_z + g_y_0_z_0_zz_xxxxyz[k];

                g_y_0_z_0_zzz_xxxxz[k] = -g_y_0_z_0_zz_xxxxz[k] * ab_z + g_y_0_z_0_zz_xxxxzz[k];

                g_y_0_z_0_zzz_xxxyy[k] = -g_y_0_z_0_zz_xxxyy[k] * ab_z + g_y_0_z_0_zz_xxxyyz[k];

                g_y_0_z_0_zzz_xxxyz[k] = -g_y_0_z_0_zz_xxxyz[k] * ab_z + g_y_0_z_0_zz_xxxyzz[k];

                g_y_0_z_0_zzz_xxxzz[k] = -g_y_0_z_0_zz_xxxzz[k] * ab_z + g_y_0_z_0_zz_xxxzzz[k];

                g_y_0_z_0_zzz_xxyyy[k] = -g_y_0_z_0_zz_xxyyy[k] * ab_z + g_y_0_z_0_zz_xxyyyz[k];

                g_y_0_z_0_zzz_xxyyz[k] = -g_y_0_z_0_zz_xxyyz[k] * ab_z + g_y_0_z_0_zz_xxyyzz[k];

                g_y_0_z_0_zzz_xxyzz[k] = -g_y_0_z_0_zz_xxyzz[k] * ab_z + g_y_0_z_0_zz_xxyzzz[k];

                g_y_0_z_0_zzz_xxzzz[k] = -g_y_0_z_0_zz_xxzzz[k] * ab_z + g_y_0_z_0_zz_xxzzzz[k];

                g_y_0_z_0_zzz_xyyyy[k] = -g_y_0_z_0_zz_xyyyy[k] * ab_z + g_y_0_z_0_zz_xyyyyz[k];

                g_y_0_z_0_zzz_xyyyz[k] = -g_y_0_z_0_zz_xyyyz[k] * ab_z + g_y_0_z_0_zz_xyyyzz[k];

                g_y_0_z_0_zzz_xyyzz[k] = -g_y_0_z_0_zz_xyyzz[k] * ab_z + g_y_0_z_0_zz_xyyzzz[k];

                g_y_0_z_0_zzz_xyzzz[k] = -g_y_0_z_0_zz_xyzzz[k] * ab_z + g_y_0_z_0_zz_xyzzzz[k];

                g_y_0_z_0_zzz_xzzzz[k] = -g_y_0_z_0_zz_xzzzz[k] * ab_z + g_y_0_z_0_zz_xzzzzz[k];

                g_y_0_z_0_zzz_yyyyy[k] = -g_y_0_z_0_zz_yyyyy[k] * ab_z + g_y_0_z_0_zz_yyyyyz[k];

                g_y_0_z_0_zzz_yyyyz[k] = -g_y_0_z_0_zz_yyyyz[k] * ab_z + g_y_0_z_0_zz_yyyyzz[k];

                g_y_0_z_0_zzz_yyyzz[k] = -g_y_0_z_0_zz_yyyzz[k] * ab_z + g_y_0_z_0_zz_yyyzzz[k];

                g_y_0_z_0_zzz_yyzzz[k] = -g_y_0_z_0_zz_yyzzz[k] * ab_z + g_y_0_z_0_zz_yyzzzz[k];

                g_y_0_z_0_zzz_yzzzz[k] = -g_y_0_z_0_zz_yzzzz[k] * ab_z + g_y_0_z_0_zz_yzzzzz[k];

                g_y_0_z_0_zzz_zzzzz[k] = -g_y_0_z_0_zz_zzzzz[k] * ab_z + g_y_0_z_0_zz_zzzzzz[k];
            }

            /// Set up 1260-1281 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 1260 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 1261 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 1262 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 1263 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 1264 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 1265 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 1266 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 1267 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 1268 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 1269 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 1270 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 1271 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 1272 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 1273 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 1274 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 1275 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 1276 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 1277 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 1278 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 1279 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 1280 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xx_xxxxx, g_z_0_x_0_xx_xxxxxx, g_z_0_x_0_xx_xxxxxy, g_z_0_x_0_xx_xxxxxz, g_z_0_x_0_xx_xxxxy, g_z_0_x_0_xx_xxxxyy, g_z_0_x_0_xx_xxxxyz, g_z_0_x_0_xx_xxxxz, g_z_0_x_0_xx_xxxxzz, g_z_0_x_0_xx_xxxyy, g_z_0_x_0_xx_xxxyyy, g_z_0_x_0_xx_xxxyyz, g_z_0_x_0_xx_xxxyz, g_z_0_x_0_xx_xxxyzz, g_z_0_x_0_xx_xxxzz, g_z_0_x_0_xx_xxxzzz, g_z_0_x_0_xx_xxyyy, g_z_0_x_0_xx_xxyyyy, g_z_0_x_0_xx_xxyyyz, g_z_0_x_0_xx_xxyyz, g_z_0_x_0_xx_xxyyzz, g_z_0_x_0_xx_xxyzz, g_z_0_x_0_xx_xxyzzz, g_z_0_x_0_xx_xxzzz, g_z_0_x_0_xx_xxzzzz, g_z_0_x_0_xx_xyyyy, g_z_0_x_0_xx_xyyyyy, g_z_0_x_0_xx_xyyyyz, g_z_0_x_0_xx_xyyyz, g_z_0_x_0_xx_xyyyzz, g_z_0_x_0_xx_xyyzz, g_z_0_x_0_xx_xyyzzz, g_z_0_x_0_xx_xyzzz, g_z_0_x_0_xx_xyzzzz, g_z_0_x_0_xx_xzzzz, g_z_0_x_0_xx_xzzzzz, g_z_0_x_0_xx_yyyyy, g_z_0_x_0_xx_yyyyz, g_z_0_x_0_xx_yyyzz, g_z_0_x_0_xx_yyzzz, g_z_0_x_0_xx_yzzzz, g_z_0_x_0_xx_zzzzz, g_z_0_x_0_xxx_xxxxx, g_z_0_x_0_xxx_xxxxy, g_z_0_x_0_xxx_xxxxz, g_z_0_x_0_xxx_xxxyy, g_z_0_x_0_xxx_xxxyz, g_z_0_x_0_xxx_xxxzz, g_z_0_x_0_xxx_xxyyy, g_z_0_x_0_xxx_xxyyz, g_z_0_x_0_xxx_xxyzz, g_z_0_x_0_xxx_xxzzz, g_z_0_x_0_xxx_xyyyy, g_z_0_x_0_xxx_xyyyz, g_z_0_x_0_xxx_xyyzz, g_z_0_x_0_xxx_xyzzz, g_z_0_x_0_xxx_xzzzz, g_z_0_x_0_xxx_yyyyy, g_z_0_x_0_xxx_yyyyz, g_z_0_x_0_xxx_yyyzz, g_z_0_x_0_xxx_yyzzz, g_z_0_x_0_xxx_yzzzz, g_z_0_x_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxx_xxxxx[k] = -g_z_0_x_0_xx_xxxxx[k] * ab_x + g_z_0_x_0_xx_xxxxxx[k];

                g_z_0_x_0_xxx_xxxxy[k] = -g_z_0_x_0_xx_xxxxy[k] * ab_x + g_z_0_x_0_xx_xxxxxy[k];

                g_z_0_x_0_xxx_xxxxz[k] = -g_z_0_x_0_xx_xxxxz[k] * ab_x + g_z_0_x_0_xx_xxxxxz[k];

                g_z_0_x_0_xxx_xxxyy[k] = -g_z_0_x_0_xx_xxxyy[k] * ab_x + g_z_0_x_0_xx_xxxxyy[k];

                g_z_0_x_0_xxx_xxxyz[k] = -g_z_0_x_0_xx_xxxyz[k] * ab_x + g_z_0_x_0_xx_xxxxyz[k];

                g_z_0_x_0_xxx_xxxzz[k] = -g_z_0_x_0_xx_xxxzz[k] * ab_x + g_z_0_x_0_xx_xxxxzz[k];

                g_z_0_x_0_xxx_xxyyy[k] = -g_z_0_x_0_xx_xxyyy[k] * ab_x + g_z_0_x_0_xx_xxxyyy[k];

                g_z_0_x_0_xxx_xxyyz[k] = -g_z_0_x_0_xx_xxyyz[k] * ab_x + g_z_0_x_0_xx_xxxyyz[k];

                g_z_0_x_0_xxx_xxyzz[k] = -g_z_0_x_0_xx_xxyzz[k] * ab_x + g_z_0_x_0_xx_xxxyzz[k];

                g_z_0_x_0_xxx_xxzzz[k] = -g_z_0_x_0_xx_xxzzz[k] * ab_x + g_z_0_x_0_xx_xxxzzz[k];

                g_z_0_x_0_xxx_xyyyy[k] = -g_z_0_x_0_xx_xyyyy[k] * ab_x + g_z_0_x_0_xx_xxyyyy[k];

                g_z_0_x_0_xxx_xyyyz[k] = -g_z_0_x_0_xx_xyyyz[k] * ab_x + g_z_0_x_0_xx_xxyyyz[k];

                g_z_0_x_0_xxx_xyyzz[k] = -g_z_0_x_0_xx_xyyzz[k] * ab_x + g_z_0_x_0_xx_xxyyzz[k];

                g_z_0_x_0_xxx_xyzzz[k] = -g_z_0_x_0_xx_xyzzz[k] * ab_x + g_z_0_x_0_xx_xxyzzz[k];

                g_z_0_x_0_xxx_xzzzz[k] = -g_z_0_x_0_xx_xzzzz[k] * ab_x + g_z_0_x_0_xx_xxzzzz[k];

                g_z_0_x_0_xxx_yyyyy[k] = -g_z_0_x_0_xx_yyyyy[k] * ab_x + g_z_0_x_0_xx_xyyyyy[k];

                g_z_0_x_0_xxx_yyyyz[k] = -g_z_0_x_0_xx_yyyyz[k] * ab_x + g_z_0_x_0_xx_xyyyyz[k];

                g_z_0_x_0_xxx_yyyzz[k] = -g_z_0_x_0_xx_yyyzz[k] * ab_x + g_z_0_x_0_xx_xyyyzz[k];

                g_z_0_x_0_xxx_yyzzz[k] = -g_z_0_x_0_xx_yyzzz[k] * ab_x + g_z_0_x_0_xx_xyyzzz[k];

                g_z_0_x_0_xxx_yzzzz[k] = -g_z_0_x_0_xx_yzzzz[k] * ab_x + g_z_0_x_0_xx_xyzzzz[k];

                g_z_0_x_0_xxx_zzzzz[k] = -g_z_0_x_0_xx_zzzzz[k] * ab_x + g_z_0_x_0_xx_xzzzzz[k];
            }

            /// Set up 1281-1302 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 1281 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 1282 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 1283 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 1284 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 1285 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 1286 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 1287 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 1288 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 1289 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 1290 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 1291 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 1292 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 1293 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 1294 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 1295 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 1296 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 1297 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 1298 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 1299 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 1300 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 1301 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxy_xxxxx, g_z_0_x_0_xxy_xxxxy, g_z_0_x_0_xxy_xxxxz, g_z_0_x_0_xxy_xxxyy, g_z_0_x_0_xxy_xxxyz, g_z_0_x_0_xxy_xxxzz, g_z_0_x_0_xxy_xxyyy, g_z_0_x_0_xxy_xxyyz, g_z_0_x_0_xxy_xxyzz, g_z_0_x_0_xxy_xxzzz, g_z_0_x_0_xxy_xyyyy, g_z_0_x_0_xxy_xyyyz, g_z_0_x_0_xxy_xyyzz, g_z_0_x_0_xxy_xyzzz, g_z_0_x_0_xxy_xzzzz, g_z_0_x_0_xxy_yyyyy, g_z_0_x_0_xxy_yyyyz, g_z_0_x_0_xxy_yyyzz, g_z_0_x_0_xxy_yyzzz, g_z_0_x_0_xxy_yzzzz, g_z_0_x_0_xxy_zzzzz, g_z_0_x_0_xy_xxxxx, g_z_0_x_0_xy_xxxxxx, g_z_0_x_0_xy_xxxxxy, g_z_0_x_0_xy_xxxxxz, g_z_0_x_0_xy_xxxxy, g_z_0_x_0_xy_xxxxyy, g_z_0_x_0_xy_xxxxyz, g_z_0_x_0_xy_xxxxz, g_z_0_x_0_xy_xxxxzz, g_z_0_x_0_xy_xxxyy, g_z_0_x_0_xy_xxxyyy, g_z_0_x_0_xy_xxxyyz, g_z_0_x_0_xy_xxxyz, g_z_0_x_0_xy_xxxyzz, g_z_0_x_0_xy_xxxzz, g_z_0_x_0_xy_xxxzzz, g_z_0_x_0_xy_xxyyy, g_z_0_x_0_xy_xxyyyy, g_z_0_x_0_xy_xxyyyz, g_z_0_x_0_xy_xxyyz, g_z_0_x_0_xy_xxyyzz, g_z_0_x_0_xy_xxyzz, g_z_0_x_0_xy_xxyzzz, g_z_0_x_0_xy_xxzzz, g_z_0_x_0_xy_xxzzzz, g_z_0_x_0_xy_xyyyy, g_z_0_x_0_xy_xyyyyy, g_z_0_x_0_xy_xyyyyz, g_z_0_x_0_xy_xyyyz, g_z_0_x_0_xy_xyyyzz, g_z_0_x_0_xy_xyyzz, g_z_0_x_0_xy_xyyzzz, g_z_0_x_0_xy_xyzzz, g_z_0_x_0_xy_xyzzzz, g_z_0_x_0_xy_xzzzz, g_z_0_x_0_xy_xzzzzz, g_z_0_x_0_xy_yyyyy, g_z_0_x_0_xy_yyyyz, g_z_0_x_0_xy_yyyzz, g_z_0_x_0_xy_yyzzz, g_z_0_x_0_xy_yzzzz, g_z_0_x_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxy_xxxxx[k] = -g_z_0_x_0_xy_xxxxx[k] * ab_x + g_z_0_x_0_xy_xxxxxx[k];

                g_z_0_x_0_xxy_xxxxy[k] = -g_z_0_x_0_xy_xxxxy[k] * ab_x + g_z_0_x_0_xy_xxxxxy[k];

                g_z_0_x_0_xxy_xxxxz[k] = -g_z_0_x_0_xy_xxxxz[k] * ab_x + g_z_0_x_0_xy_xxxxxz[k];

                g_z_0_x_0_xxy_xxxyy[k] = -g_z_0_x_0_xy_xxxyy[k] * ab_x + g_z_0_x_0_xy_xxxxyy[k];

                g_z_0_x_0_xxy_xxxyz[k] = -g_z_0_x_0_xy_xxxyz[k] * ab_x + g_z_0_x_0_xy_xxxxyz[k];

                g_z_0_x_0_xxy_xxxzz[k] = -g_z_0_x_0_xy_xxxzz[k] * ab_x + g_z_0_x_0_xy_xxxxzz[k];

                g_z_0_x_0_xxy_xxyyy[k] = -g_z_0_x_0_xy_xxyyy[k] * ab_x + g_z_0_x_0_xy_xxxyyy[k];

                g_z_0_x_0_xxy_xxyyz[k] = -g_z_0_x_0_xy_xxyyz[k] * ab_x + g_z_0_x_0_xy_xxxyyz[k];

                g_z_0_x_0_xxy_xxyzz[k] = -g_z_0_x_0_xy_xxyzz[k] * ab_x + g_z_0_x_0_xy_xxxyzz[k];

                g_z_0_x_0_xxy_xxzzz[k] = -g_z_0_x_0_xy_xxzzz[k] * ab_x + g_z_0_x_0_xy_xxxzzz[k];

                g_z_0_x_0_xxy_xyyyy[k] = -g_z_0_x_0_xy_xyyyy[k] * ab_x + g_z_0_x_0_xy_xxyyyy[k];

                g_z_0_x_0_xxy_xyyyz[k] = -g_z_0_x_0_xy_xyyyz[k] * ab_x + g_z_0_x_0_xy_xxyyyz[k];

                g_z_0_x_0_xxy_xyyzz[k] = -g_z_0_x_0_xy_xyyzz[k] * ab_x + g_z_0_x_0_xy_xxyyzz[k];

                g_z_0_x_0_xxy_xyzzz[k] = -g_z_0_x_0_xy_xyzzz[k] * ab_x + g_z_0_x_0_xy_xxyzzz[k];

                g_z_0_x_0_xxy_xzzzz[k] = -g_z_0_x_0_xy_xzzzz[k] * ab_x + g_z_0_x_0_xy_xxzzzz[k];

                g_z_0_x_0_xxy_yyyyy[k] = -g_z_0_x_0_xy_yyyyy[k] * ab_x + g_z_0_x_0_xy_xyyyyy[k];

                g_z_0_x_0_xxy_yyyyz[k] = -g_z_0_x_0_xy_yyyyz[k] * ab_x + g_z_0_x_0_xy_xyyyyz[k];

                g_z_0_x_0_xxy_yyyzz[k] = -g_z_0_x_0_xy_yyyzz[k] * ab_x + g_z_0_x_0_xy_xyyyzz[k];

                g_z_0_x_0_xxy_yyzzz[k] = -g_z_0_x_0_xy_yyzzz[k] * ab_x + g_z_0_x_0_xy_xyyzzz[k];

                g_z_0_x_0_xxy_yzzzz[k] = -g_z_0_x_0_xy_yzzzz[k] * ab_x + g_z_0_x_0_xy_xyzzzz[k];

                g_z_0_x_0_xxy_zzzzz[k] = -g_z_0_x_0_xy_zzzzz[k] * ab_x + g_z_0_x_0_xy_xzzzzz[k];
            }

            /// Set up 1302-1323 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 1302 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 1303 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 1304 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 1305 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 1306 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 1307 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 1308 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 1309 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 1310 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 1311 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 1312 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 1313 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 1314 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 1315 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 1316 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 1317 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 1318 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 1319 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 1320 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 1321 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 1322 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxz_xxxxx, g_z_0_x_0_xxz_xxxxy, g_z_0_x_0_xxz_xxxxz, g_z_0_x_0_xxz_xxxyy, g_z_0_x_0_xxz_xxxyz, g_z_0_x_0_xxz_xxxzz, g_z_0_x_0_xxz_xxyyy, g_z_0_x_0_xxz_xxyyz, g_z_0_x_0_xxz_xxyzz, g_z_0_x_0_xxz_xxzzz, g_z_0_x_0_xxz_xyyyy, g_z_0_x_0_xxz_xyyyz, g_z_0_x_0_xxz_xyyzz, g_z_0_x_0_xxz_xyzzz, g_z_0_x_0_xxz_xzzzz, g_z_0_x_0_xxz_yyyyy, g_z_0_x_0_xxz_yyyyz, g_z_0_x_0_xxz_yyyzz, g_z_0_x_0_xxz_yyzzz, g_z_0_x_0_xxz_yzzzz, g_z_0_x_0_xxz_zzzzz, g_z_0_x_0_xz_xxxxx, g_z_0_x_0_xz_xxxxxx, g_z_0_x_0_xz_xxxxxy, g_z_0_x_0_xz_xxxxxz, g_z_0_x_0_xz_xxxxy, g_z_0_x_0_xz_xxxxyy, g_z_0_x_0_xz_xxxxyz, g_z_0_x_0_xz_xxxxz, g_z_0_x_0_xz_xxxxzz, g_z_0_x_0_xz_xxxyy, g_z_0_x_0_xz_xxxyyy, g_z_0_x_0_xz_xxxyyz, g_z_0_x_0_xz_xxxyz, g_z_0_x_0_xz_xxxyzz, g_z_0_x_0_xz_xxxzz, g_z_0_x_0_xz_xxxzzz, g_z_0_x_0_xz_xxyyy, g_z_0_x_0_xz_xxyyyy, g_z_0_x_0_xz_xxyyyz, g_z_0_x_0_xz_xxyyz, g_z_0_x_0_xz_xxyyzz, g_z_0_x_0_xz_xxyzz, g_z_0_x_0_xz_xxyzzz, g_z_0_x_0_xz_xxzzz, g_z_0_x_0_xz_xxzzzz, g_z_0_x_0_xz_xyyyy, g_z_0_x_0_xz_xyyyyy, g_z_0_x_0_xz_xyyyyz, g_z_0_x_0_xz_xyyyz, g_z_0_x_0_xz_xyyyzz, g_z_0_x_0_xz_xyyzz, g_z_0_x_0_xz_xyyzzz, g_z_0_x_0_xz_xyzzz, g_z_0_x_0_xz_xyzzzz, g_z_0_x_0_xz_xzzzz, g_z_0_x_0_xz_xzzzzz, g_z_0_x_0_xz_yyyyy, g_z_0_x_0_xz_yyyyz, g_z_0_x_0_xz_yyyzz, g_z_0_x_0_xz_yyzzz, g_z_0_x_0_xz_yzzzz, g_z_0_x_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxz_xxxxx[k] = -g_z_0_x_0_xz_xxxxx[k] * ab_x + g_z_0_x_0_xz_xxxxxx[k];

                g_z_0_x_0_xxz_xxxxy[k] = -g_z_0_x_0_xz_xxxxy[k] * ab_x + g_z_0_x_0_xz_xxxxxy[k];

                g_z_0_x_0_xxz_xxxxz[k] = -g_z_0_x_0_xz_xxxxz[k] * ab_x + g_z_0_x_0_xz_xxxxxz[k];

                g_z_0_x_0_xxz_xxxyy[k] = -g_z_0_x_0_xz_xxxyy[k] * ab_x + g_z_0_x_0_xz_xxxxyy[k];

                g_z_0_x_0_xxz_xxxyz[k] = -g_z_0_x_0_xz_xxxyz[k] * ab_x + g_z_0_x_0_xz_xxxxyz[k];

                g_z_0_x_0_xxz_xxxzz[k] = -g_z_0_x_0_xz_xxxzz[k] * ab_x + g_z_0_x_0_xz_xxxxzz[k];

                g_z_0_x_0_xxz_xxyyy[k] = -g_z_0_x_0_xz_xxyyy[k] * ab_x + g_z_0_x_0_xz_xxxyyy[k];

                g_z_0_x_0_xxz_xxyyz[k] = -g_z_0_x_0_xz_xxyyz[k] * ab_x + g_z_0_x_0_xz_xxxyyz[k];

                g_z_0_x_0_xxz_xxyzz[k] = -g_z_0_x_0_xz_xxyzz[k] * ab_x + g_z_0_x_0_xz_xxxyzz[k];

                g_z_0_x_0_xxz_xxzzz[k] = -g_z_0_x_0_xz_xxzzz[k] * ab_x + g_z_0_x_0_xz_xxxzzz[k];

                g_z_0_x_0_xxz_xyyyy[k] = -g_z_0_x_0_xz_xyyyy[k] * ab_x + g_z_0_x_0_xz_xxyyyy[k];

                g_z_0_x_0_xxz_xyyyz[k] = -g_z_0_x_0_xz_xyyyz[k] * ab_x + g_z_0_x_0_xz_xxyyyz[k];

                g_z_0_x_0_xxz_xyyzz[k] = -g_z_0_x_0_xz_xyyzz[k] * ab_x + g_z_0_x_0_xz_xxyyzz[k];

                g_z_0_x_0_xxz_xyzzz[k] = -g_z_0_x_0_xz_xyzzz[k] * ab_x + g_z_0_x_0_xz_xxyzzz[k];

                g_z_0_x_0_xxz_xzzzz[k] = -g_z_0_x_0_xz_xzzzz[k] * ab_x + g_z_0_x_0_xz_xxzzzz[k];

                g_z_0_x_0_xxz_yyyyy[k] = -g_z_0_x_0_xz_yyyyy[k] * ab_x + g_z_0_x_0_xz_xyyyyy[k];

                g_z_0_x_0_xxz_yyyyz[k] = -g_z_0_x_0_xz_yyyyz[k] * ab_x + g_z_0_x_0_xz_xyyyyz[k];

                g_z_0_x_0_xxz_yyyzz[k] = -g_z_0_x_0_xz_yyyzz[k] * ab_x + g_z_0_x_0_xz_xyyyzz[k];

                g_z_0_x_0_xxz_yyzzz[k] = -g_z_0_x_0_xz_yyzzz[k] * ab_x + g_z_0_x_0_xz_xyyzzz[k];

                g_z_0_x_0_xxz_yzzzz[k] = -g_z_0_x_0_xz_yzzzz[k] * ab_x + g_z_0_x_0_xz_xyzzzz[k];

                g_z_0_x_0_xxz_zzzzz[k] = -g_z_0_x_0_xz_zzzzz[k] * ab_x + g_z_0_x_0_xz_xzzzzz[k];
            }

            /// Set up 1323-1344 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1323 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1324 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1325 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1326 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1327 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1328 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1329 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1330 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1331 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1332 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1333 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1334 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1335 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1336 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1337 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1338 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1339 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1340 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1341 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1342 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyy_xxxxx, g_z_0_x_0_xyy_xxxxy, g_z_0_x_0_xyy_xxxxz, g_z_0_x_0_xyy_xxxyy, g_z_0_x_0_xyy_xxxyz, g_z_0_x_0_xyy_xxxzz, g_z_0_x_0_xyy_xxyyy, g_z_0_x_0_xyy_xxyyz, g_z_0_x_0_xyy_xxyzz, g_z_0_x_0_xyy_xxzzz, g_z_0_x_0_xyy_xyyyy, g_z_0_x_0_xyy_xyyyz, g_z_0_x_0_xyy_xyyzz, g_z_0_x_0_xyy_xyzzz, g_z_0_x_0_xyy_xzzzz, g_z_0_x_0_xyy_yyyyy, g_z_0_x_0_xyy_yyyyz, g_z_0_x_0_xyy_yyyzz, g_z_0_x_0_xyy_yyzzz, g_z_0_x_0_xyy_yzzzz, g_z_0_x_0_xyy_zzzzz, g_z_0_x_0_yy_xxxxx, g_z_0_x_0_yy_xxxxxx, g_z_0_x_0_yy_xxxxxy, g_z_0_x_0_yy_xxxxxz, g_z_0_x_0_yy_xxxxy, g_z_0_x_0_yy_xxxxyy, g_z_0_x_0_yy_xxxxyz, g_z_0_x_0_yy_xxxxz, g_z_0_x_0_yy_xxxxzz, g_z_0_x_0_yy_xxxyy, g_z_0_x_0_yy_xxxyyy, g_z_0_x_0_yy_xxxyyz, g_z_0_x_0_yy_xxxyz, g_z_0_x_0_yy_xxxyzz, g_z_0_x_0_yy_xxxzz, g_z_0_x_0_yy_xxxzzz, g_z_0_x_0_yy_xxyyy, g_z_0_x_0_yy_xxyyyy, g_z_0_x_0_yy_xxyyyz, g_z_0_x_0_yy_xxyyz, g_z_0_x_0_yy_xxyyzz, g_z_0_x_0_yy_xxyzz, g_z_0_x_0_yy_xxyzzz, g_z_0_x_0_yy_xxzzz, g_z_0_x_0_yy_xxzzzz, g_z_0_x_0_yy_xyyyy, g_z_0_x_0_yy_xyyyyy, g_z_0_x_0_yy_xyyyyz, g_z_0_x_0_yy_xyyyz, g_z_0_x_0_yy_xyyyzz, g_z_0_x_0_yy_xyyzz, g_z_0_x_0_yy_xyyzzz, g_z_0_x_0_yy_xyzzz, g_z_0_x_0_yy_xyzzzz, g_z_0_x_0_yy_xzzzz, g_z_0_x_0_yy_xzzzzz, g_z_0_x_0_yy_yyyyy, g_z_0_x_0_yy_yyyyz, g_z_0_x_0_yy_yyyzz, g_z_0_x_0_yy_yyzzz, g_z_0_x_0_yy_yzzzz, g_z_0_x_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyy_xxxxx[k] = -g_z_0_x_0_yy_xxxxx[k] * ab_x + g_z_0_x_0_yy_xxxxxx[k];

                g_z_0_x_0_xyy_xxxxy[k] = -g_z_0_x_0_yy_xxxxy[k] * ab_x + g_z_0_x_0_yy_xxxxxy[k];

                g_z_0_x_0_xyy_xxxxz[k] = -g_z_0_x_0_yy_xxxxz[k] * ab_x + g_z_0_x_0_yy_xxxxxz[k];

                g_z_0_x_0_xyy_xxxyy[k] = -g_z_0_x_0_yy_xxxyy[k] * ab_x + g_z_0_x_0_yy_xxxxyy[k];

                g_z_0_x_0_xyy_xxxyz[k] = -g_z_0_x_0_yy_xxxyz[k] * ab_x + g_z_0_x_0_yy_xxxxyz[k];

                g_z_0_x_0_xyy_xxxzz[k] = -g_z_0_x_0_yy_xxxzz[k] * ab_x + g_z_0_x_0_yy_xxxxzz[k];

                g_z_0_x_0_xyy_xxyyy[k] = -g_z_0_x_0_yy_xxyyy[k] * ab_x + g_z_0_x_0_yy_xxxyyy[k];

                g_z_0_x_0_xyy_xxyyz[k] = -g_z_0_x_0_yy_xxyyz[k] * ab_x + g_z_0_x_0_yy_xxxyyz[k];

                g_z_0_x_0_xyy_xxyzz[k] = -g_z_0_x_0_yy_xxyzz[k] * ab_x + g_z_0_x_0_yy_xxxyzz[k];

                g_z_0_x_0_xyy_xxzzz[k] = -g_z_0_x_0_yy_xxzzz[k] * ab_x + g_z_0_x_0_yy_xxxzzz[k];

                g_z_0_x_0_xyy_xyyyy[k] = -g_z_0_x_0_yy_xyyyy[k] * ab_x + g_z_0_x_0_yy_xxyyyy[k];

                g_z_0_x_0_xyy_xyyyz[k] = -g_z_0_x_0_yy_xyyyz[k] * ab_x + g_z_0_x_0_yy_xxyyyz[k];

                g_z_0_x_0_xyy_xyyzz[k] = -g_z_0_x_0_yy_xyyzz[k] * ab_x + g_z_0_x_0_yy_xxyyzz[k];

                g_z_0_x_0_xyy_xyzzz[k] = -g_z_0_x_0_yy_xyzzz[k] * ab_x + g_z_0_x_0_yy_xxyzzz[k];

                g_z_0_x_0_xyy_xzzzz[k] = -g_z_0_x_0_yy_xzzzz[k] * ab_x + g_z_0_x_0_yy_xxzzzz[k];

                g_z_0_x_0_xyy_yyyyy[k] = -g_z_0_x_0_yy_yyyyy[k] * ab_x + g_z_0_x_0_yy_xyyyyy[k];

                g_z_0_x_0_xyy_yyyyz[k] = -g_z_0_x_0_yy_yyyyz[k] * ab_x + g_z_0_x_0_yy_xyyyyz[k];

                g_z_0_x_0_xyy_yyyzz[k] = -g_z_0_x_0_yy_yyyzz[k] * ab_x + g_z_0_x_0_yy_xyyyzz[k];

                g_z_0_x_0_xyy_yyzzz[k] = -g_z_0_x_0_yy_yyzzz[k] * ab_x + g_z_0_x_0_yy_xyyzzz[k];

                g_z_0_x_0_xyy_yzzzz[k] = -g_z_0_x_0_yy_yzzzz[k] * ab_x + g_z_0_x_0_yy_xyzzzz[k];

                g_z_0_x_0_xyy_zzzzz[k] = -g_z_0_x_0_yy_zzzzz[k] * ab_x + g_z_0_x_0_yy_xzzzzz[k];
            }

            /// Set up 1344-1365 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1344 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1345 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1346 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1347 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1348 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1349 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1350 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1351 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1352 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1353 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1354 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1355 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1356 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1357 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1358 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1359 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1360 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1361 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1362 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1363 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1364 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyz_xxxxx, g_z_0_x_0_xyz_xxxxy, g_z_0_x_0_xyz_xxxxz, g_z_0_x_0_xyz_xxxyy, g_z_0_x_0_xyz_xxxyz, g_z_0_x_0_xyz_xxxzz, g_z_0_x_0_xyz_xxyyy, g_z_0_x_0_xyz_xxyyz, g_z_0_x_0_xyz_xxyzz, g_z_0_x_0_xyz_xxzzz, g_z_0_x_0_xyz_xyyyy, g_z_0_x_0_xyz_xyyyz, g_z_0_x_0_xyz_xyyzz, g_z_0_x_0_xyz_xyzzz, g_z_0_x_0_xyz_xzzzz, g_z_0_x_0_xyz_yyyyy, g_z_0_x_0_xyz_yyyyz, g_z_0_x_0_xyz_yyyzz, g_z_0_x_0_xyz_yyzzz, g_z_0_x_0_xyz_yzzzz, g_z_0_x_0_xyz_zzzzz, g_z_0_x_0_yz_xxxxx, g_z_0_x_0_yz_xxxxxx, g_z_0_x_0_yz_xxxxxy, g_z_0_x_0_yz_xxxxxz, g_z_0_x_0_yz_xxxxy, g_z_0_x_0_yz_xxxxyy, g_z_0_x_0_yz_xxxxyz, g_z_0_x_0_yz_xxxxz, g_z_0_x_0_yz_xxxxzz, g_z_0_x_0_yz_xxxyy, g_z_0_x_0_yz_xxxyyy, g_z_0_x_0_yz_xxxyyz, g_z_0_x_0_yz_xxxyz, g_z_0_x_0_yz_xxxyzz, g_z_0_x_0_yz_xxxzz, g_z_0_x_0_yz_xxxzzz, g_z_0_x_0_yz_xxyyy, g_z_0_x_0_yz_xxyyyy, g_z_0_x_0_yz_xxyyyz, g_z_0_x_0_yz_xxyyz, g_z_0_x_0_yz_xxyyzz, g_z_0_x_0_yz_xxyzz, g_z_0_x_0_yz_xxyzzz, g_z_0_x_0_yz_xxzzz, g_z_0_x_0_yz_xxzzzz, g_z_0_x_0_yz_xyyyy, g_z_0_x_0_yz_xyyyyy, g_z_0_x_0_yz_xyyyyz, g_z_0_x_0_yz_xyyyz, g_z_0_x_0_yz_xyyyzz, g_z_0_x_0_yz_xyyzz, g_z_0_x_0_yz_xyyzzz, g_z_0_x_0_yz_xyzzz, g_z_0_x_0_yz_xyzzzz, g_z_0_x_0_yz_xzzzz, g_z_0_x_0_yz_xzzzzz, g_z_0_x_0_yz_yyyyy, g_z_0_x_0_yz_yyyyz, g_z_0_x_0_yz_yyyzz, g_z_0_x_0_yz_yyzzz, g_z_0_x_0_yz_yzzzz, g_z_0_x_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyz_xxxxx[k] = -g_z_0_x_0_yz_xxxxx[k] * ab_x + g_z_0_x_0_yz_xxxxxx[k];

                g_z_0_x_0_xyz_xxxxy[k] = -g_z_0_x_0_yz_xxxxy[k] * ab_x + g_z_0_x_0_yz_xxxxxy[k];

                g_z_0_x_0_xyz_xxxxz[k] = -g_z_0_x_0_yz_xxxxz[k] * ab_x + g_z_0_x_0_yz_xxxxxz[k];

                g_z_0_x_0_xyz_xxxyy[k] = -g_z_0_x_0_yz_xxxyy[k] * ab_x + g_z_0_x_0_yz_xxxxyy[k];

                g_z_0_x_0_xyz_xxxyz[k] = -g_z_0_x_0_yz_xxxyz[k] * ab_x + g_z_0_x_0_yz_xxxxyz[k];

                g_z_0_x_0_xyz_xxxzz[k] = -g_z_0_x_0_yz_xxxzz[k] * ab_x + g_z_0_x_0_yz_xxxxzz[k];

                g_z_0_x_0_xyz_xxyyy[k] = -g_z_0_x_0_yz_xxyyy[k] * ab_x + g_z_0_x_0_yz_xxxyyy[k];

                g_z_0_x_0_xyz_xxyyz[k] = -g_z_0_x_0_yz_xxyyz[k] * ab_x + g_z_0_x_0_yz_xxxyyz[k];

                g_z_0_x_0_xyz_xxyzz[k] = -g_z_0_x_0_yz_xxyzz[k] * ab_x + g_z_0_x_0_yz_xxxyzz[k];

                g_z_0_x_0_xyz_xxzzz[k] = -g_z_0_x_0_yz_xxzzz[k] * ab_x + g_z_0_x_0_yz_xxxzzz[k];

                g_z_0_x_0_xyz_xyyyy[k] = -g_z_0_x_0_yz_xyyyy[k] * ab_x + g_z_0_x_0_yz_xxyyyy[k];

                g_z_0_x_0_xyz_xyyyz[k] = -g_z_0_x_0_yz_xyyyz[k] * ab_x + g_z_0_x_0_yz_xxyyyz[k];

                g_z_0_x_0_xyz_xyyzz[k] = -g_z_0_x_0_yz_xyyzz[k] * ab_x + g_z_0_x_0_yz_xxyyzz[k];

                g_z_0_x_0_xyz_xyzzz[k] = -g_z_0_x_0_yz_xyzzz[k] * ab_x + g_z_0_x_0_yz_xxyzzz[k];

                g_z_0_x_0_xyz_xzzzz[k] = -g_z_0_x_0_yz_xzzzz[k] * ab_x + g_z_0_x_0_yz_xxzzzz[k];

                g_z_0_x_0_xyz_yyyyy[k] = -g_z_0_x_0_yz_yyyyy[k] * ab_x + g_z_0_x_0_yz_xyyyyy[k];

                g_z_0_x_0_xyz_yyyyz[k] = -g_z_0_x_0_yz_yyyyz[k] * ab_x + g_z_0_x_0_yz_xyyyyz[k];

                g_z_0_x_0_xyz_yyyzz[k] = -g_z_0_x_0_yz_yyyzz[k] * ab_x + g_z_0_x_0_yz_xyyyzz[k];

                g_z_0_x_0_xyz_yyzzz[k] = -g_z_0_x_0_yz_yyzzz[k] * ab_x + g_z_0_x_0_yz_xyyzzz[k];

                g_z_0_x_0_xyz_yzzzz[k] = -g_z_0_x_0_yz_yzzzz[k] * ab_x + g_z_0_x_0_yz_xyzzzz[k];

                g_z_0_x_0_xyz_zzzzz[k] = -g_z_0_x_0_yz_zzzzz[k] * ab_x + g_z_0_x_0_yz_xzzzzz[k];
            }

            /// Set up 1365-1386 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1365 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1366 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1367 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1368 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1369 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1370 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1371 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1372 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1373 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1374 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1375 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1376 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1377 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1378 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1379 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1380 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1381 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1382 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1383 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1384 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1385 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzz_xxxxx, g_z_0_x_0_xzz_xxxxy, g_z_0_x_0_xzz_xxxxz, g_z_0_x_0_xzz_xxxyy, g_z_0_x_0_xzz_xxxyz, g_z_0_x_0_xzz_xxxzz, g_z_0_x_0_xzz_xxyyy, g_z_0_x_0_xzz_xxyyz, g_z_0_x_0_xzz_xxyzz, g_z_0_x_0_xzz_xxzzz, g_z_0_x_0_xzz_xyyyy, g_z_0_x_0_xzz_xyyyz, g_z_0_x_0_xzz_xyyzz, g_z_0_x_0_xzz_xyzzz, g_z_0_x_0_xzz_xzzzz, g_z_0_x_0_xzz_yyyyy, g_z_0_x_0_xzz_yyyyz, g_z_0_x_0_xzz_yyyzz, g_z_0_x_0_xzz_yyzzz, g_z_0_x_0_xzz_yzzzz, g_z_0_x_0_xzz_zzzzz, g_z_0_x_0_zz_xxxxx, g_z_0_x_0_zz_xxxxxx, g_z_0_x_0_zz_xxxxxy, g_z_0_x_0_zz_xxxxxz, g_z_0_x_0_zz_xxxxy, g_z_0_x_0_zz_xxxxyy, g_z_0_x_0_zz_xxxxyz, g_z_0_x_0_zz_xxxxz, g_z_0_x_0_zz_xxxxzz, g_z_0_x_0_zz_xxxyy, g_z_0_x_0_zz_xxxyyy, g_z_0_x_0_zz_xxxyyz, g_z_0_x_0_zz_xxxyz, g_z_0_x_0_zz_xxxyzz, g_z_0_x_0_zz_xxxzz, g_z_0_x_0_zz_xxxzzz, g_z_0_x_0_zz_xxyyy, g_z_0_x_0_zz_xxyyyy, g_z_0_x_0_zz_xxyyyz, g_z_0_x_0_zz_xxyyz, g_z_0_x_0_zz_xxyyzz, g_z_0_x_0_zz_xxyzz, g_z_0_x_0_zz_xxyzzz, g_z_0_x_0_zz_xxzzz, g_z_0_x_0_zz_xxzzzz, g_z_0_x_0_zz_xyyyy, g_z_0_x_0_zz_xyyyyy, g_z_0_x_0_zz_xyyyyz, g_z_0_x_0_zz_xyyyz, g_z_0_x_0_zz_xyyyzz, g_z_0_x_0_zz_xyyzz, g_z_0_x_0_zz_xyyzzz, g_z_0_x_0_zz_xyzzz, g_z_0_x_0_zz_xyzzzz, g_z_0_x_0_zz_xzzzz, g_z_0_x_0_zz_xzzzzz, g_z_0_x_0_zz_yyyyy, g_z_0_x_0_zz_yyyyz, g_z_0_x_0_zz_yyyzz, g_z_0_x_0_zz_yyzzz, g_z_0_x_0_zz_yzzzz, g_z_0_x_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzz_xxxxx[k] = -g_z_0_x_0_zz_xxxxx[k] * ab_x + g_z_0_x_0_zz_xxxxxx[k];

                g_z_0_x_0_xzz_xxxxy[k] = -g_z_0_x_0_zz_xxxxy[k] * ab_x + g_z_0_x_0_zz_xxxxxy[k];

                g_z_0_x_0_xzz_xxxxz[k] = -g_z_0_x_0_zz_xxxxz[k] * ab_x + g_z_0_x_0_zz_xxxxxz[k];

                g_z_0_x_0_xzz_xxxyy[k] = -g_z_0_x_0_zz_xxxyy[k] * ab_x + g_z_0_x_0_zz_xxxxyy[k];

                g_z_0_x_0_xzz_xxxyz[k] = -g_z_0_x_0_zz_xxxyz[k] * ab_x + g_z_0_x_0_zz_xxxxyz[k];

                g_z_0_x_0_xzz_xxxzz[k] = -g_z_0_x_0_zz_xxxzz[k] * ab_x + g_z_0_x_0_zz_xxxxzz[k];

                g_z_0_x_0_xzz_xxyyy[k] = -g_z_0_x_0_zz_xxyyy[k] * ab_x + g_z_0_x_0_zz_xxxyyy[k];

                g_z_0_x_0_xzz_xxyyz[k] = -g_z_0_x_0_zz_xxyyz[k] * ab_x + g_z_0_x_0_zz_xxxyyz[k];

                g_z_0_x_0_xzz_xxyzz[k] = -g_z_0_x_0_zz_xxyzz[k] * ab_x + g_z_0_x_0_zz_xxxyzz[k];

                g_z_0_x_0_xzz_xxzzz[k] = -g_z_0_x_0_zz_xxzzz[k] * ab_x + g_z_0_x_0_zz_xxxzzz[k];

                g_z_0_x_0_xzz_xyyyy[k] = -g_z_0_x_0_zz_xyyyy[k] * ab_x + g_z_0_x_0_zz_xxyyyy[k];

                g_z_0_x_0_xzz_xyyyz[k] = -g_z_0_x_0_zz_xyyyz[k] * ab_x + g_z_0_x_0_zz_xxyyyz[k];

                g_z_0_x_0_xzz_xyyzz[k] = -g_z_0_x_0_zz_xyyzz[k] * ab_x + g_z_0_x_0_zz_xxyyzz[k];

                g_z_0_x_0_xzz_xyzzz[k] = -g_z_0_x_0_zz_xyzzz[k] * ab_x + g_z_0_x_0_zz_xxyzzz[k];

                g_z_0_x_0_xzz_xzzzz[k] = -g_z_0_x_0_zz_xzzzz[k] * ab_x + g_z_0_x_0_zz_xxzzzz[k];

                g_z_0_x_0_xzz_yyyyy[k] = -g_z_0_x_0_zz_yyyyy[k] * ab_x + g_z_0_x_0_zz_xyyyyy[k];

                g_z_0_x_0_xzz_yyyyz[k] = -g_z_0_x_0_zz_yyyyz[k] * ab_x + g_z_0_x_0_zz_xyyyyz[k];

                g_z_0_x_0_xzz_yyyzz[k] = -g_z_0_x_0_zz_yyyzz[k] * ab_x + g_z_0_x_0_zz_xyyyzz[k];

                g_z_0_x_0_xzz_yyzzz[k] = -g_z_0_x_0_zz_yyzzz[k] * ab_x + g_z_0_x_0_zz_xyyzzz[k];

                g_z_0_x_0_xzz_yzzzz[k] = -g_z_0_x_0_zz_yzzzz[k] * ab_x + g_z_0_x_0_zz_xyzzzz[k];

                g_z_0_x_0_xzz_zzzzz[k] = -g_z_0_x_0_zz_zzzzz[k] * ab_x + g_z_0_x_0_zz_xzzzzz[k];
            }

            /// Set up 1386-1407 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1386 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1387 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1388 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1389 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1390 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1391 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1392 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1393 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1394 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1395 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1396 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1397 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1398 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1399 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1400 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1401 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1402 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1403 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1404 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1405 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1406 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yy_xxxxx, g_z_0_x_0_yy_xxxxxy, g_z_0_x_0_yy_xxxxy, g_z_0_x_0_yy_xxxxyy, g_z_0_x_0_yy_xxxxyz, g_z_0_x_0_yy_xxxxz, g_z_0_x_0_yy_xxxyy, g_z_0_x_0_yy_xxxyyy, g_z_0_x_0_yy_xxxyyz, g_z_0_x_0_yy_xxxyz, g_z_0_x_0_yy_xxxyzz, g_z_0_x_0_yy_xxxzz, g_z_0_x_0_yy_xxyyy, g_z_0_x_0_yy_xxyyyy, g_z_0_x_0_yy_xxyyyz, g_z_0_x_0_yy_xxyyz, g_z_0_x_0_yy_xxyyzz, g_z_0_x_0_yy_xxyzz, g_z_0_x_0_yy_xxyzzz, g_z_0_x_0_yy_xxzzz, g_z_0_x_0_yy_xyyyy, g_z_0_x_0_yy_xyyyyy, g_z_0_x_0_yy_xyyyyz, g_z_0_x_0_yy_xyyyz, g_z_0_x_0_yy_xyyyzz, g_z_0_x_0_yy_xyyzz, g_z_0_x_0_yy_xyyzzz, g_z_0_x_0_yy_xyzzz, g_z_0_x_0_yy_xyzzzz, g_z_0_x_0_yy_xzzzz, g_z_0_x_0_yy_yyyyy, g_z_0_x_0_yy_yyyyyy, g_z_0_x_0_yy_yyyyyz, g_z_0_x_0_yy_yyyyz, g_z_0_x_0_yy_yyyyzz, g_z_0_x_0_yy_yyyzz, g_z_0_x_0_yy_yyyzzz, g_z_0_x_0_yy_yyzzz, g_z_0_x_0_yy_yyzzzz, g_z_0_x_0_yy_yzzzz, g_z_0_x_0_yy_yzzzzz, g_z_0_x_0_yy_zzzzz, g_z_0_x_0_yyy_xxxxx, g_z_0_x_0_yyy_xxxxy, g_z_0_x_0_yyy_xxxxz, g_z_0_x_0_yyy_xxxyy, g_z_0_x_0_yyy_xxxyz, g_z_0_x_0_yyy_xxxzz, g_z_0_x_0_yyy_xxyyy, g_z_0_x_0_yyy_xxyyz, g_z_0_x_0_yyy_xxyzz, g_z_0_x_0_yyy_xxzzz, g_z_0_x_0_yyy_xyyyy, g_z_0_x_0_yyy_xyyyz, g_z_0_x_0_yyy_xyyzz, g_z_0_x_0_yyy_xyzzz, g_z_0_x_0_yyy_xzzzz, g_z_0_x_0_yyy_yyyyy, g_z_0_x_0_yyy_yyyyz, g_z_0_x_0_yyy_yyyzz, g_z_0_x_0_yyy_yyzzz, g_z_0_x_0_yyy_yzzzz, g_z_0_x_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyy_xxxxx[k] = -g_z_0_x_0_yy_xxxxx[k] * ab_y + g_z_0_x_0_yy_xxxxxy[k];

                g_z_0_x_0_yyy_xxxxy[k] = -g_z_0_x_0_yy_xxxxy[k] * ab_y + g_z_0_x_0_yy_xxxxyy[k];

                g_z_0_x_0_yyy_xxxxz[k] = -g_z_0_x_0_yy_xxxxz[k] * ab_y + g_z_0_x_0_yy_xxxxyz[k];

                g_z_0_x_0_yyy_xxxyy[k] = -g_z_0_x_0_yy_xxxyy[k] * ab_y + g_z_0_x_0_yy_xxxyyy[k];

                g_z_0_x_0_yyy_xxxyz[k] = -g_z_0_x_0_yy_xxxyz[k] * ab_y + g_z_0_x_0_yy_xxxyyz[k];

                g_z_0_x_0_yyy_xxxzz[k] = -g_z_0_x_0_yy_xxxzz[k] * ab_y + g_z_0_x_0_yy_xxxyzz[k];

                g_z_0_x_0_yyy_xxyyy[k] = -g_z_0_x_0_yy_xxyyy[k] * ab_y + g_z_0_x_0_yy_xxyyyy[k];

                g_z_0_x_0_yyy_xxyyz[k] = -g_z_0_x_0_yy_xxyyz[k] * ab_y + g_z_0_x_0_yy_xxyyyz[k];

                g_z_0_x_0_yyy_xxyzz[k] = -g_z_0_x_0_yy_xxyzz[k] * ab_y + g_z_0_x_0_yy_xxyyzz[k];

                g_z_0_x_0_yyy_xxzzz[k] = -g_z_0_x_0_yy_xxzzz[k] * ab_y + g_z_0_x_0_yy_xxyzzz[k];

                g_z_0_x_0_yyy_xyyyy[k] = -g_z_0_x_0_yy_xyyyy[k] * ab_y + g_z_0_x_0_yy_xyyyyy[k];

                g_z_0_x_0_yyy_xyyyz[k] = -g_z_0_x_0_yy_xyyyz[k] * ab_y + g_z_0_x_0_yy_xyyyyz[k];

                g_z_0_x_0_yyy_xyyzz[k] = -g_z_0_x_0_yy_xyyzz[k] * ab_y + g_z_0_x_0_yy_xyyyzz[k];

                g_z_0_x_0_yyy_xyzzz[k] = -g_z_0_x_0_yy_xyzzz[k] * ab_y + g_z_0_x_0_yy_xyyzzz[k];

                g_z_0_x_0_yyy_xzzzz[k] = -g_z_0_x_0_yy_xzzzz[k] * ab_y + g_z_0_x_0_yy_xyzzzz[k];

                g_z_0_x_0_yyy_yyyyy[k] = -g_z_0_x_0_yy_yyyyy[k] * ab_y + g_z_0_x_0_yy_yyyyyy[k];

                g_z_0_x_0_yyy_yyyyz[k] = -g_z_0_x_0_yy_yyyyz[k] * ab_y + g_z_0_x_0_yy_yyyyyz[k];

                g_z_0_x_0_yyy_yyyzz[k] = -g_z_0_x_0_yy_yyyzz[k] * ab_y + g_z_0_x_0_yy_yyyyzz[k];

                g_z_0_x_0_yyy_yyzzz[k] = -g_z_0_x_0_yy_yyzzz[k] * ab_y + g_z_0_x_0_yy_yyyzzz[k];

                g_z_0_x_0_yyy_yzzzz[k] = -g_z_0_x_0_yy_yzzzz[k] * ab_y + g_z_0_x_0_yy_yyzzzz[k];

                g_z_0_x_0_yyy_zzzzz[k] = -g_z_0_x_0_yy_zzzzz[k] * ab_y + g_z_0_x_0_yy_yzzzzz[k];
            }

            /// Set up 1407-1428 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1407 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1408 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1409 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1410 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1411 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1412 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1413 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1414 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1415 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1416 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1417 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1418 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1419 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1420 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1421 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1422 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1423 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1424 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1425 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1426 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyz_xxxxx, g_z_0_x_0_yyz_xxxxy, g_z_0_x_0_yyz_xxxxz, g_z_0_x_0_yyz_xxxyy, g_z_0_x_0_yyz_xxxyz, g_z_0_x_0_yyz_xxxzz, g_z_0_x_0_yyz_xxyyy, g_z_0_x_0_yyz_xxyyz, g_z_0_x_0_yyz_xxyzz, g_z_0_x_0_yyz_xxzzz, g_z_0_x_0_yyz_xyyyy, g_z_0_x_0_yyz_xyyyz, g_z_0_x_0_yyz_xyyzz, g_z_0_x_0_yyz_xyzzz, g_z_0_x_0_yyz_xzzzz, g_z_0_x_0_yyz_yyyyy, g_z_0_x_0_yyz_yyyyz, g_z_0_x_0_yyz_yyyzz, g_z_0_x_0_yyz_yyzzz, g_z_0_x_0_yyz_yzzzz, g_z_0_x_0_yyz_zzzzz, g_z_0_x_0_yz_xxxxx, g_z_0_x_0_yz_xxxxxy, g_z_0_x_0_yz_xxxxy, g_z_0_x_0_yz_xxxxyy, g_z_0_x_0_yz_xxxxyz, g_z_0_x_0_yz_xxxxz, g_z_0_x_0_yz_xxxyy, g_z_0_x_0_yz_xxxyyy, g_z_0_x_0_yz_xxxyyz, g_z_0_x_0_yz_xxxyz, g_z_0_x_0_yz_xxxyzz, g_z_0_x_0_yz_xxxzz, g_z_0_x_0_yz_xxyyy, g_z_0_x_0_yz_xxyyyy, g_z_0_x_0_yz_xxyyyz, g_z_0_x_0_yz_xxyyz, g_z_0_x_0_yz_xxyyzz, g_z_0_x_0_yz_xxyzz, g_z_0_x_0_yz_xxyzzz, g_z_0_x_0_yz_xxzzz, g_z_0_x_0_yz_xyyyy, g_z_0_x_0_yz_xyyyyy, g_z_0_x_0_yz_xyyyyz, g_z_0_x_0_yz_xyyyz, g_z_0_x_0_yz_xyyyzz, g_z_0_x_0_yz_xyyzz, g_z_0_x_0_yz_xyyzzz, g_z_0_x_0_yz_xyzzz, g_z_0_x_0_yz_xyzzzz, g_z_0_x_0_yz_xzzzz, g_z_0_x_0_yz_yyyyy, g_z_0_x_0_yz_yyyyyy, g_z_0_x_0_yz_yyyyyz, g_z_0_x_0_yz_yyyyz, g_z_0_x_0_yz_yyyyzz, g_z_0_x_0_yz_yyyzz, g_z_0_x_0_yz_yyyzzz, g_z_0_x_0_yz_yyzzz, g_z_0_x_0_yz_yyzzzz, g_z_0_x_0_yz_yzzzz, g_z_0_x_0_yz_yzzzzz, g_z_0_x_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyz_xxxxx[k] = -g_z_0_x_0_yz_xxxxx[k] * ab_y + g_z_0_x_0_yz_xxxxxy[k];

                g_z_0_x_0_yyz_xxxxy[k] = -g_z_0_x_0_yz_xxxxy[k] * ab_y + g_z_0_x_0_yz_xxxxyy[k];

                g_z_0_x_0_yyz_xxxxz[k] = -g_z_0_x_0_yz_xxxxz[k] * ab_y + g_z_0_x_0_yz_xxxxyz[k];

                g_z_0_x_0_yyz_xxxyy[k] = -g_z_0_x_0_yz_xxxyy[k] * ab_y + g_z_0_x_0_yz_xxxyyy[k];

                g_z_0_x_0_yyz_xxxyz[k] = -g_z_0_x_0_yz_xxxyz[k] * ab_y + g_z_0_x_0_yz_xxxyyz[k];

                g_z_0_x_0_yyz_xxxzz[k] = -g_z_0_x_0_yz_xxxzz[k] * ab_y + g_z_0_x_0_yz_xxxyzz[k];

                g_z_0_x_0_yyz_xxyyy[k] = -g_z_0_x_0_yz_xxyyy[k] * ab_y + g_z_0_x_0_yz_xxyyyy[k];

                g_z_0_x_0_yyz_xxyyz[k] = -g_z_0_x_0_yz_xxyyz[k] * ab_y + g_z_0_x_0_yz_xxyyyz[k];

                g_z_0_x_0_yyz_xxyzz[k] = -g_z_0_x_0_yz_xxyzz[k] * ab_y + g_z_0_x_0_yz_xxyyzz[k];

                g_z_0_x_0_yyz_xxzzz[k] = -g_z_0_x_0_yz_xxzzz[k] * ab_y + g_z_0_x_0_yz_xxyzzz[k];

                g_z_0_x_0_yyz_xyyyy[k] = -g_z_0_x_0_yz_xyyyy[k] * ab_y + g_z_0_x_0_yz_xyyyyy[k];

                g_z_0_x_0_yyz_xyyyz[k] = -g_z_0_x_0_yz_xyyyz[k] * ab_y + g_z_0_x_0_yz_xyyyyz[k];

                g_z_0_x_0_yyz_xyyzz[k] = -g_z_0_x_0_yz_xyyzz[k] * ab_y + g_z_0_x_0_yz_xyyyzz[k];

                g_z_0_x_0_yyz_xyzzz[k] = -g_z_0_x_0_yz_xyzzz[k] * ab_y + g_z_0_x_0_yz_xyyzzz[k];

                g_z_0_x_0_yyz_xzzzz[k] = -g_z_0_x_0_yz_xzzzz[k] * ab_y + g_z_0_x_0_yz_xyzzzz[k];

                g_z_0_x_0_yyz_yyyyy[k] = -g_z_0_x_0_yz_yyyyy[k] * ab_y + g_z_0_x_0_yz_yyyyyy[k];

                g_z_0_x_0_yyz_yyyyz[k] = -g_z_0_x_0_yz_yyyyz[k] * ab_y + g_z_0_x_0_yz_yyyyyz[k];

                g_z_0_x_0_yyz_yyyzz[k] = -g_z_0_x_0_yz_yyyzz[k] * ab_y + g_z_0_x_0_yz_yyyyzz[k];

                g_z_0_x_0_yyz_yyzzz[k] = -g_z_0_x_0_yz_yyzzz[k] * ab_y + g_z_0_x_0_yz_yyyzzz[k];

                g_z_0_x_0_yyz_yzzzz[k] = -g_z_0_x_0_yz_yzzzz[k] * ab_y + g_z_0_x_0_yz_yyzzzz[k];

                g_z_0_x_0_yyz_zzzzz[k] = -g_z_0_x_0_yz_zzzzz[k] * ab_y + g_z_0_x_0_yz_yzzzzz[k];
            }

            /// Set up 1428-1449 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1428 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1429 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1430 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1431 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1432 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1433 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1434 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1435 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1436 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1437 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1438 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1439 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1440 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1441 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1442 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1443 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1444 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1445 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1446 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1447 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1448 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzz_xxxxx, g_z_0_x_0_yzz_xxxxy, g_z_0_x_0_yzz_xxxxz, g_z_0_x_0_yzz_xxxyy, g_z_0_x_0_yzz_xxxyz, g_z_0_x_0_yzz_xxxzz, g_z_0_x_0_yzz_xxyyy, g_z_0_x_0_yzz_xxyyz, g_z_0_x_0_yzz_xxyzz, g_z_0_x_0_yzz_xxzzz, g_z_0_x_0_yzz_xyyyy, g_z_0_x_0_yzz_xyyyz, g_z_0_x_0_yzz_xyyzz, g_z_0_x_0_yzz_xyzzz, g_z_0_x_0_yzz_xzzzz, g_z_0_x_0_yzz_yyyyy, g_z_0_x_0_yzz_yyyyz, g_z_0_x_0_yzz_yyyzz, g_z_0_x_0_yzz_yyzzz, g_z_0_x_0_yzz_yzzzz, g_z_0_x_0_yzz_zzzzz, g_z_0_x_0_zz_xxxxx, g_z_0_x_0_zz_xxxxxy, g_z_0_x_0_zz_xxxxy, g_z_0_x_0_zz_xxxxyy, g_z_0_x_0_zz_xxxxyz, g_z_0_x_0_zz_xxxxz, g_z_0_x_0_zz_xxxyy, g_z_0_x_0_zz_xxxyyy, g_z_0_x_0_zz_xxxyyz, g_z_0_x_0_zz_xxxyz, g_z_0_x_0_zz_xxxyzz, g_z_0_x_0_zz_xxxzz, g_z_0_x_0_zz_xxyyy, g_z_0_x_0_zz_xxyyyy, g_z_0_x_0_zz_xxyyyz, g_z_0_x_0_zz_xxyyz, g_z_0_x_0_zz_xxyyzz, g_z_0_x_0_zz_xxyzz, g_z_0_x_0_zz_xxyzzz, g_z_0_x_0_zz_xxzzz, g_z_0_x_0_zz_xyyyy, g_z_0_x_0_zz_xyyyyy, g_z_0_x_0_zz_xyyyyz, g_z_0_x_0_zz_xyyyz, g_z_0_x_0_zz_xyyyzz, g_z_0_x_0_zz_xyyzz, g_z_0_x_0_zz_xyyzzz, g_z_0_x_0_zz_xyzzz, g_z_0_x_0_zz_xyzzzz, g_z_0_x_0_zz_xzzzz, g_z_0_x_0_zz_yyyyy, g_z_0_x_0_zz_yyyyyy, g_z_0_x_0_zz_yyyyyz, g_z_0_x_0_zz_yyyyz, g_z_0_x_0_zz_yyyyzz, g_z_0_x_0_zz_yyyzz, g_z_0_x_0_zz_yyyzzz, g_z_0_x_0_zz_yyzzz, g_z_0_x_0_zz_yyzzzz, g_z_0_x_0_zz_yzzzz, g_z_0_x_0_zz_yzzzzz, g_z_0_x_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzz_xxxxx[k] = -g_z_0_x_0_zz_xxxxx[k] * ab_y + g_z_0_x_0_zz_xxxxxy[k];

                g_z_0_x_0_yzz_xxxxy[k] = -g_z_0_x_0_zz_xxxxy[k] * ab_y + g_z_0_x_0_zz_xxxxyy[k];

                g_z_0_x_0_yzz_xxxxz[k] = -g_z_0_x_0_zz_xxxxz[k] * ab_y + g_z_0_x_0_zz_xxxxyz[k];

                g_z_0_x_0_yzz_xxxyy[k] = -g_z_0_x_0_zz_xxxyy[k] * ab_y + g_z_0_x_0_zz_xxxyyy[k];

                g_z_0_x_0_yzz_xxxyz[k] = -g_z_0_x_0_zz_xxxyz[k] * ab_y + g_z_0_x_0_zz_xxxyyz[k];

                g_z_0_x_0_yzz_xxxzz[k] = -g_z_0_x_0_zz_xxxzz[k] * ab_y + g_z_0_x_0_zz_xxxyzz[k];

                g_z_0_x_0_yzz_xxyyy[k] = -g_z_0_x_0_zz_xxyyy[k] * ab_y + g_z_0_x_0_zz_xxyyyy[k];

                g_z_0_x_0_yzz_xxyyz[k] = -g_z_0_x_0_zz_xxyyz[k] * ab_y + g_z_0_x_0_zz_xxyyyz[k];

                g_z_0_x_0_yzz_xxyzz[k] = -g_z_0_x_0_zz_xxyzz[k] * ab_y + g_z_0_x_0_zz_xxyyzz[k];

                g_z_0_x_0_yzz_xxzzz[k] = -g_z_0_x_0_zz_xxzzz[k] * ab_y + g_z_0_x_0_zz_xxyzzz[k];

                g_z_0_x_0_yzz_xyyyy[k] = -g_z_0_x_0_zz_xyyyy[k] * ab_y + g_z_0_x_0_zz_xyyyyy[k];

                g_z_0_x_0_yzz_xyyyz[k] = -g_z_0_x_0_zz_xyyyz[k] * ab_y + g_z_0_x_0_zz_xyyyyz[k];

                g_z_0_x_0_yzz_xyyzz[k] = -g_z_0_x_0_zz_xyyzz[k] * ab_y + g_z_0_x_0_zz_xyyyzz[k];

                g_z_0_x_0_yzz_xyzzz[k] = -g_z_0_x_0_zz_xyzzz[k] * ab_y + g_z_0_x_0_zz_xyyzzz[k];

                g_z_0_x_0_yzz_xzzzz[k] = -g_z_0_x_0_zz_xzzzz[k] * ab_y + g_z_0_x_0_zz_xyzzzz[k];

                g_z_0_x_0_yzz_yyyyy[k] = -g_z_0_x_0_zz_yyyyy[k] * ab_y + g_z_0_x_0_zz_yyyyyy[k];

                g_z_0_x_0_yzz_yyyyz[k] = -g_z_0_x_0_zz_yyyyz[k] * ab_y + g_z_0_x_0_zz_yyyyyz[k];

                g_z_0_x_0_yzz_yyyzz[k] = -g_z_0_x_0_zz_yyyzz[k] * ab_y + g_z_0_x_0_zz_yyyyzz[k];

                g_z_0_x_0_yzz_yyzzz[k] = -g_z_0_x_0_zz_yyzzz[k] * ab_y + g_z_0_x_0_zz_yyyzzz[k];

                g_z_0_x_0_yzz_yzzzz[k] = -g_z_0_x_0_zz_yzzzz[k] * ab_y + g_z_0_x_0_zz_yyzzzz[k];

                g_z_0_x_0_yzz_zzzzz[k] = -g_z_0_x_0_zz_zzzzz[k] * ab_y + g_z_0_x_0_zz_yzzzzz[k];
            }

            /// Set up 1449-1470 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1449 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1450 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1451 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1452 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1453 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1454 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1455 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1456 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1457 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1458 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1459 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1460 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1461 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1462 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1463 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1464 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1465 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1466 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1467 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1468 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zz_xxxxx, g_0_0_x_0_zz_xxxxy, g_0_0_x_0_zz_xxxxz, g_0_0_x_0_zz_xxxyy, g_0_0_x_0_zz_xxxyz, g_0_0_x_0_zz_xxxzz, g_0_0_x_0_zz_xxyyy, g_0_0_x_0_zz_xxyyz, g_0_0_x_0_zz_xxyzz, g_0_0_x_0_zz_xxzzz, g_0_0_x_0_zz_xyyyy, g_0_0_x_0_zz_xyyyz, g_0_0_x_0_zz_xyyzz, g_0_0_x_0_zz_xyzzz, g_0_0_x_0_zz_xzzzz, g_0_0_x_0_zz_yyyyy, g_0_0_x_0_zz_yyyyz, g_0_0_x_0_zz_yyyzz, g_0_0_x_0_zz_yyzzz, g_0_0_x_0_zz_yzzzz, g_0_0_x_0_zz_zzzzz, g_z_0_x_0_zz_xxxxx, g_z_0_x_0_zz_xxxxxz, g_z_0_x_0_zz_xxxxy, g_z_0_x_0_zz_xxxxyz, g_z_0_x_0_zz_xxxxz, g_z_0_x_0_zz_xxxxzz, g_z_0_x_0_zz_xxxyy, g_z_0_x_0_zz_xxxyyz, g_z_0_x_0_zz_xxxyz, g_z_0_x_0_zz_xxxyzz, g_z_0_x_0_zz_xxxzz, g_z_0_x_0_zz_xxxzzz, g_z_0_x_0_zz_xxyyy, g_z_0_x_0_zz_xxyyyz, g_z_0_x_0_zz_xxyyz, g_z_0_x_0_zz_xxyyzz, g_z_0_x_0_zz_xxyzz, g_z_0_x_0_zz_xxyzzz, g_z_0_x_0_zz_xxzzz, g_z_0_x_0_zz_xxzzzz, g_z_0_x_0_zz_xyyyy, g_z_0_x_0_zz_xyyyyz, g_z_0_x_0_zz_xyyyz, g_z_0_x_0_zz_xyyyzz, g_z_0_x_0_zz_xyyzz, g_z_0_x_0_zz_xyyzzz, g_z_0_x_0_zz_xyzzz, g_z_0_x_0_zz_xyzzzz, g_z_0_x_0_zz_xzzzz, g_z_0_x_0_zz_xzzzzz, g_z_0_x_0_zz_yyyyy, g_z_0_x_0_zz_yyyyyz, g_z_0_x_0_zz_yyyyz, g_z_0_x_0_zz_yyyyzz, g_z_0_x_0_zz_yyyzz, g_z_0_x_0_zz_yyyzzz, g_z_0_x_0_zz_yyzzz, g_z_0_x_0_zz_yyzzzz, g_z_0_x_0_zz_yzzzz, g_z_0_x_0_zz_yzzzzz, g_z_0_x_0_zz_zzzzz, g_z_0_x_0_zz_zzzzzz, g_z_0_x_0_zzz_xxxxx, g_z_0_x_0_zzz_xxxxy, g_z_0_x_0_zzz_xxxxz, g_z_0_x_0_zzz_xxxyy, g_z_0_x_0_zzz_xxxyz, g_z_0_x_0_zzz_xxxzz, g_z_0_x_0_zzz_xxyyy, g_z_0_x_0_zzz_xxyyz, g_z_0_x_0_zzz_xxyzz, g_z_0_x_0_zzz_xxzzz, g_z_0_x_0_zzz_xyyyy, g_z_0_x_0_zzz_xyyyz, g_z_0_x_0_zzz_xyyzz, g_z_0_x_0_zzz_xyzzz, g_z_0_x_0_zzz_xzzzz, g_z_0_x_0_zzz_yyyyy, g_z_0_x_0_zzz_yyyyz, g_z_0_x_0_zzz_yyyzz, g_z_0_x_0_zzz_yyzzz, g_z_0_x_0_zzz_yzzzz, g_z_0_x_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzz_xxxxx[k] = -g_0_0_x_0_zz_xxxxx[k] - g_z_0_x_0_zz_xxxxx[k] * ab_z + g_z_0_x_0_zz_xxxxxz[k];

                g_z_0_x_0_zzz_xxxxy[k] = -g_0_0_x_0_zz_xxxxy[k] - g_z_0_x_0_zz_xxxxy[k] * ab_z + g_z_0_x_0_zz_xxxxyz[k];

                g_z_0_x_0_zzz_xxxxz[k] = -g_0_0_x_0_zz_xxxxz[k] - g_z_0_x_0_zz_xxxxz[k] * ab_z + g_z_0_x_0_zz_xxxxzz[k];

                g_z_0_x_0_zzz_xxxyy[k] = -g_0_0_x_0_zz_xxxyy[k] - g_z_0_x_0_zz_xxxyy[k] * ab_z + g_z_0_x_0_zz_xxxyyz[k];

                g_z_0_x_0_zzz_xxxyz[k] = -g_0_0_x_0_zz_xxxyz[k] - g_z_0_x_0_zz_xxxyz[k] * ab_z + g_z_0_x_0_zz_xxxyzz[k];

                g_z_0_x_0_zzz_xxxzz[k] = -g_0_0_x_0_zz_xxxzz[k] - g_z_0_x_0_zz_xxxzz[k] * ab_z + g_z_0_x_0_zz_xxxzzz[k];

                g_z_0_x_0_zzz_xxyyy[k] = -g_0_0_x_0_zz_xxyyy[k] - g_z_0_x_0_zz_xxyyy[k] * ab_z + g_z_0_x_0_zz_xxyyyz[k];

                g_z_0_x_0_zzz_xxyyz[k] = -g_0_0_x_0_zz_xxyyz[k] - g_z_0_x_0_zz_xxyyz[k] * ab_z + g_z_0_x_0_zz_xxyyzz[k];

                g_z_0_x_0_zzz_xxyzz[k] = -g_0_0_x_0_zz_xxyzz[k] - g_z_0_x_0_zz_xxyzz[k] * ab_z + g_z_0_x_0_zz_xxyzzz[k];

                g_z_0_x_0_zzz_xxzzz[k] = -g_0_0_x_0_zz_xxzzz[k] - g_z_0_x_0_zz_xxzzz[k] * ab_z + g_z_0_x_0_zz_xxzzzz[k];

                g_z_0_x_0_zzz_xyyyy[k] = -g_0_0_x_0_zz_xyyyy[k] - g_z_0_x_0_zz_xyyyy[k] * ab_z + g_z_0_x_0_zz_xyyyyz[k];

                g_z_0_x_0_zzz_xyyyz[k] = -g_0_0_x_0_zz_xyyyz[k] - g_z_0_x_0_zz_xyyyz[k] * ab_z + g_z_0_x_0_zz_xyyyzz[k];

                g_z_0_x_0_zzz_xyyzz[k] = -g_0_0_x_0_zz_xyyzz[k] - g_z_0_x_0_zz_xyyzz[k] * ab_z + g_z_0_x_0_zz_xyyzzz[k];

                g_z_0_x_0_zzz_xyzzz[k] = -g_0_0_x_0_zz_xyzzz[k] - g_z_0_x_0_zz_xyzzz[k] * ab_z + g_z_0_x_0_zz_xyzzzz[k];

                g_z_0_x_0_zzz_xzzzz[k] = -g_0_0_x_0_zz_xzzzz[k] - g_z_0_x_0_zz_xzzzz[k] * ab_z + g_z_0_x_0_zz_xzzzzz[k];

                g_z_0_x_0_zzz_yyyyy[k] = -g_0_0_x_0_zz_yyyyy[k] - g_z_0_x_0_zz_yyyyy[k] * ab_z + g_z_0_x_0_zz_yyyyyz[k];

                g_z_0_x_0_zzz_yyyyz[k] = -g_0_0_x_0_zz_yyyyz[k] - g_z_0_x_0_zz_yyyyz[k] * ab_z + g_z_0_x_0_zz_yyyyzz[k];

                g_z_0_x_0_zzz_yyyzz[k] = -g_0_0_x_0_zz_yyyzz[k] - g_z_0_x_0_zz_yyyzz[k] * ab_z + g_z_0_x_0_zz_yyyzzz[k];

                g_z_0_x_0_zzz_yyzzz[k] = -g_0_0_x_0_zz_yyzzz[k] - g_z_0_x_0_zz_yyzzz[k] * ab_z + g_z_0_x_0_zz_yyzzzz[k];

                g_z_0_x_0_zzz_yzzzz[k] = -g_0_0_x_0_zz_yzzzz[k] - g_z_0_x_0_zz_yzzzz[k] * ab_z + g_z_0_x_0_zz_yzzzzz[k];

                g_z_0_x_0_zzz_zzzzz[k] = -g_0_0_x_0_zz_zzzzz[k] - g_z_0_x_0_zz_zzzzz[k] * ab_z + g_z_0_x_0_zz_zzzzzz[k];
            }

            /// Set up 1470-1491 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 1470 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 1471 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 1472 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 1473 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 1474 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 1475 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 1476 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 1477 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 1478 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 1479 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 1480 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 1481 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 1482 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 1483 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 1484 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 1485 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 1486 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 1487 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 1488 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 1489 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 1490 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xx_xxxxx, g_z_0_y_0_xx_xxxxxx, g_z_0_y_0_xx_xxxxxy, g_z_0_y_0_xx_xxxxxz, g_z_0_y_0_xx_xxxxy, g_z_0_y_0_xx_xxxxyy, g_z_0_y_0_xx_xxxxyz, g_z_0_y_0_xx_xxxxz, g_z_0_y_0_xx_xxxxzz, g_z_0_y_0_xx_xxxyy, g_z_0_y_0_xx_xxxyyy, g_z_0_y_0_xx_xxxyyz, g_z_0_y_0_xx_xxxyz, g_z_0_y_0_xx_xxxyzz, g_z_0_y_0_xx_xxxzz, g_z_0_y_0_xx_xxxzzz, g_z_0_y_0_xx_xxyyy, g_z_0_y_0_xx_xxyyyy, g_z_0_y_0_xx_xxyyyz, g_z_0_y_0_xx_xxyyz, g_z_0_y_0_xx_xxyyzz, g_z_0_y_0_xx_xxyzz, g_z_0_y_0_xx_xxyzzz, g_z_0_y_0_xx_xxzzz, g_z_0_y_0_xx_xxzzzz, g_z_0_y_0_xx_xyyyy, g_z_0_y_0_xx_xyyyyy, g_z_0_y_0_xx_xyyyyz, g_z_0_y_0_xx_xyyyz, g_z_0_y_0_xx_xyyyzz, g_z_0_y_0_xx_xyyzz, g_z_0_y_0_xx_xyyzzz, g_z_0_y_0_xx_xyzzz, g_z_0_y_0_xx_xyzzzz, g_z_0_y_0_xx_xzzzz, g_z_0_y_0_xx_xzzzzz, g_z_0_y_0_xx_yyyyy, g_z_0_y_0_xx_yyyyz, g_z_0_y_0_xx_yyyzz, g_z_0_y_0_xx_yyzzz, g_z_0_y_0_xx_yzzzz, g_z_0_y_0_xx_zzzzz, g_z_0_y_0_xxx_xxxxx, g_z_0_y_0_xxx_xxxxy, g_z_0_y_0_xxx_xxxxz, g_z_0_y_0_xxx_xxxyy, g_z_0_y_0_xxx_xxxyz, g_z_0_y_0_xxx_xxxzz, g_z_0_y_0_xxx_xxyyy, g_z_0_y_0_xxx_xxyyz, g_z_0_y_0_xxx_xxyzz, g_z_0_y_0_xxx_xxzzz, g_z_0_y_0_xxx_xyyyy, g_z_0_y_0_xxx_xyyyz, g_z_0_y_0_xxx_xyyzz, g_z_0_y_0_xxx_xyzzz, g_z_0_y_0_xxx_xzzzz, g_z_0_y_0_xxx_yyyyy, g_z_0_y_0_xxx_yyyyz, g_z_0_y_0_xxx_yyyzz, g_z_0_y_0_xxx_yyzzz, g_z_0_y_0_xxx_yzzzz, g_z_0_y_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxx_xxxxx[k] = -g_z_0_y_0_xx_xxxxx[k] * ab_x + g_z_0_y_0_xx_xxxxxx[k];

                g_z_0_y_0_xxx_xxxxy[k] = -g_z_0_y_0_xx_xxxxy[k] * ab_x + g_z_0_y_0_xx_xxxxxy[k];

                g_z_0_y_0_xxx_xxxxz[k] = -g_z_0_y_0_xx_xxxxz[k] * ab_x + g_z_0_y_0_xx_xxxxxz[k];

                g_z_0_y_0_xxx_xxxyy[k] = -g_z_0_y_0_xx_xxxyy[k] * ab_x + g_z_0_y_0_xx_xxxxyy[k];

                g_z_0_y_0_xxx_xxxyz[k] = -g_z_0_y_0_xx_xxxyz[k] * ab_x + g_z_0_y_0_xx_xxxxyz[k];

                g_z_0_y_0_xxx_xxxzz[k] = -g_z_0_y_0_xx_xxxzz[k] * ab_x + g_z_0_y_0_xx_xxxxzz[k];

                g_z_0_y_0_xxx_xxyyy[k] = -g_z_0_y_0_xx_xxyyy[k] * ab_x + g_z_0_y_0_xx_xxxyyy[k];

                g_z_0_y_0_xxx_xxyyz[k] = -g_z_0_y_0_xx_xxyyz[k] * ab_x + g_z_0_y_0_xx_xxxyyz[k];

                g_z_0_y_0_xxx_xxyzz[k] = -g_z_0_y_0_xx_xxyzz[k] * ab_x + g_z_0_y_0_xx_xxxyzz[k];

                g_z_0_y_0_xxx_xxzzz[k] = -g_z_0_y_0_xx_xxzzz[k] * ab_x + g_z_0_y_0_xx_xxxzzz[k];

                g_z_0_y_0_xxx_xyyyy[k] = -g_z_0_y_0_xx_xyyyy[k] * ab_x + g_z_0_y_0_xx_xxyyyy[k];

                g_z_0_y_0_xxx_xyyyz[k] = -g_z_0_y_0_xx_xyyyz[k] * ab_x + g_z_0_y_0_xx_xxyyyz[k];

                g_z_0_y_0_xxx_xyyzz[k] = -g_z_0_y_0_xx_xyyzz[k] * ab_x + g_z_0_y_0_xx_xxyyzz[k];

                g_z_0_y_0_xxx_xyzzz[k] = -g_z_0_y_0_xx_xyzzz[k] * ab_x + g_z_0_y_0_xx_xxyzzz[k];

                g_z_0_y_0_xxx_xzzzz[k] = -g_z_0_y_0_xx_xzzzz[k] * ab_x + g_z_0_y_0_xx_xxzzzz[k];

                g_z_0_y_0_xxx_yyyyy[k] = -g_z_0_y_0_xx_yyyyy[k] * ab_x + g_z_0_y_0_xx_xyyyyy[k];

                g_z_0_y_0_xxx_yyyyz[k] = -g_z_0_y_0_xx_yyyyz[k] * ab_x + g_z_0_y_0_xx_xyyyyz[k];

                g_z_0_y_0_xxx_yyyzz[k] = -g_z_0_y_0_xx_yyyzz[k] * ab_x + g_z_0_y_0_xx_xyyyzz[k];

                g_z_0_y_0_xxx_yyzzz[k] = -g_z_0_y_0_xx_yyzzz[k] * ab_x + g_z_0_y_0_xx_xyyzzz[k];

                g_z_0_y_0_xxx_yzzzz[k] = -g_z_0_y_0_xx_yzzzz[k] * ab_x + g_z_0_y_0_xx_xyzzzz[k];

                g_z_0_y_0_xxx_zzzzz[k] = -g_z_0_y_0_xx_zzzzz[k] * ab_x + g_z_0_y_0_xx_xzzzzz[k];
            }

            /// Set up 1491-1512 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 1491 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 1492 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 1493 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 1494 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 1495 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 1496 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 1497 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 1498 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 1499 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 1500 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 1501 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 1502 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 1503 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 1504 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 1505 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 1506 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 1507 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 1508 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 1509 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 1510 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxy_xxxxx, g_z_0_y_0_xxy_xxxxy, g_z_0_y_0_xxy_xxxxz, g_z_0_y_0_xxy_xxxyy, g_z_0_y_0_xxy_xxxyz, g_z_0_y_0_xxy_xxxzz, g_z_0_y_0_xxy_xxyyy, g_z_0_y_0_xxy_xxyyz, g_z_0_y_0_xxy_xxyzz, g_z_0_y_0_xxy_xxzzz, g_z_0_y_0_xxy_xyyyy, g_z_0_y_0_xxy_xyyyz, g_z_0_y_0_xxy_xyyzz, g_z_0_y_0_xxy_xyzzz, g_z_0_y_0_xxy_xzzzz, g_z_0_y_0_xxy_yyyyy, g_z_0_y_0_xxy_yyyyz, g_z_0_y_0_xxy_yyyzz, g_z_0_y_0_xxy_yyzzz, g_z_0_y_0_xxy_yzzzz, g_z_0_y_0_xxy_zzzzz, g_z_0_y_0_xy_xxxxx, g_z_0_y_0_xy_xxxxxx, g_z_0_y_0_xy_xxxxxy, g_z_0_y_0_xy_xxxxxz, g_z_0_y_0_xy_xxxxy, g_z_0_y_0_xy_xxxxyy, g_z_0_y_0_xy_xxxxyz, g_z_0_y_0_xy_xxxxz, g_z_0_y_0_xy_xxxxzz, g_z_0_y_0_xy_xxxyy, g_z_0_y_0_xy_xxxyyy, g_z_0_y_0_xy_xxxyyz, g_z_0_y_0_xy_xxxyz, g_z_0_y_0_xy_xxxyzz, g_z_0_y_0_xy_xxxzz, g_z_0_y_0_xy_xxxzzz, g_z_0_y_0_xy_xxyyy, g_z_0_y_0_xy_xxyyyy, g_z_0_y_0_xy_xxyyyz, g_z_0_y_0_xy_xxyyz, g_z_0_y_0_xy_xxyyzz, g_z_0_y_0_xy_xxyzz, g_z_0_y_0_xy_xxyzzz, g_z_0_y_0_xy_xxzzz, g_z_0_y_0_xy_xxzzzz, g_z_0_y_0_xy_xyyyy, g_z_0_y_0_xy_xyyyyy, g_z_0_y_0_xy_xyyyyz, g_z_0_y_0_xy_xyyyz, g_z_0_y_0_xy_xyyyzz, g_z_0_y_0_xy_xyyzz, g_z_0_y_0_xy_xyyzzz, g_z_0_y_0_xy_xyzzz, g_z_0_y_0_xy_xyzzzz, g_z_0_y_0_xy_xzzzz, g_z_0_y_0_xy_xzzzzz, g_z_0_y_0_xy_yyyyy, g_z_0_y_0_xy_yyyyz, g_z_0_y_0_xy_yyyzz, g_z_0_y_0_xy_yyzzz, g_z_0_y_0_xy_yzzzz, g_z_0_y_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxy_xxxxx[k] = -g_z_0_y_0_xy_xxxxx[k] * ab_x + g_z_0_y_0_xy_xxxxxx[k];

                g_z_0_y_0_xxy_xxxxy[k] = -g_z_0_y_0_xy_xxxxy[k] * ab_x + g_z_0_y_0_xy_xxxxxy[k];

                g_z_0_y_0_xxy_xxxxz[k] = -g_z_0_y_0_xy_xxxxz[k] * ab_x + g_z_0_y_0_xy_xxxxxz[k];

                g_z_0_y_0_xxy_xxxyy[k] = -g_z_0_y_0_xy_xxxyy[k] * ab_x + g_z_0_y_0_xy_xxxxyy[k];

                g_z_0_y_0_xxy_xxxyz[k] = -g_z_0_y_0_xy_xxxyz[k] * ab_x + g_z_0_y_0_xy_xxxxyz[k];

                g_z_0_y_0_xxy_xxxzz[k] = -g_z_0_y_0_xy_xxxzz[k] * ab_x + g_z_0_y_0_xy_xxxxzz[k];

                g_z_0_y_0_xxy_xxyyy[k] = -g_z_0_y_0_xy_xxyyy[k] * ab_x + g_z_0_y_0_xy_xxxyyy[k];

                g_z_0_y_0_xxy_xxyyz[k] = -g_z_0_y_0_xy_xxyyz[k] * ab_x + g_z_0_y_0_xy_xxxyyz[k];

                g_z_0_y_0_xxy_xxyzz[k] = -g_z_0_y_0_xy_xxyzz[k] * ab_x + g_z_0_y_0_xy_xxxyzz[k];

                g_z_0_y_0_xxy_xxzzz[k] = -g_z_0_y_0_xy_xxzzz[k] * ab_x + g_z_0_y_0_xy_xxxzzz[k];

                g_z_0_y_0_xxy_xyyyy[k] = -g_z_0_y_0_xy_xyyyy[k] * ab_x + g_z_0_y_0_xy_xxyyyy[k];

                g_z_0_y_0_xxy_xyyyz[k] = -g_z_0_y_0_xy_xyyyz[k] * ab_x + g_z_0_y_0_xy_xxyyyz[k];

                g_z_0_y_0_xxy_xyyzz[k] = -g_z_0_y_0_xy_xyyzz[k] * ab_x + g_z_0_y_0_xy_xxyyzz[k];

                g_z_0_y_0_xxy_xyzzz[k] = -g_z_0_y_0_xy_xyzzz[k] * ab_x + g_z_0_y_0_xy_xxyzzz[k];

                g_z_0_y_0_xxy_xzzzz[k] = -g_z_0_y_0_xy_xzzzz[k] * ab_x + g_z_0_y_0_xy_xxzzzz[k];

                g_z_0_y_0_xxy_yyyyy[k] = -g_z_0_y_0_xy_yyyyy[k] * ab_x + g_z_0_y_0_xy_xyyyyy[k];

                g_z_0_y_0_xxy_yyyyz[k] = -g_z_0_y_0_xy_yyyyz[k] * ab_x + g_z_0_y_0_xy_xyyyyz[k];

                g_z_0_y_0_xxy_yyyzz[k] = -g_z_0_y_0_xy_yyyzz[k] * ab_x + g_z_0_y_0_xy_xyyyzz[k];

                g_z_0_y_0_xxy_yyzzz[k] = -g_z_0_y_0_xy_yyzzz[k] * ab_x + g_z_0_y_0_xy_xyyzzz[k];

                g_z_0_y_0_xxy_yzzzz[k] = -g_z_0_y_0_xy_yzzzz[k] * ab_x + g_z_0_y_0_xy_xyzzzz[k];

                g_z_0_y_0_xxy_zzzzz[k] = -g_z_0_y_0_xy_zzzzz[k] * ab_x + g_z_0_y_0_xy_xzzzzz[k];
            }

            /// Set up 1512-1533 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 1512 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 1513 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 1514 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 1515 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 1516 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 1517 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 1518 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 1519 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 1520 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 1521 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 1522 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 1523 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 1524 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 1525 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 1526 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 1527 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 1528 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 1529 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 1530 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 1531 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 1532 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxz_xxxxx, g_z_0_y_0_xxz_xxxxy, g_z_0_y_0_xxz_xxxxz, g_z_0_y_0_xxz_xxxyy, g_z_0_y_0_xxz_xxxyz, g_z_0_y_0_xxz_xxxzz, g_z_0_y_0_xxz_xxyyy, g_z_0_y_0_xxz_xxyyz, g_z_0_y_0_xxz_xxyzz, g_z_0_y_0_xxz_xxzzz, g_z_0_y_0_xxz_xyyyy, g_z_0_y_0_xxz_xyyyz, g_z_0_y_0_xxz_xyyzz, g_z_0_y_0_xxz_xyzzz, g_z_0_y_0_xxz_xzzzz, g_z_0_y_0_xxz_yyyyy, g_z_0_y_0_xxz_yyyyz, g_z_0_y_0_xxz_yyyzz, g_z_0_y_0_xxz_yyzzz, g_z_0_y_0_xxz_yzzzz, g_z_0_y_0_xxz_zzzzz, g_z_0_y_0_xz_xxxxx, g_z_0_y_0_xz_xxxxxx, g_z_0_y_0_xz_xxxxxy, g_z_0_y_0_xz_xxxxxz, g_z_0_y_0_xz_xxxxy, g_z_0_y_0_xz_xxxxyy, g_z_0_y_0_xz_xxxxyz, g_z_0_y_0_xz_xxxxz, g_z_0_y_0_xz_xxxxzz, g_z_0_y_0_xz_xxxyy, g_z_0_y_0_xz_xxxyyy, g_z_0_y_0_xz_xxxyyz, g_z_0_y_0_xz_xxxyz, g_z_0_y_0_xz_xxxyzz, g_z_0_y_0_xz_xxxzz, g_z_0_y_0_xz_xxxzzz, g_z_0_y_0_xz_xxyyy, g_z_0_y_0_xz_xxyyyy, g_z_0_y_0_xz_xxyyyz, g_z_0_y_0_xz_xxyyz, g_z_0_y_0_xz_xxyyzz, g_z_0_y_0_xz_xxyzz, g_z_0_y_0_xz_xxyzzz, g_z_0_y_0_xz_xxzzz, g_z_0_y_0_xz_xxzzzz, g_z_0_y_0_xz_xyyyy, g_z_0_y_0_xz_xyyyyy, g_z_0_y_0_xz_xyyyyz, g_z_0_y_0_xz_xyyyz, g_z_0_y_0_xz_xyyyzz, g_z_0_y_0_xz_xyyzz, g_z_0_y_0_xz_xyyzzz, g_z_0_y_0_xz_xyzzz, g_z_0_y_0_xz_xyzzzz, g_z_0_y_0_xz_xzzzz, g_z_0_y_0_xz_xzzzzz, g_z_0_y_0_xz_yyyyy, g_z_0_y_0_xz_yyyyz, g_z_0_y_0_xz_yyyzz, g_z_0_y_0_xz_yyzzz, g_z_0_y_0_xz_yzzzz, g_z_0_y_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxz_xxxxx[k] = -g_z_0_y_0_xz_xxxxx[k] * ab_x + g_z_0_y_0_xz_xxxxxx[k];

                g_z_0_y_0_xxz_xxxxy[k] = -g_z_0_y_0_xz_xxxxy[k] * ab_x + g_z_0_y_0_xz_xxxxxy[k];

                g_z_0_y_0_xxz_xxxxz[k] = -g_z_0_y_0_xz_xxxxz[k] * ab_x + g_z_0_y_0_xz_xxxxxz[k];

                g_z_0_y_0_xxz_xxxyy[k] = -g_z_0_y_0_xz_xxxyy[k] * ab_x + g_z_0_y_0_xz_xxxxyy[k];

                g_z_0_y_0_xxz_xxxyz[k] = -g_z_0_y_0_xz_xxxyz[k] * ab_x + g_z_0_y_0_xz_xxxxyz[k];

                g_z_0_y_0_xxz_xxxzz[k] = -g_z_0_y_0_xz_xxxzz[k] * ab_x + g_z_0_y_0_xz_xxxxzz[k];

                g_z_0_y_0_xxz_xxyyy[k] = -g_z_0_y_0_xz_xxyyy[k] * ab_x + g_z_0_y_0_xz_xxxyyy[k];

                g_z_0_y_0_xxz_xxyyz[k] = -g_z_0_y_0_xz_xxyyz[k] * ab_x + g_z_0_y_0_xz_xxxyyz[k];

                g_z_0_y_0_xxz_xxyzz[k] = -g_z_0_y_0_xz_xxyzz[k] * ab_x + g_z_0_y_0_xz_xxxyzz[k];

                g_z_0_y_0_xxz_xxzzz[k] = -g_z_0_y_0_xz_xxzzz[k] * ab_x + g_z_0_y_0_xz_xxxzzz[k];

                g_z_0_y_0_xxz_xyyyy[k] = -g_z_0_y_0_xz_xyyyy[k] * ab_x + g_z_0_y_0_xz_xxyyyy[k];

                g_z_0_y_0_xxz_xyyyz[k] = -g_z_0_y_0_xz_xyyyz[k] * ab_x + g_z_0_y_0_xz_xxyyyz[k];

                g_z_0_y_0_xxz_xyyzz[k] = -g_z_0_y_0_xz_xyyzz[k] * ab_x + g_z_0_y_0_xz_xxyyzz[k];

                g_z_0_y_0_xxz_xyzzz[k] = -g_z_0_y_0_xz_xyzzz[k] * ab_x + g_z_0_y_0_xz_xxyzzz[k];

                g_z_0_y_0_xxz_xzzzz[k] = -g_z_0_y_0_xz_xzzzz[k] * ab_x + g_z_0_y_0_xz_xxzzzz[k];

                g_z_0_y_0_xxz_yyyyy[k] = -g_z_0_y_0_xz_yyyyy[k] * ab_x + g_z_0_y_0_xz_xyyyyy[k];

                g_z_0_y_0_xxz_yyyyz[k] = -g_z_0_y_0_xz_yyyyz[k] * ab_x + g_z_0_y_0_xz_xyyyyz[k];

                g_z_0_y_0_xxz_yyyzz[k] = -g_z_0_y_0_xz_yyyzz[k] * ab_x + g_z_0_y_0_xz_xyyyzz[k];

                g_z_0_y_0_xxz_yyzzz[k] = -g_z_0_y_0_xz_yyzzz[k] * ab_x + g_z_0_y_0_xz_xyyzzz[k];

                g_z_0_y_0_xxz_yzzzz[k] = -g_z_0_y_0_xz_yzzzz[k] * ab_x + g_z_0_y_0_xz_xyzzzz[k];

                g_z_0_y_0_xxz_zzzzz[k] = -g_z_0_y_0_xz_zzzzz[k] * ab_x + g_z_0_y_0_xz_xzzzzz[k];
            }

            /// Set up 1533-1554 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1533 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1534 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1535 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1536 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1537 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1538 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1539 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1540 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1541 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1542 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1543 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1544 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1545 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1546 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1547 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1548 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1549 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1550 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1551 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1552 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1553 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyy_xxxxx, g_z_0_y_0_xyy_xxxxy, g_z_0_y_0_xyy_xxxxz, g_z_0_y_0_xyy_xxxyy, g_z_0_y_0_xyy_xxxyz, g_z_0_y_0_xyy_xxxzz, g_z_0_y_0_xyy_xxyyy, g_z_0_y_0_xyy_xxyyz, g_z_0_y_0_xyy_xxyzz, g_z_0_y_0_xyy_xxzzz, g_z_0_y_0_xyy_xyyyy, g_z_0_y_0_xyy_xyyyz, g_z_0_y_0_xyy_xyyzz, g_z_0_y_0_xyy_xyzzz, g_z_0_y_0_xyy_xzzzz, g_z_0_y_0_xyy_yyyyy, g_z_0_y_0_xyy_yyyyz, g_z_0_y_0_xyy_yyyzz, g_z_0_y_0_xyy_yyzzz, g_z_0_y_0_xyy_yzzzz, g_z_0_y_0_xyy_zzzzz, g_z_0_y_0_yy_xxxxx, g_z_0_y_0_yy_xxxxxx, g_z_0_y_0_yy_xxxxxy, g_z_0_y_0_yy_xxxxxz, g_z_0_y_0_yy_xxxxy, g_z_0_y_0_yy_xxxxyy, g_z_0_y_0_yy_xxxxyz, g_z_0_y_0_yy_xxxxz, g_z_0_y_0_yy_xxxxzz, g_z_0_y_0_yy_xxxyy, g_z_0_y_0_yy_xxxyyy, g_z_0_y_0_yy_xxxyyz, g_z_0_y_0_yy_xxxyz, g_z_0_y_0_yy_xxxyzz, g_z_0_y_0_yy_xxxzz, g_z_0_y_0_yy_xxxzzz, g_z_0_y_0_yy_xxyyy, g_z_0_y_0_yy_xxyyyy, g_z_0_y_0_yy_xxyyyz, g_z_0_y_0_yy_xxyyz, g_z_0_y_0_yy_xxyyzz, g_z_0_y_0_yy_xxyzz, g_z_0_y_0_yy_xxyzzz, g_z_0_y_0_yy_xxzzz, g_z_0_y_0_yy_xxzzzz, g_z_0_y_0_yy_xyyyy, g_z_0_y_0_yy_xyyyyy, g_z_0_y_0_yy_xyyyyz, g_z_0_y_0_yy_xyyyz, g_z_0_y_0_yy_xyyyzz, g_z_0_y_0_yy_xyyzz, g_z_0_y_0_yy_xyyzzz, g_z_0_y_0_yy_xyzzz, g_z_0_y_0_yy_xyzzzz, g_z_0_y_0_yy_xzzzz, g_z_0_y_0_yy_xzzzzz, g_z_0_y_0_yy_yyyyy, g_z_0_y_0_yy_yyyyz, g_z_0_y_0_yy_yyyzz, g_z_0_y_0_yy_yyzzz, g_z_0_y_0_yy_yzzzz, g_z_0_y_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyy_xxxxx[k] = -g_z_0_y_0_yy_xxxxx[k] * ab_x + g_z_0_y_0_yy_xxxxxx[k];

                g_z_0_y_0_xyy_xxxxy[k] = -g_z_0_y_0_yy_xxxxy[k] * ab_x + g_z_0_y_0_yy_xxxxxy[k];

                g_z_0_y_0_xyy_xxxxz[k] = -g_z_0_y_0_yy_xxxxz[k] * ab_x + g_z_0_y_0_yy_xxxxxz[k];

                g_z_0_y_0_xyy_xxxyy[k] = -g_z_0_y_0_yy_xxxyy[k] * ab_x + g_z_0_y_0_yy_xxxxyy[k];

                g_z_0_y_0_xyy_xxxyz[k] = -g_z_0_y_0_yy_xxxyz[k] * ab_x + g_z_0_y_0_yy_xxxxyz[k];

                g_z_0_y_0_xyy_xxxzz[k] = -g_z_0_y_0_yy_xxxzz[k] * ab_x + g_z_0_y_0_yy_xxxxzz[k];

                g_z_0_y_0_xyy_xxyyy[k] = -g_z_0_y_0_yy_xxyyy[k] * ab_x + g_z_0_y_0_yy_xxxyyy[k];

                g_z_0_y_0_xyy_xxyyz[k] = -g_z_0_y_0_yy_xxyyz[k] * ab_x + g_z_0_y_0_yy_xxxyyz[k];

                g_z_0_y_0_xyy_xxyzz[k] = -g_z_0_y_0_yy_xxyzz[k] * ab_x + g_z_0_y_0_yy_xxxyzz[k];

                g_z_0_y_0_xyy_xxzzz[k] = -g_z_0_y_0_yy_xxzzz[k] * ab_x + g_z_0_y_0_yy_xxxzzz[k];

                g_z_0_y_0_xyy_xyyyy[k] = -g_z_0_y_0_yy_xyyyy[k] * ab_x + g_z_0_y_0_yy_xxyyyy[k];

                g_z_0_y_0_xyy_xyyyz[k] = -g_z_0_y_0_yy_xyyyz[k] * ab_x + g_z_0_y_0_yy_xxyyyz[k];

                g_z_0_y_0_xyy_xyyzz[k] = -g_z_0_y_0_yy_xyyzz[k] * ab_x + g_z_0_y_0_yy_xxyyzz[k];

                g_z_0_y_0_xyy_xyzzz[k] = -g_z_0_y_0_yy_xyzzz[k] * ab_x + g_z_0_y_0_yy_xxyzzz[k];

                g_z_0_y_0_xyy_xzzzz[k] = -g_z_0_y_0_yy_xzzzz[k] * ab_x + g_z_0_y_0_yy_xxzzzz[k];

                g_z_0_y_0_xyy_yyyyy[k] = -g_z_0_y_0_yy_yyyyy[k] * ab_x + g_z_0_y_0_yy_xyyyyy[k];

                g_z_0_y_0_xyy_yyyyz[k] = -g_z_0_y_0_yy_yyyyz[k] * ab_x + g_z_0_y_0_yy_xyyyyz[k];

                g_z_0_y_0_xyy_yyyzz[k] = -g_z_0_y_0_yy_yyyzz[k] * ab_x + g_z_0_y_0_yy_xyyyzz[k];

                g_z_0_y_0_xyy_yyzzz[k] = -g_z_0_y_0_yy_yyzzz[k] * ab_x + g_z_0_y_0_yy_xyyzzz[k];

                g_z_0_y_0_xyy_yzzzz[k] = -g_z_0_y_0_yy_yzzzz[k] * ab_x + g_z_0_y_0_yy_xyzzzz[k];

                g_z_0_y_0_xyy_zzzzz[k] = -g_z_0_y_0_yy_zzzzz[k] * ab_x + g_z_0_y_0_yy_xzzzzz[k];
            }

            /// Set up 1554-1575 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1554 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1555 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1556 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1557 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1558 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1559 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1560 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1561 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1562 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1563 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1564 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1565 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1566 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1567 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1568 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1569 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1570 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1571 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1572 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1573 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1574 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyz_xxxxx, g_z_0_y_0_xyz_xxxxy, g_z_0_y_0_xyz_xxxxz, g_z_0_y_0_xyz_xxxyy, g_z_0_y_0_xyz_xxxyz, g_z_0_y_0_xyz_xxxzz, g_z_0_y_0_xyz_xxyyy, g_z_0_y_0_xyz_xxyyz, g_z_0_y_0_xyz_xxyzz, g_z_0_y_0_xyz_xxzzz, g_z_0_y_0_xyz_xyyyy, g_z_0_y_0_xyz_xyyyz, g_z_0_y_0_xyz_xyyzz, g_z_0_y_0_xyz_xyzzz, g_z_0_y_0_xyz_xzzzz, g_z_0_y_0_xyz_yyyyy, g_z_0_y_0_xyz_yyyyz, g_z_0_y_0_xyz_yyyzz, g_z_0_y_0_xyz_yyzzz, g_z_0_y_0_xyz_yzzzz, g_z_0_y_0_xyz_zzzzz, g_z_0_y_0_yz_xxxxx, g_z_0_y_0_yz_xxxxxx, g_z_0_y_0_yz_xxxxxy, g_z_0_y_0_yz_xxxxxz, g_z_0_y_0_yz_xxxxy, g_z_0_y_0_yz_xxxxyy, g_z_0_y_0_yz_xxxxyz, g_z_0_y_0_yz_xxxxz, g_z_0_y_0_yz_xxxxzz, g_z_0_y_0_yz_xxxyy, g_z_0_y_0_yz_xxxyyy, g_z_0_y_0_yz_xxxyyz, g_z_0_y_0_yz_xxxyz, g_z_0_y_0_yz_xxxyzz, g_z_0_y_0_yz_xxxzz, g_z_0_y_0_yz_xxxzzz, g_z_0_y_0_yz_xxyyy, g_z_0_y_0_yz_xxyyyy, g_z_0_y_0_yz_xxyyyz, g_z_0_y_0_yz_xxyyz, g_z_0_y_0_yz_xxyyzz, g_z_0_y_0_yz_xxyzz, g_z_0_y_0_yz_xxyzzz, g_z_0_y_0_yz_xxzzz, g_z_0_y_0_yz_xxzzzz, g_z_0_y_0_yz_xyyyy, g_z_0_y_0_yz_xyyyyy, g_z_0_y_0_yz_xyyyyz, g_z_0_y_0_yz_xyyyz, g_z_0_y_0_yz_xyyyzz, g_z_0_y_0_yz_xyyzz, g_z_0_y_0_yz_xyyzzz, g_z_0_y_0_yz_xyzzz, g_z_0_y_0_yz_xyzzzz, g_z_0_y_0_yz_xzzzz, g_z_0_y_0_yz_xzzzzz, g_z_0_y_0_yz_yyyyy, g_z_0_y_0_yz_yyyyz, g_z_0_y_0_yz_yyyzz, g_z_0_y_0_yz_yyzzz, g_z_0_y_0_yz_yzzzz, g_z_0_y_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyz_xxxxx[k] = -g_z_0_y_0_yz_xxxxx[k] * ab_x + g_z_0_y_0_yz_xxxxxx[k];

                g_z_0_y_0_xyz_xxxxy[k] = -g_z_0_y_0_yz_xxxxy[k] * ab_x + g_z_0_y_0_yz_xxxxxy[k];

                g_z_0_y_0_xyz_xxxxz[k] = -g_z_0_y_0_yz_xxxxz[k] * ab_x + g_z_0_y_0_yz_xxxxxz[k];

                g_z_0_y_0_xyz_xxxyy[k] = -g_z_0_y_0_yz_xxxyy[k] * ab_x + g_z_0_y_0_yz_xxxxyy[k];

                g_z_0_y_0_xyz_xxxyz[k] = -g_z_0_y_0_yz_xxxyz[k] * ab_x + g_z_0_y_0_yz_xxxxyz[k];

                g_z_0_y_0_xyz_xxxzz[k] = -g_z_0_y_0_yz_xxxzz[k] * ab_x + g_z_0_y_0_yz_xxxxzz[k];

                g_z_0_y_0_xyz_xxyyy[k] = -g_z_0_y_0_yz_xxyyy[k] * ab_x + g_z_0_y_0_yz_xxxyyy[k];

                g_z_0_y_0_xyz_xxyyz[k] = -g_z_0_y_0_yz_xxyyz[k] * ab_x + g_z_0_y_0_yz_xxxyyz[k];

                g_z_0_y_0_xyz_xxyzz[k] = -g_z_0_y_0_yz_xxyzz[k] * ab_x + g_z_0_y_0_yz_xxxyzz[k];

                g_z_0_y_0_xyz_xxzzz[k] = -g_z_0_y_0_yz_xxzzz[k] * ab_x + g_z_0_y_0_yz_xxxzzz[k];

                g_z_0_y_0_xyz_xyyyy[k] = -g_z_0_y_0_yz_xyyyy[k] * ab_x + g_z_0_y_0_yz_xxyyyy[k];

                g_z_0_y_0_xyz_xyyyz[k] = -g_z_0_y_0_yz_xyyyz[k] * ab_x + g_z_0_y_0_yz_xxyyyz[k];

                g_z_0_y_0_xyz_xyyzz[k] = -g_z_0_y_0_yz_xyyzz[k] * ab_x + g_z_0_y_0_yz_xxyyzz[k];

                g_z_0_y_0_xyz_xyzzz[k] = -g_z_0_y_0_yz_xyzzz[k] * ab_x + g_z_0_y_0_yz_xxyzzz[k];

                g_z_0_y_0_xyz_xzzzz[k] = -g_z_0_y_0_yz_xzzzz[k] * ab_x + g_z_0_y_0_yz_xxzzzz[k];

                g_z_0_y_0_xyz_yyyyy[k] = -g_z_0_y_0_yz_yyyyy[k] * ab_x + g_z_0_y_0_yz_xyyyyy[k];

                g_z_0_y_0_xyz_yyyyz[k] = -g_z_0_y_0_yz_yyyyz[k] * ab_x + g_z_0_y_0_yz_xyyyyz[k];

                g_z_0_y_0_xyz_yyyzz[k] = -g_z_0_y_0_yz_yyyzz[k] * ab_x + g_z_0_y_0_yz_xyyyzz[k];

                g_z_0_y_0_xyz_yyzzz[k] = -g_z_0_y_0_yz_yyzzz[k] * ab_x + g_z_0_y_0_yz_xyyzzz[k];

                g_z_0_y_0_xyz_yzzzz[k] = -g_z_0_y_0_yz_yzzzz[k] * ab_x + g_z_0_y_0_yz_xyzzzz[k];

                g_z_0_y_0_xyz_zzzzz[k] = -g_z_0_y_0_yz_zzzzz[k] * ab_x + g_z_0_y_0_yz_xzzzzz[k];
            }

            /// Set up 1575-1596 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1575 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1576 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1577 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1578 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1579 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1580 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1581 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1582 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1583 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1584 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1585 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1586 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1587 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1588 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1589 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1590 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1591 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1592 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1593 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1594 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1595 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzz_xxxxx, g_z_0_y_0_xzz_xxxxy, g_z_0_y_0_xzz_xxxxz, g_z_0_y_0_xzz_xxxyy, g_z_0_y_0_xzz_xxxyz, g_z_0_y_0_xzz_xxxzz, g_z_0_y_0_xzz_xxyyy, g_z_0_y_0_xzz_xxyyz, g_z_0_y_0_xzz_xxyzz, g_z_0_y_0_xzz_xxzzz, g_z_0_y_0_xzz_xyyyy, g_z_0_y_0_xzz_xyyyz, g_z_0_y_0_xzz_xyyzz, g_z_0_y_0_xzz_xyzzz, g_z_0_y_0_xzz_xzzzz, g_z_0_y_0_xzz_yyyyy, g_z_0_y_0_xzz_yyyyz, g_z_0_y_0_xzz_yyyzz, g_z_0_y_0_xzz_yyzzz, g_z_0_y_0_xzz_yzzzz, g_z_0_y_0_xzz_zzzzz, g_z_0_y_0_zz_xxxxx, g_z_0_y_0_zz_xxxxxx, g_z_0_y_0_zz_xxxxxy, g_z_0_y_0_zz_xxxxxz, g_z_0_y_0_zz_xxxxy, g_z_0_y_0_zz_xxxxyy, g_z_0_y_0_zz_xxxxyz, g_z_0_y_0_zz_xxxxz, g_z_0_y_0_zz_xxxxzz, g_z_0_y_0_zz_xxxyy, g_z_0_y_0_zz_xxxyyy, g_z_0_y_0_zz_xxxyyz, g_z_0_y_0_zz_xxxyz, g_z_0_y_0_zz_xxxyzz, g_z_0_y_0_zz_xxxzz, g_z_0_y_0_zz_xxxzzz, g_z_0_y_0_zz_xxyyy, g_z_0_y_0_zz_xxyyyy, g_z_0_y_0_zz_xxyyyz, g_z_0_y_0_zz_xxyyz, g_z_0_y_0_zz_xxyyzz, g_z_0_y_0_zz_xxyzz, g_z_0_y_0_zz_xxyzzz, g_z_0_y_0_zz_xxzzz, g_z_0_y_0_zz_xxzzzz, g_z_0_y_0_zz_xyyyy, g_z_0_y_0_zz_xyyyyy, g_z_0_y_0_zz_xyyyyz, g_z_0_y_0_zz_xyyyz, g_z_0_y_0_zz_xyyyzz, g_z_0_y_0_zz_xyyzz, g_z_0_y_0_zz_xyyzzz, g_z_0_y_0_zz_xyzzz, g_z_0_y_0_zz_xyzzzz, g_z_0_y_0_zz_xzzzz, g_z_0_y_0_zz_xzzzzz, g_z_0_y_0_zz_yyyyy, g_z_0_y_0_zz_yyyyz, g_z_0_y_0_zz_yyyzz, g_z_0_y_0_zz_yyzzz, g_z_0_y_0_zz_yzzzz, g_z_0_y_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzz_xxxxx[k] = -g_z_0_y_0_zz_xxxxx[k] * ab_x + g_z_0_y_0_zz_xxxxxx[k];

                g_z_0_y_0_xzz_xxxxy[k] = -g_z_0_y_0_zz_xxxxy[k] * ab_x + g_z_0_y_0_zz_xxxxxy[k];

                g_z_0_y_0_xzz_xxxxz[k] = -g_z_0_y_0_zz_xxxxz[k] * ab_x + g_z_0_y_0_zz_xxxxxz[k];

                g_z_0_y_0_xzz_xxxyy[k] = -g_z_0_y_0_zz_xxxyy[k] * ab_x + g_z_0_y_0_zz_xxxxyy[k];

                g_z_0_y_0_xzz_xxxyz[k] = -g_z_0_y_0_zz_xxxyz[k] * ab_x + g_z_0_y_0_zz_xxxxyz[k];

                g_z_0_y_0_xzz_xxxzz[k] = -g_z_0_y_0_zz_xxxzz[k] * ab_x + g_z_0_y_0_zz_xxxxzz[k];

                g_z_0_y_0_xzz_xxyyy[k] = -g_z_0_y_0_zz_xxyyy[k] * ab_x + g_z_0_y_0_zz_xxxyyy[k];

                g_z_0_y_0_xzz_xxyyz[k] = -g_z_0_y_0_zz_xxyyz[k] * ab_x + g_z_0_y_0_zz_xxxyyz[k];

                g_z_0_y_0_xzz_xxyzz[k] = -g_z_0_y_0_zz_xxyzz[k] * ab_x + g_z_0_y_0_zz_xxxyzz[k];

                g_z_0_y_0_xzz_xxzzz[k] = -g_z_0_y_0_zz_xxzzz[k] * ab_x + g_z_0_y_0_zz_xxxzzz[k];

                g_z_0_y_0_xzz_xyyyy[k] = -g_z_0_y_0_zz_xyyyy[k] * ab_x + g_z_0_y_0_zz_xxyyyy[k];

                g_z_0_y_0_xzz_xyyyz[k] = -g_z_0_y_0_zz_xyyyz[k] * ab_x + g_z_0_y_0_zz_xxyyyz[k];

                g_z_0_y_0_xzz_xyyzz[k] = -g_z_0_y_0_zz_xyyzz[k] * ab_x + g_z_0_y_0_zz_xxyyzz[k];

                g_z_0_y_0_xzz_xyzzz[k] = -g_z_0_y_0_zz_xyzzz[k] * ab_x + g_z_0_y_0_zz_xxyzzz[k];

                g_z_0_y_0_xzz_xzzzz[k] = -g_z_0_y_0_zz_xzzzz[k] * ab_x + g_z_0_y_0_zz_xxzzzz[k];

                g_z_0_y_0_xzz_yyyyy[k] = -g_z_0_y_0_zz_yyyyy[k] * ab_x + g_z_0_y_0_zz_xyyyyy[k];

                g_z_0_y_0_xzz_yyyyz[k] = -g_z_0_y_0_zz_yyyyz[k] * ab_x + g_z_0_y_0_zz_xyyyyz[k];

                g_z_0_y_0_xzz_yyyzz[k] = -g_z_0_y_0_zz_yyyzz[k] * ab_x + g_z_0_y_0_zz_xyyyzz[k];

                g_z_0_y_0_xzz_yyzzz[k] = -g_z_0_y_0_zz_yyzzz[k] * ab_x + g_z_0_y_0_zz_xyyzzz[k];

                g_z_0_y_0_xzz_yzzzz[k] = -g_z_0_y_0_zz_yzzzz[k] * ab_x + g_z_0_y_0_zz_xyzzzz[k];

                g_z_0_y_0_xzz_zzzzz[k] = -g_z_0_y_0_zz_zzzzz[k] * ab_x + g_z_0_y_0_zz_xzzzzz[k];
            }

            /// Set up 1596-1617 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1596 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1597 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1598 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1599 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1600 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1601 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1602 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1603 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1604 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1605 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1606 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1607 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1608 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1609 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1610 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1611 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1612 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1613 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1614 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1615 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1616 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yy_xxxxx, g_z_0_y_0_yy_xxxxxy, g_z_0_y_0_yy_xxxxy, g_z_0_y_0_yy_xxxxyy, g_z_0_y_0_yy_xxxxyz, g_z_0_y_0_yy_xxxxz, g_z_0_y_0_yy_xxxyy, g_z_0_y_0_yy_xxxyyy, g_z_0_y_0_yy_xxxyyz, g_z_0_y_0_yy_xxxyz, g_z_0_y_0_yy_xxxyzz, g_z_0_y_0_yy_xxxzz, g_z_0_y_0_yy_xxyyy, g_z_0_y_0_yy_xxyyyy, g_z_0_y_0_yy_xxyyyz, g_z_0_y_0_yy_xxyyz, g_z_0_y_0_yy_xxyyzz, g_z_0_y_0_yy_xxyzz, g_z_0_y_0_yy_xxyzzz, g_z_0_y_0_yy_xxzzz, g_z_0_y_0_yy_xyyyy, g_z_0_y_0_yy_xyyyyy, g_z_0_y_0_yy_xyyyyz, g_z_0_y_0_yy_xyyyz, g_z_0_y_0_yy_xyyyzz, g_z_0_y_0_yy_xyyzz, g_z_0_y_0_yy_xyyzzz, g_z_0_y_0_yy_xyzzz, g_z_0_y_0_yy_xyzzzz, g_z_0_y_0_yy_xzzzz, g_z_0_y_0_yy_yyyyy, g_z_0_y_0_yy_yyyyyy, g_z_0_y_0_yy_yyyyyz, g_z_0_y_0_yy_yyyyz, g_z_0_y_0_yy_yyyyzz, g_z_0_y_0_yy_yyyzz, g_z_0_y_0_yy_yyyzzz, g_z_0_y_0_yy_yyzzz, g_z_0_y_0_yy_yyzzzz, g_z_0_y_0_yy_yzzzz, g_z_0_y_0_yy_yzzzzz, g_z_0_y_0_yy_zzzzz, g_z_0_y_0_yyy_xxxxx, g_z_0_y_0_yyy_xxxxy, g_z_0_y_0_yyy_xxxxz, g_z_0_y_0_yyy_xxxyy, g_z_0_y_0_yyy_xxxyz, g_z_0_y_0_yyy_xxxzz, g_z_0_y_0_yyy_xxyyy, g_z_0_y_0_yyy_xxyyz, g_z_0_y_0_yyy_xxyzz, g_z_0_y_0_yyy_xxzzz, g_z_0_y_0_yyy_xyyyy, g_z_0_y_0_yyy_xyyyz, g_z_0_y_0_yyy_xyyzz, g_z_0_y_0_yyy_xyzzz, g_z_0_y_0_yyy_xzzzz, g_z_0_y_0_yyy_yyyyy, g_z_0_y_0_yyy_yyyyz, g_z_0_y_0_yyy_yyyzz, g_z_0_y_0_yyy_yyzzz, g_z_0_y_0_yyy_yzzzz, g_z_0_y_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyy_xxxxx[k] = -g_z_0_y_0_yy_xxxxx[k] * ab_y + g_z_0_y_0_yy_xxxxxy[k];

                g_z_0_y_0_yyy_xxxxy[k] = -g_z_0_y_0_yy_xxxxy[k] * ab_y + g_z_0_y_0_yy_xxxxyy[k];

                g_z_0_y_0_yyy_xxxxz[k] = -g_z_0_y_0_yy_xxxxz[k] * ab_y + g_z_0_y_0_yy_xxxxyz[k];

                g_z_0_y_0_yyy_xxxyy[k] = -g_z_0_y_0_yy_xxxyy[k] * ab_y + g_z_0_y_0_yy_xxxyyy[k];

                g_z_0_y_0_yyy_xxxyz[k] = -g_z_0_y_0_yy_xxxyz[k] * ab_y + g_z_0_y_0_yy_xxxyyz[k];

                g_z_0_y_0_yyy_xxxzz[k] = -g_z_0_y_0_yy_xxxzz[k] * ab_y + g_z_0_y_0_yy_xxxyzz[k];

                g_z_0_y_0_yyy_xxyyy[k] = -g_z_0_y_0_yy_xxyyy[k] * ab_y + g_z_0_y_0_yy_xxyyyy[k];

                g_z_0_y_0_yyy_xxyyz[k] = -g_z_0_y_0_yy_xxyyz[k] * ab_y + g_z_0_y_0_yy_xxyyyz[k];

                g_z_0_y_0_yyy_xxyzz[k] = -g_z_0_y_0_yy_xxyzz[k] * ab_y + g_z_0_y_0_yy_xxyyzz[k];

                g_z_0_y_0_yyy_xxzzz[k] = -g_z_0_y_0_yy_xxzzz[k] * ab_y + g_z_0_y_0_yy_xxyzzz[k];

                g_z_0_y_0_yyy_xyyyy[k] = -g_z_0_y_0_yy_xyyyy[k] * ab_y + g_z_0_y_0_yy_xyyyyy[k];

                g_z_0_y_0_yyy_xyyyz[k] = -g_z_0_y_0_yy_xyyyz[k] * ab_y + g_z_0_y_0_yy_xyyyyz[k];

                g_z_0_y_0_yyy_xyyzz[k] = -g_z_0_y_0_yy_xyyzz[k] * ab_y + g_z_0_y_0_yy_xyyyzz[k];

                g_z_0_y_0_yyy_xyzzz[k] = -g_z_0_y_0_yy_xyzzz[k] * ab_y + g_z_0_y_0_yy_xyyzzz[k];

                g_z_0_y_0_yyy_xzzzz[k] = -g_z_0_y_0_yy_xzzzz[k] * ab_y + g_z_0_y_0_yy_xyzzzz[k];

                g_z_0_y_0_yyy_yyyyy[k] = -g_z_0_y_0_yy_yyyyy[k] * ab_y + g_z_0_y_0_yy_yyyyyy[k];

                g_z_0_y_0_yyy_yyyyz[k] = -g_z_0_y_0_yy_yyyyz[k] * ab_y + g_z_0_y_0_yy_yyyyyz[k];

                g_z_0_y_0_yyy_yyyzz[k] = -g_z_0_y_0_yy_yyyzz[k] * ab_y + g_z_0_y_0_yy_yyyyzz[k];

                g_z_0_y_0_yyy_yyzzz[k] = -g_z_0_y_0_yy_yyzzz[k] * ab_y + g_z_0_y_0_yy_yyyzzz[k];

                g_z_0_y_0_yyy_yzzzz[k] = -g_z_0_y_0_yy_yzzzz[k] * ab_y + g_z_0_y_0_yy_yyzzzz[k];

                g_z_0_y_0_yyy_zzzzz[k] = -g_z_0_y_0_yy_zzzzz[k] * ab_y + g_z_0_y_0_yy_yzzzzz[k];
            }

            /// Set up 1617-1638 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1617 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1618 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1619 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1620 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1621 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1622 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1623 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1624 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1625 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1626 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1627 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1628 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1629 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1630 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1631 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1632 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1633 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1634 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1635 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1636 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1637 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyz_xxxxx, g_z_0_y_0_yyz_xxxxy, g_z_0_y_0_yyz_xxxxz, g_z_0_y_0_yyz_xxxyy, g_z_0_y_0_yyz_xxxyz, g_z_0_y_0_yyz_xxxzz, g_z_0_y_0_yyz_xxyyy, g_z_0_y_0_yyz_xxyyz, g_z_0_y_0_yyz_xxyzz, g_z_0_y_0_yyz_xxzzz, g_z_0_y_0_yyz_xyyyy, g_z_0_y_0_yyz_xyyyz, g_z_0_y_0_yyz_xyyzz, g_z_0_y_0_yyz_xyzzz, g_z_0_y_0_yyz_xzzzz, g_z_0_y_0_yyz_yyyyy, g_z_0_y_0_yyz_yyyyz, g_z_0_y_0_yyz_yyyzz, g_z_0_y_0_yyz_yyzzz, g_z_0_y_0_yyz_yzzzz, g_z_0_y_0_yyz_zzzzz, g_z_0_y_0_yz_xxxxx, g_z_0_y_0_yz_xxxxxy, g_z_0_y_0_yz_xxxxy, g_z_0_y_0_yz_xxxxyy, g_z_0_y_0_yz_xxxxyz, g_z_0_y_0_yz_xxxxz, g_z_0_y_0_yz_xxxyy, g_z_0_y_0_yz_xxxyyy, g_z_0_y_0_yz_xxxyyz, g_z_0_y_0_yz_xxxyz, g_z_0_y_0_yz_xxxyzz, g_z_0_y_0_yz_xxxzz, g_z_0_y_0_yz_xxyyy, g_z_0_y_0_yz_xxyyyy, g_z_0_y_0_yz_xxyyyz, g_z_0_y_0_yz_xxyyz, g_z_0_y_0_yz_xxyyzz, g_z_0_y_0_yz_xxyzz, g_z_0_y_0_yz_xxyzzz, g_z_0_y_0_yz_xxzzz, g_z_0_y_0_yz_xyyyy, g_z_0_y_0_yz_xyyyyy, g_z_0_y_0_yz_xyyyyz, g_z_0_y_0_yz_xyyyz, g_z_0_y_0_yz_xyyyzz, g_z_0_y_0_yz_xyyzz, g_z_0_y_0_yz_xyyzzz, g_z_0_y_0_yz_xyzzz, g_z_0_y_0_yz_xyzzzz, g_z_0_y_0_yz_xzzzz, g_z_0_y_0_yz_yyyyy, g_z_0_y_0_yz_yyyyyy, g_z_0_y_0_yz_yyyyyz, g_z_0_y_0_yz_yyyyz, g_z_0_y_0_yz_yyyyzz, g_z_0_y_0_yz_yyyzz, g_z_0_y_0_yz_yyyzzz, g_z_0_y_0_yz_yyzzz, g_z_0_y_0_yz_yyzzzz, g_z_0_y_0_yz_yzzzz, g_z_0_y_0_yz_yzzzzz, g_z_0_y_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyz_xxxxx[k] = -g_z_0_y_0_yz_xxxxx[k] * ab_y + g_z_0_y_0_yz_xxxxxy[k];

                g_z_0_y_0_yyz_xxxxy[k] = -g_z_0_y_0_yz_xxxxy[k] * ab_y + g_z_0_y_0_yz_xxxxyy[k];

                g_z_0_y_0_yyz_xxxxz[k] = -g_z_0_y_0_yz_xxxxz[k] * ab_y + g_z_0_y_0_yz_xxxxyz[k];

                g_z_0_y_0_yyz_xxxyy[k] = -g_z_0_y_0_yz_xxxyy[k] * ab_y + g_z_0_y_0_yz_xxxyyy[k];

                g_z_0_y_0_yyz_xxxyz[k] = -g_z_0_y_0_yz_xxxyz[k] * ab_y + g_z_0_y_0_yz_xxxyyz[k];

                g_z_0_y_0_yyz_xxxzz[k] = -g_z_0_y_0_yz_xxxzz[k] * ab_y + g_z_0_y_0_yz_xxxyzz[k];

                g_z_0_y_0_yyz_xxyyy[k] = -g_z_0_y_0_yz_xxyyy[k] * ab_y + g_z_0_y_0_yz_xxyyyy[k];

                g_z_0_y_0_yyz_xxyyz[k] = -g_z_0_y_0_yz_xxyyz[k] * ab_y + g_z_0_y_0_yz_xxyyyz[k];

                g_z_0_y_0_yyz_xxyzz[k] = -g_z_0_y_0_yz_xxyzz[k] * ab_y + g_z_0_y_0_yz_xxyyzz[k];

                g_z_0_y_0_yyz_xxzzz[k] = -g_z_0_y_0_yz_xxzzz[k] * ab_y + g_z_0_y_0_yz_xxyzzz[k];

                g_z_0_y_0_yyz_xyyyy[k] = -g_z_0_y_0_yz_xyyyy[k] * ab_y + g_z_0_y_0_yz_xyyyyy[k];

                g_z_0_y_0_yyz_xyyyz[k] = -g_z_0_y_0_yz_xyyyz[k] * ab_y + g_z_0_y_0_yz_xyyyyz[k];

                g_z_0_y_0_yyz_xyyzz[k] = -g_z_0_y_0_yz_xyyzz[k] * ab_y + g_z_0_y_0_yz_xyyyzz[k];

                g_z_0_y_0_yyz_xyzzz[k] = -g_z_0_y_0_yz_xyzzz[k] * ab_y + g_z_0_y_0_yz_xyyzzz[k];

                g_z_0_y_0_yyz_xzzzz[k] = -g_z_0_y_0_yz_xzzzz[k] * ab_y + g_z_0_y_0_yz_xyzzzz[k];

                g_z_0_y_0_yyz_yyyyy[k] = -g_z_0_y_0_yz_yyyyy[k] * ab_y + g_z_0_y_0_yz_yyyyyy[k];

                g_z_0_y_0_yyz_yyyyz[k] = -g_z_0_y_0_yz_yyyyz[k] * ab_y + g_z_0_y_0_yz_yyyyyz[k];

                g_z_0_y_0_yyz_yyyzz[k] = -g_z_0_y_0_yz_yyyzz[k] * ab_y + g_z_0_y_0_yz_yyyyzz[k];

                g_z_0_y_0_yyz_yyzzz[k] = -g_z_0_y_0_yz_yyzzz[k] * ab_y + g_z_0_y_0_yz_yyyzzz[k];

                g_z_0_y_0_yyz_yzzzz[k] = -g_z_0_y_0_yz_yzzzz[k] * ab_y + g_z_0_y_0_yz_yyzzzz[k];

                g_z_0_y_0_yyz_zzzzz[k] = -g_z_0_y_0_yz_zzzzz[k] * ab_y + g_z_0_y_0_yz_yzzzzz[k];
            }

            /// Set up 1638-1659 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1638 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1639 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1640 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1641 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1642 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1643 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1644 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1645 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1646 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1647 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1648 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1649 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1650 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1651 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1652 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1653 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1654 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1655 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1656 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1657 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1658 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzz_xxxxx, g_z_0_y_0_yzz_xxxxy, g_z_0_y_0_yzz_xxxxz, g_z_0_y_0_yzz_xxxyy, g_z_0_y_0_yzz_xxxyz, g_z_0_y_0_yzz_xxxzz, g_z_0_y_0_yzz_xxyyy, g_z_0_y_0_yzz_xxyyz, g_z_0_y_0_yzz_xxyzz, g_z_0_y_0_yzz_xxzzz, g_z_0_y_0_yzz_xyyyy, g_z_0_y_0_yzz_xyyyz, g_z_0_y_0_yzz_xyyzz, g_z_0_y_0_yzz_xyzzz, g_z_0_y_0_yzz_xzzzz, g_z_0_y_0_yzz_yyyyy, g_z_0_y_0_yzz_yyyyz, g_z_0_y_0_yzz_yyyzz, g_z_0_y_0_yzz_yyzzz, g_z_0_y_0_yzz_yzzzz, g_z_0_y_0_yzz_zzzzz, g_z_0_y_0_zz_xxxxx, g_z_0_y_0_zz_xxxxxy, g_z_0_y_0_zz_xxxxy, g_z_0_y_0_zz_xxxxyy, g_z_0_y_0_zz_xxxxyz, g_z_0_y_0_zz_xxxxz, g_z_0_y_0_zz_xxxyy, g_z_0_y_0_zz_xxxyyy, g_z_0_y_0_zz_xxxyyz, g_z_0_y_0_zz_xxxyz, g_z_0_y_0_zz_xxxyzz, g_z_0_y_0_zz_xxxzz, g_z_0_y_0_zz_xxyyy, g_z_0_y_0_zz_xxyyyy, g_z_0_y_0_zz_xxyyyz, g_z_0_y_0_zz_xxyyz, g_z_0_y_0_zz_xxyyzz, g_z_0_y_0_zz_xxyzz, g_z_0_y_0_zz_xxyzzz, g_z_0_y_0_zz_xxzzz, g_z_0_y_0_zz_xyyyy, g_z_0_y_0_zz_xyyyyy, g_z_0_y_0_zz_xyyyyz, g_z_0_y_0_zz_xyyyz, g_z_0_y_0_zz_xyyyzz, g_z_0_y_0_zz_xyyzz, g_z_0_y_0_zz_xyyzzz, g_z_0_y_0_zz_xyzzz, g_z_0_y_0_zz_xyzzzz, g_z_0_y_0_zz_xzzzz, g_z_0_y_0_zz_yyyyy, g_z_0_y_0_zz_yyyyyy, g_z_0_y_0_zz_yyyyyz, g_z_0_y_0_zz_yyyyz, g_z_0_y_0_zz_yyyyzz, g_z_0_y_0_zz_yyyzz, g_z_0_y_0_zz_yyyzzz, g_z_0_y_0_zz_yyzzz, g_z_0_y_0_zz_yyzzzz, g_z_0_y_0_zz_yzzzz, g_z_0_y_0_zz_yzzzzz, g_z_0_y_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzz_xxxxx[k] = -g_z_0_y_0_zz_xxxxx[k] * ab_y + g_z_0_y_0_zz_xxxxxy[k];

                g_z_0_y_0_yzz_xxxxy[k] = -g_z_0_y_0_zz_xxxxy[k] * ab_y + g_z_0_y_0_zz_xxxxyy[k];

                g_z_0_y_0_yzz_xxxxz[k] = -g_z_0_y_0_zz_xxxxz[k] * ab_y + g_z_0_y_0_zz_xxxxyz[k];

                g_z_0_y_0_yzz_xxxyy[k] = -g_z_0_y_0_zz_xxxyy[k] * ab_y + g_z_0_y_0_zz_xxxyyy[k];

                g_z_0_y_0_yzz_xxxyz[k] = -g_z_0_y_0_zz_xxxyz[k] * ab_y + g_z_0_y_0_zz_xxxyyz[k];

                g_z_0_y_0_yzz_xxxzz[k] = -g_z_0_y_0_zz_xxxzz[k] * ab_y + g_z_0_y_0_zz_xxxyzz[k];

                g_z_0_y_0_yzz_xxyyy[k] = -g_z_0_y_0_zz_xxyyy[k] * ab_y + g_z_0_y_0_zz_xxyyyy[k];

                g_z_0_y_0_yzz_xxyyz[k] = -g_z_0_y_0_zz_xxyyz[k] * ab_y + g_z_0_y_0_zz_xxyyyz[k];

                g_z_0_y_0_yzz_xxyzz[k] = -g_z_0_y_0_zz_xxyzz[k] * ab_y + g_z_0_y_0_zz_xxyyzz[k];

                g_z_0_y_0_yzz_xxzzz[k] = -g_z_0_y_0_zz_xxzzz[k] * ab_y + g_z_0_y_0_zz_xxyzzz[k];

                g_z_0_y_0_yzz_xyyyy[k] = -g_z_0_y_0_zz_xyyyy[k] * ab_y + g_z_0_y_0_zz_xyyyyy[k];

                g_z_0_y_0_yzz_xyyyz[k] = -g_z_0_y_0_zz_xyyyz[k] * ab_y + g_z_0_y_0_zz_xyyyyz[k];

                g_z_0_y_0_yzz_xyyzz[k] = -g_z_0_y_0_zz_xyyzz[k] * ab_y + g_z_0_y_0_zz_xyyyzz[k];

                g_z_0_y_0_yzz_xyzzz[k] = -g_z_0_y_0_zz_xyzzz[k] * ab_y + g_z_0_y_0_zz_xyyzzz[k];

                g_z_0_y_0_yzz_xzzzz[k] = -g_z_0_y_0_zz_xzzzz[k] * ab_y + g_z_0_y_0_zz_xyzzzz[k];

                g_z_0_y_0_yzz_yyyyy[k] = -g_z_0_y_0_zz_yyyyy[k] * ab_y + g_z_0_y_0_zz_yyyyyy[k];

                g_z_0_y_0_yzz_yyyyz[k] = -g_z_0_y_0_zz_yyyyz[k] * ab_y + g_z_0_y_0_zz_yyyyyz[k];

                g_z_0_y_0_yzz_yyyzz[k] = -g_z_0_y_0_zz_yyyzz[k] * ab_y + g_z_0_y_0_zz_yyyyzz[k];

                g_z_0_y_0_yzz_yyzzz[k] = -g_z_0_y_0_zz_yyzzz[k] * ab_y + g_z_0_y_0_zz_yyyzzz[k];

                g_z_0_y_0_yzz_yzzzz[k] = -g_z_0_y_0_zz_yzzzz[k] * ab_y + g_z_0_y_0_zz_yyzzzz[k];

                g_z_0_y_0_yzz_zzzzz[k] = -g_z_0_y_0_zz_zzzzz[k] * ab_y + g_z_0_y_0_zz_yzzzzz[k];
            }

            /// Set up 1659-1680 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1659 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1660 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1661 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1662 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1663 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1664 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1665 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1666 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1667 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1668 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1669 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1670 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1671 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1672 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1673 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1674 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1675 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1676 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1677 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1678 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zz_xxxxx, g_0_0_y_0_zz_xxxxy, g_0_0_y_0_zz_xxxxz, g_0_0_y_0_zz_xxxyy, g_0_0_y_0_zz_xxxyz, g_0_0_y_0_zz_xxxzz, g_0_0_y_0_zz_xxyyy, g_0_0_y_0_zz_xxyyz, g_0_0_y_0_zz_xxyzz, g_0_0_y_0_zz_xxzzz, g_0_0_y_0_zz_xyyyy, g_0_0_y_0_zz_xyyyz, g_0_0_y_0_zz_xyyzz, g_0_0_y_0_zz_xyzzz, g_0_0_y_0_zz_xzzzz, g_0_0_y_0_zz_yyyyy, g_0_0_y_0_zz_yyyyz, g_0_0_y_0_zz_yyyzz, g_0_0_y_0_zz_yyzzz, g_0_0_y_0_zz_yzzzz, g_0_0_y_0_zz_zzzzz, g_z_0_y_0_zz_xxxxx, g_z_0_y_0_zz_xxxxxz, g_z_0_y_0_zz_xxxxy, g_z_0_y_0_zz_xxxxyz, g_z_0_y_0_zz_xxxxz, g_z_0_y_0_zz_xxxxzz, g_z_0_y_0_zz_xxxyy, g_z_0_y_0_zz_xxxyyz, g_z_0_y_0_zz_xxxyz, g_z_0_y_0_zz_xxxyzz, g_z_0_y_0_zz_xxxzz, g_z_0_y_0_zz_xxxzzz, g_z_0_y_0_zz_xxyyy, g_z_0_y_0_zz_xxyyyz, g_z_0_y_0_zz_xxyyz, g_z_0_y_0_zz_xxyyzz, g_z_0_y_0_zz_xxyzz, g_z_0_y_0_zz_xxyzzz, g_z_0_y_0_zz_xxzzz, g_z_0_y_0_zz_xxzzzz, g_z_0_y_0_zz_xyyyy, g_z_0_y_0_zz_xyyyyz, g_z_0_y_0_zz_xyyyz, g_z_0_y_0_zz_xyyyzz, g_z_0_y_0_zz_xyyzz, g_z_0_y_0_zz_xyyzzz, g_z_0_y_0_zz_xyzzz, g_z_0_y_0_zz_xyzzzz, g_z_0_y_0_zz_xzzzz, g_z_0_y_0_zz_xzzzzz, g_z_0_y_0_zz_yyyyy, g_z_0_y_0_zz_yyyyyz, g_z_0_y_0_zz_yyyyz, g_z_0_y_0_zz_yyyyzz, g_z_0_y_0_zz_yyyzz, g_z_0_y_0_zz_yyyzzz, g_z_0_y_0_zz_yyzzz, g_z_0_y_0_zz_yyzzzz, g_z_0_y_0_zz_yzzzz, g_z_0_y_0_zz_yzzzzz, g_z_0_y_0_zz_zzzzz, g_z_0_y_0_zz_zzzzzz, g_z_0_y_0_zzz_xxxxx, g_z_0_y_0_zzz_xxxxy, g_z_0_y_0_zzz_xxxxz, g_z_0_y_0_zzz_xxxyy, g_z_0_y_0_zzz_xxxyz, g_z_0_y_0_zzz_xxxzz, g_z_0_y_0_zzz_xxyyy, g_z_0_y_0_zzz_xxyyz, g_z_0_y_0_zzz_xxyzz, g_z_0_y_0_zzz_xxzzz, g_z_0_y_0_zzz_xyyyy, g_z_0_y_0_zzz_xyyyz, g_z_0_y_0_zzz_xyyzz, g_z_0_y_0_zzz_xyzzz, g_z_0_y_0_zzz_xzzzz, g_z_0_y_0_zzz_yyyyy, g_z_0_y_0_zzz_yyyyz, g_z_0_y_0_zzz_yyyzz, g_z_0_y_0_zzz_yyzzz, g_z_0_y_0_zzz_yzzzz, g_z_0_y_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzz_xxxxx[k] = -g_0_0_y_0_zz_xxxxx[k] - g_z_0_y_0_zz_xxxxx[k] * ab_z + g_z_0_y_0_zz_xxxxxz[k];

                g_z_0_y_0_zzz_xxxxy[k] = -g_0_0_y_0_zz_xxxxy[k] - g_z_0_y_0_zz_xxxxy[k] * ab_z + g_z_0_y_0_zz_xxxxyz[k];

                g_z_0_y_0_zzz_xxxxz[k] = -g_0_0_y_0_zz_xxxxz[k] - g_z_0_y_0_zz_xxxxz[k] * ab_z + g_z_0_y_0_zz_xxxxzz[k];

                g_z_0_y_0_zzz_xxxyy[k] = -g_0_0_y_0_zz_xxxyy[k] - g_z_0_y_0_zz_xxxyy[k] * ab_z + g_z_0_y_0_zz_xxxyyz[k];

                g_z_0_y_0_zzz_xxxyz[k] = -g_0_0_y_0_zz_xxxyz[k] - g_z_0_y_0_zz_xxxyz[k] * ab_z + g_z_0_y_0_zz_xxxyzz[k];

                g_z_0_y_0_zzz_xxxzz[k] = -g_0_0_y_0_zz_xxxzz[k] - g_z_0_y_0_zz_xxxzz[k] * ab_z + g_z_0_y_0_zz_xxxzzz[k];

                g_z_0_y_0_zzz_xxyyy[k] = -g_0_0_y_0_zz_xxyyy[k] - g_z_0_y_0_zz_xxyyy[k] * ab_z + g_z_0_y_0_zz_xxyyyz[k];

                g_z_0_y_0_zzz_xxyyz[k] = -g_0_0_y_0_zz_xxyyz[k] - g_z_0_y_0_zz_xxyyz[k] * ab_z + g_z_0_y_0_zz_xxyyzz[k];

                g_z_0_y_0_zzz_xxyzz[k] = -g_0_0_y_0_zz_xxyzz[k] - g_z_0_y_0_zz_xxyzz[k] * ab_z + g_z_0_y_0_zz_xxyzzz[k];

                g_z_0_y_0_zzz_xxzzz[k] = -g_0_0_y_0_zz_xxzzz[k] - g_z_0_y_0_zz_xxzzz[k] * ab_z + g_z_0_y_0_zz_xxzzzz[k];

                g_z_0_y_0_zzz_xyyyy[k] = -g_0_0_y_0_zz_xyyyy[k] - g_z_0_y_0_zz_xyyyy[k] * ab_z + g_z_0_y_0_zz_xyyyyz[k];

                g_z_0_y_0_zzz_xyyyz[k] = -g_0_0_y_0_zz_xyyyz[k] - g_z_0_y_0_zz_xyyyz[k] * ab_z + g_z_0_y_0_zz_xyyyzz[k];

                g_z_0_y_0_zzz_xyyzz[k] = -g_0_0_y_0_zz_xyyzz[k] - g_z_0_y_0_zz_xyyzz[k] * ab_z + g_z_0_y_0_zz_xyyzzz[k];

                g_z_0_y_0_zzz_xyzzz[k] = -g_0_0_y_0_zz_xyzzz[k] - g_z_0_y_0_zz_xyzzz[k] * ab_z + g_z_0_y_0_zz_xyzzzz[k];

                g_z_0_y_0_zzz_xzzzz[k] = -g_0_0_y_0_zz_xzzzz[k] - g_z_0_y_0_zz_xzzzz[k] * ab_z + g_z_0_y_0_zz_xzzzzz[k];

                g_z_0_y_0_zzz_yyyyy[k] = -g_0_0_y_0_zz_yyyyy[k] - g_z_0_y_0_zz_yyyyy[k] * ab_z + g_z_0_y_0_zz_yyyyyz[k];

                g_z_0_y_0_zzz_yyyyz[k] = -g_0_0_y_0_zz_yyyyz[k] - g_z_0_y_0_zz_yyyyz[k] * ab_z + g_z_0_y_0_zz_yyyyzz[k];

                g_z_0_y_0_zzz_yyyzz[k] = -g_0_0_y_0_zz_yyyzz[k] - g_z_0_y_0_zz_yyyzz[k] * ab_z + g_z_0_y_0_zz_yyyzzz[k];

                g_z_0_y_0_zzz_yyzzz[k] = -g_0_0_y_0_zz_yyzzz[k] - g_z_0_y_0_zz_yyzzz[k] * ab_z + g_z_0_y_0_zz_yyzzzz[k];

                g_z_0_y_0_zzz_yzzzz[k] = -g_0_0_y_0_zz_yzzzz[k] - g_z_0_y_0_zz_yzzzz[k] * ab_z + g_z_0_y_0_zz_yzzzzz[k];

                g_z_0_y_0_zzz_zzzzz[k] = -g_0_0_y_0_zz_zzzzz[k] - g_z_0_y_0_zz_zzzzz[k] * ab_z + g_z_0_y_0_zz_zzzzzz[k];
            }

            /// Set up 1680-1701 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxx_xxxxx = cbuffer.data(fh_geom_1010_off + 1680 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxxy = cbuffer.data(fh_geom_1010_off + 1681 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxxz = cbuffer.data(fh_geom_1010_off + 1682 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxyy = cbuffer.data(fh_geom_1010_off + 1683 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxyz = cbuffer.data(fh_geom_1010_off + 1684 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxxzz = cbuffer.data(fh_geom_1010_off + 1685 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxyyy = cbuffer.data(fh_geom_1010_off + 1686 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxyyz = cbuffer.data(fh_geom_1010_off + 1687 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxyzz = cbuffer.data(fh_geom_1010_off + 1688 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xxzzz = cbuffer.data(fh_geom_1010_off + 1689 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyyyy = cbuffer.data(fh_geom_1010_off + 1690 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyyyz = cbuffer.data(fh_geom_1010_off + 1691 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyyzz = cbuffer.data(fh_geom_1010_off + 1692 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xyzzz = cbuffer.data(fh_geom_1010_off + 1693 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xzzzz = cbuffer.data(fh_geom_1010_off + 1694 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyyyy = cbuffer.data(fh_geom_1010_off + 1695 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyyyz = cbuffer.data(fh_geom_1010_off + 1696 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyyzz = cbuffer.data(fh_geom_1010_off + 1697 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yyzzz = cbuffer.data(fh_geom_1010_off + 1698 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yzzzz = cbuffer.data(fh_geom_1010_off + 1699 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_zzzzz = cbuffer.data(fh_geom_1010_off + 1700 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xx_xxxxx, g_z_0_z_0_xx_xxxxxx, g_z_0_z_0_xx_xxxxxy, g_z_0_z_0_xx_xxxxxz, g_z_0_z_0_xx_xxxxy, g_z_0_z_0_xx_xxxxyy, g_z_0_z_0_xx_xxxxyz, g_z_0_z_0_xx_xxxxz, g_z_0_z_0_xx_xxxxzz, g_z_0_z_0_xx_xxxyy, g_z_0_z_0_xx_xxxyyy, g_z_0_z_0_xx_xxxyyz, g_z_0_z_0_xx_xxxyz, g_z_0_z_0_xx_xxxyzz, g_z_0_z_0_xx_xxxzz, g_z_0_z_0_xx_xxxzzz, g_z_0_z_0_xx_xxyyy, g_z_0_z_0_xx_xxyyyy, g_z_0_z_0_xx_xxyyyz, g_z_0_z_0_xx_xxyyz, g_z_0_z_0_xx_xxyyzz, g_z_0_z_0_xx_xxyzz, g_z_0_z_0_xx_xxyzzz, g_z_0_z_0_xx_xxzzz, g_z_0_z_0_xx_xxzzzz, g_z_0_z_0_xx_xyyyy, g_z_0_z_0_xx_xyyyyy, g_z_0_z_0_xx_xyyyyz, g_z_0_z_0_xx_xyyyz, g_z_0_z_0_xx_xyyyzz, g_z_0_z_0_xx_xyyzz, g_z_0_z_0_xx_xyyzzz, g_z_0_z_0_xx_xyzzz, g_z_0_z_0_xx_xyzzzz, g_z_0_z_0_xx_xzzzz, g_z_0_z_0_xx_xzzzzz, g_z_0_z_0_xx_yyyyy, g_z_0_z_0_xx_yyyyz, g_z_0_z_0_xx_yyyzz, g_z_0_z_0_xx_yyzzz, g_z_0_z_0_xx_yzzzz, g_z_0_z_0_xx_zzzzz, g_z_0_z_0_xxx_xxxxx, g_z_0_z_0_xxx_xxxxy, g_z_0_z_0_xxx_xxxxz, g_z_0_z_0_xxx_xxxyy, g_z_0_z_0_xxx_xxxyz, g_z_0_z_0_xxx_xxxzz, g_z_0_z_0_xxx_xxyyy, g_z_0_z_0_xxx_xxyyz, g_z_0_z_0_xxx_xxyzz, g_z_0_z_0_xxx_xxzzz, g_z_0_z_0_xxx_xyyyy, g_z_0_z_0_xxx_xyyyz, g_z_0_z_0_xxx_xyyzz, g_z_0_z_0_xxx_xyzzz, g_z_0_z_0_xxx_xzzzz, g_z_0_z_0_xxx_yyyyy, g_z_0_z_0_xxx_yyyyz, g_z_0_z_0_xxx_yyyzz, g_z_0_z_0_xxx_yyzzz, g_z_0_z_0_xxx_yzzzz, g_z_0_z_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxx_xxxxx[k] = -g_z_0_z_0_xx_xxxxx[k] * ab_x + g_z_0_z_0_xx_xxxxxx[k];

                g_z_0_z_0_xxx_xxxxy[k] = -g_z_0_z_0_xx_xxxxy[k] * ab_x + g_z_0_z_0_xx_xxxxxy[k];

                g_z_0_z_0_xxx_xxxxz[k] = -g_z_0_z_0_xx_xxxxz[k] * ab_x + g_z_0_z_0_xx_xxxxxz[k];

                g_z_0_z_0_xxx_xxxyy[k] = -g_z_0_z_0_xx_xxxyy[k] * ab_x + g_z_0_z_0_xx_xxxxyy[k];

                g_z_0_z_0_xxx_xxxyz[k] = -g_z_0_z_0_xx_xxxyz[k] * ab_x + g_z_0_z_0_xx_xxxxyz[k];

                g_z_0_z_0_xxx_xxxzz[k] = -g_z_0_z_0_xx_xxxzz[k] * ab_x + g_z_0_z_0_xx_xxxxzz[k];

                g_z_0_z_0_xxx_xxyyy[k] = -g_z_0_z_0_xx_xxyyy[k] * ab_x + g_z_0_z_0_xx_xxxyyy[k];

                g_z_0_z_0_xxx_xxyyz[k] = -g_z_0_z_0_xx_xxyyz[k] * ab_x + g_z_0_z_0_xx_xxxyyz[k];

                g_z_0_z_0_xxx_xxyzz[k] = -g_z_0_z_0_xx_xxyzz[k] * ab_x + g_z_0_z_0_xx_xxxyzz[k];

                g_z_0_z_0_xxx_xxzzz[k] = -g_z_0_z_0_xx_xxzzz[k] * ab_x + g_z_0_z_0_xx_xxxzzz[k];

                g_z_0_z_0_xxx_xyyyy[k] = -g_z_0_z_0_xx_xyyyy[k] * ab_x + g_z_0_z_0_xx_xxyyyy[k];

                g_z_0_z_0_xxx_xyyyz[k] = -g_z_0_z_0_xx_xyyyz[k] * ab_x + g_z_0_z_0_xx_xxyyyz[k];

                g_z_0_z_0_xxx_xyyzz[k] = -g_z_0_z_0_xx_xyyzz[k] * ab_x + g_z_0_z_0_xx_xxyyzz[k];

                g_z_0_z_0_xxx_xyzzz[k] = -g_z_0_z_0_xx_xyzzz[k] * ab_x + g_z_0_z_0_xx_xxyzzz[k];

                g_z_0_z_0_xxx_xzzzz[k] = -g_z_0_z_0_xx_xzzzz[k] * ab_x + g_z_0_z_0_xx_xxzzzz[k];

                g_z_0_z_0_xxx_yyyyy[k] = -g_z_0_z_0_xx_yyyyy[k] * ab_x + g_z_0_z_0_xx_xyyyyy[k];

                g_z_0_z_0_xxx_yyyyz[k] = -g_z_0_z_0_xx_yyyyz[k] * ab_x + g_z_0_z_0_xx_xyyyyz[k];

                g_z_0_z_0_xxx_yyyzz[k] = -g_z_0_z_0_xx_yyyzz[k] * ab_x + g_z_0_z_0_xx_xyyyzz[k];

                g_z_0_z_0_xxx_yyzzz[k] = -g_z_0_z_0_xx_yyzzz[k] * ab_x + g_z_0_z_0_xx_xyyzzz[k];

                g_z_0_z_0_xxx_yzzzz[k] = -g_z_0_z_0_xx_yzzzz[k] * ab_x + g_z_0_z_0_xx_xyzzzz[k];

                g_z_0_z_0_xxx_zzzzz[k] = -g_z_0_z_0_xx_zzzzz[k] * ab_x + g_z_0_z_0_xx_xzzzzz[k];
            }

            /// Set up 1701-1722 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxy_xxxxx = cbuffer.data(fh_geom_1010_off + 1701 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxxy = cbuffer.data(fh_geom_1010_off + 1702 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxxz = cbuffer.data(fh_geom_1010_off + 1703 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxyy = cbuffer.data(fh_geom_1010_off + 1704 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxyz = cbuffer.data(fh_geom_1010_off + 1705 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxxzz = cbuffer.data(fh_geom_1010_off + 1706 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxyyy = cbuffer.data(fh_geom_1010_off + 1707 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxyyz = cbuffer.data(fh_geom_1010_off + 1708 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxyzz = cbuffer.data(fh_geom_1010_off + 1709 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xxzzz = cbuffer.data(fh_geom_1010_off + 1710 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyyyy = cbuffer.data(fh_geom_1010_off + 1711 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyyyz = cbuffer.data(fh_geom_1010_off + 1712 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyyzz = cbuffer.data(fh_geom_1010_off + 1713 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xyzzz = cbuffer.data(fh_geom_1010_off + 1714 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xzzzz = cbuffer.data(fh_geom_1010_off + 1715 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyyyy = cbuffer.data(fh_geom_1010_off + 1716 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyyyz = cbuffer.data(fh_geom_1010_off + 1717 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyyzz = cbuffer.data(fh_geom_1010_off + 1718 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yyzzz = cbuffer.data(fh_geom_1010_off + 1719 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yzzzz = cbuffer.data(fh_geom_1010_off + 1720 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_zzzzz = cbuffer.data(fh_geom_1010_off + 1721 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxy_xxxxx, g_z_0_z_0_xxy_xxxxy, g_z_0_z_0_xxy_xxxxz, g_z_0_z_0_xxy_xxxyy, g_z_0_z_0_xxy_xxxyz, g_z_0_z_0_xxy_xxxzz, g_z_0_z_0_xxy_xxyyy, g_z_0_z_0_xxy_xxyyz, g_z_0_z_0_xxy_xxyzz, g_z_0_z_0_xxy_xxzzz, g_z_0_z_0_xxy_xyyyy, g_z_0_z_0_xxy_xyyyz, g_z_0_z_0_xxy_xyyzz, g_z_0_z_0_xxy_xyzzz, g_z_0_z_0_xxy_xzzzz, g_z_0_z_0_xxy_yyyyy, g_z_0_z_0_xxy_yyyyz, g_z_0_z_0_xxy_yyyzz, g_z_0_z_0_xxy_yyzzz, g_z_0_z_0_xxy_yzzzz, g_z_0_z_0_xxy_zzzzz, g_z_0_z_0_xy_xxxxx, g_z_0_z_0_xy_xxxxxx, g_z_0_z_0_xy_xxxxxy, g_z_0_z_0_xy_xxxxxz, g_z_0_z_0_xy_xxxxy, g_z_0_z_0_xy_xxxxyy, g_z_0_z_0_xy_xxxxyz, g_z_0_z_0_xy_xxxxz, g_z_0_z_0_xy_xxxxzz, g_z_0_z_0_xy_xxxyy, g_z_0_z_0_xy_xxxyyy, g_z_0_z_0_xy_xxxyyz, g_z_0_z_0_xy_xxxyz, g_z_0_z_0_xy_xxxyzz, g_z_0_z_0_xy_xxxzz, g_z_0_z_0_xy_xxxzzz, g_z_0_z_0_xy_xxyyy, g_z_0_z_0_xy_xxyyyy, g_z_0_z_0_xy_xxyyyz, g_z_0_z_0_xy_xxyyz, g_z_0_z_0_xy_xxyyzz, g_z_0_z_0_xy_xxyzz, g_z_0_z_0_xy_xxyzzz, g_z_0_z_0_xy_xxzzz, g_z_0_z_0_xy_xxzzzz, g_z_0_z_0_xy_xyyyy, g_z_0_z_0_xy_xyyyyy, g_z_0_z_0_xy_xyyyyz, g_z_0_z_0_xy_xyyyz, g_z_0_z_0_xy_xyyyzz, g_z_0_z_0_xy_xyyzz, g_z_0_z_0_xy_xyyzzz, g_z_0_z_0_xy_xyzzz, g_z_0_z_0_xy_xyzzzz, g_z_0_z_0_xy_xzzzz, g_z_0_z_0_xy_xzzzzz, g_z_0_z_0_xy_yyyyy, g_z_0_z_0_xy_yyyyz, g_z_0_z_0_xy_yyyzz, g_z_0_z_0_xy_yyzzz, g_z_0_z_0_xy_yzzzz, g_z_0_z_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxy_xxxxx[k] = -g_z_0_z_0_xy_xxxxx[k] * ab_x + g_z_0_z_0_xy_xxxxxx[k];

                g_z_0_z_0_xxy_xxxxy[k] = -g_z_0_z_0_xy_xxxxy[k] * ab_x + g_z_0_z_0_xy_xxxxxy[k];

                g_z_0_z_0_xxy_xxxxz[k] = -g_z_0_z_0_xy_xxxxz[k] * ab_x + g_z_0_z_0_xy_xxxxxz[k];

                g_z_0_z_0_xxy_xxxyy[k] = -g_z_0_z_0_xy_xxxyy[k] * ab_x + g_z_0_z_0_xy_xxxxyy[k];

                g_z_0_z_0_xxy_xxxyz[k] = -g_z_0_z_0_xy_xxxyz[k] * ab_x + g_z_0_z_0_xy_xxxxyz[k];

                g_z_0_z_0_xxy_xxxzz[k] = -g_z_0_z_0_xy_xxxzz[k] * ab_x + g_z_0_z_0_xy_xxxxzz[k];

                g_z_0_z_0_xxy_xxyyy[k] = -g_z_0_z_0_xy_xxyyy[k] * ab_x + g_z_0_z_0_xy_xxxyyy[k];

                g_z_0_z_0_xxy_xxyyz[k] = -g_z_0_z_0_xy_xxyyz[k] * ab_x + g_z_0_z_0_xy_xxxyyz[k];

                g_z_0_z_0_xxy_xxyzz[k] = -g_z_0_z_0_xy_xxyzz[k] * ab_x + g_z_0_z_0_xy_xxxyzz[k];

                g_z_0_z_0_xxy_xxzzz[k] = -g_z_0_z_0_xy_xxzzz[k] * ab_x + g_z_0_z_0_xy_xxxzzz[k];

                g_z_0_z_0_xxy_xyyyy[k] = -g_z_0_z_0_xy_xyyyy[k] * ab_x + g_z_0_z_0_xy_xxyyyy[k];

                g_z_0_z_0_xxy_xyyyz[k] = -g_z_0_z_0_xy_xyyyz[k] * ab_x + g_z_0_z_0_xy_xxyyyz[k];

                g_z_0_z_0_xxy_xyyzz[k] = -g_z_0_z_0_xy_xyyzz[k] * ab_x + g_z_0_z_0_xy_xxyyzz[k];

                g_z_0_z_0_xxy_xyzzz[k] = -g_z_0_z_0_xy_xyzzz[k] * ab_x + g_z_0_z_0_xy_xxyzzz[k];

                g_z_0_z_0_xxy_xzzzz[k] = -g_z_0_z_0_xy_xzzzz[k] * ab_x + g_z_0_z_0_xy_xxzzzz[k];

                g_z_0_z_0_xxy_yyyyy[k] = -g_z_0_z_0_xy_yyyyy[k] * ab_x + g_z_0_z_0_xy_xyyyyy[k];

                g_z_0_z_0_xxy_yyyyz[k] = -g_z_0_z_0_xy_yyyyz[k] * ab_x + g_z_0_z_0_xy_xyyyyz[k];

                g_z_0_z_0_xxy_yyyzz[k] = -g_z_0_z_0_xy_yyyzz[k] * ab_x + g_z_0_z_0_xy_xyyyzz[k];

                g_z_0_z_0_xxy_yyzzz[k] = -g_z_0_z_0_xy_yyzzz[k] * ab_x + g_z_0_z_0_xy_xyyzzz[k];

                g_z_0_z_0_xxy_yzzzz[k] = -g_z_0_z_0_xy_yzzzz[k] * ab_x + g_z_0_z_0_xy_xyzzzz[k];

                g_z_0_z_0_xxy_zzzzz[k] = -g_z_0_z_0_xy_zzzzz[k] * ab_x + g_z_0_z_0_xy_xzzzzz[k];
            }

            /// Set up 1722-1743 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxz_xxxxx = cbuffer.data(fh_geom_1010_off + 1722 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxxy = cbuffer.data(fh_geom_1010_off + 1723 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxxz = cbuffer.data(fh_geom_1010_off + 1724 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxyy = cbuffer.data(fh_geom_1010_off + 1725 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxyz = cbuffer.data(fh_geom_1010_off + 1726 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxxzz = cbuffer.data(fh_geom_1010_off + 1727 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxyyy = cbuffer.data(fh_geom_1010_off + 1728 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxyyz = cbuffer.data(fh_geom_1010_off + 1729 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxyzz = cbuffer.data(fh_geom_1010_off + 1730 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xxzzz = cbuffer.data(fh_geom_1010_off + 1731 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyyyy = cbuffer.data(fh_geom_1010_off + 1732 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyyyz = cbuffer.data(fh_geom_1010_off + 1733 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyyzz = cbuffer.data(fh_geom_1010_off + 1734 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xyzzz = cbuffer.data(fh_geom_1010_off + 1735 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xzzzz = cbuffer.data(fh_geom_1010_off + 1736 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyyyy = cbuffer.data(fh_geom_1010_off + 1737 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyyyz = cbuffer.data(fh_geom_1010_off + 1738 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyyzz = cbuffer.data(fh_geom_1010_off + 1739 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yyzzz = cbuffer.data(fh_geom_1010_off + 1740 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yzzzz = cbuffer.data(fh_geom_1010_off + 1741 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_zzzzz = cbuffer.data(fh_geom_1010_off + 1742 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxz_xxxxx, g_z_0_z_0_xxz_xxxxy, g_z_0_z_0_xxz_xxxxz, g_z_0_z_0_xxz_xxxyy, g_z_0_z_0_xxz_xxxyz, g_z_0_z_0_xxz_xxxzz, g_z_0_z_0_xxz_xxyyy, g_z_0_z_0_xxz_xxyyz, g_z_0_z_0_xxz_xxyzz, g_z_0_z_0_xxz_xxzzz, g_z_0_z_0_xxz_xyyyy, g_z_0_z_0_xxz_xyyyz, g_z_0_z_0_xxz_xyyzz, g_z_0_z_0_xxz_xyzzz, g_z_0_z_0_xxz_xzzzz, g_z_0_z_0_xxz_yyyyy, g_z_0_z_0_xxz_yyyyz, g_z_0_z_0_xxz_yyyzz, g_z_0_z_0_xxz_yyzzz, g_z_0_z_0_xxz_yzzzz, g_z_0_z_0_xxz_zzzzz, g_z_0_z_0_xz_xxxxx, g_z_0_z_0_xz_xxxxxx, g_z_0_z_0_xz_xxxxxy, g_z_0_z_0_xz_xxxxxz, g_z_0_z_0_xz_xxxxy, g_z_0_z_0_xz_xxxxyy, g_z_0_z_0_xz_xxxxyz, g_z_0_z_0_xz_xxxxz, g_z_0_z_0_xz_xxxxzz, g_z_0_z_0_xz_xxxyy, g_z_0_z_0_xz_xxxyyy, g_z_0_z_0_xz_xxxyyz, g_z_0_z_0_xz_xxxyz, g_z_0_z_0_xz_xxxyzz, g_z_0_z_0_xz_xxxzz, g_z_0_z_0_xz_xxxzzz, g_z_0_z_0_xz_xxyyy, g_z_0_z_0_xz_xxyyyy, g_z_0_z_0_xz_xxyyyz, g_z_0_z_0_xz_xxyyz, g_z_0_z_0_xz_xxyyzz, g_z_0_z_0_xz_xxyzz, g_z_0_z_0_xz_xxyzzz, g_z_0_z_0_xz_xxzzz, g_z_0_z_0_xz_xxzzzz, g_z_0_z_0_xz_xyyyy, g_z_0_z_0_xz_xyyyyy, g_z_0_z_0_xz_xyyyyz, g_z_0_z_0_xz_xyyyz, g_z_0_z_0_xz_xyyyzz, g_z_0_z_0_xz_xyyzz, g_z_0_z_0_xz_xyyzzz, g_z_0_z_0_xz_xyzzz, g_z_0_z_0_xz_xyzzzz, g_z_0_z_0_xz_xzzzz, g_z_0_z_0_xz_xzzzzz, g_z_0_z_0_xz_yyyyy, g_z_0_z_0_xz_yyyyz, g_z_0_z_0_xz_yyyzz, g_z_0_z_0_xz_yyzzz, g_z_0_z_0_xz_yzzzz, g_z_0_z_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxz_xxxxx[k] = -g_z_0_z_0_xz_xxxxx[k] * ab_x + g_z_0_z_0_xz_xxxxxx[k];

                g_z_0_z_0_xxz_xxxxy[k] = -g_z_0_z_0_xz_xxxxy[k] * ab_x + g_z_0_z_0_xz_xxxxxy[k];

                g_z_0_z_0_xxz_xxxxz[k] = -g_z_0_z_0_xz_xxxxz[k] * ab_x + g_z_0_z_0_xz_xxxxxz[k];

                g_z_0_z_0_xxz_xxxyy[k] = -g_z_0_z_0_xz_xxxyy[k] * ab_x + g_z_0_z_0_xz_xxxxyy[k];

                g_z_0_z_0_xxz_xxxyz[k] = -g_z_0_z_0_xz_xxxyz[k] * ab_x + g_z_0_z_0_xz_xxxxyz[k];

                g_z_0_z_0_xxz_xxxzz[k] = -g_z_0_z_0_xz_xxxzz[k] * ab_x + g_z_0_z_0_xz_xxxxzz[k];

                g_z_0_z_0_xxz_xxyyy[k] = -g_z_0_z_0_xz_xxyyy[k] * ab_x + g_z_0_z_0_xz_xxxyyy[k];

                g_z_0_z_0_xxz_xxyyz[k] = -g_z_0_z_0_xz_xxyyz[k] * ab_x + g_z_0_z_0_xz_xxxyyz[k];

                g_z_0_z_0_xxz_xxyzz[k] = -g_z_0_z_0_xz_xxyzz[k] * ab_x + g_z_0_z_0_xz_xxxyzz[k];

                g_z_0_z_0_xxz_xxzzz[k] = -g_z_0_z_0_xz_xxzzz[k] * ab_x + g_z_0_z_0_xz_xxxzzz[k];

                g_z_0_z_0_xxz_xyyyy[k] = -g_z_0_z_0_xz_xyyyy[k] * ab_x + g_z_0_z_0_xz_xxyyyy[k];

                g_z_0_z_0_xxz_xyyyz[k] = -g_z_0_z_0_xz_xyyyz[k] * ab_x + g_z_0_z_0_xz_xxyyyz[k];

                g_z_0_z_0_xxz_xyyzz[k] = -g_z_0_z_0_xz_xyyzz[k] * ab_x + g_z_0_z_0_xz_xxyyzz[k];

                g_z_0_z_0_xxz_xyzzz[k] = -g_z_0_z_0_xz_xyzzz[k] * ab_x + g_z_0_z_0_xz_xxyzzz[k];

                g_z_0_z_0_xxz_xzzzz[k] = -g_z_0_z_0_xz_xzzzz[k] * ab_x + g_z_0_z_0_xz_xxzzzz[k];

                g_z_0_z_0_xxz_yyyyy[k] = -g_z_0_z_0_xz_yyyyy[k] * ab_x + g_z_0_z_0_xz_xyyyyy[k];

                g_z_0_z_0_xxz_yyyyz[k] = -g_z_0_z_0_xz_yyyyz[k] * ab_x + g_z_0_z_0_xz_xyyyyz[k];

                g_z_0_z_0_xxz_yyyzz[k] = -g_z_0_z_0_xz_yyyzz[k] * ab_x + g_z_0_z_0_xz_xyyyzz[k];

                g_z_0_z_0_xxz_yyzzz[k] = -g_z_0_z_0_xz_yyzzz[k] * ab_x + g_z_0_z_0_xz_xyyzzz[k];

                g_z_0_z_0_xxz_yzzzz[k] = -g_z_0_z_0_xz_yzzzz[k] * ab_x + g_z_0_z_0_xz_xyzzzz[k];

                g_z_0_z_0_xxz_zzzzz[k] = -g_z_0_z_0_xz_zzzzz[k] * ab_x + g_z_0_z_0_xz_xzzzzz[k];
            }

            /// Set up 1743-1764 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1743 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1744 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1745 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1746 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1747 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1748 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1749 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1750 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1751 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1752 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1753 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1754 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1755 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1756 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1757 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1758 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1759 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1760 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1761 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1762 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1763 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyy_xxxxx, g_z_0_z_0_xyy_xxxxy, g_z_0_z_0_xyy_xxxxz, g_z_0_z_0_xyy_xxxyy, g_z_0_z_0_xyy_xxxyz, g_z_0_z_0_xyy_xxxzz, g_z_0_z_0_xyy_xxyyy, g_z_0_z_0_xyy_xxyyz, g_z_0_z_0_xyy_xxyzz, g_z_0_z_0_xyy_xxzzz, g_z_0_z_0_xyy_xyyyy, g_z_0_z_0_xyy_xyyyz, g_z_0_z_0_xyy_xyyzz, g_z_0_z_0_xyy_xyzzz, g_z_0_z_0_xyy_xzzzz, g_z_0_z_0_xyy_yyyyy, g_z_0_z_0_xyy_yyyyz, g_z_0_z_0_xyy_yyyzz, g_z_0_z_0_xyy_yyzzz, g_z_0_z_0_xyy_yzzzz, g_z_0_z_0_xyy_zzzzz, g_z_0_z_0_yy_xxxxx, g_z_0_z_0_yy_xxxxxx, g_z_0_z_0_yy_xxxxxy, g_z_0_z_0_yy_xxxxxz, g_z_0_z_0_yy_xxxxy, g_z_0_z_0_yy_xxxxyy, g_z_0_z_0_yy_xxxxyz, g_z_0_z_0_yy_xxxxz, g_z_0_z_0_yy_xxxxzz, g_z_0_z_0_yy_xxxyy, g_z_0_z_0_yy_xxxyyy, g_z_0_z_0_yy_xxxyyz, g_z_0_z_0_yy_xxxyz, g_z_0_z_0_yy_xxxyzz, g_z_0_z_0_yy_xxxzz, g_z_0_z_0_yy_xxxzzz, g_z_0_z_0_yy_xxyyy, g_z_0_z_0_yy_xxyyyy, g_z_0_z_0_yy_xxyyyz, g_z_0_z_0_yy_xxyyz, g_z_0_z_0_yy_xxyyzz, g_z_0_z_0_yy_xxyzz, g_z_0_z_0_yy_xxyzzz, g_z_0_z_0_yy_xxzzz, g_z_0_z_0_yy_xxzzzz, g_z_0_z_0_yy_xyyyy, g_z_0_z_0_yy_xyyyyy, g_z_0_z_0_yy_xyyyyz, g_z_0_z_0_yy_xyyyz, g_z_0_z_0_yy_xyyyzz, g_z_0_z_0_yy_xyyzz, g_z_0_z_0_yy_xyyzzz, g_z_0_z_0_yy_xyzzz, g_z_0_z_0_yy_xyzzzz, g_z_0_z_0_yy_xzzzz, g_z_0_z_0_yy_xzzzzz, g_z_0_z_0_yy_yyyyy, g_z_0_z_0_yy_yyyyz, g_z_0_z_0_yy_yyyzz, g_z_0_z_0_yy_yyzzz, g_z_0_z_0_yy_yzzzz, g_z_0_z_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyy_xxxxx[k] = -g_z_0_z_0_yy_xxxxx[k] * ab_x + g_z_0_z_0_yy_xxxxxx[k];

                g_z_0_z_0_xyy_xxxxy[k] = -g_z_0_z_0_yy_xxxxy[k] * ab_x + g_z_0_z_0_yy_xxxxxy[k];

                g_z_0_z_0_xyy_xxxxz[k] = -g_z_0_z_0_yy_xxxxz[k] * ab_x + g_z_0_z_0_yy_xxxxxz[k];

                g_z_0_z_0_xyy_xxxyy[k] = -g_z_0_z_0_yy_xxxyy[k] * ab_x + g_z_0_z_0_yy_xxxxyy[k];

                g_z_0_z_0_xyy_xxxyz[k] = -g_z_0_z_0_yy_xxxyz[k] * ab_x + g_z_0_z_0_yy_xxxxyz[k];

                g_z_0_z_0_xyy_xxxzz[k] = -g_z_0_z_0_yy_xxxzz[k] * ab_x + g_z_0_z_0_yy_xxxxzz[k];

                g_z_0_z_0_xyy_xxyyy[k] = -g_z_0_z_0_yy_xxyyy[k] * ab_x + g_z_0_z_0_yy_xxxyyy[k];

                g_z_0_z_0_xyy_xxyyz[k] = -g_z_0_z_0_yy_xxyyz[k] * ab_x + g_z_0_z_0_yy_xxxyyz[k];

                g_z_0_z_0_xyy_xxyzz[k] = -g_z_0_z_0_yy_xxyzz[k] * ab_x + g_z_0_z_0_yy_xxxyzz[k];

                g_z_0_z_0_xyy_xxzzz[k] = -g_z_0_z_0_yy_xxzzz[k] * ab_x + g_z_0_z_0_yy_xxxzzz[k];

                g_z_0_z_0_xyy_xyyyy[k] = -g_z_0_z_0_yy_xyyyy[k] * ab_x + g_z_0_z_0_yy_xxyyyy[k];

                g_z_0_z_0_xyy_xyyyz[k] = -g_z_0_z_0_yy_xyyyz[k] * ab_x + g_z_0_z_0_yy_xxyyyz[k];

                g_z_0_z_0_xyy_xyyzz[k] = -g_z_0_z_0_yy_xyyzz[k] * ab_x + g_z_0_z_0_yy_xxyyzz[k];

                g_z_0_z_0_xyy_xyzzz[k] = -g_z_0_z_0_yy_xyzzz[k] * ab_x + g_z_0_z_0_yy_xxyzzz[k];

                g_z_0_z_0_xyy_xzzzz[k] = -g_z_0_z_0_yy_xzzzz[k] * ab_x + g_z_0_z_0_yy_xxzzzz[k];

                g_z_0_z_0_xyy_yyyyy[k] = -g_z_0_z_0_yy_yyyyy[k] * ab_x + g_z_0_z_0_yy_xyyyyy[k];

                g_z_0_z_0_xyy_yyyyz[k] = -g_z_0_z_0_yy_yyyyz[k] * ab_x + g_z_0_z_0_yy_xyyyyz[k];

                g_z_0_z_0_xyy_yyyzz[k] = -g_z_0_z_0_yy_yyyzz[k] * ab_x + g_z_0_z_0_yy_xyyyzz[k];

                g_z_0_z_0_xyy_yyzzz[k] = -g_z_0_z_0_yy_yyzzz[k] * ab_x + g_z_0_z_0_yy_xyyzzz[k];

                g_z_0_z_0_xyy_yzzzz[k] = -g_z_0_z_0_yy_yzzzz[k] * ab_x + g_z_0_z_0_yy_xyzzzz[k];

                g_z_0_z_0_xyy_zzzzz[k] = -g_z_0_z_0_yy_zzzzz[k] * ab_x + g_z_0_z_0_yy_xzzzzz[k];
            }

            /// Set up 1764-1785 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1764 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1765 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1766 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1767 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1768 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1769 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1770 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1771 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1772 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1773 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1774 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1775 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1776 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1777 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1778 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1779 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1780 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1781 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1782 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1783 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1784 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyz_xxxxx, g_z_0_z_0_xyz_xxxxy, g_z_0_z_0_xyz_xxxxz, g_z_0_z_0_xyz_xxxyy, g_z_0_z_0_xyz_xxxyz, g_z_0_z_0_xyz_xxxzz, g_z_0_z_0_xyz_xxyyy, g_z_0_z_0_xyz_xxyyz, g_z_0_z_0_xyz_xxyzz, g_z_0_z_0_xyz_xxzzz, g_z_0_z_0_xyz_xyyyy, g_z_0_z_0_xyz_xyyyz, g_z_0_z_0_xyz_xyyzz, g_z_0_z_0_xyz_xyzzz, g_z_0_z_0_xyz_xzzzz, g_z_0_z_0_xyz_yyyyy, g_z_0_z_0_xyz_yyyyz, g_z_0_z_0_xyz_yyyzz, g_z_0_z_0_xyz_yyzzz, g_z_0_z_0_xyz_yzzzz, g_z_0_z_0_xyz_zzzzz, g_z_0_z_0_yz_xxxxx, g_z_0_z_0_yz_xxxxxx, g_z_0_z_0_yz_xxxxxy, g_z_0_z_0_yz_xxxxxz, g_z_0_z_0_yz_xxxxy, g_z_0_z_0_yz_xxxxyy, g_z_0_z_0_yz_xxxxyz, g_z_0_z_0_yz_xxxxz, g_z_0_z_0_yz_xxxxzz, g_z_0_z_0_yz_xxxyy, g_z_0_z_0_yz_xxxyyy, g_z_0_z_0_yz_xxxyyz, g_z_0_z_0_yz_xxxyz, g_z_0_z_0_yz_xxxyzz, g_z_0_z_0_yz_xxxzz, g_z_0_z_0_yz_xxxzzz, g_z_0_z_0_yz_xxyyy, g_z_0_z_0_yz_xxyyyy, g_z_0_z_0_yz_xxyyyz, g_z_0_z_0_yz_xxyyz, g_z_0_z_0_yz_xxyyzz, g_z_0_z_0_yz_xxyzz, g_z_0_z_0_yz_xxyzzz, g_z_0_z_0_yz_xxzzz, g_z_0_z_0_yz_xxzzzz, g_z_0_z_0_yz_xyyyy, g_z_0_z_0_yz_xyyyyy, g_z_0_z_0_yz_xyyyyz, g_z_0_z_0_yz_xyyyz, g_z_0_z_0_yz_xyyyzz, g_z_0_z_0_yz_xyyzz, g_z_0_z_0_yz_xyyzzz, g_z_0_z_0_yz_xyzzz, g_z_0_z_0_yz_xyzzzz, g_z_0_z_0_yz_xzzzz, g_z_0_z_0_yz_xzzzzz, g_z_0_z_0_yz_yyyyy, g_z_0_z_0_yz_yyyyz, g_z_0_z_0_yz_yyyzz, g_z_0_z_0_yz_yyzzz, g_z_0_z_0_yz_yzzzz, g_z_0_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyz_xxxxx[k] = -g_z_0_z_0_yz_xxxxx[k] * ab_x + g_z_0_z_0_yz_xxxxxx[k];

                g_z_0_z_0_xyz_xxxxy[k] = -g_z_0_z_0_yz_xxxxy[k] * ab_x + g_z_0_z_0_yz_xxxxxy[k];

                g_z_0_z_0_xyz_xxxxz[k] = -g_z_0_z_0_yz_xxxxz[k] * ab_x + g_z_0_z_0_yz_xxxxxz[k];

                g_z_0_z_0_xyz_xxxyy[k] = -g_z_0_z_0_yz_xxxyy[k] * ab_x + g_z_0_z_0_yz_xxxxyy[k];

                g_z_0_z_0_xyz_xxxyz[k] = -g_z_0_z_0_yz_xxxyz[k] * ab_x + g_z_0_z_0_yz_xxxxyz[k];

                g_z_0_z_0_xyz_xxxzz[k] = -g_z_0_z_0_yz_xxxzz[k] * ab_x + g_z_0_z_0_yz_xxxxzz[k];

                g_z_0_z_0_xyz_xxyyy[k] = -g_z_0_z_0_yz_xxyyy[k] * ab_x + g_z_0_z_0_yz_xxxyyy[k];

                g_z_0_z_0_xyz_xxyyz[k] = -g_z_0_z_0_yz_xxyyz[k] * ab_x + g_z_0_z_0_yz_xxxyyz[k];

                g_z_0_z_0_xyz_xxyzz[k] = -g_z_0_z_0_yz_xxyzz[k] * ab_x + g_z_0_z_0_yz_xxxyzz[k];

                g_z_0_z_0_xyz_xxzzz[k] = -g_z_0_z_0_yz_xxzzz[k] * ab_x + g_z_0_z_0_yz_xxxzzz[k];

                g_z_0_z_0_xyz_xyyyy[k] = -g_z_0_z_0_yz_xyyyy[k] * ab_x + g_z_0_z_0_yz_xxyyyy[k];

                g_z_0_z_0_xyz_xyyyz[k] = -g_z_0_z_0_yz_xyyyz[k] * ab_x + g_z_0_z_0_yz_xxyyyz[k];

                g_z_0_z_0_xyz_xyyzz[k] = -g_z_0_z_0_yz_xyyzz[k] * ab_x + g_z_0_z_0_yz_xxyyzz[k];

                g_z_0_z_0_xyz_xyzzz[k] = -g_z_0_z_0_yz_xyzzz[k] * ab_x + g_z_0_z_0_yz_xxyzzz[k];

                g_z_0_z_0_xyz_xzzzz[k] = -g_z_0_z_0_yz_xzzzz[k] * ab_x + g_z_0_z_0_yz_xxzzzz[k];

                g_z_0_z_0_xyz_yyyyy[k] = -g_z_0_z_0_yz_yyyyy[k] * ab_x + g_z_0_z_0_yz_xyyyyy[k];

                g_z_0_z_0_xyz_yyyyz[k] = -g_z_0_z_0_yz_yyyyz[k] * ab_x + g_z_0_z_0_yz_xyyyyz[k];

                g_z_0_z_0_xyz_yyyzz[k] = -g_z_0_z_0_yz_yyyzz[k] * ab_x + g_z_0_z_0_yz_xyyyzz[k];

                g_z_0_z_0_xyz_yyzzz[k] = -g_z_0_z_0_yz_yyzzz[k] * ab_x + g_z_0_z_0_yz_xyyzzz[k];

                g_z_0_z_0_xyz_yzzzz[k] = -g_z_0_z_0_yz_yzzzz[k] * ab_x + g_z_0_z_0_yz_xyzzzz[k];

                g_z_0_z_0_xyz_zzzzz[k] = -g_z_0_z_0_yz_zzzzz[k] * ab_x + g_z_0_z_0_yz_xzzzzz[k];
            }

            /// Set up 1785-1806 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1785 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1786 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1787 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1788 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1789 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1790 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1791 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1792 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1793 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1794 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1795 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1796 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1797 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1798 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1799 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1800 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1801 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1802 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1803 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1804 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1805 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzz_xxxxx, g_z_0_z_0_xzz_xxxxy, g_z_0_z_0_xzz_xxxxz, g_z_0_z_0_xzz_xxxyy, g_z_0_z_0_xzz_xxxyz, g_z_0_z_0_xzz_xxxzz, g_z_0_z_0_xzz_xxyyy, g_z_0_z_0_xzz_xxyyz, g_z_0_z_0_xzz_xxyzz, g_z_0_z_0_xzz_xxzzz, g_z_0_z_0_xzz_xyyyy, g_z_0_z_0_xzz_xyyyz, g_z_0_z_0_xzz_xyyzz, g_z_0_z_0_xzz_xyzzz, g_z_0_z_0_xzz_xzzzz, g_z_0_z_0_xzz_yyyyy, g_z_0_z_0_xzz_yyyyz, g_z_0_z_0_xzz_yyyzz, g_z_0_z_0_xzz_yyzzz, g_z_0_z_0_xzz_yzzzz, g_z_0_z_0_xzz_zzzzz, g_z_0_z_0_zz_xxxxx, g_z_0_z_0_zz_xxxxxx, g_z_0_z_0_zz_xxxxxy, g_z_0_z_0_zz_xxxxxz, g_z_0_z_0_zz_xxxxy, g_z_0_z_0_zz_xxxxyy, g_z_0_z_0_zz_xxxxyz, g_z_0_z_0_zz_xxxxz, g_z_0_z_0_zz_xxxxzz, g_z_0_z_0_zz_xxxyy, g_z_0_z_0_zz_xxxyyy, g_z_0_z_0_zz_xxxyyz, g_z_0_z_0_zz_xxxyz, g_z_0_z_0_zz_xxxyzz, g_z_0_z_0_zz_xxxzz, g_z_0_z_0_zz_xxxzzz, g_z_0_z_0_zz_xxyyy, g_z_0_z_0_zz_xxyyyy, g_z_0_z_0_zz_xxyyyz, g_z_0_z_0_zz_xxyyz, g_z_0_z_0_zz_xxyyzz, g_z_0_z_0_zz_xxyzz, g_z_0_z_0_zz_xxyzzz, g_z_0_z_0_zz_xxzzz, g_z_0_z_0_zz_xxzzzz, g_z_0_z_0_zz_xyyyy, g_z_0_z_0_zz_xyyyyy, g_z_0_z_0_zz_xyyyyz, g_z_0_z_0_zz_xyyyz, g_z_0_z_0_zz_xyyyzz, g_z_0_z_0_zz_xyyzz, g_z_0_z_0_zz_xyyzzz, g_z_0_z_0_zz_xyzzz, g_z_0_z_0_zz_xyzzzz, g_z_0_z_0_zz_xzzzz, g_z_0_z_0_zz_xzzzzz, g_z_0_z_0_zz_yyyyy, g_z_0_z_0_zz_yyyyz, g_z_0_z_0_zz_yyyzz, g_z_0_z_0_zz_yyzzz, g_z_0_z_0_zz_yzzzz, g_z_0_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzz_xxxxx[k] = -g_z_0_z_0_zz_xxxxx[k] * ab_x + g_z_0_z_0_zz_xxxxxx[k];

                g_z_0_z_0_xzz_xxxxy[k] = -g_z_0_z_0_zz_xxxxy[k] * ab_x + g_z_0_z_0_zz_xxxxxy[k];

                g_z_0_z_0_xzz_xxxxz[k] = -g_z_0_z_0_zz_xxxxz[k] * ab_x + g_z_0_z_0_zz_xxxxxz[k];

                g_z_0_z_0_xzz_xxxyy[k] = -g_z_0_z_0_zz_xxxyy[k] * ab_x + g_z_0_z_0_zz_xxxxyy[k];

                g_z_0_z_0_xzz_xxxyz[k] = -g_z_0_z_0_zz_xxxyz[k] * ab_x + g_z_0_z_0_zz_xxxxyz[k];

                g_z_0_z_0_xzz_xxxzz[k] = -g_z_0_z_0_zz_xxxzz[k] * ab_x + g_z_0_z_0_zz_xxxxzz[k];

                g_z_0_z_0_xzz_xxyyy[k] = -g_z_0_z_0_zz_xxyyy[k] * ab_x + g_z_0_z_0_zz_xxxyyy[k];

                g_z_0_z_0_xzz_xxyyz[k] = -g_z_0_z_0_zz_xxyyz[k] * ab_x + g_z_0_z_0_zz_xxxyyz[k];

                g_z_0_z_0_xzz_xxyzz[k] = -g_z_0_z_0_zz_xxyzz[k] * ab_x + g_z_0_z_0_zz_xxxyzz[k];

                g_z_0_z_0_xzz_xxzzz[k] = -g_z_0_z_0_zz_xxzzz[k] * ab_x + g_z_0_z_0_zz_xxxzzz[k];

                g_z_0_z_0_xzz_xyyyy[k] = -g_z_0_z_0_zz_xyyyy[k] * ab_x + g_z_0_z_0_zz_xxyyyy[k];

                g_z_0_z_0_xzz_xyyyz[k] = -g_z_0_z_0_zz_xyyyz[k] * ab_x + g_z_0_z_0_zz_xxyyyz[k];

                g_z_0_z_0_xzz_xyyzz[k] = -g_z_0_z_0_zz_xyyzz[k] * ab_x + g_z_0_z_0_zz_xxyyzz[k];

                g_z_0_z_0_xzz_xyzzz[k] = -g_z_0_z_0_zz_xyzzz[k] * ab_x + g_z_0_z_0_zz_xxyzzz[k];

                g_z_0_z_0_xzz_xzzzz[k] = -g_z_0_z_0_zz_xzzzz[k] * ab_x + g_z_0_z_0_zz_xxzzzz[k];

                g_z_0_z_0_xzz_yyyyy[k] = -g_z_0_z_0_zz_yyyyy[k] * ab_x + g_z_0_z_0_zz_xyyyyy[k];

                g_z_0_z_0_xzz_yyyyz[k] = -g_z_0_z_0_zz_yyyyz[k] * ab_x + g_z_0_z_0_zz_xyyyyz[k];

                g_z_0_z_0_xzz_yyyzz[k] = -g_z_0_z_0_zz_yyyzz[k] * ab_x + g_z_0_z_0_zz_xyyyzz[k];

                g_z_0_z_0_xzz_yyzzz[k] = -g_z_0_z_0_zz_yyzzz[k] * ab_x + g_z_0_z_0_zz_xyyzzz[k];

                g_z_0_z_0_xzz_yzzzz[k] = -g_z_0_z_0_zz_yzzzz[k] * ab_x + g_z_0_z_0_zz_xyzzzz[k];

                g_z_0_z_0_xzz_zzzzz[k] = -g_z_0_z_0_zz_zzzzz[k] * ab_x + g_z_0_z_0_zz_xzzzzz[k];
            }

            /// Set up 1806-1827 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyy_xxxxx = cbuffer.data(fh_geom_1010_off + 1806 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxxy = cbuffer.data(fh_geom_1010_off + 1807 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxxz = cbuffer.data(fh_geom_1010_off + 1808 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxyy = cbuffer.data(fh_geom_1010_off + 1809 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxyz = cbuffer.data(fh_geom_1010_off + 1810 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxxzz = cbuffer.data(fh_geom_1010_off + 1811 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxyyy = cbuffer.data(fh_geom_1010_off + 1812 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxyyz = cbuffer.data(fh_geom_1010_off + 1813 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxyzz = cbuffer.data(fh_geom_1010_off + 1814 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xxzzz = cbuffer.data(fh_geom_1010_off + 1815 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyyyy = cbuffer.data(fh_geom_1010_off + 1816 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyyyz = cbuffer.data(fh_geom_1010_off + 1817 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyyzz = cbuffer.data(fh_geom_1010_off + 1818 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xyzzz = cbuffer.data(fh_geom_1010_off + 1819 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xzzzz = cbuffer.data(fh_geom_1010_off + 1820 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyyyy = cbuffer.data(fh_geom_1010_off + 1821 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyyyz = cbuffer.data(fh_geom_1010_off + 1822 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyyzz = cbuffer.data(fh_geom_1010_off + 1823 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yyzzz = cbuffer.data(fh_geom_1010_off + 1824 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yzzzz = cbuffer.data(fh_geom_1010_off + 1825 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_zzzzz = cbuffer.data(fh_geom_1010_off + 1826 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yy_xxxxx, g_z_0_z_0_yy_xxxxxy, g_z_0_z_0_yy_xxxxy, g_z_0_z_0_yy_xxxxyy, g_z_0_z_0_yy_xxxxyz, g_z_0_z_0_yy_xxxxz, g_z_0_z_0_yy_xxxyy, g_z_0_z_0_yy_xxxyyy, g_z_0_z_0_yy_xxxyyz, g_z_0_z_0_yy_xxxyz, g_z_0_z_0_yy_xxxyzz, g_z_0_z_0_yy_xxxzz, g_z_0_z_0_yy_xxyyy, g_z_0_z_0_yy_xxyyyy, g_z_0_z_0_yy_xxyyyz, g_z_0_z_0_yy_xxyyz, g_z_0_z_0_yy_xxyyzz, g_z_0_z_0_yy_xxyzz, g_z_0_z_0_yy_xxyzzz, g_z_0_z_0_yy_xxzzz, g_z_0_z_0_yy_xyyyy, g_z_0_z_0_yy_xyyyyy, g_z_0_z_0_yy_xyyyyz, g_z_0_z_0_yy_xyyyz, g_z_0_z_0_yy_xyyyzz, g_z_0_z_0_yy_xyyzz, g_z_0_z_0_yy_xyyzzz, g_z_0_z_0_yy_xyzzz, g_z_0_z_0_yy_xyzzzz, g_z_0_z_0_yy_xzzzz, g_z_0_z_0_yy_yyyyy, g_z_0_z_0_yy_yyyyyy, g_z_0_z_0_yy_yyyyyz, g_z_0_z_0_yy_yyyyz, g_z_0_z_0_yy_yyyyzz, g_z_0_z_0_yy_yyyzz, g_z_0_z_0_yy_yyyzzz, g_z_0_z_0_yy_yyzzz, g_z_0_z_0_yy_yyzzzz, g_z_0_z_0_yy_yzzzz, g_z_0_z_0_yy_yzzzzz, g_z_0_z_0_yy_zzzzz, g_z_0_z_0_yyy_xxxxx, g_z_0_z_0_yyy_xxxxy, g_z_0_z_0_yyy_xxxxz, g_z_0_z_0_yyy_xxxyy, g_z_0_z_0_yyy_xxxyz, g_z_0_z_0_yyy_xxxzz, g_z_0_z_0_yyy_xxyyy, g_z_0_z_0_yyy_xxyyz, g_z_0_z_0_yyy_xxyzz, g_z_0_z_0_yyy_xxzzz, g_z_0_z_0_yyy_xyyyy, g_z_0_z_0_yyy_xyyyz, g_z_0_z_0_yyy_xyyzz, g_z_0_z_0_yyy_xyzzz, g_z_0_z_0_yyy_xzzzz, g_z_0_z_0_yyy_yyyyy, g_z_0_z_0_yyy_yyyyz, g_z_0_z_0_yyy_yyyzz, g_z_0_z_0_yyy_yyzzz, g_z_0_z_0_yyy_yzzzz, g_z_0_z_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyy_xxxxx[k] = -g_z_0_z_0_yy_xxxxx[k] * ab_y + g_z_0_z_0_yy_xxxxxy[k];

                g_z_0_z_0_yyy_xxxxy[k] = -g_z_0_z_0_yy_xxxxy[k] * ab_y + g_z_0_z_0_yy_xxxxyy[k];

                g_z_0_z_0_yyy_xxxxz[k] = -g_z_0_z_0_yy_xxxxz[k] * ab_y + g_z_0_z_0_yy_xxxxyz[k];

                g_z_0_z_0_yyy_xxxyy[k] = -g_z_0_z_0_yy_xxxyy[k] * ab_y + g_z_0_z_0_yy_xxxyyy[k];

                g_z_0_z_0_yyy_xxxyz[k] = -g_z_0_z_0_yy_xxxyz[k] * ab_y + g_z_0_z_0_yy_xxxyyz[k];

                g_z_0_z_0_yyy_xxxzz[k] = -g_z_0_z_0_yy_xxxzz[k] * ab_y + g_z_0_z_0_yy_xxxyzz[k];

                g_z_0_z_0_yyy_xxyyy[k] = -g_z_0_z_0_yy_xxyyy[k] * ab_y + g_z_0_z_0_yy_xxyyyy[k];

                g_z_0_z_0_yyy_xxyyz[k] = -g_z_0_z_0_yy_xxyyz[k] * ab_y + g_z_0_z_0_yy_xxyyyz[k];

                g_z_0_z_0_yyy_xxyzz[k] = -g_z_0_z_0_yy_xxyzz[k] * ab_y + g_z_0_z_0_yy_xxyyzz[k];

                g_z_0_z_0_yyy_xxzzz[k] = -g_z_0_z_0_yy_xxzzz[k] * ab_y + g_z_0_z_0_yy_xxyzzz[k];

                g_z_0_z_0_yyy_xyyyy[k] = -g_z_0_z_0_yy_xyyyy[k] * ab_y + g_z_0_z_0_yy_xyyyyy[k];

                g_z_0_z_0_yyy_xyyyz[k] = -g_z_0_z_0_yy_xyyyz[k] * ab_y + g_z_0_z_0_yy_xyyyyz[k];

                g_z_0_z_0_yyy_xyyzz[k] = -g_z_0_z_0_yy_xyyzz[k] * ab_y + g_z_0_z_0_yy_xyyyzz[k];

                g_z_0_z_0_yyy_xyzzz[k] = -g_z_0_z_0_yy_xyzzz[k] * ab_y + g_z_0_z_0_yy_xyyzzz[k];

                g_z_0_z_0_yyy_xzzzz[k] = -g_z_0_z_0_yy_xzzzz[k] * ab_y + g_z_0_z_0_yy_xyzzzz[k];

                g_z_0_z_0_yyy_yyyyy[k] = -g_z_0_z_0_yy_yyyyy[k] * ab_y + g_z_0_z_0_yy_yyyyyy[k];

                g_z_0_z_0_yyy_yyyyz[k] = -g_z_0_z_0_yy_yyyyz[k] * ab_y + g_z_0_z_0_yy_yyyyyz[k];

                g_z_0_z_0_yyy_yyyzz[k] = -g_z_0_z_0_yy_yyyzz[k] * ab_y + g_z_0_z_0_yy_yyyyzz[k];

                g_z_0_z_0_yyy_yyzzz[k] = -g_z_0_z_0_yy_yyzzz[k] * ab_y + g_z_0_z_0_yy_yyyzzz[k];

                g_z_0_z_0_yyy_yzzzz[k] = -g_z_0_z_0_yy_yzzzz[k] * ab_y + g_z_0_z_0_yy_yyzzzz[k];

                g_z_0_z_0_yyy_zzzzz[k] = -g_z_0_z_0_yy_zzzzz[k] * ab_y + g_z_0_z_0_yy_yzzzzz[k];
            }

            /// Set up 1827-1848 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyz_xxxxx = cbuffer.data(fh_geom_1010_off + 1827 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxxy = cbuffer.data(fh_geom_1010_off + 1828 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxxz = cbuffer.data(fh_geom_1010_off + 1829 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxyy = cbuffer.data(fh_geom_1010_off + 1830 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxyz = cbuffer.data(fh_geom_1010_off + 1831 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxxzz = cbuffer.data(fh_geom_1010_off + 1832 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxyyy = cbuffer.data(fh_geom_1010_off + 1833 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxyyz = cbuffer.data(fh_geom_1010_off + 1834 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxyzz = cbuffer.data(fh_geom_1010_off + 1835 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xxzzz = cbuffer.data(fh_geom_1010_off + 1836 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyyyy = cbuffer.data(fh_geom_1010_off + 1837 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyyyz = cbuffer.data(fh_geom_1010_off + 1838 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyyzz = cbuffer.data(fh_geom_1010_off + 1839 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xyzzz = cbuffer.data(fh_geom_1010_off + 1840 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xzzzz = cbuffer.data(fh_geom_1010_off + 1841 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyyyy = cbuffer.data(fh_geom_1010_off + 1842 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyyyz = cbuffer.data(fh_geom_1010_off + 1843 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyyzz = cbuffer.data(fh_geom_1010_off + 1844 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yyzzz = cbuffer.data(fh_geom_1010_off + 1845 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yzzzz = cbuffer.data(fh_geom_1010_off + 1846 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_zzzzz = cbuffer.data(fh_geom_1010_off + 1847 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyz_xxxxx, g_z_0_z_0_yyz_xxxxy, g_z_0_z_0_yyz_xxxxz, g_z_0_z_0_yyz_xxxyy, g_z_0_z_0_yyz_xxxyz, g_z_0_z_0_yyz_xxxzz, g_z_0_z_0_yyz_xxyyy, g_z_0_z_0_yyz_xxyyz, g_z_0_z_0_yyz_xxyzz, g_z_0_z_0_yyz_xxzzz, g_z_0_z_0_yyz_xyyyy, g_z_0_z_0_yyz_xyyyz, g_z_0_z_0_yyz_xyyzz, g_z_0_z_0_yyz_xyzzz, g_z_0_z_0_yyz_xzzzz, g_z_0_z_0_yyz_yyyyy, g_z_0_z_0_yyz_yyyyz, g_z_0_z_0_yyz_yyyzz, g_z_0_z_0_yyz_yyzzz, g_z_0_z_0_yyz_yzzzz, g_z_0_z_0_yyz_zzzzz, g_z_0_z_0_yz_xxxxx, g_z_0_z_0_yz_xxxxxy, g_z_0_z_0_yz_xxxxy, g_z_0_z_0_yz_xxxxyy, g_z_0_z_0_yz_xxxxyz, g_z_0_z_0_yz_xxxxz, g_z_0_z_0_yz_xxxyy, g_z_0_z_0_yz_xxxyyy, g_z_0_z_0_yz_xxxyyz, g_z_0_z_0_yz_xxxyz, g_z_0_z_0_yz_xxxyzz, g_z_0_z_0_yz_xxxzz, g_z_0_z_0_yz_xxyyy, g_z_0_z_0_yz_xxyyyy, g_z_0_z_0_yz_xxyyyz, g_z_0_z_0_yz_xxyyz, g_z_0_z_0_yz_xxyyzz, g_z_0_z_0_yz_xxyzz, g_z_0_z_0_yz_xxyzzz, g_z_0_z_0_yz_xxzzz, g_z_0_z_0_yz_xyyyy, g_z_0_z_0_yz_xyyyyy, g_z_0_z_0_yz_xyyyyz, g_z_0_z_0_yz_xyyyz, g_z_0_z_0_yz_xyyyzz, g_z_0_z_0_yz_xyyzz, g_z_0_z_0_yz_xyyzzz, g_z_0_z_0_yz_xyzzz, g_z_0_z_0_yz_xyzzzz, g_z_0_z_0_yz_xzzzz, g_z_0_z_0_yz_yyyyy, g_z_0_z_0_yz_yyyyyy, g_z_0_z_0_yz_yyyyyz, g_z_0_z_0_yz_yyyyz, g_z_0_z_0_yz_yyyyzz, g_z_0_z_0_yz_yyyzz, g_z_0_z_0_yz_yyyzzz, g_z_0_z_0_yz_yyzzz, g_z_0_z_0_yz_yyzzzz, g_z_0_z_0_yz_yzzzz, g_z_0_z_0_yz_yzzzzz, g_z_0_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyz_xxxxx[k] = -g_z_0_z_0_yz_xxxxx[k] * ab_y + g_z_0_z_0_yz_xxxxxy[k];

                g_z_0_z_0_yyz_xxxxy[k] = -g_z_0_z_0_yz_xxxxy[k] * ab_y + g_z_0_z_0_yz_xxxxyy[k];

                g_z_0_z_0_yyz_xxxxz[k] = -g_z_0_z_0_yz_xxxxz[k] * ab_y + g_z_0_z_0_yz_xxxxyz[k];

                g_z_0_z_0_yyz_xxxyy[k] = -g_z_0_z_0_yz_xxxyy[k] * ab_y + g_z_0_z_0_yz_xxxyyy[k];

                g_z_0_z_0_yyz_xxxyz[k] = -g_z_0_z_0_yz_xxxyz[k] * ab_y + g_z_0_z_0_yz_xxxyyz[k];

                g_z_0_z_0_yyz_xxxzz[k] = -g_z_0_z_0_yz_xxxzz[k] * ab_y + g_z_0_z_0_yz_xxxyzz[k];

                g_z_0_z_0_yyz_xxyyy[k] = -g_z_0_z_0_yz_xxyyy[k] * ab_y + g_z_0_z_0_yz_xxyyyy[k];

                g_z_0_z_0_yyz_xxyyz[k] = -g_z_0_z_0_yz_xxyyz[k] * ab_y + g_z_0_z_0_yz_xxyyyz[k];

                g_z_0_z_0_yyz_xxyzz[k] = -g_z_0_z_0_yz_xxyzz[k] * ab_y + g_z_0_z_0_yz_xxyyzz[k];

                g_z_0_z_0_yyz_xxzzz[k] = -g_z_0_z_0_yz_xxzzz[k] * ab_y + g_z_0_z_0_yz_xxyzzz[k];

                g_z_0_z_0_yyz_xyyyy[k] = -g_z_0_z_0_yz_xyyyy[k] * ab_y + g_z_0_z_0_yz_xyyyyy[k];

                g_z_0_z_0_yyz_xyyyz[k] = -g_z_0_z_0_yz_xyyyz[k] * ab_y + g_z_0_z_0_yz_xyyyyz[k];

                g_z_0_z_0_yyz_xyyzz[k] = -g_z_0_z_0_yz_xyyzz[k] * ab_y + g_z_0_z_0_yz_xyyyzz[k];

                g_z_0_z_0_yyz_xyzzz[k] = -g_z_0_z_0_yz_xyzzz[k] * ab_y + g_z_0_z_0_yz_xyyzzz[k];

                g_z_0_z_0_yyz_xzzzz[k] = -g_z_0_z_0_yz_xzzzz[k] * ab_y + g_z_0_z_0_yz_xyzzzz[k];

                g_z_0_z_0_yyz_yyyyy[k] = -g_z_0_z_0_yz_yyyyy[k] * ab_y + g_z_0_z_0_yz_yyyyyy[k];

                g_z_0_z_0_yyz_yyyyz[k] = -g_z_0_z_0_yz_yyyyz[k] * ab_y + g_z_0_z_0_yz_yyyyyz[k];

                g_z_0_z_0_yyz_yyyzz[k] = -g_z_0_z_0_yz_yyyzz[k] * ab_y + g_z_0_z_0_yz_yyyyzz[k];

                g_z_0_z_0_yyz_yyzzz[k] = -g_z_0_z_0_yz_yyzzz[k] * ab_y + g_z_0_z_0_yz_yyyzzz[k];

                g_z_0_z_0_yyz_yzzzz[k] = -g_z_0_z_0_yz_yzzzz[k] * ab_y + g_z_0_z_0_yz_yyzzzz[k];

                g_z_0_z_0_yyz_zzzzz[k] = -g_z_0_z_0_yz_zzzzz[k] * ab_y + g_z_0_z_0_yz_yzzzzz[k];
            }

            /// Set up 1848-1869 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1848 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1849 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1850 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1851 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1852 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1853 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1854 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1855 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1856 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1857 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1858 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1859 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1860 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1861 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1862 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1863 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1864 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1865 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1866 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1867 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1868 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzz_xxxxx, g_z_0_z_0_yzz_xxxxy, g_z_0_z_0_yzz_xxxxz, g_z_0_z_0_yzz_xxxyy, g_z_0_z_0_yzz_xxxyz, g_z_0_z_0_yzz_xxxzz, g_z_0_z_0_yzz_xxyyy, g_z_0_z_0_yzz_xxyyz, g_z_0_z_0_yzz_xxyzz, g_z_0_z_0_yzz_xxzzz, g_z_0_z_0_yzz_xyyyy, g_z_0_z_0_yzz_xyyyz, g_z_0_z_0_yzz_xyyzz, g_z_0_z_0_yzz_xyzzz, g_z_0_z_0_yzz_xzzzz, g_z_0_z_0_yzz_yyyyy, g_z_0_z_0_yzz_yyyyz, g_z_0_z_0_yzz_yyyzz, g_z_0_z_0_yzz_yyzzz, g_z_0_z_0_yzz_yzzzz, g_z_0_z_0_yzz_zzzzz, g_z_0_z_0_zz_xxxxx, g_z_0_z_0_zz_xxxxxy, g_z_0_z_0_zz_xxxxy, g_z_0_z_0_zz_xxxxyy, g_z_0_z_0_zz_xxxxyz, g_z_0_z_0_zz_xxxxz, g_z_0_z_0_zz_xxxyy, g_z_0_z_0_zz_xxxyyy, g_z_0_z_0_zz_xxxyyz, g_z_0_z_0_zz_xxxyz, g_z_0_z_0_zz_xxxyzz, g_z_0_z_0_zz_xxxzz, g_z_0_z_0_zz_xxyyy, g_z_0_z_0_zz_xxyyyy, g_z_0_z_0_zz_xxyyyz, g_z_0_z_0_zz_xxyyz, g_z_0_z_0_zz_xxyyzz, g_z_0_z_0_zz_xxyzz, g_z_0_z_0_zz_xxyzzz, g_z_0_z_0_zz_xxzzz, g_z_0_z_0_zz_xyyyy, g_z_0_z_0_zz_xyyyyy, g_z_0_z_0_zz_xyyyyz, g_z_0_z_0_zz_xyyyz, g_z_0_z_0_zz_xyyyzz, g_z_0_z_0_zz_xyyzz, g_z_0_z_0_zz_xyyzzz, g_z_0_z_0_zz_xyzzz, g_z_0_z_0_zz_xyzzzz, g_z_0_z_0_zz_xzzzz, g_z_0_z_0_zz_yyyyy, g_z_0_z_0_zz_yyyyyy, g_z_0_z_0_zz_yyyyyz, g_z_0_z_0_zz_yyyyz, g_z_0_z_0_zz_yyyyzz, g_z_0_z_0_zz_yyyzz, g_z_0_z_0_zz_yyyzzz, g_z_0_z_0_zz_yyzzz, g_z_0_z_0_zz_yyzzzz, g_z_0_z_0_zz_yzzzz, g_z_0_z_0_zz_yzzzzz, g_z_0_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzz_xxxxx[k] = -g_z_0_z_0_zz_xxxxx[k] * ab_y + g_z_0_z_0_zz_xxxxxy[k];

                g_z_0_z_0_yzz_xxxxy[k] = -g_z_0_z_0_zz_xxxxy[k] * ab_y + g_z_0_z_0_zz_xxxxyy[k];

                g_z_0_z_0_yzz_xxxxz[k] = -g_z_0_z_0_zz_xxxxz[k] * ab_y + g_z_0_z_0_zz_xxxxyz[k];

                g_z_0_z_0_yzz_xxxyy[k] = -g_z_0_z_0_zz_xxxyy[k] * ab_y + g_z_0_z_0_zz_xxxyyy[k];

                g_z_0_z_0_yzz_xxxyz[k] = -g_z_0_z_0_zz_xxxyz[k] * ab_y + g_z_0_z_0_zz_xxxyyz[k];

                g_z_0_z_0_yzz_xxxzz[k] = -g_z_0_z_0_zz_xxxzz[k] * ab_y + g_z_0_z_0_zz_xxxyzz[k];

                g_z_0_z_0_yzz_xxyyy[k] = -g_z_0_z_0_zz_xxyyy[k] * ab_y + g_z_0_z_0_zz_xxyyyy[k];

                g_z_0_z_0_yzz_xxyyz[k] = -g_z_0_z_0_zz_xxyyz[k] * ab_y + g_z_0_z_0_zz_xxyyyz[k];

                g_z_0_z_0_yzz_xxyzz[k] = -g_z_0_z_0_zz_xxyzz[k] * ab_y + g_z_0_z_0_zz_xxyyzz[k];

                g_z_0_z_0_yzz_xxzzz[k] = -g_z_0_z_0_zz_xxzzz[k] * ab_y + g_z_0_z_0_zz_xxyzzz[k];

                g_z_0_z_0_yzz_xyyyy[k] = -g_z_0_z_0_zz_xyyyy[k] * ab_y + g_z_0_z_0_zz_xyyyyy[k];

                g_z_0_z_0_yzz_xyyyz[k] = -g_z_0_z_0_zz_xyyyz[k] * ab_y + g_z_0_z_0_zz_xyyyyz[k];

                g_z_0_z_0_yzz_xyyzz[k] = -g_z_0_z_0_zz_xyyzz[k] * ab_y + g_z_0_z_0_zz_xyyyzz[k];

                g_z_0_z_0_yzz_xyzzz[k] = -g_z_0_z_0_zz_xyzzz[k] * ab_y + g_z_0_z_0_zz_xyyzzz[k];

                g_z_0_z_0_yzz_xzzzz[k] = -g_z_0_z_0_zz_xzzzz[k] * ab_y + g_z_0_z_0_zz_xyzzzz[k];

                g_z_0_z_0_yzz_yyyyy[k] = -g_z_0_z_0_zz_yyyyy[k] * ab_y + g_z_0_z_0_zz_yyyyyy[k];

                g_z_0_z_0_yzz_yyyyz[k] = -g_z_0_z_0_zz_yyyyz[k] * ab_y + g_z_0_z_0_zz_yyyyyz[k];

                g_z_0_z_0_yzz_yyyzz[k] = -g_z_0_z_0_zz_yyyzz[k] * ab_y + g_z_0_z_0_zz_yyyyzz[k];

                g_z_0_z_0_yzz_yyzzz[k] = -g_z_0_z_0_zz_yyzzz[k] * ab_y + g_z_0_z_0_zz_yyyzzz[k];

                g_z_0_z_0_yzz_yzzzz[k] = -g_z_0_z_0_zz_yzzzz[k] * ab_y + g_z_0_z_0_zz_yyzzzz[k];

                g_z_0_z_0_yzz_zzzzz[k] = -g_z_0_z_0_zz_zzzzz[k] * ab_y + g_z_0_z_0_zz_yzzzzz[k];
            }

            /// Set up 1869-1890 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzz_xxxxx = cbuffer.data(fh_geom_1010_off + 1869 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxxy = cbuffer.data(fh_geom_1010_off + 1870 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxxz = cbuffer.data(fh_geom_1010_off + 1871 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxyy = cbuffer.data(fh_geom_1010_off + 1872 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxyz = cbuffer.data(fh_geom_1010_off + 1873 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxxzz = cbuffer.data(fh_geom_1010_off + 1874 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxyyy = cbuffer.data(fh_geom_1010_off + 1875 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxyyz = cbuffer.data(fh_geom_1010_off + 1876 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxyzz = cbuffer.data(fh_geom_1010_off + 1877 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xxzzz = cbuffer.data(fh_geom_1010_off + 1878 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyyyy = cbuffer.data(fh_geom_1010_off + 1879 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyyyz = cbuffer.data(fh_geom_1010_off + 1880 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyyzz = cbuffer.data(fh_geom_1010_off + 1881 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xyzzz = cbuffer.data(fh_geom_1010_off + 1882 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xzzzz = cbuffer.data(fh_geom_1010_off + 1883 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyyyy = cbuffer.data(fh_geom_1010_off + 1884 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyyyz = cbuffer.data(fh_geom_1010_off + 1885 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyyzz = cbuffer.data(fh_geom_1010_off + 1886 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yyzzz = cbuffer.data(fh_geom_1010_off + 1887 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yzzzz = cbuffer.data(fh_geom_1010_off + 1888 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_zzzzz = cbuffer.data(fh_geom_1010_off + 1889 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zz_xxxxx, g_0_0_z_0_zz_xxxxy, g_0_0_z_0_zz_xxxxz, g_0_0_z_0_zz_xxxyy, g_0_0_z_0_zz_xxxyz, g_0_0_z_0_zz_xxxzz, g_0_0_z_0_zz_xxyyy, g_0_0_z_0_zz_xxyyz, g_0_0_z_0_zz_xxyzz, g_0_0_z_0_zz_xxzzz, g_0_0_z_0_zz_xyyyy, g_0_0_z_0_zz_xyyyz, g_0_0_z_0_zz_xyyzz, g_0_0_z_0_zz_xyzzz, g_0_0_z_0_zz_xzzzz, g_0_0_z_0_zz_yyyyy, g_0_0_z_0_zz_yyyyz, g_0_0_z_0_zz_yyyzz, g_0_0_z_0_zz_yyzzz, g_0_0_z_0_zz_yzzzz, g_0_0_z_0_zz_zzzzz, g_z_0_z_0_zz_xxxxx, g_z_0_z_0_zz_xxxxxz, g_z_0_z_0_zz_xxxxy, g_z_0_z_0_zz_xxxxyz, g_z_0_z_0_zz_xxxxz, g_z_0_z_0_zz_xxxxzz, g_z_0_z_0_zz_xxxyy, g_z_0_z_0_zz_xxxyyz, g_z_0_z_0_zz_xxxyz, g_z_0_z_0_zz_xxxyzz, g_z_0_z_0_zz_xxxzz, g_z_0_z_0_zz_xxxzzz, g_z_0_z_0_zz_xxyyy, g_z_0_z_0_zz_xxyyyz, g_z_0_z_0_zz_xxyyz, g_z_0_z_0_zz_xxyyzz, g_z_0_z_0_zz_xxyzz, g_z_0_z_0_zz_xxyzzz, g_z_0_z_0_zz_xxzzz, g_z_0_z_0_zz_xxzzzz, g_z_0_z_0_zz_xyyyy, g_z_0_z_0_zz_xyyyyz, g_z_0_z_0_zz_xyyyz, g_z_0_z_0_zz_xyyyzz, g_z_0_z_0_zz_xyyzz, g_z_0_z_0_zz_xyyzzz, g_z_0_z_0_zz_xyzzz, g_z_0_z_0_zz_xyzzzz, g_z_0_z_0_zz_xzzzz, g_z_0_z_0_zz_xzzzzz, g_z_0_z_0_zz_yyyyy, g_z_0_z_0_zz_yyyyyz, g_z_0_z_0_zz_yyyyz, g_z_0_z_0_zz_yyyyzz, g_z_0_z_0_zz_yyyzz, g_z_0_z_0_zz_yyyzzz, g_z_0_z_0_zz_yyzzz, g_z_0_z_0_zz_yyzzzz, g_z_0_z_0_zz_yzzzz, g_z_0_z_0_zz_yzzzzz, g_z_0_z_0_zz_zzzzz, g_z_0_z_0_zz_zzzzzz, g_z_0_z_0_zzz_xxxxx, g_z_0_z_0_zzz_xxxxy, g_z_0_z_0_zzz_xxxxz, g_z_0_z_0_zzz_xxxyy, g_z_0_z_0_zzz_xxxyz, g_z_0_z_0_zzz_xxxzz, g_z_0_z_0_zzz_xxyyy, g_z_0_z_0_zzz_xxyyz, g_z_0_z_0_zzz_xxyzz, g_z_0_z_0_zzz_xxzzz, g_z_0_z_0_zzz_xyyyy, g_z_0_z_0_zzz_xyyyz, g_z_0_z_0_zzz_xyyzz, g_z_0_z_0_zzz_xyzzz, g_z_0_z_0_zzz_xzzzz, g_z_0_z_0_zzz_yyyyy, g_z_0_z_0_zzz_yyyyz, g_z_0_z_0_zzz_yyyzz, g_z_0_z_0_zzz_yyzzz, g_z_0_z_0_zzz_yzzzz, g_z_0_z_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzz_xxxxx[k] = -g_0_0_z_0_zz_xxxxx[k] - g_z_0_z_0_zz_xxxxx[k] * ab_z + g_z_0_z_0_zz_xxxxxz[k];

                g_z_0_z_0_zzz_xxxxy[k] = -g_0_0_z_0_zz_xxxxy[k] - g_z_0_z_0_zz_xxxxy[k] * ab_z + g_z_0_z_0_zz_xxxxyz[k];

                g_z_0_z_0_zzz_xxxxz[k] = -g_0_0_z_0_zz_xxxxz[k] - g_z_0_z_0_zz_xxxxz[k] * ab_z + g_z_0_z_0_zz_xxxxzz[k];

                g_z_0_z_0_zzz_xxxyy[k] = -g_0_0_z_0_zz_xxxyy[k] - g_z_0_z_0_zz_xxxyy[k] * ab_z + g_z_0_z_0_zz_xxxyyz[k];

                g_z_0_z_0_zzz_xxxyz[k] = -g_0_0_z_0_zz_xxxyz[k] - g_z_0_z_0_zz_xxxyz[k] * ab_z + g_z_0_z_0_zz_xxxyzz[k];

                g_z_0_z_0_zzz_xxxzz[k] = -g_0_0_z_0_zz_xxxzz[k] - g_z_0_z_0_zz_xxxzz[k] * ab_z + g_z_0_z_0_zz_xxxzzz[k];

                g_z_0_z_0_zzz_xxyyy[k] = -g_0_0_z_0_zz_xxyyy[k] - g_z_0_z_0_zz_xxyyy[k] * ab_z + g_z_0_z_0_zz_xxyyyz[k];

                g_z_0_z_0_zzz_xxyyz[k] = -g_0_0_z_0_zz_xxyyz[k] - g_z_0_z_0_zz_xxyyz[k] * ab_z + g_z_0_z_0_zz_xxyyzz[k];

                g_z_0_z_0_zzz_xxyzz[k] = -g_0_0_z_0_zz_xxyzz[k] - g_z_0_z_0_zz_xxyzz[k] * ab_z + g_z_0_z_0_zz_xxyzzz[k];

                g_z_0_z_0_zzz_xxzzz[k] = -g_0_0_z_0_zz_xxzzz[k] - g_z_0_z_0_zz_xxzzz[k] * ab_z + g_z_0_z_0_zz_xxzzzz[k];

                g_z_0_z_0_zzz_xyyyy[k] = -g_0_0_z_0_zz_xyyyy[k] - g_z_0_z_0_zz_xyyyy[k] * ab_z + g_z_0_z_0_zz_xyyyyz[k];

                g_z_0_z_0_zzz_xyyyz[k] = -g_0_0_z_0_zz_xyyyz[k] - g_z_0_z_0_zz_xyyyz[k] * ab_z + g_z_0_z_0_zz_xyyyzz[k];

                g_z_0_z_0_zzz_xyyzz[k] = -g_0_0_z_0_zz_xyyzz[k] - g_z_0_z_0_zz_xyyzz[k] * ab_z + g_z_0_z_0_zz_xyyzzz[k];

                g_z_0_z_0_zzz_xyzzz[k] = -g_0_0_z_0_zz_xyzzz[k] - g_z_0_z_0_zz_xyzzz[k] * ab_z + g_z_0_z_0_zz_xyzzzz[k];

                g_z_0_z_0_zzz_xzzzz[k] = -g_0_0_z_0_zz_xzzzz[k] - g_z_0_z_0_zz_xzzzz[k] * ab_z + g_z_0_z_0_zz_xzzzzz[k];

                g_z_0_z_0_zzz_yyyyy[k] = -g_0_0_z_0_zz_yyyyy[k] - g_z_0_z_0_zz_yyyyy[k] * ab_z + g_z_0_z_0_zz_yyyyyz[k];

                g_z_0_z_0_zzz_yyyyz[k] = -g_0_0_z_0_zz_yyyyz[k] - g_z_0_z_0_zz_yyyyz[k] * ab_z + g_z_0_z_0_zz_yyyyzz[k];

                g_z_0_z_0_zzz_yyyzz[k] = -g_0_0_z_0_zz_yyyzz[k] - g_z_0_z_0_zz_yyyzz[k] * ab_z + g_z_0_z_0_zz_yyyzzz[k];

                g_z_0_z_0_zzz_yyzzz[k] = -g_0_0_z_0_zz_yyzzz[k] - g_z_0_z_0_zz_yyzzz[k] * ab_z + g_z_0_z_0_zz_yyzzzz[k];

                g_z_0_z_0_zzz_yzzzz[k] = -g_0_0_z_0_zz_yzzzz[k] - g_z_0_z_0_zz_yzzzz[k] * ab_z + g_z_0_z_0_zz_yzzzzz[k];

                g_z_0_z_0_zzz_zzzzz[k] = -g_0_0_z_0_zz_zzzzz[k] - g_z_0_z_0_zz_zzzzz[k] * ab_z + g_z_0_z_0_zz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

