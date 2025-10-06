#include "ElectronRepulsionGeom2000ContrRecIDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_idxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_idxx,
                                            const size_t idx_geom_10_hdxx,
                                            const size_t idx_geom_20_hdxx,
                                            const size_t idx_geom_20_hfxx,
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
            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_10_off = idx_geom_10_hdxx + i * dcomps + j;

            auto g_x_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 269 * ccomps * dcomps);

            auto g_z_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 270 * ccomps * dcomps);

            auto g_z_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 271 * ccomps * dcomps);

            auto g_z_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 272 * ccomps * dcomps);

            auto g_z_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 273 * ccomps * dcomps);

            auto g_z_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 274 * ccomps * dcomps);

            auto g_z_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 275 * ccomps * dcomps);

            auto g_z_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 276 * ccomps * dcomps);

            auto g_z_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 277 * ccomps * dcomps);

            auto g_z_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 278 * ccomps * dcomps);

            auto g_z_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 279 * ccomps * dcomps);

            auto g_z_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 280 * ccomps * dcomps);

            auto g_z_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 281 * ccomps * dcomps);

            auto g_z_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 282 * ccomps * dcomps);

            auto g_z_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 283 * ccomps * dcomps);

            auto g_z_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 284 * ccomps * dcomps);

            auto g_z_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 285 * ccomps * dcomps);

            auto g_z_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 286 * ccomps * dcomps);

            auto g_z_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 287 * ccomps * dcomps);

            auto g_z_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 288 * ccomps * dcomps);

            auto g_z_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 289 * ccomps * dcomps);

            auto g_z_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 290 * ccomps * dcomps);

            auto g_z_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 291 * ccomps * dcomps);

            auto g_z_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 292 * ccomps * dcomps);

            auto g_z_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 293 * ccomps * dcomps);

            auto g_z_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 294 * ccomps * dcomps);

            auto g_z_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 295 * ccomps * dcomps);

            auto g_z_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 296 * ccomps * dcomps);

            auto g_z_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 297 * ccomps * dcomps);

            auto g_z_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 298 * ccomps * dcomps);

            auto g_z_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 299 * ccomps * dcomps);

            auto g_z_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 300 * ccomps * dcomps);

            auto g_z_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 301 * ccomps * dcomps);

            auto g_z_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 302 * ccomps * dcomps);

            auto g_z_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 303 * ccomps * dcomps);

            auto g_z_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 304 * ccomps * dcomps);

            auto g_z_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 305 * ccomps * dcomps);

            auto g_z_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 306 * ccomps * dcomps);

            auto g_z_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 307 * ccomps * dcomps);

            auto g_z_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 308 * ccomps * dcomps);

            auto g_z_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 309 * ccomps * dcomps);

            auto g_z_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 310 * ccomps * dcomps);

            auto g_z_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 311 * ccomps * dcomps);

            auto g_z_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 312 * ccomps * dcomps);

            auto g_z_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 313 * ccomps * dcomps);

            auto g_z_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 314 * ccomps * dcomps);

            auto g_z_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 315 * ccomps * dcomps);

            auto g_z_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 316 * ccomps * dcomps);

            auto g_z_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 317 * ccomps * dcomps);

            auto g_z_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 318 * ccomps * dcomps);

            auto g_z_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 319 * ccomps * dcomps);

            auto g_z_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 320 * ccomps * dcomps);

            auto g_z_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 321 * ccomps * dcomps);

            auto g_z_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 322 * ccomps * dcomps);

            auto g_z_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 323 * ccomps * dcomps);

            auto g_z_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 324 * ccomps * dcomps);

            auto g_z_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 325 * ccomps * dcomps);

            auto g_z_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 326 * ccomps * dcomps);

            auto g_z_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 327 * ccomps * dcomps);

            auto g_z_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 328 * ccomps * dcomps);

            auto g_z_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 329 * ccomps * dcomps);

            auto g_z_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 330 * ccomps * dcomps);

            auto g_z_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 331 * ccomps * dcomps);

            auto g_z_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 332 * ccomps * dcomps);

            auto g_z_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 333 * ccomps * dcomps);

            auto g_z_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 334 * ccomps * dcomps);

            auto g_z_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 335 * ccomps * dcomps);

            auto g_z_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 357 * ccomps * dcomps);

            auto g_z_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 358 * ccomps * dcomps);

            auto g_z_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 359 * ccomps * dcomps);

            auto g_z_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 360 * ccomps * dcomps);

            auto g_z_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 361 * ccomps * dcomps);

            auto g_z_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 362 * ccomps * dcomps);

            auto g_z_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 363 * ccomps * dcomps);

            auto g_z_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 377 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_20_off = idx_geom_20_hdxx + i * dcomps + j;

            auto g_xx_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 401 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 407 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 413 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 419 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 425 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 431 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 437 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 443 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 503 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 509 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 515 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 521 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 524 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 527 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 533 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 539 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 540 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 541 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 542 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 543 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 544 * ccomps * dcomps);

            auto g_yz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 545 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 546 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 547 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 548 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 549 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 550 * ccomps * dcomps);

            auto g_yz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 551 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 552 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 553 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 554 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 555 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 556 * ccomps * dcomps);

            auto g_yz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 557 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 558 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 559 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 560 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 561 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 562 * ccomps * dcomps);

            auto g_yz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 563 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 564 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 565 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 566 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 567 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 568 * ccomps * dcomps);

            auto g_yz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 569 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 570 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 571 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 572 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 573 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 574 * ccomps * dcomps);

            auto g_yz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 575 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 576 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 577 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 578 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 579 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 580 * ccomps * dcomps);

            auto g_yz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 581 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 582 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 583 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 584 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 585 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 586 * ccomps * dcomps);

            auto g_yz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 587 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 588 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 589 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 590 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 591 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 592 * ccomps * dcomps);

            auto g_yz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 593 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 594 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 595 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 596 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 597 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 598 * ccomps * dcomps);

            auto g_yz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 629 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 635 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 641 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 647 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 648 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 649 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 650 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 651 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 652 * ccomps * dcomps);

            auto g_zz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 653 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 654 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 655 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 656 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 657 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 658 * ccomps * dcomps);

            auto g_zz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 659 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 660 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 661 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 662 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 663 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 664 * ccomps * dcomps);

            auto g_zz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 665 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 666 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 667 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 668 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 669 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 670 * ccomps * dcomps);

            auto g_zz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 671 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 672 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 673 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 674 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 675 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 676 * ccomps * dcomps);

            auto g_zz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 677 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 678 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 679 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 680 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 681 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 682 * ccomps * dcomps);

            auto g_zz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 683 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 684 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 685 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 686 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 687 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 688 * ccomps * dcomps);

            auto g_zz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 689 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 690 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 691 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 692 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 693 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 694 * ccomps * dcomps);

            auto g_zz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 695 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 696 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 697 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 698 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 699 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 700 * ccomps * dcomps);

            auto g_zz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 701 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 702 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 703 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 704 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 705 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 706 * ccomps * dcomps);

            auto g_zz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 707 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 708 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 709 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 710 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 711 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 712 * ccomps * dcomps);

            auto g_zz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 713 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 714 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 715 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 716 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 717 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 718 * ccomps * dcomps);

            auto g_zz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 719 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 720 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 721 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 722 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 723 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 724 * ccomps * dcomps);

            auto g_zz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 725 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 726 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 727 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 728 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 729 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 730 * ccomps * dcomps);

            auto g_zz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 731 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 732 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 733 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 734 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 735 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 736 * ccomps * dcomps);

            auto g_zz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 737 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 738 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 739 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 740 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 741 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 742 * ccomps * dcomps);

            auto g_zz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 743 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 744 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 745 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 746 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 747 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 748 * ccomps * dcomps);

            auto g_zz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 755 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HFSS

            const auto hf_geom_20_off = idx_geom_20_hfxx + i * dcomps + j;

            auto g_xx_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 149 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 155 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 161 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 164 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 167 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 168 * ccomps * dcomps);

            auto g_xx_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 169 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 170 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 171 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 172 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 173 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 174 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 175 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 176 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 177 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 178 * ccomps * dcomps);

            auto g_xx_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 179 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 180 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 181 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 182 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 183 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 184 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 185 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 186 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 187 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 188 * ccomps * dcomps);

            auto g_xx_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 189 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 190 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 191 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 192 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 193 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 194 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 195 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 196 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 197 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 198 * ccomps * dcomps);

            auto g_xx_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 199 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 200 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 201 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 202 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 203 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 204 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 205 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 206 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 207 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 208 * ccomps * dcomps);

            auto g_xx_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 269 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 299 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 300 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 301 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 302 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 303 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 304 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 305 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 306 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 307 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 308 * ccomps * dcomps);

            auto g_xy_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 309 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 310 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 311 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 312 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 313 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 314 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 317 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 323 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 329 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 335 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 336 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 337 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 338 * ccomps * dcomps);

            auto g_xy_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 339 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 340 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 341 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 342 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 343 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 344 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 345 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 346 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 347 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 348 * ccomps * dcomps);

            auto g_xy_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 349 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 350 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 351 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 352 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 353 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 354 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 355 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 356 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 357 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 358 * ccomps * dcomps);

            auto g_xy_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 359 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 360 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 361 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 362 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 363 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 364 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 365 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 366 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 367 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 368 * ccomps * dcomps);

            auto g_xy_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 369 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 370 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 371 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 372 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 373 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 374 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 375 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 376 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 377 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 378 * ccomps * dcomps);

            auto g_xy_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 379 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 380 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 381 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 382 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 383 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 384 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 385 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 386 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 387 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 388 * ccomps * dcomps);

            auto g_xy_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 389 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 390 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 391 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 392 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 393 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 394 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 395 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 396 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 397 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 398 * ccomps * dcomps);

            auto g_xy_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 399 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 400 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 401 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 402 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 403 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 404 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 405 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 406 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 407 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 408 * ccomps * dcomps);

            auto g_xy_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 409 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 410 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 411 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 412 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 413 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 414 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 415 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 416 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 417 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 418 * ccomps * dcomps);

            auto g_xy_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 419 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 449 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 450 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 451 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 452 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 453 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 454 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 455 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 456 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 457 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 458 * ccomps * dcomps);

            auto g_xz_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 459 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 460 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 461 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 462 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 463 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 464 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 465 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 466 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 467 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 468 * ccomps * dcomps);

            auto g_xz_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 469 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 470 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 471 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 472 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 473 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 474 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 475 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 476 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 477 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 478 * ccomps * dcomps);

            auto g_xz_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 479 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 480 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 481 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 482 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 483 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 484 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 485 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 486 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 487 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 488 * ccomps * dcomps);

            auto g_xz_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 489 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 490 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 491 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 492 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 493 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 494 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 495 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 496 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 497 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 498 * ccomps * dcomps);

            auto g_xz_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 499 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 500 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 501 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 502 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 503 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 504 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 505 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 506 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 507 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 508 * ccomps * dcomps);

            auto g_xz_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 509 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 510 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 511 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 512 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 513 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 514 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 515 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 516 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 517 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 518 * ccomps * dcomps);

            auto g_xz_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 519 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 520 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 521 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 522 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 523 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 524 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 525 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 526 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 527 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 528 * ccomps * dcomps);

            auto g_xz_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 529 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 530 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 531 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 532 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 533 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 534 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 535 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 536 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 537 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 538 * ccomps * dcomps);

            auto g_xz_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 539 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 540 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 541 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 542 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 543 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 544 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 545 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 546 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 547 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 548 * ccomps * dcomps);

            auto g_xz_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 549 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 550 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 551 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 552 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 553 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 554 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 555 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 556 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 557 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 558 * ccomps * dcomps);

            auto g_xz_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 559 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 560 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 561 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 562 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 563 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 564 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 565 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 566 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 567 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 568 * ccomps * dcomps);

            auto g_xz_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 569 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 570 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 571 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 572 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 573 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 574 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 575 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 576 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 577 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 578 * ccomps * dcomps);

            auto g_xz_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 579 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 580 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 581 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 582 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 583 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 584 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 585 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 586 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 587 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 588 * ccomps * dcomps);

            auto g_xz_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 589 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 590 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 591 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 592 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 593 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 594 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 595 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 596 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 597 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 598 * ccomps * dcomps);

            auto g_xz_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 599 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 600 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 601 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 602 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 603 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 604 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 605 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 606 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 607 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 608 * ccomps * dcomps);

            auto g_xz_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 609 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 610 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 611 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 612 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 613 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 614 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 615 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 616 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 617 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 618 * ccomps * dcomps);

            auto g_xz_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 619 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 620 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 621 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 622 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 623 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 624 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 625 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 626 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 627 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 628 * ccomps * dcomps);

            auto g_xz_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 629 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 630 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 631 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 632 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 633 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 634 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 635 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 636 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 637 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 638 * ccomps * dcomps);

            auto g_yy_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 639 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 640 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 641 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 642 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 643 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 644 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 645 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 646 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 647 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 648 * ccomps * dcomps);

            auto g_yy_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 649 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 650 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 651 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 652 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 653 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 654 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 655 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 656 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 657 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 658 * ccomps * dcomps);

            auto g_yy_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 659 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 660 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 661 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 662 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 663 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 664 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 665 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 666 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 667 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 668 * ccomps * dcomps);

            auto g_yy_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 669 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 670 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 671 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 672 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 673 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 674 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 675 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 676 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 677 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 678 * ccomps * dcomps);

            auto g_yy_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 679 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 680 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 681 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 682 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 683 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 684 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 685 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 686 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 687 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 688 * ccomps * dcomps);

            auto g_yy_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 689 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 690 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 691 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 692 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 693 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 694 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 695 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 696 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 697 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 698 * ccomps * dcomps);

            auto g_yy_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 699 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 700 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 701 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 702 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 703 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 704 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 705 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 706 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 707 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 708 * ccomps * dcomps);

            auto g_yy_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 709 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 710 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 711 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 712 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 713 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 714 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 715 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 716 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 717 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 718 * ccomps * dcomps);

            auto g_yy_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 719 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 720 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 721 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 722 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 723 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 724 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 725 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 726 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 727 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 728 * ccomps * dcomps);

            auto g_yy_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 729 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 730 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 731 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 732 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 733 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 734 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 735 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 736 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 737 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 738 * ccomps * dcomps);

            auto g_yy_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 739 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 740 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 741 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 742 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 743 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 744 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 745 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 746 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 747 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 748 * ccomps * dcomps);

            auto g_yy_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 749 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 750 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 751 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 752 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 753 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 754 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 755 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 756 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 757 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 758 * ccomps * dcomps);

            auto g_yy_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 759 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 760 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 761 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 762 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 763 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 764 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 765 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 766 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 767 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 768 * ccomps * dcomps);

            auto g_yy_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 769 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 770 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 771 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 772 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 773 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 774 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 775 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 776 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 777 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 778 * ccomps * dcomps);

            auto g_yy_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 779 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 780 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 781 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 782 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 783 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 784 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 785 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 786 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 787 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 788 * ccomps * dcomps);

            auto g_yy_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 789 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 790 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 791 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 792 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 793 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 794 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 795 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 796 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 797 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 798 * ccomps * dcomps);

            auto g_yy_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 799 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 800 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 801 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 802 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 803 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 804 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 805 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 806 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 807 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 808 * ccomps * dcomps);

            auto g_yy_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 809 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 810 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 811 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 812 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 813 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 814 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 815 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 816 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 817 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 818 * ccomps * dcomps);

            auto g_yy_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 819 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 820 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 821 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 822 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 823 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 824 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 825 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 826 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 827 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 828 * ccomps * dcomps);

            auto g_yy_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 829 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 830 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 831 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 832 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 833 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 834 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 835 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 836 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 837 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 838 * ccomps * dcomps);

            auto g_yy_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 839 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 840 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 841 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 842 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 843 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 844 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 845 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 846 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 847 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 848 * ccomps * dcomps);

            auto g_yz_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 849 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 850 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 851 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 852 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 853 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 854 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 855 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 856 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 857 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 858 * ccomps * dcomps);

            auto g_yz_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 859 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 860 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 861 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 862 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 863 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 864 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 865 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 866 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 867 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 868 * ccomps * dcomps);

            auto g_yz_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 869 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 870 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 871 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 872 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 873 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 874 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 875 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 876 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 877 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 878 * ccomps * dcomps);

            auto g_yz_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 879 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 880 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 881 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 882 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 883 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 884 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 885 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 886 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 887 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 888 * ccomps * dcomps);

            auto g_yz_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 889 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 890 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 891 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 892 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 893 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 894 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 895 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 896 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 897 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 898 * ccomps * dcomps);

            auto g_yz_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 899 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 900 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 901 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 902 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 903 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 904 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 905 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 906 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 907 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 908 * ccomps * dcomps);

            auto g_yz_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 909 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 910 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 911 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 912 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 913 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 914 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 915 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 916 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 917 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 918 * ccomps * dcomps);

            auto g_yz_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 919 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 920 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 921 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 922 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 923 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 924 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 925 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 926 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 927 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 928 * ccomps * dcomps);

            auto g_yz_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 929 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 930 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 931 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 932 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 933 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 934 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 935 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 936 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 937 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 938 * ccomps * dcomps);

            auto g_yz_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 939 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 940 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 941 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 942 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 943 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 944 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 945 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 946 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 947 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 948 * ccomps * dcomps);

            auto g_yz_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 949 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 950 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 951 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 952 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 953 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 954 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 955 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 956 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 957 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 958 * ccomps * dcomps);

            auto g_yz_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 959 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 960 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 961 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 962 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 963 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 964 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 965 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 966 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 967 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 968 * ccomps * dcomps);

            auto g_yz_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 969 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 970 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 971 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 972 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 973 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 974 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 975 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 976 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 977 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 978 * ccomps * dcomps);

            auto g_yz_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 979 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 980 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 981 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 982 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 983 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 984 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 985 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 986 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 987 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 988 * ccomps * dcomps);

            auto g_yz_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 989 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 990 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 991 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 992 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 993 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 994 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 995 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 996 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 997 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 998 * ccomps * dcomps);

            auto g_yz_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 999 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 1000 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 1001 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 1002 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 1003 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 1004 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 1005 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 1006 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 1007 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 1008 * ccomps * dcomps);

            auto g_yz_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 1009 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 1010 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 1011 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 1012 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 1013 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 1014 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 1015 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 1016 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 1017 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 1018 * ccomps * dcomps);

            auto g_yz_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 1019 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 1020 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 1021 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 1022 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 1023 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 1024 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 1025 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 1026 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 1027 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 1028 * ccomps * dcomps);

            auto g_yz_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 1029 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 1030 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 1031 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 1032 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 1033 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 1034 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 1035 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 1036 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 1037 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 1038 * ccomps * dcomps);

            auto g_yz_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 1039 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 1040 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 1041 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 1042 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 1043 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 1044 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 1045 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 1046 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 1047 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 1048 * ccomps * dcomps);

            auto g_yz_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 1049 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxx = cbuffer.data(hf_geom_20_off + 1050 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxy = cbuffer.data(hf_geom_20_off + 1051 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxz = cbuffer.data(hf_geom_20_off + 1052 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xyy = cbuffer.data(hf_geom_20_off + 1053 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xyz = cbuffer.data(hf_geom_20_off + 1054 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xzz = cbuffer.data(hf_geom_20_off + 1055 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yyy = cbuffer.data(hf_geom_20_off + 1056 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yyz = cbuffer.data(hf_geom_20_off + 1057 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yzz = cbuffer.data(hf_geom_20_off + 1058 * ccomps * dcomps);

            auto g_zz_0_xxxxx_zzz = cbuffer.data(hf_geom_20_off + 1059 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxx = cbuffer.data(hf_geom_20_off + 1060 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxy = cbuffer.data(hf_geom_20_off + 1061 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxz = cbuffer.data(hf_geom_20_off + 1062 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xyy = cbuffer.data(hf_geom_20_off + 1063 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xyz = cbuffer.data(hf_geom_20_off + 1064 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xzz = cbuffer.data(hf_geom_20_off + 1065 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yyy = cbuffer.data(hf_geom_20_off + 1066 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yyz = cbuffer.data(hf_geom_20_off + 1067 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yzz = cbuffer.data(hf_geom_20_off + 1068 * ccomps * dcomps);

            auto g_zz_0_xxxxy_zzz = cbuffer.data(hf_geom_20_off + 1069 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxx = cbuffer.data(hf_geom_20_off + 1070 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxy = cbuffer.data(hf_geom_20_off + 1071 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxz = cbuffer.data(hf_geom_20_off + 1072 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xyy = cbuffer.data(hf_geom_20_off + 1073 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xyz = cbuffer.data(hf_geom_20_off + 1074 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xzz = cbuffer.data(hf_geom_20_off + 1075 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yyy = cbuffer.data(hf_geom_20_off + 1076 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yyz = cbuffer.data(hf_geom_20_off + 1077 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yzz = cbuffer.data(hf_geom_20_off + 1078 * ccomps * dcomps);

            auto g_zz_0_xxxxz_zzz = cbuffer.data(hf_geom_20_off + 1079 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxx = cbuffer.data(hf_geom_20_off + 1080 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxy = cbuffer.data(hf_geom_20_off + 1081 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxz = cbuffer.data(hf_geom_20_off + 1082 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xyy = cbuffer.data(hf_geom_20_off + 1083 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xyz = cbuffer.data(hf_geom_20_off + 1084 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xzz = cbuffer.data(hf_geom_20_off + 1085 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yyy = cbuffer.data(hf_geom_20_off + 1086 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yyz = cbuffer.data(hf_geom_20_off + 1087 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yzz = cbuffer.data(hf_geom_20_off + 1088 * ccomps * dcomps);

            auto g_zz_0_xxxyy_zzz = cbuffer.data(hf_geom_20_off + 1089 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxx = cbuffer.data(hf_geom_20_off + 1090 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxy = cbuffer.data(hf_geom_20_off + 1091 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxz = cbuffer.data(hf_geom_20_off + 1092 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xyy = cbuffer.data(hf_geom_20_off + 1093 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xyz = cbuffer.data(hf_geom_20_off + 1094 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xzz = cbuffer.data(hf_geom_20_off + 1095 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yyy = cbuffer.data(hf_geom_20_off + 1096 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yyz = cbuffer.data(hf_geom_20_off + 1097 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yzz = cbuffer.data(hf_geom_20_off + 1098 * ccomps * dcomps);

            auto g_zz_0_xxxyz_zzz = cbuffer.data(hf_geom_20_off + 1099 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxx = cbuffer.data(hf_geom_20_off + 1100 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxy = cbuffer.data(hf_geom_20_off + 1101 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxz = cbuffer.data(hf_geom_20_off + 1102 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xyy = cbuffer.data(hf_geom_20_off + 1103 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xyz = cbuffer.data(hf_geom_20_off + 1104 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xzz = cbuffer.data(hf_geom_20_off + 1105 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yyy = cbuffer.data(hf_geom_20_off + 1106 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yyz = cbuffer.data(hf_geom_20_off + 1107 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yzz = cbuffer.data(hf_geom_20_off + 1108 * ccomps * dcomps);

            auto g_zz_0_xxxzz_zzz = cbuffer.data(hf_geom_20_off + 1109 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxx = cbuffer.data(hf_geom_20_off + 1110 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxy = cbuffer.data(hf_geom_20_off + 1111 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxz = cbuffer.data(hf_geom_20_off + 1112 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xyy = cbuffer.data(hf_geom_20_off + 1113 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xyz = cbuffer.data(hf_geom_20_off + 1114 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xzz = cbuffer.data(hf_geom_20_off + 1115 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yyy = cbuffer.data(hf_geom_20_off + 1116 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yyz = cbuffer.data(hf_geom_20_off + 1117 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yzz = cbuffer.data(hf_geom_20_off + 1118 * ccomps * dcomps);

            auto g_zz_0_xxyyy_zzz = cbuffer.data(hf_geom_20_off + 1119 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxx = cbuffer.data(hf_geom_20_off + 1120 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxy = cbuffer.data(hf_geom_20_off + 1121 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxz = cbuffer.data(hf_geom_20_off + 1122 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xyy = cbuffer.data(hf_geom_20_off + 1123 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xyz = cbuffer.data(hf_geom_20_off + 1124 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xzz = cbuffer.data(hf_geom_20_off + 1125 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yyy = cbuffer.data(hf_geom_20_off + 1126 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yyz = cbuffer.data(hf_geom_20_off + 1127 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yzz = cbuffer.data(hf_geom_20_off + 1128 * ccomps * dcomps);

            auto g_zz_0_xxyyz_zzz = cbuffer.data(hf_geom_20_off + 1129 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxx = cbuffer.data(hf_geom_20_off + 1130 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxy = cbuffer.data(hf_geom_20_off + 1131 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxz = cbuffer.data(hf_geom_20_off + 1132 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xyy = cbuffer.data(hf_geom_20_off + 1133 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xyz = cbuffer.data(hf_geom_20_off + 1134 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xzz = cbuffer.data(hf_geom_20_off + 1135 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yyy = cbuffer.data(hf_geom_20_off + 1136 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yyz = cbuffer.data(hf_geom_20_off + 1137 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yzz = cbuffer.data(hf_geom_20_off + 1138 * ccomps * dcomps);

            auto g_zz_0_xxyzz_zzz = cbuffer.data(hf_geom_20_off + 1139 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxx = cbuffer.data(hf_geom_20_off + 1140 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxy = cbuffer.data(hf_geom_20_off + 1141 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxz = cbuffer.data(hf_geom_20_off + 1142 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xyy = cbuffer.data(hf_geom_20_off + 1143 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xyz = cbuffer.data(hf_geom_20_off + 1144 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xzz = cbuffer.data(hf_geom_20_off + 1145 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yyy = cbuffer.data(hf_geom_20_off + 1146 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yyz = cbuffer.data(hf_geom_20_off + 1147 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yzz = cbuffer.data(hf_geom_20_off + 1148 * ccomps * dcomps);

            auto g_zz_0_xxzzz_zzz = cbuffer.data(hf_geom_20_off + 1149 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxx = cbuffer.data(hf_geom_20_off + 1150 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxy = cbuffer.data(hf_geom_20_off + 1151 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxz = cbuffer.data(hf_geom_20_off + 1152 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xyy = cbuffer.data(hf_geom_20_off + 1153 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xyz = cbuffer.data(hf_geom_20_off + 1154 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xzz = cbuffer.data(hf_geom_20_off + 1155 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yyy = cbuffer.data(hf_geom_20_off + 1156 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yyz = cbuffer.data(hf_geom_20_off + 1157 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yzz = cbuffer.data(hf_geom_20_off + 1158 * ccomps * dcomps);

            auto g_zz_0_xyyyy_zzz = cbuffer.data(hf_geom_20_off + 1159 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxx = cbuffer.data(hf_geom_20_off + 1160 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxy = cbuffer.data(hf_geom_20_off + 1161 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxz = cbuffer.data(hf_geom_20_off + 1162 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xyy = cbuffer.data(hf_geom_20_off + 1163 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xyz = cbuffer.data(hf_geom_20_off + 1164 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xzz = cbuffer.data(hf_geom_20_off + 1165 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yyy = cbuffer.data(hf_geom_20_off + 1166 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yyz = cbuffer.data(hf_geom_20_off + 1167 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yzz = cbuffer.data(hf_geom_20_off + 1168 * ccomps * dcomps);

            auto g_zz_0_xyyyz_zzz = cbuffer.data(hf_geom_20_off + 1169 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxx = cbuffer.data(hf_geom_20_off + 1170 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxy = cbuffer.data(hf_geom_20_off + 1171 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxz = cbuffer.data(hf_geom_20_off + 1172 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xyy = cbuffer.data(hf_geom_20_off + 1173 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xyz = cbuffer.data(hf_geom_20_off + 1174 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xzz = cbuffer.data(hf_geom_20_off + 1175 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yyy = cbuffer.data(hf_geom_20_off + 1176 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yyz = cbuffer.data(hf_geom_20_off + 1177 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yzz = cbuffer.data(hf_geom_20_off + 1178 * ccomps * dcomps);

            auto g_zz_0_xyyzz_zzz = cbuffer.data(hf_geom_20_off + 1179 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxx = cbuffer.data(hf_geom_20_off + 1180 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxy = cbuffer.data(hf_geom_20_off + 1181 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxz = cbuffer.data(hf_geom_20_off + 1182 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xyy = cbuffer.data(hf_geom_20_off + 1183 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xyz = cbuffer.data(hf_geom_20_off + 1184 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xzz = cbuffer.data(hf_geom_20_off + 1185 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yyy = cbuffer.data(hf_geom_20_off + 1186 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yyz = cbuffer.data(hf_geom_20_off + 1187 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yzz = cbuffer.data(hf_geom_20_off + 1188 * ccomps * dcomps);

            auto g_zz_0_xyzzz_zzz = cbuffer.data(hf_geom_20_off + 1189 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxx = cbuffer.data(hf_geom_20_off + 1190 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxy = cbuffer.data(hf_geom_20_off + 1191 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxz = cbuffer.data(hf_geom_20_off + 1192 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xyy = cbuffer.data(hf_geom_20_off + 1193 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xyz = cbuffer.data(hf_geom_20_off + 1194 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xzz = cbuffer.data(hf_geom_20_off + 1195 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yyy = cbuffer.data(hf_geom_20_off + 1196 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yyz = cbuffer.data(hf_geom_20_off + 1197 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yzz = cbuffer.data(hf_geom_20_off + 1198 * ccomps * dcomps);

            auto g_zz_0_xzzzz_zzz = cbuffer.data(hf_geom_20_off + 1199 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxx = cbuffer.data(hf_geom_20_off + 1200 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxy = cbuffer.data(hf_geom_20_off + 1201 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxz = cbuffer.data(hf_geom_20_off + 1202 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xyy = cbuffer.data(hf_geom_20_off + 1203 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xyz = cbuffer.data(hf_geom_20_off + 1204 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xzz = cbuffer.data(hf_geom_20_off + 1205 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yyy = cbuffer.data(hf_geom_20_off + 1206 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yyz = cbuffer.data(hf_geom_20_off + 1207 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yzz = cbuffer.data(hf_geom_20_off + 1208 * ccomps * dcomps);

            auto g_zz_0_yyyyy_zzz = cbuffer.data(hf_geom_20_off + 1209 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxx = cbuffer.data(hf_geom_20_off + 1210 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxy = cbuffer.data(hf_geom_20_off + 1211 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxz = cbuffer.data(hf_geom_20_off + 1212 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xyy = cbuffer.data(hf_geom_20_off + 1213 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xyz = cbuffer.data(hf_geom_20_off + 1214 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xzz = cbuffer.data(hf_geom_20_off + 1215 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yyy = cbuffer.data(hf_geom_20_off + 1216 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yyz = cbuffer.data(hf_geom_20_off + 1217 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yzz = cbuffer.data(hf_geom_20_off + 1218 * ccomps * dcomps);

            auto g_zz_0_yyyyz_zzz = cbuffer.data(hf_geom_20_off + 1219 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxx = cbuffer.data(hf_geom_20_off + 1220 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxy = cbuffer.data(hf_geom_20_off + 1221 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxz = cbuffer.data(hf_geom_20_off + 1222 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xyy = cbuffer.data(hf_geom_20_off + 1223 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xyz = cbuffer.data(hf_geom_20_off + 1224 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xzz = cbuffer.data(hf_geom_20_off + 1225 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yyy = cbuffer.data(hf_geom_20_off + 1226 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yyz = cbuffer.data(hf_geom_20_off + 1227 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yzz = cbuffer.data(hf_geom_20_off + 1228 * ccomps * dcomps);

            auto g_zz_0_yyyzz_zzz = cbuffer.data(hf_geom_20_off + 1229 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxx = cbuffer.data(hf_geom_20_off + 1230 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxy = cbuffer.data(hf_geom_20_off + 1231 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxz = cbuffer.data(hf_geom_20_off + 1232 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xyy = cbuffer.data(hf_geom_20_off + 1233 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xyz = cbuffer.data(hf_geom_20_off + 1234 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xzz = cbuffer.data(hf_geom_20_off + 1235 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yyy = cbuffer.data(hf_geom_20_off + 1236 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yyz = cbuffer.data(hf_geom_20_off + 1237 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yzz = cbuffer.data(hf_geom_20_off + 1238 * ccomps * dcomps);

            auto g_zz_0_yyzzz_zzz = cbuffer.data(hf_geom_20_off + 1239 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxx = cbuffer.data(hf_geom_20_off + 1240 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxy = cbuffer.data(hf_geom_20_off + 1241 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxz = cbuffer.data(hf_geom_20_off + 1242 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xyy = cbuffer.data(hf_geom_20_off + 1243 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xyz = cbuffer.data(hf_geom_20_off + 1244 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xzz = cbuffer.data(hf_geom_20_off + 1245 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yyy = cbuffer.data(hf_geom_20_off + 1246 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yyz = cbuffer.data(hf_geom_20_off + 1247 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yzz = cbuffer.data(hf_geom_20_off + 1248 * ccomps * dcomps);

            auto g_zz_0_yzzzz_zzz = cbuffer.data(hf_geom_20_off + 1249 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxx = cbuffer.data(hf_geom_20_off + 1250 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxy = cbuffer.data(hf_geom_20_off + 1251 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxz = cbuffer.data(hf_geom_20_off + 1252 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xyy = cbuffer.data(hf_geom_20_off + 1253 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xyz = cbuffer.data(hf_geom_20_off + 1254 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xzz = cbuffer.data(hf_geom_20_off + 1255 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yyy = cbuffer.data(hf_geom_20_off + 1256 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yyz = cbuffer.data(hf_geom_20_off + 1257 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yzz = cbuffer.data(hf_geom_20_off + 1258 * ccomps * dcomps);

            auto g_zz_0_zzzzz_zzz = cbuffer.data(hf_geom_20_off + 1259 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_idxx

            const auto id_geom_20_off = idx_geom_20_idxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xx, g_x_0_xxxxx_xy, g_x_0_xxxxx_xz, g_x_0_xxxxx_yy, g_x_0_xxxxx_yz, g_x_0_xxxxx_zz, g_xx_0_xxxxx_xx, g_xx_0_xxxxx_xxx, g_xx_0_xxxxx_xxy, g_xx_0_xxxxx_xxz, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_xyy, g_xx_0_xxxxx_xyz, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_xzz, g_xx_0_xxxxx_yy, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_zz, g_xx_0_xxxxxx_xx, g_xx_0_xxxxxx_xy, g_xx_0_xxxxxx_xz, g_xx_0_xxxxxx_yy, g_xx_0_xxxxxx_yz, g_xx_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxx_xx[k] = -2.0 * g_x_0_xxxxx_xx[k] - g_xx_0_xxxxx_xx[k] * ab_x + g_xx_0_xxxxx_xxx[k];

                g_xx_0_xxxxxx_xy[k] = -2.0 * g_x_0_xxxxx_xy[k] - g_xx_0_xxxxx_xy[k] * ab_x + g_xx_0_xxxxx_xxy[k];

                g_xx_0_xxxxxx_xz[k] = -2.0 * g_x_0_xxxxx_xz[k] - g_xx_0_xxxxx_xz[k] * ab_x + g_xx_0_xxxxx_xxz[k];

                g_xx_0_xxxxxx_yy[k] = -2.0 * g_x_0_xxxxx_yy[k] - g_xx_0_xxxxx_yy[k] * ab_x + g_xx_0_xxxxx_xyy[k];

                g_xx_0_xxxxxx_yz[k] = -2.0 * g_x_0_xxxxx_yz[k] - g_xx_0_xxxxx_yz[k] * ab_x + g_xx_0_xxxxx_xyz[k];

                g_xx_0_xxxxxx_zz[k] = -2.0 * g_x_0_xxxxx_zz[k] - g_xx_0_xxxxx_zz[k] * ab_x + g_xx_0_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_xx, g_xx_0_xxxxx_xxy, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_xyy, g_xx_0_xxxxx_xyz, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_yy, g_xx_0_xxxxx_yyy, g_xx_0_xxxxx_yyz, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_yzz, g_xx_0_xxxxx_zz, g_xx_0_xxxxxy_xx, g_xx_0_xxxxxy_xy, g_xx_0_xxxxxy_xz, g_xx_0_xxxxxy_yy, g_xx_0_xxxxxy_yz, g_xx_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxy_xx[k] = -g_xx_0_xxxxx_xx[k] * ab_y + g_xx_0_xxxxx_xxy[k];

                g_xx_0_xxxxxy_xy[k] = -g_xx_0_xxxxx_xy[k] * ab_y + g_xx_0_xxxxx_xyy[k];

                g_xx_0_xxxxxy_xz[k] = -g_xx_0_xxxxx_xz[k] * ab_y + g_xx_0_xxxxx_xyz[k];

                g_xx_0_xxxxxy_yy[k] = -g_xx_0_xxxxx_yy[k] * ab_y + g_xx_0_xxxxx_yyy[k];

                g_xx_0_xxxxxy_yz[k] = -g_xx_0_xxxxx_yz[k] * ab_y + g_xx_0_xxxxx_yyz[k];

                g_xx_0_xxxxxy_zz[k] = -g_xx_0_xxxxx_zz[k] * ab_y + g_xx_0_xxxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_xx, g_xx_0_xxxxx_xxz, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_xyz, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_xzz, g_xx_0_xxxxx_yy, g_xx_0_xxxxx_yyz, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_yzz, g_xx_0_xxxxx_zz, g_xx_0_xxxxx_zzz, g_xx_0_xxxxxz_xx, g_xx_0_xxxxxz_xy, g_xx_0_xxxxxz_xz, g_xx_0_xxxxxz_yy, g_xx_0_xxxxxz_yz, g_xx_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxz_xx[k] = -g_xx_0_xxxxx_xx[k] * ab_z + g_xx_0_xxxxx_xxz[k];

                g_xx_0_xxxxxz_xy[k] = -g_xx_0_xxxxx_xy[k] * ab_z + g_xx_0_xxxxx_xyz[k];

                g_xx_0_xxxxxz_xz[k] = -g_xx_0_xxxxx_xz[k] * ab_z + g_xx_0_xxxxx_xzz[k];

                g_xx_0_xxxxxz_yy[k] = -g_xx_0_xxxxx_yy[k] * ab_z + g_xx_0_xxxxx_yyz[k];

                g_xx_0_xxxxxz_yz[k] = -g_xx_0_xxxxx_yz[k] * ab_z + g_xx_0_xxxxx_yzz[k];

                g_xx_0_xxxxxz_zz[k] = -g_xx_0_xxxxx_zz[k] * ab_z + g_xx_0_xxxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxy_xx, g_xx_0_xxxxy_xxy, g_xx_0_xxxxy_xy, g_xx_0_xxxxy_xyy, g_xx_0_xxxxy_xyz, g_xx_0_xxxxy_xz, g_xx_0_xxxxy_yy, g_xx_0_xxxxy_yyy, g_xx_0_xxxxy_yyz, g_xx_0_xxxxy_yz, g_xx_0_xxxxy_yzz, g_xx_0_xxxxy_zz, g_xx_0_xxxxyy_xx, g_xx_0_xxxxyy_xy, g_xx_0_xxxxyy_xz, g_xx_0_xxxxyy_yy, g_xx_0_xxxxyy_yz, g_xx_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyy_xx[k] = -g_xx_0_xxxxy_xx[k] * ab_y + g_xx_0_xxxxy_xxy[k];

                g_xx_0_xxxxyy_xy[k] = -g_xx_0_xxxxy_xy[k] * ab_y + g_xx_0_xxxxy_xyy[k];

                g_xx_0_xxxxyy_xz[k] = -g_xx_0_xxxxy_xz[k] * ab_y + g_xx_0_xxxxy_xyz[k];

                g_xx_0_xxxxyy_yy[k] = -g_xx_0_xxxxy_yy[k] * ab_y + g_xx_0_xxxxy_yyy[k];

                g_xx_0_xxxxyy_yz[k] = -g_xx_0_xxxxy_yz[k] * ab_y + g_xx_0_xxxxy_yyz[k];

                g_xx_0_xxxxyy_zz[k] = -g_xx_0_xxxxy_zz[k] * ab_y + g_xx_0_xxxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxyz_xx, g_xx_0_xxxxyz_xy, g_xx_0_xxxxyz_xz, g_xx_0_xxxxyz_yy, g_xx_0_xxxxyz_yz, g_xx_0_xxxxyz_zz, g_xx_0_xxxxz_xx, g_xx_0_xxxxz_xxy, g_xx_0_xxxxz_xy, g_xx_0_xxxxz_xyy, g_xx_0_xxxxz_xyz, g_xx_0_xxxxz_xz, g_xx_0_xxxxz_yy, g_xx_0_xxxxz_yyy, g_xx_0_xxxxz_yyz, g_xx_0_xxxxz_yz, g_xx_0_xxxxz_yzz, g_xx_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyz_xx[k] = -g_xx_0_xxxxz_xx[k] * ab_y + g_xx_0_xxxxz_xxy[k];

                g_xx_0_xxxxyz_xy[k] = -g_xx_0_xxxxz_xy[k] * ab_y + g_xx_0_xxxxz_xyy[k];

                g_xx_0_xxxxyz_xz[k] = -g_xx_0_xxxxz_xz[k] * ab_y + g_xx_0_xxxxz_xyz[k];

                g_xx_0_xxxxyz_yy[k] = -g_xx_0_xxxxz_yy[k] * ab_y + g_xx_0_xxxxz_yyy[k];

                g_xx_0_xxxxyz_yz[k] = -g_xx_0_xxxxz_yz[k] * ab_y + g_xx_0_xxxxz_yyz[k];

                g_xx_0_xxxxyz_zz[k] = -g_xx_0_xxxxz_zz[k] * ab_y + g_xx_0_xxxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxz_xx, g_xx_0_xxxxz_xxz, g_xx_0_xxxxz_xy, g_xx_0_xxxxz_xyz, g_xx_0_xxxxz_xz, g_xx_0_xxxxz_xzz, g_xx_0_xxxxz_yy, g_xx_0_xxxxz_yyz, g_xx_0_xxxxz_yz, g_xx_0_xxxxz_yzz, g_xx_0_xxxxz_zz, g_xx_0_xxxxz_zzz, g_xx_0_xxxxzz_xx, g_xx_0_xxxxzz_xy, g_xx_0_xxxxzz_xz, g_xx_0_xxxxzz_yy, g_xx_0_xxxxzz_yz, g_xx_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxzz_xx[k] = -g_xx_0_xxxxz_xx[k] * ab_z + g_xx_0_xxxxz_xxz[k];

                g_xx_0_xxxxzz_xy[k] = -g_xx_0_xxxxz_xy[k] * ab_z + g_xx_0_xxxxz_xyz[k];

                g_xx_0_xxxxzz_xz[k] = -g_xx_0_xxxxz_xz[k] * ab_z + g_xx_0_xxxxz_xzz[k];

                g_xx_0_xxxxzz_yy[k] = -g_xx_0_xxxxz_yy[k] * ab_z + g_xx_0_xxxxz_yyz[k];

                g_xx_0_xxxxzz_yz[k] = -g_xx_0_xxxxz_yz[k] * ab_z + g_xx_0_xxxxz_yzz[k];

                g_xx_0_xxxxzz_zz[k] = -g_xx_0_xxxxz_zz[k] * ab_z + g_xx_0_xxxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyy_xx, g_xx_0_xxxyy_xxy, g_xx_0_xxxyy_xy, g_xx_0_xxxyy_xyy, g_xx_0_xxxyy_xyz, g_xx_0_xxxyy_xz, g_xx_0_xxxyy_yy, g_xx_0_xxxyy_yyy, g_xx_0_xxxyy_yyz, g_xx_0_xxxyy_yz, g_xx_0_xxxyy_yzz, g_xx_0_xxxyy_zz, g_xx_0_xxxyyy_xx, g_xx_0_xxxyyy_xy, g_xx_0_xxxyyy_xz, g_xx_0_xxxyyy_yy, g_xx_0_xxxyyy_yz, g_xx_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyy_xx[k] = -g_xx_0_xxxyy_xx[k] * ab_y + g_xx_0_xxxyy_xxy[k];

                g_xx_0_xxxyyy_xy[k] = -g_xx_0_xxxyy_xy[k] * ab_y + g_xx_0_xxxyy_xyy[k];

                g_xx_0_xxxyyy_xz[k] = -g_xx_0_xxxyy_xz[k] * ab_y + g_xx_0_xxxyy_xyz[k];

                g_xx_0_xxxyyy_yy[k] = -g_xx_0_xxxyy_yy[k] * ab_y + g_xx_0_xxxyy_yyy[k];

                g_xx_0_xxxyyy_yz[k] = -g_xx_0_xxxyy_yz[k] * ab_y + g_xx_0_xxxyy_yyz[k];

                g_xx_0_xxxyyy_zz[k] = -g_xx_0_xxxyy_zz[k] * ab_y + g_xx_0_xxxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyyz_xx, g_xx_0_xxxyyz_xy, g_xx_0_xxxyyz_xz, g_xx_0_xxxyyz_yy, g_xx_0_xxxyyz_yz, g_xx_0_xxxyyz_zz, g_xx_0_xxxyz_xx, g_xx_0_xxxyz_xxy, g_xx_0_xxxyz_xy, g_xx_0_xxxyz_xyy, g_xx_0_xxxyz_xyz, g_xx_0_xxxyz_xz, g_xx_0_xxxyz_yy, g_xx_0_xxxyz_yyy, g_xx_0_xxxyz_yyz, g_xx_0_xxxyz_yz, g_xx_0_xxxyz_yzz, g_xx_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyz_xx[k] = -g_xx_0_xxxyz_xx[k] * ab_y + g_xx_0_xxxyz_xxy[k];

                g_xx_0_xxxyyz_xy[k] = -g_xx_0_xxxyz_xy[k] * ab_y + g_xx_0_xxxyz_xyy[k];

                g_xx_0_xxxyyz_xz[k] = -g_xx_0_xxxyz_xz[k] * ab_y + g_xx_0_xxxyz_xyz[k];

                g_xx_0_xxxyyz_yy[k] = -g_xx_0_xxxyz_yy[k] * ab_y + g_xx_0_xxxyz_yyy[k];

                g_xx_0_xxxyyz_yz[k] = -g_xx_0_xxxyz_yz[k] * ab_y + g_xx_0_xxxyz_yyz[k];

                g_xx_0_xxxyyz_zz[k] = -g_xx_0_xxxyz_zz[k] * ab_y + g_xx_0_xxxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyzz_xx, g_xx_0_xxxyzz_xy, g_xx_0_xxxyzz_xz, g_xx_0_xxxyzz_yy, g_xx_0_xxxyzz_yz, g_xx_0_xxxyzz_zz, g_xx_0_xxxzz_xx, g_xx_0_xxxzz_xxy, g_xx_0_xxxzz_xy, g_xx_0_xxxzz_xyy, g_xx_0_xxxzz_xyz, g_xx_0_xxxzz_xz, g_xx_0_xxxzz_yy, g_xx_0_xxxzz_yyy, g_xx_0_xxxzz_yyz, g_xx_0_xxxzz_yz, g_xx_0_xxxzz_yzz, g_xx_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyzz_xx[k] = -g_xx_0_xxxzz_xx[k] * ab_y + g_xx_0_xxxzz_xxy[k];

                g_xx_0_xxxyzz_xy[k] = -g_xx_0_xxxzz_xy[k] * ab_y + g_xx_0_xxxzz_xyy[k];

                g_xx_0_xxxyzz_xz[k] = -g_xx_0_xxxzz_xz[k] * ab_y + g_xx_0_xxxzz_xyz[k];

                g_xx_0_xxxyzz_yy[k] = -g_xx_0_xxxzz_yy[k] * ab_y + g_xx_0_xxxzz_yyy[k];

                g_xx_0_xxxyzz_yz[k] = -g_xx_0_xxxzz_yz[k] * ab_y + g_xx_0_xxxzz_yyz[k];

                g_xx_0_xxxyzz_zz[k] = -g_xx_0_xxxzz_zz[k] * ab_y + g_xx_0_xxxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxzz_xx, g_xx_0_xxxzz_xxz, g_xx_0_xxxzz_xy, g_xx_0_xxxzz_xyz, g_xx_0_xxxzz_xz, g_xx_0_xxxzz_xzz, g_xx_0_xxxzz_yy, g_xx_0_xxxzz_yyz, g_xx_0_xxxzz_yz, g_xx_0_xxxzz_yzz, g_xx_0_xxxzz_zz, g_xx_0_xxxzz_zzz, g_xx_0_xxxzzz_xx, g_xx_0_xxxzzz_xy, g_xx_0_xxxzzz_xz, g_xx_0_xxxzzz_yy, g_xx_0_xxxzzz_yz, g_xx_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzzz_xx[k] = -g_xx_0_xxxzz_xx[k] * ab_z + g_xx_0_xxxzz_xxz[k];

                g_xx_0_xxxzzz_xy[k] = -g_xx_0_xxxzz_xy[k] * ab_z + g_xx_0_xxxzz_xyz[k];

                g_xx_0_xxxzzz_xz[k] = -g_xx_0_xxxzz_xz[k] * ab_z + g_xx_0_xxxzz_xzz[k];

                g_xx_0_xxxzzz_yy[k] = -g_xx_0_xxxzz_yy[k] * ab_z + g_xx_0_xxxzz_yyz[k];

                g_xx_0_xxxzzz_yz[k] = -g_xx_0_xxxzz_yz[k] * ab_z + g_xx_0_xxxzz_yzz[k];

                g_xx_0_xxxzzz_zz[k] = -g_xx_0_xxxzz_zz[k] * ab_z + g_xx_0_xxxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyy_xx, g_xx_0_xxyyy_xxy, g_xx_0_xxyyy_xy, g_xx_0_xxyyy_xyy, g_xx_0_xxyyy_xyz, g_xx_0_xxyyy_xz, g_xx_0_xxyyy_yy, g_xx_0_xxyyy_yyy, g_xx_0_xxyyy_yyz, g_xx_0_xxyyy_yz, g_xx_0_xxyyy_yzz, g_xx_0_xxyyy_zz, g_xx_0_xxyyyy_xx, g_xx_0_xxyyyy_xy, g_xx_0_xxyyyy_xz, g_xx_0_xxyyyy_yy, g_xx_0_xxyyyy_yz, g_xx_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyy_xx[k] = -g_xx_0_xxyyy_xx[k] * ab_y + g_xx_0_xxyyy_xxy[k];

                g_xx_0_xxyyyy_xy[k] = -g_xx_0_xxyyy_xy[k] * ab_y + g_xx_0_xxyyy_xyy[k];

                g_xx_0_xxyyyy_xz[k] = -g_xx_0_xxyyy_xz[k] * ab_y + g_xx_0_xxyyy_xyz[k];

                g_xx_0_xxyyyy_yy[k] = -g_xx_0_xxyyy_yy[k] * ab_y + g_xx_0_xxyyy_yyy[k];

                g_xx_0_xxyyyy_yz[k] = -g_xx_0_xxyyy_yz[k] * ab_y + g_xx_0_xxyyy_yyz[k];

                g_xx_0_xxyyyy_zz[k] = -g_xx_0_xxyyy_zz[k] * ab_y + g_xx_0_xxyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyyz_xx, g_xx_0_xxyyyz_xy, g_xx_0_xxyyyz_xz, g_xx_0_xxyyyz_yy, g_xx_0_xxyyyz_yz, g_xx_0_xxyyyz_zz, g_xx_0_xxyyz_xx, g_xx_0_xxyyz_xxy, g_xx_0_xxyyz_xy, g_xx_0_xxyyz_xyy, g_xx_0_xxyyz_xyz, g_xx_0_xxyyz_xz, g_xx_0_xxyyz_yy, g_xx_0_xxyyz_yyy, g_xx_0_xxyyz_yyz, g_xx_0_xxyyz_yz, g_xx_0_xxyyz_yzz, g_xx_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyz_xx[k] = -g_xx_0_xxyyz_xx[k] * ab_y + g_xx_0_xxyyz_xxy[k];

                g_xx_0_xxyyyz_xy[k] = -g_xx_0_xxyyz_xy[k] * ab_y + g_xx_0_xxyyz_xyy[k];

                g_xx_0_xxyyyz_xz[k] = -g_xx_0_xxyyz_xz[k] * ab_y + g_xx_0_xxyyz_xyz[k];

                g_xx_0_xxyyyz_yy[k] = -g_xx_0_xxyyz_yy[k] * ab_y + g_xx_0_xxyyz_yyy[k];

                g_xx_0_xxyyyz_yz[k] = -g_xx_0_xxyyz_yz[k] * ab_y + g_xx_0_xxyyz_yyz[k];

                g_xx_0_xxyyyz_zz[k] = -g_xx_0_xxyyz_zz[k] * ab_y + g_xx_0_xxyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyzz_xx, g_xx_0_xxyyzz_xy, g_xx_0_xxyyzz_xz, g_xx_0_xxyyzz_yy, g_xx_0_xxyyzz_yz, g_xx_0_xxyyzz_zz, g_xx_0_xxyzz_xx, g_xx_0_xxyzz_xxy, g_xx_0_xxyzz_xy, g_xx_0_xxyzz_xyy, g_xx_0_xxyzz_xyz, g_xx_0_xxyzz_xz, g_xx_0_xxyzz_yy, g_xx_0_xxyzz_yyy, g_xx_0_xxyzz_yyz, g_xx_0_xxyzz_yz, g_xx_0_xxyzz_yzz, g_xx_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyzz_xx[k] = -g_xx_0_xxyzz_xx[k] * ab_y + g_xx_0_xxyzz_xxy[k];

                g_xx_0_xxyyzz_xy[k] = -g_xx_0_xxyzz_xy[k] * ab_y + g_xx_0_xxyzz_xyy[k];

                g_xx_0_xxyyzz_xz[k] = -g_xx_0_xxyzz_xz[k] * ab_y + g_xx_0_xxyzz_xyz[k];

                g_xx_0_xxyyzz_yy[k] = -g_xx_0_xxyzz_yy[k] * ab_y + g_xx_0_xxyzz_yyy[k];

                g_xx_0_xxyyzz_yz[k] = -g_xx_0_xxyzz_yz[k] * ab_y + g_xx_0_xxyzz_yyz[k];

                g_xx_0_xxyyzz_zz[k] = -g_xx_0_xxyzz_zz[k] * ab_y + g_xx_0_xxyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzzz_xx, g_xx_0_xxyzzz_xy, g_xx_0_xxyzzz_xz, g_xx_0_xxyzzz_yy, g_xx_0_xxyzzz_yz, g_xx_0_xxyzzz_zz, g_xx_0_xxzzz_xx, g_xx_0_xxzzz_xxy, g_xx_0_xxzzz_xy, g_xx_0_xxzzz_xyy, g_xx_0_xxzzz_xyz, g_xx_0_xxzzz_xz, g_xx_0_xxzzz_yy, g_xx_0_xxzzz_yyy, g_xx_0_xxzzz_yyz, g_xx_0_xxzzz_yz, g_xx_0_xxzzz_yzz, g_xx_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzzz_xx[k] = -g_xx_0_xxzzz_xx[k] * ab_y + g_xx_0_xxzzz_xxy[k];

                g_xx_0_xxyzzz_xy[k] = -g_xx_0_xxzzz_xy[k] * ab_y + g_xx_0_xxzzz_xyy[k];

                g_xx_0_xxyzzz_xz[k] = -g_xx_0_xxzzz_xz[k] * ab_y + g_xx_0_xxzzz_xyz[k];

                g_xx_0_xxyzzz_yy[k] = -g_xx_0_xxzzz_yy[k] * ab_y + g_xx_0_xxzzz_yyy[k];

                g_xx_0_xxyzzz_yz[k] = -g_xx_0_xxzzz_yz[k] * ab_y + g_xx_0_xxzzz_yyz[k];

                g_xx_0_xxyzzz_zz[k] = -g_xx_0_xxzzz_zz[k] * ab_y + g_xx_0_xxzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzzz_xx, g_xx_0_xxzzz_xxz, g_xx_0_xxzzz_xy, g_xx_0_xxzzz_xyz, g_xx_0_xxzzz_xz, g_xx_0_xxzzz_xzz, g_xx_0_xxzzz_yy, g_xx_0_xxzzz_yyz, g_xx_0_xxzzz_yz, g_xx_0_xxzzz_yzz, g_xx_0_xxzzz_zz, g_xx_0_xxzzz_zzz, g_xx_0_xxzzzz_xx, g_xx_0_xxzzzz_xy, g_xx_0_xxzzzz_xz, g_xx_0_xxzzzz_yy, g_xx_0_xxzzzz_yz, g_xx_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzzz_xx[k] = -g_xx_0_xxzzz_xx[k] * ab_z + g_xx_0_xxzzz_xxz[k];

                g_xx_0_xxzzzz_xy[k] = -g_xx_0_xxzzz_xy[k] * ab_z + g_xx_0_xxzzz_xyz[k];

                g_xx_0_xxzzzz_xz[k] = -g_xx_0_xxzzz_xz[k] * ab_z + g_xx_0_xxzzz_xzz[k];

                g_xx_0_xxzzzz_yy[k] = -g_xx_0_xxzzz_yy[k] * ab_z + g_xx_0_xxzzz_yyz[k];

                g_xx_0_xxzzzz_yz[k] = -g_xx_0_xxzzz_yz[k] * ab_z + g_xx_0_xxzzz_yzz[k];

                g_xx_0_xxzzzz_zz[k] = -g_xx_0_xxzzz_zz[k] * ab_z + g_xx_0_xxzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyy_xx, g_xx_0_xyyyy_xxy, g_xx_0_xyyyy_xy, g_xx_0_xyyyy_xyy, g_xx_0_xyyyy_xyz, g_xx_0_xyyyy_xz, g_xx_0_xyyyy_yy, g_xx_0_xyyyy_yyy, g_xx_0_xyyyy_yyz, g_xx_0_xyyyy_yz, g_xx_0_xyyyy_yzz, g_xx_0_xyyyy_zz, g_xx_0_xyyyyy_xx, g_xx_0_xyyyyy_xy, g_xx_0_xyyyyy_xz, g_xx_0_xyyyyy_yy, g_xx_0_xyyyyy_yz, g_xx_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyy_xx[k] = -g_xx_0_xyyyy_xx[k] * ab_y + g_xx_0_xyyyy_xxy[k];

                g_xx_0_xyyyyy_xy[k] = -g_xx_0_xyyyy_xy[k] * ab_y + g_xx_0_xyyyy_xyy[k];

                g_xx_0_xyyyyy_xz[k] = -g_xx_0_xyyyy_xz[k] * ab_y + g_xx_0_xyyyy_xyz[k];

                g_xx_0_xyyyyy_yy[k] = -g_xx_0_xyyyy_yy[k] * ab_y + g_xx_0_xyyyy_yyy[k];

                g_xx_0_xyyyyy_yz[k] = -g_xx_0_xyyyy_yz[k] * ab_y + g_xx_0_xyyyy_yyz[k];

                g_xx_0_xyyyyy_zz[k] = -g_xx_0_xyyyy_zz[k] * ab_y + g_xx_0_xyyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyyz_xx, g_xx_0_xyyyyz_xy, g_xx_0_xyyyyz_xz, g_xx_0_xyyyyz_yy, g_xx_0_xyyyyz_yz, g_xx_0_xyyyyz_zz, g_xx_0_xyyyz_xx, g_xx_0_xyyyz_xxy, g_xx_0_xyyyz_xy, g_xx_0_xyyyz_xyy, g_xx_0_xyyyz_xyz, g_xx_0_xyyyz_xz, g_xx_0_xyyyz_yy, g_xx_0_xyyyz_yyy, g_xx_0_xyyyz_yyz, g_xx_0_xyyyz_yz, g_xx_0_xyyyz_yzz, g_xx_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyz_xx[k] = -g_xx_0_xyyyz_xx[k] * ab_y + g_xx_0_xyyyz_xxy[k];

                g_xx_0_xyyyyz_xy[k] = -g_xx_0_xyyyz_xy[k] * ab_y + g_xx_0_xyyyz_xyy[k];

                g_xx_0_xyyyyz_xz[k] = -g_xx_0_xyyyz_xz[k] * ab_y + g_xx_0_xyyyz_xyz[k];

                g_xx_0_xyyyyz_yy[k] = -g_xx_0_xyyyz_yy[k] * ab_y + g_xx_0_xyyyz_yyy[k];

                g_xx_0_xyyyyz_yz[k] = -g_xx_0_xyyyz_yz[k] * ab_y + g_xx_0_xyyyz_yyz[k];

                g_xx_0_xyyyyz_zz[k] = -g_xx_0_xyyyz_zz[k] * ab_y + g_xx_0_xyyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyzz_xx, g_xx_0_xyyyzz_xy, g_xx_0_xyyyzz_xz, g_xx_0_xyyyzz_yy, g_xx_0_xyyyzz_yz, g_xx_0_xyyyzz_zz, g_xx_0_xyyzz_xx, g_xx_0_xyyzz_xxy, g_xx_0_xyyzz_xy, g_xx_0_xyyzz_xyy, g_xx_0_xyyzz_xyz, g_xx_0_xyyzz_xz, g_xx_0_xyyzz_yy, g_xx_0_xyyzz_yyy, g_xx_0_xyyzz_yyz, g_xx_0_xyyzz_yz, g_xx_0_xyyzz_yzz, g_xx_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyzz_xx[k] = -g_xx_0_xyyzz_xx[k] * ab_y + g_xx_0_xyyzz_xxy[k];

                g_xx_0_xyyyzz_xy[k] = -g_xx_0_xyyzz_xy[k] * ab_y + g_xx_0_xyyzz_xyy[k];

                g_xx_0_xyyyzz_xz[k] = -g_xx_0_xyyzz_xz[k] * ab_y + g_xx_0_xyyzz_xyz[k];

                g_xx_0_xyyyzz_yy[k] = -g_xx_0_xyyzz_yy[k] * ab_y + g_xx_0_xyyzz_yyy[k];

                g_xx_0_xyyyzz_yz[k] = -g_xx_0_xyyzz_yz[k] * ab_y + g_xx_0_xyyzz_yyz[k];

                g_xx_0_xyyyzz_zz[k] = -g_xx_0_xyyzz_zz[k] * ab_y + g_xx_0_xyyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzzz_xx, g_xx_0_xyyzzz_xy, g_xx_0_xyyzzz_xz, g_xx_0_xyyzzz_yy, g_xx_0_xyyzzz_yz, g_xx_0_xyyzzz_zz, g_xx_0_xyzzz_xx, g_xx_0_xyzzz_xxy, g_xx_0_xyzzz_xy, g_xx_0_xyzzz_xyy, g_xx_0_xyzzz_xyz, g_xx_0_xyzzz_xz, g_xx_0_xyzzz_yy, g_xx_0_xyzzz_yyy, g_xx_0_xyzzz_yyz, g_xx_0_xyzzz_yz, g_xx_0_xyzzz_yzz, g_xx_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzzz_xx[k] = -g_xx_0_xyzzz_xx[k] * ab_y + g_xx_0_xyzzz_xxy[k];

                g_xx_0_xyyzzz_xy[k] = -g_xx_0_xyzzz_xy[k] * ab_y + g_xx_0_xyzzz_xyy[k];

                g_xx_0_xyyzzz_xz[k] = -g_xx_0_xyzzz_xz[k] * ab_y + g_xx_0_xyzzz_xyz[k];

                g_xx_0_xyyzzz_yy[k] = -g_xx_0_xyzzz_yy[k] * ab_y + g_xx_0_xyzzz_yyy[k];

                g_xx_0_xyyzzz_yz[k] = -g_xx_0_xyzzz_yz[k] * ab_y + g_xx_0_xyzzz_yyz[k];

                g_xx_0_xyyzzz_zz[k] = -g_xx_0_xyzzz_zz[k] * ab_y + g_xx_0_xyzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzzz_xx, g_xx_0_xyzzzz_xy, g_xx_0_xyzzzz_xz, g_xx_0_xyzzzz_yy, g_xx_0_xyzzzz_yz, g_xx_0_xyzzzz_zz, g_xx_0_xzzzz_xx, g_xx_0_xzzzz_xxy, g_xx_0_xzzzz_xy, g_xx_0_xzzzz_xyy, g_xx_0_xzzzz_xyz, g_xx_0_xzzzz_xz, g_xx_0_xzzzz_yy, g_xx_0_xzzzz_yyy, g_xx_0_xzzzz_yyz, g_xx_0_xzzzz_yz, g_xx_0_xzzzz_yzz, g_xx_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzzz_xx[k] = -g_xx_0_xzzzz_xx[k] * ab_y + g_xx_0_xzzzz_xxy[k];

                g_xx_0_xyzzzz_xy[k] = -g_xx_0_xzzzz_xy[k] * ab_y + g_xx_0_xzzzz_xyy[k];

                g_xx_0_xyzzzz_xz[k] = -g_xx_0_xzzzz_xz[k] * ab_y + g_xx_0_xzzzz_xyz[k];

                g_xx_0_xyzzzz_yy[k] = -g_xx_0_xzzzz_yy[k] * ab_y + g_xx_0_xzzzz_yyy[k];

                g_xx_0_xyzzzz_yz[k] = -g_xx_0_xzzzz_yz[k] * ab_y + g_xx_0_xzzzz_yyz[k];

                g_xx_0_xyzzzz_zz[k] = -g_xx_0_xzzzz_zz[k] * ab_y + g_xx_0_xzzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzzz_xx, g_xx_0_xzzzz_xxz, g_xx_0_xzzzz_xy, g_xx_0_xzzzz_xyz, g_xx_0_xzzzz_xz, g_xx_0_xzzzz_xzz, g_xx_0_xzzzz_yy, g_xx_0_xzzzz_yyz, g_xx_0_xzzzz_yz, g_xx_0_xzzzz_yzz, g_xx_0_xzzzz_zz, g_xx_0_xzzzz_zzz, g_xx_0_xzzzzz_xx, g_xx_0_xzzzzz_xy, g_xx_0_xzzzzz_xz, g_xx_0_xzzzzz_yy, g_xx_0_xzzzzz_yz, g_xx_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzzz_xx[k] = -g_xx_0_xzzzz_xx[k] * ab_z + g_xx_0_xzzzz_xxz[k];

                g_xx_0_xzzzzz_xy[k] = -g_xx_0_xzzzz_xy[k] * ab_z + g_xx_0_xzzzz_xyz[k];

                g_xx_0_xzzzzz_xz[k] = -g_xx_0_xzzzz_xz[k] * ab_z + g_xx_0_xzzzz_xzz[k];

                g_xx_0_xzzzzz_yy[k] = -g_xx_0_xzzzz_yy[k] * ab_z + g_xx_0_xzzzz_yyz[k];

                g_xx_0_xzzzzz_yz[k] = -g_xx_0_xzzzz_yz[k] * ab_z + g_xx_0_xzzzz_yzz[k];

                g_xx_0_xzzzzz_zz[k] = -g_xx_0_xzzzz_zz[k] * ab_z + g_xx_0_xzzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyy_xx, g_xx_0_yyyyy_xxy, g_xx_0_yyyyy_xy, g_xx_0_yyyyy_xyy, g_xx_0_yyyyy_xyz, g_xx_0_yyyyy_xz, g_xx_0_yyyyy_yy, g_xx_0_yyyyy_yyy, g_xx_0_yyyyy_yyz, g_xx_0_yyyyy_yz, g_xx_0_yyyyy_yzz, g_xx_0_yyyyy_zz, g_xx_0_yyyyyy_xx, g_xx_0_yyyyyy_xy, g_xx_0_yyyyyy_xz, g_xx_0_yyyyyy_yy, g_xx_0_yyyyyy_yz, g_xx_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyy_xx[k] = -g_xx_0_yyyyy_xx[k] * ab_y + g_xx_0_yyyyy_xxy[k];

                g_xx_0_yyyyyy_xy[k] = -g_xx_0_yyyyy_xy[k] * ab_y + g_xx_0_yyyyy_xyy[k];

                g_xx_0_yyyyyy_xz[k] = -g_xx_0_yyyyy_xz[k] * ab_y + g_xx_0_yyyyy_xyz[k];

                g_xx_0_yyyyyy_yy[k] = -g_xx_0_yyyyy_yy[k] * ab_y + g_xx_0_yyyyy_yyy[k];

                g_xx_0_yyyyyy_yz[k] = -g_xx_0_yyyyy_yz[k] * ab_y + g_xx_0_yyyyy_yyz[k];

                g_xx_0_yyyyyy_zz[k] = -g_xx_0_yyyyy_zz[k] * ab_y + g_xx_0_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyyz_xx, g_xx_0_yyyyyz_xy, g_xx_0_yyyyyz_xz, g_xx_0_yyyyyz_yy, g_xx_0_yyyyyz_yz, g_xx_0_yyyyyz_zz, g_xx_0_yyyyz_xx, g_xx_0_yyyyz_xxy, g_xx_0_yyyyz_xy, g_xx_0_yyyyz_xyy, g_xx_0_yyyyz_xyz, g_xx_0_yyyyz_xz, g_xx_0_yyyyz_yy, g_xx_0_yyyyz_yyy, g_xx_0_yyyyz_yyz, g_xx_0_yyyyz_yz, g_xx_0_yyyyz_yzz, g_xx_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyz_xx[k] = -g_xx_0_yyyyz_xx[k] * ab_y + g_xx_0_yyyyz_xxy[k];

                g_xx_0_yyyyyz_xy[k] = -g_xx_0_yyyyz_xy[k] * ab_y + g_xx_0_yyyyz_xyy[k];

                g_xx_0_yyyyyz_xz[k] = -g_xx_0_yyyyz_xz[k] * ab_y + g_xx_0_yyyyz_xyz[k];

                g_xx_0_yyyyyz_yy[k] = -g_xx_0_yyyyz_yy[k] * ab_y + g_xx_0_yyyyz_yyy[k];

                g_xx_0_yyyyyz_yz[k] = -g_xx_0_yyyyz_yz[k] * ab_y + g_xx_0_yyyyz_yyz[k];

                g_xx_0_yyyyyz_zz[k] = -g_xx_0_yyyyz_zz[k] * ab_y + g_xx_0_yyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyzz_xx, g_xx_0_yyyyzz_xy, g_xx_0_yyyyzz_xz, g_xx_0_yyyyzz_yy, g_xx_0_yyyyzz_yz, g_xx_0_yyyyzz_zz, g_xx_0_yyyzz_xx, g_xx_0_yyyzz_xxy, g_xx_0_yyyzz_xy, g_xx_0_yyyzz_xyy, g_xx_0_yyyzz_xyz, g_xx_0_yyyzz_xz, g_xx_0_yyyzz_yy, g_xx_0_yyyzz_yyy, g_xx_0_yyyzz_yyz, g_xx_0_yyyzz_yz, g_xx_0_yyyzz_yzz, g_xx_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyzz_xx[k] = -g_xx_0_yyyzz_xx[k] * ab_y + g_xx_0_yyyzz_xxy[k];

                g_xx_0_yyyyzz_xy[k] = -g_xx_0_yyyzz_xy[k] * ab_y + g_xx_0_yyyzz_xyy[k];

                g_xx_0_yyyyzz_xz[k] = -g_xx_0_yyyzz_xz[k] * ab_y + g_xx_0_yyyzz_xyz[k];

                g_xx_0_yyyyzz_yy[k] = -g_xx_0_yyyzz_yy[k] * ab_y + g_xx_0_yyyzz_yyy[k];

                g_xx_0_yyyyzz_yz[k] = -g_xx_0_yyyzz_yz[k] * ab_y + g_xx_0_yyyzz_yyz[k];

                g_xx_0_yyyyzz_zz[k] = -g_xx_0_yyyzz_zz[k] * ab_y + g_xx_0_yyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzzz_xx, g_xx_0_yyyzzz_xy, g_xx_0_yyyzzz_xz, g_xx_0_yyyzzz_yy, g_xx_0_yyyzzz_yz, g_xx_0_yyyzzz_zz, g_xx_0_yyzzz_xx, g_xx_0_yyzzz_xxy, g_xx_0_yyzzz_xy, g_xx_0_yyzzz_xyy, g_xx_0_yyzzz_xyz, g_xx_0_yyzzz_xz, g_xx_0_yyzzz_yy, g_xx_0_yyzzz_yyy, g_xx_0_yyzzz_yyz, g_xx_0_yyzzz_yz, g_xx_0_yyzzz_yzz, g_xx_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzzz_xx[k] = -g_xx_0_yyzzz_xx[k] * ab_y + g_xx_0_yyzzz_xxy[k];

                g_xx_0_yyyzzz_xy[k] = -g_xx_0_yyzzz_xy[k] * ab_y + g_xx_0_yyzzz_xyy[k];

                g_xx_0_yyyzzz_xz[k] = -g_xx_0_yyzzz_xz[k] * ab_y + g_xx_0_yyzzz_xyz[k];

                g_xx_0_yyyzzz_yy[k] = -g_xx_0_yyzzz_yy[k] * ab_y + g_xx_0_yyzzz_yyy[k];

                g_xx_0_yyyzzz_yz[k] = -g_xx_0_yyzzz_yz[k] * ab_y + g_xx_0_yyzzz_yyz[k];

                g_xx_0_yyyzzz_zz[k] = -g_xx_0_yyzzz_zz[k] * ab_y + g_xx_0_yyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzzz_xx, g_xx_0_yyzzzz_xy, g_xx_0_yyzzzz_xz, g_xx_0_yyzzzz_yy, g_xx_0_yyzzzz_yz, g_xx_0_yyzzzz_zz, g_xx_0_yzzzz_xx, g_xx_0_yzzzz_xxy, g_xx_0_yzzzz_xy, g_xx_0_yzzzz_xyy, g_xx_0_yzzzz_xyz, g_xx_0_yzzzz_xz, g_xx_0_yzzzz_yy, g_xx_0_yzzzz_yyy, g_xx_0_yzzzz_yyz, g_xx_0_yzzzz_yz, g_xx_0_yzzzz_yzz, g_xx_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzzz_xx[k] = -g_xx_0_yzzzz_xx[k] * ab_y + g_xx_0_yzzzz_xxy[k];

                g_xx_0_yyzzzz_xy[k] = -g_xx_0_yzzzz_xy[k] * ab_y + g_xx_0_yzzzz_xyy[k];

                g_xx_0_yyzzzz_xz[k] = -g_xx_0_yzzzz_xz[k] * ab_y + g_xx_0_yzzzz_xyz[k];

                g_xx_0_yyzzzz_yy[k] = -g_xx_0_yzzzz_yy[k] * ab_y + g_xx_0_yzzzz_yyy[k];

                g_xx_0_yyzzzz_yz[k] = -g_xx_0_yzzzz_yz[k] * ab_y + g_xx_0_yzzzz_yyz[k];

                g_xx_0_yyzzzz_zz[k] = -g_xx_0_yzzzz_zz[k] * ab_y + g_xx_0_yzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzzz_xx, g_xx_0_yzzzzz_xy, g_xx_0_yzzzzz_xz, g_xx_0_yzzzzz_yy, g_xx_0_yzzzzz_yz, g_xx_0_yzzzzz_zz, g_xx_0_zzzzz_xx, g_xx_0_zzzzz_xxy, g_xx_0_zzzzz_xy, g_xx_0_zzzzz_xyy, g_xx_0_zzzzz_xyz, g_xx_0_zzzzz_xz, g_xx_0_zzzzz_yy, g_xx_0_zzzzz_yyy, g_xx_0_zzzzz_yyz, g_xx_0_zzzzz_yz, g_xx_0_zzzzz_yzz, g_xx_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzzz_xx[k] = -g_xx_0_zzzzz_xx[k] * ab_y + g_xx_0_zzzzz_xxy[k];

                g_xx_0_yzzzzz_xy[k] = -g_xx_0_zzzzz_xy[k] * ab_y + g_xx_0_zzzzz_xyy[k];

                g_xx_0_yzzzzz_xz[k] = -g_xx_0_zzzzz_xz[k] * ab_y + g_xx_0_zzzzz_xyz[k];

                g_xx_0_yzzzzz_yy[k] = -g_xx_0_zzzzz_yy[k] * ab_y + g_xx_0_zzzzz_yyy[k];

                g_xx_0_yzzzzz_yz[k] = -g_xx_0_zzzzz_yz[k] * ab_y + g_xx_0_zzzzz_yyz[k];

                g_xx_0_yzzzzz_zz[k] = -g_xx_0_zzzzz_zz[k] * ab_y + g_xx_0_zzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 164 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzzz_xx, g_xx_0_zzzzz_xxz, g_xx_0_zzzzz_xy, g_xx_0_zzzzz_xyz, g_xx_0_zzzzz_xz, g_xx_0_zzzzz_xzz, g_xx_0_zzzzz_yy, g_xx_0_zzzzz_yyz, g_xx_0_zzzzz_yz, g_xx_0_zzzzz_yzz, g_xx_0_zzzzz_zz, g_xx_0_zzzzz_zzz, g_xx_0_zzzzzz_xx, g_xx_0_zzzzzz_xy, g_xx_0_zzzzzz_xz, g_xx_0_zzzzzz_yy, g_xx_0_zzzzzz_yz, g_xx_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzzz_xx[k] = -g_xx_0_zzzzz_xx[k] * ab_z + g_xx_0_zzzzz_xxz[k];

                g_xx_0_zzzzzz_xy[k] = -g_xx_0_zzzzz_xy[k] * ab_z + g_xx_0_zzzzz_xyz[k];

                g_xx_0_zzzzzz_xz[k] = -g_xx_0_zzzzz_xz[k] * ab_z + g_xx_0_zzzzz_xzz[k];

                g_xx_0_zzzzzz_yy[k] = -g_xx_0_zzzzz_yy[k] * ab_z + g_xx_0_zzzzz_yyz[k];

                g_xx_0_zzzzzz_yz[k] = -g_xx_0_zzzzz_yz[k] * ab_z + g_xx_0_zzzzz_yzz[k];

                g_xx_0_zzzzzz_zz[k] = -g_xx_0_zzzzz_zz[k] * ab_z + g_xx_0_zzzzz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_xx, g_xy_0_xxxxx_xxx, g_xy_0_xxxxx_xxy, g_xy_0_xxxxx_xxz, g_xy_0_xxxxx_xy, g_xy_0_xxxxx_xyy, g_xy_0_xxxxx_xyz, g_xy_0_xxxxx_xz, g_xy_0_xxxxx_xzz, g_xy_0_xxxxx_yy, g_xy_0_xxxxx_yz, g_xy_0_xxxxx_zz, g_xy_0_xxxxxx_xx, g_xy_0_xxxxxx_xy, g_xy_0_xxxxxx_xz, g_xy_0_xxxxxx_yy, g_xy_0_xxxxxx_yz, g_xy_0_xxxxxx_zz, g_y_0_xxxxx_xx, g_y_0_xxxxx_xy, g_y_0_xxxxx_xz, g_y_0_xxxxx_yy, g_y_0_xxxxx_yz, g_y_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxx_xx[k] = -g_y_0_xxxxx_xx[k] - g_xy_0_xxxxx_xx[k] * ab_x + g_xy_0_xxxxx_xxx[k];

                g_xy_0_xxxxxx_xy[k] = -g_y_0_xxxxx_xy[k] - g_xy_0_xxxxx_xy[k] * ab_x + g_xy_0_xxxxx_xxy[k];

                g_xy_0_xxxxxx_xz[k] = -g_y_0_xxxxx_xz[k] - g_xy_0_xxxxx_xz[k] * ab_x + g_xy_0_xxxxx_xxz[k];

                g_xy_0_xxxxxx_yy[k] = -g_y_0_xxxxx_yy[k] - g_xy_0_xxxxx_yy[k] * ab_x + g_xy_0_xxxxx_xyy[k];

                g_xy_0_xxxxxx_yz[k] = -g_y_0_xxxxx_yz[k] - g_xy_0_xxxxx_yz[k] * ab_x + g_xy_0_xxxxx_xyz[k];

                g_xy_0_xxxxxx_zz[k] = -g_y_0_xxxxx_zz[k] - g_xy_0_xxxxx_zz[k] * ab_x + g_xy_0_xxxxx_xzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxxy_xx, g_xy_0_xxxxxy_xy, g_xy_0_xxxxxy_xz, g_xy_0_xxxxxy_yy, g_xy_0_xxxxxy_yz, g_xy_0_xxxxxy_zz, g_xy_0_xxxxy_xx, g_xy_0_xxxxy_xxx, g_xy_0_xxxxy_xxy, g_xy_0_xxxxy_xxz, g_xy_0_xxxxy_xy, g_xy_0_xxxxy_xyy, g_xy_0_xxxxy_xyz, g_xy_0_xxxxy_xz, g_xy_0_xxxxy_xzz, g_xy_0_xxxxy_yy, g_xy_0_xxxxy_yz, g_xy_0_xxxxy_zz, g_y_0_xxxxy_xx, g_y_0_xxxxy_xy, g_y_0_xxxxy_xz, g_y_0_xxxxy_yy, g_y_0_xxxxy_yz, g_y_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxy_xx[k] = -g_y_0_xxxxy_xx[k] - g_xy_0_xxxxy_xx[k] * ab_x + g_xy_0_xxxxy_xxx[k];

                g_xy_0_xxxxxy_xy[k] = -g_y_0_xxxxy_xy[k] - g_xy_0_xxxxy_xy[k] * ab_x + g_xy_0_xxxxy_xxy[k];

                g_xy_0_xxxxxy_xz[k] = -g_y_0_xxxxy_xz[k] - g_xy_0_xxxxy_xz[k] * ab_x + g_xy_0_xxxxy_xxz[k];

                g_xy_0_xxxxxy_yy[k] = -g_y_0_xxxxy_yy[k] - g_xy_0_xxxxy_yy[k] * ab_x + g_xy_0_xxxxy_xyy[k];

                g_xy_0_xxxxxy_yz[k] = -g_y_0_xxxxy_yz[k] - g_xy_0_xxxxy_yz[k] * ab_x + g_xy_0_xxxxy_xyz[k];

                g_xy_0_xxxxxy_zz[k] = -g_y_0_xxxxy_zz[k] - g_xy_0_xxxxy_zz[k] * ab_x + g_xy_0_xxxxy_xzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_xx, g_xy_0_xxxxx_xxz, g_xy_0_xxxxx_xy, g_xy_0_xxxxx_xyz, g_xy_0_xxxxx_xz, g_xy_0_xxxxx_xzz, g_xy_0_xxxxx_yy, g_xy_0_xxxxx_yyz, g_xy_0_xxxxx_yz, g_xy_0_xxxxx_yzz, g_xy_0_xxxxx_zz, g_xy_0_xxxxx_zzz, g_xy_0_xxxxxz_xx, g_xy_0_xxxxxz_xy, g_xy_0_xxxxxz_xz, g_xy_0_xxxxxz_yy, g_xy_0_xxxxxz_yz, g_xy_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxz_xx[k] = -g_xy_0_xxxxx_xx[k] * ab_z + g_xy_0_xxxxx_xxz[k];

                g_xy_0_xxxxxz_xy[k] = -g_xy_0_xxxxx_xy[k] * ab_z + g_xy_0_xxxxx_xyz[k];

                g_xy_0_xxxxxz_xz[k] = -g_xy_0_xxxxx_xz[k] * ab_z + g_xy_0_xxxxx_xzz[k];

                g_xy_0_xxxxxz_yy[k] = -g_xy_0_xxxxx_yy[k] * ab_z + g_xy_0_xxxxx_yyz[k];

                g_xy_0_xxxxxz_yz[k] = -g_xy_0_xxxxx_yz[k] * ab_z + g_xy_0_xxxxx_yzz[k];

                g_xy_0_xxxxxz_zz[k] = -g_xy_0_xxxxx_zz[k] * ab_z + g_xy_0_xxxxx_zzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxyy_xx, g_xy_0_xxxxyy_xy, g_xy_0_xxxxyy_xz, g_xy_0_xxxxyy_yy, g_xy_0_xxxxyy_yz, g_xy_0_xxxxyy_zz, g_xy_0_xxxyy_xx, g_xy_0_xxxyy_xxx, g_xy_0_xxxyy_xxy, g_xy_0_xxxyy_xxz, g_xy_0_xxxyy_xy, g_xy_0_xxxyy_xyy, g_xy_0_xxxyy_xyz, g_xy_0_xxxyy_xz, g_xy_0_xxxyy_xzz, g_xy_0_xxxyy_yy, g_xy_0_xxxyy_yz, g_xy_0_xxxyy_zz, g_y_0_xxxyy_xx, g_y_0_xxxyy_xy, g_y_0_xxxyy_xz, g_y_0_xxxyy_yy, g_y_0_xxxyy_yz, g_y_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyy_xx[k] = -g_y_0_xxxyy_xx[k] - g_xy_0_xxxyy_xx[k] * ab_x + g_xy_0_xxxyy_xxx[k];

                g_xy_0_xxxxyy_xy[k] = -g_y_0_xxxyy_xy[k] - g_xy_0_xxxyy_xy[k] * ab_x + g_xy_0_xxxyy_xxy[k];

                g_xy_0_xxxxyy_xz[k] = -g_y_0_xxxyy_xz[k] - g_xy_0_xxxyy_xz[k] * ab_x + g_xy_0_xxxyy_xxz[k];

                g_xy_0_xxxxyy_yy[k] = -g_y_0_xxxyy_yy[k] - g_xy_0_xxxyy_yy[k] * ab_x + g_xy_0_xxxyy_xyy[k];

                g_xy_0_xxxxyy_yz[k] = -g_y_0_xxxyy_yz[k] - g_xy_0_xxxyy_yz[k] * ab_x + g_xy_0_xxxyy_xyz[k];

                g_xy_0_xxxxyy_zz[k] = -g_y_0_xxxyy_zz[k] - g_xy_0_xxxyy_zz[k] * ab_x + g_xy_0_xxxyy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_xx, g_xy_0_xxxxy_xxz, g_xy_0_xxxxy_xy, g_xy_0_xxxxy_xyz, g_xy_0_xxxxy_xz, g_xy_0_xxxxy_xzz, g_xy_0_xxxxy_yy, g_xy_0_xxxxy_yyz, g_xy_0_xxxxy_yz, g_xy_0_xxxxy_yzz, g_xy_0_xxxxy_zz, g_xy_0_xxxxy_zzz, g_xy_0_xxxxyz_xx, g_xy_0_xxxxyz_xy, g_xy_0_xxxxyz_xz, g_xy_0_xxxxyz_yy, g_xy_0_xxxxyz_yz, g_xy_0_xxxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyz_xx[k] = -g_xy_0_xxxxy_xx[k] * ab_z + g_xy_0_xxxxy_xxz[k];

                g_xy_0_xxxxyz_xy[k] = -g_xy_0_xxxxy_xy[k] * ab_z + g_xy_0_xxxxy_xyz[k];

                g_xy_0_xxxxyz_xz[k] = -g_xy_0_xxxxy_xz[k] * ab_z + g_xy_0_xxxxy_xzz[k];

                g_xy_0_xxxxyz_yy[k] = -g_xy_0_xxxxy_yy[k] * ab_z + g_xy_0_xxxxy_yyz[k];

                g_xy_0_xxxxyz_yz[k] = -g_xy_0_xxxxy_yz[k] * ab_z + g_xy_0_xxxxy_yzz[k];

                g_xy_0_xxxxyz_zz[k] = -g_xy_0_xxxxy_zz[k] * ab_z + g_xy_0_xxxxy_zzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxz_xx, g_xy_0_xxxxz_xxz, g_xy_0_xxxxz_xy, g_xy_0_xxxxz_xyz, g_xy_0_xxxxz_xz, g_xy_0_xxxxz_xzz, g_xy_0_xxxxz_yy, g_xy_0_xxxxz_yyz, g_xy_0_xxxxz_yz, g_xy_0_xxxxz_yzz, g_xy_0_xxxxz_zz, g_xy_0_xxxxz_zzz, g_xy_0_xxxxzz_xx, g_xy_0_xxxxzz_xy, g_xy_0_xxxxzz_xz, g_xy_0_xxxxzz_yy, g_xy_0_xxxxzz_yz, g_xy_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxzz_xx[k] = -g_xy_0_xxxxz_xx[k] * ab_z + g_xy_0_xxxxz_xxz[k];

                g_xy_0_xxxxzz_xy[k] = -g_xy_0_xxxxz_xy[k] * ab_z + g_xy_0_xxxxz_xyz[k];

                g_xy_0_xxxxzz_xz[k] = -g_xy_0_xxxxz_xz[k] * ab_z + g_xy_0_xxxxz_xzz[k];

                g_xy_0_xxxxzz_yy[k] = -g_xy_0_xxxxz_yy[k] * ab_z + g_xy_0_xxxxz_yyz[k];

                g_xy_0_xxxxzz_yz[k] = -g_xy_0_xxxxz_yz[k] * ab_z + g_xy_0_xxxxz_yzz[k];

                g_xy_0_xxxxzz_zz[k] = -g_xy_0_xxxxz_zz[k] * ab_z + g_xy_0_xxxxz_zzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyyy_xx, g_xy_0_xxxyyy_xy, g_xy_0_xxxyyy_xz, g_xy_0_xxxyyy_yy, g_xy_0_xxxyyy_yz, g_xy_0_xxxyyy_zz, g_xy_0_xxyyy_xx, g_xy_0_xxyyy_xxx, g_xy_0_xxyyy_xxy, g_xy_0_xxyyy_xxz, g_xy_0_xxyyy_xy, g_xy_0_xxyyy_xyy, g_xy_0_xxyyy_xyz, g_xy_0_xxyyy_xz, g_xy_0_xxyyy_xzz, g_xy_0_xxyyy_yy, g_xy_0_xxyyy_yz, g_xy_0_xxyyy_zz, g_y_0_xxyyy_xx, g_y_0_xxyyy_xy, g_y_0_xxyyy_xz, g_y_0_xxyyy_yy, g_y_0_xxyyy_yz, g_y_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyy_xx[k] = -g_y_0_xxyyy_xx[k] - g_xy_0_xxyyy_xx[k] * ab_x + g_xy_0_xxyyy_xxx[k];

                g_xy_0_xxxyyy_xy[k] = -g_y_0_xxyyy_xy[k] - g_xy_0_xxyyy_xy[k] * ab_x + g_xy_0_xxyyy_xxy[k];

                g_xy_0_xxxyyy_xz[k] = -g_y_0_xxyyy_xz[k] - g_xy_0_xxyyy_xz[k] * ab_x + g_xy_0_xxyyy_xxz[k];

                g_xy_0_xxxyyy_yy[k] = -g_y_0_xxyyy_yy[k] - g_xy_0_xxyyy_yy[k] * ab_x + g_xy_0_xxyyy_xyy[k];

                g_xy_0_xxxyyy_yz[k] = -g_y_0_xxyyy_yz[k] - g_xy_0_xxyyy_yz[k] * ab_x + g_xy_0_xxyyy_xyz[k];

                g_xy_0_xxxyyy_zz[k] = -g_y_0_xxyyy_zz[k] - g_xy_0_xxyyy_zz[k] * ab_x + g_xy_0_xxyyy_xzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_xx, g_xy_0_xxxyy_xxz, g_xy_0_xxxyy_xy, g_xy_0_xxxyy_xyz, g_xy_0_xxxyy_xz, g_xy_0_xxxyy_xzz, g_xy_0_xxxyy_yy, g_xy_0_xxxyy_yyz, g_xy_0_xxxyy_yz, g_xy_0_xxxyy_yzz, g_xy_0_xxxyy_zz, g_xy_0_xxxyy_zzz, g_xy_0_xxxyyz_xx, g_xy_0_xxxyyz_xy, g_xy_0_xxxyyz_xz, g_xy_0_xxxyyz_yy, g_xy_0_xxxyyz_yz, g_xy_0_xxxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyz_xx[k] = -g_xy_0_xxxyy_xx[k] * ab_z + g_xy_0_xxxyy_xxz[k];

                g_xy_0_xxxyyz_xy[k] = -g_xy_0_xxxyy_xy[k] * ab_z + g_xy_0_xxxyy_xyz[k];

                g_xy_0_xxxyyz_xz[k] = -g_xy_0_xxxyy_xz[k] * ab_z + g_xy_0_xxxyy_xzz[k];

                g_xy_0_xxxyyz_yy[k] = -g_xy_0_xxxyy_yy[k] * ab_z + g_xy_0_xxxyy_yyz[k];

                g_xy_0_xxxyyz_yz[k] = -g_xy_0_xxxyy_yz[k] * ab_z + g_xy_0_xxxyy_yzz[k];

                g_xy_0_xxxyyz_zz[k] = -g_xy_0_xxxyy_zz[k] * ab_z + g_xy_0_xxxyy_zzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyz_xx, g_xy_0_xxxyz_xxz, g_xy_0_xxxyz_xy, g_xy_0_xxxyz_xyz, g_xy_0_xxxyz_xz, g_xy_0_xxxyz_xzz, g_xy_0_xxxyz_yy, g_xy_0_xxxyz_yyz, g_xy_0_xxxyz_yz, g_xy_0_xxxyz_yzz, g_xy_0_xxxyz_zz, g_xy_0_xxxyz_zzz, g_xy_0_xxxyzz_xx, g_xy_0_xxxyzz_xy, g_xy_0_xxxyzz_xz, g_xy_0_xxxyzz_yy, g_xy_0_xxxyzz_yz, g_xy_0_xxxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyzz_xx[k] = -g_xy_0_xxxyz_xx[k] * ab_z + g_xy_0_xxxyz_xxz[k];

                g_xy_0_xxxyzz_xy[k] = -g_xy_0_xxxyz_xy[k] * ab_z + g_xy_0_xxxyz_xyz[k];

                g_xy_0_xxxyzz_xz[k] = -g_xy_0_xxxyz_xz[k] * ab_z + g_xy_0_xxxyz_xzz[k];

                g_xy_0_xxxyzz_yy[k] = -g_xy_0_xxxyz_yy[k] * ab_z + g_xy_0_xxxyz_yyz[k];

                g_xy_0_xxxyzz_yz[k] = -g_xy_0_xxxyz_yz[k] * ab_z + g_xy_0_xxxyz_yzz[k];

                g_xy_0_xxxyzz_zz[k] = -g_xy_0_xxxyz_zz[k] * ab_z + g_xy_0_xxxyz_zzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxzz_xx, g_xy_0_xxxzz_xxz, g_xy_0_xxxzz_xy, g_xy_0_xxxzz_xyz, g_xy_0_xxxzz_xz, g_xy_0_xxxzz_xzz, g_xy_0_xxxzz_yy, g_xy_0_xxxzz_yyz, g_xy_0_xxxzz_yz, g_xy_0_xxxzz_yzz, g_xy_0_xxxzz_zz, g_xy_0_xxxzz_zzz, g_xy_0_xxxzzz_xx, g_xy_0_xxxzzz_xy, g_xy_0_xxxzzz_xz, g_xy_0_xxxzzz_yy, g_xy_0_xxxzzz_yz, g_xy_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzzz_xx[k] = -g_xy_0_xxxzz_xx[k] * ab_z + g_xy_0_xxxzz_xxz[k];

                g_xy_0_xxxzzz_xy[k] = -g_xy_0_xxxzz_xy[k] * ab_z + g_xy_0_xxxzz_xyz[k];

                g_xy_0_xxxzzz_xz[k] = -g_xy_0_xxxzz_xz[k] * ab_z + g_xy_0_xxxzz_xzz[k];

                g_xy_0_xxxzzz_yy[k] = -g_xy_0_xxxzz_yy[k] * ab_z + g_xy_0_xxxzz_yyz[k];

                g_xy_0_xxxzzz_yz[k] = -g_xy_0_xxxzz_yz[k] * ab_z + g_xy_0_xxxzz_yzz[k];

                g_xy_0_xxxzzz_zz[k] = -g_xy_0_xxxzz_zz[k] * ab_z + g_xy_0_xxxzz_zzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyyy_xx, g_xy_0_xxyyyy_xy, g_xy_0_xxyyyy_xz, g_xy_0_xxyyyy_yy, g_xy_0_xxyyyy_yz, g_xy_0_xxyyyy_zz, g_xy_0_xyyyy_xx, g_xy_0_xyyyy_xxx, g_xy_0_xyyyy_xxy, g_xy_0_xyyyy_xxz, g_xy_0_xyyyy_xy, g_xy_0_xyyyy_xyy, g_xy_0_xyyyy_xyz, g_xy_0_xyyyy_xz, g_xy_0_xyyyy_xzz, g_xy_0_xyyyy_yy, g_xy_0_xyyyy_yz, g_xy_0_xyyyy_zz, g_y_0_xyyyy_xx, g_y_0_xyyyy_xy, g_y_0_xyyyy_xz, g_y_0_xyyyy_yy, g_y_0_xyyyy_yz, g_y_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyy_xx[k] = -g_y_0_xyyyy_xx[k] - g_xy_0_xyyyy_xx[k] * ab_x + g_xy_0_xyyyy_xxx[k];

                g_xy_0_xxyyyy_xy[k] = -g_y_0_xyyyy_xy[k] - g_xy_0_xyyyy_xy[k] * ab_x + g_xy_0_xyyyy_xxy[k];

                g_xy_0_xxyyyy_xz[k] = -g_y_0_xyyyy_xz[k] - g_xy_0_xyyyy_xz[k] * ab_x + g_xy_0_xyyyy_xxz[k];

                g_xy_0_xxyyyy_yy[k] = -g_y_0_xyyyy_yy[k] - g_xy_0_xyyyy_yy[k] * ab_x + g_xy_0_xyyyy_xyy[k];

                g_xy_0_xxyyyy_yz[k] = -g_y_0_xyyyy_yz[k] - g_xy_0_xyyyy_yz[k] * ab_x + g_xy_0_xyyyy_xyz[k];

                g_xy_0_xxyyyy_zz[k] = -g_y_0_xyyyy_zz[k] - g_xy_0_xyyyy_zz[k] * ab_x + g_xy_0_xyyyy_xzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_xx, g_xy_0_xxyyy_xxz, g_xy_0_xxyyy_xy, g_xy_0_xxyyy_xyz, g_xy_0_xxyyy_xz, g_xy_0_xxyyy_xzz, g_xy_0_xxyyy_yy, g_xy_0_xxyyy_yyz, g_xy_0_xxyyy_yz, g_xy_0_xxyyy_yzz, g_xy_0_xxyyy_zz, g_xy_0_xxyyy_zzz, g_xy_0_xxyyyz_xx, g_xy_0_xxyyyz_xy, g_xy_0_xxyyyz_xz, g_xy_0_xxyyyz_yy, g_xy_0_xxyyyz_yz, g_xy_0_xxyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyz_xx[k] = -g_xy_0_xxyyy_xx[k] * ab_z + g_xy_0_xxyyy_xxz[k];

                g_xy_0_xxyyyz_xy[k] = -g_xy_0_xxyyy_xy[k] * ab_z + g_xy_0_xxyyy_xyz[k];

                g_xy_0_xxyyyz_xz[k] = -g_xy_0_xxyyy_xz[k] * ab_z + g_xy_0_xxyyy_xzz[k];

                g_xy_0_xxyyyz_yy[k] = -g_xy_0_xxyyy_yy[k] * ab_z + g_xy_0_xxyyy_yyz[k];

                g_xy_0_xxyyyz_yz[k] = -g_xy_0_xxyyy_yz[k] * ab_z + g_xy_0_xxyyy_yzz[k];

                g_xy_0_xxyyyz_zz[k] = -g_xy_0_xxyyy_zz[k] * ab_z + g_xy_0_xxyyy_zzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyz_xx, g_xy_0_xxyyz_xxz, g_xy_0_xxyyz_xy, g_xy_0_xxyyz_xyz, g_xy_0_xxyyz_xz, g_xy_0_xxyyz_xzz, g_xy_0_xxyyz_yy, g_xy_0_xxyyz_yyz, g_xy_0_xxyyz_yz, g_xy_0_xxyyz_yzz, g_xy_0_xxyyz_zz, g_xy_0_xxyyz_zzz, g_xy_0_xxyyzz_xx, g_xy_0_xxyyzz_xy, g_xy_0_xxyyzz_xz, g_xy_0_xxyyzz_yy, g_xy_0_xxyyzz_yz, g_xy_0_xxyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyzz_xx[k] = -g_xy_0_xxyyz_xx[k] * ab_z + g_xy_0_xxyyz_xxz[k];

                g_xy_0_xxyyzz_xy[k] = -g_xy_0_xxyyz_xy[k] * ab_z + g_xy_0_xxyyz_xyz[k];

                g_xy_0_xxyyzz_xz[k] = -g_xy_0_xxyyz_xz[k] * ab_z + g_xy_0_xxyyz_xzz[k];

                g_xy_0_xxyyzz_yy[k] = -g_xy_0_xxyyz_yy[k] * ab_z + g_xy_0_xxyyz_yyz[k];

                g_xy_0_xxyyzz_yz[k] = -g_xy_0_xxyyz_yz[k] * ab_z + g_xy_0_xxyyz_yzz[k];

                g_xy_0_xxyyzz_zz[k] = -g_xy_0_xxyyz_zz[k] * ab_z + g_xy_0_xxyyz_zzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyzz_xx, g_xy_0_xxyzz_xxz, g_xy_0_xxyzz_xy, g_xy_0_xxyzz_xyz, g_xy_0_xxyzz_xz, g_xy_0_xxyzz_xzz, g_xy_0_xxyzz_yy, g_xy_0_xxyzz_yyz, g_xy_0_xxyzz_yz, g_xy_0_xxyzz_yzz, g_xy_0_xxyzz_zz, g_xy_0_xxyzz_zzz, g_xy_0_xxyzzz_xx, g_xy_0_xxyzzz_xy, g_xy_0_xxyzzz_xz, g_xy_0_xxyzzz_yy, g_xy_0_xxyzzz_yz, g_xy_0_xxyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzzz_xx[k] = -g_xy_0_xxyzz_xx[k] * ab_z + g_xy_0_xxyzz_xxz[k];

                g_xy_0_xxyzzz_xy[k] = -g_xy_0_xxyzz_xy[k] * ab_z + g_xy_0_xxyzz_xyz[k];

                g_xy_0_xxyzzz_xz[k] = -g_xy_0_xxyzz_xz[k] * ab_z + g_xy_0_xxyzz_xzz[k];

                g_xy_0_xxyzzz_yy[k] = -g_xy_0_xxyzz_yy[k] * ab_z + g_xy_0_xxyzz_yyz[k];

                g_xy_0_xxyzzz_yz[k] = -g_xy_0_xxyzz_yz[k] * ab_z + g_xy_0_xxyzz_yzz[k];

                g_xy_0_xxyzzz_zz[k] = -g_xy_0_xxyzz_zz[k] * ab_z + g_xy_0_xxyzz_zzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzzz_xx, g_xy_0_xxzzz_xxz, g_xy_0_xxzzz_xy, g_xy_0_xxzzz_xyz, g_xy_0_xxzzz_xz, g_xy_0_xxzzz_xzz, g_xy_0_xxzzz_yy, g_xy_0_xxzzz_yyz, g_xy_0_xxzzz_yz, g_xy_0_xxzzz_yzz, g_xy_0_xxzzz_zz, g_xy_0_xxzzz_zzz, g_xy_0_xxzzzz_xx, g_xy_0_xxzzzz_xy, g_xy_0_xxzzzz_xz, g_xy_0_xxzzzz_yy, g_xy_0_xxzzzz_yz, g_xy_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzzz_xx[k] = -g_xy_0_xxzzz_xx[k] * ab_z + g_xy_0_xxzzz_xxz[k];

                g_xy_0_xxzzzz_xy[k] = -g_xy_0_xxzzz_xy[k] * ab_z + g_xy_0_xxzzz_xyz[k];

                g_xy_0_xxzzzz_xz[k] = -g_xy_0_xxzzz_xz[k] * ab_z + g_xy_0_xxzzz_xzz[k];

                g_xy_0_xxzzzz_yy[k] = -g_xy_0_xxzzz_yy[k] * ab_z + g_xy_0_xxzzz_yyz[k];

                g_xy_0_xxzzzz_yz[k] = -g_xy_0_xxzzz_yz[k] * ab_z + g_xy_0_xxzzz_yzz[k];

                g_xy_0_xxzzzz_zz[k] = -g_xy_0_xxzzz_zz[k] * ab_z + g_xy_0_xxzzz_zzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyyy_xx, g_xy_0_xyyyyy_xy, g_xy_0_xyyyyy_xz, g_xy_0_xyyyyy_yy, g_xy_0_xyyyyy_yz, g_xy_0_xyyyyy_zz, g_xy_0_yyyyy_xx, g_xy_0_yyyyy_xxx, g_xy_0_yyyyy_xxy, g_xy_0_yyyyy_xxz, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_xyy, g_xy_0_yyyyy_xyz, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_xzz, g_xy_0_yyyyy_yy, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_zz, g_y_0_yyyyy_xx, g_y_0_yyyyy_xy, g_y_0_yyyyy_xz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yz, g_y_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyy_xx[k] = -g_y_0_yyyyy_xx[k] - g_xy_0_yyyyy_xx[k] * ab_x + g_xy_0_yyyyy_xxx[k];

                g_xy_0_xyyyyy_xy[k] = -g_y_0_yyyyy_xy[k] - g_xy_0_yyyyy_xy[k] * ab_x + g_xy_0_yyyyy_xxy[k];

                g_xy_0_xyyyyy_xz[k] = -g_y_0_yyyyy_xz[k] - g_xy_0_yyyyy_xz[k] * ab_x + g_xy_0_yyyyy_xxz[k];

                g_xy_0_xyyyyy_yy[k] = -g_y_0_yyyyy_yy[k] - g_xy_0_yyyyy_yy[k] * ab_x + g_xy_0_yyyyy_xyy[k];

                g_xy_0_xyyyyy_yz[k] = -g_y_0_yyyyy_yz[k] - g_xy_0_yyyyy_yz[k] * ab_x + g_xy_0_yyyyy_xyz[k];

                g_xy_0_xyyyyy_zz[k] = -g_y_0_yyyyy_zz[k] - g_xy_0_yyyyy_zz[k] * ab_x + g_xy_0_yyyyy_xzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_xx, g_xy_0_xyyyy_xxz, g_xy_0_xyyyy_xy, g_xy_0_xyyyy_xyz, g_xy_0_xyyyy_xz, g_xy_0_xyyyy_xzz, g_xy_0_xyyyy_yy, g_xy_0_xyyyy_yyz, g_xy_0_xyyyy_yz, g_xy_0_xyyyy_yzz, g_xy_0_xyyyy_zz, g_xy_0_xyyyy_zzz, g_xy_0_xyyyyz_xx, g_xy_0_xyyyyz_xy, g_xy_0_xyyyyz_xz, g_xy_0_xyyyyz_yy, g_xy_0_xyyyyz_yz, g_xy_0_xyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyz_xx[k] = -g_xy_0_xyyyy_xx[k] * ab_z + g_xy_0_xyyyy_xxz[k];

                g_xy_0_xyyyyz_xy[k] = -g_xy_0_xyyyy_xy[k] * ab_z + g_xy_0_xyyyy_xyz[k];

                g_xy_0_xyyyyz_xz[k] = -g_xy_0_xyyyy_xz[k] * ab_z + g_xy_0_xyyyy_xzz[k];

                g_xy_0_xyyyyz_yy[k] = -g_xy_0_xyyyy_yy[k] * ab_z + g_xy_0_xyyyy_yyz[k];

                g_xy_0_xyyyyz_yz[k] = -g_xy_0_xyyyy_yz[k] * ab_z + g_xy_0_xyyyy_yzz[k];

                g_xy_0_xyyyyz_zz[k] = -g_xy_0_xyyyy_zz[k] * ab_z + g_xy_0_xyyyy_zzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyz_xx, g_xy_0_xyyyz_xxz, g_xy_0_xyyyz_xy, g_xy_0_xyyyz_xyz, g_xy_0_xyyyz_xz, g_xy_0_xyyyz_xzz, g_xy_0_xyyyz_yy, g_xy_0_xyyyz_yyz, g_xy_0_xyyyz_yz, g_xy_0_xyyyz_yzz, g_xy_0_xyyyz_zz, g_xy_0_xyyyz_zzz, g_xy_0_xyyyzz_xx, g_xy_0_xyyyzz_xy, g_xy_0_xyyyzz_xz, g_xy_0_xyyyzz_yy, g_xy_0_xyyyzz_yz, g_xy_0_xyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyzz_xx[k] = -g_xy_0_xyyyz_xx[k] * ab_z + g_xy_0_xyyyz_xxz[k];

                g_xy_0_xyyyzz_xy[k] = -g_xy_0_xyyyz_xy[k] * ab_z + g_xy_0_xyyyz_xyz[k];

                g_xy_0_xyyyzz_xz[k] = -g_xy_0_xyyyz_xz[k] * ab_z + g_xy_0_xyyyz_xzz[k];

                g_xy_0_xyyyzz_yy[k] = -g_xy_0_xyyyz_yy[k] * ab_z + g_xy_0_xyyyz_yyz[k];

                g_xy_0_xyyyzz_yz[k] = -g_xy_0_xyyyz_yz[k] * ab_z + g_xy_0_xyyyz_yzz[k];

                g_xy_0_xyyyzz_zz[k] = -g_xy_0_xyyyz_zz[k] * ab_z + g_xy_0_xyyyz_zzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyzz_xx, g_xy_0_xyyzz_xxz, g_xy_0_xyyzz_xy, g_xy_0_xyyzz_xyz, g_xy_0_xyyzz_xz, g_xy_0_xyyzz_xzz, g_xy_0_xyyzz_yy, g_xy_0_xyyzz_yyz, g_xy_0_xyyzz_yz, g_xy_0_xyyzz_yzz, g_xy_0_xyyzz_zz, g_xy_0_xyyzz_zzz, g_xy_0_xyyzzz_xx, g_xy_0_xyyzzz_xy, g_xy_0_xyyzzz_xz, g_xy_0_xyyzzz_yy, g_xy_0_xyyzzz_yz, g_xy_0_xyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzzz_xx[k] = -g_xy_0_xyyzz_xx[k] * ab_z + g_xy_0_xyyzz_xxz[k];

                g_xy_0_xyyzzz_xy[k] = -g_xy_0_xyyzz_xy[k] * ab_z + g_xy_0_xyyzz_xyz[k];

                g_xy_0_xyyzzz_xz[k] = -g_xy_0_xyyzz_xz[k] * ab_z + g_xy_0_xyyzz_xzz[k];

                g_xy_0_xyyzzz_yy[k] = -g_xy_0_xyyzz_yy[k] * ab_z + g_xy_0_xyyzz_yyz[k];

                g_xy_0_xyyzzz_yz[k] = -g_xy_0_xyyzz_yz[k] * ab_z + g_xy_0_xyyzz_yzz[k];

                g_xy_0_xyyzzz_zz[k] = -g_xy_0_xyyzz_zz[k] * ab_z + g_xy_0_xyyzz_zzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzzz_xx, g_xy_0_xyzzz_xxz, g_xy_0_xyzzz_xy, g_xy_0_xyzzz_xyz, g_xy_0_xyzzz_xz, g_xy_0_xyzzz_xzz, g_xy_0_xyzzz_yy, g_xy_0_xyzzz_yyz, g_xy_0_xyzzz_yz, g_xy_0_xyzzz_yzz, g_xy_0_xyzzz_zz, g_xy_0_xyzzz_zzz, g_xy_0_xyzzzz_xx, g_xy_0_xyzzzz_xy, g_xy_0_xyzzzz_xz, g_xy_0_xyzzzz_yy, g_xy_0_xyzzzz_yz, g_xy_0_xyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzzz_xx[k] = -g_xy_0_xyzzz_xx[k] * ab_z + g_xy_0_xyzzz_xxz[k];

                g_xy_0_xyzzzz_xy[k] = -g_xy_0_xyzzz_xy[k] * ab_z + g_xy_0_xyzzz_xyz[k];

                g_xy_0_xyzzzz_xz[k] = -g_xy_0_xyzzz_xz[k] * ab_z + g_xy_0_xyzzz_xzz[k];

                g_xy_0_xyzzzz_yy[k] = -g_xy_0_xyzzz_yy[k] * ab_z + g_xy_0_xyzzz_yyz[k];

                g_xy_0_xyzzzz_yz[k] = -g_xy_0_xyzzz_yz[k] * ab_z + g_xy_0_xyzzz_yzz[k];

                g_xy_0_xyzzzz_zz[k] = -g_xy_0_xyzzz_zz[k] * ab_z + g_xy_0_xyzzz_zzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzzz_xx, g_xy_0_xzzzz_xxz, g_xy_0_xzzzz_xy, g_xy_0_xzzzz_xyz, g_xy_0_xzzzz_xz, g_xy_0_xzzzz_xzz, g_xy_0_xzzzz_yy, g_xy_0_xzzzz_yyz, g_xy_0_xzzzz_yz, g_xy_0_xzzzz_yzz, g_xy_0_xzzzz_zz, g_xy_0_xzzzz_zzz, g_xy_0_xzzzzz_xx, g_xy_0_xzzzzz_xy, g_xy_0_xzzzzz_xz, g_xy_0_xzzzzz_yy, g_xy_0_xzzzzz_yz, g_xy_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzzz_xx[k] = -g_xy_0_xzzzz_xx[k] * ab_z + g_xy_0_xzzzz_xxz[k];

                g_xy_0_xzzzzz_xy[k] = -g_xy_0_xzzzz_xy[k] * ab_z + g_xy_0_xzzzz_xyz[k];

                g_xy_0_xzzzzz_xz[k] = -g_xy_0_xzzzz_xz[k] * ab_z + g_xy_0_xzzzz_xzz[k];

                g_xy_0_xzzzzz_yy[k] = -g_xy_0_xzzzz_yy[k] * ab_z + g_xy_0_xzzzz_yyz[k];

                g_xy_0_xzzzzz_yz[k] = -g_xy_0_xzzzz_yz[k] * ab_z + g_xy_0_xzzzz_yzz[k];

                g_xy_0_xzzzzz_zz[k] = -g_xy_0_xzzzz_zz[k] * ab_z + g_xy_0_xzzzz_zzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_xx, g_x_0_yyyyy_xy, g_x_0_yyyyy_xz, g_x_0_yyyyy_yy, g_x_0_yyyyy_yz, g_x_0_yyyyy_zz, g_xy_0_yyyyy_xx, g_xy_0_yyyyy_xxy, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_xyy, g_xy_0_yyyyy_xyz, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_yy, g_xy_0_yyyyy_yyy, g_xy_0_yyyyy_yyz, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_yzz, g_xy_0_yyyyy_zz, g_xy_0_yyyyyy_xx, g_xy_0_yyyyyy_xy, g_xy_0_yyyyyy_xz, g_xy_0_yyyyyy_yy, g_xy_0_yyyyyy_yz, g_xy_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyy_xx[k] = -g_x_0_yyyyy_xx[k] - g_xy_0_yyyyy_xx[k] * ab_y + g_xy_0_yyyyy_xxy[k];

                g_xy_0_yyyyyy_xy[k] = -g_x_0_yyyyy_xy[k] - g_xy_0_yyyyy_xy[k] * ab_y + g_xy_0_yyyyy_xyy[k];

                g_xy_0_yyyyyy_xz[k] = -g_x_0_yyyyy_xz[k] - g_xy_0_yyyyy_xz[k] * ab_y + g_xy_0_yyyyy_xyz[k];

                g_xy_0_yyyyyy_yy[k] = -g_x_0_yyyyy_yy[k] - g_xy_0_yyyyy_yy[k] * ab_y + g_xy_0_yyyyy_yyy[k];

                g_xy_0_yyyyyy_yz[k] = -g_x_0_yyyyy_yz[k] - g_xy_0_yyyyy_yz[k] * ab_y + g_xy_0_yyyyy_yyz[k];

                g_xy_0_yyyyyy_zz[k] = -g_x_0_yyyyy_zz[k] - g_xy_0_yyyyy_zz[k] * ab_y + g_xy_0_yyyyy_yzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 300 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 301 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 302 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 303 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 304 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyy_xx, g_xy_0_yyyyy_xxz, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_xyz, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_xzz, g_xy_0_yyyyy_yy, g_xy_0_yyyyy_yyz, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_yzz, g_xy_0_yyyyy_zz, g_xy_0_yyyyy_zzz, g_xy_0_yyyyyz_xx, g_xy_0_yyyyyz_xy, g_xy_0_yyyyyz_xz, g_xy_0_yyyyyz_yy, g_xy_0_yyyyyz_yz, g_xy_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyz_xx[k] = -g_xy_0_yyyyy_xx[k] * ab_z + g_xy_0_yyyyy_xxz[k];

                g_xy_0_yyyyyz_xy[k] = -g_xy_0_yyyyy_xy[k] * ab_z + g_xy_0_yyyyy_xyz[k];

                g_xy_0_yyyyyz_xz[k] = -g_xy_0_yyyyy_xz[k] * ab_z + g_xy_0_yyyyy_xzz[k];

                g_xy_0_yyyyyz_yy[k] = -g_xy_0_yyyyy_yy[k] * ab_z + g_xy_0_yyyyy_yyz[k];

                g_xy_0_yyyyyz_yz[k] = -g_xy_0_yyyyy_yz[k] * ab_z + g_xy_0_yyyyy_yzz[k];

                g_xy_0_yyyyyz_zz[k] = -g_xy_0_yyyyy_zz[k] * ab_z + g_xy_0_yyyyy_zzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 306 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 307 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 308 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 309 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 310 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyz_xx, g_xy_0_yyyyz_xxz, g_xy_0_yyyyz_xy, g_xy_0_yyyyz_xyz, g_xy_0_yyyyz_xz, g_xy_0_yyyyz_xzz, g_xy_0_yyyyz_yy, g_xy_0_yyyyz_yyz, g_xy_0_yyyyz_yz, g_xy_0_yyyyz_yzz, g_xy_0_yyyyz_zz, g_xy_0_yyyyz_zzz, g_xy_0_yyyyzz_xx, g_xy_0_yyyyzz_xy, g_xy_0_yyyyzz_xz, g_xy_0_yyyyzz_yy, g_xy_0_yyyyzz_yz, g_xy_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyzz_xx[k] = -g_xy_0_yyyyz_xx[k] * ab_z + g_xy_0_yyyyz_xxz[k];

                g_xy_0_yyyyzz_xy[k] = -g_xy_0_yyyyz_xy[k] * ab_z + g_xy_0_yyyyz_xyz[k];

                g_xy_0_yyyyzz_xz[k] = -g_xy_0_yyyyz_xz[k] * ab_z + g_xy_0_yyyyz_xzz[k];

                g_xy_0_yyyyzz_yy[k] = -g_xy_0_yyyyz_yy[k] * ab_z + g_xy_0_yyyyz_yyz[k];

                g_xy_0_yyyyzz_yz[k] = -g_xy_0_yyyyz_yz[k] * ab_z + g_xy_0_yyyyz_yzz[k];

                g_xy_0_yyyyzz_zz[k] = -g_xy_0_yyyyz_zz[k] * ab_z + g_xy_0_yyyyz_zzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 312 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 313 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 314 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyzz_xx, g_xy_0_yyyzz_xxz, g_xy_0_yyyzz_xy, g_xy_0_yyyzz_xyz, g_xy_0_yyyzz_xz, g_xy_0_yyyzz_xzz, g_xy_0_yyyzz_yy, g_xy_0_yyyzz_yyz, g_xy_0_yyyzz_yz, g_xy_0_yyyzz_yzz, g_xy_0_yyyzz_zz, g_xy_0_yyyzz_zzz, g_xy_0_yyyzzz_xx, g_xy_0_yyyzzz_xy, g_xy_0_yyyzzz_xz, g_xy_0_yyyzzz_yy, g_xy_0_yyyzzz_yz, g_xy_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzzz_xx[k] = -g_xy_0_yyyzz_xx[k] * ab_z + g_xy_0_yyyzz_xxz[k];

                g_xy_0_yyyzzz_xy[k] = -g_xy_0_yyyzz_xy[k] * ab_z + g_xy_0_yyyzz_xyz[k];

                g_xy_0_yyyzzz_xz[k] = -g_xy_0_yyyzz_xz[k] * ab_z + g_xy_0_yyyzz_xzz[k];

                g_xy_0_yyyzzz_yy[k] = -g_xy_0_yyyzz_yy[k] * ab_z + g_xy_0_yyyzz_yyz[k];

                g_xy_0_yyyzzz_yz[k] = -g_xy_0_yyyzz_yz[k] * ab_z + g_xy_0_yyyzz_yzz[k];

                g_xy_0_yyyzzz_zz[k] = -g_xy_0_yyyzz_zz[k] * ab_z + g_xy_0_yyyzz_zzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzzz_xx, g_xy_0_yyzzz_xxz, g_xy_0_yyzzz_xy, g_xy_0_yyzzz_xyz, g_xy_0_yyzzz_xz, g_xy_0_yyzzz_xzz, g_xy_0_yyzzz_yy, g_xy_0_yyzzz_yyz, g_xy_0_yyzzz_yz, g_xy_0_yyzzz_yzz, g_xy_0_yyzzz_zz, g_xy_0_yyzzz_zzz, g_xy_0_yyzzzz_xx, g_xy_0_yyzzzz_xy, g_xy_0_yyzzzz_xz, g_xy_0_yyzzzz_yy, g_xy_0_yyzzzz_yz, g_xy_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzzz_xx[k] = -g_xy_0_yyzzz_xx[k] * ab_z + g_xy_0_yyzzz_xxz[k];

                g_xy_0_yyzzzz_xy[k] = -g_xy_0_yyzzz_xy[k] * ab_z + g_xy_0_yyzzz_xyz[k];

                g_xy_0_yyzzzz_xz[k] = -g_xy_0_yyzzz_xz[k] * ab_z + g_xy_0_yyzzz_xzz[k];

                g_xy_0_yyzzzz_yy[k] = -g_xy_0_yyzzz_yy[k] * ab_z + g_xy_0_yyzzz_yyz[k];

                g_xy_0_yyzzzz_yz[k] = -g_xy_0_yyzzz_yz[k] * ab_z + g_xy_0_yyzzz_yzz[k];

                g_xy_0_yyzzzz_zz[k] = -g_xy_0_yyzzz_zz[k] * ab_z + g_xy_0_yyzzz_zzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzzz_xx, g_xy_0_yzzzz_xxz, g_xy_0_yzzzz_xy, g_xy_0_yzzzz_xyz, g_xy_0_yzzzz_xz, g_xy_0_yzzzz_xzz, g_xy_0_yzzzz_yy, g_xy_0_yzzzz_yyz, g_xy_0_yzzzz_yz, g_xy_0_yzzzz_yzz, g_xy_0_yzzzz_zz, g_xy_0_yzzzz_zzz, g_xy_0_yzzzzz_xx, g_xy_0_yzzzzz_xy, g_xy_0_yzzzzz_xz, g_xy_0_yzzzzz_yy, g_xy_0_yzzzzz_yz, g_xy_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzzz_xx[k] = -g_xy_0_yzzzz_xx[k] * ab_z + g_xy_0_yzzzz_xxz[k];

                g_xy_0_yzzzzz_xy[k] = -g_xy_0_yzzzz_xy[k] * ab_z + g_xy_0_yzzzz_xyz[k];

                g_xy_0_yzzzzz_xz[k] = -g_xy_0_yzzzz_xz[k] * ab_z + g_xy_0_yzzzz_xzz[k];

                g_xy_0_yzzzzz_yy[k] = -g_xy_0_yzzzz_yy[k] * ab_z + g_xy_0_yzzzz_yyz[k];

                g_xy_0_yzzzzz_yz[k] = -g_xy_0_yzzzz_yz[k] * ab_z + g_xy_0_yzzzz_yzz[k];

                g_xy_0_yzzzzz_zz[k] = -g_xy_0_yzzzz_zz[k] * ab_z + g_xy_0_yzzzz_zzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzzz_xx, g_xy_0_zzzzz_xxz, g_xy_0_zzzzz_xy, g_xy_0_zzzzz_xyz, g_xy_0_zzzzz_xz, g_xy_0_zzzzz_xzz, g_xy_0_zzzzz_yy, g_xy_0_zzzzz_yyz, g_xy_0_zzzzz_yz, g_xy_0_zzzzz_yzz, g_xy_0_zzzzz_zz, g_xy_0_zzzzz_zzz, g_xy_0_zzzzzz_xx, g_xy_0_zzzzzz_xy, g_xy_0_zzzzzz_xz, g_xy_0_zzzzzz_yy, g_xy_0_zzzzzz_yz, g_xy_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzzz_xx[k] = -g_xy_0_zzzzz_xx[k] * ab_z + g_xy_0_zzzzz_xxz[k];

                g_xy_0_zzzzzz_xy[k] = -g_xy_0_zzzzz_xy[k] * ab_z + g_xy_0_zzzzz_xyz[k];

                g_xy_0_zzzzzz_xz[k] = -g_xy_0_zzzzz_xz[k] * ab_z + g_xy_0_zzzzz_xzz[k];

                g_xy_0_zzzzzz_yy[k] = -g_xy_0_zzzzz_yy[k] * ab_z + g_xy_0_zzzzz_yyz[k];

                g_xy_0_zzzzzz_yz[k] = -g_xy_0_zzzzz_yz[k] * ab_z + g_xy_0_zzzzz_yzz[k];

                g_xy_0_zzzzzz_zz[k] = -g_xy_0_zzzzz_zz[k] * ab_z + g_xy_0_zzzzz_zzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_xx, g_xz_0_xxxxx_xxx, g_xz_0_xxxxx_xxy, g_xz_0_xxxxx_xxz, g_xz_0_xxxxx_xy, g_xz_0_xxxxx_xyy, g_xz_0_xxxxx_xyz, g_xz_0_xxxxx_xz, g_xz_0_xxxxx_xzz, g_xz_0_xxxxx_yy, g_xz_0_xxxxx_yz, g_xz_0_xxxxx_zz, g_xz_0_xxxxxx_xx, g_xz_0_xxxxxx_xy, g_xz_0_xxxxxx_xz, g_xz_0_xxxxxx_yy, g_xz_0_xxxxxx_yz, g_xz_0_xxxxxx_zz, g_z_0_xxxxx_xx, g_z_0_xxxxx_xy, g_z_0_xxxxx_xz, g_z_0_xxxxx_yy, g_z_0_xxxxx_yz, g_z_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxx_xx[k] = -g_z_0_xxxxx_xx[k] - g_xz_0_xxxxx_xx[k] * ab_x + g_xz_0_xxxxx_xxx[k];

                g_xz_0_xxxxxx_xy[k] = -g_z_0_xxxxx_xy[k] - g_xz_0_xxxxx_xy[k] * ab_x + g_xz_0_xxxxx_xxy[k];

                g_xz_0_xxxxxx_xz[k] = -g_z_0_xxxxx_xz[k] - g_xz_0_xxxxx_xz[k] * ab_x + g_xz_0_xxxxx_xxz[k];

                g_xz_0_xxxxxx_yy[k] = -g_z_0_xxxxx_yy[k] - g_xz_0_xxxxx_yy[k] * ab_x + g_xz_0_xxxxx_xyy[k];

                g_xz_0_xxxxxx_yz[k] = -g_z_0_xxxxx_yz[k] - g_xz_0_xxxxx_yz[k] * ab_x + g_xz_0_xxxxx_xyz[k];

                g_xz_0_xxxxxx_zz[k] = -g_z_0_xxxxx_zz[k] - g_xz_0_xxxxx_zz[k] * ab_x + g_xz_0_xxxxx_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_xx, g_xz_0_xxxxx_xxy, g_xz_0_xxxxx_xy, g_xz_0_xxxxx_xyy, g_xz_0_xxxxx_xyz, g_xz_0_xxxxx_xz, g_xz_0_xxxxx_yy, g_xz_0_xxxxx_yyy, g_xz_0_xxxxx_yyz, g_xz_0_xxxxx_yz, g_xz_0_xxxxx_yzz, g_xz_0_xxxxx_zz, g_xz_0_xxxxxy_xx, g_xz_0_xxxxxy_xy, g_xz_0_xxxxxy_xz, g_xz_0_xxxxxy_yy, g_xz_0_xxxxxy_yz, g_xz_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxy_xx[k] = -g_xz_0_xxxxx_xx[k] * ab_y + g_xz_0_xxxxx_xxy[k];

                g_xz_0_xxxxxy_xy[k] = -g_xz_0_xxxxx_xy[k] * ab_y + g_xz_0_xxxxx_xyy[k];

                g_xz_0_xxxxxy_xz[k] = -g_xz_0_xxxxx_xz[k] * ab_y + g_xz_0_xxxxx_xyz[k];

                g_xz_0_xxxxxy_yy[k] = -g_xz_0_xxxxx_yy[k] * ab_y + g_xz_0_xxxxx_yyy[k];

                g_xz_0_xxxxxy_yz[k] = -g_xz_0_xxxxx_yz[k] * ab_y + g_xz_0_xxxxx_yyz[k];

                g_xz_0_xxxxxy_zz[k] = -g_xz_0_xxxxx_zz[k] * ab_y + g_xz_0_xxxxx_yzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxxz_xx, g_xz_0_xxxxxz_xy, g_xz_0_xxxxxz_xz, g_xz_0_xxxxxz_yy, g_xz_0_xxxxxz_yz, g_xz_0_xxxxxz_zz, g_xz_0_xxxxz_xx, g_xz_0_xxxxz_xxx, g_xz_0_xxxxz_xxy, g_xz_0_xxxxz_xxz, g_xz_0_xxxxz_xy, g_xz_0_xxxxz_xyy, g_xz_0_xxxxz_xyz, g_xz_0_xxxxz_xz, g_xz_0_xxxxz_xzz, g_xz_0_xxxxz_yy, g_xz_0_xxxxz_yz, g_xz_0_xxxxz_zz, g_z_0_xxxxz_xx, g_z_0_xxxxz_xy, g_z_0_xxxxz_xz, g_z_0_xxxxz_yy, g_z_0_xxxxz_yz, g_z_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxz_xx[k] = -g_z_0_xxxxz_xx[k] - g_xz_0_xxxxz_xx[k] * ab_x + g_xz_0_xxxxz_xxx[k];

                g_xz_0_xxxxxz_xy[k] = -g_z_0_xxxxz_xy[k] - g_xz_0_xxxxz_xy[k] * ab_x + g_xz_0_xxxxz_xxy[k];

                g_xz_0_xxxxxz_xz[k] = -g_z_0_xxxxz_xz[k] - g_xz_0_xxxxz_xz[k] * ab_x + g_xz_0_xxxxz_xxz[k];

                g_xz_0_xxxxxz_yy[k] = -g_z_0_xxxxz_yy[k] - g_xz_0_xxxxz_yy[k] * ab_x + g_xz_0_xxxxz_xyy[k];

                g_xz_0_xxxxxz_yz[k] = -g_z_0_xxxxz_yz[k] - g_xz_0_xxxxz_yz[k] * ab_x + g_xz_0_xxxxz_xyz[k];

                g_xz_0_xxxxxz_zz[k] = -g_z_0_xxxxz_zz[k] - g_xz_0_xxxxz_zz[k] * ab_x + g_xz_0_xxxxz_xzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxy_xx, g_xz_0_xxxxy_xxy, g_xz_0_xxxxy_xy, g_xz_0_xxxxy_xyy, g_xz_0_xxxxy_xyz, g_xz_0_xxxxy_xz, g_xz_0_xxxxy_yy, g_xz_0_xxxxy_yyy, g_xz_0_xxxxy_yyz, g_xz_0_xxxxy_yz, g_xz_0_xxxxy_yzz, g_xz_0_xxxxy_zz, g_xz_0_xxxxyy_xx, g_xz_0_xxxxyy_xy, g_xz_0_xxxxyy_xz, g_xz_0_xxxxyy_yy, g_xz_0_xxxxyy_yz, g_xz_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyy_xx[k] = -g_xz_0_xxxxy_xx[k] * ab_y + g_xz_0_xxxxy_xxy[k];

                g_xz_0_xxxxyy_xy[k] = -g_xz_0_xxxxy_xy[k] * ab_y + g_xz_0_xxxxy_xyy[k];

                g_xz_0_xxxxyy_xz[k] = -g_xz_0_xxxxy_xz[k] * ab_y + g_xz_0_xxxxy_xyz[k];

                g_xz_0_xxxxyy_yy[k] = -g_xz_0_xxxxy_yy[k] * ab_y + g_xz_0_xxxxy_yyy[k];

                g_xz_0_xxxxyy_yz[k] = -g_xz_0_xxxxy_yz[k] * ab_y + g_xz_0_xxxxy_yyz[k];

                g_xz_0_xxxxyy_zz[k] = -g_xz_0_xxxxy_zz[k] * ab_y + g_xz_0_xxxxy_yzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxyz_xx, g_xz_0_xxxxyz_xy, g_xz_0_xxxxyz_xz, g_xz_0_xxxxyz_yy, g_xz_0_xxxxyz_yz, g_xz_0_xxxxyz_zz, g_xz_0_xxxxz_xx, g_xz_0_xxxxz_xxy, g_xz_0_xxxxz_xy, g_xz_0_xxxxz_xyy, g_xz_0_xxxxz_xyz, g_xz_0_xxxxz_xz, g_xz_0_xxxxz_yy, g_xz_0_xxxxz_yyy, g_xz_0_xxxxz_yyz, g_xz_0_xxxxz_yz, g_xz_0_xxxxz_yzz, g_xz_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyz_xx[k] = -g_xz_0_xxxxz_xx[k] * ab_y + g_xz_0_xxxxz_xxy[k];

                g_xz_0_xxxxyz_xy[k] = -g_xz_0_xxxxz_xy[k] * ab_y + g_xz_0_xxxxz_xyy[k];

                g_xz_0_xxxxyz_xz[k] = -g_xz_0_xxxxz_xz[k] * ab_y + g_xz_0_xxxxz_xyz[k];

                g_xz_0_xxxxyz_yy[k] = -g_xz_0_xxxxz_yy[k] * ab_y + g_xz_0_xxxxz_yyy[k];

                g_xz_0_xxxxyz_yz[k] = -g_xz_0_xxxxz_yz[k] * ab_y + g_xz_0_xxxxz_yyz[k];

                g_xz_0_xxxxyz_zz[k] = -g_xz_0_xxxxz_zz[k] * ab_y + g_xz_0_xxxxz_yzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxzz_xx, g_xz_0_xxxxzz_xy, g_xz_0_xxxxzz_xz, g_xz_0_xxxxzz_yy, g_xz_0_xxxxzz_yz, g_xz_0_xxxxzz_zz, g_xz_0_xxxzz_xx, g_xz_0_xxxzz_xxx, g_xz_0_xxxzz_xxy, g_xz_0_xxxzz_xxz, g_xz_0_xxxzz_xy, g_xz_0_xxxzz_xyy, g_xz_0_xxxzz_xyz, g_xz_0_xxxzz_xz, g_xz_0_xxxzz_xzz, g_xz_0_xxxzz_yy, g_xz_0_xxxzz_yz, g_xz_0_xxxzz_zz, g_z_0_xxxzz_xx, g_z_0_xxxzz_xy, g_z_0_xxxzz_xz, g_z_0_xxxzz_yy, g_z_0_xxxzz_yz, g_z_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxzz_xx[k] = -g_z_0_xxxzz_xx[k] - g_xz_0_xxxzz_xx[k] * ab_x + g_xz_0_xxxzz_xxx[k];

                g_xz_0_xxxxzz_xy[k] = -g_z_0_xxxzz_xy[k] - g_xz_0_xxxzz_xy[k] * ab_x + g_xz_0_xxxzz_xxy[k];

                g_xz_0_xxxxzz_xz[k] = -g_z_0_xxxzz_xz[k] - g_xz_0_xxxzz_xz[k] * ab_x + g_xz_0_xxxzz_xxz[k];

                g_xz_0_xxxxzz_yy[k] = -g_z_0_xxxzz_yy[k] - g_xz_0_xxxzz_yy[k] * ab_x + g_xz_0_xxxzz_xyy[k];

                g_xz_0_xxxxzz_yz[k] = -g_z_0_xxxzz_yz[k] - g_xz_0_xxxzz_yz[k] * ab_x + g_xz_0_xxxzz_xyz[k];

                g_xz_0_xxxxzz_zz[k] = -g_z_0_xxxzz_zz[k] - g_xz_0_xxxzz_zz[k] * ab_x + g_xz_0_xxxzz_xzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyy_xx, g_xz_0_xxxyy_xxy, g_xz_0_xxxyy_xy, g_xz_0_xxxyy_xyy, g_xz_0_xxxyy_xyz, g_xz_0_xxxyy_xz, g_xz_0_xxxyy_yy, g_xz_0_xxxyy_yyy, g_xz_0_xxxyy_yyz, g_xz_0_xxxyy_yz, g_xz_0_xxxyy_yzz, g_xz_0_xxxyy_zz, g_xz_0_xxxyyy_xx, g_xz_0_xxxyyy_xy, g_xz_0_xxxyyy_xz, g_xz_0_xxxyyy_yy, g_xz_0_xxxyyy_yz, g_xz_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyy_xx[k] = -g_xz_0_xxxyy_xx[k] * ab_y + g_xz_0_xxxyy_xxy[k];

                g_xz_0_xxxyyy_xy[k] = -g_xz_0_xxxyy_xy[k] * ab_y + g_xz_0_xxxyy_xyy[k];

                g_xz_0_xxxyyy_xz[k] = -g_xz_0_xxxyy_xz[k] * ab_y + g_xz_0_xxxyy_xyz[k];

                g_xz_0_xxxyyy_yy[k] = -g_xz_0_xxxyy_yy[k] * ab_y + g_xz_0_xxxyy_yyy[k];

                g_xz_0_xxxyyy_yz[k] = -g_xz_0_xxxyy_yz[k] * ab_y + g_xz_0_xxxyy_yyz[k];

                g_xz_0_xxxyyy_zz[k] = -g_xz_0_xxxyy_zz[k] * ab_y + g_xz_0_xxxyy_yzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 379 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyyz_xx, g_xz_0_xxxyyz_xy, g_xz_0_xxxyyz_xz, g_xz_0_xxxyyz_yy, g_xz_0_xxxyyz_yz, g_xz_0_xxxyyz_zz, g_xz_0_xxxyz_xx, g_xz_0_xxxyz_xxy, g_xz_0_xxxyz_xy, g_xz_0_xxxyz_xyy, g_xz_0_xxxyz_xyz, g_xz_0_xxxyz_xz, g_xz_0_xxxyz_yy, g_xz_0_xxxyz_yyy, g_xz_0_xxxyz_yyz, g_xz_0_xxxyz_yz, g_xz_0_xxxyz_yzz, g_xz_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyz_xx[k] = -g_xz_0_xxxyz_xx[k] * ab_y + g_xz_0_xxxyz_xxy[k];

                g_xz_0_xxxyyz_xy[k] = -g_xz_0_xxxyz_xy[k] * ab_y + g_xz_0_xxxyz_xyy[k];

                g_xz_0_xxxyyz_xz[k] = -g_xz_0_xxxyz_xz[k] * ab_y + g_xz_0_xxxyz_xyz[k];

                g_xz_0_xxxyyz_yy[k] = -g_xz_0_xxxyz_yy[k] * ab_y + g_xz_0_xxxyz_yyy[k];

                g_xz_0_xxxyyz_yz[k] = -g_xz_0_xxxyz_yz[k] * ab_y + g_xz_0_xxxyz_yyz[k];

                g_xz_0_xxxyyz_zz[k] = -g_xz_0_xxxyz_zz[k] * ab_y + g_xz_0_xxxyz_yzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyzz_xx, g_xz_0_xxxyzz_xy, g_xz_0_xxxyzz_xz, g_xz_0_xxxyzz_yy, g_xz_0_xxxyzz_yz, g_xz_0_xxxyzz_zz, g_xz_0_xxxzz_xx, g_xz_0_xxxzz_xxy, g_xz_0_xxxzz_xy, g_xz_0_xxxzz_xyy, g_xz_0_xxxzz_xyz, g_xz_0_xxxzz_xz, g_xz_0_xxxzz_yy, g_xz_0_xxxzz_yyy, g_xz_0_xxxzz_yyz, g_xz_0_xxxzz_yz, g_xz_0_xxxzz_yzz, g_xz_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyzz_xx[k] = -g_xz_0_xxxzz_xx[k] * ab_y + g_xz_0_xxxzz_xxy[k];

                g_xz_0_xxxyzz_xy[k] = -g_xz_0_xxxzz_xy[k] * ab_y + g_xz_0_xxxzz_xyy[k];

                g_xz_0_xxxyzz_xz[k] = -g_xz_0_xxxzz_xz[k] * ab_y + g_xz_0_xxxzz_xyz[k];

                g_xz_0_xxxyzz_yy[k] = -g_xz_0_xxxzz_yy[k] * ab_y + g_xz_0_xxxzz_yyy[k];

                g_xz_0_xxxyzz_yz[k] = -g_xz_0_xxxzz_yz[k] * ab_y + g_xz_0_xxxzz_yyz[k];

                g_xz_0_xxxyzz_zz[k] = -g_xz_0_xxxzz_zz[k] * ab_y + g_xz_0_xxxzz_yzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 391 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzzz_xx, g_xz_0_xxxzzz_xy, g_xz_0_xxxzzz_xz, g_xz_0_xxxzzz_yy, g_xz_0_xxxzzz_yz, g_xz_0_xxxzzz_zz, g_xz_0_xxzzz_xx, g_xz_0_xxzzz_xxx, g_xz_0_xxzzz_xxy, g_xz_0_xxzzz_xxz, g_xz_0_xxzzz_xy, g_xz_0_xxzzz_xyy, g_xz_0_xxzzz_xyz, g_xz_0_xxzzz_xz, g_xz_0_xxzzz_xzz, g_xz_0_xxzzz_yy, g_xz_0_xxzzz_yz, g_xz_0_xxzzz_zz, g_z_0_xxzzz_xx, g_z_0_xxzzz_xy, g_z_0_xxzzz_xz, g_z_0_xxzzz_yy, g_z_0_xxzzz_yz, g_z_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzzz_xx[k] = -g_z_0_xxzzz_xx[k] - g_xz_0_xxzzz_xx[k] * ab_x + g_xz_0_xxzzz_xxx[k];

                g_xz_0_xxxzzz_xy[k] = -g_z_0_xxzzz_xy[k] - g_xz_0_xxzzz_xy[k] * ab_x + g_xz_0_xxzzz_xxy[k];

                g_xz_0_xxxzzz_xz[k] = -g_z_0_xxzzz_xz[k] - g_xz_0_xxzzz_xz[k] * ab_x + g_xz_0_xxzzz_xxz[k];

                g_xz_0_xxxzzz_yy[k] = -g_z_0_xxzzz_yy[k] - g_xz_0_xxzzz_yy[k] * ab_x + g_xz_0_xxzzz_xyy[k];

                g_xz_0_xxxzzz_yz[k] = -g_z_0_xxzzz_yz[k] - g_xz_0_xxzzz_yz[k] * ab_x + g_xz_0_xxzzz_xyz[k];

                g_xz_0_xxxzzz_zz[k] = -g_z_0_xxzzz_zz[k] - g_xz_0_xxzzz_zz[k] * ab_x + g_xz_0_xxzzz_xzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 399 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyy_xx, g_xz_0_xxyyy_xxy, g_xz_0_xxyyy_xy, g_xz_0_xxyyy_xyy, g_xz_0_xxyyy_xyz, g_xz_0_xxyyy_xz, g_xz_0_xxyyy_yy, g_xz_0_xxyyy_yyy, g_xz_0_xxyyy_yyz, g_xz_0_xxyyy_yz, g_xz_0_xxyyy_yzz, g_xz_0_xxyyy_zz, g_xz_0_xxyyyy_xx, g_xz_0_xxyyyy_xy, g_xz_0_xxyyyy_xz, g_xz_0_xxyyyy_yy, g_xz_0_xxyyyy_yz, g_xz_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyy_xx[k] = -g_xz_0_xxyyy_xx[k] * ab_y + g_xz_0_xxyyy_xxy[k];

                g_xz_0_xxyyyy_xy[k] = -g_xz_0_xxyyy_xy[k] * ab_y + g_xz_0_xxyyy_xyy[k];

                g_xz_0_xxyyyy_xz[k] = -g_xz_0_xxyyy_xz[k] * ab_y + g_xz_0_xxyyy_xyz[k];

                g_xz_0_xxyyyy_yy[k] = -g_xz_0_xxyyy_yy[k] * ab_y + g_xz_0_xxyyy_yyy[k];

                g_xz_0_xxyyyy_yz[k] = -g_xz_0_xxyyy_yz[k] * ab_y + g_xz_0_xxyyy_yyz[k];

                g_xz_0_xxyyyy_zz[k] = -g_xz_0_xxyyy_zz[k] * ab_y + g_xz_0_xxyyy_yzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 404 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyyz_xx, g_xz_0_xxyyyz_xy, g_xz_0_xxyyyz_xz, g_xz_0_xxyyyz_yy, g_xz_0_xxyyyz_yz, g_xz_0_xxyyyz_zz, g_xz_0_xxyyz_xx, g_xz_0_xxyyz_xxy, g_xz_0_xxyyz_xy, g_xz_0_xxyyz_xyy, g_xz_0_xxyyz_xyz, g_xz_0_xxyyz_xz, g_xz_0_xxyyz_yy, g_xz_0_xxyyz_yyy, g_xz_0_xxyyz_yyz, g_xz_0_xxyyz_yz, g_xz_0_xxyyz_yzz, g_xz_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyz_xx[k] = -g_xz_0_xxyyz_xx[k] * ab_y + g_xz_0_xxyyz_xxy[k];

                g_xz_0_xxyyyz_xy[k] = -g_xz_0_xxyyz_xy[k] * ab_y + g_xz_0_xxyyz_xyy[k];

                g_xz_0_xxyyyz_xz[k] = -g_xz_0_xxyyz_xz[k] * ab_y + g_xz_0_xxyyz_xyz[k];

                g_xz_0_xxyyyz_yy[k] = -g_xz_0_xxyyz_yy[k] * ab_y + g_xz_0_xxyyz_yyy[k];

                g_xz_0_xxyyyz_yz[k] = -g_xz_0_xxyyz_yz[k] * ab_y + g_xz_0_xxyyz_yyz[k];

                g_xz_0_xxyyyz_zz[k] = -g_xz_0_xxyyz_zz[k] * ab_y + g_xz_0_xxyyz_yzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 409 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyzz_xx, g_xz_0_xxyyzz_xy, g_xz_0_xxyyzz_xz, g_xz_0_xxyyzz_yy, g_xz_0_xxyyzz_yz, g_xz_0_xxyyzz_zz, g_xz_0_xxyzz_xx, g_xz_0_xxyzz_xxy, g_xz_0_xxyzz_xy, g_xz_0_xxyzz_xyy, g_xz_0_xxyzz_xyz, g_xz_0_xxyzz_xz, g_xz_0_xxyzz_yy, g_xz_0_xxyzz_yyy, g_xz_0_xxyzz_yyz, g_xz_0_xxyzz_yz, g_xz_0_xxyzz_yzz, g_xz_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyzz_xx[k] = -g_xz_0_xxyzz_xx[k] * ab_y + g_xz_0_xxyzz_xxy[k];

                g_xz_0_xxyyzz_xy[k] = -g_xz_0_xxyzz_xy[k] * ab_y + g_xz_0_xxyzz_xyy[k];

                g_xz_0_xxyyzz_xz[k] = -g_xz_0_xxyzz_xz[k] * ab_y + g_xz_0_xxyzz_xyz[k];

                g_xz_0_xxyyzz_yy[k] = -g_xz_0_xxyzz_yy[k] * ab_y + g_xz_0_xxyzz_yyy[k];

                g_xz_0_xxyyzz_yz[k] = -g_xz_0_xxyzz_yz[k] * ab_y + g_xz_0_xxyzz_yyz[k];

                g_xz_0_xxyyzz_zz[k] = -g_xz_0_xxyzz_zz[k] * ab_y + g_xz_0_xxyzz_yzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzzz_xx, g_xz_0_xxyzzz_xy, g_xz_0_xxyzzz_xz, g_xz_0_xxyzzz_yy, g_xz_0_xxyzzz_yz, g_xz_0_xxyzzz_zz, g_xz_0_xxzzz_xx, g_xz_0_xxzzz_xxy, g_xz_0_xxzzz_xy, g_xz_0_xxzzz_xyy, g_xz_0_xxzzz_xyz, g_xz_0_xxzzz_xz, g_xz_0_xxzzz_yy, g_xz_0_xxzzz_yyy, g_xz_0_xxzzz_yyz, g_xz_0_xxzzz_yz, g_xz_0_xxzzz_yzz, g_xz_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzzz_xx[k] = -g_xz_0_xxzzz_xx[k] * ab_y + g_xz_0_xxzzz_xxy[k];

                g_xz_0_xxyzzz_xy[k] = -g_xz_0_xxzzz_xy[k] * ab_y + g_xz_0_xxzzz_xyy[k];

                g_xz_0_xxyzzz_xz[k] = -g_xz_0_xxzzz_xz[k] * ab_y + g_xz_0_xxzzz_xyz[k];

                g_xz_0_xxyzzz_yy[k] = -g_xz_0_xxzzz_yy[k] * ab_y + g_xz_0_xxzzz_yyy[k];

                g_xz_0_xxyzzz_yz[k] = -g_xz_0_xxzzz_yz[k] * ab_y + g_xz_0_xxzzz_yyz[k];

                g_xz_0_xxyzzz_zz[k] = -g_xz_0_xxzzz_zz[k] * ab_y + g_xz_0_xxzzz_yzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzzz_xx, g_xz_0_xxzzzz_xy, g_xz_0_xxzzzz_xz, g_xz_0_xxzzzz_yy, g_xz_0_xxzzzz_yz, g_xz_0_xxzzzz_zz, g_xz_0_xzzzz_xx, g_xz_0_xzzzz_xxx, g_xz_0_xzzzz_xxy, g_xz_0_xzzzz_xxz, g_xz_0_xzzzz_xy, g_xz_0_xzzzz_xyy, g_xz_0_xzzzz_xyz, g_xz_0_xzzzz_xz, g_xz_0_xzzzz_xzz, g_xz_0_xzzzz_yy, g_xz_0_xzzzz_yz, g_xz_0_xzzzz_zz, g_z_0_xzzzz_xx, g_z_0_xzzzz_xy, g_z_0_xzzzz_xz, g_z_0_xzzzz_yy, g_z_0_xzzzz_yz, g_z_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzzz_xx[k] = -g_z_0_xzzzz_xx[k] - g_xz_0_xzzzz_xx[k] * ab_x + g_xz_0_xzzzz_xxx[k];

                g_xz_0_xxzzzz_xy[k] = -g_z_0_xzzzz_xy[k] - g_xz_0_xzzzz_xy[k] * ab_x + g_xz_0_xzzzz_xxy[k];

                g_xz_0_xxzzzz_xz[k] = -g_z_0_xzzzz_xz[k] - g_xz_0_xzzzz_xz[k] * ab_x + g_xz_0_xzzzz_xxz[k];

                g_xz_0_xxzzzz_yy[k] = -g_z_0_xzzzz_yy[k] - g_xz_0_xzzzz_yy[k] * ab_x + g_xz_0_xzzzz_xyy[k];

                g_xz_0_xxzzzz_yz[k] = -g_z_0_xzzzz_yz[k] - g_xz_0_xzzzz_yz[k] * ab_x + g_xz_0_xzzzz_xyz[k];

                g_xz_0_xxzzzz_zz[k] = -g_z_0_xzzzz_zz[k] - g_xz_0_xzzzz_zz[k] * ab_x + g_xz_0_xzzzz_xzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyy_xx, g_xz_0_xyyyy_xxy, g_xz_0_xyyyy_xy, g_xz_0_xyyyy_xyy, g_xz_0_xyyyy_xyz, g_xz_0_xyyyy_xz, g_xz_0_xyyyy_yy, g_xz_0_xyyyy_yyy, g_xz_0_xyyyy_yyz, g_xz_0_xyyyy_yz, g_xz_0_xyyyy_yzz, g_xz_0_xyyyy_zz, g_xz_0_xyyyyy_xx, g_xz_0_xyyyyy_xy, g_xz_0_xyyyyy_xz, g_xz_0_xyyyyy_yy, g_xz_0_xyyyyy_yz, g_xz_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyy_xx[k] = -g_xz_0_xyyyy_xx[k] * ab_y + g_xz_0_xyyyy_xxy[k];

                g_xz_0_xyyyyy_xy[k] = -g_xz_0_xyyyy_xy[k] * ab_y + g_xz_0_xyyyy_xyy[k];

                g_xz_0_xyyyyy_xz[k] = -g_xz_0_xyyyy_xz[k] * ab_y + g_xz_0_xyyyy_xyz[k];

                g_xz_0_xyyyyy_yy[k] = -g_xz_0_xyyyy_yy[k] * ab_y + g_xz_0_xyyyy_yyy[k];

                g_xz_0_xyyyyy_yz[k] = -g_xz_0_xyyyy_yz[k] * ab_y + g_xz_0_xyyyy_yyz[k];

                g_xz_0_xyyyyy_zz[k] = -g_xz_0_xyyyy_zz[k] * ab_y + g_xz_0_xyyyy_yzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyyz_xx, g_xz_0_xyyyyz_xy, g_xz_0_xyyyyz_xz, g_xz_0_xyyyyz_yy, g_xz_0_xyyyyz_yz, g_xz_0_xyyyyz_zz, g_xz_0_xyyyz_xx, g_xz_0_xyyyz_xxy, g_xz_0_xyyyz_xy, g_xz_0_xyyyz_xyy, g_xz_0_xyyyz_xyz, g_xz_0_xyyyz_xz, g_xz_0_xyyyz_yy, g_xz_0_xyyyz_yyy, g_xz_0_xyyyz_yyz, g_xz_0_xyyyz_yz, g_xz_0_xyyyz_yzz, g_xz_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyz_xx[k] = -g_xz_0_xyyyz_xx[k] * ab_y + g_xz_0_xyyyz_xxy[k];

                g_xz_0_xyyyyz_xy[k] = -g_xz_0_xyyyz_xy[k] * ab_y + g_xz_0_xyyyz_xyy[k];

                g_xz_0_xyyyyz_xz[k] = -g_xz_0_xyyyz_xz[k] * ab_y + g_xz_0_xyyyz_xyz[k];

                g_xz_0_xyyyyz_yy[k] = -g_xz_0_xyyyz_yy[k] * ab_y + g_xz_0_xyyyz_yyy[k];

                g_xz_0_xyyyyz_yz[k] = -g_xz_0_xyyyz_yz[k] * ab_y + g_xz_0_xyyyz_yyz[k];

                g_xz_0_xyyyyz_zz[k] = -g_xz_0_xyyyz_zz[k] * ab_y + g_xz_0_xyyyz_yzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyzz_xx, g_xz_0_xyyyzz_xy, g_xz_0_xyyyzz_xz, g_xz_0_xyyyzz_yy, g_xz_0_xyyyzz_yz, g_xz_0_xyyyzz_zz, g_xz_0_xyyzz_xx, g_xz_0_xyyzz_xxy, g_xz_0_xyyzz_xy, g_xz_0_xyyzz_xyy, g_xz_0_xyyzz_xyz, g_xz_0_xyyzz_xz, g_xz_0_xyyzz_yy, g_xz_0_xyyzz_yyy, g_xz_0_xyyzz_yyz, g_xz_0_xyyzz_yz, g_xz_0_xyyzz_yzz, g_xz_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyzz_xx[k] = -g_xz_0_xyyzz_xx[k] * ab_y + g_xz_0_xyyzz_xxy[k];

                g_xz_0_xyyyzz_xy[k] = -g_xz_0_xyyzz_xy[k] * ab_y + g_xz_0_xyyzz_xyy[k];

                g_xz_0_xyyyzz_xz[k] = -g_xz_0_xyyzz_xz[k] * ab_y + g_xz_0_xyyzz_xyz[k];

                g_xz_0_xyyyzz_yy[k] = -g_xz_0_xyyzz_yy[k] * ab_y + g_xz_0_xyyzz_yyy[k];

                g_xz_0_xyyyzz_yz[k] = -g_xz_0_xyyzz_yz[k] * ab_y + g_xz_0_xyyzz_yyz[k];

                g_xz_0_xyyyzz_zz[k] = -g_xz_0_xyyzz_zz[k] * ab_y + g_xz_0_xyyzz_yzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzzz_xx, g_xz_0_xyyzzz_xy, g_xz_0_xyyzzz_xz, g_xz_0_xyyzzz_yy, g_xz_0_xyyzzz_yz, g_xz_0_xyyzzz_zz, g_xz_0_xyzzz_xx, g_xz_0_xyzzz_xxy, g_xz_0_xyzzz_xy, g_xz_0_xyzzz_xyy, g_xz_0_xyzzz_xyz, g_xz_0_xyzzz_xz, g_xz_0_xyzzz_yy, g_xz_0_xyzzz_yyy, g_xz_0_xyzzz_yyz, g_xz_0_xyzzz_yz, g_xz_0_xyzzz_yzz, g_xz_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzzz_xx[k] = -g_xz_0_xyzzz_xx[k] * ab_y + g_xz_0_xyzzz_xxy[k];

                g_xz_0_xyyzzz_xy[k] = -g_xz_0_xyzzz_xy[k] * ab_y + g_xz_0_xyzzz_xyy[k];

                g_xz_0_xyyzzz_xz[k] = -g_xz_0_xyzzz_xz[k] * ab_y + g_xz_0_xyzzz_xyz[k];

                g_xz_0_xyyzzz_yy[k] = -g_xz_0_xyzzz_yy[k] * ab_y + g_xz_0_xyzzz_yyy[k];

                g_xz_0_xyyzzz_yz[k] = -g_xz_0_xyzzz_yz[k] * ab_y + g_xz_0_xyzzz_yyz[k];

                g_xz_0_xyyzzz_zz[k] = -g_xz_0_xyzzz_zz[k] * ab_y + g_xz_0_xyzzz_yzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 450 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 451 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 452 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 453 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 454 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzzz_xx, g_xz_0_xyzzzz_xy, g_xz_0_xyzzzz_xz, g_xz_0_xyzzzz_yy, g_xz_0_xyzzzz_yz, g_xz_0_xyzzzz_zz, g_xz_0_xzzzz_xx, g_xz_0_xzzzz_xxy, g_xz_0_xzzzz_xy, g_xz_0_xzzzz_xyy, g_xz_0_xzzzz_xyz, g_xz_0_xzzzz_xz, g_xz_0_xzzzz_yy, g_xz_0_xzzzz_yyy, g_xz_0_xzzzz_yyz, g_xz_0_xzzzz_yz, g_xz_0_xzzzz_yzz, g_xz_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzzz_xx[k] = -g_xz_0_xzzzz_xx[k] * ab_y + g_xz_0_xzzzz_xxy[k];

                g_xz_0_xyzzzz_xy[k] = -g_xz_0_xzzzz_xy[k] * ab_y + g_xz_0_xzzzz_xyy[k];

                g_xz_0_xyzzzz_xz[k] = -g_xz_0_xzzzz_xz[k] * ab_y + g_xz_0_xzzzz_xyz[k];

                g_xz_0_xyzzzz_yy[k] = -g_xz_0_xzzzz_yy[k] * ab_y + g_xz_0_xzzzz_yyy[k];

                g_xz_0_xyzzzz_yz[k] = -g_xz_0_xzzzz_yz[k] * ab_y + g_xz_0_xzzzz_yyz[k];

                g_xz_0_xyzzzz_zz[k] = -g_xz_0_xzzzz_zz[k] * ab_y + g_xz_0_xzzzz_yzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 456 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 457 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 458 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 459 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 460 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzzz_xx, g_xz_0_xzzzzz_xy, g_xz_0_xzzzzz_xz, g_xz_0_xzzzzz_yy, g_xz_0_xzzzzz_yz, g_xz_0_xzzzzz_zz, g_xz_0_zzzzz_xx, g_xz_0_zzzzz_xxx, g_xz_0_zzzzz_xxy, g_xz_0_zzzzz_xxz, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_xyy, g_xz_0_zzzzz_xyz, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_xzz, g_xz_0_zzzzz_yy, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_zz, g_z_0_zzzzz_xx, g_z_0_zzzzz_xy, g_z_0_zzzzz_xz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzzz_xx[k] = -g_z_0_zzzzz_xx[k] - g_xz_0_zzzzz_xx[k] * ab_x + g_xz_0_zzzzz_xxx[k];

                g_xz_0_xzzzzz_xy[k] = -g_z_0_zzzzz_xy[k] - g_xz_0_zzzzz_xy[k] * ab_x + g_xz_0_zzzzz_xxy[k];

                g_xz_0_xzzzzz_xz[k] = -g_z_0_zzzzz_xz[k] - g_xz_0_zzzzz_xz[k] * ab_x + g_xz_0_zzzzz_xxz[k];

                g_xz_0_xzzzzz_yy[k] = -g_z_0_zzzzz_yy[k] - g_xz_0_zzzzz_yy[k] * ab_x + g_xz_0_zzzzz_xyy[k];

                g_xz_0_xzzzzz_yz[k] = -g_z_0_zzzzz_yz[k] - g_xz_0_zzzzz_yz[k] * ab_x + g_xz_0_zzzzz_xyz[k];

                g_xz_0_xzzzzz_zz[k] = -g_z_0_zzzzz_zz[k] - g_xz_0_zzzzz_zz[k] * ab_x + g_xz_0_zzzzz_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 462 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 463 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 464 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 465 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 466 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyy_xx, g_xz_0_yyyyy_xxy, g_xz_0_yyyyy_xy, g_xz_0_yyyyy_xyy, g_xz_0_yyyyy_xyz, g_xz_0_yyyyy_xz, g_xz_0_yyyyy_yy, g_xz_0_yyyyy_yyy, g_xz_0_yyyyy_yyz, g_xz_0_yyyyy_yz, g_xz_0_yyyyy_yzz, g_xz_0_yyyyy_zz, g_xz_0_yyyyyy_xx, g_xz_0_yyyyyy_xy, g_xz_0_yyyyyy_xz, g_xz_0_yyyyyy_yy, g_xz_0_yyyyyy_yz, g_xz_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyy_xx[k] = -g_xz_0_yyyyy_xx[k] * ab_y + g_xz_0_yyyyy_xxy[k];

                g_xz_0_yyyyyy_xy[k] = -g_xz_0_yyyyy_xy[k] * ab_y + g_xz_0_yyyyy_xyy[k];

                g_xz_0_yyyyyy_xz[k] = -g_xz_0_yyyyy_xz[k] * ab_y + g_xz_0_yyyyy_xyz[k];

                g_xz_0_yyyyyy_yy[k] = -g_xz_0_yyyyy_yy[k] * ab_y + g_xz_0_yyyyy_yyy[k];

                g_xz_0_yyyyyy_yz[k] = -g_xz_0_yyyyy_yz[k] * ab_y + g_xz_0_yyyyy_yyz[k];

                g_xz_0_yyyyyy_zz[k] = -g_xz_0_yyyyy_zz[k] * ab_y + g_xz_0_yyyyy_yzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 468 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 469 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 470 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 471 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 472 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyyz_xx, g_xz_0_yyyyyz_xy, g_xz_0_yyyyyz_xz, g_xz_0_yyyyyz_yy, g_xz_0_yyyyyz_yz, g_xz_0_yyyyyz_zz, g_xz_0_yyyyz_xx, g_xz_0_yyyyz_xxy, g_xz_0_yyyyz_xy, g_xz_0_yyyyz_xyy, g_xz_0_yyyyz_xyz, g_xz_0_yyyyz_xz, g_xz_0_yyyyz_yy, g_xz_0_yyyyz_yyy, g_xz_0_yyyyz_yyz, g_xz_0_yyyyz_yz, g_xz_0_yyyyz_yzz, g_xz_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyz_xx[k] = -g_xz_0_yyyyz_xx[k] * ab_y + g_xz_0_yyyyz_xxy[k];

                g_xz_0_yyyyyz_xy[k] = -g_xz_0_yyyyz_xy[k] * ab_y + g_xz_0_yyyyz_xyy[k];

                g_xz_0_yyyyyz_xz[k] = -g_xz_0_yyyyz_xz[k] * ab_y + g_xz_0_yyyyz_xyz[k];

                g_xz_0_yyyyyz_yy[k] = -g_xz_0_yyyyz_yy[k] * ab_y + g_xz_0_yyyyz_yyy[k];

                g_xz_0_yyyyyz_yz[k] = -g_xz_0_yyyyz_yz[k] * ab_y + g_xz_0_yyyyz_yyz[k];

                g_xz_0_yyyyyz_zz[k] = -g_xz_0_yyyyz_zz[k] * ab_y + g_xz_0_yyyyz_yzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 474 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 475 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 476 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 477 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 478 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyzz_xx, g_xz_0_yyyyzz_xy, g_xz_0_yyyyzz_xz, g_xz_0_yyyyzz_yy, g_xz_0_yyyyzz_yz, g_xz_0_yyyyzz_zz, g_xz_0_yyyzz_xx, g_xz_0_yyyzz_xxy, g_xz_0_yyyzz_xy, g_xz_0_yyyzz_xyy, g_xz_0_yyyzz_xyz, g_xz_0_yyyzz_xz, g_xz_0_yyyzz_yy, g_xz_0_yyyzz_yyy, g_xz_0_yyyzz_yyz, g_xz_0_yyyzz_yz, g_xz_0_yyyzz_yzz, g_xz_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyzz_xx[k] = -g_xz_0_yyyzz_xx[k] * ab_y + g_xz_0_yyyzz_xxy[k];

                g_xz_0_yyyyzz_xy[k] = -g_xz_0_yyyzz_xy[k] * ab_y + g_xz_0_yyyzz_xyy[k];

                g_xz_0_yyyyzz_xz[k] = -g_xz_0_yyyzz_xz[k] * ab_y + g_xz_0_yyyzz_xyz[k];

                g_xz_0_yyyyzz_yy[k] = -g_xz_0_yyyzz_yy[k] * ab_y + g_xz_0_yyyzz_yyy[k];

                g_xz_0_yyyyzz_yz[k] = -g_xz_0_yyyzz_yz[k] * ab_y + g_xz_0_yyyzz_yyz[k];

                g_xz_0_yyyyzz_zz[k] = -g_xz_0_yyyzz_zz[k] * ab_y + g_xz_0_yyyzz_yzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 480 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 481 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 482 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 483 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 484 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzzz_xx, g_xz_0_yyyzzz_xy, g_xz_0_yyyzzz_xz, g_xz_0_yyyzzz_yy, g_xz_0_yyyzzz_yz, g_xz_0_yyyzzz_zz, g_xz_0_yyzzz_xx, g_xz_0_yyzzz_xxy, g_xz_0_yyzzz_xy, g_xz_0_yyzzz_xyy, g_xz_0_yyzzz_xyz, g_xz_0_yyzzz_xz, g_xz_0_yyzzz_yy, g_xz_0_yyzzz_yyy, g_xz_0_yyzzz_yyz, g_xz_0_yyzzz_yz, g_xz_0_yyzzz_yzz, g_xz_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzzz_xx[k] = -g_xz_0_yyzzz_xx[k] * ab_y + g_xz_0_yyzzz_xxy[k];

                g_xz_0_yyyzzz_xy[k] = -g_xz_0_yyzzz_xy[k] * ab_y + g_xz_0_yyzzz_xyy[k];

                g_xz_0_yyyzzz_xz[k] = -g_xz_0_yyzzz_xz[k] * ab_y + g_xz_0_yyzzz_xyz[k];

                g_xz_0_yyyzzz_yy[k] = -g_xz_0_yyzzz_yy[k] * ab_y + g_xz_0_yyzzz_yyy[k];

                g_xz_0_yyyzzz_yz[k] = -g_xz_0_yyzzz_yz[k] * ab_y + g_xz_0_yyzzz_yyz[k];

                g_xz_0_yyyzzz_zz[k] = -g_xz_0_yyzzz_zz[k] * ab_y + g_xz_0_yyzzz_yzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 486 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 487 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 488 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 489 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 490 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzzz_xx, g_xz_0_yyzzzz_xy, g_xz_0_yyzzzz_xz, g_xz_0_yyzzzz_yy, g_xz_0_yyzzzz_yz, g_xz_0_yyzzzz_zz, g_xz_0_yzzzz_xx, g_xz_0_yzzzz_xxy, g_xz_0_yzzzz_xy, g_xz_0_yzzzz_xyy, g_xz_0_yzzzz_xyz, g_xz_0_yzzzz_xz, g_xz_0_yzzzz_yy, g_xz_0_yzzzz_yyy, g_xz_0_yzzzz_yyz, g_xz_0_yzzzz_yz, g_xz_0_yzzzz_yzz, g_xz_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzzz_xx[k] = -g_xz_0_yzzzz_xx[k] * ab_y + g_xz_0_yzzzz_xxy[k];

                g_xz_0_yyzzzz_xy[k] = -g_xz_0_yzzzz_xy[k] * ab_y + g_xz_0_yzzzz_xyy[k];

                g_xz_0_yyzzzz_xz[k] = -g_xz_0_yzzzz_xz[k] * ab_y + g_xz_0_yzzzz_xyz[k];

                g_xz_0_yyzzzz_yy[k] = -g_xz_0_yzzzz_yy[k] * ab_y + g_xz_0_yzzzz_yyy[k];

                g_xz_0_yyzzzz_yz[k] = -g_xz_0_yzzzz_yz[k] * ab_y + g_xz_0_yzzzz_yyz[k];

                g_xz_0_yyzzzz_zz[k] = -g_xz_0_yzzzz_zz[k] * ab_y + g_xz_0_yzzzz_yzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 492 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 493 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 494 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 495 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 496 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzzz_xx, g_xz_0_yzzzzz_xy, g_xz_0_yzzzzz_xz, g_xz_0_yzzzzz_yy, g_xz_0_yzzzzz_yz, g_xz_0_yzzzzz_zz, g_xz_0_zzzzz_xx, g_xz_0_zzzzz_xxy, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_xyy, g_xz_0_zzzzz_xyz, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_yy, g_xz_0_zzzzz_yyy, g_xz_0_zzzzz_yyz, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_yzz, g_xz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzzz_xx[k] = -g_xz_0_zzzzz_xx[k] * ab_y + g_xz_0_zzzzz_xxy[k];

                g_xz_0_yzzzzz_xy[k] = -g_xz_0_zzzzz_xy[k] * ab_y + g_xz_0_zzzzz_xyy[k];

                g_xz_0_yzzzzz_xz[k] = -g_xz_0_zzzzz_xz[k] * ab_y + g_xz_0_zzzzz_xyz[k];

                g_xz_0_yzzzzz_yy[k] = -g_xz_0_zzzzz_yy[k] * ab_y + g_xz_0_zzzzz_yyy[k];

                g_xz_0_yzzzzz_yz[k] = -g_xz_0_zzzzz_yz[k] * ab_y + g_xz_0_zzzzz_yyz[k];

                g_xz_0_yzzzzz_zz[k] = -g_xz_0_zzzzz_zz[k] * ab_y + g_xz_0_zzzzz_yzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 498 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 499 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 500 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 501 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 502 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_xx, g_x_0_zzzzz_xy, g_x_0_zzzzz_xz, g_x_0_zzzzz_yy, g_x_0_zzzzz_yz, g_x_0_zzzzz_zz, g_xz_0_zzzzz_xx, g_xz_0_zzzzz_xxz, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_xyz, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_xzz, g_xz_0_zzzzz_yy, g_xz_0_zzzzz_yyz, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_yzz, g_xz_0_zzzzz_zz, g_xz_0_zzzzz_zzz, g_xz_0_zzzzzz_xx, g_xz_0_zzzzzz_xy, g_xz_0_zzzzzz_xz, g_xz_0_zzzzzz_yy, g_xz_0_zzzzzz_yz, g_xz_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzzz_xx[k] = -g_x_0_zzzzz_xx[k] - g_xz_0_zzzzz_xx[k] * ab_z + g_xz_0_zzzzz_xxz[k];

                g_xz_0_zzzzzz_xy[k] = -g_x_0_zzzzz_xy[k] - g_xz_0_zzzzz_xy[k] * ab_z + g_xz_0_zzzzz_xyz[k];

                g_xz_0_zzzzzz_xz[k] = -g_x_0_zzzzz_xz[k] - g_xz_0_zzzzz_xz[k] * ab_z + g_xz_0_zzzzz_xzz[k];

                g_xz_0_zzzzzz_yy[k] = -g_x_0_zzzzz_yy[k] - g_xz_0_zzzzz_yy[k] * ab_z + g_xz_0_zzzzz_yyz[k];

                g_xz_0_zzzzzz_yz[k] = -g_x_0_zzzzz_yz[k] - g_xz_0_zzzzz_yz[k] * ab_z + g_xz_0_zzzzz_yzz[k];

                g_xz_0_zzzzzz_zz[k] = -g_x_0_zzzzz_zz[k] - g_xz_0_zzzzz_zz[k] * ab_z + g_xz_0_zzzzz_zzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxx_xx, g_yy_0_xxxxx_xxx, g_yy_0_xxxxx_xxy, g_yy_0_xxxxx_xxz, g_yy_0_xxxxx_xy, g_yy_0_xxxxx_xyy, g_yy_0_xxxxx_xyz, g_yy_0_xxxxx_xz, g_yy_0_xxxxx_xzz, g_yy_0_xxxxx_yy, g_yy_0_xxxxx_yz, g_yy_0_xxxxx_zz, g_yy_0_xxxxxx_xx, g_yy_0_xxxxxx_xy, g_yy_0_xxxxxx_xz, g_yy_0_xxxxxx_yy, g_yy_0_xxxxxx_yz, g_yy_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxx_xx[k] = -g_yy_0_xxxxx_xx[k] * ab_x + g_yy_0_xxxxx_xxx[k];

                g_yy_0_xxxxxx_xy[k] = -g_yy_0_xxxxx_xy[k] * ab_x + g_yy_0_xxxxx_xxy[k];

                g_yy_0_xxxxxx_xz[k] = -g_yy_0_xxxxx_xz[k] * ab_x + g_yy_0_xxxxx_xxz[k];

                g_yy_0_xxxxxx_yy[k] = -g_yy_0_xxxxx_yy[k] * ab_x + g_yy_0_xxxxx_xyy[k];

                g_yy_0_xxxxxx_yz[k] = -g_yy_0_xxxxx_yz[k] * ab_x + g_yy_0_xxxxx_xyz[k];

                g_yy_0_xxxxxx_zz[k] = -g_yy_0_xxxxx_zz[k] * ab_x + g_yy_0_xxxxx_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxy_xx, g_yy_0_xxxxxy_xy, g_yy_0_xxxxxy_xz, g_yy_0_xxxxxy_yy, g_yy_0_xxxxxy_yz, g_yy_0_xxxxxy_zz, g_yy_0_xxxxy_xx, g_yy_0_xxxxy_xxx, g_yy_0_xxxxy_xxy, g_yy_0_xxxxy_xxz, g_yy_0_xxxxy_xy, g_yy_0_xxxxy_xyy, g_yy_0_xxxxy_xyz, g_yy_0_xxxxy_xz, g_yy_0_xxxxy_xzz, g_yy_0_xxxxy_yy, g_yy_0_xxxxy_yz, g_yy_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxy_xx[k] = -g_yy_0_xxxxy_xx[k] * ab_x + g_yy_0_xxxxy_xxx[k];

                g_yy_0_xxxxxy_xy[k] = -g_yy_0_xxxxy_xy[k] * ab_x + g_yy_0_xxxxy_xxy[k];

                g_yy_0_xxxxxy_xz[k] = -g_yy_0_xxxxy_xz[k] * ab_x + g_yy_0_xxxxy_xxz[k];

                g_yy_0_xxxxxy_yy[k] = -g_yy_0_xxxxy_yy[k] * ab_x + g_yy_0_xxxxy_xyy[k];

                g_yy_0_xxxxxy_yz[k] = -g_yy_0_xxxxy_yz[k] * ab_x + g_yy_0_xxxxy_xyz[k];

                g_yy_0_xxxxxy_zz[k] = -g_yy_0_xxxxy_zz[k] * ab_x + g_yy_0_xxxxy_xzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 519 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxz_xx, g_yy_0_xxxxxz_xy, g_yy_0_xxxxxz_xz, g_yy_0_xxxxxz_yy, g_yy_0_xxxxxz_yz, g_yy_0_xxxxxz_zz, g_yy_0_xxxxz_xx, g_yy_0_xxxxz_xxx, g_yy_0_xxxxz_xxy, g_yy_0_xxxxz_xxz, g_yy_0_xxxxz_xy, g_yy_0_xxxxz_xyy, g_yy_0_xxxxz_xyz, g_yy_0_xxxxz_xz, g_yy_0_xxxxz_xzz, g_yy_0_xxxxz_yy, g_yy_0_xxxxz_yz, g_yy_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxz_xx[k] = -g_yy_0_xxxxz_xx[k] * ab_x + g_yy_0_xxxxz_xxx[k];

                g_yy_0_xxxxxz_xy[k] = -g_yy_0_xxxxz_xy[k] * ab_x + g_yy_0_xxxxz_xxy[k];

                g_yy_0_xxxxxz_xz[k] = -g_yy_0_xxxxz_xz[k] * ab_x + g_yy_0_xxxxz_xxz[k];

                g_yy_0_xxxxxz_yy[k] = -g_yy_0_xxxxz_yy[k] * ab_x + g_yy_0_xxxxz_xyy[k];

                g_yy_0_xxxxxz_yz[k] = -g_yy_0_xxxxz_yz[k] * ab_x + g_yy_0_xxxxz_xyz[k];

                g_yy_0_xxxxxz_zz[k] = -g_yy_0_xxxxz_zz[k] * ab_x + g_yy_0_xxxxz_xzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 524 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyy_xx, g_yy_0_xxxxyy_xy, g_yy_0_xxxxyy_xz, g_yy_0_xxxxyy_yy, g_yy_0_xxxxyy_yz, g_yy_0_xxxxyy_zz, g_yy_0_xxxyy_xx, g_yy_0_xxxyy_xxx, g_yy_0_xxxyy_xxy, g_yy_0_xxxyy_xxz, g_yy_0_xxxyy_xy, g_yy_0_xxxyy_xyy, g_yy_0_xxxyy_xyz, g_yy_0_xxxyy_xz, g_yy_0_xxxyy_xzz, g_yy_0_xxxyy_yy, g_yy_0_xxxyy_yz, g_yy_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyy_xx[k] = -g_yy_0_xxxyy_xx[k] * ab_x + g_yy_0_xxxyy_xxx[k];

                g_yy_0_xxxxyy_xy[k] = -g_yy_0_xxxyy_xy[k] * ab_x + g_yy_0_xxxyy_xxy[k];

                g_yy_0_xxxxyy_xz[k] = -g_yy_0_xxxyy_xz[k] * ab_x + g_yy_0_xxxyy_xxz[k];

                g_yy_0_xxxxyy_yy[k] = -g_yy_0_xxxyy_yy[k] * ab_x + g_yy_0_xxxyy_xyy[k];

                g_yy_0_xxxxyy_yz[k] = -g_yy_0_xxxyy_yz[k] * ab_x + g_yy_0_xxxyy_xyz[k];

                g_yy_0_xxxxyy_zz[k] = -g_yy_0_xxxyy_zz[k] * ab_x + g_yy_0_xxxyy_xzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 529 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 531 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyz_xx, g_yy_0_xxxxyz_xy, g_yy_0_xxxxyz_xz, g_yy_0_xxxxyz_yy, g_yy_0_xxxxyz_yz, g_yy_0_xxxxyz_zz, g_yy_0_xxxyz_xx, g_yy_0_xxxyz_xxx, g_yy_0_xxxyz_xxy, g_yy_0_xxxyz_xxz, g_yy_0_xxxyz_xy, g_yy_0_xxxyz_xyy, g_yy_0_xxxyz_xyz, g_yy_0_xxxyz_xz, g_yy_0_xxxyz_xzz, g_yy_0_xxxyz_yy, g_yy_0_xxxyz_yz, g_yy_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyz_xx[k] = -g_yy_0_xxxyz_xx[k] * ab_x + g_yy_0_xxxyz_xxx[k];

                g_yy_0_xxxxyz_xy[k] = -g_yy_0_xxxyz_xy[k] * ab_x + g_yy_0_xxxyz_xxy[k];

                g_yy_0_xxxxyz_xz[k] = -g_yy_0_xxxyz_xz[k] * ab_x + g_yy_0_xxxyz_xxz[k];

                g_yy_0_xxxxyz_yy[k] = -g_yy_0_xxxyz_yy[k] * ab_x + g_yy_0_xxxyz_xyy[k];

                g_yy_0_xxxxyz_yz[k] = -g_yy_0_xxxyz_yz[k] * ab_x + g_yy_0_xxxyz_xyz[k];

                g_yy_0_xxxxyz_zz[k] = -g_yy_0_xxxyz_zz[k] * ab_x + g_yy_0_xxxyz_xzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxzz_xx, g_yy_0_xxxxzz_xy, g_yy_0_xxxxzz_xz, g_yy_0_xxxxzz_yy, g_yy_0_xxxxzz_yz, g_yy_0_xxxxzz_zz, g_yy_0_xxxzz_xx, g_yy_0_xxxzz_xxx, g_yy_0_xxxzz_xxy, g_yy_0_xxxzz_xxz, g_yy_0_xxxzz_xy, g_yy_0_xxxzz_xyy, g_yy_0_xxxzz_xyz, g_yy_0_xxxzz_xz, g_yy_0_xxxzz_xzz, g_yy_0_xxxzz_yy, g_yy_0_xxxzz_yz, g_yy_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxzz_xx[k] = -g_yy_0_xxxzz_xx[k] * ab_x + g_yy_0_xxxzz_xxx[k];

                g_yy_0_xxxxzz_xy[k] = -g_yy_0_xxxzz_xy[k] * ab_x + g_yy_0_xxxzz_xxy[k];

                g_yy_0_xxxxzz_xz[k] = -g_yy_0_xxxzz_xz[k] * ab_x + g_yy_0_xxxzz_xxz[k];

                g_yy_0_xxxxzz_yy[k] = -g_yy_0_xxxzz_yy[k] * ab_x + g_yy_0_xxxzz_xyy[k];

                g_yy_0_xxxxzz_yz[k] = -g_yy_0_xxxzz_yz[k] * ab_x + g_yy_0_xxxzz_xyz[k];

                g_yy_0_xxxxzz_zz[k] = -g_yy_0_xxxzz_zz[k] * ab_x + g_yy_0_xxxzz_xzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyy_xx, g_yy_0_xxxyyy_xy, g_yy_0_xxxyyy_xz, g_yy_0_xxxyyy_yy, g_yy_0_xxxyyy_yz, g_yy_0_xxxyyy_zz, g_yy_0_xxyyy_xx, g_yy_0_xxyyy_xxx, g_yy_0_xxyyy_xxy, g_yy_0_xxyyy_xxz, g_yy_0_xxyyy_xy, g_yy_0_xxyyy_xyy, g_yy_0_xxyyy_xyz, g_yy_0_xxyyy_xz, g_yy_0_xxyyy_xzz, g_yy_0_xxyyy_yy, g_yy_0_xxyyy_yz, g_yy_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyy_xx[k] = -g_yy_0_xxyyy_xx[k] * ab_x + g_yy_0_xxyyy_xxx[k];

                g_yy_0_xxxyyy_xy[k] = -g_yy_0_xxyyy_xy[k] * ab_x + g_yy_0_xxyyy_xxy[k];

                g_yy_0_xxxyyy_xz[k] = -g_yy_0_xxyyy_xz[k] * ab_x + g_yy_0_xxyyy_xxz[k];

                g_yy_0_xxxyyy_yy[k] = -g_yy_0_xxyyy_yy[k] * ab_x + g_yy_0_xxyyy_xyy[k];

                g_yy_0_xxxyyy_yz[k] = -g_yy_0_xxyyy_yz[k] * ab_x + g_yy_0_xxyyy_xyz[k];

                g_yy_0_xxxyyy_zz[k] = -g_yy_0_xxyyy_zz[k] * ab_x + g_yy_0_xxyyy_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 549 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyz_xx, g_yy_0_xxxyyz_xy, g_yy_0_xxxyyz_xz, g_yy_0_xxxyyz_yy, g_yy_0_xxxyyz_yz, g_yy_0_xxxyyz_zz, g_yy_0_xxyyz_xx, g_yy_0_xxyyz_xxx, g_yy_0_xxyyz_xxy, g_yy_0_xxyyz_xxz, g_yy_0_xxyyz_xy, g_yy_0_xxyyz_xyy, g_yy_0_xxyyz_xyz, g_yy_0_xxyyz_xz, g_yy_0_xxyyz_xzz, g_yy_0_xxyyz_yy, g_yy_0_xxyyz_yz, g_yy_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyz_xx[k] = -g_yy_0_xxyyz_xx[k] * ab_x + g_yy_0_xxyyz_xxx[k];

                g_yy_0_xxxyyz_xy[k] = -g_yy_0_xxyyz_xy[k] * ab_x + g_yy_0_xxyyz_xxy[k];

                g_yy_0_xxxyyz_xz[k] = -g_yy_0_xxyyz_xz[k] * ab_x + g_yy_0_xxyyz_xxz[k];

                g_yy_0_xxxyyz_yy[k] = -g_yy_0_xxyyz_yy[k] * ab_x + g_yy_0_xxyyz_xyy[k];

                g_yy_0_xxxyyz_yz[k] = -g_yy_0_xxyyz_yz[k] * ab_x + g_yy_0_xxyyz_xyz[k];

                g_yy_0_xxxyyz_zz[k] = -g_yy_0_xxyyz_zz[k] * ab_x + g_yy_0_xxyyz_xzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 554 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyzz_xx, g_yy_0_xxxyzz_xy, g_yy_0_xxxyzz_xz, g_yy_0_xxxyzz_yy, g_yy_0_xxxyzz_yz, g_yy_0_xxxyzz_zz, g_yy_0_xxyzz_xx, g_yy_0_xxyzz_xxx, g_yy_0_xxyzz_xxy, g_yy_0_xxyzz_xxz, g_yy_0_xxyzz_xy, g_yy_0_xxyzz_xyy, g_yy_0_xxyzz_xyz, g_yy_0_xxyzz_xz, g_yy_0_xxyzz_xzz, g_yy_0_xxyzz_yy, g_yy_0_xxyzz_yz, g_yy_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyzz_xx[k] = -g_yy_0_xxyzz_xx[k] * ab_x + g_yy_0_xxyzz_xxx[k];

                g_yy_0_xxxyzz_xy[k] = -g_yy_0_xxyzz_xy[k] * ab_x + g_yy_0_xxyzz_xxy[k];

                g_yy_0_xxxyzz_xz[k] = -g_yy_0_xxyzz_xz[k] * ab_x + g_yy_0_xxyzz_xxz[k];

                g_yy_0_xxxyzz_yy[k] = -g_yy_0_xxyzz_yy[k] * ab_x + g_yy_0_xxyzz_xyy[k];

                g_yy_0_xxxyzz_yz[k] = -g_yy_0_xxyzz_yz[k] * ab_x + g_yy_0_xxyzz_xyz[k];

                g_yy_0_xxxyzz_zz[k] = -g_yy_0_xxyzz_zz[k] * ab_x + g_yy_0_xxyzz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 559 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzzz_xx, g_yy_0_xxxzzz_xy, g_yy_0_xxxzzz_xz, g_yy_0_xxxzzz_yy, g_yy_0_xxxzzz_yz, g_yy_0_xxxzzz_zz, g_yy_0_xxzzz_xx, g_yy_0_xxzzz_xxx, g_yy_0_xxzzz_xxy, g_yy_0_xxzzz_xxz, g_yy_0_xxzzz_xy, g_yy_0_xxzzz_xyy, g_yy_0_xxzzz_xyz, g_yy_0_xxzzz_xz, g_yy_0_xxzzz_xzz, g_yy_0_xxzzz_yy, g_yy_0_xxzzz_yz, g_yy_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzzz_xx[k] = -g_yy_0_xxzzz_xx[k] * ab_x + g_yy_0_xxzzz_xxx[k];

                g_yy_0_xxxzzz_xy[k] = -g_yy_0_xxzzz_xy[k] * ab_x + g_yy_0_xxzzz_xxy[k];

                g_yy_0_xxxzzz_xz[k] = -g_yy_0_xxzzz_xz[k] * ab_x + g_yy_0_xxzzz_xxz[k];

                g_yy_0_xxxzzz_yy[k] = -g_yy_0_xxzzz_yy[k] * ab_x + g_yy_0_xxzzz_xyy[k];

                g_yy_0_xxxzzz_yz[k] = -g_yy_0_xxzzz_yz[k] * ab_x + g_yy_0_xxzzz_xyz[k];

                g_yy_0_xxxzzz_zz[k] = -g_yy_0_xxzzz_zz[k] * ab_x + g_yy_0_xxzzz_xzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyy_xx, g_yy_0_xxyyyy_xy, g_yy_0_xxyyyy_xz, g_yy_0_xxyyyy_yy, g_yy_0_xxyyyy_yz, g_yy_0_xxyyyy_zz, g_yy_0_xyyyy_xx, g_yy_0_xyyyy_xxx, g_yy_0_xyyyy_xxy, g_yy_0_xyyyy_xxz, g_yy_0_xyyyy_xy, g_yy_0_xyyyy_xyy, g_yy_0_xyyyy_xyz, g_yy_0_xyyyy_xz, g_yy_0_xyyyy_xzz, g_yy_0_xyyyy_yy, g_yy_0_xyyyy_yz, g_yy_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyy_xx[k] = -g_yy_0_xyyyy_xx[k] * ab_x + g_yy_0_xyyyy_xxx[k];

                g_yy_0_xxyyyy_xy[k] = -g_yy_0_xyyyy_xy[k] * ab_x + g_yy_0_xyyyy_xxy[k];

                g_yy_0_xxyyyy_xz[k] = -g_yy_0_xyyyy_xz[k] * ab_x + g_yy_0_xyyyy_xxz[k];

                g_yy_0_xxyyyy_yy[k] = -g_yy_0_xyyyy_yy[k] * ab_x + g_yy_0_xyyyy_xyy[k];

                g_yy_0_xxyyyy_yz[k] = -g_yy_0_xyyyy_yz[k] * ab_x + g_yy_0_xyyyy_xyz[k];

                g_yy_0_xxyyyy_zz[k] = -g_yy_0_xyyyy_zz[k] * ab_x + g_yy_0_xyyyy_xzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyz_xx, g_yy_0_xxyyyz_xy, g_yy_0_xxyyyz_xz, g_yy_0_xxyyyz_yy, g_yy_0_xxyyyz_yz, g_yy_0_xxyyyz_zz, g_yy_0_xyyyz_xx, g_yy_0_xyyyz_xxx, g_yy_0_xyyyz_xxy, g_yy_0_xyyyz_xxz, g_yy_0_xyyyz_xy, g_yy_0_xyyyz_xyy, g_yy_0_xyyyz_xyz, g_yy_0_xyyyz_xz, g_yy_0_xyyyz_xzz, g_yy_0_xyyyz_yy, g_yy_0_xyyyz_yz, g_yy_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyz_xx[k] = -g_yy_0_xyyyz_xx[k] * ab_x + g_yy_0_xyyyz_xxx[k];

                g_yy_0_xxyyyz_xy[k] = -g_yy_0_xyyyz_xy[k] * ab_x + g_yy_0_xyyyz_xxy[k];

                g_yy_0_xxyyyz_xz[k] = -g_yy_0_xyyyz_xz[k] * ab_x + g_yy_0_xyyyz_xxz[k];

                g_yy_0_xxyyyz_yy[k] = -g_yy_0_xyyyz_yy[k] * ab_x + g_yy_0_xyyyz_xyy[k];

                g_yy_0_xxyyyz_yz[k] = -g_yy_0_xyyyz_yz[k] * ab_x + g_yy_0_xyyyz_xyz[k];

                g_yy_0_xxyyyz_zz[k] = -g_yy_0_xyyyz_zz[k] * ab_x + g_yy_0_xyyyz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 579 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyzz_xx, g_yy_0_xxyyzz_xy, g_yy_0_xxyyzz_xz, g_yy_0_xxyyzz_yy, g_yy_0_xxyyzz_yz, g_yy_0_xxyyzz_zz, g_yy_0_xyyzz_xx, g_yy_0_xyyzz_xxx, g_yy_0_xyyzz_xxy, g_yy_0_xyyzz_xxz, g_yy_0_xyyzz_xy, g_yy_0_xyyzz_xyy, g_yy_0_xyyzz_xyz, g_yy_0_xyyzz_xz, g_yy_0_xyyzz_xzz, g_yy_0_xyyzz_yy, g_yy_0_xyyzz_yz, g_yy_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyzz_xx[k] = -g_yy_0_xyyzz_xx[k] * ab_x + g_yy_0_xyyzz_xxx[k];

                g_yy_0_xxyyzz_xy[k] = -g_yy_0_xyyzz_xy[k] * ab_x + g_yy_0_xyyzz_xxy[k];

                g_yy_0_xxyyzz_xz[k] = -g_yy_0_xyyzz_xz[k] * ab_x + g_yy_0_xyyzz_xxz[k];

                g_yy_0_xxyyzz_yy[k] = -g_yy_0_xyyzz_yy[k] * ab_x + g_yy_0_xyyzz_xyy[k];

                g_yy_0_xxyyzz_yz[k] = -g_yy_0_xyyzz_yz[k] * ab_x + g_yy_0_xyyzz_xyz[k];

                g_yy_0_xxyyzz_zz[k] = -g_yy_0_xyyzz_zz[k] * ab_x + g_yy_0_xyyzz_xzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 584 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzzz_xx, g_yy_0_xxyzzz_xy, g_yy_0_xxyzzz_xz, g_yy_0_xxyzzz_yy, g_yy_0_xxyzzz_yz, g_yy_0_xxyzzz_zz, g_yy_0_xyzzz_xx, g_yy_0_xyzzz_xxx, g_yy_0_xyzzz_xxy, g_yy_0_xyzzz_xxz, g_yy_0_xyzzz_xy, g_yy_0_xyzzz_xyy, g_yy_0_xyzzz_xyz, g_yy_0_xyzzz_xz, g_yy_0_xyzzz_xzz, g_yy_0_xyzzz_yy, g_yy_0_xyzzz_yz, g_yy_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzzz_xx[k] = -g_yy_0_xyzzz_xx[k] * ab_x + g_yy_0_xyzzz_xxx[k];

                g_yy_0_xxyzzz_xy[k] = -g_yy_0_xyzzz_xy[k] * ab_x + g_yy_0_xyzzz_xxy[k];

                g_yy_0_xxyzzz_xz[k] = -g_yy_0_xyzzz_xz[k] * ab_x + g_yy_0_xyzzz_xxz[k];

                g_yy_0_xxyzzz_yy[k] = -g_yy_0_xyzzz_yy[k] * ab_x + g_yy_0_xyzzz_xyy[k];

                g_yy_0_xxyzzz_yz[k] = -g_yy_0_xyzzz_yz[k] * ab_x + g_yy_0_xyzzz_xyz[k];

                g_yy_0_xxyzzz_zz[k] = -g_yy_0_xyzzz_zz[k] * ab_x + g_yy_0_xyzzz_xzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 589 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzzz_xx, g_yy_0_xxzzzz_xy, g_yy_0_xxzzzz_xz, g_yy_0_xxzzzz_yy, g_yy_0_xxzzzz_yz, g_yy_0_xxzzzz_zz, g_yy_0_xzzzz_xx, g_yy_0_xzzzz_xxx, g_yy_0_xzzzz_xxy, g_yy_0_xzzzz_xxz, g_yy_0_xzzzz_xy, g_yy_0_xzzzz_xyy, g_yy_0_xzzzz_xyz, g_yy_0_xzzzz_xz, g_yy_0_xzzzz_xzz, g_yy_0_xzzzz_yy, g_yy_0_xzzzz_yz, g_yy_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzzz_xx[k] = -g_yy_0_xzzzz_xx[k] * ab_x + g_yy_0_xzzzz_xxx[k];

                g_yy_0_xxzzzz_xy[k] = -g_yy_0_xzzzz_xy[k] * ab_x + g_yy_0_xzzzz_xxy[k];

                g_yy_0_xxzzzz_xz[k] = -g_yy_0_xzzzz_xz[k] * ab_x + g_yy_0_xzzzz_xxz[k];

                g_yy_0_xxzzzz_yy[k] = -g_yy_0_xzzzz_yy[k] * ab_x + g_yy_0_xzzzz_xyy[k];

                g_yy_0_xxzzzz_yz[k] = -g_yy_0_xzzzz_yz[k] * ab_x + g_yy_0_xzzzz_xyz[k];

                g_yy_0_xxzzzz_zz[k] = -g_yy_0_xzzzz_zz[k] * ab_x + g_yy_0_xzzzz_xzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyy_xx, g_yy_0_xyyyyy_xy, g_yy_0_xyyyyy_xz, g_yy_0_xyyyyy_yy, g_yy_0_xyyyyy_yz, g_yy_0_xyyyyy_zz, g_yy_0_yyyyy_xx, g_yy_0_yyyyy_xxx, g_yy_0_yyyyy_xxy, g_yy_0_yyyyy_xxz, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_xyy, g_yy_0_yyyyy_xyz, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_xzz, g_yy_0_yyyyy_yy, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyy_xx[k] = -g_yy_0_yyyyy_xx[k] * ab_x + g_yy_0_yyyyy_xxx[k];

                g_yy_0_xyyyyy_xy[k] = -g_yy_0_yyyyy_xy[k] * ab_x + g_yy_0_yyyyy_xxy[k];

                g_yy_0_xyyyyy_xz[k] = -g_yy_0_yyyyy_xz[k] * ab_x + g_yy_0_yyyyy_xxz[k];

                g_yy_0_xyyyyy_yy[k] = -g_yy_0_yyyyy_yy[k] * ab_x + g_yy_0_yyyyy_xyy[k];

                g_yy_0_xyyyyy_yz[k] = -g_yy_0_yyyyy_yz[k] * ab_x + g_yy_0_yyyyy_xyz[k];

                g_yy_0_xyyyyy_zz[k] = -g_yy_0_yyyyy_zz[k] * ab_x + g_yy_0_yyyyy_xzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 600 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 601 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 602 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 603 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 604 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyz_xx, g_yy_0_xyyyyz_xy, g_yy_0_xyyyyz_xz, g_yy_0_xyyyyz_yy, g_yy_0_xyyyyz_yz, g_yy_0_xyyyyz_zz, g_yy_0_yyyyz_xx, g_yy_0_yyyyz_xxx, g_yy_0_yyyyz_xxy, g_yy_0_yyyyz_xxz, g_yy_0_yyyyz_xy, g_yy_0_yyyyz_xyy, g_yy_0_yyyyz_xyz, g_yy_0_yyyyz_xz, g_yy_0_yyyyz_xzz, g_yy_0_yyyyz_yy, g_yy_0_yyyyz_yz, g_yy_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyz_xx[k] = -g_yy_0_yyyyz_xx[k] * ab_x + g_yy_0_yyyyz_xxx[k];

                g_yy_0_xyyyyz_xy[k] = -g_yy_0_yyyyz_xy[k] * ab_x + g_yy_0_yyyyz_xxy[k];

                g_yy_0_xyyyyz_xz[k] = -g_yy_0_yyyyz_xz[k] * ab_x + g_yy_0_yyyyz_xxz[k];

                g_yy_0_xyyyyz_yy[k] = -g_yy_0_yyyyz_yy[k] * ab_x + g_yy_0_yyyyz_xyy[k];

                g_yy_0_xyyyyz_yz[k] = -g_yy_0_yyyyz_yz[k] * ab_x + g_yy_0_yyyyz_xyz[k];

                g_yy_0_xyyyyz_zz[k] = -g_yy_0_yyyyz_zz[k] * ab_x + g_yy_0_yyyyz_xzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 606 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 607 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 608 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 609 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 610 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyzz_xx, g_yy_0_xyyyzz_xy, g_yy_0_xyyyzz_xz, g_yy_0_xyyyzz_yy, g_yy_0_xyyyzz_yz, g_yy_0_xyyyzz_zz, g_yy_0_yyyzz_xx, g_yy_0_yyyzz_xxx, g_yy_0_yyyzz_xxy, g_yy_0_yyyzz_xxz, g_yy_0_yyyzz_xy, g_yy_0_yyyzz_xyy, g_yy_0_yyyzz_xyz, g_yy_0_yyyzz_xz, g_yy_0_yyyzz_xzz, g_yy_0_yyyzz_yy, g_yy_0_yyyzz_yz, g_yy_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyzz_xx[k] = -g_yy_0_yyyzz_xx[k] * ab_x + g_yy_0_yyyzz_xxx[k];

                g_yy_0_xyyyzz_xy[k] = -g_yy_0_yyyzz_xy[k] * ab_x + g_yy_0_yyyzz_xxy[k];

                g_yy_0_xyyyzz_xz[k] = -g_yy_0_yyyzz_xz[k] * ab_x + g_yy_0_yyyzz_xxz[k];

                g_yy_0_xyyyzz_yy[k] = -g_yy_0_yyyzz_yy[k] * ab_x + g_yy_0_yyyzz_xyy[k];

                g_yy_0_xyyyzz_yz[k] = -g_yy_0_yyyzz_yz[k] * ab_x + g_yy_0_yyyzz_xyz[k];

                g_yy_0_xyyyzz_zz[k] = -g_yy_0_yyyzz_zz[k] * ab_x + g_yy_0_yyyzz_xzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 612 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 613 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 614 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 615 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 616 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzzz_xx, g_yy_0_xyyzzz_xy, g_yy_0_xyyzzz_xz, g_yy_0_xyyzzz_yy, g_yy_0_xyyzzz_yz, g_yy_0_xyyzzz_zz, g_yy_0_yyzzz_xx, g_yy_0_yyzzz_xxx, g_yy_0_yyzzz_xxy, g_yy_0_yyzzz_xxz, g_yy_0_yyzzz_xy, g_yy_0_yyzzz_xyy, g_yy_0_yyzzz_xyz, g_yy_0_yyzzz_xz, g_yy_0_yyzzz_xzz, g_yy_0_yyzzz_yy, g_yy_0_yyzzz_yz, g_yy_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzzz_xx[k] = -g_yy_0_yyzzz_xx[k] * ab_x + g_yy_0_yyzzz_xxx[k];

                g_yy_0_xyyzzz_xy[k] = -g_yy_0_yyzzz_xy[k] * ab_x + g_yy_0_yyzzz_xxy[k];

                g_yy_0_xyyzzz_xz[k] = -g_yy_0_yyzzz_xz[k] * ab_x + g_yy_0_yyzzz_xxz[k];

                g_yy_0_xyyzzz_yy[k] = -g_yy_0_yyzzz_yy[k] * ab_x + g_yy_0_yyzzz_xyy[k];

                g_yy_0_xyyzzz_yz[k] = -g_yy_0_yyzzz_yz[k] * ab_x + g_yy_0_yyzzz_xyz[k];

                g_yy_0_xyyzzz_zz[k] = -g_yy_0_yyzzz_zz[k] * ab_x + g_yy_0_yyzzz_xzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 618 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 619 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 620 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 621 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 622 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzzz_xx, g_yy_0_xyzzzz_xy, g_yy_0_xyzzzz_xz, g_yy_0_xyzzzz_yy, g_yy_0_xyzzzz_yz, g_yy_0_xyzzzz_zz, g_yy_0_yzzzz_xx, g_yy_0_yzzzz_xxx, g_yy_0_yzzzz_xxy, g_yy_0_yzzzz_xxz, g_yy_0_yzzzz_xy, g_yy_0_yzzzz_xyy, g_yy_0_yzzzz_xyz, g_yy_0_yzzzz_xz, g_yy_0_yzzzz_xzz, g_yy_0_yzzzz_yy, g_yy_0_yzzzz_yz, g_yy_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzzz_xx[k] = -g_yy_0_yzzzz_xx[k] * ab_x + g_yy_0_yzzzz_xxx[k];

                g_yy_0_xyzzzz_xy[k] = -g_yy_0_yzzzz_xy[k] * ab_x + g_yy_0_yzzzz_xxy[k];

                g_yy_0_xyzzzz_xz[k] = -g_yy_0_yzzzz_xz[k] * ab_x + g_yy_0_yzzzz_xxz[k];

                g_yy_0_xyzzzz_yy[k] = -g_yy_0_yzzzz_yy[k] * ab_x + g_yy_0_yzzzz_xyy[k];

                g_yy_0_xyzzzz_yz[k] = -g_yy_0_yzzzz_yz[k] * ab_x + g_yy_0_yzzzz_xyz[k];

                g_yy_0_xyzzzz_zz[k] = -g_yy_0_yzzzz_zz[k] * ab_x + g_yy_0_yzzzz_xzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 624 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 625 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 626 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 627 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 628 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzzz_xx, g_yy_0_xzzzzz_xy, g_yy_0_xzzzzz_xz, g_yy_0_xzzzzz_yy, g_yy_0_xzzzzz_yz, g_yy_0_xzzzzz_zz, g_yy_0_zzzzz_xx, g_yy_0_zzzzz_xxx, g_yy_0_zzzzz_xxy, g_yy_0_zzzzz_xxz, g_yy_0_zzzzz_xy, g_yy_0_zzzzz_xyy, g_yy_0_zzzzz_xyz, g_yy_0_zzzzz_xz, g_yy_0_zzzzz_xzz, g_yy_0_zzzzz_yy, g_yy_0_zzzzz_yz, g_yy_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzzz_xx[k] = -g_yy_0_zzzzz_xx[k] * ab_x + g_yy_0_zzzzz_xxx[k];

                g_yy_0_xzzzzz_xy[k] = -g_yy_0_zzzzz_xy[k] * ab_x + g_yy_0_zzzzz_xxy[k];

                g_yy_0_xzzzzz_xz[k] = -g_yy_0_zzzzz_xz[k] * ab_x + g_yy_0_zzzzz_xxz[k];

                g_yy_0_xzzzzz_yy[k] = -g_yy_0_zzzzz_yy[k] * ab_x + g_yy_0_zzzzz_xyy[k];

                g_yy_0_xzzzzz_yz[k] = -g_yy_0_zzzzz_yz[k] * ab_x + g_yy_0_zzzzz_xyz[k];

                g_yy_0_xzzzzz_zz[k] = -g_yy_0_zzzzz_zz[k] * ab_x + g_yy_0_zzzzz_xzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 630 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 631 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 632 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 633 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 634 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xx, g_y_0_yyyyy_xy, g_y_0_yyyyy_xz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yz, g_y_0_yyyyy_zz, g_yy_0_yyyyy_xx, g_yy_0_yyyyy_xxy, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_xyy, g_yy_0_yyyyy_xyz, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_yy, g_yy_0_yyyyy_yyy, g_yy_0_yyyyy_yyz, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_yzz, g_yy_0_yyyyy_zz, g_yy_0_yyyyyy_xx, g_yy_0_yyyyyy_xy, g_yy_0_yyyyyy_xz, g_yy_0_yyyyyy_yy, g_yy_0_yyyyyy_yz, g_yy_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyy_xx[k] = -2.0 * g_y_0_yyyyy_xx[k] - g_yy_0_yyyyy_xx[k] * ab_y + g_yy_0_yyyyy_xxy[k];

                g_yy_0_yyyyyy_xy[k] = -2.0 * g_y_0_yyyyy_xy[k] - g_yy_0_yyyyy_xy[k] * ab_y + g_yy_0_yyyyy_xyy[k];

                g_yy_0_yyyyyy_xz[k] = -2.0 * g_y_0_yyyyy_xz[k] - g_yy_0_yyyyy_xz[k] * ab_y + g_yy_0_yyyyy_xyz[k];

                g_yy_0_yyyyyy_yy[k] = -2.0 * g_y_0_yyyyy_yy[k] - g_yy_0_yyyyy_yy[k] * ab_y + g_yy_0_yyyyy_yyy[k];

                g_yy_0_yyyyyy_yz[k] = -2.0 * g_y_0_yyyyy_yz[k] - g_yy_0_yyyyy_yz[k] * ab_y + g_yy_0_yyyyy_yyz[k];

                g_yy_0_yyyyyy_zz[k] = -2.0 * g_y_0_yyyyy_zz[k] - g_yy_0_yyyyy_zz[k] * ab_y + g_yy_0_yyyyy_yzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 636 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 637 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 638 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 639 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 640 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyy_xx, g_yy_0_yyyyy_xxz, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_xyz, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_xzz, g_yy_0_yyyyy_yy, g_yy_0_yyyyy_yyz, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_yzz, g_yy_0_yyyyy_zz, g_yy_0_yyyyy_zzz, g_yy_0_yyyyyz_xx, g_yy_0_yyyyyz_xy, g_yy_0_yyyyyz_xz, g_yy_0_yyyyyz_yy, g_yy_0_yyyyyz_yz, g_yy_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyz_xx[k] = -g_yy_0_yyyyy_xx[k] * ab_z + g_yy_0_yyyyy_xxz[k];

                g_yy_0_yyyyyz_xy[k] = -g_yy_0_yyyyy_xy[k] * ab_z + g_yy_0_yyyyy_xyz[k];

                g_yy_0_yyyyyz_xz[k] = -g_yy_0_yyyyy_xz[k] * ab_z + g_yy_0_yyyyy_xzz[k];

                g_yy_0_yyyyyz_yy[k] = -g_yy_0_yyyyy_yy[k] * ab_z + g_yy_0_yyyyy_yyz[k];

                g_yy_0_yyyyyz_yz[k] = -g_yy_0_yyyyy_yz[k] * ab_z + g_yy_0_yyyyy_yzz[k];

                g_yy_0_yyyyyz_zz[k] = -g_yy_0_yyyyy_zz[k] * ab_z + g_yy_0_yyyyy_zzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 642 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 643 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 644 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 645 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 646 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyz_xx, g_yy_0_yyyyz_xxz, g_yy_0_yyyyz_xy, g_yy_0_yyyyz_xyz, g_yy_0_yyyyz_xz, g_yy_0_yyyyz_xzz, g_yy_0_yyyyz_yy, g_yy_0_yyyyz_yyz, g_yy_0_yyyyz_yz, g_yy_0_yyyyz_yzz, g_yy_0_yyyyz_zz, g_yy_0_yyyyz_zzz, g_yy_0_yyyyzz_xx, g_yy_0_yyyyzz_xy, g_yy_0_yyyyzz_xz, g_yy_0_yyyyzz_yy, g_yy_0_yyyyzz_yz, g_yy_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyzz_xx[k] = -g_yy_0_yyyyz_xx[k] * ab_z + g_yy_0_yyyyz_xxz[k];

                g_yy_0_yyyyzz_xy[k] = -g_yy_0_yyyyz_xy[k] * ab_z + g_yy_0_yyyyz_xyz[k];

                g_yy_0_yyyyzz_xz[k] = -g_yy_0_yyyyz_xz[k] * ab_z + g_yy_0_yyyyz_xzz[k];

                g_yy_0_yyyyzz_yy[k] = -g_yy_0_yyyyz_yy[k] * ab_z + g_yy_0_yyyyz_yyz[k];

                g_yy_0_yyyyzz_yz[k] = -g_yy_0_yyyyz_yz[k] * ab_z + g_yy_0_yyyyz_yzz[k];

                g_yy_0_yyyyzz_zz[k] = -g_yy_0_yyyyz_zz[k] * ab_z + g_yy_0_yyyyz_zzz[k];
            }

            /// Set up 648-654 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 648 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 649 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 650 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 651 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 652 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyzz_xx, g_yy_0_yyyzz_xxz, g_yy_0_yyyzz_xy, g_yy_0_yyyzz_xyz, g_yy_0_yyyzz_xz, g_yy_0_yyyzz_xzz, g_yy_0_yyyzz_yy, g_yy_0_yyyzz_yyz, g_yy_0_yyyzz_yz, g_yy_0_yyyzz_yzz, g_yy_0_yyyzz_zz, g_yy_0_yyyzz_zzz, g_yy_0_yyyzzz_xx, g_yy_0_yyyzzz_xy, g_yy_0_yyyzzz_xz, g_yy_0_yyyzzz_yy, g_yy_0_yyyzzz_yz, g_yy_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzzz_xx[k] = -g_yy_0_yyyzz_xx[k] * ab_z + g_yy_0_yyyzz_xxz[k];

                g_yy_0_yyyzzz_xy[k] = -g_yy_0_yyyzz_xy[k] * ab_z + g_yy_0_yyyzz_xyz[k];

                g_yy_0_yyyzzz_xz[k] = -g_yy_0_yyyzz_xz[k] * ab_z + g_yy_0_yyyzz_xzz[k];

                g_yy_0_yyyzzz_yy[k] = -g_yy_0_yyyzz_yy[k] * ab_z + g_yy_0_yyyzz_yyz[k];

                g_yy_0_yyyzzz_yz[k] = -g_yy_0_yyyzz_yz[k] * ab_z + g_yy_0_yyyzz_yzz[k];

                g_yy_0_yyyzzz_zz[k] = -g_yy_0_yyyzz_zz[k] * ab_z + g_yy_0_yyyzz_zzz[k];
            }

            /// Set up 654-660 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 654 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 655 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 656 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 657 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 658 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzzz_xx, g_yy_0_yyzzz_xxz, g_yy_0_yyzzz_xy, g_yy_0_yyzzz_xyz, g_yy_0_yyzzz_xz, g_yy_0_yyzzz_xzz, g_yy_0_yyzzz_yy, g_yy_0_yyzzz_yyz, g_yy_0_yyzzz_yz, g_yy_0_yyzzz_yzz, g_yy_0_yyzzz_zz, g_yy_0_yyzzz_zzz, g_yy_0_yyzzzz_xx, g_yy_0_yyzzzz_xy, g_yy_0_yyzzzz_xz, g_yy_0_yyzzzz_yy, g_yy_0_yyzzzz_yz, g_yy_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzzz_xx[k] = -g_yy_0_yyzzz_xx[k] * ab_z + g_yy_0_yyzzz_xxz[k];

                g_yy_0_yyzzzz_xy[k] = -g_yy_0_yyzzz_xy[k] * ab_z + g_yy_0_yyzzz_xyz[k];

                g_yy_0_yyzzzz_xz[k] = -g_yy_0_yyzzz_xz[k] * ab_z + g_yy_0_yyzzz_xzz[k];

                g_yy_0_yyzzzz_yy[k] = -g_yy_0_yyzzz_yy[k] * ab_z + g_yy_0_yyzzz_yyz[k];

                g_yy_0_yyzzzz_yz[k] = -g_yy_0_yyzzz_yz[k] * ab_z + g_yy_0_yyzzz_yzz[k];

                g_yy_0_yyzzzz_zz[k] = -g_yy_0_yyzzz_zz[k] * ab_z + g_yy_0_yyzzz_zzz[k];
            }

            /// Set up 660-666 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 660 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 661 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 662 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 663 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 664 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzzz_xx, g_yy_0_yzzzz_xxz, g_yy_0_yzzzz_xy, g_yy_0_yzzzz_xyz, g_yy_0_yzzzz_xz, g_yy_0_yzzzz_xzz, g_yy_0_yzzzz_yy, g_yy_0_yzzzz_yyz, g_yy_0_yzzzz_yz, g_yy_0_yzzzz_yzz, g_yy_0_yzzzz_zz, g_yy_0_yzzzz_zzz, g_yy_0_yzzzzz_xx, g_yy_0_yzzzzz_xy, g_yy_0_yzzzzz_xz, g_yy_0_yzzzzz_yy, g_yy_0_yzzzzz_yz, g_yy_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzzz_xx[k] = -g_yy_0_yzzzz_xx[k] * ab_z + g_yy_0_yzzzz_xxz[k];

                g_yy_0_yzzzzz_xy[k] = -g_yy_0_yzzzz_xy[k] * ab_z + g_yy_0_yzzzz_xyz[k];

                g_yy_0_yzzzzz_xz[k] = -g_yy_0_yzzzz_xz[k] * ab_z + g_yy_0_yzzzz_xzz[k];

                g_yy_0_yzzzzz_yy[k] = -g_yy_0_yzzzz_yy[k] * ab_z + g_yy_0_yzzzz_yyz[k];

                g_yy_0_yzzzzz_yz[k] = -g_yy_0_yzzzz_yz[k] * ab_z + g_yy_0_yzzzz_yzz[k];

                g_yy_0_yzzzzz_zz[k] = -g_yy_0_yzzzz_zz[k] * ab_z + g_yy_0_yzzzz_zzz[k];
            }

            /// Set up 666-672 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 666 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 667 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 668 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 669 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 670 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzzz_xx, g_yy_0_zzzzz_xxz, g_yy_0_zzzzz_xy, g_yy_0_zzzzz_xyz, g_yy_0_zzzzz_xz, g_yy_0_zzzzz_xzz, g_yy_0_zzzzz_yy, g_yy_0_zzzzz_yyz, g_yy_0_zzzzz_yz, g_yy_0_zzzzz_yzz, g_yy_0_zzzzz_zz, g_yy_0_zzzzz_zzz, g_yy_0_zzzzzz_xx, g_yy_0_zzzzzz_xy, g_yy_0_zzzzzz_xz, g_yy_0_zzzzzz_yy, g_yy_0_zzzzzz_yz, g_yy_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzzz_xx[k] = -g_yy_0_zzzzz_xx[k] * ab_z + g_yy_0_zzzzz_xxz[k];

                g_yy_0_zzzzzz_xy[k] = -g_yy_0_zzzzz_xy[k] * ab_z + g_yy_0_zzzzz_xyz[k];

                g_yy_0_zzzzzz_xz[k] = -g_yy_0_zzzzz_xz[k] * ab_z + g_yy_0_zzzzz_xzz[k];

                g_yy_0_zzzzzz_yy[k] = -g_yy_0_zzzzz_yy[k] * ab_z + g_yy_0_zzzzz_yyz[k];

                g_yy_0_zzzzzz_yz[k] = -g_yy_0_zzzzz_yz[k] * ab_z + g_yy_0_zzzzz_yzz[k];

                g_yy_0_zzzzzz_zz[k] = -g_yy_0_zzzzz_zz[k] * ab_z + g_yy_0_zzzzz_zzz[k];
            }

            /// Set up 672-678 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 674 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxx_xx, g_yz_0_xxxxx_xxx, g_yz_0_xxxxx_xxy, g_yz_0_xxxxx_xxz, g_yz_0_xxxxx_xy, g_yz_0_xxxxx_xyy, g_yz_0_xxxxx_xyz, g_yz_0_xxxxx_xz, g_yz_0_xxxxx_xzz, g_yz_0_xxxxx_yy, g_yz_0_xxxxx_yz, g_yz_0_xxxxx_zz, g_yz_0_xxxxxx_xx, g_yz_0_xxxxxx_xy, g_yz_0_xxxxxx_xz, g_yz_0_xxxxxx_yy, g_yz_0_xxxxxx_yz, g_yz_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxx_xx[k] = -g_yz_0_xxxxx_xx[k] * ab_x + g_yz_0_xxxxx_xxx[k];

                g_yz_0_xxxxxx_xy[k] = -g_yz_0_xxxxx_xy[k] * ab_x + g_yz_0_xxxxx_xxy[k];

                g_yz_0_xxxxxx_xz[k] = -g_yz_0_xxxxx_xz[k] * ab_x + g_yz_0_xxxxx_xxz[k];

                g_yz_0_xxxxxx_yy[k] = -g_yz_0_xxxxx_yy[k] * ab_x + g_yz_0_xxxxx_xyy[k];

                g_yz_0_xxxxxx_yz[k] = -g_yz_0_xxxxx_yz[k] * ab_x + g_yz_0_xxxxx_xyz[k];

                g_yz_0_xxxxxx_zz[k] = -g_yz_0_xxxxx_zz[k] * ab_x + g_yz_0_xxxxx_xzz[k];
            }

            /// Set up 678-684 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 679 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxy_xx, g_yz_0_xxxxxy_xy, g_yz_0_xxxxxy_xz, g_yz_0_xxxxxy_yy, g_yz_0_xxxxxy_yz, g_yz_0_xxxxxy_zz, g_yz_0_xxxxy_xx, g_yz_0_xxxxy_xxx, g_yz_0_xxxxy_xxy, g_yz_0_xxxxy_xxz, g_yz_0_xxxxy_xy, g_yz_0_xxxxy_xyy, g_yz_0_xxxxy_xyz, g_yz_0_xxxxy_xz, g_yz_0_xxxxy_xzz, g_yz_0_xxxxy_yy, g_yz_0_xxxxy_yz, g_yz_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxy_xx[k] = -g_yz_0_xxxxy_xx[k] * ab_x + g_yz_0_xxxxy_xxx[k];

                g_yz_0_xxxxxy_xy[k] = -g_yz_0_xxxxy_xy[k] * ab_x + g_yz_0_xxxxy_xxy[k];

                g_yz_0_xxxxxy_xz[k] = -g_yz_0_xxxxy_xz[k] * ab_x + g_yz_0_xxxxy_xxz[k];

                g_yz_0_xxxxxy_yy[k] = -g_yz_0_xxxxy_yy[k] * ab_x + g_yz_0_xxxxy_xyy[k];

                g_yz_0_xxxxxy_yz[k] = -g_yz_0_xxxxy_yz[k] * ab_x + g_yz_0_xxxxy_xyz[k];

                g_yz_0_xxxxxy_zz[k] = -g_yz_0_xxxxy_zz[k] * ab_x + g_yz_0_xxxxy_xzz[k];
            }

            /// Set up 684-690 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxz_xx, g_yz_0_xxxxxz_xy, g_yz_0_xxxxxz_xz, g_yz_0_xxxxxz_yy, g_yz_0_xxxxxz_yz, g_yz_0_xxxxxz_zz, g_yz_0_xxxxz_xx, g_yz_0_xxxxz_xxx, g_yz_0_xxxxz_xxy, g_yz_0_xxxxz_xxz, g_yz_0_xxxxz_xy, g_yz_0_xxxxz_xyy, g_yz_0_xxxxz_xyz, g_yz_0_xxxxz_xz, g_yz_0_xxxxz_xzz, g_yz_0_xxxxz_yy, g_yz_0_xxxxz_yz, g_yz_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxz_xx[k] = -g_yz_0_xxxxz_xx[k] * ab_x + g_yz_0_xxxxz_xxx[k];

                g_yz_0_xxxxxz_xy[k] = -g_yz_0_xxxxz_xy[k] * ab_x + g_yz_0_xxxxz_xxy[k];

                g_yz_0_xxxxxz_xz[k] = -g_yz_0_xxxxz_xz[k] * ab_x + g_yz_0_xxxxz_xxz[k];

                g_yz_0_xxxxxz_yy[k] = -g_yz_0_xxxxz_yy[k] * ab_x + g_yz_0_xxxxz_xyy[k];

                g_yz_0_xxxxxz_yz[k] = -g_yz_0_xxxxz_yz[k] * ab_x + g_yz_0_xxxxz_xyz[k];

                g_yz_0_xxxxxz_zz[k] = -g_yz_0_xxxxz_zz[k] * ab_x + g_yz_0_xxxxz_xzz[k];
            }

            /// Set up 690-696 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyy_xx, g_yz_0_xxxxyy_xy, g_yz_0_xxxxyy_xz, g_yz_0_xxxxyy_yy, g_yz_0_xxxxyy_yz, g_yz_0_xxxxyy_zz, g_yz_0_xxxyy_xx, g_yz_0_xxxyy_xxx, g_yz_0_xxxyy_xxy, g_yz_0_xxxyy_xxz, g_yz_0_xxxyy_xy, g_yz_0_xxxyy_xyy, g_yz_0_xxxyy_xyz, g_yz_0_xxxyy_xz, g_yz_0_xxxyy_xzz, g_yz_0_xxxyy_yy, g_yz_0_xxxyy_yz, g_yz_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyy_xx[k] = -g_yz_0_xxxyy_xx[k] * ab_x + g_yz_0_xxxyy_xxx[k];

                g_yz_0_xxxxyy_xy[k] = -g_yz_0_xxxyy_xy[k] * ab_x + g_yz_0_xxxyy_xxy[k];

                g_yz_0_xxxxyy_xz[k] = -g_yz_0_xxxyy_xz[k] * ab_x + g_yz_0_xxxyy_xxz[k];

                g_yz_0_xxxxyy_yy[k] = -g_yz_0_xxxyy_yy[k] * ab_x + g_yz_0_xxxyy_xyy[k];

                g_yz_0_xxxxyy_yz[k] = -g_yz_0_xxxyy_yz[k] * ab_x + g_yz_0_xxxyy_xyz[k];

                g_yz_0_xxxxyy_zz[k] = -g_yz_0_xxxyy_zz[k] * ab_x + g_yz_0_xxxyy_xzz[k];
            }

            /// Set up 696-702 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 699 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyz_xx, g_yz_0_xxxxyz_xy, g_yz_0_xxxxyz_xz, g_yz_0_xxxxyz_yy, g_yz_0_xxxxyz_yz, g_yz_0_xxxxyz_zz, g_yz_0_xxxyz_xx, g_yz_0_xxxyz_xxx, g_yz_0_xxxyz_xxy, g_yz_0_xxxyz_xxz, g_yz_0_xxxyz_xy, g_yz_0_xxxyz_xyy, g_yz_0_xxxyz_xyz, g_yz_0_xxxyz_xz, g_yz_0_xxxyz_xzz, g_yz_0_xxxyz_yy, g_yz_0_xxxyz_yz, g_yz_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyz_xx[k] = -g_yz_0_xxxyz_xx[k] * ab_x + g_yz_0_xxxyz_xxx[k];

                g_yz_0_xxxxyz_xy[k] = -g_yz_0_xxxyz_xy[k] * ab_x + g_yz_0_xxxyz_xxy[k];

                g_yz_0_xxxxyz_xz[k] = -g_yz_0_xxxyz_xz[k] * ab_x + g_yz_0_xxxyz_xxz[k];

                g_yz_0_xxxxyz_yy[k] = -g_yz_0_xxxyz_yy[k] * ab_x + g_yz_0_xxxyz_xyy[k];

                g_yz_0_xxxxyz_yz[k] = -g_yz_0_xxxyz_yz[k] * ab_x + g_yz_0_xxxyz_xyz[k];

                g_yz_0_xxxxyz_zz[k] = -g_yz_0_xxxyz_zz[k] * ab_x + g_yz_0_xxxyz_xzz[k];
            }

            /// Set up 702-708 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 704 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxzz_xx, g_yz_0_xxxxzz_xy, g_yz_0_xxxxzz_xz, g_yz_0_xxxxzz_yy, g_yz_0_xxxxzz_yz, g_yz_0_xxxxzz_zz, g_yz_0_xxxzz_xx, g_yz_0_xxxzz_xxx, g_yz_0_xxxzz_xxy, g_yz_0_xxxzz_xxz, g_yz_0_xxxzz_xy, g_yz_0_xxxzz_xyy, g_yz_0_xxxzz_xyz, g_yz_0_xxxzz_xz, g_yz_0_xxxzz_xzz, g_yz_0_xxxzz_yy, g_yz_0_xxxzz_yz, g_yz_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxzz_xx[k] = -g_yz_0_xxxzz_xx[k] * ab_x + g_yz_0_xxxzz_xxx[k];

                g_yz_0_xxxxzz_xy[k] = -g_yz_0_xxxzz_xy[k] * ab_x + g_yz_0_xxxzz_xxy[k];

                g_yz_0_xxxxzz_xz[k] = -g_yz_0_xxxzz_xz[k] * ab_x + g_yz_0_xxxzz_xxz[k];

                g_yz_0_xxxxzz_yy[k] = -g_yz_0_xxxzz_yy[k] * ab_x + g_yz_0_xxxzz_xyy[k];

                g_yz_0_xxxxzz_yz[k] = -g_yz_0_xxxzz_yz[k] * ab_x + g_yz_0_xxxzz_xyz[k];

                g_yz_0_xxxxzz_zz[k] = -g_yz_0_xxxzz_zz[k] * ab_x + g_yz_0_xxxzz_xzz[k];
            }

            /// Set up 708-714 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 709 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyy_xx, g_yz_0_xxxyyy_xy, g_yz_0_xxxyyy_xz, g_yz_0_xxxyyy_yy, g_yz_0_xxxyyy_yz, g_yz_0_xxxyyy_zz, g_yz_0_xxyyy_xx, g_yz_0_xxyyy_xxx, g_yz_0_xxyyy_xxy, g_yz_0_xxyyy_xxz, g_yz_0_xxyyy_xy, g_yz_0_xxyyy_xyy, g_yz_0_xxyyy_xyz, g_yz_0_xxyyy_xz, g_yz_0_xxyyy_xzz, g_yz_0_xxyyy_yy, g_yz_0_xxyyy_yz, g_yz_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyy_xx[k] = -g_yz_0_xxyyy_xx[k] * ab_x + g_yz_0_xxyyy_xxx[k];

                g_yz_0_xxxyyy_xy[k] = -g_yz_0_xxyyy_xy[k] * ab_x + g_yz_0_xxyyy_xxy[k];

                g_yz_0_xxxyyy_xz[k] = -g_yz_0_xxyyy_xz[k] * ab_x + g_yz_0_xxyyy_xxz[k];

                g_yz_0_xxxyyy_yy[k] = -g_yz_0_xxyyy_yy[k] * ab_x + g_yz_0_xxyyy_xyy[k];

                g_yz_0_xxxyyy_yz[k] = -g_yz_0_xxyyy_yz[k] * ab_x + g_yz_0_xxyyy_xyz[k];

                g_yz_0_xxxyyy_zz[k] = -g_yz_0_xxyyy_zz[k] * ab_x + g_yz_0_xxyyy_xzz[k];
            }

            /// Set up 714-720 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyz_xx, g_yz_0_xxxyyz_xy, g_yz_0_xxxyyz_xz, g_yz_0_xxxyyz_yy, g_yz_0_xxxyyz_yz, g_yz_0_xxxyyz_zz, g_yz_0_xxyyz_xx, g_yz_0_xxyyz_xxx, g_yz_0_xxyyz_xxy, g_yz_0_xxyyz_xxz, g_yz_0_xxyyz_xy, g_yz_0_xxyyz_xyy, g_yz_0_xxyyz_xyz, g_yz_0_xxyyz_xz, g_yz_0_xxyyz_xzz, g_yz_0_xxyyz_yy, g_yz_0_xxyyz_yz, g_yz_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyz_xx[k] = -g_yz_0_xxyyz_xx[k] * ab_x + g_yz_0_xxyyz_xxx[k];

                g_yz_0_xxxyyz_xy[k] = -g_yz_0_xxyyz_xy[k] * ab_x + g_yz_0_xxyyz_xxy[k];

                g_yz_0_xxxyyz_xz[k] = -g_yz_0_xxyyz_xz[k] * ab_x + g_yz_0_xxyyz_xxz[k];

                g_yz_0_xxxyyz_yy[k] = -g_yz_0_xxyyz_yy[k] * ab_x + g_yz_0_xxyyz_xyy[k];

                g_yz_0_xxxyyz_yz[k] = -g_yz_0_xxyyz_yz[k] * ab_x + g_yz_0_xxyyz_xyz[k];

                g_yz_0_xxxyyz_zz[k] = -g_yz_0_xxyyz_zz[k] * ab_x + g_yz_0_xxyyz_xzz[k];
            }

            /// Set up 720-726 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyzz_xx, g_yz_0_xxxyzz_xy, g_yz_0_xxxyzz_xz, g_yz_0_xxxyzz_yy, g_yz_0_xxxyzz_yz, g_yz_0_xxxyzz_zz, g_yz_0_xxyzz_xx, g_yz_0_xxyzz_xxx, g_yz_0_xxyzz_xxy, g_yz_0_xxyzz_xxz, g_yz_0_xxyzz_xy, g_yz_0_xxyzz_xyy, g_yz_0_xxyzz_xyz, g_yz_0_xxyzz_xz, g_yz_0_xxyzz_xzz, g_yz_0_xxyzz_yy, g_yz_0_xxyzz_yz, g_yz_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyzz_xx[k] = -g_yz_0_xxyzz_xx[k] * ab_x + g_yz_0_xxyzz_xxx[k];

                g_yz_0_xxxyzz_xy[k] = -g_yz_0_xxyzz_xy[k] * ab_x + g_yz_0_xxyzz_xxy[k];

                g_yz_0_xxxyzz_xz[k] = -g_yz_0_xxyzz_xz[k] * ab_x + g_yz_0_xxyzz_xxz[k];

                g_yz_0_xxxyzz_yy[k] = -g_yz_0_xxyzz_yy[k] * ab_x + g_yz_0_xxyzz_xyy[k];

                g_yz_0_xxxyzz_yz[k] = -g_yz_0_xxyzz_yz[k] * ab_x + g_yz_0_xxyzz_xyz[k];

                g_yz_0_xxxyzz_zz[k] = -g_yz_0_xxyzz_zz[k] * ab_x + g_yz_0_xxyzz_xzz[k];
            }

            /// Set up 726-732 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 727 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 729 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzzz_xx, g_yz_0_xxxzzz_xy, g_yz_0_xxxzzz_xz, g_yz_0_xxxzzz_yy, g_yz_0_xxxzzz_yz, g_yz_0_xxxzzz_zz, g_yz_0_xxzzz_xx, g_yz_0_xxzzz_xxx, g_yz_0_xxzzz_xxy, g_yz_0_xxzzz_xxz, g_yz_0_xxzzz_xy, g_yz_0_xxzzz_xyy, g_yz_0_xxzzz_xyz, g_yz_0_xxzzz_xz, g_yz_0_xxzzz_xzz, g_yz_0_xxzzz_yy, g_yz_0_xxzzz_yz, g_yz_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzzz_xx[k] = -g_yz_0_xxzzz_xx[k] * ab_x + g_yz_0_xxzzz_xxx[k];

                g_yz_0_xxxzzz_xy[k] = -g_yz_0_xxzzz_xy[k] * ab_x + g_yz_0_xxzzz_xxy[k];

                g_yz_0_xxxzzz_xz[k] = -g_yz_0_xxzzz_xz[k] * ab_x + g_yz_0_xxzzz_xxz[k];

                g_yz_0_xxxzzz_yy[k] = -g_yz_0_xxzzz_yy[k] * ab_x + g_yz_0_xxzzz_xyy[k];

                g_yz_0_xxxzzz_yz[k] = -g_yz_0_xxzzz_yz[k] * ab_x + g_yz_0_xxzzz_xyz[k];

                g_yz_0_xxxzzz_zz[k] = -g_yz_0_xxzzz_zz[k] * ab_x + g_yz_0_xxzzz_xzz[k];
            }

            /// Set up 732-738 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 734 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyy_xx, g_yz_0_xxyyyy_xy, g_yz_0_xxyyyy_xz, g_yz_0_xxyyyy_yy, g_yz_0_xxyyyy_yz, g_yz_0_xxyyyy_zz, g_yz_0_xyyyy_xx, g_yz_0_xyyyy_xxx, g_yz_0_xyyyy_xxy, g_yz_0_xyyyy_xxz, g_yz_0_xyyyy_xy, g_yz_0_xyyyy_xyy, g_yz_0_xyyyy_xyz, g_yz_0_xyyyy_xz, g_yz_0_xyyyy_xzz, g_yz_0_xyyyy_yy, g_yz_0_xyyyy_yz, g_yz_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyy_xx[k] = -g_yz_0_xyyyy_xx[k] * ab_x + g_yz_0_xyyyy_xxx[k];

                g_yz_0_xxyyyy_xy[k] = -g_yz_0_xyyyy_xy[k] * ab_x + g_yz_0_xyyyy_xxy[k];

                g_yz_0_xxyyyy_xz[k] = -g_yz_0_xyyyy_xz[k] * ab_x + g_yz_0_xyyyy_xxz[k];

                g_yz_0_xxyyyy_yy[k] = -g_yz_0_xyyyy_yy[k] * ab_x + g_yz_0_xyyyy_xyy[k];

                g_yz_0_xxyyyy_yz[k] = -g_yz_0_xyyyy_yz[k] * ab_x + g_yz_0_xyyyy_xyz[k];

                g_yz_0_xxyyyy_zz[k] = -g_yz_0_xyyyy_zz[k] * ab_x + g_yz_0_xyyyy_xzz[k];
            }

            /// Set up 738-744 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 739 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyz_xx, g_yz_0_xxyyyz_xy, g_yz_0_xxyyyz_xz, g_yz_0_xxyyyz_yy, g_yz_0_xxyyyz_yz, g_yz_0_xxyyyz_zz, g_yz_0_xyyyz_xx, g_yz_0_xyyyz_xxx, g_yz_0_xyyyz_xxy, g_yz_0_xyyyz_xxz, g_yz_0_xyyyz_xy, g_yz_0_xyyyz_xyy, g_yz_0_xyyyz_xyz, g_yz_0_xyyyz_xz, g_yz_0_xyyyz_xzz, g_yz_0_xyyyz_yy, g_yz_0_xyyyz_yz, g_yz_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyz_xx[k] = -g_yz_0_xyyyz_xx[k] * ab_x + g_yz_0_xyyyz_xxx[k];

                g_yz_0_xxyyyz_xy[k] = -g_yz_0_xyyyz_xy[k] * ab_x + g_yz_0_xyyyz_xxy[k];

                g_yz_0_xxyyyz_xz[k] = -g_yz_0_xyyyz_xz[k] * ab_x + g_yz_0_xyyyz_xxz[k];

                g_yz_0_xxyyyz_yy[k] = -g_yz_0_xyyyz_yy[k] * ab_x + g_yz_0_xyyyz_xyy[k];

                g_yz_0_xxyyyz_yz[k] = -g_yz_0_xyyyz_yz[k] * ab_x + g_yz_0_xyyyz_xyz[k];

                g_yz_0_xxyyyz_zz[k] = -g_yz_0_xyyyz_zz[k] * ab_x + g_yz_0_xyyyz_xzz[k];
            }

            /// Set up 744-750 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyzz_xx, g_yz_0_xxyyzz_xy, g_yz_0_xxyyzz_xz, g_yz_0_xxyyzz_yy, g_yz_0_xxyyzz_yz, g_yz_0_xxyyzz_zz, g_yz_0_xyyzz_xx, g_yz_0_xyyzz_xxx, g_yz_0_xyyzz_xxy, g_yz_0_xyyzz_xxz, g_yz_0_xyyzz_xy, g_yz_0_xyyzz_xyy, g_yz_0_xyyzz_xyz, g_yz_0_xyyzz_xz, g_yz_0_xyyzz_xzz, g_yz_0_xyyzz_yy, g_yz_0_xyyzz_yz, g_yz_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyzz_xx[k] = -g_yz_0_xyyzz_xx[k] * ab_x + g_yz_0_xyyzz_xxx[k];

                g_yz_0_xxyyzz_xy[k] = -g_yz_0_xyyzz_xy[k] * ab_x + g_yz_0_xyyzz_xxy[k];

                g_yz_0_xxyyzz_xz[k] = -g_yz_0_xyyzz_xz[k] * ab_x + g_yz_0_xyyzz_xxz[k];

                g_yz_0_xxyyzz_yy[k] = -g_yz_0_xyyzz_yy[k] * ab_x + g_yz_0_xyyzz_xyy[k];

                g_yz_0_xxyyzz_yz[k] = -g_yz_0_xyyzz_yz[k] * ab_x + g_yz_0_xyyzz_xyz[k];

                g_yz_0_xxyyzz_zz[k] = -g_yz_0_xyyzz_zz[k] * ab_x + g_yz_0_xyyzz_xzz[k];
            }

            /// Set up 750-756 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 750 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 751 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 752 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 753 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 754 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzzz_xx, g_yz_0_xxyzzz_xy, g_yz_0_xxyzzz_xz, g_yz_0_xxyzzz_yy, g_yz_0_xxyzzz_yz, g_yz_0_xxyzzz_zz, g_yz_0_xyzzz_xx, g_yz_0_xyzzz_xxx, g_yz_0_xyzzz_xxy, g_yz_0_xyzzz_xxz, g_yz_0_xyzzz_xy, g_yz_0_xyzzz_xyy, g_yz_0_xyzzz_xyz, g_yz_0_xyzzz_xz, g_yz_0_xyzzz_xzz, g_yz_0_xyzzz_yy, g_yz_0_xyzzz_yz, g_yz_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzzz_xx[k] = -g_yz_0_xyzzz_xx[k] * ab_x + g_yz_0_xyzzz_xxx[k];

                g_yz_0_xxyzzz_xy[k] = -g_yz_0_xyzzz_xy[k] * ab_x + g_yz_0_xyzzz_xxy[k];

                g_yz_0_xxyzzz_xz[k] = -g_yz_0_xyzzz_xz[k] * ab_x + g_yz_0_xyzzz_xxz[k];

                g_yz_0_xxyzzz_yy[k] = -g_yz_0_xyzzz_yy[k] * ab_x + g_yz_0_xyzzz_xyy[k];

                g_yz_0_xxyzzz_yz[k] = -g_yz_0_xyzzz_yz[k] * ab_x + g_yz_0_xyzzz_xyz[k];

                g_yz_0_xxyzzz_zz[k] = -g_yz_0_xyzzz_zz[k] * ab_x + g_yz_0_xyzzz_xzz[k];
            }

            /// Set up 756-762 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 756 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 757 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 758 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 759 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 760 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 761 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzzz_xx, g_yz_0_xxzzzz_xy, g_yz_0_xxzzzz_xz, g_yz_0_xxzzzz_yy, g_yz_0_xxzzzz_yz, g_yz_0_xxzzzz_zz, g_yz_0_xzzzz_xx, g_yz_0_xzzzz_xxx, g_yz_0_xzzzz_xxy, g_yz_0_xzzzz_xxz, g_yz_0_xzzzz_xy, g_yz_0_xzzzz_xyy, g_yz_0_xzzzz_xyz, g_yz_0_xzzzz_xz, g_yz_0_xzzzz_xzz, g_yz_0_xzzzz_yy, g_yz_0_xzzzz_yz, g_yz_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzzz_xx[k] = -g_yz_0_xzzzz_xx[k] * ab_x + g_yz_0_xzzzz_xxx[k];

                g_yz_0_xxzzzz_xy[k] = -g_yz_0_xzzzz_xy[k] * ab_x + g_yz_0_xzzzz_xxy[k];

                g_yz_0_xxzzzz_xz[k] = -g_yz_0_xzzzz_xz[k] * ab_x + g_yz_0_xzzzz_xxz[k];

                g_yz_0_xxzzzz_yy[k] = -g_yz_0_xzzzz_yy[k] * ab_x + g_yz_0_xzzzz_xyy[k];

                g_yz_0_xxzzzz_yz[k] = -g_yz_0_xzzzz_yz[k] * ab_x + g_yz_0_xzzzz_xyz[k];

                g_yz_0_xxzzzz_zz[k] = -g_yz_0_xzzzz_zz[k] * ab_x + g_yz_0_xzzzz_xzz[k];
            }

            /// Set up 762-768 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 762 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 763 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 764 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 765 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 766 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 767 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyy_xx, g_yz_0_xyyyyy_xy, g_yz_0_xyyyyy_xz, g_yz_0_xyyyyy_yy, g_yz_0_xyyyyy_yz, g_yz_0_xyyyyy_zz, g_yz_0_yyyyy_xx, g_yz_0_yyyyy_xxx, g_yz_0_yyyyy_xxy, g_yz_0_yyyyy_xxz, g_yz_0_yyyyy_xy, g_yz_0_yyyyy_xyy, g_yz_0_yyyyy_xyz, g_yz_0_yyyyy_xz, g_yz_0_yyyyy_xzz, g_yz_0_yyyyy_yy, g_yz_0_yyyyy_yz, g_yz_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyy_xx[k] = -g_yz_0_yyyyy_xx[k] * ab_x + g_yz_0_yyyyy_xxx[k];

                g_yz_0_xyyyyy_xy[k] = -g_yz_0_yyyyy_xy[k] * ab_x + g_yz_0_yyyyy_xxy[k];

                g_yz_0_xyyyyy_xz[k] = -g_yz_0_yyyyy_xz[k] * ab_x + g_yz_0_yyyyy_xxz[k];

                g_yz_0_xyyyyy_yy[k] = -g_yz_0_yyyyy_yy[k] * ab_x + g_yz_0_yyyyy_xyy[k];

                g_yz_0_xyyyyy_yz[k] = -g_yz_0_yyyyy_yz[k] * ab_x + g_yz_0_yyyyy_xyz[k];

                g_yz_0_xyyyyy_zz[k] = -g_yz_0_yyyyy_zz[k] * ab_x + g_yz_0_yyyyy_xzz[k];
            }

            /// Set up 768-774 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 768 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 769 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 770 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 771 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 772 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 773 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyz_xx, g_yz_0_xyyyyz_xy, g_yz_0_xyyyyz_xz, g_yz_0_xyyyyz_yy, g_yz_0_xyyyyz_yz, g_yz_0_xyyyyz_zz, g_yz_0_yyyyz_xx, g_yz_0_yyyyz_xxx, g_yz_0_yyyyz_xxy, g_yz_0_yyyyz_xxz, g_yz_0_yyyyz_xy, g_yz_0_yyyyz_xyy, g_yz_0_yyyyz_xyz, g_yz_0_yyyyz_xz, g_yz_0_yyyyz_xzz, g_yz_0_yyyyz_yy, g_yz_0_yyyyz_yz, g_yz_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyz_xx[k] = -g_yz_0_yyyyz_xx[k] * ab_x + g_yz_0_yyyyz_xxx[k];

                g_yz_0_xyyyyz_xy[k] = -g_yz_0_yyyyz_xy[k] * ab_x + g_yz_0_yyyyz_xxy[k];

                g_yz_0_xyyyyz_xz[k] = -g_yz_0_yyyyz_xz[k] * ab_x + g_yz_0_yyyyz_xxz[k];

                g_yz_0_xyyyyz_yy[k] = -g_yz_0_yyyyz_yy[k] * ab_x + g_yz_0_yyyyz_xyy[k];

                g_yz_0_xyyyyz_yz[k] = -g_yz_0_yyyyz_yz[k] * ab_x + g_yz_0_yyyyz_xyz[k];

                g_yz_0_xyyyyz_zz[k] = -g_yz_0_yyyyz_zz[k] * ab_x + g_yz_0_yyyyz_xzz[k];
            }

            /// Set up 774-780 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 774 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 775 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 776 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 777 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 778 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyzz_xx, g_yz_0_xyyyzz_xy, g_yz_0_xyyyzz_xz, g_yz_0_xyyyzz_yy, g_yz_0_xyyyzz_yz, g_yz_0_xyyyzz_zz, g_yz_0_yyyzz_xx, g_yz_0_yyyzz_xxx, g_yz_0_yyyzz_xxy, g_yz_0_yyyzz_xxz, g_yz_0_yyyzz_xy, g_yz_0_yyyzz_xyy, g_yz_0_yyyzz_xyz, g_yz_0_yyyzz_xz, g_yz_0_yyyzz_xzz, g_yz_0_yyyzz_yy, g_yz_0_yyyzz_yz, g_yz_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyzz_xx[k] = -g_yz_0_yyyzz_xx[k] * ab_x + g_yz_0_yyyzz_xxx[k];

                g_yz_0_xyyyzz_xy[k] = -g_yz_0_yyyzz_xy[k] * ab_x + g_yz_0_yyyzz_xxy[k];

                g_yz_0_xyyyzz_xz[k] = -g_yz_0_yyyzz_xz[k] * ab_x + g_yz_0_yyyzz_xxz[k];

                g_yz_0_xyyyzz_yy[k] = -g_yz_0_yyyzz_yy[k] * ab_x + g_yz_0_yyyzz_xyy[k];

                g_yz_0_xyyyzz_yz[k] = -g_yz_0_yyyzz_yz[k] * ab_x + g_yz_0_yyyzz_xyz[k];

                g_yz_0_xyyyzz_zz[k] = -g_yz_0_yyyzz_zz[k] * ab_x + g_yz_0_yyyzz_xzz[k];
            }

            /// Set up 780-786 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 780 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 781 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 782 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 783 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 784 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 785 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzzz_xx, g_yz_0_xyyzzz_xy, g_yz_0_xyyzzz_xz, g_yz_0_xyyzzz_yy, g_yz_0_xyyzzz_yz, g_yz_0_xyyzzz_zz, g_yz_0_yyzzz_xx, g_yz_0_yyzzz_xxx, g_yz_0_yyzzz_xxy, g_yz_0_yyzzz_xxz, g_yz_0_yyzzz_xy, g_yz_0_yyzzz_xyy, g_yz_0_yyzzz_xyz, g_yz_0_yyzzz_xz, g_yz_0_yyzzz_xzz, g_yz_0_yyzzz_yy, g_yz_0_yyzzz_yz, g_yz_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzzz_xx[k] = -g_yz_0_yyzzz_xx[k] * ab_x + g_yz_0_yyzzz_xxx[k];

                g_yz_0_xyyzzz_xy[k] = -g_yz_0_yyzzz_xy[k] * ab_x + g_yz_0_yyzzz_xxy[k];

                g_yz_0_xyyzzz_xz[k] = -g_yz_0_yyzzz_xz[k] * ab_x + g_yz_0_yyzzz_xxz[k];

                g_yz_0_xyyzzz_yy[k] = -g_yz_0_yyzzz_yy[k] * ab_x + g_yz_0_yyzzz_xyy[k];

                g_yz_0_xyyzzz_yz[k] = -g_yz_0_yyzzz_yz[k] * ab_x + g_yz_0_yyzzz_xyz[k];

                g_yz_0_xyyzzz_zz[k] = -g_yz_0_yyzzz_zz[k] * ab_x + g_yz_0_yyzzz_xzz[k];
            }

            /// Set up 786-792 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 786 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 787 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 788 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 789 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 790 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 791 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzzz_xx, g_yz_0_xyzzzz_xy, g_yz_0_xyzzzz_xz, g_yz_0_xyzzzz_yy, g_yz_0_xyzzzz_yz, g_yz_0_xyzzzz_zz, g_yz_0_yzzzz_xx, g_yz_0_yzzzz_xxx, g_yz_0_yzzzz_xxy, g_yz_0_yzzzz_xxz, g_yz_0_yzzzz_xy, g_yz_0_yzzzz_xyy, g_yz_0_yzzzz_xyz, g_yz_0_yzzzz_xz, g_yz_0_yzzzz_xzz, g_yz_0_yzzzz_yy, g_yz_0_yzzzz_yz, g_yz_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzzz_xx[k] = -g_yz_0_yzzzz_xx[k] * ab_x + g_yz_0_yzzzz_xxx[k];

                g_yz_0_xyzzzz_xy[k] = -g_yz_0_yzzzz_xy[k] * ab_x + g_yz_0_yzzzz_xxy[k];

                g_yz_0_xyzzzz_xz[k] = -g_yz_0_yzzzz_xz[k] * ab_x + g_yz_0_yzzzz_xxz[k];

                g_yz_0_xyzzzz_yy[k] = -g_yz_0_yzzzz_yy[k] * ab_x + g_yz_0_yzzzz_xyy[k];

                g_yz_0_xyzzzz_yz[k] = -g_yz_0_yzzzz_yz[k] * ab_x + g_yz_0_yzzzz_xyz[k];

                g_yz_0_xyzzzz_zz[k] = -g_yz_0_yzzzz_zz[k] * ab_x + g_yz_0_yzzzz_xzz[k];
            }

            /// Set up 792-798 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 792 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 793 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 794 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 795 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 796 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzzz_xx, g_yz_0_xzzzzz_xy, g_yz_0_xzzzzz_xz, g_yz_0_xzzzzz_yy, g_yz_0_xzzzzz_yz, g_yz_0_xzzzzz_zz, g_yz_0_zzzzz_xx, g_yz_0_zzzzz_xxx, g_yz_0_zzzzz_xxy, g_yz_0_zzzzz_xxz, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_xyy, g_yz_0_zzzzz_xyz, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_xzz, g_yz_0_zzzzz_yy, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzzz_xx[k] = -g_yz_0_zzzzz_xx[k] * ab_x + g_yz_0_zzzzz_xxx[k];

                g_yz_0_xzzzzz_xy[k] = -g_yz_0_zzzzz_xy[k] * ab_x + g_yz_0_zzzzz_xxy[k];

                g_yz_0_xzzzzz_xz[k] = -g_yz_0_zzzzz_xz[k] * ab_x + g_yz_0_zzzzz_xxz[k];

                g_yz_0_xzzzzz_yy[k] = -g_yz_0_zzzzz_yy[k] * ab_x + g_yz_0_zzzzz_xyy[k];

                g_yz_0_xzzzzz_yz[k] = -g_yz_0_zzzzz_yz[k] * ab_x + g_yz_0_zzzzz_xyz[k];

                g_yz_0_xzzzzz_zz[k] = -g_yz_0_zzzzz_zz[k] * ab_x + g_yz_0_zzzzz_xzz[k];
            }

            /// Set up 798-804 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 798 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 799 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 800 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 801 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 802 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 803 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyy_xx, g_yz_0_yyyyy_xxy, g_yz_0_yyyyy_xy, g_yz_0_yyyyy_xyy, g_yz_0_yyyyy_xyz, g_yz_0_yyyyy_xz, g_yz_0_yyyyy_yy, g_yz_0_yyyyy_yyy, g_yz_0_yyyyy_yyz, g_yz_0_yyyyy_yz, g_yz_0_yyyyy_yzz, g_yz_0_yyyyy_zz, g_yz_0_yyyyyy_xx, g_yz_0_yyyyyy_xy, g_yz_0_yyyyyy_xz, g_yz_0_yyyyyy_yy, g_yz_0_yyyyyy_yz, g_yz_0_yyyyyy_zz, g_z_0_yyyyy_xx, g_z_0_yyyyy_xy, g_z_0_yyyyy_xz, g_z_0_yyyyy_yy, g_z_0_yyyyy_yz, g_z_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyy_xx[k] = -g_z_0_yyyyy_xx[k] - g_yz_0_yyyyy_xx[k] * ab_y + g_yz_0_yyyyy_xxy[k];

                g_yz_0_yyyyyy_xy[k] = -g_z_0_yyyyy_xy[k] - g_yz_0_yyyyy_xy[k] * ab_y + g_yz_0_yyyyy_xyy[k];

                g_yz_0_yyyyyy_xz[k] = -g_z_0_yyyyy_xz[k] - g_yz_0_yyyyy_xz[k] * ab_y + g_yz_0_yyyyy_xyz[k];

                g_yz_0_yyyyyy_yy[k] = -g_z_0_yyyyy_yy[k] - g_yz_0_yyyyy_yy[k] * ab_y + g_yz_0_yyyyy_yyy[k];

                g_yz_0_yyyyyy_yz[k] = -g_z_0_yyyyy_yz[k] - g_yz_0_yyyyy_yz[k] * ab_y + g_yz_0_yyyyy_yyz[k];

                g_yz_0_yyyyyy_zz[k] = -g_z_0_yyyyy_zz[k] - g_yz_0_yyyyy_zz[k] * ab_y + g_yz_0_yyyyy_yzz[k];
            }

            /// Set up 804-810 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 804 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 805 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 806 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 807 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 808 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyyz_xx, g_yz_0_yyyyyz_xy, g_yz_0_yyyyyz_xz, g_yz_0_yyyyyz_yy, g_yz_0_yyyyyz_yz, g_yz_0_yyyyyz_zz, g_yz_0_yyyyz_xx, g_yz_0_yyyyz_xxy, g_yz_0_yyyyz_xy, g_yz_0_yyyyz_xyy, g_yz_0_yyyyz_xyz, g_yz_0_yyyyz_xz, g_yz_0_yyyyz_yy, g_yz_0_yyyyz_yyy, g_yz_0_yyyyz_yyz, g_yz_0_yyyyz_yz, g_yz_0_yyyyz_yzz, g_yz_0_yyyyz_zz, g_z_0_yyyyz_xx, g_z_0_yyyyz_xy, g_z_0_yyyyz_xz, g_z_0_yyyyz_yy, g_z_0_yyyyz_yz, g_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyz_xx[k] = -g_z_0_yyyyz_xx[k] - g_yz_0_yyyyz_xx[k] * ab_y + g_yz_0_yyyyz_xxy[k];

                g_yz_0_yyyyyz_xy[k] = -g_z_0_yyyyz_xy[k] - g_yz_0_yyyyz_xy[k] * ab_y + g_yz_0_yyyyz_xyy[k];

                g_yz_0_yyyyyz_xz[k] = -g_z_0_yyyyz_xz[k] - g_yz_0_yyyyz_xz[k] * ab_y + g_yz_0_yyyyz_xyz[k];

                g_yz_0_yyyyyz_yy[k] = -g_z_0_yyyyz_yy[k] - g_yz_0_yyyyz_yy[k] * ab_y + g_yz_0_yyyyz_yyy[k];

                g_yz_0_yyyyyz_yz[k] = -g_z_0_yyyyz_yz[k] - g_yz_0_yyyyz_yz[k] * ab_y + g_yz_0_yyyyz_yyz[k];

                g_yz_0_yyyyyz_zz[k] = -g_z_0_yyyyz_zz[k] - g_yz_0_yyyyz_zz[k] * ab_y + g_yz_0_yyyyz_yzz[k];
            }

            /// Set up 810-816 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 810 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 811 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 812 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 813 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 814 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 815 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyzz_xx, g_yz_0_yyyyzz_xy, g_yz_0_yyyyzz_xz, g_yz_0_yyyyzz_yy, g_yz_0_yyyyzz_yz, g_yz_0_yyyyzz_zz, g_yz_0_yyyzz_xx, g_yz_0_yyyzz_xxy, g_yz_0_yyyzz_xy, g_yz_0_yyyzz_xyy, g_yz_0_yyyzz_xyz, g_yz_0_yyyzz_xz, g_yz_0_yyyzz_yy, g_yz_0_yyyzz_yyy, g_yz_0_yyyzz_yyz, g_yz_0_yyyzz_yz, g_yz_0_yyyzz_yzz, g_yz_0_yyyzz_zz, g_z_0_yyyzz_xx, g_z_0_yyyzz_xy, g_z_0_yyyzz_xz, g_z_0_yyyzz_yy, g_z_0_yyyzz_yz, g_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyzz_xx[k] = -g_z_0_yyyzz_xx[k] - g_yz_0_yyyzz_xx[k] * ab_y + g_yz_0_yyyzz_xxy[k];

                g_yz_0_yyyyzz_xy[k] = -g_z_0_yyyzz_xy[k] - g_yz_0_yyyzz_xy[k] * ab_y + g_yz_0_yyyzz_xyy[k];

                g_yz_0_yyyyzz_xz[k] = -g_z_0_yyyzz_xz[k] - g_yz_0_yyyzz_xz[k] * ab_y + g_yz_0_yyyzz_xyz[k];

                g_yz_0_yyyyzz_yy[k] = -g_z_0_yyyzz_yy[k] - g_yz_0_yyyzz_yy[k] * ab_y + g_yz_0_yyyzz_yyy[k];

                g_yz_0_yyyyzz_yz[k] = -g_z_0_yyyzz_yz[k] - g_yz_0_yyyzz_yz[k] * ab_y + g_yz_0_yyyzz_yyz[k];

                g_yz_0_yyyyzz_zz[k] = -g_z_0_yyyzz_zz[k] - g_yz_0_yyyzz_zz[k] * ab_y + g_yz_0_yyyzz_yzz[k];
            }

            /// Set up 816-822 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 816 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 817 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 818 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 819 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 820 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 821 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzzz_xx, g_yz_0_yyyzzz_xy, g_yz_0_yyyzzz_xz, g_yz_0_yyyzzz_yy, g_yz_0_yyyzzz_yz, g_yz_0_yyyzzz_zz, g_yz_0_yyzzz_xx, g_yz_0_yyzzz_xxy, g_yz_0_yyzzz_xy, g_yz_0_yyzzz_xyy, g_yz_0_yyzzz_xyz, g_yz_0_yyzzz_xz, g_yz_0_yyzzz_yy, g_yz_0_yyzzz_yyy, g_yz_0_yyzzz_yyz, g_yz_0_yyzzz_yz, g_yz_0_yyzzz_yzz, g_yz_0_yyzzz_zz, g_z_0_yyzzz_xx, g_z_0_yyzzz_xy, g_z_0_yyzzz_xz, g_z_0_yyzzz_yy, g_z_0_yyzzz_yz, g_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzzz_xx[k] = -g_z_0_yyzzz_xx[k] - g_yz_0_yyzzz_xx[k] * ab_y + g_yz_0_yyzzz_xxy[k];

                g_yz_0_yyyzzz_xy[k] = -g_z_0_yyzzz_xy[k] - g_yz_0_yyzzz_xy[k] * ab_y + g_yz_0_yyzzz_xyy[k];

                g_yz_0_yyyzzz_xz[k] = -g_z_0_yyzzz_xz[k] - g_yz_0_yyzzz_xz[k] * ab_y + g_yz_0_yyzzz_xyz[k];

                g_yz_0_yyyzzz_yy[k] = -g_z_0_yyzzz_yy[k] - g_yz_0_yyzzz_yy[k] * ab_y + g_yz_0_yyzzz_yyy[k];

                g_yz_0_yyyzzz_yz[k] = -g_z_0_yyzzz_yz[k] - g_yz_0_yyzzz_yz[k] * ab_y + g_yz_0_yyzzz_yyz[k];

                g_yz_0_yyyzzz_zz[k] = -g_z_0_yyzzz_zz[k] - g_yz_0_yyzzz_zz[k] * ab_y + g_yz_0_yyzzz_yzz[k];
            }

            /// Set up 822-828 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 822 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 823 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 824 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 825 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 826 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 827 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzzz_xx, g_yz_0_yyzzzz_xy, g_yz_0_yyzzzz_xz, g_yz_0_yyzzzz_yy, g_yz_0_yyzzzz_yz, g_yz_0_yyzzzz_zz, g_yz_0_yzzzz_xx, g_yz_0_yzzzz_xxy, g_yz_0_yzzzz_xy, g_yz_0_yzzzz_xyy, g_yz_0_yzzzz_xyz, g_yz_0_yzzzz_xz, g_yz_0_yzzzz_yy, g_yz_0_yzzzz_yyy, g_yz_0_yzzzz_yyz, g_yz_0_yzzzz_yz, g_yz_0_yzzzz_yzz, g_yz_0_yzzzz_zz, g_z_0_yzzzz_xx, g_z_0_yzzzz_xy, g_z_0_yzzzz_xz, g_z_0_yzzzz_yy, g_z_0_yzzzz_yz, g_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzzz_xx[k] = -g_z_0_yzzzz_xx[k] - g_yz_0_yzzzz_xx[k] * ab_y + g_yz_0_yzzzz_xxy[k];

                g_yz_0_yyzzzz_xy[k] = -g_z_0_yzzzz_xy[k] - g_yz_0_yzzzz_xy[k] * ab_y + g_yz_0_yzzzz_xyy[k];

                g_yz_0_yyzzzz_xz[k] = -g_z_0_yzzzz_xz[k] - g_yz_0_yzzzz_xz[k] * ab_y + g_yz_0_yzzzz_xyz[k];

                g_yz_0_yyzzzz_yy[k] = -g_z_0_yzzzz_yy[k] - g_yz_0_yzzzz_yy[k] * ab_y + g_yz_0_yzzzz_yyy[k];

                g_yz_0_yyzzzz_yz[k] = -g_z_0_yzzzz_yz[k] - g_yz_0_yzzzz_yz[k] * ab_y + g_yz_0_yzzzz_yyz[k];

                g_yz_0_yyzzzz_zz[k] = -g_z_0_yzzzz_zz[k] - g_yz_0_yzzzz_zz[k] * ab_y + g_yz_0_yzzzz_yzz[k];
            }

            /// Set up 828-834 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 828 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 829 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 830 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 831 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 832 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 833 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzzz_xx, g_yz_0_yzzzzz_xy, g_yz_0_yzzzzz_xz, g_yz_0_yzzzzz_yy, g_yz_0_yzzzzz_yz, g_yz_0_yzzzzz_zz, g_yz_0_zzzzz_xx, g_yz_0_zzzzz_xxy, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_xyy, g_yz_0_zzzzz_xyz, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_yy, g_yz_0_zzzzz_yyy, g_yz_0_zzzzz_yyz, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_yzz, g_yz_0_zzzzz_zz, g_z_0_zzzzz_xx, g_z_0_zzzzz_xy, g_z_0_zzzzz_xz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzzz_xx[k] = -g_z_0_zzzzz_xx[k] - g_yz_0_zzzzz_xx[k] * ab_y + g_yz_0_zzzzz_xxy[k];

                g_yz_0_yzzzzz_xy[k] = -g_z_0_zzzzz_xy[k] - g_yz_0_zzzzz_xy[k] * ab_y + g_yz_0_zzzzz_xyy[k];

                g_yz_0_yzzzzz_xz[k] = -g_z_0_zzzzz_xz[k] - g_yz_0_zzzzz_xz[k] * ab_y + g_yz_0_zzzzz_xyz[k];

                g_yz_0_yzzzzz_yy[k] = -g_z_0_zzzzz_yy[k] - g_yz_0_zzzzz_yy[k] * ab_y + g_yz_0_zzzzz_yyy[k];

                g_yz_0_yzzzzz_yz[k] = -g_z_0_zzzzz_yz[k] - g_yz_0_zzzzz_yz[k] * ab_y + g_yz_0_zzzzz_yyz[k];

                g_yz_0_yzzzzz_zz[k] = -g_z_0_zzzzz_zz[k] - g_yz_0_zzzzz_zz[k] * ab_y + g_yz_0_zzzzz_yzz[k];
            }

            /// Set up 834-840 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 834 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 835 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 836 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 837 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 838 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_xx, g_y_0_zzzzz_xy, g_y_0_zzzzz_xz, g_y_0_zzzzz_yy, g_y_0_zzzzz_yz, g_y_0_zzzzz_zz, g_yz_0_zzzzz_xx, g_yz_0_zzzzz_xxz, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_xyz, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_xzz, g_yz_0_zzzzz_yy, g_yz_0_zzzzz_yyz, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_yzz, g_yz_0_zzzzz_zz, g_yz_0_zzzzz_zzz, g_yz_0_zzzzzz_xx, g_yz_0_zzzzzz_xy, g_yz_0_zzzzzz_xz, g_yz_0_zzzzzz_yy, g_yz_0_zzzzzz_yz, g_yz_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzzz_xx[k] = -g_y_0_zzzzz_xx[k] - g_yz_0_zzzzz_xx[k] * ab_z + g_yz_0_zzzzz_xxz[k];

                g_yz_0_zzzzzz_xy[k] = -g_y_0_zzzzz_xy[k] - g_yz_0_zzzzz_xy[k] * ab_z + g_yz_0_zzzzz_xyz[k];

                g_yz_0_zzzzzz_xz[k] = -g_y_0_zzzzz_xz[k] - g_yz_0_zzzzz_xz[k] * ab_z + g_yz_0_zzzzz_xzz[k];

                g_yz_0_zzzzzz_yy[k] = -g_y_0_zzzzz_yy[k] - g_yz_0_zzzzz_yy[k] * ab_z + g_yz_0_zzzzz_yyz[k];

                g_yz_0_zzzzzz_yz[k] = -g_y_0_zzzzz_yz[k] - g_yz_0_zzzzz_yz[k] * ab_z + g_yz_0_zzzzz_yzz[k];

                g_yz_0_zzzzzz_zz[k] = -g_y_0_zzzzz_zz[k] - g_yz_0_zzzzz_zz[k] * ab_z + g_yz_0_zzzzz_zzz[k];
            }

            /// Set up 840-846 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxx_xx = cbuffer.data(id_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_xy = cbuffer.data(id_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_xz = cbuffer.data(id_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_yy = cbuffer.data(id_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_yz = cbuffer.data(id_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_zz = cbuffer.data(id_geom_20_off + 845 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxx_xx, g_zz_0_xxxxx_xxx, g_zz_0_xxxxx_xxy, g_zz_0_xxxxx_xxz, g_zz_0_xxxxx_xy, g_zz_0_xxxxx_xyy, g_zz_0_xxxxx_xyz, g_zz_0_xxxxx_xz, g_zz_0_xxxxx_xzz, g_zz_0_xxxxx_yy, g_zz_0_xxxxx_yz, g_zz_0_xxxxx_zz, g_zz_0_xxxxxx_xx, g_zz_0_xxxxxx_xy, g_zz_0_xxxxxx_xz, g_zz_0_xxxxxx_yy, g_zz_0_xxxxxx_yz, g_zz_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxx_xx[k] = -g_zz_0_xxxxx_xx[k] * ab_x + g_zz_0_xxxxx_xxx[k];

                g_zz_0_xxxxxx_xy[k] = -g_zz_0_xxxxx_xy[k] * ab_x + g_zz_0_xxxxx_xxy[k];

                g_zz_0_xxxxxx_xz[k] = -g_zz_0_xxxxx_xz[k] * ab_x + g_zz_0_xxxxx_xxz[k];

                g_zz_0_xxxxxx_yy[k] = -g_zz_0_xxxxx_yy[k] * ab_x + g_zz_0_xxxxx_xyy[k];

                g_zz_0_xxxxxx_yz[k] = -g_zz_0_xxxxx_yz[k] * ab_x + g_zz_0_xxxxx_xyz[k];

                g_zz_0_xxxxxx_zz[k] = -g_zz_0_xxxxx_zz[k] * ab_x + g_zz_0_xxxxx_xzz[k];
            }

            /// Set up 846-852 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxy_xx = cbuffer.data(id_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_xy = cbuffer.data(id_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_xz = cbuffer.data(id_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_yy = cbuffer.data(id_geom_20_off + 849 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_yz = cbuffer.data(id_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_zz = cbuffer.data(id_geom_20_off + 851 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxy_xx, g_zz_0_xxxxxy_xy, g_zz_0_xxxxxy_xz, g_zz_0_xxxxxy_yy, g_zz_0_xxxxxy_yz, g_zz_0_xxxxxy_zz, g_zz_0_xxxxy_xx, g_zz_0_xxxxy_xxx, g_zz_0_xxxxy_xxy, g_zz_0_xxxxy_xxz, g_zz_0_xxxxy_xy, g_zz_0_xxxxy_xyy, g_zz_0_xxxxy_xyz, g_zz_0_xxxxy_xz, g_zz_0_xxxxy_xzz, g_zz_0_xxxxy_yy, g_zz_0_xxxxy_yz, g_zz_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxy_xx[k] = -g_zz_0_xxxxy_xx[k] * ab_x + g_zz_0_xxxxy_xxx[k];

                g_zz_0_xxxxxy_xy[k] = -g_zz_0_xxxxy_xy[k] * ab_x + g_zz_0_xxxxy_xxy[k];

                g_zz_0_xxxxxy_xz[k] = -g_zz_0_xxxxy_xz[k] * ab_x + g_zz_0_xxxxy_xxz[k];

                g_zz_0_xxxxxy_yy[k] = -g_zz_0_xxxxy_yy[k] * ab_x + g_zz_0_xxxxy_xyy[k];

                g_zz_0_xxxxxy_yz[k] = -g_zz_0_xxxxy_yz[k] * ab_x + g_zz_0_xxxxy_xyz[k];

                g_zz_0_xxxxxy_zz[k] = -g_zz_0_xxxxy_zz[k] * ab_x + g_zz_0_xxxxy_xzz[k];
            }

            /// Set up 852-858 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxz_xx = cbuffer.data(id_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_xy = cbuffer.data(id_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_xz = cbuffer.data(id_geom_20_off + 854 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_yy = cbuffer.data(id_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_yz = cbuffer.data(id_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_zz = cbuffer.data(id_geom_20_off + 857 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxz_xx, g_zz_0_xxxxxz_xy, g_zz_0_xxxxxz_xz, g_zz_0_xxxxxz_yy, g_zz_0_xxxxxz_yz, g_zz_0_xxxxxz_zz, g_zz_0_xxxxz_xx, g_zz_0_xxxxz_xxx, g_zz_0_xxxxz_xxy, g_zz_0_xxxxz_xxz, g_zz_0_xxxxz_xy, g_zz_0_xxxxz_xyy, g_zz_0_xxxxz_xyz, g_zz_0_xxxxz_xz, g_zz_0_xxxxz_xzz, g_zz_0_xxxxz_yy, g_zz_0_xxxxz_yz, g_zz_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxz_xx[k] = -g_zz_0_xxxxz_xx[k] * ab_x + g_zz_0_xxxxz_xxx[k];

                g_zz_0_xxxxxz_xy[k] = -g_zz_0_xxxxz_xy[k] * ab_x + g_zz_0_xxxxz_xxy[k];

                g_zz_0_xxxxxz_xz[k] = -g_zz_0_xxxxz_xz[k] * ab_x + g_zz_0_xxxxz_xxz[k];

                g_zz_0_xxxxxz_yy[k] = -g_zz_0_xxxxz_yy[k] * ab_x + g_zz_0_xxxxz_xyy[k];

                g_zz_0_xxxxxz_yz[k] = -g_zz_0_xxxxz_yz[k] * ab_x + g_zz_0_xxxxz_xyz[k];

                g_zz_0_xxxxxz_zz[k] = -g_zz_0_xxxxz_zz[k] * ab_x + g_zz_0_xxxxz_xzz[k];
            }

            /// Set up 858-864 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyy_xx = cbuffer.data(id_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_xy = cbuffer.data(id_geom_20_off + 859 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_xz = cbuffer.data(id_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_yy = cbuffer.data(id_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_yz = cbuffer.data(id_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_zz = cbuffer.data(id_geom_20_off + 863 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyy_xx, g_zz_0_xxxxyy_xy, g_zz_0_xxxxyy_xz, g_zz_0_xxxxyy_yy, g_zz_0_xxxxyy_yz, g_zz_0_xxxxyy_zz, g_zz_0_xxxyy_xx, g_zz_0_xxxyy_xxx, g_zz_0_xxxyy_xxy, g_zz_0_xxxyy_xxz, g_zz_0_xxxyy_xy, g_zz_0_xxxyy_xyy, g_zz_0_xxxyy_xyz, g_zz_0_xxxyy_xz, g_zz_0_xxxyy_xzz, g_zz_0_xxxyy_yy, g_zz_0_xxxyy_yz, g_zz_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyy_xx[k] = -g_zz_0_xxxyy_xx[k] * ab_x + g_zz_0_xxxyy_xxx[k];

                g_zz_0_xxxxyy_xy[k] = -g_zz_0_xxxyy_xy[k] * ab_x + g_zz_0_xxxyy_xxy[k];

                g_zz_0_xxxxyy_xz[k] = -g_zz_0_xxxyy_xz[k] * ab_x + g_zz_0_xxxyy_xxz[k];

                g_zz_0_xxxxyy_yy[k] = -g_zz_0_xxxyy_yy[k] * ab_x + g_zz_0_xxxyy_xyy[k];

                g_zz_0_xxxxyy_yz[k] = -g_zz_0_xxxyy_yz[k] * ab_x + g_zz_0_xxxyy_xyz[k];

                g_zz_0_xxxxyy_zz[k] = -g_zz_0_xxxyy_zz[k] * ab_x + g_zz_0_xxxyy_xzz[k];
            }

            /// Set up 864-870 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyz_xx = cbuffer.data(id_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_xy = cbuffer.data(id_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_xz = cbuffer.data(id_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_yy = cbuffer.data(id_geom_20_off + 867 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_yz = cbuffer.data(id_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_zz = cbuffer.data(id_geom_20_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyz_xx, g_zz_0_xxxxyz_xy, g_zz_0_xxxxyz_xz, g_zz_0_xxxxyz_yy, g_zz_0_xxxxyz_yz, g_zz_0_xxxxyz_zz, g_zz_0_xxxyz_xx, g_zz_0_xxxyz_xxx, g_zz_0_xxxyz_xxy, g_zz_0_xxxyz_xxz, g_zz_0_xxxyz_xy, g_zz_0_xxxyz_xyy, g_zz_0_xxxyz_xyz, g_zz_0_xxxyz_xz, g_zz_0_xxxyz_xzz, g_zz_0_xxxyz_yy, g_zz_0_xxxyz_yz, g_zz_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyz_xx[k] = -g_zz_0_xxxyz_xx[k] * ab_x + g_zz_0_xxxyz_xxx[k];

                g_zz_0_xxxxyz_xy[k] = -g_zz_0_xxxyz_xy[k] * ab_x + g_zz_0_xxxyz_xxy[k];

                g_zz_0_xxxxyz_xz[k] = -g_zz_0_xxxyz_xz[k] * ab_x + g_zz_0_xxxyz_xxz[k];

                g_zz_0_xxxxyz_yy[k] = -g_zz_0_xxxyz_yy[k] * ab_x + g_zz_0_xxxyz_xyy[k];

                g_zz_0_xxxxyz_yz[k] = -g_zz_0_xxxyz_yz[k] * ab_x + g_zz_0_xxxyz_xyz[k];

                g_zz_0_xxxxyz_zz[k] = -g_zz_0_xxxyz_zz[k] * ab_x + g_zz_0_xxxyz_xzz[k];
            }

            /// Set up 870-876 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxzz_xx = cbuffer.data(id_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_xy = cbuffer.data(id_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_xz = cbuffer.data(id_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_yy = cbuffer.data(id_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_yz = cbuffer.data(id_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_zz = cbuffer.data(id_geom_20_off + 875 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxzz_xx, g_zz_0_xxxxzz_xy, g_zz_0_xxxxzz_xz, g_zz_0_xxxxzz_yy, g_zz_0_xxxxzz_yz, g_zz_0_xxxxzz_zz, g_zz_0_xxxzz_xx, g_zz_0_xxxzz_xxx, g_zz_0_xxxzz_xxy, g_zz_0_xxxzz_xxz, g_zz_0_xxxzz_xy, g_zz_0_xxxzz_xyy, g_zz_0_xxxzz_xyz, g_zz_0_xxxzz_xz, g_zz_0_xxxzz_xzz, g_zz_0_xxxzz_yy, g_zz_0_xxxzz_yz, g_zz_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxzz_xx[k] = -g_zz_0_xxxzz_xx[k] * ab_x + g_zz_0_xxxzz_xxx[k];

                g_zz_0_xxxxzz_xy[k] = -g_zz_0_xxxzz_xy[k] * ab_x + g_zz_0_xxxzz_xxy[k];

                g_zz_0_xxxxzz_xz[k] = -g_zz_0_xxxzz_xz[k] * ab_x + g_zz_0_xxxzz_xxz[k];

                g_zz_0_xxxxzz_yy[k] = -g_zz_0_xxxzz_yy[k] * ab_x + g_zz_0_xxxzz_xyy[k];

                g_zz_0_xxxxzz_yz[k] = -g_zz_0_xxxzz_yz[k] * ab_x + g_zz_0_xxxzz_xyz[k];

                g_zz_0_xxxxzz_zz[k] = -g_zz_0_xxxzz_zz[k] * ab_x + g_zz_0_xxxzz_xzz[k];
            }

            /// Set up 876-882 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyy_xx = cbuffer.data(id_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_xy = cbuffer.data(id_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_xz = cbuffer.data(id_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_yy = cbuffer.data(id_geom_20_off + 879 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_yz = cbuffer.data(id_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_zz = cbuffer.data(id_geom_20_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyy_xx, g_zz_0_xxxyyy_xy, g_zz_0_xxxyyy_xz, g_zz_0_xxxyyy_yy, g_zz_0_xxxyyy_yz, g_zz_0_xxxyyy_zz, g_zz_0_xxyyy_xx, g_zz_0_xxyyy_xxx, g_zz_0_xxyyy_xxy, g_zz_0_xxyyy_xxz, g_zz_0_xxyyy_xy, g_zz_0_xxyyy_xyy, g_zz_0_xxyyy_xyz, g_zz_0_xxyyy_xz, g_zz_0_xxyyy_xzz, g_zz_0_xxyyy_yy, g_zz_0_xxyyy_yz, g_zz_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyy_xx[k] = -g_zz_0_xxyyy_xx[k] * ab_x + g_zz_0_xxyyy_xxx[k];

                g_zz_0_xxxyyy_xy[k] = -g_zz_0_xxyyy_xy[k] * ab_x + g_zz_0_xxyyy_xxy[k];

                g_zz_0_xxxyyy_xz[k] = -g_zz_0_xxyyy_xz[k] * ab_x + g_zz_0_xxyyy_xxz[k];

                g_zz_0_xxxyyy_yy[k] = -g_zz_0_xxyyy_yy[k] * ab_x + g_zz_0_xxyyy_xyy[k];

                g_zz_0_xxxyyy_yz[k] = -g_zz_0_xxyyy_yz[k] * ab_x + g_zz_0_xxyyy_xyz[k];

                g_zz_0_xxxyyy_zz[k] = -g_zz_0_xxyyy_zz[k] * ab_x + g_zz_0_xxyyy_xzz[k];
            }

            /// Set up 882-888 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyz_xx = cbuffer.data(id_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_xy = cbuffer.data(id_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_xz = cbuffer.data(id_geom_20_off + 884 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_yy = cbuffer.data(id_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_yz = cbuffer.data(id_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_zz = cbuffer.data(id_geom_20_off + 887 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyz_xx, g_zz_0_xxxyyz_xy, g_zz_0_xxxyyz_xz, g_zz_0_xxxyyz_yy, g_zz_0_xxxyyz_yz, g_zz_0_xxxyyz_zz, g_zz_0_xxyyz_xx, g_zz_0_xxyyz_xxx, g_zz_0_xxyyz_xxy, g_zz_0_xxyyz_xxz, g_zz_0_xxyyz_xy, g_zz_0_xxyyz_xyy, g_zz_0_xxyyz_xyz, g_zz_0_xxyyz_xz, g_zz_0_xxyyz_xzz, g_zz_0_xxyyz_yy, g_zz_0_xxyyz_yz, g_zz_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyz_xx[k] = -g_zz_0_xxyyz_xx[k] * ab_x + g_zz_0_xxyyz_xxx[k];

                g_zz_0_xxxyyz_xy[k] = -g_zz_0_xxyyz_xy[k] * ab_x + g_zz_0_xxyyz_xxy[k];

                g_zz_0_xxxyyz_xz[k] = -g_zz_0_xxyyz_xz[k] * ab_x + g_zz_0_xxyyz_xxz[k];

                g_zz_0_xxxyyz_yy[k] = -g_zz_0_xxyyz_yy[k] * ab_x + g_zz_0_xxyyz_xyy[k];

                g_zz_0_xxxyyz_yz[k] = -g_zz_0_xxyyz_yz[k] * ab_x + g_zz_0_xxyyz_xyz[k];

                g_zz_0_xxxyyz_zz[k] = -g_zz_0_xxyyz_zz[k] * ab_x + g_zz_0_xxyyz_xzz[k];
            }

            /// Set up 888-894 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyzz_xx = cbuffer.data(id_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_xy = cbuffer.data(id_geom_20_off + 889 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_xz = cbuffer.data(id_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_yy = cbuffer.data(id_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_yz = cbuffer.data(id_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_zz = cbuffer.data(id_geom_20_off + 893 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyzz_xx, g_zz_0_xxxyzz_xy, g_zz_0_xxxyzz_xz, g_zz_0_xxxyzz_yy, g_zz_0_xxxyzz_yz, g_zz_0_xxxyzz_zz, g_zz_0_xxyzz_xx, g_zz_0_xxyzz_xxx, g_zz_0_xxyzz_xxy, g_zz_0_xxyzz_xxz, g_zz_0_xxyzz_xy, g_zz_0_xxyzz_xyy, g_zz_0_xxyzz_xyz, g_zz_0_xxyzz_xz, g_zz_0_xxyzz_xzz, g_zz_0_xxyzz_yy, g_zz_0_xxyzz_yz, g_zz_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyzz_xx[k] = -g_zz_0_xxyzz_xx[k] * ab_x + g_zz_0_xxyzz_xxx[k];

                g_zz_0_xxxyzz_xy[k] = -g_zz_0_xxyzz_xy[k] * ab_x + g_zz_0_xxyzz_xxy[k];

                g_zz_0_xxxyzz_xz[k] = -g_zz_0_xxyzz_xz[k] * ab_x + g_zz_0_xxyzz_xxz[k];

                g_zz_0_xxxyzz_yy[k] = -g_zz_0_xxyzz_yy[k] * ab_x + g_zz_0_xxyzz_xyy[k];

                g_zz_0_xxxyzz_yz[k] = -g_zz_0_xxyzz_yz[k] * ab_x + g_zz_0_xxyzz_xyz[k];

                g_zz_0_xxxyzz_zz[k] = -g_zz_0_xxyzz_zz[k] * ab_x + g_zz_0_xxyzz_xzz[k];
            }

            /// Set up 894-900 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzzz_xx = cbuffer.data(id_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_xy = cbuffer.data(id_geom_20_off + 895 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_xz = cbuffer.data(id_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_yy = cbuffer.data(id_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_yz = cbuffer.data(id_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_zz = cbuffer.data(id_geom_20_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzzz_xx, g_zz_0_xxxzzz_xy, g_zz_0_xxxzzz_xz, g_zz_0_xxxzzz_yy, g_zz_0_xxxzzz_yz, g_zz_0_xxxzzz_zz, g_zz_0_xxzzz_xx, g_zz_0_xxzzz_xxx, g_zz_0_xxzzz_xxy, g_zz_0_xxzzz_xxz, g_zz_0_xxzzz_xy, g_zz_0_xxzzz_xyy, g_zz_0_xxzzz_xyz, g_zz_0_xxzzz_xz, g_zz_0_xxzzz_xzz, g_zz_0_xxzzz_yy, g_zz_0_xxzzz_yz, g_zz_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzzz_xx[k] = -g_zz_0_xxzzz_xx[k] * ab_x + g_zz_0_xxzzz_xxx[k];

                g_zz_0_xxxzzz_xy[k] = -g_zz_0_xxzzz_xy[k] * ab_x + g_zz_0_xxzzz_xxy[k];

                g_zz_0_xxxzzz_xz[k] = -g_zz_0_xxzzz_xz[k] * ab_x + g_zz_0_xxzzz_xxz[k];

                g_zz_0_xxxzzz_yy[k] = -g_zz_0_xxzzz_yy[k] * ab_x + g_zz_0_xxzzz_xyy[k];

                g_zz_0_xxxzzz_yz[k] = -g_zz_0_xxzzz_yz[k] * ab_x + g_zz_0_xxzzz_xyz[k];

                g_zz_0_xxxzzz_zz[k] = -g_zz_0_xxzzz_zz[k] * ab_x + g_zz_0_xxzzz_xzz[k];
            }

            /// Set up 900-906 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyy_xx = cbuffer.data(id_geom_20_off + 900 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_xy = cbuffer.data(id_geom_20_off + 901 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_xz = cbuffer.data(id_geom_20_off + 902 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_yy = cbuffer.data(id_geom_20_off + 903 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_yz = cbuffer.data(id_geom_20_off + 904 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_zz = cbuffer.data(id_geom_20_off + 905 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyy_xx, g_zz_0_xxyyyy_xy, g_zz_0_xxyyyy_xz, g_zz_0_xxyyyy_yy, g_zz_0_xxyyyy_yz, g_zz_0_xxyyyy_zz, g_zz_0_xyyyy_xx, g_zz_0_xyyyy_xxx, g_zz_0_xyyyy_xxy, g_zz_0_xyyyy_xxz, g_zz_0_xyyyy_xy, g_zz_0_xyyyy_xyy, g_zz_0_xyyyy_xyz, g_zz_0_xyyyy_xz, g_zz_0_xyyyy_xzz, g_zz_0_xyyyy_yy, g_zz_0_xyyyy_yz, g_zz_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyy_xx[k] = -g_zz_0_xyyyy_xx[k] * ab_x + g_zz_0_xyyyy_xxx[k];

                g_zz_0_xxyyyy_xy[k] = -g_zz_0_xyyyy_xy[k] * ab_x + g_zz_0_xyyyy_xxy[k];

                g_zz_0_xxyyyy_xz[k] = -g_zz_0_xyyyy_xz[k] * ab_x + g_zz_0_xyyyy_xxz[k];

                g_zz_0_xxyyyy_yy[k] = -g_zz_0_xyyyy_yy[k] * ab_x + g_zz_0_xyyyy_xyy[k];

                g_zz_0_xxyyyy_yz[k] = -g_zz_0_xyyyy_yz[k] * ab_x + g_zz_0_xyyyy_xyz[k];

                g_zz_0_xxyyyy_zz[k] = -g_zz_0_xyyyy_zz[k] * ab_x + g_zz_0_xyyyy_xzz[k];
            }

            /// Set up 906-912 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyz_xx = cbuffer.data(id_geom_20_off + 906 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_xy = cbuffer.data(id_geom_20_off + 907 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_xz = cbuffer.data(id_geom_20_off + 908 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_yy = cbuffer.data(id_geom_20_off + 909 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_yz = cbuffer.data(id_geom_20_off + 910 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_zz = cbuffer.data(id_geom_20_off + 911 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyz_xx, g_zz_0_xxyyyz_xy, g_zz_0_xxyyyz_xz, g_zz_0_xxyyyz_yy, g_zz_0_xxyyyz_yz, g_zz_0_xxyyyz_zz, g_zz_0_xyyyz_xx, g_zz_0_xyyyz_xxx, g_zz_0_xyyyz_xxy, g_zz_0_xyyyz_xxz, g_zz_0_xyyyz_xy, g_zz_0_xyyyz_xyy, g_zz_0_xyyyz_xyz, g_zz_0_xyyyz_xz, g_zz_0_xyyyz_xzz, g_zz_0_xyyyz_yy, g_zz_0_xyyyz_yz, g_zz_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyz_xx[k] = -g_zz_0_xyyyz_xx[k] * ab_x + g_zz_0_xyyyz_xxx[k];

                g_zz_0_xxyyyz_xy[k] = -g_zz_0_xyyyz_xy[k] * ab_x + g_zz_0_xyyyz_xxy[k];

                g_zz_0_xxyyyz_xz[k] = -g_zz_0_xyyyz_xz[k] * ab_x + g_zz_0_xyyyz_xxz[k];

                g_zz_0_xxyyyz_yy[k] = -g_zz_0_xyyyz_yy[k] * ab_x + g_zz_0_xyyyz_xyy[k];

                g_zz_0_xxyyyz_yz[k] = -g_zz_0_xyyyz_yz[k] * ab_x + g_zz_0_xyyyz_xyz[k];

                g_zz_0_xxyyyz_zz[k] = -g_zz_0_xyyyz_zz[k] * ab_x + g_zz_0_xyyyz_xzz[k];
            }

            /// Set up 912-918 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyzz_xx = cbuffer.data(id_geom_20_off + 912 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_xy = cbuffer.data(id_geom_20_off + 913 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_xz = cbuffer.data(id_geom_20_off + 914 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_yy = cbuffer.data(id_geom_20_off + 915 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_yz = cbuffer.data(id_geom_20_off + 916 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_zz = cbuffer.data(id_geom_20_off + 917 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyzz_xx, g_zz_0_xxyyzz_xy, g_zz_0_xxyyzz_xz, g_zz_0_xxyyzz_yy, g_zz_0_xxyyzz_yz, g_zz_0_xxyyzz_zz, g_zz_0_xyyzz_xx, g_zz_0_xyyzz_xxx, g_zz_0_xyyzz_xxy, g_zz_0_xyyzz_xxz, g_zz_0_xyyzz_xy, g_zz_0_xyyzz_xyy, g_zz_0_xyyzz_xyz, g_zz_0_xyyzz_xz, g_zz_0_xyyzz_xzz, g_zz_0_xyyzz_yy, g_zz_0_xyyzz_yz, g_zz_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyzz_xx[k] = -g_zz_0_xyyzz_xx[k] * ab_x + g_zz_0_xyyzz_xxx[k];

                g_zz_0_xxyyzz_xy[k] = -g_zz_0_xyyzz_xy[k] * ab_x + g_zz_0_xyyzz_xxy[k];

                g_zz_0_xxyyzz_xz[k] = -g_zz_0_xyyzz_xz[k] * ab_x + g_zz_0_xyyzz_xxz[k];

                g_zz_0_xxyyzz_yy[k] = -g_zz_0_xyyzz_yy[k] * ab_x + g_zz_0_xyyzz_xyy[k];

                g_zz_0_xxyyzz_yz[k] = -g_zz_0_xyyzz_yz[k] * ab_x + g_zz_0_xyyzz_xyz[k];

                g_zz_0_xxyyzz_zz[k] = -g_zz_0_xyyzz_zz[k] * ab_x + g_zz_0_xyyzz_xzz[k];
            }

            /// Set up 918-924 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzzz_xx = cbuffer.data(id_geom_20_off + 918 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_xy = cbuffer.data(id_geom_20_off + 919 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_xz = cbuffer.data(id_geom_20_off + 920 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_yy = cbuffer.data(id_geom_20_off + 921 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_yz = cbuffer.data(id_geom_20_off + 922 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_zz = cbuffer.data(id_geom_20_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzzz_xx, g_zz_0_xxyzzz_xy, g_zz_0_xxyzzz_xz, g_zz_0_xxyzzz_yy, g_zz_0_xxyzzz_yz, g_zz_0_xxyzzz_zz, g_zz_0_xyzzz_xx, g_zz_0_xyzzz_xxx, g_zz_0_xyzzz_xxy, g_zz_0_xyzzz_xxz, g_zz_0_xyzzz_xy, g_zz_0_xyzzz_xyy, g_zz_0_xyzzz_xyz, g_zz_0_xyzzz_xz, g_zz_0_xyzzz_xzz, g_zz_0_xyzzz_yy, g_zz_0_xyzzz_yz, g_zz_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzzz_xx[k] = -g_zz_0_xyzzz_xx[k] * ab_x + g_zz_0_xyzzz_xxx[k];

                g_zz_0_xxyzzz_xy[k] = -g_zz_0_xyzzz_xy[k] * ab_x + g_zz_0_xyzzz_xxy[k];

                g_zz_0_xxyzzz_xz[k] = -g_zz_0_xyzzz_xz[k] * ab_x + g_zz_0_xyzzz_xxz[k];

                g_zz_0_xxyzzz_yy[k] = -g_zz_0_xyzzz_yy[k] * ab_x + g_zz_0_xyzzz_xyy[k];

                g_zz_0_xxyzzz_yz[k] = -g_zz_0_xyzzz_yz[k] * ab_x + g_zz_0_xyzzz_xyz[k];

                g_zz_0_xxyzzz_zz[k] = -g_zz_0_xyzzz_zz[k] * ab_x + g_zz_0_xyzzz_xzz[k];
            }

            /// Set up 924-930 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzzz_xx = cbuffer.data(id_geom_20_off + 924 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_xy = cbuffer.data(id_geom_20_off + 925 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_xz = cbuffer.data(id_geom_20_off + 926 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_yy = cbuffer.data(id_geom_20_off + 927 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_yz = cbuffer.data(id_geom_20_off + 928 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_zz = cbuffer.data(id_geom_20_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzzz_xx, g_zz_0_xxzzzz_xy, g_zz_0_xxzzzz_xz, g_zz_0_xxzzzz_yy, g_zz_0_xxzzzz_yz, g_zz_0_xxzzzz_zz, g_zz_0_xzzzz_xx, g_zz_0_xzzzz_xxx, g_zz_0_xzzzz_xxy, g_zz_0_xzzzz_xxz, g_zz_0_xzzzz_xy, g_zz_0_xzzzz_xyy, g_zz_0_xzzzz_xyz, g_zz_0_xzzzz_xz, g_zz_0_xzzzz_xzz, g_zz_0_xzzzz_yy, g_zz_0_xzzzz_yz, g_zz_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzzz_xx[k] = -g_zz_0_xzzzz_xx[k] * ab_x + g_zz_0_xzzzz_xxx[k];

                g_zz_0_xxzzzz_xy[k] = -g_zz_0_xzzzz_xy[k] * ab_x + g_zz_0_xzzzz_xxy[k];

                g_zz_0_xxzzzz_xz[k] = -g_zz_0_xzzzz_xz[k] * ab_x + g_zz_0_xzzzz_xxz[k];

                g_zz_0_xxzzzz_yy[k] = -g_zz_0_xzzzz_yy[k] * ab_x + g_zz_0_xzzzz_xyy[k];

                g_zz_0_xxzzzz_yz[k] = -g_zz_0_xzzzz_yz[k] * ab_x + g_zz_0_xzzzz_xyz[k];

                g_zz_0_xxzzzz_zz[k] = -g_zz_0_xzzzz_zz[k] * ab_x + g_zz_0_xzzzz_xzz[k];
            }

            /// Set up 930-936 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyy_xx = cbuffer.data(id_geom_20_off + 930 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_xy = cbuffer.data(id_geom_20_off + 931 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_xz = cbuffer.data(id_geom_20_off + 932 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_yy = cbuffer.data(id_geom_20_off + 933 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_yz = cbuffer.data(id_geom_20_off + 934 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_zz = cbuffer.data(id_geom_20_off + 935 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyy_xx, g_zz_0_xyyyyy_xy, g_zz_0_xyyyyy_xz, g_zz_0_xyyyyy_yy, g_zz_0_xyyyyy_yz, g_zz_0_xyyyyy_zz, g_zz_0_yyyyy_xx, g_zz_0_yyyyy_xxx, g_zz_0_yyyyy_xxy, g_zz_0_yyyyy_xxz, g_zz_0_yyyyy_xy, g_zz_0_yyyyy_xyy, g_zz_0_yyyyy_xyz, g_zz_0_yyyyy_xz, g_zz_0_yyyyy_xzz, g_zz_0_yyyyy_yy, g_zz_0_yyyyy_yz, g_zz_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyy_xx[k] = -g_zz_0_yyyyy_xx[k] * ab_x + g_zz_0_yyyyy_xxx[k];

                g_zz_0_xyyyyy_xy[k] = -g_zz_0_yyyyy_xy[k] * ab_x + g_zz_0_yyyyy_xxy[k];

                g_zz_0_xyyyyy_xz[k] = -g_zz_0_yyyyy_xz[k] * ab_x + g_zz_0_yyyyy_xxz[k];

                g_zz_0_xyyyyy_yy[k] = -g_zz_0_yyyyy_yy[k] * ab_x + g_zz_0_yyyyy_xyy[k];

                g_zz_0_xyyyyy_yz[k] = -g_zz_0_yyyyy_yz[k] * ab_x + g_zz_0_yyyyy_xyz[k];

                g_zz_0_xyyyyy_zz[k] = -g_zz_0_yyyyy_zz[k] * ab_x + g_zz_0_yyyyy_xzz[k];
            }

            /// Set up 936-942 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyz_xx = cbuffer.data(id_geom_20_off + 936 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_xy = cbuffer.data(id_geom_20_off + 937 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_xz = cbuffer.data(id_geom_20_off + 938 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_yy = cbuffer.data(id_geom_20_off + 939 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_yz = cbuffer.data(id_geom_20_off + 940 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_zz = cbuffer.data(id_geom_20_off + 941 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyz_xx, g_zz_0_xyyyyz_xy, g_zz_0_xyyyyz_xz, g_zz_0_xyyyyz_yy, g_zz_0_xyyyyz_yz, g_zz_0_xyyyyz_zz, g_zz_0_yyyyz_xx, g_zz_0_yyyyz_xxx, g_zz_0_yyyyz_xxy, g_zz_0_yyyyz_xxz, g_zz_0_yyyyz_xy, g_zz_0_yyyyz_xyy, g_zz_0_yyyyz_xyz, g_zz_0_yyyyz_xz, g_zz_0_yyyyz_xzz, g_zz_0_yyyyz_yy, g_zz_0_yyyyz_yz, g_zz_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyz_xx[k] = -g_zz_0_yyyyz_xx[k] * ab_x + g_zz_0_yyyyz_xxx[k];

                g_zz_0_xyyyyz_xy[k] = -g_zz_0_yyyyz_xy[k] * ab_x + g_zz_0_yyyyz_xxy[k];

                g_zz_0_xyyyyz_xz[k] = -g_zz_0_yyyyz_xz[k] * ab_x + g_zz_0_yyyyz_xxz[k];

                g_zz_0_xyyyyz_yy[k] = -g_zz_0_yyyyz_yy[k] * ab_x + g_zz_0_yyyyz_xyy[k];

                g_zz_0_xyyyyz_yz[k] = -g_zz_0_yyyyz_yz[k] * ab_x + g_zz_0_yyyyz_xyz[k];

                g_zz_0_xyyyyz_zz[k] = -g_zz_0_yyyyz_zz[k] * ab_x + g_zz_0_yyyyz_xzz[k];
            }

            /// Set up 942-948 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyzz_xx = cbuffer.data(id_geom_20_off + 942 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_xy = cbuffer.data(id_geom_20_off + 943 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_xz = cbuffer.data(id_geom_20_off + 944 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_yy = cbuffer.data(id_geom_20_off + 945 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_yz = cbuffer.data(id_geom_20_off + 946 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_zz = cbuffer.data(id_geom_20_off + 947 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyzz_xx, g_zz_0_xyyyzz_xy, g_zz_0_xyyyzz_xz, g_zz_0_xyyyzz_yy, g_zz_0_xyyyzz_yz, g_zz_0_xyyyzz_zz, g_zz_0_yyyzz_xx, g_zz_0_yyyzz_xxx, g_zz_0_yyyzz_xxy, g_zz_0_yyyzz_xxz, g_zz_0_yyyzz_xy, g_zz_0_yyyzz_xyy, g_zz_0_yyyzz_xyz, g_zz_0_yyyzz_xz, g_zz_0_yyyzz_xzz, g_zz_0_yyyzz_yy, g_zz_0_yyyzz_yz, g_zz_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyzz_xx[k] = -g_zz_0_yyyzz_xx[k] * ab_x + g_zz_0_yyyzz_xxx[k];

                g_zz_0_xyyyzz_xy[k] = -g_zz_0_yyyzz_xy[k] * ab_x + g_zz_0_yyyzz_xxy[k];

                g_zz_0_xyyyzz_xz[k] = -g_zz_0_yyyzz_xz[k] * ab_x + g_zz_0_yyyzz_xxz[k];

                g_zz_0_xyyyzz_yy[k] = -g_zz_0_yyyzz_yy[k] * ab_x + g_zz_0_yyyzz_xyy[k];

                g_zz_0_xyyyzz_yz[k] = -g_zz_0_yyyzz_yz[k] * ab_x + g_zz_0_yyyzz_xyz[k];

                g_zz_0_xyyyzz_zz[k] = -g_zz_0_yyyzz_zz[k] * ab_x + g_zz_0_yyyzz_xzz[k];
            }

            /// Set up 948-954 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzzz_xx = cbuffer.data(id_geom_20_off + 948 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_xy = cbuffer.data(id_geom_20_off + 949 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_xz = cbuffer.data(id_geom_20_off + 950 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_yy = cbuffer.data(id_geom_20_off + 951 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_yz = cbuffer.data(id_geom_20_off + 952 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_zz = cbuffer.data(id_geom_20_off + 953 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzzz_xx, g_zz_0_xyyzzz_xy, g_zz_0_xyyzzz_xz, g_zz_0_xyyzzz_yy, g_zz_0_xyyzzz_yz, g_zz_0_xyyzzz_zz, g_zz_0_yyzzz_xx, g_zz_0_yyzzz_xxx, g_zz_0_yyzzz_xxy, g_zz_0_yyzzz_xxz, g_zz_0_yyzzz_xy, g_zz_0_yyzzz_xyy, g_zz_0_yyzzz_xyz, g_zz_0_yyzzz_xz, g_zz_0_yyzzz_xzz, g_zz_0_yyzzz_yy, g_zz_0_yyzzz_yz, g_zz_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzzz_xx[k] = -g_zz_0_yyzzz_xx[k] * ab_x + g_zz_0_yyzzz_xxx[k];

                g_zz_0_xyyzzz_xy[k] = -g_zz_0_yyzzz_xy[k] * ab_x + g_zz_0_yyzzz_xxy[k];

                g_zz_0_xyyzzz_xz[k] = -g_zz_0_yyzzz_xz[k] * ab_x + g_zz_0_yyzzz_xxz[k];

                g_zz_0_xyyzzz_yy[k] = -g_zz_0_yyzzz_yy[k] * ab_x + g_zz_0_yyzzz_xyy[k];

                g_zz_0_xyyzzz_yz[k] = -g_zz_0_yyzzz_yz[k] * ab_x + g_zz_0_yyzzz_xyz[k];

                g_zz_0_xyyzzz_zz[k] = -g_zz_0_yyzzz_zz[k] * ab_x + g_zz_0_yyzzz_xzz[k];
            }

            /// Set up 954-960 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzzz_xx = cbuffer.data(id_geom_20_off + 954 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_xy = cbuffer.data(id_geom_20_off + 955 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_xz = cbuffer.data(id_geom_20_off + 956 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_yy = cbuffer.data(id_geom_20_off + 957 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_yz = cbuffer.data(id_geom_20_off + 958 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_zz = cbuffer.data(id_geom_20_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzzz_xx, g_zz_0_xyzzzz_xy, g_zz_0_xyzzzz_xz, g_zz_0_xyzzzz_yy, g_zz_0_xyzzzz_yz, g_zz_0_xyzzzz_zz, g_zz_0_yzzzz_xx, g_zz_0_yzzzz_xxx, g_zz_0_yzzzz_xxy, g_zz_0_yzzzz_xxz, g_zz_0_yzzzz_xy, g_zz_0_yzzzz_xyy, g_zz_0_yzzzz_xyz, g_zz_0_yzzzz_xz, g_zz_0_yzzzz_xzz, g_zz_0_yzzzz_yy, g_zz_0_yzzzz_yz, g_zz_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzzz_xx[k] = -g_zz_0_yzzzz_xx[k] * ab_x + g_zz_0_yzzzz_xxx[k];

                g_zz_0_xyzzzz_xy[k] = -g_zz_0_yzzzz_xy[k] * ab_x + g_zz_0_yzzzz_xxy[k];

                g_zz_0_xyzzzz_xz[k] = -g_zz_0_yzzzz_xz[k] * ab_x + g_zz_0_yzzzz_xxz[k];

                g_zz_0_xyzzzz_yy[k] = -g_zz_0_yzzzz_yy[k] * ab_x + g_zz_0_yzzzz_xyy[k];

                g_zz_0_xyzzzz_yz[k] = -g_zz_0_yzzzz_yz[k] * ab_x + g_zz_0_yzzzz_xyz[k];

                g_zz_0_xyzzzz_zz[k] = -g_zz_0_yzzzz_zz[k] * ab_x + g_zz_0_yzzzz_xzz[k];
            }

            /// Set up 960-966 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzzz_xx = cbuffer.data(id_geom_20_off + 960 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_xy = cbuffer.data(id_geom_20_off + 961 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_xz = cbuffer.data(id_geom_20_off + 962 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_yy = cbuffer.data(id_geom_20_off + 963 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_yz = cbuffer.data(id_geom_20_off + 964 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_zz = cbuffer.data(id_geom_20_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzzz_xx, g_zz_0_xzzzzz_xy, g_zz_0_xzzzzz_xz, g_zz_0_xzzzzz_yy, g_zz_0_xzzzzz_yz, g_zz_0_xzzzzz_zz, g_zz_0_zzzzz_xx, g_zz_0_zzzzz_xxx, g_zz_0_zzzzz_xxy, g_zz_0_zzzzz_xxz, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_xyy, g_zz_0_zzzzz_xyz, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_xzz, g_zz_0_zzzzz_yy, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzzz_xx[k] = -g_zz_0_zzzzz_xx[k] * ab_x + g_zz_0_zzzzz_xxx[k];

                g_zz_0_xzzzzz_xy[k] = -g_zz_0_zzzzz_xy[k] * ab_x + g_zz_0_zzzzz_xxy[k];

                g_zz_0_xzzzzz_xz[k] = -g_zz_0_zzzzz_xz[k] * ab_x + g_zz_0_zzzzz_xxz[k];

                g_zz_0_xzzzzz_yy[k] = -g_zz_0_zzzzz_yy[k] * ab_x + g_zz_0_zzzzz_xyy[k];

                g_zz_0_xzzzzz_yz[k] = -g_zz_0_zzzzz_yz[k] * ab_x + g_zz_0_zzzzz_xyz[k];

                g_zz_0_xzzzzz_zz[k] = -g_zz_0_zzzzz_zz[k] * ab_x + g_zz_0_zzzzz_xzz[k];
            }

            /// Set up 966-972 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyy_xx = cbuffer.data(id_geom_20_off + 966 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_xy = cbuffer.data(id_geom_20_off + 967 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_xz = cbuffer.data(id_geom_20_off + 968 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_yy = cbuffer.data(id_geom_20_off + 969 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_yz = cbuffer.data(id_geom_20_off + 970 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_zz = cbuffer.data(id_geom_20_off + 971 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyy_xx, g_zz_0_yyyyy_xxy, g_zz_0_yyyyy_xy, g_zz_0_yyyyy_xyy, g_zz_0_yyyyy_xyz, g_zz_0_yyyyy_xz, g_zz_0_yyyyy_yy, g_zz_0_yyyyy_yyy, g_zz_0_yyyyy_yyz, g_zz_0_yyyyy_yz, g_zz_0_yyyyy_yzz, g_zz_0_yyyyy_zz, g_zz_0_yyyyyy_xx, g_zz_0_yyyyyy_xy, g_zz_0_yyyyyy_xz, g_zz_0_yyyyyy_yy, g_zz_0_yyyyyy_yz, g_zz_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyy_xx[k] = -g_zz_0_yyyyy_xx[k] * ab_y + g_zz_0_yyyyy_xxy[k];

                g_zz_0_yyyyyy_xy[k] = -g_zz_0_yyyyy_xy[k] * ab_y + g_zz_0_yyyyy_xyy[k];

                g_zz_0_yyyyyy_xz[k] = -g_zz_0_yyyyy_xz[k] * ab_y + g_zz_0_yyyyy_xyz[k];

                g_zz_0_yyyyyy_yy[k] = -g_zz_0_yyyyy_yy[k] * ab_y + g_zz_0_yyyyy_yyy[k];

                g_zz_0_yyyyyy_yz[k] = -g_zz_0_yyyyy_yz[k] * ab_y + g_zz_0_yyyyy_yyz[k];

                g_zz_0_yyyyyy_zz[k] = -g_zz_0_yyyyy_zz[k] * ab_y + g_zz_0_yyyyy_yzz[k];
            }

            /// Set up 972-978 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyz_xx = cbuffer.data(id_geom_20_off + 972 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_xy = cbuffer.data(id_geom_20_off + 973 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_xz = cbuffer.data(id_geom_20_off + 974 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_yy = cbuffer.data(id_geom_20_off + 975 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_yz = cbuffer.data(id_geom_20_off + 976 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_zz = cbuffer.data(id_geom_20_off + 977 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyyz_xx, g_zz_0_yyyyyz_xy, g_zz_0_yyyyyz_xz, g_zz_0_yyyyyz_yy, g_zz_0_yyyyyz_yz, g_zz_0_yyyyyz_zz, g_zz_0_yyyyz_xx, g_zz_0_yyyyz_xxy, g_zz_0_yyyyz_xy, g_zz_0_yyyyz_xyy, g_zz_0_yyyyz_xyz, g_zz_0_yyyyz_xz, g_zz_0_yyyyz_yy, g_zz_0_yyyyz_yyy, g_zz_0_yyyyz_yyz, g_zz_0_yyyyz_yz, g_zz_0_yyyyz_yzz, g_zz_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyz_xx[k] = -g_zz_0_yyyyz_xx[k] * ab_y + g_zz_0_yyyyz_xxy[k];

                g_zz_0_yyyyyz_xy[k] = -g_zz_0_yyyyz_xy[k] * ab_y + g_zz_0_yyyyz_xyy[k];

                g_zz_0_yyyyyz_xz[k] = -g_zz_0_yyyyz_xz[k] * ab_y + g_zz_0_yyyyz_xyz[k];

                g_zz_0_yyyyyz_yy[k] = -g_zz_0_yyyyz_yy[k] * ab_y + g_zz_0_yyyyz_yyy[k];

                g_zz_0_yyyyyz_yz[k] = -g_zz_0_yyyyz_yz[k] * ab_y + g_zz_0_yyyyz_yyz[k];

                g_zz_0_yyyyyz_zz[k] = -g_zz_0_yyyyz_zz[k] * ab_y + g_zz_0_yyyyz_yzz[k];
            }

            /// Set up 978-984 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyzz_xx = cbuffer.data(id_geom_20_off + 978 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_xy = cbuffer.data(id_geom_20_off + 979 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_xz = cbuffer.data(id_geom_20_off + 980 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_yy = cbuffer.data(id_geom_20_off + 981 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_yz = cbuffer.data(id_geom_20_off + 982 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_zz = cbuffer.data(id_geom_20_off + 983 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyzz_xx, g_zz_0_yyyyzz_xy, g_zz_0_yyyyzz_xz, g_zz_0_yyyyzz_yy, g_zz_0_yyyyzz_yz, g_zz_0_yyyyzz_zz, g_zz_0_yyyzz_xx, g_zz_0_yyyzz_xxy, g_zz_0_yyyzz_xy, g_zz_0_yyyzz_xyy, g_zz_0_yyyzz_xyz, g_zz_0_yyyzz_xz, g_zz_0_yyyzz_yy, g_zz_0_yyyzz_yyy, g_zz_0_yyyzz_yyz, g_zz_0_yyyzz_yz, g_zz_0_yyyzz_yzz, g_zz_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyzz_xx[k] = -g_zz_0_yyyzz_xx[k] * ab_y + g_zz_0_yyyzz_xxy[k];

                g_zz_0_yyyyzz_xy[k] = -g_zz_0_yyyzz_xy[k] * ab_y + g_zz_0_yyyzz_xyy[k];

                g_zz_0_yyyyzz_xz[k] = -g_zz_0_yyyzz_xz[k] * ab_y + g_zz_0_yyyzz_xyz[k];

                g_zz_0_yyyyzz_yy[k] = -g_zz_0_yyyzz_yy[k] * ab_y + g_zz_0_yyyzz_yyy[k];

                g_zz_0_yyyyzz_yz[k] = -g_zz_0_yyyzz_yz[k] * ab_y + g_zz_0_yyyzz_yyz[k];

                g_zz_0_yyyyzz_zz[k] = -g_zz_0_yyyzz_zz[k] * ab_y + g_zz_0_yyyzz_yzz[k];
            }

            /// Set up 984-990 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzzz_xx = cbuffer.data(id_geom_20_off + 984 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_xy = cbuffer.data(id_geom_20_off + 985 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_xz = cbuffer.data(id_geom_20_off + 986 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_yy = cbuffer.data(id_geom_20_off + 987 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_yz = cbuffer.data(id_geom_20_off + 988 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_zz = cbuffer.data(id_geom_20_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzzz_xx, g_zz_0_yyyzzz_xy, g_zz_0_yyyzzz_xz, g_zz_0_yyyzzz_yy, g_zz_0_yyyzzz_yz, g_zz_0_yyyzzz_zz, g_zz_0_yyzzz_xx, g_zz_0_yyzzz_xxy, g_zz_0_yyzzz_xy, g_zz_0_yyzzz_xyy, g_zz_0_yyzzz_xyz, g_zz_0_yyzzz_xz, g_zz_0_yyzzz_yy, g_zz_0_yyzzz_yyy, g_zz_0_yyzzz_yyz, g_zz_0_yyzzz_yz, g_zz_0_yyzzz_yzz, g_zz_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzzz_xx[k] = -g_zz_0_yyzzz_xx[k] * ab_y + g_zz_0_yyzzz_xxy[k];

                g_zz_0_yyyzzz_xy[k] = -g_zz_0_yyzzz_xy[k] * ab_y + g_zz_0_yyzzz_xyy[k];

                g_zz_0_yyyzzz_xz[k] = -g_zz_0_yyzzz_xz[k] * ab_y + g_zz_0_yyzzz_xyz[k];

                g_zz_0_yyyzzz_yy[k] = -g_zz_0_yyzzz_yy[k] * ab_y + g_zz_0_yyzzz_yyy[k];

                g_zz_0_yyyzzz_yz[k] = -g_zz_0_yyzzz_yz[k] * ab_y + g_zz_0_yyzzz_yyz[k];

                g_zz_0_yyyzzz_zz[k] = -g_zz_0_yyzzz_zz[k] * ab_y + g_zz_0_yyzzz_yzz[k];
            }

            /// Set up 990-996 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzzz_xx = cbuffer.data(id_geom_20_off + 990 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_xy = cbuffer.data(id_geom_20_off + 991 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_xz = cbuffer.data(id_geom_20_off + 992 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_yy = cbuffer.data(id_geom_20_off + 993 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_yz = cbuffer.data(id_geom_20_off + 994 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_zz = cbuffer.data(id_geom_20_off + 995 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzzz_xx, g_zz_0_yyzzzz_xy, g_zz_0_yyzzzz_xz, g_zz_0_yyzzzz_yy, g_zz_0_yyzzzz_yz, g_zz_0_yyzzzz_zz, g_zz_0_yzzzz_xx, g_zz_0_yzzzz_xxy, g_zz_0_yzzzz_xy, g_zz_0_yzzzz_xyy, g_zz_0_yzzzz_xyz, g_zz_0_yzzzz_xz, g_zz_0_yzzzz_yy, g_zz_0_yzzzz_yyy, g_zz_0_yzzzz_yyz, g_zz_0_yzzzz_yz, g_zz_0_yzzzz_yzz, g_zz_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzzz_xx[k] = -g_zz_0_yzzzz_xx[k] * ab_y + g_zz_0_yzzzz_xxy[k];

                g_zz_0_yyzzzz_xy[k] = -g_zz_0_yzzzz_xy[k] * ab_y + g_zz_0_yzzzz_xyy[k];

                g_zz_0_yyzzzz_xz[k] = -g_zz_0_yzzzz_xz[k] * ab_y + g_zz_0_yzzzz_xyz[k];

                g_zz_0_yyzzzz_yy[k] = -g_zz_0_yzzzz_yy[k] * ab_y + g_zz_0_yzzzz_yyy[k];

                g_zz_0_yyzzzz_yz[k] = -g_zz_0_yzzzz_yz[k] * ab_y + g_zz_0_yzzzz_yyz[k];

                g_zz_0_yyzzzz_zz[k] = -g_zz_0_yzzzz_zz[k] * ab_y + g_zz_0_yzzzz_yzz[k];
            }

            /// Set up 996-1002 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzzz_xx = cbuffer.data(id_geom_20_off + 996 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_xy = cbuffer.data(id_geom_20_off + 997 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_xz = cbuffer.data(id_geom_20_off + 998 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_yy = cbuffer.data(id_geom_20_off + 999 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_yz = cbuffer.data(id_geom_20_off + 1000 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_zz = cbuffer.data(id_geom_20_off + 1001 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzzz_xx, g_zz_0_yzzzzz_xy, g_zz_0_yzzzzz_xz, g_zz_0_yzzzzz_yy, g_zz_0_yzzzzz_yz, g_zz_0_yzzzzz_zz, g_zz_0_zzzzz_xx, g_zz_0_zzzzz_xxy, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_xyy, g_zz_0_zzzzz_xyz, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_yy, g_zz_0_zzzzz_yyy, g_zz_0_zzzzz_yyz, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_yzz, g_zz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzzz_xx[k] = -g_zz_0_zzzzz_xx[k] * ab_y + g_zz_0_zzzzz_xxy[k];

                g_zz_0_yzzzzz_xy[k] = -g_zz_0_zzzzz_xy[k] * ab_y + g_zz_0_zzzzz_xyy[k];

                g_zz_0_yzzzzz_xz[k] = -g_zz_0_zzzzz_xz[k] * ab_y + g_zz_0_zzzzz_xyz[k];

                g_zz_0_yzzzzz_yy[k] = -g_zz_0_zzzzz_yy[k] * ab_y + g_zz_0_zzzzz_yyy[k];

                g_zz_0_yzzzzz_yz[k] = -g_zz_0_zzzzz_yz[k] * ab_y + g_zz_0_zzzzz_yyz[k];

                g_zz_0_yzzzzz_zz[k] = -g_zz_0_zzzzz_zz[k] * ab_y + g_zz_0_zzzzz_yzz[k];
            }

            /// Set up 1002-1008 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzzz_xx = cbuffer.data(id_geom_20_off + 1002 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_xy = cbuffer.data(id_geom_20_off + 1003 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_xz = cbuffer.data(id_geom_20_off + 1004 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_yy = cbuffer.data(id_geom_20_off + 1005 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_yz = cbuffer.data(id_geom_20_off + 1006 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_zz = cbuffer.data(id_geom_20_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_xx, g_z_0_zzzzz_xy, g_z_0_zzzzz_xz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_zz, g_zz_0_zzzzz_xx, g_zz_0_zzzzz_xxz, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_xyz, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_xzz, g_zz_0_zzzzz_yy, g_zz_0_zzzzz_yyz, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_yzz, g_zz_0_zzzzz_zz, g_zz_0_zzzzz_zzz, g_zz_0_zzzzzz_xx, g_zz_0_zzzzzz_xy, g_zz_0_zzzzzz_xz, g_zz_0_zzzzzz_yy, g_zz_0_zzzzzz_yz, g_zz_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzzz_xx[k] = -2.0 * g_z_0_zzzzz_xx[k] - g_zz_0_zzzzz_xx[k] * ab_z + g_zz_0_zzzzz_xxz[k];

                g_zz_0_zzzzzz_xy[k] = -2.0 * g_z_0_zzzzz_xy[k] - g_zz_0_zzzzz_xy[k] * ab_z + g_zz_0_zzzzz_xyz[k];

                g_zz_0_zzzzzz_xz[k] = -2.0 * g_z_0_zzzzz_xz[k] - g_zz_0_zzzzz_xz[k] * ab_z + g_zz_0_zzzzz_xzz[k];

                g_zz_0_zzzzzz_yy[k] = -2.0 * g_z_0_zzzzz_yy[k] - g_zz_0_zzzzz_yy[k] * ab_z + g_zz_0_zzzzz_yyz[k];

                g_zz_0_zzzzzz_yz[k] = -2.0 * g_z_0_zzzzz_yz[k] - g_zz_0_zzzzz_yz[k] * ab_z + g_zz_0_zzzzz_yzz[k];

                g_zz_0_zzzzzz_zz[k] = -2.0 * g_z_0_zzzzz_zz[k] - g_zz_0_zzzzz_zz[k] * ab_z + g_zz_0_zzzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

